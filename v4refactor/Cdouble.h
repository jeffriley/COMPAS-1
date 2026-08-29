#pragma once

// CDOUBLE - a double wrapper that uses tolerance-aware comparisons for all
// relational operators (==, !=, <, <=, >, >=). Wrap stellar attributes that
// need tolerance-aware branch decisions (mass, radius, etc.) - the code will
// read like an ordinary numerical comparison, but the comparison ignores
// differences smaller than the tolerance.
//
// Design choices:
//
//   - Global static tolerance. All CDOUBLE instances share a single absolute
//     and relative tolerance. This is symmetric by construction: `a == b` and
//     `b == a` always agree. Set via the static accessors
//     CDOUBLE::AbsoluteTolerance(x) and CDOUBLE::RelativeTolerance(y).
//
//   - Arithmetic is CDOUBLE-in, CDOUBLE-out.  The class does not provide an
//     implicit conversion to double - intermediate results stay as CDOUBLEs.
//     Use Value() or static_cast<double>(c) to extract a raw double if required.
//
//   - Implicit construction from arithmetic types is allowed.  This is what lets
//     `mass < 1.4` work without requiring a CDOUBLE wrapper at the call site
//     (1.4 is promoted to a temporary CDOUBLE for the comparison).
//
//   - Tolerance-aware <, <=, >, >= are not transitive in general (this is a known
//     property of tolerance-based comparison).  Do *not* use CDOUBLE as a key in
//     std::set / std::map / std::sort - the strict-weak-ordering requirement of
//     those algorithms is not satisfied.
//
//   - Binary operators are defined as hidden friends so that `cdouble OP arithmetic`
//     and `arithmetic OP cdouble` both work symmetrically via implicit conversion of
//     the arithmetic operand.

#include <algorithm>
#include <cmath>
#include <ostream>
#include <type_traits>

#include "constants.h"


class CDOUBLE final {

private:

    // State

    double m_Value;                                                                     // The wrapped value

    static inline double s_AbsoluteTolerance = FLOAT_TOLERANCE_ABSOLUTE;                // Global absolute tolerance (set via AbsoluteTolerance(x))
    static inline double s_RelativeTolerance = FLOAT_TOLERANCE_RELATIVE;                // Global relative tolerance (set via RelativeTolerance(x))


    // Tolerance-aware comparison primitives
    //
    // These operate on raw doubles and use the static class tolerances.
    // Symmetric in their arguments: equals(a, b) == equals(b, a) etc.
    //
    // Note: smaller / bigger respect the tolerance band. If two values are within
    // tolerance, neither is smaller nor bigger - they are equal. This is what makes
    // operator< NON-transitive (which is fine for branch-selection use cases, but
    // do *not* use as a strict weak ordering).
    GNU_PURE static bool equals(double a, double b) {                                   // Tolerance-aware equality
        return a == b
            || std::fabs(a - b) <= std::max(s_AbsoluteTolerance,
                                            s_RelativeTolerance * std::max(std::fabs(a), std::fabs(b)));
    }
    GNU_PURE static bool smaller(double a, double b) { return !equals(a, b) && a < b; } // Tolerance-aware strict less-than
    GNU_PURE static bool bigger (double a, double b) { return !equals(a, b) && a > b; } // Tolerance-aware strict greater-than


public:

    // Constructors, special members

    // Implicit conversion from arithmetic types is allowed (no `explicit`) so that
    // `mass < 1.4` works naturally - the 1.4 becomes a temporary CDOUBLE for the
    // comparison.
    CDOUBLE(double p_Value = DEFAULT_INITIAL_DOUBLE_VALUE) : m_Value(p_Value) { }       // (Implicit) construct from a double

    // Rule-of-five: all five special members explicitly defaulted.
    //
    // CDOUBLE is a trivial value type (a single double), so the compiler-generated
    // copy/move/destroy semantics are what we want.  We still spell them out because
    // declaring any one of them (here the destructor, and the copy/move assignments)
    // makes the *implicit* generation of the others deprecated under -Wdeprecated-copy
    // or -Wdeprecated-copy-dtor.  Defaulting all five keeps the trivial semantics while
    // silencing the deprecation warnings.
    CDOUBLE(const CDOUBLE&)            = default;
    CDOUBLE(CDOUBLE&&)                 = default;
    CDOUBLE& operator=(const CDOUBLE&) = default;
    CDOUBLE& operator=(CDOUBLE&&)      = default;
    ~CDOUBLE()                         = default;


    // Value access
    GNU_PURE double Value() const noexcept               { return m_Value; }            // Raw double value
    GNU_PURE explicit operator double() const noexcept   { return m_Value; }            // Explicit conversion: static_cast<double>(c) works, `double d = c;` does *not*


    // Global tolerance (static API)
    //
    // Getters and setters share names, distinguished by signature (the setter takes
    // an argument).  Pass-through to the static members.
    // The setter ensures the tolerance is non-negative by taking the absolute value
    // of the value passed (negative tolerances would invert the comparison semantics).
    GNU_PURE static double AbsoluteTolerance()                   { return s_AbsoluteTolerance; }
    GNU_PURE static double RelativeTolerance()                   { return s_RelativeTolerance; }
             static void   AbsoluteTolerance(double p_Tolerance) { s_AbsoluteTolerance = std::abs(p_Tolerance); }
             static void   RelativeTolerance(double p_Tolerance) { s_RelativeTolerance = std::abs(p_Tolerance); }


    // Assignment from arithmetic types
    //
    // The default copy/move assignment operators (above) handle CDOUBLE-to-CDOUBLE.
    // This template handles direct assignment from any built-in arithmetic type
    // without first constructing a CDOUBLE temporary.
    template<typename T>
    CDOUBLE& operator=(const T p_Value) noexcept {
        static_assert(std::is_arithmetic_v<T>, "CDOUBLE assignment requires an arithmetic type");
        m_Value = static_cast<double>(p_Value);
        return *this;
    }


    // Compound assignment
    //
    // Take by const ref - mixed-type cases (e.g. `c += 1.0`) work via implicit construction
    // of a temporary CDOUBLE.
    CDOUBLE& operator+=(const CDOUBLE& p_Other) noexcept { m_Value += p_Other.m_Value; return *this; }
    CDOUBLE& operator-=(const CDOUBLE& p_Other) noexcept { m_Value -= p_Other.m_Value; return *this; }
    CDOUBLE& operator*=(const CDOUBLE& p_Other) noexcept { m_Value *= p_Other.m_Value; return *this; }
    CDOUBLE& operator/=(const CDOUBLE& p_Other) noexcept { m_Value /= p_Other.m_Value; return *this; }


    // Unary operators
    GNU_PURE CDOUBLE operator-() const noexcept { return CDOUBLE(-m_Value); }                    // Negation
    GNU_PURE CDOUBLE operator+() const noexcept { return *this; }                                // Unary plus (no-op)


    // Binary arithmetic operators (hidden friends)
    //
    // Defined inside the class so they have access to private members and are found via ADL.
    // These take operands by const ref.  Symmetric in the operand position because they are
    // non-members - mixed-type cases like `1.0 + c` and `c + 1.0` resolve identically (the
    // arithmetic operand is implicitly converted to a temporary CDOUBLE).
    GNU_PURE friend CDOUBLE operator+(const CDOUBLE& p_Lhs, const CDOUBLE& p_Rhs) noexcept { return CDOUBLE(p_Lhs.m_Value + p_Rhs.m_Value); }
    GNU_PURE friend CDOUBLE operator-(const CDOUBLE& p_Lhs, const CDOUBLE& p_Rhs) noexcept { return CDOUBLE(p_Lhs.m_Value - p_Rhs.m_Value); }
    GNU_PURE friend CDOUBLE operator*(const CDOUBLE& p_Lhs, const CDOUBLE& p_Rhs) noexcept { return CDOUBLE(p_Lhs.m_Value * p_Rhs.m_Value); }
    GNU_PURE friend CDOUBLE operator/(const CDOUBLE& p_Lhs, const CDOUBLE& p_Rhs) noexcept { return CDOUBLE(p_Lhs.m_Value / p_Rhs.m_Value); }


    // Comparison operators (hidden friends)
    //
    // Tolerance-aware. Symmetric in operand position because they are non-members and call
    // the (symmetric) static comparison primitives.
    GNU_PURE friend bool operator==(const CDOUBLE& p_Lhs, const CDOUBLE& p_Rhs) noexcept { return  equals (p_Lhs.m_Value, p_Rhs.m_Value); }
    GNU_PURE friend bool operator!=(const CDOUBLE& p_Lhs, const CDOUBLE& p_Rhs) noexcept { return !equals (p_Lhs.m_Value, p_Rhs.m_Value); }
    GNU_PURE friend bool operator< (const CDOUBLE& p_Lhs, const CDOUBLE& p_Rhs) noexcept { return  smaller(p_Lhs.m_Value, p_Rhs.m_Value); }
    GNU_PURE friend bool operator<=(const CDOUBLE& p_Lhs, const CDOUBLE& p_Rhs) noexcept { return !bigger (p_Lhs.m_Value, p_Rhs.m_Value); }
    GNU_PURE friend bool operator> (const CDOUBLE& p_Lhs, const CDOUBLE& p_Rhs) noexcept { return  bigger (p_Lhs.m_Value, p_Rhs.m_Value); }
    GNU_PURE friend bool operator>=(const CDOUBLE& p_Lhs, const CDOUBLE& p_Rhs) noexcept { return !smaller(p_Lhs.m_Value, p_Rhs.m_Value); }


    // Stream output
    friend std::ostream& operator<<(std::ostream& p_Os, const CDOUBLE& p_C) { return p_Os << p_C.m_Value; }
};
