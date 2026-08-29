#pragma once


#include <memory>
#include <algorithm>
#include <type_traits>

#include "constants.h"
#include "typedefs.h"
#include "Rand.h"
#include "ErrorCatalog.h"

namespace utils {

/*
 * Utility functions
 * 
 * The "utils" namespace gathers general functions that have utility across the
 * COMPAS code base.
 * 
 * These are free functions that:
 * 
 *    (a) don't really belong to any of the stellar or binary classes,
 *    (b) aren't required to run against an instantiated object, and
 *    (c) rely only on:
 *             (i) parameters passed, and/or
 *            (ii) constants.h, and/or
 *           (iii) the GLOBALS object, and/or
 *             (v) the RAND object, and/or
 *             (v) library functions, and/or
 *            (vi) other free functions (defined here)
 */


[[noreturn]] void Fail(const ERROR p_Error, const StrT p_Preamble = "");

/*
 * GetValue
 *
 * @brief
 * Returns the underlying value of a std::optional, or calls Fail() if the
 * optional has no value.
 */
template<typename T>
T GetValue(const std::optional<T>& p_Var) {
    if (p_Var.has_value()) return p_Var.value();
    Fail(ERROR::NO_VALUE);
}


GNU_PURE  SizeTVectorT BinarySearch(const DblVectorT& p_Arr, const double p_X);

GNU_CONST bool         BracketTolerance(const double p_Bracket1, const double p_Bracket2);

GNU_PURE  StrT         CentreJustify(const StrT p_Str, const SizeT p_Width);

GNU_CONST int          Compare(const double p_X, const double p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true);
GNU_PURE  int          Compare(const CDOUBLE& p_X, const CDOUBLE& p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true);
GNU_PURE  int          Compare(const CDOUBLE& p_X, const double p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true);
GNU_PURE  int          Compare(const double p_X, const CDOUBLE& p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true);
 
GNU_PURE  bool         Equals(StrT p_Str1, StrT p_Str2);

          bool         FileExists(const char *p_Filename);
          bool         FileExists(const StrT& p_Filename);

GNU_PURE  StrT         GetBoostVersion();
          StrT         GetGslVersion();
GNU_PURE  StrT         GetHDF5Version();

          StrVectorT   GetStackTrace();

GNU_CONST constexpr double IntPow(const double p_Base, const int p_Exponent);
GNU_CONST constexpr double IntPowL(const double p_Base, const int p_Exponent);

GNU_PURE  SizeT        IsBOOL(const StrT p_Str);
GNU_PURE  bool         IsDOUBLE(const StrT p_Str);
GNU_PURE  bool         IsFLOAT(const StrT p_Str);
GNU_PURE  bool         IsINT(const StrT p_Str);
GNU_PURE  bool         IsLONGDOUBLE(const StrT p_Str);
GNU_PURE  bool         IsLONGINT(const StrT p_Str);
GNU_PURE  bool         IsULONGINT(const StrT p_Str);

GNU_PURE  bool         IsOneOf(const STELLAR_TYPE p_StellarType, const StellarTypeListT p_List);

GNU_PURE  StrT         PadLeadingZeros(const StrT p_Str, const SizeT p_MaxLength);
GNU_PURE  StrT         PadTrailingSpaces(const StrT p_Str, const SizeT p_MaxLength);

          void         ShowStackTrace();
GNU_CONST SN_EVENT     SNEventType(const SN_EVENT p_SNEvent);
          StrT         SplashScreen(const bool p_Print);

          StrT&        ltrim(StrT& p_Str);
          StrT&        rtrim(StrT& p_Str);
          StrT&        trim(StrT& p_Str);

GNU_PURE  StrT         ToLower(StrT p_Str);
GNU_PURE  StrT         ToUpper(StrT p_Str);

std::tuple<ERROR, StrT, StrVectorT> CreateDirectories(const StrT p_Path);
std::tuple<ERROR, StrT, StrVectorT> RemoveDirectories(const StrVectorT p_Paths);

std::tuple<ERROR, DblVectorT>       ReadTimesteps(const StrT p_TimestepsFileName);

GNU_CONST std::tuple<ERROR, double> SolveQuadratic(const double p_A, const double p_B, double p_C);


// Template functions
// Must be in header

/*
 * Find
 *
 * @brief
 * Find an element in a vector (of any type).
 *
 * Determines if an element is contained within a vector.
 * Returns a tuple containing:
 *
 *   - a boolean indicating the result (true = found, false = not found)
 *   - the vector index of the element if found (-1 if not found)
 *
 *
 * std::tuple<bool, int> Find(const T &p_Elem, const std::vector<T> &p_Vector)
 *
 * @param       p_Elem                          The element to find in the vector
 * @param       p_Vector                        The vector in which to look for p_Elem
 * @return                                      Tuple<bool, int> indicating <found, element index>, with index = -1 if found = false
 */
template <typename T>
std::tuple<bool, int> Find(const T &p_Elem, const std::vector<T> &p_Vector) {
    auto iter = std::find(p_Vector.begin(), p_Vector.end(), p_Elem);                                                        // Find the element
    return iter != p_Vector.end() ? std::make_tuple(true, distance(p_Vector.begin(), iter)) : std::make_tuple(false, -1l);  // If found return index, otherwise -1
}


/*
 * GetMapKey
 *
 * @brief
 * Find a value in an unordered map and return the key if found, otherwise default value
 *
 * This function looks for the passed string value in an unordered map, and if the string
 * is found returns the key corresponding to the value found.  If the value is not found
 * the value passed as the default value is returned.
 *
 * The string comparison is case-insensitive.
 *
 * The default value passed must allow the functional return type to be deduced - so it
 * should be one of the values in the unordered_map.
 *
 *
 * template<typename M, typename E>
 * std::tuple<bool, E> GetMapKey(const std::string p_Value, const M& p_Map, const E& p_Default)
 *
 * @param       p_Value                         The value to find in the unordered map
 * @param       p_Map                           The unordered map in which to look for p_Value
 * @param       p_Default                       The default value to be returned if p_Value is not found
 * @return                                      Tuple<bool, E> indicationg <found, key>, where
 *                                                 E is the type of p_Default
 *                                                 found is a boolean: true if p_Value was found, false is not
 *                                                 key is:
 *                                                    the key corresponding to the value found, or
 *                                                    p_Default if the value was not found
 */
template<typename M, typename E>
std::tuple<bool, E> GetMapKey(const StrT p_Value, const M& p_Map, const E& p_Default) {
    for (auto& it: p_Map)
        if (Equals(it.second, p_Value)) return std::make_tuple(true, it.first);
    return std::make_tuple(false, p_Default);
}


/*
 * to_underlying
 *
 * Cast enum element to underlying type.
 * Name is snake_case here because C++23 has a std::to_underlying(), so once
 * we move to C++23 or later we can just switch all utils::to_underlying()
 * call so std::to_underlying().
 * 
 * template <typename E>
 * to_underlying(E e)
 * 
 * @param       e                               Enum element to per cast to underlying type
 * @return                                      static_cast<T>(E::e)
 */
template <typename E>
constexpr auto to_underlying(E e) noexcept {
    static_assert(std::is_enum_v<E>, "to_underlying requires an enum type");
    return static_cast<std::underlying_type_t<E>>(e);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 DECODE / CONVERT                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * SNEventType
 *
 * @brief
 * Return a single SN type based on the SN_EVENT parameter passed.
 * 
 * Returns (in priority order):
 *
 *    SN_EVENT::NONE    iff no bits are set
 *    SN_EVENT::CCSN    iff CCSN  bit is set and USSN bit is not set
 *    SN_EVENT::ECSN    iff ECSN  bit is set
 *    SN_EVENT::PISN    iff PISN  bit is set
 *    SN_EVENT::PPISN   iff PPISN bit is set
 *    SN_EVENT::USSN    iff USSN  bit is set
 *    SN_EVENT::AIC     iff AIC   bit is set
 *    SN_EVENT::SNIA    iff SNIA  bit is set and HeSD bit is not set
 *    SN_EVENT::HeSD    iff HeSD  bit is set
 *    SN_EVENT::UNKNOWN otherwise
 * 
 * 
 * SN_EVENT SNEventType(const SN_EVENT p_SNEvent)
 *
 * @param       p_SNEvent                       SN_EVENT mask to check for SN event type
 * @return                                      SN_EVENT
 */
inline SN_EVENT SNEventType(const SN_EVENT p_SNEvent) {

    if (p_SNEvent == SN_EVENT::NONE)                                        return SN_EVENT::NONE;

    if ((p_SNEvent & (SN_EVENT::CCSN | SN_EVENT::USSN)) == SN_EVENT::CCSN ) return SN_EVENT::CCSN;
    if ((p_SNEvent & SN_EVENT::ECSN )                   == SN_EVENT::ECSN ) return SN_EVENT::ECSN;
    if ((p_SNEvent & SN_EVENT::PISN )                   == SN_EVENT::PISN ) return SN_EVENT::PISN;
    if ((p_SNEvent & SN_EVENT::PPISN)                   == SN_EVENT::PPISN) return SN_EVENT::PPISN;
    if ((p_SNEvent & SN_EVENT::USSN )                   == SN_EVENT::USSN ) return SN_EVENT::USSN;
    if ((p_SNEvent & SN_EVENT::AIC  )                   == SN_EVENT::AIC  ) return SN_EVENT::AIC;
    if ((p_SNEvent & (SN_EVENT::SNIA | SN_EVENT::HeSD)) == SN_EVENT::SNIA ) return SN_EVENT::SNIA;
    if ((p_SNEvent & SN_EVENT::HeSD )                   == SN_EVENT::HeSD ) return SN_EVENT::HeSD;

    return SN_EVENT::UNKNOWN;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                   FIND / SEARCH                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * BinarySearch
 *
 * @brief
 * Binary search function.
 * Searches a std::vector of doubles for a given double value.  The vector of doubles is
 * expected to be sorted before calling this function (sorting here would be prohibitive
 * performance-wise - the expectation is that this function could/will be called many
 * times for the same (sorted) vector).
 * 
 * For a given double x, and a sorted array arr, the function returns a vector of (two)
 * numbers (of type SizeT) that represent the lower and upper bin edges of x in arr.
 *
 *
 * SizeTVectorT BinarySearch(const DblVectorT& p_Arr, const double p_X)
 *
 * @param       p_Arr                           Sorted array to search over
 * @param       p_X                             Value to search for
 * @return                                      Vector containing SizeT indices of the lower and upper bin edges that contain p_X.
 *                                                  If x < min(Arr), return {-1, 0}
 *                                                  If x > max(Arr), return {0, -1}
 *                                                  If x is equal to an array element, return index of that element i.e. {idx, idx}
 *                                                  Otherwise x is between two array elements, and {low, low+1} is returned (p_Arr[low] < x < p_Arr[low+1])
 */
inline SizeTVectorT BinarySearch(const DblVectorT& p_Arr, const double p_X) {

    auto it = std::ranges::lower_bound(p_Arr, p_X);                     // First element not less than p_X

    if (it == p_Arr.end()) { return {0, static_cast<SizeT>(-1)}; }      // p_X > max(p_Arr) - also covers empty p_Arr

    const SizeT idx = static_cast<SizeT>(it - p_Arr.begin());           // Index of *it

    if (*it == p_X)          { return {idx, idx}; }                     // Exact match: low = up = idx
    if (it == p_Arr.begin()) { return {static_cast<SizeT>(-1), 0}; }    // *it > p_X with nothing before -> p_X < min(p_Arr)

    return {idx - 1, idx};                                              // p_Arr[idx-1] < p_X < p_Arr[idx]
}


/*
 * IsOneOf
 *
 * @brief
 * Determine if the stellar type passed is one of a list of stellar types passed.
 *
 *
 * bool IsOneOf(const STELLAR_TYPE p_StellarType, const std::initializer_list<ST> p_List)
 *
 * @param       p_StellarType                   Stellar type to check
 * @param       p_List                          List of stellar types
 * @return                                      Boolean - true if p_StellarType is in list, false if not
 */
inline bool IsOneOf(const STELLAR_TYPE p_StellarType, const StellarTypeListT p_List) {
    for (auto elem: p_List)
        if (p_StellarType == elem) return true;
    return false;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     NUMERICS                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * Compare
 *
 * @brief
 * Tolerance-aware comparison of two floating-point values.
 *
 * Returns -1, 0, or +1 depending on whether p_X is less than, equal to (within tolerance), or
 * greater than p_Y.
 *
 * Four variant signatures are implemented here:
 * 
 * int Compare(const double p_X, const double p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true)
 * int Compare(const CDOUBLE& p_X, const double p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true)
 * int Compare(const double p_X, const CDOUBLE& p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true)
 * int Compare(const CDOUBLE& p_X, const CDOUBLE& p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true)
 * 
 * This allows the function to be called with both operands of any arithmetic type, both
 * type 'CDOUBLE', or mixed 'CDOUBLE' and any arithmetic type.
 * 
 * The global absolute and relative tolerances are declared in constants.h, and are used
 * for non-CDOUBLE arithmetic types.  The CDOUBLE absolute and relative tolerances are
 * declared inside the CDOUBLE class (in CDouble.h). The *effective* tolerances for all 
 * variants of the function are determined as follows:
 * 
 *    - if p_Tolerance >= 0.0, it overrides both the global tolerances *and* the CDOUBLE
 *      tolerances.  Whether the absolute or relative tolerances are overridden is determined
 *      by the p_Absolute selector: absolute if true, relative if false.
 * 
 *    - if p_Tolerance < 0.0 (the default):
 * 
 *       - if both operands are non-CDOUBLE arithmetic types, FLOAT_TOLERANCE_ABSOLUTE and
 *         FLOAT_TOLERANCE_RELATIVE (from constants.h) are used as the effective tolerances.
 * 
 *       - if both operands are CDOUBLE types, CDOUBLE::AbsoluteTolerance() and
 *         CDOUBLE::RelativeTolerance() (from Cdouble.h) are used as the effective tolerances.
 * 
 *       - if the operands are mixed (i.e. one is a non-CDOUBLE arithmetic type, and the other
 *         is CDOUBLE) the effective tolerances are *smaller* (element-wise, for absolute and
 *         relative) of the global tolerances (from constants.h) and the CDOUBLE tolerances
 *         from CDouble.h (i.e. CDOUBLE::AbsoluteTolerance() and CDOUBLE::RelativeTolerance()).
 *
 * Kill switch:
 * To turn the tolerance-aware comparison off (i.e. utils::Compare() just becomes an exact
 * comparison), set FLOAT_TOLERANCE_ABSOLUTE and FLOAT_TOLERANCE_RELATIVE to 0.0 in constants.h.
 * The min() means 0 dominates, so it becomes an exact equality.
 * 
 *
 * For each of the variant signatures below:
 *
 * @param       p_X                             Floating-point value to be compared
 * @param       p_Y                             Floating-point value to be compared
 * @param       p_Tolerance                     Floating-point tolerance value - if > 0.0 supersedes global tolerance
 * @param       p_Absolute                      Boolean indicating whether p_Tolerance should be treated as absolute
 *                                              tolerance (true) or relative tolerance (false)
 * @return                                      Integer indicating result of comparison:
 *                                                  -1 indicates p_X is less than p_Y
 *                                                   0 indicates equality
 *                                                  +1 indicates p_X is greater than p_Y
 */
inline int CompareImpl(double p_X, double p_Y, double p_AbsTol, double p_RelTol) {
    const double diff = std::fabs(p_X - p_Y);
    const double tol  = std::max(p_AbsTol, p_RelTol * std::max(std::fabs(p_X), std::fabs(p_Y)));

    if (diff <= tol) return 0;
    else             return p_X < p_Y ? -1 : 1; // 'else' is redundant, but I think reads better...
}

inline int Compare(const double p_X, const double p_Y, const double p_Tolerance, const bool p_Absolute) {

    const bool overrideTol = p_Tolerance >= 0.0;

    const double absTol = overrideTol ? (p_Absolute ? p_Tolerance : 0.0        ) : FLOAT_TOLERANCE_ABSOLUTE;
    const double relTol = overrideTol ? (p_Absolute ? 0.0         : p_Tolerance) : FLOAT_TOLERANCE_RELATIVE;

    return CompareImpl(p_X, p_Y, absTol, relTol);
}

inline int Compare(const CDOUBLE& p_X, const CDOUBLE& p_Y, const double p_Tolerance, const bool p_Absolute) {

    const bool overrideTol = p_Tolerance >= 0.0;

    const double absTol = overrideTol ? (p_Absolute ? p_Tolerance : 0.0        ) : CDOUBLE::AbsoluteTolerance();
    const double relTol = overrideTol ? (p_Absolute ? 0.0         : p_Tolerance) : CDOUBLE::RelativeTolerance();
        
    return CompareImpl(p_X.Value(), p_Y.Value(), absTol, relTol);
}

inline int Compare(const CDOUBLE& p_X, const double p_Y, const double p_Tolerance, const bool p_Absolute) {

    const bool overrideTol = p_Tolerance >= 0.0;

    const double absTol = overrideTol ? (p_Absolute ? p_Tolerance : 0.0        ) : std::min(CDOUBLE::AbsoluteTolerance(), FLOAT_TOLERANCE_ABSOLUTE);
    const double relTol = overrideTol ? (p_Absolute ? 0.0         : p_Tolerance) : std::min(CDOUBLE::RelativeTolerance(), FLOAT_TOLERANCE_RELATIVE);

    return CompareImpl(p_X.Value(), p_Y, absTol, relTol);
}

inline int Compare(const double p_X, const CDOUBLE& p_Y, const double p_Tolerance, const bool p_Absolute) {

    const bool overrideTol = p_Tolerance >= 0.0;

    const double absTol = overrideTol ? (p_Absolute ? p_Tolerance : 0.0        ) : std::min(CDOUBLE::AbsoluteTolerance(), FLOAT_TOLERANCE_ABSOLUTE);
    const double relTol = overrideTol ? (p_Absolute ? 0.0         : p_Tolerance) : std::min(CDOUBLE::RelativeTolerance(), FLOAT_TOLERANCE_RELATIVE);

    return CompareImpl(p_X, p_Y.Value(), absTol, relTol);
}


/*
 * IntPow
 *
 * @brief
 * Calculate x^y where x is double and y is an integer (positive or negative)
 *
 * Faster than std::pow() for integer exponent, n.  This function uses plain multiplication,
 * andcompletes in O(|n|) multiplications, and is the right choice for small |n| (|n| <= ~4).
 * For n > ~4, consider IntPowL() (below).
 *
 *
 * double IntPow(const double p_Base, const int p_Exponent)
 *
 * @param       p_Base                          Base - number to be raised to integer power
 * @param       p_Exponent                      Exponent - integer to which base should be raised
 * @return                                      Base ^ Exponent
 */
constexpr double IntPow(const double p_Base, const int p_Exponent) {
    double result = 1.0;                                    // Default (exponent = 0)
    int i = (p_Exponent < 0) ? -p_Exponent : p_Exponent;    // Cache absolute exponent (abs() not constexpr until C++23)
    while (i-- > 0) result *= p_Base;                       // Multiply
    return p_Exponent < 0 ? 1.0 / result : result;          // Invert if negative exponent
}


/*
 * IntPowL
 *
 * @brief
 * Calculate x^y where x is double and y is an integer (positive or negative)
 *
 * Faster than std::pow() for integer exponent, n.  This function uses exponentiation by squaring,
 * and completes in O(log |n|) multiplications, and is the right choice for large |n| (|n| > ~4).
 * For n <= ~4, consider IntPow() (above).
 *
 *
 * double IntPowL(const double p_Base, const int p_Exponent)
 *
 * @param       p_Base                          Base - number to be raised to integer power
 * @param       p_Exponent                      Exponent - integer to which base should be raised
 * @return                                      Base ^ Exponent
 */
constexpr double IntPowL(const double p_Base, const int p_Exponent) {
    double result = 1.0;                                    // Default (exponent = 0)
    double base   = p_Base;                                 // Working base, squared each iteration
    int i = (p_Exponent < 0) ? -p_Exponent : p_Exponent;    // Bits of |p_Exponent| still to process (abs() not constexpr until C++23)
    while (i > 0) {                                         // Loop over bits of |p_Exponent|
        if (i & 1) result *= base;                          // Current bit set: accumulate this power
        base *= base;                                       // Square base for the next bit position
        i >>= 1;                                            // Shift to next bit
    }
    return p_Exponent < 0 ? 1.0 / result : result;          // Invert if negative exponent
}


/*
 * BracketTolerance
 *
 * @brief
 * Tolerance for Boost bracket_and_solve_root()
 *
 * Determines if the brackets around the root are within the COMPAS defined tolerance.
 * 
 * 
 * bool BracketTolerance(const double p_Bracket1, const double p_Bracket2)
 * 
 * @param       p_Bracket1                      Bracket bound 1
 * @param       p_Bracket2                      Bracket bound 2
 * @return                                      Boolean indicating if the brackets bounds are within tolerance
 */
inline bool BracketTolerance(const double p_Bracket1, const double p_Bracket2) {
    const double diff = std::fabs(p_Bracket1 - p_Bracket2);                                       // Absolute value of difference
    const double min  = std::min(p_Bracket1, p_Bracket2);                                         // Minimum bracket value - could straddle 0.0
    return diff <= ROOT_ABS_TOLERANCE || std::fabs(diff / min) <= ROOT_REL_TOLERANCE;
}


} // namespace utils
