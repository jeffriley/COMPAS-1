#ifndef __Cdouble_H__
#define __Cdouble_H__

#include <type_traits>

#include "constants.h"

class CDOUBLE {

#ifdef CDOUBLE_USE_TOLERANCE
#define EQUALS(a, b) (a == b || (std::fabs(a - b) <= std::max(m_AbsoluteTolerance, m_RelativeTolerance * std::max(std::fabs(a), std::fabs(b)))))
#else
#define EQUALS(a, b) (a == b)
#endif

#define SMALLER(a, b) (EQUALS(a, b) ? false : (a < b) ? true : false)
#define BIGGER(a, b)  (EQUALS(a, b) ? false : (b < a) ? true : false)

#define AssignmentMsg "Unsupported type for CDOUBLE assignment"
#define ComparisonMsg "Unsupported type for CDOUBLE comparison"

#define IS(x) std::is_same<T, x>::value
#define AssertType(T, msg) static_assert(IS(float) || IS(double) || IS(long double) || IS(signed char) || IS(unsigned char) || IS(short int) || IS(unsigned short int) || IS(int) || IS(long unsigned int) || IS(long int) || IS(unsigned long int) || IS(long long int) || IS(unsigned long long int), msg)


private:

    double m_Value;                                                                                                                         // double value
    double m_AbsoluteTolerance;                                                                                                             // absolute tolerance applied to relational operators
    double m_RelativeTolerance;                                                                                                             // relative tolerance applied to relational operators

public:
    
    CDOUBLE(const double p_Value = DEFAULT_INITIAL_DOUBLE_VALUE, const double p_AbsoluteTolerance = FLOAT_TOLERANCE_ABSOLUTE, const double p_RelativeTolerance = FLOAT_TOLERANCE_RELATIVE) {
        m_Value             = p_Value;
        m_AbsoluteTolerance = std::abs(p_AbsoluteTolerance);
        m_RelativeTolerance = std::abs(p_RelativeTolerance);
    }

    // getters
    operator double() const          { return m_Value; }                                                                                    // implicit (default) cast operator (double)
    double Value() const             { return m_Value; }                                                                                    // value

    double AbsoluteTolerance() const { return m_AbsoluteTolerance; }                                                                        // absolute tolerance
    double RelativeTolerance() const { return m_RelativeTolerance; }                                                                        // relative tolerance

    // setters
    void AbsoluteTolerance(const double p_AbsoluteTolerance = FLOAT_TOLERANCE_ABSOLUTE) { m_AbsoluteTolerance = p_AbsoluteTolerance; }      // absolute tolerance
    void RelativeTolerance(const double p_RelativeTolerance = FLOAT_TOLERANCE_RELATIVE) { m_RelativeTolerance = p_RelativeTolerance; }      // relative tolerance
    

    // operator overloads
    // each operator has an overload that takes a CDOUBLE parameter as the other operand, and
    // a catchall overload that takes all other data types as the other operand. 
    // unsupported data types for the other operand will cause the compile to fail - deliberatley

    CDOUBLE& operator=(const CDOUBLE& p_Operand) { if (this != &p_Operand) this->m_Value = p_Operand.Value(); return *this; }               // = CDOUBLE&
    template<typename T> CDOUBLE& operator=(const T p_Operand) { AssertType(T, AssignmentMsg); this->m_Value = p_Operand; return *this; }   // = catchall

    CDOUBLE& operator+=(const CDOUBLE& p_Operand) { if (this != &p_Operand) this->m_Value += p_Operand.Value(); return *this; }             // += CDOUBLE&
    template<typename T> CDOUBLE& operator+=(const T p_Operand) { AssertType(T, AssignmentMsg); this->m_Value += p_Operand; return *this; } // += catchall

    CDOUBLE& operator-=(const CDOUBLE& p_Operand) { if (this != &p_Operand) this->m_Value -= p_Operand.Value(); return *this; }             // -= CDOUBLE&
    template<typename T> CDOUBLE& operator-=(const T p_Operand) { AssertType(T, AssignmentMsg); this->m_Value -= p_Operand; return *this; } // -= catchall

    CDOUBLE& operator*=(const CDOUBLE& p_Operand) { if (this != &p_Operand) this->m_Value *= p_Operand.Value(); return *this; }             // *= CDOUBLE&
    template<typename T> CDOUBLE& operator*=(const T p_Operand) { AssertType(T, AssignmentMsg); this->m_Value *= p_Operand; return *this; } // *= catchall

    CDOUBLE& operator/=(const CDOUBLE& p_Operand) { if (this != &p_Operand) this->m_Value /= p_Operand.Value(); return *this; }             // /= CDOUBLE&
    template<typename T> CDOUBLE& operator/=(const T p_Operand) { AssertType(T, AssignmentMsg); this->m_Value /= p_Operand; return *this; } // /= catchall

    bool operator==(const CDOUBLE& p_Operand) const { return EQUALS(m_Value, p_Operand.Value()); }                                          // == CDOUBLE&
    template<typename T> bool operator==(const T p_Operand) const { AssertType(T, ComparisonMsg); return EQUALS(m_Value, p_Operand); }      // == catchall

    bool operator!=(const CDOUBLE& p_Operand) const { return !EQUALS(m_Value, p_Operand.Value()); }                                         // != CDOUBLE&
    template<typename T> bool operator!=(const T p_Operand) const { AssertType(T, ComparisonMsg); return !EQUALS(m_Value, p_Operand); }     // != catchall

    bool operator<(const CDOUBLE& p_Operand) const { return SMALLER(m_Value, p_Operand.Value()); }                                          // < CDOUBLE&
    template<typename T> bool operator<(const T p_Operand) const { AssertType(T, ComparisonMsg); return SMALLER(m_Value, p_Operand); }      // < catchall

    bool operator<=(const CDOUBLE& p_Operand) const { return !BIGGER(m_Value, p_Operand.Value()); }                                         // <= CDOUBLE&
    template<typename T> bool operator<=(const T p_Operand) const { AssertType(T, ComparisonMsg); return !BIGGER(m_Value, p_Operand); }     // <= catchall
    
    bool operator>(const CDOUBLE& p_Operand) const { return BIGGER(m_Value, p_Operand.Value()); }                                           // > CDOUBLE&
    template<typename T> bool operator>(const T p_Operand) const { AssertType(T, ComparisonMsg); return BIGGER(m_Value, p_Operand); }       // > catchall
    
    bool operator>=(const CDOUBLE& p_Operand) const { return !SMALLER(m_Value, p_Operand.Value()); }                                        // >= CDOUBLE&
    template<typename T> bool operator>=(const T p_Operand) const { AssertType(T, ComparisonMsg); return !SMALLER(m_Value, p_Operand); }    // >= catchall
    
    CDOUBLE& operator+(const CDOUBLE& p_Operand) { if (this != &p_Operand) this->m_Value + p_Operand.Value(); return *this; }               // + CDOUBLE&
    template<typename T> CDOUBLE& operator+(const T p_Operand) { AssertType(T, AssignmentMsg); this->m_Value + p_Operand; return *this; }   // + catchall
    
    CDOUBLE& operator-(const CDOUBLE& p_Operand) { if (this != &p_Operand) this->m_Value - p_Operand.Value(); return *this; }               // - CDOUBLE&
    template<typename T> CDOUBLE& operator-(const T p_Operand) { AssertType(T, AssignmentMsg); this->m_Value - p_Operand; return *this; }   // - catchall
    
    CDOUBLE& operator*(const CDOUBLE& p_Operand) { if (this != &p_Operand) this->m_Value * p_Operand.Value(); return *this; }               // * CDOUBLE&
    template<typename T> CDOUBLE& operator*(const T p_Operand) { AssertType(T, AssignmentMsg); this->m_Value * p_Operand; return *this; }   // * catchall
    
    CDOUBLE& operator/(const CDOUBLE& p_Operand) { if (this != &p_Operand) this->m_Value / p_Operand.Value(); return *this; }               // / CDOUBLE&
    template<typename T> CDOUBLE& operator/(const T p_Operand) { AssertType(T, AssignmentMsg); this->m_Value / p_Operand; return *this; }   // / catchall
        
#undef AssertType
#undef IS
#undef ComparisonMsg
#undef AssignmentMsg
#undef BIGGER
#undef SMALLER
#undef EQUALS
};

#endif // __Cdouble_H__