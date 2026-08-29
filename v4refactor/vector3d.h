#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <utility>
#include <ostream>
#include <algorithm>

#include "constants.h"
#include "Errors.h"
#include "Cdouble.h"


class Vector3d {

private:

    double m_x = 0.0;
    double m_y = 0.0;
    double m_z = 0.0;

public:

    // Constructors
    Vector3d() = default;
    Vector3d(const double p_x, const double p_y, const double p_z) : m_x(p_x), m_y(p_y), m_z(p_z) {}
    Vector3d(const Vector3dT& p_Vec) : m_x(p_Vec[0]), m_y(p_Vec[1]), m_z(p_Vec[2]) {} // Note: implicit converting constructor


    // Getters 
    GNU_PURE DblVectorT asDblVector() const { return { m_x, m_y , m_z }; }
    GNU_PURE double     Magnitude() const   { return std::sqrt(Dot(*this, *this)); }
    GNU_PURE Vector3d   UnitVector() const  { return *this / Magnitude(); } 
    GNU_PURE double     xValue() const      { return m_x; }
    GNU_PURE double     yValue() const      { return m_y; }
    GNU_PURE double     zValue() const      { return m_z; }


    // Operator overloads
    const double& operator[](const SizeT p_i) const {
             if (p_i == 0) return m_x;
        else if (p_i == 1) return m_y;
        else if (p_i == 2) return m_z;
        else THROW_ERROR_STATIC(ERROR::INDEX_OUT_OF_RANGE);
    }
    double& operator[](const SizeT p_i) { return const_cast<double&>(std::as_const(*this)[p_i]); }

    Vector3d& operator+= (const Vector3d& p_Vec) {
        m_x += p_Vec.m_x;
        m_y += p_Vec.m_y;
        m_z += p_Vec.m_z;
        return *this;
    }
    friend Vector3d operator+ (Vector3d p_Lhs, const Vector3d& p_Rhs) {
        p_Lhs += p_Rhs;
        return p_Lhs;
    }

    Vector3d& operator-= (const Vector3d& p_Vec) {
        m_x -= p_Vec.m_x;
        m_y -= p_Vec.m_y;
        m_z -= p_Vec.m_z;
        return *this;
    }  
    friend Vector3d operator- (Vector3d p_Lhs, const Vector3d& p_Rhs) {
        p_Lhs -= p_Rhs;
        return p_Lhs;
    }
    
    Vector3d operator-() const { return Vector3d(-m_x, -m_y, -m_z); }

    Vector3d& operator*= (const double p_Scalar) {
        m_x *= p_Scalar;
        m_y *= p_Scalar;
        m_z *= p_Scalar;
        return *this;
    }
    friend Vector3d operator* (const double p_Scalar, Vector3d p_Vec) {
        p_Vec *= p_Scalar;
        return p_Vec;
    }   
    friend Vector3d operator* (Vector3d p_Vec, const double p_Scalar) {
        p_Vec *= p_Scalar;
        return p_Vec;
    }

    Vector3d& operator/= (const double p_Scalar) {
        m_x /= p_Scalar;
        m_y /= p_Scalar;
        m_z /= p_Scalar;
        return *this;
    }
    friend Vector3d operator/ (Vector3d p_Vec, const double p_Scalar) {
        p_Vec /= p_Scalar;
        return p_Vec;
    }

    friend bool operator== (const Vector3d& p_Lhs, const Vector3d& p_Rhs) {
        // We compare each component by wrapping the LHS component in CDOUBLE - this makes
        // the compares tolerance-aware by using the CDOUBLE::operator== (RHS promotes implicitly).
        // NOTE: this is non-transitive (see CDOUBLE documentation)
        return CDOUBLE(p_Lhs.m_x) == p_Rhs.m_x && CDOUBLE(p_Lhs.m_y) == p_Rhs.m_y && CDOUBLE(p_Lhs.m_z) == p_Rhs.m_z;
    }
    friend bool operator!= (const Vector3d& p_Lhs, const Vector3d& p_Rhs) { return !(p_Lhs == p_Rhs); }

    friend std::ostream &operator << (std::ostream &p_os, const Vector3d& p_Vec) { return p_os << "{" << p_Vec[0] << ", " << p_Vec[1] << ", " << p_Vec[2] << "}"; }


    // Member functions 
    GNU_PURE static double   AngleBetween(const Vector3d& p_Vec1, const Vector3d& p_Vec2);
    GNU_PURE        Vector3d ChangeBasis(const double p_ThetaE, const double p_PhiE, const double p_PsiE) const;
    GNU_PURE static Vector3d Cross(const Vector3d& p_Vec1, const Vector3d& p_Vec2);
    GNU_PURE static double   Dot(const Vector3d& p_Vec1, const Vector3d& p_Vec2);
    GNU_PURE static Vector3d MatrixMult(const Matrix3x3T& p_Matrix, const Vector3d& p_Vector);
    GNU_PURE        Vector3d RotateVectorAboutX(const double p_Theta) const;
    GNU_PURE        Vector3d RotateVectorAboutY(const double p_Theta) const;
    GNU_PURE        Vector3d RotateVectorAboutZ(const double p_Theta) const;

};


/*
 * AngleBetween
 *
 * @brief
 * Calculate the angle between two vectors.
 * 
 *
 * static double AngleBetween(const Vector3d& p_Vec1, const Vector3d& p_Vec2)
 *
 * @param       p_Vec1                          Vector 1
 * @param       p_Vec2                          Vector 2
 * @return                                      Angle between vector 1 and vector 2 (rad)
 */
inline double Vector3d::AngleBetween(const Vector3d& p_Vec1, const Vector3d& p_Vec2) {
    return std::acos(std::clamp(Dot(p_Vec1, p_Vec2) / (p_Vec1.Magnitude() * p_Vec2.Magnitude()), -1.0, 1.0));
}


/*
 * Cross
 *
 * @brief
 * Calculate the cross product of two vectors.
 * 
 *
 * static Vector3d Cross(const Vector3d& p_Vec1, const Vector3d& p_Vec2)
 *
 * @param       p_Vec1                          Vector 1
 * @param       p_Vec2                          Vector 2
 * @return                                      Cross product of vector 1 and vector 2
 */
inline Vector3d Vector3d::Cross(const Vector3d& p_Vec1, const Vector3d& p_Vec2) {
    return Vector3d(
        (p_Vec1.m_y * p_Vec2.m_z) - (p_Vec1.m_z * p_Vec2.m_y),
        (p_Vec1.m_z * p_Vec2.m_x) - (p_Vec1.m_x * p_Vec2.m_z),
        (p_Vec1.m_x * p_Vec2.m_y) - (p_Vec1.m_y * p_Vec2.m_x)
    );
}


/*
 * Dot
 *
 * @brief
 * Calculate the dot product of two vectors.
 * 
 *
 * static double Dot(const Vector3d& p_Vec1, const Vector3d& p_Vec2)
 *
 * @param       p_Vec1                          Vector 1
 * @param       p_Vec2                          Vector 2
 * @return                                      Dot product of vector 1 and vector 2
 */
inline double Vector3d::Dot(const Vector3d& p_Vec1, const Vector3d& p_Vec2) {
    return (p_Vec1.m_x * p_Vec2.m_x) + (p_Vec1.m_y * p_Vec2.m_y) + (p_Vec1.m_z * p_Vec2.m_z);
}


/*
 * ChangeBasis
 *
 * @brief
 * Redefine a vector from one coordinate basis to another using Euler angles.
 * 
 * For a vector defined in a new coordinate basis, (X',Y',Z'), we want to 
 * find it's values in a previous coordinate basis (X,Y,Z) so that we can 
 * add the new vector to vectors defined in the previous basis.
 *
 * For a change of basis from (X,Y,Z)->(X',Y',Z')
 * 
 *     ThetaE [0, pi] is the angle between Z and Z',
 *     Vector N := Z x Z' (cross product)
 *     PhiE [0, 2pi) is the angle between X and N
 *     PsiE [0, 2pi) is the angle between X' and N
 * 
 * These angles uniquely determine the change of basis, which is applied 
 * in the form of a rotation matrix R as a function of these angles.
 * 
 * V_{X,Y,Z} = R * V_{X',Y',Z'} gives the vector in the original coordinates.
 *
 * For details, see:
 * https://en.wikipedia.org/wiki/Euler_angles
 * https://en.wikipedia.org/wiki/Change_of_basis
 *
 * 
 * Vector3d::ChangeBasis(const double p_ThetaE, const double p_PhiE, const double p_PsiE) const
 * 
 * @param       p_ThetaE                        Euler angle Theta (rad) 
 * @param       p_PhiE                          Euler angle Phi   (rad) 
 * @param       p_PsiE                          Euler angle Psi   (rad) 
 * @return                                      Vector in specified basis
 */
GNU_PURE inline Vector3d Vector3d::ChangeBasis(const double p_ThetaE, const double p_PhiE, const double p_PsiE) const {

    // Evaluate the trig functions once
    const double cTheta = std::cos(p_ThetaE);
    const double cPhi   = std::cos(p_PhiE);
    const double cPsi   = std::cos(p_PsiE);
    const double sTheta = std::sin(p_ThetaE);
    const double sPhi   = std::sin(p_PhiE);
    const double sPsi   = std::sin(p_PsiE);

    // Rotation matrix - fixed 3x3, stack-allocated (no heap allocations).    
    const Matrix3x3T R = {{
        { cPhi * cPsi - sPhi * cTheta * sPsi, -cPhi * sPsi - sPhi * cTheta * cPsi,  sTheta * sPhi },
        { sPhi * cPsi + cPhi * cTheta * sPsi, -sPhi * sPsi + cPhi * cTheta * cPsi, -sTheta * cPhi },
        {                      sTheta * sPsi,                       sTheta * cPsi,  cTheta        }
    }};

    return MatrixMult(R, *this); // Apply rotation: newVector = R * (*this)
}


/*
 * MatrixMult
 *
 * @brief
 * Right-multiply a 3x3 matrix by a 3-vector.
 * 
 *
 * static Vector3d MatrixMult(const Matrix3x3T& p_Matrix, const Vector3d& p_Vec)
 *
 * @param       p_Matrix                        Matrix to be multiplied
 * @param       p_Vec                           Vector to multiply
 * @return                                      Matrix-vector product (M x v)
 */
inline Vector3d Vector3d::MatrixMult(const Matrix3x3T& p_Matrix, const Vector3d& p_Vec) {

    Vector3d result(0.0, 0.0, 0.0);

    for (SizeT row = 0; row < 3; row++)
        for (SizeT col = 0; col < 3; col++)
            result[row] += p_Matrix[row][col] * p_Vec[col];

    return result;
}


/*
 * RotateVectorAboutX
 *
 * @brief
 * Rotate a vector about the X axis.
 *
 * 
 * Vector3d RotateVectorAboutX(const double p_Theta) const
 *
 * @param       p_Theta                         Rotation angle (rad) 
 * @return                                      Rotated vector
 */
GNU_PURE inline Vector3d Vector3d::RotateVectorAboutX(const double p_Theta) const {

    // Evaluate the trig functions once
    const double cTheta = std::cos(p_Theta);
    const double sTheta = std::sin(p_Theta);

    const Matrix3x3T R = {{         // Define the rotation matrix
        { 1.0,  0.0,     0.0    },
        { 0.0,  cTheta, -sTheta },
        { 0.0,  sTheta,  cTheta }
    }};

    return MatrixMult(R, *this);    // Apply rotation and return result
}


/*
 * RotateVectorAboutY
 *
 * @brief
 * Rotate a vector about the Y axis.
 *
 * 
 * Vector3d RotateVectorAboutY( const double p_Theta) const
 *
 * @param       p_Theta                         Rotation angle (rad) 
 * @return                                      Rotated vector
 */
GNU_PURE inline Vector3d Vector3d::RotateVectorAboutY(const double p_Theta) const {

    // Evaluate the trig functions once
    const double cTheta = std::cos(p_Theta);
    const double sTheta = std::sin(p_Theta);

    const Matrix3x3T R = {{         // Define the rotation matrix
        {  cTheta, 0.0, sTheta },
        {  0.0,    1.0, 0.0    },
        { -sTheta, 0.0, cTheta }
    }};

    return MatrixMult(R, *this);    // Apply rotation and return result
}


/*
 * RotateVectorAboutZ
 *
 * @brief
 * Rotate a vector about the Z axis.
 *
 * 
 * Vector3d RotateVectorAboutZ( const double p_Theta) const
 *
 * @param       p_Theta                         Rotation angle (rad) 
 * @return                                      Rotated vector
 */
GNU_PURE inline Vector3d Vector3d::RotateVectorAboutZ(const double p_Theta) const {

    // Evaluate the trig functions once
    const double cTheta = std::cos(p_Theta);
    const double sTheta = std::sin(p_Theta);

    const Matrix3x3T R = {{         // Define the rotation matrix
        { cTheta, -sTheta,  0.0 },
        { sTheta,  cTheta,  0.0 },
        { 0.0,     0.0,     1.0 }
    }};

    return MatrixMult(R, *this);    // Apply rotation and return result
}
