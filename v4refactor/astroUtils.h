#pragma once

#include "utils.h"
#include "ZDependent.h"

#include <boost/math/distributions.hpp>      // For OStarRotationVelocityFunctor (inverse_gamma + normal)
#include <boost/math/tools/roots.hpp>        // For bracket_and_solve_root in CalculateRotationalVelocityOStar_Ramirez2013

namespace astro {

/*
 * Astrophysics utility functions
 * 
 * The "astro" namespace gathers astrophysics functions that have utility across the
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
 *            (vi) other free functions (defined here, or in utils.h)
 */

GNU_PURE    double       CalculateLifetimeToBGB_Hurley2000(const double p_Mass, const HurleyACoeffsT& p_aCoeffs);
COMPAS_PURE double       CalculateLifetimeToBGB_Hurley2000(const double p_Metallicity, const double p_Mass);
GNU_PURE    double       CalculateMinOmegaForCHE_Butler2018(const double p_Metallicity, const double p_MZAMS);
            double       CalculateRotationalVelocityOStar_Ramirez2013();
GNU_CONST   double       CalculateTemperature(const double p_Luminosity, const double p_Radius);
GNU_CONST   double       CalculateZAMSHAbundance_Pols1998(const double p_Metallicity);
GNU_CONST   double       CalculateZAMSHeAbundance_Pols1998(const double p_Metallicity);
GNU_PURE    ToutLCoeffsT CalculateZAMSLuminosityCoefficients_Tout1996(const double p_Zeta);
COMPAS_PURE double       CalculateZAMSLuminosity_Tout1996(const double p_Metallicity, const double p_MZAMS);
            double       CalculateZAMSOmega_Hurley2000(const double p_MZAMS, const double p_RZAMS);
            double       CalculateZAMSOmega(const double p_MZAMS, const double p_RZAMS);
GNU_PURE    ToutRCoeffsT CalculateZAMSRadiusCoefficients_Tout1996(const double p_Zeta);
COMPAS_PURE double       CalculateZAMSRadius_Tout1996(const double p_Metallicity, const double p_MZAMS);
            double       CalculateZAMSVrot(const double p_MZAMS);
GNU_CONST   double       ConvertPeriodInDaysToSemiMajorAxisInAU(const double p_Mass1, const double p_Mass2, const double p_Period);

GNU_CONST std::tuple<ERROR, double, double> SolveKeplersEquation(const double p_MeanAnomaly, const double p_Eccentricity);


// Template functions
// Must be in header

/*
 * OStarRotationVelocityFunctor
 *
 * @brief
 * Functor for CalculateRotationalVelocityOStar_Ramirez2013() (see below)
 *
 * The analytic cumulative distribution function (CDF) for the equatorial rotational
 * velocity of single O stars is modelled as a mixture of a gamma component and a
 * normal component, per Ramirez-Agudelo et al. 2013, eqs 1-4
 * (https://arxiv.org/abs/1309.2929).
 *
 * 
 * Constructor: initialise the class
 * template <class T> OStarRotationVelocityFunctor(double p_CDF, ERROR *p_Error)
 *
 * @param       p_CDF                           Desired CDF value
 *
 * 
 * Function: calculate the CDF of the O star rotational velocity and compare to desired value
 * T OStarRotationVelocityFunctor(double const& p_Ve)
 *
 * @param       p_Ve                            Rotational velocity, km s^-1
 * @return                                      Difference between star's Roche Lobe radius and radius after mass loss
 */
template <class T>
struct OStarRotationVelocityFunctor {
    OStarRotationVelocityFunctor(double p_CDF) { m_CDF = p_CDF; }
    T operator()(double const& p_Ve) {

        constexpr double alpha  = 4.82;
        constexpr double beta   = 1.0 / 25.0;
        constexpr double mu     = 205.0;
        constexpr double sigma  = 190.0;
        constexpr double iGamma = 0.43;

        boost::math::inverse_gamma_distribution<> gammaComponent(alpha, beta);  // (shape, scale) = (alpha, beta)
        boost::math::normal_distribution<>        normalComponent(mu, sigma);

        // Compute CDF at zero rotational velocity -- the CDF should be relative to this quantity
        const double CDFzero         = (iGamma * boost::math::cdf(gammaComponent, 0.0))  + ((1.0 - iGamma) * boost::math::cdf(normalComponent, 0.0));
        const double CDFunnormalised = (iGamma * boost::math::cdf(gammaComponent, p_Ve)) + ((1.0 - iGamma) * boost::math::cdf(normalComponent, p_Ve));

        return ((CDFunnormalised - CDFzero) / (1.0 - CDFzero)) - m_CDF;
    }
private:
    double m_CDF;
};


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ABUNDANCE                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateZAMSHAbundance_Pols1998
 *
 * @brief
 * Calculate ZAMS hydrogen abundance as a fraction of the star's mass, per Pols et al. 1998
 * 
 *
 * double CalculateZAMSHAbundance_Pols1998(const double p_Metallicity)
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @return                                      ZAMS hydrogen abundance for the star
 */
inline double CalculateZAMSHAbundance_Pols1998(const double p_Metallicity) { return 0.76 - 3.0 * p_Metallicity; }


/*
 * CalculateZAMSHeAbundance_Pols1998
 *
 * @brief
 * Calculate ZAMS helium abundance as a fraction of the star's mass, per Pols et al. 1998
 * 
 *
 * double CalculateZAMSHeAbundance_Pols1998(const double p_Metallicity)
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @return                                      ZAMS helium abundance for the star
 */
inline double CalculateZAMSHeAbundance_Pols1998(const double p_Metallicity) { return 0.24 + 2.0 * p_Metallicity; }


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateLifetimeToBGB_Hurley2000
 *
 * @brief
 * Calculate the lifetime to the Base of the Giant Branch, BGB, (end of the Hertzsprung Gap),
 * per Hurley at al. 2000, eq 4 (plotted in Hurley et al. 2000, fig 5)
 * 
 * For high mass stars, tBGB = tHeI.
 *
 *
 * double CalculateLifetimeToBGB_Hurley2000(const double p_Mass, const HurleyACoeffsT& p_aCoeffs)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_aCoeffs                       Hurley a coefficients
 * @return                                      Lifetime to the Base of the Giant Branch, tBGB (Myr)
 */
inline double CalculateLifetimeToBGB_Hurley2000(const double p_Mass, const HurleyACoeffsT& p_aCoeffs) {

    const auto& a = p_aCoeffs; // Hurley a coefficients

    const double m2 = p_Mass * p_Mass;
    const double m4 = m2 * m2;
    const double m7 = p_Mass * m2 * m4;

    return (a[1] + (a[2] * m4) + (a[3] * std::sqrt(p_Mass) * p_Mass * m4) + m7) / ((a[4] * m2) + (a[5] * m7));
}


/*
 * CalculateLifetimeToBGB_Hurley2000
 *
 * @brief
 * Calculate the lifetime to the Base of the Giant Branch, BGB, (end of the Hertzsprung Gap),
 * per Hurley at al. 2000, eq 4 (plotted in Hurley et al. 2000, fig 5)
 * 
 * For high mass stars, tBGB = tHeI.
 *
 *
 * double CalculateLifetimeToBGB_Hurley2000(const double p_Metallicity, const double p_Mass)
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Lifetime to the Base of the Giant Branch, tBGB (Myr)
 */
inline double CalculateLifetimeToBGB_Hurley2000(const double p_Metallicity, const double p_Mass) {
    return CalculateLifetimeToBGB_Hurley2000(p_Mass, ZDEP->HurleyACoeffs(p_Metallicity));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateZAMSLuminosityCoefficients_Tout1996
 *
 * @brief
 * Calculate ZAMS luminosity coefficients per Tout et al. 1996, table 1
 *
 * Luminosity coefficients depend on the star's metallicity only - so this only needs to be done at most
 * once per star (upon creation), but can also be reused if metallicity doesn't change from one star to
 * the next (e.g. in a population run).
 * 
 *
 * ToutLCoeffsT CalculateZAMSLuminosityCoefficients_Tout1996(const double p_Zeta)
 * 
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      Luminosity coefficients array
 */
inline ToutLCoeffsT CalculateZAMSLuminosityCoefficients_Tout1996(const double p_Zeta) {
    
    // Calculate some powers of zeta - for performance and readability
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;
        
    // Create and initialise luminosity coefficients array - this is the return value
    ToutLCoeffsT lCoeffs;
    
    // Populate luminosity coefficients array.
    // Iterate over luminosity coefficients constants L_COEFF (see constants.h).
    // Each row is indexed by the L_Coeff keys 'ALPHA', 'BETA', 'GAMMA', 'DELTA', 'EPSILON', 'ZETA', 'ETA',
    // and defines the coefficients of the 5 terms (LR_TCoeff coefficients 'a', 'b', 'c', 'd', 'e') 
    for (const auto& [idx, vals] : L_COEFF) {
        using enum LR_TCoeff;
        auto c = [&](LR_TCoeff p) { return vals.at(p); };
        lCoeffs[idx] = c(a) + (c(b) * p_Zeta) + (c(c) * zeta2) + (c(d) * zeta3) + (c(e) * zeta4);
    } 
    
    return lCoeffs;
}  


/*
 * CalculateZAMSLuminosity_Tout1996
 *
 * @brief
 * Calculate the luminosity of a star at ZAMS, per Tout et al. 1996, eq 1
 *
 *
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function nominally relies on the ZAMS values of one or more attributes of the star,
 * and should not be used if the ZAMS values are required but are not known.
 * 
 * The function *could* be called with something other than ZAMS values as parameters - it
 * is up to the caller what values are actually passed. In such cases the ZAMS warning does
 * not apply.
 * 
 * In particular we call this function after a MS merger to determine the luminosity of the
 * merger product - there we reset the effective intial mass of the merger product to the
 * final mass of the merger product, and use that as the ZAMS mass.
 * 
 * 
 * double CalculateZAMSLuminosity_Tout1996(const double p_Metallicity, const double p_MZAMS)
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS luminosity (Lsol)
 */
inline double CalculateZAMSLuminosity_Tout1996(const double p_Metallicity, const double p_MZAMS) {

    // calculate some powers of p_MZAMS - for performance and readability
    // pow() is slow - use multiplication where it makes sense
    const double m0_5 = std::sqrt(p_MZAMS);  // sqrt() is much faster than pow()
    const double m2   = p_MZAMS * p_MZAMS;
    const double m3   = p_MZAMS * m2;
    const double m5   = m2 * m3;
    const double m7   = m2 * m5;
    const double m8   = p_MZAMS * m7;
            
    const ToutLCoeffsT& c = ZDEP->ToutLuminosityCoeffs(p_Metallicity); // get Tout ZAMS luminosity coefficients

    const double top = (c[TOUT_L::ALPHA] * (m5 * m0_5)) + (c[TOUT_L::BETA] * (m3 * m8));    
    const double bot = (c[TOUT_L::GAMMA] + m3) + (c[TOUT_L::DELTA] * m5) + (c[TOUT_L::EPSILON] * m7) + (c[TOUT_L::ZETA] * m8) + (c[TOUT_L::ETA] * (p_MZAMS * m8 * m0_5));

    return top / bot;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       ORBIT                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * ConvertPeriodInDaysToSemiMajorAxisInAU
 *
 * @brief
 * Convert period in days to semi-major axis in AU.
 *
 *
 * double ConvertPeriodInDaysToSemiMajorAxisInAU(const double p_Mass1, const double p_Mass2, const double p_Period)
 *
 * @param       p_Mass1                         Mass of star 1 (Msol)
 * @param       p_Mass2                         Mass of star 2 (Msol)
 * @param       p_Period                        Orbital period (days)
 * @return                                      Semi-major axis (AU)
 */
inline double ConvertPeriodInDaysToSemiMajorAxisInAU(const double p_Mass1, const double p_Mass2, const double p_Period) {
    return std::cbrt((p_Mass1 + p_Mass2) * p_Period * p_Period / DAYS_IN_YEAR / DAYS_IN_YEAR);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateZAMSRadiusCoefficients_Tout1996
 *
 * @brief
 * Calculate ZAMS radius coefficients per Tout et al. 1996, table 2
 *
 * Radius coefficients depend on the star's metallicity only - so this only needs to be done at most
 * once per star (upon creation), but can also be reused if metallicity doesn't change from one star
 * to the next (e.g. in a population run).
 *
 *
 * ToutRCoeffsT CalculateZAMSRadiusCoefficients_Tout1996(const double p_Zeta)
 * 
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      Radius coefficients array
 */
inline ToutRCoeffsT CalculateZAMSRadiusCoefficients_Tout1996(const double p_Zeta) {
    
    // Calculate some powers of zeta - for performance and readability
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;
        
    // Create and inialise radius coefficients array - this is the return value
    ToutRCoeffsT rCoeffs;
    
    // Populate radius coefficients array.
    // Iterate over radius coefficients constants R_COEFF (see constants.h).
    // Each row is indexed by the R_Coeff keys 'THETA', 'IOTA', 'KAPPA', 'LAMBDA', 'MU', 'NU', 'XI', 'OMICRON', 'PI',
    // and defines the coefficients of the 5 terms (LR_TCoeff coefficients 'a', 'b', 'c', 'd', 'e') 
    for (const auto& [idx, vals] : R_COEFF) {
        using enum LR_TCoeff;
        auto c = [&](LR_TCoeff p) { return vals.at(p); };
        rCoeffs[idx] = c(a) + (c(b) * p_Zeta) + (c(c) * zeta2) + (c(d) * zeta3) + (c(e) * zeta4);
    }

    return rCoeffs;
}


/*
 * CalculateZAMSRadius_Tout1996
 *
 * @brief
 * Calculate the radius of a star at ZAMS, per Tout et al. 1996, eq 2
 *
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function nominally relies on the ZAMS values of one or more attributes of the star,
 * and should not be used if the ZAMS values are required but are not known.
 * 
 * The function *could* be called with something other than ZAMS values as parameters - it
 * is up to the caller what values are actually passed. In such cases the ZAMS warning does
 * not apply.
 * 
 * 
 * double CalculateZAMSRadius_Tout1996(const double p_Metallicity, const double p_MZAMS)
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS radius of the star (Rsol)
 */
inline double CalculateZAMSRadius_Tout1996(const double p_Metallicity, const double p_MZAMS) {

    // calculate some powers of p_MZAMS - for performance and readability
    const double m0_5  = std::sqrt(p_MZAMS);
    const double m2    = p_MZAMS * p_MZAMS;
    const double m6    = m2 * m2 * m2;
    const double m8    = m6 * m2;
    const double m11   = p_MZAMS * m8 * m2;
    const double m19   = m11 * m8;
    const double m19_5 = m19 * m0_5;
    
    const ToutRCoeffsT& c = ZDEP->ToutRadiusCoeffs(p_Metallicity); // Get Tout ZAMS radius coefficients

    const double top = (c[TOUT_R::THETA] * (m2 * m0_5)) + (c[TOUT_R::IOTA] * m6 * m0_5) + (c[TOUT_R::KAPPA] * m11) + (c[TOUT_R::LAMBDA] * m19) + (c[TOUT_R::MU] * m19_5);
    const double bot = c[TOUT_R::NU] + (c[TOUT_R::XI] * m2) + (c[TOUT_R::OMICRON] * (m8 * m0_5)) + (m6 * m6 * m6 * m0_5) + (c[TOUT_R::PI] * m19_5);
            
    return top / bot;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ROTATION                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateMinOmegaForCHE_Butler2018
 *
 * @brief
 * Calculate the minimum angular frequency (in rad yr^-1) at which CHE will occur for a star,
 * given the metallicity and ZAMS mass of the star.
 *
 * Mandel's fit from Butler 2018 (see Riley et al. 2021, appendix A
 * (https://doi.org/10.1093/mnras/stab1291))
 *
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function nominally relies on the ZAMS values of one or more attributes of the star,
 * and should not be used if the ZAMS values are required but are not known.
 * 
 * The function *could* be called with something other than ZAMS values as parameters - it
 * is up to the caller what values are actually passed. In such cases the ZAMS warning does
 * not apply.
 * 
 * 
 * double CalculateMinOmegaForCHE_Butler2018(const double p_Metallicity, const double p_MZAMS)
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      Minimum angular frequency for CHE (rad yr^-1)
 */
inline double CalculateMinOmegaForCHE_Butler2018(const double p_Metallicity, const double p_MZAMS) {

    const double mRatio = p_MZAMS;      // In Msol, so ratio is just p_MZAMS

    // Calculate omegaCHE(M, Z = 0.004).
    double omegaZ004 = 0.0;
    if (p_MZAMS <= MANDEL_BUTLER_CHE_MASS_BREAK) {
        for (SizeT i = 0; i < CHE_COEFFICIENTS.size(); i++) {
            omegaZ004 += CHE_COEFFICIENTS[i] * utils::IntPow(mRatio, i) / PPOW(mRatio, 0.4);
        }
    }
    else {
        for (SizeT i = 0; i < CHE_COEFFICIENTS.size(); i++) {
            omegaZ004 += CHE_COEFFICIENTS[i] * utils::IntPow(MANDEL_BUTLER_CHE_MASS_BREAK, i) / PPOW(mRatio, 0.4);
        }
    }

    // Calculate omegaCHE(M, Z).
    return (1.0 / ((0.09 * log(p_Metallicity / 0.004)) + 1.0) * omegaZ004) * SECONDS_IN_YEAR;
}


/*
 * CalculateZAMSOmega_Hurley2000
 *
 * @brief
 * Calculate the ZAMS angular frequency of a star, per Hurley et al. 2000, eq 108
 * 
 *
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function nominally relies on the ZAMS values of one or more attributes of the star,
 * and should not be used if the ZAMS values are required but are not known.
 * 
 * The function *could* be called with something other than ZAMS values as parameters - it
 * is up to the caller what values are actually passed. In such cases the ZAMS warning does
 * not apply.
 * 
 * 
 * double CalculateZAMSOmega_Hurley2000(const double p_MZAMS, const double p_RZAMS)
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_vRotDist                      Rotational velocity distribution
 * @return                                      ZAMS angular frequency (rad yr^-1)
 */
inline double CalculateZAMSOmega_Hurley2000(const double p_MZAMS, const double p_RZAMS) {
    double vRot = CalculateZAMSVrot(p_MZAMS);
    return 45.35 * vRot / p_RZAMS;
}


/*
 * CalculateZAMSOmega
 *
 * @brief
 * Calculate the angular frequency of a star at ZAMS.
 * 
 * Calls relevant angular frequency function based on the evolutionary mode.
 * 
 *
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function nominally relies on the ZAMS values of one or more attributes of the star,
 * and should not be used if the ZAMS values are required but are not known.
 * 
 * The function *could* be called with something other than ZAMS values as parameters - it
 * is up to the caller what values are actually passed. In such cases the ZAMS warning does
 * not apply.
 * 
 * 
 * double CalculateZAMSOmega(const double p_MZAMS, const double p_RZAMS)
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @return                                      ZAMS angular frequency (rad yr^-1)
 */
inline double CalculateZAMSOmega(const double p_MZAMS, const double p_RZAMS) {

    double omega;

    switch (OPTIONS->Mode()) {                                              // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                                    // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                                    // HURLEY BSE
            omega = CalculateZAMSOmega_Hurley2000(p_MZAMS, p_RZAMS);
            break;
        
        default:                                                            // Unexpected mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNEXPECTED_EVOLUTION_MODE);           // Throw error
    }       

    return omega;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    TEMPERATURE                                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateTemperature
 *
 * @brief
 * Calculate the effective temperature of the star, given the luminosity and radius
 * of the star, using the Stefan–Boltzmann law.
 *
 *
 * double CalculateTemperature(const double p_Luminosity, const double p_Radius)
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @return                                      Effective temperature of the star (Tsol)
 */
inline double CalculateTemperature(const double p_Luminosity, const double p_Radius) {
    return std::sqrt(std::sqrt(p_Luminosity)) / std::sqrt(p_Radius);
}

} // namespace astro
