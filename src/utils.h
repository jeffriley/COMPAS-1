#pragma once

#include <memory>
#include <iostream>
#include <stdarg.h>
#include <algorithm>
#include <cstring>
#include <iomanip>
#include <execinfo.h>
#include <cxxabi.h>
#include <type_traits>

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "changelog.h"
#include "hdf5.h"

#include "Rand.h"
#include "ErrorCatalog.h"

#include <boost/version.hpp>
#include <boost/filesystem.hpp>


/*
 * Utility functions
 * 
 * The "utils" namespace gathers functions that have utility across the COMPAS code.
 * 
 * These are free functions that:
 * 
 *    (a) don't really belong to any of the stellar or binary classes,
 *    (b) aren't required to run against an instantiated object, and
 *    (c) rely only on:
 *             (i) parameters passed, and/or
 *            (ii) constants.h, and/or
 *           (iii) the GLOBALS object, and/or
 *            (iv) the OPTIONS object, and/or
 *             (v) the RAND object, and/or
 *             (v) library functions, and/or
 *            (vi) other free functions (defined here)
 * 
 * This is a header file implementation - no .cpp file
 */

namespace utils {


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       ADMIN                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * GetBoostVersion
 *
 * @brief
 * Get Boost version string in format MM.mm.rr
 *
 * 
 * StrT GetBoostVersion()
 * 
 * @return                                    String containing Boost version in format MM.mm.rr
 */
GNU_CONST inline StrT GetBoostVersion() { 
    return std::to_string(BOOST_VERSION / 100000) + "." + std::to_string(BOOST_VERSION / 100 % 1000) + "." + std::to_string(BOOST_VERSION % 100);
}


/*
 * GetGslVersion
 *
 * @brief
 * Get gsl version string in format MM.mm.rr
 *
 * 
 * StrT GetGslVersion()
 * 
 * @return                                    String containing gsl version in format MM.mm.rr
 *                                            Will be "Not available" if not able to retrieve the actual value.
 */
GNU_CONST StrT GetGslVersion() {

    StrT versionStr = "Not available";                                                                  // Default return value

    char buffer[128];                                                                                   // Command return buffer
    std::unique_ptr<FILE, int(*)(FILE*)> pipe(popen("gsl-config --version", "r"), &pclose);             // Open pipe for command
    if (pipe) {                                                                                         // Ok?
        versionStr = "";                                                                                // Yes
        while (fgets(buffer, sizeof(buffer), pipe.get()) != nullptr) versionStr += buffer;              // Copy buffer
        if (!versionStr.empty() && versionStr[versionStr.length() - 1] == '\n') versionStr.pop_back();  // Remove trailing newline if necessary
    }

    return versionStr;
}


/*
 * GetHDF5Version
 *
 * @brief
 * Get HDF5 library version string in format MM.mm.rr
 *
 * 
 * StrT GetHDF5Version()
 * 
 * @return                                    String containing HDF5 library version in format MM.mm.rr
 *                                            Will be "Not available" if not able to retrieve the actual value.
 */
GNU_CONST StrT GetHDF5Version() {

    StrT versionStr = "Not available";                                                                              // Default return value

    unsigned majorNum, minorNum, releaseNum;
    herr_t status = H5get_libversion(&majorNum, &minorNum, &releaseNum);                                            // Retrieve HDF5 library version
    if (status >= 0)                                                                                                // Ok?
        versionStr = std::to_string(majorNum) + "." + std::to_string(minorNum) + "." + std::to_string(releaseNum);  // Yes - set version string

    return versionStr;
}


/*
 * SplashScreen
 *
 * @brief 
 * Construct COMPAS splash string.
 * Prints splash string to stdout if required.
 *
 *
 * StrT SplashScreen(const bool p_Print)
 * 
 * @param       p_Print                         Boolean indicating whether splash string should be printed.  Default is TRUE.
 * @return                                      Splash string
 */
GNU_CONST StrT SplashScreen(const bool p_Print) {

    // Set compiler name
    StrT compiler = COMPILER_ID == GCC_COMPILER ? "gcc" : (COMPILER_ID == CLANG_COMPILER ? "Clang" : (COMPILER_ID == MSC_COMPILER ? "MSCV" : "?"));

    // Construct compiler details string
    StrT major = COMPILER_MAJOR < 0 ? "?" : std::to_string(COMPILER_MAJOR);
    StrT minor = COMPILER_MINOR < 0 ? "?" : std::to_string(COMPILER_MINOR);
    StrT patch = COMPILER_PATCH < 0 ? "?" : std::to_string(COMPILER_PATCH);

    StrT compilerStr = compilerName + " " + major + "." + minor + "." + patch;

    // Construct the splash string
    StrT splashString = "\nCOMPAS v" + 
                        VERSION_STRING + 
                        " (" + compilerStr +
                        ", gsl v" + utils::GetGslVersion() + 
                        ", boost v" + utils::GetBoostVersion() + 
                        ", HDF5 v" + utils::GetHDF5Version() + ")" +
                        "\nCompact Object Mergers: Population Astrophysics and Statistics"
                        "\nby Team COMPAS (http://compas.science/index.html)"
                        "\nA binary star simulator\n"
                        "\nGo to https://compas.readthedocs.io/en/latest/index.html for the online documentation"
                        "\nCheck https://compas.readthedocs.io/en/latest/pages/whats-new.html to see what's new in the latest release\n";

    if (p_Print) std::cout << splashString << std::endl;    // Print the splash string if required

    return splashString;                                    // Return the splash string
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              ASTRO / PHYSICS / MATH                               //
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
 * double CalculateLifetimeToBGB_Hurley2000(const double p_Mass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Lifetime to the Base of the Giant Branch, tBGB (Myr)
 */
COMPAS_PURE inline double CalculateLifetimeToBGB_Hurley2000(const double p_Mass) {

    const DblVectorT a = GLOBALS->HurleyAcoefficients(); // get Hurley a coefficients

    const double m2 = p_Mass * p_Mass;
    const double m4 = m2 * m2;
    const double m7 = p_Mass * m2 * m4;

    return (a[1] + (a[2] * m4) + (a[3] * std::sqrt(p_Mass) * p_Mass * m4) + m7) / ((a[4] * m2) + (a[5] * m7));
}


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
 * double CalculateMinOmegaForCHE_Butler2018(const double p_Z, const double p_MZAMS)
 *
 * @param       p_Z                             Metallicity of the star
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      Minimum angular frequency for CHE (rad yr^-1)
 */
GNU_PURE double CalculateMinOmegaForCHE_Butler2018(const double p_Z, const double p_MZAMS) {

    const double mRatio = p_MZAMS;      // In Msol, so ratio is just p_MZAMS

    // Calculate omegaCHE(M, Z = 0.004).
    double omegaZ004 = 0.0;
    if (p_MZAMS <= MANDEL_BUTLER_CHE_MASS_BREAK) {
        for (SizeT i = 0; i < CHE_Coefficients.size(); i++) {
            omegaZ004 += CHE_Coefficients[i] * utils::intPow(mRatio, i) / PPOW(mRatio, 0.4);
        }
    }
    else {
        for (SizeT i = 0; i < CHE_Coefficients.size(); i++) {
            omegaZ004 += CHE_Coefficients[i] * utils::intPow(MANDEL_BUTLER_CHE_MASS_BREAK, i) / PPOW(mRatio, 0.4);
        }
    }

    // Calculate omegaCHE(M, Z).
    return (1.0 / ((0.09 * log(p_Z / 0.004)) + 1.0) * omegaZ004) * SECONDS_IN_YEAR;
}


/*
 * CalculateZAMSOmega
 *
 * @brief
 * Calculate the angular frequency of a star at ZAMS.
 * 
 * Calls relevant angular frequency function based on the evolutionary mode given in program options.
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
COMPAS_PURE double CalculateZAMSOmega(const double p_MZAMS, const double p_RZAMS) {

    double omega;

    Switch (OPTIONS->Mode()) {                                              // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                         // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                         // HURLEY BSE
            omega = CalculateZAMSOmega_Hurley2000(p_MZAMS, p_RZAMS);
            break;
        
        default:                                                            // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                     // throw error
    }       

    return omega;
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
 * @return                                      ZAMS angular frequency (rad yr^-1)
 */
COMPAS_PURE inline double CalculateZAMSOmega_Hurley2000(const double p_MZAMS, const double p_RZAMS) {
    double vRot = utils::CalculateZAMSVrot(p_MZAMS);
    return 45.35 * vRot / p_RZAMS;
}


/*
 * CalculateZAMSHAbundance_Pols1998
 *
 * @brief
 * Calculate ZAMS hydrogen abundance as a fraction of the star's mass, per Pols et al. 1998
 * 
 *
 * double CalculateZAMSHAbundance_Pols1998(const double p_Z)
 *
 * @param       p_Z                             Metallicity of the star
 * @return                                      ZAMS hydrogen abundance for the star
 */
GNU_CONST inline double CalculateZAMSHAbundance_Pols1998(const double p_Z) { return 0.76 - 3.0 * p_Z; }


/*
 * CalculateZAMSHeAbundance_Pols1998
 *
 * @brief
 * Calculate ZAMS helium abundance as a fraction of the star's mass, per Pols et al. 1998
 * 
 *
 * double CalculateZAMSHeAbundance_Pols1998(const double p_Z)
 * 
 * @param       p_Z                             Metallicity of the star
 * @return                                      ZAMS helium abundance for the star
 */
GNU_CONST inline double Globals::CalculateZAMSHeAbundance_Pols1998(const double p_Z) { return 0.24 + 2.0 * p_Z; }


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
 * DblVectorT CalculateZAMSLuminosityCoefficients_Tout1996(const double p_Zeta)
 * 
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      Luminosity coefficients vector
 */
GNU_PURE DblVectorT CalculateZAMSLuminosityCoefficients_Tout1996(const double p_Zeta) {
        
    // Create and inialise luminosity coefficients vector - this is the return value
    DbltVectorT lCoeffs(TOUT_L_COEFF.size(), DEFAULT_INITIAL_DOUBLE_VALUE);
    
    // Calculate some powers of zeta - for performance and readability
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;
    
    // Populate luminosity coefficients vector.
    // Iterate over luminosity coefficients constants TOUT_L_COEFF (see constants.h).
    // Each row is indexed by the TOUT_L_Coeff keys 'ALPHA', 'BETA', 'GAMMA', 'DELTA', 'EPSILON', 'ZETA', 'ETA',
    // and defines the coefficients of the 5 terms (TOUT_LR_TCoeff coefficients 'a', 'b', 'c', 'd', 'e') 
    for (const auto& [idx, vals] : TOUT_L_COEFF) {
        using enum TOUT_LR_TCoeff;
        auto coeff   = [&](TOUT_LR_TCoeff p) { return vals[static_cast<SizeT>(p)]; };
        lCoeffs[idx] = coeff(a) + (coeff(b) * p_Zeta) + (coeff(c) * zeta2) + (coeff(d) * zeta3) + (coeff(e) * zeta4);
    } 
    
    // Return the luminosity coefficients vector by value - NRVO takes care of performance/efficiency
    return lCoeffs;
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
 * double CalculateZAMSRadius_Tout1996(const double p_MZAMS) const
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS radius of the star (Rsol)
 */
COMPAS_PURE inline double CalculateZAMSRadius_Tout1996(const double p_MZAMS) const {

    // calculate some powers of p_MZAMS - for performance and readability
    const double m0_5  = std::sqrt(p_MZAMS);
    const double m2    = p_MZAMS * p_MZAMS;
    const double m6    = m2 * m2 * m2;
    const double m8    = m6 * m2;
    const double m11   = p_MZAMS * m8 * m2;
    const double m19   = m11 * m8;
    const double m19_5 = m19 * m0_5;
    
    c = GLOBALS->ToutZAMSRadiusCoefficients(); // get Tout ZAMS radius coefficients

    using enum TOUT_R_Coeff;

    const double top = (c[THETA] * (m2 * m0_5)) + (c[IOTA] * m6 * m0_5) + (c[KAPPA] * m11) + (c[LAMBDA] * m19) + (c[MU] * m19_5);
            
    return top / (c[NU] + (c[XI] * m2) + (c[OMICRON] * (m8 * m0_5)) + (m6 * m6 * m6 * m0_5) + (c[PI] * m19_5));
}


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
 * DBL_VECTOR CalculateZAMSRadiusCoefficients_Tout1996(const double p_Zeta)
 * 
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      Radius coefficients vector
 */
GNU_PURE DBL_VECTOR Globals::CalculateZAMSRadiusCoefficients_Tout1996(const double p_Zeta) {
        
    // Create and inialise radius coefficients vector - this is the return value
    DBL_VECTOR rCoeffs(TOUT_R_COEFF.size(), DEFAULT_INITIAL_DOUBLE_VALUE);
    
    // Calculate some powers of zeta - for performance and readability
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;
    
    // Populate radius coefficients vector.
    // Iterate over radius coefficients constants TOUT_R_COEFF (see constants.h).
    // Each row is indexed by the TOUT_R_Coeff keys 'THETA', 'IOTA', 'KAPPA', 'LAMBDA', 'MU', 'NU', 'XI', 'OMICRON', 'PI',
    // and defines the coefficients of the 5 terms (TOUT_LR_TCoeff coefficients 'a', 'b', 'c', 'd', 'e') 
    for (const auto& [idx, vals] : TOUT_R_COEFF) {
        using enum TOUT_LR_TCoeff;
        auto coeff   = [&](TOUT_LR_TCoeff p) { return vals[static_cast<SizeT>(p)]; };
        lCoeffs[idx] = coeff(a) + (coeff(b) * p_Zeta) + (coeff(c) * zeta2) + (coeff(d) * zeta3) + (coeff(e) * zeta4);
    }

    // Return the radius coefficients vector by value - NRVO takes care of performance/efficiency
    return rCoeffs;
}


/*
 * CalculateZAMSVrot
 *
 * @brief
 * Calculate the ZAMS equatorial rotational velocity of a star, given the ZAMS mass.
 *
 * How the rotational velocity is calculated depends on the distribution specified by the
 * `--rotational-velocity-distribution` program option (passed as a paremeter to this function).
 * 
 * Possible values of the rotational velocity distribution are:
 * 
 *    - ROTATIONAL_VELOCITY_DISTRIBUTION::HURLEY
 *         - rotational velocity calculated using Hurley et al. 2000, eq 107
 * 
 *    - ROTATIONAL_VELOCITY_DISTRIBUTION::VLTFLAMES
 *         - rotational velocity calculated using:
 * 
 *            - for O-stars:
 *                 - single stars: Ramirez-Agudelo et al. 2013 (https://arxiv.org/abs/1309.2929)
 *                 - spectroscopic binaries: Ramirez-Agudelo et al. 2015 (https://arxiv.org/abs/1507.02286)
 * 
 *            - for B-stars: Dufton et al. 2013 (https://arxiv.org/abs/1212.2424)
 * 
 *            - for lower mass stars: Hurley et al. 2000, eq 107 (for now)
 *
 *    - ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO
 *         - no rotation; rotational velocity is 0.0 
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
 * double CalculateZAMSVrot(double p_MZAMS)
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS equatorial rotational velocity (km s^-1)
 */
COMPAS_PURE double CalculateZAMSVrot(const double p_MZAMS) {

    double vRot = 0.0;

    switch (OPTIONS->RotationalVelocityDistribution()) {                // which rotational velocity distribution?

        case ROTATIONAL_VELOCITY_DISTRIBUTION::HURLEY:                  // HURLEY

            // Hurley et al. 2000, eq 107 (uses fit from Lang 1992)
            vRot = (330.0 * PPOW(p_MZAMS, 3.3)) / (15.0 + PPOW(p_MZAMS, 3.45));
            break;

         case ROTATIONAL_VELOCITY_DISTRIBUTION::VLTFLAMES:              // VLTFLAMES

            // Rotational velocity based on VLT-FLAMES survey.
            // For O-stars (taken to be above 16 Msol), use results
            // of Ramirez-Agudelo et al. (2013) https://arxiv.org/abs/1309.2929 (single stars)
            // and Ramirez-Agudelo et al. (2015) https://arxiv.org/abs/1507.02286 (spectroscopic binaries)
            // For B-stars (taken to be between 2 and 16 Msol) use results
            // of Dufton et al. (2013) https://arxiv.org/abs/1212.2424
            // For lower mass stars, default back to  Hurley et al. 2000 distribution for now

            if (utils::Compare(p_MZAMS, 16.0) >= 0) {
                vRot = CalculateRotationalVelocityOStar_Ramirez2013();
                vRot = std::max(vRot, 0.0);                             // Set to no rotation if no positive solution found; warning already raised
            }
            else if (utils::Compare(p_MZAMS, 2.0) >= 0) {
                vRot = utils::InverseSampleFromTabulatedCDF(RAND->Random(), BStarRotationalVelocityCDFTable);
            }
            else {
                // Don't know what better to use for low mass stars so for now
                // default to Hurley et al. 2000, eq 107 (uses fit from Lang 1992)
                vRot = (330.0 * PPOW(p_MZAMS, 3.3)) / (15.0 + PPOW(p_MZAMS, 3.45));
            }
            break;

        case ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO:                    // ZERO
            vRot = 0.0;
            break;

        default:                                                        // unknown prescription
            // The only way this can happen is if someone added a ROTATIONAL_VELOCITY_DISTRIBUTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_VROT_PRESCRIPTION);              // throw error
    }

    return vRot;
}


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
GNU_CONST inline double ConvertPeriodInDaysToSemiMajorAxisInAU(const double p_Mass1, const double p_Mass2, const double p_Period) {
    return std::cbrt((p_Mass1 + p_Mass2) * p_Period * p_Period / DAYS_IN_YEAR / DAYS_IN_YEAR);
}


/*
 * SolveKeplersEquation
 *
 * @brief
 * Solve Kepler's Equation using Newton-Raphson iteration.
 *
 * For a definition of the anomalies seen here:
 *
 *    https://en.wikipedia.org/wiki/Mean_anomaly
 *    https://en.wikipedia.org/wiki/True_anomaly
 *    https://en.wikipedia.org/wiki/Eccentric_anomaly
 *
 *
 * std::tuple<ERROR, double, double> SolveKeplersEquation(const double p_MeanAnomaly, const double p_Eccentricity)
 *
 * @param       p_MeanAnomaly                   The mean anomaly
 * @param       p_Eccentricity                  Eccentricity of the binary
 * @return                                      Tuple containing (in order): error value, eccentric anomaly, true anomaly
 *                                                  The error value returned will be:
 *                                                      NONE           if no error occurred
 *                                                      NO_CONVERGENCE if the Newton-Raphson iteration did not converge
 *                                                      OUT_OF_BOUNDS  if the eccentric anomaly returned is < 0 or > 2pi
 *                                                  If the error returned is not NONE, use the eccentric anomaly and 
 *                                                  true anomaly returned at your own risk
 */
GNU_PURE std::tuple<ERROR, double, double> SolveKeplersEquation(const double p_MeanAnomaly, const double p_Eccentricity) {

    ERROR error = ERROR::NONE;

    const double e = p_Eccentricity;
    const double M = p_MeanAnomaly;

    double E           = p_MeanAnomaly;                                                     // Initial guess at E is M - correct for e = 0
    double kepler      = E - (e * std::sin(E)) - M;                                         // Let f(E) = 0
    double keplerPrime = 1.0 - (e * std::cos(E));                                           // Derivative of f(E), f'(E)

    // Iterate until E is within the specified error of the true value, or max iterations exceeded
    SizeT iteration = 0;
    while (std::abs(kepler) >= NEWTON_RAPHSON_EPSILON && iteration++ < MAX_KEPLER_ITERATIONS) {
        E          -= kepler / keplerPrime;
        kepler      = E - (e * std::sin(E)) - M;                                            // Let f(E) = 0
        keplerPrime = 1.0 - (e * std::cos(E));                                              // Derivative of f(E), f'(E)
    }

    if (iteration >= MAX_KEPLER_ITERATIONS) error = ERROR::NO_CONVERGENCE;                  // No convergence

    double nu = 2.0 * std::atan((std::sqrt((1.0 + e) / (1.0 - e))) * std::tan(0.5 * E));    // Convert eccentric anomaly into true anomaly

    if (E >= M_PI && E <= _2_PI) nu += _2_PI;                                               // Add 2PI if necessary
    else if (E < 0.0 || E > _2_PI) error = ERROR::OUT_OF_BOUNDS;                            // E < 0 or E > 2pi?

    return std::make_tuple(error, E, nu);
}


/*
 * SolveQuadratic
 *
 * @brief
 * Solve quadratic Ax^2 + Bx + C.
 *
 * Returns either root, depending on discriminant will return:
 *
 *    0.0               if 0 roots
 *    root              if 1 root
 *    max(root1, root2) if 2 roots
 *
 *
 * std::tuple<ERROR, double> SolveQuadratic(const double p_A, const double p_B, double p_C)
 *
 * @param       p_A                             Coefficient of x^2
 * @param       p_B                             Coefficient of x^1
 * @param       p_C                             Coefficient of x^0 (Constant)
 * @return                                      Tuple containing (in order): error value, root found (see above)
 *                                                  The error value returned will be:
 *                                                      NONE          if no error occurred
 *                                                      NO_REAL_ROOTS if the equation has no real roots
 *                                                  If the error returned is not ERROR:NONE, use root returned at your own risk
 */
GNU_PURE std::tuple<ERROR, double> SolveQuadratic(const double p_A, const double p_B, double p_C) {

    ERROR error = ERROR::NONE;

    const double discriminant = (p_B * p_B) - (4.0 * p_A * p_C);    // d = B^2 - 4AC

    double root = 0.0;                                              // root found

    if (discriminant < 0.0) {                                       // real roots?
        error = ERROR::NO_REAL_ROOTS;                               // no real roots - set error
    }
    else if (discriminant > 0.0) {                                  // 2 real roots?
                                                                    // 2 real roots
        const double sqrtD = std::sqrt(discriminant);
        const double A2    = p_A + p_A;
        const double root1 = (-p_B + sqrtD) / A2;                   // (-B + SQRT(B^2 - 4AC)) / 2A
        const double root2 = (-p_B - sqrtD) / A2;                   // (-B - SQRT(B^2 - 4AC)) / 2A

        root = std::max(root1, root2);
    }
    else {                                                          // 1 real root
        root = -p_B / (p_A + p_A);                                  // -B / 2A,discriminant = 0.0
    }

    return std::make_tuple(error, root);
}


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
GNU_PURE inline SN_EVENT SNEventType(const SN_EVENT p_SNEvent) {

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


/*
 * to_underlying
 *
 * Cast enum element to underlying type.
 * Name is snake_case here because C++23 has a std::to_underlying(), so once
 * we move to C++23 or later we can just switch all uitils::to_underling()
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
//                                  ERROR HANDLING                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

    
/*
 * GetStackTrace
 *
 * @brief
 * Construct a vector of strings that represent the stack trace for the current thread
 * (COMPAS is single-threaded, so in our case, the current process executing COMPAS).
 * 
 * We use the gcc library functions to construct the stack trace (gcc calls this a "backtrace"):
 * 
 *     - backtrace(), which provides a list of pointers to each of the functions that make
 *       up the stack trace (the functions called, all the way from main() to the current
 *       point of execution - this is an instantaneous list, not historic).
 * 
 * 
 *     - backtrace_symbols(), translates the function pointers obtained from backtrace() into
 *       an array of strings that are the function names.  This typically only works for COMPAS
 *       functions (because we build COMPAS with debug info included).  We most likely won't have
 *       symbols for libraries (e.g. libc), so for non-COMPAS functions we insert the string
 *       "~~LIBFUNC~~" as the function name so users can identify non-COMPAS entries and handle
 *       them accordingly.
 * 
 * Returns a vector of strings representing the function names that comprise the stack trace.
 * 
 * 
 * StrVectorT GetStackTrace()
 * 
 * @return                                      Vector of strings representing the function names that comprise the stack trace
 */
StrVectorT GetStackTrace() {

    StrVectorT stackTrace = {};                                                                                 // Return vector

    void*  trace[MAX_STACK_TRACE_SIZE];                                                                         // Stack trace
    char** strings   = (char **)NULL;                                                                           // Stack trace strings
    SizeT  traceSize = 0;                                                                                       // Stack trace size
        
    traceSize = backtrace(trace, MAX_STACK_TRACE_SIZE);                                                         // Get stack trace size
    strings   = backtrace_symbols(trace, traceSize);                                                            // Get stack trace with symbols

    for (SizeT idx = 1; idx < traceSize; ++idx) {                                                               // For each stack trace entry
        // Extract function name.
        // We don't have symbols for libraries (e.g. libc), so for non-COMPAS functions we insert
        // "~~LIBFUNC~~" as the function name so the caller can identify non-COMPAS entries and
        // handle them accordingly
        SizeT start = 0;
        while (strings[idx][start] != '(' && strings[idx][start] != ' ' && strings[idx][start] != 0) ++start;   // Find function name start position
        SizeT end = start;
        while (strings[idx][end] != '+' && strings[idx][end] != 0) ++end;                                       // Find function name end position

        StrT funcName;                                                                                          // The extracted function name
        SizeT funcStrLen = end - start - 1;                                                                     // Length of (mangled) function string
        if (funcStrLen < 1) funcName = "~~LIBFUNC~~";                                                           // Library function
        else {                                                                                                  // Extract COMPAS function name
            char* funcStr = new char[funcStrLen + 1];                                                           // Allows for null terminator
            strncpy(funcStr, &strings[idx][start + 1], funcStrLen);                                             // Copy function name
            funcStr[funcStrLen] = 0;                                                                            // Make sure it is null-terminated
            funcName = StrT(funcStr);                                                                           // Function name

            int status = -1;
            char* demangledName = abi::__cxa_demangle(funcStr, NULL, NULL, &status);                            // Try to demangle the function name
            if (status == 0) funcName = StrT(demangledName);                                                    // Use the demangled name if available

            delete[] demangledName;
            delete[] funcStr;
        }

        stackTrace.push_back(funcName);                                                                         // Add function name to stacktrace

        if (funcName == "main") break;                                                                          // That's all we need
    }
    delete[] strings;

    return stackTrace;
}


/*
 * ShowStackTrace
 *
 * @brief
 * Display a stack trace, obtained by calling utils::GetStackTrace(), on stderr
 *
 * 
 * void ShowStackTrace()
 */
void ShowStackTrace() {

    StrVectorT stackTrace = utils::GetStackTrace();                     // Get stack trace

    if (!stackTrace.empty()) {                                          // Anything to show?
        std::cerr << "\nStack trace:\n";                                // Yes - display header
        for (SizeT entry = 1; entry < stackTrace.size(); entry++) {     // Ignore the eponymous entry
            std::cerr << "    " << stackTrace[entry] << "\n";           // Show stacktrace entry
        }
    }
    else {                                                              // No - nothing to show
        std::cerr << "\nNo stack trace available\n";                    // Say so
    }
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    FILE SYSTEM                                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CreateDirectory
 *
 * @brief
 * Create directories, where they don't already exist, in the path supplied
 *
 * We could just use boost::create_directories() (or FS::create_directories() for c++17 or greater), but that
 * would just create the directories as necessary and not report which directories in the path were pre-existing
 * and which were newly created.  We want to be able to clean up any directories we created but didn't use (because
 * we had an error somewhere perhaps...), and to do that we need to know which directories were actually created.
 * 
 * Returns a vector of paths created.
 * 
 * 
 * std::tuple<ERROR, StrT, StrVectorT> CreateDirectory(const StrT p_Path)
 * 
 * @param       p_Path                          Path specifying directories to be created
 * @return                                      Tuple containing error value, error string, and vector of directories created
 *                                                  The error value returned will be:
 *                                                      NONE                       if no error occurred
 *                                                      UNABLE_TO_CREATE_DIRECTORY if any non-existent directory in the path supplied could not be created
 *                                                  The error string indicates the problematic path in the case an error occurred.
 *                                                      The error string is only valid if error != NONE
 *                                                  The returned vector contains string paths for each of the directories actually created.
 *                                                      The returned vector will not contain names of directories in the path that already existed
 *                                                      The content of the returned vector for directories created is only valid if error == NONE,
 */
std::tuple<ERROR, StrT, StrVectorT> CreateDirectories(const StrT p_Path) {

    if (p_Path.empty()) return std::make_tuple(ERROR::NONE, "", StrVectorT({}));            // Nothing to do

    ERROR error             = ERROR::NONE;                                                  // Error - initially NONE
    StrT errStr             = "";                                                           // Error string - initially none
    StrVectorT pathsCreated = {};                                                           // Directories created - initially none

    // Create directories as necessary - top-down
    const char sep   = boost::filesystem::path::preferred_separator;                        // Platform-specific path separator
    StrT path       = "";                                                                   // Current path - the directory to be created
    for (const auto& name: boost::filesystem::path(p_Path)) {                               // Parse the user-supplied path string
        if (!path.empty() && path[path.length() - 1] != sep) path += sep;                   // Add separator to current path if necessary
        path += name.c_str();                                                               // Add directory name to current path
        if (!path.empty() && !boost::filesystem::exists(path)) {
            try {
                boost::system::error_code err;
                (void)boost::filesystem::create_directory(path, err);                       // Create directory - let boost throw an exception if it fails
                if (err.value() == 0) {                                                     // Ok?
                    pathsCreated.push_back(boost::filesystem::canonical(path).string());    // Yes - record creation
                }
                else  {                                                                     // Not ok...
                    error  = ERROR::UNABLE_TO_CREATE_DIRECTORY;                             // Set error
                    errStr = path;                                                          // ... and error string
                }            
            }
            catch (...) {                                                                   // Unhandled problem...
                error  = ERROR::UNABLE_TO_CREATE_DIRECTORY;                                 // Set error
                errStr = path;                                                              // ... and error string
            }                
        }

        if (error != ERROR::NONE) {                                                         // Problem?
            std::tie(error, errStr, pathsCreated) = RemoveDirectories(pathsCreated);        // Yes - clean up
            break;
        }
    }

    return std::make_tuple(error, errStr, pathsCreated);
}


/* 
 * FileExists
 *
 * @brief
 * Determine if a file with name as passed in p_Filename exists.
 * The parameter p_Filename is a character array.
 *
 *
 * bool FileExists(const char *p_Filename)
 *
 * @param       p_Filename                      Fully qualified filename (char*)
 * @return                                      Boolean indicating whether file exists
 */
inline bool FileExists(const char *p_Filename) { return (bool)(std::ifstream ifile(p_Filename)); }


/*
 * FileExists
 *
 * @brief
 * Determine if a file with name as passed in p_Filename exists.
 * The parameter p_Filename is a std::string.
 *
 *
 * bool FileExists(const std::string& p_Filename)
 *
 * @param       p_Filename                      Fully qualified filename (std::string)
 * @return                                      Boolean indicating whether file exists
 */
inline bool FileExists(const StrTg& p_Filename) { return FileExists(p_Filename.c_str()); }


/// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS - should this be in main.cpp?
/*
 * ReadTimesteps
 *
 * @brief
 * Read timesteps from timesteps file.
 *
 * Timesteps file is expected to be an ascii file with one timestep per record.
 * Timesteps must be > 0.0
 *  
 * 
 * std::tuple<ERROR, DblVectorT> ReadTimesteps(const StrT p_TimestepsFileName)
 * 
 * @param       p_TimestepsFileName             Filename to be read - should be fully qualified
 * @return                                      Tuple containing error value and timesteps vector
 *                                                  The error value returned will be:
 *                                                      NONE                                 if no error occurred
 *                                                      EMPTY_FILENAME                       if the filename provided was an empty string
 *                                                      FILE_DOES_NOT_EXIST                  if the timesteps file does not exist
 *                                                      FILE_OPEN_ERROR                      if the timesteps file exists but could not be opened
 *                                                      FILE_READ_ERROR                      if the timesteps file could not be read
 *                                                      EMPTY_FILE                           if the timesteps file contains no content
 *                                                      INVALID_VALUE_IN_FILE                if the file contains an invalid value for timestep
 *                                                      TOO_MANY_TIMESTEPS_IN_TIMESTEPS_FILE if the file contains too many timesteps (> maximum per OPTIONS)
 * 
 *                                                  If the error returned is not NONE, the content of the timesteps vector returned is not defined
 */
std::tuple<ERROR, DblVectorT> ReadTimesteps(const StrT p_TimestepsFileName) {

    ERROR error = ERROR::NONE;
 
    DblVectorT timesteps;                                                                                       // Timesteps vector

    if (p_TimestepsFileName.empty()) {                                                                          // Timesteps filename empty?
        error = ERROR::EMPTY_FILENAME;                                                                          // Yes - fail
    }
    else {

        if (!utils::FileExists(p_TimestepsFileName)) {                                                          // Timesteps file exists?
            error = ERROR::FILE_DOES_NOT_EXIST;                                                                 // No - fail
        }
        else {                                                                                                  // Yes
            std::ifstream timestepsFile;
            timestepsFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            try {
                timestepsFile.open(p_TimestepsFileName);                                                        // Open the timesteps file
                if (!timestepsFile.is_open()) {                                                                 // Open ok?
                    error = ERROR::FILE_OPEN_ERROR;                                                             // No - fail
                }
                else {                                                                                          // Yes - file open
                    std::string rec;                                                                            // Record read from file
                    SizeT numTimesteps = 0;
                    while (std::getline(timestepsFile, rec)) {                                                  // Get next record from timesteps file

                        if (rec.size() > 0 && (rec[rec.size() - 1] == '\n' || rec[rec.size() - 1] == '\r')) {   // Last character `\n` or `\r`?
                            rec.erase(rec.size() - 1);                                                          // Yes - strip it
                        }

                        rec = utils::trim(rec);                                                                 // Remove leading and trailing blanks

                        if (!(rec.empty() || rec[0] == '#')) {                                                  // Blank record or comment?                                 
                            try {                                                                               // No - process it
                                SizeT lastChar;
                                long double v = std::stold(rec, &lastChar);                                     // Try conversion
                                if (lastChar != (rec.size())) {                                                 // Conversion valid only if rec completely consumed
                                    error = ERROR::INVALID_VALUE_IN_FILE;                                       // Not a valid DOUBLE
                                    break;                                                                      // Stop processing
                                }

                                if (v < 0.0) {                                                                  // Timestep must be >= 0.0
                                    error = ERROR::INVALID_VALUE_IN_FILE;                                       // Not a valid timestep
                                    break;                                                                      // Stop processing
                                }
                                else {                                                                          // Ok - timestep >= 0.0
                                    timesteps.push_back(v);                                                     // Add timestep to timesteps vector
                                }
                    
                                numTimesteps++;                                                                 // Increment number of timesteps read
                                if (numTimesteps >= ABSOLUTE_MAXIMUM_TIMESTEPS) {                               // Number of timesteps exceeds maximum?
                                    error = ERROR::TOO_MANY_TIMESTEPS_IN_TIMESTEPS_FILE;                        // Yes - fail
                                    break;                                                                      // Stop processing
                                }
                            }
                            catch (const std::out_of_range& e) {                                                // Conversion failed
                                error = ERROR::INVALID_VALUE_IN_FILE;                                           // Not a valid DOUBLE
                                break;                                                                          // Stop processing
                            }
                            catch (const std::invalid_argument& e) {                                            // Conversion failed
                                error = ERROR::INVALID_VALUE_IN_FILE;                                           // Not a valid DOUBLE
                                break;                                                                          // Stop processing
                            }
                        }
                    }
                    try {
                        timestepsFile.close();                                                                  // Close the timesteps file
                    }
                    catch (std::ifstream::failure& e) {                                                         // Close failed
                        error = ERROR::FILE_NOT_CLOSED;                                                         // Fail
                    }
                }
            }
            catch (std::ifstream::failure& e) {                                                                 // Something was flagged...
                if (timestepsFile.eof()) {                                                                      // end-of-file?
                    if (timesteps.size() < 1) {                                                                 // Yes - at least one timestep read?
                        error = ERROR::EMPTY_FILE;                                                              // No - fail
                    }
                }
                else {                                                                                          // Not end-of-file - error
                    error = ERROR::FILE_READ_ERROR;                                                             // Fail
                }
            }

        }
    }

    return std::make_tuple(error, timesteps);
}


/*
 * RemoveDirectories
 *
 * @brief
 * Remove directories if they are empty
 *
 * Iterate in reverse order over the vector of paths passed (i.e. walk up the directory tree from the leaf) and
 * remove empty directories - but stop at the first non-empty directory (or error).  
 * 
 * Returns a vector of paths not removed.
 * 
 * 
 * std::tuple<ERROR, StrT, StrVectorT> RemoveDirectories(const StrVectorT p_Paths)
 * 
 * @param       p_Path                          Vector of paths to be removed
 * @return                                      Tuple containing error value, error string, and vector of directories not removed
 *                                                  The error value returned will be:
 *                                                      NONE                       if no error occurred
 *                                                      UNABLE_TO_CREATE_DIRECTORY if any non-existent directory in the path supplied could not be created
 *                                                  The error string indicates the problematic path in the case an error occurred.
 *                                                      The error string is only valid if error != NONE
 *                                                  The returned vector contains string paths for each of the directories from the input vector that were not removed,
 *                                                  either because they were not empty or because an error occurred (check error element of tuple returned)
 */
std::tuple<ERROR, StrT, StrVectorT> RemoveDirectories(const StrVectorT p_Paths) {

    if (p_Paths.size() < 1) return std::make_tuple(ERROR::NONE, "", StrVectorT({}));                                    // Nothing to do

    ERROR      error           = ERROR::NONE;                                                                           // Error - initially NONE
    StrT       errStr          = "";                                                                                    // Error string - initially none
    StrVectorT pathsNotRemoved = p_Paths;                                                                               // Directories not removed - initially p_Paths

    // Remove directories as necessary - bottom-up - stop if error or not empty.
    StrT path;
    boost::system::error_code err;
    for (auto iter = p_Paths.rbegin(); iter != p_Paths.rend(); ++iter) {
        path = *iter;
        if (boost::filesystem::exists(path) && boost::filesystem::is_empty(path)) {                                     // Directory exists and is empty?
            try {
                (void)boost::filesystem::remove(path, err);                                                             // Remove directory - let boost throw an exception if it fails
                if (err.value() == 0) {                                                                                 // Ok?
                                                                                                                        // Yes
                    // Erase path removed from pathsNotRemoved vector.
                    // I could just use std::find() as the parameter to erase(), but this is safer.
                    StrVectorT::iterator idx = std::find(pathsNotRemoved.begin(), pathsNotRemoved.end(), path);         // Get element index
                    if (idx != pathsNotRemoved.end()) pathsNotRemoved.erase(idx);                                       // If found, erase it
                }
                else {                                                                                                  // Not ok...
                    error  = ERROR::UNABLE_TO_REMOVE_DIRECTORY;                                                         // Set error
                    errStr = path;                                                                                      // ... and error string
                    break;
                }
            }
            catch (...) {                                                                                               // Unhandled problem...
                error  = ERROR::UNABLE_TO_REMOVE_DIRECTORY;                                                             // Set error
                errStr = path;                                                                                          // ... and error string
                break;
            }                
        }
        else {                                                                                                          // Problem...
            // Either the directory does not exist or is not empty.
            // if not empty we stop here - we don't want to delete anything other than empty directories.
            // if the directory does not exist we err on the side of caution, assume there's a problem, and stop here.
            error  = boost::filesystem::exists(path) ? ERROR::UNABLE_TO_REMOVE_DIRECTORY : ERROR::DIRECTORY_NOT_EMPTY;  // Set error
            errStr = path;                                                                                              // ... and error string
            break;
        }
    }

    return std::make_tuple(error, errStr, pathsNotRemoved);
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
 * Iterative (rather than recursive) binary search function.  Searches a std::vector of doubles
 * for a given double value.  The vector of doubles is expected to be sorted before calling this
 * function (sorting here would be prohibitive performance-wise - the expectation is that this
 * function could/will be called many times for the same (sorted) vector).
 * 
 * For a given double x, and a sorted array arr, the function returns a vector of (two) numbers
 * (of type SizeT) that represent the lower and upper bin edges of x in arr.
 *
 *
 * SizeTVectorT BinarySearch(const DblVectorT p_Arr, const double p_X)
 *
 * @param       p_Arr                           Sorted array to search over
 * @param       p_X                             Value to search for
 * @return                                      Vector containing SizeT indices of the lower and upper bin edges that contain p_X.
 *                                                  If x < min(Arr), return {-1, 0}
 *                                                  If x > max(Arr), return {0, -1}
 *                                                  If x is equal to an array element, return index of that element i.e. {idx, idx}
 */
GNU_CONST inline SizeTVectorT BinarySearch(const DblVectorT p_Arr, const double p_X) {

    SizeT low = 0;
    SizeT up  = p_Arr.size() - 1;
    SizeT mid = 0;

    // If p_X is not within p_Arr limits...
    if      (p_X < p_Arr[low]) { return {-1, 0}; }
    else if (p_X > p_Arr[up])  { return {0, -1}; }

    while(1) {                                                      // This cannot hang - eventually one of the returns must happen
        mid = roundl(0.5 * (up + low));
        if (std::abs(low - up) == 1) { return {low, low + 1}; }     // p_Arr(low) < p_X < p_Arr(up), so return low
        else if (p_X == p_Arr[low])  { return {low, low}; }         // p_Arr(low) = p_X. In this case, return low = up
        else if (p_X == p_Arr[up])   { return {up, up}; }           // p_Arr(up) = p_X. In this case, return low = up
        else if (p_X == p_Arr[mid])  { return {mid, mid}; }         // p_Arr(mid) = p_X. In this case, return low = up = mid
        else if (p_X < p_Arr[mid])   { up  = mid; }                 // Bring down upper bound
        else                         { low = mid; }                 // Bring up lower bound
    }
}


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
GNU_CONST std::tuple<bool, int> Find(const T &p_Elem, const std::vector<T> &p_Vector) {
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
GNU_CONST std::tuple<bool, E> GetMapKey(const StrT p_Value, const M& p_Map, const E& p_Default) {
    for (auto& it: p_Map)
        if (Equals(it.second, p_Value)) return std::make_tuple(true, it.first);
    return std::make_tuple(false, p_Default);
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
GNU_CONST inline bool IsOneOf(const STELLAR_TYPE p_StellarType, const StellarTypeListT p_List) {
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
GNU_CONST inline bool BracketTolerance(const double p_Bracket1, const double p_Bracket2) {
    const double diff = std::fabs(p_Bracket1 - p_Bracket2);                                       // Absolute value of difference
    const double min  = std::min(p_Bracket1, p_Bracket2);                                         // Minimum bracket value - could straddle 0.0
    return diff <= ROOT_ABS_TOLERANCE || std::fabs(diff / min) <= ROOT_REL_TOLERANCE;
}


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
namespace compare {
    inline int CompareImpl(double p_X, double p_Y, double p_AbsTol, double p_RelTol) {
        const double diff = std::fabs(p_X - p_Y);
        const double tol  = std::max(p_AbsTol, p_RelTol * std::max(std::fabs(p_X), std::fabs(p_Y)));

        if (diff <= tol) return 0;
        else             return p_X < p_Y ? -1 : 1; // 'else' is redundant, but I think reads better...
    }
}

inline int Compare(const double p_X, const double p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true) {

    const bool overrideTol = p_Tolerance >= 0.0;

    const double absTol = overrideTol ? (p_Absolute ? p_Tolerance : 0.0        ) : FLOAT_TOLERANCE_ABSOLUTE;
    const double relTol = overrideTol ? (p_Absolute ? 0.0         : p_Tolerance) : FLOAT_TOLERANCE_RELATIVE;

    return compare::CompareImpl(p_X, p_Y, absTol, relTol);
}

inline int Compare(const CDOUBLE& p_X, const CDOUBLE& p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true) {

    const bool overrideTol = p_Tolerance >= 0.0;

    const double absTol = overrideTol ? (p_Absolute ? p_Tolerance : 0.0        ) : CDOUBLE::AbsoluteTolerance();
    const double relTol = overrideTol ? (p_Absolute ? 0.0         : p_Tolerance) : CDOUBLE::RelativeTolerance();
        
    return compare::CompareImpl(p_X.Value(), p_Y.Value(), absTol, relTol);
}

inline int Compare(const CDOUBLE& p_X, const double p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true) {

    const bool overrideTol = p_Tolerance >= 0.0;

    const double absTol = overrideTol ? (p_Absolute ? p_Tolerance : 0.0        ) : std::min(CDOUBLE::AbsoluteTolerance(), FLOAT_TOLERANCE_ABSOLUTE);
    const double relTol = overrideTol ? (p_Absolute ? 0.0         : p_Tolerance) : std::min(CDOUBLE::RelativeTolerance(), FLOAT_TOLERANCE_RELATIVE);

    return compare::CompareImpl(p_X.Value(), p_Y, absTol, relTol);
}

inline int Compare(const double p_X, const CDOUBLE& p_Y, const double p_Tolerance = -1.0, const bool p_Absolute = true) {

    const bool overrideTol = p_Tolerance >= 0.0;

    const double absTol = overrideTol ? (p_Absolute ? p_Tolerance : 0.0        ) : std::min(CDOUBLE::AbsoluteTolerance(), FLOAT_TOLERANCE_ABSOLUTE);
    const double relTol = overrideTol ? (p_Absolute ? 0.0         : p_Tolerance) : std::min(CDOUBLE::RelativeTolerance(), FLOAT_TOLERANCE_RELATIVE);

    return compare::CompareImpl(p_X, p_Y.Value(), absTol, relTol);
}


/*
 * intPow
 *
 * @brief
 * Calculate x^y where x is double and y is an integer (positive or negative)
 *
 * Faster than std::pow() for integer exponent, n.  This function uses plain multiplication,
 * andcompletes in O(|n|) multiplications, and is the right choice for small |n| (|n| <= ~4).
 * For n > ~4, consider intPowL() (below).
 *
 *
 * double intPow(const double p_Base, const int p_Exponent)
 *
 * @param   [IN]    p_Base              Base - number to be raised to integer power
 * @param   [IN]    p_Exponent          Exponent - integer to which base should be raised
 * @return                              Base ^ Exponent
 */
GNU_CONST inline double intPow(const double p_Base, const int p_Exponent) {
    double result = 1.0;                            // Default (exponent = 0)
    int i = std::abs(p_Exponent);                   // Cache absolute exponent
    while (i-- > 0) result *= p_Base;               // Multiply
    return p_Exponent < 0 ? 1.0 / result : result;  // Invert if negative exponent
}


/*
 * intPowL
 *
 * @brief
 * Calculate x^y where x is double and y is an integer (positive or negative)
 *
 * Faster than std::pow() for integer exponent, n.  This function uses exponentiation by squaring,
 * and completes in O(log |n|) multiplications, and is the right choice for large |n| (|n| > ~4).
 * For n <= ~4, consider intPow() (above).
 *
 *
 * double intPowL(const double p_Base, const int p_Exponent)
 *
 * @param   [IN]    p_Base              Base - number to be raised to integer power
 * @param   [IN]    p_Exponent          Exponent - integer to which base should be raised
 * @return                              Base ^ Exponent
 */
GNU_CONST inline double intPowL(const double p_Base, const int p_Exponent) {
    double result = 1.0;                            // Default (exponent = 0)
    double base   = p_Base;                         // Working base, squared each iteration
    int i = std::abs(p_Exponent);                   // Bits of |p_Exponent| still to process
    while (i > 0) {                                 // Loop over bits of |p_Exponent|
        if (i & 1) result *= base;                  // Current bit set: accumulate this power
        base *= base;                               // Square base for the next bit position
        i >>= 1;                                    // Shift to next bit
    }
    return p_Exponent < 0 ? 1.0 / result : result;  // Invert if negative exponent
}


/*
 * IsBOOL
 *
 * @brief
 * Determine if the string passed as p_Str is a valid BOOL (as defined by Boost).
 *
 * In this context (the Boost context), a valid boolean is one of:
 * 
 *     - 0|1        ("0" or "1")
 *     - true|false ("true" or "false" - case insensitive)
 *     - yes|no     ("yes" or "no" - case insensitive)
 *     - on|off     ("on" or "off" - case insensitive)
 *
 * The function will retiurn one of {0, 1, 2, 3, 4, -1, -2, -3, -4} to indicate the result:
 * 
 *     0 = not a valid boolean
 *     1 = valid: 0|1
 *     2 = valid: true|false
 *     3 = valid: yes|no
 *     4 = valid: on|off
 * 
 * A positive return value indicates the boolean value is TRUE; a negative, FALSE.
 * 
 * 
 * SizeT IsBOOL(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - as described above
 */
GNU_CONST inline SizeT IsBOOL(const StrT p_Str) {

    if (p_Str.empty()) return 0;                    // Not valid: empty string

    if (utils::Equals(p_Str, "0")    ) return  1;   // Valid: 0|1       : TRUE
    if (utils::Equals(p_Str, "1")    ) return -1;   // Valid: 0|1       : FALSE
    if (utils::Equals(p_Str, "true") ) return  2;   // Valid: true|false: TRUE
    if (utils::Equals(p_Str, "false")) return -2;   // Valid: true|false: FALSE
    if (utils::Equals(p_Str, "yes")  ) return  3;   // Valid: yes|no    : TRUE
    if (utils::Equals(p_Str, "no")   ) return -3;   // Valid: yes|no    : FALSE
    if (utils::Equals(p_Str, "on")   ) return  4;   // Valid: on|off    : TRUE
    if (utils::Equals(p_Str, "off")  ) return -4;   // Valid: on|off    : FALSE

    return 0;                                       // Not valid
}


/*
 * IsDOUBLE
 *
 * @brief
 * Determine if the string passed as p_Str is a valid DOUBLE.
 *
 * In this context, to be a valid DOUBLE the string must convert to a double successfully
 * via the std::stod() function.
 * 
 * 
 * bool IsDOUBLE(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid DOUBLE, else FALSE
 */
GNU_CONST inline bool IsDOUBLE(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stod(p_Str, &lastChar);      // Try conversion

        result = lastChar == p_Str.size();      // Valid DOUBLE if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid DOUBLE
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid DOUBLE
    }

    return result;
}


/*
 * IsFLOAT
 *
 * @brief
 * Determine if the string passed as p_Str is a valid FLOAT.
 *
 * In this context, to be a valid FLOAT the string must convert to a float successfully
 * via the std::stof() function.
 * 
 * 
 * bool IsFLOAT(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid FLOAT, else FALSE
 */
GNU_CONST inline bool IsFLOAT(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stof(p_Str, &lastChar);      // Try conversion

        result = lastChar == p_Str.size();      // Valid FLOAT if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid FLOAT
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid FLOAT
    }

    return result;
}


/*
 * IsINT
 *
 * @brief
 * Determines if the string passed as p_Str is a valid INT.
 *
 * In this context, to be a valid INT the string must convert to an integer successfully
 * via the std::stoi() function.
 * 
 * 
 * bool IsINT(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid INT, else FALSE
 */
GNU_CONST inline bool IsINT(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stoi(p_Str, &lastChar);      // Try conversion

        result = lastChar == p_Str.size();      // Valid INT if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid INT
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid INT
    }

    return result;
}


/*
 * IsLONGDOUBLE
 *
 * @brief
 * Determine if the string passed as p_Str is a valid LONG DOUBLE.
 *
 * In this context, to be a valid LONG DOUBLE the string must convert to a long double successfully
 * via the std::stold() function.
 * 
 * 
 * bool IsLONGDOUBLE(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid LONG DOUBLE, else FALSE
 */
GNU_CONST inline bool IsLONGDOUBLE(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stold(p_Str, &lastChar);     // Try conversion

        result = lastChar == p_Str.size();      // Valid LONG DOUBLE if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid LONG DOUBLE
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid LONG DOUBLE
    }

    return result;
}


/*
 * IsLONGINT
 *
 * @brief
 * Determine if the string passed as p_Str is a valid LONG INT.
 *
 * In this context, to be a valid LONG INT the string must convert to a long integer
 * successfully via the std::stol() function.
 * 
 * 
 * bool IsLONGINT(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid LONGINT, else FALSE
 */
GNU_CONST inline bool IsLONGINT(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stol(p_Str, &lastChar);      // Try conversion

        result = lastChar == p_Str.size();      // Valid LONG INT if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid LONG INT
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid LONG INT
    }

    return result;
}


/*
 * Determine if the string passed as p_Str is a valid UNSIGNED LONG INT.
 *
 * In this context, to be a valid UNSIGNED LONG INT the string must convert to an unsigned long integer
 * successfully via the std::stoul() function.
 * 
 * 
 * bool IsULONGINT(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid UNSIGNED LONGINT, else FALSE
 */
GNU_CONST inline bool IsULONGINT(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stoul(p_Str, &lastChar);     // Try conversion

        result = lastChar == p_Str.size();      // Valid UNSIGNED LONG INT if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid UNSIGNED LONG INT
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid UNSIGNED LONG INT
    }

    return result;
}

    
///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     SAMPLING                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateCDFKroupa
 *
 * @brief
 * Calculate the value of the CDF of the Kroupa (2001) IMF at p_Mass.
 *
 * If p_Mass is outside the bounds of the IMF (< p_Min or >= p_Max), the returned CDF value will be 0.0.
 * 
 * 
 * double CalculateCDFKroupa(const double p_Mass, const double p_Max, const double p_Min)
 *
 * @param       p_Mass                          Mass at which to calculate the value of the CDF (Msol)
 * @param       p_Max                           IMF maximum
 * @param       p_Min                           IMF minimum
 * @return                                      CDF value
 */
GNU_CONST double CalculateCDFKroupa(const double p_Mass, const double p_Max, const double p_Min) {

    if ((p_Mass < p_Min) || (p_Mass >= p_Max)) return 0.0;  // Return 0.0 if mass is out of bounds of the IMF
    
    double vCDF;

    if (p_Min <= KROUPA_BREAK_1 &&
        p_Max >  KROUPA_BREAK_1 &&
        p_Max <= KROUPA_BREAK_2) {

        const double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
        const double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (PPOW(p_Max, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);

        const double C1 = 1.0 / (term1 + term2);
        const double C2 = C1 * KROUPA_BREAK_1_POWER_1_2;

        if (p_Mass >= p_Min && p_Mass < KROUPA_BREAK_1) {
            vCDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_1) - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
        }
        else if (p_Mass >= KROUPA_BREAK_1 && p_Mass < KROUPA_BREAK_2) {
            vCDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1)) +
                   ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);
        }
    }
    else if (p_Min <= KROUPA_BREAK_1 &&
             p_Max >  KROUPA_BREAK_2) {

        const double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
        const double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2);
        const double term3 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * KROUPA_BREAK_2_POWER_2_3 * (PPOW(p_Max, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);

        const double C1 = 1.0 / (term1 + term2 + term3);
        const double C2 = C1 * KROUPA_BREAK_1_POWER_1_2;
        const double C3 = C2 * KROUPA_BREAK_2_POWER_2_3;

        if (p_Mass >= p_Min && p_Mass < KROUPA_BREAK_1) {
            vCDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_1) - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
        }
        else if (p_Mass >= KROUPA_BREAK_1 && p_Mass < KROUPA_BREAK_2) {
            vCDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1)) +
                   ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);
        }
        else if (p_Mass >= KROUPA_BREAK_2 && p_Mass < p_Max) {
            vCDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1)) +
                   ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2) +
                   ONE_OVER_KROUPA_POWER_3_PLUS1 * C3 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);
        }
    }
    else if (p_Min >  KROUPA_BREAK_1 &&
             p_Min <= KROUPA_BREAK_2 &&
             p_Max >  KROUPA_BREAK_2) {

        const double term1 = ONE_OVER_KROUPA_POWER_2_PLUS1 * (KROUPA_BREAK_2_PLUS1_2 - PPOW(p_Min, KROUPA_POWER_PLUS1_2));
        const double term2 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_2_POWER_2_3 * (PPOW(p_Max, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);

        const double C2 = 1.0 / (term1 + term2);
        const double C3 = C2 * KROUPA_BREAK_2_POWER_2_3;

        if (p_Mass >= p_Min && p_Mass < KROUPA_BREAK_2) {
            vCDF = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_2) - PPOW(p_Min, KROUPA_POWER_PLUS1_2));
        }
        else if (p_Mass >= KROUPA_BREAK_2 && p_Mass < p_Max) {
            vCDF = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - PPOW(p_Min, KROUPA_POWER_PLUS1_2)) +
                   ONE_OVER_KROUPA_POWER_3_PLUS1 * C3 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);
        }
    }

    return vCDF;
}


/*
 * DrawKickDirection
 *
 * @brief
 * Draw the angular components of the supernova kick, theta and phi, from the distribution
 * specified by the by the user via program options.
 * 
 * In the current implementation, phi is always returned as RAND->Random() * _2_PI
 * 
 * 
 * Dbl_DblT DrawKickDirection(const KICK_DIRECTION_DISTRIBUTION p_Distribution, const double p_Power)
 * 
 * @param       p_Distribution                  Kick direction distribution specified by the user
 * @param       p_Power                         Exponent for power law (only used by the PWERLAW distribution)
 * @return                                      Tuple containing theta and phi
 */
Dbl_DblT DrawKickDirection(const KICK_DIRECTION_DISTRIBUTION p_Distribution, const double p_Power = 0.0) {

    constexpr double delta = 1.0 * DEGREE;                                          // Small angle in radians
    const double     phi   = RAND->Random() * _2_PI;                                // Angle in the plane - random angle in range [0, 2pi) in the plane

    double theta = 0.0;                                                             // Angle in/out of the plane (0.0 is in)
    switch (p_Distribution) {                                                       // Which kick direction distribution?

        case KICK_DIRECTION_DISTRIBUTION::ISOTROPIC:                                // ISOTROPIC
            // Draw theta isotropically
            theta = std::acos(1.0 - (2.0 * RAND->Random())) - M_PI_2;
            break;

        case KICK_DIRECTION_DISTRIBUTION::POWERLAW: {                               // POWERLAW
            // Draw theta according to a powerlaw
            // (power law power = 0 = isotropic, +infinity = kick along pole, -infinity = kick in plane)

            // Choose magnitude of power law distribution.
            // If using a negative power law that fails at 0, use a lower cutoff (currently 1E-6)
            double cosThetaMagnitude = utils::InverseSampleFromPowerLaw(p_Power, 1.0, 1E-6);
            if (p_Power < 0.0) cosThetaMagnitude = 1.0 - cosThetaMagnitude;

            // Calculating the magnitude of cos theta gives us no information about whether it was up or down,
            // so we choose a direction randomly here, then clamp cosTheta to [-1.0, 1.0]
            const double cosTheta = std::min(1.0, std::max(-1.0, (RAND->Random() < 0.5 ? -cosThetaMagnitude : cosThetaMagnitude)));

            theta = std::acos(cosTheta);                                            // Force kick to be out of the plane
            } break;

        case KICK_DIRECTION_DISTRIBUTION::INPLANE:                                  // INPLANE
            // Force the kick to be in the plane theta = 0
            theta = 0.0;                                                            // Force kick to be in the plane
            break;

        case KICK_DIRECTION_DISTRIBUTION::PERPENDICULAR:                            // PERPENDICULAR
            // Force kick to be along spin axis
            theta = RAND->Random() < 0.5 ? M_PI_2 : -M_PI_2;                        // pi/2 - choose up or down randomly
            break;

        case KICK_DIRECTION_DISTRIBUTION::POLES:                                    // POLES
            //Direct the kick in a small cone around the poles
            theta = std::fabs(RAND->RandomGaussian(delta));                         // Small random angle (could be 0.0)
            theta = RAND->Random() < 0.5 > M_PI_2 - theta : theta - M_PI_2;         // Close, or equal, to pi/2. Choose up or down randomly
            break;

        case KICK_DIRECTION_DISTRIBUTION::WEDGE:                                    // WEDGE
            // Direct kick into a wedge around the equator (theta = 0)
            theta = RAND->RandomGaussian(delta);                                    // Gaussian around 0 with a deviation delta
            break;

        default:                                                                    // Unexpected distribution
            // The only way this can happen are if someone added a kick direction distribution and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_KICK_DIRECTION_DISTRIBUTION);               // Throw error
    }

    return std::make_tuple(theta, phi);
}


/*
 * InverseSampleFromPowerLaw
 *
 * @brief
 * Draw sample from a power law distribution p(x) \propto x^(a) between x = p_Min and x = p_Max.
 * 
 * Returns p_Min if p_Min == p_Max.
 * 
 * In the general case (a != -1), the normalised CDF over [p_Min, p_Max] is
 * 
 *    F(x) = (x^(a+1) − p_Min^(a+1)) / (p_Max^(a+1) − p_Min^(a+1)) 
 * 
 * This breaks down for the special case of a == -1 for two reasons:
 * 
 *     (i) the integral \int x^a dx changes from x^(a+1)/(a+1) for the general case to ln(x) for a == -1.
 *    (ii) the general case has (a + 1) both as an exponent and a divisor, so would divide by 0 for a == -1.
 * 
 * 
 * double InverseSampleFromPowerLaw(const double p_Power, const double p_Max, const double p_Min)
 *
 * @param       p_Power                         The power for the power law
 * @param       p_Max                           Maximum of the X-interval from which to sample
 * @param       p_Min                           Minimum of the X-interval from which to sample
 * @return                                      Drawn sample
 */
double InverseSampleFromPowerLaw(const double p_Power, const double p_Max, const double p_Min) {

    // We use utils::Compare() here - if p_Min and p_Max are equal within our floating-point
    // comparison tolerance, we are unlikely to actually sample a value between them.
    if (utils::Compare(p_Min, p_Max) == 0) return p_Min;

    double result;

    // We use utils::Compare() here because if p_Power as passed is the result of a numerical
    // computation and is very close to -1.0, but not actually -1.0 (e.g. -0.9999999), a simple
    // test of (p_Power == -1.0) would fail, and the code path would default to the general case,
    // and there the result becomes neumerically unstable as (a + 1) approaches 0.0.
    if (utils::Compare(p_Power, -1.0) == 0) {                                   // p_Power very close to -1.0?
        result = std::exp(RAND->Random() * std::log(p_Max / p_Min)) * p_Min;    // Yes - special case
    }
    else {                                                                      // No - general case
        const double powerPlus1      = p_Power + 1.0;
        const double min_PowerPlus1 = PPOW(p_Min, powerPlus1);

        result = PPOW((RAND->Random() * (PPOW(p_Max, powerPlus1) - min_PowerPlus1) + min_PowerPlus1), 1.0 / powerPlus1);
    }

    return result;
}


/*
 * InverseSampleFromTabulatedCDF
 *
 * @brief
 * Inverse sample from tabulated CDF.
 *
 * Finds X given Y and an ordered map<double, double> of (X, Y) pairs.
 * Uses simple linear interpolation.
 * 
 * We expect the map passed to contain y values in the range [0.0, 1.0].
 * We assume the map passed is ordered/sorted in y.
 * 
 * Throws an error if:
 * 
 *    - the CDF table is empty
 *    - the min y value in the table is < 0.0 or the max y value in the table is > 1.0
 *    - if p_Y is outside the range [0.0, 1.0)
 *
 *
 * double InverseSampleFromTabulatedCDF(const double p_Y, const std::map<double, double> p_Table)
 *
 * @param       p_Y                             The Y value for which X is to be calculated
 * @param       p_Table                         The table to interpolate on
 * @return                                      Interpolated X value
 */
GNU_PURE double InverseSampleFromTabulatedCDF(const double p_Y, const std::map<double, double> p_Table) {

    if (p_Table.size() <= 0) THROW_ERROR(ERROR::EMPTY_MAP);                             // Bounds check - map must not be empty
    if (p_Y < 0.0 || p_Y >= 1.0) THROW_ERROR(ERROR::OUT_OF_BOUNDS);                     // Bounds check - p_Y must in [0.0, 1.0)

    double xInterp;

    std::map<double, double>::const_iterator iter;                                      // Iterator

    const double yMin = p_Table.begin()->second;                                        // Minimum y value (map is ordered)
    if (yMin < 0.0) THROW_ERROR(ERROR::OUT_OF_BOUNDS);                                  // Bounds check - min y must be >= 0.0

    const double yMax = p_Table.rbegin()->second;                                       // Maximum y value
    if (yMax > 1.0) THROW_ERROR(ERROR::OUT_OF_BOUNDS);                                  // Bounds check - max y must be <= 1.0

    const double y = yMin + (p_Y * (std::max(yMin, yMax) - std::min(yMax, yMin)));      // Normalise y - clamp to [yMin, yMax)

    for (iter = p_Table.begin(); (iter != p_Table.end() && iter->second < y); ++iter);  // Find y upper bound

    if (iter == p_Table.begin()) {                                                      // Upper bound is first entry?
        xInterp = iter->first;                                                          // Yes - return first entry x value
    }
    else {                                                                              // No - interpolate
        const double xAbove   = iter->first;                                            // Upper bound x value
        const double yAbove   = iter->second;                                           // Upper bound y value
        const double xBelow   = std::prev(iter)->first;                                 // Lower bound x value
        const double yBelow   = std::prev(iter)->second;                                // Lower bound y value
        const double gradient = (yAbove - yBelow) / (xAbove - xBelow);                  // Gradient

        xInterp = xBelow + ((y - yBelow) / gradient);                                   // Interpolate
    }

    return xInterp;
}


/*
 * SampleEccentricity
 *
 * @brief
 * Draw eccentricity from the distribution specified by the user via program options.
 * The drawn value will be in the range [p_Min, p_Max].
 * 
 * Returns p_Min if p_Min == p_Max.
 * 
 *
 * double SampleEccentricity(const ECCENTRICITY_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min)
 *
 * @param       p_Distribution                  Eccentricity distribution specified by the user
 * @param       p_Max                           Distribution maximum
 * @param       p_Min                           Distribution minimum
 * @return                                      Drawn eccentricity
 */
double SampleEccentricity(const ECCENTRICITY_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min) {

    // We use utils::Compare() here - if p_Min and p_Max are equal within our floating-point
    // comparison tolerance, we are unlikely to actually sample a value between them.
    if (utils::Compare(p_Min, p_Max) == 0) return p_Min;

    double eccentricity;

    switch (p_Distribution) {                                                       // Which eccentricity distribution?

        case ECCENTRICITY_DISTRIBUTION::ZERO:                                       // ZERO
            // All systems are initially circular i.e. have zero eccentricity
            eccentricity = 0.0;
            break;

        case ECCENTRICITY_DISTRIBUTION::FLAT:                                       // FLAT
            eccentricity = utils::InverseSampleFromPowerLaw(0.0, p_Max, p_Min);
            break;

        case ECCENTRICITY_DISTRIBUTION::THERMAL:                                    // THERMAL
            // Thermal eccentricity distribution: p(e) = 2e
            eccentricity = utils::InverseSampleFromPowerLaw(1.0, p_Max, p_Min);
            break;

        case ECCENTRICITY_DISTRIBUTION::GELLER2013:                                 // GELLER2013
            // M35 eccentricity distribution from Geller, Hurley and Mathieu 2013
            // Gaussian with mean 0.38 and sigma 0.23
            // http://iopscience.iop.org/article/10.1088/0004-6256/145/1/8/pdf
            // Sampling function taken from binpop.f in NBODY6

            // loop until we find a value between min and max - shouldn't loop forever...
            do {
                eccentricity = 0.23 * std::sqrt(-2.0 * log(RAND->Random())) * cos(_2_PI * RAND->Random()) + 0.38;
            } while (eccentricity < p_Min || eccentricity > p_Max);
            break;

        case ECCENTRICITY_DISTRIBUTION::DUQUENNOYMAYOR1991:                         // DUQUENNOYMAYOR1991
            // Eccentricity distribution from Duquennoy & Mayor (1991)
            // http://adsabs.harvard.edu/abs/1991A%26A...248..485D
            // Sampling function taken from binpop.f in NBODY6

            // loop until we find a value between min and max - shouldn't loop forever...
            do {
                eccentricity = 0.15 * std::sqrt(-2.0 * log(RAND->Random())) * cos(_2_PI * RAND->Random()) + 0.3;
            } while (eccentricity < p_Min or eccentricity > p_Max);
            break;

        case ECCENTRICITY_DISTRIBUTION::SANA2012:                                   // SANA2012
            // Eccentricity distribution from Sana et al. 2012
            // (http://science.sciencemag.org/content/sci/337/6093/444.full.pdf) distribution of eccentricities.
            // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
            // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

            eccentricity = utils::InverseSampleFromPowerLaw(-0.42, p_Max, p_Min);
            break;

        default:                                                                    // Unexpected distribution
            // The only way this can happen are if someone added an eccentricity distribution and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_ECCENTRICITY_DISTRIBUTION);               // Throw error
    }

    return eccentricity;
}


/*
 * SampleFromTabulatedCDF
 *
 * @brief
 * Sample from tabulated CDF
 *
 * Finds Y given X and an ordered map<double, double> of (X, Y) pairs
 * Uses simple linear interpolation.
 * 
 * Throws an error if the CDF table is empty or if p_X is outside the bounds of the table x value.
 * 
 * We expect the map passed to contain y values in the range [0.0, 1.0].
 * We assume the map passed is ordered/sorted in y.
 * 
 * Throws an error if:
 * 
 *    - the CDF table is empty
 *    - p_X is less than min y value in the table or greater than the max y value in the table
 * 
 *
 * double SampleFromTabulatedCDF(const double p_X, const std::map<double, double> p_Table)
 *
 * @param       p_X                             The X value for which Y is to be calculated
 * @param       p_Table                         The table to interpolate on (Ordered map<double, double> of (X, Y) pairs)
 * @return                                      Interpolated Y value
 */
GNU_PURE double SampleFromTabulatedCDF(const double p_X, const std::map<double, double> p_Table) {

    double yInterp;

    std::map<double, double>::const_iterator iter;                                      // Iterator

    const double xMin = p_Table.begin()->first;                                         // Minimum x value (map is ordered)
    const double xMax = p_Table.rbegin()->first;                                        // Maximum x value

    if (p_Table.size() <= 0) THROW_ERROR(ERROR::EMPTY_MAP);                             // Bounds check - map must not be empty
    if (p_X < xMin || p_Y > xMax) THROW_ERROR(ERROR::OUT_OF_BOUNDS);                    // Bounds check - p_X must in [xMin, xMax)

    for (iter = p_Table.begin(); (iter != p_Table.end() && iter->first < p_X); ++iter); // Find x upper bound

    if (iter == p_Table.begin()) {                                                      // Upper bound is first entry?
        yInterp = iter->second;                                                         // Yes - return first entry y value
    }
    else if (iter == p_Table.end()) {                                                   // No - last entry?
        yInterp = std::prev(iter)->second;                                              // Yes - return last entry y value
    }
    else {                                                                              // No - interpolate
        const double xAbove   = iter->first;                                            // Upper bound x value
        const double yAbove   = iter->second;                                           // Upper bound y value
        const double xBelow   = std::prev(iter)->first;                                 // Lower bound x value
        const double yBelow   = std::prev(iter)->second;                                // Lower bound y value
        const double gradient = (yAbove - yBelow) / (xAbove - xBelow);                  // Gradient

        yInterp = yBelow + ((p_X - xBelow) * gradient);                                 // Interpolate
    }

    return yInterp;
}


/*
 * SampleInitialMass
 *
 * @brief
 * Draw mass from the distribution specified by the user via program options.
 * The drawn value will be in the range [p_Min, p_Max].
 * 
 * returns p_Min if p_Min == p_Max.
 *
 *
 * double SampleInitialMass(const INITIAL_MASS_FUNCTION p_Distribution, const double p_Max, const double p_Min, const double p_Power)
 *
 * @param       p_Distribution                  The IMF to use to draw the mass
 * @param       p_Max                           IMF maximum
 * @param       p_Min                           IMF minimum
 * @param       p_Power                         IMF power (for IMF::POWERLAW)
 * @return                                      Sampled mass (Msol)
 */
double SampleInitialMass(const INITIAL_MASS_FUNCTION p_Distribution, const double p_Max, const double p_Min, const double p_Power = 0.0) {

    // We use utils::Compare() here - if p_Min and p_Max are equal within our floating-point
    // comparison tolerance, we are unlikely to actually sample a value between them.
    if (utils::Compare(p_Min, p_Max) == 0) return p_Min;

    double mass;

    switch (p_Distribution) {                                                           // Which IMF?

        case INITIAL_MASS_FUNCTION::SALPETER:                                           // SALPETER
            mass = utils::InverseSampleFromPowerLaw(SALPETER_POWER, p_Max, p_Min);      // Inverse sample from power-law
            break;

        case INITIAL_MASS_FUNCTION::POWERLAW:                                           // POWER LAW
            mass = utils::InverseSampleFromPowerLaw(p_Power, p_Max, p_Min);             // Inverse sample from power-law
            break;

        case INITIAL_MASS_FUNCTION::UNIFORM:                                            // UNIFORM
            // Convenience function for POWERLAW with slope of 0
            mass = RAND->Random(p_Min, p_Max);                                          // Random number in range [p_Min, p_Max]
            break;

        case INITIAL_MASS_FUNCTION::KROUPA:                                             // KROUPA
            // Determine where the user specified their minimum and maximum masses to generate
            if (p_Min <= KROUPA_BREAK_1 && p_Max <= KROUPA_BREAK_1) {
                mass = utils::InverseSampleFromPowerLaw(KROUPA_POWER_1, p_Max, p_Min);  // Inverse sample from power-law 
            }
            else if (p_Min > KROUPA_BREAK_1 && p_Min <= KROUPA_BREAK_2 && p_Max > KROUPA_BREAK_1 && p_Max <= KROUPA_BREAK_2) {
                mass = utils::InverseSampleFromPowerLaw(KROUPA_POWER_2, p_Max, p_Min);  // Inverse sample from power-law
            }
            else if (p_Min > KROUPA_BREAK_2 && p_Max > KROUPA_BREAK_2) {
                mass = utils::InverseSampleFromPowerLaw(KROUPA_POWER_3, p_Max, p_Min);  // Inverse sample from power-law
            }
            else if (p_Min <= KROUPA_BREAK_1 && p_Max > KROUPA_BREAK_1 && p_Max <= KROUPA_BREAK_2) {
                const double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
                const double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (PPOW(p_Max, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);
                const double C1    = 1.0 / (term1 + term2);
                const double C2    = C1 * KROUPA_BREAK_1_POWER_1_2;
                const double A     = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
                const double rand  = RAND->Random();

                if (rand < CalculateCDFKroupa(KROUPA_BREAK_1, p_Max, p_Min))
                    mass = PPOW(rand * (KROUPA_POWER_PLUS1_1 / C1) + PPOW(p_Min, KROUPA_POWER_PLUS1_1), ONE_OVER_KROUPA_POWER_1_PLUS1);
                else
                    mass = PPOW((rand - A) * (KROUPA_POWER_PLUS1_2 / C2) + KROUPA_BREAK_1_PLUS1_2, ONE_OVER_KROUPA_POWER_2_PLUS1);
            }
            else if (p_Min <= KROUPA_BREAK_1 && p_Max > KROUPA_BREAK_2_POWER_2_3) {
                const double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
                const double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2);
                const double term3 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * KROUPA_BREAK_2_POWER_2_3 * (PPOW(p_Max, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);
                const double C1    = 1.0 / (term1 + term2 + term3);
                const double C2    = C1 * KROUPA_BREAK_1_POWER_1_2;
                const double C3    = C2 * KROUPA_BREAK_2_POWER_2_3;
                const double A     = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
                const double B     = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2);
                const double rand  = RAND->Random();

                if (rand < CalculateCDFKroupa(KROUPA_BREAK_1, p_Max, p_Min))
                    mass = PPOW(rand * (KROUPA_POWER_PLUS1_1 / C1) + PPOW(p_Min, KROUPA_POWER_PLUS1_1), ONE_OVER_KROUPA_POWER_1_PLUS1);
                else if (rand < CalculateCDFKroupa(KROUPA_BREAK_2, p_Max, p_Min))
                    mass = PPOW((rand - A) * (KROUPA_POWER_PLUS1_2 / C2) + KROUPA_BREAK_1_PLUS1_2, ONE_OVER_KROUPA_POWER_2_PLUS1);
                else
                    mass = PPOW((rand - A - B) * (KROUPA_POWER_PLUS1_3 / C3) + KROUPA_BREAK_2_PLUS1_3, ONE_OVER_KROUPA_POWER_3_PLUS1);
            }
            else if (p_Min > KROUPA_BREAK_1 && p_Min <= KROUPA_BREAK_2 && p_Max > KROUPA_BREAK_2) {
                const double term1 = ONE_OVER_KROUPA_POWER_2_PLUS1 * (KROUPA_BREAK_2_PLUS1_2 - PPOW(p_Min, KROUPA_POWER_PLUS1_2));
                const double term2 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_2_POWER_2_3 * (PPOW(p_Max, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);
                const double C2    = 1.0 / (term1 + term2);
                const double C3    = C2 * KROUPA_BREAK_2_POWER_2_3;
                const double B     = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - PPOW(p_Min, KROUPA_POWER_PLUS1_2));
                const double rand  = RAND->Random();

                if (rand < CalculateCDFKroupa(KROUPA_BREAK_2, p_Max, p_Min))
                    mass = PPOW(rand * (KROUPA_POWER_PLUS1_2 / C2) + PPOW(p_Min, KROUPA_POWER_PLUS1_2), ONE_OVER_KROUPA_POWER_2_PLUS1);
                else
                    mass = PPOW((rand - B) * (KROUPA_POWER_PLUS1_3 / C3) + KROUPA_BREAK_2_PLUS1_3, ONE_OVER_KROUPA_POWER_3_PLUS1);
            }
            // No other case possible - as long as p_Min < p_Max (currently enforced in Options.cpp)
            break;

        default:                                                                        // Unexpected distribution
            // The only way this can happen are if someone added an initial mass function and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_INITIAL_MASS_FUNCTION);                       // Throw error
    }

    return mass;
}


/*
 * SampleMassRatio
 *
 * @brief
 * Draw mass ratio, q, from the distribution specified by the user via program options.
 * The drawn value will be in the range [p_Min, p_Max].
 * 
 * Returns p_Min if p_Min == p_Max.
 *
 *
 * double SampleMassRatio(const MASS_RATIO_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min)
 *
 * @param       p_Distribution                  The distribution to use to draw the ratio
 * @param       p_Max                           Distribution maximum
 * @param       p_Min                           Distribution minimum
 * @return                                      Mass ratio q
 */
double SampleMassRatio(const MASS_RATIO_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min) {

    // We use utils::Compare() here - if p_Min and p_Max are equal within our floating-point
    // comparison tolerance, we are unlikely to actually sample a value between them.
    if (utils::Compare(p_Min, p_Max) == 0) return p_Min;

    double q;

    switch (p_Distribution) {                                                           // Which q distribution?

        case MASS_RATIO_DISTRIBUTION::FLAT:                                             // FLAT
            q = utils::InverseSampleFromPowerLaw(0.0, p_Max, p_Min);
            break;

        case MASS_RATIO_DISTRIBUTION::DUQUENNOYMAYOR1991:                               // DUQUENNOYMAYOR1991                                              
            // Mass ratio distribution from Duquennoy & Mayor (1991) (http://adsabs.harvard.edu/abs/1991A%26A...248..485D)

            // loop until we find a value between min and max - shouldn't loop forever...
            do {
                q = 0.42 * std::sqrt(-2.0 * log(RAND->Random())) * cos(_2_PI * RAND->Random()) + 0.23;
            } while (q < p_Min || q > p_Max);
            break;

        case MASS_RATIO_DISTRIBUTION::SANA2012:                                         // SANA2021                                                    
            // Sana et al. 2012 (http://science.sciencemag.org/content/sci/337/6093/444.full.pdf) distribution of eccentricities.
            // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
            // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

            q = utils::InverseSampleFromPowerLaw(-0.1, p_Max, p_Min);                   // de Mink and Belczynski use min = 0.1, max = 1.0
            break;

        default:                                                                        // Unexpected distribution
            // The only way this can happen are if someone added a mass ratio distribution and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_MASS_RATIO_DISTRIBUTION);                     // Throw error
    }

    return std::min(std::max(p_Min, q), p_Max);                                         // Clamp to [min, max]
}


/*
 * SampleMetallicity
 *
 * @brief
 * Draw metallicity from the distribution specified by the user via program options.
 * The drawn value will be in the range [p_Min, p_Max].
 * 
 * p_Min and p_Max are ignored for p_Distribution == ZSOLAR.
 * 
 * Returns p_Min if p_Min == p_Max unless p_Distribution == ZSOLAR.
 *
 *
 * double SampleMetallicity(const METALLICITY_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min)
 *
 * @param       p_Distribution                  The metallicity distribution to use to draw the metallicity
 * @param       p_Max                           Distribution maximum
 * @param       p_Min                           Distribution minimum
 * @return                                      Metallicity
 */
double SampleMetallicity(const METALLICITY_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min) {

    double metallicity;

    switch (p_Distribution) {                                                           // which Z distribution?

        case METALLICITY_DISTRIBUTION::ZSOLAR:                                          // ZSOLAR
            // All systems have Z = ZSOLAR (ZSOL_ASPLUND)
            metallicity = ZSOL_ASPLUND;
            break;

        case METALLICITY_DISTRIBUTION::LOGUNIFORM: {                                    // LOGUNIFORM
            // Sample Z uniformly in log

            // We use utils::Compare() here - if p_Min and p_Max are equal within our floating-point
            // comparison tolerance, we are unlikely to actually sample a value between them.
            if (utils::Compare(p_Min, p_Max) == 0) metallicity = p_Min;
            else {
                const double logMin = std::log10(p_Min);
                const double logMax = std::log10(p_Max);
                metallicity = PPOW(10, (logMin + ((logMax - logMin) * RAND->Random())));
            }
        } break;

        default:                                                                        // Unexpected distribution
            // The only way this can happen are if someone added a metallicity distribution and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_METALLICITY_DISTRIBUTION);                    // Throw error
    }

    return metallicity;
}


/*
 * SampleOrbitalPeriod
 *
 * @brief
 * Draw orbital period from the distribution specified by the user via program options.
 * The drawn value will be in the range [p_Min, p_Max].
 * 
 * Returns p_Min if p_Min == p_Max.
 * 
 * 
 * double SampleOrbitalPeriodDistribution(const ORBITAL_PERIOD_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min)
 *
 * @param       p_Distribution                  The distribution to use to draw orbital period
 * @param       p_Max                           Distribution maximum
 * @param       p_Min                           Distribution minimum
 * @return                                      Orbital period (days)
 */
double SampleOrbitalPeriod(const ORBITAL_PERIOD_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min) {

    // We use utils::Compare() here - if p_Min and p_Max are equal within our floating-point
    // comparison tolerance, we are unlikely to actually sample a value between them.
    if (utils::Compare(p_Min, p_Max) == 0) return p_Min;

    double orbitalPeriod;

    switch (p_Distribution) {                                                           // Which distribution?

        case ORBITAL_PERIOD_DISTRIBUTION::FLATINLOG:                                    // FLATINLOG

            orbitalPeriod = utils::InverseSampleFromPowerLaw(-1.0, p_Max, p_Min);
            break;

        default:                                                                        // Unexpected distribution
            // The only way this can happen are if someone added an orbital period distribution and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_ORBITAL_PERIOD_DISTRIBUTION);                 // Throw error
    }

    return orbitalPeriod;
}


/*
 * SampleSemiMajorAxis
 *
 * @brief
 * Draw semi-major axis from the distribution specified by the user via program options.
 * The drawn value will be in the range [p_aMin, p_aMax].
 * 
 * Returns p_aMin if p_aMin == p_aMax.
 * 
 * For the DUQUENNOYMAYOR1991 distribution, the functions returns a boolean flag indicating whether
 * the search for a value converged in the allowed number of iterations (see code below).
 * 
 * 
 * std::tuple<double, bool> SampleSemiMajorAxis(
 *     const SEMI_MAJOR_AXIS_DISTRIBUTION p_Distribution, 
 *     const double                       p_aMax, 
 *     const double                       p_aMin, 
 *     const double                       p_pMax = 0.0, 
 *     const double                       p_pMin = 0.0, 
 *     const double                       p_Mass1 = 0.0, 
 *     const double                       p_Mass2 = 0.0)
 *
 * @param       p_Distribution                  The distribution to use to draw semi-major axis
 * @param       p_aMax                          Semi-major axis distribution maximum
 * @param       p_aMin                          Semi-major axis distribution minimum
 * @param       p_pMax                          Orbital period distribution maximum (only used for SANA2012 distribution)
 * @param       p_pMin                          Orbital period distribution minimum (only used for for SANA2012 distribution)
 * @param       p_Mass1                         Mass of the primary star (Msol) (only used for SANA2012 distribution)
 * @param       p_Mass2                         Mass of the secondary star (Msol) (only used for SANA2012 distribution)
 * @return                                      Tuple containing:
 *                                                  Semi-major axis (AU)
 *                                                  Booelan flag indicating convergence
 *                                                      Only applies to DUQUENNOYMAYOR1991 distribution - will be TRUE otherwise
 *                                                      If no convergence, the value returned should be used with caution
 */
std::tuple<double, bool> SampleSemiMajorAxis(
    const SEMI_MAJOR_AXIS_DISTRIBUTION p_Distribution, 
    const double                       p_aMax, 
    const double                       p_aMin, 
    const double                       p_Pmax = 0.0, 
    const double                       p_Pmin = 0.0, 
    const double                       p_Mass1 = 0.0, 
    const double                       p_Mass2 = 0.0) {

    // We use utils::Compare() here - if p_Min and p_Max are equal within our floating-point
    // comparison tolerance, we are unlikely to actually sample a value between them.
    if (utils::Compare(p_aMin, p_aMax) == 0) return p_aMin;

    double semiMajorAxis;
    double convergence = true;                                                                                      // Default value

    switch (p_Adist) {                                                                                              // Which a distribution?

        case SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG:                                                               // FLATINLOG
            semiMajorAxis = utils::InverseSampleFromPowerLaw(-1.0, p_Max, p_Min);
            break;

        case SEMI_MAJOR_AXIS_DISTRIBUTION::DUQUENNOYMAYOR1991: {                                                    // DUQUENNOYMAYOR1991
            // Duquennoy & Mayor (1991) period distribution
            // http://adsabs.harvard.edu/abs/1991A%26A...248..485D
            // See also the period distribution (Figure 1) of M35 in Geller+ 2013 https://arxiv.org/abs/1210.1575
            // See also the period distribution (Figure 13) of local solar type binaries from Raghavan et al. 2010 https://arxiv.org/abs/1007.0414
            // They have log-normal distribution with a mean of 5.03 and a standard deviation of 2.28, with a minimum period of around 0.1 days
            // Sampling function taken from binpop.f in NBODY6

            // Make sure that the drawn semi-major axis is in the range specified by the user
            SizeT samples = 0;
            do {
                const double p = PPOW(10.0, 2.3 * std::sqrt(-2.0 * log(RAND->Random())) * std::cos(_2_PI * RAND->Random()) + 4.8);
                semiMajorAxis  = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(p_Mass1, p_Mass2, p);                // Convert period (days) to semi-major axis (AU)
            } while (samples++ < SEMI_MAJOR_AXIS_SAMPLES && (semiMajorAxis < p_aMin || semiMajorAxis > p_aMax));

            if (samples >= SEMI_MAJOR_AXIS_SAMPLES) convergence = false;                                            // Signal no convergence to caller
            } break;

        case SEMI_MAJOR_AXIS_DISTRIBUTION::SANA2012: {                                                              // SANA2012
            // Sana et al. 2012
            // http://science.sciencemag.org/content/sci/337/6093/444.full.pdf
            // distribution of semi-major axes. Sana et al. fit for the orbital period, which we sample in here, before returning the semi major axis
            // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
            // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

            const double logpMin = p_pMin > 1.0 ? log(p_pMin) : 0.0;                                                // Lower bound
            const double logpMax = p_pMax > 1.0 ? log(p_pMax) : 0.0;                                                // Upper bound

            const double p = std::exp(utils::InverseSampleFromPowerLaw(-0.55, logpMax, logpMin));                   // Draw period from their distribution
            semiMajorAxis  = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(p_Mass1, p_Mass2, p);                    // Convert period (days) to semi-major axis (AU)
            } break;

        default:                                                                                                    // Unexpected distribution
            // The only way this can happen are if someone added a semi-major axis distribution and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_SEMI_MAJOR_AXIS_DISTRIBUTION);                                            // Throw error
        }

        return std::make_tuple(semiMajorAxis, convergence);
    }

    
///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      STRINGS                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CentreJustify
 *
 * @brief
 * Centre-justifies string to specified width by prepending and appending spaces.
 * Extra space will be at the end of the string if necessary.
 * Used primarily in logging code.
 *
 *
 * StrT CentreJustify(const StrT p_Str, const SizeT p_Width)
 *
 * @param       p_Str                           String to be centre-justified
 * @param       p_Width                         The required width of the resultant string
 * @return                                      String padded with leading and trailing spaces so as to (as close as possible) centre-justify p_Str
 *                                              The string returned will always be p_Width characters in length
 */
 GNU_CONST inline StrT CentreJustify(const StrT p_Str, const SizeT p_Width) {

    StrT result = p_Str;                                                        // Default is no change

    if (p_Str.length() < p_Width) {                                             // p_Str < field width?
                                                                                // Yes
        SizeT numLeadingSpaces = (p_Width - p_Str.length()) / 2;                // Number of spaces to add at start - half the deficit
        SizeT numTralingSpaces = p_Width - p_Str.length() - numLeadingSpaces;   // Number of spaces to add at the end - whatever is left
        StrT leadingSpaces(numLeadingSpaces, ' ');                              // Blank string to prepend
        StrT trailingSpaces(numTralingSpaces, ' ');                             // Blank string to append
        result = leadingSpaces + p_Str + trailingSpaces;                        // Add leading and trailing spaces to p_Str
    }

    return result;
}


/*
 * Equals
 *
 * @brief
 * Case-insensitive comparison of strings.
 * This only works with ASCII data, but I think that's all we need.
 * Note that std::string has an == operator to test for equality (actually calls std::strcmp)
 *
 *
 * bool Equals(StrT p_Str1, StrT p_Str2)
 *
 * @param       p_Str1                          String to be compared
 * @param       p_Str2                          String to be compared
 * @return                                      Boolean indicating equality (true = equal)
 */
GNU_CONST inline bool Equals(StrT p_Str1, StrT p_Str2) {
    std::transform(p_Str1.begin(), p_Str1.end(), p_Str1.begin(), ::tolower);
    std::transform(p_Str2.begin(), p_Str2.end(), p_Str2.begin(), ::tolower);
    return (std::strcmp(p_Str1.c_str(), p_Str2.c_str()) == 0);
}


/*
 * PadLeadingZeros
 *
 * @brief
 * Pads string to specified length by prepending the string with "0".
 * This only works with ASCII data, but I think that's all we need.
 *
 *
 * StrT PadLeadingZeros(const StrT p_Str, const SizeT p_MaxLength)
 *
 * @param       p_Str                           String to be padded with leading "0"s
 * @param       p_MaxLength                     The required length of the resultant string
 * @return                                      String padded with leading "0"s
 *                                              Will be unchanged from input string if length already >= required length
 */
GNU_CONST inline StrT PadLeadingZeros(const StrT p_Str, const SizeT p_MaxLength) {
    return (p_Str.length() < p_MaxLength) ? StrT(p_MaxLength - p_Str.length(), '0') + p_Str : p_Str;
}


/*
 * PadTrailingSpaces
 *
 * @brief
 * Pads string to specified length by appending the string with " ".
 * This only works with ASCII data, but I think that's all we need.
 *
 *
 * StrT PadTrailingSpaces(const StrT p_Str, const SizeT p_MaxLength)
 *
 * @param       p_Str                           String to be padded with trailing " "s
 * @param       p_MaxLength                     The required length of the resultant string
 * @return                                      String padded with leading "0"s
 *                                              Will be unchanged from input string if length already >= required length
 */
GNU_CONST inline StrT PadTrailingSpaces(const StrT p_Str, const SizeT p_MaxLength) {
    return (p_Str.length() < p_MaxLength) ? p_Str + StrT(p_MaxLength - p_Str.length(), ' ') : p_Str;
}


/*
 * ltrim
 *
 * @brief
 * Trim leading whitespace characters from a string.
 *
 *
 * StrT& ltrim(StrT& p_Str)
 *
 * @param       p_Str                           String to be trimmed of whitespace
 * @return                                      Trimmed string
 */
GNU_CONST inline StrT& ltrim(StrT& p_Str) {
    p_Str.erase(0, p_Str.find_first_not_of("\t\n\v\f\r "));
    return p_Str;
}


/*
 * rtrim
 *
 * @brief
 * Trim trailing whitespace characters from a string.
 *
 *
 * StrT& rtrim(StrT& p_Str)
 *
 * @param       p_Str                           String to be trimmed of whitespace
 * @return                                      Trimmed string
 */
GNU_CONST inline StrT& rtrim(StrT& p_Str) {
    p_Str.erase(p_Str.find_last_not_of("\t\n\v\f\r ") + 1);
    return p_Str;
}


/*
 * trim
 *
 * @brief
 * Trim both leading and trailing whitespace characters from a string.
 *
 *
 * StrT& trim(StrT& p_Str)
 *
 * @param       p_Str                           String to be trimmed of whitespace
 * @return                                      Trimmed string
 */
GNU_CONST inline StrT& trim(StrT& p_Str) { return utils::ltrim(utils::rtrim(p_Str)); }


/*
 * ToLower
 *
 * @brief
 * Downshift uppercase characters in string to lowercase.
 * This only works with ASCII data, but I think that's all we need,
 *
 *
 * StrT ToLower(StrT p_Str)
 *
 * @param       p_Str                           String to be downshifted
 * @return                                      Downshifted string
 */
GNU_CONST inline StrT ToLower(StrT p_Str) {
    std::transform(p_Str.begin(), p_Str.end(), p_Str.begin(), ::tolower);
    return p_Str;
}


/*
 * ToUpper
 *
 * @brief
 * Upshift lowercase characters in string to uppercase.
 * This only works with ASCII data, but I think that's all we need.
 *
 *
 * StrT ToUpper(StrT p_Str)
 *
 * @param       p_Str                           String to be upshifted
 * @return                                      Upshifted string
 */
GNU_CONST inline StrT ToUpper(StrT p_Str) {
    std::transform(p_Str.begin(), p_Str.end(), p_Str.begin(), ::toupper);
    return p_Str;
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< SHOULD THIS BE IN LOGGING - CALLED ANYWHERE ELSE????? JR FIX THIS
/*
 * vFormat
 *
 * @brief
 * Formats value per sprintf() and returns string.
 * Used primarily in logging code.
 *
 * From https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf/49812018#49812018
 *
 *
 * const StrT vFormat(const char* const p_zcFormat, ...)
 *
 * @param       p_zcFormat                      Format string
 * @param       ...                             Parameters to be formatted
 * @return                                      Formatted string
 */
GNU_CONST const StrT vFormat(const char* const p_zcFormat, ...) {

    // Initialize use of the variable argument array
    va_list vaArgs;
    va_start(vaArgs, p_zcFormat);

    // Reliably acquire the size from a copy of the variable argument array
    // and a functionally reliable call to mock the formatting
    va_list vaArgsCopy;
    va_copy(vaArgsCopy, vaArgs);
    const int iLen = std::vsnprintf(NULL, 0, p_zcFormat, vaArgsCopy);
    va_end(vaArgsCopy);

    // Return a formatted string without risking memory mismanagement
    // and without assuming any compiler or platform specific behavior
    std::vector<char> zc(iLen + 1);
    std::vsnprintf(zc.data(), zc.size(), p_zcFormat, vaArgs);
    va_end(vaArgs);

    return StrT(zc.data(), iLen);
}


}
