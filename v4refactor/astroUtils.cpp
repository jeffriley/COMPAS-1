#include "astroUtils.h"

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


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       ORBIT                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

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
std::tuple<ERROR, double, double> SolveKeplersEquation(const double p_MeanAnomaly, const double p_Eccentricity) {

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


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ROTATION                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateRotationalVelocityOStar_Ramirez2013
 *
 * @brief
 * Calculate rotational velocity from the analytic cumulative distribution function (CDF)
 * for the equatorial rotational velocity of single O stars, per Ramirez-Agudelo et al. 2013
 * (see https://arxiv.org/abs/1309.2929)
 *
 * Uses inverse sampling and root finding.
 *
 *
 * double CalculateRotationalVelocityOStar_Ramirez2013()
 *
 * @return                                      Rotational velocity (km s^-1)
 */
double CalculateRotationalVelocityOStar_Ramirez2013() {

    double desiredCDF            = RAND->Random();                                                  // Random desired CDF
    const boost::uintmax_t maxit = ADAPTIVE_RV_MAX_ITERATIONS;                                      // Limit to maximum iterations.
    boost::uintmax_t it          = maxit;                                                           // Initially our chosen max iterations, but updated with actual

    double guess      = 100.0;                                                                      // Guess at 100 km s^-1 (arbitrary initial guess)
    double factorFrac = ADAPTIVE_RV_SEARCH_FACTOR_FRAC;                                             // Search step size factor fractional part
    double factor     = 1.0 + factorFrac;                                                           // Factor to determine search step size (size = guess * factor)

    std::pair<double, double> root(-1.0, -1.0);                                                     // Initialise root - default return
    SizeT tries = 0;                                                                                // Number of tries
    bool done   = false;                                                                            // Finished (found root or exceed maximum tries)?
    ERROR error = ERROR::NONE;
    OStarRotationVelocityFunctor<double> func = OStarRotationVelocityFunctor<double>(desiredCDF);
    while (!done) {                                                                                 // While no error and acceptable root found

        bool isRising = true;                                                                       // Guess for direction of search; CDF increases monotonically

        try {
            error = ERROR::NONE;
            root  = boost::math::tools::bracket_and_solve_root(func, guess, factor, isRising, utils::BracketTolerance, it); // Find root
            if (error != ERROR::NONE) { SHOW_WARN_STATIC(error); }                                  // Root finder encountered an error
            else if (it >= maxit) { SHOW_WARN_STATIC(ERROR::TOO_MANY_RV_ITERATIONS); }              // Too many root finder iterations
        }
        catch(std::exception& e) {                                                                  // Catch generic boost root finding error
            if (it >= maxit) { SHOW_WARN_STATIC(ERROR::TOO_MANY_RV_ITERATIONS); }                   // Too many root finder iterations
            else             { SHOW_WARN_STATIC(ERROR::ROOT_FINDER_FAILED, e.what()); }             // Some other problem - show it as a warning
        }

        // We have a solution from the root finder - it may not be an acceptable solution
        if (std::fabs(func(root.first + (root.second - root.first) / 2.0)) <= ROOT_ABS_TOLERANCE) { // Solution within tolerance?
            done = true;
        }
        else if (std::fabs(func(root.first)) <= ROOT_ABS_TOLERANCE) {                               // Solution within tolerance at endpoint 1?
            root.second = root.first;
            done = true;
        }
        else if (std::fabs(func(root.second)) <= ROOT_ABS_TOLERANCE) {                              // Solution within tolerance at endpoint 2?
            root.first = root.second;
            done = true;
        }
        else {                                                                                      // No - try again
            factorFrac /= 2.0;
            factor      = 1.0 + factorFrac;
            tries++;
            if (tries > ADAPTIVE_RV_MAX_TRIES || std::fabs(factor - 1.0) <= ROOT_ABS_TOLERANCE) {   // Too many tries, or step size 0.0?
                root.first  = -1.0;
                root.second = -1.0;
                SHOW_WARN_STATIC(ERROR::TOO_MANY_RV_TRIES);                                         // Show warning
                done = true;
            }
        }
    }

    // Midway between brackets is our result
    return root.first + (root.second - root.first) / 2.0;
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
 * double CalculateZAMSVrot(const double p_MZAMS)
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS equatorial rotational velocity (km s^-1)
 */
double CalculateZAMSVrot(const double p_MZAMS) {

    double vRot = 0.0;

    switch (OPTIONS->RotationalVelocityDistribution()) {                // Which rotational velocity distribution?

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

            if (p_MZAMS >= 16.0) {
                vRot = CalculateRotationalVelocityOStar_Ramirez2013();
                vRot = std::max(vRot, 0.0);                             // Set to no rotation if no positive solution found; warning already raised
            }
            else if (p_MZAMS >= 2.0) {
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

        default:                                                        // Unexpected distribution
            // The only way this can happen is if someone added a ROTATIONAL_VELOCITY_DISTRIBUTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNEXPECTED_VROT_DISTRIBUTION);    // Throw error
    }

    return vRot;
}


} // namespace astro
