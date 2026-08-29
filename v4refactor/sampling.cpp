#include "sampling.h"

namespace sampling {

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

    double result = 0.0;                                                        // Default return value

    // We use utils::Compare() here because if p_Power as passed is the result of a numerical
    // computation and is very close to -1.0, but not actually -1.0 (e.g. -0.9999999), a simple
    // test of (p_Power == -1.0) would fail, and the code path would default to the general case,
    // and there the result becomes neumerically unstable as (a + 1) approaches 0.0.
    if (utils::Compare(p_Power, -1.0) == 0) {                                   // p_Power very close to -1.0?
        result = std::exp(RAND->Random() * std::log(p_Max / p_Min)) * p_Min;    // Yes - special case
    }
    else {                                                                      // No - general case
        const double powerPlus1     = p_Power + 1.0;
        const double min_PowerPlus1 = PPOW(p_Min, powerPlus1);

        result = PPOW((RAND->Random() * (PPOW(p_Max, powerPlus1) - min_PowerPlus1) + min_PowerPlus1), 1.0 / powerPlus1);
    }

    return result;
}


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
double CalculateCDFKroupa(const double p_Mass, const double p_Max, const double p_Min) {

    if ((p_Mass < p_Min) || (!(p_Mass < p_Max))) return 0.0;    // Return 0.0 if mass is out of bounds of the IMF
    
    double vCDF = 0.0;                                          // Default return value

    if (!(p_Min > KROUPA_BREAK_1) && (p_Max > KROUPA_BREAK_1) && !(p_Max > KROUPA_BREAK_2)) {

        const double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
        const double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (PPOW(p_Max, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);

        const double C1 = 1.0 / (term1 + term2);
        const double C2 = C1 * KROUPA_BREAK_1_POWER_1_2;

        if (!(p_Mass < p_Min) && (p_Mass < KROUPA_BREAK_1)) {
            vCDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_1) - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
        }
        else if (!(p_Mass < KROUPA_BREAK_1) && (p_Mass < KROUPA_BREAK_2)) {
            vCDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1)) +
                   ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);
        }
    }
    else if (!(p_Min > KROUPA_BREAK_1) && (p_Max >  KROUPA_BREAK_2)) {

        const double term1 = ONE_OVER_KROUPA_POWER_1_PLUS1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
        const double term2 = ONE_OVER_KROUPA_POWER_2_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2);
        const double term3 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_1_POWER_1_2 * KROUPA_BREAK_2_POWER_2_3 * (PPOW(p_Max, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);

        const double C1 = 1.0 / (term1 + term2 + term3);
        const double C2 = C1 * KROUPA_BREAK_1_POWER_1_2;
        const double C3 = C2 * KROUPA_BREAK_2_POWER_2_3;

        if (!(p_Mass < p_Min) && (p_Mass < KROUPA_BREAK_1)) {
            vCDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_1) - PPOW(p_Min, KROUPA_POWER_PLUS1_1));
        }
        else if (!(p_Mass < KROUPA_BREAK_1) && (p_Mass < KROUPA_BREAK_2)) {
            vCDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1)) +
                   ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_2) - KROUPA_BREAK_1_PLUS1_2);
        }
        else if (!(p_Mass < KROUPA_BREAK_2) && (p_Mass < p_Max)) {
            vCDF = ONE_OVER_KROUPA_POWER_1_PLUS1 * C1 * (KROUPA_BREAK_1_PLUS1_1 - PPOW(p_Min, KROUPA_POWER_PLUS1_1)) +
                   ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (KROUPA_BREAK_2_PLUS1_2 - KROUPA_BREAK_1_PLUS1_2) +
                   ONE_OVER_KROUPA_POWER_3_PLUS1 * C3 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);
        }
    }
    else if ((p_Min > KROUPA_BREAK_1) && !(p_Min > KROUPA_BREAK_2) && (p_Max > KROUPA_BREAK_2)) {

        const double term1 = ONE_OVER_KROUPA_POWER_2_PLUS1 * (KROUPA_BREAK_2_PLUS1_2 - PPOW(p_Min, KROUPA_POWER_PLUS1_2));
        const double term2 = ONE_OVER_KROUPA_POWER_3_PLUS1 * KROUPA_BREAK_2_POWER_2_3 * (PPOW(p_Max, KROUPA_POWER_PLUS1_3) - KROUPA_BREAK_2_PLUS1_3);

        const double C2 = 1.0 / (term1 + term2);
        const double C3 = C2 * KROUPA_BREAK_2_POWER_2_3;

        if (!(p_Mass < p_Min) && (p_Mass < KROUPA_BREAK_2)) {
            vCDF = ONE_OVER_KROUPA_POWER_2_PLUS1 * C2 * (PPOW(p_Mass, KROUPA_POWER_PLUS1_2) - PPOW(p_Min, KROUPA_POWER_PLUS1_2));
        }
        else if (!(p_Mass < KROUPA_BREAK_2) && (p_Mass < p_Max)) {
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
Dbl_DblT DrawKickDirection(const KICK_DIRECTION_DISTRIBUTION p_Distribution, const double p_Power) {

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
            theta = RAND->Random() < 0.5 ? M_PI_2 - theta : theta - M_PI_2;         // Close, or equal, to pi/2. Choose up or down randomly
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

            THROW_ERROR_STATIC(ERROR::UNEXPECTED_KICK_DIRECTION_DISTRIBUTION);      // Throw error
    }

    return std::make_tuple(theta, phi);
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
double InverseSampleFromTabulatedCDF(const double p_Y, const std::map<double, double>& p_Table) {

    if (p_Table.empty()) THROW_ERROR_STATIC(ERROR::EMPTY_MAP);                          // Bounds check - map must not be empty
    if (p_Y < 0.0 || p_Y >= 1.0) THROW_ERROR_STATIC(ERROR::OUT_OF_BOUNDS);              // Bounds check - p_Y must in [0.0, 1.0)

    double xInterp = 0.0;                                                               // default return value

    std::map<double, double>::const_iterator iter;                                      // Iterator

    const double yMin = p_Table.begin()->second;                                        // Minimum y value (map is ordered)
    if (yMin < 0.0) THROW_ERROR_STATIC(ERROR::OUT_OF_BOUNDS);                           // Bounds check - min y must be >= 0.0

    const double yMax = p_Table.rbegin()->second;                                       // Maximum y value
    if (yMax > 1.0) THROW_ERROR_STATIC(ERROR::OUT_OF_BOUNDS);                           // Bounds check - max y must be <= 1.0

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

    double eccentricity = 0.0;                                                      // Default return value

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

        case ECCENTRICITY_DISTRIBUTION::GELLER_2013:                                // GELLER_2013
            // M35 eccentricity distribution from Geller, Hurley and Mathieu 2013
            // Gaussian with mean 0.38 and sigma 0.23
            // http://iopscience.iop.org/article/10.1088/0004-6256/145/1/8/pdf
            // Sampling function taken from binpop.f in NBODY6

            // Loop until we find a value between min and max - shouldn't loop forever...
            do {
                eccentricity = 0.23 * std::sqrt(-2.0 * log(RAND->Random())) * cos(_2_PI * RAND->Random()) + 0.38;
            } while (eccentricity < p_Min || eccentricity > p_Max);
            break;

        case ECCENTRICITY_DISTRIBUTION::DUQUENNOYMAYOR1991:                         // DUQUENNOYMAYOR1991
            // Eccentricity distribution from Duquennoy & Mayor (1991)
            // http://adsabs.harvard.edu/abs/1991A%26A...248..485D
            // Sampling function taken from binpop.f in NBODY6

            // Loop until we find a value between min and max - shouldn't loop forever...
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

            THROW_ERROR_STATIC(ERROR::UNEXPECTED_ECCENTRICITY_DISTRIBUTION);        // Throw error
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
double SampleFromTabulatedCDF(const double p_X, const std::map<double, double>& p_Table) {

    double yInterp = 0.0;                                                               // Default return value

    std::map<double, double>::const_iterator iter;                                      // Iterator

    if (p_Table.empty()) THROW_ERROR_STATIC(ERROR::EMPTY_MAP);                          // Bounds check - map must not be empty

    const double xMin = p_Table.begin()->first;                                         // Minimum x value (map is ordered)
    const double xMax = p_Table.rbegin()->first;                                        // Maximum x value

    if (p_X < xMin || p_X > xMax) THROW_ERROR_STATIC(ERROR::OUT_OF_BOUNDS);             // Bounds check - p_X must in [xMin, xMax)

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
 * Returns p_Min if p_Min == p_Max.
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
double SampleInitialMass(const INITIAL_MASS_FUNCTION p_Distribution, const double p_Max, const double p_Min, const double p_Power) {

    // We use utils::Compare() here - if p_Min and p_Max are equal within our floating-point
    // comparison tolerance, we are unlikely to actually sample a value between them.
    if (utils::Compare(p_Min, p_Max) == 0) return p_Min;

    double mass = 0.0;                                                                  // Default return value

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
            else if (p_Min <= KROUPA_BREAK_1 && p_Max > KROUPA_BREAK_2) {
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

            THROW_ERROR_STATIC(ERROR::UNEXPECTED_INITIAL_MASS_FUNCTION);                // Throw error
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

    double q = 0.0;                                                         // Default return value

    switch (p_Distribution) {                                               // Which q distribution?

        case MASS_RATIO_DISTRIBUTION::FLAT:                                 // FLAT
            q = utils::InverseSampleFromPowerLaw(0.0, p_Max, p_Min);
            break;

        case MASS_RATIO_DISTRIBUTION::DUQUENNOYMAYOR1991:                   // DUQUENNOYMAYOR1991                                              
            // Mass ratio distribution from Duquennoy & Mayor (1991) (http://adsabs.harvard.edu/abs/1991A%26A...248..485D)

            // loop until we find a value between min and max - shouldn't loop forever...
            do {
                q = 0.42 * std::sqrt(-2.0 * log(RAND->Random())) * cos(_2_PI * RAND->Random()) + 0.23;
            } while (q < p_Min || q > p_Max);
            break;

        case MASS_RATIO_DISTRIBUTION::SANA2012:                             // SANA2021                                                    
            // Sana et al. 2012 (http://science.sciencemag.org/content/sci/337/6093/444.full.pdf) distribution of eccentricities.
            // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
            // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

            q = utils::InverseSampleFromPowerLaw(-0.1, p_Max, p_Min);       // de Mink and Belczynski use min = 0.1, max = 1.0
            break;

        default:                                                            // Unexpected distribution
            // The only way this can happen are if someone added a mass ratio distribution and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNEXPECTED_MASS_RATIO_DISTRIBUTION);  // Throw error
    }

    return std::min(std::max(p_Min, q), p_Max);                             // Clamp to [min, max]
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

    double metallicity = 0.0;                                               // Default return value

    switch (p_Distribution) {                                               // Which Z distribution?

        case METALLICITY_DISTRIBUTION::ZSOLAR:                              // ZSOLAR
            // All systems have Z = ZSOLAR (ZSOL_ASPLUND)
            metallicity = ZSOL_ASPLUND;
            break;

        case METALLICITY_DISTRIBUTION::LOGUNIFORM: {                        // LOGUNIFORM
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

        default:                                                            // Unexpected distribution
            // The only way this can happen are if someone added a metallicity distribution and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNEXPECTED_METALLICITY_DISTRIBUTION); // Throw error
    }

    return metallicity;
}


/*
 * SampleNSBirthMagneticField
 *
 * @brief
 * Sample (log10 of) the magnetic field (in G) for a pulsar at birth, per user-supplied
 * birth magnetic field distribution.
 * (see option --pulsar-birth-magnetic-field-distribution)
 *
 *
 * double SampleNSBirthMagneticField()
 *
 * @return                                      log10 of the birth magnetic field in G
 */
double SampleNSBirthMagneticField() {

	double log10B = 0.0;                                                                    // Default return value

    switch (OPTIONS->PulsarBirthMagneticFieldDistribution()) {                              // Which distribution?

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FLATINLOG: {                         // FLAT IN LOG
            // Distribution from Oslowski et al. 2011 https://arxiv.org/abs/0903.3538

            const double maximum = OPTIONS->PulsarBirthMagneticFieldDistributionMax();
            const double minimum = OPTIONS->PulsarBirthMagneticFieldDistributionMin();
            
            log10B = minimum + (RAND->Random() * (maximum - minimum));

            } break;

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::UNIFORM: {                           // UNIFORM
            // Flat distribution used in Kiel et al. 2008 https://arxiv.org/abs/0805.0059
            // (log10B0min = 11, log10B0max = 13.5 see section 3.4 and Table 1.)    
            const double maximum = PPOW(10.0, OPTIONS->PulsarBirthMagneticFieldDistributionMax());
            const double minimum = PPOW(10.0, OPTIONS->PulsarBirthMagneticFieldDistributionMin());

            log10B = std::log10(minimum + (RAND->Random() * (maximum - minimum)));
            } break;

        case PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::LOGNORMAL: {                         // LOG NORMAL
            // Distribution from Faucher-Giguere and Kaspi 2006 https://arxiv.org/abs/astro-ph/0512585
            const double mean  = OPTIONS->PulsarBirthMagneticFieldDistributionMean();
            const double sigma = OPTIONS->PulsarBirthMagneticFieldDistributionSigma();

            log10B = RAND->RandomGaussian(sigma) + mean;

            // Ensure magnetic field is always larger than or equal to the user-supplied value
            // of the program option --pulsar-minimum-magnetic-field
            SizeT iterations = 0;
            do { 
                log10B = RAND->RandomGaussian(sigma) + mean;
            } while (iterations++ < PULSAR_MAG_ITERATIONS && OPTIONS->PulsarLog10MinimumMagneticField() > log10B);

            if (iterations >= PULSAR_MAG_ITERATIONS) THROW_ERROR_STATIC(ERROR::TOO_MANY_PULSAR_MAG_ITERATIONS);
            } break;

        default:                                                                            // Unexpected prescription
            // The only way this can happen is if the PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION served
            // to this function is not accounted for in this code.  We should not default here, with
            // or without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR_STATIC(ERROR::UNEXPECTED_PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION); // Throw error
    }

    return log10B;
}


/*
 * SampleNSBirthSpinPeriod
 *
 * @brief
 * Sample the spin period of a Pulsar at birth per user-supplied birth spin period distribution.
 * (see option --pulsar-birth-spin-period-distribution)

 *
 * double SampleNSBirthSpinPeriod()
 *
 * @return                                      Pulsar birth spin period (s)
 */
double SampleNSBirthSpinPeriod() {

	double pSpin = 0.0;                                                                     // Default return value

    switch (OPTIONS->PulsarBirthSpinPeriodDistribution()) {                                 // Which distribution?

        case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::UNIFORM: {                              // UNIFORM
            // Distribution between minimum and maximum value as in Oslowski et al. 2011
            // https://arxiv.org/abs/0903.3538 (default Pmin = and Pmax = )
            // and also Kiel et al. 2008 https://arxiv.org/abs/0805.0059
            // (default Pmin = 10 ms and Pmax 100 ms, section 3.4)
            const double maximum = OPTIONS->PulsarBirthSpinPeriodDistributionMax();
            const double minimum = OPTIONS->PulsarBirthSpinPeriodDistributionMin();

            pSpin = minimum + (RAND->Random() * (maximum - minimum));
            } break;

        case PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::NORMAL: {                               // NORMAL
            // Distribution from Faucher-Giguere and Kaspi 2006
            // https://arxiv.org/abs/astro-ph/0512585
            const double mean  = OPTIONS->PulsarBirthSpinPeriodDistributionMean();
            const double sigma = OPTIONS->PulsarBirthSpinPeriodDistributionSigma();

            // Ensure magnetic field is always larger than 0.0
            SizeT iterations = 0;
            do { 
                pSpin = RAND->RandomGaussian(sigma) + mean;
            } while (iterations++ < PULSAR_SPIN_ITERATIONS && !(pSpin > 0.0));
            if (iterations >= PULSAR_SPIN_ITERATIONS) THROW_ERROR_STATIC(ERROR::TOO_MANY_PULSAR_SPIN_ITERATIONS);

            } break;

        default:                                                                            // Unexpected prescription
            // The only way this can happen is if the PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION served
            // to this function is not accounted for in this code.  We should not default here, with
            // or without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR_STATIC(ERROR::UNEXPECTED_PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION);    // Throw error
    }

    return pSpin * SECONDS_IN_MS;
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

    double orbitalPeriod = 0.0;                                                 // Default return value

    switch (p_Distribution) {                                                   // Which distribution?

        case ORBITAL_PERIOD_DISTRIBUTION::FLATINLOG:                            // FLATINLOG

            orbitalPeriod = utils::InverseSampleFromPowerLaw(-1.0, p_Max, p_Min);
            break;

        default:                                                                // Unexpected distribution
            // The only way this can happen are if someone added an orbital period distribution and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNEXPECTED_ORBITAL_PERIOD_DISTRIBUTION);  // Throw error
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
    const double                       p_pMax, 
    const double                       p_pMin, 
    const double                       p_Mass1, 
    const double                       p_Mass2) {

    // We use utils::Compare() here - if p_Min and p_Max are equal within our floating-point
    // comparison tolerance, we are unlikely to actually sample a value between them.
    if (utils::Compare(p_aMin, p_aMax) == 0) return {p_aMin, true};

    double semiMajorAxis = 0.0;                                                                         // Default return value
    bool   convergence   = true;                                                                        // Default return value

    switch (p_Distribution) {                                                                           // Which a distribution?

        case SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG:                                                   // FLATINLOG
            semiMajorAxis = utils::InverseSampleFromPowerLaw(-1.0, p_aMax, p_aMin);
            break;

        case SEMI_MAJOR_AXIS_DISTRIBUTION::DUQUENNOYMAYOR1991: {                                        // DUQUENNOYMAYOR1991
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
                semiMajorAxis  = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(p_Mass1, p_Mass2, p);    // Convert period (days) to semi-major axis (AU)
            } while (samples++ < SEMI_MAJOR_AXIS_SAMPLES && (semiMajorAxis < p_aMin || semiMajorAxis > p_aMax));

            if (samples >= SEMI_MAJOR_AXIS_SAMPLES) convergence = false;                                // Signal no convergence to caller
            } break;

        case SEMI_MAJOR_AXIS_DISTRIBUTION::SANA2012: {                                                  // SANA2012
            // Sana et al. 2012
            // http://science.sciencemag.org/content/sci/337/6093/444.full.pdf
            // distribution of semi-major axes. Sana et al. fit for the orbital period, which we sample in here, before returning the semi major axis
            // Taken from table S3 in http://science.sciencemag.org/content/sci/suppl/2012/07/25/337.6093.444.DC1/1223344.Sana.SM.pdf
            // See also de Mink and Belczynski 2015 http://arxiv.org/pdf/1506.03573v2.pdf

            const double logpMin = p_pMin > 1.0 ? log(p_pMin) : 0.0;                                    // Lower bound
            const double logpMax = p_pMax > 1.0 ? log(p_pMax) : 0.0;                                    // Upper bound

            const double p = std::exp(utils::InverseSampleFromPowerLaw(-0.55, logpMax, logpMin));       // Draw period from their distribution
            semiMajorAxis  = utils::ConvertPeriodInDaysToSemiMajorAxisInAU(p_Mass1, p_Mass2, p);        // Convert period (days) to semi-major axis (AU)
            } break;

        default:                                                                                        // Unexpected distribution
            // The only way this can happen are if someone added a semi-major axis distribution and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a distribution this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing distribution or, if the missing
            // distribution is superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNEXPECTED_SEMI_MAJOR_AXIS_DISTRIBUTION);                         // Throw error
    }

    return std::make_tuple(semiMajorAxis, convergence);
}

} // namespace sampling
