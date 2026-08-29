#pragma once

#include "utils.h"

namespace sampling {

/*
 * Sampling functions
 * 
 * The "sampling" namespace gathers sampling functions that have utility across the
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

            double   InverseSampleFromPowerLaw(const double p_Power, const double p_Max, const double p_Min);
GNU_CONST   double   CalculateCDFKroupa(const double p_Mass, const double p_Max, const double p_Min);
            Dbl_DblT DrawKickDirection(const KICK_DIRECTION_DISTRIBUTION p_Distribution, const double p_Power);
COMPAS_DIAG double   InverseSampleFromTabulatedCDF(const double p_Y, const std::map<double, double>& p_Table);
            double   SampleEccentricity(const ECCENTRICITY_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min);
COMPAS_DIAG double   SampleFromTabulatedCDF(const double p_X, const std::map<double, double>& p_Table);
            double   SampleInitialMass(const INITIAL_MASS_FUNCTION p_Distribution, const double p_Max, const double p_Min, const double p_Power);
            double   SampleMassRatio(const MASS_RATIO_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min);
            double   SampleMetallicity(const METALLICITY_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min);
            double   SampleNSBirthMagneticField();
            double   SampleNSBirthSpinPeriod();
            double   SampleOrbitalPeriod(const ORBITAL_PERIOD_DISTRIBUTION p_Distribution, const double p_Max, const double p_Min);

 std::tuple<double, bool> SampleSemiMajorAxis(
    const SEMI_MAJOR_AXIS_DISTRIBUTION p_Distribution, 
    const double                       p_aMax, 
    const double                       p_aMin, 
    const double                       p_pMax = 0.0, 
    const double                       p_pMin = 0.0, 
    const double                       p_Mass1 = 0.0, 
    const double                       p_Mass2 = 0.0
);


} // namespace sampling
