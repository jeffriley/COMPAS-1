#include "BH.h"

/*
 * 
 * Calculate the gravitational mass due to neutrino mass loss
 *
 * double CalculateNeutrinoMassLoss_Static(const double p_BaryonicMass)
 *
 * @param   [IN]    p_BaryonicMass              Baryonic remnant mass in Msol
 * @return                                      Gravitational mass of remnant in Msol
 */
double BH::CalculateNeutrinoMassLoss_Static(const double p_BaryonicMass) {

    double gravitationalMass = 0.0;

    switch (OPTIONS->NeutrinoMassLossAssumptionBH()) {                                                  // which assumption?

        case NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION:                                           // FIXED FRACTION
            gravitationalMass = p_BaryonicMass * (1.0 - OPTIONS->NeutrinoMassLossValueBH());
            break;

        case NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS:                                               // FIXED MASS
            gravitationalMass = p_BaryonicMass - OPTIONS->NeutrinoMassLossValueBH();
            break;
    
        default:                                                                                        // unknown prescription
            // the only way this can happen is if someone added a NEUTRINO_MASS_LOSS_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose an assumption this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing assumption or, if the missing
            // assumption is superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNKNOWN_NEUTRINO_MASS_LOSS_PRESCRIPTION);                         // throw error
    }

    return gravitationalMass;
}


/*
 * CalculateSNkickWeighting_Static
 *
 * @brief
 * Calculate the weighting to be applied to the kick magnitude for a Black Hole, BH,
 * based on the user-supplied black hole kixks mode (via the option `--black-hole-kicks-mode`).
 *
 * Option values are:
 *
 *    FALLBACK: Black holes receive a kick down-weighted by the amount of mass falling back onto them
 *    FULL    : Black holes receive the same kicks as neutron stars (weighting = 1.0)
 *    REDUCED : Black holes receive the same momentum kick as a neutron star, but down-weighted by the black hole mass
 *    ZERO    : Black holes receive zero natal kick (weighting = 0.0)
 *
 *
 * static double CalculateSNkickWeighting_Static(const double p_Mass, const double p_FallbackFraction)
 *
 * @param       p_Mass                          Mass of the remnant (BH) (Msol)
 * @param       p_FallbackFraction              Fraction of mass that falls back onto the proto-compact object [0, 1]
 * @return                                      Kick magnitude
 */
COMPAS_PURE double BH::CalculateSNkickWeighting_Static(const double p_Mass, const double p_FallbackFraction) {

    double weighting;

    switch (OPTIONS->BlackHoleKicksMode()) {                                                    // which BH kicks mode?

        case BLACK_HOLE_KICKS_MODE::FALLBACK: weighting = 1.0 - p_FallbackFraction;   break;    // using the so-called 'fallback' mode for BH kicks
        case BLACK_HOLE_KICKS_MODE::FULL    : weighting = 1.0;                        break;    // full kick - no adjustment necessary
        case BLACK_HOLE_KICKS_MODE::REDUCED : weighting = NEUTRON_STAR_MASS / p_Mass; break;    // kick is reduced by the ratio of the neutron star mass to the black hole mass
        case BLACK_HOLE_KICKS_MODE::ZERO    : weighting = 0.0;                        break;    // no kick

        default:                                                                                // unknown mode
            // the only way this can happen is if someone added a BLACK_HOLE_KICKS_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNKNOWN_BH_KICK_MODE);                                    // throw error
    }

    return weighting;
}
