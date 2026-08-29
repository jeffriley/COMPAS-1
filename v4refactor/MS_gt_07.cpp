#include "MS_gt_07.h"
#include "Star.h" 













/// MS_gt_07_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


/*
 * CalculateMTRejuvenationFactor
 *
 * @brief
 * Calculate rejuvenation factor for stellar age, based on the user-specified MT_REJUVENATION_PRESCRIPTION
 *
 *
 * double CalculateMTRejuvenationFactor(const double p_Mass, consr double p_MassPrev) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_MassPrev                      Mass of the star at the previous timestep (Msol)
 * @return                                      Rejuvenation factor
 */
double MS_gt_07_Constituent::CalculateMTRejuvenationFactor() {

    const double mass     = m_Star->Mass();
    const double massPrev = m_Star->MassPrev();

    double fRej;

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {                                  // which prescription?

        case MT_REJUVENATION_PRESCRIPTION::HURLEY:                                              // HURLEY
            fRej = 1.0;                                                                         // Hurley et al. 2000
            break;

        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                                           // STARTRACK
            // StarTrack 2008 prescription
            // Belczynski et al. 2007, section 5.6 (http://arxiv.org/pdf/astro-ph/0511811v3.pdf)
            fRej = mass < massPrev ? 1.0 : massPrev / mass;                                     // rejuvenation factor is unity for mass losing stars
            break;

        default:                                                                                // unknown prescription
            // the only way this can happen is if someone added a MT_REJUVENATION_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION);                           // throw error
    }

    return fRej;
}
