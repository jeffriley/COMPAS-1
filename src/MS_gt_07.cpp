#include "MS_gt_07.h"


/*
 * DetermineEnvelopeType
 *
 * @brief
 * Determine the star's envelope type, based on the user-specified ENVELOPE_STATE_PRESCRIPTION.
 * 
 *
 * ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Envelope type (ENVELOPE::{RADIATIVE, CONVECTIVE, REMNANT})
 */
ENVELOPE MS_gt_07::DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const {
    
    ENVELOPE envType;
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                                 // which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::CONVECTIVE_MASS_FRACTION:                 // CONVECTIVE_MASS_FRACTION
            // envelope is convective when the convective mass exceeds specified fraction of the envelope mass
            double mEnv;
            std::tie(mEnv, std::ignore) = CalculateConvectiveEnvelopeMass();
            envType = (mEnv / (p_Mass - p_CoreMass)) > OPTIONS->ConvectiveEnvelopeMassThreshold() ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:                        // FIXED_TEMPERATURE
            // envelope is radiative if temperature exceeds specified threshold, otherwise convective
            envType = (p_Temperature * TSOL) > OPTIONS->ConvectiveEnvelopeTemperatureThreshold() ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:                                   // HURLEY
            // there is some convective envelope for stars below 1.25 solar masses according to
            // Hurley et al. 2002, eq 36, but we simplify
            envType = p_Mass < 1.25 ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:                                   // COMPAS LEGACY
            envType = ENVELOPE::RADIATIVE;                                          // always radiative
            break;

        default:                                                                    // unknown prescription
            // the only way this can happen is if someone added an ENVELOPE_STATE_PRESCRIPTION and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION);                // throw error            
    }
    
    return envType;
}













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
double MS_gt_07_Constituent::CalculateMTRejuvenationFactor(const double p_Mass, consr double p_MassPrev) const {

    double fRej;

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {                                  // which prescription?

        case MT_REJUVENATION_PRESCRIPTION::HURLEY:                                              // HURLEY
            fRej = 1.0;                                                                         // Hurley et al. 2000
            break;

        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                                           // STARTRACK
            // StarTrack 2008 prescription
            // Belczynski et al. 2007, section 5.6 (http://arxiv.org/pdf/astro-ph/0511811v3.pdf)
            fRej = p_Mass < m_MassPrev ? 1.0 : p_MassPrev / p_Mass;                             // rejuvenation factor is unity for mass losing stars
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
