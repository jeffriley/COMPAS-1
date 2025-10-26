#include "MS_gt_07.h"


/*
 * CalculateMLRate_Hurley2000
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, per Hurley et al. 2000.
 *
 *
 * MASS_LOSS_T CalculateMLRate_Hurley2000(const double p_Mass,
 *                                        const double p_Radius,
 *                                        const double p_Luminosity,
 *                                        const double p_PerturbationMu,
 *                                        const double p_ZscaledHurley,
 *                                        const double p_WRfactor) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @param       p_ZscaledHurley                 Z inversely scaled by Hurley ZSOL (Z / ZSOL_HURLEY)
 * @param       p_WRfactor                      WR mass loss factor
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (coule be MASS_LOSS_TYPE::NONE)
 */
GNU_CONST MASS_LOSS_T MS_gt_07::CalculateMLRate_Hurley2000(const double p_Mass,
                                                           const double p_Radius,
                                                           const double p_Luminosity,
                                                           const double p_PerturbationMu,
                                                           const double p_ZscaledHurley,
                                                           const double p_WRfactor) const {
    const double dMdt = CalculateMLRate_NieuwenhuijzenDeJager(p_Mass, p_Radius, p_Luminosity, p_ZscaledHurley);
    return std::make_tuple(dMdt, (dMdt > 0.0 ? MASS_LOSS_TYPE::GB : MASS_LOSS_TYPE::NONE));
}


/*
 * DetermineEnvelopeType
 *
 * @brief
 * Determine the star's envelope type, based on the user-specified ENVELOPE_STATE_PRESCRIPTION.
 * 
 *
 * double DetermineEnvelopeType(const double p_Mass, const double p_Temeparture) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Envelope type (ENVELOPE::{RADIATIVE, CONVECTIVE, REMNANT})
 */
ENVELOPE MS_gt_07::DetermineEnvelopeType(const double p_Mass, const double p_Temeparture) const {
    
    ENVELOPE envelopeType = ENVELOPE::RADIATIVE;                                                    // default envelope type
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                                                 // which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:                                                   // COMPAS LEGACY
            envelopeType = ENVELOPE::RADIATIVE;                                                     // always radiative
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:                                                   // HURLEY
            // there is some convective envelope for stars below 1.25 solar masses
            // according to Hurley et al. 2002, eq 36, but we simplify
            envelopeType = p_Mass < 1.25 ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:                                        // FIXED_TEMPERATURE
            // envelope is radiative if temperature exceeds fixed threshold, otherwise convective
            envelopeType = (p_Temeparture * TSOL) > OPTIONS->ConvectiveEnvelopeTemperatureThreshold() ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::CONVECTIVE_MASS_FRACTION:                                 // CONVECTIVE_MASS_FRACTION
            // envelope is labelled convective when the convective mass exceeds a fixed fraction of the envelope mass
            double convectiveEnvelopeMass, convectiveEnvelopeMassMax;
            std::tie(convectiveEnvelopeMass, convectiveEnvelopeMassMax) = CalculateConvectiveEnvelopeMass();
            envelopeType = (convectiveEnvelopeMass / (p_Mass - m_CoreMass)) > OPTIONS->ConvectiveEnvelopeMassThreshold() ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;

        default:                                                                                    // unknown prescription
            // the only way this can happen is if someone added an ENVELOPE_STATE_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION);                                // throw error            
    }
    
    return envelopeType;
}















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
