#pragma once

#include "constants.h"
#include "BaseStar.h"
#include "MainSequence.h"

class BaseStar;
class MainSequence;

class MS_gt_07: virtual public BaseStar, public MainSequence {




public:

    // Constructors
    MS_gt_07() {};
    MS_gt_07(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) { 
        if (p_Initialise) Initialise();
    }

    // Initialise


private:

    // private member functions - alphabetically

    void Initialise() {
        CalculateTimescales(); // Initialise timescales // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // Age for MS_GT_07 is carried over from CH stars switching to MS after spinning down, so not set to 0.0 here
        
        //// JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // Initialise core mass, luminosity, radius, and temperature if Brcek core mass prescription is used
        // Only do this once - this should not be done if a CH star spins down and becomes a MS star (when using CHE_MODE::PESSIMISTIC)
        //
        if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK &&                                                          // Brcek core mass prescription?
            utils::Compare(MZAMS(), BRCEK_LOWER_MASS_LIMIT) >= 0                         &&                                                             // ZAMS mass >= BRCEK_LOWER_MASS_LIMIT?
            Time() <= 0.0) {                                                                                                                            // star not yet aged past creation?
                                                                                                                                                        // yes - initialise
            const double cnoCoreMass = CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024(MZAMS());                                                      // compute once
            m_InterimState.SetCNOProcessedCoreMass(cnoCoreMass);
            m_InterimState.SetMSCoreMass(cnoCoreMass);
            m_InterimState.SetLuminosity (MainSequence::CalculateLuminosityOnPhase(Age(), Mass0(), BaseStar::CalculateLuminosityAtZAMS(Mass0())));
            m_InterimState.SetRadius     (MainSequence::CalculateRadius());
            m_InterimState.SetTemperature(BaseStar::CalculateTemperatureOnPhase(Luminosity(), Radius()));
        }
    }


protected:

    // Member functions (not getters or setters)
    //
    // VIRTUAL FUNCTIONS may be (are expected to be) overridden by derived classes.
    // When overriding virtual functions in a derived class, use the "override" attribute.
    //
    // NON-VIRTUAL FUNCTIONS should not be overridden (declared separately) by derived classes.
    // While it is legal in C++ to declare the same (non-virtual) function in multiple classes,
    // (aka "shadowing", or "hiding"), we discourage it.  Non-virtual functions are statically
    // bound, and as such, especially with indirection, may not produce expected results.

    //////////////////////////////////////////////////
    //   ENVELOPE                                   //
    //////////////////////////////////////////////////

    COMPAS_PURE ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const override;

};


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ENVELOPE                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

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
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Envelope type (ENVELOPE::{RADIATIVE, CONVECTIVE, REMNANT})
 */
ENVELOPE MS_gt_07::DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const {
    
    ENVELOPE envType = ENVELOPE::CONVECTIVE;                                    // Default return value
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                             // Which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::CONVECTIVE_MASS_FRACTION:             // CONVECTIVE_MASS_FRACTION
            // Envelope is convective when the convective mass exceeds specified fraction of the envelope mass
            double mEnv;
            std::tie(mEnv, std::ignore) = CalculateConvectiveEnvelopeMass();
            envType = (mEnv / (p_Mass - p_CoreMass)) > OPTIONS->ConvectiveEnvelopeMassThreshold() ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:                    // FIXED_TEMPERATURE
            // Envelope is radiative if temperature exceeds specified threshold, otherwise convective
            envType = (p_Temperature * TSOL) > OPTIONS->ConvectiveEnvelopeTemperatureThreshold() ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:                               // HURLEY
            // There is some convective envelope for stars below 1.25 solar masses according to
            // Hurley et al. 2002, eq 36, but we simplify
            envType = p_Mass < 1.25 ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:                               // COMPAS LEGACY
            envType = ENVELOPE::RADIATIVE;                                      // Always radiative
            break;

        default:                                                                // Unexpected prescription
            // The only way this can happen is if the ENVELOPE_STATE_PRESCRIPTION served to this
            // function is not accounted for in this code.  We should not default here, with or
            // without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR_STATIC(ERROR::UNEXPECTED_ENVELOPE_STATE_PRESCRIPTION);  // Throw error          
    }
    
    return envType;
}









/// MS_gt_07_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//



    // MT rejuvenation
    double  CalculateMTRejuvenationFactor() override;                                                       // body in MS_gt_07.cpp

    // Critical mass ratio
    double  CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const override;       // body inline below
    double  CalculateCriticalMassRatio_Hurley2002() const override { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_GT_07; }



/*
 * CalculateCriticalMassRatio_Claeys2014
 *
 * @brief
 * Calculate the critical mass ratio, per Claeys et al. 2014
 * 
 * The critical mass ratio indicates whether the mass transfer is unstable.
 *
 * See Claeys et al. 2014, de Mink et al. 2013, and Ge et al. 2010, 2015, and 2020, for discussions.
 *
 * Assumes this star is the donor.
 * 
 * Critical mass ratio is defined as qCrit = Maccretor / Mdonor.
 * 
 *
 * double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const
 *
 * @param       p_AccretorIsDegenerate          Boolean indicating if accretor is degenerate
 * @return                                      Critical mass ratio for unstable MT 
 */
inline double MS_gt_07_Constituent::CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {                                                                                                          
    return p_AccretorIsDegenerate
            ? OPTIONS->MassTransferCriticalMassRatioMSHighMassDegenerateAccretor()       // degenerate accretor
            : OPTIONS->MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor();   // non-degenerate accretor
}

