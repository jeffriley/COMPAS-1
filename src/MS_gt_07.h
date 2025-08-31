#ifndef __MS_gt_07_h__
#define __MS_gt_07_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"
#include "MainSequence.h"

class BaseStar;
class MainSequence;

class MS_gt_07: virtual public BaseStar, public MainSequence {

public:

    MS_gt_07() { m_StellarType = STELLAR_TYPE::MS_GT_07; };
    
    MS_gt_07(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::MS_GT_07;                                                                                                         // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                 // Initialise if required
    }

    MS_gt_07* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        MS_gt_07* clone = new MS_gt_07(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static MS_gt_07* Clone(MS_gt_07& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        MS_gt_07* clone = new MS_gt_07(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                          // Initialise timescales
        // Age for MS_GT_07 is carried over from CH stars switching to MS after spinning down, so not set to 0.0 here
        
        // Initialise core mass, luminosity, radius, and temperature if Brcek core mass prescription is used
        // Only do this once - this should not be done if a CH star spins down and becomes a MS star (when using CHE_MODE::PESSIMISTIC)
        if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK &&                                                          // Brcek core mass prescription?
            utils::Compare(m_MZAMS, BRCEK_LOWER_MASS_LIMIT) >= 0                         &&                                                             // ZAMS mass >= BRCEK_LOWER_MASS_LIMIT?
            m_Time <= 0.0) {                                                                                                                            // star not yet aged past creation?
                                                                                                                                                        // yes - initialise
            m_InitialMainSequenceCoreMass = MainSequence::CalculateInitialMainSequenceCoreMass(m_MZAMS, m_InitialHeliumAbundance);
            m_MainSequenceCoreMass        = m_InitialMainSequenceCoreMass;
            m_Luminosity                  = MainSequence::CalculateLuminosityOnPhase(m_Age, m_Mass0, BaseStar::CalculateLuminosityAtZAMS(m_Mass0));
            m_Radius                      = MainSequence::CalculateRadiusOnPhase(m_Mass, m_Tau, m_RZAMS0);
            m_Temperature                 = BaseStar::CalculateTemperatureOnPhase_Static(m_Luminosity, m_Radius);
        }
    }

    // sets state zero for the star
    // Since this is a main sequence star, these attributes are ZAMS attributes
    // this function is guaranteed to be called only once during the lifetime of
    // a star - it is called from the Star class constructor, after the initial
    // switch to the required starting stellar type: it is not called when a star
    // switches stellar type (via Star::SwitchTo())
    void InitialiseState(const unsigned long int p_RandomSeed,
                         const double            p_Metallicity,
                         const double            p_Mass,
                         const KickParametersT   p_KickParameters,
                         const double            p_RotationalVelocity) {

        StarStateT state;

        state.stellarType              = STELLAR_TYPE::MS_GT_07;     // this stellar type

        state.evolutionStatus          = EVOLUTION_STATUS::CONTINUE; // initially
        state.error                    = ERROR::NONE;                // initially

        // set state attributes supplied as parameters
        state.randomSeed               = p_RandomSeed;
        state.metallicity              = p_Metallicity;
        state.mass                     = p_Mass;
        state.kickParameters           = p_KickParameters;

        // set, or calculate where possible/necessary, remaining state attributes

        state.time                     = 0.0;
        state.age                      = 0.0;
        state.tau                      = 0.0;

        state.dt                       = 0.0;

        state.dMdt                     = 0.0;
        state.dominantMassLossType     = MASS_LOSS_TYPE::NONE;

        state.luminosity               = CalculateZAMSLuminosity_Tout_1996(p_Mass, GLOBALS->LuminosityCoefficients());
        state.radius                   = CalculateZAMSRadius_Tout_1996(p_Mass, GLOBALS->RadiusCoefficients());

        state.heliumAbundance          = GLOBALS->ZAMSHeliumAbundance();
        state.heliumAbundanceCore      = heliumAbundance;
        state.heliumAbundanceSurface   = heliumAbundance;

        state.hydrogenAbundance        = GLOBALS->ZAMSHydrogenAbundance();
        state.hydrogenAbundanceCore    = hydrogenAbundance;
        state.hydrogenAbundanceSurface = hydrogenAbundance;

        state.coreMass                 = 0.0;

        // if the Brcek MS core mass prescription was specified by the user, and if the ZAMS mass
        // of the star indicates it is in the Brcek regime, initialise core mass, and recalculate
        // luminosity and radius
        if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK &&  // Brcek core mass prescription?
            utils::Compare(state.mass, BRCEK_LOWER_MASS_LIMIT) >= 0) {                          // yes - in Brcek regime?
                                                                                                // yes - (re)calculate as appropriate
            state.luminosity = MainSequence::CalculateLuminosityOnPhase(0.0, state.mass, state.luminosity);
            state.radius     = MainSequence::CalculateRadiusOnPhase(state.mass, 0.0, state.radius);

            // calculate ZAMS core mass per Shikauchi et al. 2024, eq A3
            DBL_VECTOR fmixCoeffs = std::get<1>(SHIKAUCHI_COEFFICIENTS);
            state.coreMass = state.mass * (fmixCoeffs[0] + fmixCoeffs[1] * std::exp(-state.mass / fmixCoeffs[2]));
        }

        state.coreMassMS          = state.coreMass;
        state.coreMassCO          = 0.0;
        state.coreMassHe          = 0.0;

        state.temperature         = CalculateTemperatureOnPhase(state.luminosity, state.radius);

        state.massEffective       = 0.0;
        state.radiusEffective     = 0.0;
        state.luminosityEffective = 0.0;

        state.angularFrequency = p_RotationalFrequency >= 0.0                           // valid rotational frequency passed in?
        ? _2_PI * p_RotationalFrequency                     // yes - convert from cycles/yr to rad/yr and use it
        : CalculateZAMSAngularFrequency(m_MZAMS, m_RZAMS);  // no - calculate it
m_AngularMomentum                          = CalculateMomentOfInertiaAU() * m_OmegaZAMS;


        // state (state zero) initialised - push to history
        // This will be the first entry in the state history, and will stay
        // locked at the start of the history for the lifetime of the star
        (void)m_StateHistory.Push(state);
    }

    // member functions - alphabetically

    double      CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const ;
    double      CalculateCriticalMassRatioHurleyHjellmingWebbink() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_GT_07; }
    double      CalculateMassLossRateHurley();
    double      CalculateMassTransferRejuvenationFactor();

    ENVELOPE    DetermineEnvelopeType() const;
};

#endif // __MS_gt_07_h__
