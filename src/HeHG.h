#ifndef __HeHG_h__
#define __HeHG_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "HeMS.h"


class BaseStar;
class HeMS;

class HeHG: virtual public BaseStar, public HeMS {

public:

    HeHG() { m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP; };
    
    HeHG(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeMS(p_BaseStar, false) {
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;                                                                                                                    // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                                     // Initialise if required
    }


    // member functions

    static void CalculateGBparams_Static(const double p_Mass0, const double p_Mass, const double p_LogMetallicityXi, const DBL_VECTOR &p_MassCutoffs, const DBL_VECTOR &p_AnCoefficients, const DBL_VECTOR &p_BnCoefficients, DBL_VECTOR &p_GBparams);
 


protected:

    void Initialise() {
        m_Tau = 0.0;                                                                                      // Start of phase
        CalculateTimescales();                                                                            // Initialise timescales
        // Age for HeHG is calculated before switching -
        // can get here via EvolveOneTimestep() and ResolveEnvelopeLoss(),
        // and Age is calculated differently in those cases
        
        // Update stellar properties at start of HeHG phase (since core definition changes)
        CalculateGBparams();

        EvolveOnPhase(0.0);
    }


    // on phase member functions - aphabetically
    

DBL_VECTOR CalculateTimescales_Hurley2000() const override {
    return CalculateTimescales_Hurley2000(
        m_StateHistory.CurrentState.MassEffectiveInitial(),
        m_StateHistory.CurrentState.GBparams(),
        m_StateHistory.CurrentState.TimeScales()
    );
}
COMPAS_PURE DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;

DBL_VECTOR CalculateGBparams_Hurley2000() const override {
    return CalculateGBparams_Hurley2000(m_StateHistory.CurrentState.MassEffectiveInitial(), m_StateHistory.CurrentState.GBparams());
}
COMPAS_PURE virtual DBL_VECTOR CalculateGBparams_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const;





            double          CalculateConvectiveCoreMass() const                                                     { return m_CoreMass; }

GNU_CONST inline double CalculateCoreMassAtBAGB_Hurley2000(const double p_Mass, const DBL_VECTOR& p_bN) const override { return m_Mass0; } // McBAGB = M0 (Hurley et al. 2000, discussion just before eq 89)


GNU_CONST static double CalculateConvectiveCoreRadius_Hurley2000_Static(const double p_Radius, const double p_Tau, const double p_CoreMass);
GNU_CONST static double CalculateConvectiveCoreRadius_Static(const double p_Radius, const double p_Tau);
GNU_CONST inline double CalculateConvectiveCoreRadius(const double p_Radius, const double p_Tau) const override { return CalculateConvectiveCoreRadius_Static(p_Radius, p_Tau); }


inline double CalculateRemnantRadius_Hurley2000() const override {
    return CalculateRemnantRadius_Static(m_StateHistory.CurrentState.CoreMass());
}
GNU_CONST static double CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass);


            double          CalculateCoreMassAtPhaseEnd() const                                                     { return m_CoreMass; }                                                  // NO-OP
            double          CalculateCoreMassOnPhase() const                                                        { return m_COCoreMass; }                                                // Mc(HeMS) = McCOMass





                                                     // NO-OP
    
            double          CalculateHeliumAbundanceCoreAtPhaseEnd() const                                          { return 0.0; }



            


            double          CalculateHydrogenAbundanceCoreAtPhaseEnd() const                                        { return 0.0; }




            double          CalculateLuminosityOnPhase() const;
            double          CalculateLuminosityAtPhaseEnd() const                                                   { return m_Luminosity; }                                                // NO-OP

            double          CalculateMomentOfInertia() const                                                        { return GiantBranch::CalculateMomentOfInertia(); }


GNU_CONST double CalculateHurleyPerturbationMu() const;
inline double CalculateHurleyPerturbationMuAtPhaseEnd() const { return m_Mu; }


            double          CalculateRadiusAtPhaseEnd() const                                                       { return m_Radius; }                                                    // NO-OP
   
            double          CalculateRadiusOnMassChange(double p_dM)                                                { return CalculateRadiusOnPhase(m_Mass + p_dM, m_Luminosity); }
            double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }
            double          CalculateRadiusOnPhase(double p_Mass, double p_Luminosity) const;


            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity) const;
            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                          { return CalculateRadiusAndStellarTypeOnPhase(m_Mass, m_Luminosity); }







///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inline double CalculateCOCoreMass_Hurley2000() const override {
    return HeGB::CalculateCoreMass_Hurley2000_Static(
        m_StateHistory.CurrentState.MassEffectiveInitial(),
        m_StateHistory.CurrentState.Age(),
        m_StateHistory.CurrentState.GBparams(),
        m_StateHistory.CurrentState.Timescales(tHeMS)
    );
}

inline double CalculateHeCoreMass() const override { return m_StateHistory.CurrentState.Mass(); } // McHe = Mc for HeHG stars

GNU_CONST inline double HeHG::CalculateHurleyPerturbationMu(const double p_Mass, const double p_CoreMass) const;

GNU_CONST inline double CalculateTau_Hurley2000() const override { return 0.0; } // Tau (relative age) is not used for HeHG stars in Hurley et al. 2000, so we return 0.0


    inline double CalculateLuminosity_Hurley2000() const override {
        return HeGB::CalculateLuminosity_Hurley2000_Static(
            m_StateHistory.CurrentState.CoreMass(),
            m_StateHistory.CurrentState.GBparams(static_cast<int>(HURLEY_GBP::B)),
            m_StateHistory.CurrentState.GBparams(static_cast<int>(HURLEY_GBP::D))
        );
    }

    inline double CalculateRadius_Hurley2000() const override {
        return HeGB::CalculateRadius_Hurley2000_Static(m_StateHistory.CurrentState.Mass(), m_StateHistory.CurrentState.Luminosity());
    }

    GNU_CONST inline double CalculateHAbundanceCore(const double p_Tau) const override { return 0.0; } // No hydrogen in the core for HeHG stars
    GNU_CONST inline double CalculateHeAbundanceCore(const double p_Tau) const override { return 0.0; } // No helium in the core for HeHG stars

 



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inline double CalculateCOCoreMassAtPhaseEnd_Hurley2000() const override { return CalculateCOCoreMass_Hurley2000(); }

inline double CalculateHeCoreMassAtPhaseEnd() const override { return m_StateHistory.CurrentState.Mass(); } // McHe = Mc for HeHG stars



            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return m_Temperature; }                                               // NO-OP
            double          CalculateTemperatureAtPhaseEnd() const                                                  { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }      // Use class member variables

inline double CalculateMLrateThermal() const override { return GiantBranch::CalculateMLrateThermal(); } // Skip HeMS


            double          ChooseTimestep(const double p_Time) const;


COMPAS_PURE ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const override;


            STELLAR_TYPE    EvolveToNextPhase();

            bool            IsEndOfPhase() const                                                                    { return !ShouldEvolveOnPhase(); }
            bool            IsSupernova() const;
            double          CalculateInitialSupernovaMass() const;
    
            void            PerturbLuminosityAndRadius()                                                            { GiantBranch::PerturbLuminosityAndRadius(); }                          // NO-OP

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
            void            ResolveHeliumFlash() { }                                                                                                                                        // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                                   { return m_StellarType; }                                               // NO-OP

            bool            ShouldEnvelopeBeExpelledByPulsations() const                                            { return CHeB::ShouldEnvelopeBeExpelledByPulsations(); }                // Envelope of convective star with luminosity to mass ratio beyond threshold should be expelled
            bool            ShouldEvolveOnPhase() const;
            bool            ShouldSkipPhase() const                                                                 { return false; }                                                       // Never skip HeMS phase

            void            UpdateAgeAfterMassLoss()                                                                { GiantBranch::UpdateAgeAfterMassLoss(); }                              // No action for He giants
            double          CalculateInitialMass()                                                                  { return GiantBranch::CalculateInitialMass(); }                         // Skip HeMS
};




/////// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         AGE / LIFETIME / TAU / TIMESCALES                         //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           MISC. CONSTANTS / PARAMETERS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateHurleyPerturbationMu
 *
 * @brief
 * Calculate the Hurley perturbation parameter, mu, per Hurley et al. 2000, eq 98 (&89)
 *
 *
 * double CalculateHurleyPerturbationMu(const double p_Mass, const double p_CoreMass) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Perturbation parameter, mu
 */
GNU_CONST inline double HeHG::CalculateHurleyPerturbationMu(const double p_Mass, const double p_CoreMass) const {
    // Hurley et al. 2000, eq 89 - clamped to ensure McMax >= p_CoreMass
    const double McMax = std::max(std::min(((1.45 * p_Mass) - 0.31), p_Mass), p_Coremass);
    return 5.0 * (McMax - p_CoreMass) / McMax; // Hurley et al. 2000, eq 98
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateConvectiveCoreRadius_Hurley2000_Static
 *
 * @brief
 * Calculate the convective core radius, per Hurley et al. 2000, sec 6 (last para)
 *
 * 
 * static double CalculateConvectiveCoreRadius_Hurley2000_Static(const double p_Radius, const double p_Tau, const double p_CoreMass)
 * 
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Convective core radius (Rsol)         
 */
GNU_CONST static inline double HeHG::CalculateConvectiveCoreRadius_Hurley2000_Static(const double p_Radius, const double p_Tau, const double p_CoreMass) {
    return std::min(5.0 * HeHG::CalculateRemnantRadius_Hurley2000_Static(p_CoreMass), p_Radius);
}


/*
 * CalculateConvectiveCoreRadius_Static
 *
 * @brief
 * Calculate the convective core radius.
 *
 * Calls relevant radius function based on the evolutionary mode given in program options.
 * 
 *
 * static double CalculateConvectiveCoreRadius_Static(const double p_Radius, const double p_Tau)
 * 
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      Convective core radius (Rsol) 
 */
GNU_CONST static inline double HeHG::CalculateConvectiveCoreRadius_Static(const double p_Radius, const double p_Tau) { 

    double radius;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            radius = CalculateConvectiveCoreRadius_Hurley2000_Static(p_Radius, p_Tau);
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return radius;
}


/*
 * CalculateRemnantRadius_Hurley2000_Static
 *
 * @brief
 * Calculate the radius of the remnant the star would become if it lost all of its envelope
 * immediately (i.e. M = Mc), per Hurley et al. 2000, at the end of section 6 (after eq 105).
 *
 *
 * static double CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass)
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Remnant core radius (Rsol)
 */
GNU_CONST static inline double HeHG::CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass) {
    return WhiteDwarfs::CalculateRadius_Hurley2000_Static(p_CoreMass);
}







/// HeHG_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


class HeHG_Constituent: virtual public BinaryConstituentStar, public HeHG {

public:


protected:

};

    
            double          CalculateMTRejuvenationFactor()                                               { return 1.0; }


            double          CalculateLambdaNanjingStarTrack(const double p_Mass) const;
            double          CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const { return CalculateLambdaNanjingStarTrack(0.0, 0.0); }        // 0.0 are dummy values that are not used




double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const;

double CalculateCriticalMassRatio_Hurley2002() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_HE_GIANT; }


            double          CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription)                  { return GiantBranch::CalculateZetaAdiabatic_ByEnvelopeType(p_ZetaPrescription); } // Calculate Zetas as for other giant stars (HeMS stars were an exception)
    
            double          CalculateZetaEquilibrium()                                                              { return -std::numeric_limits<double>::infinity(); }                     // Nuclear timescale MT should be impossible from HG stars that evolve on a thermal timescale


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
inline double HeHG_Constituent::CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {                                                                                                          
    return p_AccretorIsDegenerate
            ? OPTIONS->MassTransferCriticalMassRatioHeHGDegenerateAccretor()        // degenerate accretor
            : OPTIONS->MassTransferCriticalMassRatioHeHGNonDegenerateAccretor();    // non-degenerate accretor
}


#endif // __HeHG_h__
