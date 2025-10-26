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

    HeHG* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HeHG* clone = new HeHG(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static HeHG* Clone(HeHG& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HeHG* clone = new HeHG(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


    // member functions

    static void CalculateGBparams_Static(const double p_Mass0, const double p_Mass, const double p_LogMetallicityXi, const DBL_VECTOR &p_MassCutoffs, const DBL_VECTOR &p_AnCoefficients, const DBL_VECTOR &p_BnCoefficients, DBL_VECTOR &p_GBparams);
 
    double  CalculateRemnantRadius() const;

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


    // member functions - aphabetically
    
            double          CalculateConvectiveCoreMass() const                                                     { return m_CoreMass; }

GNU_CONST inline double CalculateCoreMassAtBAGB_Hurley2000(const double p_Mass, const DBL_VECTOR& p_bN) const override { return m_Mass0; } // McBAGB = M0 (Hurley et al. 2000, discussion just before eq 89)


GNU_CONST static double CalculateConvectiveCoreRadius_Hurley2000_Static(const double p_Radius, const double p_Tau, const double p_CoreMass);
GNU_CONST static double CalculateConvectiveCoreRadius_Static(const double p_Radius, const double p_Tau);
GNU_CONST inline double CalculateConvectiveCoreRadius(const double p_Radius, const double p_Tau) const override { return CalculateConvectiveCoreRadius_Static(p_Radius, p_Tau); }

GNU_CONST static double CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass);
GNU_CONST static double CalculateRemnantRadius_Static(const double p_CoreMass);
GNU_CONST inline double CalculateRemnantRadius(const double p_CoreMass) const override { return CalculateRemnantRadius_Static(p_CoreMass); }


            double          CalculateCoreMassAtPhaseEnd() const                                                     { return m_CoreMass; }                                                  // NO-OP
            double          CalculateCoreMassOnPhase() const                                                        { return m_COCoreMass; }                                                // Mc(HeMS) = McCOMass




            void            CalculateGBparams(const double p_Mass, DBL_VECTOR &p_GBparams);
            void            CalculateGBparams()                                                                     { CalculateGBparams(m_Mass0, m_GBparams); }                             // Use class member variables

            double          CalculateHeCoreMassAtPhaseEnd() const                                                   { return CalculateHeCoreMassOnPhase(); }                                // Same as on phase
            double          CalculateHeCoreMassOnPhase() const                                                      { return m_Mass; }                                                      // NO-OP
    
            double          CalculateHeliumAbundanceCoreAtPhaseEnd() const                                          { return 0.0; }



GNU_CONST inline double CalculateHAbundanceCoreOnPhase(const double p_Tau, const double p_InitialHAbundance) const override { return 0.0; }
GNU_CONST inline double CalculateHeAbundanceCoreOnPhase(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance = 0.0) const override { return 0.0; }
            


            double          CalculateHydrogenAbundanceCoreAtPhaseEnd() const                                        { return 0.0; }




            double          CalculateLuminosityOnPhase() const;
            double          CalculateLuminosityAtPhaseEnd() const                                                   { return m_Luminosity; }                                                // NO-OP

            double          CalculateMomentOfInertia() const                                                        { return GiantBranch::CalculateMomentOfInertia(); }

            double          CalculatePerturbationMu() const;
            double          CalculatePerturbationMuAtPhaseEnd() const                                               { return m_Mu; }                                                        // NO-OP

            double          CalculateRadiusAtPhaseEnd() const                                                       { return m_Radius; }                                                    // NO-OP
   
            double          CalculateRadiusOnMassChange(double p_dM)                                                { return CalculateRadiusOnPhase(m_Mass + p_dM, m_Luminosity); }
            double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }
            double          CalculateRadiusOnPhase(double p_Mass, double p_Luminosity) const;


            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity) const;
            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                          { return CalculateRadiusAndStellarTypeOnPhase(m_Mass, m_Luminosity); }







///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inline double CalculateCOCoreMass_Hurley2000() const override {
    return HeGB::CalculateCoreMass_Hurley2000_Static(m_StateHistory.CurrentState.MassEffectiveInitial(),
                                                     m_StateHistory.CurrentState.Age(),
                                                     m_StateHistory.CurrentState.GBparams(),
                                                     m_StateHistory.CurrentState.Timescales(tHeMS));
}

double CalculateTau_Hurley2000() const override { return 0.0; } // Tau (relative age) is not used for HeHG stars in Hurley et al. 2000, so we return 0.0







 



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inline double CalculateCOCoreMassAtPhaseEnd_Hurley2000() const override { return CalculateCOCoreMass_Hurley2000(); }





            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return m_Temperature; }                                               // NO-OP
            double          CalculateTemperatureAtPhaseEnd() const                                                  { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }      // Use class member variables

            double          CalculateThermalMassLossRate() const                                                    { return GiantBranch::CalculateThermalMassLossRate(); }                 // Skip HeMS

            void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
            void            CalculateTimescales()                                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables

            double          ChooseTimestep(const double p_Time) const;

            ENVELOPE        DetermineEnvelopeType() const;

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
//                             LIFETIME / AGE FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  MASS FUNCTIONS                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 RADIUS FUNCTIONS                                  //
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
 * @param       p_Tau                           HeHG-relative age of the star [0.0, 1.0]
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
 * @param       p_Tau                           HeHG-relative age of the star [0.0, 1.0]
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
 * Calculate the radius of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass), per Hurley et al. 2000, at the 
 * end of section 6 (after eq 105).
 *
 *
 * static double CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass)
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Radius of remnant core (Rsol)
 */
GNU_CONST static inline double HeHG::CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass) {
    return WhiteDwarfs::CalculateRadius_Hurley2000_Static(p_CoreMass);
}


/*
 * CalculateRemnantRadius_Static
 *
 * @brief
 * Calculate the radius of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass).
 *
 * Calls relevant radius function based on the evolutionary mode given in program options.
 * 
 *
 * static double CalculateRemnantRadius_Static(const double p_CoreMass)
 * 
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Radius of remnant core (Rsol)
 */
GNU_CONST static inline double HeHG::CalculateRemnantRadius_Static(const double p_CoreMass) { 

    double radius;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            radius = CalculateRemnantRadius_Hurley2000_Static(p_CoreMass);
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







/// HeHG_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


class HeHG_Constituent: virtual public BinaryConstituentStar, public HeHG {

public:


protected:

};

    
            double          CalculateMTRejuvenationFactor()                                               { return 1.0; }


            double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const;
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
