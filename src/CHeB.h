#ifndef __CHeB_h__
#define __CHeB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "FGB.h"


class BaseStar;
class FGB;

class CHeB: virtual public BaseStar, public FGB {

public:

    CHeB() { m_StellarType = STELLAR_TYPE::CORE_HELIUM_BURNING; };
    
    CHeB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), FGB(p_BaseStar, false) {
        m_StellarType = STELLAR_TYPE::CORE_HELIUM_BURNING;                                                                                                      // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                         // Initialise if required
    }

    CHeB* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        CHeB* clone = new CHeB(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static CHeB* Clone(CHeB& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        CHeB* clone = new CHeB(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


    // member functions - alphabetically

double CalculateMinLuminosityOnPhase_Hurley2000(const double      p_Mass,
                                                const double      p_MHeF,
                                                const double      p_MFGB,
                                                const double      p_Alpha1,
                                                const DBL_VECTOR& p_bN) const


    static double CalculateMinimumRadiusOnPhase_Static(const double      p_Mass,
                                                       const double      p_CoreMass,
                                                       const double      p_Alpha1,
                                                       const double      p_MHeF,
                                                       const double      p_MFGB,
                                                       const double      p_LuminosityOnPhase,
                                                       const DBL_VECTOR &p_BnCoefficients);


protected:

    void Initialise() {

        m_StellarType = STELLAR_TYPE::CORE_HELIUM_BURNING;                                                                                                      // Set stellar type
        CalculateTimescales();                                                                                                                                  // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tHeI)];                                                                                                // Set age appropriately
        m_MinimumLuminosityOnPhase = CalculateMinimumLuminosityOnPhase(massCutoffs(MHeF), m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), GLOBALS->HurleyBCoefficients());    // Calculate once, not many


CalculateMinLuminosityOnPhase_Hurley2000(const double      p_Mass,
                                                const double      p_MHeF,
                                                const double      p_MFGB,
                                                const double      p_Alpha1,
                                                const DBL_VECTOR& p_bN);




        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically

    double          CalculateBluePhaseFBL(const double p_Mass);


///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


GNU_CONST inline double CalculateCoreMass_Hurley2000() const override {
    return CalculateCoreMass_Hurley2000(m_StateHistory.CurrentState.Mass(), m_StateHistory.CurrentState.Tau());
}
GNU_CONST double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau) const;


GNU_CONST inline double CalculateCOCoreMass_Hurley2000() const override { return 0.0; } // CHeB stars have no CO core

GNU_CONST inline double CalculateTau_Hurley2000() const override {
    return CalculateTau_Hurley2000(m_StateHistory.CurrentState.Age(),
                                   m_StateHistory.CurrentState.Timescales(tHeI),
                                   m_StateHistory.CurrentState.Timescales(tHe));
}
GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const double p_tHeI, const double p_tHe) const;



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCoreMassAtPhaseEnd_Hurley2000() const override {
    return CalculateCoreMass_Hurley2000();                                              // per Hurley SSE code `hrdiag.f` lines 259-265
}







double CalculateRadiusOnPhase_Hurley(const double p_Mass, const double p_Luminosity, const double p_Tau) const;
double CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return GiantBranch::CalculateRadiusOnPhase(p_Mass, p_Luminosity); }
double CalculateRadiusOnPhase() const { return CalculateRadiusOnPhase(m_Mass, m_Luminosity, m_Tau); }




GNU_PURE  double CalculateCriticalMassRatio_Hurley2002() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_GT_07; }

double CalculateHeCoreMassAtPhaseEnd() const { return m_StateHistory.CurrentState.CoreMass(); }





GNU_CONST double CalculateHAbundanceCoreOnPhase(const double p_Tau, const double p_InitialHAbundance) const override { return 0.0; }
GNU_PURE  double CalculateHAbundanceSurfaceOnPhase(const double p_Tau, const double p_InitialHAbundance) const override { return GLOBALS->ZAMSHAbundance(); }
    
GNU_CONST double CalculateHeAbundanceCoreOnPhase(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance = 0.0) const override;
GNU_PURE  double CalculateHeAbundanceSurfaceOnPhase(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance) const override { return GLOBALS->ZAMSHeAbundance(); }                      
 



GNU_CONST double CalculateCELambda_Dewi(const double p_Mass,
                                        const double p_Radius,
                                        const double p_Luminosity,
                                        const double p_RZAMS,
                                        const double p_CoreMass) const;



    double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const;
    double          CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;

    double          CalculateLifetimeOnBluePhase(const double p_Mass);





double CalculatePhaseLifetime() const override;   
GNU_CONST double CalculatePhaseLifetime_Hurley2000(const double      p_Mass,
                                                   const double      p_CoreMass,
                                                   const double      p_tBGB,
                                                   const double      p_MHeF,
                                                   const DBL_VECTOR& p_bN,
                                                   const double      p_Alpha4) const;
                                                                      





    double          CalculateLuminosityAtBluePhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtBluePhaseStart(const double p_Mass) const;

    double          CalculateLuminosityAtPhaseEnd() const                       { return CalculateLuminosityAtBAGB(m_Mass0); }
    double          CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau) const;
    double          CalculateLuminosityOnPhase() const                          { return CalculateLuminosityOnPhase(m_Mass0, m_Tau); }

    double          CalculateRadiusAtBluePhaseEnd(const double p_Mass) const;
    double          CalculateRadiusAtBluePhaseStart(const double p_Mass) const;

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const;
    double          CalculateRadiusAtPhaseEnd() const                           { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }





static double CalculateMinRadiusOnPhase_Hurley2000_Static(const double      p_Mass,
                                                          const double      p_CoreMass,
                                                          const double      p_Alpha1,
                                                          const double      p_MHeF,
                                                          const double      p_MFGB,
                                                          const double      p_MinLuminosity,
                                                          const DBL_VECTOR& p_bN);



    double          CalculateRadiusRho(const double p_Mass, const double p_Tau) const;

    double          CalculateRemnantLuminosity() const;
    double          CalculateRemnantRadius() const;






GNU PURE  DBL_VECTOR CalculateTimescales_Hurley(const double p_Mass, const DBL_VECTOR& p_MassCutoffs, const DBL_VECTOR& p_Timescales) const;

    void            CalculateTimescales()                                       { CalculateTimescales(m_Mass0, m_Timescales); }                                 // Use class member variables

    double          ChooseTimestep(const double p_Time) const;

    ENVELOPE        DetermineEnvelopeType() const;

    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                        { return !ShouldEvolveOnPhase(); }                                              // Phase ends when age at or after He Burning
    bool            IsSupernova() const                                         { return false; }                                                               // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash() {  }                                                                                                                   // NO-OP

    bool            ShouldEnvelopeBeExpelledByPulsations() const { return ( OPTIONS->ExpelConvectiveEnvelopeAboveLuminosityThreshold() && DetermineEnvelopeType() == ENVELOPE::CONVECTIVE && utils::Compare( log10(m_Luminosity/m_Mass), OPTIONS->LuminosityToMassThreshold() ) >= 0 ) ; }                             // Envelope of convective star with luminosity to mass ratio beyond threshold should be expelled
    bool            ShouldEvolveOnPhase() const;
    bool            ShouldSkipPhase() const                                     { return false; }                                                               // Never skip CHeB phase

};




////// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                             LIFETIME / AGE FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculatePhaseLifetime
 *
 * @brief
 * Calculate lifetime of this phase of the evolution of the star (Core Helium Burning, tHe).
 *
 * Calls relevant lifetime function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculatePhaseLifetime()
 *
 * @return                                      CHeB lifetime, tHe (Myr)
 */
inline double CHeB::CalculatePhaseLifetime() const { 

    double lifetime;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            lifetime = CalculatePhaseLifetime_Hurley2000(m_StateHistory.CurrentState.Mass(),
                                                         m_StateHistory.CurrentState.CoreMass(),
                                                         GLOBALS->MassCutoffs(MHEF),
                                                         m_StateHistory.CurrentState.Timescales(tBGB),
                                                         GLOBALS->HurleyBCoefficients(),
                                                         GLOBALS->HurleyAlpha4());
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

    return lifetime;
}


/*
 * CalculatePhaseLifetime_Hurley2000
 *
 * @brief
 * Calculate the lifetime of Core Helium Burning, tHe, per Hurley at al. 2000, eq 57
 *
 *
 * double CalculatePhaseLifetime_Hurley2000(const double      p_Mass,
 *                                          const double      p_CoreMass,
 *                                          const double      p_tBGB,
 *                                          const double      p_MHeF,
 *                                          const DBL_VECTOR& p_bN,
 *                                          const double      p_Alpha4) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_tBGB                          Lifetime to Base of Giant Branch, tBGB (Myr)
 * @param       p_bN                            Hurley b(n) coefficients
 * @param       p_Alpha4                        Hurley alpha4 constant
 * @return                                      CHeB lifetime, tHe (Myr)
 */
GNU_CONST inline double CHeB::CalculatePhaseLifetime_Hurley2000(const double      p_Mass,
                                                                const double      p_CoreMass,
                                                                const double      p_MHeF,
                                                                const double      p_tBGB,
                                                                const DBL_VECTOR& p_bN,
                                                                const double      p_Alpha4) const {
    double tHe;

    if (p_Mass < p_MHeF) {
        const double tHeMS = HeMS::CalculatePhaseLifetime_Hurley2000_Static(p_CoreMass);
        tHe = (p_bN[39] + ((tHeMS - p_bN[39]) * PPOW((1.0 - (p_Mass / p_MHeF)), p_bN[40]))) * (1.0 + (p_Alpha4 * exp(15.0 * (p_Mass - p_MHeF))));
    }
    else {
        const double m5 = utils::IntPow(p_Mass, 5);
        tHe = p_tBGB * (((p_bN[41] * PPOW(p_Mass, p_bN[42])) + (p_bN[43] * m5)) / (p_bN[44] + m5));
    }

    return tHe;
}


/*
 * CalculateTau_Hurley2000
 *
 * @brief
 * Calculate the CHeB-relative age (fractional Core Helium Burning age) of the star,
 * per Hurley et al. 2000, just before eq 59.
 *
 * 
 * double CalculateTau_Hurley2000(const double p_Age, const double p_tHeI, const double p_tHe) const
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tHeI                          Time to helium ignition (Myr)
 * @param       p_tHe                           Time to helium burning (Myr)
 * @return                                      CHeB-relative age, [0, 1]
 */
GNU_CONST inline double CHeB::CalculateTau_Hurley2000(const double p_Age, const double p_tHeI, const double p_tHe) const {
    return std::max(0.0, std::min(1.0, (p_Age - p_tHeI) / p_tHe));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                ABUNDANCE FUNCTIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateHeAbundanceCoreOnPhase
 *
 * @brief
 * Calculate the helium abundance in the core of the star, given the phase relative age of
 * the star.  Currently just a simple linear model from the initial helium abundance to the
 * maximum helium abundance (assuming that all hydrogen is converted to helium). 
 * 
 * Should one day be updated to match detailed models.
 *
 * double CalculateHeAbundanceCoreOnPhase(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Tau                           CHeB relative age of the star [0, 1]
 * @param       p_InitialHeAbundance            Initial helium abundance of the star (not used here)
 * @return                                      Helium abundance in the core of the star
 */
GNU_CONST inline double CHeB::CalculateHeAbundanceCoreOnPhase(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance) const {
    return (1.0 - p_Metallicity) * (1.0 - p_Tau);
}





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 RADIUS FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRemnantRadius
 *
 * @brief
 * Calculate radius of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass)
 *
 * Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantRadius(const double p_Mass, const double p_Tau)
 *
 * @return                                      Radius of remnant core in Rsol
 */
GNU_CONST inline double CHeB::CalculateRemnantRadius() const {
    return HeMS::CalculateRadius_Hurley2000_Static(m_CoreMass, CalculateTauOnPhase());
}




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  MASS FUNCTIONS                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCoreMass_Hurley2000
 *
 * @brief
 * Calculate the core mass between Helium Ignition (HeI) and the Base of the
 * Asymptotic Giant Branch (BAGB), per Hurley et al. 2000, eq 67.
 *
 *
 * double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           CHeB-relative age of the star [0.0, 1.0]
 * @return                                      EAGB core mass (Msol)
 */
GNU_CONST inline double CHeB::CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau) const {
    // should become HeMS star - He mass clamped to total mass
    return std::min(((1.0 - p_Tau) * CalculateCoreMassAtHeI_Hurley2000(p_Mass)) + (p_Tau * CalculateCoreMassAtBAGB(p_Mass)), p_Mass);
}









#endif // __CHeB_h__
