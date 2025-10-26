#ifndef __MainSequence_h__
#define __MainSequence_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"

class BaseStar;

class MainSequence: virtual public BaseStar {
    
    class StarState : public BaseStar::StarState {

    #define RET_VALUE(x) { if (x.has_value()) { return x.value(); } else { std::cerr << "\nState attribute has no value: program terminated\n"; utils::ShowStackTrace(); std::exit(1); }}
            
            
    private:

        std::optional<double> m_HeAbundanceCore;       
        std::optional<double> m_HeAbundanceSurface;

    public:
            
        StarState() {}
        
        // setters
        void SetHeAbundanceCore(const double p_HeAbundanceCore)       { m_HeAbundanceCore = p_HeAbundanceCore; }
        void SetHeAbundanceSurface(const double p_HeAbundanceSurface) { m_HeAbundanceSurface = p_HeAbundanceSurface; }


        // getters
        double HeAbundanceCore() const                                { RET_VALUE(m_HeAbundanceCore); }
        double HeAbundanceSurface() const                             { RET_VALUE(m_HeAbundanceSurface); }
        
        
        // has value
        bool HeAbundanceCore_HasValue() const                         { return m_HeAbundanceCore.has_value(); }
        bool HeAbundanceSurface_HasValue() const                      { return m_HeAbundanceSurface.has_value(); }

    #undef RET_VALUE
    };


public:

    MainSequence(){};

    MainSequence(const BaseStar& p_BaseStar) : BaseStar(p_BaseStar) {}

    MT_CASE DetermineMassTransferTypeAsDonor() const                                        { return MT_CASE::A; }                                                  // Always case A
    
    const std::tuple <DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> SHIKAUCHI_COEFFICIENTS = InterpolateShikauchiCoefficients(m_Metallicity);                                 // Interpolate Shikauchi coefficients for the given metallicity

protected:
    
    // member variables
    double          m_HeliumAbundanceOutsideCore      = m_InitialHeliumAbundance;                                                                                       // Helium abundance just outside the core, used for rejuvenation calculations


    // static member functions - alphabetically

    GNU_PURE  double CalculateRadiusAtPhaseEnd_Hurley_Static(const double p_Mass, const double p_RZAMS, const DBL_VECTOR& p_aN) const;

    GNU_CONST double CalculateLuminosityAtZAMS_Tout1996(const double p_MZAMS, const DBL_VECTOR& p_LCoeffs) const;

    GNU_CONST double CalculateRadiusAtZAMS_Tout1996(const double p_MZAMS, const DBL_VECTOR& p_RCoeffs) const;



    // member functions - alphabetically

    // constants, coefficients, etc.
    GNU_CONST double CalculateHurleyAlphaL(const double p_Mass, const DBL_VECTOR& p_aN, const DBL_VECTOR& p_LConsts) const;
    GNU_CONST double CalculateHurleyAlphaR(const double p_Mass, const DBL_VECTOR& p_aN, const DBL_VECTOR& p_RConsts) const;
    GNU_CONST double CalculateHurleyBetaL(const double p_Mass, const DBL_VECTOR& p_aN, const DBL_VECTOR& p_LConsts) const;
    GNU_CONST double CalculateHurleyBetaR(const double p_Mass, const DBL_VECTOR& p_aN, const DBL_VECTOR& p_RConsts) const;
    GNU_CONST double CalculateHurleyDeltaL(const double p_Mass, const double p_HookMass, const DBL_VECTOR& p_aN, const DBL_VECTOR& p_LConsts) const;
    GNU_CONST double CalculateHurleyDeltaR(const double p_Mass, const double p_HookMass, const DBL_VECTOR& p_aN, const DBL_VECTOR& p_RConsts) const;
    GNU_CONST double CalculateHurleyEta(const double p_Metallicity, const double p_Mass) const;
    GNU_CONST double CalculateHurleyGamma(const double p_Mass, const DBL_VECTOR& p_aN, const DBL_VECTOR& p_GammaConsts) const;



///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; }            // McCO = 0.0 for MS stars



    // luminosity
    GNU_CONST double CalculateLuminosity(const EVOLUTION_MODE p_Mode,
                                         const double         p_Metallicity,
                                         const double         p_Mass,
                                         const double         p_Time,
                                         const double         p_LZAMS,
                                         const DBL_VECTOR&    p_Timescales,
                                         const DBL_VECTOR&    p_aN,
                                         const DBL_VECTOR&    p_LConsts) const;

    GNU_CONST double CalculateLuminosityAtPhaseEnd_Hurley(const double p_Mass, const DBL_VECTOR& p_aN) const;

    GNU_CONST double CalculateLuminosity_Brcek(const double      p_Metallicity,
                                               const double      p_Mass, 
                                               const double      p_Tau,
                                               const double      p_MZAMS,
                                               const double      p_LZAMS,
                                               const double      p_CoreMass,
                                               const DBL_VECTOR& p_Timescales,
                                               const DBL_VECTOR& p_aN,
                                               const DBL_VECTOR& p_ShikauchiLCoeffs) const;

    double CalculateLuminosityOnPhase();

    GNU_CONST double CalculateLuminosityOnPhase_Hurley(const double      p_Metallicity,
                                                       const double      p_Mass,
                                                       const double      p_Time,
                                                       const double      p_LZAMS,
                                                       const DBL_VECTOR& p_Timescales,                                                      
                                                       const DBL_VECTOR& p_aN,
                                                       const DBL_VECTOR& p_LConsts) const;

    GNU_CONST double CalculateLuminosity_Shikauchi(const double p_CoreMass, const double p_HeAbundanceCore, const DBL_VECTOR& p_ShikauchiLCoeffs) const;



    // radius
    GNU_CONST double CalculateRadius(const EVOLUTION_MODE p_Mode, const double p_Mass, const double p_Tau, const double p_RZAMS, const DBL_VECTOR& p_aN) const;

    GNU_CONST double CalculateRadiusAtPhaseEnd(const EVOLUTION_MODE p_Mode, const double p_Mass, const double p_RZAMS, const DBL_VECTOR& p_aN) const;

    GNU_PURE  double CalculateRadiusAtPhaseEnd_Hurley(const double p_Mass, const double p_RZAMS, const DBL_VECTOR& p_aN) const {
        return CalculateRadiusAtPhaseEnd_Hurley_Static(p_Mass, p_RZAMS, p_aN);
    }

GNU_CONST double CalculateRadiusOnPhase_Hurley(const double p_Mass, const double p_Tau, const double p_RZAMS, const DBL_VECTOR& p_aN) const;

    double CalculateRadiusOnPhase_Hurley() const {
        // use globals and state values - read only
        return CalculateRadiusOnPhase_Hurley(m_StateHistory.CurrentState().Mass(), m_StateHistory.CurrentState().Tau(), m_StateHistory.ZAMSState().Radius(), GLOBALS->HurleyACoefficients());
    }


inline double CalculateEffectiveInitialMass_Hurley2000() const override { return m_StateHistory.CurrentState.Mass(); } // per Hurley et al. 2000, section 7.1




    double          CalculateConvectiveCoreMass() const;
    double          CalculateConvectiveCoreRadius() const;
    DBL_DBL         CalculateConvectiveEnvelopeMass() const;
    double          CalculateBetaL(const double p_Mass) const;
    double          CalculateBetaR(const double p_Mass) const;

    double          CalculateDeltaR(const double p_Mass) const;


    double          CalculateGamma(const double p_Mass) const;






 



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; }  // McCO = 0.0 for MS stars






    double          CalculateCoreMassAtPhaseEnd() const                                     { return (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::MANDEL) ? std::min(MainSequenceCoreMass(), m_Mass) : 0.0; }        // Accounts for minimal core mass built up prior to mass loss through mass transfer; core mass can't exceed total mass
    double          CalculateCoreMassOnPhase() const                                        { return 0.0; }                                                         // Mc(MS) = 0.0 (Hurley et al. 2000, just before eq 28)

    double          CalculateHeCoreMassAtPhaseEnd() const                                   { return CalculateCoreMassAtPhaseEnd(); }                               // Same as He core mass
    double          CalculateHeCoreMassOnPhase() const                                      { return 0.0; }                                                         // McHe(MS) = 0.0

    double          CalculateHeliumAbundanceCoreAtPhaseEnd() const                          { return CalculateHeliumAbundanceCoreOnPhase(); }


GNU_PURE??? double CalculateHAbundanceCoreOnPhase(const double p_Tau, const double p_InitialHAbundance) const override;                                
inline double CalculateHAbundanceSurfaceOnPhase(const double p_Tau, const double p_InitialHAbundance) const override { return m_StateHistory.CurrentState.HAbundanceSurface(); }

GNU_PURE??? double CalculateHeAbundanceCoreOnPhase(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance = 0.0) const override;
inline double CalculateHeAbundanceSurfaceOnPhase(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance) const override { return m_StateHistory.CurrentState.HeAbundanceSurface(); }

    


    double          CalculateHeliumAbundanceSurfaceAtPhaseEnd() const                       { return CalculateHeliumAbundanceSurfaceOnPhase(); }
    
    double          CalculateHydrogenAbundanceCoreAtPhaseEnd() const                        { return CalculateHydrogenAbundanceCoreOnPhase(); } 
    double          CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const;                                                          
    
    double          CalculateHydrogenAbundanceSurfaceAtPhaseEnd() const                     { return CalculateHydrogenAbundanceSurfaceOnPhase(); } 

    


virtual double CalculatePhaseLifeTime() const;
GNU_CONST double CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_ZetaHurley, const double p_TBGB, const DBL_VECTOR& p_aN) const;


double CalculateTau_Hurley2000() const override {
    return CalculateTau_Hurley2000(m_StateHistory.CurrentState.Age(), m_StateHistory.CurrentState.Timescales(tMS));
}
GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const double p_tMS) const;

    double          CalculateLuminosityAtPhaseEnd() const                                   { return CalculateLuminosityAtPhaseEnd(m_Mass0); }                      // Use class member variables



    double          CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const;
    double          CalculateLuminosityOnPhase() const{ return CalculateLuminosityOnPhase(m_Age, m_Mass0, CalculateLuminosityAtZAMS(m_Mass0)); } // Use class member variables







    double          CalculateLuminosityShikauchi(const double p_CoreMass, const double p_HeliumAbundanceCore) const;
    double          CalculateLuminosityTransitionToHG(const double p_Mass, const double p_Age, double const p_LZAMS) const;
    DBL_DBL         CalculateMainSequenceCoreMassBrcek(const double p_Dt, const double p_MassLossRate);
    double          CalculateInitialMainSequenceCoreMass(const double p_Mass, const double p_HeliumAbundanceCore) const;
    double          CalculateMomentOfInertia() const                                        { return (0.1 * (m_Mass) * m_Radius * m_Radius); }                      // k2 = 0.1 as defined in Hurley et al. 2000, after eq 109

    double          CalculatePerturbationMu() const                                         { return 5.0; }                                                         // mu(MS) = 5.0 (Hurley et al. 2000, eqs 97 & 98)

    double          CalculateRadialExtentConvectiveEnvelope() const;

    double          CalculateRadiusOnMassChange(double p_dM)                                { return CalculateRadiusOnPhase(m_Mass + p_dM, m_Tau, CalculateRadiusAtZAMS(m_Mass + p_dM)); }
    

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS) const;
    double          CalculateRadiusAtPhaseEnd() const                                       { return CalculateRadiusAtPhaseEnd(m_Mass, m_RZAMS); }                  // Use class member variables


//    double CalculateRadiusOnPhase_Hurley(const double p_Mass, const double p_Tau, const double p_RZAMS) const override;
    double CalculateRadiusOnPhase_Hurley() const override;



    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Tau, const double p_RZAMS) const;
    double          CalculateRadiusOnPhase() const override;


    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const { return Radius(); }                                               // Not a meaningful function for MS stars
    double          CalculateRadiusTransitionToHG(const double p_Mass, const double p_Tau, const double p_RZAMS) const;
     



    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables

    double          CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription)  { return OPTIONS->ZetaMainSequence(); }

    double          CalculateZetaEquilibrium();
    double          ChooseTimestep(const double p_Time) const;
    STELLAR_TYPE    EvolveToNextPhase()                                                     { return STELLAR_TYPE::HERTZSPRUNG_GAP; }




double CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency = 0.0) const;

                                           


    std::tuple <DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> InterpolateShikauchiCoefficients(const double p_Metallicity) const;
    
    bool            IsEndOfPhase() const                                                    { return !ShouldEvolveOnPhase(); }                                      // Phase ends when age at or after MS timescale

    void            PerturbLuminosityAndRadius() { }                                                                                                                // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);

    bool            ShouldEvolveOnPhase() const                                             { return (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]); }    // Evolve on MS phase if age in MS timescale
    
    double          CalculateTAMSCoreMass() const;


   
    void            UpdateAfterMerger(double p_Mass, double p_HydrogenMass);
    


double CalculateAgeAfterMassLoss() const override;
GNU_CONST double CalculateAgeAfterMassLoss_Hurley(const double p_Mass, const double p_Age, const DBL_VECTOR& p_tScales, const DBL_VECTOR& p_aN) const

GNU_CONST DBL_VECTOR CH::CalculateTimescales_Hurley(const double p_Mass, const DBL_VECTOR& p_tScales, const DBL_VECTOR& p_aN) const;


    void            UpdateMainSequenceCoreMass(const double p_Dt, const double p_MassLossRate);

    void            UpdateEffectiveZAMSLandR()                                              {
                                                                                                m_LZAMS0 = CalculateLuminosityAtZAMS(m_Mass0);
                                                                                                m_RZAMS0 = CalculateRadiusAtZAMS(m_Mass0);
                                                                                            }
};



///////////////////////// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                             LIFETIME / AGE FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculatePhaseLifetime
 *
 * @brief
 * Calculate lifetime of this phase of the evolution of the star (Main Sequence, tMS).
 *
 * Calls relevant lifetime function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculatePhaseLifetime() const
 *
 * @return                                      MS lifetime, tMS (Myr)
 */
inline double MainSequence::CalculatePhaseLifetime() const { 

    double lifetime;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            lifetime = CalculatePhaseLifetime_Hurley2000(m_StateHistory.CurrentState.Mass(), 
                                                         GLOBALS->ZetaHurley(), 
                                                         m_StateHistory.CurrentState.Timescales(tBGB), 
                                                         GLOBALS->HurleyACoefficients());
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
 * Calculate lifetime of the Main Sequence, tMS, per Hurley et al. 2000, eqs 5, 6, & 7.
 *
 *
 * double CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_ZetaHurley, const double p_TBGB, const DBL_VECTOR& p_aN) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_ZetaHurley                    Hurley zeta value (log10(Z / ZSOL_HURLEY))
 * @param       p_tBGB                          Lifetime to Base of Giant Branch, tBGB (Myr)
 * @param       p_aN                            Hurley a(n) coefficients
 * @return                                      MS lifetime, tMS (Myr)
 */
GNU_CONST inline double MainSequence::CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_ZetaHurley, const double p_TBGB, const DBL_VECTOR& p_aN) const {

    const double tHook = std::max(0.5, (1.0 - (0.01 * std::max((p_aN[6] / PPOW(p_Mass, p_aN[7])), (p_aN[8] + (p_aN[9] / PPOW(p_Mass, p_aN[10]))))))) * p_tBGB;
    const double x     = std::max(0.95, std::min((0.95 - (0.03 * (p_ZetaHurley + 0.30103))), 0.99));

    return std::max(tHook, (x * p_tBGB));
}


/*
 * CalculateTau_Hurley2000
 *
 * @brief
 * Calculate the MS-relative age (fractional Main Sequence age) of the star,
 * per Hurley et al. 2000, eq 11
 *
 * 
 * double CalculateTau_Hurley2000(const double p_Age, const double p_tMS) const
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tMS                           MS lifetime (Myr)
 * @return                                      MS-relative age, [0, 1]
 */
GNU_CONST inline double MainSequence::CalculateTau_Hurley2000(const double p_Age, const double p_tMS) const {
    return std::max(0.0, std::min(1.0, p_Age / p_tMS));
}



/*
 * CalculateHurleyEta
 *
 * @brief
 * Calculate the Hurley radius exponent eta, per Hurley et al. 2000, eq 18
 * 
 * Eta is used in Hurley et al. 2000, eq 12 (see CalculateLuminosityOnPhase_Hurley())
 *
 *
 * double CalculateHurleyEta(const double p_Metallicity, const double p_Mass) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Radius exponent, eta
 */
double MainSequence::CalculateHurleyEta(const double p_Metallicity, const double p_Mass) const {
    return p_Metallicity > 0.0009 ? 10.0 : (p_Mass > 1.0 ? (p_Mass >= 1.1 ? 20.0 : (100.0 * p_Mass) - 90.0) : 10.0);
}





/*
 * CalculateAgeAfterMassLoss_Hurley2000
 *
 * @brief
 * Recalculate the star's age after mass loss, per Hurley et al. 2000, section 7.1
 *
 *
 * double CalculateAgeAfterMassLoss_Hurley2000(const double      p_Mass,
 *                                             const double      p_Age,
 *                                             const double      p_ZetaHurley,
 *                                             const double      p_tMS,
 *                                             const DBL_VECTOR& p_aN) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Current age of the star (Myr)
 * @param       p_ZetaHurley                    Hurley zeta value (log10(Z / ZSOL_HURLEY))
 * @param       p_tMS                           MS lifetime (Myr)
 * @param       p_aN                            Hurley a(n) coefficients
 * @return                                      Age of the star after mass loss (Myr)
 */
inline double MainSequence::CalculateAgeAfterMassLoss_Hurley2000(const double      p_Mass,
                                                             const double      p_Age,
                                                             const double      p_ZetaHurley,
                                                             const double      p_tMS,
                                                             const DBL_VECTOR& p_aN) const {
    return p_Age * CalculateMSLifetime_Hurley2000(p_Mass, p_ZetaHurley, CalculateLifetimeToBGB_Hurley2000(p_Mass, p_aN), p_aN) / p_tMS;
}


/*
 * CalculateAgeAfterMassLoss
 *
 * @brief
 * Recalculate the star's age after mass loss.
 * 
 * Calls relevant age function based on the evolutionary mode given in program options.
 *
 *
 * double CalculateAgeAfterMassLoss() const
 *
 * @return                                      Age of the star after mass loss (Myr)
 */
inline double MainSequence::CalculateAgeAfterMassLoss() const {
    double age;

    Switch (OPTIONS->Mode()) {                                                                                  // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                                             // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                                             // HURLEY BSE
            age = CalculateAgeAfterMassLoss_Hurley2000(m_StateHistory.CurrentState.Mass(),
                                                       m_StateHistory.CurrentState.Age(),
                                                       GLOBALS->ZetaHurley(),
                                                       m_StateHistory.CurrentState.Timescales(tMS),
                                                       GLOBALS->HurleyACoefficients());
            break;

        default:                                                                                                // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                                         // throw error
    }

    return age;
}


/*
 * CalculateLuminosity
 *
 * @brief
 * Calculate the luminosity on the main sequence.
 * Uses the method appropriate for the evolution method and as specified by program options.
 * 
 * 
 * double CalculateLuminosity(const EVOLUTION_MODE p_Mode,
 *                            const double         p_Metallicity,
 *                            const double         p_Mass,
 *                            const double         p_Time,
 *                            const double         p_LZAMS,
 *                            const DBL_VECTOR&    p_tScales,
 *                            const DBL_VECTOR&    p_aN,
 *                            const DBL_VECTOR&    p_LConsts) const
 *
 * @param       p_Mode                          Evolution mode
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Time                          Time elapsed since ZAMS (Myr)
 * @param       p_LZAMS                         ZAMS luminosity of the star (Lsol)
 * @param       p_tScales                       Hurley timescales
 * @param       p_aN                            Hurley a(n) coefficients
 * @param       p_LConsts                       Hurley luminosity constants
 * @return                                      MS luminosity (Lsol))
 */
inline double MainSequence::CalculateLuminosity(const EVOLUTION_MODE p_Mode,
                                                const double         p_Metallicity,
                                                const double         p_Mass,
                                                const double         p_Time,
                                                const double         p_LZAMS,
                                                const DBL_VECTOR&    p_Timescales,                                         
                                                const DBL_VECTOR&    p_aN,
                                                const DBL_VECTOR&    p_LConsts) const {
    double luminosity;

    Switch (p_Mode) {                                                                                           // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                                             // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                                             // HURLEY BSE
            luminosity = CalculateLuminosityOnPhase_Hurley(p_Metallicity, p_Mass, p_Time, p_LZAMS, p_tScales, p_aN, p_LConsts);
            break;

        default:                                                                                                // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                                         // throw error
    }

    return luminosity;
}




/*
 * CalculateRadius
 *
 * @brief
 * Calculate the radius on the main sequence.
 * Uses the method appropriate for the evolution method and as specified by program options.
 * 
 * 
 * double CalculateRadius(const EVOLUTION_MODE p_Mode,
 *                        const double         p_Mass,
 *                        const double         p_Tau,
 *                        const double         p_RZAMS,
 *                        const DBL_VECTOR&    p_aN) const
 *
 * @param       p_Mode                          Evolution mode
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Fractional MS age of the star
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_aN                            Hurley a(n) coefficients
 * @return                                      MS radius (Rsol)
 */
inline double MainSequence::CalculateRadius(const EVOLUTION_MODE p_Mode,
                                            const double         p_Mass,
                                            const double         p_Tau,
                                            const double         p_RZAMS,
                                            const DBL_VECTOR&    p_aN) const {
    double radius;

    Switch (p_Mode) {                                                                                           // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                                             // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                                             // HURLEY BSE
            radius = CalculateRadiusOnPhase_Hurley(p_Mass, p_Tau, p_RZAMS, p_aN);
            break;

        default:                                                                                                // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                                         // throw error
    }

    return radius;
}


/*
 * CalculateRadiusAtPhaseEnd
 *
 * @brief
 * Calculate radius at the end of the Main Sequence (TAMS).
 * Uses the method appropriate for the evolution method and as specified by program options.
 * 
 * 
 * double CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS, const DBL_VECTOR& p_aN) const
 *
 * @param       p_Mode                          Evolution mode
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_aN                            Hurley a(n) coefficients
 * @return                                      TAMS radius (Rsol)
 */
double MainSequence::CalculateRadiusAtPhaseEnd(const EVOLUTION_MODE p_Mode, const double p_Mass, const double p_RZAMS, const DBL_VECTOR& p_aN) const {

    double radius;

    Switch (p_Mode) {                                                                                           // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                                             // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                                             // HURLEY BSE
            radius = CalculateRadiusAtPhaseEnd_Hurley(p_Mass, p_RZAMS, p_aN);
            break;

        default:                                                                                                // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                                         // throw error
    }

    return radius;
}



#endif // __MainSequence_h__
