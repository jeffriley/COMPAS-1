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
    GNU_CONST double CalculateZAMSLuminosity_Tout_Static(const double p_MZAMZ, const DBL_VECTOR& p_LCoefficients) const;
    GNU_CONST double CalculateZAMSRadius_Tout_Static(const double p_MZAMS, const DBL_VECTOR& p_RCoefficients) const;


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
    
    // luminosity
    GNU_CONST double CalculateLuminosity(const double     p_Metallicity,
                                         const double      p_Mass,
                                         const double      p_Time,
                                         const double      p_LZAMS,
                                         const DBL_VECTOR& p_Timescales,
                                         const DBL_VECTOR& p_aN,
                                         const DBL_VECTOR& p_LConsts) const;

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
    GNU_CONST double CalculateRadius(const double p_Mass, const double p_Tau, const double p_RZAMS, const DBL_VECTOR& p_aN) const;

    GNU_CONST double CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS, const DBL_VECTOR& p_aN) const;

    GNU_PURE  double CalculateRadiusAtPhaseEnd_Hurley(const double p_Mass, const double p_RZAMS, const DBL_VECTOR& p_aN) const {
        return CalculateRadiusAtPhaseEnd_Hurley_Static(p_Mass, p_RZAMS, p_aN);
    }

    GNU_CONST double CalculateRadiusOnPhase_Hurley(const double p_Mass, const double p_Tau, const double p_RZAMS, const DBL_VECTOR& p_aN) const;

    double CalculateRadiusOnPhase_Hurley() const {
        // use globals and state values - read only
        return CalculateRadiusOnPhase_Hurley(m_StateHistory.CurrentState().Mass(), m_StateHistory.CurrentState().Tau(), m_StateHistory.ZAMSState().Radius(), GLOBALS->HurleyACoefficients());
    }




    double          CalculateConvectiveCoreMass() const;
    double          CalculateConvectiveCoreRadius() const;
    DBL_DBL         CalculateConvectiveEnvelopeMass() const;
    double          CalculateBetaL(const double p_Mass) const;
    double          CalculateBetaR(const double p_Mass) const;

    double          CalculateDeltaR(const double p_Mass) const;


    double          CalculateGamma(const double p_Mass) const;

    double          CalculateCOCoreMassAtPhaseEnd() const                                   { return CalculateCOCoreMassOnPhase(); }                                // Same as on phase
    double          CalculateCOCoreMassOnPhase() const                                      { return 0.0; }                                                         // McCO(MS) = 0.0

    double          CalculateCoreMassAtPhaseEnd() const                                     { return (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::MANDEL) ? std::min(MainSequenceCoreMass(), m_Mass) : 0.0; }        // Accounts for minimal core mass built up prior to mass loss through mass transfer; core mass can't exceed total mass
    double          CalculateCoreMassOnPhase() const                                        { return 0.0; }                                                         // Mc(MS) = 0.0 (Hurley et al. 2000, just before eq 28)

    double          CalculateHeCoreMassAtPhaseEnd() const                                   { return CalculateCoreMassAtPhaseEnd(); }                               // Same as He core mass
    double          CalculateHeCoreMassOnPhase() const                                      { return 0.0; }                                                         // McHe(MS) = 0.0

    double          CalculateHeliumAbundanceCoreAtPhaseEnd() const                          { return CalculateHeliumAbundanceCoreOnPhase(); }
    double          CalculateHeliumAbundanceCoreOnPhase(const double p_Tau) const;                                         
    double          CalculateHeliumAbundanceCoreOnPhase() const                             { return CalculateHeliumAbundanceCoreOnPhase(m_Tau); }                  // Use class member variables                                       
    
    double          CalculateHeliumAbundanceSurfaceAtPhaseEnd() const                       { return CalculateHeliumAbundanceSurfaceOnPhase(); }
    double          CalculateHeliumAbundanceSurfaceOnPhase() const                          { return m_HeliumAbundanceSurface; }                                    // Use class member variables
    
    double          CalculateHydrogenAbundanceCoreAtPhaseEnd() const                        { return CalculateHydrogenAbundanceCoreOnPhase(); } 
    double          CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const;                                                          
    double          CalculateHydrogenAbundanceCoreOnPhase() const                           { return CalculateHydrogenAbundanceCoreOnPhase(m_Tau); }                // Use class member variables                                 
    
    double          CalculateHydrogenAbundanceSurfaceAtPhaseEnd() const                     { return CalculateHydrogenAbundanceSurfaceOnPhase(); } 
    double          CalculateHydrogenAbundanceSurfaceOnPhase() const                        { return m_HydrogenAbundanceSurface; }                                  // Use class member variables
    
    double          CalculateLifetimeOnPhase(const double p_Mass, const double p_TBGB) const;

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
     
    double          CalculateTauAtPhaseEnd() const                                          { return 1.0; }                                                         // tau = 1.0 at end of MS
    double          CalculateTauOnPhase() const;

    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                                   { CalculateTimescales(m_Mass0, m_Timescales); }                         // Use class member variables

    double          CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription)  { return OPTIONS->ZetaMainSequence(); }

    double          CalculateZetaEquilibrium();
    double          ChooseTimestep(const double p_Time) const;
    STELLAR_TYPE    EvolveToNextPhase()                                                     { return STELLAR_TYPE::HERTZSPRUNG_GAP; }

    double          InterpolateGeEtAlQCrit(const QCRIT_PRESCRIPTION p_qCritPrescription, 
                                           const double p_massTransferEfficiencyBeta); // RTW do I need a const here?
    
    std::tuple <DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> InterpolateShikauchiCoefficients(const double p_Metallicity) const;
    
    bool            IsEndOfPhase() const                                                    { return !ShouldEvolveOnPhase(); }                                      // Phase ends when age at or after MS timescale

    void            PerturbLuminosityAndRadius() { }                                                                                                                // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);

    bool            ShouldEvolveOnPhase() const                                             { return (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]); }    // Evolve on MS phase if age in MS timescale
    
    double          CalculateTAMSCoreMass() const;

    double          CalculateEffectiveInitialMass()                                                     { return m_Mass; }                                  // Per Hurley et al. 2000, section 7.1
   
    void            UpdateAfterMerger(double p_Mass, double p_HydrogenMass);
    
    void            UpdateAgeAfterMassLoss();                                                                                                                       // Per Hurley et al. 2000, section 7.1

    void            UpdateMainSequenceCoreMass(const double p_Dt, const double p_MassLossRate);

    void            UpdateEffectiveZAMSLandR()                                              {
                                                                                                m_LZAMS0 = CalculateLuminosityAtZAMS(m_Mass0);
                                                                                                m_RZAMS0 = CalculateRadiusAtZAMS(m_Mass0);
                                                                                            }
};

#endif // __MainSequence_h__
