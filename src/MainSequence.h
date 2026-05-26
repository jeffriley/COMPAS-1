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
    


    


    ////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES   //
    ////////////////////////////////////////

    COMPAS_PURE inline double CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_tBGB) const;

    inline double CalculateTau_Hurley2000() const override { return CalculateTau_Hurley2000(Age(), m_InterimState.Timescale(TIMESCALE::tMS)); }
    GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const double p_tMS) const;

    inline DblVectorT CalculateTimescales_Hurley2000() const override { return CalculateTimescales_Hurley2000(Mass0(), m_InterimState.TimeScales()); }
    COMPAS_PURE inline DblVectorT CalculateTimescales_Hurley2000(const double p_Mass, const DblVectorT& p_tScales) const;


    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////

    // LZAMS() will exist for MS  <<<<<<<<<<<<<<<<< CHECK IF CALLED FROM OTHER STELLAR TYPES!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<
    inline double CalculateLuminosity_Hurley2000() const override { return CalculateLuminosity_Hurley2000(Mass(), Age(), LZAMS(), m_InterimState.Timescale(TIMESCALE::tMS), m_InterimState.Timescale(TIMESCALE::tBGB)); }
    COMPAS_PURE double CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Age, const double p_LZAMS, const double p_tMS, const double p_tBGB) const;

GNU_CONST double CalculateLuminosityAtPhaseEnd_Hurley(const double p_Mass, const DBL_VECTOR& p_aN) const;

GNU_CONST double CalculateLuminosity_Brcek(const double      p_Mass, 
                                               const double      p_Tau,
                                               const double      p_MZAMS,
                                               const double      p_LZAMS,
                                               const double      p_CoreMass,
                                               const DBL_VECTOR& p_Timescales,
                                               const DBL_VECTOR& p_aN,
                                               const DBL_VECTOR& p_ShikauchiLCoeffs) const;


GNU_CONST double CalculateLuminosity_Shikauchi(const double p_CoreMass, const double p_HeAbundanceCore, const DBL_VECTOR& p_ShikauchiLCoeffs) const;

        
    ////////////////////////////////////////
    //   MASS                             //
    ////////////////////////////////////////

    // MZAMS() will exist for MS if we started evolving on the MS, otherwise will fail

    // Here we use the mass loss rate from the interim state
    inline double CalculateCoreMass() const override { return CalculateCoreMass(Mass(), Luminosity(), MZAMS(), CoreMass(), HeAbundanceCore(), dt(), m_InterinState.dMdt()); }


GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; } // McCO = 0.0 for MS stars

GNU_CONST inline double CalculateHeCoreMass() const { return 0.0; } // McHe = 0.0 for MS stars

inline Dbl_DblT CalculateConvectiveEnvelopeMass() const override { return CalculateConvectiveEnvelopeMass(Mass(), Tau()); }
GNU_CONST inline Dbl_DblT CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_Tau) const override { return CalculateConvectiveEnvelopeMass_Hurley2000(p_Mass, p_Tau); }
GNU_CONST Dbl_DblT CalculateConvectiveEnvelopeMass_Hurley2000(const double p_Mass, const double p_Tau) const;

inline double CalculateEffectiveInitialMass_Hurley2000() const override { return Mass(); } // per Hurley et al. 2000, section 7.1


Dbl_DblT         CalculateConvectiveEnvelopeMass() const;

double          CalculateConvectiveCoreMass() const;


    
    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////

    // RZAMS() will exist for MS if we started evolving on the MS, otherwise will fail

                inline double CalculateRadius_Hurley2000() const override { return CalculateRadius_Hurley2000(Mass(), Tau(), RZAMS(), Timescale(TIMESCALE::tBGB)); }
    COMPAS_PURE        double CalculateRadius_Hurley2000(const double p_Mass, const double p_Tau, const double p_RZAMS, const double p_tBGB) const;


                inline double CalculateRadiusAtPhaseEnd_Hurley2000() const override { return CalculateRadiusAtPhaseEnd_Hurley2000_Static(Mass(), RZAMS()); }
    COMPAS_PURE static double CalculateRadiusAtPhaseEnd_Hurley2000_Static(const double p_Mass, const double p_RZAMS) const;



double          CalculateConvectiveCoreRadius() const;







    










 



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; }  // McCO = 0.0 for MS stars

GNU_CONST inline double CalculateHeCoreMassAtPhaseEnd() const { return 0.0; } // McHe = 0.0 for MS stars





    double          CalculateCoreMassAtPhaseEnd() const                                     { return (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::MANDEL) ? std::min(MainSequenceCoreMass(), m_Mass) : 0.0; }        // Accounts for minimal core mass built up prior to mass loss through mass transfer; core mass can't exceed total mass
    double          CalculateCoreMassOnPhase() const                                        { return 0.0; }                                                         // Mc(MS) = 0.0 (Hurley et al. 2000, just before eq 28)

    double          CalculateHeliumAbundanceCoreAtPhaseEnd() const                          { return CalculateHeliumAbundanceCoreOnPhase(); }


COMPAS_PURE double CalculateHAbundanceCoreOnPhase(const double p_Tau) const override;                                
inline double CalculateHAbundanceSurfaceOnPhase(const double p_Tau) const override { return HAbundanceSurface(); }

COMPAS_PURE double CalculateHeAbundanceCoreOnPhase(const double p_Tau) const override;
inline double CalculateHeAbundanceSurfaceOnPhase(const double p_Tau) const override { return HeAbundanceSurface(); }

    


    double          CalculateHeliumAbundanceSurfaceAtPhaseEnd() const                       { return CalculateHeliumAbundanceSurfaceOnPhase(); }
    
    double          CalculateHydrogenAbundanceCoreAtPhaseEnd() const                        { return CalculateHydrogenAbundanceCoreOnPhase(); } 
    double          CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const;                                                          
    
    double          CalculateHydrogenAbundanceSurfaceAtPhaseEnd() const                     { return CalculateHydrogenAbundanceSurfaceOnPhase(); } 

    






    double          CalculateLuminosityAtPhaseEnd() const                                   { return CalculateLuminosityAtPhaseEnd(m_Mass0); }                      // Use class member variables



    double          CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const;
    double          CalculateLuminosityOnPhase() const{ return CalculateLuminosityOnPhase(m_Age, m_Mass0, CalculateLuminosityAtZAMS(m_Mass0)); } // Use class member variables







    double          CalculateLuminosityShikauchi(const double p_CoreMass, const double p_HeliumAbundanceCore) const;
    double          CalculateLuminosityTransitionToHG(const double p_Mass, const double p_Age, double const p_LZAMS) const;
    Dbl_DblT         CalculateMainSequenceCoreMassBrcek(const double p_Dt, const double p_MassLossRate);


COMPAS_PURE double CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024(const double p_MZAMS) const;
COMPAS_PURE double CalculateCNOprocessedCoreMass_Brcek2025(const double p_Mass, const double p_HeAbundanceCore) const;


    double          CalculateMomentOfInertia() const                                        { return (0.1 * (m_Mass) * m_Radius * m_Radius); }                      // k2 = 0.1 as defined in Hurley et al. 2000, after eq 109


GNU_CONST inline double CalculateHurleyPerturbationMu() const { return 5.0; } // Hurley et al. 2000, eqs 97 & 98


    double CalculateConvectiveEnvelopeRadialExtent() const override { return CalculateConvectiveEnvelopeRadialExtent_Hurley2002(Mass(), Radius()); }
    COMPAS_PURE double CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Mass, const double p_Radius) const override;

    double          CalculateRadiusOnMassChange(double p_dM)                                { return CalculateRadius(m_Mass + p_dM, m_Tau, CalculateRadiusAtZAMS(m_Mass + p_dM)); }
    
    double          CalculateRadiusTransitionToHG(const double p_Mass, const double p_Tau, const double p_RZAMS) const;
     





    double          CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription)  { return OPTIONS->ZetaMainSequence(); }

    double          CalculateZetaEquilibrium();


    GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const override;


    STELLAR_TYPE    EvolveToNextPhase()                                                     { return STELLAR_TYPE::HERTZSPRUNG_GAP; }




double CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency = 0.0) const;

                                           


    std::tuple <DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> InterpolateShikauchiCoefficients() const;
    
    bool            IsEndOfPhase() const                                                    { return !ShouldEvolveOnPhase(); }                                      // Phase ends when age at or after MS timescale

    void            PerturbLuminosityAndRadius() { }                                                                                                                // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);

    bool            ShouldEvolveOnPhase() const                                             { return (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]); }    // Evolve on MS phase if age in MS timescale
    

GNU_CONST double CalculateCoreMassAtTAMS_Hurely2000() const;


   
    void            UpdateAfterMerger(double p_Mass, double p_HydrogenMass);
    


double CalculateAgeAfterMassLoss() const override;
GNU_CONST double CalculateAgeAfterMassLoss_Hurley(const double p_Mass, const double p_Age, const DBL_VECTOR& p_tScales, const DBL_VECTOR& p_aN) const


    void            UpdateMainSequenceCoreMass(const double p_Dt, const double p_MassLossRate);

    void            UpdateEffectiveZAMSLandR()                                              {
                                                                                                m_LZAMS0 = CalculateLuminosityAtZAMS(m_Mass0);
                                                                                                m_RZAMS0 = CalculateRadiusAtZAMS(m_Mass0);
                                                                                            }
};



///////////////////////// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculatePhaseLifetime_Hurley2000
 *
 * @brief
 * Calculate lifetime of the Main Sequence, tMS, per Hurley et al. 2000, eqs 5, 6, & 7.
 *
 *
 * double CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_TBGB) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (per Hurley timescales) (Myr)
 * @return                                      MS lifetime, tMS (Myr)
 */
inline double MainSequence::CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_tBGB) const {

    const DBL_VECTOR a = GLOBALS->HurleyAcoefficients(); // get Hurley a coefficients
    const double tHook = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10]))))))) * p_tBGB;

    return std::max(tHook, (std::max(0.95, std::min((0.95 - (0.03 * (GLOBALS->ZetaHurley() + 0.30103))), 0.99)) * p_tBGB));
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
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tMS                           MS lifetime, tMS (per Hurley timescales) (Myr)
 * @return                                      MS-relative age, [0, 1]
 */
inline double MainSequence::CalculateTau_Hurley2000(const double p_Age, const double p_tMS) const {
    return std::max(0.0, std::min(1.0, p_Age / p_tMS));
}


/*
 * CalculateTimescales_Hurley2000
 *
 * @brief
 * (Re)calculate timescales given the mass of the star, per Hurley at al. 2000.
 * 
 * Since timescales depend on a star's mass, they need to be calculated whenever
 * the mass of the star changes (probably every timestep).
 *
 *
 * DblVectorT CalculateTimescales_Hurley2000(const double p_Mass, const DblVectorT& p_tScales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
inline DblVectorT MainSequence::CalculateTimescales_Hurley2000(const double p_Mass, const DblVectorT& p_tScales) const {

    DblVectorT tScales = p_tScales; // copy given timescales

    // (re)calculate tBGB and tMS
    tScales[HURLEY_TS::BGB] = utils::CalculateLifetimeToBGB_Hurley2000(p_Mass);
    tScales[HURLEY_TS::MS]  = CalculatePhaseLifetime_Hurley2000(p_Mass, tScales[HURLEY_TS::BGB]);

    // return timescales vector by value - NRVO takes care of performance/efficiency
    return tScales;
}


/*
 * CalculateAgeAfterMassLoss_Hurley2000
 *
 * @brief
 * Recalculate the star's age after mass loss, per Hurley et al. 2000, section 7.1
 *
 *
 * double CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tMS) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tMS                           MS lifetime, tMS (per Hurley timescales) (Myr)
 * @return                                      Age of the star after mass loss (Myr)
 */
inline double MainSequence::CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tMS) const {
    // We call utils::CalculateLifetimeToBGB_Hurley2000() here rather than use timescales
    // because p_Mass may not be the same mass used to calculate timescales[HURLEY_TS::BGB]
    return p_Age * CalculatePhaseLifetime_Hurley2000(p_Mass, utils::CalculateLifetimeToBGB_Hurley2000(p_Mass)) / p_tMS;
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

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            age = CalculateAgeAfterMassLoss_Hurley2000();
            break;

        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }

    return age;
}


/*
 * ChooseTimestep_Hurley2000
 *
 * @brief
 * Choose timestep for evolution
 * See the discussion in Hurley et al. 2000, p21
 * The returned value will be clamped to minimum NUCLEAR_MINIMUM_TIMESTEP
 *
 *
 * double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales)
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Phase timescales (Myr)
 * @return                                      Suggested timestep (Myr)
 */
GNU_CONST inline double MainSequence::ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const {
#define tScales(x) p_tScales[static_cast<int>(TIMESCALE::x)]

          double dtk = 1.0E-2 * tScales(tMS);                       // stellar type specific dt (Hurley sse uses 0.05)
    const double dte = tScales(tMS) - p_Age;                        // time to end of phase (change of stellar type)

    // check that dtk is short enough to resolve the hook at the
    // end of the MS for HM stars? (JAR: why not check for HM star?)
    // if not, go an order-of-magnitude shorter
    if (dte < dtk) dtk /= 10.0;

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);  // clamp to minimum NUCLEAR_MINIMUM_TIMESTEP

#undef tScales
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLuminosityAtZAMS
 *
 * @brief
 * Calculate the luminosity of a star at ZAMS.
 * 
 * Calls relevant luminosity function based on the evolutionary mode given in program options.
 *
 *
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function ostensibly relies on the values of the ZAMS mass of the star,
 * and should not be used if the ZAMS mass is not known.
 * 
 * 
 * double CalculateLuminosityAtZAMS(const double p_MZAMS, const DBL_VECTOR& p_LCoeffs) const
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_LCoeffs                       Tout luminosity coefficients
 * @return                                      ZAMS luminosity (Lsol)
 */
 double MainSequence::CalculateLuminosityAtZAMS(const double p_MZAMS, const DBL_VECTOR& p_LCoeffs) const {  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< lcoeffs????????????????????

    double luminosity;

    Switch (OPTIONS->Mode()) {                                              // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                         // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                         // HURLEY BSE
            luminosity = CalculateLuminosityAtZAMS_Tout1996(p_MZAMS, p_LCoeffs);
            break;
        
        default:                                                            // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                     // throw error
    }       

    return luminosity;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024
 *
 * @brief
 * Calculate the CNO-processed core mass of a main sequence star at ZAMS,
 * per Shikauchi et al. 2024
 *
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function ostensibly relies on the value of the ZAMS mass of the star,
 * and should not be used if the ZAMS mass is not known.
 * 
 * 
 * double CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024(const double p_MZAMS) const
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS CNO-processed core mass (Msol)
 */
COMPAS_PURE inline double MainSequence::CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024(const double p_MZAMS) const {
    DblVectorT fCoeffs = GLOBALS->ShikauskifMixCoefficients(); // Get Shikauchi fMix coefficients from GLOBALS
    return p_MZAMS * (fCoeffs[0] + fCoeffs[1] * std::exp(-p_MZAMS / fCoeffs[2]));
}


/*
 * CalculateCNOprocessedCoreMass_Brcek2025
 *
 * @brief
 * Calculate the CNO-processed core mass of a main sequence star after full mixing (due to
 * merger or CHE) for an arbitrary central helium fraction, per Brcek et al. 2025
 * 
 *
 * double CalculateCNOprocessedCoreMass_Brcek2025(const double p_Mass, const double p_HeAbundanceCore)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core of the star 
 * @return                                      CNO-processed core mass (Msol)
 */
COMPAS_PURE inline double MainSequence::CalculateCNOprocessedCoreMass_Brcek2025(const double p_Mass, const double p_HeAbundanceCore) const {
    DbVectorT fCoeffs = BRCEK_FMIX_COEFFICIENTS; // Get Brcek fMix coefficients from constants.h
    const double mh   = p_Mass * PPOW(10.0, p_HeAbundanceCore * (p_HeAbundanceCore + 2.0) / 4.0);
    return p_Mass * (fCoeffs[0] + fCoeffs[1] * std::exp(-mh / fCoeffs[2])) * PPOW(1.0 - fCoeffs[4] / mh, fCoeffs[3]);
}


/*
 * CalculateCoreMassAtTAMS_Hurely2000
 *
 * @brief
 * Calculate the (expected) core mass at terminal age main sequence, per Hurley et al. 2000.
 * We do this by calculating the core mass at the start of the HG phase.
 *
 * double calculate CalculateCoreMassAtTAMS_Hurely2000() const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAMS core mass (Msol)
 *
 */
GNU_CONST inline double MainSequence::CalculateCoreMassAtTAMS_Hurely2000(const double p_Mass) const {
    
    // We need TAMS core mass, which is just the core mass at the start of the HG phase.
    // Since we are on the main sequence here, we can clone this object as an HG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HG phase),
    // we can  query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.

    std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::HG, OBJECT_PERSISTENCE::EPHEMERAL);

    return clone->CalculateCoreMass();
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


#endif // __MainSequence_h__
