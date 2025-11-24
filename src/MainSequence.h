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




    // member functions - alphabetically

        
    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////


    
    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////

    inline double CalculateRadius_Hurley2000() const override { return CalculateRadius_Hurley2000(Mass(), Tau(), RZAMS(), Timescale(TIMESCALE::tBGB)); } // RZAMS() will exist for MS - if called with MS having existed, will fail

    COMPAS_PURE double CalculateRadius_Hurley2000(const double p_Mass, const double p_Tau, const double p_RZAMS, const double p_tBGB) const;

    inline double CalculateRadiusAtPhaseEnd_Hurley2000() const override { return CalculateRadiusAtPhaseEnd_Hurley2000_Static(Mass(), RZAMS()); } // RZAMS()will exist for MS - if called with MS having existed, will fail

    COMPAS_PURE static double CalculateRadiusAtPhaseEnd_Hurley2000_Static(const double p_Mass, const double p_RZAMS);









///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; } // McCO = 0.0 for MS stars

GNU_CONST inline double CalculateHeCoreMass() const { return 0.0; } // McHe = 0.0 for MS stars

inline double CalculateLuminosity_Hurley2000() const override { return CalculateLuminosity_Hurley2000(Mass(), Age(), LZAMS(), Timescale(TIMESCALE::tMS), Timescale(TIMESCALE::tBGB)); } // LZAMS() will exist for MS  <<<<<<<<<<<<<<<<< CHECK IF CALLED FROM OTHER STELLAR TYPES!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<

COMPAS_PURE double CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Age, const double p_LZAMS, const double p_tMS, const double p_tBGB) const;





inline DBL_VECTOR CalculateTimescales_Hurley2000() const override { return CalculateTimescales_Hurley2000(Mass0(), TimeScales()); }

GNU_CONST inline DBL_VECTOR HeMS::CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_tScales) const;





inline DBL_DBL CalculateConvectiveEnvelopeMass() const override { return CalculateConvectiveEnvelopeMass_Hurley2000(Mass(), Tau()); }

GNU_CONST DBL_DBL CalculateConvectiveEnvelopeMass_Hurley2000(const double p_Mass, const double p_Tau) const;




    // luminosity



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



    // radius



inline double CalculateEffectiveInitialMass_Hurley2000() const override { return Mass(); } // per Hurley et al. 2000, section 7.1




    double          CalculateConvectiveCoreMass() const;
    
double          CalculateConvectiveCoreRadius() const;


    DBL_DBL         CalculateConvectiveEnvelopeMass() const;
    double          CalculateBetaL(const double p_Mass) const;
    double          CalculateBetaR(const double p_Mass) const;

    double          CalculateDeltaR(const double p_Mass) const;


    double          CalculateGamma(const double p_Mass) const;






 



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

    


COMPAS_PURE inline double CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_TBGB) const;


inline double CalculateTau_Hurley2000() const override { return CalculateTau_Hurley2000(Age(), Timescale(TIMESCALE::tMS)); }
GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const double p_tMS) const;

    double          CalculateLuminosityAtPhaseEnd() const                                   { return CalculateLuminosityAtPhaseEnd(m_Mass0); }                      // Use class member variables



    double          CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const;
    double          CalculateLuminosityOnPhase() const{ return CalculateLuminosityOnPhase(m_Age, m_Mass0, CalculateLuminosityAtZAMS(m_Mass0)); } // Use class member variables







    double          CalculateLuminosityShikauchi(const double p_CoreMass, const double p_HeliumAbundanceCore) const;
    double          CalculateLuminosityTransitionToHG(const double p_Mass, const double p_Age, double const p_LZAMS) const;
    DBL_DBL         CalculateMainSequenceCoreMassBrcek(const double p_Dt, const double p_MassLossRate);


COMPAS_PURE double CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024(const double p_MZAMS) const;
COMPAS_PURE double CalculateCNOprocessedCoreMass_Brcek2025(const double p_Mass, const double p_HeAbundanceCore) const;


    double          CalculateMomentOfInertia() const                                        { return (0.1 * (m_Mass) * m_Radius * m_Radius); }                      // k2 = 0.1 as defined in Hurley et al. 2000, after eq 109


GNU_CONST inline double CalculateHurleyPerturbationMu() const { return 5.0; } // Hurley et al. 2000, eqs 97 & 98


    double          CalculateConvectiveEnvelopeRadialExtent() const;

    double          CalculateRadiusOnMassChange(double p_dM)                                { return CalculateRadiusOnPhase(m_Mass + p_dM, m_Tau, CalculateRadiusAtZAMS(m_Mass + p_dM)); }
    

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS) const;
    double          CalculateRadiusAtPhaseEnd() const                                       { return CalculateRadiusAtPhaseEnd(m_Mass, m_RZAMS); }                  // Use class member variables


//    double CalculateRadiusOnPhase_Hurley(const double p_Mass, const double p_Tau, const double p_RZAMS) const override;
    double CalculateRadiusOnPhase_Hurley() const override;



    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Tau, const double p_RZAMS) const;
    double          CalculateRadiusOnPhase() const override;


    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const { return Radius(); }                                               // Not a meaningful function for MS stars
    double          CalculateRadiusTransitionToHG(const double p_Mass, const double p_Tau, const double p_RZAMS) const;
     





    double          CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription)  { return OPTIONS->ZetaMainSequence(); }

    double          CalculateZetaEquilibrium();
    double          ChooseTimestep(const double p_Time) const;
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
//                         AGE / LIFETIME / TAU / TIMESCALES                         //
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
COMPAS_PURE inline double MainSequence::CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_TBGB) const {

    const DBL_VECTOR a = GLOBALS->HurleyAcoefficients(); // get Hurley a coefficients
    const double tHook = std::max(0.5, (1.0 - (0.01 * std::max((p_aN[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10]))))))) * p_tBGB;

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
GNU_CONST inline double MainSequence::CalculateTau_Hurley2000(const double p_Age, const double p_tMS) const {
    return std::max(0.0, std::min(1.0, p_Age / p_tMS));
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
GNU_CONST inline double MainSequence::CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tMS) const {
                    // FIX THIS AFTER TALKING TO ILYA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return p_Age * CalculatePhaseLifetime_Hurley2000(p_Mass, BaseStar::CalculateLifetimeToBGB_Hurley2000_Static(p_Mass)) / p_tMS;
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
 * This function relies on the values of the ZAMS mass and ZAMS radius of the star,
 * and should not be used if the ZAMS mass or the ZAMS radius is not known.
 * 
 * 
 * double CalculateLuminosityAtZAMS(const double p_MZAMS, const DBL_VECTOR& p_LCoeffs) const
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_LCoeffs                       Tout luminosity coefficients
 * @return                                      ZAMS luminosity (Lsol)
 */
 double MainSequence::CalculateLuminosityAtZAMS(const double p_MZAMS, const DBL_VECTOR& p_LCoeffs) const {

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
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



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
 * This function relies on the value of the ZAMS mass of the star, and should not be used
 * if the ZAMS mass is not known.
 * 
 * 
 * double CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024(const double p_MZAMS) const
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS CNO-processed core mass (Msol)
 */
COMPAS_PURE inline double MainSequence::CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024(const double p_MZAMS) const {
    DBL_VECTOR fCoeffs = GLOBALS->ShikauskifMixCoefficients(); // get Shikauchi fMix coefficients from GLOBALS
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
    DBL_VECTOR fCoeffs = BRCEK_FMIX_COEFFICIENTS; // get Brcek fMix coefficients from constants.h
    const double mh    = p_Mass * PPOW(10.0, p_HeAbundanceCore * (p_HeAbundanceCore + 2.0) / 4.0);
    return p_Mass * (fCoeffs[0] + fCoeffs[1] * std::exp(-mh / fCoeffs[2])) * PPOW(1.0 - fCoeffs[4] / mh, fCoeffs[3]);
}


/*
 * CalculateCoreMassAtTAMS_Hurely2000
 *
 * @brief
 * Calculate the (expected) core mass at terminal age main sequence, per Hurley et al. 2000.
 * We do this by calculating the core mass at the start of the HG phase.
 *
 * static double calculate CalculateCoreMassAtTAMS_Hurely2000() const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAMS core mass (Msol)
 *
 */
GNU_CONST inline double MainSequence::CalculateCoreMassAtTAMS_Hurely2000(const double p_Mass) const {
    return HG::CalculateCoreMass_Hurley2000_Unconstrained_Static(const double p_Mass, 0.0);
}


#endif // __MainSequence_h__
