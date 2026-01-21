#ifndef __HeMS_h__
#define __HeMS_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "TPAGB.h"


// JR: todo: revisit this one day - sometimes HeMS is MS, sometimes it is GiantBranch...
// Right now it is GiantBranch - figure out which is more appropriate

class BaseStar;
class TPAGB;

class HeMS: virtual public BaseStar, public TPAGB {

public:

    HeMS() { m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS; };
    
    HeMS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), TPAGB(p_BaseStar, false) {
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                                                                                                                                         // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                                             // Initialise if required
    }


    // member functions - alphabetically














    static DBL_DBL  CalculateRadiusAtPhaseEnd_Static(const double p_Mass, const double p_Luminosity);




    



protected:

    void Initialise() {
        // initialise surface abundances
        m_HydrogenAbundanceSurface = 0.0;
        m_HeliumAbundanceSurface   = 1.0 - m_Metallicity;
        
        CalculateTimescales();
        // JR: Age for HeMS is partially calculated before switching -
        // can get here from various places in ResolveEnvelopeLoss(),
        // and Age is calculated differently in those cases
        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically
GNU_CONST inline DBL_DBL CalculateConvectiveEnvelopeMass() const override { return std::make_tuple(0.0, 0.0); } // No convective envelope for HeMS stars


            double          CalculateConvectiveCoreMass() const;
            double          CalculateConvectiveCoreRadius() const;
            double          CalculateCoreMassAtPhaseEnd() const                                                     { return CalculateHeCoreMassOnPhase(); } // Same as on phase /*ILYA*/ To fix, not everything will become CO core
            double          CalculateCoreMassOnPhase() const                                                        { return 0.0; }                                                                 // Mc(HeMS) = 0.0





///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


DBL_VECTOR CalculateTimescales_Hurley2000() const override {
    return CalculateTimescales_Hurley2000(m_StateHistory.CurrentState.MassEffectiveInitial(), m_StateHistory.CurrentState.TimeScales());
}
GNU_CONST inline DBL_VECTOR HeMS::CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_tScales) const;

DBL_VECTOR CalculateGBparams_Hurley2000() const override {
    return CalculateGBparams_Hurley2000(m_StateHistory.CurrentState.MassEffectiveInitial(), m_StateHistory.CurrentState.GBparams());
}
COMPAS_PURE virtual DBL_VECTOR CalculateGBparams_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const;



inline double CalculateRemnantRadius() const override { return m_StateHistory.CurrentState.Radius(); }




GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; } // McCO = 0.0 for HeMS stars

inline double CalculateHeCoreMass() const override { return Mass(); } // McHe(HeMS) = Mass

inline double CalculateEffectiveInitialMass_Hurley2000() const override { return Mass(); } // per Hurley et al. 2000, section 7.1





inline double CalculateLuminosity() const override;
GNU_CONST static double CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Tau);


// MAYBE BREAK THE DESCRIPTION OF THESE OUT TO INLINE BELOW
GNU_CONST inline double CalculateCoreMass_Luminosity_B_Hurley2000(const double p_Mass) const override { return 4.1E4; } // Hurley et al. 2000, just after eq 84  THIS COULD BE A CONSTANT <<<<<<<<<<<<<<<<
GNU_CONST inline double CalculateCoreMass_Luminosity_D_Hurley2000(const double p_Mass, const double p_ZetaHurley = 0.0, const double p_MHeF = 0.0) const override { return 5.5E4 / (1.0 + (0.4 * utils::IntPow(p_Mass, 4))); } // Hurley et al. 2000, just after eq 84
GNU_CONST inline double CalculateCoreMass_Luminosity_p_Hurley2000(const double p_Mass, const double p_MHeF) const override { return 5.0; }
GNU_CONST inline double CalculateCoreMass_Luminosity_q_Hurley2000(const double p_Mass, const double p_MHeF) const override { return 3.0; }
//------------------


double CalculateAgeAfterMassLoss() const override;
GNU_CONST double HeMS::CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tHeMS) const;


GNU_CONST static double CalculatePhaseLifetime_Hurley2000_Static(const double p_Mass);


GNU_CONST static double CalculateLuminosityAtPhaseEnd_Hurley2000_Static(const double p_Mass);
GNU_CONST static double CalculateLuminosityAtZAHeMS_Hurley2000_Static(const double p_Mass);




    COMPAS_PURE static MASS_LOSS_T CalculateMLrate_Belczynski2010_Static(const double p_Luminosity);

    COMPAS_PURE MASS_LOSS_T CalculateMLrate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_PerturbationMu = 0.0) const override;

COMPAS_PURE static MASS_LOSS_T CalculateMLrate_Merritt2025_Static(const double p_Luminosity, const double p_Temperature);


GNU_CONST static MASS_LOSS_T CalculateMLrate_Vink2017_Static(const double, p_Luminosity, const double p_ZetaAnders);

GNU_CONST static MASS_LOSS_T CalculateMLrateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature, const double p_SigmaHurley);


GNU_CONST static double CalculateRadiusAtZAHeMS_Hurley2000_Static(const double p_Mass);

GNU_CONST inline double CalculateRadius_Hurley2000() const override { 
    return CalculateRadius_Hurley2000_Static(m_StateHistory.CurrentState.Mass(), m_StateHistory.CurrentState.Tau());
}
GNU_CONST static double CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Tau);






 



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; } // McCO = 0.0 for HeMS stars

inline double CalculateHeCoreMassAtPhaseEnd() const override { return m_StateHistory.CurrentState.Mass(); } // McHe(HeMS) = Mass


GNU_CONST inline double CalculateRadiusAtPhaseEnd_Hurley2000() const override {
    return CalculateRadius_Hurley2000_Static(m_StateHistory.CurrentState.Mass(), 1.0); // tau = 1.0 at end of phase
}





            // Abundances



GNU_CONST inline double CalculateHAbundanceCore(const double p_Tau) const override { return 0.0; } // No hydrogen in the core for HeMS stars
double inline CalculateHAbundanceSurface(const double p_Tau) const override { return HAbundanceSurface(); }

COMPAS_PURE double CalculateHeAbundanceCore(const double p_Tau) const override;
double inline CalculateHeAbundanceSurface(const double p_Tau) const override { return HeAbundanceSurface(); }


            double          CalculateHeliumAbundanceCoreAtPhaseEnd() const                                          { return CalculateHeliumAbundanceCoreOnPhase(); }

            double          CalculateHeliumAbundanceSurfaceAtPhaseEnd() const                                       { return CalculateHeliumAbundanceSurfaceOnPhase(); }

            double          CalculateHydrogenAbundanceCoreAtPhaseEnd() const                                        { return CalculateHydrogenAbundanceCoreOnPhase(); } 

            double          CalculateHydrogenAbundanceSurfaceAtPhaseEnd() const                                     { return CalculateHydrogenAbundanceSurfaceOnPhase(); } 

            double          CalculateInitialSupernovaMass() const                                                   { return GiantBranch::CalculateInitialSupernovaMass(); }                        // Use GiantBranch
            double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const                                { return CalculateLuminosityAtPhaseEnd_Static(p_Mass); }
            double          CalculateLuminosityAtPhaseEnd() const                                                   { return CalculateLuminosityAtPhaseEnd(m_Mass); }                               // Use class member variables


            double          CalculateMassLossRateHurley();







            double          CalculateMomentOfInertia() const                                                        { return MainSequence::CalculateMomentOfInertia(); }


GNU_CONST inline double CalculateHurleyPerturbationMu() const { return 5.0; } // Hurley et al. 2000, eqs 97 & 98


            double          CalculateRadiusAtPhaseEnd(const double p_Mass) const                                    { return CalculateRadiusAtPhaseEnd_Static(p_Mass); }
            double          CalculateRadiusAtPhaseEnd() const                                                       { return CalculateRadiusAtPhaseEnd(m_Mass); }                                   // Use class member variables
    static  double          CalculateRadiusAtPhaseEnd_Static(const double p_Mass);
            double          CalculateRadiusOnMassChange(double p_dM)                                                { return CalculateRadiusOnPhaseTau(m_Mass + p_dM, m_Tau); }
            double          CalculateRadiusOnPhaseTau(const double p_Mass, const double p_Tau) const                { return CalculateRadiusOnPhase_Static(p_Mass, p_Tau); }
            double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhaseTau(m_Mass, m_Tau); }                            // Use class member variables
            double          CalculateRadiusOnPhase(double p_Mass, double p_Luminosity) const                        { return CalculateRadiusOnPhase(); }        // not a meaningful calculation for HeMS star, ignore arguments



inline double CalculateTau_Hurley2000() const override { return CalculateTau_Hurley2000(Age(), Timescales(TIMESCALE::tHeMS)); }
GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const double p_tHeMS) const;



            double          CalculateTemperatureAtPhaseEnd() const                                                  { return BaseStar::CalculateTemperatureAtPhaseEnd(); }
            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return CalculateTemperatureOnPhase(p_Luminosity, p_Radius); }                 // Same as on phase

inline double CalculateMLrateThermal() const override { return BaseStar::CalculateMLrateThermal(); } // Use BaseStar


            double          ChooseTimestep(const double p_Time) const;

            STELLAR_TYPE    EvolveToNextPhase();

GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::RADIATIVE; } // Always RADIATIVE for HeMS stars







            bool            IsEndOfPhase() const                                                                    { return !ShouldEvolveOnPhase(); }
            bool            IsSupernova() const                                                                     { return false; }                                                               // Not here

            void            PerturbLuminosityAndRadius() { }                                                                                                                                        // NO-OP

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
            void            ResolveHeliumFlash() { }                                                                                                                                                // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                                   { return m_StellarType; }                                                       // NO-OP

            void            SetSNHydrogenContent()                                                                  { m_SupernovaDetails.isHydrogenPoor = true; }                                   // Always true

            bool            ShouldEnvelopeBeExpelledByPulsations() const                                            { return false; }                                                               // No envelope to lose by pulsations
            bool            ShouldEvolveOnPhase() const                                                             { return (utils::Compare(m_Tau, 0.0) >= 0 && utils::Compare(m_Tau, 1.0) < 0); } // Evolve on HeMS phase if 0 <= tau < 1.0
            bool            ShouldSkipPhase() const                                                                 { return false; }                                                               // Never skip HeMS phase








            void            UpdateAgeAfterMassLoss();                                                                                                                                               // Per Hurley et al. 2000, section 7.1
};





/// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateHeAbundanceCore
 *
 * @brief
 * Calculate the helium abundance in the core of the star, given the phase relative age of
 * the star.  Currently just a simple linear model from the initial helium abundance to the
 * maximum helium abundance (assuming that all hydrogen is converted to helium). 
 * 
 * Should one day be updated to match detailed models.
 *
 * double CalculateHeAbundanceCore(const double p_Tau)
 * 
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      Helium abundance in the core of the star
 */
COMPAS_PURE inline double HeMS::CalculateHeAbundanceCore(const double p_Tau) const {
    return (1.0 - GLOBALS->Metallicity()) * (1.0 - p_Tau);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LIFETIME / AGE FUNCIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
inline double HeMS::CalculateAgeAfterMassLoss() const {

    double age;

    Switch (OPTIONS->Mode()) {                                                                                  // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                                             // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                                             // HURLEY BSE
            age = CalculateAgeAfterMassLoss_Hurley2000(m_StateHistory.CurrentState.Mass(),
                                                       m_StateHistory.CurrentState.Age(),
                                                       m_StateHistory.CurrentState.Timescales(tHeMS));
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
 * CalculateAgeAfterMassLoss_Hurley2000
 *
 * @brief
 * Recalculate the star's age after mass loss, per Hurley et al. 2000, section 7.1
 *
 *
 * double CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tHeMS) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tHeMS                         HeMS lifetime (Myr)
 * @return                                      Age of the star after mass loss (Myr)
 */
GNU_CONST inline double HeMS::CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tHeMS) const {
    return p_Age * CalculatePhaseLifetime_Hurley2000_Static(p_Mass) / p_tHeMS;
}


/*
 * CalculatePhaseLifetime_Hurley2000_Static
 *
 * @brief
 * Calculate lifetime of the Helium Main Sequence, tHeMS, per Hurley et al. 2000, eq 79
 *
 *
 * static double CalculatePhaseLifetime_Hurley2000_Static(const double p_Mass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      HeMS lifetime, tHeMS (Myr)
 */
 GNU_CONST inline double HeMS::CalculatePhaseLifetime_Hurley2000_Static(const double p_Mass) {

    // pow() is slow - use multiplication (sqrt() is much faster than pow())
    const double m4   = utils::IntPow(p_Mass, 4);
    const double m6   = p_Mass * p_Mass * m4;
    const double m6_5 = m6 * std::sqrt(p_Mass);

    return (0.4129 + (18.81 * m4) + (1.853 * m6)) / m6_5;
}


/*
 * CalculateTau_Hurley2000
 *
 * @brief
 * Calculate the HeNS-relative age (fractional Helium Main Sequence age) of the star,
 * per Hurley et al. 2000, just after eq 81.
 *
 * 
 * double CalculateTau_Hurley2000(const double p_Age, const double p_tHeMS) const
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tHeMS                         HeMS lifetime (Myr)
 * @return                                      HeMS-relative age, [0, 1]
 */
GNU_CONST inline double HeMS::CalculateTau_Hurley2000(const double p_Age, const double p_tHeMS) const {
    return std::max(0.0, std::min(1.0, p_Age / p_tHeMS));
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
 * DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_tScales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
COMPAS_PURE inline DBL_VECTOR HeMS::CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_tScales) const {

    DBL_VECTOR tScales = p_tScales; // copy given timescales

    tScales[static_cast<int>(TIMESCALE::tHeMS)] = CalculateLifetimeOnPhase_Hurley2000_Static(p_Mass);

    // return timescales vector by value - NRVO takes care of performance/efficiency
    return tScales;
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               LUMINOSITY FUNCTIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateLuminosity
 *
 * @brief
 * Calculate luminosity for a Helium Main Sequence star (during central He burning).
 *
 * Calls relevant luminosity function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateLuminosity() const
 *
 * @return                                      HeMS-relative age, [0, 1]
 */
inline double HeMS::CalculateLuminosity() const { 

    double luminosity;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            luminosity = CalculateLuminosity_Hurley2000_Static(m_StateHistory.CurrentState.Age(), m_StateHistory.CurrentState.Tau());
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

    return luminosity;
}


/*
 * CalculateLuminosity_Hurley2000_Static
 *
 * @brief
 * Calculate luminosity for a Helium Main Sequence star (during central He burning),
 * per Hurley et al. 2000, eqs 80 & 82
 *
 *
 * static double CalculateLuminosity_Hurley2000_Static(const double Mass, const double p_Tau)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      HeMS luminosity (Lsol)
 */
static GNU_CONST inline double HeMS::CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Tau) {
    return CalculateLuminosityAtZAHeMS_Hurley2000_Static(p_Mass) * (1.0 + (p_Tau * (0.45 + (std::max(0.0, 0.85 - 0.08 * p_Mass) * p_Tau))));
}


/*
 * CalculateLuminosityAtPhaseEnd_Hurley2000_Static
 *
 * @brief
 * Calculate luminosity at the end of the Helium Main Sequence, per Hurley et al. 2000, eq 80
 *
 *
 * static double CalculateLuminosityAtPhaseEnd_Hurley2000_Static(const double p_Mass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAHeMS luminosity (Lsol)
 */
GNU_CONST static inline double HeMS::CalculateLuminosityAtPhaseEnd_Hurley2000_Static(const double p_Mass) {
    return CalculateLuminosityAtZAHeMS_Hurley2000_Static(p_Mass) * (1.0 + 0.45 + std::max(0.0, 0.85 - (0.08 * p_Mass)));
}


/*
 * CalculateLuminosityAtZAHeMS_Hurley2000_Static
 *
 * @brief
 * Calculate luminosity at zero-age Helium Main Sequence, per Hurley et al. 2000, eq 77
 * 
 * 
 * static double CalculateLuminosityAtZAHeMS_Hurley2000_Static(const double p_MZAHeMS)
 *
 * @param       p_MZAHeMS                       ZAHeMS mass of the star (Msol)
 * @return                                      ZAHeMS luminosity (Lsol)
 */
GNU_CONST static inline double HeMS::CalculateLuminosityAtZAHeMS_Hurley2000_Static(const double p_MZAHeMS) {

    // pow() is slow - use multiplication (sqrt() is much faster than pow())
    const double m0_5   = std::sqrt(p_MZAHeMS);
    const double m3     = p_MZAHeMS * p_MZAHeMS * p_MZAHeMS;
    const double m6     = m3 * m3;
    const double m7_5   = p_MZAHeMS * m6 * m0_5;
    const double m9     = m6 * m3;
    const double m10_25 = p_MZAHeMS * m9 * std::sqrt(m0_5);

    return (15262.0 * m10_25) / (m9 + (29.54 * m7_5) + (31.18 * m6) + 0.0469);
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 RADIUS FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRadiusAtZAHeMS_Hurley2000_Static
 *
 * @brief
 * Calculate radius at zero-age Helium Main Sequence, per Hurley et al. 2000, eq 78
 *
 *
 * static double CalculateRadiusAtZAHeMS_Hurley2000_Static(const double p_MZAHeMS)
 *
 * @param       p_MZAHeMS                       ZAHeMS mass of the star (Msol)
 * @return                                      ZAHeMS radius (Rsol)
 */
GNU_CONST static inline double HeMS::CalculateRadiusAtZAHeMS_Hurley2000_Static(const double p_MZAHeMS) {
    const double m3 = p_MZAHeMS * p_MZAHeMS * p_MZAHeMS; // pow() is slow - use multiplication
    return (0.2391 * PPOW(p_MZAHeMS, 4.6)) / (p_MZAHeMS * m3 + (0.162 * m3) + 0.0065);
}


/*
 * CalculateRadius_Hurley2000_Static
 *
 * @brief
 * Calculate radius for a Helium Main Sequence star (during central He burning),
 * per Hurley et al. 2000, eq 81
 *
 *
 * static double CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Tau)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      Radius for a Helium Main Sequence star (Rsol)
 */
GNU_CONST static inline double HeMS::CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Tau) {
    return p_Mass <= 0.0 // sanity check
            ? 0.0
            : CalculateRadiusAtZAHeMS_Hurley2000_Static(p_Mass) * (1.0 + (std::max(0.0, 0.4 - 0.22 * log10(p_Mass)) * (p_Tau - utils::IntPow(p_Tau, 6))));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  MASS FUNCTIONS                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                             MASS LOSS RATE FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateMLrate_Belczynski2010_Static
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, per Belczynski 2010
 * (as implemented in StarTrack - courtesy Chris Belczynski).
 *
 * Since this is a static function and is not derived from BaseStar::CalculateMLrate_Belczynski2010(),
 * there is no requirement that the parameter list match that of BaseStar::CalculateMLrate_Belczynski2010()
 * 
 *
 * static MASS_LOSS_T CalculateMLrate_Belczynski2010_Static(const double p_Luminosity)
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
GNU_CONST static inline MASS_LOSS_T HeMS::CalculateMLrate_Belczynski2010_Static(const double p_Luminosity) {
    return BaseStar::CalculateMLrateWR_ZDependent_Static(p_Luminosity, 0.0);
}


/*
 * CalculateMLrate_Vink2017_Static
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for low-mass helium stars,
 * per Vink 2017 eq 1.
 * 
 * (See https://ui.adsabs.harvard.edu/abs/2017A%26A...607L...8V/abstract)
 * 
 * 
 * static MASS_LOSS_T CalculateMLrate_Vink2017_Static(const double p_Luminosity)
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
COMPAS_PURE static inline MASS_LOSS_T HeMS::CalculateMLrate_Vink2017_Static(const double p_Luminosity) {
    return std::make_tuple(PPOW(10.0, -13.3 + 1.36 * log10(p_Luminosity) + 0.61 * GLOBALS->ZetaAnders()), MASS_LOSS_TYPE:WR);
}


/*
 * CalculateMLrateWR_Shenar2019_Static
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for Wolf-Rayet stars,
 * per Shenar et al. 2019, eq 6, tbl 5.
 * 
 * (See https://ui.adsabs.harvard.edu/abs/2019A%26A...627A.151S/abstract)
 * 
 * Here we use the fitting coefficients for hydrogen poor WR stars (X_H < 0.05).
 * The C4 (X_He) term is = 0 and is omitted.
 *  
 * 
 * static MASS_LOSS_T CalculateMLrateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature)
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
COMPAS_PURE static inline MASS_LOSS_T HeMS::CalculateMLrateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature) {

    constexpr double C1 = -7.99;
    constexpr double C2 =  0.97;
    constexpr double C3 = -0.07;
    constexpr double C5 =  0.89;

    return std::make_tuple(PPOW(10.0, C1 + (C2 * log10(p_Luminosity)) + (C3 * log10(p_Temperature * TSOL)) + (C5 * GLOBALS->SigmaHurley())), MASS_LOSS_TYPE::WR);
}



/// HeMS_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


class HeMS_Constituent: virtual public BinaryConstituentStar, public HeMS {

public:


protected:

};



    MT_CASE         DetermineMassTransferTypeAsDonor() const                                                        { return MT_CASE::OTHER; }                                                      // Not A, B, C, or NONE

            double          CalculateMTRejuvenationFactor();


double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const;

double CalculateCriticalMassRatio_Hurley2002() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_GT_07; }




            double          CalculateCELambda_Dewi() const                                                             { return 0.5; }
            double          CalculateLambdaLoveridge(const double p_EnvMass, const bool p_IsMassLoss = false) const { return BaseStar::CalculateLambdaLoveridge(p_EnvMass, p_IsMassLoss); }   // Not supported - use BaseStar
            double          CalculateLambdaNanjingStarTrack(const double p_Mass) const  { return BaseStar::CalculateLambdaNanjingStarTrack(0.0, 0.0); }                 // Not supported - use BaseStar (0.0 are dummy values)

            

            double          CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription)                  { return OPTIONS->ZetaMainSequence(); }                                                                             // A HeMS star is treated as any other MS star for Zeta calculation purposes
            double          CalculateZetaEquilibrium()                                                              { return MainSequence::CalculateZetaEquilibrium(); }                           // A HeMS star is treated as any other MS star for Zeta calculation purposes



double CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency = 0.0) const;


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
inline double HeMS_Constituent::CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {                                                                                                          
    return p_AccretorIsDegenerate
            ? OPTIONS->MassTransferCriticalMassRatioHeMSDegenerateAccretor()        // degenerate accretor
            : OPTIONS->MassTransferCriticalMassRatioHeMSNonDegenerateAccretor();    // non-degenerate accretor
                                                                                                                        
    return qCrit;
}



#endif // __HeMS_h__
