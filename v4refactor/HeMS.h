#pragma once

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

    HeMS() { m_InterimState.SetStellarType(STELLAR_TYPE::NAKED_HELIUM_STAR_MS); };
    
    HeMS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), TPAGB(p_BaseStar, false) {
        m_InterimState.SetStellarType(STELLAR_TYPE::NAKED_HELIUM_STAR_MS);                                                                                                                                         // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                                             // Initialise if required
    }


    // member functions - alphabetically














    static Dbl_DblT  CalculateRadiusAtPhaseEnd_Static(const double p_Mass, const double p_Luminosity);




    



protected:

    void Initialise() {
        // initialise surface abundances
        m_InterimState.SetHAbundanceSurface(0.0);
        m_InterimState.SetHeAbundanceSurface(1.0 - Metallicity());
        
        CalculateTimescales();
        // JR: Age for HeMS is partially calculated before switching -
        // can get here from various places in ResolveEnvelopeLoss(),
        // and Age is calculated differently in those cases
        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically
GNU_CONST inline Dbl_DblT CalculateConvectiveEnvelopeMass() const override { return std::make_tuple(0.0, 0.0); } // No convective envelope for HeMS stars


            double          CalculateConvectiveCoreMass() const;
            double          CalculateConvectiveCoreRadius() const {  } // JR FIX THIS: placeholder
            double          CalculateConvectiveCoreRadius(const double p_Radius, const double p_Tau) const;





///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


void CalculateTimescales_Hurley2000() override { m_InterimState.SetHurleyTimescales(CalculateTimescales_Hurley2000(MassEffectiveInitial(), m_InterimState.HurleyTimescalesOrDefault())); }
GNU_CONST inline TimescalesT CalculateTimescales_Hurley2000(const double p_Mass, const TimescalesT& p_tScales) const;

inline DblVectorT CalculateGBParams_Hurley2000(const double p_Mass, DblVectorT p_GBParams) const override {
    return GiantBranch::CalculateGBParams_Hurley2000(m_InterimSate.Mass0(), m_InterimSate.GBParams());
}



inline double CalculateRemnantRadius() const override { return Radius(); }









inline double CalculateLuminosity() const override;


// MAYBE BREAK THE DESCRIPTION OF THESE OUT TO INLINE BELOW
GNU_CONST inline double CalculateCoreMass_Luminosity_B_Hurley2000(const double p_Mass) const override { return 4.1E4; } // Hurley et al. 2000, just after eq 84  THIS COULD BE A CONSTANT <<<<<<<<<<<<<<<<
GNU_CONST inline double CalculateCoreMass_Luminosity_D_Hurley2000(const double p_Mass) const override { return 5.5E4 / (1.0 + (0.4 * utils::intPow(p_Mass, 4))); } // Hurley et al. 2000, just after eq 84
GNU_CONST inline double CalculateCoreMass_Luminosity_p_Hurley2000(const double p_Mass) const override { return 5.0; }
GNU_CONST inline double CalculateCoreMass_Luminosity_q_Hurley2000(const double p_Mass) const override { return 3.0; }
//------------------





    //////////////////////////////////////////////////
    //   ABUNDANCE                                  //
    //////////////////////////////////////////////////



GNU_CONST inline double CalculateHAbundanceCore(const double p_Tau) const override { return 0.0; } // No hydrogen in the core for HeMS stars
double inline CalculateHAbundanceSurface(const double p_Tau) const override { return HAbundanceSurface(); }

COMPAS_PURE double CalculateHeAbundanceCore(const double p_Tau) const override;
double inline CalculateHeAbundanceSurface(const double p_Tau) const override { return HeAbundanceSurface(); }


            double          CalculateHeAbundanceCoreAtPhaseEnd() const                                          { return CalculateHeAbundanceCoreOnPhase(Tau()); }

            double          CalculateHeAbundanceSurfaceAtPhaseEnd() const                                       { return CalculateHeAbundanceSurfaceOnPhase(Tau()); }

            double          CalculateHAbundanceCoreAtPhaseEnd() const                                        { return CalculateHAbundanceCoreOnPhase(Tau()); } 

            double          CalculateHAbundanceSurfaceAtPhaseEnd() const                                     { return CalculateHAbundanceSurfaceOnPhase(Tau()); } 




    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////

    inline double CalculateAgeAfterMassLoss_Hurley2000() const override {
        return CalculateAgeAfterMassLoss_Hurley2000(Mass(), Age(), Timescales(TS::HeMS));
    }
    GNU_CONST double CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tHeMS) const;


    GNU_CONST static double CalculatePhaseLifetime_Hurley2000(const double p_Mass); // JR FIX THIS: DONE




    GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const TimescalesT& p_tScales); // JR FIX THIS: DONE
    GNU_CONST double CalculateTauAtPhaseEnd_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const override { return 1.0; } // 1.0 at phase end JR FIX THIS: DONE



GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const override;






    //////////////////////////////////////////////////
    //   ENVELOPE                                   //
    //////////////////////////////////////////////////

GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::RADIATIVE; } // Always RADIATIVE for HeMS stars







    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    GNU_CONST static double CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Tau); // JR FIX THIS: DONE
    GNU_CONST double CalculateLuminosity_Hurley2000() const override { return CalculateLuminosity_Hurley2000(Mass(), Tau()); } // JR FIX THIS: DONE


    GNU_CONST static double CalculateLuminosityAtPhaseEnd_Hurley2000(const double p_Mass); // JR FIX THIS: DONE
    GNU_CONST double CalculateLuminosityAtPhaseEnd_Hurley2000() const override { return CalculateLuminosityAtPhaseEnd_Hurley2000(Mass()); } // JR FIX THIS: DONE


    GNU_CONST static double CalculateLuminosityAtZAHeMS_Hurley2000(const double p_Mass); // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////

    GNU_CONST inline double CalculateCoreMass_Hurley2000() const override { return 0.0; } // Mc(HeMS) = 0.0 JR FIX THIS: DONE
    GNU_PURE  inline double CalculateCoreMassAtPhaseEnd_Hurley2000() const override { return CalculateHeCoreMassAtPhaseEnd(); } // /*Ilya*/ To fix, not everything will become CO core // JR FIX THIS: DONE

    GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; } // McCO(HeMS) = 0.0 // JR FIX THIS: DONE
    GNU_PURE  inline double CalculateCOCoreMassAtPhaseEnd() const override { return CalculateCOCoreMass(); } // Same as on phase // JR FIX THIS: DONE

    GNU_PURE  inline double CalculateHeCoreMass() const override { return Mass(); } // McHe(HeMS) = Mass // JR FIX THIS: DONE
    GNU_PURE  inline double CalculateHeCoreMassAtPhaseEnd() const override { return CalculateHeCoreMass(); } // Same as on phase // JR FIX THIS: DONE




    GNU_PURE  inline double CalculateEffectiveInitialMass_Hurley2000() const override { return Mass(); } // per Hurley et al. 2000, section 7.1 // JR FIX THIS: DONE



GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; } // McCO = 0.0 for HeMS stars

inline double CalculateHeCoreMassAtPhaseEnd() const override { return Mass(); } // McHe(HeMS) = Mass


            double          CalculateInitialSupernovaMass() const                                                   { return GiantBranch::CalculateInitialSupernovaMass(); }                        // Use GiantBranch





    //////////////////////////////////////////////////
    //   MASS LOSS / ACCRETION                      //
    //////////////////////////////////////////////////


COMPAS_PURE MassLossT CalculateMLRate_Belczynski2010(const double p_Metallicity, const double p_Luminosity, double, double, double, double, double, double, double, LBV_ML_PRESCRIPTION) const override;

COMPAS_PURE MassLossT CalculateMLRate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, double) const override;

GNU_PURE static MassLossT CalculateMLRate_Merritt2025(
    const double p_Metallicity,
    const double p_Mass,
    const double p_Radius,
    const double p_Luminosity,
    const double p_Temperature,
    const double p_HeAbundanceSurface,
    const double p_Perturb,
    const double p_CoolWindMultiplier,
    const double p_LBVFactor,
    const double p_mStart
);


    GNU_CONST static MassLossT CalculateMLRate_Vink2017_Static(const double p_Luminosity, const double p_ZetaAnders);

COMPAS_PURE MassLossT CalculateMLRateWR_Shenar2019(const double p_Metallicity, const double p_Luminosity, const double p_Temperature) const override;



inline double CalculateMLRateThermal() const override { return BaseStar::CalculateMLRateThermal(); } // Use BaseStar






    //////////////////////////////////////////////////
    //   MISC. CONSTANTS / PARAMETERS               //
    //////////////////////////////////////////////////

GNU_CONST inline double CalculateHurleyPerturbationMu() const { return 5.0; } // Hurley et al. 2000, eqs 97 & 98






    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////



GNU_CONST inline double CalculateRadius_Hurley2000() const override { 
    return CalculateRadius_Hurley2000_Static(Mass(), Tau());
}


    GNU_CONST static double CalculateRadiusAtZAHeMS_Hurley2000(const double p_Mass); // JR FIX THIS: DONE

    GNU_CONST static double         CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Tau);

GNU_CONST inline double CalculateRadiusAtPhaseEnd_Hurley2000() const override { return CalculateRadius_Hurley2000_Static(Mass(), 1.0); } // tau = 1.0 at end of phase

double          CalculateRadiusAtPhaseEnd(const double p_Mass) const                                    { return CalculateRadiusAtPhaseEnd_Static(p_Mass); }
            double          CalculateRadiusAtPhaseEnd() const                                                       { return CalculateRadiusAtPhaseEnd(Mass()); }
    static  double          CalculateRadiusAtPhaseEnd_Static(const double p_Mass) {  } // JR FIX THIS: placeholder
            double          CalculateRadiusOnMassChange(double p_dM)                                                { return CalculateRadiusOnPhaseTau(Mass() + p_dM, Tau()); }
            double          CalculateRadiusOnPhaseTau(const double p_Mass, const double p_Tau) const                { return CalculateRadius_Hurley2000_Static(p_Mass, p_Tau); }
            double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhaseTau(Mass(), Tau()); }
            double          CalculateRadiusOnPhase(double p_Mass, double p_Luminosity) const                        { return CalculateRadiusOnPhase(); }        // not a meaningful calculation for HeMS star, ignore arguments


 







    //////////////////////////////////////////////////
    //   ROTATION                                   //
    //////////////////////////////////////////////////






            double          CalculateMomentOfInertia() const                                                        { return MainSequence::CalculateMomentOfInertia(); }





            


 







    //////////////////////////////////////////////////
    //   TEMPERATURE                                //
    //////////////////////////////////////////////////




            double          CalculateTemperatureAtPhaseEnd() const                                                  { return BaseStar::CalculateTemperatureAtPhaseEnd(); }
            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return CalculateTemperatureOnPhase(p_Luminosity, p_Radius); }                 // Same as on phase




            STELLAR_TYPE    EvolveToNextPhase();








            bool            IsEndOfPhase() const                                                                    { return !ShouldEvolveOnPhase(); }
            bool            IsSupernova() const                                                                     { return false; }                                                               // Not here

            void            PerturbLuminosityAndRadius() { }                                                                                                                                        // NO-OP
    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
            void            ResolveHeliumFlash() { }                                                                                                                                                // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                                   { return StellarType(); }

            void            SetSNHydrogenContent()                                                                  { m_InterimState.SetSupernovaDetailsIsHydrogenPoor(true); }                                   // Always true

            bool            ShouldEnvelopeBeExpelledByPulsations() const                                            { return false; }                                                               // No envelope to lose by pulsations
            bool            ShouldEvolveOnPhase() const                                                             { return (utils::Compare(InterimState().Tau(), 0.0) >= 0 && utils::Compare(InterimState().Tau(), 1.0) < 0); } // Evolve on HeMS phase if 0 <= new tau < 1.0 (transition phase reads InterimState)
            bool            ShouldSkipPhase() const                                                                 { return false; }                                                               // Never skip HeMS phase








            void            UpdateAgeAfterMassLoss() {  } // JR FIX THIS <<<<<<<<<<<<<<<<<<<< "Calculate" // Per Hurley et al. 2000, section 7.1
};





/// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ABUNDANCE                                     //
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
    return (1.0 - Metallicity()) * (1.0 - p_Tau);
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tHeMS                         HeMS lifetime (Myr)
 * @return                                      Age of the star after mass loss (Myr)
 */
inline double HeMS::CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tHeMS) const {
    return p_Age * CalculatePhaseLifetime_Hurley2000(p_Mass) / p_tHeMS;
}


/*
 * CalculatePhaseLifetime_Hurley2000
 *
 * @brief
 * Calculate lifetime of the Helium Main Sequence, tHeMS, per Hurley et al. 2000, eq 79
 *
 *
 * double CalculatePhaseLifetime_Hurley2000(const double p_Mass) const 
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      HeMS lifetime, tHeMS (Myr)
 */
 inline double HeMS::CalculatePhaseLifetime_Hurley2000(const double p_Mass) const {

    // pow() is slow - use multiplication (sqrt() is much faster than pow())
    const double m4   = utils::intPow(p_Mass, 4);
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
 * double CalculateTau_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      HeMS-relative age of the star, [0, 1]
 */
GNU_CONST inline double HeMS::CalculateTau_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const {
    return std::clamp(p_Age / p_tScales[TS::HeMS], 0.0, 1.0);
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
 * DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const TimescalesT& p_tScales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
COMPAS_PURE inline DBL_VECTOR HeMS::CalculateTimescales_Hurley2000(const double p_Mass, const TimescalesT& p_tScales) const {

    DBL_VECTOR tScales = p_tScales; // copy given timescales

    tScales[TS::HeMS] = CalculatePhaseLifetime_Hurley2000_Static(p_Mass);

    // Return timescales vector by value - NRVO takes care of performance/efficiency
    return tScales;
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
 * double ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales)
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Suggested timestep (Myr)
 */
GNU_CONST inline double HeMS::ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const {
    const double dtk = 0.05 * p_tScales[TS::HeMS];                  // Stellar type specific dt
    const double dte = p_tScales[TS::HeMS] - p_Age;                 // Time to end of phase (change of stellar type)

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);  // clamp to minimum NUCLEAR_MINIMUM_TIMESTEP
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateLuminosity_Hurley2000
 *
 * @brief
 * Calculate luminosity for a Helium Main Sequence star (during central He burning),
 * per Hurley et al. 2000, eqs 80 & 82
 *
 *
 * static double CalculateLuminosity_Hurley2000(const double Mass, const double p_Tau)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      HeMS luminosity of the star (Lsol)
 */
inline double HeMS::CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Tau) {
    return CalculateLuminosityAtZAHeMS_Hurley2000(p_Mass) * (1.0 + (p_Tau * (0.45 + (std::max(0.0, 0.85 - 0.08 * p_Mass) * p_Tau))));
}


/*
 * CalculateLuminosityAtPhaseEnd_Hurley2000
 *
 * @brief
 * Calculate luminosity at the end of the Helium Main Sequence,
 * per Hurley et al. 2000, eq 80
 *
 *
 * static double CalculateLuminosityAtPhaseEnd_Hurley2000(const double p_Mass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAHeMS luminosity (Lsol)
 */
inline double HeMS::CalculateLuminosityAtPhaseEnd_Hurley2000(const double p_Mass) {
    return CalculateLuminosityAtZAHeMS_Hurley2000(p_Mass) * (1.0 + 0.45 + std::max(0.0, 0.85 - (0.08 * p_Mass)));
}


/*
 * CalculateLuminosityAtZAHeMS_Hurley2000
 *
 * @brief
 * Calculate luminosity at zero-age Helium Main Sequence, per Hurley et al. 2000, eq 77
 * 
 * 
 * static double CalculateLuminosityAtZAHeMS_Hurley2000(const double p_MZAHeMS)
 *
 * @param       p_MZAHeMS                       ZAHeMS mass of the star (Msol)
 * @return                                      ZAHeMS luminosity of the star (Lsol)
 */
inline double HeMS::CalculateLuminosityAtZAHeMS_Hurley2000(const double p_MZAHeMS) {

    // pow() is slow - use multiplication
    // sqrt() is much faster than pow()
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
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               MASS LOSS / ACCRETION                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateMLRate_Belczynski2010
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, per Belczynski 2010
 * (as implemented in StarTrack - courtesy Chris Belczynski).
 * 
 *
 * MassLossT CalculateMLRate_Belczynski2010(const double p_Metallicity, double, double, const double p_Luminosity, double, double, double, double, double, LBV_ML_PRESCRIPTION) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * << 7 unnamed double parameters, and one unnamed LBV_ML_PRESCRIPTION parameter, to match BaseStar's virtual signature >>
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::WR)
 */
inline MassLossT HeMS::CalculateMLRate_Belczynski2010(const double p_Metallicity, double, double, const double p_Luminosity, double, double, double, double, double, LBV_ML_PRESCRIPTION) const {
    return CalculateMLRateWR_ZDependent(p_Metallicity, p_Luminosity, 0.0);
}


/*
 * CalculateMLRate_Vink2017_Static
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for low-mass helium stars,
 * per Vink 2017 eq 1.
 * 
 * (See https://ui.adsabs.harvard.edu/abs/2017A%26A...607L...8V/abstract)
 * 
 * 
 * MassLossT CalculateMLRate_Vink2017_Static(const double p_Luminosity)
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_ZetaAnders                    ZDEP->ZetaAnders(Metallicity()) at call site
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::WR)
 */
// JR FIX THIS <<<<<<<<<<<<<<<<< ZDEP???
COMPAS_PURE inline MassLossT HeMS::CalculateMLRate_Vink2017_Static(const double p_Luminosity, const double p_ZetaAnders) {
    return std::make_tuple(PPOW(10.0, -13.3 + 1.36 * std::log10(p_Luminosity) + 0.61 * p_ZetaAnders), ML_TYPE::WR);
}


/*
 * CalculateMLRateWR_Shenar2019_Static
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
 * MassLossT CalculateMLRateWR_Shenar2019(const double p_Metallicity, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                   double  WR mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::WR)
 */
inline MassLossT HeMS::CalculateMLRateWR_Shenar2019(const double p_Metallicity, const double p_Luminosity, const double p_Temperature) const {

    constexpr double C1 = -7.99;
    constexpr double C2 =  0.97;
    constexpr double C3 = -0.07;
    constexpr double C5 =  0.89;

    return std::make_tuple(PPOW(10.0, C1 + (C2 * std::log10(p_Luminosity)) + (C3 * std::log10(p_Temperature * TSOL)) + (C5 * ZDEP->SigmaHurley(p_Metallicity))), ML_TYPE::WR);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           MISC. CONSTANTS / PARAMETERS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateRadiusAtZAHeMS_Hurley2000
 *
 * @brief
 * Calculate radius at zero-age Helium Main Sequence, per Hurley et al. 2000, eq 78
 *
 *
 * static double CalculateRadiusAtZAHeMS_Hurley2000(const double p_MZAHeMS)
 *
 * @param       p_MZAHeMS                       ZAHeMS mass of the star (Msol)
 * @return                                      ZAHeMS radius of the star (Rsol)
 */
inline double HeMS::CalculateRadiusAtZAHeMS_Hurley2000(const double p_MZAHeMS) {
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
GNU_CONST inline double HeMS::CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Tau) {
    return p_Mass <= 0.0 // sanity check
            ? 0.0
            : CalculateRadiusAtZAHeMS_Hurley2000(p_Mass) * (1.0 + (std::max(0.0, 0.4 - 0.22 * std::log10(p_Mass)) * (p_Tau - utils::intPow(p_Tau, 6))));
}




/// HeMS_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//

#include "ConstituentStar.h"

class HeMS_Constituent: public ConstituentStar {

public:

    explicit HeMS_Constituent(Star* p_Star = nullptr) : ConstituentStar(p_Star) {}

    // MT classification
    MT_CASE DetermineMassTransferTypeAsDonor() const override { return MT_CASE::OTHER; }                       // Not A, B, C, or NONE

    // MT rejuvenation
    double  CalculateMTRejuvenationFactor() override;                                                          // body in HeMS.cpp

    // Critical mass ratio
    double  CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const override;           // body inline below
    double  CalculateCriticalMassRatio_Hurley2002() const override { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_GT_07; }
    double  CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency = 0.0) const override {  } // JR FIX THIS: placeholder

    // CE lambda family (HeMS-specific: not really supported, use safe defaults)
    double  CalculateCELambda_Dewi() const override { return 0.5; }
    double  CalculateLambdaLoveridge(const double p_EnvMass, const bool p_IsMassLoss = false) const override { (void)p_EnvMass; (void)p_IsMassLoss; return 1.0; }      // Not supported - safe default (was: BaseStar::CalculateLambdaLoveridge fallback)
    COMPAS_PURE double  CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const override { (void)p_Mass; (void)p_Radius; (void)p_CoreMass; return 1.0; }                                                          // Not supported - safe default

    // Zeta (HeMS treated as MS-like)
    double  CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription) override { (void)p_ZetaPrescription; return OPTIONS->ZetaMainSequence(); }      // A HeMS star is treated as any other MS star for Zeta calculation purposes

    double  CalculateZetaEquilibrium() override;

};


/*
 * CalculateCriticalMassRatio_Claeys2014
 *
 * @brief
 * Calculate the critical mass ratio, per Claeys et al. 2014
 * 
 * @param       p_AccretorIsDegenerate          Boolean indicating if accretor is degenerate
 * @return                                      Critical mass ratio for unstable MT 
 */
inline double HeMS_Constituent::CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {                                                                                                          
    return p_AccretorIsDegenerate
            ? OPTIONS->MassTransferCriticalMassRatioHeliumMSDegenerateAccretor()        // degenerate accretor
            : OPTIONS->MassTransferCriticalMassRatioHeliumMSNonDegenerateAccretor();    // non-degenerate accretor
}



