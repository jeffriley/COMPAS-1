#ifndef __Remnants_h__
#define __Remnants_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "HeGB.h"


class BaseStar;
class HeGB;

class Remnants: virtual public BaseStar, public HeGB {

public:

    Remnants(){};
    
    Remnants(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeGB(p_BaseStar, false) { }

    
    // member functions
    
    



protected:


    // member functions - alphabetically





///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

DBL_VECTOR CalculateTimescales_Hurley2000() const override { return TPAGB::CalculateTimescales_Hurley2000(MassEffectiveInitial(), m_InterimState.GBparams(), m_InterimState.TimeScales()); }

inline DBL_VECTOR CalculateGBparams_Hurley2000() const override { return m_InterimState.GBparams(); } // gb params not used beyond helium stars


inline double CalculateCOCoreMass() const override { return Mass(); } // McCO = M for remnants

inline double CalculateHeCoreMass() const override { return Mass(); } // McHe = M for remnants

    double          CalculateConvectiveCoreRadius() const                                                       { return m_Radius; }                                                    // All core

inline Dbl_DblT CalculateConvectiveEnvelopeMass() const override { return std::make_tuple(0.0, 0.0); } // No envelope for remnants

    double          CalculateCoreMassOnPhase() const                                                            { return m_Mass; }                                                      // Return m_Mass


GNU_CONST inline Dbl_DblT CalculateRemnantMass(const double p_COCoreMass) { return CalculateRemnantMass_Static(const double p_COCoreMass); }
GNU_CONST static Dbl_DblT CalculateRemnantMass_Static(const double p_COCoreMass);
GNU_CONST static Dbl_DblT CalculateRemnantMass_Hurley2000_Static(const double p_COCoreMass);


inline double CalculateRemnantRadius() const override { return Radius(); }


GNU_CONST inline double CalculateHAbundanceCore(const double p_Tau) const override       { return 0.0; }; // No hydrogen in the core for remnants
GNU_CONST inline double CalculateHAbundanceSurface(const double p_Tau) const override    { return 0.0; }; // No hydrogen on the surface for remnants
GNU_CONST inline double CalculateHeAbundanceCore(const double p_Tau) const override    { return 0.0; }; // No helium in the core for remnants (except HeWD)
GNU_CONST inline double CalculateHeAbundanceSurface(const double p_Tau) const override { return 0.0; }; // No helium on the surface for remnants (except HeWD)
    
    


    double          CalculateInitialSupernovaMass() const                                                       { return GiantBranch::CalculateInitialSupernovaMass(); }                // Use GiantBranch


    GNU_CONST double          CalculateMassLossRateHurley()                                                               { m_DominantMassLossRate = MASS_LOSS_TYPE::NONE ; return 0.0; }

    GNU_CONST MASS_LOSS_T CalculateMassLossRateBelczynski2010(const double p_Luminosity, const double p_HeAbundanceSurface) { return make_tuple(0.0, MASS_LOSS_TYPE::NONE); }

    GNU_CONST double          CalculateMassLossRateMerritt2025()                                                          { m_DominantMassLossRate = MASS_LOSS_TYPE::NONE ; return 0.0; }                                                         // 


inline double CalculateHurleyPerturbationMu() const { return m_Mu; }


GNU_CONST double CalculateConvectiveEnvelopeRadialExtent() const override { return 0.0; } // WD stars don't have a convective envelope



    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                              { return BaseStar::CalculateRadiusAndStellarTypeOnPhase(); }
   

inline double CalculateTimescale_Thermal() const override { return CalculateTimescale_Dynamical(); }

inline double CalculateMLRateThermal() const override{ return BaseStar::CalculateMLRateThermal(); } // Set thermal mass gain rate to be effectively infinite, using dynamical timescale (in practice, will be Eddington limited), avoid division by zero



    GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const override;

    
GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::REMNANT; } // Always REMNANT for remnants


virtual double CalculateEddingtonCriticalRate() { return CalculateEddingtonCriticalRate(Radius()); }
virtual double CalculateEddingtonCriticalRate(const double p_Radius) const { return 2.08E-3 / 1.7 * p_Radius * MYR_TO_YEAR * OPTIONS->EddingtonAccretionFactor() ; } // Hurley+, 2002, Eq. (67)





///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inline double CalculateCOCoreMassAtPhaseEnd() const override { return Mass(); } // McCO = M for remnants

inline double CalculateHeCoreMassAtPhaseEnd() const override { return Mass(); } // McHe = M for remnants









    STELLAR_TYPE    EvolveToNextPhase()                                                                         { return BaseStar::EvolveToNextPhase(); }                               // Default to BaseStar

    bool            IsDegenerate() const                                                                        { return true; }                                                        // White Dwarfs, NS and BH are degenerate

    bool            IsEndOfPhase() const                                                                        { return !ShouldEvolveOnPhase(); }                                      // Phase ends when envelope loss or going supernova

    bool            IsSupernova() const                                                                         { return false; }                                                       // Default

    void            PerturbLuminosityAndRadius() { }                                                                                                                                    // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false)                                                   { return BaseStar::ResolveEnvelopeLoss(p_Force); }                      // Default to BaseStar

    void            ResolveEnvelopeMassAtPhaseEnd(const double p_Tau) const                                     { ResolveEnvelopeMassOnPhase(p_Tau); }                                  // Same as on phase
    void            ResolveEnvelopeMassOnPhase(const double p_Tau) const { }                                                                                                            // NO-OP

    STELLAR_TYPE    ResolveMassLoss(const double p_dt)                                                          { return m_StellarType; }                                               // NO-OP                                                                                                                                    // NO-OP
    STELLAR_TYPE    ResolveMassLossHurley(const double p_dt)                                                    { return ResolveMassLoss(p_Dt); }                                                                                                                                   // NO-OP

    STELLAR_TYPE    ResolveSkippedPhase()                                                                       { return BaseStar::ResolveSkippedPhase(); }                             // Default to BaseStar
                                                                                                                                                                                        //
    void            ResolveShellChange(const double p_AccretedMass) { }                                                                                                                 // NO-OP 
                                                                                                                                                                                        //
    STELLAR_TYPE    ResolveSupernova()                                                                          { return BaseStar::ResolveSupernova(); }                                // Default to BaseStar

    void            SetPulsarParameters() const { }                                                                                                                                     // NO-OP

    bool            ShouldEnvelopeBeExpelledByPulsations() const                                                { return false; }                                                       // No envelope to lose by pulsations
    bool            ShouldEvolveOnPhase() const                                                                 { return true; }                                                        // Default
    bool            ShouldSkipPhase() const                                                                     { return false; }                                                       // Don't skip WD phase

};



/// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



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
 * @param       p_tScales                       Phase timescales (Myr) (not used here)
 * @return                                      Suggested timestep (Myr)
 */
GNU_CONST inline double Remnants::ChooseTimestep_Hurley2000(const double p_Age, [[maybe_unused]] const DBL_VECTOR& p_tScales) const {

    double dtk = std::min(std::max(1.0, p_Age), 500.0); // modfied from Hurley
                                                                                    
    // time to end of phase (change of stellar type, dte) not used here - how to calculate?
    // clamp to minimum NUCLEAR_MINIMUM_TIMESTEP
    return std::max(dtk, NUCLEAR_MINIMUM_TIMESTEP);
}



/*
 * CalculateRemnantMass_Hurley2000_Static
 *
 * @brief
 * Calculate remnant mass, per Hurley et al. 2000, eq 92.
 * 
 * Hurley et al. 2000 indicates a star with a CO core > 7.0 Msol leads to a collapse into a BH,
 * but is ambiguous about the BH mass in this case - we assume here that it's just the CO core.
 *
 * 
 * static Dbl_DblT CalculateRemnantMass_Hurley2000_Static(const double p_COCoreMass)
 *
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
GNU_CONST static inline Dbl_DblT Remnants::CalculateRemnantMass_Hurley2000_Static(const double p_COCoreMass) {
    return std::make_tuple((p_COCoreMass > 7.0 ? p_COCoreMass : 1.17 + (0.09 * p_COCoreMass)), 0.0); // fallback fraction not defined by Hurley
}


/*
 * CalculateRemnantMass_Static
 *
 * @brief
 * Calculate remnant mass.
 * 
 * Calls relevant remnant mass function based on the evolutionary mode given in program options.
 * 
 * 
 * static Dbl_DblT CalculateRemnantMass_Static(const double p_COCoreMass)
 *
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
GNU_CONST static inline Dbl_DblT Remnants::CalculateRemnantMass_Static(const double p_COCoreMass) {

    double mass;
    double fallbackFraction;

    Switch (OPTIONS->Mode()) {                                              // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                         // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                         // HURLEY BSE
            std::tie(mass, fallbackFraction) = CalculateRemnantMass_Hurley2000_Static(const double p_COCoreMass);
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

    return std::make_tuple(mass, fallbackFraction);
}






/// Remnants_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


class Remnants_Constituent: virtual public BinaryConstituentStar, public Remnants {

public:


protected:

};


    





    double          CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription)                      { return 0.0; }                                                         // Should never be called...
    
    double          CalculateZetaEquilibrium()                                      { return -std::numeric_limits<double>::infinity(); }                                         // Nuclear timescale MT should be impossible from remnant stars; should never be called



    double          CalculateCELambda_Dewi() const                                                                 { return BaseStar::CalculateCELambda_Dewi(); }
    double          CalculateLambdaNanjingStarTrack(const double p_Mass) const      { return BaseStar::CalculateLambdaNanjingStarTrack(0.0, 0.0); }

    Dbl_DblT         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const double p_AccretorMassRate);
    Dbl_DblT         CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                const double p_AccretorMassRate,
                                                const bool   p_IsHeRich)                                        { return CalculateMassAcceptanceRate(p_DonorMassRate, p_AccretorMassRate); } // Ignore the He content for non-WDs


                                                
    Dbl_Dbl_Dbl_DblT CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const   { return std::make_tuple(0.0, 0.0, 0.0, 0.0); }              // Default is no tidal response
    Dbl_Dbl_Dbl_DblT CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const { return std::make_tuple(0.0, 0.0, 0.0, 0.0); }              // Default is no tidal response


double CalculateCriticalMassRatio(const double p_Mass,
                                  const double p_Radius,
                                  const double p_CoreMass,
                                  const bool   p_AccretorIsDegenerate,
                                  const double p_MTefficiency = 0.0) const { return 0.0; } // Should not be called (but if it is, mass transfer from a neutron star always treated as stable)

double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const { return 0.0; }

double CalculateCriticalMassRatio_Ge2020(const double p_MTefficiency) { return 0.0; }

double CalculateCriticalMassRatio_Hurley2002() const { return 0.0; }



#endif // __Remnants_h__
