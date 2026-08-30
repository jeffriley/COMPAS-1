#pragma once

#include "constants.h"
#include "typedefs.h"
#include "utils.h"

#include "HeGB.h"


class BaseStar;
class HeGB;

class Remnants: virtual public BaseStar, public HeGB {


public:

    Remnants(){};
    
    Remnants(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeGB(p_BaseStar, false) { }


protected:

    // Member functions (not getters or setters)
    //
    // VIRTUAL FUNCTIONS may be (are expected to be) overridden by derived classes.
    // When overriding virtual functions in a derived class, use the "override" attribute.
    //
    // NON-VIRTUAL FUNCTIONS should not be overridden (declared separately) by derived classes.
    // While it is legal in C++ to declare the same (non-virtual) function in multiple classes,
    // (aka "shadowing", or "hiding"), we discourage it.  Non-virtual functions are statically
    // bound, and as such, especially with indirection, may not produce expected results.

    //////////////////////////////////////////////////
    //   ABUNDANCE                                  //
    //////////////////////////////////////////////////

    GNU_CONST double CalculateHAbundanceCore(const double p_Tau) const override { // JR FIX THIS: DONE
        return 0.0; // No hydrogen in the core for remnants
    }

    GNU_CONST double CalculateHAbundanceSurface(const double p_Tau) const override { // JR FIX THIS: DONE
        return 0.0; // No hydrogen on the surface for remnants
    }


    GNU_CONST double CalculateHeAbundanceCore(const double p_Tau) const override { // JR FIX THIS: DONE
        return 0.0; // No helium in the core for remnants (except HeWD)
    }

    GNU_CONST double CalculateHeAbundanceSurface(const double p_Tau) const override { // JR FIX THIS: DONE
        return 0.0; // No helium on the surface for remnants (except HeWD)
    }
    

    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////

    
    GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const override { // JR FIX THIS: DONE
        return 0.0; // Tau (relative age) is not used for remnants in Hurley et al. 2000, so we return 0.0
    };
    GNU_CONST double CalculateTauAtPhaseEnd_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const override { // JR FIX THIS: DONE
        return 0.0; // Tau (relative age) is not used for remnants in Hurley et al. 2000, so we return 0.0
    };
    

TimescalesT CalculateTimescales_Hurley2000() const override {
    return TPAGB::CalculateTimescales_Hurley2000(Mass0(), m_InterimState.HurleyGBParams(), m_InterimState.HurleyTimescales())); }

inline double CalculateTimescale_Thermal() const override { return CalculateTimescale_Dynamical(); }



    GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const override;



inline double CalculateHurleyPerturbationMu() const { return m_InterimState.Mu(); }



    //////////////////////////////////////////////////
    //   ENVELOPE                                   //
    //////////////////////////////////////////////////

    GNU_CONST double CalculateConvectiveEnvelopeRadialExtent() const override { return 0.0; } // No convective envelope for remnants JR FIX THIS: DONE

GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::REMNANT; } // Always REMNANT for remnants

inline Dbl_DblT CalculateConvectiveEnvelopeMass() const override { return std::make_tuple(0.0, 0.0); } // No envelope for remnants



    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////

inline double CalculateCOCoreMassAtPhaseEnd() const override { return Mass(); } // McCO = M for remnants

inline double CalculateHeCoreMassAtPhaseEnd() const override { return Mass(); } // McHe = M for remnants



GNU_CONST static Dbl_DblT CalculateRemnantMass(const double p_COCoreMass);
GNU_CONST static Dbl_DblT CalculateRemnantMass_Hurley2000(const double p_COCoreMass);


    inline double CalculateCOCoreMass() const override { return Mass(); } // McCO = M for remnants

inline double CalculateHeCoreMass() const override { return Mass(); } // McHe = M for remnants

    double CalculateConvectiveCoreRadius() const { return Radius(); }


    double CalculateCoreMassOnPhase() const { return Mass(); }







    //////////////////////////////////////////////////
    //   MASS LOSS / ACCRETION                      //
    //////////////////////////////////////////////////

inline double CalculateMLRateThermal() const override{ return BaseStar::CalculateMLRateThermal(); } // Set thermal mass gain rate to be effectively infinite, using dynamical timescale (in practice, will be Eddington limited), avoid division by zero





    //////////////////////////////////////////////////
    //   MISC. CONSTANTS / PARAMETERS               //
    //////////////////////////////////////////////////

    GNU_PURE GBParamsT CalculateGBParams_Hurley2000(const double p_Mass, const GBParamsT& p_GBParams) const override { // JR FIX THIS: DONE
        // Not meaningful for BH, so we just return current GB parameters
        return GBParams();
    }





    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////

inline double CalculateRemnantRadius() const override { return Radius(); }

    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                              { return BaseStar::CalculateRadiusAndStellarTypeOnPhase(); }
   




    //////////////////////////////////////////////////
    //   SUPERNOVAE                                 //
    //////////////////////////////////////////////////


    double          CalculateInitialSupernovaMass() const                                                       { return GiantBranch::CalculateInitialSupernovaMass(); }                // Use GiantBranch


    















    STELLAR_TYPE    EvolveToNextPhase()                                                                         { return BaseStar::EvolveToNextPhase(); }                               // Default to BaseStar

    bool            IsDegenerate() const                                                                        { return true; }                                                        // White Dwarfs, NS and BH are degenerate

    bool            IsEndOfPhase() const                                                                        { return !ShouldEvolveOnPhase(); }                                      // Phase ends when envelope loss or going supernova

    bool            IsSupernova() const                                                                         { return false; }                                                       // Default

    void            PerturbLuminosityAndRadius() { }                                                                                                                                    // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false)                                                   { return BaseStar::ResolveEnvelopeLoss(p_Force); }                      // Default to BaseStar

    void            ResolveEnvelopeMassAtPhaseEnd(const double p_Tau) const                                     { ResolveEnvelopeMassOnPhase(p_Tau); }                                  // Same as on phase
    void            ResolveEnvelopeMassOnPhase(const double p_Tau) const { }                                                                                                            // NO-OP

    STELLAR_TYPE    ResolveMassLoss(const double p_dt)                                                          { return StellarType(); }                                               // NO-OP
    STELLAR_TYPE    ResolveMassLossHurley(const double p_dt)                                                    { return ResolveMassLoss(p_dt); }                                                                                                                                   // NO-OP

    STELLAR_TYPE    ResolveSkippedPhase()                                                                       { return BaseStar::ResolveSkippedPhase(); }                             // Default to BaseStar
                                                                                                                                                                                        //
    void            ResolveShellChange(const double p_AccretedMass) { }                                                                                                                 // NO-OP 
                                                                                                                                                                                        //
    STELLAR_TYPE    ResolveSupernova()                                                                          { return BaseStar::ResolveSupernova(); }                                // Default to BaseStar

    void            SetPulsarParameters() const { }                                                                                                                                     // NO-OP

    bool            ShouldEnvelopeBeExpelledByPulsations() const                                                { return false; }                                                       // No envelope to lose by pulsations
    bool            ShouldEvolveOnPhase() const                                                                 { return true; }                                                        // Default
    bool            ShouldSkipPhase() const                                                                     { return false; }                                                       // Don't skip WD phase

    Dbl_DblT CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate);

};



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

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
 * @param       p_tScales                       Timescales (Myr) (not used here)
 * @return                                      Timestep (Myr)
 */
inline double Remnants::ChooseTimestep_Hurley2000(const double p_Age, [[maybe_unused]] const TimescalesT& p_tScales) const {

    double dtk = std::min(std::max(1.0, p_Age), 500.0); // Modfied from Hurley
                                                                                    
    // Time to end of phase (change of stellar type, dte) not used here - how to calculate?  // JR FIX THIS ??? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // Clamp to minimum NUCLEAR_MINIMUM_TIMESTEP
    return std::max(dtk, NUCLEAR_MINIMUM_TIMESTEP);
}


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

/*
 * CalculateRemnantMass
 *
 * @brief
 * Calculate remnant mass.
 * 
 * Calls relevant remnant mass function based on the evolutionary mode given in program options.
 * 
 * 
 * static Dbl_DblT CalculateRemnantMass(const double p_COCoreMass)
 *
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
inline Dbl_DblT Remnants::CalculateRemnantMass(const double p_COCoreMass) {

    double mass = 0.0;                                              // Default return value                               
    double fallbackFraction = 0.0;                                  // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            std::tie(mass, fallbackFraction) = CalculateRemnantMass_Hurley2000(p_COCoreMass);
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNEXPECTED_EVOLUTION_MODE);   // Throw error
    }       

    return std::make_tuple(mass, fallbackFraction);
}

/*
 * CalculateRemnantMass_Hurley2000
 *
 * @brief
 * Calculate remnant mass, per Hurley et al. 2000, eq 92.
 * 
 * Hurley et al. 2000 indicates a star with a CO core > 7.0 Msol leads to a collapse into a BH,
 * but is ambiguous about the BH mass in this case - we assume here that it's just the CO core.
 *
 * 
 * static Dbl_DblT CalculateRemnantMass_Hurley2000(const double p_COCoreMass)
 *
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
inline Dbl_DblT Remnants::CalculateRemnantMass_Hurley2000(const double p_COCoreMass) {
    return std::make_tuple((p_COCoreMass > 7.0 ? p_COCoreMass : 1.17 + (0.09 * p_COCoreMass)), 0.0); // Fallback fraction not defined by Hurley
}











/// Remnants_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//



    // Zeta
    double  CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription) override { (void)p_ZetaPrescription; return 0.0; }                                            // Should never be called
    double  CalculateZetaEquilibrium() override { return -std::numeric_limits<double>::infinity(); }                                                                                 // Nuclear timescale MT should be impossible from remnant stars; should never be called

    // CE lambda variants
    double  CalculateCELambda_Dewi() const override { return 1.0; }                                                                                                                  // Was: BaseStar:: default
    COMPAS_PURE double  CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const override { (void)p_Mass; (void)p_Radius; (void)p_CoreMass; return 1.0; }                                                                        // Was: BaseStar:: default with dummy args

    // Mass acceptance
    Dbl_DblT CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                         const double p_AccretorMassRate,
                                         const bool   p_IsHeRich) override;

    // Tidal response
    Dbl_Dbl_Dbl_DblT CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const override   { (void)p_Omega; (void)p_SemiMajorAxis; (void)p_M2; return std::make_tuple(0.0, 0.0, 0.0, 0.0); }              // No tidal response
    Dbl_Dbl_Dbl_DblT CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const override { (void)p_Omega; (void)p_SemiMajorAxis; (void)p_M2; return std::make_tuple(0.0, 0.0, 0.0, 0.0); }              // No tidal response

    // Critical mass ratio
    double  CalculateCriticalMassRatio(const double p_Mass,
                                       const double p_Radius,
                                       const double p_CoreMass,
                                       const bool   p_AccretorIsDegenerate,
                                       const double p_MTefficiency = 0.0) const override {
        (void)p_Mass; (void)p_Radius; (void)p_CoreMass; (void)p_AccretorIsDegenerate; (void)p_MTefficiency;
        return 0.0;     // Should not be called; if it is, MT from a remnant always treated as stable
    }
    double  CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const override { (void)p_AccretorIsDegenerate; return 0.0; }
    double  CalculateCriticalMassRatio_Ge2020(const double p_MTefficiency) override { (void)p_MTefficiency; return 0.0; }
    double  CalculateCriticalMassRatio_Hurley2002() const override { return 0.0; }









virtual double CalculateEddingtonCriticalRate() const { return CalculateEddingtonCriticalRate(Radius()); }
virtual double CalculateEddingtonCriticalRate(const double p_Radius) const { return 2.08E-3 / 1.7 * p_Radius * MYR_TO_YEAR * OPTIONS->EddingtonAccretionFactor() ; } // Hurley+, 2002, Eq. (67)









