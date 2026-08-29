#pragma once

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "HeMS.h"


class BaseStar;
class HeMS;

class HeHG: virtual public BaseStar, public HeMS {




public:

    HeHG() { m_InterimState.SetStellarType(STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP); };
    
    HeHG(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeMS(p_BaseStar, false) {
        m_InterimState.SetStellarType(STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP);                                                                                                          // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                                     // Initialise if required
    }



private:

    void Initialise() {
        m_InterimState.SetTau(0.0);                                                                       // Start of phase
        CalculateTimescales();                                                                            // Initialise timescales
        // Age for HeHG is calculated before switching -
        // can get here via EvolveOneTimestep() and ResolveEnvelopeLoss(),
        // and Age is calculated differently in those cases
        
        // Update stellar properties at start of HeHG phase (since core definition changes)
        CalculateGBParams();

        EvolveOnPhase(0.0);
    }



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

    GNU_CONST double CalculateHAbundanceCore([[maybe_unused]]const double p_Tau) const override { // JR FIX THIS: DONE
        return 0.0; // No hydrogen in the core for HeHG stars
    }
    GNU_CONST double CalculateHAbundanceCoreAtPhaseEnd() const override { // JR FIX THIS: DONE
        return 0.0; // No hydrogen in the core for HeHG stars
    }

    GNU_CONST double CalculateHeAbundanceCore([[maybe_unused]]const double p_Tau) const override { // JR FIX THIS: DONE
        return 0.0; // No helium in the core for HeHG stars
    }
    GNU_CONST double CalculateHeAbundanceCoreAtPhaseEnd() const override { // JR FIX THIS: DONE
        return 0.0; // No helium in the core for HeHG stars
    }


    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateAgeAfterMassLoss_Hurley2000() const override { // JR FIX THIS: DONE
        return GiantBranch::CalculateAgeAfterMassLoss_Hurley2000(); // NO-OP
    }


    GNU_CONST double CalculateTau_Hurley2000() const override { // JR FIX THIS: DONE
        return 0.0; // Tau (relative age) is not used for HeHG stars in Hurley et al. 2000, so we return 0.0
    }
    GNU_CONST double CalculateTauAtPhaseEnd_Hurley2000([[maybe_unused]]const double p_Age, [[maybe_unused]]const TimescalesT& p_tScales) const override { // JR FIX THIS: DONE
        return 0.0; // Tau (relative age) is not used for HeHG stars in Hurley et al. 2000, so we return 0.0
    };


    GNU_CONST TimescalesT CalculateTimescales_Hurley2000(const double p_Mass, const GBParamsT& p_GBParams, const TimescalesT& p_tScales) const; // JR FIX THIS: DONE
    GNU_PURE  TimescalesT CalculateTimescales_Hurley2000() { return CalculateTimescales_Hurley2000(Mass0(), GBParams(), Timescales()); } // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateLuminosityAtPhaseEnd_Hurley2000() const override { // JR FIX THIS: DONE
        return Luminosity(); // NO-OP
    }
    
    GNU_PURE double CalculateLuminosity_Hurley2000() const override { // JR FIX THIS: DONE
        return HeGB::CalculateLuminosity_Hurley2000(COCoreMass(), GBParams()); // Hurley et al. Section 6.1
    }


    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////
   
    GNU_PURE double CalculateCOCoreMassAtPhaseEnd_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateCOCoreMass_Hurley2000(); // Same as on phase
    }

    GNU_PURE double CalculateCOCoreMass_Hurley2000() const override { // JR FIX THIS: DONE
        return HeGB::CalculateCoreMass_Hurley2000(Mass0(), Age(), GBParams(), Timescales(TS::HeMS)); // Hurley et al. Section 6.1
    }


    GNU_PURE double CalculateCoreMassAtBAGB_Hurley2000([[maybe_unused]]const double p_Mass) const override { // JR FIX THIS: DONE
        return Mass0(); // McBAGB = M0 (Hurley et al. 2000, discussion just before eq 89) 
    }


    GNU_PURE double CalculateCoreMassAtPhaseEnd_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateCoreMass_Hurley2000(); // Same as on phase
    }

    GNU_PURE double CalculateCoreMass_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateCOCoreMass_Hurley2000(); // Mc(HeHG) = MCOc
    }


    GNU_PURE inline double CalculateHeCoreMassAtPhaseEnd_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateHeCoreMass_Hurley2000(); // Same as on phase
    }

    GNU_PURE inline double CalculateHeCoreMass_Hurley2000() const override { // JR FIX THIS: DONE
        return Mass(); // McHe = M for HeHG stars
    }


            double          CalculateConvectiveCoreMass() const                                                     { return CalculateConvectiveCoreMass(CoreMass()); }
            double          CalculateConvectiveCoreMass(const double p_CoreMass) const                              { return p_CoreMass; }






    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateRadiusAtPhaseEnd_Hurley2000() const override { // JR FIX THIS: DONE
        return Radius(); // NO-OP
    }
    
    GNU_PURE double CalculateRadius_Hurley2000() const override { // JR FIX THIS: DONE
        return HeGB::CalculateLuminosity_Hurley2000(COCoreMass(), GBParams()); // Hurley et al. Section 6.1
    }




            double          CalculateRadiusAtPhaseEnd() const                                                       { return CalculateRadiusAtPhaseEnd(Radius()); }                                      // NO-OP
            double          CalculateRadiusAtPhaseEnd(const double p_Radius) const                                  { return p_Radius; }
   
            double          CalculateRadiusOnMassChange(double p_dM)                                                { return CalculateRadiusOnPhase(Mass() + p_dM, Luminosity()); }
            double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhase(Mass(), Luminosity()); }
            double          CalculateRadiusOnPhase(double p_Mass, double p_Luminosity) const { } // JR FIX THIS: placeholder


            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity) const { } // JR FIX THIS: placeholder
            std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                          { return CalculateRadiusAndStellarTypeOnPhase(Mass(), Luminosity()); }

GNU_CONST static double CalculateConvectiveCoreRadius_Hurley2000_Static(const double p_Radius, const double p_Tau, const double p_CoreMass);
GNU_CONST static double CalculateConvectiveCoreRadius_Static(const double p_Radius, const double p_Tau, const double p_CoreMass);
GNU_CONST inline double CalculateConvectiveCoreRadius(const double p_Radius, const double p_Tau, const double p_CoreMass) const { return CalculateConvectiveCoreRadius_Static(p_Radius, p_Tau, p_CoreMass); }




inline double CalculateRemnantRadius_Hurley2000() const override {
    return CalculateRemnantRadius_Hurley2000_Static(CoreMass());
}
GNU_CONST static double CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass);














            double          CalculateMomentOfInertia() const                                                        { return GiantBranch::CalculateMomentOfInertia(Mass(), Radius(), CoreMass()); }


GNU_CONST double CalculateSmallEnvPerturbationParm_Hurley2000(const double p_Mass, const double p_CoreMass, [[maybe_unused]] const double p_Luminosity) const override;

inline double CalculateHurleyPerturbationMuAtPhaseEnd() const                                                       { return CalculateHurleyPerturbationMuAtPhaseEnd(CurrentState().Mu()); }
inline double CalculateHurleyPerturbationMuAtPhaseEnd(const double p_Mu) const                                      { return p_Mu; }








///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<






    double CalculateLuminosity_Hurley2000() const override;

    double CalculateRadius_Hurley2000() const override;

 



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return Temperature(); }                                                           // NO-OP (deliberately ignores inputs; returns current temperature)
            double          CalculateTemperatureAtPhaseEnd() const                                                  { return CalculateTemperatureAtPhaseEnd(Luminosity(), Radius()); }                   // Use state-routed values

inline double CalculateMLRateThermal() const override { return GiantBranch::CalculateMLRateThermal(); } // Skip HeMS


GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const override;


COMPAS_PURE ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const override;


            STELLAR_TYPE    EvolveToNextPhase();

            bool            IsEndOfPhase() const                                                                    { return !ShouldEvolveOnPhase(); }
            bool            IsSupernova() const;
            double          CalculateInitialSupernovaMass() const;
    
            void            PerturbLuminosityAndRadius()                                                            { GiantBranch::PerturbLuminosityAndRadius(); }                          // NO-OP

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
            void            ResolveHeliumFlash() { }                                                                                                                                        // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                                   { return ResolveSkippedPhase(StellarType()); }                                       // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase(const STELLAR_TYPE p_StellarType)                                   { return p_StellarType; }

            bool            ShouldEnvelopeBeExpelledByPulsations() const                                            { return CHeB::ShouldEnvelopeBeExpelledByPulsations(); }                // Envelope of convective star with luminosity to mass ratio beyond threshold should be expelled
            bool            ShouldEvolveOnPhase() const;
            bool            ShouldSkipPhase() const                                                                 { return false; }                                                       // Never skip HeMS phase

            void            UpdateAgeAfterMassLoss()                                                                { GiantBranch::UpdateAgeAfterMassLoss(); }                              // No action for He giants
            double          CalculateInitialMass()                                                                  { return MassEffectiveInitial(); } // JR FIX THIS: check <<<<<<<<<<<<<<<<<<<<<<<<<
};




/////// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




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



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

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
 * TimescalesT CalculateTimescales_Hurley2000(const double p_Mass, const GBParamsT& p_GBParams, const TimescalesT& p_tScales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBParams                      GB parameters
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
TimescalesT HeHG::CalculateTimescales_Hurley2000(const double p_Mass, const GBParamsT& p_GBParams, const TimescalesT& p_tScales) const {

    const double AHe = p_GBParams[HURLEY_GBP::AHe];
    const double B   = p_GBParams[HURLEY_GBP::B];
    const double D   = p_GBParams[HURLEY_GBP::D];
    const double p   = p_GBParams[HURLEY_GBP::P];
    const double q   = p_GBParams[HURLEY_GBP::Q];
    const double Lx  = p_GBParams[HURLEY_GBP::LX];
    const double Mx  = p_GBParams[HURLEY_GBP::MX];

    const double p1    = p - 1.0;
    const double p1_p  = p1 / p;
    const double q1    = q - 1.0;

    const double lTHe = HeMS::CalculateLuminosityAtPhaseEnd(p_Mass);

    TimescalesT tScales = p_tScales; // Copy given timescales

    tScales[TS::HeMS]      = HeMS::CalculatePhaseLifetime_Hurley2000(p_Mass);
    tScales[TS::Inf1_HeGB] = tScales[TS::HeMS] + (1.0 / ((p1 * AHe * D)) * PPOW((D / lTHe), p1_p));
    tScales[TS::Mx_HeGB]   = tScales[TS::Inf1_HeGB] - (tScales[TS::Inf1_HeGB] - tScales[TS::HeMS]) * PPOW((lTHe / Lx), p1_p);
    tScales[TS::Inf2_HeGB] = tScales[TS::Mx_HeGB] + ((1.0 / (q1 * AHe * B)) * PPOW((B / Lx), q1 / q));

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
GNU_CONST inline double HeHG::ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const {
#define tScales(x) p_tScales[static_cast<int>(HURLEY_TS::x)]

    const double dtk = 0.02 * ((p_Age <= tScales(Mx_HeGB) ? tScales(Inf1_HeGB) : tScales(Inf2_HeGB)) - p_Age);    // stellar type specific dt
                                                                                    
    // time to end of phase (change of stellar type, dte) not used here - how to calculate?
    // clamp to minimum NUCLEAR_MINIMUM_TIMESTEP
    return std::max(dtk, NUCLEAR_MINIMUM_TIMESTEP);

#undef tScales
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




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           MISC. CONSTANTS / PARAMETERS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateSmallEnvPerturbationParm_Hurley2000
 *
 * @brief
 * Calculate the small envelope perturbation parameter, per Hurley et al. 2000, eq 98 (&89)
 *
 *
 * double CalculateSmallEnvPerturbationParm_Hurley2000(const double p_Mass, const double p_CoreMass, [[maybe_unused]] const double p_Luminosity) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol) (not used here)
 * @return                                      Perturbation parameter
 */
inline double HeHG::CalculateSmallEnvPerturbationParm_Hurley2000(const double p_Mass, const double p_CoreMass, [[maybe_unused]] const double p_Luminosity) const {
    // Hurley et al. 2000, eq 89 - clamped to ensure McMax >= p_CoreMass
    const double McMax = std::max(std::min(((1.45 * p_Mass) - 0.31), p_Mass), p_CoreMass);
    return 5.0 * (McMax - p_CoreMass) / McMax; // Hurley et al. 2000, eq 98
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
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
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Convective core radius (Rsol)         
 */
GNU_CONST inline double HeHG::CalculateConvectiveCoreRadius_Hurley2000_Static(const double p_Radius, const double p_Tau, const double p_CoreMass) {
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
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Convective core radius (Rsol) 
 */
GNU_CONST inline double HeHG::CalculateConvectiveCoreRadius_Static(const double p_Radius, const double p_Tau, const double p_CoreMass) {

    double radius;

    switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            radius = CalculateConvectiveCoreRadius_Hurley2000_Static(p_Radius, p_Tau, p_CoreMass);
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR_STATIC(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return radius;
}


/*
 * CalculateRemnantRadius_Hurley2000_Static
 *
 * @brief
 * Calculate the radius of the remnant the star would become if it lost all of its envelope
 * immediately (i.e. M = Mc), per Hurley et al. 2000, at the end of section 6 (after eq 105).
 *
 *
 * static double CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass)
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Remnant core radius (Rsol)
 */
// JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
double HeHG::CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass) const {
    return HeWD::CalculateRadiusOnPhase_Static(p_CoreMass);
}









/// HeHG_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//

    // MT rejuvenation
    double  CalculateMTRejuvenationFactor() override { return 1.0; }

    // CE lambda Nanjing variants
    COMPAS_PURE double  CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const override;                                    // body in HeHG.cpp

    // Critical mass ratio
    double  CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const override;                // body inline below
    double  CalculateCriticalMassRatio_Hurley2002() const override { return HURLEY_HJELLMING_WEBBINK_QCRIT_HE_GIANT; }

    // Zeta (HeHG treated as giant-branch-like, not MS-like)
    double  CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription) override;

    double  CalculateZetaEquilibrium() override { return -std::numeric_limits<double>::infinity(); }                // Nuclear timescale MT impossible from HG-stage He stars



/*
 * CalculateCriticalMassRatio_Claeys2014
 *
 * @brief
 * Calculate the critical mass ratio, per Claeys et al. 2014
 * 
 * @param       p_AccretorIsDegenerate          Boolean indicating if accretor is degenerate
 * @return                                      Critical mass ratio for unstable MT 
 */
inline double HeHG_Constituent::CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {                                                                                                          
    return p_AccretorIsDegenerate
            ? OPTIONS->MassTransferCriticalMassRatioHeliumHGDegenerateAccretor()        // degenerate accretor
            : OPTIONS->MassTransferCriticalMassRatioHeliumHGNonDegenerateAccretor();    // non-degenerate accretor
}
