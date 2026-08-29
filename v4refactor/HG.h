#pragma once

#include "constants.h"
#include "GiantBranch.h"
#include "ConstituentStar.h"   // for HG_Constituent (declared at the bottom of this file)

#include <boost/math/tools/roots.hpp>


class BaseStar;
class GiantBranch;

class HG: virtual public BaseStar, public GiantBranch {


public:
    
    HG() {};
    HG(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), GiantBranch(p_BaseStar) { if (p_Initialise) Initialise(); }
    
    
    
private:

    void Initialise() {
        
        m_InterimState.SetTau(0.0);                                                                                                                                                         // Start of phase
        
        // update stellar properties at start of HG phase (since core definition changes)
        CalculateGBParams();
        CalculateTimescales();
        // Initialise timescales
        m_InterimState.SetAge(m_InterimState.HurleyTimescales(TS::MS)); // Set age appropriately
        
        // update effective "initial" mass (m_Mass0) so that core mass matches main sequence core mass
        // (only relevant if MANDEL or BRCEK main sequence core mass prescription is used)
        if (utils::Compare(CalculateCoreMassOnPhase(Mass0(), Age()), std::min(Mass(), MainSequenceCoreMass())) < 0 ||
            (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK && utils::Compare(MZAMS(), BRCEK_LOWER_MASS_LIMIT) >= 0)) {
            double desiredCoreMass = std::min(Mass(), MainSequenceCoreMass());                                                                                                              // desired core mass
            double newMass0        = std::max(Mass0ToMatchDesiredCoreMass(this, desiredCoreMass), std::min(Mass(), MZAMS()));                                                                // use root finder to find new core mass estimate, newMass0 should not be lower than Mass() unless star gained mass
            if (newMass0 <= 0.0) {                                                                                                                                                          // no root found - no solution for estimated core mass
                newMass0 = Mass();                                                                                                                                                          // if no root found we keep newMass0 equal to the total mass
            }
            m_InterimState.SetMassEffectiveInitial(newMass0);                                                                                                                               // commit
            CalculateGBParams();
            CalculateTimescales();
            m_InterimState.SetAge(m_InterimState.HurleyTimescales(TS::MS));
        }
        EvolveOnPhase(0.0);  // <<<<< forces calculation of attributes - we need radius for BRCEK prescription when we clone an HG star...
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

    GNU_CONST inline double CalculateHAbundanceCore() const override { return 0.0; } // No hydrogen in the core for HG stars JR FIX THIS: DONE
    GNU_CONST inline double CalculateHAbundanceCoreAtPhaseEnd() const override { CalculateHAbundanceCore(); } // Same as on phase JR FIX THIS: DONE

    GNU_PURE  inline double CalculateHAbundanceSurface() const override { return HAbundanceSurface(); } // JR FIX THIS: DONE
    GNU_PURE  inline double CalculateHAbundanceSurfaceAtPhaseEnd() const override { return CalculateHAbundanceSurface(); } // Same as on phase JR FIX THIS: DONE
    

    GNU_PURE  inline double CalculateHeAbundanceCore() const override { return 1.0 - Metallicity(); } // JR FIX THIS: DONE
    GNU_PURE  inline double CalculateHeAbundanceCoreAtPhaseEnd() const override { return CalculateHeAbundanceCore(); } // Same as on phase JR FIX THIS: DONE

    GNU_PURE  inline double CalculateHeAbundanceSurface(const double p_Tau) const override { return HeAbundanceSurface(); } // JR FIX THIS: DONE
    GNU_PURE  inline double CalculateHeAbundanceSurfaceAtPhaseEnd(const double p_Tau) const override { return CalculateHeAbundanceSurface(); } // Same as on phase JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateAgeAfterMassLoss_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const TimescalesT& p_tScales) const; // JR FIX THIS: DONE
    double CalculateAgeAfterMassLoss_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateAgeAfterMassLoss_Hurley2000(Metallicity(), Mass(), Age(), Timescales());
    }

    GNU_CONST inline double CalculateTauAtPhaseEnd_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const override { return 1.0; } // 1.0 at phase end JR FIX THIS: DONE
    GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const override; // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   ENVELOPE                                   //
    //////////////////////////////////////////////////

    COMPAS_PURE ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const override;


    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    

inline double CalculateLuminosity_Hurley2000() const override { return CalculateLuminosity_Hurley2000_Static(MassEffectiveInitial(), Tau()); } 
COMPAS_PURE inline double CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Tau) const;





    COMPAS_PURE static double CalculateLuminosityAtPhaseEnd_Hurley2000_Static(const double p_Metallicity, const double p_Mass);

    double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const {  } // JR FIX THIS: placeholder
    double          CalculateLuminosityAtPhaseEnd() const                           { return CalculateLuminosityAtPhaseEnd(Mass0());}                                            // Use state-routed values
                                    // Use class member variables


    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////

    GNU_PURE inline double CalculateCoreMass_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateCoreMass_Hurley2000(Mass(), Tau(), CoreMass(), GBParams());
    }
    GNU_PURE double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau, const double p_CoreMass, const GBParamsT& p_GBParams) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateCoreMass_Hurley2000_Unconstrained(const double p_Mass, const double p_Tau, const DblVectorT& p_GBParams) const; // JR FIX THIS: DONE

    GNU_PURE inline double CalculateCoreMassAtPhaseEnd_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateCoreMassAtPhaseEnd_Hurley2000(Mass(), GBParams());
    }
    GNU_PURE double CalculateCoreMassAtPhaseEnd_Hurley2000(const double p_Mass, const GBParamsT& p_GBParams) const; // JR FIX THIS: DONE


    GNU_PURE inline double CalculateCoreMassAt2ndDredgeUp_Hurley2000() const override { return GBParams(HURLEY_GBP::McDU); } // JR FIX THIS: DONE


    GNU_CONST inline double CalculateCOCoreMass_Hurley2000() const override { return 0.0; } // McCO(HG) = 0.0 JR FIX THIS: DONE
    GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd_Hurley2000() const override { return CalculateCOCoreMass_Hurley2000(); } // Same as on phase JR FIX THIS: DONE


    GNU_PURE inline double CalculateHeCoreMass_Hurley2000() const override { return CoreMass(); } // McHe(HG) = Mc JR FIX THIS: DONE
    GNU_PURE inline double CalculateHeCoreMassAtPhaseEnd_Hurley2000() const override { return CalculateHeCoreMass_Hurley2000(); } // Same as on phase JR FIX THIS: DONE


    GNU_PURE inline double CalculateEffectiveInitialMass_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateEffectiveInitialMass_Hurley2000(Mass(), Mass0(), Tau(), CoreMass(), GBParams());
    }
    GNU_PURE double CalculateEffectiveInitialMass_Hurley2000(const double p_Mass, const double p_MassEffectiveInitial, const double p_Tau, const double p_Coremass, const GBParamsT& p_GBParams) const; // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   MISC. CONSTANTS / PARAMETERS               //
    //////////////////////////////////////////////////

    GNU_CONST double CalculateRho_Hurley2000(const double p_Mass) const; // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////


GNU_CONST double CalculateRadiusOnPhase_Hurley(const double p_Mass, const double p_Tau, const double p_RZAMS, const DblVectorT& p_bN) const {  } // JR FIX THIS: placeholder

    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Tau, const double p_RZAMS) const {  } // JR FIX THIS: placeholder
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return GiantBranch::CalculateRadius_Hurley2000(p_Mass, p_Luminosity); } // Treats HG stars as GB stars
    
double CalculateRadiusOnPhase() const override { return CalculateRadiusOnPhase(Mass0(), Tau(), RadiusZAMSEffective()); }

    



    double          CalculateRadiusAtPhaseEnd(const double p_Mass) const {  } // JR FIX THIS: placeholder
    double          CalculateRadiusAtPhaseEnd() const                               { return CalculateRadiusAtPhaseEnd(Mass()); }      

COMPAS_PURE double CalculateRadiusAtPhaseEnd_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity) const;


















///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<











GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const override;




    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                            { return !ShouldEvolveOnPhase(); }                                                          // Phase ends when age at or after Base Giant Branch MS timescale
    bool            IsSupernova() const                                             { return false; }                                                                           // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash() {  }                                                                                                                                   // NO-OP
    STELLAR_TYPE    ResolveSkippedPhase()                                           { return ResolveSkippedPhase(StellarType()); }                                                                              // NO-OP
    STELLAR_TYPE    ResolveSkippedPhase(const STELLAR_TYPE p_StellarType)            { return p_StellarType; }

    bool            ShouldEvolveOnPhase() const                                     { return (utils::Compare(InterimState().Age(), InterimState().HurleyTimescales(TS::BGB)) < 0); }    // Evolve on HG phase if new age < Base Giant Branch timescale (transition phase reads InterimState)
    bool            ShouldSkipPhase() const                                         { return false; }                                                                           // Never skip HG phase

////    double          CalculateTAMSCoreMass() const                                   { return 0.0; }








    COMPAS_PURE double CalculateRadiusOnPhase_Hurley2000(const double p_Mass, const double p_Luminosity, const double p_Tau, const double p_RZAMS, const double p_MHeF, const double p_MFGB, const double p_Alpha1, const DBL_VECTOR& p_bN) const;
    COMPAS_PURE double CalculateCELambda_Nanjing_Enhanced(const double p_Mass, const double p_Radius, const double p_CoreMass, const SizeT p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;





    void            UpdateEffectiveZAMSLandR()                                      { BaseStar::UpdateEffectiveZAMSLandR(); }                                                   // Skip MainSequence
       
    /*
     * Functor for Mass0ToMatchDesiredCoreMass()
     *
     *
     * Constructor: initialise the class
     * template <class T> Mass0YieldsDesiredCoreMassFunctor(HG *p_Star, double p_DesiredCoreMass)
     *
     * @param   [IN]    p_Star                      (Pointer to) The star under examination
     * @param   [IN]    p_DesiredCoreMass           The desired core mass
     * 
     * Function: core mass difference after mass loss
     * T RadiusEqualsRocheLobeFunctor(double const& p_dM)
     * 
     * @param   [IN]    p_GuessMass0                Guess for Mass0
     * @return                                      Difference between estimated core mass (based on p_GuessMass0) and desired core mass (p_DesiredCoreMass)
     */
    template <class T>
    struct Mass0YieldsDesiredCoreMassFunctor {
        Mass0YieldsDesiredCoreMassFunctor(HG *p_Star, double p_DesiredCoreMass) {
            m_Star            = p_Star;
            m_DesiredCoreMass = p_DesiredCoreMass;
        }
        T operator()(double const& p_GuessMass0) {
        
            // We need an estimate of the core mass of the star so we clone the star without
            // initialisation (i.e. we leave it where it is on the phase) so we can calculate
            // and query its core mass.
            //
            // To ensure the clone does not participate in logging, we set its persistence to EPHEMERAL.

            auto clone = m_Star->CloneAs(STELLAR_TYPE::HERTZSPRUNG_GAP, OBJECT_PERSISTENCE::EPHEMERAL, false);
            clone->AdvanceOneTimestep(0.0, 0.0, p_GuessMass0 - clone->Mass0());
            HG* hgClone = dynamic_cast<HG*>(clone.get());
            double coreMassEstimate = hgClone->CalculateCoreMassOnPhase(p_GuessMass0, hgClone->Age());      // calculate clone's core mass

            return (coreMassEstimate - m_DesiredCoreMass);
        }
    private:
        HG *m_Star;
        double m_DesiredCoreMass;
    };
    
    
    /*
     * Root solver to determine "initial" mass (m_Mass0) based on desired core mass
     *
     * Uses boost::math::tools::bracket_and_solve_root()
     *
     *
     * double Mass0ToMatchDesiredCoreMass(HG *p_Star, double p_DesiredCoreMass)
     *
     * @param   [IN]    p_Star                      (Pointer to) The star under examination
     * @param   [IN]    p_DesiredCoreMass           The desired core mass
     * @return                                      Root found: will be -1.0 if no acceptable real root found
     */
    double Mass0ToMatchDesiredCoreMass(HG *p_Star, double p_DesiredCoreMass) {

        const boost::uintmax_t maxit = ADAPTIVE_MASS0_MAX_ITERATIONS;                                       // Limit to maximum iterations.
        boost::uintmax_t it          = maxit;                                                               // Initially our chosen max iterations, but updated with actual.

        // Find root.
        // We use an iterative algorithm to find the root here:
        //    - if the root finder throws an exception, we stop and return a negative value for the root (indicating no root found)
        //    - if the root finder reaches the maximum number of (internal) iterations, we stop and return a negative value for the root (indicating no root found)
        //    - if the root finder returns a solution, we check that func(solution) = 0.0 +/ ROOT_ABS_TOLERANCE
        //       - if the solution is acceptable, we stop and return the solution
        //       - if the solution is not acceptable, we reduce the search step size and try again
        //       - if we reach the maximum number of search step reduction iterations, or the search step factor reduces to 1.0 (so search step size = 0.0),
        //         we stop and return a negative value for the root (indicating no root found)

        double guess      = p_Star->Mass();                                                                 // Rough guess at solution
        
        double factorFrac = ADAPTIVE_MASS0_SEARCH_FACTOR_FRAC;                                              // Search step size factor fractional part
        double factor     = 1.0 + factorFrac;                                                               // Factor to determine search step size (size = guess * factor)

        std::pair<double, double> root(p_DesiredCoreMass, 0.0);                                             // Initialise root - default return
        SizeT tries = 0;                                                                                    // Number of tries
        bool done   = false;                                                                                // Finished (found root or exceed maximum tries)?
        Mass0YieldsDesiredCoreMassFunctor<double> func = Mass0YieldsDesiredCoreMassFunctor<double>(p_Star, p_DesiredCoreMass);
        while (!done) {                                                                                     // While no acceptable root found
        
            bool isRising = func((double)guess) >= func((double)guess * factor) ? false : true;

            // Run the root finder.
            // Regardless of any exceptions or errors, display any problems as a warning, then
            // check if the root returned is within tolerance - so even if the root finder
            // bumped up against the maximum iterations, or couldn't bracket the root, use
            // whatever value it ended with and check if it's good enough for us - not finding
            // an acceptable root should be the exception rather than the rule, so this strategy
            // shouldn't cause undue performance issues.
            try {
                root = boost::math::tools::bracket_and_solve_root(func, guess, factor, isRising, utils::BracketTolerance, it); // Find root
                // Root finder returned without raising an exception
                if (it >= maxit) { SHOW_WARN(ERROR::TOO_MANY_MASS0_ITERATIONS); }                           // Too many root finder iterations
            }
            catch(std::exception& e) {                                                                      // Catch generic boost root finding error
                // Root finder exception.
                // Could be too many iterations, or unable to bracket root - it may not
                // be a hard error - so no matter what the reason is that we are here,
                // we'll just emit a warning and keep trying
                if (it >= maxit) { SHOW_WARN(ERROR::TOO_MANY_MASS0_ITERATIONS); }                           // Too many root finder iterations
                else             { SHOW_WARN(ERROR::ROOT_FINDER_FAILED, e.what()); }                        // Some other problem - show it as a warning
            }

            // wW have a solution from the root finder - it may not be an acceptable solution
            // so we check if it is within our preferred tolerance
            if (std::fabs(func(root.first + (root.second - root.first) / 2.0)) <= ROOT_ABS_TOLERANCE) {     // Solution within tolerance?
                done = true;                                                                                // Yes - we're done
            }
            else if (std::fabs(func(root.first)) <= ROOT_ABS_TOLERANCE) {                                   // Solution within tolerance at endpoint 1?
                root.second=root.first;
                done = true;                                                                                // Yes - we're done
            }
            else if (std::fabs(func(root.second)) <= ROOT_ABS_TOLERANCE) {                                  // Solution within tolerance at endpoint 2?
                root.first=root.second;
                done = true;                                                                                // Yes - we're done
            }
            else {                                                                                          // No - try again
                // We don't have an acceptable solution - reduce search step size and try again
                factorFrac /= 2.0;                                                                          // Reduce fractional part of factor
                factor      = 1.0 + factorFrac;                                                             // New search step size
                tries++;                                                                                    // Increment number of tries
                if (tries > ADAPTIVE_MASS0_MAX_TRIES || std::fabs(factor - 1.0) <= ROOT_ABS_TOLERANCE) {    // Too many tries, or step size 0.0?
                    // we've tried as much as we can - fail here with -ve return value
                    root.first  = -1.0;                                                                     // Yes - set error return
                    root.second = -1.0;
                    SHOW_WARN(ERROR::TOO_MANY_MASS0_TRIES);                                                 // Show warning
                    done = true;                                                                            // We're done
                }
            }
        }
        
        // Midway between brackets is our result.
        // If necessary we could return the result as an interval here.
        return root.first + (root.second - root.first) / 2.0;
    }
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
 * CalculateAgeAfterMassLoss_Hurley2000
 *
 * @brief
 * Recalculate the star's age after mass loss, per Hurley et al. 2000, section 7.1
 *
 *
 * double CalculateAgeAfterMassLoss_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const TimescalesT& p_tScales) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Age of the star after mass loss (Myr)
 */
inline double HG::CalculateAgeAfterMassLoss_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const TimescalesT& p_tScales) const {
    const double tMS       = p_tScales[TS::MS];
    const double tBGB      = p_tScales[TS::BGB];
    const double tBGBprime = astro::CalculateLifetimeToBGB_Hurley2000(p_Metallicity, p_Mass);
    const double tMSprime  = MainSequence::CalculatePhaseLifetime_Hurley2000(p_Metallicity, p_Mass, tBGB);
    return tMSprime + (((tBGBprime - tMSprime) / (tBGB - tMS)) * (p_Age - tMS));
}


/*
 * CalculateTau_Hurley2000
 *
 * @brief
 * Calculate the HG-relative age (fractional Hertzsprung Gap age) of the star,
 * per Hurley et al. 2000, eq 25
 *
 * 
 * double CalculateTau_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      HG-relative age of the star, [0, 1]
 */
inline double HG::CalculateTau_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const {
    const double tMS = p_tScales[TS::MS]
    return std::clamp((p_Age - tMS) / (p_tScales[TS::BGB] - tMS), 0.0, 1.0);
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
inline double HG::ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const {
    const double tBGB = p_tScales[TS::BGB]
    const double dtk  = 0.05 * (tBGB - p_tScales[TS::MS]);          // Stellar type specific dt
    const double dte  = tBGB - p_Age;                               // Time to end of phase (change of stellar type)
    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);  // Clamp to minimum NUCLEAR_MINIMUM_TIMESTEP
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ENVELOPE                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * DetermineEnvelopeType
 *
 * @brief
 * Determine the star's envelope type, based on the user-specified ENVELOPE_STATE_PRESCRIPTION.
 * 
 *
 * ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Envelope type (ENVELOPE::{RADIATIVE, CONVECTIVE, REMNANT})
 */
ENVELOPE HG::DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const {
 
    ENVELOPE envType;
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                             // Which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::CONVECTIVE_MASS_FRACTION:             // CONVECTIVE_MASS_FRACTION
            // Envelope is convective when the convective mass exceeds specified fraction of the envelope mass
            double mEnv;
            std::tie(mEnv, std::ignore) = CalculateConvectiveEnvelopeMass();
            envType = mEnv / (p_Mass - p_CoreMass) > OPTIONS->ConvectiveEnvelopeMassThreshold() ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:                    // FIXED_TEMPERATURE
            // Envelope is radiative if temperature exceeds specified threshold, otherwise convective
            envType = (p_Temperature * TSOL) > OPTIONS->ConvectiveEnvelopeTemperatureThreshold() ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:                               // HURLEY
            // Hurley et al. 2002 eqs 39 & 40, and Hurley et al. 2000 end of section 7.2 describe gradual
            // growth of convective envelope over HG, but we approximate it as already convective here
            envType = ENVELOPE::CONVECTIVE;                                     // Always convective
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:                               // COMPAS LEGACY
            envType = ENVELOPE::RADIATIVE;                                      // Always radiative
            break;

        default:                                                                // Unexpected prescription
            // The only way this can happen is if the ENVELOPE_STATE_PRESCRIPTION served to this
            // function is not accounted for in this code.  We should not default here, with or
            // without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR_STATIC(ERROR::UNEXPECTED_ENVELOPE_STATE_PRESCRIPTION);  // Throw error               
    }
    
    return envType;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLuminosityAtPhaseEnd_Hurley2000_Static
 *
 * @brief
 * Calculate luminosity at the end of the Hertzsprung Gap, TAHG,
 * per Hurley et al. 2000, just before eq 8
 *
 *
 * static double CalculateLuminosityAtPhaseEnd_Hurley2000_Static(const double p_Metallicity, const double p_Mass)
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAHG luminosity (Lsol)
 */
inline double HG::CalculateLuminosityAtPhaseEnd_Hurley2000_Static(const double p_Metallicity, const double p_Mass) {
    return (p_Mass < ZDEP->HurleyMCOFGB(p_Metallicity))
            ? GiantBranch::CalculateLuminosityAtBGB_Hurley2000_Static(p_Metallicity, p_Mass)
            : GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(p_Metallicity, p_Mass);
}


/*
 * CalculateLuminosity_Hurley2000_Static
 *
 * @brief
 * Calculate luminosity on the Hertzsprung Gap, per Hurley et al. 2000, eq 26
 *
 *
 * static double CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Tau)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      HG luminosity (Lsol)
 */
COMPAS_PURE inline double HG::CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Tau) const {
    const double lTMS = MainSequence::CalculateLuminosityAtPhaseEnd_Hurley(p_Mass, ZDEP->HurleyACoefficients(Metallicity()));
    return lTMS * PPOW((CalculateLuminosityAtPhaseEnd_Hurley2000_Static(p_Mass, ZDEP->HurleyACoefficients(Metallicity()), ZDEP->HurleyMassCutoffs(Metallicity(), HURLEY_MCO::FGB), ZDEP->HurleyBCoefficients(Metallicity()), ZDEP->HurleyMassCutoffs(Metallicity(), HURLEY_MCO::HeF), ZDEP->HurleyAlpha1(Metallicity())) / lTMS), p_Tau);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateCoreMassAtPhaseEnd_Hurley2000
 *
 * @brief
 * Calculate core mass at the end of the Hertzsprung Gap, per Hurley et al. 2000, eq 28.
 *
 *
 * double CalculateCoreMassAtPhaseEnd_Hurley2000(const double p_Metallicity, const double p_Mass, const GBParamsT& p_GBParams) const {
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBParams                      Hurley GB parameters
 * @return                                      Core mass of the star at TAHG (Base of the Giant Branch) (Msol)
 */
inline double HG::CalculateCoreMassAtPhaseEnd_Hurley2000(const double p_Metallicity, const double p_Mass, const GBParamsT& p_GBParams) const {

    double McTAHG = 0.0; // Default return value

    if (p_Mass < ZDEP->HurleyMCOHeF(p_Metallicity)) {
        const double lum = GiantBranch::CalculateLuminosityAtBGB_Hurley2000(p_Metallicity, p_Mass);
        McTAHG = BaseStar::CalculateCoreMass_Hurley2000(lum, p_GBParams);
    }
    else if (p_Mass < ZDEP->HurleyMCOFGB(p_Metallicity)) {
        McTAHG = p_GBParams[HURLEY_GBP::MCBGB];
    }
    else {
        McTAHG = CalculateCoreMassAtHeI_Hurley2000(p_Metallicity, p_Mass, p_GBParams);
    }

    return McTAHG;
}


/*
 * CalculateCoreMass_Hurley2000
 *
 * @brief
 * Calculate core mass on the Hertzsprung Gap, per Hurley et al. 2000, eq 30
 * (see Section 7
 *
 * The core mass returned is clamped to a minimum value of the current value of the
 * core mass of the star (passed as a parameter) - if the star is losing mass we don't
 * reduce the size of the core (per Hurley).
 *
 *
 * double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau, const double p_CoreMass, const GBParamsT& p_GBParams) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_CoreMass                      Current core mass of the star (Msol)
 * @param       p_GBParams                      Hurley GB parameters
 * @return                                      HG core mass of the star (Msol)
 */
inline double HG::CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau, const double p_CoreMass, const GBParamsT& p_GBParams) const {
    return std::max(CalculateCoreMass_Hurley2000_Unconstrained(p_Mass, p_Tau, p_GBParams), p_CoreMass);
}


/*
 * CalculateCoreMass_Hurley2000_Unconstrained
 *
 * @brief
 * Calculate core mass on the Hertzsprung Gap without accounting for previous core mass,
 * based on Hurley et al. 2000 (eq 30).
 *
 * This ignores the previous core mass constraint (see section 7 of Hurley et al. 2000)
 * when computing the expected core mass, and just follows eq. 30.  This is useful for
 * asking what the core mass would be for the given mass without considering that the
 * core mass should not be allowed to drop.
 *
 *
 * double CalculateCoreMass_Hurley2000_Unconstrained(const double p_Mass, const double p_Tau, const GBParamsT& p_GBParams) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_GBParams                      Hurley GB parameters
 * @return                                      (Unconstrained) HG core mass of the star (Msol)
 */
inline double HG::CalculateCoreMass_Hurley2000_Unconstrained(const double p_Mass, const double p_Tau, const GBParamsT& p_GBParams) const {
    return (((1.0 - p_Tau) * CalculateRho_Hurley2000(p_Mass)) + p_Tau) * CalculateCoreMassAtPhaseEnd_Hurley2000(p_Mass, p_GBParams);
}


/*
 * CalculateEffectiveInitialMass_Hurley2000
 *
 * @brief
 * Calculate effective initial mass of an HG star, per Hurley et al. 2000, section 7.1.
 * 
 * On the HG, the effective initial mass on the HG tracks the stellar mass, unless it would
 * yield an unphysical decrease in the core mass.
 * 
 * Only update effective initial mass on mass loss if the current mass would yield a core
 * mass larger than or equal to the current core mass (i.e. no unphysical core mass decrease
 * would ensue).
 * 
 * Do not update effective initial mass on mass gain on the HG - there is no instruction for
 * doing so in Hurley.  This check also avoids difficulties for the BRCEK rejuvenation
 * prescription, when effective initial mass may be set to enforce a core mass that is lower
 * than would be expected for the current mass value according to the Hurley prescription).
 *
 * 
 * double CalculateEffectiveInitialMass_Hurley2000(const double p_Mass, const double p_MassEffectiveInitial, const double p_Tau, const double p_Coremass, const GBParamsT& p_GBParams) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_MassEffectiveInitial          Current effective initial mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_GBParams                      Hurley GB parameters
 * @return                                      HG effective initial mass (Msol) (may be unchanged)
 */
inline double HG::CalculateEffectiveInitialMass_Hurley2000(const double p_Mass, const double p_MassEffectiveInitial, const double p_Tau, const double p_Coremass, const GBParamsT& p_GBParams) const {
    return ((p_MassEffectiveInitial > p_Mass) && (p_Coremass <= CalculateCoreMass_Hurley2000_Unconstrained(p_Mass, p_Tau, p_GBParams))) ? p_Mass : p_MassEffectiveInitial;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           MISC. CONSTANTS / PARAMETERS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateRho_Hurley2000
 *
 * @brief
 * Calculate the parameter rho for the Hertzsprung Gap, per Hurley et al. 2000, eq 29.
 *
 * Rho is calculated such that McTMS = Rho * McEHG
 * 
 *
 * double CalculateRho_Hurley2000(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      HG rho
 */
inline double HG::CalculateRho_Hurley2000(const double p_Mass) const {
    const double m5_25 = utils::intPow(p_Mass, 5) * std::sqrt(std::sqrt(p_Mass)); // sqrt() is much faster than pow()
    return (1.586 + m5_25) / (2.434 + (1.02 * m5_25));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRadiusAtPhaseEnd_Hurley2000
 *
 * @brief
 * Calculate radius at the end of the Hertzsprung Gap, TAHG,
 * per Hurley et al. 2000, eqs 7 & 8
 *
 *
 * double CalculateRadiusAtPhaseEnd_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass, const double p_MinLuminosity) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MinLuminosity                 Minimum luminosity on phase (Lsol)
 * @return                                      TAHG radius (Rsol)
 */
inline double HG::CalculateRadiusAtPhaseEnd_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass, const double p_MinLuminosity) const {
    return (p_Mass < ZDEP->HurleyMCOFGB(p_Metallicity))
            ? GiantBranch::CalculateRadius_Hurley2000(p_Mass, GiantBranch::CalculateLuminosityAtBGB_Hurley2000_Static(p_Metallicity, p_Mass));
            : CalculateRadiusAtHeI_Hurley2000(p_Mass, p_CoreMass, p_MinLuminosity);
}



/// HG_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//


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
inline double HG_Constituent::CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {                                                                                                          
    return p_AccretorIsDegenerate
            ? OPTIONS->MassTransferCriticalMassRatioHGDegenerateAccretor()      // degenerate accretor
            : OPTIONS->MassTransferCriticalMassRatioHGNonDegenerateAccretor();  // non-degenerate accretor
}


