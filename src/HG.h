#ifndef __HG_h__
#define __HG_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include <boost/math/tools/roots.hpp>

#include "GiantBranch.h"


class BaseStar;
class GiantBranch;

class HG: virtual public BaseStar, public GiantBranch {
    
public:
    
    HG() { m_StellarType = STELLAR_TYPE::HERTZSPRUNG_GAP; };
    
    HG(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), GiantBranch(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::HERTZSPRUNG_GAP;                                                                                                                          // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                         // Initialise if required
    }
    
    HG* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HG* clone = new HG(*this, p_Initialise);
        clone->SetPersistence(p_Persistence);
        return clone;
    }
    
    static HG* Clone(HG& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HG* clone = new HG(p_Star, p_Initialise);
        clone->SetPersistence(p_Persistence);
        return clone;
    }
    
    
private:

COMPAS_PURE double CalculateCoreMassAtPhaseEnd_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const;

GNU_CONST double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau, const double p_CoreMass) const;
static GNU_CONST double CalculateCoreMass_Hurley2000_Unconstrained_Static(const double p_Mass, const double p_Tau) const;
GNU_CONST double CalculateEffectiveInitialMass_Hurley2000(const double p_Mass, const double p_Tau, const double p_Coremass, const double p_MassEffectiveInitial) const;
static GNU_CONST double CalculateRho_Hurley2000_Static(const double p_Mass) const;
GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const double p_tMS, const double p_tBGB) const;



protected:
    
    void Initialise() {
        
        m_Tau = 0.0;                                                                                                                                                            // Start of phase
        
        // update stellar properties at start of HG phase (since core definition changes)
        CalculateGBparams();
        CalculateTimescales();
        // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tMS)];                                                                                                                 // Set age appropriately
        
        // update effective "initial" mass (m_Mass0) so that core mass matches main sequence core mass
        // (only relevant if MANDEL or BRCEK main sequence core mass prescription is used)
        if (utils::Compare(CalculateCoreMassOnPhase(m_Mass0, m_Age), std::min(m_Mass, MainSequenceCoreMass())) < 0 ||
            (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK && utils::Compare(m_MZAMS, BRCEK_LOWER_MASS_LIMIT) >= 0)) {
            double desiredCoreMass = std::min(m_Mass, MainSequenceCoreMass());                                                                                                  // desired core mass
            m_Mass0                = std::max(Mass0ToMatchDesiredCoreMass(this, desiredCoreMass), std::min(m_Mass, m_MZAMS));                                                   // use root finder to find new core mass estimate, m_Mass0 should not be lower than m_Mass unless star gained mass
            if (m_Mass0 <= 0.0) {                                                                                                                                               // no root found - no solution for estimated core mass
                m_Mass0 = m_Mass;                                                                                                                                               // if no root found we keep m_Mass0 equal to the total mass
            }
            CalculateGBparams();
            CalculateTimescales();
            m_Age = m_Timescales[static_cast<int>(TIMESCALE::tMS)];
        }
        EvolveOnPhase(0.0);  // <<<<< forces calculation of attributes - we need radius for BRCEK prescription when we clone an HG star...
    }
    
    
    // member functions - alphabetically





                                                                   // McHe(HG) = Core Mass
    
    double          CalculateHeliumAbundanceCoreAtPhaseEnd() const                  { return 1.0 - m_Metallicity; }
    
    double          CalculateHeliumAbundanceSurfaceAtPhaseEnd() const               { return CalculateHeliumAbundanceSurfaceOnPhase(); }
    
    double          CalculateHydrogenAbundanceCoreAtPhaseEnd() const                { return CalculateHydrogenAbundanceCoreOnPhase(); }
    
    double          CalculateHydrogenAbundanceSurfaceAtPhaseEnd() const             { return CalculateHydrogenAbundanceSurfaceOnPhase(); }

    

inline double CalculateLuminosity_Hurley2000() const override { 
    return CalculateLuminosity_Hurley2000_Static(m_StateHistory.CurrentState.MassEffectiveInitial(), m_StateHistory.CurrentState.Tau());
} 
GNU_CONST inline double CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Tau) const;





    double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtPhaseEnd() const                           { return CalculateLuminosityAtPhaseEnd(m_Mass0);}                                           // Use class member variables
                                    // Use class member variables
    
    // Radius

GNU_CONST double CalculateRadiusOnPhase_Hurley(const double p_Mass, const double p_Tau, const double p_RZAMS, const DBL_VECTOR& p_bN) const;

    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Tau, const double p_RZAMS) const;
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return GiantBranch::CalculateRadiusOnPhase(p_Mass, p_Luminosity); }                                // Treats HG stars as GB stars
    
double CalculateRadiusOnPhase() const override; { return CalculateRadiusOnPhase(m_Mass0, m_Tau, m_RZAMS0); } // Use class member variables

    



    double          CalculateRadiusAtPhaseEnd(const double p_Mass) const;
    double          CalculateRadiusAtPhaseEnd() const                               { return CalculateRadiusAtPhaseEnd(m_Mass); }      









///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


inline double CalculateCoreMass_Hurley2000() const override {
    return CalculateCoreMass_Hurley2000(m_StateHistory.CurrentState.Mass(),
                                        m_StateHistory.CurrentState.Tau(),
                                        m_StateHistory.CurrentState.CoreMass());
}

GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; } // McCO = 0.0 for HG stars

inline double CalculateHeCoreMass() const override { return m_StateHistory.CurrentState.Coremass(); } // McHe = Mc for HG stars

inline double CalculateTau_Hurley2000() const override {
    return CalculateTau_Hurley2000(m_StateHistory.CurrentState.Age(),
                                   m_StateHistory.CurrentState.Timescales(tMS),
                                   m_StateHistory.CurrentState.Timescales(tBGB));
}

inline double CalculateEffectiveInitialMass_Hurley2000() const override {
    return CalculateEffectiveInitialMass_Hurley2000(m_StateHistory.CurrentState.Mass(),
                                                    m_StateHistory.CurrentState.Tau(),
                                                    m_StateHistory.CurrentState.Coremass(),
                                                    m_StateHistory.CurrentState.MassEffectiveInitial());
}





GNU_CONST inline double CalculateHAbundanceCore() const override { return 0.0; } // No hydrogen in the core for HG stars
inline double CalculateHAbundanceSurface() const override { return m_StateHistory.CurrentState.HAbundanceSurface(); }
    
GNU_CONST inline double CalculateHeAbundanceCore() const override { return 1.0 - m_StateHistory.CurrentState.Metallicity(); }
inline double CalculateHeAbundanceSurface() const override { return m_StateHistory.CurrentState.HeAbundanceSurface(); }

 



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


double CalculateCoreMassAtPhaseEnd_Hurley2000() const override {
    return CalculateCoreMassAtPhaseEnd_Hurley2000(m_StateHistory.CurrentState.MassEffectiveInitial(), m_StateHistory.CurrentState.GBparams());
}

GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; } // McCO = 0.0 for HG stars

inline double CalculateHeCoreMassAtPhaseEnd() const override { return m_StateHistory.CurrentState.Coremass(); } // McHe = Mc for HG stars



GNU_CONST inline double CalculateHAbundanceCoreAtPhaseEnd() const override { return CalculateHAbundanceCore(); }
inline double CalculateHAbundanceSurfaceAtPhaseEnd() const override { return CalculateHAbundanceSurface(); }

GNU_CONST inline double CalculateHeAbundanceCoreAtPhaseEnd() const override { return CalculateHeAbundanceCore(); }
inline double CalculateHeAbundanceSurfaceAtPhaseEnd() const override { return CalculateHeAbundanceSurface(); }





    double          ChooseTimestep(const double p_Time) const;

    ENVELOPE        DetermineEnvelopeType() const;

    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                            { return !ShouldEvolveOnPhase(); }                                                          // Phase ends when age at or after Base Giant Branch MS timescale
    bool            IsSupernova() const                                             { return false; }                                                                           // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash() {  }                                                                                                                                   // NO-OP
    STELLAR_TYPE    ResolveSkippedPhase()                                           { return m_StellarType; }                                                                   // NO-OP

    bool            ShouldEvolveOnPhase() const                                     { return (utils::Compare(m_Age, m_Timescales[static_cast<int>(TIMESCALE::tBGB)]) < 0); }    // Evolve on HG phase if age < Base Giant Branch timescale
    bool            ShouldSkipPhase() const                                         { return false; }                                                                           // Never skip HG phase

    double          CalculateTAMSCoreMass() const                                   { return 0.0; }


double CalculateAgeAfterMassLoss() const override;
GNU_PURE  double CalculateAgeAfterMassLoss_Hurley(const double p_Mass, const double p_Age, const DBL_VECTOR& p_tScales, const DBL_VECTOR& p_aN) const;





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

            HG *clone = m_Star->Clone(OBJECT_PERSISTENCE::EPHEMERAL, false);
            clone->EvolveOneTimestep(0.0, p_GuessMass0 - clone->Mass0(), 0.0, true);                        // update clone's mass and age it one timestep 
            double coreMassEstimate = clone->CalculateCoreMassOnPhase(p_GuessMass0, clone->Age());          // calculate clone's core mass
            delete clone; clone = nullptr;                                                                  // return the memory allocated for the clone

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

        const boost::uintmax_t maxit = ADAPTIVE_MASS0_MAX_ITERATIONS;                                       // limit to maximum iterations.
        boost::uintmax_t it          = maxit;                                                               // initially our chosen max iterations, but updated with actual.

        // find root
        // we use an iterative algorithm to find the root here:
        //    - if the root finder throws an exception, we stop and return a negative value for the root (indicating no root found)
        //    - if the root finder reaches the maximum number of (internal) iterations, we stop and return a negative value for the root (indicating no root found)
        //    - if the root finder returns a solution, we check that func(solution) = 0.0 +/ ROOT_ABS_TOLERANCE
        //       - if the solution is acceptable, we stop and return the solution
        //       - if the solution is not acceptable, we reduce the search step size and try again
        //       - if we reach the maximum number of search step reduction iterations, or the search step factor reduces to 1.0 (so search step size = 0.0),
        //         we stop and return a negative value for the root (indicating no root found)

        double guess      = p_Star->Mass();                                                                 // rough guess at solution
        
        double factorFrac = ADAPTIVE_MASS0_SEARCH_FACTOR_FRAC;                                              // search step size factor fractional part
        double factor     = 1.0 + factorFrac;                                                               // factor to determine search step size (size = guess * factor)

        std::pair<double, double> root(p_DesiredCoreMass, 0.0);                                             // initialise root - default return
        std::size_t tries = 0;                                                                              // number of tries
        bool done         = false;                                                                          // finished (found root or exceed maximum tries)?
        Mass0YieldsDesiredCoreMassFunctor<double> func = Mass0YieldsDesiredCoreMassFunctor<double>(p_Star, p_DesiredCoreMass);
        while (!done) {                                                                                     // while no acceptable root found
        
            bool isRising = func((const double)guess) >= func((const double)guess * factor) ? false : true; // gradient direction from guess to upper search increment

            // run the root finder
            // regardless of any exceptions or errors, display any problems as a warning, then
            // check if the root returned is within tolerance - so even if the root finder
            // bumped up against the maximum iterations, or couldn't bracket the root, use
            // whatever value it ended with and check if it's good enough for us - not finding
            // an acceptable root should be the exception rather than the rule, so this strategy
            // shouldn't cause undue performance issues.
            try {
                root = boost::math::tools::bracket_and_solve_root(func, guess, factor, isRising, utils::BracketTolerance, it); // find root
                // root finder returned without raising an exception
                if (it >= maxit) { SHOW_WARN(ERROR::TOO_MANY_MASS0_ITERATIONS); }                           // too many root finder iterations
            }
            catch(std::exception& e) {                                                                      // catch generic boost root finding error
                // root finder exception
                // could be too many iterations, or unable to bracket root - it may not
                // be a hard error - so no matter what the reason is that we are here,
                // we'll just emit a warning and keep trying
                if (it >= maxit) { SHOW_WARN(ERROR::TOO_MANY_MASS0_ITERATIONS); }                           // too many root finder iterations
                else             { SHOW_WARN(ERROR::ROOT_FINDER_FAILED, e.what()); }                        // some other problem - show it as a warning
            }

            // we have a solution from the root finder - it may not be an acceptable solution
            // so we check if it is within our preferred tolerance
            if (fabs(func(root.first + (root.second - root.first) / 2.0)) <= ROOT_ABS_TOLERANCE) {          // solution within tolerance?
                done = true;                                                                                // yes - we're done
            }
            else if (fabs(func(root.first)) <= ROOT_ABS_TOLERANCE) {                                        // solution within tolerance at endpoint 1?
                root.second=root.first;
                done = true;                                                                                // yes - we're done
            }
            else if (fabs(func(root.second)) <= ROOT_ABS_TOLERANCE) {                                       // solution within tolerance at endpoint 2?
                root.first=root.second;
                done = true;                                                                                // yes - we're done
            }
            else {                                                                                          // no - try again
                // we don't have an acceptable solution - reduce search step size and try again
                factorFrac /= 2.0;                                                                          // reduce fractional part of factor
                factor      = 1.0 + factorFrac;                                                             // new search step size
                tries++;                                                                                    // increment number of tries
                if (tries > ADAPTIVE_MASS0_MAX_TRIES || fabs(factor - 1.0) <= ROOT_ABS_TOLERANCE) {         // too many tries, or step size 0.0?
                    // we've tried as much as we can - fail here with -ve return value
                    root.first  = -1.0;                                                                     // yes - set error return
                    root.second = -1.0;
                    SHOW_WARN(ERROR::TOO_MANY_MASS0_TRIES);                                                 // show warning
                    done = true;                                                                            // we're done
                }
            }
        }
        
        return root.first + (root.second - root.first) / 2.0;                                               // midway between brackets is our result, if necessary we could return the result as an interval here.
    }
};



///// INLINE CANDIDATES <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           COEFFICIENTS / CONSTANTS etc.                           //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRho_Hurley2000_Static
 *
 * @brief
 * Calculate the parameter rho for the Hertzsprung Gap, per Hurley et al. 2000, eq 29.
 *
 * Rho is calculated such that McTMS = Rho * McEHG
 * 
 *
 * static double CalculateRho_Hurley2000_Static(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      HG rho
 */
GNU_CONST inline double HG::CalculateRho_Hurley2000_Static(const double p_Mass) const {
    const double m5_25 = utils::intPow(p_Mass, 5) * std::sqrt(std::sqrt(p_Mass)); // sqrt() is much faster than pow()
    return (1.586 + m5_25) / (2.434 + (1.02 * m5_25));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         AGE / LIFETIME / TAU / TIMESCALES                         //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateAgeAfterMassLoss_Hurley2000
 *
 * @brief
 * Recalculate the star's age after mass loss, per Hurley et al. 2000, section 7.1
 *
 *
 * double CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tMS, const double p_tBGB) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tMS                           MS lifetime, tMS (per Hurley timescales) (Myr)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (per Hurley timescales) (Myr)
 * @return                                      Age of the star after mass loss (Myr)
 */
COMPAS_PURE inline double HG::CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tMS, const double p_tBGB) const {

    const double tBGBprime = BaseStar::CalculateLifetimeToBGB_Hurley2000_Static(p_Mass);
    const double tMSprime  = MainSequence::CalculatePhaseLifetime_Hurley2000(p_Mass, tBGBprime);

    return tMSprime + (((tBGBprime - tMSprime) / (p_tBGB - p_tMS)) * (p_Age - p_tMS));
}


/*
 * CalculateTau_Hurley2000
 *
 * @brief
 * Calculate the HG-relative age (fractional Hertzsprung Gap age) of the star,
 * per Hurley et al. 2000, eq 25
 *
 * 
 * double CalculateTau_Hurley2000(const double p_Age, const double p_tMS, const double p_tBGB) const
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tMS                           MS lifetime, tMS (per Hurley timescales) (Myr)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (per Hurley timescales) (Myr)
 * @return                                      HG-relative age, [0, 1]
 */
GNU_CONST inline double HG::CalculateTau_Hurley2000(const double p_Age, const double p_tMS, const double p_tBGB) const {
    return std::max(0.0, std::min(1.0, (p_Age - p_tMS) / (p_tBGB - p_tMS)));
}


/*
 * CalculateAgeAfterMassLoss  <<<<<<<<<<<< BASESTAR???????????????????? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
inline double HG::CalculateAgeAfterMassLoss() const {

    double age;

    Switch (OPTIONS->Mode()) {                                                                                  // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                                             // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                                             // HURLEY BSE
            age = CalculateAgeAfterMassLoss_Hurley2000(m_StateHistory.CurrentState.MassEffectiveInitial(),
                                                       m_StateHistory.CurrentState.Age(),
                                                       m_StateHistory.CurrentState.Timescales(tMS),
                                                       m_StateHistory.CurrentState.Timescales(tBGB)
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
 * double CalculateCoreMassAtPhaseEnd_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const {
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      TAHG core mass (Base of the Giant Branch) (Msol)
 */
COMPAS_PURE inline double HG::CalculateCoreMassAtPhaseEnd_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const {

    double McTAHG;

    if (p_Mass < GLOBALS->HurleyMassCutoffs(static_cast<int>(HURLEY_MASS_CUTOFF::MHeF))) {
        McTAHG = BaseStar::CalculateCoreMass_Hurley2000_Static(GiantBranch::CalculateLuminosityAtBGB_Hurley2000_Static(p_Mass), p_GBparams);
    }
    else if (p_Mass < GLOBALS->HurleyMassCutoffs(static_cast<int>(HURLEY_MASS_CUTOFF::MFGB))) {
        McTAHG = p_GBparams(static_cast<int>(HURLEY_GBP:::McBGB));
    }
    else {
        McTAHG = CalculateCoreMassAtHeI_Hurley2000(p_Mass);
    }

    return McTAHG;
}


/*
 * CalculateCoreMass_Hurley2000
 *
 * Calculate core mass on the Hertzsprung Gap, per Hurley et al. 2000, eq 30
 * (see Section 7
 *
 * The core mass returned is clamped to a minimum value of the current value of the
 * core mass of the star (passed as a parameter) - if the star is losing mass we don't
 * reduce the size of the core (per Hurley).
 *
 *
 * double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_CoreMass                      Current core mass of the star (Msol)
 * @return                                      HG core mass (Msol)
 */
GNU_CONST inline double HG::CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau, const double p_CoreMass) const {
    return std::max(CalculateCoreMass_Hurley2000_Unconstrained_Static(p_Mass, p_Tau), p_CoreMass);
}


/*
 * CalculateCoreMass_Hurley2000_Unconstrained_Static
 *
 * Calculate core mass on the Hertzsprung Gap without accounting for previous core mass.
 *
 * This ignores the previous core mass constraint (see section 7 of Hurley et al. 2000)
 * when computing the expected core mass, and just follows eq. 30.  This is useful for
 * asking what the core mass would be for the given mass without considering that the
 * core mass should not be allowed to drop.
 *
 *
 * static double CalculateCoreMass_Hurley2000_Unconstrained_Static(const double p_Mass, const double p_Tau) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      HG core mass (Msol)
 */
GNU_CONST inline double HG::CalculateCoreMass_Hurley2000_Unconstrained_Static(const double p_Mass, const double p_Tau) const {
    const double McEHG = CalculateCoreMassAtPhaseEnd_Hurley2000(p_Mass, p_GBparams, p_MHeF, p_MFGB, p_aN);
    return (((1.0 - p_Tau) * CalculateRho_Hurley2000_Static(p_Mass)) + p_Tau) * McEHG;
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
 * double CalculateEffectiveInitialMass_Hurley2000(const double p_Mass, const double p_Tau, const double p_Coremass, const double p_MassEffectiveInitia) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MassEffectiveInitial          Current effective initial mass of the star (Msol)
 * @return                                      HG effective initial mass (Msol) (may be unchanged)
 */
GNU_CONST inline double HG::CalculateEffectiveInitialMass_Hurley2000(const double p_Mass, const double p_Tau, const double p_Coremass, const double p_MassEffectiveInitial) const {
    return ((p_MassEffectiveInitial > p_Mass) && (p_Coremass <= CalculateCoreMass_Hurley2000_Unconstrained_Static(p_Mass, p_Tau))) ? p_Mass : p_MassEffectiveInitial;
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
 * static double CalculateLuminosityAtPhaseEnd_Hurley2000_Static(const double p_Mass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAHG luminosity (Lsol)
 */
COMPAS_PURE inline double HG::CalculateLuminosityAtPhaseEnd_Hurley2000_Static(const double p_Mass) const {
    return p_Mass < GLOBALS->HurleyMassCutoffs(static_cast<int>(HURLEY_MASS_CUTOFF::MFGB))
            ? GiantBranch::CalculateLuminosityAtBGB_Hurley2000_Static(p_Mass)
            : GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass,);
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
 * @param       p_Tau                           Phase-relative age of the star (Myr)
 * @return                                      HG luminosity (Lsol)
 */
GNU_CONST inline double HG::CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Tau) const {
    const double lTMS = MainSequence::CalculateLuminosityAtPhaseEnd_Hurley2000(p_Mass);
    return lTMS * PPOW((CalculateLuminosityAtPhaseEnd_Hurley2000_Static(p_Mass) / lTMS), p_Tau);
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
 * double CalculateRadiusAtPhaseEnd_Hurley2000(const double p_Mass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAHG radius (Rsol)
 */
COMPAS_PURE inline double HG::CalculateRadiusAtPhaseEnd_Hurley2000(const double p_Mass, const double p_MFGB, const DBL_VECTOR& p_aN) const {
    const double mFGB = GLOBALS->HurleyMassCutoffs(static_cast<int>(HURLEY_MASS_CUTOFF::MFGB));
    return p_Mass < mFGB
            ? GiantBranch::CalculateRadius_Hurley2000_Static(p_Mass, GiantBranch::CalculateLuminosityAtBGB_Hurley2000_Static(p_Mass))
            : CalculateRadiusAtHeIgnition_Hurley2000(p_Mass);
}



/// HG_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


class HG_Constituent: virtual public BinaryConstituentStar, public HG {

public:


protected:

};

     
    MT_CASE DetermineMassTransferTypeAsDonor() const { return MT_CASE::B; }                                                                                                     // Always case B

    

    double          CalculateMTRejuvenationFactor()                       { return 1.0; }


    void            UpdateAfterMerger(double p_Mass, double p_HydrogenMass) { }                                                                                                 // Nothing to do for stars beyond the Main Sequence for now

    
    double          CalculateZetaEquilibrium()                                      { return -std::numeric_limits<double>::infinity(); }                                         // Nuclear timescale MT should be impossible from HG stars that evolve on a thermal timescale


    
    double          CalculateCELambda_Dewi() const;
    double          CalculateLambdaLoveridge(const double p_EnvMass, const bool p_IsMassLoss = false) const;
    double          CalculateLambdaNanjingStarTrack(const double p_Mass) const;
    double          CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;




double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const;

double CalculateCriticalMassRatio_Hurley2002() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_HG; }



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


#endif // __HG_h__
