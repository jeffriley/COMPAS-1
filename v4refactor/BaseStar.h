#pragma once

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"
#include "vector3d.h"

#include "Rand.h"
#include "Options.h"
#include "Log.h"
#include "Errors.h"

#include "State.h"
#include "Globals.h"
#include "ZDependent.h"

#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_cdf.h>

#include <boost/math/distributions.hpp>

class Star;


class BaseStar {
    
protected:

    // FixedAttributes
    //
    // Inner class holding the attributes of a star that are set once at construction
    // and never change thereafter. These are NOT state - they are fixed properties of
    // the star itself, independent of the star's evolution.  The test is whether the
    // value of an attribute can mutate over evolution - if it can, it should be in
    // state; if not, it should be here.  The values recorded here are typically
    // specified by the user at run time, or calculated from values specified by the
    // user at run time.
    //
    // The "fixed" guarantee is enforced by the class having no setters exposed - all
    // values are populated at construction and are immutable thereafter.
    //
    // RandomSeed could conceivably be part of GLOBALS, but if we keep it here that opens
    // the possibility that binary stars could have their own random seed (that might
    // differ). That may not be somethig we ever implement, but I like to keep things as
    // flexible as possible.

    class FixedAttributes {
    
    private:
        
        // Always present
        ULongT       m_ObjectId;                // Object's unique object id - always present
        ULongT       m_RandomSeed;              // Random seed for the star - always present
        STELLAR_TYPE m_StartingStellarType;     // Star's starting stellar type - always present

        // Baryonic mass for which the gravitational remnant mass equals the max NS mass.
        const double maxMassNS                 = OPTIONS->MaximumNeutronStarMass();
        const double m_BaryonicMassOfMaxMassNS = (0.075 * maxMassNS * maxMassNS) + maxMassNS;

        // Optional - may not have a value assigned
        OptDblT      m_UserKickMagnitude;       // User supplied kick magnitude (km s^-1) - may not contain a value
        OptDblT      m_UserKickMagnitudeRandom; // User supplied kick random number U(0,1) - may not contain a value

        OptDblT      m_RotationalFrequency;     // User supplied rotational frequency (Hz) - may not contain a value
        OptDblT      m_CHEOmega;                // Minimum angular frequency for CHE (rad yr^-1) - may not contain a value

        // Pulsar constants (optional).
        // These are calculated from program options and don't change throughout evolution.
        //
        // Since they are only used when the star is a NS, we don't bother computing the
        // values until (and if) that happens - they are computed at first request (of any
        // of them).  Until then they don't contain values.
        //
        // The struct pulsarConstants will not contain a value unless the pulsar constants
        // have been calculated.  Callers should check pulsarConstants.has_value() if they
        // are unsure if the values have been calculated. If the struct exists, then each
        // of the constitunets will have been calculated.
        //
        // The has_value() check for std::optional types is extemely fast - typically a single
        // machine instruction.

        std::optional<PulsarConstants> m_PulsarConstants;


        /* 
         * CalculatePulsarConstants
         *
         * Brief
         * Calculates the NS fixed attributes from option values.
         *
         * This function mutates FixedAttributes member variables directly.  It's private,
         * as are the member variables.
         * 
         * 
         * void CalculatePulsarConstants()
         */
        void CalculatePulsarConstants() {
        
            PulsarConstants pConsts = p_pConsts;

            // NS decay constants
            pConsts.decayMassScale = OPTIONS->PulsarMagneticFieldDecayMassscale() * MSOL_TO_G;
            pConsts.decayTimeScale = OPTIONS->PulsarMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR;

            // NS minimum magnetic field
            pConsts.log10MinMagField = OPTIONS->PulsarLog10MinimumMagneticField();
            pConsts.minMagField      = PPOW(10.0, pConsts.log10MinMagField);

            return pConsts;
        }


    public:

        FixedAttributes(const ULongT       p_RandomSeed,
                        const STELLAR_TYPE p_StartingStellarType,
                        const double       p_Metallicity,
                        const double       p_Mass,
                        const OptDblT      p_UserKickMagnitude,
                        const OptDblT      p_UserKickMagnitudeRandom,
                        const OptDblT      p_RotationalFrequency) {

            // Always present
            m_ObjectId            = GLOBALS->NextObjectId();        // Always present
            m_RandomSeed          = p_RandomSeed;                   // Always present
            m_StartingStellarType = p_StartingStellarType;          // Always present

            // Optional - may not have a value assigned - the user may not have specified them
            if (p_UserKickMagnitude.has_value())       m_UserKickMagnitude       = p_UserKickMagnitude;
            if (p_UserKickMagnitudeRandom.has_value()) m_UserKickMagnitudeRandom = p_UserKickMagnitudeRandom;
            if (p_RotationalFrequency.has_value())     m_RotationalFrequency     = p_RotationalFrequency;

            // Only set these values if the star starts evolution on the main sequence
            if (utils::IsOneOf(p_StartingStellarType, MAIN_SEQUENCE)) {
                m_CHEOmega = astro::CalculateMinOmegaForCHE_Butler2018(p_Metallicity, p_Mass);
            }
        }
   

        // Getters

        // Always present
        ULongT       ObjectId() const                   { return m_ObjectId; }
        ULongT       RandomSeed() const                 { return m_RandomSeed; }
        STELLAR_TYPE StartingStellarType() const        { return m_StartingStellarType; }

        double       BaryonicMassOfMaxMassNS()          { return m_BaryonicMassOfMaxMassNS; }

        // Optional - may not have an assigned value.
        // Getters for these always return std::optional<double>
        OptDblT      CHEOmega() const                   { return m_CHEOmega; }
        OptDblT      RotationalFrequency() const        { return m_RotationalFrequency; }
        OptDblT      UserKickMagnitude() const          { return m_UserKickMagnitude; }
        OptDblT      UserKickMagnitudeRandom() const    { return m_UserKickMagnitudeRandom; }

        // Pulsar constant getters.
        // Each causes all pulsar constants to be calculated if necessary.
        // Getters for the individual values always return a double - they are
        // calculated on demand, so are guaranteed to have a value assigned.

        PulsarConstants PulsarConstants() const {
            if (!m_PulsarConstants.has_value()) m_PulsarConstants = CalculatePulsarConstants();
            return m_PulsarConstants;
        }

        double PulsarDecayMassScale() const {
            if (!m_PulsarConstants.has_value()) m_PulsarConstants = CalculatePulsarConstants();
            return m_PulsarConstants.decayMassScale.value();
        }

        double PulsarDecayTimeScale() const {
            if (!m_PulsarConstants.has_value()) m_PulsarConstants = CalculatePulsarConstants();
            return m_PulsarConstants.decayTimeScale.value();
        }

        double PulsarMinMagField() const {
            if (!m_PulsarConstants.has_value()) m_PulsarConstants = CalculatePulsarConstants();
            return m_PulsarConstants.minMagField.value();
        }

        double PulsarLog10MinMagField() const {
            if (!m_PulsarConstants.has_value()) m_PulsarConstants = CalculatePulsarConstants();
            return m_PulsarConstants.log10MinMagField.value();
        }
    };


    Star*                   m_Parent = nullptr;                             // This object's parent object

    OBJECT_PERSISTENCE      m_ObjectPersistence;                            // Instantiated object's persistence (permanent or ephemeral)

    FixedAttributes         m_FixedAttributes;                              // Fixed star attributes

    StateHistory<StarState> m_StateHistory;                                 // Empty start state already present
    StarState               m_InterimState;                                 // Interim state being built during a timestep (see AdvanceOneTimestep)



    // The interim state being built during the current timestep (see AdvanceOneTimestep).
    // Read this for THIS-timestep values of dt, dMdt, GBParams, timescales.
    // read CurrentState() for frozen start-of-timestep attributes.
    const StarState& InterimState()  const { return m_InterimState; }

public:

    BaseStar(const Star*              p_Parent,
             const ULongT             p_RandomSeed,
             const STELLAR_TYPE       p_StartingStellarType,
             const double             p_Metallicity,
             const double             p_Mass,
             const StellarKickParmsT& p_KickParameters,
             const OptDblT            p_RotationalFrequency = std::nullopt) {

        // m_Parent is the parent Star object that encapsulates this BaseStar object.
        // It is here to facilitate calls to Star::Clone(), Star::CloneAs(), and Star::MakeStar()
        m_Parent = p_Parent;

        // m_ObjectPersistence is the *expected* persistence - PERMANENT or EPHEMERAL. The value
        // given here is just an indication, set by the caller. Logging code uses it to decide
        // whether the object participates in logging. Cloned objects typically change this to
        // EPHEMERAL so they don't participate in logging.
        m_ObjectPersistence = OBJECT_PERSISTENCE::PERMANENT;

        // Populate fixed attributes
        m_FixedAttributes = FixedAttributes(
            p_RandomSeed,
            p_StartingStellarType,
            p_Metallicity,                                                      // Not stored in FixedAttributes - used to calculate CHEOmega
            p_Mass,                                                             // Not stored in FixedAttributes - used to calculate CHEOmega
            p_KickParameters.userMagnitude,
            p_KickParameters.userMagnitudeRandom,
            p_RotationalFrequency
        );

        StarState s;                                                            // Start state

        // Set state attributes supplied as parameters
        s.SetMetallicity(p_Metallicity);                                        // Can mutate over evolution
        s.SetMass(p_Mass);                                                      // Also ZAMS mass if starting on MS

        // Calculate and set angular frequency if use specified rotational frequency
        if (p_RotationalFrequency.has_value()) s.SetAngularFrequency(SECONDS_IN_YEAR * _2_PI * p_RotationalFrequency);

        // Values set only if starting on the Main Sequence (this is then also the
        // ZAMS state). Cannot be CHE here - CHE is decided later by size + rotation.
        if (utils::IsOneOf(p_StellarType, MAIN_SEQUENCE)) {

            s.SetMass0(p_Mass);

            s.SetAge(0.0);
            s.SetTau(0.0);
            s.SetTime(0.0);
           
            s.SetHAbundanceCore(ZDEP->ZAMSHAbundance(p_Metallicity));
            s.SetHAbundanceSurface(s.HAbundanceCore());                         // = core at ZAMS
            s.SetHeAbundanceCore(ZDEP->ZAMSHeAbundance(p_Metallicity));
            s.SetHeAbundanceCoreOut(s.HeAbundanceCore());                       // Helium abundance just outside the core (= core at ZAMS)
            s.SetHeAbundanceSurface(s.HeAbundanceCore());                       // = core at ZAMS

            s.SetLuminosity(astro::CalculateZAMSLuminosity_Tout1996(m_Metallicity, m_Mass));
            s.SetLuminosity0(s.Luminosity());

            s.SetRadius(astro::CalculateZAMSRadius_Tout1996(m_Metallicity, m_Mass));
            s.SetRadius0(s.Radius());

            s.SetTemperatureZAMSEffective(astro::CalculateTemperature(s.Luminosity(), s.Radius()));
            s.SetTemperature0(s.Temperature());

            s.SetCNOProcessedCoreMass(0.0);                                     // Default for MS (changed in MS_gt_07 constructor)
            s.SetCOCoreMass(0.0);                                               // Default for MS
            s.SetCoreMass(0.0);                                                 // Default for MS (changed in MS_gt_07 constructor)
            s.SetHeCoreMass(0.0);                                               // Default for MS

            // Since we know the star is starting evolution on the main sequence, we can calculate
            // the ZAMS angular frequency if necessary.  We only do this if the user did not provide
            // a value for rotational frequency - if they did, we've already used that to calculate
            // the angular frequency (see above).
            if (!p_RotationalFrequency.has_value()) {                           // Rotational frequency supplied?
                                                                                // No - calculate ZAMS omega
                s.SetAngularFrequency(astro::CalculateZAMSOmega(OPTIONS->Mode(), p_Mass, m_FixedAttributes.RZAMS()));
            }

            // Since we know angular frequency here we can calculate angular momentum
            // *Ilya* is this still valid when Brcek core mass is used?
            s.SetAngularMomentum(CalculateMomentOfInertiaAU() * s.AngularFrequency());
        }

        // Initial values for evolution variables
        s.SetError(ERROR::NONE);
        s.SetEvolutionStatus(EVOLUTION_STATUS::CONTINUE);

        // Commit the start state as the star's start (and current, and possibly ZAMS) state
        m_StateHistory.SetStartState(s);
    }
        
    BaseStar(const BaseStar& p_Star) = default;

    virtual ~BaseStar() {}

 

    
    
    
    // Member functions
    //
    // Declarations are here, but some implementations might be also:
    //
    // Short (one or two line) implementations that don't require a lot of documentation (i.e. that won't clutter here),
    // are in here (and are implied inline).
    //
    // Short implementations that can still be inlined, but that require more documentation (and would clutter here),
    // are below (after class declaration).
    //
    // Longer implementations that can still be inlined (that would clutter here) are below (after class declaration).
    //
    // All other implementations are in the cpp file.
    //
    // Function names should indicate the purpose of the function (and should genarally begin with a verb).
    //
    // e.g., functions beginning with:
    //
    //    - "Calculate" should calculate and return value(s) 
    //    - "Determine" should determine e.g. current state, envelope type, etc., and return value(s)
    //    - "Draw" should draw a value, or values, from a distribution and return value(s)
    //
    // Almost all non getter/setter functions should do some operation ("Calculate", "Determine", "Draw", etc.),
    // and should return the result(s) of that operation - they should not modify class member or state variables.
    // The exceptions are setters, and there should be very few of those.  Class member and/or state variables
    // should be modified in as few places in the code as possible - this helps:
    //
    //     (a) to keep the code clean,
    //     (b) developers to be clear about what functions do and any side-effects they may have,
    //     (c) to facilitate our "one (major) thing at a time" paradigm
    //
    // Getter function names are generally just the name of (or descriptive version of) the variable for which
    // the value is required.
    //
    // Setter function names should begin with "Set", followed by the name of (or descriptive version of) the
    // variable being modified.
    //
    // Functions that do not modify class member variables should be declared const.
    //
    // Where appropriate, the attributes GNU_CONST, GNU_PURE, and COMPAS_PURE should be applied to function
    // declarations.  See <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    
    // Object identifiers - all classes have these
                OBJECT_ID           ObjectId() const                            { return m_FixedAttributes.ObjectId(); }
                ULongT              RandomSeed() const                          { return m_FixedAttributes.RandomSeed(); }
                OBJECT_PERSISTENCE  ObjectPersistence() const                   { return m_ObjectPersistence; }
    GNU_CONST   OBJECT_TYPE         ObjectType() const                          { return OBJECT_TYPE::BASE_STAR; }


    // JR FIX THIS - DESCRIBE THE CLONE FUNCTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    std::unique_ptr<BaseStar> Clone(const std::unique_ptr<BaseStar> p_Star, const STELLAR_TYPE p_StellarType, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) const;
    std::unique_ptr<BaseStar> CloneAs(const STELLAR_TYPE p_StellarType, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) const;
    std::unique_ptr<BaseStar> MakeStar(const STELLAR_TYPE p_StellarType, const OBJECT_PERSISTENCE p_Persistence) const;


    // Getters - alphabetically
    //
    // Most of these are just convenience functions to retrieve state variables.
    // Getter implementations should generally be very short, so many will be implemented here.
    //
    //
    // Most of the state variable getters here will get current-state variables, but some will get
    // interim-state variables.  The rationale for this is, and for each interim-state variable,
    // is given below.
    //
    // The iterative method used by COMPAS pre v4.00.00 was an "asynchronous" update-ordering
    // method - meaning that the attributes of a system were updated in-place, and then those
    // already-updated attributes were used to update other attributes *in the same timestep*.  
    // Depending on the code path taken, an attributes dependencies (i.e. the attribute values
    // used to calculate the value of the attribute), may or may not have yet been updated from
    // the previous timestep - so different code paths could have resulted in different
    // attribute values simply by virtue of the code path (not necessarily for any physical
    // reason).  This meant that the starting point was not consisent for all attribute updates.
    // The value of some attributes was calculated purely from the previous timestep values,
    // whereas for aother attributes the value was calculated using atrribute values that had
    // already been updated in the current timestep.
    //
    // To ensure that the values of all attribues are updated from a consistent starting point,
    // from v4.00.00 onwards COMPAS behaves as a transition function applied by an integrator
    // (a bit like a state machine, but technically not quite).  At each timestep, a pure function
    // maps the current state to a new state.  To do that, at the start of each timestep, we first
    // create an "interim" state from the current state:
    //
    //     interimState = f(currentState),
    //
    // where
    //
    //     interimState.Metallicity = metallicity to be used for the timestep  <<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS <<<<<<<<<<<<<<<<<<
    //     interimState.GBParams    = f(interimState, interimState.GBParams)
    //     interimState.Timescales  = f(interimState, interimState.Timescales)
    //     interimState.dt          = f(interimState)  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS - p_dt may be supplied <<<<<<<<<<<<<<<<<<<
    //     interimState.dMdt        = f(interimState)  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS - d_MT may be supplied <<<<<<<<<<<<<<<<<<<
    
    // — newState = f(currentState, interimState) — and the evolution loop (the integrator) just
    // applies that function over and over to march the system forward through time. 
    // The organising principle is continuous state being integrated forward, not a graph of discrete states.
    // That's why "state machine" doesn't quite fit: the thing being advanced is continuous
    // (mass, radius, luminosity, age…), and you're numerically integrating it, not hopping between discrete modes.

    // Some getters are virtual.
    //
    // VIRTUAL FUNCTIONS may be (are expected to be) overridden by derived classes.
    // When overriding virtual functions in a derived class, use the "override" attribute.
    //
    // NON-VIRTUAL FUNCTIONS should not be overridden (declared separately) by derived classes.
    // While it is legal in C++ to declare the same (non-virtual) function in multiple classes,
    // (aka "shadowing", or "hiding"), we discourage it.  Non-virtual functions are statically
    // bound, and as such, especially with indirection, may not produce expected results.

    // Virtual getters are grouped here
    virtual bool      IsDegenerate() const                     { return false; }   // Virtual - default is not degenerate - White Dwarfs, NS and BH are degenerate
    virtual bool      IsSupernova() const                      { return false; }   // Virtual - default value


    // Non-virtual getters
    double            Age() const                              { return m_InterimState.Age(); }

    double            AngularFrequency() const                 { return m_StateHistory.CurrentState().AngularFrequency(); }
    double            Omega() const                            { return AngularFrequency(); }
    OptDblT           ZAMSOmega() const                        { return m_StateHistory.HaveZAMSState() ? m_StateHistory.ZAMSState().AngularFrequency() : std::nullopt; }

    double            AngularMomentum() const                  { return m_StateHistory.CurrentState().AngularMomentum(); }

    bool              EvolvedCHE() const                       { return m_StateHistory.CurrentState().EvolvedCHE(); }

    double            CNOProcessedCoreMass() const             { return m_StateHistory.CurrentState().CNOProcessedCoreMass(); }
    OptDblT           ZAMSCNOProcessedCoreMass() const         { return m_StateHistory.HaveZAMSState() ? m_StateHistory.ZAMSState().CNOProcessedCoreMass() : std::nullopt; }

    double            CHEOmega() const                         { return m_FixedAttributes.CHEOmega(); }

    double            COCoreMass() const                       { return m_StateHistory.CurrentState().COCoreMass(); }
    double            CoreMass() const                         { return m_StateHistory.CurrentState().CoreMass(); }

    ML_TYPE           DominantMassLossType() const             { return m_StateHistory.CurrentState().DominantMassLossType(); }
   
    double            dMdt() const                             { return m_InterimState.dMdt(); }
    double            dt() const                               { return m_InterimState.dt(); }

    OptDblT           dtPrev() const                           { return m_StateHistory.HavePreviousState() ? m_StateHistory.PreviousState().dt() : std::nullopt; }

    bool              EnvExpelledByPulsations() const          { return m_StateHistory.CurrentState().EnvExpelledByPulsations(); }

    ERROR             Error() const                            { return m_StateHistory.CurrentState().Error(); }

    EVOLUTION_STATUS  EvolutionStatus() const                  { return m_StateHistory.CurrentState().EvolutionStatus(); }

    bool              ExperiencedAIC() const                   { return (m_StateHistory.CurrentState().SNEvents().past & SN_EVENT::AIC)   == SN_EVENT::AIC; }
    bool              ExperiencedCCSN() const                  { return (m_StateHistory.CurrentState().SNEvents().past & SN_EVENT::CCSN)  == SN_EVENT::CCSN; }
    bool              ExperiencedECSN() const                  { return (m_StateHistory.CurrentState().SNEvents().past & SN_EVENT::ECSN)  == SN_EVENT::ECSN; }
    bool              ExperiencedHeSD() const                  { return (m_StateHistory.CurrentState().SNEvents().past & SN_EVENT::HeSD)  == SN_EVENT::HeSD; }
    bool              ExperiencedPISN() const                  { return (m_StateHistory.CurrentState().SNEvents().past & SN_EVENT::PISN)  == SN_EVENT::PISN; }
    bool              ExperiencedPPISN() const                 { return (m_StateHistory.CurrentState().SNEvents().past & SN_EVENT::PPISN) == SN_EVENT::PPISN; }
    bool              ExperiencedSNIA() const                  { return (m_StateHistory.CurrentState().SNEvents().past & SN_EVENT::SNIA)  == SN_EVENT::SNIA; }
    bool              ExperiencedSN_NONE() const               { return  m_StateHistory.CurrentState().SNEvents().past                    == SN_EVENT::NONE; }
    SN_EVENT          ExperiencedSN_Type() const               { return utils::SNEventType(m_StateHistory.CurrentState().SNEvents().past); }
    bool              ExperiencedUSSN() const                  { return (m_StateHistory.CurrentState().SNEvents().past & SN_EVENT::USSN)  == SN_EVENT::USSN; }

    OptDblT           GBParam(const SizeT p_GBParam) const     { return m_InterimState.GBParams(p_GBParam, get::asOptionalT); }
    GBParamsT         GBParams() const                         { return m_InterimState.GBParams(); }

    double            HAbundanceCore() const                   { return m_StateHistory.CurrentState().HAbundanceCore(); }
    double            HAbundanceSurface() const                { return m_StateHistory.CurrentState().HAbundanceSurface(); }
    double            HeAbundanceCore() const                  { return m_StateHistory.CurrentState().HeAbundanceCore(); }
    double            HeAbundanceSurface() const               { return m_StateHistory.CurrentState().HeAbundanceSurface();}
    double            HeAbundanceCoreOut() const               { return m_StateHistory.CurrentState().HeAbundanceCoreOut(); }

    double            HeCoreMass() const                       { return m_StateHistory.CurrentState().HeCoreMass(); }

    double            HurleyMu() const                         { return m_StateHistory.CurrentState().HurleyMu(); }

    double            InitialHAbundance() const                { return ZDEP->ZAMSHAbundance(Metallicity()); }
    double            InitialHeAbundance() const               { return ZDEP->ZAMSHeAbundance(Metallicity()); }

    bool              IsAIC() const                            { return (m_StateHistory.CurrentState().SNEvents().current & SN_EVENT::AIC)   == SN_EVENT::AIC; }
    bool              IsCCSN() const                           { return (m_StateHistory.CurrentState().SNEvents().current & SN_EVENT::CCSN)  == SN_EVENT::CCSN; }
    bool              IsECSN() const                           { return (m_StateHistory.CurrentState().SNEvents().current & SN_EVENT::ECSN)  == SN_EVENT::ECSN; }
    bool              IsHeSD() const                           { return (m_StateHistory.CurrentState().SNEvents().current & SN_EVENT::HeSD)  == SN_EVENT::HeSD; }
    bool              IsPISN() const                           { return (m_StateHistory.CurrentState().SNEvents().current & SN_EVENT::PISN)  == SN_EVENT::PISN; }
    bool              IsPPISN() const                          { return (m_StateHistory.CurrentState().SNEvents().current & SN_EVENT::PPISN) == SN_EVENT::PPISN; }
    bool              IsSN_NONE() const                        { return  m_StateHistory.CurrentState().SNEvents().current                    == SN_EVENT::NONE; }
    SN_EVENT          IsSN_Type() const                        { return utils::SNEventType(m_StateHistory.CurrentState().SNEvents().current); }
    bool              IsSNIA() const                           { return (m_StateHistory.CurrentState().SNEvents().current & SN_EVENT::SNIA)  == SN_EVENT::SNIA; }
    bool              IsUSSN() const                           { return (m_StateHistory.CurrentState().SNEvents().current & SN_EVENT::USSN)  == SN_EVENT::USSN; }

    double            Luminosity() const                       { return m_StateHistory.CurrentState().Luminosity(); }
    OptDblT           LuminosityPrev() const                   { return m_StateHistory.HavePreviousState() ? m_StateHistory.PreviousState().Luminosity() : std::nullopt; }
    double            LuminosityStart() const                  { return m_StateHistory.StartState().Luminosity(); }
    double            LuminosityZAMSEffective() const          { return m_StateHistory.CurrentState().LuminosityZAMSEffective(); }
    double            Luminosity0() const                      { return m_StateHistory.CurrentState().Luminosity0(); }
    OptDblT           LZAMS() const                            { return m_StateHistory.HaveZAMSState() ? m_StateHistory.ZAMSState().Luminosity() : std::nullopt; }

    double            Mass() const                             { return m_StateHistory.CurrentState().Mass(); }
    double            MassEffectiveInitial() const             { return m_StateHistory.CurrentState().MassEffectiveInitial(); }
    double            Mass0() const                            { return m_StateHistory.CurrentState().Mass0(); }
    OptDblT           MassPrev() const                         { return m_StateHistory.HavePreviousState() ? m_StateHistory.PreviousState().Mass() : std::nullopt; }
    double            MassStart() const                        { return m_StateHistory.StartState().Mass(); }
    OptDblT           MZAMS() const                            { return m_StateHistory.HaveZAMSState() ? m_StateHistory.ZAMSState().Mass() : std::nullopt; }

    double            Metallicity() const                      { return m_InterimState.Metallicity(); }

    OptDblT           PulsarMagneticField() const              { PulsarDetailsT& p = m_StateHistory.CurrentState().PulsarDetails(); return p.has_value() ? p.magneticField : std::nullopt; }
    OptDblT           PulsarSpinPeriod() const                 { PulsarDetailsT& p = m_StateHistory.CurrentState().PulsarDetails(); return p.has_value() ? p.spinPeriod : std::nullopt; }
    OptDblT           PulsarSpinFrequency() const              { PulsarDetailsT& p = m_StateHistory.CurrentState().PulsarDetails(); return p.has_value() ? p.spinFrequency : std::nullopt; }
    OptDblT           PulsarSpinDownRate() const               { PulsarDetailsT& p = m_StateHistory.CurrentState().PulsarDetails(); return p.has_value() ? p.spinDownRate : std::nullopt; }
    OptDblT           PulsarBirthPeriod() const                { PulsarDetailsT& p = m_StateHistory.CurrentState().PulsarDetails(); return p.has_value() ? p.birthPeriod : std::nullopt }
    OptDblT           PulsarBirthSpinDownRate() const          { PulsarDetailsT& p = m_StateHistory.CurrentState().PulsarDetails(); return p.has_value() ? p.birthSpinDownRate : std::nullopt; }

    double            Radius() const                           { return m_StateHistory.CurrentState().Radius(); }
    OptDblT           RadiusPrev() const                       { return m_StateHistory.HavePreviousState() ? m_StateHistory.PreviousState().Radius() : std::nullopt; }
    double            RadiusStart() const                      { return m_StateHistory.StartState().Radius(); }
    double            RadiusZAMSEffective() const              { return m_StateHistory.CurrentState().RadiusZAMSEffective(); }
    double            Radius0() const                          { return m_StateHistory.CurrentState().Radius0(); }
    OptDblT           RZAMS() const                            { return m_StateHistory.HaveZAMSState() ? m_StateHistory.ZAMSState().Radius() : std::nullopt; }

    ULongT            RandomSeed() const                       { return m_FixedAttributes.RandomSeed(); }

    STELLAR_TYPE      StellarType() const                      { return m_StateHistory.CurrentState().StellarType(); }
    std::optional<STELLAR_TYPE> StellarTypePrev() const        { return m_StateHistory.HavePreviousState() ? m_StateHistory.PreviousState().StellarType() : std::nullopt; }
    STELLAR_TYPE      StellarTypeStart() const                 { return m_StateHistory.StartState().StellarType(); }

    StellarSNDetailsT SupernovaDetails() const                 { return m_StateHistory.CurrentState.SupernovaDetails(); }
    StellarSNDetailsT SN_Details() const                       { return SupernovaDetails(); }

    double            Tau() const                              { return m_InterimState.Tau(); }

    double            Temperature() const                      { return m_StateHistory.CurrentState().Temperature(); }
    OptDblT           TemperaturePrev() const                  { return m_StateHistory.HavePreviousState() ? m_StateHistory.PreviousState().Temperature() : std::nullopt; }
    double            TemperatureStart() const                 { return m_StateHistory.StartState().Temperature(); }
    double            TemperatureZAMSEffective() const         { return m_StateHistory.CurrentState().TemperatureZAMSEffective(); }
    double            Temperature0() const                     { return m_StateHistory.CurrentState().Temperature0(); }
    OptDblT           TZAMS() const                            { return m_StateHistory.HaveZAMSState() ? m_StateHistory.ZAMSState().Temperature() : std::nullopt; }

    double            Time() const                             { return m_InterimState.Time(); }

    OptDblT           Timescale(const SizeT p_Timescale) const { return m_InterimState.Timescales(p_Timescale, get::asOptionalT); }
    TimescalesT       Timescales() const                       { return m_InterimState.Timescales(); }



    virtual double   XExponent() const { return 0.0; }  // <<<<< JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<









 
//////  FOR BINARY STAR CLASS (CONSTITUENT)    StrT                       MassTransferDonorHistory() const; // implementation in cpp file

    bool                        LBV_PhaseFlag() const                                   { return CurrentState().LBV_PhaseFlag(); }
//////  FOR BINARY STAR CLASS (CONSTITUENT)    ST_VECTOR                   MassTransferDonorHistory() const                        { return m_MassTransferDonorHistory; }
    double                      Mdot() const                                            { return CurrentState().Mdot(); }
    GNU_CONST double            CalculateOmegaBreak(const double p_Mass, const double p_Radius) const;
    double                      OmegaBreak() const                                      { return CalculateOmegaBreak(Mass(), Radius()); }
    COMPAS_VARIABLE             PropertyValue(const T_ANY_PROPERTY p_Property) const; // implementation in cpp file

    double                      SN_CoreMassAtCOFormation() const                        { return CurrentState().COFormationDetails().coreMass; }
    double                      SN_CoreRadiusAtCOFormation() const                      { return CurrentState().COFormationDetails().coreRadius; }
    double                      SN_COCoreMassAtCOFormation() const                      { return CurrentState().COFormationDetails().COCoreMass; }
    double                      SN_DrawnKickMagnitude() const                           { return CurrentState().CurrentKickParams().magnitudeDrawn.value_or(0.0); }
    double                      SN_EccentricAnomaly() const                             { return 0.0; }
    double                      SN_FallbackFraction() const                             { return CurrentState().SupernovaDetails().fallbackFraction.value_or(0.0); }
    double                      SN_HeCoreMassAtCOFormation() const                      { return CurrentState().COFormationDetails().HeCoreMass; }
    bool                        SN_IsHydrogenPoor() const                               { return CurrentState().SupernovaDetails().isHydrogenPoor.value_or(false); }
    double                      SN_KickMagnitude() const                                { return CurrentState().CurrentKickParams().magnitude.value_or(0.0); }
    double                      SN_MeanAnomaly() const                                  { return 0.0; }
    double                      SN_Phi() const                                          { return 0.0; }
    double                      SN_RocketKickMagnitude() const                          { return 0.0; }
    double                      SN_RocketKickPhi() const                                { return 0.0; }
    double                      SN_RocketKickTheta() const                              { return 0.0; }
    double                      SN_TotalMassAtCOFormation() const                       { return CurrentState().COFormationDetails().mass; }
    double                      SN_TotalRadiusAtCOFormation() const                     { return CurrentState().COFormationDetails().radius; }
    double                      SN_TrueAnomaly() const                                  { return 0.0; }
    double                      SN_Theta() const                                        { return 0.0; }
    SN_EVENT                    SN_Type() const                                         { return utils::SNEventType(CurrentState().SupernovaEvents().current); }
    double                      SN_KickMagnitudeRandom() const                          { return CurrentState().CurrentKickParams().magnitudeRandom.value_or(0.0); }
    virtual double              Speed() const                                           { return 0.0; }
    COMPAS_VARIABLE             StellarPropertyValue(const T_ANY_PROPERTY p_Property) const; // implmentation in cpp file
    virtual double              TotalMassLossRate() const                               { return 0.0; }
//////  FOR BINARY STAR CLASS (CONSTITUENT)    double                      VelocityX() const                                       { return m_ComponentVelocity.xValue(); }
//////  FOR BINARY STAR CLASS (CONSTITUENT)    double                      VelocityY() const                                       { return m_ComponentVelocity.yValue(); }
//////  FOR BINARY STAR CLASS (CONSTITUENT)    double                      VelocityZ() const                                       { return m_ComponentVelocity.zValue(); }
    virtual double              VelocityX() const                                       { return 0.0; }
    virtual double              VelocityY() const                                       { return 0.0; } 
    virtual double              VelocityZ() const                                       { return 0.0; } 
    virtual ACCRETION_REGIME    WhiteDwarfAccretionRegime() const                       { return ACCRETION_REGIME::ZERO; } // virtual - default value
    
    
    // Setters - alphabetically
    // Setter implementations should generally be very short, so many will be implemented here.

    inline void             SetPersistence(const OBJECT_PERSISTENCE p_Persistence)  { m_ObjectPersistence = p_Persistence; }
    inline void             SetStellarType(const STELLAR_TYPE p_StellarType)       { m_InterimState.SetStellarType(p_StellarType); }
    inline void             SetStartingType(const STELLAR_TYPE p_StellarType)       {
        auto opt = m_StateHistory.StartState();
        if (opt.has_value()) {
            StarState updated = opt.value();
            updated.SetStellarType(p_StellarType);
            m_StateHistory.SetStartState(updated);
        }
    }



    void                        AdvanceAgeAndTime(const double p_dt)                            {
                                                                                                    if (p_dt > 0.0) {                                                               // Only if delta > 0.0 (don't use utils::Compare() here)
                                                                                                        m_InterimState.SetAge(CurrentState().Age()   + p_dt);                       // Advance age of star
                                                                                                        m_InterimState.SetTime(CurrentState().Time() + p_dt);                       // Advance simulation time
                                                                                                    }
                                                                                                };
    void                        AdvanceAgeAndTime()                                             { AdvanceAgeAndTime(CurrentState().dt()); }                                        // Use committed dt from state



    void                        SetAngularMomentum(double p_AngularMomentum)                    { m_InterimState.SetAngularMomentum(std::max(p_AngularMomentum, 0.0)); }
    void                        SetError(const ERROR p_Error)                                   { m_InterimState.SetError(p_Error); }
//    void                        SetObjectId(const OBJECT_ID p_ObjectId)                         { m_ObjectId = p_ObjectId; }
    
    void                        SetOmega(double p_Omega)                                        { SetAngularMomentum(CalculateMomentOfInertiaAU() * p_Omega); }
    

GNU_CONST inline SN_EVENT AddSNEvent(const SN_EVENT p_SNEvent, const SN_EVENT p_ExistingSNEvent) { return static_cast<SN_EVENT>(static_cast<int>(p_ExistingSNEvent) | static_cast<int>(p_SNEvent)); }
GNU_CONST inline SN_EVENT ClearSNEvent(const SN_EVENT p_SNEvent) { return SN_EVENT::NONE; }
    

    void                        SetEvolutionStatus(const EVOLUTION_STATUS p_EvolutionStatus)    { m_InterimState.SetEvolutionStatus(p_EvolutionStatus); }           // Set evolution status (typically final outcome) for star

    void                        SetDt(const double p_dt)                                        { m_InterimState.SetDt(std::max(0.0, p_dt)); }                      // Set timestep - ensure >= 0.0

    void                        SetMass(const double p_Mass)                                    { m_InterimState.SetMass(std::max(0.0, p_Mass)); }                  // Direct mass setter - non-negative clamp
    void                        SetAge(const double p_Age)                                      { m_InterimState.SetAge(std::max(0.0, p_Age)); }                    // Direct age setter - non-negative clamp
  
    virtual void                UpdateEffectiveZAMSLandR(){}                                                                                                        // Virtual - default is NO-OP
    









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

    virtual double CalculateHAbundanceCore() const                  { return HAbundanceCore(); } // JR FIX THIS: DONE
    virtual double CalculateHAbundanceCoreAtPhaseEnd() const        { return CalculateHAbundanceCore(); } // Same as on phase JR FIX THIS: DONE

    virtual double CalculateHAbundanceSurface() const               { return HAbundanceSurface(); } // JR FIX THIS: DONE
    virtual double CalculateHAbundanceSurfaceAtPhaseEnd() const     { return CalculateHAbundanceSurface(); } // Same as on phase JR FIX THIS: DONE


    virtual double CalculateHeAbundanceCore() const                 { return HeAbundanceCore(); } // JR FIX THIS: DONE
    virtual double CalculateHeAbundanceCoreAtPhaseEnd() const       { return CalculateHeAbundanceCore(); } // Same as on phase JR FIX THIS: DONE

    virtual double CalculateHeAbundanceCoreOut() const              { return HeAbundanceCoreOut(); } // JR FIX THIS: DONE
    virtual double CalculateHeAbundanceCoreOutAtPhaseEnd() const    { return CalculateHeAbundanceCoreOut(); } // Same as on phase JR FIX THIS: DONE

    virtual double CalculateHeAbundanceSurface() const              { return HeAbundanceSurface(); } // JR FIX THIS: DONE
    virtual double CalculateHeAbundanceSurfaceAtPhaseEnd() const    { return CalculateHeAbundanceSurface(); } // Same as on phase JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateAgeAfterMassLoss() const; // JR FIX THIS: DONE
    virtual  double CalculateAgeAfterMassLoss_Hurley2000() const { return Age(); } // JR FIX THIS: DONE


    GNU_CONST double CalculateLifetimeToBAGB_Hurley2000(const double p_tHeI, const double p_tHe) const; // JR FIX THIS: DONE


    GNU_PURE double CalculatePhaseLifetime() const; // JR FIX THIS: DONE
    virtual  double CalculatePhaseLifetime_Hurley2000() const { return 0.0; } // JR FIX THIS: DONE


    GNU_PURE double CalculateTau(const double p_Age, const TimescalesT& p_tScales) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateTau() const { // JR FIX THIS: DONE
        return CalculateTau(Age(), Timescales());
    }
    GNU_PURE double CalculateTau(const double p_Age) const { // JR FIX THIS: DONE
        return CalculateTau(p_Age, Timescales());
    }
    virtual double CalculateTau_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const { return Tau(); } // JR FIX THIS: DONE


    GNU_PURE double CalculateTauAtPhaseEnd(const double p_Age, const TimescalesT& p_tScales) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateTauAtPhaseEnd() const { // JR FIX THIS: DONE
        return CalculateTauAtPhaseEnd(Age(), Timescales());
    }
    GNU_PURE double CalculateTauAtPhaseEnd(const double p_Age) const { // JR FIX THIS: DONE
        return CalculateTauAtPhaseEnd(p_Age, Timescales());
    }
    virtual double CalculateTauAtPhaseEnd_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const { return Tau(); } // JR FIX THIS: DONE


    GNU_PURE TimescalesT CalculateTimescales(const double p_Mass) const { // JR FIX THIS: DONE
        return CalculateTimescales(Metallicity(), p_Mass, Timescales());
    }
    virtual  TimescalesT CalculateTimescales() const { // JR FIX THIS: DONE
        return CalculateTimescales(Metallicity(), Mass(), Timescales());
    };
    virtual  TimescalesT CalculateTimescales_Hurley2000(const double p_Metallicity, const double p_Mass, TimescalesT& p_tScales) const { return Timescales(); } // JR FIX THIS: DONE


    GNU_PURE  double CalculateTimescale_Dynamical() const { // JR FIX THIS: DONE
        return CalculateTimescale_Dynamical_Kalogera1996(Mass(), Radius());
    }
    GNU_CONST double CalculateTimescale_Dynamical_Kalogera1996(const double p_Mass, const double p_Radius) const; // JR FIX THIS: DONE


    double CalculateTimescale_EddyTurnover() const { // JR FIX THIS: DONE
        return CalculateTimescale_EddyTurnover_Hurley2002(Mass(), Radius(), Luminosity(), CoreMass());
    }
    GNU_PURE double CalculateTimescale_EddyTurnover_Hurley2002(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const; // JR FIX THIS: DONE


    GNU_CONST OptDblT CalculateTimescale_MassChange(const STELLAR_TYPE p_StellarType, const STELLAR_TYPE p_StellarTypePrev, const double p_Mass, const double p_MassPrev, const double p_dtPrev) const; // JR FIX THIS: DONE


    GNU_CONST OptDblT CalculateTimescale_RadialExpansion(const STELLAR_TYPE p_StellarType, const STELLAR_TYPE p_StellarTypePrev, const double p_Radius, const double p_RadiusPrev, const double p_dtPrev) const; // JR FIX THIS: DONE


    virtual   double CalculateTimescale_Thermal() const { // JR FIX THIS: DONE
        return CalculateTimescale_Thermal_Kalogera1996(Mass(), Radius(), Luminosity(), CoreMass());
    }
    virtual double CalculateTimescale_Thermal(const double p_Radius) const { // JR FIX THIS: DONE
        return CalculateTimescale_Thermal_Kalogera1996(Mass(), p_Radius, Luminosity(), CoreMass());
    }
    GNU_CONST double CalculateTimescale_Thermal_Kalogera1996(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const; // JR FIX THIS: DONE


    GNU_PURE double CalculateTimestep( // JR FIX THIS: DONE
        const STELLAR_TYPE p_StellarType,
        const double       p_Mass,
        const double       p_Radius,
        const STELLAR_TYPE p_StellarTypePrev,
        const double       p_MassPrev,
        const double       p_RadiuPrev,
        const double       p_dtPrev
    ) const;
    GNU_PURE double CalculateTimestep() const { // JR FIX THIS: DONE 
        return CalculateTimestep(StellarType(), Mass() Radius(), StellarTypePrev(), MassPrev(), RadiusPrev(), dtPrev());
    }


    GNU_PURE double ChooseTimestep() const { // JR FIX THIS: DONE
        return ChooseTimestep(Age, Timescales());
    };
    virtual  double ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const { return NUCLEAR_MINIMUM_TIMESTEP; } // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   ENVELOPE                                   //
    //////////////////////////////////////////////////

    GNU_CONST double CalculateBindingEnergy(const double p_Radius, const double p_CoreMass, const double p_EnvMass, const double p_Lambda) const; // JR FIX THIS: DONE
    inline    double CalculateBindingEnergy(const double p_Lambda) const { // JR FIX THIS: DONE
        return CalculateBindingEnergy(Radius(), CoreMass(), Mass() - CoreMass(), p_Lambda);
    }


    GNU_CONST double CalculateConvectiveEnvelopeBindingEnergy_Hirai2022(const double p_Mass, const double p_Radius, const double p_EnvMass, const double p_Lambda) const; // JR FIX THIS: DONE
    GNU_CONST double CalculateConvectiveEnvelopeBindingEnergy_Hirai2022(const double p_Lambda) const { // JR FIX THIS: DONE
        return CalculateConvectiveEnvelopeBindingEnergy_Hirai2022(Mass(), Radius(), std::get<0>(CalculateConvectiveEnvelopeMass()), p_Lambda);
    }


    GNU_PURE double CalculateConvectiveEnvelopeLambda_Picker2024(const double p_Metallicity, const double p_Mass, const double p_EnvMass, const double p_EnvMassMax) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateConvectiveEnvelopeLambda_Picker2024(const Dbl_DblT p_EnvMass_Max) const { // JR FIX THIS: DONE
        return CalculateConvectiveEnvelopeLambda_Picker2024(Metallicity(), Mass(), std::get<0>(p_EnvMass_Max), std::get<1>(p_EnvMass_Max));
    }


    virtual Dbl_DblT CalculateConvectiveEnvelopeMass() const { return CalculateConvectiveEnvelopeMass(Mass(), Tau()); } // JR FIX THIS: DONE
    virtual Dbl_DblT CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_Tau) const { return {0.0, 0.0}; } // JR FIX THIS: DONE


    virtual double CalculateConvectiveEnvelopeRadialExtent() const { return 0.0; }     // Default for stars with no convective envelope // JR FIX THIS: DONE
    virtual double CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Mass, const double p_Radius) const { return 0.0; } // JR FIX THIS: DONE


    virtual ENVELOPE DetermineEnvelopeType() const { return ENVELOPE::REMNANT; }         // Default is REMNANT - but should never be called
    virtual ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const { return ENVELOPE::REMNANT; }


    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    virtual double CalculateLuminosity() const;
    virtual double CalculateLuminosityAtPhaseEnd() const;
    virtual double CalculateLuminosityAtPhaseEnd_Hurley2000() const { return Luminosity(); } 
    virtual double CalculateLuminosity_Hurley2000() const { return Luminosity(); } 
    virtual double CalculateLuminosity_NoBrcek() const { return CalculateLuminosity(); }



    double CalculateLuminosityAtZAMS(const double p_Metallicity, const double p_MZAMS) const {
        switch (OPTIONS->Mode()) {
            case EVOLUTION_MODE::SSE_HURLEY:
            case EVOLUTION_MODE::BSE_HURLEY: return CalculateLuminosityAtZAMS_Tout1996(p_Metallicity, p_MZAMS);
            default: THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE); return 0.0;
        }
    }
    COMPAS_PURE double              CalculateLuminosityAtZAMS_Tout1996(const double p_Metallicity, const double p_MZAMS) const;
    GNU_CONST   double              CalculateLuminosity_Hurley2000(const double p_CoreMass, const GBParamsT& p_GBParams) const;




    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////

    double CalculateCoreMass() const { return CalculateCoreMass(CoreMass()); };
    double CalculateCoreMass(const double p_CoreMass) const;
    double CalculateCoreMass_Hurley2000(const double p_Luminosity, const GBParamsT& p_GBParams);
    virtual double CalculateCoreMass_Hurley2000() const { return CalculateCoreMass_Hurley2000(Luminosity(), GBParams()); };

    double CalculateCoreMassAtPhaseEnd() const { return CalculateCoreMassAtPhaseEnd(CoreMass()); };
    double CalculateCoreMassAtPhaseEnd(const double p_CoreMass) const;
    virtual double CalculateCoreMassAtPhaseEnd_Hurley2000() { return CalculateCoreMass_Hurley2000(); }


    virtual double CalculateCOCoreMass() const { return CalculateCoreMass(COCoreMass()); };
    double CalculateCOCoreMass(const double p_COCoreMass) const;
    virtual double CalculateCOCoreMass_Hurley2000() const { return COCoreMass(); } 

    virtual double CalculateCOCoreMassAtPhaseEnd() const { return CalculateCoreMassAtPhaseEnd(COCoreMass()); };
    double CalculateCOCoreMassAtPhaseEnd(const double p_COCoreMass) const;
    virtual double CalculateCOCoreMassAtPhaseEnd_Hurley2000() const { return COCoreMass(); } 


    virtual double CalculateConvectiveCoreMass() const { return 0.0; }


    virtual double CalculateMass0() const { return CalculateEffectiveInitialMass(); }
    virtual double CalculateEffectiveInitialMass() const { return CalculateEffectiveInitialMass(MassEffectiveInitial()); };
    double CalculateEffectiveInitialMass(const double p_EffectiveInitialMass) const;
    virtual double CalculateEffectiveInitialMass_Hurley2000() const { return MassEffectiveInitial(); }

    virtual double CalculateEffectiveInitialMassAtPhaseEnd() const { return CalculateEffectiveInitialMassAtPhaseEnd(MassEffectiveInitial()); };
    double CalculateEffectiveInitialMassAtPhaseEnd(const double p_EffectiveInitialMass) const;
    virtual double CalculateEffectiveInitialMassAtPhaseEnd_Hurley2000() const { return MassEffectiveInitial(); }


    virtual double CalculateHeCoreMass() const { return CalculateHeCoreMass(HeCoreMass()); };
    double CalculateHeCoreMass(const double p_HeCoreMass) const;
    virtual double CalculateHeCoreMass_Hurley2000() const { return HeCoreMass(); } 

    virtual double CalculateHeCoreMassAtPhaseEnd() const { return CalculateHeCoreMassAtPhaseEnd(HeCoreMass()); };
    double CalculateHeCoreMassAtPhaseEnd(const double p_HeCoreMass) const;
    virtual double CalculateHeCoreMassAtPhaseEnd_Hurley2000() const { return HeCoreMass(); } 




static double CalculateCoreMass_Hurley2000_Static(const double p_Luminosity, const GBParamsT& p_GBParams);



    //////////////////////////////////////////////////
    //   MASS LOSS / ACCRETION                      //
    //////////////////////////////////////////////////





    // mass loss LBV
    COMPAS_PURE MassLossT CalculateMLRateLBV(const double p_Radius, const double p_Luminosity, const double p_LBVFactor, const LBV_ML_PRESCRIPTION p_MLPrescription) const;

    // mass loss OB
    COMPAS_PURE MassLossT CalculateMLRateOB(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE MassLossT CalculateMLRateOB_Bjorklund2022(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE MassLossT CalculateMLRateOB_Krticka2018(const double p_Metallicity, const double p_Luminosity) const;
    COMPAS_PURE MassLossT CalculateMLRateOB_Vink2001(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE MassLossT CalculateMLRateOB_VinkSander2021(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const;

    // mass loss PRESCRIPTIONS
    GNU_CONST   MassLossT CalculateMLRate_NieuwenhuijzenDeJager1990(const double p_Mass, const double p_Radius, const double p_Luminosity) const;
    GNU_CONST   MassLossT CalculateMLRate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const;
    GNU_CONST   MassLossT CalculateMLRate_KudritzkiReimers1978(const double p_Mass, const double p_Radius, const double p_Luminosity) const;

    // mass loss RSG
    GNU_PURE    MassLossT       CalculateMLRateRSG(
        const double p_Metallicity,
        const double p_Mass,
        const double p_Radius,
        const double p_Luminosity,
        const double p_Temperature,
        const double p_mStart
    ) const;
    GNU_CONST   MassLossT CalculateMLRateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const;
    GNU_CONST   MassLossT CalculateMLRateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const;
    GNU_CONST   MassLossT CalculateMLRateRSG_Kee2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    GNU_CONST   MassLossT CalculateMLRateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const;
    GNU_CONST   MassLossT CalculateMLRateRSG_Yang2023(const double p_Luminosity) const;

    // mass loss VMS
    COMPAS_PURE MassLossT CalculateMLRateVMS(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    GNU_CONST   MassLossT CalculateMLRateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const;
    COMPAS_PURE MassLossT CalculateMLRateVMS_Sabhahit2023(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE MassLossT CalculateMLRateVMS_Vink2011(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const;

    // mass loss WR
    GNU_CONST   MassLossT CalculateMLRateWR_Hurley2000(const double p_Luminosity, const double p_Perturb) const;
    COMPAS_PURE MassLossT CalculateMLRateWR_SanderVink2020(const double p_Metallicity, const double p_Luminosity, const double p_Perturb) const;
    COMPAS_PURE MassLossT CalculateMLRateWR_Shenar2019(const double p_Metallicity, const double p_Luminosity, const double p_Temperature) const;
    GNU_CONST   MassLossT CalculateMLRateWR_ZDependent(const double p_Metallicity, const double p_Luminosity, const double p_Perturb);





    virtual MassLossT CalculateMassLossRate() const {
        //Mu() = calculate perturbation mu
        return CalculateMassLossRate(Mass(), StartMass(), Radius(), Luminosity(), Temperature(), Mu(), HeAbundanceSurface(), OPTIONS->CoolWindMassLossMultiplier())
    };


    virtual MassLossT CalculateMLRate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_Perturb) const;

    virtual MassLossT CalculateMLRate_Merritt2025(
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
    ) const;

    virtual MassLossT CalculateMLRate_Belczynski2010(
        const double              p_Metallicity,
        const double              p_Mass,
        const double              p_Radius,
        const double              p_Luminosity,
        const double              p_Temperature,
        const double              p_HeAbundanceSurface,
        const double              p_Perturb,
        const double              p_CoolWinMultiplier,
        const double              p_LBVFactor,
        const LBV_ML_PRESCRIPTION p_LBVMLPrescription
    ) const;

    virtual double CalculateMLRateThermal() const { return Mass() / CalculateTimescale_Thermal(); }



MassLossT CalculateMassLossValues(double p_Dt, const bool p_UpdateMDot = false);

    
//// Never called    double                      CalculateNuclearMassLossRate()                                                  { return m_Mass / CalculateTimescale_RadialExpansion_DuringMT(); }




    //////////////////////////////////////////////////
    //   MISC. CONSTANTS / PARAMETERS               //
    //////////////////////////////////////////////////

    GBParamsT CalculateGBParams() const { return CalculateGBParams(Mass(), GBParams()); };
    GBParamsT CalculateGBParams(const double p_Mass) const { return CalculateGBParams(p_Mass, GBParams()); }
    GBParamsT CalculateGBParams(const double p_Mass, GBParamsT& p_GBParams) const;

    virtual GBParamsT CalculateGBParams_Hurley2000(const double p_Mass, GBParamsT& p_GBParams) const { return GBParams(); }


    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////

    virtual double CalculateConvectiveCoreRadius(const double p_Radius, const double p_Tau) const { return 0.0; }
    virtual double CalculateConvectiveCoreRadius() const { return 0.0; }

    virtual double CalculateRadius() const;
    virtual double CalculateRadius_Hurley2000() const { return RZAMS(); } 

    virtual double CalculateRadiusAtPhaseEnd() const;
    virtual double CalculateRadiusAtPhaseEnd_Hurley2000() const { return RZAMS(); } 

    virtual double CalculateRadiusOnMassChange(const double p_dM) const { return Radius(); } // Default is NO-OP // JR FIX THIS: DONE

    virtual double CalculateRadiusOnPhase() const { return Radius(); }


    GNU_PURE double CalculateRemnantRadius() const;
    virtual double CalculateRemnantRadius_Hurley2000() const { return Radius(); }




    double CalculateRadiusAtZAMS(const double p_MZAMS) const {
        switch (OPTIONS->Mode()) {
            case EVOLUTION_MODE::SSE_HURLEY:
            case EVOLUTION_MODE::BSE_HURLEY: return CalculateRadiusAtZAMS_Tout1996(p_MZAMS);
            default: THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE); return 0.0;
        }
    }


//// Never called   double                      CalculateRadialChange() const                                                   { return (utils::Compare(m_RadiusPrev,0) <= 0)? 0 : std::abs(m_Radius - m_RadiusPrev) / m_RadiusPrev; } // Return fractional radial change (if previous radius is negative or zero, return 0 to avoid NaN




    //////////////////////////////////////////////////
    //   ROTATION                                   //
    //////////////////////////////////////////////////

    // WHAT SHOULD THESE DEFAULT TO IF NOT MS START?  CHECK ALL OTHER STs
    virtual double CalculateMomentOfInertia_Hurley2000() const { const double mass = CurrentState().Mass(); const double radius = CurrentState().Radius(); return (0.1 * mass * radius * radius); } // Defaults to MS. k2 = 0.1 as defined in Hurley et al. 2000, after eq 109
    virtual double CalculateMomentOfInertiaAU() const { return CalculateMomentOfInertia_Hurley2000() * RSOL_TO_AU * RSOL_TO_AU; }




COMPAS_PURE static double       CalculateAngularFrequency_CHE_Static(const double p_MZAMS);

// COMPAS_PURE double              CalculateRotationalVelocityOStar_Ramirez2013() const; // MOVED to utils::CalculateRotationalVelocityOStar_Ramirez2013 (free function)





    //////////////////////////////////////////////////
    //   SUPERNOVAE                                 //
    //////////////////////////////////////////////////

    COMPAS_PURE StellarSNDetailsT CalculateSNkickMagnitude(
        const STELLAR_TYPE       p_RemnantType,
        const double             p_Mass,
        const double             p_EjectaMass,
        const double             p_RemnantMass,
        const StellarSNDetailsT& p_SNdetails
    ) const;

    COMPAS_PURE double              DrawSNkickMagnitude(
        const double   p_COCoreMass,
        const double   p_EjectaMass,
        const double   p_RemnantMass,
        const double   p_Sigma,
        const double   p_Rand,
        const SN_EVENT p_SNevent
    ) const;

    GNU_CONST   double              DrawSNkickMagnitude_BrayEldridge2018(const double p_EjectaMass, const double p_RemnantMass, const double p_Alpha, const double p_Beta) const;
    GNU_CONST   double              DrawSNkickMagnitude_Flat(const double p_KickMax, const double p_Rand) const;
    GNU_CONST   double              DrawSNkickMagnitude_Maxwellian(const double p_Sigma, const double p_Rand) const;
    GNU_CONST   double              DrawSNkickMagnitude_Muller2016(const double p_COCoreMass) const;
    COMPAS_PURE double              DrawSNkickMagnitude_MullerMandel2020(const double p_COCoreMass, const double p_RemnantMass, const double p_Rand, const SN_EVENT p_SNevent) const;






    //////////////////////////////////////////////////
    //   TEMPERATURE                                //
    //////////////////////////////////////////////////

    virtual double CalculateTemperature() const { return CalculateTemperature(Luminosity(), Radius()); };
    COMPAS_PURE double CalculateTemperature(const double p_Luminosity, const double p_Radius) const;

























    
    
    














    
    



        








 


        


    // AdvanceOneTimestep is the state-machine transition function (see BaseStar.cpp).
    // It supersedes EvolveOneTimestep. It reads the frozen current state, builds the
    // next state in m_InterimState (two phases: predictors, then attributes), commits
    // it to the history via Push(), and returns the next stellar type for the caller
    // (the Star wrapper) to switch to if it differs from the current type.
    //
    // Optional inputs (Fork B, the unified per-step entry point):
    //   p_Dt          - timestep override (Myr). Any non-negative value is used as-is
    //                   (caller has applied any TimestepMultiplier/quantisation); zero
    //                   is legitimate ("no time advance", used by mass-only BSE updates).
    //                   The default sentinel -1.0 (or any negative value) means "compute
    //                   internally via CalculateTimestep()" (the old behaviour).
    //   p_DeltaMass   - externally-supplied mass change (Msol), e.g. BSE mass transfer. When
    //                   non-zero, wind mass loss is skipped this step and the supplied value is
    //                   used directly. Default 0.0 means "no external delta, compute wind loss".
    //   p_DeltaMass0  - externally-supplied effective-initial-mass change. When non-zero, used
    //                   directly; default 0.0 means "recompute via Hurley".
    STELLAR_TYPE                AdvanceOneTimestep(const double p_Dt         = -1.0,
                                                   const double p_DeltaMass  =  0.0,
                                                   const double p_DeltaMass0 =  0.0);


    void HaltWinds() { m_InterimState.SetMdot(0.0); }   // Disable wind mass loss in current time step


    void ResetEnvelopeExpulsationByPulsations() { m_InterimState.ResetEnvelopeExpelledByPulsations(); }





    virtual STELLAR_TYPE ResolveEnvelopeLoss(bool p_Force = false)                                        { return CurrentState().StellarType(); }

    virtual STELLAR_TYPE        ResolveMassLoss(const double p_dt) { return StellarType(); }
    virtual STELLAR_TYPE        ResolveMassLossHurley(const double p_dt);
    virtual void                ResolveShellChange(const double p_AccretedMass) { }                                                                                                 // Default does nothing, use inheritance for WDs.
    virtual STELLAR_TYPE        ResolveSupernova(bool p_HasMTdonorHistory = false)                              { return CurrentState().StellarType(); }                            // Default is NO-OP.  
       
    void                        SetStellarTypePrev(const STELLAR_TYPE p_StellarTypePrev)                        { }
    
    virtual bool                ShouldEnvelopeBeExpelledByPulsations() const                                    { return false; }                                                   // Default is that there is no envelope expulsion by pulsations

    virtual void                SpinDownIsolatedPulsar(const double p_Stepsize) { }                                                                                                 // Default is NO-OP
   
    
    
    virtual void                UpdateMagneticFieldAndSpin(const bool p_CommonEnvelope, const bool p_RecyclesNS, const double p_Stepsize, const double p_MassGainPerTimeStep, const double p_Epsilon) { } // Default is NO-OP
    virtual void                UpdateMainSequenceCoreMass(const double p_Dt, const double p_TotalMassLossRate) { }                                                                 // Set core mass for Main Sequence stars; default is NO-OP
    virtual void                UpdateTotalMassLossRate(const double p_MassLossRate)                            { (void)p_MassLossRate; }
    
    void                        UpdateSNAnomalies(const double p_Eccentricity);


    // printing functions
    bool PrintDetailedOutput(const int p_Id, const SSE_DETAILED_RECORD_TYPE p_RecordType) const { 
        return OPTIONS->DetailedOutput() ? LOGGING->LogSSEDetailedOutput(this, p_Id, p_RecordType) : true;                                                                          // Write record to SSE Detailed Output log file
    }

    bool PrintPulsarEvolutionParameters(const SSE_PULSAR_RECORD_TYPE p_RecordType) const {
        return OPTIONS->EvolvePulsars() ? LOGGING->LogSSEPulsarEvolutionParameters(this, p_RecordType) : true;
    }

    bool PrintSupernovaDetails(const SSE_SN_RECORD_TYPE p_RecordType = SSE_SN_RECORD_TYPE::DEFAULT) const {
        return LOGGING->LogSSESupernovaDetails(this, p_RecordType);                                                                                                                 // Write record to SSE Supernovae log file
    }

    bool PrintSwitchLog() const { 
        return OPTIONS->SwitchLog() ? (LOGGING->ObjectSwitchingPersistence() == OBJECT_PERSISTENCE::PERMANENT ? LOGGING->LogSSESwitchLog(this) : true) : true;                      // Write record to SSE Switchlog log file
    }

    bool PrintSystemSnapshotLog(const SSE_SYSTEM_SNAPSHOT_RECORD_TYPE p_RecordType = SSE_SYSTEM_SNAPSHOT_RECORD_TYPE::DEFAULT) const {
        return LOGGING->LogSSESystemSnapshotLog(this, p_RecordType);                                                                                                                // Write record to SSE System Parameters file
    }

    bool PrintSystemParameters(const SSE_SYSPARMS_RECORD_TYPE p_RecordType = SSE_SYSPARMS_RECORD_TYPE::DEFAULT) const {
        return LOGGING->LogSSESystemParameters(this, p_RecordType);                                                                                                                 // Write record to SSE System Parameters file
    }





















// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// removed in PR 1445    
//GNU_CONST double inline BaseStar::CalculateMLfractionWR(const double p_HeAbundanceSurface) const;




    void                CalculateLCoefficients(const double p_LogMetallicityXi, DblVectorT &p_LCoefficients) const;






COMPAS_PURE double CalculateLuminosityAtBAGB_Hurley2000(const double p_Metallicity, const double p_Mass) const;







    
    void                CalculateMassCutoffs(DblVectorT& p_MassCutoffs) const;





























   
    
    
    
    
    



    double              CalculateMassLossRateWolfRayetTemperatureCorrectionSander2023(const double p_Mdot) const;



// JR : DON'T FORGET "OVERRIDE" ON DERIVED FUNCTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



virtual double CalculateMomentOfInertia() const                                                                     { return CalculateMomentOfInertia(Mass(), Radius()); }
virtual double CalculateMomentOfInertia(const double p_Mass, const double p_Radius) const                           { return 0.1 * p_Mass * p_Radius * p_Radius; } // Defaults to MS. k2 = 0.1 as defined in Hurley et al. 2000, after eq 109
    // Pre-pass note: duplicate CalculateMomentOfInertiaAU declaration removed here -- canonical at line 655.


    GNU_CONST double CalculateEddingtonLuminosity(const double p_Mass, const double p_HeAbundanceSurface) const;













// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ********************************************************* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>























virtual inline double CalculateHurleyPerturbationMu() const { return CurrentState().Mu(); }
virtual inline double CalculateHurleyPerturbationMuAtPhaseEnd() const { return CalculateHurleyPerturbationMu(); }

GNU_CONST inline double CalculateHurleyPerturbationB(const double p_Mass) const { return 0.002 * std::max(1.0, (2.5 / p_Mass)); } // Hurley et al. 2000, eq 103

GNU_CONST inline double CalculateHurleyPerturbationC(double p_Mass) const { return 0.006 * std::max(1.0, (2.5 / p_Mass)); } // Hurley et al. 2000, eq 104

GNU_CONST inline double CalculatePerturbationQ(const double p_Radius, const double p_Rc) const { return log(p_Radius / p_Rc); } // Hurley et al. 2000, eq 105 (really is natural log)

GNU_CONST double CalculateHurleyPerturbationR(const double p_Mass, const double p_Radius, const double p_Mu, const double p_Rc) const;

GNU_CONST double CalculateHurleyPerturbationS(const double p_Mass, const double p_Mu) const;











    virtual std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                                              { return CalculateRadiusAndStellarTypeOnPhase(StellarType()); }
    virtual std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const STELLAR_TYPE p_StellarType) const              { return std::make_tuple(CalculateRadiusOnPhase(), p_StellarType); }

    void                CalculateRCoefficients(const double p_LogMetallicityXi, DblVectorT& p_RCoefficients) const;


GNU_PURE double CalculateAngularFrequencyAtZAMS(const double p_MZAMS, const double p_RZAMS) const;
GNU_PURE double CalculateAngularFrequencyAtZAMS_Hurley2000(const double p_MZAMS, const double p_RZAMS) const;



  






























///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<








    virtual double      CalculateTemperatureAtPhaseEnd() const                                                          { return CalculateTemperatureAtPhaseEnd(Luminosity(), Radius()); }
    virtual double      CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const          { return CalculateTemperatureOnPhase(p_Luminosity, p_Radius); }             // Same as on phase
    double              CalculateTemperatureKelvinOnPhase(const double p_Luminosity, const double p_Radius) const;


////    virtual double      CalculateTemperatureOnPhase() const                                                             { return CalculateTemperatureOnPhase(m_Luminosity, m_Radius); }



    virtual double      CalculateTemperatureOnPhase(const double p_Luminosity, const double p_Radius) const;







    STELLAR_TYPE        EvolveOnPhase(const double p_DeltaTime);

    virtual STELLAR_TYPE EvolveToNextPhase()                                                                            { return EvolveToNextPhase(StellarType()); }
    virtual STELLAR_TYPE EvolveToNextPhase(const STELLAR_TYPE p_StellarType)                                            { return p_StellarType; }




    virtual bool        IsEndOfPhase() const                                                                            { return false; }

    /*
     * Perturb Luminosity and Radius
     *
     * See Hurley et al. 2000, section 6.3
     *
     * The default is no perturbation - this function does nothing and is called
     * only if the stellar class doesn't define its own perturbation function.
     * See the stellar class perturbation functions for perturbation details specific
     * to the stellar class.
     *
     * Perturbation is disabled by default when DEBUG is enabled - except when
     * DEBUG_PERTURB is defined (see below).  The stellar class perturbation
     * functions are defined away if DEBUG is defined - so this generic Star
     * function is called (and does nothing).
     *
     * If DEBUG_PERTURB is defined then perturbation is not disabled while debbuging.
     * To enable perturbation while DEBUG is enabled, define DEBUG_PERTURB.
     */
    virtual void        PerturbLuminosityAndRadius() { }                                                                                                                                            // NO-OP
    virtual void        PerturbLuminosityAndRadiusAtPhaseEnd()                                                          { PerturbLuminosityAndRadiusOnPhase(); }                                    // Same as on phase
    virtual void        PerturbLuminosityAndRadiusOnPhase()                                                             { PerturbLuminosityAndRadius(); }

    STELLAR_TYPE        ResolveEndOfPhase();
    virtual void        ResolveHeliumFlash() { }
    virtual STELLAR_TYPE ResolveSkippedPhase()                                                                          { return EvolveToNextPhase(); }                                             // Default is evolve to next phase

    double              ReweightSupernovaKickByMass(const double p_vK,
                                                    const double p_FallbackFraction,
                                                    const double p_BlackHoleMass)                                       { return p_vK; }                                                            // Default is not to re-weight, except for black holes where the --black-hole-kicks-prescription option is relevant

    
    virtual void        SetSNHydrogenContent()                                                                          { m_InterimState.SetSupernovaDetailsIsHydrogenPoor(false); }                              // Default is false

    void                SetSNCurrentEvent(const SN_EVENT p_SNEvent)                                                     { m_InterimState.AddSNCurrentEventBit(p_SNEvent); }
    void                SetSNPastEvent   (const SN_EVENT p_SNEvent)                                                     { m_InterimState.AddSNPastEventBit   (p_SNEvent); }

    bool         ShouldBeMasslessRemnant() const                                                                 { return (CurrentState().Mass() <= 0.0 || CurrentState().StellarType() == STELLAR_TYPE::MASSLESS_REMNANT); }
    virtual bool        ShouldEvolveOnPhase() const                                                                     { return true; }
    virtual bool        ShouldSkipPhase() const                                                                         { return false; }                                                           // Default is false



    virtual std::unique_ptr<BaseStar> CloneAs(STELLAR_TYPE p_StellarType, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) const { return nullptr; }
    virtual std::unique_ptr<BaseStar> MakeStar(STELLAR_TYPE p_StellarType, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) const { return nullptr; }







    /* OStarRotationVelocityFunctor - MOVED to namespace utils in utils.h (was nested template struct here) */


















//// <<<<<<<<<<<<<<<<<<<<<<<<<<<<< constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    virtual double      CalculateMTRejuvenationFactor()                                                       { return 1.0; }

    virtual MT_CASE             DetermineMassTransferTypeAsDonor() const                                        { return MT_CASE::OTHER; } // Not A, B, C, or NONE


    void                        ApplyMassTransferRejuvenationFactor()                                           { m_InterimState.SetAge(CurrentState().Age() * CalculateMTRejuvenationFactor()); } // Apply age rejuvenation factor

    virtual void                UpdateComponentVelocity(const Vector3d p_NewVelocity)           { (void)p_NewVelocity; }

    // JR FIX THIS - placeholders - populate
    virtual void                UpdateInitialMass()                                              {  }
    virtual void                ClearCurrentSNEvent()                                            {  }
    virtual double              CalculateEddyTurnoverTimescale() const                           { return 1.0; }
    virtual double              CalculateRadialExpansionTimescale() const                        { return 1.0; }
    virtual double              CalculateTAMSCoreMass() const                                    { return 0.0; }
    virtual double              CalculateConvectiveEnvelopeBindingEnergy(const double p_Lambda) const  { return CalculateConvectiveEnvelopeBindingEnergy_Hirai2022(p_Lambda); }
    virtual double              CalculateNuclearMassLossRate() const                             { return 0.0; }
    virtual double              InitialMainSequenceCoreMass() const                              { return CNOProcessedCoreMass(); } 
    bool                HasMTdonorHistory() const                                                { return CurrentState().HasMTdonorHistory(); }
    virtual StrT        GetMassTransferDonorHistoryString() const                                { return CurrentState().MTdonorHistoryString(); }

    // JR FIX THIS ...
    // BSE-only writers used by BinaryConstituentStar::UpdateMassTransferDonorHistory to
    // mirror the binary-side state into per-star reach.  SSE never calls these.
    void                SetHasMTdonorHistory(const bool p_HasHistory)                            { m_InterimState.SetHasMTdonorHistory(p_HasHistory); }
    void                SetMTdonorHistoryString(const StrT& p_Str)                               { m_InterimState.SetMTdonorHistoryString(p_Str); }
  

    
virtual double CalculateCriticalMassRatio(const double p_Mass,
                                          const double p_Radius,
                                          const double p_CoreMass,
                                          const bool   p_AccretorIsDegenerate,
                                          const double p_MTefficiency = 0.0) const { return 0.0; } // JR FIX THIS - placeholder
    GNU_PURE double CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate, const double p_MTefficiency) const;

virtual double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const { return 0.0; } // Default is 0.0

virtual double CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency = 0.0) const { return 0.0; } 

virtual double CalculateCriticalMassRatio_Ge2020(const double p_MTefficiency) const {
    return CalculateCriticalMassRatio_Ge2020_Interpolate(Mass(), Radius(), p_MTefficiency);
}

virtual double CalculateCriticalMassRatio_Hurley2002() const { return 0.0; } // Default is 0.0




inline COMPAS_PURE double CalculateZetaAdiabatic_Hurley2002(const double p_Mass, const double p_CoreMass) const;

GNU_CONST double CalculateZetaAdiabatic_Soberman1997(const double p_Mass, const double p_CoreMass) const;




    
    virtual Dbl_Dbl_Dbl_DblT     CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const;
    virtual Dbl_Dbl_Dbl_DblT     CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const ;
    virtual Dbl_Dbl_Dbl_DblT     CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const;
    

    virtual double              CalculateCELambda_Dewi() const                                                     { return 1.0; }           
    double                      CalculateCELambda_Kruckow(const double p_Radius, const double p_Alpha) const;

                                              // Default for stellar types with no LamdaDewi definitions - 1.0 is benign

GNU_CONST double CalculateCELambdaKruckow(const double p_Radius, const double p_Alpha) const { (void)p_Radius; (void)p_Alpha; return 0.0; }  // JR FIX THIS - placeholder
    double CalculateCELambdaKruckow() const { return CalculateCELambdaKruckow(Radius(), OPTIONS->CommonEnvelopeSlopeKruckow()); }
    
    virtual double              CalculateLambdaLoveridge(const double p_EnvMass, const bool p_IsMassLoss = false) const { return 1.0; } // Default for non giant branch stars - 1.0 is benign
    double                      CalculateLambdaLoveridge() const                                                { return CalculateLambdaLoveridge(Mass() - CoreMass(), false); }


//    double                      CalculateCELambda_Nanjing() const;
    


    Dbl_DblT                     CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate); 
    virtual Dbl_DblT             CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                            const double p_AccretorMassRate,
                                                            const bool   p_IsHeRich)                            { return CalculateMassAcceptanceRate(p_DonorMassRate, p_AccretorMassRate); } // Ignore the He content for non-WDs
    double                      CalculateMassAccretedForCO(const double p_Mass, const double p_CompanionMass, const double p_CompanionRadius, const double p_CompanionEnvelope) const;


COMPAS_PURE double CalculateCELambda_Nanjing(const double p_Mass, const double p_Radius, const double p_CoreMass) const;

virtual double CalculateCELambda_Nanjing_StarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const { return 1.0; } // Default for stellar types with no LamdaNanjing definitions - 1.0 is benign




virtual double CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
                                                            const double             p_Radius,
                                                            const double             p_CoreMass,
                                                            const SizeT              p_MassIndex,
                                                            const STELLAR_POPULATION p_StellarPop) const { return 1.0; } // Default for stellar types with no LamdaNanjing definitions - 1.0 is benign



    double                      CalculateThermalMassAcceptanceRate(const double p_Radius);
    double                      CalculateThermalMassAcceptanceRate()                                            { return CalculateThermalMassAcceptanceRate(Radius()); }


    
    double                      CalculateZetaAdiabatic() const;
    virtual double              CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription) const   { return 0.0; } // Use inheritance hierarchy
    virtual double              CalculateZetaEquilibrium()                                                      { return 0.0; }



    virtual ACCRETION_REGIME    DetermineAccretionRegime(const double p_DonorThermalMassLossRate, 
                                                         const bool p_HeRich)                                   { return ACCRETION_REGIME::ZERO; } // Placeholder, use inheritance for WDs

    virtual bool                IsMassAboveChandrasekhar() const                                                { return false; } // JR FIX THIS - placeholder
    virtual STELLAR_TYPE        ResolveAIC()                                                                    { return ResolveAIC(StellarType()); } // JR FIX THIS - placeholder
    virtual STELLAR_TYPE        ResolveAIC(const STELLAR_TYPE p_StellarType)                                    { return p_StellarType; } // JR FIX THIS - placeholder
    virtual STELLAR_TYPE        ResolveSNIa()                                                                   { return ResolveSNIa(StellarType()); } // JR FIX THIS - placeholder
    virtual STELLAR_TYPE        ResolveSNIa(const STELLAR_TYPE p_StellarType)                                   { return p_StellarType; } // JR FIX THIS - placeholder
    virtual STELLAR_TYPE        ResolveHeSD()                                                                   { return ResolveHeSD(StellarType()); } // JR FIX THIS - placeholder
    virtual STELLAR_TYPE        ResolveHeSD(const STELLAR_TYPE p_StellarType)                                   { return p_StellarType; } // JR FIX THIS - placeholder

    virtual double              HeShell()      const                                                            { return 0.0; } // JR FIX THIS - placeholder
    virtual double              HShell()       const                                                            { return 0.0; } // JR FIX THIS - placeholder
    virtual double              L0Ritter()     const                                                            { return 0.0; } // JR FIX THIS - placeholder
    virtual double              LambdaRitter() const                                                            { return 0.0; } // JR FIX THIS - placeholder
    virtual double              XRitter()      const                                                            { return 0.0; } // JR FIX THIS - placeholder
    virtual void                SetHeShell(const double p_Value)                                                {  } // JR FIX THIS - placeholder
    virtual void                SetHShell(const double p_Value)                                                 {  } // JR FIX THIS - placeholder
    virtual double              CalculateEtaH(const double p_MassIntakeRate)                                    { return 0.0; } // JR FIX THIS - placeholder
    virtual double              CalculateEtaHe(const double p_MassIntakeRate)                                   { return 0.0; } // JR FIX THIS - placeholder



    void                        ResolveAccretion(const double p_AccretionMass)                                  { m_InterimState.SetMass(std::max(0.0, Mass() + p_AccretionMass)); } // Handles donation and accretion - won't let mass go negative
    virtual void                ResolveAccretionRegime(const ACCRETION_REGIME p_Regime, const double p_DonorThermalMassLossRate) { }                                                // Default does nothing, only works for WDs.
    virtual double              ResolveCommonEnvelopeAccretion(const double p_FinalMass,
                                                               const double p_CompanionMass     = 0.0,
                                                               const double p_CompanionRadius   = 0.0,
                                                               const double p_CompanionEnvelope = 0.0)          { return p_FinalMass - Mass(); }                                    // Overwritten in NS.h; for now, no accretion on stars other than compact objects during CE



    virtual void                UpdateAfterMerger(double p_Mass, double p_HydrogenMass) { }                                                                                         // Default is NO-OP
    virtual void                UpdateAgeAfterMassLoss() { }  //  // Default is NO-OP




GNU_CONST double  CalculateCELambda_Nanjing_MassInterpolated(const double p_Mass, const double p_Radius, const double p_CoreMass, const STELLAR_POPULATION stellarPop) const;



GNU_PURE  double CalculateCELambda_Nanjing_MassInterpolated(const double p_Mass, const STELLAR_POPULATION p_StellarPop) const;
GNU_PURE  double CalculateCELambda_Nanjing_ZInterpolated(const double p_Mass, const double p_Radius, const double p_CoreMass, const SizeT p_MassIndex) const;
    GNU_PURE  double CalculateCELambda_Nanjing_MassAndZInterpolated(const double p_Mass, const double p_Radius, const double p_CoreMass) const;


};  // end of class BaseStar


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
 * CalculateAgeAfterMassLoss
 *
 * @brief
 * Calculate the age of a star after mass loss.
 * 
 * Calls relevant function based on the evolutionary mode given in program options.
 *
 *
 * double CalculateAgeAfterMassLoss() const
 *
 * @return                                      Age of the star after mass loss (Myr)
 */
inline double BaseStar::CalculateAgeAfterMassLoss() const {

    double age = 0.0;                                               // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            age = CalculateAgeAfterMassLoss_Hurley2000();
            break;

        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }

    return age;
}


/*
 * CalculateLifetimeToBAGB_Hurley2000
 *
 * @brief
 * Calculate the lifetime to the Base of the Asymptotic Giant Branch, BAGB,
 * per Hurley et al. 2000, just before eq 69
 * 
 * tBAGB = tHeI + tHe
 *
 *
 * double CalculateLifetimeToBAGB_Hurley2000(const double p_tHeI, const double p_tHe)
 *
 * @param       p_tHeI                          Time to helium ignition, tHeI (Myr)
 * @param       p_tHe                           Time to helium burning, tHe (Myr)
 * @return                                      Lifetime to the BAGB (Myr)
 */
inline double BaseStar::CalculateLifetimeToBAGB_Hurley2000(const double p_tHeI, const double p_tHe) const {
    return p_tHeI + p_tHe;
}


/*
 * CalculatePhaseLifetime
 *
 * @brief
 * Calculate the lifetime of the current phase.
 * 
 * Calls relevant function based on the evolutionary mode given in program options.
 *
 *
 * double CalculatePhaseLifetime() const
 *
 * @return                                      Lifetime of the phase (Myr)
 */
inline double BaseStar::CalculatePhaseLifetime() const {

    double lifetime = 0.0;                                          // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            lifetime = CalculatePhaseLifetime_Hurley2000();
            break;

        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }

    return lifetime;
}


/*
 * CalculateTau
 *
 * @brief
 * Calculate the phase-relative age of the star (fractional age on the current
 * evolutionary phase).
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateTau(const double p_Age, const TimescalesT& p_tScales) const
 * 
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Phase-relative age of the star, [0, 1]
 */
inline double BaseStar::CalculateTau(const double p_Age, const TimescalesT& p_tScales) const { 

    double tau = 0.0;                                               // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            tau = CalculateTau_Hurley2000(p_Age, p_tScales);
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return tau;
}


/*
 * CalculateTauAtPhaseEnd
 *
 * @brief
 * Calculate the phase-relative age of the star (fractional age on the current
 * evolutionary phase) at the end of the phase.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateTau(const double p_Age, const TimescalesT& p_tScales) const
 * 
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Phase-relative age of the star, [0, 1]
 */
inline double BaseStar::CalculateTauAtPhaseEnd(const double p_Age, const TimescalesT& p_tScales) const { 

    double tau = 0.0;                                               // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            tau = CalculateTauAtPhaseEnd_Hurley2000(p_Age, p_tScales);
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return tau;
}


/*
 * CalculateTimescales
 *
 * @brief
 * Calculate timescales.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 * Since timescales depend on a star's mass, they need to be calculated whenever the mass of the
 * star changes (probably every timestep).
 *
 * p_tScales is copied once and returned by value.  NRVO constructs the returned array directly
 * in the caller's return slot, so no copy on the return.
 *
 * 
 * TimescalesT CalculateTimescales(const double p_Metallicity, const double p_Mass, TimescalesT& p_tScales) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
inline TimescalesT BaseStar::CalculateTimescales(const double p_Metallicity, const double p_Mass, TimescalesT& p_tScales) const { 

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            p_tScales = CalculateTimescales_Hurley2000(p_Metallicity, p_Mass, p_tScales);
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }
    
    return p_tScales;
}


/*
 * CalculateTimescale_Dynamical_Kalogera1996
 *
 * @brief
 * Calculate dynamical timescale, per Kalogera & Webbink 1996, eq 1
 *
 *
 * double CalculateTimescale_Dynamical_Kalogera1996(const double p_Mass, const double p_Radius) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @return                                      Dynamical timescale (Myr)
 */
inline double BaseStar::CalculateTimescale_Dynamical_Kalogera1996(const double p_Mass, const double p_Radius) const {
    return 5.0 * 1.0E-5 * p_Radius * std::sqrt(p_Radius) * YEAR_TO_MYR / std::sqrt(p_Mass);
}


/*
 * CalculateTimescale_EddyTurnover_Hurley2002
 *
 * @brief
 * Calculate the eddy turnover timescale, per Hurley et al. 2002, sec 2.3,
 * particularly eq 31 of subsec 2.3.1
 *
 *
 * double CalculateTimescale_EddyTurnover_Hurley2002(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Eddy turnover timescale (yr)   // *Ilya* why is this in yr, not Myr?  Seems inconsistent.
 */
inline double BaseStar::CalculateTimescale_EddyTurnover_Hurley2002(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const {
	const double rEnv = CalculateConvectiveEnvelopeRadialExtent();
    const double mEnv = std::get<0>(CalculateConvectiveEnvelopeMass(p_Mass, p_CoreMass));
    return 0.4311 * std::cbrt((mEnv * rEnv * (p_Radius - (0.5 * rEnv))) / (3.0 * p_Luminosity));
}


/*
 * CalculateTimescale_MassChange
 *
 * @brief
 * Calculate the mass change timescale.
 *
 * If the mass change timescale can't be calculated (because, for example,
 * the stellar type has changed, or the mass has not changed from the previous
 * timestep), the return value will be std::nollopt.
 * 
 * 
 * OptDblT CalculateTimescale_MassChange(const STELLAR_TYPE p_StellarType, const STELLAR_TYPE p_StellarTypePrev, const double p_Mass, const double p_MassPrev, const double p_dtPrev) const
 *
 * @param       p_StellarType                   Stellar type of the star
 * @param       p_StellarTypePrev               Previous stellar type of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_MassPrev                      Previous mass of the star (Msol)
 * @param       p_dtPrev                        Previous timestep (Myr)
 * @return                                      Mass change timescale (Myr)
 *                                              Wrapped in std::optional<> and returns std::nullopt if the mass change
 *                                              timescale can't be calculated (because, for example, the stellar type has
 *                                              changed, or the mass has not changed)
 */
inline OptDblT BaseStar::CalculateTimescale_MassChange(const STELLAR_TYPE p_StellarType, const STELLAR_TYPE p_StellarTypePrev, const double p_Mass, const double p_MassPrev, const double p_dtPrev) const {
    return ((p_StellarTypePrev == p_StellarType) && (CDOUBLE(p_MassPrev) != p_Mass))
                ? (p_dtPrev * p_MassPrev) / std::fabs(p_Mass - p_MassPrev)
                : std::nullopt;
}


/*
 * CalculateTimescale_RadialExpansion
 *
 * @brief
 * Calculate the radial expansion timescale.
 *
 * If the radial expansion timescale can't be calculated (because, for example,
 * the stellar type has changed, or the radius has not changed, from the previous
 * timestep), the return value will be std::nullopt.
 * 
 *
 * OptDblT CalculateTimescale_RadialExpansion(const STELLAR_TYPE p_StellarType, const STELLAR_TYPE p_StellarTypePrev, const double p_Radius, const double p_RadiusPrev, const double p_dtPrev) const
 *
 * @param       p_StellarType                   Stellar type of the star
 * @param       p_StellarTypePrev               Previous stellar type of the star
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_RadiusPrev                    Previous radius of the star (Rsol)
 * @param       p_dtPrev                        Previous timestep (Myr)
 * @return                                      Radial expansion timescale (Myr)
 *                                              Wrapped in std::optional<> and returns std::nullopt if the radial expansion
 *                                              timescale can't be calculated (because, for example, the stellar type has
 *                                              changed, or the radius has not changed)
 */
inline OptDblT BaseStar::CalculateTimescale_RadialExpansion(const STELLAR_TYPE p_StellarType, const STELLAR_TYPE p_StellarTypePrev, const double p_Radius, const double p_RadiusPrev, const double p_dtPrev) const {
    return ((p_StellarTypePrev == p_StellarType) && (CDOUBLE(p_RadiusPrev) != p_Radius))
                ? (p_dtPrev * p_RadiusPrev) / std::fabs(p_Radius - p_RadiusPrev)
                : std::nullopt;
}


/*
 * CalculateTimescale_Thermal_Kalogera1996
 *
 * @brief
 * Calculate thermal timescale
 *
 * Pre-factor from Kalogera & Webbink 1996 (https://arxiv.org/abs/astro-ph/9508072), equation 2, 
 * combined with p_Mass * p_EnvMass case from equation 61 from https://arxiv.org/abs/astro-ph/0201220 for k in {2,3,4,5,6,8,9}
 * [note that equation 61 of BSE (https://arxiv.org/abs/astro-ph/0201220) approximates this with a value a factor of 3 smaller]
 * 
 * 
 * double CalculateTimescale_Thermal_Kalogera1996(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Thermal timescale (Myr)
 */
inline double BaseStar::CalculateTimescale_Thermal_Kalogera1996(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const {   
    // G*Msol^2/(Lsol*Rsol) ~ 31.4 Myr (~ 30 Myr in Kalogera & Webbink)
    return 31.4 * p_Mass * ((CDOUBLE(p_Mass) == p_CoreMass) ? p_Mass : p_Mass - p_CoreMass) / (p_Radius * p_Luminosity);
}


/*
 * CalculateTimestep
 *
 * @brief
 * Calculate next timestep for stellar evolution, based on stellar attributes.
 *
 *
 * double CalculateTimestep(
 *     const STELLAR_TYPE p_StellarType,
 *     const double       p_Mass,
 *     const double       p_Radius,
 *     const STELLAR_TYPE p_StellarTypePrev,
 *     const double       p_MassPrev,
 *     const double       p_RadiusPrev,
 *     const double       p_dtPrev
 * ) const
 * 
 * @param       p_StellarType                   Stellar type of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_StellarTypePrev               Previous stellar type of the star
 * @param       p_MassPrev                      Previous mass of the star (Msol)
 * @param       p_RadiusPrev                    Previous radius of the star (Rsol)
 * @param       p_dtPrev                        Previous timestep (Myr)
 * @return                                      Timestep (Myr)
 */
inline double BaseStar::CalculateTimestep(
    const STELLAR_TYPE p_StellarType,
    const double       p_Mass,
    const double       p_Radius,
    const STELLAR_TYPE p_StellarTypePrev,
    const double       p_MassPrev,
    const double       p_RadiusPrev,
    const double       p_dtPrev
) const {

    double dt = 0.0;                                            // Default return value

    if (OPTIONS->MassChangeFraction() > 0.0) {                  // User specified mass change fraction?
                                                                // Yes
        // Calculate the mass change timescale.
        // No value means it could not be computed (e.g., just after stellar type change)
        OptDblT tsMassChange = CalculateTimescale_MassChange(p_StellarType, p_StellarTypePrev, p_Mass, p_MassPrev, p_dtPrev);
        if (tsMassChange.has_value()) dt = p_MassChangeFraction * tsMassChange;
    }

    if (OPTIONS->RadialChangeFraction() > 0.0) {                // User specified radial change fraction?
                                                                // Yes
        // Calculate the radial expansion timescale.
        // No value means it could not be computed (e.g., just after stellar type change)
        OptDblT tsRadialExp = CalculateTimescale_RadialExpansion(p_StellarType, p_StellarTypePrev, p_Radius, p_RadiusPrev, p_dtPrev);
        if (tsRadialExp.has_value()) {
            const double thisDt = p_RadialChangeFraction * tsRadialExpansion;
            dt = (dt > 0.0) ? std::min(dt, thisDt) : thisDt;
        }
    }

    dt = (dt > 0.0) ? std::min(dt, ChooseTimestep()) : ChooseTimestep();

    // There is a chance that mass loss from winds is much faster than previously estimated
    // if, say, LBV winds have turned on.  We therefore precompute the mass loss rate to avoid
    // taking an overly long timestep, despite the extra computational costs.
    auto mlvWinds = CalculateMassLossValues(dt, false);
    double massChangeWinds = p_Mass - std::get<double>(mlvWinds);
    if (CDOUBLE(massChangeWinds) != 0.0) dt = std::min(dt, OPTIONS->MassChangeFraction() * (dt * p_Mass / std::fabs(massChangeWinds)));

    dt = std::max(QUANTISE_DT(dt), NUCLEAR_MINIMUM_TIMESTEP);   // Quantised, and clamped to NUCLEAR_MINIMUM_TIMESTEP

    return dt;
}


/*
 * ChooseTimestep
 *
 * @brief
 * Choose timestep for evolution.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double ChooseTimestep(const double p_Age, const TimescalesT& p_tScales) const
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Timestep (Myr)
 */
inline double BaseStar::ChooseTimestep(const double p_Age, const TimescalesT& p_tScales) const { 

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            return ChooseTimestep_Hurley2000(p_Age, p_tScales);
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ENVELOPE                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateBindingEnergy
 *
 * @brief
 * Calculate the absolute value of the binding energy of the envelope of the star
 *
 *
 * double CalculateBindingEnergy(const double p_CoreMass, const double p_Radius, const double p_EnvMass, const double p_Lambda) const
 *
 * @param           p_CoreMass                  Core mass of the star (Msol)
 * @param           p_Radius                    Radius of the star (Rsol)
 * @param           p_EnvMass                   Mass of the star's envelope (Msol)
 * @param           p_Lambda                    Dimensionless parameter defining the binding energy
 * @return                                      Binding energy of the envelope (erg)
 */
inline double BaseStar::CalculateBindingEnergy(const double p_CoreMass, const double p_Radius, const double p_EnvMass, const double p_Lambda) const {

    // Convert to CGS
    const double Mc   = p_CoreMass * MSOL_TO_G;
    const double Menv = p_EnvMass * MSOL_TO_G;

	return G_CGS * (Mc + Menv) * Menv / (p_Lambda * p_Radius * RSOL_TO_CM);
}


/*
 * CalculateConvectiveEnvelopeBindingEnergy_Hirai2022
 *
 * @brief
 * Calculate the convective envelope binding energy for the two-stage common envelope
 * formalism of Hirai & Mandel, 2022.  Located here because it is a stellar attribute.
 *
 *
 * double CalculateConvectiveEnvelopeBindingEnergy_Hirai2022(const double p_Mass, const double p_Radius, const double p_EnvMass, const double p_Lambda) const
 *
 * @param           p_Mass                      Mass of the star (Msol)
 * @param           p_Radius                    Radius of the star (Rsol)
 * @param           p_EnvMass                   Mass of the star's convective outer envelope (Msol)
 * @param           p_Lambda                    Lambda binding energy parameter for the convective envelope
 * @return                                      Binding energy (erg)
 */
inline double BaseStar::CalculateConvectiveEnvelopeBindingEnergy_Hirai2022(const double p_Mass, const double p_Radius, const double p_EnvMass, const double p_Lambda) const {
    return CalculateBindingEnergy(p_Mass - p_EnvMass, p_Radius, p_EnvMass, p_Lambda);
}


/*
 * CalculateConvectiveEnvelopeLambda_Picker2024
 *
 * @brief
 * Calculates (an approximate value of) the lambda binding energy parameter of the outer convective envelope.
 * Follows the fits of Picker et al., 2024, for lambda_He (https://arxiv.org/abs/2402.13180).
 *
 * This is required for the Hirai & Mandel, 2022, two-stage common envelope formalism, but located
 * here because it is a stellar attribute.
 *
 *
 * double CalculateConvectiveEnvelopeLambda_Picker2024(const double p_Metallicity, const double p_Mass, const double p_EnvMass, const double p_EnvMassMax) const
 *
 * @param           p_Metallicity               Metallicity of the star
 * @param           p_Mass                      Mass of the star (Msol)
 * @param           p_EnvMass                   Mass of the star's convective outer envelope (Msol)
 * @param           p_EnvMassMax                Maximum mass of the star's convective outer envelope (Msol)
 * @return                                      Lambda binding energy parameter for the convective outer envelope
 */
inline double BaseStar::CalculateConvectiveEnvelopeLambda_Picker2024(const double p_Metallicity, const double p_Mass, const double p_EnvMass, const double p_EnvMassMax) const {
    
    const double logZ = ZDEP->SigmaHurley(p_Metallicity);
    const double m2   = 0.0023 * logZ * logZ + 0.0088 * logZ + 0.013;   // Picker et al., 2024, eq 12, and Table 1
    const double b1   = m2 * p_Mass - 0.23;                             // ibid., eq 11
    
    return std::exp((((p_EnvMass / p_EnvMassMax) > 0.3) ? 0.42 * p_EnvMass / p_EnvMassMax + b1 : 0.3 * 0.42 + b1));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           MISC. CONSTANTS / PARAMETERS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateGBParams
 *
 * @brief
 * Calculate the Giant Branch parameters.
 * 
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 * Since Giant Branch parameters depend on a star's mass, they need to be calculated whenever
 * the mass of the star changes (probably every timestep).
 * 
 * p_GBParams is passed by value, modified in place, and returned by value.  No copy on the call,
 * and NRVO constructs GBParams directly in the caller's return slot, so no copy on the return.
 * 
 * 
 * GBParamsT CalculateGBParams(const double p_Mass, GBParamsT& p_GBParams) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBParams                      GB parameters
 * @return                                      Mutated GB parameters
 */
inline GBParamsT BaseStar::CalculateGBParams(const double p_Mass, GBParamsT& p_GBParams) { 

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            p_GBParams = CalculateGBParams_Hurley2000(p_Mass, p_GBParams);
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }
    
    return p_GBParams;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRadius
 *
 * @brief
 * Calculate the radius of the star.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateRadius() const
 *
 * @return                                      Radius (Rsol)
 */
inline double BaseStar::CalculateRadius() const { 

    double radius;

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            radius = CalculateRadius_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return radius;
}


/*
 * CalculateRadiusAtPhaseEnd
 *
 * @brief
 * Calculate the radius of the star at the end of the phase.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateRadiusAtPhaseEnd() const
 *
 * @return                                      Radius (Rsol)
 */
inline double BaseStar::CalculateRadiusAtPhaseEnd() const { 

    double radius;

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            radius = CalculateRadiusAtPhaseEnd_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return radius;
}





/*
 * CalculateRemnantRadius
 *
 * @brief
 * Calculate radius of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass)
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateRemnantRadius() const
 *
 * @return                                      Radius (Rsol)
 */
inline double BaseStar::CalculateRemnantRadius() const { 

    double radius;

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            radius = CalculateRemnantRadius_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return radius;
}




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ROTATION                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateOmegaBreak
 *
 * @brief
 * Calculate the break up angular velocity of a star, in rad yr^1,
 * where [G] = 4*pi^2 AU^3 yr^-2 Msol^-1
 * 
 * 
 * double CalculateOmegaBreak(const double p_Mass, const double p_Radius) const
 *
 * @return                                      Break up angular velocity (rad yr^-1)
 */
inline double BaseStar::CalculateOmegaBreak(const double p_Mass, const double p_Radius) const {
 return _2_PI * std::sqrt(p_Mass / (RSOL_TO_AU * RSOL_TO_AU * RSOL_TO_AU * p_Radius * p_Radius * p_Radius));
}




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    SUPERNOVAE                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * DrawSNkickMagnitude_BrayEldridge2018
 *
 * @brief
 * Draw a kick magnitude, per Bray & Eldridge 2016, 2018
 *
 * See:
 *    https://arxiv.org/abs/1605.09529
 *    https://arxiv.org/abs/1804.04414
 *
 *
 * double DrawSNkickMagnitude_BrayEldridge2018(const double p_EjectaMass, const double p_RemnantMass, const double p_Alpha, const double p_Beta) const
 *
 * @param       p_EjectaMass                    Change in mass of the exploding star (i.e. mass of the ejecta) (Msol)
 * @param       p_RemnantMass                   Mass of the remnant (Msol)
 * @param       p_Alpha                         Fitting coefficient (see Bray & Eldridge 2016, 2018)
 * @param       p_Beta                          Fitting coefficient (see Bray & Eldridge 2016, 2018)
 * @return                                      Drawn kick magnitude (km s^-1)
 */
inline double BaseStar::DrawSNkickMagnitude_BrayEldridge2018(const double p_EjectaMass, const double p_RemnantMass, const double p_Alpha, const double p_Beta) const {
    return p_Alpha * (p_EjectaMass / p_RemnantMass) + p_Beta;
}


/*
 * DrawSNkickMagnitude_Flat
 *
 * @brief
 * Draw a kick magnitude from a uniform distribution between 0 and parameter p_MaxVK.
 *
 *
 * double DrawSNkickMagnitude_Flat(const double p_KickMax, const double p_Rand) const
 *
 * @param       p_KickMax                       Maximum kick magnitude to draw (km s^-1)
 * @param       p_Rand                          Random number for drawing from the distribution [0, 1]
 * @return                                      Drawn kick magnitude (km s^-1)
 */
inline double BaseStar::DrawSNkickMagnitude_Flat(const double p_KickMax, const double p_Rand) const {
    return p_Rand * p_KickMax;
}


/*
 * DrawSNkickMagnitude_Maxwellian
 *
 * @brief
 * Draw a kick magnitude from a Maxwellian distribution.  A Maxwellian is a chi
 * distribution with three degrees of freedom - uses GSL gsl_cdf_chisq_Pinv().
 *
 *
 * double DrawSNkickMagnitude_Maxwellian(const double p_Sigma, const double p_Rand) const
 *
 * @param       p_Sigma                         Distribution scale parameter - affects the spread of the distribution
 * @param       p_Rand                          Random number for drawing from the inverse CDF of the Maxwellian [0, 1]
 * @return                                      Drawn kick magnitude (km s^-1)
 */
inline double BaseStar::DrawSNkickMagnitude_Maxwellian(const double p_Sigma, const double p_Rand) const {
    return p_Sigma * std::sqrt(gsl_cdf_chisq_Pinv(p_Rand, 3));
}













/*
 * CalculateHurleyPerturbationS
 *
 * Calculate the Hurley perturbation parameter, s, per Hurley et al. 2000, eq 101
 *
 *
 * double CalculateHurleyPerturbationS(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Mu                            Hurley perturbation parameter, mu
 * @return                                      Hurley perturbation parameter, s
 */
inline double BaseStar::CalculateHurleyPerturbationS(const double p_Mass, const double p_Mu) const {

    const double b      = CalculateHurleyPerturbationB(p_Mass);
    const double b3     = b * b * b;                // pow() is slow - use multiplication
    const double mu3_b3 = p_Mu * p_Mu * p_Mu / b3;  // calculate once

    return ((1.0 + b3) * mu3_b3) / (1.0 + mu3_b3);
}


/*
 * CalculateHurleyPerturbationR
 *
 * Calculate the Hurley perturbation parameter, r, per Hurley et al. 2000, eq 102
 *
 *
 * double CalculateHurleyPerturbationR(const double p_Mass, const double p_Radius, const double p_Mu, const double p_RemnantRadius) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Mu                            Hurley perturbation parameter, mu
 * @param       p_Rc                            Radius that the remnant would have if the star immediately lost its envelope (Rsol)
 * @return                                      Perturbation parameter, r
 */
inline double BaseStar::CalculateHurleyPerturbationR(const double p_Mass, const double p_Radius, const double p_Mu, const double p_Rc) const {

    double r = 0.0;

    // we only perturb the radius if mu > 0 and radius is larger than core radius,
    // otherwise r = 0 and perturbed radius = core radius
    if (p_Mu > 0.0 && p_Radius > p_Rc) {

        const double c      = CalculateHurleyPerturbationC(p_Mass);
        const double c3     = c * c * c;                    // pow() is slow - use multiplication
        const double mu3_c3 = p_Mu * p_Mu * p_Mu / c3;      // calculate once

        // Hurley et al. 2000 is just 0.1 / q here, but the Hurley sse code does this (`rpertf()` in `zfuncs.f`) - no explanation.
        const double exponent = std::min((0.1 / CalculatePerturbationQ(p_Radius, p_Rc)), (-14.0 / std::log10(p_Mu)));

        r = ((1.0 + c3) * mu3_c3 * PPOW((p_Mu), exponent)) / ((1.0 + mu3_c3));
    }

    return r;
}


/*
 * CalculateLuminosity_Hurley2000
 *
 * @brief
 * Calculate luminosity, given the core mass of the star, per Hurley et al. 2000, eq 37
 *
 *
 * double CalculateLuminosity_Hurley2000(const double p_CoreMass, const GBParamsT& p_GBParams) const
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_GBParams                      GB parameters
 * @return                                      AGB luminosity (Lsol)  <<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS - SHOULD MOVE TO FGB  <<<<<<<<<<<<<<<<<<<<<
 */
inline double BaseStar::CalculateLuminosity_Hurley2000(const double p_CoreMass, const GBParamsT& p_GBParams) const {
    return std::min((p_GBParams[HURLEY_GBP::B] * PPOW(p_CoreMass, p_GBParams[HURLEY_GBP::Q])), (p_GBParams[HURLEY_GBP::D] * PPOW(p_CoreMass, p_GBParams[HURLEY_GBP::P])));
}


/*
 * CalculateLuminosityAtBAGB_Hurley2000
 *
 * @brief
 * Calculate luminosity at the base of the Asymptotic Giant Branch,
 * per Hurley et al. 2000, eq 56
 *
 *
 * double CalculateLuminosityAtBAGB_Hurley2000(const double p_Metallicity, const double p_Mass) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Effective initial mass (Msol)
 * @return                                      BAGB luminosity (Lsol)
 */
inline double BaseStar::CalculateLuminosityAtBAGB_Hurley2000(const double p_Metallicity, const double p_Mass) const {
    const DblVectorT b  = ZDEP->HurleyBCoefficients(p_Metallicity);
    const double MHeF = ZDEP->HurleyMassCutoffs(p_Metallicity, HURLEY_MCO::HeF); 
    return p_Mass < MHeF
            ? (b[29] * PPOW(p_Mass, b[30])) / (1.0 + (ZDEP->HurleyAlpha3(p_Metallicity) * std::exp(15.0 * (p_Mass - MHeF))))
            : (b[31] + (b[32] * PPOW(p_Mass, (b[33] + 1.8)))) / (b[34] + PPOW(p_Mass, b[33]));
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    TEMPERATURE                                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////









////////////////////////// DO WE REALLY NEED THIS???????????????????????????????????????????????
/*
 * Calculate the effective temperature of the star in Kelvin, given the luminosity of the
 * star (in Lsol) and the radius of the star (in Rsol)
 *
 *
 * double CalculateTemperatureKelvinOnPhase(const double p_Luminosity, const double p_Radius)
 *
 * @param   [IN]    p_Luminosity                Luminosity of the star (Lsol)
 * @param   [IN]    p_Radius                    Radius of the star (Rsol)
 * @return                                      Effective temperature of the star (Kelvin)
 */
//double BaseStar::CalculateTemperatureKelvinOnPhase(const double p_Luminosity, const double p_Radius) const {
//    return CalculateTemperatureOnPhase(p_Luminosity, p_Radius) * TSOL;
//}




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCoreMass
 *
 * @brief
 * Calculate the core mass of the star.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateCoreMass(const double p_CoreMass) const
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Core mass (Msol)
 */
inline double BaseStar::CalculateCoreMass(const double p_CoreMass) const { 

    double coreMass = p_CoreMass;                                   // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            coreMass = CalculateCoreMass_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return coreMass;
}


/*
 * CalculateCoreMassAtPhaseEnd
 *
 * @brief
 * Calculate the core mass of the star at the end of the current evolutionary phase.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateCoreMassAtPhaseEnd(const double p_CoreMass) const
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Core mass (Msol)
 */
inline double BaseStar::CalculateCoreMassAtPhaseEnd(const double p_CoreMass) const { 

    double coreMass = p_CoreMass;                                   // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            coreMass = CalculateCoreMassAtPhaseEnd_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return coreMass;
}


/*
 * CalculateCoreMass_Hurley2000
 *
 * @brief
 * Calculate the core mass of the star, per Hurley et al. 2000, eqs 37 & 38.
 * 
 * The Hurley et al. 2000 Mc-L relation is described in section 5.2 of the paper (just before eq 37).
 * 
 * 
 * double CalculateCoreMass_Hurley2000(const double p_Luminosity, const GBParamsT& p_GBParams)
 * 
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_GBParams                      GB parameters
 * @return                                      Core mass of the star (Msol)
 */
inline double BaseStar::CalculateCoreMass_Hurley2000(const double p_Luminosity, const GBParamsT& p_GBParams) {
    return (p_Luminosity > p_GBParams[HURLEY_GBP::Lx])
            ? PPOW((p_Luminosity / p_GBParams[HURLEY_GBP::B]), (1.0 / p_GBParams[HURLEY_GBP::q]))
            : PPOW((p_Luminosity / p_GBParams[HURLEY_GBP::D]), (1.0 / p_GBParams[HURLEY_GBP::p]));
}


/*
 * CalculateCOCoreMass
 *
 * @brief
 * Calculate the carbon-oxygen (CO) core mass of the star.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateCOCoreMass(const double p_COCoreMass) const
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @return                                      CO core mass (Msol)
 */
inline double BaseStar::CalculateCOCoreMass(const double p_COCoreMass) const { 

    double COCoreMass = p_COCoreMass;                               // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            COCoreMass = CalculateCOCoreMass_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return COCoreMass;
}


/*
 * CalculateCOCoreMassAtPhaseEnd
 *
 * @brief
 * Calculate the carbon-oxygen (CO) core mass of the star at the end of the current evolutionary
 * phase.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateCOCoreMassAtPhaseEnd(const double p_COCoreMass) const
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @return                                      CO core mass (Msol)
 */
inline double BaseStar::CalculateCOCoreMassAtPhaseEnd(const double p_COCoreMass) const { 

    double COCoreMass = p_COCoreMass;                               // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            COCoreMass = CalculateCOCoreMassAtPhaseEnd_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return COCoreMass;
}


/*
 * CalculateEffectiveInitialMass
 *
 * @brief
 * Calculate the effective initial mass of the star.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateEffectiveInitialMass(const double p_EffectiveInitialMass) const
 *
 * @param       p_EffectiveInitialMass          Effective initial mass of the star (Msol)
 * @return                                      CO core mass (Msol)
 */
inline double BaseStar::CalculateEffectiveInitialMass(const double p_EffectiveInitialMass) const { 

    double effectiveInitialMass = p_EffectiveInitialMass;           // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            effectiveInitialMass = CalculateEffectiveInitialMass_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return effectiveInitialMass;
}


/*
 * CalculateHeCoreMass
 *
 * @brief
 * Calculate the Helium (He) core mass of the star.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateHeCoreMass(const double p_HeCoreMass) const
 *
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @return                                      He core mass (Msol)
 */
inline double BaseStar::CalculateHeCoreMass(const double p_HeCoreMass) const { 

    double HeCoreMass = p_HeCoreMass;                               // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            HeCoreMass = CalculateHeCoreMass_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return HeCoreMass;
}


/*
 * CalculateHeCoreMassAtPhaseEnd
 *
 * @brief
 * Calculate the Helium (He) core mass of the star at the end of the current evolutionary
 * phase.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateHeCoreMassAtPhaseEnd(const double p_HeCoreMass) const
 *
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @return                                      He core mass (Msol)
 */
inline double BaseStar::CalculateHeCoreMassAtPhaseEnd(const double p_HeCoreMass) const { 

    double HeCoreMass = p_HeCoreMass;                               // Default return value

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            HeCoreMass = CalculateHeCoreMassAtPhaseEnd_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return HeCoreMass;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               MASS LOSS / ACCRETION                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// removed in PR 1445
/*
 * CalculateMLfractionWR
 *
 * @brief
 * Calculate the fraction of mass loss attributable to WR mass loss, per Yoon et al. 2006
 *
 * The model described in Yoon et al. 2006 (also Szecsi et al. 2015) uses OB mass loss while the
 * He surface abundance is below 0.55, WR mass loss when the surface He abundance is above 0.7,
 * and linearly interpolate when the He surface abundance is between those limits.
 *
 * This function calculates the fraction of mass loss attributable to WR mass loss, based on
 * the He surface abundance and the abundance limits described in Yoon et al. 2006.  The value
 * returned will be 1.0 if 100% of the mass loss is attributable to WR mass lass, 0.0 if 100% of
 * the mass loss is attributable to OB mass loss, and in the range (0.0, 1.0) if the mass loss is
 * a mix of WR and OB.
 *
 *
 * double CalculateMLfractionWR(const double p_HeAbundanceSurface) const
 *
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Fraction of mass loss attributable to WR mass loss
 //
inline double BaseStar::CalculateMLfractionWR(const double p_HeAbundanceSurface) const {

    constexpr double limOB = 0.55;                                          // per Yoon et al. 2006
    constexpr double limWR = 0.70;                                          // per Yoon et al. 2006

    return std::min(1.0, std::max (0.0, (p_HeAbundanceSurface - limOB) / (limWR - limOB)));
}
*/


/*
 * CalculateMLRate_Hurley2000
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, per Hurley et al. 2000.
 * 
 * 
 * MassLossT CalculateMLRate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, double) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * << 1 unnamed double parameter to match the virtual signature >>
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (may be ML_TYPE::NONE)
 */
inline MassLossT BaseStar::CalculateMLRate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, double) const {
    return CalculateMLRate_NieuwenhuijzenDeJager1990(p_Mass, p_Radius, p_Luminosity);
}


/*
 * CalculateMLRate_KudritzkiReimers1978
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, on the GB and beyond,
 * per Hurley et al. 2000, eq 106
 * 
 * Based on a prescription taken from Kudritzki and Reimers 1978
 *
 *
 * MassLossT CalculateMLRate_KudritzkiReimers1978(const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::GB)
 * @return                                      Kudritzki and Reimers mass loss rate (in Msol yr^{-1})
 */
inline MassLossT BaseStar::CalculateMLRate_KudritzkiReimers1978(const double p_Mass, const double p_Radius, const double p_Luminosity) const {
    // Hurley et al. 2000 has eta^2 - that's wrong per Hurley SSE code
    return std::make_tuple(4.0E-13 * (MASS_LOSS_ETA * p_Luminosity * p_Radius / p_Mass), ML_TYPE::GB);
}


/* 
 * CalculateMLRate_VassiliadisWood1993
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for stars on the AGB,
 * based on the Mira pulsation period (P0), per Vassiliadis and Wood 1993.
 *
 * See Hurley et al. 2000, just after eq 106
 * Note: in the Hurley fortran code, P0 is taken to be min(p0, 2000.0) - implemented here as a minimum power
 *
 *
 * MassLossT CalculateMLRate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate on the AGB (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::GB)
 */
inline MassLossT BaseStar::CalculateMLRate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const {
    const double P0   = PPOW(10.0, std::min(3.3, (-2.07 - (0.9 * std::log10(p_Mass)) + (1.94 * std::log10(p_Radius)))));
    const double dMdt = PPOW(10.0, (-11.4 + (0.0125 * (P0 - 100.0 * std::max((p_Mass - 2.5), 0.0)))));

    return std::make_tuple(std::min(dMdt, (1.36E-9 * p_Luminosity)), ML_TYPE::GB);
}


/*
 * CalculateMLRateOB_Krticka2018
 *
 * @brief
 * Calculate mass loss rate for massive OB stars, per Krticka et al. 2018
 * (See https://arxiv.org/pdf/1712.03321.pdf)
 *
 * 
 * MassLossT CalculateMLRateOB_Krticka2018(const double p_Luminosity) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate for hot OB stars (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::OB)
 */
COMPAS_PURE inline MassLossT BaseStar::CalculateMLRateOB_Krticka2018(const double p_Metallicity, const double p_Luminosity) const {
    const double zeta = ZDEP->ZetaAsplund(p_Metallicity);
    return std::make_tuple(PPOW(10.0, -5.70 + 0.50 * zeta + (1.61 - 0.12 * zeta) * std::log10(p_Luminosity / 1.0E6)), ML_TYPE::OB);
}


/*
 * CalculateMLRateRSG_Beasor2020
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for RSG stars (Red Supergiants),
 * per Beasor et al. 2020, eq 4.
 * 
 * See https://arxiv.org/pdf/2001.07222.pdf
 * 
 * fit corrected slightly in Decin 2023, eq E.1
 * (see https://arxiv.org/pdf/2303.09385.pdf)
 * 
 * corrected again by Beasor et al. 2023
 * (see https://ui.adsabs.harvard.edu/abs/2023MNRAS.524.2460B/abstract)
 * 
 * 
 * MassLossT CalculateMLRateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const
 *
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   double  RSG mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::RSG)
 */
inline MassLossT BaseStar::CalculateMLRateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const {
    return std::make_tuple(PPOW(10.0, -21.5 - 0.15 * p_mStart + 3.6 * std::log10(p_Luminosity)), ML_TYPE::RSG);
}


/*
 * CalculateMLRateRSG_Decin2023
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for RSG stars (Red Supergiants),
 * per Decin 2023, eq 6.
 *
 * See https://arxiv.org/pdf/2303.09385.pdf
 * 
 *  
 * MassLossT CalculateMLRateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const
 *
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   double  RSG mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::RSG)
 */
inline MassLossT BaseStar::CalculateMLRateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const {
    return std::make_tuple(PPOW(10.0, -20.63 - 0.16 * p_mStart + 3.47 * std::log10(p_Luminosity)), ML_TYPE::RSG);
}


/*
 * CalculateMLRateRSG_VinkSabhahit2023
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for Red Supergiant (RSG) stars,
 * per Vink and Sabhahit 2023, eqs 1 and 2.
 * 
 * A kinked function of L and M
 * 
 * See https://arxiv.org/pdf/2309.08657.pdf
 *
 * 
 * MassLossT CalculateMLRateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   double  RSG mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::RSG)
 */
inline MassLossT BaseStar::CalculateMLRateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const {
    const double logL = std::log10(p_Luminosity);
    const double exp  = logL < 4.6 ? -8.0 + 0.7 * logL - 0.7 * std::log10(p_Mass) : -24.0 + 4.77 * logL - 3.99 * std::log10(p_Mass);
    return std::make_tuple(PPOW(10.0, exp), ML_TYPE::RSG);
}


/*
 * CalculateMLRateRSG_Yang2023
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for RSG stars (Red Supergiants),
 * per Yang 2023, eq 6.
 * 
 * See https://arxiv.org/pdf/2303.09385.pdf
 *
 *  
 * MassLossT CalculateMLRateRSG_Yang2023(const double p_Luminosity) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   double  RSG mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::RSG)
 */
inline MassLossT BaseStar::CalculateMLRateRSG_Yang2023(const double p_Luminosity) const {

    const double logL  = std::log10(p_Luminosity);
    const double logL2 = logL * logL;

    return std::make_tuple(PPOW(10.0, 0.45 * logL2 * logL - 5.26 * logL2 + 20.93 * logL - 34.56), ML_TYPE::RSG);
}


/*
 * CalculateMLRateVMS_Bestenlehner2020
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for very massive (> 100 Msol) OB stars,
 * per Bestenlehner 2020
 * 
 * See https://arxiv.org/pdf/2002.05168.pdf
 * 
 * 
 * MassLossT CalculateMLRateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate for very massive stars (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::VMS)
 */
inline MassLossT BaseStar::CalculateMLRateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const {

    constexpr double alpha       = 0.39;                                        // CAK force multiplier
    constexpr double logMdotZero = -4.78;                                       // From substituting LogMdotTrans and Gamma_e trans into eq 12. 

    const double gamma = EDDINGTON_PARAMETER_FACTOR * p_Luminosity / p_Mass;    // Eddington Parameter, not metallicity specific as in the publication

    return std::make_tuple(PPOW(10.0, logMdotZero + ((1.0 / alpha) + 0.5) * std::log10(gamma) - (((1.0 - alpha) / alpha) + 2.0) * std::log10(1.0 - gamma)), ML_TYPE::VMS);
}


/*
 * CalculateMLRateWR_Hurley2000
 *
 * @brief
 * Calculate the Wolf-Rayet like mass loss rate for small hydrogen-envelope mass (when mu < 1.0),
 * per Hurley et al. 2000, just after eq 106
 * 
 * Taken from Hamann, Koesterke & Wessolowski 1995, Hamann & Koesterke 1998.
 * 
 * In the Hurley SSE code there is a parameter 'hewind' which by default is 1.0, but it can be set
 * to zero to disable this particular part of winds. We instead opt for all winds on or off.
 *
 * Note that the reduction of this formula is imposed to match the observed number of black holes in
 * binaries (Hurley et al. 2000)
 * 
 *
 * MassLossT CalculateMLRateWR_Hurley2000(const double p_Luminosity, const double p_Perturb) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Perturb                       Small envelope perturbation parameter
 * @return                                      Tuple containing:
 *                                                   double  WR mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::WR)
 */
inline MassLossT BaseStar::CalculateMLRateWR_Hurley2000(const double p_Luminosity, const double p_Perturb) const {
    return std::make_tuple((p_Perturb >= 1.0 ? 0.0 : PPOW(p_Luminosity, 1.5) * (1.0 - p_Perturb) * 1.0E-13), ML_TYPE::WR);
} 


/*
 * CalculateMLRateWR_Shenar2019
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, forWolf-Rayet stars,
 * per Shenar et al. 2019, eq 6, tbl 5.
 * 
 * See (https://ui.adsabs.harvard.edu/abs/2019A%26A...627A.151S/abstract)
 * 
 * We use the fitting coefficients for hydrogen rich WR stars (e.g., WNh).
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
COMPAS_PURE inline MassLossT BaseStar::CalculateMLRateWR_Shenar2019(const double p_Metallicity, const double p_Luminosity, const double p_Temperature) const {

    // For H-rich WR stars (X_H > 0.4)
    constexpr double C1 = -6.78;
    constexpr double C2 =  0.66;
    constexpr double C3 = -0.12;
    constexpr double C5 =  0.74;

    return std::make_tuple(PPOW(10.0, C1 + (C2 * std::log10(p_Luminosity)) + (C3 * std::log10(p_Temperature * TSOL)) + (C5 * ZDEP->SigmaHurley(p_Metallicity))), ML_TYPE::WR);
}


/*
 * CalculateMLRateWR_ZDependent
 *
 * @brief
 * Calculate the Wolf-Rayet like mass loss rate for small hydrogen-envelope mass (mu < 1.0),
 * per Belczynski et al. 2010, eq 9
 * (taken from Hamann, Koesterke & Wessolowski 1995, Hamann & Koesterke 1998)
 *
 * Note that the reduction of this formula is imposed to match the observed number of black holes
 * in binaries (Hurley et al. 2000)
 *
 *
 * MassLossT CalculateMLRateWR_ZDependent(const double p_Metallicity, const double p_Luminosity, const double p_Perturb)
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Perturb                       Small envelope perturbation parameter
 * @return                                      Tuple containing:
 *                                                   double  WR mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (will be ML_TYPE::WR)
 */
inline MassLossT BaseStar::CalculateMLRateWR_ZDependent(const double p_Metallicity, const double p_Luminosity, const double p_Perturb) {
    // *Ilya* comments below
    // StarTrack may still do something different here.
    // There are references to Hamann & Koesterke 1998 and Vink and de Koter 2005.
    // TW - Haven't seen StarTrack but I think H&K gives the original equation and V&dK gives the Z dependence

    // Calculate p_Metallicity / ZSOL_ANDERS dierctly rather than use ZDEP->ZScaledAnders() 
    const double dMdt = (p_Perturb >= 1.0)
                            ? 0.0
                            : 1.0E-13 * p_Luminosity * std::sqrt(p_Luminosity) * PPOW(p_Metallicity / ZSOL_ANDERS, 0.86) * (1.0 - p_Perturb);

    return std::make_tuple(dMdt, ML_TYPE::WR);
}    







///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateEddingtonLuminosity
 *
 * @brief
 * Calculate the Eddington luminosity for the star
 * See e.g., above Equation 6 in Bjorklund et al. 2022 (https://arxiv.org/abs/2203.08218)
 * 
 * 
 * double CalculateEddingtonLuminosity(const double p_Mass, const double p_HeliumAbundanceSurface) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Eddington luminosity (Lsol)
 */
inline double BaseStar::CalculateEddingtonLuminosity(const double p_Mass, const double p_HeAbundanceSurface) const {
    const double HeAbundance2 = p_HeAbundanceSurface + p_HeAbundanceSurface;   
    return (4.0 * M_PI * G * C * p_Mass * MSOL_TO_KG) / ((0.4 * (1.0 + 2.0 * HeAbundance2) / (1.0 + 4.0 * (HeAbundance2 + HeAbundance2))) * OPACITY_CGS_TO_SI);
}


/*
 * CalculateLuminosity
 *
 * @brief
 * Calculate the luminosity of the star.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateLuminosity() const
 *
 * @return                                      Luminosity (Lsol)
 */
inline double BaseStar::CalculateLuminosity() const { 

    double luminosity;

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            luminosity = CalculateLuminosity_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return luminosity;
}


/*
 * CalculateLuminosityAtPhaseEnd
 *
 * @brief
 * Calculate the luminosity of the star at the end of the current phase.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateLuminosityAtPhaseEnd() const
 *
 * @return                                      Luminosity (Lsol)
 */
inline double BaseStar::CalculateLuminosityAtPhaseEnd() const { 

    double luminosity;

    switch (OPTIONS->Mode()) {                                      // Which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                            // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                            // HURLEY BSE
            luminosity = CalculateLuminosityAtPhaseEnd_Hurley2000();
            break;
        
        default:                                                    // Unexpected mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);          // Throw error
    }       

    return luminosity;
}
















// <<<<<<<<<<<<<<<<<<<<<<<< CONSTITUENT CLASS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


class BaseStar_Constituent {
    
protected:


public:



};











/*
 * CalculateTimescale_RadialExpansion_DuringMT
 *
 * @brief       // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  JR FIX THIS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 * Calculate the radial expansion timescale in the mass transfer regime.
 * We do not use CalculateTimescale_RadialExpansion() because in the process of
 * mass transfer the previous radius is determined by binary evolution, not nuclear timescale evolution
 *
 *
 * double CalculateTimescale_RadialExpansion_DuringMT()
 *
 * @return                                      Radial expansion timescale
 */
inline double BaseStar::CalculateTimescale_RadialExpansion_DuringMT() { // <<<<<<<<<<<<<<<<<<  JR FIX THIS - PASS STELLAR TYPE
    
    // We create a clone of this star, with the same stellar type, and age it slightly to
    // determine how the radius will change.
    // To be sure the clone does not participate in logging, we set its persistence to EPHEMERAL.
    std::unique_ptr<BaseStar> clone = CloneAs(StellarType(), OBJECT_PERSISTENCE::EPHEMERAL, false); // Do not re-initialise the clone

    const double timestep = std::max(1000.0 * NUCLEAR_MINIMUM_TIMESTEP, Age() / 1.0E6);             // Timestep to age clone (Myr)
    (void)clone->AdvanceOneTimestep(timestep);                                                      // Age clone by timestep - wind mass loss applied internally

    return timestep * Radius() / std::fabs(Radius() - clone->Radius());
}












///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 ZETA CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateZetaAdiabatic_Hurley2002                <<< ONLY CALLED (INDIRECTLY) FROM BASEBINARYSTAR   -  INLINE THIS
 *
 * @brief
 * Calculate the adiabatic exponent (donor radial response to mass loss), zeta,
 * per Hurley et al. 2002
 *
 *
 * double CalculateZetaAdiabatic_Hurley2002(const double p_Mass, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Adiabatic exponent, zeta
 */
inline COMPAS_PURE double BaseStar::CalculateZetaAdiabatic_Hurley2002(const double p_Mass, const double p_CoreMass) const {    

    // JR FIX THIS: check this <<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return (p_CoreMass >= p_Mass) ? 1.0 : -ZDEP->HurleyGBRadiusXExponent(Metallicity()) + (2.0 * utils::intPow(p_CoreMass / p_Mass, 5));
}








///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                        COMMON ENVELOPE LAMBDA CALCULATIONS                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/* 
 * CalculateCELambda_Kruckow
 *
 * @brief
 * Calculate the common envelope lambda parameter from Kruckow et al., 2016, fig 1, 
 * per Vigna-Gomez et al., 2018 
 *
 * Spectrum fit to the region bounded by the upper and lower limits as shown in Kruckow et al., 2016
 * (https://arxiv.org/abs/1610.04417); fit as presented in Vigna-Gomez et al. 2018 (https://arxiv.org/abs/1805.07974)
 *
 *
 * double CalculateCELambda_Kruckow(const double p_Radius, const double p_Alpha) const
 *
 * @param           p_Radius                    Radius of the star (Rsol)
 * @param           p_Alpha                     Common envelope slope (see `--common-envelope-slope-kruckow`)
 * @return                                      Common envelope lambda parameter
 */
inline double BaseStar::CalculateCELambda_Kruckow(const double p_Radius, const double p_Alpha) const {
	const double beta = std::max(-2.0 / 3.0, std::min(-1.0, p_Alpha));    // beta in Vigna-Gomez et al. 2018, sec 3.2.8; clamp to [-1.0, -2.0/3.0]
	return 1600.0 * PPOW(0.00125, -beta) * PPOW(p_Radius, beta);
}


