#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <list>
#include <memory>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <variant>

class BinaryConstituentStar;

#include "constants.h"
#include "ErrorCatalog.h"
#include "ErrorsMacros.h"
#include "utils.h"
#include "typealiases.h"
#include "vector3d.h"
#include "iterators.h"

// STATE_CHUNK_SIZE and DEFAULT_STATE_HISTORY_CAPACITY are defined in constants.h




// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS DOCUMENTATION


// We implement a state stack to record the state hsistory of systems (stars, and binaries).
//
// The single star state history stack contains state snapshots of type StarState, and the
// binary state history stack contains state snapshots of type StellarBinaryState, which in
// itself contains simultaneous snapshots of the (two) constituents stars of the binary, of
// type BinaryStarState (which includes a pointer to the StarState ). <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS
//
//

/*---



implementation is header file only

state history of a star or binary consists of the the starting state (aka the "zero state"),
the current state, and at least one intervening state.  The state immediately prior to the
current state is referred to as the previous state (from the perspective of the star or 
binary in its current state).

explain extra interveninging states - for trend analysis

the state history is stored and maintained as a LIFO doubly-linked list - with a minor twist
- explain zero state never falls off the end of the stack, and is never popped or - oldest non-zero state is removed


capacity must be 0 or >= 3 (zero state, current state, previous state); 0 indicates unbounded (no elements pop when new element pushed)
current state is last element (back of list)
previous state is second-last element (element immediately prior to back of list)
zero state is locked at first element (front of list)
acts like LIFO list, except that first element is locked - pushing a new element pops second element (elsement immediately after front of list)




Stack operations

Any stack operation function that returns a stack element, returns that element as
an object (of type T, and the object contains a copy of the element), rather than a
pointer to the actual stack element.  Returning a pointer (even a const pointer to
const data) to the actual stack element would cause problems if the stack element
was popped or removed due to a push when stack full (the pointer held by the caller
would no longer be valid), so a copy of the stack element is returned).

Sometimes, a stack operation function may not have a (copy of a) stack element to
return.  For example, in the case of the Peek() function, if the caller specifies
a p_Count value that is outside the range of the stack (i.e the caller asks to peek
at a stack element that does not exist), the Peek() function cannot return a copy
of the stack element requested.  Since the function cannot return a null pointer
(because it doesn't return a pointer, for the reasons stated above), the returned
object (nominally a copy of a stack element of type T) is wrapped in std::optional
(i.e. the object returned is of type std::optional<T>).  In the cases that the
requested element does not exist on the stack, the returned object will not be
assigned a value, and the caller should check for that circumstance.

Standard std::optional<> semantics apply.  The caller should use (something like):

      auto elem = Peek(n);

      to retrieve stack element n, measured from TOS, and

      if (elem)
          or 
      if (elem.has_value())

      to check if Peek(n) returned a non-empty object, and

      elem.value().Mass()
              or
      (*elem).Mass()
              or
      elem->Mass()

      to get the value of (e.g.) the mass variable from the stack element

The construct 'elem.value_or(default_value)' is also available - this returns the
value "default_value" in the case that elem has no value assigned.

Note that attempting to query the value of a std::optional object that has not had a
value assigned will result in undefined bahaviour.

---*/













// The clean generalization is event channels: each channel pairs a detection predicate with its own list of iterators into the main history. 
// The owner registers the channels relevant to its state type, and Push runs every registered predicate. 
// Adding a future event type is then just "register another channel."
// One enum spanning all event kinds. Each StateHistory only registers the channels that apply to its state type — StateHistory<StarState> 
// never registers COMMON_ENVELOPE, and that's fine.
enum class EventType {
    SWITCH,            // stellar-type change (single star / constituent)
    SUPERNOVA,         // SN event (single star / constituent)
    MASS_TRANSFER,     // mass-transfer event (constituent)
    COMMON_ENVELOPE,   // CE event (binary)
    // add future event types here
};






// StarState class - state for single stars
class StarState {       

private:

    std::optional<ERROR>             m_Error;

    std::optional<EVOLUTION_STATUS>  m_EvolutionStatus;

    std::optional<STELLAR_TYPE>      m_StellarType;

    OptDblT                          m_Metallicity;                     // Can vary over evolution

    OptDblT                          m_dt;                              // Timestep - delta time (Myr)
    OptDblT                          m_Age;                             // Stellar age (Myr)
    OptDblT                          m_Tau;                             // Phase-relative stellar age
    OptDblT                          m_Time;                            // Simulation time (Myr)       

    OptDblT                          m_CNOProcessedCoreMass;            // CNO-processed core mass (Msol).
    OptDblT                          m_COCoreMass;
    OptDblT                          m_CoreMass;
    OptDblT                          m_HeCoreMass;
    OptDblT                          m_Mass;
    OptDblT                          m_MassEffectiveInitial;            // Effective initial mass of the star (Msol) (Hurley et al., 2000)

    OptDblT                          m_dMdt;

    std::optional<ML_TYPE>    m_DominantMassLossType;

    OptDblT                          m_HAbundanceCore;
    OptDblT                          m_HAbundanceSurface;
    OptDblT                          m_HeAbundanceCore;
    OptDblT                          m_HeAbundanceCoreOut;              // Helium abundance in the region OUTSIDE the core (distinct from core/surface). Brcek/Shikauchi MS core tracking.
    OptDblT                          m_HeAbundanceSurface;

    OptDblT                          m_Luminosity;
    OptDblT                          m_LuminosityZAMSEffective;

    OptDblT                          m_Radius;
    OptDblT                          m_RadiusZAMSEffective;

    OptDblT                          m_Temperature;
    OptDblT                          m_TemperatureZAMSEffective;

    OptDblT                          m_AngularFrequency;                // Angular frequency (rad/yr)
    OptDblT                          m_AngularMomentum;                 // Angular momentum (Msol*AU^2/yr)
    OptDblT                          m_MomentOfInertia;                 // Moment of inertia (Msol*AU^2)


    OptBoolT                         m_EvolvedCHE;                      // CHE flag - true if the star spent entire MS as a CH star; false if evolved CH->MS
    OptBoolT                         m_LBVphaseFlag; 
    OptBoolT                         m_EnvelopeExpelledByPulsations;    // Records whether the convective envelope has been expelled by pulsations

    std::optional<PulsarDetailsT>    m_PulsarDetails;
              
    std::optional<StellarSNDetailsT> m_SupernovaDetails;
       
    // Timescales and GBParams are each constructed as a std::array (much better for
    // optimising than std::vector, especially when being passed as parameters). 
    //
    // A small wrinkle is that the size of a std::array must be fixed, and known at
    // compile time.  It may be that for different evolution modes, the sizes of the
    // timescales and GBParams arrays differ (if they're even used), and we don't
    // know the evolution mode until run time (and the options are parsed).  To get
    // around this we calculate the size of the arrays of the maximum size required
    // (see below).  We might waste a few bytes (realistically the sizes for different
    // evolutionary modes (i.e. the number of timescales, and the number of GBParams)
    // shouldn't be too different)), but the payoff in performance is likely to be
    // significant.
    //
    // Each array is constructed here with the maximim number of elements.  Each array
    // element is wrapped in std::optional<>, and initialised to std::nullopt.  The
    // arrays are guaranteed to exist, but individual elements won't have values until
    // they are calculated.
    //
    // Sizing the arrays:
    //
    // Assuming we have implemented two evolutionary modes, HURLEY and METISSE (both
    // SSE and BSE modes if applicable, and assuming both modes actually use timescales
    // and GBParams), then we would calculate the maximum sizes of the arrays as:
    //    
    // constexpr TIMESCALES_COUNT = std::max(HURLEY_TIMESCALES::COUNT, METISSE_TIMESCALES::COUNT);
    // constexpr GBPARAMS_COUNT   = std::max(HURLEY_GB_PARAMETERS::COUNT, METISSE_GB_PARAMETERS::COUNT);
    //
    // For now, since we have only implemented the HURLEY evolution modes, we use
    //
    // constexpr TIMESCALES_COUNT = HURLEY_TIMESCALES::COUNT;
    // constexpr GBPARAMS_COUNT   = HURLEY_GB_PARAMETERS::COUNT;
    //
    // The constexpr statements are in constants.h so that the constants can be used
    // throughout the code

    TimescalesT m_Timescales;   
    GBParamsT   m_GBParams;


public:

    StarState() { }


    // Getters
    //
    // Note that the getters for variables wrapped in std::optional<> have two forms, and each has
    // a parameter that informs the compiler (so the decision is made at compile time, not run time)
    // which variant should be called.
    //
    // The default form (which has the extra parameter set to a default value), returns the *value*
    // of the variable requested.
    //
    // For the default form, if the variable being requested is *not* wrapped in std::optional<>,
    // the value of the variable is returned. If the variable *is* wrapped in std::optional<>,
    // then:
    //
    //    (a) if the variable has a value, the value will be returned, but
    //    (b) if the variable *does not* have a value, an error message and stacktrace will be
    //        displayed, and the program terminated.
    //
    // The reasoning for (b) is that for most circumstances, developers will know when state
    // variables have, or are expected to have, defined values, so can just call the getters
    // and get the value of the variable - and if the variable does not have a value that's 
    // probably either a coding problem, or a corruption/timing issue, that needs to be addressed.
    //
    // For the non-default (optional) form (see below for calling syntax), the variable requested
    // is returned wrapped in std::optional<>.  The caller should then check via *.has_value() if
    // the returned variable has a value assigned.  If the variable being requested is *not*
    // actually declared as std::optional<>, the getter will wrap it in std::optional<> before
    // returning it.  Such variables will *always* have a value when returned. Variables that *are*
    // declared as std::optional<> will be returned as is, and may or may not have a value assigned.
    //
    //
    // Array variables
    // ----------------
    //
    // Some state variables are stored as arrays - e.g. Timescales.  All array variables can be
    // retrieved as an array, or as individual elements, but there is overhead in retrieving the
    // individual elements - range checking is performed on the specified index (where appropriate).
    // In some cases it might be more appropriate to retrieve the array, and let the caller do the
    // range checking as approriate.  See below for calling syntax.
    //
    // Getters for state values that are arrays come in two forms:
    //
    //    - One form takes a single parameter that is an index into the array, and returns the
    //      value of the array element at that index.  For this form of the getter, the index is
    //      checked against the bounds of the array, and if it is out of bounds the getter fails
    //      and the program is terminated with ERROR::OUT_OF_BOUNDS.
    //
    //    - The other form takes no parameters and returns a const reference to the array in the
    //      state stack (of type appropriate for the getter called).  This is a pointer into the
    //      state - so that there is no copying of the array (so saving compute cycles).  Since
    //      it is a const reference, the caller can read the values in the stack but cannot change
    //      them.
    //
    //      There is an *important caveat* for this form of the getter: the reference to the array
    //      that is returned could become stale.  The reference is valid only while the state
    //      object is alive, and the state slot hasn't been evicted.  The reference returned is
    //      essentially a pointer to an area of memory occupied by the array in the state slot,
    //      and the best usage pattern for it is "get the reference, read the values required, use
    //      the values, then forget the reference".  The reason is that the state is mutable - an
    //      entry occupying a slot may be evicted to make room if a new state is pushed onto the
    //      state stack when the stack is full (only happens when the stack is limited - for an
    //      unlimited stack nothing is evicted).  So if a caller stashes a reference and then a
    //      new state is pushed onto the stack, the slot behind the reference could be evicted
    //      from the stack, and the stashed reference now reads the wrong data (or, worse, is
    //      pointing at values unrelated to the state stack).
    //
    //      Given the paragraph above, the best thing to do is not stash a reference to a state
    //      array away and reuse it - there is no wat to determine if a stashed reference is stale.
    //
    //
    // Usage and calling syntax
    // ------------------------
    //
    // The namespace "get" is declared to facilitate calling the different getter variants
    // (see typeAliases.h).  // <<<<<<<<<<<<<<<< JR FIX THIS when the final home for namespace get is decided <<<<<<<<<<<<<<<<<<<<<
    //
    // If the caller expects a state variable to have a value assigned (whether it be a variable
    // declared as std::optional<> or not), they should use the default getter variant. For example:
    //
    //     double mass = state.Mass()
    //
    // will return the value of mass as a non-optional double value.
    //
    // The State class member variable m_Mass is declared as std::optional<double>, and there is a
    // window as COMPAS is initialising where it may not have yet been assigned a value.  If the
    // call to state.Mass() is made before the value of m_Mass has been set, the call above will
    // fail as described above (error message, stack trace, and program terminated).  If, however,
    // the call is made after the value of m_Mass has been set, the call above will return the value
    // of m_Mass as a variable of type double.
    //
    // Note that using:
    //
    //     double mass = state.Mass(get::unwrappedT)
    //
    // is the same as using 
    //
    //     double mass = state.Mass()
    //
    // (the extra parameter defaults to "get::unwrappedT")
    //
    // On the other hand, if the caller is uncertain if a state variable has a value assigned,
    // they should use the optional getter variant. For example:
    //
    //     std::optional<double> mass = state.Mass(get::asOptionalT)
    //
    // which will return a std::optional variable that can be interrogated as follows:
    //
    //     std::optional<double> mass = state.Mass(get::asOptionalT);
    //     if (mass.has_value()) {
    //         std::cout << "The value of mass is " << mass.value() << "\n";
    //     }
    //     else {
    //         // take some appropriate action here if the variable has no value assigned
    //     }
    //
    // Note that the extra parameter must be supplied to distinguish the call from the default
    // variant of the getter ("get::asOptionalT" must appear as the last parameter in the
    // parameter list).
    //
    // For getters that take a parameter (e.g. getters that return a specific element from
    // an array), the call syntax is e.g.:
    //
    //    double McBGB = GBParams(HURLEY_GBP::McBGB)
    //
    // for the default getter variant, getting the element for HURLEY_GB_PARAMETERS::McBGB, and
    //
    //    OptDblT McBGB = GBParams(HURLEY_GBP::McBGB, get::asOptionalT)
    //
    // for the optional getter variant, getting the element for HURLEY_GB_PARAMETERS::McBGB.

    // Convenience macros for getters
    // Two-form getter for a simple optional field: NAME() / NAME(get::asOptionalT)
    #define STATE_GETTER(NAME, VALTYPE, OPTTYPE, MEMBER)                              \
        VALTYPE NAME(get::unwrappedT = {}) const { return utils::GetValue(MEMBER); }  \
        OPTTYPE NAME(get::asOptionalT   ) const { return MEMBER; }

    #define STATE_GETTER_DBL(NAME, MEMBER)  STATE_GETTER(NAME, double, OptDblT,  MEMBER)
    #define STATE_GETTER_BOOL(NAME, MEMBER) STATE_GETTER(NAME, bool,   OptBoolT, MEMBER)

    STATE_GETTER_DBL (AngularFrequency,             m_AngularFrequency)
    STATE_GETTER_DBL (AngularMomentum,              m_AngularMomentum)
    STATE_GETTER_DBL (Age,                          m_Age)
    STATE_GETTER_DBL (CNOProcessedCoreMass,         m_CNOProcessedCoreMass)
    STATE_GETTER_DBL (COCoreMass,                   m_COCoreMass)
    STATE_GETTER_DBL (CoreMass,                     m_CoreMass)
    STATE_GETTER_DBL (dMdt,                         m_dMdt)
    STATE_GETTER_DBL (dt,                           m_dt)
    STATE_GETTER_BOOL(EnvelopeExpelledByPulsations, m_EnvelopeExpelledByPulsations)
    STATE_GETTER_BOOL(EvolvedCHE,                   m_EvolvedCHE)
    STATE_GETTER_DBL (HAbundanceCore,               m_HAbundanceCore)
    STATE_GETTER_DBL (HAbundanceSurface,            m_HAbundanceSurface)
    STATE_GETTER_DBL (HeAbundanceCore,              m_HeAbundanceCore)
    STATE_GETTER_DBL (HeAbundanceCoreOut,           m_HeAbundanceCoreOut)
    STATE_GETTER_DBL (HeAbundanceSurface,           m_HeAbundanceSurface)
    STATE_GETTER_DBL (HeCoreMass,                   m_HeCoreMass)
    STATE_GETTER_BOOL(LBV_PhaseFlag,                m_LBVphaseFlag)
    STATE_GETTER_DBL (Luminosity,                   m_Luminosity)
    STATE_GETTER_DBL (LuminosityZAMSEffective,      m_LuminosityZAMSEffective)
    STATE_GETTER_DBL (Luminosity0,                  m_LuminosityZAMSEffective)
    STATE_GETTER_DBL (Mass,                         m_Mass)
    STATE_GETTER_DBL (MassEffectiveInitial,         m_MassEffectiveInitial)
    STATE_GETTER_DBL (Mass0,                        m_MassEffectiveInitial)
    STATE_GETTER_DBL (Metallicity,                  m_Metallicity)
    STATE_GETTER_DBL (Radius,                       m_Radius)
    STATE_GETTER_DBL (RadiusZAMSEffective,          m_RadiusZAMSEffective)
    STATE_GETTER_DBL (Radius0,                      m_RadiusZAMSEffective)
    STATE_GETTER_DBL (Tau,                          m_Tau)
    STATE_GETTER_DBL (Temperature,                  m_Temperature)
    STATE_GETTER_DBL (TemperatureZAMSEffective,     m_TemperatureZAMSEffective)
    STATE_GETTER_DBL (Temperature0,                 m_TemperatureZAMSEffective)
    STATE_GETTER_DBL (Time,                         m_Time)

    STATE_GETTER(DominantMassLossType, ML_TYPE,    std::optional<ML_TYPE>,     m_DominantMassLossType)
    STATE_GETTER(Error,                ERROR,             std::optional<ERROR>,              m_Error)
    STATE_GETTER(EvolutionStatus,      EVOLUTION_STATUS,  std::optional<EVOLUTION_STATUS>,   m_EvolutionStatus)
    STATE_GETTER(PulsarDetails,        PulsarDetailsT,    std::optional<PulsarDetailsT>&,    m_PulsarDetails)
    STATE_GETTER(SNDetails,            StellarSNDetailsT, std::optional<StellarSNDetailsT>&, m_SupernovaDetails)
    STATE_GETTER(StellarType,          STELLAR_TYPE,      std::optional<STELLAR_TYPE>,       m_StellarType)

    UnpackedGBParamsT GBParams(get::unwrappedT = {}) const {                            // There is overhead here - use sparingly
        UnpackedGBParamsT unpacked;                                                     // This version of the getter unpacks the optional values
        for (SizeT i = 0; i < GBPARAMS_COUNT; i++) {                                    // For each GB param
            if (!m_GBParams[i].has_value()) utils::Fail(ERROR::NO_VALUE);               // Fail if m_GBParams[i] has no value
            unpacked[i] = m_GBParams[i].value();                                        // Unpack value
        }
        return unpacked;
    }
    const GBParamsT& GBParams(get::asOptionalT) const { return m_GBParams; }
    double GBParams(const SizeT p_Idx, get::unwrappedT = {}) const {                    // Indexed getter
        if (p_Idx < m_GBParams.size()) return utils::GetValue(m_GBParams[p_Idx]);       // Return value at requested index iff index is in bounds
        utils::Fail(ERROR::OUT_OF_BOUNDS);                                              // Index out of bounds - fail
    }
    OptDblT GBParams(const SizeT p_Idx, get::asOptionalT) const {                       // Indexed getter
        if (p_Idx < m_GBParams.size()) return m_GBParams[p_Idx];                        // Return value at requested index iff index is in bounds
        utils::Fail(ERROR::OUT_OF_BOUNDS);                                              // Index out of bounds - fail
    }

    SupernovaEventsT SNEvents(get::unwrappedT = {}) const {
        if (!m_SupernovaDetails.has_value()) utils::Fail(ERROR::NO_VALUE);              // Fail if m_SupernovaDetails has no value
        return m_SupernovaDetails->events;                                              // Guaranteed to exist here
    }
    std::optional<SupernovaEventsT> SNEvents(get::asOptionalT) const {
        if (!m_SupernovaDetails.has_value()) return std::nullopt;                       // No value returns std::nullopt
        return std::optional<SupernovaEventsT>(m_SupernovaDetails->events);             // Optional value (but guaranteed to exist here)
    }

    double SNFallbackFraction(get::unwrappedT = {}) const {
        if (!m_SupernovaDetails.has_value()) utils::Fail(ERROR::NO_VALUE);              // Fail if m_SupernovaDetails has no value
        return m_SupernovaDetails->fallbackFraction;                                    // Guaranteed to exist here
    }
    OptDblT SNFallbackFraction(get::asOptionalT) const {
        if (!m_SupernovaDetails.has_value()) return std::nullopt;                       // No value returns std::nullopt
        return std::optional<double>(m_SupernovaDetails->fallbackFraction);             // Optional value (but guaranteed to exist here)
    }

    bool SNIsHydrogenPoor(get::unwrappedT = {}) const {
        if (!m_SupernovaDetails.has_value()) utils::Fail(ERROR::NO_VALUE);              // Fail if m_SupernovaDetails has no value
        return m_SupernovaDetails->isHydrogenPoor;                                      // Guaranteed to exist here
    }
    OptBoolT SNIsHydrogenPoor(get::asOptionalT) const {
        if (!m_SupernovaDetails.has_value()) return std::nullopt;                       // No value returns std::nullopt
        return std::optional<bool>(m_SupernovaDetails->isHydrogenPoor);                 // Optional value (but guaranteed to exist here)
    }

    StellarKickParmsT SNKickParameters(get::unwrappedT = {}) const {
        if (!m_SupernovaDetails.has_value()) utils::Fail(ERROR::NO_VALUE);              // Fail if m_SupernovaDetails has no value
        return m_SupernovaDetails->kickParameters;                                      // Guaranteed to exist here
    }
    std::optional<StellarKickParmsT> SNKickParameters(get::asOptionalT) const {
        if (!m_SupernovaDetails.has_value()) return std::nullopt;                       // No value returns std::nullopt
        return std::optional<StellarKickParmsT>(m_SupernovaDetails->kickParameters);    // Optional value (but guaranteed to exist here)
    }

    STELLAR_TYPE SNStellarType(get::unwrappedT = {}) const {
        if (!m_SupernovaDetails.has_value()) utils::Fail(ERROR::NO_VALUE);              // Fail if m_SupernovaDetails has no value
        return m_SupernovaDetails->stellarType;                                         // Guaranteed to exist here
    }
    std::optional<STELLAR_TYPE> SNStellarType(get::asOptionalT) const {
        if (!m_SupernovaDetails.has_value()) return std::nullopt;                       // No value returns std::nullopt
        return std::optional<STELLAR_TYPE>(m_SupernovaDetails->stellarType);            // Optional value (but guaranteed to exist here)
    }

    UnpackedTimescalesT Timescales(get::unwrappedT = {}) const {                        // There is overhead here - use sparingly
        UnpackedTimescalesT unpacked;                                                   // This version of the getter unpacks the optional values
        for (SizeT i = 0; i < TIMESCALES_COUNT; i++) {                                  // For each Timescale
            if (!m_Timescales[i].has_value()) utils::Fail(ERROR::NO_VALUE);             // Fail if m_Timescales[i] has no value
            unpacked[i] = m_Timescales[i].value();                                      // Unpack value
        }
        return unpacked;
    }
    const TimescalesT& Timescales(get::asOptionalT) const   { return m_Timescales; }
    double Timescales(const SizeT p_Idx, get::unwrappedT = {}) const {                  // Indexed getter
        if (p_Idx < m_Timescales.size()) return utils::GetValue(m_Timescales[p_Idx]);   // Return value at requested index iff index is in bounds
        utils::Fail(ERROR::OUT_OF_BOUNDS);                                              // Index out of bounds - fail
    }
    OptDblT Timescales(const SizeT p_Idx, get::asOptionalT) const {                     // Indexed getter
        if (p_Idx < m_Timescales.size()) return m_Timescales[p_Idx];                    // Return value at requested index iff index is in bounds
        utils::Fail(ERROR::OUT_OF_BOUNDS);                                              // Index out of bounds - fail
    }

    #undef STATE_GETTER_BOOL
    #undef STATE_GETTER_DBL
    #undef STATE_GETTER


    // Setters & resetters

    // Convenience macros for setters and resetters
    #define STATE_SETTER(NAME, VALTYPE, MEMBER)                          \
        void Reset##NAME()                    { MEMBER = std::nullopt; } \
        void Set##NAME(const VALTYPE p_Value) { MEMBER = p_Value; }

    #define STATE_SETTER_DBL(NAME, MEMBER)  STATE_SETTER(NAME, double, MEMBER)
    #define STATE_SETTER_BOOL(NAME, MEMBER) STATE_SETTER(NAME, bool,   MEMBER)

    STATE_SETTER_DBL (AngularFrequency,             m_AngularFrequency)
    STATE_SETTER_DBL (AngularMomentum,              m_AngularMomentum)
    STATE_SETTER_DBL (Age,                          m_Age)
    STATE_SETTER_BOOL(EvolvedCHE,                   m_EvolvedCHE)
    STATE_SETTER_DBL (CNOProcessedCoreMass,         m_CNOProcessedCoreMass)   
    STATE_SETTER_DBL (COCoreMass,                   m_COCoreMass)
    STATE_SETTER_DBL (CoreMass,                     m_CoreMass)
    STATE_SETTER_DBL (HeCoreMass,                   m_HeCoreMass)
    STATE_SETTER_DBL (DMdt,                         m_dMdt)
    STATE_SETTER_DBL (Dt,                           m_dt)
    STATE_SETTER_BOOL(EnvelopeExpelledByPulsations, m_EnvelopeExpelledByPulsations)
    STATE_SETTER_DBL (HAbundanceCore,               m_HAbundanceCore)
    STATE_SETTER_DBL (HAbundanceSurface,            m_HAbundanceSurface)
    STATE_SETTER_DBL (HeAbundanceCore,              m_HeAbundanceCore)
    STATE_SETTER_DBL (HeAbundanceCoreOut,           m_HeAbundanceCoreOut)
    STATE_SETTER_DBL (HeAbundanceSurface,           m_HeAbundanceSurface)
    STATE_SETTER_BOOL(LBV_PhaseFlag,                m_LBVphaseFlag)
    STATE_SETTER_DBL (Luminosity,                   m_Luminosity)
    STATE_SETTER_DBL (LuminosityZAMSEffective,      m_LuminosityZAMSEffective)
    STATE_SETTER_DBL (Luminosity0,                  m_LuminosityZAMSEffective)
    STATE_SETTER_DBL (Mass,                         m_Mass)
    STATE_SETTER_DBL (MassEffectiveInitial,         m_MassEffectiveInitial)
    STATE_SETTER_DBL (Mass0,                        m_MassEffectiveInitial)
    STATE_SETTER_DBL (Metallicity,                  m_Metallicity)
    STATE_SETTER_DBL (Radius,                       m_Radius)
    STATE_SETTER_DBL (RadiusZAMSEffective,          m_RadiusZAMSEffective)
    STATE_SETTER_DBL (Radius0,                      m_RadiusZAMSEffective)
    STATE_SETTER_DBL (Tau,                          m_Tau)
    STATE_SETTER_DBL (Temperature,                  m_Temperature)
    STATE_SETTER_DBL (TemperatureZAMSEffective,     m_TemperatureZAMSEffective)
    STATE_SETTER_DBL (Temperature0,                 m_TemperatureZAMSEffective)
    STATE_SETTER_DBL (Time,                         m_Time) 

    STATE_SETTER(DominantMassLossType, ML_TYPE,     m_DominantMassLossType)
    STATE_SETTER(Error,                ERROR,              m_Error)
    STATE_SETTER(EvolutionStatus,      EVOLUTION_STATUS,   m_EvolutionStatus)
    STATE_SETTER(PulsarDetails,        PulsarDetailsT&,    m_PulsarDetails)
    STATE_SETTER(StellarType,          STELLAR_TYPE,       m_StellarType)
    STATE_SETTER(SNDetails,            StellarSNDetailsT&, m_SupernovaDetails)

    void ResetGBParams()                                                        { for (SizeT idx = 0; idx < m_GBParams.size(); idx++) m_GBParams[idx] = std::nullopt; }
    void ResetGBParams(const SizeT p_Idx) {                                     // Indexed resetter
        if (p_Idx < m_GBParams.size()) m_GBParams[p_Idx] = std::nullopt;        // Reset value at requested index iff index is in bounds
        else utils::Fail(ERROR::OUT_OF_BOUNDS);                                 // Index out of bounds - fail
    }
    void SetGBParams(const OptDblVectorT& p_GBParams)                           { m_GBParams = p_GBParams; }
    void SetGBParams(const SizeT p_Idx, const double p_Value) {                 // Indexed setter
        if (p_Idx < m_GBParams.size()) m_GBParams[p_Idx] = p_Value;             // Set value at requested index iff index is in bounds
        else utils::Fail(ERROR::OUT_OF_BOUNDS);                                 // Index out of bounds - fail
    }

    void ResetSNEvents() {
        if (m_SupernovaDetails.has_value()) {                                   // SN details exist?
            m_SupernovaDetails->events = {SN_EVENT::NONE, SN_EVENT::NONE};      // Yes, reset events
        }                                                                       // No need if no SN details
    }
    void SetSNEvents(const SupernovaEventsT& p_SNEvents)                        { EnsureSNDetails().events = p_SNEvents; }

    void ResetSNCurrentEvent() {
        if (m_SupernovaDetails.has_value()) {                                   // SN details exist?
            m_SupernovaDetails->events.current = SN_EVENT::NONE;                // Yes, reset current event
        }                                                                       // No need if no SN details
    }
    void ClearSNCurrentEvent()                                                  { ResetSNCurrentEvent(); }
    void SetSNCurrentEvent(const SN_EVENT p_SNEvent)                            { EnsureSNDetails().events.current = p_SNEvent; }

    void ResetSNPastEvent() {
        if (m_SupernovaDetails.has_value()) {                                   // SN details exist?
            m_SupernovaDetails->events.past = SN_EVENT::NONE;                   // Yes, reset past event
        }                                                                       // No need if no SN details
    }
    void ClearSNPastEvent()                                                     { ResetSNPastEvent(); }
    void SetSNPastEvent(const SN_EVENT p_SNEvent)                               { EnsureSNDetails().events.past = p_SNEvent; }

    void ResetSNFallbackFraction() {
        if (m_SupernovaDetails.has_value()) {                                   // SN details exist?
            m_SupernovaDetails->fallbackFraction = {};                          // Yes, reset fallback fraction
        }                                                                       // No need if no SN details
    }
    void SetSNFallbackFraction(const double p_FallbackFraction)                 { EnsureSNDetails().fallbackFraction = p_FallbackFraction; }

    void ResetSNIsHydrogenPoor() {
        if (m_SupernovaDetails.has_value()) {                                   // SN details exist?
            m_SupernovaDetails->isHydrogenPoor = {};                            // Yes, reset hydrogen poor flag
        }                                                                       // No need if no SN details
    }
    void SetSNIsHydrogenPoor(const bool p_IsHydrogenPoor)                       { EnsureSNDetails().isHydrogenPoor = p_IsHydrogenPoor; }

    void ResetSNKickParameters() {
        if (m_SupernovaDetails.has_value()) {                                   // SN details exist?
            m_SupernovaDetails->kickParameters = {};                            // Yes, reset kick parameters
        }                                                                       // No need if no SN details
    }
    void SetSNKickParameters(const StellarKickParmsT& p_KickParameters)         { EnsureSNDetails().kickParameters = p_KickParameters; }

    void ResetSNStellarType() {
        if (m_SupernovaDetails.has_value()) {                                   // SN details exist?
            m_SupernovaDetails->stellarType = {};                               // Yes, reset stellar type
        }                                                                       // No need if no SN details
    }
    void SetSNStellarType(const STELLAR_TYPE p_StellarType)                     { EnsureSNDetails().stellarType = p_StellarType; }

    void ResetTimescales()                                                      { for (SizeT idx = 0; idx < m_Timescales.size(); idx++) m_Timescales[idx] = std::nullopt; }
    void ResetTimescales(const SizeT p_Idx) {                                   // Indexed resetter
        if (p_Idx < m_Timescales.size()) m_Timescales[p_Idx] = std::nullopt;    // Reset value at requested index iff index is in bounds
        else utils::Fail(ERROR::OUT_OF_BOUNDS);                                 // Index out of bounds - fail
    }
    void SetTimescales(const OptDblVectorT& p_Timescales)                       { m_Timescales = p_Timescales; }
    void SetTimescales(const SizeT p_Idx, const double p_Value) {               // Indexed setter
        if (p_Idx < m_Timescales.size()) m_Timescales[p_Idx] = p_Value;         // Set value at requested index iff index is in bounds
        else utils::Fail(ERROR::OUT_OF_BOUNDS);                                 // Index out of bounds - fail
    }

    #undef STATE_SETTER_BOOL
    #undef STATE_SETTER_DBL
    #undef STATE_SETTER


    // Member functions

    /*
     * ResetState
     *
     * @brief
     * Resets all state member variables to the uninitialised state
     * 
     * void ResetState()
     */
    void ResetState() {

        ResetAge();

        ResetAngularFrequency();
        ResetAngularMomentum();

        ResetCNOProcessedCoreMass();
        ResetCOCoreMass();
        ResetCoreMass();
        ResetHeCoreMass();
        
        ResetDMdt();
        
        ResetDominantMassLossType();
        
        ResetDt();

        ResetEnvelopeExpelledByPulsations();

        ResetError();
        ResetEvolutionStatus();
        ResetEvolvedCHE();
        
        ResetGBParams();

        ResetHAbundanceCore();
        ResetHAbundanceSurface();
        ResetHeAbundanceCore();
        ResetHeAbundanceCoreOut();
        ResetHeAbundanceSurface();
        
        ResetLBV_PhaseFlag();

        ResetLuminosity();
        ResetLuminosityZAMSEffective();
        
        ResetMass();
        ResetMassEffectiveInitial();

        ResetMetallicity();

        ResetPulsarDetails();

        ResetRadius();
        ResetRadiusZAMSEffective();
        
        ResetStellarType();

        ResetSNDetails();

        ResetTau();
        
        ResetTemperature();
        ResetTemperatureZAMSEffective();
        
        ResetTime();

        ResetTimescales();
    }
};




// ConstituentState class - state for an individual constituent star (i.e. a star
// that is part of a binary or higher-order system).
// 
// Parallel to StarState (single-star physics state).  ConstituentState carries
// binary-context state - things that don't exist outside the context of a pairing:
// mass-transfer bookkeeping, common-envelope snapshots, RLOF flags, orbital-frame
// velocity, accretion-driven WD outcome flags, etc.
// 
//
// Members:
//   - History of MT events as donor (m_MTDonorHistory)
//   - Common-envelope snapshots (m_CEDetails)
//   - RLOF details (m_RLOFDetails)
//   - Misc flags (m_Flags: recycledNS, ...)
//   - First-MT episode marker
//   - Mass/MT differentials, orbital energy snapshots, tidal omega diff
//   - Orbital-frame velocity (m_ComponentVelocity) - was in BaseStar
//   - Total mass-loss rate (m_TotalMassLossRate) - was in BaseStar; includes MT
//   - WD-specific accretion-driven outcome state (m_HeShellDetonation,
//     m_OffCentreIgnition, m_AccretionRegime, m_LambdaRitter) - was in
//     WhiteDwarfs and derivatives; these are accretion outcomes, not
//     intrinsic WD properties, so they live in the constituent
//
// Members NOT held here:
//   - m_Companion: relationship pointer, set once at construction, lives as
//     a direct member of BinaryConstituentStar (not state, not history-tracked)
class ConstituentState {       
        
private:
        
    // MT donor history
    // Vector of historical MT events as the donor.
    // The vector will always exist - empty vector indicates no events

    std::vector<MTEventT>   m_MTDonorHistory;   // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR - DO WE NEED THIS WITH STATE EVENTS?????


    // CE / RLOF snapshots
    StellarCEDetailsT       m_CEDetails;                            // Common envelope details
    StellarRLOFDetailsT     m_RLOFDetails;                          // Roche-lobe overflow details

    // MT bookkeeping
    bool                    m_FirstMassTransferEpisode;             // Activated for the initial Mass Transfer Episode

    struct FLAGS {                                                  // Miscellaneous flags
        bool recycledNS;                                            // Indicate whether the accretor was a recycled neutron star
    }                       m_Flags;

    double                  m_MassLossDiff;
    double                  m_MassTransferDiff;

    double                  m_OmegaTidesIndividualDiff;

    double                  m_OrbitalEnergyPostSN;
    double                  m_OrbitalEnergyPreSN;

    // Migrated from BaseStar (orbital-frame / MT-aware single-star quantities)
    Vector3d                m_ComponentVelocity;                    // Component velocity in the orbital frame
    double                  m_TotalMassLossRate;                    // Includes wind + MT (per JR's call)

    // WD-specific accretion-driven outcome state
    //
    // These were members of WhiteDwarfs/HeWD/COWD/ONeWD but they're set by
    // binary-context accretion, not intrinsic WD physics.  Moved here so
    // single-star (SSE) WD evolution never reads "binary-set" state.
    //
    bool                    m_HeShellDetonation;                    // Triggered HeSD via accreted He shell
    bool                    m_OffCentreIgnition;                    // Off-centre C ignition (COWD -> ONeWD path)
    ACCRETION_REGIME        m_AccretionRegime;                      // Current accretion regime (WD family)
    bool                    m_ShouldRejuvenate;                     // HeWD rejuvenation triggered (HeWD -> HeMS path)
    bool                    m_IsSubChandrasekharTypeIa;             // HeWD SubCh Type Ia trigger



public:

    // ConstituentState class constructor.
    ConstituentState() :
        m_MTDonorHistory{},
        m_CEDetails{},
        m_RLOFDetails{},
        m_FirstMassTransferEpisode(false),
        m_Flags{false},
        m_MassLossDiff(DEFAULT_INITIAL_DOUBLE_VALUE),
        m_MassTransferDiff(DEFAULT_INITIAL_DOUBLE_VALUE),
        m_OmegaTidesIndividualDiff(DEFAULT_INITIAL_DOUBLE_VALUE),
        m_OrbitalEnergyPostSN(DEFAULT_INITIAL_DOUBLE_VALUE),
        m_OrbitalEnergyPreSN(DEFAULT_INITIAL_DOUBLE_VALUE),
        m_ComponentVelocity{},
        m_TotalMassLossRate(0.0),
        m_HeShellDetonation(false),
        m_OffCentreIgnition(false),
        m_AccretionRegime(ACCRETION_REGIME::ZERO),
        m_ShouldRejuvenate(false),
        m_IsSubChandrasekharTypeIa(false)
    {}
        

    // Setters

    // MT donor event - add to history.
    // Event only added if first event or different from latest donor stellar type.
    void AddMtdonorEvent(const STELLAR_TYPE p_DonorST, const STELLAR_TYPE p_AccretorST, const MT_CASE p_MTcase) { 

        if (m_MTDonorHistory.empty()) {                                                     // first MT donor event?
            m_MTDonorHistory.push_back({ p_DonorST, p_AccretorST, p_MTcase });              // yes - set first event
        }                                                                                   // no - have existing events
        else if (!utils::IsOneOf(p_DonorST, { m_MTDonorHistory.back().donorST })) {         // first MT donor event as donor stellar type?
            m_MTDonorHistory.push_back({ p_DonorST, p_AccretorST, p_MTcase });              // yes - add new event
        }
    }

    void SetCEDetails(const StellarCEDetailsT& p_CEDetails)               { m_CEDetails = p_CEDetails; }
    void SetRLOFDetails(const StellarRLOFDetailsT& p_RLOFDetails)         { m_RLOFDetails = p_RLOFDetails; }
    void SetRecycledNS(const bool p_Value)                                { m_Flags.recycledNS = p_Value; }
    void SetFirstMassTransferEpisode(const bool p_Value)                  { m_FirstMassTransferEpisode = p_Value; }
    void SetMassLossDiff(const double p_Value)                            { m_MassLossDiff = p_Value; }
    void SetMassTransferDiff(const double p_Value)                        { m_MassTransferDiff = p_Value; }
    void SetOmegaTidesIndividualDiff(const double p_Value)                { m_OmegaTidesIndividualDiff = p_Value; }
    void SetOrbitalEnergyPostSN(const double p_Value)                     { m_OrbitalEnergyPostSN = p_Value; }
    void SetOrbitalEnergyPreSN(const double p_Value)                      { m_OrbitalEnergyPreSN = p_Value; }
    void SetComponentVelocity(const Vector3d& p_Value)                    { m_ComponentVelocity = p_Value; }
    void UpdateComponentVelocity(const Vector3d& p_Delta)                 { m_ComponentVelocity += p_Delta; }
    void SetTotalMassLossRate(const double p_Value)                       { m_TotalMassLossRate = p_Value; }
    void SetHeShellDetonation(const bool p_Value)                         { m_HeShellDetonation = p_Value; }
    void SetOffCentreIgnition(const bool p_Value)                         { m_OffCentreIgnition = p_Value; }
    void SetAccretionRegime(const ACCRETION_REGIME p_Value)               { m_AccretionRegime = p_Value; }
    void SetShouldRejuvenate(const bool p_Value)                          { m_ShouldRejuvenate = p_Value; }
    void SetIsSubChandrasekharTypeIa(const bool p_Value)                  { m_IsSubChandrasekharTypeIa = p_Value; }


    // Getters

    StrT                                MTDonorHistoryString(get::unwrappedT = {}) const            { return utils::GetValue(m_MTDonorHistoryString); }
    OptStrT                             MTDonorHistoryString(get::asOptionalT    ) const            { return m_MTDonorHistoryString; }

    // MT donor history
    std::vector<MTEventT> MTDonorHistory() const { return m_MTDonorHistory; }

    MTEventT FirstMTDonorEvent(get::unwrappedT = {}) const { if (m_MTDonorHistory.empty()) utils::Fail(ERROR::NO_VALUE); return m_MTDonorHistory.front(); }
    std::optional<MTEventT> FirstMTDonorEvent(get::asOptionalT) const { if (m_MTDonorHistory.empty()) return std::nullopt; return std::make_optional(m_MTDonorHistory.front()); }

    MTEventT LatestMTDonorEvent(get::unwrappedT = {}) const { if (m_MTDonorHistory.empty()) utils::Fail(ERROR::NO_VALUE); return m_MTDonorHistory.back(); }
    std::optional<MTEventT> LatestMTDonorEvent(get::asOptionalT) const { if (m_MTDonorHistory.empty()) return std::nullopt; return std::make_optional(m_MTDonorHistory.back()); }

    // CE and RLOF compound state - return by const ref to avoid copying the whole struct
    const StellarCEDetailsT&   CEDetails()   const { return m_CEDetails; }
    const StellarRLOFDetailsT& RLOFDetails() const { return m_RLOFDetails; }
    // Mutable overloads for per-field writes (m_CEDetails.preCEE.X = ...)
          StellarCEDetailsT&   CEDetails()         { return m_CEDetails; }
          StellarRLOFDetailsT& RLOFDetails()       { return m_RLOFDetails; }

    // Misc flags
    bool   RecycledNS()              const { return m_Flags.recycledNS; }

    // MT bookkeeping
    bool   FirstMassTransferEpisode() const { return m_FirstMassTransferEpisode; }
    double MassLossDiff()            const { return m_MassLossDiff; }
    double MassTransferDiff()        const { return m_MassTransferDiff; }
    double OmegaTidesIndividualDiff() const { return m_OmegaTidesIndividualDiff; }

    // Orbital energy snapshots
    double OrbitalEnergyPostSN()     const { return m_OrbitalEnergyPostSN; }
    double OrbitalEnergyPreSN()      const { return m_OrbitalEnergyPreSN; }

    // Migrated from BaseStar
    const Vector3d& ComponentVelocity() const { return m_ComponentVelocity; }
    double  Speed()       const { return m_ComponentVelocity.Magnitude(); }
    double  VelocityX()   const { return m_ComponentVelocity.xValue(); }
    double  VelocityY()   const { return m_ComponentVelocity.yValue(); }
    double  VelocityZ()   const { return m_ComponentVelocity.zValue(); }
    double  TotalMassLossRate() const { return m_TotalMassLossRate; }

    // WD-specific accretion outcomes
    bool    HeShellDetonation()       const { return m_HeShellDetonation; }
    bool    OffCentreIgnition()       const { return m_OffCentreIgnition; }
    ACCRETION_REGIME AccretionRegime() const { return m_AccretionRegime; }
    bool    ShouldRejuvenate()        const { return m_ShouldRejuvenate; }
    bool    IsSubChandrasekharTypeIa() const { return m_IsSubChandrasekharTypeIa; }


    // Reset
    void ResetMTDonorHistory()             { m_MTDonorHistory.clear(); }
    void ResetCEDetails()                  { m_CEDetails = StellarCEDetailsT{}; }
    void ResetRLOFDetails()                { m_RLOFDetails = StellarRLOFDetailsT{}; }
    void ResetFlags()                      { m_Flags.recycledNS = false; }
    void ResetFirstMassTransferEpisode()   { m_FirstMassTransferEpisode = false; }
    void ResetMassLossDiff()               { m_MassLossDiff = DEFAULT_INITIAL_DOUBLE_VALUE; }
    void ResetMassTransferDiff()           { m_MassTransferDiff = DEFAULT_INITIAL_DOUBLE_VALUE; }
    void ResetOmegaTidesIndividualDiff()   { m_OmegaTidesIndividualDiff = DEFAULT_INITIAL_DOUBLE_VALUE; }
    void ResetOrbitalEnergyPostSN()        { m_OrbitalEnergyPostSN = DEFAULT_INITIAL_DOUBLE_VALUE; }
    void ResetOrbitalEnergyPreSN()         { m_OrbitalEnergyPreSN = DEFAULT_INITIAL_DOUBLE_VALUE; }
    void ResetComponentVelocity()          { m_ComponentVelocity = Vector3d{}; }
    void ResetTotalMassLossRate()          { m_TotalMassLossRate = 0.0; }
    void ResetHeShellDetonation()          { m_HeShellDetonation = false; }
    void ResetOffCentreIgnition()          { m_OffCentreIgnition = false; }
    void ResetAccretionRegime()            { m_AccretionRegime = ACCRETION_REGIME::ZERO; }
    void ResetShouldRejuvenate()           { m_ShouldRejuvenate = false; }
    void ResetIsSubChandrasekharTypeIa()   { m_IsSubChandrasekharTypeIa = false; }

    void ResetHasMTDonorHistory()                                               { m_HasMTDonorHistory = std::nullopt; }
    void SetHasMTDonorHistory(const bool p_HasHistory)                          { m_HasMTDonorHistory = p_HasHistory; }
    void ResetMTDonorHistoryString()                                            { m_MTDonorHistoryString = std::nullopt; }
    void SetMTDonorHistoryString(const StrT& p_Str)                             { m_MTDonorHistoryString = p_Str; }

    void ResetState() {
        ResetMTDonorHistory();

        ResetHasMTDonorHistory();
        ResetMTDonorHistoryString();
        ResetCEDetails();
        ResetRLOFDetails();
        ResetFlags();
        ResetFirstMassTransferEpisode();
        ResetMassLossDiff();
        ResetMassTransferDiff();
        ResetOmegaTidesIndividualDiff();
        ResetOrbitalEnergyPostSN();
        ResetOrbitalEnergyPreSN();
        ResetComponentVelocity();
        ResetTotalMassLossRate();
        ResetHeShellDetonation();
        ResetOffCentreIgnition();
        ResetAccretionRegime();
        ResetShouldRejuvenate();
        ResetIsSubChandrasekharTypeIa();
    }

};


// StellarBinaryState class - state for stellar binary systems
class StellarBinaryState {       
        
private:
        



};






// ===========================================================================
// State history storage
//
// A system's state history is logically an append-only sequence with a fixed
// start state (BOS) at index 0 and the current state (TOS) at the top. It is
// NOT a linked list of individually-allocated nodes: it is backed by fixed-size,
// address-stable chunks drawn from a shared pool. Existing chunks never move, so
// raw pointers into their slots stay valid for the life of the borrow - which is
// what the event channels rely on.
// ===========================================================================

// A StateChunk is a fixed run of state slots. The slots are default-constructed
// once (when the chunk is allocated) and reused via assignment thereafter.
template <typename T>
struct StateChunk {
    static constexpr SizeT SIZE = STATE_CHUNK_SIZE;
    T slots[SIZE];
};


// StateStorage - per-T singleton pool of reusable chunks.
//
// Stars/binaries are new-ed and delete-d per system, but their storage is
// recycled across the whole run: one malloc per chunk EVER, reused thereafter.
// Acquire() lends a chunk (allocating only when the free list is empty);
// Release() returns it to the free list (never frees mid-run). Every chunk ever
// allocated is tracked in m_All and freed at program exit (track-and-free).
//
// Single-threaded by design - COMPAS scales by running multiple processes, not
// threads - so a plain static singleton is safe.
template <typename T>
class StateStorage {

    std::vector<StateChunk<T>*> m_All;    // every chunk ever allocated (owns them)
    std::vector<StateChunk<T>*> m_Free;   // chunks currently available to lend

    StateStorage() = default;
    StateStorage(const StateStorage&)            = delete;
    StateStorage& operator=(const StateStorage&) = delete;

public:

    static StateStorage& Instance() { static StateStorage self; return self; }

    StateChunk<T>* Acquire() {
        if (!m_Free.empty()) { StateChunk<T>* c = m_Free.back(); m_Free.pop_back(); return c; }
        StateChunk<T>* c = new StateChunk<T>();                                             // the only place we allocate
        m_All.push_back(c);
        return c;
    }

    void Release(StateChunk<T>* p_Chunk) { m_Free.push_back(p_Chunk); }                     // reclaim - never frees mid-run

    ~StateStorage() { for (StateChunk<T>* c : m_All) delete c; }                            // track-and-free at program exit
};


template <typename T>
class StateHistory {

public:

    // Detection predicate: given the incoming state and the previous state
    // (prev is nullptr on the first push), return true if this push is an
    // event of this channel's type.
    using EventPredicate = std::function<bool(const T& p_New, const T* p_Prev)>;

private:

    struct EventChannel {
        EventType       type;
        EventPredicate  predicate;
        std::vector<T*> entries;         // pointers into chunk slots - stable; empty until first matching event
    };

    std::vector<StateChunk<T>*> m_Chunks;                                    // borrowed backing chunks, in logical order
    SizeT                       m_Size     = 0;                              // number of live states (BOS at logical 0)
    SizeT                       m_Capacity = DEFAULT_STATE_HISTORY_CAPACITY; // fail-stop ceiling (runaway guard); 0 = unbounded
    std::vector<EventChannel>   m_Channels;

    // Address of the state at logical index i. O(1). Never call with i >= m_Size.
    T&       At(SizeT i)       { return m_Chunks[i / StateChunk<T>::SIZE]->slots[i % StateChunk<T>::SIZE]; }
    const T& At(SizeT i) const { return m_Chunks[i / StateChunk<T>::SIZE]->slots[i % StateChunk<T>::SIZE]; }

    // Find a channel by type (few channels -> linear scan is fine)
    const EventChannel* Channel(EventType p_Type) const {
        for (const auto& ch : m_Channels) if (ch.type == p_Type) return &ch;
        return nullptr;
    }

public:

    // Read-only forward iterator over the history (BOS at 0 -> TOS at m_Size-1).
    class ConstIter {
        const StateHistory* m_H;
        SizeT               m_Idx;
    public:
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type        = T;
        using reference         = const T&;
        using pointer           = const T*;
        using difference_type   = std::ptrdiff_t;

        ConstIter(const StateHistory* p_H, SizeT p_Idx) : m_H(p_H), m_Idx(p_Idx) {}
        const T& operator*()  const { return m_H->At(m_Idx); }
        const T* operator->() const { return &m_H->At(m_Idx); }
        ConstIter& operator++() { ++m_Idx; return *this; }
        ConstIter& operator--() { --m_Idx; return *this; }
        bool operator==(const ConstIter& o) const { return m_H == o.m_H && m_Idx == o.m_Idx; }
        bool operator!=(const ConstIter& o) const { return !(*this == o); }
    };

    // Read-only reverse iterator (TOS -> BOS). rEnd() uses SizeT(-1) as the
    // "one before BOS" sentinel: decrementing the index at 0 wraps to SizeT(-1).
    class ConstRevIter {
        const StateHistory* m_H;
        SizeT               m_Idx;
    public:
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type        = T;
        using reference         = const T&;
        using pointer           = const T*;
        using difference_type   = std::ptrdiff_t;

        ConstRevIter(const StateHistory* p_H, SizeT p_Idx) : m_H(p_H), m_Idx(p_Idx) {}
        const T& operator*()  const { return m_H->At(m_Idx); }
        const T* operator->() const { return &m_H->At(m_Idx); }
        ConstRevIter& operator++() { --m_Idx; return *this; }                               // reverse: advance == step back
        ConstRevIter& operator--() { ++m_Idx; return *this; }
        bool operator==(const ConstRevIter& o) const { return m_H == o.m_H && m_Idx == o.m_Idx; }
        bool operator!=(const ConstRevIter& o) const { return !(*this == o); }
    };


    void RegisterEvent(EventType p_Type, EventPredicate p_Predicate) {
        m_Channels.push_back({p_Type, std::move(p_Predicate), {}});
    }


    // Constructors
    //
    // Every history borrows a warm chunk from the pool and seeds a single state
    // at BOS (logical 0), so the history is never empty. On destruction all
    // borrowed chunks are returned to the pool (no frees) - which is how storage
    // recycles across systems (the star is new-ed/delete-d per system, its
    // StateHistory member with it).

    StateHistory() {
        m_Chunks.push_back(StateStorage<T>::Instance().Acquire());
        At(0) = T{};                                                                        // seed one blank start state
        m_Size = 1;
    }

    StateHistory(const T& p_StartState, const SizeT p_Capacity = DEFAULT_STATE_HISTORY_CAPACITY)
        : m_Capacity(p_Capacity) {                                                          // 0 = unbounded; otherwise fail-stop ceiling
        m_Chunks.push_back(StateStorage<T>::Instance().Acquire());
        At(0) = p_StartState;
        m_Size = 1;
    }

    StateHistory(const StateHistory& o) : m_Capacity(o.m_Capacity) {
        // Deep-copy states into freshly-borrowed chunks, recording where each
        // source state landed so event-channel entries can be remapped.
        std::unordered_map<const T*, T*> remap;
        m_Chunks.push_back(StateStorage<T>::Instance().Acquire());
        for (SizeT i = 0; i < o.m_Size; ++i) {
            if (i == m_Chunks.size() * StateChunk<T>::SIZE)                                 // current chunk full - borrow another
                m_Chunks.push_back(StateStorage<T>::Instance().Acquire());
            At(i) = o.At(i);
            remap.emplace(&o.At(i), &At(i));
        }
        m_Size = o.m_Size;

        // Copy channel definitions, translating entries to the new slots.
        for (const auto& ch : o.m_Channels) {
            EventChannel nc{ ch.type, ch.predicate, {} };                                  // predicate is std::function - copyable
            nc.entries.reserve(ch.entries.size());
            for (const T* e : ch.entries) nc.entries.push_back(remap.at(e));               // source slot ptr -> new slot ptr
            m_Channels.push_back(std::move(nc));
        }
    }

    ~StateHistory() {
        for (StateChunk<T>* c : m_Chunks) StateStorage<T>::Instance().Release(c);          // return storage to the pool
    }

    friend void swap(StateHistory& a, StateHistory& b) noexcept {
        using std::swap;
        swap(a.m_Chunks,   b.m_Chunks);
        swap(a.m_Size,     b.m_Size);
        swap(a.m_Capacity, b.m_Capacity);
        swap(a.m_Channels, b.m_Channels);
    }

    StateHistory& operator=(StateHistory o) {   // BY VALUE - copy-and-swap (o built by copy or move ctor)
        swap(*this, o);
        return *this;
    }

    StateHistory(StateHistory&&) noexcept = default;


    // SetStartState - overwrite the start (BOS) state in place. Intended for an
    // object's constructor: build the initial state in a buffer, then commit it
    // here. The history is never empty (a start state is always seeded), so BOS
    // always exists.
    void SetStartState(const T& p_State) { At(0) = p_State; }


    // Push - append p_State as the new TOS.
    //
    // Pure append: no eviction, no window. If a fail-stop capacity is set and
    // reached, the system is aborted (THROW_ERROR_STATIC throws an int caught
    // up-stack; the run continues with the next system). Storage grows by
    // borrowing another chunk from the pool only when the current one is full.
    void Push(const T& p_State) {
        if (m_Capacity != 0 && m_Size >= m_Capacity)
            THROW_ERROR_STATIC(ERROR::STATE_HISTORY_CAPACITY_EXCEEDED);                     // runaway guard - aborts this system

        if (m_Size == m_Chunks.size() * StateChunk<T>::SIZE)                               // current chunk full?
            m_Chunks.push_back(StateStorage<T>::Instance().Acquire());                     // borrow another

        const T* prev = &At(m_Size - 1);                                                   // current TOS (history never empty)
        T*       slot = &At(m_Size);
        *slot = p_State;                                                                   // one copy into recycled storage
        ++m_Size;

        for (auto& ch : m_Channels)                                                        // event detection (append-only)
            if (ch.predicate(*slot, prev)) ch.entries.push_back(slot);
    }

    void AddState(const T& p_State) { Push(p_State); }                                     // AddState() is an alias for Push()


    // Pop - remove the TOS state and return a copy of it. The start state (BOS)
    // is never removed. Returns std::nullopt if only BOS remains. The trailing
    // chunk is retained (not returned to the pool) for reuse by the next Push.
    std::optional<T> Pop() {
        if (m_Size <= 1) return std::nullopt;                                              // never remove BOS

        std::optional<T> retVal = At(m_Size - 1);                                          // copy TOS for caller
        T* victim = &At(m_Size - 1);
        for (auto& ch : m_Channels)                                                        // scrub entries pointing at the removed slot
            ch.entries.erase(std::remove(ch.entries.begin(), ch.entries.end(), victim), ch.entries.end());
        --m_Size;

        return retVal;
    }

    std::optional<T> RemoveState() { return Pop(); }                                       // RemoveState() is an alias for Pop()


    // Peek - return a copy of the state p_Count elements below TOS (0 = TOS).
    // Returns std::nullopt if the requested element does not exist. O(1).
    std::optional<T> Peek(const SizeT p_Count = 0) const {
        if (p_Count < m_Size) return At(m_Size - 1 - p_Count);
        return std::nullopt;
    }

    std::optional<T> LookBack(const SizeT p_Count = 0) const { return Peek(p_Count); }     // LookBack() is an alias for Peek()


    // Bookkeeping
    SizeT Size() const  { return m_Size; }
    SizeT Count() const { return Size(); }


    // Iterators - read only - forward and reverse
    ConstIter    Begin()  const { return ConstIter(this, 0); }
    ConstIter    End()    const { return ConstIter(this, m_Size); }
    ConstRevIter rBegin() const { return ConstRevIter(this, m_Size - 1); }
    ConstRevIter rEnd()   const { return ConstRevIter(this, static_cast<SizeT>(-1)); }


    // Sentinel getters
    const T& TOS() const { return At(m_Size - 1); }                      // Top of stack - always exists (history never empty)
    const T& BOS() const { return At(0); }                              // Bottom of stack - always exists (history never empty)


    // Functions to determine if defined states exist
    // StartState and CurrentState are guaranteed to exist (history never empty)
    bool HaveZAMSState() const {
        const auto st = StartState().StellarType(get::asOptionalT{});
        return st.has_value() && utils::IsOneOf(*st, MAIN_SEQUENCE);
    }
    bool HavePreviousState() const { return m_Size > 1; }


    // By-reference sentinel accessors
    //
    // These return a const reference directly to the stored state, for callers
    // that read a field and discard the reference within the same expression
    // (e.g. wrapping getters: 'return CurrentState().Mass();'). This avoids
    // copying a whole StarState just to read one field.
    //
    // IMPORTANT: do NOT store the returned reference. It refers to a slot inside
    // a backing chunk and is invalidated by a Pop that removes that element (or
    // by the history's destruction). Transient use within a single expression is
    // safe; holding onto it across a state transition is not.
    //
    // CurrentState()/StartState() always have a value (history never empty).
    // ZAMSState()/PreviousState() will call utils::Fail() (which terminates) if
    // the requested state does not exist.

    const T& CurrentState() const { return TOS(); }                      // TOS - always exists

    const T& StartState()   const { return BOS(); }                      // BOS - always exists

    const T& ZAMSState() const {
        if (!HaveZAMSState()) utils::Fail(ERROR::NO_STATE);
        return BOS();
    }

    const T& PreviousState() const {
        if (!HavePreviousState()) utils::Fail(ERROR::NO_STATE);
        return At(m_Size - 2);
    }


    // Event queries
    std::size_t EventCount(EventType p_Type) const {
        const EventChannel* ch = Channel(p_Type);
        return ch ? ch->entries.size() : 0;
    }

    const T* LastEvent(EventType p_Type) const {
        const EventChannel* ch = Channel(p_Type);
        return (ch && !ch->entries.empty()) ? ch->entries.back() : nullptr;
    }

    template <typename Fn> void ForEachEvent(EventType p_Type, Fn&& p_Fn) const {
        const EventChannel* ch = Channel(p_Type);
        if (ch) for (const T* p : ch->entries) p_Fn(*p);                                   // p -> const T&
    }
};
