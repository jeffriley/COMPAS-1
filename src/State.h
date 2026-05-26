#pragma once

#include <list>
#include <optional>
#include <variant>

#include "constants.h"

//#include "utils.h"

constexpr SizeT DEFAULT_STATE_HISTORY_STACK_SIZE = 0; 

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
an object (of type T, and the object contains a copy of the element),rather than a
pointer to the actual stack element.an object.  Returning a pointer (even a const
pointer to const data) to the actual stack element would cause problems if the stack
element was popped or removed due to a push when stack full (the pointer held by the
caller would no longer be valid), so a copy of the stack element is returned).

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

      elem.value().RandomSeed()
              or
      (*elem).RandomSeed()
              or
      elem->RandomSeed()

      to get the value of (e.g.) the randomseed variable from the stack element

The construct 'elem.value_or(default_value)' is also available - this returns the
value "default_value" in the case that elem has no value assigned.

Note that attempting to query the value of a std::optional object that has not had a
value assigned will result in undefined bahaviour.

---*/



    /*
     * Fail
     *
     * @brief
     * Prints an error message, shows a stack trace, and terminates the program.
     *
     * 
     * void Fail(const ERROR p_Error, const StrT p_Preamble = "") const
     *  
     * @param       p_Error                         The error to be printed (index into error catalog)
     * @param       p_Preamble                      Preamble string to be printed prior to the error string
     *                                              This replaces the default preamble ("\nError retrieving state value\n")
     */
    void Fail(const ERROR p_Error, const StrT p_Preamble = "") const {

        StrT errStr = "";

        if (p_Preamble.empty()) {                               // Preamble supplied?
            errStr = "\nError retrieving state value\n";        // No - default preamble
        }
        else {                                                  // Yes
            if (p_Preamble[0] != '\n') errStr += "\n";          // Add leading newline if necessary
            errStr += p_Preamble;                               // Add supplied preamble
            if (errStr.back() != '\n') errStr += "\n";          // Add trailing newline if necessary
        }

        // Look for error in error catalog
        std::unordered_map<ERROR, std::tuple<ERROR_SCOPE, std::string>>::const_iterator iter = ERROR_CATALOG.find(p_Error);
        if (iter != ERROR_CATALOG.end()) {                      // Found?
            errStr += std::get<1>(staticIter->second) + "\n";   // Yes - add message to error string
        }
        else {                                                  // No
            errStr += "Uknown error\n";                         // Add generic message
        }

        std::cerr << errStr;                                    // Announce error
        utils::ShowStackTrace();                                // Show stack trace
        std::exit(1);                                           // Terminate program
    }


    /*
     * GetValue
     *
     * @brief
     * Gets the underlying value of a variable wrapped in std::optional<>.
     * If the valiable passed does not have an assigned value, this function calls
     * Fail() (see below) to terminate the program.
     * 
     * 
     * template<typename T>
     * T GetValue(const std::optional<T>& p_Var)
     */
    template<typename T>
    T GetValue(const std::optional<T>& p_Var) {
        if (p_Var.has_value()) return p_Var.value();            // If variable has a value, return it unwrapped
        else Fail(ERROR::NO_VALUE);                             // Otherwise fail
    }



// StarState class - state for single stars
class StarState {       

private:
        
    OptBoolT                     m_CHE; // CHE flag - true if the star spent entire MS as a CH star; false if evolved CH->MS

    OptDblT                   m_CNOProcessedCoreMass;   // Only applies to main sequence stars
    OptDblT                   m_COCoreMass;
    OptDblT                   m_CoreMass;
        
    OptDblT                   m_dMdt;
        
    std::optional<MASS_LOSS_TYPE>           m_DominantMassLossType;
        
    OptDblT                   m_dt;
        
    std::optional<ERROR>                    m_Error;
    std::optional<EVOLUTION_STATUS>         m_EvolutionStatus;
        
    OptDblT                   m_HAbundanceCore;
    OptDblT                   m_HAbundanceSurface;
    OptDblT                   m_HeAbundanceCore;
    OptDblT                   m_HeAbundanceSurface;

    OptDblT                   m_HeCoreMass;

    // For timescales and GB parameters vectors (see Hurley et al. 2000), each
    // vector is constructed here with the correct number of elements, and each
    // element is wrapped in std::optional<>, and initialised to std::nullopt.
    // The vector is guaranteed to exist, but each element is optional.

    OptDblVectorT                           m_HurleyTimescales;
    OptDblVectorT                           m_HurleyGBparams;

    OptDblT                   m_Luminosity;
    OptDblT                   m_LuminosityZAMSEffective;
        
    OptDblT                   m_Mass;
    OptDblT                   m_MassEffectiveInitial; // Effective initial mass of the star (Msol) (Hurley et al., 2000)

    OptDblT                   m_OmegaZAMS;          // Angular frequency at ZAMS

    std::optional<PulsarDetailsT>           m_PulsarDetails;
        
    OptDblT                   m_Radius;
    OptDblT                   m_RadiusZAMSEffective;
        
    std::optional<ULongT>                   m_RandomSeed;
        
    std::optional<STELLAR_TYPE>             m_StellarType;
        
    std::optional<StellarSupernovaDetailsT> m_SupernovaDetails;
        
    OptDblT                   m_Tau;
        
    OptDblT                   m_Temperature;
        
    OptDblT                   m_Time;       

        
public:

    // StarState class
    // constructor
    StarState() {

        // Create GB parameters and timescales vectors.
        // These vectors are guaranteed to exist, but each element is wrapped in std::optional<>.
        // Construct vectors of required sizes and initialise each element to std::nullopt.
        for (SizeT idx = 0; idx < utils::to_underlying(HURLEY_GBP::COUNT); idx++) m_HurleyGBparams.push_back(std::nullopt);
        for (SizeT idx = 0; idx < utils::to_underlying(TIMESCALE::COUNT); idx++) m_HurleyTimescales.push_back(std::nullopt);
    }
        



               
        


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
    // the value of the variable is returned. If the variable *is* wrapped in std::optional<>, then:
    //    (a) if the variable has a value, the value will be returned, but
    //    (b) if the variable *does not* have a value, an error message and stacktrace will be
    //        displayed, and the program terminated.
    // The reasoning for (b) is that for most circumstances, developers will know when state
    // variables have, or are expected to have, defined values, so can just call the getters
    // and get the value of the variable - and if the variable does not have a value that's 
    // probably either a coding problem, or a corruption/timing issue, that needs to be addressed.
    //
    // For the non-default (optional) form (see below for calling syntax), the variable requested
    // is returned wrapped in std::optional<>.  The caller should then check via *.has_value() if
    // the returned variable has a value assigned.  If the variable being requested is *not* actually
    // declared as std::optional<>, the getter will wrap it in std::optional<> before returning it.
    // Such variables will *always* have a value when returned. Variables that *are* declared as
    // std::optional<> will be returned as is, and may or may not have a value assigned.
    //
    //
    // Vector variables
    // ----------------
    //
    // Some state variables are stored as vectors - e.g. the Hurley 'a' and 'b' coefficients.
    // All vector variables can be retrieved as a vector, or as individual elements, but there
    // is overhead in retrieving the individual elements - range checking is performed on the
    // specified index (where appropriate).  In some cases it will be more appropriate to retrieve
    // the vector, and let the caller do the range checking as approriate
    // (the Hurley a and b coefficients, for example).  See below for calling syntax.
    //
    // When a state vector variable is retrieved as a vector, the default form of the getter will
    // return the vector as type 'std::vector<T>', where 'T' is the type of the underlying variable
    // (i.e. double, int, etc.), whereas the non-default (optional) form of the getter will return
    // the vector as type 'std::optional<std::vector<T>>', where 'T' is the type of the underlying
    // variable.
    // 
    // To facilitate retrieval of single, indexed, elements of vector variables, an alternative
    // form of the getters for vector variables is provided.  This alternative form takes a
    // parameter that is the index into the vector of the single element required (see below).
    // For getters of this form, the single element of the vector is returned as type 'T' (the
    // underlying type, i.e. double, int, etc.) for the default form of the getter, and as type
    // 'std::optional<T>' for the non-default (optional) form of the getter.
    //
    //
    // Usage and calling syntax
    // ------------------------
    //
    // The namespace "get" is declared below to facilitate calling the different getter variants.
    //
    // If the developer expects a state variable to have a value assigned (whether it be a
    // variable declated as std::optional<> or not), they should use the default getter variant.
    // For example:
    //
    //     double mass = state.Mass()
    //
    // will return the value of mass as a non-optional double value.
    // The State class member variable m_Mass is declared as std::optional<double>, and there
    // is a window as COMPAS is initialising where it may not have yet been assigned a value.
    // If the call to state.Mass() is made before the value of m_Mass has been
    // set, the call above will fail as described above (error message, stack trace, and program
    // terminated).  If, however, the call is made after the value of m_Mass has been set, the call
    // above will return the value of m_Mass as a variable of type double.
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
    // On the other hand, if the developer is uncertain if a state variable has a value assigned,
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
    // variant of the getter ("get::asOptionalT" must appear as the last parameter in the parameter
    // list).
    //
    // For getters that take a parameter (e.g. getters that return a specific element from a vector),
    // the call syntax is e.g.:
    //
    //    double McBGB = GBparams(HURLEY_GBP::MCBGB)
    //
    // for the default getter variant, getting the element for HURLEY_GBP::MCBGB, and
    //
    //    std::optional<double> McBGB = GBparams(HURLEY_GBP::MCBGB, get::asOptionalT)
    //
    // for the optional getter variant, getting the element for HURLEY_GBP::MCBGB.


    
    DblVectorT GBparams(get::unwrappedT = {}) const { return GetValue(m_HurleyGBparams); }
    Optional_DblVecT GBparams(get::asOptionalT) const { return m_HurleyGBparams; }
    double GBparams(const HURLEY_GBP p_GBParam, get::unwrappedT = {}) const {
        if (m_HurleyGBparams.has_value()) {                                                         // Vector exists?
                return m_HurleyGBparams[static_cast<int>(p_GBParam)];                               // Yes - return unwrapped value
        }
        else Fail(ERROR::NO_VALUE);                                                                 // No - vector does not exist - fail
    }
    OptDblT GBparams(const HURLEY_GBP p_GBParam, get::asOptionalT) const {
        if (m_HurleyGBparams.has_value()) {                                                         // Vector exists?
                return std::optional<double>(m_HurleyGBparams[static_cast<int>(p_GBParam)]);        // Yes - return wrapped value
        }
        else return std::nullopt;                                                                   // No - vector does not exist - return no value
    }

    
    DblVectorT Timescales(get::unwrappedT = {}) const { return GetValue(m_HurleyTimescales); }
    Optional_DblVecT Timescales(get::asOptionalT) const { return m_HurleyTimescales; }
    double Timescales(const TIMESCALE p_Timescale, get::unwrappedT = {}) const {
        if (m_HurleyTimescales.has_value()) {                                                       // Vector exists?
                return m_HurleyTimescales[static_cast<int>(p_Timescale)];                           // Yes - return unwrapped value
        }
        else Fail(ERROR::NO_VALUE);                                                                 // No - vector does not exist - fail
    }
    OptDblT Timescales(const TIMESCALE p_Timescale, get::asOptionalT) const {
        if (m_HurleyTimescales.has_value()) {                                                       // Vector exists?
                return std::optional<double>(m_HurleyTimescales[static_cast<int>(p_Timescale)]);    // Yes - return wrapped value
        }
        else return std::nullopt;                                                                   // No - vector does not exist - return no value
    }


    double   CHE(get::unwrappedT = {}) const { return GetValue(m_CHE); }
    OptBoolT CHE(get::asOptionalT    ) const { return m_CHE; }

    double   CNOProcessedCoreMass(get::unwrappedT = {}) const { return GetValue(m_CNOProcessedCoreMass); }
    OptDblT  CNOProcessedCoreMass(get::asOptionalT    ) const { return m_CNOProcessedCoreMass; }

    double   COCoreMass(get::unwrappedT = {}) const { return GetValue(m_COCoreMass); }
    OptDblT  COCoreMass(get::asOptionalT    ) const { return m_COCoreMass; }

    double   CoreMass(get::unwrappedT = {}) const { return GetValue(m_COCoreMass); }
    OptDblT  CoreMass(get::asOptionalT    ) const { return m_COCoreMass; }

    double   dMdt(get::unwrappedT = {}) const { return GetValue(m_dMdt); }
    OptDblT  dMdt(get::asOptionalT    ) const { return m_dMdt; }
        
    MASS_LOSS_TYPE DominantMassLossType(get::unwrappedT = {}) const { return GetValue(m_DominantMassLossType); }
    std::optional<MASS_LOSS_TYPE> DominantMassLossType(get::asOptionalT) const { return m_DominantMassLossType; }

    double   dt(get::unwrappedT = {}) const { return GetValue(m_dt); }
    OptDblT  dt(get::asOptionalT    ) const { return m_dt; }
 
    ERROR Error(get::unwrappedT = {}) const { return GetValue(m_Error); }
    std::optional<ERROR> Error(get::asOptionalT) const { return m_Error; }
 
    EVOLUTION_STATUS EvolutionStatus(get::unwrappedT = {}) const { return GetValue(m_EvolutionStatus); }
    std::optional<EVOLUTION_STATUS> EvolutionStatus(get::asOptionalT) const { return m_EvolutionStatus; }

    double   HAbundanceCore(get::unwrappedT = {}    ) const { return GetValue(m_HAbundanceCore); }
    OptDblT  HAbundanceCore(get::asOptionalT        ) const { return m_HAbundanceCore; }
    double   HAbundanceSurface(get::unwrappedT = {} ) const { return GetValue(m_HAbundanceSurface); }
    OptDblT  HAbundanceSurface(get::asOptionalT     ) const { return m_HAbundanceSurface; }
    double   HeAbundanceCore(get::unwrappedT = {}   ) const { return GetValue(m_HeAbundanceCore); }
    OptDblT  HeAbundanceCore(get::asOptionalT       ) const { return m_HeAbundanceCore; }
    double   HeAbundanceSurface(get::unwrappedT = {}) const { return GetValue(m_HeAbundanceSurface); }
    OptDblT  HeAbundanceSurface(get::asOptionalT    ) const { return m_HeAbundanceSurface; }

    double   HeCoreMass(get::unwrappedT = {}) const { return GetValue(m_HeCoreMass); }
    OptDblT  HeCoreMass(get::asOptionalT    ) const { return m_HeCoreMass; }
        
    double   Luminosity(get::unwrappedT = {}         ) const { return GetValue(m_Luminosity); }
    OptDblT  Luminosity(get::asOptionalT             ) const { return m_Luminosity; }
    double   Luminosity0(get::unwrappedT = {}        ) const { return LuminosityZAMSEffective(); }
    OptDblT  Luminosity0(get::asOptionalT            ) const { return LuminosityZAMSEffective(get::asOptionalT); }
    double   LuminosityZAMSEffective(get::unwrappedT = {}) const { return GetValue(m_LuminosityZAMSEffective); }
    OptDblT  LuminosityZAMSEffective(get::asOptionalT    ) const { return m_LuminosityZAMSEffective; }
        
    double   Mass(get::unwrappedT = {}         ) const { return GetValue(m_Mass); }
    OptDblT  Mass(get::asOptionalT             ) const { return m_Mass; }
    double   Mass0(get::unwrappedT = {}        ) const { return MassEffectiveInitial(); }
    OptDblT  Mass0(get::asOptionalT            ) const { return MassEffectiveInitial(get::asOptional); }
    double   MassEffectiveInitial(get::unwrappedT = {}) const { return GetValue(m_MassEffectiveInitial); }
    OptDblT  MassEffectiveInitial(get::asOptional     ) const { return m_MassEffectiveInitial; }

    double   OmegaZAMS(get::unwrappedT = {}) const { return GetValue(m_OmegaZAMS); }
    OptDblT  OmegaZAMS(get::asOptionalT    ) const { return m_OmegaZAMS; }
        
    PulsarDetailsT PulsarDetails(get::unwrappedT = {}) const { return GetValue(m_PulsarDetails); }
    std::optional<PulsarDetailsT> PulsarDetails(get::asOptional) const { return m_PulsarDetails; }

    double   Radius(get::unwrappedT = {}         ) const { return GetValue(m_Radius); }
    OptDblT  Radius(get::asOptional              ) const { return m_Radius; }
    double   Radius0(get::unwrappedT = {}        ) const { return RadiusZAMSEffective(); }
    OptDblT  Radius0(get::asOptional             ) const { return RadiusZAMSEffective(get::asOptional); }
    double   RadiusZAMSEffective(get::unwrappedT = {}) const { return GetValue(m_RadiusZAMSEffective); }
    OptDblT  RadiusZAMSEffective(get::asOptional     ) const { return m_RadiusZAMSEffective; }
        
    ULongT RandomSeed(get::unwrappedT = {}) const { return GetValue(m_RandomSeed); }
    std::optional<ULongT> RandomSeed(get::asOptional) return const { m_RandomSeed; }
        
    STELLAR_TYPE StellarType(get::unwrappedT = {}) const { return GetValue(m_StellarType); }
    std::optional<STELLAR_TYPE> StellarType(get::asOptional) return const { m_StellarType; }

    StellarKickParmsT CurrentKickParameters(get::unwrappedT = {}) const { return GetValue(m_SupernovaDetails.currentKickParams); }
    std::optional<StellarKickParmsT> CurrentKickParameters(get::asOptional) const { return m_SupernovaDetails.currentKickParams; }

    StellarKickParmsT SuppliedKickParameters(get::unwrappedT = {}) const { return GetValue(m_SupernovaDetails.suppliedKickParams); }
    std::optional<StellarKickParmsT> SuppliedKickParameters(get::asOptional) const { return m_SupernovaDetails.suppliedKickParams; }

    StellarSupernovaDetailsT SupernovaDetails(get::unwrappedT = {}) const { return GetValue(m_SupernovaDetails); }
    std::optional<StellarSupernovaDetailsT> SupernovaDetails(get::asOptional) const { return m_SupernovaDetails; }

    double   Tau(get::unwrappedT = {}) const { return GetValue(m_Tau); }
    OptDblT  Tau(get::asOptional     ) const { return m_Tau; }
        
    double   Temperature(get::unwrappedT = {}) const { return GetValue(m_Temperature); }
    OptDblT  Temperature(get::asOptional     ) const { return m_Temperature; }
        
    double   Time(get::unwrappedT = {}) const { return GetValue(m_Time); }   
    OptDblT  Time(get::asOptional     ) const { return m_Time; }   



    // setters
        
    void ResetCHE()                                                                 { m_CHE = std::nullopt; }
    void SetCHE(const bool p_CHE)                                                   { m_CHE = p_CHE; }

    void ResetCNOProcessedCoreMass()                                                { m_CNOProcessedCoreMass = std::nullopt; }
    void SetCNOProcessedCoreMass(const double p_CNOProcessedCoreMass)               { m_CNOProcessedCoreMass = p_CNOProcessedCoreMass; }    
    void ResetCOCoreMass()                                                          { m_COCoreMass = std::nullopt; }
    void SetCOCoreMass(const double p_COCoreMass)                                   { m_COCoreMass = p_COCoreMass; }
    void ResetCoreMass()                                                            { m_CoreMass = std::nullopt; }
    void SetCoreMass(const double p_CoreMass)                                       { m_CoreMass = p_CoreMass; }
    void ResetHeCoreMass()                                                          { m_HeCoreMass = std::nullopt; }
    void SetHeCoreMass(const double p_HeCoreMass)                                   { m_HeCoreMass = p_HeCoreMass; }
        
    void ResetDMdt()                                                                { m_dMdt = std::nullopt; }
    void SetDMdt(const double p_dMdt)                                               { m_dMdt = p_dMdt; }
        
    void ResetDominantMassLossType()                                                { m_DominantMassLossType = std::nullopt; }
    void SetDominantMassLossType(const MASS_LOSS_TYPE p_DominantMassLossType)       { m_DominantMassLossType = p_DominantMassLossType; }
        
    void ResetDt()                                                                  { m_dt = std::nullopt; }
    void SetDt(const double p_dt)                                                   { m_dt = p_dt; }
        
    void ResetError()                                                               { m_Error = std::nullopt; }
    void SetError(const ERROR p_Error)                                              { m_Error = p_Error; }
    void ResetEvolutionStatus()                                                     { m_EvolutionStatus = std::nullopt; }
    void SetEvolutionStatus(const EVOLUTION_STATUS p_EvolutionStatus)               { m_EvolutionStatus = p_EvolutionStatus; }
        
    void ResetGBparams()                                                            { for (SizeT idx = 0; idx < utils::to_underlying(HURLEY_GBP::COUNT); idx++) m_HurleyGBparams[idx] = std::nullopt; }
    void SetGBparams(const OptDblVectorT& p_GBparams)                               { m_HurleyGBparams = p_GBparams; }
    void ResetGBparams(const HURLEY_GBP p_GBParam)                                  { m_HurleyGBparams[utils::to_underlying(p_GBParam)] = std::nullopt; }
    void SetGBparams(const HURLEY_GBP p_GBParam, const double p_Value)              { m_HurleyGBparams[utils::to_underlying(p_GBParam)] = p_Value; }

    void ResetHAbundanceCore()                                                      { m_HAbundanceCore = std::nullopt; }
    void SetHAbundanceCore(const double p_HAbundance)                               { m_HAbundanceCore = p_HAbundance; }
    void ResetHAbundanceSurface()                                                   { m_HAbundanceSurface = std::nullopt; }
    void SetHAbundanceSurface(const double p_HAbundance)                            { m_HAbundanceSurface = p_HAbundance; }
    void ResetHeAbundanceCore()                                                     { m_HeAbundanceCore = std::nullopt; }
    void SetHeAbundanceCore(const double p_HeAbundance)                             { m_HeAbundanceCore = p_HeAbundance; }
    void ResetHeAbundanceSurface()                                                  { m_HeAbundanceSurface = std::nullopt; }
    void SetHeAbundanceSurface(const double p_HeAbundance)                          { m_HeAbundanceSurface = p_HeAbundance; }

    void ResetLuminosity()                                                          { m_Luminosity = std::nullopt; }
    void SetLuminosity(const double p_Luminosity)                                   { m_Luminosity = p_Luminosity; }
    void ResetLuminosity0()                                                         { ResetLuminosityZAMSEffective(); }
    void ResetLuminosityZAMSEffective()                                             { m_LuminosityZAMSEffective = std::nullopt; }
    void SetLuminosity0(const double p_Luminosity)                                  { SetLuminosityZAMSEffective(const double p_Luminosity); }
    void SetLuminosityZAMSEffective(const double p_Luminosity)                      { m_LuminosityZAMSEffective = p_Luminosity; }
        
    void ResetMass()                                                                { m_Mass = std::nullopt; }
    void SetMass(const double p_Mass)                                               { m_Mass = p_Mass; }
    void ResetMass0()                                                               { ResetMassEffectiveInitial(); }
    void ResetMassEffectiveInitial()                                                { m_MassEffectiveInitial = std::nullopt; }
    void SetMass0(const double p_Mass)                                              { SetMassEffectiveInitial(const double p_Mass); }
    void SetMassEffectiveInitial(const double p_Mass)                               { m_MassEffectiveInitial = p_Mass; }

    void ResetOmegaZAMS()                                                           { m_OmegaZAMS = std::nullopt; }
    void SetOmegaZAMS(const double p_OmegaZAMS)                                     { m_OmegaZAMS = p_OmegaZAMS; }

    void ResetPulsarDetails()                                                       { m_PulsarDetails = std::nullopt; }
    void SetPulsarDetails(const PulsarDetailsT& p_PulsarDetails)                    { m_PulsarDetails = p_PulsarDetails; }
        
    void ResetRadius()                                                              { m_Radius = std::nullopt; }
    void SetRadius(const double p_Radius)                                           { m_Radius = p_Radius; }
    void ResetRadius0()                                                             { ResetRadiusZAMSEffective(); }
    void ResetRadiusZAMSEffective()                                                 { m_RadiusZAMSEffective = std::nullopt; }
    void SetRadius0(const double p_Radius)                                          { SetRadiusZAMSEffective(const double p_Radius); }
    void SetRadiusZAMSEffective(const double p_Radius)                              { m_RadiusZAMSEffective = p_Radius; }
        
    void ResetRandomSeed()                                                          { m_RandomSeed = std::nullopt; }
    void SetRandomSeed(const ULongT p_RandomSeed)                                   { m_RandomSeed = p_RandomSeed; }
        
    void ResetStellarType()                                                         { m_StellarType = std::nullopt; }
    void SetStellarType(const STELLAR_TYPE p_StellarType)                           { m_StellarType = p_StellarType; }
        
    void ResetCurrentKickParameters()                                               { m_SupernovaDetails.currentKickParams = std::nullopt; }
    void SetCurrentKickParameters(const StellarKickParmsT& p_KickParams)            { m_SupernovaDetails.currentKickParams = p_KickParams; }
        
    void ResetSuppliedKickParameters()                                              { m_SupernovaDetails.suppliedKickParams = std::nullopt; }
    void SetSuppliedKickParameters(const StellarKickParmsT& p_KickParams)           { m_SupernovaDetails.suppliedKickParams = p_KickParams; }

    void ResetSupernovaDetails()                                                    { m_SupernovaDetails = std::nullopt; }
    void SetSupernovaDetails(const StellarSupernovaDetailsT& p_SupernovaDetails)    { m_SupernovaDetails = p_SupernovaDetails; }
        
    void ResetTau()                                                                 { m_Tau = std::nullopt; }
    void SetTau(const double p_Tau)                                                 { m_Tau = p_Tau; }
        
    void ResetTemperature()                                                         { m_Temperature = std::nullopt; }
    void SetTemperature(const double p_Temperature)                                 { m_Temperature = p_Temperature; }
        
    void ResetTime()                                                                { m_Time = std::nullopt; }
    void SetTime(const double p_Time)                                               { m_Time = p_Time; }

    void ResetTimescales()                                                          { for (SizeT idx = 0; idx < utils::to_underlying(TIMESCALE::COUNT); idx++) m_HurleyTimescales[idx] = std::nullopt; }
    void SetTimescales(const OptDblVectorT& p_Timescales)                           { m_HurleyTimescales = p_Timescales; }
    void ResetTimescales(const TIMESCALE p_Timescale)                               { m_HurleyTimescales[utils::to_underlying(p_Timescale)] = std::nullopt; }
    void SetTimescales(const TIMESCALE p_Timescale, const double p_Value)           { m_HurleyTimescales[utils::to_underlying(p_Timescale)] = p_Value; }


    // member functions

    /*
     * UnPackOptDblVector
     *
     * @brief
     * Unpacks an OptDblVectorT to a DblVectorT.
     * Converts the vector of std::optional<<double>> values to a vector of double values.
     * 
     * If any of the values in the std::optional<<double>> do not have a value assigned,
     * this function will fail and terminate the program. 
     * 
     * 
     * DblVectorT UnPackOptDblVector(const OptDblVectorT p_Vec) const
     * 
     * @param       p_Vec                           The vector to be unpacked
     * @return                                      Unpacked vector
     */
    GNU_CONST inline DblVectorT UnPackOptDblVector(const OptDblVectorT p_Vec) const {
        DblVectorT unpacked;
        for (SizeT idx = 0; idx < p_Vec.size(); idx++) unpacked.push_back(GetValue(p_Vec[idx].value()));
        return unpacked;
    }


    /*
     * ResetState
     *
     * @brief
     * Resets all state member variables to the uninitialised state
     * 
     * void ResetState()
     */
    void ResetState() {

        m_CHE                        = ResetCHE();

        m_CNOProcessedCoreMass       = ResetCNOProcessedCoreMass();
        m_COCoreMass                 = ResetCOCoreMas();
        m_CoreMass                   = ResetCoreMass();
        m_HeCoreMass                 = ResetHeCoreMass();
        
        m_dMdt                       = ResetDMdt();
        
        m_DominantMassLossType       = ResetDominantMassLossType();
        
        m_dt                         = ResetDt();
        
        m_Error                      = ResetError();
        m_EvolutionStatus            = ResetEvolutionStatus();
        
        m_HurleyGBparams             = ResetHurleyGBparams();

        m_HAbundanceCore             = ResetHAbundanceCore();
        m_HAbundanceSurface          = ResetHAbundanceSurface();
        m_HeAbundanceCore            = ResetHeAbundanceCore();
        m_HeAbundanceSurface         = ResetHeAbundanceSurface();
        
        m_Luminosity                 = ResetLuminosity();
        m_LuminosityZAMSEffective    = ResetLuminosityZAMSEffective();
        
        m_Mass                       = ResetMass();
        m_MassEffectiveInitial       = ResetMassEffectiveInitial();

        m_OmegaZAMS                  = ResetOmegaZAMS();       

        m_PulsarDetails              = ResetPulsarDetails();

        m_Radius                     = ResetRadius();
        m_RadiusZAMSInitial          = ResetRadiusZAMSInitial();
        
        m_RandomSeed                 = ResetRandomSeed();
        
        m_StellarType                = ResetStellarType();

        m_SupernovaDetails           = ResetSupernovaDetails();

        m_Tau                        = ResetTau();
        
        m_Temperature                = ResetTemperature();
        
        m_Time                       = ResetTime();

        m_HurleyTimescales           = ResetHurleyTimescales();
    }
};


// BinaryStarState class - state for individual binary stars
class BinaryStarState {       
        
private:
        
    // Vector of historical MT events as the donor.
    // The vector will always exist - empty vector indicates no events
    std::vector<MTEventT> m_MTdonorHistory;



    // member variables - alphabetically

    StellarCEDetailsT       m_CEDetails;                            // Common envelope details

    bool                    m_FirstMassTransferEpisode;             // Activated for the initial Mass Transfer Episode

    struct FLAGS {                                                  // Miscellaneous flags
        bool recycledNS;                                            // Indicate whether the accretor was a recycled neutron star
    }                       m_Flags;

    double                  m_MassLossDiff;
    double                  m_MassTransferDiff;

    double                  m_OmegaTidesIndividualDiff;

    double                  m_OrbitalEnergyPostSN;
    double                  m_OrbitalEnergyPreSN;

    StellarRLOFDetailsT     m_RLOFDetails;


    // the companion - set by calling SetCompanion()
    BinaryConstituentStar  *m_Companion;






public:

    // BinaryStarState class
    // constructor
    BinaryStarState() {}
        

    // setters

    // Add a mass trasfer even to m_MTdonorHistory.
    // Event only added if the first event as donor stellar type.
    void AddMtdonorEvent(const STELLAR_TYPE p_DonorST, const STELLAR_TYPE p_AccretorST, const MT_CASE p_MTcase) { 

        if (m_MTdonorHistory.empty()) {                                                     // first MT donor event?
            m_MTdonorHistory.push_back({ p_DonorST, p_AccretorST, p_MTcase });              // yes - set first event
        }                                                                                   // no - have existing events
        else if (!utils::IsOneOf(p_DonorST, { m_MTdonorHistory.back().donorST })) {         // first MT donor event as donor stellar type?
            m_MTdonorHistory.push_back({ p_DonorST, p_AccretorST, p_MTcase });              // yes - add new event
        }
    }


    // getters

    std::vector<MTEventT> MTdonorHistory() const { return m_MTdonorHistory; }
    
    MTEventT FirstMTdonorEvent(get::unwrappedT = {}) const { if (m_MTdonorHistory.empty()) Fail(ERROR::NO_VALUE); return m_MTdonorHistory.front(); }
    std::optional<MTEventT> FirstMTdonorEvent(get::asOptionalT) const { if (m_MTdonorHistory.empty()) return std::nullopt; return std::optional<m_MTdonorHistory.front()>; }
    
    MTEventT LatestMTdonorEvent() const { if (m_MTdonorHistory.empty()) Fail(ERROR::NO_VALUE); return m_MTdonorHistory.back(); }
    std::optional<MTEventT> LatestMTdonorEvent() const { if (m_MTdonorHistory.empty()) return std::nullopt; return std::optional<m_MTdonorHistory.back()>; }


};


// StellarBinaryState class - state for stellar binary systems
class StellarBinaryState {       
        
private:
        



};


template <typename T>
class StateHistory {

private:

    std::list<std::unique_ptr<T>> m_History;                                                // State history - doubly-linked list
    SizeT                         m_Capacity;                                               // History capacity (maximum number of elements)

public:

    // StateHistory class
    // constructor
    StateHistory(const T& p_StartState, const SizeT p_Capacity = DEFAULT_STATE_HISTORY_STACK_SIZE) {
        m_Capacity = p_Capacity == 0 ? 0 : std::max(3, static_cast<int>(p_Capacity));       // 0, or minimum of 3 - 0 means unbounded
        m_History.clear();
        m_History.push_back(std::make_unique<T>(p_StartState));                             // Current state 
    }

    StateHistory() : StateHistory(T{}) {}


    // member functions

    /*
     * Push
     *
     * @brief
     * Adds an element at the top of the stack (pushes the element onto the stack).
     * This results in existing elements on the stack being moved down (in rank, from
     * the top-of-stack (TOS)).  If the capacity of the stack is not bounded, all prior
     * existing elements are just moved down (i.e. the new element becomes TOS, the
     * element that was previously TOS is now at TOS-1, etc.).  If the capacity of the
     * stack is bounded, and the stack is already full, an element needs to be removed
     * from the stack to make space for the element being pushed.  Since the element at
     * the bottom of the stack (the initial state) is fixed and cannot be removed, the
     * element immediately above BOS is removed so that all elements above it can be
     * pushed down the stack to make room for the new TOS element.
     * 
     * If an element is removed from the stack to make room for the element being pushed,
     * a copy of the removed element is returned to the caller.  If no element is
     * removed, the returned object has no value assigned (see std::optional/std::nullopt).
     * The caller can ignore the returned object with (void).
     * 
     * 
     * template <typename T>
     * std::optional<T> Push(const T p_State)
     * 
     * @param       p_State                         Element to add to the stack (at TOS)
     * @return                                      std::optional<T> object containing removed element (if it exits)
     */
    std::optional<T> Push(const T p_State) {
        std::optional<T> retVal = std::nullopt;                                             // default return value has no value assigned

        if (m_Capacity > 0 && m_History.size() >= m_Capacity) {                             // need to make room?
                                                                                            // yes
            typename std::list<std::unique_ptr<T>>::iterator it = ++(m_History.begin());    // element to be erased
            retVal = **it;                                                                  // copy element to be erased
            m_History.erase(it);                                                            // erase element
        }

        m_History.push_back(std::make_unique<T>(p_State));                                  // add new element to end (TOS)

        return retVal;
    }
    std::optional<T> AddState(const T p_State) { return Push(p_State); }                    // AddState() is an alias for Push()


    /*
     * Pop
     *
     * @brief
     * Removes the element at top-of-stack (TOS) (pops the element off the stack), if
     * appropriate (the initial state, the element at the bottom of the stack, will never
     * be removed by this function).  This results in the remaining existing elements on
     * the stack being moved up (in rank, towards the top-of-stack (TOS)).
     * 
     * If an element is removed from the stack, a copy of the removed element is returned
     * to the caller.  If no element is removed, the returned object has no value assigned
     * (see std::optional/std::nullopt). The caller can ignore the returned object with (void).
     * 
     * 
     * template <typename T>
     * std::optional<T> Pop()
     * 
     * @return                                      std::optional<T> object containing removed element (if it exits)
     */
    std::optional<T> Pop() {
        std::optional<T> retVal = std::nullopt;                                             // default return value has no value assigned

        if (m_History.size() > 1) {                                                         // prevent removal of state zero
            retVal = *(m_History.back());                                                   // copy last element (TOS)
            m_History.pop_back();                                                           // remove last element (TOS)
        }

        return retVal;
    }
    std::optional<T> RemoveState() { return Pop(); }                                        // RemoveState() is an alias for Pop()

    
    /*
     * Peek
     *
     * @brief
     * Returns a copy of the stack element specified by the p_Count parameter.  The p_Count
     * parameter specifies the number of elements down the stack (below top-of-stack (TOS))
     * to look to find the element to return, with p_Count = 0 indicating the element at TOS.
     * When p_Count is not specified by the caller, p_Count defaults to 0, and the element at
     * TOS is returned.
     * 
     * If the stack element requested does not exist, the returned object will have no value
     * assigned (see std::optional/std::nullopt)
     * 
     * 
     * template <typename T>
     * std::optional<T> Peek(const SizeT p_Count = 0) const
     * 
     * @param       p_Count                         Element to retrieve (the number of elements below top-of-stack)
     * @return                                      std::optional<T> object containing element (if it exits)
     */
    std::optional<T> Peek(const SizeT p_Count = 0) const {
        std::optional<T> retVal = std::nullopt;                                             // default return value has no value assigned

        if (p_Count < m_History.size()) {                                                   // requested element exists?
            if (p_Count == 0) retVal = *(m_History.back());                                 // yes - TOS?
            else {                                                                          // not TOS - iterate to requested element
                SizeT count = 1;                                                            // start at TOS - 1
                for (typename std::list<std::unique_ptr<T>>::const_reverse_iterator it = ++(m_History.rbegin()); it != m_History.rend(); ++it) {
                    if (count++ == p_Count) {                                               // have requested element?
                        retVal = **it;                                                      // yes
                        break;                                                              // have requested element - stop iterating
                    }
                }
            }
        }

        return retVal;
    }
    std::optional<T> LookBack(const SizeT p_Count = 0) const { return Peek(p_Count); }      // LookBack() is an alias for Peek()


    // bookkeeping
    SizeT Size() const  { return m_History.size(); }
    SizeT Count() const { return Size(); }


    // iterators - read only - forward and reverse
    typename std::list<std::unique_ptr<T>>::const_iterator Begin() const          { return m_History.begin(); }
    typename std::list<std::unique_ptr<T>>::const_iterator End() const            { return m_History.end(); }
    typename std::list<std::unique_ptr<T>>::const_reverse_iterator rBegin() const { return m_History.rbegin(); }
    typename std::list<std::unique_ptr<T>>::const_reverse_iterator rEnd() const   { return m_History.rend(); }


    // sentinel getters
    T TOS() const { return *(m_History.back()); }                                           // top of stack (TOS) - will always have a value assigned (see constructor)
    T BOS() const { return *(m_History.front()); }                                          // bottom of stack (BOS) - will always have a value assigned (see constructor)


    // functions to determine if defined states exist
    // StartState and CurrentState are guaranteed to exist (see constructor)
    bool HaveZAMSState() const     { return utils::IsOneOf(StartState().StellarType(), MAIN_SEQUENCE); }
    bool HavePreviousState() const { return m_History.size() > 1; }


    // defined state getters
    // 
    //    - StartState and CurrentState are guaranteed to exist (see constructor)
    // 
    //    - ZAMSState may not exist (if star being simulated did not start on the MS),and the
    //      getter will fail if called when the state does not exist.  If there is doubt,
    //      callers should use HaveZAMSState() to check if ZAMSState exists.
    //
    //     - PreviousState may not exist (should only happen on first timestep), and the getter
    //       will fail if called when the state does not exist.  If there is doubt, callers
    //       should use HavePreviousState() to check if PreviousState exists.
    std::optional<T> StartState() const    { return BOS(); }
    std::optional<T> CurrentState() const  { return TOS(); }

    std::optional<T> ZAMSState() const     { if (!HaveZAMSState()) Fail(ERROR::NO_STATE) else return StartState(); }
    std::optional<T> PreviousState() const { if (!HavePreviousState()) Fail(ERROR::NO_STATE) else return Peek(1); }
};
