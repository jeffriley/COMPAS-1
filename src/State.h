#include <list>
#include <optional>

#include "constants.h"

//#include "utils.h"

//constexpr std::size_t DEFAULT_STATE_HISTORY_STACK_SIZE = 0; 

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





// StarState class - state for isolated single stars
class StarState {       

    #define FAIL         { std::cerr << "\nState attribute has no value: program terminated\n"; utils::ShowStackTrace(); std::exit(1); }
    #define RET_VALUE(x) { if (x.has_value()) { return x.value(); } else FAIL }
        
private:
        
    std::optional<double>            m_AngularFrequency;
        
    std::optional<double>            m_COCoreMass;
    std::optional<double>            m_CoreMass;
    std::optional<double>            m_CoreMassEffective;
        
    std::optional<double>            m_dMdt;
        
    std::optional<MASS_LOSS_TYPE>    m_DominantMassLossType;
        
    std::optional<double>            m_dt;
        
    std::optional<ERROR>             m_Error;
    std::optional<EVOLUTION_STATUS>  m_EvolutionStatus;
        
    std::optional<double>            m_HAbundanceCore;
    std::optional<double>            m_HAbundanceSurface;
    std::optional<double>            m_HeAbundanceCore;
    std::optional<double>            m_HeAbundanceSurface;

    std::optional<double>            m_HeCoreMass;
        
    std::optional<KickParametersT>   m_KickParameters;

    std::optional<double>            m_Luminosity;
    std::optional<double>            m_LuminosityEffectiveInitial;
        
    std::optional<double>            m_Mass;
    std::optional<double>            m_MassEffectiveInitial;
        
    std::optional<double>            m_Radius;
    std::optional<double>            m_RadiusEffectiveInitial;
        
    std::optional<unsigned long int> m_RandomSeed;
        
    std::optional<STELLAR_TYPE>      m_StellarType;
        
    std::optional<double>            m_Tau;
        
    std::optional<double>            m_Temperature;
        
    std::optional<double>            m_Time;       

    // timescales and GB parameters vectors (see Hurley et al. 2000)
    // Each vector is constructed here with the correct number of elements, and each
    // element is wrapped in std::optional<>, and initialised to std::nullopt.
    // The vector is guaranteed to exist, but each element is optional.

    OPT_DBL_VECTOR m_HurleyTimescales;
    OPT_DBL_VECTOR m_HurleyGBparams;

        
public:

    // StarState class
    // constructor
    StarState() {

        // create GB parameters and timescales vectors
        // construct vectors of required sizes and initialise each element to std::nullopt
        for (size_t idx = 0; idx < static_cast<int>(HURLEY_GBP:::COUNT); idx++) m_HurleyGBparams.push_back(std::nullopt);
        for (size_t idx = 0; idx < static_cast<int>(TIMESCALE::COUNT); idx++) m_HurleyTimescales.push_back(std::nullopt);
    }
        

    // setters

    void SetAngularFrequency(const double p_AngularFrequency)                 { m_AngularFrequency = p_AngularFrequency; }
        
    void SetCOCoreMass(const double p_COCoreMass)                             { m_COCoreMass = p_COCoreMass; }
    void SetCoreMass(const double p_CoreMass)                                 { m_CoreMass = p_CoreMass; }
    void SetCoreMassEffective(const double p_CoreMass)                        { m_CoreMassEffective = p_CoreMass; }
    void SetHeCoreMass(const double p_HeCoreMass)                             { m_HeCoreMass = p_HeCoreMass; }
        
    void SetDMdt(const double p_dMdt)                                         { m_dMdt = p_dMdt; }
        
    void SetDominantMassLossType(const MASS_LOSS_TYPE p_DominantMassLossType) { m_DominantMassLossType = p_DominantMassLossType; }
        
    void SetDt(const double p_dt)                                             { m_dt = p_dt; }
        
    void SetError(const ERROR p_Error)                                        { m_Error = p_Error; }
    void SetEvolutionStatus(const EVOLUTION_STATUS p_EvolutionStatus)         { m_EvolutionStatus = p_EvolutionStatus; }
        
    void SetGBparams(const OPT_DBL_VECTOR p_GBparams)                         { m_HurleyGBparams = p_GBparams; }
    void SetGBparams(const GBP p_GBParam, const double p_Value)               { m_HurleyGBparams[static_cast<int>(p_GBParam)] = p_Value; }

    void SetHAbundanceCore(const double p_HAbundance)                         { m_HAbundanceCore = p_HAbundance; }
    void SetHAbundanceSurface(const double p_HAbundance)                      { m_HAbundanceSurface = p_HAbundance; }
    void SetHeAbundanceCore(const double p_HeAbundance)                       { m_HeAbundanceCore = p_HeAbundance; }
    void SetHeAbundanceSurface(const double p_HeAbundance)                    { m_HeAbundanceSurface = p_HeAbundance; }
        
    void SetKickParameters(const KickParametersT p_KickParameters)            { m_KickParameters = p_KickParameters; }
        
    void SetLuminosity(const double p_Luminosity)                             { m_Luminosity = p_Luminosity; }
    void SetLuminosityEffectiveInitial(const double p_Luminosity)             { m_LuminosityEffectiveInitial = p_Luminosity; }
        
    void SetMass(const double p_Mass)                                         { m_Mass = p_Mass; }
    void SetMassEffectiveInitial(const double p_Mass)                         { m_MassEffectiveInitial = p_Mass; }
        
    void SetRadius(const double p_Radius)                                     { m_Radius = p_Radius; }
    void SetRadiusEffectiveInitial(const double p_Radius)                     { m_RadiusEffectiveInitial = p_Radius; }
        
    void SetRandomSeed(const unsigned long int p_RandomSeed)                  { m_RandomSeed = p_RandomSeed; }
        
    void SetStellarType(const STELLAR_TYPE p_StellarType)                     { m_StellarType = p_StellarType; }
        
    void SetTau(const double p_Tau)                                           { m_Tau = p_Tau; }
        
    void SetTemperature(const double p_Temperature)                           { m_Temperature = p_Temperature; }
        
    void SetTime(const double p_Time)                                         { m_Time = p_Time; }

    void SetTimescales(const OPT_DBL_VECTOR p_Timescales)                     { m_HurleyTimescales = p_Timescales; }
    void SetTimescales(const TIMESCALE p_Timescale, const double p_Value)     { m_HurleyTimescales[static_cast<int>(p_Timescale)] = p_Value; }
               
        
    // getters   <<<<<<<<<<<<<<<<<<<<<<< REVISIT <<<<<<<<<<<<<<<<<<< from current state??? Ah, no, to be used with state qualifier - need to document

    double AngularFrequency() const                      { RET_VALUE(m_AngularFrequency); }
        
    double COCoreMass() const                            { RET_VALUE(m_COCoreMass); }
    double CoreMass() const                              { RET_VALUE(m_CoreMass); }
    double CoreMassEffective() const                     { RET_VALUE(m_CoreMassEffective); }
        
    double dMdt() const                                  { RET_VALUE(m_dMdt); }
        
    MASS_LOSS_TYPE DominantMassLossType() const          { RET_VALUE(m_DominantMassLossType); }
        
    double dt() const                                    { RET_VALUE(m_dt); }
        
    ERROR Error() const                                  { RET_VALUE(m_Error); }
    EVOLUTION_STATUS EvolutionStatus() const             { RET_VALUE(m_EvolutionStatus); }

    DBL_VECTOR GBparams() const                          { return UnPackOptDblVector(m_HurleyGBparams); }
    double GBparams(const GBP p_GBParam) const           { RET_VALUE(m_HurleyGBparams[static_cast<int>(p_GBParam)]); }
        
    double HAbundanceCore() const                        { RET_VALUE(m_HAbundanceCore); }
    double HAbundanceSurface() const                     { RET_VALUE(m_HAbundanceSurface); }
    double HeAbundanceCore() const                       { RET_VALUE(m_HeAbundanceCore); }
    double HeAbundanceSurface() const                    { RET_VALUE(m_HeAbundanceSurface); }

    double HeCoreMass() const                            { RET_VALUE(m_HeCoreMass); }
        
    KickParametersT KickParameters() const               { RET_VALUE(m_KickParameters); }
        
    double Luminosity() const                            { RET_VALUE(m_Luminosity); }
    double LuminosityEffectiveInitial() const            { RET_VALUE(m_LuminosityEffectiveInitial); }
        
    double Mass() const                                  { RET_VALUE(m_Mass); }
    double MassEffectiveInitial() const                  { RET_VALUE(m_MassEffectiveInitial); }
        
    double Radius() const                                { RET_VALUE(m_Radius); }
    double RadiusEffectiveInitial() const                { RET_VALUE(m_RadiusEffectiveInitial); }
        
    unsigned long int RandomSeed() const                 { RET_VALUE(m_RandomSeed); }
        
    STELLAR_TYPE StellarType() const                     { RET_VALUE(m_StellarType); }
        
    double Tau() const                                   { RET_VALUE(m_Tau); }
        
    double Temperature() const                           { RET_VALUE(m_Temperature); }
        
    double Time() const                                  { RET_VALUE(m_Time); }   

    DBL_VECTOR Timescales() const                        { return UnPackOptDblVector(m_HurleyTimescales); }
    double Timescales(const TIMESCALE p_Timescale) const { RET_VALUE(m_HurleyTimescales[static_cast<int>(p_Timescale)]); }


    // has value

    bool AngularFrequency_HasValue() const                      { return m_AngularFrequency.has_value(); }
        
    bool COCoreMass_HasValue() const                            { return m_COCoreMass.has_value(); }
    bool CoreMass_HasValue() const                              { return m_CoreMass.has_value(); }
    bool CoreMassEffective_HasValue() const                     { return m_CoreMassEffective.has_value(); }
        
    bool dMdt_HasValue() const                                  { return m_dMdt.has_value(); }
        
    bool DominantMassLossType_HasValue() const                  { return m_DominantMassLossType.has_value(); }
        
    bool dt_HasValue() const                                    { return m_dt.has_value(); }
        
    bool Error_HasValue() const                                 { return m_Error.has_value(); }
    bool EvolutionStatus_HasValue() const                       { return m_EvolutionStatus.has_value(); }

    bool GBparams_HasValue() const                              { // returns true iff all elements have values, else false
                                                                    for (size_t idx = 0; idx < m_HurleyGBparams.size(); idx++) {
                                                                        if (!m_HurleyGBparams[idx].has_value()) return false;
                                                                    }
                                                                    return true;
                                                                }
    bool GBparams_HasValue(const GBP p_GBParam) const           { return m_HurleyGBparams[static_cast<int>(p_GBParam)].has_value(); }
        
    bool HAbundanceCore_HasValue() const                        { return m_HAbundanceCore.has_value(); }
    bool HAbundanceSurface_HasValue() const                     { return m_HAbundanceCore.has_value(); }
    bool HeAbundanceCore_HasValue() const                       { return m_HeAbundanceCore.has_value(); }
    bool HeAbundanceSurface_HasValue() const                    { return m_HeAbundanceCore.has_value(); }

    bool HeCoreMass_HasValue() const                            { return m_HeCoreMass.has_value(); }
        
    bool KickParameters_HasValue() const                        { return m_KickParameters.has_value(); }
        
    bool Luminosity_HasValue() const                            { return m_Luminosity.has_value(); }
    bool LuminosityEffectiveInitial_HasValue() const            { return m_LuminosityEffectiveInitial.has_value(); }
        
    bool Mass_HasValue() const                                  { return m_Mass.has_value(); }
    bool MassEffectiveInitial_HasValue() const                  { return m_MassEffectiveInitial.has_value(); }
        
    bool Radius_HasValue() const                                { return m_Radius.has_value(); }
    bool RadiusEffectiveInitial_HasValue() const                { return m_RadiusEffectiveInitial.has_value(); }
        
    bool RandomSeed_HasValue() const                            { return m_RandomSeed.has_value(); }
        
    bool StellarType_HasValue() const                           { return m_StellarType.has_value(); }
        
    bool Tau_HasValue() const                                   { return m_Tau.has_value(); }
        
    bool Temperature_HasValue() const                           { return m_Temperature.has_value(); }
        
    bool Time_HasValue() const                                  { return m_Time.has_value(); }

    bool Timescales_HasValue() const                            { // returns true iff all elements have values, else false
                                                                    for (size_t idx = 0; idx < m_HurleyTimescales.size(); idx++) {
                                                                        if (!m_HurleyTimescales[idx].has_value()) return false;
                                                                    }
                                                                    return true;
                                                                }
    bool Timescales_HasValue(const TIMESCALE p_Timescale) const { return m_HurleyTimescales[static_cast<int>(p_Timescale)].has_value(); }


    // member functions

    /*
     * UnPackOptDblVector
     *
     * @brief
     * Unpacks an OPT_DB_VECTOR to a DBL_VECTOR.
     * Converts the vector of std::optional<<double>> values to a vector of double values.
     * 
     * If any of the values in the std::optional<<double>> do not have a value assigned,
     * this function will fail and terminate the program. 
     * 
     * 
     * std::optional<T> Push(const T p_State)
     * 
     * @param       p_State                         Element to add to the stack (at TOS)
     * @return                                      std::optional<T> object containing removed element (if it exits)
     */
    DBL_VECTOR UnPackOptDblVector(const OPT_DBL_VECTOR p_Vec) {
        DBL_VECTOR unpacked;
        for (size_t idx = 0; idx < p_Vec.size(); idx++) unpacked.push_back(RET_VALUE(p_Vec[idx].value()));
        return unpacked;
    }

    #undef RET_VALUE
    #undef FAIL
};


// BinaryStarState class - state for individual binary stars
class BinaryStarState {       

    #define FAIL         { std::cerr << "\nState attribute has no value: program terminated\n"; utils::ShowStackTrace(); std::exit(1); }
    #define RET_VALUE(x) { if (x.has_value()) { return x.value(); } else FAIL }
        
private:
        
    std::vector<MTEventT> m_MTdonorHistory;     // historical list of MT events as the donor - empty vector indicates no events

public:

    // BinaryStarState class
    // constructor
    BinaryStarState() {}
        

    // setters

    void AddMtdonorEvent(const STELLAR_TYPE p_DonorST, const STELLAR_TYPE p_AccretorST, const MT_CASE p_MTcase) { 

        if (!m_MTdonorHistory.has_value() || m_MTdonorHistory.empty()) {                    // first MT donor event?
            m_MTdonorHistory = { p_DonorST, p_AccretorST, p_MTcase };                       // yes - set first event
        }                                                                                   // no - have existing events
        else if (!utils::IsOneOf(p_DonorST, { m_MTdonorHistory.back().m_StellarType })) {   // first MT donor event as donor stellar type?
            m_MTdonorHistory.push_back({ p_DonorST, p_AccretorST, p_MTcase });              // yes - add new event
        }
    }


    // getters

    std::vector<MTEventT> MTdonorHistory() const { return m_MTdonorHistory; }
    
    MTEventT LatestMTdonorEvent() const { return m_MTdonorHistory.empty() ? { STELLAR_TYPE::NONE, STELLAR_TYPE::NONE, MT_CASE::NONE } : m_MTdonorHistory.back(); }


    // has value

    bool MTdonorHistory_HasValue() const { return !m_MTdonorHistory.empty(); } // true iff vector is not empty



};


// StellarBinaryState class - state for stellar binary systems
class StellarBinaryState {       

    #define FAIL         { std::cerr << "\nState attribute has no value: program terminated\n"; utils::ShowStackTrace(); std::exit(1); }
    #define RET_VALUE(x) { if (x.has_value()) { return x.value(); } else FAIL }
        
private:
        
    std::optional<double>            m_AngularFrequency;



};


template <typename T>
class StateHistory {

private:

    std::list<std::unique_ptr<T>> m_History;                                                // state history - doubly-linked list
    size_t                        m_Capacity;                                               // history capacity (maximum number of elements)

public:

    // StateHistory class
    // constructor
    StateHistory(const T p_StartState, const size_t p_Capacity = DEFAULT_STATE_HISTORY_STACK_SIZE) {
        m_Capacity = p_Capacity == 0 ? 0 : std::max(3, static_cast<int>(p_Capacity));
        m_History.clear();
        auto state = std::make_unique<T>(p_StartState);                                     // allocate entry for list
        m_History.push_back(std::make_unique<T>(p_StartState)); 
    }


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
     * std::optional<T> Peek(const size_t p_Count = 0) const
     * 
     * @param       p_Count                         Element to retrieve (the number of elements below top-of-stack)
     * @return                                      std::optional<T> object containing element (if it exits)
     */
    std::optional<T> Peek(const size_t p_Count = 0) const {
        std::optional<T> retVal = std::nullopt;                                             // default return value has no value assigned

        if (p_Count < m_History.size()) {                                                   // requested element exists?
            if (p_Count == 0) retVal = *(m_History.back());                                 // yes - TOS?
            else {                                                                          // not TOS - iterate to requested element
                size_t count = 1;                                                           // start at TOS - 1
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
    std::optional<T> LookBack(const size_t p_Count = 0) const { return Peek(p_Count); }     // LookBack() is an alias for Peek()


    // bookkeeping
    size_t Size() const                                          { return m_History.size(); }
    size_t Count() const                                         { return Size(); }


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
    #define FAIL { std::cerr << "\nState does not exist: program terminated\n"; utils::ShowStackTrace(); std::exit(1); }

    T StartState() const    { return BOS(); }
    T CurrentState() const  { return TOS(); }

    T ZAMSState() const     { if (!HaveZAMSState()) FAIL else return StartState(); }
    T PreviousState() const { if (!HavePreviousState()) FAIL else return Peek(1); }

    #undef FAIL
};
