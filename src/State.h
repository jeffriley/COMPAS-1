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





// StarState class - state for individual stars
class StarState {       
        
private:
        
    std::optional<double>            m_AngularFrequency;
        
    std::optional<double>            m_CoreMass;
    std::optional<double>            m_CoreMassCO;
    std::optional<double>            m_CoreMassEffective;
    std::optional<double>            m_CoreMassHe;
        
    std::optional<double>            m_dMdt;
        
    std::optional<MASS_LOSS_TYPE>    m_DominantMassLossType;
        
    std::optional<double>            m_dt;
        
    std::optional<ERROR>             m_Error;
    std::optional<EVOLUTION_STATUS>  m_EvolutionStatus;
        
    std::optional<double>            m_HAbundance;
    std::optional<double>            m_HeAbundance;
        
    std::optional<KickParametersT>   m_KickParameters;
        
    std::optional<double>            m_Luminosity;
    std::optional<double>            m_LuminosityEffective;
        
    std::optional<double>            m_Mass;
    std::optional<double>            m_MassEffective;
        
    std::optional<double>            m_Metallicity;
        
    std::optional<double>            m_Radius;
    std::optional<double>            m_RadiusEffective;
        
    std::optional<unsigned long int> m_RandomSeed;
        
    std::optional<STELLAR_TYPE>      m_StellarType;
        
    std::optional<double>            m_Tau;
        
    std::optional<double>            m_Temperature;
        
    std::optional<double>            m_Time;       
        
        
public:
            
    StarState() {}
        
    // setters
    void SetAngularFrequency(const double p_AngularFrequency)                 { m_AngularFrequency = p_AngularFrequency; }
        
    void SetCoreMass(const double p_CoreMass)                                 { m_CoreMass = p_CoreMass; }
    void SetCoreMassCO(const double p_CoreMassCO)                             { m_CoreMassCO = p_CoreMassCO; }
    void SetCoreMassEffective(const double p_CoreMass)                        { m_CoreMassEffective = p_CoreMass; }
    void SetCoreMassHe(const double p_CoreMassHe)                             { m_CoreMassHe = p_CoreMassHe; }
        
    void SetDMdt(const double p_dMdt)                                         { m_dMdt = p_dMdt; }
        
    void SetDominantMassLossType(const MASS_LOSS_TYPE p_DominantMassLossType) { m_DominantMassLossType = p_DominantMassLossType; }
        
    void SetDt(const double p_dt)                                             { m_dt = p_dt; }
        
    void SetError(const ERROR p_Error)                                        { m_Error = p_Error; }
    void SetEvolutionStatus(const EVOLUTION_STATUS p_EvolutionStatus)         { m_EvolutionStatus = p_EvolutionStatus; }
        
    void SetHAbundance(const double p_HAbundance)                             { m_HAbundance = p_HAbundance; }
    void SetHeAbundance(const double p_HeAbundance)                           { m_HeAbundance = p_HeAbundance; }
        
    void SetKickParameters(const KickParametersT p_KickParameters)            { m_KickParameters = p_KickParameters; }
        
    void SetLuminosity(const double p_Luminosity)                             { m_Luminosity = p_Luminosity; }
    void SetLuminosityEffective(const double p_LuminosityEffective)           { m_LuminosityEffective = p_LuminosityEffective; }
        
    void SetMass(const double p_Mass)                                         { m_Mass = p_Mass; }
    void SetMassEffective(const double p_MassEffective)                       { m_MassEffective = p_MassEffective; }
        
    void SetMetallicity(const double p_Metallicity)                           { m_Metallicity = p_Metallicity; }
        
    void SetRadius(const double p_Radius)                                     { m_Radius = p_Radius; }
    void SetRadiusEffective(const double p_RadiusEffective)                   { m_RadiusEffective = p_RadiusEffective; }
        
    void SetRandomSeed(const unsigned long int p_RandomSeed)                  { m_RandomSeed = p_RandomSeed; }
        
    void SetStellarType(const STELLAR_TYPE p_StellarType)                     { m_StellarType = p_StellarType; }
        
    void SetTau(const double p_Tau)                                           { m_Tau = p_Tau; }
        
    void SetTemperature(const double p_Temperature)                           { m_Temperature = p_Temperature; }
        
    void SetTime(const double p_Time)                                         { m_Time = p_Time; }
               
        
    // getters

    #define FAIL         { std::cerr << "\nState attribute has no value: program terminated\n"; utils::ShowStackTrace(); std::exit(1); }
    #define RET_VALUE(x) { if (x.has_value()) { return x.value(); } else FAIL }

    double AngularFrequency() const             { RET_VALUE(m_AngularFrequency); }
        
    double CoreMass() const                     { RET_VALUE(m_CoreMass); }
    double CoreMassCO() const                   { RET_VALUE(m_CoreMassCO); }
    double CoreMassEffective() const            { RET_VALUE(m_CoreMassEffective); }
    double CoreMassHe() const                   { RET_VALUE(m_CoreMassHe); }
        
    double dMdt() const                         { RET_VALUE(m_dMdt); }
        
    MASS_LOSS_TYPE DominantMassLossType() const { RET_VALUE(m_DominantMassLossType); }
        
    double dt() const                           { RET_VALUE(m_dt); }
        
    ERROR Error() const                         { RET_VALUE(m_Error); }
    EVOLUTION_STATUS EvolutionStatus() const    { RET_VALUE(m_EvolutionStatus); }
        
    double HAbundance() const                   { RET_VALUE(m_HAbundance); }
    double HeAbundance() const                  { RET_VALUE(m_HeAbundance); }
        
    KickParametersT KickParameters() const      { RET_VALUE(m_KickParameters); }
        
    double Luminosity() const                   { RET_VALUE(m_Luminosity); }
    double LuminosityEffective() const          { RET_VALUE(m_LuminosityEffective); }
        
    double Mass() const                         { RET_VALUE(m_Mass); }
    double MassEffective() const                { RET_VALUE(m_MassEffective); }
        
    double Metallicity() const                  { RET_VALUE(m_Metallicity); }
        
    double Radius() const                       { RET_VALUE(m_Radius); }
    double RadiusEffective() const              { RET_VALUE(m_RadiusEffective); }
        
    unsigned long int RandomSeed() const        { RET_VALUE(m_RandomSeed); }
        
    STELLAR_TYPE StellarType() const            { RET_VALUE(m_StellarType); }
        
    double Tau() const                          { RET_VALUE(m_Tau); }
        
    double Temperature() const                  { RET_VALUE(m_Temperature); }
        
    double Time() const                         { RET_VALUE(m_Time); }   

    #undef RET_VALUE
    #undef FAIL


    // has value
    bool AngularFrequency_HasValue() const     { return m_AngularFrequency.has_value(); }
        
    bool CoreMass_HasValue() const             { return m_CoreMass.has_value(); }
    bool CoreMassCO_HasValue() const           { return m_CoreMassCO.has_value(); }
    bool CoreMassEffective_HasValue() const    { return m_CoreMassEffective.has_value(); }
    bool CoreMassHe_HasValue() const           { return m_CoreMassHe.has_value(); }
        
    bool dMdt_HasValue() const                 { return m_dMdt.has_value(); }
        
    bool DominantMassLossType_HasValue() const { return m_DominantMassLossType.has_value(); }
        
    bool dt_HasValue() const                   { return m_dt.has_value(); }
        
    bool Error_HasValue() const                { return m_Error.has_value(); }
    bool EvolutionStatus_HasValue() const      { return m_EvolutionStatus.has_value(); }
        
    bool HAbundance_HasValue() const           { return m_HAbundance.has_value(); }
    bool HeAbundance_HasValue() const          { return m_HeAbundance.has_value(); }
        
    bool KickParameters_HasValue() const       { return m_KickParameters.has_value(); }
        
    bool Luminosity_HasValue() const           { return m_Luminosity.has_value(); }
    bool LuminosityEffective_HasValue() const  { return m_LuminosityEffective.has_value(); }
        
    bool Mass_HasValue() const                 { return m_Mass.has_value(); }
    bool MassEffective_HasValue() const        { return m_MassEffective.has_value(); }
        
    bool Metallicity_HasValue() const          { return m_Metallicity.has_value(); }
        
    bool Radius_HasValue() const               { return m_Radius.has_value(); }
    bool RadiusEffective_HasValue() const      { return m_RadiusEffective.has_value(); }
        
    bool RandomSeed_HasValue() const           { return m_RandomSeed.has_value(); }
        
    bool StellarType_HasValue() const          { return m_StellarType.has_value(); }
        
    bool Tau_HasValue() const                  { return m_Tau.has_value(); }
        
    bool Temperature_HasValue() const          { return m_Temperature.has_value(); }
        
    bool Time_HasValue() const                 { return m_Time.has_value(); }
};



template <typename T>
class StateHistory {

private:

    std::list<std::unique_ptr<T>> m_History;                                                // state history - doubly-linked list
    size_t                        m_Capacity;                                               // history capacity (maximum number of elements)

public:



    // StateHistrory class
    // constructor
    StateHistory(const T p_FirstState, const size_t p_Capacity = DEFAULT_STATE_HISTORY_STACK_SIZE) {
        m_Capacity = p_Capacity == 0 ? 0 : std::max(3, static_cast<int>(p_Capacity));
        m_History.clear();
        auto state = std::make_unique<T>(p_FirstState);                                     // allocate entry for list
        m_History.push_back(std::make_unique<T>(p_FirstState)); 
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
    // FirstState and CurrentState are guaranteed to exist (see constructor)
    bool HaveZAMSState() const     { return utils::IsOneOf(FirstState().StellarType(), MAIN_SEQUENCE); }
    bool HavePreviousState() const { return m_History.size() > 1; }


    // defined state getters
    // 
    //    - FirstState and CurrentState are guaranteed to exist (see constructor)
    // 
    //    - ZAMSState may not exist (if star being simulated did not start on the MS),and the
    //      getter will fail if called when the state does not exist.  If there is doubt,
    //      callers should use HaveZAMSState() to check if ZAMSState exists.
    //
    //     - PreviousState may not exist (should only happen on first timestep), and the getter
    //       will fail if called when the state does not exist.  If there is doubt, callers
    //       should use HavePreviousState() to check if PreviousState exists.
    #define FAIL { std::cerr << "\nState does not exist: program terminated\n"; utils::ShowStackTrace(); std::exit(1); }

    T FirstState() const    { return BOS(); }
    T CurrentState() const  { return TOS(); }

    T ZAMSState() const     { if (!HaveZAMSState()) FAIL else return FirstState(); }
    T PreviousState() const { if (!HavePreviousState()) FAIL else return Peek(1); }

    #undef FAIL
};
