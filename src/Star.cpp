// star class is a container holding the underlying object of the required stellar type
// needs to be done this way so that underlying object can be deleted and reconstructed
// as the required class for the stellar type
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#include <algorithm>
#include <csignal>
#include <fenv.h>

#include "Star.h"


// Default constructor
//Star::Star() : m_Star(new BaseStar()) {
//
//    m_ObjectId          = globalObjectId++;                                                                             // set object id
//    m_ObjectPersistence = OBJECT_PERSISTENCE::PERMANENT;                                                                // set object persistence
//}


Star::Star(const STARTING_STELLAR_TYPE p_StartingStellarType,
           const unsigned long int     p_RandomSeed,
           const double                p_Metallicity, 
           const double                p_Mass,
           const KickParameters        p_KickParameters,
           const double                p_RotationalFrequency) {

    m_ObjectId          = globalObjectId++;                                                                             // set object id
    m_ObjectPersistence = OBJECT_PERSISTENCE::PERMANENT;                                                                // set object persistence

    double angularFrequency = (p_RotationalFrequency >= 0.0) ? SECONDS_IN_YEAR * _2_PI * p_RotationalFrequency : 0.0;   // use rotational frequency if supplied

    // We need to convert the starting stellar type passed in to an actual stellar type so we can construct
    // a star of that type.  For most stellar types the mapping from starting stellar type to actual stellar
    // type is 1:1, but for the MS starting stellar type we need to deterine if the actual stellar type will
    // be MS_LTE_07, MS_GT_07, or CH (Chemically Homeogeneous).
    //
    // If CHE is enabled (OPTIONS->CHEMode(), program option `--chemically-homogeneous-evolution-mode`), and
    // if the ZAMS angular frequency of the star is greater than the threshold for CHE to occur, the stellar
    // type of the star will be CH, otherwise MS_LTE_07 or MS_GT_07 depending upon the ZAMS mass of the star.

    STELLAR_TYPE startingStellarType = static_cast<STELLAR_TYPE>(static_cast<int>(p_StartingStellarType));              // fix up starting stellar type

    if (utils::IsOneOf(startingStellarType, MAIN_SEQUENCE)) {                                                           // starting stellar type MS (won't be CH here)?
        if (OPTIONS->CHEMode() != CHE_MODE::NONE) {                                                                     // yes - CHE enabled?
                                                                                                                        // yes
            if (p_RotationalFrequency < 0.0) {                                                                          // rotational frequency supplied?
                const double RZAMS = CalculateZAMSRadius_Tout_Static(p_Mass, GLOBALS->ToutRadiusCoefficients());        // no
                angularFrequency   = BaseStar::CalculateZAMSAngularFrequency_Hurley_Static(p_Mass, RZAMS);              // calculate ZAMS omega
            }
            
            if (angularFrequency >= BaseStar::CalculateCHEAngularFrequency_Static(p_Mass, p_Metallicity)) {             // rotating fast enough to be CH?
                startingStellarType = STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS;                                             // yes - set starting stellar type to CH
            }
        }

        if (startingStellarType != STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS) {                                              // CH?
            if (p_Mass <= 0.7) {                                                                                        // no - mass <= 0.7?
                startingStellarType = STELLAR_TYPE::MS_LTE_07;                                                          // yes, set starting stellar type to MS_LTE_07
            }
            else {                                                                                                      // no
                startingStellarType = STELLAR_TYPE::MS_GT_07;                                                           // set starting stellar type to MS_GT_07
            }
        }
    }

    // construct a BaseStar object
    m_Star = new BaseStar(startingStellarType, p_RandomSeed, p_Metallicity, p_Mass, p_KickParameters, angularFrequency);                                                                                            // create underlying BaseStar object
    
    // switch to specified starting stellar type
    (void)SwitchTo(startingStellarType, true);

    // thresholds flags for system snapshot output file
    if (OPTIONS->SystemSnapshotAgeThresholds().size()  > 0) m_SystemSnapshotAgeFlags.assign(OPTIONS->SystemSnapshotAgeThresholds().size(), -1.0);
    if (OPTIONS->SystemSnapshotTimeThresholds().size() > 0) m_SystemSnapshotTimeFlags.assign(OPTIONS->SystemSnapshotTimeThresholds().size(), false);
}


// Copy constructor - deep copy so dynamic variables are also copied
Star::Star(const Star& p_Star) {

    m_ObjectId          = globalObjectId++;                                                                                             // set object id
    m_ObjectPersistence = p_Star.ObjectPersistence();                                                                                   // set object persistence

    m_Star = p_Star.m_Star ? static_cast<BaseStar*>(p_Star.m_Star->Clone(OBJECT_PERSISTENCE::PERMANENT, false)) : nullptr;              // copy underlying BaseStar object
}


/*
 * Switch to required star type
 *
 * Instantiates new object of required class, deletes existing pointer to star object and
 * replaces it with pointer to newly instantiated object
 *
 *
 * STELLAR_TYPE SwitchTo(const STELLAR_TYPE p_StellarType, bool p_SetInitialState)
 *
 * @param   [IN]    p_StellarType               StellarType to switch to
 * @param   [IN]    p_SetInitialType            Indicates whether the initial stellar type of the star should be set to p_StellarType
 *                                              (optional, default = false)
 * @return                                      Stellar type of star before switch (previous stellar type)
 */
STELLAR_TYPE Star::SwitchTo(const STELLAR_TYPE p_StellarType, bool p_SetInitialType) {

    STELLAR_TYPE stellarTypePrev = StellarType();

    // don't switch if stellarTypePrev == p_StellarType
    if (p_StellarType != StellarType()) {
        BaseStar *ptr = nullptr;

        switch (p_StellarType) {
            case STELLAR_TYPE::MS_LTE_07                                : {ptr = new MS_lte_07(*m_Star);} break;
            case STELLAR_TYPE::MS_GT_07                                 : {ptr = new MS_gt_07(*m_Star);} break;
            case STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS                   : {ptr = new CH(*m_Star);} break;
            case STELLAR_TYPE::HERTZSPRUNG_GAP                          : {ptr = new HG(*m_Star);} break;
            case STELLAR_TYPE::FIRST_GIANT_BRANCH                       : {ptr = new FGB(*m_Star);} break;
            case STELLAR_TYPE::CORE_HELIUM_BURNING                      : {ptr = new CHeB(*m_Star);} break;
            case STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH            : {ptr = new EAGB(*m_Star);} break;
            case STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH: {ptr = new TPAGB(*m_Star);} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_MS                     : {ptr = new HeMS(*m_Star);} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP        : {ptr = new HeHG(*m_Star);} break;
            case STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH           : {ptr = new HeGB(*m_Star);} break;
            case STELLAR_TYPE::HELIUM_WHITE_DWARF                       : {ptr = new HeWD(*m_Star);} break;
            case STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF                : {ptr = new COWD(*m_Star);} break;
            case STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF                  : {ptr = new ONeWD(*m_Star);} break;
            case STELLAR_TYPE::NEUTRON_STAR                             : {ptr = new NS(*m_Star);} break;
            case STELLAR_TYPE::BLACK_HOLE                               : {ptr = new BH(*m_Star);} break;
            case STELLAR_TYPE::MASSLESS_REMNANT                         : {ptr = new MR(*m_Star);} break;

            default:                                                                                                        // unknown stellar type
                // the only ways this can happen are if someone added a STELLAR_TYPE
                // and it isn't accounted for in this code, or if there is a defect in the code that causes
                // this function to be called with a bad parameter.  We should not default here, with or without
                // a warning.
                // We are here because the user chose a prescription this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing prescription or, if the missing
                // prescription is superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNKNOWN_STELLAR_TYPE);                                                                   // throw error
        }

        if (ptr) {
            delete m_Star;
            m_Star = ptr;

            if (p_SetInitialType) m_Star->SetInitialType(p_StellarType);   // CHECK THIS - NOT INITIAL, JUST TYPE????????????????????????????????
        }

        // write to switch log file if required
        // star should be evolving from one of the evolvable types (we don't want the initial switch from Star->MS).
        // check is not necessary for BSE (handled differently), but no harm
        if (utils::IsOneOf(stellarTypePrev, EVOLVABLE_TYPES) && OPTIONS->SwitchLog()) {

            LOGGING->SetSwitchParameters(m_ObjectId, ObjectType(), m_ObjectPersistence, stellarTypePrev, p_StellarType);    // store switch details to LOGGING service
            if (OPTIONS->EvolutionMode() == EVOLUTION_MODE::BSE) {                                                          // BSE?
                raise(SIGUSR1);                                                                                             // signal to BSE that switch is occurring
            }
            else {                                                                                                          // SSE
                (void)m_Star->PrintSwitchLog();                                                                             // no need for the BSE signal shenanigans - just call the function
            }
        }
    }

    return stellarTypePrev;
}


/*
 *
 * Resolve the loss of an envelope, including switching the stellar type.
 *
 * void ResolveEnvelopeLossAndSwitch()
 *
 */
void Star::ResolveEnvelopeLossAndSwitch() {
    (void)SwitchTo(m_Star->ResolveEnvelopeLoss(true));
}


/*
 * Apply mass changes if required, age the star one timestep, advance the simulation time, and update the
 * attributes of the star.
 *
 * The star's attributes (Age, Radius, Luminosity etc.) are calculated and updated as required.
 *
 * Free parameters in the update process are the star's mass (m_Mass), initial mass (m_Mass0), the star's age
 * (m_Age) and the simulation time attribute (m_Time):
 *
 *    - if required, the star's mass is changed by the amount passed as the p_DeltaMass parameter before other
 *      attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to the
 *      star's mass before the attributes of the star are calculated.
 *
 *    - if required, the star's initial mass is changed by the amount passed as the p_DeltaMass0 parameter before
 *      other attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to
 *      the star's mass before the attributes of the star are calculated.  This should be used infrequently, and
 *      is really a kludge because the Mass0 attribute in Hurley et al. 2000 was overloaded after the introduction
 *      of mass loss (see section 7.1).  We should really separate the different uses of Mass0 in the code and
 *      use a different variable - initial mass shouldn't change (other than to initial mass upon entering a
 *      stellar phase - it doesn't make a lot of sense for initial mass to change during evolution through the
 *      phase).         JR: todo
 *
 *    - if required, the star is aged by the amount passed as the p_DeltaTime parameter, and the simulation time is
 *      advanced by the same amount, before other attributes are updated.  The p_deltaTime parameter may be zero,
 *      in which case no change is made to the star's age or the physical time attribute.
 *
 *
 * Checks whether the star:
 *    - is a massless remnant (checked after applying p_DeltaMass and p_DeltaMass0, but before applying p_DeltaTime)
 *    - has become a supernova (checked after applying p_DeltaMass and p_DeltaMass0, but before applying p_DeltaTime)
 *    - should skip this phase for this timestep (checked after applying p_DeltaMass, p_DeltaMass0 and p_DeltaTime)
 *
 * If none of the above are true the star evolves on phase for the specified timestep (which may be 0, in which case
 * the star's attributes other than age are re-calculated), then the need to evolve the star off phase is checked.
 *
 * If p_DeltaMass, p_DeltaMass0 and p_DeltaTime are all passed as zero the checks for massless remnant and supernova
 * are performed (and consequential changes made), but no other changes to the star's attributes are made - unless
 * the p_ForceRecalculate parameter is set true.
 *
 * The functional return is the stellar type to which the star should evolve.  The returned stellar type is just the
 * stellar type of the star upon entry if it should remain on phase.
 *
 * If the parameter p_Switch is true the star will switch to the new stellar type before returning.
 *
 *
 * STELLAR_TYPE UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
 *                                                const double p_DeltaMass0,
 *                                                const double p_DeltaTime,
 *                                                const bool   p_Switch,
 *                                                const bool   p_ForceRecalculate)
 *
 * @param   [IN]    p_DeltaMass                 The change in mass to apply in Msol
 * @param   [IN]    p_DeltaMass0                The change in mass0 to apply in Msol
 * @param   [IN]    p_DeltaTime                 The timestep to evolve in Myr
 * @param   [IN]    p_Switch                    Specifies whether the star should switch to new stellar type before returning
 *                                              (optional, default = true)
 * @param   [IN]    p_ForceRecalculate          Specifies whether the star's attributes should be recalculated even if the three deltas are 0.0
 *                                              (optional, default = false)
 * @return                                      New stellar type for star
 */
STELLAR_TYPE Star::UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass,
                                                     const double p_DeltaMass0,
                                                     const double p_DeltaTime,
                                                     const bool   p_Switch,
                                                     const bool   p_ForceRecalculate) {

    STELLAR_TYPE stellarType = m_Star->EvolveOneTimestep(p_DeltaMass, p_DeltaMass0, p_DeltaTime, p_ForceRecalculate);

/*
    if (p_Switch && (stellarType != m_Star->StellarType())) {                               // switch to new stellar type if necessary?
        STELLAR_TYPE stellarTypePrev = SwitchTo(stellarType);                               // yes - switch
        m_Star->SetStellarTypePrev(stellarTypePrev);                                        // record previous stellar type

        // recalculate stellar attributes after switching if necessary - transition may not be continuous
        // this is a bit of a kludge just for CH -> HeMS  JR: should revisit the best way to do this
        if (stellarTypePrev == STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS && 
                stellarType == STELLAR_TYPE::NAKED_HELIUM_STAR_MS) {                        // discontinuous transition?
            if (UpdateAttributes(0.0, 0.0, true) != stellarType) {                          // yes - recalculate stellar attributes
                // JR: need to revisit this - should we actually switch here
                // (or maybe queue a switch - not sure that's even possible...)?
                SHOW_WARN(ERROR::SWITCH_NOT_TAKEN);                                         // show warning if we think we should switch again - really for diagnostics/stats (how often does this happen?)
            }
        }
    }
*/

    return stellarType;                                                                     // return new stellar type
}


/*
 * Update the mass attributes of the star (m_Mass, m_Mass0) without ageing the star, advancing the simulation
 * time, or switching the stellar type of the star.
 *
 * Apply mass changes as required and update the attributes of the star.
 *
 * The star's attributes (Radius, Luminosity etc. (not Age or simulation time)) are calculated and updated as
 * required.
 *
 * Free parameters in the update process are the star's mass (m_Mass) and initial mass (m_Mass0):
 *
 *    - if required, the star's mass is changed by the amount passed as the p_DeltaMass parameter before other
 *      attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to the
 *      star's mass before the attributes of the star are calculated.
 *
 *    - if required, the star's initial mass is changed by the amount passed as the p_DeltaMass parameter before
 *      other attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to
 *      the star's mass before the attributes of the star are calculated.  This should be used infrequently, and
 *      is really a kludge because the Mass0 attribute in Hurley et al. 2000 was overloaded after the introduction
 *      of mass loss (see section 7.1).  We should really separate the different uses of Mass0 in the code and
 *      use a different variable - initial mass shouldn't change (other than to initial mass upon entering a
 *      stellar phase - it doesn't make a lot of sense for initial mass to change during evolution through the
 *      phase).         JR: todo
 *
 *
 * See BaseStar::EvolveOneTimestep() for details of operation.
 *
 * The functional return is the stellar type to which the star should switch as a result of the update.
 * The returned stellar type is just the stellar type of the star upon entry if it remained on the current phase.
 *
 *
 * STELLAR_TYPE UpdateAttributes(const double p_DeltaMass, const double p_DeltaMass0, const bool p_ForceRecalculate)
 *
 * @param   [IN]    p_DeltaMass                 The change in mass to apply in Msol
 * @param   [IN]    p_DeltaMass0                The change in mass0 to apply in Msol
 * @param   [IN]    p_ForceRecalculate          Specifies whether the star's attributes should be recalculated even if the two deltas are 0.0
 *                                              (optional, default = false)
 * @return                                      New stellar type for star
 */
STELLAR_TYPE Star::UpdateAttributes(const double p_DeltaMass, const double p_DeltaMass0, const bool p_ForceRecalculate) {
   return m_Star->EvolveOneTimestep(p_DeltaMass, p_DeltaMass0, 0.0, p_ForceRecalculate);
}


/*
 * Age the star a single timestep - timestep is provided as parameter
 *
 * Checks whether the star:
 *    - is a massless remnant
 *    - has become a supernova
 *    - should skip this phase for this timestep
 *
 * If none of the above are true, the star evolves on phase for this timestep,
 * then the need to evolve the star off phase is checked.
 *
 * The functional return is the stellar type to which the star should evolve.  The returned stellar
 * type is just the stellar type of the star upon entry if it should remain on phase.
 *
 * If the parameter p_Switch is true the star will switch to the new stellar type before returning.
 * 
 * This function is called from SSE code (Star::EvolveOneTimeStep()), and BSE code (BaseBinaryStar::EvolveOneTimestep())
 *
 *
 * STELLAR_TYPESTELLAR_TYPE AgeOneTimestep(const double p_Dt, bool p_Switch)
 *
 * @param   [IN]    p_DeltaTime                 The timestep to evolve in Myr
 * @param   [IN]    p_Switch                    Specifies whether the star should switch to new stellar type before returning
 *                                              (optional, default = true)
 * @return                                      New stellar type for star
 */
STELLAR_TYPE Star::AgeOneTimestep(const double p_DeltaTime, bool p_Switch) {
    return m_Star->EvolveOneTimestep(0.0, 0.0, p_DeltaTime, false);
}


/*
 * Evolve the star a single timestep (SSE mode only)
 *
 * - Log pre-mass-loss details to SSE detailed output file
 * - Resolve wind mass loss
 * - Age the star a single timestep (mass changes have already been made)
 * - Return the stellar type to which the star should switch as a result of the evolution
 * - Log post-mass-loss details to SSE detailed output file
 *
 *
 * STELLAR_TYPE EvolveOneTimestep(const double p_Dt)
 *
 * @param   [IN]    p_dt                        The timestep duration (MYr)
 * @return                                      New stellar type for star
 */
STELLAR_TYPE Star::EvolveOneTimestep(const double p_dt) {
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "Star::EvolveOneTimestep(@entry), p_dt = " << p_dt << ", m_Star->Mass() = " << m_Star->Mass() << ", m_Star->Radius() = " << m_Star->Radius() << "\n";

    STELLAR_TYPE nextStellarType;

//    m_Star->UpdateDt(p_Dt);                                                                                     // advance star dt, age, and simulation time
    
    (void)m_Star->PrintDetailedOutput(m_Id, SSE_DETAILED_RECORD_TYPE::PRE_MASS_LOSS);                           // log record - pre mass loss
    
    nextStellarType = m_Star->ResolveMassLoss(p_dt);                                                            // apply wind mass loss if required
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "Star::EvolveOneTimestep(@1), p_dt = " << p_dt << ", m_Star->Mass() = " << m_Star->Mass() << ", m_Star->Radius() = " << m_Star->Radius() << "\n";

    if (nextStellarType == m_Star->StellarType())                                                               // need to switch stellar type?
        nextStellarType = m_Star->EvolveOneTimestep(0.0, 0.0, p_dt, false);                                     // no - age the star one time step - modify stellar attributes as appropriate
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "Star::EvolveOneTimestep(@2), p_dt = " << p_dt << ", m_Star->Mass() = " << m_Star->Mass() << ", m_Star->Radius() = " << m_Star->Radius() << "\n";

    (void)m_Star->PrintDetailedOutput(m_Id, SSE_DETAILED_RECORD_TYPE::POST_MASS_LOSS);                          // log record - post mass loss

    return nextStellarType;
}


/*
 * Evolve the star through its entire lifetime
 *
 *
 * EVOLUTION_STATUS Evolve(const long int p_Id)
 *
 * @param   [IN]    p_Id                        The id for this star - can be used to name logfiles for this star
 * @return                                      Status
 */
EVOLUTION_STATUS Star::Evolve(const long int p_Id) {

    EVOLUTION_STATUS evolutionStatus = EVOLUTION_STATUS::CONTINUE;                                              // default status
    STELLAR_TYPE     nextStellarType = StellarType();                                                           // next stellar type (defult is current)

    try {

        m_Id = p_Id;                                                                                            // store the id

        // evolve the star

        m_Star->CalculateGBParams();                                                                            // calculate giant branch parameters - in case for some reason star is initially not MS

        double dt = 0.0;

        (void)m_Star->PrintDetailedOutput(m_Id, SSE_DETAILED_RECORD_TYPE::INITIAL_STATE);                       // log detailed output record 

        bool usingProvidedTimesteps = false;                                                                    // using user-provided timesteps?
        DBL_VECTOR timesteps;
        if (!OPTIONS->TimestepsFileName().empty()) {                                                            // have timesteps filename?
                                                                                                                // yes
            ERROR error;
            std::tie(error, timesteps) = utils::ReadTimesteps(OPTIONS->TimestepsFileName());                    // read timesteps from file
            if (error != ERROR::NONE) {                                                                         // ok?
                THROW_ERROR(error, ERR_MSG(ERROR::NO_TIMESTEPS_READ));                                          // no - throw error - this is not what the user asked for
            }
            else usingProvidedTimesteps = true;                                                                 // have user-provided timesteps
        }

        // the evolution of a star in COMPAS is effectively represented as a finite-state machine, though we don't know
        // the number of states a priori (the number of states is not infinite, though it could  << REWRITE - IT IS INFINTE, BUT WE DON'T VISI THEM ALL
        // be intractably large, but even so, not an infinite-state machine - more like a transition system).  We do have an
        // upper limit on the number of states - set by option values for the maximum evolution
        // time and the maximum number of timesteps.  Even without those, given that we have an
        // absolute minimum timestep, and a star can't be older than the universe, we have a
        // theoretical maximum number of states.  But the actual number of state is infinite!
        // We just don't visit them all in any single run of COMPAS...  <<<<<<<<<<<<<<<<<<<<<  NEED TO DOCUMENT THIS SOMEWHERE

        // set the initial state of our state machine

        unsigned long int stepNum = 0;                                                                          // initialise step number
        while (evolutionStatus == EVOLUTION_STATUS::CONTINUE) {
            if (StellarType() == STELLAR_TYPE::MASSLESS_REMNANT) {                                              // massless remnant?
                evolutionStatus = EVOLUTION_STATUS::MASSLESS_REMNANT;                                           // set status
            }
            if (m_Star->Time() > OPTIONS->MaxEvolutionTime()) {                                                 // out of time?
                evolutionStatus = EVOLUTION_STATUS::TIMES_UP;                                                   // set status
            }
            else if (stepNum >= OPTIONS->MaxNumberOfTimestepIterations()) {                                     // out of timesteps?
                evolutionStatus = EVOLUTION_STATUS::STEPS_UP;                                                   // set status
            }
            else if (m_Star->IsOneOf(WHITE_DWARFS) || StellarType() == STELLAR_TYPE::BLACK_HOLE || 
                    (StellarType() == STELLAR_TYPE::NEUTRON_STAR && !OPTIONS->EvolvePulsars())) {               // If a WD or BH, or an NS but not evolving pulsars, we're done
                evolutionStatus = EVOLUTION_STATUS::DONE;
            }
            else if (usingProvidedTimesteps && stepNum >= timesteps.size()) {                                   // using user-provided timesteps and all consumed?
                evolutionStatus = EVOLUTION_STATUS::TIMESTEPS_EXHAUSTED;                                        // yes - set status
                SHOW_WARN(ERROR::TIMESTEPS_EXHAUSTED);                                                          // show warning
            }
            else {                                                                                              // evolve one timestep

                m_Star->SetPrevDt(dt);
            
                // check for, and process, immediate events, in the following order:
                //
                // 1. supernova
                // 2. stellar type switch
                // 3. common envelope

                // supernova
                if (IsSupernova()) {                                                                            // is star about to go supernova?
if (OPTIONS->DebugLevel() > 0) std::cout << "Processing supernova\n";
                                                                                                                // yes
                    m_Star->SetDt(ABSOLUTE_MINIMUM_TIMESTEP);                                                   // set dt
                    m_Star->AdvanceAgeAndTime();                                                                // advance age of star and simulation time
                    nextStellarType = m_Star->ResolveSupernova();                                               // resolve the supernova event
                    if (nextStellarType != StellarType()) {                                                     // stellar type change?
                        (void)SwitchTo(nextStellarType, false);                                                 // yes - switch stellar type
                        UpdateAttributes(0.0, 0.0, true);                                                       // update stellar attributes

                        // log SN details to the SSE Supernova log (BSE does its own SN printing)
                        // - only if not an ephemeral clone
                        if (OPTIONS->EvolutionMode() == EVOLUTION_MODE::SSE && m_ObjectPersistence == OBJECT_PERSISTENCE::PERMANENT) {
                            m_Star->PrintSupernovaDetails();
                        }
                    }
                    else {
                        // THIS IS AN ERROR!
if (OPTIONS->DebugLevel() > 0) std::cout << "This should not happen!\n";
                    }
                }

                // stellar type change
                else if (nextStellarType != StellarType()) {                                                    // stellar type change?
if (OPTIONS->DebugLevel() > 0) std::cout << "Processing stellar type change\n";
                                                                                                                // yes
                    m_Star->SetDt(ABSOLUTE_MINIMUM_TIMESTEP);                                                   // set dt
                    m_Star->AdvanceAgeAndTime();                                                                // advance age of star and simulation time
                    (void)SwitchTo(nextStellarType, false);                                                     // switch stellar type
                    UpdateAttributes(0.0, 0.0, true);                                                           // update stellar attributes
                }


                // ordinary timestep
                else {
if (OPTIONS->DebugLevel() > 0) std::cout << "Processing ordinary timestep\n";


                    if (usingProvidedTimesteps) {                                                                   // user-provided timesteps
                        // get new timestep
                        //   - don't quantise
                        //   - don't apply timestep multiplier
                        // (we assume user wants the timesteps in the file)
                        dt = timesteps[stepNum];
                    }
                    else {                                                                                          // not using user-provided timesteps
                        dt = QUANTISE_DT(m_Star->CalculateTimestep() * OPTIONS->TimestepMultiplier() * OPTIONS->TimestepMultipliers(static_cast<int>(StellarType()))); // new timestep; quantised
if (OPTIONS->DebugLevel() > 0) std::cout << "Star::Evolve(@1), dt = " << dt << "\n";   
                    }
                    stepNum++;                                                                                      // increment step number                                                      

if (OPTIONS->DebugLevel() > 0) std::cout << "Star::Evolve(@2), dt = " << dt << "\n";   
                    nextStellarType = EvolveOneTimestep(dt);                                                        // evolve for timestep
                    UpdateAttributes(0.0, 0.0, true);                                                               // keeps SSE in sync with BSE
if (OPTIONS->DebugLevel() > 0) std::cout << "Star::Evolve(@3), m_Star->Dt() = " << m_Star->Dt() << "\n";   

    //(void)SwitchTo(stellarType);                                                                                // switch phase if required

    //if(OPTIONS->EvolvePulsars() && m_Star->StellarType() == STELLAR_TYPE::NEUTRON_STAR){                        // if star is a neutron star and we are evolving pulsars
    //    (void)m_Star->SpinDownIsolatedPulsar(p_Dt * MYR_TO_YEAR * SECONDS_IN_YEAR);                             // update pulsar parameters due to spin down as an isolated pulsar; convert timestep to seconds for this function (uses cgs units)
    //}
                }

                if (StellarType() == STELLAR_TYPE::NEUTRON_STAR && OPTIONS->EvolvePulsars()) {                      // Pulsar output if star is a neutron star and user wants pulsar output
                    (void)m_Star->PrintPulsarEvolutionParameters(SSE_PULSAR_RECORD_TYPE::TIMESTEP_COMPLETED);       // log pulsar evolution parameters
                } 

                (void)m_Star->PrintDetailedOutput(m_Id, SSE_DETAILED_RECORD_TYPE::TIMESTEP_COMPLETED);              // log detailed output record 

                // check thresholds for system detailed output printing
                // don't use utils::Compare() here - not for time/age

                bool printSystemSnapshotRec = false;                                                                // so we only print this timestep once

                // age threshold
                // we print a record each timestep the star crosses the threshold from below
                // notes:
                //    (a) the age of individual stars can drop for various reasons (phase change, rejuvenation, winds/mass transfer, etc.),
                //        and if the age of the star drops below an age threshold, we will log another record if the star then ages beyond
                //        the same threshold (so we might log several records for the star crossing the same threshold if the age of the
                //        star oscillates around the threshold)
                //    (b) we will print multiple records for exceeding the age threshold if the constituent stars exceed the age threshold
                //        at different timesteps (likely)
                for (size_t threshold = 0; threshold < OPTIONS->SystemSnapshotAgeThresholds().size(); threshold++) {// for each system detailed output age threshold

                    double thresholdValue = OPTIONS->SystemSnapshotAgeThresholds(threshold);                        // this threshold value
      
                    // flag need to print (log) system snapshot record
                    // we don't want to print multiple records for the same timestep, so we flag need rather than print here
                    printSystemSnapshotRec |= m_SystemSnapshotAgeFlags[threshold] < 0.0 && m_Star->Age() >= thresholdValue;

                    // record the current age of the star in the threshold flag - this is how we check for re-crossing a threshold
                    // if the age of the star has dropped below the threshold value, we reset the threshold flag for the star
                    // the check will fail if the star hasn't crossed the threshold already, but the flag will be -1.0 anyway
                    m_SystemSnapshotAgeFlags[threshold] = (m_Star->Age() < thresholdValue) ? -1.0 : m_Star->Age();
                }

                // time threshold
                // we print a record at the first timestep that the simulation time exceeds the time threshold
                for (size_t threshold = 0; threshold < OPTIONS->SystemSnapshotTimeThresholds().size(); threshold++) { // for each system snapshott time threshold
                    if (!m_SystemSnapshotTimeFlags[threshold] && m_Star->Time() >= OPTIONS->SystemSnapshotTimeThresholds(threshold)) { // need to action?
                        m_SystemSnapshotTimeFlags[threshold] = true;                                                // yes, flag action taken
                        printSystemSnapshotRec               = true;                                                // flag need to print (log) system snapshot record
                    }
                }

                if (printSystemSnapshotRec) (void)m_Star->PrintSystemSnapshotLog();                                 // print (log) system record record if necessary
                
if (OPTIONS->DebugLevel() > 0) { std::cout << "Continue?\n"; std::string tmp; std::cin >> tmp; }
            }
        }




        if (usingProvidedTimesteps && timesteps.size() > stepNum) {                                             // all user-defined timesteps consumed?
            evolutionStatus = EVOLUTION_STATUS::TIMESTEPS_NOT_CONSUMED;                                         // no - set status
            SHOW_WARN(ERROR::TIMESTEPS_NOT_CONSUMED);                                                           // show warning
        }

        (void)m_Star->PrintDetailedOutput(m_Id, SSE_DETAILED_RECORD_TYPE::FINAL_STATE);                         // log detailed output record 

        // if we trapped a floating-point error we set the star's error value to indicate a
        // floating-point error occurred, but we don't terminate evolution (we can only have
        // floating-point errors trapped here if the user has not activated the floating-point
        // error instrumentation.  i.e --fp-error-mode OFF)
        // Set the error here so that users know that a floating-point error occurred, even though
        // the evolution of the star was not terminated because an error occurred.

        if (fetestexcept(FE_DIVBYZERO) ||
            fetestexcept(FE_INVALID)   ||
            fetestexcept(FE_OVERFLOW)  ||
            fetestexcept(FE_UNDERFLOW)) m_Star->SetError(ERROR::FLOATING_POINT_ERROR);                          // floating-point error

        feclearexcept(FE_ALL_EXCEPT);                                                                           // clear all FE traps
    }
    catch (const std::runtime_error& e) {                                                                       // catch runtime exceptions
        // anything we catch here should not already have been displayed to the user,
        // so set the error value, display the error, and flag termination (do not rethrow the error)
        if (std::string(e.what()) == "FPE") m_Star->SetError(ERROR::FLOATING_POINT_ERROR);                      // floating-point error
        else                                m_Star->SetError(ERROR::ERROR);                                     // unspecified error
        SHOW_ERROR(m_Star->Error());                                                                            // display error (don't throw here - handled by returning status)
        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                              // evolution terminated
    }
    catch (int e) {
        // anything we catch here should already have been displayed to the user,
        // so just ensure error value is set and flag termination (do not rethrow the error)
        if (e != static_cast<int>(ERROR::NONE)) m_Star->SetError(static_cast<ERROR>(e));                        // specified errpr
        else                                    m_Star->SetError(ERROR::ERROR);                                 // unspecified error
        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                              // evolution terminated
    }
    catch (...) {
        // anything we catch here should not already have been displayed to the user,
        // so set the error value, display the error, and flag termination (do not rethrow the error)
        m_Star->SetError(ERROR::ERROR);                                                                         // unspecified error
        SHOW_ERROR(m_Star->Error());                                                                            // display error (don't throw here - handled by returning status)
        evolutionStatus = EVOLUTION_STATUS::ERROR;                                                              // evolution terminated
    }

    m_Star->SetEvolutionStatus(evolutionStatus);                                                                // set evolution final outcome for star

    (void)m_Star->PrintSystemParameters();                                                                      // log system parameters record

    return evolutionStatus;
}
