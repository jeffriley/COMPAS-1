#ifndef __BaseStar_h__
#define __BaseStar_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"
#include "vector3d.h"

#include "Rand.h"
#include "Options.h"
#include "Log.h"
#include "Errors.h"

#include "BH.h"

#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>

#include <boost/math/distributions.hpp>


class BaseStar {
    

    class StarState: public State::StarState;
        
public:
    
//    BaseStar();


    BaseStar::BaseStar(const unsigned long int p_RandomSeed, 
                       const double            p_Metallicity,
                       const double            p_Mass,
                       const KickParameters    p_KickParameters,
                       const double            p_AngularFrequency) {

                        
    // RECALC GLOBALS HERE IF NECESSARY <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


        // initialise as much of the first state of the star as possible
        // The state history stack has already been initialised - the single
        // entry is the first state of the star

        // set state attributes supplied as parameters
        m_State.SetAngularFrequency(p_AngularFreqency);
        m_State.SetMass(p_Mass);
        m_State.SetMassEffective(p_Mass);
        m_State.SetMetallicity(p_Metallicity);
        m_State.SetKickParameters(p_KickParameters);
        m_State.SetRandomSeed(p_RandomSeed);
        
        ////////////////m_State.SetStellarType(p_StellarType);

        // set, or calculate where possible/necessary, remaining state attributes

        m_State.SetError(ERROR::NONE);
        m_State.SetEvolutionStatus(EVOLUTION_STATUS::CONTINUE);

        m_State.SetLuminosity(CalculateLuminosityAtZAMS(p_Mass));
        m_State.SetLuminosityEffective(m_State.Luminosity());

        m_State.SetRadius(CalculateRadiusAtZAMS(p_Mass));
        m_State.SetRadiusEffective(m_State.Radius());

        m_State.SetTemperature(CalculateTemperatureOnPhase_Static(m_State.Luminosity(), m_State.Radius()));

        m_State.SetTime(0.0);
        m_State.SetTau(0.0);


//m_AngularMomentum                          = CalculateMomentOfInertiaAU() * m_OmegaZAMS;





    //stateZero.coreMass            = 0.0;                        // to be calculated
    //stateZero.coreMassMS          = 0.0;                        // to be calculated
    //stateZero.coreMassCO          = 0.0;                        // to be calculated
    //stateZero.coreMassHe          = 0.0;                        // to be calculated
    
  


    //stateZero.dt                  = 0.0;                        // initially
    //stateZero.dMdt                = 0.0;                        // initially


    }
        
    virtual ~BaseStar() {}

 

    
    
    // object identifiers - all classes have these
    OBJECT_ID                   ObjectId() const                                        { return m_ObjectId; }
    OBJECT_TYPE                 ObjectType() const                                      { return OBJECT_TYPE::BASE_STAR; }
    OBJECT_PERSISTENCE          ObjectPersistence() const                               { return m_ObjectPersistence; }
    
    
    // Getters - alphabetically
    //
    // Short (one or two line) implementations that don't require a lot of documentation (i.e. that won't clutter here), are in here (and are implied inline)
    // Short implementations that can still be inlined, but that require more documentation (and would clutter here), are below (after class declaration)
    // Longer implementations that can still be inlined (that would clutter here) are below (after class declaration)
    // All other implementations are in the cpp file
    //
    // *** getter implementations should generally be very short



    inline double           Age() const                                 { return m_StateHistory.CurrentState.Age(); }
    inline double           AngularFrequency() const                    { return m_StateHistory.CurrentState.AngularFrequency(); }
    inline double           COCoreMass() const                          { return m_StateHistory.CurrentState.COCoreMass(); }
    inline double           CoreMass() const                            { return m_StateHistory.CurrentState.CoreMass(); }
    inline double           dt() const                                  { return m_StateHistory.CurrentState.dt(); }
    inline double           dtPrev() const                              { return m_StateHistory.PreviousState.dt(); }
    inline double           HAbundanceCore() const                      { return m_StateHistory.CurrentState.HAbundanceCore(); }
    inline double           HAbundanceSurface() const                   { return m_StateHistory.CurrentState.HAbundanceSurface(); }
    inline double           HeAbundanceCore() const                     { return m_StateHistory.CurrentState.HeAbundanceCore(); }
    inline double           HeAbundanceSurface() const                  { return m_StateHistory.CurrentState.HeAbundanceSurface();}
    inline double           HeCoreMass() const                          { return m_StateHistory.CurrentState.HeCoreMass(); }
    inline double           Luminosity() const                          { return m_StateHistory.CurrentState.Luminosity(); }
    inline double           LuminosityPrev() const                      { return m_StateHistory.PreviousState.Luminosity(); }
    inline double           LZAMS() const                               { return m_StateHistory.ZAMSState.Luminosity(); }
    inline double           Mass() const                                { return m_StateHistory.CurrentState.Mass(); }
    inline double           Mass0() const                               { return m_StateHistory.CurrentState.MassEffectiveInitial(); }
    inline double           MassPrev() const                            { return m_StateHistory.PreviousState.Mass(); }
    inline double           MZAMS() const                               { return m_StateHistory.ZAMSState.Mass(); }
    inline double           Radius() const                              { return m_StateHistory.CurrentState.Radius(); }
    inline double           RadiusPrev() const                          { return m_StateHistory.PreviousState.Radius(); }
    inline double           RZAMS() const                               { return m_StateHistory.ZAMSState.Radius(); }
    inline double           StellarType() const                         { return m_StateHistory.CurrentState.StellarType(); }
    inline double           StellarTypePrev() const                     { return m_StateHistory.Previoustate.StellarType(); }
    inline double           Tau() const                                 { return m_StateHistory.CurrentState.Tau(); }
    inline double           Temperature() const                         { return m_StateHistory.CurrentState.Temperature(); }
    inline double           Time() const                                { return m_StateHistory.CurrentState.Time(); }
    inline double           TZAMS() const                               { return m_StateHistory.ZAMSState.Temperature(); }

    inline double           GBparam(HURLEY_GBP p_GBparam) const         { return m_StateHistory.CurrentState.GBparams(p_GBparam); }
    inline DBL_VECTOR       GBparams() const                            { return m_StateHistory.CurrentState.GBparams(); }
    inline double           Timescale(TIMESCALE p_Timescale) const      { return m_StateHistory.CurrentState.Timescales(p_Timescale); }
    inline DBL_VECTOR       Timescales() const                          { return m_StateHistory.CurrentState.Timescales(); }


    inline EVOLUTION_STATUS EvolutionStatus() const                                 { return m_EvolutionStatus; }



    inline void             SetPersistence(const OBJECT_PERSISTENCE p_Persistence)  { m_ObjectPersistence = p_Persistence; }
    inline void             SetStellarType(const STELLAR_TYPE p_StellarType) const  { m_StateHistory.CurrentState.SetStellarType(p_StellarType); }
    inline void             SetStartingType(const STELLAR_TYPE p_StellarType) const { m_StateHistory.StartState.SetStellarType(p_StellarType); }





    double                      AngularMomentum() const                                 { return m_AngularMomentum; }
    bool                        CHonMS() const                                          { return m_CHE; }
    int                         DominantMassLossType() const                            { return static_cast<int>(m_DominantMassLossType); }
    bool                        EnvelopeJustExpelledByPulsations() const                { return m_EnvelopeJustExpelledByPulsations; }
    ERROR                       Error() const                                           { return m_Error; }
    bool                        ExperiencedAIC() const                                  { return utils::SNEventType(m_SupernovaDetails.events.past) == SN_EVENT::AIC; }
    bool                        ExperiencedCCSN() const                                 { return utils::SNEventType(m_SupernovaDetails.events.past) == SN_EVENT::CCSN; }
    bool                        ExperiencedHeSD() const                                 { return utils::SNEventType(m_SupernovaDetails.events.past) == SN_EVENT::HeSD; }
    bool                        ExperiencedECSN() const                                 { return utils::SNEventType(m_SupernovaDetails.events.past) == SN_EVENT::ECSN; }
    bool                        ExperiencedPISN() const                                 { return utils::SNEventType(m_SupernovaDetails.events.past) == SN_EVENT::PISN; }
    bool                        ExperiencedPPISN() const                                { return utils::SNEventType(m_SupernovaDetails.events.past) == SN_EVENT::PPISN; }
    bool                        ExperiencedSNIA() const                                 { return utils::SNEventType(m_SupernovaDetails.events.past) == SN_EVENT::SNIA; }
    SN_EVENT                    ExperiencedSN_Type() const                              { return utils::SNEventType(m_SupernovaDetails.events.past); }
    bool                        ExperiencedUSSN() const                                 { return (m_SupernovaDetails.events.past & SN_EVENT::USSN) == SN_EVENT::USSN; }
    std::string                 MassTransferDonorHistory() const; // implementation in cpp file

    double                      InitialHeliumAbundance() const                          { return m_InitialHeliumAbundance; }
    double                      InitialHydrogenAbundance() const                        { return m_InitialHydrogenAbundance; }
    double                      InitialMainSequenceCoreMass() const                     { return m_InitialMainSequenceCoreMass; }
    bool                        IsAIC() const                                           { return (m_SupernovaDetails.events.current & SN_EVENT::AIC) == SN_EVENT::AIC; }
    bool                        IsCCSN() const                                          { return (m_SupernovaDetails.events.current & SN_EVENT::CCSN) == SN_EVENT::CCSN; }
    virtual bool                IsDegenerate() const                                    { return false; }                                                           // virtual - default is not degenerate - White Dwarfs, NS and BH are degenerate
    bool                        IsECSN() const                                          { return (m_SupernovaDetails.events.current & SN_EVENT::ECSN) == SN_EVENT::ECSN; }
    bool                        IsHeSD() const                                          { return (m_SupernovaDetails.events.current & SN_EVENT::HeSD) == SN_EVENT::HeSD; }
    bool                        IsOneOf(const STELLAR_TYPE_LIST p_List) const;  // inline implementation below
    bool                        IsPISN() const                                          { return (m_SupernovaDetails.events.current & SN_EVENT::PISN) == SN_EVENT::PISN; }
    bool                        IsPPISN() const                                         { return (m_SupernovaDetails.events.current & SN_EVENT::PPISN) == SN_EVENT::PPISN; }
    bool                        IsSN_NONE() const                                       { return m_SupernovaDetails.events.current == SN_EVENT::NONE; }
    bool                        IsSNIA() const                                          { return (m_SupernovaDetails.events.current & SN_EVENT::SNIA) == SN_EVENT::SNIA; }
    virtual bool                IsSupernova() const                                     { return false; }                                                           // virtual - default value
    bool                        IsUSSN() const                                          { return (m_SupernovaDetails.events.current & SN_EVENT::USSN) == SN_EVENT::USSN; }
    bool                        LBV_PhaseFlag() const                                   { return m_LBVphaseFlag; }
    double                      MainSequenceCoreMass() const                            { return m_MainSequenceCoreMass; }
    ST_VECTOR                   MassTransferDonorHistory() const                        { return m_MassTransferDonorHistory; }
    double                      Mdot() const                                            { return m_Mdot; }
    double                      Omega() const                                           { return m_AngularMomentum / CalculateMomentOfInertiaAU(); }
    double                      OmegaCHE() const                                        { return m_OmegaCHE; }
    double                      OmegaBreak() const                                      { return CalculateOmegaBreak(); }
    double                      OmegaZAMS() const                                       { return m_OmegaZAMS; }
    COMPAS_VARIABLE             PropertyValue(const T_ANY_PROPERTY p_Property) const; // implmentation in cpp file
    double                      PulsarMagneticField() const                             { return m_PulsarDetails.magneticField; }
    double                      PulsarSpinPeriod() const                                { return m_PulsarDetails.spinPeriod; }
    double                      PulsarSpinFrequency() const                             { return m_PulsarDetails.spinFrequency; }
    double                      PulsarSpinDownRate() const                              { return m_PulsarDetails.spinDownRate; }
    double                      PulsarBirthPeriod() const                               { return m_PulsarDetails.birthPeriod; }
    double                      PulsarBirthSpinDownRate() const                         { return m_PulsarDetails.birthSpinDownRate; }
    unsigned long int           RandomSeed() const                                      { return m_RandomSeed; }
    double                      SN_CoreMassAtCOFormation() const                        { return m_SupernovaDetails.coreMassAtCOFormation; }
    double                      SN_CoreRadiusAtCOFormation() const                      { return m_SupernovaDetails.coreRadiusAtCOFormation; }
    double                      SN_COCoreMassAtCOFormation() const                      { return m_SupernovaDetails.COCoreMassAtCOFormation; }
    SupernovaDetailsT           SN_Details() const                                      { return m_SupernovaDetails; }
    double                      SN_DrawnKickMagnitude() const                           { return m_SupernovaDetails.drawnKickMagnitude; }
    double                      SN_EccentricAnomaly() const                             { return m_SupernovaDetails.eccentricAnomaly; }
    double                      SN_FallbackFraction() const                             { return m_SupernovaDetails.fallbackFraction; }
    double                      SN_HeCoreMassAtCOFormation() const                      { return m_SupernovaDetails.HeCoreMassAtCOFormation; }
    bool                        SN_IsHydrogenPoor() const                               { return m_SupernovaDetails.isHydrogenPoor; }
    double                      SN_KickMagnitude() const                                { return m_SupernovaDetails.kickMagnitude; }
    double                      SN_MeanAnomaly() const                                  { return m_SupernovaDetails.meanAnomaly; }
    double                      SN_Phi() const                                          { return m_SupernovaDetails.phi; }
    double                      SN_RocketKickMagnitude() const                          { return m_SupernovaDetails.rocketKickMagnitude; }
    double                      SN_RocketKickPhi() const                                { return m_SupernovaDetails.rocketKickPhi; }
    double                      SN_RocketKickTheta() const                              { return m_SupernovaDetails.rocketKickTheta; }
    double                      SN_TotalMassAtCOFormation() const                       { return m_SupernovaDetails.totalMassAtCOFormation; }
    double                      SN_TotalRadiusAtCOFormation() const                     { return m_SupernovaDetails.totalRadiusAtCOFormation; }
    double                      SN_TrueAnomaly() const                                  { return m_SupernovaDetails.trueAnomaly; }
    double                      SN_Theta() const                                        { return m_SupernovaDetails.theta; }
    SN_EVENT                    SN_Type() const                                         { return utils::SNEventType(m_SupernovaDetails.events.current); }
    double                      SN_KickMagnitudeRandom() const                          { return m_SupernovaDetails.kickMagnitudeRandom; }
    double                      Speed() const                                           { return m_ComponentVelocity.Magnitude(); }
    COMPAS_VARIABLE             StellarPropertyValue(const T_ANY_PROPERTY p_Property) const; // implmentation in cpp file
    double                      TotalMassLossRate() const                               { return m_TotalMassLossRate; }
    double                      VelocityX() const                                       { return m_ComponentVelocity.xValue(); }
    double                      VelocityY() const                                       { return m_ComponentVelocity.yValue(); }
    double                      VelocityZ() const                                       { return m_ComponentVelocity.zValue(); }
    virtual ACCRETION_REGIME    WhiteDwarfAccretionRegime() const                       { return ACCRETION_REGIME::ZERO; }                                          // virtual - default value
    double                      XExponent() const                                       { return m_XExponent; }
    
    
    // Setters - alphabetically
    //
    // Short (one or two line) implementations that don't require a lot of documentation (i.e. that won't clutter here), are in here (and are implied inline)
    // Short implementations that can still be inlined, but that require more documentation (and would clutter here), are below (after class declaration)
    // Longer implementations that can still be inlined (that would clutter here) are below (after class declaration)
    // All other implementations are in the cpp file
    //
    // *** setter implementations should generally be very short


 

    void                        AdvanceAgeAndTime(const double p_dt)                            {
                                                                                                    if (p_dt > 0.0) {                                                               // Only if delta > 0.0 (don't use utils::Compare() here)
                                                                                                        m_Age  += p_dt;                                                             // Advance age of star
                                                                                                        m_Time += p_dt;                                                             // Advance simulation time
                                                                                                    }
                                                                                                };
    void                        AdvanceAgeAndTime()                                             { AdvanceAgeAndTime(m_dt); }                                        // Use class member variable

    void                        SetAngularMomentum(double p_AngularMomentum)                    { m_AngularMomentum = std::max(p_AngularMomentum, 0.0); }
//    void                        SetError(const ERROR p_Error)                                   { m_Error = p_Error; }
//    void                        SetObjectId(const OBJECT_ID p_ObjectId)                         { m_ObjectId = p_ObjectId; }
    
    void                        SetOmega(double p_Omega)                                        { SetAngularMomentum(CalculateMomentOfInertiaAU() * p_Omega); }
    

GNU_CONST inline SN_EVENT AddSNEvent(const SN_EVENT p_SNEvent, const SN_EVENT p_ExistingSNEvent) { return p_ExistingSNEvent |= p_SNEvent; }
GNU_CONST inline SN_EVENT ClearSNEvent(const SN_EVENT p_SNEvent) { return SN_EVENT::NONE; }
    

    void                        SetEvolutionStatus(const EVOLUTION_STATUS p_EvolutionStatus)    { m_EvolutionStatus = p_EvolutionStatus; }                          // Set evolution status (typically final outcome) for star

    void                        SetDt(const double p_dt)                                        { m_dt = std::max(0.0, p_dt); }                                     // Set timestep - ensure >= 0.0
    void                        SetPrevDt(const double p_dt)                                    { m_dtPrev = p_dt; }                                                // Set previous timestep - don't clamp - preserve actuality
  
    virtual void                UpdateEffectiveZAMSLandR(){}                                                                                                        // Virtual - default is NO-OP
    

    // Class member functions
    //
    // Short (one or two line) implementations that don't require a lot of documentation (i.e. that won't clutter here),
    // are in here (and are implied inline).
    //
    // Short implementations that can still be inlined, but that require more documentation (and would clutter here),
    // are below (after class declaration).
    //
    // Longer implementations that can still be inlined (that would clutter here) are below (after class declaration).
    //
    // All other implementations are in the cpp file
    //
    //
    // Function names should indicate the purpose of the function (and should genarally begin with a verb).
    //
    // Functions beginning with:
    //
    //    - "Calculate" should calculate and return value(s) 
    //    - "Determine" should determine e.g. current state, envelope type, etc., and return value(s)
    //    - "Draw" should draw a value, or values, from a distribution and return value(s)
    //
    // Most (almost all) functions should do some operation ("Calculate", "Determine", "Draw", etc.),
    // and should return the result(s) of that operation - they should not modify class member or state
    // variables.  class member and/or state variables should be modified in as few places in the code
    // as possible - this helps:
    //
    //     (a) to keep the code clean,
    //     (b) developers to be clear about what functions do and any side-effects they may have,
    //     (c) to facilitate our "one (major) things at a time" paradigm
    //
    // Functions that do not update (modify) class member variables should be declared const
    

GNU CONST double CalculateBindingEnergy(const double p_Radius, const double p_CoreMass, const double p_EnvMass, const double p_Lambda) const;
    double       CalculateBindingEnergy(const double p_Lambda) const { return CalculateBindingEnergy(m_Radius, m_CoreMass, m_Mass - m_CoreMass, p_Lambda); }

GNU_CONST virtual inline double CalculateConvectiveCoreMass() const {return 0.0;}

GNU_CONST virtual inline double CalculateConvectiveCoreRadius(const double p_Radius, const double p_Tau) const {return 0.0;}



GNU_CONST CalculateConvectiveEnvelopeBindingEnergy(const double p_Mass, const double p_Radius, const double p_EnvMass, const double p_Lambda) const
GNU_PURE  double CalculateConvectiveEnvelopeLambda_Picker(const double p_Mass, const double p_EnvMass, const double p_EnvMassMax) const;
    
    
    
virtual DBL_DBL CalculateConvectiveEnvelopeMass() const { return std::make_tuple(0.0, 0.0); }





    
    
    double                      CalculateEddyTurnoverTimescale() const;
    



    double                      CalculateMassLossValues(double p_Dt, const bool p_UpdateMDot = false);

    
    double                      CalculateNuclearMassLossRate()                                                  { return m_Mass / CalculateTimescale_RadialExpansion_DuringMT(); }
        
    COMPAS_PURE static double CalculateAngularFrequencyCHE_Static(const double p_MZAMS) const;

    double                      CalculateRadialChange() const                                                   { return (utils::Compare(m_RadiusPrev,0) <= 0)? 0 : std::abs(m_Radius - m_RadiusPrev) / m_RadiusPrev; } // Return fractional radial change (if previous radius is negative or zero, return 0 to avoid NaN

    double                      CalculateTimescale_RadialExpansion_DuringMT();

    virtual double              CalculateConvectiveEnvelopeRadialExtent() const                                 { return 0.0; }                                                     // Default for stars with no convective envelope
    virtual double              CalculateRadiusOnMassChange(double p_dM)                                        { return Radius(); }                                                // NO-OP



    double                      CalculateSNkickMagnitude(const double p_EjectaMass, const double p_RemnantMass, const STELLAR_TYPE p_StellarType);


 

    virtual double              CalculateTAMSCoreMass() const                                                   { return 0.0; }                                                     // Except MS stars

    double                      CalculateTimestep();



    virtual inline ENVELOPE DetermineEnvelopeType() const { return ENVELOPE::REMNANT; } // Default is REMNANT - but should never be called

    STELLAR_TYPE                EvolveOneTimestep(const double p_dM, const double p_dM0, const double p_dt, const bool p_ForceRecalculate = false);
    
    void                        HaltWinds()                                                                     { m_Mdot = 0.0; }                                                   // Disable wind mass loss in current time step (e.g., if star is a donor or accretor in a RLOF episode)










    void                        ResetEnvelopeExpulsationByPulsations()                                          { m_EnvelopeJustExpelledByPulsations = false; }


    virtual STELLAR_TYPE        ResolveEnvelopeLoss(bool p_Force = false)                                       { return m_StellarType; }
    virtual STELLAR_TYPE        ResolveMassLoss(const double p_dt);
    virtual STELLAR_TYPE        ResolveMassLossHurley(const double p_dt);
    virtual void                ResolveShellChange(const double p_AccretedMass) { }                                                                                                 // Default does nothing, use inheritance for WDs.
    virtual STELLAR_TYPE        ResolveSupernova()                                                              { return m_StellarType; }                                           // Default is NO-OP
       
    void                        SetStellarTypePrev(const STELLAR_TYPE p_StellarTypePrev)                        { m_StellarTypePrev = p_StellarTypePrev; }
    
    virtual bool                ShouldEnvelopeBeExpelledByPulsations() const                                    { return false; }                                                   // Default is that there is no envelope expulsion by pulsations

    virtual void                SpinDownIsolatedPulsar(const double p_Stepsize) { }                                                                                                 // Default is NO-OP
    
    virtual void                UpdateGBparams(const double p_Mass, DBL_VECTOR &p_GBparams)                  { }                                                                                  // Default is NO-OP
    virtual void                UpdateGBparams()                                                             { UpdateGBparams(m_Mass0, m_GBparams); }                         // Use class member variables

    


    
    
    
    
    // Default is NO-OP
    virtual void                UpdateMagneticFieldAndSpin(const bool p_CommonEnvelope, const bool p_RecyclesNS, const double p_Stepsize, const double p_MassGainPerTimeStep, const double p_Epsilon) { } // Default is NO-OP
    virtual void                UpdateMainSequenceCoreMass(const double p_Dt, const double p_TotalMassLossRate) { }                                                                 // Set core mass for Main Sequence stars; default is NO-OP
    virtual void                UpdateTotalMassLossRate(const double p_MassLossRate)                            { m_TotalMassLossRate = p_MassLossRate; }                           // m_TotalMassLossRate = -m_Mdot in SSE, during a mass transfer episode m_TotalMassLossRate = m_MassLossRateInRLOF
    
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

protected:

    // Protected member variables - grouped, then alphabetical

    // Object identifiers - all classes have these
    OBJECT_ID               m_ObjectId;                                 // Instantiated object's unique object id
    OBJECT_PERSISTENCE      m_ObjectPersistence;                        // Instantiated object's persistence (permanent or ephemeral)

    StarState               m_State;
    StateHistory            m_StateHistory(m_State);




    // Protected member function prototypes - alphabetically (implementations are in the cpp file)
    // See below (after class declaration) for inline functions declared in this header file

    double              CalculateAlpha1() const;                        // inline implementation below
    double              CalculateAlpha3() const;                        // inline implementation below
    double              CalculateAlpha4() const;                        // inline implementation below

    void                CalculateAnCoefficients(DBL_VECTOR &p_An, DBL_VECTOR &p_LConstants, DBL_VECTOR &p_RConstants, DBL_VECTOR &p_GammaConstants) const;



    void                CalculateBnCoefficients(DBL_VECTOR &p_Bn) const;



    double       CalculateConvectiveEnvelopeBindingEnergy(const double p_Lambda) const { 
 double convectiveEnvMass;
 std::tie(convectiveEnvMass, std::ignore) = CalculateConvectiveEnvelopeMass();
 return CalculateConvectiveEnvelopeBindingEnergy(m_Mass, m_Radius, convectiveEnvMass, p_Lambda);
 }

                                                     // Default is NO-OP


    virtual double      CalculateEddingtonCriticalRate() const                                                          { return 2.08E-3 / 1.7 * m_Radius * MYR_TO_YEAR * OPTIONS->EddingtonAccretionFactor() ; } // Hurley+, 2002, Eq. (67)






    double              CalculateGBRadiusXExponent() const;

                                                 // Default is NO-OP




    
GNU_CONST double inline BaseStar::CalculateMLfractionWR(const double p_HeAbundanceSurface) const;




    void                CalculateLCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_LCoefficients) const;






COMPAS_PURE double CalculateLuminosityAtBAGB_Hurley2000(double p_Mass) const;



    double              CalculateLuminosityAtBAGB(double p_Mass) const;         // inline implementation below

                                                 // Default is NO-OP



GNU_CONST inline double BaseStar::CalculateTimescale_MassChange(
    const STELLAR_TYPE p_StellarType,
    const STELLAR_TYPE p_StellarTypePrev,
    const double       p_Mass,
    const double       p_MassPrev,
    const double       p_dtPrev
) const;


    
    void                CalculateMassCutoffs(DBL_VECTOR &p_MassCutoffs) const;





























   
    
    
    
    
    



    double              CalculateMassLossRateWolfRayetTemperatureCorrectionSander2023(const double p_Mdot) const;



// JR : DON'T FORGET "OVERRIDE" ON DERIVED FUNCTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



virtual double CalculateMomentOfInertia() const { return (0.1 * (m_Mass) * m_Radius * m_Radius); } // Defaults to MS. k2 = 0.1 as defined in Hurley et al. 2000, after eq 109
    virtual double              CalculateMomentOfInertiaAU() const                                              { return CalculateMomentOfInertia() * RSOL_TO_AU * RSOL_TO_AU; }


    GNU_CONST double CalculateEddingtonLuminosity(const double p_Mass, const double p_HeAbundanceSurface) const;








GNU_CONST double CalculateTemperatureOnPhase() const;
GNU_CONST double CalculateTemperatureOnPhase(const double p_Luminosity, const double p_Radius) const;





// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ********************************************************* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


    ////////////////////////////////////////
    //   ABUNDANCE                        //
    ////////////////////////////////////////

    virtual double CalculateHAbundanceCore(const double p_Tau) const     { return HAbundanceCore(); }
    virtual double CalculateHAbundanceSurface(const double p_Tau) const  { return HAbundanceSurface(); }
    virtual double CalculateHeAbundanceCore(const double p_Tau) const    { return HeAbundanceCore(); }
    virtual double CalculateHeAbundanceSurface(const double p_Tau) const { return HeAbundanceSurface(); }


    ////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES   //
    ////////////////////////////////////////

    virtual double CalculateAgeAfterMassLoss() const;
    virtual inline double CalculateAgeAfterMassLoss_Hurley2000() const { return Age(); } 

    GNU_CONST   double CalculateLifetimeToBAGB_Hurley2000(const double p_tHeI, const double p_tHe) const;
    COMPAS_PURE double CalculateLifetimeToBGB_Hurley2000_Static(const double p_Mass) const;

    virtual double CalculateTau() const;
    virtual inline double CalculateTau_Hurley2000() const { return Tau(); }

    virtual DBL_VECTOR CalculateTimescales() const;
    virtual DBL_VECTOR CalculateTimescales_Hurley2000() const { return Timescales(); }

    inline double CalculateTimescale_Dynamical() const { 
        return CalculateTimescale_Dynamical_Kalogera1996(.Mass(), Radius());
    }

    GNU_CONST double CalculateTimescale_Dynamical_Kalogera1996(const double p_Mass, const double p_Radius) const;

    virtual inline double CalculateTimescale_Thermal() const { return CalculateTimescale_Thermal(Mass(), Radius(), Luminosity(), CoreMass()); }

    GNU_CONST double CalculateTimescale_Thermal(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const;


    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////

    virtual double CalculateLuminosity() const;
    virtual double CalculateLuminosityAtPhaseEnd() const;

    COMPAS_PURE double CalculateLuminosityAtZAMS_Tout1996(const double p_MZAMS) const;

    virtual inline double CalculateLuminosityAtPhaseEnd_Hurley2000() const { return Luminosity(); } 
    virtual inline double CalculateLuminosity_Hurley2000() const { return Luminosity(); } 

    GNU_CONST double CalculateLuminosity_Hurley2000(const double p_CoreMass, const DBL_VECTOR& p_GBparams) const;


    ////////////////////////////////////////
    //   MASS                             //
    ////////////////////////////////////////

    virtual double CalculateCoreMass() const;
    virtual inline double CalculateCoreMass_Hurley2000() const { return CoreMass(); } 

    GNU_CONST static double CalculateCoreMass_Hurley2000_Static(const double p_Luminosity, const DBL_VECTOR& p_GBparams); // <<<<< WHAT'S THIS FOR?????

    virtual double CalculateCOCoreMass() const;
    virtual inline double CalculateCOCoreMass_Hurley2000() const { return COCoreMass(); } 

    virtual double CalculateMass0() const;
    virtual inline double CalculateMass0_Hurley2000() const { return Mass0(); }

    virtual double CalculateHeCoreMass() const;
    virtual inline double CalculateHeCoreMass_Hurley2000() const { return HeCoreMass(); } 


    ////////////////////////////////////////
    //   MASS LOSS                        //
    ////////////////////////////////////////

    virtual MASS_LOSS_T CalculateMassLossRate() const;

    // mass loss LBV
    COMPAS_PURE MASS_LOSS_T CalculateMLrateLBV(const double p_Radius, const double p_Luminosity, const LBV_MASS_LOSS_PRESCRIPTION p_LBV_MLprescription) const;
    COMPAS_PURE MASS_LOSS_T CalculateMLrateLBV_Belczynski2010() const;
    GNU_CONST   MASS_LOSS_T CalculateMLrateLBV_Hurley2000(const double p_Luminosity, const double p_HDlimitfactor) const;

    // mass loss OB
    COMPAS_PURE MASS_LOSS_T CalculateMLrateOB(
        const double                    p_Mass,
        const double                    p_Luminosity,
        const double                    p_Temperature,
        const OB_MASS_LOSS_PRESCRIPTION p_OB_MLprescription
    ) const;
    COMPAS_PURE MASS_LOSS_T CalculateMLrateOB_Bjorklund2022(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE MASS_LOSS_T CalculateMLrateOB_Krticka2018(const double p_Luminosity) const;
    COMPAS_PURE MASS_LOSS_T CalculateMLrateOB_Vink2001(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE MASS_LOSS_T CalculateMLrateOB_VinkSander2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;

    // mass loss PRESCRIPTIONS
    COMPAS_PURE virtual MASS_LOSS_T CalculateMLrate_Belczynski2010(
        const double p_Mass,
        const double p_Radius,
        const double p_Luminosity,
        const double p_Temperature,
        const double p_PerturbationMu,
        const double p_HeAbundanceSurface
    ) const;

    virtual COMPAS_PURE MASS_LOSS_T CalculateMLrate_Hurley2000(
        const double p_Mass,
        const double p_Radius,
        const double p_Luminosity,
        const double p_PerturbationMu = 0.0
    ) const;

    virtual COMPAS_PURE MASS_LOSS_T CalculateMLrate_Merritt2025(
        const double p_Mass,
        const double p_Radius,
        const double p_Luminosity,
        const double p_Temperature,
        const double p_PerturbationMu,
        const double p_mStart,
        const double p_HeAbundanceSurface
    ) const;

    COMPAS_PURE MASS_LOSS_T CalculateMLrate_NieuwenhuijzenDeJager1990(const double p_Mass, const double p_Radius, const double p_Luminosity) const;
    
    GNU_CONST   MASS_LOSS_T CalculateMLrate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const;

    COMPAS_PURE double BaseStar::CalculateMLrate_WRenhanced(
        const double                p_Luminosity, 
        const double                p_Temperature, 
        const double                p_HeAbundanceSurface, 
        const double                p_dMdtOther,
        const std::optional<double> p_dMdtWR
    ) const;

    // mass loss RSG
    COMPAS_PURE MASS_LOSS_T CalculateMLrateRSG(
        const double                     p_Mass,
        const double                     p_Radius,
        const double                     p_Luminosity,
        const double                     p_Temperature,
        const double                     p_mStart,
        const RSG_MASS_LOSS_PRESCRIPTION p_RSG_MLprescription
    ) const;
    GNU_CONST   MASS_LOSS_T CalculateMLrateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const;
    GNU_CONST   MASS_LOSS_T CalculateMLrateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const;
    GNU_CONST   MASS_LOSS_T CalculateMLrateRSG_Kee2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    GNU_CONST   MASS_LOSS_T CalculateMLrateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const;
    GNU_CONST   MASS_LOSS_T CalculateMLrateRSG_Yang2023(const double p_Luminosity) const;

    // mass loss THERMAL
    virtual double CalculateMLrateThermal() const { return Mass() / CalculateTimescale_Thermal(); }

    // mass loss VMS
    COMPAS_PURE MASS_LOSS_T CalculateMLrateVMS(
        const double                     p_Mass,
        const double                     p_Luminosity,
        const double                     p_Temperature,
        const VMS_MASS_LOSS_PRESCRIPTION p_VMS_MLprescription
    ) const;
    GNU_CONST   MASS_LOSS_T CalculateMLrateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const;
    COMPAS_PURE MASS_LOSS_T CalculateMLrateVMS_Sabhahit2023(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE MASS_LOSS_T CalculateMLrateVMS_Vink2011(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;

    // mass loss WR
    GNU_CONST   MASS_LOSS_T CalculateMLrateWR_Hurley2000(const double p_Luminosity, const double p_PerturbationMu) const;
    COMPAS_PURE static MASS_LOSS_T CalculateMLrateWR_SanderVink2020_Static(const double p_Luminosity, const double p_PerturbationMu) const;
    COMPAS_PURE static MASS_LOSS_T CalculateMLrateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE static MASS_LOSS_T CalculateMLrateWR_ZDependent_Static(const double p_Luminosity, const double p_PerturbationMu) const;

    
    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////

    virtual double CalculateRadius() const;

    COMPAS_PURE double CalculateRadiusAtZAMS_Tout1996(const double p_MZAMS) const;

    virtual inline double CalculateRadius_Hurley2000() const { return RZAMS(); } 

    virtual double CalculateRemnantRadius() const;
    virtual inline double CalculateRemnantRadius_Hurley2000() const { return Radius(); } // relevant for MS stars, over-written for GB stars


    ////////////////////////////////////////
    //   ROTATION                         //
    ////////////////////////////////////////

    COMPAS_PURE double BaseStar::CalculateRotationalVelocityOStar_Ramirez2013() const;


    ////////////////////////////////////////
    //   SUPERNOVAE                       //
    ////////////////////////////////////////

    COMPAS_PURE StellarSNDetailsT CalculateSNkickMagnitude(
        const STELLAR_TYPE       p_RemnantType,
        const double             p_Mass,
        const double             p_EjectaMass,
        const double             p_RemnantMass,
        const StellarSNDetailsT& p_SNdetails
    ) const;

    COMPAS_PURE double DrawSNkickMagnitude(
        const double   p_COCoreMass,
        const double   p_EjectaMass,
        const double   p_RemnantMass,
        const double   p_Sigma,
        const double   p_Rand,
        const SN_EVENT p_SNevent
    ) const;

    GNU_CONST   double DrawSNkickMagnitude_Muller2016(const double p_COCoreMass) const;
    COMPAS_PURE double DrawSNkickMagnitude_MullerMandel2020(const double p_COCoreMass, const double p_RemnantMass, const double p_Rand, const SN_EVENT p_SNevent) const;

    GNU_CONST   double DrawSNkickMagnitude_BrayEldridge2018(const double p_EjectaMass, const double p_RemnantMass, const double p_Alpha, const double p_Beta) const;
    GNU_CONST   double DrawSNkickMagnitude_Flat(const double p_KickMax, const double p_Rand) const;
    GNU_CONST   double DrawSNkickMagnitude_Maxwellian(const double p_Sigma, const double p_Rand) const;
    GNU_CONST   double DrawSNkickMagnitude_Muller2016(const double p_COCoreMass) const;
    COMPAS_PURE double DrawSNkickMagnitude_MullerMandel2020(const double p_COCoreMass, const double p_RemnantMass, const double p_Rand, const SN_EVENT p_SNevent) const;




























virtual inline double CalculateHurleyPerturbationMu() const { return m_Mu; }
virtual inline double CalculateHurleyPerturbationMuAtPhaseEnd() const { return CalculateHurleyPerturbationMu(); }

GNU_CONST inline double CalculateHurleyPerturbationB(const double p_Mass) const { return 0.002 * std::max(1.0, (2.5 / p_Mass)); } // Hurley et al. 2000, eq 103

GNU_CONST inline double CalculateHurleyPerturbationC(double p_Mass) const { return 0.006 * std::max(1.0, (2.5 / p_Mass)); } // Hurley et al. 2000, eq 104

GNU_CONST inline double CalculatePerturbationQ(const double p_Radius, const double p_Rc) const { return log(p_Radius / p_Rc); } // Hurley et al. 2000, eq 105 (really is natural log)

GNU_CONST double CalculateHurleyPerturbationR(const double p_Mass, const double p_Radius, const double p_Mu, const double p_Rc) const;

GNU_CONST double CalculateHurleyPerturbationS(const double p_Mass, const double p_Mu) const;



GNU_CONST inline double BaseStar::CalculateTimescale_RadialExpansion(
    const STELLAR_TYPE p_StellarType,
    const STELLAR_TYPE p_StellarTypePrev,
    const double       p_Radius,
    const double       p_RadiusPrev,
    const double       p_dtPrev
) const;








    virtual std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                              { return std::make_tuple(CalculateRadiusOnPhase(), m_StellarType); }

    void                CalculateRCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_RCoefficients) const;


GNU_PURE double CalculateAngularFrequencyAtZAMS_Hurley2000(const double p_MZAMS, const double p_RZAMS);
GNU_PURE double CalculateRotationalVelocityAtZAMS(double p_MZAMS);





///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

virtual DBL_VECTOR CalculateGBparams() const;
virtual DBL_VECTOR CalculateGBparams_Hurley2000() const { return GBparams(); } // only required for GB stars






























///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

virtual double CalculateCoreMassAtPhaseEnd() const;
virtual inline double CalculateCoreMassAtPhaseEnd_Hurley2000() const { return CoreMass(); }


virtual double CalculateCOCoreMassAtPhaseEnd() const;
virtual inline double CalculateCOCoreMassAtPhaseEnd_Hurley2000() const { COCoreMass(); }

virtual double CalculateHeCoreMassAtPhaseEnd() const;
virtual inline double CalculateHeCoreMassAtPhaseEnd_Hurley2000() const { HeCoreMass(); }

virtual double CalculateRadiusAtPhaseEnd() const;
virtual inline double CalculateRadiusAtPhaseEnd_Hurley2000() const { Radius(); }







    virtual double      CalculateTemperatureAtPhaseEnd() const                                                          { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }
    virtual double      CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const          { return CalculateTemperatureOnPhase(p_Luminosity, p_Radius); }             // Same as on phase
    double              CalculateTemperatureKelvinOnPhase(const double p_Luminosity, const double p_Radius) const;


////    virtual double      CalculateTemperatureOnPhase() const                                                             { return CalculateTemperatureOnPhase(m_Luminosity, m_Radius); }



    virtual double      CalculateTemperatureOnPhase(const double p_Luminosity, const double p_Radius) const;






    virtual double      ChooseTimestep(const double p_Time) const                                                       { return m_dt; }


    STELLAR_TYPE        EvolveOnPhase(const double p_DeltaTime);

    virtual STELLAR_TYPE EvolveToNextPhase()                                                                            { return m_StellarType; }




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

    
    virtual void        SetSNHydrogenContent()                                                                          { m_SupernovaDetails.isHydrogenPoor = false; }                              // Default is false

    bool         ShouldBeMasslessRemnant() const                                                                 { return (m_Mass <= 0.0 || m_StellarType == STELLAR_TYPE::MASSLESS_REMNANT); }
    virtual bool        ShouldEvolveOnPhase() const                                                                     { return true; }
    virtual bool        ShouldSkipPhase() const                                                                         { return false; }                                                           // Default is false



    /*
     * OStarRotationVelocityFunctor
     *
     * @brief
     * Functor for CalculateRotationalVelocityOStar_Ramirez2013()
     *
     * Find the ifference between star's Roche Lobe radius and radius after mass loss.
     * 
     * The analytic cumulative distribution function (CDF) for the equatorial rotational velocity
     * of single O stars is modelled as a mixture of a gamma component and a normal component, and
     * computed per Ramirez-Agudelo et al. 2013, eqs 1-4 (see https://arxiv.org/abs/1309.2929).
     *
     * 
     * Constructor: initialise the class
     * template <class T> OStarRotationVelocityFunctor(double p_CDF, ERROR *p_Error)
     *
     * @param       p_CDF                           Desired CDF value
     *
     * 
     * Function: calculate the CDF of the O star rotational velocity and compare to desired value
     * T OStarRotationVelocityFunctor(double const& p_Ve)
     *
     * @param       p_Ve                            Rotational velocity, km s^-1
     * @return                                      Difference between star's Roche Lobe radius and radius after mass loss
     */
    template <class T>
    struct OStarRotationVelocityFunctor {
        OStarRotationVelocityFunctor(double p_CDF) { m_CDF = p_CDF; }
        T operator()(double const& p_Ve) { 

            constexpr double alpha  = 4.82;
            constexpr double beta   = 1.0 / 25.0;
            constexpr double mu     = 205.0;
            constexpr double sigma  = 190.0;
            constexpr double iGamma = 0.43;

            boost::math::inverse_gamma_distribution<> gammaComponent(alpha, beta); // (shape, scale) = (alpha, beta)
            boost::math::normal_distribution<> normalComponent(mu, sigma);
    
            // Compute CDF at zero rotational velocity -- the CDF should be relative to this quantity
            const double CDFzero         = (iGamma * boost::math::cdf(gammaComponent, 0.0)) + ((1.0 - iGamma) * boost::math::cdf(normalComponent, 0.0));
            const double CDFunnormalised = (iGamma * boost::math::cdf(gammaComponent, p_Ve)) + ((1.0 - iGamma) * boost::math::cdf(normalComponent, p_Ve));
    
            return ((CDFunnormalised - CDFzero) / (1.0 - CDFzero)) - m_CDF;
        }
    private:
        double m_CDF;
    };


};






///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         AGE / LIFETIME / TAU / TIMESCALES                         //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateAgeAfterMassLoss
 *
 * @brief
 * Recalculate the age of a star after mass loss.
 * 
 * Calls relevant age function based on the evolutionary mode given in program options.
 *
 *
 * double CalculateAgeAfterMassLoss() const
 *
 * @return                                      Age of the star after mass loss (Myr)
 */
inline double BaseStar::CalculateAgeAfterMassLoss() const {

    double age;

    Switch (OPTIONS->Mode()) {                                                                                  // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                                             // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                                             // HURLEY BSE
            age = CalculateAgeAfterMassLoss_Hurley2000();
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
 * @param       p_tHeI                          Time to helium ignition, tHeI (per Hurley timescales) (Myr)
 * @param       p_tHe                           Time to helium burning, tHe (per Hurley timescales) (Myr)
 * @return                                      Lifetime to the BAGB (Myr)
 */
GNU_CONST inline double BaseStar::CalculateLifetimeToBAGB_Hurley2000(const double p_tHeI, const double p_tHe) const {
    return p_tHeI + p_tHe;
}


/*
 * CalculateLifetimeToBGB_Hurley2000_Static
 *
 * @brief
 * Calculate the lifetime to the Base of the Giant Branch, BGB, (end of the Hertzsprung Gap),
 * per Hurley at al. 2000, eq 4 (plotted in Hurley et al. 2000, fig 5)
 * 
 * For high mass stars, tBGB = tHeI.
 *
 *
 * static double CalculateLifetimeToBGB_Hurley2000_Static(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Lifetime to the Base of the Giant Branch, tBGB (Myr)
 */
COMPAS_PURE inline double BaseStar::CalculateLifetimeToBGB_Hurley2000_Static(const double p_Mass) const {

    const DBL_VECTOR a = GLOBALS->HurleyAcoefficients(); // get Hurley a coefficients

    const double m2 = p_Mass * p_Mass;
    const double m4 = m2 * m2;
    const double m7 = p_Mass * m2 * m4;

    return (a[1] + (a[2] * m4) + (a[3] * std::sqrt(p_Mass) * p_Mass * m4) + m7) / ((a[4] * m2) + (a[5] * m7));
}


/*
 * CalculateTau
 *
 * @brief
 * Calculate the phase-relative age of the star (fractional age on the current
 * evolutionary phase).
 *
 * Calls relevant tau function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateTau() const
 *
 * @return                                      MS-relative age, [0, 1]
 */
inline double BaseStar::CalculateTau() const { 

    double tau;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            tau = CalculateTau_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return tau;
}


/*
 * CalculateTimescales
 *
 * @brief
 * (Re)calculate timescales.
 * 
 * Timescales should be recalculated every timestep.
 *
 * Calls relevant timescales function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateTimescales() const
 *
 * @return                                      Timescales (Myr)
 */
inline DBL_VECTOR BaseStar::CalculateTimescales() const { 

    DBL_VECTOR tScales;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            tScales = CalculateTimescales_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return tScales;
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
GNU_CONST inline double BaseStar::CalculateTimescale_Dynamical_Kalogera1996(const double p_Mass, const double p_Radius) const {
    return 5.0 * 1.0E-5 * p_Radius * std::sqrt(p_Radius) * YEAR_TO_MYR / std::sqrt(p_Mass);
}


/*
 * CalculateTimescale_Thermal
 *
 * @brief
 * Calculate thermal timescale
 *
 * Pre-factor from Kalogera & Webbink 1996 (https://arxiv.org/abs/astro-ph/9508072), equation 2, 
 * combined with p_Mass * p_EnvMass case from equation 61 from https://arxiv.org/abs/astro-ph/0201220 for k in {2,3,4,5,6,8,9}
 * [note that equation 61 of BSE (https://arxiv.org/abs/astro-ph/0201220) approximates this with a value a factor of 3 smaller]
 * 
 * 
 * double CalculateTimescale_Thermal(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Thermal timescale (Myr)
 */
GNU_CONST inline double BaseStar::CalculateTimescale_Thermal(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const {   
    // G*Msol^2/(Lsol*Rsol) ~ 31.4 Myr (~ 30 Myr in Kalogera & Webbink)
    return 31.4 * p_Mass * (p_Mass == p_CoreMass ? p_Mass : p_Mass - p_CoreMass) / (p_Radius * m_Luminosity);
}


/*
 * CalculateTimescale_RadialExpansion
 *
 * @brief
 * Calculate radial expansion timescale
 *
 * If the radial expansion timescale can't be calculated (because, for example,
 * the stellar type has changed, or the radius has not changed), the radial
 * expansion timescale is returned as -1.0.
 * 
 *
 * double CalculateTimescale_RadialExpansion(
 *     const STELLAR_TYPE p_StellarType,
 *     const STELLAR_TYPE p_StellarTypePrev,
 *     const double       p_Radius,
 *     const double       p_RadiusPrev,
 *     const double       p_dtPrev
 * ) const
 *
 * @param       p_StellarType                   Current stellar type of the star
 * @param       p_StellarTypePrev               Previous stellar type of the star
 * @param       p_Radius                        Current radius of the star (Rsol)
 * @param       p_RadiusPrev                    Previous radius of the star (Rsol)
 * @param       p_dtPrev                        Previous timestep (Myr)
 * @return                                      Radial expansion timescale (Myr)
 *                                              Note: returns -1.0 if radial expansion timescale can't be calculated
 *                                              (because, for example, the stellar type has changed, or the radius has not changed)
 */
GNU_CONST inline double BaseStar::CalculateTimescale_RadialExpansion(
    const STELLAR_TYPE p_StellarType,
    const STELLAR_TYPE p_StellarTypePrev,
    const double       p_Radius,
    const double       p_RadiusPrev,
    const double       p_dtPrev
) const {
    return (p_StellarTypePrev == p_StellarType && p_RadiusPrev != p_Radius) ? (p_dtPrev * p_RadiusPrev) / std::fabs(p_Radius - p_RadiusPrev) : -1.0;
}


/*
 * CalculateTimescale_RadialExpansion_DuringMT
 *
 * @brief       // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 * Calculate the radial expansion timescale in the mass transfer regime
 * We do not use CalculateRadialExpansionTimescale(), since in the process of mass transfer the previous radius
 * is determined by binary evolution, not nuclear timescale evolution
 *
 *
 * double CalculateTimescale_RadialExpansion_DuringMT()
 *
 * @return                                      Radial expansion timescale
 */
double BaseStar::CalculateTimescale_RCalculateTimescale_RadialExpansion_DuringMTadialExpansionDuringMT() {
    
    // We create a clone of this star, with the same stellar type, and age it slightly to
    // determine how the radius will change.
    // To be sure the clone does not participate in logging, we set its persistence to EPHEMERAL.
    std::unique_ptr<BaseStar> clone = CloneAs(StellarType(), OBJECT_PERSISTENCE::EPHEMERAL, false); // do not re-initialise the clone

    const double timestep = std::max(1000.0 * NUCLEAR_MINIMUM_TIMESTEP, Age() / 1.0E6);   // timestep to age clone (Myr)
    (void)clone->EvolveOneTimestep(0.0, 0.0, timestep, true); // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    return timestep * Radius() / std::fabs(Radius() - clone->Radius());
}


/*
 * CalculateTimescale_MassChange
 *
 * @brief
 * Calculate mass change timescale
 *
 * If the mass change timescale can't be calculated (because, for example,
 * the stellar type has changed, or the mass has not changed), the mass
 * change timescale is returned as -1.0.
 * 
 * 
 * double CalculateTimescale_MassChange(
 *     const STELLAR_TYPE p_StellarType,
 *     const STELLAR_TYPE p_StellarTypePrev,
 *     const double       p_Mass,
 *     const double       p_MassPrev,
 *     const double       p_dtPrev
 * ) const
 *
 * @param       p_StellarType                   Current stellar type of the star
 * @param       p_StellarTypePrev               Previous stellar type of the star
 * @param       p_Mass                          Current mass of the star (Msol)
 * @param       p_MassPrev                      Previous mass of the star (Msol)
 * @param       p_dtPrev                        Previous timestep (Myr)
 * @return                                      Mass change timescale (Myr)
 *                                              Note: returns -1.0 if mass change timescale can't be calculated
 *                                              (because, for example, the stellar type has changed, or the mass has not changed)
 */
GNU_CONST inline double BaseStar::CalculateTimescale_MassChange(
    const STELLAR_TYPE p_StellarType,
    const STELLAR_TYPE p_StellarTypePrev,
    const double       p_Mass,
    const double       p_MassPrev,
    const double       p_dtPrev
) const {
    return p_StellarTypePrev == p_StellarType && p_MassPrev != p_Mass ? (p_dtPrev * p_MassPrev) / std::fabs(p_Mass - p_MassPrev) : -1.0;
}



/*
 * Calculate the eddy turnover timescale
 * Hurley et al. 2002, sec 2.3, particularly eq 31 of subsec 2.3.1
 *
 *
 * double CalculateEddyTurnoverTimescale()
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Mass change timescale (Myr)
 * @return                                      eddy turnover timescale (yr)   /// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 */
double BaseStar::CalculateTimescale_EddyTurnover_Hurley2002(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const {

	const double rEnv = CalculateConvectiveEnvelopeRadialExtent();
    double mEnv;
    std::tie(mEnv, std::ignore) = CalculateConvectiveEnvelopeMass(p_Mass, p_CoreMass);
    return 0.4311 * std::cbrt((mEnv * rEnv * (p_Radius - (0.5 * rEnv))) / (3.0 * p_Luminosity));
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           MISC. CONSTANTS / PARAMETERS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateGBparams
 *
 * @brief
 * (Re)calculate the Giant Branch parameters.
 * 
 * Giant Branch parameters should be recalculated every timestep.
 *
 * Calls relevant GBparams function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateGBparams() const
 *
 * @return                                      Timescales (Myr)
 */
GNU_CONST inline DBL_VECTOR BaseStar::CalculateGBparams() const { 

    DBL_VECTOR GBparams;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            GBparams = CalculateGBparams_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return GBparams;
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
 * Calls relevant radius function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateRadius() const
 *
 * @return                                      Radius (Rsol)
 */
inline double BaseStar::CalculateRadius() const { 

    double radius;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            radius = CalculateRadius_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
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
 * Calls relevant radius function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateRemnantRadius() const
 *
 * @return                                      Radius (Rsol)
 */
inline double BaseStar::CalculateRemnantRadius() const { 

    double radius;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            radius = CalculateRemnantRadius_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
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
GNU_CONST inline double CalculateOmegaBreak(const double p_Mass, const double p_Radius) const {
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
GNU_CONST inline double BaseStar::DrawSNkickMagnitude_BrayEldridge2018(const double p_EjectaMass, const double p_RemnantMass, const double p_Alpha, const double p_Beta) const {
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
GNU_CONST inline double BaseStar::DrawSNkickMagnitude_Flat(const double p_KickMax, const double p_Rand) const {
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
GNU_CONST inline double BaseStar::DrawSNkickMagnitude_Maxwellian(const double p_Sigma, const double p_Rand) const {
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
double BaseStar::CalculateHurleyPerturbationS(const double p_Mass, const double p_Mu) const {

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
GNU_CONST inline double BaseStar::CalculateHurleyPerturbationR(const double p_Mass, const double p_Radius, const double p_Mu, const double p_Rc) const {

    double r = 0.0;

    // we only perturb the radius if mu > 0 and radius is larger than core radius,
    // otherwise r = 0 and perturbed radius = core radius
    if (p_Mu > 0.0 && p_Radius > p_Rc) {

        const double c      = CalculateHurleyPerturbationC(p_Mass);
        const double c3     = c * c * c;                    // pow() is slow - use multiplication
        const double mu3_c3 = p_Mu * p_Mu * p_Mu / c3;      // calculate once

        // Hurley et al. 2000 is just 0.1 / q here, but the Hurley sse code does this (`rpertf()` in `zfuncs.f`) - no explanation.
        const double exponent = min((0.1 / CalculatePerturbationQ(p_Radius, p_Rc)), (-14.0 / log10(p_Mu)));

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
 * double CalculateLuminosity_Hurley2000(const double p_CoreMass, const DBL_VECTOR& p_GBparams) const
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      AGB luminosity (Lsol)
 */
GNU_CONST inline double BaseStar::CalculateLuminosity_Hurley2000(const double p_CoreMass, const DBL_VECTOR& p_GBparams) const {
#define GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)] // for convenience and readability - undefined at end of function
    return std::min((GBparams(B) * PPOW(p_CoreMass, GBparams(q))), (GBparams(D) * PPOW(p_CoreMass, GBparams(p))));
#undef GBparams
}


/*
 * CalculateLuminosityAtBAGB_Hurley2000
 *
 * @brief
 * Calculate luminosity at the base of the Asymptotic Giant Branch,
 * per Hurley et al. 2000, eq 56
 *
 *
 * double CalculateLuminosityAtBAGB_Hurley2000(double p_Mass) const
 *
 * @param       p_Mass                          Effective initial mass (Msol)
 * @return                                      BAGB luminosity (Lsol)
 */
COMPAS_PURE inline double BaseStar::CalculateLuminosityAtBAGB_Hurley2000(double p_Mass) const {
#define b(x) GLOBALS->HurleyBCoefficients(x) // for convenience and readability - undefined at end of function

    return p_Mass < p_MHeF
            ? (b(29) * PPOW(p_Mass, b(30))) / (1.0 + (GLOBALS->HurleyAlpha3() * exp(15.0 * (p_Mass - GLOBALS->HurleyMassCutoffs(static_cast<int>(MASS_CUTOFF::MHeF))))))
            : (b(31) + (b(32) * PPOW(p_Mass, (b(33) + 1.8)))) / (b(34) + PPOW(p_Mass, b(33)));

#undef bN
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               TEMPERATURE FUNCTIONS                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateTemperatureOnPhase
 *
 * @brief
 * Calculate the temperature of the star at the current evolutionary phase.
 *
 * Calls relevant temperature function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateTemperatureOnPhase()
 *
 * @return                                      Temperature of the star (Tsol)
 */
inline double BaseStar::CalculateTemperatureOnPhase() const { 

    double temperature;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            // retrieve luminosity and radius from current state and call relevant function
            temperature = CalculateTemperatureOnPhase(Luminosity(), Radius());
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return temperature;
}


/*
 * CalculateTemperatureOnPhase
 *
 * @brief
 * Calculate the effective temperature of the star at the current evolutionary phase,
 * given the luminosity and radius of the star, using the Stefan–Boltzmann law.
 *
 *
 * double CalculateTemperatureOnPhase(const double p_Luminosity, const double p_Radius) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @return                                      Effective temperature of the star (Tsol)
 */
GNU_CONST inline double BaseStar::CalculateTemperatureOnPhase(const double p_Luminosity, const double p_Radius) const {
    return std::sqrt(std::sqrt(p_Luminosity)) / std::sqrt(p_Radius);
}



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
 * Calls relevant core mass function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateCoreMass() const
 *
 * @return                                      Core mass (Msol)
 */
inline double BaseStar::CalculateCoreMass() const { 

    double coreMass;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            coreMass = CalculateCoreMass_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       


    switch (OPTIONS->MainSequenceCoreMassPrescription()) {              // which MS core mass prescription?
        
        case MS_CORE_MASS_PRESCRIPTION::BRCEK:                          // BRCEK


    return coreMass;
}


/*
 * CalculateCoreMassAtPhaseEnd
 *
 * @brief
 * Calculate the core mass of the star at the end of the current evolutionary phase.
 *
 * Calls relevant core mass function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateCoreMassAtPhaseEnd() const
 *
 * @return                                      Core mass (Msol)
 */
inline double BaseStar::CalculateCoreMassAtPhaseEnd() const { 

    double coreMass;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            coreMass = CalculateCoreMassAtPhaseEnd_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return coreMass;
}


/*
 * CalculateCoreMass_Hurley2000_Static
 *
 * @brief
 * Calculate the core mass of the star, per Hurley et al. 2000, eqs 37 & 38.
 * 
 * The Hurley et al. 2000 Mc-L relation is described in section 5.2 of the paper (just before eq 37).
 * 
 * 
 * static double CalculateCoreMass_Hurley2000_Static(const double p_Luminosity, const DBL_VECTOR& p_GBparams)
 * 
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Core mass (Msol)
 */
GNU_CONST static inline double CalculateCoreMass_Hurley2000_Static(const double p_Luminosity, const DBL_VECTOR& p_GBparams) {
    return p_Luminosity > p_GBparams[static_cast<int>(HURLEY_GBP:::Lx)]
            ? PPOW((p_Luminosity / p_GBparams[static_cast<int>(HURLEY_GBP:::B)]), (1.0 / p_GBparams[static_cast<int>(HURLEY_GBP:::q)]))
            : PPOW((p_Luminosity / p_GBparams[static_cast<int>(HURLEY_GBP:::D)]), (1.0 / p_GBparams[static_cast<int>(HURLEY_GBP:::p)]));
}


/*
 * CalculateCOCoreMass
 *
 * @brief
 * Calculate the carbon-oxygen (CO) core mass of the star.
 *
 * Calls relevant CO core mass function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateCOCoreMass() const
 *
 * @return                                      CO core mass (Msol)
 */
inline double BaseStar::CalculateCOCoreMass() const { 

    double COcoreMass;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            COcoreMass = CalculateCOCoreMass_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return COcoreMass;
}


/*
 * CalculateCOCoreMassAtPhaseEnd
 *
 * @brief
 * Calculate the carbon-oxygen (CO) core mass of the star at the end of the current evolutionary
 * phase.
 *
 * Calls relevant CO core mass function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateCOCoreMassAtPhaseEnd() const
 *
 * @return                                      CO core mass (Msol)
 */
inline double BaseStar::CalculateCOCoreMassAtPhaseEnd() const { 

    double COcoreMass;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            COcoreMass = CalculateCOCoreMassAtPhaseEnd_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return COcoreMass;
}


/*
 * CalculateHeCoreMass
 *
 * @brief
 * Calculate the Helium (He) core mass of the star.
 *
 * Calls relevant He core mass function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateHeCoreMass() const
 *
 * @return                                      He core mass (Msol)
 */
inline double BaseStar::CalculateHeCoreMass() const { 

    double HecoreMass;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            HecoreMass = CalculateHeCoreMass_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return HecoreMass;
}


/*
 * CalculateHeCoreMassAtPhaseEnd
 *
 * @brief
 * Calculate the Helium (He) core mass of the star at the end of the current evolutionary
 * phase.
 *
 * Calls relevant He core mass function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateHeCoreMassAtPhaseEnd() const
 *
 * @return                                      He core mass (Msol)
 */
inline double BaseStar::CalculateHeCoreMassAtPhaseEnd() const { 

    double HecoreMass;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            HecoreMass = CalculateHeCoreMassAtPhaseEnd_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return HecoreMass;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     MASS LOSS                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
 */
GNU_CONST double inline BaseStar::CalculateMLfractionWR(const double p_HeAbundanceSurface) const {

    constexpr double limOB = 0.55;                                          // per Yoon et al. 2006
    constexpr double limWR = 0.70;                                          // per Yoon et al. 2006

    return std::min(1.0, std::max (0.0, (p_HeAbundanceSurface - limOB) / (limWR - limOB)));
}


/*
 * CalculateMLrate_Hurley2000
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, per Hurley et al. 2000.
 * 
 * 
 * MASS_LOSS_T CalculateMLrate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_PerturbationMu) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu (optional parameter)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (may be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE inline MASS_LOSS_T BaseStar::CalculateMLrate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_PerturbationMu) const {
    return CalculateMLrate_NieuwenhuijzenDeJager1990(p_Mass, p_Radius, p_Luminosity);
}


/*
 * CalculateMLrate_KudritzkiReimers1978
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, on the GB and beyond,
 * per Hurley et al. 2000, eq 106
 * 
 * Based on a prescription taken from Kudritzki and Reimers 1978
 *
 *
 * MASS_LOSS_T CalculateMLrate_KudritzkiReimers1978(const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::GB)
 * @return                                      Kudritzki and Reimers mass loss rate (in Msol yr^{-1})
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLrate_KudritzkiReimers1978(const double p_Mass, const double p_Radius, const double p_Luminosity) const {
    // Hurley et al. 2000 has eta^2 - that's wrong per Hurley SSE code
    return std::make_tuple(4.0E-13 * (MASS_LOSS_ETA * p_Luminosity * p_Radius / p_Mass), MASS_LOSS_TYPE::GB);
}


/* 
 * CalculateMLrate_VassiliadisWood1993
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for stars on the AGB,
 * based on the Mira pulsation period (P0), per Vassiliadis and Wood 1993.
 *
 * See Hurley et al. 2000, just after eq 106
 * Note: in the Hurley fortran code, P0 is taken to be min(p0, 2000.0) - implemented here as a minimum power
 *
 *
 * MASS_LOSS_T CalculateMLrate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate on the AGB (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::GB)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLrate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const {
    const double P0   = PPOW(10.0, std::min(3.3, (-2.07 - (0.9 * log10(p_Mass)) + (1.94 * log10(p_Radius)))));
    const double dMdt = PPOW(10.0, (-11.4 + (0.0125 * (P0 - 100.0 * std::max((p_Mass - 2.5), 0.0)))));

    return std::make_tuple(std::min(dMdt, (1.36E-9 * p_Luminosity)), MASS_LOSS_TYPE::GB);
}


/*
 * CalculateMLrateLBV_Belczynski2010
 *
 * @brief
 * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit
 * (Humphreys & Davidson 1994), per Belczynski et al. 2010, eq 8 
 * 
 *
 * MASS_LOSS_T CalculateMLrateLBV_Belczynski2010() const
 *
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for LBV stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::LBV)
 */    
COMPAS_PURE inline MASS_LOSS_T BaseStar::CalculateMLrateLBV_Belczynski2010() const {
    return std::make_tuple(OPTIONS->LuminousBlueVariableFactor() * 1.0E-4, MASS_LOSS_TYPE::LBV);
} 


/*
 * CalculateMLrateLBV_Hurley2000
 *
 * @brief
 * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit
 * (Humphreys & Davidson 1994), per Hurley et al. 2000, sec 7.1, unlabelled equation
 * a few equations after eq 106
 *  
 *
 * MASS_LOSS_T CalculateMLrateLBV_Hurley2000(const double p_Luminosity, const double p_HDlimitfactor) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_HDlimitfactor                 Factor by which star is above Humphreys-Davidson limit
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for LBV stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::LBV)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLrateLBV_Hurley2000(const double p_Luminosity, const double p_HDlimitfactor) const {
    return std::make_tuple(0.1 * utils::IntPow(p_HDlimitfactor - 1.0, 3) * ((p_Luminosity / 6.0E5) - 1.0), MASS_LOSS_TYPE::LBV);
}


/*
 * CalculateMLrateOB_Krticka2018
 *
 * @brief
 * Calculate mass loss rate for massive OB stars, per Krticka et al. 2018
 * (See https://arxiv.org/pdf/1712.03321.pdf)
 *
 * 
 * MASS_LOSS_T CalculateMLrateOB_Krticka2018(const double p_Luminosity) const
 * 
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for hot OB stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::OB)
 */
COMPAS_PURE inline MASS_LOSS_T BaseStar::CalculateMLrateOB_Krticka2018(const double p_Luminosity) const {
    return std::make_tuple(PPOW(10.0, -5.70 + 0.50 * GLOBALS->ZetaAsplund()() + (1.61 - 0.12 * GLOBALS->ZetaAsplund()()) * log10(p_Luminosity / 1.0E6)), MASS_LOSS_TYPE::OB);
}


/*
 * CalculateMLrateRSG_Beasor2020
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
 * MASS_LOSS_T CalculateMLrateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const
 *
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLrateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const {
    return std::make_tuple(PPOW(10.0, -21.5 - 0.15 * p_mStart + 3.6 * log10(p_Luminosity)), MASS_LOSS_TYPE::RSG);
}


/*
 * CalculateMLrateRSG_Decin2023
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for RSG stars (Red Supergiants),
 * per Decin 2023, eq 6.
 *
 * See https://arxiv.org/pdf/2303.09385.pdf
 * 
 *  
 * MASS_LOSS_T CalculateMLrateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const
 *
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLrateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const {
    return std::make_tuple(PPOW(10.0, -20.63 - 0.16 * p_mStart + 3.47 * log10(p_Luminosity)), MASS_LOSS_TYPE::RSG);
}


/*
 * CalculateMLrateRSG_VinkSabhahit2023
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
 * MASS_LOSS_T CalculateMLrateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLrateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const {
    const double exp = logL < 4.6 ? -8.0 + 0.7 * log10(p_Luminosity) - 0.7 * log10(p_Mass) : -24.0 + 4.77 * log10(p_Luminosity) - 3.99 * log10(p_Mass);
    return std::make_tuple(PPOW(10.0, exp), MASS_LOSS_TYPE::RSG);
}


/*
 * CalculateMLrateRSG_Yang2023
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for RSG stars (Red Supergiants),
 * per Yang 2023, eq 6.
 * 
 * See https://arxiv.org/pdf/2303.09385.pdf
 *
 *  
 * GNU_CONST MASS_LOSS_T CalculateMLrateRSG_Yang2023(const double p_Luminosity) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLrateRSG_Yang2023(const double p_Luminosity) const {

    const double logL  = log10(p_Luminosity);
    const double logL2 = logL * logL;

    return std::make_tuple(PPOW(10.0, 0.45 * logL2 * logL - 5.26 * logL2 + 20.93 * logL - 34.56), MASS_LOSS_TYPE::RSG);
}


/*
 * CalculateMLrateVMS_Bestenlehner2020
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for very massive (> 100 Msol) OB stars,
 * per Bestenlehner 2020
 * 
 * See https://arxiv.org/pdf/2002.05168.pdf
 * 
 * 
 * MASS_LOSS_T CalculateMLrateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for very massive stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::VMS)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLrateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const {

    constexpr double alpha       = 0.39;                                        // CAK force multiplier
    constexpr double logMdotZero = -4.78;                                       // from substituting LogMdotTrans and Gamma_e trans into eq 12. 

    const double gamma = EDDINGTON_PARAMETER_FACTOR * p_Luminosity / p_Mass;    // Eddington Parameter, not metallicity specific as in the publication

    return std::make_tuple(PPOW(10.0, logMdotZero + ((1.0 / alpha) + 0.5) * log10(gamma) - (((1.0 - alpha) / alpha) + 2.0) * log10(1.0 - gamma)), MASS_LOSS_TYPE::VMS);
}


/*
 * CalculateMLrateWR_Hurley2000
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
 * MASS_LOSS_T CalculateMLrateWR_Hurley2000(const double p_Luminosity, const double p_PerturbationMu) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLrateWR_Hurley2000(const double p_Luminosity, const double p_PerturbationMu) const {
    return std::make_tuple((p_PerturbationMu >= 1.0 ? 0.0 : PPOW(p_Luminosity, 1.5) * (1.0 - p_PerturbationMu) * 1.0E-13), MASS_LOSS_TYPE::WR);
} 


/*
 * CalculateMLrateWR_Shenar2019_Static
 *
 * Calculate mass loss rate, and the dominant mass loss type, forWolf-Rayet stars,
 * per Shenar et al. 2019, eq 6, tbl 5.
 * 
 * See (https://ui.adsabs.harvard.edu/abs/2019A%26A...627A.151S/abstract)
 * 
 * We use the fitting coefficients for hydrogen rich WR stars (e.g., WNh).
 * The C4 (X_He) term is = 0 and is omitted.
 *  
 * 
 * static MASS_LOSS_T CalculateMLrateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
COMPAS_PURE inline MASS_LOSS_T BaseStar::CalculateMLrateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature) const {

    // For H-rich WR stars (X_H > 0.4)
    constexpr double C1 = -6.78;
    constexpr double C2 =  0.66;
    constexpr double C3 = -0.12;
    constexpr double C5 =  0.74;

    return std::make_tuple(PPOW(10.0, C1 + (C2 * log10(p_Luminosity)) + (C3 * log10(p_Temperature * TSOL)) + (C5 * GLOBALS->SigmaHurley())), MASS_LOSS_TYPE::WR);
}


/*
 * CalculateMLrateWR_ZDependent_Static
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
 * static MASS_LOSS_T CalculateMLrateWR_ZDependent_Static(const double p_Luminosity, const double p_PerturbationMu) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
COMPAS_PURE inline MASS_LOSS_T BaseStar::CalculateMLrateWR_ZDependent_Static(const double p_Luminosity, const double p_PerturbationMu) const {
    // StarTrack may still do something different here.
    // There are references to Hamann & Koesterke 1998 and Vink and de Koter 2005.
    // TW - Haven't seen StarTrack but I think H&K gives the original equation and V&dK gives the Z dependence
    const double dMdt = p_PerturbationMu >= 1.0 ? 0.0 : 1.0E-13 * p_Luminosity * std::sqrt(p_Luminosity) * PPOW(GLOBALS->Metallicity() / ZSOL_ANDERS, 0.86) * (1.0 - p_PerturbationMu);
    return std::make_tuple(dMdt, MASS_LOSS_TYPE::WR);
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
GNU_CONST inline double BaseStar::CalculateEddingtonLuminosity(const double p_Mass, const double p_HeAbundanceSurface) const {
    const double HeAbundance2 = p_HeAbundanceSurface + p_HeAbundanceSurface;   
    return (4.0 * M_PI * G * C * p_Mass * MSOL_TO_KG) / ((0.4 * (1.0 + 2.0 * HeAbundance2) / (1.0 + 4.0 * (HeAbundance2 + HeAbundance2))) * OPACITY_CGS_TO_S);
}


/*
 * CalculateLuminosity
 *
 * @brief
 * Calculate the luminosity of the star.
 *
 * Calls relevant luminosity function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateLuminosity() const
 *
 * @return                                      Luminosity (Lsol)
 */
inline double BaseStar::CalculateLuminosity() const { 

    double luminosity;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            luminosity = CalculateLuminosity_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return luminosity;
}


/*
 * CalculateLuminosityAtPhaseEnd
 *
 * @brief
 * Calculate the luminosity of the star at the end of the current phase.
 *
 * Calls relevant luminosity function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateLuminosityAtPhaseEnd() const
 *
 * @return                                      Luminosity (Lsol)
 */
inline double BaseStar::CalculateLuminosityAtPhaseEnd() const { 

    double luminosity;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            luminosity = CalculateLuminosityAtPhaseEnd_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return luminosity;
}












/*
 * IsOneOf
 *
 * @brief
 * Determines if the current stellar type of the star is one of the stellar types
 * in the list of stellar types passed as p_List.
 * 
 * 
 * bool IsOneOf(const STELLAR_TYPE_LIST p_List) const
 * 
 * @param
 * @param       p_List                          List of stellar types (std::initializer_list<STELLAR_TYPE>)
 * @return                                      True if the star's stellar type is one of the stellar types in p_List, otherwise False
 */
inline bool BaseStar::IsOneOf(const STELLAR_TYPE_LIST p_List) const {
    for (auto elem: p_List) {
        if (StellarType() == elem) return true;
    }
    return false;
};












//// <<<<<<<<<<<<<<<<<<<<<<<<<<<<< constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    virtual double      CalculateMTRejuvenationFactor()                                                       { return 1.0; }

    virtual MT_CASE             DetermineMassTransferTypeAsDonor() const                                        { return MT_CASE::OTHER; }                                          // Not A, B, C, or NONE


    void                        ApplyMassTransferRejuvenationFactor()                                           { m_Age *= CalculateMTRejuvenationFactor(); }             // Apply age rejuvenation factor

    void                        UpdateComponentVelocity(const Vector3d p_NewVelocity)           { m_ComponentVelocity += p_NewVelocity; }
    void                        UpdateMassTransferDonorHistory();
  

    
virtual double CalculateCriticalMassRatio(const double p_Mass,
                                          const double p_Radius,
                                          const double p_CoreMass,
                                          const bool   p_AccretorIsDegenerate,
                                          const double p_MTefficiency = 0.0) const;

virtual double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const { return 0.0; } // Default is 0.0

virtual double CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency = 0.0) const { return 0.0; } 

virtual double CalculateCriticalMassRatio_Ge2020(const double p_MTefficiency) {
    return CalculateCriticalMassRatio_Ge2020_Interpolate(Mass(), Radius(), p_MTefficiency);
}

virtual double CalculateCriticalMassRatio_Hurley2002() const { return 0.0; } // Default is 0.0




inline COMPAS_PURE double CalculateZetaAdiabatic_Hurley2002(const double p_Mass, const double p_CoreMass) const;

GNU_CONST double CalculateZetaAdiabatic_Soberman1997(const double p_CoreMass) const;




    
    virtual DBL_DBL_DBL_DBL     CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const;
    virtual DBL_DBL_DBL_DBL     CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const ;
    virtual DBL_DBL_DBL_DBL     CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const;
    

    virtual double              CalculateCELambda_Dewi() const                                                     { return 1.0; }           

                                              // Default for stellar types with no LamdaDewi definitions - 1.0 is benign

GNU_CONST double CalculateCELambdaKruckow(const double p_Radius, const double p_Alpha) const;
    double CalculateCELambdaKruckow() const { return CalculateCELambdaKruckow(m_Radius, OPTIONS->CommonEnvelopeSlopeKruckow()); }
    
    virtual double              CalculateLambdaLoveridge(const double p_EnvMass, const bool p_IsMassLoss = false) const { return 1.0; }                                             // Default for non giant branch stars - 1.0 is benign
    double                      CalculateLambdaLoveridge() const                                                { return CalculateLambdaLoveridge(m_Mass - m_CoreMass, false); }


    double                      CalculateCELambda_Nanjing() const;
    


    DBL_DBL                     CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate); 
    virtual DBL_DBL             CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                            const double p_AccretorMassRate,
                                                            const bool   p_IsHeRich)                            { return CalculateMassAcceptanceRate(p_DonorMassRate, p_AccretorMassRate); } // Ignore the He content for non-WDs
    double                      CalculateMassAccretedForCO(const double p_Mass, const double p_CompanionMass, const double p_CompanionRadius, const double p_CompanionEnvelope) const;




COMPAS_PURE  virtual double CalculateCELambda_Nanjing_StarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const { return 1.0; } // Default for stellar types with no LamdaNanjing definitions - 1.0 is benign




GNU_PURE  virtual double CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
                                                            const double             p_Radius,
                                                            const double             p_CoreMass,
                                                            const size_t             p_MassIndex,
                                                            const STELLAR_POPULATION p_StellarPop) const { return 1.0; } // Default for stellar types with no LamdaNanjing definitions - 1.0 is benign



    double                      CalculateThermalMassAcceptanceRate(const double p_Radius);
    double                      CalculateThermalMassAcceptanceRate()                                            { return CalculateThermalMassAcceptanceRate(m_Radius); }


    
    double                      CalculateZetaAdiabatic() const;
    virtual double              CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription)          { return 0.0; }                                                     // Use inheritance hierarchy
    virtual double              CalculateZetaEquilibrium()                                                      { return 0.0; }



    virtual ACCRETION_REGIME    DetermineAccretionRegime(const double p_DonorThermalMassLossRate, 
                                                         const bool p_HeRich)                                   { return ACCRETION_REGIME::ZERO; }                                  // Placeholder, use inheritance for WDs



    void                        ResolveAccretion(const double p_AccretionMass)                                  { m_Mass = std::max(0.0, m_Mass + p_AccretionMass); }               // Handles donation and accretion - won't let mass go negative
    virtual void                ResolveAccretionRegime(const ACCRETION_REGIME p_Regime, const double p_DonorThermalMassLossRate) { }                                                // Default does nothing, only works for WDs.
    virtual double              ResolveCommonEnvelopeAccretion(const double p_FinalMass,
                                                               const double p_CompanionMass     = 0.0,
                                                               const double p_CompanionRadius   = 0.0,
                                                               const double p_CompanionEnvelope = 0.0)          { return p_FinalMass - Mass(); }                                    // Overwritten in NS.h; for now, no accretion on stars other than compact objects during CE



    virtual void                UpdateAfterMerger(double p_Mass, double p_HydrogenMass) { }                                                                                         // Default is NO-OP
//    virtual void                UpdateAgeAfterMassLoss() { }                                                                                                                        // Default is NO-OP




GNU_CONST double  CalculateCELambda_Nanjing_MassInterpolated(const double p_Mass, const double p_Radius, const double p_CoreMass, const STELLAR_POPULATION stellarPop) const;



GNU_PURE  double CalculateCELambda_Nanjing_MassInterpolated(const double p_Mass, const STELLAR_POPULATION p_StellarPop) const;
GNU_PURE  double CalculateCELambda_Nanjing_ZInterpolated(const double p_Z, const int p_MassIndex) const;



/*
 * CalculateConvectiveEnvelopeBindingEnergy
 *
 * @brief
 * Calculate the convective envelope binding energy for the two-stage common envelope formalism of Hirai & Mandel, 2022
 *
 *
 * double CalculateConvectiveEnvelopeBindingEnergy(const double p_Mass, const double p_Radius, const double p_EnvMass, const double p_Lambda) const
 *
 * @param           p_Mass                      Mass (total) of the star (Msol)
 * @param           p_Radius                    Radius of the star (Rsol)
 * @param           p_EnvMass                   Mass of the convective outer envelope  (Msol)
 * @param           p_Lambda                    Lambda parameter for the convective envelope
 * @return                                      Binding energy (erg)
 */
inline double BaseStar::CalculateConvectiveEnvelopeBindingEnergy(const double p_Mass, const double p_Radius, const double p_EnvMass, const double p_Lambda) const {
    return CalculateBindingEnergy(p_Mass - p_EnvMass, p_Radius, p_EnvMass, p_Lambda);
}


/*
 * CalculateConvectiveEnvelopeLambda_Picker
 *
 * @brief
 * Calculates (an approximate value of) the lambda binding energy parameter of the outer convective envelope.
 * Follows the fits of Picker et al., 2024, for lambda_He (https://arxiv.org/abs/2402.13180).
 *
 * This is required for the Hirai & Mandel, 2022, two-stage CE formalism.
 *
 *
 * double BaseStar::CalculateConvectiveEnvelopeLambda_Picker(const double p_Mass, const double p_EnvMass, const double p_EnvMassMax) const
 *
 * @param           p_Mass                      Mass of the star (Msol)
 * @param           p_EnvMass                   Mass of the star's convective outer envelope (Msol)
 * @param           p_EnvMassMax                Maximum mass of the star's convective outer envelope (Msol)
 * @return                                      Lambda binding energy parameter for the convective outer envelope
 */
inline double BaseStar::CalculateConvectiveEnvelopeLambda_Picker(const double p_Mass, const double p_EnvMass, const double p_EnvMassMax) const {
    
    const double m2 = 0.0023 * GLOBALS->SigmaHurley() * GLOBALS->SigmaHurley() + 0.0088 * GLOBALS->SigmaHurley() + 0.013;   // Picker et al., 2024, eq 12, and Table 1
    const double b1 = m2 * p_Mass - 0.23;                                                                                   // ibid., eq 11
    
    return std::exp((p_EnvMass / p_EnvMassMax > 0.3 ? 0.42 * p_EnvMass / p_EnvMassMax + b1 : 0.3 * 0.42 + b1));
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
    return (p_CoreMass >= m_Mass) ? : -GLOBALS->HurleyRadiusXexponent() + (2.0 * utils::IntPow(p_CoreMass / m_Mass, 5));
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
	const double beta = max(-2.0 / 3.0, min(-1.0, p_Alpha));    // beta in Vigna-Gomez et al. 2018, sec 3.2.8; clamp to [-1.0, -2.0/3.0]
	return 1600.0 * PPOW(0.00125, -beta) * PPOW(p_Radius, beta);
}


#endif // __BaseStar_h__
