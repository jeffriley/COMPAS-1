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

    







    /*
     * The following Clone() functions should be used to clone a star - any steller type, including BaseStar.
     * The BaseStar functions will never actually be executed - the derived functions will be called as required.
     * The static function is provided for cases where a clone is required (e.g. for hypothesis testing), but
     * no object of the correct class exists.
     *
     * Important:
     *
     * The Clone() functions return a pointer, created by the 'new' operator.  The 'new' operator dynamically
     * allocates memory on the heap, not the stack, which is why it is available to the caller of this function
     * after this function has exited and its stack frame collapsed.  It is the responsibility of the caller of
     * this function to delete the pointer returned when it is no longer required so that the allocated memory
     * is return to the pool of available memory - failing to do so will cause a memory leak and the program will
     * eventually exhaust available memory and fail.  The preferred usage pattern is:
     *
     *     T* ptr = Clone(obj, persistence)
     *     ...
     *     ...
     *     delete ptr; ptr = nullptr;
     *
     *
     * BaseStar* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true)
     *
     * @param   [IN]    p_Star                      (address of) The star to be cloned
     *                                              Must be BaseStar or one of the stellar type classes derived from the BaseStar class
     * @param   [IN]    p_Persistence               Specifies the object persistence to be assigned to the cloned star.
     *                                              If the cloned star is intended to be used temporarily (e.g. for hypothesis testing),
     *                                              persistence should be EPHEMERAL (the default), otherwise PERMANENT.
     * @param   [IN]    p_Initialise                Specifies whether the clone should be initialised via the class Initialise() function.
     *                                              (Default is `true`)
     * @return                                      (pointer to) The cloned star
     */
    
    virtual BaseStar* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) { return this; }
    static  BaseStar* Clone(BaseStar* p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) { return p_Star; }
    

    
    
    // object identifiers - all classes have these
    OBJECT_ID                   ObjectId() const                                        { return m_ObjectId; }
    OBJECT_TYPE                 ObjectType() const                                      { return OBJECT_TYPE::BASE_STAR; }
    OBJECT_PERSISTENCE          ObjectPersistence() const                               { return m_ObjectPersistence; }



//    STELLAR_TYPE                StartingStellarType() const                             { return m_StartingStellarType; }
//    STELLAR_TYPE                StellarType() const                                     { return m_StellarType; }
//    STELLAR_TYPE                StellarTypePrev() const                                 { return m_StellarTypePrev; }
    
    
    // Getters - alphabetically
    //
    // Short (one or two line) implementations that don't require a lot of documentation (i.e. that won't clutter here), are in here (and are implied inline)
    // Short implementations that can still be inlined, but that require more documentation (and would clutter here), are below (after class declaration)
    // Longer implementations that can still be inlined (that would clutter here) are below (after class declaration)
    // All other implementations are in the cpp file
    //
    // *** getter implementations should generally be very short

    double                      Age() const                                             { return m_Age; }
    double                      AngularMomentum() const                                 { return m_AngularMomentum; }
    bool                        CHonMS() const                                          { return m_CHE; }
    double                      COCoreMass() const                                      { return m_COCoreMass; }
    double                      CoreMass() const                                        { return m_CoreMass; }
    int                         DominantMassLossType() const                            { return static_cast<int>(m_DominantMassLossType); }
    double                      Dt() const                                              { return m_dt; }
    double                      DtPrev() const                                          { return m_dtPrev; }
    bool                        EnvelopeJustExpelledByPulsations() const                { return m_EnvelopeJustExpelledByPulsations; }
    ERROR                       Error() const                                           { return m_Error; }
    EVOLUTION_STATUS            EvolutionStatus() const                                 { return m_EvolutionStatus; }
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
    double                      HeCoreMass() const                                      { return m_HeCoreMass; }
    double                      HeliumAbundanceCore() const                             { return m_HeliumAbundanceCore; }
    double                      HeliumAbundanceSurface() const                          { return m_HeliumAbundanceSurface;}
    double                      HydrogenAbundanceCore() const                           { return m_HydrogenAbundanceCore; }
    double                      HydrogenAbundanceSurface() const                        { return m_HydrogenAbundanceSurface; }
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
    double                      LogMetallicityRho() const                               { return LogMetallicityXiHurley() + 1.0; }                                  // rho in Hurley+ 2000
    double                      LogMetallicitySigma() const                             { return m_Log10Metallicity; }                                              // sigma in Hurley+ 2000
    double                      LogMetallicityXiHurley() const                          { return m_Log10Metallicity - LOG10_ZSOL_HURLEY; }                          // xi in Hurley+ 2000
    double                      LogMetallicityXiAnders() const                          { return m_Log10Metallicity - LOG10_ZSOL_ANDERS; }                          // log10(Z / ZSOL_ANDERS)
    double                      LogMetallicityXiAsplund() const                         { return m_Log10Metallicity - LOG10_ZSOL_ASPLUND; }                         // log10(Z / ZSOL_ASPLUND)
    double                      Luminosity() const                                      { return m_Luminosity; }
    double                      MainSequenceCoreMass() const                            { return m_MainSequenceCoreMass; }
    double                      Mass() const                                            { return m_Mass; }
    double                      Mass0() const                                           { return m_Mass0; }
    double                      MassPrev() const                                        { return m_MassPrev; }
    ST_VECTOR                   MassTransferDonorHistory() const                        { return m_MassTransferDonorHistory; }
    double                      Mdot() const                                            { return m_Mdot; }
    double                      Metallicity() const                                     { return m_Metallicity; }
    double                      MZAMS() const                                           { return m_MZAMS; }
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
    double                      Radius() const                                          { return m_Radius; }
    double                      RadiusPrev() const                                      { return m_RadiusPrev; }
    unsigned long int           RandomSeed() const                                      { return m_RandomSeed; }
    double                      RZAMS() const                                           { return m_RZAMS; }
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
    double                      Tau() const                                             { return m_Tau; }
    double                      Temperature() const                                     { return m_Temperature; }
    double                      Time() const                                            { return m_Time; }
    double                      Timescale(TIMESCALE p_Timescale) const                  { return m_Timescales[static_cast<int>(p_Timescale)]; }
    double                      TotalMassLossRate() const                               { return m_TotalMassLossRate; }
    double                      TZAMS() const                                           { return m_TZAMS; }
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
    void                        SetError(const ERROR p_Error)                                   { m_Error = p_Error; }
    void                        SetObjectId(const OBJECT_ID p_ObjectId)                         { m_ObjectId = p_ObjectId; }
    void                        SetPersistence(const OBJECT_PERSISTENCE p_Persistence)          { m_ObjectPersistence = p_Persistence; }
    
    void                        SetOmega(double p_Omega)                                        { SetAngularMomentum(CalculateMomentOfInertiaAU() * p_Omega); }
    

GNU_CONST inline SN_EVENT AddSNEvent(const SN_EVENT p_SNEvent, const SN_EVENT p_ExistingSNEvent) { return p_ExistingSNEvent |= p_SNEvent; }
GNU_CONST inline SN_EVENT ClearSNEvent(const SN_EVENT p_SNEvent) { return SN_EVENT::NONE; }
    

    void                        SetEvolutionStatus(const EVOLUTION_STATUS p_EvolutionStatus)    { m_EvolutionStatus = p_EvolutionStatus; }                          // Set evolution status (typically final outcome) for star

    void                        SetDt(const double p_dt)                                        { m_dt = std::max(0.0, p_dt); }                                     // Set timestep - ensure >= 0.0
    void                        SetPrevDt(const double p_dt)                                    { m_dtPrev = p_dt; }                                                // Set previous timestep - don't clamp - preserve actuality
  
    virtual void                UpdateEffectiveZAMSLandR(){}                                                                                                        // Virtual - default is NO-OP
    

    // Public member function prototypes - alphabetically
    //
    // Short (one or two line) implementations that don't require a lot of documentation (i.e. that won't clutter here), are in here (and are implied inline)
    // Short implementations that can still be inlined, but that require more documentation (and would clutter here), are below (after class declaration)
    // Longer implementations that can still be inlined (that would clutter here) are below (after class declaration)
    // All other implementations are in the cpp file

    // Function names should indicate the purpose of the function (and should genarally begin with a verb)
    //
    // Functions beginning with:
    //
    //    - "Calculate" should calculate and return value(s) - they should not update (modify) class member variables
    //    - "Determine" should determine e.g. current state etc., and return value(s) - they should not update (modify) class member variables
    //
    //    - "Apply"     should update (modify) class member variables, and may return value(s)  /// NO!!!!!! ... maybe <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //    - "Update"    should update (modify) class member variables, and may return value(s)  /// NO!!!!!! ... maybe <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    //
    // Functions that do not update (modify) class member variables should be declared const
    

GNU CONST double CalculateBindingEnergy(const double p_CoreMass, const double p_EnvMass, const double p_Radius, const double p_Lambda) const;
    double       CalculateBindingEnergy(const double p_Lambda) const                                             { return CalculateBindingEnergy(m_CoreMass, m_Mass - m_CoreMass, m_Radius, p_Lambda); }

    virtual double              CalculateConvectiveCoreMass() const                                             {return 0.0;}



GNU_CONST virtual inline double CalculateConvectiveCoreRadius(const double p_Radius, const double p_Tau) const {return 0.0;}

virtual double              CalculateRemnantRadius() const                                                  { return Radius(); }                                                // Relevant for MS stars, over-written for GB stars


GNU_CONST CalculateConvectiveEnvelopeBindingEnergy(const double p_Mass, const double p_EnvMass, const double p_Radius, const double p_Lambda) const
GNU_PURE  double CalculateConvectiveEnvelopeLambda_Picker(const double p_Mass, const double p_EnvMass, const double p_EnvMassMax) const;
    
    
    
    virtual DBL_DBL             CalculateConvectiveEnvelopeMass() const                                         { return std::tuple<double, double> (0.0, 0.0); }





    double                      CalculateDynamicalTimescale() const                                             { return CalculateDynamicalTimescale_Static(m_Mass, m_Radius); }    // Use class member variables
    
    
    double                      CalculateEddyTurnoverTimescale() const;
    



    double                      CalculateMassChangeTimescale() const                                            { return CalculateMassChangeTimescale_Static(m_StellarType, m_StellarTypePrev, m_Mass, m_MassPrev, m_dtPrev); }  // Use class member variables
    double                      CalculateMassLossValues(double p_Dt, const bool p_UpdateMDot = false);

    virtual double              CalculateMomentOfInertia() const                                                { return (0.1 * (m_Mass) * m_Radius * m_Radius); }                  // Defaults to MS. k2 = 0.1 as defined in Hurley et al. 2000, after eq 109
    virtual double              CalculateMomentOfInertiaAU() const                                              { return CalculateMomentOfInertia() * RSOL_TO_AU * RSOL_TO_AU; }
    
    double                      CalculateNuclearMassLossRate()                                                  { return m_Mass / CalculateRadialExpansionTimescaleDuringMassTransfer(); }
        
    GNU_PURE double CalculateCHEAngularFrequency_Static(const double p_MZAMS, const double p_Metallicity) const;

    double                      CalculateRadialChange() const                                                   { return (utils::Compare(m_RadiusPrev,0) <= 0)? 0 : std::abs(m_Radius - m_RadiusPrev) / m_RadiusPrev; } // Return fractional radial change (if previous radius is negative or zero, return 0 to avoid NaN
    double                      CalculateRadialExpansionTimescale() const                                       { return CalculateRadialExpansionTimescale_Static(m_StellarType, m_StellarTypePrev, m_Radius, m_RadiusPrev, m_dtPrev); } // Use class member variables
    double                      CalculateRadialExpansionTimescaleDuringMassTransfer();
    virtual double              CalculateRadialExtentConvectiveEnvelope() const                                 { return 0.0; }                                                     // Default for stars with no convective envelope
    virtual double              CalculateRadiusOnMassChange(double p_dM)                                        { return Radius(); }                                                // NO-OP



    double                      CalculateSNKickMagnitude(const double p_RemnantMass, const double p_EjectaMass, const STELLAR_TYPE p_StellarType);


    virtual double              CalculateThermalMassLossRate() const                                            { return m_Mass / CalculateThermalTimescale(); }                    // Use class member variables - and inheritance hierarchy
    virtual double              CalculateThermalTimescale(const double p_Radius) const;                                                                                             // Use inheritance hierarchy
    virtual double              CalculateThermalTimescale() const                                               { return CalculateThermalTimescale(m_Radius); }                     // Use inheritance hierarchy

    virtual double              CalculateTAMSCoreMass() const                                                   { return 0.0; }                                                     // Except MS stars

    double                      CalculateTimestep();



    virtual ENVELOPE            DetermineEnvelopeType() const                                                   { return ENVELOPE::REMNANT; }                                       // Default is REMNANT - but should never be called

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



/*
    STELLAR_TYPE            m_StartingStellarType;                      // Stellar type at start of simulation, defined in Hurley et al. 2000
    STELLAR_TYPE            m_StellarType;                              // Stellar type defined in Hurley et al. 2000
    ERROR                   m_Error;                                    // Records most recent error encountered for this star
*/
    
//    StarStateT stateZero;
//    StarStateT stateZAMS;  // << ---- flag this as not set - until it is set (if it is set)

/*

    // Stellar variables
    bool                    m_CHE;                                      // CHE flag - true if the star spent entire MS as a CH star; false if evolved CH->MS
    EVOLUTION_STATUS        m_EvolutionStatus;                          // Status of evolution for this star (typically final outcome e.g. DONE, TIMES_UP, etc.)
    unsigned long int       m_RandomSeed;                               // Seeds the random number generator for this star

    // Zero Age Main Sequence
    double                  m_InitialHeliumAbundance;                   // Initial helium abundance (Y)
    double                  m_InitialHydrogenAbundance;                 // Initial hydrogen abundance (X)
    double                  m_InitialMainSequenceCoreMass;              // Initial main sequence core mass (used in BRCEK core mass prescription)
    double                  m_LZAMS;                                    // ZAMS Luminosity
    double                  m_MZAMS;                                    // ZAMS Mass
    double                  m_OmegaZAMS;                                // ZAMS Angular Frequency
    double                  m_OmegaCHE;                                 // Minimum angular frequency at which CHE will occur (calculated at ZAMS)
    double                  m_RZAMS;                                    // ZAMS Radius
    double                  m_TZAMS;                                    // ZAMS Temperature

    // Effective Zero Age Main Sequence
    double                  m_LZAMS0;                                   // Effective ZAMS Luminosity
    double                  m_RZAMS0;                                   // Effective ZAMS Radius

    // Current timestep variables
    double                  m_Age;                                      // Current effective age (changes with mass loss/gain) (Myr)
    double                  m_AngularMomentum;                          // Angular Momentum (Msol * AU^2 / yr)
    double                  m_COCoreMass;                               // Current CO core mass (Msol)
    double                  m_CoreMass;                                 // Current core mass (Msol)
    double                  m_dt;                                       // Size of current timestep (Myr)
    bool                    m_EnvelopeJustExpelledByPulsations;         // Flag to know if the convective envelope has just been expelled by pulsations
    double                  m_HeCoreMass;                               // Current He core mass (Msol)
    double                  m_HeliumAbundanceCore;                      // Helium abundance in the core
    double                  m_HeliumAbundanceSurface;                   // Helium abundance at the surface
    double                  m_HydrogenAbundanceCore;                    // Hydrogen abundance in the core
    double                  m_HydrogenAbundanceSurface;                 // Hydrogen abundance at the surface
    bool                    m_LBVphaseFlag;                             // Flag to know if the star satisfied the conditions, at any point in its evolution, to be considered a Luminous Blue Variable (LBV)
    double                  m_Luminosity;                               // Current luminosity (Lsol)
    double                  m_MainSequenceCoreMass;                     // Core mass of main sequence stars (Msol)
    double                  m_Mass;                                     // Current mass (Msol)
    double                  m_Mass0;                                    // Current effective initial mass (Msol)
    double                  m_MinimumLuminosityOnPhase;                 // Only required for CHeB stars, but only needs to be calculated once per star
    double                  m_Mdot;                                     // Current mass loss rate in winds (Msol per yr)
    MASS_LOSS_TYPE          m_DominantMassLossRate;                     // Current dominant type of wind mass loss

    double                  m_Mu;                                       // Current small envelope parameter mu
    double                  m_Radius;                                   // Current radius (Rsol)
    double                  m_Tau;                                      // Relative time
    double                  m_Temperature;                              // Current temperature (Tsol)
    double                  m_Time;                                     // Current physical time the star has been evolved (Myr)
    double                  m_TotalMassLossRate;                        // Current mass loss/gain rate from mass transfer or winds (Msol per yr)

    // Previous timestep variables
    double                  m_dtPrev;                                   // Previous timestep
    double                  m_MassPrev;                                 // Previous mass (Msol)
    double                  m_RadiusPrev;                               // Previous radius (Rsol)
    STELLAR_TYPE            m_StellarTypePrev;                          // Stellar type at previous timestep

    // Metallicity variables
    double                  m_Metallicity;                              // Metallicity
    double                  m_Log10Metallicity;                         // log10(Metallicity) - for performance
*/

/*
    // Metallicity dependent constants
    double                  m_Alpha1;                                   // Alpha1 in Hurley et al. 2000, just after eq 49
    double                  m_Alpha3;                                   // Alpha3 in Hurley et al. 2000, just after eq 56
    double                  m_Alpha4;                                   // Alpha4 in Hurley et al. 2000, just after eq 57
    double                  m_XExponent;                                // Exponent to which R depends on M - 'x' in Hurley et al. 2000, eq 47

    // Constants only calculated once
    double                  m_BaryonicMassOfMaximumNeutronStarMass;     // Baryonic mass of MaximumNeutronStarMass 

    // The following vectors were initially implemented as unordered_maps.  The code worked
    // quite well, except for one small problem - access times (presumably due to hashing)
    // were prohibitive when accessed hundreds of thousands, and in some cases, millions,
    // of times as we evolve the star.  So vectors are now used instead - the code is not as
    // elegant, but performance is better by an order of magnitude

    // Timescales, Giant Branch parameters, mass cutoffs
    DBL_VECTOR              m_GBparams;                                 // Giant Branch Parameters
    DBL_VECTOR              m_MassCutoffs;                              // Mass cutoffs
    DBL_VECTOR              m_Timescales;                               // Timescales

    // Luminosity, Radius, a(n) and b(n) coefficients
    DBL_VECTOR              m_AnCoefficients;                           // a(n) coefficients
    DBL_VECTOR              m_BnCoefficients;                           // b(n) coefficients
    DBL_VECTOR              m_LCoefficients;                            // Luminosity coefficients
    DBL_VECTOR              m_RCoefficients;                            // Radius coefficients

    // Luminosity, Radius and Gamma constants
    // These are calculated in CalculateAnCoefficients()
    // Calculating the a(n) coefficients requires one of the R constants, and the L, R and Gamma
    // constants are calculated using the a(n) coefficients
    DBL_VECTOR              m_GammaConstants;                           // Gamma constants
    DBL_VECTOR              m_LConstants;                               // Luminosity constants
    DBL_VECTOR              m_RConstants;                               // Radius constants
*/

/*
    // Stellar details squirrelled away...
    SupernovaDetailsT       m_SupernovaDetails;                         // Supernova attributes
    PulsarDetailsT          m_PulsarDetails;                            // Pulsar attributes

    // Star vector velocity 
    Vector3d                m_ComponentVelocity;                        // Isolated star velocity vector (binary's centre-of-mass velocity for bound binary)

    // Star mass transfer history 
    ST_VECTOR               m_MassTransferDonorHistory;                 // List of MT donor stellar types - mostly relevant for binary stars
*/

    // Protected member function prototypes - alphabetically (implementations are in the cpp file)
    // See below (after class declaration) for inline functions declared in this header file

    double              CalculateAlpha1() const;                        // inline implementation below
    double              CalculateAlpha3() const;                        // inline implementation below
    double              CalculateAlpha4() const;                        // inline implementation below

    void                CalculateAnCoefficients(DBL_VECTOR &p_An, DBL_VECTOR &p_LConstants, DBL_VECTOR &p_RConstants, DBL_VECTOR &p_GammaConstants) const;



    void                CalculateBnCoefficients(DBL_VECTOR &p_Bn) const;



    double       CalculateConvectiveEnvelopeBindingEnergy(const double p_Lambda) const                           { 
                                                                                                                            double convectiveEnvMass;
                                                                                                                            std::tie(convectiveEnvMass, std::ignore) = CalculateConvectiveEnvelopeMass();
                                                                                                                            return CalculateConvectiveEnvelopeBindingEnergy(m_Mass, convectiveEnvMass, m_Radius, p_Lambda);
                                                                                                                        }

                                                     // Default is NO-OP

    static  double      CalculateDynamicalTimescale_Static(const double p_Mass, const double p_Radius);

    virtual double      CalculateEddingtonCriticalRate() const                                                          { return 2.08E-3 / 1.7 * m_Radius * MYR_TO_YEAR * OPTIONS->EddingtonAccretionFactor() ; } // Hurley+, 2002, Eq. (67)






    double              CalculateGBRadiusXExponent() const;

    virtual double      CalculateHeCoreMassAtPhaseEnd() const                                                           { return m_HeCoreMass; }                                                    // Default is NO-OP
    virtual double      CalculateHeCoreMassOnPhase() const                                                              { return m_HeCoreMass; }                                                    // Default is NO-OP









    void                CalculateLCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_LCoefficients) const;


GNU_CONST double CalculateLifetimeToBAGB_Hurley2000(const double p_tHeI, const double p_tHe) const;
GNU_CONST double CalculateLifetimeToBGB_Hurley2000(const double p_Mass, const DBL_VECTOR& p_aN) const;





    double              CalculateLuminosityAtBAGB(double p_Mass) const;         // inline implementation below

    virtual double      CalculateLuminosityAtPhaseEnd() const                                                           { return m_Luminosity; }                                                    // Default is NO-OP



    
    
    


    virtual double      CalculateLuminosityOnPhase() const                                                              { return m_Luminosity; }                                                    // Default is NO-OP




    static  double      CalculateMassChangeTimescale_Static(const STELLAR_TYPE p_StellarType,
                                                            const STELLAR_TYPE p_StellarTypePrev,
                                                            const double       p_Mass,
                                                            const double       p_MassPrev,
                                                            const double       p_DtPrev);
    
    void                CalculateMassCutoffs(const double p_Metallicity, const double p_LogMetallicityXi, DBL_VECTOR &p_MassCutoffs) const;





























   
    
    
    
    
    



    double              CalculateMassLossRateWolfRayetTemperatureCorrectionSander2023(const double p_Mdot) const;



// JR : DON'T FORGET "OVERRIDE" ON DERIVED FUNCTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





GNU_CONST double CalculateEddingtonLuminosity(const double p_Mass, const double p_HeAbundanceSurface) const;



virtual double CalculateRadiusOnPhase() const { return m_StateHistory.CurrentState.Radius(); }




GNU_CONST double CalculateTemperatureOnPhase() const;
GNU_CONST double CalculateTemperatureOnPhase(const double p_Luminosity, const double p_Radius) const;



GNU_PURE virtual MASS_LOSS_T CalculateMLRate_Belczynski2010(const double                     p_Metallicity,
                                                            const double                     p_Mass,
                                                            const double                     p_Radius,
                                                            const double                     p_Luminosity,
                                                            const double                     p_Temperature,
                                                            const double                     p_PerturbationMu,
                                                            const double                     p_ZscaledHurley,
                                                            const double                     p_HeAbundanceSurface,
                                                            const double                     p_CoolWindsMultiplier,
                                                            const double                     p_WRfactor,
                                                            const bool                       p_ScaleWithSurfaceHe,
                                                            const LBV_MASS_LOSS_PRESCRIPTION p_LBVprescription) const;

GNU_CONST virtual MASS_LOSS_T CalculateMLRate_Hurley2000(const double p_Mass,
                                                         const double p_Radius,
                                                         const double p_Luminosity,
                                                         const double p_PerturbationMu,
                                                         const double p_ZscaledHurley,
                                                         const double p_WRfactor) const;

COMPAS_PURE virtual MASS_LOSS_T CalculateMLRate_Merritt2025(const double p_Metallicity,
                                                            const double p_Mass,
                                                            const double p_Radius,
                                                            const double p_Luminosity,
                                                            const double p_Temperature,
                                                            const double p_PerturbationMu,
                                                            const double p_mStart,
                                                            const double p_SigmaHurley,
                                                            const double p_ZetaAnders,
                                                            const double p_ZetaAsplund,
                                                            const double p_ZscaledHurley,
                                                            const double p_HeAbundanceSurface,
                                                            const double p_WRfactor,
                                                            const double p_TerminalWindScalePower) const;

GNU_CONST MASS_LOSS_T CalculateMLRate_NieuwenhuijzenDeJager1990(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_ZscaledHurley) const;

GNU_CONST MASS_LOSS_T CalculateMLRate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const


GNU_CONST MASS_LOSS_T CalculateMLRateLBV(const double p_Radius, const double p_Luminosity, const LBV_MASS_LOSS_PRESCRIPTION p_LBVprescription, const double p_LBVfactor) const;
GNU_CONST MASS_LOSS_T CalculateMLRateLBV_Belczynski2010(const double p_LBVfactor) const;
GNU_CONST MASS_LOSS_T CalculateMLRateLBV_Hurley2000(const double p_Luminosity, const double p_HDlimitfactor) const;


GNU_CONST MASS_LOSS_T CalculateMLRateOB(const double                    p_Metallicity,
                                        const double                    p_Mass,
                                        const double                    p_Luminosity,
                                        const double                    p_Temperature,
                                        const double                    p_ZetaAnders,
                                        const double                    p_ZetaAsplund,
                                        const double                    p_TerminalWindScalePower,
                                        const OB_MASS_LOSS_PRESCRIPTION p_MassLossPrescription) const;
GNU_CONST MASS_LOSS_T CalculateMLRateOB_Bjorklund2022(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
GNU_CONST MASS_LOSS_T CalculateMLRateOB_Krticka2018(const double p_Luminosity, const double p_ZetaAsplund) const;
GNU_CONST MASS_LOSS_T CalculateMLRateOB_Vink2001(const double p_Mass,
                                                 const double p_Luminosity,
                                                 const double p_Temperature,
                                                 const double p_ZetaAnders,
                                                 const double p_TerminalWindScalePower) const;
GNU_CONST MASS_LOSS_T CalculateMLRateOB_VinkSander2021(const double p_Mass, const double p_Luminosity, const double p_Temperature, const double p_ZetaAnders) const;

GNU_CONST MASS_LOSS_T CalculateMLRateRSG(const double                     p_Mass,
                                         const double                     p_Radius,
                                         const double                     p_Luminosity,
                                         const double                     p_mStart,
                                         const double                     p_Temperature,
                                         const double                     p_ZscaledHurley,
                                         const RSG_MASS_LOSS_PRESCRIPTION p_MassLossPrescription) const;
GNU_CONST MASS_LOSS_T CalculateMLRateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const;
GNU_CONST MASS_LOSS_T CalculateMLRateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const;
GNU_CONST MASS_LOSS_T CalculateMLRateRSG_Kee2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
GNU_CONST MASS_LOSS_T CalculateMLRateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const;
GNU_CONST MASS_LOSS_T CalculateMLRateRSG_Yang2023(const double p_Luminosity) const;

GNU_CONST MASS_LOSS_T CalculateMLRateVMS(const double                     p_Metallicity,
                                         const double                     p_Mass,
                                         const double                     p_Luminosity,
                                         const double                     p_Temperature,
                                         const double                     p_SigmaHurley,
                                         const double                     p_ZetaAnders,
                                         const double                     p_ZetaAsplund,
                                         const double                     p_TerminalWindScalePower,
                                         const OB_MASS_LOSS_PRESCRIPTION  p_OBprescription,
                                         const VMS_MASS_LOSS_PRESCRIPTION p_VMSprescription) const;
GNU_CONST MASS_LOSS_T CalculateMLRateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const;
GNU_CONST MASS_LOSS_T CalculateMLRateVMS_Sabhahit2023(const double                    p_Metallicity,
                                                      const double                    p_Mass,
                                                      const double                    p_Luminosity,
                                                      const double                    p_Temperature,
                                                      const double                    p_SigmaHurley,
                                                      const double                    p_ZetaAnders,
                                                      const double                    p_ZetaAsplund,
                                                      const OB_MASS_LOSS_PRESCRIPTION p_OBprescription) const;
GNU_CONST MASS_LOSS_T CalculateMLRateVMS_Vink2011(const double p_Mass,
                                                  const double p_Luminosity,
                                                  const double p_Temperature,
                                                  const double p_ZetaAnders,
                                                  const double p_TerminalWindScalePower) const;

GNU_CONST MASS_LOSS_T CalculateMLRateWR_Hurley2000(const double p_Luminosity, const double p_PerturbationMu) const;
static GNU_CONST MASS_LOSS_T CalculateMLRateWR_SanderVink2020_Static(const double p_Luminosity, const double p_PerturbationMu, const double p_ZetaAnders) const;
static GNU_CONST MASS_LOSS_T CalculateMLRateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature, const double p_SigmaHurley) const;
static GNU_CONST MASS_LOSS_T CalculateMLRateWR_ZDependent_Static(const double p_Metallicity, const double p_Luminosity, const double p_PerturbationMu) const;




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 RADIUS FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////







virtual double      CalculateRadiusAtPhaseEnd() const                                                               { return m_Radius; }                                                        // Default is NO-OP






    double CalculateRadiusOnPhase(const double p_Mass, const double p_Tau, const double p_RZAMS) const;














// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<















    


    virtual double      CalculateMTRejuvenationFactor()                                                       { return 1.0; }

    double       CalculateMaximumCoreMass(double p_Mass) const                                                   { return std::min(((1.45 * p_Mass) - 0.31), p_Mass); }                      // Hurley et al. 2000, eq 89

    /*
     * CalculateOmegaBreak
     *
     * Calculate the break up angular velocity of a star, in rad yr^1, where [G] = 4*pi^2 AU^3 yr^-2 Msol^-1
     *
     * Uses current values of:
     * 
     *    - m_Radius
     * 
     * 
     * double CalculateOmegaBreak() const
     *
     * @return                                      Break up angular velocity (rad yr^-1)
     */
    double CalculateOmegaBreak() const {
	    return _2_PI * std::sqrt(m_Mass / (RSOL_TO_AU * RSOL_TO_AU * RSOL_TO_AU * m_Radius * m_Radius * m_Radius));
    }


    static  double      CalculateOStarRotationalVelocityAnalyticCDF_Static(const double p_Ve);
    static  double      CalculateOStarRotationalVelocityAnalyticCDFInverse_Static(double p_Ve, void *p_Params);
    static  double      CalculateOStarRotationalVelocity_Static(const double p_Xmin, const double p_Xmax);

    double       CalculatePerturbationB(const double p_Mass) const                                               { return 0.002 * std::max(1.0, (2.5 / p_Mass)); }                           // Hurley et al. 2000, eq 103
    double       CalculatePerturbationC(double p_Mass) const                                                     { return 0.006 * std::max(1.0, (2.5 / p_Mass)); }                           // Hurley et al. 2000, eq 104
    virtual double      CalculatePerturbationMu() const                                                                 { return m_Mu; }                                                            // Default is NO-OP
    virtual double      CalculatePerturbationMuAtPhaseEnd() const                                                       { return CalculatePerturbationMuOnPhase(); }                                // Same as on phase
    virtual double      CalculatePerturbationMuOnPhase() const                                                          { return CalculatePerturbationMu(); }
    double       CalculatePerturbationQ(const double p_Radius, const double p_Rc) const                          { return log(p_Radius / p_Rc); }                                            // Hurley et al. 2000, eq 105 (really is natural log)
    double              CalculatePerturbationR(const double p_Mu, const double p_Mass, const double p_Radius, const double p_Rc) const;
    double              CalculatePerturbationS(const double p_Mu, const double p_Mass) const;

    static  double      CalculateRadialExpansionTimescale_Static(const STELLAR_TYPE p_StellarType,
                                                                 const STELLAR_TYPE p_StellarTypePrev,
                                                                 const double       p_Radius,
                                                                 const double       p_RadiusPrev,
                                                                 const double       p_DtPrev);








    virtual std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                              { return std::make_tuple(CalculateRadiusOnPhase(), m_StellarType); }

    void                CalculateRCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_RCoefficients) const;

GNU_PURE double CalculateZAMSAngularFrequency_Hurley_Static(const double p_MZAMS, const double p_RZAMS);
GNU_PURE double CalculateZAMSRotationalVelocity_Static(double p_MZAMS);





///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

virtual double CalculateAgeAfterMassLoss() const;
virtual inline double CalculateAgeAfterMassLoss_Hurley2000() const { return m_StateHistory.CurrentState.Age(); } 

virtual double CalculateCoreMass() const;
virtual inline double CalculateCoreMass_Hurley2000() const { return m_StateHistory.CurrentState.CoreMass(); } 

GNU_CONST static double CalculateCoreMass_Hurley2000_Static(const double p_Luminosity, const DBL_VECTOR& p_GBparams); // <<<<< WHAT'S THIS FOR?????

virtual double CalculateCOCoreMass() const;
virtual inline double CalculateCOCoreMass_Hurley2000() const { return m_StateHistory.CurrentState.COCoreMass(); } 

virtual double CalculateEffectiveInitialMass() const;
virtual inline double CalculateEffectiveInitialMass_Hurley2000() const { return m_StateHistory.CurrentState.MassEffectiveInitial(); }

virtual MASS_LOSS_T CalculateMassLossRate() const;

virtual double CalculateRadius() const;
virtual inline double CalculateRadius_Hurley2000() const { return m_StateHistory.ZAMSState().Radius(); }  

virtual double CalculateTau() const;
virtual inline double CalculateTau_Hurley2000() const { return m_StateHistory.CurrentState.Tau(); }


virtual double CalculateHAbundanceCore(const double p_Tau, const double p_InitialHAbundance) const { return m_StateHistory.CurrentState.HAbundanceCore(); }
virtual double CalculateHAbundanceSurface(const double p_Tau, const double p_InitialHAbundance) const { return m_StateHistory.CurrentState.HAbundanceSurface(); }
virtual double CalculateHeAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance = 0.0) const { return m_StateHistory.CurrentState.HeAbundanceCore(); }
virtual double CalculateHeAbundanceSurface(const double p_Metallicity, const double p_Tau, const double p_InitialHAbundance) const { return m_StateHistory.CurrentState.HeAbundanceSurface(); }


GNU_CONST double CalculateLuminosityGivenCoreMass_Hurley2000(const double p_CoreMass, const DBL_VECTOR& p_GBparams) const;










///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

virtual double CalculateCoreMassAtPhaseEnd() const;
virtual inline double CalculateCoreMassAtPhaseEnd_Hurley2000() const { return m_StateHistory.CurrentState.CoreMass(); }


virtual double CalculateCOCoreMassAtPhaseEnd() const;
virtual inline double CalculateCOCoreMassAtPhaseEnd_Hurley2000() const { m_StateHistory.CurrentState.COCoreMass(); }








    virtual double      CalculateTemperatureAtPhaseEnd() const                                                          { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }
    virtual double      CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const          { return CalculateTemperatureOnPhase(p_Luminosity, p_Radius); }             // Same as on phase
    double              CalculateTemperatureKelvinOnPhase(const double p_Luminosity, const double p_Radius) const;


////    virtual double      CalculateTemperatureOnPhase() const                                                             { return CalculateTemperatureOnPhase(m_Luminosity, m_Radius); }



    virtual double      CalculateTemperatureOnPhase(const double p_Luminosity, const double p_Radius) const;



    virtual void        CalculateTimescales()                                                                           { CalculateTimescales(m_Mass0, m_Timescales); }                             // Use class member variables
    virtual void        CalculateTimescales(const double p_Mass, DBL_VECTOR &p_tScales) { }                                                                                                      // Default is NO-OP


    virtual double      ChooseTimestep(const double p_Time) const                                                       { return m_dt; }

    double              DrawKickMagnitudeBrayEldridge(const double p_EjectaMass, const double p_RemnantMass, const double p_Alpha, const double p_Beta) const;
    double              DrawKickMagnitudeDistributionFlat(const double p_MaxVK, const double p_Rand) const;
    double              DrawKickMagnitudeDistributionMaxwell(const double p_Sigma, const double p_Rand) const;
    double              DrawRemnantKickMuller(const double p_COCoreMass) const;
    double              DrawRemnantKickMullerMandel(const double p_COCoreMass, const double p_Rand, const double p_RemnantMass) const;
    double              DrawSNKickMagnitude(const double p_Sigma, const double p_COCoreMass, const double p_Rand, const double p_EjectaMass, const double p_RemnantMass) const;

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

};






///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 RADIUS FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRadiusOnPhase
 *
 * @brief
 * Calculate the radius of the star at the current evolutionary phase.
 *
 * Calls relevant radius function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateRadiusOnPhase()
 *
 * @return                                      Radius of the star (Rsol)
 */
GNU_CONST inline double BaseStar::CalculateRadiusOnPhase() const { 

    double radius;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            radius = ();
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
//                            LIFETIME / AGE CALCULATIONS                            //
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
inline double BASESTAR::CalculateAgeAfterMassLoss() const {

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
GNU_CONST inline double BaseStar::CalculateTau() const { 

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
 * CalculateLifetimeToBAGB_Hurley2000
 *
 * @brief
 * Calculate the lifetime to the Base of the Asymptotic Giant Branch, per Hurley et al. 2000,
 * just before eq 69
 * 
 * tBAGB = tHeI + tHe
 *
 *
 * double CalculateLifetimeToBAGB_Hurley2000(const double p_tHeI, const double p_tHe)
 *
 * @param       p_tHeI                          Time to helium ignition (Myr)
 * @param       p_tHe                           Time to helium burning (Myr)
 * @return                                      Lifetime to the Base of the Asymptotic Giant Branch (Myr)
 */
inline GNU_CONST double BaseStar::CalculateLifetimeToBAGB_Hurley2000(const double p_tHeI, const double p_tHe) const {
    return p_tHeI + p_tHe;
}


/*
 * CalculateLifetimeToBGB_Hurley2000
 *
 * @brief
 * Calculate the lifetime to the base of the Giant Branch (end of the Hertzsprung Gap),
 * per Hurley at al. 2000, eq 4 (plotted in Hurley et al. 2000, fig 5)
 * 
 * For high mass stars, t_BGB = t_HeI.
 *
 *
 * double CalculateLifetimeToBGB_Hurley2000(const double p_Mass, const DBL_VECTOR& p_aN) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_aN                            Hurley a(n) coefficients
 * @return                                      Lifetime to the Base of the Giant Branch (Myr)
 */
inline double BaseStar::CalculateLifetimeToBGB_Hurley2000(const double p_Mass, const DBL_VECTOR& p_aN) const {
    // pow() is slow - use multiplication (sqrt() is much faster than pow())
    const double m2 = p_Mass * p_Mass;
    const double m4 = m2 * m2;
    const double m7 = p_Mass * m2 * m4;

    return (p_aN[1] + (p_aN[2] * m4) + (p_aN[3] * std::sqrt(p_Mass) * p_Mass * m4) + m7) / ((p_aN[4] * m2) + (p_aN[5] * m7));
}


/*
 * Calculate dynamical timescale
 *
 * Kalogera & Webbink 1996, eq 1
 *
 *
 * double CalculateDynamicalTimescale_Static(const double p_Mass, const double p_Radius)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @return                                      Dynamical timescale in Myr
 */
double BaseStar::CalculateDynamicalTimescale_Static(const double p_Mass, const double p_Radius) {
    return 5.0 * 1.0E-5 * p_Radius * std::sqrt(p_Radius) * YEAR_TO_MYR / std::sqrt(p_Mass);   // sqrt() is much faster than pow()
}


/*
 * Calculate thermal timescale
 *
 * pre-factor from Kalogera & Webbink 1996 (https://arxiv.org/abs/astro-ph/9508072), equation 2, 
 * combined with p_Mass * p_EnvMass case from equation 61 from https://arxiv.org/abs/astro-ph/0201220 for k in {2,3,4,5,6,8,9}
 * [note that equation 61 of BSE (https://arxiv.org/abs/astro-ph/0201220) approximates this with a value a factor of 3 smaller]
 * 
 * 
 * double CalculateThermalTimescale(const double p_Radius) const
 *
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @return                                      Thermal timescale in Myr
 *
 * The p_Radius parameter is to accommodate the call (of this function) in BaseBinaryStar::CalculateMassTransfer()
*/
double BaseStar::CalculateThermalTimescale(const double p_Radius) const {   
    return 31.4 * m_Mass * (m_Mass == m_CoreMass ? m_Mass : m_Mass - m_CoreMass) / (p_Radius * m_Luminosity); // G*Msol^2/(Lsol*Rsol) ~ 31.4 Myr (~ 30 Myr in Kalogera & Webbink)
}


/*
 * Calculate radial expansion timescale
 *
 *
 * double CalculateRadialExpansionTimescale_Static(const STELLAR_TYPE p_StellarType,
 *                                                 const STELLAR_TYPE p_StellarTypePrev,
 *                                                 const double       p_Radius,
 *                                                 const double       p_RadiusPrev,
 *                                                 const double       p_dtPrev)
 *
 * @param   [IN]    p_StellarType               Current stellar type of star
 * @param   [IN]    p_StellarTypePrev           Previous stellar type of star
 * @param   [IN]    p_Radius                    Current radius of star in Rsol
 * @param   [IN]    p_RadiusPrev                Previous radius of star in Rsol
 * @param   [IN]    p_dtPrev                    Previous timestep in Myr
 * @return                                      Radial expansion timescale in Myr
 *                                              Returns -1.0 if radial expansion timescale can't be calculated
 *                                              (i.e. stellar type has changed or radius has not changed)
 */
double BaseStar::CalculateRadialExpansionTimescale_Static(const STELLAR_TYPE p_StellarType,
                                                          const STELLAR_TYPE p_StellarTypePrev,
                                                          const double       p_Radius,
                                                          const double       p_RadiusPrev,
                                                          const double       p_dtPrev) {
    return (p_StellarTypePrev == p_StellarType && p_RadiusPrev != p_Radius) ? (p_dtPrev * p_RadiusPrev) / fabs(p_Radius - p_RadiusPrev) : -1.0;
}


/*
 * Calculate the radial expansion timescale in the mass transfer regime
 * We do not use CalculateRadialExpansionTimescale(), since in the process of mass transfer the previous radius
 * is determined by binary evolution, not nuclear timescale evolution
 *
 *
 * double CalculateRadialExpansionTimescaleDuringMassTransfer()
 *
 * @return                                      Radial expansion timescale
 */
double BaseStar::CalculateRadialExpansionTimescaleDuringMassTransfer() {
    
    // We create and age it slightly to determine how the radius will change.
    // To be sure the clone does not participate in logging, we set its persistence to EPHEMERAL.
    BaseStar *clone = Clone(OBJECT_PERSISTENCE::EPHEMERAL, false);                              // do not re-initialise the clone

    double timestep = std::max(1000.0 * NUCLEAR_MINIMUM_TIMESTEP, m_Age / 1.0E6);
    (void)clone->EvolveOneTimestep(0.0, 0.0, timestep, true);
    double radiusAfterAging = clone->Radius();
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone

    return timestep * m_Radius / fabs(m_Radius - radiusAfterAging);
}


/*
 * Calculate mass change timescale
 *
 *
 * double CalculateMassChangeTimescale_Static(const STELLAR_TYPE p_StellarType,
 *                                            const STELLAR_TYPE p_StellarTypePrev,
 *                                            const double       p_Mass,
 *                                            const double       p_MassPrev,
 *                                            const double       p_dtPrev)
 *
 * @param   [IN]    p_StellarType               Current stellar type of star
 * @param   [IN]    p_StellarTypePrev           Previous stellar type of star
 * @param   [IN]    p_Mass                      Current mass of star in Msol
 * @param   [IN]    p_MassPrev                  Previous radius of star in Msol
 * @param   [IN]    p_dtPrev                    Previous timestep in Myr
 * @return                                      Mass change timescale in Myr
 *                                              Returns -1.0 if mass change timescale can't be calculated
 *                                              (i.e. stellar type has changed or mass has not changed)
 */
double BaseStar::CalculateMassChangeTimescale_Static(const STELLAR_TYPE p_StellarType,
                                                     const STELLAR_TYPE p_StellarTypePrev,
                                                     const double       p_Mass,
                                                     const double       p_MassPrev,
                                                     const double       p_dtPrev) {

    return p_StellarTypePrev == p_StellarType && utils::Compare(p_MassPrev, p_Mass) != 0
            ? (p_dtPrev * p_MassPrev) / fabs(p_Mass - p_MassPrev)
            : -1.0;
}



/*
 * Calculate the eddy turnover timescale
 * Hurley+2002, sec. 2.3, particularly eq. 31 of subsec. 2.3.1
 *
 *
 * double CalculateEddyTurnoverTimescale()
 *
 * @return                                      eddy turnover timescale (yr)
 */
double BaseStar::CalculateEddyTurnoverTimescale() const {

	double rEnv	= CalculateRadialExtentConvectiveEnvelope();
    double mEnv, mEnvmax;
    std::tie(mEnv, mEnvmax) = CalculateConvectiveEnvelopeMass();
    return 0.4311 * cbrt((mEnv * rEnv * (m_Radius - (0.5 * rEnv))) / (3.0 * m_Luminosity));
}




/*
 * CalculateHurleyAlpha1
 *
 * @brief
 * Calculate the constant alpha1, per Hurley et al, 2000, just after eq 49
 *
 * Alpha1 depends on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 *
 * double CalculateHurleyAlpha1(const double p_MHeF, const DBL_VECTOR& p_bN) const
 *
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      Hurley constant alpha1
 */
inline GNU_CONST double BaseStar::CalculateHurleyAlpha1(const double p_MHeF, const DBL_VECTOR& p_bN) const {
    const double LHeI_MHeF = (p_bN[11] + (p_bN[12] * PPOW(p_MHeF, 3.8))) / (p_bN[13] + (p_MHeF * p_MHeF));
    return ((p_bN[9] * PPOW(p_MHeF, p_bN[10])) - LHeI_MHeF) / LHeI_MHeF;
}


/*
 * CalculateHurleyAlpha3
 *
 * @brief
 * Calculate the constant alpha3, per Hurley et al. 2000, just after eq 56
 *
 * Alpha3 depends on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 *
 * double CalculateHurleyAlpha3(const double p_MHeF, const DBL_VECTOR& p_bN) const
 * 
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      Hurley constant alpha3
 */
inline GNU_CONST double BaseStar::CalculateHurleyAlpha3(const double p_MHeF, const DBL_VECTOR& p_bN) const {
    const double LBAGB = (p_bN[31] + (p_bN[32] * PPOW(p_MHeF, (p_bN[33] + 1.8)))) / (p_bN[34] + PPOW(p_MHeF, p_bN[33]));
    return ((p_bN[29] * PPOW(p_MHeF, p_bN[30])) - LBAGB) / LBAGB;
}


/*
 * CalculateHurleyAlpha4
 *
 * @brief
 * Calculate the constant alpha4, per Hurley et al. 2000, just after eq 57
 *
 * Alpha4 depends on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 *
 * double CalculateHurleyAlpha4(const double p_MHeF, const DBL_VECTOR& p_aN, const DBL_VECTOR& p_bN) const
 *
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_aN                            Hurley a(n) coefficients
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      Hurley constant alpha4
 */
inline GNU_CONST double BaseStar::CalculateHurleyAlpha4(const double p_MHeF, const DBL_VECTOR& p_aN, const DBL_VECTOR& p_bN) const {
    const double MHeF5 = p_MHeF * p_MHeF * p_MHeF * p_MHeF * p_MHeF;    // pow() is slow - use multiplication   
    return (((CalculateLifetimeToBGB_Hurley(p_MHeF, p_aN) * (p_bN[41] * PPOW(p_MHeF, p_bN[42]) + p_bN[43] * MHeF5) / (p_bN[44] + MHeF5)) - p_bN[39]) / p_bN[39]);
}


/*
 * CalculatePerturbationS
 *
 * Calculate the perturbation parameter s
 *
 * Hurley et al. 2000, eq 101
 *
 *
 * double CalculatePerturbationS(const double p_Mass) const
 *
 * @param   [IN]    p_Mu                        Perturbation parameter mu
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Perturbation parameter s
 */
double BaseStar::CalculatePerturbationS(const double p_Mu, const double p_Mass) const {

    double b      = CalculatePerturbationB(p_Mass);
    double b_3    = b * b * b;                      // pow() is slow - use multiplication
    double mu_b_3 = p_Mu * p_Mu * p_Mu / b_3;       // calculate once, use many times...

    return ((1.0 + b_3) * mu_b_3) / (1.0 + mu_b_3);
}


/*
 * CalculatePerturbationR
 *
 * Calculate the perturbation parameter r
 *
 * Hurley et al. 2000, eq 102
 *
 *
 * double CalculatePerturbationR(const double p_Mu, const double p_Mass, const double p_Radius, const double p_Rc) const
 *
 * @param   [IN]    p_Mu                        Perturbation parameter mu
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Radius                    Radius in Rsol
 * @param   [IN]    p_Rc                        Radius that the remnant would have if the star immediately lost its envelope (in Rsol)
 * @return                                      Perturbation parameter r
 */
inline double BaseStar::CalculatePerturbationR(const double p_Mu, const double p_Mass, const double p_Radius, const double p_Rc) const {

    double r = 0.0;

    if (utils::Compare(p_Mu, 0.0) > 0 && utils::Compare(p_Radius, p_Rc) > 0) {  // only if mu > 0 and radius is larger than core radius, otherwise r = 0 and perturbed radius = core radius

        double c        = CalculatePerturbationC(p_Mass);
        double c_3      = c * c * c;                                            // pow() is slow - use multiplication
        double mu_c_3   = p_Mu * p_Mu * p_Mu / c_3;                             // calculate once

        double q        = CalculatePerturbationQ(p_Radius, p_Rc);
        double exponent = min((0.1 / q), (-14.0 / log10(p_Mu)));                // Hurley et al. 2000 is just 0.1 / q, but the Hurley sse code does this (`rpertf()` in `zfuncs.f`) - no explanation.

        r = ((1.0 + c_3) * mu_c_3 * PPOW((p_Mu), exponent)) / ((1.0 + mu_c_3));
    }

    return r;
}


/*
 * CalculateLuminosityGivenCoreMass_Hurley2000
 *
 * @brief
 * Calculate luminosity given the core mass of the star, per Hurley et al. 2000, eq 37
 * (core mass - luminosity relationship)
 *
 *
 * double CalculateLuminosityGivenCoreMass_Hurley2000(const double p_CoreMass, const DBL_VECTOR& p_GBparams) const
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Luminosity at the Base of the Asymptotic Giant Branch (Lsol)
 */
GNU_CONST inline double BaseStar::CalculateLuminosityGivenCoreMass_Hurley2000(const double p_CoreMass, const DBL_VECTOR& p_GBparams) const {
#define GBparams(x) p_GBparams[static_cast<int>(GBP::x)] // for convenience and readability - undefined at end of function
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
 * double CalculateLuminosityAtBAGB_Hurley2000(double p_Mass, double p_MHeF, const double p_Alpha3, const DBL_VECTOR& p_bN) const
 *
 * @param       p_Mass                          Effective initial mass (Msol)
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_Alpha3                        Hurley alpha3 constant
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      BAGB luminosity (Lsol)
 */
GNU_CONST inline double BaseStar::CalculateLuminosityAtBAGB_Hurley2000(double p_Mass, double p_MHeF, const double p_Alpha3, const DBL_VECTOR& p_bN) const {
    return p_Mass < p_MHeF
            ? (p_bN[29] * PPOW(p_Mass, p_bN[30])) / (1.0 + (p_Alpha3 * exp(15.0 * (p_Mass - p_MHeF))))
            : (p_bN[31] + (p_bN[32] * PPOW(p_Mass, (p_bN[33] + 1.8)))) / (p_bN[34] + PPOW(p_Mass, p_bN[33]));
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
GNU_CONST inline double BaseStar::CalculateTemperatureOnPhase() const { 

    double temperature;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            // retrieve luminosity and radius from current state and call relevant function
            temperature = CalculateTemperatureOnPhase(m_StateHistory.CurrentState.Luminosity(), m_StateHistory.CurrentState.Radius());
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
//                                  MASS FUNCTIONS                                   //
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
GNU_CONST inline double BaseStar::CalculateCoreMass() const { 

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
GNU_CONST inline double BaseStar::CalculateCoreMassAtPhaseEnd() const { 

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
    return p_Luminosity > p_GBparams[static_cast<int>(GBP::Lx)]
            ? PPOW((p_Luminosity / p_GBparams[static_cast<int>(GBP::B)]), (1.0 / p_GBparams[static_cast<int>(GBP::q)]))
            : PPOW((p_Luminosity / p_GBparams[static_cast<int>(GBP::D)]), (1.0 / p_GBparams[static_cast<int>(GBP::p)]));
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
GNU_CONST inline double BaseStar::CalculateCOCoreMass() const { 

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
GNU_CONST inline double BaseStar::CalculateCOCoreMassAtPhaseEnd() const { 

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
 * CalculateEffectiveInitialMass
 *
 * @brief
 * Calculate the effective initial mass of the star.
 *
 * Calls relevant effective initial mass function based on the evolutionary mode given in program options.
 * 
 * The effective initial mass of the star is the variable M0 in Hurley et al. 2000.  M0 is introduced in
 * section 7 of Hurley et al. 2000 and is described there as the "initial mass". M0 is the phase-specific
 * initial mass - i.e., it is the initial mass for the current evolutionary phase of the star, but can 
 * track Mt (mass at time t), depending upon the phase, so (as stated in Hurley et al. 2000), it is really
 * the "effective initial mass".
 * 
 *
 * double CalculateEffectiveInitialMass() const
 *
 * @return                                      Effective initial mass (Msol)
 */
GNU_CONST inline double BaseStar::CalculateEffectiveInitialMass() const { 

    double mass;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            mass = CalculateEffectiveInitialMass_Hurley2000();
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

    return mass;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                             MASS LOSS RATE FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////
    

/*
 * CalculateMLRate_Hurley2000
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, per Hurley et al. 2000.
 * 
 * 
 * MASS_LOSS_T CalculateMLRate_Hurley2000(const double p_Mass,
 *                                        const double p_Radius,
 *                                        const double p_Luminosity,
 *                                        const double p_PerturbationMu,
 *                                        const double p_ZscaledHurley,
 *                                        const double p_WRfactor) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @param       p_ZscaledHurley                 Z inversely scaled by Hurley ZSOL (Z / ZSOL_HURLEY)
 * @param       p_WRfactor                      WR mass loss factor
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (may be MASS_LOSS_TYPE::NONE)
 */
GNU_CONST inline MASS_LOSS_T CalculateMLRate_Hurley2000(const double p_Mass,
                                                        const double p_Radius,
                                                        const double p_Luminosity,
                                                        const double p_PerturbationMu,
                                                        const double p_ZscaledHurley,
                                                        const double p_WRfactor) const {
    return CalculateMLRate_NieuwenhuijzenDeJager1990(p_Mass, p_Radius, p_Luminosity, p_ZscaledHurley);
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
 * MASS_LOSS_T CalculateMLRate_KudritzkiReimers1978(const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::GB)
 * @return                                      Kudritzki and Reimers mass loss rate (in Msol yr^{-1})
 */
GNU_CONST MASS_LOSS_T BaseStar::CalculateMLRate_KudritzkiReimers1978(const double p_Mass, const double p_Radius, const double p_Luminosity) const {
    // Hurley et al. 2000 has eta^2 - that's wrong per Hurley SSE code
    return std::make_tuple(4.0E-13 * (MASS_LOSS_ETA * p_Luminosity * p_Radius / p_Mass), MASS_LOSS_TYPE::GB);
}


/* 
 * CalculateMLRate_VassiliadisWood1993
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for stars on the AGB,
 *  based on the Mira pulsation period (P0), per Vassiliadis and Wood 1993.
 *
 * See Hurley et al. 2000, just after eq 106
 * Note: in the Hurley fortran code, P0 is taken to be min(p0, 2000.0) - implemented here as a minimum power
 *
 *
 * MASS_LOSS_T CalculateMLRate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate on the AGB (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::GB)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const {
    const double P0   = PPOW(10.0, std::min(3.3, (-2.07 - (0.9 * log10(p_Mass)) + (1.94 * log10(p_Radius)))));
    const double dMdt = PPOW(10.0, (-11.4 + (0.0125 * (P0 - 100.0 * std::max((p_Mass - 2.5), 0.0)))));

    return std::make_tuple(std::min(dMdt, (1.36E-9 * p_Luminosity)), MASS_LOSS_TYPE::GB);
}


/*
 * CalculateMLRateLBV_Belczynski2010
 *
 * @brief
 * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit
 * (Humphreys & Davidson 1994), per Belczynski et al. 2010, eq 8 
 * 
 *
 * MASS_LOSS_T CalculateMLRateLBV_Belczynski2010(const double p_LBVfactor) const
 *
 * @param       p_LBVfactor                     Luminous Blue Variable factor
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for LBV stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::LBV)
 */    
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateLBV_Belczynski2010(const double p_LBVfactor) const {
    return std::make_tuple(p_LBVfactor * 1.0E-4, MASS_LOSS_TYPE::LBV);
} 


/*
 * CalculateMLRateLBV_Hurley2000
 *
 * @brief
 * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit
 * (Humphreys & Davidson 1994), per Hurley et al. 2000, sec 7.1, unlabelled equation
 * a few equations after eq 106
 *  
 *
 * MASS_LOSS_T CalculateMLRateLBV_Hurley2000(const double p_Luminosity, const double p_HDlimitfactor) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_HDlimitfactor                 Factor by which star is above Humphreys-Davidson limit
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for LBV stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::LBV)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateLBV_Hurley2000(const double p_Luminosity, const double p_HDlimitfactor) const {
    return std::make_tuple(0.1 * utils::IntPow(p_HDlimitfactor - 1.0, 3) * ((p_Luminosity / 6.0E5) - 1.0), MASS_LOSS_TYPE::LBV);
}


/*
 * CalculateMLRateOB_Krticka2018
 *
 * Calculate mass loss rate for massive OB stars, per Krticka et al 2018
 * https://arxiv.org/pdf/1712.03321.pdf
 *
 * Uses current values of:
 * 
 *    - m_Luminosity
 *
 * 
 * MASS_LOSS_T CalculateMLRateOB_Krticka2018(const double p_Luminosity, const double p_ZetaAsplund) const
 * 
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_ZetaAsplund                   Asplund zeta value (log10(Z / ZSOL_ASPLUND))
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for hot OB stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::OB)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateOB_Krticka2018(const double p_Luminosity, const double p_ZetaAsplund) const {
    return std::make_tuple(PPOW(10.0, -5.70 + 0.50 * p_ZetaAsplund() + (1.61 - 0.12 * p_ZetaAsplund()) * log10(p_Luminosity / 1.0E6)), MASS_LOSS_TYPE::OB);
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
 * MASS_LOSS_T CalculateMLRateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const
 *
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const {
    return std::make_tuple(PPOW(10.0, -21.5 - 0.15 * p_mStart + 3.6 * log10(p_Luminosity)), MASS_LOSS_TYPE::RSG);
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
 * MASS_LOSS_T CalculateMLRateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const
 *
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const {
    return std::make_tuple(PPOW(10.0, -20.63 - 0.16 * p_mStart + 3.47 * log10(p_Luminosity)), MASS_LOSS_TYPE::RSG);
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
 * MASS_LOSS_T CalculateMLRateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const {
    const double exp = logL < 4.6 ? -8.0 + 0.7 * log10(p_Luminosity) - 0.7 * log10(p_Mass) : -24.0 + 4.77 * log10(p_Luminosity) - 3.99 * log10(p_Mass);
    return std::make_tuple(PPOW(10.0, exp), MASS_LOSS_TYPE::RSG);
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
 * GNU_CONST MASS_LOSS_T CalculateMLRateRSG_Yang2023(const double p_Luminosity) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateRSG_Yang2023(const double p_Luminosity) const {

    const double logL  = log10(p_Luminosity);
    const double logL2 = logL * logL;

    return std::make_tuple(PPOW(10.0, 0.45 * logL2 * logL - 5.26 * logL2 + 20.93 * logL - 34.56), MASS_LOSS_TYPE::RSG);
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
 * MASS_LOSS_T CalculateMLRateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for very massive stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::VMS)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const {

    constexpr double alpha       = 0.39;                                        // CAK force multiplier
    constexpr double logMdotZero = -4.78;                                       // from substituting LogMdotTrans and Gamma_e trans into eq 12. 

    const double gamma = EDDINGTON_PARAMETER_FACTOR * p_Luminosity / p_Mass;    // Eddington Parameter, not metallicity specific as in the publication

    return std::make_tuple(PPOW(10.0, logMdotZero + ((1.0 / alpha) + 0.5) * log10(gamma) - (((1.0 - alpha) / alpha) + 2.0) * log10(1.0 - gamma)), MASS_LOSS_TYPE::VMS);
}


/*
 * CalculateMLRateWR_Hurley2000
 *
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
 * MASS_LOSS_T CalculateMLRateWR_Hurley2000(const double p_Luminosity, const double p_PerturbationMu) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateWR_Hurley2000(const double p_Luminosity, const double p_PerturbationMu) const {
    return std::make_tuple((p_PerturbationMu >= 1.0 ? 0.0 : PPOW(p_Luminosity, 1.5) * (1.0 - p_PerturbationMu) * 1.0E-13), MASS_LOSS_TYPE::WR);
} 


/*
 * CalculateMLRateWR_Shenar2019_Static
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
 * MASS_LOSS_T CalculateMLRateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature, const double p_SigmaHurley) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_SigmaHurley                   Sigma from Hurley et al. 2000 p24, sigma = log10(Z)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature, const double p_SigmaHurley) const {

    // For H-rich WR stars (X_H > 0.4)
    constexpr double C1 = -6.78;
    constexpr double C2 =  0.66;
    constexpr double C3 = -0.12;
    constexpr double C5 =  0.74;

    return std::make_tuple(PPOW(10.0, C1 + (C2 * log10(p_Luminosity)) + (C3 * log10(p_Temperature * TSOL)) + (C5 * p_SigmaHurley)), MASS_LOSS_TYPE::WR);
}


/*
 * CalculateMLRateWR_ZDependent_Static
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
 * MASS_LOSS_T CalculateMLRateWR_ZDependent_Static(const double p_Metallicity, const double p_Luminosity, const double p_PerturbationMu) const
 *
 * @param       p_Metallicity                   (Fractional) wetallicity of the star
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
GNU_CONST inline MASS_LOSS_T BaseStar::CalculateMLRateWR_ZDependent_Static(const double p_Metallicity, const double p_Luminosity, const double p_PerturbationMu) const {
    // StarTrack may still do something different here.
    // There are references to Hamann & Koesterke 1998 and Vink and de Koter 2005.
    // TW - Haven't seen StarTrack but I think H&K gives the original equation and V&dK gives the Z dependence
    const double dMdt = p_PerturbationMu >= 1.0 ? 0.0 : 1.0E-13 * p_Luminosity * std::sqrt(p_Luminosity) * PPOW(p_Metallicity / ZSOL_ANDERS, 0.86) * (1.0 - p_PerturbationMu);
    return std::make_tuple(dMdt, MASS_LOSS_TYPE::WR);
}    











/*
 * CalculateConvectiveEnvelopeBindingEnergy
 *
 * @brief
 * Calculate the convective envelope binding energy for the two-stage common envelope formalism of Hirai & Mandel, 2022
 *
 *
 * double CalculateConvectiveEnvelopeBindingEnergy(const double p_Mass, const double p_EnvMass, const double p_Radius, const double p_Lambda) const
 *
 * @param           p_Mass                      Mass (total) of the star (Msol)
 * @param           p_EnvMass                   Mass of the convective outer envelope  (Msol)
 * @param           p_Radius                    Radius of the star (Rsol)
 * @param           p_Lambda                    Lambda parameter for the convective envelope
 * @return                                      Binding energy (erg)
 */
inline double BaseStar::CalculateConvectiveEnvelopeBindingEnergy(const double p_Mass, const double p_EnvMass, const double p_Radius, const double p_Lambda) const {
    return CalculateBindingEnergy(p_Mass - p_EnvMass, p_EnvMass, p_Radius, p_Lambda);
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




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 RADIUS FUNCTIONS                                  //
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
GNU_CONST inline double BaseStar::CalculateRadius() const { 

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
        if (m_StateHistory.CurrentState().StellarType() == elem) return true;
    }
    return false;
};












//// <<<<<<<<<<<<<<<<<<<<<<<<<<<<< constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


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
    return CalculateCriticalMassRatio_Ge2020_Interpolate(m_StateHistory.CurrentState().Mass(), m_StateHistory.CurrentState().Radius(), p_MTefficiency);
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




GNU_PURE  virtual double CalculateCELambda_Nanjing_StarTrack(const double p_Metallicity,
                                                             const double p_Mass,
                                                             const double p_Radius,
                                                             const double p_CoreMass) const { return 1.0; } // Default for stellar types with no LamdaNanjing definitions - 1.0 is benign




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
