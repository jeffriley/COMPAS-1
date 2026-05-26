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
    
protected:

    // Object identifiers - all classes have these
    OBJECT_ID               m_ObjectId;                                 // Instantiated object's unique object id
    OBJECT_PERSISTENCE      m_ObjectPersistence;                        // Instantiated object's persistence (permanent or ephemeral)

    StateHistory<StarState> m_StateHistory{StarState{}};                // State history stack - current state is empty
    StarState               m_InterimState;                             // Interim state of the star

    // Convenience aliases for star state.
    // Each of these is wrapped in std::optional<>, and may not have an assigned values
    // (at least initially).
    auto& CurrentState = m_StateHistory.CurrentState();                 // The most recently recorded state of the star
    auto& StartState   = m_StateHistory.StartState();                   // The first state of the star - may or may not be the ZAMS state
    auto& ZAMSState    = m_StateHistory.ZAMSState();                    // The ZAMS state of the star (if known)
    auto& PrevState    = m_StateHistory.PreviousState();                // The previous state of the star

public:
    
    BaseStar::BaseStar(const UlongT             p_RandomSeed,
                       const STELLAR_TYPE       p_StellarType, 
                       const double             p_Metallicity,
                       const double             p_Mass,
                       const StellarKickParmsT& p_KickParameters,
                       const double             p_RotationalFrequency = -1.0) {

                        
        // Initialise as much of the first state of the star as possible.
        // The state history stack has already been initialised - the single
        // entry is the first state of the star, but it is empty.


        // Set state attributes supplied as parameters
        CurrentState.SetRandomSeed(p_RandomSeed);
        CurrentState.SetStellarType(p_StellarType);
        CurrentState.SetMetallicity(p_Metallicity);
        CurrentState.SetMass(p_Mass);
        CurrentState.SetMassEffectiveInitial(p_Mass);
        CurrentState.SetSuppliedKickParameters(p_KickParameters);


        // Set, or calculate where possible/necessary, remaining state attributes.
        // Any state attribute for which we either don't know, or can't calculate, the
        // correct initial value, we leave uninitialised - the state getters will return
        // a value wrapped in std::optional<> set to std::nullopt (no value assigned).

        CurrentState.SetError(ERROR::NONE);
        CurrentState.SetEvolutionStatus(EVOLUTION_STATUS::CONTINUE);

        CurrentState.SetLuminosity(CalculateLuminosityAtZAMS_Tout1996(p_Mass));
        CurrentState.LuminosityZAMSEffective(state.Luminosity());

        CurrentState.SetRadius(CalculateRadiusAtZAMS_Tout1996(p_Mass));
        CurrentState.SetRadiusZAMSEffective(state.Radius());

        if (p_RotationalFrequency < 0.0) {                                                      // Rotational frequency passed in?
            if (utils::isOneOf(p_StellarType, MAIN_SEQUENCE)) {                                 // Starting on Main Sequence (can't be CHE here)?
                CurrentState.SetOmegaZAMS(utils::CalculateZAMSOmega(p_Mass, state.Radius()));   // Yes - calculate Omega at ZAMS and set state
            }                                                                                   // N0 - leave state value uninitialised
        }
        else CurrentState.SetOmegaZAMS(p_RotationalFrequency * _2_PI);                          // No - use the rotational frequency value passed in



        m_AngularMomentum                          = CalculateMomentOfInertiaAU() * m_OmegaZAMS;    // Angular Momentum (Msol * AU^2 / yr)




        CurrentState.SetTemperature(CalculateTemperature(state.Luminosity(), state.Radius()));

        CurrentState.SetAge(0.0);
        CurrentState.SetTau(0.0);
        CurrentState.SetTime(0.0);

        CurrentState.CNOProcessedCoreMass(0.0);    // Default for MS: changed in MS_gt_07 constructor
        CurrentState.COCoreMass(0.0);              // Default for MS
        CurrentState.CoreMass(0.0);                // Default for MS: changed in MS_gt_07 constructor
        CurrentState.HeCoreMass(0.0);              // Default for MS

        CurrentState.SetCHE(false);                // Initially
    }
        
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
    // declarations.  See <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS

    
    // Object identifiers - all classes have these
                OBJECT_ID           ObjectId() const                            { return m_ObjectId; }
                OBJECT_PERSISTENCE  ObjectPersistence() const                   { return m_ObjectPersistence; }
    GNU_CONST   OBJECT_TYPE         ObjectType() const                          { return OBJECT_TYPE::BASE_STAR; }


    // Getters - alphabetically
    // Getter implementations should generally be very short, so many will be implemented here.

                double              Age() const                                 { return m_StateHistory.CurrentState.Age(); }
                double              AngularFrequency() const                    { return m_StateHistory.CurrentState.AngularFrequency(); }
                bool                CHE() const                                 { return m_StateHistory.CurrentState.CHE(); }
                double              COCoreMass() const                          { return m_StateHistory.CurrentState.COCoreMass(); }
                double              CoreMass() const                            { return m_StateHistory.CurrentState.CoreMass(); }
                MASS_LOSS_TYPE      DominantMassLossType() const                { return m_StateHistory.CurrentState.DominantMassLossType(); }
                double              dt() const                                  { return m_StateHistory.CurrentState.dt(); }
                double              dtPrev() const                              { return m_StateHistory.PreviousState.dt(); }
                ERROR               Error() const                               { return m_StateHistory.Error(); }
                EVOLUTION_STATUS    EvolutionStatus() const                     { return m_StateHistory.CurrentState.EvolutionStatus(); }
                bool                ExperiencedAIC() const                      { return utils::SNEventType(m_StateHistory.CurrentState.SupernovaDetails().events.past) == SN_EVENT::AIC; }
                bool                ExperiencedCCSN() const                     { return utils::SNEventType(m_StateHistory.CurrentState.SupernovaDetails().events.past) == SN_EVENT::CCSN; }
                bool                ExperiencedHeSD() const                     { return utils::SNEventType(m_StateHistory.CurrentState.SupernovaDetails().events.past) == SN_EVENT::HeSD; }
                bool                ExperiencedECSN() const                     { return utils::SNEventType(m_StateHistory.CurrentState.SupernovaDetails().events.past) == SN_EVENT::ECSN; }
                bool                ExperiencedPISN() const                     { return utils::SNEventType(m_StateHistory.CurrentState.SupernovaDetails().events.past) == SN_EVENT::PISN; }
                bool                ExperiencedPPISN() const                    { return utils::SNEventType(m_StateHistory.CurrentState.SupernovaDetails().events.past) == SN_EVENT::PPISN; }
                bool                ExperiencedSNIA() const                     { return utils::SNEventType(m_StateHistory.CurrentState.SupernovaDetails().events.past) == SN_EVENT::SNIA; }
                SN_EVENT            ExperiencedSN_Type() const                  { return utils::SNEventType(m_StateHistory.CurrentState.SupernovaDetails().events.past); }
                bool                ExperiencedUSSN() const                     { return (m_StateHistory.CurrentState.SupernovaDetails().events.past & SN_EVENT::USSN) == SN_EVENT::USSN; }
                double              GBparam(HURLEY_GBP p_GBparam) const         { return m_StateHistory.CurrentState.GBparams(p_GBparam); }
                DblVectorT          GBparams() const                            { return m_StateHistory.CurrentState.GBparams(); }
                double              HAbundanceCore() const                      { return m_StateHistory.CurrentState.HAbundanceCore(); }
                double              HAbundanceSurface() const                   { return m_StateHistory.CurrentState.HAbundanceSurface(); }
                double              HeAbundanceCore() const                     { return m_StateHistory.CurrentState.HeAbundanceCore(); }
                double              HeAbundanceSurface() const                  { return m_StateHistory.CurrentState.HeAbundanceSurface();}
                double              HeCoreMass() const                          { return m_StateHistory.CurrentState.HeCoreMass(); }
                bool                IsAIC() const                               { return (m_StateHistory.CurrentState.SupernovaDetails().events.current & SN_EVENT::AIC) == SN_EVENT::AIC; }
                bool                IsCCSN() const                              { return (m_StateHistory.CurrentState.SupernovaDetails().events.current & SN_EVENT::CCSN) == SN_EVENT::CCSN; }
        virtual bool                IsDegenerate() const                        { return false; }   // virtual - default is not degenerate - White Dwarfs, NS and BH are degenerate
                bool                IsECSN() const                              { return (m_StateHistory.CurrentState.SupernovaDetails().events.current & SN_EVENT::ECSN) == SN_EVENT::ECSN; }
                bool                IsHeSD() const                              { return (m_StateHistory.CurrentState.SupernovaDetails().events.current & SN_EVENT::HeSD) == SN_EVENT::HeSD; }
                bool                IsPISN() const                              { return (m_StateHistory.CurrentState.SupernovaDetails().events.current & SN_EVENT::PISN) == SN_EVENT::PISN; }
                bool                IsPPISN() const                             { return (m_StateHistory.CurrentState.SupernovaDetails().events.current & SN_EVENT::PPISN) == SN_EVENT::PPISN; }
                bool                IsSN_NONE() const                           { return m_StateHistory.CurrentState.SupernovaDetails().events.current == SN_EVENT::NONE; }
                bool                IsSNIA() const                              { return (m_StateHistory.CurrentState.SupernovaDetails().events.current & SN_EVENT::SNIA) == SN_EVENT::SNIA; }
        virtual bool                IsSupernova() const                         { return false; }   // virtual - default value
                bool                IsUSSN() const                              { return (m_SupernovaDetails.events.current & SN_EVENT::USSN) == SN_EVENT::USSN; }
                double              Luminosity() const                          { return m_StateHistory.CurrentState.Luminosity(); }
                double              LuminosityPrev() const                      { return m_StateHistory.PreviousState.Luminosity(); }
                double              LZAMS() const                               { return m_StateHistory.ZAMSState.Luminosity(); }
                double              Mass() const                                { return m_StateHistory.CurrentState.Mass(); }
                double              Mass0() const                               { return m_StateHistory.CurrentState.MassEffectiveInitial(); }
                double              MassEffectiveInitial() const                { return m_StateHistory.CurrentState.MassEffectiveInitial(); }
                double              MassPrev() const                            { return m_StateHistory.PreviousState.Mass(); }
                double              MZAMS() const                               { return m_StateHistory.ZAMSState.Mass(); }
                double              PulsarMagneticField() const                 { return m_StateHistory.PulsarDetails().magneticField; }
                double              PulsarSpinPeriod() const                    { return m_StateHistory.PulsarDetails().spinPeriod; }
                double              PulsarSpinFrequency() const                 { return m_StateHistory.PulsarDetails().spinFrequency; }
                double              PulsarSpinDownRate() const                  { return m_StateHistory.PulsarDetails().spinDownRate; }
                double              PulsarBirthPeriod() const                   { return m_StateHistory.PulsarDetails().birthPeriod; }
                double              PulsarBirthSpinDownRate() const             { return m_StateHistory.PulsarDetails().birthSpinDownRate; }
                double              Radius() const                              { return m_StateHistory.CurrentState.Radius(); }
                double              RadiusPrev() const                          { return m_StateHistory.PreviousState.Radius(); }
                ULongT              RandomSeed() const                          { return m_StateHistory.RandomSeed(); }
                double              RZAMS() const                               { return m_StateHistory.ZAMSState.Radius(); }
                double              StellarType() const                         { return m_StateHistory.CurrentState.StellarType(); }
                double              StellarTypePrev() const                     { return m_StateHistory.Previoustate.StellarType(); }
                double              Tau() const                                 { return m_StateHistory.CurrentState.Tau(); }
                double              Temperature() const                         { return m_StateHistory.CurrentState.Temperature(); }
                double              Time() const                                { return m_StateHistory.CurrentState.Time(); }
                double              Timescale(TIMESCALE p_Timescale) const      { return m_StateHistory.CurrentState.Timescales(p_Timescale); }
                DblVectorT          Timescales() const                          { return m_StateHistory.CurrentState.Timescales(); }
                double              TZAMS() const                               { return m_StateHistory.ZAMSState.Temperature(); }










    ////double                      AngularMomentum() const                                 { return m_AngularMomentum; }
    bool                        EnvelopeJustExpelledByPulsations() const                { return m_EnvelopeJustExpelledByPulsations; }
    _STR_                       MassTransferDonorHistory() const; // implementation in cpp file

    double                      InitialHeliumAbundance() const                          { return m_InitialHeliumAbundance; }
    double                      InitialHydrogenAbundance() const                        { return m_InitialHydrogenAbundance; }
    double                      InitialMainSequenceCoreMass() const                     { return m_InitialMainSequenceCoreMass; }

    bool                        LBV_PhaseFlag() const                                   { return m_LBVphaseFlag; }
    double                      MainSequenceCoreMass() const                            { return m_MainSequenceCoreMass; }
    ST_VECTOR                   MassTransferDonorHistory() const                        { return m_MassTransferDonorHistory; }
    double                      Mdot() const                                            { return m_Mdot; }
    double                      Omega() const                                           { return m_AngularMomentum / CalculateMomentOfInertiaAU(); }
    double                      OmegaCHE() const                                        { return m_OmegaCHE; }
    double                      OmegaBreak() const                                      { return CalculateOmegaBreak(); }
    double                      OmegaZAMS() const                                       { return m_OmegaZAMS; }
    COMPAS_VARIABLE             PropertyValue(const T_ANY_PROPERTY p_Property) const; // implmentation in cpp file

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
    // Setter implementations should generally be very short, so many will be implemented here.

    inline void             SetPersistence(const OBJECT_PERSISTENCE p_Persistence)  { m_ObjectPersistence = p_Persistence; }
    inline void             SetStellarType(const STELLAR_TYPE p_StellarType) const  { m_StateHistory.CurrentState.SetStellarType(p_StellarType); }
    inline void             SetStartingType(const STELLAR_TYPE p_StellarType) const { m_StateHistory.StartState.SetStellarType(p_StellarType); }



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
    





    // Member functions (not getters or setters)

    ////////////////////////////////////////
    //   VIRTUAL FUNCTIONS                //
    ////////////////////////////////////////
     
    // By function, then alphabetical.
    // May be (are expected to be) overridden by derived classes.
    // When overriding virtual functions in a derived class, use the "override" attribute.


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

    virtual double     CalculateAgeAfterMassLoss() const;
    virtual double     CalculateAgeAfterMassLoss_Hurley2000() const { return Age(); } 

    virtual double     CalculateTau() const;
    virtual double     CalculateTau_Hurley2000() const { return Tau(); }

    virtual DblVectorT CalculateTimescales() const;
    virtual DblVectorT CalculateTimescales_Hurley2000() const { return m_InterimState.Timescales(); }

    virtual double     CalculateTimescale_Thermal() const { return CalculateTimescale_Thermal(Mass(), Radius(), Luminosity(), CoreMass()); }


    ////////////////////////////////////////
    //   ENVELOPE                         //
    ////////////////////////////////////////

    virtual DblVectorT CalculateConvectiveEnvelopeMass() const { return std::make_tuple(0.0, 0.0); }
    virtual DblVectorT CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_Tau) const override { return std::make_tuple(0.0, 0.0); }

    virtual double     CalculateConvectiveEnvelopeRadialExtent() const { return 0.0; }     // Default for stars with no convective envelope

    virtual ENVELOPE   DetermineEnvelopeType() const { return ENVELOPE::REMNANT; }         // Default is REMNANT - but should never be called


    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////

    virtual double CalculateLuminosity() const;
    virtual double CalculateLuminosityAtPhaseEnd() const;
    virtual double CalculateLuminosityAtPhaseEnd_Hurley2000() const { return Luminosity(); } 
    virtual double CalculateLuminosity_Hurley2000() const { return Luminosity(); } 


    ////////////////////////////////////////
    //   MASS                             //
    ////////////////////////////////////////

    virtual double CalculateCoreMass() const;
    virtual double CalculateCoreMass_Hurley2000() const { return CoreMass(); } 

    virtual double CalculateCOCoreMass() const;
    virtual double CalculateCOCoreMass_Hurley2000() const { return COCoreMass(); } 

    virtual double CalculateConvectiveCoreMass() const { return 0.0; }

    virtual double CalculateEffectiveInitialMass_Hurley2000() const { return Mass(); }

    virtual double CalculateHeCoreMass() const;
    virtual double CalculateHeCoreMass_Hurley2000() const { return HeCoreMass(); } 

    virtual double CalculateMass0() const;
    virtual double CalculateMass0_Hurley2000() const { return Mass0(); }


    ////////////////////////////////////////
    //   MASS LOSS                        //
    ////////////////////////////////////////

    virtual MassLossT CalculateMassLossRate() const;

    virtual MassLossT CalculateMLrate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_PerturbationMu = 0.0) const;

    virtual MassLossT CalculateMLrate_Merritt2025(
        const double p_Mass,
        const double p_Radius,
        const double p_Luminosity,
        const double p_Temperature,
        const double p_PerturbationMu,
        const double p_mStart,
        const double p_HeAbundanceSurface
    ) const;

    // mass loss PRESCRIPTIONS
    virtual MassLossT CalculateMLrate_Belczynski2010(
        const double p_Mass,
        const double p_Radius,
        const double p_Luminosity,
        const double p_Temperature,
        const double p_PerturbationMu,
        const double p_HeAbundanceSurface
    ) const;

    // mass loss THERMAL
    virtual double CalculateMLRateThermal() const { return Mass() / CalculateTimescale_Thermal(); }


    ////////////////////////////////////////
    //   MISC. CONSTANTS / PARAMETERS     //
    ////////////////////////////////////////

    virtual DblVectorT CalculateGBparams() const;
    virtual DblVectorT CalculateGBparams_Hurley2000() const { return m_InterimState.GBparams(); } // only required for GB stars


    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////

    virtual double CalculateConvectiveCoreRadius(const double p_Radius, const double p_Tau) const { return 0.0; }

    virtual double CalculateRadius() const;
    virtual double CalculateRadius_Hurley2000() const { return RZAMS(); } 

    virtual double CalculateRadiusAtPhaseEnd() const;
    virtual double CalculateRadiusAtPhaseEnd_Hurley2000() const { return RZAMS(); } 

    virtual double CalculateRadiusOnMassChange(double p_dM) { return Radius(); }   // Default is NO-OP


    virtual double CalculateRemnantRadius() const;
    virtual double CalculateRemnantRadius_Hurley2000() const { return Radius(); }  // Relevant for MS stars, over-written for GB stars


    ////////////////////////////////////////
    //   ROTATION                         //
    ////////////////////////////////////////

    // WHAT SHOULD THESE DEFAULT TO IF NOT MS START?  CHECK ALL OTHER STs
    virtual double CalculateMomentOfInertia() const { return (0.1 * (m_Mass) * m_Radius * m_Radius); } // Defaults to MS. k2 = 0.1 as defined in Hurley et al. 2000, after eq 109 . *Ilya* is this valid for Brcek MS core?
    virtual double CalculateMomentOfInertiaAU() const { return CalculateMomentOfInertia() * RSOL_TO_AU * RSOL_TO_AU; }

    ////////////////////////////////////////
    //   TEMPERATURE                      //
    ////////////////////////////////////////

    virtual double CalculateTemperature() const { return CalculateTemperature(Luminosity(), Radius()); };
    COMPAS_PURE double CalculateTemperature(const double p_Luminosity, const double p_Radius) const;



    // Non-virtual functions - by function, then alphabetical

    ////////////////////////////////////////
    //   NON-VIRTUAL FUNCTIONS            //
    ////////////////////////////////////////

    // Should not be overridden (declared separately) by derived classes.
    // While it is legal in C++ to declare the same (non-virtual) function in multiple classes,
    // (aka "shadowing", or "hiding"), we discourage it.  Non-virtual functions are statically
    // linked, and as such, especially with indirection, may not produce expected results.


    ////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES   //
    ////////////////////////////////////////

    GNU_CONST   double              CalculateLifetimeToBAGB_Hurley2000(const double p_tHeI, const double p_tHe) const;

    inline      double              CalculateTimescale_Dynamical() const { return CalculateTimescale_Dynamical_Kalogera1996(Mass(), Radius()); }

    GNU_CONST   double              CalculateTimescale_Dynamical_Kalogera1996(const double p_Mass, const double p_Radius) const;

    GNU_CONST   double              CalculateTimescale_Thermal(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const;

    inline      double              CalculateTimescale_EddyTurnover() const { return CalculateTimescale_EddyTurnover_Hurley2002(Mass(), Radius(), Luminosity(), CoreMass()); }
    COMPAS_PURE inline double       CalculateTimescale_EddyTurnover_Hurley2002(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const;


    ////////////////////////////////////////
    //   ENERGY                           //
    ////////////////////////////////////////

    inline double                   CalculateBindingEnergy(const double p_Lambda) const { return CalculateBindingEnergy(Radius(), CoreMass(), Mass() - CoreMass(), p_Lambda); }
    GNU_CONST double                CalculateBindingEnergy(const double p_Radius, const double p_CoreMass, const double p_EnvMass, const double p_Lambda) const;

    inline double                   CalculateConvectiveEnvelopeBindingEnergy_Hirai2022(const double p_Lambda) const { 
        double convectiveEnvMass;
        std::tie(convectiveEnvMass, std::ignore) = CalculateConvectiveEnvelopeMass();
        return CalculateConvectiveEnvelopeBindingEnergy_Hirai2022(Mass(), Radius(), convectiveEnvMass, p_Lambda);
    }
    GNU_CONST double                CalculateConvectiveEnvelopeBindingEnergy_Hirai2022(const double p_Mass, const double p_Radius, const double p_EnvMass, const double p_Lambda) const;


    ////////////////////////////////////////
    //   ENVELOPE                         //
    ////////////////////////////////////////

    COMPAS_PURE double              CalculateConvectiveEnvelopeLambda_Picker2024(const double p_Mass, const double p_EnvMass, const double p_EnvMassMax) const;


    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////

    COMPAS_PURE double              CalculateLuminosityAtZAMS_Tout1996(const double p_MZAMS) const;
    GNU_CONST   double              CalculateLuminosity_Hurley2000(const double p_CoreMass, const DblVectorT& p_GBparams) const;


    ////////////////////////////////////////
    //   MASS                             //
    ////////////////////////////////////////

    GNU_CONST double CalculateCoreMass_Hurley2000(const double p_Luminosity, const DblVectorT& p_GBparams);


    ////////////////////////////////////////
    //   MASS LOSS                        //
    ////////////////////////////////////////

    // mass loss LBV
    COMPAS_PURE MassLossT       CalculateMLrateLBV(const double p_Radius, const double p_Luminosity, const LBV_MASS_LOSS_PRESCRIPTION p_LBV_MLprescription) const;
    COMPAS_PURE MassLossT       CalculateMLrateLBV_Belczynski2010() const;
    GNU_CONST   MassLossT       CalculateMLrateLBV_Hurley2000(const double p_Luminosity, const double p_HDlimitfactor) const;

    // mass loss OB
    COMPAS_PURE MassLossT       CalculateMLrateOB(const double p_Mass, const double p_Luminosity, const double p_Temperature, const OB_MASS_LOSS_PRESCRIPTION p_OB_MLprescription) const;
    COMPAS_PURE MassLossT       CalculateMLrateOB_Bjorklund2022(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE MassLossT       CalculateMLrateOB_Krticka2018(const double p_Luminosity) const;
    COMPAS_PURE MassLossT       CalculateMLrateOB_Vink2001(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE MassLossT       CalculateMLrateOB_VinkSander2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;

    // mass loss PRESCRIPTIONS
    COMPAS_PURE MassLossT       CalculateMLrate_NieuwenhuijzenDeJager1990(const double p_Mass, const double p_Radius, const double p_Luminosity) const;
    GNU_CONST   MassLossT       CalculateMLrate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const;

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// removed in PR 1445
/*
    COMPAS_PURE double BaseStar::CalculateMLrate_WRenhanced(
        const double                p_Luminosity, 
        const double                p_Temperature, 
        const double                p_HeAbundanceSurface, 
        const double                p_dMdtOther,
        const std::optional<double> p_dMdtWR
    ) const;
*/

    // mass loss RSG
    COMPAS_PURE MassLossT       CalculateMLrateRSG(
        const double                     p_Mass,
        const double                     p_Radius,
        const double                     p_Luminosity,
        const double                     p_Temperature,
        const double                     p_mStart,
        const RSG_MASS_LOSS_PRESCRIPTION p_RSG_MLprescription
    ) const;
    GNU_CONST   MassLossT       CalculateMLrateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const;
    GNU_CONST   MassLossT       CalculateMLrateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const;
    GNU_CONST   MassLossT       CalculateMLrateRSG_Kee2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    GNU_CONST   MassLossT       CalculateMLrateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const;
    GNU_CONST   MassLossT       CalculateMLrateRSG_Yang2023(const double p_Luminosity) const;

    // mass loss VMS
    COMPAS_PURE MassLossT       CalculateMLrateVMS(const double p_Mass, const double p_Luminosity, const double p_Temperature, const VMS_MASS_LOSS_PRESCRIPTION p_VMS_MLprescription) const;
    GNU_CONST   MassLossT       CalculateMLrateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const;
    COMPAS_PURE MassLossT       CalculateMLrateVMS_Sabhahit2023(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE MassLossT       CalculateMLrateVMS_Vink2011(const double p_Mass, const double p_Luminosity, const double p_Temperature) const;

    // mass loss WR
    GNU_CONST   MassLossT       CalculateMLrateWR_Hurley2000(const double p_Luminosity, const double p_PerturbationMu) const;
    COMPAS_PURE static MassLossT CalculateMLrateWR_SanderVink2020_Static(const double p_Luminosity, const double p_PerturbationMu) const;
    COMPAS_PURE static MassLossT CalculateMLrateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature) const;
    COMPAS_PURE static MassLossT CalculateMLrateWR_ZDependent_Static(const double p_Luminosity, const double p_PerturbationMu) const;

    
    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////
    COMPAS_PURE double              CalculateRadiusAtZAMS_Tout1996(const double p_MZAMS) const;


    ////////////////////////////////////////
    //   ROTATION                         //
    ////////////////////////////////////////

    COMPAS_PURE static double       CalculateAngularFrequency_CHE_Static(const double p_MZAMS) const;

    COMPAS_PURE double              CalculateRotationalVelocityOStar_Ramirez2013() const;


    ////////////////////////////////////////
    //   SUPERNOVAE                       //
    ////////////////////////////////////////

    COMPAS_PURE StellarSNDetailsT   CalculateSNkickMagnitude(
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












    
    
    














    
    



    double                      CalculateMassLossValues(double p_Dt, const bool p_UpdateMDot = false);

    
//// Never called    double                      CalculateNuclearMassLossRate()                                                  { return m_Mass / CalculateTimescale_RadialExpansion_DuringMT(); }
        

//// Never called   double                      CalculateRadialChange() const                                                   { return (utils::Compare(m_RadiusPrev,0) <= 0)? 0 : std::abs(m_Radius - m_RadiusPrev) / m_RadiusPrev; } // Return fractional radial change (if previous radius is negative or zero, return 0 to avoid NaN

    double                      CalculateTimescale_RadialExpansion_DuringMT();






 
    double                      CalculateTimestep();




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
    
    virtual void                UpdateGBparams(const double p_Mass, DblVectorT &p_GBparams)                  { }                                                                                  // Default is NO-OP
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

















    double              CalculateGBRadiusXExponent() const;

                                                 // Default is NO-OP




// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// removed in PR 1445    
//GNU_CONST double inline BaseStar::CalculateMLfractionWR(const double p_HeAbundanceSurface) const;




    void                CalculateLCoefficients(const double p_LogMetallicityXi, DblVectorT &p_LCoefficients) const;






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


    
    void                CalculateMassCutoffs(DblVectorT& p_MassCutoffs) const;





























   
    
    
    
    
    



    double              CalculateMassLossRateWolfRayetTemperatureCorrectionSander2023(const double p_Mdot) const;



// JR : DON'T FORGET "OVERRIDE" ON DERIVED FUNCTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



virtual double CalculateMomentOfInertia() const { return (0.1 * (m_Mass) * m_Radius * m_Radius); } // Defaults to MS. k2 = 0.1 as defined in Hurley et al. 2000, after eq 109
    virtual double              CalculateMomentOfInertiaAU() const                                              { return CalculateMomentOfInertia() * RSOL_TO_AU * RSOL_TO_AU; }


    GNU_CONST double CalculateEddingtonLuminosity(const double p_Mass, const double p_HeAbundanceSurface) const;













// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ********************************************************* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>























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

    void                CalculateRCoefficients(const double p_LogMetallicityXi, DblVectorT& p_RCoefficients) const;


GNU_PURE double CalculateAngularFrequencyAtZAMS(const double p_MZAMS, const double p_RZAMS) const;
GNU_PURE double CalculateAngularFrequencyAtZAMS_Hurley2000(const double p_MZAMS, const double p_RZAMS) const;
GNU_PURE double CalculateRotationalVelocityAtZAMS(double p_MZAMS) const;


/////  CALCULATE THIS?? - in UTILS  JR FIX THIS    GNU_CONST double    CalculateMinOmegaForCHE_Butler2018(const double p_RefZ, const double p_MZAMS) const

/////  CALCULATE THIS?? - in UTILS  JR FIX THIS    COMPAS_PURE inline double CalculateZAMSOmega(const double p_MZAMS, const double p_RZAMS)

/////  CALCULATE THIS?? - in UTILS  JR FIX THIS    double CalculateZAMSVrot(const double p_MZAMS) const;


///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



  






























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





    virtual double ChooseTimestep() const;


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
 * Recalculate the age of a star after mass loss.
 * 
 * Calls relevant function based on the evolutionary mode given in program options.
 *
 *
 * double CalculateAgeAfterMassLoss() const
 *
 * @return                                      Age of the star after mass loss (Myr)
 */
inline double BaseStar::CalculateAgeAfterMassLoss() const {

    double age;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            age = CalculateAgeAfterMassLoss_Hurley2000();
            break;

        default:                                                                        // unknown mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
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
 * CalculateTau
 *
 * @brief
 * Calculate the phase-relative age of the star (fractional age on the current
 * evolutionary phase).
 *
 * Calls relevant function based on the evolutionary mode given in program options.
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * DblVectorT CalculateTimescales() const
 *
 * @return                                      Timescales (Myr)
 */
inline DblVectorT BaseStar::CalculateTimescales() const { 

    DblVectorT tScales;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            tScales = CalculateTimescales_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
 * the stellar type has changed, or the radius has not changed, from the previous
 * timestep), the radial expansion timescale is returned as -1.0.
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
double BaseStar::CalculateTimescale_RadialExpansion_DuringMT() {
    
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
 * @return                                      Mass change timescale (Myr)
 * @return                                      eddy turnover timescale (yr)   // *Ilya* why is this in yr, not Myr?  Seems inconsistent.
 */
double BaseStar::CalculateTimescale_EddyTurnover_Hurley2002(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_CoreMass) const {

	const double rEnv = CalculateConvectiveEnvelopeRadialExtent();
    double mEnv;
    std::tie(mEnv, std::ignore) = CalculateConvectiveEnvelopeMass(p_Mass, p_CoreMass);
    return 0.4311 * std::cbrt((mEnv * rEnv * (p_Radius - (0.5 * rEnv))) / (3.0 * p_Luminosity));
}


/*
 * ChooseTimestep
 *
 * @brief
 * Choose timestep for evolution
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double ChooseTimestep() const
 *
 * @return                                      Suggested timestep (Myr)
 */
inline double BaseStar::ChooseTimestep() const { 

    DblVectorT GBparams;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            GBparams = ChooseTimestep_Hurley2000(m_InterimState.Age(), m_InterimState.Timescales());
            break;
        
        default:                                                                        // unknown mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
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
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * DblVectorT CalculateGBparams() const
 *
 * @return                                      Timescales (Myr)
 */
inline DblVectorT BaseStar::CalculateGBparams() const { 

    DblVectorT GBparams;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            GBparams = CalculateGBparams_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
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
 * Calls relevant function based on the evolutionary mode given in program options.
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
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

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            radius = CalculateRadiusAtPhaseEnd_Hurley2000();
            break;
        
        default:                                                                        // unknown mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
 * Calls relevant function based on the evolutionary mode given in program options.
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
        const double exponent = min((0.1 / CalculatePerturbationQ(p_Radius, p_Rc)), (-14.0 / std::log10(p_Mu)));

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
 * double CalculateLuminosity_Hurley2000(const double p_CoreMass, const DblVectorT& p_GBparams) const
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      AGB luminosity (Lsol)
 */
GNU_CONST inline double BaseStar::CalculateLuminosity_Hurley2000(const double p_CoreMass, const DblVectorT& p_GBparams) const {
#define GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)]
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
#define b(x) GLOBALS->HurleyBCoefficients(x)

    return p_Mass < p_MHeF
            ? (b(29) * PPOW(p_Mass, b(30))) / (1.0 + (GLOBALS->HurleyAlpha3() * std::exp(15.0 * (p_Mass - GLOBALS->HurleyMassCutoffs(static_cast<int>(MASS_CUTOFF::MHeF))))))
            : (b(31) + (b(32) * PPOW(p_Mass, (b(33) + 1.8)))) / (b(34) + PPOW(p_Mass, b(33)));

#undef bN
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    TEMPERATURE                                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateTemperature
 *
 * @brief
 * Calculate the temperature of the star at the current evolutionary phase.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateTemperature(const double p_Luminosity, const double p_Radius) const
 *
 * @return                                      Temperature of the star (Tsol)
 */
inline double BaseStar::CalculateTemperature(const double p_Luminosity, const double p_Radius) const { 

    double temperature;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            // retrieve luminosity and radius from current state and call relevant function
            temperature = CalculateTemperature_Hurley2000(p_Luminosity, p_Radius);
            break;
        
        default:                                                                        // unknown mode
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
 * Calls relevant function based on the evolutionary mode given in program options.
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
 * Calls relevant function based on the evolutionary mode given in program options.
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
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
 * double CalculateCoreMass_Hurley2000(const double p_Luminosity, const DblVectorT& p_GBparams)
 * 
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Core mass (Msol)
 */
inline double CalculateCoreMass_Hurley2000(const double p_Luminosity, const DblVectorT& p_GBparams) {
    return p_Luminosity > p_GBparams[HURLEY_GBP::LX]
            ? PPOW((p_Luminosity / p_GBparams[HURLEY_GBP::B]), (1.0 / p_GBparams[HURLEY_GBP::Q]))
            : PPOW((p_Luminosity / p_GBparams[HURLEY_GBP::D]), (1.0 / p_GBparams[HURLEY_GBP::P]));
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
 * Calls relevant function based on the evolutionary mode given in program options.
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
 * Calls relevant function based on the evolutionary mode given in program options.
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
 * Calls relevant function based on the evolutionary mode given in program options.
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
GNU_CONST inline double BaseStar::CalculateMLfractionWR(const double p_HeAbundanceSurface) const {

    constexpr double limOB = 0.55;                                          // per Yoon et al. 2006
    constexpr double limWR = 0.70;                                          // per Yoon et al. 2006

    return std::min(1.0, std::max (0.0, (p_HeAbundanceSurface - limOB) / (limWR - limOB)));
}
*/


/*
 * CalculateMLrate_Hurley2000
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, per Hurley et al. 2000.
 * 
 * 
 * MassLossT CalculateMLrate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_PerturbationMu) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu (optional parameter)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (may be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE inline MassLossT BaseStar::CalculateMLrate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_PerturbationMu) const {
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
 * MassLossT CalculateMLrate_KudritzkiReimers1978(const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::GB)
 * @return                                      Kudritzki and Reimers mass loss rate (in Msol yr^{-1})
 */
GNU_CONST inline MassLossT BaseStar::CalculateMLrate_KudritzkiReimers1978(const double p_Mass, const double p_Radius, const double p_Luminosity) const {
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
 * MassLossT CalculateMLrate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate on the AGB (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::GB)
 */
GNU_CONST inline MassLossT BaseStar::CalculateMLrate_VassiliadisWood1993(const double p_Mass, const double p_Radius, const double p_Luminosity) const {
    const double P0   = PPOW(10.0, std::min(3.3, (-2.07 - (0.9 * std::log10(p_Mass)) + (1.94 * std::log10(p_Radius)))));
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
 * MassLossT CalculateMLrateLBV_Belczynski2010() const
 *
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for LBV stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::LBV)
 */    
COMPAS_PURE inline MassLossT BaseStar::CalculateMLrateLBV_Belczynski2010() const {
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
 * MassLossT CalculateMLrateLBV_Hurley2000(const double p_Luminosity, const double p_HDlimitfactor) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_HDlimitfactor                 Factor by which star is above Humphreys-Davidson limit
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for LBV stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::LBV)
 */
GNU_CONST inline MassLossT BaseStar::CalculateMLrateLBV_Hurley2000(const double p_Luminosity, const double p_HDlimitfactor) const {
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
 * MassLossT CalculateMLrateOB_Krticka2018(const double p_Luminosity) const
 * 
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for hot OB stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::OB)
 */
COMPAS_PURE inline MassLossT BaseStar::CalculateMLrateOB_Krticka2018(const double p_Luminosity) const {
    return std::make_tuple(PPOW(10.0, -5.70 + 0.50 * GLOBALS->ZetaAsplund()() + (1.61 - 0.12 * GLOBALS->ZetaAsplund()()) * std::log10(p_Luminosity / 1.0E6)), MASS_LOSS_TYPE::OB);
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
 * MassLossT CalculateMLrateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const
 *
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MassLossT BaseStar::CalculateMLrateRSG_Beasor2020(const double p_mStart, const double p_Luminosity) const {
    return std::make_tuple(PPOW(10.0, -21.5 - 0.15 * p_mStart + 3.6 * std::log10(p_Luminosity)), MASS_LOSS_TYPE::RSG);
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
 * MassLossT CalculateMLrateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const
 *
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MassLossT BaseStar::CalculateMLrateRSG_Decin2023(const double p_mStart, const double p_Luminosity) const {
    return std::make_tuple(PPOW(10.0, -20.63 - 0.16 * p_mStart + 3.47 * std::log10(p_Luminosity)), MASS_LOSS_TYPE::RSG);
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
 * MassLossT CalculateMLrateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MassLossT BaseStar::CalculateMLrateRSG_VinkSabhahit2023(const double p_Mass, const double p_Luminosity) const {
    const double exp = logL < 4.6 ? -8.0 + 0.7 * std::log10(p_Luminosity) - 0.7 * std::log10(p_Mass) : -24.0 + 4.77 * std::log10(p_Luminosity) - 3.99 * std::log10(p_Mass);
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
 * GNU_CONST MassLossT CalculateMLrateRSG_Yang2023(const double p_Luminosity) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST inline MassLossT BaseStar::CalculateMLrateRSG_Yang2023(const double p_Luminosity) const {

    const double logL  = std::log10(p_Luminosity);
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
 * MassLossT CalculateMLrateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for very massive stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::VMS)
 */
GNU_CONST inline MassLossT BaseStar::CalculateMLrateVMS_Bestenlehner2020(const double p_Mass, const double p_Luminosity) const {

    constexpr double alpha       = 0.39;                                        // CAK force multiplier
    constexpr double logMdotZero = -4.78;                                       // from substituting LogMdotTrans and Gamma_e trans into eq 12. 

    const double gamma = EDDINGTON_PARAMETER_FACTOR * p_Luminosity / p_Mass;    // Eddington Parameter, not metallicity specific as in the publication

    return std::make_tuple(PPOW(10.0, logMdotZero + ((1.0 / alpha) + 0.5) * std::log10(gamma) - (((1.0 - alpha) / alpha) + 2.0) * std::log10(1.0 - gamma)), MASS_LOSS_TYPE::VMS);
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
 * MassLossT CalculateMLrateWR_Hurley2000(const double p_Luminosity, const double p_PerturbationMu) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
GNU_CONST inline MassLossT BaseStar::CalculateMLrateWR_Hurley2000(const double p_Luminosity, const double p_PerturbationMu) const {
    return std::make_tuple((p_PerturbationMu >= 1.0 ? 0.0 : PPOW(p_Luminosity, 1.5) * (1.0 - p_PerturbationMu) * 1.0E-13), MASS_LOSS_TYPE::WR);
} 


/*
 * CalculateMLrateWR_Shenar2019_Static
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
 * MassLossT CalculateMLrateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
COMPAS_PURE inline MassLossT BaseStar::CalculateMLrateWR_Shenar2019_Static(const double p_Luminosity, const double p_Temperature) const {

    // For H-rich WR stars (X_H > 0.4)
    constexpr double C1 = -6.78;
    constexpr double C2 =  0.66;
    constexpr double C3 = -0.12;
    constexpr double C5 =  0.74;

    return std::make_tuple(PPOW(10.0, C1 + (C2 * std::log10(p_Luminosity)) + (C3 * std::log10(p_Temperature * TSOL)) + (C5 * GLOBALS->SigmaHurley())), MASS_LOSS_TYPE::WR);
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
 * static MassLossT CalculateMLrateWR_ZDependent_Static(const double p_Luminosity, const double p_PerturbationMu)
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
COMPAS_PURE inline MassLossT BaseStar::CalculateMLrateWR_ZDependent_Static(const double p_Luminosity, const double p_PerturbationMu) {
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
 * Calls relevant function based on the evolutionary mode given in program options.
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
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
 * Calls relevant function based on the evolutionary mode given in program options.
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
            // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return luminosity;
}























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




    
    virtual Dbl_Dbl_Dbl_DblT     CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const;
    virtual Dbl_Dbl_Dbl_DblT     CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const ;
    virtual Dbl_Dbl_Dbl_DblT     CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const;
    

    virtual double              CalculateCELambda_Dewi() const                                                     { return 1.0; }           

                                              // Default for stellar types with no LamdaDewi definitions - 1.0 is benign

GNU_CONST double CalculateCELambdaKruckow(const double p_Radius, const double p_Alpha) const;
    double CalculateCELambdaKruckow() const { return CalculateCELambdaKruckow(m_Radius, OPTIONS->CommonEnvelopeSlopeKruckow()); }
    
    virtual double              CalculateLambdaLoveridge(const double p_EnvMass, const bool p_IsMassLoss = false) const { return 1.0; }                                             // Default for non giant branch stars - 1.0 is benign
    double                      CalculateLambdaLoveridge() const                                                { return CalculateLambdaLoveridge(m_Mass - m_CoreMass, false); }


//    double                      CalculateCELambda_Nanjing() const;
    


    Dbl_DblT                     CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate); 
    virtual Dbl_DblT             CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                            const double p_AccretorMassRate,
                                                            const bool   p_IsHeRich)                            { return CalculateMassAcceptanceRate(p_DonorMassRate, p_AccretorMassRate); } // Ignore the He content for non-WDs
    double                      CalculateMassAccretedForCO(const double p_Mass, const double p_CompanionMass, const double p_CompanionRadius, const double p_CompanionEnvelope) const;


COMPAS_PURE double BaseStar::CalculateCELambda_Nanjing(const double p_Mass, const double p_Radius, const double p_CoreMass) const;

COMPAS_PURE  virtual double CalculateCELambda_Nanjing_StarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const { return 1.0; } // Default for stellar types with no LamdaNanjing definitions - 1.0 is benign




GNU_PURE  virtual double CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
                                                            const double             p_Radius,
                                                            const double             p_CoreMass,
                                                            const SizeT              p_MassIndex,
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
 * CalculateConvectiveEnvelopeBindingEnergy_Hirai2022
 *
 * @brief
 * Calculate the convective envelope binding energy for the two-stage common envelope
 * formalism of Hirai & Mandel, 2022
 *
 *
 * double CalculateConvectiveEnvelopeBindingEnergy_Hirai2022(const double p_Mass, const double p_Radius, const double p_EnvMass, const double p_Lambda) const
 *
 * @param           p_Mass                      Mass (total) of the star (Msol)
 * @param           p_Radius                    Radius of the star (Rsol)
 * @param           p_EnvMass                   Mass of the star's convective outer envelope (Msol)
 * @param           p_Lambda                    Lambda parameter for the convective envelope
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
 * This is required for the Hirai & Mandel, 2022, two-stage CE formalism.
 *
 *
 * double BaseStar::CalculateConvectiveEnvelopeLambda_Picker2024(const double p_Mass, const double p_EnvMass, const double p_EnvMassMax) const
 *
 * @param           p_Mass                      Mass of the star (Msol)
 * @param           p_EnvMass                   Mass of the star's convective outer envelope (Msol)
 * @param           p_EnvMassMax                Maximum mass of the star's convective outer envelope (Msol)
 * @return                                      Lambda binding energy parameter for the convective outer envelope
 */
inline double BaseStar::CalculateConvectiveEnvelopeLambda_Picker2024(const double p_Mass, const double p_EnvMass, const double p_EnvMassMax) const {
    
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
