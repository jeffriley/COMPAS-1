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
    
public:
    
    BaseStar();
    
    BaseStar(const unsigned long int p_RandomSeed,
             const double            p_MZAMS,
             const double            p_Metallicity,
             const KickParameters    p_KickParameters,
             const double            p_RotationalVelocity = -1.0);
    
    
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
    
    virtual ~BaseStar() {}
    
    
    // object identifiers - all classes have these
    OBJECT_ID                   ObjectId() const                                        { return m_ObjectId; }
    OBJECT_TYPE                 ObjectType() const                                      { return OBJECT_TYPE::BASE_STAR; }
    OBJECT_PERSISTENCE          ObjectPersistence() const                               { return m_ObjectPersistence; }
    STELLAR_TYPE                InitialStellarType() const                              { return m_InitialStellarType; }
    STELLAR_TYPE                StellarType() const                                     { return m_StellarType; }
    STELLAR_TYPE                StellarTypePrev() const                                 { return m_StellarTypePrev; }
    
    
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
    void                        SetInitialType(const STELLAR_TYPE p_InitialType)                { m_InitialStellarType = p_InitialType; }
    void                        SetError(const ERROR p_Error)                                   { m_Error = p_Error; }
    void                        SetObjectId(const OBJECT_ID p_ObjectId)                         { m_ObjectId = p_ObjectId; }
    void                        SetPersistence(const OBJECT_PERSISTENCE p_Persistence)          { m_ObjectPersistence = p_Persistence; }
    
    void                        SetOmega(double p_Omega)                                        { SetAngularMomentum(CalculateMomentOfInertiaAU() * p_Omega); }
    
    void                        SetSNCurrentEvent(const SN_EVENT p_SNEvent)                     { m_SupernovaDetails.events.current |= p_SNEvent; }                 // Set supernova primary event/state for current timestep
    void                        SetSNPastEvent(const SN_EVENT p_SNEvent)                        { m_SupernovaDetails.events.past |= p_SNEvent; }                    // Set supernova primary event/state for any past timestep
    void                        ClearCurrentSNEvent()                                           { m_SupernovaDetails.events.current = SN_EVENT::NONE; }             // Clear supernova event/state for current timestep
    
    void                        SetEvolutionStatus(const EVOLUTION_STATUS p_EvolutionStatus)    { m_EvolutionStatus = p_EvolutionStatus; }                          // Set evolution status (typically final outcome) for star

    void                        SetDt(const double p_dt)                                        { m_dt = std::max(0.0, p_dt); }                                     // Set timestep - ensure >= 0.0
    void                        SetPrevDt(const double p_dt)                                    { m_dtPrev = p_dt; }                                                // Set previous timestep - don't clamp - preserve actuality

    void                        UpdateComponentVelocity(const Vector3d p_NewVelocity)           { m_ComponentVelocity += p_NewVelocity; }
    void                        UpdateMassTransferDonorHistory();
    
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
    //    - "Apply"     should update (modify) class member variables, and may return value(s)
    //    - "Update"    should update (modify) class member variables, and may return value(s)
    //
    // Functions that do not update (modify) class member variables should be declared const

    void                        ApplyMassTransferRejuvenationFactor()                                           { m_Age *= CalculateMassTransferRejuvenationFactor(); }             // Apply age rejuvenation factor
    
    double                      CalculateBindingEnergy(const double p_CoreMass, const double p_EnvMass, const double p_Radius, const double p_Lambda) const;

    virtual double              CalculateConvectiveCoreMass() const                                             {return 0.0;}
    virtual double              CalculateConvectiveCoreRadius() const                                           {return 0.0;}
    double                      CalculateConvectiveEnvelopeBindingEnergy(const double p_TotalMass, const double p_ConvectiveEnvelopeMass, const double p_Radius, const double p_Lambda) const;
    double                      CalculateConvectiveEnvelopeLambdaPicker(const DBL_DBL p_ConvectiveEnvelopeMass) const;
    virtual DBL_DBL             CalculateConvectiveEnvelopeMass() const                                         { return std::tuple<double, double> (0.0, 0.0); }
    
    virtual double              CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate, const double p_massTransferEfficiencyBeta) const;
    virtual double              CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const     { return 0.0; }                                                     // Default is 0.0
    virtual double              CalculateCriticalMassRatioGeEtAl(const QCRIT_PRESCRIPTION p_qCritPrescription, const double p_massTransferEfficiencyBeta) {
                                    return InterpolateGeEtAlQCrit(p_qCritPrescription, p_massTransferEfficiencyBeta);
                                }
    virtual double              CalculateCriticalMassRatioHurleyHjellmingWebbink() const                        { return 0.0; }                                                     // Default is 0.0
    
    double                      CalculateDynamicalTimescale() const                                             { return CalculateDynamicalTimescale_Static(m_Mass, m_Radius); }    // Use class member variables
    
    double                      CalculateEddingtonLuminosity() const                                            { return CalculateEddingtonLuminosity_Static(m_Mass, m_HeliumAbundanceSurface); } // Use class member variables
    
    double                      CalculateEddyTurnoverTimescale() const;
    
    double                      CalculateInitialHeliumAbundance() const                                         { return 0.24 + 2.0 * m_Metallicity; }                              // Pols et al. 1998
    double                      CalculateInitialHydrogenAbundance() const                                       { return 0.76 - 3.0 * m_Metallicity; }                              // Pols et al. 1998
    
    virtual DBL_DBL_DBL_DBL     CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const;
    virtual DBL_DBL_DBL_DBL     CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const ;
    virtual DBL_DBL_DBL_DBL     CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const;
    
    virtual double              CalculateLambdaDewi() const                                                     { return 1.0; }                                                     // Default for stellar types with no LamdaDewi definitions - 1.0 is benign
    double                      CalculateLambdaKruckow(const double p_Radius, const double p_Alpha) const;
    double                      CalculateLambdaKruckow() const                                                  { return CalculateLambdaKruckow(m_Radius, OPTIONS->CommonEnvelopeSlopeKruckow()); }
    virtual double              CalculateLambdaLoveridge(const double p_EnvMass, const bool p_IsMassLoss = false) const { return 1.0; }                                             // Default for non giant branch stars - 1.0 is benign
    double                      CalculateLambdaLoveridge() const                                                { return CalculateLambdaLoveridge(m_Mass - m_CoreMass, false); }
    double                      CalculateLambdaNanjing() const;
    
    DBL_DBL                     CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate); 
    virtual DBL_DBL             CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                            const double p_AccretorMassRate,
                                                            const bool   p_IsHeRich)                            { return CalculateMassAcceptanceRate(p_DonorMassRate, p_AccretorMassRate); } // Ignore the He content for non-WDs
    double                      CalculateMassAccretedForCO(const double p_Mass, const double p_CompanionMass, const double p_CompanionRadius, const double p_CompanionEnvelope) const;
    double                      CalculateMassChangeTimescale() const                                            { return CalculateMassChangeTimescale_Static(m_StellarType, m_StellarTypePrev, m_Mass, m_MassPrev, m_dtPrev); }  // Use class member variables
    double                      CalculateMassLossValues(double p_Dt, const bool p_UpdateMDot = false);

    virtual double              CalculateMomentOfInertia() const                                                { return (0.1 * (m_Mass) * m_Radius * m_Radius); }                  // Defaults to MS. k2 = 0.1 as defined in Hurley et al. 2000, after eq 109
    virtual double              CalculateMomentOfInertiaAU() const                                              { return CalculateMomentOfInertia() * RSOL_TO_AU * RSOL_TO_AU; }
    
    double                      CalculateNuclearMassLossRate()                                                  { return m_Mass / CalculateRadialExpansionTimescaleDuringMassTransfer(); }
        
    double                      CalculateOmegaCHE(const double p_MZAMS, const double p_Metallicity) const;

    double                      CalculateRadialChange() const                                                   { return (utils::Compare(m_RadiusPrev,0) <= 0)? 0 : std::abs(m_Radius - m_RadiusPrev) / m_RadiusPrev; } // Return fractional radial change (if previous radius is negative or zero, return 0 to avoid NaN
    double                      CalculateRadialExpansionTimescale() const                                       { return CalculateRadialExpansionTimescale_Static(m_StellarType, m_StellarTypePrev, m_Radius, m_RadiusPrev, m_dtPrev); } // Use class member variables
    double                      CalculateRadialExpansionTimescaleDuringMassTransfer();
    virtual double              CalculateRadialExtentConvectiveEnvelope() const                                 { return 0.0; }                                                     // Default for stars with no convective envelope
    virtual double              CalculateRadiusOnMassChange(double p_dM)                                        { return Radius(); }                                                // NO-OP
    virtual double              CalculateRemnantRadius() const                                                  { return Radius(); }                                                // Relevant for MS stars, over-written for GB stars

    double                      CalculateSNKickMagnitude(const double p_RemnantMass, const double p_EjectaMass, const STELLAR_TYPE p_StellarType);

    double                      CalculateThermalMassAcceptanceRate(const double p_Radius);
    double                      CalculateThermalMassAcceptanceRate()                                            { return CalculateThermalMassAcceptanceRate(m_Radius); }
    virtual double              CalculateThermalMassLossRate() const                                            { return m_Mass / CalculateThermalTimescale(); }                    // Use class member variables - and inheritance hierarchy
    virtual double              CalculateThermalTimescale(const double p_Radius) const;                                                                                             // Use inheritance hierarchy
    virtual double              CalculateThermalTimescale() const                                               { return CalculateThermalTimescale(m_Radius); }                     // Use inheritance hierarchy

    virtual double              CalculateTAMSCoreMass() const                                                   { return 0.0; }                                                     // Except MS stars

    double                      CalculateTimestep();
    
    double                      CalculateZetaAdiabatic() const;
    virtual double              CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription)          { return 0.0; }                                                     // Use inheritance hierarchy
    virtual double              CalculateZetaEquilibrium()                                                      { return 0.0; }



    virtual ACCRETION_REGIME    DetermineAccretionRegime(const double p_DonorThermalMassLossRate, 
                                                         const bool p_HeRich)                                   { return ACCRETION_REGIME::ZERO; }                                  // Placeholder, use inheritance for WDs

    virtual ENVELOPE            DetermineEnvelopeType() const                                                   { return ENVELOPE::REMNANT; }                                       // Default is REMNANT - but should never be called
    virtual MT_CASE             DetermineMassTransferTypeAsDonor() const                                        { return MT_CASE::OTHER; }                                          // Not A, B, C, or NONE

    STELLAR_TYPE                EvolveOneTimestep(const double p_dM, const double p_dM0, const double p_dt, const bool p_ForceRecalculate = false);
    
    void                        HaltWinds()                                                                     { m_Mdot = 0.0; }                                                   // Disable wind mass loss in current time step (e.g., if star is a donor or accretor in a RLOF episode)

    virtual double              InterpolateGeEtAlQCrit(const QCRIT_PRESCRIPTION p_qCritPrescription, 
                                                       const double p_massTransferEfficiencyBeta)               { return 0.0; }                                                     // Placeholder, use interpolator for either H-rich or H-poor stars


    void                        ResetEnvelopeExpulsationByPulsations()                                          { m_EnvelopeJustExpelledByPulsations = false; }

    void                        ResolveAccretion(const double p_AccretionMass)                                  { m_Mass = std::max(0.0, m_Mass + p_AccretionMass); }               // Handles donation and accretion - won't let mass go negative
    virtual void                ResolveAccretionRegime(const ACCRETION_REGIME p_Regime, const double p_DonorThermalMassLossRate) { }                                                // Default does nothing, only works for WDs.
    virtual double              ResolveCommonEnvelopeAccretion(const double p_FinalMass,
                                                               const double p_CompanionMass     = 0.0,
                                                               const double p_CompanionRadius   = 0.0,
                                                               const double p_CompanionEnvelope = 0.0)          { return p_FinalMass - Mass(); }                                    // Overwritten in NS.h; for now, no accretion on stars other than compact objects during CE
    virtual STELLAR_TYPE        ResolveEnvelopeLoss(bool p_Force = false)                                       { return m_StellarType; }
    virtual STELLAR_TYPE        ResolveMassLoss(const double p_dt);
    virtual STELLAR_TYPE        ResolveMassLossHurley(const double p_dt);
    virtual void                ResolveShellChange(const double p_AccretedMass) { }                                                                                                 // Default does nothing, use inheritance for WDs.
    virtual STELLAR_TYPE        ResolveSupernova()                                                              { return m_StellarType; }                                           // Default is NO-OP
       
    void                        SetStellarTypePrev(const STELLAR_TYPE p_StellarTypePrev)                        { m_StellarTypePrev = p_StellarTypePrev; }
    
    virtual bool                ShouldEnvelopeBeExpelledByPulsations() const                                    { return false; }                                                   // Default is that there is no envelope expulsion by pulsations

    virtual void                SpinDownIsolatedPulsar(const double p_Stepsize) { }                                                                                                 // Default is NO-OP
    
    virtual void                UpdateGBParams(const double p_Mass, DBL_VECTOR &p_GBParams)                  { }                                                                                  // Default is NO-OP
    virtual void                UpdateGBParams()                                                             { UpdateGBParams(m_Mass0, m_GBParams); }                         // Use class member variables
    
    virtual void                UpdateAfterMerger(double p_Mass, double p_HydrogenMass) { }                                                                                         // Default is NO-OP
    virtual void                UpdateAgeAfterMassLoss() { }                                                                                                                        // Default is NO-OP
    virtual double              CalculateEffectiveInitialMass()                                                 { return m_Mass0; }                                                 // Default is NO-OP
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
    STELLAR_TYPE            m_InitialStellarType;                       // Stellar type at birth, defined in Hurley et al. 2000
    STELLAR_TYPE            m_StellarType;                              // Stellar type defined in Hurley et al. 2000
    ERROR                   m_Error;                                    // Records most recent error encountered for this star
  
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
    DBL_VECTOR              m_GBParams;                                 // Giant Branch Parameters
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

    // Stellar details squirrelled away...
    SupernovaDetailsT       m_SupernovaDetails;                         // Supernova attributes
    PulsarDetailsT          m_PulsarDetails;                            // Pulsar attributes

    // Star vector velocity 
    Vector3d                m_ComponentVelocity;                        // Isolated star velocity vector (binary's centre-of-mass velocity for bound binary)

    // Star mass transfer history 
    ST_VECTOR               m_MassTransferDonorHistory;                 // List of MT donor stellar types - mostly relevant for binary stars


    // Protected member function prototypes - alphabetically (implementations are in the cpp file)
    // See below (after class declaration) for inline functions declared in this header file

    double              CalculateAlpha1() const;                        // inline implementation below
    double              CalculateAlpha3() const;                        // inline implementation below
    double              CalculateAlpha4() const;                        // inline implementation below

    void                CalculateAnCoefficients(DBL_VECTOR &p_An, DBL_VECTOR &p_LConstants, DBL_VECTOR &p_RConstants, DBL_VECTOR &p_GammaConstants) const;

    double       CalculateBindingEnergy(const double p_Lambda) const                                             { return CalculateBindingEnergy(m_CoreMass, m_Mass - m_CoreMass, m_Radius, p_Lambda); }

    void                CalculateBnCoefficients(DBL_VECTOR &p_Bn) const;

    virtual double      CalculateCOCoreMassAtPhaseEnd() const                                                           { return m_COCoreMass; }                                                    // Default is NO-OP
    virtual double      CalculateCOCoreMassOnPhase() const                                                              { return m_COCoreMass; }                                                    // Default is NO-OP

    double       CalculateConvectiveEnvelopeBindingEnergy(const double p_Lambda) const                           { 
                                                                                                                            double convectiveEnvMass;
                                                                                                                            std::tie(convectiveEnvMass, std::ignore) = CalculateConvectiveEnvelopeMass();
                                                                                                                            return CalculateConvectiveEnvelopeBindingEnergy(m_Mass, convectiveEnvMass, m_Radius, p_Lambda);
                                                                                                                        }

    virtual double      CalculateCoreMassAtPhaseEnd() const                                                             { return m_CoreMass; }                                                      // Default is NO-OP
    static  double      CalculateCoreMassGivenLuminosity_Static(const double p_Luminosity, const DBL_VECTOR &p_GBParams){                                                                           // Hurley et al. 2000, eqs 37 & 38
                                                                                                                            return (utils::Compare(p_Luminosity, p_GBParams[static_cast<int>(GBP::Lx)]) > 0)
                                                                                                                                ? PPOW((p_Luminosity / p_GBParams[static_cast<int>(GBP::B)]), (1.0 / p_GBParams[static_cast<int>(GBP::q)]))
                                                                                                                                : PPOW((p_Luminosity / p_GBParams[static_cast<int>(GBP::D)]), (1.0 / p_GBParams[static_cast<int>(GBP::p)]));
                                                                                                                        }
    virtual double      CalculateCoreMassOnPhase() const                                                                { return m_CoreMass; }                                                      // Default is NO-OP

    static  double      CalculateDynamicalTimescale_Static(const double p_Mass, const double p_Radius);

    virtual double      CalculateEddingtonCriticalRate() const                                                          { return 2.08E-3 / 1.7 * m_Radius * MYR_TO_YEAR * OPTIONS->EddingtonAccretionFactor() ; } // Hurley+, 2002, Eq. (67)
    static  double      CalculateEddingtonLuminosity_Static(const double p_Mass, const double p_HeliumAbundanceSurface);

    double              CalculateGBRadiusXExponent() const;

    virtual double      CalculateHeCoreMassAtPhaseEnd() const                                                           { return m_HeCoreMass; }                                                    // Default is NO-OP
    virtual double      CalculateHeCoreMassOnPhase() const                                                              { return m_HeCoreMass; }                                                    // Default is NO-OP

    virtual double      CalculateHeliumAbundanceCoreOnPhase() const                                                     { return m_HeliumAbundanceCore; }                                           // Default is NO-OP
    virtual double      CalculateHeliumAbundanceSurfaceOnPhase() const                                                  { return m_HeliumAbundanceSurface; }                                        // Default is NO-OP

    virtual double      CalculateHydrogenAbundanceCoreOnPhase() const                                                   { return m_HydrogenAbundanceCore; }                                         // Default is NO-OP
    virtual double      CalculateHydrogenAbundanceSurfaceOnPhase() const                                                { return m_HydrogenAbundanceSurface; }                                      // Default is NO-OP    

    static  double      CalculateHeRateConstant_Static()                                                                { return HE_RATE_CONSTANT; }                                                // Only >= CHeB stars need AHe, but no drama if other stars calculate (retrieve it) - it's only a constant (we could just use the constant inline...)
    static  double      CalculateHHeRateConstant_Static()                                                               { return HHE_RATE_CONSTANT; }                                               // Only TPAGB stars need AHHe, but no drama if other stars calculate (retrieve it) - it's only a constant (we could just use the constant inline...)

    virtual double      CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const          { return 1.0; }                                                             // Default for stellar types with no LamdaNanjing definitions - 1.0 is benign
    virtual double      CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const { return 1.0; }                                                          // Default for stellar types with no LamdaNanjing definitions - 1.0 is benign

    void                CalculateLCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_LCoefficients) const;

    double              CalculateLifetimeToBAGB(const double p_tHeI, const double p_tHe) const;
    double              CalculateLifetimeToBGB(const double p_Mass) const;

    double              CalculateLuminosityAtBAGB(double p_Mass) const;         // inline implementation below

    virtual double      CalculateLuminosityAtPhaseEnd() const                                                           { return m_Luminosity; }                                                    // Default is NO-OP
    double              CalculateLuminosityAtZAMS(const double p_MZAMS) const;
    double       CalculateLuminosityGivenCoreMass(const double p_CoreMass) const                                 { return std::min((m_GBParams[static_cast<int>(GBP::B)] * PPOW(p_CoreMass, m_GBParams[static_cast<int>(GBP::q)])), (m_GBParams[static_cast<int>(GBP::D)] * PPOW(p_CoreMass, m_GBParams[static_cast<int>(GBP::p)]))); } // Hurley et al. 2000, eq 37

    virtual double      CalculateLuminosityOnPhase() const                                                              { return m_Luminosity; }                                                    // Default is NO-OP

    double              CalculateMassAndZInterpolatedLambdaNanjing(const double p_Mass, const double p_Z) const;
    double              CalculateMassInterpolatedLambdaNanjing(const double p_Mass, const STELLAR_POPULATION p_StellarPop) const;
    double              CalculateZInterpolatedLambdaNanjing(const double p_Z, const int p_MassIndex) const;

    static  double      CalculateMassChangeTimescale_Static(const STELLAR_TYPE p_StellarType,
                                                            const STELLAR_TYPE p_StellarTypePrev,
                                                            const double       p_Mass,
                                                            const double       p_MassPrev,
                                                            const double       p_DtPrev);
    
    void                CalculateMassCutoffs(const double p_Metallicity, const double p_LogMetallicityXi, DBL_VECTOR &p_MassCutoffs) const;


    virtual std::tuple<double, MASS_LOSS_TYPE>  CalculateMassLossRate() const;
    virtual double                              CalculateMassLossRateBelczynski2010() const;

    /*
     * CalculateMassLossRateEnhancementForRotatingStars
     *
     * Calculate mass loss rate enhancement for rapidly rotating stars
     *
     * Langer 1998 (https://ui.adsabs.harvard.edu/abs/1998A%26A...329..551L/abstract) eq 3
     * 
     * The exponent originally comes from Bjorkman & Cassinelli 1993 (https://ui.adsabs.harvard.edu/abs/1993ApJ...409..429B/abstract),
     * based on a fit to data from Friend & Abbott 1986 (https://ui.adsabs.harvard.edu/abs/1986ApJ...311..701F/abstract) 
     *
     * Uses no class member variables (directly).
     * 
     * 
     * double CalculateMassLossRateEnhancementForRotatingStars() const
     *
     * @return                                      Mass loss enhancement factor for rapidly rotating stars
     */
    double CalculateMassLossRateEnhancementForRotatingStars() const {
        return OPTIONS->EnableRotationallyEnhancedMassLoss() ? PPOW((1.0 - Omega() / OmegaBreak()), -0.43) : 1.0;
    }

    /*
     * CalculateMassLossRateHeliumStarVink2017
     *
     * Calculate the mass-loss rate for helium stars per Vink 2017
     * https://ui.adsabs.harvard.edu/abs/2017A%26A...607L...8V/abstract eq. 1
     *
     * Uses current values of:
     * 
     *    - m_Luminosity
     * 
     * 
     * double CalculateMassLossRateHeliumStarVink2017() const
     *
     * @return                                      Mass loss rate (Msol yr^-1)
     */
    double CalculateMassLossRateHeliumStarVink2017() const {
        return PPOW(10.0, -13.3 + 1.36 * log10(m_Luminosity) + 0.61 * LogMetallicityXiAnders());
    }


    /*
     * CalculateMassLossRateHurley
     *
     * Calculate the dominant mass loss mechanism and associated rate for the star
     * at the current evolutionary phase, per Hurley et al. 2000
     *
     * 
     * double CalculateMassLossRateHurley()
     *
     * @return                                      Mass loss rate in Msol per year
     */
    virtual std::tuple<double, MASS_LOSS_TYPE> CalculateMassLossRateHurley() {
        return CalculateMassLossRateNieuwenhuijzenDeJager();
    }




    double       CalculateMassLossRateKudritzkiReimers() const                                                   { return 4.0E-13 * (MASS_LOSS_ETA * m_Luminosity * m_Radius / m_Mass); }    // Hurley et al. 2000, eq 106 (based on a prescription taken from Kudritzki and Reimers 1978). Note: shouldn't be eta squared like in paper!}




    std::tuple<double, MASS_LOSS_TYPE> CalculateMassLossRateLBV(const LBV_MASS_LOSS_PRESCRIPTION p_LBVprescription);


    /*
     * CalculateMassLossRateLBVBelczynski
     *
     * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994)
     * per Belczynski et al. 2010, eq 8
     *
     * Does not use any class member variables 
     * 
     *
     * double CalculateMassLossRateLBVBelczynski() const
     *
     * @return                                      LBV-like mass loss rate (Msol yr^-1)
     */    
    double CalculateMassLossRateLBVBelczynski() const {
        return OPTIONS->LuminousBlueVariableFactor() * 1.0E-4;
    } 

    /*
     * CalculateMassLossRateLBVHurley
     *
     * Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994)
     * per Hurley+ 2000, section 7.1, a few equations after Eq. 106 (Equation not labelled)
     *
     * Uses current values of:
     * 
     *    - m_Luminosity
     *  
     *
     * double CalculateMassLossRateLBVHurley(const double p_HDlimitfactor) const
     *
     * @param   [IN]    p_HDlimitfactor             Factor by which star is above Humphreys-Davidson limit
     * @return                                      LBV-like mass loss rate (Msol yr^-1)
     */
    double CalculateMassLossRateLBVHurley(const double p_HDlimitfactor) const {
        double v = p_HDlimitfactor - 1.0;
        return 0.1 * v * v * v * ((m_Luminosity / 6.0E5) - 1.0);
    }


    virtual double      CalculateMassLossRateMerritt2025();
    double              CalculateMassLossRateNieuwenhuijzenDeJager() const;
    std::tuple<double, MASS_LOSS_TYPE> CalculateMassLossRateOB(const OB_MASS_LOSS_PRESCRIPTION p_OBMassLossPrescription) const;
    std::tuple<double, MASS_LOSS_TYPE> CalculateMassLossRateOBBjorklund2022() const;
    std::tuple<double, MASS_LOSS_TYPE> CalculateMassLossRateOBKrticka2018() const;



    /*
     * CalculateMassLossRateOBKrticka2018
     *
     * Calculate mass loss rate for massive OB stars per Krticka+ 2018
     * https://arxiv.org/pdf/1712.03321.pdf
     *
     * Uses current values of:
     * 
     *    - m_Luminosity
     *
     * 
     * std::tuple<double, MASS_LOSS_TYPE> CalculateMassLossRateOBKrticka2018() const
     * 
     * @return                                      Tuple containing:
     *                                                   DOUBLE         Mass loss rate for hot OB stars (Msol yr^-1)
     *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::OB)
     */
    std::tuple<double, MASS_LOSS_TYPE> CalculateMassLossRateOBKrticka2018() const {
        return std::make_tuple(PPOW(10.0, -5.70 + 0.50 * LogMetallicityXiAsplund() + (1.61 - 0.12 * LogMetallicityXiAsplund()) * log10(m_Luminosity / 1.0E6)), MASS_LOSS_TYPE::OB);
    }

    std::tuple<double, MASS_LOSS_TYPE>              CalculateMassLossRateOBVink2001() const;
    std::tuple<double, MASS_LOSS_TYPE>              CalculateMassLossRateOBVinkSander2021() const;
    std::tuple<double, MASS_LOSS_TYPE>              CalculateMassLossRateRSG(const RSG_MASS_LOSS_PRESCRIPTION p_RSG_mass_loss) const;

    /*
     * CalculateMassLossRateRSGBeasor2020
     *
     * Calculate mass loss rate for RSG stars per Beasor+2020
     * https://arxiv.org/pdf/2001.07222.pdf eq 4.
     * 
     * fit corrected slightly in Decin 2023, eq E.1 
     * https://arxiv.org/pdf/2303.09385.pdf
     * 
     * corrected again by Beasor+2023
     * https://ui.adsabs.harvard.edu/abs/2023MNRAS.524.2460B/abstract
     *
     * Uses current values of:
     * 
     *    - m_Luminosity
     *    - m_MZAMS
     * 
     *  
     * double CalculateMassLossRateRSGBeasor2020() const
     *
     * @return                                      Tuple containing:
     *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
     *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
     */
    std::tuple<double, MASS_LOSS_TYPE> CalculateMassLossRateRSGBeasor2020() const {
        return std::make_tuple(PPOW(10.0, -21.5 - 0.15 * m_MZAMS + 3.6 * log10(m_Luminosity)), MASS_LOSS_TYPE::RSG);
    }

    /*
     * CalculateMassLossRateRSGDecin2023
     *
     * Calculate mass loss rate for RSG stars per Decin 2023
     * https://arxiv.org/pdf/2303.09385.pdf eq 6.
     *
     * Uses current values of:
     * 
     *    - m_Luminosity
     *    - m_MZAMS
     *
     * 
     * double CalculateMassLossRateRSGDecin2023()c onst
     *
     * @return                                      Mass loss rate for RSG stars (Msol yr^-1)
     */

    double CalculateMassLossRateRSGDecin2023() const {
        return PPOW(10.0, -20.63 - 0.16 * m_MZAMS + 3.47 * log10(m_Luminosity));
    }

    double              CalculateMassLossRateRSGKee2021() const;
    double              CalculateMassLossRateRSGVinkSabhahit2023() const;
    double              CalculateMassLossRateRSGYang2023() const;
    double              CalculateMassLossRateVassiliadisWood() const;
    double              CalculateMassLossRateVMS(const VMS_MASS_LOSS_PRESCRIPTION p_VMS_mass_loss);
    double              CalculateMassLossRateVMSBestenlehner2020() const;
    double              CalculateMassLossRateVMSSabhahit2023() const;
    double              CalculateMassLossRateVMSVink2011() const;




    /*
     * CalculateMassLossRateWolfRayet
     *
     * Calculate the Wolf-Rayet like mass loss rate for small hydrogen-envelope mass (when mu < 1.0).
     *
     * Hurley et al. 2000, just after eq 106 (taken from Hamann, Koesterke & Wessolowski 1995, Hamann & Koesterke 1998)
     * 
     * In Hurley's fortran code there is a parameter 'hewind' which by default is 1.0, but it can be set to zero to
      * disable this particular part of winds. We instead opt for all winds on or off.
     *
     * Note that the reduction of this formula is imposed to match the observed number of black holes in binaries (Hurley et al. 2000)
     *
     * Uses current values of:
     * 
     *    - m_Luminosity
     * 
     *
     * double CalculateMassLossRateWolfRayet(const double p_Mu) const
     *
     * @param   [IN]    p_Mu                        Small envelope parameter (see Hurley et al. 2000, eq 97 & 98)
     * @return                                      Mass loss rate (Msol yr^-1)
     */
    double BaseStar::CalculateMassLossRateWolfRayet(const double p_Mu) const {
        return p_Mu >= 1.0 ? 0.0 : PPOW(m_Luminosity, 1.5) * (1.0 - p_Mu) * 1.0E-13;    // don't use utils::Compare() for thresholds
    }    
    
    
    
    
    
    double              CalculateMassLossRateWolfRayetSanderVink2020(const double p_Mu) const;
    virtual double      CalculateMassLossRateWolfRayetShenar2019() const;
    double              CalculateMassLossRateWolfRayetTemperatureCorrectionSander2023(const double p_Mdot) const;


    /*
     * CalculateMassLossRateWolfRayetZDependent
     *
     * Calculate the Wolf-Rayet like mass loss rate for small hydrogen-envelope mass (when mu < 1.0).
     * Belczynski et al. 2010, eq 9 (taken from Hamann, Koesterke & Wessolowski 1995, Hamann & Koesterke 1998)
     *
     * Note that the reduction of this formula is imposed to match the observed number of black holes in binaries (Hurley et al. 2000)
     *
     * Uses current values of:
     * 
     *    - m_Luminosity
     *    - m_Metallicity
     *
     *
     * double CalculateMassLossRateWolfRayetZDependent(const double p_Mu) const
     *
     * @param   [IN]    p_Mu                        Small envelope parameter (see Hurley et al. 2000, eq 97 & 98)
     * @return                                      Mass loss rate (Msol yr^-1)
     */
    double BaseStar::CalculateMassLossRateWolfRayetZDependent(const double p_Mu) const {
    // StarTrack may still do something different here.  There are references to Hamann & Koesterke 1998
    // and Vink and de Koter 2005.
    // TW - Haven't seen StarTrack but I think H&K gives the original equation and V&dK gives the Z dependence
        return p_Mu >= 1.0 ? 0.0 : 1.0E-13 * PPOW(m_Luminosity, 1.5) * PPOW(m_Metallicity / ZSOL_ANDERS, 0.86) * (1.0 - p_Mu);  // don't use utils::Compare() for thresholds
    }    



    virtual double      CalculateMassTransferRejuvenationFactor()                                                       { return 1.0; }

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

    virtual double      CalculateRadiusAtPhaseEnd() const                                                               { return m_Radius; }                                                        // Default is NO-OP
    double              CalculateRadiusAtZAMS(const double p_MZAMS) const;
    virtual double      CalculateRadiusOnPhase() const                                                                  { return m_Radius; }                                                        // Default is NO-OP
    virtual std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                              { return std::make_tuple(CalculateRadiusOnPhase(), m_StellarType); }

    void                CalculateRCoefficients(const double p_LogMetallicityXi, DBL_VECTOR &p_RCoefficients) const;

    double              CalculateRotationalVelocity(double p_MZAMS);

    virtual double      CalculateTauOnPhase() const                                                                     { return m_Tau; }                                                           // Default is NO-OP
    virtual double      CalculateTauAtPhaseEnd() const                                                                  { return m_Tau; }                                                           // Default is NO-OP

    virtual double      CalculateTemperatureAtPhaseEnd() const                                                          { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }
    virtual double      CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const          { return CalculateTemperatureOnPhase(p_Luminosity, p_Radius); }             // Same as on phase
    double              CalculateTemperatureKelvinOnPhase(const double p_Luminosity, const double p_Radius) const;
    virtual double      CalculateTemperatureOnPhase() const                                                             { return CalculateTemperatureOnPhase(m_Luminosity, m_Radius); }
    virtual double      CalculateTemperatureOnPhase(const double p_Luminosity, const double p_Radius) const             { return CalculateTemperatureOnPhase_Static(p_Luminosity, p_Radius); }
    static  double      CalculateTemperatureOnPhase_Static(const double p_Luminosity, const double p_Radius);

    virtual void        CalculateTimescales()                                                                           { CalculateTimescales(m_Mass0, m_Timescales); }                             // Use class member variables
    virtual void        CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) { }                                                                                                      // Default is NO-OP

    double              CalculateZAMSAngularFrequency(const double p_MZAMS, const double p_RZAMS);

    double              CalculateZetaAdiabaticHurley2002(const double p_CoreMass) const;
    double              CalculateZetaAdiabaticSPH(const double p_CoreMass) const;

    virtual double      ChooseTimestep(const double p_Time) const                                                       { return m_dt; }

    double              DrawKickMagnitudeBrayEldridge(const double p_EjectaMass, const double p_RemnantMass, const double p_Alpha, const double p_Beta) const;
    double              DrawKickMagnitudeDistributionFlat(const double p_MaxVK, const double p_Rand) const;
    double              DrawKickMagnitudeDistributionMaxwell(const double p_Sigma, const double p_Rand) const;
    double              DrawRemnantKickMuller(const double p_COCoreMass) const;
    double              DrawRemnantKickMullerMandel(const double p_COCoreMass, const double p_Rand, const double p_RemnantMass) const;
    double              DrawSNKickMagnitude(const double p_Sigma, const double p_COCoreMass, const double p_Rand, const double p_EjectaMass, const double p_RemnantMass) const;

    STELLAR_TYPE        EvolveOnPhase(const double p_DeltaTime);

    virtual STELLAR_TYPE EvolveToNextPhase()                                                                            { return m_StellarType; }

    double              FindLambdaNanjingNearestMassIndex(const double p_Mass) const;

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



/*
 * CalculateAlpha1
 *
 * Calculate the constant alpha1
 *
 * Hurley et al, 2000, just after eq 49
 *
 * Alpha1 depends on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 *
 * double CalculateAlpha1() const
 *
 * @return                                      Metallicity dependent constant alpha1
 */
inline double BaseStar::CalculateAlpha1() const {
    double LHeI_MHeF = (m_Bn[11] + (m_Bn[12] * PPOW(m_MassCutoffs(MHeF), 3.8))) / (m_Bn[13] + (m_MassCutoffs(MHeF) * m_MassCutoffs(MHeF)));
    return ((m_Bn[9] * PPOW(m_MassCutoffs(MHeF), m_Bn[10])) - LHeI_MHeF) / LHeI_MHeF;
}


/*
 * CalculateAlpha3
 *
 * Calculate the constant alpha3
 *
 * Hurley et al. 2000, just after eq 56
 *
 * Alpha3 depends on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 *
 * double CalculateAlpha3() const
 *
 * @return                                      Metallicity dependent constant alpha3
 */
inline double BaseStar::CalculateAlpha3() const {
    double LBAGB = (m_Bn[31] + (m_Bn[32] * PPOW(m_MassCutoffs(MHeF), (m_Bn[33] + 1.8)))) / (m_Bn[34] + PPOW(m_MassCutoffs(MHeF), m_Bn[33]));
    return ((m_Bn[29] * PPOW(m_MassCutoffs(MHeF), m_Bn[30])) - LBAGB) / LBAGB;
}


/*
 * CalculateAlpha4
 *
 * Calculate the constant alpha4
 *
 * Hurley et al. 2000, just after eq 57
 *
 * Alpha4 depends on a star's metallicity only - so this only needs to be done once per star (upon creation)
 *
 *
 * double CalculateAlpha4() const
 *
 * @return                                      Metallicity dependent constant alpha4
 */
inline double BaseStar::CalculateAlpha4() const {

    double MHeF      = m_MassCutoffs(MHeF);
    double MHeF_5    = MHeF * MHeF * MHeF * MHeF * MHeF;    // pow() is slow - use multiplication
    double tBGB_MHeF = CalculateLifetimeToBGB(MHeF);        // tBGB for mass M = MHeF
    double tHe_MHeF  = tBGB_MHeF * (m_Bn[41] * PPOW(MHeF, m_Bn[42]) + m_Bn[43] * MHeF_5) / (m_Bn[44] + MHeF_5);
    
    return ((tHe_MHeF - m_Bn[39]) / m_Bn[39]);
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
 * CalculateLuminosityAtBAGB
 *
 * Calculate luminosity at the base of the Asymptotic Giant Branch
 *
 * Hurley et al. 2000, eq 56
 *
 *
 * double CalculateLuminosityAtBAGB(double p_Mass) const
 *
 * @param   [IN]    p_Mass                      (Effective) mass in Msol
 * @return                                      Luminosity at BAGB in Lsol
 */
inline double BaseStar::CalculateLuminosityAtBAGB(double p_Mass) const {
    return (utils::Compare(p_Mass, m_MassCutoffs(MHeF)) < 0)
            ? (m_Bn[29] * PPOW(p_Mass, m_Bn[30])) / (1.0 + (m_Alpha3 * exp(15.0 * (p_Mass - m_MassCutoffs(MHeF)))))
            : (m_Bn[31] + (m_Bn[32] * PPOW(p_Mass, (m_Bn[33] + 1.8)))) / (m_Bn[34] + PPOW(p_Mass, m_Bn[33]));
}


/*
 * CalculateMassLossRateWolfRayetShenar2019
 *
 * Calculate the mass-loss rate for Wolf--Rayet stars according to the
 * prescription of Shenar et al. 2019 (https://ui.adsabs.harvard.edu/abs/2019A%26A...627A.151S/abstract)
 * 
 * See their Eq. 6 and Table 5
 * 
 * We use the fitting coefficients for hydrogen rich WR stars (e.g., WNh)
 * The C4 (X_He) term is = 0 and is omitted
 *
 * Uses current values of:
 * 
 *    - m_Log10Metallicity
 *    - m_Luminosity
 *    - m_Temperature
 *  
 * 
 * double CalculateMassLossRateWolfRayetShenar2019()
 *
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::WR)
 */
std::tuple<double, MASS_LOSS_TYPE> BaseStar::CalculateMassLossRateWolfRayetShenar2019() const {

    // For H-rich WR stars (X_H > 0.4)
    constexpr double C1 = -6.78;
    constexpr double C2 =  0.66;
    constexpr double C3 = -0.12;
    constexpr double C5 =  0.74;

    return std::make_tuple(PPOW(10.0, C1 + (C2 * log10(m_Luminosity)) + (C3 * log10(m_Temperature * TSOL)) + (C5 * m_Log10Metallicity)), MASS_LOSS_TYPE::WR);
}


/*
 * CalculateMassLossRateRSGYang2023
 *
 * Calculate mass loss rate for RSG stars (Red Supergiant) per Yang 2023
 * https://arxiv.org/pdf/2303.09385.pdf eq 6.
 *
 *
 * Uses current values of:
 * 
 *    - m_Luminosity
 * 
 *  
 * std::tuple<double, MASS_LOSS_TYPE> CalculateMassLossRateRSGYang2023() const
 *
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
std::tuple<double, MASS_LOSS_TYPE> BaseStar::CalculateMassLossRateRSGYang2023() const {

    double logL   = log10(m_Luminosity);
    double logL_2 = logL * logL;

    return std::make_tuple(PPOW(10.0, 0.45 * logL_2 * logL - 5.26 * logL_2 + 20.93 * logL - 34.56), MASS_LOSS_TYPE::RSG);
}


/*
 * CalculateMassLossRateVMSBestenlehner2020
 *
 * Calculate mass loss rate for very massive (>100 Msol) OB stars per Bestenlehner 2020
 * https://arxiv.org/pdf/2002.05168.pdf
 *
 * Uses current values of:
 * 
 *    - m_Luminosity
 *    - m_Mass
 * 
 * 
 * std::tuple<double, MASS_LOSS_TYPE> CalculateMassLossRateVMSBestenlehner2020() const
 *
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for very massive stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::VMS)
 */
std::tuple<double, MASS_LOSS_TYPE> BaseStar::CalculateMassLossRateVMSBestenlehner2020() const {

    constexpr double alpha       = 0.39;                                // CAK force multiplier
    constexpr double logMdotZero = -4.78;                               // from substituting LogMdotTrans and Gamma_e trans into eq 12. 

    double gamma = EDDINGTON_PARAMETER_FACTOR * m_Luminosity / m_Mass;  // Eddington Parameter, not metallicity specific as in the publication

    return std::make_tuple(PPOW(10.0, logMdotZero + ((1.0 / alpha) + 0.5) * log10(gamma) - (((1.0 - alpha) / alpha) + 2.0) * log10(1.0 - gamma)), MASS_LOSS_TYPE::VMS);
}


/*
 * IsOneOf
 *
 */
inline bool BaseStar::IsOneOf(const STELLAR_TYPE_LIST p_List) const {
    for (auto elem: p_List) {
        if (m_StellarType == elem) return true;
    }
    return false;
};

#endif // __BaseStar_h__
