#ifndef __Star_h__
#define __Star_h__

#include <fstream>

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"
#include "MS_lte_07.h"
#include "MS_gt_07.h"
#include "CH.h"
#include "HG.h"
#include "FGB.h"
#include "CHeB.h"
#include "EAGB.h"
#include "TPAGB.h"
#include "HeMS.h"
#include "HeHG.h"
#include "HeGB.h"
#include "HeWD.h"
#include "COWD.h"
#include "ONeWD.h"
#include "NS.h"
#include "BH.h"
#include "MR.h"

class BaseStar;
class MS_lte_07;
class MS_gt_07;
class CH;
class HG;
class FGB;
class CHeB;
class EAGB;
class TPAGB;
class HeMS;
class HeHG;
class HeGB;
class HeWD;
class COWD;
class ONeWD;
class NS;
class BH;
class MR;


class Star {

public:

//    Star();

    Star(const STARTING_STELLAR_TYPE p_StartingStellarType,
         const unsigned long int     p_RandomSeed, 
         const double                p_Metallicity, 
         const double                p_Mass, 
         const KickParameters        p_KickParameters,
         const double                p_RotationalFrequency = -1.0); 

    Star(const Star& p_Star);

    virtual ~Star() { delete m_Star; }


    // object identifiers - all classes have these
    inline OBJECT_ID            ObjectId() const                                                            { return m_ObjectId; }
    inline OBJECT_TYPE          ObjectType() const                                                          { return OBJECT_TYPE::STAR; }
    inline OBJECT_PERSISTENCE   ObjectPersistence() const                                                   { return m_ObjectPersistence; }



    inline STELLAR_TYPE         StartingStellarType() const                                                 { return m_Star->StartingStellarType(); }
    inline STELLAR_TYPE         StellarType() const                                                         { return m_Star->StellarType(); }

    // getters - alphabetically
    inline double               Age() const                                                                 { return m_Star->Age(); }
    inline double               AngularMomentum() const                                                     { return m_Star->AngularMomentum(); }
    inline double               CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate, 
                                                   const double p_massTransferEfficiencyBeta) const         { return m_Star->CalculateCriticalMassRatio(p_AccretorIsDegenerate, p_massTransferEfficiencyBeta); }
    inline double               CalculateCriticalMassRatioClaeys14(const bool p_AccretorIsDegenerate) const { return m_Star->CalculateCriticalMassRatioClaeys14(p_AccretorIsDegenerate); }
    inline double               CalculateCriticalMassRatioHurleyHjellmingWebbink() const                    { return m_Star->CalculateCriticalMassRatioHurleyHjellmingWebbink(); }
    inline double               CalculateDynamicalTimescale() const                                         { return m_Star->CalculateDynamicalTimescale(); }
    inline double               CalculateLambdaLoveridge() const                                            { return m_Star->CalculateLambdaLoveridge(); }
    inline double               CalculateLambdaNanjing() const                                              { return m_Star->CalculateLambdaNanjing(); }
    inline double               CalculateLambdaKruckow() const                                              { return m_Star->CalculateLambdaKruckow(); }
    inline double               CalculateLambdaDewi() const                                                 { return m_Star->CalculateLambdaDewi(); }
    inline double               CalculateRadialExpansionTimescale() const                                   { return m_Star->CalculateRadialExpansionTimescale(); }
    inline double               CalculateThermalTimescale() const                                           { return m_Star->CalculateThermalTimescale(); }
    inline double               COCoreMass() const                                                          { return m_Star->COCoreMass(); }
    inline double               CoreMass() const                                                            { return m_Star->CoreMass(); }
    inline double               Dt() const                                                                  { return m_Star->Dt(); }
    inline bool                 EnvelopeJustExpelledByPulsations() const                                    { return m_Star->EnvelopeJustExpelledByPulsations(); }
    inline bool                 ExperiencedAIC() const                                                      { return m_Star->ExperiencedAIC(); }
    inline bool                 ExperiencedCCSN() const                                                     { return m_Star->ExperiencedCCSN(); }
    inline bool                 ExperiencedECSN() const                                                     { return m_Star->ExperiencedECSN(); }
    inline bool                 ExperiencedPISN() const                                                     { return m_Star->ExperiencedPISN() ; }
    inline bool                 ExperiencedPPISN() const                                                    { return m_Star->ExperiencedPPISN(); }
    inline bool                 ExperiencedUSSN() const                                                     { return m_Star->ExperiencedUSSN(); }
    inline double               HeCoreMass() const                                                          { return m_Star->HeCoreMass(); }
    inline double               HeliumAbundanceCore() const                                                 { return m_Star->HeliumAbundanceCore(); }
    inline double               HeliumAbundanceSurface() const                                              { return m_Star->HeliumAbundanceSurface();} 
    inline double               HydrogenAbundanceCore() const                                               { return m_Star->HydrogenAbundanceCore(); }
    inline double               HydrogenAbundanceSurface() const                                            { return m_Star->HydrogenAbundanceSurface(); }
    inline double               InitialHeliumAbundance() const                                              { return m_Star->InitialHeliumAbundance(); }
    inline double               InitialHydrogenAbundance() const                                            { return m_Star->InitialHydrogenAbundance(); }
    inline double               InitialMainSequenceCoreMass() const                                         { return m_Star->InitialMainSequenceCoreMass(); }
    inline bool                 IsAIC() const                                                               { return m_Star->IsAIC(); }
    inline bool                 IsCCSN() const                                                              { return m_Star->IsCCSN(); }
    inline bool                 IsDegenerate() const                                                        { return m_Star->IsDegenerate(); }
    inline bool                 IsECSN() const                                                              { return m_Star->IsECSN(); }
    inline bool                 IsHeSD() const                                                              { return m_Star->IsHeSD(); }
    inline bool                 IsOneOf(STELLAR_TYPE_LIST p_List) const                                     { return m_Star->IsOneOf(p_List); }
    inline bool                 IsPISN() const                                                              { return m_Star->IsPISN(); }
    inline bool                 IsPPISN() const                                                             { return m_Star->IsPPISN(); }
    inline bool                 IsSNIA() const                                                              { return m_Star->IsSNIA(); }
    inline bool                 IsSupernova() const                                                         { return m_Star->IsSupernova(); }
    inline bool                 IsUSSN() const                                                              { return m_Star->IsUSSN(); }
    inline double               Luminosity() const                                                          { return m_Star->Luminosity(); }
    inline double               MainSequenceCoreMass() const                                                { return m_Star->MainSequenceCoreMass(); }
    inline double               Mass() const                                                                { return m_Star->Mass(); }
    inline double               Mass0() const                                                               { return m_Star->Mass0(); }
    inline double               MassPrev() const                                                            { return m_Star->MassPrev(); }
    inline double               Metallicity() const                                                         { return m_Star->Metallicity(); }
    inline double               MZAMS() const                                                               { return m_Star->MZAMS(); }
    inline double               Omega() const                                                               { return m_Star->Omega(); }
    inline double               OmegaBreak() const                                                          { return m_Star->OmegaBreak(); }
    inline double               OmegaCHE() const                                                            { return m_Star->OmegaCHE(); }
    inline double               Radius() const                                                              { return m_Star->Radius(); }
    inline double               RadiusPrev() const                                                          { return m_Star->RadiusPrev(); }
    inline unsigned long int    RandomSeed() const                                                          { return m_Star->RandomSeed(); }
    inline double               RZAMS() const                                                               { return m_Star->RZAMS(); }
    inline SupernovaDetailsT    SN_Details() const                                                          { return m_Star->SN_Details(); }
    inline double               SN_EccentricAnomaly() const                                                 { return m_Star->SN_EccentricAnomaly(); }
    inline double               SN_KickMagnitude()                                                          { return m_Star->SN_KickMagnitude() ; }
    inline double               SN_Phi() const                                                              { return m_Star->SN_Phi(); }
    inline double               SN_RocketKickMagnitude()                                                    { return m_Star->SN_RocketKickMagnitude(); }
    inline double               SN_RocketKickPhi()                                                          { return m_Star->SN_RocketKickPhi(); }
    inline double               SN_RocketKickTheta()                                                        { return m_Star->SN_RocketKickTheta(); }
    inline double               SN_Theta() const                                                            { return m_Star->SN_Theta(); }
    inline double               SN_TotalMassAtCOFormation() const                                           { return m_Star->SN_TotalMassAtCOFormation(); }
    inline double               SN_TotalRadiusAtCOFormation() const                                         { return m_Star->SN_TotalRadiusAtCOFormation(); }
    inline double               SN_TrueAnomaly() const                                                      { return m_Star->SN_TrueAnomaly(); }
    inline SN_EVENT             SN_Type() const                                                             { return m_Star->SN_Type(); }
    inline double               Speed() const                                                               { return m_Star->Speed(); }
    inline COMPAS_VARIABLE      StellarPropertyValue(const T_ANY_PROPERTY p_Property) const                 { return m_Star->StellarPropertyValue(p_Property); }
    inline STELLAR_TYPE         StellarTypePrev() const                                                     { return m_Star->StellarTypePrev(); }
    inline double               Tau() const                                                                 { return m_Star->Tau(); }
    inline double               Temperature() const                                                         { return m_Star->Temperature(); }
    inline double               Timescale(TIMESCALE p_Timescale) const                                      { return m_Star->Timescale(p_Timescale); }
    inline double               TotalMassLossRate() const                                                   { return m_Star->TotalMassLossRate(); }    
    inline ACCRETION_REGIME     WhiteDwarfAccretionRegime() const                                           { return m_Star->WhiteDwarfAccretionRegime(); }
    inline double               XExponent() const                                                           { return m_Star->XExponent(); }

    
    // setters
    inline void                 SetObjectId(const OBJECT_ID p_ObjectId)                                     { m_ObjectId = p_ObjectId; }
    inline void                 SetPersistence(const OBJECT_PERSISTENCE p_Persistence)                      { m_ObjectPersistence = p_Persistence; }

    inline void                 ClearCurrentSNEvent()                                                       { m_Star->ClearCurrentSNEvent(); }

    inline void                 ResetEnvelopeExpulsationByPulsations()                                      { m_Star->ResetEnvelopeExpulsationByPulsations(); }

    inline void                 SetAngularMomentum(double p_AngularMomentum)                                { m_Star->SetAngularMomentum(p_AngularMomentum); }
    inline void                 SetDt(const double p_dt)                                                    { m_Star->SetDt(p_dt); }
    inline void                 SetOmega(double p_Omega)                                                    { m_Star->SetOmega(p_Omega); }
    inline void                 SetPrevDt(const double p_dt)                                                { m_Star->SetPrevDt(p_dt); }
    inline void                 SetSNCurrentEvent(const SN_EVENT p_SNEvent)                                 { m_Star->SetSNCurrentEvent(p_SNEvent); }
    inline void                 SetSNPastEvent(const SN_EVENT p_SNEvent)                                    { m_Star->SetSNPastEvent(p_SNEvent); }

    inline void                 UpdateMassTransferDonorHistory()                                            { m_Star->UpdateMassTransferDonorHistory(); }


    // member functions - alphabetically
    STELLAR_TYPE                AgeOneTimestep(const double p_Dt, bool p_Switch = true);

    inline void                 ApplyMassTransferRejuvenationFactor()                                                               { m_Star->ApplyMassTransferRejuvenationFactor(); }

    inline double               CalculateBindingEnergy(const double p_CoreMass, const double p_EnvMass, const double p_Radius, const double p_Lambda) const {
                                    return m_Star->CalculateBindingEnergy(p_CoreMass, p_EnvMass, p_Radius, p_Lambda);
                                }

    inline double               CalculateConvectiveCoreMass()                                                                       { return m_Star->CalculateConvectiveCoreMass(); }
    inline double               CalculateConvectiveCoreRadius()                                                                     { return m_Star->CalculateConvectiveCoreRadius(); }

    inline double               CalculateConvectiveEnvelopeBindingEnergy(const double p_TotalMass, const double p_ConvectiveEnvelopeMass, const double p_Radius, const double p_Lambda) {
                                    return m_Star->CalculateConvectiveEnvelopeBindingEnergy(p_TotalMass, p_ConvectiveEnvelopeMass, p_Radius, p_Lambda);
                                }

    inline double               CalculateConvectiveEnvelopeLambdaPicker(const DBL_DBL p_convectiveEnvelopeMass) const               { return m_Star->CalculateConvectiveEnvelopeLambdaPicker(p_convectiveEnvelopeMass); }
    
    inline DBL_DBL              CalculateConvectiveEnvelopeMass()                                                                   { return m_Star->CalculateConvectiveEnvelopeMass(); }
    
    inline double               CalculateEddyTurnoverTimescale()                                                                    { return m_Star->CalculateEddyTurnoverTimescale(); }
    
    inline DBL_DBL_DBL_DBL      CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)      { return m_Star->CalculateImKnmDynamical(p_Omega, p_SemiMajorAxis, p_M2); }
    inline DBL_DBL_DBL_DBL      CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)    { return m_Star->CalculateImKnmEquilibrium(p_Omega, p_SemiMajorAxis, p_M2); }
    inline DBL_DBL_DBL_DBL      CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)          { return m_Star->CalculateImKnmTidal(p_Omega, p_SemiMajorAxis, p_M2); }

    inline DBL_DBL              CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate, const bool p_IsHeRich) {
                                    return m_Star->CalculateMassAcceptanceRate(p_DonorMassRate, p_AccretorMassRate, p_IsHeRich);
                                }

    inline double               CalculateMassLossValues(double p_Dt, const bool p_UpdateMDot = false)                               { return m_Star->CalculateMassLossValues(p_Dt, p_UpdateMDot); }

    inline double               CalculateMomentOfInertia() const                                                                    { return m_Star->CalculateMomentOfInertia(); }
    inline double               CalculateMomentOfInertiaAU() const                                                                  { return m_Star->CalculateMomentOfInertiaAU(); }
    
    inline double               CalculateNuclearMassLossRate()                                                                      { return m_Star->CalculateNuclearMassLossRate(); }
    
    inline double               CalculateRadialExpansionTimescaleDuringMassTransfer()                                               { return m_Star->CalculateRadialExpansionTimescaleDuringMassTransfer(); }
    
    inline double               CalculateRadialExtentConvectiveEnvelope()                                                           { return m_Star->CalculateRadialExtentConvectiveEnvelope(); }

    inline double               CalculateRadiusOnMassChange(double p_dM)                                                            { return m_Star->CalculateRadiusOnMassChange(p_dM); }
    
    inline double               CalculateRemnantRadius()                                                                            { return m_Star->CalculateRemnantRadius(); }
    
    inline void                 CalculateSNAnomalies(const double p_Eccentricity)                                                   { m_Star->CalculateSNAnomalies(p_Eccentricity); }
    
    inline double               CalculateSNKickMagnitude(const double p_RemnantMass, const double p_EjectaMass, const STELLAR_TYPE p_StellarType) {
                                    return m_Star->CalculateSNKickMagnitude(p_RemnantMass, p_EjectaMass, p_StellarType);
                                }

    inline double               CalculateThermalMassAcceptanceRate(const double p_Radius)                                           { return m_Star->CalculateThermalMassAcceptanceRate(p_Radius); }
    inline double               CalculateThermalMassAcceptanceRate()                                                                { return m_Star->CalculateThermalMassAcceptanceRate(); }

    inline double               CalculateThermalMassLossRate() const                                                                { return m_Star->CalculateThermalMassLossRate(); }

    inline double               CalculateThermalTimescale(const double p_Radius) const                                              { return m_Star->CalculateThermalTimescale(p_Radius); }

    inline double               CalculateTimestep()                                                                                 { return m_Star->CalculateTimestep(); }

    inline double               CalculateZetaAdiabatic() const                                                                      { return m_Star->CalculateZetaAdiabatic(); }
    inline double               CalculateZetaConstantsByEnvelope(ZETA_PRESCRIPTION p_ZetaPrescription)                              { return m_Star->CalculateZetaConstantsByEnvelope(p_ZetaPrescription); }
    inline double               CalculateZetaEquilibrium()                                                                          { return m_Star->CalculateZetaEquilibrium(); }

    inline ACCRETION_REGIME     DetermineAccretionRegime(const double p_DonorThermalMassLossRate, const bool p_HeRich)              { return m_Star->DetermineAccretionRegime(p_DonorThermalMassLossRate, p_HeRich); }

    inline ENVELOPE             DetermineEnvelopeType() const                                                                       { return m_Star->DetermineEnvelopeType(); }

    EVOLUTION_STATUS            Evolve(const long int p_Id);

    STELLAR_TYPE                EvolveOneTimestep(const double p_Dt);

    inline double               InterpolateGeEtAlQCrit(const QCRIT_PRESCRIPTION p_qCritPrescription, const double p_massTransferEfficiencyBeta) {
                                    return m_Star->InterpolateGeEtAlQCrit(p_qCritPrescription, p_massTransferEfficiencyBeta);
                                }

    inline void                 HaltWinds()                                                                                         { m_Star->HaltWinds(); }

    inline void                 ResolveAccretion(const double p_AccretionMass)                                                      { m_Star->ResolveAccretion(p_AccretionMass); }

    inline void                 ResolveAccretionRegime(const ACCRETION_REGIME p_Regime, const double p_DonorThermalMassLossRate)    { m_Star->ResolveAccretionRegime(p_Regime, p_DonorThermalMassLossRate); }
    
    inline double               ResolveCommonEnvelopeAccretion(const double p_FinalMass, const double p_CompanionMass, const double p_CompanionRadius, const double p_CompanionEnvelope) {
                                    return m_Star->ResolveCommonEnvelopeAccretion(p_FinalMass, p_CompanionMass, p_CompanionRadius, p_CompanionEnvelope);
                                } 

    inline STELLAR_TYPE         ResolveEnvelopeLoss(bool p_Force = false)                                                           { return m_Star->ResolveEnvelopeLoss(p_Force); } 
    void                        ResolveEnvelopeLossAndSwitch();

    inline void                 ResolveShellChange(const double p_AccretedMass)                                                     { m_Star->ResolveShellChange(p_AccretedMass); }

    bool                        RevertState();
    void                        SaveState();

    
    inline void                 SpinDownIsolatedPulsar(const double p_Stepsize)                                                     { m_Star->SpinDownIsolatedPulsar(p_Stepsize); }

    STELLAR_TYPE                SwitchTo(const STELLAR_TYPE p_StellarType, bool p_SetInitialType = false);
    
    inline double               CalculateTAMSCoreMass() const                                                                       { return m_Star->CalculateTAMSCoreMass(); }

    inline void                 UpdateAfterMerger(double p_Mass, double p_HydrogenMass)                                             { m_Star->UpdateAfterMerger(p_Mass, p_HydrogenMass); }
    inline void                 UpdateAgeAfterMassLoss()                                                                            { m_Star->UpdateAgeAfterMassLoss(); }

    inline void                 UpdateAttributes()                                                                                  { (void)UpdateAttributes(0.0, 0.0, true); }
    STELLAR_TYPE                UpdateAttributes(const double p_DeltaMass, const double p_DeltaMass0, const bool p_ForceRecalculate = false);

    STELLAR_TYPE                UpdateAttributesAndAgeOneTimestep(const double p_DeltaMass, const double p_DeltaMass0, const double p_DeltaTime, const bool p_Switch = true, const bool p_ForceRecalculate = false);

    inline void                 UpdateComponentVelocity(const Vector3d p_newVelocity)                                               { m_Star->UpdateComponentVelocity(p_newVelocity); }

    inline void                 UpdateInitialMass()                                                                                 { m_Star->UpdateInitialMass(); }

    inline void                 UpdateMagneticFieldAndSpin(const bool p_CommonEnvelope, const bool p_RecycledNS, const double p_Stepsize, const double p_MassGainPerTimeStep, const double p_Epsilon) {
                                    m_Star->UpdateMagneticFieldAndSpin(p_CommonEnvelope, p_RecycledNS, p_Stepsize, p_MassGainPerTimeStep, p_Epsilon);
                                }

    inline void                 UpdateMainSequenceCoreMass(const double p_Dt, const double p_TotalMassLossRate)                     { m_Star->UpdateMainSequenceCoreMass(p_Dt, p_TotalMassLossRate); }
    
    inline void                 UpdateTotalMassLossRate(const double p_MassLossRate)                                                { m_Star->UpdateTotalMassLossRate(p_MassLossRate); }


private:

    long int  m_Id;                             // id used to name output files - uses p_Id as passed (usually the step number of multiple single stars being produced)

    BaseStar *m_Star;                           // pointer to star

    std::vector<double> m_Timesteps;            // timesteps vector - for debugging/testing

    // thresholds flags for system snapshot file
    DBL_VECTOR  m_SystemSnapshotAgeFlags;
    BOOL_VECTOR m_SystemSnapshotTimeFlags;

protected:

    OBJECT_ID          m_ObjectId;              // instantiated object's unique object id
    OBJECT_PERSISTENCE m_ObjectPersistence;     // instantiated object's persistence

};

#endif // __Star_h__
