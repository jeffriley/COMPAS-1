#ifndef __NS_h__
#define __NS_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Rand.h"
#include "Remnants.h"
#include "BH.h"


class BaseStar;
class Remnants;

class NS: virtual public BaseStar, public Remnants {

public:

    NS() {
        m_StellarType = STELLAR_TYPE::NEUTRON_STAR;                                                                                                             // Set stellar type

        // set NS values based on options provided (so we don't need to do this every timestep)
        // JR: these are good candidates for the "globals/constants" singleton...
        NS_MAG_FIELD_LOWER_LIMIT = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField());
        NS_DECAY_MASS_SCALE      = OPTIONS->PulsarMagneticFieldDecayMassscale() * MSOL_TO_G;
        NS_DECAY_TIME_SCALE      = OPTIONS->PulsarMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR; 
    };
    
    NS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::NEUTRON_STAR;                                                                                                             // Set stellar type

        // set NS values based on options provided (so we don't need to do this every timestep)
        // JR: these are good candidates for the "globals/constants" singleton...
        NS_MAG_FIELD_LOWER_LIMIT = PPOW(10.0, OPTIONS->PulsarLog10MinimumMagneticField());
        NS_DECAY_MASS_SCALE      = OPTIONS->PulsarMagneticFieldDecayMassscale() * MSOL_TO_G;
        NS_DECAY_TIME_SCALE      = OPTIONS->PulsarMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR; 

        if (p_Initialise) Initialise();                                                                                                                         // Initialise if required
    }



   
    inline double CalculateLuminosity_Hurley2000() const override { return CalculateLuminosity_Hurley2000_Static(Mass(), Age()); }
    GNU_CONST inline static double CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Age) {
        const double t = std::max(p_Age, 0.1);
        return 0.02 * PPOW(p_Mass, 2.0 / 3.0) / (t * t); // Hurley et al. 2000, eq 93
    }


    inline double CalculateRadius() const override { return CalculateRadiusInKM_Static(Mass() * KM_TO_RSOL); }  // MSol
    COMPAS_PURE double CalculateRadiusInKM_Static(const double p_Mass);                                         // km





    static double       DeltaJByAccretion_Static(const double p_Mass, const double p_Radius_6, const double p_MagField, const double p_SpinFrequency, const double p_mDot, const double p_Epsilon);

    MT_CASE             DetermineMassTransferTypeAsDonor() const                { return MT_CASE::NONE; }                                                       // Always NONE


protected:
    
    void Initialise() {
        
        // set internal properties to zero to avoid meaningless values
        m_Age        = 0.0;
        m_COCoreMass = 0.0;
        m_HeCoreMass = 0.0;
        m_CoreMass   = 0.0;
        m_Mass0      = 0.0;
        
        EvolveOnPhase(0.0);

        CalculateAndSetPulsarParameters();
    }

    // member variables

    double m_AngularMomentum_CGS;                                                                                                                               // Current angular momentum in CGS - only required in NS class
    double m_MomentOfInertia_CGS;                                                                                                                               // MoI in CGS - only required in NS class


    // member functions - alphabetically

inline DBL_VECTOR CalculateTimescales_Hurley2000() const override { return m_InterimState.Timescales(); }     // not meaningful for NS and BH  <<<<<<<<<<<<<<<<<<<<<<<< IS THIS NECESSARY????


            void            CalculateAndSetPulsarParameters();

            double          CalculateBirthMagneticField();
            double          CalculateBirthSpinPeriod();
    
 double CalculateCriticalMassRatio_Hurley2002() const { return 0.0; }

    
    static  double          CalculateMomentOfInertiaCGS_Static(const double p_Mass, const double p_Radius);                                                     // MoI in CGS            
            double          CalculateMomentOfInertiaCGS() const                 { return CalculateMomentOfInertiaCGS_Static(m_Mass * MSOL_TO_G, m_Radius * RSOL_TO_CM); } // MOI in CGS - use member variables
            double          CalculateMomentOfInertia() const                    { return CalculateMomentOfInertiaCGS() / MSOL_TO_G / RSOL_TO_CM / RSOL_TO_CM; } // MoI (default is solar units)

    
            double          CalculateRadiusOnPhase(double p_Mass, double p_Luminosity) const { return CalculateRadiusOnPhase(); }                               // not a meaningful calculation for NS, ignore arguments

            double          CalculateSpinDownRate(const double p_Omega, const double p_MomentOfInteria, const double p_MagField, const double p_Radius) const;
  

            GNU_CONST inline double ChooseTimestep() const override { return ChooseTimestep(m_InterimState.Age()); }
            double ChooseTimestep(const double p_Age) const;

            STELLAR_TYPE    EvolveToNextPhase()                                 { return STELLAR_TYPE::BLACK_HOLE; }

            double          ResolveCommonEnvelopeAccretion(const double p_FinalMass,
                                                           const double p_CompanionMass     = 0.0,
                                                           const double p_CompanionRadius   = 0.0,
                                                           const double p_CompanionEnvelope = 0.0);
    
            bool            ShouldEvolveOnPhase() const                         { return (m_Mass <= OPTIONS->MaximumNeutronStarMass()); }                       // Evolve as a neutron star unless mass > maximum neutron star mass (e.g. through accretion)
            void            SpinDownIsolatedPulsar(const double p_Stepsize);
            void            UpdateMagneticFieldAndSpin(const bool   p_CommonEnvelope,
                                                       const bool   p_RecycledNS,
                                                       const double p_Stepsize,
                                                       const double p_MassGainPerTimeStep,
                                                       const double p_Epsilon);






///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline MASS_LOSS_T CalculateMassLossRate() const override { return std::make_tuple(0.0, MASS_LOSS_TYPE::NONE); } // Ensure that NSs don't lose mass in winds








///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



};

#endif // __NS_h__
