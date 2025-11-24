#ifndef __BH_h__
#define __BH_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Remnants.h"

class BaseStar;
class NS;

class BH: virtual public BaseStar, public Remnants {
    
public:
    
    BH() { m_StellarType = STELLAR_TYPE::BLACK_HOLE; };
    
    BH(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::BLACK_HOLE;                                                                                                       // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                 // Initialise if required
    }
    
    
    // member functions - alphabetically
inline static double CalculateLuminosityOnPhase_Hurley2000() const override { return CalculateLuminosityOnPhase_Hurley2000_Static(); }
inline static double CalculateLuminosityOnPhase_Hurley2000_Static() { return 1.0E-10; } // Hurley et al. 2000, eq 96



static  double      CalculateNeutrinoMassLoss_Static(const double p_BaryonicMass);




COMPAS_PURE static double CalculateSNkickWeighting_Static(const double p_Mass, const double p_FallbackFraction);
   
    
protected:
    
    void Initialise() {

        // set internal properties to zero to avoid meaningless values
        m_Age        = 0.0;
        m_COCoreMass = 0.0;
        m_HeCoreMass = 0.0;
        m_CoreMass   = 0.0;
        m_Mass0      = 0.0;
        
        EvolveOnPhase(0.0);
    }
    
    
    // member functions - alphabetically
    double  CalculateConvergedMassStepZetaNuclear() const                           { return 0.0; }
    double  CalculateEddingtonCriticalRate() const                                  { return 2.6E-8 * m_Mass * MYR_TO_YEAR; }                           // E.g., Marchant+, 2017, Eq. 3, assuming accretion efficiency of 10%
    double  CalculateLuminosityOnPhase() const                                      { return CalculateLuminosityOnPhase_Static(); }

    



///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline MASS_LOSS_T CalculateMassLossRate() const override { return std::make_tuple(0.0, MASS_LOSS_TYPE::NONE); } // Ensure that NSs don't lose mass in winds

inline double CalculateRadius_Hurley2000() const override { return CalculateRadius_Hurley2000_Static(Mass()); }

GNU_CONST static inline double CalculateRadius_Hurley2000_Static(const double p_Mass) { return CalculateSchwarzschildRadius_Static(p_Mass); } // Hurley et al. 2000, eq 94

GNU_CONST static inline double CalculateSchwarzschildRadius_Static(const double p_Mass) { return 4.24E-6 * p_Mass; } // Schwarzschild radius of black hole

inline double CalculateMomentOfInertia() const override { return (2.0 / 5.0) * Mass() * Radius() * Radius(); } // MoI for solid sphere *ILYA* JR: that's not really right, is it?





///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




    
};

#endif // __BH_h__
