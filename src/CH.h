#ifndef __CH_h__
#define __CH_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"
#include "MS_gt_07.h"
#include "HeMS.h"

class BaseStar;
class MS_gt_07;

class CH: virtual public BaseStar, public MS_gt_07 {

public:

    CH() { m_StellarType = STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS; };
    
    CH(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MS_gt_07(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS;                                                                                                                       // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                             // Initialise if required
    }

    CH* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        CH* clone = new CH(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static CH* Clone(CH& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        CH* clone = new CH(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                                                      // Initialise timescales
        m_Age = 0.0;                                                                                                                                                                // Set age appropriately
        m_CHE = true;                                                                                                                                                               // initially for CH stars                                                                                                                                                            // Set age appropriately
    }

    // member functions

    // Abundances
    double          CalculateHeliumAbundanceCoreOnPhase(const double p_Tau) const;
    double          CalculateHeliumAbundanceCoreOnPhase() const         { return CalculateHeliumAbundanceCoreOnPhase(m_Tau); };
    double          CalculateHeliumAbundanceSurfaceOnPhase(const double p_Tau) const;
    double          CalculateHeliumAbundanceSurfaceOnPhase() const      { return CalculateHeliumAbundanceSurfaceOnPhase(m_Tau); };

    double          CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const;
    double          CalculateHydrogenAbundanceCoreOnPhase() const       { return CalculateHydrogenAbundanceCoreOnPhase(m_Tau); };
    double          CalculateHydrogenAbundanceSurfaceOnPhase(const double p_Tau) const;
    double          CalculateHydrogenAbundanceSurfaceOnPhase() const    { return CalculateHydrogenAbundanceSurfaceOnPhase(m_Tau); };
    
    // Lifetime
    double          CalculateLogLifetimeRatio(const double p_Mass) const;
    double          CalculateLifetimeRatio(const double p_Mass) const;

    // Luminosity
    double          CalculateLogLuminosityRatio(const double p_Mass, const double p_Tau) const;

    double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtPhaseEnd() const               { return CalculateLuminosityAtPhaseEnd(m_Mass0); }                                                          // Use class member variables


double CalculateLuminosityOnPhase() const; // Uses globals and state variables
double CalculateLuminosityOnPhase(m_StateHistory.CurrentState().Mass(),
                                  m_StateHistory.ZAMSState().Luminosity(),
                                  m_StateHistory.Timescales(),
                                  m_StateHistory.CurrentState().Time(),
                                  GLOBALS->ReferenceMetallicity(),
                                  GLOBALS->HurleyACoefficients(),
                                  GLOBALS->HurleyLuminosityConstants());

    
    // Mass loss rate
std::tuple<double, MASS_LOSS_TYPE> CalculateMLRateBelczynski2010(const double p_Metallicity, const double p_Luminosity, const double p_HeAbundanceSurface) const
std::tuple<double, MASS_LOSS_TYPE> CalculateMLRateMerritt2025(const double p_Metallicity, const double p_Luminosity, const double p_HeAbundanceSurface) const


    double          CalculateMassLossRateWeightOB(const double p_HeliumAbundanceSurface);
    
    // Radius
    double          CalculateRadiusOnPhase() const                      { return m_RZAMS; }                                                                                         // Constant from birth
    double          CalculateRadiusAtPhaseEnd() const                   { return CalculateRadiusOnPhase(); }                                                                        // Same as on phase

    // Timescales
DBL_VECTOR CalculateTimescales(const double p_Mass, const DBL_VECTOR& p_Timescales) const;


    void            CalculateTimescales()                               { return BaseStar::CalculateTimescales(); };  /// ?????????????????????????????????????????

    STELLAR_TYPE    EvolveToNextPhase();

    bool            ShouldEvolveOnPhase() const                         { return m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)] && (OPTIONS->OptimisticCHE() || Omega() >= m_OmegaCHE); } // Evolve on CHE phase if age in MS timescale and spinning at least as fast as CHE threshold

    void            UpdateAgeAfterMassLoss();
    
    void            UpdateMainSequenceCoreMass(const double p_Dt, const double p_TotalMassLossRate)    { };                                                                         // Do not use core mass calculations during CHE phase

};




/*
 * CalculateHeAbundanceCore
 *
 * @brief
 * Calculate the helium abundance in the core of the star, given the fractional
 * MS age of the star.  Currently just a simple linear model from the initial
 * helium abundance to the maximum helium abundance (assuming that all hydrogen
 * is converted to helium). 
 * 
 * Should one day be updated to match detailed models.
 * 
 *
 * double CalculateHeAbundanceCore(const double p_Tau, const double p_InitialHeAbundance) const
 * 
 * @param       p_Tau                           Fractional MS age of the star
 * @param       p_InitialHeAbundance            Initial helium abundance of the star
 * @return                                      Helium abundance in the core of the star
 */
inline double CH::CalculateHeAbundanceCore(const double p_Tau, const double p_InitialHeAbundance) const {
    return ((1.0 - m_Metallicity - p_InitialHeAbundance) * p_Tau) + p_InitialHeAbundance;
}


/*
 * CalculateHeAbundanceSurface
 *
 * @brief
 * Calculate the helium abundance at the surface of the star, given the fractional
 * MS age of the star.  Since the star is checmically homogeneous, the core and the
 * surface have the same abundances.
 * 
 * 
 * double CalculateHeAbundanceSurface(const double p_Tau, const double p_InitialHeAbundance) const
 * 
 * @param       p_Tau                           Fractional MS age of the star
 * @param       p_InitialHeAbundance            Initial helium abundance of the star
 * @return                                      Helium abundance at the surface of the star
 */
inline double CH::CalculateHeAbundanceSurface(const double p_Tau, const double p_InitialHeAbundance) const {
    return CalculateHeAbundanceCore(p_Tau, p_InitialHeAbundance);
}


/*
 * CalculateHAbundanceCore
 *
 * @brief
 * Calculate the hydrogen abundance in the core of the star, given the fractional MS age
 * of the star.  Currently just a simple linear model that assumes that hydrogen in the
 * core of the star is burned to helium at a constant rate throughout the star's lifetime.
 * 
 * Should one day be updated to match detailed models.
 * 
 *
 * double CalculateHAbundanceCore(const double p_Tau, const double p_InitialHAbundance) const
 * 
 * @param       p_Tau                           Fractional MS age of the star
 * @param       p_InitialHAbundance             Initial hydrogen abundance of the star
 * @return                                      Hydrogen abundance in the core of the star
 */
inline double CH::CalculateHAbundanceCore(const double p_Tau, const double p_InitialHAbundance) const {
    return p_InitialHAbundance * (1.0 - p_Tau);
}


/*
 * CalculateHAbundanceSurface
 *
 * @brief
 * Calculate the hydrogen abundance at the surface of the star, given the fractional
 * MS age of the star.  Since the star is chemically homogeneous, the core and the
 * surface have the same abundances.
 * 
 *
 * double CalculateHAbundanceSurface(const double p_Tau, const double p_InitialHeAbundance) const
 * 
 * @param       p_Tau                           Fractional MS age of the star
 * @param       p_InitialHAbundance             Initial hydrogen abundance of the star
 * @return                                      Hydrogen abundance at the surface of the star
 */
inline double CH::CalculateHAbundanceSurface(const double p_Tau, const double p_InitialHAbundance) const {
    return CalculateHAbundanceCore(p_Tau, p_InitialHAbundance);
}


/*
 * CalculateLifetimesRatio
 *
 * @brief
 * Calculate the ratio of the lifetimes of a CH star and a non-CH MS star of the same mass
 * (t_CHE/t_MS).  Uses a polynomial fit derived using IZw18 and IZw18CHE BoOST models from
 * Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203).
 *
 * 
 * double CalculateLifetimesRatio(const double p_Mass) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Lifetimes ratio
 */
inline double CH::CalculateLifetimesRatio(const double p_Mass) const {

    const double x  = log10(p_Mass);
    const double x2 = x * x;

    return PPOW(10.0, -0.15929168474199387 + (1.050069750483549 * x) + (-0.8233601359988406 * x2) + (0.17772610259473764 * x * x2));
}


/*
 * CalculateLogLuminositiesRatio
 *
 * @brief
 * Calculate the ratio of the log luminosities of a CH star and a non-CH MS of the same mass
 * (log(L_CH) / log(L_MS)).  Uses a polynomial fit derived using IZw18 and IZw18CHE BoOST 
 * models from Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203).
 *
 * 
 * double CalculateLogLuminositiesRatio(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Log luminosities ratio
 */
inline double CH::CalculateLogLuminositiesRatio(const double p_Mass) const {
         
    const double x  = log10(p_Mass);
    const double x2 = x * x;

    return 1.8261540986808193 + (-1.1636822407341694  * x) + (0.5876329884434304  * x2) + (-0.10236336828026288 * x * x2);
}


/*
 * CalculateMassLossFractionOB
 *
 * @brief
 * Calculate the fraction of mass loss attributable to OB mass loss, per Yoon et al. 2006
 * 
 * The model described in Yoon et al. 2006 (also Szecsi et al. 2015) uses OB mass loss while the
 * He surface abundance is below 0.55, WR mass loss when the surface He abundance is above 0.7,
 * and linearly interpolate when the He surface abundance is between those limits.
 * 
 * This function calculates the fraction of mass loss attributable to OB mass loss, based on
 * the He surface abundance and the abundance limits described in Yoon et al. 2006.  The value
 * returned will be 1.0 if 100% of the mass loss is attributable to OB mass lass, 0.0 if 100% of
 * the mass loss is attributable to WR mass loss, and in the range (0.0, 1.0) if the mass loss is
 * a mix of OB and WR.
 * 
 *
 * double CalculateMassLossFractionOB(const double p_HeAbundanceSurface) const
 *
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Fraction of mass loss attributable to OB mass loss
 */
inline double CH::CalculateMassLossFractionOB(const double p_HeAbundanceSurface) const {

    constexpr double limOB = 0.55;      // per Yoon et al. 2006
    constexpr double limWR = 0.70;      // per Yoon et al. 2006

    return p_HeAbundanceSurface <= limOB ? 1.0 : (p_HeAbundanceSurface >= limWR ? 0.0 : (limWR - p_HeAbundanceSurface) / (limWR - limOB));
}

#endif // __CH_h__
