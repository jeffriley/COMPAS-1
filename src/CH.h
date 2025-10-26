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


///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_PURE double CalculateAgeAfterMassLoss_Hurley2000(const double      p_Age,
                                                     const double      p_ZetaHurley,
                                                     const double      p_tMS,
                                                     const DBL_VECTOR& p_aN,
                                                     const bool        p_EnhanceLifetime) const override;


inline double CalculateRadius() const override { return m_StateHistory.ZAMSState().Radius(); }              // CH radius is fixed throughout lifetime

GNU_CONST double CalculateHAbundanceCore(const double p_Tau, const double p_InitialHAbundance) const override;
GNU_CONST double CalculateHAbundanceSurface(const double p_Tau, const double p_InitialHAbundance) const override;
GNU_CONST double CalculateHeAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance = 0.0) const override;
GNU_CONST double CalculateHeAbundanceSurface(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance) const override;



GNU_CONST double CalculateLifetimesRatio(const double p_Mass) const;

GNU_CONST double CalculateLogLuminositiesRatio(const double p_Mass) const;

GNU_CONST double CalculateMLFractionOB(const double p_HeAbundanceSurface) const;

GNU_CONST DBL_VECTOR CalculateTimescales_Hurley(const double p_Mass, const DBL_VECTOR& p_tScales, const DBL_VECTOR& p_aN) const;


 



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inline double CalculateRadiusAtPhaseEnd() const override { return m_StateHistory.ZAMSState().Radius(); }    // CH radius is fixed throughout lifetime 







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
MASS_LOSS_T CalculateMLRate_Belczynski2010(const double                     p_Metallicity,
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
                                           const LBV_MASS_LOSS_PRESCRIPTION p_LBVprescription) const

COMPAS_PURE MASS_LOSS_T CH::CalculateMLRate_Merritt2025(const double p_Metallicity,
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
                                                        const double p_TerminalWindScalePower) const override;


    double          CalculateMassLossRateWeightOB(const double p_HeliumAbundanceSurface);

    
// Radius
                                                     // Same as on phase







    STELLAR_TYPE    EvolveToNextPhase();

    bool            ShouldEvolveOnPhase() const                         { return m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)] && (OPTIONS->OptimisticCHE() || Omega() >= m_OmegaCHE); } // Evolve on CHE phase if age in MS timescale and spinning at least as fast as CHE threshold











    void            UpdateMainSequenceCoreMass(const double p_Dt, const double p_TotalMassLossRate)    { };                                                                         // Do not use core mass calculations during CHE phase

};





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                ABUNDANCE FUNCTIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateHAbundanceCore
 *
 * @brief
 * Calculate the hydrogen abundance in the core of a CH star, given the phase-relative age
 * of the star.  Currently just a simple linear model that assumes that hydrogen in the
 * core of the star is burned to helium at a constant rate throughout the star's lifetime.
 * 
 * Should one day be updated to match detailed models.
 * 
 *
 * double CalculateHAbundanceCore(const double p_Tau, const double p_InitialHAbundance) const
 * 
 * @param       p_Tau                           CH (MS) relative age of the star [0, 1]
 * @param       p_InitialHAbundance             Initial hydrogen abundance of the star
 * @return                                      Hydrogen abundance in the core of the star
 */
GNU_CONST inline double CH::CalculateHAbundanceCore(const double p_Tau, const double p_InitialHAbundance) const {
    return p_InitialHAbundance * (1.0 - p_Tau);
}


/*
 * CalculateHAbundanceSurface
 *
 * @brief
 * Calculate the hydrogen abundance at the surface of a CH star, given the phase-relative age
 * of the star.  Since the star is chemically homogeneous, the core and the surface have the
 * same abundances.
 * 
 *
 * double CalculateHAbundanceSurface(const double p_Tau, const double p_InitialHAbundance) const
 * 
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_InitialHAbundance             Initial hydrogen abundance of the star
 * @return                                      Hydrogen abundance at the surface of the star
 */
GNU_CONST inline double CH::CalculateHAbundanceSurfaceOnPhase(const double p_Tau, const double p_InitialHAbundance) const {
    return CalculateHAbundanceCore(p_Tau, p_InitialHAbundance);
}


/*
 * CalculateHeAbundanceCore
 *
 * @brief
 * Calculate the helium abundance in the core of a CH star, given the phase-relative age
 * of the star.  Currently just a simple linear model from the initial helium abundance
 * to the maximum helium abundance (assuming that all hydrogen is converted to helium). 
 * 
 * Should one day be updated to match detailed models.
 * 
 *
 * double CalculateHeAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_InitialHeAbundance            Initial helium abundance of the star
 * @return                                      Helium abundance in the core of the star
 */
GNU_CONST inline double CH::CalculateHeAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance) const {
    return ((1.0 - p_Metallicity - p_InitialHeAbundance) * p_Tau) + p_InitialHeAbundance;
}


/*
 * CalculateHeAbundanceSurface
 *
 * @brief
 * Calculate the helium abundance at the surface of a CH star, given the phase-relative age
 * of the star.  Since the star is checmically homogeneous, the core and the surface have the
 * same abundances.
 * 
 * 
 * double CalculateHeAbundanceSurface(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance) const
 * 
 * @param       p_Metallicity                   Metallicity of the star (not used here)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_InitialHeAbundance            Initial helium abundance of the star
 * @return                                      Helium abundance at the surface of the star
 */
GNU_CONST inline double CH::CalculateHeAbundanceSurface(const double p_Metallicity, const double p_Tau, const double p_InitialHeAbundance) const {
    return CalculateHeAbundanceCore(p_Tau, p_InitialHeAbundance);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                             LIFETIME / AGE FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateAgeAfterMassLoss_Hurley2000
 *
 * @brief
 * Recalculate the age of a CH star after mass loss, per Hurley et al. 2000, section 7.1.
 * 
 * The age will be further modified if option `--enhance-CHE-lifetimes-luminosities`
 * was specified.
 *
 *
 * double CalculateAgeAfterMassLoss_Hurley2000(const double      p_Mass,
 *                                             const double      p_Age,
 *                                             const double      p_ZetaHurley,
 *                                             const double      p_tMS,
 *                                             const DBL_VECTOR& p_aN,
 *                                             const bool        p_EnhanceLifetime) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Current age of the star (Myr)
 * @param       p_ZetaHurley                    Hurley zeta value (log10(Z / ZSOL_HURLEY))
 * @param       p_tMS                           MS lifetime (Myr)
 * @param       p_aN                            Hurley a(n) coefficients
 * @param       p_EnhanceLifetime               User-specified flag indicating whether to enhance CH lifetime
 *                                              (value of option `--enhance-CHE-lifetimes-luminosities`)
 * @return                                      Age of the star after mass loss (Myr)
 */
inline double CH::CalculateAgeAfterMassLoss_Hurley2000(const double      p_Mass,
                                                       const double      p_Age,
                                                       const double      p_ZetaHurley,
                                                       const double      p_tMS,
                                                       const DBL_VECTOR& p_aN,
                                                       const bool        p_EnhanceLifetime) const {
    const double age = p_Age * CalculateMSLifetime_Hurley2000(p_Mass, p_ZetaHurley, CalculateLifetimeToBGB_Hurley(p_Mass, p_aN), p_aN) / p_tMS;   
    return p_EnhanceLifetime ? age *= CalculateLifetimesRatio(p_Mass) : age;    // enhance lifetime of star if required
}


/*
 * CalculateLifetimesRatio_Szecsi2020
 *
 * @brief
 * Calculate the ratio of the lifetimes of a CH star and a non-CH MS star of the same mass
 * (t_CHE/t_MS).  Uses a polynomial fit derived using IZw18 and IZw18CHE BoOST models from
 * Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203).
 *
 * 
 * double CalculateLifetimesRatio_Szecsi2020(const double p_Mass) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Lifetimes ratio
 */
inline double CH::CalculateLifetimesRatio_Szecsi2020(const double p_Mass) const {

    const double x  = log10(p_Mass);
    const double x2 = x * x;

    return PPOW(10.0, -0.15929168474199387 + (1.050069750483549 * x) + (-0.8233601359988406 * x2) + (0.17772610259473764 * x * x2));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            MASS / MASS LOSS FUNCTIONS                             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateMLRateRotationEnhancement_Langer1998
 *
 * @brief
 * Calculate mass loss rate enhancement for rapidly rotating stars, per Langer 1998, eq 3
 * (See https://ui.adsabs.harvard.edu/abs/1998A%26A...329..551L/abstract)
 * 
 * The exponent originally comes from Bjorkman & Cassinelli 1993
 * (see https://ui.adsabs.harvard.edu/abs/1993ApJ...409..429B/abstract), based on a fit to data
 * from Friend & Abbott 1986 (see https://ui.adsabs.harvard.edu/abs/1986ApJ...311..701F/abstract) 
 *
 * 
 * double CalculateMLRateRotationEnhancement_Langer1998() const
 *
 * @return                                      Mass loss enhancement factor for rapidly rotating stars
 */
GNU_CONST inline double CalculateMLRateRotationEnhancement_Langer1998() const {
    return PPOW((1.0 - Omega() / OmegaBreak()), -0.43);
}


/*
 * CalculateMLFractionOB
 *
 * @brief
 * Calculate the fraction of mass loss attributable to OB mass loss, for a CH star,
 * per Yoon et al. 2006
 * 
 * The model described in Yoon et al. 2006 (also Szecsi et al. 2015) uses OB mass loss while the
 * He surface abundance is below 0.55, WR mass loss when the surface He abundance is above 0.7,
 * and linearly interpolates when the He surface abundance is between those limits.
 * 
 * This function calculates the fraction of mass loss attributable to OB mass loss, based on
 * the He surface abundance and the abundance limits described in Yoon et al. 2006.  The value
 * returned will be 1.0 if 100% of the mass loss is attributable to OB mass lass, 0.0 if 100% of
 * the mass loss is attributable to WR mass loss, and in the range (0.0, 1.0) if the mass loss is
 * a mix of OB and WR.
 * 
 *
 * double CalculateMLFractionOB(const double p_HeAbundanceSurface) const
 *
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Fraction of mass loss attributable to OB mass loss
 */
GNU_CONST inline double CH::CalculateMLFractionOB(const double p_HeAbundanceSurface) const {

    constexpr double limOB = 0.55;      // per Yoon et al. 2006
    constexpr double limWR = 0.70;      // per Yoon et al. 2006

    return std::min(1.0, std::max (0.0, (limWR - p_HeAbundanceSurface) / (limWR - limOB)));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LUMINOSITY CALCULATIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
 * CalculateLuminosityAtPhaseEnd
 *
 * @brief
 * Calculate the luminosity of a CH star at the end of the (CH) MS.  The luminosity
 * will be enhanced if option `--enhance-CHE-lifetimes-luminosities` was specified.
 * 
 * 
 * double CalculateLuminosityAtPhaseEnd(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAMS CH luminosity (Lsol)
 */
inline double CH::CalculateLuminosityAtPhaseEnd(const double p_Mass) const {

    // unenhanced CH luminosity at the end of the MS is just MS luminosity at the end of the MS
    const double luminosity = MainSequence::CalculateLuminosityAtPhaseEnd(p_Mass);
    
    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {                       // enhance luminosity of CH stars?
                                                                            // yes
        // apply enhancement, which at TAMS is just the ratio log(L_CH) / log(L_MS)
        // enhancement should not reduce luminosity, so ratio is clamped to a minimum of +1.0
        luminosity = PPOW(10.0, log10(luminosity) * std::max(CalculateLogLuminositiesRatio(p_Mass), 1.0));
    }

    return luminosity;
}



#endif // __CH_h__
