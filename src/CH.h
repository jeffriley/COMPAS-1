#ifndef __CH_h__
#define __CH_h__

#include "constants.h"
//// #include "typedefs.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "profiling.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "utils.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#include "BaseStar.h"
#include "MS_gt_07.h"

class BaseStar;
class MS_gt_07;

class CH: virtual public BaseStar, public MS_gt_07 {

public:

    // constructors

    CH() { m_StellarType = STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS; };
    
    CH(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MS_gt_07(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS; // Set stellar type // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if (p_Initialise) Initialise();  // Initialise if required
    }


private: 

    void Initialise() {
        CalculateTimescales(); // Initialise timescales
        m_Age = 0.0; // Set age appropriately   // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        m_CHE = true; // initially for CH stars  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                                                                                                                                                          // Set age appropriately
    }


protected:

    // member functions


///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // abundances
    GNU_CONST   double CalculateHAbundanceCore(const double p_Tau, const double p_InitialHAbundance) const override;
    GNU_CONST   double CalculateHAbundanceSurface(const double p_Tau, const double p_InitialHAbundance) const override;
    COMPAS_PURE double CalculateHeAbundanceCore(const double p_Tau, const double p_InitialHeAbundance = 0.0) const override;
    COMPAS_PURE double CalculateHeAbundanceSurface(const double p_Tau, const double p_InitialHeAbundance) const override;


    // age, lifetime, tau, timescales

    COMPAS_PURE double CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tMS, const bool p_EnhanceLifetime) const override;

    GNU_CONST double CalculateLifetimesRatio_Szecsi2020(const double p_Mass) const;

    inline DBL_VECTOR CalculateTimescales_Hurley2000() const override {
        return CalculateTimescales_Hurley2000(m_StateHistory.CurrentState.MassEffectiveInitial(), m_StateHistory.CurrentState.GBparams(), m_StateHistory.CurrentState.TimeScales());
    }

    COMPAS_PURE DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;


    // luminosity

    GNU_CONST double CalculateLogLuminositiesRatio(const double p_Mass) const;

    COMPAS_PURE double CH::CalculateLuminosity_Hurley2000(
        const double p_Mass,
        const double p_Tau,
        const double p_Time,
        const double p_LZAMS,
        const double p_tMS,
        const double p_tBGB
    ) const;


    // radius
    inline double CalculateRadius() const override { return m_StateHistory.ZAMSState().Radius(); }              // CH radius is fixed throughout lifetime
    inline double CalculateRadiusAtPhaseEnd() const override { return m_StateHistory.ZAMSState().Radius(); }    // CH radius is fixed throughout lifetime 















///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<








    // Lifetime
    double          CalculateLogLifetimeRatio(const double p_Mass) const;
    double          CalculateLifetimeRatio(const double p_Mass) const;

    // Luminosity
    double          CalculateLogLuminosityRatio(const double p_Mass, const double p_Tau) const;

    double          CalculateLuminosityAtPhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtPhaseEnd() const               { return CalculateLuminosityAtPhaseEnd(m_Mass0); }                                                          // Use class member variables


double CalculateLuminosityOnPhase() const; // Uses globals and state variables
double CalculateLuminosityOnPhase(
    m_StateHistory.CurrentState().Mass(),
    m_StateHistory.ZAMSState().Luminosity(),
    m_StateHistory.Timescales(),
    m_StateHistory.CurrentState().Time(),
    GLOBALS->ReferenceMetallicity(),
    GLOBALS->HurleyACoefficients(),
    GLOBALS->HurleyLuminosityConstants()
);

    
    // Mass loss rate
    COMPAS_PURE MASS_LOSS_T CalculateMLrate_Belczynski2010(
        const double p_Mass,
        const double p_Radius,
        const double p_Luminosity,
        const double p_Temperature,
        const double p_PerturbationMu,
        const double p_HeAbundanceSurface
    ) const override;

    COMPAS_PURE MASS_LOSS_T CH::CalculateMLrate_Merritt2025(
        const double p_Mass,
        const double p_Radius,
        const double p_Luminosity,
        const double p_Temperature,
        const double p_PerturbationMu,
        const double p_mStart,
        const double p_HeAbundanceSurface,
    ) const override;


    double          CalculateMassLossRateBelczynski2010()               { return BaseStar::CalculateMassLossRateBelczynski2010() * CalculateMassLossRateEnhancementRotation(); }
    double          CalculateMassLossRateMerritt2025()                  { return BaseStar::CalculateMassLossRateBelczynski2010() * CalculateMassLossRateEnhancementRotation(); }


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
//                                     ABUNDANCE                                     //
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
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
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
GNU_CONST inline double CH::CalculateHAbundanceSurface(const double p_Tau, const double p_InitialHAbundance) const {
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
 * double CalculateHeAbundanceCore(const double p_Tau, const double p_InitialHeAbundance) const
 * 
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_InitialHeAbundance            Initial helium abundance of the star
 * @return                                      Helium abundance in the core of the star
 */
COMPAS_PURE inline double CH::CalculateHeAbundanceCore(const double p_Tau, const double p_InitialHeAbundance) const {
    return ((1.0 - GLOBALS->ReferenceMetallicity() - p_InitialHeAbundance) * p_Tau) + p_InitialHeAbundance;
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
 * double CalculateHeAbundanceSurface(const double p_Tau, const double p_InitialHeAbundance) const
 * 
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_InitialHeAbundance            Initial helium abundance of the star
 * @return                                      Helium abundance at the surface of the star
 */
COMPAS_PURE inline double CH::CalculateHeAbundanceSurface(const double p_Tau, const double p_InitialHeAbundance) const {
    return CalculateHeAbundanceCore(p_Tau, p_InitialHeAbundance);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         AGE / LIFETIME / TAU / TIMESCALES                         //
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
 * double CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tMS, const bool p_EnhanceLifetime) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tMS                           MS lifetime, tMS (per Hurley timescales) (Myr)
 * @param       p_EnhanceLifetime               User-specified flag indicating whether to enhance CH lifetime
 *                                              (value of option `--enhance-CHE-lifetimes-luminosities`)
 * @return                                      Age of the star after mass loss (Myr)
 */
COMPAS_PURE inline double CH::CalculateAgeAfterMassLoss_Hurley2000(const double p_Mass, const double p_Age, const double p_tMS, const bool p_EnhanceLifetime) const {

    // JR FIX THIS AFTER TALKING TO ILYA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                                        // why not timescales tBGB ???????????????????? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    const double age = p_Age * MainSequence::CalculatePhaseLifetime_Hurley2000(p_Mass, BaseStar::CalculateLifetimeToBGB_Hurley2000_Static(p_Mass)) / p_tMS;   

    return p_EnhanceLifetime ? age *= CalculateLifetimesRatio_Szecsi2020(p_Mass) : age; // enhance lifetime of star if required
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
GNU_CONST inline double CH::CalculateLifetimesRatio_Szecsi2020(const double p_Mass) const {

    const double x  = log10(p_Mass);
    const double x2 = x * x;

    return PPOW(10.0, -0.15929168474199387 + (1.050069750483549 * x) + (-0.8233601359988406 * x2) + (0.17772610259473764 * x * x2));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     MASS LOSS                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateMLrateRotationEnhancement_Langer1998
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
 * double CalculateMLrateRotationEnhancement_Langer1998() const
 *
 * @return                                      Mass loss enhancement factor for rapidly rotating stars
 */
GNU_CONST inline double CalculateMLrateRotationEnhancement_Langer1998() const {
    return PPOW((1.0 - Omega() / OmegaBreak()), -0.43);
}


/*
 * CalculateMLrate_Belczynski2010
 *
 * @brief
 * Calculate the mass loss rate, per Belczynski et al. 2010, based on the StarTrack
 * implementation, and modified for CH stars.
 * 
 * If option `--enable-rotationally-enhanced-mass-loss` was specified, the  mass
 * loss rate will be enhanced for rotation.
 * 
 * 
 * MASS_LOSS_T CalculateMLrate_Belczynski2010(
 *     const double p_Mass,
 *     const double p_Radius,
 *     const double p_Luminosity,
 *     const double p_Temperature,
 *     const double p_PerturbationMu,
 *     const double p_HeAbundanceSurface,
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type
 *                                                                  (will be MASS_LOSS_TYPE::WR or MASS_LOSS_TYPE::OB)
 */
COMPAS_PURE inline MASS_LOSS_T CH::CalculateMLrate_Belczynski2010(
    const double p_Mass,
    const double p_Radius,
    const double p_Luminosity,
    const double p_Temperature,
    const double p_PerturbationMu,
    const double p_HeAbundanceSurface
) const {

    MASS_LOSS_TYPE dominantMLtype;
    double dMdt;
    std::tie(dMdt, dominantMLtype) = BaseStar::CalculateMLrate_Belczynski2010(
        p_Mass, 
        p_Radius, 
        p_Luminosity, 
        p_Temperature, 
        p_PerturbationMu, 
        p_HeAbundanceSurface
    );

    // return mass loss rate, enhanced for rotation if required, and dominant mass loss type
    return std::make_tuple(dMdt * (OPTIONS->EnableRotationallyEnhancedMassLoss() ? CalculateMLrateRotationEnhancement_Langer1998() : 1.0), dominantMLtype);
}


/*
 * CalculateMLrate_Merritt2025
 *
 * Calculate mass loss rate, and dominant mass loss type, at the current evolutionary phase,
 * per Merritt et al., 2025.
 * Calculate the mass loss rate, per Merritt et al. 2024, modified for CH stars.
 * 
 * If option `--enable-rotationally-enhanced-mass-loss` was specified, the mass loss rate will
 * be enhanced for rotation.
 * 
 * 
 * MASS_LOSS_T CalculateMLrate_Merritt2025(
 *     const double p_Mass,
 *     const double p_Radius,
 *     const double p_Luminosity,
 *     const double p_Temperature,
 *     const double p_PerturbationMu,
 *     const double p_mStart,
 *     const double p_HeAbundanceSurface
 * ) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T inline CH::CalculateMLrate_Merritt2025(
    const double p_Mass,
    const double p_Radius,
    const double p_Luminosity,
    const double p_Temperature,
    const double p_PerturbationMu,
    const double p_mStart,
    const double p_HeAbundanceSurface
) const {

    MASS_LOSS_TYPE dominantMLtype;
    double dMdt;
    std::tie(dMdt, dominantMLtype) = BaseStar::CalculateMLrate_Merritt2025(
        p_Mass, 
        p_Radius, 
        p_Luminosity, 
        p_Temperature, 
        p_PerturbationMu, 
        p_mStart, 
        p_HeAbundanceSurface
    );

    // return mass loss rate, enhanced for rotation if required, and dominant mass loss type
    return std::make_tuple(dMdt * (OPTIONS->EnableRotationallyEnhancedMassLoss() ? CalculateMLrateRotationEnhancement_Langer1998() : 1.0), dominantMLtype);
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
GNU_CONST inline double CH::CalculateLogLuminositiesRatio(const double p_Mass) const {
         
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
