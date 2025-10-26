#include "CH.h"







/*
 * CalculateLuminosityOnPhase
 *
 * @brief
 * Calculate the luminosity of a CH star on the (CH) MS.  The luminosity will be enhanced
 * if option `--enhance-CHE-lifetimes-luminosities` was specified.
 * 
 * 
 * double CalculateLuminosityOnPhase(const double      p_Metallicity,
 *                                   const double      p_Mass,
 *                                   const double      p_Time,
 *                                   const double      p_LZAMS,
 *                                   const DBL_VECTOR& p_tScales,
 *                                   const DBL_VECTOR& p_aN,
 *                                   const DBL_VECTOR& p_LConstants) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Time                          Time elapsed since ZAMS (Myr)
 * @param       p_LZAMS                         ZAMS luminosity of the star (Lsol)
 * @param       p_tScales                       Phase timescales
 * @param       p_aN                            Hurley a(n) coefficients
 * @param       p_LConstants                    Hurley luminosity constants
 * @return                                      CH luminosity (Lsol)
 */
double CH::CalculateLuminosityOnPhase(const double      p_Metallicity,
                                      const double      p_Mass,
                                      const double      p_Time,
                                      const double      p_LZAMS,
                                      const DBL_VECTOR& p_tScales,
                                      const DBL_VECTOR& p_aN,
                                      const DBL_VECTOR& p_LConstants) const {

    // unenhanced CH luminosity is just MS luminosity
    const double LZAMS = m_StateHistory.ZAMSState().Luminosity();
    const double mass  = m_StateHistory.CurrentState().Mass();
    const double time  = m_StateHistory.CurrentState().Time();

    double luminosity = MainSequence::CalculateLuminosity(p_Metallicity, p_Mass, p_Time, p_LZAMS, p_tScales, p_aN, p_LConstants);

    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {                           // enhance luminosity of CH stars?
                                                                                // yes
        const double tau = m_StateHistory.CurrentState().Tau();

        // enhancement should not reduce luminosity, so ratio is clamped to a minimum of +1.0
        // enhancement amount grows from 1 to logLuminosityRatio over main-sequence
        const double enhancement = 1.0 + (std::max(CalculateLogLuminositiesRatio(mass), 1.0) - 1.0) * tau * tau;

        luminosity = PPOW(10.0, log10(luminosity) * enhancement);               // apply enhancement
    }

    return luminosity;
}


/*
 * CalculateTimescales_Hurley2000
 *
 * @brief
 * (Re)calculate timescales given the mass of the star, per Hurley at al. 2000.
 * 
 * Since timescales depend on a star's mass, they need to be calculated whenever
 * the mass of the star changes (probably every timestep).
 *
 *
 * DBL_VECTOR CalculateTimescales_Hurley2000(const double      p_Mass,
 *                                           const double      p_ZetaHurley,
 *                                           const DBL_VECTOR& p_GBparams,
 *                                           const DBL_VECTOR& p_MassCutoffs,
 *                                           const DBL_VECTOR& p_Timescales,
 *                                           const double      p_Alpha3,
 *                                           const DBL_VECTOR& p_aN,
 *                                           const DBL_VECTOR& p_bN) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_ZetaHurley                    Hurley zeta value (log10(Z / ZSOL_HURLEY))
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_MassCutoffs                   Hurley mass cutoffs (Msol)
 * @param       p_tScales                       Hurley timescales (Myr)
 * @param       p_Alpha3                        Hurley alpha3 constant
 * @param       p_aN                            Hurley a(n) coefficients
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      Mutated timescales (Myr)
 */
DBL_VECTOR CH::CalculateTimescales_Hurley2000(const double      p_Mass,
                                              const double      p_ZetaHurley,
                                              const DBL_VECTOR& p_GBparams,
                                              const DBL_VECTOR& p_MassCutoffs,
                                              const DBL_VECTOR& p_tScales,
                                              const double      p_Alpha3,
                                              const DBL_VECTOR& p_aN,
                                              const DBL_VECTOR& p_bN) const {

    DBL_VECTOR tScales = p_tScales;                     // copy given timescales

    // (re)calculate MS timescales
    tScales = MainSequence::CalculateTimescales_Hurley2000(p_Mass, p_ZetaHurley, p_GBparams, p_MassCutoffs, tScales, p_Alpha3, p_aN, p_bN);

    // enhamce lifetimes as appropriate
    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {   // enhance lifetime of CH stars?
                                                        // yes
        const double lifetimesRatio = CalculateLifetimesRatio_Szecsi2020(p_Mass);

        tScales[static_cast<int>(TIMESCALE::tBGB)] *= lifetimesRatio;
        tScales[static_cast<int>(TIMESCALE::tMS)]  *= lifetimesRatio;
    }

    // return timescales vector by value - NRVO takes care of performance/efficiency
    return tScales;
}


/*
 * CalculateMLRate_Belczynski2010
 *
 * @brief
 * Calculate the mass loss rate, per Belczynski et al. 2010, based on the StarTrack
 * implementation, and modified for CH stars.
 *
 * If option `--scale-CHE-mass-loss-with-surface-helium-abundance` was specified,
 * the mass loss rate will be scaled with the surface helium abundance.
 * 
 * If option `--enable-rotationally-enhanced-mass-loss` was specified, the  mass
 * loss rate will be enhanced for rotation.
 * 
 * 
 * MASS_LOSS_T CalculateMLRate_Belczynski2010(const double                     p_Metallicity,
 *                                            const double                     p_Mass,
 *                                            const double                     p_Radius,
 *                                            const double                     p_Luminosity,
 *                                            const double                     p_Temperature,
 *                                            const double                     p_PerturbationMu,
 *                                            const double                     p_ZscaledHurley,
 *                                            const double                     p_HeAbundanceSurface,
 *                                            const double                     p_CoolWindsMultiplier,
 *                                            const double                     p_WRfactor,
 *                                            const bool                       p_ScaleWithSurfaceHe,
 *                                            const LBV_MASS_LOSS_PRESCRIPTION p_LBVprescription) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @param       p_ZscaledHurley                 Z inversely scaled by Hurley ZSOL (Z / ZSOL_HURLEY)
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @param       p_CoolWindsMultiplier           Cool winds mass loss multiplier
 * @param       p_WRfactor                      WR mass loss factor
 * @param       p_ScaleWithSurfaceHelium        Indicates whether mass loss should be scaled with surface He abundance
 * @param       p_LBVprescription               LBV mass loss prescription to use
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type
 *                                                                  (will be MASS_LOSS_TYPE::WR or MASS_LOSS_TYPE::OB)
 */
COMPAS_PURE MASS_LOSS_T CH::CalculateMLRate_Belczynski2010(const double                     p_Metallicity,
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
                                                           const LBV_MASS_LOSS_PRESCRIPTION p_LBVprescription) const {
    // set defaults
    MASS_LOSS_TYPE dominantMLType = MASS_LOSS_TYPE::OB;
    double dMdt = BaseStar::CalculateMLRate_Belczynski2010(p_Metallicity,
                                                           p_Mass,
                                                           p_Radius,
                                                           p_Luminosity,
                                                           p_Temperature,
                                                           p_PerturbationMu,
                                                           p_ZscaledHurley,
                                                           p_HeAbundanceSurface,
                                                           p_CoolWindsMultiplier,
                                                           p_WRfactor,
                                                           p_ScaleWithSurfaceHe,
                                                           p_LBVprescription);

    // scale mass loss with the surface helium abundance if necessary
    if (p_ScaleWithSurfaceHe) {                                                                         // transition between OB and WR mass loss rates?
                                                                                                        // yes
        const double dMdtWR     = BaseStar::CalculateMLRateWR_ZDependent_Static(p_Luminosity, p_Metallicity, 0.0); // WR mass loss rate
        const double fractionOB = CalculateMLFractionOB(p_HeAbundanceSurface);                          // mass loss fraction attributable to OB mass loss
        
        if (((1.0 - fractionOB) * dMdtWR) > (fractionOB * dMdt)) dominantMLType = MASS_LOSS_TYPE::WR;   // dominant mass loss type

        dMdt = (fractionOB * dMdt) + ((1.0 - fractionOB) * dMdtWR);                                     // combined mass loss rate
    }

    // calculate mass loss enhancement due to rotation, and dominant mass loss type
    const double rotEnhancement = OPTIONS->EnableRotationallyEnhancedMassLoss() ? CalculateMLRateRotationEnhancement_Langer1998() : 1.0;

    return std::make_tuple(dMdt * rotEnhancement, dominantMLType);
}


/*
 * CalculateMLRate_Merritt2025
 *
 * @brief
 * Calculate the mass loss rate, per Merritt et al. 2024, modified for CH stars.
 *
 * If option `--scale-CHE-mass-loss-with-surface-helium-abundance` was specified,
 * the mass loss rate will be scaled with the surface helium abundance.
 * 
 * If option `--enable-rotationally-enhanced-mass-loss` was specified, the  mass
 * loss rate will be enhanced for rotation.
 * 
 * 
 * MASS_LOSS_T CalculateMLRate_Merritt2025(const double p_Metallicity, const double p_Luminosity, const double p_HeAbundanceSurface) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_SigmaHurley                   Hurley sigma value (log10(Z))
 * @param       p_ZetaAnders                    Anders zeta value (log10(Z / ZSOL_ANDERS))
 * @param       p_ZetaAsplund                   Asplund zeta value (log10(Z / ZSOL_ASPLUND))
 * @param       p_ZscaledHurley                 Z inversely scaled by Hurley ZSOL (Z / ZSOL_HURLEY)
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @param       p_WRfactor                      WR mass loss factor
 * @param       p_TerminalWindScalePower        Power with which to scale terminal wind velocity with metallicity
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type
 *                                                                  (will be MASS_LOSS_TYPE::WR or MASS_LOSS_TYPE::OB)
 */
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
                                                        const double p_TerminalWindScalePower) const {

    // set defaults
    MASS_LOSS_TYPE dominantMLType = MASS_LOSS_TYPE::OB;
    double dMdt = BaseStar::CalculateMLRateOB(p_Metallicity, p_Mass, p_Luminosity, p_Temperature, p_ZetaAnders, p_ZetaAsplund, p_TerminalWindScalePower, OPTIONS->OBMassLossPrescription());

    // scale mass loss with the surface helium abundance if necessary
    if (OPTIONS->ScaleCHEMassLossWithSurfaceHeliumAbundance()) {                                        // transition between OB and WR mass loss rates?
                                                                                                        // yes
        const double dMdtWR     = HeMS::CalculateMLRateMerritt2025_Static(p_Metallicity, p_Luminosity, p_Temperature, p_SigmaHurley, p_ZetaAnders); // WR mass loss rate
        const double fractionOB = CalculateMLFractionOB(p_HeAbundanceSurface);                          // mass loss fraction attributable to OB mass loss

        if ((1.0 - fractionOB) * dMdtWR > fractionOB * dMdt) dominantMLType = MASS_LOSS_TYPE::WR;       // dominant mass loss type

        dMdt = (fractionOB * dMdt) + ((1.0 - fractionOB) * dMdtWR);                                     // combined mass loss rate
    }

    // calculate mass loss enhancement due to rotation, and dominant mass loss type
    const double rotEnhancement = OPTIONS->EnableRotationallyEnhancedMassLoss() ? CalculateMLRateRotationEnhancement_Langer1998() : 1.0;

    return std::make_tuple(dMdt * rotEnhancement, dominantMLType);
}





STELLAR_TYPE CH::EvolveToNextPhase() {

    STELLAR_TYPE stellarType = STELLAR_TYPE::MS_GT_07;

    if (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]) {           // evolving off because of age?
        stellarType = STELLAR_TYPE::MS_GT_07;                               // no - must have spun down - evolve as MS star now
        m_CHE       = false;                                                // evolved CH->MS
        
        // if BRCEK core mass calculations enabled, initialise the core mass based on current mass and central helium fraction
        if ((OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK) && (utils::Compare(m_MZAMS, BRCEK_LOWER_MASS_LIMIT) >= 0))
            m_MainSequenceCoreMass = MainSequence::CalculateInitialMainSequenceCoreMass(m_Mass, m_HeliumAbundanceCore);
    }
    else {                                                                  // yes
        stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                   // evolve as HeMS star now
        m_Age       = 0.0;                                                  // can't use Hurley et al. 2000, eq 76 here - timescales(tHe) not calculated yet
        m_Tau       = 0.0;
        m_CHE       = true;                                                 // stayed on MS as CH
    }

    return stellarType;
}

