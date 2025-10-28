#include "CH.h"







/*
 * CalculateLuminosityOnPhase
 *
 * @brief
 * Calculate the luminosity of a CH star on the (CH) MS.  The luminosity will be
 * enhanced if option `--enhance-CHE-lifetimes-luminosities` was specified.
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
 * (Re)calculate CH timescales.
 * CH timescales are MS timescales, per Hurley at al. 2000, optionally enhanced
 * if the user specified the `--enhance-CHE-lifetimes-luminosities` option.
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

