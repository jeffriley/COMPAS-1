#include "CH.h"




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         AGE / LIFETIME / TAU / TIMESCALES                         //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
 * DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_Timescales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
COMPAS_PURE DBL_VECTOR CH::CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const {

    DBL_VECTOR tScales = p_tScales;                     // copy given timescales

    // (re)calculate MS timescales
    tScales = MainSequence::CalculateTimescales_Hurley2000(p_Mass, p_GBparams, tScales);

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


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLuminosity_Hurley2000
 *
 * @brief
 * Calculate the luminosity of a CH star on the (CH) MS, per Hurley et al. 2000, eq 12
 * 
 * The luminosity will be enhanced if option `--enhance-CHE-lifetimes-luminosities` was specified.
 * 
 * 
 * double CalculateLuminosity_Hurley2000(
 *     const double p_Mass,
 *     const double p_Tau,
 *     const double p_Age,
 *     const double p_LZAMS,
 *     const double p_tMS,
 *     const double p_tBGB
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_Age                           Time elapsed since ZAMS (Myr)
 * @param       p_LZAMS                         ZAMS luminosity of the star (Lsol)
 * @param       p_tMS                           MS lifetime, tMS (per Hurley timescales) (Myr)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (per Hurley timescales) (Myr)
 * @return                                      CH luminosity (Lsol)
 */
COMPAS_PURE double CH::CalculateLuminosity_Hurley2000(
    const double p_Mass,
    const double p_Tau,
    const double p_Age,
    const double p_LZAMS,
    const double p_tMS,
    const double p_tBGB
) const {

    // unenhanced CH luminosity is just MS luminosity
    double luminosity = MainSequence::CalculateLuminosity_Hurley2000(p_Mass, p_Time, p_LZAMS, p_tMS, p_tBGB);

    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {               // enhance luminosity of CH stars?
                                                                    // yes
        // enhancement should not reduce luminosity, so ratio is clamped to a minimum of +1.0
        // enhancement amount grows from 1 to logLuminosityRatio over main-sequence
        const double enhancement = 1.0 + (std::max(CalculateLogLuminositiesRatio(mass), 1.0) - 1.0) * p_Tau * p_Tau;

        luminosity = PPOW(10.0, log10(luminosity) * enhancement);   // apply enhancement
    }

    return luminosity;
}










STELLAR_TYPE CH::EvolveToNextPhase() {

    STELLAR_TYPE stellarType = STELLAR_TYPE::MS_GT_07;

    if (m_Age < m_Timescales[static_cast<int>(TIMESCALE::tMS)]) {           // evolving off because of age?
        stellarType = STELLAR_TYPE::MS_GT_07;                               // no - must have spun down - evolve as MS star now
        m_CHE       = false;                                                // evolved CH->MS
        
        // if BRCEK core mass calculations enabled, initialise the core mass based on current mass and central helium fraction
        if ((OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK) && (utils::Compare(m_MZAMS, BRCEK_LOWER_MASS_LIMIT) >= 0))
            m_MainSequenceCoreMass = MainSequence::CalculateCNOprocessedCoreMass_Brcek2025(m_Mass, m_HeliumAbundanceCore);
    }
    else {                                                                  // yes
        stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                   // evolve as HeMS star now
        m_Age       = 0.0;                                                  // can't use Hurley et al. 2000, eq 76 here - timescales(tHe) not calculated yet
        m_Tau       = 0.0;
        m_CHE       = true;                                                 // stayed on MS as CH
    }

    return stellarType;
}

