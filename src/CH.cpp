#include "CH.h"


/*
 * CalculateAgeAfterMassLoss_Hurley
 *
 * @brief
 * Recalculate the star's age after mass loss, per Hurley et al. 2000, section 7.1
 * 
 * Note that the CH timescales should be recalculated after recalculating the star's age.
 *
 *
 * CalculateAgeAfterMassLoss_Hurley(const double p_Mass, const double p_Age) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Current age of the star (Myr)
 * @param       p_tMS                           MS timescale (Myr)
 * @return                                      Age of the star after mass loss (Myr)
 */
double CH::CalculateAgeAfterMassLoss_Hurley(const double p_Mass, const double p_Age, const double p_tMS) const {

    // calculate tMS for mass as passed (tMS')
    double tMSprime = MainSequence::CalculateLifetimeOnPhase(p_Mass, CalculateLifetimeToBGB(p_Mass));

    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {   // enhance lifetime of CH stars?                      
        tMSprime *= CalculateLifetimesRatio(p_Mass);    // yes
    }

    return p_Age * tMSprime / p_tMS;
}


/*
 * CalculateLuminosityAtPhaseEnd
 *
 * @brief
 * Calculate the luminosity of a CH star at the end of the (CH) MS.  The luminosity will be
 * enhanced if option `--enhance-CHE-lifetimes-luminosities` was specified.
 * 
 * 
 * double CalculateLuminosityAtPhaseEnd(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAMS CH luminosity (Lsol)
 */
double CH::CalculateLuminosityAtPhaseEnd(const double p_Mass) const {

    // unenhanced CH luminosity at the end of the MS is just MS luminosity at the end of the MS
    const double luminosity = MainSequence::CalculateLuminosityAtPhaseEnd(p_Mass);
    
    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {                       // enhance luminosity of CH stars?
                                                                            // yes
        const double ratio = CalculateLogLuminositiesRatio(p_Mass);         // log(L_CH) / log(L_MS)

        // apply enhancement, which at TAMS is just the ratio log(L_CH) / log(L_MS)
        // enhancement should not reduce luminosity, so ratio is clamped to a minimum of +1.0
        luminosity = PPOW(10.0, log10(luminosity) * std::max(ratio, 1.0));
    }

    return luminosity;
}


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
 *                                   const DBL_VECTOR& p_Timescales,
 *                                   const DBL_VECTOR& p_aCoefficients,
 *                                   const DBL_VECTOR& p_LConstants) const
 *
 * @param       p_Metallicity                   (Fractional) metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Time                          Time elapsed since ZAMS (Myr)
 * @param       p_LZAMS                         ZAMS luminosity of the star (Lsol)
 * @param       p_Timescales                    Hurley timescales
 * @param       p_aCoefficients                 Hurley a(n) coefficients
 * @param       p_LConstants                    Hurley luminosity constants
 * @return                                      CH luminosity (Lsol)
 */
double CH::CalculateLuminosityOnPhase(const double      p_Metallicity,
                                      const double      p_Mass,
                                      const double      p_Time,
                                      const double      p_LZAMS,
                                      const DBL_VECTOR& p_Timescales,
                                      const DBL_VECTOR& p_aCoefficients,
                                      const DBL_VECTOR& p_LConstants) const {

    // unenhanced CH luminosity is just MS luminosity
    const double LZAMS = m_StateHistory.ZAMSState().Luminosity();
    const double mass  = m_StateHistory.CurrentState().Mass();
    const double time  = m_StateHistory.CurrentState().Time();

    double luminosity = MainSequence::CalculateLuminosity(p_Metallicity, p_Mass, p_Time, p_LZAMS, p_Timescales, p_aCoefficients, p_LConstants);

    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {                           // enhance luminosity of CH stars?
                                                                                // yes
        const double tau   = m_StateHistory.CurrentState().Tau();
        const double ratio = CalculateLogLuminositiesRatio(mass);               // log(L_CH) / log(L_MS)

        // enhancement should not reduce luminosity, so ratio is clamped to a minimum of +1.0
        // enhancement amount grows from 1 to logLuminosityRatio over main-sequence
        const double enhancement = 1.0 + (std::max(ratio, 1.0) - 1.0) * tau * tau;

        luminosity = PPOW(10.0, log10(luminosity) * enhancement);               // apply enhancement
    }

    return luminosity;
}


/*
 * CalculateTimescales
 *
 * @brief
 * (Re)calculate the CH timescales given the mass of the star.  Only timescales
 * relevant to CH are modified.  Since timescales depend on a star's mass, they
 * need to be calculated whenever the mass of the star changes (at least at each
 * timestep).
 *
 *
 * DBL_VECTOR CalculateTimescales(const double p_Mass, const DBL_VECTOR& p_Timescales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Timescales                    Timescales vector
 * @return                                      Mutated timescales vector
 */
DBL_VECTOR CH::CalculateTimescales(const double p_Mass, const DBL_VECTOR& p_Timescales) const {

    DBL_VECTOR timescales = p_Timescales;                               // copy given timescales

    // (re)calculate tBGB and tMS
    timescales[static_cast<int>(TIMESCALE::tBGB)] = CalculateLifetimeToBGB(p_Mass);
    timescales[static_cast<int>(TIMESCALE::tMS)]  = CalculateLifetimeOnPhase(p_Mass, timescales[static_cast<int>(TIMESCALE::tBGB)]);

    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {                   // enhance lifetime of CH stars?
                                                                        // yes
        const double lifetimesRatio = CalculateLifetimesRatio(p_Mass);

        timescales[static_cast<int>(TIMESCALE::tBGB)] *= lifetimesRatio;
        timescales[static_cast<int>(TIMESCALE::tMS)]  *= lifetimesRatio;
    }

    // return timescales vector by value - NRVO takes care of performance/efficiency
    return timescales;
}


/*
 * CalculateMLRateBelczynski2010
 *
 * @brief
 * Calculate the mass loss rate, per Belczynski et al. 2010, based on the StarTrack
 * implementation, and modified for CH stars.
 * 
 * 
 * std::tuple<double, MASS_LOSS_TYPE> CalculateMLRateBelczynski2010(const double p_Metallicity, const double p_Luminosity, const double p_HeAbundanceSurface) const
 *
 * @param       p_Metallicity                   (Fractional) metallicity of the star
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type
 *                                                                  (will be MASS_LOSS_TYPE::WR or MASS_LOSS_TYPE::OB)
 */
std::tuple<double, MASS_LOSS_TYPE> CH::CalculateMLRateBelczynski2010(const double p_Metallicity, const double p_Luminosity, const double p_HeAbundanceSurface) const {

    MASS_LOSS_TYPE dominantMLType = MASS_LOSS_TYPE::OB;                                                 // set default dominant mass loss type
    double dMdt                   = BaseStar::CalculateMLRateBelczynski2010();                          // set default mass loss rate (OB)

    if (OPTIONS->ScaleCHEMassLossWithSurfaceHeliumAbundance()) {                                        // transition between OB and WR mass loss rates?
                                                                                                        // yes
        const double dMdtWR     = BaseStar::CalculateMLRateWRZDependent_Static(p_Luminosity, p_Metallicity, 0.0); // WR mass loss rate
        const double OBfraction = CalculateMassLossFractionOB(p_HeAbundanceSurface);                    // mass loss fraction attributable to OB mass loss
        
        dMdt = (OBfraction * dMdt) + ((1.0 - OBfraction) * dMdtWR);                                     // combined mass loss rate

        if (((1.0 - fractionOB) * MdotWR) > (fractionOB * MdotOB)) dominantMLType = MASS_LOSS_TYPE::WR; // determine dominant mass loss rate
    }

    // return mass loss rate, enhanced due to rotation, and dominant mass loss type
    return std::make_tuple(dMdt * CalculateMLRateEnhancementRotation(), dominantMLType);
}


/*
 * CalculateMLRateMerritt2025
 *
 * @brief
 * Calculate the mass loss rate, per Merritt et al. 2024, modified for CH stars.
 *
 * 
 * std::tuple<double, MASS_LOSS_TYPE> CalculateMLRateMerritt2025(const double p_Metallicity, const double p_Luminosity, const double p_HeAbundanceSurface) const
 * 
 * @param       p_Metallicity                   (Fractional) metallicity of the star
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type
 *                                                                  (will be MASS_LOSS_TYPE::WR or MASS_LOSS_TYPE::OB)
 */
std::tuple<double, MASS_LOSS_TYPE> CH::CalculateMLRateMerritt2025(const double p_Metallicity, const double p_Luminosity, const double p_HeAbundanceSurface) const {

    MASS_LOSS_TYPE dominantMLType = MASS_LOSS_TYPE::OB;                                                 // set default dominant mass loss type
    double dMdt                   = BaseStar::CalculateMLRateOB(OPTIONS->OBMassLossPrescription());     // set default mass loss rate (OB)

    if (OPTIONS->ScaleCHEMassLossWithSurfaceHeliumAbundance()) {                                        // transition between OB and WR mass loss rates?
                                                                                                        // yes
        const double dMdtWR     = HeMS::CalculateMLRateMerritt2025_Static(p_Metallicity, p_Luminosity, p_Temperature, p_SigmaHurley, p_ZetaAnders); // WR mass loss rate
        const double OBfraction = CalculateMassLossFractionOB(p_HeAbundanceSurface);                    // mass loss fraction attributable to OB mass loss

        dMdt = (OBfraction * dMdt) + ((1.0 - OBfraction) * dMdtWR);                                     // combined mass loss rate

        if (OBfraction <= 0.5) dominantMLType = MASS_LOSS_TYPE::WR;                                     // determine dominant mass loss rate
    }

    // return mass loss rate, enhanced due to rotation, and dominant mass loss type
    return std::make_tuple(dMdt * CalculateMLRateEnhancementRotation(), dominantMLType);
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

