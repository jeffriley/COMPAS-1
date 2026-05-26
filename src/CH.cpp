#include "CH.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
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
 * DblVectorT CalculateTimescales_Hurley2000(const double p_Mass, const DblVectorT& p_GBparams, const DblVectorT& p_Timescales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
DblVectorT CH::CalculateTimescales_Hurley2000(const double p_Mass, const DblVectorT& p_GBparams, const DblVectorT& p_tScales) const {

    DblVectorT tScales = p_tScales;                     // copy given timescales

    // (re)calculate MS timescales
    tScales = MainSequence::CalculateTimescales_Hurley2000(p_Mass, p_GBparams, tScales);

    // enhamce lifetimes as appropriate
    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {   // enhance lifetime of CH stars?
                                                        // yes
        const double lifetimesRatio = CalculateLifetimesRatio_Szecsi2020(p_Mass);

        tScales[HURLEY_TS::BGB] *= lifetimesRatio;
        tScales[HURLEY_TS::MS]  *= lifetimesRatio;
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
double CH::CalculateLuminosity_Hurley2000(
    const double p_Mass,
    const double p_Tau,
    const double p_Age,
    const double p_LZAMS,
    const double p_tMS,
    const double p_tBGB
) const {

    // unenhanced CH luminosity is just MS luminosity
    double luminosity = MainSequence::CalculateLuminosity_Hurley2000(p_Mass, p_Time, p_LZAMS, p_tMS, p_tBGB);

    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {                   // enhance luminosity of CH stars?
                                                                        // yes
        // enhancement should not reduce luminosity, so ratio is clamped to a minimum of +1.0
        // enhancement amount grows from 1 to logLuminosityRatio over main-sequence
        const double enhancement = 1.0 + (std::max(CalculateLogLuminositiesRatio(mass), 1.0) - 1.0) * p_Tau * p_Tau;

        luminosity = PPOW(10.0, std::log10(luminosity) * enhancement);  // apply enhancement
    }

    return luminosity;
}


/*
 * CalculateMLrate_Belczynski2010
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, per Belczynski 2010
 * (as implemented in StarTrack - courtesy Chris Belczynski).
 *
 * If option `--scale-mass-loss-with-surface-helium-abundance` was specified,
 * the mass loss rate will be scaled with the surface helium abundance.
 *
 * If option `--enable-rotationally-enhanced-mass-loss` was specified,
 * the mass loss rate will be enhance due to rotation.
 * 
 * 
 * MassLossT CalculateMLrate_Belczynski2010(
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
MassLossT CH::CalculateMLrate_Belczynski2010(
    const double p_Mass,
    const double p_Radius,
    const double p_Luminosity,
    const double p_Temperature,
    const double p_PerturbationMu,
    const double p_HeAbundanceSurface
) const {

    // calculate OB mass loss rate
    double dMdt;
    MASS_LOSS_TYPE dominantMLtype;
    std::tie(dMdt, dominantMLtype) = BaseStar::CalculateMLrate_Belczynski2010(p_Mass, p_Radius, p_Luminosity, p_Temperature, p_PerturbationMu, p_HeAbundanceSurface);

    // scale mass loss with the surface helium abundance if required
    // (transition between OB and WR mass loss rates)
    if (OPTIONS->ScaleCHEMassLossWithSurfaceHeliumAbundance()) {

        double fractionOB = CalculateMLfractionOB(p_HeAbundanceSurface);

        if (fractionOB < 1.0) {
            double dMdtWR = 0.0;  
            MASS_LOSS_TYPE dominantMLtypeWR;
            std::tie(dMdtWR, dominantMLTypeWR) = CalculateMLrateWR_ZDependent_Static(p_Luminosity, 0.0); // *Ilya* should we use p_PerturbationMu here (since we have it)?
        
            dMdt   *= fractionOB;
            dMdtWR *= (1.0 - fractionOB);

            if (dMdtWR > dMdt) dominantMLtype = dominantMLTypeWR;   // dominant mass loss type
        
            dMdt += dMdtWR;                                         // scaled mass loss rate
        }
    }

    // enhance mass loss rate due to rotation if required
    if (OPTIONS->EnableRotationallyEnhancedMassLoss()) dMdt *= CalculateMLrateRotationEnhancement_Langer1998();

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrate_Merritt2025
 *
 * Calculate mass loss rate, and dominant mass loss type, at the current evolutionary phase,
 * per Merritt et al., 2025.
 * 
 * If option `--scale-mass-loss-with-surface-helium-abundance` was specified,
 * the mass loss rate will be scaled with the surface helium abundance.
 *
 * If option `--enable-rotationally-enhanced-mass-loss` was specified,
 * the mass loss rate will be enhance due to rotation.
 * 
 * 
 * MassLossT CalculateMLrate_Merritt2025(
 *     const double p_Mass,
 *     const double p_Radius,
 *     const double p_Luminosity,
 *     const double p_Temperature,
 *     [[maybe_unused]] const double p_PerturbationMu,
 *     [[maybe_unused]] const double p_mStart,
 *     const double p_HeAbundanceSurface
 * ) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu (not used here)
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol) (not used here)
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
MassLossT CH::CalculateMLrate_Merritt2025(
    const double p_Mass,
    const double p_Radius,
    const double p_Luminosity,
    const double p_Temperature,
    [[maybe_unused]] const double p_PerturbationMu,
    [[maybe_unused]] const double p_mStart,
    const double p_HeAbundanceSurface
) const {

    // calculate OB mass loss rate 
    double dMdt;
    MASS_LOSS_TYPE dominantMLtype;
    std::tie(dMdt, dominantMLtype) = BaseStar::CalculateMLrateOB(p_Mass, p_Luminosity, p_Temperature, OPTIONS->OBMassLossPrescription());

    // scale mass loss with the surface helium abundance if required
    // (transition between OB and WR mass loss rates)
    if (OPTIONS->ScaleCHEMassLossWithSurfaceHeliumAbundance()) {

        double fractionOB = CalculateMassLossFractionOB(p_HeAbundanceSurface);

        if (fractionOB < 1.0) {
            double dMdtWR = 0.0;  
            MASS_LOSS_TYPE dominantMLtypeWR;
            std::tie(dMdtWR, dominantMLTypeWR) = HeMS::CalculateMLrate_Merritt2025_Static(const double p_Luminosity, const double p_Temperature);

            dMdt   *= fractionOB;
            dMdtWR *= (1.0 - fractionOB);

            if (dMdtWR > dMdt) dominantMLtype = dominantMLTypeWR;   // dominant mass loss type
        
            dMdt += dMdtWR;                                         // scaled mass loss rate
        }
    }

    // enhance mass loss rate due to rotation if required
    if (OPTIONS->EnableRotationallyEnhancedMassLoss()) dMdt *= CalculateMLrateRotationEnhancement_Langer1998();

    return std::make_tuple(dMdt, dominantMLtype);
}










/////////////////////////////////////////////////////////////////// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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

