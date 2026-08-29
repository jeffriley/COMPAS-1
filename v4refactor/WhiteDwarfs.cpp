#include "WhiteDwarfs.h"
#include "NS.h"
#include "Star.h"








/// WhiteDwarfs Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//




/* Calculate eta_hydrogen from Claeys+ 2014, appendix B. This parameter depends 
 * on three regimes for the mass transfer rate, which here are distinguished by the 
 * thresholds logMdotUppH and logMdotLowH. In Claeys+ 2014, the mass transfer rate is
 * \dot{M}_{tr} and the thresholds are \dot{M}_{cr,H} and \dot{M}_{cr,H}/8, respectively. 
 *
 * However, we have used improved thresholds from Nomoto+ 2007 in which the 
 * lower boundary is \dot{M}_{stable} and the upper boundary is \dot{M}_{RG}. 
 * More precisely, we implemented quadratic fits to the values in Nomoto+ 2007,
 * table 5, as described in the second COMPAS methods paper (in prep). 
 *
 * double CalculateEtaH(const double p_MassTransferRate)
 *
 * @param   [IN]    p_MassTransferRate     Mass transfer rate onto the WD surface (Msun/Myr)
 * @return                                 Hydrogen accretion efficiency
 */
double WhiteDwarfs::CalculateEtaH(const double p_MassTransferRate) {

    double etaH = 0.0;                                      // default return value

    double logMassTransferRate = std::log10(p_MassTransferRate / MYR_TO_YEAR);
    double m_Mass_2            = Mass() * Mass();

    // The following coefficients come from quadratic fits to Nomoto+ 2007 results (table 5) in Mass vs log10 Mdot space, to cover the low-mass end.
    double logMdotUppH = WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_0 + WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_1 * Mass() + WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_2 * m_Mass_2; 
    double logMdotLowH = WD_LOG_MT_LIMIT_NOMOTO_STABLE_0   + WD_LOG_MT_LIMIT_NOMOTO_STABLE_1   * Mass() + WD_LOG_MT_LIMIT_NOMOTO_STABLE_2   * m_Mass_2;
    
    if (utils::Compare(logMassTransferRate, logMdotUppH) >= 0) {
        etaH = PPOW(10, logMdotUppH - logMassTransferRate);
    } 
    else if (utils::Compare(logMassTransferRate, logMdotLowH) >= 0) {
        etaH = 1.0;
    } 

    return etaH;
}


/* Calculate eta_helium from Claeys+ 2014, appendix B. Similarly to CalculateEtaH
 * above, this parameter depends on four regimes for the mass transfer rate, distinguished
 * here by logMdotUppHe, logMdotMidHe, and logMdotLowHe. In Claeys+ 2014, these thresholds
 * are \dot{M}_{up}, \dot{M}_{cr,He}, and \dot{M}_{low}, respectively. 
 *
 * However, we have again updated the thresholds to those described in Piersanti+ 2014,
 * table A1. The thresholds here are named by the boundaries RG/SS, SS/MF, and SF/Dt, 
 * respectively (see text for details). Note that the different flashes regimes from 
 * Piersanti+ 2014 have been merged into one, i.e we omit the MF/SF boundary, and 
 * the accumulation regime has been changed so we can get double detonations. Finally, 
 * eta_KH04 has also been updated with the accretion efficiency values from Piersanti+ 2014.
 *
 * double CalculateEtaHe(const double p_MassTransferRate)
 *
 * @param   [IN]    p_MassTransferRate     Mass transfer rate onto the WD surface (Msun/Myr)
 * @return                                 Helium accretion efficiency
 */
double WhiteDwarfs::CalculateEtaHe(const double p_MassTransferRate) {

    double etaHe = 1.0;                                     // default return value - so we can have double detonations
    
    double logMassTransferRate = std::log10(p_MassTransferRate / MYR_TO_YEAR);

    // The following coefficients in massTransfer limits come from table A1 in Piersanti+ 2014.
    double logMdotUppHe = WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_0 + WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_1 * Mass();
    double logMdotMidHe = WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_0 + WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_1 * Mass();
    double logMdotLowHe = WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_0 + WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_1 * Mass();

    if (utils::Compare(logMassTransferRate, logMdotUppHe) >= 0) {
        etaHe = PPOW(10, logMdotUppHe - logMassTransferRate);
    } 
    else if (utils::Compare(logMassTransferRate, logMdotMidHe) >= 0) {
        etaHe = 1.0;
    } 
    else if (utils::Compare(logMassTransferRate, logMdotLowHe) >= 0) {
        etaHe = CalculateEtaPTY(p_MassTransferRate);
    } 

    return etaHe;
}


/* Calculate accretion efficiency as indicated in Piersanti+ 2014, section A3. Their recipe works
 * for specific mass and Mdot values, so a better implementation would require interpolation and
 * extrapolation (specially towards the low-mass end). Right now, we just adopt a
 * piece-wise approach. Note that the authors also specify that this is based on the first
 * strong flash only, but we use it for all episodes.
 *
 * double CalculateEtaPTY(const double p_MassTransferRate)
 *
 * @param   [IN]    p_MassTransferRate     Mass transfer rate onto the WD surface (Msun/Myr)
 * @return                                 Accretion efficiency during the first stron helium flash, Piersanti+ 2014
 */
double WhiteDwarfs::CalculateEtaPTY(const double p_MassTransferRate) {

    double etaPTY = 0.0;                            // default return value

    double massRate   = p_MassTransferRate * 100;   // Piersanti + 2014 assumes the rate is in units of 10^-8 Msun / yr (see Eq. A3)
    double massRate_2 = massRate * massRate;
    double massRate_3 = massRate_2 * massRate;

    // Limits on each conditional statement come from masses from each model in Piersanti+ 2014. The final etaPTY value is based on table A3.
    if (utils::Compare(Mass(), 0.6) <= 0) {
        etaPTY = WD_PIERSANTI_M060_G0 + WD_PIERSANTI_M060_G1 * massRate + WD_PIERSANTI_M060_G2 * massRate_2 - WD_PIERSANTI_M060_G3 * massRate_3;
    } 
    else if  (utils::Compare(Mass(), 0.7) <= 0) {
        etaPTY = -WD_PIERSANTI_M070_G0 + WD_PIERSANTI_M070_G1 * massRate - WD_PIERSANTI_M070_G2 * massRate_2 + WD_PIERSANTI_M070_G3 * massRate_3;
    } 
    else if (utils::Compare(Mass(), 0.81) <= 0) {
        etaPTY = WD_PIERSANTI_M081_G0 + WD_PIERSANTI_M081_G1 * massRate + WD_PIERSANTI_M081_G2 * massRate_2 - WD_PIERSANTI_M081_G3 * massRate_3;
    } 
    else if (utils::Compare(Mass(), 0.92) <= 0) { 
        etaPTY = -WD_PIERSANTI_M092_G0 + WD_PIERSANTI_M092_G1 * massRate + WD_PIERSANTI_M092_G2 * massRate_2 - WD_PIERSANTI_M092_G3 * massRate_3;
    } 
    else {
        etaPTY = -WD_PIERSANTI_M102_G0 + WD_PIERSANTI_M102_G1 * massRate - WD_PIERSANTI_M102_G2 * massRate_2 + WD_PIERSANTI_M102_G3 * massRate_3;
    }

    return std::clamp(etaPTY, 0.0, 1.0);    // Clamp eta to [0.0, 1.0]
}










/// WhiteDwarfs_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//



/* Calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the retention efficiency parameter
 *
 * Currently used for COWDs and ONeWDs
 *
 * For a given mass transfer rate, this function computes the amount of mass a WD would retain after
 * flashes, as given by appendix B of Claeys+ 2014. 
 * https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..83C/abstract 
 *
 *
 * Dbl_DblT CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate from the donor
 * @param   [IN]    p_IsHeRich                  Material is He-rich or not
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate (Msun/yr) and Retention Efficiency Parameter
 */
Dbl_DblT WhiteDwarfs_Constituent::CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                              const double p_AccretorMassRate,
                                                              const bool   p_IsHeRich) {

    const ACCRETION_REGIME regime = DetermineAccretionRegime(p_DonorMassRate, p_IsHeRich);
    CurrentState().SetAccretionRegime(regime);

    double acceptanceRate   = 0.0;                                                       // acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                       // accretion fraction - default = 0.0

    acceptanceRate = p_DonorMassRate * m_Star->CalculateEtaHe(p_DonorMassRate);
    if (!p_IsHeRich) acceptanceRate *= m_Star->CalculateEtaH(p_DonorMassRate);

    fractionAccreted = acceptanceRate / p_DonorMassRate;

    return std::make_tuple(acceptanceRate, fractionAccreted);
}





/* 
 * Determine the WD accretion regime based on the MT rate and whether the donor is He rich. Also,
 * initialize He-Shell detonation or Off-centre ignition when necessary, by setting
 * CurrentState().SetHeShellDetonation() or CurrentState().SetOffCentreIgnition() (respectively).
 *
 * The accretion regime is one of the options listed in enum ACCRETION_REGIME (constants.h)
 *
 * Note that we have merged the different flashes regimes from Piersanti+ 2014 into a single regime.
 *
 * ACCRETION_REGIME DetermineAccretionRegime(const double p_DonorMassLossRate, const bool p_HeRich) 
 *
 * @param   [IN]    p_DonorMassLossRate      Donor mass loss rate, in units of Msol / Myr
 * @param   [IN]    p_HeRich                 Whether the accreted material is helium-rich or not
 * @return                                   Current WD accretion regime
 */
ACCRETION_REGIME WhiteDwarfs_Constituent::DetermineAccretionRegime(const double p_DonorMassLossRate, const bool p_HeRich) {

    const double mass    = m_Star->Mass();
    const double heShell = m_Star->HeShell();

    double logMdot          = std::log10(p_DonorMassLossRate / MYR_TO_YEAR);                                                // logarithm of the accreted mass (M_sun/yr)
    ACCRETION_REGIME regime = ACCRETION_REGIME::ZERO;

    if (p_HeRich) {
        // The following coefficients in logMassTransfer limits come from table A1 in Piersanti+ 2014.
        double logMassTransferCrit       = WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_0 + WD_LOG_MT_LIMIT_PIERSANTI_RG_SS_1 * mass;
        double logMassTransferStable     = WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_0 + WD_LOG_MT_LIMIT_PIERSANTI_SS_MF_1 * mass;    // Piersanti+2014 has several Flashes regimes. Here we group them into one.
        double logMassTransferDetonation = WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_0 + WD_LOG_MT_LIMIT_PIERSANTI_SF_Dt_1 * mass;    // critical value for double detonation regime in Piersanti+ 2014
        if (utils::Compare(logMdot, logMassTransferStable) < 0) {
            if (utils::Compare(logMdot, logMassTransferDetonation) > 0) {
                regime = ACCRETION_REGIME::HELIUM_FLASHES;
            } 
            else {
                regime = ACCRETION_REGIME::HELIUM_ACCUMULATION;
                if ((utils::Compare(mass, MASS_DOUBLE_DETONATION_CO) >= 0) && (utils::Compare(heShell, WD_HE_SHELL_MCRIT_DETONATION) >= 0)) {
                    CurrentState().SetHeShellDetonation(true);
                }
            }
        } 
        else if (utils::Compare(logMdot, logMassTransferCrit) > 0) {
            regime = ACCRETION_REGIME::HELIUM_OPT_THICK_WINDS;
        } 
        else {
            regime = ACCRETION_REGIME::HELIUM_STABLE_BURNING;
            if ((utils::Compare(logMdot, COWD_LOG_MDOT_MIN_OFF_CENTER_IGNITION) > 0) && (utils::Compare(mass, COWD_MASS_MIN_OFF_CENTER_IGNITION) > 0)) {  // JR FIX THIS - add British spelling!!! <<<<<<<<<<<<<<<<<<<
                CurrentState().SetOffCentreIgnition(true);
            }
        }
    } 
    else {
        // The following coefficients in logMassTransfer limits come from quadratic fits to Nomoto+ 2007 results (table 5) in Mass vs log10 Mdot space, to cover the low-mass end.
        double mass_2 = mass * mass;
        double logMassTransferCrit   = WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_0 + WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_1 * mass + WD_LOG_MT_LIMIT_NOMOTO_REDGIANT_2 * mass_2;
        double logMassTransferStable = WD_LOG_MT_LIMIT_NOMOTO_STABLE_0   + WD_LOG_MT_LIMIT_NOMOTO_STABLE_1   * mass + WD_LOG_MT_LIMIT_NOMOTO_STABLE_2   * mass_2;

        if (utils::Compare(logMdot, logMassTransferStable) < 0) {
            regime = ACCRETION_REGIME::HYDROGEN_FLASHES;
        } 
        else if (utils::Compare(logMdot, logMassTransferCrit) > 0) {
            regime = ACCRETION_REGIME::HYDROGEN_OPT_THICK_WINDS;
        } 
        else {
            regime = ACCRETION_REGIME::HYDROGEN_STABLE_BURNING;
        }
    }

    return regime;
}


/* 
 * Increase shell size after mass transfer episode. Hydrogen and helium shells are kept separately.
 * Only applies the full mass increase from accretion to one of the shells. Does not account for, e.g,
 * the H layer burning and building up the He layer, which may be desired in the future. - RTW 9/14/22
 *
 * void ResolveShellChange(const double p_AccretedMass)
 *
 * @param   [IN]    p_AccretedMass              Mass accreted
 */
void WhiteDwarfs_Constituent::ResolveShellChange(const double p_AccretedMass) {
    

    switch (CurrentState().AccretionRegime()) {

        case ACCRETION_REGIME::HELIUM_ACCUMULATION:
        case ACCRETION_REGIME::HELIUM_FLASHES:
        case ACCRETION_REGIME::HELIUM_STABLE_BURNING:
        case ACCRETION_REGIME::HELIUM_OPT_THICK_WINDS:
        case ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR:
        case ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION:
	        m_Star->SetHeShell(m_Star->HeShell() + p_AccretedMass);
            break;

        case ACCRETION_REGIME::HYDROGEN_FLASHES:
        case ACCRETION_REGIME::HYDROGEN_STABLE_BURNING:
        case ACCRETION_REGIME::HYDROGEN_OPT_THICK_WINDS:
        case ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES:
        case ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION:
	        m_Star->SetHShell(m_Star->HShell() + p_AccretedMass);
            break;

        case ACCRETION_REGIME::NONE:    // DEPRECATED June 2024 - remove end 2024 
        case ACCRETION_REGIME::ZERO:
            SHOW_WARN(ERROR::UNEXPECTED_ACCRETION_REGIME, "No mass added to shell");        // show warning
            break;

        default:                                                                            // unknown stellar population
            THROW_ERROR(ERROR::UNKNOWN_ACCRETION_REGIME);                                   // throw error
    }
}




/*
 * Resolve Accretion-Induced Collapse of a WD
 *
 * Following Hurley et al. 2000, Section 6.2.1
 *
 * An AIC of a WD results in a NS, which we are 
 * here assuming to have a low mass equal to the ECSN
 * remnant NS mass, and no natal kick. 
 *
 * STELLAR_TYPE ResolveAIC() 
 *
 * @return                                      Stellar type of remnant (STELLAR_TYPE::NEUTRON_STAR if SN, otherwise current type)
 */
STELLAR_TYPE WhiteDwarfs::ResolveAIC() { 

    if (!IsSupernova()) return m_StellarType;                                           // shouldn't be here if no SN

    m_SupernovaDetails.totalMassAtCOFormation  = m_Mass;
    m_SupernovaDetails.HeCoreMassAtCOFormation = m_HeCoreMass;
    m_SupernovaDetails.COCoreMassAtCOFormation = m_COCoreMass;
    m_SupernovaDetails.coreMassAtCOFormation   = m_CoreMass;
    SetSNHydrogenContent();                                                             // SN to be H-poor. 

    m_Mass                                = MECS_REM;                                   // defined in constants.h
    
    m_SupernovaDetails.drawnKickMagnitude = 0.0;
    m_SupernovaDetails.kickMagnitude      = 0.0;

    SetSNCurrentEvent(SN_EVENT::AIC);                                                   // AIC happening now
    SetSNPastEvent(SN_EVENT::AIC);                                                      // ... and will be a past event

    return STELLAR_TYPE::NEUTRON_STAR;
}


/*
 * Resolve Type 1a Supernova 
 *
 * A Type 1a SN results in a massless remnant.
 *
 * STELLAR_TYPE ResolveSNIa() 
 *
 * @return                                      Stellar type of remnant (STELLAR_TYPE::MASSLESS_REMNANT if SN, otherwise current type)
 */
STELLAR_TYPE WhiteDwarfs::ResolveSNIa() { 

    if (!IsSupernova()) return m_StellarType;                                           // shouldn't be here if no SN

    m_SupernovaDetails.totalMassAtCOFormation  = m_Mass;
    m_SupernovaDetails.HeCoreMassAtCOFormation = m_HeCoreMass;
    m_SupernovaDetails.COCoreMassAtCOFormation = m_COCoreMass;
    m_SupernovaDetails.coreMassAtCOFormation   = m_CoreMass;
    SetSNHydrogenContent();                                                             // SN to be H-poor. 
        
    m_Mass       = 0.0;
    m_Radius     = 0.0;
    m_Luminosity = 0.0;
    m_Age        = 0.0;

    m_SupernovaDetails.drawnKickMagnitude = 0.0;
    m_SupernovaDetails.kickMagnitude      = 0.0;

    SetSNCurrentEvent(SN_EVENT::SNIA);                                                  // SN Type Ia happening now
    SetSNPastEvent(SN_EVENT::SNIA);                                                     // ... and will be a past event

    return STELLAR_TYPE::MASSLESS_REMNANT;
}


/*
 * Resolve Double Detonation     
 *
 * A double detonation results in a massless remnant.
 *
 * STELLAR_TYPE ResolveHeSD() 
 *
 * @return                                      Stellar type of remnant (STELLAR_TYPE::MASSLESS_REMNANT if SN, otherwise current type)
 */
STELLAR_TYPE WhiteDwarfs::ResolveHeSD() { 

    if (!IsSupernova()) return m_StellarType;                                           // shouldn't be here if no SN

    m_Mass       = 0.0;
    m_Radius     = 0.0;
    m_Luminosity = 0.0;
    m_Age        = 0.0;
    
    m_SupernovaDetails.drawnKickMagnitude = 0.0;
    m_SupernovaDetails.kickMagnitude      = 0.0;

    SetSNCurrentEvent(SN_EVENT::HeSD);                                                  // SN Type Ia (HeSD) happening now
    SetSNPastEvent(SN_EVENT::HeSD);                                                     // ... and will be a past event

    return STELLAR_TYPE::MASSLESS_REMNANT;
}
