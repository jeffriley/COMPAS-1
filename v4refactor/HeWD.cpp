#include "HeWD.h"
#include "Star.h"
#include "HeMS.h"







///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




/* For HeWD, calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the retention efficiency parameter
 *
 *
 * For He WDs, we calculate the mass accretion rate following the 
 * StarTrack prescription (Belczynski+ 2008, sect 5.7.1).
 * https://ui.adsabs.harvard.edu/abs/2008ApJS..174..223B/abstract
 *
 *
 * Dbl_DblT CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate of the donor
 * @param   [IN]    p_IsHeRich                  Material is He-rich or not
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate (Msun/yr) and Retention Efficiency Parameter
 */
Dbl_DblT HeWD_Constituent::CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                                       const double p_AccretorMassRate,
                                                       const bool   p_IsHeRich) {

    const ACCRETION_REGIME regime = DetermineAccretionRegime(p_DonorMassRate, p_IsHeRich);
    CurrentState().SetAccretionRegime(regime);
                                                                               
    double acceptanceRate   = 0.0;                                                                                      // acceptance mass rate - default = 0.0
    double fractionAccreted = regime == ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES ? 0.0 : 1.0;              // accretion fraction - default = 1.0, but flashes restrict accumulation

    return std::make_tuple(acceptanceRate, fractionAccreted);
}


Dbl_DblT HeWD::CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich) {
    (void)p_DonorMassRate; (void)p_IsHeRich;
    return std::make_tuple(0.0, 0.0);
}


/* 
 * Determine the WD accretion regime based on the MT rate and whether the donor is He rich. Also,
 * initialize Sub-Chandrasekhar SN Ia or rejuvenation (evolution into HeMS) when necessary, by 
 * setting CurrentState().SetIsSubChandrasekharTypeIa() or SetShouldRejuvenate() (respectively). 
 *
 * The accretion regime is one of the options listed in enum ACCRETION_REGIME (constants.h)
 *
 * ACCRETION_REGIME DetermineAccretionRegime(const double p_DonorMassLossRate, const bool p_HeRich) 
 *
 * @param   [IN]    p_DonorMassRate      Donor mass loss rate, in units of Msol / Myr
 * @param   [IN]    p_HeRich             Whether the accreted material is helium-rich or not
 * @return                               WD accretion regime
 */
ACCRETION_REGIME HeWD_Constituent::DetermineAccretionRegime(const double p_DonorMassRate, const bool p_HeRich) {


    const double mass         = m_Star->Mass();
    const double heShell      = m_Star->HeShell();
    const double l0Ritter     = m_Star->L0Ritter();
    const double lambdaRitter = m_Star->LambdaRitter();
    const double xRitter      = m_Star->XRitter();

    double Mdot = p_DonorMassRate / MYR_TO_YEAR;                                                        // Accreted mass rate (M_sun/yr)

    ACCRETION_REGIME regime;
    if (p_HeRich) {
        if (utils::Compare(Mdot, HEWD_HE_MDOT_CRIT) <= 0) {
            regime           = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR;           // Could lead to Sub-Chandrasekhar SN Ia
            double massSubCh = WD_BELCZYNSKI_SN_CONSTANT - WD_BELCZYNSKI_SN_LINEAR * Mdot;              // Minimum mass for Sub-Chandrasekhar Mass detonation. Eq 62, Belczynski+ 2008.
            if (utils::Compare(mass, massSubCh) >= 0 ) {
                CurrentState().SetIsSubChandrasekharTypeIa(true);
            }
        } 
        else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION;                              // Could lift degeneracy and evolve into He MS. Requires minimum mass ! on top of the shell size
            if (utils::Compare(mass, HEWD_MINIMUM_MASS_IGNITION) >= 0) {
                if (utils::Compare(Mdot, WD_BELCZYNSKI_IMMEDIATE_FLASH) < 0) {                          // Accretion limit from eq 61, Belczynski+ 2008.
                    double mCritHeShell = WD_BELCZYNSKI_MINIMUM_HE_CONSTANT - WD_BELCZYNSKI_MINIMUM_HE_LINEAR * Mdot; // Minimum shell mass of He for ignition. Eq 61, Belczynski+ 2008. This helium should not be burnt, but not implemented this yet. Ruiter+ 2014.
                    if (utils::Compare(heShell, mCritHeShell) >= 0) {
                        CurrentState().SetShouldRejuvenate(true);
                    }
                } 
                else {
                    CurrentState().SetShouldRejuvenate(true);
                }
            }
        }
    } 
    else {
        double Mcrit = l0Ritter * PPOW(mass, lambdaRitter) / (xRitter * Q_HYDROGEN_BURNING);            // Eq. 60 in Belczynski+ 2008. 6e18 is the energy yield of H burning in ergs/g.
        if (utils::Compare(Mdot, Mcrit) <= 0) {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES;                             // Flashes restrict accumulation
        } 
        else {
            regime = ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION;                        // Material piles up on the WD. Leads to merger or CEE.
        }
    }

    return regime;
}


ACCRETION_REGIME HeWD::DetermineAccretionRegime(const double p_DonorMassRate, const bool p_HeRich) {
    (void)p_DonorMassRate; (void)p_HeRich;
    return ACCRETION_REGIME::ZERO;
}


/*
 * Specifies next stage, if the star changes its phase.
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                               Stellar Type for next phase
 */
STELLAR_TYPE HeWD::EvolveToNextPhase() {

    const double coreMass = Mass();
    m_InterimState.SetCoreMass  (coreMass);
    m_InterimState.SetRadius    (HeMS::CalculateRadiusAtZAHeMS_Hurley2000(coreMass));
    m_InterimState.SetLuminosity(HeMS::CalculateLuminosityAtZAHeMS_Hurley2000_Static(coreMass));
    m_InterimState.SetTau       (0.0);

    return STELLAR_TYPE::NAKED_HELIUM_STAR_MS;
}







