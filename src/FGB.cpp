#include "FGB.h"
#include "HeMS.h"
#include "HeWD.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCoreMass_Hurley2000
 *
 * @brief
 * Calculate core mass on the First Giant Branch, per Hurley et al. 2000, eqs 39 & 45
 *
 *
 * double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Age, const double p_Tau, const DblVectorT& p_GBparams, const DblVectorT& p_tScales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_Tau                           Phase-relative age of the star (Myr)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      FGB core mass (Msol)
 */
double FGB::CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Age, const double p_Tau, const DblVectorT& p_GBparams, const DblVectorT& p_tScales) const {

    double McGB;
    if (p_Mass < GLOBALS->HurleyMassCutoffs(static_cast<int>(HURLEY_MCO::MHeF))) {
        const double AH = p_GBparams[HURLEY_GBP::AH];
        const double B  = p_GBparams[HURLEY_GBP::B];
        const double D  = p_GBparams[HURLEY_GBP::D];
        const double p  = p_GBparams[HURLEY_GBP::P];
        const double q  = p_GBparams[HURLEY_GBP::Q];

        McGB = p_Age <= p_tScales[HURLEY_TS::MX_FGB]
                ? PPOW(((p - 1.0) * AH * D * (p_tScales[HURLEY_TS::INF1_FGB] - p_Age)), (1.0 / (1.0 - p)))
                : PPOW(((q - 1.0) * AH * B * (p_tScales[HURLEY_TS::INF2_FGB] - p_Age)), (1.0 / (1.0 - q)));
    }
    else {
        const double McBGB = p_GBparams[HURLEY_GBP::MCBGB];
        McGB = McBGB + ((CalculateCoreMassAtHeI_Hurley2000(p_Mass) - McBGB) * p_Tau);
    }

    return McGB;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////







///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Modify the star after it loses its envelope
 *
 * Hurley et al. 2000, section 6 just before eq 76 (see also after Eq. 105)
 *
 * Where necessary updates attributes of star (depending upon stellar type):
 *
 *     - m_StellarType
 *     - m_Timescales
 *     - m_GBparams
 *     - m_Luminosity
 *     - m_Radius
 *     - m_Mass
 *     - m_Mass0
 *     - m_CoreMass
 *     - m_HeCoreMass
 *     - m_COCoreMass
 *     - m_Age
 *
 * STELLAR_TYPE ResolveEnvelopeLoss(bool p_Force)
 *
 * @param   [IN]    p_Force                     Boolean to indicate whether the resolution of the loss of the envelope should be performed
 *                                              without checking the precondition(s).
 *                                              Default is false.
 *
 * @return                                      Stellar Type to which star should evolve after losing envelope
 */
STELLAR_TYPE FGB::ResolveEnvelopeLoss(bool p_Force) {

    STELLAR_TYPE stellarType = m_StellarType;

    if (p_Force || utils::Compare(m_CoreMass, m_Mass) >= 0) {                                       // Envelope loss

        m_Mass       = std::min(m_CoreMass, m_Mass);
        m_CoreMass   = m_HeCoreMass;
        m_Mass       = m_CoreMass;
        m_COCoreMass = 0.0;
        
        if (utils::Compare(m_Mass0, massCutoffs(MHeF)) < 0) {                                       // Star evolves to Helium White Dwarf

            stellarType = STELLAR_TYPE::HELIUM_WHITE_DWARF;

            m_Age       = 0.0;
            m_Radius    = WhiteDwarfs::CalculateRadius_Hurley2000_Static(m_Mass);
        }
        else {                                                                                      // Star evolves to Zero age Naked Helium Main Star

            stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;

            m_Mass0     = m_Mass;
            m_Age       = 0.0;
            m_Radius    = HeMS::CalculateRadiusAtZAHeMS_Hurley2000_Static(m_Mass);
        }
    }

    return stellarType;
}


/*
 * Modify the star due to (possible) helium flash
 *
 *
 * void ResolveHeliumFlash()
 *
 * Deletermine if Helium Flash occurs, and if so set m_Mass0 equal to current mass as described in Hurley+ (2000), last paragraph before start of 7.1.1.
 */
void FGB::ResolveHeliumFlash() {

    if (utils::Compare(m_Mass0, massCutoffs(MHeF)) < 0) {               // Helium flash if initial mass < Helium Flash cutoff
        m_Mass0 = m_Mass;                                               // for LM star at ZAHB (end of GB/beginning of CHeB) due to helium flash when doing mass loss
    }
}


/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE FGB::EvolveToNextPhase() {
    ResolveHeliumFlash();   // ...before we evolve

    return STELLAR_TYPE::CORE_HELIUM_BURNING;
}
