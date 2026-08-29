#include "HeHG.h"
#include "WhiteDwarfs.h"
#include "Star.h"
//#include "HeGB.h" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#include "COWD.h"
//#include "HeWD.h" <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              CONSTANTS / PARAMETERS                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                   MISCELLANEOUS                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * DetermineEnvelopeType
 *
 * @brief
 * Determine the star's envelope type, based on the user-specified ENVELOPE_STATE_PRESCRIPTION.
 * 
 *
 * ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Envelope type (ENVELOPE::{RADIATIVE, CONVECTIVE, REMNANT})
 */
COMPAS_PURE ENVELOPE HeHG::DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const {
    
    ENVELOPE envType;
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                             // which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::CONVECTIVE_MASS_FRACTION:             // CONVECTIVE_MASS_FRACTION
            // envelope is convective when the convective mass exceeds specified fraction of the envelope mass
            double mEnv;
            std::tie(mEnv, std::ignore) = CalculateConvectiveEnvelopeMass();
            envType = mEnv / (p_Mass - p_CoreMass) > OPTIONS->ConvectiveEnvelopeMassThreshold() ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:                    // FIXED_TEMPERATURE
            // envelope is radiative if temperature exceeds specified threshold, otherwise convective
            envType = (p_Temperature * TSOL) > OPTIONS->ConvectiveEnvelopeTemperatureThreshold() ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:                               // HURLEY
            // Hurley et al. 2002 eqs 39 & 40, and Hurley et al. 2000 end of section 7.2 describe gradual
            // growth of convective envelope over HG, but we approximate it as already convective here
            envType = ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:                               // COMPAS LEGACY
            envType = ENVELOPE::RADIATIVE;                                      // always radiative
            break;

        default:                                                                // unknown prescription
            // the only way this can happen is if someone added an ENVELOPE_STATE_PRESCRIPTION and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION);            // throw error                
    }
    
    return envType;
}







/*
 * Determine if evolution should continue on this phase, or whether evolution
 * on this phase should end (and so evolve to next phase)
 *
 *
 * bool ShouldEvolveOnPhase()
 *
 * @return                                      Boolean flag: true if evolution on this phase should continue, false if not
 */
bool HeHG::ShouldEvolveOnPhase() const {

    double McMax = CalculateMaximumCoreMass(m_Mass);
    double McSN  = CalculateCoreMassAtSupernova_Static(MECS, m_GBParams[static_cast<int>(GBP::McBAGB)]);
    return ((utils::Compare(COCoreMass(), McMax) <= 0 || utils::Compare(McMax, McSN) >= 0) && !ShouldEnvelopeBeExpelledByPulsations());    // Evolve on HeHG phase if McCO <= McMax or McMax >= McSN and envelope is not ejected by pulsations
}


/*
 * Modify the star after it loses its envelope
 *
 * Hurley et al. 2000, section 6 just before eq 76 and after Eq. 105
 *
 * Where necessary updates attributes of star (depending upon stellar type):
 *
 *     - m_StellarType
 *     - m_Timescales
 *     - m_GBParams
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
STELLAR_TYPE HeHG::ResolveEnvelopeLoss(bool p_Force) {

    STELLAR_TYPE stellarType = m_StellarType;
    
    if (ShouldEnvelopeBeExpelledByPulsations()) m_EnvelopeJustExpelledByPulsations = true;
    
    if (p_Force || utils::Compare(m_CoreMass, m_Mass) >= 0 || m_EnvelopeJustExpelledByPulsations) {         // Envelope lost - determine what type of star to form

        m_Mass       = std::min(m_CoreMass, m_Mass);
        m_CoreMass   = m_Mass;
        m_HeCoreMass = m_Mass;
        m_COCoreMass = m_Mass;
        m_Mass0      = m_Mass;
        m_Radius     = COWD::CalculateRadiusOnPhase_Static(m_Mass);
        m_Age        = 0.0;
        if (!IsSupernova()) {
            // note that this uses the CO core mass, rather than the core mass at base of AGB or He mass at He star birth suggested by Hurley+, 2000
            stellarType = (utils::Compare(m_COCoreMass, OPTIONS->MCBUR1() ) < 0) ? STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF : STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
        }
    }
    return stellarType;
}


/*
 * Determine if star should continue evolution as a supernova
 *
 *
 * bool IsSupernova()
 *
 * @return                                      Boolean flag: true if star has gone Supernova, false if not
 */
bool HeHG::IsSupernova() const {
    if (utils::Compare(CoreMass(), Mass()) == 0) {      // special case of ultra-stripped-star -- go SN immediately if over ECSN limit
        return (utils::Compare(Mass(), MECS) > 0);
    }
        
    return (utils::Compare(COCoreMass(), CalculateCoreMassAtSN(MECS, CurrentState().HurleyGBParams(HURLEY_GB_PARAMETERS::McBAGB))) >= 0); // Go supernova if CO core mass large enough
}

/*
 * Assistant function for determining the supernova explosion type
 *
 *
 * double       CalculateInitialSupernovaMass()
 *
 * @return                                      double: Initial supernova supernova mass variable
 */
double HeHG::CalculateInitialSupernovaMass() const {
    if (utils::Compare(CoreMass(), Mass()) == 0) {      // special case of ultra-stripped-star -- use current mass
        return std::max(Mass(), CurrentState().HurleyGBParams(HURLEY_GB_PARAMETERS::McBAGB));
    }
    return GiantBranch::CalculateInitialSupernovaMass();
}

/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE HeHG::EvolveToNextPhase() {
    return STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF;
}





///////////////////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                LAMBDA CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the common envelope lambda parameter using the "Nanjing" prescription
 * from X.-J. Xu and X.-D. Li arXiv:1004.4957 (v1, 28Apr2010) as implemented in STARTRACK
 *
 * This implementation adapted from the STARTRACK implementation (STARTRACK courtesy Chris Belczynski)
 *
 * This function is for HeHG and HeGB stars (for Helium stars: always use Natasha Ivanova's fit)
 *
 *
 * double CalculateLambdaNanjingStarTrack(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass
 * 
 * @return                                      Nanjing lambda for use in common envelope
 */
COMPAS_PURE double HeHG_Constituent::CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const {
    (void)p_Mass;                                       // unused: HeHG StarTrack lambda is purely radius-driven
    (void)p_CoreMass;                                   // unused: same as above

    const double radius = p_Radius;

    double rMin = 0.25;                              // minimum considered radius: Natasha
	double rMax = 120.0;                             // maximum considered radius: Natasha

	double rMinLambda = 0.3 * PPOW(rMin, -0.8);
	double rMaxLambda = 0.3 * PPOW(rMax, -0.8);

	return radius < rMin ? rMinLambda : (radius > rMax ? rMaxLambda : 0.3 * PPOW(radius, -0.8));
}


