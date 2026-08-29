#include "GiantBranch.h"
#include "HeMS.h"
#include "WhiteDwarfs.h"
#include "ONeWD.h"
#include "NS.h"
#include "BH.h"
#include "Star.h"



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                     COEFFICIENT AND CONSTANT CALCULATIONS ETC.                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
 * DblVectorT CalculateTimescales_Hurley2000(const double p_Mass, const DblVectorT& p_GBParams, const DblVectorT& p_Timescales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBParams                      Hurley GB parameters
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
COMPAS_PURE DBL_VECTOR GiantBranch::CalculateTimescales_Hurley2000(const double p_Mass, const DblVectorT& p_GBParams, const DblVectorT& p_tScales) const {

// #defines for convenience and readability - undefined at end of function
#define GBParams(x) p_GBParams[static_cast<int>(HURLEY_GB_PARAMETERS::x)]
#define tScales(x) tScales[static_cast<int>(HURLEY_TIMESCALES::x)]

    const double p1   = GBParams(p) - 1.0;
    const double q1   = GBParams(q) - 1.0;
    const double p1_p = p1 / GBParams(p);
    const double q1_q = q1 / GBParams(q);

    const double lBGB = GiantBranch::CalculateLuminosityAtBGB_Hurley2000_Static(p_Metallicity, p_Mass);

    DBL_VECTOR tScales = p_tScales; // copy given timescales


    tScales = MainSequence::CalculateTimescales_Hurley2000(p_Mass, tScales);

    tScales[static_cast<int>(HURLEY_TIMESCALES::Inf1_FGB)] = tScales(BGB) + ((1.0 / (p1 * GBParams(AH) * GBParams(D))) * PPOW((GBParams(D) / lBGB), p1_p));
    tScales[static_cast<int>(HURLEY_TIMESCALES::Mx_FGB)]   = tScales(Inf1_FGB) - ((tScales(Inf1_FGB) - tScales(BGB)) * PPOW((lBGB / GBParams(Lx)), p1_p));
    tScales[static_cast<int>(HURLEY_TIMESCALES::Inf2_FGB)] = tScales(Mx_FGB) + ((1.0 / (q1 * GBParams(AH) * GBParams(B))) * PPOW((GBParams(B) / GBParams(Lx)), q1_q));

    tScales[static_cast<int>(HURLEY_TIMESCALES::HeI)]      = CalculateLifetimeToHeI_Hurley2000(p_Mass, tScales(Inf1_FGB), tScales(Inf2_FGB), p_GBParams); // JR FIX THIS: parameter order????
    tScales[static_cast<int>(HURLEY_TIMESCALES::HeMS)]     = HeMS::CalculatePhaseLifetime_Hurley2000_Static(p_Mass);

    // return timescales vector by value - NRVO takes care of performance/efficiency
    return tScales;

#undef tScales
#undef GBParams
}



/*
 * CalculateGBParams_Hurley2000
 *
 * @brief
 * Calculate Giant Branch (GB) parameters, per Hurley et al. 2000.
 *
 * Since Giant Branch parameters depend on a star's mass, they need to be calculated whenever
 * the mass of the star changes (probably every timestep).
 *
 * p_GBParams is passed by value, modified in place, and returned by value.  No copy on the call,
 * and NRVO constructs GBParams directly in the caller's return slot, so no copy on the return.
 * 
 * 
 * DblVectorT CalculateGBParams_Hurley2000(const double p_Mass, DblVectorT p_GBParams) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBParams                      Hurley GB parameters
 * @return                                      Mutated GB parameters
 */
DblVectorT GiantBranch::CalculateGBParams_Hurley2000(const double p_Mass, DblVectorT& p_GBParams) const {

    p_GBParams[HURLEY_GBP::AH]     = CalculateHRateConstant_Hurley2000(p_Mass);
    p_GBParams[HURLEY_GBP::AHHe]   = HHE_RATE_CONSTANT_HURLEY2000;
    p_GBParams[HURLEY_GBP::AHe]    = HE_RATE_CONSTANT_HURLEY2000;

    p_GBParams[HURLEY_GBP::B]      = CalculateCoreMass_Luminosity_B_Hurley2000(p_Mass);
    p_GBParams[HURLEY_GBP::D]      = CalculateCoreMass_Luminosity_D_Hurley2000(p_Mass);

    p_GBParams[HURLEY_GBP::p]      = CalculateCoreMass_Luminosity_p_Hurley2000(p_Mass);
    p_GBParams[HURLEY_GBP::q]      = CalculateCoreMass_Luminosity_q_Hurley2000(p_Mass);

    p_GBParams[HURLEY_GBP::Mx]     = CalculateCoreMass_Luminosity_Mx_Hurley2000(p_GBParams);
    p_GBParams[HURLEY_GBP::Lx]     = CalculateCoreMass_Luminosity_Lx_Hurley2000(p_GBParams);

    p_GBParams[HURLEY_GBP::McBAGB] = CalculateCoreMassAtBAGB_Hurley2000(p_Mass);
    p_GBParams[HURLEY_GBP::McDU]   = CalculateCoreMassAt2ndDredgeUp_Hurley2000(p_GBParams[HURLEY_GBP::McBAGB]);
    p_GBParams[HURLEY_GBP::McBGB]  = CalculateCoreMassAtBGB_Hurley2000(p_Mass, p_GBParams[HURLEY_GBP::McBAGB]);

    return p_GBParams;
}


/*
 * Calculate Giant Branch (GB) parameters per Hurley et al. 2000
 *
 * Giant Branch Parameters depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 *
 * This function exists to facilitate the calculation of gbParams in EAGB::ResolveEnvelopeLoss() for the
 * stellar type to which the star will evolve.  The calculations of some of the stellar attributes there
 * depend on new gbParams.  
 * JR: This really needs to be revisited one day - these calculations should really be performed after
 *     switching to the new stellar type, but other calculations are done (in the legacy code) before the 
 *     switch (see evolveOneTimestep() in star.cpp for EAGB stars in the legacy code).
 *
 *
 * void CalculateGBParams_Static(const double      p_Mass, 
 *                               const double      p_LogMetallicityXi, 
 *                               const DBL_VECTOR &p_MassCutoffs, 
 *                               const DBL_VECTOR &p_An, 
 *                               const DBL_VECTOR &p_Bn, 
 *                                     DBL_VECTOR &p_GBParams) const
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN]        p_LogMetallicityXi      log10(Metallicity / Zsol) - called xi in Hurley et al. 2000
 * @param   [IN]        p_MassCutoffs           Mass cutoffs
 * @param   [IN]        p_An                    a(n) coefficients
 * @param   [IN]        p_Bn                    b(n) coefficients
 * @param   [IN/OUT]    p_GBParams              Giant Branch Parameters - calculated here
 */
void GiantBranch::CalculateGBParams_Static(const double      p_Mass, 
                                           const double      p_LogMetallicityXi, 
                                           const DBL_VECTOR &p_MassCutoffs, 
                                           const DBL_VECTOR &p_An, 
                                           const DBL_VECTOR &p_Bn, 
                                                 DBL_VECTOR &p_GBParams) const {
#define p_GBParams(x) p_GBParams[static_cast<int>(HURLEY_GB_PARAMETERS::x)]

    p_GBParams(AH)     = CalculateHRateConstant_Hurley2000(p_Mass);
    p_GBParams(AHHe)   = HHE_RATE_CONSTANT_HURLEY2000;
    p_GBParams(AHe)    = HE_RATE_CONSTANT_HURLEY2000;

    p_GBParams(B)      = CalculateCoreMass_Luminosity_B_Hurley2000(p_Mass);
    p_GBParams(D)      = CalculateCoreMass_Luminosity_D_Hurley2000(p_Mass);

    p_GBParams(p)      = CalculateCoreMass_Luminosity_p_Hurley2000(p_Mass);
    p_GBParams(q)      = CalculateCoreMass_Luminosity_q_Hurley2000(p_Mass);

    p_GBParams(Mx)     = CalculateCoreMass_Luminosity_Mx_Hurley2000(p_GBParams);
    p_GBParams(Lx)     = CalculateCoreMass_Luminosity_Lx_Hurley2000(p_GBParams);

    p_GBParams(McDU)   = CalculateCoreMassAt2ndDredgeUp_Hurley2000(p_GBParams(McBAGB));
    p_GBParams(McBAGB) = CalculateCoreMassAtBAGB_Hurley2000(p_Mass, p_Bn);
    p_GBParams(McBGB)  = CalculateCoreMassAtBGB_Hurley2000(p_Mass, p_GBParams(McBAGB));

#undef p_GBParams
}












/// JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

/*
 * Perturb Luminosity and Radius
 *
 * Perturbs Luminosity and Radius per Hurley et al. 2000, section 6.3
 * The attributes of the star are updated.
 *
 * Perturbation is disabled by default when DEBUG is enabled - except when
 * DEBUG_PERTURB is defined (see below).  The stellar class perturbation
 * functions are defined away if DEBUG is defined - so the generic Star
 * function is called (and does nothing). (So far FGB is the only class that
 * defines this function where it actually does anything)
 *
 * If DEBUG_PERTURB is defined then perturbation is not disabled while debugging.
 * To enable perturbation while DEBUG is enabled, define DEBUG_PERTURB.
 *
 *
 * void PerturbLuminosityAndRadius()
 */
#if !defined(DEBUG) || defined(DEBUG_PERTURB)   // don't perturb if debugging and DEBUG_PERTURB not defined
void GiantBranch::PerturbLuminosityAndRadius() {

    if (utils::Compare(CurrentState().Mu(), 1.0) < 0) {   // perturb only if mu < 1.0

        double Lc = CalculateRemnantLuminosity();
        double Rc = CalculateRemnantRadius();

        double s = CalculateHurleyPerturbationS(Mass(), CurrentState().Mu());
        double r = CalculateHurleyPerturbationR(Mass(), Radius(), CurrentState().Mu(), Rc);

        m_InterimState.SetLuminosity(Lc * PPOW((Luminosity() / Lc), s));        
        m_InterimState.SetRadius(Rc * PPOW((Radius() / Rc), r));
    }
}
#else
void GiantBranch::PerturbLuminosityAndRadius() { }
#endif






///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////













///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////







///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              SUPERNOVA CALCULATIONS                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



/*
 * CalculateRemnantMass_Schneider2020
 *
 * @brief
 * Calculate remnant mass for isolated single stars, per Schneider et al. 2020.
 * (See arxiv:2008.08599).
 * 
 * Function provides for standard and alternative fits (single stars) from Schneider et al. 2020.
 * 
 *
 * Dbl_DblT CalculateRemnantMass_Schneider2020(const double p_Mass,
 *                                             const double p_COCoreMass,
 *                                             const double p_HeCoreMassPreSN,
 *                                             const double p_MaxNSMass,
 *                                             const REMNANT_MASS_PRESCRIPTION p_RemnantMassPrescription) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMassPreSN               Pre-SN (at CO formation) He core mass (Msol)
 * @param       p_MaxNSMass                     User-specified maximum mass of a neutron star (`--maximum-neutron-star-mass`) (Msol)
 * @param       p_RemnantMassPrescription       The remnant mass prescription to use (should be SCHNEIDER2020 or SCHNEIDER2020ALT) 
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
Dbl_DblT GiantBranch::CalculateRemnantMass_Schneider2020(const double                    p_Mass,
                                                         const double                    p_COCoreMass,
                                                         const double                    p_HeCoreMassPreSN,
                                                         const double                    p_MaxNSMass,
                                                         const REMNANT_MASS_PRESCRIPTION p_RemnantMassPrescription) const {
    double logRemnantMass;
    switch (p_RemnantMassPrescription) {                                // which REMNANT_MASS_PRESCRIPTION?

        case REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020:

                 if (p_COCoreMass <  6.357) { logRemnantMass = std::log10(0.03357 * p_COCoreMass + 1.31780); }
            else if (p_COCoreMass <  7.311) { logRemnantMass = -0.02466 * p_COCoreMass + 1.28070; }
            else if (p_COCoreMass < 12.925) { logRemnantMass = std::log10(0.03357 * p_COCoreMass + 1.31780); }
            else                            { logRemnantMass = 0.01940 * p_COCoreMass + 0.98462; }
            break;

        case REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020ALT:

                 if (p_COCoreMass <  6.357) { logRemnantMass = std::log10(0.04199 * p_COCoreMass + 1.28128); }
            else if (p_COCoreMass <  7.311) { logRemnantMass = -0.02466 * p_COCoreMass + 1.28070; }
            else if (p_COCoreMass < 12.925) { logRemnantMass = std::log10(0.04701 * (p_COCoreMass * p_COCoreMass) - 0.91403 * p_COCoreMass + 5.93380); }
            else                            { logRemnantMass = 0.01940*p_COCoreMass + 0.98462; }
            break;

        default:                                                        // unexpected prescription
            // the only way this can happen is if the REMNANT_MASS_PRESCRIPTION passed to this function
            // is not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.

            THROW_ERROR(ERROR::UNEXPECTED_REMNANT_MASS_PRESCRIPTION);   // throw error  
    }

    // convert remnant mass to linear value, and limit to the pre-SN He Core mass
    const double remnantMass = std::min(PPOW(10.0, logRemnantMass), p_HeCoreMassPreSN);

    // fallback fraction of mass beyond proto-neutron-star for BH formation and kicks
    double fallbackFraction = remnantMass > p_MaxNSMass ? (remnantMass - NEUTRON_STAR_MASS) / (p_Mass - NEUTRON_STAR_MASS) : 0.0;

    return std::make_tuple(remnantMass, std::clamp(fallbackFraction, 0.0, 1.0));
}


/*
 * CalculateRemnantMass_Maltsev2025
 *
 * @brief
 * Calculate remnant mass for isolated single stars, per Maltsev et al. 2025.
 * (See https://doi.org/10.1051/0004-6361/202554931)
 * 
 *
 * Dbl_DblT CalculateRemnantMass_Maltsev2025(const double       p_COCoreMass,
 *                                         const double       p_HeCoreMass,
 *                                         const double       p_ZetaAsplund, 
 *                                         const double       p_MaltsevFallbackFraction,
 *                                         const MALTSEV_MODE p_MaltsevMode) const
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @param       p_ZetaAsplund                   Asplund zeta value (log10(Z / ZSOL_ASPLUND))
 * @param       p_MaltsevFallbackFraction       User-specfifid Maltsev fallback fraction (`--maltsev-fallback-fraction`) [0.0, 1.0]
 * @param       p_MaltsevMode                   User-specified Maltsev mode (`--maltsev-mode`))
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
Dbl_DblT GiantBranch::CalculateRemnantMass_Maltsev2025(const double p_COCoreMass, const double p_HeCoreMass) const {
    double fallbackFraction;
    double remnantMass;

    if (p_COCoreMass < MALTSEV2024_MMIN) {                                                      // CO core mass < Maltsev lower threshold?
                                                                                                // Yes - NS formation regardless of metallicity
        fallbackFraction = 0.0;
        remnantMass      = NEUTRON_STAR_MASS;
    }
    else if (p_COCoreMass > MALTSEV2024_MMAX) {                                                 // No - CO core mass > Maltsev upper threshold?
                                                                                                // Yes - BH formation regardless of metallicity
        fallbackFraction = 1.0;
        remnantMass      = p_HeCoreMass;
    }
    else {                                                                                      // No - calculate remnant mass

        double Zbound;                                                                      
        switch (OPTIONS->MaltsevMode()) {                                                       // Which Maltsev mode?                                                                                 

            case MALTSEV_MODE::OPTIMISTIC:                                                      // OPTIMISTIC               
                Zbound = ZDEP->ZetaAsplund(Metallicity());                                                // Asplund zeta (log10(Z / ZSOL_ASPLUND))
                break;

            case MALTSEV_MODE::BALANCED:                                                        // BALANCED                        
                Zbound = std::min(std::max(ZDEP->ZetaAsplund(Metallicity()), -1.6989700043360185), 0.0);  // Asplund zeta clamped to [log10(1/50), log10(1)]
                break;

            case MALTSEV_MODE::PESSIMISTIC:                                                     // Pessimistic                          
                Zbound = std::min(std::max(ZDEP->ZetaAsplund(Metallicity()), -1.0), 0.0);                 // Asplund zeta clamped to [log10(1/10), log10(1)]
                break;

            default:                                                                            // Unexpected mode
                // The only way this can happen is if the MALTSEV_MODE passed to this function
                // isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the code passed a mode that this function doesn't account for, and
                // that should be flagged as an error and result in termination of the evolution of the star
                // or binary.
                // The correct fix for this is to add code to this function for the missing mode, or fix the
                // calling code to pass a mode that is handled by this function.
                THROW_ERROR(ERROR::UNEXPECTED_MALTSEV_MODE);                                    // Throw error
        }
  
        const double M1 = MALTSEV2024_M1S + (MALTSEV2024_M1S - MALTSEV2024_M1SZ01) * Zbound;
        const double M2 = MALTSEV2024_M2S + (MALTSEV2024_M2S - MALTSEV2024_M2SZ01) * Zbound;
        const double M3 = MALTSEV2024_M3S + (MALTSEV2024_M3S - MALTSEV2024_M3SZ01) * Zbound;

        if (p_COCoreMass >= M3 || (p_COCoreMass >= M1 && p_COCoreMass <= M2) ) {                // Complete fallback into BH
            fallbackFraction = 1.0;
            remnantMass      = p_HeCoreMass;
        }
        else if (p_COCoreMass > M2 && p_COCoreMass < M3 && RAND->Random(0, 1) <= 0.1) {         // Partial fallback - BH formation
            fallbackFraction = OPTIONS->MaltsevFallbackFraction();
            remnantMass      = (p_HeCoreMass - NEUTRON_STAR_MASS) * fallbackFraction + NEUTRON_STAR_MASS;
        }
        else {                                                                                  // No fallback - NS formation
            fallbackFraction = 0.0;
            remnantMass      = NEUTRON_STAR_MASS;
        }
    }

    return std::make_tuple(remnantMass, fallbackFraction);  // JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
}



/*
 * CalculateRemnantMass_MullerMandel2020
 *
 * @brief
 * Calculate remnant mass, per Mandel & Mueller 2020
 *
 *
 * Dbl_DblT CalculateRemnantMass_MullerMandel2020 (const double p_COCoreMass, const double p_HeCoreMass)
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
Dbl_DblT GiantBranch::CalculateRemnantMass_MullerMandel2020(const double p_COCoreMass, const double p_HeCoreMass) {

    double remnantMass       = 0.0;   
    double pBH               = 0.0;
    double pCompleteCollapse = 0.0;
   
    if (utils::Compare(p_COCoreMass, MULLERMANDEL_M1) < 0 || utils::Compare(p_HeCoreMass, OPTIONS->MaximumNeutronStarMass()) <= 0 )
	    pBH = 0.0;
    else if (utils::Compare(p_COCoreMass, MULLERMANDEL_M3) < 0)
    	pBH = 1.0 / (MULLERMANDEL_M3-MULLERMANDEL_M1) * (p_COCoreMass-MULLERMANDEL_M1);
    else
	    pBH = 1.0;

    if (utils::Compare(RAND->Random(0, 1), pBH) < 0) {  // this is a BH
        if (utils::Compare(p_COCoreMass, MULLERMANDEL_M4) < 0)
		    pCompleteCollapse = 1.0 / (MULLERMANDEL_M4 - MULLERMANDEL_M1) * (p_COCoreMass - MULLERMANDEL_M1);
        else
		    pCompleteCollapse = 1.0;

	    if (utils::Compare(RAND->Random(0, 1), pCompleteCollapse) < 0)
		    remnantMass = p_HeCoreMass;
	    else
            remnantMass = CalculateBHMassAfterFallback_MullerMandel2020(p_COCoreMass, p_HeCoreMass);
    }
    else                                               // this is an NS
        remnantMass = CalculateRemnantNSMassMullerMandel(p_COCoreMass, p_HeCoreMass);

    // fallback fraction = 0.0 - no subsequent kick adjustment by fallback fraction needed
    return std::make_tuple(remnantMass, 0.0);
}


/*
 * Calculate NS remnant mass given COCoreMass and HeCoreMass
 *
 * Mandel & Mueller, 2020; also used by Maltsev2025, so separated out into helper function
 *
 *
 * double CalculateRemnantNSMassByMullerMandel (const double p_COCoreMass, const double p_HeCoreMass)
 *
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @param   [IN]    p_HeCoreMass                HeCoreMass in Msol
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantNSMassMullerMandel(const double p_COCoreMass, const double p_HeCoreMass) {

    double remnantMass = 0.0;
    SizeT  iterations  = 0;

    if (utils::Compare(p_COCoreMass, MULLERMANDEL_M1) < 0) {
        while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS            &&
               (utils::Compare(remnantMass, MULLERMANDEL_MINNS) < 0                ||
                utils::Compare(remnantMass, OPTIONS->MaximumNeutronStarMass()) > 0 ||
                utils::Compare(remnantMass, p_COCoreMass) > 0)) {
            remnantMass = MULLERMANDEL_MU1 + RAND->RandomGaussian(MULLERMANDEL_SIGMA1);
        }
        if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS)
            // Failure to find a solution implies a narrow range; just pick a midpoint in this case
            remnantMass = (std::min(OPTIONS->MaximumNeutronStarMass(), p_COCoreMass) + MULLERMANDEL_MINNS) / 2.0;
    }
    else if (utils::Compare(p_COCoreMass, MULLERMANDEL_M2) < 0) {
        while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS            &&
               (utils::Compare(remnantMass, MULLERMANDEL_MINNS) < 0                ||
                utils::Compare(remnantMass, OPTIONS->MaximumNeutronStarMass()) > 0 ||
                utils::Compare(remnantMass, p_COCoreMass) > 0)) {
            remnantMass = MULLERMANDEL_MU2A + MULLERMANDEL_MU2B / (MULLERMANDEL_M2 - MULLERMANDEL_M1) * (p_COCoreMass - MULLERMANDEL_M1) + RAND->RandomGaussian(MULLERMANDEL_SIGMA2);
        }
        if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS)
            // Failure to find a solution implies a narrow range; just pick a midpoint in this case
            remnantMass = (std::min(OPTIONS->MaximumNeutronStarMass(), p_COCoreMass) + MULLERMANDEL_MINNS) / 2.0;
    }
    else {
        while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS            &&
               (utils::Compare(remnantMass, MULLERMANDEL_MINNS) < 0                ||
                utils::Compare(remnantMass, OPTIONS->MaximumNeutronStarMass()) > 0 ||
                utils::Compare(remnantMass, p_COCoreMass) > 0)) {
            remnantMass = MULLERMANDEL_MU3A + MULLERMANDEL_MU3B / (MULLERMANDEL_M3 - MULLERMANDEL_M2) * (p_COCoreMass - MULLERMANDEL_M2) + RAND->RandomGaussian(MULLERMANDEL_SIGMA3);
        }
        if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS)
            // Failure to find a solution implies a narrow range; just pick a midpoint in this case
            remnantMass = (std::min(OPTIONS->MaximumNeutronStarMass(), p_COCoreMass) + MULLERMANDEL_MINNS) / 2.0;
    }
    return remnantMass;
}




/*
 * CalculateGravitationalRemnantMass_Fryer2012
 *
 * @brief
 * Calculate the gravitational mass of the remnant, per
 *
 *    - Fryer et al. 2012, eq 13 for Neutron Stars
 *    - An enhanced version of Fryer et al. 2012, eq 14 for Black Holes
 *
 *
 * double CalculateGravitationalRemnantMass_Fryer2012(const double p_BaryonicRemnantMass)
 *
 * @param       p_BaryonicRemnantMass           Baryonic remnant mass (Msol)
 * @param       p_NSmaxBaryonicMass             Maximum baryonic mass for NS maximum mass (Msol)
 * @return                                      Gravitational mass of the remnant (Msol)
 */
double GiantBranch::CalculateGravitationalRemnantMass_Fryer2012(const double p_BaryonicRemnantMass, const double p_NSmaxBaryonicMass) const {

    double mass;
    
    if (p_BaryonicRemnantMass < p_NSmaxBaryonicMass) {                                      // below NS threshold?
        ERROR error;                                                                        // yes - neutron star
        std::tie(error, mass) = utils::SolveQuadratic(0.075, 1.0, -p_BaryonicRemnantMass);  // find gravitational mass
        if (error != ERROR::NONE) THROW_ERROR(error);                                       // throw error if no root found
    } 
    else {                                                                                  // no - black hole
        mass = BH::CalculateNeutrinoMassLoss_Static(p_BaryonicRemnantMass);                 // convert baryonic mass to gravitational mass
    }

    return mass;
}


/*
 * CalculateFallbackFraction_Fryer2012_Rapid
 *
 * @brief
 * Calculate the fraction of mass falling back onto the proto compact object,
 * using the rapid supernove mechanism, per Fryer et al. 2012, eq 15s & 16.
 *
 *
 * double CalculateFallbackFraction_Fryer2012_Rapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const
 *
 * @param       p_PreSNMass                     Pre-SN mass of the star (Msol)
 * @param       p_ProtoMass                     Fe/Ni proto compact object mass (Msol)
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
double GiantBranch::CalculateFallbackFraction_Fryer2012_Rapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const {

    double fallbackFraction;

    if (p_COCoreMass < 2.5) {
        fallbackFraction = 0.2 / (p_PreSNMass - p_ProtoMass);
    }
    else if (p_COCoreMass < 6.0) {
        fallbackFraction = ((0.286 * p_COCoreMass) - 0.514) / (p_PreSNMass - p_ProtoMass);
    }
    else if (p_COCoreMass < 7.0) {
        fallbackFraction = 1.0;
    }
    else if (p_COCoreMass < 11.0) {
        const double a1 = 0.25 - (1.275 / (p_PreSNMass - p_ProtoMass));
        fallbackFraction = (a1 * p_COCoreMass) + ((-11.0 * a1) + 1.0);
    }
    else {
        fallbackFraction = 1.0;
    }

    return std::clamp(fallbackFraction, 0,0, 1.0);
}



/*
 * CalculateFallbackFraction_Fryer2012_Delayed
 *
 * @brief
 * Calculate the fraction of mass falling back onto the proto compact object,
 * using the delayed supernova mechanis, per Fryer et al. 2012, eq 19.
 *
 *
 * double CalculateFallbackFraction_Fryer2012_Delayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const
 *
 * @param       p_PreSNMass                     Pre-SN mass of the star (Msol)
 * @param       p_ProtoMass                     Fe/Ni proto compact object mass (Msol)
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
double GiantBranch::CalculateFallbackFraction_Fryer2012_Delayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const {

    double fallbackFraction;

    if (p_COCoreMass < 2.5) {
        fallbackFraction = 0.2 / (p_PreSNMass - p_ProtoMass);
    }
    else if (p_COCoreMass < 3.5) {
        fallbackFraction = ((0.5 * p_COCoreMass) - 1.05) / (p_PreSNMass - p_ProtoMass);
    }
    else if (p_COCoreMass < 11.0) {
        const double a2 = 0.133 - (0.093 / (p_PreSNMass - p_ProtoMass));
        fallbackFraction = (a2 * p_COCoreMass) + ((-11.0 * a2) + 1.0);
    }
    else {
        fallbackFraction = 1.0;
    }

    return std::clamp(fallbackFraction, 0.0, 1.0);
}


/*
 * CalculateRemnantMass_Fryer2012
 *
 * @brief
 * Calculate remnant mass, per Fryer et al. 2012.
 *
 *
 * Dbl_DblT CalculateRemnantMass_Fryer2012(const double p_Mass, const double p_COCoreMass) const
 *
 * @param       p_Mass                          Pre-SN mass of the star (Msol)
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
Dbl_DblT GiantBranch::CalculateRemnantMass_Fryer2012(const double p_Mass, const double p_COCoreMass) const {

    double baryonicMass;
    double fallbackFraction;
    double gravitationalMass;
    double mProto;

    switch (OPTIONS->FryerSupernovaEngine()) {                              // which Fryer SN mechanism?

        case SN_ENGINE::DELAYED:                                            // DELAYED

            mProto            = CalculateProtoCoreMass_Fryer2012_Delayed(p_COCoreMass);
            fallbackFraction  = CalculateFallbackFraction_Fryer2012_Delayed(p_Mass, mProto, p_COCoreMass);
            baryonicMass      = CalculateBaryonicRemnantMass_Fryer2012(mProto, CalculateFallbackMass_Fryer2012(p_Mass, mProto, fallbackFraction));
            gravitationalMass = CalculateGravitationalRemnantMass_Fryer2012(baryonicMass, GLOBALS->NSMaxBaryonicMass());
            break;

        case SN_ENGINE::RAPID:                                              // RAPID

            mProto            = FRYER_PROTO_CORE_MASS_RAPID;
            fallbackFraction  = CalculateFallbackFractionRapid(p_Mass, mProto, p_COCoreMass);
            baryonicMass      = CalculateBaryonicRemnantMass_Fryer2012(mProto, CalculateFallbackMass_Fryer2012(p_Mass, mProto, fallbackFraction));
            gravitationalMass = CalculateGravitationalRemnantMass_Fryer2012(baryonicMass, GLOBALS->NSMaxBaryonicMass());
            break;
    
        default:                                                            // unexpected prescription
            // the only way this can happen is if someone added an SN_ENGINE and it isn't accounted for
            // in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.
            THROW_ERROR(ERROR::UNEXPECTED_SN_ENGINE);                       // throw error
    }

    return std::make_tuple(gravitationalMass, std::clamp(fallbackFraction, 0.0, 1.0));
}


/*
 * CalculateRemnantMass_Fryer2022
 *
 * @brief
 * Calculate the remnant mass, per Fryer et al. 2022, eq 5.
 *
 *
 * Dbl_DblT CalculateRemnantMass_Fryer2022(const double p_Mass, const double p_COCoreMass, const double p_fMix, const double p_mCrit) const
 *
 * @param       p_Mass                          Pre-SN mass of the star (Msol)
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @param       p_fMix                          Fryer 2022 mixing growth time
 * @param       p_mCrit                         Fryer 2022 critical mass for BH formation (Msol)        
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
Dbl_DblT GiantBranch::CalculateRemnantMass_Fryer2022(const double p_Mass, const double p_COCoreMass, const double p_fMix, const double p_mCrit) const {

    // calculate baryonic mass, clamped to total mass
    const double f = p_COCoreMass / p_fMix;
    double baryonicMass = std::min(1.2 + 0.05 * p_fMix + 0.01 * f * f + std::exp(p_fMix * (p_COCoreMass - p_mCrit)), p_Mass);

    // calculate fallback fraction and gravitational mass
    double mProto;
    [[maybe_unused]] double fallbackMass;
    double fallbackFraction;
    double gravitationalMass;

    switch (OPTIONS->FryerSupernovaEngine()) {                              // which Fryer SN mechanism?

        case SN_ENGINE::DELAYED:                                            // DELAYED  

            mProto            = CalculateProtoCoreMass_Fryer2012_Delayed(p_COCoreMass);
            fallbackFraction  = std::clamp(std::max(0.0, baryonicMass - mProto) / (p_Mass - mProto), 0.0, 1.0);
            gravitationalMass = CalculateGravitationalRemnantMass_Fryer2012(baryonicMass, m_FixedAttributes->BaryonicMassOfMaxMassNS()); // JR FIX THIS - don't use class member <<<<<<<<<<<<<<<<<<<<<<<<<<<
            break;

        case SN_ENGINE::RAPID:  

            mProto            = FRYER_PROTO_CORE_MASS_RAPID;
            fallbackFraction  = std::clamp(std::max(0.0, baryonicMass - mProto) / (p_Mass - mProto), 0.0, 1.0);
            gravitationalMass = CalculateGravitationalRemnantMass_Fryer2012(baryonicMass, m_FixedAttributes->BaryonicMassOfMaxMassNS()); // JR FIX THIS - don't use class member <<<<<<<<<<<<<<<<<<<<<<<<<<<
            break;
    
        default:                                                            // unexpected prescription
            // the only way this can happen is if someone added an SN_ENGINE and it isn't accounted for
            // in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.
            THROW_ERROR(ERROR::UNEXPECTED_SN_ENGINE);                       // throw error
    }
                                   
    return std::make_tuple(gravitationalMass, std::clamp(fallbackFraction, 0.0, 1.0));
}





/*
 * Driver function for Core Collapse Supernovas
 *
 * This function determines which prescription is used for the core collapse SN (via program options)
 *
 * The function calls prescription functions that update the following parameters:
 *      Mass, stellarType, drawnKickMagnitude, kickMagnitude
 *
 * At the end of this function we set the following parameters which are (so far) independent of the
 * ccSN prescriptions (but do depend on the parameters above):
 *      Luminosity, Radius, Temperature, supernova events: current = SN, past = CCSN
 *
 *
 * STELLAR_TYPE ResolveCoreCollapseSN()
 *
 * @return                                      The stellar type to which the star should evolve
 */
STELLAR_TYPE GiantBranch::ResolveCoreCollapseSN() {

    STELLAR_TYPE stellarType = m_StellarType;
    double mass              = m_Mass;                                                                      // initial mass

    switch (OPTIONS->RemnantMassPrescription()) {                                                           // which prescription?

        case REMNANT_MASS_PRESCRIPTION::HURLEY2000:                                                         // Hurley 2000

            m_SupernovaDetails.fallbackFraction = 0.0;                                                      // not defined
            m_Mass                              = Remnants::CalculateRemnantMass_Static(m_COCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::BELCZYNSKI2002:                                                     // Belczynski 2002

            m_SupernovaDetails.fallbackFraction = CalculateFallbackByBelczynski2002(m_CoreMass);
            m_Mass                              = CalculateRemnantMassByBelczynski2002(m_Mass, m_CoreMass, m_SupernovaDetails.fallbackFraction);
            break;

        case REMNANT_MASS_PRESCRIPTION::FRYER2012:                                                          // Fryer 2012

            std::tie(m_Mass, m_SupernovaDetails.fallbackFraction) = CalculateRemnantMassByFryer2012(m_Mass, m_COCoreMass);
            break;


        case REMNANT_MASS_PRESCRIPTION::FRYER2022:                                                          // Fryer 2022

            std::tie(m_Mass, m_SupernovaDetails.fallbackFraction) = CalculateRemnantMassByFryer2022(m_Mass, m_COCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::MULLER2016:                                                         // Muller 2016

            m_SupernovaDetails.fallbackFraction = 0.0;                                                      // no subsequent kick adjustment by fallback fraction needed
            m_Mass                              = CalculateRemnantMassByMuller2016(m_Mass, m_COCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::MULLERMANDEL:                                                       // Mandel & Mueller, 2020

            m_SupernovaDetails.fallbackFraction = 0.0;                                                      // no subsequent kick adjustment by fallback fraction needed
            m_Mass                              = CalculateRemnantMassByMullerMandel(m_COCoreMass, m_HeCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020:                                                      // Schneider 2020

            m_Mass                              = CalculateRemnantMassBySchneider2020(m_COCoreMass);
            m_SupernovaDetails.fallbackFraction = utils::Compare(m_Mass, OPTIONS->MaximumNeutronStarMass() ) > 0 ? (m_Mass - NEUTRON_STAR_MASS) / (mass - NEUTRON_STAR_MASS) : 0.0;		// fallback fraction of mass beyond proto-neutron-star for BH formation and kicks
            break;

        case REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020ALT:                                                   // Schneider 2020, alternative

            m_Mass                              = CalculateRemnantMassBySchneider2020Alt(m_COCoreMass);
            m_SupernovaDetails.fallbackFraction = utils::Compare(m_Mass, OPTIONS->MaximumNeutronStarMass() ) > 0 ? (m_Mass - NEUTRON_STAR_MASS) / (mass - NEUTRON_STAR_MASS) : 0.0;		// fallback fraction of mass beyond proto-neutron-star for BH formation and kicks
            break;
        
        case REMNANT_MASS_PRESCRIPTION::MALTSEV2024:                                                        // Maltsev+ 2024

            m_Mass                              = CalculateRemnantMassByMaltsev2024(m_COCoreMass, m_HeCoreMass);																		// fallback fraction determined internally
            break;
            
    
        default:                                                                                            // unknown prescription
            // the only way this can happen is if someone added a REMNANT_MASS_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_REMNANT_MASS_PRESCRIPTION);                                          // throw error
    }
    
    // Set the stellar type to which the star should evolve (either use prescription or MAXIMUM_NS_MSS)
    if (OPTIONS->RemnantMassPrescription() == REMNANT_MASS_PRESCRIPTION::MULLER2016) {
        stellarType = CalculateRemnantTypeByMuller2016(m_COCoreMass);
    }
    else if (OPTIONS->RemnantMassPrescription() == REMNANT_MASS_PRESCRIPTION::MULLERMANDEL) {
        if (utils::Compare(m_Mass, OPTIONS->MaximumNeutronStarMass() ) > 0)
            stellarType = STELLAR_TYPE::BLACK_HOLE;
        else
            stellarType = STELLAR_TYPE::NEUTRON_STAR;
    }
    else if (OPTIONS->RemnantMassPrescription() == REMNANT_MASS_PRESCRIPTION::HURLEY2000) {
        stellarType = (utils::Compare(m_Mass, 1.8 ) > 0)
                        ? STELLAR_TYPE::BLACK_HOLE
                        : STELLAR_TYPE::NEUTRON_STAR;                                                       // Hurley+ 2000, Eq. (92)
    }
    else if (utils::Compare(m_Mass, OPTIONS->MaximumNeutronStarMass()) > 0) {
        std::tie(m_Luminosity, m_Radius, m_Temperature) = BH::CalculateCoreCollapseSNParams_Static(m_Mass);
        stellarType = STELLAR_TYPE::BLACK_HOLE;
    }
    else {
        std::tie(m_Luminosity, m_Radius, m_Temperature) = NS::CalculateCoreCollapseSNParams_Static(m_Mass);
        stellarType = STELLAR_TYPE::NEUTRON_STAR;
    }

    if (utils::Compare(mass,m_CoreMass) == 0 && utils::Compare(m_HeCoreMass, m_COCoreMass) == 0) {          // entire star is CO core, so this is a USSN
        SetSNCurrentEvent(SN_EVENT::USSN);                                                                  // flag ultra-stripped SN happening now
        SetSNPastEvent(SN_EVENT::USSN);                                                                     // ... and will be a past event
    }

    SetSNCurrentEvent(SN_EVENT::CCSN);                                                                      // flag core-collapse SN happening now
    SetSNPastEvent(SN_EVENT::CCSN);                                                                         // ... and will be a past event

    return stellarType;
}


/*
 * Resolve Electron Capture Supernova
 *
 * Calculate the mass of the remnant and set remnant type - always a Neutron Star
 * Updates attributes of star; sets SN flags
 *
 *
 * Short hand wavy story. The core ignites ONeMg and collapses through electron
 * capture (e.g., Nomoto 1984 for thorough discussion and a summary in Nomoto 1987).
 * The explosion is likely accompanied by a low natal kick
 *
 *
 * STELLAR_TYPE ResolveElectronCaptureSN()
 *
 * @return                                      Stellar type of remnant (always STELLAR_TYPE::NEUTRON_STAR)
 */
STELLAR_TYPE GiantBranch::ResolveElectronCaptureSN() {
    m_Mass       = MECS_REM;                                                                // defined in constants.h
    m_CoreMass   = m_Mass;
    m_HeCoreMass = m_Mass;
    m_COCoreMass = m_Mass;
    m_Mass0      = m_Mass;
        
    SetSNCurrentEvent(SN_EVENT::ECSN);                                                      // electron capture SN happening now
    SetSNPastEvent(SN_EVENT::ECSN);                                                         // ... and will be a past event

    return STELLAR_TYPE::NEUTRON_STAR;
}


/*
 * Resolve Pair-Instability Supernova
 *
 * Calculate the mass of the remnant and set remnant type according to mass
 * Updates attributes of star; sets SN events
 *
 *
 * Short handwavy story.  The core is hot and massive enough that there is a significant
 * amount of high-energetic gamma rays which can create an electron-positron pair.
 * When a significant amounts of pairs are created the photons taken away
 * reduce the radiative pressure. The core contracts becomes hotter creating
 * more pairs. This runaway process will end in a SN that explodes the entire core without
 * leaving a remnant.
 *
 *
 * STELLAR_TYPE ResolvePairInstabilitySN()
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch::ResolvePairInstabilitySN() {

    m_Luminosity  = 0.0;
    m_Radius      = 0.0;
    m_Temperature = 0.0;
    m_Mass        = 0.0;
    m_CoreMass    = 0.0;
    m_COCoreMass  = 0.0;
    m_HeCoreMass  = 0.0;

    m_SupernovaDetails.drawnKickMagnitude = 0.0;
    m_SupernovaDetails.kickMagnitude      = 0.0;
    m_SupernovaDetails.fallbackFraction   = 0.0;

    SetSNCurrentEvent(SN_EVENT::PISN);                                                                  // pair instability SN happening now
    SetSNPastEvent(SN_EVENT::PISN);                                                                     // ... and will be a past event

    return STELLAR_TYPE::MASSLESS_REMNANT;
}


/*
 * Resolve Pulsational Pair-Instability Supernova
 *
 * Calculate the mass of the remnant and set remnant type according to mass
 * Updates attributes of star; sets SN events
 *
 *
 * STELLAR_TYPE ResolvePulsationalPairInstabilitySN()
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch::ResolvePulsationalPairInstabilitySN() {

    STELLAR_TYPE stellarType = m_StellarType;

    double baryonicMass;
    switch (OPTIONS->PulsationalPairInstabilityPrescription()) {                                        // which prescription?

        case PPI_PRESCRIPTION::COMPAS:                                                                  // deprecated Feb 2025, to be removed
        case PPI_PRESCRIPTION::WOOSLEY:                                                                 // Woosley 2017 https://arxiv.org/abs/1608.08939
            baryonicMass = m_HeCoreMass;                                                                // strip off the hydrogen envelope if any was left
            m_Mass       = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                          // convert to gravitational mass due to neutrino mass loss
            break;

        case PPI_PRESCRIPTION::STARTRACK:                                                               // Belczynski et al. 2016 https://arxiv.org/abs/1607.03116
            baryonicMass = std::min(m_HeCoreMass, STARTRACK_PPISN_HE_CORE_MASS);                        // strip off the hydrogen envelope if any was left, limit to STARTRACK_PPISN_HE_CORE_MASS (default 45 Msun)
            m_Mass       = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                          // convert to gravitational mass due to neutrino mass loss
            break;

        case PPI_PRESCRIPTION::MARCHANT: {                                                              // Marchant et al. 2018 https://arxiv.org/abs/1810.13412

            // pow() is slow - use multiplication
            double HeCoreMass_2 = m_HeCoreMass * m_HeCoreMass;
            double HeCoreMass_3 = HeCoreMass_2 * m_HeCoreMass;
            double HeCoreMass_4 = HeCoreMass_2 * HeCoreMass_2;
            double HeCoreMass_5 = HeCoreMass_3 * HeCoreMass_2;
            double HeCoreMass_6 = HeCoreMass_3 * HeCoreMass_3;
            double HeCoreMass_7 = HeCoreMass_6 * m_HeCoreMass;

            double ratioOfRemnantToHeCoreMass = std::max(0.0, std::min(1.0, (-1.63057326E-08 * HeCoreMass_7) +
                                                                            ( 5.36316755E-06 * HeCoreMass_6) +
                                                                            (-7.52206933E-04 * HeCoreMass_5) +
                                                                            ( 5.83107626E-02 * HeCoreMass_4) +
                                                                            (-2.69801221E+00 * HeCoreMass_3) +
                                                                            ( 7.45060098E+01 * HeCoreMass_2) +
                                                                            (-1.13694590E+03 * m_HeCoreMass) +
                                                                              7.39643451E+03));

            baryonicMass = ratioOfRemnantToHeCoreMass * m_HeCoreMass;                                   // strip off the hydrogen envelope if any was left
            m_Mass       = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                          // convert to gravitational mass due to neutrino mass loss
            } break;

        case PPI_PRESCRIPTION::FARMER: {                                                                // Farmer et al. 2019 http://dx.doi.org/10.3847/1538-4357/ab518b
            double totalMassPrePPISN = m_Mass;                                                          // save the total stellar mass 
                                                                                                        // three cases:
            if (utils::Compare(m_COCoreMass, FARMER_PPISN_UPP_LIM_LIN_REGIME) < 0) {
                m_Mass = m_COCoreMass + 4.0;                                                            // a linear relation below CO core masses of 38 Msun
            }
            else if (utils::Compare(m_COCoreMass, FARMER_PPISN_UPP_LIM_QUAD_REGIME) < 0) {              // a quadratic relation in CO core mass for 38 =< CO_core < 60  JR: shouldn't have constants here (in the comment) - they may change
                const double a1 = -0.096;
                const double a2 = 8.564;
                const double a3 = -2.07;
                const double a4 = -152.97;

                m_Mass = a1 * PPOW(m_COCoreMass, 2.0) + a2 * m_COCoreMass + a3 * m_Log10Metallicity + a4;
            }
            else if (utils::Compare(m_COCoreMass, FARMER_PPISN_UPP_LIM_INSTABILLITY) < 0) {             // no remnant between 60 - 140 Msun  JR: shouldn't have constants here (in the comment) - they may change
                m_Mass = 0.0;
            }
            else {                                                                                      // BH mass becomes CO-core mass above the PISN gap
                m_Mass = m_COCoreMass;
            }

            m_Mass = std::min(totalMassPrePPISN, m_Mass);                                               // check if remnant mass is bigger than total mass    
            } break;
    
        case PPI_PRESCRIPTION::HENDRIKS: {    
            // Prescription from Hendriks et al. 2023 (https://arxiv.org/abs/2309.09339)
            // Based on Renzo et al. 2022 (https://iopscience.iop.org/article/10.3847/2515-5172/ac503e)
            // 
            // Suggest using --PPI-upper-limit 80.0 and --PISN-lower-limit 80.0

            double DeltaMPPICOShift = OPTIONS->PulsationalPairInstabilityCOCoreShiftHendriks();
            double DeltaMPPIExtraML = 0.0; 								// Make an option? Currently does nothing

            // Equation (6) of Hendricks et al. 2023			
            double PPIOnset        = m_COCoreMass - DeltaMPPICOShift - 34.8;
            double PPIOnsetSquared = PPIOnset * PPIOnset;
            double PPIOnsetCubed   = PPIOnsetSquared * PPIOnset;
            double firstTerm  = (0.0006 * m_Log10Metallicity + 0.0054) * PPIOnsetCubed;
            double secondTerm = 0.0013 * PPIOnsetSquared;
            double DeltaMPPI  = firstTerm - secondTerm + DeltaMPPIExtraML;
            
            DeltaMPPI = std::max(DeltaMPPI, 0.0);						// DeltaMPPI, the amount of the He core that's lost in pulsations, is non-negative
            m_Mass = std::max(m_HeCoreMass - DeltaMPPI, 0.0);			// Remnant mass should be non-negative		
            m_Mass = m_Mass > 10.0 ? m_Mass : 0.0;                      // If the predicted remnant mass is below 10 Msun, set it equal to 0 (assume a PISN)

        } break;

        default:                                                                                        // unknown prescription
            // the only way this can happen is if someone added a REMNANT_MASS_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_PPI_PRESCRIPTION);                                               // throw error
    }

    if (utils::Compare(m_Mass, 0.0) <= 0) {                                                             // remnant mass <= 0?
        stellarType = ResolvePairInstabilitySN();                                                       // yes - PISN rather than PPISN
    }
    else {                                                                                              // no - PPISN
        SetSNCurrentEvent(SN_EVENT::PPISN);                                                             // pulsational pair instability SN happening now
        SetSNPastEvent(SN_EVENT::PPISN);                                                                // ... and will be a past event

        stellarType   = STELLAR_TYPE::BLACK_HOLE;                                                       // -> black hole
        
        m_Luminosity  = BH::CalculateLuminosityOnPhase_Static();                                        // black hole luminosity
        m_Radius      = BH::CalculateRadiusOnPhase_Static(m_Mass);                                      // Schwarzschild radius (not correct for rotating BH)
        m_Temperature = CalculateTemperatureOnPhase(m_Luminosity, m_Radius);

        m_SupernovaDetails.fallbackFraction = 1.0;                                                      // fraction of mass that falls back
    }

    return stellarType;
}


/*
 * The main supernova function
 *
 * This function determines the type of the supernova and calls the appropriate functions
 * to calculate attributes correctly, and to determine the type of remnant to which the
 * star should evolve.
 *
 *
 * STELLAR_TYPE ResolveSupernova()
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch::ResolveSupernova() {

    STELLAR_TYPE stellarType = m_StellarType;

    if (IsSupernova()) {                                                                            // is going supernova
                                                                                                    // yes - resolve new supernova event
        // squirrel away some attributes before they get changed...
        m_SupernovaDetails.totalMassAtCOFormation  = m_Mass;
        m_SupernovaDetails.totalRadiusAtCOFormation= m_Radius;
        m_SupernovaDetails.HeCoreMassAtCOFormation = m_HeCoreMass;
        m_SupernovaDetails.COCoreMassAtCOFormation = m_COCoreMass;
        m_SupernovaDetails.coreMassAtCOFormation   = m_CoreMass;
        m_SupernovaDetails.coreRadiusAtCOFormation = CalculateConvectiveCoreRadius();

        double snMass = CalculateInitialSupernovaMass();                                            // calculate SN initial mass
        
        SetSNHydrogenContent();                                                                     // set H-rich or H-poor  JR: why don't we do this when we initialise the star at change of stellar type?

        if (                             OPTIONS->UsePulsationalPairInstability()              &&
            utils::Compare(m_HeCoreMass, OPTIONS->PulsationalPairInstabilityLowerLimit()) >= 0 &&
            utils::Compare(m_HeCoreMass, OPTIONS->PulsationalPairInstabilityUpperLimit()) <= 0) {   // Pulsational Pair Instability Supernova

            stellarType = ResolvePulsationalPairInstabilitySN();                                    // BH or MR
        }
        else if (                        OPTIONS->UsePairInstabilitySupernovae()    &&
            utils::Compare(m_HeCoreMass, OPTIONS->PairInstabilityLowerLimit()) >= 0 &&
            utils::Compare(m_HeCoreMass, OPTIONS->PairInstabilityUpperLimit()) <= 0) {              // Pair Instability Supernova

            stellarType = ResolvePairInstabilitySN();                                               // MR
        }
        else if (utils::Compare(snMass, MCBUR2) < 0 && (!m_MassTransferDonorHistory.empty() || OPTIONS->AllowNonStrippedECSN())) {
            stellarType = ResolveElectronCaptureSN();                                               // electron capture SN; requires progenitor to have been a MT donor unless non-stripped ECSN are allowed; forms NS
        }
        else {                                                                                      // Core Collapse Supernova
            stellarType = ResolveCoreCollapseSN();                                                  // BH or NS
        }
        
        // check if the SN actually happened
        if (utils::IsOneOf(stellarType, { STELLAR_TYPE::NEUTRON_STAR, STELLAR_TYPE::BLACK_HOLE, STELLAR_TYPE::MASSLESS_REMNANT })) {
                                                                                                    // SN happened
            if (utils::SNEventType(m_SupernovaDetails.events.current) != SN_EVENT::PISN && !utils::IsOneOf(stellarType, { STELLAR_TYPE::MASSLESS_REMNANT }))
                CalculateSNKickMagnitude(m_Mass, m_SupernovaDetails.totalMassAtCOFormation - m_Mass, stellarType);
        
            if (!utils::IsOneOf(stellarType, { STELLAR_TYPE::NEUTRON_STAR }))
                m_SupernovaDetails.rocketKickMagnitude = 0;                                         // only NSs can get rocket kicks

            // Print SN details to the SSE Supernova log.
            // Only if SSE (BSE does its own SN printing), and only if not an ephemeral clone
            if (OPTIONS->EvolutionMode() == EVOLUTION_MODE::SSE && m_ObjectPersistence == OBJECT_PERSISTENCE::PERMANENT) {
                PrintSupernovaDetails();
            }
       }
    }

    return stellarType;
}














///////////////////////////////////////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<







/*
 * CalculateRemnantMass_Schneider2020
 *
 * @brief
 * Calculate remnant mass for binary stars, per Schneider et al. 2020.
 * (See arxiv:2008.08599).
 * 
 * Function provides for standard and alternative fits (single stars) from Schneider et al. 2020.
 *
 * Remnant mass is based on the donor type at the latest (most recent) MT event.
 * 
 * Note that Schneider only prescribes remnant masses for the simple cases of single episode mass
 * transfer, so some of the double episode cases here are a bit uncertain, and may need to be refined
 * at a later date.
 * 
 *
 * Dbl_DblT CalculateRemnantMass_Schneider2020(const double                    p_Mass,
 *                                             const double                    p_COCoreMass,
 *                                             const double                    p_HeCoreMassPreSN,
 *                                             const double                    p_MaxNSMass,
 *                                             const MT_CASE                   p_MTcase,
 *                                             const REMNANT_MASS_PRESCRIPTION p_RemnantMassPrescription) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMassPreSN               Pre-SN He core mass of the star (Msol)
 * @param       p_MaxNSMass                     User-specified maximum mass of a neutron star (`--maximum-neutron-star-mass`) (Msol)
 * @param       p_MTcase                        Mass transfer case of latest MT event (MT_CASE::NONE if no MT history) // m_StateHistory.CurrentState.LatestMTdonorEvent().MTcase <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 * @param       p_RemnantMassPrescription       Remnant mass prescription to use
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
Dbl_DblT GiantBranch_Constituent::CalculateRemnantMass_Schneider2020(const double                    p_Mass,
                                                                     const double                    p_COCoreMass,
                                                                     const double                    p_HeCoreMassPreSN,
                                                                     const double                    p_MaxNSMass,
                                                                     const MT_CASE                   p_MTcase,
                                                                     const REMNANT_MASS_PRESCRIPTION p_RemnantMassPrescription) const {
    double remnantMass;
    double fallbackFraction;

    // if no history of MT, or if MT happens from naked He stars, WDs, etc.,
    // assume that the core properties are not affected
    if (p_MTcase == MT_CASE::NONE || p_MTcase == MT_CASE::OTHER) {
        remnantMass = 0.0;
        fallbackFraction = 0.0;
    }
    else {                                                                          // handle MT cases
        double logRemnantMass = 0.0;
        switch (p_MTcase) {

            case MT_CASE::A:                                                        // case A MT

                     if (p_COCoreMass <  7.064) { logRemnantMass = std::log10(0.02128 * p_COCoreMass + 1.35349); }
                else if (p_COCoreMass <  8.615) { logRemnantMass = 0.03866 * p_COCoreMass + 0.64417; }
                else if (p_COCoreMass < 15.187) { logRemnantMass = std::log10(0.02128 * p_COCoreMass + 1.35349); }
                else                            { logRemnantMass = 0.02573 * p_COCoreMass + 0.79027; }
                break;

            case MT_CASE::B:                                                        // case B MT

                     if (p_COCoreMass <  7.548) { logRemnantMass = std::log10(0.01909 * p_COCoreMass + 1.34529); }
                else if (p_COCoreMass <  8.491) { logRemnantMass = 0.03306 * p_COCoreMass + 0.68978; }
                else if (p_COCoreMass < 15.144) { logRemnantMass = std::log10(0.01909 * p_COCoreMass + 1.34529); }
                else                            { logRemnantMass = 0.02477 * p_COCoreMass + 0.80614; }
                break;

            case MT_CASE::C:                                                        // case C MT

                     if (p_COCoreMass <  6.357) { logRemnantMass = std::log10(0.03781 * p_COCoreMass + 1.36363); }
                else if (p_COCoreMass <  7.311) { logRemnantMass = 0.05264 * p_COCoreMass + 0.58531; }
                else if (p_COCoreMass < 14.008) { logRemnantMass = std::log10(0.03781 * p_COCoreMass + 1.36363); }
                else                            { logRemnantMass = 0.01795 * p_COCoreMass + 0.98797; }
                break;

            default:                                                              // unexpected MT_CASE
                // the only way this can happen is if the MT_CASE passed to this function is not accounted
                // for in this code.  We should not default here, with or without a warning.
                // We are here because the code passed a mass transfer case that this function doesn't account
                // for, and that should be flagged as an error and result in termination of the evolution of
                // binary.
                // The correct fix for this is to add code to this function for the missing mass transfer case,
                // or fix the calling code to pass a mass transfer case that is handled by this function.

                THROW_ERROR(ERROR::UNEXPECTED_MT_CASE);                             // throw error
        }

        // convert log mass to linear and clamp to pre-SN He core mass of the star
        remnantMass = std::min(PPOW(10.0, logRemnantMass), p_HeCoreMassPreSN);
    }

    // fallback fraction of mass beyond proto-neutron-star for BH formation and kicks
    fallbackFraction = remnantMass > p_MaxNSMass ? (remnantMass - NEUTRON_STAR_MASS) / (p_Mass - NEUTRON_STAR_MASS) : 0.0;

    return std::make_tuple(remnantMass, std::clamp(fallbackFraction, 0.0, 1.0));
}



/*
 * CalculateRemnantMass_Maltsev2025
 *
 * Calculate remnant mass for binary stars, per Maltsev et al. 2025.
 * (See https://doi.org/10.1051/0004-6361/202554931)
 *
 * Remnant mass is based on the donor type at the first MT event.
 * 
 *
 * Dbl_DblT CalculateRemnantMass_Maltsev2025(const double       p_COCoreMass,
 *                                           const double       p_HeCoreMass,
 *                                           const double       p_ZetaAsplund, 
 *                                           const double       p_MaltsevFallbackFraction,
 *                                           const MT_CASE      p_MTcase,
 *                                           const MALTSEV_MODE p_MaltsevMode) const
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @param       p_ZetaAsplund                   Asplund zeta value (log10(Z / ZSOL_ASPLUND))
 * @param       p_MaltsevFallbackFraction       User-specfifid Maltsev fallback fraction (`--maltsev-fallback-fraction`) [0.0, 1.0]
 * @param       p_MTcase                        Mass transfer case of first MT event (MT_CASE::NONE if no MT history) // m_StateHistory.CurrentState.FirstMTdonorEvent().MTcase <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 * @param       p_MaltsevMode                   User-specified Maltsev mode (`--maltsev-mode`))
 * @return                                      Remnant mass in Msol
 */
Dbl_DblT GiantBranch_Constituent::CalculateRemnantMass_Maltsev2025(const double       p_COCoreMass,
                                                                   const double       p_HeCoreMass,
                                                                   const double       p_ZetaAsplund, 
                                                                   const double       p_MaltsevFallbackFraction, //OPTIONS->MaltsevFallbackFraction()  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                                                   const MT_CASE      p_MTcase,
                                                                   const MALTSEV_MODE p_MaltsevMode) const { // OPTIONS->MaltsevMode() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    double fallbackFraction;
    double remnantMass;

    if (p_COCoreMass < MALTSEV2024_MMIN) {                                                  // CO core mass < Maltsev lower threshold?
                                                                                            // yes - NS formation regardless of metallicity
        fallbackFraction = 0.0;
        remnantMass = NEUTRON_STAR_MASS;
    }
    else if (p_COCoreMass > MALTSEV2024_MMAX) {                                             // no - CO core mass > Maltsev upper threshold?
                                                                                            // yes - BH formation regardless of metallicity
        fallbackFraction = 1.0;
        remnantMass = p_HeCoreMass;
    }
    else {                                                                                  // no - calculate remnant mass

        double Zbound;                                                                      
        switch (p_MaltsevMode) {                                                            // which Maltsev mode?                                                                                 

            case MALTSEV_MODE::OPTIMISTIC:                                                  // OPTIMISTIC               
                Zbound = p_ZetaAsplund;                                                     // Asplund zeta (log10(Z / ZSOL_ASPLUND))
                break;

            case MALTSEV_MODE::BALANCED:                                                    // BALANCED                        
                Zbound = std::min(std::max(p_ZetaAsplund, -1.6989700043360185), 0.0);       // Asplund zeta clamped to [log10(1/50), log10(1)]
                break;

            case MALTSEV_MODE::PESSIMISTIC:                                                 // Pessimistic                          
                Zbound = std::min(std::max(p_ZetaAsplund, -1.0), 0.0);                      // Asplund zeta clamped to [log10(1/10), log10(1)]
                break;

            default:                                                                        // unexpected mode
                // the only way this can happen is if the MALTSEV_MODE passed to this function
                // isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the code passed a mode that this function doesn't account for, and
                // that should be flagged as an error and result in termination of the evolution of the star
                // or binary.
                // The correct fix for this is to add code to this function for the missing mode, or fix the
                // calling code to pass a mode that is handled by this function.
                THROW_ERROR(ERROR::UNEXPECTED_MALTSEV_MODE);                                // throw error
        }
       
        MT_CASE MTcase = p_MTcase;                                                          // copy given mass transfer case

        // if the star was stripped by winds, treat it as if it experienced Case B mass transfer
        if (MTcase == MT_CASE::NONE && m_Star->HydrogenAbundanceSurface() == 0.0) MTcase = MT_CASE::B;

        // if no history of MT, or if MT happens from naked He stars, WDs, etc.,
        // assume that the core properties are not affected
        if (p_MTcase == MT_CASE::NONE || p_MTcase == MT_CASE::OTHER) {
            // JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        }
        else {                                                                              // handle MT cases
            double M1 = 0.0, M2 = 0.0, M3 = 0.0;
            switch (MTcase) {

                case MT_CASE::A:                                                            // case A MT
                    M1 = MALTSEV2024_M1A + (MALTSEV2024_M1A - MALTSEV2024_M1AZ01) * Zbound;
                    M2 = MALTSEV2024_M2A + (MALTSEV2024_M2A - MALTSEV2024_M2AZ01) * Zbound;
                    M3 = MALTSEV2024_M3A + (MALTSEV2024_M3A - MALTSEV2024_M3AZ01) * Zbound;
                    break;

                case MT_CASE::B:                                                            // case B MT
                    M1 = MALTSEV2024_M1B + (MALTSEV2024_M1B - MALTSEV2024_M1BZ01) * Zbound;
                    M2 = MALTSEV2024_M2B + (MALTSEV2024_M2B - MALTSEV2024_M2BZ01) * Zbound;
                    M3 = MALTSEV2024_M3B + (MALTSEV2024_M3B - MALTSEV2024_M3BZ01) * Zbound;
                    break;

                case MT_CASE::C:                                                            // case C MT
                    M1 = MALTSEV2024_M1C + (MALTSEV2024_M1C - MALTSEV2024_M1CZ01) * Zbound;
                    M2 = MALTSEV2024_M2C + (MALTSEV2024_M2C - MALTSEV2024_M2CZ01) * Zbound;
                    M3 = MALTSEV2024_M3C + (MALTSEV2024_M3C - MALTSEV2024_M3CZ01) * Zbound;
                    break;

                default:                                                                    // unexpected MT_CASE
                    // the only way this can happen is if the MT_CASE passed to this function is not accounted
                    // for in this code.  We should not default here, with or without a warning.
                    // We are here because the code passed a mass transfer case that this function doesn't account
                    // for, and that should be flagged as an error and result in termination of the evolution of
                    // binary.
                    // The correct fix for this is to add code to this function for the missing mass transfer case,
                    // or fix the calling code to pass a mass transfer case that is handled by this function.

                    THROW_ERROR(ERROR::UNEXPECTED_MT_CASE);                                 // throw error 
            }
        

            if (p_COCoreMass >= M3 || (p_COCoreMass >= M1 && p_COCoreMass <= M2)) {         // complete fallback into BH
                fallbackFraction = 1.0;
                remnantMass = p_HeCoreMass;
            }
            else if (p_COCoreMass > M2 && p_COCoreMass < M3 && RAND->Random(0, 1) <= 0.1) { // partial fallback - BH formation
                fallbackFraction = p_MaltsevFallbackFraction;
                remnantMass = (p_HeCoreMass - NEUTRON_STAR_MASS) * fallbackFraction + NEUTRON_STAR_MASS;
            }
            else {                                                                          // no fallback - NS formation
                fallbackFraction = 0.0;
                remnantMass = NEUTRON_STAR_MASS;
            }
        }
    }

    return std::make_tuple(remnantMass, std::clamp(fallbackFraction, 0.0, 1.0));
}






/*
 * CalculateCriticalMassRatio_Claeys2014
 *
 * @brief
 * Calculate the critical mass ratio for unstable mass transfer, per Claeys et al. 2014.
 *
 * See e.g de Mink et al. 2013, Claeys et al. 2014, and Ge et al. 2010, 2015, 2020 for discussions.
 *
 * Assumes this star is the donor.
 * 
 * Critical mass ratio is defined as qCrit = Maccretor / Mdonor.
 * 
 *
 * double CalculateCriticalMassRatio_Claeys2014(const double p_Mass, const double p_HeCoreMass, const bool p_AccretorIsDegenerate) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_HeCoreMass                    Helium core mass of the star (Msol)
 * @param       p_AccretorIsDegenerate          Boolean indicating if accretor is degenerate
 * @return                                      Critical mass ratio for unstable MT 
 */
double GiantBranch_Constituent::CalculateCriticalMassRatio_Claeys2014(const double p_Mass, const double p_HeCoreMass, const bool p_AccretorIsDegenerate) const {

    double qCrit;
                                                                                                                            
    if (p_AccretorIsDegenerate) {                                                                           // degenerate accretor?
        qCrit = OPTIONS->MassTransferCriticalMassRatioGiantDegenerateAccretor();                            // yes - use option value
    }
    else {                                                                                                  // no - non-degenerate accretor 
        if (OPTIONS->OptionDefaulted("critical-mass-ratio-giant-non-degenerate-accretor")) {
                                                                                                            // no - calculate it per Claeys et al. 2014
            if (p_HeCoreMass > 0.0) {
                const double mRatio = p_HeCoreMass / p_Mass;
                qCrit = 2.13 / ( 1.67 - ZDEP->HurleyGBRadiusXExponent(m_Star->Metallicity()) + 2.0 * utils::intPow(mRatio, 5)) // Claeys et al. 2014, Table 2
            }
            else qCrit = 0.0;
        }
        else {                                                                                              // yes - user sepcified value
            qCrit = OPTIONS->MassTransferCriticalMassRatioGiantNonDegenerateAccretor();                     // use the option value specified
        }
    }

    return qCrit;
}









///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 ZETA CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateZetaAdiabatic_ConvectiveEnvelopeGiant
 *
 * @brief
 * Calculate the adiabatic exponent (donor radial response to mass loss), zeta, for
 * convective-envelope giant-like stars
 *
 *
 * double CalculateZetaAdiabatic_ConvectiveEnvelopeGiant(ZETA_PRESCRIPTION p_ZetaPrescription) const
 *
 * @param       p_ZetaPrescription              User-specified (or default) ZETA_PRESCRIPTION
 * @return                                      Adiabatic exponent, zeta
 */
double GiantBranch_Constituent::CalculateZetaAdiabatic_ConvectiveEnvelopeGiant(ZETA_PRESCRIPTION p_ZetaPrescription) {
    
    double zeta;

    switch (p_ZetaPrescription) {                               // which ZETA_PRESCRIPTION?
            
        case ZETA_PRESCRIPTION::HURLEY:                         // HURLEY
            // JR REVIEW: BaseStar::CalculateZetaAdiabatic_Hurley2002 takes (p_Mass, p_CoreMass).
            // The pre-refactor body was passing only (m_CoreMass) - missing arg.  Fixed here
            // by passing both Mass and CoreMass from m_Star.
            zeta = m_Star->CalculateZetaAdiabatic_Hurley2002(m_Star->Mass(), m_Star->CoreMass());
            break;
            
        case ZETA_PRESCRIPTION::ARBITRARY:                      // ARBITRARY
            zeta = OPTIONS->ZetaAdiabaticArbitrary();
            break;

        case ZETA_PRESCRIPTION::SOBERMAN:                       // SOBERMAN
            zeta = m_Star->CalculateZetaAdiabatic_Soberman1997(m_Star->Mass(), m_Star->CoreMass());
            break;

        default:                                                // unexpected prescription
            // the only way this can happen is if the ZETA_PRESCRIPTION passed to this function is
            // not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.

            THROW_ERROR(ERROR::UNEXPECTED_ZETA_PRESCRIPTION);   // throw error             
    }
    
    return zeta;
}



/*
 * CalculateZetaAdiabatic_ByEnvelopeType
 * 
 * @brief
 * Calculate the adiabatic exponent (donor radial response to mass loss), zeta,
 * based on the envelope type of the star.
 *
 * Hurley et al. 2000, eqs 97 & 98
 *
 *
 * double CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription)
 *
 * @param   [IN]    p_ZetaPrescription          ZETA_PRESCRIPTION
 * @return                                      mass-radius response exponent Zeta
 */
double GiantBranch_Constituent::CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription) {
    
    double zeta;
    
    switch (m_Star->DetermineEnvelopeType()) {                              // which envelope type?
        case ENVELOPE::RADIATIVE:                                           // RADIATIVE
            zeta = OPTIONS->ZetaRadiativeEnvelopeGiant();                   // use option value
            break;
            
        case ENVELOPE::CONVECTIVE:                                          // CONVECTIVE
            zeta = CalculateZetaAdiabatic_ConvectiveEnvelopeGiant(p_ZetaPrescription);  // same-class sibling
            break;

        case ENVELOPE::REMNANT:                                             // REMNANT
            zeta = 0.0;                                                     // always 0.0
            break;

        default:                                                            // unknown envelope type
            // the only way this can happen is if someone added a new envelope type (to ENVELOPE)
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because DetermineEnvelopeType() returned an envelope type this code doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution of the
            // star or binary.
            // The correct fix for this is to add code for the missing envelope type or, if the missing envelope
            // type is incorrect/superfluous, remove it from ENVELOPE.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_TYPE);                      // throw error
    }
    
    return zeta;
}










