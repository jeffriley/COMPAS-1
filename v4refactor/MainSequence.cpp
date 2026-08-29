#include "MainSequence.h"
#include "MS_gt_07.h"
#include "HG.h"

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateLuminosity_Hurley2000
 *
 * @brief
 * Calculate luminosity on the Main Sequence, per Hurley et al. 2000, eq 12
 *
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function nominally relies on the ZAMS values of one or more attributes of the star,
 * and should not be used if the ZAMS values are required but are not known.
 * 
 * The function *could* be called with something other than ZAMS values as parameters - it
 * is up to the caller what values are actually passed. In such cases the ZAMS warning does
 * not apply.
 * 
 *
 * double CalculateLuminosity_Hurley2000(
 *     const double       p_Metallicity,
 *     const double       p_Mass,
 *     const double       p_Age,
 *     const double       p_LZAMS,
 *     const double       p_MHook,
 *     const TimescalesT& p_tScales
 * ) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_LZAMS                         ZAMS luminosity of the star (Lsol)
 * @param       p_MHook                         Mass above which hook appears on MS (Msol)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      MS luminosity (Lsol))
 */
double MainSequence::CalculateLuminosity_Hurley2000(
    const double       p_Metallicity,
    const double       p_Mass,
    const double       p_Age,
    const double       p_LZAMS,
    const double       p_MHook,
    const TimescalesT& p_tScales
) const {

    auto& a = ZDEP->HurleyACoeffs(p_Metallicity);   // Hurley a coefficients
    auto& l = ZDEP->HurleyLConsts(p_Metallicity);   // Hurley luminosity constants

    // Calculate the Hurley luminosity alpha coefficient, alphaL, per Hurley et al. 2000, eqs 19a & 19b
    double alphaL;
         if (p_Mass <   0.5) alphaL = a[49];
    else if (p_Mass <   0.7) alphaL = a[49] + (5.0 * (0.3 - a[49]) * (p_Mass - 0.5));
    else if (p_Mass < a[52]) alphaL = 0.3 + ((a[50] - 0.3) * (p_Mass - 0.7) / (a[52] - 0.7));
    else if (p_Mass < a[53]) alphaL = a[50] + ((a[51] - a[50]) * (p_Mass - a[52]) / (a[53] - a[52]));
    else if (p_Mass <   2.0) alphaL = a[51] + ((p_LConsts[HURLEY_L::B_ALPHA_L] - a[51]) * (p_Mass - a[53]) / (2.0 - a[53]));
    else                     alphaL = (a[45] + (a[46] * PPOW(p_Mass, a[48]))) / (PPOW(p_Mass, 0.4) + (a[47] * PPOW(p_Mass, 1.9)));
    
    // Calculate the Hurley luminosity beta coefficient, betaL, per Hurley et al. 2000, eq 20
    double betaL = std::max(0.0, (a[54] - (a[55] * PPOW(p_Mass, a[56]))));
    if (p_Mass > a[57] && betaL > 0.0) {
        betaL = std::max(0.0, (l[HURLEY_L::B_BETA_L] - 10.0 * (p_Mass - a[57]) * l[HURLEY_L::B_BETA_L]));
    }

    // Calculate the Hurley luminosity perturbation value, deltaL, per Hurley et al. 2000, eq 16
    double deltaL;
         if (p_Mass <= p_MHook) deltaL = 0.0; // 0.0 in BSE Fortran code
    else if (p_Mass <  a[33]  ) deltaL = l[HURLEY_L::B_DELTA_L] * PPOW(((p_Mass - p_MHook) / (a[33] - p_MHook)), 0.4);
    else                        deltaL = std::min((a[34] / PPOW(p_Mass, a[35])), (a[36] / PPOW(p_Mass, a[37])));

    // Calculate the Hurley radius exponent eta, per Hurley et al. 2000, eq 18
    const double eta = p_Metallicity > 0.0009 ? 10.0 : (p_Mass > 1.0 ? (p_Mass >= 1.1 ? 20.0 : (100.0 * p_Mass) - 90.0) : 10.0);

    // Calculate luminosity

    constexpr double epsilon = 0.01;

    const double lTAMS = CalculateLuminosityAtPhaseEnd_Hurley2000(p_Mass);
    const double mu    = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10])))))));  // Hurley et al. 2000, eq 7
    const double tHook = mu * p_tScales[TS::BGB];                                                                                       // ibid., just after eq 5
    const double t     = p_Age / p_tScales[TS::MS];                                                                                     // ibid., eq 11
    const double t1    = std::min(1.0, (p_Age / tHook));                                                                                // ibid., eq 14
    const double t2    = std::clamp((p_Age - ((1.0 - epsilon) * tHook)) / (epsilon * tHook), 0.0, 1.0);                                 // ibid., eq 15

    // ibid., eq 12
    return p_LZAMS * PPOW(10.0, (alphaL * t) + (betaL * PPOW(t, eta)) + ((std::log10(lTAMS / p_LZAMS) - alphaL - betaL) * t * t) - (deltaL * ((t1 * t1) - (t2 * t2))));
}


/*
 * CalculateLuminosity_Shikauchi2024
 *
 * @brief
 * Calculate luminosity on the Main Sequence, as a function of current core mass and central
 * helium fraction, per Shikauchi et al. 2024, eq A5.
 * 
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function nominally relies on the ZAMS values of one or more attributes of the star,
 * and should not be used if the ZAMS values are required but are not known.
 * 
 * The function *could* be called with something other than ZAMS values as parameters - it
 * is up to the caller what values are actually passed. In such cases the ZAMS warning does
 * not apply.
 * 
 * 
 * Only valid during core hydrogen burning, and for:
 * 
 *     MZAMS (Msol) >= SHIKAUCHI_LOWER_MASS_LIMIT, or 
 *     MZAMS (Msol) >= max(SHIKAUCHI_LOWER_MASS_LIMIT, BRCEK_LOWER_MASS_LIMIT) when the
 *     BRCEK core mass prescription is specified.
 * 
 * This function returns a std::optional<double> value.  If the star is in the Shikauchi regime, the
 * function will calculate the MS luminosity per Shikauchi et al. 2024 and return the luminosity
 * (a .has_value() check of the returned value will return true), otherwise std::nullopt is returned
 * (a .has_value() check of the returned value will return false).
 * 
 * 
 * OptDblT CalculateLuminosity_Shikauchi2024(const double p_Metallicity, const double p_MZAMS, const double p_CoreMass, const double p_HeAbundanceCore) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core of the star [0, 1]
 * @return                                      Luminosity (Lsol)
 */
OptDblT double MainSequence::CalculateLuminosity_Shikauchi2024(const double p_Metallicity, const double p_MZAMS, const double p_CoreMass, const double p_HeAbundanceCore) const {

    if (p_HeAbundanceCore > 0.0) return std::nullopt;                       // Helium abundance not in Shikauchi regime

    // Determine the lower mass limit (see description above)
    const double lowerMassLimit = (p_MSCoreMassPrescription == MS_CORE_MASS_PRESCRIPTION::BRCEK)
                                    ? std::max(SHIKAUCHI_LOWER_MASS_LIMIT, BRCEK_LOWER_MASS_LIMIT)
                                    : SHIKAUCHI_LOWER_MASS_LIMIT;

    if (p_MZAMS >= lowerMassLimit) return std::nullopt;                     // MZAMS not in Shikauchi regime


    // In Shikauschi regime - calculate luminosity
    
    // Common factors
    const double logMixingCoreMass      = std::log10(p_CoreMass);
    const double logMixingCoreMass2     = logMixingCoreMass * logMixingCoreMass;
    const double logMixingCoreMass3     = logMixingCoreMass * logMixingCoreMass2;
    const double logMixingCoreMass4     = logMixingCoreMass * logMixingCoreMass3;
    
    const double HeAbundanceCore        = p_HeAbundanceCore;
    const double HeAbundanceCore2       = p_HeAbundanceCore * p_HeAbundanceCore;
    const double HeAbundanceCore3       = p_HeAbundanceCore * HeAbundanceCore2;
    const double HeAbundanceCore4       = p_HeAbundanceCore * HeAbundanceCore3;

    const double logMixingCoreMass_He   = logMixingCoreMass * p_HeAbundanceCore;
    const double logMixingCoreMass_He2  = logMixingCoreMass * HeAbundanceCore2;
    const double logMixingCoreMass_He3  = logMixingCoreMass * HeAbundanceCore3;

    const double He_LogMixingCoreMass2  = p_HeAbundanceCore * logMixingCoreMass2;
    const double He_LogMixingCoreMass3  = p_HeAbundanceCore * logMixingCoreMass3;

    const double logMixingCoreMass2_He2 = logMixingCoreMass2 * HeAbundanceCore2;

    // Log luminosity 
    auto& lCoeffs = ZDEP->ShikauchiLuminosityCoefficients(p_Metallicity);   // Shikauchi L coefficients
    const double logL = lCoeffs[0]  * logMixingCoreMass      + 
                        lCoeffs[1]  * p_HeAbundanceCore      + 
                        lCoeffs[2]  * logMixingCoreMass_He   + 
                        lCoeffs[3]  * logMixingCoreMass2     + 
                        lCoeffs[4]  * HeAbundanceCore2       + 
                        lCoeffs[5]  * logMixingCoreMass3     + 
                        lCoeffs[6]  * HeAbundanceCore3       + 
                        lCoeffs[7]  * He_LogMixingCoreMass2  + 
                        lCoeffs[8]  * logMixingCoreMass_He2  + 
                        lCoeffs[9]  * logMixingCoreMass4     + 
                        lCoeffs[10] * HeAbundanceCore4       + 
                        lCoeffs[11] * logMixingCoreMass_He3  + 
                        lCoeffs[12] * logMixingCoreMass2_He2 + 
                        lCoeffs[13] * He_LogMixingCoreMass3  + 
                        lCoeffs[14];
    
    return PPOW(10.0, logL);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRadiusAtPhaseEnd_Hurley2000
 *
 * @brief
 * Calculate radius at the end of the Main Sequence (TAMS),
 * per Hurley et al. 2000, eqs 9a & 9b
 *
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function nominally relies on the ZAMS values of one or more attributes of the star,
 * and should not be used if the ZAMS values are required but are not known.
 * 
 * The function *could* be called with something other than ZAMS values as parameters - it
 * is up to the caller what values are actually passed. In such cases the ZAMS warning does
 * not apply.
 * 
 *
 * double CalculateRadiusAtPhaseEnd_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_RZAMS) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @return                                      TAMS radius (Rsol)
 */
double MainSequence::CalculateRadiusAtPhaseEnd_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_RZAMS) const {
    
    auto& a = ZDEP->HurleyACoeffs(p_Metallicity);       // Hurley a coefficients

    const double mA = a[17] + 0.1;

    double rTAMS = 0.0;                                 // Default return value

    if (p_Mass <= a[17]) {
        rTAMS = (a[18] + (a[19] * PPOW(p_Mass, a[21]))) / (a[20] + PPOW(p_Mass, a[22]));
        if (p_Mass < 0.5) rTAMS = std::max(rTAMS, 1.5 * p_RZAMS);
    }
    else if (p_Mass >= mA) {     
        const double m3 = p_Mass * p_Mass * p_Mass;     // pow() is slow - use multiplication
        rTAMS = ((C_COEFF.at(1) * m3) + (a[23] * PPOW(p_Mass, a[26])) + (a[24] * PPOW(p_Mass, a[26] + 1.5))) / (a[25] + m3 * p_Mass * p_Mass);
    }
    else {                                              // Interpolate between the end points
        const double mA3 = mA * mA * mA;                // pow() is slow - use multiplication

        const double y2 = ((C_COEFF.at(1) * mA3) + (a[23] * PPOW(mA, a[26])) + (a[24] * PPOW(mA, a[26] + 1.5))) / (a[25] + mA3 * mA * mA);
        const double y1 = (a[18] + (a[19] * PPOW(a[17], a[21]))) / (a[20] + PPOW(a[17], a[22]));

        const double gradient  = (y2 - y1) / 0.1;
        const double intercept = y1 - (gradient * a[17]);

        rTAMS = (gradient * p_Mass) + intercept;
    }

    return rTAMS;
}


/*
 * CalculateRadius_Hurley2000
 *
 * @brief
 * Calculate radius on the Main Sequence per Hurley et al. 2000, eq 13
 *
 *
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function nominally relies on the ZAMS values of one or more attributes of the star,
 * and should not be used if the ZAMS values are required but are not known.
 * 
 * The function *could* be called with something other than ZAMS values as parameters - it
 * is up to the caller what values are actually passed. In such cases the ZAMS warning does
 * not apply.
 * 
 * 
 * double CalculateRadius_Hurley2000(
 *     const double p_Metallicity,
 *     const double p_Mass,
 *     const double p_Tau,
 *     const double p_RZAMS,
 *     const double p_MHook,
 *     const double p_tBGB
 * ) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_MHook                         Mass above which hook appears on MS (Msol)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (Myr)
 * @return                                      MS radius (Rsol))
 */
double MainSequence::CalculateRadius_Hurley2000(
    const double p_Metallicity,
    const double p_Mass,
    const double p_Tau,
    const double p_RZAMS,
    const double p_MHook,
    const double p_tBGB
) const {

    auto& a = ZDEP->HurleyACoeffs(p_Metallicity);   // Hurley a coefficients
    auto& g = ZDEP->HurleyGConsts(Metallicity());   // Hurley gamma constants
    auto& r = ZDEP->HurleyRConsts(Metallicity());   // Hurley radius constants

    // Calculate the Hurley radius constant, alphaR, per Hurley et al. 2000, eqs 21a & 21b   
    double alphaR;
         if (p_Mass <    0.5) alphaR = a[62];
    else if (p_Mass <   0.65) alphaR = a[62] + (a[63] - a[62]) * (p_Mass - 0.5) / 0.15;
    else if (p_Mass <  a[68]) alphaR = a[63] + (a[64] - a[63]) * (p_Mass - 0.65) / (a[68] - 0.65);
    else if (p_Mass <  a[66]) alphaR = a[64] + (r[HURLEY_R::B_ALPHA_R] - a[64]) * (p_Mass - a[68]) / (a[66] - a[68]);
    else if (p_Mass <= a[67]) alphaR = a[58] * PPOW(p_Mass, a[60]) / (a[59] + PPOW(p_Mass, a[61]));
    else                      alphaR = r[HURLEY_R::C_ALPHA_R] + a[65] * (p_Mass - a[67]);

    // Calculate the Hurley radius constant, betaR, per Hurley et al. 2000, eqs 22a & 22b
    double betaR;
         if (p_Mass <=  1.0) betaR = 1.06;
    else if (p_Mass < a[74]) betaR = 1.06 + (a[72] - 1.06) * (p_Mass - 1.0) / (a[74] - 1.06);
    else if (p_Mass <   2.0) betaR = a[72] + (r[HURLEY_R::B_BETA_R] - a[72]) * (p_Mass - a[74]) / (2.0 - a[74]);
    else if (p_Mass <= 16.0) betaR = (a[69] * p_Mass * p_Mass * p_Mass * std::sqrt(p_Mass)) / (a[70] + PPOW(p_Mass, a[71]));
    else                     betaR = r[HURLEY_R::C_BETA_R] + a[73] * (p_Mass - 16.0);

    betaR -= 1.0;

    // Calculate the Hurley radius perturbation value, deltaR, per Hurley et al. 2000, eq 17
    double deltaR;
         if (p_Mass <= p_MHook) deltaR = 0.0; // 0.0 in BSE Fortran code
    else if (p_Mass <= a[42]  ) deltaR = a[43] * std::sqrt((p_Mass - p_MHook) / (a[42] - p_MHook));
    else if (p_Mass <    2.0  ) deltaR = a[43] + ((r[HURLEY_R::B_DELTA_R] - a[43]) * PPOW(((p_Mass - a[42]) / (2.0 - a[42])), a[44]));
    else                        deltaR = (a[38] + (a[39] * p_Mass * p_Mass * p_Mass * std::sqrt(p_Mass)) / ((a[40] * p_Mass * p_Mass * p_Mass) + PPOW(p_Mass, a[41]))) - 1.0;

    // Calculate the Hurley radius coefficient gamma, per Hurley et al. 2000, eq 23
    double gamma;
    if (p_Mass <= 1.0) {
        // BSE Fortran code has abs()
        gamma = a[76] + (a[77] * PPOW(std::abs(p_Mass - a[78]), a[79]));
    }
    else if (p_Mass <= a[75]) {
        gamma = g[HURLEY_G::B_GAMMA] + (a[80] - g[HURLEY_G::B_GAMMA]) * PPOW((p_Mass - 1.0) / (a[75] - 1.0), a[81]);
    }
    else if (p_Mass <= (a[75] + 0.1)) {
        // See discussion just prior to eq 23 - the end point is wrong in the arXiv version of
        // Hurley et al. 2000 (should be 0.1, not 1.0) - confirmed in BSE Fortran code
        gamma = g[HURLEY_G::C_GAMMA] - (10.0 * (p_Mass - a[75]) * g[HURLEY_G::C_GAMMA]);
    }
    else {
        // See discussion just prior to eq 23 (confirmed in BSE Fortran code)
        gamma = 0.0;
    }
    
    // See discussion following eq 23 (confirmed in BSE Fortran code)
    gamma = std::max(0.0, gamma); 

    // Calculate radius per Hurley et al. 2000
    const double tBGB   = astro::CalculateLifetimeToBGB_Hurley2000(p_Metallicity, p_Mass);
    const double tMS    = CalculatePhaseLifetime_Hurley2000(p_Metallicity, p_Mass, p_tBGB);
    const double rTMS   = CalculateRadiusAtPhaseEnd_Hurley2000(p_Metallicity, p_Mass, p_RZAMS);
    const double mu     = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10])))))));
    const double tHook  = mu * p_tBGB;                                                                  // Just after eq 5
    const double time   = tMS * p_Tau;
    const double tau1   = std::min(1.0, (time / tHook));                                                // Eq 14
    const double tau2   = std::clamp((time - ((1.0 - 0.01) * tHook)) / (0.01 * tHook), 0.0, 1.0);       // Eq 15, epsilon= 0.01
    const double tau_3  = p_Tau * p_Tau * p_Tau;                                                        // pow() is slow - use multiplication where it makes sense
    const double tau_10 = p_Tau < FLOAT_TOLERANCE_ABSOLUTE ? 0.0: p_Tau * tau_3 * tau_3 * tau_3;        // Tolerance comparison to avoid underflow
    const double tau_40 = tau_10 < FLOAT_TOLERANCE_ABSOLUTE ? 0.0: tau_10 * tau_10 * tau_10 * tau_10;   // Tolerance comparison to avoid underflow
    const double tau1_3 = tau1 * tau1 * tau1;
    const double tau2_3 = tau2 * tau2 * tau2;

    // Eq 13
    const double logRMS_RZAMS = alphaR * p_Tau + betaR * tau_10 + gamma * tau_40 + (std::log10(rTMS / p_RZAMS) - alphaR - betaR - gamma) * tau_3 - deltaR * (tau1_3 - tau2_3);

    return p_RZAMS * PPOW(10.0, logRMS_RZAMS);
}

    






// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< use below to calculate abundances <<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< then delete this !!!
/*
 * CalculateRadiusOnPhase_Brcek
 *
 * @brief
 * Calculate radius on the Main Sequence per Brcek et al. 2025
 *
 * Calculate radius to smoothly connect the beginning of MS hook and the beginning of HG (@TAMS)
 * The MS hook starts at tau = 0.99 on the main sequence
 *
 * double CalculateRadiusOnPhase_Brcek(const double p_Mass, const double p_Tau, const double p_RZAMS, const double p_Radius, const DBL_VECTOR& p_bCoefficients) const
 *
 * @param       p_Mass                          Current mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_Radius                        Current radius of the star (Rsol)
 * @param       p_bCoefficients                 Hurley b(n) coefficients
 * @return                                      Main sequence radius of the star (Rsol)
 */
double MainSequence::CalculateRadiusOnPhase_Brcek(const double p_Mass, const double p_Tau, const double p_RZAMS, const double p_Radius, const DBL_VECTOR& p_bCoefficients) const {

    const double radiusTAMS        = HG::CalculateRadiusOnPhase_Hurley2000(p_Mass, p_Tau, p_RZAMS, p_bCoefficients); // JR FIX THIS: check <<<<<<<<<<<<<<<<<<
    const double radiusAtHookStart = CalculateRadius_Hurley2000(p_Mass, 0.99, p_RZAMS, 0.0);  // JR FIX THIS: check <<<<<<<<<<<<<<<<<<
        
    double radius = (radiusAtHookStart * (1.0 - p_Tau) + radiusTAMS * (p_Tau - 0.99)) / 0.01;   // linear interpolation

    // if the star has been stripped below its initial core mass we need to adjust the radius
  



    double heliumAbundanceSurface = HeAbundanceSurface();
        if (p_Mass < CNOProcessedCoreMass())
            // By tracing the helium profile in the star, calculate how the surface helium abundance changes if mass drops below the initial core mass
            heliumAbundanceSurface = HeAbundanceCoreOut() + (p_Mass - MainSequenceCoreMass()) * (HeAbundanceSurface() - HeAbundanceCoreOut()) / (CNOProcessedCoreMass() - MainSequenceCoreMass());
        
        // Factor that scales radius based on surface helium abundance
        double surfaceAbundanceFactor = (utils::Compare(HeAbundanceCore(), InitialHeAbundance()) != 0) ? (heliumAbundanceSurface - InitialHeAbundance()) / (HeAbundanceCore() - InitialHeAbundance()) : 0.0;
        
        radius += (p_RZAMS - radius) * surfaceAbundanceFactor;






    return radius;
}















///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////












// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< use below to calculate abundances <<<<<<<<<<<<<<<<<<<<<<
// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< then delete this ???
/*
 * CalculateCoreMass_Brcek
 *
 * @brief
 * Calculate the convective core mass and central helium fraction of a main sequence star that
 * loses mass, per Shikauchi et al. 2024
 *
 * ///////////////////////////// do this next bit in a separate function
 * This function also accounts for mass gain by modelling rejuvenation and updates the initial mixing core mass
 * and the helium abundance just outside the core. 
 * 
 * If star undergoes significant stripping and the total mass
 * reaches the initial convective core mass, surface helium and hydrogen abundances are also updated
 * ////////////////////////////////
 * 
 *
 * double CalculateCoreMass_Brcek(const double p_dt, const double p_dMdt) const
 *
 * @param       p_Mass              Current mass of star (Msol)
 * @param       p_dt                Time step (Myr)
 * @param       p_dMdt              Mass loss rate in (Msol yr^-1)   //// >>>>>> CHECK THIS (negative for mass loss, positive for mass gain) <<<<<<<<<<<<<<<<
 * @param       p_CoreMass          Current MS core mass of star (Msol)
 * @return                          MS convective core mass
 */
double MainSequence::CalculateCoreMass_Brcek(
                                            const double p_Mass,
                                             const double p_dt,
                                             const double p_dMdt,
                                             const double p_CoreMass,
                                            const double p_CNOprocessedCoreMass,
                                            const double p_MSCoreMass,
                                            const double p_InitialHeAbundance,
                                            const double p_HeAbundanceCore,
                                            const double p_HeAbundanceOutsideCore,
                                            const double p_HeAbundanceSurface
                                            ) const { 
// check dt = 0 - do nothing

    DBL_VECTOR aCoeffs = ZDEP->ShikauchiAlphaCoefficients(Metallicity());                                               // get Shikauchi alpha coefficients from GLOBALS
    DBL_VECTOR fCoeffs = ZDEP->ShikauchifMixCoefficients(Metallicity());                                                // get Shikauchi fMix coefficients from GLOBALS
    DBL_VECTOR dCoeffs = SHIKAUCHI_DELTA_COEFFICIENTS;                                                                  // get Shikauchi delta coefficients from constants.h


    auto fmix = [&](double mass) { return fCoeffs[0] + fCoeffs[1] * std::exp(-mass / fCoeffs[2]); };                    // Shikauchi et al. 2024, eq A3


    const double updatedMass = p_Mass + (p_dMdt * p_dt * MYR_TO_YEAR);                                                 // mass adjusted for loss/gain

    double delta;                                                                                                       // mixing core inertia
    if (p_dMdt <= 0.0) {                                                                                                // mass loss?
                                                                                                                        // yes
        const double g  = dCoeffs[1] * p_MSCoreMass + dCoeffs[2];                                                       // ibid., eq A7
        const double Yc = -dCoeffs[0] * (p_HeAbundanceCore - p_InitialHeAbundance) / (1.0 - p_InitialHeAbundance - Metallicity());
        delta = std::min(PPOW(10.0, Yc + g), 1.0);                                                                      // ibid., eq A6
    }
    else {                                                                                                              // no - mass gain
        delta = PPOW(2.0, -(p_HeAbundanceCore - p_InitialHeAbundance) / (1.0 - p_InitialHeAbundance - Metallicity()));  // updated prescription for mass gain
    }

    const double deltaCoreMassML      = p_MSCoreMass * delta * (updatedMass * fmix(updatedMass) / (p_Mass * fmix(p_Mass)) - 1.0); // core mass delta due to mass loss/gain
    const double alpha                = PPOW(10.0, std::max(-2.0, aCoeffs[1] * p_MSCoreMass + aCoeffs[2])) + aCoeffs[0];// ibid., eq (A2)
    const double deltaYc              = CalculateLuminosityOnPhase() / (Q_CNO * p_MSCoreMass) * p_dt;                   // ibid., eq (12) (central helium fraction delta)
    const double deltaCoreMassNatural = -alpha / (1 - alpha * p_HeAbundanceCore) * deltaYc * p_MSCoreMass;              // ibid., eq (4)  (core mass delta due to natural decay)
    const double deltaCoreMass        = deltaCoreMassNatural + deltaCoreMassML;                                         // total core mass delta
    
    const double mixingCoreMass   = std::min(p_MSCoreMass + deltaCoreMass, BRCEK_CORE_MASS_TO_MASS_RATIO_LIMIT * updatedMass); // mixing core mass, must be <= total mass
    double centralHeFraction      = std::min(p_HeAbundanceCore + deltaYc, 1.0 - Metallicity());                         // central helium fraction, clamped to a maximum of (1-Z)
    double HeAbundanceOutsideCore = p_HeAbundanceOutsideCore;                                                           // He abundance outside the core
    [[maybe_unused]] double CNOprocessedCoreMass   = p_CNOprocessedCoreMass;                                                             // CNO processed core mass

    if (deltaCoreMass > 0.0) {                                                                                          // core mass increased?
                                                                                                                        // yes - need to account for rejuvenation
        if (mixingCoreMass < p_CNOprocessedCoreMass) {                                                                  // mixing core mass < CNO processed core mass?
                                                                                                                        // yes
            const double f1 = p_HeAbundanceOutsideCore - p_InitialHeAbundance;
            const double f2 = p_MSCoreMass - p_CNOprocessedCoreMass;
            const double f3 = p_MSCoreMass + deltaCoreMass;

            // calculate the change in core helium abundance, assuming a linear profile between the
            // central helium fraction (Yc) and the initial helium fraction (Y0), and that the the
            // accreted gas has helium fraction Y0
            const double deltaY = (p_HeAbundanceOutsideCore - p_HeAbundanceCore) / f3 * deltaCoreMass + 0.5 / f3 * f1 / f2 * deltaCoreMass * deltaCoreMass;

            centralHeFraction       = p_HeAbundanceCore + deltaY;
            HeAbundanceOutsideCore += f1 / f2 * deltaCoreMass;
        }
        else {                                                                                                          // no - mixing core mass >= CNO processed core mass
            const double deltaCoreMass1 = p_CNOprocessedCoreMass - p_MSCoreMass;                                        // mass accreted up to the CNO processed core mass
            const double deltaCoreMass2 = deltaCoreMass - deltaCoreMass1;                                               // remaining accreted mass

            centralHeFraction      = (p_MSCoreMass * p_HeAbundanceCore + 0.5 * (p_HeAbundanceOutsideCore + p_InitialHeAbundance) * deltaCoreMass1 + deltaCoreMass2 * p_InitialHeAbundance) / (p_MSCoreMass + deltaCoreMass);
            HeAbundanceOutsideCore = p_InitialHeAbundance;
            CNOprocessedCoreMass   = mixingCoreMass;
        }
    }
    else {                                                                                                              // no - core mass decreased (core decayed)
        // if total mass dropped below the CNO processed core mass, partially processed material
        // is exposed and surface abundance needs to be adjusted
        if (updatedMass < p_CNOprocessedCoreMass) {                                                                     // updated total mass < CNO processed core mass?
                                                                                                                        // yes
            // calculate surface helium and hydrogen abundances
            const double HeAbundanceSurface = p_HeAbundanceOutsideCore + (updatedMass - p_MSCoreMass) * (p_HeAbundanceSurface - p_HeAbundanceOutsideCore) / (p_CNOprocessedCoreMass - p_MSCoreMass);
            [[maybe_unused]] const double HAbundanceSurface  = 1.0 - Metallicity() - HeAbundanceSurface;

            CNOprocessedCoreMass = updatedMass;
        }

        HeAbundanceOutsideCore = centralHeFraction;
    }
    
    // JR FIX THIS: check vvvvvvvvvvvvvvvvvvvvvvvvvvv  <<<<<<<<<<<<<<<<<<<<<<<<<<
    // also return HeAbundanceOutsideCore, CNOprocessedCoreMass, HeAbundanceSurface, and HAbundanceSurface
    return mixingCoreMass;
}





/*
 * CalculateCoreMass_Brcek2025
 *
 * @brief
 * Calculate the convective core mass of a main sequence star, per Brcek et al. 2025
 * (follows Shikauchi et al. 2024)
 *
 *
 * double CalculateCoreMass_Brcek2025(
 *     const double p_Metallicity,
 *     const double p_Mass,
 *     const double p_Luminosity,
 *     const double p_MZAMS,
 *     const double p_CoreMass,
 *     const double p_HeAbundanceCore,
 *     const double p_dt,
 *     const double p_dMdt,
 * ) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core of the star [0, 1]
 * @param       p_dt                            Time step (Myr)
 * @param       p_dMdt                          Mass loss rate (Msol yr^-1)
 * @return                                      MS convective core mass (Msol)
 */
COMPAS_PURE double MainSequence::CalculateCoreMass_Brcek2025(
    const double p_Metallicity,
    const double p_Mass,
    const double p_Luminosity,
    const double p_MZAMS,
    const double p_CoreMass,
    const double p_HeAbundanceCore,
    const double p_dt,
    const double p_dMdt
) const {
     
    double coreMass = p_Coremass;                                                                                           // Default is no change

    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS - IS THIS CHECK REALLY REQUIRED??? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (p_dt <= 0.0) return coreMass;                                                                                       // Sanity check - no change if 0 timestep

    if (p_MZAMS >= BRCEK_LOWER_MASS_LIMIT) {                                                                                // In BRCEK regime?
                                                                                                                            // Yes
        if (p_HeAbundanceCore < (1.0 - p_Metallicity)) {                                                                    // In MS hook?
                                                                                                                            // Yes
            auto& aCoeffs = std::get<const DblVectorT&>(ZDEP->ShikauchiAlphaCoefficients(p_Metallicity));                   // Shikauchi alpha coefficients
            auto& fCoeffs = std::get<const DblVectorT&>(ZDEP->ShikauchifMixCoefficients(p_Metallicity));                    // Shikauchi fMix coefficients
            auto& dCoeffs = SHIKAUCHI_DELTA_COEFFICIENTS;                                                                   // Shikauchi delta coefficients

            const double ZAMSHeAbundance = ZDEP->ZAMSHeAbundance(p_Metallicity);                                            // ZAMS helium abundance

            auto fmix = [&](double mass) { return fCoeffs[0] + fCoeffs[1] * std::exp(-mass / fCoeffs[2]); };                // Shikauchi et al. 2024, eq A3

            double delta;                                                                                                   // Mixing core inertia
            if (p_dMdt <= 0.0) {                                                                                            // Mass loss?
                                                                                                                            // Yes
                const double g  = dCoeffs[1] * p_CoreMass + dCoeffs[2];                                                     // ibid., eq A7
                const double Yc = -dCoeffs[0] * (p_HeAbundanceCore - ZAMSHeAbundance) / (1.0 - ZAMSHeAbundance - p_Metallicity);
                delta = std::min(PPOW(10.0, Yc + g), 1.0);                                                                  // ibid., eq A6
            }
            else {                                                                                                          // No - mass gain
                delta = PPOW(2.0, -(p_HeAbundanceCore - ZAMSHeAbundance) / (1.0 - ZAMSHeAbundance - p_Metallicity));
            }

            const double mAdjusted    = p_Mass + (p_dMdt * p_dt * MYR_TO_YEAR);                                             // Mass adjusted for loss/gain
            const double deltaMcML    = p_CoreMass * delta * (mAdjusted * fmix(mAdjusted) / (p_Mass * fmix(p_Mass)) - 1.0); // Core mass delta due to mass loss/gain
            const double alpha        = PPOW(10.0, std::max(-2.0, aCoeffs[1] * p_CoreMass + aCoeffs[2])) + aCoeffs[0];      // ibid., eq A2
            const double deltaYc      = p_Luminosity / (Q_CNO * p_CoreMass) * p_dt;                                         // ibid., eq 12 (central helium fraction delta)
            const double deltaMcDecay = -alpha / (1 - alpha * p_HeAbundanceCore) * deltaYc * p_CoreMass;                    // ibid., eq 4  (core mass delta due to natural decay)
            const double deltaMc      = deltaMcDecay + deltaMcML;                                                           // Total core mass delta
    
            coreMass = std::min(p_CoreMass + deltaMc, BRCEK_CORE_MASS_TO_MASS_RATIO_LIMIT * mAdjusted);                     // Mixing core mass, clamped to total mass
        }
    }

    return coreMass;
}













/*
 * Calculate the Equilibrium zeta
 *
 *
 * double CalculateZetaEquilibrium
 *
 * @return                                      Equilibrium zeta
 */
double MainSequence_Constituent::CalculateZetaEquilibrium() {
    const double mass   = m_Star->Mass();
    const double radius = m_Star->Radius();
    double deltaMass           = -mass / 1.0E5;
    double radiusAfterMassGain = m_Star->CalculateRadiusOnMassChange(deltaMass);
    double zetaEquilibrium     = (radiusAfterMassGain - radius) / deltaMass * mass / radius;        // dlnR / dlnM

    return zetaEquilibrium;
}
    






///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              MISCELLANEOUS FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////













///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



/*
 * Resolve changes to the remnant after the star loses its envelope
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
 *
 * STELLAR_TYPE ResolveEnvelopeLoss(bool p_Force)
 *
 * @param   [IN]    p_Force                     Boolean to indicate whether the resolution of the loss of the envelope should be performed
 *                                              without checking the precondition(s).
 *                                              Default is false.
 *
 * @return                                      Stellar type to which star should evolve
 */
STELLAR_TYPE MainSequence::ResolveEnvelopeLoss(bool p_Force) {

    STELLAR_TYPE stellarType = m_StellarType;
    
    if (p_Force || utils::Compare(m_Mass, 0.0) <= 0) {      // envelope loss
        stellarType = STELLAR_TYPE::MASSLESS_REMNANT;
        m_Radius    = 0.0;
        m_Mass      = 0.0;
    }
    
    return stellarType;
}


/*
 * Sets the mass and age of a merger product of two main sequence stars
 * (note: treats merger products as main-sequence stars, not CHE, and does not check for MS_lte_07)
 *
 * Uses prescription of Wang et al., 2022, https://www.nature.com/articles/s41550-021-01597-5
 *
 * void UpdateAfterMerger(double p_Mass, double p_HydrogenMass)
 *
 * @param   [IN]    p_Mass                      New value of stellar mass
 * @param   [IN]    p_HydrogenMass              Desired value of hydrogen mass of merger remnant
 *
 */
void MainSequence::UpdateAfterMerger(double p_Mass, double p_HydrogenMass) {
    #define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function

    m_Mass                 = p_Mass;
    m_Mass0                = m_Mass;
    m_MainSequenceCoreMass = 0.0;
    
    double initialHydrogenFraction = m_InitialHydrogenAbundance;
    
    CalculateTimescales();
    CalculateGBParams();
            
    m_Tau = (initialHydrogenFraction - p_HydrogenMass / m_Mass) / initialHydrogenFraction;      // assumes uniformly mixed merger product and a uniform rate of H fusion on main sequence
    
    m_Age = m_Tau * timescales(tMS);
    
    m_HeliumAbundanceCore   = 1.0 - m_Metallicity - p_HydrogenMass / p_Mass;
    m_HydrogenAbundanceCore = 1.0 - m_Metallicity - m_HeliumAbundanceCore;
    
    m_HeliumAbundanceSurface   = m_HeliumAbundanceCore;                                         // abundances are the same throughout the star, assuming uniform mixing after merger
    m_HydrogenAbundanceSurface = m_HydrogenAbundanceCore;
    
    if ((OPTIONS->MainSequenceCoreMassPrescription() == CORE_MASS_PRESCRIPTION::BRCEK) && (utils::Compare(m_MZAMS, BRCEK_LOWER_MASS_LIMIT) >= 0)) {
        m_InitialMainSequenceCoreMass = CalculateInitialMainSequenceCoreMass(p_Mass, m_HeliumAbundanceCore);           // update initial mixing core mass
        m_MainSequenceCoreMass        = m_InitialMainSequenceCoreMass;                                                 // update core mass
    }
    
    UpdateAttributesAndAgeOneTimestep(0.0, 0.0, 0.0, true);
    
    #undef timescales
}

















////////////////////////////////////////////////// MainSequence_Constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
















/* 
 * CalculateCriticalMassRatio_Ge2020_Interpolate    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FIX DESCRIPTION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 *
 * @brief
 * Interpolate Ge+ critical mass ratio (aka qCrit), for H-rich stars, per Ge et al. 2020
 * 
 * Function takes input QCRIT_PRESCRIPTION, currently either of the prescriptions for critical mass ratios
 * from Ge et al. (2020), GE or GE_IC. The first is the full adiabatic response, the second assumes
 * artificially isentropic envelopes. From private communication with Ge, we have an updated datatable that
 * includes qCrit for fully conservative and fully non-conservative MT, so we now interpolate on those as well.
 *
 * Interpolation is done linearly in logM, logR, and logZ
 * 
 * double CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_MTefficiency                  Mass transfer accretion efficiency
 * @return                                      Critical mass ratio (aka qCrit)
 */ 
double MainSequence_Constituent::CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency) const {

    const double logMass   = std::log10(p_Mass);   
    const double logRadius = std::log10(p_Radius);

    // Iterate over the (currently 2) QCRIT_GE tables to get the qcrits at each metallicity
    std::vector<const GE_QCRIT_TABLE*> qCritTables = { &QCRIT_GE_LOW_Z, &QCRIT_GE_HIGH_Z };
    double qCritPerZ[qCritTables.size()];

    for (SizeT tblIdx = 0; tblIdx < qCritTables.size(); tblIdx++) {     // for each qCrit table

        const DBL_VECTOR massVec = std::get<0>(*(qCritTables[tblIdx])); // vector of masses from qCrit table

        SizeTVectorT indices = utils::BinarySearch(massVec, p_Mass);
        int lowMIdx = indices[0];                                       // bin lower bound index
        int uppMIdx = indices[1];                                       // bin upper bound index
    
        // if masses are out of range, set to endpoints
        if (lowMIdx == -1) {                                            // below minimum?
            lowMIdx = 0;                                                // yes - lower edge of lowest bin 
            uppMIdx = 1;                                                // upper edge of lowest bin
        } 
        else if (uppMIdx == -1) {                                       // no - above maximum? 
            lowMIdx = massVec.size() - 2;                               // yes - lower edge of highest bin 
            uppMIdx = massVec.size() - 1;                               // upper edge of highest bin
        } 

        // vector of radii and qCrits from qCrit table
        const GE_QCRIT_RADII_QCRIT_VECTOR& radiiQCritsVec = std::get<1>(*(qCritTables[tblIdx]));

        // (log) radii vectors for mass bin bounds
        const DBL_VECTOR& radiiVecLow = std::get<0>(radiiQCritsVec[lowMIdx]);
        const DBL_VECTOR& radiiVecUpp = std::get<0>(radiiQCritsVec[uppMIdx]);

        // qCrits vectors for mass bin bounds -- adapter from 7-tuple form to vec-of-vecs (indexed by GE_QCRIT_MODEL)
        const auto& tupLow = radiiQCritsVec[lowMIdx];
        const auto& tupUpp = radiiQCritsVec[uppMIdx];
        const std::vector<DBL_VECTOR> qCritsVecLow = {
            std::get<1>(tupLow), std::get<2>(tupLow), std::get<3>(tupLow),
            std::get<4>(tupLow), std::get<5>(tupLow), std::get<6>(tupLow)
        };  // order matches GE_QCRIT_MODEL enum (ST_FULL=0, ST_HALF=1, ST_NON_C=2, IC_FULL=3, IC_HALF=4, IC_NON_C=5)
        const std::vector<DBL_VECTOR> qCritsVecUpp = {
            std::get<1>(tupUpp), std::get<2>(tupUpp), std::get<3>(tupUpp),
            std::get<4>(tupUpp), std::get<5>(tupUpp), std::get<6>(tupUpp)
        };
 
        std::vector<double> qCritVecUppELowM;
        std::vector<double> qCritVecUppEUppM;
        std::vector<double> qCritVecLowELowM;
        std::vector<double> qCritVecLowEUppM;
        
        // set the appropriate qCrit vector based on QCRIT_PRESCRIPTION and MT efficiency
        switch (OPTIONS->QCritPrescription()) {                         // which QCRIT_PRESCRIPTION?

            case QCRIT_PRESCRIPTION::GE:
                if (p_MTefficiency > 0.5) {
                    qCritVecUppELowM = qCritsVecLow[static_cast<int>(GE_QCRIT_MODEL::ST_FULL)];
                    qCritVecUppEUppM = qCritsVecUpp[static_cast<int>(GE_QCRIT_MODEL::ST_FULL)];
                    qCritVecLowELowM = qCritsVecLow[static_cast<int>(GE_QCRIT_MODEL::ST_HALF)];
                    qCritVecLowEUppM = qCritsVecUpp[static_cast<int>(GE_QCRIT_MODEL::ST_HALF)];
                }
                else {
                    qCritVecUppELowM = qCritsVecLow[static_cast<int>(GE_QCRIT_MODEL::ST_HALF)];
                    qCritVecUppEUppM = qCritsVecUpp[static_cast<int>(GE_QCRIT_MODEL::ST_HALF)];
                    qCritVecLowELowM = qCritsVecLow[static_cast<int>(GE_QCRIT_MODEL::ST_NON_C)];
                    qCritVecLowEUppM = qCritsVecUpp[static_cast<int>(GE_QCRIT_MODEL::ST_NON_C)];
                }
                break;

            case QCRIT_PRESCRIPTION::GE_IC:
                if (p_MTefficiency > 0.5) {
                    qCritVecUppELowM = qCritsVecLow[static_cast<int>(GE_QCRIT_MODEL::IC_FULL)];
                    qCritVecUppEUppM = qCritsVecUpp[static_cast<int>(GE_QCRIT_MODEL::IC_FULL)];
                    qCritVecLowELowM = qCritsVecLow[static_cast<int>(GE_QCRIT_MODEL::IC_HALF)];
                    qCritVecLowEUppM = qCritsVecUpp[static_cast<int>(GE_QCRIT_MODEL::IC_HALF)];
                }
                else {
                    qCritVecUppELowM = qCritsVecLow[static_cast<int>(GE_QCRIT_MODEL::IC_HALF)];
                    qCritVecUppEUppM = qCritsVecUpp[static_cast<int>(GE_QCRIT_MODEL::IC_HALF)];
                    qCritVecLowELowM = qCritsVecLow[static_cast<int>(GE_QCRIT_MODEL::IC_NON_C)];
                    qCritVecLowEUppM = qCritsVecUpp[static_cast<int>(GE_QCRIT_MODEL::IC_NON_C)];
                }
                break;

            default:                                                    // unexpected prescription
                // the only way this can happen is if the QCRIT_PRESCRIPTION specified by the user is not
                // accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user specified a prescription that this function doesn't account
                // for, and that should be flagged as an error and result in termination of the evolution
                // of the star or binary.  We know this function does't account for all QCRIT_PRESCRIPTIONs,
                // but it *should* account for all Ge-related QCRIT_PRESCRIPTIONs.
                // The correct fix for this is to add code to this function for the missing prescription,
                // or fix the calling code to ensure this function is not called when the user specifies
                // a QCRIT_PRESCRIPTION that is not handled by this function.

                THROW_ERROR(ERROR::UNEXPECTED_QCRIT_PRESCRIPTION);      // throw error
        }

        indices = utils::BinarySearch(radiiVecLow, logRadius);          // find lower mass bound radii bin for p_Radius
        int lowRLowMIdx = indices[0];                                   // bin lower bound index
        int uppRLowMIdx = indices[1];                                   // bin upper bound index
    
        // if radii are out of range, set to endpoints
        if (lowRLowMIdx == -1) {                                        // below minimum?                                      
            lowRLowMIdx = 0;                                            // yes - lower edge of lowest bin  
            uppRLowMIdx = 1;                                            // upper edge of lowest bin 
        }
        else if (uppRLowMIdx == -1) {                                   // no - above maximum?                                             
            lowRLowMIdx = radiiVecLow.size() - 2;                       // yes - lower edge of highest bin 
            uppRLowMIdx = radiiVecLow.size() - 1;                       // upper edge of highest bin 
        }
    
        indices = utils::BinarySearch(radiiVecUpp, logRadius);          // find upper mass bound radii bin for p_Radius
        int lowRUppMIdx = indices[0];
        int uppRUppMIdx = indices[1];
    
        // if radii are out of range, set to endpoints
        if (lowRUppMIdx == -1) {                                        // below minimum? 
            lowRUppMIdx = 0;                                            // yes - lower edge of lowest bin  
            uppRUppMIdx = 1;                                            // upper edge of lowest bin 
        }
        else if (uppRUppMIdx == -1) {                                   // no - above maximum? 
            lowRUppMIdx = radiiVecUpp.size() - 2;                       // yes - lower edge of highest bin
            uppRUppMIdx = radiiVecUpp.size() - 1;                       // upper edge of highest bin
        }
    
        // Set the boundary points for the 2D interpolation
        const double rLowMLowR = radiiVecLow[lowRLowMIdx];
        const double rLowMUppR = radiiVecLow[uppRLowMIdx];
        const double rUppMlowR = radiiVecUpp[lowRUppMIdx];
        const double rUppMUppR = radiiVecUpp[uppRUppMIdx];

        const double rSlope1 = (rLowMUppR - logRadius) / (rLowMUppR - rLowMLowR);
        const double rSlope2 = (rUppMUppR - logRadius) / (rUppMUppR - rUppMlowR);
    
        const double qUppLowLow = qCritVecUppELowM[lowRLowMIdx];
        const double qUppLowUpp = qCritVecUppELowM[uppRLowMIdx];
        const double qUppUppLow = qCritVecUppEUppM[lowRUppMIdx];
        const double qUppUppUpp = qCritVecUppEUppM[uppRUppMIdx];
        const double qLowLowLow = qCritVecLowELowM[lowRLowMIdx];
        const double qLowLowUpp = qCritVecLowELowM[uppRLowMIdx];
        const double qLowUppLow = qCritVecLowEUppM[lowRUppMIdx];
        const double qLowUppUpp = qCritVecLowEUppM[uppRUppMIdx];
    
        // Interpolate on logR first, then logM, then on efficiency, using nearest neighbour for extrapolation

        // qCrit upper efficiency, lower mass
        double qCritUppELowM;
             if (logRadius < rLowMLowR) qCritUppELowM = qUppLowLow;
        else if (logRadius > rLowMUppR) qCritUppELowM = qUppLowUpp;
        else                            qCritUppELowM = qUppLowLow + rSlope1 * (qUppLowUpp - qUppLowLow);

        // qCrit upper efficiency, upper mass
        double qCritUppEUppM;
             if (logRadius < rUppMlowR) qCritUppEUppM = qUppUppLow;
        else if (logRadius > rUppMUppR) qCritUppEUppM = qUppUppUpp;
        else                            qCritUppEUppM = qUppUppLow + rSlope2 * (qUppUppUpp - qUppUppLow);

        // qCrit lower efficiency, lower mass
        double qCritLowELowM;
             if (logRadius < rLowMLowR) qCritLowELowM = qLowLowLow;
        else if (logRadius > rLowMUppR) qCritLowELowM = qLowLowUpp;
        else                            qCritLowELowM = qLowLowLow + rSlope1 * (qLowLowUpp - qLowLowLow);

        // qCrit lower efficiency, upper mass
        double qCritLowEUppM;
             if (logRadius < rUppMlowR) qCritLowEUppM = qLowUppLow;
        else if (logRadius > rUppMUppR) qCritLowEUppM = qLowUppUpp;
        else                            qCritLowEUppM = qLowUppLow + rSlope2 * (qLowUppUpp - qLowUppLow);
    
        const double logLowerMass = std::log10(massVec[lowMIdx]);
        const double logUpperMass = std::log10(massVec[uppMIdx]);
        const double mSlope       = (logUpperMass - logMass) / (logUpperMass - logLowerMass);

        double qCritUppE;
             if (logMass < logLowerMass) qCritUppE = qCritUppELowM;
        else if (logMass > logUpperMass) qCritUppE = qCritUppEUppM;
        else                             qCritUppE = qCritUppELowM + mSlope * (qCritUppEUppM - qCritUppELowM);

        double qCritLowE;
             if (logMass < logLowerMass) qCritLowE = qCritLowELowM;
        else if (logMass > logUpperMass) qCritLowE = qCritLowEUppM;
        else                             qCritLowE = qCritLowELowM + mSlope * (qCritLowEUppM - qCritLowELowM);
    
        qCritPerZ[tblIdx] = p_MTefficiency * qCritUppE + (1.0 - p_MTefficiency) * qCritLowE;
    }

    constexpr double logZlo = -3;                                       // log10(0.001)
    constexpr double logZhi = LOG10_ZSOL_HURLEY;                        // log10(0.02)  // JR FIX THIS <<<<<<<<<<<<<<<<<
    
    return qCritPerZ[1] + (ZDEP->SigmaHurley(m_Star->Metallicity()) - logZhi) * (qCritPerZ[1] - qCritPerZ[0]) / (logZhi - logZlo);
}
