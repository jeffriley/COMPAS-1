#include "MainSequence.h"
#include "MS_gt_07.h"
#include "HG.h"



// Here:
//
// Luminosity functions
// Mass functions
// Miscellaneous functions
//    -
//    -
// Radius functions








///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLuminosityAtPhaseEnd_Hurley2000
 *
 * @brief
 * Calculate luminosity at the end of the Main Sequence (TAMS),
 * per Hurley et al. 2000, eq 8
 *
 *
 * double CalculateLuminosityAtPhaseEnd_Hurley2000(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAMS luminosity (Lsol)
 */
COMPAS_PURE double MainSequence::CalculateLuminosityAtPhaseEnd_Hurley2000(const double p_Mass) const {
   
    const DBL_VECTOR a = GLOBALS->HurleyAcoefficients(); // get Hurley a coefficients

    const double m2 = p_Mass * p_Mass;
    const double m3 = p_Mass * m2;

    return (a[11] * m3) + (a[12] * m2 * m2) + (a[13] * PPOW(p_Mass, (a[16] + 1.8))) / (a[14] + (a[15] * m2 * m3) + PPOW(p_Mass, a[16]));
}


/*
 * CalculateLuminosity_Brcek2025
 *
 * @brief
 * Calculate luminosity on the Main Sequence per Brcek et al. 2025.
 * 
 * Uses Shikauchi et al. 2024 if possible, otherwise defaults to the method
 * appropriate for the evolution method and as specified by program options.
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
 * not apply.  This function exists in the main code rather than in the GLOBALS module in
 * recognition of the usefulness of such cases.
 * 
 *
 * double CalculateLuminosity_Brcek2025(
 *    const double p_Mass, 
 *    const double p_Tau,
 *    const double p_Time,
 *    const double p_MZAMS,
 *    const double p_LZAMS,
 *    const double p_CoreMass,
 *    const double p_HeAbundanceCore,
 *    const double p_tMS,
 *    const double p_tBGB
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_Time                          Time elapsed since ZAMS (Myr)
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_LZAMS                         ZAMS luminosity of the star (Lsol)
 * @param       p_CoreMass                      MS core mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core of the star 
 * @param       p_tMS                           MS lifetime, tMS (per Hurley timescales) (Myr)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (per Hurley timescales) (Myr)
 * @return                                      MS luminosity (Lsol)
 */
double MainSequence::CalculateLuminosity_Brcek2025(
    const double p_Mass, 
    const double p_Tau,
    const double p_Time,
    const double p_MZAMS,
    const double p_LZAMS,
    const double p_CoreMass,
    const double p_HeAbundanceCore,
    const double p_tMS,
    const double p_tBGB
) const {

    double luminosity;

    const double tHook = 0.99 * p_tMS;                                                      // MS hook start time

    if (p_Tau > tHook) {                                                                    // star in MS hook?
                                                                                            // yes
        // calculate luminosity on the transition from MS to HG - interpolate to smoothly
        // connect the beginning of MS hook and the beginning of HG (TAMS)

        double lHook;                                                                       // luminosity at the start of the MS hook
        if (p_MZAMS >= std::max(SHIKAUCHI_LOWER_MASS_LIMIT, BRCEK_LOWER_MASS_LIMIT)) {      // mass in Brcek/Shikauchi regime?
                                                                                            // yes - in the hook
            // core helium abundance fixed at 1-Z and core mass is not changing
            lHook = CalculateLuminosity_Shikauchi(p_CoreMass, p_HeAbundanceCore);
        }
        else {                                                                              // no - can't use Shikauchi here
                                                                                            // use default method
            lHook = CalculateLuminosity_Hurley2000(p_Mass, tHook, p_LZAMS, p_tMs, p_tBGB); // <<<<<<<<<<<<<<<< HURLEY -> generic <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        }
    
        const double lTAMS = CalculateLuminosityAtPhaseEnd_Hurley2000(p_Mass); // <<<<<<<<<<<<<<<< HURLEY -> generic <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        luminosity = (lHook * (tMS - p_Age) + lTAMS * (p_Age - tHook)) / (tMS - tHook);     // interpolate
    }
    else {                                                                                  // not in the MS hook
        if (utils::Compare(m_MZAMS, 15.0) >= 0) {                                           // MZAMS >= 15 Msol?
            luminosity = CalculateLuminosityShikauchi(p_CoreMass, p_HeAbundanceCore);       // yes, use Shikauchi
        }
        else {                                                                              // no, use default method
            luminosity = CalculateLuminosity_Hurley2000(p_Mass, p_Time, p_LZAMS, p_tMS, p_tBGB); // <<<<<<<<<<<<<<<< HURLEY -> generic <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        }
    }

    return luminosity;
}


/*
 * CalculateLuminosity
 *
 * @brief
 * Calculate luminosity on the Main Sequence.
 *
 * If the Brcek MS core mass prescription was specified by the user, and the star is in the 
 * Brcek regime, calculate the MS luminosity per Brcek et al. 2025, otherwise uses the method
 * appropriate for the evolution method and as specified by program options.
 * 
 * Uses globals and state variables - read only.
 * 
 * 
 * double CalculateLuminositye() const
 *
 * @return                                      MS luminosity (Lsol)
 */
double MainSequence::CalculateLuminosity() const {
        
    // common variables
    const double mZAMS          = MZAMS();                                                            // ZAMS mass of the star

    double luminosity;
    if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK &&      // Brcek MS core mass prescription specified, and ...
        mZAMS >= BRCEK_LOWER_MASS_LIMIT) {                                                      // ... in Brcek regime?
                                                                                                                                // yes - use Brcek
    const double Z              = GLOBALS->Metallicity();                                                              // metallicity of the star
    const double mass           = Mass();                                                         // current mass of the star
    const double time           = Time();                                                         // time elapsed since ZAMS
        const double tau                  = Tau();                                                // MS fractional age of the star
        const double coreMass             = CoreMass();                                           // current core mass of the star
        const double HeAbundanceCore      = HeAbundanceCore();                                    // core Helium abundance
        const DBL_VECTOR aCoeffs          = GLOBALS->HurleyACoefficients();                                                     // Hurley a(n) coefficients
        const DBL_VECTOR ShikauchiLCoeffs = GLOBALS->ShikauchiLCoefficients();                                                  // Shikauchi luminosity coefficients
    const double lZAMS          = LZAMS();                                                      // ZAMS luminosity of the star
    const DBL_VECTOR timescales = GLOBALS->Timescales();                                                                        // Phase timescales

        luminosity = CalculateLuminosity_Brcek(GLOBALS->Metallicity(),mass, tau, time, mZAMS, lZAMS, coreMass, HeAbundanceCore, timescales, aCoeffs, lConstants, ShikauchiLCoeffs);
    }
    else {                                                                                                                      // no - use default method
        luminosity = BaseStar::CalculateLuminosity();

        ////const DBL_VECTOR aCoeffs    = GLOBALS->HurleyACoefficients();                                                           // Hurley a(n) coefficients
        ////const DBL_VECTOR lConstants = GLOBALS->HurleyLConstants();                                                              // Hurley luminosity constants
        ////luminosity = CalculateLuminosity(Z, mass, time, lZAMS, timescales, aCoeffs, lConstants);
    }

    return luminosity;
}


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
 * double CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Age, const double p_LZAMS, const double p_tMS, const double p_tBGB) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_LZAMS                         ZAMS luminosity of the star (Lsol)
 * @param       p_tMS                           MS lifetime, tMS (per Hurley timescales) (Myr)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (per Hurley timescales) (Myr)
 * @return                                      MS luminosity (Lsol))
 */
double MainSequence::CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Age, const double p_LZAMS, const double p_tMS, const double p_tBGB) const {

    const DBL_VECTOR a = GLOBALS->HurleyAcoefficients(); // get Hurley a coefficients

    // Calculate the Hurley luminosity alpha coefficient, alphaL, per Hurley et al. 2000, eqs 19a & 19b
    double alphaL;
         if (p_Mass <   0.5) alphaL = a[49];
    else if (p_Mass <   0.7) alphaL = a[49] + (5.0 * (0.3 - a[49]) * (p_Mass - 0.5));
    else if (p_Mass < a[52]) alphaL = 0.3 + ((a[50] - 0.3) * (p_Mass - 0.7) / (a[52] - 0.7));
    else if (p_Mass < a[53]) alphaL = a[50] + ((a[51] - a[50]) * (p_Mass - a[52]) / (a[53] - a[52]));
    else if (p_Mass <   2.0) alphaL = a[51] + ((GLOBALS->HurleyLuminosityConstants(static_cast<int>(HURLEY_L_CONSTANTS::B_ALPHA_L)) - a[51]) * (p_Mass - a[53]) / (2.0 - a[53]));
    else                     alphaL = (a[45] + (a[46] * PPOW(p_Mass, a[48]))) / (PPOW(p_Mass, 0.4) + (a[47] * PPOW(p_Mass, 1.9)));
    
    // Calculate the Hurley luminosity beta coefficient, betaL, per Hurley et al. 2000, eq 20
    double betaL = std::max(0.0, (a[54] - (a[55] * PPOW(p_Mass, a[56]))));
    if (p_Mass > a[57] && betaL > 0.0) {
        const double bBetaL = GLOBALS->HurleyLuminosityConstants(static_cast<int>(HURLEY_L_CONSTANTS::B_BETA_L));
        betaL = std::max(0.0, (bBetaL - 10.0 * (p_Mass - a[57]) * bBetaL));
    }

    // Calculate the Hurley luminosity perturbation value, deltaL, per Hurley et al. 2000, eq 16
    const double mHook = GLOBALS->HurleyMassCutoffs(static_cast<int>(MASS_CUTOFF::MHook));

    double deltaL;
         if (p_Mass <= mHook) deltaL = 0.0; // 0.0 in BSE Fortran code
    else if (p_Mass <  a[33]) deltaL = GLOBALS->HurleyLuminosityConstants(static_cast<int>(HURLEY_L_CONSTANTS::B_DELTA_L)) * PPOW(((p_Mass - mHook) / (a[33] - mHook)), 0.4);
    else                      deltaL = std::min((a[34] / PPOW(p_Mass, a[35])), (a[36] / PPOW(p_Mass, a[37])));

    // Calculate the Hurley radius exponent eta, per Hurley et al. 2000, eq 18
    const double eta = GLOBALS->Metallicity() > 0.0009 ? 10.0 : (p_Mass > 1.0 ? (p_Mass >= 1.1 ? 20.0 : (100.0 * p_Mass) - 90.0) : 10.0);

    // Calculate luminosity

    constexpr double epsilon = 0.01;

    const double lTAMS = CalculateLuminosityAtPhaseEnd_Hurley2000(p_Mass);
    const double mu    = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10])))))));  // Hurley et al. 2000, eq 7
    const double tHook = mu * p_tBGB;                                                                                                   // ibid., just after eq 5
    const double t     = p_Age / p_tMS;                                                                                                 // ibid., eq 11
    const double t1    = std::min(1.0, (p_Age / tHook));                                                                                // ibid., eq 14
    const double t2    = std::max(0.0, std::min(1.0, (p_Age - ((1.0 - epsilon) * tHook)) / (epsilon * tHook)));                         // ibid., eq 15

    // ibid., eq 12
    return p_LZAMS * PPOW(10.0, (alphaL * t) + (betaL * PPOW(t, eta)) + ((std::log10(lTAMS / p_LZAMS) - alphaL - betaL) * t * t) - (deltaL * ((t1 * t1) - (t2 * t2))));
}


/*
 * CalculateLuminosity_Shikauchi2024
 *
 * @brief
 * Calculate luminosity on the Main Sequence, as a function of current core mass and central
 * helium fraction, per Shikauchi et al. 2024, eq A5.  Only valid during core hydrogen burning,
 * and for stars with MZAMS >= SHIKAUCHI_LOWER_MASS_LIMIT (Msol) (see constants.h).
 *
 * 
 * double CalculateLuminosity_Shikauchi2024(const double p_CoreMass, const double p_HeAbundanceCore) const
 *
 * @param       p_CoreMass                      MS core mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core of the star
 * @return                                      MS luminosity (Lsol)
 */
COMPAS_PURE double MainSequence::CalculateLuminosity_Shikauchi2024(const double p_CoreMass, const double p_HeAbundanceCore) const {

    const DBL_VECTOR ShikauchiLCoeffs = GLOBALS->ShikauchiLCoefficients(); // get Shikauchi L coefficients
    
    // common factors
    const double logMixingCoreMass      = std::log10(p_CoreMass);
    const double logMixingCoreMass2     = logMixingCoreMass * logMixingCoreMass;
    const double logMixingCoreMass3     = logMixingCoreMass * logMixingCoreMass2;
    const double logMixingCoreMass4     = logMixingCoreMass * logMixingCoreMass3;
    
    const double heliumAbundanceCore    = p_HeAbundanceCore;
    const double heliumAbundanceCore2   = heliumAbundanceCore * heliumAbundanceCore;
    const double heliumAbundanceCore3   = heliumAbundanceCore * heliumAbundanceCore2;
    const double heliumAbundanceCore4   = heliumAbundanceCore * heliumAbundanceCore3;

    const double logMixingCoreMass_He   = logMixingCoreMass * heliumAbundanceCore;
    const double logMixingCoreMass_He2  = logMixingCoreMass * heliumAbundanceCore2;
    const double logMixingCoreMass_He3  = logMixingCoreMass * heliumAbundanceCore3;

    const double He_LogMixingCoreMass2  = heliumAbundanceCore * logMixingCoreMass2;
    const double He_LogMixingCoreMass3  = heliumAbundanceCore * logMixingCoreMass3;

    const double logMixingCoreMass2_He2 = logMixingCoreMass2 * heliumAbundanceCore2;

    // log luminosity 
    const double logL = ShikauchiLCoeffs[0]  * logMixingCoreMass      + 
                        ShikauchiLCoeffs[1]  * heliumAbundanceCore    + 
                        ShikauchiLCoeffs[2]  * logMixingCoreMass_He   + 
                        ShikauchiLCoeffs[3]  * logMixingCoreMass2     + 
                        ShikauchiLCoeffs[4]  * heliumAbundanceCore2   + 
                        ShikauchiLCoeffs[5]  * logMixingCoreMass3     + 
                        ShikauchiLCoeffs[6]  * heliumAbundanceCore3   + 
                        ShikauchiLCoeffs[7]  * He_LogMixingCoreMass2  + 
                        ShikauchiLCoeffs[8]  * logMixingCoreMass_He2  + 
                        ShikauchiLCoeffs[9]  * logMixingCoreMass4     + 
                        ShikauchiLCoeffs[10] * heliumAbundanceCore4   + 
                        ShikauchiLCoeffs[11] * logMixingCoreMass_He3  + 
                        ShikauchiLCoeffs[12] * logMixingCoreMass2_He2 + 
                        ShikauchiLCoeffs[13] * He_LogMixingCoreMass3  + 
                        ShikauchiLCoeffs[14];
    
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
 * This function ostensibly relies on the value of the ZAMS radius of the star,
 * and should not be used if the ZAMS radius is not known.
 * 
 *
 * static double CalculateRadiusAtPhaseEnd_Hurley2000(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @return                                      TAMS radius (Rsol)
 */
double MainSequence::CalculateRadiusAtPhaseEnd_Hurley2000(const double p_Mass, const double p_RZAMS) const {
    
    const DBL_VECTOR a = GLOBALS->HurleyAcoefficients();    // get Hurley a coefficients

    const double mA = a[17] + 0.1;

    double rTAMS;

    if (p_Mass <= a[17]) {
        rTAMS = (a[18] + (a[19] * PPOW(p_Mass, a[21]))) / (a[20] + PPOW(p_Mass, a[22]));
        if (p_Mass < 0.5) rTAMS = std::max(rTAMS, 1.5 * p_RZAMS);
    }
    else if (p_Mass >= mA) {     
        const double m3 = p_Mass * p_Mass * p_Mass;         // pow() is slow - use multiplication
        rTAMS = ((HURLEY_C_COEFF[1] * m3) + (a[23] * PPOW(p_Mass, a[26])) + (a[24] * PPOW(p_Mass, a[26] + 1.5))) / (a[25] + m3 * p_Mass * p_Mass);
    }
    else {                                                  // interpolate between the end points
        const double mA3 = mA * mA * mA;                    // pow() is slow - use multiplication

        const double y2 = ((HURLEY_C_COEFF[1] * mA3) + (a[23] * PPOW(mA, a[26])) + (a[24] * PPOW(mA, a[26] + 1.5))) / (a[25] + mA3 * mA * mA);
        const double y1 = (a[18] + (a[19] * PPOW(a[17], a[21]))) / (a[20] + PPOW(a[17], a[22]));

        const double gradient  = (y2 - y1) / 0.1;
        const double intercept = y1 - (gradient * a[17]);

        rTAMS = (gradient * p_Mass) + intercept;
    }

    return rTAMS;
}


/*
 * CalculateRadius_Brcek2025
 *
 * @brief
 * Calculate radius on the Main Sequence per Brcek et al. 2025
 * Uses Shikauchi et al. 2024 if possible, otherwise defaults to the method
 * appropriate for the evolution method and as specified by program options.
 * 
 *
 * double CalculateRadius_Brcek2025(const double      p_Mass, 
 *                              const double      p_Tau,
 *                              const double      p_Age,
 *                              const double      p_MZAMS,
 *                              const double      p_LZAMS,
 *                              const double      p_CoreMass,
 *                              const double      p_HeAbundanceCore,
 *                              const DBL_VECTOR& p_Timescales,
 *                              const DBL_VECTOR& p_aN,
 *                              const DBL_VECTOR& p_LConsts,
 *                              const DBL_VECTOR& p_ShikauchiLCoeffs) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_LZAMS                         ZAMS luminosity of the star (Lsol)
 * @param       p_CoreMass                      MS core mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core of the star 
 * @param       p_Timescales                    Phase timescales
 * @param       p_aN                            Hurley a(n) coefficients
 * @param       p_LConsts                       Hurley luminosity constants
 * @param       p_ShikauchiLCoeffs              Shikauchi luminosity coefficients
 * @return                                      MS radius (Rsol)
 */
double MainSequence::CalculateRadius_Brcek2025(const double      p_Mass, 
                                           const double      p_Tau,
                                           const double      p_Age,
                                           const double      p_MZAMS,
                                           const double      p_LZAMS,
                                           const double      p_CoreMass,
                                           const double      p_HeAbundanceCore,
                                           const DBL_VECTOR& p_Timescales,
                                           const DBL_VECTOR& p_aN,
                                           const DBL_VECTOR& p_LConsts,
                                           const DBL_VECTOR& p_ShikauchiLCoeffs) const {

    double radius;

    // The core mass prescription from Shikauchi et al. (2024) cannot be used beyond the MS hook
    // (beyond Tau = 0.99), so if the star is in the MS hook, we return the radius that smoothly
    // connects the beginning of MS hook and the beginning of the HG.
    if (p_Tau > 0.99) {                                                                         // In MS hook?
                                                                                                // Yes
        // We need the TAMS radius, which is just the radius at the start of the HG phase.
        // Since we are on the main sequence here, we can clone this object as an HG object and,
        // as long as it is initialised (to correctly set Tau to 0.0 on the HG phase), we can 
        // query the cloned object for its radius.
        //
        // The clone should not evolve, and so should not log anything, but to be sure the
        // clone does not participate in logging, we set its persistence to EPHEMERAL.

        std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::HG, OBJECT_PERSISTENCE::EPHEMERAL);

        // Select radius at TAMS from the HG clone or current radius (whichever is smaller), relevant for
        // stars that were significantly stripped as this prevents radius expansion during the hook, and
        // delays possible mass transfer to the start of HG.
        const double radiusAtTAMS      = std::min(clone->Radius(), Radius());
        const double radiusAtHookStart = CalculateRadiusOnPhase(p_Mass, 0.99, p_RZAMS);         // Hook starts at Tau = 0.99
    
        radius = (radiusAtHookStart * (1.0 - p_Tau) + radiusAtTAMS * (p_Tau - 0.99)) / 0.01;    // Linear interpolation
    }
    else {                                                                                      // No - not in MS hook
                                                                                           
        radius = CalculateRadius(Mass(), Time(), LZAMS());  // Calculate radius using the default method <<<<<<<<<<<<<<<<<< JR FIX THIS

        // The star was stripped below its initial core mass, so the radius may need to be adjusted
        if (utils::Compare(m_HeAbundanceCore, m_InitialHeAbundance) != 0) {                     // Core helium abindance changed?
                                                                                                // Yes - scale radius based on surface helium abundance
            double HeAbundanceSurface;
            if (p_Mass < m_InitialMainSequenceCoreMass) {                                       // Mass dropped below the initial core mass?
                                                                                                // Yes - calculate new surface helium abundance
                HeAbundanceSurface = m_HeAbundanceCoreOut + (p_Mass - m_MSCoreMass) * (m_HeAbundanceSurface - m_HeAbundanceCoreOut) / (m_InitialMSCoreMass - m_MSCoreMass);
            }
            else HeAbundanceSurface = m_HeAbundanceSurface;                                     // No - use current surface helium abundance

            // Ajust radius
            radius += (p_RZAMS - radius) * (HeAbundanceSurface - m_InitialHeAbundance) / (m_HeAbundanceCore - m_InitialHeAbundance);
        }
    }
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
 * This function ostensibly relies on the value of the ZAMS radius of the star,
 * and should not be used if the ZAMS radius is not known.
 * 
 * 
 * double CalculateRadius_Hurley2000(const double p_Mass, const double p_Tau, const double p_RZAMS, const double p_tBGB) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (per Hurley timescales) (Myr)
 * @return                                      MS radius (Rsol)
 */
COMPAS_PURE double MainSequence::CalculateRadius_Hurley2000(const double p_Mass, const double p_Tau, const double p_RZAMS, const double p_tBGB) const {

    const DBL_VECTOR a = GLOBALS->HurleyAcoefficients(); // get Hurley a coefficients

    // calculate the Hurley radius constant, alphaR, per Hurley et al. 2000, eqs 21a & 21b   
    double alphaR;
         if (p_Mass <    0.5) alphaR = a[62];
    else if (p_Mass <   0.65) alphaR = a[62] + (a[63] - a[62]) * (p_Mass - 0.5) / 0.15;
    else if (p_Mass <  a[68]) alphaR = a[63] + (a[64] - a[63]) * (p_Mass - 0.65) / (a[68] - 0.65);
    else if (p_Mass <  a[66]) alphaR = a[64] + (GLOBALS->HurleyRadiusConstants(HURLEY_R::B_ALPHA_R) - a[64]) * (p_Mass - a[68]) / (a[66] - a[68]);
    else if (p_Mass <= a[67]) alphaR = a[58] * PPOW(p_Mass, a[60]) / (a[59] + PPOW(p_Mass, a[61]));
    else                      alphaR = GLOBALS->HurleyRadiusConstants(HURLEY_R::C_ALPHA_R) + a[65] * (p_Mass - a[67]);

    // calculate the Hurley radius constant, betaR, per Hurley et al. 2000, eqs 22a & 22b
    double betaR;
         if (p_Mass <=  1.0) betaR = 1.06;
    else if (p_Mass < a[74]) betaR = 1.06 + (a[72] - 1.06) * (p_Mass - 1.0) / (a[74] - 1.06);
    else if (p_Mass <   2.0) betaR = a[72] + (GLOBALS->HurleyRadiusConstants(HURLEY_R::B_BETA_R) - a[72]) * (p_Mass - a[74]) / (2.0 - a[74]);
    else if (p_Mass <= 16.0) betaR = (a[69] * p_Mass * p_Mass * p_Mass * std::sqrt(p_Mass)) / (a[70] + PPOW(p_Mass, a[71]));
    else                     betaR = GLOBALS->HurleyRadiusConstants(HURLEY_R::C_BETA_R) + a[73] * (p_Mass - 16.0);

    betaR -= 1.0;

    // calculate the Hurley radius perturbation value, deltaR, per Hurley et al. 2000, eq 17
    double deltaR;
         if (p_Mass <= mHook) deltaR = 0.0; // 0.0 in BSE Fortran code
    else if (p_Mass <= a[42]) deltaR = a[43] * std::sqrt((p_Mass - mHook) / (a[42] - mHook));
    else if (p_Mass <    2.0) deltaR = a[43] + ((GLOBALS->HurleyRadiusConstants(HURLEY_R::B_DELTA_R) - a[43]) * PPOW(((p_Mass - a[42]) / (2.0 - a[42])), a[44]));
    else                      deltaR = (a[38] + (a[39] * p_Mass * p_Mass * p_Mass * std::sqrt(p_Mass)) / ((a[40] * p_Mass * p_Mass * p_Mass) + PPOW(p_Mass, a[41]))) - 1.0;

    // calculate the Hurley radius coefficient gamma, per Hurley et al. 2000, eq 23
    double gamma;
    if (p_Mass <= 1.0) {
        // BSE Fortran code has abs()
        gamma = a[76] + (a[77] * PPOW(std::abs(p_Mass - a[78]), a[79]));
    }
    else if (p_Mass <= a[75]) {
        const double bGamma = GLOBALS->HurleyGammaConstants(HURLEY_G::B_GAMMA);
        gamma = bGamma + (a[80] - bGamma) * PPOW((p_Mass - 1.0) / (a[75] - 1.0), a[81]);
    }
    else if (p_Mass <= (a[75] + 0.1)) {
        // see discussion just prior to eq 23 - the end point is wrong in the arXiv version of
        // Hurley et al. 2000 (should be 0.1, not 1.0) - confirmed in BSE Fortran code
        const double cGamma = GLOBALS->HurleyGammaConstants(HURLEY_G::C_GAMMA);
        gamma = cGamma - (10.0 * (p_Mass - a[75]) * cGamma);
    }
    else {
        // see discussion just prior to eq 23
        // confirmed in BSE Fortran code
        gamma = 0.0;
    }
    
    gamma = std::max(0.0, gamma); // see discussion following eq 23 - confirmed in BSE Fortran code


    // calculate radius per Hurley et al. 2000
    const double tBGB  = utils::CalculateLifetimeToBGB_Hurley2000(p_Mass);
    const double tMS   = CalculatePhaseLifetime_Hurley2000(p_Mass, p_tBGB);
    const double rTMS  = CalculateRadiusAtPhaseEnd_Hurley2000(p_Mass, p_RZAMS);
    const double mu    = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, p_aN[7])), (p_aN[8] + (p_aN[9] / PPOW(p_Mass, p_aN[10]))))))); // eq 7
    const double tHook = mu * p_tBGB;                                                                       // just after eq 5
    const double time  = tMS * p_Tau;
    const double tau1  = std::min(1.0, (time / tHook));                                                     // eq 14
    const double tau2  = std::max(0.0, std::min(1.0, (time - ((1.0 - 0.01) * tHook)) / (0.01 * tHook)));    // eq 15, epsilon= 0.01
    const double tau_3  = p_Tau * p_Tau * p_Tau;                                                            // pow() is slow - use multiplication where it makes sense
    const double tau_10 = p_Tau < FLOAT_TOLERANCE_ABSOLUTE ? 0.0: p_Tau * tau_3 * tau_3 * tau_3;            // tolerance comparison to avoid underflow
    const double tau_40 = tau_10 < FLOAT_TOLERANCE_ABSOLUTE ? 0.0: tau10 * tau10 * tau10 * tau10;           // tolerance comparison to avoid underflow
    const double tau1_3 = tau1 * tau1 * tau1;
    const double tau2_3 = tau2 * tau2 * tau2;

    // eq 13
    const double logRMS_RZAMS = alphaR * p_Tau + betaR * tau_10 + gamma * tau_40 + (std::log10(rTMS / p_RZAMS) - alphaR - betaR - gamma) * tau_3 - deltaR * (tau1_3 - tau2_3);

    return p_RZAMS * PPOW(10.0, logRMS_RZAMS);
}

    







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

    const double radiusTAMS        = std::min(HG::CalculateRadiusOnPhase_Static(p_Mass, p_Tau, p_RZAMS, p_bCoefficients), p_Radius);
    const double radiusAtHookStart = CalculateRadiusOnPhase_Hurley(mass, 0.99, p_RZAMS); 
        
    double radius = (radiusAtHookStart * (1.0 - p_Tau) + radiusTAMS * (p_Tau - 0.99)) / 0.01;   // linear interpolation

    // if the star has been stripped below its initial core mass we need to adjust the radius
  



    double heliumAbundanceSurface = m_HeliumAbundanceSurface;
        if (p_Mass < m_InitialMainSequenceCoreMass)
            // By tracing the helium profile in the star, calculate how the surface helium abundance changes if mass drops below the initial core mass
            heliumAbundanceSurface = m_HeliumAbundanceOutsideCore + (p_Mass - m_MainSequenceCoreMass) * (m_HeliumAbundanceSurface - m_HeliumAbundanceOutsideCore) / (m_InitialMainSequenceCoreMass - m_MainSequenceCoreMass);
        
        // Factor that scales radius based on surface helium abundance
        double surfaceAbundanceFactor = (utils::Compare(m_HeliumAbundanceCore, m_InitialHeliumAbundance) != 0) ? (heliumAbundanceSurface - m_InitialHeliumAbundance) / (m_HeliumAbundanceCore - m_InitialHeliumAbundance) : 0.0;
        
        radius += (p_RZAMS - radius) * surfaceAbundanceFactor;






    return radius;
}






/*
 * CalculateRadius
 *
 * @brief
 * Calculate radius on the Main Sequence.
 *
 * If:
 *    (a) the Brcek MS core mass prescription was specified by the user,
 *    (b) the star is in the Brcek regime (mass >= BRCEK_LOWER_MASS_LIMIT), and
 *    (c) the star is in the MS hook (tau > 0.99),
 * 
 * calculate the MS radius per Brcek et al. 2025, otherwise uses the method appropriate
 * for the evolution method and as specified by program options.
 * 
 * 
 * double CalculateRadius() const
 *
 * @return                                      MS radius (Rsol)
 */
double MainSequence::CalculateRadius() const {
        
    double radius;

    if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK &&      // Brcek prescription specified?
        mZAMS >= BRCEK_LOWER_MASS_LIMIT) {                                                      // Yes - in Brcek regime?    
                                                                                                // Yes - use Brcek
        radius = CalculateRadius_Brcek2025(Mass(), Tau(), Time(), MZAMS(), LZAMS(), CoreMass(), HeAbundanceCore());
    }
    else {                                                                                      // No - use default method
        radius = CalculateRadius(Mass(), Time(), LZAMS());
    }

    return radius;
}



/*
 * CalculateRadiusOnPhase
 *
 * @brief
 * Calculate radius on the Main Sequence
 *
 *
 * double CalculateRadiusOnPhase() const
 *
 * @return                                      Main Sequence radius of the star (Rsol)
 */
double MainSequence::CalculateRadiusOnPhase() const { 

    const double Z              = GLOBALS->Metallicity();                                                              // metallicity of the star
    const double mass           = Mass();                                                         // current mass of the star
    const double tau           = Tau();                                                         // time elapsed since ZAMS
    const double rZAMS  = RZAMS();

    const double currentRadius = Radius();

    double radius;

    // if the user specified the BRCEK MS core mass prescription, and the star is in the 
    // BRCEK MS core mass prescription regine, and is on the MS hook, calculate the MS radius
    // per Brcek et al., 2025


    if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK && MZAMS() >= BRCEK_LOWER_MASS_LIMIT  && // yes - in Brcek MS core mass prescription regime?
        Tau() > 0.99) {                                                                       // yes - on MS hook?
                                                                                            // yes
        radius = CalculateRadiusOnPhase_Brcek(Mass(), Tau(), RZAMS(), Radius());
    }
    else {
        Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

            EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
            EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
                radius = CalculateRadiusOnPhase_Hurley();
                break;
        
            default:                                                                        // unknown mode
                // the only way this can happen is if someone added an EVOLUTION_MODE
                // and it isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user chose a mode this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing mode or, if the missing mode is
                // superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
        }       
    }

    return radius;
}


/*
 * CalculateRadiusMStoHG
 *
 * @brief
 * Calculate radius on the transition from the Main Sequence to the Hertzsprung Gap.
 * 
 * Used when BRCEK core mass prescription is used.
 *
 * Core mass prescription from Shikauchi et al. (2024) cannot be used beyond the MS hook
 * (beyond Tau = 0.99), and this function smoothly connects the radius between the beginning
 * of the MS hook and the beginning of the Hertzsprung Gap.
 *
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function ostensibly relies on the value of the ZAMS radius of the star,
 * and should not be used if the ZAMS radius is not known.
 * 
 *
 * double CalculateRadiusMStoHG(const double p_Mass, const double p_Tau)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @return                                      MS radius, for 0.99 <= Tau <= 1.0 (Rsol)
 */
double MainSequence::CalculateRadiusMStoHG(const double p_Mass, const double p_Tau, const double p_RZAMS) const {
    HG *clone = HG::Clone(static_cast<HG&>(const_cast<MainSequence&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    // Select radius at TAMS from the HG clone or current radius (whichever is smaller), relevant for stars that were significantly
    // stripped as this prevents radius expansion during the hook, and delays possible mass transfer to the start of HG
    double radiusTAMS = std::min(clone->Radius(), m_Radius);                                                                    // Get radius from clone (with updated Mass0)
    delete clone; clone = nullptr;                                                      // Return the memory allocated for the clone
    
    double radiusAtHookStart = CalculateRadiusOnPhase(p_Mass, 0.99, p_RZAMS);           // Hook starts at Tau = 0.99  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    return (radiusAtHookStart * (1.0 - p_Tau) + radiusTAMS * (p_Tau - 0.99)) / 0.01;    // Linear interpolation
}


/*
 * CalculateConvectiveEnvelopeRadialExtentHurley2002
 *
 * Calculate the radial extent of the star's convective envelope (if it has one),
 * per Hurley et al. 2002, sec. 2.3, particularly subsec. 2.3.1, eqs 36-38
 *
 * Uses radius of a 0.35 solar mass star at ZAMS rather than at fractional age Tau,
 * but such low-mass stars only grow by a maximum factor of 1.5
 * [just above Eq. (10) in Hurley, Pols, Tout (2000)], so this is a reasonable approximation
 * 
 *
 * double CalculateConvectiveEnvelopeRadialExtentHurley2002(const double p_Mass, const double p_Radius) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Msol)
 * @return                                      Radial extent of the star's convective envelope (Rsol)
 */
double MainSequence::CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Mass, const double p_Radius) const {

    double radiusEnvelope;

         if (p_Mass >= 1.25) radiusEnvelope = 0.0;
    else if (p_Mass > 0.35 ) radiusEnvelope = CalculateRadiusAtZAMS_Tout1996(0.35) * std::sqrt((1.25 - p_Mass) / 0.9);
    else                     radiusEnvelope = p_Radius;

    return radiusEnvelope * std::sqrt(std::sqrt(1.0 - m_Tau));
}


/*
 * CalculateConvectiveCoreRadius
 *
 * Calculate the radial extent of the star's convective core (if it has one)
 *
 * Uses preliminary fit from Minori Shikauchi @ ZAMS, then a smooth interpolation to the HG
 *
 *
 * double CalculateConvectiveCoreRadius()
 *
 * @return                                      Radial extent of the star's convective core in Rsol
 */
double MainSequence::CalculateConvectiveCoreRadius(const double p_Mass, const double p_Tau) const {

    double radius = 0.0;                                                                            // Default return value

    if (p_Mass >= 1.25) {                                                                           // Star has convective core?
                                                                                                    // Yes
        double ZAMSconvectiveCoreRadius = p_Mass * (0.06 + 0.05 * std::exp(-p_Mass / 61.57));
    
        // We need TAMS core radius, which is just the core radius at the start of the HG phase.
        // Since we are on the main sequence here, we can clone this object as an HG object and,
        // as long as it is initialised (to correctly set Tau to 0.0 on the HG phase), we can 
        // query the cloned object for its core radius.
        //
        // The clone should not evolve, and so should not log anything, but to be sure the
        // clone does not participate in logging, we set its persistence to EPHEMERAL.

        std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::HG, OBJECT_PERSISTENCE::EPHEMERAL);

        radius = ZAMSconvectiveCoreRadius - p_Tau * (ZAMSconvectiveCoreRadius - clone->CalculateRemnantRadius());
    }

    return radius;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCoreMass_Shikauchi2024
 *
 * @brief
 * Calculate the mass of the convective core, per Shikauchi, Hirai, Mandel, 2024
 * (see https://arxiv.org/abs/2409.00460)
 * 
 * Core mass shrinks to 60% of initial value over the course of the MS
 *
 *
 * double CalculateCoreMass_Shikauchi2024(const double p_Tau) const
 * 
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 *
 * @return                                      Mass of convective core (Msol)
 */
double MainSequence::CalculateCoreMass_Shikauchi2024(const double p_Tau) const {
    const double finalConvectiveCoreMass   = CalculateTAMSCoreMass();   // core mass at TAMS
    const double initialConvectiveCoreMass = finalConvectiveCoreMass / 0.6;

    return (initialConvectiveCoreMass - p_Tau * (initialConvectiveCoreMass - finalConvectiveCoreMass));
}


/*
 * CalculateConvectiveEnvelopeMass_Hurley2000
 *
 * @brief
 * Calculate the convective envelope mass for Main Sequence, MS, stars.
 *
 * Based on Hurley et al. 2000, section 7.2, after eq 111
 *
 *
 * Dbl_DblT CalculateConvectiveEnvelopeMass_Hurley2000(const double p_Mass, const double p_Tau) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      Tuple containing:
 *                                                   MS convective envelope mass (Msol)
 *                                                   Maximum MS (ZAMS) convective envelope mass (Msol)
 */
Dbl_DblT MainSequence::CalculateConvectiveEnvelopeMass_Hurley2000(const double p_Mass, const double p_Tau) const {

    double mEnv;
    double mEnv0;

    if (p_Mass <= 1.25) {   // no convective envelope above 1.25 Msol
        mEnv0 = p_Mass > 0.35 ? 0.35 * (1.25 - p_Mass) * (1.25 - p_Mass) / 0.81 : p_Mass;
        mEnv  = mEnv0 * std::sqrt(std::sqrt(1.0 - p_Tau));
    }
    else {
        mEnv0 = 0.0;
        mEnv  = 0.0;
    }

    return std::make_tuple(mEnv, mEnv0);
}


// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< use below to calculate abundances <<<<<<<<<<<<<<<<<<<<<<
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
                                            const double p_CNOprocessedCoreMass,   /// CNOprocessedCoreMass

                                            const double p_MSCoreMass,

                                            const double p_InitialHeAbundance,      // initial??????  better name <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                            const double p_HeAbundanceCore,
                                            const double p_HeAbundanceOutsideCore,
                                            const double p_HeAbundanceSurface,

                                            ) const { 
// check dt = 0 - do nothing

    DBL_VECTOR aCoeffs = GLOBALS->ShikauskiAlphaCoefficients();                                                         // get Shikauchi alpha coefficients from GLOBALS
    DBL_VECTOR fCoeffs = GLOBALS->ShikauskifMixCoefficients();                                                          // get Shikauchi fMix coefficients from GLOBALS
    DBL_VECTOR dCoeffs = SHIKAUCHI_DELTA_COEFFICIENTS;                                                                  // get Shikauchi delta coefficients from constants.h


    auto fmix = [&](double mass) { return fCoeffs[0] + fCoeffs[1] * std::exp(-mass / fCoeffs[2]); }                     // Shikauchi et al. 2024, eq A3


    const double updatedMass = p_Mass + (p_dMdt * p_Dt * MYR_TO_YEAR);                                                  // mass adjusted for loss/gain

    double delta;                                                                                                       // mixing core inertia
    if (p_dMdt <= 0.0) {                                                                                                // mass loss?
                                                                                                                        // yes
        const double g  = dCoeffs[1] * p_MSCoreMass + dCoeffs[2];                                                       // ibid., eq A7
        const double Yc = -dCoeffs[0] * (p_HeAbundanceCore - p_InitialHeAbundance) / (1.0 - p_InitialHeAbundance - GLOBALS->Metallicity());
        delta = std::min(PPOW(10.0, Yc + g), 1.0);                                                                      // ibid., eq A6
    }
    else {                                                                                                              // no - mass gain
        delta = PPOW(2.0, -(p_HeAbundanceCore - p_InitialHeAbundance) / (1.0 - p_InitialHeAbundance - GLOBALS->Metallicity()));  // updated prescription for mass gain
    }

    const double deltaCoreMassML      = p_MSCoreMass * delta * (updatedMass * fmix(updatedMass) / (p_Mass * fmix(p_Mass)) - 1.0); // core mass delta due to mass loss/gain
    const double alpha                = PPOW(10.0, std::max(-2.0, aCoeffs[1] * p_MSCoreMass + aCoeffs[2])) + aCoeffs[0];// ibid., eq (A2)
    const double deltaYc              = CalculateLuminosityOnPhase() / (Q_CNO * p_MSCoreMass) * p_Dt;                   // ibid., eq (12) (central helium fraction delta)
    const double deltaCoreMassNatural = -alpha / (1 - alpha * p_HeAbundanceCore) * deltaYc * p_MSCoreMass;              // ibid., eq (4)  (core mass delta due to natural decay)
    const double deltaCoreMass        = deltaCoreMassNatural + deltaCoreMassML;                                         // total core mass delta
    
    const double mixingCoreMass   = std::min(p_MSCoreMass + deltaCoreMass, BRCEK_CORE_MASS_TO_MASS_RATIO_LIMIT * updatedMass); // mixing core mass, must be <= total mass
    double centralHeFraction      = std::min(p_HeAbundanceCore + deltaYc, 1.0 - GLOBALS->Metallicity());                         // central helium fraction, clamped to a maximum of (1-Z)
    double HeAbundanceOutsideCore = p_HeAbundanceOutsideCore;                                                           // He abundance outside the core
    double CNOprocessedCoreMass   = p_CNOprocessedCoreMass;                                                             // CNO processed core mass

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
            const double HAbundanceSurface  = 1.0 - GLOBALS->Metallicity() - HeAbundanceSurface;

            CNOprocessedCoreMass = updatedMass;
        }

        HeAbundanceOutsideCore = centralHeFraction;
    }
    
    // also return HeAbundanceOutsideCore, CNOprocessedCoreMass, HeAbundanceSurface, and HAbundanceSurface <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    return std::make_tuple(mixingCoreMass, std::min(centralHeFraction, 1.0 - GLOBALS->Metallicity()));
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
 *     const double p_Mass,
 *     const double p_Luminosity,
 *     const double p_CoreMass,
 *     const double p_HeAbundanceCore,
 *     const double p_dt,
 *     const double p_dMdt,
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core of the star 
 * @param       p_dt                            Time step (Myr)
 * @param       p_dMdt                          Mass loss rate (Msol yr^-1)
 * @return                                      MS convective core mass (Msol)
 */
COMPAS_PURE double MainSequence::CalculateCoreMass_Brcek2025(
    const double p_Mass,
    const double p_Luminosity,
    const double p_CoreMass,
    const double p_HeAbundanceCore,
    const double p_dt,
    const double p_dMdt
) const {
     
    double coreMass = p_CoreMass;                                                                                       // default is no change

    if (p_dt > 0.0) {                                                                                                   // only update core mass if timestep > 0

        DBL_VECTOR aCoeffs = GLOBALS->ShikauskiAlphaCoefficients();                                                     // get Shikauchi alpha coefficients from GLOBALS
        DBL_VECTOR fCoeffs = GLOBALS->ShikauskifMixCoefficients();                                                      // get Shikauchi fMix coefficients from GLOBALS
        DBL_VECTOR dCoeffs = SHIKAUCHI_DELTA_COEFFICIENTS;                                                              // get Shikauchi delta coefficients from constants.h

        auto fmix = [&](double mass) { return fCoeffs[0] + fCoeffs[1] * std::exp(-mass / fCoeffs[2]); }                 // Shikauchi et al. 2024, eq A3

        double delta;                                                                                                   // mixing core inertia
        if (p_dMdt <= 0.0) {                                                                                            // mass loss?
                                                                                                                        // yes
            const double g  = dCoeffs[1] * p_CoreMass + dCoeffs[2];                                                     // ibid., eq A7
            const double Yc = -dCoeffs[0] * (p_HeAbundanceCore - GLOBALS->ZAMSHeAbundance()) / (1.0 - GLOBALS->ZAMSHeAbundance() - GLOBALS->Metallicity());
            delta = std::min(PPOW(10.0, Yc + g), 1.0);                                                                  // ibid., eq A6
        }
        else {                                                                                                          // no - mass gain
            delta = PPOW(2.0, -(p_HeAbundanceCore - GLOBALS->ZAMSHeAbundance()) / (1.0 - GLOBALS->ZAMSHeAbundance() - GLOBALS->Metallicity()));
        }

        const double mAdjusted    = p_Mass + (p_dMdt * p_Dt * MYR_TO_YEAR);                                             // mass adjusted for loss/gain
        const double deltaMcML    = p_CoreMass * delta * (mAdjusted * fmix(mAdjusted) / (p_Mass * fmix(p_Mass)) - 1.0); // core mass delta due to mass loss/gain
        const double alpha        = PPOW(10.0, std::max(-2.0, aCoeffs[1] * p_CoreMass + aCoeffs[2])) + aCoeffs[0];      // ibid., eq A2
        const double deltaYc      = p_Luminosity / (Q_CNO * p_CoreMass) * p_Dt;                                         // ibid., eq 12 (central helium fraction delta)
        const double deltaMcDecay = -alpha / (1 - alpha * p_HeAbundanceCore) * deltaYc * p_CoreMass;                    // ibid., eq 4  (core mass delta due to natural decay)
        const double deltaMc      = deltaMcDecay + deltaMcML;                                                           // total core mass delta
    
        coreMass = std::min(p_CoreMass + deltaMc, BRCEK_CORE_MASS_TO_MASS_RATIO_LIMIT * mAdjusted);                     // mixing core mass, clamped to total mass
    }

    return coreMass;
}


/*
 * CalculateCoreMass   <<<<<<<<<<<<<<<<<<<<<<<<<<<< needs constituent version <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 *
 * @brief
 * Calculate the convective core mass of a main sequence star, based on the
 * MS core mass prescription specified by the user.
 * 
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function ostensibly relies on the value of the ZAMS mass of the star,
 * and should not be used if the ZAMS mass is not known.
 * 
 *
 * double CalculateCoreMass(
 *     const double p_Mass,
 *     const double p_Luminosity,
 *     const double p_MZAMS,
 *     const double p_CoreMass,
 *     const double p_HeAbundanceCore,
 *     const double p_dt, 
 *     const double p_dMdt
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core of the star 
 * @param       p_dt                            Time step (Myr)
 * @param       p_dMdt                          Mass loss rate (Msol yr^-1) 
 */
COMPAS_PURE double MainSequence::CalculateCoreMass(
    const double p_Mass,
    const double p_Luminosity,
    const double p_MZAMS,
    const double p_CoreMass,
    const double p_HeAbundanceCore,
    const double p_dt, 
    const double p_dMdt
) const {

    double coreMass = p_CoreMass;                                               // default is unchanged

    switch (OPTIONS->MainSequenceCoreMassPrescription()) {                      // which MS core mass prescription?
        
        case MS_CORE_MASS_PRESCRIPTION::BRCEK:                                  // BRCEK
            // calculate MS core mass per Brcek et al. 2025, following Shikauchi et al. 2024
            // account for rejuvenation if core grows
            if (p_MZAMS >= BRCEK_LOWER_MASS_LIMIT) {                            // in BRCEK regime?
                                                                                // yes
                if (p_HeAbundanceCore < (1.0 - GLOBALS->Metallicity())) {       // in MS hook?
                                                                                // no - proceed
                    coreMass = CalculateCoreMass_Brcek2025(p_Mass, p_Luminosity, p_CoreMass, p_HeAbundanceCore, p_dt, p_dMdt);
                }
            }
            break;

        case MS_CORE_MASS_PRESCRIPTION::HURLEY:                                 // HURLEY
            coreMass = 0.0;                                                     // no MS core in Hurley et al. 2000
            break;
        
        case MS_CORE_MASS_PRESCRIPTION::MANDEL:                                 // MANDEL
            // calculate MS core mass per Brcek et al. 2025, following Shikauchi et al. 2024
            // account for rejuvenation if core grows
            if (p_MZAMS >= BRCEK_LOWER_MASS_LIMIT) {                            // in BRCEK regime?
                                                                                // yes
                if (p_HeAbundanceCore < (1.0 - GLOBALS->Metallicity())) {       // in MS hook?
                                                                                // no - proceed
                    coreMass = CalculateCoreMass_Brcek2025(p_Mass, p_Luminosity, p_CoreMass, p_HeAbundanceCore, p_dt, p_dMdt);
                }
            }
            break;
        
        default:                                                                // unknown prescription
            // the only way this can happen is if someone added a MS_CORE_MASS_PRESCRIPTION and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MS_CORE_MASS_PRESCRIPTION);              // throw error
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
double MainSequence::CalculateZetaEquilibrium() {
    double deltaMass           = -m_Mass / 1.0E5;
    double radiusAfterMassGain = CalculateRadiusOnMassChange(deltaMass);
    double zetaEquilibrium     = (radiusAfterMassGain - m_Radius) / deltaMass * m_Mass / m_Radius;      // dlnR / dlnM

    return zetaEquilibrium;
}
    






///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              MISCELLANEOUS FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateHeAbundanceCoreOnPhase
 *
 * @brief
 * Calculate the helium abundance in the core of the star, given tau
 * 
 * Simple linear model from the initial helium abundance to the maximum helium
 * abundance (assuming that all hydrogen is converted to helium).
 * 
 * 
 
 * 
 * Should be updated to match detailed models.
 * 
 *
 * double CalculateHeAbundanceCoreOnPhase(const double p_Tau)
 * 
 * @param   [IN]    p_Tau                       Fraction of main sequence lifetime
 * @return                                      Helium abundance in the core of the star
 */
double MainSequence::CalculateHeAbundanceCoreOnPhase(const double p_Tau, const double p_HeAbundanceCore) const {
    
    // If BRCEK core mass prescription is used, core helium abundance is calculated with the core mass   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FIX THIS
    return (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK && m_MZAMS >= BRCEK_LOWER_MASS_LIMIT)
            ? p_HeAbundanceCore
            : ((1.0 - GLOBALS->Metallicity() - GLOBALS->ZAMSHeAbundance()) * p_Tau) + GLOBALS->ZAMSHeAbundance();
}


/*
 * Calculate the hydrogen abundance in the core of the star
 * 
 * Currently just a simple linear model. Assumes that hydrogen in the core of 
 * the star is burned to helium at a constant rate throughout the lifetime. 
 * 
 * 
 * Calculate the hydrogen abundance in the core of the star, given the phase-relative age
 * of the star.  Currently just a simple linear model that assumes that hydrogen in the
 * core of the star is burned to helium at a constant rate throughout the star's lifetime.
 * 
 * Should be updated to match detailed models.
 *
 * double CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau, const double p_InitialHAbundance)
 * 
 * @param   [IN]    p_Tau                       Fraction of main sequence lifetime
 * @return                                      Hydrogen abundance in the core (X_c)
 */
double MainSequence::CalculateHAbundanceCoreOnPhase(const double p_Tau, const double p_InitialHAbundance) const {
    
    // If BRCEK core mass prescription is used, core helium abundance is calculated with the core mass   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FIX THIS
    return OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK && m_MZAMS >= BRCEK_LOWER_MASS_LIMIT
            ? 1.0 - m_HeliumAbundanceCore - m_Metallicity
            : m_InitialHydrogenAbundance * (1.0 - p_Tau);
}











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

    m_Mass                 = p_Mass;
    m_Mass0                = m_Mass;
    m_MainSequenceCoreMass = 0.0;
    
    double initialHydrogenFraction = m_InitialHydrogenAbundance;
    
    CalculateTimescales();      /// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    CalculateGBparams();        /// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
    m_Tau = (initialHydrogenFraction - p_HydrogenMass / m_Mass) / initialHydrogenFraction;      // assumes uniformly mixed merger product and a uniform rate of H fusion on main sequence
    
    m_Age = m_Tau * timescales(tMS);
    
    m_HeliumAbundanceCore   = 1.0 - m_Metallicity - p_HydrogenMass / p_Mass;
    m_HydrogenAbundanceCore = 1.0 - m_Metallicity - m_HeliumAbundanceCore;
    
    m_HeliumAbundanceSurface   = m_HeliumAbundanceCore;                                         // abundances are the same throughout the star, assuming uniform mixing after merger
    m_HydrogenAbundanceSurface = m_HydrogenAbundanceCore;
    
    if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK && m_MZAMS >= BRCEK_LOWER_MASS_LIMIT) {
        m_InitialMainSequenceCoreMass = CalculateCNOprocessedCoreMass_Brcek2025(p_Mass, m_HeliumAbundanceCore);           // update initial mixing core mass  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        m_MainSequenceCoreMass        = m_InitialMainSequenceCoreMass;                                                 // update core mass
    }
    
    EvolveOneTimestep(0.0, 0.0, 0.0, true);
}
















////////////////////////////////////////////////// MainSequence_Constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<







/*
 * CalculateCoreMass
 *
 * @brief
 * Calculate the core mass on the main sequence star
 * /////// that loses mass through winds or Case A mass transfer
 * When BRCEK core prescription is used, also update the core helium abundance and effective age  ///////
 *
 *
 * void UpdateMainSequenceCoreMass(const double p_Dt, const double p_MassLossRate)
 * double CalculateCoreMass(const double p_Dt, const double p_MassLossRate)
 *
 * @param   [IN]      p_Dt                      Current timestep in Myr
 * @param   [IN]      p_MassLossRate            Mass loss rate either from stellar winds or mass transfer in Msol yr-1
 */
double MainSequence_Constituent::CalculateCoreMass(
    const double p_Dt, 
    const double p_MZAMS,
    const double p_CoreMass,
    const double p_HAbundanceInit,
    const double p_HeAbundanceInit,
    const double p_HeAbundanceCore,
    const double p_dMdt,
    const double p_tMS,
    const double p_MassLossRate) {

    double mainSequenceCoreMass = m_MainSequenceCoreMass;                                                                               // default is no change
    double HeAbundanceCore  = p_HeAbundanceCore;                                                                                // default is no change
    double age                  = m_Age;                                                                                                // default is no change

    const double oneMinusZ = 1.0 - GLOBALS->Metallicity();
    double coreMass = p_CoreMass;  // <<<<<<<<<<<<<<<<<<<<<<<<<< need assignment?

    switch (OPTIONS->MainSequenceCoreMassPrescription()) {              // which MS core mass prescription?
        
        case MS_CORE_MASS_PRESCRIPTION::BRCEK:                          // BRCEK
            // calculate MS core mass following Shikauchi et al. 2024
            // account for rejuvenation if core grows
            if (p_MZAMS >= BRCEK_LOWER_MASS_LIMIT) {                                    // in BRCEK regime?
                                                                                        // yes
                // Only proceed with calculations if star is not in MS hook (Yc < 1-Z) and time step is not zero
                if (p_HeAbundanceCore < oneMinusZ && p_Dt > 0.0) {  /////  the p_Dt check shouldn't be necessary  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< THIS SHOULDN'T BE NECESSARY - CHECK LATER <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR
                    // Update the core mass and central helium fraction only if the mass loss rate argument is equal
                    // to the total mass loss rate (i.e. total mass loss rate was updated, this prevents the calculation
                    // in SSE if it was executed as part of BSE for the same time step)


                    if (utils::Compare(p_MassLossRate, m_TotalMassLossRate) == 0) {   ////// BSE & SSE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        // Calculate and update the core mass and central helium fraction
                        std::tie(coreMass, HeAbundanceCore) = CalculateCoreMass_Brcek(p_Dt, p_MassLossRate);
                        // Update effective age here only if core hydrogen was exhausted
                        age = HeAbundanceCore == oneMinusZ ? 0.99 * p_tMS : age;
                    }
                    // Update effective age only when stars are aged in SSE (when p_MassLossRate = -Mdot)
                    if (p_MassLossRate == -p_dMdt)
                        // Update the effective age based on central helium fraction
                        age = (p_HeAbundanceCore - p_HeAbundanceInit) / p_HAbundanceInit * 0.99 * p_tMS;  //// SSE <<<<<<<<<<<<<<<<<<<<<<
                }
            }
            else {                                                                                                                      // no - MANDEL prescription used
                // Only applied to donors as part of binary evolution, not applied to SSE
                if (p_MassLossRate < 0.0 && p_MassLossRate != -p_dMdt)
                    coreMass = std::max(p_CoreMass, CalculateTauOnPhase() * CalculateCoreMassAtTAMS_Hurely2000());
            }
            break;

            case MS_CORE_MASS_PRESCRIPTION::HURLEY:                                             // HURLEY
            coreMass = 0.0;                                                                     // no MS core in Hurley et al. 2000
            break;
        
        case MS_CORE_MASS_PRESCRIPTION::MANDEL:                                                 // MANDEL
            // Calculate the minimum core mass of a main sequence star that loses mass through 
            // Case A mass transfer as the core mass of a TAMS star, scaled by the fractional age.
            // Only applied to donors as part of binary evolution, not applied to SSE
            if ((OPTIONS->RetainCoreMassDuringCaseAMassTransfer()) && (p_MassLossRate < 0.0) && (utils::Compare(p_MassLossRate, -m_Mdot) != 0))
                coreMass = std::max(m_MainSequenceCoreMass, CalculateTauOnPhase() * CalculateCoreMassAtTAMS_Hurely2000());
            break;

        default: 
            // the only way this can happen is if someone added a MS_CORE_MASS_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MS_CORE_MASS_PRESCRIPTION);                                                                      // throw error
    }

    m_MainSequenceCoreMass = mainSequenceCoreMass;                                                                                      // update core mass
    m_HeliumAbundanceCore  = heliumAbundanceCore;                                                                                       // update core helium abundance
    m_Age                  = age;                                                                                                       // update age
}







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
    std::vector<GE_QCRIT_TABLE*> qCritTables = { &QCRIT_GE_LOW_Z, &QCRIT_GE_HIGH_Z };
    double qCritPerZ[qCritTables.size()];

    for (int tblIdx = 0; tblIdx < qCritTables.size(); tblIdx++) {       // for each qCrit table

        const DBL_VECTOR massVec = std::get<0>(*(qCritTables[tblIdx])); // vector of masses from qCrit table

        INT_VECTOR indices = utils::BinarySearch(massVec, p_Mass);      // find mass bin for p_Mass
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
        const std::vector<std::tuple<DBL_VECTOR, std::vector<DBL_VECTOR>>> radiiQCritsVec = std::get<1>(*(qCritTables[tblIdx]));

        // (log) radii vectors for mass bin bounds
        const DBL_VECTOR radiiVecLow = std::get<0>(radiiQCritsVec[lowMIdx]);
        const DBL_VECTOR radiiVecUpp = std::get<0>(radiiQCritsVec[uppMIdx]);

        // qCrits vectors for mass bin bounds
        const DBL_VECTOR qCritsVecLow = std::get<1>(radiiQCritsVec[lowMIdx]);
        const DBL_VECTOR qCritsVecUpp = std::get<1>(radiiQCritsVec[uppMIdx]);
 
        std::vector<double> qCritVecUppELowM;
        std::vector<double> qCritVecUppEUppM;
        std::vector<double> qCritVecLowELowM;
        std::vector<double> qCritVecLowEUppM;
        
        // set the appropriate qCrit vector based on QCRIT_PRESCRIPTION and MT efficiency
        switch (OPTIONS->QCritPrescription()) {                         // which QCRIT_PRESCRIPTION?

            QCRIT_PRESCRIPTION::GE:                                     // GE
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

            QCRIT_PRESCRIPTION::GE_IC:                                  // GE_IC
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
        const double rLowMLowR = radiiVecLower[lowRLowMIdx];
        const double rLowMUppR = radiiVecLower[uppRLowMIdx];
        const double rUppMlowR = radiiVecUpper[lowRUppMIdx];
        const double rUppMUppR = radiiVecUpper[uppRUppMIdx];

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
             if (logRadius < rLowMLowR) ? qCritUppELowM = qUppLowLow;
        else if (logRadius > rLowMUppR) ? qCritUppELowM = qUppLowUpp;
        else                              qCritUppELowM = qUppLowLow + rSlope1 * (qUppLowUpp - qUppLowLow);

        // qCrit upper efficiency, upper mass
        double qCritUppEUppM;
             if (logRadius < rUppMLowR) qCritUppEUppM = qUppUppLow;
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
    constexpr double logZhi = LOG10_ZSOL_HURLEY;                        // log10(0.02)  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FIX
    
    return qCritPerZ[1] + (m_Log10Metallicity - logZhi) * (qCritPerZ[1] - qCritPerZ[0]) / (logZhi - logZlo);  // <<<<<<<<<<<<<<<<<<<<<<<<< FIX m_Log10Metallicity
}
