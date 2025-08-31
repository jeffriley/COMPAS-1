#include "MainSequence.h"
#include "MS_gt_07.h"
#include "HG.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateHeliumAbundanceCoreOnPhase
 *
 * @brief
 * Calculate the helium abundance in the core of the star
 * 
 * Currently just a simple linear model from the initial helium abundance to 
 * the maximum helium abundance (assuming that all hydrogen is converted to
 * helium). 
 * 
 * When tau = 0, heliumAbundanceCore = m_InitialHeliumAbundance
 * When tau = 1, heliumAbundanceCore = heliumAbundanceCoreMax = 1.0 - m_Metallicity
 * 
 * Should be updated to match detailed models.
 *
 * double CalculateHeliumAbundanceCoreOnPhase(const double p_Tau)
 * 
 * @param   [IN]    p_Tau                       Fraction of main sequence lifetime
 * @return                                      Helium abundance in the core
 */
double MainSequence::CalculateHeliumAbundanceCoreOnPhase(const double p_Tau) const {
    
    // If BRCEK core mass prescription is used, core helium abundance is calculated with the core mass
    return BRCEK_MS_CORE_MASS_REGIME(m_MZAMS) ? m_HeliumAbundanceCore : ((1.0 - m_Metallicity - m_InitialHeliumAbundance) * p_Tau) + m_InitialHeliumAbundance;
}


/*
 * Calculate the hydrogen abundance in the core of the star
 * 
 * Currently just a simple linear model. Assumes that hydrogen in the core of 
 * the star is burned to helium at a constant rate throughout the lifetime. 
 * 
 * Should be updated to match detailed models.
 *
 * double CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau)
 * 
 * @param   [IN]    p_Tau                       Fraction of main sequence lifetime
 * @return                                      Hydrogen abundance in the core (X_c)
 */
double MainSequence::CalculateHydrogenAbundanceCoreOnPhase(const double p_Tau) const {
    
    // If BRCEK core mass prescription is used, core helium abundance is calculated with the core mass
    return BRCEK_MS_CORE_MASS_REGIME(m_MZAMS) ? 1.0 - m_HeliumAbundanceCore - m_Metallicity : m_InitialHydrogenAbundance * (1.0 - p_Tau);
}


/*
 * Calculate timescales in units of Myr
 *
 * Timescales depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 *
 * void CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales)
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN/OUT]    p_Timescales            Timescales
 */
void MainSequence::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
    timescales(tBGB) = CalculateLifetimeToBGB(p_Mass);
    timescales(tMS)  = CalculateLifetimeOnPhase(p_Mass, timescales(tBGB));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LUMINOSITY CALCULATIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the luminosity perturbation delta_L
 *
 * Hurley et al. 2000, eq 16
 *
 *
 * double CalculateDeltaL(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity perturbation (delta_L in Hurley et al. 2000)
 */
double MainSequence::CalculateDeltaL(const double p_Mass) const {
#define an m_AnCoefficients                                              // for convenience and readability - undefined at end of function

    double deltaL;

    if (utils::Compare(p_Mass, massCutoffs(MHook)) <= 0) {              // per Hurley et al. 2000, eq 16
        deltaL = 0.0;                                                   // this really is supposed to be zero
    }
    else if (utils::Compare(p_Mass, an[33]) < 0) {
        double top    = p_Mass - massCutoffs(MHook);
        double bottom = an[33] - massCutoffs(MHook);
        deltaL        = m_LConstants[static_cast<int>(L_CONSTANTS::B_DELTA_L)] * PPOW((top / bottom), 0.4);
    }
    else {
        deltaL = std::min((an[34] / PPOW(p_Mass, an[35])), (an[36] / PPOW(p_Mass, an[37])));
    }

    return deltaL;

#undef an
}


/*
 * Calculate the luminosity beta coefficient
 *
 * Hurley et al. 2000, eq 20
 *
 *
 * double CalculateBetaL(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity beta coefficient (beta_L in Hurley et al. 2000)
 */
double MainSequence::CalculateBetaL(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    double betaL  = std::max(0.0, (an[54] - (an[55] * PPOW(p_Mass, an[56]))));
    if ((utils::Compare(p_Mass, an[57]) > 0) && (utils::Compare(betaL, 0.0) > 0)) {
        double bBetaL = m_LConstants[static_cast<int>(L_CONSTANTS::B_BETA_L)];

        betaL = std::max(0.0, (bBetaL - 10.0 * (p_Mass - an[57]) * bBetaL));
    }

    return betaL;

#undef a
}


/*
 * Calculate the luminosity alpha constant alpha_L
 *
 * Hurley et al. 2000, eqs 19a & 19b
 *
 *
 * double CalculateAlphaL(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity alpha constant (alpha_L in Hurley et al. 2000)
 */
double MainSequence::CalculateAlphaL(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    // You might find that functions you give Metallicity to as an argument don't actually need it -- metallicity dependence is in an/MFGB etc.
    // Also, if this is likely to be called in a loop, try to precompute it (only depends on initial values of mass/metallicity right?
    // Of course the only problem with this kind of thing is it makes it less flexible if you have to change one of those)

    double alphaL  = 0.0;

         if (utils::Compare(p_Mass, 0.5)   < 0) alphaL = an[49];
    else if (utils::Compare(p_Mass, 0.7)   < 0) alphaL = an[49] + (5.0 * (0.3 - an[49]) * (p_Mass - 0.5));
    else if (utils::Compare(p_Mass, an[52]) < 0) alphaL = 0.3 + ((an[50] - 0.3) * (p_Mass - 0.7) / (an[52] - 0.7));
    else if (utils::Compare(p_Mass, an[53]) < 0) alphaL = an[50] + ((an[51] - an[50]) * (p_Mass - an[52]) / (an[53] - an[52]));
    else if (utils::Compare(p_Mass, 2.0)   < 0) alphaL = an[51] + ((m_LConstants[static_cast<int>(L_CONSTANTS::B_ALPHA_L)] - an[51]) * (p_Mass - an[53]) / (2.0 - an[53]));
    else                                        alphaL = (an[45] + (an[46] * PPOW(p_Mass, an[48]))) / (PPOW(p_Mass, 0.4) + (an[47] * PPOW(p_Mass, 1.9)));

    return alphaL;

#undef a
}


/*
 * Calculate the exponent eta (for Hurley et al. 2000, eq 12)
 *
 * Hurley et al. 2000, eq 18
 *
 *
 * double CalculateEta(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      The exponent eta (for Hurley et al. 2000, eq 12)
 */
double MainSequence::CalculateEta(const double p_Mass) const {

    double eta = 10.0;

    if (utils::Compare(m_Metallicity, 0.0009) <= 0) {
        if (utils::Compare(p_Mass, 1.1) >= 0) {
            eta = 20.0;
        }
        else if (utils::Compare(p_Mass, 1.0) > 0) {
            eta = (100.0 * p_Mass) - 90.0;  // linear interpolation between end points
        }
    }

    return eta;
}


/*
 * Calculate the gamma coefficient
 *
 * Hurley et al. 2000, eq 23
 *
 *
 * double CalculateGamma(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      The exponent eta (for Hurley et al. 2000, eq 12)
 */
double MainSequence::CalculateGamma(const double p_Mass) const {
#define a m_AnCoefficients                                                      // for convenience and readability - undefined at end of function
#define B_GAMMA m_GammaConstants[static_cast<int>(GAMMA_CONSTANTS::B_GAMMA)]    // for convenience and readability - undefined at end of function
#define C_GAMMA m_GammaConstants[static_cast<int>(GAMMA_CONSTANTS::C_GAMMA)]    // for convenience and readability - undefined at end of function

    double gamma = 0.0;                                                                                                         // default return value

         if (utils::Compare(p_Mass,  1.0)          <= 0) gamma = a[76] + (a[77] * PPOW(std::abs(p_Mass - a[78]), a[79]));       // BSE Fortran code has abs()
    else if (utils::Compare(p_Mass,  a[75])        <= 0) gamma = B_GAMMA + (a[80] - B_GAMMA) * PPOW((p_Mass - 1.0) / (a[75] - 1.0), a[81]);
    else if (utils::Compare(p_Mass, (a[75] + 0.1)) <= 0) gamma = C_GAMMA - (10.0 * (p_Mass - a[75]) * C_GAMMA);                 // see discussion just prior to eq 23 - the end point is wrong in the arxiv version of Hurley et al. 2000 (should be 0.1, not 1.0) - confirmed in BSE Fortran code
    else                                                 gamma = 0.0;                                                           // see discussion just prior to eq 23 - confirmed in BSE Fortran code

    return std::max(0.0, gamma);                                                                                                // see discussion following eq 23 - confirmed in BSE Fortran code

#undef C_GAMMA
#undef B_GAMMA
#undef a
}


/*
 * CalculateZAMSLuminosity_Tout_1996
 *
 * @brief
 * Calculate ZAMS luminosity (in Lsol) per Tout et al. 1996, eq 1
 *
 *
 * double CalculateZAMSLuminosity_Tout_1996(const double p_MZAMS, const DBL_VECTOR& p_LuminosityCoefficients) const
 * 
 * @param       p_MZAMS                         Zero age main sequence mass (Msol)
 * @param       p_LuminosityCoefficients        Tout luminosity coefficients
 * @return                                      ZAMS luminosity (Lsol)
 */
double MainSequence::CalculateZAMSLuminosity_Tout_1996(const double p_MZAMS, const DBL_VECTOR& p_LuminosityCoefficients) const {
    #define lCoeffs(x) p_LuminosityCoefficients[static_cast<int>(TOUT_L_Coeff::x)] // for convenience and readability - undefined at end of function
        
        // calculate some powers of p_MZAMS - for performance and readability
        // this function is only called once per star, and at most twice per binary (but probably once), so not too onerous
        // pow() is slow - use multiplication where it makes sense
        const double M0_5 = std::sqrt(p_MZAMS);  // sqrt() is much faster than pow()
        const double M2   = p_MZAMS * p_MZAMS;
        const double M3   = p_MZAMS * M2;
        const double M5   = M2 * M3;
        const double M7   = M2 * M5;
        const double M8   = p_MZAMS * M7;
            
        double top = (lCoeffs(ALPHA) * (M5 * M0_5)) + (lCoeffs(BETA) * (M3 * M8));
            
        return top / (lCoeffs(GAMMA) + M3) + (lCoeffs(DELTA) * M5) + (lCoeffs(EPSILON) * M7) + (lCoeffs(ZETA) * M8) + (lCoeffs(ETA) * ( M8 * p_MZAMS * M0_5));
            
    #undef lCoeffs
    }
    
    
    /*
     * CalculateZAMSRadius_Tout_1996
     *
     * Calculate radius at ZAMS (in Rsol) per Tout et al. 1996, eq 2
     *
     *
     * double CalculateZAMSRadius_Tout_1996(const double p_MZAMS, const DBL_VECTOR& p_RadiusCoefficients)
     *
     * @param       p_MZAMS                         Zero age main sequence mass (Msol)
     * @param       p_RadiusCoefficients            Tout radius coefficients
     * @return                                      ZAMS radius (Rsol)
     */
    double MainSequence::CalculateZAMSRadius_Tout_1996(const double p_MZAMS, const DBL_VECTOR& p_RadiusCoefficients) const {
    #define rCoeffs(x) p_RadiusCoefficients[static_cast<int>(TOUT_R_Coeff::x)] // for convenience and readability - undefined at end of function
        
        // calculate some powers of p_MZAMS - for performance and readability
        // this function is only called once per star, and at most twice per binary (but probably once), so not too onerous
        // pow() is slow - use multiplication where it makes sense
        const double M0_5  = std::sqrt(p_MZAMS);
        const double M2    = p_MZAMS * p_MZAMS;
        const double M6    = M2 * M2 * M2;
        const double M6_5  = M6 * M0_5;
        const double M8    = M6 * M2;
        const double M11   = M8 * M2 * p_MZAMS;
        const double M19   = M11 * M8;
        const double M19_5 = M19 * M0_5;
        
        double top = (rCoeffs(THETA) * (M2 * M0_5)) + (rCoeffs(IOTA) * M6_5) + (rCoeffs(KAPPA) * M11) + (rCoeffs(LAMBDA) * M19) + (rCoeffs(MU) * M19_5);
        
        return top / rCoeffs(NU) + (rCoeffs(XI) * M2) + (rCoeffs(OMICRON) * (M8 * M0_5)) + (M6 * M6 * M6_5) + (rCoeffs(PI) * M19_5);
        
    #undef coeff
    }



/*
 * Calculate luminosity at the end of the Main Sequence
 *
 * Hurley et al. 2000, eq 8
 *
 *
 * double CalculateLuminosityAtPhaseEnd(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Luminosity at the end of the Main Sequence in Lsol
 */
double MainSequence::CalculateLuminosityAtPhaseEnd(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    // pow() is slow - use multiplication
    double m_3 = p_Mass * p_Mass * p_Mass;
    double m_4 = m_3 * p_Mass;
    double m_5 = m_4 * p_Mass;

    double top    = (an[11] * m_3) + (an[12] * m_4) + (an[13] * PPOW(p_Mass, (an[16] + 1.8)));
    double bottom = an[14] + (an[15] * m_5) + PPOW(p_Mass, an[16]);

    return top / bottom;

#undef a
}


/*
 * Calculate luminosity on the Main Sequence
 *
 * Hurley et al. 2000, eq 12 OR BRCEK <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Fix this
 *
 *
 * double CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS)
 *
 * @param   [IN]    p_Time                      Time (after ZAMS) in Myr
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_LZAMS                     Zero Age Main Sequence (ZAMS) Luminosity
 * @return                                      Luminosity on the Main Sequence as a function of time
 */
double MainSequence::CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const {
    
    double luminosity = 0.0;                                                                                                            // default return value

    // if BRCEK core prescription is used, return luminosity from Shikauchi et al. (2024) during core
    // hydrogen burning (valid for MZAMS >= 15 Msol) or luminosity that smoothly connects MS and HG
    // during MS hook (valid for MZAMS >= BRCEK_LOWER_MASS_LIMIT).
    // do not use Shikauchi luminosity prescription during CHE
    if (BRCEK_MS_CORE_MASS_REGIME(m_MZAMS) && 
       !m_CHE) {
        if (utils::Compare(p_Time, 0.99 * timescales(tMS)) > 0)                                                                         // star in MS hook?
        luminosity = CalculateLuminosityTransitionToHG(p_Mass, p_Time, p_LZAMS);
        else {
            if (utils::Compare(m_MZAMS, 15.0) >= 0)                                                                                     // use Shikauchi luminosity if MZAMS >= 15 Msun
            luminosity = CalculateLuminosityShikauchi(m_MainSequenceCoreMass, m_HeliumAbundanceCore);
        }
    }
    else {
    const double epsilon = 0.01;

    double LTMS   = CalculateLuminosityAtPhaseEnd(p_Mass);
    double alphaL = CalculateAlphaL(p_Mass);
    double betaL  = CalculateBetaL(p_Mass);
    double deltaL = CalculateDeltaL(p_Mass);
    double eta    = CalculateEta(p_Mass);

    double mu     = std::max(0.5, (1.0 - (0.01 * std::max((an[6] / PPOW(p_Mass, an[7])), (an[8] + (an[9] / PPOW(p_Mass, an[10])))))));  // Hurley et al. 2000, eq 7
    double tHook  = mu * timescales(tBGB);                                                                                              // Hurley et al. 2000, just after eq 5
    double tau    = p_Time / timescales(tMS);                                                                                           // Hurley et al. 2000, eq 11
    double tau1   = std::min(1.0, (p_Time / tHook));                                                                                    // Hurley et al. 2000, eq 14
    double tau2   = std::max(0.0, std::min(1.0, (p_Time - ((1.0 - epsilon) * tHook)) / (epsilon * tHook)));                             // Hurley et al. 2000, eq 15

    // pow() is slow - use multipliaction where it makes sense
    double logLMS_LZAMS  = alphaL * tau;                                                                                                // Hurley et al. 2000, eq 12, part 1
           logLMS_LZAMS += betaL * PPOW(tau, eta);                                                                                      // Hurley et al. 2000, eq 12, part 2
           logLMS_LZAMS += (log10(LTMS / p_LZAMS) - alphaL - betaL) * tau * tau;                                                        // Hurley et al. 2000, eq 12, part 3
           logLMS_LZAMS -= deltaL * ((tau1 * tau1) - (tau2 * tau2));                                                                    // Hurley et al. 2000, eq 12, part 4

    luminosity = p_LZAMS * PPOW(10.0, logLMS_LZAMS);                                                                                          // rewrite Hurley et al. 2000, eq 12 for L(t)
    }

    return luminosity;
}


/*
 * Calculate luminosity on the Main Sequence when BRCEK core mass prescription is used
 *
 * During core hydrogen burning uses eq (A5) from Shikauchi et al. (2024), valid for stars with MZAMS >= 15 Msol
 *
 * double CalculateLuminosityShikauchi(const double p_CoreMass, const double p_HeliumAbundanceCore)
 *
 * @param   [IN]    p_CoreMass                  Main sequence core mass in Msol
 * @param   [IN]    p_HeliumAbundanceCore       Central helium fraction
 * @return                                      Luminosity on the Main Sequence as a function of current core mass and central helium fraction
 */
double MainSequence::CalculateLuminosityShikauchi(const double p_CoreMass, const double p_HeliumAbundanceCore) const {
    DBL_VECTOR L_COEFFICIENTS = std::get<2>(SHIKAUCHI_COEFFICIENTS);
    
    // common factors
    double logMixingCoreMass   = std::log10(p_CoreMass);
    double logMixingCoreMass_2 = logMixingCoreMass * logMixingCoreMass;
    double logMixingCoreMass_3 = logMixingCoreMass_2 * logMixingCoreMass;
    
    double heliumAbundanceCore_2 = p_HeliumAbundanceCore * p_HeliumAbundanceCore;
    double heliumAbundanceCore_3 = heliumAbundanceCore_2 * p_HeliumAbundanceCore;
    
    double logL = L_COEFFICIENTS[0] * logMixingCoreMass + L_COEFFICIENTS[1] * p_HeliumAbundanceCore + L_COEFFICIENTS[2] * logMixingCoreMass * p_HeliumAbundanceCore + L_COEFFICIENTS[3] * logMixingCoreMass_2 + L_COEFFICIENTS[4] * heliumAbundanceCore_2 + L_COEFFICIENTS[5] * logMixingCoreMass_3 + L_COEFFICIENTS[6] * heliumAbundanceCore_3 + L_COEFFICIENTS[7] * logMixingCoreMass_2 * p_HeliumAbundanceCore + L_COEFFICIENTS[8] * logMixingCoreMass * heliumAbundanceCore_2 + L_COEFFICIENTS[9] * logMixingCoreMass_3 * logMixingCoreMass + L_COEFFICIENTS[10] * heliumAbundanceCore_3 * p_HeliumAbundanceCore + L_COEFFICIENTS[11] * logMixingCoreMass * heliumAbundanceCore_3 + L_COEFFICIENTS[12] * logMixingCoreMass_2 * heliumAbundanceCore_2 + L_COEFFICIENTS[13] * logMixingCoreMass_3 * p_HeliumAbundanceCore + L_COEFFICIENTS[14];
    
    return PPOW(10.0, logL);
}


/*
 * Calculate luminosity on the Main Sequence
 *
 * Hurley et al. 2000, eq 12 OR BRCEK <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Fix this
 *
 *
 * double CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS)
 *
 * @param   [IN]    p_Time                      Time (after ZAMS) in Myr
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_LZAMS                     Zero Age Main Sequence (ZAMS) Luminosity
 * @return                                      Luminosity on the Main Sequence as a function of time
 */
double MainSequence::CalculateLuminosityOnPhase(const double p_Time, const double p_Mass, const double p_LZAMS) const {

    // if BRCEK core prescription is used, return luminosity from Shikauchi et al. (2024) during core
    // hydrogen burning (valid for MZAMS >= 15 Msol) or luminosity that smoothly connects MS and HG
    // during MS hook (valid for MZAMS >= BRCEK_LOWER_MASS_LIMIT).
    // do not use Shikauchi luminosity prescription during CHE
    if (m_StateHistory.CurrentState().stellarType != STELLER_TYPE::CH                   &&  // CH star?
        BRCEK_MS_CORE_MASS_REGIME(m_StateHistory.StateZero().mass)) {                    &&  // yes - in Brcek regime?
                                                                                            // yes
        if (p_Time > (0.99 * timescales(tMS))) {                                            // in MS hook?                                                                                       
            luminosity = CalculateLuminosityTransitionToHG(p_Mass, p_Time, p_LZAMS);    // yes -  DESCRIBE <<<<<<<<<<<<<<<<<<<<<<<<<<<< 
        }
        else {                                                                              // no - no in MS hook
            if (utils::Compare(m_MZAMS, 15.0) >= 0)                                                                                     // use Shikauchi luminosity if MZAMS >= 15 Msun
            luminosity = CalculateLuminosityShikauchi(m_MainSequenceCoreMass, m_HeliumAbundanceCore);
        }
    }


    // if the Brcek MS core mass prescription was specified by the user, and the star is in the 
    // Brcek regime, and is on the MS hook, calculate the MS luminosity per Brcek et al. 2025
    if (BRCEK_MS_CORE_MASS_REGIME(p_MZAMS)                                               && // yes - in Brcek regime?
        tau > 0.99) {                                                                       // yes - on MS hook?
                                                                                            // yes
        retVal = CalculateRadiusOnPhase_Brcek(mass, tau, m_StateHistory.StateZero().radius, radius, GLOBALS->bCoefficients);
    }
    else {
        Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

            EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
            EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
                retVal = CalculateRadiusOnPhase_Hurley();
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


}


/*
 * CalculateLuminosityOnPhase_Brcek
 *
 * @brief
 * Calculate luminosity on the Main Sequence per Brcek et al. 2025
 *
 * Calculate luminosity to smoothly connect the beginning of MS hook and the beginning of HG (@TAMS)
 * The MS hook starts at tau = 0.99 on the main sequence
 *
 * GNU_CONST double CalculateLuminosityOnPhase_Brcek(const double p_Mass, const double p_Tau, const double p_RZAMS, const double p_Radius, const DBL_VECTOR p_bCoefficients) const
 *
 * @param       p_Mass                          Current mass of the star (Msol)
 * @param       p_Tau                           Main sequence fractional age of the star
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_Radius                        Current radius of the star (Rsol)
 * @param       p_bCoefficients                 Hurley b(n) coefficients
 * @return                                      Main sequence radius of the star (Rsol)
 */
double MainSequence::CalculateLuminosityOnPhase_Brcek(const double p_Mass, const double p_Tau, const double p_RZAMS, const double p_Radius, const DBL_VECTOR p_bCoefficients) const {

    const double luminosityTAMS        = std::min(HG::CalculateLuminosityOnPhase_Static(p_Mass, p_Tau, p_RZAMS, p_bCoefficients), p_Radius);
    const double luminosityAtHookStart = CalculateRadiusOnPhase_Hurley(mass, 0.99, p_RZAMS); 
        
    double radius = (radiusAtHookStart * (1.0 - p_Tau) + radiusTAMS * (p_Tau - 0.99)) / 0.01;   // linear interpolation

    // if the star has been stripped below its initial core mass we need to adjust the radius

/*
 * Calculate luminosity on the transition from the Main Sequence to the HG when BRCEK core mass prescription is used
 *
 * Luminosity prescription from Shikauchi et al. (2024) cannot be used beyond the MS hook (beyond age 0.99 * tMS), and this
 * function smoothly connects the luminosity between the beginning of the hook and the beginning of the HG
 *
 *
 * double CalculateLuminosityTransitionToHG(const double p_Mass, const double p_Age, double const p_LZAMS)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Age                       Age in Myr
 * @param   [IN]    p_LZAMS                     Zero Age Main Sequence (ZAMS) luminosity
 * @return                                      Luminosity on the Main Sequence (for age between tHook and tMS)
 */
double MainSequence::CalculateLuminosityTransitionToHG(const double p_Mass, const double p_Age, double const p_LZAMS) const {

    HG *clone = HG::Clone(static_cast<HG&>(const_cast<MainSequence&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    double luminosityTAMS = clone->Luminosity();                                                                                // Get luminosity from clone (with updated Mass0)
    delete clone; clone = nullptr;                                                                                              // Return the memory allocated for the clone
    
    double ageAtHookStart = 0.99 * timescales(tMS);
    
    double luminosityAtHookStart;
    if (utils::Compare(m_MZAMS, std::max(15.0, BRCEK_LOWER_MASS_LIMIT)) >= 0)
        luminosityAtHookStart = CalculateLuminosityShikauchi(m_MainSequenceCoreMass, m_HeliumAbundanceCore);                    // In the hook, core helium abundance fixed at 1-Z and core mass is not changing
    else
        luminosityAtHookStart = CalculateLuminosityOnPhase(ageAtHookStart, p_Mass, p_LZAMS);                                    // Do not use Shikauchi luminosity for MZAMS < 15 Msun
    
    // Linear interpolation
    return (luminosityAtHookStart * (timescales(tMS) - p_Age) + luminosityTAMS * (p_Age - ageAtHookStart)) / (timescales(tMS) - ageAtHookStart);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////




/*
 * Calculate the radius constant alpha_R
 * Hurley et al. 2000, eqs 21a & 21b
 *
 *
 * double CalculateAlphaR(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius constant alpha_R
 */
double MainSequence::CalculateAlphaR(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    double alphaR = 0.0;

         if (utils::Compare(p_Mass,   0.5) <  0) alphaR = an[62];
    else if (utils::Compare(p_Mass,  0.65) <  0) alphaR = an[62] + (an[63] - an[62]) * (p_Mass - 0.5) / 0.15;
    else if (utils::Compare(p_Mass, an[68]) <  0) alphaR = an[63] + (an[64] - an[63]) * (p_Mass - 0.65) / (an[68] - 0.65);
    else if (utils::Compare(p_Mass, an[66]) <  0) alphaR = an[64] + (m_RConstants[static_cast<int>(R_CONSTANTS::B_ALPHA_R)] - an[64]) * (p_Mass - an[68]) / (an[66] - an[68]);
    else if (utils::Compare(p_Mass, an[67]) <= 0) alphaR = an[58] * PPOW(p_Mass, an[60]) / (an[59] + PPOW(p_Mass, an[61]));
    else                                         alphaR = m_RConstants[static_cast<int>(R_CONSTANTS::C_ALPHA_R)] + an[65] * (p_Mass - an[67]);

    return alphaR;

#undef a
}


/*
 * Calculate the radius constant beta_R
 *
 * Hurley et al. 2000, eqs 22a & 22b
 *
 *
 * double Star::CalculateBetaR(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius constant beta_R
 */
double MainSequence::CalculateBetaR(const double p_Mass) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    double betaRPrime = 0.0;

         if (utils::Compare(p_Mass, 1.0)   <= 0) betaRPrime = 1.06;
    else if (utils::Compare(p_Mass, an[74]) <  0) betaRPrime = 1.06 + (an[72] - 1.06) * (p_Mass - 1.0) / (an[74] - 1.06);
    else if (utils::Compare(p_Mass, 2.0)   <  0) betaRPrime = an[72] + (m_RConstants[static_cast<int>(R_CONSTANTS::B_BETA_R)] - an[72]) * (p_Mass - an[74]) / (2.0 - an[74]);
    else if (utils::Compare(p_Mass, 16.0)  <= 0) betaRPrime = (an[69] * p_Mass * p_Mass * p_Mass * std::sqrt(p_Mass)) / (an[70] + PPOW(p_Mass, an[71]));  // pow()is slow - use multiplication (sqrt() is faster than pow())
    else                                         betaRPrime = m_RConstants[static_cast<int>(R_CONSTANTS::C_BETA_R)] + an[73] * (p_Mass - 16.0);

    return betaRPrime - 1.0;

#undef a
}


/*
 * Calculate the value of the radius perturbation DeltaR
 *
 * Hurley et al. 2000, eq 17
 *
 *
 * double CalculateDeltaR(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      The radius perturbation DeltaR
 */
double MainSequence::CalculateDeltaR(const double p_Mass) const {
#define an m_AnCoefficients                                              // for convenience and readability - undefined at end of function

    double deltaR;

    if (utils::Compare(p_Mass, massCutoffs(MHook)) <= 0) deltaR = 0.0;   // this really is supposed to be 0
    else if (utils::Compare(p_Mass, an[42])        <= 0) deltaR = an[43] * std::sqrt((p_Mass - massCutoffs(MHook)) / (an[42] - massCutoffs(MHook)));
    else if (utils::Compare(p_Mass, 2.0)           <  0) deltaR = an[43] + ((m_RConstants[static_cast<int>(R_CONSTANTS::B_DELTA_R)] - an[43]) * PPOW(((p_Mass - an[42]) / (2.0 - an[42])), an[44]));
    else {
        // pow() is slow - use multiplication (sqrt() is faster than pow())
        double top    = an[38] + (an[39] * p_Mass * p_Mass * p_Mass * std::sqrt(p_Mass));
        double bottom = (an[40] * p_Mass * p_Mass * p_Mass) + PPOW(p_Mass, an[41]);
        deltaR = (top / bottom) - 1.0;
    }

    return deltaR;

#undef an
}


/*
 * Calculate radius at the end of the Main Sequence
 *
 * Hurley et al. 2000, eqs 9a & 9b
 *
 *
 * double CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS)
 *
 * @param   [IN]    p_Mass                      Stellar mass (Msol)
 * @param   [IN]    p_RZAMS                     Zero Age Main Sequence (ZAMS) Radius
 * @return                                      Radius at the end of the Main Sequence in Rsol
 */
double MainSequence::CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_RZAMS) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    double RTMS;
    double mAsterisk = an[17] + 0.1;

    if (utils::Compare(p_Mass, an[17]) <= 0) {
        RTMS = (an[18] + (an[19] * PPOW(p_Mass, an[21]))) / (an[20] + PPOW(p_Mass, an[22]));

        if (utils::Compare(p_Mass, 0.5) < 0) {
            RTMS = std::max(RTMS, 1.5 * p_RZAMS);
        }
    }
    else if (utils::Compare(p_Mass, mAsterisk) >= 0) {
        // pow() is slow - use multiplication
        double m_3 = p_Mass * p_Mass * p_Mass;
        double m_5 = m_3 * p_Mass * p_Mass;

        RTMS = ((C_COEFF.at(1) * m_3) + (an[23] * PPOW(p_Mass, an[26])) + (an[24] * PPOW(p_Mass, an[26] + 1.5))) / (an[25] + m_5);
    }
    else {
        // for stars with masses between a17, a17 + 0.1 interpolate between the end points (y = mx + c)

        // pow() is slow - use multiplication
        double mA_3 = mAsterisk * mAsterisk * mAsterisk;
        double mA_5 = mA_3 * mAsterisk * mAsterisk;

        double y2   = ((C_COEFF.at(1) * mA_3) + (an[23] * PPOW(mAsterisk, an[26])) + (an[24] * PPOW(mAsterisk, an[26] + 1.5))) / (an[25] + mA_5);    // RTMS(mAsterisk)
        double y1   = (an[18] + (an[19] * PPOW(an[17], an[21]))) / (an[20] + PPOW(an[17], an[22]));                                                    // RTMS(a17)

        double gradient  = (y2 - y1) / 0.1;
        double intercept = y1 - (gradient * an[17]);

        RTMS = (gradient * p_Mass) + intercept;
    }

    return RTMS;

#undef a
}


/*
 * CalculateRadiusOnPhase_Hurley
 *
 * @brief
 * Calculate radius on the Main Sequence per Hurley et al. 2000
 *
 *
 * double CalculateRadiusOnPhase_Hurley(const double p_Mass, const double p_Tau, const double p_RZAMS, const DBL_VECTOR p_aCoefficients)
 *
 * @param       p_Mass                          Current mass of the star (Msol)
 * @param       p_Tau                           Main Sequence fractional age of the star
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_aCoefficients                 Hurley a(n) coefficients
 * @return                                      Main Sequence radius of the star (Rsol)
 */
double MainSequence::CalculateRadiusOnPhase_Hurley(const double p_Mass, const double p_Tau, const double p_RZAMS, const DBL_VECTOR p_aCoefficients) const {
#define a p_aCoefficients // for convenience and readability - undefined at end of function
    
    const double tBGB    = CalculateLifetimeToBGB(p_Mass);
    const double tMS     = CalculateLifetimeOnPhase(p_Mass, tBGB);

    const double RTMS    = CalculateRadiusAtPhaseEnd(p_Mass, p_RZAMS);
    const double alphaR  = CalculateAlphaR(p_Mass);
    const double betaR   = CalculateBetaR(p_Mass);
    const double deltaR  = CalculateDeltaR(p_Mass);
    const double gamma   = CalculateGamma(p_Mass);

    const double mu      = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10])))))));    // Hurley et al. 2000, eq 7
    const double tHook   = mu * tBGB;                                                                                                       // ibid, just after eq 5
    const double time    = tMS * p_Tau;
    const double tau_1   = std::min(1.0, (time / tHook));                                                                                   // ibid, eq 14
    const double tau_2   = std::max(0.0, std::min(1.0, (time - ((1.0 - 0.01) * tHook)) / (0.01 * tHook)));                                  // ibid, eq 15, epsilon= 0.01

    // pow() is slow - use multiplication where it makes sense
    const double tau3   = p_Tau * p_Tau * p_Tau;
    const double tau_13 = tau_1 * tau_1 * tau_1;
    const double tau_23 = tau_2 * tau_2 * tau_2;
    const double tau10  = p_Tau < FLOAT_TOLERANCE_ABSOLUTE ? 0.0: tau3 * tau3 * tau3 * p_Tau;                                               // tolerance comparison to avoid underflow
    const double tau40  = tau_10 < FLOAT_TOLERANCE_ABSOLUTE ? 0.0: tau10 * tau10 * tau10 * tau10;                                           // tolerance comparison to avoid underflow

    const double logRMS_RZAMS = alphaR * p_Tau + betaR * tau10 + gamma * tau40 + (log10(RTMS / p_RZAMS) - alphaR - betaR - gamma) * tau3 - deltaR * (tau_13 - tau_23); // ibid, eq 13

    return p_RZAMS * PPOW(10.0, logRMS_RZAMS);                                                                                              // rewrite Hurley et al. 2000, eq 13 for R(t)

#undef a
}

double CalculateRadiusOnPhase_Hurley() const {
    return CalculateRadiusOnPhase(m_StateHistory.CurrentState().mass, m_StateHistory.CurrentState().tau, m_StateHistory.StateZero().mass);
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
 * GNU_CONST double CalculateRadiusOnPhase_Brcek(const double p_Mass, const double p_Tau, const double p_RZAMS, const double p_Radius, const DBL_VECTOR p_bCoefficients) const
 *
 * @param       p_Mass                          Current mass of the star (Msol)
 * @param       p_Tau                           Main sequence fractional age of the star
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_Radius                        Current radius of the star (Rsol)
 * @param       p_bCoefficients                 Hurley b(n) coefficients
 * @return                                      Main sequence radius of the star (Rsol)
 */
double MainSequence::CalculateRadiusOnPhase_Brcek(const double p_Mass, const double p_Tau, const double p_RZAMS, const double p_Radius, const DBL_VECTOR p_bCoefficients) const {

    const double radiusTAMS        = std::min(HG::CalculateRadiusOnPhase_Static(p_Mass, p_Tau, p_RZAMS, p_bCoefficients), p_Radius);
    const double radiusAtHookStart = CalculateRadiusOnPhase_Hurley(mass, 0.99, p_RZAMS); 
        
    double radius = (radiusAtHookStart * (1.0 - p_Tau) + radiusTAMS * (p_Tau - 0.99)) / 0.01;   // linear interpolation

    // if the star has been stripped below its initial core mass we need to adjust the radius
  



    double heliumAbundanceSurface = m_HeliumAbundanceSurface;
        if (p_Mass < m_InitialMainSequenceCoreMass)
            // By tracing the helium profile in the star, calculate how the surface helium abundance changes if mass drops below the initial core mass
            heliumAbundanceSurface = m_HeliumAbundanceCoreOut + (p_Mass - m_MainSequenceCoreMass) * (m_HeliumAbundanceSurface - m_HeliumAbundanceCoreOut) / (m_InitialMainSequenceCoreMass - m_MainSequenceCoreMass);
        
        // Factor that scales radius based on surface helium abundance
        double surfaceAbundanceFactor = (utils::Compare(m_HeliumAbundanceCore, m_InitialHeliumAbundance) != 0) ? (heliumAbundanceSurface - m_InitialHeliumAbundance) / (m_HeliumAbundanceCore - m_InitialHeliumAbundance) : 0.0;
        
        radius += (p_RZAMS - radius) * surfaceAbundanceFactor;






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

    const double MZAMS  = m_StateHistory.StateZero().mass;
    const double tau    = m_StateHistory.CurrentState().tau;
    const double radius = m_StateHistory.CurrentState().radius;

    double retVal = radius;                                                                 // calculated radius - default is unchanged

    // if the user specified the BRCEK MS core mass prescription, and the star is in the 
    // BRCEK MS core mass prescription regine, and is on the MS hook, calculate the MS radius
    // per Brcek et al., 2025
    if (BRCEK_MS_CORE_MASS_REGIME(MZAMS)  && // yes - in Brcek MS core mass prescription regime?
        tau > 0.99) {                                                                       // yes - on MS hook?
                                                                                            // yes
        retVal = CalculateRadiusOnPhase_Brcek(mass, tau, m_StateHistory.StateZero().radius, radius, GLOBALS->bCoefficients);
    }
    else {
        Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

            EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
            EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
                retVal = CalculateRadiusOnPhase_Hurley();
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

    return retVal;
}


/*
 * Calculate radius on the transition from the Main Sequence to the HG when BRCEK core mass prescription is used
 *
 * Core mass prescription from Shikauchi et al. (2024) cannot be used beyond the MS hook (beyond Tau = 0.99), and this
 * function smoothly connects the radius between the beginning of the hook and the beginning of the HG
 *
 *
 * double CalculateRadiusTransitionToHG(const double p_Mass, const double p_Tau)
 
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       Fractional age on Main Sequence
 * @param   [IN]    p_RZAMS                     Zero Age Main Sequence (ZAMS) Radius
 * @return                                      Radius on the Main Sequence (for Tau between 0.99 and 1)
 */
double MainSequence::CalculateRadiusTransitionToHG(const double p_Mass, const double p_Tau, const double p_RZAMS) const {
    HG *clone = HG::Clone(static_cast<HG&>(const_cast<MainSequence&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    // Select radius at TAMS from the HG clone or current radius (whichever is smaller), relevant for stars that were significantly
    // stripped as this prevents radius expansion during the hook, and delays possible mass transfer to the start of HG
    double radiusTAMS = std::min(clone->Radius(), m_Radius);                                                                    // Get radius from clone (with updated Mass0)
    delete clone; clone = nullptr;                                                                                              // Return the memory allocated for the clone
    
    double radiusAtHookStart = CalculateRadiusOnPhase(p_Mass, 0.99, p_RZAMS);                                                   // Hook starts at Tau = 0.99
    
    return (radiusAtHookStart * (1.0 - p_Tau) + radiusTAMS * (p_Tau - 0.99)) / 0.01;                                            // Linear interpolation
}


/*
 * Calculate the radial extent of the star's convective envelope (if it has one)
 *
 * Hurley et al. 2002, sec. 2.3, particularly subsec. 2.3.1, eqs 36-38
 *
 *
 * double CalculateRadialExtentConvectiveEnvelope()
 *
 * @return                                      Radial extent of the star's convective envelope in Rsol
 */
double MainSequence::CalculateRadialExtentConvectiveEnvelope() const {
    double radiusEnvelope0 = m_Radius;

    if ( utils::Compare(m_Mass, 1.25) >= 0)
        radiusEnvelope0 = 0.0;
    else if (utils::Compare(m_Mass, 0.35) > 0) {
        // uses radius of a 0.35 solar mass star at ZAMS rather than at fractional age Tau,
        // but such low-mass stars only grow by a maximum factor of 1.5
        // [just above Eq. (10) in Hurley, Pols, Tout (2000)], so this is a reasonable approximation
        radiusEnvelope0 = CalculateRadiusAtZAMS(0.35) * std::sqrt((1.25 - m_Mass) / 0.9);
    }

    return radiusEnvelope0 * std::sqrt(std::sqrt(1.0 - m_Tau));
}


/*
 * Calculate the radial extent of the star's convective core (if it has one)
 *
 * Uses preliminary fit from Minori Shikauchi @ ZAMS, then a smooth interpolation to the HG
 *
 *
 * double CalculateRadialExtentConvectiveEnvelope()
 *
 * @return                                      Radial extent of the star's convective core in Rsol
 */
double MainSequence::CalculateConvectiveCoreRadius() const {
    if(utils::Compare(m_Mass, 1.25) < 0) return 0.0;                                            // low-mass star with a radiative core
       
    double convectiveCoreRadiusZAMS = m_Mass * (0.06 + 0.05 * exp(-m_Mass / 61.57));
    
    // We need TAMSCoreRadius, which is just the core radius at the start of the HG phase.
    // Since we are on the main sequence here, we can clone this object as an HG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.

    HG *clone = HG::Clone(static_cast<HG&>(const_cast<MainSequence&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    double TAMSCoreRadius = clone->CalculateRemnantRadius();                                    // get core radius from clone
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone

    return (convectiveCoreRadiusZAMS - m_Tau * (convectiveCoreRadiusZAMS - TAMSCoreRadius));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the mass of the convective core
 *
 * Based on Shikauchi, Hirai, Mandel (2024), core mass shrinks to 60% of initial value over the course of the MS
 *
 *
 * double CalculateConvectiveCoreMass() const
 *
 * @return                                      Mass of convective core in Msol
 */
double MainSequence::CalculateConvectiveCoreMass() const {
    double finalConvectiveCoreMass   = CalculateTAMSCoreMass();                 // core mass at TAMS
    double initialConvectiveCoreMass = finalConvectiveCoreMass / 0.6;
    return (initialConvectiveCoreMass - m_Tau * (initialConvectiveCoreMass - finalConvectiveCoreMass));
}


/*
 * Calculate the mass of the convective envelope
 *
 * Based on section 7.2 (after Eq. 111) of Hurley, Pols, Tout (2000)
 *
 *
 * DBL_DBL CalculateConvectiveEnvelopeMass() const
 *
 * @return                                      Tuple containing current mass of convective envelope and mass at ZAMS in Msol
 */
DBL_DBL MainSequence::CalculateConvectiveEnvelopeMass() const {
    if (utils::Compare(m_Mass, 1.25) > 0) return std::tuple<double, double> (0.0, 0.0);

    double massEnvelope0 = m_Mass;
    if (utils::Compare(m_Mass, 0.35) > 0) massEnvelope0 = 0.35 * (1.25 - m_Mass) * (1.25 - m_Mass) / 0.81;
    
    double massEnvelope  = massEnvelope0 * sqrt(sqrt(1.0 - m_Tau));
    
    return std::tuple<double, double> (massEnvelope, massEnvelope0);
}


/*
 * CalculateCoreMass_Brcek
 *
 * @brief
 * Calculate the convective core mass and central helium fraction of a main sequence star that
 * loses mass (e.g. through winds, case A mass transfer, or both), per Shikauchi et al. 2024
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
 * GNU_PURE double CalculateCoreMass_Brcek(const double p_dt, const double p_dMdt) const
 *
 * @param       p_Mass              Current mass of star (Msol)
 * @param       p_dt                Time step (Myr)
 * @param       p_dMdt              Mass loss rate in (Msol yr^-1)   //// >>>>>> CHECK THIS (negative for mass loss, positive for mass gain) <<<<<<<<<<<<<<<<
 * @param       p_CoreMass          Current MS core mass of star (Msol)
 * @return                          MS convective core mass
 */
double MainSequence::CalculateCoreMass_Brcek(const double p_Mass, const double p_dt, const double p_dMdt, const double p_CoreMass) const { 

    // get Shikauski coefficients from GLOBALS.  These (should) have  been calculated based
    // on the metallicity of the star.  If the coefficients haven't been calculated, throw
    // an appropriate error (this is a coding issue).

    std::optional<DBL_VECTOR> alphaCoeffs = GLOBALS->ShikauskiAlphaCoefficients();                                          // get Shikauchi alpha coefficients from GLOBALS
    if (!alphaCoeffs.has_value()) {                                                                                         // Shikauski alpha coefficients exist?
        THROW_ERROR(ERROR::NO_SHIKAUSKI_ALPHA_COEFFS);                                                                      // no - throw error
    }

    std::optional<DBL_VECTOR> fMixCoeffs = GLOBALS->ShikauskifMixCoefficients();                                            // get Shikauchi fMix coefficients from GLOBALS
    if (!fMixCoeffs.has_value()) {                                                                                          // Shikauski fMix coefficients exist?
        THROW_ERROR(ERROR::NO_SHIKAUSKI_FMIX_COEFFS);                                                                       // no - throw error
    }

    deltaCoeffs = SHIKAUCHI_DELTA_COEFFICIENTS;                                                                             // from constants.h

    auto fmix = [&](double mass) { return fMixCoeffs[0] + fMixCoeffs[1] * std::exp(-mass / fMixCoeffs[2]); };               // Shikauchi et al. 2024, eq A3


    double alpha = PPOW(10.0, std::max(-2.0, alphaCoeffs[1] * m_MainSequenceCoreMass + alphaCoeffs[2])) + alphaCoeffs[0];   // ibid, eq (A2)
    double g     = deltaCoeffs[1] * m_MainSequenceCoreMass + deltaCoeffs[2];                                                // ibid, eq (A7)
    
    double delta;
    if (p_dMdt <= 0.0)
        delta = std::min(PPOW(10.0, -deltaCoeffs[0] * (m_HeliumAbundanceCore - m_InitialHeliumAbundance) / (1.0 - m_InitialHeliumAbundance - m_Metallicity) + g), 1.0);          // ibid, eq (A6)
    else
        delta = PPOW(2.0, -(m_HeliumAbundanceCore - m_InitialHeliumAbundance) / (1.0 - m_InitialHeliumAbundance - m_Metallicity));                              // updated prescription for mass gain
    
    double deltaYc              = CalculateLuminosityOnPhase() / (Q_CNO * m_MainSequenceCoreMass) * p_Dt;                                                       // Change in central helium fraction; ibid, eq (12)
    double deltaMass            = p_dMdt * p_Dt * MYR_TO_YEAR;                                                                                          // Total mass lost/gained
    double deltaCoreMassML      = m_MainSequenceCoreMass * delta * ((p_Mass + deltaMass) * fmix(p_Mass + deltaMass) / (p_Mass * fmix(p_Mass)) - 1);             // Change in core mass due to mass loss/gain
    double deltaCoreMassNatural = -alpha / (1 - alpha * m_HeliumAbundanceCore) * deltaYc * m_MainSequenceCoreMass;                                              // Change in core mass due to natural decay; ibid, eq (4)
    double deltaCoreMass        = deltaCoreMassNatural + deltaCoreMassML;                                                                                       // Total difference in core mass
    


    double newMixingCoreMass        = std::min(m_MainSequenceCoreMass + deltaCoreMass, BRCEK_CORE_MASS_TO_MASS_RATIO_LIMIT * (p_Mass + deltaMass));             // New mixing core mass, always has to be smaller than the total mass







    double newCentralHeliumFraction = std::min(m_HeliumAbundanceCore + deltaYc, 1.0 - m_Metallicity);                                                           // New central helium fraction, capped at 1-Z

    if (deltaCoreMass > 0.0) {                                                                                                                                  // If the core grows, we need to account for rejuvenation
        if (utils::Compare(newMixingCoreMass, m_InitialMainSequenceCoreMass) < 0) {                                                                             // New core mass less than initial core mass?
            // Common factors
            double f1 = m_HeliumAbundanceCoreOut - m_InitialHeliumAbundance;
            double f2 = m_MainSequenceCoreMass - m_InitialMainSequenceCoreMass;
            double f3 = m_MainSequenceCoreMass + deltaCoreMass;

            // Calculate change in helium abundance just outside the core
            double deltaYout = f1 / f2 * deltaCoreMass;

            // Calculate the change in core helium abundance, assuming linear profile between Yc and Y0, and that the the accreted gas has helium fraction Y0
            double deltaY             = (m_HeliumAbundanceCoreOut - m_HeliumAbundanceCore) / f3 * deltaCoreMass + 0.5 / f3 * f1 / f2 * deltaCoreMass * deltaCoreMass;
            newCentralHeliumFraction  = m_HeliumAbundanceCore + deltaY;
            m_HeliumAbundanceCoreOut += deltaYout;
        }
        else {                                                                                                                                                  // New core mass greater or equal to the initial core mass?
            double deltaCoreMass1         = m_InitialMainSequenceCoreMass - m_MainSequenceCoreMass;                                                             // Mass accreted up to the initial core mass
            double deltaCoreMass2         = deltaCoreMass - deltaCoreMass1;                                                                                     // Remaining accreted mass
            newCentralHeliumFraction      = (m_MainSequenceCoreMass * m_HeliumAbundanceCore + 0.5 * (m_HeliumAbundanceCoreOut + m_InitialHeliumAbundance) * deltaCoreMass1 + deltaCoreMass2 * m_InitialHeliumAbundance) / (m_MainSequenceCoreMass + deltaCoreMass);
            m_HeliumAbundanceCoreOut      = m_InitialHeliumAbundance;
            m_InitialMainSequenceCoreMass = newMixingCoreMass;
        }
    }
    else {                                                                                                                                                      // Core decayed
        // If total mass dropped below the initial core mass, partially processed material is exposed and surface abundance needs to be adjusted
        if (utils::Compare(p_Mass + deltaMass, m_InitialMainSequenceCoreMass) < 0) {
            // Set surface helium abundance following the helium abundance profile in the star
            m_HeliumAbundanceSurface      = m_HeliumAbundanceCoreOut + (p_Mass + deltaMass - m_MainSequenceCoreMass) * (m_HeliumAbundanceSurface - m_HeliumAbundanceCoreOut) / (m_InitialMainSequenceCoreMass - m_MainSequenceCoreMass);
            m_HydrogenAbundanceSurface    = 1.0 - m_Metallicity - m_HeliumAbundanceSurface;
            m_InitialMainSequenceCoreMass = p_Mass + deltaMass;                                                                                                 // Update the initial core mass
        }
        m_HeliumAbundanceCoreOut = newCentralHeliumFraction;                                                                                                    // If core did not grow, Y_out = Y_c
    }
    
    return std::tuple<double, double> (newMixingCoreMass, std::min(newCentralHeliumFraction, 1.0 - m_Metallicity));
}


/*
 * Calculate the initial convective core mass of a main sequence star after full mixing (due to merger or CHE)
 * for an arbitrary central helium fraction using the approach described in Brcek et al. (2025)
 *
 * double CalculateInitialMainSequenceCoreMass(const double p_Mass, const double p_HeliumAbundanceCore)
 *
 * @param   [IN]    p_Mass                      Mass after merger or after spin down of CH star in Msol
 * @param   [IN]    p_HeliumAbundanceCore       Central helium fraction
 * @return                                      Mass of the convective core at ZAMS or after merger in Msol
 */
double MainSequence::CalculateInitialMainSequenceCoreMass(const double p_Mass, const double p_HeliumAbundanceCore) const {
    
    // After full mixing not at ZAMS, use the approach from Brcek+ (2025)
    double h = PPOW(10.0, p_HeliumAbundanceCore * (p_HeliumAbundanceCore + 2.0) / 4.0);
    double fmix = (BRCEK_FMIX_COEFFICIENTS[0] + BRCEK_FMIX_COEFFICIENTS[1] * std::exp(-p_Mass * h / BRCEK_FMIX_COEFFICIENTS[2])) * PPOW(1.0 - BRCEK_FMIX_COEFFICIENTS[4] / (p_Mass * h), BRCEK_FMIX_COEFFICIENTS[3]);

    return fmix * p_Mass;
}


/*
 * CalculateCoreMassOnPhase
 *
 * @brief
 * Calculate the core mass on the main sequence star
 * /////// that loses mass through winds or Case A mass transfer
 * When BRCEK core prescription is used, also update the core helium abundance and effective age  ///////
 *
 *
 * void UpdateMainSequenceCoreMass(const double p_Dt, const double p_MassLossRate)
 *
 * @param   [IN]      p_Dt                      Current timestep in Myr
 * @param   [IN]      p_MassLossRate            Mass loss rate either from stellar winds or mass transfer in Msol yr-1
 */
void MainSequence::UpdateMainSequenceCoreMass(const double p_Dt, const double p_MassLossRate) {

    double mainSequenceCoreMass = m_MainSequenceCoreMass;                                                                               // default is no change
    double heliumAbundanceCore  = m_HeliumAbundanceCore;                                                                                // default is no change
    double age                  = m_Age;                                                                                                // default is no change

    switch (OPTIONS->MainSequenceCoreMassPrescription()) {
        case MS_CORE_MASS_PRESCRIPTION::ZERO: 
            mainSequenceCoreMass = 0.0;
            break;
        
        case MS_CORE_MASS_PRESCRIPTION::MANDEL: 
            // Calculate the minimum core mass of a main sequence star that loses mass through 
            // Case A mass transfer as the core mass of a TAMS star, scaled by the fractional age.
            // Only applied to donors as part of binary evolution, not applied to SSE
            if ((OPTIONS->RetainCoreMassDuringCaseAMassTransfer()) && (p_MassLossRate < 0.0) && (utils::Compare(p_MassLossRate, -m_Mdot) != 0))
                mainSequenceCoreMass = std::max(m_MainSequenceCoreMass, CalculateTauOnPhase() * CalculateTAMSCoreMass());
            break;
        
        case MS_CORE_MASS_PRESCRIPTION::BRCEK:
            // Set core mass following Shikauchi et al. (2024) and account for rejuvenation if core grows
            if (utils::Compare(m_MZAMS, BRCEK_LOWER_MASS_LIMIT) >= 0) {                                                                 // BRCEK prescription valid?
                // Only proceed with calculations if star is not in MS hook (Yc < 1-Z) and time step is not zero
                if ((utils::Compare(m_HeliumAbundanceCore, 1.0 - m_Metallicity) < 0) && (utils::Compare(p_Dt, 0.0) != 0)) {
                    // Update the core mass and central helium fraction only if the mass loss rate argument is equal
                    // to the total mass loss rate (i.e. total mass loss rate was updated, this prevents the calculation
                    // in SSE if it was executed as part of BSE for the same time step)
                    if (utils::Compare(p_MassLossRate, m_TotalMassLossRate) == 0) {
                        // Calculate and update the core mass and central helium fraction
                        std::tie(mainSequenceCoreMass, heliumAbundanceCore) = CalculateCoreMass_Brcek(p_Dt, p_MassLossRate);
                        // Update effective age here only if core hydrogen was exhausted
                        age = heliumAbundanceCore == 1.0 - m_Metallicity ? 0.99 * timescales(tMS) : age;
                    }
                    // Update effective age only when stars are aged in SSE (when p_MassLossRate = -Mdot)
                    if (utils::Compare(p_MassLossRate, -m_Mdot) == 0)
                        // Update the effective age based on central helium fraction
                        age = (heliumAbundanceCore - m_InitialHeliumAbundance) / m_InitialHydrogenAbundance * 0.99 * timescales(tMS);
                }
            }
            else {                                                                                                                      // no - MANDEL prescription used
                // Only applied to donors as part of binary evolution, not applied to SSE
                if ((p_MassLossRate < 0.0) && (utils::Compare(p_MassLossRate, -m_Mdot) != 0))
                    mainSequenceCoreMass = std::max(m_MainSequenceCoreMass, CalculateTauOnPhase() * CalculateTAMSCoreMass());
            }
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


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate relative age on the Main Sequence
 *
 * Hurley et al. 2000, eq 11
 * Naturally bounded by [0, 1], but clamp here anyway
 *
 * double CalculateTauOnPhase()
 *
 * @return                                      MS relative age, clamped to [0, 1]
 */
double MainSequence::CalculateTauOnPhase() const {
    return std::max(0.0, std::min(1.0, m_Age / timescales(tMS)));
}


/*
 * Calculate lifetime of Main Sequence
 *
 * Hurley et al. 2000, eq 5
 *
 *
 * double CalculatePhaseLifetime(const double p_Mass, const double p_TBGB)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_TBGB                      Lifetime to Base of Giant Branch
 * @return                                      Lifetime of Main Sequence in Myr
 */
double MainSequence::CalculateLifetimeOnPhase(const double p_Mass, const double p_TBGB) const {
#define a m_AnCoefficients    // for convenience and readability - undefined at end of function

    // Calculate time to Hook
    // Hurley et al. 2000, eqs 5, 6 & 7
    double mu    = std::max(0.5, (1.0 - (0.01 * std::max((an[6] / PPOW(p_Mass, an[7])), (an[8] + (an[9] / PPOW(p_Mass, an[10])))))));
    double tHook = mu * p_TBGB;

    // For mass < Mhook, x > mu (i.e. for stars without a hook)
    double x = std::max(0.95, std::min((0.95 - (0.03 * (LogMetallicityXiHurley() + 0.30103))), 0.99));

    return std::max(tHook, (x * p_TBGB));

#undef a
}


/*
 * Recalculates the star's age after mass loss
 *
 * Hurley et al. 2000, section 7.1
 *
 * Modifies attribute m_Age
 *
 *
 * UpdateAgeAfterMassLoss()
 *
 */
void MainSequence::UpdateAgeAfterMassLoss() {

    double tMS       = m_Timescales[static_cast<int>(TIMESCALE::tMS)];
    double tBGBprime = CalculateLifetimeToBGB(m_Mass);
    double tMSprime  = MainSequence::CalculateLifetimeOnPhase(m_Mass, tBGBprime);
    
    m_Age *= tMSprime / tMS;
    CalculateTimescales(m_Mass, m_Timescales);                                      // must update timescales
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
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Choose timestep for evolution
 *
 * Given in the discussion in Hurley et al. 2000
 *
 *
 * ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double MainSequence::ChooseTimestep(const double p_Time) const {

    double dtk = 1.0E-2 * timescales(tMS);  // 0.01 of MS timescale (sse uses 0.05)
    double dte = timescales(tMS) - p_Time;  // time remaining on MS

    if (utils::Compare(dte, dtk) < 0) {     // short enough to resolve the hook at the end of the MS for HM stars? JAR: why not check for HM star?
        dtk /= 10.0;                        // no - go an order-of-magnitude shorter
    }

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);
}


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
 * Calculate the expected core mass at terminal age main sequence, i.e., at the start of the HG phase
 *
 * double calculate TAMSCoreMass() const
 *
 *
 * @return                                      TAMS core Mass (Msol)
 *
 */
double MainSequence::CalculateTAMSCoreMass() const {
    // Since we are on the main sequence here, we can clone this object as an HG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.
    
    HG *clone = HG::Clone(static_cast<HG&>(const_cast<MainSequence&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    double TAMSCoreMass = clone->CoreMass();                                                    // get core mass from clone
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone
    
    return TAMSCoreMass;
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
    
    CalculateTimescales();
    CalculateGBParams();
            
    m_Tau = (initialHydrogenFraction - p_HydrogenMass / m_Mass) / initialHydrogenFraction;      // assumes uniformly mixed merger product and a uniform rate of H fusion on main sequence
    
    m_Age = m_Tau * timescales(tMS);
    
    m_HeliumAbundanceCore   = 1.0 - m_Metallicity - p_HydrogenMass / p_Mass;
    m_HydrogenAbundanceCore = 1.0 - m_Metallicity - m_HeliumAbundanceCore;
    
    m_HeliumAbundanceSurface   = m_HeliumAbundanceCore;                                         // abundances are the same throughout the star, assuming uniform mixing after merger
    m_HydrogenAbundanceSurface = m_HydrogenAbundanceCore;
    
    if (BRCEK_MS_CORE_MASS_REGIME(m_MZAMS)) {
        m_InitialMainSequenceCoreMass = CalculateInitialMainSequenceCoreMass(p_Mass, m_HeliumAbundanceCore);           // update initial mixing core mass
        m_MainSequenceCoreMass        = m_InitialMainSequenceCoreMass;                                                 // update core mass
    }
    
    EvolveOneTimestep(0.0, 0.0, 0.0, true);
}


/* 
 * Interpolate Ge+ Critical Mass Ratios, for H-rich stars
 * 
 * Function takes input QCRIT_PRESCRIPTION, currently either of the prescriptions for critical mass ratios
 * from Ge et al. (2020), GE or GE_IC. The first is the full adiabatic response, the second assumes
 * artificially isentropic envelopes. From private communication with Ge, we have an updated datatable that
 * includes qCrit for fully conservative and fully non-conservative MT, so we now interpolate on those as well.
 *
 * Interpolation is done linearly in logM, logR, and logZ
 * 
 * double BaseStar::InterpolateGeEtAlQCrit(const QCRIT_PRESCRIPTION p_qCritPrescription, const double p_massTransferEfficiencyBeta) 
 * 
 * @param   [IN]    p_qCritPrescription          Adopted critical mass ratio prescription
 * @param   [IN]    p_massTransferEfficiencyBeta Mass transfer accretion efficiency
 * @return                                       Interpolated value of either the critical mass ratio or zeta for given stellar mass / radius
 */ 
double MainSequence::InterpolateGeEtAlQCrit(const QCRIT_PRESCRIPTION p_qCritPrescription, const double p_massTransferEfficiencyBeta) {

    // Iterate over the two QCRIT_GE tables to get the qcrits at each metallicity
    double qCritPerMetallicity[2];
    std::vector<GE_QCRIT_TABLE> qCritTables = { QCRIT_GE_LOW_Z, QCRIT_GE_HIGH_Z };

    for (int ii=0; ii<2; ii++) { // iterate over the vector of tables to store qCrits per metallicity

        // Get vector of masses from qCritTable
        GE_QCRIT_TABLE &qCritTable = qCritTables[ii];
        DBL_VECTOR massesFromQCritTable = std::get<0>(qCritTable);
        GE_QCRIT_RADII_QCRIT_VECTOR radiiQCritsFromQCritTable = std::get<1>(qCritTable);

        INT_VECTOR indices = utils::BinarySearch(massesFromQCritTable, m_Mass);
        int lowerMassIndex = indices[0];
        int upperMassIndex = indices[1];
    
        if (lowerMassIndex == -1) {                                                   // if masses are out of range, set to endpoints
            lowerMassIndex = 0; 
            upperMassIndex = 1;
        } 
        else if (upperMassIndex == -1) { 
            lowerMassIndex = massesFromQCritTable.size() - 2; 
            upperMassIndex = massesFromQCritTable.size() - 1;
        } 
    
        // Get vector of radii from qCritTable for the lower and upper mass indices
        std::vector<double> logRadiusVectorLowerMass = std::get<0>(radiiQCritsFromQCritTable[lowerMassIndex]);
        std::vector<double> logRadiusVectorUpperMass = std::get<0>(radiiQCritsFromQCritTable[upperMassIndex]);
    
        // Get the qCrit vector for the lower and upper mass bounds 
        std::vector<double> qCritVectorUpperEffLowerMass;
        std::vector<double> qCritVectorUpperEffUpperMass;
        std::vector<double> qCritVectorLowerEffLowerMass;
        std::vector<double> qCritVectorLowerEffUpperMass;
        
        // Set the appropriate qCrit vector, depends on MT eff and whether you use GE STD or IC
        if (p_qCritPrescription == QCRIT_PRESCRIPTION::GE) {
            if (p_massTransferEfficiencyBeta > 0.5) {
                qCritVectorUpperEffLowerMass = std::get<1>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorUpperEffUpperMass = std::get<1>(radiiQCritsFromQCritTable[upperMassIndex]);
                qCritVectorLowerEffLowerMass = std::get<2>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorLowerEffUpperMass = std::get<2>(radiiQCritsFromQCritTable[upperMassIndex]);
            }
            else {
                qCritVectorUpperEffLowerMass = std::get<2>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorUpperEffUpperMass = std::get<2>(radiiQCritsFromQCritTable[upperMassIndex]);
                qCritVectorLowerEffLowerMass = std::get<3>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorLowerEffUpperMass = std::get<3>(radiiQCritsFromQCritTable[upperMassIndex]);

            }
        }
        else if (p_qCritPrescription == QCRIT_PRESCRIPTION::GE_IC) {
            if (p_massTransferEfficiencyBeta > 0.5) {
                qCritVectorUpperEffLowerMass = std::get<4>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorUpperEffUpperMass = std::get<4>(radiiQCritsFromQCritTable[upperMassIndex]);
                qCritVectorLowerEffLowerMass = std::get<5>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorLowerEffUpperMass = std::get<5>(radiiQCritsFromQCritTable[upperMassIndex]);
            }
            else {
                qCritVectorUpperEffLowerMass = std::get<5>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorUpperEffUpperMass = std::get<5>(radiiQCritsFromQCritTable[upperMassIndex]);
                qCritVectorLowerEffLowerMass = std::get<6>(radiiQCritsFromQCritTable[lowerMassIndex]);
                qCritVectorLowerEffUpperMass = std::get<6>(radiiQCritsFromQCritTable[upperMassIndex]);

            }
        }
    
        // Get vector of radii from qCritTable for both lower and upper masses
        INT_VECTOR indicesR0          = utils::BinarySearch(logRadiusVectorLowerMass, log10(m_Radius));
        int lowerRadiusLowerMassIndex = indicesR0[0];
        int upperRadiusLowerMassIndex = indicesR0[1];
    
        if (lowerRadiusLowerMassIndex == -1) {                                        // if radii are out of range, set to endpoints
            lowerRadiusLowerMassIndex = 0; 
            upperRadiusLowerMassIndex = 1; 
        }
        else if (upperRadiusLowerMassIndex == -1) {                                                   
            lowerRadiusLowerMassIndex = logRadiusVectorLowerMass.size() - 2; 
            upperRadiusLowerMassIndex = logRadiusVectorLowerMass.size() - 1; 
        }
    
        INT_VECTOR indicesR1          = utils::BinarySearch(logRadiusVectorUpperMass, log10(m_Radius));
        int lowerRadiusUpperMassIndex = indicesR1[0];
        int upperRadiusUpperMassIndex = indicesR1[1];
    
        if (lowerRadiusUpperMassIndex == -1) {                                        // if radii are out of range, set to endpoints
            lowerRadiusUpperMassIndex = 0; 
            upperRadiusUpperMassIndex = 1; 
        }
        else if (upperRadiusUpperMassIndex == -1) {                                                   
            lowerRadiusUpperMassIndex = logRadiusVectorUpperMass.size() - 2; 
            upperRadiusUpperMassIndex = logRadiusVectorUpperMass.size() - 1; 
        }
    
        // Set the 4 boundary points for the 2D interpolation
        double qUppLowLow = qCritVectorUpperEffLowerMass[lowerRadiusLowerMassIndex];
        double qUppLowUpp = qCritVectorUpperEffLowerMass[upperRadiusLowerMassIndex];
        double qUppUppLow = qCritVectorUpperEffUpperMass[lowerRadiusUpperMassIndex];
        double qUppUppUpp = qCritVectorUpperEffUpperMass[upperRadiusUpperMassIndex];
        double qLowLowLow = qCritVectorLowerEffLowerMass[lowerRadiusLowerMassIndex];
        double qLowLowUpp = qCritVectorLowerEffLowerMass[upperRadiusLowerMassIndex];
        double qLowUppLow = qCritVectorLowerEffUpperMass[lowerRadiusUpperMassIndex];
        double qLowUppUpp = qCritVectorLowerEffUpperMass[upperRadiusUpperMassIndex];
    
        double lowerLogRadiusLowerMass = logRadiusVectorLowerMass[lowerRadiusLowerMassIndex];
        double upperLogRadiusLowerMass = logRadiusVectorLowerMass[upperRadiusLowerMassIndex];
        double lowerLogRadiusUpperMass = logRadiusVectorUpperMass[lowerRadiusUpperMassIndex];
        double upperLogRadiusUpperMass = logRadiusVectorUpperMass[upperRadiusUpperMassIndex];

        double logLowerMass   = log10(massesFromQCritTable[lowerMassIndex]);
        double logUpperMass   = log10(massesFromQCritTable[upperMassIndex]);
    
        // Interpolate on logR first, then logM, then on the efficiency, using nearest neighbor for extrapolation
        double logRadius = log10(m_Radius);
        double qCritUpperEffLowerMass = (logRadius < lowerLogRadiusLowerMass) ? qUppLowLow
                                      : (logRadius > upperLogRadiusLowerMass) ? qUppLowUpp
                                      : qUppLowLow + (upperLogRadiusLowerMass - logRadius) / (upperLogRadiusLowerMass - lowerLogRadiusLowerMass) * (qUppLowUpp - qUppLowLow);
        double qCritUpperEffUpperMass = (logRadius < lowerLogRadiusUpperMass) ? qUppUppLow
                                      : (logRadius > upperLogRadiusUpperMass) ? qUppUppUpp
                                      : qUppUppLow + (upperLogRadiusUpperMass - logRadius) / (upperLogRadiusUpperMass - lowerLogRadiusUpperMass) * (qUppUppUpp - qUppUppLow);
        double qCritLowerEffLowerMass = (logRadius < lowerLogRadiusLowerMass) ? qLowLowLow
                                      : (logRadius > upperLogRadiusLowerMass) ? qLowLowUpp
                                      : qLowLowLow + (upperLogRadiusLowerMass - logRadius) / (upperLogRadiusLowerMass - lowerLogRadiusLowerMass) * (qLowLowUpp - qLowLowLow);
        double qCritLowerEffUpperMass = (logRadius < lowerLogRadiusUpperMass) ? qLowUppLow
                                      : (logRadius > upperLogRadiusUpperMass) ? qLowUppUpp
                                      : qLowUppLow + (upperLogRadiusUpperMass - logRadius) / (upperLogRadiusUpperMass - lowerLogRadiusUpperMass) * (qLowUppUpp - qLowUppLow);
    
        double logMass = log10(m_Mass);
        double interpolatedQCritUpperEff = (logMass < logLowerMass) ? qCritUpperEffLowerMass
                                         : (logMass > logUpperMass) ? qCritUpperEffUpperMass
                                         : qCritUpperEffLowerMass + (logUpperMass - logMass) / (logUpperMass - logLowerMass) * (qCritUpperEffUpperMass - qCritUpperEffLowerMass);
        double interpolatedQCritLowerEff = (logMass < logLowerMass) ? qCritLowerEffLowerMass
                                         : (logMass > logUpperMass) ? qCritLowerEffUpperMass
                                         : qCritLowerEffLowerMass + (logUpperMass - logMass) / (logUpperMass - logLowerMass) * (qCritLowerEffUpperMass - qCritLowerEffLowerMass);
    
        double interpolatedQCritForZ = p_massTransferEfficiencyBeta * interpolatedQCritUpperEff + (1.0 - p_massTransferEfficiencyBeta) * interpolatedQCritLowerEff;                 // Don't need to use nearest neighbor for this, beta is always between 0 and 1
        qCritPerMetallicity[ii] = interpolatedQCritForZ;
    }
    double logZlo = -3;                // log10(0.001)
    double logZhi = LOG10_ZSOL_HURLEY; // log10(0.02)
    
    return qCritPerMetallicity[1] + (m_Log10Metallicity - logZhi)*(qCritPerMetallicity[1] - qCritPerMetallicity[0])/(logZhi - logZlo);
}
