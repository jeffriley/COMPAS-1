#include "Globals.h"

using std::max;
using std::min;

// initialise instance pointer
Globals* Globals::m_Instance = nullptr;


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         HURLEY COEFFICIENTS AND CONSTANTS                         //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateACoefficients_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 'a' coefficients.
 * 
 * Hurley a coefficients depend on the star's metallicity only, so this only needs to be done
 * at most once per star (upon creation), but can also be reused if metallicity doesn't change
 * from one star to the next (e.g. in a population run).
 * 
 * Values are from table in Appendix A of Hurley et al. 2000.
 * 
 * 
 * DblVectorT CalculateACoefficients_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta) const
 * 
 * @param       p_Z                             Metallicity
 * @param       p_Sigma                         Sigma from Hurley et al. 2000 sigma, p24
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      a coefficients vector
 */
DblVectorT Globals::CalculateACoefficients_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta) const {

    // Create and inialise a coefficients vector - this is the return value.
    DblVectorT a(HURLEY_A_COEFF.size() + 1, DEFAULT_INITIAL_DOUBLE_VALUE);

    // Calculate some powers of zeta - for performance and readability.
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;

    // Populate a coefficients vector with initial values for a coefficients.
    // Iterate over a coefficients constants HURLEY_A_COEFF (see constants.h).
    // Each row (indexed by coeff.first (int)) defines the coefficients of the 5 terms
    // (HURLEY_AB_TCoeff coefficients 'ALPHA', 'BETA', 'GAMMA', 'ETA', 'MU').
    // The first entry (index 0) is a dummy entry - so our index is the same as that in
    // Hurley et al. 2000 (we just ignore the zeroeth entry).
    for (const auto& [idx, vals] : HURLEY_A_COEFF) {
        using enum AB_TCoeff;
        auto coeff = [&](AB_TCoeff p) { return vals[static_cast<SizeT>(p)]; };
        a[idx] = coeff(ALPHA) + coeff(BETA) * p_Zeta + coeff(GAMMA) * zeta2 + coeff(ETA) * zeta3 + coeff(MU) * zeta4;
    }    

    // Special cases - see Hurley et al. 2000.
    a[11] *= a[14];
    a[12] *= a[14];
    a[17]  = PPOW(10.0, max((0.097 - (0.1072 * (p_Sigma + 3.0))), max(0.097, min(0.1461, (0.1461 + (0.1237 * (p_Sigma + 2.0)))))));
    a[18] *= a[20];
    a[19] *= a[20];
    a[29]  = PPOW(a[29], (a[32]));
    a[33]  = max(0.6355 - (0.4192 * p_Zeta), max(1.25, min(1.4, 1.5135 + (0.3769 * p_Zeta))));
    a[42]  = min(1.25, max(1.1, a[42]));
    a[44]  = min(1.3, max(0.45, a[44]));
    a[49]  = max(a[49], 0.145);
    a[50]  = min(a[50], (0.306 + (0.053 * p_Zeta)));
    a[51]  = min(a[51], (0.3625 + (0.062 * p_Zeta)));
    a[52]  = max(a[52], 0.9);;
    a[52]  = p_Z > 0.01 ? min(a[52], 1.0) : a[52];
    a[53]  = max(a[53], 1.0);
    a[53]  = p_Z > 0.01 ? min(a[53], 1.1) : a[53];
    a[57]  = max((0.6355 - (0.4192 * p_Zeta)), max(1.25, min(1.4, a[57])));
    a[62]  = max(0.065, a[62]);
    a[63]  = p_Z < 0.004 ? min(0.055, a[63]) : a[63];
    a[66]  = max(0.8, min(0.8 - (2.0 * p_Zeta), max(a[66], min(1.6, -0.308 - (1.046 * p_Zeta)))));
    a[68]  = max(0.9, min(a[68], 1.0));
    a[72]  = p_Z > 0.01 ? max(a[72], 0.95) : a[72];
    a[74]  = max(1.4, min(a[74], 1.6));
    a[75]  = max(1.0, min(a[75], 1.27));
    a[75]  = max(a[75], 0.6355 - (0.4192 * p_Zeta));
    a[76]  = max(a[76], -0.1015564 - (0.2161264 * p_Zeta) - (0.05182516 * zeta2));
    a[77]  = max((-0.3868776 - (0.5457078 * p_Zeta) - (0.1463472 * zeta2)), min(0.0, a[77]));
    a[78]  = max(0.0, min(a[78], 7.454 + (9.046 * p_Zeta)));
    a[79]  = min(a[79], max(2.0, -13.3 - (18.6 * p_Zeta)));
    a[80]  = max(0.0585542, a[80]);
    a[81]  = min(1.5, max(0.4, a[81]));

    // (Re)compute special cases that depend on radius constants.
    // Uses alpahR = (a[58] * PPOW(a[66], a[60])) / (a[59] + PPOW(a[66], a[61])) as given in
    // Hurley et al. 2000, eq 21a (wrong in the arxiv version - says = a59*M**(a61))
    a[64] = a[68] > a[66] ? (a[58] * PPOW(a[66], a[60])) / (a[59] + PPOW(a[66], a[61])) : max(0.091, min(0.121, a[64]));
    a[68] = min(a[68], a[66]);  // must calculate after calculating a[64]

    // Return the a coefficients vector by value - NRVO takes care of performance/efficiency.
    return a;
}
    
    
/*
 * CalculateBCoefficients_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 'b' coefficients.
 * 
 * Hurley b coefficients depend on the star's metallicity only, so this only needs to be done
 * at most once per star (upon creation), but can also be reused if metallicity doesn't change
 * from one star to the next (e.g. in a population run).
 * 
 * Values are from table in Appendix A of Hurley et al. 2000.
 *
 *
 * DblVectorT CalculateBCoefficients_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta, const double p_Rho) const
 * 
 * @param       p_Z                             Metallicity
 * @param       p_Sigma                         Sigma from Hurley et al. 2000 sigma, p24
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @param       p_Rho                           Rho from Hurley et al. 2000 sigma, p24
 * @return                                      b coefficients vector
 */
DblVectorT Globals::CalculateBCoefficients_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta, const double p_Rho) const {

    // Create and inialise b coefficients vector - this is the return value.
    DblVectorT b(HURLEY_B_COEFF.size() + 1, DEFAULT_INITIAL_DOUBLE_VALUE);

    // Calculate some powers of zeta and rho - for performance and readability.
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;
    const double zeta5 = p_Zeta * zeta4;
    const double rho2  = p_Rho * p_Rho;
    const double rho3  = p_Rho * rho2;

    // Populate b coefficients vector with initial values for b coefficients.
    // Iterate over b coefficients constants HURLEY_B_COEFF (see constants.h)
    // each row (indexed by coeff.first (int)) defines the coefficients of the 5 terms
    // (HURLEY_AB_TCoeff coefficients 'ALPHA', 'BETA', 'GAMMA', 'ETA', 'MU').
    // The first entry (index 0) is a dummy entry - so our index is the same as that in
    // Hurley et al. 2000 (we just ignore the zeroeth entry).
    for (const auto& [idx, vals] : HURLEY_B_COEFF) {
        using enum AB_TCoeff;
        auto coeff = [&](AB_TCoeff p) { return vals[static_cast<SizeT>(p)]; };
        b[idx] = coeff(ALPHA) + coeff(BETA) * p_Zeta + coeff(GAMMA) * zeta2 + coeff(ETA) * zeta3 + coeff(MU) * zeta4;
    }

    // Special cases - see Hurley et al. 2000.
    b[ 1]  = min(0.54, b[1]);
    b[ 2]  = min(max(PPOW(10.0, (-4.6739 - (0.9394 * p_Sigma))), (-0.04167 + (55.67 * p_Z))), (0.4771 - (9329.21 * PPOW(p_Z, 2.94))));
    b[ 3]  = PPOW(10.0, max(-0.1451, (-2.2794 - (1.5175 * p_Sigma) - (0.254 * p_Sigma * p_Sigma))));
    b[ 3]  = p_Z > 0.004 ? max(b[3], 0.7307 + (14265.1 * PPOW(p_Z, 3.395))) : b[3];
    b[ 4] += 0.1231572 * zeta5;
    b[ 6] += 0.01640687 * zeta5;
    b[11] *= b[11];
    b[13] *= b[13];
    b[14]  = PPOW(b[14], b[15]);
    b[16]  = PPOW(b[16], b[15]);
    b[17]  = xi > -1.0 ? 1.0 - (0.3880523 * PPOW(p_Rho, 2.862149)) : 1.0;
    b[24]  = PPOW(b[24], b[28]);
    b[26]  = 5.0 - (0.09138012 * PPOW(p_Z, -0.3671407));
    b[27]  = PPOW(b[27], (2.0 * b[28]));
    b[31]  = PPOW(b[31], b[33]);
    b[34]  = PPOW(b[34], b[33]);
    b[36] *= b[36] * b[36] * b[36];
    b[37] *= 4.0;
    b[38] *= b[38] * b[38] * b[38];
    b[40]  = max(b[40], 1.0);
    b[41]  = PPOW(b[41], b[42]);
    b[44] *= b[44] * b[44] * b[44] * b[44];
    b[45]  = p_Rho <= 0.0 ? 1.0 : 1.0 - ((2.47162 * p_Rho) - (5.401682 * rho2) + (3.247361 * rho3));
    b[46]  = -1.0 * b[46] * std::log10(MassCutoffs(HURLEY_MCO::HEF) / MassCutoffs(HURLEY_MCO::FGB));
    b[47]  = (1.127733 * p_Rho) + (0.2344416 * rho2) - (0.3793726 * rho3);
    b[51] -= 0.1343798 * zeta5;
    b[53] += 0.4426929 * zeta5;
    b[55]  = min((0.99164 - (743.123 * PPOW(p_Z, 2.83))), b[55]);
    b[56] += 0.1140142 * zeta5;
    b[57] -= 0.01308728 * zeta5;

    // Return the b coefficients vector by value - NRVO takes care of performance/efficiency.
    return b;
}


/*
 * CalculateAlphas_Hurley2000
 *
 * @brief
 * Calculate Hurley constants alpha1, alpha3, and alpha4
 *
 *   - alpha1 per Hurley et al, 2000, p12, just after eq 49
 *   - alpha3 per Hurley et al, 2000, p12, just after eq 56
 *   - alpha4 per Hurley et al, 2000, p12, just after eq 57
 *
 * alpha1, alpha3, and alpha4 depend on the star's metallicity only - so this only needs to be done once
 * per star (upon creation)
 * 
 * alpha2, per Hurley et al, 2000, p12, just after eq 53, depends on the star's core mass, so is not constant
 * and is not calculated here.
 *
 *
 * DblVectorT CalculateAlphas_Hurley2000(const double p_MHeF, const DblVectorT& p_bCoeffs) const
 * 
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_bCoeffs                       Hurley b coefficients
 * @return                                      Alpha values vector
 */
DblVectorT Globals::CalculateAlphas_Hurley2000(const double p_MHeF, const DblVectorT& p_bCoeffs) const {

    // Create and inialise alphas vector - this is the return value.
    DblVectorT alphas(4, DEFAULT_INITIAL_DOUBLE_VALUE);

    const auto& b = p_bCoeffs;
    
    // alpha1
    const double lHeI_MHeF = (b[11] + (b[12] * PPOW(p_MHeF, 3.8))) / (b[13] + (p_MHeF * p_MHeF));

    alphas[0] = ((b[9] * PPOW(p_MHeF, b[10])) - lHeI_MHeF) / lHeI_MHeF;

    // alpha3
    const double lBAGB = (b[31] + (b[32] * PPOW(p_MHeF, (b[33] + 1.8)))) / (b[34] + PPOW(p_MHeF, b[33]));

    alphas[1] = ((b[29] * PPOW(p_MHeF, b[30])) - lBAGB) / lBAGB;

    // alpha4
    const double tBGB_MHeF = utils::CalculateLifetimeToBGB_Hurley2000(p_MHeF); // tBGB for mass M = MHeF
    const double mHeF5     = utils::intPow(p_MHeF, 5);
    const double tHe_MHeF  = tBGB_MHeF * (b[41] * PPOW(p_MHeF, b[42]) + b[43] * mHeF5) / (b[44] + mHeF5); 

    alphas[2] = ((tHe_MHeF - b[39]) / b[39]);

    // Return the alphas vector by value - NRVO takes care of performance/efficiency
    return alphas;
}    


/*
 * CalculateGammaConstants_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 gamma constants
 * 
 * Gamma constants depend on the star's metallicity only - so this only needs to be done at
 * most once per star (upon creation), but can also be reused if metallicity doesn't change 
 * from one star to the next (e.g. in a population run).
 * 
 * The gamma constant values are calculated using Hurley et al. 2000:
 * 
 *    - B_GAMMA: eq 23, with M = 1.0 (see discussion immediately following eq 23)
 *    - C_GAMMA: see discussion immediately following eq 23
 * 
 * 
 * DblVectorT CalculateGammaConstants_Hurley2000(const DblVectorT& p_aCoeffs) const
 * 
 * @param       p_aCoeffs                       Hurley a coefficients
 * @return                                      Gamma constants vector
 */
DblVectorT Globals::CalculateGammaConstants_Hurley2000(const DblVectorT& p_aCoeffs) const {
               
    using enum HURLEY_GAMMA_CONSTANT;

    // Create and inialise gamma constants vector - this is the return value.
    DblVectorT gammas(utils::to_underlying(COUNT), DEFAULT_INITIAL_DOUBLE_VALUE);
        
    // Populate gamma constants vector.
    // Hurley et al. 2000, eq 23 and discussion following.
    const auto& a = p_aCoeffs;
    double bGamma = max(0.0, a[76] + (a[77] * PPOW((1.0 - a[78]), a[79])));
    gammas[utils::to_underlying(B_GAMMA)] = bGamma;
    gammas[utils::to_underlying(C_GAMMA)] = a[75] <= 1.0 ? bGamma : a[80];
        
    // Return the gamma constants vector by value - NRVO takes care of performance/efficiency.
    return gammas;
}


/*
 * CalculateLuminosityConstants_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 luminosity constants
 * 
 * Luminosity constants depend on the star's metallicity only - so this only needs to be done
 * at most once per star (upon creation), but can also be reused if metallicity doesn't change 
 * from one star to the next (e.g. in a population run).
 * 
 * The luminosity constant values are calculated using Hurley et al. 2000:
 * 
 *    - B_DELTA_L: eq 16, with M = a(33) (see discussion immediately following eq 16)
 *    - B_ALPHA_L: eq 19, with M = 2.0 (see discussion immediately following eq 19a)
 *    - B_BETA_L : eq 20
 * 
 * 
 * DblVectorT CalculateLuminosityConstants_Hurley2000(const DblVectorT& p_aCoefficients) const
 * 
 * @param       p_aCoeffs                       Hurley a coefficients
 * @return                                      Luminosity constants vector
 */
DblVectorT Globals::CalculateLuminosityConstants_Hurley2000(const DblVectorT& p_aCoeffs) const {
   
    using enum HURLEY_L_CONSTANT;

    // Create and inialise luminosity constants vector - this is the return value.
    DblVectorT lums(utils::to_underlying(COUNT), DEFAULT_INITIAL_DOUBLE_VALUE);

    // Populate luminosity constants vector.
    const auto& a = p_aCoeffs;
    lums[utils::to_underlying(B_DELTA_L)] = min((a[34] / PPOW(a[33], a[35])), (a[36] / PPOW(a[33], a[37])));                    // Hurley et al. 2000, eq 16
    lums[utils::to_underlying(B_ALPHA_L)] = (a[45] + (a[46] * PPOW(2.0, a[48]))) / (PPOW(2.0, 0.4) + (a[47] * PPOW(2.0, 1.9))); // ibid., eq 19
    lums[utils::to_underlying(B_BETA_L)]  = max(0.0, (a[54] - (a[55] * PPOW(a[57], a[56]))));                                   // ibid., eq 20

    // Return the luminosity constants vector by value - NRVO takes care of performance/efficiency.
    return lums;
}


/*
 * CalculateMassCutoffs_Hurley2000
 *
 * @brief
 * Calculate mass cutoffs per Hurley et al. 2000:
 *
 *   Hook: the mass above which a hook appears on the MS (Hurley et al. 2000, eq 1)
 *   HeF : the maximum initial mass for which He ignites degenerately in the He Flash (Hurley et al. 2000, eq 2)
 *   FGB : the maximum mass at which He ignites degenerately on the First Giant Branch (FGB) (Hurley et al. 2000, eq 3)
 *
 * Hurley mass cutoffs depend on the star's metallicity only - so this only needs to be done at most once per star
 * (upon creation), but can also be reused if metallicity doesn't change from one star to the next (e.g. in a population run).
 *
 *
 * DblVectorT CalculateMassCutoffs_Hurley2000(const double p_Z, const double p_Zeta) const
 * 
 * @param       p_Z                             Metallicity
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      Mass cutoffs vector
 */
DblVectorT Globals::CalculateMassCutoffs_Hurley2000(const double p_Z, const double p_Zeta) const {
   
    using enum HURLEY_MASS_CUTOFFS;

    // Create and inialise mass cutoffs vector - this is the return value.
    DblVectorT massCutoffs(utils::to_underlying(COUNT), DEFAULT_INITIAL_DOUBLE_VALUE);

    // Calculate some powers of zeta - for performance and readability.
    const double zeta2 = p_Zeta * p_Zeta;

    // Populate mass cutoffs vector.
    massCutoffs[utils::to_underlying(Hook)] = 1.0185 + (0.16015 * p_Zeta) + (0.0892 * zeta2);
    massCutoffs[utils::to_underlying(HeF)]  = 1.995 + (0.25 * p_Zeta) + (0.087 * zeta2);  
    massCutoffs[utils::to_underlying(FGB)]  = 13.048 * PPOW((p_Z / ZSOL_HURLEY), 0.06) / (1.0 + (0.0012 * PPOW((ZSOL_HURLEY / p_Z), 1.27)));    

    // Return the mass cutoffs vector by value - NRVO takes care of performance/efficiency.
    return massCutoffs;
}


/*
 * CalculateRadiusConstants_Hurley2000
 *
 * @brief
 * Calculate Hurley et al. 2000 radius constants
 * 
 * Radius constants depend on the star's metallicity only - so this only needs to be done at
 * most once per star (upon creation), but can also be reused if metallicity doesn't change 
 * from one star to the next (e.g. in a population run).
 * 
 * The radius constant values are calculated using Hurley et al. 2000:
 * 
 *    - B_DELTA_R: eq 17, with M = 2.0 (see discussion immediately following eq 17)
 *    - B_ALPHA_R: eq 21, with M = a(66) (see discussion immediately following eq 21a)
 *    - C_ALPHA_R: eq 21, with M = a(67) (see discussion immediately following eq 21a)
 *    - B_BETA_R : eq 22, with M = 2.0 (see discussion immediately following eq 22a)
 *    - C_BETA_R : eq 22, with M = 16.0 (see discussion immediately following eq 22a)
 * 
 * 
 * DblVectorT CalculateRadiusConstants_Hurley2000(const DblVectorT& p_aCoeffs) const
 * 
 * @param       p_aCoeffs                       Hurley a coefficients
 * @return                                      Radius constants vector
 */
DblVectorT Globals::CalculateRadiusConstants_Hurley2000(const DblVectorT& p_aCoeffs) const {
     
    using enum HURLEY_R_CONSTANT;

    // Create and inialise radius constants vector - this is the return value.
    DblVectorT rads(utils::to_underlying(COUNT), DEFAULT_INITIAL_DOUBLE_VALUE);
    
    // Populate radius constants vector.
    const auto& a = p_aCoeffs;
    rads[utils::to_underlying(B_DELTA_R)] = (a[38] + a[39] * 8.0 * M_SQRT2) / (a[40] * 8.0 + PPOW(2.0, a[41])) - 1.0;   // Hurley et al. 2000, eq 17
    rads[utils::to_underlying(B_ALPHA_R)] = (a[58] * PPOW(a[66], a[60])) / (a[59] + PPOW(a[66], a[61]));                // ibid., eq 21 (wrong in the arxiv version - says = a59*M**(a61))
    rads[utils::to_underlying(C_ALPHA_R)] = (a[58] * PPOW(a[67], a[60])) / (a[59] + PPOW(a[67], a[61]));                // ibid., eq 21 (wrong in the arxiv version - as above)
    rads[utils::to_underlying(B_BETA_R)]  = (a[69] * 8.0 * M_SQRT2) / (a[70] + PPOW(2.0, a[71]));                       // ibid., eq 22
    rads[utils::to_underlying(C_BETA_R)]  = (a[69] * 16384.0) / (a[70] + PPOW(16.0, a[71]));                            // ibid., eq 22
    
    // Return the radius constants vector by value - NRVO takes care of performance/efficiency.
    return rads;
}


/*
 * CalculateGBRadiusXexponent_Hurley2000
 *
 * @brief
 * Calculate the parameter x for the Giant Branch ('x' exponent to which Radius depends on Mass
 * (at constant Luminosity), per Hurley et al. 2000, eq 47)
 * 
 * Hybrid of b5 and b7 from Hurley et al. 2000
 *
 * 'x' depends on the star's metallicity only - so this only needs to be done once per star (upon creation),
 * but can also be reused if metallicity doesn't change from one star to the next (e.g. in a population run).
 *
 *
 * double CalculateGBRadiusXexponent_Hurley2000(const double p_Zeta) const
 * 
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      Giant Branch radius 'x' exponent
 */
double Globals::CalculateGBRadiusXexponent_Hurley2000(const double p_Zeta) const {

    // Calculate some powers of zeta - for performance and readability.
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;

    return 0.30406 + (0.0805 * p_Zeta) + (0.0897 * zeta2) + (0.0878 * zeta3) + (0.0222 * zeta4);
}


/*
 * CalculateZdependentValues_Hurley2000
 *
 * @brief
 * Calculate Hurley Z-dependent values.
 * 
 * 
 * struct HurleyZdependentT CalculateZdependentValues_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta, const HurleyZdependentT& p_Zvalues) const
 * 
 * @param       p_Z                             Metallicity value to be used calculate Hurley Z-dependent values
 * @param       p_Sigma                         Sigma from Hurley et al. 2000 p24, sigma = log10(Z)
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @param       p_Zvalues                       Struct containing Hurley Z-dependent values
 * @return                                      Struct containing (possibly updated) Hurley Z-dependent values
 */
HurleyZdependentT CalculateZdependentValues_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta, const HurleyZdependentT& p_Zvalues) const {

    HurleyZdependentT values   = p_Zvalues;         // Return value - default is unchanged

    // Calculate and set member values.
    values.rho                 = p_Zeta + 1.0;      // Hurley et al. 2000 p24, rho = zeta + 1.0

    values.aCoefficients       = CalculateACoefficients_Hurley2000(p_Z, p_Sigma, p_Zeta);
    values.massCutoffs         = CalculateMassCutoffs_Hurley2000(p_Z, p_Zeta);
    values.bCoefficients       = CalculateBCoefficients_Hurley2000(p_Z, p_Sigma, p_Zeta, values.rho, values.massCutoffs);

    values.gammaConstants      = CalculateGammaConstants_Hurley2000(values.aCoefficients);
    values.luminosityConstants = CalculateLuminosityConstants_Hurley2000(values.aCoefficients);
    values.radiusConstants     = CalculateRadiusConstants_Hurley2000(values.aCoefficients);

    values.radiusXexponent     = CalculateGBRadiusXexponent_Hurley2000();

    values.alphas              = CalculateAlphas_Hurley2000(values.massCutoffs[utils::to_underlying(HURLEY_MCO::HeF)], values.bCoefficients);

    // Return the values struct by value - NRVO takes care of performance/efficiency.
    return values;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         OTHER COEFFICIENTS AND CONSTANTS                          //
//                                                                                   //
/////////////////////////////////////////////////////////////////////////////////////// 


/*
 * CalculateCoeffs_Shikauchi2024
 *
 * @brief
 * Calculate metallicity-dependent coefficients per Shikauchi et al. 2024.
 * Shikauchi et al. 2024 gives values for the coefficients of relationships for:
 * 
 *    - the natural decline rate of the mixing core mass in the absence of mass loss (alpha),
 *    - the fraction of the mass contained in the mixing core (fMix),
 *    - luminosity
 * 
 * The relationships described by Shikauchi et al. 2024 are metallicity-dependent.  This function
 * calculates the values for coefficents at a specified metallicity either by extending the extremities
 * of the relationships (assuming a constant value equal to the bound), or linear interpolation between
 * the bounds defined by Shikauchi et al. 2024.
 * 
 *
 * std::tuple<DblVectorT, DblVectorT, DblVectorT> CalculateCoeffs_Shikauchi2024(const double p_Z) const
 *
 * @param       p_Z                             Metallicity of the star
 * @return                                      Tuple containing vectors of coefficients (alpha, fMix, luminosity)
 */
std::tuple<DblVectorT, DblVectorT, DblVectorT> Globals::CalculateCoeffs_Shikauchi2024(const double p_Z) const {

    using fMixLower  = SHIKAUCHI_FMIX_COEFF[SHIKAUCHI_Coeff::ONE_TENTH_Z_SOL];
    using fMixMid    = SHIKAUCHI_FMIX_COEFF[SHIKAUCHI_Coeff::ONE_THIRD_Z_SOL];
    using fMixUpper  = SHIKAUCHI_FMIX_COEFF[SHIKAUCHI_Coeff::Z_SOL];

    using alphaLower = SHIKAUCHI_ALPHA_COEFF[SHIKAUCHI_Coeff::ONE_TENTH_Z_SOL];
    using alphaMid   = SHIKAUCHI_ALPHA_COEFF[SHIKAUCHI_Coeff::ONE_THIRD_Z_SOL];
    using alphaUpper = SHIKAUCHI_ALPHA_COEFF[SHIKAUCHI_Coeff::Z_SOL];

    using lumLower   = SHIKAUCHI_LUMINOSITY_COEFF[SHIKAUCHI_Coeff::ONE_TENTH_Z_SOL];
    using lumMid     = SHIKAUCHI_LUMINOSITY_COEFF[SHIKAUCHI_Coeff::ONE_THIRD_Z_SOL];
    using lumUpper   = SHIKAUCHI_LUMINOSITY_COEFF[SHIKAUCHI_Coeff::Z_SOL];

    // Create and initialise return values - vectors of coefficients.
    DblVectorT alphaCoeff(fMixLower.size(), 0.0);                           // Alpha coefficients
    DblVectorT fMixCoeff(SalphaLower.size(), 0.0);                          // fMix coefficients
    DblVectorT luminosityCoeff(lumLower.size(), 0.0);                       // Luminosity coefficients
       
    // Common factors for each of the metallicities defined by Shikauchi et al. 2024.
    const double logZ              = std::log10(p_Z);

    const double logLowerZ         = std::log10(0.1 * ZSOL_HURLEY);         // Z lower bound: SHIKAUCHI_Coeff::ONE_TENTH_Z_SOL
    const double logMidZ           = std::log10(1.0 / 3.0 * ZSOL_HURLEY);   // Z mid-range:   SHIKAUCHI_Coeff::ONE_THIRD_Z_SOL
    const double logUpperZ         = LOG10_ZSOL_HURLEY;                     // Z upper bound: SHIKAUCHI_Coeff::Z_SOL

    const double logZ_logLowerZ    = p_logZ - logLowerZ;
    const double logZ_logMidZ      = p_logZ - logMidZ;
    const double logMidZ_logZ      = logMidZ  - p_logZ;
    const double logMidZ_logLowerZ = logMidZ  - logLowerZ;
    const double logUpperZ_logZ    = logUpperZ - p_logZ;
    const double logUpperZ_logMidZ = logUpperZ - logMidZ;
    
    // Calculate coefficients for specified metallicity.
    if (p_logZ <= logLowestZ) {                                             // p_Z at or below lower metallicity bound?
                                                                            // Yes, clamp coefficients to lower bound values
        alphaCoeff      = alphaLower;
        fMixCoeff       = fMixLower;
        luminosityCoeff = lumLower;
    }
    else if (p_logZ <= middle) {                                            // p_Z in lower to middle metallicity band?
                                                                            // Yes, interpolate
        for (SizeT i = 0; i < alphaCoeff.size(); i++)
            alphaCoeff[i] = (alphaLower[i] * logMidZ_logZ + alphaMid[i] * logZ_logLowerZ) / logMidZ_logLowerZ;
        for (SizeT i = 0; i < fMixCoeff.size(); i++)
            fMixCoeff[i]  = (fMixLower[i] * logMidZ_logZ + fMixMid[i] * logZ_logLowerZ) / logMidZ_logLowerZ;
        for (SizeT i = 0; i < luminosityCoeff.size; i++)
            luminosityCoeff[i] = (lumLower[i] * logMidZ_logZ + lumMid[i] * logZ_logLowerZ) / logMidZ_logLowerZ;
    }
    else if (p_logZ < high) {                                               // p_Z in middle to upper metallicity band?
                                                                            // Yes, interpolate
        for (SizeT i = 0; i < alphaCoeff.size(); i++)
            alphaCoeff[i] = (alphaMid[i] * logUpperZ_logZ + alphaUpper[i] * logZ_logMidZ) / logUpperZ_logMidZ;
        for (SizeT i = 0; i < fMixCoeff.size(); i++)
            fMixCoeff[i]  = (fMixMid[i] * logUpperZ_logZ + fMixUpper[i] * logZ_logMidZ) / logUpperZ_logMidZ;
        for (SizeT i = 0; i < luminosityCoeff.size; i++)
            luminosityCoeff[i] = (lumMid[i] * logUpperZ_logZ + lumUpper[i] * logZ_logMidZ) / logUpperZ_logMidZ;
    }
    else {                                                                  // p_Z at or above upper metallicity bound
                                                                            // Clamp coefficients to upper bound values
        alphaCoeff      = alphaUpper;
        fMixCoeff       = fMixUpper;
        luminosityCoeff = lumUpper;
    }
    
    return std::make_tuple(alphaCoeff, fMixCoeff, luminosityCoeff);
}
