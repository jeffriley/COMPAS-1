#include "Globals.h"

using std::max;
using std::min;

// initialise instance pointer
Globals* Globals::m_Instance = nullptr;

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
// Some of the functions here are small enough to be candidates for inlining, but    //
// since we call these functions at most once per system (star or binary), there's   //
// no real advantage in them being inlined, and we don't bother splitting candidates //
// out to the header file.                                                           //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////




// ADMIN-ish


void Initialise() {

    // options are already initialised here


    // regardles of evolution maode (SSE or BSE), we will need to record details
    // for at east one star
    StarDetails* m_Star1 = new StarDetails;                 // star 1

    // assumption: metallicity is same for both stars in a binary.  The code here could easily be changed to have metallicity different for each star, but COMPAS options currently don't allow that, so for now we use that as a simplification.

    m_RefZ                   = OPTIONS->Metallicity();                                                      // reference metallicity

    m_ZAMSheliumAbundance    = CalculateZAMSHeliumAbundance_Pols_1998(m_RefZ);                              // ZAMS helium abundance
    m_ZAMShydrogenAbundance  = CalculateZAMSHydrogenAbundance_Pols_1998(m_RefZ);                            // ZAMS hydrogen abundance

    m_HurleyZdependentValues = CalculateHurleyZdependentValues(m_RefZ, m_HurleyZdependentValues);           // Hurley Z-dependent values

    m_LuminosityCoefficients = CalculateLuminosityCoefficients_Tout_1996(m_HurleyZdependentValues.zeta);
    m_RadiusCoefficients     = CalculateRadiusCoefficients_Tout_1996(m_HurleyZdependentValues.zeta);





    // first, determine what mode the user requested
    EVOLUTION_MODE evolutionMode = OPTIONS->EvolutionMode();
    
    if (evolutionMode == EVOLUTION_MODE::BSE_HURLEY) {      // binary mode?
                                                            // yes
        StarDetails* m_Star2 = new StarDetails;             // star 2

        // star 2 z-dependent values are the same as star 1
        m_Star2->refZ = m_Star1->refZ;

        m_Star2->ZAMSheliumAbundance   = m_Star1->ZAMSheliumAbundance;
        m_Star2->ZAMShydrogenAbundance = m_Star1->ZAMShydrogenAbundance;
            
    }


}





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                 (RE)CALCULATE HURLEY METALLICITY-DEPENDENT VALUES                 //
//                                                                                   // 
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         HURLEY COEFFICIENTS AND CONSTANTS                         //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateHurleyACoefficients
 *
 * @brief
 * Calculate Hurley et al. 2000 a(n) coefficients
 * 
 * a(n) coefficients depend on the star's metallicity only - so this only needs to be done at
 * most once per star (upon creation), but can also be reused if metallicity doesn't change 
 * from one star to the next (e.g. in a population run).
 * 
 * The a(n) values are from table in Appendix A of Hurley et al. 2000.
 * 
 * 
 * DBL_VECTOR CalculateHurleyACoefficients(const double p_Z, const double p_Sigma, const double p_Zeta) const
 * 
 * @param       p_Z                             Metallicity
 * @param       p_Sigma                         Sigma from Hurley et al. 2000 sigma, p24
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      a(n) coefficients vector
 */
DBL_VECTOR Globals::CalculateHurleyACoefficients(const double p_Z, const double p_Sigma, const double p_Zeta) const {
// macros for convenience and readability - undefined at end of function
#define index    coeff.first
#define coeff(x) coeff.second[AB_TCoeff::x]

    // create and inialise a(n) coefficients vector - this is the return value
    DBL_VECTOR a(HURLEY_A_COEFF.size() + 1, DEFAULT_INITIAL_DOUBLE_VALUE);

    // calculate some powers of zeta - for performance and readability
    // this function is only called once per star, and at most twice per binary (but probably once), so not too onerous
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;

    // populate a(n) coefficients vector with initial values for a(n) coefficients
    // iterate over a(n) coefficients constants HURLEY_A_COEFF (see constants.h)
    // each row (indexed by coeff.first (int)) defines the coefficients of the 5 terms (HURLEY_AB_TCoeff coefficients 'ALPHA', 'BETA', 'GAMMA', 'ETA', 'MU') 
    // first entry (index 0) is a dummy entry - so our index is the same as that in Hurley et al. 2000 (we just ignore the zeroeth entry)
    for (auto coeff: HURLEY_A_COEFF) a[index] = coeff(ALPHA) + (coeff(BETA) * p_Zeta) + (coeff(GAMMA) * zeta2) + (coeff(ETA) * zeta3) + (coeff(MU) * zeta4);
    
    // special cases - see Hurley et al. 2000
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

    // (re)compute a(n) special cases that depend on radius constants
    // uses alpahR = (a[58] * PPOW(a[66], a[60])) / (a[59] + PPOW(a[66], a[61])) as given in
    // Hurley et al. 2000, eq 21a (wrong in the arxiv version - says = a59*M**(a61))
    a[64] = a[68] > a[66] ? (a[58] * PPOW(a[66], a[60])) / (a[59] + PPOW(a[66], a[61])) : max(0.091, min(0.121, a[64]));
    a[68] = min(a[68], a[66]);  // must calculate after calculating a[64]

    // return the a(n) coefficients vector by value - NRVO takes care of performance/efficiency
    return a;

#undef coeff
#undef index
}
    
    
/*
 * CalculateHurleyBCoefficients
 *
 * @brief
 * Calculate Hurley et al. 2000 b(n) coefficients
 * 
 * b(n) coefficients depend on metallicity-dependent mass cutoffs, and the star's metallicity only - so
 * this only needs to be done at most once per star (upon creation), but can also be reused if metallicity
 * doesn't change from one star to the next (e.g. in a population run).
 * 
 * Values are from table in Appendix A of Hurley et al. 2000
 *
 *
 * DBL_VECTOR CalculateHurleyBCoefficients(const double p_Z, const double p_Sigma, const double p_Zeta, const double p_Rho, const DBL_VECTOR& p_MassCutoffs) const
 * 
 * @param       p_Z                             Metallicity
 * @param       p_Sigma                         Sigma from Hurley et al. 2000 sigma, p24
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @param       p_Rho                           Rho from Hurley et al. 2000 sigma, p24
 * @param       p_MassCutoffs                   Hurley mass cutoffs vector
 * @return                                      b(n) coefficients vector
 */
DBL_VECTOR Globals::CalculateHurleyBCoefficients(const double p_Z, const double p_Sigma, const double p_Zeta, const double p_Rho, const DBL_VECTOR& p_MassCutoffs) const {
// macros for convenience and readability - undefined at end of function
#define index    coeff.first
#define coeff(x) coeff.second[AB_TCoeff::x]
    
    // create and inialise b(n) coefficients vector - this is the return value
    DBL_VECTOR b(HURLEY_B_COEFF.size() + 1, DEFAULT_INITIAL_DOUBLE_VALUE);

    // calculate some powers of zeta and rho - for performance and readability
    // this function is only called once per star, and at most twice per binary (but probably once), so not too onerous
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;
    const double zeta5 = p_Zeta * zeta4;
    const double rho2  = p_Rho * p_Rho;
    const double rho3  = p_Rho * rho2;

    // populate b(n) coefficients vector with initial values for b(n) coefficients
    // iterate over b(n) coefficients constants HURLEY_B_COEFF (see constants.h)
    // each row (indexed by coeff.first (int)) defines the coefficients of the 5 terms (HURLEY_AB_TCoeff coefficients 'ALPHA', 'BETA', 'GAMMA', 'ETA', 'MU') 
    // first entry (index 0) is a dummy entry - so our index is the same as that in Hurley et al. 2000 (we just ignore the zeroeth entry)
    for (auto coeff: HURLEY_B_COEFF) b[index] = coeff(ALPHA) + (coeff(BETA) * p_Zeta) + (coeff(GAMMA) * zeta2) + (coeff(ETA) * zeta3) + (coeff(MU) * zeta4);

    // special cases - see Hurley et al. 2000
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
    b[46]  = -1.0 * b[46] * log10(p_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)] / p_MassCutoffs[static_cast<int>(MASS_CUTOFF::MFGB)]);
    b[47]  = (1.127733 * p_Rho) + (0.2344416 * rho2) - (0.3793726 * rho3);
    b[51] -= 0.1343798 * zeta5;
    b[53] += 0.4426929 * zeta5;
    b[55]  = min((0.99164 - (743.123 * PPOW(p_Z, 2.83))), b[55]);
    b[56] += 0.1140142 * zeta5;
    b[57] -= 0.01308728 * zeta5;

    // return the b(n) coefficients vector by value - NRVO takes care of performance/efficiency
    return b;

#undef coeff
#undef index
}


/*
 * CalculateHurleyLuminosityConstants
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
 * DBL_VECTOR CalculateHurleyLuminosityConstants(const DBL_VECTOR& p_aCoefficients) const
 * 
 * @param       p_aCoefficients                 Hurley a(n) coefficients
 * @return                                      Luminosity constants vector
 */
DBL_VECTOR Globals::CalculateHurleyLuminosityConstants(const DBL_VECTOR& p_aCoefficients) const {
#define a p_aCoefficients // for convenience and readability - undefined at end of function
   
    // create and inialise luminosity constants vector - this is the return value
    DBL_VECTOR lums(HURLEY_L_CONSTANT::COUNT, DEFAULT_INITIAL_DOUBLE_VALUE);

    // populate luminosity constants vector
    lums[static_cast<int>(HURLEY_L_CONSTANT::B_DELTA_L)] = min((a[34] / PPOW(a[33], a[35])), (a[36] / PPOW(a[33], a[37])));                     // Hurley et al. 2000, eq 16
    lums[static_cast<int>(HURLEY_L_CONSTANT::B_ALPHA_L)] = (a[45] + (a[46] * PPOW(2.0, a[48]))) / (PPOW(2.0, 0.4) + (a[47] * PPOW(2.0, 1.9)));  // Hurley et al. 2000, eq 19
    lums[static_cast<int>(HURLEY_L_CONSTANT::B_BETA_L)]  = max(0.0, (a[54] - (a[55] * PPOW(a[57], a[56]))));                                    // Hurley et al. 2000, eq 20

    // return the luminosity constants vector by value - NRVO takes care of performance/efficiency
    return lums;

#undef a
}


/*
 * CalculateHurleyRadiusConstants
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
 * DBL_VECTOR CalculateHurleyRadiusConstants(const DBL_VECTOR& p_aCoefficients) const
 * 
 * @param       p_aCoefficients                 Hurley a(n) coefficients
 * @return                                      Radius constants vector
 */
DBL_VECTOR Globals::CalculateHurleyRadiusConstants(const DBL_VECTOR& p_aCoefficients) const {
#define a p_aCoefficients // for convenience and readability - undefined at end of function
           
    // create and inialise radius constants vector - this is the return value
    DBL_VECTOR rads(HURLEY_R_CONSTANT::COUNT, DEFAULT_INITIAL_DOUBLE_VALUE);
    
    // populate radius constants vector
    rads[static_cast<int>(HURLEY_R_CONSTANT::B_DELTA_R)] = (a[38] + a[39] * 8.0 * M_SQRT2) / (a[40] * 8.0 + PPOW(2.0, a[41])) - 1.0;    // Hurley et al. 2000, eq 17
    rads[static_cast<int>(HURLEY_R_CONSTANT::B_ALPHA_R)] = (a[58] * PPOW(a[66], a[60])) / (a[59] + PPOW(a[66], a[61]));                 // Hurley et al. 2000, eq 21 (wrong in the arxiv version - says = a59*M**(a61))
    rads[static_cast<int>(HURLEY_R_CONSTANT::C_ALPHA_R)] = (a[58] * PPOW(a[67], a[60])) / (a[59] + PPOW(a[67], a[61]));                 // Hurley et al. 2000, eq 21 (wrong in the arxiv version - as above)
    rads[static_cast<int>(HURLEY_R_CONSTANT::B_BETA_R)]  = (a[69] * 8.0 * M_SQRT2) / (a[70] + PPOW(2.0, a[71]));                        // Hurley et al. 2000, eq 22
    rads[static_cast<int>(HURLEY_R_CONSTANT::C_BETA_R)]  = (a[69] * 16384.0) / (a[70] + PPOW(16.0, a[71]));                             // Hurley et al. 2000, eq 22
    
    // return the radius constants vector by value - NRVO takes care of performance/efficiency
    return rads;

#undef a
}


/*
 * CalculateHurleyGammaConstants
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
 * DBL_VECTOR CalculateHurleyGammaConstants(const DBL_VECTOR& p_aCoefficients) const
 * 
 * @param       p_aCoefficients                 Hurley a(n) coefficients
 * @return                                      Gamma constants vector
 */
DBL_VECTOR Globals::CalculateHurleyGammaConstants(const DBL_VECTOR& p_aCoefficients) const {
#define a p_aCoefficients // for convenience and readability - undefined at end of function
               
    // create and inialise gamma constants vector - this is the return value
    DBL_VECTOR gammas(HURLEY_GAMMA_CONSTANT::COUNT, DEFAULT_INITIAL_DOUBLE_VALUE);
        
    // populate gamma constants vector
    double bGamma = max(0.0, a[76] + (a[77] * PPOW((1.0 - a[78]), a[79])));                     // Hurley et al. 2000, eq 23 and discussion following
    gammas[static_cast<int>(HURLEY_GAMMA_CONSTANT::B_GAMMA)] = bGamma;
    gammas[static_cast<int>(HURLEY_GAMMA_CONSTANT::C_GAMMA)] = a[75] <= 1.0 ? bGamma : a[80];   // Hurley et al. 2000, eq 23 and discussion following  
        
    // return the gamma constants vector by value - NRVO takes care of performance/efficiency
    return gammas;

#undef a
}


/*
 * CalculateHurleyMassCutoffs
 *
 * @brief
 * Calculate mass cutoffs per Hurley et al. 2000:
 *
 *   MHook: the mass above which a hook appears on the MS (Hurley et al. 2000, eq 1)
 *   MHeF : the maximum initial mass for which He ignites degenerately in the He Flash (Hurley et al. 2000, eq 2)
 *   MFGB : the maximum mass at which He ignites degenerately on the First Giant Branch (FGB) (Hurley et al. 2000, eq 3)
 *
 * Hurley mass cutoffs depend on the star's metallicity only - so this only needs to be done at most once per star
 * (upon creation), but can also be reused if metallicity doesn't change from one star to the next (e.g. in a population run).
 *
 *
 * DBL_VECTOR CalculateHurleyMassCutoffs(const double p_Z, const double p_Zeta) const
 * 
 * @param       p_Z                             Metallicity
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      Mass cutoffs vector
 */
DBL_VECTOR Globals::CalculateHurleyMassCutoffs(const double p_Z, const double p_Zeta) const {
   
    // create and inialise mass cutoffs vector - this is the return value
    DBL_VECTOR massCutoffs(MASS_CUTOFF::COUNT, DEFAULT_INITIAL_DOUBLE_VALUE);

    // calculate some powers of zeta - for performance and readability
    // this function is only called once per star, and at most twice per binary (but probably once), so not too onerous
    const double zeta2 = p_Zeta * p_Zeta;

    // populate mass cutoffs vector
    massCutoffs[static_cast<int>MASS_CUTOFF::MHook] = 1.0185 + (0.16015 * p_Zeta) + (0.0892 * zeta2);
    massCutoffs[static_cast<int>MASS_CUTOFF::MHeF]  = 1.995 + (0.25 * p_Zeta) + (0.087 * zeta2);  
    massCutoffs[static_cast<int>MASS_CUTOFF::MFGB]  = 13.048 * PPOW((p_Z / ZSOL_HURLEY), 0.06) / (1.0 + (0.0012 * PPOW((ZSOL_HURLEY / p_Z), 1.27)));    

    // return the mass cutoffs vector by value - NRVO takes care of performance/efficiency
    return massCutoffs;
}


/*
 * CalculateHurleyAlphas
 *
 * @brief
 * Calculate Hurley constants alpha1, alpha2, and alpha3
 *
 *   - alpha1 per Hurley et al, 2000, p12, just after eq 49
 *   - alpha3 per Hurley et al, 2000, p12, just after eq 56
 *   - alpha4 per Hurley et al, 2000, p12, just after eq 57
 *
 * alpha1, alpha3, and alpha4 depend on the star's metallicity only - so this only needs to be done once
 * per star (upon creation)
 * 
 * alpha2, per Hurley et al, 2000, p12, just after eq 53, depends on the star's core mass, so is not constant
 * and is not calculated here - but we leave an empty slot in the returned vector so that alpha1, alpha3, and
 * alpha4 can be indexed directly (index is alpha ordinal-1 (e.g. alpha1 is at vec[0]; alpha3 at vec[2], etc.)).
 *
 *
 * DBL_VECTOR CalculateHurleyAlphas(const DBL_VECTOR& p_bCoefficients, const DBL_VECTOR& p_MassCutoffs) const)
 * 
 * @param       p_bCoefficients                 Hurley b(n) coefficients
 * @param       p_MassCutoffs                   Hurley mass cutoffs vector
 * @return                                      Alpha values vector
 */
DBL_VECTOR Globals::CalculateHurleyAlphas(const DBL_VECTOR& p_bCoefficients, const DBL_VECTOR& p_MassCutoffs) const {
// macros for convenience and readability - undefined at end of function
#define massCutoffs(x) p_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)] 
#define b              p_bCoefficients
   
    // create and inialise alphas vector - this is the return value
    DBL_VECTOR alphas(4, DEFAULT_INITIAL_DOUBLE_VALUE);
    
    // alpha1
    const double LHeI_MHeF = (b[11] + (b[12] * PPOW(massCutoffs(MHeF), 3.8))) / (b[13] + (massCutoffs(MHeF) * massCutoffs(MHeF)));
    alphas[0]              = ((p_bCoeffs[9] * PPOW(massCutoffs(MHeF), b[10])) - LHeI_MHeF) / LHeI_MHeF;

    // alpha2 - empty slot at alphas[1] (value will be DEFAULT_INITIAL_DOUBLE_VALUE)

    // alpha3
    const double LBAGB = (b[31] + (b[32] * PPOW(massCutoffs(MHeF), (b[33] + 1.8)))) / (b[34] + PPOW(massCutoffs(MHeF), b[33]));
    alphas[2]          = ((b[29] * PPOW(massCutoffs(MHeF), b[30])) - LBAGB) / LBAGB;

    // alpha4
    const double MHeF5     = massCutoffs(MHeF) * massCutoffs(MHeF) * massCutoffs(MHeF) * massCutoffs(MHeF) * massCutoffs(MHeF); // pow() is slow - use multiplication
    const double tBGB_MHeF = BaseStar::CalculateLifetimeToBGB_Static(massCutoffs(MHeF));                                        // tBGB for mass M = MHeF
    const double tHe_MHeF  = tBGB_MHeF * (b[41] * PPOW(massCutoffs(MHeF), b[42]) + b[43] * massCutoffs(MHeF)) / (b[44] + massCutoffs(MHeF)); 
    alphas[3]              = ((tHe_MHeF - b[39]) / b[39]);

    // return the alphas vector by value - NRVO takes care of performance/efficiency
    return alphas;

#undef b
#undef massCutoffs
}    
    

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            MISCELLANEOUS HURLEY VALUES                            //
//                                                                                   //
/////////////////////////////////////////////////////////////////////////////////////// 

/*
 * CalculateHurleyGBRadiusXexponent
 *
 * @brief
 * Calculate the parameter x for the Giant Branch ('x' exponent to which Radius depends on Mass
 * (at constant Luminosity) - 'x' in Hurley et al. 2000, eq 47)
 * Hybrid of b5 and b7 from Hurley et al. 2000
 * Hurley et al. 2000, eq 47
 *
 * 'x' depends on the star's metallicity only - so this only needs to be done once per star (upon creation),
 * but can also be reused if metallicity doesn't change from one star to the next (e.g. in a population run).
 *
 *
 * double CalculateHurleyGBRadiusXexponent(const double p_Zeta) const
 * 
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      Giant Branch radius 'x' exponent
 */
double Globals::CalculateHurleyGBRadiusXexponent(const double p_Zeta) const {

    // calculate some powers of zeta - for performance and readability
    // this function is only called once per star, and at most twice per binary (but probably once), so not too onerous
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;

    return 0.30406 + (0.0805 * p_Zeta) + (0.0897 * zeta2) + (0.0878 * zeta3) + (0.0222 * zeta4);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//            DRIVER TO (RE)CALCULATE HURLEY METALLICITY-DEPENDENT VALUES            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateHurleyZdependentValues
 *
 * @brief
 * Calculate and set all Hurley Z-dependent values if necessary - if the reference metallicity passed
 * is different from the reference metallicity recorded in the struct parameter passed, then all values
 * will be calculated and set, otherwise the struct passed will be returned unchanged.
 * 
 * Note that the reference metallicity passed to this function will be range checked and clamped
 * to [MINIMUM_METALLICITY_HURLEY, MAXIMUM_METALLICITY_HURLEY] before anything else is done.
 * (This shouldn't happen - options code should already have caught this, but defensive...)
 * 
 * 
 * struct HurleyZdependentValues CalculateHurleyZdependentValues(const double p_Z, const HurleyZdependentValues& p_HurleyZdependentValues) const
 * 
 * @param       p_Z                             Metallicity value to be used calculate Hurley Z-dependent values
 * @param       p_HurleyZdependentValues        Struct containing Hurley Z-dependent values
 * @return                                      Struct containing (possibly updated) Hurley Z-dependent values
 */
HurleyZdependentT CalculateHurleyZdependentValues(const double p_Z, const HurleyZdependentValues& p_HurleyZdependentValues) const {

    HurleyZdependentValues values = p_HurleyZdependentValues;                               // return value - default is unchanged

    double Z = p_Z;

    // range check metallicty passed in and clamp to [MINIMUM_METALLICITY_HURLEY, MAXIMUM_METALLICITY_HURLEY]
    if (Z < MINIMUM_METALLICITY_HURLEY || Z > MAXIMUM_METALLICITY_HURLEY) {                 // reference metallicity outside range?
                                                                                            // yes
        Z = std::max(MAXIMUM_METALLICITY_HURLEY, std::min(Z, MINIMUM_METALLICITY_HURLEY));  // clamp it
                                                                                            // issue warning
        // ISSUE WARNING & CLAMP TO [MINIMUM, MAXIMUM]  
    }

    if (Z != p_HurleyZdependentValues.refZ) {                                               // reference metallicity changed?
                                                                                            // yes
        // (re)calculate and set class member values
        values.refZ                = Z;                                                     // Hurley reference metallicity
        values.sigma               = log10(values.refZ);                                    // Hurley et al. 2000 p24, sigma = log10(Z)
        values.zeta                = values.sigma - LOG10_ZSOL_HURLEY;                      // Hurley et al. 2000 p5, just before eq 1, zeta = log10(Z/0.02) = log10(Z) - log10(0.02) = sigma - LOG10_ZSOL_HURLEY
        values.zetaAnders          = values.sigma - LOG10_ZSOL_ANDERS;                      // log10(Z / ZSOL_ANDERS)
        values.zetaAsplund         = values.sigma - LOG10_ZSOL_ASPLUND;                     // log10(Z / ZSOL_ASPLUND)
        values.rho                 = values.zeta + 1.0;                                     // Hurley et al. 2000 p24, rho = zeta + 1.0

        values.aCoefficients       = CalculateHurleyACoefficients(Z, values.sigma, values.zeta);
        values.massCutoffs         = CalculateHurleyMassCutoffs(Z, values.zeta);
        values.bCoefficients       = CalculateHurleyBCoefficients(Z, values.sigma, values.zeta, values.rho, values.massCutoffs);

        values.gammaConstants      = CalculateHurleyGammaConstants(values.aCoefficients);
        values.luminosityConstants = CalculateHurleyLuminosityConstants(values.aCoefficients);
        values.radiusConstants     = CalculateHurleyRadiusConstants(values.aCoefficients);

        values.xExponent           = CalculateHurleyGBRadiusXexponent();

        values.alphas              = CalculateHurleyAlphas(values.bCoefficients, values.massCutoffs);
    }

    // return the values struct by value - NRVO takes care of performance/efficiency
    return values;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//               (RE)CALCULATE GLOBAL METALLICITY-DEPENDENT VARIABLES                //
//                                                                                   // 
//                                 !!! Be Aware !!!                                  //
//                                                                                   //
// None of the "CalculateAndSet...()" functions below are const or pure - all modify //
// modify members of the Globals class.                                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////
    
    
/*
 * CalculateLuminosityCoefficients_Tout_1996
 *
 * @brief
 * Calculate luminosity coefficients per Tout et al. 1996, table 1
 *
 * Luminosity coefficients depend on the star's metallicity only - so this only needs to be done at most
 * once per star (upon creation), but can also be reused if metallicity doesn't change from one star to
 * the next (e.g. in a population run).
 * 
 *
 * DBL_VECTOR CalculateLuminosityCoefficients_Tout_1996(const double p_Zeta) const
 * 
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      Luminosity coefficients vector
 */
DBL_VECTOR Globals::CalculateLuminosityCoefficients_Tout_1996(const double p_Zeta) const {
// macros for convenience and readability - undefined at end of function
#define index    static_cast<int>(coeff.first)
#define coeff(x) coeff.second[LR_TCoeff::x]
        
    // create and inialise luminosity coefficients vector - this is the return value
    DBL_VECTOR lCoeffs(TOUT_L_COEFF.size(), DEFAULT_INITIAL_DOUBLE_VALUE);
    
    // calculate some powers of zeta - for performance and readability
    // this function is only called once per star, and at most twice per binary (but probably once), so not too onerous
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;
    
    // populate luminosity coefficients vector
    // iterate over luminosity coefficients constants TOUT_L_COEFF (see constants.h)
    // each row is indexed by the TOUT_L_Coeff keys 'ALPHA', 'BETA', 'GAMMA', 'DELTA', 'EPSILON', 'ZETA', 'ETA',
    // and defines the coefficients of the 5 terms (TOUT_LR_TCoeff coefficients 'a', 'b', 'c', 'd', 'e') 
    for (auto coeff: TOUT_L_COEFF) lCoeffs[index] = coeff(a) + (coeff(b) * p_Zeta) + (coeff(c) * zeta2) + (coeff(d) * zeta3) + (coeff(e) * zeta4);
    
    // return the luminosity coefficients vector by value - NRVO takes care of performance/efficiency
    return lCoeffs;
    
#undef coeff
#undef index
}
        
        
/*
 * CalculateRadiusCoefficients_Tout_1996
 *
 * @brief
 * Calculate radius coefficients per Tout et al. 1996, table 2
 *
 * Radius coefficients depend on the star's metallicity only - so this only needs to be done at most
 * once per star (upon creation), but can also be reused if metallicity doesn't change from one star
 * to the next (e.g. in a population run).
 * 
 *
 * DBL_VECTOR CalculateRadiusCoefficients_Tout_1996(const double p_Zeta) const
 * 
 * @param       p_Zeta                          Zeta from Hurley et al. 2000, p5, just before eq 1
 * @return                                      Radius coefficients vector
 */
DBL_VECTOR Globals::CalculateRadiusCoefficients_Tout_1996(const double p_Zeta) const {
// macros for convenience and readability - undefined at end of function
#define index    static_cast<int>(coeff.first)
#define coeff(x) coeff.second[TOUT_LR_TCoeff::x]
        
    // create and inialise radius coefficients vector - this is the return value
    DBL_VECTOR rCoeffs(TOUT_R_COEFF.size(), DEFAULT_INITIAL_DOUBLE_VALUE);
    
    // calculate some powers of zeta - for performance and readability
    // this function is only called once per star, and at most twice per binary (but probably once), so not too onerous
    const double zeta2 = p_Zeta * p_Zeta;
    const double zeta3 = p_Zeta * zeta2;
    const double zeta4 = p_Zeta * zeta3;
    
    // populate radius coefficients vector
    // iterate over radius coefficients constants TOUT_R_COEFF (see constants.h)
    // each row is indexed by the TOUT_R_Coeff keys 'THETA', 'IOTA', 'KAPPA', 'LAMBDA', 'MU', 'NU', 'XI', 'OMICRON', 'PI',
    // and defines the coefficients of the 5 terms (TOUT_LR_TCoeff coefficients 'a', 'b', 'c', 'd', 'e') 
    for (auto coeff: TOUT_R_COEFF) rCoeffs[index] = coeff(a) + (coeff(b) * p_Zeta) + (coeff(c) * zeta2) + (coeff(d) * zeta3) + (coeff(e) * zeta4);
    
    // return the radius coefficients vector by value - NRVO takes care of performance/efficiency
    return rCoeffs;
    
#undef coeff
#undef index
}


/*
 * CalculateZAMSHeliumAbundance_Pols_1998
 *
 * @brief
 * Calculate ZAMS helium abundance as a fraction of the star's mass, per Pols et al. 1998
 *
 *
 * double CalculateZAMSHeliumAbundance_Pols_1998(const double p_Z) const
 * 
 * @param       p_Z                             Metallicity
 * @return                                      ZAMS helium abundance
 */
double Globals::CalculateZAMSHeliumAbundance_Pols_1998(const double p_Z) const {
    return 0.24 + 2.0 * p_Z;
}


/*
 * CalculateZAMSHydrogenAbundance_Pols_1998
 *
 * @brief
 * Calculate ZAMS hydrogen abundance as a fraction of the star's mass, per Pols et al. 1998
 *
 *
 * double CalculateZAMSHydrogenAbundance_Pols_1998(const double p_Z) const
 *
 * @param       p_Z                             Metallicity
 * @return                                      ZAMS hydrogen abundance
 */
double Globals::CalculateZAMSHydrogenAbundance_Pols_1998(const double p_Z) const {
    return 0.76 - 3.0 * p_Z;
}


/*
 * CalculateShikauchiCoefficients
 *
 * Calculate metallicity-dependent coefficients per Shikauchi et al. 2024
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
 * std::tuple<DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> InterpolateShikauchiCoefficients(const double p_Metallicity)
 *
 * @param       p_Z                             Metallicity
 * @return                                      Tuple containing vectors of coefficients (alpha, fMix, luminosity)
 */
std::tuple<DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> Globals::CalculateShikauchiCoefficients(const double p_Z) const {
// macros for convenience and readability - undefined at end of function
#define lower SHIKAUCHI_Coeff::ONE_TENTH_Z_SOL
#define mid   SHIKAUCHI_Coeff::ONE_THIRD_Z_SOL
#define upper SHIKAUCHI_Coeff::Z_SOL

    // create and initialise return values - vectors of coefficients
    DBL_VECTOR alphaCoeff(SHIKAUCHI_FMIX_COEFF[lower].size(), 0.0);             // alpha coefficients
    DBL_VECTOR fMixCoeff(SHIKAUCHI_ALPHA_COEFF[lower].size(), 0.0);             // fMix coefficients
    DBL_VECTOR luminosityCoeff(SHIKAUCHI_LUMINOSITY_COEFF[lower].size(), 0.0);  // luminosity coefficients
       
    // common factors for each of the metallicities defined by Shikauchi et al. 2024
    const double logLowerZ         = std::log10(0.1 * ZSOL_HURLEY);             // Z lower bound: SHIKAUCHI_Coeff::ONE_TENTH_Z_SOL
    const double logMidZ           = std::log10(1.0 / 3.0 * ZSOL_HURLEY);       // Z mid-range:   SHIKAUCHI_Coeff::ONE_THIRD_Z_SOL
    const double logUpperZ         = std::log10(ZSOL_HURLEY);                   // Z upper bound: SHIKAUCHI_Coeff::Z_SOL

    const double logZ              = std::log10(p_Z);
    const double logZ_logLowerZ    = logZ - logLowerZ;
    const double logZ_logMidZ      = logZ - logMidZ;
    const double logMidZ_logZ      = logMidZ  - logZ;
    const double logMidZ_logLowerZ = logMidZ  - logLowerZ;
    const double logUpperZ_logZ    = logUpperZ - logZ;
    const double logUpperZ_logMidZ = logUpperZ - logMidZ;
    
    // calculate coefficients for specified metallicity
    if (logZ <= logLowestZ) {                                                   // p_Z at or below lower metallicity bound?
                                                                                // yes, clamp coefficients to lower bound values
        alphaCoeff      = SHIKAUCHI_ALPHA_COEFF[lower];
        fMixCoeff       = SHIKAUCHI_FMIX_COEFF[lower];
        luminosityCoeff = SHIKAUCHI_LUMINOSITY_COEFF[lower];
    }
    else if (logZ <= middle) {                                                  // p_Z in lower to middle metallicity band?
                                                                                // yes, interpolate
        for (size_t i = 0; i < alphaCoeff.size(); i++)
            alphaCoeff[i] = (SHIKAUCHI_ALPHA_COEFF[lower][i] * logMidZ_logZ + SHIKAUCHI_ALPHA_COEFF[mid][i] * logZ_logLowerZ) / logMidZ_logLowerZ;
        for (size_t i = 0; i < fMixCoeff.size(); i++)
            fMixCoeff[i]  = (SHIKAUCHI_FMIX_COEFF[lower][i] * logMidZ_logZ + SHIKAUCHI_FMIX_COEFF[mid][i] * logZ_logLowerZ) / logMidZ_logLowerZ;
        for (size_t i = 0; i < luminosityCoeff.size; i++)
            luminosityCoeff[i] = (SHIKAUCHI_LUMINOSITY_COEFF[lower][i] * logMidZ_logZ + SHIKAUCHI_LUMINOSITY_COEFF[mid][i] * logZ_logLowerZ) / logMidZ_logLowerZ;
    }
    else if (logZ < high) {                                                     // p_Z in middle to upper metallicity band?
                                                                                // yes, interpolate
            for (size_t i = 0; i < alphaCoeff.size(); i++)
                alphaCoeff[i] = (SHIKAUCHI_ALPHA_COEFF[mid][i] * logUpperZ_logZ + SHIKAUCHI_ALPHA_COEFF[upper][i] * logZ_logMidZ) / logUpperZ_logMidZ;
            for (size_t i = 0; i < fMixCoeff.size(); i++)
                fMixCoeff[i]  = (SHIKAUCHI_FMIX_COEFF[mid][i] * logUpperZ_logZ + SHIKAUCHI_FMIX_COEFF[upper][i] * logZ_logMidZ) / logUpperZ_logMidZ;
            for (size_t i = 0; i < luminosityCoeff.size; i++)
                luminosityCoeff[i] = (SHIKAUCHI_LUMINOSITY_COEFF[mid][i] * logUpperZ_logZ + SHIKAUCHI_LUMINOSITY_COEFF[upper][i] * logZ_logMidZ) / logUpperZ_logMidZ;
    }
    else {                                                                      // p_Z at or above upper metallicity bound
                                                                                // clamp coefficients to upper bound values
        alphaCoeff      = SHIKAUCHI_ALPHA_COEFF[upper];
        fMixCoeff       = SHIKAUCHI_FMIX_COEFF[upper];
        luminosityCoeff = SHIKAUCHI_LUMINOSITY_COEFF[upper];
    }
    
    return std::make_tuple(alphaCoeff, fMixCoeff, luminosityCoeff);

#undef upper
#undef mid
#undef lower
}
