#include "HeGB.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         AGE / LIFETIME / TAU / TIMESCALES                         //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateAge_Static
 *
 * @brief
 * Calculate the age of an evolved Helium Giant Branch, HeGB, star,
 * per Hurley et al. 2000, eq 39 (using eqs 40-42), modified as described after eq 84
 *
 *
 * static double CalculateAge_Static(const double p_Mass, const double p_CoreMass, const double p_tHeMS, const DBL_VECTOR& p_GBparams)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_tHeMS                         HeMs lifetime (per Hurley timescales) (Myr)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      HeGB age (Myr)
 */
GNU_CONST double HeGB::CalculateAge_Static(const double p_Mass, const double p_CoreMass, const double p_tHeMS, const DBL_VECTOR& p_GBparams) {
#define GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)] // for convenience and readability - undefined at end of function

    double age;

    const double p1    = GBparams(p) - 1.0;
    const double p1_p  = p1 / GBparams(p);

    const double ltHe  = HeMS::CalculateLuminosityAtPhaseEnd_Hurley2000_Static(p_Mass);
    const double tinf1 = p_tHeMS + ((1.0 / (p1 * GBparams(AHe) * GBparams(D))) * PPOW(GBparams(D) / ltHe, p1_p));                       // eq 40, modified

    if (p_CoreMass > GBparams(Mx)) {
        const double q1    = GBparams(q) - 1.0;
        const double tx    = tinf1 - ((tinf1 - p_tHeMS) * PPOW(ltHe / GBparams(Lx), p1_p));                                             // eq 41, modified
        const double tinf2 = tx + ((1.0 / (q1 * GBparams(AHe) * GBparams(B))) * PPOW(GBparams(B) / GBparams(Lx), q1 / GBparams(q)));    // eq 42, modfied

        age = tinf2 - (PPOW(p_CoreMass, 1.0 - GBparams(q)) / (q1 * GBparams(AHe) * GBparams(B)));                                       // eq 39, modified                 
    }
    else {
        age = tinf1 - (PPOW(p_CoreMass, 1.0 - GBparams(p)) / (p1 * GBparams(AHe) * GBparams(D)));                                       // eq 39, modified 
    }

    return age;

#undef GBparams
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCoreMass_Hurley2000_Static
 *
 * @brief
 * Calculate the core mass on the Helium Giant Branch, HeGB,
 * per Hurley et al. 2000, eq 39 (using eqs 40-42), modified as described after eq 84
 *
 *
 * static double CalculateCoreMass_Hurley2000_Static(const double p_Mass, const double p_Age, const DBL_VECTOR& p_GBparams, const double p_tHeMS)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tHeMS                         HeMs lifetime (per Hurley timescales) (Myr)
 * @return                                      HeGB core mass (Msol)
 */
GNU_CONST double HeGB::CalculateCoreMass_Hurley2000_Static(const double p_Mass, const double p_Age, const DBL_VECTOR& p_GBparams, const double p_tHeMS) {
#define GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)] // for convenience and readability - undefined at end of function

    double coreMass;

    const double p1    = GBparams(p) - 1.0;
    const double p1_p  = p1 / GBparams(p);

    const double ltHe  = HeMS::CalculateLuminosityAtPhaseEnd_Hurley2000_Static(p_Mass);
    const double tinf1 = p_tHeMS + ((1.0 / (p1 * GBparams(AHe) * GBparams(D))) * PPOW(GBparams(D) / ltHe, p1_p));                       // eq 40, modified
    const double tx    = tinf1 - (tinf1 - p_tHeMS) * PPOW((ltHe / GBparams(Lx)), p1_p);                                                 // eq 41, modified
    
    if (p_Age > tx) {
        const double q1    = GBparams(q) - 1.0;
        const double tinf2 = tx + ((1.0 / (q1 * GBparams(AHe) * GBparams(B))) * PPOW(GBparams(B) / GBparams(Lx), q1 / GBparams(q)));    // eq 42, modfied

        coreMass = PPOW(q1 * GBparams(AHe) * GBparams(B) * (tinf2 - p_Age), 1.0 / (1.0 - GBparams(q)));                                 // eq 39, modfied
    }
    else {
        coreMass = PPOW(p1 * GBparams(AHe) * GBparams(D) * (tinf1 - p_Age), 1.0 / (1.0 - GBparams(p)));                                 // eq 39, modfied
    }

    return coreMass;

#undef GBparams
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRadius_Hurley2000_Static
 *
 * @brief
 * Calculate the Helium Giant Branch, HeGB, radius, per Hurley et al. 2000
 * 
 * Per discussion around eqs 85-88, returns the minimum of eq 86 (using eqs 85 & 87),
 * and eq 88 (via HeGB::CalculateRadius_Hurley2000_HeHayashi_Static())
 *
 * 
 * May need to change this according to section 2.1.2 of http://iopscience.iop.org/article/10.1086/340304/pdf 
 * which talks about updated helium star evolution
 * Or replace with helium star tracks from binary_c
 * Especially important for low mass helium stars, BNS progenitors
 * This paper suggest mass below which envelope is convective is Mconv = 4.5 Msol, they leave it as an uncertain
 * Rapid expansion given by the second term in eq 86
 *
 *
 * static double CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Luminosity)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      HeGB (post-HeMS) radius (Rsol)
 */
GNU_CONST double HeGB::CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Luminosity) {

    // sanity check for mass and luminosity - just return 0.0 if mass or luminosity <= 0
    if (p_Mass <= 0.0 || p_Luminosity <= 0.0) return 0.0;

    const double rZHe = HeMS::CalculateRadiusAtZAHeMS_Hurley2000_Static(p_Mass);
    const double lTHe = HeMS::CalculateLuminosityAtPhaseEnd_Hurley2000_Static(p_Mass);
    
    const double m2   = p_Mass * p_Mass;            // pow() is slow - use multiplication
    const double m2_5 = m_2 * std::sqrt(p_Mass);    // srqt() is much faster than pow()

    const double lamda = 500.0 * (2.0 + (m_2_5 * m_2_5)) / m2_5;

    const double r1 = rZHe * PPOW((p_Luminosity / lTHe), 0.2) + (0.02 * (exp(p_Luminosity / lamda) - exp(lTHe / lamda)));
    const double r2 = CalculateRadius_Hurley2000_HeHayashi_Static(p_Luminosity);

    return std::min(r1, r2);
}
