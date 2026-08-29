#pragma once

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "HeHG.h"


class BaseStar;
class HeHG;

class HeGB: virtual public BaseStar, public HeHG {



public:

    HeGB() { m_InterimState.SetStellarType(STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH); };
    
    HeGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeHG(p_BaseStar, false) {
        STELLAR_TYPE currentStellarType = StellarType();
        m_InterimState.SetStellarType(STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH); // Set stellar type
        if (p_Initialise) Initialise(currentStellarType); // Initialise if required
    }



private:

    void Initialise(const STELLAR_TYPE p_PreviousStellarType) {
        // If we're using one of the HURLEY evolution modes, and the star is
        // not evolving from HeHG then we need to recalculate the age of the star.
        // (Currently this is only done for Hurley evolution modes - revisit if also true for other evolution modes)
        if (p_PreviousStellarType != STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP && 
           (OPTIONS->Mode() == EVOLUTION_MODE::SSE_HURLEY || OPTIONS->Mode() == EVOLUTION_MODE::BSE_HURLEY)) {
            m_InterimState.SetAge(CalculateAge_Hurley2000(Mass(), COCoreMass(), GBParams(), Timescales(TS::HeMS)));
        }
    }
    
    
protected:

    // Member functions (not getters or setters)
    //
    // VIRTUAL FUNCTIONS may be (are expected to be) overridden by derived classes.
    // When overriding virtual functions in a derived class, use the "override" attribute.
    //
    // NON-VIRTUAL FUNCTIONS should not be overridden (declared separately) by derived classes.
    // While it is legal in C++ to declare the same (non-virtual) function in multiple classes,
    // (aka "shadowing", or "hiding"), we discourage it.  Non-virtual functions are statically
    // bound, and as such, especially with indirection, may not produce expected results.

    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////

    GNU_CONST static double CalculateAge_Hurley2000(const double p_Mass, const double p_CoreMass, const GBParamsT& p_GBParams, const double p_tHeMS); // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   ENVELOPE                                   //
    //////////////////////////////////////////////////

    GNU_CONST ENVELOPE DetermineEnvelopeType() const override { // JR FIX THIS: DONE
        return ENVELOPE::CONVECTIVE; // Always CONVECTIVE for HeGB stars
    }


    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    GNU_CONST static double CalculateLuminosity_Hurley2000(const double p_CoreMass, const GBParamsT& p_GBParams); // JR FIX THIS: DONE
    GNU_PURE double CalculateLuminosity_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateLuminosity_Hurley2000(CoreMass(), GBParams());
    }


    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////

    GNU_CONST static double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Age, const GBParamsT& p_GBParams, const double p_tHeMS); // JR FIX THIS: DONE
    GNU_PURE double CalculateCoreMass_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateCoreMass_Hurley2000(Mass(), Age(), GBParams(), Timescales(TS::tHeMS));
    }


    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////

    GNU_CONST static double CalculateRadius_Hurley2000(const double p_Mass, const double p_Luminosity); // JR FIX THIS: DONE
    double CalculateRadius_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateRadius_Hurley2000(Mass(), Luminosity());
    }


    GNU_CONST static double CalculateRadius_Hurley2000_HeHayashi(const double p_Luminosity); // JR FIX THIS: DONE

};


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateAge_Hurley2000
 *
 * @brief
 * Calculate the age of an evolved Helium Giant Branch, HeGB, star,
 * per Hurley et al. 2000, eq 39 (using eqs 40-42), modified as described after eq 84
 *
 *
 * static double CalculateAge_Hurley2000(const double p_Mass, const double p_CoreMass, const GBParamsT& p_GBParams, const double p_tHeMS)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_GBParams                      GB parameters
 * @param       p_tHeMS                         HeMs lifetime (Myr)
 * @return                                      HeGB age (Myr)
 */
double HeGB::CalculateAge_Hurley2000(const double p_Mass, const double p_CoreMass, const GBParamsT& p_GBParams, const double p_tHeMS) {

    double age = 0.0;                                                               // Default return value

    const double AHe = p_GBParams[GBP::AHe];
    const double B   = p_GBParams[GBP::B];
    const double D   = p_GBParams[GBP::D];
    const double p   = p_GBParams[GBP::P];
    const double q   = p_GBParams[GBP::Q];
    const double Lx  = p_GBParams[GBP::LX];
    const double Mx  = p_GBParams[GBP::MX];

    const double p1    = p - 1.0;
    const double p1_p  = p1 / p;
    const double ltHe  = HeMS::CalculateLuminosityAtPhaseEnd_Hurley2000(p_Mass);
    const double tinf1 = p_tHeMS + ((1.0 / (p1 * AHe * D)) * PPOW(D / ltHe, p1_p)); // Eq 40, modified

    if (p_CoreMass > Mx) {
        const double q1    = q - 1.0;
        const double tx    = tinf1 - ((tinf1 - p_tHeMS) * PPOW(ltHe / Lx, p1_p));   // Eq 41, modified
        const double tinf2 = tx + ((1.0 / (q1 * AHe * B)) * PPOW(B / Lx, q1 / q));  // Eq 42, modfied

        age = tinf2 - (PPOW(p_CoreMass, 1.0 - q) / (q1 * AHe * B));                 // Eq 39, modified                 
    }
    else {
        age = tinf1 - (PPOW(p_CoreMass, 1.0 - p) / (p1 * AHe * D));                 // Eq 39, modified 
    }

    return age;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLuminosity_Hurley2000
 *
 * @brief
 * Calculate luminosity on the Helium Giant Branch (HeGB), per Hurley et al. 2000, eq 84
 *
 *
 * static double CalculateLuminosity_Hurley2000(const double p_CoreMass, const GBParamsT& p_GBParams)
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_GBParams                      GB parameters
 * @return                                      HeGB (post-HeMS) luminosity (Lsol)
 *
 * p_GBPB and p_GBPD passed as parameters so function can be declared static
 */
double HeGB::CalculateLuminosity_Hurley2000(const double p_CoreMass, const GBParamsT& p_GBParams) {
    const double Mc3 = p_CoreMass * p_CoreMass * p_CoreMass;
    return std::min((p_GBParams[GBP::B] * Mc3), (p_GBParams[GBP::D] * p_CoreMass * p_CoreMass * Mc3));
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateCoreMass_Hurley2000
 *
 * @brief
 * Calculate the core mass on the Helium Giant Branch, HeGB,
 * per Hurley et al. 2000, eq 39 (using eqs 40-42), modified as described after eq 84
 *
 *
 * static double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Age, const GBParamsT& p_GBParams, const double p_tHeMS)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_GBParams                      GB parameters
 * @param       p_tHeMS                         HeMs lifetime (Myr)
 * @return                                      HeGB core mass (Msol)
 */
double HeGB::CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Age, const GBParamsT& p_GBParams, const double p_tHeMS) {

    double coreMass = 0.0;                                                          // Default return value

    const double AHe = p_GBParams[GBP::AHe];
    const double B   = p_GBParams[GBP::B];
    const double D   = p_GBParams[GBP::D];
    const double p   = p_GBParams[GBP::P];
    const double q   = p_GBParams[GBP::Q];
    const double Lx  = p_GBParams[GBP::LX];

    const double p1    = p - 1.0;
    const double p1_p  = p1 / p;
    const double ltHe  = HeMS::CalculateLuminosityAtPhaseEnd_Hurley2000(p_Mass);
    const double tinf1 = p_tHeMS + ((1.0 / (p1 * AHe * D)) * PPOW(D / ltHe, p1_p)); // Eq 40, modified
    const double tx    = tinf1 - (tinf1 - p_tHeMS) * PPOW((ltHe / Lx), p1_p);       // Eq 41, modified
    
    if (p_Age > tx) {
        const double q1    = q - 1.0;
        const double tinf2 = tx + ((1.0 / (q1 * AHe * B)) * PPOW(B / Lx, q1 / q));  // Eq 42, modfied

        coreMass = PPOW(q1 * AHe * B * (tinf2 - p_Age), 1.0 / (1.0 - q));           // Eq 39, modfied
    }
    else {
        coreMass = PPOW(p1 * AHe * D * (tinf1 - p_Age), 1.0 / (1.0 - p));           // Eq 39, modfied
    }

    return coreMass;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateRadius_Hurley2000
 *
 * @brief
 * Calculate the Helium Giant Branch, HeGB, radius, per Hurley et al. 2000
 * 
 * Per discussion around eqs 85-88, returns the minimum of eq 86 (using eqs 85 & 87),
 * and eq 88 (via HeGB::CalculateRadius_Hurley2000_HeHayashi())
 *
 * 
 * Note:
 * May need to change this according to section 2.1.2 of http://iopscience.iop.org/article/10.1086/340304/pdf 
 * which talks about updated helium star evolution.
 * Or replace with helium star tracks from binary_c.
 * Especially important for low mass helium stars, BNS progenitors.
 * This paper suggest mass below which envelope is convective is Mconv = 4.5 Msol, they leave it as an uncertain.
 * Rapid expansion given by the second term in eq 86.
 *
 *
 * static double CalculateRadius_Hurley2000(const double p_Mass, const double p_Luminosity)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      HeGB (post-HeMS) radius (Rsol)
 */
inline double HeGB::CalculateRadius_Hurley2000(const double p_Mass, const double p_Luminosity) {

    // Sanity check for mass and luminosity - just return 0.0 if mass or luminosity <= 0
    if (!(p_Mass > 0.0) || !(p_Luminosity > 0.0)) return 0.0;

    const double rZHe = HeMS::CalculateRadiusAtZAHeMS_Hurley2000(p_Mass);
    const double lTHe = HeMS::CalculateLuminosityAtPhaseEnd_Hurley2000(p_Mass);
    
    const double m2   = p_Mass * p_Mass;            // pow() is slow - use multiplication
    const double m2_5 = m2 * std::sqrt(p_Mass);     // srqt() is much faster than pow()

    const double lamda = 500.0 * (2.0 + (m2_5 * m2_5)) / m2_5;

    const double r1 = rZHe * PPOW((p_Luminosity / lTHe), 0.2) + (0.02 * (std::exp(p_Luminosity / lamda) - std::exp(lTHe / lamda)));
    const double r2 = CalculateRadius_Hurley2000_HeHayashi(p_Luminosity);

    return std::min(r1, r2);
}


/*
 * CalculateRadius_Hurley2000_HeHayashi
 *
 * @brief
 * Calculate the giant branch radius for a helium star, per Hurley et al. 2000 eq 88
 * (mimics the Hayashi track)
 * 
 * static double CalculateRadius_Hurley2000_HeHayashi(const double p_Luminosity)
 * 
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      HeGB (post-HeMS) radius (mimicking the Hayashi track) (Rsol)
 */
inline double HeGB::CalculateRadius_Hurley2000_HeHayashi(const double p_Luminosity) {
    return 0.08 * PPOW(p_Luminosity, 0.75);
}












/// HeGB_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    double  CalculateZetaEquilibrium() override { return 0.0; }                                                  // At lowest order, giants with a convective envelope have radii that are insensitive to mass loss

    double  CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const override;            // body inline below
    double  CalculateCriticalMassRatio_Hurley2002() const override { return HURLEY_HJELLMING_WEBBINK_QCRIT_HE_GIANT; }



/*
 * CalculateCriticalMassRatio_Claeys2014
 *
 * @brief
 * Calculate the critical mass ratio, per Claeys et al. 2014
 * 
 * @param       p_AccretorIsDegenerate          Boolean indicating if accretor is degenerate
 * @return                                      Critical mass ratio for unstable MT 
 */
inline double HeGB_Constituent::CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {                                                                                                          
    return p_AccretorIsDegenerate
            ? OPTIONS->MassTransferCriticalMassRatioHeliumGiantDegenerateAccretor()     // degenerate accretor
            : OPTIONS->MassTransferCriticalMassRatioHeliumGiantNonDegenerateAccretor(); // non-degenerate accretor
}

