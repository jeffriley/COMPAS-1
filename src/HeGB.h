#ifndef __HeGB_h__
#define __HeGB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "HeHG.h"


class BaseStar;
class HeHG;

class HeGB: virtual public BaseStar, public HeHG {

public:

    HeGB() { m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH; };
    
    HeGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HeHG(p_BaseStar, false) {
        STELLAR_TYPE currentStellarType = m_StellarType;                                                                                                                // Stellar type evolving from
        m_StellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH;                                                                                                   // Set stellar type
        if (p_Initialise) Initialise(currentStellarType);                                                                                                               // Initialise if required
    }


    // member functions - alphabetically

    GNU_CONST static double CalculateAge_Static(const double p_Mass, const double p_CoreMass, const double p_tHeMS, const DBL_VECTOR& p_GBparams);

    GNU_CONST static double CalculateCoreMass_Hurley2000_Static(const double p_Mass, const double p_Age, const DBL_VECTOR& p_GBparams, const double p_tHeMS);

    inline double CalculateLuminosity_Hurley2000() const override {
        return CalculateLuminosity_Hurley2000_Static(
            m_StateHistory.CurrentState.CoreMass(),
            m_StateHistory.CurrentState.GBparams(static_cast<int>(HURLEY_GBP::B)),
            m_StateHistory.CurrentState.GBparams(static_cast<int>(HURLEY_GBP::D))
        );
    }

    GNU_CONST static double CalculateLuminosity_Hurley2000_Static(const double p_CoreMass, const double p_GBparamB, const double p_GBparamD);

    inline double CalculateRadius_Hurley2000() const override {
        return CalculateRadius_Hurley2000_Static(m_StateHistory.CurrentState.Mass(), m_StateHistory.CurrentState.Luminosity());
    }

    GNU_CONST static double CalculateRadius_Hurley2000_HeHayashi_Static(const double p_Luminosity);

    GNU_CONST static double CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Luminosity);

    GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::CONVECTIVE; } // Always CONVECTIVE for HeGB stars


private:

    void Initialise(const STELLAR_TYPE p_PreviousStellarType) {
        CalculateTimescales();                                                                                                                                          // Initialise timescales
        if (p_PreviousStellarType != STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP)                                                                                   // If not evolving from HeHG...
            m_Age = CalculateAgeOnPhase_Static(m_Mass, m_COCoreMass, m_Timescales[static_cast<int>(TIMESCALE::tHeMS)], m_GBparams);                                 // ... Set age appropriately
        EvolveOnPhase(0.0);
    }           
};


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLuminosity_Hurley2000_Static
 *
 * @brief
 * Calculate luminosity on the Helium Giant Branch (HeGB), per Hurley et al. 2000, eq 84
 *
 *
 * static double CalculateLuminosity_Hurley2000_Static(const double p_CoreMass, const double p_GBparamB, const double p_GBparamD)
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_GBparamB                      Hurley GB parameter B
 * @param       p_GBparamD                      Hurley GB parameter D
 * @return                                      HeGB (post-HeMS) luminosity (Lsol)
 *
 * p_GBPB and p_GBPD passed as parameters so function can be declared static
 */
GNU_CONST inline double HeGB::CalculateLuminosity_Hurley2000_Static(const double p_CoreMass, const double p_GBparamB, const double p_GBparamD) {
    const double Mc3 = p_CoreMass * p_CoreMass * p_CoreMass;
    return std::min((p_GBparamB * Mc3), (p_GBparamD * p_CoreMass * p_CoreMass * Mc3));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRadius_Hurley2000_HeHayashi_Static
 *
 * @brief
 * Calculate the giant branch radius for a helium star, per Hurley et al. 2000 eq 88
 * (mimics the Hayashi track)
 * 
 * static double CalculateRadius_Hurley2000_HeHayashi_Static(const double p_Luminosity)
 * 
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      HeGB (post-HeMS) radius (mimicking the Hayashi track) (Rsol)
 */
GNU_CONST inline double HeGB::CalculateRadius_Hurley2000_HeHayashi_Static(const double p_Luminosity) {
    return = 0.08 * PPOW(p_Luminosity, 0.75);
}






/// HeGB_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    

class HeGB_Constituent: virtual public BinaryConstituentStar, public HeGB {

public:


protected:

};


    double          CalculateZetaEquilibrium()                                                                      { return 0.0; }                                         // At lowest order, giants with a convective envelope have radii that are insensitive to mass loss 

double      CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const;

double CalculateCriticalMassRatio_Hurley2002() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_HE_GIANT; }




/*
 * CalculateCriticalMassRatio_Claeys2014
 *
 * @brief
 * Calculate the critical mass ratio, per Claeys et al. 2014
 * 
 * The critical mass ratio indicates whether the mass transfer is unstable.
 *
 * See Claeys et al. 2014, de Mink et al. 2013, and Ge et al. 2010, 2015, and 2020, for discussions.
 *
 * Assumes this star is the donor.
 * 
 * Critical mass ratio is defined as qCrit = Maccretor / Mdonor.
 * 
 *
 * double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const
 *
 * @param       p_AccretorIsDegenerate          Boolean indicating if accretor is degenerate
 * @return                                      Critical mass ratio for unstable MT 
 */
inline double HeGB_Constituent::CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {                                                                                                          
    return p_AccretorIsDegenerate
            ? OPTIONS->MassTransferCriticalMassRatioHeGiantDegenerateAccretor()     // degenerate accretor
            : OPTIONS->MassTransferCriticalMassRatioHeGiantNonDegenerateAccretor(); // non-degenerate accretor
}


#endif // __HeGB_h__
