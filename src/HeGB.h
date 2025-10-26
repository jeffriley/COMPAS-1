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

    HeGB* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HeGB* clone = new HeGB(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static HeGB* Clone(HeGB& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        HeGB* clone = new HeGB(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


    // member functions - alphabetically
    static  double  CalculateAgeOnPhase_Static(const double p_Mass, const double p_CoreMass, const double p_tHeMS, const DBL_VECTOR &p_GBparams);



GNU_CONST static double CalculateCoreMass_Hurley2000_Static(const double p_Mass, const double p_Age, const DBL_VECTOR &p_GBparams, const double p_tHeMS);



    static  double  CalculateLuminosityOnPhase_Static(const double p_CoreMass, const double p_GBPB, const double p_GBPD);

    static  DBL_DBL CalculateRadiusOnPhase_Static(const double p_Mass, const double p_Luminosity);


protected:

    void Initialise(const STELLAR_TYPE p_PreviousStellarType) {
        CalculateTimescales();                                                                                                                                          // Initialise timescales
        if (p_PreviousStellarType != STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP)                                                                                   // If not evolving from HeHG...
            m_Age = CalculateAgeOnPhase_Static(m_Mass, m_COCoreMass, m_Timescales[static_cast<int>(TIMESCALE::tHeMS)], m_GBparams);                                 // ... Set age appropriately
        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically




    double      CalculateLuminosityOnPhase(const double p_CoreMass, const double p_GBPB, const double p_GBPD) const { return CalculateLuminosityOnPhase_Static(p_CoreMass, p_GBPB, p_GBPD); }
    double      CalculateLuminosityOnPhase() const                                                                  { return CalculateLuminosityOnPhase(m_CoreMass, m_GBparams[static_cast<int>(GBP::B)], m_GBparams[static_cast<int>(GBP::D)]); }

    double      CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const;
    double      CalculateRadiusOnPhase() const                                                                      { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }

    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase(const double p_Mass, const double p_Luminosity) const;
    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                                  { return CalculateRadiusAndStellarTypeOnPhase(m_Mass, m_Luminosity); }
            
    ENVELOPE    DetermineEnvelopeType() const                                                                       { return ENVELOPE::CONVECTIVE; }                        // Always CONVECTIVE    
};






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
