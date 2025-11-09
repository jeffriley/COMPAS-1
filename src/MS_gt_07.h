#ifndef __MS_gt_07_h__
#define __MS_gt_07_h__

#include "constants.h"
//// #include "typedefs.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "profiling.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "utils.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#include "BaseStar.h"
#include "MainSequence.h"

class BaseStar;
class MainSequence;

class MS_gt_07: virtual public BaseStar, public MainSequence {

public:

    // constructors

    MS_gt_07() { m_StellarType = STELLAR_TYPE::MS_GT_07; };
    
    MS_gt_07(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::MS_GT_07; // Set stellar type // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if (p_Initialise) Initialise();                                     // Initialise if required
    }

    // member functions - alphabetically

    COMPAS_PURE ENVELOPE DetermineEnvelopeType() const override;


private:

    void Initialise() {
        CalculateTimescales(); // Initialise timescales // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // Age for MS_GT_07 is carried over from CH stars switching to MS after spinning down, so not set to 0.0 here
        
        //// JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        // Initialise core mass, luminosity, radius, and temperature if Brcek core mass prescription is used
        // Only do this once - this should not be done if a CH star spins down and becomes a MS star (when using CHE_MODE::PESSIMISTIC)
        if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK &&                                                          // Brcek core mass prescription?
            utils::Compare(m_MZAMS, BRCEK_LOWER_MASS_LIMIT) >= 0                         &&                                                             // ZAMS mass >= BRCEK_LOWER_MASS_LIMIT?
            m_Time <= 0.0) {                                                                                                                            // star not yet aged past creation?
                                                                                                                                                        // yes - initialise
            m_InitialMainSequenceCoreMass = MainSequence::CalculateInitialMainSequenceCoreMass(m_MZAMS, m_InitialHeliumAbundance);
            m_MainSequenceCoreMass        = m_InitialMainSequenceCoreMass;
            m_Luminosity                  = MainSequence::CalculateLuminosityOnPhase(m_Age, m_Mass0, BaseStar::CalculateLuminosityAtZAMS(m_Mass0));
            m_Radius                      = MainSequence::CalculateRadiusOnPhase(m_Mass, m_Tau, m_RZAMS0);
            m_Temperature                 = BaseStar::CalculateTemperatureOnPhase(m_Luminosity, m_Radius);
        }
    }

};






/// MS_gt_07_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


class MS_gt_07_Constituent: virtual public BinaryConstituentStar, public MS_gt_07 {

public:


protected:

};


    double      CalculateMTRejuvenationFactor();


double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const ;
double CalculateCriticalMassRatio_Hurley2002() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_GT_07; }


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
inline double MS_gt_07_Constituent::CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {                                                                                                          
    return p_AccretorIsDegenerate
            ? OPTIONS->MassTransferCriticalMassRatioMSHighMassDegenerateAccretor()       // degenerate accretor
            : OPTIONS->MassTransferCriticalMassRatioMSHighMassNonDegenerateAccretor();   // non-degenerate accretor
}


#endif // __MS_gt_07_h__
