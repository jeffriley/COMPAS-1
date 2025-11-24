#ifndef __MS_lte_07_h__
#define __MS_lte_07_h__

#include "constants.h"
//// #include "typedefs.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "profiling.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "utils.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#include "BaseStar.h"
#include "MainSequence.h"

class BaseStar;
class MainSequence;

class MS_lte_07: virtual public BaseStar, public MainSequence {

public:

    // constructors

    MS_lte_07() {}
    MS_lte_07(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) { if (p_Initialise) Initialise(); }


    // public member functions - alphabetically

    GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::CONVECTIVE; }   // Always CONVECTIVE for MS_lte_07 stars


private:

    // private member functions - alphabetically

    inline void Initialise() {
        CalculateTimescales(); // Initialise timescales - need to do here for cloning // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        m_StateHistory.CurrentState.SetAge(0.0); // Set age appropriately - again, for cloning // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    }

};









/// MS_lte_07_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


class MS_lte_07_Constituent: virtual public BinaryConstituentStar, public MS_lte_07 {

public:


protected:

};



    double      CalculateMTRejuvenationFactor()                                   { return 1.0; }


double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const ;
double CalculateCriticalMassRatio_Hurley2002() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_LTE_07; }



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
inline double MS_lte_07_Constituent::CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {                                                                                                          
    return p_AccretorIsDegenerate
            ? OPTIONS->MassTransferCriticalMassRatioMSLowMassDegenerateAccretor()       // degenerate accretor
            : OPTIONS->MassTransferCriticalMassRatioMSLowMassNonDegenerateAccretor();   // non-degenerate accretor
}


#endif // __MS_lte_07_h__
