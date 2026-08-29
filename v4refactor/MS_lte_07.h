#progma once

#include "constants.h"
#include "BaseStar.h"
#include "MainSequence.h"

class BaseStar;
class MainSequence;

class MS_lte_07: virtual public BaseStar, public MainSequence {

public:

    // Constructors
    MS_lte_07() {}
    MS_lte_07(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) { 
        if (p_Initialise) Initialise();
    }


    private:

    // Initialise
    inline void Initialise() {
        CalculateTimescales(); // Initialise timescales - need to do here for cloning // <<<<<<<<< JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        m_StateHistory.CurrentStateRef().SetAge(0.0); // Set age appropriately - again, for cloning // <<<<<<<<< JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
    //   ENVELOPE                                   //
    //////////////////////////////////////////////////

    GNU_CONST ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::CONVECTIVE; } // Always CONVECTIVE for MS_lte_07 stars

};




/// MS_lte_07_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//



    // MT rejuvenation
    double  CalculateMTRejuvenationFactor() override { return 1.0; }

    // Critical mass ratio
    double  CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const override;       // body inline below
    double  CalculateCriticalMassRatio_Hurley2002() const override { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_LTE_07; }



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

