#ifndef __MR_h__
#define __MR_h__

#include "constants.h"
#include "typedefs.h"
#include "limits.h"

#include "Remnants.h"

class BaseStar;
class Remnants;

class MR: virtual public BaseStar, public Remnants {

public:

    MR() {};
    
    MR(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {
        if (p_Initialise) Initialise();
    }


protected:

    void Initialise() {
        // ensure it's a massless remnant...
        m_InterimState.SetAge(0.0);
        m_InterimState.SetMass(0.0);
        m_InterimState.SetCOCoreMass(0.0);
        m_InterimState.SetHeCoreMass(0.0);
        m_InterimState.SetCoreMass(0.0);
        m_InterimState.SetMassEffectiveInitial(0.0);
        m_InterimState.SetLuminosity(0.0);
        m_InterimState.SetRadius(0.0);
        m_InterimState.SetTemperature(0.0);
    }


    // member functions
   	 double     CalculateMomentOfInertia() const        { return 0.0; }                                     // No moment of inertia for massless remnants - use 0.0
   	 double     CalculateMomentOfInertiaAU() const      { return 0.0; }                                     // No moment of inertia for massless remnants - use 0.0
    
     GNU_CONST inline double  ChooseTimestep() const override { return std::numeric_limits<double>::max(); }// Can take arbitrarily long time steps for massless remnants -- nothing is happening

     void       SetPulsarParameters() const { }                                                             // NO-OP

     bool       ShouldEvolveOnPhase() const             { return true; }                                    // Always
     bool       ShouldSkipPhase() const                 { return false; }                                   // Don't skip
};

#endif // __MR_h__
