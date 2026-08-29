#pragma once

#include "constants.h"
#include "typedefs.h"

#include "WhiteDwarfs.h"


class BaseStar;
class WhiteDwarfs;

class ONeWD: virtual public BaseStar, public WhiteDwarfs {

public:

    ONeWD() { m_InterimState.SetStellarType(STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF); };
    
    ONeWD(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), WhiteDwarfs(p_BaseStar) {
        m_InterimState.SetStellarType(STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF); // Set stellar type
        if (p_Initialise) Initialise(); // Initialise if required
    }


private:

    void Initialise() {
        m_InterimState.SetAge(0.0);
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
    //   EVOLUTION                                  //
    //////////////////////////////////////////////////
    
    GNU_PURE  STELLAR_TYPE EvolveToNextPhase() const override { return StellarType(); } // Single-star: no AIC, stay ONeWD // JR FIX THIS: DONE
    GNU_CONST bool ShouldEvolveOnPhase() const override { return true; }                // Single-star: no AIC, keep cooling // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateLuminosity() const override { // JR FIX THIS: DONE
        return CalculateLuminosity_Hurley2000(Metallicity(), Mass(), Age(), WD_Baryon_Number.at(STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF));
    }


    //////////////////////////////////////////////////
    //   SUPERNOVAE                                 //
    //////////////////////////////////////////////////
    
    GNU_CONST bool IsSupernova() const override { return false; }                        // Single-star: no AIC // JR FIX THIS: DONE
    GNU_PURE  STELLAR_TYPE ResolveSupernova() const override { return StellarType(); }   // Single-star: no AIC, NO-OP) // JR FIX THIS: DONE
                                                           
};














/// ONeWD_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//




