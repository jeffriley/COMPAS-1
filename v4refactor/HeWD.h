#pragma once

#include "constants.h"
#include "typedefs.h"

#include "WhiteDwarfs.h"


class BaseStar;
class WhiteDwarfs;

class HeWD: virtual public BaseStar, public WhiteDwarfs {


public:

    HeWD() { m_InterimState.SetStellarType(STELLAR_TYPE::HELIUM_WHITE_DWARF); };
    
    HeWD(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), WhiteDwarfs(p_BaseStar) {
        m_InterimState.SetStellarType(STELLAR_TYPE::HELIUM_WHITE_DWARF);               // Set stellar type
        if (p_Initialise) Initialise();                                 // Initialise if required
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
    //   ABUNDANCE                                  //
    //////////////////////////////////////////////////

    GNU_PURE inline double CalculateHeAbundanceCore(const double p_Tau) const override { return 1.0 - Metallicity(); }
    GNU_PURE inline double CalculateHeAbundanceSurface(const double p_Tau) const override { return 1.0 - Metallicity(); }

    //////////////////////////////////////////////////
    //   EVOLUTION                                  //
    //////////////////////////////////////////////////
    
    GNU_CONST inline bool ShouldEvolveOnPhase() const override { return true; }                 // Single-star: no SNIa or HeSD, keep cooling
    GNU_CONST inline STELLAR_TYPE EvolveToNextPhase() const override { return StellarType(); }  // Single-star: no SNIa or HeSD, stay HeWD


    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateLuminosity_Hurley2000() { // JR FIX THIS: DONE
        return CalculateLuminosity_Hurley2000(Metallicity(), Mass(), Age(), WD_Baryon_Number.at(STELLAR_TYPE::HELIUM_WHITE_DWARF));
    }


    //////////////////////////////////////////////////
    //   SUPERNOVAE                                 //
    //////////////////////////////////////////////////
    
    GNU_CONST inline bool IsSupernova() const override { return false; }                        // Single-star: no SNIa or HeSD
    GNU_CONST inline STELLAR_TYPE ResolveSupernova() const override { return StellarType(); }   // Single-star: no SNIa or HeSD, NO-OP)


    

////    double          CalculateRadiusOnPhase(const double p_Mass) const                                       { return CalculateRadiusOnPhase_Static(p_Mass); }
////    double          CalculateRadiusOnPhase(double p_Mass, double p_Luminosity) const                        { return CalculateRadiusOnPhase(p_Mass); }                                      // ignore luminosity argument for WDs
////    double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhase(m_Mass); }                                      // Use class member variables


    std::tuple <double, STELLAR_TYPE> CalculateRadiusAndStellarTypeOnPhase() const                          { return BaseStar::CalculateRadiusAndStellarTypeOnPhase(); }

};











/// HeWD_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//


