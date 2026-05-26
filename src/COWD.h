#ifndef __COWD_h__
#define __COWD_h__

#include "constants.h"
#include "typedefs.h"

#include "WhiteDwarfs.h"


class BaseStar;
class WhiteDwarfs;

class COWD: virtual public BaseStar, public WhiteDwarfs {

public:

    COWD() { m_StellarType = STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF; };
    
    COWD(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), WhiteDwarfs(p_BaseStar) {
        m_StellarType = STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF;                                                                                                // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                         // Initialise if required
    }


    // member functions

    static double CalculateLuminosityOnPhase_Static(const double p_Mass, const double p_Time) { 
        return WhiteDwarfs::CalculateLuminosityOnPhase_Static(p_Mass, p_Time, WD_Baryon_Number.at(STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF));
    }

protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                                  // Initialise timescales
        m_HShell            = 0.0;                                                                                                                              // Initialize hydrogen shell
        m_HeShell           = 0.0;                                                                                                                              // Initialize helium shell
        m_HeShellDetonation = false;
        m_OffCentreIgnition = false;
        m_AccretionRegime   = ACCRETION_REGIME::ZERO;
        
        EvolveOnPhase(0.0);
    }







    double CalculateLuminosityOnPhase(const double p_Mass, const double p_Time) const {
        return CalculateLuminosityOnPhase_Static(p_Mass, p_Time);
    }
    double          CalculateLuminosityOnPhase() const                              { return CalculateLuminosityOnPhase(m_Mass, m_Age, m_Metallicity); }        // Use class member variables

    
    STELLAR_TYPE    EvolveToNextPhase()                                             { return m_OffCentreIgnition ? STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF : m_HeShellDetonation ? ResolveHeSD() : ResolveSNIa(); }
    bool            IsSupernova() const                                             { return m_HeShellDetonation || IsMassAboveChandrasekhar(); };                                             
    STELLAR_TYPE    ResolveSupernova()                                              { return m_HeShellDetonation ? ResolveHeSD() : ResolveSNIa(); }
    bool            ShouldEvolveOnPhase() const                                     { return m_OffCentreIgnition ? false : !IsSupernova(); };                   // From https://ui.adsabs.harvard.edu/abs/2017MNRAS.472.1593W/abstract around the end of section 3.2. Also, allows SN.                                               

};

#endif // __COWD_h__
