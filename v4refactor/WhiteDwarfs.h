#pragma once

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Remnants.h"
#include "ConstituentStar.h"   // for WhiteDwarfs_Constituent (bottom of file)


class BaseStar;
class Remnants;

class WhiteDwarfs: virtual public BaseStar, public Remnants {

public:

    WhiteDwarfs(){};

    WhiteDwarfs(const BaseStar &p_BaseStar) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {}





    Dbl_DblT         CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich);
    Dbl_DblT         CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate, const bool p_IsHeRich) { return CalculateMassAcceptanceRate(p_DonorMassRate, p_IsHeRich); }

    double           HeShell()        const override { return m_HeShell; }
    double           HShell()         const override { return m_HShell; }
    double           L0Ritter()       const override { return m_L0Ritter; }
    double           LambdaRitter()   const override { return m_LambdaRitter; }
    double           XRitter()        const override { return m_XRitter; }
    void             SetHeShell(const double p_Value)      override { m_HeShell = p_Value; }
    void             SetHShell(const double p_Value)       override { m_HShell = p_Value; }
    double           CalculateEtaH(const double p_MassIntakeRate)  override;        // bodies in WhiteDwarfs.cpp; promoted from protected to public
    double           CalculateEtaHe(const double p_MassIntakeRate) override;


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

    GNU_CONST ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::CONVECTIVE; }  // Always CONVECTIVE for White Dwarfs // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    COMPAS_PURE static double CalculateLuminosity_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const double p_WDBaryonNumber); // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateCOCoreMass() const override { return COCoreMass(); }               // McCO constant for WDs // JR FIX THIS: DONE
    GNU_PURE double CalculateCOCoreMassAtPhaseEnd() const override { return COCoreMass(); }     // McCO constant for WDs // JR FIX THIS: DONE


    GNU_PURE double CalculateHeCoreMass() const override { return HeCoreMass(); }               // McHe constant for WDs // JR FIX THIS: DONE
    GNU_PURE double CalculateHeCoreMassAtPhaseEnd() const override { return HeCoreMass(); }     // McHe constant for WDs // JR FIX THIS: DONE


    GNU_PURE bool   IsMassAboveChandrasekhar() const override { return Mass() > MCH; } // JR FIX THIS: DOME


    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////

    GNU_PURE  double CalculateRadius() const override { return CalculateRadius_Marsh2004(Mass()); } // JR FIX THIS: DONE


    GNU_CONST double CalculateRadius_Marsh2004(const double p_Mass) const; // JR FIX THIS: DONE


    







            


            double           CalculateEtaPTY(const double p_MassIntakeRate);

            double           Calculatel0Ritter() const { return (utils::Compare(Metallicity(), 0.01) > 0) ? L0_RITTER_HIGH_Z : L0_RITTER_LOW_Z; }

            double           CalculateXRitter() const { return utils::Compare(Metallicity(), 0.01) > 0 ? 0.7 : 0.8; }

            double           CalculateLambdaRitter() const { return utils::Compare(Metallicity(), 0.01) > 0 ? 8.0 : 5.0; }

            double           CalculateInitialSupernovaMass() const { return 0.0; }





















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
 * CalculateLuminosity_Hurley2000
 *
 * @brief
 * Calculate the luminosity of a White Dwarf as it cools, per Hurley et al. 2000, eq 90
 *
 *
 * static double CalculateLuminosity_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const double p_WDBaryonNumber)
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Age since White Dwarf formation (Myr)
 * @param       p_WDBaryonNumber                Whitedwarf effective baryon number
 * @return                                      Luminosity of the star (Lsol)
 */
inline double WhiteDwarfs::CalculateLuminosity_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const double p_WDBaryonNumber) {
    return (635.0 * p_Mass * PPOW(p_Metallicity, 0.4)) / PPOW(p_WDBaryonNumber * (p_Age + 0.1), 1.4);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateRadius_Marsh2004
 *
 * @brief
 * Calculate the radius of a white dwarf, per Marsh et al. 2004 eq 24.
 * (See https://academic.oup.com/mnras/article/350/1/113/986306)
 *
 * Originally from Eggleton 1986, quoted in Verbunt & Rappaport 1988, equation
 * from Marsh et al. 2004.
 * 
 * Compared to the Hurley et al. 2000 prescription, the additional factor that
 * includes WD_MP allows for the change to a constant density configuration at
 * low masses (e.g., Zapolsky & Salpeter 1969) after mass loss episodes.
 * 
 * Since a WD is ~Earth-size, expect the returned WD mass to be ~0.009 for 
 * p_Mass > 0.0.  If p_Mass is above the Chandrasekhar mass, the returned WD
 * radius will be the radius of a nuetron star (NEUTRON_STAR_RADIUS).  Zero
 * (or below!) mass means zero radius.
 *
 * 
 * double CalculateRadius_Marsh2004(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Radius of the star (Rsol)
 */
inline double WhiteDwarfs::CalculateRadius_Marsh2004(const double p_Mass) const {

    double radius = 0.0;                    // Default return value - only if mass <= 0.0

    if (p_Mass > 0.0) {                     // Mass of star > 0.0?
        if (p_Mass >= MCH) {                // Yes - mass > Chandrasekhar mass?
            radius = NEUTRON_STAR_RADIUS;   // Yes - giant star core/remnant
        }
        else {                              // Use Marsh et al. 2004 eq 24
            const double MCH_Mass_1_3 = std::cbrt(MCH / p_Mass); 
            const double MCH_Mass_2_3 = MCH_Mass_1_3 * MCH_Mass_1_3;
            const double MP_Mass      = WD_MP / p_Mass;
            const double MP_Mass_2_3  = MP_Mass / std::cbrt(WD_MP / p_Mass); 
            const double f            = 1.0 + 3.5 * MP_Mass_2_3 + MP_Mass;

            radius = std::max(NEUTRON_STAR_RADIUS, 0.0114 * std::sqrt((MCH_Mass_2_3 - 1.0 / MCH_Mass_2_3)) * std::cbrt(f) / f);
        }
    }

    return radius;
}







/// WhiteDwarfs_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//
class WhiteDwarfs_Constituent: public Remnants_Constituent {

public:






            STELLAR_TYPE     ResolveAIC()  override;  
            STELLAR_TYPE     ResolveSNIa() override;  
            STELLAR_TYPE     ResolveHeSD() override;  
            //STELLAR_TYPE     ResolveSupernova()                                                                 { return EvolveToNextPhase(); }                                 // SNe for WDs are handled internally to each WD type




    // member variables
    //
 
            double           m_HeShell;                                                                                                                                         // Current WD He-shell size (Msol). Increases through accretion.
            double           m_HShell;                                                                                                                                          // Current WD H-shell size (Msol). Increases through accretion.
            double           m_L0Ritter;                                                                                                                                        // Parameter from numerical calculations, see Ritter 1999, section 3. Eqs 10 and 12, as well as table 2. Corresponds to L0.
            double           m_LambdaRitter;                                                                                                                                    // Parameter from numerical calculations, see Ritter 1999, section 3. Eqs 10 and 12, as well as table 2.  Cached metallicity constant (not binary state).
            double           m_XRitter;                                                                                                                                         // Assumed hydrogen-mass fraction of material being accreted by He WD, as in Ritter 1999, table 2.
            
            // member functions - alphabetically
            double           CalculateAccretionRegime(const bool   p_DonorIsHeRich,
                                                      const bool   p_DonorIsGiant,
                                                      const double p_DonorThermalMassLossRate,
                                                      const double p_MassLostByDonor);



double CalculateCriticalMassRatio(const double p_Mass,
                                  const double p_Radius,
                                  const double p_CoreMass,
                                  const bool   p_AccretorIsDegenerate,
                                  const double p_MTefficiency = 0.0) const { return CalculateCriticalMassRatio_Hurley2002(); }

double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {
    return CalculateCriticalMassRatio_Hurley2002();
}

double CalculateCriticalMassRatio_Ge2020(const double p_MTefficiency) {
    return CalculateCriticalMassRatio_Hurley2002();
}

double CalculateCriticalMassRatio_Hurley2002() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_WD; }






















    MT_CASE DetermineMassTransferTypeAsDonor() const override { return MT_CASE::OTHER; }                       // Not A, B, C, or NONE

    double  CalculateCriticalMassRatio(const double p_Mass,
                                       const double p_Radius,
                                       const double p_CoreMass,
                                       const bool   p_AccretorIsDegenerate,
                                       const double p_MTefficiency = 0.0) const override {
        (void)p_Mass; (void)p_Radius; (void)p_CoreMass; (void)p_AccretorIsDegenerate; (void)p_MTefficiency;
        return CalculateCriticalMassRatio_Hurley2002();
    }
    double  CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const override { (void)p_AccretorIsDegenerate; return CalculateCriticalMassRatio_Hurley2002(); }
    double  CalculateCriticalMassRatio_Ge2020(const double p_MTefficiency) override { (void)p_MTefficiency; return CalculateCriticalMassRatio_Hurley2002(); }
    double  CalculateCriticalMassRatio_Hurley2002() const override { return HURLEY_HJELLMING_WEBBINK_QCRIT_WD; }

    double  CalculateLambdaRitter() const override;

    ACCRETION_REGIME WhiteDwarfAccretionRegime() const override { return CurrentState().AccretionRegime(); }
    ACCRETION_REGIME DetermineAccretionRegime(const double p_DonorThermalMassLossRate, const bool p_HeRich) override;     // body in WhiteDwarfs.cpp
    void             ResolveShellChange(const double p_AccretedMass) override;                                            // body in WhiteDwarfs.cpp

    Dbl_DblT CalculateMassAcceptanceRate(const double p_DonorMassRate,
                                         const double p_AccretorMassRate,
                                         const bool   p_IsHeRich) override;                                               // body in WhiteDwarfs.cpp
};



    

    //////////////////////////////////////////////////
    //   MASS TRANSFER                              //
    //////////////////////////////////////////////////

    GNU_CONST MT_CASE DetermineMassTransferTypeAsDonor() const { return MT_CASE::OTHER; }       // Not A, B, C, or NONE // JR FIX THIS: DONE


