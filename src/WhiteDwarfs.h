#ifndef __WhiteDwarfs_h__
#define __WhiteDwarfs_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Remnants.h"


class BaseStar;
class Remnants;

class WhiteDwarfs: virtual public BaseStar, public Remnants {

public:

    WhiteDwarfs(){};

    WhiteDwarfs(const BaseStar &p_BaseStar) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {}


    // member functions


    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////


    COMPAS_PURE static double CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Age);




    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////





inline double CalculateRadius_Hurley2000() const override { return CalculateRadius_Hurley2000_Static(Mass()); }
GNU_CONST static inline double CalculateRadius_Hurley2000_Static(const double p_Mass) { return CalculateRadius_Marsh2004_Static(p_Mass); }
GNU_CONST static double CalculateRadius_Marsh2004_Static(const double p_Mass);


    MT_CASE          DetermineMassTransferTypeAsDonor() const                                                   { return MT_CASE::OTHER; }                                      // Not A, B, C, or NONE

    ACCRETION_REGIME DetermineAccretionRegime(const double p_DonorThermalMassLossRate, const bool p_HeRich);                                                                    // Get the current accretion regime. Can also change m_HeShellDetonation and m_OffCentreIgnition flags.
    
    void             ResolveShellChange(const double p_AccretedMass);


protected:
    // member variables

            bool             m_HeShellDetonation;                                                                                                                               // Flag to initialize He-Shell detonation (i.e. as described in Wang. 2018, sect 5 2018RAA....18...49W)
            double           m_HeShell;                                                                                                                                         // Current WD He-shell size (Msol). Increases through accretion.
            double           m_HShell;                                                                                                                                          // Current WD H-shell size (Msol). Increases through accretion.
            double           m_L0Ritter;                                                                                                                                        // Parameter from numerical calculations, see Ritter 1999, section 3. Eqs 10 and 12, as well as table 2. Corresponds to L0.
            double           m_LambdaRitter;                                                                                                                                    // Parameter from numerical calculations, see Ritter 1999, section 3. Eqs 10 and 12, as well as table 2.
            bool             m_OffCentreIgnition;                                                                                                                               // Flag for CO WD evolution into ONe WD
            bool             m_ShouldRejuvenate;                                                                                                                                // Flag for evolution of HeWD back into HeMS
            bool             m_IsSubChandrasekharTypeIa;                                                                                                                        // Flag for SubCh SN of HeWD
            double           m_XRitter;                                                                                                                                         // Assumed hydrogen-mass fraction of material being accreted by He WD, as in Ritter 1999, table 2.
            ACCRETION_REGIME m_AccretionRegime;
            
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






///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<        


inline double CalculateCOCoreMass() const override { return COCoreMass(); } // McCO constant for WDs

inline double CalculateHeCoreMass() const override { return HeCoreMass(); } // McHe constant for WDs


            

            double           CalculateEtaH(const double p_MassIntakeRate);

            double           CalculateEtaHe(const double p_MassIntakeRate);

            double           CalculateEtaPTY(const double p_MassIntakeRate);

            double           Calculatel0Ritter() const                                                          { return (utils::Compare(m_Metallicity, 0.01) > 0) ? L0_RITTER_HIGH_Z : L0_RITTER_LOW_Z; }

            double           CalculateXRitter() const                                                           { return utils::Compare(m_Metallicity, 0.01) > 0 ? 0.7 : 0.8; } // Assumed Hydrogen-mass fraction

            double           CalculateLambdaRitter() const                                                      { return utils::Compare(m_Metallicity, 0.01) > 0 ? 8.0 : 5.0; } // Exponent for the assumed core-mass and luminosity relationship in Ritter 1999

            double           CalculateInitialSupernovaMass() const                                              { return 0.0; }








GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::CONVECTIVE; } // Always CONVECTIVE for White Dwarfs


            bool             IsMassAboveChandrasekhar() const                                                   { return (utils::Compare(m_Mass, MCH) > 0); }                   // Mass exceeds Chandrasekhar limit 







///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


inline double CalculateCOCoreMassAtPhaseEnd() const override { return COCoreMass(); } // McCO constant for WDs

inline double CalculateCOCoreMassAtPhaseEnd() const override { return HeCoreMass(); } // McHe constant for WDs




            STELLAR_TYPE     ResolveAIC();  
            STELLAR_TYPE     ResolveSNIa();  
            STELLAR_TYPE     ResolveHeSD();  
            //STELLAR_TYPE     ResolveSupernova()                                                                 { return EvolveToNextPhase(); }                                 // SNe for WDs are handled internally to each WD type


};



///////////// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
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
 * Calculate the luminosity of a White Dwarf as it cools, per Hurley et al. 2000, eq 90
 *
 *
 * double CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Age)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Time                          Age since White Dwarf formation (Myr)
 * @return                                      Luminosity of the White Dwarf (Lsol)
 */
inline double WhiteDwarfs::CalculateLuminosity_Hurley2000_Static(const double p_Mass, const double p_Age) {
    return (635.0 * p_Mass * PPOW(GLOBALS->Metallicity(), 0.4)) / PPOW(WD_Baryon_Number.at(STELLAR_TYPE::HELIUM_WHITE_DWARF) * (p_Age + 0.1), 1.4);
}




/// WhiteDwarfs_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


class WhiteDwarfs_Constituent: virtual public BinaryConstituentStar, public WhiteDwarfs {

public:


protected:

};


            ACCRETION_REGIME WhiteDwarfAccretionRegime() const                                                  { return m_AccretionRegime; }
    

            Dbl_DblT CalculateMassAcceptanceRate(const double p_DonorMassRate, const bool p_IsHeRich);          
            Dbl_DblT CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate, const bool p_IsHeRich) { return CalculateMassAcceptanceRate(p_DonorMassRate, p_IsHeRich); }




#endif // __WhiteDwarfs_h__
