#ifndef __CHeB_h__
#define __CHeB_h__

#include "constants.h"
//// #include "typedefs.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "profiling.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "utils.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#include "FGB.h"


class BaseStar;
class FGB;

class CHeB: virtual public BaseStar, public FGB {

public:

    CHeB() { m_StellarType = STELLAR_TYPE::CORE_HELIUM_BURNING; };
    
    CHeB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), FGB(p_BaseStar, false) {
        m_StellarType = STELLAR_TYPE::CORE_HELIUM_BURNING;                                                                                                      // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                         // Initialise if required
    }


    // member functions - alphabetically

double CalculateMinLuminosityOnPhase_Hurley2000(const double      p_Mass,
                                                const double      p_MHeF,
                                                const double      p_MFGB,
                                                const double      p_Alpha1,
                                                const DBL_VECTOR& p_bN) const




protected:

    void Initialise() {

        m_StellarType = STELLAR_TYPE::CORE_HELIUM_BURNING;                                                                                                      // Set stellar type
        CalculateTimescales();                                                                                                                                  // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tHeI)];                                                                                                // Set age appropriately
        m_MinimumLuminosityOnPhase = CalculateMinimumLuminosityOnPhase(massCutoffs(MHeF), m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), GLOBALS->HurleyBCoefficients());    // Calculate once, not many
   // SET GLOBAL VALUE <<<<<<<<<<<<<<<<<<<<<<  m_MinLuminosity_CHeB             = CalculateMinLuminosity_Hurley2000_Static(m_HurleyZdependentValues.massCutoffs(static_cast<int>(HURLEY_MCO::MHeF)));

CalculateMinLuminosityOnPhase_Hurley2000(const double      p_Mass,
                                                const double      p_MHeF,
                                                const double      p_MFGB,
                                                const double      p_Alpha1,
                                                const DBL_VECTOR& p_bN);




        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically


///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    
    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////


    COMPAS_PURE inline double CalculateRadiusAtPhaseEnd_Hurley2000() const override { return EAGB::CalculateRadiusOnPhase_Hurley2000_Static(Mass(), Luminosity()); }
    
    COMPAS_PURE static double CalculateMinRadius_Hurley2000_Static(const double p_Mass, const double p_CoreMass);



double CalculateRadius_Hurley2000(const double p_Mass, const double p_Luminosity, const double p_Tau) const;
double CalculateRadius(const double p_Mass, const double p_Luminosity) const    { return GiantBranch::CalculateRadius(p_Mass, p_Luminosity); }
double CalculateRadius() const { return CalculateRadius(m_Mass, m_Luminosity, m_Tau); }







COMPAS_PURE double CalculateBluePhaseFBL_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity) const;

COMPAS_PURE double CalculateLifetimeOnBluePhase_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity) const;

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR CHECK THESE PARAMETERS BELOW
DblVectorT CalculateTimescales_Hurley2000() const override { return CalculateTimescales_Hurley2000(MassEffectiveInitial(), m_InterimState.CoreMass(), GBparams(), m_InterimState.TimeScales()); }
COMPAS_PURE DblVectorT CalculateTimescales_Hurley2000(const double p_Mass, const double p_CoreMass, const DblVectorT& p_GBparams, const DblVectorT& p_tScales) const;


GNU_CONST inline double CalculateCoreMass_Hurley2000() const override { return CalculateCoreMass_Hurley2000(Mass(), Tau()); }
GNU_CONST double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau) const;

GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; } // McCO = 0.0 for CHeB stars

inline double CalculateHeCoreMass() const override { return CalculateCoreMass(); } // McHe = Mc for CHeB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed

inline double CalculateTau_Hurley2000() const override { return CalculateTau_Hurley2000(Age(), m_InterimState.Timescales()); }
GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const;


GNU_CONST double CalculateRemnantLuminosity_Hurley2000(const double p_Age, const double p_CoreMass, const DblVectorT& p_tScales) const;

inline double CalculateRemnantRadius_Hurley2000() const override { return CalculateRemnantRadius_Hurley2000(CoreMass(), Tau()); }
GNU_CONST double CalculateRemnantRadius_Hurley2000(const double p_CoreMass, const double p_Tau) const;




COMPAS_PURE static double CalculateMinLuminosity_Hurley2000_Static(const double p_Mass);


COMPAS_PURE double CalculateLuminosityAtBluePhaseStart_Hurley2000(const double p_Mass, const double p_CoreMass) const;
COMPAS_PURE double CalculateLuminosityAtBluePhaseEnd_Hurley2000(const double p_Mass, const double p_CoreMass, const DblVectorT& p_tScales) const;


///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCoreMassAtPhaseEnd_Hurley2000() const override { return CalculateCoreMass_Hurley2000(); } // per Hurley SSE code `hrdiag.f` lines 259-265 (tau is calculated, not necessarily 1.0)

GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; } // McCO = 0.0 for CHeB stars

inline double CalculateHeCoreMassAtPhaseEnd() const override { return CalculateCoreMass(); } // McHe = Mc for CHeB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed





GNU_PURE  double CalculateCriticalMassRatio_Hurley2002() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_GT_07; }






GNU_CONST   double CalculateHAbundanceCore(const double p_Tau) const override { return 0.0; } // No hydrogen in the core for CHeB stars
COMPAS_PURE double CalculateHAbundanceSurface(const double p_Tau) const override { return GLOBALS->ZAMSHAbundance(); }
    
COMPAS_PURE double CalculateHeAbundanceCore(const double p_Tau) const override;
COMPAS_PURE double CalculateHeAbundanceSurface(const double p_Tau) const override { return GLOBALS->ZAMSHeAbundance(); }                      
 



GNU_CONST double CalculateCELambda_Dewi(const double p_Mass,
                                        const double p_Radius,
                                        const double p_Luminosity,
                                        const double p_RZAMS,
                                        const double p_CoreMass) const;



    double          CalculateLambdaNanjingStarTrack(const double p_Mass) const;
    double          CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;

    double          CalculateLifetimeOnBluePhase(const double p_Mass);






COMPAS_PURE double CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_CoreMass, const DblVectorT& p_tScales) const;
                                                                      





    double          CalculateLuminosityAtBluePhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtBluePhaseStart(const double p_Mass) const;

    double          CalculateLuminosityAtPhaseEnd() const                       { return CalculateLuminosityAtBAGB(m_Mass0); }
    double          CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau) const;
    double          CalculateLuminosityOnPhase() const                          { return CalculateLuminosityOnPhase(m_Mass0, m_Tau); }

    double          CalculateRadiusAtBluePhaseEnd(const double p_Mass) const;
    double          CalculateRadiusAtBluePhaseStart(const double p_Mass) const;

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const;
    double          CalculateRadiusAtPhaseEnd() const                           { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }





COMPAS_PURE static double CalculateMinRadiusOnPhase_Hurley2000_Static(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity);



COMPAS_PURE double CalculateRadiusRho(const double p_Mass, const double p_Tau, const double p_CoreMass, const DblVectorT& p_tScales) const;








GNU PURE  DBL_VECTOR CalculateTimescales_Hurley(const double p_Mass, const DBL_VECTOR& p_MassCutoffs, const DBL_VECTOR& p_Timescales) const; // override??? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const override;

COMPAS_PURE ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const override;


    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                        { return !ShouldEvolveOnPhase(); }                                              // Phase ends when age at or after He Burning
    bool            IsSupernova() const                                         { return false; }                                                               // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash() {  }                                                                                                                   // NO-OP

    bool            ShouldEnvelopeBeExpelledByPulsations() const { return ( OPTIONS->ExpelConvectiveEnvelopeAboveLuminosityThreshold() && DetermineEnvelopeType() == ENVELOPE::CONVECTIVE && utils::Compare( std::log10(m_Luminosity / m_Mass), OPTIONS->LuminosityToMassThreshold() ) >= 0 ) ; }                             // Envelope of convective star with luminosity to mass ratio beyond threshold should be expelled
    bool            ShouldEvolveOnPhase() const;
    bool            ShouldSkipPhase() const                                     { return false; }                                                               // Never skip CHeB phase

};




////// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ABUNDANCE                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateHeAbundanceCore
 *
 * @brief
 * Calculate the helium abundance in the core of the star, given the phase relative age of
 * the star.  Currently just a simple linear model from the initial helium abundance to the
 * maximum helium abundance (assuming that all hydrogen is converted to helium). 
 * 
 * Should one day be updated to match detailed models.
 *
 * double CalculateHeAbundanceCore(const double p_Tau) const
 * 
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      Helium abundance in the core of the star
 */
inline double CHeB::CalculateHeAbundanceCore(const double p_Tau) const {
    return (1.0 - GLOBALS->Metallicity()) * (1.0 - p_Tau);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculatePhaseLifetime_Hurley2000
 *
 * @brief
 * Calculate the lifetime of Core Helium Burning, tHe, per Hurley at al. 2000, eq 57
 *
 *
 * double CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_CoreMass, const DblVectorT& p_tScales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_tScales                       Phase timescales (Myr)
 * @return                                      CHeB lifetime, tHe (Myr)
 */
inline double CHeB::CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_CoreMass, const DblVectorT& p_tScales) const {

    double tHe;

    const DblVectorT b = GLOBALS->HurleyBcoefficients();   
    const double MHeF  = GLOBALS->HurleyMassCutoffs(HURLEY_MCO::HEF);

    if (p_Mass < MHeF) {
        const double tHeMS = HeMS::CalculatePhaseLifetime_Hurley2000_Static(p_CoreMass);
        tHe = (b(39) + ((tHeMS - b(39)) * PPOW((1.0 - (p_Mass / MHeF)), b(40)))) * (1.0 + (GLOBALS->HurleyAlpha4() * std::exp(15.0 * (p_Mass - MHeF))));
    }
    else {
        const double M5 = utils::IntPow(p_Mass, 5);
        tHe = p_tScales[HURLEY_TS::BGB] * (((b(41) * PPOW(p_Mass, b(42))) + (b(43) * M5)) / (b(44) + M5));
    }

    return tHe;
}


/*
 * CalculateTau_Hurley2000
 *
 * @brief
 * Calculate the CHeB-relative age (fractional Core Helium Burning age) of the star,
 * per Hurley et al. 2000, just before eq 59.
 *
 * 
 * double CalculateTau_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tScales                       Phase timescales (Myr)
 * @return                                      CHeB-relative age, [0, 1]
 */
inline double CHeB::CalculateTau_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const {
    return std::max(0.0, std::min(1.0, (p_Age - p_tScales[HURLEY_TS::HEI]) / p_tScales[HURLEY_TS::HE]));
}


/*
 * ChooseTimestep_Hurley2000
 *
 * @brief
 * Choose timestep for evolution
 * See discussion in Hurley et al. 2000, p21
 * The returned value will be clamped to minimum NUCLEAR_MINIMUM_TIMESTEP
 *
 *
 * double ChooseTimestep_Hurley2000(const double p_Age, const DblVectorT& p_tScales)
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tScales                       Phase timescales (Myr)
 * @return                                      Suggested timestep (Myr)
 */
inline double CHeB::ChooseTimestep_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const {

    const double dtk = 2.0E-3 * p_tScales[HURLEY_TS::HE];                      // stellar type specific dt (JAR: check 2E-3 vs 2E-2)
    const double dte = p_tScales[_T_HEI] + p_tScales[HURLEY_TS::HE] - p_Age;   // time to end of phase (change of stellar type)

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);      // clamp to minimum NUCLEAR_MINIMUM_TIMESTEP
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateMinLuminosity_Hurley2000_Static
 *
 * @brief
 * Calculate minimum luminosity during Core Helium Burning for intermediate mass (IM) stars,
 * per Hurley et al. 2000, eq 51
 * 
 * **Ilya**
 * JR: is there a check for IM stars? Should we enforce IM, or remove the caveat/qualifier?
 * Hurley defines LM, IM, and HM at the end of Hurley et al. 2000, sec 5.  Note that the
 * Hurley definition of HM stars differs from the COMPAS constant HIGH_MASS_THRESHOLD
 * introduced (I think) by JD with the winds update.
 *
 *
 * static double CalculateMinLuminosity_Hurley2000_Static(const double p_Mass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      CHeB minimum luminosity (Lsol)
 */
inline double CHeB::CalculateMinLuminosity_Hurley2000_Static(const double p_Mass) {

    const DblVectorT b = GLOBALS->HurleyBcoefficients();   
    const double MFGB  = GLOBALS->HurleyMassCutoffs(_M_FGB_);
            
    const double c    = (b[17] / PPOW(MFGB, 0.1)) + (((b[16] * b[17]) - b[14]) / (PPOW(MFGB, (b[15] + 0.1))));
    const double lHeI = GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass);

    return lHeI * ((b[14] + (c * PPOW(p_Mass, (b[15] + 0.1)))) / (b[16] + PPOW(p_Mass, b[15])));
}


/*
 * CalculateRemnantLuminosity_Hurley2000
 *
 * @brief
 * Calculate luminosity of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass), per Hurley et al. 2000, just after eq 105
 *
 * This function calls CalculateTau_Hurley2000() - we can't assume we can use the (perhaps)
 * already calculated value of Tau (from the interim state) becaue we can't control what
 * the caller passes as values for p_Age and p_tScales.
 * 
 *
 * double CalculateRemnantLuminosity_Hurley2000(const double p_Age, const double p_CoreMass, const DblVectorT& p_tScales)
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_tScales                       Phase timescales (Myr)
 * @return                                      Remnant core luminosity (Lsol)
 */
inline double CHeB::CalculateRemnantLuminosity_Hurley2000(const double p_Age, const double p_CoreMass, const DblVectorT& p_tScales) const {
    return HeMS::CalculateLuminosityOnPhase_Static(p_CoreMass, CalculateTau_Hurley2000(p_Age, p_tScales));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCoreMass_Hurley2000
 *
 * @brief
 * Calculate the core mass between Helium Ignition (HeI) and the Base of the
 * Asymptotic Giant Branch (BAGB), per Hurley et al. 2000, eq 67.
 *
 *
 * double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      CHeB core mass (Msol)
 */
inline double CHeB::CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau) const {
    // should become HeMS star - He mass clamped to total mass
    return std::min(((1.0 - p_Tau) * CalculateCoreMassAtHeI_Hurley2000(p_Mass)) + (p_Tau * CalculateCoreMassAtBAGB(p_Mass)), p_Mass);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRemnantRadius_Hurley2000
 *
 * @brief
 * Calculate radius of the remnant the star would become if it lost all of its envelope
 * immediately (i.e. M = Mc), per Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantRadius_Hurley2000(const double p_CoreMass, const double p_Tau)
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      Remnant core radius (Rsol)
 */
inline double CHeB::CalculateRemnantRadius_Hurley2000(const double p_CoreMass, const double p_Tau) const {
    return HeMS::CalculateRadius_Hurley2000_Static(p_CoreMass, p_Tau);
}


#endif // __CHeB_h__
