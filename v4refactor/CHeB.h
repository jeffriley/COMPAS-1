#pragma once

#include "constants.h"
//// #include "typedefs.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "profiling.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "utils.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#include "FGB.h"


class BaseStar;
class FGB;

class CHeB: virtual public BaseStar, public FGB {

public:

    CHeB() { m_InterimState.SetStellarType(STELLAR_TYPE::CORE_HELIUM_BURNING); };
    
    CHeB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), FGB(p_BaseStar, false) {
        m_InterimState.SetStellarType(STELLAR_TYPE::CORE_HELIUM_BURNING);                                                                                                      // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                         // Initialise if required
    }


    // member functions - alphabetically


double CalculateMinLuminosityOnPhase_Hurley2000(const double p_Mass, const double p_MHeF, const double p_MFGB, const double p_Alpha1, const DBL_VECTOR& p_bN) const; // JR FIX THIS <<<<<<<<<<<<<<<<<<<<




protected:

    void Initialise() {

        m_InterimState.SetStellarType(STELLAR_TYPE::CORE_HELIUM_BURNING);                                                                                                      // Set stellar type
        CalculateTimescales();                                                                                                                                  // Initialise timescales
        m_InterimState.SetAge(m_InterimState.Timescales(TS::HeI));

        EvolveOnPhase(0.0);
    }


    // Member functions (not getters or setters)
    //
    // VIRTUAL FUNCTIONS may be (are expected to be) overridden by derived classes.
    // When overriding virtual functions in a derived class, use the "override" attribute.
    //
    // NON-VIRTUAL FUNCTIONS should not be overridden (declared separately) by derived classes.
    // While it is legal in C++ to declare the same (non-virtual) function in multiple classes,
    // (aka "shadowing", or "hiding"), we discourage it.  Non-virtual functions are statically
    // bound, and as such, especially with indirection, may not produce expected results.



    ////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES   //
    ////////////////////////////////////////

    GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const override; // JR FIX THIS: DONE
    GNU_CONST double CalculateTauAtPhasEnd_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const override { // JR FIX THIS: DONE
        return CalculateTau_Hurley2000(p_Age, p_tScales); // Same as on phase
    }

    double CalculatePhaseLifetime_Hurley2000(
        const double      p_Mass,
        const double      p_CoreMass,
        const double      p_Alpha4,
        const double      p_mHeF,
        const double      p_tBGB,
        const DblVectorT& p_bCoeffs
    ) const;




    ////////////////////////////////////////
    //   MASS                             //
    ////////////////////////////////////////


GNU_PURE inline double CalculateCoreMass_Hurley2000() const override { return CalculateCoreMass_Hurley2000(Mass(), Tau(), GBParams()); } // JR FIX THIS: DONE
GNU_PURE double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau, const GBParamsT& p_GBParams) const; // JR FIX THIS: DONE
GNU_PURE inline double CalculateCoreMassAtPhaseEnd_Hurley2000() const override { return CalculateCoreMass_Hurley2000(); } // Per Hurley SSE code `hrdiag.f` lines 259-265 (tau is calculated, not necessarily 1.0) JR FIX THIS: DONE





GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; } // McCO = 0.0 for CHeB stars

inline double CalculateHeCoreMass() const override { return CalculateCoreMass(); } // McHe = Mc for CHeB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed




    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////


COMPAS_PURE double CalculateRadiusAtPhaseEnd_Hurley2000() const override;
    
    COMPAS_PURE static double CalculateMinRadius_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass); // JR FIX THIS: DONE


double CalculateRadius_Hurley2000(const double p_Mass, const double p_Luminosity, const double p_Tau) const;
double CalculateRadius(const double p_Mass, const double p_Luminosity) const    { return GiantBranch::CalculateRadius(p_Mass, p_Luminosity); }
double CalculateRadius() const { return CalculateRadius(Mass(), Luminosity()); }







COMPAS_PURE double CalculateBluePhaseFBL_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity) const;

COMPAS_PURE double CalculateLifetimeOnBluePhase_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity) const;

void CalculateTimescales_Hurley2000() override { m_InterimState.SetHurleyTimescales(CalculateTimescales_Hurley2000(MassEffectiveInitial(), m_InterimState.CoreMass(), m_MinimumLuminosityOnPhase, m_InterimState.HurleyGBParamsOrDefault(), m_InterimState.HurleyTimescalesOrDefault())); }
COMPAS_PURE DblVectorT CalculateTimescales_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity, const DblVectorT& p_GBParams, const DblVectorT& p_tScales) const;

void CalculateTimescales(const double p_Mass, DblVectorT &p_Timescales) override;








GNU_CONST double CalculateRemnantLuminosity_Hurley2000(const double p_Age, const double p_CoreMass, const DblVectorT& p_tScales) const;

inline double CalculateRemnantRadius_Hurley2000() const override { return CalculateRemnantRadius_Hurley2000(CoreMass(), Tau()); }
GNU_CONST double CalculateRemnantRadius_Hurley2000(const double p_CoreMass, const double p_Tau) const;




COMPAS_PURE static double CalculateMinLuminosity_Hurley2000(const double p_Mass, const double p_MHeF, const double p_MFGB, const DblVectorT& p_bCoeffs);


COMPAS_PURE double CalculateLuminosityAtBluePhaseStart_Hurley2000(const double p_Mass, const double p_CoreMass) const;
COMPAS_PURE double CalculateLuminosityAtBluePhaseEnd_Hurley2000(const double p_Mass, const double p_CoreMass, const DblVectorT& p_tScales) const;
COMPAS_PURE double CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Tau, const double p_CoreMass, const DblVectorT& p_tScales) const;
COMPAS_PURE double CalculateRadius_Hurley2000(const double p_Mass, const double p_Luminosity, const double p_Tau, const double p_CoreMass, const DblVectorT& p_Timescales) const;
COMPAS_PURE double CalculateCELambda_Nanjing_Enhanced(const double p_Mass, const double p_Radius, const double p_CoreMass, const SizeT p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;


///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; } // McCO = 0.0 for CHeB stars

inline double CalculateHeCoreMassAtPhaseEnd() const override { return CalculateCoreMass(); } // McHe = Mc for CHeB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed













GNU_CONST   double CalculateHAbundanceCore(const double p_Tau) const override { return 0.0; } // No hydrogen in the core for CHeB stars
COMPAS_PURE double CalculateHAbundanceSurface(const double p_Tau) const override { return ZDEP->ZAMSHAbundance(Metallicity()); }
    
COMPAS_PURE double CalculateHeAbundanceCore(const double p_Tau) const override;
COMPAS_PURE double CalculateHeAbundanceSurface(const double p_Tau) const override { return ZDEP->ZAMSHeAbundance(Metallicity()); }                      
 




    double          CalculateLifetimeOnBluePhase(const double p_Mass);







                                                                      





    double          CalculateLuminosityAtBluePhaseEnd(const double p_Mass) const {  } // JR FIX THIS: placeholder
    double          CalculateLuminosityAtBluePhaseStart(const double p_Mass) const;

    double          CalculateLuminosityAtPhaseEnd() const                       { return CalculateLuminosityAtBAGB(Mass0()); }
    double          CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau) const;
    double          CalculateLuminosityOnPhase() const                          { return CalculateLuminosityOnPhase(Mass0(), Tau()); }

    double          CalculateRadiusAtBluePhaseEnd(const double p_Mass) const;
    double          CalculateRadiusAtBluePhaseStart(const double p_Mass) const {  } // JR FIX THIS: placeholder
    double          CalculateRadiusAtBluePhaseStart(const double p_Mass, const double p_CoreMass) const;

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const {  } // JR FIX THIS: placeholder
    double          CalculateRadiusAtPhaseEnd() const { return CalculateRadiusAtPhaseEnd(Mass(), Luminosity()); }





COMPAS_PURE static double CalculateMinRadiusOnPhase_Hurley2000_Static(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity);  // JR FIX THIS <<<<<<<<<<<<<<<<<



COMPAS_PURE double CalculateRadiusRho(const double p_Mass, const double p_Tau, const double p_CoreMass, const DblVectorT& p_tScales) const;








GNU_PURE DBL_VECTOR CalculateTimescales_Hurley(const double p_Mass, const DBL_VECTOR& p_MassCutoffs, const DBL_VECTOR& p_Timescales) const;  // JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<< override???



GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const override;

COMPAS_PURE ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const override;


    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                        { return !ShouldEvolveOnPhase(); }                                              // Phase ends when age at or after He Burning
    bool            IsSupernova() const                                         { return false; }                                                               // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash() {  }                                                                                                                   // NO-OP

    bool            ShouldEnvelopeBeExpelledByPulsations() const { return ( OPTIONS->ExpelConvectiveEnvelopeAboveLuminosityThreshold() && DetermineEnvelopeType() == ENVELOPE::CONVECTIVE && utils::Compare( std::log10(Luminosity() / Mass()), OPTIONS->LuminosityToMassThreshold() ) >= 0 ) ; } // JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<
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
    return (1.0 - Metallicity()) * (1.0 - p_Tau);
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
 * double CalculatePhaseLifetime_Hurley2000(
 *     const double      p_Mass,
 *     const double      p_CoreMass,
 *     const double      p_Alpha4,
 *     const double      p_mHeF,
 *     const double      p_tBGB,
 *     const DblVectorT& p_bCoeffs
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_Alpha4                        Hurley alpha4
 * @param       p_MHeF                          Maximum initial mass for which helium ignites degenerately in a Helium Flash (HeF)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (Myr)
 * @param       p_bCoeffs                       Hurley b coefficients
 * @return                                      CHeB lifetime, tHe (Myr)
 */
double CHeB::CalculatePhaseLifetime_Hurley2000(
    const double      p_Mass,
    const double      p_CoreMass,
    const double      p_Alpha4,
    const double      p_mHeF,
    const double      p_tBGB,
    const DblVectorT& p_bCoeffs
) const {

    auto& b = p_bCoeffs;

    double tHe;

    if (p_Mass < p_mHeF) {
        const double tHeMS = HeMS::CalculatePhaseLifetime_Hurley2000_Static(p_CoreMass);
        tHe = (b[39] + ((tHeMS - b[39]) * PPOW((1.0 - (p_Mass / p_mHeF)), b[40]))) * (1.0 + (p_Alpha4 * std::exp(15.0 * (p_Mass - p_mHeF))));
    }
    else {
        const double M5 = utils::intPow(p_Mass, 5);
        tHe = p_tBGB * (((b[41] * PPOW(p_Mass, b[42])) + (b[43] * M5)) / (b[44] + M5));
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
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      CHeB-relative age of the star, [0, 1]
 */
inline double CHeB::CalculateTau_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const {
    return std::clamp((p_Age - p_tScales[TS::HeI]) / p_tScales[TS::He], 0.0, 1.0);
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
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Suggested timestep (Myr)
 */
inline double CHeB::ChooseTimestep_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const {

    const double dtk = 2.0E-3 * p_tScales[TS::HE];                      // Stellar type specific dt (JAR: check 2E-3 vs 2E-2)
    const double dte = p_tScales[TS::HEI] + p_tScales[TS::HE] - p_Age;  // Time to end of phase (change of stellar type)

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);      // Clamp to minimum NUCLEAR_MINIMUM_TIMESTEP
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateMinLuminosity_Hurley2000
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
 * static double CalculateMinLuminosity_Hurley2000_Static(const double p_Mass, const double p_MHeF, const double p_MFGB, const double p_Alpha1, const DblVectorT& p_bCoeffs)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_MHeF                          Maximum initial mass for which helium ignites degenerately in a Helium Flash (HeF)
 * @param       p_MFGB                          Maximum initial mass for which helium ignites on the First Giant Branch (FGB)
 * @param       p_Alpha1                        Hurley alpha1
 * @param       p_bCoeffs                       Hurley b coefficients
 * @return                                      CHeB minimum luminosity (Lsol)
 */
static inline double CHeB::CalculateMinLuminosity_Hurley2000(const double p_Mass, const double p_MHeF, const double p_MFGB, const double p_Alpha1, const DblVectorT& p_bCoeffs) {

    const double c    = (p_bCoeffs[17] / PPOW(p_MFGB, 0.1)) + (((p_bCoeffs[16] * p_bCoeffs[17]) - p_bCoeffs[14]) / (PPOW(p_MFGB, (p_bCoeffs[15] + 0.1))));
    const double lHeI = GiantBranch::CalculateLuminosityAtHeI_Hurley2000(p_Mass, p_MHeF, p_Alpha1, p_bCoeffs);

    return lHeI * ((p_bCoeffs[14] + (c * PPOW(p_Mass, (p_bCoeffs[15] + 0.1)))) / (p_bCoeffs[16] + PPOW(p_Mass, p_bCoeffs[15])));
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
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Remnant core luminosity (Lsol)
 */
double CHeB::CalculateRemnantLuminosity_Hurley2000() const {
    // JR FIX THIS: check original................
    return HeMS::CalculateLuminosityOnPhase_Static(m_CoreMass, CalculateTauOnPhase());
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
 * Calculate the core mass between helium ignition (HeI) and the Base of the
 * Asymptotic Giant Branch (BAGB), per Hurley et al. 2000, eq 67.
 *
 *
 * double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      CHeB core mass of the star (Msol)
 */
inline double CHeB::CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau, const GBParamsT& p_GBParams) const {
    // Should become HeMS star - He mass clamped to total mass
    return std::min(((1.0 - p_Tau) * CalculateCoreMassAtHeI_Hurley2000(p_Metallicity, p_Mass, p_GBParams)) + (p_Tau * CalculateCoreMassAtBAGB_Hurley2000(p_Metallicity, p_Mass)), p_Mass);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateMinRadius_Hurley2000
 *
 * @brief
 * Calculate the minimum radius during Core Helium Burning (on the blue loop),
 * per Hurley et al. 2000, eq 55
 *
 *
 * static double CalculateMinRadius_Hurley2000(const doubl p_Metallicity, const double p_Mass, const double p_CoreMass) 
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core Mass of the star (Msol)
 * @return                                      Minimum radius during Core Helium Burning (Rsol)
 */
double CHeB::CalculateMinRadius_Hurley2000(const doubl p_Metallicity, const double p_Mass, const double p_CoreMass) {

    double radius = 0.0;                                    // Default return value

    const double lMin = ZDEP->HurleyMinLuminosityOnCHeB(p_Metallicity);
    const double mHeF = ZDEP->HurleyMCOHeF(p_Metallicity);
    const auto&  b    = ZDEP->HurleyBCoeffs(p_Metallicity); // Hurley b coefficients

    if (mHeF < p_Mass) {
        const double m_b28 = PPOW(p_Mass, b[28]);           // pow() is slow - do it once only
        radius = ((b[24] * p_Mass) + (PPOW((b[25] * p_Mass), b[26]) * m_b28)) / (b[27] + m_b28);
    }
    else {
        const double lZAHB_mHeF = CalculateLuminosityOnZAHB_Hurley2000(p_Metallicity, mHeF, p_CoreMass);
        const double lZAHB      = CalculateLuminosityOnZAHB_Hurley2000(p_Metallicity p_Mass, p_CoreMass);
        const double mHeF_b28   = PPOW(mHeF, b[28]);        // pow() is slow - do it once only
        const double top        = ((b[24] * mHeF) + (PPOW((b[25] * mHeF), b[26]) * mHeF_b28)) / (b[27] + mHeF_b28);
        const double bottom     = GiantBranch::CalculateRadius_Hurley2000(p_Metallicity, mHeF, lZAHB_mHeF);

        radius = GiantBranch::CalculateRadius_Hurley2000(p_Metallicity, p_Mass, lZAHB) * PPOW(top / bottom, p_Mass / mHeF);
    }

    return radius;
}


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
double CHeB::CalculateRemnantRadius_Hurley2000() const {
    // JR CHECK THIS: check original..........
    return HeMS::CalculateRadiusOnPhase_Static(m_CoreMass, CalculateTauOnPhase());
}





///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



