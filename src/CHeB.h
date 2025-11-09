#ifndef __CHeB_h__
#define __CHeB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

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

    CHeB* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        CHeB* clone = new CHeB(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static CHeB* Clone(CHeB& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        CHeB* clone = new CHeB(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


    // member functions - alphabetically

double CalculateMinLuminosityOnPhase_Hurley2000(const double      p_Mass,
                                                const double      p_MHeF,
                                                const double      p_MFGB,
                                                const double      p_Alpha1,
                                                const DBL_VECTOR& p_bN) const


    static double CalculateMinimumRadiusOnPhase_Static(const double      p_Mass,
                                                       const double      p_CoreMass,
                                                       const double      p_Alpha1,
                                                       const double      p_MHeF,
                                                       const double      p_MFGB,
                                                       const double      p_LuminosityOnPhase,
                                                       const DBL_VECTOR &p_BnCoefficients);


protected:

    void Initialise() {

        m_StellarType = STELLAR_TYPE::CORE_HELIUM_BURNING;                                                                                                      // Set stellar type
        CalculateTimescales();                                                                                                                                  // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tHeI)];                                                                                                // Set age appropriately
        m_MinimumLuminosityOnPhase = CalculateMinimumLuminosityOnPhase(massCutoffs(MHeF), m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), GLOBALS->HurleyBCoefficients());    // Calculate once, not many
   // SET GLOBAL VALUE <<<<<<<<<<<<<<<<<<<<<<  m_MinLuminosity_CHeB             = CalculateMinLuminosity_Hurley2000_Static(m_HurleyZdependentValues.massCutoffs(static_cast<int>(HURLEY_MASS_CUTOFFS::MHeF)));

CalculateMinLuminosityOnPhase_Hurley2000(const double      p_Mass,
                                                const double      p_MHeF,
                                                const double      p_MFGB,
                                                const double      p_Alpha1,
                                                const DBL_VECTOR& p_bN);




        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically

    double          CalculateBluePhaseFBL(const double p_Mass);


///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

DBL_VECTOR CalculateTimescales_Hurley2000() const override {
    return CalculateTimescales_Hurley2000(
        m_StateHistory.CurrentState.MassEffectiveInitial(),
        m_StateHistory.CurrentState.CoreMass(),
        m_StateHistory.CurrentState.GBparams(),
        m_StateHistory.CurrentState.TimeScales()
    );
}
COMPAS_PURE DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const double p_CoreMass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;


GNU_CONST inline double CalculateCoreMass_Hurley2000() const override {
    return CalculateCoreMass_Hurley2000(m_StateHistory.CurrentState.Mass(), m_StateHistory.CurrentState.Tau());
}
GNU_CONST double CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau) const;

GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; } // McCO = 0.0 for CHeB stars

inline double CalculateHeCoreMass() const override { return CalculateCoreMass(); } // McHe = Mc for CHeB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed

GNU_CONST inline double CalculateTau_Hurley2000() const override {
    return CalculateTau_Hurley2000(
        m_StateHistory.CurrentState.Age(),
        m_StateHistory.CurrentState.Timescales(tHeI),
        m_StateHistory.CurrentState.Timescales(tHe)
    );
}
GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const double p_tHeI, const double p_tHe) const;


GNU_CONST double CalculateRemnantLuminosity_Hurley2000(const double p_Age, const double p_CoreMass, const double p_tHeI, const double p_tHe) const;

inline double CalculateRemnantRadius_Hurley2000() const override {
    return CalculateRemnantRadius_Hurley2000(m_StateHistory.CurrentState.CoreMass(), m_StateHistory.CurrentState.Tau());
}
GNU_CONST double CalculateRemnantRadius_Hurley2000(const double p_CoreMass, const double p_Tau) const;




GNU_CONST static double CalculateMinLuminosity_Hurley2000_Static(
    const double      p_Mass,
    const double      p_MHeF,
    const double      p_MFGB,
    const double      p_Alpha1,
    const DBL_VECTOR& p_bCoeffs
);

COMPAS_PURE double CHeB::CalculateLuminosityAtBluePhaseStart_Hurley2000(const double p_Mass, const double p_CoreMass) const;


///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCoreMassAtPhaseEnd_Hurley2000() const override {
    return CalculateCoreMass_Hurley2000(); // per Hurley SSE code `hrdiag.f` lines 259-265 (tau is calculated, not necessarily 1.0)
}

GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; } // McCO = 0.0 for CHeB stars

inline double CalculateHeCoreMassAtPhaseEnd() const override { return CalculateCoreMass(); } // McHe = Mc for CHeB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed


GNU CONST inline double CalculateRadiusAtPhaseEnd_Hurley2000() const override {
    return EAGB::CalculateRadiusOnPhase_Hurley2000_Static(m_StateHistory.CurrentState.Mass(), m_StateHistory.CurrentState.Luminosity());
}


double CalculateRadiusOnPhase_Hurley(const double p_Mass, const double p_Luminosity, const double p_Tau) const;
double CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return GiantBranch::CalculateRadiusOnPhase(p_Mass, p_Luminosity); }
double CalculateRadiusOnPhase() const { return CalculateRadiusOnPhase(m_Mass, m_Luminosity, m_Tau); }




GNU_PURE  double CalculateCriticalMassRatio_Hurley2002() const { return HURLEY_HJELLMING_WEBBINK_QCRIT_MS_GT_07; }






GNU_CONST double CalculateHAbundanceCore(const double p_Tau, const double p_InitialHAbundance) const override { return 0.0; } // No hydrogen in the core for CHeB stars
GNU_PURE  double CalculateHAbundanceSurface(const double p_Tau, const double p_InitialHAbundance) const override { return GLOBALS->ZAMSHAbundance(); }
    
GNU_CONST double CalculateHeAbundanceCore(const double p_Tau, const double p_InitialHeAbundance = 0.0) const override;
GNU_PURE  double CalculateHeAbundanceSurface(const double p_Tau, const double p_InitialHeAbundance) const override { return GLOBALS->ZAMSHeAbundance(); }                      
 



GNU_CONST double CalculateCELambda_Dewi(const double p_Mass,
                                        const double p_Radius,
                                        const double p_Luminosity,
                                        const double p_RZAMS,
                                        const double p_CoreMass) const;



    double          CalculateLambdaNanjingStarTrack(const double p_Mass) const;
    double          CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;

    double          CalculateLifetimeOnBluePhase(const double p_Mass);






COMPAS_PURE double CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_tBGB) const;
                                                                      





    double          CalculateLuminosityAtBluePhaseEnd(const double p_Mass) const;
    double          CalculateLuminosityAtBluePhaseStart(const double p_Mass) const;

    double          CalculateLuminosityAtPhaseEnd() const                       { return CalculateLuminosityAtBAGB(m_Mass0); }
    double          CalculateLuminosityOnPhase(const double p_Mass, const double p_Tau) const;
    double          CalculateLuminosityOnPhase() const                          { return CalculateLuminosityOnPhase(m_Mass0, m_Tau); }

    double          CalculateRadiusAtBluePhaseEnd(const double p_Mass) const;
    double          CalculateRadiusAtBluePhaseStart(const double p_Mass) const;

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const;
    double          CalculateRadiusAtPhaseEnd() const                           { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }





static double CalculateMinRadiusOnPhase_Hurley2000_Static(const double      p_Mass,
                                                          const double      p_CoreMass,
                                                          const double      p_Alpha1,
                                                          const double      p_MHeF,
                                                          const double      p_MFGB,
                                                          const double      p_MinLuminosity,
                                                          const DBL_VECTOR& p_bN);



    double          CalculateRadiusRho(const double p_Mass, const double p_Tau) const;








GNU PURE  DBL_VECTOR CalculateTimescales_Hurley(const double p_Mass, const DBL_VECTOR& p_MassCutoffs, const DBL_VECTOR& p_Timescales) const;


    double          ChooseTimestep(const double p_Time) const;

    ENVELOPE        DetermineEnvelopeType() const;

    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                        { return !ShouldEvolveOnPhase(); }                                              // Phase ends when age at or after He Burning
    bool            IsSupernova() const                                         { return false; }                                                               // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash() {  }                                                                                                                   // NO-OP

    bool            ShouldEnvelopeBeExpelledByPulsations() const { return ( OPTIONS->ExpelConvectiveEnvelopeAboveLuminosityThreshold() && DetermineEnvelopeType() == ENVELOPE::CONVECTIVE && utils::Compare( log10(m_Luminosity/m_Mass), OPTIONS->LuminosityToMassThreshold() ) >= 0 ) ; }                             // Envelope of convective star with luminosity to mass ratio beyond threshold should be expelled
    bool            ShouldEvolveOnPhase() const;
    bool            ShouldSkipPhase() const                                     { return false; }                                                               // Never skip CHeB phase

};




////// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         AGE / LIFETIME / TAU / TIMESCALES                         //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculatePhaseLifetime_Hurley2000
 *
 * @brief
 * Calculate the lifetime of Core Helium Burning, tHe, per Hurley at al. 2000, eq 57
 *
 *
 * double CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_tBGB) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (per Hurley timescales) (Myr)
 * @return                                      CHeB lifetime, tHe (Myr)
 */
COMPAS_PURE inline double CHeB::CalculatePhaseLifetime_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_tBGB) const {
#define b(x) GLOBALS->HurleyBCoefficients(x) // for convenience and readability - undefined at end of function

    double tHe;

    const double mHeF = GLOBALS->HurleyMassCutoffs(static_cast<int>(HURLEY_MASS_CUTOFF:MHeF));

    if (p_Mass < mHeF) {
        const double tHeMS = HeMS::CalculatePhaseLifetime_Hurley2000_Static(p_CoreMass);
        tHe = (b(39) + ((tHeMS - b(39)) * PPOW((1.0 - (p_Mass / mHeF)), b(40)))) * (1.0 + (GLOBALS->HurleyAlpha4() * exp(15.0 * (p_Mass - mHeF))));
    }
    else {
        const double m5 = utils::IntPow(p_Mass, 5);
        tHe = p_tBGB * (((b(41) * PPOW(p_Mass, bN(42))) + (b(43) * m5)) / (b(44) + m5));
    }

    return tHe;

#undef bN
}


/*
 * CalculateTau_Hurley2000
 *
 * @brief
 * Calculate the CHeB-relative age (fractional Core Helium Burning age) of the star,
 * per Hurley et al. 2000, just before eq 59.
 *
 * 
 * double CalculateTau_Hurley2000(const double p_Age, const double p_tHeI, const double p_tHe) const
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tHeI                          Time to helium ignition (per Hurley timescales) (Myr)
 * @param       p_tHe                           Time to helium burning (per Hurley timescales) (Myr)
 * @return                                      CHeB-relative age, [0, 1]
 */
GNU_CONST inline double CHeB::CalculateTau_Hurley2000(const double p_Age, const double p_tHeI, const double p_tHe) const {
    return std::max(0.0, std::min(1.0, (p_Age - p_tHeI) / p_tHe));
}


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
 * double CalculateHeAbundanceCore(const double p_Tau, const double p_InitialHeAbundance) const
 * 
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_InitialHeAbundance            Initial helium abundance of the star (not used here)
 * @return                                      Helium abundance in the core of the star
 */
COMPAS_PURE inline double CHeB::CalculateHeAbundanceCore(const double p_Tau, const double p_InitialHeAbundance) const {
    return (1.0 - GLOBALS->ReferenceMetallicity()) * (1.0 - p_Tau);
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
GNU_CONST inline double CHeB::CalculateRemnantRadius_Hurley2000(const double p_CoreMass, const double p_Tau) const {
    return HeMS::CalculateRadius_Hurley2000_Static(p_CoreMass, p_Tau);
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
 * @return                                      EAGB core mass (Msol)
 */
GNU_CONST inline double CHeB::CalculateCoreMass_Hurley2000(const double p_Mass, const double p_Tau) const {
    // should become HeMS star - He mass clamped to total mass
    return std::min(((1.0 - p_Tau) * CalculateCoreMassAtHeI_Hurley2000(p_Mass)) + (p_Tau * CalculateCoreMassAtBAGB(p_Mass)), p_Mass);
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
 * static double CalculateMinLuminosity_Hurley2000_Static(
 *     const double p_Mass, 
 *     const double p_MHeF, 
 *     const double p_MFGB, 
 *     const double p_Alpha1, 
 *     const DBL_VECTOR& p_bCoeffs
 * )
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_MFGB                          Maximum initial mass at helium ignition on the FGB (Hurley masscutoffs[MFGB]) (Msol)
 * @param       p_Alpha1                        Hurley alpha1 constant
 * @param       p_bCoeffs                       Hurley b coefficients
 * @return                                      CHeB minimum luminosity (Lsol)
 */
GNU_CONST inline double CHeB::CalculateMinLuminosity_Hurley2000_Static(
    const double      p_Mass,
    const double      p_MHeF,
    const double      p_MFGB,
    const double      p_Alpha1,
    const DBL_VECTOR& p_bCoeffs
) {
    const double c    = (p_bCoeffs[17] / PPOW(p_MFGB, 0.1)) + (((p_bCoeffs[16] * p_bCoeffs[17]) - p_bCoeffs[14]) / (PPOW(p_MFGB, (p_bCoeffs[15] + 0.1))));
    const double lHeI = GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass, p_MHeF, p_Alpha1, p_bCoeffs);

    return lHeI * ((p_bCoeffs[14] + (c * PPOW(p_Mass, (p_bCoeffs[15] + 0.1)))) / (p_bCoeffs[16] + PPOW(p_Mass, p_bCoeffs[15])));
}


/*
 * CalculateRemnantLuminosity_Hurley2000
 *
 * @brief
 * Calculate luminosity of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass), per Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantLuminosity_Hurley2000(const double p_Age, const double p_CoreMass, const double p_tHeI, const double p_tHe)
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_tHeI                          Time to helium ignition (per Hurley timescales) (Myr)
 * @param       p_tHe                           Time to helium burning (per Hurley timescales) (Myr)
 * @return                                      Remant core luminosity (Lsol)
 */
GNU_CONST inline double CHeB::CalculateRemnantLuminosity_Hurley2000(const double p_Age, const double p_CoreMass, const double p_tHeI, const double p_tHe) const {
    return HeMS::CalculateLuminosityOnPhase_Static(p_CoreMass, CalculateTau_Hurley2000(p_Age, p_tHeI, p_tHe));
}


#endif // __CHeB_h__
