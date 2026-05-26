#ifndef __EAGB_h__
#define __EAGB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "CHeB.h"


class BaseStar;
class CHeB;

class EAGB: virtual public BaseStar, public CHeB {

public:

    EAGB() { m_StellarType = STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH; };
    
    EAGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), CHeB(p_BaseStar, false) {
        m_StellarType = STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH;                                                                                                    // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                 // Initialise if required
    }





protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                                          // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tHeI)] + m_Timescales[static_cast<int>(TIMESCALE::tHe)];                                                       // Set age appropriately
        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically






///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


DBL_VECTOR CalculateTimescales_Hurley2000() const override { return CalculateTimescales_Hurley2000(MassEffectiveInitial(), m_InterimState.GBparams(), m_InterimState.TimeScales()); }
COMPAS_PURE DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;



inline double CalculateCOCoreMass_Hurley2000() const override { return CalculateCOCoreMass_Hurley2000(Age(), m_InterimState.GBparams(), m_InterimState.Timescales()); }
GNU_CONST double CalculateCOCoreMass_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;

inline double CalculateHeCoreMass() const override { return HeCoreMass(); } // McHe is constant for EAGB stars



inline double CalculateLuminosity_Hurley2000() const override {
    // use Hurley core mass - luminosity relationship, per Hurley et al. eq 37
    return CalculateLuminosityGivenCoreMass_Hurley2000(CoreMass(), GBparams());
}




    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////

    COMPAS_PURE static double CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Luminosity);
GNU_CONST   inline double EAGB::CalculateRemnantRadius_Hurley2000(const double p_HeCoreMass) const;


GNU_CONST inline double CalculateTau_Hurley2000() const override { return 0.0; }; // Tau (relative age) is not used for EAGB stars in Hurley et al. 2000, so we return 0.0




///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd_Hurley2000() const override { return CurrentState.GBparams(McDU); }    // McCO = McDU at phase end for EAGB stars, per Hurley et al. 2000, section 5.4

inline double CalculateHeCoreMassAtPhaseEnd() const override { return HeCoreMass(); } // McHe is constant for EAGB stars


    double          CalculateCoreMassAtPhaseEnd() const                                             { GBparams[static_cast<int>(HURLEY_GBP:::McDU)]; }                 // Mc(EAGB) = McDU at phase end (Hurley et al. 2000, section 5.4)
    double          CalculateCoreMassOnPhase() const                                                { GBparams[static_cast<int>(HURLEY_GBP:::McBAGB)]; }               // Mc(EAGB) = McHe(EAGB) = McBAGB on phase (Hurley et al. 2000, section 5.4)

    double          CalculateInitialSupernovaMass() const                                           { GBparams[static_cast<int>(HURLEY_GBP:::McBAGB)]; }               // For EAGB & TPAGB we use the mass at Base Asymptotic Giant Branch to determine SN type

    double          CalculateLambdaNanjingStarTrack(const double p_Mass) const;
    double          CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;

    double          CalculateLifetimeTo2ndDredgeUp(const double p_Tinf1_FAGB, const double p_Tinf2_FAGB) const;

    double          CalculateLuminosityAtPhaseEnd(const double p_CoreMass) const                    { return CalculateLuminosityOnPhase(p_CoreMass); }                  // Same as on phase
    double          CalculateLuminosityAtPhaseEnd() const                                           { return CalculateLuminosityAtPhaseEnd(m_COCoreMass); }             // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_CoreMass) const;
    double          CalculateLuminosityOnPhase() const                                              { return CalculateLuminosityOnPhase(m_COCoreMass); }



    COMPAS_PURE MASS_LOSS_T CalculateMLrate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_PerturbationMu = 0.0) const override;



    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const { return CalculateRadiusOnPhase(p_Mass, p_Luminosity); }            // Same as on phase
    double          CalculateRadiusAtPhaseEnd() const                                               { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }         // Use class member variables
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, m_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)], m_BnCoefficients); }
    double          CalculateRadiusOnPhase() const                                                  { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }            // Use class member variables

    double          CalculateRemnantLuminosity() const;
    double          CalculateRemnantRadius() const;







GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const override;


GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::CONVECTIVE; } // Always CONVECTIVE for EAGB stars


    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                                            { return !ShouldEvolveOnPhase(); }                                  // Phase ends when age at or after DU timescale, and no TPAGB
    bool            IsSupernova() const;

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash() {  }                                                                                                                           // NO-OP
    STELLAR_TYPE    ResolveSkippedPhase()                                                           { return m_StellarType; }                                           // NO-OP

    bool            ShouldEvolveOnPhase() const;
    bool            ShouldSkipPhase() const;

};




///////////// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLifetimeTo2ndDU_Hurley2000
 *
 * @brief
 * Calculate the lifetime to second dredge up (n.b. there is no explicit first),
 * per Hurley et al. 2000, eqs 70, 71 & 72
 * 
 * This is the time of transition between the EAGB and the TPAGB
 *
 * 
 * double CalculateLifetimeTo2ndDU_Hurley2000(const DBL_VECTOR& p_GBparams, const double p_Tinf1_FAGB, const double p_Tinf2_FAGB)
 *
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tInf1_FAGB                    EAGB integration constant 1, tinf1_FAGB (per Hurley timescales) (Myr)
 * @param       p_tInf2_FAGB                    EAGB integration constant 2, tinf2_FAGB (per Hurley timescales) (Myr)
 * @return                                      Lifetime to second dredge up (tDU) (Myr)
 */
double EAGB::CalculateLifetimeTo2ndDU_Hurley2000(const DBL_VECTOR& p_GBparams, const double p_Tinf1_FAGB, const double p_Tinf2_FAGB) const {
#define GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)]

    const double lDU = CalculateLuminosityGivenCoreMass_Hurley2000(GBparams(McDU), p_GBparams);
    const double p1  = GBparams(p) - 1.0;
    const double q1  = GBparams(q) - 1.0;

    return lDU <= GBparams(Lx)
            ? p_Tinf1_FAGB - (1.0 / (p1 * GBparams(AHe) * GBparams(D))) * PPOW((GBparams(D) / lDU), (p1 / GBparams(p)))
            : p_Tinf2_FAGB - (1.0 / (q1 * GBparams(AHe) * GBparams(B))) * PPOW((GBparams(B) / lDU), (q1 / GBparams(q)));

#undef GBparams
}


/*
 * ChooseTimestep_Hurley2000
 *
 * @brief
 * Choose timestep for evolution
 * See the discussion in Hurley et al. 2000, p21
 * The returned value will be clamped to minimum NUCLEAR_MINIMUM_TIMESTEP
 *
 *
 * double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales)
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Phase timescales (Myr)
 * @return                                      Suggested timestep (Myr)
 */
GNU_CONST inline double EAGB::ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const {
#define tScales(x) p_tScales[static_cast<int>(TIMESCALE::x)]

    const double dtk = 0.02 * ((p_Age <= tScales(tMx_FAGB) ? tScales(tinf1_FAGB) : tScales(tinf2_FAGB)) - p_Age);   // stellar type specific dt

    // time to end of phase (change of stellar type, dte) not used here - how to calculate?
    // clamp to minimum NUCLEAR_MINIMUM_TIMESTEP
    return std::max(dtk, NUCLEAR_MINIMUM_TIMESTEP);

#undef tScales
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRemnantLuminosity_Hurley2000
 *
 * @brief
 * Calculate luminosity of the remnant the star would become if it lost all of its envelope
 * immediately (i.e. M = Mc), per Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantLuminosity_Hurley2000(const double p_COCoreMass, const DBL_VECTOR& p_GBparams)
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Remnant core luminosity (Lsol)
 */
GNU_CONST inline double EAGB::CalculateRemnantLuminosity_Hurley2000(const double p_COCoreMass, const double p_GBparamB, const double p_GBparamD) const {
    return HeGB::CalculateLuminosity_Hurley2000_Static(m_COCoreMass, p_GBparamB, const double p_GBparamD);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCOCoreMass_Hurley2000
 *
 * @brief
 * Calculate CO core mass on the Early Asymptotic Giant Branch,
 * per Hurley et al. 2000, eq 39, modified as described in Section 5.4
 *
 *
 * double CalculateCOCoreMass_Hurley2000(const double p_Time, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales)
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      EAGB CO core mass (Msol)
 */
GNU_CONST inline double EAGB::CalculateCOCoreMass_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const {
// macros for convenience and readability - undefined at end of function
#define GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)]
#define tScales(x) p_tScales[static_cast<int>(TIMESCALE::x)]

    return p_Age <= tScales(tMx_FAGB)
            ? PPOW((GBparams(p) - 1.0) * GBparams(AHe) * GBparams(D) * (tScales(tinf1_FAGB) - p_Age), 1.0 / (1.0 - GBparams(p)))
            : PPOW((GBparams(q) - 1.0) * GBparams(AHe) * GBparams(B) * (tScales(tinf2_FAGB) - p_Age), 1.0 / (1.0 - GBparams(q)));

#undef tScales
#undef GBparams
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
 * double CalculateRemnantRadius_Hurley2000(const double p_HeCoreMass)
 *
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @return                                      Remnant core radius (Rsol)
 */
GNU_CONST inline double EAGB::CalculateRemnantRadius_Hurley2000(const double p_HeCoreMass) const {
    double R1, R2;
    std::tie(R1, R2) = HeGB::CalculateRadius_Hurley2000_Static(p_HeCoreMass, CalculateRemnantLuminosity_Hurley2000());

    return std::min(m_Radius, std::min(R1, R2));
}





///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    MT_CASE         DetermineMassTransferTypeAsDonor() const                                        { return MT_CASE::C; }                                              // Always case C

double CalculateCriticalMassRatio_Hurley2002() const { return GiantBranch_Constituent::CalculateCriticalMassRatio_Hurley2002(); }



#endif // __EAGB_h__
