#pragma once

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "CHeB.h"


class BaseStar;
class CHeB;

class EAGB: virtual public BaseStar, public CHeB {



public:

    EAGB() { m_InterimState.SetStellarType(STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH); };
    
    EAGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), CHeB(p_BaseStar, false) {
        m_InterimState.SetStellarType(STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH);                                                                                                    // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                 // Initialise if required
    }






protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                                          // Initialise timescales
        m_InterimState.SetAge(m_InterimState.HurleyTimescales(HURLEY_TIMESCALES::HeI) + m_InterimState.HurleyTimescales(HURLEY_TIMESCALES::He));                       // Set age appropriately
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






inline double CalculateLuminosity_Hurley2000() const override {
    // use Hurley core mass - luminosity relationship, per Hurley et al. eq 37
    return CalculateLuminosity_Hurley2000(CoreMass(), m_InterimState.HurleyGBParams());
}
    using BaseStar::CalculateLuminosity_Hurley2000;



    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////

    // Tau (relative age) is not used for EAGB stars in Hurley et al. 2000, so we return 0.0
    GNU_CONST inline double CalculateTau_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const override { return 0.0; }; // JR FIX THIS: DONE
    GNU_CONST inline double CalculateTauAtPhaseEnd_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const override { return 0.0; }; // JR FIX THIS: DONE


void EAGB::CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) {
    p_Timescales = CalculateTimescales_Hurley2000(p_Mass, m_InterimState.HurleyGBParams(), p_Timescales);
}
void CalculateTimescales_Hurley2000() override { m_InterimState.SetHurleyTimescales(CalculateTimescales_Hurley2000(MassEffectiveInitial(), m_InterimState.HurleyGBParamsOrDefault(), m_InterimState.HurleyTimescalesOrDefault())); }
COMPAS_PURE DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const;

void CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) override;




    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    GNU_CONST double CalculateRemnantLuminosity_Hurley2000(const double p_COCoreMass, const GBParamsT& p_GBParams) const; // JR FIX THIS: DONE
    GNU_CONST double CalculateRemnantLuminosity_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateRemnantLuminosity_Hurley2000(COCoreMass(), GBParams());
    }




    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////

    GNU_PURE inline double CalculateCoreMass_Hurley2000() const override{ return GBParams(HURLEY_GBP::McBAGB); } // Mc(EAGB) = McHe(EAGB) = McBAGB (Hurley et al. 2000, section 5.4) JR FIX THIS: DONE
    GNU_PURE inline double CalculateCoreMassAtPhaseEnd_Hurley2000() const override{ return GBParams(HURLEY_GBP::McDU); } // Mc(EAGB) = McDU (Hurley et al. 2000, section 5.4) JR FOX THIS: DONE




inline double CalculateCOCoreMass_Hurley2000() const override { return CalculateCOCoreMass_Hurley2000(Age(), m_InterimState.HurleyGBParams(), m_InterimState.HurleyTimescales()); }
GNU_CONST double CalculateCOCoreMass_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const;

inline double CalculateHeCoreMass() const override { return HeCoreMass(); } // McHe is constant for EAGB stars


GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd_Hurley2000() const override { return CurrentState().HurleyGBParams(HURLEY_GB_PARAMETERS::McDU); }    // McCO = McDU at phase end for EAGB stars, per Hurley et al. 2000, section 5.4

inline double CalculateHeCoreMassAtPhaseEnd() const override { return HeCoreMass(); } // McHe is constant for EAGB stars






    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////

    GNU_PURE static double CalculateRadius_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Luminosity);




    //////////////////////////////////////////////////
    //   REMNANT                                    //
    //////////////////////////////////////////////////

    GNU_CONST double CalculateRemnantRadius_Hurley2000(const double p_Radius, const double p_HeCoreMass) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateRemnantRadius_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateRemnantRadius_Hurley2000(Radius(), HeCoreMass());
    }





///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    double          CalculateInitialSupernovaMass() const                                           { return m_InterimState.HurleyGBParams(HURLEY_GB_PARAMETERS::McBAGB); }

    double          CalculateLifetimeTo2ndDredgeUp(const double p_Tinf1_FAGB, const double p_Tinf2_FAGB) const { } // JR FIX THIS: placeholder

    double          CalculateLuminosityAtPhaseEnd(const double p_CoreMass) const                    { return CalculateLuminosityOnPhase(p_CoreMass); }                  // Same as on phase
    double          CalculateLuminosityAtPhaseEnd() const                                           { return CalculateLuminosityAtPhaseEnd(COCoreMass()); }
    double          CalculateLuminosityOnPhase(const double p_CoreMass) const { } // JR FIX THIS: placeholder
    double          CalculateLuminosityOnPhase() const                                              { return CalculateLuminosityOnPhase(COCoreMass()); }



COMPAS_PURE MassLossT CalculateMLRate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_Perturb) const override;



    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const { return CalculateRadiusOnPhase(p_Mass, p_Luminosity); }            // Same as on phase
    double          CalculateRadiusAtPhaseEnd() const                                               { return CalculateRadiusAtPhaseEnd(Mass(), Luminosity()); }
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return CalculateRadius_Hurley2000_Static(p_Mass, p_Luminosity, ZDEP->HurleyBCoefficients(Metallicity()), ZDEP->HurleyMassCutoffs(Metallicity(), HURLEY_MCO::HeF)); }
    double          CalculateRadiusOnPhase() const                                                  { return CalculateRadiusOnPhase(Mass(), Luminosity()); }

    double          CalculateRemnantLuminosity() const;
    double          CalculateRemnantRadius() const {  } // JR FIX THIS: placeholder







GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const override;


GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::CONVECTIVE; } // Always CONVECTIVE for EAGB stars


    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                                            { return !ShouldEvolveOnPhase(); }                                  // Phase ends when age at or after DU timescale, and no TPAGB
    bool            IsSupernova() const;

 
    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash() {  }  // JR FIX THIS: placeholder

    bool            ShouldEvolveOnPhase() const;
    bool            ShouldSkipPhase() const;

    COMPAS_PURE double CalculateLifetimeTo2ndDU_Hurley2000(const DBL_VECTOR& p_GBParams, const double p_Tinf1_FAGB, const double p_Tinf2_FAGB) const;

    COMPAS_PURE double CalculateCELambda_Nanjing_Enhanced(const double p_Mass, const double p_Radius, const double p_CoreMass, const SizeT p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;

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
 * double CalculateLifetimeTo2ndDU_Hurley2000(const DBL_VECTOR& p_GBParams, const double p_Tinf1_FAGB, const double p_Tinf2_FAGB)
 *
 * @param       p_GBParams                      Hurley GB parameters
 * @param       p_tInf1_FAGB                    EAGB integration constant 1, tinf1_FAGB (Myr)
 * @param       p_tInf2_FAGB                    EAGB integration constant 2, tinf2_FAGB (Myr)
 * @return                                      Lifetime to second dredge up (tDU) (Myr)
 */
inline double EAGB::CalculateLifetimeTo2ndDU_Hurley2000(const DBL_VECTOR& p_GBParams, const double p_Tinf1_FAGB, const double p_Tinf2_FAGB) const {
#define GBParams(x) p_GBParams[static_cast<int>(HURLEY_GB_PARAMETERS::x)]

    const double lDU = CalculateLuminosity_Hurley2000(GBParams(McDU), p_GBParams);
    const double p1  = GBParams(p) - 1.0;
    const double q1  = GBParams(q) - 1.0;

    return lDU <= GBParams(Lx)
            ? p_Tinf1_FAGB - (1.0 / (p1 * GBParams(AHe) * GBParams(D))) * PPOW((GBParams(D) / lDU), (p1 / GBParams(p)))
            : p_Tinf2_FAGB - (1.0 / (q1 * GBParams(AHe) * GBParams(B))) * PPOW((GBParams(B) / lDU), (q1 / GBParams(q)));

#undef GBParams
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
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Suggested timestep (Myr)
 */
GNU_CONST inline double EAGB::ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const {
#define tScales(x) p_tScales[static_cast<int>(HURLEY_TIMESCALES::x)]

    const double dtk = 0.02 * ((p_Age <= tScales(Mx_FAGB) ? tScales(Inf1_FAGB) : tScales(Inf2_FAGB)) - p_Age);   // stellar type specific dt

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
 * double CalculateRemnantLuminosity_Hurley2000(const double p_COCoreMass, const GBParamsT& p_GBParams) const
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_GBParams                      GB parameters
 * @return                                      Remnant core luminosity (Lsol)
 */
inline double CalculateRemnantLuminosity_Hurley2000(const double p_COCoreMass, const GBParamsT& p_GBParams) const {
    return HeGB::CalculateLuminosityOnPhase(p_COCoreMass, p_GBParams);
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
 * double CalculateCOCoreMass_Hurley2000(const double p_Time, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales)
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_GBParams                      Hurley GB parameters
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      EAGB CO core mass (Msol)
 */
GNU_CONST inline double EAGB::CalculateCOCoreMass_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const {
// macros for convenience and readability - undefined at end of function
#define GBParams(x) p_GBParams[static_cast<int>(HURLEY_GB_PARAMETERS::x)]
#define tScales(x) p_tScales[static_cast<int>(HURLEY_TIMESCALES::x)]

    return p_Age <= tScales(Mx_FAGB)
            ? PPOW((GBParams(p) - 1.0) * GBParams(AHe) * GBParams(D) * (tScales(Inf1_FAGB) - p_Age), 1.0 / (1.0 - GBParams(p)))
            : PPOW((GBParams(q) - 1.0) * GBParams(AHe) * GBParams(B) * (tScales(Inf2_FAGB) - p_Age), 1.0 / (1.0 - GBParams(q)));

#undef tScales
#undef GBParams
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
 * double CalculateRemnantRadius_Hurley2000(const double p_Radius, const double p_HeCoreMass) const
 *
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @return                                      Remnant core radius (Rsol)
 */
GNU_CONST double EAGB::CalculateRemnantRadius_Hurley2000(const double p_Radius, const double p_HeCoreMass) const {
    const double lum = CalculateRemnantLuminosity_Hurley2000();
    return std::min(p_Radius, HeGB::CalculateRadius_Hurley2000(p_HeCoreMass, lum));
}



















///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    MT_CASE DetermineMassTransferTypeAsDonor() const override { return MT_CASE::C; }                                                                                                                                            // Always case C
    double  CalculateCriticalMassRatio_Hurley2002() const override { return GiantBranch_Constituent::CalculateCriticalMassRatio_Hurley2002(); }

    COMPAS_PURE double  CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const override;                 // body in EAGB.cpp



