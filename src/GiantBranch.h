#ifndef __GiantBranch_h__
#define __GiantBranch_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Rand.h"
#include "MainSequence.h"


class BaseStar;
class MainSequence;

class GiantBranch: virtual public BaseStar, public MainSequence {

public:

    GiantBranch(){};
    
    GiantBranch(const BaseStar &p_BaseStar) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) {}
    


protected:













    

    virtual double          CalculateInitialSupernovaMass() const                                           { return m_Mass; }                                                  // Use class member variables

   







            double          CalculateMassLossRateHurley();






            
            
            
            double          CalculateFallbackFraction_Fryer2012_Delayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass);


            double          CalculateFallbackFractionRapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass);



            
            
            double          CalculateGravitationalRemnantMass(const double p_BaryonicRemnantMass);




            double          CalculateRemnantMassByMullerMandel(const double p_COCoreMass, const double p_HeCoreMass);
            
            double          CalculateRemnantMassBySchneider2020(const double p_COCoreMass, const bool p_UseSchneiderAlt = false);
            double          CalculateRemnantMassBySchneider2020Alt(const double p_COCoreMass)               { return CalculateRemnantMassBySchneider2020(p_COCoreMass, true); }
            double          CalculateRemnantNSMassMullerMandel(const double p_COCoreMass, const double p_HeCoreMass);








            double          CalculateConvectiveEnvelopeRadialExtent() const;
            double GiantBranch::CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Radius) const






















	





            double          CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription);
            double          CalculateZetaConvectiveEnvelopeGiant(ZETA_PRESCRIPTION p_ZetaPrescription);

    virtual void            PerturbLuminosityAndRadius();

            void            ResolveSupernovaPreamble();
            STELLAR_TYPE    ResolveSupernova();
            STELLAR_TYPE    ResolveCoreCollapseSN();
            STELLAR_TYPE    ResolveElectronCaptureSN();
            STELLAR_TYPE    ResolvePairInstabilitySN();
            STELLAR_TYPE    ProcessPPISN();
    
            void            UpdateAgeAfterMassLoss() { }                                                                                                                        // NO-OP for most stellar types


}
    
            void            UpdateMainSequenceCoreMass(const double p_Dt, const double p_TotalMassLossRate) { }                                                                 // NO-OP for most stellar types




///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


DBL_VECTOR CalculateTimescales_Hurley2000() const override { return CalculateTimescales_Hurley2000(MassEffectiveInitial(), m_InterimState.GBparams(), m_InterimState.TimeScales()); }
COMPAS_PURE DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;

DBL_VECTOR CalculateGBparams_Hurley2000() const override { return CalculateGBparams_Hurley2000(MassEffectiveInitial(), m_InterimState.GBparams()); }
COMPAS_PURE virtual DBL_VECTOR CalculateGBparams_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const;








 







    // Member functions (not getters or setters)

    ////////////////////////////////////////
    //   VIRTUAL FUNCTIONS                //
    ////////////////////////////////////////
     
    // By function, then alphabetical.
    // May be (are expected to be) overridden by derived classes.
    // When overriding virtual functions in a derived class, use the "override" attribute.

    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////

    virtual double CalculateRemnantLuminosity_Hurley2000(const double p_Mass, const double p_CoreMass) const;


    ////////////////////////////////////////
    //   MISC. CONSTANTS / PARAMETERS     //
    ////////////////////////////////////////

    virtual double CalculateCoreMass_Luminosity_B_Hurley2000(const double p_Mass) const;
    virtual double CalculateCoreMass_Luminosity_D_Hurley2000(const double p_Mass) const;
    virtual double CalculateCoreMass_Luminosity_p_Hurley2000(const double p_Mass) const;
    virtual double CalculateCoreMass_Luminosity_q_Hurley2000(const double p_Mass) const;










    ////////////////////////////////////////
    //   NON-VIRTUAL FUNCTIONS            //
    ////////////////////////////////////////

    // Should not be overridden (declared separately) by derived classes.
    // While it is legal in C++ to declare the same (non-virtual) function in multiple classes,
    // (aka "shadowing", or "hiding"), we discourage it.  Non-virtual functions are statically
    // linked, and as such, especially with indirection, may not produce expected results.


    ////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES   //
    ////////////////////////////////////////

    double CalculateLifetimeToHeI_Hurley2000(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB) const;


    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////

    
    COMPAS_PURE static double CalculateLuminosityAtBGB_Hurley2000_Static(const double p_Mass);
    COMPAS_PURE static double CalculateLuminosityAtHeI_Hurley2000_Static(const double p_Mass);

    COMPAS_PURE double CalculateLuminosityOnZAHB_Hurley2000(const double p_Mass, const double p_CoreMass) const;


    ////////////////////////////////////////
    //   MASS                             //
    ////////////////////////////////////////

    GNU_CONST   double   CalculateBaryonicRemnantMass_Fryer2012(const double p_ProtoMass, double p_FallbackMass) const;

                double   CalculateBHMassAfterFallback_MullerMandel2020(const double p_COCoreMass, const double p_HeCoreMass) const;

    inline      double   CalculateConvectiveCoreMass() const override { return CoreMass(); }

    inline      Dbl_DblT CalculateConvectiveEnvelopeMass() const override { return CalculateConvectiveEnvelopeMass_Picker2024(Mass(), Mass0(), Temperature()); }
    COMPAS_PURE Dbl_DblT CalculateConvectiveEnvelopeMass_Picker2024(const double p_Mass, const double p_MassEffectiveInitial, const double p_Temperature) const;

    GNU_CONST   double   CalculateCoreMassAt2ndDredgeUp_Hurley2000(const double p_McBAGB) const;
    COMPAS_PURE double   CalculateCoreMassAtBAGB_Hurley2000(const double p_Mass) const;
    COMPAS_PURE double   CalculateCoreMassAtBGB_Hurley2000(const double p_Mass, const DblVectorT& p_GBparams) const;
    COMPAS_PURE double   CalculateCoreMassAtHeI_Hurley2000(const double p_Mass) const;
    GNU_CONST   double   CalculateCoreMassAtSN(const double p_Mthreshold, const double p_McBAGB);

    inline      double   CalculateEffectiveInitialMass_Hurley2000() const override { return BaseStar::CalculateEffectiveInitialMass_Hurley2000(); } // per Hurley et al. 2000, section 7.1

    GNU_CONST   double   CalculateFallbackMass_Fryer2012(const double p_PreSNMass, const double p_ProtoMass, const double p_FallbackFraction) const;

    GNU_CONST   double   CalculateProtoCoreMass_Fryer2012_Delayed(const double p_COCoreMass) const;

    GNU_CONST   Dbl_DblT CalculateRemnantMass_Belczynski2002(const double p_Mass, const double p_COCoreMass) const;
    COMPAS_PURE Dbl_DblT CalculateRemnantMass_Fryer2012(const double p_Mass, const double p_COCoreMassn) const;
    COMPAS_PURE Dbl_DblT CalculateRemnantMass_Fryer2022(const double p_Mass, const double p_COCoreMass, const double p_fMix, const double p_mCrit) const;
                Dbl_DblT CalculateRemnantMass_Maltsev2025(const double p_COCoreMass, const double p_HeCoreMass);
    GNU_CONST   double   CalculateRemnantMass_Muller2016(const double p_Mass, const double p_COCoreMass) const;


    ////////////////////////////////////////
    //   MASS LOSS                        //
    ////////////////////////////////////////

    inline double CalculateMLRateThermal() const { return (Mass() - CoreMass()) / CalculateTimescale_Thermal(); }


    ////////////////////////////////////////
    //   MISC. CONSTANTS / PARAMETERS     //
    ////////////////////////////////////////

    GNU_CONST double CalculateCoreMass_Luminosity_Lx_Hurley2000(const DblVectorT& p_GBparams) const;
    GNU_CONST double CalculateCoreMass_Luminosity_Mx_Hurley2000(const DblVectorT& p_GBparams) const;

    GNU_CONST double CalculateHRateConstant_Hurley2000(const double p_Mass) const;
    GNU_CONST double CalculatePerturbationMu_Hurley2000(const double p_Mass, const double p_Luminosity, const double p_CoreMass) const;


    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////

    inline double CalculateConvectiveCoreRadius() const override { return std::min(CalculateRemnantRadius(), Radius()); } // Last paragraph of section 6 of Hurley+ 2000

    double CalculateRadius(const double p_Mass, const double p_Luminosity) const;

    COMPAS_PURE double CalculateRadiusAtHeI_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity) const;

    inline double CalculateRadiusOnMassChange(double p_dM) const override { return CalculateRadius(Mass() + p_dM, Luminosity()); }

    COMPAS_PURE CalculateRadiusOnZAHB_Hurley2000(const double p_Mass, const double p_CoreMass) const;

    inline double CalculateRadius_Hurley2000() const override { return CalculateRadius_Hurley2000(Mass(), Luminosity()); }
    COMPAS_PURE double CalculateRadius_Hurley2000(const double p_Mass, const double p_Luminosity) const;


    ////////////////////////////////////////
    //   REMNANTS                         //
    ////////////////////////////////////////

    inline double CalculateRemnantRadius_Hurley2000() const override { return CalculateRemnantRadius_Hurley2000(Mass0(), CoreMass()); }
    COMPAS_PURE double CalculateRemnantRadius_Hurley2000(const double p_Mass, const double p_CoreMass) const;

    GNU_CONST STELLAR_TYPE GiantBranch::CalculateRemnantType_Muller2016(const double p_COCoreMass) const;


    ////////////////////////////////////////
    //   ROTATION                         //
    ////////////////////////////////////////

    COMPAS_PURE double CalculateMomentOfInertia(const double p_Mass, const double p_Radius, const double p_CoreMass) const;

















///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<






};


// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLifetimeToHeI_Hurley2000
 *
 * @brief
 * Calculate the time to Helium ignition tHeI, per Hurley et al. 2000, eq 43
 *
 *
 * double CalculateLifetimeToHeI_Hurley2000(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB) const
 *
 *
 * @param       p_Mass                      Mass of the star (Msol)
 * @param       p_Tinf1_FGB                 tinf1_FGB (Timescales[TIMESCALE::tinf1_FGB]) - First Giant Branch tinf1
 * @param       p_Tinf2_FGB                 tinf2_FGB (Timescales[TIMESCALE::tinf2_FGB]) - First Giant Branch tinf2
 * @param       p_GBparams                      Hurley GB parameters

 * @return                                      Lifetime to He ignition (tHeI)
 */
inline double GiantBranch::CalculateLifetimeToHeI_Hurley2000(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB, const DblVectorT& p_GBparams) const {

const DblVectorT& p_GBparams) const {

    const double AH = p_GBparams[HURLEY_GBP::AH];
    const double B  = p_GBparams[HURLEY_GBP::B];
    const double D  = p_GBparams[HURLEY_GBP::D];
    const double p  = p_GBparams[HURLEY_GBP::P];
    const double q  = p_GBparams[HURLEY_GBP::Q];
    const double p1 = p - 1.0;
    const double q1 = q - 1.0;

    const double lHeI = CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass);

    return lHeI <= gbParams(Lx)
            ? p_Tinf1_FGB - (1.0 / (p1 * AH * D)) * PPOW((D / lHeI), (p1 / p))
            : p_Tinf2_FGB - (1.0 / (q1 * AH * B)) * PPOW((B / lHeI), (q1 / q));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateBHMassAfterFallback_MullerMandel2020
 *
 * @brief
 * Calculate the mass of the black hole created in a fallback supernova, per Mandel & Mueller, 2020
 * 
 * Also used by Maltsev2025.
 *
 *
 * double CalculateBHMassAfterFallback_MullerMandel2020(const double p_COCoreMass, const double p_HeCoreMass) const
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @return                                      Remnant mass (Msol)
 */
inline double GiantBranch::CalculateBHMassAfterFallback_MullerMandel2020(const double p_COCoreMass, const double p_HeCoreMass) const {
    
    double remnantMass;

    SizeT iterations = 0;
    while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS &&
          (remnantMass < OPTIONS->MaximumNeutronStarMass() || remnantMass > p_HeCoreMass)) {
        remnantMass = MULLERMANDEL_MUBH * p_COCoreMass + RAND->RandomGaussian(MULLERMANDEL_SIGMABH);
    }
    
    if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS) {   // Converged?
                                                                    // No
        // Failure to find a solution implies a narrow range - just pick a midpoint in this case
        remnantMass = (OPTIONS->MaximumNeutronStarMass() + p_HeCoreMass) / 2.0;
    }
    
    return remnantMass;
}


/*
 * CalculateConvectiveEnvelopeMass_Picker2024
 *
 * @brief
 * Calculate the mass of the convective envelope.
 * 
 * Approximates the mass of the outer convective envelope.
 * Follows the fits of Picker, Hirai, Mandel 2024 (see https://arxiv.org/pdf/2402.13180)
 * 
 * See also Mandel, Hirai, Picker 2024 (see https://arxiv.org/pdf/2412.10691)
 *
 *
 * Dbl_DblT CalculateConvectiveEnvelopeMass_Picker2024(const double p_Mass, const double p_MassEffectiveInitial, const double p_Temperature) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_MassEffectiveInitial          Effective initial mass of the star (Msol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                   GB convective envelope mass (Msol)
 *                                                   Maximum GB convective envelope mass (Msol)
 */
inline Dbl_DblT GiantBranch::CalculateConvectiveEnvelopeMass_Picker2024(const double p_Mass, const double p_MassEffectiveInitial, const double p_Temperature) const {
    
    // Ratio of final intershell mass to final core mass - Picker at al. 2024 eq 8.
    const double MiFinalMcFinal = -0.023 * GLOBALS->SigmaHurley() - 0.0023;

    // We need the temperature of the star just after BAGB, which is the temperature at the
    // start of the EAGB phase.  Since we are on the giant branch here, we can clone this
    // object as an EAGB object and, as long as it is initialised (to the start of the phase),
    // we can query the cloned object for its temperature.
    //
    // To ensure the clone does not participate in logging, we set its persistence to EPHEMERAL.
    std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::EAGB, OBJECT_PERSISTENCE::EPHEMERAL);                   // Clone
    const double tMin = clone->Temperature();                                                                       // Get temperature of clone

    // Use Eq. 6 of Mandel et al. 2024 eq 6 rather than Picker at al. 2024 eq 6 for tOnset
    // to avoid issues caused by differences between temperatures in MESA models (used in
    // fits from Picker et al. 2024) and Pols models (used in Hurley SSE tracks).
    const double tOnset  = tMin / std::min(0.695 - 0.057 * GLOBALS->SigmaHurley(), 0.95);                           // Mandel et al. 2024, eq 6
    const double McFinal = CalculateCoreMassAtBAGB_Hurley2000(p_MassEffectiveInitial);
    
    double mEnvMax = std::max(p_Mass - McFinal * (1.0 + MiFinalMcFinal), 0.0);                                      // Picker at al. 2024, eq 9

    // Picker+ 2024 fits were only made for stars above 8.0 solar masses, with runs down to 5.0 Msol,
    // so we use the final core mass as an approximate threshold of validity.
    // Unlike massive stars, intermediate-mass stars have almost no radiative intershell at maximum
    // convective envelope extent.
    if (McFinal < 1.5) mEnvMax = std::max(p_Mass - mCoreFinal, 0.0);
    
    const double mEnv = mEnvMax / (1.0 + std::exp(4.6 * (tMin + tOnset - 2.0 * p_Temperature) / (tMin - tOnset)));  // Picker at al. 2024, eq 7
    
    return std::make_tuple(mEnv, mEnvMax);
}


/*
 * CalculateCoreMassAt2ndDredgeUp_Hurley2000
 *
 * @brief
 * Calculate the mass of the core during second dredge up, per Hurley et al. 2000, eq 69
 * (not numbered in the paper, but between eqs 68 & 70)
 *
 *
 * double CalculateCoreMassAt2ndDredgeUp_Hurley2000(const double p_McBAGB) const
 *
 * @param       p_McBAGB                        Core mass at the Base of the Asymptotic Giant Branch (Msol)
 * @return                                      Core Mass at second dredge up (Msol)
 */
inline double GiantBranch::CalculateCoreMassAt2ndDredgeUp_Hurley2000(const double p_McBAGB) const {
    return p_McBAGB >= 0.8 ? ((0.44 * p_McBAGB) + 0.448) : p_McBAGB;
}


/*
 * CalculateCoreMassAtBAGB_Hurley2000
 *
 * @brief
 * Calculate core mass at the Base of the Asymptotic Giant Branch,
 * per Hurley et al. 2000, eq 66
 *
 *
 * double CalculateCoreMassAtBAGB_Hurley2000(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      BAGB core mass (Msol)
 */
inline double GiantBranch::CalculateCoreMassAtBAGB_Hurley2000(const double p_Mass) const {
    const DBL_VECTOR b = GLOBALS->HurleyBcoefficients();                    // get Hurley b coefficients
    return std::sqrt(std::sqrt((b[36] * PPOW(p_Mass, b[37])) + b[38]));     // sqrt() is much faster than PPOW()
}


/*
 * CalculateCoreMassAtBGB_Hurley2000
 *
 * @brief
 * Calculate core mass at the Base of the Giant Branch, per Hurley et al. 2000, eq 44
 *
 * For large enough M, we have McBGB ~ 0.098 * Mass^1.35
 *
 *
 * double CalculateCoreMassAtBGB_Hurley2000(const double p_Mass, const double p_McBAGB) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_McBAGB                        BAGB core mass (Msol)
 * @return                                      BGB core mass (Msol)
 */
double GiantBranch::CalculateCoreMassAtBGB_Hurley2000(const double p_Mass, const double p_McBAGB) const {

    double McBGB = 0.0;     // default value for McBGB (see comment below)

    const double mHeF = GLOBALS-HurelyMassCutoffs(static_cast<int>(MHeF));

    // no McBGB for stars with mass below the helium flash threshold
    // See Hurley at al. 2000, text immediately prior to eq 44
    if (p_Mass >= mHeF) { 
        const double luminosity = GiantBranch::CalculateLuminosityAtBGB_Hurley2000(mHeF);
        const double Mc_MHeF    = BaseStar::CalculateCoreMass_Hurley2000_Static(luminosity, p_GBparams);
        const double c          = utils::IntPow(Mc_MHeF, 4) - MC_L_C1 * PPOW(mHeF, MC_L_C2);
    
        McBGB = std::min(0.95 * p_GBparams(McBAGB), std::sqrt(std::sqrt(c + MC_L_C1 * PPOW(p_Mass, MC_L_C2))));
    }

    return McBGB;
}


/*
 * CalculateCoreMassAtHeI_Hurley2000
 *
 * @brief
 * Calculate core mass at Helium Ignition (HeI), per:
 * 
 *    - Hurley et al. 2000, eq 37 for low mass stars (p_Mass < p_MHef, see p_MHeF below)
 *    - Hurley et al. 2000, eq 44 otherwise, replacing Mc(LBGB(MHeF)) with Mc(LHeI(MHeF)),
 *                                           per Hurley et al. 2000 end of section 5.3
 *
 *
 * double CalculateCoreMassAtHeI_Hurley2000(const double p_Mass, const DblVectorT& p_GBparams) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Core mass at HeI (Msol)
 */
inline double GiantBranch::CalculateCoreMassAtHeI_Hurley2000(const double p_Mass, const DblVectorT& p_GBparams) const {

    double coreMass;

    if (p_Mass < p_MHef) {
        const double luminosity = CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass);
        coreMass = BaseStar::CalculateCoreMass_Hurley2000_Static(luminosity, p_GBparams);
    }
    else {
        const double mHeF    = GLOBALS->HurleyMassCutoffs(HURLY_MCO::HEF);
        const double lMHeF   = CalculateLuminosityAtHeI_Hurley2000_Static(mHeF);
        const double Mc_MHeF = BaseStar::CalculateCoreMass_Hurley2000(lMHeF, p_GBparams);
        const double McBAGB  = CalculateCoreMassAtBAGB_Hurley2000(p_Mass);
        const double c       = (utils::intPow(Mc_MHeF, 4)) - (MC_L_C1 * PPOW(mHeF, MC_L_C2));

        coreMass = std::min((0.95 * McBAGB), std::sqrt(std::sqrt(c + (MC_L_C1 * PPOW(p_Mass, MC_L_C2))))); // sqrt() is much faster than PPOW()
    }

    return coreMass;
}


/*
 * CalculateCoreMassAtSN
 *
 * @brief
 * Calculate the core mass at which the Asymptotic Giant Branch phase
 * is terminated in a SN/loss of envelope, per Hurley et al. 2000, eq 75
 * 
 * The mass returned is clamped to a minimum CO core mass (passed as p_mThreshold), which
 * is MCH in Hurley et al. eq 75, but we deviate from Hurley for stars that may explode in
 * ECSNe and use MECS (instead of MCH) in those cases.
 *
 *
 * static double CalculateCoreMassAtSN(const double p_mThreshold, const double p_McBAGB)
 *
 * @param       p_mThreshold                    Threshold mass (see notes above) (Msol)
 * @param       p_McBAGB                        Core mass at the Base of the Asymptotic Giant Branch (Msol)
 * @return                                      Maximum core mass pre-SN on the Asymptotic Giant Branch (Msol)
 */
inline double GiantBranch::CalculateCoreMassAtSN(const double p_mThreshold, const double p_McBAGB) {
    return std::max(p_mThreshold, (0.773 * p_McBAGB) - 0.35);
}


/*
 * CalculateFallbackMass_Fryer2012
 *
 * @brief
 * Calculate the mass falling back onto the proto compact object,
 * per Fryer et al. 2012, eq 11
 *
 *
 * double CalculateFallbackMass_Fryer2012(const double p_PreSNMass, const double p_ProtoMass, const double p_Fallback) const
 *
 * @param       p_PreSNMass                     Pre-SN mass of the star (Msol)
 * @param       p_ProtoMass                     Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @param       p_FallbackFraction              Fraction of mass falling back onto proto object [0.0, 1.0]
 * @return                                      Mass falling back onto proto object (Msol)
 */
inline double GiantBranch::CalculateFallbackMass_Fryer2012(const double p_PreSNMass, const double p_ProtoMass, const double p_FallbackFraction) const {
    return p_FallbackFraction * (p_PreSNMass - p_ProtoMass);
}


/*
 * CalculateProtoCoreMass_Fryer2012_Delayed
 *
 * @brief
 * Calculate the mass of the proto core, using the delayed supernova mechanism,
 * per Fryer et al. 2012, eq 18.
 *
 *
 * double CalculateProtoCoreMass_Fryer2012_Delayed(const double p_COCoreMass) const
 *
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @return                                      Mass of the Fe/Ni proto core (Msol)
 */
inline double GiantBranch::CalculateProtoCoreMass_Fryer2012_Delayed(const double p_COCoreMass) const {
    return p_COCoreMass < 3.5 ? 1.2 : (p_COCoreMass < 6.0 ? 1.3 : (p_COCoreMass < 11.0 ? 1.4 : 1.6));
}


/*
 * CalculateRemnantMass_Belczynski2002
 *
 * @brief
 * Calculate remnant mass, per Belczynski et al. 2002
 *
 * This is also used in Hurley SSE code
 * (equation from Belczynski et al. 2002, not documented in Hurley et al. 2000)
 *
 *
 * Dbl_DblT CalculateRemnantMass_Belczynski2002(const double p_Mass, const double p_COCoreMass) const
 *
 * @param       p_Mass                          Pre-SN mass of the star (Msol)
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
inline Dbl_DblT GiantBranch::CalculateRemnantMass_Belczynski2002(const double p_Mass, const double p_COCoreMass) const {
    const double fallbackFraction = p_COCoreMass <= 5.0 ? 0.0 : (p_COCoreMass < 7.6 ? (p_COCoreMass - 5.0) / 2.6 : 1.0);
    const double McFeNi           = p_COCoreMass < 2.5 ? (0.161767 * p_COCoreMass) + 1.067055 : (0.314154 * p_COCoreMass) + 0.686088; // iron core mass
    return std::make_tuple(McFeNi + (fallbackFraction * (p_Mass - McFeNi)), fallbackFraction);
}


/*
 * CalculateBaryonicRemnantMass_Fryer2012
 *
 * @brief
 * Calculate the baryonic mass of the remnant, per Fryer et al. 2012, eq 12
 *
 *
 * double CalculateBaryonicRemnantMass_Fryer2012(const double p_ProtoMass, double p_FallbackMass) const
 *
 * @param       p_ProtoMass                     Mass of proto compact object (Msol)
 * @param       p_FallbackMass                  Mass falling back onto proto compact object (Msol)
 * @return                                      Baryonic mass of the remnant (Msol)
 */
inline double GiantBranch::CalculateBaryonicRemnantMass_Fryer2012(const double p_ProtoMass, double p_FallbackMass) const {
    return p_ProtoMass + p_FallbackMass;
}


/*
 * CalculateRemnantMass_Muller2016
 *
 * @brief
 * Calculate remnant mass, per Muller et al. 2016
 * (as presented in eq B4 of Vigna-Gomez et al. 2018 (See arXiv:1805.07974))
 *
 * 
 * double CalculateRemnantMass_Muller2016(const double p_Mass, const double p_COCoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
inline double GiantBranch::CalculateRemnantMass_Muller2016(const double p_Mass, const double p_COCoreMass) const {

    double	remnantMass;

    if (p_COCoreMass < 1.372) {
        // Not explicitly pointed out in Appendix B of Vigna-Gomez et al. 2018,
        // but assumed for continuity and simplicity.
        // Muller et al. 2016 didn't go as low as this in CO Core mass (see Figure A1)
        remnantMass = 1.21;
    }                                      
	else if (p_COCoreMass < 1.49) remnantMass = 1.21 - (0.4  * (p_COCoreMass - 1.372));
	else if (p_COCoreMass < 1.65) remnantMass = 1.16;
    else if (p_COCoreMass < 2.4 ) remnantMass = 1.32 + (0.3  * (p_COCoreMass - 1.65));
    else if (p_COCoreMass < 3.2 ) remnantMass = 1.42 + (0.7  * (p_COCoreMass - 2.4));
    else if (p_COCoreMass < 3.6 ) remnantMass = 1.32 + (0.25 * (p_COCoreMass - 3.2));
    else if (p_COCoreMass < 4.05) remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; // -> Black Hole
    else if (p_COCoreMass < 4.6 ) remnantMass = 1.5;
    else if (p_COCoreMass < 5.7 ) remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; // -> Black Hole
    else if (p_COCoreMass < 6.0 ) remnantMass = 1.64 - (0.2  * (p_COCoreMass - 5.7));
    else                          remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; // -> Black Hole

    return remnantMass;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           MISC. CONSTANTS / PARAMETERS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateHRateConstant_Hurley2000
 *
 * @brief
 * Calculate the Hydrogen rate constant AH', per Hurley et al. 2000,
 * just after eq 43 (before eq 44)
 *
 *
 * double CalculateHRateConstant_Hurley2000(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Hydrogen rate constant (Msol Lsol^-1 Myr^-1)
 */
inline double GiantBranch::CalculateHRateConstant_Hurley2000(const double p_Mass) const {
    return PPOW(10.0, std::max(-4.8, std::min((-5.7 + (0.8 * p_Mass)), (-4.1 + (0.14 * p_Mass)))));
}


/*
 * CalculatePerturbationMu_Hurley2000
 *
 * @brief
 * Calculate the Hurley small envelope perturbation parameter, mu,
 * per Hurley et al. 2000 eqs 97 & 98
 *
 *
 * double CalculatePerturbationMu_Hurley2000(const double p_Mass, const double p_Luminosity, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Small envelope perturbation parameter, mu
 */
inline double GiantBranch::CalculatePerturbationMu_Hurley2000(const double p_Mass, const double p_Luminosity, const double p_CoreMass) const {
    constexpr double kappa = -0.5;
    constexpr double L0    = 7.0E4;
    return ((p_Mass - p_CoreMass) / p_Mass) * (std::min(5.0, std::max(1.2, PPOW((p_Luminosity / L0), kappa))));
}


/*
 * CalculateCoreMass_Luminosity_Lx_Hurley2000
 *
 * @brief
 * Calculate First Giant Branch (FGB) Core mass - Luminosity relation parameter, Lx,
 * per Hurley et al. 2000, eq 37
 *
 *
 * double CalculateCoreMass_Luminosity_Lx_Hurley2000(const DblVectorT& p_GBparams) const
 *
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Core mass - Luminosity relation parameter, Lx
 */
inline double GiantBranch::CalculateCoreMass_Luminosity_Lx_Hurley2000(const DblVectorT& p_GBparams) const {

    const double B  = p_GBparams[HURLEY_GBP::B];
    const double D  = p_GBparams[HURLEY_GBP::D];
    const double Mx = p_GBparams[HURLEY_GBP::MX];
    const double p  = p_GBparams[HURLEY_GBP::P];
    const double q  = p_GBparams[HURLEY_GBP::Q];

    // since the mass used here is the mass at crossover (Mx), these
    // should give the same answer - but we'll take the minimum anyway
    return std::min((B * PPOW(Mx, q)), (D * PPOW(Mx, p)));
}


/*
 * CalculateCoreMass_Luminosity_Mx_Hurley2000
 *
 * @brief
 * Calculate the Core mass - Luminosity relation parameter, Mx,
 * per Hurley et al. 2000, eq 38
 *
 * Mx is the point at which the low- and high-luminosity approximations cross
 *
 *
 * double CalculateCoreMass_Luminosity_Mx_Hurley2000(const DblVectorT& p_GBparams) const
 *
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Core mass - Luminosity relation parameter, Mx
 */
inline double GiantBranch::CalculateCoreMass_Luminosity_Mx_Hurley2000(const DblVectorT &p_GBparams) const {

    const double B  = p_GBparams[HURLEY_GBP::B];
    const double D  = p_GBparams[HURLEY_GBP::D];
    const double p  = p_GBparams[HURLEY_GBP::P];
    const double q  = p_GBparams[HURLEY_GBP::Q];

    return PPOW(B / D, (1.0 / (p - q)));
}


/*
 * CalculateCoreMass_Luminosity_B_Hurley2000
 *
 * @brief
 * Calculate the Core mass - Luminosity relation parameter, B, 
 * per Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_B_Hurley2000(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Core mass - Luminosity relation parameter, B
 */
inline double GiantBranch::CalculateCoreMass_Luminosity_B_Hurley2000(const double p_Mass) const {
    return std::max(3.0E4, (500.0 + (1.75E4 * PPOW(p_Mass, 0.6))));
}


/*
 * CalculateCoreMass_Luminosity_D_Hurley2000
 *
 * @brief
 * Calculate the Core mass - Luminosity relation parameter D, per Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_D_Hurley2000(onst double p_Mass) const 
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Core mass - Luminosity relation parameter D
 */
inline double GiantBranch::CalculateCoreMass_Luminosity_D_Hurley2000(const double p_Mass) const {

    const double D0 = 5.37 + (0.135 * GLOBALS->ZetaHurley());
    const double D1 = (0.975 * D0) - (0.18 * p_Mass);

    double logD = D0;   // default value

    const double mHeF = GLOBALS->HurleyMassCutoffs(HURLEY_MCO::HEF);

    if (p_Mass > mHeF) {
        if (p_Mass >= 2.5) {
            logD = std::max(std::max(-1.0, D1), (0.5 * D0) - (0.06 * p_Mass));
        }
        else {          // Linear interpolation between end points
            const double gradient  = (D0 - D1) / (mHeF - 2.5);
            const double intercept = D0 - (mHeF * gradient);
            logD = (gradient * p_Mass) + intercept;
        }
    }

    return PPOW(10.0, logD);
}


/*
 * CalculateCoreMass_Luminosity_p_Hurley2000
 *
 * @brief
 * Calculate the Core mass - Luminosity relation parameter, p, per Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_p_Hurley2000(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Core mass - Luminosity relation parameter, p
 */
inline double GiantBranch::CalculateCoreMass_Luminosity_p_Hurley2000(const double p_Mass) const {

    double p = 6.0;                                             // default value

    const double mHeF = GLOBALS->HurleyMassCutoffs(HURLEY_MCO::HEF);

    if (p_Mass > mHeF) {
        if (p_Mass >= 2.5) {
            p = 5.0;
        }
        else {                                                  // linear interpolation between end points
            const double gradient  = 1.0 / (mHeF - 2.5);        // will be negative
            const double intercept = 5.0 - (2.5 * gradient);
            p = (gradient * p_Mass) + intercept;
        }
    }

    return p;
}


/*
 * CalculateCoreMass_Luminosity_q_Hurley2000
 *
 * @brief
 * Calculate the Core mass - Luminosity relation parameter, q, per Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_q_Hurley2000(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Core mass - Luminosity relation parameter q
 */
inline double GiantBranch::CalculateCoreMass_Luminosity_q_Hurley2000(const double p_Mass) const {

    double q = 3.0;     // default value

    const double mHeF = GLOBALS->HurleyMassCutoffs(HURLEY_MCO::HEF);

    if (p_Mass > mHeF) {
        if (p_Mass >= 2.5) {
            q = 2.0;
        }
        else {          // linear interpolation between end points
            double gradient  = 1.0 / (mHeF - 2.5);
            double intercept = 2.0 - (2.5 * gradient);
            q = (gradient * p_Mass) + intercept;
        }
    }

    return q;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRadius
 *
 * @brief
 * Calculate the radius of the star.
 *
 * Calls relevant function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateRadius(const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Msol)
 * @return                                      Radius (Rsol)
 */
inline double GiantBranch::CalculateRadius(const double p_Mass, const double p_Luminosity) const { 

    double radius;

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            radius = CalculateRadius_Hurley2000(const double p_Mass, const double p_Luminosity);
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't accounted
            // for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return radius;
}


/*
 * CalculateRadiusAtHeI_Hurley2000
 *
 * @brief
 * Calculate radius at Helium Ignition, per Hurley et al. 2000, eq 50
 *
 *
 * double CalculateRadiusAtHeI_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MinLuminosity                 Minimum luminosity on phase (Lsol)
 * @return                                      Radius at Helium Ignition (Rsol)
 */
inline double GiantBranch::CalculateRadiusAtHeI_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MinLuminosity) const {

    double radius;

    const double mFGB = GLOBALS->HurleyMassCutoffs(HURLEY_MCO::FGB);
    const double lHeI = CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass);

    if (p_Mass <= mFGB) {
        radius = CalculateRadius_Hurley2000(p_Mass, lHeI);
    }
    else {
        const double rmHe = CHeB::CalculateMinRadius_Hurley2000_Static(p_Mass, p_CoreMass, p_MinLuminosity);
    
        if (p_Mass >= std::max(mFGB, HIGH_MASS_THRESHOLD)) {
            radius = std::min(rmHe, EAGB::CalculateRadius_Hurley2000_Static(p_Mass, lHeI));
        }
        else {
            radius = rmHe * PPOW(CalculateRadius_Hurley2000(p_Mass, lHeI) / rmHe, std::log10(p_Mass / HIGH_MASS_THRESHOLD) / std::log10(mFGB / HIGH_MASS_THRESHOLD));
        }
    }

    return radius;
}


/*
 * CalculateRadiusOnZAHB_Hurley2000
 *
 * @brief
 * Calculate radius on the Zero Age Horizontal Branch, per Hurley et al. 2000, eq 54
 *
 *
 * double CalculateRadiusOnZAHB_Hurley2000(onst double p_Mass, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Radius on the Zero Age Horizontal Branch (Rsol)
 */
inline double GiantBranch::CalculateRadiusOnZAHB_Hurley2000(const double p_Mass, const double p_CoreMass) const {

    const DblVectorT b = GLOBALS->HurleyBcoefficients(); // get Hurley b coefficients
    const double mHeF  = GLOBALS->HurleyMassCutoffs(HURLEY_MCO::HEF);
    const double rGB   = CalculateRadius_Hurley2000(p_Mass, CalculateLuminosityOnZAHB_Hurley2000(p_Mass, p_CoreMass));
    const double f     = ((1.0 + b[21]) * PPOW((p_Mass - p_CoreMass) / (mHeF - p_CoreMass), b[22])) / (1.0 + b[21] * PPOW(mu, b[23]));

    return ((1.0 - f)) * HeMS::CalculateRadiusAtZAHeMS_Hurley2000_Static(p_CoreMass) + (f * rGB);
}


/*
 * CalculateRadius_Hurley2000
 *
 * @brief
 * Calculate radius on the Giant Branch, per Hurley et al. 2000, eq 46
 * (see also just before eq 47)
 *
 *
 * double CalculateRadius_Hurley2000(const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      GB radius (Rsol)
 */
inline double GiantBranch::CalculateRadius_Hurley2000(const double p_Mass, const double p_Luminosity) const {
    const DblVectorT b = GLOBALS->HurleyBcoefficients(); // get Hurley b coefficients
    return std::min((b[4] * PPOW(p_Mass, -b[5])), (b[6] * PPOW(p_Mass, -b[7]))) * (PPOW(p_Luminosity, b[1]) + (b[2] * PPOW(p_Luminosity, b[3])));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     REMNANTS                                      //
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
 * double CalculateRemnantRadius_Hurley2000(const double p_Mass, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol) (typically effective mass, mass0)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Remnant core radius (Rsol)
 */
inline double GiantBranch::CalculateRemnantRadius_Hurley2000(const double p_Mass, const double p_CoreMass) const {
    return p_Mass > GLOBALS->HurleyMassCutoffs(HURLEY_MCO::HEF)
            ? HeMS::CalculateRadiusAtZAHeMS_Hurley2000_Static(p_CoreMass)
            : WhiteDwarfs::CalculateRadius_Marsh2004_Static(p_CoreMass);
}


/*
 * CalculateRemnantType_Muller2016
 *
 * @brief
 * Calculate remnant type given, per Muller et al. 2016.
 * 
 * As presented in Vigna-Gomez et al. 2018, appendix B (See arXiv:1805.07974)
 *
 *
 * STELLAR_TYPE CalculateRemnantType_Muller2016(const double p_COCoreMass) const
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @return                                      Remnant stellar type
 */
inline STELLAR_TYPE GiantBranch::CalculateRemnantType_Muller2016(const double p_COCoreMass) const {

    STELLAR_TYPE stellarType;

         if (p_COCoreMass < 3.6 ) stellarType = STELLAR_TYPE::NEUTRON_STAR;
    else if (p_COCoreMass < 4.05) stellarType = STELLAR_TYPE::BLACK_HOLE;
    else if (p_COCoreMass < 4.6 ) stellarType = STELLAR_TYPE::NEUTRON_STAR;
    else if (p_COCoreMass < 5.7 ) stellarType = STELLAR_TYPE::BLACK_HOLE;
    else if (p_COCoreMass < 6.0 ) stellarType = STELLAR_TYPE::NEUTRON_STAR;
    else                          stellarType = STELLAR_TYPE::BLACK_HOLE;

    return stellarType;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ROTATION                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateMomentOfInertia_Hurley2000
 *
 * Calculate moment of inertia, per Hurley et al., 2000,
 * paragraph immediately following eq 109 
 *
 * 
 * double CalculateMomentOfInertia()
 * 
 * @return                                      Moment of inertia (Msol AU^2)
 */
inline double GiantBranch::CalculateMomentOfInertia(const double p_Mass, const double p_Radius, const double p_CoreMass) const {
    const double Rc = CalculateRemnantRadius_Hurley2000();
    return (0.1 * (p_Mass - p_CoreMass) * p_Radius * p_Radius) + (0.21 * p_CoreMass * Rc * Rc);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLuminosityAtBGB_Hurley2000_Static
 *
 * @brief
 * Calculate luminosity at the base of the giant branch (BGB),
 * per Hurley et al. 2000, eq 10
 *
 *
 * static double CalculateLuminosityAtBGB_Hurley2000_Static(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      BGB luminosity (Lsol)
 */
inline double GiantBranch::CalculateLuminosityAtBGB_Hurley2000_Static(const double p_Mass) {

    const DBL_VECTOR a = GLOBALS->HurleyAcoefficients(); // get Hurley a coefficients

    return (a[27] * PPOW(p_Mass, a[31])) + (a[28] * PPOW(p_Mass, HURLEY_C_COEFF[2])) / a[29] + (a[30] * PPOW(p_Mass, HURLEY_C_COEFF[3])) + PPOW(p_Mass, a[32]);
}


/*
 * CalculateLuminosityAtHeI_Hurley2000_Static
 *
 * @brief
 * Calculate luminosity at Helium Ignition, HeI, per Hurley et al. 2000, eq 49
 *
 *
 * static double CalculateLuminosityAtHeI_Hurley2000_Static(const double p_Mass)
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      HeI luminosity (Lsol)
 */
static inline double GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(const double p_Mass) {

    const DblVectorT b = GLOBALS->HurleyBcoefficients();
    const double MHeF  = GLOBALS->HurleyMassCutoffs(HURLEY_MCO::HEF);
    
    return p_Mass < MHeF
            ? (b[9] * PPOW(p_Mass, b[10])) / (1.0 + (GLOBALS->HurleyAlpha1() * std::exp(15.0 * (p_Mass - MHeF))))
            : (b[11] + (b[12] * PPOW(p_Mass, 3.8))) / (b[13] + (p_Mass * p_Mass));
}


/*
 * CalculateLuminosityOnZAHB_Hurley2000
 *
 * @brief
 * Calculate luminosity on the Zero Age Horizontal Branch, for Low Mass stars,
 * per Hurley et al. 2000, eq 53
 *
 * **Ilya**
 * JR: is there a check for LM stars? Should we enforce LM, or remove the caveat/qualifier?
 * Hurley defines LM, IM, and HM at the end of Hurley et al. 2000, sec 5.  JD introduced the
 * HIGH_MASS_THRESHOLD constant, but I don't think we have a LM equivalent (or IM for that matter).
 * 
 *
 * double CalculateLuminosityOnZAHB_Hurley2000(const double p_Mass, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      ZAHB luminosity (Lsol)
 */
inline double GiantBranch::CalculateLuminosityOnZAHB_Hurley2000(const double p_Mass, const double p_CoreMass) const {

    const DblVectorT b = GLOBALS->HurleyBcoefficients(); // get Hurley b coefficients
    const double MHeF  = GLOBALS->HurleyMassCutoffs(HURLEY::MCO::HEF);
    const double lCHeB = GLOBALS->MinLuminosity_CHeB();
    const double lZHe  = HeMS::CalculateLuminosityAtZAHeMS_Hurley2000_Static(p_CoreMass);

    const double mu    = (p_Mass - p_CoreMass) / (MHeF - p_CoreMass);
    const double tmp   = (b[18] * PPOW(mu, b[19])) / (1.0 + (b[18] + lZHe - lCHeB) / (lCHeB - lZHe) * std::exp(15.0 * (p_Mass - MHeF)));

    return lZHe + ((1.0 + b[20]) / (1.0 + (b[20] * PPOW(mu, 1.6479))) * tmp);
}


/*
 * CalculateRemnantLuminosity_Hurley2000
 *
 * @brief
 * Calculate the luminosity of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. mass = coreMass), per Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantLuminosity_Hurley2000(const double p_Mass, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol) (typically effective mass, mass0)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Luminosity of remnant core (Lsol)
 */
inline double GiantBranch::CalculateRemnantLuminosity_Hurley2000(const double p_Mass, const double p_CoreMass) const {
    return p_Mass > GLOBALS->HurleyMassCutoffs(HURLEY_MCO::HEF)
            ? HeMS::CalculateLuminosityAtZAHeMS_Hurley2000_Static(p_CoreMass)
            : WhiteDwarfs::CalculateLuminosity_Hurley2000_Static(p_CoreMass, 0.0);
}
























/// GiantBranch_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


class GiantBranch_Constituent: virtual public BinaryConstituentStar, public GiantBranch {

public:


protected:

};



COMPAS_PURE double CalculateCriticalMassRatio_Claeys2014(const double p_Mass, const double p_HeCoreMass, const bool p_AccretorIsDegenerate) const; 
double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const {
    return CalculateCriticalMassRatio_Claeys2014(m_StateHistory.CurrentState().Mass(), m_StateHistory.CurrentState().HeCoreMass(), p_AccretorIsDegenerate);
}



GNU_CONST double CalculateCriticalMassRatio_Hurley2002(const double p_Mass, const double p_CoreMass) const; 
double CalculateCriticalMassRatio_Hurley2002() const {
    return CalculateCriticalMassRatio_Hurley2002(m_StateHistory.CurrentState().Mass(), m_StateHistory.CurrentState().CoreMass());
}; 



/*
 * CalculateCriticalMassRatio_Hurley2002
 *
 * @brief
 * Calculate the critical mass ratio for unstable mass transfer, per Hurley et al. 2002 sect. 2.6.1.
 *
 * See Hurley et al. 2002 sect. 2.6.1, and Hjellming & Webbink 1987.
 *
 * Assumes this star is the donor.
 * 
 * Critical mass ratio is defined as qCrit = Maccretor / Mdonor.
 * Note: mass ratio is defined as Mdonor / Maccretor in Hurley et al. 2002, so need to invert here.
 *
 *
 * double CalculateCriticalMassRatio_Hurley2002(const double p_Mass, const double p_CoreMass) const 
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Critical mass ratio for unstable MT 
 */
inline GNU_CONST double GiantBranch_Constituent::CalculateCriticalMassRatio_Hurley2002(const double p_Mass, const double p_CoreMass) const {
    return 1.0 / (0.362 + 1.0 / (3.0 * (1.0 - p_CoreMass / p_Mass)));   // Hurley et al. 2002, just after eq. 57
}










#endif // __GiantBranch_h__
