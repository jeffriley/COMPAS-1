#pragma once

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "Rand.h"
#include "MainSequence.h"
#include "ConstituentStar.h"   // for GiantBranch_Constituent (declared at the bottom of this file)


class BaseStar;
class MainSequence;

class GiantBranch: virtual public BaseStar, public MainSequence {



public:
  
    GiantBranch(const BaseStar &p_BaseStar) : BaseStar(p_BaseStar), MainSequence(p_BaseStar) {}  


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

    ////////////////////////////////////////
    //   ABUNDANCE                        //
    ////////////////////////////////////////











    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateAgeAfterMassLoss_Hurley2000() const override { return Age(); } // NO-OP // JR FIX THIS: DONE

    double CalculateLifetimeToHeI_Hurley2000(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB, const DblVectorT& p_GBParams) const;






    ////////////////////////////////////////
    //   ENVELOPE                         //
    ////////////////////////////////////////

    


   







    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////







        
    ////////////////////////////////////////
    //   MASS                             //
    ////////////////////////////////////////


    virtual double CalculateInitialSupernovaMass() const { return Mass(); }

    






    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////

            
            
            double      CalculateFallbackFraction_Fryer2012_Delayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const;
            double      CalculateFallbackFraction_Fryer2012_Rapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const;
            double      CalculateGravitationalRemnantMass_Fryer2012(const double p_BaryonicRemnantMass, const double p_NSmaxBaryonicMass) const;
            Dbl_DblT    CalculateRemnantMass_Schneider2020(const double p_Mass, const double p_COCoreMass, const double p_HeCoreMassPreSN, const double p_MaxNSMass, const REMNANT_MASS_PRESCRIPTION p_RemnantMassPrescription) const;
            Dbl_DblT    CalculateRemnantMass_MullerMandel2020(const double p_COCoreMass, const double p_HeCoreMass);


            double          CalculateFallbackFractionRapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const {  } // JR FIX THIS: placeholder



            
            
            double          CalculateGravitationalRemnantMass(const double p_BaryonicRemnantMass);




            double          CalculateRemnantMassByMullerMandel(const double p_COCoreMass, const double p_HeCoreMass);
            
            double          CalculateRemnantMassBySchneider2020(const double p_COCoreMass, const bool p_UseSchneiderAlt = false);
            double          CalculateRemnantMassBySchneider2020Alt(const double p_COCoreMass)               { return CalculateRemnantMassBySchneider2020(p_COCoreMass, true); }
            double          CalculateRemnantNSMassMullerMandel(const double p_COCoreMass, const double p_HeCoreMass);





























	





            double          CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription);
            double          CalculateZetaConvectiveEnvelopeGiant(ZETA_PRESCRIPTION p_ZetaPrescription);

    virtual void            PerturbLuminosityAndRadius();

    
            void            UpdateAgeAfterMassLoss() { }                                                                                                                        // NO-OP for most stellar types


            void            UpdateMainSequenceCoreMass(const double p_Dt, const double p_TotalMassLossRate) { }                                                                 // NO-OP for most stellar types




///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


void CalculateTimescales_Hurley2000() override { m_InterimState.SetHurleyTimescales(CalculateTimescales_Hurley2000(MassEffectiveInitial(), m_InterimState.HurleyGBParamsOrDefault(), m_InterimState.HurleyTimescalesOrDefault())); }
COMPAS_PURE DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const;










 







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
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////







    //////////////////////////////////////////////////
    //   ENVELOPE                                   //
    //////////////////////////////////////////////////

    COMPAS_PURE Dbl_DblT CalculateConvectiveEnvelopeMass_Picker2024(const double p_Metallicity, const double p_Mass, const double p_Mass0, const double p_Temperature) const; // JR FIX THIS: DONE
    COMPAS_PURE Dbl_DblT CalculateConvectiveEnvelopeMass() const override { // JR FIX THIS: DONE
        return CalculateConvectiveEnvelopeMass_Picker2024(Metallicity(), Mass(), Mass0(), Temperature());
    }



    double          CalculateConvectiveEnvelopeRadialExtent() const { return 0.0; }
    COMPAS_PURE double CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Radius) const;




    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

       
    GNU_PURE double CalculateLuminosityAtBGB_Hurley2000(const double p_Metallicity, const double p_Mass) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateLuminosityAtHeI_Hurley2000(const double p_Metallicity, const double p_Mass) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateLuminosityOnZAHB_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const; // JR FIX THIS: DONE

    virtual double CalculateRemnantLuminosity_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const;




    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////

    GNU_CONST   double   CalculateBaryonicRemnantMass_Fryer2012(const double p_ProtoMass, double p_FallbackMass) const;

                double   CalculateBHMassAfterFallback_MullerMandel2020(const double p_COCoreMass, const double p_HeCoreMass) const;

    GNU_PURE double CalculateConvectiveCoreMass() const override { return CoreMass(); }


    GNU_CONST double CalculateCoreMassAt2ndDredgeUp_Hurley2000(const double p_McBAGB) const; // JR FIX THIS: DONE
    GNU_PURE  double CalculateCoreMassAtBAGB_Hurley2000(const double p_Metallicity, const double p_Mass) const; // JR FIX THIS: DONE
    GNU_PURE  double CalculateCoreMassAtBGB_Hurley2000(const double p_Metallicity, const double p_Mass, const GBParamsT& p_GBParams) const; // JR FIX THIS: DONE
    GNU_PURE  double CalculateCoreMassAtHeI_Hurley2000(const double p_Metallicity, const double p_Mass, const GBParamsT& p_GBParams) const; // JR FIX THIS: DONE
    GNU_CONST double CalculateCoreMassAtSN(const double p_mThreshold, const double p_McBAGB) const; // JR FIX THIS: DONE


    double   CalculateEffectiveInitialMass_Hurley2000() const override { return BaseStar::CalculateEffectiveInitialMass_Hurley2000(); } // per Hurley et al. 2000, section 7.1

    GNU_CONST   double   CalculateFallbackMass_Fryer2012(const double p_PreSNMass, const double p_ProtoMass, const double p_FallbackFraction) const;

    GNU_CONST   double   CalculateProtoCoreMass_Fryer2012_Delayed(const double p_COCoreMass) const;

    GNU_CONST   Dbl_DblT CalculateRemnantMass_Belczynski2002(const double p_Mass, const double p_COCoreMass) const;
    COMPAS_PURE Dbl_DblT CalculateRemnantMass_Fryer2012(const double p_Mass, const double p_COCoreMassn) const;
    COMPAS_PURE Dbl_DblT CalculateRemnantMass_Fryer2022(const double p_Mass, const double p_COCoreMass, const double p_fMix, const double p_mCrit) const;
                double   CalculateRemnantMass_Maltsev2025(const double p_COCoreMass, const double p_HeCoreMass) const;
    GNU_CONST   double   CalculateRemnantMass_Muller2016(const double p_Mass, const double p_COCoreMass) const;


    //////////////////////////////////////////////////
    //   MASS LOSS / ACCRETION                      //
    //////////////////////////////////////////////////

    double CalculateMLRateThermal() const { return (Mass() - CoreMass()) / CalculateTimescale_Thermal(); }


    //////////////////////////////////////////////////
    //   MISC. CONSTANTS / PARAMETERS               //
    //////////////////////////////////////////////////

    virtual double CalculateCoreMass_Luminosity_B_Hurley2000(const double p_Mass) const;
    virtual double CalculateCoreMass_Luminosity_D_Hurley2000(const double p_Mass) const;
    virtual double CalculateCoreMass_Luminosity_p_Hurley2000(const double p_Mass) const;
    virtual double CalculateCoreMass_Luminosity_q_Hurley2000(const double p_Mass) const;


    DblVectorT CalculateGBParams_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBParams) const override;

    GNU_CONST double CalculateCoreMass_Luminosity_Lx_Hurley2000(const DblVectorT& p_GBParams) const;
    GNU_CONST double CalculateCoreMass_Luminosity_Mx_Hurley2000(const DblVectorT& p_GBParams) const;

    GNU_CONST double CalculateHRateConstant_Hurley2000(const double p_Mass) const;
    GNU_CONST double CalculateSmallEnvPerturbationParm_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_Luminosity = 0.0) const override;


    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////

    double CalculateConvectiveCoreRadius() const override { // JR FIX THIS: DONE
        return std::min(CalculateRemnantRadius(), Radius()); // Last paragraph of section 6 of Hurley et al. 2000
    }

    COMPAS_PURE double CalculateRadiusAtHeI_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const; // JR FIX THIS: DONE
    }

    double CalculateRadiusOnMassChange(double p_dM) const override { return CalculateRadius(Mass() + p_dM, Luminosity()); } // JR FIX THIS: DONE

    COMPAS_PURE double CalculateRadiusOnZAHB_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const; // JR FIX THIS: DONE

    COMPAS_PURE double CalculateRadius_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Luminosity) const; // JR FIX THIS: DONE
    double CalculateRadius_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateRadius_Hurley2000(Metallicity(), Mass(), Luminosity());
    }

    COMPAS_PURE double CalculateRemnantRadius_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const; // JR FIX THIS: DONE
    COMPAS_PURE double CalculateRemnantRadius_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateRemnantRadius_Hurley2000(Metallicity(), Mass0(), CoreMass());
    }


    //////////////////////////////////////////////////
    //   REMNANTS                                   //
    //////////////////////////////////////////////////


    GNU_CONST STELLAR_TYPE CalculateRemnantType_Muller2016(const double p_COCoreMass) const;


    //////////////////////////////////////////////////
    //   ROTATION                                   //
    //////////////////////////////////////////////////

    COMPAS_PURE double CalculateMomentOfInertia(const double p_Mass, const double p_Radius, const double p_CoreMass) const;

















///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<






};




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
 * CalculateLifetimeToHeI_Hurley2000
 *
 * @brief
 * Calculate the time to Helium ignition, tHeI, per Hurley et al. 2000, eq 43
 *
 *
 * double CalculateLifetimeToHeI_Hurley2000(const double p_Metallicity, const double p_Mass, const GBParamsT& p_GBParams, const TimescalesT& p_tScales) const
 *
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBParams                      GB parameters
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Lifetime to He ignition (tHeI)
 */
inline double GiantBranch::CalculateLifetimeToHeI_Hurley2000(const double p_Metallicity, const double p_Mass, const GBParamsT& p_GBParams, const TimescalesT& p_tScales) const {

    const double AH = p_GBParams[HURLEY_GBP::AH];
    const double B  = p_GBParams[HURLEY_GBP::B];
    const double D  = p_GBParams[HURLEY_GBP::D];
    const double p  = p_GBParams[HURLEY_GBP::P];
    const double q  = p_GBParams[HURLEY_GBP::Q];
    const double Lx = p_GBParams[HURLEY_GBP::LX];
    const double p1 = p - 1.0;
    const double q1 = q - 1.0;

    const double mHeF    = ZDEP->HurleyMCOHeF(p_Metallicity);
    const double alpha1  = ZDEP->HurleyAlpha1(p_Metallicity);
    const auto&  bCoeffs = ZDEP->HurleyBCoeffs(p_Metallicity); // Hurley b coefficients

    const double lHeI = CalculateLuminosityAtHeI_Hurley2000(p_Mass, mHeF, alpha1, bCoeffs);

    return (lHeI > Lx)
            ? p_tScales[TS::INF2_FGB] - (1.0 / (q1 * AH * B)) * PPOW((B / lHeI), (q1 / q))
            : p_tScales[TS::INF1_FGB] - (1.0 / (p1 * AH * D)) * PPOW((D / lHeI), (p1 / p));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ENVELOPE                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

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
 * Dbl_DblT CalculateConvectiveEnvelopeMass_Picker2024(const double p_Metallicity, const double p_Mass, const double p_Mass0, const double p_Temperature) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Mass0                         Effective initial mass of the star (Msol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                   GB convective envelope mass (Msol)
 *                                                   Maximum GB convective envelope mass (Msol)
 */
inline Dbl_DblT GiantBranch::CalculateConvectiveEnvelopeMass_Picker2024(const double p_Metallicity, const double p_Mass, const double p_Mass0, const double p_Temperature) const {
    
    const double sigmaHurley = ZDEP->SigmaHurley(p_Metallicity);

    // We need the temperature of the star just after BAGB, which is the temperature at the
    // start of the EAGB phase.  Since we are on the giant branch here, we can clone this
    // object as an EAGB object and, as long as it is initialised (to the start of the phase),
    // we can query the cloned object for its temperature.
    //
    // To ensure the clone does not participate in logging, we set its persistence to EPHEMERAL.
    std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH, OBJECT_PERSISTENCE::EPHEMERAL);
    const double tMin = clone->Temperature();
    
    // Use Eq. 6 of Mandel et al. 2024 eq 6 rather than Picker at al. 2024 eq 6 for tOnset
    // to avoid issues caused by differences between temperatures in MESA models (used in
    // fits from Picker et al. 2024) and Pols models (used in Hurley SSE tracks).
    const double tOnset = tMin / std::min(0.695 - 0.057 * sigmaHurley, 0.95);                                       // Mandel et al. 2024, eq 6

    const double McFinal = CalculateCoreMassAtBAGB_Hurley2000(p_Metallicity, p_Mass0);

    // Ratio of final intershell mass to final core mass 
    const double MiFinalMcFinal = -0.023 * sigmaHurley - 0.0023;                                                    // Picker at al. 2024 eq 8.

    // Picker+ 2024 fits were only made for stars above 8.0 solar masses, with runs down to 5.0 Msol,
    // so we use the final core mass as an approximate threshold of validity.
    // Unlike massive stars, intermediate-mass stars have almost no radiative intershell at maximum
    // convective envelope extent.
    double mEnvMax = (McFinal < 1.5) 
                        ? std::max(p_Mass - McFinal, 0.0)
                        : std::max(p_Mass - McFinal * (1.0 + MiFinalMcFinal), 0.0);                                 // Picker at al. 2024, eq 9

    const double mEnv = mEnvMax / (1.0 + std::exp(4.6 * (tMin + tOnset - 2.0 * p_Temperature) / (tMin - tOnset)));  // Picker at al. 2024, eq 7
    
    return std::make_tuple(mEnv, mEnvMax);
}


/*
 * CalculateConvectiveEnvelopeRadialExtent_Hurley2002
 *
 * @brief
 * Calculate the radial extent of the convective outer envelope,
 * per Hurley et al. 2002
 *
 * Implementation is a combination of:
 * 
 *    Hurley et al. 2000, end of sec 7.2, and
 *    Hurley et al. 2002, sec 2.3, particularly subsec 2.3.1, eqs 39-40
 *
 *
 * double CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Radius) const
 *
 * @param       p_Radius                        Radius of the star (Rsol)
 * @return                                      Radial extent of the convective outer envelope (Rsol)
 */
inline double GiantBranch::CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Radius) const {

    double Menv;
    double MZAMSenv;
    std::tie(Menv, MZAMSenv) = CalculateConvectiveEnvelopeMass();
    
    return (!(Menv > 0.0) || !(MZAMSenv > 0.0)) ? 0.0 : std::sqrt(Menv / MZAMSenv) * (p_Radius - CalculateConvectiveCoreRadius());
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
 * @return                                      Core mass of the star at second dredge up (Msol)
 */
inline double GiantBranch::CalculateCoreMassAt2ndDredgeUp_Hurley2000(const double p_McBAGB) const {
    return (p_McBAGB >= 0.8) ? ((0.44 * p_McBAGB) + 0.448) : p_McBAGB;
}


/*
 * CalculateCoreMassAtBAGB_Hurley2000
 *
 * @brief
 * Calculate core mass at the Base of the Asymptotic Giant Branch, per Hurley et al. 2000, eq 66
 *
 *
 * double CalculateCoreMassAtBAGB_Hurley2000(const double p_Metallicity, const double p_Mass) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Core mass of the star at BAGB (Msol)
 */
inline double GiantBranch::CalculateCoreMassAtBAGB_Hurley2000(const double p_Metallicity, const double p_Mass) const {
    const auto& b = ZDEP->HurleyBCoefficients(p_Metallicity);           // Hurley b coefficients
    return std::sqrt(std::sqrt((b[36] * PPOW(p_Mass, b[37])) + b[38])); // sqrt() is much faster than pow()
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
 * double CalculateCoreMassAtBGB_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_McBAGB, const GBParamsT& p_GBParams) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBParams                      GB parameters
 * @return                                      Core mass of the star at BGB (Msol)
 */
inline double GiantBranch::CalculateCoreMassAtBGB_Hurley2000(const double p_Metallicity, const double p_Mass, const GBParamsT& p_GBParams) const {

    double McBGB = 0.0;                                                                                 // Default return value

    const double mHeF = ZDEP->HurleyMCOHeF(p_Metallicity);

    // No McBGB for stars with mass below the helium flash threshold.
    // See Hurley at al. 2000, text immediately prior to eq 44.
    if (!(mHeF < p_Mass)) {
        double luminosity = CalculateLuminosityAtBGB_Hurley2000(p_Metallicity, mHeF);
        double Mc_MHeF    = CalculateCoreMass_Hurley2000(luminosity, p_GBParams);
        double c          = (Mc_MHeF * Mc_MHeF * Mc_MHeF * Mc_MHeF) - (MC_L_C1 * PPOW(mHeF, MC_L_C2));  // pow() is slow - use multiplication   
        
        McBGB = std::min((0.95 * p_GBParams[HURLEY_GBP::McBAGB]), std::sqrt(std::sqrt(c + (MC_L_C1 * PPOW(p_Mass, MC_L_C2))))); // sqrt is much faster than pow()
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
 * double CalculateCoreMassAtHeI_Hurley2000(const double p_Metallicity, const double p_Mass, const GBParamsT& p_GBParams) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBParams                      GB parameters
 * @return                                      Core mass of the star at HeI (Msol)
 */
inline double GiantBranch::CalculateCoreMassAtHeI_Hurley2000(const double p_Metallicity, const double p_Mass, const GBParamsT& p_GBParams) const {

    double coreMass = p_Mass;                                                                               // Default return value

    const double mHeF    = ZDEP->HurleyMCOHeF(p_Metallicity);
    const double alpha1  = ZDEP->HurleyAlpha1(p_Metallicity);
    const auto&  bCoeffs = ZDEP->HurleyBCoeffs(p_Metallicity);                                              // Hurley b coefficients

    const double lHeI = CalculateLuminosityAtHeI_Hurley2000(p_Mass, mHeF, alpha1, bCoeffs);

    if (p_Mass < mHeF) {
        coreMass = CalculateCoreMass_Hurley2000(lHeI, p_GBParams);
    }
    else {
        const double Mc_MHeF = CalculateCoreMass_Hurley2000(lHeI, p_GBParams);
        const double McBAGB  = CalculateCoreMassAtBAGB_Hurley2000(p_Metallicity, p_Mass);
        const double c       = (Mc_MHeF * Mc_MHeF * Mc_MHeF * Mc_MHeF) - (MC_L_C1 * PPOW(mHeF, MC_L_C2));   // pow() is slow - use multiplication

        coreMass = std::min((0.95 * McBAGB), std::sqrt(std::sqrt(c + (MC_L_C1 * PPOW(p_Mass, MC_L_C2)))));  // sqrt() is much faster than pow()
    }

    return coreMass;
}


/*
 * CalculateCoreMassAtSN
 *
 * @brief
 * Calculate the core mass at which the Asymptotic Giant Branch phase
 * is terminated in a SN/loss of envelope, per Hurley et al. 2000, eq 75.
 * 
 * The mass returned is clamped to a minimum CO core mass (passed as p_mThreshold), which
 * is MCH in Hurley et al. eq 75, but we deviate from Hurley for stars that may explode in
 * ECSNe and pass MECS (instead of MCH) in those cases.
 *
 *
 * static double CalculateCoreMassAtSN(const double p_mThreshold, const double p_McBAGB)
 *
 * @param       p_mThreshold                    Threshold mass (see notes above) (Msol)
 * @param       p_McBAGB                        Core mass at the Base of the Asymptotic Giant Branch (Msol)
 * @return                                      Maximum core mass pre-SN on the Asymptotic Giant Branch (Msol)
 */
inline double GiantBranch::CalculateCoreMassAtSN(const double p_mThreshold, const double p_McBAGB) const {
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
 * CalculateSmallEnvPerturbationParm_Hurley2000
 *
 * @brief
 * Calculate the small envelope perturbation parameter, per Hurley et al. 2000 eqs 97 & 98
 *
 *
 * double CalculateSmallEnvPerturbationParm_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Small envelope perturbation parameter
 */
inline double GiantBranch::CalculateSmallEnvPerturbationParm_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_Luminosity) const {
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
 * double CalculateCoreMass_Luminosity_Lx_Hurley2000(const DblVectorT& p_GBParams) const
 *
 * @param       p_GBParams                      GB parameters
 * @return                                      Core mass - Luminosity relation parameter, Lx
 */
inline double GiantBranch::CalculateCoreMass_Luminosity_Lx_Hurley2000(const DblVectorT& p_GBParams) const {

    const double B  = p_GBParams[HURLEY_GBP::B];
    const double D  = p_GBParams[HURLEY_GBP::D];
    const double Mx = p_GBParams[HURLEY_GBP::MX];
    const double p  = p_GBParams[HURLEY_GBP::P];
    const double q  = p_GBParams[HURLEY_GBP::Q];

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
 * double CalculateCoreMass_Luminosity_Mx_Hurley2000(const DblVectorT& p_GBParams) const
 *
 * @param       p_GBParams                      GB parameters
 * @return                                      Core mass - Luminosity relation parameter, Mx
 */
inline double GiantBranch::CalculateCoreMass_Luminosity_Mx_Hurley2000(const DblVectorT &p_GBParams) const {

    const double B  = p_GBParams[HURLEY_GBP::B];
    const double D  = p_GBParams[HURLEY_GBP::D];
    const double p  = p_GBParams[HURLEY_GBP::P];
    const double q  = p_GBParams[HURLEY_GBP::Q];

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

    const double D0 = 5.37 + (0.135 * ZDEP->ZetaHurley(Metallicity()));
    const double D1 = (0.975 * D0) - (0.18 * p_Mass);

    double logD = D0;   // default value

    const double mHeF = ZDEP->HurleyMassCutoffs(Metallicity(), HURLEY_MCO::HeF);

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

    const double mHeF = ZDEP->HurleyMassCutoffs(Metallicity(), HURLEY_MCO::HeF);

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

    const double mHeF = ZDEP->HurleyMassCutoffs(Metallicity(), HURLEY_MCO::HeF);

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
 * double CalculateRadius(const double p_Metallicity, const double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Msol)
 * @return                                      Radius (Rsol)
 */
inline double GiantBranch::CalculateRadius(const double p_Metallicity, const double p_Mass, const double p_Luminosity) const { 

    double radius;

    switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            radius = CalculateRadius_Hurley2000(p_Metallicity, p_Mass, p_Luminosity);
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
 * double CalculateRadiusAtHeI_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Radius at Helium Ignition (Rsol)
 */
inline double GiantBranch::CalculateRadiusAtHeI_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const {

    double radius = 0.0;                                        // Default return value

    const double mHeF    = ZDEP->HurleyMCOHeF(p_Metallicity);
    const double alpha1  = ZDEP->HurleyAlpha1(p_Metallicity);
    const auto&  bCoeffs = ZDEP->HurleyBCoeffs(p_Metallicity);  // Hurley b coefficients

    const double lHeI = CalculateLuminosityAtHeI_Hurley2000(p_Mass, mHeF, alpha1, bCoeffs);

    const double mFGB = ZDEP->HurleyMCOFGB(p_Metallicity);

    if (mFGB > p_Mass) {
        radius = CalculateRadius_Hurley2000(p_Metallicity, p_Mass, lHeI);
    }
    else {
        const double rmHe = CHeB::CalculateMinRadius_Hurley2000(p_Metallicity, p_Mass, p_CoreMass);
        if (p_Mass < std::max(mFGB, HIGH_MASS_THRESHOLD)) {
            radius = rmHe * PPOW(CalculateRadius_Hurley2000(p_Metallicity, p_Mass, lHeI) / rmHe, std::log10(p_Mass / HIGH_MASS_THRESHOLD) / std::log10(mFGB / HIGH_MASS_THRESHOLD));
        }
        else {
            radius = std::min(rmHe, EAGB::CalculateRadius_Hurley2000(p_Metallicity, p_Mass, lHeI));
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
 * double CalculateRadiusOnZAHB_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Radius on the Zero Age Horizontal Branch (Rsol)
 */
double GiantBranch::CalculateRadiusOnZAHB_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const {

    const double mHeF = ZDEP->HurleyMCOHeF(p_Metallicity);
    const auto&  b    = ZDEP->HurleyBCoeffs(p_Metallicity); // Hurley b coefficients

    const double rZAHeMS = HeMS::CalculateRadiusAtZAHeMS_Hurley2000(p_CoreMass);
    const double lZAHB   = CalculateLuminosityOnZAHB(p_Metallicity, p_Mass, p_CoreMass);
    const double rGB     = CalculateRadius_Hurley2000(p_Metallicity, p_Mass, lZAHB);

    const double mu    = (p_Mass - p_CoreMass) / (mHeF - p_CoreMass);
    const double f     = ((1.0 + b[21]) * PPOW(mu, b[22])) / (1.0 + b[21] * PPOW(mu, b[23]));

    return ((1.0 - f)) * rZAHeMS + (f * rGB);
}



/*
 * CalculateRadius_Hurley2000
 *
 * @brief
 * Calculate radius on the Giant Branch, per Hurley et al. 2000, eq 46
 * (see also just before eq 47)
 *
 *
 * double CalculateRadius_Hurley2000(const double p_Metallicity, onst double p_Mass, const double p_Luminosity) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      GB radius (Rsol)
 */
inline double GiantBranch::CalculateRadius_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Luminosity) const {
    const auto& b = ZDEP->HurleyBCoeffs(p_Metallicity); // Hurley b coefficients
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
 * immediately (i.e. M = Mc), per Hurley et al. 2000, just after eq 105.
 * 
 *
 * double CalculateRemnantRadius_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol) (typically effective mass, mass0)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Remnant core radius (Rsol)
 */
inline double GiantBranch::CalculateRemnantRadius_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const {
    return (p_Mass > ZDEP->HurleyMCOHeF(p_Metallicity))
            ? HeMS::CalculateRadiusAtZAHeMS_Hurley2000(p_CoreMass)
            : WhiteDwarfs::CalculateRadius_Marsh2004(p_CoreMass);
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
 * double CalculateMomentOfInertia(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_CoreMass)
 * 
 * @return                                      Moment of inertia (Msol AU^2)
 */
inline double GiantBranch::CalculateMomentOfInertia(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_CoreMass) const {
    const double Rc = CalculateRemnantRadius_Hurley2000(p_Metallicity, p_Mass, p_CoreMass);
    return (0.1 * (p_Mass - p_CoreMass) * p_Radius * p_Radius) + (0.21 * p_CoreMass * Rc * Rc);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLuminosityAtBGB_Hurley2000
 *
 * @brief
 * Calculate luminosity at the base of the giant branch (BGB), per Hurley et al. 2000, eq 10
 *
 *
 * static double CalculateLuminosityAtBGB_Hurley2000(const double p_Metallicity, const double p_Mass) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Luminosity of the etsra at BGB (Lsol)
 */
inline double GiantBranch::CalculateLuminosityAtBGB_Hurley2000(const double p_Metallicity, const double p_Mass) const {

    const auto&  a      = ZDEP->HurleyACoefficients(p_Metallicity); // Hurley a coefficients
    const double top    = (a[27] * PPOW(p_Mass, a[31])) + (a[28] * PPOW(p_Mass, C_COEFF.at(2)));
    const double bottom = a[29] + (a[30] * PPOW(p_Mass, C_COEFF.at(3))) + PPOW(p_Mass, a[32]);

    return top / bottom;
}


/*
 * CalculateLuminosityAtHeI_Hurley2000
 *
 * @brief
 * Calculate luminosity at helium ignition, HeI, per Hurley et al. 2000, eq 49
 *
 *
 * double CalculateLuminosityAtHeI_Hurley2000(const double p_Metallicity, const double p_Mass) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Luminosity of the star at HeI (Lsol)
 */
inline double GiantBranch::CalculateLuminosityAtHeI_Hurley2000(const double p_Metallicity, const double p_Mass) const {

    const auto&  b    = ZDEP->HurleyBCoefficients(p_Metallicity); // Hurley b coefficients
    const double mHeF = ZDEP->HurleyMCOHeF(p_Metallicity);

    return (p_Mass < mHeF)
            ? (b[9] * PPOW(p_Mass, b[10])) / (1.0 + (ZDEP->HurleyAlpha1(p_Metallicity) * std::exp(15.0 * (p_Mass - mHeF))))
            : (b[11] + (b[12] * PPOW(p_Mass, 3.8))) / (b[13] + (p_Mass * p_Mass));
}


/*
 * CalculateLuminosityOnZAHB_Hurley2000
 *
 * @brief
 * Calculate luminosity on the Zero Age Horizontal Branch, for low mass stars,
 * per Hurley et al. 2000, eq 53
 *
 * *Ilya*
 * Is there a check for LM stars? Should we enforce LM, or remove the caveat/qualifier?
 * Hurley defines LM, IM, and HM at the end of Hurley et al. 2000, sec 5.  JD introduced the
 * HIGH_MASS_THRESHOLD constant, but I don't think we have a LM equivalent (or IM for that matter).
 * 
 *
 * double CalculateLuminosityOnZAHB_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      ZAHB luminosity of the star (Lsol)
 */
double GiantBranch::CalculateLuminosityOnZAHB_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const {

    const double lMin = ZDEP->HurleyMinLuminosityOnCHeB(p_Metallicity);
    const double mHeF = ZDEP->HurleyMCOHeF(p_Metallicity);
    const auto&  b    = ZDEP->HurleyBCoefficients(p_Metallicity); // Hurley b coefficients

    const double lZHe = HeMS::CalculateLuminosityAtZAHeMS_Hurley2000(p_CoreMass);
    const double mu   = (p_Mass - p_CoreMass) / (mHeF - p_CoreMass);
    const double tmp  = (b[18] * PPOW(mu, b[19])) / (1.0 + (b[18] + lZHe - lMin) / (lMin - lZHe) * std::exp(15.0 * (p_Mass - mHeF)));

    return lZHe + (((1.0 + b[20]) / (1.0 + (b[20] * PPOW(mu, 1.6479)))) * tmp);
}


/*
 * CalculateRemnantLuminosity_Hurley2000
 *
 * @brief
 * Calculate the luminosity of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. mass = coreMass), per Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantLuminosity_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const
 *
 * @param       p_Metallisity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol) (typically effective mass, mass0)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Luminosity of remnant core (Lsol)
 */
double GiantBranch::CalculateRemnantLuminosity_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_CoreMass) const {
    return (p_Mass > ZDEP->HurleyMCOHeF(p_Metallicity))
                ? HeMS::CalculateLuminosityAtZAHeMS_Hurley2000(p_CoreMass)
                : WhiteDwarfs::CalculateLuminosity_Hurley2000_Static(p_CoreMass, 0.0, Metallicity());
}
























/// GiantBranch_Constituent <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


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
inline double GiantBranch_Constituent::CalculateCriticalMassRatio_Hurley2002(const double p_Mass, const double p_CoreMass) const {
    return 1.0 / (0.362 + 1.0 / (3.0 * (1.0 - p_CoreMass / p_Mass)));   // Hurley et al. 2002, just after eq. 57
}

