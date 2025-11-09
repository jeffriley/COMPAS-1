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


    // member functions - alphabetically (sort of - some are grouped by functionality)
            double          CalculateConvectiveCoreMass() const                                             { return m_CoreMass; }
            
inline DBL_DBL CalculateConvectiveEnvelopeMass() const override {
    return 
};







 double CalculateCoreMassAtHeI_Hurley2000(const double p_Mass) const;



    static  double          CalculateCoreMassAtSN_Static(const double p_Mthreshold, const double p_McBAGB);





// radius


double          CalculateConvectiveCoreRadius() const { return std::min(CalculateRemnantRadius(), m_Radius); } // Last paragraph of section 6 of Hurley+ 2000








    

    virtual double          CalculateInitialSupernovaMass() const                                           { return m_Mass; }                                                  // Use class member variables

            double          CalculateLifetimeToHeIgnition(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB);







            double          CalculateMassLossRateHurley();






            
            
            
            double          CalculateFallbackFraction_Fryer2012_Delayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass);


            double          CalculateFallbackFractionRapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass);



            
            
            double          CalculateGravitationalRemnantMass(const double p_BaryonicRemnantMass);




            double          CalculateRemnantMassByMullerMandel(const double p_COCoreMass, const double p_HeCoreMass);
            
            double          CalculateRemnantMassBySchneider2020(const double p_COCoreMass, const bool p_UseSchneiderAlt = false);
            double          CalculateRemnantMassBySchneider2020Alt(const double p_COCoreMass)               { return CalculateRemnantMassBySchneider2020(p_COCoreMass, true); }
            double          CalculateRemnantNSMassMullerMandel(const double p_COCoreMass, const double p_HeCoreMass);


            double          CalculateMomentOfInertia() const;





            double          CalculateRadialExtentConvectiveEnvelope() const;







            double          CalculateRadiusOnMassChange(double p_dM)                                        { return CalculateRadiusOnPhase(m_Mass + p_dM, m_Luminosity); }




    virtual double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, m_BnCoefficients); }



            double          CalculateRadiusOnPhase() const                                                  { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }











	



inline double CalculateMLrateThermal() const { 
    return (m_StateHistory.CurrentState.Mass() - m_StateHistory.CurrentState.m_CoreMass()) / CalculateTimescale_Thermal();
}


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


DBL_VECTOR CalculateTimescales_Hurley2000() const override {
    return CalculateTimescales_Hurley2000(
        m_StateHistory.CurrentState.MassEffectiveInitial(),
        m_StateHistory.CurrentState.GBparams(),
        m_StateHistory.CurrentState.TimeScales()
    );
}
COMPAS_PURE DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;

DBL_VECTOR CalculateGBparams_Hurley2000() const override {
    return CalculateGBparams_Hurley2000(m_StateHistory.CurrentState.MassEffectiveInitial(), m_StateHistory.CurrentState.GBparams());
}
COMPAS_PURE virtual DBL_VECTOR CalculateGBparams_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const;






inline double CalculateEffectiveInitialMass_Hurley2000() const override { return BaseStar::CalculateEffectiveInitialMass_Hurley2000(); } // per Hurley et al. 2000, section 7.1


 
GNU_CONST inline STELLAR_TYPE GiantBranch::CalculateRemnantType_Muller2016(const double p_COCoreMass) const;

GNU_CONST inline double CalculatePerturbationMu_Hurley2000(const double p_Mass, const double p_Luminosity, const double p_CoreMass) const;


GNU_CONST double CalculateProtoCoreMass_Fryer2012_Delayed(const double p_COCoreMass) const;

GNU_CONST DBL_DBL CalculateRemnantMass_Belczynski2002(const double p_Mass, const double p_COCoreMass) const;

GNU_CONST DBL_DBL CalculateRemnantMass_Fryer2012(const double p_Mass, const double p_COCoreMass, const double p_NSmaxBaryonicMass, const SN_ENGINE p_SNenginePrescription) const;

GNU_CONST DBL_DBL CalculateRemnantMass_Fryer2022(const double    p_Mass,
                                                 const double    p_COCoreMass,
                                                 const double    p_NSmaxBaryonicMass,
                                                 const double    p_fMix,
                                                 const double    p_mCrit
                                                 const SN_ENGINE p_SNenginePrescription) const;




GNU_CONST DBL_DBL CalculateRemnantMass_Maltsev2025(const double p_COCoreMass, const double p_HeCoreMass);




GNU_CONST double CalculateRemnantMass_Muller2016(const double p_Mass, const double p_COCoreMass) const;


GNU_CONST double CalculateBaryonicRemnantMass_Fryer2012(const double p_ProtoMass, double p_FallbackMass) const;


GNU_CONST double CalculateFallbackMass_Fryer2012(const double p_PreSNMass, const double p_ProtoMass, const double p_FallbackFraction) const;


double CalculateBHMassAfterFallback_MullerMandel2020(const double p_COCoreMass, const double p_HeCoreMass) const;


GNU_CONST double GiantBranch::CalculateHRateConstant_Hurley2000(const double p_Mass) const;





COMPAS_PURE static double CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Luminosity);

GNU_CONST   static double CalculateRadiusOnZAHB_Hurley2000_Static(const double      p_Mass,
                                                                const double      p_CoreMass,
                                                                const double      p_MHeF,
                                                                const double      p_MinLuminosity,
                                                                const double      p_Alpha1,
                                                                const DBL_VECTOR& p_bN);

inline double CalculateRemnantRadius_Hurley2000() const override {
    return CalculateRemnantRadius_Hurley2000(m_StateHistory.CurrentState.MassEffectiveInitial(), m_StateHistory.CurrentState.CoreMass());
}
COMPAS_PURE double CalculateRemnantRadius_Hurley2000(const double p_Mass, const double p_CoreMass) const;


GNU_CONST double CalculateRadiusAtHeIgnition_Hurley2000(const double      p_Mass,
                                                        const double      p_CoreMass,
                                                        const double      p_MHeF,
                                                        const double      p_MFGB,
                                                        const double      p_MinLuminosity,
                                                        const double      p_Alpha1,
                                                        const DBL_VECTOR& p_bN) const;





// radius
CalculateRadiusOnPhase() const { return CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Luminosity); }



// Luminosity
COMPAS_PURE static double CalculateLuminosityAtBGB_Hurley2000_Static(const double p_Mass) const;
COMPAS_PURE static double CalculateLuminosityAtHeI_Hurley2000_Static(const double p_Mass);


GNU_CONST static double CalculateLuminosityOnZAHB_Hurley2000_Static(const double p_Mass, const double p_CoreMass);
GNU_CONST virtual double CalculateRemnantLuminosity_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MHeF) const;



// constants etc.
GNU_CONST virtual double CalculateCoreMass_Luminosity_B_Hurley2000(const double p_Mass) const;
GNU_CONST virtual double CalculateCoreMass_Luminosity_D_Hurley2000(const double p_Mass, const double p_ZetaHurley = 0.0, const double p_MHeF = 0.0) const;
GNU_CONST virtual double CalculateCoreMass_Luminosity_p_Hurley2000(const double p_Mass, const double p_MHeF) const;
GNU_CONST virtual double CalculateCoreMass_Luminosity_q_Hurley2000(const double p_Mass, const double p_MHeF) const;


GNU_CONST double CalculateCoreMass_Luminosity_Lx_Hurley2000(const DBL_VECTOR& p_GBparams) const;
GNU_CONST double CalculateCoreMass_Luminosity_Mx_Hurley2000(const DBL_VECTOR& p_GBparams) const;


// mass
GNU_CONST   double CalculateCoreMassAt2ndDredgeUp_Hurley2000(const double p_McBAGB) const;
COMPAS_PURE double CalculateCoreMassAtBAGB_Hurley2000(const double p_Mass) const;
GNU_CONST double CalculateCoreMassAtBGB_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const double p_MHef, const DBL_VECTOR& p_aN) const;



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<






};


// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


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
GNU_CONST inline double GiantBranch::CalculateHRateConstant_Hurley2000(const double p_Mass) const {
    return PPOW(10.0, std::max(-4.8, std::min((-5.7 + (0.8 * p_Mass)), (-4.1 + (0.14 * p_Mass)))));
}


/*
 * CalculatePerturbationMu_Hurley2000
 *
 * @brief
 * Calculate the small envelope perturbation parameter, mu, per Hurley et al. 2000,
 * eqs 97 & 98
 *
 *
 * double CalculatePerturbationMu_Hurley2000(const double p_Mass, const double p_Luminosity, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Small envelope perturbation parameter, mu
 */
GNU_CONST inline double GiantBranch::CalculatePerturbationMu_Hurley2000(const double p_Mass, const double p_Luminosity, const double p_CoreMass) const {
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
 * double CalculateCoreMass_Luminosity_Lx_Hurley2000(const DBL_VECTOR& p_GBparams) const
 *
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Core mass - Luminosity relation parameter, Lx
 */
GNU_CONST inline double GiantBranch::CalculateCoreMass_Luminosity_Lx_Hurley2000(const DBL_VECTOR& p_GBparams) const {
#define GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)] // for convenience and readability - undefined at end of function
    // since the mass used here is the mass at crossover (Mx), these
    // should give the same answer - but we'll take the minimum anyway
    return std::min((GBparams(B) * PPOW(GBparams(Mx), GBparams(q))), (GBparams(D) * PPOW(GBparams(Mx), GBparams(p))));
#undef GBparams
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
 * double CalculateCoreMass_Luminosity_Mx_Hurley2000(const DBL_VECTOR& p_GBparams) const
 *
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Core mass - Luminosity relation parameter, Mx
 */
GNU_CONST inline double GiantBranch::CalculateCoreMass_Luminosity_Mx_Hurley2000(const DBL_VECTOR &p_GBparams) const {
#define GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)] // for convenience and readability - undefined at end of function
    return PPOW(GBparams(B) / GBparams(D), (1.0 / (GBparams(p) - GBparams(q))));
#undef GBparams
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
GNU_CONST inline double GiantBranch::CalculateCoreMass_Luminosity_B_Hurley2000(const double p_Mass) const {
    return std::max(3.0E4, (500.0 + (1.75E4 * PPOW(p_Mass, 0.6))));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
GNU_CONST inline double GiantBranch::CalculateCoreMassAt2ndDredgeUp_Hurley2000(const double p_McBAGB) const {
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
COMPAS_PURE inline double GiantBranch::CalculateCoreMassAtBAGB_Hurley2000(const double p_Mass) const {
    const DBL_VECTOR b = GLOBALS->HurleyBcoefficients();                    // get Hurley b coefficients
    return std::sqrt(std::sqrt((b[36] * PPOW(p_Mass, b[37])) + b[38]));     // sqrt() is much faster than PPOW()
}


/*
 * CalculateCoreMassAtSN_Static
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
 * static double CalculateCoreMassAtSN_Static(const double p_mThreshold, const double p_McBAGB)
 *
 * @param       p_mThreshold                    Threshold mass (see notes above) (Msol)
 * @param       p_McBAGB                        Core mass at the Base of the Asymptotic Giant Branch (Msol)
 * @return                                      Maximum core mass pre-SN on the Asymptotic Giant Branch (Msol)
 */
GNU_CONST static inline double GiantBranch::CalculateCoreMassAtSN_Static(const double p_mThreshold, const double p_McBAGB) {
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
GNU_CONST inline double GiantBranch::CalculateFallbackMass_Fryer2012(const double p_PreSNMass, const double p_ProtoMass, const double p_FallbackFraction) const {
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
GNU_CONST inline double GiantBranch::CalculateProtoCoreMass_Fryer2012_Delayed(const double p_COCoreMass) const {
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
 * DBL_DBL CalculateRemnantMass_Belczynski2002(const double p_Mass, const double p_COCoreMass) const
 *
 * @param       p_Mass                          Pre-SN mass of the star (Msol)
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
GNU_CONST inline DBL_DBL GiantBranch::CalculateRemnantMass_Belczynski2002(const double p_Mass, const double p_COCoreMasson) const {
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
GNU_CONST inline double GiantBranch::CalculateBaryonicRemnantMass_Fryer2012(const double p_ProtoMass, double p_FallbackMass) const {
    return p_ProtoMass + p_FallbackMass;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 RADIUS FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRadius_Hurley2000_Static
 *
 * @brief
 * Calculate radius on the Giant Branch, per Hurley et al. 2000, eq 46
 *
 *
 * static double CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Luminosity)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      GB radius (Rsol)
 */
COMPAS_PURE inline double GiantBranch::CalculateRadius_Hurley2000_Static(const double p_Mass, const double p_Luminosity) {
    const DBL_VECTOR b = GLOBALS->HurleyBcoefficients(); // get Hurley b coefficients
    return std::min((b[4] * PPOW(p_Mass, -b[5])), (b[6] * PPOW(p_Mass, -b[7]))) * (PPOW(p_Luminosity, b[1]) + (b[2] * PPOW(p_Luminosity, b[3])));
}


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
COMPAS_PURE inline double GiantBranch::CalculateRemnantRadius_Hurley2000(const double p_Mass, const double p_CoreMass) const {
    const double mHeF = GLOBALS->HurleyMassCutoffs(static_cast<int>(MHeF));
    return p_Mass > mHeF ? HeMS::CalculateRadiusAtZAHeMS_Hurley2000_Static(p_CoreMass) : WhiteDwarfs::CalculateRadiusOnPhase_Marsh2004_Static(p_CoreMass);
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
double GiantBranch::CalculateMomentOfInertia(const double p_Mass, const double p_Radius, const double p_CoreMass) const {
    const double Rc = CalculateRemnantRadius_Hurley2000();
    return (0.1 * (p_Mass - p_CoreMass) * p_Radius * p_Radius) + (0.21 * m_CoreMass * Rc * Rc);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               LUMINOSITY FUNCTIONS                                //
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
COMPAS_PURE inline double GiantBranch::CalculateLuminosityAtBGB_Hurley2000_Static(const double p_Mass) const {

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
COMPAS_PURE static inline double GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(const double p_Mass) {

    const DBL_VECTOR b = GLOBALS->HurleyBcoefficients(); // get Hurley b coefficients
    const double mHeF  = GLOBALS->HurleyMassCutoffs(static_cast<int>(MHeF));
    
    return p_Mass < mHeF
            ? (b[9] * PPOW(p_Mass, b[10])) / (1.0 + (GLOBALS->HurleyAlpha1() * exp(15.0 * (p_Mass - mHeF))))
            : (b[11] + (b[12] * PPOW(p_Mass, 3.8))) / (b[13] + (p_Mass * p_Mass));
}


/*
 * CalculateLuminosityOnZAHB_Hurley2000_Static
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
 * static double CalculateLuminosityOnZAHB_Hurley2000_Static(const double p_Mass, const double p_CoreMass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      ZAHB luminosity (Lsol)
 */
GNU_CONST static inline double GiantBranch::CalculateLuminosityOnZAHB_Hurley2000_Static(const double p_Mass, const double p_CoreMass) {
#define b(x) GLOBALS->HurleyBCoefficients(x) // for convenience and readability - undefined at end of function

    const double mHef = GLOBALS->HurleyMassCutoffs(static_cast<int>(MHeF));
    const double mu   = (p_Mass - p_CoreMass) / (mHef - p_CoreMass);
    const double lZHe = HeMS::CalculateLuminosityAtZAHeMS_Hurley2000_Static(p_CoreMass);
    const double tmp  = (b(18) * PPOW(mu, b(19))) / (1.0 + (b(18) + lZHe - GLOBALS->MinLuminosity_CHeB()) / (GLOBALS->MinLuminosity_CHeB() - lZHe) * exp(15.0 * (p_Mass - mHef)));

    return lZHe + ((1.0 + b(20)) / (1.0 + (b(20) * PPOW(mu, 1.6479))) * tmp);

#undef bN
}


/*
 * CalculateRemnantLuminosity_Hurley2000
 *
 * @brief
 * Calculate the luminosity of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. mass = coreMass), per Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantLuminosity_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MHeF) const
 *
 * @param       p_Mass                          Mass of the star (Msol) (typically effective mass, mass0)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @return                                      Luminosity of remnant core (Lsol)
 */
COMPAS_PURE inline double GiantBranch::CalculateRemnantLuminosity_Hurley2000(const double p_Mass, const double p_CoreMass, const double p_MHeF) const {
    return p_Mass > p_MHeF
            ? HeMS::CalculateLuminosityAtZAHeMS_Hurley2000_Static(p_CoreMass)
            : WhiteDwarfs::CalculateLuminosityOnPhase_Hurley2000_Static(p_CoreMass, 0.0, GLOBALS->ReferenceMetallicity(), WD_Baryon_Number.at(STELLAR_TYPE::HELIUM_WHITE_DWARF));
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
 * @return                                      Remnant stellar type (STELLAR_TYPE)
 */
GNU_CONST inline STELLAR_TYPE GiantBranch::CalculateRemnantType_Muller2016(const double p_COCoreMass) const {

    STELLAR_TYPE stellarType;

         if (p_COCoreMass < 3.6 ) stellarType = STELLAR_TYPE::NEUTRON_STAR;
    else if (p_COCoreMass < 4.05) stellarType = STELLAR_TYPE::BLACK_HOLE;
    else if (p_COCoreMass < 4.6 ) stellarType = STELLAR_TYPE::NEUTRON_STAR;
    else if (p_COCoreMass < 5.7 ) stellarType = STELLAR_TYPE::BLACK_HOLE;
    else if (p_COCoreMass < 6.0 ) stellarType = STELLAR_TYPE::NEUTRON_STAR;
    else                          stellarType = STELLAR_TYPE::BLACK_HOLE;

    return stellarType;
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
