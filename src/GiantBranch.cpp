#include "GiantBranch.h"
#include "HeMS.h"
#include "WhiteDwarfs.h"
#include "ONeWD.h"
#include "NS.h"
#include "BH.h"

#include "stellarUtils.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                     COEFFICIENT AND CONSTANT CALCULATIONS ETC.                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateTimescales_Hurley2000
 *
 * @brief
 * (Re)calculate timescales given the mass of the star, per Hurley at al. 2000.
 * 
 * Since timescales depend on a star's mass, they need to be calculated whenever
 * the mass of the star changes (probably every timestep).
 *
 *
 * DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_Timescales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
COMPAS_PURE DBL_VECTOR GiantBranch::CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const {

// #defines for convenience and readability - undefined at end of function
#define GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)]
#define tScales(x) tScales[static_cast<int>(TIMESCALE::x)]

    const double p1   = GBparams(p) - 1.0;
    const double q1   = GBparams(q) - 1.0;
    const double p1_p = p1 / GBparams(p);
    const double q1_q = q1 / GBparams(q);

    const double lBGB = GiantBranch::CalculateLuminosityAtBGB_Hurley2000(p_Mass);

    DBL_VECTOR tScales = p_tScales; // copy given timescales

    // (re)calculate MS timescales
    DBL_VECTOR tScales = MainSequence::CalculateTimescales_Hurley2000(p_Mass, p_GBparams, tScales);

    tScales[static_cast<int>(TIMESCALE::tinf1_FGB)] = tScales(tBGB) + ((1.0 / (p1 * GBparams(AH) * GBparams(D))) * PPOW((GBparams(D) / lBGB), p1_p));
    tScales[static_cast<int>(TIMESCALE::tMx_FGB)]   = tScales(tinf1_FGB) - ((tScales(tinf1_FGB) - tScales(tBGB)) * PPOW((lBGB / GBparams(Lx)), p1_p));
    tScales[static_cast<int>(TIMESCALE::tinf2_FGB)] = tScales(tMx_FGB) + ((1.0 / (q1 * GBparams(AH) * GBparams(B))) * PPOW((GBparams(B) / GBparams(Lx)), q1_q));

    tScales[static_cast<int>(TIMESCALE::tHeI)]      = CalculateLifetimeToHeI_Hurley2000(p_Mass, tScales(tinf1_FGB), tScales(tinf2_FGB));
    tScales[static_cast<int>(TIMESCALE::tHeMS)]     = HeMS::CalculateLifetimeOnPhase_Hurley2000_Static(p_Mass);

    // return timescales vector by value - NRVO takes care of performance/efficiency
    return tScales;

#undef tScales
#undef GBparams
}



/*
 * CalculateGBparams_Hurley2000
 *
 * @brief
 * Calculate Giant Branch (GB) parameters, per Hurley et al. 2000
 *
 * Since Giant Branch parameters depend on a star's mass, they need to be calculated
 * whenever the mass of the star changes (probably every timestep).
 *
 *
 * DBL_VECTOR CalculateGBparams_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Mutated GB parameters
 */
COMPAS_PURE DBL_VECTOR GiantBranch::CalculateGBparams_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const {

    DBL_VECTOR GBparams = p_GBparams;   // copy given GBparams

    GBparams[static_cast<int>(HURLEY_GBP:::AH)]     = CalculateHRateConstant_Hurley2000(p_Mass);
    GBparams[static_cast<int>(HURLEY_GBP:::AHHe)]   = HHE_RATE_CONSTANT_HURLEY2000;
    GBparams[static_cast<int>(HURLEY_GBP:::AHe)]    = HE_RATE_CONSTANT_HURLEY2000;

    GBparams[static_cast<int>(HURLEY_GBP:::B)]      = CalculateCoreMass_Luminosity_B_Hurley2000(p_Mass);
    GBparams[static_cast<int>(HURLEY_GBP:::D)]      = CalculateCoreMass_Luminosity_D_Hurley2000(p_Mass);

    GBparams[static_cast<int>(HURLEY_GBP:::p)]      = CalculateCoreMass_Luminosity_p_Static_Hurley2000(p_Mass);
    GBparams[static_cast<int>(HURLEY_GBP:::q)]      = CalculateCoreMass_Luminosity_q_Static_Hurley2000(p_Mass);

    GBparams[static_cast<int>(HURLEY_GBP:::Mx)]     = CalculateCoreMass_Luminosity_Mx_Hurley2000(GBparams);
    GBparams[static_cast<int>(HURLEY_GBP:::Lx)]     = CalculateCoreMass_Luminosity_Lx_Hurley2000(GBparams);

    GBparams[static_cast<int>(HURLEY_GBP:::McBAGB)] = CalculateCoreMassAtBAGB_Hurley2000(p_Mass);
    GBparams[static_cast<int>(HURLEY_GBP:::McDU)]   = CalculateCoreMassAt2ndDredgeUp_Hurley2000(GBparams[static_cast<int>(HURLEY_GBP:::McBAGB)]);
    GBparams[static_cast<int>(HURLEY_GBP:::McBGB)]  = CalculateCoreMassAtBGB_Hurley2000(p_Mass, GBparams[static_cast<int>(HURLEY_GBP:::McBAGB)]);

    // return GB parameters vector by value - NRVO takes care of performance/efficiency
    return GBparams;
}


/*
 * Calculate Giant Branch (GB) parameters per Hurley et al. 2000
 *
 * Giant Branch Parameters depend on a star's mass, so this needs to be called at least each timestep
 *
 * Vectors are passed by reference here for performance - preference would be to pass const& and
 * pass modified value back by functional return, but this way is faster - and this function is
 * called many, many times.
 *
 *
 * This function exists to facilitate the calculation of gbParams in EAGB::ResolveEnvelopeLoss() for the
 * stellar type to which the star will evolve.  The calculations of some of the stellar attributes there
 * depend on new gbParams.  
 * JR: This really needs to be revisited one day - these calculations should really be performed after
 *     switching to the new stellar type, but other calculations are done (in the legacy code) before the 
 *     switch (see evolveOneTimestep() in star.cpp for EAGB stars in the legacy code).
 *
 *
 * void CalculateGBparams_Static(const double      p_Mass, 
 *                               const double      p_LogMetallicityXi, 
 *                               const DBL_VECTOR &p_MassCutoffs, 
 *                               const DBL_VECTOR &p_An, 
 *                               const DBL_VECTOR &p_Bn, 
 *                                     DBL_VECTOR &p_GBparams) const
 *
 * @param   [IN]        p_Mass                  Mass in Msol
 * @param   [IN]        p_LogMetallicityXi      log10(Metallicity / Zsol) - called xi in Hurley et al. 2000
 * @param   [IN]        p_MassCutoffs           Mass cutoffs
 * @param   [IN]        p_An                    a(n) coefficients
 * @param   [IN]        p_Bn                    b(n) coefficients
 * @param   [IN/OUT]    p_GBparams              Giant Branch Parameters - calculated here
 */
void GiantBranch::CalculateGBparams_Static(const double      p_Mass, 
                                           const double      p_LogMetallicityXi, 
                                           const DBL_VECTOR &p_MassCutoffs, 
                                           const DBL_VECTOR &p_An, 
                                           const DBL_VECTOR &p_Bn, 
                                                 DBL_VECTOR &p_GBparams) const {
#define p_GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)]    // for convenience and readability - undefined at end of function

    p_GBparams(AH)     = CalculateHRateConstant_Hurley2000_Static(p_Mass);
    p_GBparams(AHHe)   = HHE_RATE_CONSTANT_HURLEY2000;
    p_GBparams(AHe)    = HE_RATE_CONSTANT_HURLEY2000;

    p_GBparams(B)      = CalculateCoreMass_Luminosity_B_Hurley2000_Static(p_Mass);
    p_GBparams(D)      = CalculateCoreMass_Luminosity_D_Hurley2000_Static(p_Mass, p_LogMetallicityXi, p_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)]);

    p_GBparams(p)      = CalculateCoreMass_Luminosity_p_Static(p_Mass, p_MassCutoffs);
    p_GBparams(q)      = CalculateCoreMass_Luminosity_q_Static(p_Mass, p_MassCutoffs);

    p_GBparams(Mx)     = CalculateCoreMass_Luminosity_Mx_Hurley2000_Static(p_GBparams);
    p_GBparams(Lx)     = CalculateCoreMass_Luminosity_Lx_Hurley2000_Static(p_GBparams);

    p_GBparams(McDU)   = CalculateCoreMassAt2ndDredgeUp_Static(p_GBparams(McBAGB));
    p_GBparams(McBAGB) = CalculateCoreMassAtBAGB_Static(p_Mass, p_Bn);
    p_GBparams(McBGB)  = CalculateCoreMassAtBGB_Static(p_Mass, p_MassCutoffs, p_An, p_GBparams);

#undef p_GBparams
}





/*
 * CalculateCoreMass_Luminosity_D_Hurley2000
 *
 * @brief
 * Calculate the Core mass - Luminosity relation parameter D, per Hurley et al. 2000, eqs 31 - 38
 *
 *
 * double CalculateCoreMass_Luminosity_D_Hurley2000(onst double p_Mass, const double p_ZetaHurley = 0.0, const double p_MHeF = 0.0) const 
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_ZetaHurley                    Hurley zeta value (log10(Z / ZSOL_HURLEY))
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @return                                      Core mass - Luminosity relation parameter D
 */
double GiantBranch::CalculateCoreMass_Luminosity_D_Hurley2000(const double p_Mass, const double p_ZetaHurley = 0.0, const double p_MHeF = 0.0) const {

    const double D0 = 5.37 + (0.135 * p_ZetaHurley);
    const double D1 = (0.975 * D0) - (0.18 * p_Mass);

    double logD = D0;                                           // default value

    if (p_Mass > p_MHeF) {
        if (p_Mass >= 2.5) {
            logD = std::max(std::max(-1.0, D1), (0.5 * D0) - (0.06 * p_Mass));
        }
        else {                                                  // Linear interpolation between end points
            const double gradient  = (D0 - D1) / (p_MHeF - 2.5);
            const double intercept = D0 - (p_MHeF * gradient);
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
 * double CalculateCoreMass_Luminosity_p_Hurley2000(const double p_Mass, const double p_MHeF) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @return                                      Core mass - Luminosity relation parameter, p
 */
double GiantBranch::CalculateCoreMass_Luminosity_p_Hurley2000(const double p_Mass, const double p_MHeF) const {

    double p = 6.0;                                             // default value

    if (p_Mass > p_MHeF) {
        if (p_Mass >= 2.5) {
            p = 5.0;
        }
        else {                                                  // linear interpolation between end points
            const double gradient  = 1.0 / (p_MHeF - 2.5);      // will be negative
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
 * double CalculateCoreMass_Luminosity_q_Hurley2000(const double p_Mass, const double p_MHeF) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @return                                      Core mass - Luminosity relation parameter q
 */
double GiantBranch::CalculateCoreMass_Luminosity_q_Hurley2000(const double p_Mass, const double p_MHeF) const {

    double q = 3.0;                                             // default value

    if (p_Mass > p_MHeF) {
        if (p_Mass >= 2.5) {
            q = 2.0;
        }
        else {                                                  // linear interpolation between end points
            double gradient  = 1.0 / (p_MHeF - 2.5);
            double intercept = 2.0 - (2.5 * gradient);
            q = (gradient * p_Mass) + intercept;
        }
    }

    return q;
}







/// JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

/*
 * Perturb Luminosity and Radius
 *
 * Perturbs Luminosity and Radius per Hurley et al. 2000, section 6.3
 * The attributes of the star are updated.
 *
 * Perturbation is disabled by default when DEBUG is enabled - except when
 * DEBUG_PERTURB is defined (see below).  The stellar class perturbation
 * functions are defined away if DEBUG is defined - so the generic Star
 * function is called (and does nothing). (So far FGB is the only class that
 * defines this function where it actually does anything)
 *
 * If DEBUG_PERTURB is defined then perturbation is not disabled while debugging.
 * To enable perturbation while DEBUG is enabled, define DEBUG_PERTURB.
 *
 *
 * void PerturbLuminosityAndRadius()
 */
#if !defined(DEBUG) || defined(DEBUG_PERTURB)   // don't perturb if debugging and DEBUG_PERTURB not defined
void GiantBranch::PerturbLuminosityAndRadius() {

    if (utils::Compare(m_Mu, 1.0) < 0) {   // perturb only if mu < 1.0

        double Lc = CalculateRemnantLuminosity();
        double Rc = CalculateRemnantRadius();

        double s = CalculateHurleyPerturbationS(m_Mass, m_Mu);
        double r = CalculateHurleyPerturbationR(m_Mass, m_Radius, m_Mu, Rc);

        m_Luminosity = Lc * PPOW((m_Luminosity / Lc), s);        
        m_Radius     = Rc * PPOW((m_Radius / Rc), r);
    }
}
#else
void GiantBranch::PerturbLuminosityAndRadius() { }
#endif






///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
 * double CalculateConvectiveEnvelopeRadialExtent_Hurley2002() const
 *
 * @param       p_Radius                        Radius of the star (Rsol)
 * @return                                      Radial extent of the convective outer envelope (Rsol)
 */
double GiantBranch::CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Radius) const {

    double envMass;
    double envZAMSMass;
    std::tie(envMass, envZAMSMass) = CalculateConvectiveEnvelopeMass();
    
    return envMass <= 0.0 || envZAMSMass <= 0.0 ? 0.0 : std::sqrt(envMass / envZAMSMass) * (p_Radius - CalculateConvectiveCoreRadius());          
}


/*
 * CalculateRadiusAtHeIgnition_Hurley2000
 *
 * @brief
 * Calculate radius at Helium Ignition, per Hurley et al. 2000, eq 50
 *
 *
 * double CalculateRadiusAtHeIgnition_Hurley2000(
 *     const double p_Mass,
 *     const double p_CoreMass,
 *     const double p_MHeF,
 *     const double p_MFGB,
 *     const double p_MinLuminosity
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_MFGB                          Maximum initial mass at helium ignition on the FGB (Hurley masscutoffs[MFGB]) (Msol)
 * @param       p_MinLuminosity                 Minimum luminosity on phase (Lsol)
 * @param       p_Alpha1                        Hurley alpha1 constant
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      Radius at Helium Ignition (Rsol)
 */
double GiantBranch::CalculateRadiusAtHeIgnition_Hurley2000(
    const double p_Mass,
    const double p_CoreMass,
    const double p_MHeF,
    const double p_MFGB,
    const double p_MinLuminosity
) const {

    double radius;

    const double lHeI = CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass);

    if (p_Mass <= p_MFGB) {
        radius = CalculateRadiusOnPhase(p_Mass, lHeI);
    }
    else if (p_Mass >= std::max(p_MFGB, HIGH_MASS_THRESHOLD)) {
        const double rmHe = CHeB::CalculateMinimumRadiusOnPhase_Static(p_Mass, p_CoreMass, p_Alpha1, p_MHeF, p_MFGB, p_MinLuminosity);
        radius = std::min(rmHe, EAGB::CalculateRadiusOnPhase_Static(p_Mass, lHeI, p_MHeF));
    }
    else {
        const double rmHe = CHeB::CalculateMinimumRadiusOnPhase_Static(p_Mass, p_CoreMass, p_Alpha1, p_MHeF, p_MFGB, p_MinLuminosity);
        radius = rmHe * PPOW((CalculateRadiusOnPhase(p_Mass, lHeI); / rmHe), log10(p_Mass / HIGH_MASS_THRESHOLD) / log10(p_MFGB / HIGH_MASS_THRESHOLD));
    }

    return radius;
}


/*
 * CalculateRadiusOnZAHB_Hurley2000_Static
 *
 * @brief
 * Calculate radius on the Zero Age Horizontal Branch, per Hurley et al. 2000, eq 54
 *
 *
 * static double CalculateRadiusOnZAHB_Hurley2000_Static(const double      p_Mass,
 *                                                       const double      p_CoreMass,
 *                                                       const double      p_MHeF,
 *                                                       const double      p_MinLuminosity,
 *                                                       const double      p_Alpha1,
 *                                                       const DBL_VECTOR& p_bN)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param   [IN]    p_MinLuminosity  Minimum luminosity on phase (only required for CHeB stars) - calculated once per star <<<<<<<<<<<<<<<<<<<<<<<<<<< ???????????????
 * @param       p_MinLuminosity                 Minimum luminosity on phase (Lsol)
 * @param       p_Alpha1                        Hurley alpha1 constant
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      Radius on the Zero Age Horizontal Branch (Rsol)
 */
double GiantBranch::CalculateRadiusOnZAHB_Hurley2000_Static(const double      p_Mass,
                                                            const double      p_CoreMass,
                                                            const double      p_MHeF) {

    const double rZHe  = HeMS::CalculateRadiusAtZAHeMS_Hurley2000_Static(p_CoreMass);
    const double lZAHB = CalculateLuminosityOnZAHB_Hurley2000_Static(p_Mass, p_CoreMass);
    const double rGB   = CalculateRadiusOnPhase_Hurley2000_Static(p_Mass, lZAHB, p_bN);
    const double f     = ((1.0 + p_bN[21]) * PPOW((p_Mass - p_CoreMass) / (p_MHeF - p_CoreMass), p_bN[22])) / (1.0 + p_bN[21] * PPOW(mu, p_bN[23]));

    return ((1.0 - f)) * rZHe + (f * rGB);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


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
COMPAS_PURE double GiantBranch::CalculateCoreMassAtBGB_Hurley2000(const double p_Mass, const double p_McBAGB) const {

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
 * double CalculateCoreMassAtHeI_Hurley2000(const double      p_Mass,
 *                                          const DBL_VECTOR& p_GBparams,
 *                                          const double      p_MHef,
 *                                          const double      p_Alpha1,
 *                                          const DBL_VECTOR& p_bN) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_MHef                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_Alpha1                        Hurley alpha1 constant
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      Core mass at HeI (Msol)
 */
double GiantBranch::CalculateCoreMassAtHeI_Hurley2000(const double      p_Mass,
                                                      const DBL_VECTOR& p_GBparams,
                                                      const double      p_MHef,
                                                      const double      p_Alpha1,
                                                      const DBL_VECTOR& p_bN) const {
    double coreMass;

    if (p_Mass < p_MHef) {
        const double luminosity = CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass);
        coreMass = BaseStar::CalculateCoreMass_Hurley2000_Static(luminosity, p_GBparams);
    }
    else {
        const double lMHeF   = CalculateLuminosityAtHeI_Hurley2000_Static(p_MHef);
        const double Mc_MHeF = BaseStar::CalculateCoreMass_Hurley2000_Static(lMHeF, p_GBparams);
        const double McBAGB  = CalculateCoreMassAtBAGB_Hurley2000(p_Mass, p_bN);
        const double c       = (utils::intPow(Mc_MHeF, 4)) - (MC_L_C1 * PPOW(p_MHef, MC_L_C2));

        coreMass = std::min((0.95 * McBAGB), std::sqrt(std::sqrt(c + (MC_L_C1 * PPOW(p_Mass, MC_L_C2))))); // sqrt() is much faster than PPOW()
    }

    return coreMass;
}




/*
 * Calculate the dominant mass loss mechanism and associated rate for the star
 * at the current evolutionary phase.
 *
 * According to Hurley et al. 2000
 *
 * double CalculateMassLossRateHurley()
 *
 * @return                                      Mass loss rate in Msol per year
 */
double GiantBranch::CalculateMassLossRateHurley() {

    double rateNJ = CalculateMLrate_NieuwenhuijzenDeJager1990(p_Mass, p_Radius, p_Luminosity);
    double rateKR = CalculateMLrate_KudritzkiReimers1978(p_Mass, p_Radius, p_Luminosity);
    double rateWR = CalculateMLrateWR_Hurley2000(p_Luminosity, p_Mu);

    m_DominantMassLossRate = MASS_LOSS_TYPE::GB;
    double dominantRate    = std::max(rateNJ, rateKR);

    if (utils::Compare(rateWR, dominantRate) > 0) {
        dominantRate           = rateWR;
        m_DominantMassLossRate = MASS_LOSS_TYPE::WR;
    }

    return dominantRate;
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLifetimeToHeI_Hurley2000
 *
 * @brief
 * Calculate the time to Helium ignition tHeI, per Hurley et al. 2000, eq 43
 *
 *
 * double CalculateLifetimeToHeI_Hurley2000(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB)
 *
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tinf1_FGB                 tinf1_FGB (Timescales[TIMESCALE::tinf1_FGB]) - First Giant Branch tinf1
 * @param   [IN]    p_Tinf2_FGB                 tinf2_FGB (Timescales[TIMESCALE::tinf2_FGB]) - First Giant Branch tinf2
 * @return                                      Lifetime to He ignition (tHeI)
 */
double GiantBranch::CalculateLifetimeToHeI_Hurley2000(const double p_Mass, const double p_Tinf1_FGB, const double p_Tinf2_FGB) {

    double LHeI = CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass);
    double p1   = gbParams(p) - 1.0;
    double q1   = gbParams(q) - 1.0;

    return LHeI <= gbParams(Lx)
            ? p_Tinf1_FGB - (1.0 / (p1 * gbParams(AH) * gbParams(D))) * PPOW((gbParams(D) / LHeI), (p1 / gbParams(p)))
            : p_Tinf2_FGB - (1.0 / (q1 * gbParams(AH) * gbParams(B))) * PPOW((gbParams(B) / LHeI), (q1 / gbParams(q)));
}






///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              SUPERNOVA CALCULATIONS                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



/*
 * CalculateRemnantMass_Schneider2020
 *
 * @brief
 * Calculate remnant mass for isolated single stars, per Schneider et al. 2020.
 * (See arxiv:2008.08599).
 * 
 * Function provides for standard and alternative fits (single stars) from Schneider et al. 2020.
 * 
 *
 * DBL_DBL CalculateRemnantMass_Schneider2020(const double p_Mass,
 *                                            const double p_COCoreMass,
 *                                            const double p_HeCoreMassPreSN,
 *                                            const double p_MaxNSMass,
 *                                            const REMNANT_MASS_PRESCRIPTION p_RemnantMassPrescription) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMassPreSN               Pre-SN (at CO formation) He core mass (Msol)
 * @param       p_MaxNSMass                     User-specified maximum mass of a neutron star (`--maximum-neutron-star-mass`) (Msol)
 * @param       p_RemnantMassPrescription       The remnant mass prescription to use (should be SCHNEIDER2020 or SCHNEIDER2020ALT) 
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
DBL_DBL GiantBranch::CalculateRemnantMass_Schneider2020(const double                    p_Mass,
                                                        const double                    p_COCoreMass,
                                                        const double                    p_HeCoreMassPreSN,
                                                        const double                    p_MaxNSMass,
                                                        const REMNANT_MASS_PRESCRIPTION p_RemnantMassPrescription) const {
    double logRemnantMass;
    switch (p_RemnantMassPrescription) {                                // which REMNANT_MASS_PRESCRIPTION?

        REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020:                       // SCHNEIDER2020 standard prescription

                 if (p_COCoreMass <  6.357) { logRemnantMass = log10(0.03357 * p_COCoreMass + 1.31780); }
            else if (p_COCoreMass <  7.311) { logRemnantMass = -0.02466 * p_COCoreMass + 1.28070; }
            else if (p_COCoreMass < 12.925) { logRemnantMass = log10(0.03357 * p_COCoreMass + 1.31780); }
            else                            { logRemnantMass = 0.01940 * p_COCoreMass + 0.98462; }
            break;

        REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020ALT:                    // SCHNEIDER2020 alternative prescription

                 if (p_COCoreMass <  6.357) { logRemnantMass = log10(0.04199 * p_COCoreMass + 1.28128); }
            else if (p_COCoreMass <  7.311) { logRemnantMass = -0.02466 * p_COCoreMass + 1.28070; }
            else if (p_COCoreMass < 12.925) { logRemnantMass = log10(0.04701 * (p_COCoreMass * p_COCoreMass) - 0.91403 * p_COCoreMass + 5.93380); }
            else                            { logRemnantMass = 0.01940*p_COCoreMass + 0.98462; }
            break;

        default:                                                        // unexpected prescription
            // the only way this can happen is if the REMNANT_MASS_PRESCRIPTION passed to this function
            // is not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.

            THROW_ERROR(ERROR::UNEXPECTED_REMNANT_MASS_PRESCRIPTION);   // throw error  
    }

    // convert remnant mass to linear value, and limit to the pre-SN He Core mass
    const double remnantMass = std::min(PPOW(10.0, logRemnantMass), p_HeCoreMassPreSN);

    // fallback fraction of mass beyond proto-neutron-star for BH formation and kicks
    const double fallbackFraction = remnantMass > p_MaxNSMass ? (remnantMass - NEUTRON_STAR_MASS) / (p_Mass - NEUTRON_STAR_MASS) : 0.0;		

    return std::make_tuple(remnantMass, std::max(0.0, std::min(1.0, fallbackFraction)));
}


/*
 * CalculateRemnantMass_Maltsev2025  *** JR: should we rename this? *Ilya*
 *
 * @brief
 * Calculate remnant mass for isolated single stars, per Maltsev et al. 2025.
 * (See https://doi.org/10.1051/0004-6361/202554931)
 * 
 *
 * double CalculateRemnantMass_Maltsev2025(const double       p_COCoreMass,
 *                                         const double       p_HeCoreMass,
 *                                         const double       p_ZetaAsplund, 
 *                                         const double       p_MaltsevFallbackFraction,
 *                                         const MALTSEV_MODE p_MaltsevMode) const
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @param       p_ZetaAsplund                   Asplund zeta value (log10(Z / ZSOL_ASPLUND))
 * @param       p_MaltsevFallbackFraction       User-specfifid Maltsev fallback fraction (`--maltsev-fallback-fraction`) [0.0, 1.0]
 * @param       p_MaltsevMode                   User-specified Maltsev mode (`--maltsev-mode`))
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
double GiantBranch::CalculateRemnantMass_Maltsev2025(const double       p_COCoreMass,
                                                     const double       p_HeCoreMass,
                                                     const double       p_ZetaAsplund, 
                                                     const double       p_MaltsevFallbackFraction, //OPTIONS->MaltsevFallbackFraction()  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                                     const MALTSEV_MODE p_MaltsevMode) const { // OPTIONS->MaltsevMode() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    double fallbackFraction;
    double remnantMass;

    if (p_COCoreMass < MALTSEV2025_MMIN) {                                              // CO core mass < Maltsev lower threshold?
                                                                                        // yes - NS formation regardless of metallicity
        fallbackFraction = 0.0;
        remnantMass = NEUTRON_STAR_MASS;
    }
    else if (p_COCoreMass > MALTSEV2025_MMAX) {                                         // no - CO core mass > Maltsev upper threshold?
                                                                                        // yes - BH formation regardless of metallicity
        fallbackFraction = 1.0;
        remnantMass = p_HeCoreMass;
    }
    else {                                                                              // no - calculate remnant mass

        double Zbound;                                                                      
        switch (p_MaltsevMode) {                                                        // which Maltsev mode?                                                                                 

            case MALTSEV_MODE::OPTIMISTIC:                                              // OPTIMISTIC               
                Zbound = p_ZetaAsplund;                                                 // Asplund zeta (log10(Z / ZSOL_ASPLUND))
                break;

            case MALTSEV_MODE::BALANCED:                                                // BALANCED                        
                Zbound = std::min(std::max(p_ZetaAsplund, -1.6989700043360185), 0.0);   // Asplund zeta clamped to [log10(1/50), log10(1)]
                break;

            case MALTSEV_MODE::PESSIMISTIC:                                             // Pessimistic                          
                Zbound = std::min(std::max(p_ZetaAsplund, -1.0), 0.0);                  // Asplund zeta clamped to [log10(1/10), log10(1)]
                break;

            default:                                                                    // unexpected mode
                // the only way this can happen is if the MALTSEV_MODE passed to this function
                // isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the code passed a mode that this function doesn't account for, and
                // that should be flagged as an error and result in termination of the evolution of the star
                // or binary.
                // The correct fix for this is to add code to this function for the missing mode, or fix the
                // calling code to pass a mode that is handled by this function.
                THROW_ERROR(ERROR::UNEXPECTED_MALTSEV_MODE);                            // throw error
        }
  
        const double M1 = MALTSEV2025_M1S + (MALTSEV2025_M1S - MALTSEV2025_M1SZ01) * Zbound;
        const double M2 = MALTSEV2025_M2S + (MALTSEV2025_M2S - MALTSEV2025_M2SZ01) * Zbound;
        const double M3 = MALTSEV2025_M3S + (MALTSEV2025_M3S - MALTSEV2025_M3SZ01) * Zbound;

        if (p_COCoreMass >= M3 || (p_COCoreMass >= M1 && p_COCoreMass <= M2) ) {        // complete fallback into BH
            fallbackFraction = 1.0;
            remnantMass = p_HeCoreMass;
        }
        else if (p_COCoreMass > M2 && p_COCoreMass < M3 && RAND->Random(0, 1) <= 0.1) { // partial fallback - BH formation
            fallbackFraction = p_MaltsevFallbackFraction;
            remnantMass = (p_HeCoreMass - NEUTRON_STAR_MASS) * fallbackFraction + NEUTRON_STAR_MASS;
        }
        else {                                                                          // no fallback - NS formation
            fallbackFraction = 0.0;
            remnantMass = NEUTRON_STAR_MASS;
        }
    }

    return std::make_tuple(remnantMass, std::max(0.0, std::min(1.0, fallbackFraction)));
}



/*
 * CalculateRemnantMass_MullerMandel2020
 *
 * @brief
 * Calculate remnant mass, per Mandel & Mueller 2020
 *
 *
 * DBL_DBL CalculateRemnantMass_MullerMandel2020 (const double p_COCoreMass, const double p_HeCoreMass)
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
DBL_DBL GiantBranch::CalculateRemnantMass_MullerMandel2020(const double p_COCoreMass, const double p_HeCoreMass) {

    double remnantMass       = 0.0;   
    double pBH               = 0.0;
    double pCompleteCollapse = 0.0;
   
    if (utils::Compare(p_COCoreMass, MULLERMANDEL_M1) < 0 || utils::Compare(p_HeCoreMass, OPTIONS->MaximumNeutronStarMass()) <= 0 )
	    pBH = 0.0;
    else if (utils::Compare(p_COCoreMass, MULLERMANDEL_M3) < 0)
    	pBH = 1.0 / (MULLERMANDEL_M3-MULLERMANDEL_M1) * (p_COCoreMass-MULLERMANDEL_M1);
    else
	    pBH = 1.0;

    if (utils::Compare(RAND->Random(0, 1), pBH) < 0) {  // this is a BH
        if (utils::Compare(p_COCoreMass, MULLERMANDEL_M4) < 0)
		    pCompleteCollapse = 1.0 / (MULLERMANDEL_M4 - MULLERMANDEL_M1) * (p_COCoreMass - MULLERMANDEL_M1);
        else
		    pCompleteCollapse = 1.0;

	    if (utils::Compare(RAND->Random(0, 1), pCompleteCollapse) < 0)
		    remnantMass = p_HeCoreMass;
	    else
            remnantMass = CalculateBHMassAfterFallback_MullerMandel2020(p_COCoreMass, p_HeCoreMass);
    }
    else                                               // this is an NS
        remnantMass = CalculateRemnantNSMassMullerMandel(p_COCoreMass, p_HeCoreMass);

    // fallback fraction = 0.0 - no subsequent kick adjustment by fallback fraction needed
    return std::make_tuple(remnantMass, 0.0);
}


/*
 * Calculate NS remnant mass given COCoreMass and HeCoreMass
 *
 * Mandel & Mueller, 2020; also used by Maltsev2025, so separated out into helper function
 *
 *
 * double CalculateRemnantNSMassByMullerMandel (const double p_COCoreMass, const double p_HeCoreMass)
 *
 * @param   [IN]    p_COCoreMass                COCoreMass in Msol
 * @param   [IN]    p_HeCoreMass                HeCoreMass in Msol
 * @return                                      Remnant mass in Msol
 */
double GiantBranch::CalculateRemnantNSMassMullerMandel(const double p_COCoreMass, const double p_HeCoreMass) {

    double remnantMass      = 0.0;
    std::size_t iterations  = 0;

    if (utils::Compare(p_COCoreMass, MULLERMANDEL_M1) < 0) {
        while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS            &&
               (utils::Compare(remnantMass, MULLERMANDEL_MINNS) < 0                ||
                utils::Compare(remnantMass, OPTIONS->MaximumNeutronStarMass()) > 0 ||
                utils::Compare(remnantMass, p_HeCoreMass) > 0)) {
            remnantMass = MULLERMANDEL_MU1 + RAND->RandomGaussian(MULLERMANDEL_SIGMA1);
        }
        if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS) // failure to find a solution implies a narrow range; just pick a midpoint in this case
            remnantMass = (std::min(OPTIONS->MaximumNeutronStarMass(), p_HeCoreMass) + MULLERMANDEL_MINNS) / 2.0;
    }
    else if (utils::Compare(p_COCoreMass, MULLERMANDEL_M2) < 0) {
        while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS            &&
               (utils::Compare(remnantMass, MULLERMANDEL_MINNS) < 0                ||
                utils::Compare(remnantMass, OPTIONS->MaximumNeutronStarMass()) > 0 ||
                utils::Compare(remnantMass, p_HeCoreMass) > 0)) {
            remnantMass = MULLERMANDEL_MU2A + MULLERMANDEL_MU2B / (MULLERMANDEL_M2 - MULLERMANDEL_M1) * (p_COCoreMass - MULLERMANDEL_M1) + RAND->RandomGaussian(MULLERMANDEL_SIGMA2);
        }
        if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS) // failure to find a solution implies a narrow range; just pick a midpoint in this case
            remnantMass = (std::min(OPTIONS->MaximumNeutronStarMass(), p_HeCoreMass) + MULLERMANDEL_MINNS) / 2.0;
    }
    else {
        while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS            &&
               (utils::Compare(remnantMass, MULLERMANDEL_MINNS) < 0                ||
                utils::Compare(remnantMass, OPTIONS->MaximumNeutronStarMass()) > 0 ||
                utils::Compare(remnantMass, p_HeCoreMass) > 0)) {
            remnantMass = MULLERMANDEL_MU3A + MULLERMANDEL_MU3B / (MULLERMANDEL_M3 - MULLERMANDEL_M2) * (p_COCoreMass - MULLERMANDEL_M2) + RAND->RandomGaussian(MULLERMANDEL_SIGMA3);
        }
        if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS) // failure to find a solution implies a narrow range; just pick a midpoint in this case
            remnantMass = (std::min(OPTIONS->MaximumNeutronStarMass(), p_HeCoreMass) + MULLERMANDEL_MINNS) / 2.0;
    }
    return remnantMass;
}


/*
 * CalculateBHMassAfterFallback_MullerMandel2020
 *
 * @brief
 * Calculate the mass of the black hole created in a fallback supernova, per Mandel & Mueller, 2020
 * 
 * Also used by Maltsev2025
 *
 *
 * double CalculateBHMassAfterFallback_MullerMandel2020(const double p_COCoreMass, const double p_HeCoreMass) const
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @return                                      Remnant mass (Msol)
 */
double GiantBranch::CalculateBHMassAfterFallback_MullerMandel2020(const double p_COCoreMass, const double p_HeCoreMass) const {
    
    double remnantMass;

    std::size_t iterations = 0;    
    while (iterations++ < MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS &&
          (remnantMass < OPTIONS->MaximumNeutronStarMass() || remnantMass > p_HeCoreMass)) {
        remnantMass = MULLERMANDEL_MUBH * p_COCoreMass + RAND->RandomGaussian(MULLERMANDEL_SIGMABH);
    }
    
    if (iterations >= MULLERMANDEL_REMNANT_MASS_MAX_ITERATIONS) {   // no convergence?
        // failure to find a solution implies a narrow range; just pick a midpoint in this case
        remnantMass = (OPTIONS->MaximumNeutronStarMass() + p_HeCoreMass) / 2.0;
    }
    
    return remnantMass;
}


/*
 * CalculateRemnantMass_Muller2016
 *
 * @brief
 * Calculate remnant mass, per Muller et al. 2016
 * (as presented in eq B4 of Vigna-Gomez et al. 2018 (See arXiv:1805.07974))
 *
 * 
 * DBL_DBL CalculateRemnantMass_Muller2016(const double p_Mass, const double p_COCoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
DBL_DBL GiantBranch::CalculateRemnantMass_Muller2016(const double p_Mass, const double p_COCoreMass) const {

    double	remnantMass;

    if (p_COCoreMass < 1.372) {
        // not explicitly pointed out in Appendix B of Vigna-Gomez et al. 2018, but assumed for continuity and simplicity
        // Muller et al. 2016 didn't go as low as this in CO Core mass (see Figure A1)
        remnantMass = 1.21;
    }                                      
	else if (p_COCoreMass < 1.49) remnantMass = 1.21 - (0.4  * (p_COCoreMass - 1.372));
	else if (p_COCoreMass < 1.65) remnantMass = 1.16;
    else if (p_COCoreMass < 2.4 ) remnantMass = 1.32 + (0.3  * (p_COCoreMass - 1.65));
    else if (p_COCoreMass < 3.2 ) remnantMass = 1.42 + (0.7  * (p_COCoreMass - 2.4));
    else if (p_COCoreMass < 3.6 ) remnantMass = 1.32 + (0.25 * (p_COCoreMass - 3.2));
    else if (p_COCoreMass < 4.05) remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; // going to be a Black Hole
    else if (p_COCoreMass < 4.6 ) remnantMass = 1.5;
    else if (p_COCoreMass < 5.7 ) remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; // going to be a Black Hole
    else if (p_COCoreMass < 6.0 ) remnantMass = 1.64 - (0.2  * (p_COCoreMass - 5.7));
    else                          remnantMass = p_Mass * NEUTRINO_LOSS_FALLBACK_FACTOR; // Going to be a Black Hole

    // fallback fraction = 0.0 - no subsequent kick adjustment by fallback fraction needed
    return std::make_tuple(remnantMass, 0.0);
}


/*
 * CalculateGravitationalRemnantMass_Fryer2012
 *
 * @brief
 * Calculate the gravitational mass of the remnant, per
 *
 *    - Fryer et al. 2012, eq 13 for Neutron Stars
 *    - An enhanced version of Fryer et al. 2012, eq 14 for Black Holes
 *
 *
 * double CalculateGravitationalRemnantMass_Fryer2012(const double p_BaryonicRemnantMass)
 *
 * @param       p_BaryonicRemnantMass           Baryonic remnant mass (Msol)
 * @param       p_NSmaxBaryonicMass             Maximum baryonic mass for NS maximum mass (Msol)
 * @return                                      Gravitational mass of the remnant (Msol)
 */
double GiantBranch::CalculateGravitationalRemnantMass_Fryer2012(const double p_BaryonicRemnantMass, const double p_NSmaxBaryonicMass) {

    double mass;
    
    if (p_BaryonicRemnantMass < p_NSmaxBaryonicMass) {                                      // below NS threshold?
        ERROR error;                                                                        // yes - neutron star
        std::tie(error, mass) = utils::SolveQuadratic(0.075, 1.0, -p_BaryonicRemnantMass);  // find gravitational mass
        if (error != ERROR::NONE) THROW_ERROR(error);                                       // throw error if no root found
    } 
    else {                                                                                  // no - black hole
        mass = BH::CalculateNeutrinoMassLoss_Static(p_BaryonicRemnantMass);                 // convert baryonic mass to gravitational mass
    }

    return mass;
}


/*
 * CalculateFallbackFraction_Fryer2012_Rapid
 *
 * @brief
 * Calculate the fraction of mass falling back onto the proto compact object,
 * using the rapid supernove mechanism, per Fryer et al. 2012, eq 15s & 16.
 *
 *
 * double CalculateFallbackFraction_Fryer2012_Rapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const
 *
 * @param       p_PreSNMass                     Pre-SN mass of the star (Msol)
 * @param       p_ProtoMass                     Fe/Ni proto compact object mass (Msol)
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
double GiantBranch::CalculateFallbackFraction_Fryer2012_Rapid(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const {

    double fallbackFraction;

    if (p_COCoreMass < 2.5) {
        fallbackFraction = 0.2 / (p_PreSNMass - p_ProtoMass);
    }
    else if (p_COCoreMass < 6.0) {
        fallbackFraction = ((0.286 * p_COCoreMass) - 0.514) / (p_PreSNMass - p_ProtoMass);
    }
    else if (p_COCoreMass < 7.0) {
        fallbackFraction = 1.0;
    }
    else if (p_COCoreMass < 11.0) {
        const double a1 = 0.25 - (1.275 / (p_PreSNMass - p_ProtoMass));
        fallbackFraction = (a1 * p_COCoreMass) + ((-11.0 * a1) + 1.0);
    }
    else {
        fallbackFraction = 1.0;
    }

    return std::max(0.0, std::min(1.0, fallbackFraction));
}



/*
 * CalculateFallbackFraction_Fryer2012_Delayed
 *
 * @brief
 * Calculate the fraction of mass falling back onto the proto compact object,
 * using the delayed supernova mechanis, per Fryer et al. 2012, eq 19.
 *
 *
 * double CalculateFallbackFraction_Fryer2012_Delayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const
 *
 * @param       p_PreSNMass                     Pre-SN mass of the star (Msol)
 * @param       p_ProtoMass                     Fe/Ni proto compact object mass (Msol)
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
double GiantBranch::CalculateFallbackFraction_Fryer2012_Delayed(const double p_PreSNMass, const double p_ProtoMass, const double p_COCoreMass) const {

    double fallbackFraction;

    if (p_COCoreMass < 2.5) {
        fallbackFraction = 0.2 / (p_PreSNMass - p_ProtoMass);
    }
    else if (p_COCoreMass < 3.5) {
        fallbackFraction = ((0.5 * p_COCoreMass) - 1.05) / (p_PreSNMass - p_ProtoMass);
    }
    else if (p_COCoreMass < 11.0) {
        const double a2 = 0.133 - (0.093 / (p_PreSNMass - p_ProtoMass));
        fallbackFraction = (a2 * p_COCoreMass) + ((-11.0 * a2) + 1.0);
    }
    else {
        fallbackFraction = 1.0;
    }

    return std::max(0.0, std::min(1.0, fallbackFraction));
}


/*
 * CalculateRemnantMass_Fryer2012
 *
 * @brief
 * Calculate remnant mass, per Fryer et al. 2012.
 *
 *
 * DBL_DBL CalculateRemnantMass_Fryer2012(const double p_Mass, const double p_COCoreMass, const double p_NSmaxBaryonicMass, const SN_ENGINE p_SNenginePrescription) const
 *
 * @param       p_Mass                          Pre-SN mass of the star (Msol)
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @param       p_NSmaxBaryonicMass             Maximum baryonic mass for NS maximum mass (Msol)
 * @param       p_SNenginePrescription          Fryer2012 supernova mechanism (DELAYED, RAPID)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
DBL_DBL GiantBranch::CalculateRemnantMass_Fryer2012(const double p_Mass, const double p_COCoreMass, const double p_NSmaxBaryonicMass, const SN_ENGINE p_SNenginePrescription) const {

    double baryonicMass;
    double fallbackFraction;
    double gravitationalMass;
    double mProto;

    switch (p_SNenginePrescription) {                                       // which Fryer SN mechanism?

        case SN_ENGINE::DELAYED:                                            // DELAYED

            mProto            = CalculateProtoCoreMass_Fryer2012_Delayed(p_COCoreMass);
            fallbackFraction  = CalculateFallbackFraction_Fryer2012_Delayed(p_Mass, mProto, p_COCoreMass);
            baryonicMass      = CalculateBaryonicRemnantMass_Fryer2012(mProto, CalculateFallbackMass_Fryer2012(p_Mass, mProto, fallbackFraction));
            gravitationalMass = CalculateGravitationalRemnantMass_Fryer2012(baryonicMass, p_NSmaxBaryonicMass);
            break;

        case SN_ENGINE::RAPID:                                              // RAPID

            mProto            = FRYER_PROTO_CORE_MASS_RAPID;
            fallbackFraction  = CalculateFallbackFractionRapid(p_Mass, mProto, p_COCoreMass);
            baryonicMass      = CalculateBaryonicRemnantMass_Fryer2012(mProto, CalculateFallbackMass_Fryer2012(p_Mass, mProto, fallbackFraction));
            gravitationalMass = CalculateGravitationalRemnantMass_Fryer2012(baryonicMass, p_NSmaxBaryonicMass);
            break;
    
        default:                                                            // unexpected prescription
            // the only way this can happen is if someone added an SN_ENGINE and it isn't accounted for
            // in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.
            THROW_ERROR(ERROR::UNEXPECTED_SN_ENGINE);                       // throw error
    }

    return std::make_tuple(gravitationalMass, std::max(0.0, std::min(1.0, fallbackFraction)));
}


/*
 * CalculateRemnantMass_Fryer2022
 *
 * @brief
 * Calculate the remnant mass, per Fryer et al. 2022, eq 5.
 *
 *
 * DBL_DBL CalculateRemnantMass_Fryer2022(const double    p_Mass,
 *                                        const double    p_COCoreMass,
 *                                        const double    p_NSmaxBaryonicMass,
 *                                        const double    p_fMix,
 *                                        const double    p_mCrit
 *                                        const SN_ENGINE p_SNenginePrescription) const
 *
 * @param       p_Mass                          Pre-SN mass of the star (Msol)
 * @param       p_COCoreMass                    Pre-SN Carbon Oxygen (CO) core mass of the star (Msol)
 * @param       p_NSmaxBaryonicMass             Maximum baryonic mass for NS maximum mass (Msol)
 * @param       p_fMix                          Fryer 2022 mixing growth time
 * @param       p_mCrit                         Fryer 2022 critical mass for BH formation (Msol)        
 * @param       p_SNenginePrescription          Fryer2012 supernova mechanism (DELAYED, RAPID)
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
DBL_DBL GiantBranch::CalculateRemnantMass_Fryer2022(const double    p_Mass,
                                                    const double    p_COCoreMass,
                                                    const double    p_NSmaxBaryonicMass,
                                                    const double    p_fMix,
                                                    const double    p_mCrit
                                                    const SN_ENGINE p_SNenginePrescription) const {

    // calculate baryonic mass, clamped to total mass
    const double f = p_COCoreMass / p_fMix;
    double baryonicMass = std::min(1.2 + 0.05 * p_fMix + 0.01 * f * f + exp(p_fMix * (p_COCoreMass - p_mCrit)), p_Mass);

    // calculate fallback fraction and gravitational mass
    double mProto;
    double fallbackMass;
    double fallbackFraction;
    double gravitationalMass;

    switch (p_SNenginePrescription) {                                       // which Fryer SN mechanism?

        case SN_ENGINE::DELAYED:                                            // DELAYED  

            mProto            = CalculateProtoCoreMass_Fryer2012_Delayed(p_COCoreMass);
            fallbackFraction  = std::max(0.0, baryonicMass - mProto) / (p_Mass - mProto);
            fallbackFraction  = std::max(0.0, std::min(1.0, fallbackFraction));
            gravitationalMass = CalculateGravitationalRemnantMass_Fryer2012(baryonicMass, p_NSmaxBaryonicMass);
            break;

        case SN_ENGINE::RAPID:  

            mProto            = FRYER_PROTO_CORE_MASS_RAPID;
            fallbackFraction  = std::max(0.0, baryonicMass - mProto) / (p_Mass - mProto);
            fallbackFraction  = std::max(0.0, std::min(1.0, fallbackFraction));
            gravitationalMass = CalculateGravitationalRemnantMass_Fryer2012(baryonicMass, p_NSmaxBaryonicMass);
            break;
    
        default:                                                            // unexpected prescription
            // the only way this can happen is if someone added an SN_ENGINE and it isn't accounted for
            // in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.
            THROW_ERROR(ERROR::UNEXPECTED_SN_ENGINE);                       // throw error
    }
                                   
    return std::make_tuple(gravitationalMass, std::max(0.0, std::min(1.0, fallbackFraction)));
}





/*
 * ProcessCoreCollapseSN
 *
 * Driver function for Core Collapse Supernovae
 *
 * This function determines which prescription is used for the core collapse SN (via program options)
 *
 * The function calls prescription functions that update the following parameters:
 *      Mass, stellarType, drawnKickMagnitude, kickMagnitude
 *
 * At the end of this function we set the following parameters which are (so far) independent of the
 * ccSN prescriptions (but do depend on the parameters above):
 *      Luminosity, Radius, Temperature, supernova events: current = SN, past = CCSN
 *
 *
 * STELLAR_TYPE ProcessCoreCollapseSN(const StellarSNDetailsT& p_SNdetails)
 *
 * @return                                      The stellar type to which the star should evolve
 */
STELLAR_TYPE GiantBranch::ProcessCoreCollapseSN(const STELLAR_TYPE p_StellarType,
                                                const double p_Mass,
                                                const double p_CoreMass,
                                                const StellarSNDetailsT& p_SNdetails,
                                                const REMNANT_MASS_PRESCRIPTION p_RemnantMassPrescription) {

    StellarSNDetailsT SNdetails = p_SNdetails;      // copy given supernova details

    SNdetails.stellarTypePreSN = stellarType;  // stellar type pre compact object formation

    double mass;

    switch (p_RemnantMassPrescription) {                                                           // which remnant mass prescription?

        case REMNANT_MASS_PRESCRIPTION::BELCZYNSKI2002:                                                     // Belczynski 2002

            std::tie(mass, SNdetails.fallbackFraction) = CalculateRemnantMass_Belczynski2002(p_Mass, p_CoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::FRYER2012:                                                          // Fryer 2012
            std::tie(mass, SNdetails.fallbackFraction) = CalculateRemnantMass_Fryer2012(p_Mass, m_COCoreMass, p_NSmaxBaryonicMass, OPTIONS->FryerSupernovaEngine());
            break;

        case REMNANT_MASS_PRESCRIPTION::FRYER2022:                                                          // Fryer 2022

            std::tie(mass, SNdetails.fallbackFraction) = CalculateRemnantMass_Fryer2022(p_Mass,
                                                                                        p_COCoreMass,
                                                                                        p_NSmaxBaryonicMass,
                                                                                        OPTIONS->Fryer2022fmix(),
                                                                                        OPTIONS->Fryer2022mCrit(),
                                                                                        OPTIONS->FryerSupernovaEngine());
            break;

        case REMNANT_MASS_PRESCRIPTION::HURLEY2000:                                                         // Hurley 2000

            std::tie(mass, SNdetails.fallbackFraction) = Remnants::CalculateRemnantMass_Hurley2000_Static(p_COCoreMass);
            break;
        
        case REMNANT_MASS_PRESCRIPTION::MALTSEV2025:                                                        // Maltsev+ 2024

            std::tie(mass, SNdetails.fallbackFraction) = CalculateRemnantMass_Maltsev2025(m_COCoreMass, m_HeCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::MULLER2016:                                                         // Muller 2016

            std::tie(mass, SNdetails.fallbackFraction) = CalculateRemnantMass_Muller2016(p_Mass, m_COCoreMass);
            SNdetails.fallbackFraction = 0.0;
            break;

        case REMNANT_MASS_PRESCRIPTION::MULLERMANDEL:                                                       // Mandel & Mueller 2020

            std::tie(mass, SNdetails.fallbackFraction) = CalculateRemnantMass_MullerMandel2020(m_COCoreMass, m_HeCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020:                                                      // Schneider 2020
        case REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020ALT:                                                   // Schneider 2020, alternative

            std::tie(mass, SNdetails.fallbackFraction) = CalculateRemnantMass_Schneider2020(p_Mass,
                                                                                            p_COCoreMass,
                                                                                            p_HeCoreMassPreSN,
                                                                                            OPTIONS->MaximumNeutronStarMass(),
                                                                                            p_RemnantMassPrescription);
            break;
            
        default:                                                        // unexpected prescription
            // the only way this can happen is if the REMNANT_MASS_PRESCRIPTION passed to this function
            // is not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.

            THROW_ERROR(ERROR::UNEXPECTED_REMNANT_MASS_PRESCRIPTION);   // throw error  
    }

    // determine the stellar type to which the star should evolve (either use prescription or MAXIMUM_NS_MSS)
    switch (p_RemnantMassPrescription) {                                                           // which remnant mass prescription?

        case REMNANT_MASS_PRESCRIPTION::MULLER2016:                                                         // Muller 2016
            SNdetails.stellarTypePostSN = CalculateRemnantType_Muller2016(m_COCoreMass);
            break;

        case REMNANT_MASS_PRESCRIPTION::MULLERMANDEL:                                                       // Mandel & Mueller 2020
            SNdetails.stellarTypePostSN = mass > OPTIONS->MaximumNeutronStarMass() ? STELLAR_TYPE::BLACK_HOLE : STELLAR_TYPE::NEUTRON_STAR;
            break;

        case REMNANT_MASS_PRESCRIPTION::HURLEY2000:                                                         // Hurley 2000
            SNdetails.stellarTypePostSN = mass > 1.8 ? STELLAR_TYPE::BLACK_HOLE : STELLAR_TYPE::NEUTRON_STAR;               // Hurley et al. 2000, eq 92
            break;

        default:
            // default is ok here - unknown/unexpected remnant mass prescription already dealt with above

            if (mass > OPTIONS->MaximumNeutronStarMass()) {                                                  // collapse to BH
                m_Luminosity  = BH::CalculateLuminosityOnPhase_Static();                                            // luminosity of BH
                m_Radius      = BH::CalculateSchwarzschildRadius_Static(mass);                                    // Schwarzschild radius (not correct for rotating BH)
                m_Temperature = BaseStar::CalculateTemperatureOnPhase(m_Luminosity, m_Radius);               // temperature of BH
                SNdetails.stellarTypePostSN = STELLAR_TYPE::BLACK_HOLE;
            }
            else {                                                                                                  // collapse to NS
                m_Luminosity  = NS::CalculateLuminosityOnPhase_Static(mass, 0.0);                                 // luminosity of NS as it cools
                m_Radius      = NS::CalculateRadiusOnPhase_Static(mass);                                          // radius of NS
                m_Temperature = BaseStar::CalculateTemperatureOnPhase(m_Luminosity, m_Radius);               // temperature of NS
                SNdetails.stellarTypePostSN = STELLAR_TYPE::NEUTRON_STAR;
            }
    }

    if (p_Mass == p_CoreMass && m_HeCoreMass == m_COCoreMass) {                                     // entire star CO core?
                                                                                                    // yes - this is a USSN
        SNdetails.currentEvent = AddSNEvent(SN_EVENT::USSN, const SN_EVENT SNdetails.currentEvent); // flag ultra-stripped SN happening now, and ...
        SNdetails.pastEvent    = AddSNEvent(SN_EVENT::USSN, const SN_EVENT SNdetails.pastEvent);    // ... will be a past event
    }

    SNdetails.currentEvent = AddSNEvent(SN_EVENT::CCSN, const SN_EVENT SNdetails.currentEvent);     // flag core-collapse SN happening now, and ...
    SNdetails.pastEvent    = AddSNEvent(SN_EVENT::CCSN, const SN_EVENT SNdetails.pastEvent);        // ... will be a past event

    // return SN details struct by value - NRVO takes care of performance/efficiency
    return SNdetails;
}


/*
 * ResolveElectronCaptureSN
 *
 * Calculate the mass of the remnant and set remnant type - always a Neutron Star
 * Updates attributes of star; sets SN flags
 *
 *
 * Short hand wavy story. The core ignites ONeMg and collapses through electron
 * capture (e.g., Nomoto 1984 for thorough discussion and a summary in Nomoto 1987).
 * The explosion is likely accompanied by a low natal kick
 *
 *
 * STELLAR_TYPE ResolveElectronCaptureSN()
 *
 * @return                                      Stellar type of remnant (always STELLAR_TYPE::NEUTRON_STAR)
 */
STELLAR_TYPE GiantBranch::ResolveElectronCaptureSN() {
    m_Mass       = MECS_REM;                                                                // defined in constants.h
    m_CoreMass   = m_Mass;
    m_HeCoreMass = m_Mass;
    m_COCoreMass = m_Mass;
    m_Mass0      = m_Mass;
        
    SetSNCurrentEvent(SN_EVENT::ECSN);                                                      // electron capture SN happening now
    SetSNPastEvent(SN_EVENT::ECSN);                                                         // ... and will be a past event

    return STELLAR_TYPE::NEUTRON_STAR;
}


/*
 * ResolvePairInstabilitySN
 *
 * Calculate the mass of the remnant and set remnant type according to mass
 * Updates attributes of star; sets SN events
 *
 *
 * Short handwavy story.  The core is hot and massive enough that there is a significant
 * amount of high-energetic gamma rays which can create an electron-positron pair.
 * When a significant amounts of pairs are created the photons taken away
 * reduce the radiative pressure. The core contracts becomes hotter creating
 * more pairs. This runaway process will end in a SN that explodes the entire core without
 * leaving a remnant.
 *
 *
 * STELLAR_TYPE ResolvePairInstabilitySN()
 *
 * @return                                      Stellar type of remnant (always STELLAR_TYPE::MASSLESS_REMNANT)
 */
STELLAR_TYPE GiantBranch::ResolvePairInstabilitySN() {

    m_Luminosity  = 0.0;
    m_Radius      = 0.0;
    m_Temperature = 0.0;
    m_Mass        = 0.0;
    m_CoreMass    = 0.0;
    m_COCoreMass  = 0.0;
    m_HeCoreMass  = 0.0;

    m_SupernovaDetails.drawnKickMagnitude = 0.0;
    m_SupernovaDetails.kickMagnitude      = 0.0;
    m_SupernovaDetails.fallbackFraction   = 0.0;

    SetSNCurrentEvent(SN_EVENT::PISN);                                                                  // pair instability SN happening now
    SetSNPastEvent(SN_EVENT::PISN);                                                                     // ... and will be a past event

    return STELLAR_TYPE::MASSLESS_REMNANT;
}


/*
 * ProcessPPISN
 *
 * @brief
 * Process a Pulsational Pair Instability Supernova event.
 * 
 * Calculate the mass of the remnant and set remnant type according to mass
 * Updates attributes of star; sets SN events
 *
 *
 * STELLAR_TYPE ProcessPPISN(const StellarSNDetailsT& p_SNdetails) const
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch::ProcessPPISN(const STELLAR_TYPE       p_StellarType,
                                       const double             p_Mass,
                                       const double             p_Radius,
                                       const double             p_CoreMass,
                                       const double             p_COCoreMass,
                                       const double             p_HeCoreMass,
                                       const PPI_PRESCRIPTION   p_PPIprescription,
                                       const StellarSNDetailsT& p_SNdetails) const {

    StellarSNDetailsT SNdetails = p_SNdetails;      // copy given supernova details

    // update some SNdetails attributes - assumes PPISN
    // these may be overwritten if PISN (see Hendriks prescription below)
    SNdetails.stellarTypePreCOFormation = p_StellarType;
    SNdetails.stellarTypePreCOFormation = STELLAR_TYPE::BLACK_HOLE;     // remnant is a black hole

    SNdetails.massAtCOFormation       = p_Mass;
    SNdetails.coreMassAtCOFormation   = p_CoreMass;
    SNdetails.COCoreMassAtCOFormation = p_COCoreMass;
    SNdetails.HeCoreMassAtCOFormation = p_HeCoreMass;

    SNdetails.radiusAtCOFormation     = p_Radius;
    SNdetails.coreRadiusAtCOFormation = CalculateConvectiveCoreRadius();


    SNdetails.fallbackFraction = 1.0;                                              // fraction of mass that falls back

    SNdetails.currentEvent = SN_EVENT::PPISN;                                   // PPISN happening now, and ...
    SNdetails.pastEvent    = SN_EVENT::PPISN;                                           // ... will be a past event
        
        m_Luminosity  = BH::CalculateLuminosityOnPhase_Static();                                // black hole luminosity
        m_Radius      = BH::CalculateRadiusOnPhase_Static(m_Mass);                              // Schwarzschild radius (not correct for rotating BH)
        m_Temperature = CalculateTemperatureOnPhase(m_Luminosity, m_Radius);

    }


    double mass = p_Mass;


    double baryonicMass;
    switch (p_PPIprescription) {                                // which prescription?

        case PPI_PRESCRIPTION::FARMER:                                                      // FARMER
            // Farmer et al. 2019 (See http://dx.doi.org/10.3847/1538-4357/ab518b)

            if (p_COCoreMass < FARMER_PPISN_UPP_LIM_LIN_REGIME)) {            // CO core mass < FARMER_PPISN_UPP_LIM_LIN_REGIME threshold?
                mass = p_COCoreMass + 4.0;                                  // yes -linear relation below threshold
            }
            else if (p_COCoreMass < FARMER_PPISN_UPP_LIM_QUAD_REGIME) {     // no - CO core mass < FARMER_PPISN_UPP_LIM_QUAD_REGIME threshold?
                                                                            // yes - quadratic relation below threshold
                mass = -0.096 * p_COCoreMass * p_COCoreMass + 8.564 * p_COCoreMass + -2.07 * m_Log10Metallicity + -152.97;
            }
            else if (p_COCoreMass < FARMER_PPISN_UPP_LIM_INSTABILLITY) {     // no - CO core mass < FARMER_PPISN_UPP_LIM_INSTABILLITY threshold?
                mass = 0.0;                                                     // yes - zero BH mass
            }
            else {                                                              // otherwise
                mass = p_COCoreMass;                                // BH mass becomes CO-core mass above the PISN gap
            }

            mass = std::max(std::min(mass, p_Mass), 0.0);           // clamp remnant mass to [0.0, total mass pre PPISN]

            break;
    
        case PPI_PRESCRIPTION::HENDRIKS: {                                                      // HENDRIKS
            // Prescription from Hendriks et al. 2023 eq 6 (See (https://arxiv.org/abs/2309.09339))
            // Based on Renzo et al. 2022 (See https://iopscience.iop.org/article/10.3847/2515-5172/ac503e)
            // 
            // Suggest using --PPI-upper-limit 80.0 and --PISN-lower-limit 80.0
            //
            // Hendriks et al. 2023 eq 6 calculates the amount of the He core that's lost in pulsations

            const double deltaMPPICOShift = OPTIONS->PulsationalPairInstabilityCOCoreShiftHendriks();
            const double deltaMPPIExtraML = OPTIONS->PulsationalPairInstabilityExtraMLHendriks(); // make an option - default 0.0 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            const double PPIOnset  = p_COCoreMass - deltaMPPICOShift - 34.8;
            const double PPIOnset2 = PPIOnset * PPIOnset;
            const double deltaMPPI = std::max(((0.0006 * m_Log10Metallicity + 0.0054) * PPIOnset * PPIOnset2) - (0.0013 * PPIOnset2) + deltaMPPIExtraML, 0.0);
            
            mass = std::max(p_HeCoreMass - deltaMPPI, 0.0);                                // remnant mass should be non-negative		
            mass = mass >= 10.0 ? mass : 0.0;                                           // if the remnant mass is below 10 Msol, assume PISN

            } break;

        case PPI_PRESCRIPTION::MARCHANT: {                                                      // MARCHANT
            // Marchant et al. 2018 (See https://arxiv.org/abs/1810.13412)

            // pow() is slow - use multiplication
            const double HeCoreMass2 = p_HeCoreMass * p_HeCoreMass;
            const double HeCoreMass3 = p_HeCoreMass * HeCoreMass2;
            const double HeCoreMass4 = HeCoreMass2  * HeCoreMass2;
            const double HeCoreMass5 = HeCoreMass3  * HeCoreMass2;
            const double HeCoreMass6 = HeCoreMass3  * HeCoreMass3;
            const double HeCoreMass7 = p_HeCoreMass * HeCoreMass6;

            const double ratioOfRemnantToHeCoreMass = std::max(0.0, std::min(1.0, (-1.63057326E-08 * HeCoreMass_7) +
                                                                                  ( 5.36316755E-06 * HeCoreMass_6) +
                                                                                  (-7.52206933E-04 * HeCoreMass_5) +
                                                                                  ( 5.83107626E-02 * HeCoreMass_4) +
                                                                                  (-2.69801221E+00 * HeCoreMass_3) +
                                                                                  ( 7.45060098E+01 * HeCoreMass_2) +
                                                                                  (-1.13694590E+03 * p_HeCoreMass) +
                                                                                    7.39643451E+03));

            baryonicMass = ratioOfRemnantToHeCoreMass * p_HeCoreMass;                           // strip off the hydrogen envelope if any was left
            m_Mass       = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                  // convert to gravitational mass due to neutrino mass loss
            } break;

        case PPI_PRESCRIPTION::STARTRACK:                                                       // STATRACK
            // Belczynski et al. 2016 https://arxiv.org/abs/1607.03116
            baryonicMass = std::min(p_HeCoreMass, STARTRACK_PPISN_HE_CORE_MASS);                // strip off the hydrogen envelope if any was left, limit to STARTRACK_PPISN_HE_CORE_MASS (default 45 Msun)
            m_Mass       = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                  // convert to gravitational mass due to neutrino mass loss
            break;

        case PPI_PRESCRIPTION::WOOSLEY:                                                         // Woosley 2017 https://arxiv.org/abs/1608.08939
            baryonicMass = p_HeCoreMass;                                                        // strip off the hydrogen envelope if any was left
            m_Mass       = BH::CalculateNeutrinoMassLoss_Static(baryonicMass);                  // convert to gravitational mass due to neutrino mass loss
            break;

        default:                                                                                // unknown prescription
            // the only way this can happen is if someone added a REMNANT_MASS_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_PPI_PRESCRIPTION);                                       // throw error
    }

    if (utils::Compare(m_Mass, 0.0) <= 0) {                                                     // remnant mass <= 0?
        stellarType = ResolvePairInstabilitySN();                                               // yes - PISN rather than PPISN
    }
    else {                                                                                      // no - PPISN
        SetSNCurrentEvent(SN_EVENT::PPISN);                                                     // pulsational pair instability SN happening now
        SetSNPastEvent(SN_EVENT::PPISN);                                                        // ... and will be a past event

        stellarType   = STELLAR_TYPE::BLACK_HOLE;                                               // -> black hole
        
        m_Luminosity  = BH::CalculateLuminosityOnPhase_Static();                                // black hole luminosity
        m_Radius      = BH::CalculateRadiusOnPhase_Static(m_Mass);                              // Schwarzschild radius (not correct for rotating BH)
        m_Temperature = CalculateTemperatureOnPhase(m_Luminosity, m_Radius);

        m_SupernovaDetails.fallbackFraction = 1.0;                                              // fraction of mass that falls back
    }

    return stellarType;
}





/*
 * ProcessSupernova
 *
 * This function checks if the star is about to undergo a SN event, and if so determines
 * the type of the supernova and calls the appropriate functions to calculate attributes
 * correctly, and to determine the type of remnant to which the star should evolve.
 * 
 * The WhiteDwarfs class also has a ProcessSupernova() function that handles WD SNe.
 *
 * DRIVER, so USES state, globals, options etc.  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 * 
 * 
 * STELLAR_TYPE ProcessSupernova()
 *
 * @return                                      Stellar type of remnant
 */
StellarSNDetailsT GiantBranch::ProcessSupernova(const bool p_AllowECSN, const StellarSNDetailsT& p_SNdetails) {
    // electron capture SN; requires progenitor to have been a MT donor unless non-stripped ECSN are allowed
    // parm should be TRUE for SSE if OPTIONS->AllowNonStrippedECSN() is TRUE
    // parm should be true for call from constituent if !m_MassTransferDonorHistory.empty() || OPTIONS->AllowNonStrippedECSN()

    StellarSNDetailsT SNdetails = p_SNdetails;      // copy given supernova details

    // set preSN values
    SNdetails.COCoreMassPreSN  = m_StateHistory.CurrentState.COCoreMass();
    SNdetails.coreMassPreSN    = m_StateHistory.CurrentState.CoreMass();
    SNdetails.coreRadiusPreSN  = CalculateConvectiveCoreRadius();
    SNdetails.HeCoreMassPreSN  = m_StateHistory.CurrentState.HeCoreMass();
    SNdetails.massPreSN        = m_StateHistory.CurrentState.Mass();
    SNdetails.radiusPreSN      = m_StateHistory.CurrentState.Radius();
    SNdetails.stellarTypePreSN = m_StateHistory.CurrentState.StellarType();

    SetSNHydrogenContent();





                                                                 // store pre-SN stellar attributes; set H content

    if (                OPTIONS->UsePulsationalPairInstability()        &&  // Pulsational Pair Instability enabled, and ...
        m_HeCoreMass >= OPTIONS->PulsationalPairInstabilityLowerLimit() &&  // ... He core mass at or above lower PPI threshold, and ...
        m_HeCoreMass <= OPTIONS->PulsationalPairInstabilityUpperLimit()) {  // ... He core mass at or below upper PPI threshold?
                                                                            // yes
        stellarType = ProcessPPISN();                // Pulsational Pair Instability Supernova (PPISN)
    }
    else if (                OPTIONS->UsePairInstabilitySupernovae() &&     // Pair Instability enabled, and ...
             m_HeCoreMass >= OPTIONS->PairInstabilityLowerLimit()    &&     // ... He core mass at or above lower PI threshold, and ...
             m_HeCoreMass <= OPTIONS->PairInstabilityUpperLimit()) {        // ... He core mass at or below upper PI threshold?
                                                                            // yes
        stellarType = ResolvePairInstabilitySN();                           // Pair Instability Supernova (PISN)
    }
    else if (CalculateInitialSupernovaMass() < MCBUR2 &&                    // initial SN mass below non-degenerate threshold, and ...
             p_AllowECSN) {                                                 // ... ECSN allowed (see notes above)?
                                                                            // yes
        stellarType = ResolveElectronCaptureSN();                           // Electron Capture Supernova (ECSN)
    }
    else {                                                                  // otherwise
        stellarType = ResolveCoreCollapseSN(p_SNdetails);                              // Core Collapse Supernova (CCSN)
    }
        
        // check if the SN will actually happen (i.e. new stellar type indicates a switch)
    if (utils::IsOneOf(stellarType, SN_REMNANTS)) {
        if (utils::SNEventType(m_SupernovaDetails.events.current) != SN_EVENT::PISN &&          // NOT a PISN event, and ...
           !utils::IsOneOf(stellarType, {STELLAR_TYPE::MASSLESS_REMNANT})) {           // NOT a massless remnant?
            CalculateSNkickMagnitude(m_SupernovaDetails.totalMassAtCOFormation - p_Mass, m_Mass, stellarType);  // yes - calculate kick magnitude
        }
    }
    else {
            // THROW ERROR HERE - NO SN HAPPENED!!!
    }



SNdetails.stellarTypePostCOFormation = stellarType;

    // return SN details struct by value - NRVO takes care of performance/efficiency
    return SNdetails;
}





/*
 * CalculateConvectiveEnvelopeMass
 *
 * @brief
 * Calculate the mass of the convective envelope.
 * 
 * Approximates the mass of the outer convective envelope.
 * Follows the fits of Picker, Hirai, Mandel 2024 (see https://arxiv.org/pdf/2402.13180)
 * 
 * Mandel, Hirai, Picker 2024 (see https://arxiv.org/pdf/2412.10691)
 *
 *
 * DBL_DBL CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_CoreMass) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_MassEffectiveInitial          Effective initial mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   GB convective envelope mass (Msol)
 *                                                   Maximum GB convective envelope mass (Msol)
 */
COMPAS_PURE DBL_DBL GiantBranch::CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_MassEffectiveInitial) const {
    
    // ratio of final intershell mass to final core mass - Picker at al. 2024 eq 8
    double MiFinalMcFinal = -0.023 * GLOBALS->SigmaHurley() - 0.0023;

    // We need the temperature of the star just after BAGB, which is the temperature at the
    // start of the EAGB phase.  Since we are on the giant branch here, we can clone this
    // object as an EAGB object and, as long as it is initialised (to the start of the phase),
    // we can query the cloned object for its temperature.
    //
    // To ensure the clone does not participate in logging, we set its persistence to EPHEMERAL.
    //
    // Furthermore, 'this' is const in this function, so we first remove its const-ness (required
    // to call Clone()) via the use of const_cast<>().  Since we don't know what class the
    // underlying object is, we cast it to EAGB&.

    EAGB *clone = EAGB::Clone(static_cast<EAGB&>(const_cast<GiantBranch&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    clone->EvolveOneTimestep(0.0, 0.0, 0.0, true);                                                                          // otherwise, temperature not updated
    double tMin = clone->Temperature();                                                                                     // get temperature of clone
    delete clone; clone = nullptr;                                                                                          // return the memory allocated for the clone
    
    // Use Eq. 6 of Mandel et al. 2024 eq 6 rather than Picker at al. 2024 eq 6 for tOnset
    // to avoid issues caused by differences between temperatures in MESA models (used in
    // fits from Picker et al. 2024) and Pols models (used in Hurley SSE tracks).
    double tOnset  = tMin / std::min(0.695 - 0.057 * GLOBALS->SigmaHurley(), 0.95);                                        // Mandel et al. 2024, eq 6
    double McFinal = CalculateCoreMassAtBAGB_Hurley2000(p_MassEffectiveInitial);
    double mEnvMax = std::max(p_Mass - mCoreFinal * (1.0 + MiFinalMcFinal), 0.0);                                         // Picker at al. 2024, eq 9

    // Picker+ 2024 fits were only made for stars above 8.0 solar masses, with runs down to 5.0 solar masses, 
    // so using the final core mass as an approximate threshold of validity
    if(utils::Compare(mCoreFinal, 1.5) < 0) mEnvMax = std::max(p_Mass - mCoreFinal, 0.0);                                  // unlike massive stars, intermediate-mass stars have almost no radiative intershell at maximum convective envelope extent
    
    double mEnv = mConvMax / (1.0 + exp(4.6 * (tMin + tOnset - 2.0 * m_Temperature) / (tMin - tOnset)));  // Picker at al. 2024, eq 7
    
    return std::make_tuple(mEnv, mEnvMax);
}









///////////////////////////////////////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



/*
 * ResolveSupernova
 *
 * This function checks if the star is about to undergo a SN event, and if so determines
 * the type of the supernova and calls the appropriate functions to calculate attributes
 * correctly, and to determine the type of remnant to which the star should evolve.
 * 
 * The WhiteDwarfs class also has a ResolveSupernova() function that handles WD SNe.
 *
 *
 * STELLAR_TYPE ResolveSupernova()
 *
 * @return                                      Stellar type of remnant
 */
STELLAR_TYPE GiantBranch_Constituent::ResolveSupernova(pass MT history here - or flag?) {

    STELLAR_TYPE stellarType = GiantBranch::ResolveSupernova(!m_MassTransferDonorHistory.empty());

    // only NSs can get rocket kicks
    if (!utils::IsOneOf(stellarType, { STELLAR_TYPE::NEUTRON_STAR })) m_SupernovaDetails.rocketKickMagnitude = 0.0;
}



/*
 * CalculateRemnantMass_Schneider2020
 *
 * @brief
 * Calculate remnant mass for binary stars, per Schneider et al. 2020.
 * (See arxiv:2008.08599).
 * 
 * Function provides for standard and alternative fits (single stars) from Schneider et al. 2020.
 *
 * Remnant mass is based on the donor type at the latest (most recent) MT event.
 * 
 * Note that Schneider only prescribes remnant masses for the simple cases of single episode mass
 * transfer, so some of the double episode cases here are a bit uncertain, and may need to be refined
 * at a later date.
 * 
 *
 * DBL_DBL CalculateRemnantMass_Schneider2020(const double                    p_Mass,
 *                                            const double                    p_COCoreMass,
 *                                            const double                    p_HeCoreMassPreSN,
 *                                            const double                    p_MaxNSMass,
 *                                            const MT_CASE                   p_MTcase,
 *                                            const REMNANT_MASS_PRESCRIPTION p_RemnantMassPrescription) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMassPreSN               Pre-SN He core mass of the star (Msol)
 * @param       p_MaxNSMass                     User-specified maximum mass of a neutron star (`--maximum-neutron-star-mass`) (Msol)
 * @param       p_MTcase                        Mass transfer case of latest MT event (MT_CASE::NONE if no MT history) // m_StateHistory.CurrentState.LatestMTdonorEvent().MTcase <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 * @param       p_RemnantMassPrescription       Remnant mass prescription to use
 * @return                                      Tuple containing:
 *                                                   DOUBLE Remnant mass (Msol)
 *                                                   DOUBLE Fraction of mass falling back onto compact object [0.0, 1.0]
 */
DBL_DBL GiantBranch_Constituent::CalculateRemnantMass_Schneider2020(const double                    p_Mass,
                                                                    const double                    p_COCoreMass,
                                                                    const double                    p_HeCoreMassPreSN,
                                                                    const double                    p_MaxNSMass,
                                                                    const MT_CASE                   p_MTcase,
                                                                    const REMNANT_MASS_PRESCRIPTION p_RemnantMassPrescription) const {
    double remnantMass;
    double fallbackFraction;

    // if no history of MT, or if MT happens from naked He stars, WDs, etc.,
    // assume that the core properties are not affected
    if (MTcase == MT_CASE::NONE || MTcase == MT_CASE::OTHER) {
        std::tie(remnantMass, fallbackFraction) = GiantBranch::CalculateRemnantMass_Schneider2020(p_Mass,
                                                                                                  p_COCoreMass,
                                                                                                  p_HeCoreMassPreSN,
                                                                                                  p_MaxNSMass,
                                                                                                  p_RemnantMassPrescription);
    }
    else {                                                                          // handle MT cases
        switch (MTcase) {                                                           // which MT_CASE?

            double logRemnantMass;
            case MT_CASE::A:                                                        // case A MT

                     if (p_COCoreMass <  7.064) { logRemnantMass = log10(0.02128 * p_COCoreMass + 1.35349); }
                else if (p_COCoreMass <  8.615) { logRemnantMass = 0.03866 * p_COCoreMass + 0.64417; }
                else if (p_COCoreMass < 15.187) { logRemnantMass = log10(0.02128 * p_COCoreMass + 1.35349); }
                else                            { logRemnantMass = 0.02573 * p_COCoreMass + 0.79027; }
                break;

            case MT_CASE::B:                                                        // case B MT

                     if (p_COCoreMass <  7.548) { logRemnantMass = log10(0.01909 * p_COCoreMass + 1.34529); }
                else if (p_COCoreMass <  8.491) { logRemnantMass = 0.03306 * p_COCoreMass + 0.68978; }
                else if (p_COCoreMass < 15.144) { logRemnantMass = log10(0.01909 * p_COCoreMass + 1.34529); }
                else                            { logRemnantMass = 0.02477 * p_COCoreMass + 0.80614; }
                break;

            case MT_CASE::C:                                                        // case C MT

                     if (p_COCoreMass <  6.357) { logRemnantMass = log10(0.03781 * p_COCoreMass + 1.36363); }
                else if (p_COCoreMass <  7.311) { logRemnantMass = 0.05264 * p_COCoreMass + 0.58531; }
                else if (p_COCoreMass < 14.008) { logRemnantMass = log10(0.03781 * p_COCoreMass + 1.36363); }
                else                            { logRemnantMass = 0.01795 * p_COCoreMass + 0.98797; }
                break;

            default:                                                              // unexpected MT_CASE
                // the only way this can happen is if the MT_CASE passed to this function is not accounted
                // for in this code.  We should not default here, with or without a warning.
                // We are here because the code passed a mass transfer case that this function doesn't account
                // for, and that should be flagged as an error and result in termination of the evolution of
                // binary.
                // The correct fix for this is to add code to this function for the missing mass transfer case,
                // or fix the calling code to pass a mass transfer case that is handled by this function.

                THROW_ERROR(ERROR::UNEXPECTED_MT_CASE);                             // throw error    
        }

        // convert log mass to linear and clamp to pre-SN He core mass of the star
        remnantMass = std::min(PPOW(10.0, logRemnantMass), p_HeCoreMassPreSN);
    }

    // fallback fraction of mass beyond proto-neutron-star for BH formation and kicks
    const double fallbackFraction = remnantMass > p_MaxNSMass ? (remnantMass - NEUTRON_STAR_MASS) / (p_Mass - NEUTRON_STAR_MASS) : 0.0;		

    return std::make_tuple(remnantMass, std::max(0.0, std::min(1.0, fallbackFraction)));
}



/*
 * CalculateRemnantMass_Maltsev2025  *** JR: should we rename this? *Ilya*
 *
 * Calculate remnant mass for binary stars, per Maltsev et al. 2024.
 * (See https://doi.org/10.1051/0004-6361/202554931)
 *
 * Remnant mass is based on the donor type at the first MT event.
 * 
 *
 * DBL_DBL CalculateRemnantMass_Maltsev2025(const double       p_COCoreMass,
 *                                          const double       p_HeCoreMass,
 *                                          const double       p_ZetaAsplund, 
 *                                          const double       p_MaltsevFallbackFraction,
 *                                          const MT_CASE      p_MTcase,
 *                                          const MALTSEV_MODE p_MaltsevMode) const
 *
 * @param       p_COCoreMass                    CO core mass of the star (Msol)
 * @param       p_HeCoreMass                    He core mass of the star (Msol)
 * @param       p_ZetaAsplund                   Asplund zeta value (log10(Z / ZSOL_ASPLUND))
 * @param       p_MaltsevFallbackFraction       User-specfifid Maltsev fallback fraction (`--maltsev-fallback-fraction`) [0.0, 1.0]
 * @param       p_MTcase                        Mass transfer case of first MT event (MT_CASE::NONE if no MT history) // m_StateHistory.CurrentState.FirstMTdonorEvent().MTcase <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 * @param       p_MaltsevMode                   User-specified Maltsev mode (`--maltsev-mode`))
 * @return                                      Remnant mass in Msol
 */
DBL_DBL GiantBranch_Constituent::CalculateRemnantMass_Maltsev2025(const double       p_COCoreMass,
                                                                  const double       p_HeCoreMass,
                                                                  const double       p_ZetaAsplund, 
                                                                  const double       p_MaltsevFallbackFraction, //OPTIONS->MaltsevFallbackFraction()  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                                                  const MT_CASE      p_MTcase,
                                                                  const MALTSEV_MODE p_MaltsevMode) const { // OPTIONS->MaltsevMode() <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    double fallbackFraction;
    double remnantMass;

    if (p_COCoreMass < MALTSEV2025_MMIN) {                                                  // CO core mass < Maltsev lower threshold?
                                                                                            // yes - NS formation regardless of metallicity
        fallbackFraction = 0.0;
        remnantMass = NEUTRON_STAR_MASS;
    }
    else if (p_COCoreMass > MALTSEV2025_MMAX) {                                             // no - CO core mass > Maltsev upper threshold?
                                                                                            // yes - BH formation regardless of metallicity
        fallbackFraction = 1.0;
        remnantMass = p_HeCoreMass;
    }
    else {                                                                                  // no - calculate remnant mass

        double Zbound;                                                                      
        switch (p_MaltsevMode) {                                                            // which Maltsev mode?                                                                                 

            case MALTSEV_MODE::OPTIMISTIC:                                                  // OPTIMISTIC               
                Zbound = p_ZetaAsplund;                                                     // Asplund zeta (log10(Z / ZSOL_ASPLUND))
                break;

            case MALTSEV_MODE::BALANCED:                                                    // BALANCED                        
                Zbound = std::min(std::max(p_ZetaAsplund, -1.6989700043360185), 0.0);       // Asplund zeta clamped to [log10(1/50), log10(1)]
                break;

            case MALTSEV_MODE::PESSIMISTIC:                                                 // Pessimistic                          
                Zbound = std::min(std::max(p_ZetaAsplund, -1.0), 0.0);                      // Asplund zeta clamped to [log10(1/10), log10(1)]
                break;

            default:                                                                        // unexpected mode
                // the only way this can happen is if the MALTSEV_MODE passed to this function
                // isn't accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the code passed a mode that this function doesn't account for, and
                // that should be flagged as an error and result in termination of the evolution of the star
                // or binary.
                // The correct fix for this is to add code to this function for the missing mode, or fix the
                // calling code to pass a mode that is handled by this function.
                THROW_ERROR(ERROR::UNEXPECTED_MALTSEV_MODE);                                // throw error
        }
       
        MT_CASE MTcase = p_MTcase;                                                          // copy given mass transfer case

        // if the star was stripped by winds, treat it as if it experienced Case B mass transfer
        if (MTcase == MT_CASE::NONE && HydrogenAbundanceSurface() == 0.0) MTcase = MT_CASE::B;

        // if no history of MT, or if MT happens from naked He stars, WDs, etc.,
        // assume that the core properties are not affected
        if (MTcase == MT_CASE::NONE || MTcase == MT_CASE::OTHER) {
            std::tie(remnantMass, fallbackFraction) = GiantBranch::CalculateRemnantMass_Maltsev2025(const double       p_COCoreMass,
                                                                                                    const double       p_HeCoreMass,
                                                                                                    const double       p_ZetaAsplund, 
                                                                                                    const double       p_MaltsevFallbackFraction,
                                                                                                    const MALTSEV_MODE p_MaltsevMode);
        }
        else {                                                                              // handle MT cases
            switch (massTransferCase) {                                                     // which MT_CASE?

                case MT_CASE::A:                                                            // case A MT
                    M1 = MALTSEV2025_M1A + (MALTSEV2025_M1A - MALTSEV2025_M1AZ01) * Zbound;
                    M2 = MALTSEV2025_M2A + (MALTSEV2025_M2A - MALTSEV2025_M2AZ01) * Zbound;
                    M3 = MALTSEV2025_M3A + (MALTSEV2025_M3A - MALTSEV2025_M3AZ01) * Zbound;
                    break;

                case MT_CASE::B:                                                            // case B MT
                    M1 = MALTSEV2025_M1B + (MALTSEV2025_M1B - MALTSEV2025_M1BZ01) * Zbound;
                    M2 = MALTSEV2025_M2B + (MALTSEV2025_M2B - MALTSEV2025_M2BZ01) * Zbound;
                    M3 = MALTSEV2025_M3B + (MALTSEV2025_M3B - MALTSEV2025_M3BZ01) * Zbound;
                    break;

                case MT_CASE::C:                                                            // case C MT
                    M1 = MALTSEV2025_M1C + (MALTSEV2025_M1C - MALTSEV2025_M1CZ01) * Zbound;
                    M2 = MALTSEV2025_M2C + (MALTSEV2025_M2C - MALTSEV2025_M2CZ01) * Zbound;
                    M3 = MALTSEV2025_M3C + (MALTSEV2025_M3C - MALTSEV2025_M3CZ01) * Zbound;
                    break;

                default:                                                                    // unexpected MT_CASE
                    // the only way this can happen is if the MT_CASE passed to this function is not accounted
                    // for in this code.  We should not default here, with or without a warning.
                    // We are here because the code passed a mass transfer case that this function doesn't account
                    // for, and that should be flagged as an error and result in termination of the evolution of
                    // binary.
                    // The correct fix for this is to add code to this function for the missing mass transfer case,
                    // or fix the calling code to pass a mass transfer case that is handled by this function.

                    THROW_ERROR(ERROR::UNEXPECTED_MT_CASE);                                 // throw error 
            }
        

            if (p_COCoreMass >= M3 || (p_COCoreMass >= M1 && p_COCoreMass <= M2)) {         // complete fallback into BH
                fallbackFraction = 1.0;
                remnantMass = p_HeCoreMass;
            }
            else if (p_COCoreMass > M2 && p_COCoreMass < M3 && RAND->Random(0, 1) <= 0.1) { // partial fallback - BH formation
                fallbackFraction = p_MaltsevFallbackFraction;
                remnantMass = (p_HeCoreMass - NEUTRON_STAR_MASS) * fallbackFraction + NEUTRON_STAR_MASS;
            }
            else {                                                                          // no fallback - NS formation
                fallbackFraction = 0.0;
                remnantMass = NEUTRON_STAR_MASS;
            }
        }
    }

    return std::make_tuple(remnantMass, std::max(0.0, std::min(1.0, fallbackFraction)));
}






/*
 * CalculateCriticalMassRatio_Claeys2014
 *
 * @brief
 * Calculate the critical mass ratio for unstable mass transfer, per Claeys et al. 2014.
 *
 * See e.g de Mink et al. 2013, Claeys et al. 2014, and Ge et al. 2010, 2015, 2020 for discussions.
 *
 * Assumes this star is the donor.
 * 
 * Critical mass ratio is defined as qCrit = Maccretor / Mdonor.
 * 
 *
 * double CalculateCriticalMassRatio_Claeys2014(const double p_Mass, const double p_HeCoreMass, const bool p_AccretorIsDegenerate) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_HeCoreMass                    Helium core mass of the star (Msol)
 * @param       p_AccretorIsDegenerate          Boolean indicating if accretor is degenerate
 * @return                                      Critical mass ratio for unstable MT 
 */
double GiantBranch_Constituent::CalculateCriticalMassRatio_Claeys2014(const double p_Mass, const double p_HeCoreMass, const bool p_AccretorIsDegenerate) const {

    double qCrit;
                                                                                                                            
    if (p_AccretorIsDegenerate) {                                                                           // degenerate accretor?
        qCrit = OPTIONS->MassTransferCriticalMassRatioGiantDegenerateAccretor();                            // yes - use option value
    }
    else {                                                                                                  // no - non-degenerate accretor 
        if OPTIONS->OptionDefaulted("critical-mass-ratio-giant-non-degenerate-accretor") {                  // user specified qCrit value?
                                                                                                            // no - calculate it per Claeys et al. 2014
            if (p_HeCoreMass > 0.0) {
                const double mRatio = p_HeCoreMass / p_Mass;
                qCrit = 2.13 / ( 1.67 - GLOBALS->HurleyRadiusXexponent() + 2.0 * utils::IntPow(mRatio, 5)); // Claeys et al. 2014, Table 2
            }
            else qCrit = 0.0;
        }
        else {                                                                                              // yes - user sepcified value
            qCrit = OPTIONS->MassTransferCriticalMassRatioGiantNonDegenerateAccretor();                     // use the option value specified
        }
    }

    return qCrit;
}









///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 ZETA CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateZetaAdiabatic_ConvectiveEnvelopeGiant
 *
 * @brief
 * Calculate the adiabatic exponent (donor radial response to mass loss), zeta, for
 * convective-envelope giant-like stars
 *
 *
 * double CalculateZetaAdiabatic_ConvectiveEnvelopeGiant(ZETA_PRESCRIPTION p_ZetaPrescription) const
 *
 * @param       p_ZetaPrescription              User-specified (or default) ZETA_PRESCRIPTION
 * @return                                      Adiabatic exponent, zeta
 */
double GiantBranch_Constituent::CalculateZetaAdiabatic_ConvectiveEnvelopeGiant(ZETA_PRESCRIPTION p_ZetaPrescription) {
    
    double zeta;

    switch (p_ZetaPrescription) {                               // which ZETA_PRESCRIPTION?
            
        case ZETA_PRESCRIPTION::HURLEY:                         // HURLEY
            zeta = CalculateZetaAdiabatic_Hurley2002(m_CoreMass);
            break;
            
        case ZETA_PRESCRIPTION::ARBITRARY:                      // ARBITRARY
            zeta = OPTIONS->ZetaAdiabaticArbitrary();
            break;

        case ZETA_PRESCRIPTION::SOBERMAN:                       // SOBERMAN
            zeta = CalculateZetaAdiabatic_Soberman1997(m_CoreMass);
            break;

        default:                                                // unexpected prescription
            // the only way this can happen is if the ZETA_PRESCRIPTION passed to this function is
            // not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.

            THROW_ERROR(ERROR::UNEXPECTED_ZETA_PRESCRIPTION);   // throw error             
    }
    
    return zeta;
}



/*
 * CalculateZetaAdiabatic_ByEnvelopeType
 * 
 * @brief
 * Calculate the adiabatic exponent (donor radial response to mass loss), zeta,
 * based on the envelope type of the star.
 *
 * Hurley et al. 2000, eqs 97 & 98
 *
 *
 * double CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription)
 *
 * @param   [IN]    p_ZetaPrescription          ZETA_PRESCRIPTION
 * @return                                      mass-radius response exponent Zeta
 */
double GiantBranch_Constituent::CalculateZetaAdiabatic_ByEnvelopeType(ZETA_PRESCRIPTION p_ZetaPrescription) {
    
    double zeta;
    
    switch (DetermineEnvelopeType()) {                                      // which envelope type?
        case ENVELOPE::RADIATIVE:                                           // RADIATIVE
            zeta = OPTIONS->ZetaRadiativeEnvelopeGiant();                   // use option value
            break;
            
        case ENVELOPE::CONVECTIVE:                                          // CONVECTIVE
            zeta = CalculateZetaAdiabatic_ConvectiveEnvelopeGiant(p_ZetaPrescription);
            break;

        case ENVELOPE::REMNANT:                                             // REMNANT
            zeta = 0.0;                                                     // always 0.0
            break;

        default:                                                            // unknown envelope type
            // the only way this can happen is if someone added a new envelope type (to ENVELOPE)
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because DetermineEnvelopeType() returned an envelope type this code doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution of the
            // star or binary.
            // The correct fix for this is to add code for the missing envelope type or, if the missing envelope
            // type is incorrect/superfluous, remove it from ENVELOPE.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_TYPE);                      // throw error
    }
    
    return zeta;
}



