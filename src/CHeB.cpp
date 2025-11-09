#include "CHeB.h"
#include "EAGB.h"
#include "HeMS.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         AGE / LIFETIME / TAU / TIMESCALES                         //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateBluePhaseFBL_Hurley2000
 *
 * @brief
 * Calculate the Hurley blue phase parameter, fbl, per Hurley et al. 2000, just after eq 58
 *
 *
 * double CalculateBluePhaseFBL_Hurley2000(const double p_Mass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Blue phase parameter, fbl
 */
double CHeB::CalculateBluePhaseFBL_Hurley2000(const double p_Mass) {
#define b(x) GLOBALS->HurleyBCoefficients(x) // for convenience and readability - undefined at end of function

    // calculate RmHe for M > MFGB > MHeF
    const double m_b28 = PPOW(p_Mass, b(28));
    const double top   = ((bN(24) * p_Mass) + (PPOW((b(25) * p_Mass), b(26)) * m_b28)) / (b(27) + m_b28);

    // might be that we are supposed to use min(RmHe, Rx = RHeI)
    const double rHeI = CalculateRadiusAtHeI(p_Mass);
    const double lHeI = GiantBranch::CalculateLuminosityAtHeI_Static(p_Mass);

    // calculate RAGB(LHeI(M)) for M > MFGB > MHeF
    double brackets = 1.0 - (std::min(top, rHeI) / EAGB::CalculateRadiusOnPhase_Static(p_Mass, lHeI));

    return PPOW(p_Mass, b(48)) * PPOW(brackets, b(49));

#undef bN
}


/*
 * CalculateLifetimeOnBluePhase_Hurley2000
 *
 * @brief
 * Calculate tHe relative lifetime of the blue phase of Core Helium Burning, tbl,
 * per Hurley et al. 2000, eq 58
 *
 *
 * double CalculateLifetimeOnBluePhase_Hurley2000(double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      CHeB relative lifetime of blue phase, tbl [0.0, 1.0]
 */
double CHeB::CalculateLifetimeOnBluePhase_Hurley2000(const double p_Mass) const {                                          
// macros for convenience and readability - undefined at end of function
#define massCutoffs(x) GLOBALS->HurleyMassCutoffs(static_cast<int>(x))

    double tbl;

    if (p_Mass < massCutoffs(MHeF)) {
        tbl = 1.0;
    }
    else if (p_Mass <= massCutoffs(MFGB)) {
        const double mass_MFGB = p_Mass / massCutoffs(MFGB);
        const double term1     = (GLOBALS->HurleyBCoefficients(45) * PPOW(mass_MFGB, 0.414));
        tbl = term1 + ((1.0 - firstTerm) * PPOW((log10(mass_MFGB) / log10(massCutoffs(MHeF) / massCutoffs(MFGB))), GLOBALS->HurleyBCoefficients(46)));
    }
    else {
        const double fblM    = CalculateBluePhaseFBL(p_Mass);
        const double fblMFGB = CalculateBluePhaseFBL(massCutoffs(MFGB));

        tbl = (1.0 - GLOBALS->HurleyBCoefficients(47)) * (fblM / fblMFGB);
    }

    return std::min(1.0, std::max(0.0, tbl));

#undef massCutoffs
}


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
 * DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const double p_CoreMass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
COMPAS_PURE DBL_VECTOR CHeB::CalculateTimescales_Hurley2000(const double p_Mass, const double p_CoreMass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const {

    DBL_VECTOR tScales = p_tScales;                                                                         // copy given timescales

    // (re)calculate GB timescales
    tScales = GiantBranch::CalculateTimescales_Hurley2000(p_Mass, p_GBparams, tScales);

    tScales[static_cast<int>(TIMESCALE::tHe)] = CalculatePhaseLifetime_Hurley2000(p_Mass, p_CoreMass, tScales[static_cast<int>(TIMESCALE::p_tBGB)]);

	tScales[static_cast<int>(TIMESCALE::tau_BL)] = CalculateLifetimeOnBluePhase_Hurley2000(p_Mass);

    // JR: The blue loop on CHeB can be 0-length/duration (see Hurley et al., 2000, section 5.3, 
    // particularly eq 58 and beyond).  COMPAS does not allow for a 0-length blue loop - some of 
    // the equations used (e.g. to calculate Radius) result in nan or inf.  As a temporary workaround 
    // until we work out how to skip the blue loop (when it is 0-length) we will set the length of a 
    // 0-length blue loop to the absolute minimum timestep (currently ~100 seconds).
    //
    // Note that this works around a long-standing problem, which was worked around in legacy COMPAS
    // by the following code in calculateBluePhaseFBL() in star.cpp:
	//
    //    if(brackets ==0){brackets = 1e-12;}  //If zero gives R=NaN Coen Neijssel 10-01-2017
    //
    // The workaround implemented here is closer to the source of the problem (the blue loop does not
    // actually exist for some stars), and maybe a bit more meaningful (we're just using a very short 
    // duration blue loop instead of a no duration (non-existent) one)

    if (tScales[static_cast<int>(TIMESCALE::tau_BL)] <= 0.0) tScales[static_cast<int>(TIMESCALE::tau_BL)] = ABSOLUTE_MINIMUM_TIMESTEP;

    // calculate the relative age at the start of the blue phase of Core Helium Burning
    // Hurley et al. 2000, just before eq 59
    // Naturally clamped to [0, 1]
	if (p_Mass >= p_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)] && p_Mass < p_MassCutoffs[static_cast<int>(MASS_CUTOFF::MFGB)]) {
        tScales[static_cast<int>(TIMESCALE::tauX_BL)] = 1.0 - tScales[static_cast<int>(TIMESCALE::tau_BL)]; // intermediate mass stars
    }
    else {
        tScales[static_cast<int>(TIMESCALE::tauX_BL)] = 0.0;                                                // low and high mass stars
    }

    // calculate the relative age at the end of the blue phase of Core Helium Burning
    // Hurley et al. 2000, just before eq 64
    // Naturally clamped to [0, 1]
	if (p_Mass >= p_MassCutoffs[static_cast<int>(MASS_CUTOFF::MFGB)]) {
        tScales[static_cast<int>(TIMESCALE::tauY_BL)] = tScales[static_cast<int>(TIMESCALE::tau_BL)];       // high mass stars
    }
    else {
        tScales[static_cast<int>(TIMESCALE::tauY_BL)] = 1.0;                                                // intermediate and low mass stars
    }

    // return timescales vector by value - NRVO takes care of performance/efficiency
    return tScales;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////






/*
 * CalculateLuminosityAtBluePhaseStart_Hurley2000
 *
 * @brief
 * Calculate luminosity at the start of the blue phase of Core Helium Burning,
 * per Hurley et al. 2000, eq 59
 *
 *
 * double CalculateLuminosityAtBluePhaseStart_Hurley2000(const double p_Mass, const double p_CoreMass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Luminosity at the start of CHeB blue phase (Lsol)
 */
COMPAS_PURE double CHeB::CalculateLuminosityAtBluePhaseStart_Hurley2000(const double p_Mass, const double p_CoreMass) const {

    double Lx;
    if (p_Mass < GLOBALS->HurleyMassCutoffs(static_cast<int>(MHeF))) {
        Lx = GiantBranch::CalculateLuminosityOnZAHB_Hurley2000_Static(p_Mass, p_CoreMass);
    }
    else if (p_Mass < GLOBALS->HurleyMassCutoffs(static_cast<int>(MFGB))) {
        Lx = CalculateMinLuminosity_Hurley2000_Static(
            p_Mass, 
            GLOBALS->HurleyMassCutoffs(static_cast<int>(MHeF)), 
            GLOBALS->HurleyMassCutoffs(static_cast<int>(MFGB)), 
            GLOBALS->HurleyAlpha1(), 
            GLOBALS->HurleyCoefficients()
        );
    }
    else {
        Lx = GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(
            p_Mass,
            GLOBALS->HurleyMassCutoffs(static_cast<int>(MHeF)), 
            GLOBALS->HurleyAlpha1(), 
            GLOBALS->HurleyCoefficients()
        );
    }

    return Lx;
}


/*
 * CalculateLuminosityAtBluePhaseEnd_Hurley2000
 *
 * @brief
 * Calculate luminosity at the end of the blue phase of Core Helium Burning,
 * per Hurley et al. 2000, eq 61 (see discussion just before eq 64)
 *
 *
 * double CalculateLuminosityAtBluePhaseEnd_Hurley2000(const double p_Mass, const double p_CoreMass, const DBL_VECTOR& p_tScales)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      Luminosity at the end of CHeB blue phase (Lsol)
 */
double CHeB::CalculateLuminosityAtBluePhaseEnd_Hurley2000(const double p_Mass, const double p_CoreMass, const DBL_VECTOR& p_tScales) const {

    const double tx = p_tScales(tauX_BL);
    const double ty = p_tScales(tauY_BL);
    const double lx = CalculateLuminosityAtBluePhaseStart_Hurley2000(p_Mass, p_CoreMass);

    double ly;
    if (ty >= tx) {
        const double rMinHe = CalculateMinimumRadiusOnPhase_Static(
            p_Mass,
            p_CoreMass,
            GLOBALS->HurleyAlpha1(),
            GLOBALS->HurleyMassCutoffs(static_cast<int>(MHeF)),
            GLOBALS->HurleyMassCutoffs(static_cast<int>(MFGB)),
            GLOBALS->MinLuminosity_CHeB(),
            GLOBALS->HurleyBcoefficients()
        );

        const double epsilon = std::min(2.5, std::max(0.4, rMinHe / CalculateRadiusAtBluePhaseStart(p_Mass, pCoreMass)));
        const double lambda  = (ty == tx || tx == 1.0) ? 0.0 : PPOW((ty - tx) / (1.0 - tx), epsilon);   // tx can be 1.0 here - if so, lambda = 0.0

        ly = lx * PPOW(CalculateLuminosityAtBAGB(p_Mass) / lx, lambda);
    }
    else {
        const double lHeI = GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(
            p_Mass,
            GLOBALS->HurleyMassCutoffs(static_cast<int>(MHeF)),
            GLOBALS->HurleyAlpha1(),
            GLOBALS->HurleyBcoefficients()
        );

        ly = lx * PPOW(lHeI / lx, utils::intPow((tx - ty) / tx, 3));                                    // tx cannot be 0.0 here - so safe (tx > ty, ty = [0, 1])
    }

    return ly;
}


/*
 * CalculateLuminosity_Hurley2000
 *
 * @brief
 * Calculate luminosity during Core Helium Burning, per Hurley et al. 2000, eqs 61, 62 & 63
 *
 *
 * double CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Tau, const double p_CoreMass) const
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       CHeB relative age
 * @return                                      Luminosity during Core Helium Burning in Lsol
 */
double CHeB::CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Tau, const double p_CoreMass) const {

    double lCHeB;

    double tx = timescales(tauX_BL);                                                                                                        // 0 for LM and HM stars, non-zero for IM stars
    double Lx = CalculateLuminosityAtBluePhaseStart(p_Mass);

    if (utils::Compare(p_Tau, tx) >= 0) {                                                                                                   // on the blue loop
        double Rx      = CalculateRadiusAtBluePhaseStart(p_Mass);
        double RmHe    = CalculateMinimumRadiusOnPhase_Static(p_Mass, p_CoreMass, GLOBALS->HurleyAlpha1(), massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
        double LBAGB   = CalculateLuminosityAtBAGB(p_Mass);

        // the following check for high mass stars was added to match the Hurley sse code
        // - see Hurley sse `hrdiag.f` line 297
        if (p_Mass > HIGH_MASS_THRESHOLD) {
            Rx = RmHe;
        }

        double epsilon = std::min(2.5, std::max(0.4, (RmHe / Rx)));
        double lambda  = (utils::Compare(p_Tau, tx) == 0) ? 0.0 : PPOW((p_Tau - tx) / (1.0 - tx), epsilon);                                 // tx can be 1.0 here - if so, lambda = 0.0
        lCHeB          = Lx * PPOW(LBAGB / Lx, lambda);
    }
    else {                                                                                                                                  // before the blue loop
        CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass, p_MHeF, m_Alpha1, p_bN)
        double LHeI        = GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass, m_Alpha1, massCutoffs(MHeF), m_BnCoefficients);    // pow() is slow - use multiplication
        double tmp         = (tx - p_Tau) / tx;                                                                                             // tx cannot be 0.0 here, so safe (tx > tau, tau = [0, 1])
        double lambdaPrime = tmp * tmp * tmp;
        lCHeB              = Lx * PPOW((LHeI / Lx), lambdaPrime);
    }

    return lCHeB;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateMinRadius_Hurley2000_Static
 *
 * @brief
 * Calculate the minimum radius during Core Helium Burning (on the blue loop),
 * per Hurley et al. 2000, eq 55
 *
 *
 * static double CalculateMinRadius_Hurley2000_Static(const double      p_Mass,
 *                                                           const double      p_CoreMass,
 *                                                           const double      p_Alpha1,
 *                                                           const double      p_MHeF,
 *                                                           const double      p_MFGB,
 *                                                           const double      p_MinLuminosity,
 *                                                           const DBL_VECTOR &p_bN)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core Mass of the star (Msol)
 * @param       p_Alpha1                        Hurley alpha1 constant
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_MFGB                          Maximum initial mass at helium ignition on the FGB (Hurley masscutoffs[MFGB]) (Msol)
 * @param       p_MinLuminosity                 Minimum luminosity during Core Helium Burning (Lsol)
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      Minimum radius during Core Helium Burning (Rsol)
 */
double CHeB::CalculateMinRadius_Hurley2000_Static(const double      p_Mass,
                                                         const double      p_CoreMass,
                                                         const double      p_Alpha1,
                                                         const double      p_MHeF,
                                                         const double      p_MFGB,
                                                         const double      p_MinLuminosity,
                                                         const DBL_VECTOR& p_bN) {
    double radius;

    const DBL_VECTOR b = GLOBALS->HurleyBcoefficients();        // get Hurley b coefficients
    const double mHeF  = p_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)];

    if (mHeF < p_Mass) {
        const double Mb28 = PPOW(p_Mass, b[28]);                // pow() is slow - do it once only
        radius = ((b[24] * p_Mass) + (PPOW((b[25] * p_Mass), b[26]) * Mb28)) / (b[27] + Mb28);
    }
    else {
        const double lZAHB_mHeF = GiantBranch::CalculateLuminosityOnZAHB_Hurley2000_Static(mHeF, p_CoreMass);
        const double lZAHB      = GiantBranch::CalculateLuminosityOnZAHB_Hurley2000_Static(p_Mass, p_CoreMass);
        const double mHeF_b28   = PPOW(mHeF, b[28]);            // pow() is slow - do it once only
        const double top        = ((b[24] * mHeF) + (PPOW((b[25] * mHeF), b[26]) * mHeF_b28)) / (b[27] + mHeF_b28);
        const double bottom     = GiantBranch::CalculateRadiusOnPhase_Static(mHeF, lZAHB_mHeF);

        radius = GiantBranch::CalculateRadiusOnPhase_Hurley2000_Static(p_Mass, lZAHB) * PPOW(top / bottom, p_Mass / mHeF);
    }

    return radius;
}


/*
 * Calculate the radius at the start of the blue phase of Core Helium Burning
 *
 * Hurley et al. 2000, eq 60
 *
 *
 * double CalculateRadiusAtBluePhaseStart(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at the start of the blue phase of Core Helium Burning in Rsol
 */
double CHeB::CalculateRadiusAtBluePhaseStart(const double p_Mass) const {

    double Rx;

    if (utils::Compare(p_Mass, massCutoffs(MHeF)) < 0) {
        Rx = GiantBranch::CalculateRadiusOnZAHB_Static(p_Mass, m_CoreMass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
    }
    else if (utils::Compare(p_Mass, massCutoffs(MFGB)) < 0) {
        double luminosity = CalculateMinimumLuminosityOnPhase(p_Mass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_BnCoefficients);
        Rx = GiantBranch::CalculateRadiusOnPhase(p_Mass, luminosity);
    }
    else {
        Rx = CalculateRadiusAtHeIgnition(p_Mass);
    }

    return Rx;
}


/*
 * Calculate the radius at the end of the blue phase of Core Helium Burning
 *
 * Hurley et al. 2000, just before eq 64
 *
 *
 * double CalculateRadiusAtBluePhaseEnd(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @return                                      Radius at the end of the blue phase of Core Helium Burning in Rsol
 */
double CHeB::CalculateRadiusAtBluePhaseEnd(const double p_Mass) const {
    return EAGB::CalculateRadiusOnPhase_Static(p_Mass, CalculateLuminosityAtBluePhaseEnd(m_Mass0), massCutoffs(MHeF), m_BnCoefficients);
}


/*
 * Calculate the parameter Rho for radius during Core Helium Burning
 *
 * Hurley et al. 2000, eq 65
 *
 *
 * double CalculateRadiusRho(const double p_Mass, const double p_Tau)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Tau                       CHeB relative age
 * @return                                      Rho
 */
double CHeB::CalculateRadiusRho(const double p_Mass, const double p_Tau) const {

    double tx    = timescales(tauX_BL);
    double ty    = timescales(tauY_BL);

    double Rx    = CalculateRadiusAtBluePhaseStart(p_Mass);
    double Ry    = CalculateRadiusAtBluePhaseEnd(p_Mass);
    double RmHe  = CalculateMinimumRadiusOnPhase_Static(p_Mass, m_CoreMass, m_Alpha1, massCutoffs(MHeF), massCutoffs(MFGB), m_MinimumLuminosityOnPhase, m_BnCoefficients);
    double Rmin  = std::min(RmHe, Rx);

    double ty_tx = ty - tx;

    double one   = std::cbrt(log((Ry / Rmin)));
    double two   = (p_Tau - tx) / ty_tx;
    double three = std::cbrt(log((Rx / Rmin)));
    double four  = (ty - p_Tau) / ty_tx;

    return (one * two) - (three * four);
}


/*
 * CalculateRadiusOnPhase_Hurley2000
 *
 * @brief
 * Calculate the radius during Core Helium Burning, per Hurley et al. 2000, eq 64
 *
 *
 * double CalculateRadiusOnPhase_Hurley2000(const double      p_Mass,
 *                                          const double      p_Luminosity,
 *                                          const double      p_Tau,
 *                                          const double      p_CoreMass,
 *                                          const DBL_VECTOR& p_MassCutoffs,
 *                                          const DBL_VECTOR& p_Timescales,
 *                                          const DBL_VECTOR& p_aN) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MassCutoffs                   Hurley mass cutoffs (Msol)
 * @param       p_Timescales                    Phase timescales
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      Radius of the star (Rsol)
 */
double CHeB::CalculateRadiusOnPhase_Hurley2000(const double      p_Mass,
                                               const double      p_Luminosity,
                                               const double      p_Tau,
                                               const double      p_CoreMass,
                                               const DBL_VECTOR& p_MassCutoffs,
                                               const DBL_VECTOR& p_Timescales,
                                               const DBL_VECTOR& p_bN) const {
    double radius;

    if (p_Tau < p_Timescales[static_cast<int>(TIMESCALE::tauX_BL)]) {
        radius = GiantBranch::CalculateRadiusOnPhase_Hurley2000_Static(p_Mass, p_Luminosity, p_bN);
    }
    else if (p_Tau > p_Timescales[static_cast<int>(TIMESCALE::tauY_BL)]) {
        radius = EAGB::CalculateRadiusOnPhase_Hurley2000_Static(p_Mass, p_Luminosity, p_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)], p_bN);
    }
    else  {
        const double RmHe = CalculateMinRadiusOnPhase_Hurley2000_Static(p_Mass, p_CoreMass, m_Alpha1, p_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)], p_MassCutoffs[static_cast<int>(MASS_CUTOFF::MFGB)], m_MinimumLuminosityOnPhase, p_bN);
        const double rho  = std::abs(CalculateRadiusRho(p_Mass, p_Tau));

        radius = std::min(RmHe, CalculateRadiusAtBluePhaseStart(p_Mass)) * exp(rho * rho * rho);
    }

    return radius;
}












///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////







/*
 * Determine whether star should continue to evolve on phase
 *
 *
 * bool            ShouldEvolveOnPhase()
 *
 * @return         true if evolution should continue on phase, false otherwise
 */
bool CHeB::ShouldEvolveOnPhase() const {

    bool afterHeIgnition      = m_Age >= timescales(tHeI);
    bool beforeEndOfHeBurning = m_Age < (timescales(tHeI) + timescales(tHe));
    bool coreIsNotTooMassive  = utils::Compare(m_HeCoreMass, m_Mass) < 0;

    // Evolve on CHeB phase if age after He Ign and while He Burning and He core mass does not exceed total mass (could happen due to mass loss)
    return (afterHeIgnition && beforeEndOfHeBurning && coreIsNotTooMassive && !ShouldEnvelopeBeExpelledByPulsations());
}




/*
 * Determine the star's envelope type.
 *
 *
 *
 * ENVELOPE DetermineEnvelopeType()
 *
 * @return                                      ENVELOPE::{ RADIATIVE, CONVECTIVE, REMNANT }
 */
ENVELOPE CHeB::DetermineEnvelopeType() const {
    
    ENVELOPE envelope = ENVELOPE::CONVECTIVE;                                                       // default envelope type
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                                                 // which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:
            envelope = ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:
            // envelope is radiative if temperature exceeds fixed threshold, otherwise convective
            envelope = utils::Compare(Temperature() * TSOL, OPTIONS->ConvectiveEnvelopeTemperatureThreshold()) > 0 ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;
            break;

        case ENVELOPE_STATE_PRESCRIPTION::CONVECTIVE_MASS_FRACTION:
            // envelope is labeled convective when the convective mass exceeds a fixed fraction of the envelope mass
            double convectiveEnvelopeMass, convectiveEnvelopeMassMax;
            std::tie(convectiveEnvelopeMass, convectiveEnvelopeMassMax) = CalculateConvectiveEnvelopeMass();
            envelope = utils::Compare(convectiveEnvelopeMass / (m_Mass - m_CoreMass), OPTIONS->ConvectiveEnvelopeMassThreshold()) > 0 ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        default:                                                                                    // unknown prescription
            // the only way this can happen is if someone added an ENVELOPE_STATE_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION);                                // throw error             
    }
    
    return envelope;
}


/*
 * Choose timestep for evolution
 *
 * Given in the discussion in Hurley et al. 2000
 *
 *
 * ChooseTimestep(const double p_Time)
 *
 * @param   [IN]    p_Time                      Current age of star in Myr
 * @return                                      Suggested timestep (dt)
 */
double CHeB::ChooseTimestep(const double p_Time) const {

    double dtk = 2.0E-3 * timescales(tHe);
    double dte = timescales(tHeI) + timescales(tHe) - p_Time;

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);
}


/*
 * Modify the star after it loses its envelope
 *
 * Hurley et al. 2000, section 6 just before eq 76 and after Eq. 105
 *
 * Where necessary updates attributes of star (depending upon stellar type):
 *
 *     - m_StellarType
 *     - m_Timescales
 *     - m_GBparams
 *     - m_Luminosity
 *     - m_Radius
 *     - m_Mass
 *     - m_Mass0
 *     - m_CoreMass
 *     - m_HeCoreMass
 *     - m_COCoreMass
 *     - m_Age
 *
 *
 * STELLAR_TYPE ResolveEnvelopeLoss(bool p_Force)
 *
 * @param   [IN]    p_Force                     Boolean to indicate whether the resolution of the loss of the envelope should be performed
 *                                              without checking the precondition(s).
 *                                              Default is false.
 *
 * @return                                      Stellar Type to which star should evolve after losing envelope
 */
STELLAR_TYPE CHeB::ResolveEnvelopeLoss(bool p_Force) {

    STELLAR_TYPE stellarType = m_StellarType;
    
    if (ShouldEnvelopeBeExpelledByPulsations()) m_EnvelopeJustExpelledByPulsations = true;

    if (p_Force || utils::Compare(m_CoreMass, m_Mass) >= 0 || m_EnvelopeJustExpelledByPulsations ) {    // Envelope loss

        m_Mass       = std::min(m_CoreMass, m_Mass);
        m_CoreMass   = m_Mass;
        m_Mass0      = m_Mass;
        m_HeCoreMass = m_CoreMass;
        m_COCoreMass = 0.0;

 		// set evolved time for naked helium star since already has some core mass.
		double tHeIPrime = timescales(tHeI);
		double tHePrime  = timescales(tHe);

        m_Tau = (m_Age - tHeIPrime) / tHePrime;                                                         // Hurley et al. 2000, just after eq 81
        m_Age = m_Tau * HeMS::CalculatePhaseLifetime_Hurley2000_Static(m_Mass0);    // Hurley et al. 2000, eq 76 and following discussion

        CalculateTimescales(m_Mass0, m_Timescales);
        CalculateGBparams(m_Mass0, m_GBparams);

        m_Luminosity = HeMS::CalculateLuminosityOnPhase_Static(m_Mass, m_Tau);
        m_Radius     = HeMS::CalculateRadius_Hurley2000_Static(m_Mass, m_Tau);
        stellarType  = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;                                              // will evolve to an evolved helium star
    }

    return stellarType;
}


/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE CHeB::EvolveToNextPhase() {
    return STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH;
}







////////////////////////////////////////////////// CHeB_Constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<








///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 LAMBDA FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCELambda_Dewi
 *
 * @brief
 * Calculate the common envelope lambda parameter, per Dewi and Tauris, 2000
 * (https://arxiv.org/abs/astro-ph/0007034)
 * 
 * Uses the fit from Claeys et al., 2014, Appendix A (https://arxiv.org/abs/1401.2895)
 *
 * ALEJANDRO - 17/05/2017 - Not fully tested, nor fully coded, nor fully trusted.
 * Any other \lambda prescription is personally preferred.
 * Missing (A.6),(A.7),(A.8),(A.9),(A.10),(A.11) and (A.12), which have to do with ionization energy
 *
 *
 * double CalculateCELambda_Dewi(const double p_Mass,
 *                               const double p_Radius,
 *                               const double p_Luminosity,
 *                               const double p_RZAMS,
 *                               const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Common envelope lambda parameter
 */
double CHeB::CalculateCELambda_Dewi(const double p_Mass,
                                    const double p_Radius,
                                    const double p_Luminosity,
                                    const double p_RZAMS,
                                    const double p_CoreMass) const {

    const double log10Luminosity = log10(p_Luminosity);                                                 // log is expensive - do it just once

    double lambda3 = std::min(-0.9, 0.58 + (0.75 * log10(p_Mass))) - (0.08 * log10Luminosity);          // Claeys et al., 2014, A.4
    double lambda1 = std::min(lambda3, std::min(0.8, 1.25 - (0.15 * log10Luminosity)));                 // ibid., A.5, top
	double lambda2 = 0.42 * PPOW(p_RZAMS / p_Radius, 0.4);                                              // ibid., A.2
	double envMass = p_CoreMass > 0.0 && p_Mass > p_CoreMass ? p_Mass - p_CoreMass : 0.0;

    double lambdaCE;
         if (envMass >= 1.0) lambdaCE = 2.0 * lambda1;                                                  // ibid., A.1, bottom
	else if (envMass >  0.0) lambdaCE = 2.0 * (lambda2 + (std::sqrt(envMass) * (lambda1 - lambda2)));   // ibid., A.1, mid
	else                     lambdaCE = 2.0 * lambda2;                                                  // ibid., A.1, top

	return lambdaCE;
}


/*
 * CalculateCELambda_Nanjing_Enhanced
 *
 * @brief
 * Calculate the common envelope lambda parameter, per Xu & Li, 2010 (Nanjing - for the university)
 * (https://arxiv.org/abs/1004.4957, v1, 28Apr2010)
 * (https://iopscience.iop.org/article/10.1088/0004-637X/716/1/114)
 *
 * This function good for CHeB stars.
 *
 *
 * double CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
 *                                           const double             p_Radius,
 *                                           const double             p_CoreMass,
 *                                           const size_t             p_MassIndex,
 *                                           const STELLAR_POPULATION p_StellarPop) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MassIndex                     Index of mass bin in NANJING_MASSES (see constants.h)
 * @param       p_StellarPop                    Stellar population (POP I or POP II)
 * @return                                      Common envelope lambda parameter
 */
double CHeB::CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
                                                const double             p_Radius,
                                                const double             p_CoreMass,
                                                const size_t             p_MassIndex,
                                                const STELLAR_POPULATION p_StellarPop) const {

    constexpr size_t evolStage = 2;                                                         // CHeB evolutionary stage from Xu & Li, 2010

    size_t coeffsBGidx = 0;                                                                 // index into coefficients array
    bool   useLambdas  = false;                                                             // flag - use lambdas defined in the paper

    if (p_StellarPop == STELLAR_POPULATION::POPULATION_I) {                                 // pop I
        if (p_MassIndex == 1 && p_Radius > 8.5 && p_Radius < 60.0) useLambdas = true;  
        else if (p_MassIndex == 10) {
                 if (p_Radius >   0.0 && p_Radius <= 350.0) coeffsBGidx = 0;
            else if (p_Radius > 350.0 && p_Radius <= 600.0) coeffsBGidx = 1;
            else                                            coeffsBGidx = 2;
        }
    }
    else {                                                                                  // pop II
             if (p_MassIndex ==  1 && p_Radius >   6.0 && p_Radius <  50.0) useLambdas = true;
        else if (p_MassIndex ==  2 && p_Radius >  36.0 && p_Radius <  53.0) useLambdas = true;
        else if (p_MassIndex ==  3 && p_Radius >  19.0 && p_Radius <  85.0) useLambdas = true;
        else if (p_MassIndex ==  4 && p_Radius >  85.0 && p_Radius < 120.0) useLambdas = true;
        else if (p_MassIndex ==  5 && p_Radius > 115.0 && p_Radius < 165.0) useLambdas = true;
        else if (p_MassIndex ==  6 && p_Radius > 150.0 && p_Radius < 210.0) useLambdas = true;
        else if (p_MassIndex ==  7 && p_Radius > 190.0 && p_Radius < 260.0) useLambdas = true;
        else if (p_MassIndex ==  8 && p_Radius > 180.0 && p_Radius < 300.0) useLambdas = true;
        else if (p_MassIndex == 10 && p_Radius > 200.0 && p_Radius < 410.0) useLambdas = true;
        else if (p_MassIndex == 11 && p_Radius > 250.0 && p_Radius < 490.0) useLambdas = true;
        else if (p_MassIndex == 12 && p_Radius > 200.0 && p_Radius < 570.0) useLambdas = true;
        else if (p_MassIndex == 13 && p_Radius > 230.0 && p_Radius < 755.0) useLambdas = true;
    }

    // get limits and (defined) lambdas
    std::tuple<NANJING_POP_LIMITS_LAMBDAS, NANJING_POP_LIMITS_LAMBDAS> evolStageLimitsLambdas = NANJING_LIMITS_LAMBDAS_ENHANCED[evolStage - 1];
    NANJING_POP_LIMITS_LAMBDAS                                         popLimitsLambdas       = p_StellarPop == STELLAR_POPULATION::POPULATION_I ? std::get<0>(evolStageValues) : std::get<1>(evolStageValues);
    std::tuple<NANJING_LIMITS_ENHANCED, NANJING_LAMBDAS>               limitsLambdas          = popValues[p_MassIndex];

    std::tuple<double, double, double> maxBGR = std::get<0>(limitsLambdas)[0];              // {maxB, maxG, maxR}

    double lambdaB;
    double lambdaG;
    if (useLambdas) {                                                                       // use lambdas defined in the paper?
                                                                                            // yes
        NANJING_LAMBDAS lambdaBG = std::get<1>(limitsLambdas)[0];                           // defined {lambdaB, lambdaG}

        lambdaB = std::get<0>(lambdaBG);
        lambdaG = std::get<1>(lambdaBG);
    }
    else {                                                                                  // no - calculate lambdas (per paper)

        // get B & G coefficients vector
        std::tuple<NANJING_POP_COEFFICIENTS, NANJING_POP_COEFFICIENTS> evolStageCoeffs = NANJING_COEFFICIENTS[evolStage - 1];
        NANJING_POP_COEFFICIENTS                                       popCoeffs       = p_StellarPop == STELLAR_POPULATION::POPULATION_I ? std::get<0>(evolStageCoeffs) : std::get<1>(evolStageCoeffs);
        std::tuple<DBL_VECTOR, DBL_VECTOR>                             BGcoeffs        = popCoeffs[p_MassIndex][coeffsBGidx];

        DBL_VECTOR Bcoeffs = std::get<0>(BGcoeffs);
        DBL_VECTOR Gcoeffs = std::get<1>(BGcoeffs);
        
        if (p_StellarPop == STELLAR_POPULATION::POPULATION_I && p_MassIndex == 11 && p_Radius > 69.0 && p_Radius < 126.0) { 
            lambdaB = 0.5 - (p_Radius * 8.77E-04);
            lambdaG = 0.18;
        }
        else {
            double Rin;
            if (p_StellarPop == STELLAR_POPULATION::POPULATION_I && p_MassIndex == 0)
                Rin = (p_Mass - p_CoreMass) / p_Mass;
            else
                Rin = std::min(p_Radius, std::get<2>(maxBGR));                              // clamp to maximum allowed radius (maxR) to prevent exceeding domain of the polynomial fits
            
            const double Rin2 = Rin  * Rin;
            const double Rin3 = Rin  * Rin2;
            const double Rin4 = Rin2 * Rin2;
            const double Rin5 = Rin2 * Rin3;

            lambdaB = Bcoeffs[0] + (Bcoeffs[1] * Rin) + (Bcoeffs[2] * Rin2) + (Bcoeffs[3] * Rin3) + (Bcoeffs[4] * Rin4) + (Bcoeffs[5] * Rin5);
            lambdaG = Gcoeffs[0] + (Gcoeffs[1] * Rin) + (Gcoeffs[2] * Rin2) + (Gcoeffs[3] * Rin3) + (Gcoeffs[4] * Rin4) + (Gcoeffs[5] * Rin5);

            if (p_StellarPop == STELLAR_POPULATION::POPULATION_I && p_MassIndex == 0) {
                lambdaB = 1.0 / lambdaB;
                lambdaG = 1.0 / lambdaG;                
            }
        }   
    }

    // Limit lambda to some 'reasonable' range
    lambdaG = std::min(std::max(0.05, lambdaG), std::min(1.0, std::get<1>(maxBGR)));        // clamp lambda G to [0.05, min(1, maxG)]
    lambdaB = std::max(std::min(lambdaB, std::get<0>(maxBGR)), std::max(0.05, lambdaG));    // clamp lambda B to [max(0.05, lambdaG), maxB]

    // Calculate lambda as some combination of lambdaB and lambdaG by
    // lambda = alpha_th • lambdaB + (1-alpha_th) • lambdaG
    // STARTRACK uses alpha_th = 1/2
    return (OPTIONS->CommonEnvelopeAlphaThermal() * lambdaB) + ((1.0 - OPTIONS->CommonEnvelopeAlphaThermal()) * lambdaG);
}


/*
 * CalculateLambdaNanjingStarTrack
 *
 * @brief
 * Calculate the common envelope lambda parameter, per Belczynski
 * (as implemented in StarTrack - courtesy Chris Belczynski)
 *
 * This function good for HG and FGB stars.
 *
 *
 * double CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Common envelope lambda parameter
 */
double CHeB::CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const {

    constexpr size_t evolStage = 2;                                                         // CHeB evolutionary stage from Xu & Li, 2010
                                           
    size_t coeffsBGidx = 0;                                                                 // index into coefficients vector
    size_t limitBGRidx = 0;                                                                 // index into limits vector
    size_t lambdaBGidx = -1;                                                                // index into lambdas vector: -ve indicates calculate lambdas

    // determine mass index based on p_Mass
    auto it = std::upper_bound(NANJING_MASSES_MIDPOINTS.begin(), NANJING_MASSES_MIDPOINTS.end(), p_Mass);
    const size_t massIndex = it != arr.end() ? std::distance(NANJING_MASSES_MIDPOINTS.begin(), it) : NANJING_MASSES_MIDPOINTS.size();

    if (GLOBALS->ReferenceMetallicity() > LAMBDA_NANJING_ZLIMIT_STARTRACK) {                // Z > LAMBDA_NANJING_ZLIMIT_STARTRACK?
                                                                                            // yes
             if (massIndex == 0 && p_Radius > 200.0) lambdaBGidx = 0;
        else if (massIndex == 1) {
                 if (p_Radius > 340.0)                  lambdaBGidx = 0;
            else if (p_Radius > 8.5 && p_Radius < 60.0) lambdaBGidx = 1;
        }
        else if (massIndex == 2) {
            if (p_Radius > 400.0) lambdaBGidx = 0;
            else                  limitBGRidx = 1;
        }
        else if (massIndex == 3) {
            if (p_Radius > 410.0) lambdaBGidx = 0;
            else                  limitBGRidx = 1;
        }
        else if (massIndex == 4 && p_Radius > 430.0) lambdaBGidx = 0;
        else if (massIndex == 5 && p_Radius > 440.0) lambdaBGidx = 0;
        else if (massIndex == 6 && p_Radius > 420.0) lambdaBGidx = 0;
        else if (massIndex == 7) {
            if (p_Radius > 490.0) lambdaBGidx = 0;
            else                  limitBGRidx = 1;
        }
        else if (massIndex == 8 && p_Radius > 530.0) lambdaBGidx = 0;
        else if (massIndex == 9 && p_Radius > 600.0) lambdaBGidx = 0;
        else if (p_MassIndex == 10) {
                 if (p_Radius > 850.0)                      lambdaBGidx = 0;
            else if (p_Radius >   0.0 && p_Radius <= 350.0) coeffsBGidx = 0;
            else if (p_Radius > 350.0 && p_Radius <= 600.0) coeffsBGidx = 1;
            else                                            coeffsBGidx = 2;
        }
        else if (massIndex == 11) {
                 if (p_Radius > 1000.0)                   lambdaBGidx = 0;
            else if (p_Radius > 69.0 && p_Radius < 126.0) lambdaBGidx = 1;
        }
        else if (massIndex == 12 && p_Radius > 1050.0) lambdaBGidx = 0;
        else if (massIndex == 13 && p_Radius > 1200.0) lambdaBGidx = 0;
    }
    else {                                                                                  // no - Z <= LAMBDA_NANJING_ZLIMIT_STARTRACK
             if (massIndex == 0 && p_Radius > 160.0) lambdaBGidx = 0;
        else if (massIndex == 1) {
                 if (p_Radius > 350.0)                  lambdaBGidx = 0;
            else if (p_Radius > 6.0 && p_Radius < 50.0) lambdaBGidx = 2;
        }
        else if (massIndex == 2) {
                 if (p_Radius > 400.0)                   lambdaBGidx = 0;
            else if (p_Radius > 36.0 && p_Radius < 53.0) lambdaBGidx = 1;
        }
        else if (massIndex ==  3) {
                 if (p_Radius > 410.0)                   lambdaBGidx = 0;
            else if (p_Radius > 19.0 && p_Radius < 85.0) lambdaBGidx = 1;
        }
        else if (massIndex == 4) {
                 if (p_Radius > 320.0)                    lambdaBGidx = 0;
            else if (p_Radius > 85.0 && p_Radius < 120.0) lambdaBGidx = 1;
        }
        else if (massIndex == 5) {
                 if (p_Radius > 330.0)                     lambdaBGidx = 0;
            else if (p_Radius > 115.0 && p_Radius < 165.0) lambdaBGidx = 1;
        }
        else if (massIndex == 6) {
                 if (p_Radius > 360.0)                     lambdaBGidx = 0;
            else if (p_Radius > 150.0 && p_Radius < 210.0) lambdaBGidx = 1;
        }
        else if (massIndex == 7) {
                 if (p_Radius > 400.0)                     lambdaBGidx = 0;
            else if (p_Radius > 190.0 && p_Radius < 260.0) lambdaBGidx = 1;
        }
        else if (massIndex == 8) {
                 if (p_Radius > 440.0)                     lambdaBGidx = 0;
            else if (p_Radius > 180.0 && p_Radius < 300.0) lambdaBGidx = 1;
        }
        else if (massIndex ==  9 && p_Radius > 500.0) lambdaBGidx = 0;
        else if (massIndex == 10) {
                 if (p_Radius > 600.0)                     lambdaBGidx = 0;
            else if (p_Radius > 200.0 && p_Radius < 410.0) lambdaBGidx = 1;
        }
        else if (massIndex == 11) {
                 if (p_Radius > 650.0)                     lambdaBGidx = 0;
            else if (p_Radius > 250.0 && p_Radius < 490.0) lambdaBGidx = 1;
        }
        else if (massIndex == 12) {
                 if (p_Radius > 750.0)                     lambdaBGidx = 0;
            else if (p_Radius > 200.0 && p_Radius < 570.0) lambdaBGidx = 1;
        }
        else if (massIndex == 13) {
                 if (p_Radius > 900.0)                     lambdaBGidx = 0;
            else if (p_Radius > 230.0 && p_Radius < 755.0) lambdaBGidx = 1;
        }
    }

    // get limits and (defined) lambdas
    NANJING_Z_LIMITS_LAMBDAS                              ZlimitsLambdas = GLOBALS->ReferenceMetallicity() < LAMBDA_NANJING_ZLIMIT ? std::get<0>(NANJING_LIMITS_LAMBDAS_STARTRACK) : std::get<1>(NANJING_LIMITS_LAMBDAS_STARTRACK);
    std::tuple<NANJING_LIMITS_STARTRACK, NANJING_LAMBDAS> limitsLambdas  = ZlimitsLambdas[p_MassIndex];

    std::tuple<double, double> maxBG = std::get<0>(limitsLambdas)[limitBGidx];              // {maxB, maxG}

    double lambdaB;
    double lambdaG;
    if (useLambdas) {                                                                       // use lambdas defined by StarTrack?
                                                                                            // yes
        std::tuple<double, double> lambdaBG = std::get<1>(limitsLambdas)[lambdaBGidx];      // defined {lambdaB, lambdaG}

        lambdaB = std::get<0>(lambdaBG);
        lambdaG = std::get<1>(lambdaBG);
    }
    else {                                                                                  // no - calculate lambdas (per StarTrack)

        // get B & G coefficients vector
        std::tuple<NANJING_POP_COEFFICIENTS, NANJING_POP_COEFFICIENTS> evolStageCoeffs = NANJING_COEFFICIENTS[evolStage - 1];
        NANJING_POP_COEFFICIENTS                                       ZCoeffs         = GLOBALS->ReferenceMetallicity() < LAMBDA_NANJING_ZLIMIT ? std::get<0>(evolStageCoeffs) : std::get<1>(evolStageCoeffs);
        std::tuple<DBL_VECTOR, DBL_VECTOR>                             BGcoeffs        = ZCoeffs[p_MassIndex][coeffsBGidx];

        DBL_VECTOR Bcoeffs = std::get<0>(BGcoeffs);
        DBL_VECTOR Gcoeffs = std::get<1>(BGcoeffs);
        
        double Rin = p_Radius;

        if (GLOBALS->ReferenceMetallicity() < LAMBDA_NANJING_ZLIMIT && p_MassIndex == 0 && p_Radius > 2.7) {
            lambdaB = 2.33 - (Rin * 9.18E-03);
            lambdaG = 1.12 - (Rin * 4.59E-03);
        }
        else if (GLOBALS->ReferenceMetallicity() < LAMBDA_NANJING_ZLIMIT && p_MassIndex == 13) {
            lambdaB = 1.2 * exp(-Rin / 90.0);
            lambdaG = 0.55 * exp(-Rin / 160.0);
        }
        else if (GLOBALS->ReferenceMetallicity() >= LAMBDA_NANJING_ZLIMIT && p_MassIndex == 0 && p_Radius > 12.0) {
            lambdaB = 1.8 * exp(-Rin / 80.0);
            lambdaG = exp(-Rin / 45.0);
        }
        else if (GLOBALS->ReferenceMetallicity() >= LAMBDA_NANJING_ZLIMITI && p_MassIndex == 9) {
            const double tmp = exp(-Rin / 35.0);
            lambdaB = 1.75 * tmp;
            lambdaG = 0.9 * tmp;
        }
        else {
            if (GLOBALS->ReferenceMetallicity() < LAMBDA_NANJING_ZLIMIT && p_MassIndex == 0) Rin = (p_Mass - p_CoreMass) / p_Mass;
            
            const double Rin2 = Rin  * Rin;
            const double Rin3 = Rin  * Rin2;
            const double Rin4 = Rin2 * Rin2;
            const double Rin5 = Rin2 * Rin3;

            lambdaB = Bcoeffs[0] + (Bcoeffs[1] * Rin) + (Bcoeffs[2] * Rin2) + (Bcoeffs[3] * Rin3) + (Bcoeffs[4] * Rin4) + (Bcoeffs[5] * Rin5);
            lambdaG = Gcoeffs[0] + (Gcoeffs[1] * Rin) + (Gcoeffs[2] * Rin2) + (Gcoeffs[3] * Rin3) + (Gcoeffs[4] * Rin4) + (Gcoeffs[5] * Rin5);

            if (GLOBALS->ReferenceMetallicity() < LAMBDA_NANJING_ZLIMIT && p_MassIndex == 0) {
                lambdaB = 1.0 / lambdaB;
                lambdaG = 1.0 / lambdaG;                
            }
        }       
    }

    // Limit lambda to some 'reasonable' range
    lambdaG = std::min(std::max(0.05, lambdaG), std::min(1.0, std::get<1>(maxBGR)));        // clamp lambda G to [0.05, min(1, maxG)]
    lambdaB = std::max(std::min(lambdaB, std::get<0>(maxBGR)), std::max(0.05, lambdaG));    // clamp lambda B to [max(0.05, lambdaG), maxB]

    // Calculate lambda as some combination of lambdaB and lambdaG by
    // lambda = alpha_th • lambdaB + (1-alpha_th) • lambdaG
    // STARTRACK uses alpha_th = 1/2
    return (OPTIONS->CommonEnvelopeAlphaThermal() * lambdaB) + ((1.0 - OPTIONS->CommonEnvelopeAlphaThermal()) * lambdaG);
}