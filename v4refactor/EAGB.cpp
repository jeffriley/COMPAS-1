#include "EAGB.h"
#include "Star.h"
#include "HeMS.h"
#include "HeGB.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
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
 * DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Metallicity, const double p_Mass, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_Timescales) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBParams                      Hurley GB parameters
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
COMPAS_PURE DBL_VECTOR EAGB::CalculateTimescales_Hurley2000(const double p_Metallicity, const double p_Mass, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const {

// #defines for convenience and readability - undefined at end of function
#define GBParams(x) p_GBParams[static_cast<int>(HURLEY_GB_PARAMETERS::x)]
#define tScales(x) p_tScales[static_cast<int>(HURLEY_TIMESCALES::x)]

// JR: CHECK SSE - why don't we recalculate GB timescales here?????????? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    const double p1   = GBParams(p) - 1.0;
    const double q1   = GBParams(q) - 1.0;
    const double p1_p = p1 / GBParams(p);
    const double q1_q = q1 / GBParams(q);

    const double tBAGB = CalculateLifetimeToBAGB_Hurley2000(tScales(HeI), tScales(He));
    const double lBAGB = CalculateLuminosityAtBAGB_Hurley2000(p_Metallicity, p_Mass);

    DBL_VECTOR tScales = p_tScales; // copy given timescales

    tScales[static_cast<int>(HURLEY_TIMESCALES::Inf1_FAGB)] = tBAGB + ((1.0 / (p1 * GBParams(AHe) * GBParams(D))) * PPOW((GBParams(D) / lBAGB), p1_p));
    tScales[static_cast<int>(HURLEY_TIMESCALES::Mx_FAGB)]   = tScales(Inf1_FAGB) - ((tScales(Inf1_FAGB) - tBAGB) * PPOW((lBAGB / GBParams(Lx)), p1_p));
    tScales[static_cast<int>(HURLEY_TIMESCALES::Inf2_FAGB)] = tScales(Mx_FAGB) + ((1.0 / (q1 * GBParams(AHe) * GBParams(B))) * PPOW((GBParams(B) / GBParams(Lx)), q1_q));

    // return timescales vector by value - NRVO takes care of performance/efficiency
    return tScales;

#undef timescales
#undef GBParams    
}


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
 * CalculateRadius_Hurley2000
 *
 * @brief
 * Calculate radius on the Early Asymptotic Giant Branch, per Hurley et al. 2000, eq 74
 *
 *
 * static double CalculateRadius_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Luminosity)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_BCoefficients                 ZDEP->HurleyBCoefficients(Metallicity()) at call site
 * @param       p_MHeF                          Maximum initial mass for which helium ignites degenerately in a Helium Flash (HeF)
 * @return                                      EAGB radius (Rsol)
 */
double EAGB::CalculateRadius_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Luminosity) {

    // sanity check for mass and luminosity - just return 0.0 if mass or luminosity <= 0
    if (!(p_Mass > 0.0) || !(p_Luminosity > 0.0)) return 0.0;

    const double mHeF = ZDEP->HurleyMCOHeF(p_Metallicity);
    const auto&  b    = ZDEP->HurleyBCoeffs(p_Metallicity); // Hurley b coefficients

    // calculate radius constant A (Hurley et al. 2000, eq 74), and coefficient b(50)
    double A;
    double b50;

    if (mHeF < p_Mass) {                                    // Mass above He flash threshold?
                                                            // Yes
        A   = std::min((b[51] * PPOW(p_Mass, -b[52])), (b[53] * PPOW(p_Mass, -b[54])));
        b50 = b[55] * b[3];
    }
    else if ((mHeF - 0.2) > p_Mass) {                       // No - mass below He flash threshold?
                                                            // Yes
        A   = b[56] + (b[57] * p_Mass);   
        b50 = b[3];
    }
    else {                                                  // No - linear interpolation between end points
        const double x1    = mHeF - 0.2;
        const double x2    = mHeF;
        const double x2_x1 = x2 - x1;

        double y1        = b[56] + (b[57] * x1);
        double y2        = std::min((b[51] * PPOW(x2, -b[52])), (b[53] * PPOW(x2, -b[54])));
        double gradient  = (y2 - y1) / x2_x1;
        double intercept = y2 - (gradient * x2);

               A         = (gradient * p_Mass) + intercept;

               y1        = b[3];
               y2        = b[55] * b[3];
               gradient  = (y2 - y1) / x2_x1;
               intercept = y2 - (gradient * x2);

               b50       = (gradient * p_Mass) + intercept;
    }

    return A * (PPOW(p_Luminosity, b[1]) + (b[2] * PPOW(p_Luminosity, b50)));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               MASS LOSS / ACCRETION                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateMLRate_Hurley2000
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for EAGB stars,
 * per Hurley et al. 2000.
 * 
 * 
 * MassLossT CalculateMLRate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_Perturb) const
 *     
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Perturb                       Small envelope perturbation parameter
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (may be ML_TYPE::NONE)
 */
COMPAS_PURE MassLossT EAGB::CalculateMLRate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_Perturb) const {
   
    // Calculate max GB mass loss rate - default rate

    const double dMdtWR = std::get<0>(CalculateMLRate_Hurley2000(p_Mass, p_Radius, p_Luminosity, p_Perturb));
    const double dMdtKR = std::get<0>(CalculateMLRate_KudritzkiReimers1978(p_Mass, p_Radius, p_Luminosity));
    const double dMdtNJ = std::get<0>(CalculateMLRate_NieuwenhuijzenDeJager1990(p_Mass, p_Radius, p_Luminosity));
    const double dMdtVW = std::get<0>(CalculateMLRate_VassiliadisWood1993(p_Mass, p_Radius, p_Luminosity));

    double dMdt = std::max(dMdtVW, std::max(dMdtNJ, dMdtKR));   // Max GB mass loss rate

    ML_TYPE DominantMassLossType = ML_TYPE::GB;                 // Default dominant mass loss type is GB

    if (dMdtWR > dMdt) {                                        // WR rate > max GB rate?
        dominantMassLossType = ML_TYPE::WR;                     // Yes - set dominant type to WR
        dMdt = dMdtWR;                                          // Use WR ML rate
    }

    return std::make_tuple(dMdt, dominantMassLossType);
}







///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Modify the star after it loses its envelope
 *
 * Hurley et al. 2000, section 6 just before eq 76 and after Eq. 105
 *
 * Where necessary updates attributes of star (depending upon stellar type):
 *
 *     - m_StellarType
 *     - m_Timescales
 *     - m_GBParams
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
STELLAR_TYPE EAGB::ResolveEnvelopeLoss(bool p_Force) {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]  // for convenience and readability - undefined at end of function
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]            // for convenience and readability - undefined at end of function

    STELLAR_TYPE stellarType = m_StellarType;

    if (ShouldEnvelopeBeExpelledByPulsations()) m_EnvelopeJustExpelledByPulsations = true;

    if (p_Force || utils::Compare(m_CoreMass, m_Mass) >= 0 || m_EnvelopeJustExpelledByPulsations) {                                         // Envelope lost, form an evolved naked helium giant

        m_Mass       = std::min(m_CoreMass, m_Mass);
        m_HeCoreMass = m_Mass;
        m_Mass0      = m_Mass;
        m_CoreMass   = m_COCoreMass;

        double p1    = gbParams(p) - 1.0;
        double q1    = gbParams(q) - 1.0;
        double p1_p  = p1 / gbParams(p);
        double q1_q  = q1 / gbParams(q);

        timescales(tHeMS) = HeMS::CalculateLifetimeOnPhase_Static(m_Mass);                                                                  // calculate common values

        double LTHe = HeMS::CalculateLuminosityAtPhaseEnd_Static(m_Mass);

        timescales(tinf1_HeGB) = timescales(tHeMS) + (1.0 / ((p1 * gbParams(AHe) * gbParams(D))) * PPOW((gbParams(D) / LTHe), p1_p));
        timescales(tx_HeGB)    = timescales(tinf1_HeGB) - (timescales(tinf1_HeGB) - timescales(tHeMS)) * PPOW((LTHe / gbParams(Lx)), p1_p);
        timescales(tinf2_HeGB) = timescales(tx_HeGB) + ((1.0 / (q1 * gbParams(AHe) * gbParams(B))) * PPOW((gbParams(B) / gbParams(Lx)), q1_q));

        m_Age = HeGB::CalculateAgeOnPhase_Static(m_Mass, m_COCoreMass, timescales(tHeMS), m_GBParams);

        HeHG::CalculateGBParams_Static(m_Mass0, m_Mass, LogMetallicityXiHurley(), m_MassCutoffs, m_AnCoefficients, m_BnCoefficients, m_GBParams); 
        m_Luminosity = HeGB::CalculateLuminosityOnPhase_Static(m_COCoreMass, gbParams(B), gbParams(D));

        double R1, R2;
        std::tie(R1, R2) = HeGB::CalculateRadiusOnPhase_Static(m_Mass, m_Luminosity);
        if (utils::Compare(R1, R2) < 0) {
            m_Radius    = R1;
            stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;
        }
        else {
            m_Radius    = R2;
            stellarType = STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH;                                 // Has a deep convective envelope
        }
    }

    return stellarType;

#undef gbParams
#undef timescales
}


/*
 * Determine if this phase should be skipped entirely
 *
 *
 * bool ShouldSkipPhase()
 *
 * @return                                      Boolean flag: true if this phase should be skipped, false if not
 */
bool EAGB::ShouldSkipPhase() const {

    double McCOBAGB = CalculateCOCoreMass_Hurley2000(Timescale(HURLEY_TIMESCALES::HeI) + Timescale(HURLEY_TIMESCALES::He), CurrentState().HurleyGBParams(), CurrentState().HurleyTimescales());
    double McSN     = std::max(CalculateCoreMassAtSN(MCH, CurrentState().HurleyGBParams(HURLEY_GB_PARAMETERS::McBAGB)), 1.05 * McCOBAGB);                       // hack from Hurley fortran code, doesn't seem to be in the paper

    return (utils::Compare(McSN, COCoreMass()) < 0);                                            // skip phase if core is heavy enough to go supernova
}


/*
 * Determine if evolution should continue on this phase, or whether evolution
 * on this phase should end (and so evolve to next phase)
 *
 *
 * bool ShouldEvolveOnPhase()
 *
 * @return                                      Boolean flag: true if evolution on this phase should continue, false if not
 */
bool EAGB::ShouldEvolveOnPhase() const {

	double tDU = CalculateLifetimeTo2ndDredgeUp(Timescale(HURLEY_TIMESCALES::Inf1_FAGB), Timescale(HURLEY_TIMESCALES::Inf2_FAGB));
    return ((utils::Compare(Age(), tDU) < 0 || utils::Compare(CurrentState().HurleyGBParams(HURLEY_GB_PARAMETERS::McBAGB), MCBUR2) >= 0) && !ShouldEnvelopeBeExpelledByPulsations());
}


/*
 * Determine if star should continue evolution as a Supernova
 *
 *
 * bool IsSupernova()
 *
 * @return                                      Boolean flag: true if star has gone Supernova, false if not
 */
bool EAGB::IsSupernova() const {

    double McCOBAGB = CalculateCOCoreMass_Hurley2000(Timescale(HURLEY_TIMESCALES::HeI) + Timescale(HURLEY_TIMESCALES::He), CurrentState().HurleyGBParams(), CurrentState().HurleyTimescales());
    double McSN     = std::max(CalculateCoreMassAtSN(MCH, CurrentState().HurleyGBParams(HURLEY_GB_PARAMETERS::McBAGB)), 1.05 * McCOBAGB);                                // hack from Hurley fortran code, doesn't seem to be in the paper   JR: do we know why? **Ilya**

    return (utils::Compare(McSN, COCoreMass()) <= 0);                                           // core is heavy enough to go Supernova
}


/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE EAGB::EvolveToNextPhase() {
    return STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH;
}




/// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                LAMBDA CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCELambda_Nanjing_Enhanced
 *
 * @brief
 * Calculate the common envelope lambda parameter, per Xu & Li, 2010 (Nanjing - for the university)
 * (https://arxiv.org/abs/1004.4957, v1, 28Apr2010)
 * (https://iopscience.iop.org/article/10.1088/0004-637X/716/1/114)
 *
 * This function good for EAGB stars.
 *
 *
 * double CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
 *                                           const double             p_Radius,
 *                                           const double             p_CoreMass,
 *                                           const SizeT              p_MassIndex,
 *                                           const STELLAR_POPULATION p_StellarPop) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MassIndex                     Index of mass bin in NANJING_MASSES (see constants.h)
 * @param       p_StellarPop                    Stellar population (POP I or POP II)
 * @return                                      Common envelope lambda parameter
 */
double EAGB::CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
                                                const double             p_Radius,
                                                const double             p_CoreMass,
                                                const SizeT              p_MassIndex,
                                                const STELLAR_POPULATION p_StellarPop) const {

    [[maybe_unused]] constexpr SizeT evolStage = 3;                                                          // EAGB evolutionary stage from Xu & Li, 2010

    SizeT coeffsBGidx = 0;                                                                  // Index into coefficients array
    bool  useLambdas  = false;                                                              // Flag - use lambdas defined in the paper

    ApplyNanjingEnhancedDecisions(
        p_StellarPop == STELLAR_POPULATION::POPULATION_I ? NANJING_DECISIONS_EAGB_POPI : NANJING_DECISIONS_EAGB_POPII,
        p_MassIndex, p_Radius,
        useLambdas, coeffsBGidx);

    // Get limits and (defined) lambdas
    const auto& enData   = (p_StellarPop == STELLAR_POPULATION::POPULATION_I) ? NANJING_DATA_EAGB_POPI : NANJING_DATA_EAGB_POPII;
    const auto& massData = enData[p_MassIndex];
    std::tuple<double, double, double> maxBGR = std::make_tuple(massData.maxB, massData.maxG, massData.maxR);   // {maxB, maxG, maxR}

    double lambdaB;
    double lambdaG;
    if (useLambdas) {                                                                       // Use lambdas defined in the paper?
                                                                                            // Yes
        lambdaB = massData.definedLambdaB;
        lambdaG = massData.definedLambdaG;
    }
    else {                                                                                  // No - calculate lambdas (per paper)

        // Get B & G coefficients vector
        const auto&       coeffsTable = (p_StellarPop == STELLAR_POPULATION::POPULATION_I) ? NANJING_COEFFS_EAGB_POPI : NANJING_COEFFS_EAGB_POPII;
        const auto&       BGcoeffs    = coeffsTable[p_MassIndex][coeffsBGidx];
        const DBL_VECTOR& Bcoeffs     = std::get<0>(BGcoeffs);
        const DBL_VECTOR& Gcoeffs     = std::get<1>(BGcoeffs);

        const double rinClamped  = std::min(p_Radius, std::get<2>(maxBGR));
        const double envMassFrac = (p_Mass - p_CoreMass) / p_Mass;                          // not used by any EAGB rule, but supplied for helper-signature uniformity

        const auto& branchRules = (p_StellarPop == STELLAR_POPULATION::POPULATION_I) ? NANJING_POLY_BRANCH_EAGB_POPI : NANJING_POLY_BRANCH_EAGB_POPII;

        if (!ApplyNanjingPolyBranchOverride(branchRules, p_MassIndex, p_Radius, rinClamped, envMassFrac,
                                            Bcoeffs, Gcoeffs, lambdaB, lambdaG)) {
            // No rule matched -- default 5th-order polynomial on RIN_CLAMPED
            const double Rin  = rinClamped;
            const double Rin2 = Rin  * Rin;
            const double Rin3 = Rin  * Rin2;
            const double Rin4 = Rin2 * Rin2;
            const double Rin5 = Rin2 * Rin3;

            lambdaB = Bcoeffs[0] + (Bcoeffs[1] * Rin) + (Bcoeffs[2] * Rin2) + (Bcoeffs[3] * Rin3) + (Bcoeffs[4] * Rin4) + (Bcoeffs[5] * Rin5);
            lambdaG = Gcoeffs[0] + (Gcoeffs[1] * Rin) + (Gcoeffs[2] * Rin2) + (Gcoeffs[3] * Rin3) + (Gcoeffs[4] * Rin4) + (Gcoeffs[5] * Rin5);
        }
    }

    // Limit lambda to some 'reasonable' range
    lambdaG = std::min(std::max(0.05, lambdaG), std::min(1.0, std::get<1>(maxBGR)));        // Clamp lambda G to [0.05, min(1, maxG)]
    lambdaB = std::max(std::min(lambdaB, std::get<0>(maxBGR)), std::max(0.05, lambdaG));    // Clamp lambda B to [max(0.05, lambdaG), maxB]

    // Calculate lambda as some combination of lambdaB and lambdaG by
    // lambda = alpha_th • lambdaB + (1-alpha_th) • lambdaG
    // STARTRACK uses alpha_th = 1/2
    return (OPTIONS->CommonEnvelopeAlphaThermal() * lambdaB) + ((1.0 - OPTIONS->CommonEnvelopeAlphaThermal()) * lambdaG);
}


/*
 * Calculate the common envelope lambda parameter using the "Nanjing" prescription
 * from X.-J. Xu and X.-D. Li arXiv:1004.4957 (v1, 28Apr2010) as implemented in STARTRACK
 *
 * This implementation adapted from the STARTRACK implementation (STARTRACK courtesy Chris Belczynski)
 *
 *
 * JR: todo: the coefficients and factors here are hard-coded until I figure out an
 * efficient way of putting them in constants.h.  Because they are indexed by a few
 * things: stellar type, metallicity and ZAMS mass the easiet thing would be to  put
 * them in a map - but because they can be re-calculated at every timestep the hashing
 * overhead becomes a performance concern.  Vectors are a good alternative, but I need
 * to figure out how best to structure them for reasonable (and inuitive) access.
 *
 * This function good for EAGB stars.
 *
 *
 * double CalculateLambdaNanjingStarTrack(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass
 * 
 * @return                                      Nanjing lambda for use in common envelope
 */
COMPAS_PURE double EAGB_Constituent::CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const {

    const bool  isHighZ = utils::Compare(m_Star->Metallicity(), LAMBDA_NANJING_ZLIMIT_STARTRACK) > 0;
    const auto& bins    = isHighZ ? NANJING_STARTRACK_BINS_EAGB_POPI : NANJING_STARTRACK_BINS_EAGB_POPII;
    const auto& bin     = bins[FindNanjingStarTrackMassBin(bins, p_Mass)];

    DBL_VECTOR lambdaBG{};
    SizeT      coeffsIdx = 0;
    const bool gotLambda = ApplyNanjingStarTrackBinWindows(bin, p_Radius, lambdaBG, coeffsIdx);

    if (!gotLambda) {
        DBL_VECTOR a, b;

        const double x  = (bin.polyMode == NanjingPolyMode::INVERSE) ? (m_Star->Mass() - p_CoreMass) / m_Star->Mass() : p_Radius;
        const double x2 = x  * x;
        const double x3 = x2 * x;
        const double x4 = x2 * x2;
        const double x5 = x3 * x2;
        const double y1 = a[0] + (a[1] * x) + (a[2] * x2) + (a[3] * x3) + (a[4] * x4) + (a[5] * x5);
        const double y2 = b[0] + (b[1] * x) + (b[2] * x2) + (b[3] * x3) + (b[4] * x4) + (b[5] * x5);

        switch (bin.polyMode) {
            case NanjingPolyMode::INVERSE:  lambdaBG = { 1.0 / y1, 1.0 / y2 };  break;
            case NanjingPolyMode::EXP10:    lambdaBG = { pow(10.0, y1), y2 };   break;
            case NanjingPolyMode::IDENTITY: lambdaBG = { y1, y2 };              break;
        }
    }

    // Clamp lambda to [0.05, maxBG]
    lambdaBG[0] = std::min(std::max(0.05, lambdaBG[0]), std::get<0>(bin.maxBG));
    lambdaBG[1] = std::min(std::max(0.05, lambdaBG[1]), std::get<1>(bin.maxBG));

    // Calculate lambda as alpha-thermal-weighted combination of lambda_b and lambda_g.
    // Note that this is different from STARTRACK.
    return (OPTIONS->CommonEnvelopeAlphaThermal() * lambdaBG[0])
         + ((1.0 - OPTIONS->CommonEnvelopeAlphaThermal()) * lambdaBG[1]);
}








