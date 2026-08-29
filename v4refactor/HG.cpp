#include "HG.h"
#include "HeMS.h"
#include "HeWD.h"
#include "Star.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                      COEFFICIENTS, CONSTANTS etc. FUNCTIONS                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////





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
 * Calculate the radius on the Hertzsprung Gap, using a modified version of
 * Hurley et al. 2000, eq 27
 * 
 * See Hurley SSE code `hrdiag.f` lines 92, 188-203.  Here we replace the numerator, REHG, with the
 * GB radius if mass is below the threshold for He ignition, and a calculated value if mass is above
 * the threshold for He ignition (see code below)
 *
 *
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function nominally relies on the ZAMS values of one or more attributes of the star,
 * and should not be used if the ZAMS values are required but are not known.
 * 
 * The function *could* be called with something other than ZAMS values as parameters - it
 * is up to the caller what values are actually passed. In such cases the ZAMS warning does
 * not apply.  This function exists in the main code rather than in the GLOBALS module in
 * recognition of the usefulness of such cases.
 * 
 *
 * double CalculateRadiusOnPhase_Hurley2000(const double      p_Mass,
 *                                          const double      p_Luminosity,
 *                                          const double      p_Tau,
 *                                          const double      p_RZAMS,
 *                                          const double      p_MHeF,
 *                                          const double      p_MFGB,
 *                                          const double      p_Alpha1,
 *                                          const DBL_VECTOR& p_bN) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_MHeF                          Maximum initial mass at Helium Flash (Hurley masscutoffs[MHeF]) (Msol)
 * @param       p_MFGB                          Maximum initial mass at helium ignition on the FGB (Hurley masscutoffs[MFGB]) (Msol)
 * @param       p_Alpha1                        Hurley alpha1 constant
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      HG radius (Rsol)
 */
double HG::CalculateRadiusOnPhase_Hurley2000(const double      p_Mass,
                                             const double      p_Luminosity, 
                                             const double      p_Tau,
                                             const double      p_RZAMS,
                                             const double      p_MHeF,
                                             const double      p_MFGB,
                                             const double      p_Alpha1,
                                             const DBL_VECTOR& p_bN) const {

                                                    // p_Alpha1

    // FIX THIS - REMOVE BRCEK CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    [[maybe_unused]] double RTMS;  
    if ((OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK) && (utils::Compare(MZAMS(), BRCEK_LOWER_MASS_LIMIT) >= 0))
        // p_Mass generally has the value of m_Mass0, but since m_Mass is used for radius calculations on the MS and m_Mass0
        // is updated to a new value when BRCEK prescription is used, we need to use m_Mass here to keep radius continuous
        RTMS = MainSequence::CalculateRadiusAtPhaseEnd_Hurley2000(Mass(), p_RZAMS);
    else
        RTMS = MainSequence::CalculateRadiusAtPhaseEnd_Hurley2000(p_Mass, p_RZAMS);
    // FIX THIS - REMOVE BRCEK CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



    const double rTAMS = MainSequence::CalculateRadiusAtPhaseEnd_Hurley2000(p_Mass, p_RZAMS);
    const double rGB   = GiantBranch::CalculateRadius_Hurley2000(p_Mass, p_Luminosity);

    double rx = rGB;                                                                                                   // rx in Hurley SSE Fortran code

    if (p_Mass > p_MFGB) {                                                                      // mass above threshold for He ignition?
                                                                                                                        // yes
        // rMinHe is Hurley et al. 2000, eq 55 - first part (M >= MHeF)
        const double Mb28   = PPOW(p_Mass, p_bN[28]);                                                                   // pow() is slow - do it once only
        const double rMinHe = ((p_bN[24] * p_Mass) + (PPOW((p_bN[25] * p_Mass), p_bN[26]) * Mb28)) / (p_bN[27] + Mb28); // rmin in Hurley SSE Fortran code
        const double lum    = GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass, ZDEP->HurleyBCoefficients(Metallicity()), ZDEP->HurleyMassCutoffs(Metallicity(), HURLEY_MCO::HeF), ZDEP->HurleyAlpha1(Metallicity()));

        // mass is used here in the Hurley SSE Fortran code, mass0 everywhere else
        // ry in Hurley SSE Fortran code
        const double ry = EAGB::CalculateRadius_Hurley2000_Static(p_Mass, lum, p_bN, p_MHeF);

        // calculate radius at He ignition for MFGB < p_Mass < HM
        // Hurley et al. 2000, eq 50
        
        rx = std::min(rMinHe, ry);
        
        if (p_Mass < HIGH_MASS_THRESHOLD) {
            rx = rMinHe * PPOW(rGB / rMinHe, std::log10(p_Mass / HIGH_MASS_THRESHOLD) / std::log10(p_MFGB / HIGH_MASS_THRESHOLD));
        }

        // this piece of code resets rx if the blue loop is relatively short
        // see Hurley SSE Fortran code, function tblf() in `zfuncs1.f`
        // JR: I suspect this is where the check came from in the discussion re blue loop in CHeB::CalculateTimescales() - I don't like this much... **Ilya**
        const double r1 = std::max(1.0 - rMinHe / ry, 1.0E-12); 

        double tblf = (1.0 - p_bN[47]) * PPOW(p_Mass, p_bN[48]) * PPOW(r1, p_bN[49]);                                        // calculate blue-loop fraction of He-burning
        tblf = std::min(1.0, std::max(0.0, tblf));                                                                  // clamp to [0.0, 1.0]

        if (tblf < MINIMUM_BLUE_LOOP_FRACTION) rx = ry;                                                             // reset rx if short blue loop
    }

    return rTAMS * PPOW(rx / rTAMS, p_Tau);
}



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////










///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                            LIFETIME / AGE CALCULATIONS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////




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
 * STELLAR_TYPE ResolveEnvelopeLoss(bool p_Force)
 *
 * @param   [IN]    p_Force                     Boolean to indicate whether the resolution of the loss of the envelope should be performed
 *                                              without checking the precondition(s).
 *                                              Default is false.
 *
 * @return                                      Stellar Type to which star should evolve after losing envelope
 */
STELLAR_TYPE HG::ResolveEnvelopeLoss(bool p_Force) {
#define timescales(x) m_Timescales[static_cast<int>(TIMESCALE::x)]              // for convenience and readability - undefined at end of function
#define massCutoffs(x) m_MassCutoffs[static_cast<int>(MASS_CUTOFF::x)]          // for convenience and readability - undefined at end of function

    STELLAR_TYPE stellarType = m_StellarType;

    if (p_Force || utils::Compare(m_CoreMass, m_Mass) >= 0) {                   // envelope loss

        m_Mass = std::min(m_CoreMass, m_Mass);

        if (utils::Compare(m_Mass0, massCutoffs(MHeF)) < 0) {                   // star evolves to Helium White Dwarf

            stellarType  = STELLAR_TYPE::HELIUM_WHITE_DWARF;

            m_Radius     = HeWD::CalculateRadiusOnPhase_Static(m_Mass);
            m_Age        = 0.0;                                                 // see Hurley et al. 2000, discussion after eq 76
        }
        else {                                                                  // star evolves to Zero age Naked Helium Main Star

            stellarType  = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;

            m_Mass0      = m_Mass;
            m_Radius     = HeMS::CalculateRadiusAtZAMS_Static(m_Mass);          
            m_Luminosity = HeMS::CalculateLuminosityAtZAMS_Static(m_Mass);
            m_Age        = 0.0;                                                 // can't use Hurley et al. 2000, eq 76 here - timescales(tHe) not calculated yet
        }
    }

    return stellarType;

#undef massCutoffs
#undef timescales
}


/*
 * Set parameters for evolution to next phase and return Stellar Type for next phase
 *
 *
 * STELLAR_TYPE EvolveToNextPhase()
 *
 * @return                                      Stellar Type for next phase
 */
STELLAR_TYPE HG::EvolveToNextPhase() {
#define massCutoffs(x) ZDEP->HurleyMassCutoffs(Metallicity(), HURLEY_MCO::x)

    STELLAR_TYPE stellarType;

    if (utils::Compare(Mass0(), massCutoffs(FGB)) < 0) {
        stellarType = STELLAR_TYPE::FIRST_GIANT_BRANCH;
    }
    else {
        stellarType = STELLAR_TYPE::CORE_HELIUM_BURNING;
    }    

    return stellarType;
#undef massCutoffs
}





//// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                LAMBDA CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate the common envelope lambda parameter using the "Dewi" prescription.
 *
 * Fit from Appendix A in Claeys+2014, based on the method described in Dewi & Tauris 2000
 * arXiv:1401.2895 for Claeys+2014
 * arXiv:0007034 for Dewi and Tauris 2000
 *
 * ALEJANDRO - 17/05/2017 - Not fully tested, nor fully coded, nor fully trusted. Any other \lambda prescription is personally preferred.
 * Missing (A.6),(A.7),(A.8),(A.9),(A.10),(A.11) and (A.12), which have to do with ionization energy
 *
 *
 * double CalculateCELambda_Dewi()
 *
 * @return                                      Dewi lambda for use in common envelope
 */
double HG_Constituent::CalculateCELambda_Dewi() const {

    const double m       = m_Star->Mass();
    const double r       = m_Star->Radius();
    const double rZAMS   = m_Star->RZAMS();
    const double lum     = m_Star->Luminosity();
    const double cmass   = m_Star->CoreMass();

	double lambda1 = std::min(0.80, (3.0 / (2.4 + PPOW(m,-3.0 / 2.0))) - (0.15 * std::log10(lum)));                           // (A.3) Claeys+2014
	double lambda2 = 0.42 * PPOW(rZAMS / r, 0.4);                                                                            // (A.2) Claeys+2014
	double envMass = utils::Compare(cmass, 0.0) > 0 && utils::Compare(m, cmass) > 0 ? m - cmass : 0.0;

    double lambdaCE;

         if (utils::Compare(envMass, 1.0) >= 0) lambdaCE = 2.0 * lambda1;                                                   // (A.1) Bottom, Claeys+2014
	else if (utils::Compare(envMass, 0.0) >  0) lambdaCE = 2.0 * (lambda2 + (std::sqrt(envMass) * (lambda1 - lambda2)));    // (A.1) Mid, Claeys+2014
	else                                        lambdaCE = 2.0 * lambda2;			                                        // (A.1) Top, Claeys+2014

	return	lambdaCE;
}

/*
 * Calculata the lambda parameter using the Loveridge prescription
 *
 * Binding energy from detailed models (Loveridge et al. 2011) is given in [E]=erg, so use cgs
 *
 *
 * double CalculateLambdaLoveridge(const double p_EnvMass, const bool p_IsMassLoss)
 *
 * @param   [IN]    p_EnvMass                   Envelope mass (Msol)
 * @param   [IN]    p_IsMassLoss                Boolean indicating whether mass-loss correction should be applied
 * @return                                      Common envelope lambda parameter
 */
double HG_Constituent::CalculateLambdaLoveridge(const double p_EnvMass, const bool p_IsMassLoss) const {

    const double mass        = m_Star->Mass();
    const double mzams       = m_Star->MZAMS();
    const double radius      = m_Star->Radius();
    const double metallicity = m_Star->Metallicity();
    const double coCoreMass  = m_Star->COCoreMass();

    // find closest metallicity covered by Loveridge et al. 2011
    // (see LOVERIDGE_METALLICITY and LOVERIDGE_METALLICITYValue)

    int lMetallicity = 0;
    double minDiff   = std::numeric_limits<double>::max();

    // initialise m_MassCutoffs vector - so we have the right number of entries
    for (int i = 0; i < static_cast<int>(LOVERIDGE_METALLICITY::COUNT); i++) {
        double thisDiff = std::abs(metallicity - std::get<1>(LOVERIDGE_METALLICITY_VALUE[i]));
        if (utils::Compare(thisDiff, minDiff) < 0) {
            lMetallicity = i;
            minDiff      = thisDiff;
        }
    }

    // Determine the evolutionary stage of the star (see LOVERIDGE_GROUP)

    LOVERIDGE_GROUP lGroup;

    if (utils::Compare(mass, LOVERIDGE_LM_HM_CUTOFFS[lMetallicity]) > 0) {                  // mass > low mass / high mass cutoff?
        lGroup = LOVERIDGE_GROUP::HM;                                                       // yes, group is HM - High Mass
    }
    else {                                                                                  // no - low mass
        if (utils::Compare(coCoreMass, 0.0) > 0) {                                          // CO core exists?
            lGroup = LOVERIDGE_GROUP::LMA;                                                  // yes, group is LMA - Low mass on the AGB
        }
        else {                                                                              // no - low mass star on RGB

            // calculate early / late cutoff for low mass RGB stars
            constexpr double deltaM   = 1.0E-5;
                      double cutOff   = 0.0;
                      int    exponent = 0;
            for (auto const& aCoefficient: LOVERIDGE_LM1_LM2_CUTOFFS[lMetallicity]) {
                cutOff += aCoefficient * utils::intPow(std::log10(mass + deltaM), exponent++);
            }

            // set evolutionary stage based on cutoff
            lGroup = utils::Compare(std::log10(radius), cutOff) > 0 ? LOVERIDGE_GROUP::LMR2 : LOVERIDGE_GROUP::LMR1;
        }
    }

    // calculate log10(binding energy)
    constexpr double deltaR           = 1.0E-5;
              double logBindingEnergy = 0.0;
    for (auto const& lCoefficients: LOVERIDGE_COEFFICIENTS[lMetallicity][static_cast<int>(lGroup)]) {
        logBindingEnergy += lCoefficients.alpha_mr * utils::intPow(std::log10(mass), lCoefficients.m) * utils::intPow(std::log10(radius + deltaR), lCoefficients.r);
    }

    double MZAMS_Mass = (mzams - mass) / mzams;                                             // should m_ZAMS really be Mass0() (i.e., account for change in effective mass through mass loss in winds, MS mass transfer?)
    logBindingEnergy *= p_IsMassLoss ? 1.0 + (0.25 * MZAMS_Mass * MZAMS_Mass) : 1.0;        // apply mass-loss correction factor (lambda)

    logBindingEnergy += 33.29866;                                                           // + logBE0
    double bindingEnergy = PPOW(10.0, logBindingEnergy);
    
    double lambda = utils::Compare(bindingEnergy, 0.0) > 0 && utils::Compare(1.0 / bindingEnergy, 0.0) > 0 && utils::Compare(p_EnvMass, MAXIMUM_MASS_LOSS_FRACTION * mass) > 0
            ? (G_CGS * mass * MSOL_TO_G * p_EnvMass * MSOL_TO_G) / (radius * RSOL_TO_AU * AU_TO_CM * bindingEnergy)
            : 1.0;                                                                          // default to 1.0 (usual lambda default) if binding energy is not sensible [sometimes can be infinite if logBindingEnergy is too high] or if envelope mass is too low to reliably evaluate lambda [can be zero]
    return utils::Compare(lambda, 0.0) > 0 ? lambda : 1.0;                                  // final check to avoid returning zero lambda
}


/*
 * CalculateCELambda_Nanjing_Enhanced
 *
 * @brief
 * Calculate the common envelope lambda parameter, per Xu & Li, 2010 (Nanjing - for the university)
 * (https://arxiv.org/abs/1004.4957, v1, 28Apr2010)
 * (https://iopscience.iop.org/article/10.1088/0004-637X/716/1/114)
 *
 * This function good for HG and FGB stars.
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
double HG::CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
                                              const double             p_Radius,
                                              const double             p_CoreMass,
                                              const SizeT              p_MassIndex,
                                              const STELLAR_POPULATION p_StellarPop) const {

    [[maybe_unused]] constexpr SizeT evolStage = 1;                                         // HG evolutionary stage from Xu & Li, 2010

    SizeT coeffsBGidx = 0;                                                                  // Index into coefficients vector
    bool  useLambdas  = false;                                                              // Flag - use lambdas defined in the paper

    ApplyNanjingEnhancedDecisions(
        p_StellarPop == STELLAR_POPULATION::POPULATION_I ? NANJING_DECISIONS_HG_POPI : NANJING_DECISIONS_HG_POPII,
        p_MassIndex, p_Radius,
        useLambdas, coeffsBGidx);

    // Get limits and (defined) lambdas
    const auto& enData   = (p_StellarPop == STELLAR_POPULATION::POPULATION_I) ? NANJING_DATA_HG_POPI : NANJING_DATA_HG_POPII;
    const auto& massData = enData[p_MassIndex];
    Dbl_Dbl_DblT maxBGR  = std::make_tuple(massData.maxB, massData.maxG, massData.maxR);    // {maxB, maxG, maxR}

    double lambdaB;
    double lambdaG;
    if (useLambdas) {                                                                       // Use lambdas defined in the paper?
                                                                                            // Yes
        lambdaB = massData.definedLambdaB;
        lambdaG = massData.definedLambdaG;
    }
    else {                                                                                  // No - calculate lambdas (per paper)

        // Get B & G coefficients vector
        const auto&       coeffsTable = (p_StellarPop == STELLAR_POPULATION::POPULATION_I) ? NANJING_COEFFS_HG_POPI : NANJING_COEFFS_HG_POPII;
        const auto&       BGcoeffs    = coeffsTable[p_MassIndex][coeffsBGidx];
        const DBL_VECTOR& Bcoeffs     = std::get<0>(BGcoeffs);
        const DBL_VECTOR& Gcoeffs     = std::get<1>(BGcoeffs);

        const double rinClamped  = std::min(p_Radius, std::get<2>(maxBGR));                 // clamp to maxR to stay in polynomial-fit domain
        const double envMassFrac = (p_Mass - p_CoreMass) / p_Mass;                          // input for POLY_INVERSE at massIdx==0

        const auto& branchRules = (p_StellarPop == STELLAR_POPULATION::POPULATION_I) ? NANJING_POLY_BRANCH_HG_POPI : NANJING_POLY_BRANCH_HG_POPII;

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
COMPAS_PURE double HG_Constituent::CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const {

    [[maybe_unused]] constexpr SizeT evolStage = 1;                                         // HG evolutionary stage from Xu & Li, 2010
                                           
    SizeT coeffsBGidx = 0;                                                                  // Index into coefficients vector
    SizeT limitBGRidx = 0;                                                                  // Index into limits vector
    SizeT lambdaBGidx = 0;                                                                  // Index into lambdas vector (only used if useLambdas)
    bool  useLambdas  = false;                                                              // Flag - use defined lambda set (set by helper)

    // determine mass index based on p_Mass
    auto it = std::upper_bound(NANJING_MASSES_MIDPOINTS.begin(), NANJING_MASSES_MIDPOINTS.end(), p_Mass);
    const SizeT massIndex = it != NANJING_MASSES_MIDPOINTS.end() ? std::distance(NANJING_MASSES_MIDPOINTS.begin(), it) : NANJING_MASSES_MIDPOINTS.size();

    ApplyNanjingStarTrackDecisions(
        m_Star->Metallicity() > LAMBDA_NANJING_ZLIMIT_STARTRACK ? NANJING_STARTRACK_DECISIONS_HG_HIGHZ : NANJING_STARTRACK_DECISIONS_HG_LOWZ,
        massIndex, p_Radius,
        useLambdas, lambdaBGidx, limitBGRidx, coeffsBGidx);

    // get limits and (defined) lambdas
    const auto& stData     = (m_Star->Metallicity() > LAMBDA_NANJING_ZLIMIT) ? NANJING_DATA_HG_STARTRACK_HIGHZ : NANJING_DATA_HG_STARTRACK_LOWZ;
    const auto& stMassData = stData[massIndex];
    Dbl_DblT maxBG = stMassData.limits[limitBGRidx];                                        // {maxB, maxG}

    double lambdaB;
    double lambdaG;
    if (useLambdas) {                                                                       // Use lambdas defined by StarTrack?
                                                                                            // Yes
        Dbl_DblT lambdaBG = stMassData.lambdas[lambdaBGidx];                                // Defined {lambdaB, lambdaG}

        lambdaB = std::get<0>(lambdaBG);
        lambdaG = std::get<1>(lambdaBG);
    }
    else {                                                                                  // No - calculate lambdas (per StarTrack)

        // Get B & G coefficients vector
        const auto&       coeffsTable = (m_Star->Metallicity() > LAMBDA_NANJING_ZLIMIT) ? NANJING_COEFFS_HG_POPI : NANJING_COEFFS_HG_POPII;
        const auto&       BGcoeffs    = coeffsTable[massIndex][coeffsBGidx];
        const DBL_VECTOR& Bcoeffs     = std::get<0>(BGcoeffs);
        const DBL_VECTOR& Gcoeffs     = std::get<1>(BGcoeffs);

        const double rinClamped  = p_Radius;
        const double envMassFrac = (p_Mass - p_CoreMass) / p_Mass;

        const auto& branchRules = (m_Star->Metallicity() > LAMBDA_NANJING_ZLIMIT) ? NANJING_POLY_BRANCH_HG_POPI : NANJING_POLY_BRANCH_HG_POPII;

        if (!ApplyNanjingPolyBranchOverride(branchRules, massIndex, p_Radius, rinClamped, envMassFrac,
                                            Bcoeffs, Gcoeffs, lambdaB, lambdaG)) {
            // No rule matched -- default 5th-order polynomial
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
    lambdaG = std::min(std::max(0.05, lambdaG), std::min(1.0, std::get<1>(maxBG)));         // Clamp lambda G to [0.05, min(1, maxG)]
    lambdaB = std::max(std::min(lambdaB, std::get<0>(maxBG)), std::max(0.05, lambdaG));     // Clamp lambda B to [max(0.05, lambdaG), maxB]

    // Calculate lambda as some combination of lambdaB and lambdaG by
    // lambda = alpha_th • lambdaB + (1-alpha_th) • lambdaG
    // STARTRACK uses alpha_th = 1/2
    return (OPTIONS->CommonEnvelopeAlphaThermal() * lambdaB) + ((1.0 - OPTIONS->CommonEnvelopeAlphaThermal()) * lambdaG);
}
