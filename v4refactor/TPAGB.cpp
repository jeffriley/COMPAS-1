#include "TPAGB.h"
#include "Star.h"
#include "WhiteDwarfs.h"
#include "HeWD.h"
#include "COWD.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////





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
double TPAGB_Constituent::CalculateCELambda_Dewi() const {

    const double mass    = m_Star->Mass();
    const double radius  = m_Star->Radius();
    const double rZAMS   = m_Star->RZAMS();
    const double lum     = m_Star->Luminosity();
    const double cmass   = m_Star->CoreMass();

    const double log10mass = std::log10(mass);

    double lambda3 = std::min(0.9, 0.58 + (0.75 * log10mass)) - (0.08 * std::log10(lum));                                   // (A.4) Claeys+2014 with corrected typo (see e.g. Appendix E.1 in Marchant+2021)
    double lambda1 = std::max(1.0, std::max(lambda3, -3.5 - (0.75 * log10mass) + std::log10(lum)));                         // (A.5) Bottom, Claeys+2014
	double lambda2 = 0.42 * PPOW(rZAMS / radius, 0.4);                                                                       // (A.2) Claeys+2014
	double envMass = utils::Compare(cmass, 0.0) > 0 && utils::Compare(mass, cmass) > 0 ? mass - cmass : 0.0;

    double lambdaCE;

         if (utils::Compare(envMass, 1.0) >= 0) lambdaCE = 2.0 * lambda1;                                                   // (A.1) Bottom, Claeys+2014
	else if (utils::Compare(envMass, 0.0) >  0) lambdaCE = 2.0 * (lambda2 + (std::sqrt(envMass) * (lambda1 - lambda2)));    // (A.1) Mid, Claeys+2014
	else                                        lambdaCE = 2.0 * lambda2;                                                   // (A.1) Top, Claeys+2014

	return	lambdaCE;
}


/*
 * CalculateCELambda_Nanjing_Enhanced
 *
 * @brief
 * Calculate the common envelope lambda parameter, per Xu & Li, 2010 (Nanjing - for the university)
 * (https://arxiv.org/abs/1004.4957, v1, 28Apr2010)
 * (https://iopscience.iop.org/article/10.1088/0004-637X/716/1/114)
 *
 * This function good for TPGB stars.
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
double TPAGB::CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
                                                 const double             p_Radius,
                                                 const double             p_CoreMass,
                                                 const SizeT              p_MassIndex,
                                                 const STELLAR_POPULATION p_StellarPop) const {

    [[maybe_unused]] constexpr SizeT evolStage   = 3;                                                        // TPAGB evolutionary stage from Xu & Li, 2010

    SizeT coeffsBGidx = 0;                                                                  // index into coefficients array
    bool useLambdas   = false;                                                              // flag - use lambdas defined in the paper

    ApplyNanjingEnhancedDecisions(
        p_StellarPop == STELLAR_POPULATION::POPULATION_I ? NANJING_DECISIONS_TPAGB_POPI : NANJING_DECISIONS_TPAGB_POPII,
        p_MassIndex, p_Radius,
        useLambdas, coeffsBGidx);

    [[maybe_unused]] const SizeT popIdx = static_cast<SizeT>(p_StellarPop);                                  // set pop index - pop I or pop II

    // get limits and (defined) lambdas
    const auto& enData   = (p_StellarPop == STELLAR_POPULATION::POPULATION_I) ? NANJING_DATA_TPAGB_POPI : NANJING_DATA_TPAGB_POPII;
    const auto& massData = enData[p_MassIndex];
    std::tuple<double, double, double> maxBGR = std::make_tuple(massData.maxB, massData.maxG, massData.maxR);   // {maxB, maxG, maxR}

    double lambdaB;
    double lambdaG;
    if (useLambdas) {                                                                       // use lambdas defined in the paper?
                                                                                            // yes
        lambdaB = massData.definedLambdaB;
        lambdaG = massData.definedLambdaG;
    }
    else {                                                                                  // no - calculate lambdas (per paper)

        // get B & G coefficients vector
        const auto&       coeffsTable = (p_StellarPop == STELLAR_POPULATION::POPULATION_I) ? NANJING_COEFFS_TPAGB_POPI : NANJING_COEFFS_TPAGB_POPII;
        const auto&       BGcoeffs    = coeffsTable[p_MassIndex][coeffsBGidx];
        const DBL_VECTOR& Bcoeffs     = std::get<0>(BGcoeffs);
        const DBL_VECTOR& Gcoeffs     = std::get<1>(BGcoeffs);

        const double rinClamped  = std::min(p_Radius, std::get<2>(maxBGR));
        const double envMassFrac = (p_Mass - p_CoreMass) / p_Mass;                          // supplied for helper-signature uniformity

        const auto& branchRules = (p_StellarPop == STELLAR_POPULATION::POPULATION_I) ? NANJING_POLY_BRANCH_TPAGB_POPI : NANJING_POLY_BRANCH_TPAGB_POPII;

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
    lambdaG = std::min(std::max(0.05, lambdaG), std::min(1.0, std::get<1>(maxBGR)));        // clamp lambda G to [0.05, min(1, maxG)]
    lambdaB = std::max(std::min(lambdaB, std::get<0>(maxBGR)), std::max(0.05, lambdaG));    // clamp lambda B to [max(0.05, lambdaG), maxB]

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
 * This function good for TPAGB stars.
 *
 *
 * double CalculateLambdaNanjingStarTrack(const double p_Mass)
 *
 * @param   [IN]    p_Mass                      Mass
 * 
 * @return                                      Nanjing lambda for use in common envelope
 */
COMPAS_PURE double TPAGB_Constituent::CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const {


    const bool  isHighZ = utils::Compare(m_Star->Metallicity(), LAMBDA_NANJING_ZLIMIT_STARTRACK) > 0;
    const auto& bins    = isHighZ ? NANJING_STARTRACK_BINS_TPAGB_POPI : NANJING_STARTRACK_BINS_TPAGB_POPII;
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



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               LUMINOSITY FUNCTIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                RADIUS CALCULATIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate radius on the Thermally Pulsing Asymptotic Giant Branch
 *
 * Uses EAGB calculation - Hurley et al. 2000, eq 74
 *
 *
 * double CalculateRadiusOnPhase_Static(const double      p_Mass,
 *                                      const double      p_Luminosity,
 *                                      const double      p_MHeF,
 *                                      const DBL_VECTOR &p_BnCoefficients)
 *
 * @param   [IN]    p_Mass                      Mass in Msol
 * @param   [IN]    p_Luminosity                Luminosity in Lsol
 * @param   [IN]    p_MHeF                      Maximum initial mass for which helium ignites degenerately in a Helium Flash
 * @param   [IN]    p_BnCoefficients            b(n) coefficients
 * @return                                      Radius on the Thermally Pulsing Asymptotic Giant Branch in Rsol
 *
 * p_MHeF, p_MFGB and p_BnCoefficients passed as parameters so function can be declared static
 */
double TPAGB::CalculateRadiusOnPhase_Static(const double      p_Mass,
                                            const double      p_Luminosity,
                                            const double      p_MHeF,
                                            const DBL_VECTOR &p_BnCoefficients) {
    return EAGB::CalculateRadius_Hurley2000_Static(p_Mass, p_Luminosity, p_BnCoefficients, p_MHeF);
}




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
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
 *
 * STELLAR_TYPE ResolveEnvelopeLoss(bool p_Force)
 *
 * @param   [IN]    p_Force                     Boolean to indicate whether the resolution of the loss of the envelope should be performed
 *                                              without checking the precondition(s).
 *                                              Default is false.
 *
 * @return                                      Stellar Type to which star should evolve after losing envelope
 */
STELLAR_TYPE TPAGB::ResolveEnvelopeLoss(bool p_Force) {
#define gbParams(x) m_GBParams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

    STELLAR_TYPE stellarType = m_StellarType;

    if (ShouldEnvelopeBeExpelledByPulsations()) m_EnvelopeJustExpelledByPulsations = true;

    if (p_Force || (utils::Compare(m_CoreMass, m_Mass)) >= 0 || m_EnvelopeJustExpelledByPulsations) {   // envelope loss
                
        m_Mass       = std::min(m_CoreMass, m_Mass);
        m_CoreMass   = m_Mass;
        m_HeCoreMass = m_Mass;
        m_COCoreMass = m_Mass;
        m_Mass0      = m_Mass;
        m_Radius     = COWD::CalculateRadiusOnPhase_Static(m_Mass);
        m_Age        = 0.0;
        stellarType  = (utils::Compare(gbParams(McBAGB), OPTIONS->MCBUR1() ) < 0) ? STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF : STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF;
    }

    return stellarType;

#undef gbParams
}



/*
 * Determine if star should continue evolution as a Supernova
 *
 *
 * bool IsSupernova()
 *
 * @return                                      Boolean flag: true if star has gone Supernova, false if not
 */
bool TPAGB::IsSupernova() const {
    double snMass = CalculateInitialSupernovaMass();
    bool isCCSN = utils::Compare(COCoreMass(), CalculateCoreMassAtSN(MCH, CurrentState().HurleyGBParams(HURLEY_GB_PARAMETERS::McBAGB))) >= 0 &&
                  utils::Compare(snMass, OPTIONS->MCBUR1()) >= 0 && utils::Compare(COCoreMass(), Mass()) < 0;
    bool isECSN = utils::Compare(snMass, MCBUR2) < 0 && (HasMTdonorHistory() || OPTIONS->AllowNonStrippedECSN()) &&
                  utils::Compare(COCoreMass(), CalculateCoreMassAtSN(MECS, CurrentState().HurleyGBParams(HURLEY_GB_PARAMETERS::McBAGB))) >= 0 &&
                  utils::Compare(snMass, OPTIONS->MCBUR1()) >= 0 && utils::Compare(COCoreMass(), Mass()) < 0;
    return isCCSN || isECSN;
}


GNU_CONST double TPAGB::CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass) {
    return WhiteDwarfs::CalculateRadius_Marsh2004(p_CoreMass);
}

