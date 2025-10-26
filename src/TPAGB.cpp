#include "TPAGB.h"
#include "HeWD.h"
#include "COWD.h"


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
 * DBL_VECTOR CalculateTimescales_Hurley2000(const double      p_Mass,
 *                                           const double      p_ZetaHurley,
 *                                           const DBL_VECTOR& p_GBparams,
 *                                           const DBL_VECTOR& p_MassCutoffs,
 *                                           const DBL_VECTOR& p_Timescales,
 *                                           const double      p_Alpha3,
 *                                           const DBL_VECTOR& p_aN,
 *                                           const DBL_VECTOR& p_bN) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_ZetaHurley                    Hurley zeta value (log10(Z / ZSOL_HURLEY))
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_MassCutoffs                   Hurley mass cutoffs (Msol)
 * @param       p_tScales                       Hurley timescales (Myr)
 * @param       p_Alpha3                        Hurley alpha3 constant
 * @param       p_aN                            Hurley a(n) coefficients
 * @param       p_bN                            Hurley b(n) coefficients
 * @return                                      Mutated timescales (Myr)
 */
DBL_VECTOR TPAGB::CalculateTimescales_Hurley2000(const double      p_Mass,
                                                 const double      p_ZetaHurley,
                                                 const DBL_VECTOR& p_GBparams,
                                                 const DBL_VECTOR& p_MassCutoffs,
                                                 const DBL_VECTOR& p_tScales,
                                                 const double      p_Alpha3,
                                                 const DBL_VECTOR& p_aN,
                                                 const DBL_VECTOR& p_bN) const {

// #defines for convenience and readability - undefined at end of function
#define GBparams(x) p_GBparams[static_cast<int>(GBP::x)]
#define tScales(x) tScales[static_cast<int>(TIMESCALE::x)]

    double p1   = GBparams(p) - 1.0;
    double q1   = GBparams(q) - 1.0;
    double p1_p = p1 / GBparams(p);
    double q1_q = q1 / GBparams(q);

    double lDU  = CalculateLuminosityGivenCoreMass(GBparams(McDU));

    DBL_VECTOR tScales = p_tScales; // copy given timescales

    // (re)calculate EAGB timescales (Note: EAGB does not recalculate earlier timescales) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    tScales = EAGB::CalculateTimescales_Hurley(p_Mass, p_ZetaHurley, p_GBparams, p_MassCutoffs, tScales, p_Alpha3, p_aN, p_bN);

    tScales[static_cast<int>(TIMESCALE::tP)] = CalculateLifetimeTo2ndDredgeUp(tScales(tinf1_FAGB), tScales(tinf2_FAGB));

    if (lDU > GBparams(Lx)) {
        tScales[static_cast<int>(TIMESCALE::tinf1_SAGB)] = tScales(tinf1_FAGB);
        tScales[static_cast<int>(TIMESCALE::tMx_SAGB)]   = tScales(tMx_FAGB);
        tScales[static_cast<int>(TIMESCALE::tinf2_SAGB)] = tScales(tP) + ((1.0 / (q1 * GBparams(AHHe) * GBparams(B))) * PPOW((GBparams(B) / LDU), q1_q));
    }
    else {
        tScales[static_cast<int>(TIMESCALE::tinf1_SAGB)] = tScales(tP) + ((1.0 / (p1 * GBparams(AHHe) * GBparams(D) )) * PPOW((GBparams(D) / LDU), p1_p));
        tScales[static_cast<int>(TIMESCALE::tMx_SAGB)]   = tScales(tinf1_SAGB) - ((tScales(tinf1_SAGB) - tScales(tP)) * PPOW((LDU / GBparams(Lx)), p1_p));
        tScales[static_cast<int>(TIMESCALE::tinf2_SAGB)] = tScales(tMx_SAGB) + ((1.0 / (q1 * GBparams(AHHe) * GBparams(B))) * PPOW((GBparams(B) / GBparams(Lx)), q1_q));
    }

    // return timescales vector by value - NRVO takes care of performance/efficiency
    return tScales;

#undef timescales
#undef GBparams
}


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
double TPAGB::CalculateCELambda_Dewi() const {

    double lambda3 = std::min(-0.9, 0.58 + (0.75 * log10(m_Mass))) - (0.08 * log10(m_Luminosity));                          // (A.4) Claeys+2014
    double lambda1 = std::max(1.0, std::max(lambda3, -3.5 - (0.75 * log10(m_Mass)) + log10(m_Luminosity)));                 // (A.5) Bottom, Claeys+2014
	double lambda2 = 0.42 * PPOW(m_RZAMS / m_Radius, 0.4);                                                                  // (A.2) Claeys+2014
	double envMass = utils::Compare(m_CoreMass, 0.0) > 0 && utils::Compare(m_Mass, m_CoreMass) > 0 ? m_Mass - m_CoreMass : 0.0;

    double lambdaCE;

         if (utils::Compare(envMass, 1.0) >= 0) lambdaCE = 2.0 * lambda1;                                                   // (A.1) Bottom, Claeys+2014
	else if (utils::Compare(envMass, 0.0) >  0) lambdaCE = 2.0 * (lambda2 + (std::sqrt(envMass) * (lambda1 - lambda2)));         // (A.1) Mid, Claeys+2014
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
double TPAGB::CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
                                                 const double             p_Radius,
                                                 const double             p_CoreMass,
                                                 const size_t             p_MassIndex,
                                                 const STELLAR_POPULATION p_StellarPop) const {

    constexpr size_t evolStage   = 3;                                                       // TPAGB evolutionary stage from Xu & Li, 2010

    size_t           coeffsBGidx = 0;                                                       // index into coefficients array
    bool             useLambdas  = false;                                                   // flag - use lambdas defined in the paper

    if (p_StellarPop == STELLAR_POPULATION::POPULATION_I) {                                 // pop I
        if (p_MassIndex == 10) {
                 if (p_Radius >   0.0 && p_Radius <= 350.0) coeffsBGidx = 0;
            else if (p_Radius > 350.0 && p_Radius <= 600.0) coeffsBGidx = 1;
            else                                            coeffsBGidx = 2;
        }
    }
    else {                                                                                  // pop II
        if (p_MassIndex == 2 && p_Radius > 36.0 && p_Radius < 53.0) useLambdas = true;
    }

    const size_t popIdx = static_cast<size_t>(p_StellarPop);                                // set pop index - pop I or pop II

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
       
        const double Rin = std::min(p_Radius, std::get<2>(maxBGR));                         // clamp to maximum allowed radius (maxR) to prevent exceeding domain of the polynomial fits

        if (p_StellarPop == STELLAR_POPULATION::POPULATION_I && p_MassIndex == 0) {
            const double tmp = 0.1 - (Rin * 3.57E-04);
            lambdaB = tmp;
            lambdaG = tmp;
        }
        else {
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
 * double CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity)
 *
 * @param   [IN]    p_Mass                      Mass
 * @param   [IN]    p_Metallicity               Metallicity
 * 
 * @return                                      Nanjing lambda for use in common envelope
 */
double TPAGB::CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const {

	DBL_VECTOR maxBG    = {};                                                       // [0] = maxB, [1] = maxG
	DBL_VECTOR lambdaBG = {};                                                       // [0] = lambdaB, [1] = lambdaG
	DBL_VECTOR a        = {};                                                       // 0..5 a_coefficients
	DBL_VECTOR b        = {};                                                       // 0..5 b_coefficients

    if (utils::Compare(p_Metallicity, LAMBDA_NANJING_ZLIMIT_STARTRACK) > 0) {                 // Z>0.5 Zsun: popI
        if (utils::Compare(p_Mass, 1.5) < 0) {
            maxBG = { 2.5, 1.5 };
            if (utils::Compare(m_Radius, 200.0) > 0) lambdaBG = { 0.05, 0.05 };
            else  {
                double tmp = 0.1 - (m_Radius * 3.57E-04);
                lambdaBG   = { tmp, tmp };
            }
        }
        else if (utils::Compare(p_Mass, 2.5) < 0) {
            maxBG = { 4.0, 2.0 };
            if (utils::Compare(m_Radius, 340.0) > 0) lambdaBG = { 3.589970, 0.514132 };
            else {
                a = { 0.88954, 0.0098 , -3.1411E-05 , 7.66979E-08,  0.0       , 0.0 };
                b = { 0.48271, 0.00584, -6.22051E-05, 2.41531E-07, -3.1872E-10, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 3.5) < 0) {
            maxBG = { 500.0, 10.0 };
            if (utils::Compare(m_Radius, 400.0) > 0) lambdaBG = { 116.935557, 0.848808 };
            else {
                a = { -0.04669, 0.00764, -4.32726E-05, 9.31942E-08, 0.0        ,  0.0 };
                b = {  0.44889, 0.01102, -6.46629E-05, 5.66857E-09, 7.21818E-10, -1.2201E-12 };
            }
        }
        else if (utils::Compare(p_Mass, 4.5) < 0) {
            maxBG = { 1000.0, 8.0 };
            if (utils::Compare(m_Radius, 410.0) > 0) lambdaBG = { 52.980056, 1.109736 };
            else {
                a = { -0.37322, 0.00943, -3.26033E-05, 5.37823E-08, 0.0, 0.0 };
                b = {  0.13153, 0.00984, -2.89832E-05, 2.63519E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 5.5) < 0) {
            maxBG = { 1000.0, 8.0 };
            if (utils::Compare(m_Radius, 430.0) > 0) lambdaBG = { 109.593522, 1.324248 };
            else {
                a = { -0.80011, 0.00992, -3.03247E-05,  5.26235E-08, 0.0, 0.0 };
                b = { -0.00456, 0.00426,  4.71117E-06, -1.72858E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 6.5) < 0) {
            maxBG = { 25.5, 5.0 };
            if (utils::Compare(m_Radius, 440.0) > 0) lambdaBG = { 16.279603, 1.352166 };
            else {
                a = { -2.7714 ,  0.06467, -4.01537E-04,  7.98466E-07, 0.0, 0.0 };
                b = {  0.23083, -0.00266,  2.21788E-05, -2.35696E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 7.5) < 0) {
            maxBG = { 9.0, 3.0 };
            if (utils::Compare(m_Radius, 420.0) > 0) lambdaBG = { 5.133959, 1.004036 };
            else {
                a = { -0.63266,  0.02054, -1.3646E-04 ,  2.8661E-07 , 0.0, 0.0 };
                b = {  0.26294, -0.00253,  1.32272E-05, -7.12205E-09, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 8.5) < 0) {
            maxBG = { 7.0, 3.0 };
            if (utils::Compare(m_Radius, 490.0) > 0) lambdaBG = { 4.342985, 0.934659 };
            else {
                a = { -0.1288 ,  0.0099 , -6.71455E-05,  1.33568E-07, 0.0, 0.0 };
                b = {  0.26956, -0.00219,  7.97743E-06, -1.53296E-09, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 9.5) < 0) {
            maxBG = { 4.0, 2.0 };
            if (utils::Compare(m_Radius, 530.0) > 0) lambdaBG = { 2.441672, 0.702310 };
            else {
                a = { 1.19804, -0.01961, 1.28222E-04, -3.41278E-07, 3.35614E-10, 0.0 };
                b = { 0.40587, -0.0051 , 2.73866E-05, -5.74476E-08, 4.90218E-11, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 11.0) < 0) {
            maxBG = { 3.0, 1.5 };
            if (utils::Compare(m_Radius, 600.0) > 0) lambdaBG = { 1.842314, 0.593854 };
            else {
                a = { 0.3707 ,  2.67221E-04, -9.86464E-06, 2.26185E-08, 0.0, 0.0 };
                b = { 0.25549, -0.00152    ,  3.35239E-06, 2.24224E-10, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 13.0) < 0) {
            maxBG = { 1.5, 1.0 };
                 if (utils::Compare(m_Radius, 850.0) > 0) lambdaBG = { 0.392470, 0.176660 };
            else if (utils::Compare(m_Radius, 0.0) > 0 && utils::Compare(m_Radius, 350.0) <= 0) {
                a = { 1.28593, -0.02209, 1.79764E-04, -6.21556E-07, 7.59444E-10, 0.0 };
                b = { 0.68544, -0.01394, 1.20845E-04, -4.29071E-07, 5.29169E-10, 0.0 };
            }
            else if (utils::Compare(m_Radius, 350.0) > 0 && utils::Compare(m_Radius, 600.0) <= 0) {
                a = { -11.99537,  0.0992, -2.8981E-04,  3.62751E-07, -1.65585E-10, 0.0 };
                b = {   0.46156, -0.0066,  3.9625E-05, -9.98667E-08, -8.84134E-11, 0.0 };
            }
            else {
                a = { -58.03732, 0.23633, -3.20535E-04, 1.45129E-07, 0.0, 0.0 };
                b = { -15.11672, 0.06331, -8.81542E-05, 4.0982E-08 , 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 15.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 1000.0) > 0) lambdaBG = { 0.414200, 0.189008 };
            else {
                a = { -106.90553, 0.36469, -4.1472E-04 , 1.57349E-07, 0.0, 0.0 };
                b = {  -39.93089, 0.13667, -1.55958E-04, 5.94076E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 18.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 1050.0) > 0) lambdaBG = { 0.2, 0.1 };
            else {
                a = { -154.70559, 0.46718, -4.70169E-04, 1.57773E-07, 0.0, 0.0 };
                b = {  -65.39602, 0.19763, -1.99078E-04, 6.68766E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 35.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 1200.0) > 0) lambdaBG = { 0.05, 0.05 };
            else {
                a = { -260484.85724, 4.26759E+06, -2.33016E+07, 4.24102E+07, 0.0, 0.0 };
                b = { -480055.67991, 7.87484E+06, -4.30546E+07, 7.84699E+07, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 75.0) < 0) {
            maxBG = { 1.0, 0.5 };
            a     = { 0.31321, -7.50384E-04, 5.38545E-07, -1.16946E-10, 0.0, 0.0 };
            b     = { 0.159  , -3.94451E-04, 2.88452E-07, -6.35132E-11, 0.0, 0.0 };
        }
        else {
            maxBG = { 1.0, 0.5 };
            a     = { 0.376 , -0.0018 , 2.81083E-06, -1.67386E-09, 3.35056E-13, 0.0 };
            b     = { 0.2466, -0.00121, 1.89029E-06, -1.12066E-09, 2.2258E-13 , 0.0 };
        }
    }
    else {                                                                          // Z<=0.5 Zsun: popI and popII
        if (utils::Compare(p_Mass, 1.5) < 0) {
            maxBG = { 2.0, 1.5 };
            if (utils::Compare(m_Radius, 160.0) > 0) lambdaBG = { 0.05, 0.05 };
            else {
                a = { 0.24012, -0.01907, 6.09529E-04, -8.17819E-06, 4.83789E-08, -1.04568e-10 };
                b = { 0.15504, -0.01238, 3.96633E-04, -5.3329E-06 , 3.16052E-08, -6.84288e-11 };
            }
        }
        else if (utils::Compare(p_Mass, 2.5) < 0) {
            maxBG = { 4.0, 2.0 };
            if (utils::Compare(m_Radius, 350.0) > 0) lambdaBG = { 2.868539, 0.389991 };
            else {
                a = { 0.5452 ,  0.00212    , 6.42941E-05, -1.46783E-07, 0.0       ,  0.0 };
                b = { 0.30594, -9.58858E-04, 1.12174E-04, -1.04079E-06, 3.4564E-09, -3.91536e-12 };
            }
        }
        else if (utils::Compare(p_Mass, 3.5) < 0) {
            maxBG = { 600.0, 2.0 };
            if (utils::Compare(m_Radius, 400.0) > 0)                                           lambdaBG = { 398.126442, 0.648560 };
            else if (utils::Compare(m_Radius, 36.0) > 0 && utils::Compare(m_Radius, 53.0) < 0) lambdaBG = { 1.0, 1.0 };
            else {
                a = { -0.475  , -0.00328, 1.31101E-04, -6.03669E-07, 8.49549E-10, 0.0 };
                b = {  0.05434,  0.0039 , 9.44609E-06, -3.87278E-08, 0.0        , 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 4.5) < 0) {
            maxBG = { 600.0, 2.0 };
            if (utils::Compare(m_Radius, 410.0) > 0) lambdaBG = { 91.579093, 1.032432 };
            else {
                a = { -0.2106 , -0.01574, 2.01107E-04, -6.90334E-07, 7.92713E-10, 0.0 };
                b = {  0.36779, -0.00991, 1.19411E-04, -3.59574E-07, 3.33957E-10, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 5.5) < 0) {
            maxBG = { 10.0, 3.0 };
            if (utils::Compare(m_Radius, 320.0) > 0) lambdaBG = { 7.618019, 1.257919 };
            else {
                a = { -0.12027,  0.01981, -2.27908E-04,  7.55556E-07, 0.0, 0.0 };
                b = {  0.31252, -0.00527,  3.60348E-05, -3.22445E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 6.5) < 0) {
            maxBG = { 4.0, 1.5 };
            if (utils::Compare(m_Radius, 330.0) > 0) lambdaBG = { 2.390575, 0.772091 };
            else {
                a = { 0.26578,  0.00494, -7.02203E-05, 2.25289E-07, 0.0, 0.0 };
                b = { 0.26802, -0.00248,  6.45229E-06, 1.69609E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 7.5) < 0) {
            maxBG = { 2.5, 1.0 };
            if (utils::Compare(m_Radius, 360.0) > 0) lambdaBG = { 1.878174, 0.646353 };
            else {
                a = { 0.8158 , -0.01633, 1.46552E-04, -5.75308E-07, 8.77711E-10, 0.0 };
                b = { 0.26883, -0.00219, 4.12941E-06,  1.33138E-08, 0.0        , 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 8.5) < 0) {
            maxBG = { 2.0, 1.0 };
            if (utils::Compare(m_Radius, 400.0) > 0) lambdaBG = { 1.517662, 0.553169 };
            else {
                a = { 0.74924, -0.01233, 9.55715E-05, -3.37117E-07, 4.67367E-10, 0.0 };
                b = { 0.25249, -0.00161, 8.35478E-07,  1.25999E-08, 0.0        , 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 9.5) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 440.0) > 0) lambdaBG = { 1.136394, 0.478963 };
            else {
                a = { 0.73147, -0.01076, 7.54308E-05, -2.4114E-07 , 2.95543E-10, 0.0 };
                b = { 0.31951, -0.00392, 2.31815E-05, -6.59418E-08, 7.99575E-11, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 11.0) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 500.0) > 0) lambdaBG = { 1.068300, 0.424706 };
            else {
                a = { -9.26519,  0.08064, -2.30952E-04, 2.21986E-07, 0.0, 0.0 };
                b = {  0.81491, -0.00161, -8.13352E-06, 1.95775E-08, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 13.0) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 600.0) > 0) lambdaBG = { 0.537155, 0.211105 };
            else {
                a = { -51.15252, 0.30238, -5.95397E-04, 3.91798E-07, 0.0, 0.0 };
                b = { -13.44   , 0.08141, -1.641E-04  , 1.106E-07  , 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 15.0) < 0) {
            maxBG = { 1.6, 1.0 };
            if (utils::Compare(m_Radius, 650.0) > 0) lambdaBG = { 0.3, 0.160696 };
            else {
                a = { -140.0   , 0.7126 , -0.00121    , 6.846E-07  , 0.0, 0.0 };
                b = {  -44.1964, 0.22592, -3.85124E-04, 2.19324E-07, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 18.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 750.0) > 0) lambdaBG = { 0.5, 0.204092 };
            else {
                a = { -358.4    , 1.599  , -0.00238   , 1.178E-06  , 0.0, 0.0 };
                b = { -118.13757, 0.52737, -7.8479E-04, 3.89585E-07, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 35.0) < 0) {
            maxBG = { 1.5, 1.0 };
            if (utils::Compare(m_Radius, 900.0) > 0) lambdaBG = { 0.2, 0.107914 };
            else {
                a = { -436.00777, 1.41375, -0.00153    , 5.47573E-07, 0.0, 0.0 };
                b = { -144.53456, 0.46579, -4.99197E-04, 1.78027E-07, 0.0, 0.0 };
            }
        }
        else if (utils::Compare(p_Mass, 75.0) < 0) {
            maxBG = { 20.0, 3.0 };
            a     = { 0.821  , -0.00669, 1.57665E-05, -1.3427E-08 , 3.74204E-12, 0.0 };
            b     = { 0.49287, -0.00439, 1.06766E-05, -9.22015E-09, 2.58926E-12, 0.0 };
        }
        else {
            maxBG = { 4.0, 2.0 };
            a     = { 1.25332, -0.02065, 1.3107E-04 , -3.67006E-07, 4.58792E-10, -2.09069E-13 };
            b     = { 0.81716, -0.01436, 9.31143E-05, -2.6539E-07 , 3.30773E-10, -1.51207E-13 };
        }
    }

    if (lambdaBG.empty()) {                                                         // calculate lambda B & G - not approximated by hand
         if (utils::Compare(p_Metallicity, LAMBDA_NANJING_ZLIMIT_STARTRACK) > 0 &&
            (utils::Compare(p_Mass, 1.5) < 0 || (utils::Compare(p_Mass, 25.0) < 0 && utils::Compare(p_Mass, 18.0) >= 0)) ) {
            double x  = (m_Mass - m_CoreMass) / m_Mass;
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x2 * x2;
            double x5 = x3 * x2;

            double y1 = an[0] + (an[1] * x) + (an[2] * x2) + (an[3] * x3) + (an[4] * x4) + (an[5] * x5);
            double y2 = b[0] + (b[1] * x) + (b[2] * x2) + (b[3] * x3) + (b[4] * x4) + (b[5] * x5);

            lambdaBG = { 1.0 / y1, 1.0 / y2 };
        }
        else if ( (utils::Compare(p_Metallicity, LAMBDA_NANJING_ZLIMIT_STARTRACK) > 0  && utils::Compare(p_Mass, 2.5) >= 0 && utils::Compare(p_Mass, 5.5) < 0) ||
                  (utils::Compare(p_Metallicity, LAMBDA_NANJING_ZLIMIT_STARTRACK) <= 0 && utils::Compare(p_Mass, 2.5) >= 0 && utils::Compare(p_Mass, 4.5) < 0)) {
            double x  = m_Radius;
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x2 * x2;
            double x5 = x3 * x2;

            double y1 = an[0] + (an[1] * x) + (an[2] * x2) + (an[3] * x3) + (an[4] * x4) + (an[5] * x5);
            double y2 = b[0] + (b[1] * x) + (b[2] * x2) + (b[3] * x3) + (b[4] * x4) + (b[5] * x5);

            lambdaBG = { pow(10.0, y1), y2 };
        }
        else {
            double x  = m_Radius;
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x2 * x2;
            double x5 = x3 * x2;

            double y1 = an[0] + (an[1] * x) + (an[2] * x2) + (an[3] * x3) + (an[4] * x4) + (an[5] * x5);
            double y2 = b[0] + (b[1] * x) + (b[2] * x2) + (b[3] * x3) + (b[4] * x4) + (b[5] * x5);

            lambdaBG = { y1, y2 };
        }
    }

    // Limit lambda to some 'reasonable' range
    lambdaBG[0] = std::min(std::max(0.05, lambdaBG[0]), maxBG[0]);                  // clamp lambda B to [0.05, maxB]
    lambdaBG[1] = std::min(std::max(0.05, lambdaBG[1]), maxBG[1]);                  // clamp lambda G to [0.05, maxG]

    // Calculate lambda as some combination of lambda_b and lambda_g by
    // lambda = alpha_th • lambda_b + (1-alpha_th) • lambda_g
    // Note that this is different from STARTRACK
    return (OPTIONS->CommonEnvelopeAlphaThermal() * lambdaBG[0]) + ((1.0 - OPTIONS->CommonEnvelopeAlphaThermal()) * lambdaBG[1]);
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
    return EAGB::CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, p_MHeF, p_BnCoefficients);
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
STELLAR_TYPE TPAGB::ResolveEnvelopeLoss(bool p_Force) {
#define gbParams(x) m_GBparams[static_cast<int>(GBP::x)]    // for convenience and readability - undefined at end of function

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
double TPAGB::ChooseTimestep(const double p_Time) const {

    double dtk = utils::Compare(p_Time, timescales(tMx_SAGB)) <= 0
                    ? 0.02 * (timescales(tinf1_SAGB) - p_Time)
                    : 0.02 * (timescales(tinf2_SAGB) - p_Time);

    double dte = 5.0E-3;

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);
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
    bool isCCSN = utils::Compare(m_COCoreMass, CalculateCoreMassAtSN_Static(MCH, m_GBparams[static_cast<int>(GBP::McBAGB)])) >= 0 &&
        utils::Compare(snMass, OPTIONS->MCBUR1()) >= 0 && utils::Compare(m_COCoreMass, m_Mass) < 0;
    bool isECSN = utils::Compare(snMass, MCBUR2) < 0 && (!m_MassTransferDonorHistory.empty() || OPTIONS->AllowNonStrippedECSN()) &&
        utils::Compare(m_COCoreMass, CalculateCoreMassAtSN_Static(MECS, m_GBparams[static_cast<int>(GBP::McBAGB)])) >= 0 &&
        utils::Compare(snMass, OPTIONS->MCBUR1()) >= 0 && utils::Compare(m_COCoreMass, m_Mass) < 0;
    return isCCSN || isECSN;
}
