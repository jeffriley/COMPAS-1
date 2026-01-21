#include "HG.h"
#include "HeMS.h"
#include "HeWD.h"


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
 * This function ostensibly relies on the value of the ZAMS radius of the star,
 * and should not be used if the ZAMS radius is not known.
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
    double RTMS;  
    if ((OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK) && (utils::Compare(m_MZAMS, BRCEK_LOWER_MASS_LIMIT) >= 0))
        // p_Mass generally has the value of m_Mass0, but since m_Mass is used for radius calculations on the MS and m_Mass0
        // is updated to a new value when BRCEK prescription is used, we need to use m_Mass here to keep radius continuous
        RTMS = MainSequence::CalculateRadiusAtPhaseEnd(m_Mass, p_RZAMS);
    else
        RTMS = MainSequence::CalculateRadiusAtPhaseEnd(p_Mass, p_RZAMS);
    // FIX THIS - REMOVE BRCEK CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



    const double rTAMS = MainSequence::CalculateRadiusAtPhaseEnd(p_Mass, p_RZAMS);
    const double rGB   = GiantBranch::CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, p_bN);

    double rx = rGB;                                                                                                   // rx in Hurley SSE Fortran code

    if (p_Mass > p_MFGB) {                                                                      // mass above threshold for He ignition?
                                                                                                                        // yes
        // rMinHe is Hurley et al. 2000, eq 55 - first part (M >= MHeF)
        const double Mb28   = PPOW(p_Mass, p_bN[28]);                                                                   // pow() is slow - do it once only
        const double rMinHe = ((p_bN[24] * p_Mass) + (PPOW((p_bN[25] * p_Mass), p_bN[26]) * Mb28)) / (p_bN[27] + Mb28); // rmin in Hurley SSE Fortran code
        const double lum    = GiantBranch::CalculateLuminosityAtHeI_Hurley2000_Static(p_Mass);

        // mass is used here in the Hurley SSE Fortran code, mass0 everywhere else
        // ry in Hurley SSE Fortran code
        const double ry = EAGB::CalculateRadiusOnPhase_Static(p_Mass, lum, p_MHeF, p_bN);

        // calculate radius at He ignition for MFGB < p_Mass < HM
        // Hurley et al. 2000, eq 50
        
        rx = std::min(rMinHe, ry);
        
        if (p_Mass < HIGH_MASS_THRESHOLD) {
            rx = rMinHe * PPOW(rGB / rMinHe, log10(p_Mass / HIGH_MASS_THRESHOLD) / log10(p_MFGB / HIGH_MASS_THRESHOLD));
        }

        // this piece of code resets rx if the blue loop is relatively short
        // see Hurley SSE Fortran code, function tblf() in `zfuncs1.f`
        // JR: I suspect this is where the check came from in the dicussion re blue loop in CHeB::CalculateTimescales() - I don't like this much... **Ilya**
        const double r1 = std::max(1.0 - rMinHe / ry, 1.0E-12); 

        double tblf = (1.0 - b[47]) * PPOW(p_Mass, b[48]) * PPOW(r1, b[49]);                                        // calculate blue-loop fraction of He-burning
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
/*
 * DetermineEnvelopeType
 *
 * @brief
 * Determine the star's envelope type, based on the user-specified ENVELOPE_STATE_PRESCRIPTION.
 * 
 *
 * ENVELOPE DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Envelope type (ENVELOPE::{RADIATIVE, CONVECTIVE, REMNANT})
 */
COMPAS_PURE ENVELOPE HG::DetermineEnvelopeType(const double p_Mass, const double p_Temperature, const double p_CoreMass) const {
 
    ENVELOPE envType;
    
    switch (OPTIONS->EnvelopeStatePrescription()) {                                 // which envelope prescription?
            
        case ENVELOPE_STATE_PRESCRIPTION::CONVECTIVE_MASS_FRACTION:                 // CONVECTIVE_MASS_FRACTION
            // envelope is convective when the convective mass exceeds specified fraction of the envelope mass
            double mEnv;
            std::tie(mEnv, std::ignore) = CalculateConvectiveEnvelopeMass();
            envType = mEnv / (p_Mass - p_CoreMass) > OPTIONS->ConvectiveEnvelopeMassThreshold() ? ENVELOPE::CONVECTIVE : ENVELOPE::RADIATIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE:                        // FIXED_TEMPERATURE
            // envelope is radiative if temperature exceeds specified threshold, otherwise convective
            envType = (p_Temperature * TSOL) > OPTIONS->ConvectiveEnvelopeTemperatureThreshold() ? ENVELOPE::RADIATIVE : ENVELOPE::CONVECTIVE;
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::HURLEY:                                   // HURLEY
            // Hurley et al. 2002 eqs 39 & 40, and Hurley et al. 2000 end of section 7.2 describe gradual
            // growth of convective envelope over HG, but we approximate it as already convective here
            envType = ENVELOPE::CONVECTIVE;                                         // always convective
            break;
            
        case ENVELOPE_STATE_PRESCRIPTION::LEGACY:                                   // COMPAS LEGACY
            envType = ENVELOPE::RADIATIVE;                                          // always radiative
            break;

        default:                                                                    // unknown prescription
            // the only way this can happen is if someone added an ENVELOPE_STATE_PRESCRIPTION and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_ENVELOPE_STATE_PRESCRIPTION);                // throw error               
    }
    
    return envType;
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
double HG::ChooseTimestep(const double p_Time) const {

    double dtk = 0.05 * (timescales(tBGB) - timescales(tMS));
    double dte = timescales(tBGB) - p_Time;    

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
 * STELLAR_TYPE ResolveEnvelopeLoss(bool p_Force)
 *
 * @param   [IN]    p_Force                     Boolean to indicate whether the resolution of the loss of the envelope should be performed
 *                                              without checking the precondition(s).
 *                                              Default is false.
 *
 * @return                                      Stellar Type to which star should evolve after losing envelope
 */
STELLAR_TYPE HG::ResolveEnvelopeLoss(bool p_Force) {

    STELLAR_TYPE stellarType = m_StellarType;

    if (p_Force || utils::Compare(m_CoreMass, m_Mass) >= 0) {                   // envelope loss

        m_Mass = std::min(m_CoreMass, m_Mass);

        if (utils::Compare(m_Mass0, massCutoffs(MHeF)) < 0) {                   // star evolves to Helium White Dwarf

            stellarType  = STELLAR_TYPE::HELIUM_WHITE_DWARF;

            m_Radius     = WhiteDwarfs::CalculateRadius_Hurley2000_Static(m_Mass);
            m_Age        = 0.0;                                                 // see Hurley et al. 2000, discussion after eq 76
        }
        else {                                                                  // star evolves to Zero age Naked Helium Main Star

            stellarType  = STELLAR_TYPE::NAKED_HELIUM_STAR_MS;

            m_Mass0      = m_Mass;
            m_Radius     = HeMS::CalculateRadiusAtZAHeMS_Hurley2000_Static(m_Mass);          
            m_Luminosity = HeMS::CalculateLuminosityAtZAHeMS_Hurley2000_Static(m_Mass);
            m_Age        = 0.0;                                                 // can't use Hurley et al. 2000, eq 76 here - timescales(tHe) not calculated yet
        }
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
STELLAR_TYPE HG::EvolveToNextPhase() {

    STELLAR_TYPE stellarType;

    if (utils::Compare(m_Mass0, massCutoffs(MFGB)) < 0) {
        stellarType = STELLAR_TYPE::FIRST_GIANT_BRANCH;
    }
    else {
        stellarType = STELLAR_TYPE::CORE_HELIUM_BURNING;
    }    

    return stellarType;
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
double HG::CalculateCELambda_Dewi() const {

	double lambda1 = std::min(0.80, (3.0 / (2.4 + PPOW(m_Mass,-3.0 / 2.0))) - (0.15 * log10(m_Luminosity)));                // (A.3) Claeys+2014
	double lambda2 = 0.42 * PPOW(m_RZAMS / m_Radius, 0.4);                                                                  // (A.2) Claeys+2014
	double envMass = utils::Compare(m_CoreMass, 0.0) > 0 && utils::Compare(m_Mass, m_CoreMass) > 0 ? m_Mass - m_CoreMass : 0.0;

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
double HG::CalculateLambdaLoveridge(const double p_EnvMass, const bool p_IsMassLoss) const {
    
    // find closest metallicity covered by Loveridge et al. 2011
    // (see LOVERIDGE_METALLICITY and LOVERIDGE_METALLICITYValue)

    int lMetallicity = 0;
    double minDiff   = std::numeric_limits<double>::max();

    // initialise m_MassCutoffs vector - so we have the right number of entries
    for (int i = 0; i < static_cast<int>(LOVERIDGE_METALLICITY::COUNT); i++) {
        double thisDiff = std::abs(m_Metallicity - std::get<1>(LOVERIDGE_METALLICITY_VALUE[i]));
        if (utils::Compare(thisDiff, minDiff) < 0) {
            lMetallicity = i;
            minDiff      = thisDiff;
        }
    }

    // Determine the evolutionary stage of the star (see LOVERIDGE_GROUP)

    LOVERIDGE_GROUP lGroup;

    if (utils::Compare(m_Mass, LOVERIDGE_LM_HM_CUTOFFS[lMetallicity]) > 0) {                // mass > low mass / high mass cutoff?
        lGroup = LOVERIDGE_GROUP::HM;                                                       // yes, group is HM - High Mass
    }
    else {                                                                                  // no - low mass
        if (utils::Compare(m_COCoreMass, 0.0) > 0) {                                        // CO core exists?
            lGroup = LOVERIDGE_GROUP::LMA;                                                  // yes, group is LMA - Low mass on the AGB
        }
        else {                                                                              // no - low mass star on RGB

            // calculate early / late cutoff for low mass RGB stars
            constexpr double deltaM   = 1.0E-5;
                      double cutOff   = 0.0;
                      int    exponent = 0;
            for (auto const& aCoefficient: LOVERIDGE_LM1_LM2_CUTOFFS[lMetallicity]) {
                cutOff += aCoefficient * utils::intPow(log10(m_Mass + deltaM), exponent++);
            }

            // set evolutionary stage based on cutoff
            lGroup = utils::Compare(log10(m_Radius), cutOff) > 0 ? LOVERIDGE_GROUP::LMR2 : LOVERIDGE_GROUP::LMR1;
        }
    }

    // calculate log10(binding energy)
    constexpr double deltaR           = 1.0E-5;
              double logBindingEnergy = 0.0;
    for (auto const& lCoefficients: LOVERIDGE_COEFFICIENTS[lMetallicity][static_cast<int>(lGroup)]) {
        logBindingEnergy += lCoefficients.alpha_mr * utils::intPow(log10(m_Mass), lCoefficients.m) * utils::intPow(log10(m_Radius + deltaR), lCoefficients.r);
    }

    double MZAMS_Mass = (m_MZAMS - m_Mass) / m_MZAMS;                                       // should m_ZAMS really be m_Mass0 (i.e., account for change in effective mass through mass loss in winds, MS mass transfer?)
    logBindingEnergy *= p_IsMassLoss ? 1.0 + (0.25 * MZAMS_Mass * MZAMS_Mass) : 1.0;        // apply mass-loss correction factor (lambda)

    logBindingEnergy += 33.29866;                                                           // + logBE0
    double bindingEnergy = PPOW(10.0, logBindingEnergy);
    
    double lambda = utils::Compare(bindingEnergy, 0.0) > 0 && utils::Compare(1.0 / bindingEnergy, 0.0) > 0 && utils::Compare(p_EnvMass, MAXIMUM_MASS_LOSS_FRACTION * m_Mass) > 0
            ? (G_CGS * m_Mass * MSOL_TO_G * p_EnvMass * MSOL_TO_G) / (m_Radius * RSOL_TO_AU * AU_TO_CM * bindingEnergy)
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
double HG::CalculateCELambda_Nanjing_Enhanced(const double             p_Mass,
                                              const double             p_Radius,
                                              const double             p_CoreMass,
                                              const size_t             p_MassIndex,
                                              const STELLAR_POPULATION p_StellarPop) const {

    constexpr size_t evolStage = 1;                                                         // HG evolutionary stage from Xu & Li, 2010

    size_t coeffsBGidx = 0;                                                                 // index into coefficients vector
    bool   useLambdas  = false;                                                             // flag - use lambdas defined in the paper

    if (p_StellarPop == STELLAR_POPULATION::POPULATION_I) {                                 // pop I
             if (p_MassIndex ==  1 && p_Radius >   8.5 && p_Radius <  60.0) useLambdas = true;
        else if (p_MassIndex == 11 && p_Radius > 190.0 && p_Radius < 600.0) useLambdas = true;   
        else if (p_MassIndex == 12 && p_Radius > 120.0 && p_Radius < 170.0) useLambdas = true;   
        else if (p_MassIndex == 10) {
                 if (p_Radius >   0.0 && p_Radius <= 350.0) coeffsBGidx = 0;
            else if (p_Radius > 350.0 && p_Radius <= 600.0) coeffsBGidx = 1;
            else                                            coeffsBGidx = 2;
        }
    }
    else {                                                                                  // pop II
        if (p_MassIndex ==  1 && p_Radius >  22.0 && p_Radius <  87.0) useLambdas = true;
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
        
        double Rin = std::min(p_Radius, std::get<2>(maxBGR));                               // clamp to maximum allowed radius (maxR) to prevent exceeding domain of the polynomial fits

        if (p_StellarPop == STELLAR_POPULATION::POPULATION_I && p_MassIndex == 0 && Rin > 2.7) {
            lambdaB = 2.33 - (Rin * 9.18E-03);
            lambdaG = 1.12 - (Rin * 4.59E-03);
        }
        else if (p_StellarPop == STELLAR_POPULATION::POPULATION_I && p_MassIndex == 13) {
            lambdaB = 1.2 * exp(-Rin / 90.0);
            lambdaG = 0.55 * exp(-Rin / 160.0);
        }
        else if (p_StellarPop == STELLAR_POPULATION::POPULATION_II && p_MassIndex == 0 && Rin > 12.0) {
            lambdaB = 1.8 * exp(-Rin / 80.0);
            lambdaG = exp(-Rin / 45.0);
        }
        else if (p_StellarPop == STELLAR_POPULATION::POPULATION_II && p_MassIndex == 9) {
            const double tmp = exp(-Rin / 35.0);
            lambdaB = 1.75 * tmp;
            lambdaG = 0.9 * tmp;
        }
        else {
            if (p_StellarPop == STELLAR_POPULATION::POPULATION_I && p_MassIndex == 0) Rin = (p_Mass - p_CoreMass) / p_Mass;
            
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
double HG::CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Radius, const double p_CoreMass) const {

    constexpr size_t evolStage = 1;                                                         // HG evolutionary stage from Xu & Li, 2010
                                           
    size_t coeffsBGidx = 0;                                                                 // index into coefficients vector
    size_t limitBGRidx = 0;                                                                 // index into limits vector
    size_t lambdaBGidx = -1;                                                                // index into lambdas vector: -ve indicates calculate lambdas

    // determine mass index based on p_Mass
    auto it = std::upper_bound(NANJING_MASSES_MIDPOINTS.begin(), NANJING_MASSES_MIDPOINTS.end(), p_Mass);
    const size_t massIndex = it != arr.end() ? std::distance(NANJING_MASSES_MIDPOINTS.begin(), it) : NANJING_MASSES_MIDPOINTS.size();

    if (GLOBALS->Metallicity() > LAMBDA_NANJING_ZLIMIT_STARTRACK) {                // Z > LAMBDA_NANJING_ZLIMIT_STARTRACK?
                                                                                            // yes
             if (massIndex == 0 && p_Radius > 200.0) lambdaBGidx = 0;
        else if (massIndex == 1 && p_Radius > 340.0) lambdaBGidx = 0;
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
        else if (massIndex == 7 && p_Radius > 490.0) lambdaBGidx = 0;
        else if (massIndex == 8 && p_Radius > 530.0) lambdaBGidx = 0;
        else if (massIndex == 9) {
            if (p_Radius > 600.0) lambdaBGidx = 0;
            else                  limitBGRidx = 1;
        }
        else if (p_MassIndex == 10) {
                 if (p_Radius > 850.0)                      lambdaBGidx = 0;
            else if (p_Radius >   0.0 && p_Radius <= 350.0) coeffsBGidx = 0;
            else if (p_Radius > 350.0 && p_Radius <= 600.0) coeffsBGidx = 1;
            else                                            coeffsBGidx = 2;
        }
        else if (massIndex == 11) {
                 if (p_Radius > 1000.0)                    lambdaBGidx = 0;
            else if (p_Radius > 190.0 && p_Radius < 600.0) lambdaBGidx = 1;
        }
        else if (massIndex == 12) {
                 if (p_Radius > 1050.0)                    lambdaBGidx = 0;
            else if (p_Radius > 120.0 && p_Radius < 170.0) lambdaBGidx = 1;
        }
        else if (massIndex == 13 && p_Radius > 1200.0) lambdaBGidx = 0;
    }
    else {                                                                                  // no - Z <= LAMBDA_NANJING_ZLIMIT_STARTRACK
             if (massIndex == 0 && p_Radius > 160.0) lambdaBGidx = 0;
        else if (massIndex == 1) {
                 if (p_Radius > 350.0)                   lambdaBGidx = 0;
            else if (p_Radius > 22.0 && p_Raidus < 87.0) lambdaBGidx = 1;
        }
        else if (massIndex ==  2 && p_Radius > 400.0) lambdaBGidx = 0;
        else if (massIndex ==  3 && p_Radius > 410.0) lambdaBGidx = 0;
        else if (massIndex ==  4 && p_Radius > 320.0) lambdaBGidx = 0;
        else if (massIndex ==  5 && p_Radius > 330.0) lambdaBGidx = 0;
        else if (massIndex ==  6 && p_Radius > 360.0) lambdaBGidx = 0;
        else if (massIndex ==  7 && p_Radius > 400.0) lambdaBGidx = 0;
        else if (massIndex ==  8 && p_Radius > 440.0) lambdaBGidx = 0;
        else if (massIndex ==  9 && p_Radius > 500.0) lambdaBGidx = 0;
        else if (massIndex == 10 && p_Radius > 600.0) lambdaBGidx = 0;
        else if (massIndex == 11 && p_Radius > 650.0) lambdaBGidx = 0;
        else if (massIndex == 12 && p_Radius > 750.0) lambdaBGidx = 0;
        else if (massIndex == 13 && p_Radius > 900.0) lambdaBGidx = 0;
    }

    // get limits and (defined) lambdas
    NANJING_Z_LIMITS_LAMBDAS                              ZlimitsLambdas = GLOBALS->Metallicity() < LAMBDA_NANJING_ZLIMIT ? std::get<0>(NANJING_LIMITS_LAMBDAS_STARTRACK) : std::get<1>(NANJING_LIMITS_LAMBDAS_STARTRACK);
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
        NANJING_POP_COEFFICIENTS                                       ZCoeffs         = GLOBALS->Metallicity() < LAMBDA_NANJING_ZLIMIT ? std::get<0>(evolStageCoeffs) : std::get<1>(evolStageCoeffs);
        std::tuple<DBL_VECTOR, DBL_VECTOR>                             BGcoeffs        = ZCoeffs[p_MassIndex][coeffsBGidx];

        DBL_VECTOR Bcoeffs = std::get<0>(BGcoeffs);
        DBL_VECTOR Gcoeffs = std::get<1>(BGcoeffs);
        
        double Rin = p_Radius;

        if (GLOBALS->Metallicity() < LAMBDA_NANJING_ZLIMIT && p_MassIndex == 0 && p_Radius > 2.7) {
            lambdaB = 2.33 - (Rin * 9.18E-03);
            lambdaG = 1.12 - (Rin * 4.59E-03);
        }
        else if (GLOBALS->Metallicity() < LAMBDA_NANJING_ZLIMIT && p_MassIndex == 13) {
            lambdaB = 1.2 * exp(-Rin / 90.0);
            lambdaG = 0.55 * exp(-Rin / 160.0);
        }
        else if (GLOBALS->Metallicity() >= LAMBDA_NANJING_ZLIMIT && p_MassIndex == 0 && p_Radius > 12.0) {
            lambdaB = 1.8 * exp(-Rin / 80.0);
            lambdaG = exp(-Rin / 45.0);
        }
        else if (GLOBALS->Metallicity() >= LAMBDA_NANJING_ZLIMITI && p_MassIndex == 9) {
            const double tmp = exp(-Rin / 35.0);
            lambdaB = 1.75 * tmp;
            lambdaG = 0.9 * tmp;
        }
        else {
            if (GLOBALS->Metallicity() < LAMBDA_NANJING_ZLIMIT && p_MassIndex == 0) Rin = (p_Mass - p_CoreMass) / p_Mass;
            
            const double Rin2 = Rin  * Rin;
            const double Rin3 = Rin  * Rin2;
            const double Rin4 = Rin2 * Rin2;
            const double Rin5 = Rin2 * Rin3;

            lambdaB = Bcoeffs[0] + (Bcoeffs[1] * Rin) + (Bcoeffs[2] * Rin2) + (Bcoeffs[3] * Rin3) + (Bcoeffs[4] * Rin4) + (Bcoeffs[5] * Rin5);
            lambdaG = Gcoeffs[0] + (Gcoeffs[1] * Rin) + (Gcoeffs[2] * Rin2) + (Gcoeffs[3] * Rin3) + (Gcoeffs[4] * Rin4) + (Gcoeffs[5] * Rin5);

            if (GLOBALS->Metallicity() < LAMBDA_NANJING_ZLIMIT && p_MassIndex == 0) {
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
