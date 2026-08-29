#include "HeMS.h"
#include "Star.h"
#include "HeWD.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////










///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 RADIUS FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * Calculate convective core radius
 *
 * Assume equal to total radius at start (for continuity with stripped CHeB or HG star), continuity with HeHG at end of phase, linear growth
 *
 *
 * double CalculateConvectiveCoreRadius()
 *
 * @return                                      Convective core radius (solar radii)
 */
double HeMS::CalculateConvectiveCoreRadius(const double p_Radius, const double p_Tau) const {

    // We need core radius at end of phase, which is just the core radius at the start of the HeHG phase.
    // Since we are on the He main sequence here, we can clone this object as an HeHG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HeHG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.     
    std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP, OBJECT_PERSISTENCE::EPHEMERAL, true);
    HeHG* hgClone = dynamic_cast<HeHG*>(clone.get());
    double finalConvectiveCoreRadius = HeHG::CalculateConvectiveCoreRadius_Static(hgClone->Radius(), hgClone->Tau(), hgClone->CoreMass());

    double initialConvectiveCoreRadius = p_Radius;
    return (initialConvectiveCoreRadius - p_Tau * (initialConvectiveCoreRadius - finalConvectiveCoreRadius));
}











///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 MASS CALCULATIONS                                 //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate convective core mass
 *
 * Assume equal to total mass at start (for continuity with stripped CHeB or HG star), continuity with HeHG at end of phase, linear growth
 *
 *
 * double CalculateConvectiveCoreMass()
 *
 * @return                                      Convective core mass (solar masses)
 */
double HeMS::CalculateConvectiveCoreMass() const {

    // We need core mass at end of phase, which is just the core mass at the start of the HeHG phase.
    // Since we are on the He main sequence here, we can clone this object as an HeHG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HeHG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.      
    std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP, OBJECT_PERSISTENCE::EPHEMERAL, true);
    HeHG* hgClone = dynamic_cast<HeHG*>(clone.get());
    double finalConvectiveCoreMass = hgClone->CoreMass(); // get core mass from clone
    // unique_ptr handles cleanup

    double initialConvectiveCoreMass = Mass();
    return (initialConvectiveCoreMass - Tau() * (initialConvectiveCoreMass - finalConvectiveCoreMass));
}


/*
 * Calculate rejuvenation factor for stellar age based on mass lost/gained during mass transfer
 *
 * Description?  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 *
 *
 * double CalculateMTRejuvenationFactor()
 *
 * @return                                      Rejuvenation factor
 */
double HeMS_Constituent::CalculateMTRejuvenationFactor() {

    const double mass     = m_Star->Mass();
    const double massPrev = m_Star->MassPrev();

    double fRej = 1.0;                                                                          // default value

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {

        case MT_REJUVENATION_PRESCRIPTION::HURLEY:                                              // use default Hurley et al. 2000 prescription = 1.0
            break;

        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                                           // StarTrack 2008 prescription - section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf
            fRej = utils::Compare(mass, massPrev) <= 0 ? 1.0 : massPrev / mass;                 // rejuvenation factor is unity for mass losing stars
            break;

        default:                                                                                // unknown prescription
            THROW_ERROR(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION);                           // throw error
    }

    return fRej;
}


/*
 * CalculateMLRate_Hurley2000
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for HeMS stars,
 * per Hurley et al. 2000.
 *
 * 
 * MassLossT CalculateMLRate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, double) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * << 1 unnamed double parameter to match BaseStar's virtual signature >>
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT HeMS::CalculateMLRate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, double) const {

    const double dMdtWR = std::get<0>(CalculateMLRateWR_Hurley2000(p_Luminosity, 0.0)) * OPTIONS->WolfRayetFactor();
    const double dMdtKR = std::get<0>(CalculateMLRate_KudritzkiReimers1978(p_Mass, p_Radius, p_Luminosity));
    const double dMdtNJ = std::get<0>(CalculateMLRate_NieuwenhuijzenDeJager1990(p_Mass, p_Radius, p_Luminosity));

    double dMdt = std::max(dMdtNJ, dMdtKR);         // Default mass loss rate
    ML_TYPE DominantMLtype = ML_TYPE::GB;           // Default dominant mass loss type

    if (dMdtWR > dMdt) {                            // WR ML rate bigger?
        dMdt = dMdtWR;                              // Yes - use WR rate
        dominantMLtype = ML_TYPE::WR;
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRate_Merritt2025
 *
 * @brief
 * Calculate mass loss rate, and dominant mass loss type, for HeMS stars,
 * per Merritt et al., 2025.
 *
 * 
 * MassLossT CalculateMLRate_Merritt2025(
 *     const double p_Metallicity,
 *     const double p_Mass,
 *     const double p_Radius,
 *     const double p_Luminosity,
 *     const double p_Temperature,
 *     const double p_HeAbundanceSurface,
 *     const double p_Perturb,
 *     const double p_CoolWindMultiplier,
 *     const double p_LBVFactor,
 *     double
 * ) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @param       p_Perturb                       Small envelope perturbation parameter
 * @param       p_CoolWindMultiplier            Cool wind ML multiplier (see --cool-wind-mass-loss-multiplier)
 * @param       p_LBVFactor                     LBV factor (see --luminous-blue-variable-multiplier)
 * @param       p_mStart                        Not used here - unnamed to suppress compiler warning
 * @return                                      Tuple containing:
 *                                                   double  WR mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT HeMS::CalculateMLRate_Merritt2025(
    const double p_Metallicity,
    const double p_Mass,
    const double p_Radius,
    const double p_Luminosity,
    const double p_Temperature,
    const double p_HeAbundanceSurface,
    const double p_Perturb,
    const double p_CoolWindMultiplier,
    const double p_LBVFactor,
    double
) const {

    ML_TYPE DominantMLtype;
    double dMdt;

    switch (OPTIONS->WRMassLossPrescription()) {                                            // Which WR mass loss prescription?

        case WR_ML_PRESCRIPTION::BELCZYNSKI2010:                                            // BELCZYNSKI2010
            std::tie(dMdt, dominantMLtype) = CalculateMLRate_Belczynski2010(
                p_Metallicity,
                p_Mass,
                p_Radius,
                p_Luminosity,
                p_Temperature,
                p_HeAbundanceSurface,
                p_Perturb,
                p_CoolWindMultiplier,
                p_LBVFactor,
                OPTIONS->LBVMassLossPrescription()
            );
            break;

        case WR_ML_PRESCRIPTION::SANDERVINK2023: {                                          // SANDERVINK2023

            // Start with Sander & Vink 2020
            std::tie(dMdt, dominantMLtype) = CalculateMLRateWR_SanderVink2020(p_Metalliciy, p_Luminosity, 0.0);

            // Apply the Sander et al. 2023 temperature correction to the Sander & Vink 2020
            // rate if necessary - gives the Sander & Vink 2023 rate
            // https://arxiv.org/abs/2301.01785
            // 
            // Use the correction given in eq 18, with the effective temperature
            // (what they refer to as T_\star in eq 1) as T_eff,crit
            if (dMdt > 0.0) {                                                               // Only apply the correction for positive mass loss rates
                constexpr double teffMin = 100.0E3;                                         // Minimum effective temperature (K) for which correction applies
                constexpr double teffRef = 141.0E3;                                         // Reference effective temperature in Kelvin
                const double     teff    = p_Temperature * TSOL;                            // Effective temperature in Kelvin
                if (teff > teffMin) {                                                       // Correction applicable?
                    dMdt = PPOW(10.0, std::log10(dMdt) - 6.0 * std::log10(teff / teffRef)); // Yes, apply correction - gives Sander & Vink 2023 mass loss rate
                }
            }

            // Compare Sander & Vink 2020/2023 mass loss rate to Vink 2017, and clamp to a minimum of
            // Vink 2017 - will typically result in Vink 2017 mass loss rate for low mass or luminosity,
            // and Sander & Vink 2020/2023 for high mass or luminosity
            const auto [dMdtVink2017, dominantMLtypeVink2017] = CalculateMLRate_Vink2017(p_Metallicity, p_Luminosity);
            if (dMdtVink2017 > dMdt) {
                dMdt           = dMdtVink2017;
                dominantMLtype = dominantMLtypeVink2017;
            }

            } break;

        case WR_ML_PRESCRIPTION::SHENAR2019: {                                              // SHENAR2019

            // Start with Shenar+ 2019
            std::tie(dMdt, dominantMLtype) = CalculateMLRateWR_Shenar2019(p_Metallicity, p_Luminosity, p_Temperature);

            // Compare Shenar+ 2019 mass loss rate to Vink 2017, and clamp to a minimum
            // of Vink 2017 to avoid extrapolating to low luminosity
            const auto [dMdtVink2017, dominantMLtypeVink2017] = CalculateMLRate_Vink2017(p_Metallicity, p_Luminosity);
            if (dMdtVink2017 > dMdt) {
                dMdt           = dMdtVink2017;
                dominantMLtype = dominantMLtypeVink2017;
            }

            } break;

        case WR_ML_PRESCRIPTION::ZERO:                                                      // ZERO
            dMdt           = 0.0;                                                           // No mass loss
            dominantMLtype = ML_TYPE::NONE;
            break;

        default:                                                                            // Unexpected prescription
            // The only way this can happen is if the WR_ML_PRESCRIPTION served to this function is
            // not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR_STATIC(ERROR::UNEXPECTED_WR_ML_PRESCRIPTION);                       // Throw error
    }

    // Return mass loss rate with user supplied WR factor applied
    return std::make_tuple(dMdt * OPTIONS->WolfRayetFactor(), dominantMLtype);
}





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              LIFETIME / AGE FUNCIONS                              //
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
STELLAR_TYPE HeMS::ResolveEnvelopeLoss(bool p_Force) {

    STELLAR_TYPE stellarType = m_StellarType;

    if (p_Force || utils::Compare(m_Mass, 0.0) <= 0) {
        stellarType = STELLAR_TYPE::MASSLESS_REMNANT;
        m_Radius    = 0.0;
        m_Mass      = 0.0;
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
STELLAR_TYPE HeMS::EvolveToNextPhase() {
#define timescales(x) m_InterimState.HurleyTimescales(TS::x)

    m_InterimState.SetAge(timescales(HeMS));

    return STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;
#undef timescales
}
























//////////////////////////////// HeMS_Constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

/* 
 * CalculateCriticalMassRatio_Ge2020_Interpolate
 *
 * @brief
 * Calculate critical mass ratios for H-poor stars, per Ge et al. DATE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 * 
 * Interpolate in log mass and log radius to calculate the stellar response of a He star to mass loss.
 *
 * The existing data table only applies for fully conservative mass transfer and the Ge et al. <<<<<DATE>>>>> <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 * fully adiabatic response, not the artificially isentropic one. Also only for Z = Zsol.   <<<< WHICH ZSOL?????? <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 *
 * 
 * double CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_MTefficiency                  Mass transfer accretion efficiency (ignored by this function for now)
 * @return                                      Critical mass ratio (aka qCrit)
 */ 
inline double MainSequence_Constituent::CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency) const {

    const double logMass   = std::log10(p_Mass);   
    const double logRadius = std::log10(p_Radius);

    const DBL_VECTOR massVec = std::get<0>(QCRIT_GE_HE_STAR);   // vector of masses from QCRIT_GE_HE_STAR

    SizeTVectorT indices = utils::BinarySearch(massVec, p_Mass);
    SizeT lowMIdx = indices[0];
    SizeT uppMIdx = indices[1];
    
    // if masses are out of range, set to endpoints
    if (lowMIdx == SIZE_MAX) {                                  // below minimum?
        lowMIdx = 0;                                            // yes - lower edge of lowest bin 
        uppMIdx = 1;                                            // upper edge of lowest bin
    } 
    else if (uppMIdx == SIZE_MAX) {                             // no - above maximum?
        lowMIdx = massVec.size() - 2;                           // yes - lower edge of highest bin 
        uppMIdx = massVec.size() - 1;                           // upper edge of highest bin
    } 
    
    // vector of radii and qCrits from QCRIT_GE_HE_STAR
    const std::vector<std::tuple<DBL_VECTOR, DBL_VECTOR>> radiiQCritsVec = std::get<1>(QCRIT_GE_HE_STAR);

    // (log) radii vectors for mass bin bounds
    const DBL_VECTOR radiiVecLow = std::get<0>(radiiQCritsVec[lowMIdx]);
    const DBL_VECTOR radiiVecUpp = std::get<0>(radiiQCritsVec[uppMIdx]);

    // qCrit vectors for mass bin bounds
    const DBL_VECTOR qCritVecLow = std::get<1>(radiiQCritsVec[lowMIdx]);
    const DBL_VECTOR qCritVecUpp = std::get<1>(radiiQCritsVec[uppMIdx]);

    indices = utils::BinarySearch(radiiVecLow, logRadius);      // find lower mass bound radii bin for p_Radius
    int lowRLowMIdx = indices[0];                               // bin lower bound index
    int uppRLowMIdx = indices[1];                               // bin upper bound index
    
    // if radii are out of range, set to endpoints
    if (lowRLowMIdx == -1) {                                    // below minimum?                                      
        lowRLowMIdx = 0;                                        // yes - lower edge of lowest bin  
        uppRLowMIdx = 1;                                        // upper edge of lowest bin 
    }
    else if (uppRLowMIdx == -1) {                               // no - above maximum?                                             
        lowRLowMIdx = radiiVecLow.size() - 2;                   // yes - lower edge of highest bin 
        uppRLowMIdx = radiiVecLow.size() - 1;                   // upper edge of highest bin 
    }
    
    indices = utils::BinarySearch(radiiVecUpp, logRadius);      // find upper mass bound radii bin for p_Radius
    int lowRUppMIdx = indices[0];
    int uppRUppMIdx = indices[1];
    
    // if radii are out of range, set to endpoints
    if (lowRUppMIdx == -1) {                                    // below minimum? 
        lowRUppMIdx = 0;                                        // yes - lower edge of lowest bin  
        uppRUppMIdx = 1;                                        // upper edge of lowest bin 
    }
    else if (uppRUppMIdx == -1) {                               // no - above maximum? 
        lowRUppMIdx = radiiVecUpp.size() - 2;                   // yes - lower edge of highest bin
        uppRUppMIdx = radiiVecUpp.size() - 1;                   // upper edge of highest bin
    }
    
    // set the boundary points for the 2D interpolation
    const double rLowMLowR = radiiVecLow[lowRLowMIdx];
    const double rLowMUppR = radiiVecLow[uppRLowMIdx];
    const double rUppMlowR = radiiVecUpp[lowRUppMIdx];
    const double rUppMUppR = radiiVecUpp[uppRUppMIdx];
    
    const double qLowMLowR = qCritVecLow[lowRLowMIdx];
    const double qLowMUppR = qCritVecLow[uppRLowMIdx];
    const double qUppMLowR = qCritVecUpp[lowRUppMIdx];
    const double qUppMUppR = qCritVecUpp[uppRUppMIdx];
    
    // interpolate on logR first, then logM, using nearest neighbour for extrapolation
    // qCrit lower mass
    double qCritLowM;
         if (logRadius < rLowMLowR) qCritLowM = qLowMLowR;      // below lower bound?  Use lower bound
    else if (logRadius > rLowMUppR) qCritLowM = qLowMUppR;      // no - above upper bound?  Use upper bound
    else                            qCritLowM = qLowMLowR + (rLowMUppR - logRadius) / (rLowMUppR - rLowMLowR) * (qLowMUppR - qLowMLowR); // no - interpolate

    // qCrit upper mass
    double qCritUppM;
         if (logRadius < rUppMlowR) qCritUppM = qUppMLowR;      // below lower bound?  Use lower bound
    else if (logRadius > rUppMUppR) qCritUppM = qUppMUppR;      // no - above upper bound?  Use upper bound
    else                            qCritUppM = qUppMLowR + (rUppMUppR - logRadius) / (rUppMUppR - rUppMlowR) * (qUppMUppR - qUppMLowR); // no - interpolate

    // qCrit
    double qCrit;

    const double logLowM = std::log10(massVec[lowMIdx]);
    const double logUppM = std::log10(massVec[uppMIdx]);

         if (logMass < logLowM) qCrit = qCritLowM;        // below lower bound?
    else if (logMass > logUppM) qCrit = qCritUppM;        // no - above upper bound?  Use upper bound
    else                             qCrit = qCritLowM + (logUppM - logMass) / (logUppM - logLowM) * (qCritUppM - qCritLowM); // no - interpolate

    return qCrit;
}



/*
 * HeMS_Constituent::CalculateZetaEquilibrium()
 *
 */
double HeMS_Constituent::CalculateZetaEquilibrium() {
    const double mass   = m_Star->Mass();
    const double radius = m_Star->Radius();
    const double deltaMass           = -mass / 1.0E5;
    const double radiusAfterMassGain = m_Star->CalculateRadiusOnMassChange(deltaMass);
    return (radiusAfterMassGain - radius) / deltaMass * mass / radius;
}
