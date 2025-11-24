#include "HeMS.h"
#include "HeWD.h"


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//             PARAMETERS, MISCELLANEOUS CALCULATIONS AND FUNCTIONS ETC.             //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////







/*
 * CalculateGBparams_Hurley2000
 *
 * @brief
 * Calculate Giant Branch (GB) parameters per Hurley et al. 2000
 *
 * Since Giant Branch parameters depend on a star's mass, they need to be calculated
 * whenever the mass of the star changes (probably every timestep).
 *
 *
 * DBL_VECTOR CalculateGBparams_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @return                                      Mutated GB parameters (Myr)
 */
COMPAS_PURE DBL_VECTOR HeMS::CalculateGBparams_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams) const {

    DBL_VECTOR GBparams = p_GBparams;   // copy given GBparams

	GBparams[static_cast<int>(HURLEY_GBP:::B)] = CalculateCoreMass_Luminosity_B_Hurley2000(p_Mass);
	GBparams[static_cast<int>(HURLEY_GBP:::D)] = CalculateCoreMass_Luminosity_D_Hurley2000(p_Mass);

    GBparams[static_cast<int>(HURLEY_GBP:::p)] = CalculateCoreMass_Luminosity_p_Hurley2000(p_Mass);
    GBparams[static_cast<int>(HURLEY_GBP:::q)] = CalculateCoreMass_Luminosity_q_Hurley2000(p_Mass);
    
    GBparams[static_cast<int>(HURLEY_GBP:::Mx)] = GiantBranch::CalculateCoreMass_Luminosity_Mx_Hurley2000(GBparams);

    // return GB parameters vector by value - NRVO takes care of performance/efficiency
    return GBparams;
}








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
      
    HeHG *clone = HeHG::Clone(static_cast<HeHG&>(const_cast<HeMS&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL);
    double finalConvectiveCoreRadius = clone->CalculateConvectiveCoreRadius();                  // get core radius from clone
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone

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
      
    HeHG *clone = HeHG::Clone(static_cast<HeHG&>(const_cast<HeMS&>(*this)), OBJECT_PERSISTENCE::EPHEMERAL, true);
    double finalConvectiveCoreMass = clone->CoreMass();                                         // get core mass from clone
    delete clone; clone = nullptr;                                                              // return the memory allocated for the clone

    double initialConvectiveCoreMass = m_Mass;
    return (initialConvectiveCoreMass - m_Tau * (initialConvectiveCoreMass - finalConvectiveCoreMass));
}


/*
 * Calculate rejuvenation factor for stellar age based on mass lost/gained during mass transfer
 *
 * Description?
 *
 *
 * double CalculateMTRejuvenationFactor()
 *
 * @return                                      Rejuvenation factor
 */
double HeMS::CalculateMTRejuvenationFactor() {

    double fRej = 1.0;                                                                          // default value

    switch (OPTIONS->MassTransferRejuvenationPrescription()) {

        case MT_REJUVENATION_PRESCRIPTION::HURLEY:                                              // use default Hurley et al. 2000 prescription = 1.0
            break;

        case MT_REJUVENATION_PRESCRIPTION::STARTRACK:                                           // StarTrack 2008 prescription - section 5.6 of http://arxiv.org/pdf/astro-ph/0511811v3.pdf
            fRej = utils::Compare(m_Mass, m_MassPrev) <= 0 ? 1.0 : m_MassPrev / m_Mass;         // rejuvenation factor is unity for mass losing stars
            break;

        default:                                                                                // unknown prescription
            // the only way this can happen is if someone added a MT_REJUVENATION_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MT_REJUVENATION_PRESCRIPTION);                           // throw error
    }

    return fRej;
}


/*
 * CalculateMLrate_Hurley2000
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for HeMS stars,
 * per Hurley et al. 2000.
 *
 * 
 * MASS_LOSS_T CalculateMLrate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_PerturbationMu) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for massive stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T HeMS::CalculateMLrate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_PerturbationMu) const {

    const double dMdtWR = CalculateMLrateWR_Hurley2000(p_Luminosity, 0.0) * OPTIONS->WolfRayetFactor();
    const double dMdtKR = CalculateMLrate_KudritzkiReimers1978(p_Mass, p_Radius, p_Luminosity);
    const double dMdtNJ = CalculateMLrate_NieuwenhuijzenDeJager1990(p_Mass, p_Radius, p_Luminosity);

    double dMdt = std::max(dMdtNJ, dMdtKR);                 // default mass loss rate
    MASS_LOSS_TYPE dominantMLtype = MASS_LOSS_TYPE::GB;     // default dominant mass loss type

    if (dMdtWR > dMdt) {                                    // WR ML rate bigger?
        dMdt = dMdtWR;                                      // yes - use WR rate
        dominantMLtype = MASS_LOSS_TYPE::WR;
    }

    // return mass loss rate with user supplied WR factor applied
    return std::make_tuple(dMdt, dominantMLtype);
}








/*
 * CalculateMLrate_Merritt2025_Static
 *
 * @brief
 * Calculate mass loss rate, and dominant mass loss type, for HeMS stars,
 * per Merritt et al., 2025.
 *
 * 
 * static MASS_LOSS_T CalculateMLrate_Merritt2025_Static(const double p_Luminosity, const double p_Temperature)
 * 
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE static MASS_LOSS_T HeMS::CalculateMLrate_Merritt2025_Static(const double p_Luminosity, const double p_Temperature) {


    MASS_LOSS_TYPE dominantMLtype;
    double dMdt;

    switch (OPTIONS->WRMassLossPrescription()) {                                                     // which WR mass loss prescription?

        case WR_MASS_LOSS_PRESCRIPTION::BELCZYNSKI2010:                             // BELCZYNSKI2010
            std::tie(dMdt, dominantMLtype) = HeMS::CalculateMLrate_Belczynski2010_Static(GLOBALS->Metallicity(), p_Luminosity);
            break;

        case WR_MASS_LOSS_PRESCRIPTION::SANDERVINK2023:                             // SANDERVINK2023

            // start with Sander & Vink 2020
            std::tie(dMdt, dominantMLtype) = BaseStar::CalculateMLrateWR_SanderVink2020_Static(p_Luminosity, 0.0, GLOBALS->ZetaAnders());

            // apply the Sander et al. 2023 temperature correction to the Sander & Vink 2020
            // rate if necessary - gives the Sander & Vink 2023 rate
            // https://arxiv.org/abs/2301.01785
            // 
            // use the correction given in eq 18, with the effective temperature
            // (what they refer to as T_\star in eq 1) as T_eff,crit
             
            if (dMdt > 0.0) {                                                       // only apply the correction for positive mass loss rates
                constexpr double teffMin = 100.0E3;                                 // minimum effective temperature (K) for which correction applies
                constexpr double teffRef = 141.0E3;                                 // reference effective temperature in Kelvin
                const double     teff    = p_Temperature * TSOL;                    // effective temperature in Kelvin
                if (teff > teffMin) {                                               // correction applicable?
                    dMdt = PPOW(10.0, log10(dMdt) - 6.0 * log10(teff / teffRef));   // yes, apply correction - gives Sander & Vink 2023 mass loss rate
                }
            }

            // compare Sander & Vink 2020/2023 mass loss rate to Vink 2017, and clamp to a minimum of
            // Vink 2017 - will typically result in Vink 2017 mass loss rate for low mass or luminosity,
            // and Sander & Vink 2020/2023 for high mass or luminosity

            const double dMdtVink2017;
            const double dominantMLtypeVink2017;
            std::tie(dMdtVink2017, dominantMLtypeVink2017) = HeMS::CalculateMLrate_Vink2017_Static(p_Luminosity, GLOBALS->ZetaAnders());
            if (dMdtVink2017 > dMdt) {
                dMdt           = dMdtVink2017;
                dominantMLtype = dominantMLtypeVink2017;
            }
            break;

        case WR_MASS_LOSS_PRESCRIPTION::SHENAR2019:                                 // SHENAR2019

            // start with Shenar+ 2019
            std::tie(dMdt, dominantMLtype) = HeMS::CalculateMLrateWR_Shenar2019_Static(p_Luminosity, p_Temperature, GLOBALS->SigmaHurley());

            // compare Shenar+ 2019 mass loss rate to Vink 2017, and clamp to a minimum
            // of Vink 2017 to avoid extrapolating to low luminosity

            const double dMdtVink2017;
            const double dominantMLtypeVink2017;
            std::tie(dMdtVink2017, dominantMLtypeVink2017) = HeMS::CalculateMLrate_Vink2017_Static(p_Luminosity, GLOBALS->ZetaAnders());
            if (dMdtVink2017 > dMdt) {
                dMdt           = dMdtVink2017;
                dominantMLtype = dominantMLtypeVink2017;
            }
            break;

        case WR_MASS_LOSS_PRESCRIPTION::ZERO:                                       // ZERO
            dMdt           = 0.0;                                                   // no mass loss
            dominantMLtype = MASS_LOSS_TYPE::NONE;
            break;

        default:                                                                    // unexpected prescription
            // the only way this can happen is if the WR_MASS_LOSS_PRESCRIPTION passed to this function
            // is not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.
            THROW_ERROR(ERROR::UNEXPECTED_WR_MASS_LOSS_PRESCRIPTION);               // throw error
    }

    // return mass loss rate with user supplied WR factor applied
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
double HeMS::ChooseTimestep(const double p_Time) const {

    double dtk = 0.05 * timescales(tHeMS);
    double dte = timescales(tHeMS) - p_Time;

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);
}


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

    m_Age = timescales(tHeMS);

    return STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP;
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
double MainSequence_Constituent::CalculateCriticalMassRatio_Ge2020_Interpolate(const double p_Mass, const double p_Radius, const double p_MTefficiency) const {

    const double logMass   = log10(p_Mass);   
    const double logRadius = log10(p_Radius);

    const DBL_VECTOR massVec = std::get<0>(QCRIT_GE_HE_STAR);   // vector of masses from QCRIT_GE_HE_STAR

    INT_VECTOR indices = utils::BinarySearch(massVec, p_Mass);  // find mass bin for p_Mass
    int lowMIdx = indices[0];                                   // bin lower bound index
    int uppMIdx = indices[1];                                   // bin upper bound index
    
    // if masses are out of range, set to endpoints
    if (lowMIdx == -1) {                                        // below minimum?
        lowMIdx = 0;                                            // yes - lower edge of lowest bin 
        uppMIdx = 1;                                            // upper edge of lowest bin
    } 
    else if (uppMIdx == -1) {                                   // no - above maximum? 
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
         if (logRadius < rUppMLowR) qCritUppM = qUppMLowR;      // below lower bound?  Use lower bound
    else if (logRadius > rUppMUppR) qCritUppM = qUppMUppR;      // no - above upper bound?  Use upper bound
    else                            qCritUppM = qUppMLowR + (rUppMUppR - logRadius) / (rUppMUppR - rUppMlowR) * (qUppMUppR - qUppMLowR); // no - interpolate

    // qCrit
    double qCrit;

    const double logLowM = log10(massVec[lowMIdx]);
    const double logUppM = log10(massVec[uppMIdx]);

         if (logMass < logLowerMass) qCrit = qCritLowM;         // below lower bound?
    else if (logMass > logUpperMass) qCrit = qCritUppM;         // no - above upper bound?  Use upper bound
    else                             qCrit = qCritLowM + (logUppM - logMass) / (logUppM - logLowM) * (qCritUppM - qCritLowM); // no - interpolate

    return qCrit;
}

