#include "BaseStar.h"







    










///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  CLASS FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * StellarPropertyValue
 *
 * Determine the value of the requested property of the constituent star (parameter p_Property)
 *
 * The property is a boost variant variable, and is one of the following types:
 *
 *      STAR_PROPERTY           - any individual star property
 *      STAR_1_PROPERTY         - property of the primary (m_Star1)
 *      STAR_2_PROPERTY         - property of the secondary (m_Star2)
 *      SUPERNOVA_PROPERTY      - property of the star that has gone supernova
 *      COMPANION_PROPERTY      - property of the companion to the supernova
 *      BINARY_PROPERTY         - property of the binary
 *      PROGRAM_OPTION          - program option
 *
 * This function handles properties of type:
 * 
 *    STAR_PROPERTY, STAR_1_PROPERTY, STAR_2_PROPERTY, SUPERNOVA_PROPERTY, COMPANION_PROPERTY
 * 
 * only - anything else will result in an error being thrown and the evolution of the star (or binary)
 * terminated.
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.  
 *
 *
 * COMPAS_VARIABLE StellarPropertyValue(const T_ANY_PROPERTY p_Property)
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BaseStar::StellarPropertyValue(const T_ANY_PROPERTY p_Property) const {

    COMPAS_VARIABLE value;

    ANY_STAR_PROPERTY property;

    switch (boost::apply_visitor(VariantPropertyType(), p_Property)) {

        case ANY_PROPERTY_TYPE::T_STAR_PROPERTY     : { STAR_PROPERTY      prop = boost::get<STAR_PROPERTY>(p_Property);      property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_STAR_1_PROPERTY   : { STAR_1_PROPERTY    prop = boost::get<STAR_1_PROPERTY>(p_Property);    property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_STAR_2_PROPERTY   : { STAR_2_PROPERTY    prop = boost::get<STAR_2_PROPERTY>(p_Property);    property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_SUPERNOVA_PROPERTY: { SUPERNOVA_PROPERTY prop = boost::get<SUPERNOVA_PROPERTY>(p_Property); property = (ANY_STAR_PROPERTY)prop; } break;
        case ANY_PROPERTY_TYPE::T_COMPANION_PROPERTY: { COMPANION_PROPERTY prop = boost::get<COMPANION_PROPERTY>(p_Property); property = (ANY_STAR_PROPERTY)prop; } break;

        default:                                                                                                        // unexpected stellar property type
            // the only ways this can happen are if someone added a stellar type property (into ANY_PROPERTY_TYPE)
            // and it isn't accounted for in this code, or if there is a defect in the code that causes
            // this function to be called with a bad parameter.  We should not default here, with or without a
            // warning - this is a code defect, so we flag it as an error and that will result in termination of
            // the evolution of the star or binary.
            // The correct fix for this is to add code for the missing property type or find and fix the code defect.

            THROW_ERROR(ERROR::UNEXPECTED_STELLAR_PROPERTY_TYPE);                                                       // throw error
    }

    switch (property) {
        case ANY_STAR_PROPERTY::AGE:                                                value = Age();                                                  break;
        case ANY_STAR_PROPERTY::ANGULAR_MOMENTUM:                                   value = AngularMomentum();                                      break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_FIXED:                               value = CalculateBindingEnergy(OPTIONS->CommonEnvelopeLambda()); break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_NANJING:                             value = CalculateBindingEnergy(CalculateCELambda_Nanjing());    break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE:                           value = CalculateBindingEnergy(CalculateLambdaLoveridge());     break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS:                     value = CalculateBindingEnergy(CalculateLambdaLoveridge(m_Mass - m_CoreMass, true)); break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_KRUCKOW:                             value = CalculateBindingEnergy(CalculateCELambdaKruckow());     break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_CONVECTIVE_ENVELOPE:                 value = CalculateConvectiveEnvelopeBindingEnergy(CalculateConvectiveEnvelopeLambdaPicker(CalculateConvectiveEnvelopeMass(Mass(), CoreMass()))); break;
        case ANY_STAR_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE:               value = CHonMS();                                               break;
        case ANY_STAR_PROPERTY::CO_CORE_MASS:                                       value = COCoreMass();                                           break;
        case ANY_STAR_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION:           value = SN_COCoreMassAtCOFormation();                           break;
        case ANY_STAR_PROPERTY::CONVECTIVE_ENV_MASS:                                std::tie(value, std::ignore) = CalculateConvectiveEnvelopeMass(Mass(), CoreMass());  break;
        case ANY_STAR_PROPERTY::CORE_MASS:                                          value = CoreMass();                                             break;
        case ANY_STAR_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION:              value = SN_CoreMassAtCOFormation();                             break;
        case ANY_STAR_PROPERTY::CORE_RADIUS_AT_COMPACT_OBJECT_FORMATION:            value = SN_CoreRadiusAtCOFormation();                           break; 
        case ANY_STAR_PROPERTY::DRAWN_KICK_MAGNITUDE:                               value = SN_DrawnKickMagnitude();                                break;
        case ANY_STAR_PROPERTY::DOMINANT_MASS_LOSS_TYPE:                            value = DominantMassLossType();                                 break; /* was _RATE <<<<<<*/
        case ANY_STAR_PROPERTY::DT:                                                 value = Dt();                                                   break;
        case ANY_STAR_PROPERTY::DYNAMICAL_TIMESCALE:                                value = CalculateTimescale_Dynamical();                          break;
        case ANY_STAR_PROPERTY::ECCENTRIC_ANOMALY:                                  value = SN_EccentricAnomaly();                                  break;
        case ANY_STAR_PROPERTY::ENV_MASS:                                           value = Mass() - CoreMass();                                    break;
        case ANY_STAR_PROPERTY::ERROR:                                              value = Error();                                                break;
        case ANY_STAR_PROPERTY::EVOL_STATUS:                                        value = EvolutionStatus();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_AIC:                                    value = ExperiencedAIC();                                       break;
        case ANY_STAR_PROPERTY::EXPERIENCED_HeSD:                                   value = ExperiencedHeSD();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_CCSN:                                   value = ExperiencedCCSN();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_ECSN:                                   value = ExperiencedECSN();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_PISN:                                   value = ExperiencedPISN();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_PPISN:                                  value = ExperiencedPPISN();                                     break;
        case ANY_STAR_PROPERTY::EXPERIENCED_SNIA:                                   value = ExperiencedSNIA();                                      break;
        case ANY_STAR_PROPERTY::EXPERIENCED_SN_TYPE:                                value = ExperiencedSN_Type();                                   break;
        case ANY_STAR_PROPERTY::EXPERIENCED_USSN:                                   value = ExperiencedUSSN();                                      break;
        case ANY_STAR_PROPERTY::FALLBACK_FRACTION:                                  value = SN_FallbackFraction();                                  break;
        case ANY_STAR_PROPERTY::MASS_TRANSFER_DONOR_HISTORY:                        value = GetMassTransferDonorHistoryString();                    break;
        case ANY_STAR_PROPERTY::HE_CORE_MASS:                                       value = HeCoreMass();                                           break;
        case ANY_STAR_PROPERTY::HE_CORE_MASS_AT_COMPACT_OBJECT_FORMATION:           value = SN_HeCoreMassAtCOFormation();                           break;
        case ANY_STAR_PROPERTY::HELIUM_ABUNDANCE_CORE:                              value = HeliumAbundanceCore();                                  break;
        case ANY_STAR_PROPERTY::HELIUM_ABUNDANCE_SURFACE:                           value = HeliumAbundanceSurface();                               break;
        case ANY_STAR_PROPERTY::HYDROGEN_ABUNDANCE_CORE:                            value = HydrogenAbundanceCore();                                break;
        case ANY_STAR_PROPERTY::HYDROGEN_ABUNDANCE_SURFACE:                         value = HydrogenAbundanceSurface();                             break;
        case ANY_STAR_PROPERTY::IS_HYDROGEN_POOR:                                   value = SN_IsHydrogenPoor();                                    break;
        case ANY_STAR_PROPERTY::ID:                                                 value = ObjectId();                                             break;
        case ANY_STAR_PROPERTY::INITIAL_HELIUM_ABUNDANCE:                           value = CalculateInitialHeliumAbundance();                      break;
        case ANY_STAR_PROPERTY::INITIAL_HYDROGEN_ABUNDANCE:                         value = CalculateInitialHydrogenAbundance();                    break;
        case ANY_STAR_PROPERTY::INITIAL_STELLAR_TYPE:                               value = InitialStellarType();                                   break;
        case ANY_STAR_PROPERTY::INITIAL_STELLAR_TYPE_NAME:                          value = STELLAR_TYPE_LABEL.at(InitialStellarType());            break;
        case ANY_STAR_PROPERTY::IS_AIC:                                             value = IsAIC();                                                break;
        case ANY_STAR_PROPERTY::IS_CCSN:                                            value = IsCCSN();                                               break;
        case ANY_STAR_PROPERTY::IS_HeSD:                                            value = IsHeSD();                                               break;
        case ANY_STAR_PROPERTY::IS_ECSN:                                            value = IsECSN();                                               break;
        case ANY_STAR_PROPERTY::IS_PISN:                                            value = IsPISN();                                               break;
        case ANY_STAR_PROPERTY::IS_PPISN:                                           value = IsPPISN();                                              break;
        case ANY_STAR_PROPERTY::IS_SNIA:                                            value = IsSNIA();                                               break;
        case ANY_STAR_PROPERTY::IS_USSN:                                            value = IsUSSN();                                               break;
        case ANY_STAR_PROPERTY::KICK_MAGNITUDE:                                     value = SN_KickMagnitude();                                     break;
        case ANY_STAR_PROPERTY::LAMBDA_CONVECTIVE_ENVELOPE:                         value = CalculateConvectiveEnvelopeLambdaPicker(CalculateConvectiveEnvelopeMass(Mass(), CoreMass())); break;
        case ANY_STAR_PROPERTY::LAMBDA_DEWI:                                        value = CalculateCELambda_Dewi();                                  break;
        case ANY_STAR_PROPERTY::LAMBDA_FIXED:                                       value = OPTIONS->CommonEnvelopeLambda();                        break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW:                                     value = CalculateCELambdaKruckow();                             break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_BOTTOM:                              value = CalculateCELambdaKruckow(m_Radius, -1.0);               break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_MIDDLE:                              value = CalculateCELambdaKruckow(m_Radius, -4.0 / 5.0);         break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_TOP:                                 value = CalculateCELambdaKruckow(m_Radius, -2.0 / 3.0);         break;
        case ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE:                                   value = CalculateLambdaLoveridge(m_Mass - m_CoreMass, false);   break;
        case ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE_WINDS:                             value = CalculateLambdaLoveridge(m_Mass - m_CoreMass, true);    break;
        case ANY_STAR_PROPERTY::LAMBDA_NANJING:                                     value = CalculateCELambda_Nanjing();                            break;
        case ANY_STAR_PROPERTY::LBV_PHASE_FLAG:                                     value = LBV_PhaseFlag();                                        break;
        case ANY_STAR_PROPERTY::LUMINOSITY:                                         value = Luminosity();                                           break;
        case ANY_STAR_PROPERTY::MASS:                                               value = Mass();                                                 break;
        case ANY_STAR_PROPERTY::MASS_0:                                             value = Mass0();                                                break;
        case ANY_STAR_PROPERTY::MDOT:                                               value = Mdot();                                                 break;
        case ANY_STAR_PROPERTY::MEAN_ANOMALY:                                       value = SN_MeanAnomaly();                                       break;
        case ANY_STAR_PROPERTY::METALLICITY:                                        value = GLOBALS->Metallicity();                                          break;
        case ANY_STAR_PROPERTY::MOMENT_OF_INERTIA:                                  value = CalculateMomentOfInertia();                             break;
        case ANY_STAR_PROPERTY::MZAMS:                                              value = MZAMS();                                                break;
        case ANY_STAR_PROPERTY::OMEGA:                                              value = Omega() / SECONDS_IN_YEAR;                              break;
        case ANY_STAR_PROPERTY::OMEGA_BREAK:                                        value = OmegaBreak() / SECONDS_IN_YEAR;                         break;
        case ANY_STAR_PROPERTY::OMEGA_ZAMS:                                         value = OmegaZAMS() / SECONDS_IN_YEAR;                          break;
        case ANY_STAR_PROPERTY::PULSAR_MAGNETIC_FIELD:                              value = PulsarMagneticField();                                  break;
        case ANY_STAR_PROPERTY::PULSAR_SPIN_DOWN_RATE:                              value = PulsarSpinDownRate();                                   break;
        case ANY_STAR_PROPERTY::PULSAR_SPIN_FREQUENCY:                              value = PulsarSpinFrequency();                                  break;
        case ANY_STAR_PROPERTY::PULSAR_SPIN_PERIOD:                                 value = PulsarSpinPeriod();                                     break;
        case ANY_STAR_PROPERTY::PULSAR_BIRTH_PERIOD:                                value = PulsarBirthPeriod();                                    break;
        case ANY_STAR_PROPERTY::PULSAR_BIRTH_SPIN_DOWN_RATE:                        value = PulsarBirthSpinDownRate();                              break;
        case ANY_STAR_PROPERTY::RADIAL_EXPANSION_TIMESCALE:                         value = CalculateRadialExpansionTimescale();                    break;
        case ANY_STAR_PROPERTY::RADIUS:                                             value = Radius();                                               break;
        case ANY_STAR_PROPERTY::RANDOM_SEED:                                        value = RandomSeed();                                           break;
        case ANY_STAR_PROPERTY::ROCKET_KICK_MAGNITUDE:                              value = SN_RocketKickMagnitude();                               break;
        case ANY_STAR_PROPERTY::ROCKET_KICK_PHI:                                    value = SN_RocketKickPhi();                                     break;
        case ANY_STAR_PROPERTY::ROCKET_KICK_THETA:                                  value = SN_RocketKickTheta();                                   break;
        case ANY_STAR_PROPERTY::RZAMS:                                              value = RZAMS();                                                break;
        case ANY_STAR_PROPERTY::SN_TYPE:                                            value = SN_Type();                                              break;
        case ANY_STAR_PROPERTY::SPEED:                                              value = Speed();												break;
        case ANY_STAR_PROPERTY::STELLAR_TYPE:                                       value = StellarType();                                          break;
        case ANY_STAR_PROPERTY::STELLAR_TYPE_NAME:                                  value = STELLAR_TYPE_LABEL.at(StellarType());                   break;
        case ANY_STAR_PROPERTY::STELLAR_TYPE_PREV:                                  value = StellarTypePrev();                                      break;
        case ANY_STAR_PROPERTY::STELLAR_TYPE_PREV_NAME:                             value = STELLAR_TYPE_LABEL.at(StellarTypePrev());               break;
        case ANY_STAR_PROPERTY::SUPERNOVA_KICK_MAGNITUDE_RANDOM_NUMBER:             value = SN_KickMagnitudeRandom();                               break;
        case ANY_STAR_PROPERTY::SUPERNOVA_PHI:                                      value = SN_Phi();                                               break;
        case ANY_STAR_PROPERTY::SUPERNOVA_THETA:                                    value = SN_Theta();                                             break;
        case ANY_STAR_PROPERTY::TEMPERATURE:                                        value = Temperature() * TSOL;                                   break;
        case ANY_STAR_PROPERTY::THERMAL_TIMESCALE:                                  value = CalculateTimescale_Thermal();                            break;
        case ANY_STAR_PROPERTY::TIME:                                               value = Time();                                                 break;
        case ANY_STAR_PROPERTY::TIMESCALE_MS:                                       value = Timescale(TIMESCALE::tMS);                              break;
        case ANY_STAR_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION:             value = SN_TotalMassAtCOFormation();                            break;
        case ANY_STAR_PROPERTY::TOTAL_RADIUS_AT_COMPACT_OBJECT_FORMATION:           value = SN_TotalRadiusAtCOFormation();                          break;
        case ANY_STAR_PROPERTY::TRUE_ANOMALY:                                       value = SN_TrueAnomaly();                                       break;
        case ANY_STAR_PROPERTY::TZAMS:                                              value = TZAMS() * TSOL;                                         break;
        case ANY_STAR_PROPERTY::VELOCITY_X:                                         value = VelocityX();											break;
        case ANY_STAR_PROPERTY::VELOCITY_Y:                                         value = VelocityY();											break;
        case ANY_STAR_PROPERTY::VELOCITY_Z:                                         value = VelocityZ();											break;
        case ANY_STAR_PROPERTY::ZETA_HURLEY:                                        value = CalculateZetaAdiabatic_Hurley2002(m_CoreMass);           break;
        case ANY_STAR_PROPERTY::ZETA_HURLEY_HE:                                     value = CalculateZetaAdiabatic_Hurley2002(m_HeCoreMass);         break;
        case ANY_STAR_PROPERTY::ZETA_SOBERMAN:                                      value = CalculateZetaAdiabatic_Soberman1997(m_CoreMass);                  break;
        case ANY_STAR_PROPERTY::ZETA_SOBERMAN_HE:                                   value = CalculateZetaAdiabatic_Soberman1997(m_HeCoreMass);                break;

        default:                                                                                                        // unexpected stellar property
            // the only ways this can happen are if someone added a stellar property (into ANY_STAR_PROPERTY),
            // or allowed users to specify a stellar property (via the logfile definitions file), and it isn't
            // accounted for in this code.  We should not default here, with or without a warning - this is a
            // code defect, so we flag it as an error and that will result in termination of the evolution of
            // the star or binary.
            // The correct fix for this is to add code for the missing property, or prevent it from being 
            // specified in the logfile definitions file.

            THROW_ERROR(ERROR::UNEXPECTED_STELLAR_PROPERTY);                                                            // throw error
    }

    return value;
}


/*
 * PropertyValue
 *
 * Determine the value of the requested property of the star (parameter p_Property)
 *
 * The property is a boost variant variable, and is one of the following types:
 *
 *      STAR_PROPERTY           - any individual star property
 *      STAR_1_PROPERTY         - property of the primary (m_Star1)
 *      STAR_2_PROPERTY         - property of the secondary (m_Star2)
 *      SUPERNOVA_PROPERTY      - property of the star that has gone supernova
 *      COMPANION_PROPERTY      - property of the companion to the supernova
 *      BINARY_PROPERTY         - property of the binary
 *      PROGRAM_OPTION          - program option
 *
 * This function calls the appropriate helper function to retrieve the value.
 * 
 * This function handles properties of type:
 * 
 *    STAR_PROPERTY, PROGRAM_OPTION
 * 
 * only - anything else will result in an error being thrown and the evolution of the star (or binary)
 * terminated.
 *
 * This is the function used to retrieve values for properties required to be printed.
 * This allows the composition of the log records to be dynamically modified - this is
 * how we allow users to specify what properties they want recorded in log files.
 *
 * The functional return is the value of the property requested.  
 *
 *
 * COMPAS_VARIABLE PropertyValue(const T_ANY_PROPERTY p_Property) const
 *
 * @param   [IN]    p_Property                  The property for which the value is required
 * @return                                      The value of the requested property
 */
COMPAS_VARIABLE BaseStar::PropertyValue(const T_ANY_PROPERTY p_Property) const {

    COMPAS_VARIABLE value = 0.0;                                                                                                // default property value

    switch (boost::apply_visitor(VariantPropertyType(), p_Property)) {                                                          // which property type?

        case ANY_PROPERTY_TYPE::T_STAR_PROPERTY:  value = StellarPropertyValue(p_Property); break;                              // star property
        case ANY_PROPERTY_TYPE::T_PROGRAM_OPTION: value = OPTIONS->OptionValue(p_Property); break;                              // program option

        default:                                                                                                                // unexpected property type
            // the only ways this can happen are if someone added a stellar type property (into ANY_PROPERTY_TYPE)
            // and it isn't accounted for in this code, or if there is a defect in the code that causes
            // this function to be called with a bad parameter.  We should not default here, with or without a
            // warning - this is a code defect, so we flag it as an error and that will result in termination of
            // the evolution of the star or binary.
            // The correct fix for this is to add code for the missing property type or find and fix the code defect.

            THROW_ERROR(ERROR::UNEXPECTED_STELLAR_PROPERTY_TYPE);                                                               // throw error
    }

    return value;
}



















///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         AGE / LIFETIME / TAU / TIMESCALES                         //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      ENERGY                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateBindingEnergy
 *
 * @brief
 * Calculate the absolute value of the binding energy of the envelope of the star
 *
 *
 * double CalculateBindingEnergy(const double p_Radius, const double p_CoreMass, const double p_EnvMass, const double p_Lambda)
 *
 * @param           p_CoreMass                  Mass of the star's core (Msol)
 * @param           p_Radius                    Radius of the star (Rsol)
 * @param           p_EnvMass                   Mass of the star's envelope (Msol)
 * @param           p_Lambda                    Dimensionless parameter defining the binding energy
 * @return                                      Binding energy (erg)
 */
double BaseStar::CalculateBindingEnergy(const double p_Radius, const double p_CoreMass, const double p_EnvMass, const double p_Lambda) const {

    // convert to CGS
    const double Mc   = p_CoreMass * MSOL_TO_G;
    const double Menv = p_EnvMass * MSOL_TO_G;

	return G_CGS * (Mc + Menv) * Menv / (p_Lambda * p_Radius * RSOL_TO_CM);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLuminosityAtZAMS_Tout1996
 *
 * @brief
 * Calculate the luminosity of a star at ZAMS, per Tout et al. 1996, eq 1
 *
 *
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function relies on the values of the ZAMS mass and ZAMS radius of the star,
 * and should not be used if the ZAMS mass or the ZAMS radius is not known.
 * 
 * *Caveat*: we call this function after a MS merger to determine the luminosity of the
 * merger product - there we reset the effective intial mass of the merger product to
 * the final mass of the merger product, and use that as the ZAMS mass.
 * 
 * double CalculateLuminosityAtZAMS_Tout1996(const double p_MZAMSs) const
 * 
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS luminosity (Lsol)
 */
COMPAS_PURE double BaseStar::CalculateLuminosityAtZAMS_Tout1996(const double p_MZAMS) const {
#define lCoeffs(x) lCoeffs[static_cast<int>(TOUT_L_Coeff::x)] // for convenience and readability - undefined at end of function
        
    lCoeffs = GLOBALS->ToutZAMSLuminosityCoefficients(); // get Tout ZAMS luminosity coefficients

    // calculate some powers of p_MZAMS - for performance and readability
    // pow() is slow - use multiplication where it makes sense
    const double m0_5 = std::sqrt(p_MZAMS);  // sqrt() is much faster than pow()
    const double m2   = p_MZAMS * p_MZAMS;
    const double m3   = p_MZAMS * m2;
    const double m5   = m2 * m3;
    const double m7   = m2 * m5;
    const double m8   = p_MZAMS * m7;
            
    const double top  = (lCoeffs(ALPHA) * (m5 * m0_5)) + (lCoeffs(BETA) * (m3 * m8));
            
    return top / (lCoeffs(GAMMA) + m3) + (lCoeffs(DELTA) * m5) + (lCoeffs(EPSILON) * m7) + (lCoeffs(ZETA) * m8) + (lCoeffs(ETA) * (p_MZAMS * m8 * m0_5));
            
#undef lCoeffs
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     MASS LOSS                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateMassLossRate
 *
 * @brief
 * Calculate mass loss rate, and dominant mass loss type, at the current evolutionary phase.
 *
 * Calls relevant mass loss function based on mass loss prescription given in program options.
 * 
 * The first level of program option used to determine the process flow for mass loss calculations
 * is `--mass-loss-prescription`.  Based on the value specified by the user, we then drill-down to
 * more specific options to drive the mass loss process flow (e.g. `--LBV-mass-loss-prescription`,
 * `--VMS-mass-loss-prescription`, etc.).
 * 
 * Since this is the driver function for winds mass loss rate calculations, we gather the state
 * variables, option values, and whatever else is needed here and pass them to functions called
 * from here, so that this function is telling the called functions what to do, rather than the
 * called functions interrogating state and options to determine what they should do.  The only
 * exception to this is the `CalculateMLrate_Merritt2025()` function which we pass variables to
 * as parameters, but we allow the function itself to determine which option values will be used
 * (because it is in some ways a driver function for winds mass loss rate calculations).
 * 
 * 
 * MASS_LOSS_T CalculateMassLossRate() const
 *
 * @return                                      Tuple containing:
 *                                                  DOUBLE         mass loss rate (Msol yr^-1)
 *                                                                 will be clamped to [0.0, MAXIMUM_WIND_MASS_LOSS_RATE]
 *                                                  MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
MASS_LOSS_T BaseStar::CalculateMassLossRate() const {

    // *Ilya* check this please - mStart is used in:
    //    - BaseStar::CalculateMLrate_Merritt2025()
    //    - BaseStar::CalculateMLrateRSG()
    //    - CalculateMLrateRSG_Beasor2020()
    //    - CalculateMassLossRateRSG_Decin2023()
    // All were originally ZAMS mass
    // We can easily get ZAMS mass if it exists (i.e. we started on MS), but we need to manage if it doesn't

    const double mass               = Mass();                 // current mass of the star
    const double mStart             = m_StateHistory.StartState().Mass();                   // mass of the star at the start of the simulation
    const double radius             = Radius();               // current radius of the star
    const double luminosity         = Luminosity();           // current luminosity of the star
    const double temperature        = Temperature();          // current temperature of the star
    const double perturbationMu     = PerturbationMu();       // current small envelope perturbation parameter
    const double HeAbundanceSurface = HeAbundanceSurface();   // He abundance on the surface of the star

    double dMdt;
    MASS_LOSS_TYPE dominantMLtype;

    switch (OPTIONS->MassLossPrescription()) {                                              // which mass loss prescription?

        case MASS_LOSS_PRESCRIPTION::BELCZYNSKI2010:                                        // BELCZYNSKI2010
            std::tie(dMdt, dominantMLtype) = CalculateMLrate_Belczynski2010(mass, radius, luminosity, temperature, perturbationMu, HeAbundanceSurface);
            break;

        case MASS_LOSS_PRESCRIPTION::HURLEY:                                                // HURLEY
            std::tie(dMdt, dominantMLtype) = CalculateMLrate_Hurley2000(mass, radius, luminosity, perturbationMu);

            double dMdtLBV;
            MASS_LOSS_TYPE dominantMLtypeLBV;
            std::tie(dMdtLBV, dominantMLtypeLBV) = CalculateMLrateLBV(radius, luminosity, LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD);

            if (dMdtLBV > dMdt) dominantMLtype = dominantMLtypeLBV;                         // dominant ML type
            dMdt += dMdtLBV;                                                                // sum rates
            break;

        case MASS_LOSS_PRESCRIPTION::MERRITT2025:                                           // MERRITT2025
            std::tie(dMdt, dominantMLtype) = CalculateMLrate_Merritt2025(mass, radius, luminosity, temperature, perturbationMu, mStart, HeAbundanceSurface);
            break;

        case MASS_LOSS_PRESCRIPTION::ZERO:                                                  // ZERO
            dMdt = 0.0;                                                                     // no mass loss
            dominantMLtype = MASS_LOSS_TYPE::NONE;
            break;

        default:                                                                            // unknown prescription
            // the only way this can happen is if someone added a MASS_LOSS_PRESCRIPTION/ and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MASS_LOSS_PRESCRIPTION);                             // throw error
    }

    // apply overall wind mass loss multiplier and clamp winds to [0.0, MAXIMUM_WIND_MASS_LOSS_RATE]
    // to avoid convergence issues (maximum is typically 0.1 solar masses per year)    
    return std::make_tuple(std::max(std::min(dMdt * OPTIONS->OverallWindMassLossMultiplier(), MAXIMUM_WIND_MASS_LOSS_RATE), 0.0), dominantMLtype);
}


/*
 * CalculateMLrateLBV
 *
 * @brief
 * Calculate the LBV-like mass loss rate, and the dominant mass loss type, for stars
 * beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994), based on the
 * LBV_MASS_LOSS_PRESCRIPTION passed in p_LBVprescription. 
 * 
 * The LBV mass loss prescription to be used is passed as a parameter, rather than this
 * function using any option value specified by the user, so that the caller can choose
 * the prescription to be used.
 *  
 *  
 * MASS_LOSS_T CalculateMLrateLBV(const double p_Radius, const double p_Luminosity, const LBV_MASS_LOSS_PRESCRIPTION p_LBV_MLprescription) const
 *
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_LBV_MLprescription            LBV mass loss prescription to use
 * @return                                      Tuple containing:
 *                                                  DOUBLE         LBV-like mass loss rate (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrateLBV(const double p_Radius, const double p_Luminosity, const LBV_MASS_LOSS_PRESCRIPTION p_LBV_MLprescription) const {

    double dMdt = 0.0;                                                          // default mass loss rate
    MASS_LOSS_TYPE dominantMLtype = MASS_LOSS_TYPE::NONE;                       // default dominant mass loss type

    const double HDlimitfactor = p_Radius * std::sqrt(p_Luminosity) * 1.0E-5;   // factor by which the star is above the HD limit
    if (p_Luminosity > LBV_LUMINOSITY_LIMIT_STARTRACK && HDlimitfactor > 1.0) { // LBV?
        
        switch (p_LBV_MLprescription) {                                         // which LBV mass loss prescription?
            
            case LBV_MASS_LOSS_PRESCRIPTION::BELCZYNSKI:                        // BELCZYNSKI
                std::tie(dMdt, dominantMLtype) = CalculateMLrateLBV_Belczynski2010();
                break;

            case LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD:                        // HURLEY_ADD
            case LBV_MASS_LOSS_PRESCRIPTION::HURLEY:                            // HURLEY
                std::tie(dMdt, dominantMLtype) = CalculateMLrateLBV_Hurley2000(p_Luminosity, HDlimitfactor);
                break;

            case LBV_MASS_LOSS_PRESCRIPTION::ZERO:                              // ZERO
                dMdt = 0.0;                                                     // no mass loss
                dominantMLtype = MASS_LOSS_TYPE::NONE;                
            break;

            default:                                                            // unexpected prescription
                // the only way this can happen is if the LBV_MASS_LOSS_PRESCRIPTION passed to this function
                // is not accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the code passed a prescription that this function doesn't account
                // for, and that should be flagged as an error and result in termination of the evolution
                // of the star or binary.
                // The correct fix for this is to add code to this function for the missing prescription,
                // or fix the calling code to pass a prescription that is handled by this function.
                THROW_ERROR(ERROR::UNEXPECTED_LBV_MASS_LOSS_PRESCRIPTION);      // throw error
        }
    }

    // NOTE: CALLER SHOULD SET m_LBVphaseFlag BASED ON dominantMLtype - LBV vs NONE <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrateOB
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for main sequence stars,
 * based on the OB_MASS_LOSS_PRESCRIPTION passed in p_MassLossPrescription. 
 * 
 * The OB mass loss prescription to be used is passed as a parameter, rather than this
 * function using any option value specified by the user, so that the caller can choose
 * the prescription to be used.
 * 
 * 
 * MASS_LOSS_T CalculateMLrateOB(const double p_Mass, const double p_Luminosity, const double p_Temperature, const OB_MASS_LOSS_PRESCRIPTION p_OB_MLprescription) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_OB_MLprescription             OB mass loss prescription to use          
 * @return                                      Tuple containing:
 *                                                  DOUBLE         OB mass loss rate (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrateOB(const double p_Mass, const double p_Luminosity, const double p_Temperature, const OB_MASS_LOSS_PRESCRIPTION p_OB_MLprescription) const {

    double dMdt;                      
    MASS_LOSS_TYPE dominantMLtype;
    
    switch (p_OB_MLprescription) {                                      // which OB mass loss prescription?

        case OB_MASS_LOSS_PRESCRIPTION::BJORKLUND2022:                  // BJORKLUND2022
            std::tie(dMdt, dominantMLtype) = CalculateMLrateOB_Bjorklund2022(p_Mass, p_Luminosity, p_Temperature);
            break;

        case OB_MASS_LOSS_PRESCRIPTION::KRTICKA2018:                    // KRTICKA2018
            std::tie(dMdt, dominantMLtype) = CalculateMLrateOB_Krticka2018(p_Luminosity);
            break;

        case OB_MASS_LOSS_PRESCRIPTION::VINK2001:                       // VINK2001
            std::tie(dMdt, dominantMLtype) = CalculateMLrateOB_Vink2001(p_Mass, p_Luminosity, p_Temperature);
            break;

        case OB_MASS_LOSS_PRESCRIPTION::VINK2021:                       // VINK2021 (Vink & Sander 2021)
            std::tie(dMdt, dominantMLtype) = CalculateMLrateOB_VinkSander2021(p_Mass, p_Luminosity, p_Temperature);
            break;

        case OB_MASS_LOSS_PRESCRIPTION::ZERO:                           // ZERO
            dMdt = 0.0;                                                 // no mass loss
            dominantMLtype = MASS_LOSS_TYPE::NONE;                
            break;

        default:                                                        // unexpected prescription
            // the only way this can happen is if the OB_MASS_LOSS_PRESCRIPTION passed to this function
            // is not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.
            THROW_ERROR(ERROR::UNEXPECTED_OB_MASS_LOSS_PRESCRIPTION);   // throw error
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrateOB_Bjorklund2022
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for massive OB stars,
 * per Bjorklund et al. 2022.
 * 
 * See eq 7 and surrounding text in https://arxiv.org/abs/2203.08218
 * 
 * This prescription is calibrated to the following ranges and values:
 * 
 *    10^4.5 < L / Lsol < 10^6
 *    15,000 < Teff / K < 50,000
 *    15 < M / Msol < 80
 * 
 *    Zsol = 0.014
 *    Zlmc = Zsol / 2.0
 *    Zsmc = Zsol / 5.0
 *
 * 
 * MASS_LOSS_T CalculateMLrateOB_Bjorklund2022(const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  DOUBLE         Mass loss rate for massive stars (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::OB)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrateOB_Bjorklund2022(const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    const double gamma   = (p_Luminosity * LSOLW) / CalculateEddingtonLuminosity(p_Mass, 0.1);  // Bjorklund et al. 2022, para 3, assumes He abundance = 0.1
    const double logZ    = log10(GLOBALS->Metallicity() / 0.014);                      // Bjorklund et al. 2022 uses 0.014
    const double logL    = log10(p_Luminosity / 1.0E6);
    const double logTeff = log10(p_Temperature * TSOL / 45000.0);           
    const double logMeff = log10(p_Mass * (1.0 - gamma) / 45.0);

    // eq 7 in Bjorklund et al. 2022
    return std::make_tuple(PPOW(10.0, -5.52 + (2.39 * logL) + (-1.48 * logMeff) + (2.12 * logTeff) + ((0.75 - (1.87 * logTeff)) * logZ)), MASS_LOSS_TYPE::OB);
}


/*
 * CalculateMLrateOB_Vink2001
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for massive OB stars,
 * per Vink et al. 2001.
 * 
 * See:
 *    - Vink et al. 2001, eqs 24 & 25
 *    - Belczynski et al. 2010, eqs 6 & 7
 *
 * 
 * MASS_LOSS_T CalculateMLrateOB_Vink2001(const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  DOUBLE         Mass loss rate for hot OB stars (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrateOB_Vink2001(const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    double dMdt;
    MASS_LOSS_TYPE dominantMLtype;

    const double teff = p_Temperature * TSOL;                                                                           // Kelvin

    if (teff >= VINK_MASS_LOSS_MINIMUM_TEMP && teff <= VINK_MASS_LOSS_BISTABILITY_TEMP) {
        const double v = 1.3 * PPOW(GLOBALS->ZetaAnders(), OPTIONS->ScaleTerminalWindVelocityWithMetallicityPower());   // v_inf / v_esc, scaled with metallicity

        dMdt = PPOW(10.0, -6.688 + (2.210 * log10(p_Luminosity / 1.0E5)) - (1.339 * log10(p_Mass / 30.0)) - (1.601 * log10(v / 2.0)) + (0.85 * GLOBALS->ZetaAnders()) + (1.07 * log10(teff / 20000.0)));
        dominantMLtype = MASS_LOSS_TYPE::OB; 
    }
    else if (teff > VINK_MASS_LOSS_BISTABILITY_TEMP) {
        const double v1 = 2.6 * PPOW(GLOBALS->ZetaAnders(), OPTIONS->ScaleTerminalWindVelocityWithMetallicityPower());  // v_inf / v_esc, scaled with metallicity
        const double v2 = log10(teff / 40000.0)

        dMdt = PPOW(10.0, -6.697 + (2.194 * log10(p_Luminosity / 1.0E5)) - (1.313 * log10(p_Mass / 30.0)) - (1.226 * log10(v1 / 2.0)) + (0.85 * GLOBALS->ZetaAnders()) + (0.933 * v2) - (10.92 * v2 * v2));
        dominantMLtype = MASS_LOSS_TYPE::OB; 

        SHOW_WARN_IF(teff > VINK_MASS_LOSS_MAXIMUM_TEMP, ERROR::HIGH_TEFF_WINDS);                                       // show warning if winds being used outside comfort zone
    }
    else {                                                                                                              // too cold to use winds
        dMdt = 0.0;                                                                                                     // turn winds off
        dominantMLtype = MASS_LOSS_TYPE::NONE;                                                                          // set dominant type (NONE)
        SHOW_WARN(ERROR::LOW_TEFF_WINDS, "Mass Loss Rate = 0.0");                                                       // show warning
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrateOB_VinkSander2021
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for massive OB stars,
 * per Vink and Sander 2021.
 * 
 * See https://arxiv.org/pdf/2103.12736.pdf
 * 
 * Features two bistability jumps, at T1 and T2
 * offset = {"cold":-5.99, "inter":-6.688, "hot":-6.697}
 *
 * 
 * MASS_LOSS_T CalculateMLrateOB_VinkSander2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  DOUBLE         Mass loss rate for hot OB stars (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrateOB_VinkSander2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    double dMdt;
    MASS_LOSS_TYPE dominantMLtype;

    const double teff = p_Temperature * TSOL;                                           // Kelvin

    if (teff < VINK_MASS_LOSS_MINIMUM_TEMP) {                                           // temp below Vink minimum?
                                                                                        // yes - too cold to use winds
        dMdt = 0.0;                                                                     // turn winds off
        dominantMLtype = MASS_LOSS_TYPE::NONE;

        SHOW_WARN(ERROR::LOW_TEFF_WINDS, "Mass Loss Rate = 0.0");                       // show warning
    }
    else {                                                                              // no - temp is at or above minimum
        dominantMLtype = MASS_LOSS_TYPE::OB;                                            // OB winds

        const double logL5  = log10(p_Luminosity / 1.0E5);                              // common value used below
        const double logM30 = log10(p_Mass / 30.0);                                     // common value used below

        constexpr double Zexp2001 = 0.85;                                               // Vink et al. 2001
    
        const double gamma = EDDINGTON_PARAMETER_FACTOR * p_Luminosity / p_Mass;
        const double rho   = -14.94 + (3.1857 * gamma) + (zExp * GLOBALS->ZetaAnders());// characteristic density
        const double T1    = ( 100.0 + (6.0 * rho) ) * 1000.0;                          // bistability jump 1: typically around 20000.0, has similar behavior when fixed

        if (teff <= T1) {                                                               // temp at or below jump 1?
                                                                                        // yes
            constexpr double v = 0.7;                                                   // v_inf / v_esc

            dMdt = PPOW(10.0, -5.99 + (2.210 * logL5) - (1.339 * logM30) - (1.601 * log10(v / 2.0)) + (Zexp2001 * GLOBALS->ZetaAnders()) + (1.07 * log10(teff / 20000.0)));
        }
        else {                                                                          // temp is above jump 1
            SHOW_WARN_IF(teff > VINK_MASS_LOSS_MAXIMUM_TEMP, ERROR::HIGH_TEFF_WINDS);   // show warning if temp above Vink maximum

            const double T2    = ( 61.2 + (2.59 * rho) ) * 1000.0;                      // bistability jump 2: typically around 25000.0, higher jump first as in Vink python recipe

            if (teff <= T2) {                                                           // temp at or below jump 2?
                                                                                        // yes
                constexpr double v = 1.3;                                               // v_inf / v_esc

                dMdt = PPOW(10.0, -6.688 + (2.210 * logL5) - (1.339 * logM30) - (1.601 * log10(v / 2.0)) + (Zexp2001 * GLOBALS->ZetaAnders()) + (1.07 * log10(teff / 20000.0)));
            }
            else {                                                                      // temp is above jump 2                      
                constexpr double Zexp2021 = 0.42;                                       // Vink and Sander 2021
                constexpr double v        = 2.6;                                        // v_inf / v_esc
                const     double logT40   = log10(teff / 40000.0);

                dMdt = PPOW(10.0, -6.697 + (2.194 * logL5) - (1.313 * logM30) - (1.226 * log10(v / 2.0)) + (Zexp2021 * GLOBALS->ZetaAnders()) + (0.933 * logT40) - (10.92 * logT40 * logT40));
            }
        }
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrateRSG
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for RSG stars (Red Supergiants),
 * based on the RSG_MASS_LOSS_PRESCRIPTION passed in p_MassLossPrescription. 
 * 
 * The RSG mass loss prescription to be used is passed as a parameter, rather than this
 * function using any option value specified by the user, so that the caller can choose
 * the prescription to be used.
 * 
 * 
 * MASS_LOSS_T CalculateMLrateRSG(
 *     const double                     p_Mass,
 *     const double                     p_Radius,
 *     const double                     p_Luminosity,
 *     const double                     p_Temperature,
 *     const double                     p_mStart,
 *     const RSG_MASS_LOSS_PRESCRIPTION p_RSG_MLprescription
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_RSG_MLprescription            RSG Mass loss prescription
 * @return                                      Tuple containing:
 *                                                  DOUBLE         RSG mass loss rate (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrateRSG(
    const double                     p_Mass,
    const double                     p_Radius,
    const double                     p_Luminosity,
    const double                     p_Temperature,
    const double                     p_mStart,
    const RSG_MASS_LOSS_PRESCRIPTION p_RSG_MLprescription
) const {

    double dMdt;                      
    MASS_LOSS_TYPE dominantMLtype;

    switch (p_RSG_MLprescription) {                                     // which RSG mass loss prescription?

        case RSG_MASS_LOSS_PRESCRIPTION::BEASOR2020:                    // BEASOR2020
            std::tie(dMdt, dominantMLtype) = CalculateMLrateRSG_Beasor2020(p_mStart, p_Luminosity);
            break;

        case RSG_MASS_LOSS_PRESCRIPTION::DECIN2023:                     // DECIN2023
            std::tie(dMdt, dominantMLtype) = CalculateMLrateRSG_Decin2023(p_mStart, p_Luminosity);
            break;

        case RSG_MASS_LOSS_PRESCRIPTION::KEE2021:                       // KEE2021
            std::tie(dMdt, dominantMLtype) = CalculateMLrateRSG_Kee2021(p_Mass, p_Luminosity, p_Temperature);
            break;

        case RSG_MASS_LOSS_PRESCRIPTION::NJ90:                          // NJ90
            std::tie(dMdt, dominantMLtype) = CalculateMLrate_NieuwenhuijzenDeJager1990(p_Mass, p_Radius, p_Luminosity);
            break;

        case RSG_MASS_LOSS_PRESCRIPTION::VINKSABHAHIT2023:              // VINKSABHAHIT2023
            std::tie(dMdt, dominantMLtype) = CalculateMLrateRSG_VinkSabhahit2023(p_Mass, p_Luminosity);
            break;

        case RSG_MASS_LOSS_PRESCRIPTION::YANG2023:                      // YANG2023
            std::tie(dMdt, dominantMLtype) = CalculateMLrateRSG_Yang2023(p_Luminosity);
            break;    

        case RSG_MASS_LOSS_PRESCRIPTION::ZERO:                          // ZERO
            dMdt = 0.0;                                                 // no mass loss
            dominantMLtype = MASS_LOSS_TYPE::NONE;
            break;

        default:                                                        // unexpected prescription
            // the only way this can happen is if the RSG_MASS_LOSS_PRESCRIPTION passed to this function
            // is not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.
            THROW_ERROR(ERROR::UNEXPECTED_RSG_MASS_LOSS_PRESCRIPTION);  // throw error
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrateRSG_Kee2021
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for Red Supergiant (RSG) stars,
 * per Kee et al. 2021, eqs 5, 13, 14, 25.
 * 
 * See https://arxiv.org/pdf/2101.03070.pdf  
 * 
 * 
 * MASS_LOSS_T CalculateMLrateRSG_Kee2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  DOUBLE         Mass loss rate for RSG stars (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::RSG)
 */
GNU_CONST MASS_LOSS_T BaseStar::CalculateMLrateRSG_Kee2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    constexpr double vTurb  = 1.5E4;                                                                        // turbulent velocity, ms^-1, for a typical RSG
    constexpr double kBoltz = 1.38E-23;                                                                     // Boltzmann Constant in J K^-1
    constexpr double sigma  = 5.67E-8;                                                                      // Stefan Boltzmann constant W m^-2 K^-4
    constexpr double mH     = 1.67E-27;                                                                     // mass of hydrogen in Kg
    constexpr double kappa  = 0.01 * OPACITY_CGS_TO_SI;                                                     // given after eq 16 

    const double tEff  = TSOL * p_Temperature;                                                              // Kelvin
    
    const double rSI   = std::sqrt((p_Luminosity * LSOLW) / (4.0 * M_PI * sigma * intPow(teff, 4)));        // SI units
    const double mSI   = p_Mass * MSOL_TO_KG;                                                               // SI units
    const double cs    = std::sqrt(kBoltz * tEff / mH);
    const double gamma = (kappa * p_Luminosity * LSOLW) / (4.0 * M_PI * G * C * mSI);
    const double vEsc  = std::sqrt(2.0 * G * (mSI) / (rSI));                                                // ms^-1, not vEsc,eff

    const double rPmod = G * (M_SI) * (1.0 - gamma) / (2.0 * ((cs * cs) + (vTurb * vTurb)));                // modified parker radius, in m
    const double rho   = (4.0 / 3.0) * (rPmod / (kappa * (rSI) * (rSI))) * (std::exp(-(2.0 * rPmod / (rSI)) + (3.0 / 2.0))) / (1.0 - std::exp(-2.0 * rPmod / (rSI)));

    const double mDotAnalytical = 4.0 * M_PI * rho * std::sqrt(cs * cs + vTurb * vTurb) * rPmod * rPmod;    // kg s^-1

    return std::make_tuple(PPOW(((vTurb / 17000.0) / (vEsc / 60000.0)), 1.30) * mDotAnalytical * SECONDS_IN_YEAR / MSOL_TO_KG, MASS_LOSS_TYPE::RSG);
}   


/*
 * CalculateMLrateVMS
 *
 * Calculate mass loss rate, and the dominant mass loss type, for very massive main sequence
 * stars (> 100 Msol), based on the VMS_MASS_LOSS_PRESCRIPTION passed in p_MassLossPrescription. 
 * 
 * The VMS mass loss prescription to be used is passed as a parameter, rather than this
 * function using any option value specified by the user, so that the caller can choose
 * the prescription to be used.
 * 
 * 
 * MASS_LOSS_T CalculateMLrateVMS(const double p_Mass, const double p_Luminosity, const double p_Temperature, const VMS_MASS_LOSS_PRESCRIPTION p_VMS_MLprescription) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_VMS_MLprescription            VMS Mass loss prescription to use
 * @return                                      Tuple containing:
 *                                                  DOUBLE         VMS mass loss rate (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrateVMS(const double p_Mass, const double p_Luminosity, const double p_Temperature, const VMS_MASS_LOSS_PRESCRIPTION p_VMS_MLprescription) const {

    double dMdt;                         
    MASS_LOSS_TYPE dominantMLtype;

    switch (p_VMS_MLprescription) {                                     // which VMS mass loss prescription?

        case VMS_MASS_LOSS_PRESCRIPTION::BESTENLEHNER2020:              // BESTENLEHNER2020
            std::tie(dMdt, dominantMLtype) = CalculateMLrateVMS_Bestenlehner2020(p_Mass, p_Luminosity);
            break;
        
        case VMS_MASS_LOSS_PRESCRIPTION::SABHAHIT2023:                  // SABHAHIT2023
            std::tie(dMdt, dominantMLtype) = CalculateMLrateVMS_Sabhahit2023(p_Mass, p_Luminosity, p_Temperature);
            break;
        
        case VMS_MASS_LOSS_PRESCRIPTION::VINK2011:                      // VINK2011
            std::tie(dMdt, dominantMLtype) = CalculateMLrateVMS_Vink2011(p_Mass, p_Luminosity, p_Temperature);
            break;

        case VMS_MASS_LOSS_PRESCRIPTION::ZERO:                          // ZERO
            dMdt = 0.0;                                                 // no mass loss
            dominantMLtype = MASS_LOSS_TYPE::NONE;
        break;

        default:                                                        // unexpected prescription
            // the only way this can happen is if the VMS_MASS_LOSS_PRESCRIPTION passed to this function
            // is not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.
            THROW_ERROR(ERROR::UNEXPECTED_VMS_MASS_LOSS_PRESCRIPTION);  // throw error
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrateVMS_Sabhahit2023
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for Very Massive Stars (VMS),
 * per Sabhahit 2023.
 * 
 * See https://arxiv.org/pdf/2306.11785.pdf
 *
 * 
 * MASS_LOSS_T CalculateMLrateVMS_Sabhahit2023(const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  DOUBLE         Mass loss rate for very massive stars (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::VMS or MASS_LOSS_TYPE::OB)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrateVMS_Sabhahit2023(const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    double dMdt;                        
    MASS_LOSS_TYPE dominantMLtype;

    const double mSwitch = PPOW(GLOBALS->Metallicity(), -1.574) * 0.0615 + 18.10;  // obtained from a powerlaw fit to Sabhahit 2023, table 2, given teff = 45kK
    const double lSwitch = PPOW(10.0, (-1.91 * GLOBALS->SigmaHurley() + 2.36));             // loglinear fits to Sabhahit 2023, table 2 

    if ((p_Luminosity / p_Mass) > (lSwitch / mSwitch)) {                                    // in the VMS regime according to Sabhahit+ 2023?
                                                                                            // yes
        dMdt = PPOW(10.0, (-1.86 * GLOBALS->SigmaHurley() - 8.90)) * PPOW(p_Luminosity / lSwitch , 4.77) * PPOW(p_Mass / mSwitch, -3.99);
        dominantMLtype = MASS_LOSS_TYPE::VMS;                                               // dominant mass loss type is VMS

    }
    else {                                                                                  // no, fall back to default OB mass loss prescription
        std::tie(dMdt, dominantMLtype) = CalculateMLrateOB(GLOBALS->Metallicity(), p_Mass, p_Luminosity, p_Temperature, OPTIONS->OBMassLossPrescription());
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrateVMS_Vink2011
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for Very Massive (OB) Stars (VMS),
 * per Vink 2011.
 * 
 * See https://arxiv.org/pdf/1105.0556.pdf
 *
 * 
 * MASS_LOSS_T CalculateMLrateVMS_Vink2011(p_Mass, p_Luminosity, p_Temperature) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  DOUBLE         Mass loss rate for very massive stars (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (will be MASS_LOSS_TYPE::VMS or MASS_LOSS_TYPE::OB)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrateVMS_Vink2011(const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    double dMdt;
    MASS_LOSS_TYPE dominantMLtype;

    // start with rate for massive OB stars, per Vink et al. 2001
    std::tie(dMdt, dominantMLtype) = CalculateMLrateOB_Vink2001(p_Mass, p_Luminosity, p_Temperature);

    const double gamma = EDDINGTON_PARAMETER_FACTOR * p_Luminosity / p_Mass;    // Eddington parameter, independent of surface composition

    if (gamma > 0.5) {                                                          // apply correction to high gamma only
        dMdt           = PPOW(10.0, (0.04468 + (0.3091 * gamma) + (0.2434 * gamma * gamma) + log10(dMdt)));
        dominantMLtype = MASS_LOSS_TYPE::VMS;
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrate_Belczynski2010
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, per Belczynski 2010
 * (as implemented in StarTrack - courtesy Chris Belczynski).
 *
 * If option `--scale-mass-loss-with-surface-helium-abundance` was specified,
 * the mass loss rate will be scaled with the surface helium abundance.
 * 
 * 
 * MASS_LOSS_T CalculateMLrate_Belczynski2010(
 *     const double p_Mass,
 *     const double p_Radius,
 *     const double p_Luminosity,
 *     const double p_Temperature,
 *     const double p_PerturbationMu,
 *     const double p_HeAbundanceSurface,
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Tuple containing:
 *                                                  DOUBLE         mass loss rate (Msol yr^-1)
 *                                                  MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrate_Belczynski2010(
    const double p_Mass,
    const double p_Radius,
    const double p_Luminosity,
    const double p_Temperature,
    const double p_PerturbationMu,
    const double p_HeAbundanceSurface
) const {

    double dMdt;  
    MASS_LOSS_TYPE dominantMLtype;

    // start with LBV winds (can be, and is often, 0.0)
    std::tie(dMdt, dominantMLtype) = CalculateMLrateLBV(p_Radius, p_Luminosity, OPTIONS->LBVMassLossPrescription());

    // other winds - if not in LBV regime, or user specified HURLEY_ADD 
    if (dominantMLtype != MASS_LOSS_TYPE::LBV || OPTIONS->LBVMassLossPrescription() == LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD) {

        double dMdtOther;
        MASS_LOSS_TYPE dominantMLtypeOther;

        if (p_Temperature * TSOL < VINK_MASS_LOSS_MINIMUM_TEMP) {                   // cool star?
                                                                                    // yes
            //add Hurley et al. 2000 winds and apply cool wind mass loss multiplier
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLrate_Hurley2000(p_Mass, p_Radius, p_Luminosity, p_PerturbationMu) * OPTIONS->CoolWindMassLossMultiplier();
        }
        else  {                                                                     // no - hot star
            // add Vink et al. 2001 winds (ignoring bistability jump)
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLrateOB_Vink2001(p_Mass, p_Luminosity, p_Temperature);

            // scale mass loss with the surface helium abundance if necessary
            // (transition between OB and WR mass loss rates)
            if (OPTIONS->ScaleMassLossWithSurfaceHeliumAbundance()) {
                double dMdtWR;
                std::tie(dMdtWR, std::ignore) = CalculateMLrateWR_ZDependent_Static(GLOBALS->Metallicity(), p_Luminosity, 0.0);
                MASS_LOSS_TYPE thisDominantMLType;
                std::tie(dMdtOther, thisDominantMLType) = CalculateMLrate_WRenhanced(GLOBALS->Metallicity(), p_Luminosity, p_Temperature, p_HeAbundanceSurface, dMdtOther, dMdtWR);
                if (thisDominantMLType != MASS_LOSS_TYPE::NONE) dominantMLtypeOther = thisDominantMLType;
            }
        }

        if (dMdtOther > dMdt) dominantMLtype = dominantMLtypeOther;                 // dominant ML type
        dMdt += dMdtOther;                                                          // sum rates
    }

    // Note: BSE and StarTrack have some multiplier they apply here
    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrate_WRenhanced
 *
 * @brief
 * Calculate mass loss rate with enhancement due to WR winds (see CalculateMLfractionWR)
 *
 * 
 * double CalculateMLrate_WRenhanced(
 *     const double                p_Luminosity, 
 *     const double                p_Temperature, 
 *     const double                p_HeAbundanceSurface, 
 *     const double                p_dMdtOther,
 *     const std::optional<double> p_dMdtWR
 * ) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @param       p_dMdtOther                     Wind mass loss rate due to OB or VMS winds
 * @param       p_dMdtWR                        Wind mass loss rate due to WR winds - optional
 *                                              if p_dMdtWR has a valid value, use it, otherwise calculate it
 * @return                                      Tuple containing:
 *                                                  DOUBLE         mass loss rate (Msol yr^-1)
 *                                                                 will be clamped to [0.0, MAXIMUM_WIND_MASS_LOSS_RATE]
 *                                                  MASS_LOSS_TYPE dominant mass loss type
 *                                                                 will be MASS_LOSS_TYPE::WR if WR winds are dominant,
 *                                                                 otherwise MASS_LOSS_TYPE::NONE, indicating no change in dominance)
 */
COMPAS_PURE double BaseStar::CalculateMLrate_WRenhanced(
    const double                p_Luminosity, 
    const double                p_Temperature, 
    const double                p_HeAbundanceSurface, 
    const double                p_dMdtOther,
    const std::optional<double> p_dMdtWR
) const {
    
    const double fractionWR = CalculateMLfractionWR(p_HeAbundanceSurface);

    // set defaults
    double dMdt = (1.0 - fractionWR) * p_dMdtOther;
    MASS_LOSS_TYPE dominantMLtype = MASS_LOSS_TYPE::NONE;

    // determine WR winds impact    
    if (fractionWR > 0.0) {                             // non-zero impact?
                                                        // yes
        double dMdtWRfactor = 0.0;

        if (pdMdtWR.has_value()) dMdtWRfactor = fractionWR * pdMdtWR.value();
        else                     dMdtWRfactor = fractionWR * HeMS::CalculateMLrate_Merritt2025_Static(GLOBALS->Metallicity(), p_Luminosity, p_Temperature);
        
        if (dMdtWRfactor > dMdt) {                      // WR winds dominant?
            m_DominantMassLossRate =MASS_LOSS_TYPE::WR; // yes 
        }

        dMdt += dMdtWRfactor;                           // total rate
    }
    
    // clamp winds to [0.0, MAXIMUM_WIND_MASS_LOSS_RATE]
    return std::make_tuple(std::max(std::min(dMdt, MAXIMUM_WIND_MASS_LOSS_RATE), 0.0), dominantMLtype);
}


/*
 * CalculateMLrate_NieuwenhuijzenDeJager1990
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for massive stars (L > 4000 Lsol),
 * per Nieuwenhuijzen & de Jager 1990, modified by a metallicity dependent factor (Kudritzki et al. 1989).
 *
 * See Hurley et al. 2000, just after eq 106
 *
 *
 * MASS_LOSS_T CalculateMLrate_NieuwenhuijzenDeJager1990(const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         Mass loss rate for massive stars (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrate_NieuwenhuijzenDeJager1990(const double p_Mass, const double p_Radius, const double p_Luminosity) const {

    double dMdt = 0.0;                                                          // default mass loss rate
    MASS_LOSS_TYPE dominantMLtype = MASS_LOSS_TYPE::NONE;                       // default dominant mass loss type

    // mass loss only if star's luminosity is above minimum for Nieuwenhuijzen & de Jager 1990
    if (p_Luminosity > NJ_MINIMUM_LUMINOSITY) {
        dominantMLtype = MASS_LOSS_TYPE::GB;

        const double taper = std::min(1.0, (p_Luminosity - 4000.0) / 500.0);    // smooth taper between no mass loss and mass loss       
        dMdt = std::min(std::sqrt(GLOBALS->ZscaledHurley()) * taper * 9.6E-15 * PPOW(p_Radius, 0.81) * PPOW(p_Luminosity, 1.24) * PPOW(p_Mass, 0.16), (1.36E-9 * p_Luminosity));
    }
    
    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrateWR_SanderVink2020_Static
 *
 * @brief
 * Calculate the mass loss rate for Wolf-Rayet stars, per Sander & Vink 2020, eq 13
 * https://arxiv.org/abs/2009.01849
 * 
 * 
 * static MASS_LOSS_T CalculateMLrateWR_SanderVink2020_Static(const double p_Luminosity, const double p_PerturbationMu) const
 *
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrateWR_SanderVink2020_Static(const double p_Luminosity, const double p_PerturbationMu) const {

    double dMdt = 0.0;                                                                                  // default mass loss rate                          
    MASS_LOSS_TYPE dominantMLtype = MASS_LOSS_TYPE::NONE;                                               // default dominant mass loss type

    if (p_Mu < 1.0) {                                                                                   // small envelope?
                                                                                                        // yes
        const double logL      = log10(p_Luminosity);
        const double logL0     = -0.87 * GLOBALS->ZetaAnders() + 5.06;                                  // Sander & Vink 2020, eq 19
        const double alpha     =  0.32 * GLOBALS->ZetaAnders() + 1.4;                                   // ibid., eq 18
        const double logMdot10 = -0.75 * GLOBALS->ZetaAnders() - 4.06;                                  // ibid., eq 20

        if (logL0 <= logL) {                                                                            // no mass loss for L < L0     
            dMdt = PPOW(10.0, alpha * log10(logL - logL0) + 0.75 * (logL - logL0 - 1.0) + logMdot10);   // ibid., eq 13
            dominantMLtype = MASS_LOSS_TYPE::WR;
        }
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLrate_Merritt2025
 *
 * Calculate mass loss rate, and dominant mass loss type, at the current evolutionary phase,
 * per Merritt et al., 2025.
 *
 * Mass loss rates are divided into several classes:
 *      RSG winds,
 *      cool star winds,
 *      very massive star (VMS) winds,
 *      OB star winds
 * 
 * If option `--scale-mass-loss-with-surface-helium-abundance` was specified,
 * the mass loss rate will be scaled with the surface helium abundance.
 * 
 * 
 * MASS_LOSS_T CalculateMLrate_Merritt2025(
 *     const double p_Mass,
 *     const double p_Radius,
 *     const double p_Luminosity,
 *     const double p_Temperature,
 *     const double p_PerturbationMu,
 *     const double p_mStart,
 *     const double p_HeAbundanceSurface
 * ) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_PerturbationMu                Small envelope perturbation parameter, mu
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Tuple containing:
 *                                                   DOUBLE         WR mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)
 */
COMPAS_PURE MASS_LOSS_T BaseStar::CalculateMLrate_Merritt2025(
    const double p_Mass,
    const double p_Radius,
    const double p_Luminosity,
    const double p_Temperature,
    const double p_PerturbationMu,
    const double p_mStart,
    const double p_HeAbundanceSurface
) const {

    double dMdt;
    MASS_LOSS_TYPE dominantMLtype;

    // start with LBV winds (can be, and is often, 0.0)
    std::tie(dMdt, dominantMLtype) = CalculateMLrateLBV(p_Radius, p_Luminosity, OPTIONS->LBVMassLossPrescription());

    // calculate other winds rate if necessary
    // We add other winds to the LBV winds if the use specified LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD,
    // or if the winds are not in LBV regime
    if (OPTIONS->LBVMassLossPrescription() == LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD || dominantMLtype != MASS_LOSS_TYPE::LBV) { 

        double teff = p_Temperature * TSOL;                                 // Kelvin

        double dMdtOther;
        MASS_LOSS_TYPE dominantMLtypeOther;

        // RSG winds regime, massive star, and core helium burning giant (CHeB, FGB, EAGB, TPAGB), or HG?
        if (teff < RSG_MAXIMUM_TEMP && p_mStart >= MASSIVE_THRESHOLD && (IsOneOf(GIANTS) || IsOneOf({STELLAR_TYPE::HERTZSPRUNG_GAP}))) {
            // yes - RSG mass loss rate
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLrateRSG(p_Mass, p_Radius, p_Luminosity, p_Temperature, p_mStart, OPTIONS->RSGMassLossPrescription());
        }

        // no - cool star?
        else if (teff < VINK_MASS_LOSS_MINIMUM_TEMP) {
            // yes - HURLEY mass loss rate
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLrate_Hurley2000(p_Mass, p_Radius, p_Luminosity, p_PerturbationMu);
            dMdt *= OPTIONS->CoolWindMassLossMultiplier();                  // apply cool wind mass loss multiplier
        }

        // no - VMS winds regime?
        else if (p_Mass >= VMS_MASS_THRESHOLD) {
            // yes - VMS mass loss rate
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLrateVMS(GLOBALS->Metallicity(), p_Mass, p_Luminosity, p_Temperature, GLOBALS->ZetaAnders(), OPTIONS->ScaleTerminalWindVelocityWithMetallicityPower());

            // scale mass loss with the surface helium abundance if necessary
            // (transition between OB and WR mass loss rates)
            if (OPTIONS->ScaleMassLossWithSurfaceHeliumAbundance()) {
                MASS_LOSS_TYPE thisDominantMLType;
                std::tie(dMdtOther, thisDominantMLType) = CalculateMLrate_WRenhanced(GLOBALS->Metallicity(), p_Luminosity, p_Temperature, p_HeAbundanceSurface, dMdtOther, std::nullopt);
                if (thisDominantMLType != MASS_LOSS_TYPE::NONE) dominantMLtypeOther = thisDominantMLType;
            }
        }

        // otherwise...
        else {
            // OB mass loss rate
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLrateOB(GLOBALS->Metallicity(), p_Mass, p_Luminosity, p_Temperature, OPTIONS->OBMassLossPrescription());

            // scale mass loss with the surface helium abundance if necessary
            // (transition between OB and WR mass loss rates)
            if (OPTIONS->ScaleMassLossWithSurfaceHeliumAbundance()) {
                MASS_LOSS_TYPE thisDominantMLType;
                std::tie(dMdtOther, thisDominantMLType) = CalculateMLrate_WRenhanced(GLOBALS->Metallicity(), p_Luminosity, p_Temperature, p_HeAbundanceSurface, dMdtOther, std::nullopt);
                if (thisDominantMLType != MASS_LOSS_TYPE::NONE) dominantMLtypeOther = thisDominantMLType;
            }
        }

        if (dMdtOther > dMdt) dominantMLtype = dominantMLtypeOther;         // dominant ML type
        dMdt += dMdtOther;                                                  // sum rates
    }

    return std::make_tuple(dMdt, dominantMLtype);
}




////////////////////////////// What's this for ?????????????????????????????????????????????????????????????????????????????????????????????????
/*
 * CalculateMassLossValues
 *
 * Calculate values for mDot and mass assuming mass loss is applied
 *
 * - calculates mass loss
 * - calculates new mass loss rate (mDot) to match (possibly limited) mass loss
 * - calculates new mass (mass) based on (possibly limited) mass loss
 * - returns existing value for mass if mass loss not being used (program option)
 *
 *
 * double CalculateMassLossValues(double p_dt)
 *
 * @param   [IN]    p_dt                        time step (Myr)
 * @return                                      Tuple containing:
 *                                                   DOUBLE         mass loss rate (Msol yr^-1)
 *                                                   MASS_LOSS_TYPE dominant mass loss type (could be MASS_LOSS_TYPE::NONE)

 */
double BaseStar::CalculateMassLossValues(double p_dt) { // DONT'T NEED p_Dt HERE

// PLACEHOLDER U NTIL CODE CLEANED UP - DO SOMETHING WITH PHOTON TIRING LIMIT

//        if (OPTIONS->CheckPhotonTiringLimit()) {
//            double lim = m_Luminosity / (G_SOLAR_YEAR * m_Mass / m_Radius);     // calculate the photon tiring limit in Msol yr^-1 using Owocki & Gayley 1997, equation slightly clearer in Owocki+2004 Eq. 20
//            massLoss   = std::min(massLoss, lim);                               // limit mass loss to the photon tiring limit
//            if (p_UpdateMDot) m_Mdot = massLoss / p_dt / 1.0E6;                 // update class member variable if necessary
//        }


    return CalculateMassLossRate();
}


/*
 * ResolveMassLossHurley
 *
 * @brief
 * Calculate winds mass loss per Hurley et al., 2000
 *
 * - calculate mass loss rate
 * - calculates (and limits) mass loss
 * - resets mass loss rate (m_Mdot) to match (possibly limited) mass loss
 * - calculates and sets new mass (m_Mass) based on (possibly limited) mass loss
 * - applies mass rejuvenation factor and calculates new age
 * - updates angular momentum of mass-losing star
 *
 *
 * STELLAR_TYPE ResolveMassLossHurley(const double p_dt)
 *  
 * @param           p_dt                        time step (Myr)
 * @return                                      New stellar type for star
 */
STELLAR_TYPE BaseStar::ResolveMassLossHurley(const double p_dt) {

    STELLAR_TYPE nextStellarType = m_StellarType;                                                   // next stellar type - defaults to current

    if (OPTIONS->MassLossPrescriptionHurley() != MASS_LOSS_PRESCRIPTION::ZERO) {                    // mass loss enabled for Hurley?
                                                                                                    // yes

        double Mdot = CalculateMassLossRate();                                                      // mass loss rate (Msol yr^-1)

/// DO THIS ELSEWHERE        double angularMomentumChange = (2.0 / 3.0) * (mass - m_Mass) * m_Radius * RSOL_TO_AU * m_Radius * RSOL_TO_AU * Omega();
          
        nextStellarType = EvolveOneTimestep(mass - m_Mass, 0.0, 0.0, false);    // apply mass loss, but don't age the star yet JR FIX THIS - NEW FUNC/NAME


//        UpdateInitialMass();                                                                        // update effective initial mass (MS, HG & HeMS)
//        UpdateAgeAfterMassLoss();                                                                   // update age (MS, HG & HeMS)
//        ApplyMassTransferRejuvenationFactor();                                                      // apply age rejuvenation factor
//        SetAngularMomentum(m_AngularMomentum + angularMomentumChange);
    }

    return nextStellarType;
}























///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRadiusAtZAMS_Tout1996
 *
 * @brief
 * Calculate the radius of a star at ZAMS, per Tout et al. 1996, eq 2
 *
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function relies on the value of the ZAMS mass of the star, and should not be used
 * if the ZAMS mass is not known.
 * 
 * 
 * double CalculateRadiusAtZAMS_Tout1996(const double p_MZAMS)
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS radius of the star (Rsol)
 */
COMPAS_PURE double BaseStar::CalculateRadiusAtZAMS_Tout1996(const double p_MZAMS) const {
#define rCoeffs(x) rCoeffs[static_cast<int>(TOUT_R_Coeff::x)] // for convenience and readability - undefined at end of function

    rCoeffs = GLOBALS->ToutZAMSRadiusCoefficients(); // get Tout ZAMS radius coefficients

    // calculate some powers of p_MZAMS - for performance and readability
    // pow() is slow - use multiplication where it makes sense
    const double m0_5  = std::sqrt(p_MZAMS);
    const double m2    = p_MZAMS * p_MZAMS;
    const double m6    = m2 * m2 * m2;
    const double m8    = m6 * m2;
    const double m11   = p_MZAMS * m8 * m2;
    const double m19   = m11 * m8;
    const double m19_5 = m19 * m0_5;

    const double top = (rCoeffs(THETA) * (m2 * m0_5)) + (rCoeffs(IOTA) * m6 * m0_5) + (rCoeffs(KAPPA) * m11) + (rCoeffs(LAMBDA) * m19) + (rCoeffs(MU) * m19_5);
            
    return top / (rCoeffs(NU) + (rCoeffs(XI) * m2) + (rCoeffs(OMICRON) * (m8 * m0_5)) + (m6 * m6 * m6 * m0_5) + (rCoeffs(PI) * m19_5));

#undef rCoeffs
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ROTATION                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRotationalVelocityOStar_Ramirez2013
 *
 * @brief
 * Calculate rotational velocity from the analytic cumulative distribution function (CDF)
 * for the equatorial rotational velocity of single O stars, per Ramirez-Agudelo et al. 2013
 * (see https://arxiv.org/abs/1309.2929)
 *
 * Uses inverse sampling and root finding.
 *
 * 
 * double CalculateRotationalVelocityOStar_Ramirez2013() const
 *
 * @return                                      Rotational velocity (km s^-1)
 */
COMPAS_PURE double BaseStar::CalculateRotationalVelocityOStar_Ramirez2013() const {

    double desiredCDF            = RAND->Random();                                                  // random desired CDF

    const boost::uintmax_t maxit = ADAPTIVE_RV_MAX_ITERATIONS;                                      // limit to maximum iterations.
    boost::uintmax_t it          = maxit;                                                           // initially our chosen max iterations, but updated with actual

    // find root
    // we use an iterative algorithm to find the root here:
    //    - if the root finder throws an exception, we stop and return a negative value for the root (indicating no root found)
    //    - if the root finder reaches the maximum number of (internal) iterations, we stop and return a negative value for the root (indicating no root found)
    //    - if the root finder returns a solution, we check that func(solution) = 0.0 +/ ROOT_ABS_TOLERANCE
    //       - if the solution is acceptable, we stop and return the solution
    //       - if the solution is not acceptable, we reduce the search step size and try again
    //       - if we reach the maximum number of search step reduction iterations, or the search step factor reduces to 1.0 (so search step size = 0.0),
    //         we stop and return a negative value for the root (indicating no root found)
   
    double guess      = 100.0;                                                                      // guess at 100 km s^-1 (arbitrary initial guess)

    double factorFrac = ADAPTIVE_RV_SEARCH_FACTOR_FRAC;                                             // search step size factor fractional part
    double factor     = 1.0 + factorFrac;                                                           // factor to determine search step size (size = guess * factor)
    
    std::pair<double, double> root(-1.0, -1.0);                                                     // initialise root - default return
    std::size_t tries = 0;                                                                          // number of tries
    bool done         = false;                                                                      // finished (found root or exceed maximum tries)?
    ERROR error       = ERROR::NONE;
    OStarRotationVelocityFunctor<double> func = OStarRotationVelocityFunctor<double>(desiredCDF);
    while (!done) {                                                                                 // while no error and acceptable root found

        bool isRising = true;                                                                       // guess for direction of search; CDF increases monotonically

        // run the root finder
        // regardless of any exceptions or errors, display any problems as a warning, then
        // check if the root returned is within tolerance - so even if the root finder
        // bumped up against the maximum iterations, or couldn't bracket the root, use
        // whatever value it ended with and check if it's good enough for us - not finding
        // an acceptable root should be the exception rather than the rule, so this strategy
        // shouldn't cause undue performance issues.
        try {
            error = ERROR::NONE;
            root  = boost::math::tools::bracket_and_solve_root(func, guess, factor, isRising, utils::BracketTolerance, it); // find root
            // root finder returned without raising an exception
            if (error != ERROR::NONE) { SHOW_WARN(error); }                                         // root finder encountered an error
            else if (it >= maxit) { SHOW_WARN(ERROR::TOO_MANY_RV_ITERATIONS); }                     // too many root finder iterations
        }
        catch(std::exception& e) {                                                                  // catch generic boost root finding error
            // root finder exception
            // could be too many iterations, or unable to bracket root - it may not
            // be a hard error - so no matter what the reason is that we are here,
            // we'll just emit a warning and keep trying
            if (it >= maxit) { SHOW_WARN(ERROR::TOO_MANY_RV_ITERATIONS); }                          // too many root finder iterations
            else             { SHOW_WARN(ERROR::ROOT_FINDER_FAILED, e.what()); }                    // some other problem - show it as a warning
        }

        // we have a solution from the root finder - it may not be an acceptable solution
        // so we check if it is within our preferred tolerance
        if (std::fabs(func(root.first + (root.second - root.first) / 2.0)) <= ROOT_ABS_TOLERANCE) { // solution within tolerance?
            done = true;                                                                            // yes - we're done
        }
        else if (std::fabs(func(root.first)) <= ROOT_ABS_TOLERANCE) {                               // solution within tolerance at endpoint 1?
            root.second=root.first;
            done = true;                                                                            // yes - we're done
        }
        else if (std::fabs(func(root.second)) <= ROOT_ABS_TOLERANCE) {                              // solution within tolerance at endpoint 2?
            root.first=root.second;
            done = true;                                                                            // yes - we're done
        }
        else {                                                                                      // no - try again
            // we don't have an acceptable solution - reduce search step size and try again
            factorFrac /= 2.0;                                                                      // reduce fractional part of factor
            factor      = 1.0 + factorFrac;                                                         // new search step size
            tries++;                                                                                // increment number of tries
            if (tries > ADAPTIVE_RV_MAX_TRIES || std::fabs(factor - 1.0) <= ROOT_ABS_TOLERANCE) {   // too many tries, or step size 0.0?
                // we've tried as much as we can - fail here with -ve return value
                root.first  = -1.0;                                                                 // yes - set error return
                root.second = -1.0;
                SHOW_WARN(ERROR::TOO_MANY_RV_TRIES);                                                // show warning
                done = true;                                                                        // we're done
            }
        }
    }
    
    // midway between brackets is our result
    // if necessary we could return the result as an interval here
    return root.first + (root.second - root.first) / 2.0;
}


/*
 * CalculateRotationalVelocityAtZAMS
 *
 * @brief
 * Calculate the ZAMS equatorial rotational velocity of a star.
 *
 * How the rotational velocity is calculated depends on the distribution specified by the
 * `--rotational-velocity-distribution` program option (passed as a paremeter to this function).
 * 
 * Possible values of the rotational velocity distribution are:
 * 
 *    - ROTATIONAL_VELOCITY_DISTRIBUTION::HURLEY
 *         - rotational velocity calculated using Hurley et al. 2000, eq 107
 * 
 *    - ROTATIONAL_VELOCITY_DISTRIBUTION::VLTFLAMES
 *         - rotational velocity calculated using:
 * 
 *            - for O-stars:
 *                 - single stars: Ramirez-Agudelo et al. 2013 (https://arxiv.org/abs/1309.2929)
 *                 - spectroscopic binaries: Ramirez-Agudelo et al. 2015 (https://arxiv.org/abs/1507.02286)
 * 
 *            - for B-stars: Dufton et al. 2013 (https://arxiv.org/abs/1212.2424)
 * 
 *            - for lower mass stars: Hurley et al. 2000, eq 107 (for now)
 *
 *    - ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO
 *         - no rotation; rotational velocity is 0.0 
 *
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function relies on the value of the ZAMS mass of the star, and should not be used
 * if the ZAMS mass is not known.
 * 
 *
 * double CalculateRotationalVelocityAtZAMS(double p_MZAMS) const
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS equatorial rotational velocity (km s^-1)
 */
COMPAS_PURE double BaseStar::CalculateRotationalVelocityAtZAMS(const double p_MZAMS) const {

    double vRot = 0.0;

    switch (OPTIONS->RotationalVelocityDistribution()) {                // which rotational velocity distribution?

        case ROTATIONAL_VELOCITY_DISTRIBUTION::HURLEY:                  // HURLEY

            // Hurley et al. 2000, eq 107 (uses fit from Lang 1992)
            vRot = (330.0 * PPOW(p_MZAMS, 3.3)) / (15.0 + PPOW(p_MZAMS, 3.45));
            break;

         case ROTATIONAL_VELOCITY_DISTRIBUTION::VLTFLAMES:              // VLTFLAMES

            // Rotational velocity based on VLT-FLAMES survey.
            // For O-stars (taken to be above 16 Msol), use results
            // of Ramirez-Agudelo et al. (2013) https://arxiv.org/abs/1309.2929 (single stars)
            // and Ramirez-Agudelo et al. (2015) https://arxiv.org/abs/1507.02286 (spectroscopic binaries)
            // For B-stars (taken to be between 2 and 16 Msol) use results
            // of Dufton et al. (2013) https://arxiv.org/abs/1212.2424
            // For lower mass stars, default back to  Hurley et al. 2000 distribution for now

            if (utils::Compare(p_MZAMS, 16.0) >= 0) {
                vRot = CalculateRotationalVelocityOStar_Ramirez2013();
                vRot = std::max(vRot, 0.0);                             // Set to no rotation if no positive solution found; warning already raised
            }
            else if (utils::Compare(p_MZAMS, 2.0) >= 0) {
                vRot = utils::InverseSampleFromTabulatedCDF(RAND->Random(), BStarRotationalVelocityCDFTable);
            }
            else {
                // Don't know what better to use for low mass stars so for now
                // default to Hurley et al. 2000, eq 107 (uses fit from Lang 1992)
                vRot = (330.0 * PPOW(p_MZAMS, 3.3)) / (15.0 + PPOW(p_MZAMS, 3.45));
            }
            break;

        case ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO:                    // ZERO
            vRot = 0.0;
            break;

        default:                                                        // unknown prescription
            // the only way this can happen is if someone added a ROTATIONAL_VELOCITY_DISTRIBUTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_VROT_PRESCRIPTION);              // throw error
    }

    return vRot;
}


/*
 * CalculateAngularFrequencyAtZAMS
 *
 * @brief
 * Calculate the angular frequency of a star at ZAMS.
 * 
 * Calls relevant angular frequency function based on the evolutionary mode given in program options.
 *
 *
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function relies on the values of the ZAMS mass and ZAMS radius of the star,
 * and should not be used if the ZAMS mass or the ZAMS radius is not known.
 * 
 * 
 * double CalculateAngularFrequencyAtZAMS(const double p_MZAMS, const double p_RZAMS) const
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @return                                      ZAMS angular frequency (rad yr^-1)
 */
 COMPAS_PURE double CalculateAngularFrequencyAtZAMS(const double p_MZAMS, const double p_RZAMS) const {

    double omega;

    Switch (OPTIONS->Mode()) {                                              // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                         // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                         // HURLEY BSE
            omega = CalculateAngularFrequencyAtZAMS_Hurley2000(p_MZAMS, p_RZAMS);
            break;
        
        default:                                                            // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                     // throw error
    }       

    return omega;
}


/*
 * CalculateAngularFrequencyAtZAMS_Hurley2000
 *
 * @brief
 * Calculate the ZAMS angular frequency of a star, per Hurley et al. 2000, eq 108
 *
 *
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function relies on the values of the ZAMS mass and ZAMS radius of the star,
 * and should not be used if the ZAMS mass or the ZAMS radius is not known.
 * 
 * 
 * double CalculateAngularFrequencyAtZAMS_Hurley2000(const double p_MZAMS, const double p_RZAMS)
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_RZAMS                         ZAMS radius of the star (Rsol)
 * @param       p_vRotDistribution              The rotational velocity distribution to use
 * @return                                      ZAMS angular frequency (rad yr^-1)
 */
COMPAS_PURE double BaseStar::CalculateAngularFrequencyAtZAMS_Hurley2000(const double p_MZAMS, const double p_RZAMS) {
    double vRot = CalculateRotationalVelocityAtZAMS(p_MZAMS);
    return 45.35 * vRot / p_RZAMS;
}


/*
 * CalculateAngularFrequencyCHE_Static
 *
 * @brief
 * Calculate the minimum angular frequency (in rad yr^-1) at which CHE will occur for a star,
 * given the ZAMS mass and the metallicity of the star
 *
 * Mandel's fit from Butler 2018 (see Riley et al. 2021, appendix A (https://doi.org/10.1093/mnras/stab1291))
 *
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function relies on the values of the ZAMS mass and ZAMS radius of the star,
 * and should not be used if the ZAMS mass or the ZAMS radius is not known.
 * 
 *
 * static double CalculateAngularFrequencyCHE_Static(const double p_MZAMS)
 *
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      Minimum angular frequency for CHE (rad yr^-1)
 */
COMPAS_PURE double BaseStar::CalculateAngularFrequencyCHE_Static(const double p_MZAMS) {

    constexpr double mRatio = p_MZAMS; // in MSol, so ratio is just p_MZAMS

    // calculate omegaCHE(M, Z = 0.004)
    double omegaZ004 = 0.0;
    if (p_MZAMS <= MANDEL_BUTLER_CHE_MASS_BREAK) {
        for (std::size_t i = 0; i < CHE_Coefficients.size(); i++) {
            omegaZ004 += CHE_Coefficients[i] * utils::intPow(mRatio, i) / PPOW(mRatio, 0.4);
        }
    }
    else {
        for (std::size_t i = 0; i < CHE_Coefficients.size(); i++) {
            omegaZ004 += CHE_Coefficients[i] * utils::intPow(MANDEL_BUTLER_CHE_MASS_BREAK, i) / PPOW(mRatio, 0.4);
        }
    }

    // calculate omegaCHE(M, Z)
    return (1.0 / ((0.09 * log(GLOBALS->Metallicity() / 0.004)) + 1.0) * omegaZ004) * SECONDS_IN_YEAR;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    SUPERNOVAE                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateSNkickMagnitude
 *
 * @brief
 * Calculate supernova kick magnitude, based on the user-specified kick magnitude distribution
 * (via option `--kick-magnitude-distribution`), and the current supernova event type.
 *
 *
 * double BaseStar::CalculateSNkickMagnitude(
 *     const STELLAR_TYPE       p_RemnantType,
 *     const double             p_Mass,
 *     const double             p_EjectaMass,
 *     const double             p_RemnantMass,
 *     const StellarSNDetailsT& p_SNdetails
 * ) const
 *
 * @param       p_RemnantType                   Expected stellar type of the remnant (must be NS or BH)
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_EjectaMass                    Change in mass of the exploding star (i.e. mass of the ejecta) (Msol)
 * @param       p_RemnantMass                   The mass of the remnant (Msol)
 * @param       p_SNdetails                     Supernova datails object
 * @return                                      Kick magnitude (km s^-1)
 */
COMPAS_PURE StellarSNDetailsT BaseStar::CalculateSNkickMagnitude(
    const STELLAR_TYPE       p_RemnantType,
    const double             p_Mass,
    const double             p_EjectaMass,
    const double             p_RemnantMass,
    const StellarSNDetailsT& p_SNdetails
) const {

    if (!utils::IsOneOf(p_RemnantType, {STELLAR_TYPE::NS, STELLAR_TYPE::BLACK_HOLE})) { // NS or BS expected?
        THROW_ERROR(ERROR::UNEXPECTED_REMNANT_TYPE);                                    // no - throw error - this is a code defect
    }

    StellarSNDetailsT SNdetails;                                                        // to be returned

	double kick = 0.0;                                                                  // default kick magnitude
    ERROR error = ERROR::NONE;                                                          // no error presumed

    if (!p_SNdetails.suppliedKickParamaters.magnitude.has_value() ||                    // user did not supply kick magnitude, or
         p_SNdetails.suppliedKickParamaters.magnitudeRandom.has_value()) {              // ... wants to draw magnitude using supplied random number

        double sigma;
        switch (p_SNdetails.currentEvent.p_SNevent) {                                   // what type of supernova event happening now?

		    case SN_EVENT::AIC:                                                         // no kick for AIC 
		    case SN_EVENT::SNIA:                                                        // no kick for SNIA 
		    case SN_EVENT::HeSD:                                                        // no kick for HeSD 
			    sigma = 0.0;
                break;

		    case SN_EVENT::PPISN:                                                       // PPISN

                if (!OPTIONS->NatalKickForPPISN()) {                                    // user specified natal kicks for PPISN?
                    sigma = 0.0;                                                        // no - no natal kicks for PPISN
                    break;                                                              // exit switch for PPISN and no kick
                }                                                                       // otherwise fall thorugh

            case SN_EVENT::CCSN:                                                        // CCSN (PPISN with kicks falls through)                                    

                // draw a random kick magnitude from the user-specified distribution
                // based on the expected stellar type of the remnant (NS or BH)

                if (p_RemnantType == STELLAR_TYPE::NEUTRON_STAR)                        // neutron star?
                    sigma = OPTIONS->KickMagnitudeDistributionSigmaCCSN_NS();           // yes
                else if (p_RemnantType == STELLAR_TYPE::BLACK_HOLE)                     // no - black hole?
                    sigma = OPTIONS->KickMagnitudeDistributionSigmaCCSN_BH();           // yes
                else                                                                    // unexpected stellar type - shouldn't happen                                        
                    error = ERROR::UNEXPECTED_STELLAR_TYPE;                             // set error value

                break;

		    case SN_EVENT::ECSN:                                                        // ECSN may have a separate kick prescription
			    sigma = OPTIONS->KickMagnitudeDistributionSigmaForECSN();
                break;

		    case SN_EVENT::USSN:                                                        // USSN may have a separate kick prescription
			    sigma = OPTIONS->KickMagnitudeDistributionSigmaForUSSN();
                break;

            default:                                                                    // unexpected SN event type
                // the only way this can happen is if the SN_EVENT passed to this function is not accounted
                // for in this code.  We should not default here, with or without a warning.
                // We are here because the code passed an event type that this function doesn't account for,
                // and that should be flagged as an error and result in termination of the evolution of the
                // star or binary.
                // The correct fix for this is to add code to this function for the missing event type,
                // or fix the calling code to pass an event type that is handled by this function.

                error = ERROR::UNEXPECTED_SN_EVENT;                                     // set error value
	    }
    
	    if (error == ERROR::NONE) {                                                     // check for errors
                                                                                        // no errors - draw kick magnitude
            kick = DrawSNkickMagnitude(p_SNdetails.COCoreMassAtCOFormation, p_EjectaMass, p_RemnantMass, sigma, p_KickRandom, p_SNdetails.currentEvent.p_SNevent);
        }
    }
    else {                                                                              // user supplied kick parameters, and wants to use supplied kick magnitude, so ...
        kick = p_KickMagnitude;                                                         // ... use it 
    }

	if (error == ERROR::NONE) {                                                         // check for errors
                                                                                        // no errors
        SNdetails = p_SNdetails;                                                        // copy given SN details

        SNdetails.currentKickParameters.magnitudeDrawn = kick;                          // set drawn kick magnitude

        // if there was a core-collapse supernova event this timestep, and the expected remnant is a black hole,
        // re-weight the (drawn) kick by the mass of the remnant according to the user-specified black hole kicks
        // option, if relevant (default is no reweighting)
        if (p_SNevent == SN_EVENT::CCSN && p_RemnantType == STELLAR_TYPE::BLACK_HOLE) {
            kick *= BH::CalculateSNkickWeighting_Static(p_Mass, p_SNdetails.fallbackFraction);
        }
        else {                                                                          // otherwise
            SNdetails.fallbackFraction = 0.0;                                           // set fallback fraction to zero
        }

        SNdetails.currentKickParameters.magnitude = kick;                               // set updated kick magnitude
    }
    else {                                                                              // error occurred
        THROW_ERROR(error);                                                             // throw error
    }

    // return SN details object by value - NRVO takes care of performance/efficiency
    return SNdetails;
}


/*
 * DrawSNkickMagnitude
 *
 * @brief
 * Draw a kick magnitude.  Calls appropriate functions based on the user-specified
 * kick magnitude distribution (via option `--kick-magnitude-distribution`)
 *
 *
 * double DrawSNkickMagnitude(
 *     const double   p_COCoreMass,
 *     const double   p_EjectaMass,
 *     const double   p_RemnantMass,
 *     const double   p_Sigma,
 *     const double   p_Rand,
 *     const SN_EVENT p_SNevent
 * ) const
 *
 * @param       p_COCoreMass                    CO core mass of exploding star (Msol)
 * @param       p_EjectaMass                    Change in mass of the exploding star (i.e. mass of the ejecta) (Msol)
 * @param       p_RemnantMass                   Mass of the remnant (Msol)
 * @param       p_Sigma                         Distribution scale parameter - affects the spread of the distribution
 * @param       p_Rand                          Random number for drawing from the distribution [0, 1]
 * @param       p_SNevent                       Current SN event type                      
 * @return                                      Drawn kick magnitude (km s^-1)
 */
COMPAS_PURE double BaseStar::DrawSNkickMagnitude(
    const double   p_COCoreMass,
    const double   p_EjectaMass,
    const double   p_RemnantMass,
    const double   p_Sigma,
    const double   p_Rand,
    const SN_EVENT p_SNevent
) const {

	double kick;

    switch (OPTIONS->KickMagnitudeDistribution()) {                                 // which kick magnitude distribution?

        case KICK_MAGNITUDE_DISTRIBUTION::BRAYELDRIDGE:                             // BRAY ELDRIDGE
            kick = DrawSNkickMagnitude_BrayEldridge2018(
                p_EjectaMass,
                p_RemnantMass,
                BRAY_ELDRIDGE_CONSTANT_VALUES.at(BRAY_ELDRIDGE_CONSTANT::ALPHA),
                BRAY_ELDRIDGE_CONSTANT_VALUES.at(BRAY_ELDRIDGE_CONSTANT::BETA)
            );
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::FIXED:                                    // FIXED
            kick = p_Sigma;
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::FLAT:                                     // FLAT
            kick = DrawSNkickMagnitude_Flat(OPTIONS->KickMagnitudeDistributionMaximum(), p_Rand);
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::LOGNORMAL:                                // LOGNORMAL
            // only draw Disberg & Mandel 2025 kicks for CCSN or PPISN (if they receive a kick)
            // use Maxwellians with the value of sigma set in CalculateSNkickMagnitude() for other SN types
            if (p_SNevent == SN_EVENT::CCSN || (p_SNevent == SN_EVENT::PPISN && OPTIONS->NatalKickForPPISN())) {
                // maximum kick of DISBERG_MANDEL_MAX_KICK = 1000 km/s, following Disberg & Mandel 2025
                // normalise the draw so that the kick is uniformly drawn from the inverse CDF below this maximum
                const double cdfMax = gsl_cdf_lognormal_P(DISBERG_MANDEL_MAX_KICK, DISBERG_MANDEL_MU, DISBERG_MANDEL_SIGMA);
                kick = gsl_cdf_lognormal_Pinv(p_Rand * cdfMax, DISBERG_MANDEL_MU, DISBERG_MANDEL_SIGMA);
            }
            else {
                kick = DrawSNkickMagnitude_Maxwellian(p_Sigma, p_Rand);
            }
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::MAXWELLIAN:                               // MAXWELLIAN
            kick = DrawSNkickMagnitude_Maxwellian(p_Sigma, p_Rand);
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::MULLER2016:                               // MULLER2016
            kick = DrawSNkickMagnitude_Muller2016(p_COCoreMass);
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::MULLER2016MAXWELLIAN: {                   // MULLER2016-MAXWELLIAN

            double mullerSigma = DrawSNkickMagnitude_Muller2016(p_COCoreMass) / std::sqrt(3.0);

            kick = DrawSNkickMagnitude_Maxwellian(mullerSigma, p_Rand);
            } break;

        case  KICK_MAGNITUDE_DISTRIBUTION::MULLERMANDEL:                            // MULLERMANDEL
            kick = DrawSNkickMagnitude_MullerMandel2020(p_COCoreMass, p_RemnantMass, p_Rand, p_SNevent);
            break;

        case KICK_MAGNITUDE_DISTRIBUTION::ZERO:                                     // ZERO
            kick = 0.0;                                                             // no kick
            break;
            
        default:                                                                    // unknown prescription
            // the only way this can happen is if someone added a KICK_MAGNITUDE_DISTRIBUTION and it
            // isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_KICK_MAGNITUDE_DISTRIBUTION);                // throw error
    }

    return kickMagnitude / OPTIONS->KickScalingFactor();
}


/*
 * DrawSNkickMagnitude_Muller2016
 *
 * @brief
 * Draw kick magnitude, per Muller et al. 2016,
 * as presented in eq. B5 of Vigna-Gomez et al. 2018 (see arXiv:1805.07974)
 * 
 * Follows Vigna-Gomez et al. 2018, using 35 km/s as the peak of a low-kick
 * Maxwellian (e.g. USSN, ECSN)
 * 
 * Black Holess do not get natal kicks.
 * 
 *
 * double DrawSNkickMagnitude_Muller2016(const double p_COCoreMass) const
 *
 * @param       p_COCoreMass                    CO core mass of the (exploding) star (Msol)
 * @return                                      Drawn kick magnitude (km s^-1)
 */
GNU_CONST double BaseStar::DrawSNkickMagnitude_Muller2016(const double p_COCoreMass) const {

    double kick;

	     if (p_COCoreMass < 1.372) kick = 0.0;
	else if (p_COCoreMass < 1.49 ) kick = 35.0  + (1000.0 * (p_COCoreMass - 1.372));    // following Vigna-Gomez et al. 2018
    else if (p_COCoreMass < 1.65 ) kick = 90.0  + ( 650.0 * (p_COCoreMass - 1.49));
	else if (p_COCoreMass < 2.4  ) kick = 100.0 + ( 175.0 * (p_COCoreMass - 1.65));
    else if (p_COCoreMass < 3.2  ) kick = 200.0 + ( 550.0 * (p_COCoreMass - 2.4));
    else if (p_COCoreMass < 3.6  ) kick = 80.0  + ( 120.0 * (p_COCoreMass - 3.2));
    else if (p_COCoreMass < 4.05 ) kick = 0.0;                                          // will be a Black Hole
    else if (p_COCoreMass < 4.6  ) kick = 350.0 + (  50.0 * (p_COCoreMass - 4.05));
    else if (p_COCoreMass < 5.7  ) kick = 0.0;                                          // will be a Black Hole
    else if (p_COCoreMass < 6.0  ) kick = 275.0 - ( 300.0 * (p_COCoreMass - 5.7));
    else                           kick = 0.0;                                          // will be a Black Hole

    return kick;
}


/*
 * DrawSNkickMagnitude_MullerMandel2020
 *
 * @brief
 * Draw kick magnitude, per Mandel & Mueller 2020
 * 
 *
 * double DrawSNkickMagnitude_MullerMandel2020(const double p_COCoreMass, const double p_RemnantMass, const double p_Rand) const
 * 
 * @param       p_COCoreMass                    CO core mass of the (exploding) star (Msol)
 * @param       p_RemnantMass                   Mass of the remnant (Msol)
 * @param       p_Rand                          Random number for drawing from the distribution [0, 1]
 * @param       p_SNevent                       Current SN event type
 * @return                                      Drawn kick magnitude (km s^-1)
 */
COMPAS_PURE double BaseStar::DrawSNkickMagnitude_MullerMandel2020(const double p_COCoreMass, const double p_RemnantMass, const double p_Rand, const SN_EVENT p_SNevent) const {
    					
	double kick;

    // Mandel & Mueller 2020 call for USSN kicks to be treated in the same way as CCSN kicks,
    // but if the override flag is set (via option `--USSN-kicks-override-mandel-muller`), 
    // we set the USSN kick to be equal to the user-provided magnitude
    if (OPTIONS->USSNKicksOverrideMandelMuller() && p_SNevent == SN_EVENT::USSN) {  // use user-supplied kick magnitude? // utils::SNEventType(m_SupernovaDetails.events.current) == SN_EVENT::USSN ) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        kick = OPTIONS->KickMagnitudeDistributionSigmaForUSSN();                    // yes
    }
    else {                                                                          // no - calculate kick magnitude
	    double muKick;
        double sigmaKick;

	    if (p_RemnantMass < OPTIONS->MaximumNeutronStarMass()) {
		    muKick    = std::max(OPTIONS->MullerMandelKickMultiplierNS() * (p_COCoreMass - p_RemnantMass) / p_RemnantMass, 0.0);
            sigmaKick = OPTIONS->MullerMandelSigmaKickNS();
	    }
	    else {
		    muKick    = std::max(OPTIONS->MullerMandelKickMultiplierBH() * (p_COCoreMass - p_RemnantMass) / p_RemnantMass, 0.0);
            sigmaKick = OPTIONS->MullerMandelSigmaKickBH();
	    }

        // quantile of -1 in the Gaussian CDF
        // the goal is to draw from the cut-off Gaussian since the kick must exceed 0
        const double quantile0 = gsl_cdf_gaussian_P(-1.0, sigmaKick);  
        const double rand      = quantile0 + p_Rand * (1.0 - quantile0);

        kick = muKick * (1.0 + gsl_cdf_gaussian_Pinv(rand, sigmaKick));
    }

	return kick;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    TEMPERATURE                                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////











//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//.............................................................................................................................................................
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//.............................................................................................................................................................
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//.............................................................................................................................................................
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//.............................................................................................................................................................
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//.............................................................................................................................................................
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//.............................................................................................................................................................
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//.............................................................................................................................................................
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~











///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    MISCELLANEOUS FUNCTIONS / CONTROL FUNCTIONS                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Calculate next timestep for stellar evolution
 *
 * Timestep based on stellar type, age, etc.
 *
 *
 * double CalculateTimestep()
 *
 * @return                                      Timestep
 */
double BaseStar::CalculateTimestep() {
    
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::CalculateTimestep(@start), m_StellarType = " << (int)m_StellarType << ", m_StellarTypePrev = " << (int)m_StellarTypePrev << ", m_Radius = " << m_Radius << ", m_RadiusPrev = " << m_RadiusPrev << ", m_dtPrev = " << m_dtPrev << "\n";  
    double radialExpansionTimescale = CalculateRadialExpansionTimescale();
    double massChangeTimescale      = CalculateMassChangeTimescale();
    double dt                       = 0.0;
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::CalculateTimestep(@0), massChangeTimescale = " << massChangeTimescale << ", radialExpansionTimescale = " << radialExpansionTimescale << "\n";  

    if (massChangeTimescale > 0.0)                                                                           // non-positive means it could not be computed (e.g., just after stellar type change)
        dt = OPTIONS->MassChangeFraction() * massChangeTimescale;
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::CalculateTimestep(@1), dt = " << dt << "\n";  

    if (radialExpansionTimescale > 0.0)                                                                      // non-positive means it could not be computed (e.g., just after stellar type change)
        dt = dt <= 0.0 ? OPTIONS->RadialChangeFraction() * radialExpansionTimescale : std::min(dt, OPTIONS->RadialChangeFraction() * radialExpansionTimescale);
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::CalculateTimestep(@2), dt = " << dt << "\n";  

    // the GBparams and Timescale calculations need to be done
    // before the timestep calculation - since the binary code
    // calls this functiom, the GBparams and Timescale functions
    // are called here
    CalculateGBparams();                                                                                    // calculate giant branch parameters
    CalculateTimescales();                                                                                  // calculate timescales

    dt = dt <= 0.0 ? ChooseTimestep(m_Age) : std::min(dt, ChooseTimestep(m_Age));
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::CalculateTimestep(@3), dt = " << dt << "\n";  

    // there is a chance that mass loss from winds is much faster than previously estimated if, say, LBV winds have turned on
    // we therefore precompute the mass loss rate to avoid taking an overly long timestep, despite the extra computational costs
    double massChangeWinds = m_Mass - CalculateMassLossValues(dt, false);
    if (utils::Compare(massChangeWinds, 0.0) != 0) dt = std::min(dt, OPTIONS->MassChangeFraction() * (dt * m_Mass / std::fabs(massChangeWinds)));
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::CalculateTimestep(@4), dt = " << dt << "\n";  

    dt = std::max(QUANTISE_DT(dt), NUCLEAR_MINIMUM_TIMESTEP);                                               // quantised; not less than nuclear minimum

if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << std::boolalpha << std::setprecision(15) << "BaseStar::CalculateTimestep(), returning dt = " << dt << "\n";    
    return dt;
}











/*
 * EvolveOneTimestep
 *
 * @brief
 * Apply mass changes (dM, dM0) if required, and evolve the star one timestep (dt)  <<<< mass changes for when called from BSE only
 * 
 * Evolving the star one timestep advances the simulation time and the star's age by dt, and recalculates the
 * attributes of the star (Age, Radius, Luminosity etc.) given the mass changes and dt.  If as a result of the
 * changes to the attributes of the star it should change stellar type, the new stellar type is returned as the
 * functional return (if no change to stellar type is required the star's current stellar type is returned). 
 * 
 * No change to stellar type is effected here, or as a result of the call to this function.  The caller of this
 * function is expected to check the stellar type returned and manage any stellar type switch required.
 *
 * Free parameters in the update process are the star's mass (m_Mass), initial mass (m_Mass0), the star's age
 * (m_Age) and the simulation time attribute (m_Time):
 *
 *    - if required, the star's mass is changed by the amount passed as the p_DeltaMass parameter before other
 *      attributes are updated.  The p_DeltaMass parameter may be zero, in which case no change is made to the
 *      star's mass before the attributes of the star are calculated.
 *
 *    - if required, the star's effective initial mass is changed by the amount passed as the p_DeltaMass0
 *      parameter before other attributes are updated.  The p_DeltaMass0 parameter may be zero, in which case
 *      no change is made to the star's effective initial mass before the attributes of the star are calculated.
 *      Mass0 in Hurley et al. 2000 is overloaded by the introduction of mass loss (see section 7.1).
 *
 *    - if required, the star is aged by the amount passed as the p_DeltaTime parameter, and the simulation time
 *      is advanced by the same amount, before other attributes are updated.  The p_dt parameter may be zero, in 
 *      which case no change is made to the star's age or the simulation time attribute.
 *
 *
 * Before updating attributes we check whether the star:
 *    - is due to become a supernova - if so we return with no change (supernovae are handled elsewhere)
 *    - is (or should be) a massless remnant - if so we return twith no change
 *    - should skip this phase for this timestep (checked after applying p_DeltaMass, p_DeltaMass0 and p_DeltaTime)
 *
 * If none of the above are true the star evolves on phase for the specified timestep (which may be 0, in which case
 * the star's attributes other than age are re-calculated), then the need to evolve the star off phase is checked.
 *
 * If p_DeltaMass, p_DeltaMass0 and p_DeltaTime are all passed as zero the checks for massless remnant and supernova
 * are performed, but no other changes to the star's attributes are made - unless the p_ForceRecalculate parameter is
 * set true.
 *
 *
 * STELLAR_TYPE EvolveOneTimestep(const double p_DeltaMass, const double p_DeltaMass0, const double p_DeltaTime, const bool p_ForceRecalculate)
 *
 * @param           p_dt                        The timestep to take to change state in (Myr)
 * @param   [IN]    p_dM                        The change in mass to apply in Msol     <<<< mass changes for when called from BSE only
 * @param   [IN]    p_dM0                       The change in mass0 to apply in Msol    <<<< mass changes for when called from BSE only
 * @param   [IN]    p_ForceRecalculate          Specifies whether the star's attributes should be recalculated even if the three deltas are 0.0  REVISIT THIS<<<<<<<<<<<<<<<<<<
 *                                              (optional, default = false)
 * @return                                      Stellar type to which star should evolve
 */
STELLAR_TYPE BaseStar::EvolveOneTimestep(const double p_dM, const double p_dM0, const double p_dt, const bool p_ForceRecalculate) {
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOneTimestep(@entry), p_dM = " << p_dM << ", p_dM0 = " << p_dM0 << ", p_dt = " << p_dt << ", p_ForceRecalculate = " << p_ForceRecalculate << ", m_Mass = " << m_Mass << ", m_Radius = " << m_Radius << "\n";

    STELLAR_TYPE nextStellarType = m_StellarType;                                               // next stellar type - defaults to current

    // THIS IS NEW HERE - WAS IN Star.cpp <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    nextStellarType = m_Star->ResolveMassLoss(p_dt);                                                            // apply wind mass loss if required


    if (ShouldBeMasslessRemnant()) return STELLAR_TYPE::MASSLESS_REMNANT;                       // do not update the star if it has lost all of its mass

    if (IsSupernova()) return m_StellarType;                                                    // do nothing if supernova is pending (handled elsewhere)

    // advance the simulation of this star by one timestep
    bool recalc = p_ForceRecalculate;                                                           // need to recalculate attribute values?
    
    // update mass as required (only change if delta != 0) and prevent -ve
    if (utils::Compare(p_dM,  0.0) != 0) { 
        m_Mass = std::max(0.0, m_Mass + p_dM);
        recalc = true;
    }
    
    // update mass0 as required (only change if delta != 0) and prevent -ve
    if (utils::Compare(p_dM0, 0.0) != 0) {
        m_Mass0 = std::max(0.0, m_Mass0 + p_dM0);
        recalc  = true;
    }
    
    // GBparams and Timescale calculations need to be done before taking the timestep
    UpdateGBparams();
    UpdateTimescales();

    // record some current values before they are (possibly) changed by evolution
    // since these will be previous timestep values we only record them if dt > 0
    // (i.e. we are actually taking a timestep)
    if (p_dt > 0.0) {                                                                           // don't use utils::Compare() here
        m_StellarTypePrev = m_StellarType;
        m_MassPrev        = m_Mass;
        m_RadiusPrev      = m_Radius;
        recalc            = true;
    }

    // update attributes if necessary
    if (recalc) {                                                                               // need to update attributes?
                                                                                                // yes
        SetDt(p_dt);                                                                            // set timestep
        UpdateEffectiveZAMSLandR();                                                             // update effective ZAMS luminosity and radius if necessary

        // evolve the star one timestep
        if (ShouldSkipPhase()) nextStellarType = ResolveSkippedPhase();                         // skip phase if required - per stellar type
        else {                                                                                  // phase not skipped
            nextStellarType = EvolveOnPhase(m_dt);                                              // evolve on phase
            if (nextStellarType == m_StellarType) {                                             // need to switch to new stellar type?
                nextStellarType = ResolveEndOfPhase();                                          // no - check for need to move off phase
            }   
        }
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOneTimestep(@1), m_Mass = " << m_Mass << ", m_Radius = " << m_Radius << "\n";
    }

    return nextStellarType;                                                                     // stellar type to which star should evolve
}












/*
 * Evolve the star on its current phase - take one timestep on the current phase
 *
 * If as a result of the evolution of the star it should change stellar type, the new stellar type is returned
 * as the functional return (if no change to stellar type is required the star's current stellar type is returned). 
 * 
 * No change to stellar type is effected here, or as a result of the call to this function.  The caller of this
 * function is expected to check the stellar type returned and manage any stellar type switch required.
 *
 * STELLAR_TYPE EvolveOnPhase(const double p_dt)
 *
 * @param   [IN]    p_dt                        Timestep in Myr
 * @return                                      Stellar Type to which star should evolve
 */
STELLAR_TYPE BaseStar::EvolveOnPhase(const double p_dt) {
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@entry), p_dt = " << p_dt << ", m_Age= " << m_Age << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << "\n";
    
    STELLAR_TYPE nextStellarType = m_StellarType;                                       // next stellar type - defaults to current

    if (ShouldEvolveOnPhase()) {                                                        // should evolve timestep on phase?
                                                                                        // yes
        UpdateMainSequenceCoreMass(p_dt, -m_Mdot);                                      // update core mass, relevant for MS stars

if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@1), p_dt = " << p_dt << ", m_Age= " << m_Age << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";
        m_Tau        = CalculateTauOnPhase();
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@2), p_dt = " << p_dt << ", m_Age= " << m_Age << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";

        m_COCoreMass = CalculateCOCoreMassOnPhase();
        m_CoreMass   = CalculateCoreMassOnPhase();
        m_HeCoreMass = CalculateHeCoreMassOnPhase();

        m_Luminosity = CalculateLuminosityOnPhase();

        // Calculate abundances
        m_HydrogenAbundanceCore    = CalculateHAbundanceCore();
        m_HydrogenAbundanceSurface = CalculateHAbundanceSurface();  
        m_HeliumAbundanceCore      = CalculateHeAbundanceCore();
        m_HeliumAbundanceSurface   = CalculateHeAbundanceSurface();
       
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@3), p_dt = " << p_dt << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";
        std::tie(m_Radius, nextStellarType) = CalculateRadiusAndStellarTypeOnPhase();   // radius and possibly new stellar type
if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@4), p_dt = " << p_dt << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";

        m_Mu = CalculateHurleyPerturbationMu();

        PerturbLuminosityAndRadiusOnPhase();

        m_Temperature = CalculateTemperatureOnPhase();

        if (p_dt > 0.0) {
            STELLAR_TYPE thisStellarType = ResolveEnvelopeLoss();                       // resolve envelope loss if it occurs - possibly new stellar type
            if (thisStellarType != m_StellarType) {                                     // thisStellarType overrides stellarType (from CalculateRadiusAndStellarTypeOnPhase())
                nextStellarType = thisStellarType;
            }
            AdvanceAgeAndTime(p_dt);                                                    // advance age of star and simulation time
        }
    }

    return nextStellarType;
}
















// SHOULD THIS BE IN STAR class?????????????????????????????????????????????????????????????????  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// p_dt = 0.0 means do nothing
StarState BaseStar::AdvanceOneTimestep(const double p_dt, const StarState& p_State) {
    
    if (p_dt <= 0.0) return p_State;                            // do nothing if timestep <= 0.  State is returned unchanged.


    StarState interimState = p_State;                           // copy given state to interim (local) state

    StarState currentState = m_StateHistory.CurrentState;       // current state

    // unpack variables from the current state - so we call the getters only once
    
    const double age        = Age();
    const double luminosity = Luminosity();
    const double mass       = Mass();
    const double radius     = Radius();
    const double tau        = Tau();
    const double time       = Time();

    

    // these are just sanity checks - massless remnant and SNe should have
    // been handled outside this function, but if we're here then just return
////////    if (ShouldBeMasslessRemnant()) return STELLAR_TYPE::MASSLESS_REMNANT;                       // do not update the star if it has lost all of its mass
////////    if (IsSupernova()) return m_StellarType;                                                    // do nothing if supernova is pending (handled elsewhere)    



    // (Re)calculate attributes of star and update interim state object
    // (we do this to allow for diagnostic printing of interim state)
    //
    // Since all (re)calculations of state variables (stellar attributes) in a
    // timestep are performed using the *current* state values, the order of the
    // (re)calculations of state variables is not significant - interim stellar
    // attributes (local state variables) should not be used to (re)calculate any
    // state variables.
    //
    // We calculate delta values here, and only apply those deltas to the state
    // at the end of this timestep (before starting the next timestep).  Ideally
    // we would calculate the rate-of-change of each attribute using derivatives,
    // but for now we settle on deltas (we don't bother changing the deltas to
    // rate-of-change over the timestep becuase we would just multiply by the
    // timestep at the end of the timestep anyway, so just waste CPU cycles.)


    // age, tau, time

    double deltaAge = p_dt;                                                     // age delta is just the size of the timestep
    interimState.SetAge(age + deltaAge);                                        // update interim state

    double deltaTau = CalculateTau() - tau;                                     // <<<<<<<<<< DONE >>>>>>>>>>  phase-relative age delta is just the size of the timestep
    interimState.SetTau(tau + deltaTau);                                        // update interim state

    double deltaTime = p_dt;                                                    // simulated evolution time delta is just the size of the timestep
    interimState.SetTime(time + deltaTime);                                     // update interim state


    // masses

    // effective initial mass
    //
    // The effective initial mass of the star is the variable M0 in Hurley et al. 2000.  M0 is introduced in
    // section 7 of Hurley et al. 2000 and is described there as the "initial mass". M0 is the phase-specific
    // initial mass - i.e., it is the initial mass for the current evolutionary phase of the star, but can 
    // track Mt (mass at time t), depending upon the phase, so (as stated in Hurley et al. 2000), it is really
    // the "effective initial mass".
    //
    // The effective initial mass in Hurley et al. 2000 is relavant to MS, HG, and HeMS stars only.
    //
    // Since M0, or the effective initial mass, is a Hurley construct, we only calculate it here if the evolution
    // mode is a Hurley evolution mode, and we only provide the Hurley method for calculating it.  If the addition
    // of other evolution modes require a similar "effective intial mass" construct, we can add functionality later
    // to calculate M0 as required for different modes of evolution.
    //  
    if (OPTIONS->Mode() == EVOLUTION_MODE::SSE_HURLEY || OPTIONS->Mode() == EVOLUTION_MODE::BSE_HURLEY) {   // Hurley evolution mode?
                                                                                                            // yes
        interimState.SetEffectiveInitialMass(CalculateEffectiveInitialMass_Hurley2000());                   // update interim state  <<<<<<<<<< DONE >>>>>>>>>> 
    }

    // calculate mass loss for dt, clamped to [0.0, photon tiring limit]
    //
    // the photon tiring limit is (effectively) infinity unless the user specified
    // the `--check-photon-tiring-limit` option

    double dMdt;
    MASS_LOSS_TYPE dominantMLtype;
    std::tie(dMdt, dominantMLtype) = CalculateMassLossRate();                   // <<<<<<<<<< DONE >>>>>>>>>>  Msol yr^1, [0, 1]
    
    dMdt = dt > 0.0 ? dMdt * 1.0E6 : 0.0;                                       // dMdt (Msol Myr^-1, [0, 1])

    const double upperBound = OPTIONS->CheckPhotonTiringLimit() ? luminosity / (G_SOLAR_YEAR * mass / radius) : DBL_MAX;
    const double deltaMass  = -std::min(std::max(0.0, dMdt * dt), upperBound);  // mass loss (Msol) (-ve for mass loss)

    // calculate masses

    interimState.SetMass(mass + deltaMass);                                     // update interim state


    UpdateMainSequenceCoreMass(p_dt, -m_Mdot);                                      // update core mass, relevant for MS stars
    double coreMass   = CalculateCoreMass();  /// <<<< BRCEK IN HERE  MSCoreMass??????


    double COcoreMass = CalculateCOCoreMass();                                  // <<<<<<<<<< DONE >>>>>>>>>>

    double HecoreMass = CalculateHeCoreMass();                                  // <<<<<<<<<< DONE >>>>>>>>>>


////////////////////////// FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    // calculate new radius and check need for stellar type change
    // only He stars (specifically HeHG and HeGB) might need to change stellar type
    // in response to new radius
    double       radius;
    STELLAR_TYPE stellarType;
    std::tie(radius, stellarType) = CalculateRadiusOnPhase();

    

    double luminosity = CalculateLuminosityOnPhase();


    

    
    
            // Calculate abundances
            m_HydrogenAbundanceCore    = CalculateHAbundanceCore();         /// <<< see BRCEK
            m_HydrogenAbundanceSurface = CalculateHAbundanceSurface();          /// <<< see BRCEK
            m_HeliumAbundanceCore      = CalculateHeAbundanceCore();         /// <<< see BRCEK
            m_HeliumAbundanceSurface   = CalculateHeAbundanceSurface();         /// <<< see BRCEK

    
            UpdateAgeAfterMassLoss(); // update age (MS, HG & HeMS)

            ApplyMassTransferRejuvenationFactor(); // apply age rejuvenation factor

            SetAngularMomentum(m_AngularMomentum + angularMomentumChange);

            


    // update state

    // copy some current values to prev values before they are (possibly) changed
    // since these will be previous timestep values we only record them if dt > 0
    // (i.e. we are actually taking a timestep)
    if (dt > 0.0) {                                                                             // don't use utils::Compare() for thresholds
        m_StellarTypePrev = m_StellarType;
        m_MassPrev        = m_Mass;
        m_RadiusPrev      = m_Radius;
    }




        // update mass as required (only change if delta != 0) and prevent -ve
        if (utils::Compare(p_dM,  0.0) != 0) { 
            m_Mass = std::max(0.0, m_Mass + p_dM);
            recalc = true;
        }
        
        // update mass0 as required (only change if delta != 0) and prevent -ve
        if (utils::Compare(p_dM0, 0.0) != 0) {
            m_Mass0 = std::max(0.0, m_Mass0 + p_dM0);
            recalc  = true;
        }


        
        
    
        // update attributes if necessary
        if (recalc) {                                                                               // need to update attributes?
                                                                                                    // yes
            SetDt(p_dt);                                                                            // set timestep
            UpdateEffectiveZAMSLandR();                                                             // update effective ZAMS luminosity and radius if necessary
    
            // evolve the star one timestep
            if (ShouldSkipPhase()) nextStellarType = ResolveSkippedPhase();                         // skip phase if required - per stellar type
            else {                                                                                  // phase not skipped
                nextStellarType = EvolveOnPhase(m_dt);                                              // evolve on phase
                if (nextStellarType == m_StellarType) {                                             // need to switch to new stellar type?
                    nextStellarType = ResolveEndOfPhase();                                          // no - check for need to move off phase
                }   
            }
    if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOneTimestep(@1), m_Mass = " << m_Mass << ", m_Radius = " << m_Radius << "\n";
        }

        
    

    // timescales and giantbranch parameters (as necessary)

    interimState.SetTimescales(CalculateTimescales());                          // <<<<<<<<<< DONE >>>>>>>>>>
    interimState.SetGBparams(CalculateGBparams());                              // <<<<<<<<<< DONE >>>>>>>>>>

        return nextStellarType;                                                                     // stellar type to which star should evolve
    }
    
    














    /*
     * Evolve the star on its current phase - take one timestep on the current phase
     *
     * If as a result of the evolution of the star it should change stellar type, the new stellar type is returned
     * as the functional return (if no change to stellar type is required the star's current stellar type is returned). 
     * 
     * No change to stellar type is effected here, or as a result of the call to this function.  The caller of this
     * function is expected to check the stellar type returned and manage any stellar type switch required.
     *
     * STELLAR_TYPE EvolveOnPhase(const double p_dt)
     *
     * @param   [IN]    p_dt                        Timestep in Myr
     * @return                                      Stellar Type to which star should evolve
     */
    STELLAR_TYPE BaseStar::EvolveOnPhase(const double p_dt) {
    if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@entry), p_dt = " << p_dt << ", m_Age= " << m_Age << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << "\n";
        
        STELLAR_TYPE nextStellarType = m_StellarType;                                       // next stellar type - defaults to current
    
        if (ShouldEvolveOnPhase()) {                                                        // should evolve timestep on phase?
                                                                                            // yes
            UpdateMainSequenceCoreMass(p_dt, -m_Mdot);                                      // update core mass, relevant for MS stars
    
    if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@1), p_dt = " << p_dt << ", m_Age= " << m_Age << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";
            m_Tau        = CalculateTauOnPhase();
    if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@2), p_dt = " << p_dt << ", m_Age= " << m_Age << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";
    
            m_COCoreMass = CalculateCOCoreMassOnPhase();
            m_CoreMass   = CalculateCoreMassOnPhase();
            m_HeCoreMass = CalculateHeCoreMassOnPhase();
    
            m_Luminosity = CalculateLuminosityOnPhase();
    
            // Calculate abundances
            m_HydrogenAbundanceCore    = CalculateHAbundanceCore();
            m_HydrogenAbundanceSurface = CalculateHAbundanceSurface();  
            m_HeliumAbundanceCore      = CalculateHeAbundanceCore();
            m_HeliumAbundanceSurface   = CalculateHeAbundanceSurface();
           
    if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@3), p_dt = " << p_dt << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";
            std::tie(m_Radius, nextStellarType) = CalculateRadiusAndStellarTypeOnPhase();   // radius and possibly new stellar type
    if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@4), p_dt = " << p_dt << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";
    
            m_Mu = CalculateHurleyPerturbationMu();
    
            PerturbLuminosityAndRadiusOnPhase();
    
            m_Temperature = CalculateTemperatureOnPhase();
    
            if (p_dt > 0.0) {
                STELLAR_TYPE thisStellarType = ResolveEnvelopeLoss();                       // resolve envelope loss if it occurs - possibly new stellar type
                if (thisStellarType != m_StellarType) {                                     // thisStellarType overrides stellarType (from CalculateRadiusAndStellarTypeOnPhase())
                    nextStellarType = thisStellarType;
                }
                AdvanceAgeAndTime(p_dt);                                                    // advance age of star and simulation time
            }
        }
    
        return nextStellarType;
    }

    
    STELLAR_TYPE BaseStar::EvolveOnPhaseNew(const double p_dt) {
        if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@entry), p_dt = " << p_dt << ", m_Age= " << m_Age << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << "\n";
            
            STELLAR_TYPE nextStellarType = m_StellarType;                                       // next stellar type - defaults to current
        
            if (ShouldEvolveOnPhase()) {                                                        // should evolve timestep on phase?
                                                                                                // yes
                UpdateMainSequenceCoreMass(p_dt, -m_Mdot);                                      // update core mass, relevant for MS stars
        
        if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@1), p_dt = " << p_dt << ", m_Age= " << m_Age << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";
                m_Tau        = CalculateTauOnPhase();
        if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@2), p_dt = " << p_dt << ", m_Age= " << m_Age << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";
        
                m_COCoreMass = CalculateCOCoreMassOnPhase();
                m_CoreMass   = CalculateCoreMassOnPhase();
                m_HeCoreMass = CalculateHeCoreMassOnPhase();
        
                m_Luminosity = CalculateLuminosityOnPhase();
        
                // Calculate abundances
                m_HydrogenAbundanceCore    = CalculateHAbundanceCore(p_Tau, p_InitialHAbundance);
                m_HydrogenAbundanceSurface = CalculateHAbundanceSurface(p_Tau, p_InitialHAbundance);  
                m_HeliumAbundanceCore      = CalculateHeAbundanceCore(GLOBALS->Metallicity(), p_Tau, p_InitialHeAbundance);
                m_HeliumAbundanceSurface   = CalculateHeAbundanceSurface(p_Tau, p_InitialHAbundance);
               
        if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@3), p_dt = " << p_dt << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";
                std::tie(m_Radius, nextStellarType) = CalculateRadiusAndStellarTypeOnPhase();   // radius and possibly new stellar type
        if (OPTIONS->DebugLevel() > 0) std::cout << std::boolalpha << std::setprecision(15) << "BaseStar::EvolveOnPhase(@4), p_dt = " << p_dt << ", m_Mass = " << m_Mass << ", m_radius = " << m_Radius << ", m_Tau = " << m_Tau << "\n";
        
                m_Mu = CalculateHurleyPerturbationMu();
        
                PerturbLuminosityAndRadiusOnPhase();
        
                m_Temperature = CalculateTemperatureOnPhase();
        
                if (p_dt > 0.0) {
                    STELLAR_TYPE thisStellarType = ResolveEnvelopeLoss();                       // resolve envelope loss if it occurs - possibly new stellar type
                    if (thisStellarType != m_StellarType) {                                     // thisStellarType overrides stellarType (from CalculateRadiusAndStellarTypeOnPhase())
                        nextStellarType = thisStellarType;
                    }
                    AdvanceAgeAndTime(p_dt);                                                    // advance age of star and simulation time
                }
            }
        
            return nextStellarType;
        }



    

/*
 * Evolve the star onto the next phase if necessary
 *
 * If the star is at the end of its current phase, recalculate stellar attributes and determine the stellar type
 * for the next phase.  If the star is not at the end of its current phase, do nothing.
 *
 * 
 * If as a result of the evolution of the star it should change stellar type, the new stellar type is returned
 * as the functional return (if no change to stellar type is required the star's current stellar type is returned). 
 * 
 * No change to stellar type is effected here, or as a result of the call to this function.  The caller of this
 * function is expected to check the stellar type returned and manage any stellar type switch required.
 *
 * STELLAR_TYPE ResolveEndOfPhase()
 *
 * @return                                      Stellar Type to which star should evolve
 */
STELLAR_TYPE BaseStar::ResolveEndOfPhase() {
    
    STELLAR_TYPE nextStellarType = m_StellarType;                       // next stellar type - defaults to current

    if (IsEndOfPhase()) {                                               // end of phase?
                                                                        // yes
        m_Tau           = CalculateTauAtPhaseEnd();

        m_COCoreMass    = CalculateCOCoreMassAtPhaseEnd();           //       <<<<<<<<<< DONE >>>>>>>>>>
        m_CoreMass      = CalculateCoreMassAtPhaseEnd();
        m_HeCoreMass    = CalculateHeCoreMassAtPhaseEnd();

        m_Luminosity    = CalculateLuminosityAtPhaseEnd();

        m_Radius        = CalculateRadiusAtPhaseEnd();

        m_Mu            = CalculateHurleyPerturbationMu();

        PerturbLuminosityAndRadiusAtPhaseEnd();

        m_Temperature   = CalculateTemperatureAtPhaseEnd();

        nextStellarType = EvolveToNextPhase();                          // determine the stellar type to which the star should evolve
    }

    return nextStellarType;
}







///////////////////// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// BASE CONSTITUENT STAR FUNCTIONS 

///////////////////// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         COMMON ENVELOPE LAMBDA FUNCTIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/* 
 * CalculateCELambda_Nanjing
 *
 * @brief
 * Calculate the common envelope lambda parameter, per Xu & Li, 2010 (Nanjing - for the university)
 * 
 * Should only be called if the user specified the LAMBDA_NANJING prescription via the commandline
 * option `--common-envelope-lambda-prescription`.
 * 
 * 
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * The assumed ZAMS mass passed to this function should be the actual ZAMS mass of the star, unless
 * the user specified the `--common-envelope-lambda-nanjing-use-rejuvenated-mass` commandline option,
 * in which case the assumed mass should be the "rejuvenated" (or effective) mass of the star.  The
 * corollary is that this function should not be used if the ZAMS mass is not known, unless the user
 * specified the `--common-envelope-lambda-nanjing-use-rejuvenated-mass` commandline option.
 * 
 * 
 * double BaseStar::CalculateCELambda_Nanjing(const double p_Mass, const double p_Radius, const double p_CoreMass) const
 * 
 * @param       p_Mass                          Assumed ZAMS mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Common envelope lambda parameter
 */ 
COMPAS_PURE double BaseStar::CalculateCELambda_Nanjing(const double p_Mass, const double p_Radius, const double p_CoreMass) const {
    double lambda;

    if (OPTIONS->CommonEnvelopeLambdaNanjingEnhanced()) {                                               // use enhanced Nanjing prescription?
                                                                                                        // yes
        // set stellar population based on metallicity
        STELLAR_POPULATION stellarPop = GLOBALS->Metallicity() < LAMBDA_NANJING_ZLIMIT ? STELLAR_POPULATION::POPULATION_II : STELLAR_POPULATION::POPULATION_I;

        if (OPTIONS->CommonEnvelopeLambdaNanjingInterpolateInMass()) {                                  // interpolate across mass models?
            if (OPTIONS->CommonEnvelopeLambdaNanjingInterpolateInMetallicity()) {                       // yes - interpolate across stellar populations?
                lambda = BaseStar::CalculateCELambda_Nanjing_MassAndZInterpolated(GLOBALS->Metallicity(), p_Mass, p_Radius, p_CoreMass);   // yes
            }
            else {                                                                                      // no - not interpolating across stellar populations
                lambda = BaseStar::CalculateCELambda_Nanjing_MassInterpolated(p_Mass, p_Radius, p_CoreMass, stellarPop);
            }
        }
        else {                                                                                          // no - not interpolating across mass models
            // calculate index into NANJING_MASSES (see constants.h) for p_Mass
                 if (p_Mass < NANJING_MASSES_MIDPOINTS[0])      massIndex = 0;                          // use lambda for the 1 Msun model
            else if (p_Mass >= NANJING_MASSES_MIDPOINTS.back()) massIndex = NANJING_MASSES.size() - 1;  // use lambda for the 100 Msun model
            else                                                massIndex = utils::BinarySearch(NANJING_MASSES_MIDPOINTS, p_Mass)[1]; // bin edge indices - use upper

            if (OPTIONS->CommonEnvelopeLambdaNanjingInterpolateInMetallicity()) {                       // interpolate across stellar populations?
                lambda = BaseStar::CalculateCELambda_Nanjing_ZInterpolated(GLOBALS->Metallicity(), p_Mass, p_Radius, p_CoreMass, GLOBALS->SigmaHurley(), massIndex); // yes
            }
            else {                                                                                      // no - not interpolating across stellar populations
                lambda = BaseStar::CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, massIndex, stellarPop);
            }
        }
    }
    else {                                                                                              // no - not enhanced Nanjing prescription: use StarTrack prescription
        lambda = CalculateLambdaNanjingStarTrack(p_Mass, m_Metallicity);
    }

    return lambda;
}


/* 
 * CalculateCELambda_Nanjing_MassAndZInterpolated
 *
 * @brief
 * Calculate CE lambda, per Xu & Li 2010, by interpolating across different mass models,
 * and metallicities (stellar populations), where necessary.
 * 
 * 
 * double CalculateCELambda_Nanjing_MassAndZInterpolated(const double p_Mass, const double p_Radius, const double p_CoreMass) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Common envelope lambda parameter
 */ 
double BaseStar::CalculateCELambda_Nanjing_MassAndZInterpolated(const double p_Mass, const double p_Radius, const double p_CoreMass) const {

    double lambda;

         if (GLOBALS->SigmaHurley() <= LAMBDA_NANJING_POPII_Z) lambda = BaseStar::CalculateCELambda_Nanjing_MassInterpolated(p_Mass, p_Radius, p_CoreMass, STELLAR_POPULATION::POPULATION_II);   // use lambda for pop. II metallicity
    else if (GLOBALS->SigmaHurley() >= LAMBDA_NANJING_POPI_Z)  lambda = BaseStar::CalculateCELambda_Nanjing_MassInterpolated(p_Mass, p_Radius, p_CoreMass, STELLAR_POPULATION::POPULATION_I);    // use lambda for pop. I metallicity
    else {
        // linear interpolation in logZ between pop. I and pop. II metallicities
        double lambdaLow = BaseStar::CalculateCELambda_Nanjing_MassInterpolated(p_Mass, p_Radius, p_CoreMass, STELLAR_POPULATION::POPULATION_II);
        double lambdaUp  = BaseStar::CalculateCELambda_Nanjing_MassInterpolated(p_Mass, p_Radius, p_CoreMass, STELLAR_POPULATION::POPULATION_I);

        lambda = lambdaLow + (m_Log10Metallicity - LAMBDA_NANJING_POPII_LOGZ) / (LAMBDA_NANJING_POPI_LOGZ - LAMBDA_NANJING_POPII_LOGZ) * (lambdaUp - lambdaLow);
    }

    return lambda;
}


/* 
 * CalculateCELambda_Nanjing_MassInterpolated
 *
 * @brief
 * Calculate CE lambda, per Xu & Li 2010, by interpolating across different mass models,
 * where necessary, for a given metallicity (stellar population).
 * 
 * 
 * double CalculateCELambda_Nanjing_MassInterpolated(const double p_Mass, const double p_Radius, const double p_CoreMass, const STELLAR_POPULATION p_StellarPop) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_StellarPop                    Stellar population (POP I or POP II)
 * @return                                      Common envelope lambda parameter
 */ 
double BaseStar::CalculateCELambda_Nanjing_MassInterpolated(const double p_Mass, const double p_Radius, const double p_CoreMass, const STELLAR_POPULATION p_StellarPop) const {

    INT_VECTOR indices = utils::BinarySearch(NANJING_MASSES, p_Mass);                                   // find mass in NANJING_MASSES
    size_t lower = indices[0];                                                                          // bin lower bound
    size_t upper = indices[1];                                                                          // bin upper bound

    double lambda;
    if ((lower < 0) && (upper >= 0)) {                                                                  // mass below range calculated by Xu & Li 2010?
        lambda = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, 0, p_StellarPop);     // yes - return lambda for minimum mass
    }
    else if ((lower >= 0) && (upper < 0)) {                                                             // mass above range calculated by Xu & Li 2010?
        lambda = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, NANJING_MASSES.size() - 1, p_StellarPop); // yes - return lambda for maximum mass
    }
    else if (lower == upper) {                                                                          // mass is exactly equal to the mass of a model from Xu & Li 2010?
        lambda = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, lower, p_StellarPop); // yes - return lambda for mass
    }
    else {                                                                                              // mass between bins
        // linear interpolation between upper and lower mass bins
        const double lambdaLower = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, lower, p_StellarPop);
        const double lambdaUpper = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, upper, p_StellarPop);

        lambda = lambdaLower + (p_Mass - NANJING_MASSES[lower]) / (NANJING_MASSES[upper] - NANJING_MASSES[lower]) * (lambdaUpper - lambdaLower);
    }

    return lambda;
}


/* 
 * CalculateCELambda_Nanjing_ZInterpolated
 *
 * @brief
 * Calculate CE lambda, per Xu & Li 2010, by interpolating across different metallicities
 * (stellar populations), where necessary, for a given mass.
 * 
 * double CalculateCELambda_Nanjing_ZInterpolated(const double p_Mass, const double p_Radius, const double p_CoreMass, const size_t p_MassIndex) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MassIndex                     Index of mass bin in NANJING_MASSES (see constants.h)
 * @return                                      Common envelope lambda parameter
 */ 
double BaseStar::CalculateCELambda_Nanjing_ZInterpolated(const double p_Mass, const double p_Radius, const double p_CoreMass, const size_t p_MassIndex) const {
    double lambda;
    
         if (GLOBALS->Metallicity() < LAMBDA_NANJING_POPII_Z) lambda = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, p_MassIndex, STELLAR_POPULATION::POPULATION_II);   // lambda for pop. II metallicity
    else if (GLOBALS->Metallicity() > LAMBDA_NANJING_POPI_Z)  lambda = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, p_MassIndex, STELLAR_POPULATION::POPULATION_I);    // lambda for pop. I metallicity
    else {
        // linear interpolation in logZ between pop. I and pop. II metallicities
        const double lambdaLow = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, p_MassIndex, STELLAR_POPULATION::POPULATION_II);
        const double lambdaUp  = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, p_MassIndex, STELLAR_POPULATION::POPULATION_I);

        lambda = lambdaLow + (GLOBALS->SigmaHurley() - LAMBDA_NANJING_POPII_LOGZ) / (LAMBDA_NANJING_POPI_LOGZ - LAMBDA_NANJING_POPII_LOGZ) * (lambdaUp - lambdaLow);
    }

    return lambda;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  ZETA FUNCTIONS                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateZetaAdiabatic
 *
 * @brief
 * Calculate the adiabatic exponent (donor radial response to mass loss), zeta
 * 
 *
 * double BaseStar::CalculateZetaAdiabatic() const
 *
 * @return                                      Adiabatic exponent, zeta
 */
double BaseStar::CalculateZetaAdiabatic() const { 
                                                                                
    double zeta;

    switch (OPTIONS->StellarZetaPrescription()) {

        case ZETA_PRESCRIPTION::SOBERMAN: 
        case ZETA_PRESCRIPTION::HURLEY:   
        case ZETA_PRESCRIPTION::ARBITRARY:
            zeta = CalculateZetaAdiabatic_ByEnvelopeType(OPTIONS->StellarZetaPrescription());
            break;

        default:                                                                        // unknown prescription
            // the only way this can happen is if someone added a ZETA_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_ZETA_PRESCRIPTION);                              // throw error
    }

    return zeta;
}


/*
 * CalculateZetaAdiabatic_Soberman1997
 *
 * @brief
 * Calculate the adiabatic exponent (donor radial response to mass loss), zeta,
 * per Soberman, Phinney, van den Heuvel 1997, eq 61
 *
 *
 * double CalculateZetaAdiabatic_Soberman1997(const double p_Mass, const double p_CoreMasss) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Adiabatic exponent, zeta
 */
double BaseStar::CalculateZetaAdiabatic_Soberman1997(const double p_Mass, const double p_CoreMass) const {
    
    if (p_CoreMass >= p_Mass) return 0.0;               // if the object is all core, the calculation is meaningless

    const double mRatio    = p_CoreMass / m_Mass;       // Soberman, Phinney, van den Heuvel 1997, eq 57
    const double oneMinusM = 1.0 - mRatio;

    // Soberman, Phinney, van den Heuvel 1997, eq 61
    return ((2.0 / 3.0) * mRatio / oneMinusM) - ((1.0 / 3.0) * (oneMinusM / (1.0 + (mRatio + mRatio)))) - (0.03 * mRatio) + (0.2 * mRatio / (1.0 + (1.0 / utils::IntPow(oneMinusM, 6)))); 
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           CRITICAL MASS RATIO FUNCTIONS                           //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCriticalMassRatio
 *
 * @brief
 * Calculate the critical mass ratio for unstable mass transfer, based on user-specified
 * commandline option `--critical-mass-ratio-prescription`.
 *
 * 
 * double BaseStar::CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate, const double p_MTefficiency) const
 *
 * @param       p_AccretorIsDegenerate          Flag to indicate if the accretor is a degenerate star
 * @param       p_MTefficiency                  Mass transfer accretion efficiency
 * @return                                      Critical mass ratio
 */
double BaseStar::CalculateCriticalMassRatio(const bool p_AccretorIsDegenerate, const double p_MTefficiency) const {
    
    double qCrit = 0.0;                                                                 // default return value

    switch (OPTIONS->QCritPrescription()) {                                             // which QCRIT_PRESCRIPTION?
                
        case QCRIT_PRESCRIPTION::GE:                                                    // Ge (adiabatic assumption)
        case QCRIT_PRESCRIPTION::GE_IC:                                                 // Ge (isentropic envelope assumption)
            qCrit = CalculateCriticalMassRatio_Ge2020(p_MTefficiency);       // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NOT PURE/CONST  uses state (BaseStar)
            break;

        case QCRIT_PRESCRIPTION::CLAEYS:                                                // CLAEYS
            qCrit = CalculateCriticalMassRatio_Claeys2014(p_AccretorIsDegenerate);    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NOT PURE/CONST  uses state (GiantBranch)
            break;

        case QCRIT_PRESCRIPTION::HURLEY_HJELLMING_WEBBINK:                              // HURLEY_HJELLMING_WEBBINK
            qCrit = CalculateCriticalMassRatio_Hurley2002();    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NOT PURE/CONST  uses state (GiantBranch)
            break;
        
        default:                                                                        // unexpected prescription
            // the only way this can happen is if the QCRIT_PRESCRIPTION specified by the user is not
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user specified a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.  We know this function does't account for all QCRIT_PRESCRIPTIONs,
            // but it (currently) *should* account for all except QCRIT_PRESCRIPTION::NONE.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to ensure this function is not called when the user specifies
            // a QCRIT_PRESCRIPTION that is not handled by this function.

            THROW_ERROR(ERROR::UNEXPECTED_QCRIT_PRESCRIPTION);                          // throw error
    }

    return qCrit;
}




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                              MASS TRANSFER FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateMassAcceptanceRate
 *
 * @brief
 * Calculate:
 *
 *     (a) the maximum mass acceptance rate of this star, as the accretor, during mass transfer, and
 *     (b) the accretion efficiency parameter
 *
 * The maximum acceptance rate of the accretor star during mass transfer is based on stellar type: this function
 * is for main sequence stars, or stars that have evolved off of the main sequence but are not yet remnants.
 *
 * Mass transfer is assumed Eddington limited for BHs and NSs.  The formalism of Nomoto/Claeys is used for WDs.
 *
 * For non compact objects:
 *
 *    1) Kelvin-Helmholtz (thermal) timescale if THERMAL (thermally limited) mass transfer efficiency
 *    2) Choose a fraction of the mass rate that will be effectively accreted for FIXED fraction mass transfer (as in StarTrack)
 *
 * Assumes this star is the accretor.
 * 
 *
 * DBL_DBL CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate of the donor
 * @param   [IN]    p_AccretorMassRate          Thermal mass transfer rate of the accretor (this star)
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate and the Accretion Efficiency Parameter
 */
DBL_DBL BaseStar::CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate) {

    double acceptanceRate   = 0.0;                                                                      // acceptance mass rate - default = 0.0
    double fractionAccreted = 0.0;                                                                      // accretion fraction - default  = 0.0

    switch (OPTIONS->MassTransferAccretionEfficiencyPrescription()) {

        case MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED:                                   // thermally limited mass transfer
            acceptanceRate   = std::min(OPTIONS->MassTransferCParameter() * p_AccretorMassRate, p_DonorMassRate);
            fractionAccreted = acceptanceRate / p_DonorMassRate;
            break;

        case MT_ACCRETION_EFFICIENCY_PRESCRIPTION::HAMSTARS:                                            // thermally limited mass transfer, following Lau+, 2024
            // the mass transfer C parameter is fit using the data from Figure (1) of Lau+, 2024
            acceptanceRate   = std::min(PPOW(10.0, (4.0 / PPOW((Mass() + 0.2), 0.3) - 0.6)) * p_AccretorMassRate, p_DonorMassRate);
            fractionAccreted = acceptanceRate / p_DonorMassRate;
            break;

        case MT_ACCRETION_EFFICIENCY_PRESCRIPTION::FIXED_FRACTION:                                      // fixed fraction of mass accreted, as in StarTrack
            fractionAccreted = OPTIONS->MassTransferFractionAccreted();
            acceptanceRate   = std::min(p_DonorMassRate, fractionAccreted * p_DonorMassRate);
            break;

        default:                                                                                        // unknown prescription
            // the only way this can happen is if someone added an MT_ACCRETION_EFFICIENCY_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MT_ACCRETION_EFFICIENCY_PRESCRIPTION);                           // throw error
    }

    return std::make_tuple(acceptanceRate, fractionAccreted);
}


/*
 * CalculateThermalMassAcceptanceRate
 *
 * Calculate thermal mass acceptance rate
 *
 *
 * double CalculateThermalMassAcceptanceRate(const double p_Radius)
 *
 * @param   [IN]    p_Radius                    Radius of the accretor (Rsol) [typically called with Roche Lobe radius]
 * @return                                      Thermal mass acceptance rate
 */
double BaseStar::CalculateThermalMassAcceptanceRate(const double p_Radius) {
    double acceptanceRate = 0.0;                                                                    // thermal mass acceptance rate

    switch (OPTIONS->MassTransferThermallyLimitedVariation()) {
        case MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE:
            acceptanceRate = (m_Mass - m_CoreMass) / CalculateTimescale_Thermal(p_Mass, p_Radius, p_CoreMass);           // uses provided accretor radius (should be Roche lobe radius in practice)
            break;
        case MT_THERMALLY_LIMITED_VARIATION::C_FACTOR:
            acceptanceRate = CalculateMLrateThermal();
            break;

        default:                                                                                    // unknown prescription
            // the only way this can happen is if someone added an MT_THERMALLY_LIMITED_VARIATION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_MT_THERMALLY_LIMITED_VARIATION);                             // throw error
    }

    return acceptanceRate;
}


/*
 * CalculateMassAccretedForCO
 *
 * Calculate the mass accreted by a Neutron Star given mass and radius of companion
 *
 * JR: todo: flesh-out this documentation
 * JR: is this just for NS?  Maybe a change of name... Or move it to NS.cpp? *Ilya* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 *
 *
 * double CalculateMassAccretedForCO(const double p_Mass,
 *                                   const double p_CompanionMass,
 *                                   const double p_CompanionRadius,
 *                                   const double p_CompanionEnvelope)
 *
 * @param   [IN]    p_Mass                      The mass of the accreting star (Msol)
 * @param   [IN]    p_CompanionMass             The mass of the companion star (Msol)
 * @param   [IN]    p_CompanionRadius           The radius of the companion star (Rsol)
 * @param   [IN]    p_CompanionEnvelope         Envelope of the companion pre-CE
 * @return                                      Mass accreted by the Neutron Star (Msol)
 */
double BaseStar::CalculateMassAccretedForCO(const double p_Mass,
                                            const double p_CompanionMass,
                                            const double p_CompanionRadius,
                                            const double p_CompanionEnvelope) const {

     double deltaMass;

     switch (OPTIONS->CommonEnvelopeMassAccretionPrescription()) {                                              // which prescription?

        case CE_ACCRETION_PRESCRIPTION::ZERO:                                                                   // ZERO
            deltaMass = 0.0;
            break;

        case CE_ACCRETION_PRESCRIPTION::CONSTANT:                                                               // CONSTANT
            deltaMass = OPTIONS->CommonEnvelopeMassAccretionConstant();                                         // use program option
            break;

        case CE_ACCRETION_PRESCRIPTION::UNIFORM:                                                                // UNIFORM
            deltaMass = RAND->Random(OPTIONS->CommonEnvelopeMassAccretionMin(), OPTIONS->CommonEnvelopeMassAccretionMax()); // uniform random distribution - Oslowski+ (2011)
            break;

        case CE_ACCRETION_PRESCRIPTION::MACLEOD: {                                                              // MACLEOD
                                                                                                                // linear regression estimated from Macleod+ (2015)
            double mm = -1.0714285714285712E-05;                                                                // gradient of the linear fit for gradient
            double cm =  0.00012057142857142856;                                                                // intercept of the linear fit for gradient
            double mc =  0.01588571428571428;                                                                   // gradient of the linear fit for intercept
            double cc = -0.15462857142857137;                                                                   // intercept of the linear fir for intercept
            double m  = mm * p_CompanionMass + cm;                                                              // gradient of linear fit for mass
            double c  = mc * p_CompanionMass + cc;                                                              // intercept of linear fit for mass

            // calculate mass to accrete and clamp to minimum and maximum from program options
            deltaMass = std::min(OPTIONS->CommonEnvelopeMassAccretionMax(), std::max(OPTIONS->CommonEnvelopeMassAccretionMin(), m * p_CompanionRadius + c));
            } break;

        case CE_ACCRETION_PRESCRIPTION::CHEVALIER:                                                              // CHEVALIER
                                                                                                                // Model 2 from van Son et al. 2020
            deltaMass = (p_Mass * p_CompanionMass) / (2.0 * (p_Mass + p_CompanionMass)) ;                       // Hoyle-Lyttleton accretion rate times inspiral time
            break;

        default:                                                                                                // unknown prescription
            // the only way this can happen is if someone added a CE_ACCRETION_PRESCRIPTION
            // and it isn't accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_CE_ACCRETION_PRESCRIPTION);                                              // throw error
    }

    deltaMass = std::min(p_CompanionEnvelope, deltaMass);                                                       // clamp the mass accretion to be no more than the envelope of the companion pre CE

    return deltaMass;
}


/*
 * MassTransferDonorHistory
 *
 * Construct string representing the Mass Transfer Donor History vector
 *
 * This is so that a string is passed to the output, not a vector of stellar types.
 *
 * std::string BaseStar::MassTransferDonorHistory() 
 *
 * @return                              string of dash-separated Hurley stellar type numbers
 */
std::string BaseStar::MassTransferDonorHistory() const {

    ST_VECTOR   mtHistVec = m_MassTransferDonorHistory;      
    std::string mtHistStr = "";

    if (mtHistVec.empty()) {                                                            // this star was never a donor for MT
        mtHistStr = "NA";
    }
    else {                                                                              // this star was a donor, return the stellar type string
        for (size_t ii = 0; ii < mtHistVec.size(); ii++) {
            mtHistStr += std::to_string(static_cast<int>(mtHistVec[ii])) + "-";         // create string of stellar type followed by dash
        }
        mtHistStr.pop_back();                                                           // remove final dash
    }
    return mtHistStr;
}


/*
 * UpdateMassTransferDonorHistory
 *
 * Add new MT event to event history - only for donor stars
 *
 * void BaseStar::UpdateMassTransferDonorHistory()
 *
 */
void BaseStar::UpdateMassTransferDonorHistory() {

    if (m_MassTransferDonorHistory.empty()) {                                           // no history?
        m_MassTransferDonorHistory.push_back(m_StellarType);                            // yes - first event
    }
    else if (!utils::IsOneOf(m_StellarType, { m_MassTransferDonorHistory.back() })) {   // first MT as current stellar type?
        m_MassTransferDonorHistory.push_back(m_StellarType);                            // yes - new event
    }
}
















/*
 * Calculate the Dynamical tides contribution to the l=2, (n,m) = [(1,0), (1,2), (2,2), (3,2)] imaginary components of the 
 * potential tidal Love number
 *
 * Gravity Waves, Core boundary:
 * Zahn, 1977, Eq. (5.5) , with the value of E_2 coming from Kushnir et al., 2017, by comparing Eq. (8) to Eq. (1)
 * 
 * Gravity Waves, Envelope Boundary:
 * Ahuir et al, 2021, Eq. (131)
 * 
 * Inertial Waves, Convective Envelope:
 * Ogilvie, 2013, Eq. (B3)
 * 
 * DBL_DBL_DBL_DBL CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)
 *
 * @param   [IN]    p_Omega                     Orbital angular frequency (1/yr)
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of binary (AU)
 * @param   [IN]    p_M2                        Mass of companion star (Msol)
 * @return                                      [(1,0), (1,2), (2,2), (3,2)] Imaginary components of the 
 *                                              potential tidal Love number, Dynamical tides only (unitless)
 */
DBL_DBL_DBL_DBL BaseStar::CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const {
    
    double coreMass = CalculateConvectiveCoreMass();

    double envMass, envMassMax;
    std::tie(envMass, envMassMax) = CalculateConvectiveEnvelopeMass();
    
    double radIntershellMass = m_Mass - coreMass - envMass;                                             // refers to the combined mass of non-convective layers

    double radiusAU              = m_Radius * RSOL_TO_AU;
    double coreRadiusAU          = CalculateConvectiveCoreRadius() * RSOL_TO_AU;
    double convectiveEnvRadiusAU = CalculateConvectiveEnvelopeRadialExtent() * RSOL_TO_AU;
    double radiusIntershellAU    = radiusAU - convectiveEnvRadiusAU;                                    // Outer radial coordinate of radiative intershell

    // There should be no Dynamical tides if the entire star is convective, i.e. if there are no convective-radiative boundaries. 
    // If so, return 0.0 for all dynamical components of ImKnm.
    // This condition should be true for low-mass MS stars (<= 0.35 Msol) at ZAMS.
    if (utils::Compare(radIntershellMass/m_Mass, TIDES_MINIMUM_FRACTIONAL_EXTENT) <= 0 || utils::Compare(radiusIntershellAU, coreRadiusAU) <= 0) {
        return std::make_tuple(0.0, 0.0, 0.0, 0.0);                           
    }
 
    double R_3            = radiusAU * radiusAU * radiusAU;
    double R3OverG_M      = (R_3 / G_AU_Msol_yr / m_Mass);
    double sqrtR3OverG_M  = std::sqrt(R3OverG_M);

    double k10GravityCore = 0.0;                                                                        // gravity Wave dissipation, core boundary
    double k12GravityCore = 0.0;
    double k22GravityCore = 0.0;
    double k32GravityCore = 0.0;

    double k10GravityEnv  = 0.0;                                                                        // gravity Wave dissipation, envelope boundary
    double k12GravityEnv  = 0.0;
    double k22GravityEnv  = 0.0;
    double k32GravityEnv  = 0.0;

    double k22InertialEnv = 0.0;                                                                        // inertial Wave dissipation, envelope
    
    double omegaSpin      = Omega();
    double twoOmegaSpin   = omegaSpin + omegaSpin;

    double w10 = p_Omega;
    double w12 = ((p_Omega) - twoOmegaSpin);
    double w22 = ((p_Omega + p_Omega) - twoOmegaSpin);
    double w32 = ((p_Omega + p_Omega + p_Omega) - twoOmegaSpin);
        
    // Assume that GW dissipation from core boundary is only efficient if the radiative region extends to the surface, i.e. there is no convective envelope.
    if (utils::Compare(coreRadiusAU/radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0 && utils::Compare(coreMass/m_Mass, TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0 && utils::Compare(convectiveEnvRadiusAU/radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) < 0 && utils::Compare(envMass/m_Mass, TIDES_MINIMUM_FRACTIONAL_EXTENT) < 0) {                   
        constexpr double beta2Dynamical         = 1.0;
        constexpr double rhoFactorDynamcial     = 0.1;
        double coreRadiusOverRadius   = coreRadiusAU / radiusAU;
        double coreRadiusOverRadius_3 = coreRadiusOverRadius * coreRadiusOverRadius * coreRadiusOverRadius;
        double coreRadiusOverRadius_9 = coreRadiusOverRadius_3 * coreRadiusOverRadius_3 * coreRadiusOverRadius_3;
        double massOverCoreMass       = m_Mass / coreMass;
        double E2Dynamical            = (2.0 / 3.0) * coreRadiusOverRadius_9 * massOverCoreMass * std::cbrt(massOverCoreMass) * beta2Dynamical * rhoFactorDynamcial;

        // (l=2, n=1, m=0), Gravity Wave dissipation from core boundary
        double s10     = w10 * sqrtR3OverG_M;
        double s10_4_3 = s10 * std::cbrt(s10);
        double s10_8_3 = s10_4_3 * s10_4_3;
        k10GravityCore = E2Dynamical *  std::copysign(s10_8_3, w10);
        if (std::isnan(k10GravityCore)) k10GravityCore = 0.0;

        // (l=2, n=1, m=2), Gravity Wave dissipation from core boundary
        double s12     = w12 * sqrtR3OverG_M;
        double s12_4_3 = s12 * std::cbrt(s12);
        double s12_8_3 = s12_4_3 * s12_4_3;
        k12GravityCore = E2Dynamical * std::copysign(s12_8_3, w12);
        if (std::isnan(k12GravityCore)) k12GravityCore = 0.0;

        // (l=2, n=2, m=2), Gravity Wave dissipation from core boundary
        double s22     = w22 * sqrtR3OverG_M;
        double s22_4_3 = s22 * std::cbrt(s22);
        double s22_8_3 = s22_4_3 * s22_4_3;
        k22GravityCore = E2Dynamical * std::copysign(s22_8_3, w22);
        if (std::isnan(k22GravityCore)) k22GravityCore = 0.0;

        // (l=2, n=3, m=2), Gravity Wave dissipation from core boundary
        double s32     = w32 * sqrtR3OverG_M;
        double s32_4_3 = s32 * std::cbrt(s32);
        double s32_8_3 = s32_4_3 * s32_4_3;
        k32GravityCore = E2Dynamical * std::copysign(s32_8_3, w32);
        if (std::isnan(k32GravityCore)) k32GravityCore = 0.0;    
    }

    double rint_3 = radiusIntershellAU * radiusIntershellAU * radiusIntershellAU;
    double rc_3   = coreRadiusAU * coreRadiusAU * coreRadiusAU;
    double gamma  = (envMass / (R_3 - rint_3)) / (radIntershellMass / (rint_3 - rc_3));

    // There is no GW or IW dissipation from the envelope boundary if no convective envelope
    if ((utils::Compare(convectiveEnvRadiusAU / radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0) || (utils::Compare(envMass / m_Mass, TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0)) {    

        constexpr double dynPrefactor     = 3.207452512782476;                                                        // 3^(11/3) * Gamma(1/3)^2 / 40 PI
        constexpr double m_l_factor_22    = 0.183440402716368;                                                        // m * (l(l+1))^{-4/3}
        double cbrtdNdlnr       = std::cbrt(G_AU_Msol_yr * radIntershellMass / radiusIntershellAU / (radiusAU - radiusIntershellAU) / (radiusAU - radiusIntershellAU));
        
        double alpha            = radiusIntershellAU / radiusAU;
        double oneMinusAlpha    = 1.0 - alpha;
        double beta             = radIntershellMass / m_Mass;

        double alpha_2          = alpha * alpha;
        double alpha_3          = alpha_2 * alpha;
        double alpha_5          = alpha_3 * alpha_2;
        double alpha_11         = alpha_5 * alpha_5 * alpha;
        double oneMinusAlpha_2  = oneMinusAlpha * oneMinusAlpha;
        double oneMinusAlpha_3  = 1.0 - alpha_3;
        double beta_2           = beta * beta;

    
        double oneMinusGamma    = 1.0 - gamma;
        double oneMinusGamma_2  = oneMinusGamma * oneMinusGamma;
        double alpha_2_3Minus_1 = (alpha * 2.0 / 3.0) - 1.0;

        // Assume GW dissipation from the envelope boundary only acts if the radiative zone extends to the core, i.e. if there is no convective core.
        if (utils::Compare(coreRadiusAU/radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) < 0) {                              
            double Epsilon       = alpha_11 * envMass / m_Mass * oneMinusGamma_2 * alpha_2_3Minus_1 * alpha_2_3Minus_1 / beta_2 / oneMinusAlpha_3 / oneMinusAlpha_2;

            // (l=2, n=1, m=0), Gravity Wave dissipation from envelope boundary is always 0.0 since m * (l(l+1))^{-4/3} = 0

            // (l=2, n=1, m=2), Gravity Wave dissipation from envelope boundary
            double w12_4_3       = w12 * std::cbrt(w12);
            double w12_8_3       = w12_4_3 * w12_4_3;
            k12GravityEnv        = dynPrefactor * m_l_factor_22 * std::copysign(w12_8_3, w12) * R3OverG_M * Epsilon / cbrtdNdlnr;
            if (std::isnan(k12GravityEnv)) k12GravityEnv = 0.0;  

            // (l=2, n=2, m=2), Gravity Wave dissipation from envelope boundary
            double w22_4_3       = w22 * std::cbrt(w22);
            double w22_8_3       = w22_4_3 * w22_4_3;
            k22GravityEnv        = dynPrefactor * m_l_factor_22 * std::copysign(w22_8_3, w22)* R3OverG_M * Epsilon / cbrtdNdlnr;
            if (std::isnan(k22GravityEnv)) k22GravityEnv = 0.0;  

            // (l=2, n=3, m=2), Gravity Wave dissipation from envelope boundary
            double w32_4_3       = w32 * std::cbrt(w32);
            double w32_8_3       = w32_4_3 * w32_4_3;
            k32GravityEnv        = dynPrefactor * m_l_factor_22 * std::copysign(w32_8_3, w32) * R3OverG_M * Epsilon / cbrtdNdlnr;
            if (std::isnan(k32GravityEnv)) k32GravityEnv = 0.0;  
        }

        // (l=2, n=2, m=2), Inertial Wave dissipation, convective envelope
        // IW dissipation is only efficient for highly spinning stars, as in Esseldeurs, et al., 2024 
        if (utils::Compare(twoOmegaSpin, p_Omega) >= 0) {                                                                            
            double epsilonIW_2       = omegaSpin * omegaSpin * R3OverG_M;
            double oneMinusAlpha_4 = oneMinusAlpha_2 * oneMinusAlpha_2;
            double bracket1          = 1.0 + (2.0 * alpha) + (3.0 * alpha_2) + (3.0 * alpha_3 / 2.0);
            double bracket2          = 1.0 + (oneMinusGamma / gamma) * alpha_3;
            double bracket3          = 1.0 + (3.0 * gamma / 2.0) + (5.0 * alpha_3 / (2.0 * gamma) * (1.0 + (gamma / 2.0) - (3.0* gamma * gamma / 2.0))) - (9.0 / 4.0 * oneMinusGamma * alpha_5);
            k22InertialEnv           = (100.0 * M_PI / 63.0) * epsilonIW_2 * (alpha_5 / (1.0 - alpha_5)) * oneMinusGamma_2 * oneMinusAlpha_4 * bracket1 * bracket1 * bracket2 / bracket3 / bracket3;
            k22InertialEnv           = std::copysign(k22InertialEnv, w22);
            if (std::isnan(k22InertialEnv)) k22InertialEnv = 0.0;  
        }
    }

    // return ImKnmDynamical
    return std::make_tuple(k10GravityCore + k10GravityEnv, k12GravityCore + k12GravityEnv, k22GravityCore + k22GravityEnv + k22InertialEnv, k32GravityCore + k32GravityEnv);
}


/*
 * Calculate the Equilibrium tides contribution to the l=2, (n,m) = [(1,0), (1,2), (2,2), (3,2)] imaginary components of the 
 * potential tidal Love number
 * 
 * Barker (2020), Eqs. (20) to (27), (l=2, n=2, m=2 mode only).
 *
 * DBL_DBL_DBL_DBL CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)
 *
 * @param   [IN]    p_Omega                     Orbital angular frequency (1/yr)
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of binary (AU)
 * @param   [IN]    p_M2                        Mass of companion star (Msol)
 * @return                                      [(1,0), (1,2), (2,2), (3,2)] Imaginary components of the 
 *                                              potential tidal Love number, Equilibrium tides only (unitless)
 */
DBL_DBL_DBL_DBL BaseStar::CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const {

    // Viscous dissipation
    // No contribution from convective core; only convective envelope.

    double rEnvAU = CalculateConvectiveEnvelopeRadialExtent() * RSOL_TO_AU;
    double envMass, envMassMax;
    std::tie(envMass, envMassMax) = CalculateConvectiveEnvelopeMass();
    
    double rOutAU = m_Radius * RSOL_TO_AU;                                                      // outer boundary of convective envelope
    double rInAU  = (rOutAU - rEnvAU);                                                          // inner boundary of convective envelope
    
    if (utils::Compare(rEnvAU/rOutAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) <= 0 || utils::Compare(envMass/m_Mass, TIDES_MINIMUM_FRACTIONAL_EXTENT) <= 0 || std::isnan(envMass)) return std::make_tuple(0.0, 0.0, 0.0, 0.0);           // skip calculations if there is no convective envelope (to avoid Imk22 = NaN)

    double rOut_2  = rOutAU * rOutAU;
    double rOut_3  = rOut_2 * rOutAU;
    double rOut_5  = rOut_2 * rOut_3;
    double rOut_7  = rOut_2 * rOut_5;
    double rOut_9  = rOut_2 * rOut_7;

    double rIn_2  = rInAU * rInAU;
    double rIn_3  = rIn_2 * rInAU;
    double rIn_5  = rIn_2 * rIn_3;
    double rIn_7  = rIn_2 * rIn_5;
    double rIn_9  = rIn_2 * rIn_7;

    double omegaSpin      = Omega();
    double twoOmegaSpin   = omegaSpin + omegaSpin;

    double rhoConv        = envMass / (4.0 * M_PI * (rOut_3 - rIn_3) / 3.0);
    double lConv          = rEnvAU;                                                              // set length scale to height of convective envelope
    double tConv          = CalculateEddyTurnoverTimescale();
    double vConv          = lConv / tConv;
    double omegaConv      = 1.0 / tConv;                                                         // absent factor of 2*PI, following Barker (2020)
    double vl             = vConv * lConv;
    double M_2            = m_Mass * m_Mass;

    double vl_5           = 5.0 * vl;
    double vl25OverRoot20 = vl * (25.0 / std::sqrt(20.0));
    double vlOver2        = 0.5 * vl;

    double w10 = p_Omega;
    double w12 = ((p_Omega) - (twoOmegaSpin));
    double w22 = ((p_Omega + p_Omega) - (twoOmegaSpin));
    double w32 = ((p_Omega + p_Omega + p_Omega) - (twoOmegaSpin));

    double k2_prefactor   = (224.0 * M_PI / 15.0) * (rOut_9 - rIn_9) * rhoConv / G_AU_Msol_yr / M_2 / rOut_5;

    // (l=2, n=1, m=0), Viscous dissipation, convective envelope
    double omega_t_10            = std::abs(w10);                                               
    double omega_tOverOmega_c_10 = omega_t_10 / omegaConv;
    double nuTidal10             = vl_5;
    if (utils::Compare(omega_tOverOmega_c_10, 5.0) > 0) {             
        nuTidal10 = vl25OverRoot20 / omega_tOverOmega_c_10 / omega_tOverOmega_c_10;
    }
    else if (utils::Compare(omega_tOverOmega_c_10, 0.01) > 0) {
        nuTidal10 = vlOver2 / std::sqrt(omega_tOverOmega_c_10);    
    }
    double k10Equilibrium    = k2_prefactor * nuTidal10 * omega_t_10;
    if (std::isnan(k10Equilibrium)) k10Equilibrium = 0.0;
    if (w10 < 0.0) k10Equilibrium = -std::abs(k10Equilibrium);


    // (l=2, n=1, m=2), Viscous dissipation, convective envelope
    double omega_t_12            = std::abs(w12);                                               
    double omega_tOverOmega_c_12 = omega_t_12 / omegaConv;
    double nuTidal12             = vl_5;
    if (utils::Compare(omega_tOverOmega_c_12, 5.0) > 0) {             
        nuTidal12 = vl25OverRoot20 / omega_tOverOmega_c_12 / omega_tOverOmega_c_12;
    }
    else if (utils::Compare(omega_tOverOmega_c_12, 0.01) > 0) {
        nuTidal12 = vlOver2 / std::sqrt(omega_tOverOmega_c_12);    
    }
    double k12Equilibrium    = k2_prefactor * nuTidal12 * omega_t_12;
    if (std::isnan(k12Equilibrium)) k12Equilibrium = 0.0;
    if (w12 < 0) k12Equilibrium = -std::abs(k12Equilibrium);


    // (l=2, n=2, m=2), Viscous dissipation, convective envelope
    double omega_t_22            = std::abs(w22);                                               
    double omega_tOverOmega_c_22 = omega_t_22 / omegaConv;
    double nuTidal22             = vl_5;
    if (utils::Compare(omega_tOverOmega_c_22, 5.0) > 0) {             
        nuTidal22 = vl25OverRoot20 / omega_tOverOmega_c_22 / omega_tOverOmega_c_22;
    }
    else if (utils::Compare(omega_tOverOmega_c_22, 0.01) > 0) {
        nuTidal22 = vlOver2 / std::sqrt(omega_tOverOmega_c_22);    
    }
    double k22Equilibrium    = k2_prefactor * nuTidal22 * omega_t_22;
    if (std::isnan(k22Equilibrium)) k22Equilibrium = 0.0;
    if (w22 < 0.0) k22Equilibrium = -std::abs(k22Equilibrium);


    // (l=2, n=3, m=2), Viscous dissipation, convective envelope
    double omega_t_32              = std::abs(w32);                                               
    double omega_t_over_omega_c_32 = omega_t_32 / omegaConv;
    double nuTidal32               = vl_5;
    if (utils::Compare(omega_t_over_omega_c_32, 5.0) > 0) {             
        nuTidal32 = vl25OverRoot20 / omega_t_over_omega_c_32 / omega_t_over_omega_c_32;
    }
    else if (utils::Compare(omega_t_over_omega_c_32, 0.01) > 0) {
        nuTidal32 = vlOver2 / std::sqrt(omega_t_over_omega_c_32);    
    }
    double k32Equilibrium    = k2_prefactor * nuTidal32 * omega_t_32;
    if (std::isnan(k32Equilibrium)) k32Equilibrium = 0.0;
    if (w32 < 0.0) k32Equilibrium = -std::abs(k32Equilibrium);

    // return ImKnmEquilibrium
    return std::make_tuple(k10Equilibrium, k12Equilibrium, k22Equilibrium, k32Equilibrium);
}


/*
 * Calculate the l=2, (n,m) = [(1,0), (1,2), (2,2), (3,2)] imaginary components of the potential tidal Love number 
 * by combining Equilibrium and Dynamical tidal contributions.
 *
 * DBL_DBL_DBL_DBL CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)
 *
 * @param   [IN]    p_Omega                     Orbital angular frequency (1/yr)
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of binary (AU)
 * @param   [IN]    p_M2                        Mass of companion star (Msol)
 * @return                                      [(1,0), (1,2), (2,2), (3,2)] Imaginary components of the 
 *                                              potential tidal Love number (unitless)
 */
DBL_DBL_DBL_DBL BaseStar::CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const {
    
    double Imk10Dynamical, Imk12Dynamical, Imk22Dynamical, Imk32Dynamical;
    std::tie(Imk10Dynamical, Imk12Dynamical, Imk22Dynamical, Imk32Dynamical) = CalculateImKnmDynamical(p_Omega, p_SemiMajorAxis, p_M2);

    double Imk10Equilibrium, Imk12Equilibrium, Imk22Equilibrium, Imk32Equilibrium;
    std::tie(Imk10Equilibrium, Imk12Equilibrium, Imk22Equilibrium, Imk32Equilibrium) = CalculateImKnmEquilibrium(p_Omega, p_SemiMajorAxis, p_M2);
    
    // return combined ImKnm terms;
    return std::make_tuple(Imk10Dynamical + Imk10Equilibrium, Imk12Dynamical + Imk12Equilibrium, Imk22Dynamical + Imk22Equilibrium, Imk32Dynamical + Imk32Equilibrium);
}




// THIS SHOULD BE IN either BINARYCONSTITUENTSTAR or BASEBINARYSTAR - NOT HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/*
 * Calculate eccentric anomaly and true anomaly - uses kepler's equation
 *
 * Modifies class member variables m_SupernovaDetails.eccentricAnomaly and m_SupernovaDetails.trueAnomaly
 *
 *
 * void CalculateSNAnomalies(const double p_Eccentricity)
 *
 * @param   [IN]    p_Eccentricity              Eccentricity of the binary
 */
void BaseStar::CalculateSNAnomalies(const double p_Eccentricity) {

    ERROR error = ERROR::NONE;

    std::tie(error, m_SupernovaDetails.eccentricAnomaly, m_SupernovaDetails.trueAnomaly) = utils::SolveKeplersEquation(m_SupernovaDetails.meanAnomaly, p_Eccentricity);

         if (error == ERROR::NO_CONVERGENCE) { THROW_ERROR(error, "Solving Kepler's equation"); }       // no convergence - throw error   
    else if (error == ERROR::OUT_OF_BOUNDS ) { THROW_ERROR(error, "Eccentric anomaly"); }               // eccentric anomaly out of bounds - throw error

    return;
}