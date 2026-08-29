#include "BaseStar.h"
#include "BH.h"          // needed for BH::CalculateSNkickWeighting_Static (was previously transitive via BaseStar.h)

#include <Star.h>


std::unique_ptr<BaseStar> BaseStar::Clone(const std::unique_ptr<BaseStar> p_Star, const STELLAR_TYPE p_StellarType, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) const {
    return m_Parent->Clone(p_Star, p_StellarType, p_Persistence, p_Initialise);
}

std::unique_ptr<BaseStar> BaseStar::CloneAs(const STELLAR_TYPE p_StellarType, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise) const {
    return m_Parent->CloneAs(p_StellarType, p_Persistence, p_Initialise);
}

std::unique_ptr<BaseStar> BaseStar::MakeStar(const STELLAR_TYPE p_StellarType, const OBJECT_PERSISTENCE p_Persistence) const {
    return m_Parent->MakeStar(p_StellarType, p_Persistence);
}
    










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
        case ANY_STAR_PROPERTY::BINDING_ENERGY_NANJING:                             value = CalculateBindingEnergy(CalculateCELambda_Nanjing(Mass(), Radius(), CoreMass()));    break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE:                           value = CalculateBindingEnergy(CalculateLambdaLoveridge());     break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_LOVERIDGE_WINDS:                     value = CalculateBindingEnergy(CalculateLambdaLoveridge(Mass() - CoreMass(), true)); break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_KRUCKOW:                             value = CalculateBindingEnergy(CalculateCELambdaKruckow());     break;
        case ANY_STAR_PROPERTY::BINDING_ENERGY_CONVECTIVE_ENVELOPE:                 value = CalculateConvectiveEnvelopeBindingEnergy(CalculateConvectiveEnvelopeLambda_Picker2024(CalculateConvectiveEnvelopeMass(Mass(), CoreMass()))); break;
        case ANY_STAR_PROPERTY::CHEMICALLY_HOMOGENEOUS_MAIN_SEQUENCE:               value = CHE() && utils::IsOneOf(StellarType(), MAIN_SEQUENCE);                                              break;
        case ANY_STAR_PROPERTY::CO_CORE_MASS:                                       value = COCoreMass();                                           break;
        case ANY_STAR_PROPERTY::CO_CORE_MASS_AT_COMPACT_OBJECT_FORMATION:           value = SN_COCoreMassAtCOFormation();                           break;
        case ANY_STAR_PROPERTY::CONVECTIVE_ENV_MASS:                                std::tie(value, std::ignore) = CalculateConvectiveEnvelopeMass(Mass(), CoreMass());  break;
        case ANY_STAR_PROPERTY::CORE_MASS:                                          value = CoreMass();                                             break;
        case ANY_STAR_PROPERTY::CORE_MASS_AT_COMPACT_OBJECT_FORMATION:              value = SN_CoreMassAtCOFormation();                             break;
        case ANY_STAR_PROPERTY::CORE_RADIUS_AT_COMPACT_OBJECT_FORMATION:            value = SN_CoreRadiusAtCOFormation();                           break; 
        case ANY_STAR_PROPERTY::DRAWN_KICK_MAGNITUDE:                               value = SN_DrawnKickMagnitude();                                break;
        case ANY_STAR_PROPERTY::DOMINANT_MASS_LOSS_RATE:                            value = static_cast<int>(DominantMassLossType());               break; /* was _RATE <<<<<<*/
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
        case ANY_STAR_PROPERTY::HELIUM_ABUNDANCE_CORE:                              value = HeAbundanceCore();                                  break;
        case ANY_STAR_PROPERTY::HELIUM_ABUNDANCE_SURFACE:                           value = HeAbundanceSurface();                               break;
        case ANY_STAR_PROPERTY::HYDROGEN_ABUNDANCE_CORE:                            value = HAbundanceCore();                                break;
        case ANY_STAR_PROPERTY::HYDROGEN_ABUNDANCE_SURFACE:                         value = HAbundanceSurface();                             break;
        case ANY_STAR_PROPERTY::IS_HYDROGEN_POOR:                                   value = SN_IsHydrogenPoor();                                    break;
        case ANY_STAR_PROPERTY::ID:                                                 value = ObjectId();                                             break;
        case ANY_STAR_PROPERTY::INITIAL_HELIUM_ABUNDANCE:                           value = InitialHeAbundance();                      break;
        case ANY_STAR_PROPERTY::INITIAL_HYDROGEN_ABUNDANCE:                         value = InitialHAbundance();                    break;
        case ANY_STAR_PROPERTY::INITIAL_STELLAR_TYPE:                               value = StartingStellarType();                                   break;
        case ANY_STAR_PROPERTY::INITIAL_STELLAR_TYPE_NAME:                          value = STELLAR_TYPE_LABEL.at(StartingStellarType());            break;
        case ANY_STAR_PROPERTY::IS_AIC:                                             value = IsAIC();                                                break;
        case ANY_STAR_PROPERTY::IS_CCSN:                                            value = IsCCSN();                                               break;
        case ANY_STAR_PROPERTY::IS_HeSD:                                            value = IsHeSD();                                               break;
        case ANY_STAR_PROPERTY::IS_ECSN:                                            value = IsECSN();                                               break;
        case ANY_STAR_PROPERTY::IS_PISN:                                            value = IsPISN();                                               break;
        case ANY_STAR_PROPERTY::IS_PPISN:                                           value = IsPPISN();                                              break;
        case ANY_STAR_PROPERTY::IS_SNIA:                                            value = IsSNIA();                                               break;
        case ANY_STAR_PROPERTY::IS_USSN:                                            value = IsUSSN();                                               break;
        case ANY_STAR_PROPERTY::KICK_MAGNITUDE:                                     value = SN_KickMagnitude();                                     break;
        case ANY_STAR_PROPERTY::LAMBDA_CONVECTIVE_ENVELOPE:                         value = CalculateConvectiveEnvelopeLambda_Picker2024(CalculateConvectiveEnvelopeMass(Mass(), CoreMass())); break;
        case ANY_STAR_PROPERTY::LAMBDA_DEWI:                                        value = CalculateCELambda_Dewi();                                  break;
        case ANY_STAR_PROPERTY::LAMBDA_FIXED:                                       value = OPTIONS->CommonEnvelopeLambda();                        break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW:                                     value = CalculateCELambdaKruckow();                             break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_BOTTOM:                              value = CalculateCELambdaKruckow(Radius(), -1.0);               break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_MIDDLE:                              value = CalculateCELambdaKruckow(Radius(), -4.0 / 5.0);         break;
        case ANY_STAR_PROPERTY::LAMBDA_KRUCKOW_TOP:                                 value = CalculateCELambdaKruckow(Radius(), -2.0 / 3.0);         break;
        case ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE:                                   value = CalculateLambdaLoveridge(Mass() - CoreMass(), false);   break;
        case ANY_STAR_PROPERTY::LAMBDA_LOVERIDGE_WINDS:                             value = CalculateLambdaLoveridge(Mass() - CoreMass(), true);    break;
        case ANY_STAR_PROPERTY::LAMBDA_NANJING:                                     value = CalculateCELambda_Nanjing(Mass(), Radius(), CoreMass()) break;
        case ANY_STAR_PROPERTY::LBV_PHASE_FLAG:                                     value = LBV_PhaseFlag();                                        break;
        case ANY_STAR_PROPERTY::LUMINOSITY:                                         value = Luminosity();                                           break;
        case ANY_STAR_PROPERTY::MASS:                                               value = Mass();                                                 break;
        case ANY_STAR_PROPERTY::MASS_0:                                             value = Mass0();                                                break;
        case ANY_STAR_PROPERTY::MDOT:                                               value = Mdot();                                                 break;
        case ANY_STAR_PROPERTY::MEAN_ANOMALY:                                       value = SN_MeanAnomaly();                                       break;
        case ANY_STAR_PROPERTY::METALLICITY:                                        value = Metallicity();                                          break;
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
        case ANY_STAR_PROPERTY::THERMAL_TIMESCALE:                                  value = CalculateTimescale_Thermal();                           break;
        case ANY_STAR_PROPERTY::TIME:                                               value = Time();                                                 break;
        case ANY_STAR_PROPERTY::TIMESCALE_MS:                                       value = Timescale(TS::MS);                                      break;
        case ANY_STAR_PROPERTY::TOTAL_MASS_AT_COMPACT_OBJECT_FORMATION:             value = SN_TotalMassAtCOFormation();                            break;
        case ANY_STAR_PROPERTY::TOTAL_RADIUS_AT_COMPACT_OBJECT_FORMATION:           value = SN_TotalRadiusAtCOFormation();                          break;
        case ANY_STAR_PROPERTY::TRUE_ANOMALY:                                       value = SN_TrueAnomaly();                                       break;
        case ANY_STAR_PROPERTY::TZAMS:                                              value = TZAMS() * TSOL;                                         break;
        case ANY_STAR_PROPERTY::VELOCITY_X:                                         value = VelocityX();											break;
        case ANY_STAR_PROPERTY::VELOCITY_Y:                                         value = VelocityY();											break;
        case ANY_STAR_PROPERTY::VELOCITY_Z:                                         value = VelocityZ();											break;
        case ANY_STAR_PROPERTY::ZETA_HURLEY:                                        value = CalculateZetaAdiabatic_Hurley2002(Mass(), CoreMass()); break;
        case ANY_STAR_PROPERTY::ZETA_HURLEY_HE:                                     value = CalculateZetaAdiabatic_Hurley2002(Mass(), HeCoreMass()); break;
        case ANY_STAR_PROPERTY::ZETA_SOBERMAN:                                      value = CalculateZetaAdiabatic_Soberman1997(Mass(), CoreMass()); break;
        case ANY_STAR_PROPERTY::ZETA_SOBERMAN_HE:                                   value = CalculateZetaAdiabatic_Soberman1997(Mass(), HeCoreMass()); break;

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
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      ENERGY                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               MASS LOSS / ACCRETION                               //
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
 * variables, some option values, and whatever else is needed here and pass them to functions
 * called from here, so that this function is telling the called functions what to do, rather
 * than the called functions interrogating state and (some) options to determine what they should
 * do.  We *do* allow called functions to interrogate option values, *except* where we might
 * want to choose whether the function uses an option value or a value we want to provide - these
 * instances are uncommon (the usual practice is to let the function determine what it uses to
 * calculate its returns value, but there are occasions where we want to override that).
 * The corallary is that most winds ML functions are GNU_CONST or GNU_PURE, but some are
 * COMPAS_PURE
 * 
 * 
 * The only
 * exception to this is the `CalculateMLRate_Merritt2025()` function which we pass variables to
 * as parameters, but we allow the function itself to determine which option values will be used
 * (because it is in some ways itself a driver function for winds mass loss rate calculations).
 * 
 * 
 * MassLossT CalculateMassLossRate() const
 *
 * @return                                      Tuple containing:
 *                                                  double  Mass loss rate (Msol yr^-1)
 *                                                          Will be clamped to [0.0, MAXIMUM_WIND_MASS_LOSS_RATE]
 *                                                  ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMassLossRate() const {

    // *Ilya* check this please - mStart is used in:
    //    - BaseStar::CalculateMLRate_Merritt2025()
    //    - BaseStar::CalculateMLRateRSG()
    //    - CalculateMLRateRSG_Beasor2020()
    //    - CalculateMassLossRateRSG_Decin2023()
    // All were originally ZAMS mass
    // We can easily get ZAMS mass if it exists (i.e. we started on MS), but we need to manage if it doesn't

    const double metallicity        = Metallicity();                            // Metallicity of the star
    const double mass               = Mass();                                   // Mass of the star
    const double mStart             = MassStart();                              // Mass of the star at the start of the simulation
    const double radius             = Radius();                                 // Radius of the star
    const double luminosity         = Luminosity();                             // Luminosity of the star
    const double temperature        = Temperature();                            // Temperature of the star
    const double HeAbundanceSurface = HeAbundanceSurface();                     // Surface helium abundance of the star

    const double perturb = Mu();                                     // Small envelope perturbation parameter

    double dMdt;
    ML_TYPE DominantMLtype;

    switch (OPTIONS->MassLossPrescription()) {                                  // Which mass loss prescription?

        case MASS_LOSS_PRESCRIPTION::BELCZYNSKI2010:                            // BELCZYNSKI2010
            std::tie(dMdt, dominantMLtype) = CalculateMLRate_Belczynski2010(
                mass,
                radius,
                luminosity,
                temperature,
                perturb,
                HeAbundanceSurface,
                OPTIONS->CoolWindMassLossMultiplier(),
                OPTIONS->LuminousBlueVariableFactor(),
                OPTIONS->LBVMassLossPrescription()
            );
            break;

        case MASS_LOSS_PRESCRIPTION::HURLEY:                                    // HURLEY
            std::tie(dMdt, dominantMLtype) = CalculateMLRate_Hurley2000(mass, radius, luminosity, perturb);

            const auto [dMdtLBV, dominantMLtypeLBV] = CalculateMLRateLBV(radius, luminosity, OPTIONS->LuminousBlueVariableFactor(), LBV_ML_PRESCRIPTION::HURLEY_ADD);

            if (dMdtLBV > dMdt) dominantMLtype = dominantMLtypeLBV;             // Dominant ML type
            dMdt += dMdtLBV;                                                    // Sum rates
            break;

        case MASS_LOSS_PRESCRIPTION::MERRITT2025:                               // MERRITT2025
            std::tie(dMdt, dominantMLtype) = CalculateMLRate_Merritt2025(
                metallicity,
                mass,
                radius,
                luminosity,
                temperature,
                HeAbundanceSurface,
                perturb,
                OPTIONS->CoolWindMassLossMultiplier(),
                OPTIONS->LuminousBlueVariableFactor(),
                mStart
            );
            break;

        case MASS_LOSS_PRESCRIPTION::ZERO:                                      // ZERO
            dMdt = 0.0;                                                         // No mass loss
            dominantMLtype = ML_TYPE::NONE;
            break;

        default:                                                                // Unexpected prescription
            // The only way this can happen is if someone added a MASS_LOSS_PRESCRIPTION and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a prescription this code doesn't account for, and that
            // should be flagged as an error and result in termination of the evolution of the star or
            // binary.
            // The correct fix for this is to add code for the missing prescription or, if the missing
            // prescription is superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNEXPECTED_MASS_LOSS_PRESCRIPTION);              // Throw error
    }

    // Apply overall wind mass loss multiplier and clamp winds to [0.0, MAXIMUM_WIND_MASS_LOSS_RATE]
    // to avoid convergence issues (maximum is typically 0.1 solar masses per year)    
    return std::make_tuple(std::clamp(dMdt * OPTIONS->OverallWindMassLossMultiplier(), 0.0, MAXIMUM_WIND_MASS_LOSS_RATE), dominantMLtype);
}


/*
 * CalculateMLRateLBV
 *
 * @brief
 * Calculate the LBV-like mass loss rate, and the dominant mass loss type, for stars
 * beyond the Humphreys-Davidson limit (Humphreys & Davidson 1994), based on the
 * LBV_ML_PRESCRIPTION passed in p_LBVprescription. 
 * 
 * The LBV mass loss prescription to be used is passed as a parameter, rather than this
 * function using any option value specified by the user, so that the caller can choose
 * the prescription to be used.
 *  
 *  
 * MassLossT CalculateMLRateLBV(const double p_Radius, const double p_Luminosity, const double p_LBVFactor, const LBV_ML_PRESCRIPTION p_MLPrescription) const
 *
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_LBVFactor                     LBV factor (see OPTIONS->LuminousBlueVariableFactor())
 * @param       p_MLPrescription                LBV mass loss prescription to use
 * @return                                      Tuple containing:
 *                                                  double  LBV-like mass loss rate (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMLRateLBV(const double p_Radius, const double p_Luminosity, const double p_LBVFactor, const LBV_ML_PRESCRIPTION p_MLPrescription) const {

    double dMdt = 0.0;                                                          // Default mass loss rate
    ML_TYPE DominantMLtype = ML_TYPE::NONE;                                     // Default dominant mass loss type

    const double HDlimitfactor = p_Radius * std::sqrt(p_Luminosity) * 1.0E-5;   // Factor by which the star is above the HD limit
    if (p_Luminosity > LBV_LUMINOSITY_LIMIT_STARTRACK && HDlimitfactor > 1.0) { // LBV?
        
        switch (p_MLPrescription) {                                             // Which LBV mass loss prescription?
            
            case LBV_ML_PRESCRIPTION::BELCZYNSKI:                               // BELCZYNSKI
                // Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit,
                // (Humphreys & Davidson 1994), per Belczynski et al. 2010, eq 8
                dMdt = p_LBVFactor * 1.0E-4;
                dominantMLtype = ML_TYPE::LBV;
                break;

            case LBV_ML_PRESCRIPTION::HURLEY_ADD:                               // HURLEY_ADD
            case LBV_ML_PRESCRIPTION::HURLEY:                                   // HURLEY
                // Calculate LBV-like mass loss rate for stars beyond the Humphreys-Davidson limit
                // (Humphreys & Davidson 1994), per Hurley et al. 2000, sec 7.1
                // (unlabelled equation a few equations after eq 106)
                dMdt = 0.1 * utils::intPow(HDlimitfactor - 1.0, 3) * ((p_Luminosity / 6.0E5) - 1.0);
                dominantMLtype = ML_TYPE::LBV;
                break;

            case LBV_ML_PRESCRIPTION::ZERO:                                     // ZERO
                dMdt = 0.0;                                                     // No mass loss
                dominantMLtype = ML_TYPE::NONE;                
            break;

            default:                                                            // Unexpected prescription
                // The only way this can happen is if the LBV_ML_PRESCRIPTION passed to this function is
                // not accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the code passed a prescription that this function doesn't account
                // for, and that should be flagged as an error and result in termination of the evolution
                // of the star or binary.
                // The correct fix for this is to add code to this function for the missing prescription,
                // or fix the calling code to pass a prescription that is handled by this function.
                THROW_ERROR(ERROR::UNEXPECTED_LBV_ML_PRESCRIPTION);             // Throw error
        }
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRateOB
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for main sequence stars,
 * based on the OB_ML_PRESCRIPTION. 
 * 
 * 
 * MassLossT CalculateMLRateOB(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  double  OB mass loss rate (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMLRateOB(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    double dMdt;                      
    ML_TYPE DominantMLtype;
    
    switch (OPTIONS->OBMassLossPrescription()) {                // Which OB mass loss prescription?

        case OB_ML_PRESCRIPTION::BJORKLUND2022:                  // BJORKLUND2022
            std::tie(dMdt, dominantMLtype) = CalculateMLRateOB_Bjorklund2022(p_Metallicity, p_Mass, p_Luminosity, p_Temperature);
            break;

        case OB_ML_PRESCRIPTION::KRTICKA2018:                    // KRTICKA2018
            std::tie(dMdt, dominantMLtype) = CalculateMLRateOB_Krticka2018(p_Metallicity, p_Luminosity);
            break;

        case OB_ML_PRESCRIPTION::VINK2001:                       // VINK2001
            std::tie(dMdt, dominantMLtype) = CalculateMLRateOB_Vink2001(p_Metallicity, p_Mass, p_Luminosity, p_Temperature);
            break;

        case OB_ML_PRESCRIPTION::VINK2021:                       // VINK2021 (Vink & Sander 2021)
            std::tie(dMdt, dominantMLtype) = CalculateMLRateOB_VinkSander2021(p_Mass, p_Luminosity, p_Temperature);
            break;

        case OB_ML_PRESCRIPTION::ZERO:                           // ZERO
            dMdt = 0.0;                                          // No mass loss
            dominantMLtype = ML_TYPE::NONE;                
            break;

        default:                                                 // Unexpected prescription
            // The only way this can happen is if the OB_ML_PRESCRIPTION served to this function is
            // not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR(ERROR::UNEXPECTED_OB_ML_PRESCRIPTION);  // Throw error
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRateOB_Bjorklund2022
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
 * MassLossT CalculateMLRateOB_Bjorklund2022(const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  double  Mass loss rate for massive stars (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (will be ML_TYPE::OB)
 */
MassLossT BaseStar::CalculateMLRateOB_Bjorklund2022(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    const double gamma   = p_Luminosity * LSOLW / CalculateEddingtonLuminosity(p_Mass, 0.1);    // Bjorklund et al. 2022, para 3, assumes He abundance = 0.1
    const double logZ    = ZDEP->SigmaHurley(p_Metallicity) + 1.853871964321762;                // log10(Z / 0.014) = log10(Z) - log10(0.014) (Bjorklund et al. 2022 uses 0.014)
    const double logL    = std::log10(p_Luminosity / 1.0E6);
    const double logTeff = std::log10(p_Temperature * TSOL / 45000.0);           
    const double logMeff = std::log10(p_Mass * (1.0 - gamma) / 45.0);

    // Eq 7 in Bjorklund et al. 2022
    return std::make_tuple(PPOW(10.0, -5.52 + (2.39 * logL) + (-1.48 * logMeff) + (2.12 * logTeff) + ((0.75 - (1.87 * logTeff)) * logZ)), ML_TYPE::OB);
}


/*
 * CalculateMLRateOB_Vink2001
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
 * MassLossT CalculateMLRateOB_Vink2001(const double p_Metallicity, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  double  Mass loss rate for hot OB stars (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMLRateOB_Vink2001(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    double dMdt;
    ML_TYPE DominantMLtype;

    const double teff = p_Temperature * TSOL;                                                                                       // Kelvin

    if (teff >= VINK_MASS_LOSS_MINIMUM_TEMP && teff <= VINK_MASS_LOSS_BISTABILITY_TEMP) {
        const double zeta = ZDEP->ZetaAnders(p_Metallicity);
        const double v = 1.3 * PPOW(zeta, OPTIONS->ScaleTerminalWindVelocityWithMetallicityPower());                                // v_inf / v_esc, scaled with metallicity

        dMdt = PPOW(10.0, -6.688 + (2.210 * std::log10(p_Luminosity / 1.0E5)) - (1.339 * std::log10(p_Mass / 30.0)) - (1.601 * std::log10(v / 2.0)) + (0.85 * zeta) + (1.07 * std::log10(teff / 20000.0)));
        dominantMLtype = ML_TYPE::OB; 
    }
    else if (teff > VINK_MASS_LOSS_BISTABILITY_TEMP) {
        const double v1 = 2.6 * PPOW(ZDEP->ZetaAnders(p_Metallicity), OPTIONS->ScaleTerminalWindVelocityWithMetallicityPower());    // v_inf / v_esc, scaled with metallicity
        const double v2 = std::log10(teff / 40000.0);

        dMdt = PPOW(10.0, -6.697 + (2.194 * std::log10(p_Luminosity / 1.0E5)) - (1.313 * std::log10(p_Mass / 30.0)) - (1.226 * std::log10(v1 / 2.0)) + (0.85 * ZDEP->ZetaAnders(Metallicity())) + (0.933 * v2) - (10.92 * v2 * v2));
        dominantMLtype = ML_TYPE::OB; 

        SHOW_WARN_IF(teff > VINK_MASS_LOSS_MAXIMUM_TEMP, ERROR::HIGH_TEFF_WINDS);                                                   // Show warning if winds being used outside comfort zone
    }
    else {                                                                                                                          // Too cold to use winds
        dMdt = 0.0;                                                                                                                 // Turn winds off
        dominantMLtype = ML_TYPE::NONE;                                                                                             // Set dominant type (NONE)
        SHOW_WARN(ERROR::LOW_TEFF_WINDS, "Mass Loss Rate = 0.0");                                                                   // Show warning
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRateOB_VinkSander2021
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
 * MassLossT CalculateMLRateOB_VinkSander2021(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  double  Mass loss rate for hot OB stars (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMLRateOB_VinkSander2021(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    double dMdt;
    ML_TYPE DominantMLtype;

    const double teff = p_Temperature * TSOL;                                           // Kelvin

    if (teff < VINK_MASS_LOSS_MINIMUM_TEMP) {                                           // Temp below Vink minimum?
                                                                                        // Yes - too cold to use winds
        dMdt = 0.0;                                                                     // Turn winds off
        dominantMLtype = ML_TYPE::NONE;

        SHOW_WARN(ERROR::LOW_TEFF_WINDS, "Mass Loss Rate = 0.0");                       // Show warning
    }
    else {                                                                              // No - temp is at or above minimum
        dominantMLtype = ML_TYPE::OB;                                                   // OB winds

        const double logL5  = std::log10(p_Luminosity / 1.0E5);
        const double logM30 = std::log10(p_Mass / 30.0);

        constexpr double Zexp2001 = 0.85;                                               // Vink et al. 2001

        const double zeta  = ZDEP->ZetaAnders(p_Metallicity);
        const double gamma = EDDINGTON_PARAMETER_FACTOR * p_Luminosity / p_Mass;
        const double rho   = -14.94 + (3.1857 * gamma) + (Zexp2001 * zeta);
        const double T1    = ( 100.0 + (6.0 * rho) ) * 1000.0;                          // Bistability jump 1: typically around 20000.0, has similar behavior when fixed

        if (teff <= T1) {                                                               // Temp at or below jump 1?
                                                                                        // Yes
            constexpr double v = 0.7;                                                   // v_inf / v_esc

            dMdt = PPOW(10.0, -5.99 + (2.210 * logL5) - (1.339 * logM30) - (1.601 * std::log10(v / 2.0)) + (Zexp2001 * zeta) + (1.07 * std::log10(teff / 20000.0)));
        }
        else {                                                                          // Temp is above jump 1
            SHOW_WARN_IF(teff > VINK_MASS_LOSS_MAXIMUM_TEMP, ERROR::HIGH_TEFF_WINDS);   // Show warning if temp above Vink maximum

            const double T2    = ( 61.2 + (2.59 * rho) ) * 1000.0;                      // Bistability jump 2: typically around 25000.0, higher jump first as in Vink python recipe

            if (teff <= T2) {                                                           // Temp at or below jump 2?
                                                                                        // Yes
                constexpr double v = 1.3;                                               // v_inf / v_esc

                dMdt = PPOW(10.0, -6.688 + (2.210 * logL5) - (1.339 * logM30) - (1.601 * std::log10(v / 2.0)) + (Zexp2001 * zeta) + (1.07 * std::log10(teff / 20000.0)));
            }
            else {                                                                      // Temp is above jump 2                      
                constexpr double Zexp2021 = 0.42;                                       // Vink and Sander 2021
                constexpr double v        = 2.6;                                        // v_inf / v_esc
                const     double logT40   = std::log10(teff / 40000.0);

                dMdt = PPOW(10.0, -6.697 + (2.194 * logL5) - (1.313 * logM30) - (1.226 * std::log10(v / 2.0)) + (Zexp2021 * zeta) + (0.933 * logT40) - (10.92 * logT40 * logT40));
            }
        }
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRateRSG
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for RSG stars (Red Supergiants),
 * based on the RSG_ML_PRESCRIPTION passed in p_MassLossPrescription. 
 * 
 * The RSG mass loss prescription to be used is passed as a parameter, rather than this
 * function using any option value specified by the user, so that the caller can choose
 * the prescription to be used.
 * 
 * 
 * MassLossT CalculateMLRateRSG(
 *     const double p_Metallicity,
 *     const double p_Mass,
 *     const double p_Radius,
 *     const double p_Luminosity,
 *     const double p_Temperature,
 *     const double p_mStart
 * ) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @return                                      Tuple containing:
 *                                                  double  RSG mass loss rate (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMLRateRSG(
    const double p_Metallicity,
    const double p_Mass,
    const double p_Radius,
    const double p_Luminosity,
    const double p_Temperature,
    const double p_mStart
) const {

    double dMdt;                      
    ML_TYPE DominantMLtype;

    switch (OPTIONS->RSGMassLossPrescription()) {               // Which RSG mass loss prescription?

        case RSG_ML_PRESCRIPTION::BEASOR2020:                   // BEASOR2020
            std::tie(dMdt, dominantMLtype) = CalculateMLRateRSG_Beasor2020(p_mStart, p_Luminosity);
            break;

        case RSG_ML_PRESCRIPTION::DECIN2023:                    // DECIN2023
            std::tie(dMdt, dominantMLtype) = CalculateMLRateRSG_Decin2023(p_mStart, p_Luminosity);
            break;

        case RSG_ML_PRESCRIPTION::KEE2021:                      // KEE2021
            std::tie(dMdt, dominantMLtype) = CalculateMLRateRSG_Kee2021(p_Mass, p_Luminosity, p_Temperature);
            break;

        case RSG_ML_PRESCRIPTION::NJ90:                         // NJ90
            std::tie(dMdt, dominantMLtype) = CalculateMLRate_NieuwenhuijzenDeJager1990(p_Metallicity, p_Mass, p_Radius, p_Luminosity);
            break;

        case RSG_ML_PRESCRIPTION::VINKSABHAHIT2023:             // VINKSABHAHIT2023
            std::tie(dMdt, dominantMLtype) = CalculateMLRateRSG_VinkSabhahit2023(p_Mass, p_Luminosity);
            break;

        case RSG_ML_PRESCRIPTION::YANG2023:                     // YANG2023
            std::tie(dMdt, dominantMLtype) = CalculateMLRateRSG_Yang2023(p_Luminosity);
            break;    

        case RSG_ML_PRESCRIPTION::ZERO:                         // ZERO
            dMdt = 0.0;                                         // No mass loss
            dominantMLtype = ML_TYPE::NONE;
            break;

        default:                                                // Unexpected prescription
            // The only way this can happen is if the RSG_ML_PRESCRIPTION served to this function is
            // not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR(ERROR::UNKNOWN_RSG_ML_PRESCRIPTION);    // Throw error
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRateRSG_Kee2021
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for Red Supergiant (RSG) stars,
 * per Kee et al. 2021, eqs 5, 13, 14, 25.
 * 
 * See https://arxiv.org/pdf/2101.03070.pdf  
 * 
 * 
 * MassLossT CalculateMLRateRSG_Kee2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  double  RSG mass loss rate (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (will be ML_TYPE::RSG)
 */
MassLossT BaseStar::CalculateMLRateRSG_Kee2021(const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    constexpr double vTurb  = 1.5E4;                                                                        // Turbulent velocity, ms^-1, for a typical RSG
    constexpr double kBoltz = 1.38E-23;                                                                     // Boltzmann Constant in J K^-1
    constexpr double sigma  = 5.67E-8;                                                                      // Stefan Boltzmann constant W m^-2 K^-4
    constexpr double mH     = 1.67E-27;                                                                     // Mass of hydrogen in Kg
    constexpr double kappa  = 0.01 * OPACITY_CGS_TO_SI;                                                     // Given after eq 16 

    const double tEff  = TSOL * p_Temperature;                                                              // Kelvin
    
    const double rSI   = std::sqrt((p_Luminosity * LSOLW) / (4.0 * M_PI * sigma * utils::intPow(tEff, 4))); // SI units
    const double mSI   = p_Mass * MSOL_TO_KG;                                                               // SI units
    const double cs    = std::sqrt(kBoltz * tEff / mH);
    const double gamma = (kappa * p_Luminosity * LSOLW) / (4.0 * M_PI * G * C * mSI);
    const double vEsc  = std::sqrt(2.0 * G * (mSI) / (rSI));                                                // ms^-1, not vEsc,eff

    const double rPmod = G * (mSI) * (1.0 - gamma) / (2.0 * ((cs * cs) + (vTurb * vTurb)));                 // Modified parker radius, in m
    const double rho   = (4.0 / 3.0) * (rPmod / (kappa * (rSI) * (rSI))) * (std::exp(-(2.0 * rPmod / (rSI)) + (3.0 / 2.0))) / (1.0 - std::exp(-2.0 * rPmod / (rSI)));

    const double mDotAnalytical = 4.0 * M_PI * rho * std::sqrt(cs * cs + vTurb * vTurb) * rPmod * rPmod;    // kg s^-1

    return std::make_tuple(PPOW(((vTurb / 17000.0) / (vEsc / 60000.0)), 1.30) * mDotAnalytical * SECONDS_IN_YEAR / MSOL_TO_KG, ML_TYPE::RSG);
}   


/*
 * CalculateMLRateVMS
 *
 * Calculate mass loss rate, and the dominant mass loss type, for very massive main sequence
 * stars (> 100 Msol), based on the VMS_ML_PRESCRIPTION. 
 * 
 * 
 * MassLossT CalculateMLRateVMS(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  double  VMS mass loss rate (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMLRateVMS(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    double dMdt;                         
    ML_TYPE DominantMLtype;

    switch (OPTIONS->VMSMassLossPrescription()) {               // Which VMS mass loss prescription?

        case VMS_ML_PRESCRIPTION::BESTENLEHNER2020:             // BESTENLEHNER2020
            std::tie(dMdt, dominantMLtype) = CalculateMLRateVMS_Bestenlehner2020(p_Mass, p_Luminosity);
            break;
        
        case VMS_ML_PRESCRIPTION::SABHAHIT2023:                 // SABHAHIT2023
            std::tie(dMdt, dominantMLtype) = CalculateMLRateVMS_Sabhahit2023(p_Metallicity, p_Mass, p_Luminosity, p_Temperature);
            break;
        
        case VMS_ML_PRESCRIPTION::VINK2011:                     // VINK2011
            std::tie(dMdt, dominantMLtype) = CalculateMLRateVMS_Vink2011(p_Metallicity, p_Mass, p_Luminosity, p_Temperature);
            break;

        case VMS_ML_PRESCRIPTION::ZERO:                         // ZERO
            dMdt = 0.0;                                         // No mass loss
            dominantMLtype = ML_TYPE::NONE;
        break;

        default:                                                // Unexpected prescription
            // The only way this can happen is if the VMS_ML_PRESCRIPTION served to this function is
            // not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR(ERROR::UNEXPECTED_VMS_ML_PRESCRIPTION); // Throw error
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRateVMS_Sabhahit2023
 *
 * @brief
 * Calculate mass loss rate, and the dominant mass loss type, for Very Massive Stars (VMS),
 * per Sabhahit 2023.
 * 
 * See https://arxiv.org/pdf/2306.11785.pdf
 *
 * 
 * MassLossT CalculateMLRateVMS_Sabhahit2023(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  double  Mass loss rate for very massive stars (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (will be ML_TYPE::VMS or ML_TYPE::OB)
 */
MassLossT BaseStar::CalculateMLRateVMS_Sabhahit2023(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    double dMdt;                        
    ML_TYPE DominantMLtype;

    const double zeta   = ZDEP->SigmaHurley(p_Metallicity);
    const double mSwitch = PPOW(Metallicity(), -1.574) * 0.0615 + 18.10;    // Obtained from a powerlaw fit to Sabhahit 2023, table 2, given teff = 45kK
    const double lSwitch = PPOW(10.0, (-1.91 * zeta + 2.36));               // Loglinear fits to Sabhahit 2023, table 2 

    if ((p_Luminosity / p_Mass) > (lSwitch / mSwitch)) {                    // In the VMS regime according to Sabhahit+ 2023?
                                                                            // Yes
        dMdt = PPOW(10.0, (-1.86 * zeta - 8.90)) * PPOW(p_Luminosity / lSwitch , 4.77) * PPOW(p_Mass / mSwitch, -3.99);
        dominantMLtype = ML_TYPE::VMS;                                      // Dominant mass loss type is VMS

    }
    else {                                                                  // No, fall back to default OB mass loss prescription
        std::tie(dMdt, dominantMLtype) = CalculateMLRateOB(p_Metallicity, p_Mass, p_Luminosity, p_Temperature);
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRateVMS_Vink2011
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for Very Massive (OB) Stars (VMS),
 * per Vink 2011.
 * 
 * See https://arxiv.org/pdf/1105.0556.pdf
 *
 * 
 * MassLossT CalculateMLRateVMS_Vink2011(const double p_Metallicity, p_Mass, p_Luminosity, p_Temperature) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @return                                      Tuple containing:
 *                                                  double  Mass loss rate for very massive stars (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (will be ML_TYPE::VMS or ML_TYPE::OB)
 */
MassLossT BaseStar::CalculateMLRateVMS_Vink2011(const double p_Metallicity, const double p_Mass, const double p_Luminosity, const double p_Temperature) const {

    double dMdt;
    ML_TYPE DominantMLtype;

    // Start with rate for massive OB stars, per Vink et al. 2001
    std::tie(dMdt, dominantMLtype) = CalculateMLRateOB_Vink2001(p_Metallicity, p_Mass, p_Luminosity, p_Temperature);

    const double gamma = EDDINGTON_PARAMETER_FACTOR * p_Luminosity / p_Mass;    // Eddington parameter, independent of surface composition

    if (gamma > 0.5) {                                                          // Apply correction to high gamma only
        dMdt           = PPOW(10.0, (0.04468 + (0.3091 * gamma) + (0.2434 * gamma * gamma) + std::log10(dMdt)));
        dominantMLtype = ML_TYPE::VMS;
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRate_Belczynski2010
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, per Belczynski 2010
 * (as implemented in StarTrack - courtesy Chris Belczynski).
 *
 * If option `--scale-mass-loss-with-surface-helium-abundance` was specified,
 * the mass loss rate will be scaled with the surface helium abundance.
 * 
 * 
 * MassLossT CalculateMLRate_Belczynski2010(
 *     const double                     p_Metallicity,
 *     const double                     p_Mass,
 *     const double                     p_Radius,
 *     const double                     p_Luminosity,
 *     const double                     p_Temperature,
 *     const double                     p_Perturb,
 *     const double                     p_HeAbundanceSurface,
 *     const double                     p_CoolWinMultiplier,
 *     const double                     p_LBVFactor,
 *     const LBV_ML_PRESCRIPTION p_LBVMLPrescription
 * ) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_Perturb                       Small envelope perturbation parameter
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @param       p_CoolWindMultiplier            Cool wind ML multiplier
 * @param       p_LBVFactor                     LBV factor (see OPTIONS->LuminousBlueVariableFactor())
 * @param       p_MLPrescription                LBV ML prescription
 * @return                                      Tuple containing:
 *                                                  double  Mass loss rate (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMLRate_Belczynski2010(
    const double              p_Metallicity,
    const double              p_Mass,
    const double              p_Radius,
    const double              p_Luminosity,
    const double              p_Temperature,
    const double              p_Perturb,
    const double              p_HeAbundanceSurface,
    const double              p_CoolWinMultiplier,
    const double              p_LBVFactor,
    const LBV_ML_PRESCRIPTION p_MLPrescription
) const {

    double dMdt;  
    ML_TYPE DominantMLtype;
   
    // Start with LBV winds (can be, and is often, 0.0)
    std::tie(dMdt, dominantMLtype) = CalculateMLRateLBV(p_Radius, p_Luminosity, p_LBVFactor, p_MLPrescription);

    // Other winds - if not in LBV regime, or user specified HURLEY_ADD 
    if (dominantMLtype != ML_TYPE::LBV || p_MLPrescription == LBV_ML_PRESCRIPTION::HURLEY_ADD) {

        double dMdtOther;
        ML_TYPE DominantMLtypeOther;

        if (p_Temperature * TSOL < VINK_MASS_LOSS_MINIMUM_TEMP) {                   // Cool star?
                                                                                    // Yes
            // Add Hurley et al. 2000 winds and apply cool wind mass loss multiplier
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLRate_Hurley2000(p_Mass, p_Radius, p_Luminosity, p_Perturb);
            dMdtOther *= p_CoolWinMultiplier;
        }
        else  {                                                                     // No - hot star
            // Add Vink et al. 2001 winds (ignoring bistability jump)
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLRateOB_Vink2001(p_Metallicity, p_Mass, p_Luminosity, p_Temperature);
        }

        if (dMdtOther > dMdt) dominantMLtype = dominantMLtypeOther;                 // Dominant ML type
        dMdt += dMdtOther;                                                          // Sum rates
    }

    // Note: BSE and StarTrack have a multiplier they apply here
    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRate_NieuwenhuijzenDeJager1990
 *
 * @brief
 * Calculate the mass loss rate, and the dominant mass loss type, for massive stars (L > 4000 Lsol),
 * per Nieuwenhuijzen & de Jager 1990, modified by a metallicity dependent factor (Kudritzki et al. 1989).
 *
 * See Hurley et al. 2000, just after eq 106
 *
 *
 * MassLossT CalculateMLRate_NieuwenhuijzenDeJager1990(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_Luminosity) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @return                                      Tuple containing:
 *                                                   double  RSG mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMLRate_NieuwenhuijzenDeJager1990(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_Luminosity) const {

    double dMdt = 0.0;                                                          // Default mass loss rate
    ML_TYPE DominantMLtype = ML_TYPE::NONE;                                     // Default dominant mass loss type

    // Mass loss only if star's luminosity is above minimum for Nieuwenhuijzen & de Jager 1990
    if (p_Luminosity > NJ_MINIMUM_LUMINOSITY) {
        dominantMLtype = ML_TYPE::GB;

        const double taper = std::min(1.0, (p_Luminosity - 4000.0) / 500.0);    // Smooth taper between no mass loss and mass loss
        // Calculate p_Metallicity / ZSOL_HURLEY dierctly rather than use ZDEP->ZScaledHurley()    
        dMdt = std::sqrt(p_Metallicity / ZSOL_HURLEY) * taper * 9.6E-15 * PPOW(p_Radius, 0.81) * PPOW(p_Luminosity, 1.24) * PPOW(p_Mass, 0.16);
    }
    
    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRateWR_SanderVink2020
 *
 * @brief
 * Calculate the mass loss rate for Wolf-Rayet stars, per Sander & Vink 2020, eq 13
 * https://arxiv.org/abs/2009.01849
 * 
 * 
 * static MassLossT CalculateMLRateWR_SanderVink2020(const double p_Metallicity, const double p_Luminosity, const double p_Perturb) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Perturb                       Small envelope perturbation parameter
 * @return                                      Tuple containing:
 *                                                   double  WR mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMLRateWR_SanderVink2020(const double p_Metallicity, const double p_Luminosity, const double p_Perturb) const {

    double dMdt = 0.0;                                                                                      // Default mass loss rate                          
    ML_TYPE DominantMLtype = ML_TYPE::NONE;                                                                 // Default dominant mass loss type

    const double zeta = ZDEP->ZetaAnders(p_Metallicity);

    if (p_Perturb < 1.0) {                                                                                  // Small envelope?
                                                                                                            // Yes
        const double logL      = std::log10(p_Luminosity);
        const double logL0     = -0.87 * zeta + 5.06;                                                       // Sander & Vink 2020, eq 19
        const double alpha     =  0.32 * zeta + 1.4;                                                        // ibid., eq 18
        const double logMdot10 = -0.75 * zeta - 4.06;                                                       // ibid., eq 20

        if (logL0 <= logL) {                                                                                // No mass loss for L < L0     
            dMdt = PPOW(10.0, alpha * std::log10(logL - logL0) + 0.75 * (logL - logL0 - 1.0) + logMdot10);  // ibid., eq 13
            dominantMLtype = ML_TYPE::WR;
        }
    }

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRate_Merritt2025
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
 *     const double p_mStart
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
 * @param       p_mStart                        Mass of the star at the start of the simulation (first state) (Msol)
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT BaseStar::CalculateMLRate_Merritt2025(
    const double p_Metallicity,
    const double p_Mass,
    const double p_Radius,
    const double p_Luminosity,
    const double p_Temperature,
    const double p_HeAbundanceSurface,
    const double p_Perturb,
    const double p_CoolWindMultiplier,
    const double p_LBVFactor,
    const double p_mStart
) const {

    double dMdt;
    ML_TYPE DominantMLtype;

    LBV_ML_PRESCRIPTION LBV_MLPrescription = OPTIONS->LBVMassLossPrescription();

    // Start with LBV winds (can be, and is often, 0.0)
    std::tie(dMdt, dominantMLtype) = CalculateMLRateLBV(p_Radius, p_Luminosity, p_LBVFactor, LBV_MLPrescription);

    // Calculate other winds rate if necessary.
    // We add other winds to the LBV winds if the use specified LBV_ML_PRESCRIPTION::HURLEY_ADD,
    // or if the winds are not in LBV regime
    if (LBV_MLPrescription == LBV_ML_PRESCRIPTION::HURLEY_ADD || dominantMLtype != ML_TYPE::LBV) { 

        double teff = p_Temperature * TSOL;                                 // Kelvin

        double dMdtOther;
        ML_TYPE DominantMLtypeOther;

        // RSG winds regime, massive star, and core helium burning giant (CHeB, FGB, EAGB, TPAGB), or HG?
        if (teff < RSG_MAXIMUM_TEMP && p_mStart >= MASSIVE_THRESHOLD && 
           (utils::IsOneOf(StellarType(), GIANTS) || utils::IsOneOf(StellarType(), {STELLAR_TYPE::HERTZSPRUNG_GAP}))) {
            // Yes - RSG mass loss rate
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLRateRSG(p_Mass, p_Radius, p_Luminosity, p_Temperature, p_mStart, OPTIONS->RSGMassLossPrescription());
        }

        // No - cool star?
        else if (teff < VINK_MASS_LOSS_MINIMUM_TEMP) {                      // Yes - HURLEY mass loss rate
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLRate_Hurley2000(p_Mass, p_Radius, p_Luminosity, p_Perturb);
            dMdt *= p_CoolWindMultiplier;                                   // Apply cool wind mass loss multiplier
        }

                                                                            // No - VMS winds regime?
        else if (p_Mass >= VMS_MASS_THRESHOLD) {                            // Yes - VMS mass loss rate
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLRateVMS(p_Metallicity, p_Mass, p_Luminosity, p_Temperature);
        }

        // Otherwise...
        else {                                                              // OB mass loss rate
            std::tie(dMdtOther, dominantMLtypeOther) = CalculateMLRateOB(p_Mass, p_Luminosity, p_Temperature, OPTIONS->OBMassLossPrescription());
        }

        if (dMdtOther > dMdt) dominantMLtype = dominantMLtypeOther;         // Dominant ML type
        dMdt += dMdtOther;                                                  // Sum rates
    }

    return std::make_tuple(dMdt, dominantMLtype);
}





/*
 * Calculate values for mDot and mass assuming mass loss is applied
 *
 * Class member variable m_Mdot is updated directly by this function if required (see parameters)
 * Class member variable m_Mass is not updated directly by this function - the calculated mass is returned as the functional return
 *
 * - calculates mass loss
 * - calculates new mass loss rate (mDot) to match (possibly limited) mass loss
 * - calculates new mass (mass) based on (possibly limited) mass loss
 * - returns existing value for mass if mass loss not being used (program option)
 *
 *
 * double CalculateMassLossValues(double p_Dt, const bool p_UpdateMDot)
 *
 * @param   [IN]    p_Dt                        time step (Myr)
 * @param   [IN]    p_UpdateMDot                flag to indicate whether the class member variable m_Mdot should be updated (default is false)
 * @return                                      calculated mass (mSol)
 */
double BaseStar::CalculateMassLossValues(double p_Dt, const bool p_UpdateMDot) {

    double mass = m_Mass;

    if (OPTIONS->MassLossPrescription() != MASS_LOSS_PRESCRIPTION::ZERO) {      // mass loss enabled?
                                                                                // yes
        double mDot     = CalculateMassLossRate();                              // calculate mass loss rate
        double massLoss = max(0.0, mDot * p_Dt * 1.0E6);                        // calculate mass loss; mass loss rate given in Msol per year, times are in Myr so need to multiply by 10^6
        if (p_UpdateMDot) m_Mdot = mDot;                                        // update class member variable if necessary

        if (OPTIONS->CheckPhotonTiringLimit()) {
            double lim = m_Luminosity / (G_SOLAR_YEAR * m_Mass / m_Radius);     // calculate the photon tiring limit in Msol yr^-1 using Owocki & Gayley 1997, equation slightly clearer in Owocki+2004 Eq. 20
            massLoss   = std::min(massLoss, lim);                               // limit mass loss to the photon tiring limit
            if (p_UpdateMDot) m_Mdot = massLoss / p_Dt / 1.0E6;                 // update class member variable if necessary
        }

        mass -= massLoss;                                                       // new mass based on mass loss
    }

    return mass;
}


/*
 * ResolveMassLossHurley
 *
 * @brief
 * Calculate winds mass loss per Hurley et al., 2000
 *
 * - calculate mass loss rate
 * - calculates (and limits) mass loss
 * - resets mass loss rate (CurrentState().Mdot()) to match (possibly limited) mass loss
 * - calculates and sets new mass (Mass()) based on (possibly limited) mass loss
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

    STELLAR_TYPE nextStellarType = StellarType();                                                   // next stellar type - defaults to current

    if (OPTIONS->MassLossPrescription() != MASS_LOSS_PRESCRIPTION::ZERO) {                          // mass loss enabled for Hurley?
                                                                                                    // yes

        double Mdot;
        ML_TYPE DominantMLtype;
        std::tie(Mdot, dominantMLtype) = CalculateMassLossRate();                                   // mass loss rate (Msol yr^-1)

/// DO THIS ELSEWHERE        double angularMomentumChange = (2.0 / 3.0) * (mass - m_Mass) * m_Radius * RSOL_TO_AU * m_Radius * RSOL_TO_AU * Omega();
          
        //nextStellarType = EvolveOneTimestep(mass - m_Mass, 0.0, 0.0, false);    // OLD - was: apply mass loss, but don't age the star yet


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



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ROTATION                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


// CalculateRotationalVelocityOStar_Ramirez2013 - MOVED to utils::CalculateRotationalVelocityOStar_Ramirez2013 (utils.h, free function)



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
 * StellarSNDetailsT BaseStar::CalculateSNkickMagnitude(
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
 * @param       p_SNDetails                     Supernova details object
 * @return                                      Mutated supernova details object
 */
StellarSNDetailsT BaseStar::CalculateSNkickMagnitude(
    const STELLAR_TYPE       p_RemnantType,
    const double             p_Mass,
    const double             p_EjectaMass,
    const double             p_RemnantMass,
    const StellarSNDetailsT& p_SNDetails
) const {

    if (!p_SNDetails.events.current.has_value()) {                                                                      // Current SN event?
        THROW_ERROR(ERROR::EXPECTED_SN_EVENT);                                                                          // No - throw error - this is a code defect
    }

    if (!utils::IsOneOf(p_RemnantType, {STELLAR_TYPE::NEUTRON_STAR, STELLAR_TYPE::BLACK_HOLE})) {                       // NS or BS expected?
        THROW_ERROR(ERROR::UNEXPECTED_STELLAR_TYPE);                                                                    // No - throw error - this is a code defect
    }


	double kick = 0.0;                                                                                                  // Default kick magnitude
    ERROR error = ERROR::NONE;                                                                                          // No error presumed

    StellarKickParmsT kickParameters = p_SNDetails.kickParameters;                  // Struct guaranteed to exist - members may not


    if (kickParameters.userMagnitude.has_value() && !kickParameters.userMagnitudeRandom.has_value()) {  // User supplied kick magnitude, but not kick random number?
        kick = userKickParams.magnitude.value();                                                        // Yes - use user supplied kick magnitude 
    }
    else {                                                                                              // No
        double sigma;
        switch (p_SNDetails.events.current) {                                       // What type of supernova event happening now?

		    case SN_EVENT::AIC:                                                     // AIC 
		    case SN_EVENT::SNIA:                                                    // SNIA 
		    case SN_EVENT::HeSD:                                                    // HeSD 
			    sigma = 0.0;                                                        // No kick
                break;

		    case SN_EVENT::PPISN:                                                   // PPISN
                if (!OPTIONS->NatalKickForPPISN()) {                                // User specified natal kicks for PPISN?
                    sigma = 0.0;                                                    // No - no natal kicks for PPISN
                    break;                                                          // Exit switch for PPISN and no kick
                }
                [[fallthrough]];                                                    // Otherwise fall through to CCSN

            case SN_EVENT::CCSN:                                                    // CCSN (PPISN with kicks falls through)                                    
                // Draw a random kick magnitude from the user-specified distribution
                // based on the expected stellar type of the remnant (NS or BH)
                sigma = p_RemnantType == STELLAR_TYPE::NEUTRON_STAR                 // Which remnant?
                            ? OPTIONS->KickMagnitudeDistributionSigmaCCSN_NS();     // Neutron star
                            : OPTIONS->KickMagnitudeDistributionSigmaCCSN_BH();     // Black hole
                break;

		    case SN_EVENT::ECSN:                                                    // ECSN may have a separate kick prescription
			    sigma = OPTIONS->KickMagnitudeDistributionSigmaForECSN();
                break;

		    case SN_EVENT::USSN:                                                    // USSN may have a separate kick prescription
			    sigma = OPTIONS->KickMagnitudeDistributionSigmaForUSSN();
                break;

            default:                                                                // Unexpected SN event type
                // the only way this can happen is if the SN_EVENT passed to this function is not accounted
                // for in this code.  We should not default here, with or without a warning.
                // We are here because the code passed an event type that this function doesn't account for,
                // and that should be flagged as an error and result in termination of the evolution of the
                // star or binary.
                // The correct fix for this is to add code to this function for the missing event type,
                // or fix the calling code to pass an event type that is handled by this function.

                error = ERROR::UNEXPECTED_SN_EVENT;                                 // Set error value (throw later)
	    }
    
	    if (error == ERROR::NONE) {                                                 // Check for errors
            // Sample the kick magnitude
            // kickParams.magnitudeRandom will have a valid value - either the user-specified value or the default value
            kick = DrawSNkickMagnitude(p_SNDetails.preSN.COCoreMass, p_EjectaMass, p_RemnantMass, sigma, kickParams.magnitudeRandom, p_SNdetails.events.current);
        }
    }
	if (error != ERROR::NONE) THROW_ERROR(error);                                   // Throw error if one occurred

    // No errors - continue
    StellarSupernovaDetailsT SNDetails = p_SNDetails;                               // Copy given SN details

    if (!SNdetails.currentKickParams.has_value()) SNdetails.currentKickParams = StellarKickParmsT{}; // just in case
    SNdetails.currentKickParams->magnitudeDrawn = kick;                             // set drawn kick magnitude

    // if there was a core-collapse supernova event this timestep, and the expected remnant is a black hole,
    // re-weight the (drawn) kick by the mass of the remnant according to the user-specified black hole kicks
    // option, if relevant (default is no reweighting)
    if ((p_SNDetails.events.has_value() ? p_SNDetails.events->current : SN_EVENT::NONE) == SN_EVENT::CCSN && p_RemnantType == STELLAR_TYPE::BLACK_HOLE) {
        kick *= BH::CalculateSNkickWeighting_Static(p_Mass, p_SNDetails.fallbackFraction.value_or(0.0));
    }
    else {                                                                          // otherwise
        SNDetails.fallbackFraction = OptDblT(0.0);                                  // set fallback fraction to zero
    }

    if (!SNDetails.currentKickParams.has_value()) SNDetails.currentKickParams = StellarKickParmsT{};
    SNDetails.currentKickParams->magnitude = kick;                                  // set updated kick magnitude

    // return SN details object by value - NRVO takes care of performance/efficiency
    return SNDetails;
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
double BaseStar::DrawSNkickMagnitude(
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

    return kick / OPTIONS->KickScalingFactor();
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
double BaseStar::DrawSNkickMagnitude_Muller2016(const double p_COCoreMass) const {

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
double BaseStar::DrawSNkickMagnitude_MullerMandel2020(const double p_COCoreMass, const double p_RemnantMass, const double p_Rand, const SN_EVENT p_SNevent) const {
    					
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
    
    STELLAR_TYPE nextStellarType = StellarType();                                       // next stellar type - defaults to current

    if (ShouldEvolveOnPhase()) {                                                        // should evolve timestep on phase?
                                                                                        // yes
        UpdateMainSequenceCoreMass(p_dt, -CurrentState().Mdot());                       // update core mass, relevant for MS stars

        m_InterimState.SetTau(CalculateTau());

        m_InterimState.SetCOCoreMass(CalculateCOCoreMass());
        m_InterimState.SetCoreMass(CalculateCoreMass());
        m_InterimState.SetHeCoreMass(CalculateHeCoreMass());

        m_InterimState.SetLuminosity(CalculateLuminosity());

        // Calculate abundances
        m_InterimState.SetHAbundanceCore(CalculateHAbundanceCore(Tau()));
        m_InterimState.SetHAbundanceSurface(CalculateHAbundanceSurface(Tau()));  
        m_InterimState.SetHeAbundanceCore(CalculateHeAbundanceCore(Tau()));
        m_InterimState.SetHeAbundanceSurface(CalculateHeAbundanceSurface(Tau()));
       
        double newRadius;
        std::tie(newRadius, nextStellarType) = CalculateRadiusAndStellarTypeOnPhase();  // radius and possibly new stellar type
        m_InterimState.SetRadius(newRadius);

        m_InterimState.SetMu(CalculateHurleyPerturbationMu());

        PerturbLuminosityAndRadiusOnPhase();

        m_InterimState.SetTemperature(CalculateTemperature());

        if (p_dt > 0.0) {
            STELLAR_TYPE thisStellarType = ResolveEnvelopeLoss();
            if (thisStellarType != StellarType()) {                                     // thisStellarType overrides stellarType
                nextStellarType = thisStellarType;
            }

            AdvanceAgeAndTime(p_dt);                                                    // advance age of star and simulation time
        }
    }

    return nextStellarType;
}













// AdvanceOneTimestep() notes:
//
// - This is the function that takes the current state of the star and evolves it one quantum of time (one timestep)
// - The assumption on entry to this function is that the star is fully described by the current state object
// - In general:
//    (a) the "firststate" state object describes the star at the start of evolution - this is not guaranteed to be at ZAMS. The "firststate" state object is guaranteed to exist.
//    (b) the "ZAMSstate" state object describes the star at ZAMS, but will only exist if evolution started at ZAMS.  Existence is testable.
//    (c) the timescales and GB parameters values specified in a state object are the values calculated before calculating dt and stellar attribute values (e.g. mass).
//    (d) per (c), the calculated (interim) values of timescales and GB parameters are used in the current timestep to calculate the new state of the star.
//    (e) dt (the timestep duration) is calculated after the new timescales and GB parameters are calculated.
//    (f) the attribute values in a state object describe the state of the star at the time/age/tau specified in that state object.
//    (g) the time/age/tau values specified in a state object are the values calculated after applying dt to the current state (i.e. they represent the end of the timestep that determined that state).
//    (h) the value of dt (the length of the timestep) in a state object is the value used to calculate state values in that state object (i.e. in that timestep).
//    (i) the value of dMdt (the rate of change of mass) in a state object is the value used to calculate state values in that state object (i.e. in that timestep).
//    (j) per (g), (h), and (i), the values of dt and dMdt are calculated using the current star state, and those calculated (interim) values are used in the current timestep to calculate the new state of the star.
// 
//    
//
//   



// AdvanceOneTimestep lives in BaseStar (not Star): it operates on the star's own
// state. The Star wrapper calls it and, after it returns, performs any stellar-type
// switch (SwitchTo) the returned type implies.


/*
 * AdvanceOneTimestep
 *
 * @brief
 * The transition function of the state machine: takes the star's frozen current
 * state and produces the next state, one quantum of time (one timestep) later.
 *
 *
 * The new state is built in m_InterimState and atomically committed to the state
 * history via Push() at the end. The next stellar type is returned; the caller
 * (the Star wrapper) switches type after the commit if it differs from the current
 * type.
 *
 * Massless-remnant and supernova handling occur OUTSIDE this function (before it is
 * called) - if we are here, the star is evolving normally.
 *
 *
 *
 * STELLAR_TYPE AdvanceOneTimestep()
 *
 * @return                                      Stellar type to which the star should evolve
 *                                              (== current type if no change). The new state has
 *                                              already been committed to the history.
 */
STELLAR_TYPE BaseStar::AdvanceOneTimestep(const double p_Metallicity, const double p_Dt, const double p_DeltaMass, const double p_DeltaMass0) {

    // Seed the interim state with a copy of the current state. Attributes we do not
    // explicitly recalculate this timestep therefore carry forward unchanged, and the
    // current GBParams / timescales are in place for the recalculation below.
    m_InterimState = CurrentState();


    // Set new metallicity first - all subsequent calculations use new metallicity
    m_InterimState.SetMetallicity(p_Metallicity);


    // GB parameters and timescales. CalculateGBParams()/CalculateTimescales() update
    // m_InterimState internally. These must be in place before CalculateTimestep()
    // reads them.
    m_InterimState.SetGBParams(CalculateGBParams());
    m_InterimState.SetTimescales(CalculateTimescales());


    // Timestep duration. If the caller supplied an explicit non-negative p_Dt it is
    // always honoured (zero is legitimate - it means "apply external mass change(s)   <<<<<<<<<<<<<<<<<<<< JR FIX THIS - 0 = legitimate for BSE ???
    // with no time advance", the BSE mass-only path). Only an internally-computed
    // dt <= 0 short-circuits the step (e.g. star already at end of life).
    double dt;
    if (p_Dt >= 0.0) {
        dt = p_Dt;                                                                                      // Explicit; zero OK
    }
    else {
        dt = CalculateTimestep();
        if (dt <= 0.0) {                                                                                // Computed step ran out
            m_InterimState.SetDt(dt);
            return StellarType(); // <<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS - AND COMMENT ABOVE
        }
    }
    m_InterimState.SetDt(dt);


    // Time / age / tau. The time and age deltas are just the timestep length.
    m_InterimState.SetTime(m_StateHistory.CurrentState.Time() + dt);
    m_InterimState.SetAge(m_StateHistory.CurrentState.Age() + dt);
    m_InterimState.SetTau(CalculateTau()); // JR FIX THIS: DONE


    // Unpack current-state inputs once (cheaper than repeated getter calls, and
    // makes explicit that these are the frozen start-of-timestep values). <<<<<<<<<<<<<<<< JR FIX THIS <<<<<<<<<<<<<<
    const double age    = Age();
    const double time   = Time();
    const double mass   = Mass();
    const double radius = Radius();
    const double lum    = Luminosity();


    // Mass change for this timestep. If the caller supplied an external delta
    // (BSE mass transfer), use it and skip wind mass loss. Otherwise compute the
    // wind rate and clamp to the photon-tiring limit (effectively infinite unless
    // --check-photon-tiring-limit is set).
    double deltaMass = 0.0;
    double dMdt      = 0.0;
    ML_TYPE DominantMLtype = ML_TYPE::NONE;

    if (CDOUBLE(p_DeltaMass) == 0.0) {                                          // Externally supplied (e.g. BSE mass transfer)
                                                                                // No
        if (OPTIONS->MassLossPrescription() != MASS_LOSS_PRESCRIPTION::ZERO) {  // Mass loss enabled?
                                                                                // Yes
            // Calculate maximum mass loss - check photon tiring limit if required
            const double maxMassLoss = OPTIONS->CheckPhotonTiringLimit() ? lum / (G_SOLAR_YEAR * mass / radius) : DBL_MAX;

            // Calculate mass loss rate
            std::tie(dMdt, dominantMLtype) = CalculateMassLossRate();           // Mass loss rate, dominant mass-loss type
            deltaMass = std::min(dMdt * dt * 1.0E6, maxMassLoss);               // Limit mass loss to maximum
            dMdt      = deltaMass / dt / 1.0E6;                                 // Actual mass loss rate
        }
    }
    else {                                                                      // Yes - externally supplied
        deltaMass = p_DeltaMass;                                                // Use it
    }

// <<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS - check that deltaMass is negative for winds mass loss (see std::min() above)!!
    m_InterimState.SetDMdt(dMdt);                                               // Set mass loss rate
    m_InterimState.SetDominantMLType(dominantMLtype);                           // Set dominant mass-loss type
// <<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS - if (deltaMass > 0.0) ??????????????????????????????????  shouldn't deltamass be < 0.0 for wind mass loss?????
    if (deltaMass > 0.0) m_InterimState.SetMass(mass + deltaMass);              // Set new mass


    // Effective initial mass (Mass0, Hurley et al. 2000 section 7)
    m_InterimState.SetMass0((CDOUBLE(p_DeltaMass0) == 0.0) ? CalculateMass0() : Mass0() + p_DeltaMass0);


    // Core masses
    m_InterimState.SetCoreMass(CalculateCoreMass());
    m_InterimState.SetCOCoreMass(CalculateCOCoreMass());
    m_InterimState.SetHeCoreMass(CalculateHeCoreMass());

    m_InterimState.SetCNOProcessedCoreMass(b.cnoProcessedCoreMass);  // <<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    // Abundances
    m_InterimState.SetHAbundanceCore(CalculateHAbundanceCore());
    m_InterimState.SetHAbundanceSurface(CalculateHAbundanceSurface()); 
    m_InterimState.SetHeAbundanceCore(CalculateHeAbundanceCore());
    m_InterimState.SetHeAbundanceCoreOut(m_InterimState.); // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    m_InterimState.SetHeAbundanceSurface(CalculateHeAbundanceSurface());


    // Radius.
    //
    // For the Hurley evolutionary model (EVOLUTION_MODE == SSE_HURLEY or BSE_HURLEY):
    //
    // In the Hurley model, HeHG and HeGB may change stellar type in response to a change
    // in radius (Hurley eqs 85-88, the min(R1, R2) construct).
    //
    // R1 is calculated per Hurley et al., 2000 eqs 86 & 87.
    // R2 mimics the Hayashi track for HeGB stars, per Hurley et al., 2000 eq 88.
    //
    // The CalculateRadius() function for HeHG and HeGB stars returns the minimum of
    // R1 and R2 (per Hurley et al., 2000 eqs 85-88), where R2 is calculated by the 
    // CalculateRadius_Hurley2000_Hayashi() function, per Hurley et al., 2000 eq 88,
    // as implemented in the HeHG and HeGB classes.  For classes other than HeHG and
    // HeGB, CalculateRadius_Hurley2000_Hayashi() returns DBL_MAX (to facilitate the
    // comparions later in this function).
    //    
    // Since the (R1, R2) construct is a Hurley construct (it may show up in other
    // evolutionary models, and if it does we'll need to address it then), we only
    // change the stellar type expectation for HeHG and HeGB stars if we are evolving
    // using one of the Hurley evolution modes.

    m_InterimState.SetRadius(CalculateRadius());                // Radius

    m_InterimState.SetLuminosity(CalculateLuminosity());        // Luminosity

    // Abundances (see BRCEK for the core-mass-prescription interaction).   <<<<<<<<<<<<<<<<<<<<<< JR FIX THIS
    m_InterimState.SetHAbundanceCore(CalculateHAbundanceCore(Tau()));
    m_InterimState.SetHAbundanceSurface(CalculateHAbundanceSurface(Tau()));
    m_InterimState.SetHeAbundanceCore(CalculateHeAbundanceCore(Tau()));
    m_InterimState.SetHeAbundanceSurface(CalculateHeAbundanceSurface(Tau()));






    STELLAR_TYPE nextStellarType = StellarType();
    if (ShouldSkipPhase()) {
        nextStellarType = ResolveSkippedPhase();
    }
    else {
        nextStellarType = EvolveOnPhase(dt);
        if (nextStellarType == StellarType()) nextStellarType = ResolveEndOfPhase();                   // still on type - check for end of phase
    }
    m_InterimState.SetStellarType(nextStellarType);





    // -------------------------------------------------------------------------
    // Commit
    // -------------------------------------------------------------------------
    // Atomically commit the new state to the history. The previous state is now
    // accessible via PreviousState(); there is no separate m_XxxPrev bookkeeping.
    (void)m_StateHistory.Push(m_InterimState);

    return nextStellarType;                                                                            // caller switches type after commit if changed
}
    
    














/*
 * Evolve the star on it's current phase - take one timestep on the current phase
 *
 *
 * STELLAR_TYPE EvolveOnPhase(const double p_DeltaTime)
 *
 * @param   [IN]    p_DeltaTime                 Timestep in Myr
 * @return                                      Stellar Type to which star should evolve - unchanged if not moving off current phase
 */
STELLAR_TYPE BaseStar::EvolveOnPhase(const double p_DeltaTime) {

    STELLAR_TYPE stellarType = m_StellarType;

    if (ShouldEvolveOnPhase()) {                                                    // evolve timestep on phase
        
        UpdateMainSequenceCoreMass(p_DeltaTime, -m_Mdot);                           // update core mass, relevant for MS stars

        m_Tau        = CalculateTauOnPhase();

        m_COCoreMass = CalculateCOCoreMassOnPhase();
        m_CoreMass   = CalculateCoreMassOnPhase();
        m_HeCoreMass = CalculateHeCoreMassOnPhase();

        m_Luminosity = CalculateLuminosityOnPhase();

        // Calculate abundances
        m_HeliumAbundanceCore      = CalculateHeliumAbundanceCoreOnPhase();
        m_HeliumAbundanceSurface   = CalculateHeliumAbundanceSurfaceOnPhase();
        m_HydrogenAbundanceCore    = CalculateHydrogenAbundanceCoreOnPhase();
        m_HydrogenAbundanceSurface = CalculateHydrogenAbundanceSurfaceOnPhase();  
        
        std::tie(m_Radius, stellarType) = CalculateRadiusAndStellarTypeOnPhase();   // radius and possibly new stellar type

        m_Mu = CalculatePerturbationMuOnPhase();

        PerturbLuminosityAndRadiusOnPhase();

        m_Temperature = CalculateTemperatureOnPhase();

        if (p_DeltaTime > 0.0) {
            STELLAR_TYPE thisStellarType = ResolveEnvelopeLoss();                   // resolve envelope loss if it occurs - possibly new stellar type
            if (thisStellarType != m_StellarType) {                                 // thisStellarType overrides stellarType (from CalculateRadiusAndStellarTypeOnPhase())
                stellarType = thisStellarType;
            }
        }
    }

    return stellarType;
}


/*
 * Evolve the star onto the next phase if necessary - take one timestep at the end of the current phase
 *
 *
 * STELLAR_TYPE ResolveEndOfPhase(const bool p_ResolveEnvelopeLoss)
 *
 * @param   [IN]    p_ResolveEnvelopeLoss       Specifies whether envelope loss should be resolved here
 *                                              (optional, default = true)
 *                                              JR: this is a bit of a kludge to resolve problems introduced by modifying stellar attributes in
 *                                                  anticipation of switching stellar type, but using those attributes before the actual switch
 *                                                  to the new stellar type - we need to resolve those situations in the code.
 * @return                                      Stellar Type to which star should evolve - unchanged if not moving off current phase
 */
STELLAR_TYPE BaseStar::ResolveEndOfPhase(const bool p_ResolveEnvelopeLoss) {

    STELLAR_TYPE stellarType = m_StellarType;

    if (IsEndOfPhase()) {                                                       // end of phase

        if (p_ResolveEnvelopeLoss) stellarType = ResolveEnvelopeLoss();         // if required, resolve envelope loss if it occurs

        if (stellarType == m_StellarType) {                                     // staying on phase?
            
            m_Tau         = CalculateTauAtPhaseEnd();

            m_COCoreMass  = CalculateCOCoreMassAtPhaseEnd();
            m_CoreMass    = CalculateCoreMassAtPhaseEnd();
            m_HeCoreMass  = CalculateHeCoreMassAtPhaseEnd();

            m_Luminosity  = CalculateLuminosityAtPhaseEnd();

            m_Radius      = CalculateRadiusAtPhaseEnd();

            m_Mu          = CalculatePerturbationMuAtPhaseEnd();

            PerturbLuminosityAndRadiusAtPhaseEnd();

            m_Temperature = CalculateTemperatureAtPhaseEnd();

            stellarType   = EvolveToNextPhase();                                // determine the stellar type to which the star should evolve
        }
    }

    return stellarType;
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
 * The assumed ZAMS mass passed to this function should be the actual ZAMS mass of the star, unless
 * the user specified the `--common-envelope-lambda-nanjing-use-rejuvenated-mass` commandline option,
 * in which case the assumed mass should be the "rejuvenated" (or effective) mass of the star.  The
 * corollary is that this function should not be used if the ZAMS mass is not known, unless the user
 * specified the `--common-envelope-lambda-nanjing-use-rejuvenated-mass` commandline option.
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
 * not apply.
 * 
 * 
 * double CalculateCELambda_Nanjing(const double p_MZAMS, const double p_Radius, const double p_CoreMass) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Common envelope lambda parameter
 */ 
double BaseStar_Constituent::CalculateCELambda_Nanjing(const double p_Metallicity, const double p_MZAMS, const double p_Radius, const double p_CoreMass) const {

    double lambda = 0.0;                                                            // Default return value

    if (OPTIONS->CommonEnvelopeLambdaNanjingEnhanced()) {                           // Use enhanced Nanjing prescription?
                                                                                    // Yes
        // Set stellar population based on metallicity
        STELLAR_POPULATION pop = p_Metallicity < LAMBDA_NANJING_ZLIMIT ? STELLAR_POPULATION::POPULATION_II : STELLAR_POPULATION::POPULATION_I;

        if (OPTIONS->CommonEnvelopeLambdaNanjingInterpolateInMass()) {              // Interpolate across mass models?
            if (OPTIONS->CommonEnvelopeLambdaNanjingInterpolateInMetallicity()) {   // Yes - interpolate across metallicities?
                lambda = BaseStar::CalculateCELambda_Nanjing_MassAndZInterpolated(p_Metallicity, p_MZAMS, p_Radius, p_CoreMass); // Yes
            }
            else {                                                                  // No - not interpolating across metallicities
                lambda = BaseStar::CalculateCELambda_Nanjing_MassInterpolated(p_MZAMS, p_Radius, p_CoreMass, pop);
            }
        }
        else {                                                                      // No - not interpolating across mass models
            // calculate index into NANJING_MASSES (see constants.h) for p_Mass
            std::size_t massIndex;
                 if (p_MZAMS < NANJING_MASSES_MIDPOINTS[0])      massIndex = 0;     // Use lambda for the 1 Msun model
            else if (p_MZAMS >= NANJING_MASSES_MIDPOINTS.back()) massIndex = NANJING_MASSES.size() - 1; // Use lambda for the 100 Msun model
            else                                                 massIndex = utils::BinarySearch(NANJING_MASSES_MIDPOINTS, p_MZAMS)[1]; // Bin edge indices - use upper

            if (OPTIONS->CommonEnvelopeLambdaNanjingInterpolateInMetallicity()) {   // Interpolate across stellar populations?
                lambda = BaseStar::CalculateCELambda_Nanjing_ZInterpolated(p_MZAMS, p_Radius, p_CoreMass, massIndex); // Yes
            }
            else {                                                                  // No - not interpolating across stellar populations
                lambda = BaseStar::CalculateCELambda_Nanjing_Enhanced(p_MZAMS, p_Radius, p_CoreMass, massIndex, pop);
            }
        }
    }
    else {                                                                          // No - not enhanced Nanjing prescription: use StarTrack prescription
        lambda = CalculateCELambda_Nanjing_StarTrack(p_MZAMS, p_Radius, p_CoreMass);
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
 * double CalculateCELambda_Nanjing_MassAndZInterpolated(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_CoreMass) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Common envelope lambda parameter
 */ 
double BaseStar::CalculateCELambda_Nanjing_MassAndZInterpolated(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_CoreMass) const {

    double lambda = 0.0; // Default return value

    // Use lambda for pop. I or pop. II metallicity, otherwise interpolate
         if (p_Metallicity < LAMBDA_NANJING_POPII_Z) lambda = BaseStar::CalculateCELambda_Nanjing_MassInterpolated(p_Mass, p_Radius, p_CoreMass, STELLAR_POPULATION::POPULATION_II);
    else if (p_Metallicity > LAMBDA_NANJING_POPI_Z)  lambda = BaseStar::CalculateCELambda_Nanjing_MassInterpolated(p_Mass, p_Radius, p_CoreMass, STELLAR_POPULATION::POPULATION_I);
    else {
        // linear interpolation in logZ between pop. I and pop. II metallicities
        const double lambdaLow = BaseStar::CalculateCELambda_Nanjing_MassInterpolated(p_Mass, p_Radius, p_CoreMass, STELLAR_POPULATION::POPULATION_II);
        const double lambdaUp  = BaseStar::CalculateCELambda_Nanjing_MassInterpolated(p_Mass, p_Radius, p_CoreMass, STELLAR_POPULATION::POPULATION_I);

        lambda = lambdaLow + (ZDEP->SigmaHurley(p_Metallicity) - LAMBDA_NANJING_POPII_LOGZ) / (LAMBDA_NANJING_POPI_LOGZ - LAMBDA_NANJING_POPII_LOGZ) * (lambdaUp - lambdaLow);
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

    SizeTVectorT indices = utils::BinarySearch(NANJING_MASSES, p_Mass);                                 // find mass in NANJING_MASSES

    SizeT lower = indices[0];                                                                           // bin lower bound
    SizeT upper = indices[1];                                                                           // bin upper bound

    double lambda;
    if ((lower == SIZE_MAX) && (upper != SIZE_MAX)) {                                                   // mass below range calculated by Xu & Li 2010?
        lambda = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, 0, p_StellarPop);     // yes - return lambda for minimum mass
    }
    else if ((lower != SIZE_MAX) && (upper == SIZE_MAX)) {                                              // mass above range calculated by Xu & Li 2010?
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
 * double CalculateCELambda_Nanjing_ZInterpolated(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_CoreMass, const SizeT p_MassIndex) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_MassIndex                     Index of mass bin in NANJING_MASSES (see constants.h)
 * @return                                      Common envelope lambda parameter
 */ 
double BaseStar::CalculateCELambda_Nanjing_ZInterpolated(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_CoreMass, const SizeT p_MassIndex) const {
    
    double lambda = 0.0; // Default return value
    
    // Use lambda for pop. I or pop. II metallicity, otherwise interpolate
         if (p_Metallicity < LAMBDA_NANJING_POPII_Z) lambda = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, p_MassIndex, STELLAR_POPULATION::POPULATION_II);
    else if (p_Metallicity > LAMBDA_NANJING_POPI_Z)  lambda = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, p_MassIndex, STELLAR_POPULATION::POPULATION_I);
    else {
        // Linear interpolation in logZ between pop. I and pop. II metallicities
        const double lambdaLow = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, p_MassIndex, STELLAR_POPULATION::POPULATION_II);
        const double lambdaUp  = CalculateCELambda_Nanjing_Enhanced(p_Mass, p_Radius, p_CoreMass, p_MassIndex, STELLAR_POPULATION::POPULATION_I);

        lambda = lambdaLow + (ZDEP->SigmaHurley(p_Metallicity) - LAMBDA_NANJING_POPII_LOGZ) / (LAMBDA_NANJING_POPI_LOGZ - LAMBDA_NANJING_POPII_LOGZ) * (lambdaUp - lambdaLow);
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

    const double mRatio    = p_CoreMass / Mass();       // Soberman, Phinney, van den Heuvel 1997, eq 57
    const double oneMinusM = 1.0 - mRatio;

    // Soberman, Phinney, van den Heuvel 1997, eq 61
    return ((2.0 / 3.0) * mRatio / oneMinusM) - ((1.0 / 3.0) * (oneMinusM / (1.0 + (mRatio + mRatio)))) - (0.03 * mRatio) + (0.2 * mRatio / (1.0 + (1.0 / utils::intPow(oneMinusM, 6)))); 
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
            qCrit = CalculateCriticalMassRatio_Ge2020(p_MTefficiency);
            break;

        case QCRIT_PRESCRIPTION::CLAEYS:                                                // CLAEYS
            qCrit = CalculateCriticalMassRatio_Claeys2014(p_AccretorIsDegenerate);
            break;

        case QCRIT_PRESCRIPTION::HURLEY_HJELLMING_WEBBINK:                              // HURLEY_HJELLMING_WEBBINK
            qCrit = CalculateCriticalMassRatio_Hurley2002();
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
 * Dbl_DblT CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate)
 *
 * @param   [IN]    p_DonorMassRate             Mass transfer rate of the donor
 * @param   [IN]    p_AccretorMassRate          Thermal mass transfer rate of the accretor (this star)
 * @return                                      Tuple containing the Maximum Mass Acceptance Rate and the Accretion Efficiency Parameter
 */
Dbl_DblT BaseStar::CalculateMassAcceptanceRate(const double p_DonorMassRate, const double p_AccretorMassRate) {

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
            acceptanceRate = (Mass() - CoreMass()) / CalculateTimescale_Thermal(Mass(), p_Radius, Luminosity(), CoreMass());
            break;
        case MT_THERMALLY_LIMITED_VARIATION::C_FACTOR:
            acceptanceRate = CalculateMLRateThermal();
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
















// JR: why is this in basestar?  a single star shouldn't know about tides.... <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
 * Dbl_Dbl_Dbl_DblT CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)
 *
 * @param   [IN]    p_Omega                     Orbital angular frequency (1/yr)
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of binary (AU)
 * @param   [IN]    p_M2                        Mass of companion star (Msol)
 * @return                                      [(1,0), (1,2), (2,2), (3,2)] Imaginary components of the 
 *                                              potential tidal Love number, Dynamical tides only (unitless)
 */
Dbl_Dbl_Dbl_DblT BaseStar::CalculateImKnmDynamical(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const {
    
    double coreMass = CalculateConvectiveCoreMass();

    double envMass, envMassMax;
    std::tie(envMass, envMassMax) = CalculateConvectiveEnvelopeMass();
    
    double radIntershellMass = Mass() - coreMass - envMass;                                             // refers to the combined mass of non-convective layers

    double radiusAU              = Radius() * RSOL_TO_AU;
    double coreRadiusAU          = CalculateConvectiveCoreRadius() * RSOL_TO_AU;
    double convectiveEnvRadiusAU = CalculateConvectiveEnvelopeRadialExtent() * RSOL_TO_AU;
    double radiusIntershellAU    = radiusAU - convectiveEnvRadiusAU;                                    // Outer radial coordinate of radiative intershell, or inner radial coordinate of convective envelope

    // There should be no Dynamical tides if the entire star is convective, i.e. if there are no convective-radiative boundaries. 
    // If so, return 0.0 for all dynamical components of ImKnm.
    // This condition should be true for low-mass MS stars (<= 0.35 Msol) at ZAMS.
    if (utils::Compare(radIntershellMass/Mass(), TIDES_MINIMUM_FRACTIONAL_EXTENT) <= 0 || utils::Compare(radiusIntershellAU, coreRadiusAU) <= 0) {
        return std::make_tuple(0.0, 0.0, 0.0, 0.0);                           
    }
 
    double R_3            = radiusAU * radiusAU * radiusAU;
    double R3OverG_M      = (R_3 / G_AU_Msol_yr / Mass());
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
    if (utils::Compare(coreRadiusAU/radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0 && utils::Compare(coreMass/Mass(), TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0 && utils::Compare(convectiveEnvRadiusAU/radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) < 0 && utils::Compare(envMass/Mass(), TIDES_MINIMUM_FRACTIONAL_EXTENT) < 0) {                   
        constexpr double beta2Dynamical         = 1.0;
        constexpr double rhoFactorDynamcial     = 0.1;
        double coreRadiusOverRadius   = coreRadiusAU / radiusAU;
        double coreRadiusOverRadius_3 = coreRadiusOverRadius * coreRadiusOverRadius * coreRadiusOverRadius;
        double coreRadiusOverRadius_9 = coreRadiusOverRadius_3 * coreRadiusOverRadius_3 * coreRadiusOverRadius_3;
        double massOverCoreMass       = Mass() / coreMass;
        double E2Core                 = (2.0 / 3.0) * coreRadiusOverRadius_9 * massOverCoreMass * std::cbrt(massOverCoreMass) * beta2Dynamical * rhoFactorDynamcial;

        // (l=2, n=1, m=0), Gravity Wave dissipation from core boundary
        double s10     = w10 * sqrtR3OverG_M;
        double s10_4_3 = s10 * std::cbrt(s10);
        double s10_8_3 = s10_4_3 * s10_4_3;
        k10GravityCore = E2Core *  std::copysign(s10_8_3, w10);
        if (std::isnan(k10GravityCore)) k10GravityCore = 0.0;

        // (l=2, n=1, m=2), Gravity Wave dissipation from core boundary
        double s12     = w12 * sqrtR3OverG_M;
        double s12_4_3 = s12 * std::cbrt(s12);
        double s12_8_3 = s12_4_3 * s12_4_3;
        k12GravityCore = E2Core * std::copysign(s12_8_3, w12);
        if (std::isnan(k12GravityCore)) k12GravityCore = 0.0;

        // (l=2, n=2, m=2), Gravity Wave dissipation from core boundary
        double s22     = w22 * sqrtR3OverG_M;
        double s22_4_3 = s22 * std::cbrt(s22);
        double s22_8_3 = s22_4_3 * s22_4_3;
        k22GravityCore = E2Core * std::copysign(s22_8_3, w22);
        if (std::isnan(k22GravityCore)) k22GravityCore = 0.0;

        // (l=2, n=3, m=2), Gravity Wave dissipation from core boundary
        double s32     = w32 * sqrtR3OverG_M;
        double s32_4_3 = s32 * std::cbrt(s32);
        double s32_8_3 = s32_4_3 * s32_4_3;
        k32GravityCore = E2Core * std::copysign(s32_8_3, w32);
        if (std::isnan(k32GravityCore)) k32GravityCore = 0.0;    
    }

    double rint_3 = radiusIntershellAU * radiusIntershellAU * radiusIntershellAU;
    double rc_3   = coreRadiusAU * coreRadiusAU * coreRadiusAU;
    double gamma  = (envMass / (R_3 - rint_3)) / (radIntershellMass / (rint_3 - rc_3));

    // There is no GW or IW dissipation from the envelope boundary if no convective envelope
    if ((utils::Compare(convectiveEnvRadiusAU / radiusAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0) || (utils::Compare(envMass / Mass(), TIDES_MINIMUM_FRACTIONAL_EXTENT) > 0)) {    

        constexpr double dynPrefactor     = 3.207452512782476;                                                        // 3^(11/3) * Gamma(1/3)^2 / 40 PI
        constexpr double m_l_factor_22    = 0.091720201358184;                                                        // (l(l+1))^{-4/3}, assuming l=2

        double cbrtdNdlnr       = std::cbrt(G_AU_Msol_yr * radIntershellMass / radiusIntershellAU / radiusIntershellAU / (radiusAU - radiusIntershellAU));
        
        double alpha            = radiusIntershellAU / radiusAU;
        double oneMinusAlpha    = 1.0 - alpha;
        double beta             = radIntershellMass / Mass();

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
            double Epsilon       = alpha_11 * envMass / Mass() * oneMinusGamma_2 * alpha_2_3Minus_1 * alpha_2_3Minus_1 / beta_2 / oneMinusAlpha_3 / oneMinusAlpha_2;
            double R3_Epsilon_dNdlnr_factor = R3OverG_M * Epsilon / cbrtdNdlnr;
            double E2Envelope    = dynPrefactor * m_l_factor_22 * R3_Epsilon_dNdlnr_factor;

            // (l=2, n=1, m=0), Gravity Wave dissipation from envelope boundary
            double w10_4_3       = w10 * std::cbrt(w10);
            double w10_8_3       = w10_4_3 * w10_4_3;
            k10GravityEnv        = E2Envelope * std::copysign(w10_8_3, w10);
            if (std::isnan(k10GravityEnv)) k10GravityEnv = 0.0;  

            // (l=2, n=1, m=2), Gravity Wave dissipation from envelope boundary
            double w12_4_3       = w12 * std::cbrt(w12);
            double w12_8_3       = w12_4_3 * w12_4_3;
            k12GravityEnv        = E2Envelope * std::copysign(w12_8_3, w12);
            if (std::isnan(k12GravityEnv)) k12GravityEnv = 0.0;  

            // (l=2, n=2, m=2), Gravity Wave dissipation from envelope boundary
            double w22_4_3       = w22 * std::cbrt(w22);
            double w22_8_3       = w22_4_3 * w22_4_3;
            k22GravityEnv        = E2Envelope * std::copysign(w22_8_3, w22);
            if (std::isnan(k22GravityEnv)) k22GravityEnv = 0.0;  

            // (l=2, n=3, m=2), Gravity Wave dissipation from envelope boundary
            double w32_4_3       = w32 * std::cbrt(w32);
            double w32_8_3       = w32_4_3 * w32_4_3;
            k32GravityEnv        = E2Envelope * std::copysign(w32_8_3, w32);
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
 * Dbl_Dbl_Dbl_DblT CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)
 *
 * @param   [IN]    p_Omega                     Orbital angular frequency (1/yr)
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of binary (AU)
 * @param   [IN]    p_M2                        Mass of companion star (Msol)
 * @return                                      [(1,0), (1,2), (2,2), (3,2)] Imaginary components of the 
 *                                              potential tidal Love number, Equilibrium tides only (unitless)
 */
Dbl_Dbl_Dbl_DblT BaseStar::CalculateImKnmEquilibrium(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const {

    // Viscous dissipation
    // No contribution from convective core; only convective envelope.

    double rEnvAU = CalculateConvectiveEnvelopeRadialExtent() * RSOL_TO_AU;
    double envMass, envMassMax;
    std::tie(envMass, envMassMax) = CalculateConvectiveEnvelopeMass();
    
    double rOutAU = Radius() * RSOL_TO_AU;                                                      // outer boundary of convective envelope
    double rInAU  = (rOutAU - rEnvAU);                                                          // inner boundary of convective envelope
    
    if (utils::Compare(rEnvAU/rOutAU, TIDES_MINIMUM_FRACTIONAL_EXTENT) <= 0 || utils::Compare(envMass/Mass(), TIDES_MINIMUM_FRACTIONAL_EXTENT) <= 0 || std::isnan(envMass)) return std::make_tuple(0.0, 0.0, 0.0, 0.0);           // skip calculations if there is no convective envelope (to avoid Imk22 = NaN)

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
    double lConv          = rEnvAU / 2.0;                                                              // set length scale to height of convective envelope
    double tConv          = CalculateEddyTurnoverTimescale();
    double vConv          = lConv / tConv;
    double omegaConv      = 1.0 / tConv;                                                         // absent factor of 2*PI, following Barker (2020)
    double vl             = vConv * lConv;
    double M_2            = Mass() * Mass();

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
 * Dbl_Dbl_Dbl_DblT CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2)
 *
 * @param   [IN]    p_Omega                     Orbital angular frequency (1/yr)
 * @param   [IN]    p_SemiMajorAxis             Semi-major axis of binary (AU)
 * @param   [IN]    p_M2                        Mass of companion star (Msol)
 * @return                                      [(1,0), (1,2), (2,2), (3,2)] Imaginary components of the 
 *                                              potential tidal Love number (unitless)
 */
Dbl_Dbl_Dbl_DblT BaseStar::CalculateImKnmTidal(const double p_Omega, const double p_SemiMajorAxis, const double p_M2) const {
    
    double Imk10Dynamical, Imk12Dynamical, Imk22Dynamical, Imk32Dynamical;
    std::tie(Imk10Dynamical, Imk12Dynamical, Imk22Dynamical, Imk32Dynamical) = CalculateImKnmDynamical(p_Omega, p_SemiMajorAxis, p_M2);

    double Imk10Equilibrium, Imk12Equilibrium, Imk22Equilibrium, Imk32Equilibrium;
    std::tie(Imk10Equilibrium, Imk12Equilibrium, Imk22Equilibrium, Imk32Equilibrium) = CalculateImKnmEquilibrium(p_Omega, p_SemiMajorAxis, p_M2);
    
    // return combined ImKnm terms;
    return std::make_tuple(Imk10Dynamical + Imk10Equilibrium, Imk12Dynamical + Imk12Equilibrium, Imk22Dynamical + Imk22Equilibrium, Imk32Dynamical + Imk32Equilibrium);
}




