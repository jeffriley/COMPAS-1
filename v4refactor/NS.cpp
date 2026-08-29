#include "NS.h"






///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


/*
 * CalculateMagneticFieldAndSpin
 *
 * @brief
 * Calculate the magnetic field and spin of a neutron star, in an interacting binary system,
 * that is undergoing mass-transfer-induced spin change due to Roche lobe overflow or a common
 * envelope event.
 * 
 * Note: this function operates in CGS units: both input parameters and the returned value are in CGS units.
 * *Ilya* why do we do this - can't we be consistent with the rest of the code?
 * 
 * Note: when the star is not spinning (frequency or magnetic field == 0), the spin period is set to infinity,
 *       and all other parameters to 0.
 * 
 *
 * PulsarDetailsT CalculateMagneticFieldAndSpin(
 *     const double            p_Mass,
 *     const double            p_Radius,
 *     const bool              p_CommonEnvelope,
 *     const bool              p_RecycledNS,
 *     const double            p_Stepsize,
 *     const double            p_MassGain,
 *     const double            p_Epsilon,
 *     const PulsarDetailsT&   p_pulsarDetails,
 *     const PulsarConstantsT& p_pulsarConstants 
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_CommonEnvelope                Boolean flag indicating whether there there is a common envelope
 * @param       p_RecycledNS                    Boolean flag indicating whether this star is/was a recycled NS
 * @param       p_Stepsize                      Timestep size for integration (seconds)
 * @param       p_MassGain                      Required accretor mass gain for timestep (g)
 * @param       p_Epsilon                       Efficiency factor allowing for uncertainties of coupling magnetic field and matter
 * @param       p_PulsarDetails                 Struct containing current pulsar parameters
 * @param       p_PulsarConstants               Struct containing pulsar constants
 * @return                                      Struct containing calculated pulsar parameters
 */
PulsarDetailsT NS_Constituent::CalculateMagneticFieldAndSpin(
    const double            p_Mass,
    const double            p_Radius,
    const bool              p_CommonEnvelope,
    const bool              p_RecycledNS,
    const double            p_Stepsize,
    const double            p_MassGain,
    const double            p_Epsilon,
    const PulsarDetailsT&   p_PulsarDetails,
    const PulsarConstantsT& p_PulsarConstants
) const {

    PulsarDetailsT pDet = p_PulsarDetails;                                                                                          // Local copy of current pulsar details

    if (p_RecycledNS || (p_CommonEnvelope && CDOUBLE(p_MassGain) != 0.0)) {                                                         // Isolated pulsar?
                                                                                                                                    // No
        const double massG     = Mass() * MSOL_TO_G;                                                                                // Mass in g
        const double radiusCM  = Radius() * RSOL_TO_CM;                                                                             // Radius in cm

        if (p_CommonEnvelope && (OPTIONS->NeutronStarAccretionInCE() == NS_ACCRETION_IN_CE::SURFACE)) {                             // Mass transfer through CE when accretion happens at the surface of the NS

            pDet.momentOfInertia   = NS::CalculateMomentOfInertia(massG, radiusCM);
            const double jAcc      = std::sqrt(G_CGS * massG * radiusCM) * p_MassGain;    
            pDet.angularMomentum += jAcc;                                                                                           // Angular momentum of the accreted material as it falls onto the surface of the NS
        } 
        else if ((p_MassGain > 0.0)        &&
                 (pDet.spinPeriod > 0.001) && 
                 (!p_CommonEnvelope || (p_CommonEnvelope && OPTIONS->NeutronStarAccretionInCE() == NS_ACCRETION_IN_CE::DISK))) {    // Pulsar recycling through accretion

            // Recycling happens in a binary system with mass transfer for pulsars
            // with spin period larger than 1 ms.
            // The pulsar being recycled is either in a common envelope, or should
            // have started the recycling process in previous time steps.
       
            // Solve for the angular momentum of the NS after accretion.
            // The accretor will gain p_MassGain g over p_Stepsize seconds
            const double massFinal = massG + p_MassGain;                                                                            // Required final mass of NS (after accretion) in g
            const double radius_12 = utils::IntPow(radiusCM, 12);                                                                   // For performance - do it once
            const double mDot_2    = utils::IntPow(p_MassGain / p_Stepsize, 2);                                                     // (Required mass transfer rate (g s^-1))^2 - for performance - do it once

            // Calculate initial mass slice size for integration
            const double jAcc      = CalculateDeltaJByAccretion(massG, radius_12, pDet.magField, pDet.spinFrequency, mDot_2, p_Epsilon);
            const double massSlice = std::fabs(p_PulsarDetails.angularMomentum / 1000.0 / jAcc);                                    // abs(Jx10^-3 / dJdM)

            // Use the boost ODE solver for speed and accuracy

            // Define the ODE

            // Initial state
            state_type x(1);
            x[0] = p_PulsarDetails.angularMomentum;                                                                                 // Angular momentum

            // ODE
            struct ode {
                double p_Mass, p_Radius, p_Radius_12, p_MagField, p_Mdot_2, p_Epsilon;
                ode(double mass, double radius, double radius_12, double magField, double mDot_2, double epsilon) :
                    p_Mass(mass), p_Radius(radiusCM), p_Radius_12(radius_12), p_MagField(magField), p_Mdot_2(mDot_2), p_Epsilon(epsilon) { }

                // x is the current state of the ODE (x[0] = angular momentum J)
                // dxdm is the change of state wrt mass (dxdm[0] = dJdm)
                // p_MassDelta is the cumulative change in mass of the NS
                void operator () (const state_type& x, state_type& dxdm, double p_MassDelta ) const {
                    const double m = p_Mass + p_MassDelta;
                    const double B = (p_MagField - p_PulsarConstants.minMagField) * std::exp(-p_MassDelta / p_PulsarConstants.decayMassScale) + p_PulsarConstants.minMagField;
                    const double f = x[0] / NS::CalculateMomentOfInertia(m, p_Radius);

                    dxdm[0] = CalculateDeltaJByAccretion(m, p_Radius_12, B, f, p_Mdot_2, p_Epsilon);                                    
                }
            };

            // Integrate
            controlled_stepper_type stepper;
            (void)integrate_adaptive(stepper, ode{ massG, radiusCM, radius_12, pDet.magField, mDot_2, p_Epsilon }, x, 0.0, p_MassGain, massSlice);
                
            // Final values
            pDet.angularMomentum = x[0];
            pDet.momentOfInertia = NS::CalculateMomentOfInertia(massFinal, radiusCM);
        }      

        if (pDet.magField < p_PulsarConstants.minMagField) {                                                                        // Magnetic field < lower limit?
            pDet.magField = p_PulsarConstants.minMagField;                                                                          // Yes - set it to the lower limit
        }
        else {                                                                                                                      // No - recalculate it
            pDet.magField = (pDet.magField - p_PulsarConstants.minMagField) * std::exp(-p_MassGain / p_PulsarConstants.decayMassScale) + p_PulsarConstants.minMagField; // eq. 12 in arxiv:1912.02415 
        }
        
        pDet.spinFrequency = pDet.angularMomentum / pDet.momentOfInertia;
        pDet.spinPeriod    = _2_PI / pDet.spinFrequency;
        const double fDot  = (pDet.spinFrequency - p_Pulsardetails.spinFrequency) / p_Stepsize;       
        pDet.spinDownRate  = -fDot * pDet.spinPeriod * pDet.spinPeriod / _2_PI;

    }
    else {                                                                                                                          // Otherwise...    
        pDet = NS::CalculateSpinDownValues(p_Radius, p_Stepsize, p_PulsarDetails, p_PulsarConstants);                               // ...treat the pulsar as isolated - spin down
    }

    return pDet;
}


/*
 * CalculateMassAccretedForCEE
 *
 * Calculate the mass accreted by a neutron star or black hole during common envelope evolution,
 * per the user-specified prescription (see --common-envelope-mass-accretion-prescription).
 * 
 * 
 * double CalculateMassAccretedForCEE(const double p_Mass, const double p_CompanionMass, const double p_CompanionRadius, const double p_CompanionEnvelope) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CompanionMass                 Mass of the companion star (Msol)
 * @param       p_CompanionRadius               Radius of the companion star (Rsol)
 * @param       p_CompanionEnvelope             Envelope of the companion pre-CE (Msol)
 * @return                                      Mass accreted by the star (Msol)
 */
double NS_Constituent::CalculateMassAccretedForCEE(const double p_Mass, const double p_CompanionMass, const double p_CompanionRadius, const double p_CompanionEnvelope) const {

     double deltaMass = 0.0;                                                                // Default return value

     switch (OPTIONS->CommonEnvelopeMassAccretionPrescription()) {                          // Which prescription?

        case CE_ACCRETION_PRESCRIPTION::ZERO:                                               // ZERO  //<<<<<<<<<<<<< JR FIX THIS - NONE? <<<<<<<<<<<<<
            deltaMass = 0.0;                                                                // No acceretion
            break;

        case CE_ACCRETION_PRESCRIPTION::CONSTANT:                                           // CONSTANT
            // User-supplied value (see --common-envelope-mass-accretion-constant)
            deltaMass = OPTIONS->CommonEnvelopeMassAccretionConstant();
            break;

        case CE_ACCRETION_PRESCRIPTION::UNIFORM:                                            // UNIFORM
            // Uniform random distribution - Oslowski+ (2011)
            deltaMass = RAND->Random(OPTIONS->CommonEnvelopeMassAccretionMin(), OPTIONS->CommonEnvelopeMassAccretionMax());
            break;

        case CE_ACCRETION_PRESCRIPTION::MACLEOD: {                                          // MACLEOD
            // Linear regression estimated from Macleod+ (2015)
            constexpr double mm = -1.0714285714285712E-05;                                  // Gradient of the linear fit for gradient
            constexpr double cm =  0.00012057142857142856;                                  // Intercept of the linear fit for gradient
            constexpr double mc =  0.01588571428571428;                                     // Gradient of the linear fit for intercept
            constexpr double cc = -0.15462857142857137;                                     // Intercept of the linear fit for intercept
            const     double m  = mm * p_CompanionMass + cm;                                // Gradient of linear fit for mass
            const     double c  = mc * p_CompanionMass + cc;                                // Intercept of linear fit for mass

            // Calculate mass to accrete and clamp to minimum and maximum from program options
            deltaMass = std::clamp(m * p_CompanionRadius + c, OPTIONS->CommonEnvelopeMassAccretionMin(), OPTIONS->CommonEnvelopeMassAccretionMax());

            } break;

        case CE_ACCRETION_PRESCRIPTION::CHEVALIER:                                          // CHEVALIER
            // Model 2 from van Son et al. 2020
            deltaMass = (p_Mass * p_CompanionMass) / (2.0 * (p_Mass + p_CompanionMass));    // Hoyle-Lyttleton accretion rate times inspiral time
            break;

        default:                                                                            // Unexpected prescription
            // The only way this can happen is if the CE_ACCRETION_PRESCRIPTION served to this function
            // is not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR(ERROR::UNEXPECTED_CE_ACCRETION_PRESCRIPTION);                       // Throw error
    }

    // Clamp the mass accretion to be no more than the envelope of the companion pre-CE
    deltaMass = std::min(p_CompanionEnvelope, deltaMass);

    return deltaMass;
}
