#pragma once

#include "constants.h"
#include "typedefs.h"
#include "sampling.h"

#include "Rand.h"
#include "Remnants.h"
#include "ConstituentStar.h"
#include "BH.h"


class BaseStar;
class Remnants;

class NS: virtual public BaseStar, public Remnants {

public:

    NS() { m_InterimState.SetStellarType(STELLAR_TYPE::NEUTRON_STAR); } // Set stellar type  <<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS
    
    NS(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {
        m_InterimState.SetStellarType(STELLAR_TYPE::NEUTRON_STAR); // Set stellar type


        if (p_Initialise) Initialise(); // Initialise if required
    }


private:
    
    void Initialise() {
        
        // Set attributes to zero to avoid meaningless values
        m_InterimState.SetAge(0.0);
        m_InterimState.SetCOCoreMass(0.0);
        m_InterimState.SetHeCoreMass(0.0);
        m_InterimState.SetCoreMass(0.0);
        m_InterimState.SetMassEffectiveInitial(0.0);
        
        // Calculate pulsar details at birth of pulsar
        m_InterimState.SetPulsarDetails(CalculateBirthParameters(Mass(), Radius()));
    }


protected:


    // Member functions (not getters or setters)
    //
    // VIRTUAL FUNCTIONS may be (are expected to be) overridden by derived classes.
    // When overriding virtual functions in a derived class, use the "override" attribute.
    //
    // NON-VIRTUAL FUNCTIONS should not be overridden (declared separately) by derived classes.
    // While it is legal in C++ to declare the same (non-virtual) function in multiple classes,
    // (aka "shadowing", or "hiding"), we discourage it.  Non-virtual functions are statically
    // bound, and as such, especially with indirection, may not produce expected results.

    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////

    GNU_PURE  double ChooseTimestep() const override { return ChooseTimestep(Age()); } // JR FIX THIS: DONE
    GNU_CONST double ChooseTimestep(const double p_Age) const; // JR FIX THIS: DONE

    GNU_PURE  TimescalesT CalculateTimescales(const double p_Mass, const TimescalesT& p_tScales) const override { // JR FIX THIS: DONE
        // Not meaningful for NS, so we just return current timescales
        return Timescales();
    }

    
    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    GNU_PURE  double CalculateLuminosity_Hurley2000() const override { return CalculateLuminosity_Hurley2000(Mass(), Age()); } // JR FIX THIS: DONE
    GNU_CONST double CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Age) const; // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   MASS LOSS / ACCRETION                      //
    //////////////////////////////////////////////////

    GNU_CONST MassLossT CalculateMassLossRate() const override { return std::make_tuple(0.0, ML_TYPE::NONE); } // Ensure BHs don't lose mass via winds // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   MISC. CONSTANTS / PARAMETERS               //
    //////////////////////////////////////////////////

    PulsarDetailsT CalculateBirthParameters(const double p_Mass, const double p_Radius) const; // JR FIX THIS: DONE

    
    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////

    double CalculateRadius() const override { return DetermineRadiusInKM(Mass()) * KM_TO_RSOL;} // JR FIX THIS: DONE


    static double DetermineRadiusInKM(const double p_Mass); // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   SPIN PARAMETERS                            //
    //////////////////////////////////////////////////

    GNU_CONST static double CalculateMomentOfInertia(const double p_Mass, const double p_Radius); // JR FIX THIS: DONE
    GNU_PURE         double CalculateMomentOfInertia() const override { // JR FIX THIS: DONE
        return CalculateMomentOfInertia(const double Mass(), const double Radius());
    }


    GNU_CONST static double CalculateSpinDownRate(const double p_Radius, const double p_SpinPeriod, const double p_MoI, const double p_MagField) const; // JR FIX THIS: DONE


    GNU_CONST static PulsarDetailsT CalculateSpinDownValues(const double p_Radius, const double p_Stepsize, const PulsarDetailsT& p_PulsarDetails, const PulsarConstantsT& p_PulsarConstants) const; // JR FIX THIS: DONE



 




            STELLAR_TYPE    EvolveToNextPhase()                                 { return STELLAR_TYPE::BLACK_HOLE; }

    
            bool            ShouldEvolveOnPhase() const                         { return (Mass() <= OPTIONS->MaximumNeutronStarMass()); }

};


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * ChooseTimestep
 *
 * @brief
 * Choose timestep for evolution.
 * 
 * Pulsars evolve very fast when they are first born, and slower as they age, so the
 * timestep is chosen to be small when the pulsar is young, and slowly increases as
 * the pulsar ages.
 * 
 * The returned value will be clamped to minimum NUCLEAR_MINIMUM_TIMESTEP
 *
 *
 * double ChooseTimestep(const double p_Age)
 *
 * @param       p_Age                           Age of the star (Myr)
 * @return                                      Suggested timestep (Myr)
 */
inline double NS::ChooseTimestep(const double p_Age) const {

    double dt = 500.0;                                      // Default return value

         if (p_Age < 0.01 ) dt = 0.001;
    else if (p_Age < 0.1  ) dt = 0.01;
    else if (p_Age < 1.0  ) dt = 0.1;
    else if (p_Age < 10.0 ) dt = 1.0;
    else if (p_Age < 500.0) {                               // Interporlate
        // Slope = log10(500.0) / (log10(500.0) - 1.0) = 1.58859191006
        // Step  = log10(age) - 1.0
        dt = PPOW(10.0, 1.58859191006 * (std::log10(p_Age) - 1.0));
    }

    return std::clamp(dt, 0.0, NUCLEAR_MINIMUM_TIMESTEP);   // Clamp to minimum NUCLEAR_MINIMUM_TIMESTEP (and 0.0)
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                           MISC. CONSTANTS / PARAMETERS                            //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateBirthParameters
 *
 * @brief
 * Calculate pulsar parameters at birth (of the pulsar).
 * 
 * 
 * PulsarDetailsT CalculateBirthParameters(const double p_Mass, const double p_Radius) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @return                                      Struct containing calculated pulsar parameters
 */
inline PulsarDetailsT NS::CalculateBirthParameters(const double p_Mass, const double p_Radius) const {

    PulsarDetailsT pDet;

    pDet.magneticField     = PPOW(10.0, sampling::SampleBirthMagneticField());  // Gauss
    pDet.spinPeriod        = sampling::SampleBirthSpinPeriod();                 // Seconds (s)
    pDet.birthSpinPeriod   = pDet.spinPeriod;                                   // At birth of pulsar, not birth of star. Seconds (s)
    pDet.spinFrequency     = _2_PI / pDet.spinPeriod;                           // s^-1
    pDet.momentOfInertia   = CalculateMomentOfInertia(p_Mass * MSOL_TO_G, p_Radius * RSOL_TO_CM); // g cm^2
    pDet.spinDownRate      = CalculateSpinDownRate(p_Radius * RSOL_TO_CM, pDet.spinPeriod, pDet.momentOfInertia, pDet.magneticField); // s^-2  
    pDet.birthSpinDownRate = pDet.spinDownRate;                                 // At birth of pulsar, not birth of star. s^-2 
    pDet.angularMomentum   = pDet.momentOfInertia * pDet.spinFrequency;         // g cm^2 s^-1

    return pDet;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateLuminosity_Hurley2000
 *
 * @brief
 * Calculate the luminosity of a neutron star, per Hurley et al. 2000, eq 93
 * 
 * 
 * double CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Age) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Age of the star (Myr)
 * @return                                      Luminosity of the star (Lsol)
*/
inline double CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Age) const {
    const double t = std::max(p_Age, 0.1);
    return 0.02 * PPOW(p_Mass, 2.0 / 3.0) / (t * t); // Hurley et al. 2000, eq 93
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * DetermineRadiusInKM
 *
 * @brief
 * Determine radius of a neutron star per user-supplied NS equation-of-state
 * (see option --neutron-star-equation-of-state).  The radius could be sampled
 * from table ARP3_MASS_RADIUS_RELATION defined in constants.h.
 * 
 * 
 * static double DetermineRadiusInKM(const double p_Mass)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Radius of the star (km)
 */
inline double NS::DetermineRadiusInKM(const double p_Mass) {

    double radius = 0.0;                                            // Default return value

    switch (OPTIONS->NeutronStarEquationOfState()) {                // Which equation-of-state?

        case NS_EOS::SSE:                                           // SSE
            radius = 10.0;
            break;

        case NS_EOS::ARP3: {                                        // ARP3

            // Use table ARP3_MASS_RADIUS_RELATION defined in constants.h.
            // Sample masses inside table limits, but don't extrapolate outside table limits
            // (masses outside table are just set to table extremities)

            std::map<double, double>::const_iterator iter;          // Map iterator
            
            iter = ARP3_MASS_RADIUS_RELATION.begin();               // First element
            const double minMass     = iter->first;                 // Mass
            const double radiusAtMin = iter->second;                // Radius

            iter = ARP3_MASS_RADIUS_RELATION.end();                 // Last + 1 element
            const double maxMass     = std::prev(iter)->first;      // Last element mass
            const double radiusAtMax = std::prev(iter)->second;     // Last element radius

                 if (p_Mass < minMass) radius = radiusAtMin;        // Clamp to lower bound
            else if (p_Mass > maxMass) radius = radiusAtMax;        // Clamp to upper bound
            else                       radius = utils::SampleFromTabulatedCDF(p_Mass, ARP3_MASS_RADIUS_RELATION); // Sample
        } break;

        default:                                                    // Unexpected prescription
            // The only way this can happen is if the NS_EOS served to this function is not accounted
            // for in this code.  We should not default here, with or without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR_STATIC(ERROR::UNEXPECTED_NS_EOS);           // Throw error
	}

	return radius;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  SPIN PARAMETERS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateMomentOfInertia
 *
 * @brief
 * Calculate the moment of inertia for a neutron star using a model independent relation between
 * the moment of inertia, mass and radius of a neutron star, per Raithel et al. 2016, eq 8.
 * (see https://arxiv.org/abs/1603.06594)
 *
 * Note: this function operates in CGS units: both input parameters and the returned value are in CGS units.
 * *Ilya* why do we do this - can't we be consistent with the rest of the code?
 * 
 * This function is required to be static for the Boost ODE (used in NS_Constituent::CalculateMagneticFieldAndSpin()).
 * 
 * 
 * static double CalculateMomentOfInertia(const double p_Mass, const double p_Radius)
 *
 * @param       p_Mass                          Mass of the star (g)
 * @param       p_Radius                        Radius of the star (cm)
 * @return                                      Moment of inertia (g cm^2)
 */
inline double NS::CalculateMomentOfInertia(const double p_Mass, const double p_Radius) {
    const double m_r = (p_Mass / MSOL_TO_G) / (p_Radius / KM_TO_CM);
    return 0.237 * p_Mass * p_Radius * p_Radius * (1.0 + (4.2 * m_r) + 90.0 * m_r * m_r * m_r * m_r);
}


/*
 * CalculateSpinDownRate
 *
 * @brief
 * Calculate the spin down rate, in s^-2, for an isolated neutron star.
 *
 * See Equation 5 in https://arxiv.org/abs/2406.11428
 * Note that magnetic and rotational axes are orthogonal, leading to sin^2(alpha) = 1 in this equation. 
 * A model with evolving alpha will be implemented in a future version. 
 * 
 * Evolution of the inclination between pulsar magnetic and rotational axes will be considered in a future version. 
 * 
 * Note: this function operates in CGS units: both input parameters and the returned value are in CGS units.
 * *Ilya* why do we do this - can't we be consistent with the rest of the code?
 *
 *
 * static double CalculateSpinDownRate(const double p_Radius, const double p_SpinPeriod, const double p_MoI, const double p_MagField)
 *
 * @param       p_Radius                        Radius of the star (cm)
 * @param       p_SpinPeriod                    Spin period of the star (s). 
 * @param       p_MoI                           Moment of Inertia of the star (g cm^2)
 * @param       p_MagField                      Magnetic field of the star (Gauss)
 * @return                                      Spin down rate (spin period derivative) (s^-2)
 */
inline double NS::CalculateSpinDownRate(const double p_Radius, const double p_SpinPeriod, const double p_MoI, const double p_MagField) {

   constexpr double _8_PI_2    = 8.0 * PI_2;
   constexpr double _3_C_3     = 3.0E6 * C * C * C; // 3.0 * (C * 100.0) * (C * 100.0) * (C * 100.0)
   const     double pDotTop    = _8_PI_2 * utils::intPow(p_Radius, 6) * p_MagField * p_MagField;
   const     double pDotBottom = _3_C_3 * p_MoI * p_SpinPeriod;
   
   return pDotTop / pDotBottom;                                                                         
}


/*
 * CalculateSpinDownValues
 *
 * @brief
 * Calculate attribute values for an isolated pulsar spinning down.
 * If the pulsar is not spinning, nothing changes.
 *  
 * We assume the rotational and magnetic axis are orthogonal and are not evolved
 * in the current model.
 * 
 * A model with evolving alpha will be implemented in a future version.
 * *Ilya* is the above statement likely to be true, and do we need it here?
 *
 *
 * static PulsarDetailsT CalculateSpinDownValues(const double p_Radius, const double p_Stepsize, const PulsarDetailsT& p_PulsarDetails, const PulsarConstantsT& p_PulsarConstants)
 *
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Stepsize                      Timestep size for integration (seconds)
 * @param       p_PulsarDetails                 Struct containing current pulsar parameters
 * @param       p_PulsarConstants               Struct containing pulsar constants
 * @return                                      Struct containing calculated pulsar parameters
 */
inline PulsarDetailsT NS::CalculateSpinDownValues(const double p_Radius, const double p_Stepsize, const PulsarDetailsT& p_PulsarDetails, const PulsarConstantsT& p_PulsarConstants) {

    PulsarDetailsT pDet = p_Pulsardetails;          // Local copy of current details

    // Calculate the decay of magnetic field for an isolated neutron star
    // (see Eq 6 in arXiv:0903.3538v2)
    const double magField = pDet.magField;          // Start with value at beginning of the timestep
    if (p_PulsarConstants.minMagField < magField) { // magfield >= lower limit?
                                                    // Yes - update magfield
        magField = p_PulsarConstants.minMagField + (pDet.magField - p_PulsarConstants.minMagField) * std::exp(-p_Stepsize / p_PulsarConstants.decayTimeScale);
    }

    // Calculate the spin down values for isolated neutron stars (see Eq 3 in arxiv:2406.11428).
    // Note that magnetic and rotational axes are orthogonal, leading to sin^2(alpha) = 1 in this equation. 

    constexpr double pi_2_8 = 8.0 * PI_2;
    constexpr double C_3_3  = 3.0E6 * C * C * C;    // 3.0 * (C * 100.0) * (C * 100.0) * (C * 100.0)

    const double radius   = p_Radius * RSOL_TO_CM;
    const double c2       = (pi_2_8 * utils::intPow(radius, 6)) / (C_3_3 * pDet.momentOfInertia);
    const double term1    = p_PulsarConstants.minMagField * p_PulsarConstants.minMagField * p_Stepsize;
    const double term2    = p_PulsarConstants.decayTimeScale * p_PulsarConstants.minMagField * ( magField  - pDet.magField);
    const double term3    = (p_PulsarConstants.decayTimeScale / 2.0) * ((magField * magField) - (pDet.magField * pDet.magField));
    const double Psquared = 2.0 * c2 * (term1 - term2 - term3) + (pDet.spinPeriod * pDet.spinPeriod);
    
    pDet.spinPeriod       = std::sqrt(Psquared);
    pDet.spinFrequency    = _2_PI / pDet.spinPeriod;
    pDet.spinDownRate     = CalculateSpinDownRate(radius, pDet.spinPeriod, pDet.momentOfInertia, magField); 
    pDet.angularMomentum  = pDet.spinFrequency * pDet.momentOfInertia;
    pDet.magField         = magField;

    return pDet;
}





///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



class NS_Constituent: virtual public BaseBinaryStar {

public:

    NS_Constituent() {

        m_InterimState.SetStellarType(STELLAR_TYPE::NEUTRON_STAR); // Set stellar type  <<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS

    };
    
    NS_Constituent(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseBinaryStar(p_BaseBinaryStar) {
        m_InterimState.SetStellarType(STELLAR_TYPE::NEUTRON_STAR); // Set stellar type


        if (p_Initialise) Initialise(); // Initialise if required
    }


protected:
    
    void Initialise() {
        
    }



    GNU_CONST inline double CalculateCriticalMassRatio_Hurley2002() const override { return 0.0; }

    GNU_CONST static double CalculateDeltaJByAccretion(
        const double p_Mass,
        const double p_Radius_12,
        const double p_MagField,
        const double p_SpinFrequency,
        const double p_mDot_2,
        const double p_Epsilon
    );

    GNU_CONST PulsarDetailsT CalculateMagneticFieldAndSpin(
        const double            p_Mass,
        const double            p_Radius,
        const bool              p_CommonEnvelope,
        const bool              p_RecycledNS,
        const double            p_Stepsize,
        const double            p_MassGain,
        const double            p_Epsilon,
        const PulsarDetailsT&   p_pulsarDetails,
        const PulsarConstantsT& p_pulsarConstants 
    ) const;

    GNU_PURE double CalculateMassAccretedForCEE(const double p_Mass, const double p_CompanionMass, const double p_CompanionRadius, const double p_CompanionEnvelope) const;

    GNU_PURE inline MT_CASE DetermineMTTypeAsDonor() const override { return MT_CASE::NONE; } // NSs are never donors

    GNU_PURE  inline double ResolveCommonEnvelopeAccretion() const override { return ResolveCommonEnvelopeAccretion(Mass()); }
    GNU_CONST double ChooseTimestep(const double p_Age) const;

            double ResolveCommonEnvelopeAccretion(const double p_Mass, const double p_CompanionMass, const double p_CompanionRadius, const double p_CompanionEnvelope) const;

/* 
 * Resolve common envelope accretion
 *
 * For stellar types other than Black hole or Neutron Star just set the star's mass to the parameter passed
 * For Black holes or Neutron Stars calculate the mass accreted during a CE
 *
 *
 * double ResolveCommonEnvelopeAccretion(const double p_FinalMass, 
 *                                       const double p_CompanionMass, 
 *                                       const double p_CompanionRadius, 
 *                                       const double p_CompanionEnvelope)
 *
 * @param   [IN]    p_FinalMass                 Mass of the accreting object post mass transfer (Msol) (not used here)
 * @param   [IN]    p_CompanionMass             Mass of the companion
 * @param   [IN]    p_CompanionRadius           Radius of the companion
 * @param   [IN]    p_CompanionEnvelope         Envelope of the companion pre-CE
 * @return                                      Mass delta                                      
 * 
 */
double NS::ResolveCommonEnvelopeAccretion(const double p_FinalMass,
                                          const double p_CompanionMass,
                                          const double p_CompanionRadius,
                                          const double p_CompanionEnvelope) {
    return CalculateMassAccretedForCEE(Mass(), p_CompanionMass, p_CompanionRadius, p_CompanionEnvelope);
}


};


/*
 * CalculateDeltaJByAccretion
 *
 * @brief
 * Calculate the change in angular momentum wrt mass (dJ/dM) of a neutron star
 * when it accretes mass through RLOF.
 * 
 * Note: this function operates in CGS units: both input parameters and the returned value are in CGS units.
 * *Ilya* why do we do this - can't we be consistent with the rest of the code?
 * 
 * See sec. 2.2.1 in arxiv:1912.02415
 * 
 * This function is required to be static for the Boost ODE (used in NS_Constituent::CalculateMagneticFieldAndSpin()).
 * 
 *
 * static double CalculateDeltaJByAccretion(
 *     const double p_Mass,
 *     const double p_Radius_12,
 *     const double p_MagField,
 *     const double p_SpinFrequency,
 *     const double p_mDot_2,
 *     const double p_Epsilon
 * )
 * 
 * @param       p_Mass                          Initial mass of the NS (g)
 * @param       p_Radius_12                     (Radius of the NS (cm))^12 (for performance - so it isn't calculated at every integration step)
 * @param       p_MagField                      NS magnetic field strength at the beginning of accretion (Gauss)
 * @param       p_SpinFrequency                 Angular frequency for the NS at the beginning of accretion (rad/s)
 * @param       p_mDot_2                        (Mass transfer rate (g s^-1))^2  (for performance - so it isn't calculated at every integration step)
 * @param       p_Epsilon                       Efficiency factor allowing for uncertainties of coupling magnetic field and matter
 * @return                                      Change in angular momentum wrt mass (dJ/dM) of NS due to accretion
 */
inline double NS_Constituent::CalculateDeltaJByAccretion(
    const double p_Mass,
    const double p_Radius_12,
    const double p_MagField,
    const double p_SpinFrequency,
    const double p_mDot_2,
    const double p_Epsilon
) {
    // Calculate the Alfven radius for an accreting neutron star.
    // See eq 10 in arxiv:1912.02415.
    const double p              = p_Radius_12 / (p_mDot_2 * p_Mass);
    const double q              = PPOW(p, 1.0 / 7.0);
    const double magneticRadius = ALFVEN_CONST * q * PPOW(p_MagField, 4.0 / 7.0) / 2.0; // Alfven radius / 2.0 (cm)
    
    // Calculate the difference in the keplerian angular velocity at the
    // magnetic radius and surface angular velocity of the NS.
    // See eq 2 in 1994MNRAS.269..455J / eq 9 in arxiv:1912.02415.
    const double omegaK = std::sqrt(G_CGS * p_Mass / magneticRadius) / magneticRadius;  // rad/s
    const double vDiff  = omegaK - p_SpinFrequency;                                     // rad/s

    return p_Epsilon * vDiff * magneticRadius * magneticRadius;                         // Eq 12 in arXiv:0805.0059 / eq 8 in arxiv:1912.02415
}
