#pragma once

#include "constants.h"
#include "typedefs.h"
#include "Remnants.h"


class BaseStar;

class BH: virtual public BaseStar, public Remnants {
    
public:
    
    BH() { m_InterimState.SetStellarType(STELLAR_TYPE::BLACK_HOLE); };
    
    BH(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), Remnants(p_BaseStar) {
        m_InterimState.SetStellarType(STELLAR_TYPE::BLACK_HOLE);    // Set stellar type
        if (p_Initialise) Initialise(); // Initialise if required
    }

    
private:
    
    void Initialise() {

        // set internal properties to zero to avoid meaningless values
        m_InterimState.SetAge(0.0);
        m_InterimState.SetCOCoreMass(0.0);
        m_InterimState.SetHeCoreMass(0.0);
        m_InterimState.SetCoreMass(0.0);
        m_InterimState.SetMassEffectiveInitial(0.0);
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

    GNU_PURE TimescalesT CalculateTimescales(const double p_Mass, const TimescalesT& p_tScales) const override { // JR FIX THIS: DONE
        // Not meaningful for BH, so we just return current timescales
        return Timescales();
    }

    
    //////////////////////////////////////////////////
    //   KICKS                                      //
    //////////////////////////////////////////////////

    COMPAS_PURE static double CalculateSNKickWeighting(const double p_Mass, const double p_FallbackFraction); // JR FIX THIS: DONE

    
    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////
   
    GNU_CONST double CalculateLuminosity_Hurley2000() const override { return 1.0E-10; }                                            // Hurley et al. 2000, eq 96 // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   MASS LOSS / ACCRETION                      //
    //////////////////////////////////////////////////

    GNU_CONST double CalculateEddingtonCriticalRate(const double p_Mass) const override { return 2.6E-8 * p_Mass * MYR_TO_YEAR; }   // e.g., Marchant+, 2017, Eq. 3, assuming accretion efficiency of 10% // JR FIX THIS: DONE
    GNU_PURE  double CalculateEddingtonCriticalRate() const override { return CalculateEddingtonCriticalRate(Mass()); } // JR FIX THIS: DONE


    GNU_CONST MassLossT CalculateMassLossRate() const override { return std::make_tuple(0.0, ML_TYPE::NONE); }                      // Ensure BHs don't lose mass via winds // JR FIX THIS: DONE

    GNU_PURE  static double CalculateNeutrinoMassLoss(const double p_BaryonicMass); // JR FIX THIS: DONE

    
    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////

    // Radius of Schwarzschild black hole, per Hurley et al. 2000, eq 94.
    // Note that the Schwarzschild radius is not correct for a rotating BH.
    GNU_CONST double CalculateRadius_Hurley2000(const double p_Mass) const override { return 4.24E-6 * p_Mass; } // JR FIX THIS: DONE
    GNU_PURE  double CalculateRadius_Hurley2000() const override { return CalculateRadius_Hurley2000(Mass()); } // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   SPIN PARAMETERS                            //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateMomentOfInertia() const override { return (2.0 / 5.0) * Mass() * Radius() * Radius(); }                // MoI for solid sphere *Ilya* JR: that's not really right, is it? // JR FIX THIS: DONE
    
};


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       KICKS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateSNKickWeighting
 *
 * @brief
 * Calculate the weighting to be applied to the kick magnitude for a Black Hole, BH,
 * based on the user-supplied black hole kicks mode (via the option `--black-hole-kicks-mode`).
 *
 * Option values are:
 *
 *    FALLBACK: Black holes receive a kick down-weighted by the amount of mass falling back onto them
 *    FULL    : Black holes receive the same kicks as neutron stars (weighting = 1.0)
 *    REDUCED : Black holes receive the same momentum kick as a neutron star, but down-weighted by the black hole mass
 *    ZERO    : Black holes receive zero natal kick (weighting = 0.0)
 *
 *
 * double CalculateSNKickWeighting(const double p_Mass, const double p_FallbackFraction)
 *
 * @param       p_Mass                          Mass of the remnant (Msol)
 * @param       p_FallbackFraction              Fraction of mass that falls back onto the proto-compact object [0, 1]
 * @return                                      Kick magnitude
 */
inline double BH::CalculateSNKickWeighting(const double p_Mass, const double p_FallbackFraction) {

    double weighting = 0.0;                                                                     // Default return value

    switch (OPTIONS->BlackHoleKicksMode()) {                                                    // Which BH kicks mode?

        case BLACK_HOLE_KICKS_MODE::FALLBACK: weighting = 1.0 - p_FallbackFraction;   break;    // Using the so-called 'fallback' mode for BH kicks
        case BLACK_HOLE_KICKS_MODE::FULL    : weighting = 1.0;                        break;    // Full kick - no adjustment necessary
        case BLACK_HOLE_KICKS_MODE::REDUCED : weighting = NEUTRON_STAR_MASS / p_Mass; break;    // Kick is reduced by the ratio of the neutron star mass to the black hole mass
        case BLACK_HOLE_KICKS_MODE::ZERO    : weighting = 0.0;                        break;    // No kick
    
        default:                                                                                // Unexpected prescription
            // The only way this can happen is if the BLACK_HOLE_KICKS_MODE served to this function is
            // not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR_STATIC(ERROR::UNEXPECTED_BH_KICKS_MODE);                                // Throw error
    }

    return weighting;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               MASS LOSS / ACCRETION                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateNeutrinoMassLoss
 *
 * @brief
 * Calculate the gravitational mass due to neutrino mass loss.
 * 
 *
 * double CalculateNeutrinoMassLoss(const double p_BaryonicMass)
 *
 * @param       p_BaryonicMass                  Baryonic remnant mass (Msol)
 * @return                                      Gravitational mass of remnant (Msol)
 */
inline double BH::CalculateNeutrinoMassLoss(const double p_BaryonicMass) {

    double gravitationalMass = 0.0;                                         // Default return value

    switch (OPTIONS->NeutrinoMassLossAssumptionBH()) {                      // Which prescription?

        case NEUTRINO_ML_PRESCRIPTION::FIXED_FRACTION:                      // FIXED FRACTION
            gravitationalMass = p_BaryonicMass * (1.0 - OPTIONS->NeutrinoMassLossValueBH());
            break;

        case NEUTRINO_ML_PRESCRIPTION::FIXED_MASS:                          // FIXED MASS
            gravitationalMass = p_BaryonicMass - OPTIONS->NeutrinoMassLossValueBH();
            break;
    
        default:                                                            // Unexpected prescription
            // The only way this can happen is if the NEUTRINO_ML_PRESCRIPTION served to this function
            // is not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because OPTIONS served a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the OPTIONS code to prevent it serving a prescription that is not handled by this
            // function.
            THROW_ERROR_STATIC(ERROR::UNEXPECTED_NEUTRINO_ML_PRESCRIPTION); // Throw error
    }

    return gravitationalMass;
}
