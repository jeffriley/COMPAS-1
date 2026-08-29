#pragma once

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "EAGB.h"


class BaseStar;
class EAGB;

class TPAGB: virtual public BaseStar, public EAGB {

public:

    TPAGB() { m_InterimState.SetStellarType(STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH); };
    
    TPAGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), EAGB(p_BaseStar, false) {
        m_InterimState.SetStellarType(STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH);                                                                                                                  // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                                                     // Initialise if required
    }

    TPAGB* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        TPAGB* clone = new TPAGB(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static TPAGB* Clone(TPAGB& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        TPAGB* clone = new TPAGB(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                                                                              // Initialise timescales
        m_InterimState.SetAge(m_InterimState.HurleyTimescales(HURLEY_TIMESCALES::DU2));                                                                                                                     // Set age appropriately (was: m_Age = m_Timescales[HURLEY_TIMESCALES::DU2], where tP -> DU2 under current HURLEY_TIMESCALES enum)
        
        EvolveOnPhase(0.0);
   }


   // member functions - alphabetically




            double          CalculateConvectiveCoreRadius() const                                                   { return CalculateConvectiveCoreRadius(Radius()); }                                     // Last paragraph of section 6 of Hurley+ 2000
            double          CalculateConvectiveCoreRadius(const double p_Radius) const                              { return std::min(5.0 * CalculateRemnantRadius(), p_Radius); }
                                                                // NO-OP


            double          CalculateLuminosityOnPhase(const double p_Time) const;
            double          CalculateLuminosityOnPhase() const                                                      { return CalculateLuminosityOnPhase(Age()); }                                                                                                                                                                                       // Use state-routed values
            double          CalculateLuminosityAtPhaseEnd() const                                                   { return CalculateLuminosityAtPhaseEnd(Luminosity()); }                                                                                                                                                                              // NO-OP
            double          CalculateLuminosityAtPhaseEnd(const double p_Luminosity) const                          { return p_Luminosity; }


            double          CalculateRadiusAtPhaseEnd() const                                                       { return CalculateRadiusAtPhaseEnd(Radius()); }                                                                                                                                                                                      // NO-OP
            double          CalculateRadiusAtPhaseEnd(const double p_Radius) const                                  { return p_Radius; }
            double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const            { return CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, ZDEP->HurleyMassCutoffs(Metallicity(), HURLEY_MCO::HeF), ZDEP->HurleyBCoefficients(Metallicity())); }                            // metallicity-dependent lookups via ZDEPENDENT (was: m_MassCutoffs[MHeF], m_BnCoefficients)
            double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhase(Mass(), Luminosity()); }                                                                                                                                                                            // Use state-routed values
    static  double          CalculateRadiusOnPhase_Static(const double      p_Mass,
                                                          const double      p_Luminosity,
                                                          const double      p_MHeF,
                                                          const DBL_VECTOR &p_BnCoefficients);

            double          CalculateRemnantLuminosity() const;













///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       



    ////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES   //
    ////////////////////////////////////////

    // Tau (relative age) is not used for TPAGB stars in Hurley et al. 2000, so we return 0.0
    GNU_CONST inline double CalculateTau_Hurley2000() const override { return 0.0; } // JR FIX THIS: DONE
    GNU_CONST inline double CalculateTauAtPhaseEnd_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const override { return 0.0; }; // JR FIX THIS: DONE


void CalculateTimescales_Hurley2000() override { m_InterimState.SetHurleyTimescales(CalculateTimescales_Hurley2000(MassEffectiveInitial(), m_InterimState.HurleyGBParamsOrDefault(), m_InterimState.HurleyTimescalesOrDefault())); }
COMPAS_PURE DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const;

void CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales) override;






    



    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////










    ////////////////////////////////////////
    //   MASS                             //
    ////////////////////////////////////////


inline Dbl_DblT CalculateConvectiveEnvelopeMass() const override { return CalculateConvectiveEnvelopeMass(Mass(), CoreMass()); }
GNU_CONST Dbl_DblT CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_CoreMass) const;


        
double CalculateCoreMass() const override;
GNU_CONST double CalculateCoreMass_Hurley2000(const double p_Age, const double p_MassEffectiveInitial, const DBL_VECTOR& p_GBParams, const double p_McDU, const DBL_VECTOR& p_tScales) const;
GNU_CONST double CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const;
using BaseStar::CalculateLuminosity_Hurley2000;
GNU_CONST double CalculateCoreMassPrime_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const;

inline double CalculateCOCoreMass() const override { return CalculateCoreMass(); } // McCO = Mc for TPAGB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed
inline double CalculateHeCoreMass() const override { return CalculateCoreMass(); } // McHe = Mc for TPAGB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed

GNU_CONST double CalculateMcPrime_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const;


inline double CalculateRemnantRadius_Hurley2000() const override { return CalculateRemnantRadius_Hurley2000_Static(CoreMass()); }
GNU_CONST static double CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass);



            double          CalculateCoreMassAtPhaseEnd() const                                                     { return CalculateCoreMassAtPhaseEnd(CoreMass()); }                                                                          // NO-OP
            double          CalculateCoreMassAtPhaseEnd(const double p_CoreMass) const                              { return p_CoreMass; }
            double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Time) const;
            double          CalculateCoreMassOnPhase() const                                                        { return CalculateCoreMassOnPhase(Mass0(), Age()); }                                    // Use state-routed values





///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inline double CalculateCOCoreMassAtPhaseEnd() const override { return std::min(CalculateCoreMass(), Mass()); } // McCO = Mc for TPAGB stars; McCO should be <= M  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed

inline double CalculateHeCoreMassAtPhaseEnd() const override { return CalculateCoreMass(); } // McHe = Mc for TPAGB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed






    



    ////////////////////////////////////////
    //   RADIUS                           //
    ////////////////////////////////////////




            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return Temperature(); }                                                                                                                                                                                                           // NO-OP (deliberately ignores inputs; returns current temperature)
            double          CalculateTemperatureAtPhaseEnd() const                                                  { return CalculateTemperatureAtPhaseEnd(Luminosity(), Radius()); }                                                                                                                                                                   // Use state-routed values



GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const override;

GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::CONVECTIVE; } // Always CONVECTIVE for TPAGB stars

            STELLAR_TYPE    EvolveToNextPhase()                                                                     { return EvolveToNextPhase(StellarType()); }                                                                                                                                                                                        // NO-OP
            STELLAR_TYPE    EvolveToNextPhase(const STELLAR_TYPE p_StellarType)                                     { return p_StellarType; }

            bool            IsEndOfPhase() const                                                                    { return !ShouldEvolveOnPhase(); }                                                      // Phase ends when envelope loss or going supernova
            bool            IsSupernova() const;

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
            void            ResolveHeliumFlash() { }                                                                                                                                                        // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                                   { return ResolveSkippedPhase(StellarType()); }                                                                                                                                                                                      // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase(const STELLAR_TYPE p_StellarType)                                   { return p_StellarType; }

            bool            ShouldEvolveOnPhase() const                                                             { return (utils::Compare(COCoreMass(), Mass()) < 0 && !IsSupernova() && !ShouldEnvelopeBeExpelledByPulsations()); } // Evolve on TPAGB phase if envelope is not lost and not going supernova
            bool            ShouldSkipPhase() const                                                                 { return false; }                                                                       // Never skip TPAGB phase

    COMPAS_PURE double CalculateCELambda_Nanjing_Enhanced(const double p_Mass, const double p_Radius, const double p_CoreMass, const SizeT p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;

};



//// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 

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
 * ChooseTimestep_Hurley2000
 *
 * @brief
 * Choose timestep for evolution
 * See the discussion in Hurley et al. 2000, p21
 * The returned value will be clamped to minimum NUCLEAR_MINIMUM_TIMESTEP
 *
 *
 * double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales)
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Suggested timestep (Myr)
 */
GNU_CONST inline double TPAGB::ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const {
#define tScales(x) p_tScales[static_cast<int>(HURLEY_TIMESCALES::x)]

    const double dtk = 0.02 * ((p_Age <= tScales(Mx_SAGB) ? tScales(Inf1_SAGB) : tScales(Inf2_SAGB)) - p_Age);   // stellar type specific dt
    const double dte = 5.0E-3;                                                                                      // (artificial) time to end of phase (change of stellar type)

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);                                                  // clamp to minimum NUCLEAR_MINIMUM_TIMESTEP

#undef tScales
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
 * Calculate luminosity on the Thermally Pulsing Asymptotic Giant Branch,
 * per Hurley et al. 2000, eq 37
 *
 *
 * double CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_GBParams                      Hurley GB parameters
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      TPAGB luminosity (Lsol)
 */
GNU_CONST inline double TPAGB::CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const {
    return CalculateLuminosity_Hurley2000(CalculateCoreMassPrime_Hurley2000(p_Age, p_GBParams, p_tScales), p_GBParams);
}


/*
 * Calculate luminosity of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass)
 *
 * Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantLuminosity()
 *
 * @return                                      Luminosity of remnant core in Lsol
 */
// JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
double TPAGB::CalculateRemnantLuminosity() const {
    return COWD::CalculateLuminosityOnPhase_Static(m_CoreMass, 0.0, m_Metallicity);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

      
/*
 * CalculateConvectiveEnvelopeMass
 *
 * @brief
 * Calculate the convective envelope mass for Thermally Pulsing Asymptotic Giant Branch, TPAGB, stars.
 * We assume the entire envelope is convective for TPAGB stars.
 *
 * Dbl_DblT CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_CoreMass) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   TPAGB convective envelope mass (Msol)
 *                                                   Maximum TPAGB convective envelope mass (Msol)
 */
inline Dbl_DblT TPAGB::CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_CoreMass) const {
    const double mEnv = p_Mass - p_CoreMass;
    return std::make_tuple(mEnv, mEnv);
}
   

/*
 * CalculateCoreMass
 *
 * @brief
 * Calculate core mass on this phase of the evolution of the star
 * (Thermally Pulsing Asymptotic Giant Branch).
 *
 * Calls relevant core mass function based on the evolutionary mode given in program options.
 * 
 *
 * double CalculateCoreMass() const
 *
 * @return                                      EAGB CO core mass (Msol)
 */
inline double TPAGB::CalculateCoreMass() const { 

    double coreMass;

    switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            coreMass = CalculateCoreMass_Hurley2000(Age(), MassEffectiveInitial(), m_InterimState.HurleyGBParams(), m_InterimState.GBParams(GBP::McDU), m_InterimState.Timescales());
            break;
        
        default:                                                                        // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
    }       

    return coreMass;
}


/*
 * CalculateCoreMass_Hurley2000
 *
 * @brief
 * Calculate core mass on the Thermally Pulsing Asymptotic Giant Branch (TPAGB),
 * per Hurley et al. 2000, just after eq 73
 *
 *
 * double CalculateCoreMass_Hurley2000(
 *     const double      p_Age,
 *     const double      p_MassEffectiveInitial,
 *     const DBL_VECTOR& p_GBParams,
 *     const double      p_McDU,
 *     const DBL_VECTOR& p_tScales
 * ) const
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_MassEffectiveInitial          Effictive initial Mass of the star (Msol)
 * @param       p_GBParams                      Hurley GB parameters
 * @param       p_McDU                          Core mass at second dredge up (Msol) (see Hurley et al. 2000, eq 69)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      TPAGB core mass of the star (Msol)
 */
GNU_CONST inline double TPAGB::CalculateCoreMass_Hurley2000(
    const double      p_Age,
    const double      p_MassEffectiveInitial,
    const DBL_VECTOR& p_GBParams,
    const double      p_McDU,
    const DBL_VECTOR& p_tScales
) const {

    const double lambda = std::min(0.9, 0.3 + (0.001 * utils::intPow(p_MassEffectiveInitial, 5))); // Hurley et al. 2000, eq 73
    // clamp core to maximum p_MassEffectiveInitial - core should not exceed total mass
    return std::min((p_McDU +  ((1.0 - lambda) * (CalculateCoreMassPrime_Hurley2000(p_Age, p_GBParams, p_tScales) - p_McDU))), p_MassEffectiveInitial);
}


/*
 * CalculateCoreMassPrime_Hurley2000
 *
 * @brief
 * Calculate M'c (core mass used for calculating luminosity and core mass on TPAGB),
 * per Hurley et al. 2000, just after eq 73.
 * 
 * Calculated using Hurley et al. 2000, eq 39, modified as described in Section 5.4
 *
 *
 * double CalculateCoreMassPrime_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_GBParams                      Hurley GB parameters
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      M'c (Msol)
 *
 * JR: todo: this is used in two places - work on calculating it only once.
 * I have done it this way to allow me to have a generic EvolveOneTimestep() function,
 * but there must be an elegant way of calculating once and using twice...
 */
GNU_CONST inline double TPAGB::CalculateCoreMassPrime_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const {
// #defines for convenience and readability - undefined at end of function
#define GBParams(x) p_GBParams[static_cast<int>(HURLEY_GB_PARAMETERS::x)]
#define tScales(x) p_tScales[static_cast<int>(HURLEY_TIMESCALES::x)]

    return p_Age <= tScales(Mx_SAGB)
            ? PPOW((GBParams(p) - 1.0) * GBParams(AHHe) * GBParams(D) * (tScales(Inf1_SAGB) - p_Age), 1.0 / (1.0 - GBParams(p)))
            : PPOW((GBParams(q) - 1.0) * GBParams(AHHe) * GBParams(B) * (tScales(Inf2_SAGB) - p_Age), 1.0 / (1.0 - GBParams(q)));

#undef tScales
#undef GBParams
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateRemnantRadius_Hurley2000_Static
 *
 * @brief
 * Calculate the radius of the remnant the star would become if it lost all of its envelope
 * immediately (i.e. M = Mc), per Hurley et al. 2000, at the end of section 6 (after eq 105).
 *
 *
 * static double CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass)
 *
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Remnant core radius (Rsol)
 */
double TPAGB::CalculateRemnantRadius() const {
    return HeWD::CalculateRadiusOnPhase_Static(m_CoreMass);
}






/// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//



