#ifndef __TPAGB_h__
#define __TPAGB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "EAGB.h"


class BaseStar;
class EAGB;

class TPAGB: virtual public BaseStar, public EAGB {

public:

    TPAGB() { m_StellarType = STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH; };
    
    TPAGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), EAGB(p_BaseStar, false) {
        m_StellarType = STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH;                                                                                                                            // Set stellar type
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
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tP)];                                                                                                                                              // Set age appropriately
        
        EvolveOnPhase(0.0);
   }


   // member functions - alphabetically




            double          CalculateConvectiveCoreRadius() const                                                   { return std::min(5.0 * CalculateRemnantRadius(), m_Radius); }       // Last paragraph of section 6 of Hurley+ 2000

            double          CalculateCoreMassAtPhaseEnd() const                                                     { return m_CoreMass; }                                                                  // NO-OP
            double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Time) const;
            double          CalculateCoreMassOnPhase() const                                                        { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                    // Use class member variables

                                                                // NO-OP

            double          CalculateCELambda_Dewi() const;
            double          CalculateLambdaNanjingStarTrack(const double p_Mass) const;
            double          CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;

            double          CalculateLuminosityOnPhase(const double p_Time) const;
            double          CalculateLuminosityOnPhase() const                                                      { return CalculateLuminosityOnPhase(m_Age); }                                           // Use class member variables
            double          CalculateLuminosityAtPhaseEnd() const                                                   { return m_Luminosity; }                                                                // NO-OP


            double          CalculateRadiusAtPhaseEnd() const                                                       { return m_Radius; }                                                                    // NO-OP
            double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const            { return CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, m_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)], m_BnCoefficients); }
            double          CalculateRadiusOnPhase() const                                                          { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }                                // Use class member variables
    static  double          CalculateRadiusOnPhase_Static(const double      p_Mass,
                                                          const double      p_Luminosity,
                                                          const double      p_MHeF,
                                                          const DBL_VECTOR &p_BnCoefficients);

            double          CalculateRemnantLuminosity() const;













///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<       


DBL_VECTOR CalculateTimescales_Hurley2000() const override {
    return CalculateTimescales_Hurley2000(
        m_StateHistory.CurrentState.MassEffectiveInitial(),
        m_StateHistory.CurrentState.GBparams(),
        m_StateHistory.CurrentState.TimeScales()
    );
}
COMPAS_PURE DBL_VECTOR CalculateTimescales_Hurley2000(const double p_Mass, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;

inline DBL_DBL CalculateConvectiveEnvelopeMass() const override {
    return CalculateConvectiveEnvelopeMass(m_StateHistory.CurrentState.Mass(), m_StateHistory.CurrentState.CoreMass());
}
GNU_CONST DBL_DBL CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_CoreMass) const;


        
double CalculateCoreMass() const override;
GNU_CONST double CalculateCoreMass_Hurley2000(const double p_MassEffectiveInitial, const double p_Age, const double p_McDU) const;

inline double CalculateCOCoreMass() const override { return CalculateCoreMass(); } // McCO = Mc for TPAGB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed
inline double CalculateHeCoreMass() const override { return CalculateCoreMass(); } // McHe = Mc for TPAGB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed

double CalculateTau_Hurley2000() const override { return 0.0; } // Tau (relative age) is not used for TPAGB stars in Hurley et al. 2000, so we return 0.0

GNU_CONST double CalculateMcPrime_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;


inine double CalculateRemnantRadius_Hurley2000() const override {
    return CalculateRemnantRadius_Hurley2000_Static(m_StateHistory.CurrentState.CoreMass());
}
GNU_CONST static double CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass);






///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inline double CalculateCOCoreMassAtPhaseEnd() const override {
    return std::min(CalculateCoreMass(), m_StateHistory.CurrentState.Mass()); // McCO = Mc for TPAGB stars; McCO should be <= M  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed
}

inline double CalculateHeCoreMassAtPhaseEnd() const override { return CalculateCoreMass(); } // McHe = Mc for TPAGB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed




            double          CalculateTemperatureAtPhaseEnd(const double p_Luminosity, const double p_Radius) const  { return m_Temperature; }                                                               // NO-OP
            double          CalculateTemperatureAtPhaseEnd() const                                                  { return CalculateTemperatureAtPhaseEnd(m_Luminosity, m_Radius); }                      // Use class member variables



            double          ChooseTimestep(const double p_Time) const;

            ENVELOPE        DetermineEnvelopeType() const                                                           { return ENVELOPE::CONVECTIVE; }                                                        // Always CONVECTIVE

            STELLAR_TYPE    EvolveToNextPhase()                                                                     { return m_StellarType; }                                                               // NO-OP

            bool            IsEndOfPhase() const                                                                    { return !ShouldEvolveOnPhase(); }                                                      // Phase ends when envelope loss or going supernova
            bool            IsSupernova() const;

            STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
            void            ResolveHeliumFlash() { }                                                                                                                                                        // NO-OP
            STELLAR_TYPE    ResolveSkippedPhase()                                                                   { return m_StellarType; }                                                               // NO-OP

            bool            ShouldEvolveOnPhase() const                                                             { return (utils::Compare(m_COCoreMass, m_Mass) < 0 && !IsSupernova() && !ShouldEnvelopeBeExpelledByPulsations()); }                        // Evolve on TPAGB phase if envelope is not lost and not going supernova
            bool            ShouldSkipPhase() const                                                                 { return false; }                                                                       // Never skip TPAGB phase

};



//// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  MASS FUNCTIONS                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

      
/*
 * CalculateConvectiveEnvelopeMass
 *
 * @brief
 * Calculate the convective envelope mass for Thermally Pulsing Asymptotic Giant Bramch, TPAGB, stars.
 * We assume entire envelope is convective for TPAGB stars.
 *
 * DBL_DBL CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_CoreMass) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @return                                      Tuple containing:
 *                                                   TPAGB convective envelope mass (Msol)
 *                                                   Maximum TPAGB convective envelope mass (Msol)
 */
GNU_CONST inline DBL_DBL CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_CoreMass) const {
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

    Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

        EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
        EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            coreMass = CalculateCoreMass_Hurley2000(m_StateHistory.CurrentState.MassEffectiveInitial(),
                                                    m_StateHistory.CurrentState.Age(),
                                                    m_StateHistory.CurrentState.GBparams(McDU));
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
 * Calculate core mass on the Thermally Pulsing Asymptotic Giant Branch,
 * per Hurley et al. 2000, just after eq 73
 *
 *
 * double CalculateCoreMass_Hurley2000(
 *     const double      p_Age,
 *     const double      p_MassEffectiveInitial,
 *     const DBL_VECTOR& p_GBparams,
 *     const double      p_McDU,
 *     const DBL_VECTOR& p_tScales
 * ) const
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_MassEffectiveInitial          Effictive initial Mass of the star (Msol)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_McDU                          Core mass at second dredge up (Msol) (see Hurley et al. 2000, eq 69)
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      Core mass on the Thermally Pulsing Asymptotic Giant Branch in Msol
 */
GNU_CONST inline double TPAGB::CalculateCoreMass_Hurley2000(
    const double      p_Age,
    const double      p_MassEffectiveInitial,
    const DBL_VECTOR& p_GBparams,
    const double      p_McDU,
    const DBL_VECTOR& p_tScales
) const {

    const double lambda = std::min(0.9, 0.3 + (0.001 * utils::intPow(p_MassEffectiveInitial, 5))); // Hurley et al. 2000, eq 73
    // Clam core to maximum p_MassEffectiveInitial - core should not exceed total mass
    return std::min((p_McDU +  ((1.0 - lambda) * (CalculateCoreMassPrime_Hurley2000(p_Age, p_GBparams, p_tScales) - p_McDU))), p_MassEffectiveInitial);
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
 * double CalculateCoreMassPrime_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      M'c (Msol)
 *
 * JR: todo: this is used in two places - work on calculating it only once.
 * I have done it this way to allow me to have a generic EvolveOneTimestep() function,
 * but there must be an elegant way of calculating once and using twice...
 */
GNU_CONST inline double TPAGB::CalculateCoreMassPrime_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const {
// #defines for convenience and readability - undefined at end of function
#define GBparams(x) p_GBparams[static_cast<int>(HURLEY_GBP:::x)]
#define tScales(x) tScales[static_cast<int>(TIMESCALE::x)]

    return p_Age <= tScales(tMx_SAGB)
            ? PPOW((GBparams(p) - 1.0) * GBparams(AHHe) * GBparams(D) * (tScales(tinf1_SAGB) - p_Age), 1.0 / (1.0 - GBparams(p)))
            : PPOW((GBparams(q) - 1.0) * GBparams(AHHe) * GBparams(B) * (tScales(tinf2_SAGB) - p_Age), 1.0 / (1.0 - GBparams(q)));

#undef tScales
#undef GBparams
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                             LIFETIME / AGE FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               LUMINOSITY FUNCTIONS                                //
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
 * double CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      TPAGB luminosity (Lsol)
 */
GNU_CONST inline double TPAGB::CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const {
    return CalculateLuminosityGivenCoreMass_Hurley2000(CalculateCoreMassPrime_Hurley2000(p_Age, p_GBparams, p_tScales), p_GBparams);
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
double TPAGB::CalculateRemnantLuminosity() const {
    return COWD::CalculateLuminosityOnPhase_Static(m_CoreMass, 0.0, m_Metallicity);
}




///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 RADIUS FUNCTIONS                                  //
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
GNU_CONST static inline double TPAGB::CalculateRemnantRadius_Hurley2000_Static(const double p_CoreMass) {
    return WhiteDwarfs::CalculateRadiusOnPhase_Marsh2004_Static(p_CoreMass);
}


/// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


#endif // __TPAGB_h__
