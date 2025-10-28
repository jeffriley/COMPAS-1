#ifndef __EAGB_h__
#define __EAGB_h__

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "CHeB.h"


class BaseStar;
class CHeB;

class EAGB: virtual public BaseStar, public CHeB {

public:

    EAGB() { m_StellarType = STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH; };
    
    EAGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), CHeB(p_BaseStar, false) {
        m_StellarType = STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH;                                                                                                    // Set stellar type
        if (p_Initialise) Initialise();                                                                                                                                 // Initialise if required
    }

    EAGB* Clone(const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        EAGB* clone = new EAGB(*this, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }

    static EAGB* Clone(EAGB& p_Star, const OBJECT_PERSISTENCE p_Persistence, const bool p_Initialise = true) {
        EAGB* clone = new EAGB(p_Star, p_Initialise); 
        clone->SetPersistence(p_Persistence); 
        return clone; 
    }


    // member functions
    static double   CalculateRadiusOnPhase_Static(const double      p_Mass,
                                                  const double      p_Luminosity,
                                                  const double      p_MHeF,
                                                  const DBL_VECTOR &p_BnCoefficients);


protected:

    void Initialise() {
        CalculateTimescales();                                                                                                                                          // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tHeI)] + m_Timescales[static_cast<int>(TIMESCALE::tHe)];                                                       // Set age appropriately
        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically






///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inline double CalculateCOCoreMass_Hurley2000() const override {
    return CalculateCOCoreMass_Hurley2000(m_StateHistory.CurrentState.Age(),
                                          m_StateHistory.CurrentState.GBparams(),
                                          m_StateHistory.CurrentState.Timescales());
}
GNU_CONST double CalculateCOCoreMass_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;


inline double CalculateLuminosityOnPhase_Hurley2000() const override {
    // use Hurley core mass - luminosity relationship, per Hurley et al. eq 37
    return CalculateLuminosityGivenCoreMass_Hurley2000(m_StateHistory.CurrentState.CoreMass(), m_StateHistory.CurrentState.GBparams());
}

double CalculateTau_Hurley2000() const override { return 0.0; }; // Tau (relative age) is not used for EAGB stars in Hurley et al. 2000, so we return 0.0







///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd_Hurley2000() const override {
    return m_StateHistory.CurrentState.GBparams(McDU); }    // McCO = McDU at phase end for EAGB stars, per Hurley et al. 2000, section 5.4
}



    double          CalculateCoreMassAtPhaseEnd() const                                             { return m_GBparams[static_cast<int>(GBP::McDU)]; }                 // Mc(EAGB) = McDU at phase end (Hurley et al. 2000, section 5.4)
    double          CalculateCoreMassOnPhase() const                                                { return m_GBparams[static_cast<int>(GBP::McBAGB)]; }               // Mc(EAGB) = McHe(EAGB) = McBAGB on phase (Hurley et al. 2000, section 5.4)

    double          CalculateHeCoreMassAtPhaseEnd() const                                           { return CalculateHeCoreMassOnPhase(); }                            // Same as on phase
    double          CalculateHeCoreMassOnPhase() const                                              { return m_HeCoreMass; }                                            // NO-OP

    double          CalculateInitialSupernovaMass() const                                           { return m_GBparams[static_cast<int>(GBP::McBAGB)]; }               // For EAGB & TPAGB we use the mass at Base Asymptotic Giant Branch to determine SN type

    double          CalculateLambdaNanjingStarTrack(const double p_Mass, const double p_Metallicity) const;
    double          CalculateLambdaNanjingEnhanced(const int p_MassIndex, const STELLAR_POPULATION p_StellarPop) const;

    double          CalculateLifetimeTo2ndDredgeUp(const double p_Tinf1_FAGB, const double p_Tinf2_FAGB) const;

    double          CalculateLuminosityAtPhaseEnd(const double p_CoreMass) const                    { return CalculateLuminosityOnPhase(p_CoreMass); }                  // Same as on phase
    double          CalculateLuminosityAtPhaseEnd() const                                           { return CalculateLuminosityAtPhaseEnd(m_COCoreMass); }             // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_CoreMass) const;
    double          CalculateLuminosityOnPhase() const                                              { return CalculateLuminosityOnPhase(m_COCoreMass); }



    COMPAS_PURE MASS_LOSS_T CalculateMLrate_Hurley2000(const double p_Mass, const double p_Radius, const double p_Luminosity, const double p_PerturbationMu) const override;



    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const { return CalculateRadiusOnPhase(p_Mass, p_Luminosity); }            // Same as on phase
    double          CalculateRadiusAtPhaseEnd() const                                               { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }         // Use class member variables
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return CalculateRadiusOnPhase_Static(p_Mass, p_Luminosity, m_MassCutoffs[static_cast<int>(MASS_CUTOFF::MHeF)], m_BnCoefficients); }
    double          CalculateRadiusOnPhase() const                                                  { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }            // Use class member variables

    double          CalculateRemnantLuminosity() const;
    double          CalculateRemnantRadius() const;






    void            CalculateTimescales(const double p_Mass, DBL_VECTOR &p_Timescales);
    void            CalculateTimescales()                                                           { CalculateTimescales(m_Mass0, m_Timescales); }                     // Use class member variables

    double          ChooseTimestep(const double p_Time) const;

    ENVELOPE        DetermineEnvelopeType() const                                                   { return ENVELOPE::CONVECTIVE; }                                    // Always CONVECTIVE

    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                                            { return !ShouldEvolveOnPhase(); }                                  // Phase ends when age at or after DU timescale, and no TPAGB
    bool            IsSupernova() const;

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash() {  }                                                                                                                           // NO-OP
    STELLAR_TYPE    ResolveSkippedPhase()                                                           { return m_StellarType; }                                           // NO-OP

    bool            ShouldEvolveOnPhase() const;
    bool            ShouldSkipPhase() const;

};




///////////// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                             LIFETIME / AGE FUNCTIONS                              //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  MASS FUNCTIONS                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateCOCoreMass_Hurley2000
 *
 * @brief
 * Calculate CO core mass on the Early Asymptotic Giant Branch,
 * per Hurley et al. 2000, eq 39, modified as described in Section 5.4
 *
 *
 * double CalculateCOCoreMass_Hurley2000(const double p_Time, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales)
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      EAGB CO core mass (Msol)
 */
GNU_CONST inline double EAGB::CalculateCOCoreMass_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const {
// #defines for convenience and readability - undefined at end of function
#define GBparams(x) p_GBparams[static_cast<int>(GBP::x)]
#define tScales(x) p_tScales[static_cast<int>(TIMESCALE::x)]

    return p_Age <= tScales(tMx_FAGB)
            ? PPOW((GBparams(p) - 1.0) * GBparams(AHe) * GBparams(D) * (tScales(tinf1_FAGB) - p_Age), 1.0 / (1.0 - GBparams(p)))
            : PPOW((GBparams(q) - 1.0) * GBparams(AHe) * GBparams(B) * (tScales(tinf2_FAGB) - p_Age), 1.0 / (1.0 - GBparams(q)));

#undef tScales
#undef GBparams
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                 RADIUS FUNCTIONS                                  //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * Calculate radius of the remnant the star would become if it lost all of its
 * envelope immediately (i.e. M = Mc, coreMass)
 *
 * Hurley et al. 2000, just after eq 105
 *
 *
 * double CalculateRemnantRadius()
 *
 * @return                                      Radius of remnant core in Rsol
 */
double EAGB::CalculateRemnantRadius() const {
    double R1, R2;
    std::tie(R1, R2) = HeGB::CalculateRadiusOnPhase_Static(m_HeCoreMass, CalculateRemnantLuminosity());

    return std::min(m_Radius, std::min(R1, R2));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               LUMINOSITY FUNCTIONS                                //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////





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
double EAGB::CalculateRemnantLuminosity() const {
    return HeGB::CalculateLuminosityOnPhase_Static(m_COCoreMass, gbParams(B), gbParams(D));
}




///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    MT_CASE         DetermineMassTransferTypeAsDonor() const                                        { return MT_CASE::C; }                                              // Always case C

double CalculateCriticalMassRatio_Hurley2002() const { return GiantBranch_Constituent::CalculateCriticalMassRatio_Hurley2002(); }



#endif // __EAGB_h__
