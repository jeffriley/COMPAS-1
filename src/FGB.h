#ifndef __FGB_h__
#define __FGB_h__

#include "constants.h"
//// #include "typedefs.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "profiling.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//// #include "utils.h"   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#include "HG.h"


class BaseStar;
class HG;

class FGB: virtual public BaseStar, public HG {

public:

    FGB() {};
    FGB(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), HG(p_BaseStar, false) { if (p_Initialise) Initialise(); }


protected:

    inline void Initialise() {
        CalculateTimescales(); // Initialise timescales
        m_Age = m_Timescales[static_cast<int>(TIMESCALE::tBGB)]; // Set age appropriately
        
        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically




///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

inline double CalculateCoreMass_Hurley2000() const override { 
    return CalculateCoreMass_Hurley2000(MassEffective(), Age(), Tau(), GBparams(), TimeScales());
}
COMPAS_PURE double CalculateCoreMass_Hurley2000(
    const double      p_Mass,
    const double      p_Age,
    const double      p_Tau,
    const DBL_VECTOR& p_GBparams,
    const DBL_VECTOR& p_tScales
) const;


inline double CalculateLuminosity_Hurley2000() const override { 
    return CalculateLuminosity_Hurley2000(Age(), GBparams(), TimeScales());
}
GNU_CONST inline double CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const;


GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; }            // McCO = 0.0 for FGB stars

inline double CalculateHeCoreMass() const override { return CalculateCoreMass(); } // McHe = Mc for FGB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed


inline double CalculateTau_Hurley2000() const override { return CalculateTau_Hurley2000(Age(), Timescales(TIMESCALE::tBGB), Timescales(TIMESCALE::tHeI)); }
GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const double p_tBGB, const double p_tHeI) const;


inline double CalculateEffectiveInitialMass_Hurley2000() const override { return GiantBranch::CalculateEffectiveInitialMass_Hurley2000(); } // per Hurley et al. 2000, section 7.1

 



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; }  // McCO = 0.0 for FGB stars

inline double CalculateHeCoreMassAtPhaseEnd() const override { return CalculateCoreMass(); } // McHe = Mc for FGB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed








    double          CalculateCoreMassAtPhaseEnd(const double p_Mass, const double p_Time) const     { return CalculateCoreMassOnPhase(p_Mass, p_Time); }                                            // Same as on phase
    double          CalculateCoreMassAtPhaseEnd() const                                             { return CalculateCoreMassAtPhaseEnd(m_Mass0, m_Age); }                                         // Use class member variables
    double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Time) const;
    double          CalculateCoreMassOnPhase() const                                                { return CalculateCoreMassOnPhase(m_Mass0, m_Age); }                                            // Use class member variables


    double          CalculateLuminosityAtPhaseEnd(const double p_Time) const                        { return CalculateLuminosityOnPhase(p_Time); }                                                  // Same as on phase
    double          CalculateLuminosityAtPhaseEnd() const                                           { return CalculateLuminosityAtPhaseEnd(m_Age); }                                                // Use class member variables
    double          CalculateLuminosityOnPhase(const double p_Time) const;
    double          CalculateLuminosityOnPhase() const                                              { return CalculateLuminosityOnPhase(m_Age); }                                                   // Use class member variables

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const { return GiantBranch::CalculateRadiusOnPhase(p_Mass, p_Luminosity); }                  // Skip HG - same as on phase
    double          CalculateRadiusAtPhaseEnd() const                                               { return CalculateRadiusAtPhaseEnd(m_Mass, m_Luminosity); }                                     // Use class member variables
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return GiantBranch::CalculateRadiusOnPhase(p_Mass, p_Luminosity); }
    double          CalculateRadiusOnPhase() const                                                  { return CalculateRadiusOnPhase(m_Mass, m_Luminosity); }                                        // Use class member variables







GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const override;


GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::CONVECTIVE; } // Always CONVECTIVE for FGB stars


    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                                            { return !ShouldEvolveOnPhase(); }                                                              // Phase ends when age at or after He ignition timescale
    bool            IsSupernova() const                                                             { return false; }                                                                               // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash();
    STELLAR_TYPE    ResolveSkippedPhase()                                                           { return STELLAR_TYPE::CORE_HELIUM_BURNING; }                                                   // Evolve to CHeB if phase is skipped

    bool            ShouldEvolveOnPhase() const                                                     { return (utils::Compare(m_Age, m_Timescales[static_cast<int>(TIMESCALE::tHeI)]) < 0); }        // Evolve on FGB phase if age < He ignition timescale
    bool            ShouldSkipPhase() const                                                         { return (utils::Compare(m_Mass0, m_MassCutoffs[static_cast<int>(MASS_CUTOFF::MFGB)]) >= 0); }  // Skip phase if mass >= FGB mass cutoff

    void            UpdateAgeAfterMassLoss()                                                        { GiantBranch::UpdateAgeAfterMassLoss(); }                                                      // Skip HG





};






//// inline candidates <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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
 * CalculateTau_Hurley2000
 *
 * @brief
 * Calculate the FGB-relative age (fractional First Giant Branch age) of the star,
 * per Hurley et al. 2000, just after eq 45.
 *
 * 
 * double CalculateTau_Hurley2000(const double p_Age, const double p_tBGB, const double p_tHeI) const
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (per Hurley timescales) (Myr)
 * @param       p_tHeI                          Time to helium ignition, tHeI (per Hurley timescales) (Myr)
 * @return                                      FGB-relative age, [0, 1]
 */
GNU_CONST inline double FGB::CalculateTau_Hurley2000(const double p_Age, const double p_tBGB, const double p_tHeI) const {
    return std::max(0.0, std::min(1.0, (p_Age - p_tBGB) / (p_tHeI - p_tBGB)));
}


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
 * @param       p_tScales                       Phase timescales (Myr)
 * @return                                      Suggested timestep (Myr)
 */
GNU_CONST inline double FGB::ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const {

    const double dtk = 0.02 * ((p_Age <= p_tScales[HURLEY_TS::MX_FGB] ? p_tScales[HURLEY_TS::INF1_FGB] : p_tScales[HURLEY_TS::INF2_FGB]) - p_Age); // stellar type specific dt
    const double dte = p_tScales[HURLEY_TS::HEI] - p_Age;           // time to end of phase (change of stellar type)

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);  // clamp to minimum NUCLEAR_MINIMUM_TIMESTEP
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
 * Calculate luminosity on the First Giant Branch,
 * per Hurley et al. 2000, eqs 37 & 39
 *
 *
 * double CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales)
 *
 * @param       p_Age                           Effective age of the star (Myr)
 * @param       p_GBparams                      Hurley GB parameters
 * @param       p_tScales                       Hurley timescales (Myr)
 * @return                                      FGB luminosity (Lsol)
 */
GNU_CONST inline double FGB::CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBparams, const DBL_VECTOR& p_tScales) const {

    const double AH = p_GBparams[HURLEY_GBP::AH];
    const double B  = p_GBparams[HURLEY_GBP::B];
    const double D  = p_GBparams[HURLEY_GBP::D];
    const double p  = p_GBparams[HURLEY_GBP::P];
    const double q  = p_GBparams[HURLEY_GBP::Q];

    // Calculate the core mass according to Hurley et al. 2000, eq 39, regardless
    // of whether it is the correct expression to use given the star's mass
    double McGB = p_Age < p_tScales[HURLEY_TS::MX_FGB] <= 0
                    ? PPOW(((p - 1.0) * AH * D * (p_tScales[HURLEY_TS::INF1_FGB] - p_Age)), (1.0 / (1.0 - p)))
                    : PPOW(((q - 1.0) * AH * B * (p_tScales[HURLEY_TS::INF2_FGB] - p_Age)), (1.0 / (1.0 - q)));

    return std::min((B * PPOW(McGB, q)), (D * PPOW(McGB, p))); // Hurley at al. 2000, eq 37
}





///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    
double          CalculateZetaEquilibrium()                                                      { return 0.0; }                                                     // At lowest order, giants with a convective envelope have radii that are insensitive to mass loss (but see Hurley+ 2002, Eq. 56 and Hurley+ 2000, Eq. 47)


double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const { return GiantBranch_Constituent::CalculateCriticalMassRatio_Claeys2014(p_AccretorIsDegenerate); } // Skip HG 

double CalculateCriticalMassRatio_Hurley2002() const { return GiantBranch_Constituent::CalculateCriticalMassRatio_Hurley2002(); }


#endif // __FGB_h__
