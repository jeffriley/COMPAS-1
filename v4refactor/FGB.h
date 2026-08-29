#pragma once

#include "constants.h"
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
        m_InterimState.SetAge(m_InterimState.HurleyTimescales(TS::BGB));
        
        EvolveOnPhase(0.0);
    }


    // member functions - alphabetically




///// ON PHASE FUNCTIONS   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




    ////////////////////////////////////////
    //   ABUNDANCE                        //
    ////////////////////////////////////////








    ////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES   //
    ////////////////////////////////////////

    inline double CalculateAgeAfterMassLoss_Hurley2000() const override { return GiantBranch::CalculateAgeAfterMassLoss_Hurley2000(); }




    GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const DblVectortT& p_tScales) const override; // JR FIX THIS: DONE
    GNU_CONST double CalculateTauAtPhaseEnd_Hurley2000(const double p_Age, const DblVectortT& p_tScales) const override { // JR FIX THIS: DONE
        return CalculateTau_Hurley2000(p_Age, p_tScales); // Same as on phase
    }



    ////////////////////////////////////////
    //   LUMINOSITY                       //
    ////////////////////////////////////////



inline double CalculateLuminosity_Hurley2000() const override { 
    return CalculateLuminosity_Hurley2000(Age(), m_InterimState.HurleyGBParams(), m_InterimState.HurleyTimescales());
}
GNU_CONST inline double CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const;






    ////////////////////////////////////////
    //   MASS                             //
    ////////////////////////////////////////

    GNU_PURE inline double CalculateCoreMass_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateCoreMass_Hurley2000(Metallicity(), MassEffectiveInitial(), Age(), Tau(), GBParams(), Timescales());
    }
    GNU_PURE double CalculateCoreMass_Hurley2000( // JR FIX THIS: DONE
        const double      p_Metallicity,
        const double      p_Mass,
        const double      p_Age,
        const double      p_Tau,
        const DBL_VECTOR& p_GBParams,
        const DBL_VECTOR& p_tScales
    ) const;

    GNU_PURE inline double CalculateCoreMassAtPhaseEnd_Hurley2000() const override { return CalculateCoreMass_Hurley2000(); } // Same as on phase JR FIX THIS: DONE





GNU_CONST inline double CalculateCOCoreMass() const override { return 0.0; }            // McCO = 0.0 for FGB stars

inline double CalculateHeCoreMass() const override { return CalculateCoreMass(); } // McHe = Mc for FGB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed




inline double CalculateEffectiveInitialMass_Hurley2000() const override { return GiantBranch::CalculateEffectiveInitialMass_Hurley2000(); } // per Hurley et al. 2000, section 7.1

 

    double          CalculateCoreMassAtPhaseEnd(const double p_Mass, const double p_Time) const     { return CalculateCoreMassOnPhase(p_Mass, p_Time); }                                            // Same as on phase
    double          CalculateCoreMassAtPhaseEnd() const                                             { return CalculateCoreMassAtPhaseEnd(Mass0(), Age()); }
    double          CalculateCoreMassOnPhase(const double p_Mass, const double p_Time) const {  } // JR FIX THIS: placeholder
    double          CalculateCoreMassOnPhase() const                                                { return CalculateCoreMassOnPhase(Mass0(), Age()); }




///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

GNU_CONST inline double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; }  // McCO = 0.0 for FGB stars

inline double CalculateHeCoreMassAtPhaseEnd() const override { return CalculateCoreMass(); } // McHe = Mc for FGB stars  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< check for already computed








    double          CalculateLuminosityAtPhaseEnd(const double p_Time) const                        { return CalculateLuminosityOnPhase(p_Time); }                                                  // Same as on phase
    double          CalculateLuminosityAtPhaseEnd() const                                           { return CalculateLuminosityAtPhaseEnd(Age()); }
    double          CalculateLuminosityOnPhase(const double p_Time) const {  } // JR FIX THIS: placeholder
    double          CalculateLuminosityOnPhase() const                                              { return CalculateLuminosityOnPhase(Age()); }

    double          CalculateRadiusAtPhaseEnd(const double p_Mass, const double p_Luminosity) const { return GiantBranch::CalculateRadius(p_Mass, p_Luminosity); }                  // Skip HG - same as on phase
    double          CalculateRadiusAtPhaseEnd() const                                               { return CalculateRadiusAtPhaseEnd(Mass(), Luminosity()); }
    double          CalculateRadiusOnPhase(const double p_Mass, const double p_Luminosity) const    { return GiantBranch::CalculateRadius(p_Mass, p_Luminosity); }
    double          CalculateRadiusOnPhase() const                                                  { return CalculateRadiusOnPhase(Mass(), Luminosity()); }







GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const override;


GNU_CONST inline ENVELOPE DetermineEnvelopeType() const override { return ENVELOPE::CONVECTIVE; } // Always CONVECTIVE for FGB stars


    STELLAR_TYPE    EvolveToNextPhase();

    bool            IsEndOfPhase() const                                                            { return !ShouldEvolveOnPhase(); }                                                              // Phase ends when age at or after He ignition timescale
    bool            IsSupernova() const                                                             { return false; }                                                                               // Not here

    STELLAR_TYPE    ResolveEnvelopeLoss(bool p_Force = false);
    void            ResolveHeliumFlash();
    STELLAR_TYPE    ResolveSkippedPhase()                                                           { return STELLAR_TYPE::CORE_HELIUM_BURNING; }                                                   // Evolve to CHeB if phase is skipped

    bool            ShouldEvolveOnPhase() const                                                     { return (utils::Compare(Age(), m_InterimState.HurleyTimescales(TS::HeI)) < 0); }
    bool            ShouldSkipPhase() const                                                         { return (utils::Compare(Mass0(), ZDEP->HurleyMassCutoffs(Metallicity(), MCO::FGB)) >= 0); }

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
 * double CalculateTau_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      FGB-relative age of the star, [0, 1]
 */
inline double FGB::CalculateTau_Hurley2000(const double p_Age, const DblVectorT& p_tScales) const {
    const double tBGB = p_tScales[TIMESCALES::BGB];
    return std::clamp((p_Age - tBGB) / (p_tScales[TIMESCALES::HEI] - tBGB), 0.0, 1.0);
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
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Suggested timestep (Myr)
 */
GNU_CONST inline double FGB::ChooseTimestep_Hurley2000(const double p_Age, const DBL_VECTOR& p_tScales) const {

    const double dtk = 0.02 * ((p_Age <= p_tScales[TIMESCALES::MX_FGB] ? p_tScales[TIMESCALES::INF1_FGB] : p_tScales[TIMESCALES::INF2_FGB]) - p_Age); // stellar type specific dt
    const double dte = p_tScales[TIMESCALES::HEI] - p_Age;           // time to end of phase (change of stellar type)

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
 * double CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales)
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_GBParams                      Hurley GB parameters
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      FGB luminosity (Lsol)
 */
GNU_CONST inline double FGB::CalculateLuminosity_Hurley2000(const double p_Age, const DBL_VECTOR& p_GBParams, const DBL_VECTOR& p_tScales) const {

    const double AH = p_GBParams[HURLEY_GBP::AH];
    const double B  = p_GBParams[HURLEY_GBP::B];
    const double D  = p_GBParams[HURLEY_GBP::D];
    const double p  = p_GBParams[HURLEY_GBP::P];
    const double q  = p_GBParams[HURLEY_GBP::Q];

    // Calculate the core mass according to Hurley et al. 2000, eq 39, regardless
    // of whether it is the correct expression to use given the star's mass
    double McGB = ((p_Age < p_tScales[TS::MX_FGB]) <= 0)
                    ? PPOW(((p - 1.0) * AH * D * (p_tScales[TIMESCALES::INF1_FGB] - p_Age)), (1.0 / (1.0 - p)))
                    : PPOW(((q - 1.0) * AH * B * (p_tScales[TIMESCALES::INF2_FGB] - p_Age)), (1.0 / (1.0 - q)));

    return std::min((B * PPOW(McGB, q)), (D * PPOW(McGB, p))); // Hurley at al. 2000, eq 37
}

















///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    double CalculateZetaEquilibrium() override { return 0.0; }                                          // At lowest order, giants with a convective envelope have radii that are insensitive to mass loss (but see Hurley+ 2002, Eq. 56 and Hurley+ 2000, Eq. 47)

    // Skip HG's HG-specific qcrit -- delegate explicitly to GiantBranch_Constituent
    double CalculateCriticalMassRatio_Claeys2014(const bool p_AccretorIsDegenerate) const override { return GiantBranch_Constituent::CalculateCriticalMassRatio_Claeys2014(p_AccretorIsDegenerate); }
    double CalculateCriticalMassRatio_Hurley2002() const override { return GiantBranch_Constituent::CalculateCriticalMassRatio_Hurley2002(); }

