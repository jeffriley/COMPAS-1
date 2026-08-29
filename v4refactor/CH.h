#pragma once

#include "constants.h"
#include "BaseStar.h"
#include "MS_gt_07.h"

class BaseStar;
class MS_gt_07;

class CH: virtual public BaseStar, public MS_gt_07 {


public:

    // Constructors
    CH() {};
    CH(const BaseStar &p_BaseStar, const bool p_Initialise = true) : BaseStar(p_BaseStar), MS_gt_07(p_BaseStar) {
        if (p_Initialise) Initialise();
    }


private: 

    // Initialise
    inline void Initialise() {
        // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS --- CHECK THIS IS STILL OK <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        CalculateTimescales();                       // Initialise timescales
        m_InterimState.SetAge(0.0);
        m_InterimState.SetCHE(true);
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
    //   ABUNDANCE                                  //
    //////////////////////////////////////////////////

    GNU_CONST double CalculateHAbundanceCore(const double p_Tau, const double p_ZAMSHAbundance) const; // JR FIX THIS: DONE
    GNU_PURE  double CalculateHAbundanceCore() const override { // JR FIX THIS: DONE
        return CalculateHAbundanceCore(Tau(), ZDEP->ZAMSHAbundance(Metallicty()));
    }

    GNU_PURE double CalculateHAbundanceSurface(const double p_Tau) const override { // JR FIX THIS: DONE
        return CalculateHAbundanceCore(p_Tau, ZDEP->ZAMSHAbundance(Metallicty()));
    }


    GNU_CONST double CalculateHeAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_ZAMSHeAbundance) const; // JR FIX THIS: DONE
    GNU_PURE  double CalculateHeAbundanceCore() const override { // JR FIX THIS: DONE
        return CalculateHeAbundanceCore(Metallicty(), Tau(), ZDEP->ZAMSHeAbundance(Metallicty()));
    }

    inline double CalculateHeAbundanceSurface(const double p_Tau) const override { // JR FIX THIS: DONE
        return CalculateHeAbundanceCore(Metallicty(), p_Tau, ZDEP->ZAMSHeAbundance(Metallicty()));
    }


    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////
    
    GNU_PURE double CalculateAgeAfterMassLoss_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const double p_tMS) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateAgeAfterMassLoss_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateAgeAfterMassLoss_Hurley2000(Metallicity(), Mass(), Age(), Timescales(TS::MS));
    }


    GNU_CONST double CalculateLifetimesRatio_Szecsi2020(const double p_Mass) const; // JR FIX THIS: DONE


    GNU_PURE TimescalesT CalculateTimescales_Hurley2000(const double p_Metallicity, const double p_Mass, const TimescalesT& p_tScales) const override; // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////

    GNU_CONST double CalculateLogLuminositiesRatio_Szecsi2020(const double p_Mass) const; // JR FIX THIS: DONE

    inline double CalculateLuminosity_Hurley2000() const override {
        return CalculateLuminosity_Hurley2000(
            Mass(),
            Tau(),
            Age(),
            m_StateHistory.ZAMSStateRef().Luminosity(),  // <<<<<<<<<<<<<<<<<<< JR FIX THIS <<<<<<<<<<<<<<<<<<<<<<
            Timescales(TS::MS),
            Timescales(TS::BGB)
        );
    }
    COMPAS_PURE double CalculateLuminosity_Hurley2000(const double p_Mass, const double p_Tau, const double p_Age, const double p_LZAMS, const double p_tMS, const double p_tBGB) const;

    inline double CalculateLuminosityAtPhaseEnd() const override { return CalculateLuminosityAtPhaseEnd(Mass0()); }                                                          // Use class member variables
    GNU_PURE double CalculateLuminosityAtPhaseEnd(const double p_Mass) const;






    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////







    //////////////////////////////////////////////////
    //   MASS LOSS / ACCRETION                      //
    //////////////////////////////////////////////////



    // mass loss rate
    COMPAS_PURE MassLossT CalculateMLRate_Belczynski2010(
        const double p_Metallicity,
        const double p_Mass,
        const double p_Radius,
        const double p_Luminosity,
        const double p_Temperature,
        const double p_PerturbationMu,
        const double p_HeAbundanceSurface
    ) const override;

    GNU_CONST   double CalculateMLRateRotationEnhancement_Langer1998() const;
    GNU_CONST   double CalculateMLfractionOB(const double p_HeAbundanceSurface) const;

    GNU_PURE MassLossT CalculateMLRate_Merritt2025(
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
    ) const override;



    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////


    inline double CalculateRadius() const override { return ZAMSState().Radius(); }           // CH radius is fixed throughout lifetime
    inline double CalculateRadiusAtPhaseEnd() const override { return CalculateRadius(); }  // CH radius is fixed throughout lifetime 













    STELLAR_TYPE    EvolveToNextPhase() override;

//    bool            ShouldEvolveOnPhase() const                         { return m_Age < m_Timescales[static_cast<int>(TS::MS)] && (OPTIONS->OptimisticCHE() || Omega() >= m_OmegaCHE); } // Evolve on CHE phase if age in MS timescale and spinning at least as fast as CHE threshold

};









///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                         INLINE CANDIDATE IMPLEMENTATIONS                          //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ABUNDANCE                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateHAbundanceCore
 *
 * @brief
 * Calculate the hydrogen abundance in the core of a CH star, given the phase-relative age
 * of the star.  Currently just a simple linear model that assumes that hydrogen in the
 * core of the star is burned to helium at a constant rate throughout the star's lifetime.
 * 
 * Should one day be updated to match detailed models.
 * 
 *
 * double CalculateHAbundanceCore(const double p_Metallicity, const double p_Tau) const
 * 
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_ZAMSHAbundance                ZAMS hydrogen abundance of the star [0, 1]
 * @return                                      Hydrogen abundance in the core of the star [0, 1]
 */
inline double CH::CalculateHAbundanceCore(const double p_Tau, const double p_ZAMSHAbundance) const {
    return p_ZAMSHAbundance * (1.0 - p_Tau);
}


/*
 * CalculateHeAbundanceCore
 *
 * @brief
 * Calculate the helium abundance in the core of a CH star, given the phase-relative age
 * of the star.  Currently just a simple linear model from the initial helium abundance
 * to the maximum helium abundance (assuming that all hydrogen is converted to helium). 
 * 
 * Should one day be updated to match detailed models.
 * 
 *
 * double CalculateHeAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_ZAMSHeAbundance) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_ZAMSHeAbundance               ZAMS helium abundance of the star [0, 1]
 * @return                                      Helium abundance in the core of the star [0, 1]
 */
inline double CH::CalculateHeAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_ZAMSHeAbundance) const {
    return ((1.0 - p_Metallicity - p_ZAMSHeAbundance) * p_Tau) + p_ZAMSHeAbundance;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                    AGE / LIFETIME / TAU / TIMESCALES / TIMESTEP                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateAgeAfterMassLoss_Hurley2000
 *
 * @brief
 * Recalculate the age of a CH star after mass loss, per Hurley et al. 2000, section 7.1.
 * 
 * The age will be further modified if option `--enhance-CHE-lifetimes-luminosities`
 * was specified.
 *
 *
 * double CalculateAgeAfterMassLoss_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const double p_tMS) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tMS                           MS lifetime, tMS (Myr)
 * @return                                      Age of the star after mass loss (Myr)
 */
inline double CH::CalculateAgeAfterMassLoss_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const double p_tMS) const {
    // We call astro::CalculateLifetimeToBGB_Hurley2000() here rather than use timescales
    // because p_Mass may not be the same mass used to calculate timescales[TS::BGB]
    const double tBGBprime = astro::CalculateLifetimeToBGB_Hurley2000(p_Metallicity, p_Mass);
    double age = p_Age * MainSequence::CalculatePhaseLifetime_Hurley2000(p_Metallicity, p_Mass, tBGBprime) / p_tMS;   

    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {       // Enhance lifetime of CH stars?
        age *= CalculateLifetimesRatio_Szecsi2020(p_Mass);  // Yes
    }

    return age;
}


/*
 * CalculateLifetimesRatio_Szecsi2020
 *
 * @brief
 * Calculate the ratio of the lifetimes of a CH star and a non-CH MS star of the same mass
 * (t_CHE/t_MS).  Uses a polynomial fit derived using IZw18 and IZw18CHE BoOST models from
 * Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203).
 *
 * 
 * double CalculateLifetimesRatio_Szecsi2020(const double p_Mass) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Lifetimes ratio
 */
inline double CH::CalculateLifetimesRatio_Szecsi2020(const double p_Mass) const {

    const double x  = std::log10(p_Mass);
    const double x2 = x * x;

    return PPOW(10.0, -0.15929168474199387 + (1.050069750483549 * x) + (-0.8233601359988406 * x2) + (0.17772610259473764 * x * x2));
}


/*
 * CalculateTimescales_Hurley2000
 *
 * @brief
 * Calculate CH timescales.
 * CH timescales are MS timescales, per Hurley at al. 2000, optionally enhanced if the user
 * specified the `--enhance-CHE-lifetimes-luminosities` option.
 * 
 * Since timescales depend on a star's mass, they need to be calculated whenever the mass of the
 * star changes (probably every timestep).
 *
 * p_tScales is copied once and returned by value.  NRVO constructs the returned array directly
 * in the caller's return slot, so no copy on the return.
 *
 *
 * DblVectorT CalculateTimescales_Hurley2000(const double p_Metallicity, const double p_Mass, const TimescaleT& p_tScales) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
inline DblVectorT CH::CalculateTimescales_Hurley2000(const double p_Metallicity, const double p_Mass, const TimescaleT& p_tScales) const {

    TimescaleT<TIMESCALES_COUNT> tScales = std::move(p_tScales);

    tScales[TS::BGB] = astro::CalculateLifetimeToBGB_Hurley2000(p_Metallicity, p_Mass);
    tScales[TS::MS]  = MainSequence::CalculatePhaseLifetime_Hurley2000(p_Metallicity, p_Mass, tScales[TS::BGB]);

    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {   // Enhance lifetime of CH stars?
                                                        // Yes
        const double lifetimesRatio = CalculateLifetimesRatio_Szecsi2020(p_Mass);

        tScales[TS::BGB] *= lifetimesRatio;
        tScales[TS::MS]  *= lifetimesRatio;
    }

    return tScales;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                               MASS LOSS / ACCRETION                               //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateMLfractionOB
 *
 * @brief
 * Calculate the fraction of mass loss attributable to OB mass loss, per Yoon et al. 2006
 *
 * The model described in Yoon et al. 2006 (also Szecsi et al. 2015) uses OB mass loss while the
 * He surface abundance is below 0.55, WR mass loss when the surface He abundance is above 0.7,
 * and linearly interpolate when the He surface abundance is between those limits.
 *
 * This function calculates the fraction of mass loss attributable to OB mass loss, based on
 * the He surface abundance and the abundance limits described in Yoon et al. 2006.  The value
 * returned will be 1.0 if 100% of the mass loss is attributable to OB mass lass, 0.0 if 100% of
 * the mass loss is attributable to WR mass loss, and in the range (0.0, 1.0) if the mass loss is
 * a mix of OB and WR.
 *
 *
 * double CalculateMLfractionOB(const double p_HeAbundanceSurface) const
 *
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Fraction of mass loss attributable to OB mass loss
 */
GNU_CONST inline double CH::CalculateMLfractionOB(const double p_HeAbundanceSurface) const {

    constexpr double limOB = 0.55;  // per Yoon et al. 2006
    constexpr double limWR = 0.70;  // per Yoon et al. 2006

    return std::min(1.0, std::max (0.0, (limWR - p_HeAbundanceSurface) / (limWR - limOB)));
}


/*
 * CalculateMLRateRotationEnhancement_Langer1998
 *
 * @brief
 * Calculate mass loss rate enhancement for rapidly rotating stars, per Langer 1998, eq 3
 * (See https://ui.adsabs.harvard.edu/abs/1998A%26A...329..551L/abstract)
 * 
 * The exponent originally comes from Bjorkman & Cassinelli 1993
 * (see https://ui.adsabs.harvard.edu/abs/1993ApJ...409..429B/abstract), based on a fit to data
 * from Friend & Abbott 1986 (see https://ui.adsabs.harvard.edu/abs/1986ApJ...311..701F/abstract) 
 *
 * 
 * double CalculateMLRateRotationEnhancement_Langer1998() const
 *
 * @return                                      Mass loss enhancement factor for rapidly rotating stars
 */
GNU_CONST inline double CH::CalculateMLRateRotationEnhancement_Langer1998() const {
    return PPOW((1.0 - Omega() / OmegaBreak()), -0.43);
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
 * If option `--enable-rotationally-enhanced-mass-loss` was specified,
 * the mass loss rate will be enhance due to rotation.
 * 
 * 
 * MassLossT CalculateMLRate_Belczynski2010(
 *     const double p_Metallicity,
 *     const double p_Mass,
 *     const double p_Radius,
 *     const double p_Luminosity,
 *     const double p_Temperature,
 *     const double p_Perturb,
 *     const double p_HeAbundanceSurface,
 * ) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_Temperature                   Temperature of the star (Tsol)
 * @param       p_Perturb                       Small envelope perturbation parameter
 * @param       p_HeAbundanceSurface            Helium abundance at the surface of the star
 * @return                                      Tuple containing:
 *                                                   double  Mass loss rate (Msol yr^-1)
 *                                                   ML_TYPE Dominant mass loss type
 *                                                                  (will be ML_TYPE::WR or ML_TYPE::OB)
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
 * @param       p_LBVMLPrescription             LBV ML prescription
 * @return                                      Tuple containing:
 *                                                  double  Mass loss rate (Msol yr^-1)
 *                                                  ML_TYPE Dominant mass loss type (could be ML_TYPE::NONE)
 */
MassLossT CH::CalculateMLRate_Belczynski2010(
    const double                     p_Metallicity,
    const double                     p_Mass,
    const double                     p_Radius,
    const double                     p_Luminosity,
    const double                     p_Temperature,
    const double                     p_Perturb,
    const double                     p_HeAbundanceSurface,
    const double                     p_CoolWinMultiplier,
    const double                     p_LBVFactor,
    const LBV_ML_PRESCRIPTION p_LBVMLPrescription
) const {

    // Calculate OB mass loss rate
    double dMdt;
    ML_TYPE DominantMLtype;
    std::tie(dMdt, dominantMLtype) = BaseStar::CalculateMLRate_Belczynski2010(p_Mass, p_Radius, p_Luminosity, p_Temperature, p_Perturb, p_HeAbundanceSurface);

    // Scale mass loss with the surface helium abundance if required
    // (transition between OB and WR mass loss rates)
    if (OPTIONS->ScaleCHEMassLossWithSurfaceHeliumAbundance()) {

        double fractionOB = CalculateMLfractionOB(p_HeAbundanceSurface);

        if (fractionOB < 1.0) {
            double dMdtWR = 0.0;  
            ML_TYPE DominantMLtypeWR;
            std::tie(dMdtWR, dominantMLtypeWR) = CalculateMLRateWR_ZDependent(p_Metallicity, p_Luminosity, 0.0); // *Ilya* should we use p_Perturb here (since we have it)?
        
            dMdt   *= fractionOB;
            dMdtWR *= (1.0 - fractionOB);

            if (dMdtWR > dMdt) dominantMLtype = dominantMLtypeWR;   // Dominant mass loss type
        
            dMdt += dMdtWR;                                         // Scaled mass loss rate
        }
    }

    // Enhance mass loss rate due to rotation if required
    if (OPTIONS->EnableRotationallyEnhancedMassLoss()) dMdt *= CalculateMLRateRotationEnhancement_Langer1998();

    return std::make_tuple(dMdt, dominantMLtype);
}


/*
 * CalculateMLRate_Merritt2025
 *
 * Calculate mass loss rate, and dominant mass loss type, at the current evolutionary phase,
 * per Merritt et al., 2025.
 * 
 * If option `--scale-mass-loss-with-surface-helium-abundance` was specified,
 * the mass loss rate will be scaled with the surface helium abundance.
 *
 * If option `--enable-rotationally-enhanced-mass-loss` was specified, the
 * mass loss rate will be enhance due to rotation.
 * 
 * 
 * MassLossT CalculateMLRate_Merritt2025(double, const double p_Mass, double, const double p_Luminosity, const double p_Temperature, const double p_HeAbundanceSurface, double, double, double, double) const
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
MassLossT CH::CalculateMLRate_Merritt2025(
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
    
    // Calculate OB mass loss rate 
    double dMdt;
    ML_TYPE DominantMLtype;
    std::tie(dMdt, dominantMLtype) = CalculateMLRateOB(p_Mass, p_Luminosity, p_Temperature);

    // Scale mass loss with the surface helium abundance if required
    // (transition between OB and WR mass loss rates)
    if (OPTIONS->ScaleCHEMassLossWithSurfaceHeliumAbundance()) {

        // Calculate weight for combining OB and WR mass loss rates into total mass loss rate
        const double fractionOB = CalculateMLfractionOB(p_HeAbundanceSurface);

        if (fractionOB < 1.0) {

            // We need to know what the mass loss rate would be if the star was a HeMS star.
            // We can clone this object as an HG object and have the clone calculate the mass
            // loss rate for a HeMS star.
            //
            // The clone should not evolve, and so should not log anything, but to be sure the
            // clone does not participate in logging, we set its persistence to EPHEMERAL.
            std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::NAKED_HELIUM_STAR_MS, OBJECT_PERSISTENCE::EPHEMERAL);
            double dMdtWR = clone->CalculateMassLossRateMerritt2025(    // Calculate WR mass loss rate  
                p_Metallicity,
                p_Mass,
                p_Radius,
                p_Luminosity,
                p_Temperature,
                p_HeAbundanceSurface,
                p_Perturb,
                p_CoolWindMultiplier,
                p_LBVFactor,
                p_mStart
            );            

            dMdt   *= fractionOB;
            dMdtWR *= (1.0 - fractionOB);

            if (dMdtWR > dMdt) dominantMLtype = dominantMLtypeWR;       // Dominant mass loss type
        
            dMdt += dMdtWR;                                             // Scaled mass loss rate
        }
    }

    // Enhance mass loss rate due to rotation if required
    if (OPTIONS->EnableRotationallyEnhancedMassLoss()) dMdt *= CalculateMLRateRotationEnhancement_Langer1998();

    return std::make_tuple(dMdt, dominantMLtype);
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLogLuminositiesRatio_Szecsi2020
 *
 * @brief
 * Calculate the ratio of the log luminosities of a CH star and a non-CH MS of the same mass
 * (log(L_CH) / log(L_MS)).  Uses a polynomial fit derived using IZw18 and IZw18CHE BoOST 
 * models from Szecsi et al. 2020 (https://arxiv.org/abs/2004.08203).
 *
 * 
 * double CalculateLogLuminositiesRatio_Szecsi2020(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      Log luminosities ratio
 */
GNU_CONST inline double CH::CalculateLogLuminositiesRatio_Szecsi2020(const double p_Mass) const {
         
    const double x  = std::log10(p_Mass);
    const double x2 = x * x;

    return 1.8261540986808193 + (-1.1636822407341694  * x) + (0.5876329884434304  * x2) + (-0.10236336828026288 * x * x2);
}


/*
 * CalculateLuminosityAtPhaseEnd
 *
 * @brief
 * Calculate the luminosity of a CH star at the end of the (CH) MS.  The luminosity
 * will be enhanced if option `--enhance-CHE-lifetimes-luminosities` was specified.
 * 
 * 
 * double CalculateLuminosityAtPhaseEnd(const double p_Mass) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAMS CH luminosity (Lsol)
 */
GNU_PURE inline double CH::CalculateLuminosityAtPhaseEnd(const double p_Mass) const {

    // unenhanced CH luminosity at the end of the MS is just MS luminosity at the end of the MS
    double luminosity = MainSequence::CalculateLuminosityAtPhaseEnd_Hurley2000(p_Mass);
    
    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {                       // enhance luminosity of CH stars?
                                                                            // yes
        // apply enhancement, which at TAMS is just the ratio log(L_CH) / log(L_MS)
        // enhancement should not reduce luminosity, so ratio is clamped to a minimum of +1.0
        luminosity = PPOW(10.0, std::log10(luminosity) * std::max(CalculateLogLuminositiesRatio_Szecsi2020(p_Mass), 1.0));
    }

    return luminosity;
}


/*
 * CalculateLuminosity_Hurley2000
 *
 * @brief
 * Calculate the luminosity of a CH star on the (CH) MS, per Hurley et al. 2000, eq 12
 * 
 * The luminosity will be enhanced if option `--enhance-CHE-lifetimes-luminosities` was specified.
 * 
 * 
 * double CalculateLuminosity_Hurley2000(
 *     const double p_Mass,
 *     const double p_Tau,
 *     const double p_Age,
 *     const double p_LZAMS,
 *     const double p_tMS,
 *     const double p_tBGB
 * ) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_Age                           Time elapsed since ZAMS (Myr)
 * @param       p_LZAMS                         ZAMS luminosity of the star (Lsol)
 * @param       p_tMS                           MS lifetime, tMS (Myr)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (Myr)
 * @return                                      CH luminosity (Lsol)
 */
double CH::CalculateLuminosity_Hurley2000(
    const double p_Mass,
    const double p_Tau,
    const double p_Age,
    const double p_LZAMS,
    const double p_tMS,
    const double p_tBGB
) const {

    // unenhanced CH luminosity is just MS luminosity
    double luminosity = MainSequence::CalculateLuminosity_Hurley2000(p_Mass, p_Age, p_LZAMS, p_tMS, p_tBGB);

    if (OPTIONS->EnhanceCHELifetimesLuminosities()) {                   // enhance luminosity of CH stars?
                                                                        // yes
        // enhancement should not reduce luminosity, so ratio is clamped to a minimum of +1.0
        // enhancement amount grows from 1 to logLuminosityRatio over main-sequence
        const double enhancement = 1.0 + (std::max(CalculateLogLuminositiesRatio(p_Mass), 1.0) - 1.0) * p_Tau * p_Tau; 

        luminosity = PPOW(10.0, std::log10(luminosity) * enhancement);  // apply enhancement
    }

    return luminosity;
}
