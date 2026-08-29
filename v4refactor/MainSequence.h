#pragma once

#include "constants.h"
#include "typedefs.h"
#include "profiling.h"
#include "utils.h"

#include "BaseStar.h"

class BaseStar;



class MainSequence: virtual public BaseStar {

public:

    MainSequence(const BaseStar& p_BaseStar) : BaseStar(p_BaseStar) {}


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

    GNU_PURE double CalculateHAbundanceCore() const override { return CalculateHAbundanceCore(Metallicty(), Tau(), MZAMS(), HeAbundanceCore()); } // JR FIX THIS: DONE
    GNU_PURE double CalculateHAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_MZAMS, const double p_HeAbundanceCore) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateHAbundanceCoreAtPhaseEnd() const override { return CalculateHAbundanceCore(); } // Same as on phase JR FIX THIS: DONE

    GNU_PURE double CalculateHAbundanceSurface() const override { return HAbundanceSurface(); } // JR FIX THIS: DONE
    GNU_PURE double CalculateHAbundanceSurfaceAtPhaseEnd() const override { return CalculateHAbundanceSurface(); } // Same as on phase JR FIX THIS: DONE

    
    GNU_PURE double CalculateHeAbundanceCore() const override { return CalculateHeAbundanceCore(Metallicty(), Tau(), MZAMS(), HeAbundanceCore()); } // JR FIX THIS: DONE
    GNU_PURE double CalculateHeAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_MZAMS, const double p_HeAbundanceCore) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateHeAbundanceCoreAtPhaseEnd() const override { return CalculateHeAbundanceCore(); } // Same as on phase JR FIX THIS: DONE

    GNU_PURE double CalculateHeAbundanceSurface() const override { return HeAbundanceSurface(); } // JR FIX THIS: DONE   
    GNU_PURE double CalculateHeAbundanceSurfaceAtPhaseEnd() const override { return CalculateHeAbundanceSurface(); } // Same as on phase JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   AGE, LIFETIME, TAU, TIMESCALES, TIMESTEP   //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateAgeAfterMassLoss_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const double p_tMS) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateAgeAfterMassLoss_Hurley2000() const override { // JR FIX THIS: DONE 
        return CalculateAgeAfterMassLoss_Hurley2000(Metallicity(), Mass(), Age(), Timescales(TS::MS));
    }

    
    GNU_PURE double CalculatePhaseLifetime_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_tBGB) const; // JR FIX THIS: DONE
    GNU_PURE double CalculatePhaseLifetime_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculatePhaseLifetime_Hurley2000(Metallicity(), Mass(), Timescales(TS::BGB));
    }


    GNU_CONST double CalculateTauAtPhaseEnd_Hurley2000(const double p_Age, const DblVectortT& p_tScales) const override { return 1.0; } // 1.0 at phase end JR FIX THIS: DONE


    GNU_CONST double CalculateTau_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const override; // JR FIX THIS: DONE


    GNU_PURE TimescalesT CalculateTimescales_Hurley2000(const double p_Metallicity, const double p_Mass, const TimescalesT& p_tScales) const override; // JR FIX THIS: DONE
    GNU_PURE TimescalesT CalculateTimescales_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateTimescales_Hurley2000(Metallicity(), Mass0(), Timescales());
    }


    //////////////////////////////////////////////////
    //   ENVELOPE                                   //
    //////////////////////////////////////////////////

    GNU_CONST Dbl_DblT CalculateConvectiveEnvelopeMass(const double p_Mass, const double p_Tau) const override { // JR FIX THIS: DONE
    GNU_PURE Dbl_DblT CalculateConvectiveEnvelopeMass() const override { // JR FIX THIS: DONE
        return CalculateConvectiveEnvelopeMass(Mass(), Tau());
    }


    GNU_CONST Dbl_DblT CalculateConvectiveEnvelopeMass_Hurley2000(const double p_Mass, const double p_Tau) const; // JR FIX THIS: DONE
        return CalculateConvectiveEnvelopeMass_Hurley2000(p_Mass, p_Tau);
    }
    

    GNU_PURE inline double CalculateConvectiveEnvelopeRadialExtent() const override { // JR FIX THIS: DONE
        return CalculateConvectiveEnvelopeRadialExtent_Hurley2002(Metallicity(), Mass(), Radius(), Tau());
    }


    GNU_PURE double CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_Tau) const override; // JR FIX THIS: DONE


    //////////////////////////////////////////////////
    //   LUMINOSITY                                 //
    //////////////////////////////////////////////////
 
    GNU_PURE double CalculateLuminosity( // JR FIX THIS: DONE
        const double p_Metallicity,
        const double p_Mass,
        const double p_Age,
        const double p_MZAMS,
        const double p_CoreMass,
        const double p_tMS,
        const bool   p_TryBrcek = true
    ) const;
    GNU_PURE double CalculateLuminosity() const override { // JR FIX THIS: DONE
        return CalculateLuminosity(Metallicity(), Mass(), Age(), MZAMS(), CoreMass(), Timescales(TS::MS));
    }

    
    GNU_PURE double CalculateLuminosityAtPhaseEnd_Hurley2000(const double p_Metallicity, const double p_Mass) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateLuminosityAtPhaseEnd_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateLuminosityAtPhaseEnd_Hurley2000(Metallicity(), Mass0());
    }


    GNU_PURE OptDblT CalculateLuminosity_Brcek2025( // JR FIX THIS: DONE
        const double p_Metallicity,
        const double p_Mass,
        const double p_Age,
        const double p_MZAMS,
        const double p_CoreMass,
        const double p_HeAbundanceCore,
        const double p_tMS
    ) const;


    GNU_PURE double CalculateLuminosity_Hurley2000( // JR FIX THIS: DONE
        const double       p_Metallicity,
        const double       p_Mass,
        const double       p_Age,
        const double       p_LZAMS,
        const double       p_MHook,
        const TimescalesT& p_tScales
    ) const;
    GNU_PURE double CalculateLuminosity_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateLuminosity_Hurley2000(Metallicity(), Mass(), Age(), LZAMS(), ZDEP->HurleMCOHook(Metallicity()), Timescales());
    }


    GNU_PURE double CalculateLuminosity_NoBrcek() const override { // JR FIX THIS: DONE
        return CalculateLuminosity(Metallicity(), Mass(), Age(), MZAMS(), CoreMass(), Timescales(TS::MS), false);
    }


    GNU_PURE OptDblT CalculateLuminosity_Shikauchi2024(const double p_Metallicity, const double p_MZAMS, const double p_CoreMass, const double p_HeAbundanceCore) const; // JR FIX THIS: DONE

        
    //////////////////////////////////////////////////
    //   MASS                                       //
    //////////////////////////////////////////////////

    GNU_PURE  double CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024(const double p_MZAMS) const; // JR FIX THIS: DONE
    GNU_CONST double CalculateCNOprocessedCoreMass_Brcek2025(const double p_Mass, const double p_HeAbundanceCore) const; // JR FIX THIS: DONE


    double CalculateConvectiveCoreMass_Shikauchi2024(const double p_Tau) const; // JR FIX THIS: DONE
    double CalculateConvectiveCoreMass() const override { // JR FIX THIS: DONE
        return CalculateConvectiveCoreMass_Shikauchi2024(Tau());
    }


    GNU_CONST double CalculateCOCoreMass() const override { return 0.0; } // McCO = 0.0 for MS stars
    GNU_CONST double CalculateCOCoreMassAtPhaseEnd() const override { return 0.0; } // McCO = 0.0 for MS stars


    inline double CalculateCoreMass() const override { // JR FIX THIS: DONE
        return CalculateCoreMass(Metallicity(), Mass(), Luminosity(), MZAMS(), CoreMass(), HeAbundanceCore(), dt(), dMdt());
    }
    COMPAS_PURE double CalculateCoreMass( // JR FIX THIS: DONE
        const double p_Metallicity,
        const double p_Mass,
        const double p_Luminosity,
        const double p_MZAMS,
        const double p_CoreMass,
        const double p_HeAbundanceCore,
        const double p_dt, 
        const double p_dMdt
    ) const;

    GNU_CONST double CalculateCoreMassAtPhaseEnd() const { // JR FIX THIS: DONE
        // Accounts for minimal core mass built up prior to mass loss through mass transfer.
        // Core mass can't exceed total mass
        return (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::MANDEL) ? std::min(CoreMass(), Mass()) : 0.0;
    }

    COMPAS_PURE double CalculateCoreMass_Brcek2025( // JR FIX THIS: DONE - EXCEPT FOR CHECK IN FUNCTION BODY
        const double p_Metallicity,
        const double p_Mass,
        const double p_Luminosity,
        const double p_MZAMS,
        const double p_CoreMass,
        const double p_HeAbundanceCore,
        const double p_dt,
        const double p_dMdt,
    ) const;


    GNU_PURE double CalculateEffectiveInitialMass_Hurley2000() const override { return Mass(); } // Per Hurley et al. 2000, section 7.1


    GNU_CONST double CalculateHeCoreMass() const { return 0.0; } // McHe = 0.0 for MS stars
    GNU_CONST double CalculateHeCoreMassAtPhaseEnd() const { return 0.0; } // McHe = 0.0 for MS stars


    //////////////////////////////////////////////////
    //   RADIUS                                     //
    //////////////////////////////////////////////////

    GNU_PURE double CalculateConvectiveCoreRadius(const double p_Mass, const double p_Tau) const; // JR FIX THIS: DONE


    GNU_PURE double CalculateRadius(const double p_Mass, const double p_Radius, const double p_Age, const double p_tMS, const bool p_TryBrcek = true) const; // JR FIX THIS: DONE
    GNU_PURE double CalculateRadius() const override { // JR FIX THIS: DONE 
        return CalculateRadius(Mass(), Radius(), Age(), Timescales(TS::MS));
    }


    GNU_PURE double CalculateRadiusAtPhaseEnd_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_RZAMS) const;
    GNU_PURE double CalculateRadiusAtPhaseEnd_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateRadiusAtPhaseEnd_Hurley2000_Static(Metallicity(), Mass(), RZAMS());
    }


    GNU_PURE double CalculateRadiusOnMassChange(double p_dM) const override {
        return CalculateRadius(Mass() + p_dM, astro::CalculateZAMSRadius_Tout1996(Metallicity(), Mass() + p_dM), Age(), Timescales(TS::MS));
    }


    COMPAS_PURE double CalculateRadius_Brcek2025(const double p_Mass, const double p_Radius, const double p_Age, const double p_tMS) const; // JR FIX THIS: DONE


    GNU_PURE double CalculateRadius_Hurley2000(
        const double p_Metallicity,
        const double p_Mass,
        const double p_Tau,
        const double p_RZAMS,
        const double p_MHook,
        const double p_tBGB
    ) const {
    GNU_PURE double CalculateRadius_Hurley2000() const override { // JR FIX THIS: DONE
        return CalculateRadius_Hurley2000(Metallicity(), Mass(), Tau(), RZAMS(), ZDEP->HurleyMCOHook(Metallicity()), Timescales(TS::BGB));
    }
    GNU_PURE double CalculateRadius_Hurley2000(const double p_Mass) const override { // JR FIX THIS: DONE
        return CalculateRadius_Hurley2000(Metallicity(), p_Mass, Tau(), RZAMS(), ZDEP->HurleyMCOHook(Metallicity()), Timescales(TS::BGB));
    }


    GNU_PURE double CalculateRadius_NoBrcek() const override { // JR FIX THIS: DONE 
        return CalculateRadius(Mass(), Radius(), Age(), Timescales(TS::MS), false);
    }
    
     























    







    










 



///// PHASE END, ETC.      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<







    



    



















    double          CalculateMomentOfInertia() const                                        { return CalculateMomentOfInertia(Mass(), Radius()); }

    // JR FIX THIS: is this required??? see baseStar
    double          CalculateMomentOfInertia(const double p_Mass, const double p_Radius) const { return 0.1 * p_Mass * p_Radius * p_Radius; } // k2 = 0.1 as defined in Hurley et al. 2000, after eq 109.


GNU_CONST inline double CalculateHurleyPerturbationMu() const { return 5.0; } // Hurley et al. 2000, eqs 97 & 98








    GNU_CONST double ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const override;


    STELLAR_TYPE    EvolveToNextPhase()                                                     { return STELLAR_TYPE::HERTZSPRUNG_GAP; }





                                           



   
    



                                                                                            }

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
 * Calculate the hydrogen abundance in the core of the star, given the phase-relative age
 * of the star.  Currently just a simple linear model that assumes that hydrogen in the
 * core of the star is burned to helium at a constant rate throughout the star's lifetime.
 * 
 * Should be updated to match detailed models.
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
 * double CalculateHAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_MZAMS, const double p_HeAbundanceCore) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core [0, 1]
 * @return                                      Hydrogen abundance in the core [0, 1]
 */
inline double MainSequence::CalculateHAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_MZAMS, const double p_HeAbundanceCore) const {   
    return ((OPTIONS->MSCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK) && (BRCEK_LOWER_MASS_LIMIT < p_MZAMS))
            ? 1.0 - p_HeAbundanceCore - p_Metallicity
            : ZDEP->ZAMSHeAbundance(p_Metallicity) * (1.0 - p_Tau);
}


/*
 * CalculateHeAbundanceCore
 *
 * @brief
 * Calculate the helium abundance in the core of the star, given the phase-relative age
 * of the star.  Currently just a simple linear model that assumes that hydrogen in the
 * core of the star is burned to helium at a constant rate throughout the star's lifetime.
 * 
 * Should be updated to match detailed models.
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
 * double CalculateHeAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_MZAMS, const double p_HeAbundanceCore) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core [0, 1]
 * @return                                      Helium abundance in the core [0, 1]
 */
inline double MainSequence::CalculateHeAbundanceCore(const double p_Metallicity, const double p_Tau, const double p_MZAMS, const double p_HeAbundanceCore) const {
    return ((OPTIONS->MSCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK) && (BRCEK_LOWER_MASS_LIMIT < p_MZAMS))
            ? p_HeAbundanceCore
            : ((1.0 - p_Metallicity - ZDEP->ZAMSHeAbundance(p_Metallicity)) * p_Tau) + ZDEP->ZAMSHeAbundance(p_Metallicity);
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
 * Calculate the star's age after mass loss, per Hurley et al. 2000, section 7.1
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
inline double MainSequence::CalculateAgeAfterMassLoss_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_Age, const double p_tMS) const {
    // We call astro::CalculateLifetimeToBGB_Hurley2000() here rather than use timescales
    // because p_Mass may not be the same mass used to calculate timescales[TS::BGB]
    const double tBGBprime = astro::CalculateLifetimeToBGB_Hurley2000(p_Metallicity, p_Mass);
    return p_Age * CalculatePhaseLifetime_Hurley2000(p_Metallicity, p_Mass, tBGBprime) / p_tMS;
}


/*
 * CalculatePhaseLifetime_Hurley2000
 *
 * @brief
 * Calculate lifetime of the Main Sequence, tMS, per Hurley et al. 2000, eqs 5, 6, & 7.
 *
 *
 * double CalculatePhaseLifetime_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_tBGB) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_tBGB                          Time to Base of Giant Branch, tBGB (Myr)
 * @return                                      MS lifetime, tMS (Myr)
 */
inline double MainSequence::CalculatePhaseLifetime_Hurley2000(const double p_Metallicity, const double p_Mass, const double p_tBGB) const {
    const auto& a = ZDEP->HurleyACoeffs(p_Metallicity);   // Hurley a coefficients
    const double tHook = std::max(0.5, (1.0 - (0.01 * std::max((a[6] / PPOW(p_Mass, a[7])), (a[8] + (a[9] / PPOW(p_Mass, a[10]))))))) * p_tBGB;
    return std::max(tHook, std::clamp(0.95 - (0.03 * (ZDEP->ZetaHurley(p_Metallicity), + 0.30103)), 0.95, 0.99) * p_tBGB);
}


/*
 * CalculateTau_Hurley2000
 *
 * @brief
 * Calculate the MS-relative age (fractional Main Sequence age) of the star,
 * per Hurley et al. 2000, eq 11
 *
 * 
 * double CalculateTau_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      MS-relative age of the star, [0, 1]
 */
inline double MainSequence::CalculateTau_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const {
    return std::clamp(p_Age / p_tScales[TS::MS], 0.0, 1.0);
}


/*
 * CalculateTimescales_Hurley2000
 *
 * @brief
 * Calculate timescales, per Hurley at al. 2000.
 * 
 * Since timescales depend on a star's mass, they need to be calculated whenever the mass of the
 * star changes (probably every timestep).
 *
 * p_tScales is copied once and returned by value.  NRVO constructs the returned array directly
 * in the caller's return slot, so no copy on the return.
 *
 * 
 * TimescalesT CalculateTimescales_Hurley2000(const double p_Metallicity, const double p_Mass, const TimescalesT& p_tScales) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Mutated timescales (Myr)
 */
inline TimescalesT MainSequence::CalculateTimescales_Hurley2000(const double p_Metallicity, const double p_Mass, const TimescalesT& p_tScales) const {

    TimescalesT tScales = std::move(p_tScales);

    tScales[TS::BGB] = astro::CalculateLifetimeToBGB_Hurley2000(p_Metallicity, p_Mass);
    tScales[TS::MS]  = CalculatePhaseLifetime_Hurley2000(p_Metallicity, p_Mass, tScales[TS::BGB]);

    return tScales;
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
 * double ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales)
 *
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tScales                       Timescales (Myr)
 * @return                                      Suggested timestep (Myr)
 */
inline double MainSequence::ChooseTimestep_Hurley2000(const double p_Age, const TimescalesT& p_tScales) const {

          double dtk = 1.0E-2 * p_tScales[TS::MS];                  // Stellar type specific dt (Hurley sse uses 0.05)
    const double dte = p_tScales[TS::MS] - p_Age;                   // Time to end of phase

    // Check that dtk is short enough to resolve the hook at the end of the MS for HM stars.  // *Ilya*  why not check for HM star here?
    // If not, go an order-of-magnitude shorter
    if (dte < dtk) dtk /= 10.0; // *Ilya* why stop at one order of magnitude smaller?  Why not repeat the check until dtk <= dte?  We just clamp to dte anyway...

    return std::max(std::min(dtk, dte), NUCLEAR_MINIMUM_TIMESTEP);  // Clamp to minimum NUCLEAR_MINIMUM_TIMESTEP
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     ENVELOPE                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateConvectiveEnvelopeMass_Hurley2000
 *
 * @brief
 * Calculate the convective envelope mass for Main Sequence, MS, stars.
 *
 * Based on Hurley et al. 2000, section 7.2, after eq 111
 *
 *
 * Dbl_DblT CalculateConvectiveEnvelopeMass_Hurley2000(const double p_Mass, const double p_Tau) const
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      Tuple containing:
 *                                                   MS convective envelope mass (Msol)
 *                                                   Maximum MS (ZAMS) convective envelope mass (Msol)
 */
Dbl_DblT MainSequence::CalculateConvectiveEnvelopeMass_Hurley2000(const double p_Mass, const double p_Tau) const {

    double mEnv  = 0.0; // Default return value
    double mEnv0 = 0.0; // Default return value

    if (p_Mass <= 1.25) {   // No convective envelope above 1.25 Msol
        mEnv0 = (p_Mass > 0.35) ? 0.35 * (1.25 - p_Mass) * (1.25 - p_Mass) / 0.81 : p_Mass;
        mEnv  = mEnv0 * std::sqrt(std::sqrt(1.0 - p_Tau));
    }

    return std::make_tuple(mEnv, mEnv0);
}


/*
 * CalculateConvectiveEnvelopeRadialExtent_Hurley2002
 *
 * @brief
 * Calculate the radial extent of the star's convective envelope (if it has one),
 * per Hurley et al. 2002, sec. 2.3, particularly subsec. 2.3.1, eqs 36-38.
 *
 * Uses radius of a 0.35 solar mass star at ZAMS rather than at fractional age Tau,
 * but such low-mass stars only grow by a maximum factor of 1.5
 * [just above Eq. (10) in Hurley, Pols, Tout (2000)], so this is a reasonable approximation.
 * 
 *
 * double CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_Tau) const
 * 
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      Radial extent of the star's convective envelope (Rsol)
 */
inline double MainSequence::CalculateConvectiveEnvelopeRadialExtent_Hurley2002(const double p_Metallicity, const double p_Mass, const double p_Radius, const double p_Tau) const {

    double radiusEnvelope = p_Radius; // Default return value

         if (1.25 < p_Mass) radiusEnvelope = 0.0;
    else if (p_Mass > 0.35) radiusEnvelope = CalculateZAMSRadius_Tout1996(p_Metallicity, 0.35) * std::sqrt((1.25 - p_Mass) / 0.9);

    return radiusEnvelope * std::sqrt(std::sqrt(1.0 - p_Tau));
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    LUMINOSITY                                     //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * CalculateLuminosity
 *
 * @brief
 * Calculate luminosity on the main sequence.
 *
 * If the Brcek MS core mass prescription was specified by the user, and the star is in the 
 * Brcek regime, calculate the MS luminosity per Brcek et al. 2025, otherwise use the method
 * appropriate for the evolution mode.
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
 * double CalculateLuminosity(
 *     const double p_Metallicity,
 *     const double p_Mass,
 *     const double p_Age,
 *     const double p_MZAMS,
 *     const double p_CoreMass,
 *     const double p_tMS,
 *     const bool   p_TryBrcek
 * ) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_tMS                           MS lifetime, tMS (Myr)
 * @param       p_TryBrcek                      Flag to indicate if Brcek luminosity calculation should be tried first (default is true)
 * @return                                      Luminosity of the star (Lsol)
 */
inline double MainSequence::CalculateLuminosity(
    const double p_Metallicity,
    const double p_Mass,
    const double p_Age,
    const double p_MZAMS,
    const double p_CoreMass,
    const double p_tMS,
    const bool   p_TryBrcek
) const {

    OptDblT luminosity;                                             // No default here - set later

    // Try Brcek - std::nullopt if not in Brcek regime.
    // (Returns immediately if user did not specify Brcek core mass prescription)
    if (p_TryBrcek) luminosity = CalculateLuminosity_Brcek2025(p_Metallicity, p_Mass, p_Age, p_MZAMS, p_CoreMass, ZDEP->HeAbundanceCore(), p_tMS);

    if (!luminosity.has_value()) {                                  // Have luminosity from Brceek?
        luminosity = 0.0;                                           // No - set default return value, and use method appropriate for evolution mode

        switch (OPTIONS->Mode()) {                                  // Which evolution mode?

            case EVOLUTION_MODE::SSE_HURLEY:                        // HURLEY SSE
            case EVOLUTION_MODE::BSE_HURLEY:                        // HURLEY BSE
                luminosity = CalculateLuminosity_Hurley2000();
                break;

            default:                                                // Unexpected mode
                // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
                // accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user chose a mode this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing mode or, if the missing mode is
                // superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);      // Throw error
        }
    }

    return luminosity;
}


/*
 * CalculateLuminosityAtPhaseEnd_Hurley2000
 *
 * @brief
 * Calculate luminosity at end of Main Sequence, per Hurley et al. 2000, eq 8
 *
 *
 * double CalculateLuminosityAtPhaseEnd_Hurley2000(const double p_Metallicity, const double p_Mass) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @return                                      TAMS luminosity (Lsol)
 */
inline double MainSequence::CalculateLuminosityAtPhaseEnd_Hurley2000(const double p_Metallicity, const double p_Mass) const {
   
    auto& a = ZDEP->HurleyACoeffs(p_Metallicity);   // Hurley a coefficients

    const double m2 = p_Mass * p_Mass;
    const double m3 = p_Mass * m2;

    return (a[11] * m3) + (a[12] * m2 * m2) + (a[13] * PPOW(p_Mass, (a[16] + 1.8))) / (a[14] + (a15 * m2 * m3) + PPOW(p_Mass, a[16]));
}


/*
 * CalculateLuminosityAtZAMS
 *
 * @brief
 * Calculate the luminosity of a star at ZAMS.
 * 
 * Calls relevant luminosity function based on the evolutionary mode given in program options.
 *
 *
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * !*!*!*!*! ZAMS attribute warning *!*!*!*!*!
 * !*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!
 * 
 * This function ostensibly relies on the values of the ZAMS mass of the star,
 * and should not be used if the ZAMS mass is not known.
 * 
 * 
 * double CalculateLuminosityAtZAMS(const double p_Metallicity, const double p_MZAMS) const
 *
 * @param       p_MZAMS                         Metallicity of the star
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS luminosity (Lsol)
 */
inline double MainSequence::CalculateLuminosityAtZAMS(const double p_Metallicity,const double p_MZAMS) const {

    double luminosity;

    switch (OPTIONS->Mode()) {                                              // which evolution mode?

        case EVOLUTION_MODE::SSE_HURLEY:                                         // HURLEY SSE
        case EVOLUTION_MODE::BSE_HURLEY:                                         // HURLEY BSE
            luminosity = CalculateLuminosityAtZAMS_Tout1996(p_Metallicity, p_MZAMS);
            break;
        
        default:                                                            // unknown mode
            // the only way this can happen is if someone added an EVOLUTION_MODE and it isn't
            // accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the user chose a mode this code doesn't account for, and that should
            // be flagged as an error and result in termination of the evolution of the star or binary.
            // The correct fix for this is to add code for the missing mode or, if the missing mode is
            // superfluous, remove it from the option.

            THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                     // throw error
    }       

    return luminosity;
}


/*
 * CalculateLuminosity_Brcek2025
 *
 * @brief
 * Calculate luminosity on the main sequence.
 *
 * This function returns a std::optional<double> value.  If the Brcek MS core mass prescription
 * was specified by the user, and the star is in the Brcek regime, the function will calculate
 * the MS luminosity per Brcek et al. 2025 and return the luminosity (a .has_value() check of
 * the returned value will return true), otherwise std::nullopt is returned (a .has_value()
 * check of the returned value will return false).
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
 * OptDblT CalculateLuminosity_Brcek2025(
 *     const double p_Metallicity,
 *     const double p_Mass,
 *     const double p_Age,
 *     const double p_MZAMS,
 *     const double p_CoreMass,
 *     const double p_HeAbundanceCore,
 *     const double p_tMS,
 * ) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core [0, 1]
 * @param       p_tMS                           MS lifetime, tMS (Myr)
 * @return                                      Luminosity of the star (Lsol) or std::nullopt
 */
inline OptDblT MainSequence::CalculateLuminosity_Brcek2025(
    const double p_Metallicity,
    const double p_Mass,
    const double p_Age,
    const double p_MZAMS,
    const double p_CoreMass,
    const double p_HeAbundanceCore,
    const double p_tMS
) const {

    OptDbl luminosity;                                                              // Default return value is std::nullopt

    if (OPTIONS->MSCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK) {    // User specified BRCEK core mass prescription?
        if (p_MZAMS >= BRCEK_LOWER_MASS_LIMIT) {                                    // Yes - in BRCEK regime?
                                                                                    // Yes
            const double ageAtHook = 0.99 * p_tMS;                                  // Age at hook start
            if (p_Age > ageAtHook) {                                                // In MS hook?
                                                                                    // Yes
                // Calculate the luminosity at hook start.
                // In the hook, core helium abundance is fixed at 1-Z and core mass is not changing.
                // Try Shikauchi first - if not in Shikauchi regime, use the method appropriate for
                // the evolution mode.
                OptDbl lHook = CalculateLuminosity_Shikauchi2024(p_Metallicity, p_MZAMS, p_CoreMass, p_HeAbundanceCore);
                if (!lHook.has_value()) lHook = CalculateLuminosity_NoBrcek(p_Metallicity, p_Mass, ageAtHook, p_MZAMS, p_CoreMass, p_tMS);

                // We need the TAMS luminosity, which is just the luminosity at the start of the HG phase.
                // Since we are on the main sequence here, we can clone this object as an HG object and, as
                // long as it is initialised (to correctly set Tau to 0.0 on the HG phase), we can query the
                // cloned object for its luminosity.
                //
                // The clone should not evolve, and so should not log anything, but to be sure the
                // clone does not participate in logging, we set its persistence to EPHEMERAL.
                std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::HERTZSPRUNG_GAP, OBJECT_PERSISTENCE::EPHEMERAL);

                // Linear interpolation
                luminosity = (lHook.value() * (p_tMS - p_Age) + clone->Luminosity() * (p_Age - ageAtHook)) / (p_tMS - ageAtHook);
            }
            else {                                                                  // No, not in MS hook
                // Try Shikauchi - std::nullopt if not in Shikauchi regime
                luminosity = CalculateLuminosity_Shikauchi2024(p_Metallicity, p_MZAMS, p_CoreMass, p_HeAbundanceCore);
            }
        }
    }

    return luminosity;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       MASS                                        //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024
 *
 * @brief
 * Calculate the CNO-processed core mass of a main sequence star at ZAMS,
 * per Shikauchi et al. 2024
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
 * double CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024(const double p_Metallicity, const double p_MZAMS) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @return                                      ZAMS CNO-processed core mass (Msol)
 */
inline double MainSequence::CalculateCNOprocessedCoreMassAtZAMS_Shikauchi2024(const double p_Metallicity, const double p_MZAMS) const {
    DblVectorT fCoeffs = ZDEP->ShikauchiFMixCoefficients(p_Metallicity); // Shikauchi fMix coefficients
    return p_MZAMS * (fCoeffs[0] + fCoeffs[1] * std::exp(-p_MZAMS / fCoeffs[2]));
}


/*
 * CalculateCNOprocessedCoreMass_Brcek2025
 *
 * @brief
 * Calculate the CNO-processed core mass of a main sequence star after full mixing (due to
 * merger or CHE) for an arbitrary central helium fraction, per Brcek et al. 2025
 * 
 *
 * double CalculateCNOprocessedCoreMass_Brcek2025(const double p_Mass, const double p_HeAbundanceCore)
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core [0, 1]
 * @return                                      CNO-processed core mass (Msol)
 */
inline double MainSequence::CalculateCNOprocessedCoreMass_Brcek2025(const double p_Mass, const double p_HeAbundanceCore) const {
    DblVectorT fCoeffs = BRCEK_FMIX_COEFFICIENTS; // Brcek fMix coefficients (from constants.h)
    const double mh    = p_Mass * PPOW(10.0, p_HeAbundanceCore * (p_HeAbundanceCore + 2.0) / 4.0);
    return p_Mass * (fCoeffs[0] + fCoeffs[1] * std::exp(-mh / fCoeffs[2])) * PPOW(1.0 - fCoeffs[4] / mh, fCoeffs[3]);
}


/*
 * CalculateConvectiveCoreMass_Shikauchi2024
 *
 * @brief
 * Calculate the mass of the convective core on the main sequence.
 *
 * Based on Shikauchi, Hirai, Mandel (2024), core mass shrinks to 60% of
 * initial value over the course of the MS.
 *
 *
 * double CalculateConvectiveCoreMass_Shikauchi2024(const double p_Tau) const
 * 
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      Mass of convective core (Msol)
 */
inline double MainSequence::CalculateConvectiveCoreMass_Shikauchi2024(const double p_Tau) const {
    double finalConvectiveCoreMass   = DetermineTAMSCoreMass();  // Core mass at TAMS
    double initialConvectiveCoreMass = finalConvectiveCoreMass / 0.6;
    return initialConvectiveCoreMass - p_Tau * (initialConvectiveCoreMass - finalConvectiveCoreMass);
}


/*
 * CalculateCoreMass
 *
 * @brief
 * Calculate the convective core mass of a main sequence star, based on the
 * MS core mass prescription specified by the user.
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
 * double CalculateCoreMass(
 *     const double p_Metallicity,
 *     const double p_Mass,
 *     const double p_Luminosity,
 *     const double p_MZAMS,
 *     const double p_CoreMass,
 *     const double p_HeAbundanceCore,
 *     const double p_dt, 
 *     const double p_dMdt
 * ) const
 *
 * @param       p_Metallicity                   Metallicity of the star
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Luminosity                    Luminosity of the star (Lsol)
 * @param       p_MZAMS                         ZAMS mass of the star (Msol)
 * @param       p_CoreMass                      Core mass of the star (Msol)
 * @param       p_HeAbundanceCore               Helium abundance in the core [0, 1] 
 * @param       p_dt                            Time step (Myr)
 * @param       p_dMdt                          Mass loss rate (Msol yr^-1) 
 * @return                                      Core mass of the star (Msol)
 */
inline double MainSequence::CalculateCoreMass(
    const double p_Metallicity,
    const double p_Mass,
    const double p_Luminosity,
    const double p_MZAMS,
    const double p_CoreMass,
    const double p_HeAbundanceCore,
    const double p_dt, 
    const double p_dMdt
) const {

    double coreMass = p_CoreMass;                                       // Default return value

    switch (OPTIONS->MainSequenceCoreMassPrescription()) {              // Which MS core mass prescription?
        
        case MS_CORE_MASS_PRESCRIPTION::BRCEK:                          // BRCEK
            // Calculate MS core mass per Brcek et al. 2025, following Shikauchi et al. 2024.
            coreMass = CalculateCoreMass_Brcek2025(p_Metallicity, p_Mass, p_Luminosity, p_MZAMS, p_CoreMass, p_HeAbundanceCore, p_dt, p_dMdt);
            break;

        case MS_CORE_MASS_PRESCRIPTION::HURLEY:                         // HURLEY
            coreMass = 0.0;                                             // No MS core in Hurley et al. 2000
            break;
        
        case MS_CORE_MASS_PRESCRIPTION::MANDEL:                         // MANDEL
            coreMass = p_CoreMass;                                      // No change to core mass for SSE
            break;
        
        default:                                                        // Unexpected prescription
            // The only way this can happen is if the MS_CORE_MASS_PRESCRIPTION passed to this function
            // is not accounted for in this code.  We should not default here, with or without a warning.
            // We are here because the code passed a prescription that this function doesn't account
            // for, and that should be flagged as an error and result in termination of the evolution
            // of the star or binary.
            // The correct fix for this is to add code to this function for the missing prescription,
            // or fix the calling code to pass a prescription that is handled by this function.

            THROW_ERROR(ERROR::UNEXPECTED_MS_CORE_MASS_PRESCRIPTION);   // Throw error
    }

    return coreMass;
}


/*
 * DetermineTAMSCoreMass
 *
 * @brief
 * Determine the expected core mass at terminal age main sequence (TAMS),
 * i.e., at the start of the HG phase.
 *
 * double DetermineTAMSCoreMass() const
 *
 *
 * @return                                      TAMS core Mass (Msol)
 *
 */
inline double MainSequence::DetermineTAMSCoreMass() const {
    // Since we are on the main sequence here, we can clone this object as an HG object
    // and, as long as it is initialised (to correctly set Tau to 0.0 on the HG phase),
    // we can query the cloned object for its core mass.
    //
    // The clone should not evolve, and so should not log anything, but to be sure the
    // clone does not participate in logging, we set its persistence to EPHEMERAL.
    std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::HERTZSPRUNG_GAP, OBJECT_PERSISTENCE::EPHEMERAL);
    return clone->CoreMass(); // Return core mass of clone
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      RADIUS                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CalculateConvectiveCoreRadius
 *
 * @brief
 * Calculate the radial extent of the star's convective core (if it has one)
 *
 * Uses preliminary fit from Minori Shikauchi @ ZAMS, then a smooth interpolation to the HG
 *
 *
 * double CalculateConvectiveCoreRadius()
 * 
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Tau                           Phase-relative age of the star [0, 1]
 * @return                                      Radial extent of the star's convective core (Rsol)
 */
inline double MainSequence::CalculateConvectiveCoreRadius(const double p_Mass, const double p_Tau) const {

    double radius = 0.0;        // Default return value (no convective core)

    if (1.25 < p_Mass) {        // Star has convective core?
                                // Yes
        const double ZAMSconvectiveCoreRadius = p_Mass * (0.06 + 0.05 * std::exp(-p_Mass / 61.57));
    
        // We need the TAMS core radius, which is just the core radius at the start of the
        // HG phase. Since we are on the main sequence here, we can clone this object as an
        // HG object and, as long as it is initialised (to correctly set Tau to 0.0 on the
        // HG phase), we can  query the cloned object for its core radius.
        //
        // The clone should not evolve, and so should not log anything, but to be sure the
        // clone does not participate in logging, we set its persistence to EPHEMERAL.

        std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::HERTZSPRUNG_GAP, OBJECT_PERSISTENCE::EPHEMERAL);

        radius = ZAMSconvectiveCoreRadius - p_Tau * (ZAMSconvectiveCoreRadius - clone->CalculateRemnantRadius());
    }

    return radius;
}


/*
 * CalculateRadius
 *
 * @brief
 * Calculate radius on the main sequence.
 *
 * If the Brcek MS core mass prescription was specified by the user, and the star is in the 
 * Brcek regime, calculate the MS radius per Brcek et al. 2025, otherwise use the method
 * appropriate for the evolution mode.
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
 * double CalculateRadius(const double p_Mass, const double p_Radius, const double p_Age, const double p_tMS, const bool p_TryBrcek) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tMS                           MS lifetime, tMS (Myr)
 * @param       p_TryBrcek                      Flag to indicate if Brcek radius calculation should be tried first (default is true)
 * @return                                      Radius of the star (Rsol)
 */
inline double MainSequence::CalculateRadius(const double p_Mass, const double p_Radius, const double p_Age, const double p_tMS, const bool p_TryBrcek) const {

    OptDblT radius;                                                 // No default here - set later

    // Try Brcek if required - std::nullopt if not in Brcek regime.
    // (Returns immediately if user did not specify Brcek core mass prescription)
    if (p_TryBrcek) radius = CalculateRadius_Brcek2025(p_Mass, p_Radius, p_Age, p_tMS);

    if (!radius.has_value()) {                                  // Have radius from Brceek?
        radius = p_Radius;                                      // No - set default return value, and use method appropriate for evolution mode

        switch (OPTIONS->Mode()) {                              // Which evolution mode?

            case EVOLUTION_MODE::SSE_HURLEY:                    // HURLEY SSE
            case EVOLUTION_MODE::BSE_HURLEY:                    // HURLEY BSE
                radius = CalculateRadius_Hurley2000(p_Mass);
                break;

            default:                                            // Unexpected mode
                // The only way this can happen is if someone added an EVOLUTION_MODE and it isn't
                // accounted for in this code.  We should not default here, with or without a warning.
                // We are here because the user chose a mode this code doesn't account for, and that should
                // be flagged as an error and result in termination of the evolution of the star or binary.
                // The correct fix for this is to add code for the missing mode or, if the missing mode is
                // superfluous, remove it from the option.

                THROW_ERROR(ERROR::UNEXPECTED_EVOLUTION_MODE);  // Throw error
        }
    }

    return radius.value();
}


/*
 * CalculateRadius_Brcek2025
 *
 * @brief
 * Calculate radius on the main sequence, per Brcek et al. 2025
 *
 * This function returns a std::optional<double> value.  If the Brcek MS core mass prescription
 * was specified by the user, and the star is in the Brcek regime, the function will calculate
 * the MS radius per Brcek et al. 2025 and return the radius (a .has_value() check of the
 * returned value will return true), otherwise std::nullopt is returned (a .has_value() check
 * of the returned value will return false).
 * 
 *
 * OptDblT CalculateRadius_Brcek2025(const double p_Mass, const double p_Radius, const double p_Age, const double p_tMS) const
 *
 * @param       p_Mass                          Mass of the star (Msol)
 * @param       p_Radius                        Radius of the star (Rsol)
 * @param       p_Age                           Age of the star (Myr)
 * @param       p_tMS                           MS lifetime, tMS (Myr)
 * @return                                      Radius of the star (Rsol) or std::nullopt
 */
inline OptDblT MainSequence::CalculateRadius_Brcek2025(const double p_Mass, const double p_Radius, const double p_Age, const double p_tMS) const {

    OptDbl radius;                                                                  // No default return value - optional

    if (OPTIONS->MSCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK) {    // User specified BRCEK core mass prescription?
        if (p_MZAMS >= BRCEK_LOWER_MASS_LIMIT) {                                    // Yes - in BRCEK regime?
                                                                                    // Yes
            const double ageAtHook = 0.99 * p_tMS;                                  // Age at hook start
            if (p_Age > ageAtHook) {                                                // In MS hook?
                                                                                    // Yes
                // Calculate the radius at hook start.
                // Use the method appropriate for the evolution mode.
                rHook = CalculateRadius_NoBrcek(p_Mass, p_Radius, ageAtHook, p_tMS);

                // We need the TAMS radius, which is just the radius at the start of the HG phase.
                // Since we are on the main sequence here, we can clone this object as an HG object and, as
                // long as it is initialised (to correctly set Tau to 0.0 on the HG phase), we can query the
                // cloned object for its radius.
                //
                // The clone should not evolve, and so should not log anything, but to be sure the
                // clone does not participate in logging, we set its persistence to EPHEMERAL.
                std::unique_ptr<BaseStar> clone = CloneAs(STELLAR_TYPE::HERTZSPRUNG_GAP, OBJECT_PERSISTENCE::EPHEMERAL);

                // Select radius at TAMS from the HG clone or current radius (whichever is smaller).
                // This is for stars that were significantly stripped as this prevents radius expansion
                // during the hook, and delays possible mass transfer to the start of HG
                const double rTAMS = std::min(clone->Radius(), p_Radius);

                // Linear interpolation
                radius = (rHook * (p_tMS - p_Age) + rTAMS * (p_Age - ageAtHook)) / (p_tMS - ageAtHook);
            }
        }
    }

    return radius;
}









///////// constituent functions <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//



