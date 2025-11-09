#ifndef __Globals_H__
#define __Globals_H__

#include "constants.h"


#define GLOBALS Globals::Instance()


/*
 * Globals Singleton - "global" variables
 * 
 * This singleton class holds variables that we really only need to calculate or
 * initialise once, or that might be applicable for multiple stars/binaries so 
 * need to be calculated infrequently (i.e. not for every star/binary).
 *
 * Singletons and global variables are sometimes frowned-upon, but doing it this
 * way means the objects don't need to be passed around to all and sundry.
 * I think convenience and clarity sometimes trump dogma.
 */

class Globals {

private:

    ////////////////////////////////////////////////////////////////////////////////////////
    //                                                                                    //
    // NOTE: this code assumes that metallicity is same for both stars in a binary.       //
    //                                                                                    //
    // The code here could easily be changed to have metallicity different for each star, //
    // but COMPAS options currently don't allow that (and may never), so for now we use   //
    // that as a  simplification.                                                         //
    //                                                                                    //
    ////////////////////////////////////////////////////////////////////////////////////////



    // Values of the variables in this struct depend on the star's metallicity only - so these values only need to
    // be calculated once per star (upon creation), but can also be reused if metallicity doesn't change from one
    // star to the next (e.g. in a population run).
    
    // This struct contains variables to be (re)calculated for each star
    
    // This struct contains variables to be calculated per Hurley et al. 2000 & 2002
    // Values of the variables in this struct depend on the star's metallicity only - so these values only need to
    // be calculated once per star (upon creation), but can also be reused if metallicity doesn't change from one
    // star to the next (e.g. in a population run).
    typdef struct HurleyZdependent {

        std::optional<double>      rho;                     // Hurley et al. 2000 rho, p24
        std::optional<double>      radiusXexponent;         // Hurley et al. 2000 'x' exponent to which radius depends on Mass (at constant luminosity)

        // The following variables are implemented as vectors rather than maps due to prohibitive
        // map access times (presumably due to hashing) when accessed hundreds of thousands, and
        // in some cases, millions, of times as we evolve the star/binary.  The code is not as
        // elegant, but performance is better by an order of magnitude.
        //
        // The vectors below are wrapped in std::optional.  In each case, if the vector exists,
        // all elements will exist (so no need for vector to be OPT_DBL_VECTOR).
    
        std::optional<DBL_VECTOR>  aCoefficients;           // Hurley et al. 2000 a(n) coefficients
        std::optional<DBL_VECTOR>  bCoefficients;           // Hurley et al. 2000 b(n) coefficients
    
        std::optional<DBL_VECTOR>  gammaConstants;          // Hurley et al. 2000 gamma constants
        std::optional<DBL_VECTOR>  luminosityConstants;     // Hurley et al. 2000 luminosity constants
        std::optional<DBL_VECTOR>  radiusConstants;         // Hurley et al. 2000 radius constants

        std::optional<DBL_VECTOR>  massCutoffs;             // Hurley et al. 2000 mass cutoffs

        std::optional<DBL_VECTOR>  alphas;                  // Hurley at al. 2000 alpha values (alpha1, alpha3, and alpha4; alpha2 is not constant, so not calculated here)

    } HurleyZdependentT;


    Globals() {

        // calculate the baryonic mass for which the gravitational remnant mass will be equal
        // to the maximum NS mass (see option '--maximum-neutron-star-mass')
        // calculated once per run
        const double maxNSMass    = OPTIONS->MaximumNeutronStarMass();
        m_BaryonicMassOfMaxNSMass = (0.075 * maxNSMass * maxNSMass) + maxNSMass;
    } 
    
    Globals(Globals const&) = delete;
    Globals& operator = (Globals const&) = delete;

    static Globals* m_Instance;


    // reference metallicity - used to calculate metallicity-dependent globals
    // The reference metallicity may not exist - it will not exist prior to the
    // initialisation of the globals, and/or until a valid value is passed via
    // SetReferenceMetallicity().
    //
    // The globals Initialise() function uses the value specified via the 
    // `--metallicity` option as the reference metallicity, and that should
    // have been range checked by the options code, so the reference metallicity
    // should exist and be valid after globals initialisation.
    std::optional<double> m_RefZ;


    // non Z-dependent values
    // these are calculated once per run only (in the constructor), and are guaranteed to
    // exist only if m_RefZ has a value, and then will be correct for the value of m_RefZ 
    
    std::optional<double> m_BaryonicMassOfMaxNSMass;                    // baryonic mass for which the gravitational remnant mass will be equal to the max NS mass


    // Z-dependent values
    // these are calculated whenever the reference metallicity changes, and are guaranteed
    // to exist only if m_RefZ has a value, and then will be correct for the value of m_RefZ 

    // although the following two values are ZAMS values, they can be calculated for all stars because
    // they only depend on the ZAMS metallicity of the star (COMPAS does not change the metallicity of
    // the star throughout its lifetime, and the metallicity given for a star is the ZAMS metallicity)
    std::optional<double> m_ZAMSHeAbundance;                            // Z-dependent ZAMS helium abundance fraction
    std::optional<double> m_ZAMSHAbundance;                             // Z-dependent ZAMS hydrogen abundance fraction

    std::optional<DBL_VECTOR> m_ToutZAMSLuminosityCoefficients;         // Tout luminosity coefficients - if the vector exists, all elements will exist
    std::optional<DBL_VECTOR> m_ToutZAMSRadiusCoefficients;             // Tout radius coefficients - if the vector exists, all elements will exist

   
    // variables defined in Hurley at al. 2000, but possibly used more widely

    HurleyZdependentT m_HurleyZdependentValues;                         // Hurley Z-dependent values - struct always exists; members are wrapped in std::optional<>

    std::optional<double> m_SigmaHurley;                                // Hurley et al. 2000 p24, sigma = log10(Z)

    std::optional<double> m_ZetaHurley;                                 // Hurley et al. 2000 zeta, p5, just before eq 1 (log10(refZ / ZSOL_HURLEY))
    std::optional<double> m_ZetaAnders;                                 // Anders zeta (log10(refZ / ZSOL_ANDERS))
    std::optional<double> m_ZetaAsplund;                                // Asplund zeta (log10(refZ / ZSOL_ASPLUND))

    std::optional<double> m_ZscaledAnders;                              // Z scaled inversely by Anders ZSOL (refZ / ZSOL_ANDERS)
    std::optional<double> m_ZscaledAsplund;                             // Z scaled inversely by Asplund ZSOL (refZ / ZSOL_ASPLUND)
    std::optional<double> m_ZscaledHurley;                              // Z scaled inversely by Hurley ZSOL (refZ / ZSOL_HURLEY)

    
    // ZAMS Z-dependent values
    // these are calculated whenever the reference metallicity changes, and are guaranteed to exist
    // only if m_RefZ has a value *and* the star starts on the main sequence, and then will be correct
    // ZAMS values for the value of m_RefZ 

    // Shikauchi et al. 2024 coefficients - only present if required:
    //    - if the BRCEK MS core mass prescription was specified, and 
    //    - star starts on the main sequence
    // The vector is not guaranteed to exist (only if the above conditions
    // are met), hence the std::optional wrapper, but, in each case, if the
    // vector exists, all elements will exist (so no need for vector to be
    // OPT_DBL_VECTOR).
    std::optional<DBL_VECTOR> m_ShikauchiACoeffs;                       // alpha: the natural decline rate of fMix
    std::optional<DBL_VECTOR> m_ShikauchiFCoeffs;                       // fMix: fraction of the mass contained in themixing core mass at ZAMS
    std::optional<DBL_VECTOR> m_ShikauchiLCoeffs;                       // luminosity


    std::optional<double> m_MinimumLuminosity_CHeB;                     // minimum luminosity on the Core Helium Burning phase


    struct StarDetails {

        double      refZ = -1.0;            // reference metallicty, initially undefined


    };




    // these are per binary - should be in BaseStar  ???? Huh?  That's the point of the globals - per binary or per run <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    double m_ZAMSluminosity;


public:

    static Globals* Instance() {
        if (!m_Instance) m_Instance = new Globals();
        return m_Instance;
    }
    
    void Initialise();


    // getters

    // caller should use:
    //
    // auto x = Peek(n)
    // if (x) or if (x.has_value()) // "x.value().randomSeed" or "x->randomSeed" or "(*x).randomSeed"
    // return value may not have a value (depends on stack size and lookback count) - caler needs to deal with this possibility

    #define FAIL         { std::cerr << "\nGlobals attribute has no value: program terminated\n"; utils::ShowStackTrace(); std::exit(1); }
    #define RET_VALUE(x) { if (x.has_value()) { return x.value(); } else FAIL }

    // all vector variables can be retrieved as a vector, or as individual elements,
    // but there is overhead in retrieving the individual elements - range checking is
    // performed on the specified index (where appropriate).  In some cases it will be
    // more appropriate to retrieve the vector, and let the caller do the range checking
    // as approriate (the Hurley a and b coefficients, for example).

    double     HurleyACoefficients(const size_t p_Index) const { if (p_Index < m_HurleyZdependentValues.aCoefficients.size()) RET_VALUE(m_HurleyZdependentValues.aCoefficients[p_Index]); else FAIL; }
    DBL_VECTOR HurleyACoefficients() const { RET_VALUE(m_HurleyZdependentValues.aCoefficients); }

    double     HurleyAlpha1() const { RET_VALUE(m_HurleyZdependentValues.alphas[0]); }
    double     HurleyAlpha3() const { RET_VALUE(m_HurleyZdependentValues.alphas[1]); }
    double     HurleyAlpha4() const { RET_VALUE(m_HurleyZdependentValues.alphas[2]); }

    double     HurleyBCoefficients(const size_t p_Index) const { if (p_Index < m_HurleyZdependentValues.bCoefficients.size()) RET_VALUE(m_HurleyZdependentValues.bCoefficients[p_Index]); else FAIL; }
    DBL_VECTOR HurleyBCoefficients() const { RET_VALUE(m_HurleyZdependentValues.bCoefficients); }

    double     HurleyGammaConstants(const size_t p_Index) const { if (p_Index < m_HurleyZdependentValues.gammaConstants.size()) RET_VALUE(m_HurleyZdependentValues.gammaConstants[p_Index]); else FAIL; }
    DBL_VECTOR HurleyGammaConstants() const { RET_VALUE(m_HurleyZdependentValues.gammaConstants); }

    double     HurleyLuminosityConstants(const size_t p_Index) const { if (p_Index < m_HurleyZdependentValues.luminosityConstants.size()) RET_VALUE(m_HurleyZdependentValues.luminosityConstants[p_Index]); else FAIL; }
    DBL_VECTOR HurleyLuminosityConstants() const { RET_VALUE(m_HurleyZdependentValues.luminosityConstants); }

    double     HurleyMassCutoffs(const size_t p_Index) const { if (p_Index < m_HurleyZdependentValues.massCutoffs.size()) RET_VALUE(m_HurleyZdependentValues.massCutoffs[p_Index]); else FAIL; }
    DBL_VECTOR HurleyMassCutoffs() const { RET_VALUE(m_HurleyZdependentValues.massCutoffs); }

    double     HurleyRadiusConstants(const size_t p_Index) const { if (p_Index < m_HurleyZdependentValues.radiusConstants.size()) RET_VALUE(m_HurleyZdependentValues.radiusConstants[p_Index]); else FAIL; }
    DBL_VECTOR HurleyRadiusConstants() const { RET_VALUE(m_HurleyZdependentValues.radiusConstants); }

    double     HurleyRadiusXexponent() const { RET_VALUE(m_HurleyZdependentValues.radiusXexponent); }

    double     HurleyRho() const { RET_VALUE(m_HurleyZdependentValues.rho); }

    double     ReferenceMetallicity() const { return RET_VALUE(m_RefZ); } 

    double     ShikauchiACoefficients(const size_t p_Index) const { if (p_Index < m_ShikauchiACoeffs.size()) RET_VALUE(m_ShikauchiACoeffs[p_Index]); else FAIL; }
    DBL_VECTOR ShikauchiACoefficients() const { return UnPackOptDblVector(m_ShikauchiACoeffs); }

    double     ShikauchiFCoefficients(const size_t p_Index) const { if (p_Index < m_ShikauchiFCoeffs.size()) RET_VALUE(m_ShikauchiFCoeffs[p_Index]); else FAIL; }
    DBL_VECTOR ShikauchiFCoefficients() const { return UnPackOptDblVector(m_ShikauchiFCoeffs); }

    double     ShikauchiLCoefficients(const size_t p_Index) const { if (p_Index < m_ShikauchiLCoeffs.size()) RET_VALUE(m_ShikauchiLCoeffs[p_Index]); else FAIL; }
    DBL_VECTOR ShikauchiLCoefficients() const { return UnPackOptDblVector(m_ShikauchiLCoeffs); }

    double     SigmaHurley() const { return RET_VALUE(m_SigmaHurley); }

    double     ToutZAMSLuminosityCoefficients(const size_t p_Index) const { if (p_Index < m_ToutZAMSLuminosityCoefficients.size()) RET_VALUE(m_ToutZAMSLuminosityCoefficients[p_Index]); else FAIL; }
    DBL_VECTOR ToutZAMSLuminosityCoefficients() const { RET_VALUE(m_ToutZAMSLuminosityCoefficients); }

    double     ToutZAMSRadiusCoefficients(const size_t p_Index) const { if (p_Index < m_ToutZAMSRadiusCoefficients.size()) RET_VALUE(m_ToutZAMSRadiusCoefficients[p_Index]); else FAIL; }
    DBL_VECTOR ToutZAMSRadiusCoefficients() const { RET_VALUE(m_ToutZAMSRadiusCoefficients); }

    double     ZAMSHAbundance() const  { return RET_VALUE(m_ZAMSHAbundance); }     // allowed because ZAMS H abundance depends only on ZAMS metallicity (per Pols)
    double     ZAMSHeAbundance() const { return RET_VALUE(m_ZAMSHeAbundance); }    // allowed because ZAMS He abundance depends only on ZAMS metallicity (per Pols)

    double     ZetaAnders() const  { return RET_VALUE(m_ZetaAnders); }
    double     ZetaAsplund() const { return RET_VALUE(m_ZetaAsplund); }
    double     ZetaHurley() const  { return RET_VALUE(m_ZetaHurley); }

    double     ZscaledAnders() const  { return RET_VALUE(m_ZscaledAnders); }
    double     ZscaledAsplund() const { return RET_VALUE(m_ZscaledAsplund); }
    double     ZscaledHurley() const  { return RET_VALUE(m_ZscaledHurley); }

    #undef RET_VALUE


    // setters
    void SetReferenceMetallicity(const double p_ReferenceMetallicity);



    void SetReferenceMetallicity(const double p_ReferenceMetallicity) {
        
        double refZ = p_ReferenceMetallicity;

        // range check and clamp to [MINIMUM_METALLICITY, MAXIMUM_METALLICITY] before anything else
        // this shouldn't happen - options code should already have caught this, but defensive...
        if (refZ < MINIMUM_METALLICITY || refZ > MAXIMUM_METALLICITY) {                         // reference metallicity outside range?
                                                                                                // yes
            refZ = std::max(MAXIMUM_METALLICITY, std::min(refZ, MAXIMUM_METALLICITY));          // clamp it
                                                                                                // issue warning
            // ISSUE WARNING & CLAMP TO [MINIMUM, MAXIMUM]  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        }
        
        if (refZ != m_RefZ) {                                                                   // reference metallicity changed?
                                                                                                
            m_RefZ = refZ;                                                                      // yes - set new reference metallicity

            // (re)calculate metallicity-dependent values

            // the following are evolution-mode independent, or may be used for different evolution modes,
            // so are re calculated whenever the reference metallicity changes

            m_SigmaHurley = log10(m_RefZ);

            m_ZetaAnders  = m_SigmaHurley - LOG10_ZSOL_ANDERS;
            m_ZetaAsplund = m_SigmaHurley - LOG10_ZSOL_ASPLUND;
            m_ZetaHurley  = m_SigmaHurley - LOG10_ZSOL_HURLEY;

            m_ZscaledAnders  = m_RefZ / ZSOL_ANDERS;
            m_ZscaledAsplund = m_RefZ / ZSOL_ASPLUND;
            m_ZscaledHurley  = m_RefZ / ZSOL_HURLEY;

            m_ZAMSHeAbundance = CalculateZAMSHeAbundance_Pols(m_RefZ);
            m_ZAMSHAbundance  = CalculateZAMSHAbundance_Pols(m_RefZ);
                         
            // ZAMS coefficients, per Tout.
            // These vectors are wrapped in std::optional<>, but elements are not
            m_ToutZAMSLuminosityCoefficients = CalculateZAMSLuminosityCoefficients_Tout1996(m_ZetaHurley);
            m_ToutZAMSRadiusCoefficients     = CalculateZAMSRadiusCoefficients_Tout1996(m_ZetaHurley);

            // Shikauchi et al. 2024 coefficients are only calculated if the BRCEK MS core mass
            // prescription was specified, and star starts on the main sequence.
            // These vectors are wrapped in std::optional<>, but elements are not
            if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK && utils::IsOneOf(OPTIONS->StellarType(), MAIN_SEQUENCE)) {
                std::tie(m_ShikauchiACoeffs, m_ShikauchiFCoeffs, m_ShikauchiLyCoeffs) = CalculateShikauchiCoefficients(m_RefZ, m_SigmaHurley);
            }

            // the following are evolution-mode dependent, and are calculated only if required
            Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

                EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
                EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
                    m_HurleyZdependentValues = CalculateHurleyZdependentValues(m_RefZ, m_SigmaHurley, m_ZetaHurley, m_HurleyZdependentValues);
                    break;
        
                default:                                                                        // unknown mode
                    // the only way this can happen is if someone added an EVOLUTION_MODE
                    // and it isn't accounted for in this code.  We should not default here, with or without a warning.
                    // We are here because the user chose a mode this code doesn't account for, and that should
                    // be flagged as an error and result in termination of the evolution of the star or binary.
                    // The correct fix for this is to add code for the missing mode or, if the missing mode is
                    // superfluous, remove it from the option.
    
                    THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                 // throw error
            }
        }
    }



    // member functions


    GNU_PURE  DBL_VECTOR CalculateZAMSLuminosityCoefficients_Tout1996(const double p_Zeta) const;
    GNU_PURE  DBL_VECTOR CalculateZAMSRadiusCoefficients_Tout1996(const double p_Zeta) const;

    GNU_CONST double     CalculateZAMSHeAbundance_Pols(const double p_RefZ) const;
    GNU_CONST double     CalculateZAMSHAbundance_Pols(const double p_RefZ) const;
    
    GNU_PURE  HurleyZdependentT CalculateHurleyZdependentValues(const double p_Z, const double p_Sigma, const double p_Zeta, const HurleyZdependentValues& p_HurleyZdependentValues) const;

    GNU_PURE  DBL_VECTOR CalculateHurleyACoefficients(const double p_RefZ, const double p_Sigma, const double p_Zeta) const;
    GNU_CONST DBL_VECTOR CalculateHurleyAlphas(const DBL_VECTOR& p_MassCutoffs, const DBL_VECTOR& p_aCoeffs, const DBL_VECTOR& p_bCoeffs) const;
    GNU_PURE  DBL_VECTOR CalculateHurleyBCoefficients(
        const double      p_Z,
        const double      p_Sigma,
        const double      p_Zeta,
        const double      p_Rho,
        const DBL_VECTOR& p_MassCutoffs
    ) const;
    GNU_CONST DBL_VECTOR CalculateHurleyGammaConstants(const DBL_VECTOR& p_aCoeffs) const;
    GNU_CONST double     CalculateGBRadiusXexponent(const double p_Zeta) const;
    GNU_CONST DBL_VECTOR CalculateHurleyLuminosityConstants(const DBL_VECTOR& p_aCoeffs) const;
    GNU_CONST DBL_VECTOR CalculateHurleyMassCutoffs(const double p_RefZ, const double p_Zeta) const;
    GNU_CONST DBL_VECTOR CalculateHurleyRadiusConstants(const DBL_VECTOR& p_aCoeffs) const;



    GNU_PURE std::tuple<DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> CalculateShikauchiCoefficients(const double p_Z, const double p_logZ) const;

    /*
     * UnPackOptDblVector
     *
     * @brief
     * Unpacks an OPT_DB_VECTOR to a DBL_VECTOR.
     * Converts the vector of std::optional<<double>> values to a vector of double values.
     * 
     * If any of the values in the std::optional<<double>> do not have a value assigned,
     * this function will fail and terminate the program. 
     * 
     * 
     * std::optional<T> Push(const T p_State)
     * 
     * @param       p_State                         Element to add to the stack (at TOS)
     * @return                                      std::optional<T> object containing removed element (if it exits)
     */
    DBL_VECTOR UnPackOptDblVector(const OPT_DBL_VECTOR p_Vec) {
        DBL_VECTOR unpacked;
        for (size_t idx = 0; idx < p_Vec.size(); idx++) unpacked.push_back(RET_VALUE(p_Vec[idx].value()));
        return unpacked;
    }
};

#endif // __Globals_H__
