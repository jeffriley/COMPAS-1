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

        double      refZ = -1.0;            // reference metallicty, initially undefined

        double      sigma;                  // Hurley et al. 2000 sigma, p24
        double      zeta;                   // Hurley et al. 2000 zeta, p5, just before eq 1 (log10(refZ / ZSOL_HURLEY), ZSOL_HURLEY = 0.02)
        double      zetaAnders;             // zeta using ZSOL_ANDERS instead of ZSOL_HURLEY (i.e. log10(refZ / ZSOL_ANDERS))
        double      zetaAsplund;            // zeta using ZSOL_ASPLUND instead of ZSOL_HURLEY (i.e. log10(refZ / ZSOL_ASPLUND))
        double      rho;                    // Hurley et al. 2000 rho, p24

        double      xExponent;              // Hurley et al. 2000 'x' exponent to which radius depends on Mass (at constant luminosity)

        // The following variables are implemented as vectors rather than maps due to prohibitive
        // map access times (presumably due to hashing) when accessed hundreds of thousands, and
        // in some cases, millions, of times as we evolve the star/binary.  The code is not as
        // elegant, but performance is better by an order of magnitude.
   
        DBL_VECTOR  aCoefficients;          // Hurley et al. 2000 a(n) coefficients
        DBL_VECTOR  bCoefficients;          // Hurley et al. 2000 b(n) coefficients
    
        DBL_VECTOR  gammaConstants;         // Hurley et al. 2000 gamma constants
        DBL_VECTOR  luminosityConstants;    // Hurley et al. 2000 luminosity constants
        DBL_VECTOR  radiusConstants;        // Hurley et al. 2000 radius constants

        DBL_VECTOR  massCutoffs;            // Hurley et al. 2000 mass cutoffs

        DBL_VECTOR  alphas;                 // Hurley at al. 2000 alpha values (alpha1, alpha3, and alpha4; alpha2 is not constant, so not calculated here)

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
    std::optional<double> m_RefZ;


    // non Z-dependent values
    // these are calculated once per run only (in the constructor), and are guaranteed to exist
    
    double m_BaryonicMassOfMaxNSMass;                                   // baryonic mass for which the gravitational remnant mass will be equal to the max NS mass


    // Z-dependent values
    // these are calculated whenever the reference metallicity changes, and are guaranteed to exist
    // only if m_RefZ has a value, and then will be correct for the value of m_RefZ 

    // although the following two values are ZAMS values, they can be calculated for all stars because
    // they only depend on the ZAMS metallicity of the star (COMPAS does not change the metallicity of
    // the star throughout its lifetime, and the metallicity given for a star is the ZAMS metallicity)
    double m_ZAMSheliumAbundance;                                       // Z-dependent ZAMS helium abundance fraction
    double m_ZAMShydrogenAbundance;                                     // Z-dependent ZAMS hydrogen abundance fraction

    DBL_VECTOR m_LuminosityCoefficients;                                // luminosity coefficients
    DBL_VECTOR m_RadiusCoefficients;                                    // radius coefficients

    HurleyZdependentT m_HurleyZdependent;                               // Hurley Z-dependent values
   
    // ZAMS Z-dependent values
    // these are calculated whenever the reference metallicity changes, and are guaranteed to exist
    // only if m_RefZ has a value *and* the star starts on the main sequence, and then will be correct
    // ZAMS values for the value of m_RefZ 

    // Shikauchi et al. 2024 coefficients - only present if required:
    // if the BRCEK MS core mass prescription was specified, and star starts on the main sequence
    std::optional<DBL_VECTOR> m_ShikauchiACoeffs;                       // alpha: the natural decline rate of fMix
    std::optional<DBL_VECTOR> m_ShikauchiFCoeffs;                       // fMix: fraction of the mass contained in themixing core mass at ZAMS
    std::optional<DBL_VECTOR> m_ShikauchiLCoeffs;                       // luminosity


    
    struct StarDetails {

        double      refZ = -1.0;            // reference metallicty, initially undefined


    };




    // these are per binary - should be in BaseStar
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


    double        ReferenceMetallicity() const  { return m_RefZ; } 

    DBL_VECTOR    LuminosityCoefficients()      { return m_LuminosityCoefficients; }

    std::optional<DBL_VECTOR> ShikauskiAlphaCoefficients()      { return m_ShikauchiACoeffs; }
    std::optional<DBL_VECTOR> ShikauskifMixCoefficients()       { return m_ShikauchiFCoeffs; }
    std::optional<DBL_VECTOR> ShikauskiLuminosityCoefficients() { return m_ShikauchiLCoeffs; }
    

    // setters
    void            SetReferenceMetallicity(const double p_ReferenceMetallicity);



    void SetReferenceMetallicity(const double p_ReferenceMetallicity) {
        
        double refZ = p_ReferenceMetallicity;

        // range check and clamp to [MINIMUM_METALLICITY, MAXIMUM_METALLICITY] before anything else
        // this shouldn't happen - options code should already have caught this, but defensive...
        if (refZ < MINIMUM_METALLICITY || refZ > MAXIMUM_METALLICITY) {                         // reference metallicity outside range?
                                                                                                // yes
            refZ = std::max(MAXIMUM_METALLICITY, std::min(refZ, MAXIMUM_METALLICITY));          // clamp it
                                                                                                // issue warning
            // ISSUE WARNING & CLAMP TO [MINIMUM, MAXIMUM]  
        }
        
        if (refZ != m_RefZ) {                                                                   // reference metallicity changed?
                                                                                                
            m_RefZ = refZ;                                                                      // yes - set new reference metallicity

            // (re)calculate metallicity-dependent values


            // the following are evolution-mode independent, and are always recalculated
            // when the reference metallicity changes
            m_ZAMSheliumAbundance    = CalculateZAMSHeliumAbundance_Pols_1998(m_RefZ);
            m_ZAMShydrogenAbundance  = CalculateZAMSHydrogenAbundance_Pols_1998(m_RefZ);
                                                                                            
            m_LuminosityCoefficients = CalculateLuminosityCoefficients_Tout_1996(m_HurleyZdependentValues.zeta);
            m_RadiusCoefficients     = CalculateRadiusCoefficients_Tout_1996(m_HurleyZdependentValues.zeta);

            // Shikauchi et al. 2024 coefficients are only calculated if the BRCEK MS core mass
            // prescription was specified, and star starts on the main sequence
            if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK &&
                utils::IsOneOf(OPTIONS->StellarType(), STELLAR_TYPE_LIST::MAIN_SEQUENCE)) {
                std::tie(m_ShikauchiACoeffs, m_ShikauchiFCoeffs, m_ShikauchiLyCoeffs) = CalculateShikauchiCoefficients(m_RefZ);
            }


            // the following are evolution-mode dependent, and are calculated only if required
            Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

                EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
                EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
                    m_HurleyZdependentValues = CalculateHurleyZdependentValues(m_RefZ, m_HurleyZdependentValues);
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


    GNU_PURE  DBL_VECTOR CalculateLuminosityCoefficients_Tout_1996(const double p_Zeta) const;
    GNU_PURE  DBL_VECTOR CalculateRadiusCoefficients_Tout_1996(const double p_Zeta) const;

    GNU_CONST double     CalculateZAMSHeliumAbundance_Pols_1998(const double p_RefZ) const;
    GNU_CONST double     CalculateZAMSHydrogenAbundance_Pols_1998(const double p_RefZ) const;
    
    GNU_PURE  HurleyZdependentT CalculateHurleyZdependentValues(const double p_RefZ);

    GNU_PURE  DBL_VECTOR CalculateHurleyACoefficients(const double p_RefZ, const double p_Sigma, const double p_Zeta) const;
    GNU_CONST DBL_VECTOR CalculateHurleyAlphas(const DBL_VECTOR& p_bCoefficients, const DBL_VECTOR& p_MassCutoffs) const;
    GNU_PURE  DBL_VECTOR CalculateHurleyBCoefficients(const double p_RefZ, const double p_Sigma, const double p_Zeta, const double p_Rho, const DBL_VECTOR& p_MassCutoffs) const;
    GNU_CONST DBL_VECTOR CalculateHurleyGammaConstants(const DBL_VECTOR& p_aCoefficients) const;
    GNU_CONST double     CalculateHurleyGBRadiusXexponent(const double p_Zeta) const;
    GNU_CONST DBL_VECTOR CalculateHurleyLuminosityConstants(const DBL_VECTOR& p_aCoefficients) const;
    GNU_CONST DBL_VECTOR CalculateHurleyMassCutoffs(const double p_RefZ, const double p_Zeta) const;
    GNU_CONST DBL_VECTOR CalculateHurleyRadiusConstants(const DBL_VECTOR& p_aCoefficients) const;

    GNU_PURE std::tuple<DBL_VECTOR, DBL_VECTOR, DBL_VECTOR> CalculateShikauchiCoefficients(const double p_Z) const;
};

#endif // __Globals_H__



    //////////////////////////////// DBL_VECTOR  GBParams;               // Hurley et al. 2000 Giant Branch parameters       // NOT IN GLOBALS!!! Per timestep
    //////////////////////////////// DBL_VECTOR  timescales;             // Hurley et al. 2000 timescales                    // NOT IN GLOBALS!!! Per timestep

        // construct vectors of required sizes and initialise each element to default value
        //for (size_t i = 0; i < static_cast<int>(TIMESCALE::COUNT); i++)       m_TimescalesHurley.push_back(DEFAULT_INITIAL_DOUBLE_VALUE);   // NOT IN GLOBALS!!! Per timestep
        //for (size_t i = 0; i < static_cast<int>(GBP::COUNT); i++)             m_GBParamsHurley.push_back(DEFAULT_INITIAL_DOUBLE_VALUE);     // NOT IN GLOBALS!!! Per timestep