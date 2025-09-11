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
    std::optional<double> m_ZAMSheliumAbundance;                        // Z-dependent ZAMS helium abundance fraction
    std::optional<double> m_ZAMShydrogenAbundance;                      // Z-dependent ZAMS hydrogen abundance fraction

    std::optional<DBL_VECTOR> m_ToutZAMSLuminosityCoefficients;         // Tout luminosity coefficients
    std::optional<DBL_VECTOR> m_ToutZAMSRadiusCoefficients;             // Tout radius coefficients

    HurleyZdependentT m_HurleyZdependentValues;                         // Hurley Z-dependent values

    
    // variables defined in Hurley at al. 2000, but possibly used more widely
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

    #define RET_VALUE(x) { if (x.has_value()) { return x.value(); } else { std::cerr << "\nGlobals attribute has no value: program terminated\n"; /*utils::ShowStackTrace();*/ std::exit(1); }}

    double ReferenceMetallicity() const { return RET_VALUE(m_RefZ); } 

    double SigmaHurley() const          { return RET_VALUE(m_SigmaHurley); }

    double ZetaAnders() const           { return RET_VALUE(m_ZetaAnders); }
    double ZetaAsplund() const          { return RET_VALUE(m_ZetaAsplund); }
    double ZetaHurley() const           { return RET_VALUE(m_ZetaHurley); }

    double ZscaledAnders() const        { return RET)VALUE(m_ZscaledAnders); }
    double ZscaledAsplund() const       { return RET)VALUE(m_ZscaledAsplund); }
    double ZscaledHurley() const        { return RET)VALUE(m_ZscaledHurley); }

    DBL_VECTOR ShikauschACoefficients() const      { RET_VALUE(m_ShikauchiACoeffs); }
    DBL_VECTOR ShikauchiFCoefficients() const       { RET_VALUE(m_ShikauchiFCoeffs); }
    DBL_VECTOR ShikauchiLCoefficients() const { RET_VALUE(m_ShikauchiLCoeffs); }
    
    std::optional<DBL_VECTOR> ToutZAMSLuminosityCoefficients() const  { return m_ToutZAMSLuminosityCoefficients; }
    std::optional<DBL_VECTOR> ToutZAMSRadiusCoefficients() const      { return m_ToutZAMSRadiusCoefficients; }

    std::optional<DBL_VECTOR> HurleyACoefficients() const             { return m_HurleyZdependentValues.aCoefficients; }
    std::optional<DBL_VECTOR> HurleyBCoefficients() const             { return m_HurleyZdependentValues.aCoefficients; }
    std::optional<DBL_VECTOR> HurleyGammaConstants() const            { return m_HurleyZdependentValues.gammaConstants; }
    std::optional<DBL_VECTOR> HurleyLuminosityConstants() const       { return m_HurleyZdependentValues.luminosityConstants; }
    std::optional<DBL_VECTOR> HurleyMassCutoffs() const               { return m_HurleyZdependentValues.massCutoffs; }
    std::optional<DBL_VECTOR> HurleyRadiusConstants() const           { return m_HurleyZdependentValues.radiusConstants; }


    #undev RET_VALUE
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

            m_Sigma          = log10(m_RefZ);

            m_ZetaAnders     = m_Sigma - LOG10_ZSOL_ANDERS;
            m_ZetaAsplund    = m_Sigma - LOG10_ZSOL_ASPLUND;
            m_ZetaHurley     = m_Sigma - LOG10_ZSOL_HURLEY;

            m_ZscaledAnders  = m_RefZ / ZSOL_ANDERS;
            m_ZscaledAsplund = m_RefZ / ZSOL_ASPLUND;
            m_ZscaledHurley  = m_RefZ / ZSOL_HURLEY;

            // The following values are evolution-mode independent, and are always recalculated
            // when the reference metallicity changes.
            // These values are not wrapped in std::optional<>
            m_ZAMSheliumAbundance   = CalculateZAMSHeliumAbundance_Pols(m_RefZ);
            m_ZAMShydrogenAbundance = CalculateZAMSHydrogenAbundance_Pols(m_RefZ);
                                                                                            
            m_ToutZAMSLuminosityCoefficients = CalculateZAMSLuminosityCoefficients_Tout(m_ZetaHurley);
            m_ToutZAMSRadiusCoefficients     = CalculateZAMSRadiusCoefficients_Tout(m_ZetaHurley);

            // Shikauchi et al. 2024 coefficients are only calculated if the BRCEK MS core mass
            // prescription was specified, and star starts on the main sequence
            // These values are wrapped in std::optional<>
            if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK &&
                utils::IsOneOf(OPTIONS->StellarType(), MAIN_SEQUENCE) {
                std::tie(m_ShikauchiACoeffs, m_ShikauchiFCoeffs, m_ShikauchiLyCoeffs) = CalculateShikauchiCoefficients(m_RefZ);
            }


            // the following are evolution-mode dependent, and are calculated only if required
            Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

                EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
                EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
                    m_HurleyZdependentValues = CalculateHurleyZdependentValues(m_RefZ, m_Sigma, m_ZetaHurley, m_HurleyZdependentValues);
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


    GNU_PURE  DBL_VECTOR CalculateZAMSLuminosityCoefficients_Tout(const double p_Zeta) const;
    GNU_PURE  DBL_VECTOR CalculateZAMSRadiusCoefficients_Tout(const double p_Zeta) const;

    GNU_CONST double     CalculateZAMSHeliumAbundance_Pols(const double p_RefZ) const;
    GNU_CONST double     CalculateZAMSHydrogenAbundance_Pols(const double p_RefZ) const;
    
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