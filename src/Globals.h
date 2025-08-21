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

    Globals() { m_RefZ = -1.0; }            // reference metallicty, initially undefined   
    
    Globals(Globals const&) = delete;
    Globals& operator = (Globals const&) = delete;

    static Globals* m_Instance;

    double          m_RefZ;                 // reference metallicity - used to calculate metallicity-dependent globals



    // Values of the variables in this struct depend on the star's metallicity only - so these values only need to
    // be calculated once per star (upon creation), but can also be reused if metallicity doesn't change from one
    // star to the next (e.g. in a population run).
    
    // This struct contains variables to be (re)calculated for each star
    
    // This struct contains variables to be calculated per Hurley et al. 2000 & 2002
    // Values of the variables in this struct depend on the star's metallicity only - so these values only need to
    // be calculated once per star (upon creation), but can also be reused if metallicity doesn't change from one
    // star to the next (e.g. in a population run).
    typedef struct {

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
        DBL_VECTOR  luminosityCoefficients;          // Hurley et al. 2000 luminosity coefficients
        DBL_VECTOR  radiusCoefficients;          // Hurley et al. 2000 radius coefficients
    
        DBL_VECTOR  gammaConstants;             // Hurley et al. 2000 gamma constants
        DBL_VECTOR  luminosityConstants;             // Hurley et al. 2000 luminosity constants
        DBL_VECTOR  radiusConstants;             // Hurley et al. 2000 radius constants

        DBL_VECTOR  massCutoffs;            // Hurley et al. 2000 mass cutoffs

        DBL_VECTOR  alphas;                 // Hurley at al. 2000 alpha values (alpha1, alpha3, and alpha4; alpha2 is not constant, so not calculated here)
    } m_HurleyZdependent;
    

    double m_ZAMSluminosity;
    double m_ZAMSheliumAbundance;
    double m_ZAMShydrogenAbundance;

public:

    static Globals* Instance() {
        if (!m_Instance) m_Instance = new Globals();
        return m_Instance;
    }
    
    void            Initialise();
    void            Free();

    // getters
    double          ReferenceMetallicity()                                  { return m_RefZ; } 
    
    // setters
    void            SetReferenceMetallicity(const double p_ReferenceMetallicity);



    void SetReferenceMetallicity(const double p_ReferenceMetallicity) {
        
        double refZ = p_ReferenceMetallicity;

        // range check and clamp to [MINIMUM_METALLICITY, MAXIMUM_METALLICITY] before anything else
        // this shouldn't happen - options code should already have caught this, but defensive...
        if (refZ < MINIMUM_METALLICITY || refZ > MAXIMUM_METALLICITY) {                             // reference metallicity outside range?
                                                                                                    // yes
            refZ = std::max(MAXIMUM_METALLICITY, std::min(refZ, MAXIMUM_METALLICITY));              // clamp it
                                                                                                    // issue warning
            // ISSUE WARNING & CLAMP TO [MINIMUM, MAXIMUM]  
        }
        
        if (refZ == m_RefZ) return;                                                                 // no change - nothing to do

        // reference metallicity has changed - recalculate metallicity-dependent values

        Switch (OPTIONS->Mode()) {                                                          // which evolution mode?

            EVOLUTION_MODE::SSE_HURLEY:                                                     // HURLEY SSE
            EVOLUTION_MODE::BSE_HURLEY:                                                     // HURLEY BSE
            InitialiseHurleyZdependentValues();
                break;

            default: // THIS IS AN ERROR!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        }
    }



    // member functions

    // calculate and set Hurley Z-dependent values
    void CalculateAndSetHurleyMassCutoffs();                                            // mass cutoffs

    void CalculateAndSetHurleyACoefficientsAndConstants();                              // a(n) coefficients, gamma, luminosity, and radius constants
    void CalculateAndSetHurleyBCoefficients();                                          // b(n) coefficients
    void CalculateAndSetHurleyLCoefficients();                                          // luminosity coefficients
    void CalculateAndSetHurleyRCoefficients();                                          // radius coefficients

    void CalculateAndSetHurleyGBRadiusXexponent();                                      // GB radius "x" exponent

    void CalculateAndSetHurleyAlpha1();                                                 // alpha1
    void CalculateAndSetHurleyAlpha3();                                                 // alpha3
    void CalculateAndSetHurleyAlpha4();                                                 // alpha4

    void CalculateAndSetHurleyZdependentValues(const double p_ReferenceMetallicity);    // driver function
    
    // calculate and set ZAMS Z-dependent values
    void CalculateAndSetZAMSHeliumAbundance();                                          // helium abundance
    void CalculateAndSetZAMSHydrogenAbundance();                                        // hydrogen abundance
    void CalculateAndSetZAMSLuminosity();                                               // luminosity
    
};

#endif // __Globals_H__



    //////////////////////////////// DBL_VECTOR  GBParams;               // Hurley et al. 2000 Giant Branch parameters       // NOT IN GLOBALS!!! Per timestep
    //////////////////////////////// DBL_VECTOR  timescales;             // Hurley et al. 2000 timescales                    // NOT IN GLOBALS!!! Per timestep

        // construct vectors of required sizes and initialise each element to default value
        //for (size_t i = 0; i < static_cast<int>(TIMESCALE::COUNT); i++)       m_TimescalesHurley.push_back(DEFAULT_INITIAL_DOUBLE_VALUE);   // NOT IN GLOBALS!!! Per timestep
        //for (size_t i = 0; i < static_cast<int>(GBP::COUNT); i++)             m_GBParamsHurley.push_back(DEFAULT_INITIAL_DOUBLE_VALUE);     // NOT IN GLOBALS!!! Per timestep