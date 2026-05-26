#pragma once

#include "constants.h"


////////////////////////////////////////////////////////////////////////////////////////
//                                                                                    //
//                              GLOBALS Singleton Object                              //
//                                                                                    //
////////////////////////////////////////////////////////////////////////////////////////

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


    Globals() { // Constructor

        // We calculate anything that only needs to be calculated once per run here.

/// <<<< CHECK  THESE - DEFINITION OF ONCE PER RUN!!! THESE DEPENDS ON OPTION VALUES WHICH MAY CHANGE IN GRID FILE  (OR RANGES/SETS) !!!!!!!!!!!!   JR FIX THIS

    } 
    
    Globals() { };
    Globals(Globals const&) = delete;
    Globals& operator = (Globals const&) = delete;

    // Meyers-static singleton state flag. C++17 inline static so the
    // definition lives in this (header-only) class without a separate
    // .cpp file. See IsInstantiated() below.
    inline static bool m_Instantiated = false;



    using HurleyZdependentT = struct {

        // This struct contains variables to be calculated per Hurley et al. 2000 & 2002.
        // The values of the variables in this struct depend on the star's metallicity only, so
        // these values only need to be calculated once per star (upon creation), but can also
        // be reused if metallicity doesn't change from one star to the next (e.g. in a
        // population run).

        OptDblT radiusXexponent;                // Hurley et al. 2000 'x' exponent to which radius depends on Mass (at constant luminosity)
        OptDblT rho;                            // Hurley et al. 2000 rho, p24

        // The following variables are implemented as vectors rather than maps due to prohibitive
        // map access times (presumably due to hashing) when accessed hundreds of thousands, and
        // in some cases, millions, of times as we evolve the star/binary.  The code is not as
        // elegant, but performance is better by at least an order of magnitude.
        //
        // The vectors below are wrapped in std::optional.  In each case, if the vector exists,
        // all elements will exist (so no need for vector to be OptDblVectorT).
    
        Optional_DblVecT alphas;                // Hurley at al. 2000 alpha values (alpha1, alpha3, and alpha4; alpha2 is not constant, so not calculated here)

        Optional_DblVecT aCoefficients;         // Hurley et al. 2000 a(n) coefficients
        Optional_DblVecT bCoefficients;         // Hurley et al. 2000 b(n) coefficients
    
        Optional_DblVecT gammaConstants;        // Hurley et al. 2000 gamma constants
        Optional_DblVecT luminosityConstants;   // Hurley et al. 2000 luminosity constants
        Optional_DblVecT radiusConstants;       // Hurley et al. 2000 radius constants

        Optional_DblVecT massCutoffs;           // Hurley et al. 2000 mass cutoffs
    };


    using StarDetailsT = struct {

        double refZ = -1.0;            // reference metallicty, initially undefined


    };



    // Class member variables

    // Object count.
    // Object ids are used to uniquely identify instantiated objects - used primarily for
    // error printing.
    // This variable keeps a count of served object ids (see getter ObjectId()) - each time
    // an object id is served, the counter is incremented.
    ObjectIdT objectCount;                                  // Instantiated object count

    // Evolution mode.
    // The evolution made (currently SSE_HURLEY or BSE_HURLEY) selected by the user.  This
    // may have been explicitly supplied via the '--mode' option, or left as the default
    // (BSE_HURLEY).
    //
    // Either way, we use m_Mode as a sentinal to determine of the Globals object has been
    // initialised: if the m_Mode variable does not have an assigned value, the Globals object
    // has not been initialised, whereas if m_Mode has an assigned value then the Globals object
    // has been initialised.  Callers can either use the convenience 'Initialised()' getter, or
    // just check the whether the mode has an assigned value via the Mode() getter.
    std::optional<EVOLUTION_MODE> m_Mode;                   // Evolution mode

    // Reference metallicity - used to calculate metallicity-dependent globals.
    // The reference metallicity may not exist - it will not exist prior to the initialisation
    // of the globals, and/or until a valid value is passed via SetReferenceMetallicity().
    //
    // The globals Initialise() function uses the value specified via the  `--metallicity`
    // option as the reference metallicity, and that should have been range checked by the
    // options code, so the reference metallicity should exist and be valid after globals
    // initialisation.
    //
    // Assumption: metallicity is same for both stars in a binary.  
    // The code here could easily be changed to have metallicity different for each star, but
    // COMPAS options currently don't allow that, so for now we use that as a simplification.
    // We may never want to do that, but the point is, should we decide we want to (because
    // someone wants to play), we could.
    OptDblT m_RefZ;                                         // Reference metallicity

    // Non Z-dependent values.
    // These are guaranteed to exist at all times.
    // These depend on the value(s) of various program options (if they didn't they'd be constants
    // calculated in constants.h), and are calculated whenever the values of those options change
    // (mostly - in some cases the variable calculation is small enough that checking to determine
    // if the option values changed is (alomst) as computationally demanding as just recalculating
    // the value of the variable, so in those cases we just recalculate the value of the variable
    // without bothering to check whether the option values have changed).

    // Baryonic mass for which the gravitational remnant mass will be equal to the max NS mass
    // See option '--maximum-neutron-star-mass'
    OptDblT m_BaryonicMassOfMaxMassNS;

    // NS constants.
    OptDblT m_NSDecayMassScale;                             // Mass scale on which magnetic field decays during accretion (Msol) (see --pulsar-magnetic-field-decay-massscale)
    OptDblT m_NSDecayTimeScale;                             // Timescale on which magnetic field decays (Myr) (see --pulsar-magnetic-field-decay-timescale)
    OptDblT m_NSLog10MinMagField;                           // User supplied value for option '--pulsar-minimum-magnetic-field' 
    OptDblT m_NSMagFieldLowerLimit;                         // Minimum pulsar magnetic field (log10(Gauss)) (see --pulsar-minimum-magnetic-field)

    // Z-dependent values.
    // These are calculated whenever the reference metallicity changes, and are guaranteed to
    // exist only if m_RefZ has a value, and then will be correct for the value of m_RefZ.

    // Although the following two values are ZAMS values, they can be calculated for all stars
    // because they only depend on the ZAMS metallicity of the star (COMPAS does not change the
    // metallicity of the star throughout its lifetime, and the metallicity given for a star is
    // the ZAMS metallicity).
    OptDblT m_ZAMSHeAbundance;                              // Z-dependent ZAMS helium abundance fraction
    OptDblT m_ZAMSHAbundance;                               // Z-dependent ZAMS hydrogen abundance fraction

    // Tout luminosity and radius coefficients - if the vector exists, all elements will exist
    // (so no need for each elelemnt to be wrapped in std::optional<>).
    Optional_DblVecT m_ToutZAMSLuminosityCoefficients;      // Tout luminosity coefficients     
    Optional_DblVecT m_ToutZAMSRadiusCoefficients;          // Tout radius coefficients  

    // Variables defined in Hurley at al. 2000 and related/dependent variables.

    // Hurley Z-dependent values struct.
    // The struct always exists - members are wrapped in std::optional<>
    HurleyZdependentT m_HurleyZdependentValues;

    OptDblT m_SigmaHurley;                                  // Hurley et al. 2000 p24, sigma = log10(Z)

    OptDblT m_ZetaAnders;                                   // Anders zeta (log10(refZ / ZSOL_ANDERS))
    OptDblT m_ZetaAsplund;                                  // Asplund zeta (log10(refZ / ZSOL_ASPLUND))
    OptDblT m_ZetaHurley;                                   // Hurley et al. 2000 zeta, p5, just before eq 1 (log10(refZ / ZSOL_HURLEY))

    OptDblT m_ZscaledAnders;                                // Z scaled inversely by Anders ZSOL (refZ / ZSOL_ANDERS)
    OptDblT m_ZscaledAsplund;                               // Z scaled inversely by Asplund ZSOL (refZ / ZSOL_ASPLUND)
    OptDblT m_ZscaledHurley;                                // Z scaled inversely by Hurley ZSOL (refZ / ZSOL_HURLEY)

    
    // ZAMS Z-dependent values.
    // These are calculated whenever the reference metallicity changes, and are guaranteed to
    // exist only if m_RefZ has a value *and* the star starts on the main sequence, and then
    // will be correct ZAMS values for the value of m_RefZ.
    OptDblT m_MinimumLuminosity_CHeB;                       // Minimum luminosity on the Core Helium Burning phase  <<<<<< ???? CHECK ZAMS Z-dependence  JR FIX THIS
    OptDblT m_MinimumAngularFrequency_CHE;      // Minimum angular frequency for Chemically Homogeneous Evolution to occur
     
    // Shikauchi et al. 2024 coefficients - only present if required:
    //    - if the BRCEK MS core mass prescription was specified, *and* 
    //    - the star starts on the main sequence
    // The vector is not guaranteed to exist (only if the above conditions are met), hence the
    // std::optional<> wrapper, but, in each case, if the vector exists, all elements will exist
    // (so no need for each elelemnt to be wrapped in std::optional<>).
    Optional_DblVecT m_ShikauchiACoeffs;                    // alpha: the natural decline rate of fMix
    Optional_DblVecT m_ShikauchiFCoeffs;                    // fMix: fraction of the mass contained in themixing core mass at ZAMS
    Optional_DblVecT m_ShikauchiLCoeffs;                    // Luminosity



    // these are per binary - should be in BaseStar  <<<<<  JR FIX THIS
    OptDblT              m_ZAMSluminosity;


    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< do we need these??????????  JR FIX THIS
    std::optional<StarDetailsT> m_Star1;
    std::optional<StarDetailsT> m_Star2;




    // Member functions

    GNU_PURE  DblVectorT        CalculateACoeffs_Hurley2000(const double p_RefZ, const double p_Sigma, const double p_Zeta) const;
    GNU_CONST DblVectorT        CalculateAlphas_Hurley2000(const double p_MHeF, const DblVectorT& p_bCoeffs) const;
    GNU_PURE  DblVectorT        CalculateBCoeffs_Hurley2000(const double p_Z, const double p_Sigma, const double p_Zeta, const double p_Rho, const DblVectorT& p_MassCutoffs) const;
    GNU_PURE  std::tuple<DblVectorT, DblVectorT, DblVectorT> CalculateCoeffs_Shikauchi2024(const double p_Z, const double p_logZ) const;
    GNU_CONST DblVectorT        CalculateGammaConstants_Hurley2000(const DblVectorT& p_aCoeffs) const;
    GNU_CONST double            CalculateGBRadiusXexponent_Hurley2000(const double p_Zeta) const;
    GNU_PURE  HurleyZdependentT CalculateHurleyZdependentValues(const double p_Z, const double p_Sigma, const double p_Zeta, const HurleyZdependentT& p_Zvalues) const;
    GNU_CONST DblVectorT        CalculateLuminosityConstants_Hurley2000(const DblVectorT& p_aCoeffs) const;
    GNU_CONST DblVectorT        CalculateMassCutoffs_Hurley2000(const double p_RefZ, const double p_Zeta) const;
    GNU_CONST DblVectorT        CalculateRadiusConstants_Hurley2000(const DblVectorT& p_aCoeffs) const;


    /*
     * Fail
     *
     * @brief
     * Prints an error message, shows a stack trace, and terminates the program.
     *
     * 
     * void Fail(const ERROR p_Error, const StrT p_Preamble = "") const
     *  
     * @param       p_Error                         The error to be printed (index into error catalog)
     * @param       p_Preamble                      Preamble string to be printed prior to the error string
     *                                              This replaces the default preamble ("\nError retrieving global value\n")
     */
    void Fail(const ERROR p_Error, const StrT p_Preamble = "") const {

        StrT errStr = "";

        if (p_Preamble.empty()) {                               // Preamble supplied?
            errStr = "\nError retrieving global value\n";       // No - default preamble
        }
        else {                                                  // Yes
            if (p_Preamble[0] != '\n') errStr += "\n";          // Add leading newline if necessary
            errStr += p_Preamble;                               // Add supplied preamble
            if (errStr.back() != '\n') errStr += "\n";          // Add trailing newline if necessary
        }

        // Look for error in error catalog
        std::unordered_map<ERROR, std::tuple<ERROR_SCOPE, std::string>>::const_iterator iter = ERROR_CATALOG.find(p_Error);
        if (iter != ERROR_CATALOG.end()) {                      // Found?
            errStr += std::get<1>(staticIter->second) + "\n";   // Yes - add message to error string
        }
        else {                                                  // No
            errStr += "Uknown error\n";                         // Add generic message
        }

        std::cerr << errStr;                                    // Announce error
        utils::ShowStackTrace();                                // Show stack trace
        std::exit(1);                                           // Terminate program
    }


    /*
     * GetValue
     *
     * @brief
     * Gets the underlying value of a variable wrapped in std::optional<>.
     * If the valiable passed does not have an assigned value, this function calls
     * Fail() (see below) to terminate the program.
     * 
     * 
     * template<typename T>
     * T GetValue(const std::optional<T>& p_Var)
     */
    template<typename T>
    T GetValue(const std::optional<T>& p_Var) {
        if (p_Var.has_value()) return p_Var.value();            // If variable has a value, return it unwrapped
        else Fail(ERROR::NO_VALUE);                             // Otherwise fail
    }


    /*
     * UnPackOptDblVector
     *
     * @brief
     * Unpacks an OptDblVectorT to a DblVectorT.
     * Converts the vector of std::optional<<double>> values to a vector of double values.
     * 
     * If any of the values in the std::optional<<double>> do not have a value assigned,
     * this function will fail and terminate the program. 
     * 
     * 
     * DblVectorT UnPackOptDblVector(const OptDblVectorT p_Vec) const
     * 
     * @param       p_Vec                           The vector to be unpacked
     * @return                                      Unpacked vector
     */
    GNU_CONST inline DblVectorT UnPackOptDblVector(const OptDblVectorT p_Vec) const {
        DblVectorT unpacked;
        for (SizeT idx = 0; idx < p_Vec.size(); idx++) unpacked.push_back(GetValue(p_Vec[idx].value()));
        return unpacked;
    }


public:

    ~Globals() {
        // Meyers-static singleton cleans itself up at program exit:
        // nothing to do here besides clearing m_Instantiated so
        // callers can detect we're past destruction. 
        m_Instantiated = false;
    }

    static Globals* Instance() {

        static Globals self;
        m_Instantiated = true;     // set AFTER `self` is fully constructed
        return &self;
    }


    /*
     * Initialise
     *
     * @brief
     * Initialises member variables as appropriate.
     * 
     * This function initialises the Globals structure, but does not set the values
     * of all member variables.  This function only assumes that the program options
     * have been initialised, and sets the values of member variables that can be
     * determined from option values - with the exception of any program otions that
     * require sampling in order to determine their values (e.g. possibly metallicity
     * and mass depending upon the option values).
     * 
     * This function should be called after program options are initialised and prior
     * to evolution of any star or binary.
     *
     * 
     * void Initialise(const EVOLUTION_MODE p_Mode, const CDOUBLE p_RefZ)
     * void Initialise(const EVOLUTION_MODE p_Mode, const double  p_RefZ)
     * 
     */
    void Initialise(const EVOLUTION_MODE p_Mode, const CDOUBLE p_RefZ) { return Initialise(p_Mode, p_RefZ.Value()); }
    void Initialise(const EVOLUTION_MODE p_Mode, const double  p_RefZ) {

        // Options are expected to have been initialised before this function is called.
        // Initialising the Globals singleton class prior to options being initialised
        // will produce unexpected, and probably, erroneous results.   ///  <<<<<<<<<<<<<< JR MAYBE ADD OPTIONS->Initialised() getter so we can check this?????  JR FIX THIS

        m_Mode = p_Mode;

        // Set reference metallicity and calculate Z-dependent global variables.
        SetReferenceMetallicity(p_RefZ);

        // Calculate the baryonic mass for which the gravitational remnant mass will be equal
        // to the maximum NS mass (see option '--maximum-neutron-star-mass').
        // Calculation is simple enough not to bother checking if OPTIONS->MaximumNeutronStarMass() has changed.
        const double maxMassNS    = OPTIONS->MaximumNeutronStarMass();
        m_BaryonicMassOfMaxMassNS = (0.075 * maxMassNS * maxMassNS) + maxMassNS;

        // Calculate NS constants.
        // See options:
        //    '--pulsar-magnetic-field-decay-massscale'
        //    '--pulsar-magnetic-field-decay-timescale'
        //    '--pulsar-minimum-magnetic-field'

        // Calculations for these two are simple enough not to bother checking if
        // OPTIONS->PulsarMagneticFieldDecayMassscale() or OPTIONS->PulsarMagneticFieldDecayTimescale() have changed.
        m_NSDecayMassScale = OPTIONS->PulsarMagneticFieldDecayMassscale() * MSOL_TO_G;
        m_NSDecayTimeScale = OPTIONS->PulsarMagneticFieldDecayTimescale() * MYR_TO_YEAR * SECONDS_IN_YEAR;

        // Benefit here in checking if OPTIONS->PulsarLog10MinimumMagneticField() has changed.
        if (utils::Compare(OPTIONS->PulsarLog10MinimumMagneticField(), m_NSLog10MinMagField) != 0) {    // Changed materially?
                                                                                                        // Yes
            m_NSLog10MinMagField   = OPTIONS->PulsarLog10MinimumMagneticField();                        // Record new value
            m_NSMagFieldLowerLimit = PPOW(10.0, m_NSLog10MinMagField);                                  // Calculate new lower limit
        }
    }


    // Getters
    //
    // Note that the getters for variables wrapped in std::optional<> have two forms, and each has
    // a parameter that informs the compiler (so the decision is made at compile time, not run time)
    // which variant should be called.
    //
    // The default form (which has the extra parameter set to a default value), returns the *value*
    // of the variable requested.
    //
    // For the default form, if the variable being requested is *not* wrapped in std::optional<>,
    // the value of the variable is returned. If the variable *is* wrapped in std::optional<>, then:
    //    (a) if the variable has a value, the value will be returned, but
    //    (b) if the variable *does not* have a value, an error message and stacktrace will be
    //        displayed, and the program terminated.
    // The reasoning for (b) is that for most circumstances, developers will know when global
    // variables have, or are expected to have, defined values, so can just call the getters
    // and get the value of the variable - and if the variable does not have a value that's 
    // probably either a coding problem, or a corruption/timing issue, that needs to be addressed.
    //
    // For the non-default (optional) form (see below for calling syntax), the variable requested
    // is returned wrapped in std::optional<>.  The caller should then check via *.has_value() if
    // the returned variable has a value assigned.  If the variable being requested is *not* actually
    // declared as std::optional<>, the getter will wrap it in std::optional<> before returning it.
    // Such variables will *always* have a value when returned. Variables that *are* declared as
    // std::optional<> will be returned as is, and may or may not have a value assigned.
    //
    //
    // Vector variables
    // ----------------
    //
    // Some global variables are stored as vectors - e.g. the Hurley 'a' and 'b' coefficients.
    // All vector variables can be retrieved as a vector, or as individual elements, but there
    // is overhead in retrieving the individual elements - range checking is performed on the
    // specified index (where appropriate).  In some cases it will be more appropriate to retrieve
    // the vector, and let the caller do the range checking as approriate
    // (the Hurley a and b coefficients, for example).  See below for calling syntax.
    //
    // When a global vector variable is retrieved as a vector, the default form of the getter will
    // return the vector as type 'std::vector<T>', where 'T' is the type of the underlying variable
    // (i.e. double, int, etc.), whereas the non-default (optional) form of the getter will return
    // the vector as type 'std::optional<std::vector<T>>', where 'T' is the type of the underlying
    // variable.
    // 
    // To facilitate retrieval of single, indexed, elements of vector variables, an alternative
    // form of the getters for vector variables is provided.  This alternative form takes a
    // parameter that is the index into the vector of the single element required (see below).
    // For getters of this form, the single element of the vector is returned as type 'T' (the
    // underlying type, i.e. double, int, etc.) for the default form of the getter, and as type
    // 'std::optional<T>' for the non-default (optional) form of the getter.
    //
    //
    // Usage and calling syntax
    // ------------------------
    //
    // The namespace "get" is declared below to facilitate calling the different getter variants.
    //
    // If the developer expects a global variable to have a value assigned (whether it be a
    // variable declated as std::optional<> or not), they should use the default getter variant.
    // For example:
    //
    //     double refZ = GLOBALS->ReferenceMetallicity()
    //
    // will return the value of the reference metallicty as a non-optional double value.
    // The Globals class member variable m_RefZ is declared as std::optional<double>, and there
    // is a small window as COMPAS is initialising where it hasn't yet been assigned a value.
    // If the call to GLOBALS->ReferenceMetallicity() is made before the value of m_RefZ has been
    // set, the call above will fail as described above (error message, stack trace, and program
    // terminated).  If, however, the call is made after the value of m_RefZ has been set, the call
    // above will return the value of m_RefZ as a variable of type double.
    // Note that using:
    //
    //     double refZ = GLOBALS->ReferenceMetallicity(get::unwrapped)
    //
    // is the same as using 
    //
    //     double refZ = GLOBALS->ReferenceMetallicity()
    //
    // (the extra parameter defaults to "get::unwrapped")
    //
    // On the other hand, if the developer is uncertain if a global variable has a value assigned,
    // they should use the optional getter variant. For example:
    //
    //     std::optional<double> refZ = GLOBALS->ReferenceMetallicity(asOptional)
    //
    // which will return a std::optional variable that can be interrogated as follows:
    //
    //     std::optional<double> refZ = GLOBALS->ReferenceMetallicity(get::asOptional);
    //     if (refZ.has_value()) {
    //         std::cout << "The value of refZ is " << refZ.value() << "\n";
    //     }
    //     else {
    //         // take some appropriate action here if the variable has no value assigned
    //     }
    //
    // Note that the extra parameter must be supplied to distinguish the call from the default
    // variant of the getter ("get::asOptional" must appear as the last parameter in the parameter
    // list).
    //
    // For getters that take a parameter (e.g. getters that return a specific element from a vector),
    // the call syntax is e.g.:
    //
    //    double a21 = HurleyACoefficients(21)
    //
    // for the default getter variant, getting element 21 of the vector, and
    //
    //    std::optional<double> a21 = HurleyACoefficients(21, get::asOptional)
    //
    // for the optional getter variant, getting element 21 of the vector.

    
    static bool IsInstantiated() { return m_Instantiated; }                                             // Instantiated flag

    double HurleyAlpha1(get::unwrappedT = {}) const {
        if (m_HurleyZdependentValues.alphas.has_value())                                                // Alphas vector exists?
            return m_HurleyZdependentValues.alphas[0];                                                  // Yes - all values will exist - return first value unwrapped
        else                                                                                            // No - vector does not exist
            Fail(ERROR::NO_VALUE);                                                                      // Fail - no value assigned
    }
    OptDblT HurleyAlpha1(get::asOptionalT) const {
        if (m_HurleyZdependentValues.alphas.has_value())                                                // Alphas vector exists?
            return std::optional<double>(m_HurleyZdependentValues.alphas[0]);                           // Yes - all values will exist - return first value wrapped
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }
    double HurleyAlpha3(get::unwrappedT = {}) const {
        if (m_HurleyZdependentValues.alphas.has_value())                                                // Alphas vector exists?
            return m_HurleyZdependentValues.alphas[1];                                                  // Yes - all values will exist - return first value unwrapped
        else                                                                                            // No - vector does not exist
            Fail(ERROR::NO_VALUE);                                                                      // Fail - no value assigned
    }
    OptDblT HurleyAlpha3(get::asOptionalT) const {
        if (m_HurleyZdependentValues.alphas.has_value())                                                // Alphas vector exists?
            return std::optional<double>(m_HurleyZdependentValues.alphas[1]);                           // Yes - all values will exist - return first value wrapped
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }
    double HurleyAlpha4(get::unwrappedT = {}) const {
        if (m_HurleyZdependentValues.alphas.has_value())                                                // Alphas vector exists?
            return m_HurleyZdependentValues.alphas[2];                                                  // Yes - all values will exist - return first value unwrapped
        else                                                                                            // No - vector does not exist
            Fail(ERROR::NO_VALUE);                                                                      // Fail - no value assigned
    }
    OptDblT HurleyAlpha4(get::asOptionalT) const {
        if (m_HurleyZdependentValues.alphas.has_value())                                                // Alphas vector exists?
            return std::optional<double>(m_HurleyZdependentValues.alphas[2]);                           // Yes - all values will exist - return first value wrapped
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT HurleyACoefficients(get::unwrappedT = {}) const { return GetValue(m_HurleyZdependentValues.aCoefficients); }
    Optional_DblVecT HurleyACoefficients(get::asOptionalT) const { return m_HurleyZdependentValues.aCoefficients; }
    double HurleyACoefficients(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_HurleyZdependentValues.aCoefficients.has_value()) {                                       // Vector exists?
            if (p_Index < m_HurleyZdependentValues.aCoefficients.size())                                // Yes - index in range?
                return m_HurleyZdependentValues.aCoefficients[p_Index];                                 // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT HurleyACoefficients(const SizeT p_Index, get::asOptionalT) const {
        if (m_HurleyZdependentValues.aCoefficients.has_value()) {                                       // Vector exists?
            if (p_Index < m_HurleyZdependentValues.aCoefficients.size())                                // Yes - index in range?
                return std::optional<double>(m_HurleyZdependentValues.aCoefficients[p_Index]);          // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT HurleyBCoefficients(get::unwrappedT = {}) const { return GetValue(m_HurleyZdependentValues.bCoefficients); }
    Optional_DblVecT HurleyBCoefficients(get::asOptionalT) const { return m_HurleyZdependentValues.bCoefficients; }
    double HurleyBCoefficients(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_HurleyZdependentValues.bCoefficients.has_value()) {                                       // Vector exists?
            if (p_Index < m_HurleyZdependentValues.bCoefficients.size())                                // Yes - index in range?
                return m_HurleyZdependentValues.bCoefficients[p_Index];                                 // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT HurleyBCoefficients(const SizeT p_Index, get::asOptionalT) const {
        if (m_HurleyZdependentValues.bCoefficients.has_value()) {                                       // Vector exists?
            if (p_Index < m_HurleyZdependentValues.bCoefficients.size())                                // Yes - index in range?
                return std::optional<double>(m_HurleyZdependentValues.bCoefficients[p_Index]);          // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT HurleyGammaConstants(get::unwrappedT = {}) const { return GetValue(m_HurleyZdependentValues.gammaConstants); }
    Optional_DblVecT HurleyGammaConstants(get::asOptionalT) const { return m_HurleyZdependentValues.gammaConstants; }
    double HurleyGammaConstants(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_HurleyZdependentValues.gammaConstants.has_value()) {                                      // Vector exists?
            if (p_Index < m_HurleyZdependentValues.gammaConstants.size())                               // Yes - index in range?
                return m_HurleyZdependentValues.gammaConstants[p_Index];                                // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT HurleyGammaConstants(const SizeT p_Index, get::asOptionalT) const {
        if (m_HurleyZdependentValues.gammaConstants.has_value()) {                                      // Vector exists?
            if (p_Index < m_HurleyZdependentValues.gammaConstants.size())                               // Yes - index in range?
                return std::optional<double>(m_HurleyZdependentValues.gammaConstants[p_Index]);         // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT HurleyLuminosityConstants(get::unwrappedT = {}) const { return GetValue(m_HurleyZdependentValues.luminosityConstants); }
    Optional_DblVecT HurleyLuminosityConstants(get::asOptionalT) const { return m_HurleyZdependentValues.luminosityConstants; }
    double HurleyLuminosityConstants(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_HurleyZdependentValues.luminosityConstants.has_value()) {                                 // Vector exists?
            if (p_Index < m_HurleyZdependentValues.luminosityConstants.size())                          // Yes - index in range?
                return m_HurleyZdependentValues.luminosityConstants[p_Index];                           // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT HurleyLuminosityConstants(const SizeT p_Index, get::asOptionalT) const {
        if (m_HurleyZdependentValues.luminosityConstants.has_value()) {                                 // Vector exists?
            if (p_Index < m_HurleyZdependentValues.luminosityConstants.size())                          // Yes - index in range?
                return std::optional<double>(m_HurleyZdependentValues.luminosityConstants[p_Index]);    // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT HurleyMassCutoffs(get::unwrappedT = {}) const { return GetValue(m_HurleyZdependentValues.massCutoffs); }
    Optional_DblVecT HurleyMassCutoffs(get::asOptionalT) const { return m_HurleyZdependentValues.massCutoffs; }
    double HurleyMassCutoffs(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_HurleyZdependentValues.massCutoffs.has_value()) {                                         // Vector exists?
            if (p_Index < m_HurleyZdependentValues.massCutoffs.size())                                  // Yes - index in range?
                return m_HurleyZdependentValues.massCutoffs[p_Index];                                   // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT HurleyMassCutoffs(const SizeT p_Index, get::asOptionalT) const {
        if (m_HurleyZdependentValues.massCutoffs.has_value()) {                                         // Vector exists?
            if (p_Index < m_HurleyZdependentValues.massCutoffs.size())                                  // Yes - index in range?
                return std::optional<double>(m_HurleyZdependentValues.massCutoffs[p_Index]);            // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT HurleyRadiusConstants(get::unwrappedT = {}) const { return GetValue(m_HurleyZdependentValues.radiusConstants); }
    Optional_DblVecT HurleyRadiusConstants(get::asOptionalT) const { return m_HurleyZdependentValues.radiusConstants; }
    double HurleyRadiusConstants(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_HurleyZdependentValues.radiusConstants.has_value()) {                                     // Vector exists?
            if (p_Index < m_HurleyZdependentValues.radiusConstants.size())                              // Yes - index in range?
                return m_HurleyZdependentValues.radiusConstants[p_Index];                               // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT HurleyRadiusConstants(const SizeT p_Index, get::asOptionalT) const {
        if (m_HurleyZdependentValues.radiusConstants.has_value()) {                                     // Vector exists?
            if (p_Index < m_HurleyZdependentValues.radiusConstants.size())                              // Yes - index in range?
                return std::optional<double>(m_HurleyZdependentValues.radiusConstants[p_Index]);        // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT ShikauchiACoefficients(get::unwrappedT = {}) const { return GetValue(m_ShikauchiACoeffs); }
    Optional_DblVecT ShikauchiACoefficients(get::asOptionalT) const { return m_ShikauchiACoeffs; }
    double ShikauchiACoefficients(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_ShikauchiACoeffs.has_value()) {                                                           // Vector exists?
            if (p_Index < m_ShikauchiACoeffs.size())                                                    // Yes - index in range?
                return m_ShikauchiACoeffs[p_Index];                                                     // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT ShikauchiACoefficients(const SizeT p_Index, get::asOptionalT) const {
        if (m_ShikauchiACoeffs.has_value()) {                                                           // Vector exists?
            if (p_Index < m_ShikauchiACoeffs.size())                                                    // Yes - index in range?
                return std::optional<double>(m_ShikauchiACoeffs[p_Index]);                              // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT ShikauchiFCoefficients(get::unwrappedT = {}) const { return GetValue(m_ShikauchiFCoeffs); }
    Optional_DblVecT ShikauchiFCoefficients(get::asOptionalT) const { return m_ShikauchiFCoeffs; }
    double ShikauchiFCoefficients(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_ShikauchiFCoeffs.has_value()) {                                                           // Vector exists?
            if (p_Index < m_ShikauchiFCoeffs.size())                                                    // Yes - index in range?
                return m_ShikauchiFCoeffs[p_Index];                                                     // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT ShikauchiFCoefficients(const SizeT p_Index, get::asOptionalT) const {
        if (m_ShikauchiFCoeffs.has_value()) {                                                           // Vector exists?
            if (p_Index < m_ShikauchiFCoeffs.size())                                                    // Yes - index in range?
                return std::optional<double>(m_ShikauchiFCoeffs[p_Index]);                              // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT ShikauchiLCoefficients(get::unwrappedT = {}) const { return GetValue(m_ShikauchiLCoeffs); }
    Optional_DblVecT ShikauchiLCoefficients(get::asOptionalT) const { return m_ShikauchiLCoeffs; }
    double ShikauchiLCoefficients(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_ShikauchiLCoeffs.has_value()) {                                                           // Vector exists?
            if (p_Index < m_ShikauchiLCoeffs.size())                                                    // Yes - index in range?
                return m_ShikauchiLCoeffs[p_Index];                                                     // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT ShikauchiLCoefficients(const SizeT p_Index, get::asOptionalT) const {
        if (m_ShikauchiLCoeffs.has_value()) {                                                           // Vector exists?
            if (p_Index < m_ShikauchiLCoeffs.size())                                                    // Yes - index in range?
                return std::optional<double>(m_ShikauchiLCoeffs[p_Index]);                              // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT ToutZAMSLuminosityCoefficients(get::unwrappedT = {}) const { return GetValue(m_ToutZAMSLuminosityCoefficients); }
    Optional_DblVecT ToutZAMSLuminosityCoefficients(get::asOptionalT) const { return m_ToutZAMSLuminosityCoefficients; }
    double ToutZAMSLuminosityCoefficients(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_ToutZAMSLuminosityCoefficients.has_value()) {                                             // Vector exists?
            if (p_Index < m_ToutZAMSLuminosityCoefficients.size())                                      // Yes - index in range?
                return m_ToutZAMSLuminosityCoefficients[p_Index];                                       // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT ToutZAMSLuminosityCoefficients(const SizeT p_Index, get::asOptionalT) const {
        if (m_ToutZAMSLuminosityCoefficients.has_value()) {                                             // Vector exists?
            if (p_Index < m_ToutZAMSLuminosityCoefficients.size())                                      // Yes - index in range?
                return std::optional<double>(m_ToutZAMSLuminosityCoefficients[p_Index]);                // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    DblVectorT ToutZAMSRadiusCoefficients(get::unwrappedT = {}) const { return GetValue(m_ToutZAMSRadiusCoefficients); }
    Optional_DblVecT ToutZAMSRadiusCoefficients(get::asOptionalT) const { return m_ToutZAMSRadiusCoefficients; }
    double ToutZAMSRadiusCoefficients(const SizeT p_Index, get::unwrappedT = {}) const {
        if (m_ToutZAMSRadiusCoefficients.has_value()) {                                                 // Vector exists?
            if (p_Index < m_ToutZAMSRadiusCoefficients.size())                                          // Yes - index in range?
                return m_ToutZAMSRadiusCoefficients[p_Index];                                           // Yes - return unwrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else Fail(ERROR::NO_VALUE);                                                                     // No - vector does not exist - fail
    }
    OptDblT ToutZAMSRadiusCoefficients(const SizeT p_Index, get::asOptionalT) const {
        if (m_ToutZAMSRadiusCoefficients.has_value()) {                                                 // Vector exists?
            if (p_Index < m_ToutZAMSRadiusCoefficients.size())                                          // Yes - index in range?
                return std::optional<double>(m_ToutZAMSRadiusCoefficients[p_Index]);                    // Yes - return wrapped value
            else Fail(ERROR::INDEX_OUT_OF_RANGE);                                                       // No - index out of range - fail
        }
        else return std::nullopt;                                                                       // No - vector does not exist - return no value
    }


    double  HurleyRadiusXexponent(get::unwrappedT = {}) const { return GetValue(m_HurleyZdependentValues.radiusXexponent); }
    OptDblT HurleyRadiusXexponent(get::asOptionalT    ) const { return m_HurleyZdependentValues.radiusXexponent; }

    double  HurleyRho(get::unwrappedT = {}) const { return GetValue(m_HurleyZdependentValues.rho); }
    OptDblT HurleyRho(get::asOptionalT    ) const { return m_HurleyZdependentValues.rho; }

    double  NSDecayMassScale(get::unwrappedT = {}) const { return GetValue(m_NSDecayMassScale); }
    OptDblT NSDecayMassScale(get::asOptionalT    ) const { return m_NSDecayMassScale; }

    double  NSDecayTimeScale(get::unwrappedT = {}) const { return GetValue(m_NSDecayTimeScale); }
    OptDblT NSDecayTimeScale(get::asOptionalT    ) const { return m_NSDecayTimeScale; }

    double  NSMagFieldLowerLimit(get::unwrappedT = {}) const { return GetValue(m_NSMagFieldLowerLimit); }
    OptDblT NSMagFieldLowerLimit(get::asOptionalT    ) const { return m_NSMagFieldLowerLimit; }

    double  NSMaxBaryonicMass(get::unwrappedT = {}) const { return GetValue(m_BaryonicMassOfMaxMassNS); }
    OptDblT NSMaxBaryonicMass(get::asOptionalT    ) const { return m_BaryonicMassOfMaxMassNS; }
  
    double  ReferenceMetallicity(get::unwrappedT = {}) const { return GetValue(m_RefZ); } 
    OptDblT ReferenceMetallicity(get::asOptionalT    ) const { return m_RefZ; } 
    double  RefZ(get::unwrappedT = {}) const { return GetValue(m_RefZ); } 
    OptDblT RefZ(get::asOptionalT    ) const { return m_RefZ; } 

    double  SigmaHurley(get::unwrappedT = {}) const { return GetValue(m_SigmaHurley); }
    OptDblT SigmaHurley(get::asOptionalT    ) const { return m_SigmaHurley; }

    double  ZAMSHAbundance(get::unwrappedT = {}) const { return GetValue(m_ZAMSHAbundance); }
    OptDblT ZAMSHAbundance(get::asOptionalT    ) const { return m_ZAMSHAbundance; }

    double  ZAMSHeAbundance(get::unwrappedT = {}) const { return GetValue(m_ZAMSHeAbundance); }
    OptDblT ZAMSHeAbundance(get::asOptionalT    ) const { return m_ZAMSHeAbundance; }

    double  ZetaAnders(get::unwrappedT = {} ) const { return GetValue(m_ZetaAnders); }
    OptDblT ZetaAnders(get::asOptionalT     ) const { return m_ZetaAnders; }
    double  ZetaAsplund(get::unwrappedT = {}) const { return GetValue(m_ZetaAsplund); }
    OptDblT ZetaAsplund(get::asOptionalT    ) const { return m_ZetaAsplund; }
    double  ZetaHurley(get::unwrappedT = {} ) const { return GetValue(m_ZetaHurley); }
    OptDblT ZetaHurley(get::asOptionalT     ) const { return m_ZetaHurley; }

    double  ZscaledAnders(get::unwrappedT = {} ) const { return GetValue(m_ZscaledAnders); }
    OptDblT ZscaledAnders(get::asOptionalT     ) const { return m_ZscaledAnders; }
    double  ZscaledAsplund(get::unwrappedT = {}) const { return GetValue(m_ZscaledAsplund); }
    OptDblT ZscaledAsplund(get::asOptionalT    ) const { return m_ZscaledAsplund; }
    double  ZscaledHurley(get::unwrappedT = {} ) const { return GetValue(m_ZscaledHurley); }
    OptDblT ZscaledHurley(get::asOptionalT     ) const { return m_ZscaledHurley; }


    // Setters

    void SetReferenceMetallicity(const CDOUBLE p_RefZ) { return SetZ(p_RefZ.Value()); }
    void SetReferenceMetallicity(const double  p_RefZ) { return SetZ(p_RefZ); }
    void SetZ(const CDOUBLE p_RefZ) { return SetZ(p_RefZ.Value()); }
    void SetZ(const double  p_RefZ) {

        // Check whether the Globals object has been initialised.  If not, fail here.
        if (!m_Mode.has_value(ERROR::NO_VALUE, "Globals object not initialised")) Fail();

        // If the Globals object has been initialised, then we have a value for m_RefZ.
        // We don't need to do anything if the reference metallicity hasn't (materially) changed.
        // We use utils::Compare() here, with an absolute tolerance of 1E-7.

        if (utils::Compare(p_RefZ, m_RefZ, 1E-7, true) != 0) {                                          // Reference metallicity changed?
                                                                                                        // Yes - set the new value
            // Range check the supplied reference metallicity to [MINIMUM_METALLICITY, MAXIMUM_METALLICITY]
            // before doing anything else.  It *should* be in range - options code should already have done
            // the range check, but we need metallicity to be in range.  If it is out of range here, that's
            // a coding problem (with the options code probably), so we'll fail and terminate. 
            if (p_RefZ < MINIMUM_METALLICITY || p_RefZ > MAXIMUM_METALLICITY) {                         // Reference metallicity outside range?
                                                                                                        // Yes
                Global::Fail(ERROR::OUT_OF_BOUNDS, "Reference metallicity outside allowable range");    // Fail
            }
                                                                                                        // No - proceed             
            m_RefZ = p_RefZ;                                                                            // Set new reference metallicity

            // Since the reference metallicity has changed, we need to recalculate metallicity-dependent
            // values.

            // The following are evolution-mode independent, or may be used for different evolution modes,
            // so are re calculated whenever the reference metallicity changes,

            // Hurley sigma and related variables.
            m_SigmaHurley     = std::log10(m_RefZ);

            m_ZetaAnders      = m_SigmaHurley - LOG10_ZSOL_ANDERS;
            m_ZetaAsplund     = m_SigmaHurley - LOG10_ZSOL_ASPLUND;
            m_ZetaHurley      = m_SigmaHurley - LOG10_ZSOL_HURLEY;

            m_ZscaledAnders   = m_RefZ / ZSOL_ANDERS;
            m_ZscaledAsplund  = m_RefZ / ZSOL_ASPLUND;
            m_ZscaledHurley   = m_RefZ / ZSOL_HURLEY;

            // ZAMS abundances, per Pols.
            m_ZAMSHeAbundance = utils::CalculateZAMSHeAbundance_Pols1998(m_RefZ);
            m_ZAMSHAbundance  = utils::CalculateZAMSHAbundance_Pols1998(m_RefZ);
                         
            // ZAMS coefficients, per Tout.
            // These vectors are wrapped in std::optional<>, but elements are not.
            m_ToutZAMSLuminosityCoefficients = utils::CalculateZAMSLuminosityCoefficients_Tout1996(m_ZetaHurley);
            m_ToutZAMSRadiusCoefficients     = utils::CalculateZAMSRadiusCoefficients_Tout1996(m_ZetaHurley);

            // Shikauchi et al. 2024 coefficients are only calculated if the BRCEK MS core mass
            // prescription was specified, and the star starts on the main sequence.  These vectors
            // are wrapped in std::optional<>, but individual elements are not.
            if (OPTIONS->MainSequenceCoreMassPrescription() == MS_CORE_MASS_PRESCRIPTION::BRCEK) {                                                  // BRCEK MS core mass?
                if ((m_Mode == EVOLUTION_MODE::SSE_HURLEY && utils::IsOneOf(OPTIONS->StellarType(), MAIN_SEQUENCE))  ||                             // Yes - SSE & star MS start?, or
                    (m_Mode == EVOLUTION_MODE::BSE_HURLEY && (utils::IsOneOf(OPTIONS->StellarType1(), MAIN_SEQUENCE) ||                             // BSE & star1 MS start, or
                                                              utils::IsOneOf(OPTIONS->StellarType2(), MAIN_SEQUENCE)))) {                           //       star2 MS start?
                    std::tie(m_ShikauchiACoeffs, m_ShikauchiFCoeffs, m_ShikauchiLyCoeffs) = CalculateCoeffs_Shikauchi2024(m_RefZ, m_SigmaHurley);   // Yes - calculate Shikauchi coeffs
                }
            }


            // The following are evolution-mode dependent, and are calculated only if required.
            Switch (OPTIONS->Mode()) {                                                                  // Which evolution mode?

                EVOLUTION_MODE::SSE_HURLEY:                                                             // HURLEY SSE
                EVOLUTION_MODE::BSE_HURLEY:                                                             // HURLEY BSE
                    m_HurleyZdependentValues = CalculateZdependentValues_Hurley2000(m_RefZ, m_SigmaHurley, m_ZetaHurley, m_HurleyZdependentValues);
                    break;
        
                default:                                                                                // Unknown mode
                    // the only way this can happen is if someone added an EVOLUTION_MODE
                    // and it isn't accounted for in this code.  We should not default here, with or without a warning.
                    // We are here because the user chose a mode this code doesn't account for, and that should
                    // be flagged as an error and result in termination of the evolution of the star or binary.
                    // The correct fix for this is to add code for the missing mode or, if the missing mode is
                    // superfluous, remove it from the option.
    
                    THROW_ERROR(ERROR::UNKNOWN_EVOLUTION_MODE);                                         // Throw error
            }
        }
    }
};
