#ifndef __typedefs_h__
#define __typedefs_h__


// common type definitions
// easiest way of making them available globally is to put them here

using StrT              = std::string;
using SizeT             = std::size_t;
using ShortT            = short int;
using LongT             = long int;
using LongLongT         = long long int;
using UIntT             = unsigned int;
using UShortT           = unsigned short int;
using ULongT            = unsigned long int;
using ULongLongT        = unsigned long long int;
using LongDblT          = long double;

using OptStrT           = std::optional<StrT>;
using OptIntT           = std::optional<int>;
using OptBoolT          = std::optional<bool>;
using OptDblT           = std::optional<double>;

using Dbl_DblT          = std::tuple<double, double>;
using Dbl_Dbl_DblT      = std::tuple<double, double, double>;
using Dbl_Dbl_Dbl_DblT  = std::tuple<double, double, double, double>;

using Int_IntT          = std::tuple<int, int>;

using Str_StrT          = std::tuple<StrT, StrT>;
using Str_Str_StrT      = std::tuple<StrT, StrT, StrT>;
using Str_Str_Str_StrT  = std::tuple<StrT, StrT, StrT, StrT>;

using StrVectorT        = std::vector<StrT>;
using DblVectorT        = std::vector<double>;
using IntVectorT        = std::vector<int>;
using SizeTVectorT      = std::vector<SizeT>;
using BoolVectorT       = std::vector<bool>;
using OptDblVectorT     = std::vector<OptDblT>;
using Optional_DblVecT  = std::optional<DblVectorT>; // Hopefully this is sufficiently different so as not to confuse...

using ObjectIdT = ULongT;


// This is where developer-defined types are defined - except for types that pertain directly to
// the COMPAS logging functionality (including the definition of the default record composition
// for the various log files) - those are listed in LogTypedefs.h

// #include <boost/math/tools/roots.hpp>
// #include <boost/numeric/odeint.hpp>

//#include "EnumHash.h"
#include "LogTypedefs.h"
#include "ErrorCatalog.h"



// JR: todo: clean this up and document it better


// Bitwise operators for Enum Class - |, |=, &, &=, ^, ^=, ~ only
// from http://blog.bitwigglers.org/using-enum-classes-as-type-safe-bitmasks/
#define ENABLE_BITMASK_OPERATORS(x)     \
template<>                              \
struct EnableBitMaskOperators<x> {      \
    static const bool enable = true;    \
};

template<typename Enum>  
struct EnableBitMaskOperators {
    static const bool enable = false;
};

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type  
operator |(Enum lhs, Enum rhs) {
    return static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) |
        static_cast<typename std::underlying_type<Enum>::type>(rhs)
    );
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type  
operator |=(Enum &lhs, Enum rhs) {
    lhs = static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) |
        static_cast<typename std::underlying_type<Enum>::type>(rhs)           
    );

    return lhs;
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type 
operator &(Enum lhs, Enum rhs) {
    return static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) &
        static_cast<typename std::underlying_type<Enum>::type>(rhs)
    );
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type 
operator &=(Enum &lhs, Enum rhs) {
    lhs = static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) &
        static_cast<typename std::underlying_type<Enum>::type>(rhs)           
    );

    return lhs;
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type 
operator ^(Enum lhs, Enum rhs) {
    return static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) ^
        static_cast<typename std::underlying_type<Enum>::type>(rhs)
    );
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type 
operator ^=(Enum &lhs, Enum rhs) {
    lhs = static_cast<Enum> (
        static_cast<typename std::underlying_type<Enum>::type>(lhs) ^
        static_cast<typename std::underlying_type<Enum>::type>(rhs)           
    );

    return lhs;
}

template<typename Enum>  
typename std::enable_if<EnableBitMaskOperators<Enum>::enable, Enum>::type 
operator ~(Enum rhs) {
    return static_cast<Enum> (
        ~static_cast<typename std::underlying_type<Enum>::type>(rhs)
    );
}





// put these somewhere appropriate <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 






// enum class types
// ================
//
// Listed alphabetically (with the exception of stellar types - listed first).
// Categories might work, but for now alphabetically makes things easy to find.
//
// The COMPAS nameing convention for enum class is SCREAMING_SNAKE_CASE.
//
// Some enum class have associated maps allowing lookup of description, numerical velue, etc. associated
// with the enum class entries - these maps are define here and listed with the enum class.
//
// The order of entries in most enum classes is not significant - where entries are not set to a specific
// integer value, the integer value for the entry is its ordinal position.  Code should not rely on these
// being/remaining in any specific order.
//
// The value of entries in an enum class must be unique.
//
// Some enum classes are defined in constants.h because they are needed there for constants definition.
//
// some typedefs are listed after enum classes (because they may use enum classes)  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< JR FIX THIS


// We list stellar types and associated initializer lists first so they are grouped and easy to find.
// These are symolic names for the stellar types from Hurley et al. 2000.
//
// The integer value for any stellar type is its ordinal position.  Note that the order of entries is not
// significant - the code should not rely on these being in any order (and so the stellar type symbols
// having any value or order.  e.g. It is not guaranteed that the integer value for stellar type NEUTRON_STAR
// will be > the integer value for stellar type HERTZSPRUNG_GAP).
//
// Initializer lists group stellar types for various tasks performed in COMPAS.
// StellarTypeListT is a type alias for "std::initializer_list<STELLAR_TYPE>"

enum class STELLAR_TYPE: int {               // Hurley
    MS_LTE_07                                 = 0,
    MS_GT_07                                  = 1,
    HERTZSPRUNG_GAP                           = 2,
    HG                                        = 2,
    FIRST_GIANT_BRANCH                        = 3,
    FGB                                       = 3,
    CORE_HELIUM_BURNING                       = 4,
    CHeB                                      = 4,
    EARLY_ASYMPTOTIC_GIANT_BRANCH             = 5,
    EAGB                                      = 5,
    THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH = 6,
    TPAGB                                     = 6,
    NAKED_HELIUM_STAR_MS                      = 7,
    HeMS                                      = 7,
    NAKED_HELIUM_STAR_HERTZSPRUNG_GAP         = 8,
    HeHG                                      = 8,
    NAKED_HELIUM_STAR_GIANT_BRANCH            = 9,
    HeGB                                      = 9,
    HELIUM_WHITE_DWARF                        = 10,
    HeWD                                      = 10,
    CARBON_OXYGEN_WHITE_DWARF                 = 11,
    COWD                                      = 11,
    OXYGEN_NEON_WHITE_DWARF                   = 12,
    ONeWD                                     = 12,
    NEUTRON_STAR                              = 13,
    NS                                        = 13,
    BLACK_HOLE                                = 14,
    BH                                        = 14,
    MASSLESS_REMNANT                          = 15,
    MR                                        = 15,
    CHEMICALLY_HOMOGENEOUS                    = 16, // this is here to preserve the Hurley type numbers, but note that Hurley type number progression doesn't necessarily indicate class inheritance
    CH                                        = 16,

    STAR                                      = 17, // star is created this way, then switches as required (down here so stellar types consistent with Hurley et al. 2000)
    BINARY_STAR                               = 18, // mainly for diagnostics
    STELLAR_BINARY                            = 19, // maily for diagnostics
    NONE                                      = 20, // here mainly for diagnostics

    STELLAR_TYPE_COUNT                        = 17  // count of unique evolvable stellar types
};

const std::unordered_map<STELLAR_TYPE, std::string> STELLAR_TYPE_LABEL = {
    { STELLAR_TYPE::MS_LTE_07,                                 "Main_Sequence_<=_0.7" },
    { STELLAR_TYPE::MS_GT_07,                                  "Main_Sequence_>_0.7" },
    { STELLAR_TYPE::HERTZSPRUNG_GAP,                           "Hertzsprung_Gap" },
    { STELLAR_TYPE::HG,                                        "Hertzsprung_Gap" },
    { STELLAR_TYPE::FIRST_GIANT_BRANCH,                        "First_Giant_Branch" },
    { STELLAR_TYPE::FGB,                                       "First_Giant_Branch" },
    { STELLAR_TYPE::CORE_HELIUM_BURNING,                       "Core_Helium_Burning" },
    { STELLAR_TYPE::CHeB,                                      "Core_Helium_Burning" },
    { STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,             "Early_Asymptotic_Giant_Branch" },
    { STELLAR_TYPE::EAGB,                                      "Early_Asymptotic_Giant_Branch" },
    { STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH, "Thermally_Pulsing_Asymptotic_Giant_Branch" },
    { STELLAR_TYPE::TPAGB,                                     "Thermally_Pulsing_Asymptotic_Giant_Branch" },
    { STELLAR_TYPE::NAKED_HELIUM_STAR_MS,                      "Naked_Helium_Star_MS" },
    { STELLAR_TYPE::HeMS,                                      "Naked_Helium_Star_MS" },
    { STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,         "Naked_Helium_Star_Hertzsprung_Gap" },
    { STELLAR_TYPE::HeHG,                                      "Naked_Helium_Star_Hertzsprung_Gap" },
    { STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH,            "Naked_Helium_Star_Giant_Branch" },
    { STELLAR_TYPE::HeGB,                                      "Naked_Helium_Star_Giant_Branch" },
    { STELLAR_TYPE::HELIUM_WHITE_DWARF,                        "Helium_White_Dwarf" },
    { STELLAR_TYPE::HeWD,                                      "Helium_White_Dwarf" },
    { STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,                 "Carbon-Oxygen_White_Dwarf" },
    { STELLAR_TYPE::COWD,                                      "Carbon-Oxygen_White_Dwarf" },
    { STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF,                   "Oxygen-Neon_White_Dwarf" },
    { STELLAR_TYPE::ONeWD,                                     "Oxygen-Neon_White_Dwarf" },
    { STELLAR_TYPE::NEUTRON_STAR,                              "Neutron_Star" },
    { STELLAR_TYPE::NS,                                        "Neutron_Star" },
    { STELLAR_TYPE::BLACK_HOLE,                                "Black_Hole" },
    { STELLAR_TYPE::BH,                                        "Black_Hole" },
    { STELLAR_TYPE::MASSLESS_REMNANT,                          "Massless_Remnant" },
    { STELLAR_TYPE::MR,                                        "Massless_Remnant" },
    { STELLAR_TYPE::CH,                                        "Chemically_Homogeneous" },
    { STELLAR_TYPE::CH,                                        "Chemically_Homogeneous" },
    { STELLAR_TYPE::STAR,                                      "Star" },
    { STELLAR_TYPE::BINARY_STAR,                               "Binary_Star" },
    { STELLAR_TYPE::NONE,                                      "Not_a_Star!" }
};


enum class STARTING_STELLAR_TYPE: int { // Hurley
    ms,                                 //   0 - preserves Hurley numbering - resolves to MS
    MS,                                 //   1
    HG,                                 //   2
    FGB,                                //   3
    CHeB,                               //   4
    EAGB,                               //   5
    TPAGB,                              //   6
    HeMS,                               //   7
    HeHG,                               //   8
    HeGB,                               //   9
    HeWD,                               //  10
    COWD,                               //  11
    ONeWD,                              //  12
    NS,                                 //  13
    BH                                  //  14
};

const std::unordered_map<STARTING_STELLAR_TYPE, std::string> STARTING_STELLAR_TYPE_LABEL = {
    { STARTING_STELLAR_TYPE::ms,    "MS" },
    { STARTING_STELLAR_TYPE::MS,    "MS" },
    { STARTING_STELLAR_TYPE::HG,    "HG" },
    { STARTING_STELLAR_TYPE::FGB,   "FGB" },
    { STARTING_STELLAR_TYPE::CHeB,  "CHeB" },
    { STARTING_STELLAR_TYPE::EAGB,  "EAGB" },
    { STARTING_STELLAR_TYPE::TPAGB, "TPAGB" },
    { STARTING_STELLAR_TYPE::HeMS,  "HeMS" },
    { STARTING_STELLAR_TYPE::HeHG,  "HeHG" },
    { STARTING_STELLAR_TYPE::HeGB,  "HeGB" },
    { STARTING_STELLAR_TYPE::HeWD,  "HeWD" },
    { STARTING_STELLAR_TYPE::COWD,  "COWD" },
    { STARTING_STELLAR_TYPE::ONeWD, "ONeWD" },
    { STARTING_STELLAR_TYPE::NS,    "NS" },
    { STARTING_STELLAR_TYPE::BH,    "BH" }
};


// (convenience) initializer list for "evolvable" stellar types
// i.e. not STAR, BINARY_STAR, or NONE
const StellarTypeListT EVOLVABLE_TYPES = {
    STELLAR_TYPE::MS_LTE_07,
    STELLAR_TYPE::MS_GT_07,
    STELLAR_TYPE::HERTZSPRUNG_GAP,
    STELLAR_TYPE::FIRST_GIANT_BRANCH,
    STELLAR_TYPE::CORE_HELIUM_BURNING,
    STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::NAKED_HELIUM_STAR_MS,
    STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,
    STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH,
    STELLAR_TYPE::HELIUM_WHITE_DWARF,
    STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,
    STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF,
    STELLAR_TYPE::NEUTRON_STAR,
    STELLAR_TYPE::BLACK_HOLE,
    STELLAR_TYPE::MASSLESS_REMNANT,
    STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS
};


// (convenience) initializer list for MAIN SEQUENCE stars
// (does not include NAKED_HELIUM_STAR_MS)
const StellarTypeListT MAIN_SEQUENCE = {
    STELLAR_TYPE::MS_LTE_07,
    STELLAR_TYPE::MS_GT_07,
    STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS
};


// (convenience) initializer list for ALL MAIN SEQUENCE stars
// (includes NAKED_HELIUM_STAR_MS)
const StellarTypeListT ALL_MAIN_SEQUENCE = {
    STELLAR_TYPE::MS_LTE_07,
    STELLAR_TYPE::MS_GT_07,
    STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS,
    STELLAR_TYPE::NAKED_HELIUM_STAR_MS
};


// (convenience) initializer list for ALL HERTZSPRUNG GAP
// (includes NAKED_HELIUM_STAR_HERTZSPRUNG_GAP)
const StellarTypeListT ALL_HERTZSPRUNG_GAP = {
    STELLAR_TYPE::HERTZSPRUNG_GAP,
    STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP
};


// (convenience) initializer list for non-COMPACT OBJECTS
const StellarTypeListT NON_COMPACT_OBJECTS = {
    STELLAR_TYPE::MS_LTE_07,
    STELLAR_TYPE::MS_GT_07,
    STELLAR_TYPE::CHEMICALLY_HOMOGENEOUS,
    STELLAR_TYPE::HERTZSPRUNG_GAP,
    STELLAR_TYPE::FIRST_GIANT_BRANCH,
    STELLAR_TYPE::CORE_HELIUM_BURNING,
    STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::NAKED_HELIUM_STAR_MS,
    STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,
    STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH
};


// (convenience) initializer list for COMPACT OBJECTS
const StellarTypeListT COMPACT_OBJECTS = {
    STELLAR_TYPE::HELIUM_WHITE_DWARF,
    STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,
    STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF,
    STELLAR_TYPE::NEUTRON_STAR,
    STELLAR_TYPE::BLACK_HOLE,
    STELLAR_TYPE::MASSLESS_REMNANT
};


// (convenience) initializer list for SN REMNANTS
const StellarTypeListT SN_REMNANTS = {
    STELLAR_TYPE::NEUTRON_STAR,
    STELLAR_TYPE::BLACK_HOLE,
    STELLAR_TYPE::MASSLESS_REMNANT
};


// (convenience) initializer list for GIANTS
const StellarTypeListT GIANTS = {
    STELLAR_TYPE::FIRST_GIANT_BRANCH,
    STELLAR_TYPE::CORE_HELIUM_BURNING,
    STELLAR_TYPE::EARLY_ASYMPTOTIC_GIANT_BRANCH,
    STELLAR_TYPE::THERMALLY_PULSING_ASYMPTOTIC_GIANT_BRANCH
};


// (convenience) initializer list for WHITE DWARFS
const StellarTypeListT WHITE_DWARFS = {
    STELLAR_TYPE::HELIUM_WHITE_DWARF,
    STELLAR_TYPE::CARBON_OXYGEN_WHITE_DWARF,
    STELLAR_TYPE::OXYGEN_NEON_WHITE_DWARF
};


// (convenience) initializer list for He rich stellar types
const StellarTypeListT He_RICH_TYPES = {
    STELLAR_TYPE::NAKED_HELIUM_STAR_MS,
    STELLAR_TYPE::NAKED_HELIUM_STAR_HERTZSPRUNG_GAP,
    STELLAR_TYPE::NAKED_HELIUM_STAR_GIANT_BRANCH,
    STELLAR_TYPE::HELIUM_WHITE_DWARF
};


// Start of alphabetical listing

// Accretion regimes
// Symbolic names for WD accretion regimes
enum class ACCRETION_REGIME: int {
    ZERO,
    HELIUM_ACCUMULATION,
    HELIUM_FLASHES,
    HELIUM_STABLE_BURNING,
    HELIUM_OPT_THICK_WINDS,
    HYDROGEN_FLASHES,
    HYDROGEN_STABLE_BURNING,
    HYDROGEN_OPT_THICK_WINDS,
    HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR,
    HELIUM_WHITE_DWARF_HELIUM_IGNITION,
    HELIUM_WHITE_DWARF_HYDROGEN_FLASHES,
    HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION
};

const std::unordered_map<ACCRETION_REGIME, std::string> ACCRETION_REGIME_LABEL = {
    { ACCRETION_REGIME::ZERO,                                        "No accretion" },
    { ACCRETION_REGIME::HELIUM_ACCUMULATION,                         "Helium piles up without burning, full efficiency" },
    { ACCRETION_REGIME::HELIUM_FLASHES,                              "Helium ignites in flashes, partial accretion efficiency" },
    { ACCRETION_REGIME::HELIUM_STABLE_BURNING,                       "Helium is burnt without flashes, full efficiency" },
    { ACCRETION_REGIME::HELIUM_OPT_THICK_WINDS,                      "Helium is being accreted at a high rate, producing winds and limiting the accretion efficiency to a critical value" },
    { ACCRETION_REGIME::HYDROGEN_FLASHES,                            "Hydrogen ignites in flashes, partial accretion efficiency" },
    { ACCRETION_REGIME::HYDROGEN_STABLE_BURNING,                     "Hydrogen is burnt without flashes, full efficiency" },
    { ACCRETION_REGIME::HYDROGEN_OPT_THICK_WINDS,                    "Hydrogen is being accreted at a high rate, producing winds and limiting the accretion efficiency to a critical value" },
    { ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_SUB_CHANDRASEKHAR, "Full accretion leads to transient, but it would not make enough radioactive Ni-56 to be classified as a SN Ia" },
    { ACCRETION_REGIME::HELIUM_WHITE_DWARF_HELIUM_IGNITION,          "Full accretion until material is ignited in a flash and degeneracy is lifted" },
    { ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_FLASHES,         "Unstable hydrogen flashes lead to net accretion being zero" },
    { ACCRETION_REGIME::HELIUM_WHITE_DWARF_HYDROGEN_ACCUMULATION,    "Material piles up. Depending on the companion, it could lead to a CE episode or merger" }
};


// Options to add program option columns to [BSE/SSE] SYSPARMS file
enum class ADD_OPTIONS_TO_SYSPARMS: int { ALWAYS, GRID, NEVER };
const std::unordered_map<ADD_OPTIONS_TO_SYSPARMS, std::string> ADD_OPTIONS_TO_SYSPARMS_LABEL = {
    { ADD_OPTIONS_TO_SYSPARMS::ALWAYS, "ALWAYS" },
    { ADD_OPTIONS_TO_SYSPARMS::GRID,   "GRID" },
    { ADD_OPTIONS_TO_SYSPARMS::NEVER,  "NEVER" }
};


// Black hole kick options
enum class BLACK_HOLE_KICKS_MODE: int { FULL, REDUCED, ZERO, FALLBACK };
const std::unordered_map<BLACK_HOLE_KICKS_MODE, std::string> BLACK_HOLE_KICKS_MODE_LABEL = {
    { BLACK_HOLE_KICKS_MODE::FULL,     "FULL" },
    { BLACK_HOLE_KICKS_MODE::REDUCED,  "REDUCED" },
    { BLACK_HOLE_KICKS_MODE::ZERO,     "ZERO" },
    { BLACK_HOLE_KICKS_MODE::FALLBACK, "FALLBACK" }
};


// BOOST map update options for program options (see options code for use)
enum class BOOST_MAP: int { UPDATE, NO_UPDATE };


// Kick magnitude distributions from Bray & Eldridge 2016,2018
enum class BRAY_ELDRIDGE_CONSTANT: int { ALPHA, BETA };
const std::unordered_map<BRAY_ELDRIDGE_CONSTANT, double> BRAY_ELDRIDGE_CONSTANT_VALUES = {
    { BRAY_ELDRIDGE_CONSTANT::ALPHA, 100.0 },
    { BRAY_ELDRIDGE_CONSTANT::BETA, -170.0 }
};


// Case BB mass transfer stability prescriptions
enum class CASE_BB_STABILITY_PRESCRIPTION: int{ ALWAYS_STABLE, ALWAYS_STABLE_ONTO_NSBH, TREAT_AS_OTHER_MT, ALWAYS_UNSTABLE };
const std::unordered_map<CASE_BB_STABILITY_PRESCRIPTION, std::string> CASE_BB_STABILITY_PRESCRIPTION_LABEL = {
    { CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE,           "ALWAYS_STABLE" },
    { CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_STABLE_ONTO_NSBH, "ALWAYS_STABLE_ONTO_NSBH" },
    { CASE_BB_STABILITY_PRESCRIPTION::TREAT_AS_OTHER_MT,       "TREAT_AS_OTHER_MT" },
    { CASE_BB_STABILITY_PRESCRIPTION::ALWAYS_UNSTABLE,         "ALWAYS_UNSTABLE" }
};


// Common envelope cccretion prescriptions
enum class CE_ACCRETION_PRESCRIPTION: int { ZERO, CONSTANT, UNIFORM, MACLEOD, CHEVALIER };
const std::unordered_map<CE_ACCRETION_PRESCRIPTION, std::string> CE_ACCRETION_PRESCRIPTION_LABEL = {
    { CE_ACCRETION_PRESCRIPTION::ZERO,      "ZERO" },
    { CE_ACCRETION_PRESCRIPTION::CONSTANT,  "CONSTANT" },
    { CE_ACCRETION_PRESCRIPTION::UNIFORM,   "UNIFORM" },
    { CE_ACCRETION_PRESCRIPTION::MACLEOD,   "MACLEOD" },
    { CE_ACCRETION_PRESCRIPTION::CHEVALIER, "CHEVALIER" }
};
    

// Common envelope formalisms
enum class CE_FORMALISM: int { ENERGY, TWO_STAGE };
const std::unordered_map<CE_FORMALISM, std::string> CE_FORMALISM_LABEL = {
    { CE_FORMALISM::ENERGY,    "ENERGY" },
    { CE_FORMALISM::TWO_STAGE, "TWO_STAGE" }
};


// Common envelope lambda prescriptions
enum class CE_LAMBDA_PRESCRIPTION: int { FIXED, LOVERIDGE, NANJING, KRUCKOW, DEWI };
const std::unordered_map<CE_LAMBDA_PRESCRIPTION, std::string> CE_LAMBDA_PRESCRIPTION_LABEL = {
    { CE_LAMBDA_PRESCRIPTION::FIXED,     "LAMBDA_FIXED" },
    { CE_LAMBDA_PRESCRIPTION::LOVERIDGE, "LAMBDA_LOVERIDGE" },
    { CE_LAMBDA_PRESCRIPTION::NANJING,   "LAMBDA_NANJING" },
    { CE_LAMBDA_PRESCRIPTION::KRUCKOW,   "LAMBDA_KRUCKOW" },
    { CE_LAMBDA_PRESCRIPTION::DEWI,      "LAMBDA_DEWI" }
};   


// CHE (Chemically Homogeneous Evolution) Options
enum class CHE_MODE: int { NONE, OPTIMISTIC, PESSIMISTIC };
const std::unordered_map<CHE_MODE, std::string> CHE_MODE_LABEL = {
    { CHE_MODE::NONE,        "NONE" },
    { CHE_MODE::OPTIMISTIC,  "OPTIMISTIC" },
    { CHE_MODE::PESSIMISTIC, "PESSIMISTIC" }
};


// Main sequence core mass prescription
enum class MS_CORE_MASS_PRESCRIPTION: int { BRCEK, HURLEY, MANDEL };
const std::unordered_map<MS_CORE_MASS_PRESCRIPTION, std::string> MS_CORE_MASS_PRESCRIPTION_LABEL = {
    { MS_CORE_MASS_PRESCRIPTION::BRCEK,  "BRCEK" },
    { MS_CORE_MASS_PRESCRIPTION::HURLEY, "HURLEY" },
    { MS_CORE_MASS_PRESCRIPTION::MANDEL, "MANDEL" }
};


// Logfile delimiters
enum class DELIMITER: int { TAB, SPACE, COMMA };
const std::unordered_map<DELIMITER, std::string> DELIMITERLabel = { // labels
    { DELIMITER::TAB,   "TAB" },
    { DELIMITER::SPACE, "SPACE" },
    { DELIMITER::COMMA, "COMMA" }
};
const std::unordered_map<DELIMITER, std::string> DELIMITERValue = { // values
    { DELIMITER::TAB,   "\t" },
    { DELIMITER::SPACE, " " },
    { DELIMITER::COMMA, "," }
};


// Eccentricity distributions
enum class ECCENTRICITY_DISTRIBUTION: int { ZERO, FLAT, THERMAL, GELLER_2013, DUQUENNOYMAYOR1991, SANA2012 };
const std::unordered_map<ECCENTRICITY_DISTRIBUTION, std::string> ECCENTRICITY_DISTRIBUTION_LABEL = {
    { ECCENTRICITY_DISTRIBUTION::ZERO,               "ZERO" },
    { ECCENTRICITY_DISTRIBUTION::FLAT,               "FLAT" },
    { ECCENTRICITY_DISTRIBUTION::THERMAL,            "THERMAL" },
    { ECCENTRICITY_DISTRIBUTION::GELLER_2013,        "GELLER+2013" },
    { ECCENTRICITY_DISTRIBUTION::DUQUENNOYMAYOR1991, "DUQUENNOYMAYOR1991" },
    { ECCENTRICITY_DISTRIBUTION::SANA2012,           "SANA2012"}
};


// Envelope types
enum class ENVELOPE: int { RADIATIVE, CONVECTIVE, REMNANT };
const std::unordered_map<ENVELOPE, std::string> ENVELOPE_LABEL = {
    { ENVELOPE::RADIATIVE,  "RADIATIVE" },
    { ENVELOPE::CONVECTIVE, "CONVECTIVE" },
    { ENVELOPE::REMNANT,    "REMNANT" }
};


// Envelope state prescriptions
enum class ENVELOPE_STATE_PRESCRIPTION: int { LEGACY, HURLEY, FIXED_TEMPERATURE, CONVECTIVE_MASS_FRACTION };
const std::unordered_map<ENVELOPE_STATE_PRESCRIPTION, std::string> ENVELOPE_STATE_PRESCRIPTION_LABEL = {
    { ENVELOPE_STATE_PRESCRIPTION::LEGACY,                   "LEGACY" },
    { ENVELOPE_STATE_PRESCRIPTION::HURLEY,                   "HURLEY" },
    { ENVELOPE_STATE_PRESCRIPTION::FIXED_TEMPERATURE,        "FIXED_TEMPERATURE" },
    { ENVELOPE_STATE_PRESCRIPTION::CONVECTIVE_MASS_FRACTION, "CONVECTIVE_MASS_FRACTION"}
};


// Evolution status constants
enum class EVOLUTION_STATUS: int {
    CONTINUE,
    DONE,
    ERROR,
    TIMES_UP,
    STEPS_UP,
    NO_TIMESTEPS_READ,
    TIMESTEPS_EXHAUSTED,
    TIMESTEPS_NOT_CONSUMED,
    SSE_ERROR,
    BINARY_ERROR,
    DCO_MERGER_TIME,
    STARS_TOUCHING,
    STELLAR_MERGER,
    STELLAR_MERGER_AT_BIRTH,
    DCO,
    WD_WD,
    MASSLESS_REMNANT,
    UNBOUND,
    NOT_STARTED,
    STARTED
};
// These descritions are deliberately succinct (as much as possible) so running status doesn't scroll off the page...
const std::unordered_map<EVOLUTION_STATUS, std::string> EVOLUTION_STATUS_LABEL = {
    { EVOLUTION_STATUS::CONTINUE,                "Continue evolution" },
    { EVOLUTION_STATUS::DONE,                    "Simulation completed" },
    { EVOLUTION_STATUS::ERROR,                   "Evolution stopped because an error occurred" },
    { EVOLUTION_STATUS::TIMES_UP,                "Allowed time exceeded" },
    { EVOLUTION_STATUS::STEPS_UP,                "Allowed timesteps exceeded" },
    { EVOLUTION_STATUS::NO_TIMESTEPS_READ,       "No user-provided timesteps read" },
    { EVOLUTION_STATUS::TIMESTEPS_EXHAUSTED,     "User-provided timesteps exhausted" },
    { EVOLUTION_STATUS::TIMESTEPS_NOT_CONSUMED,  "User-provided timesteps not consumed" },
    { EVOLUTION_STATUS::SSE_ERROR,               "SSE error for one of the constituent stars" },
    { EVOLUTION_STATUS::BINARY_ERROR,            "Error evolving binary" },
    { EVOLUTION_STATUS::DCO_MERGER_TIME,         "Time exceeded DCO merger (formation + coalescence) time" },
    { EVOLUTION_STATUS::STARS_TOUCHING,          "Stars touching" },
    { EVOLUTION_STATUS::STELLAR_MERGER,          "Stars merged" },
    { EVOLUTION_STATUS::STELLAR_MERGER_AT_BIRTH, "Stars merged at birth" },
    { EVOLUTION_STATUS::DCO,                     "DCO formed" },
    { EVOLUTION_STATUS::WD_WD,                   "Double White Dwarf formed" },
    { EVOLUTION_STATUS::MASSLESS_REMNANT,        "Massless Remnant formed" },
    { EVOLUTION_STATUS::UNBOUND,                 "Unbound binary" },
    { EVOLUTION_STATUS::NOT_STARTED,             "Simulation not started" },
    { EVOLUTION_STATUS::STARTED,                 "Simulation started" }
};


// Evolution mode (SSE_HURLEY or BSE_HURLEY)
enum class EVOLUTION_MODE: int { SSE_HURLEY, BSE_HURLEY };
const std::unordered_map<EVOLUTION_MODE, std::string> EVOLUTION_MODE_LABEL = {
    { EVOLUTION_MODE::SSE_HURLEY, "SSE_HURLEY" },
    { EVOLUTION_MODE::BSE_HURLEY, "BSE_HURLEY" }
};


// Floating-point error handling mode
// OFF   specifies that no floating-point error checking is performed
// ON    specifies that the current star/binary will be terminted if a floating-point error occurs
// DEBUG specifies that a stack trace will be printed and the program halted if a floating-point error occurs
enum class FP_ERROR_MODE: int { OFF, ON, DEBUG };
const std::unordered_map<FP_ERROR_MODE, std::string> FP_ERROR_MODE_LABEL = {
    { FP_ERROR_MODE::OFF,   "OFF" },
    { FP_ERROR_MODE::ON,    "ON" },
    { FP_ERROR_MODE::DEBUG, "DEBUG" }
};


// Ge critical mass ratio models
enum class GE_QCRIT_MODEL: int { ST_FULL, ST_HALF, ST_NON_C, IC_FULL, IC_HALF, IC_NON_C};
const std::unordered_map<GE_QCRIT_MODEL, std::string> GE_QCRIT_MODEL_LABEL = {
    { GE_QCRIT_MODEL::ST_FULL,  "ST_FULL" },
    { GE_QCRIT_MODEL::ST_HALF,  "ST_HALF" },
    { GE_QCRIT_MODEL::ST_NON_C, "ST_NON_C" },
    { GE_QCRIT_MODEL::IC_FULL,  "IC_FULL" },
    { GE_QCRIT_MODEL::IC_HALF,  "IC_HALF" },
    { GE_QCRIT_MODEL::IC_NON_C, "IC_NON_C" }
};


// Symbolic names for the Hurley gamma constants.
// These must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end.
// It's a bit of a hack, but it lets us calculate the number of HURLEY_GAMMA_CONSTANTS.
enum class HURLEY_GAMMA_CONSTANTS: int { B_GAMMA, C_GAMMA, COUNT };


// Symbolic names for Giant Branch Parameters.
// These must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end.
// It's a bit of a hack, but it lets us calculate the number of GB parameters.
enum class HURLEY_GB_PARAMETERS: int {
    AH,                     // Hydrogen rate constant.  Hurley et al. 2000, p553
    AHHe,                   // Effective combined rate constant for both hydrogen and helium shell burning.  Hurley et al. 2000, eq 71
    AHe,                    // Helium rate constant.  Hurley et al. 2000, eq 68
    B,                      // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    D,                      // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    p,                      // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    q,                      // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    Lx,                     // Luminosity parameter on the first giant branch (FGB) Lx as a function of the core mass (really a function of Mx).
    Mx,                     // Crosover point of high-luminosity and low-luminosity in core mass - luminosity relation. Hurley et al. 2000, p552, eq38
    McBGB,                  // Core mass at BGB (Base of Giant Branch)
    McBAGB,                 // Core mass at BAGB (Base of Asymptotic Giant Branch).  Hurley et al. 2000, eq 66 (also see eq 75 and discussion)
    McDU,                   // Core mass at second dredge up.  Hurley et al. 2000, eq 69

    COUNT                   // Sentinel for entry count
};


// Symbolic names for the Hurley luminosity constants.
// These must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end.
// It's a bit of a hack, but it lets us calculate the number of HURLEY_L_CONSTANTS.
enum class HURLEY_L_CONSTANTS: int { B_ALPHA_L, B_BETA_L, B_DELTA_L, COUNT };


// Symbolic names for the Hurley mass cutoffs.
// These must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end.
// It's a bit of a hack, but it lets us calculate the number of mass cutoffs.
enum class HURLEY_MASS_CUTOFFS: int {
    Hook,   // Mass above which hook appears on MS (in Msol)
    HeF,    // Maximum initial mass for which helium ignites degenerately in a Helium Flash (HeF)
    FGB,    // Maximum initial mass for which helium ignites on the First Giant Branch (FGB)

    COUNT   // Sentinel for entry count
};


// Symbolic names for the Hurley radius constants.
// These must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end.
// It's a bit of a hack, but it lets us calculate the number of HURLEY_R_CONSTANTS.
enum class HURLEY_R_CONSTANTS: int { B_ALPHA_R, C_ALPHA_R, B_BETA_R, C_BETA_R, B_DELTA_R, COUNT };


// Symbolic names for Hurley timescales.
// These must be left as default values - their order can be changed with the caveat that the sentinel "COUNT" must stay at the end.
// It's a bit of a hack, but it lets us calculate the number of Timescales.
enum class HURLEY_TIMESCALES: int {
    MS,         // Main sequence
    BGB,        // Base of Giant Branch
    HeI,        // Helium ignition
    He,         // Helium burning

                // First Giant Branch (FGB)
    Inf1_FGB,   // First Giant Branch tinf1 (integration constant)
    Inf2_FGB,   // First Giant Branch tinf2 (integration constant)
    Mx_FGB,     // First Giant Branch t(Mx)

                // Early Asymptotic Giant Branch (EAGB) (FAGB in Hurley's sse)
    Inf1_FAGB,  // Early Asymptotic Giant Branch tinf1 (integration constant)
    Inf2_FAGB,  // Early Asymptotic Giant Branch tinf2 (integration constant)
    Mx_FAGB,    // Early Asymptotic Giant Branch t(Mx)

                // Thermally Pulsating Asymptotic Giant Branch (TPAGB) (SAGB in Hurley's sse)
    Inf1_SAGB,  // Thermally Pulsating Asymptotic Giant Branch tinf1 (integration constant)
    Inf2_SAGB,  // Thermally Pulsating Asymptotic Giant Branch tinf2 (integration constant)
    Mx_SAGB,    // Thermally Pulsating Asymptotic Giant Branch t(Mx)

    DU2,        // Lifetime to second dredge-up

                // Helium Giant Branch
    HeMS,       // Naked Helium Star central helium burning lifetime (HeMs)
    Inf1_HeGB,  // Helium Giant Branch tinf1 (integration constant)
    Inf2_HeGB,  // Helium Giant Branch tinf2 (integration constant)
    Mx_HeGB,    // Helium Giant Branch tx (is this t(Mx)?)

    Tau_BL,     // Relative duration of blue loop taubl
    TauX_BL,    // Relative start of blue loop taux
    TauY_BL,    // Relative end of blue loop tauy

    COUNT       // Sentinel for entry count
};


// Immediate events.
// These are events that must be processed immediately, on their own, and in a minimum timestep.
enum class IMMEDIATE_EVENT: int { NONE, RLOF, SWITCH, SUPERNOVA };
const std::unordered_map<IMMEDIATE_EVENT, std::string> IMMEDIATE_EVENT_LABEL = {
    { IMMEDIATE_EVENT::NONE,      "NONE" },
    { IMMEDIATE_EVENT::RLOF,      "ROCHE_LOBE_OVERFLOW" },
    { IMMEDIATE_EVENT::SWITCH,    "STELLAR_TYPE_SWITCH" },
    { IMMEDIATE_EVENT::SUPERNOVA, "SUPERNOVA" }
};


// Initial mass functions
enum class INITIAL_MASS_FUNCTION: int { SALPETER, POWERLAW, UNIFORM, KROUPA };
const std::unordered_map<INITIAL_MASS_FUNCTION, std::string> INITIAL_MASS_FUNCTION_LABEL = {
    { INITIAL_MASS_FUNCTION::SALPETER, "SALPETER" },
    { INITIAL_MASS_FUNCTION::POWERLAW, "POWERLAW" },
    { INITIAL_MASS_FUNCTION::UNIFORM,  "UNIFORM" },
    { INITIAL_MASS_FUNCTION::KROUPA,   "KROUPA" }
};


// Kick magnitude distributions
enum class KICK_MAGNITUDE_DISTRIBUTION: int { ZERO, FIXED, FLAT, MAXWELLIAN, BRAYELDRIDGE, MULLER2016, MULLER2016MAXWELLIAN, MULLERMANDEL, LOGNORMAL};
const std::unordered_map<KICK_MAGNITUDE_DISTRIBUTION, std::string> KICK_MAGNITUDE_DISTRIBUTION_LABEL = {
    { KICK_MAGNITUDE_DISTRIBUTION::ZERO,                 "ZERO" },
    { KICK_MAGNITUDE_DISTRIBUTION::FIXED,                "FIXED" },
    { KICK_MAGNITUDE_DISTRIBUTION::FLAT,                 "FLAT" },
    { KICK_MAGNITUDE_DISTRIBUTION::MAXWELLIAN,           "MAXWELLIAN" },
    { KICK_MAGNITUDE_DISTRIBUTION::BRAYELDRIDGE,         "BRAYELDRIDGE" },
    { KICK_MAGNITUDE_DISTRIBUTION::MULLER2016,           "MULLER2016" },
    { KICK_MAGNITUDE_DISTRIBUTION::MULLER2016MAXWELLIAN, "MULLER2016MAXWELLIAN" },
    { KICK_MAGNITUDE_DISTRIBUTION::MULLERMANDEL,         "MULLERMANDEL" },
    { KICK_MAGNITUDE_DISTRIBUTION::LOGNORMAL,            "LOGNORMAL" }
};


// Kick direction distributions
enum class KICK_DIRECTION_DISTRIBUTION: int { ISOTROPIC, INPLANE, PERPENDICULAR, POWERLAW, WEDGE, POLES };
const std::unordered_map<KICK_DIRECTION_DISTRIBUTION, std::string> KICK_DIRECTION_DISTRIBUTION_LABEL = {
    { KICK_DIRECTION_DISTRIBUTION::ISOTROPIC,     "ISOTROPIC" },
    { KICK_DIRECTION_DISTRIBUTION::INPLANE,       "INPLANE" },
    { KICK_DIRECTION_DISTRIBUTION::PERPENDICULAR, "PERPENDICULAR" },
    { KICK_DIRECTION_DISTRIBUTION::POWERLAW,      "POWERLAW" },
    { KICK_DIRECTION_DISTRIBUTION::WEDGE,         "WEDGE" },
    { KICK_DIRECTION_DISTRIBUTION::POLES,         "POLES" }
};


// LBV mass loss prescriptions
enum class LBV_MASS_LOSS_PRESCRIPTION: int { ZERO, HURLEY_ADD, HURLEY, BELCZYNSKI };
const std::unordered_map<LBV_MASS_LOSS_PRESCRIPTION, std::string> LBV_MASS_LOSS_PRESCRIPTION_LABEL = {
    { LBV_MASS_LOSS_PRESCRIPTION::ZERO,       "ZERO" },
    { LBV_MASS_LOSS_PRESCRIPTION::HURLEY_ADD, "HURLEY_ADD" },
    { LBV_MASS_LOSS_PRESCRIPTION::HURLEY,     "HURLEY" },
    { LBV_MASS_LOSS_PRESCRIPTION::BELCZYNSKI, "BELCZYNSKI" }
};


// Symbolic names for GB groups described in Loveridge et al., 2011.
// These are used as indices into the loveridgeCoefficients multi-dimensional vector (described below).
enum class LOVERIDGE_GROUP: int { LMR1, LMR2, LMA, HM, RECOM };
const std::unordered_map<LOVERIDGE_GROUP, std::string> LOVERIDGE_GROUP_LABEL = {
    { LOVERIDGE_GROUP::LMR1,  "Low mass early Red Giant Branch (RGB) (before dredge-up)" },
    { LOVERIDGE_GROUP::LMR2,  "Low mass late Red Giant Branch (RGB) (after dredge-up)" },
    { LOVERIDGE_GROUP::LMA,   "Low mass Asymptotic Giant Branch (AGB)" },
    { LOVERIDGE_GROUP::HM,    "High mass" },
    { LOVERIDGE_GROUP::RECOM, "Recombination energy" }
};


// Maltsev remnant mass prescription variant
enum class MALTSEV_MODE: int { OPTIMISTIC, BALANCED, PESSIMISTIC };
const std::unordered_map<MALTSEV_MODE, std::string> MALTSEV_MODE_LABEL = {
    { MALTSEV_MODE::OPTIMISTIC,  "OPTIMISTIC" },
    { MALTSEV_MODE::BALANCED,    "BALANCED"   },
    { MALTSEV_MODE::PESSIMISTIC, "PESSIMISTIC"},
};


// Mass loss prescriptions
enum class MASS_LOSS_PRESCRIPTION: int { ZERO, HURLEY, BELCZYNSKI2010, MERRITT2025 };
const std::unordered_map<MASS_LOSS_PRESCRIPTION, std::string> MASS_LOSS_PRESCRIPTION_LABEL = {
    { MASS_LOSS_PRESCRIPTION::ZERO,           "ZERO" },
    { MASS_LOSS_PRESCRIPTION::HURLEY,         "HURLEY" },
    { MASS_LOSS_PRESCRIPTION::BELCZYNSKI2010, "BELCZYNSKI2010" },
    { MASS_LOSS_PRESCRIPTION::MERRITT2025,    "MERRITT2025" }
};


// Mass loss type
enum class MASS_LOSS_TYPE: int { NONE, GB, LBV, OB, RSG, VMS, WR};
const std::unordered_map<MASS_LOSS_TYPE, std::string> MASS_LOSS_TYPE_LABEL = {
    { MASS_LOSS_TYPE::NONE, "NONE" },
    { MASS_LOSS_TYPE::GB,   "GB" },
    { MASS_LOSS_TYPE::LBV,  "LBV" },
    { MASS_LOSS_TYPE::OB,   "OB" },
    { MASS_LOSS_TYPE::RSG,  "RSG" },
    { MASS_LOSS_TYPE::VMS,  "VMS" },
    { MASS_LOSS_TYPE::WR,   "WR" }
};


// Mass ratio distributions
enum class MASS_RATIO_DISTRIBUTION: int { FLAT, DUQUENNOYMAYOR1991, SANA2012 };
const std::unordered_map<MASS_RATIO_DISTRIBUTION, std::string> MASS_RATIO_DISTRIBUTION_LABEL = {
    { MASS_RATIO_DISTRIBUTION::FLAT,               "FLAT" },
    { MASS_RATIO_DISTRIBUTION::DUQUENNOYMAYOR1991, "DUQUENNOYMAYOR1991" },
    { MASS_RATIO_DISTRIBUTION::SANA2012,           "SANA2012" }
};


// Mass transfer timescale types
enum class MT_TIMESCALE: int { NONE, NUCLEAR, THERMAL, CE };
const std::unordered_map<MT_TIMESCALE, std::string> MT_TIMESCALE_LABEL = {
    { MT_TIMESCALE::NONE,    "NONE" },
    { MT_TIMESCALE::NUCLEAR, "NUCLEAR" },
    { MT_TIMESCALE::THERMAL, "THERMAL" },
    { MT_TIMESCALE::CE,      "CE" }
};


// Metallicity distributions
enum class METALLICITY_DISTRIBUTION: int { ZSOLAR, LOGUNIFORM };
const std::unordered_map<METALLICITY_DISTRIBUTION, std::string> METALLICITY_DISTRIBUTION_LABEL = {
    { METALLICITY_DISTRIBUTION::ZSOLAR,     "ZSOLAR" },
    { METALLICITY_DISTRIBUTION::LOGUNIFORM, "LOGUNIFORM" }
};


// Mass transfer accretion efficiency prescriptions
enum class MT_ACCRETION_EFFICIENCY_PRESCRIPTION: int { THERMALLY_LIMITED, FIXED_FRACTION, HAMSTARS };
const std::unordered_map<MT_ACCRETION_EFFICIENCY_PRESCRIPTION, std::string> MT_ACCRETION_EFFICIENCY_PRESCRIPTION_LABEL = {
    { MT_ACCRETION_EFFICIENCY_PRESCRIPTION::THERMALLY_LIMITED, "THERMAL" },
    { MT_ACCRETION_EFFICIENCY_PRESCRIPTION::FIXED_FRACTION,    "FIXED" },
    { MT_ACCRETION_EFFICIENCY_PRESCRIPTION::HAMSTARS,          "HAMSTARS"}
};


// Mass transfer angular momentum loss prescriptions
enum class MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION: int { JEANS, ISOTROPIC_RE_EMISSION, CIRCUMBINARY_RING, MACLEOD_LINEAR, KLENCKI_LINEAR, ARBITRARY };
const std::unordered_map<MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION, std::string> MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION_LABEL = {
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::JEANS,                 "JEANS" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ISOTROPIC_RE_EMISSION, "ISOTROPIC" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::CIRCUMBINARY_RING,     "CIRCUMBINARY" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::MACLEOD_LINEAR,        "MACLEOD_LINEAR" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::KLENCKI_LINEAR,        "KLENCKI_LINEAR" },
    { MT_ANGULAR_MOMENTUM_LOSS_PRESCRIPTION::ARBITRARY,             "ARBITRARY" }
};


// Mass transfer cases
enum class MT_CASE: int { NONE, A, B, C, OTHER };
const std::unordered_map<MT_CASE, std::string> MT_CASE_LABEL = {
    { MT_CASE::NONE, "Mass Transfer CASE NONE: No Mass Transfer" },
    { MT_CASE::A,    "Mass Transfer CASE A" },                          // mass transfer while donor is on main sequence
    { MT_CASE::B,    "Mass Transfer CASE B" },                          // donor star is in (or evolving to) Red Giant phase
    { MT_CASE::C,    "Mass Transfer CASE C" },                          // SuperGiant phase
    { MT_CASE::OTHER,"Mass Transfer CASE OTHER: Multiple MT events" }   // default value, or multiple MT events
};


// Mass transfer rejuvenation prescriptions
enum class MT_REJUVENATION_PRESCRIPTION: int { HURLEY, STARTRACK };
const std::unordered_map<MT_REJUVENATION_PRESCRIPTION, std::string> MT_REJUVENATION_PRESCRIPTION_LABEL = {
    { MT_REJUVENATION_PRESCRIPTION::HURLEY,    "HURLEY" },
    { MT_REJUVENATION_PRESCRIPTION::STARTRACK, "STARTRACK" }
};


// Mass transfer thermally limited variation options
enum class MT_THERMALLY_LIMITED_VARIATION: int { C_FACTOR, RADIUS_TO_ROCHELOBE };
const std::unordered_map<MT_THERMALLY_LIMITED_VARIATION, std::string> MT_THERMALLY_LIMITED_VARIATION_LABEL = {
    { MT_THERMALLY_LIMITED_VARIATION::C_FACTOR,            "CFACTOR" },
    { MT_THERMALLY_LIMITED_VARIATION::RADIUS_TO_ROCHELOBE, "ROCHELOBE" }
};


// Mass transfer timing options for writing to BSE_RLOF file
enum class MT_TIMING: int { PRE_MT, POST_MT };
const std::unordered_map<MT_TIMING, std::string> MT_TIMING_LABEL = {
    { MT_TIMING::PRE_MT,  "PRE_MT" },
    { MT_TIMING::POST_MT, "POST_MT" }
};


// Symbolic names for the mass transfer tracking constants
enum class MT_TRACKING: int { NO_MASS_TRANSFER, STABLE_1_TO_2_SURV, STABLE_2_TO_1_SURV, CE_1_TO_2_SURV, CE_2_TO_1_SURV, CE_DOUBLE_SURV, MERGER }; 
const std::unordered_map<MT_TRACKING, std::string> MT_TRACKING_LABEL = {
    { MT_TRACKING::NO_MASS_TRANSFER,   "NO MASS TRANSFER" },
    { MT_TRACKING::STABLE_1_TO_2_SURV, "MASS TRANSFER STABLE STAR1 -> STAR2" },
    { MT_TRACKING::STABLE_2_TO_1_SURV, "MASS TRANSFER STABLE STAR2 -> STAR1" },
    { MT_TRACKING::CE_1_TO_2_SURV,     "MASS TRANSFER COMMON ENVELOPE STAR1 -> STAR2" },
    { MT_TRACKING::CE_2_TO_1_SURV,     "MASS TRANSFER COMMON ENVELOPE STAR2 -> STAR1" },
    { MT_TRACKING::CE_DOUBLE_SURV,     "MASS TRANSFER COMMON ENVELOPE DOUBLE CORE" },
    { MT_TRACKING::MERGER,             "MASS TRANSFER -> MERGER" }
};


// Neutrino mass loss BH formation prescriptions
enum class NEUTRINO_MASS_LOSS_PRESCRIPTION: int { FIXED_FRACTION, FIXED_MASS };
const std::unordered_map<NEUTRINO_MASS_LOSS_PRESCRIPTION, std::string> NEUTRINO_MASS_LOSS_PRESCRIPTION_LABEL = {
    { NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_FRACTION, "FIXED_FRACTION" },
    { NEUTRINO_MASS_LOSS_PRESCRIPTION::FIXED_MASS,     "FIXED_MASS" }
};


// Neutron star accretion scenario under common envelope
enum class NS_ACCRETION_IN_CE: int { ZERO, SURFACE, DISK };
const std::unordered_map<NS_ACCRETION_IN_CE, std::string> NS_ACCRETION_IN_CE_LABEL = {
    { NS_ACCRETION_IN_CE::ZERO,    "ZERO" },
    { NS_ACCRETION_IN_CE::SURFACE, "SURFACE" },
    { NS_ACCRETION_IN_CE::DISK,    "DISK" },
};


// Neutron star equations of state
enum class NS_EOS: int { SSE, ARP3 };
const std::unordered_map<NS_EOS, std::string> NS_EOS_LABEL = {
    { NS_EOS::SSE,  "SSE" },
    { NS_EOS::ARP3, "ARP3" }
};


// OB (main sequence) mass loss prescriptions
enum class OB_MASS_LOSS_PRESCRIPTION: int { ZERO, VINK2001, VINK2021, BJORKLUND2022, KRTICKA2018};
const std::unordered_map<OB_MASS_LOSS_PRESCRIPTION, std::string> OB_MASS_LOSS_PRESCRIPTION_LABEL = {
    { OB_MASS_LOSS_PRESCRIPTION::ZERO,          "ZERO" },
    { OB_MASS_LOSS_PRESCRIPTION::VINK2001,      "VINK2001" },
    { OB_MASS_LOSS_PRESCRIPTION::VINK2021,      "VINK2021" },
    { OB_MASS_LOSS_PRESCRIPTION::BJORKLUND2022, "BJORKLUND2022" },
    { OB_MASS_LOSS_PRESCRIPTION::KRTICKA2018,   "KRTICKA2018" }
};


// Object persistence.
// Specifies whether an object is permanent or ephemaeral.
// EPHEMERAL is typically used for clones so that they don't participate in logging, etc.
enum class OBJECT_PERSISTENCE: int { PERMANENT, EPHEMERAL };
const std::unordered_map<OBJECT_PERSISTENCE, std::string> OBJECT_PERSISTENCE_LABEL = {
    { OBJECT_PERSISTENCE::PERMANENT, "Permanent" },
    { OBJECT_PERSISTENCE::EPHEMERAL, "Ephemeral" }
};


// Object types.
// Identifies the type of an object.
// If BASE_STAR, check STELLAR_TYPE.    
enum class OBJECT_TYPE: int { NONE, MAIN, PROFILING, UTILS, STAR, BASE_STAR, BINARY_STAR, BASE_BINARY_STAR, BINARY_CONSTITUENT_STAR };
const std::unordered_map<OBJECT_TYPE, std::string> OBJECT_TYPE_LABEL = {
    { OBJECT_TYPE::NONE,                    "Not_an_Object!" },
    { OBJECT_TYPE::MAIN,                    "Main" },
    { OBJECT_TYPE::PROFILING,               "Profiling" },
    { OBJECT_TYPE::UTILS,                   "Utils" },
    { OBJECT_TYPE::STAR,                    "Star" },
    { OBJECT_TYPE::BASE_STAR,               "BaseStar" },
    { OBJECT_TYPE::BINARY_STAR,             "BinaryStar" },
    { OBJECT_TYPE::BASE_BINARY_STAR,        "BaseBinaryStar" },
    { OBJECT_TYPE::BINARY_CONSTITUENT_STAR, "BinaryConstituentStar" }
};


// Program options origin indicator (command line or gridfile line)
enum class OPTIONS_ORIGIN: int { CMDLINE, GRIDFILE };


// Orbital period distributions
enum class ORBITAL_PERIOD_DISTRIBUTION: int { FLATINLOG };
const std::unordered_map<ORBITAL_PERIOD_DISTRIBUTION, std::string> ORBITAL_PERIOD_DISTRIBUTION_LABEL = {
    { ORBITAL_PERIOD_DISTRIBUTION::FLATINLOG, "FLATINLOG" },
};


// Pulsational pair instability prescriptions
enum class PPI_PRESCRIPTION: int { WOOSLEY, STARTRACK, MARCHANT, FARMER, HENDRIKS };
const std::unordered_map<PPI_PRESCRIPTION, std::string> PPI_PRESCRIPTION_LABEL = {
    { PPI_PRESCRIPTION::WOOSLEY,   "WOOSLEY" },
    { PPI_PRESCRIPTION::STARTRACK, "STARTRACK" },
    { PPI_PRESCRIPTION::MARCHANT,  "MARCHANT" },
    { PPI_PRESCRIPTION::FARMER,    "FARMER" },
    { PPI_PRESCRIPTION::HENDRIKS,  "HENDRIKS" }
};


// Program status
enum class PROGRAM_STATUS: int { SUCCESS, CONTINUE, STOPPED, ERROR_IN_COMMAND_LINE, LOGGING_FAILED, ERROR_UNHANDLED_EXCEPTION };


// Pulsar birth magnetic field distributions
enum class PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION: int { FLATINLOG, UNIFORM, LOGNORMAL };
const std::unordered_map<PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION, std::string> PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION_LABEL = {
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::FLATINLOG, "FLATINLOG" },
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::UNIFORM,   "UNIFORM" },
    { PULSAR_BIRTH_MAGNETIC_FIELD_DISTRIBUTION::LOGNORMAL, "LOGNORMAL" }
};


// Pulsar birth spin period distributions
enum class PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION: int { UNIFORM, NORMAL };
const std::unordered_map<PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION, std::string> PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION_LABEL = {
    { PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::UNIFORM, "UNIFORM" },
    { PULSAR_BIRTH_SPIN_PERIOD_DISTRIBUTION::NORMAL,  "NORMAL" }
};


// Critical mass ratio prescriptions
enum class QCRIT_PRESCRIPTION: int { NONE, CLAEYS, GE, GE_IC, HURLEY_HJELLMING_WEBBINK};
const std::unordered_map<QCRIT_PRESCRIPTION, std::string> QCRIT_PRESCRIPTION_LABEL = {
    { QCRIT_PRESCRIPTION::NONE,                     "NONE" },
    { QCRIT_PRESCRIPTION::CLAEYS,                   "CLAEYS" },
    { QCRIT_PRESCRIPTION::GE_IC,                    "GE_IC" },
    { QCRIT_PRESCRIPTION::GE,                       "GE" },
    { QCRIT_PRESCRIPTION::HURLEY_HJELLMING_WEBBINK, "HURLEY_HJELLMING_WEBBINK" }
};


// remnant mass prescriptions
enum class REMNANT_MASS_PRESCRIPTION: int { HURLEY2000, BELCZYNSKI2002, FRYER2012, FRYER2022, MULLER2016, MULLERMANDEL, SCHNEIDER2020, SCHNEIDER2020ALT, MALTSEV2025};
const std::unordered_map<REMNANT_MASS_PRESCRIPTION, std::string> REMNANT_MASS_PRESCRIPTION_LABEL = {
    { REMNANT_MASS_PRESCRIPTION::HURLEY2000,       "HURLEY2000" },
    { REMNANT_MASS_PRESCRIPTION::BELCZYNSKI2002,   "BELCZYNSKI2002" },
    { REMNANT_MASS_PRESCRIPTION::FRYER2012,        "FRYER2012" },
    { REMNANT_MASS_PRESCRIPTION::FRYER2022,        "FRYER2022" },
    { REMNANT_MASS_PRESCRIPTION::MULLER2016,       "MULLER2016" },
    { REMNANT_MASS_PRESCRIPTION::MULLERMANDEL,     "MULLERMANDEL" },
    { REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020,    "SCHNEIDER2020" },
    { REMNANT_MASS_PRESCRIPTION::SCHNEIDER2020ALT, "SCHNEIDER2020ALT" },
    { REMNANT_MASS_PRESCRIPTION::MALTSEV2025,      "MALTSEV2025" }
};


// Response of star to spin-up beyond the Keplerian frequency
enum class RESPONSE_TO_SPIN_UP: int { TRANSFER_TO_ORBIT, KEPLERIAN_LIMIT, NO_LIMIT };
const std::unordered_map<RESPONSE_TO_SPIN_UP, std::string> RESPONSE_TO_SPIN_UP_LABEL = {
    { RESPONSE_TO_SPIN_UP::TRANSFER_TO_ORBIT, "TRANSFER_TO_ORBIT" },
    { RESPONSE_TO_SPIN_UP::KEPLERIAN_LIMIT,   "KEPLERIAN_LIMIT" },
    { RESPONSE_TO_SPIN_UP::NO_LIMIT,          "NO_LIMIT"}
};


// Rotational velocity distributions
enum class ROTATIONAL_VELOCITY_DISTRIBUTION: int { ZERO, HURLEY, VLTFLAMES };
const std::unordered_map<ROTATIONAL_VELOCITY_DISTRIBUTION, std::string> ROTATIONAL_VELOCITY_DISTRIBUTION_LABEL = {
    { ROTATIONAL_VELOCITY_DISTRIBUTION::ZERO,      "ZERO" },
    { ROTATIONAL_VELOCITY_DISTRIBUTION::HURLEY,    "HURLEY" },
    { ROTATIONAL_VELOCITY_DISTRIBUTION::VLTFLAMES, "VLTFLAMES" }
};


// RSG mass loss prescriptions
enum class RSG_MASS_LOSS_PRESCRIPTION: int { ZERO, VINKSABHAHIT2023, BEASOR2020, DECIN2023, YANG2023, KEE2021, NJ90};
const std::unordered_map<RSG_MASS_LOSS_PRESCRIPTION, std::string> RSG_MASS_LOSS_PRESCRIPTION_LABEL = {
    { RSG_MASS_LOSS_PRESCRIPTION::ZERO,             "ZERO" },
    { RSG_MASS_LOSS_PRESCRIPTION::VINKSABHAHIT2023, "VINKSABHAHIT2023" },
    { RSG_MASS_LOSS_PRESCRIPTION::BEASOR2020,       "BEASOR2020" },
    { RSG_MASS_LOSS_PRESCRIPTION::DECIN2023,        "DECIN2023" },
    { RSG_MASS_LOSS_PRESCRIPTION::YANG2023,         "YANG2023" },
    { RSG_MASS_LOSS_PRESCRIPTION::KEE2021,          "KEE2021" },
    { RSG_MASS_LOSS_PRESCRIPTION::NJ90,             "NJ90" }
};


// Semi-major axis distributions
enum class SEMI_MAJOR_AXIS_DISTRIBUTION: int { FLATINLOG, DUQUENNOYMAYOR1991, SANA2012 };
const std::unordered_map<SEMI_MAJOR_AXIS_DISTRIBUTION, std::string> SEMI_MAJOR_AXIS_DISTRIBUTION_LABEL = {
    { SEMI_MAJOR_AXIS_DISTRIBUTION::FLATINLOG,          "FLATINLOG" },
    { SEMI_MAJOR_AXIS_DISTRIBUTION::DUQUENNOYMAYOR1991, "DUQUENNOYMAYOR1991" },
    { SEMI_MAJOR_AXIS_DISTRIBUTION::SANA2012,           "SANA2012" }
};


// Fryer 2012 supernova engines
enum class SN_ENGINE: int { RAPID, DELAYED };
const std::unordered_map<SN_ENGINE, std::string> SN_ENGINE_LABEL = {
    { SN_ENGINE::RAPID,   "RAPID" },
    { SN_ENGINE::DELAYED, "DELAYED" }
};


// Supernova events/states
//
// The values here for SN_EVENT are powers of 2 so that they can be used in a bit map
// and manipulated with bit-wise logical operators
//
// Ordinarily we might expect that an SN event could be only one of
//
//    NONE, CCSN, ECSN, PISN, PPISN, USSN, AIC, SNIA, or HeSD
//
// Note that the CCSN value here replaces the SN value in the legacy code
// The legacy code implemented these values as boolean flags, and the SN flag was always set when
// the USSN flag was set (but not the converse).  In the legacy code when the ECSN flag was set 
// the SN flag was not set.  In the legacy code the PISN and PPISN flags were used to track history
// and we only set for the "experienced" condition (I think).
//
// To match the legacy code usage of these flags, here the "is" and "experienced" conditions 
// ("current" and "past" SN events) are implemented as bit maps - different values can be
// ORed or ANDed into the bit map (that way the USSN and CCSN flags can be set at the same
// time - necessary for the code flow (from the legacy code) - which we should probably one
// day look at and rewrite).
//
// HeSD stands for helium-shell detonation
//
// A convenience function has been provided in utils.cpp to interpret the bit map (utils::SNEventType()).
// Given an SN_EVENT bitmap (current or past), it returns (in priority order):
//     
//    SN_EVENT::NONE    iff no bits are set
//    SN_EVENT::CCSN    iff CCSN  bit is set and USSN bit is not set
//    SN_EVENT::ECSN    iff ECSN  bit is set
//    SN_EVENT::PISN    iff PISN  bit is set
//    SN_EVENT::PPISN   iff PPISN bit is set
//    SN_EVENT::USSN    iff USSN  bit is set
//    SN_EVENT::AIC     iff AIC   bit is set
//    SN_EVENT::SNIA    iff SNIA  bit is set and HeSD bit is not set
//    SN_EVENT::HeSD    iff HeSD  bit is set
//    SN_EVENT::UNKNOWN otherwise
//
enum class SN_EVENT: int { 
    NONE         = 0, 
    CCSN         = 1, 
    ECSN         = 2, 
    PISN         = 4, 
    PPISN        = 8, 
    USSN         = 16,
    AIC          = 32,
    SNIA         = 64,
    HeSD         = 128,
    UNKNOWN      = 32768 // doesn't really matter what this is because the value is never used (as long as it's > sum of all the others)
};
const std::unordered_map<SN_EVENT, std::string> SN_EVENT_LABEL = {
    { SN_EVENT::NONE,    "No Supernova" },
    { SN_EVENT::CCSN,    "Core Collapse Supernova" },
    { SN_EVENT::ECSN,    "Electron Capture Supernova" },
    { SN_EVENT::PISN,    "Pair Instability Supernova" },
    { SN_EVENT::PPISN,   "Pulsational Pair Instability Supernova" },
    { SN_EVENT::USSN,    "Ultra Stripped Supernova" },
    { SN_EVENT::AIC,     "Accretion-Induced Collapse" }, 
    { SN_EVENT::SNIA,    "Supernova Type Ia" }, 
    { SN_EVENT::HeSD,    "Helium-shell detonation" }, 
    { SN_EVENT::UNKNOWN, "Unknown Supernova Type" }
};
ENABLE_BITMASK_OPERATORS(SN_EVENT);


// Supernova states
enum class SN_STATE: int { NONE, STAR1, STAR2, BOTH };
const std::unordered_map<SN_STATE, std::string> SN_STATE_LABEL = {
    { SN_STATE::NONE,  "No Supernova" },
    { SN_STATE::STAR1, "Star1 only" },
    { SN_STATE::STAR2, "Star2 only" },
    { SN_STATE::BOTH,  "Both stars" }
};


// Stellar populations
enum class STELLAR_POPULATION: int { POPULATION_I, POPULATION_II };
const std::unordered_map<STELLAR_POPULATION, std::string> STELLAR_POPULATION_LABEL = {
    { STELLAR_POPULATION::POPULATION_I,  "POPULATION_I" },
    { STELLAR_POPULATION::POPULATION_II, "POPULATION_II" }
};


// Tides prescriptions
enum class TIDES_PRESCRIPTION: int { NONE, PERFECT, KAPIL2026, ZAHN1977 };
const std::unordered_map<TIDES_PRESCRIPTION, std::string> TIDES_PRESCRIPTION_LABEL = {
    { TIDES_PRESCRIPTION::NONE,      "NONE" },
    { TIDES_PRESCRIPTION::PERFECT,   "PERFECT" },
    { TIDES_PRESCRIPTION::KAPIL2026, "KAPIL2026" },
    { TIDES_PRESCRIPTION::ZAHN1977,  "ZAHN1977" }
};


// VMS (very massive stars) mass loss prescriptions
enum class VMS_MASS_LOSS_PRESCRIPTION: int { ZERO, VINK2011, BESTENLEHNER2020, SABHAHIT2023};
const std::unordered_map<VMS_MASS_LOSS_PRESCRIPTION, std::string> VMS_MASS_LOSS_PRESCRIPTION_LABEL = {
    { VMS_MASS_LOSS_PRESCRIPTION::ZERO,             "ZERO" },
    { VMS_MASS_LOSS_PRESCRIPTION::VINK2011,         "VINK2011" },
    { VMS_MASS_LOSS_PRESCRIPTION::BESTENLEHNER2020, "BESTENLEHNER2020" },
    { VMS_MASS_LOSS_PRESCRIPTION::SABHAHIT2023,     "SABHAHIT2023" }
};

// WR mass loss prescriptions
enum class WR_MASS_LOSS_PRESCRIPTION: int { ZERO, BELCZYNSKI2010, SANDERVINK2023, SHENAR2019 };
const std::unordered_map<WR_MASS_LOSS_PRESCRIPTION, std::string> WR_MASS_LOSS_PRESCRIPTION_LABEL = {
    { WR_MASS_LOSS_PRESCRIPTION::ZERO,           "ZERO"},
    { WR_MASS_LOSS_PRESCRIPTION::BELCZYNSKI2010, "BELCZYNSKI2010" },
    { WR_MASS_LOSS_PRESCRIPTION::SANDERVINK2023, "SANDERVINK2023" },
    { WR_MASS_LOSS_PRESCRIPTION::SHENAR2019,     "SHENAR2019" }
};

// common envelope zeta prescriptions
enum class ZETA_PRESCRIPTION: int { SOBERMAN, HURLEY, ARBITRARY };
const std::unordered_map<ZETA_PRESCRIPTION, std::string> ZETA_PRESCRIPTION_LABEL = {
    { ZETA_PRESCRIPTION::SOBERMAN,  "SOBERMAN" },
    { ZETA_PRESCRIPTION::HURLEY,    "HURLEY" },
    { ZETA_PRESCRIPTION::ARBITRARY, "ARBITRARY" }
};

#endif // __typedefs_h__
