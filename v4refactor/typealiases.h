#pragma once

#include "typedefs.h"
#include <type_traits>


    namespace get {
        struct unwrappedT  {};
        struct asOptionalT {};
        inline constexpr unwrappedT  unwrapped{};
        inline constexpr asOptionalT asOptional{};
    }

// stellar type list initializer
using StellarTypeListT = std::initializer_list<STELLAR_TYPE>;


// immediate events
// "declared using unnamed type" issue with using-alias-of-anonymous-struct pattern)
struct ImmediateEventT_struct {
    IMMEDIATE_EVENT eventType;
    STELLAR_TYPE    stellarType1;
    STELLAR_TYPE    stellarType2;
};
using ImmediateEventT = ImmediateEventT_struct;


// mass loss
using MassLossT = std::tuple<double, ML_TYPE>;



// Hurley timescales
namespace HURLEY_TS {
    inline constexpr auto MS        = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::MS)>>(HURLEY_TIMESCALES::MS);           // Main sequence
    inline constexpr auto BGB       = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::BGB)>>(HURLEY_TIMESCALES::BGB);          // Base of Giant Branch
    inline constexpr auto HeI       = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::HeI)>>(HURLEY_TIMESCALES::HeI);          // Helium ignition
    inline constexpr auto He        = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::He)>>(HURLEY_TIMESCALES::He);           // Helium burning

                                                                                        // First Giant Branch (FGB)
    inline constexpr auto Inf1_FGB  = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Inf1_FGB)>>(HURLEY_TIMESCALES::Inf1_FGB);     // First Giant Branch tinf1 (integration constant)
    inline constexpr auto Inf2_FGB  = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Inf2_FGB)>>(HURLEY_TIMESCALES::Inf2_FGB);     // First Giant Branch tinf2 (integration constant)
    inline constexpr auto Mx_FGB    = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Mx_FGB)>>(HURLEY_TIMESCALES::Mx_FGB);       // First Giant Branch t(Mx)

                                                                                        // Early Asymptotic Giant Branch (EAGB) (FAGB in Hurley's sse)
    inline constexpr auto Inf1_FAGB = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Inf1_FAGB)>>(HURLEY_TIMESCALES::Inf1_FAGB);    // Early Asymptotic Giant Branch tinf1 (integration constant)
    inline constexpr auto Inf2_FAGB = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Inf2_FAGB)>>(HURLEY_TIMESCALES::Inf2_FAGB);    // Early Asymptotic Giant Branch tinf2 (integration constant)
    inline constexpr auto Mx_FAGB   = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Mx_FAGB)>>(HURLEY_TIMESCALES::Mx_FAGB);      // Early Asymptotic Giant Branch t(Mx)

                                                                                        // Thermally Pulsating Asymptotic Giant Branch (TPAGB) (SAGB in Hurley's sse)
    inline constexpr auto Inf1_SAGB = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Inf1_SAGB)>>(HURLEY_TIMESCALES::Inf1_SAGB);    // Thermally Pulsating Asymptotic Giant Branch tinf1 (integration constant)
    inline constexpr auto Inf2_SAGB = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Inf2_SAGB)>>(HURLEY_TIMESCALES::Inf2_SAGB);    // Thermally Pulsating Asymptotic Giant Branch tinf2 (integration constant)
    inline constexpr auto Mx_SAGB   = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Mx_SAGB)>>(HURLEY_TIMESCALES::Mx_SAGB);      // Thermally Pulsating Asymptotic Giant Branch t(Mx)

    inline constexpr auto DU2       = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::DU2)>>(HURLEY_TIMESCALES::DU2);          // Lifetime to second dredge-up

                                                                                        // Helium Giant Branch
    inline constexpr auto HeMS      = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::HeMS)>>(HURLEY_TIMESCALES::HeMS);         // Naked Helium Star central helium burning lifetime (HeMs)
    inline constexpr auto Inf1_HeGB = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Inf1_HeGB)>>(HURLEY_TIMESCALES::Inf1_HeGB);    // Helium Giant Branch tinf1 (integration constant)
    inline constexpr auto Inf2_HeGB = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Inf2_HeGB)>>(HURLEY_TIMESCALES::Inf2_HeGB);    // Helium Giant Branch tinf2 (integration constant)
    inline constexpr auto Mx_HeGB   = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Mx_HeGB)>>(HURLEY_TIMESCALES::Mx_HeGB);      // Helium Giant Branch tx (is this t(Mx)?)

    inline constexpr auto Tau_BL    = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::Tau_BL)>>(HURLEY_TIMESCALES::Tau_BL);       // Relative duration of blue loop taubl
    inline constexpr auto TauX_BL   = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::TauX_BL)>>(HURLEY_TIMESCALES::TauX_BL);      // Relative start of blue loop taux
    inline constexpr auto TauY_BL   = static_cast<std::underlying_type_t<decltype(HURLEY_TIMESCALES::TauY_BL)>>(HURLEY_TIMESCALES::TauY_BL);      // Relative end of blue loop tauy
}

namespace TIMESCALES = HURLEY_TS;   // For now - we can deal with new timescales when they happen.
namespace TS = HURLEY_TS;           // For now - we can deal with new timescales when they happen.

// Hurley Giantbranch parameters
namespace HURLEY_GBP {
    inline constexpr auto AH     = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::AH);               // Hydrogen rate constant.  Hurley et al. 2000, p553
    inline constexpr auto AHHe   = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::AHHe);             // Effective combined rate constant for both hydrogen and helium shell burning.  Hurley et al. 2000, eq 71
    inline constexpr auto AHe    = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::AHe);              // Helium rate constant.  Hurley et al. 2000, eq 68
    inline constexpr auto B      = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::B);                // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    inline constexpr auto D      = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::D);                // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    inline constexpr auto p      = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::p);                // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    inline constexpr auto q      = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::q);                // Hurley et al. 2000, p552, eq38 (does this represent something physical?  If so, what?  How should this be described?)
    inline constexpr auto Lx     = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::Lx);               // Luminosity parameter on the first giant branch (FGB) Lx as a function of the core mass (really a function of Mx).
    inline constexpr auto Mx     = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::Mx);               // Crosover point of high-luminosity and low-luminosity in core mass - luminosity relation. Hurley et al. 2000, p552, eq38
    inline constexpr auto McBGB  = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::McBGB);            // Core mass at BGB (Base of Giant Branch)
    inline constexpr auto McBAGB = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::McBAGB);           // Core mass at BAGB (Base of Asymptotic Giant Branch).  Hurley et al. 2000, eq 66 (also see eq 75 and discussion)
    inline constexpr auto McDU   = static_cast<std::underlying_type_t<HURLEY_GB_PARAMETERS>>(HURLEY_GB_PARAMETERS::McDU);             // Core mass at second dredge up.  Hurley et al. 2000, eq 69
}

namespace GBPARAMS = HURLEY_GBP;    // For now - we can deal with new GB parameters when they happen.
namespace GBP = HURLEY_GBP;         // For now - we can deal with new GB parameters when they happen.


// Hurley mass cutoffs
namespace HURLEY_MCO {
    inline constexpr auto Hook = static_cast<std::underlying_type_t<HURLEY_MASS_CUTOFFS>>(HURLEY_MASS_CUTOFFS::Hook);                // Mass above which hook appears on MS (in Msol)
    inline constexpr auto HeF  = static_cast<std::underlying_type_t<HURLEY_MASS_CUTOFFS>>(HURLEY_MASS_CUTOFFS::HeF);                 // Maximum initial mass for which helium ignites degenerately in a Helium Flash (HeF)
    inline constexpr auto FGB  = static_cast<std::underlying_type_t<HURLEY_MASS_CUTOFFS>>(HURLEY_MASS_CUTOFFS::FGB);                 // Maximum initial mass for which helium ignites on the First Giant Branch (FGB)
}


// Hurley R constants
namespace HURLEY_R {
    inline constexpr auto B_ALPHA_R = static_cast<std::underlying_type_t<HURLEY_R_CONSTANTS>>(HURLEY_R_CONSTANTS::B_ALPHA_R);
    inline constexpr auto C_ALPHA_R = static_cast<std::underlying_type_t<HURLEY_R_CONSTANTS>>(HURLEY_R_CONSTANTS::C_ALPHA_R);
    inline constexpr auto B_BETA_R  = static_cast<std::underlying_type_t<HURLEY_R_CONSTANTS>>(HURLEY_R_CONSTANTS::B_BETA_R);
    inline constexpr auto C_BETA_R  = static_cast<std::underlying_type_t<HURLEY_R_CONSTANTS>>(HURLEY_R_CONSTANTS::C_BETA_R);
    inline constexpr auto B_DELTA_R = static_cast<std::underlying_type_t<HURLEY_R_CONSTANTS>>(HURLEY_R_CONSTANTS::B_DELTA_R);
}


// Hurley L constants
namespace HURLEY_L {
    inline constexpr auto B_ALPHA_L = static_cast<std::underlying_type_t<HURLEY_L_CONSTANTS>>(HURLEY_L_CONSTANTS::B_ALPHA_L);
    inline constexpr auto B_BETA_L  = static_cast<std::underlying_type_t<HURLEY_L_CONSTANTS>>(HURLEY_L_CONSTANTS::B_BETA_L);
    inline constexpr auto B_DELTA_L = static_cast<std::underlying_type_t<HURLEY_L_CONSTANTS>>(HURLEY_L_CONSTANTS::B_DELTA_L);
}


// Hurley GAMMA constants
namespace HURLEY_G {
    inline constexpr auto B_GAMMA = static_cast<std::underlying_type_t<HURLEY_GAMMA_CONSTANTS>>(HURLEY_GAMMA_CONSTANTS::B_GAMMA);
    inline constexpr auto C_GAMMA = static_cast<std::underlying_type_t<HURLEY_GAMMA_CONSTANTS>>(HURLEY_GAMMA_CONSTANTS::C_GAMMA);
}


// Hurley A & B coefficients
namespace HURLEY_AB {
    inline constexpr auto ALPHA = static_cast<std::underlying_type_t<AB_TCoeff>>(AB_TCoeff::ALPHA);
    inline constexpr auto BETA  = static_cast<std::underlying_type_t<AB_TCoeff>>(AB_TCoeff::BETA);
    inline constexpr auto GAMMA = static_cast<std::underlying_type_t<AB_TCoeff>>(AB_TCoeff::GAMMA);
    inline constexpr auto ETA   = static_cast<std::underlying_type_t<AB_TCoeff>>(AB_TCoeff::ETA);
    inline constexpr auto MU    = static_cast<std::underlying_type_t<AB_TCoeff>>(AB_TCoeff::MU);
}


// Tout L coefficients
namespace TOUT_L {
    inline constexpr auto ALPHA   = static_cast<std::underlying_type_t<L_Coeff>>(L_Coeff::ALPHA);
    inline constexpr auto BETA    = static_cast<std::underlying_type_t<L_Coeff>>(L_Coeff::BETA);
    inline constexpr auto GAMMA   = static_cast<std::underlying_type_t<L_Coeff>>(L_Coeff::GAMMA);
    inline constexpr auto DELTA   = static_cast<std::underlying_type_t<L_Coeff>>(L_Coeff::DELTA);
    inline constexpr auto EPSILON = static_cast<std::underlying_type_t<L_Coeff>>(L_Coeff::EPSILON);
    inline constexpr auto ZETA    = static_cast<std::underlying_type_t<L_Coeff>>(L_Coeff::ZETA);
    inline constexpr auto ETA     = static_cast<std::underlying_type_t<L_Coeff>>(L_Coeff::ETA);
}


// Tout R coefficients
namespace TOUT_R {
    inline constexpr auto THETA   = static_cast<std::underlying_type_t<R_Coeff>>(R_Coeff::THETA);
    inline constexpr auto IOTA    = static_cast<std::underlying_type_t<R_Coeff>>(R_Coeff::IOTA);
    inline constexpr auto KAPPA   = static_cast<std::underlying_type_t<R_Coeff>>(R_Coeff::KAPPA);
    inline constexpr auto LAMBDA  = static_cast<std::underlying_type_t<R_Coeff>>(R_Coeff::LAMBDA);
    inline constexpr auto MU      = static_cast<std::underlying_type_t<R_Coeff>>(R_Coeff::MU);
    inline constexpr auto NU      = static_cast<std::underlying_type_t<R_Coeff>>(R_Coeff::NU);
    inline constexpr auto XI      = static_cast<std::underlying_type_t<R_Coeff>>(R_Coeff::XI);
    inline constexpr auto OMICRON = static_cast<std::underlying_type_t<R_Coeff>>(R_Coeff::OMICRON);
    inline constexpr auto PI      = static_cast<std::underlying_type_t<R_Coeff>>(R_Coeff::PI);
}


// Tout LR coefficients
namespace TOUT_LR {
    inline constexpr auto a = static_cast<std::underlying_type_t<LR_TCoeff>>(LR_TCoeff::a);
    inline constexpr auto b = static_cast<std::underlying_type_t<LR_TCoeff>>(LR_TCoeff::b);
    inline constexpr auto c = static_cast<std::underlying_type_t<LR_TCoeff>>(LR_TCoeff::c);
    inline constexpr auto d = static_cast<std::underlying_type_t<LR_TCoeff>>(LR_TCoeff::d);
    inline constexpr auto e = static_cast<std::underlying_type_t<LR_TCoeff>>(LR_TCoeff::e);
}















struct SupernovaEventsT {               // Supernova events - current and past
    SN_EVENT current;    // At the current timestep: SN_EVENT::NONE if no supernova event happening
    SN_EVENT past;       // At any past timestep   : SN_EVENT::NONE if no supernova event happened in any past timestep
};


struct StellarKickParmsT {          // Stellar SN kick parameters
    OptDblT magnitude;              // Kick magnitude (km s^-1)
    OptDblT magnitudeDrawn;         // Kick magnitude drawn from user-specified distribution (km s^-1) (may be reweighted - see `magnitude`)
    OptDblT magnitudeRandom;        // Random number U(0,1) for sampling the supernova kick magnitude
};

struct StellarBinaryKickParmsT {    // Stellar Binary SN kick parameters
    OptDblT eccentricAnomaly;       // Eccentric anomaly at instataneous time of the SN
    OptDblT trueAnomaly;            // True anomaly at instantaneous time of the SN
    OptDblT meanAnomaly;            // Mean anomaly at instantaneous time of the SN - uniform in [0, 2pi]
    OptDblT phi;                    // Kick angle in the orbital plane, defined CCW from the radial vector pointed away from the Companion (rad) [0, 2pi)
    OptDblT theta;                  // Kick angle out of the orbital plane, toward the orbital angular momentum axis (rad) [-pi/2, pi/2]

    OptDblT rocketKickMagnitude;    // Rocket kick magnitude the system received after the supernova (km s^-1)
    OptDblT rocketKickPhi;          // Rocket kick azimuthal angle phi the system received after the supernova
    OptDblT rocketKickTheta;        // Rocket kick polar angle theta the system received after the supernova
};





struct StellarSNDetailsT {                              // Stellar SN flags and attributes
    STELLAR_TYPE      stellarType;        // Stellar type immediately prior to CO formation
    StellarKickParmsT kickParameters;     // Kick parameters 
    SupernovaEventsT  events;             // Supernova events - current and past
    double            fallbackFraction;   // Fallback fraction during a supernova event
    bool              isHydrogenPoor;     // Flag to indicate if exploding star is hydrogen-poor. We consider an H-rich star all SN progenitors that have an H envelope, otherwise H-poor
};



struct BinarySNDetailsT {                            // Binary star SN flags and attributes
    StellarBinaryKickParmsT kickParameters;    // User-supplied initial kick parameters - if present used in place of drawing randomly/from distributions

    SN_STATE                  supernovaState;        // Indicates which star (or stars) are undergoing / have undergone a supernova event
};








// BOOST variant definition for allowed data types
// used for variable specification to define logfile records
using CompasVariableT = boost::variant<
    bool,
    ShortT,
    int,
    LongT,
    LongLongT,
    UShortT,
    UIntT,
    ULongT,
    ULongLongT,
    float,
    double,
    LongDblT,
    DblVectorT,
    StrT,
    StrVectorT,
    ERROR,
    STELLAR_TYPE,
    MT_CASE,
    MT_TRACKING,
    MT_TIMESCALE,
    SN_EVENT,
    SN_STATE,
    EVOLUTION_STATUS
>;


using _STELLAR_TYPE_VEC_    = std::vector<STELLAR_TYPE>;
using _COMPAS_VARIABLE_VEC_ = std::vector<CompasVariableT>;









using GridfileDetailsT = struct {           // Gridfile details
    StrT           filename;               // filename for grid file
    ERROR           error;                  // status - ERROR::NONE if no problem, otherwise an error number
    std::ifstream   handle;                 // the file handle
    std::streamsize startLine;              // the first line of the grid file to process (0-based)
    std::streamsize currentLine;            // the grid line currently being processed
    std::streamsize linesProcessed;         // the number of grid lines processed so far in this run
    std::streamsize linesToProcess;         // the number of grid lines to process (from start line)
};







using MassTransferEventT = struct {     // Details of a mass transfer event
    STELLAR_TYPE donorST;               // Stellar type of the donor star
    STELLAR_TYPE accretorST;            // Stellar type of the accretor star
    MT_CASE      MTcase;                // MT_CASE for this evene
};























// RLOF properties
// JR: add descriptive comments

// JR: add descriptive comments

// JR: add descriptive comments

// Common Envelope properties
// JR: add descriptive comments

// JR: add descriptive comments


// JR: add descriptive comments

// JR: add descriptive comments


// For boost ODE integrators, see https://www.boost.org/doc/libs/1_83_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/tutorial/harmonic_oscillator.html
typedef DblVectorT state_type;
typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type> error_stepper_type;
typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;

// MTEventT - record of a mass-transfer event when this star was the donor
// (used by ConstituentState::m_MTdonorHistory).
// Fields populated via aggregate-init: { donorST, accretorST, MTcase }.
struct MTEventT {
    STELLAR_TYPE donorST;
    STELLAR_TYPE accretorST;
    MT_CASE      MTcase;
};

