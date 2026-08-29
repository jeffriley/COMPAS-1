#include <iostream>
#include <stdarg.h>
#include <fstream>
#include <string>
#include <iomanip>

#include "Options.h"
#include "yaml.h"
#include "utils.h"


/*
 * YAML file creation
 *
 * Code to create a new YAML file based on current COMPAS defaults and any user-input option values.
 * See yaml.h for description of the YAML template.
 * 
 * There is no doubt that some of the code here is not the most efficient way of doing this, but realistically
 * this will be run once or twice every not very often - and it runs in a fraction of a second - so I'm not going
 * to engage in premature optimisation...
 */

namespace yaml {
    const std::vector<std::string> yamlTemplate {

            "",
            "booleanChoices:",
            "",
            "    ### LOGISTICS",
            "    --debug-to-file",
            "    --detailed-output                                               # WARNING! this creates a data heavy file",
            "    --emit-gravitational-radiation",
            "    --enable-warnings                                               # option to enable/disable warning messages",
            "    --errors-to-file",
            "    --evolve-main-sequence-mergers",
            "    --evolve-unbound-systems",
            "    --population-data-printing",
            "    --print-bool-as-string",
            "    --quiet",
            "    --rlof-printing",
            "    --store-input-files",
            "    --switch-log",
            "",
            "    ### STELLAR PROPERTIES",
            "    --check-photon-tiring-limit",
            "    --enable-rotationally-enhanced-mass-loss",
            "    --enhance-CHE-lifetimes-luminosities",
            "    --expel-convective-envelope-above-luminosity-threshold",
            "    --natal-kick-for-PPISN",
            "    --scale-CHE-mass-loss-with-surface-helium-abundance",
            "",
            "    ### BINARY PROPERTIES",
            "    --allow-touching-at-birth                                       # record binaries that have stars touching at birth in output files",
            "    --include-WD-binaries-as-DCO                                    # record WD binaries in the BSE DCO file",
            "",
            "    ### MASS TRANSFER",
            "    --angular-momentum-conservation-during-circularisation",
            "    --allow-rlof-at-birth                                           # allow binaries that have one or both stars in RLOF at birth to evolve, particularly useful in the context of CHE binaries",
            "    --circularise-binary-during-mass-transfer",
            "    --hmxr-binaries",
            "    --use-mass-transfer",
            "    --retain-core-mass-during-caseA-mass-transfer",
            "",
            "    ### COMMON ENVELOPE",
            "    --common-envelope-allow-immediate-RLOF-post-CE-survive",
            "    --common-envelope-allow-main-sequence-survive                   # Allow main sequence stars to survive CE",
            "    --common-envelope-allow-radiative-envelope-survive",
            "    --common-envelope-lambda-nanjing-enhanced",
            "    --common-envelope-lambda-nanjing-interpolate-in-mass",
            "    --common-envelope-lambda-nanjing-interpolate-in-metallicity",
            "    --common-envelope-lambda-nanjing-use-rejuvenated-mass",
            "    --revised-energy-formalism-nandez-ivanova",
            "",
            "    ### SUPERNOVAE, KICKS AND REMNANTS",
            "    --allow-non-stripped-ECSN",
            "    --pair-instability-supernovae",
            "    --pulsational-pair-instability",
            "    --USSN-kicks-override-mandel-muller",
            "",
            "    ### PULSAR PARAMETERS",
            "    --evolve-pulsars",
            "",
            "   ### WHITE DWARF PARAMETERS",
            "    --evolve-double-white-dwarfs",
            "",
            "",
            "numericalChoices:",
            "",
            "    ### LOGISTICS",
            "    --debug-level",
            "    --logfile-common-envelopes-record-types",
            "    --logfile-detailed-output-record-types",
            "    --logfile-double-compact-objects-record-types",
            "    --logfile-pulsar-evolution-record-types",
            "    --logfile-rlof-parameters-record-types",
            "    --logfile-supernovae-record-types",
            "    --logfile-system-snapshot-log-record-types",
            "    --logfile-system-parameters-record-types",
            "    --grid-lines-to-process",
            "    --grid-start-line",
            "    --hdf5-chunk-size",
            "    --hdf5-buffer-size",
            "    --log-level",
            "    --mass-change-fraction                                          # approximate desired fractional changes in stellar mass per timestep",
            "    --maximum-evolution-time                                        # maximum physical time a system can be evolved [Myr]",
            "    --maximum-number-timestep-iterations",
            "    --number-of-systems                                             # number of systems per batch",
            "    --radial-change-fraction                                        # approximate desired fractional changes in stellar radius per timestep",
            "    --system-snapshot-age-thresholds                                # age thresholds for the system snapshot log",
            "    --system-snapshot-time-thresholds                               # time thresholds for the system snapshot log",
            "    --timestep-multiplier                                           # optional multiplier relative to default time step duration",
            "    --timestep-multipliers                                          # optional phase-dependent multipliers relative to default time step duration",
            "",
            "    ### STELLAR PROPERTIES",
            "    --cool-wind-mass-loss-multiplier",
            "    --initial-mass                                                  # initial mass for SSE",
            "    --initial-mass-function-min                                     # use 5.0 for DCOs [Msol]",
            "    --initial-mass-function-max                                     # stellar tracks extrapolated above 50 Msol (Hurley+2000) [Msol]",
            "    --initial-mass-function-power",
            "    --luminosity-to-mass-threshold",
            "    --metallicity                                                   # metallicity for both SSE and BSE - Solar metallicity Asplund+2010",
            "    --metallicity-min",
            "    --metallicity-max",
            "    --luminous-blue-variable-multiplier",
            "    --overall-wind-mass-loss-multiplier",
            "    --random-seed",
            "    --rotational-frequency",
            "    --rotational-frequency-1",
            "    --rotational-frequency-2",
            "    --scale-terminal-wind-velocity-with-metallicity-power",
            "    --wolf-rayet-multiplier",
            "",
            "    ### BINARY PROPERTIES",
            "    --eccentricity                                                  # eccentricity for BSE",
            "    --eccentricity-min",
            "    --eccentricity-max",
            "    --initial-mass-1                                                # primary initial mass for BSE",
            "    --initial-mass-2                                                # secondary initial mass for BSE",
            "    --mass-ratio",
            "    --mass-ratio-min",
            "    --mass-ratio-max",
            "    --minimum-sampled-secondary-mass                                # Brown dwarf limit [Msol]",
            "    --orbital-period                                                # orbital period for BSE",
            "    --orbital-period-min                                            # [days]",
            "    --orbital-period-max                                            # [days]",
            "    --semi-major-axis                                               # semi-major axis for BSE",
            "    --semi-major-axis-min                                           # [AU]",
            "    --semi-major-axis-max                                           # [AU]",
            "",
            "    ### MASS TRANSFER",
            "    --convective-envelope-mass-threshold                            # Only if using envelope-state-prescription = 'CONVECTIVE_MASS_FRACTION'",
            "    --convective-envelope-temperature-threshold                     # Only if using envelope-state-prescription = 'FIXED_TEMPERATURE'",
            "    --critical-mass-ratio-HG-degenerate-accretor",
            "    --critical-mass-ratio-HG-non-degenerate-accretor",
            "    --critical-mass-ratio-MS-high-mass-degenerate-accretor",
            "    --critical-mass-ratio-MS-high-mass-non-degenerate-accretor",
            "    --critical-mass-ratio-MS-low-mass-degenerate-accretor",
            "    --critical-mass-ratio-MS-low-mass-non-degenerate-accretor",
            "    --critical-mass-ratio-giant-degenerate-accretor",
            "    --critical-mass-ratio-giant-non-degenerate-accretor",
            "    --critical-mass-ratio-helium-HG-degenerate-accretor",
            "    --critical-mass-ratio-helium-HG-non-degenerate-accretor",
            "    --critical-mass-ratio-helium-MS-degenerate-accretor",
            "    --critical-mass-ratio-helium-giant-degenerate-accretor",
            "    --critical-mass-ratio-helium-MS-non-degenerate-accretor",
            "    --critical-mass-ratio-helium-giant-non-degenerate-accretor",
            "    --critical-mass-ratio-white-dwarf-degenerate-accretor",
            "    --critical-mass-ratio-white-dwarf-non-degenerate-accretor",
            "    --mass-transfer-fa                                              # Only if using mass-transfer-accretion-efficiency-prescription = 'FIXED'",
            "    --mass-transfer-jloss                                           # Only if using mass-transfer-angular-momentum-loss-prescription = 'FIXED'",
            "    --mass-transfer-jloss-linear-fraction-degen",
            "    --mass-transfer-jloss-linear-fraction-non-degen",
            "    --mass-transfer-thermal-limit-C",
            "    --zeta-adiabatic-arbitrary",
            "    --zeta-main-sequence",
            "    --zeta-radiative-envelope-giant",
            "",
            "    ### COMMON ENVELOPE",
            "    --common-envelope-alpha",
            "    --common-envelope-alpha-thermal                                 # lambda = alpha_th*lambda_b + (1-alpha_th)*lambda_g",
            "    --common-envelope-lambda                                        # Only if using 'LAMBDA_FIXED'",
            "    --common-envelope-lambda-multiplier                             # Multiply common envelope lambda by some constant",
            "    --common-envelope-mass-accretion-constant",
            "    --common-envelope-mass-accretion-max                            # For 'MACLEOD+2014' [Msol]",
            "    --common-envelope-mass-accretion-min                            # For 'MACLEOD+2014' [Msol]",
            "    --common-envelope-recombination-energy-density",
            "    --common-envelope-second-stage-beta",
            "    --common-envelope-slope-kruckow",
            "    --maximum-mass-donor-nandez-ivanova",
            "",
            "    ### SUPERNOVAE, KICKS AND REMNANTS",
            "    --eddington-accretion-factor                                    # multiplication Factor for eddington accretion onto NS&BH",
            "    --fix-dimensionless-kick-magnitude",
            "    --fryer-22-fmix                                                 # parameter describing mixing growth time when using the 'FRYER2022' remnant mass prescription",
            "    --fryer-22-mcrit                                                # critical mass for BH formation when using the 'FRYER2022' remnant mass prescription",
            "    --kick-direction-power",
            "    --kick-magnitude-sigma-CCSN-NS                                  # [km/s]",
            "    --kick-magnitude-sigma-CCSN-BH                                  # [km/s]",
            "    --kick-magnitude-max",
            "    --kick-magnitude-random                                         # (SSE) used to draw the kick magnitude for the star should it undergo a supernova event",
            "    --kick-magnitude                                                # (SSE) (drawn) kick magnitude for the star should it undergo a supernova event [km/s]",
            "    --kick-magnitude-random-1                                       # (BSE) used to draw the kick magnitude for the primary star should it undergo a supernova event",
            "    --kick-magnitude-1                                              # (BSE) (drawn) kick magnitude for the primary star should it undergo a supernova event [km/s]",
            "    --kick-theta-1                                                  # (BSE) angle between the orbital plane and the 'z' axis of the supernova vector for the primary star should it undergo a supernova event [radians]",
            "    --kick-phi-1                                                    # (BSE) angle between 'x' and 'y', both in the orbital plane of the supernova vector, for the primary star should it undergo a supernova event [radians]",
            "    --kick-mean-anomaly-1                                           # (BSE) mean anomaly at the instant of the supernova for the primary star should it undergo a supernova event - should be uniform in [0, 2pi) [radians]",
            "    --kick-magnitude-random-2                                       # (BSE) used to draw the kick velocity for the secondary star should it undergo a supernova event",
            "    --kick-magnitude-2                                              # (BSE) (drawn) kick magnitude for the secondary star should it undergo a supernova event [km/s]",
            "    --kick-theta-2                                                  # (BSE) angle between the orbital plane and the 'z' axis of the supernova vector for the secondary star should it undergo a supernova event [radians]",
            "    --kick-phi-2                                                    # (BSE) angle between 'x' and 'y', both in the orbital plane of the supernova vector, for the secondary star should it undergo a supernova event [radians]",
            "    --kick-mean-anomaly-2                                           # (BSE) mean anomaly at the instant of the supernova for the secondary star should it undergo a supernova event - should be uniform in [0, 2pi) [radians]",
            "    --kick-magnitude-sigma-ECSN                                     # [km/s]",
            "    --kick-magnitude-sigma-USSN                                     # [km/s]",
            "    --kick-scaling-factor",
            "    --maximum-neutron-star-mass",
            "    --mcbur1",
            "    --muller-mandel-kick-multiplier-BH                              # scaling prefactor for BH kicks when using the 'MULLERMANDEL' kick magnitude distribution",
            "    --muller-mandel-kick-multiplier-NS                              # scaling prefactor for NS kicks when using the 'MULLERMANDEL' kick magnitude distribution",
            "    --muller-mandel-sigma-kick-BH                                   # BH kick scatter when using the 'MULLERMANDEL' kick magnitude distribution",
            "    --muller-mandel-sigma-kick-NS                                   # NS kick scatter when using the 'MULLERMANDEL' kick magnitude distribution",
            "    --neutrino-mass-loss-BH-formation-value",
            "    --pisn-lower-limit                                              # Minimum core mass for PISN [Msol]",
            "    --pisn-upper-limit                                              # Maximum core mass for PISN [Msol]",
            "    --ppi-lower-limit                                               # Minimum core mass for PPI [Msol]",
            "    --ppi-upper-limit                                               # Maximum core mass for PPI [Msol]",
            "    --ppi-co-core-shift-hendriks                                    # CO Core shift for PPI [Msol] from Hendriks+23",
            "    --rocket-kick-magnitude-1                                       # (BSE) rocket kick magnitude for the primary should it undergo a supernova event [km/s]",
            "    --rocket-kick-magnitude-2                                       # (BSE) rocket kick magnitude for the secondary should it undergo a supernova event [km/s]",
            "    --rocket-kick-theta-1                                           # (BSE) angle between the orbital plane and the 'z' axis of the rocket kick for the primary star [radians]",
            "    --rocket-kick-phi-1                                             # (BSE) angle between 'x' and 'y', both in the orbital plane of the rocket kick , for the primary star [radians]",
            "    --rocket-kick-theta-2                                           # (BSE) angle between the orbital plane and the 'z' axis of the rocket kick for the secondary star [radians]",
            "    --rocket-kick-phi-2                                             # (BSE) angle between 'x' and 'y', both in the orbital plane of the rocket kick , for the secondary star [radians]",
            "",
            "    ### PULSAR PARAMETERS",
            "    --pulsar-birth-magnetic-field-distribution-min                  # [log10(B/G)]",
            "    --pulsar-birth-magnetic-field-distribution-max                  # [log10(B/G)]",
            "    --pulsar-birth-magnetic-field-distribution-mean                 # [log10(B/G)]",
            "    --pulsar-birth-magnetic-field-distribution-sigma                # [log10(B/G)]",
            "    --pulsar-birth-spin-period-distribution-min                     # [ms]",
            "    --pulsar-birth-spin-period-distribution-max                     # [ms]",
            "    --pulsar-birth-spin-period-distribution-mean                    # [ms]",
            "    --pulsar-birth-spin-period-distribution-sigma                   # [ms]",
            "    --pulsar-magnetic-field-decay-timescale                         # [Myr]",
            "    --pulsar-magnetic-field-decay-massscale                         # [Msol]",
            "    --pulsar-minimum-magnetic-field                                 # [log10(B/G)]",
            "",
            "",
            "stringChoices:",
            "",
            "    ### LOGISTICS",
            "    --add-options-to-sysparms",
            "    --fp-error-mode                                                 # specifies floating-point error handling mode",
            "    --grid                                                          # grid file name (e.g. 'mygrid.txt')",
            "    --mode                                                          # evolving single (SSE) or binary stars (BSE)",
            "    --notes",
            "    --notes-hdrs",
            "    --output-container",
            "    --timesteps-filename",
            "",
            "    ### STELLAR PROPERTIES",
            "    --chemically-homogeneous-evolution-mode                         # chemically homogeneous evolution mode",
            "    --envelope-state-prescription",
            "    --initial-mass-function",
            "    --LBV-mass-loss-prescription",
            "    --main-sequence-core-mass-prescription",
            "    --mass-loss-prescription",
            "    --OB-mass-loss-prescription",
            "    --RSG-mass-loss-prescription",
            "    --VMS-mass-loss-prescription",
            "    --WR-mass-loss-prescription",
            "    --metallicity-distribution",
            "    --pulsational-pair-instability-prescription",
            "",
            "    ### BINARY PROPERTIES",
            "    --eccentricity-distribution",
            "    --mass-ratio-distribution",
            "    --orbital-period-distribution",
            "    --rotational-velocity-distribution",
            "    --semi-major-axis-distribution",
            "",
            "    ### MASS TRANSFER",
            "    --case-BB-stability-prescription",
            "    --critical-mass-ratio-prescription",
            "    --stellar-zeta-prescription",
            "    --mass-transfer-angular-momentum-loss-prescription",
            "    --mass-transfer-accretion-efficiency-prescription",
            "    --mass-transfer-rejuvenation-prescription",
            "    --mass-transfer-thermal-limit-accretor-multiplier",
            "    --response-to-spin-up",
            "",
            "    ### COMMON ENVELOPE",
            "    --common-envelope-formalism",
            "    --common-envelope-lambda-prescription                           # Xu & Li 2010",
            "    --common-envelope-mass-accretion-prescription",
            "    --common-envelope-second-stage-gamma-prescription",
            "",
            "    ### TIDES",
            "    --tides-prescription",
            "",
            "    ### SUPERNOVAE, KICKS AND REMNANTS",
            "    --black-hole-kicks-mode",
            "    --fryer-supernova-engine",
            "    --kick-magnitude-distribution",
            "    --kick-direction-distribution",
            "    --maltsev-fallback",
            "    --maltsev-mode",
            "    --neutron-star-accretion-in-ce",
            "    --neutron-star-equation-of-state",
            "    --neutrino-mass-loss-BH-formation",
            "    --pulsar-birth-magnetic-field-distribution",
            "    --pulsar-birth-spin-period-distribution",
            "    --remnant-mass-prescription",
            "",
            "    ### LOGFILES AND OUTPUTS ",
            "    --logfile-type",
            "    --logfile-name-prefix",
            "    --logfile-definitions",
            "    --logfile-common-envelopes",
            "    --logfile-detailed-output",
            "    --logfile-double-compact-objects",
            "    --logfile-pulsar-evolution",
            "    --logfile-rlof-parameters",
            "    --logfile-supernovae",
            "    --logfile-switch-log",
            "    --logfile-system-snapshot-log",
            "    --logfile-system-parameters",
            "    --output-path",
            "",
            "",
            "listChoices: ",
            "",
            "    --log-classes",
            "    --debug-classes",
        };


    // ----- D7: named constants (centralised to ensure the writer and the parser agree) ----------

    // Marker prefix for COMPAS-generated header records in the YAML file. Lines
    // beginning with this prefix are regenerated on each MakeYAMLfile call
    // (i.e. not preserved from the template).
    constexpr const char* COMPAS_HEADER_MARKER  = "##~!!~##";

    // YAML comment field labels. These are written by COMPAS and not preserved
    // from any template - the template's value strings beginning with these
    // labels (up to the next '#' or end of line) are stripped.
    constexpr const char* DEFAULT_LABEL         = "# Default: ";
    constexpr const char* OPTIONS_LABEL         = "# Options: ";

    // Top-level category headers in the YAML file structure.
    constexpr const char* CATEGORY_BOOLEAN      = "booleanChoices";
    constexpr const char* CATEGORY_NUMERICAL    = "numericalChoices";
    constexpr const char* CATEGORY_STRING       = "stringChoices";
    constexpr const char* CATEGORY_LIST         = "listChoices";

    // Option-line prefixes used when writing the YAML file.
    // (the 5-space indent makes the option lines align under the category headers)
    constexpr const char* OPT_LINE_USER_SET     = "     --";   // user-supplied option (uncommented)
    constexpr const char* OPT_LINE_DEFAULT      = "#    --";   // option using its COMPAS default (commented)

    // ----- D2/D3/D4: YamlRecord struct replaces the six parallel vectors and SET_STRINGS macro ----------
    //
    // A single YamlRecord captures everything needed to emit one line of the
    // YAML file: the preamble (header text or 'set'/'unset' marker), the
    // option name, the option's current value, allowed values, default value,
    // and any trailing comment from the template. `isOption` distinguishes
    // option records (which get aligned default/allowed/comment columns) from
    // header/blank records (which are emitted as-is).
    //
    // This replaces the 3-char string prefix tags ("HDR"/"OPT") that used
    // to be embedded in the output buffer (D4).
    struct YamlRecord {
        std::string preamble;
        std::string option;
        std::string value;
        std::string allowed;
        std::string defaultValue;
        std::string comment;
        bool        isOption = false;
    };


    /*
     * Read content from a YAML template file.
     *
     * Returns a pair (status, content):
     *   - status indicates the outcome (see READ_TEMPLATE_RESULT in yaml.h)
     *   - content is the file's contents (one std::string per line); empty
     *     unless status == OK
     *
     * The caller decides what to do with the content - typically use it to
     * replace the default template for one MakeYAMLfile call. ReadYAMLtemplate
     * itself has no side effects on any global state (B4 fix: previously this
     * function mutated a global default-template vector).
     */
    std::pair<READ_TEMPLATE_RESULT, std::vector<std::string>>
    ReadYAMLtemplate(const std::string p_YAMLtemplateName) {
        std::vector<std::string> content;

        if (!utils::FileExists(p_YAMLtemplateName)) {                                                           // template file exists?
            return {READ_TEMPLATE_RESULT::FILE_NOT_FOUND, content};                                             // no
        }

        std::ifstream yamlTemplateFile(p_YAMLtemplateName);                                                     // open the file
        if (!yamlTemplateFile.is_open()) {                                                                      // open ok?
            return {READ_TEMPLATE_RESULT::IO_ERROR, content};                                                   // no - IO error
        }

        std::string rec;
        while (std::getline(yamlTemplateFile, rec)) {                                                           // for all records in file
            content.push_back(rec);                                                                             // add record to content vector
        }

        if (content.empty()) {                                                                                  // file has any content?
            return {READ_TEMPLATE_RESULT::EMPTY, content};                                                      // no
        }
        return {READ_TEMPLATE_RESULT::OK, content};                                                             // yes - return content
    }


    /*
     * Write YAML content to a YAML file.
     *
     * If the file already exists, prompt the user for permission to overwrite
     * (single Y/N response on stdin). Returns a WRITE_FILE_RESULT indicating
     * the outcome.
     */
    WRITE_FILE_RESULT WriteYAMLfile(const std::string p_YAMLname, const std::vector<std::string> p_YAMLcontent) {

        if (utils::FileExists(p_YAMLname)) {                                                                    // file exists?
                                                                                                                // yes - ask user whether to overwrite
            const std::string prompt = "YAML file '" + p_YAMLname + "' already exists - overwrite (Y/N)?";
            std::string response;
            do {
                std::cout << prompt << std::flush;
                std::string line;
                if (!std::getline(std::cin, line)) {                                                            // EOF or stream error
                    return WRITE_FILE_RESULT::EXISTS_NO_OVERWRITE;
                }
                line = utils::trim(line);                                                                       // strip surrounding whitespace
                if (line.length() == 1) response = utils::ToLower(line);                                        // accept single-char Y/N answers
            } while (response != "y" && response != "n");
            if (response == "n") return WRITE_FILE_RESULT::EXISTS_NO_OVERWRITE;
        }

        std::ofstream yamlFile(p_YAMLname);                                                                     // open the file
        if (!yamlFile.is_open()) return WRITE_FILE_RESULT::IO_ERROR;                                            // open failed

        for (const std::string& rec : p_YAMLcontent) {                                                          // for each record
            if (!yamlFile.write(rec.c_str(), rec.length())) {
                return WRITE_FILE_RESULT::IO_ERROR;                                                             // mid-write IO error
            }
        }
        return WRITE_FILE_RESULT::OK;
    }


    /*
     * Format an option value or default string for output to the YAML file.
     *
     * - floats/doubles: fixed precision (6 dp if |v|<100, else 2 dp)
     * - bools: "True" / "False" (Python-compatible)
     * - strings: add surrounding '' quotes if not already present
     * - vector-of-strings: parse on whitespace, emit as {'a', 'b', 'c'}
     * - other arithmetic types: pass through unchanged
     *
     * Pure function on its arguments; no captures.
     */
    static std::string FormatString(const std::string& p_Str,                                                               // string to be formatted
                                    const TYPENAME     p_ShortType,                                                         // short data type
                                    const std::string& p_DetailedType) {                                                    // detailed data type
        std::stringstream ss;
        switch (p_ShortType) {
            case TYPENAME::FLOAT     : {       float v = std::stof(p_Str);          ss << std::fixed << std::setprecision(v < 100 ? 6 : 2) << v;     } break;
            case TYPENAME::DOUBLE    : {      double v = std::stod(p_Str);          ss << std::fixed << std::setprecision(v < 100 ? 6 : 2) << v;     } break;
            case TYPENAME::CDOUBLE   : {      double v = std::stod(p_Str);          ss << std::fixed << std::setprecision(v < 100 ? 6 : 2) << v;     } break;   // mirrors DOUBLE; defensive
            case TYPENAME::LONGDOUBLE: { long double v = std::stold(p_Str);         ss << std::fixed << std::setprecision(v < 100 ? 6 : 2) << v;     } break;
            case TYPENAME::BOOL      : { std::string s = p_Str; s = utils::trim(s); ss << (s == "1" || utils::Equals(s, "true") ? "True" : "False"); } break;
            case TYPENAME::STRING    : {
                bool vecType = (p_DetailedType == "VECTOR<STRING>");                                                        // vector of strings?
                size_t len = p_Str.length();
                if (len == 0) {                                                                                             // empty str?
                    ss << (vecType ? "{ }" : "''");                                                                       // yes - just quotes/braces
                }
                else {                                                                                                      // no - str not empty
                    // add quotes (or braces) if not already present
                    if (!vecType) {                                                                                         // scalar string
                        if (p_Str[0] != '\'' || p_Str[len - 1] != '\'') { ss << "'"; ss << p_Str;  ss << "'"; }
                        else ss << p_Str;
                    }
                    else {                                                                                                  // vector of strings
                        std::string str = p_Str;
                        (void)utils::trim(str);                                                                             // call for the in-place side effect; cast-to-void suppresses [-Wunused-value]
                        if (str[0] != '{' || str[len - 1] != '}') {                                                         // missing braces?
                            ss << "{ ";
                            size_t start = str[0] == '{' ? 1 : 0;
                            size_t end   = 0;
                            bool   first = true;
                            while ((end = p_Str.find(' ', start)) != std::string::npos) {
                                std::string v = p_Str.substr(start, end - start);
                                if (v.length() > 0) {
                                    if (!first) ss << ", ";
                                    if (v[0] != '\'' || v[len - 1] != '\'') { ss << "'"; ss << v; ss << "'"; }
                                    else ss << v;
                                }
                                start = end + 1;
                                first = false;
                            }
                            std::string v = p_Str.substr(start, end - start);
                            if (v.length() > 0) {
                                if (!first) ss << ", ";
                                if (v[0] != '\'' || v[len - 1] != '\'') { ss << "'"; ss << v; ss << "'"; }
                                else ss << v;
                            }
                            if (str[0] != '{' || str[len - 1] != '}') ss << " }";
                        }
                        else ss << p_Str;                                                                                   // has braces - just pass through
                    }
                }
            } break;
            default: ss << p_Str;                                                                                           // just pass it through
        }
        return ss.str();
    }

    /*
     * Construct YAML file content and write to file passed in p_YAMLfilename
     *
     * Constructs YAML file content based on YAML template read from file passed in
     * p_YAMLtemplateName.  If no filename passed in p_YAMLtemplateName, the default
     * template from yaml.h is used.
     * 
     * Constructed content is written to file passed in p_YAMLfilename.
     * 
     * 
     * void MakeYAMLfile(const std::string p_YAMLfilename, const std::string p_YAMLtemplateName)
     * 
     * @param   [IN]    p_YAMLfilename          Filename to be written - should be fully qualified
     * @param   [IN]    p_YAMLtemplateName      Template filename to be read - should be fully qualified
     */

    void MakeYAMLfile(const std::string p_YAMLfilename, const std::string p_YAMLtemplateName) {
        // The single vector of YamlRecords replaces the previous SET_STRINGS macro
        // and six parallel string vectors. A local lambda `addRecord` keeps the
        // call sites short.
        std::vector<YamlRecord> records;
        auto addHeader = [&records](const std::string& p_Preamble) {
            records.push_back({p_Preamble, "", "", "", "", "", /*isOption=*/false});
        };
        auto addOption = [&records](const std::string& p_Preamble,
                                    const std::string& p_Option,
                                    const std::string& p_Value,
                                    const std::string& p_Allowed,
                                    const std::string& p_Default,
                                    const std::string& p_Comment) {
            records.push_back({p_Preamble, p_Option, p_Value, p_Allowed, p_Default, p_Comment, /*isOption=*/true});
        };


        // process the option
        // - format strings as required and populate the yaml vectors (declared above)
        auto ProcessOption = [&] (const std::string              p_OptionStr,                                               // option name string
                                  const std::string              p_ValueStr,                                                // option value string
                                  const std::vector<std::string> p_AllowedStr,                                              // option allowed value strings
                                  const std::string              p_DefaultStr,                                              // option default string
                                  const std::string              p_CommentStr,                                              // comment string
                                  const bool                     p_UserSpecified,                                           // user specified option value?
                                  const TYPENAME                 p_DataType,                                                // short data type
                                  const std::string              p_TypeStr) {                                               // detailed data type


            std::string valueStr = FormatString(p_ValueStr, p_DataType, p_TypeStr);
            std::string defaultStr = FormatString(p_DefaultStr, p_DataType, p_TypeStr);

            const std::string preamble = p_UserSpecified ? OPT_LINE_USER_SET : OPT_LINE_DEFAULT;                            // preamble - all option records are commented, except options that are user supplied

            if (p_AllowedStr.size() == 0) {                                                                                 // multiple option values?
                addOption(preamble, p_OptionStr, valueStr, "", defaultStr, p_CommentStr);                                 // no - record for YAML file output - value set
            }
            else {                                                                                                          // yes - multiple option values
                // We list the allowed values in the order they are stored in the COMPASUnorderedMap
                // in constants.h.  We could instead list in alphabetical order, but the assumption
                // is that the order in constants.h was deliberate, so we'll maintain it - if that's
                // not true then we could just sort alphabetically here
                std::string allowedValuesStr = "Options: [";                                                                // (note: this is the actual stored value; DEFAULT_LABEL/OPTIONS_LABEL are the prefixes used during parsing)
                for (size_t idx = 0; idx < p_AllowedStr.size(); idx++) {                                                    // for each allowed value
                    allowedValuesStr += p_AllowedStr[idx];                                                                  // show allowed value
                    if (idx < p_AllowedStr.size() - 1) allowedValuesStr += ",";                                             // add delimiter if necessary
                }
                allowedValuesStr += "]";                                                                                    // close bracket
                addOption(preamble, p_OptionStr, valueStr, allowedValuesStr, defaultStr, p_CommentStr);                   // option record for YAML file output
            }
        };


        // begin MakeYAMLfile() body

        // get COMPAS option details
        std::vector<OptionDetailsT> optionDetails = OPTIONS->CmdLineOptionsDetails();

        if (optionDetails.empty()) {                                                                                        // have COMPAS option details?
            std::cerr << "Unable to access COMPAS option details - no YAML file created.\n";                                // no - announce error
        }
        else {                                                                                                              // yes - have options
                                                                                                                            // process options/create YAML file
            auto wallTime           = std::chrono::system_clock::now();                                                     // get wall time
            std::time_t timeNow     = std::chrono::system_clock::to_time_t(wallTime);                                       // current time and date ...
            std::string currentTime = std::string(std::ctime(&timeNow));                                                    // ... as a string ...
            currentTime             = utils::trim(currentTime);                                                             // ... trimmed of whitespace (including newline etc.)

            // read YAML template file if supplied

            // The 'active template' for this MakeYAMLfile call. Either populated
            // from the file passed in p_YAMLtemplateName (if successfully read)
            // or, otherwise, the default in-header `activeTemplate`. The default
            // is no longer mutated (B4 fix).
            std::vector<std::string> activeTemplate = yamlTemplate;
            if (p_YAMLtemplateName.length() > 0) {                                                                          // have YAML template filename?
                auto [readResult, content] = ReadYAMLtemplate(p_YAMLtemplateName);                                          // yes - read YAML template file
                const std::string preamble = "*WARNING* File '" + p_YAMLtemplateName + "' ";
                std::string warn;
                switch (readResult) {
                    case READ_TEMPLATE_RESULT::IO_ERROR:        warn = "not read: IO error.  COMPAS default template will be used.";        break;
                    case READ_TEMPLATE_RESULT::FILE_NOT_FOUND:  warn = "not read: file does not exist.  COMPAS default template will be used."; break;
                    case READ_TEMPLATE_RESULT::EMPTY:           warn = "is empty.  COMPAS default template will be used."; break;
                    case READ_TEMPLATE_RESULT::OK:              activeTemplate = std::move(content);                       break;
                }
                if (!warn.empty()) std::cout << preamble << warn << "\n";
            }

            // add default records (file headers)

            addHeader(std::string(COMPAS_HEADER_MARKER) + " COMPAS option values");                                               // header record for YAML file output

            std::string s = std::string(COMPAS_HEADER_MARKER) + " File Created " + currentTime + " by COMPAS v" + VERSION_STRING;                       // time and version stamp
            if (s[30] == ' ') s[30] = '0';                                                                                  // add leading zero for day if necessary (ctime format 'Www Mmm dd ...', day digit at offset 9 of ctime, +22 prefix length = 30)
            addHeader(s);                                                                             // time and version stamp record for YAML file output

            addHeader(COMPAS_HEADER_MARKER);                                                                          // Blank header line
            addHeader(std::string(COMPAS_HEADER_MARKER) + " The default COMPAS YAML file (``compasConfigDefault.yaml``), as distributed, has");   // Notice regarding commented line in default YAML                                                                            // time and version stamp record for YAML file output
            addHeader(std::string(COMPAS_HEADER_MARKER) + " all COMPAS option entries commented so that the COMPAS default value for the");       // Notice regarding commented line in default YAML                                                                            // time and version stamp record for YAML file output
            addHeader(std::string(COMPAS_HEADER_MARKER) + " option is used by default. To use a value other than the COMPAS default value,");     // Notice regarding commented line in default YAML                                                                            // time and version stamp record for YAML file output
            addHeader(std::string(COMPAS_HEADER_MARKER) + " users must uncomment the entry and change the option value to the desired value.");   // Notice regarding commented line in default YAML                                                                            // time and version stamp record for YAML file output

            size_t numTemplateRecords = activeTemplate.size();                                                                // number of records in the YAML template

            // create YAML content

            if (numTemplateRecords == 0) {                                                                                  // have YAML template?
                                                                                                                            // no...
                // we have no template - this should never happen (we've got a default template in the code,
                // but we need to deal with that being missing just in case) - so we'll just write the COMPAS
                // options alphabetically within the datatype categories (booleanChoices, numericalChoices, 
                // stringChoices, listChoices)

                std::cerr << "*WARNING* No YAML template - will write COMPAS options alphabetically.\n";                    // announce warning

                // this is not pretty - we're going to iterate over the COMPAS options 4 times - once
                // for each of the datatype categories mentioned above (see note about premature optimisation
                // at the top of this file)

                const std::vector<std::string> categories { CATEGORY_BOOLEAN, CATEGORY_NUMERICAL, CATEGORY_STRING, CATEGORY_LIST }; // category names for headers

                for (size_t idx = 0; idx < categories.size(); idx++) {                                                      // for each category

                    std::string category = categories[idx];                                                                 // which category?
                    addHeader("");                                                                    // blank record for YAML file output
                    addHeader(category + ":");                                                        // category header for YAML file output

                    for (size_t optionIdx = 0; optionIdx < optionDetails.size(); optionIdx++) {                             // for each COMPAS option

                        if (optionDetails[optionIdx].dataType == TYPENAME::BOOL) {                                          // boolean option?
                            if (category != CATEGORY_BOOLEAN) continue;                                                     // skip if this category is not booleanChoices
                        }
                        else if (optionDetails[optionIdx].dataType == TYPENAME::STRING) {                                   // string, or vector of strings, option?
                            if (optionDetails[optionIdx].typeStr == "VECTOR<STRING>") {                                     // ... vector of strings
                                if (category != CATEGORY_LIST) continue;                                                    // skip if this category is not listChoices
                            }
                            else {                                                                                          // ... string
                                if (category != CATEGORY_STRING) continue;                                                  // skip if this category is not stringChoices
                            }
                        }
                        else {                                                                                              // assume numerical option
                            if (category != CATEGORY_NUMERICAL) continue;                                                  // skip if this category is not numericalChoices
                        }

                        // process the option
                        ProcessOption(optionDetails[optionIdx].optionStr,                                                   // option name string
                                      optionDetails[optionIdx].valueStr,                                                    // option value string
                                      optionDetails[optionIdx].allowedStr,                                                  // option allowed value strings
                                      optionDetails[optionIdx].defaultStr,                                                  // option default string
                                      "",                                                                                   // no comment
                                      utils::Equals(optionDetails[optionIdx].sourceStr, "user_supplied"),                   // user specified option value?
                                      optionDetails[optionIdx].dataType,                                                    // option (short) data type  
                                      optionDetails[optionIdx].typeStr);                                                    // option (detailed) data type  
                    }
                }
            }
            else {                                                                                                          // yes - have YAML template
                                                                                                                            // use template to format YAML file
                for (size_t rec = 0; rec < numTemplateRecords; rec++) {                                                     // for each template record

                    std::string thisRec = activeTemplate[rec];                                                                // template record ...
                    thisRec = utils::trim(thisRec);                                                                         // ... trimmed of whitespace

                    // need to allow for commented options in template, especially since the template
                    // might be an existing YAML file.  Since the template may have the hash character
                    // ("#") throughout, I'll refer to the hash at the start of the record (or as the
                    // first non-whitespace character) as the "leading hash" to differentiate it from
                    // other instances of the hash character in a template record
                    bool   optRec = false;                                                                                  // option template record?
                    size_t startPos = 0;                                                                                    // where to start looking for non-leading hash

                    if (thisRec.length() > 3) {                                                                             // record long enough to be option record?
                        if (thisRec.substr(0, 2) == "--") {                                                                 // yes - has leading "--"?
                            optRec = true;                                                                                  // yes - assume option record
                        }
                        else if (thisRec[0] == '#') {                                                                       // no - not leading "--": leading hash?
                            optRec = false;                                                                                 // for now...
                            if (thisRec.length() > 4) {                                                                     // enough characters to be option record?
                                std::string s = thisRec.substr(1);                                                          // yes... strip leading hash
                                s = utils::ltrim(s);                                                                        // trim whitespace from start
                                if (s.length() > 2 && s.substr(0, 2) == "--") {                                             // more than 2 characters, and option indicater?
                                    optRec = true;                                                                          // yes - looks like it might be an option record
                                    startPos = 1;                                                                           // skip over leading hash in later searches
                                }
                            }
                        }
                    }

                    // parse template record
                    if (!optRec) {                                                                                          // option record?
                        if (activeTemplate[rec].length() < 8 || activeTemplate[rec].substr(0, 8) != COMPAS_HEADER_MARKER) {     // no - skip COMPAS generated headers, and ...
                            addHeader(activeTemplate[rec]);                                             // ... output the record as is
                        }
                    }
                    else {                                                                                                  // yes - option record

                        // look for a comment on the option record and, if present, preserve it.
                        // I refer to the hash indication comment on the option record as the
                        // "comment hash" to differentiate it from all the other hashes that might
                        // be there...  Because the template might be an existing YAML file, one or
                        // both of the (special) strings "# Default:" and "# Options:" may be present,
                        // and if so we just skip over them looking for any comment on the option record

                        size_t hashStart = startPos;                                                                        // where to start looking for the comment hash
                        size_t p;                                                                                           // position in string
                        if ((p = activeTemplate[rec].find(DEFAULT_LABEL, hashStart)) != std::string::npos) {                   // "# Default: " present?
                            if (p > hashStart) hashStart = p + 1;                                                           // yes - adjust starting position for comment hash search
                        }

                        if ((p = activeTemplate[rec].find(OPTIONS_LABEL, hashStart)) != std::string::npos) {                   // "# Options: " present?
                            if (p > hashStart) hashStart = p + 1;                                                           // yes - adjust starting position for comment hash search
                        }

                        std::string optionStr;                                                                              // string to hold option name and value
                        std::string commentStr;                                                                             // string to hold optional comment
                        size_t hashLoc = activeTemplate[rec].find_first_of("#", hashStart);                                   // find the comment if there is one
                        if (hashLoc == std::string::npos) {                                                                 // found hash?
                            optionStr = activeTemplate[rec].substr(startPos);                                                 // no - no comment present - option string is the record (from startPos)
                            commentStr = "";                                                                                // no comment
                        }
                        else {                                                                                              // yes - hash (comment present)
                            optionStr = activeTemplate[rec].substr(startPos, hashLoc);                                        // option string is the record up to hash
                            commentStr = activeTemplate[rec].substr(hashLoc);                                                 // comment string (includes "#")
                        }

                        optionStr = utils::trim(optionStr);                                                                 // trim whitespace from both ends of option string

                        // need to allow for option values in template, especially since the template might
                        // be an existing YAML file. Values in the template are not preserved - we only use
                        // the template to preserve header records, blank records, option groupings, and comments.
                        size_t colonLoc = optionStr.find_first_of(":");                                                     // find the first colon if there is one
                        if (colonLoc != std::string::npos) {                                                                // found colon?
                            optionStr = optionStr.substr(0, colonLoc);                                                      // option string is everything prior to the ":""
                        }

                        // Find the value for the option (COMPAS value, not any value present in the template)
                        // - note: could be user supplied
                        size_t pos = optionStr.find_first_not_of("\t -");                                                   // find option indicator
                        std::string templateOptionStr = optionStr.substr(pos);                                              // strip "--" from template option string

                        bool found = false;                                                                                 // matching COMPAS option string found?  Not initially...
                        for (std::size_t optionIdx = 0; optionIdx < optionDetails.size(); optionIdx++) {                    // for each COMPAS option

                            std::string compasOptionStr = optionDetails[optionIdx].optionStr;                               // COMPAS option name string
                            if (utils::Equals(optionDetails[optionIdx].sourceStr, "calculated")) continue;                  // ignore calculated options                                         
                            if (utils::Equals(compasOptionStr, templateOptionStr)) {                                        // match?

                                found = true;                                                                               // yes
                                ProcessOption(templateOptionStr,                                                            // option name string (as written in the YAML template)
                                              optionDetails[optionIdx].valueStr,                                            // option (COMPAS) value string
                                              optionDetails[optionIdx].allowedStr,                                          // option allowed value strings
                                              optionDetails[optionIdx].defaultStr,                                          // option default string
                                              commentStr,                                                                   // comment string (as written in the YAML template)
                                              utils::Equals(optionDetails[optionIdx].sourceStr, "user_supplied"),           // user specified option value?
                                              optionDetails[optionIdx].dataType,                                            // option (short) data type  
                                              optionDetails[optionIdx].typeStr);                                            // option (detailed) data type
                                break;                                                                                      // we're done
                            }
                        }

                        if (!found) {                                                                                       // matching COMPAS option string found?
                            std::cerr << "*WARNING* Option '" 
                                      << templateOptionStr 
                                      << "' in YAML template is not a valid COMPAS option - ignored.\n";                    // no - announce warning
                        }
                    }
                }
            }

            // add any COMPAS options not in the YAML template - alphabetically at end

            bool extra = false;                                                                                             // extra records?  Initially false
            for (std::size_t optionIdx = 0; optionIdx < optionDetails.size(); optionIdx++) {                                // for each COMPAS option

                std::string compasOptionStr = optionDetails[optionIdx].optionStr;                                           // option name string

                if (utils::Equals(optionDetails[optionIdx].sourceStr, "calculated")) continue;                              // ignore calculated options    
                if (utils::Equals(compasOptionStr, "help")) continue;                                                       // ignore 'help' option          
                if (utils::Equals(compasOptionStr, "version")) continue;                                                    // ignore 'version' option          
                if (utils::Equals(compasOptionStr, "create-YAML-file")) continue;                                           // ignore 'create-yaml-file' option          
                if (utils::Equals(compasOptionStr, "YAML-template")) continue;                                              // ignore 'YAML-template' option          

                bool match = false;
                for (const YamlRecord& yr : records) {                                                                      // for each existing record
                    if (yr.isOption && utils::Equals(compasOptionStr, yr.option)) {                                         // match on option name?
                        match = true;
                        break;
                    }
                }

                if (!match) {                                                                                               // COMPAS option in YAML options?
                    if (!extra) {                                                                                           // no - first extra record?
                        extra = true;                                                                                       // yes - flag it
                        addHeader("\n\n### Additional COMPAS options not found in YAML template ###");// header record for YAML file output                     
                    }
                    ProcessOption(compasOptionStr,                                                                          // option name string (as written in the COMPAS OPTIONS class)
                                  optionDetails[optionIdx].valueStr,                                                        // option value string
                                  optionDetails[optionIdx].allowedStr,                                                      // option allowed value strings
                                  optionDetails[optionIdx].defaultStr,                                                      // option default string
                                  "",                                                                                       // no comment
                                  utils::Equals(optionDetails[optionIdx].sourceStr, "user_supplied"),                       // user specified option value?
                                  optionDetails[optionIdx].dataType,                                                        // option (short) data type  
                                  optionDetails[optionIdx].typeStr);                                                        // option (detailed) data type
                }
            }

            // Build the final output strings from the YamlRecord vector.
            // Two passes:
            //   1. compute the maximum widths of option, default, and allowed
            //      strings so the columns line up
            //   2. emit each record into a std::string, padding as required

            // pass 1: maximum widths
            size_t maxOptionStrLen  = 0;
            size_t maxDefaultStrLen = 0;
            size_t maxAllowedStrLen = 0;
            for (const YamlRecord& r : records) {
                if (r.isOption) {
                    const size_t optLen = r.preamble.length() + r.option.length() + 2 + r.value.length();   // "+2" for ": "
                    if (optLen > maxOptionStrLen) maxOptionStrLen = optLen;
                }
                if (r.defaultValue.length() > maxDefaultStrLen) maxDefaultStrLen = r.defaultValue.length();
                if (r.allowed.length()      > maxAllowedStrLen) maxAllowedStrLen = r.allowed.length();
            }

            // pass 2: emit padded strings
            const size_t defaultPos = maxOptionStrLen + 2;
            const size_t allowedPos = defaultPos + maxDefaultStrLen + 13;
            const size_t commentPos = allowedPos + maxAllowedStrLen + 4;
            std::vector<std::string> yamlRecords;
            yamlRecords.reserve(records.size());
            for (const YamlRecord& r : records) {
                std::string s;
                if (!r.isOption) {
                    s = r.preamble;                                                                                         // header/blank record - emit as-is
                }
                else {
                    s = r.preamble + r.option + ": " + r.value;                                                            // option record
                    if (s.length() < defaultPos) s.append(defaultPos - s.length(), ' ');                                    // pad to default column
                    s += "# Default: " + r.defaultValue;
                    if (!r.allowed.empty()) {                                                                               // allowed values?
                        if (s.length() < allowedPos) s.append(allowedPos - s.length(), ' ');
                        s += "# " + r.allowed;
                    }
                    if (!r.comment.empty()) {                                                                               // trailing comment?
                        if (s.length() < commentPos) s.append(commentPos - s.length(), ' ');
                        s += "  " + r.comment;
                    }
                }
                yamlRecords.push_back(s + "\n");
            }

            // write new YAML file

            const WRITE_FILE_RESULT writeResult = WriteYAMLfile(p_YAMLfilename, yamlRecords);                               // write records to file
            switch (writeResult) {
                case WRITE_FILE_RESULT::IO_ERROR:             std::cerr << "*ERROR* File '" << p_YAMLfilename << "' not written: IO error.\n";              break;
                case WRITE_FILE_RESULT::EXISTS_NO_OVERWRITE:  std::cout << "File '" << p_YAMLfilename << "' not written - file already exists.\n";          break;
                case WRITE_FILE_RESULT::OK:                   std::cout << "File '" << p_YAMLfilename << "' written.\n";                                   break;
            }
        }
    }
}
