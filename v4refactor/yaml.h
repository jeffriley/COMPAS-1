#pragma once

#include <algorithm>
#include <chrono>
#include <utility>
#include <vector>
#include <string>

namespace yaml {

// YAML template rules (in no particular order):
//
//  1. The following records will be automatically written to the start of YAML file:
//         ##~!!~## COMPAS option values
//         ##~!!~## Created at ddd MMM DD HH:MM:SS YYYY by COMPAS vxx.yy.zz
//         ##~!!~## 
//         ##~!!~## The default COMPAS YAML file (``compasConfigDefault.yaml``), as distributed, has
//         ##~!!~## all COMPAS option entries commented so that the COMPAS default value for the
//         ##~!!~## option is used by default. To use a value other than the COMPAS default value,
//         ##~!!~## users must uncomment the entry and change the option value to the desired value.
//
//  2. Lines in the template beginning with "##~!!~##"" will not be preserved (these are assumed to be COMPAS generated headers, and will be rewritten by COMPAS).
//  3. Leading '#' characters on option definition lines in the template will not be preserved (but they may be rewritten by COMPAS).
//  4. Option comments in the template must be preceded by "# " or they will not be preserved.
//  5. Strings in the template beginning with "# Default: " and up to (but not including) the next '#' (or end of line if no #) will not be preserved.
//  6. Strings in the template beginning with "# Options: " and up to (but not including) the next '#' (or end of line if no #) will not be preserved.
//  7. Blank lines in the template will be preserved.
//  8. Option values in the template will not be preserved (but they may be rewritten by COMPAS).
//  9. Option values written by COMPAS will be the option default values unless COMPAS was run with command-line options set - if the user executed COMPAS and
//     specified options on the command line, the user-specified values will be written to the YAML file, and those option records in the YAML file will not be
//     commented.  This gives users the option of creating project-specific YAML files via this method.
// 10. Options present in the template that are not valid COMPAS options will be ignored and not written to the YAML file.
// 11. Any COMPAS options that are not present in the template will be written in alphabetical order at the end of the YAML file.
//
// In the following example template:
//
// 0001     ##~!!~## COMPAS option values
// 0002     ##~!!~## File Created Tue Feb 14 13:09:06 2023 by COMPAS v02.34.06
// 0003     ##~!!~## 
// 0004     ##~!!~## The default COMPAS YAML file (``compasConfigDefault.yaml``), as distributed, has
// 0005     ##~!!~## all COMPAS option entries commented so that the COMPAS default value for the
// 0006     ##~!!~## option is used by default. To use a value other than the COMPAS default value,
// 0007     ##~!!~## users must uncomment the entry and change the option value to the desired value.
// 0008
// 0009     # first comment
// 0010
// 0011     booleanChoices:
// 0012         ### BINARY PROPERTIES
// 0013     #    --allow-touching-at-birth          # Default: False                                        # second comment
// 0014
// 0015         ### STELLAR PROPERTIES
// 0016         --mass-loss-prescription: 'HURLEY'  # Default: 'VINK'  # Options: ['VINK','HURLEY','NONE']    third comment
//
// Lines 0001 - 0007 will not be preserved (but will be replaced by new COMPAS headers).
// The blank line at line 0008 will be preserved.
// The comment "first comment" (on line 0009) will be preserved.
// The blank line at line 0010 will be preserved.
// The header "booleanChoices:" on line 0011 will be preserved.
// The header "### BINARY PROPERTIES" on line 0012 will be preserved.
// The leading '#' on line 0013 will not be preserved (but may be rewritten by COMPAS if the option is set to default).
// The string beginning with "# Default: " and extending to the next '#' on line 0013 will not be preserved (but will be replaced by COMPAS).
// The comment "second comment" on line 0013 will be preserved.
// The blank line at line 0014 will be preserved.
// The header "### STELLAR PROPERTIES" on line 0015 will be preserved.
// The string beginning with "# Default: " and extending to the next '#' on line 0016 will not be preserved (but will be replaced by COMPAS).
// The string beginning with "# Options: " and extending to the next '#' (or, in this case because there is no subsequent #, the end of the 
// line) on line 0016 will not be preserved (but will be replaced by COMPAS).
// The comment "third comment" on line 0016 will not be preserved - there is no "# " prefix, so it will be subsumed by the "# Options: " string
// (which extends from "# Options: " to the end of the line).


    // The default COMPAS YAML template is defined in yaml.cpp.
    // (Previously a literal vector inside an anonymous namespace in this header,
    // which gave every TU its own copy; moved to the .cpp for proper single-definition
    // semantics.)
    extern const std::vector<std::string> yamlTemplate;


    // namespace functions

    // Status enums for the internal functions. MakeYAMLfile is the only
    // public entry point; ReadYAMLtemplate and WriteYAMLfile are exposed
    // for testing and unusual uses, but most callers should use
    // MakeYAMLfile.

    enum class READ_TEMPLATE_RESULT {
        IO_ERROR,           // file present but could not be opened/read
        FILE_NOT_FOUND,     // file does not exist
        EMPTY,              // file read ok but contains no records
        OK                  // file read ok; content returned
    };

    enum class WRITE_FILE_RESULT {
        IO_ERROR,           // file could not be opened, or write failed mid-stream
        EXISTS_NO_OVERWRITE,// file exists and user declined to overwrite
        OK                  // file written successfully (possibly overwriting)
    };

    void                                         MakeYAMLfile(const std::string p_YAMLfilename, const std::string p_YAMLtemplate);
    std::pair<READ_TEMPLATE_RESULT, std::vector<std::string>>
                                                 ReadYAMLtemplate(const std::string p_YAMLtemplateName);
    WRITE_FILE_RESULT                            WriteYAMLfile(const std::string p_YAMLname, const std::vector<std::string> p_YAMLcontent);

}
