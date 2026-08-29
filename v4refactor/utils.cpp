#include "utils.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <execinfo.h>
#include <cxxabi.h>

#include "changelog.h"
#include "hdf5.h"

#include <boost/version.hpp>
#include <boost/filesystem.hpp>

namespace utils {

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                       ADMIN                                       //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * Fail
 *
 * @brief
 * Prints an error message, shows a stack trace, and terminates the program.
 *
 *
 * void Fail(const ERROR p_Error, const StrT p_Preamble = "")
 *
 * @param   p_Error     The error to be printed (index into error catalog)
 * @param   p_Preamble  Preamble string to be printed prior to the error string.
 *                      Replaces the default preamble ("\nError retrieving state value\n").
 */
void Fail(const ERROR p_Error, const StrT p_Preamble) {

    StrT errStr;

    if (p_Preamble.empty()) {
        errStr = "\nError retrieving optional value\n";
    }
    else {
        if (p_Preamble.front() != '\n') errStr += "\n";
        errStr += p_Preamble;
        if (errStr.back() != '\n')      errStr += "\n";
    }

    const auto iter = ERROR_CATALOG.find(p_Error);
    if (iter != ERROR_CATALOG.end()) {
        errStr += iter->second.text + "\n";
    }
    else {
        errStr += "Unknown error\n";
    }

    std::cerr << errStr;
    utils::ShowStackTrace();
    std::exit(1);
}


/*
 * GetBoostVersion
 *
 * @brief
 * Get Boost version string in format MM.mm.rr
 *
 * 
 * StrT GetBoostVersion()
 * 
 * @return                                      String containing Boost version in format MM.mm.rr
 */
StrT GetBoostVersion() { 
    return std::to_string(BOOST_VERSION / 100000) + "." + std::to_string(BOOST_VERSION / 100 % 1000) + "." + std::to_string(BOOST_VERSION % 100);
}


/*
 * GetGslVersion
 *
 * @brief
 * Get gsl version string in format MM.mm.rr
 *
 * 
 * StrT GetGslVersion()
 * 
 * @return                                      String containing gsl version in format MM.mm.rr
 *                                              Will be "Not available" if not able to retrieve the actual value.
 */
StrT GetGslVersion() {

    StrT versionStr = "Not available";                                                                  // Default return value

    char buffer[128];                                                                                   // Command return buffer
    std::unique_ptr<FILE, int(*)(FILE*)> pipe(popen("gsl-config --version", "r"), &pclose);             // Open pipe for command
    if (pipe) {                                                                                         // Ok?
        versionStr = "";                                                                                // Yes
        while (fgets(buffer, sizeof(buffer), pipe.get()) != nullptr) versionStr += buffer;              // Copy buffer
        if (!versionStr.empty() && versionStr[versionStr.length() - 1] == '\n') versionStr.pop_back();  // Remove trailing newline if necessary
    }

    return versionStr;
}


/*
 * GetHDF5Version
 *
 * @brief
 * Get HDF5 library version string in format MM.mm.rr
 *
 * 
 * StrT GetHDF5Version()
 * 
 * @return                                      String containing HDF5 library version in format MM.mm.rr
 *                                              Will be "Not available" if not able to retrieve the actual value.
 */
StrT GetHDF5Version() {

    StrT versionStr = "Not available";                                                                              // Default return value

    unsigned majorNum, minorNum, releaseNum;
    herr_t status = H5get_libversion(&majorNum, &minorNum, &releaseNum);                                            // Retrieve HDF5 library version
    if (status >= 0)                                                                                                // Ok?
        versionStr = std::to_string(majorNum) + "." + std::to_string(minorNum) + "." + std::to_string(releaseNum);  // Yes - set version string

    return versionStr;
}


/*
 * SplashScreen
 *
 * @brief 
 * Construct COMPAS splash string.
 * Prints splash string to stdout if required.
 *
 *
 * StrT SplashScreen(const bool p_Print)
 * 
 * @param       p_Print                         Boolean indicating whether splash string should be printed.  Default is TRUE.
 * @return                                      Splash string
 */
StrT SplashScreen(const bool p_Print) {

    // Set compiler name
    StrT compilerName = COMPILER_ID == GCC_COMPILER ? "gcc" : (COMPILER_ID == CLANG_COMPILER ? "Clang" : (COMPILER_ID == MSC_COMPILER ? "MSCV" : "?"));

    // Construct compiler details string
    StrT major = COMPILER_MAJOR < 0 ? "?" : std::to_string(COMPILER_MAJOR);
    StrT minor = COMPILER_MINOR < 0 ? "?" : std::to_string(COMPILER_MINOR);
    StrT patch = COMPILER_PATCH < 0 ? "?" : std::to_string(COMPILER_PATCH);

    StrT compilerStr = compilerName + " " + major + "." + minor + "." + patch;

    // Construct the splash string
    StrT splashString = "\nCOMPAS v" + 
                        VERSION_STRING + 
                        " (" + compilerStr +
                        ", gsl v" + utils::GetGslVersion() + 
                        ", boost v" + utils::GetBoostVersion() + 
                        ", HDF5 v" + utils::GetHDF5Version() + ")" +
                        "\nCompact Object Mergers: Population Astrophysics and Statistics"
                        "\nby Team COMPAS (http://compas.science/index.html)"
                        "\nA binary star simulator\n"
                        "\nGo to https://compas.readthedocs.io/en/latest/index.html for the online documentation"
                        "\nCheck https://compas.readthedocs.io/en/latest/pages/whats-new.html to see what's new in the latest release\n";

    if (p_Print) std::cout << splashString << std::endl;    // Print the splash string if required

    return splashString;                                    // Return the splash string
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                  ERROR HANDLING                                   //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

    
/*
 * GetStackTrace
 *
 * @brief
 * Construct a vector of strings that represent the stack trace for the current thread
 * (COMPAS is single-threaded, so in our case, the current process executing COMPAS).
 * 
 * We use the gcc library functions to construct the stack trace (gcc calls this a "backtrace"):
 * 
 *     - backtrace(), which provides a list of pointers to each of the functions that make
 *       up the stack trace (the functions called, all the way from main() to the current
 *       point of execution - this is an instantaneous list, not historic).
 * 
 * 
 *     - backtrace_symbols(), translates the function pointers obtained from backtrace() into
 *       an array of strings that are the function names.  This typically only works for COMPAS
 *       functions (because we build COMPAS with debug info included).  We most likely won't have
 *       symbols for libraries (e.g. libc), so for non-COMPAS functions we insert the string
 *       "~~LIBFUNC~~" as the function name so users can identify non-COMPAS entries and handle
 *       them accordingly.
 * 
 * Returns a vector of strings representing the function names that comprise the stack trace.
 * 
 * 
 * StrVectorT GetStackTrace()
 * 
 * @return                                      Vector of strings representing the function names that comprise the stack trace
 */
StrVectorT GetStackTrace() {

    StrVectorT stackTrace = {};                                                                                 // Return vector

    void*  trace[MAX_STACK_TRACE_SIZE];                                                                         // Stack trace
    char** strings   = (char **)NULL;                                                                           // Stack trace strings
    SizeT  traceSize = 0;                                                                                       // Stack trace size
        
    traceSize = backtrace(trace, MAX_STACK_TRACE_SIZE);                                                         // Get stack trace size
    strings   = backtrace_symbols(trace, traceSize);                                                            // Get stack trace with symbols

    for (SizeT idx = 1; idx < traceSize; ++idx) {                                                               // For each stack trace entry
        // Extract function name.
        // We don't have symbols for libraries (e.g. libc), so for non-COMPAS functions we insert
        // "~~LIBFUNC~~" as the function name so the caller can identify non-COMPAS entries and
        // handle them accordingly
        SizeT start = 0;
        while (strings[idx][start] != '(' && strings[idx][start] != ' ' && strings[idx][start] != 0) ++start;   // Find function name start position
        SizeT end = start;
        while (strings[idx][end] != '+' && strings[idx][end] != 0) ++end;                                       // Find function name end position

        StrT funcName;                                                                                          // The extracted function name
        SizeT funcStrLen = end - start - 1;                                                                     // Length of (mangled) function string
        if (funcStrLen < 1) funcName = "~~LIBFUNC~~";                                                           // Library function
        else {                                                                                                  // Extract COMPAS function name
            char* funcStr = new char[funcStrLen + 1];                                                           // Allows for null terminator
            strncpy(funcStr, &strings[idx][start + 1], funcStrLen);                                             // Copy function name
            funcStr[funcStrLen] = 0;                                                                            // Make sure it is null-terminated
            funcName = StrT(funcStr);                                                                           // Function name

            int status = -1;
            char* demangledName = abi::__cxa_demangle(funcStr, NULL, NULL, &status);                            // Try to demangle the function name
            if (status == 0) funcName = StrT(demangledName);                                                    // Use the demangled name if available

            delete[] demangledName;
            delete[] funcStr;
        }

        stackTrace.push_back(funcName);                                                                         // Add function name to stacktrace

        if (funcName == "main") break;                                                                          // That's all we need
    }
    delete[] strings;

    return stackTrace;
}


/*
 * ShowStackTrace
 *
 * @brief
 * Display a stack trace, obtained by calling utils::GetStackTrace(), on stderr
 *
 * 
 * void ShowStackTrace()
 */
void ShowStackTrace() {

    StrVectorT stackTrace = utils::GetStackTrace();                     // Get stack trace

    if (!stackTrace.empty()) {                                          // Anything to show?
        std::cerr << "\nStack trace:\n";                                // Yes - display header
        for (SizeT entry = 1; entry < stackTrace.size(); entry++) {     // Ignore the eponymous entry
            std::cerr << "    " << stackTrace[entry] << "\n";           // Show stacktrace entry
        }
    }
    else {                                                              // No - nothing to show
        std::cerr << "\nNo stack trace available\n";                    // Say so
    }
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                    FILE SYSTEM                                    //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CreateDirectory
 *
 * @brief
 * Create directories, where they don't already exist, in the path supplied
 *
 * We could just use boost::create_directories() (or FS::create_directories() for c++17 or greater), but that
 * would just create the directories as necessary and not report which directories in the path were pre-existing
 * and which were newly created.  We want to be able to clean up any directories we created but didn't use (because
 * we had an error somewhere perhaps...), and to do that we need to know which directories were actually created.
 * 
 * Returns a vector of paths created.
 * 
 * 
 * std::tuple<ERROR, StrT, StrVectorT> CreateDirectory(const StrT p_Path)
 * 
 * @param       p_Path                          Path specifying directories to be created
 * @return                                      Tuple containing error value, error string, and vector of directories created
 *                                                  The error value returned will be:
 *                                                      NONE                       if no error occurred
 *                                                      UNABLE_TO_CREATE_DIRECTORY if any non-existent directory in the path supplied could not be created
 *                                                  The error string indicates the problematic path in the case an error occurred.
 *                                                      The error string is only valid if error != NONE
 *                                                  The returned vector contains string paths for each of the directories actually created.
 *                                                      The returned vector will not contain names of directories in the path that already existed
 *                                                      The content of the returned vector for directories created is only valid if error == NONE,
 */
std::tuple<ERROR, StrT, StrVectorT> CreateDirectories(const StrT p_Path) {

    if (p_Path.empty()) return std::make_tuple(ERROR::NONE, "", StrVectorT({}));            // Nothing to do

    ERROR error             = ERROR::NONE;                                                  // Error - initially NONE
    StrT errStr             = "";                                                           // Error string - initially none
    StrVectorT pathsCreated = {};                                                           // Directories created - initially none

    // Create directories as necessary - top-down
    const char sep   = boost::filesystem::path::preferred_separator;                        // Platform-specific path separator
    StrT path       = "";                                                                   // Current path - the directory to be created
    for (const auto& name: boost::filesystem::path(p_Path)) {                               // Parse the user-supplied path string
        if (!path.empty() && path[path.length() - 1] != sep) path += sep;                   // Add separator to current path if necessary
        path += name.c_str();                                                               // Add directory name to current path
        if (!path.empty() && !boost::filesystem::exists(path)) {
            try {
                boost::system::error_code err;
                (void)boost::filesystem::create_directory(path, err);                       // Create directory - let boost throw an exception if it fails
                if (err.value() == 0) {                                                     // Ok?
                    pathsCreated.push_back(boost::filesystem::canonical(path).string());    // Yes - record creation
                }
                else  {                                                                     // Not ok...
                    error  = ERROR::UNABLE_TO_CREATE_DIRECTORY;                             // Set error
                    errStr = path;                                                          // ... and error string
                }            
            }
            catch (...) {                                                                   // Unhandled problem...
                error  = ERROR::UNABLE_TO_CREATE_DIRECTORY;                                 // Set error
                errStr = path;                                                              // ... and error string
            }                
        }

        if (error != ERROR::NONE) {                                                         // Problem?
            std::tie(error, errStr, pathsCreated) = RemoveDirectories(pathsCreated);        // Yes - clean up
            break;
        }
    }

    return std::make_tuple(error, errStr, pathsCreated);
}


/* 
 * FileExists
 *
 * @brief
 * Determine if a file with name as passed in p_Filename exists.
 * The parameter p_Filename is a character array.
 *
 *
 * bool FileExists(const char *p_Filename)
 *
 * @param       p_Filename                      Fully qualified filename (char*)
 * @return                                      Boolean indicating whether file exists
 */
bool FileExists(const char *p_Filename) { return (bool)std::ifstream(p_Filename); }


/*
 * FileExists
 *
 * @brief
 * Determine if a file with name as passed in p_Filename exists.
 * The parameter p_Filename is a std::string.
 *
 *
 * bool FileExists(const std::string& p_Filename)
 *
 * @param       p_Filename                      Fully qualified filename (std::string)
 * @return                                      Boolean indicating whether file exists
 */
bool FileExists(const StrT& p_Filename) { return FileExists(p_Filename.c_str()); }


/*
 * ReadTimesteps
 *
 * @brief
 * Read timesteps from timesteps file.
 *
 * Timesteps file is expected to be an ascii file with one timestep per record.
 * Timesteps must be > 0.0
 *  
 * 
 * std::tuple<ERROR, DblVectorT> ReadTimesteps(const StrT p_TimestepsFileName)
 * 
 * @param       p_TimestepsFileName             Filename to be read - should be fully qualified
 * @return                                      Tuple containing error value and timesteps vector
 *                                                  The error value returned will be:
 *                                                      NONE                                 if no error occurred
 *                                                      EMPTY_FILENAME                       if the filename provided was an empty string
 *                                                      FILE_DOES_NOT_EXIST                  if the timesteps file does not exist
 *                                                      FILE_OPEN_ERROR                      if the timesteps file exists but could not be opened
 *                                                      FILE_READ_ERROR                      if the timesteps file could not be read
 *                                                      EMPTY_FILE                           if the timesteps file contains no content
 *                                                      INVALID_VALUE_IN_FILE                if the file contains an invalid value for timestep
 *                                                      TOO_MANY_TIMESTEPS_IN_TIMESTEPS_FILE if the file contains too many timesteps (> maximum per OPTIONS)
 * 
 *                                                  If the error returned is not NONE, the content of the timesteps vector returned is not defined
 */
std::tuple<ERROR, DblVectorT> ReadTimesteps(const StrT p_TimestepsFileName) {

    ERROR error = ERROR::NONE;
 
    DblVectorT timesteps;                                                                                       // Timesteps vector

    if (p_TimestepsFileName.empty()) {                                                                          // Timesteps filename empty?
        error = ERROR::EMPTY_FILENAME;                                                                          // Yes - fail
    }
    else {

        if (!utils::FileExists(p_TimestepsFileName)) {                                                          // Timesteps file exists?
            error = ERROR::FILE_DOES_NOT_EXIST;                                                                 // No - fail
        }
        else {                                                                                                  // Yes
            std::ifstream timestepsFile;
            timestepsFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            try {
                timestepsFile.open(p_TimestepsFileName);                                                        // Open the timesteps file
                if (!timestepsFile.is_open()) {                                                                 // Open ok?
                    error = ERROR::FILE_OPEN_ERROR;                                                             // No - fail
                }
                else {                                                                                          // Yes - file open
                    std::string rec;                                                                            // Record read from file
                    SizeT numTimesteps = 0;
                    while (std::getline(timestepsFile, rec)) {                                                  // Get next record from timesteps file

                        if (rec.size() > 0 && (rec[rec.size() - 1] == '\n' || rec[rec.size() - 1] == '\r')) {   // Last character `\n` or `\r`?
                            rec.erase(rec.size() - 1);                                                          // Yes - strip it
                        }

                        utils::trim(rec);                                                                 // Remove leading and trailing blanks

                        if (!(rec.empty() || rec[0] == '#')) {                                                  // Blank record or comment?                                 
                            try {                                                                               // No - process it
                                SizeT lastChar;
                                long double v = std::stold(rec, &lastChar);                                     // Try conversion
                                if (lastChar != (rec.size())) {                                                 // Conversion valid only if rec completely consumed
                                    error = ERROR::INVALID_VALUE_IN_FILE;                                       // Not a valid DOUBLE
                                    break;                                                                      // Stop processing
                                }

                                if (v < 0.0) {                                                                  // Timestep must be >= 0.0
                                    error = ERROR::INVALID_VALUE_IN_FILE;                                       // Not a valid timestep
                                    break;                                                                      // Stop processing
                                }
                                else {                                                                          // Ok - timestep >= 0.0
                                    timesteps.push_back(v);                                                     // Add timestep to timesteps vector
                                }
                    
                                numTimesteps++;                                                                 // Increment number of timesteps read
                                if (numTimesteps >= ABSOLUTE_MAXIMUM_TIMESTEPS) {                               // Number of timesteps exceeds maximum?
                                    error = ERROR::TOO_MANY_TIMESTEPS_IN_TIMESTEPS_FILE;                        // Yes - fail
                                    break;                                                                      // Stop processing
                                }
                            }
                            catch (const std::out_of_range& e) {                                                // Conversion failed
                                error = ERROR::INVALID_VALUE_IN_FILE;                                           // Not a valid DOUBLE
                                break;                                                                          // Stop processing
                            }
                            catch (const std::invalid_argument& e) {                                            // Conversion failed
                                error = ERROR::INVALID_VALUE_IN_FILE;                                           // Not a valid DOUBLE
                                break;                                                                          // Stop processing
                            }
                        }
                    }
                    try {
                        timestepsFile.close();                                                                  // Close the timesteps file
                    }
                    catch (std::ifstream::failure& e) {                                                         // Close failed
                        error = ERROR::FILE_NOT_CLOSED;                                                         // Fail
                    }
                }
            }
            catch (std::ifstream::failure& e) {                                                                 // Something was flagged...
                if (timestepsFile.eof()) {                                                                      // end-of-file?
                    if (timesteps.size() < 1) {                                                                 // Yes - at least one timestep read?
                        error = ERROR::EMPTY_FILE;                                                              // No - fail
                    }
                }
                else {                                                                                          // Not end-of-file - error
                    error = ERROR::FILE_READ_ERROR;                                                             // Fail
                }
            }

        }
    }

    return std::make_tuple(error, timesteps);
}


/*
 * RemoveDirectories
 *
 * @brief
 * Remove directories if they are empty
 *
 * Iterate in reverse order over the vector of paths passed (i.e. walk up the directory tree from the leaf) and
 * remove empty directories - but stop at the first non-empty directory (or error).  
 * 
 * Returns a vector of paths not removed.
 * 
 * 
 * std::tuple<ERROR, StrT, StrVectorT> RemoveDirectories(const StrVectorT p_Paths)
 * 
 * @param       p_Path                          Vector of paths to be removed
 * @return                                      Tuple containing error value, error string, and vector of directories not removed
 *                                                  The error value returned will be:
 *                                                      NONE                       if no error occurred
 *                                                      UNABLE_TO_CREATE_DIRECTORY if any non-existent directory in the path supplied could not be created
 *                                                  The error string indicates the problematic path in the case an error occurred.
 *                                                      The error string is only valid if error != NONE
 *                                                  The returned vector contains string paths for each of the directories from the input vector that were not removed,
 *                                                  either because they were not empty or because an error occurred (check error element of tuple returned)
 */
std::tuple<ERROR, StrT, StrVectorT> RemoveDirectories(const StrVectorT p_Paths) {

    if (p_Paths.size() < 1) return std::make_tuple(ERROR::NONE, "", StrVectorT({}));                                    // Nothing to do

    ERROR      error           = ERROR::NONE;                                                                           // Error - initially NONE
    StrT       errStr          = "";                                                                                    // Error string - initially none
    StrVectorT pathsNotRemoved = p_Paths;                                                                               // Directories not removed - initially p_Paths

    // Remove directories as necessary - bottom-up - stop if error or not empty.
    StrT path;
    boost::system::error_code err;
    for (auto iter = p_Paths.rbegin(); iter != p_Paths.rend(); ++iter) {
        path = *iter;
        if (boost::filesystem::exists(path) && boost::filesystem::is_empty(path)) {                                     // Directory exists and is empty?
            try {
                (void)boost::filesystem::remove(path, err);                                                             // Remove directory - let boost throw an exception if it fails
                if (err.value() == 0) {                                                                                 // Ok?
                                                                                                                        // Yes
                    // Erase path removed from pathsNotRemoved vector.
                    // I could just use std::find() as the parameter to erase(), but this is safer.
                    StrVectorT::iterator idx = std::find(pathsNotRemoved.begin(), pathsNotRemoved.end(), path);         // Get element index
                    if (idx != pathsNotRemoved.end()) pathsNotRemoved.erase(idx);                                       // If found, erase it
                }
                else {                                                                                                  // Not ok...
                    error  = ERROR::UNABLE_TO_REMOVE_DIRECTORY;                                                         // Set error
                    errStr = path;                                                                                      // ... and error string
                    break;
                }
            }
            catch (...) {                                                                                               // Unhandled problem...
                error  = ERROR::UNABLE_TO_REMOVE_DIRECTORY;                                                             // Set error
                errStr = path;                                                                                          // ... and error string
                break;
            }                
        }
        else {                                                                                                          // Problem...
            // Either the directory does not exist or is not empty.
            // if not empty we stop here - we don't want to delete anything other than empty directories.
            // if the directory does not exist we err on the side of caution, assume there's a problem, and stop here.
            error  = boost::filesystem::exists(path) ? ERROR::UNABLE_TO_REMOVE_DIRECTORY : ERROR::DIRECTORY_NOT_EMPTY;  // Set error
            errStr = path;                                                                                              // ... and error string
            break;
        }
    }

    return std::make_tuple(error, errStr, pathsNotRemoved);
}
 

///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                     NUMERICS                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * IsBOOL
 *
 * @brief
 * Determine if the string passed as p_Str is a valid BOOL (as defined by Boost).
 *
 * In this context (the Boost context), a valid boolean is one of:
 * 
 *     - 0|1        ("0" or "1")
 *     - true|false ("true" or "false" - case insensitive)
 *     - yes|no     ("yes" or "no" - case insensitive)
 *     - on|off     ("on" or "off" - case insensitive)
 *
 * The function will retiurn one of {0, 1, 2, 3, 4, -1, -2, -3, -4} to indicate the result:
 * 
 *     0 = not a valid boolean
 *     1 = valid: 0|1
 *     2 = valid: true|false
 *     3 = valid: yes|no
 *     4 = valid: on|off
 * 
 * A positive return value indicates the boolean value is TRUE; a negative, FALSE.
 * 
 * 
 * SizeT IsBOOL(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - as described above
 */
SizeT IsBOOL(const StrT p_Str) {

    if (p_Str.empty()) return 0;                    // Not valid: empty string

    if (utils::Equals(p_Str, "1")    ) return  1;   // Valid: 0|1       : TRUE
    if (utils::Equals(p_Str, "0")    ) return -1;   // Valid: 0|1       : FALSE
    if (utils::Equals(p_Str, "true") ) return  2;   // Valid: true|false: TRUE
    if (utils::Equals(p_Str, "false")) return -2;   // Valid: true|false: FALSE
    if (utils::Equals(p_Str, "yes")  ) return  3;   // Valid: yes|no    : TRUE
    if (utils::Equals(p_Str, "no")   ) return -3;   // Valid: yes|no    : FALSE
    if (utils::Equals(p_Str, "on")   ) return  4;   // Valid: on|off    : TRUE
    if (utils::Equals(p_Str, "off")  ) return -4;   // Valid: on|off    : FALSE

    return 0;                                       // Not valid
}


/*
 * IsDOUBLE
 *
 * @brief
 * Determine if the string passed as p_Str is a valid DOUBLE.
 *
 * In this context, to be a valid DOUBLE the string must convert to a double successfully
 * via the std::stod() function.
 * 
 * 
 * bool IsDOUBLE(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid DOUBLE, else FALSE
 */
bool IsDOUBLE(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stod(p_Str, &lastChar);      // Try conversion

        result = lastChar == p_Str.size();      // Valid DOUBLE if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid DOUBLE
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid DOUBLE
    }

    return result;
}


/*
 * IsFLOAT
 *
 * @brief
 * Determine if the string passed as p_Str is a valid FLOAT.
 *
 * In this context, to be a valid FLOAT the string must convert to a float successfully
 * via the std::stof() function.
 * 
 * 
 * bool IsFLOAT(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid FLOAT, else FALSE
 */
bool IsFLOAT(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stof(p_Str, &lastChar);      // Try conversion

        result = lastChar == p_Str.size();      // Valid FLOAT if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid FLOAT
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid FLOAT
    }

    return result;
}


/*
 * IsINT
 *
 * @brief
 * Determines if the string passed as p_Str is a valid INT.
 *
 * In this context, to be a valid INT the string must convert to an integer successfully
 * via the std::stoi() function.
 * 
 * 
 * bool IsINT(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid INT, else FALSE
 */
bool IsINT(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stoi(p_Str, &lastChar);      // Try conversion

        result = lastChar == p_Str.size();      // Valid INT if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid INT
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid INT
    }

    return result;
}


/*
 * IsLONGDOUBLE
 *
 * @brief
 * Determine if the string passed as p_Str is a valid LONG DOUBLE.
 *
 * In this context, to be a valid LONG DOUBLE the string must convert to a long double successfully
 * via the std::stold() function.
 * 
 * 
 * bool IsLONGDOUBLE(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid LONG DOUBLE, else FALSE
 */
bool IsLONGDOUBLE(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stold(p_Str, &lastChar);     // Try conversion

        result = lastChar == p_Str.size();      // Valid LONG DOUBLE if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid LONG DOUBLE
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid LONG DOUBLE
    }

    return result;
}


/*
 * IsLONGINT
 *
 * @brief
 * Determine if the string passed as p_Str is a valid LONG INT.
 *
 * In this context, to be a valid LONG INT the string must convert to a long integer
 * successfully via the std::stol() function.
 * 
 * 
 * bool IsLONGINT(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid LONGINT, else FALSE
 */
bool IsLONGINT(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stol(p_Str, &lastChar);      // Try conversion

        result = lastChar == p_Str.size();      // Valid LONG INT if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid LONG INT
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid LONG INT
    }

    return result;
}


/*
 * Determine if the string passed as p_Str is a valid UNSIGNED LONG INT.
 *
 * In this context, to be a valid UNSIGNED LONG INT the string must convert to an unsigned long integer
 * successfully via the std::stoul() function.
 * 
 * 
 * bool IsULONGINT(const StrT p_Str)
 *
 * @param       p_Str                           String to check
 * @return                                      Result - TRUE if string is a valid UNSIGNED LONGINT, else FALSE
 */
bool IsULONGINT(const StrT p_Str) {

    bool result = false;                        // Default result

    try {
        SizeT lastChar;
        (void)std::stoul(p_Str, &lastChar);     // Try conversion

        result = lastChar == p_Str.size();      // Valid UNSIGNED LONG INT if p_Str completely consumed
    }
    catch (const std::out_of_range& e) {        // Conversion failed
        result = false;                         // Not a valid UNSIGNED LONG INT
    }
    catch (const std::invalid_argument& e) {    // Conversion failed
        result = false;                         // Not a valid UNSIGNED LONG INT
    }

    return result;
}


///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      SOLVERS                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////


/*
 * SolveQuadratic
 *
 * @brief
 * Solve quadratic Ax^2 + Bx + C.
 *
 * Returns either root, depending on discriminant will return:
 *
 *    0.0               if 0 roots
 *    root              if 1 root
 *    max(root1, root2) if 2 roots
 *
 *
 * std::tuple<ERROR, double> SolveQuadratic(const double p_A, const double p_B, double p_C)
 *
 * @param       p_A                             Coefficient of x^2
 * @param       p_B                             Coefficient of x^1
 * @param       p_C                             Coefficient of x^0 (Constant)
 * @return                                      Tuple containing (in order): error value, root found (see above)
 *                                                  The error value returned will be:
 *                                                      NONE          if no error occurred
 *                                                      NO_REAL_ROOTS if the equation has no real roots
 *                                                  If the error returned is not ERROR:NONE, use root returned at your own risk
 */
std::tuple<ERROR, double> SolveQuadratic(const double p_A, const double p_B, double p_C) {

    ERROR error = ERROR::NONE;

    const double discriminant = (p_B * p_B) - (4.0 * p_A * p_C);    // d = B^2 - 4AC

    double root = 0.0;                                              // Root found

    if (discriminant < 0.0) {                                       // Real roots?
        error = ERROR::NO_REAL_ROOTS;                               // No real roots - set error
    }
    else if (discriminant > 0.0) {                                  // 2 real roots?
                                                                    // 2 real roots
        const double sqrtD = std::sqrt(discriminant);
        const double A2    = p_A + p_A;
        const double root1 = (-p_B + sqrtD) / A2;                   // (-B + SQRT(B^2 - 4AC)) / 2A
        const double root2 = (-p_B - sqrtD) / A2;                   // (-B - SQRT(B^2 - 4AC)) / 2A

        root = std::max(root1, root2);
    }
    else {                                                          // 1 real root
        root = -p_B / (p_A + p_A);                                  // -B / 2A,discriminant = 0.0
    }

    return std::make_tuple(error, root);
}

    
///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
//                                      STRINGS                                      //
//                                                                                   //
///////////////////////////////////////////////////////////////////////////////////////

/*
 * CentreJustify
 *
 * @brief
 * Centre-justifies string to specified width by prepending and appending spaces.
 * Extra space will be at the end of the string if necessary.
 * Used primarily in logging code.
 *
 *
 * StrT CentreJustify(const StrT p_Str, const SizeT p_Width)
 *
 * @param       p_Str                           String to be centre-justified
 * @param       p_Width                         The required width of the resultant string
 * @return                                      String padded with leading and trailing spaces so as to (as close as possible) centre-justify p_Str
 *                                              The string returned will always be p_Width characters in length
 */
StrT CentreJustify(const StrT p_Str, const SizeT p_Width) {

    StrT result = p_Str;                                                        // Default is no change

    if (p_Str.length() < p_Width) {                                             // p_Str < field width?
                                                                                // Yes
        SizeT numLeadingSpaces = (p_Width - p_Str.length()) / 2;                // Number of spaces to add at start - half the deficit
        SizeT numTralingSpaces = p_Width - p_Str.length() - numLeadingSpaces;   // Number of spaces to add at the end - whatever is left
        StrT leadingSpaces(numLeadingSpaces, ' ');                              // Blank string to prepend
        StrT trailingSpaces(numTralingSpaces, ' ');                             // Blank string to append
        result = leadingSpaces + p_Str + trailingSpaces;                        // Add leading and trailing spaces to p_Str
    }

    return result;
}


/*
 * Equals
 *
 * @brief
 * Case-insensitive comparison of strings.
 * This only works with ASCII data, but I think that's all we need.
 * Note that std::string has an == operator to test for equality (actually calls std::strcmp)
 *
 *
 * bool Equals(StrT p_Str1, StrT p_Str2)
 *
 * @param       p_Str1                          String to be compared
 * @param       p_Str2                          String to be compared
 * @return                                      Boolean indicating equality (true = equal)
 */
bool Equals(StrT p_Str1, StrT p_Str2) {
    std::transform(p_Str1.begin(), p_Str1.end(), p_Str1.begin(), ::tolower);
    std::transform(p_Str2.begin(), p_Str2.end(), p_Str2.begin(), ::tolower);
    return (std::strcmp(p_Str1.c_str(), p_Str2.c_str()) == 0);
}


/*
 * PadLeadingZeros
 *
 * @brief
 * Pads string to specified length by prepending the string with "0".
 * This only works with ASCII data, but I think that's all we need.
 *
 *
 * StrT PadLeadingZeros(const StrT p_Str, const SizeT p_MaxLength)
 *
 * @param       p_Str                           String to be padded with leading "0"s
 * @param       p_MaxLength                     The required length of the resultant string
 * @return                                      String padded with leading "0"s
 *                                              Will be unchanged from input string if length already >= required length
 */
StrT PadLeadingZeros(const StrT& p_Str, const SizeT p_MaxLength) {
    return (p_Str.length() < p_MaxLength) ? StrT(p_MaxLength - p_Str.length(), '0') + p_Str : p_Str;
}


/*
 * PadTrailingSpaces
 *
 * @brief
 * Pads string to specified length by appending the string with " ".
 * This only works with ASCII data, but I think that's all we need.
 *
 *
 * StrT PadTrailingSpaces(const StrT p_Str, const SizeT p_MaxLength)
 *
 * @param       p_Str                           String to be padded with trailing " "s
 * @param       p_MaxLength                     The required length of the resultant string
 * @return                                      String padded with trailing " "s
 *                                              Will be unchanged from input string if length already >= required length
 */
StrT PadTrailingSpaces(const StrT& p_Str, const SizeT p_MaxLength) {
    return (p_Str.length() < p_MaxLength) ? p_Str + StrT(p_MaxLength - p_Str.length(), ' ') : p_Str;
}


/*
 * ltrim
 *
 * @brief
 * Trim leading whitespace characters from a string.
 * 
 * Usage is either of:
 * 
 *    s = trim(s)
 *    trim(s)
 * 
 * ('s' is trimmed in both cases)
 *
 *
 * StrT& ltrim(StrT& p_Str)
 *
 * @param       p_Str                           String to be trimmed of whitespace
 * @return                                      Trimmed string
 */
StrT& ltrim(StrT& p_Str) {
    p_Str.erase(0, p_Str.find_first_not_of("\t\n\v\f\r "));
    return p_Str;
}


/*
 * rtrim
 *
 * @brief
 * Trim trailing whitespace characters from a string.
 *
 * Usage is either of:
 * 
 *    s = trim(s)
 *    trim(s)
 * 
 * ('s' is trimmed in both cases)
 * 
 *
 * StrT& rtrim(StrT& p_Str)
 *
 * @param       p_Str                           String to be trimmed of whitespace
 * @return                                      Trimmed string
 */
StrT& rtrim(StrT& p_Str) {
    p_Str.erase(p_Str.find_last_not_of("\t\n\v\f\r ") + 1);
    return p_Str;
}


/*
 * trim
 *
 * @brief
 * Trim both leading and trailing whitespace characters from a string.
 *
 * Usage is either of:
 * 
 *    s = trim(s)
 *    trim(s)
 * 
 * ('s' is trimmed in both cases)
 * 
 *
 * StrT& trim(StrT& p_Str)
 *
 * @param       p_Str                           String to be trimmed of whitespace
 * @return                                      Trimmed string
 */
StrT& trim(StrT& p_Str) { return utils::ltrim(utils::rtrim(p_Str)); }


/*
 * ToLower
 *
 * @brief
 * Downshift uppercase characters in string to lowercase.
 * This only works with ASCII data, but I think that's all we need,
 *
 *
 * StrT ToLower(StrT p_Str)
 *
 * @param       p_Str                           String to be downshifted
 * @return                                      Downshifted string
 */
StrT ToLower(StrT p_Str) {
    std::transform(p_Str.begin(), p_Str.end(), p_Str.begin(), ::tolower);
    return p_Str;
}


/*
 * ToUpper
 *
 * @brief
 * Upshift lowercase characters in string to uppercase.
 * This only works with ASCII data, but I think that's all we need.
 *
 *
 * StrT ToUpper(StrT p_Str)
 *
 * @param       p_Str                           String to be upshifted
 * @return                                      Upshifted string
 */
StrT ToUpper(StrT p_Str) {
    std::transform(p_Str.begin(), p_Str.end(), p_Str.begin(), ::toupper);
    return p_Str;
}

} // namespace utils
