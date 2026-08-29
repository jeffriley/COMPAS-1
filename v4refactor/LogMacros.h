#pragma once

#include <sstream>

// =============================================================================
//
// The macros below are thin wrappers around variadic-template functions in the
// log_helpers:: namespace (see end of Log.h). The wrapper layer preserves the
// original compile-time-strip behaviour:
//
//   - DBG family is #ifdef DEBUG gated. When DEBUG is off, every DBG macro
//     expands to `((void)0)`; the call site and its argument expressions
//     vanish entirely from the source (no calls, no string formatting, no
//     side-effects from increment/decrement etc. - the production binary is
//     fully stripped of debug code).
//
//   - DBG_WARN family is #ifdef DEBUG_WARNINGS gated, same pattern.
//
//   - SAY / LOG / LOGV families are always enabled (they're production
//     logging, not debug-mode-only).
//
// Equivalent function forms in log_helpers:: are also available for callers
// that want call-evaluation regardless of the compile-time flag (the templates
// gate the BODY via `if constexpr` but the call and its arg evaluation still
// happen). Use the macro form for production debug instrumentation, the
// function form for callable contexts or template code.
// =============================================================================


//=============================================================================
// Compile-time switches.
//=============================================================================

//#define DEBUG // comment this line out, or #undef DEBUG, to build production executable (i.e. no DEBUG code)

#define DEBUG_WARNINGS  // comment this line out, or #undef DEBUG_WARNINGS, to build executable without WARNing statements

#define DEBUG_PERTURB // see discussion around FGB::PerturbLuminosityAndRadius()


//=============================================================================
// SAY family - thin wrapper around log_helpers::say.
//
// The macro is just a one-line wrapper that forwards to log_helpers::say
// (which is a variadic template - see end of Log.h).
//=============================================================================

#define SAY(...) do { std::stringstream _ss_say; _ss_say << __VA_ARGS__; ::log_helpers::say(_ss_say.str()); } while (0)


//=============================================================================
// DBG family - gated on #ifdef DEBUG.
//
// When DEBUG is on:   macros expand to log_helpers::xxx(...) calls.
// When DEBUG is off:  macros expand to ((void)0) - call site and arguments
//                     are completely removed from the source.
//=============================================================================

#ifdef DEBUG
    #define DBG(...)                              ::log_helpers::dbg(__VA_ARGS__)
    #define DBG_WITH(cls, lvl, ...)               ::log_helpers::dbg_with(cls, lvl __VA_OPT__(,) __VA_ARGS__)
    #define DBG_IF(cond, ...)                     ::log_helpers::dbg_if(cond __VA_OPT__(,) __VA_ARGS__)
    #define DBG_IF_WITH(cond, cls, lvl, ...)      ::log_helpers::dbg_if_with(cond, cls, lvl __VA_OPT__(,) __VA_ARGS__)
    #define DBG_WAIT(...)                         ::log_helpers::dbg_wait(__VA_ARGS__)
    #define DBG_WAIT_WITH(cls, lvl, ...)          ::log_helpers::dbg_wait_with(cls, lvl __VA_OPT__(,) __VA_ARGS__)
    #define DBG_WAIT_IF(cond, ...)                ::log_helpers::dbg_wait_if(cond __VA_OPT__(,) __VA_ARGS__)
    #define DBG_ID(...)                           ::log_helpers::dbg_id_impl(__PRETTY_FUNCTION__ __VA_OPT__(,) __VA_ARGS__)
    #define DBG_ID_IF(cond, ...)                  ::log_helpers::dbg_id_if_impl(__PRETTY_FUNCTION__, cond __VA_OPT__(,) __VA_ARGS__)
#else
    #define DBG(...)                              ((void)0)
    #define DBG_WITH(...)                         ((void)0)
    #define DBG_IF(...)                           ((void)0)
    #define DBG_IF_WITH(...)                      ((void)0)
    #define DBG_WAIT(...)                         ((void)0)
    #define DBG_WAIT_WITH(...)                    ((void)0)
    #define DBG_WAIT_IF(...)                      ((void)0)
    #define DBG_ID(...)                           ((void)0)
    #define DBG_ID_IF(...)                        ((void)0)
#endif


//=============================================================================
// DBG_WARN family - gated on #ifdef DEBUG_WARNINGS.
//=============================================================================

#ifdef DEBUG_WARNINGS
    #define DBG_WARN(...)                         ::log_helpers::dbg_warn(__VA_ARGS__)
    #define DBG_WARN_WITH(cls, lvl, ...)          ::log_helpers::dbg_warn_with(cls, lvl __VA_OPT__(,) __VA_ARGS__)
    #define DBG_WARN_IF(cond, ...)                ::log_helpers::dbg_warn_if(cond __VA_OPT__(,) __VA_ARGS__)
    #define DBG_WARN_ID(...)                      ::log_helpers::dbg_warn_id_impl(__PRETTY_FUNCTION__ __VA_OPT__(,) __VA_ARGS__)
    #define DBG_WARN_ID_IF(cond, ...)             ::log_helpers::dbg_warn_id_if_impl(__PRETTY_FUNCTION__, cond __VA_OPT__(,) __VA_ARGS__)
#else
    #define DBG_WARN(...)                         ((void)0)
    #define DBG_WARN_WITH(...)                    ((void)0)
    #define DBG_WARN_IF(...)                      ((void)0)
    #define DBG_WARN_ID(...)                      ((void)0)
    #define DBG_WARN_ID_IF(...)                   ((void)0)
#endif


//=============================================================================
// SAY auxiliary macros - always enabled (production code).
//=============================================================================

#define SAY_WITH(cls, lvl, ...)                   ::log_helpers::say_with(cls, lvl __VA_OPT__(,) __VA_ARGS__)
#define SAY_IF(cond, ...)                         ::log_helpers::say_if(cond __VA_OPT__(,) __VA_ARGS__)
#define SAY_ID(...)                               ::log_helpers::say_id_impl(__PRETTY_FUNCTION__ __VA_OPT__(,) __VA_ARGS__)
#define SAY_ID_IF(cond, ...)                      ::log_helpers::say_id_if_impl(__PRETTY_FUNCTION__, cond __VA_OPT__(,) __VA_ARGS__)


//=============================================================================
// LOG family - always enabled. Writes to a logfile only (no stdout echo).
//=============================================================================

#define LOG(logfileId, ...)                       ::log_helpers::log_to(logfileId __VA_OPT__(,) __VA_ARGS__)
#define LOG_WITH(logfileId, cls, lvl, ...)        ::log_helpers::log_to_with(logfileId, cls, lvl __VA_OPT__(,) __VA_ARGS__)
#define LOG_IF(logfileId, cond, ...)              ::log_helpers::log_if(logfileId, cond __VA_OPT__(,) __VA_ARGS__)
#define LOG_ID(logfileId, ...)                    ::log_helpers::log_id_impl(__PRETTY_FUNCTION__, logfileId __VA_OPT__(,) __VA_ARGS__)
#define LOG_ID_IF(logfileId, cond, ...)           ::log_helpers::log_id_if_impl(__PRETTY_FUNCTION__, logfileId, cond __VA_OPT__(,) __VA_ARGS__)


//=============================================================================
// LOGV family - always enabled. Writes to a logfile AND echoes to stdout.
//=============================================================================

#define LOGV(logfileId, ...)                      ::log_helpers::logv(logfileId __VA_OPT__(,) __VA_ARGS__)
#define LOGV_WITH(logfileId, cls, lvl, ...)       ::log_helpers::logv_with(logfileId, cls, lvl __VA_OPT__(,) __VA_ARGS__)
#define LOGV_IF(logfileId, cond, ...)             ::log_helpers::logv_if(logfileId, cond __VA_OPT__(,) __VA_ARGS__)
#define LOGV_ID(logfileId, ...)                   ::log_helpers::logv_id_impl(__PRETTY_FUNCTION__, logfileId __VA_OPT__(,) __VA_ARGS__)
#define LOGV_ID_IF(logfileId, cond, ...)          ::log_helpers::logv_id_if_impl(__PRETTY_FUNCTION__, logfileId, cond __VA_OPT__(,) __VA_ARGS__)
