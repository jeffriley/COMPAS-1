#pragma once

#include "ErrorCatalog.h"   // ERROR enum
#include "typedefs.h"       // OBJECT_ID, OBJECT_TYPE, STELLAR_TYPE
#include <string>


#define ERROR_PREFIX   "ERROR: "
#define WARNING_PREFIX "WARNING: "

namespace error_helpers {
    void show_error(ERROR e,                          OBJECT_ID, OBJECT_TYPE, STELLAR_TYPE, const char*);
    void show_error(ERROR e, const std::string& qual, OBJECT_ID, OBJECT_TYPE, STELLAR_TYPE, const char*);
    void show_warn (ERROR e,                          OBJECT_ID, OBJECT_TYPE, STELLAR_TYPE, const char*);
    void show_warn (ERROR e, const std::string& qual, OBJECT_ID, OBJECT_TYPE, STELLAR_TYPE, const char*);
    void throw_error(ERROR e,                          OBJECT_ID, OBJECT_TYPE, STELLAR_TYPE, const char*);
    void throw_error(ERROR e, const std::string& qual, OBJECT_ID, OBJECT_TYPE, STELLAR_TYPE, const char*);
}

// ----------------------------------------------------------------
// Instance-bound macros (use the calling object's id/type/stellar-type)
// ----------------------------------------------------------------

#define SHOW_ERROR(error, ...)                                                  \
    { ::error_helpers::show_error(error __VA_OPT__(,) __VA_ARGS__,              \
                                  ObjectId(), ObjectType(), StellarType(),      \
                                  __PRETTY_FUNCTION__); }

#define SHOW_ERROR_IF(cond, ...)                                                \
    { if (cond) SHOW_ERROR(__VA_ARGS__) }

#define THROW_ERROR(error, ...)                                                 \
    { ::error_helpers::throw_error(error __VA_OPT__(,) __VA_ARGS__,             \
                                   ObjectId(), ObjectType(), StellarType(),     \
                                   __PRETTY_FUNCTION__); }

#define THROW_ERROR_IF(cond, ...)                                               \
    { if (cond) THROW_ERROR(__VA_ARGS__) }

#define SHOW_WARN(error, ...)                                                   \
    { ::error_helpers::show_warn(error __VA_OPT__(,) __VA_ARGS__,               \
                                 ObjectId(), ObjectType(), StellarType(),       \
                                 __PRETTY_FUNCTION__); }

#define SHOW_WARN_IF(cond, ...)                                                 \
    { if (cond) SHOW_WARN(__VA_ARGS__) }


// ----------------------------------------------------------------
// Static-context macros (for use outside member functions: main,
// free functions, static methods - any context where ObjectId() and
// friends are not in scope). Hardcode -1L / NONE / NONE for context.
// ----------------------------------------------------------------

#define SHOW_ERROR_STATIC(error, ...)                                           \
    { ::error_helpers::show_error(error __VA_OPT__(,) __VA_ARGS__,              \
                                  -1L, OBJECT_TYPE::NONE, STELLAR_TYPE::NONE,   \
                                  __PRETTY_FUNCTION__); }

#define SHOW_ERROR_IF_STATIC(cond, ...)                                         \
    { if (cond) SHOW_ERROR_STATIC(__VA_ARGS__) }

#define THROW_ERROR_STATIC(error, ...)                                          \
    { ::error_helpers::throw_error(error __VA_OPT__(,) __VA_ARGS__,             \
                                   -1L, OBJECT_TYPE::NONE, STELLAR_TYPE::NONE,  \
                                   __PRETTY_FUNCTION__); }

#define THROW_ERROR_IF_STATIC(cond, ...)                                        \
    { if (cond) THROW_ERROR_STATIC(__VA_ARGS__) }

#define SHOW_WARN_STATIC(error, ...)                                            \
    { ::error_helpers::show_warn(error __VA_OPT__(,) __VA_ARGS__,               \
                                 -1L, OBJECT_TYPE::NONE, STELLAR_TYPE::NONE,    \
                                 __PRETTY_FUNCTION__); }

#define SHOW_WARN_IF_STATIC(cond, ...)                                          \
    { if (cond) SHOW_WARN_STATIC(__VA_ARGS__) }



