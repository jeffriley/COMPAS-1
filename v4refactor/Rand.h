#pragma once

#define RAND Rand::Instance()

#include <algorithm>
#include <cstdlib>
#include <ctime>

#include <gsl/gsl_rng.h>                                    // GSL random number generator
#include <gsl/gsl_randist.h>

#include "Errors.h"


/*
 * Rand Singleton - interface to GSL rng (random number generator)
 *
 * Singletons and global variables are sometimes frowned-upon, but doing it this
 * way means the objects don't need to be passed around to all and sundry.
 * I think convenience and clarity sometimes trump dogma.
 *
 * This is a header-only implementation.
 */

class Rand {

private:

    Rand() = default;
    Rand(Rand const&)            = delete;
    Rand& operator=(Rand const&) = delete;

    // Meyers-static singleton state flag - inline static so the definition
    // lives in this (header-only) class without a separate .cpp file.
    // See IsInstantiated() below.
    inline static bool m_Instantiated = false;

    gsl_rng* m_Rng  = nullptr; // GSL random number generator
    ULongT   m_Seed = 0;


public:

    ~Rand() {
        m_Instantiated = false;
        Free();
    }

    static Rand* Instance() {
        static Rand self;
        m_Instantiated = true; // Set after `self` is fully constructed
        return &self;
    }

    static bool IsInstantiated() { return m_Instantiated; }


    /*
     * Initialise
     *
     * @brief
     * Initialise the random number generator.
     *
     * void Initialise()
     */
    void Initialise() {

        if (m_Rng) return; // Already initialised - just return silently

        // Set up the gsl random number generator
        gsl_rng_env_setup();

        // Seed the random number generator.
        // Preferably use environment GSL_RNG_SEED, otherwise time(nullptr)
        if (!std::getenv("GSL_RNG_SEED")) {
            gsl_rng_default_seed = std::time(nullptr);
        }

        m_Seed = gsl_rng_default_seed;

        m_Rng = gsl_rng_alloc(gsl_rng_default);
        if (!m_Rng) THROW_ERROR_STATIC(ERROR::NO_RNG);
    }


    /*
     * Free
     *
     * @brief
     * Free the dynamically allocated memory.
     *
     * Idempotent: safe to call multiple times. After freeing, m_Rng is set to
     * nullptr so a subsequent call (e.g. from ~Rand() after main has already
     * invoked Free()) is a no-op rather than a double-free.
     *
     * void Free()
     */
    void Free() {
        if (m_Rng) {
            gsl_rng_free(m_Rng);
            m_Rng = nullptr;
        }
    }


    GNU_PURE ULongT CurrentSeed() const { return m_Seed; }
    GNU_PURE ULongT DefaultSeed() const { return gsl_rng_default_seed; }
    ULongT Seed(const ULongT p_Seed) { 
        if (!m_Rng) THROW_ERROR_STATIC(ERROR::NO_RNG);
        gsl_rng_set(m_Rng, p_Seed); m_Seed = p_Seed;
        return p_Seed;
    }


    /*
     * Random
     *
     * @brief
     * Return a random floating point number uniformly distributed in the range [0.0, 1.0)
     *
     * double Random()
     *
     * @return                                      Random floating point number uniformly distributed in the range [0.0, 1.0)
     */
    double Random() {
        if (!m_Rng) THROW_ERROR_STATIC(ERROR::NO_RNG);
        return gsl_rng_uniform(m_Rng);
    }


    /*
     * Random
     *
     * @brief
     * Return a random floating point number uniformly distributed in the range [p_Lower, p_Upper),
     * where p_Lower <= p_Upper (p_Lower and p_Upper will be swapped if p_Lower > p_Upper as passed)
     *
     * double Random(const double p_Lower, const double p_Upper)
     *
     * @param       p_Lower                         Inclusive lower bound of range of distribution
     * @param       p_Upper                         Exclusive upper bound of range of distribution
     * @return                                      Random floating point number uniformly distributed in the range [p_Lower, p_Upper)
     */
    double Random(const double p_Lower, const double p_Upper) {
        if (!m_Rng) THROW_ERROR_STATIC(ERROR::NO_RNG);
        const double lower = std::min(p_Lower, p_Upper);
        const double upper = std::max(p_Lower, p_Upper);
        return (gsl_rng_uniform(m_Rng) * (upper - lower)) + lower;
    }


    /*
     * RandomInt
     *
     * @brief
     * Return a random integer number uniformly distributed in the range [p_Lower, p_Upper),
     * where p_Lower <= p_Upper (p_Lower and p_Upper will be swapped if p_Lower > p_Upper as passed)
     *
     * int RandomInt(const int p_Lower, const int p_Upper)
     *
     * @param       p_Lower                         Inclusive lower bound of range of distribution
     * @param       p_Upper                         Exclusive upper bound of range of distribution
     * @return                                      Random integer number uniformly distributed in the range [p_Lower, p_Upper)
     *                                              (Will be p_Lower if p_Lower == p_Upper)
     */
    int RandomInt(const int p_Lower, const int p_Upper) {
        if (!m_Rng) THROW_ERROR_STATIC(ERROR::NO_RNG);
        const int lower = std::min(p_Lower, p_Upper);
        const int upper = std::max(p_Lower, p_Upper);
        return (upper == lower) ? lower : gsl_rng_uniform_int(m_Rng, (upper - lower)) + lower;
    }
    int RandomInt(const int p_Upper) { return p_Upper < 0 ? 0 : RandomInt(0, p_Upper); }


    /*
     * RandomGaussian
     *
     * @brief
     * Return a Gaussian random variate, with mean zero and standard deviation p_Sigma
     *
     * double RandomGaussian(const double p_Sigma)
     *
     * @param       p_Sigma                         Standard deviation of distribution
     * @return                                      Gaussian random variate, with mean zero and standard deviation p_Sigma
     */
    double RandomGaussian(const double p_Sigma) {
        if (!m_Rng) THROW_ERROR_STATIC(ERROR::NO_RNG);
        return gsl_ran_gaussian(m_Rng, p_Sigma);
    }
};
