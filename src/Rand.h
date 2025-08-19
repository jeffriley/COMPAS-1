#ifndef __Rand_H__
#define __Rand_H__

#define RAND Rand::Instance()

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <algorithm>
#include <time.h>
#include <iostream>


/*
 * Rand Singleton - interface to GSL rng (random number generator)
 *
 * Singletons and global variables are sometimes frowned-upon, but doing it this
 * way means the objects don't need to be passed around to all and sundry.
 * I think convenience and clarity sometimes trump dogma.
 */

class Rand {

private:

    Rand() {
        // create instance of gsl random number generator
        m_Rng = gsl_rng_alloc(gsl_rng_default);     
        if (m_Rng) {                                    // ok?
                                                        // yes
            // set up the gsl random number generator
            gsl_rng_env_setup();

            // seed the random number generator
            // use environment GSL_RNG_SEED if set, otherwise time(0)
            if (!getenv("GSL_RNG_SEED")) gsl_rng_default_seed = time(NULL);
            m_Seed = gsl_rng_default_seed;
        }
    };

    Rand(Rand const&) = delete;
    Rand& operator = (Rand const&) = delete;

    unsigned long int m_Seed;
    gsl_rng*          m_Rng;    // GSL random number generator

    static Rand*      m_Instance;

public:

    static Rand* Instance() {
        if (!m_Instance) m_Instance = new Rand();
        return m_Instance;
    }

    void                Free()                              { gsl_rng_free(m_Rng); }                            // free dynamically allocated memory

    unsigned long int   CurrentSeed() const                 { return m_Seed; }                                  // return current seed value
    unsigned long int   DefaultSeed() const                 { return gsl_rng_default_seed; }                    // return default seed value
    void                Seed(const unsigned long p_Seed)    { gsl_rng_set(m_Rng, p_Seed); m_Seed = p_Seed; }    // set current seed value

    // see implementations below
    double              Random() const;
    double              Random(const double p_Lower, const double p_Upper) const;
    int                 RandomInt(const int p_Lower, const int p_Upper) const;
    int                 RandomInt(const int p_Upper) const;
    double              RandomGaussian(const double p_Sigma) const;
};


// initialise instance pointer
Rand* Rand::m_Instance = nullptr;


/*
 * Random
 * Return a random floating point number uniformly distributed in the range [0.0, 1.0)
 *
 *
 * double Random() const
 *
 * @return                                      Random floating point number uniformly distributed in the range [0.0, 1.0)
 */
inline double Rand::Random() const  { return gsl_rng_uniform(m_Rng); }


/*
 * Random
 * Return a random floating point number uniformly distributed in the range [p_Lower, p_Upper), where p_Lower <= p_Upper
 * (p_Lower and p_Upper will be swapped if p_Lower > p_Upper as passed)
 *
 *
 * double Random(const double p_Lower, const double p_Upper) const
 *
 * @param   [IN]    p_Lower                     Inclusive lower bound of range of distribution
 * @param   [IN]    p_Upper                     Exclusive upper bound of range of distribution
 * @return                                      Random floating point number uniformly distributed in the range [p_Lower, p_Upper)
 */
inline double Rand::Random(const double p_Lower, const double p_Upper) const {
    double lower = std::min(p_Lower, p_Upper);
    return (gsl_rng_uniform(m_Rng) * (std::max(p_Lower, p_Upper) - lower)) + lower;
}


/*
 * RandomInt
 *
 * Return a random integer number uniformly distributed in the range [p_Lower, p_Upper), where p_Lower <= p_Upper
 * (p_Lower and p_Upper will be swapped if p_Lower > p_Upper as passed)
 *
 *
 * int RandomInt(const int p_Lower, const int p_Upper) const
 *
 * @param   [IN]    p_Lower                     Inclusive lower bound of range of distribution
 * @param   [IN]    p_Upper                     Exclusive upper bound of range of distribution
 * @return                                      Random integer number uniformly distributed in the range [p_Lower, p_Upper)
 */
inline int Rand::RandomInt(const int p_Lower, const int p_Upper) const {
    int lower = std::min(p_Lower, p_Upper);
    return gsl_rng_uniform_int(m_Rng, (std::max(p_Lower, p_Upper) - lower)) + lower;
}


/*
 * RandomInt
 * Return a random integer number uniformly distributed in the range [0, p_Upper)
 *
 *
 * int RandomInt(const int p_Upper) const
 *
 * @param   [IN]    p_Upper                     Exclusive upper bound of range of distribution
 * @return                                      Random integer number uniformly distributed in the range [0, p_Upper)
 *                                              Will be 0 if p_Upper is < 0
 */
inline int Rand::RandomInt(const int p_Upper) const { return p_Upper < 0 ? 0 : RandomInt(0, p_Upper); }


/*
 * RandomGaussian
 * Return a Gaussian random variate, with mean zero and standard deviation p_Sigma
 *
 *
 * double RandomGaussian(const double p_Sigma) const
 *
 * @param   [IN]    p_Sigma                     Standard deviation of distribution
 * @return                                      Gaussian random variate, with mean zero and standard deviation p_Sigma
 */
inline double Rand::RandomGaussian(const double p_Sigma) const { return gsl_ran_gaussian(m_Rng, p_Sigma); }

#endif // __Rand_H__
