#ifndef OPENMC_RANDOM_LCG_H
#define OPENMC_RANDOM_LCG_H

#include <cstdint>


namespace openmc {

//==============================================================================
// Module constants.
//==============================================================================

constexpr int N_STREAMS = 6;
extern "C" const int STREAM_TRACKING;
extern "C" const int STREAM_TALLIES;
extern "C" const int STREAM_SOURCE;
extern "C" const int STREAM_URR_PTABLE;
extern "C" const int STREAM_VOLUME;
extern "C" const int STREAM_PHOTON;
constexpr int64_t DEFAULT_SEED = 1;

//==============================================================================
//! Generate a pseudo-random number using a linear congruential generator.
//! @param prn_seeds Pseudorandom number seed array
//! @param stream Pseudorandom number stream index
//! @return A random number between 0 and 1
//==============================================================================

extern "C" double prn(uint64_t * seeds, int stream);

//==============================================================================
//! Generate a random number which is 'n' times ahead from the current seed.
//!
//! The result of this function will be the same as the result from calling
//! `prn()` 'n' times.
//! @param n The number of RNG seeds to skip ahead by
//! @param prn_seeds Pseudorandom number seed array
//! @param stream Pseudorandom number stream index
//! @return A random number between 0 and 1
//==============================================================================

extern "C" double future_prn(int64_t n, uint64_t * prn_seeds, int stream);

//==============================================================================
//! Set the RNG seeds to unique values based on the ID of the particle.
//! @param prn_seeds Pseudorandom number seed array
//! @param id The particle ID
//==============================================================================

extern "C" void set_particle_seed(int64_t id, uint64_t * prn_seeds );

//==============================================================================
//! Advance the random number seed 'n' times from the current seed.
//! @param prn_seeds Pseudorandom number seed array
//! @param stream Pseudorandom number stream index
//! @param n The number of RNG seeds to skip ahead by
//==============================================================================

extern "C" void advance_prn_seed(int64_t n, uint64_t * prn_seeds, int stream);

//==============================================================================
//! Advance a random number seed 'n' times.
//!
//! This is usually used to skip a fixed number of random numbers (the stride)
//! so that a given particle always has the same starting seed regardless of
//! how many processors are used.
//! @param n The number of RNG seeds to skip ahead by
//! @param seed The starting to seed to advance from
//==============================================================================

uint64_t future_seed(uint64_t n, uint64_t seed);

//==============================================================================
//                               API FUNCTIONS
//==============================================================================

//==============================================================================
//! Get OpenMC's master seed.
//==============================================================================

extern "C" int64_t openmc_get_seed();

//==============================================================================
//! Set OpenMC's master seed.
//! @param new_seed The master seed. All other seeds will be derived from this
//! one.
//==============================================================================

extern "C" void openmc_set_seed(int64_t new_seed);

} // namespace openmc
#endif // OPENMC_RANDOM_LCG_H
