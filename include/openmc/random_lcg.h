#ifndef OPENMC_RANDOM_LCG_H
#define OPENMC_RANDOM_LCG_H

#include <cstdint>


namespace openmc {

//==============================================================================
// Module constants.
//==============================================================================

constexpr int N_STREAMS         {6};
constexpr int STREAM_TRACKING   {0};
constexpr int STREAM_TALLIES    {1};
constexpr int STREAM_SOURCE     {2};
constexpr int STREAM_URR_PTABLE {3};
constexpr int STREAM_VOLUME     {4};
constexpr int STREAM_PHOTON     {5};
constexpr int64_t DEFAULT_SEED  {1};

//==============================================================================
//! Generate a pseudo-random number using a linear congruential generator.
//! @param seed Pseudorandom number seed pointer
//! @return A random number between 0 and 1
//==============================================================================

double prn(uint64_t* seed);

//==============================================================================
//! Generate a random number which is 'n' times ahead from the current seed.
//!
//! The result of this function will be the same as the result from calling
//! `prn()` 'n' times, though without the side effect of altering the RNG
//! state.
//! @param n The number of RNG seeds to skip ahead by
//! @param seed Pseudorandom number seed
//! @return A random number between 0 and 1
//==============================================================================

double future_prn(int64_t n, uint64_t seed);

//==============================================================================
//! Set a RNG seed to a unique value based on a unique particle ID by striding
//! the seed.
//! @param id The particle ID
//! @param offset The offset from the master seed to be used (e.g., for creating
//! different streams)
//! @return The initialized seed value
//==============================================================================

uint64_t init_seed(int64_t id, int offset);

//==============================================================================
//! Set the RNG seeds to unique values based on the ID of the particle. This
//! function initializes the seeds for all RNG streams of the particle via
//! striding.
//! @param seeds Pseudorandom number seed array
//! @param id The particle ID
//==============================================================================

void init_particle_seeds(int64_t id, uint64_t* seeds);

//==============================================================================
//! Advance the random number seed 'n' times from the current seed. This
//! differs from the future_prn() function in that this function does alter
//! the RNG state.
//! @param seed Pseudorandom number seed pointer
//! @param n The number of RNG seeds to skip ahead by
//==============================================================================

void advance_prn_seed(int64_t n, uint64_t* seed);

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
