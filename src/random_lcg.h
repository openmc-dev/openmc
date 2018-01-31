#ifndef RANDOM_LCG_H
#define RANDOM_LCG_H

#include <cstdint>


namespace openmc {

//==============================================================================
// Module constants.
//==============================================================================

extern "C" const int N_STREAMS;
extern "C" const int STREAM_TRACKING;
extern "C" const int STREAM_TALLIES;
extern "C" const int STREAM_SOURCE;
extern "C" const int STREAM_URR_PTABLE;
extern "C" const int STREAM_VOLUME;

//==============================================================================
//! Generate a pseudo-random number using a linear congruential generator.
//! @return A random number between 0 and 1
//==============================================================================

extern "C" double prn();

//==============================================================================
//! Generate a random number which is 'n' times ahead from the current seed.
//!
//! The result of this function will be the same as the result from calling
//! `prn()` 'n' times.
//! @param n The number of RNG seeds to skip ahead by
//! @return A random number between 0 and 1
//==============================================================================

extern "C" double future_prn(int64_t n);

//==============================================================================
//! Set the RNG seed to a unique value based on the ID of the particle.
//! @param id The particle ID
//==============================================================================

extern "C" void set_particle_seed(int64_t id);

//==============================================================================
//! Advance the random number seed 'n' times from the current seed.
//! @param n The number of RNG seeds to skip ahead by
//==============================================================================

extern "C" void advance_prn_seed(int64_t n);

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
//! Switch the RNG to a different stream of random numbers.
//!
//! If random numbers are needed in routines not used directly for tracking
//! (e.g. physics), this allows the numbers to be generated without affecting
//! reproducibility of the physics.
//! @param n The RNG stream to switch to. Use the constants such as
//! `STREAM_TRACKING` and `STREAM_TALLIES` for this argument.
//==============================================================================

extern "C" void prn_set_stream(int n);

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
#endif // RANDOM_LCG_H
