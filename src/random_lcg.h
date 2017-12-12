#ifndef RANDOM_LCG_H
#define RANDOM_LCG_H

#include <cstdint>

//==============================================================================
// Module constants.
//==============================================================================

extern "C" const int N_STREAMS         = 5;
extern "C" const int STREAM_TRACKING   = 0;
extern "C" const int STREAM_TALLIES    = 1;
extern "C" const int STREAM_SOURCE     = 2;
extern "C" const int STREAM_URR_PTABLE = 3;
extern "C" const int STREAM_VOLUME     = 4;

//==============================================================================
// PRN generates a pseudo-random number using a linear congruential generator.
//==============================================================================

extern "C" double prn();

//==============================================================================
// FUTURE_PRN generates a pseudo-random number which is 'n' times ahead from the
// current seed.
//==============================================================================

extern "C" double future_prn(int64_t n);

//==============================================================================
// SET_PARTICLE_SEED sets the seed to a unique value based on the ID of the
// particle.
//==============================================================================

extern "C" void set_particle_seed(int64_t id);

//==============================================================================
// ADVANCE_PRN_SEED advances the random number seed 'n' times from the current
// seed.
//==============================================================================

extern "C" void advance_prn_seed(int64_t n);

//==============================================================================
// FUTURE_SEED advances the random number seed 'skip' times. This is usually
// used to skip a fixed number of random numbers (the stride) so that a given
// particle always has the same starting seed regardless of how many processors
// are used.
//==============================================================================

uint64_t future_seed(uint64_t n, uint64_t seed);

//==============================================================================
// PRN_SET_STREAM changes the random number stream. If random numbers are needed
// in routines not used directly for tracking (e.g. physics), this allows the
// numbers to be generated without affecting reproducibility of the physics.
//==============================================================================

extern "C" void prn_set_stream(int n);

//==============================================================================
//                               API FUNCTIONS
//==============================================================================

extern "C" int openmc_set_seed(int64_t new_seed);

#endif // RANDOM_LCG_H
