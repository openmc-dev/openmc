#include "random_lcg.h"
#include <cmath>


// Starting seed
int64_t seed = 1;

// LCG parameters
const int64_t prn_mult   = 2806196910506780709LL;   // multiplication factor, g
const int64_t prn_add    = 1;                       // additive factor, c
const int64_t prn_mod    = -0x8000000000000000;     // -2^63
const int64_t prn_mask   = 0x7fffffffffffffff;      // 2^63 - 1
const int64_t prn_stride = 152917LL;                // stride between particles
const double  prn_norm = pow(2, -63);               // 2^-63

// Module constants
const int N_STREAMS = 5;
const int STREAM_TRACKING   = 0;
const int STREAM_TALLIES    = 1;
const int STREAM_SOURCE     = 2;
const int STREAM_URR_PTABLE = 3;
const int STREAM_VOLUME     = 4;

// Current PRNG state
int64_t prn_seed[N_STREAMS];  // current seed
int     stream;               // current RNG stream
#pragma omp threadprivate(prn_seed, stream)


//==============================================================================
// PRN
//==============================================================================

extern "C" double
prn()
{
  // This algorithm uses bit-masking to find the next integer(8) value to be
  // used to calculate the random number.
  prn_seed[stream] = (prn_mult*prn_seed[stream] + prn_add) & prn_mask;

  // Once the integer is calculated, we just need to divide by 2**m,
  // represented here as multiplying by a pre-calculated factor
  double pseudo_rn = prn_seed[stream] * prn_norm;
  return pseudo_rn;
}

//==============================================================================
// FUTURE_PRN
//==============================================================================

extern "C" double
future_prn(int64_t n)
{
  double pseudo_rn = future_seed(n, prn_seed[stream]) * prn_norm;
  return pseudo_rn;
}

//==============================================================================
// SET_PARTICLE_SEED
//==============================================================================

extern "C" void
set_particle_seed(int64_t id)
{
  for (int i = 0; i < N_STREAMS; i++) {
    prn_seed[i] = future_seed(id * prn_stride, seed + i);
  }
}

//==============================================================================
// ADVANCE_PRN_SEED
//==============================================================================

extern "C" void
advance_prn_seed(int64_t n)
{
  prn_seed[stream] = future_seed(n, prn_seed[stream]);
}

//==============================================================================
// FUTURE_SEED
//==============================================================================

extern "C" int64_t
future_seed(int64_t n, int64_t seed)
{
  // In cases where we want to skip backwards, we add the period of the random
  // number generator until the number of PRNs to skip is positive since
  // skipping ahead that much is the same as skipping backwards by the original
  // amount.

  int64_t nskip = n;
  while (nskip < 0) nskip += prn_mod;

  // Make sure nskip is less than 2^M.
  nskip &= prn_mask;

  // The algorithm here to determine the parameters used to skip ahead is
  // described in F. Brown, "Random Number Generation with Arbitrary Stride,"
  // Trans. Am. Nucl. Soc. (Nov. 1994). This algorithm is able to skip ahead in
  // O(log2(N)) operations instead of O(N). Basically, it computes parameters G
  // and C which can then be used to find x_N = G*x_0 + C mod 2^M.

  // Initialize constants
  int64_t g     = prn_mult;
  int64_t c     = prn_add;
  int64_t g_new = 1;
  int64_t c_new = 0;

  while (nskip > 0) {
    // Check if the least significant bit is 1.
    if (nskip & 1) {
      g_new = (g_new * g) & prn_mask;
      c_new = (c_new * g + c) & prn_mask;
    }
    c = ((g + 1) * c) & prn_mask;
    g = (g * g) & prn_mask;

    // Move bits right, dropping least significant bit.
    nskip >>= 1;
  }

  // With G and C, we can now find the new seed.
  int64_t new_seed = (g_new * seed + c_new) & prn_mask;
  return new_seed;
}

//==============================================================================
// PRN_SET_STREAM
//==============================================================================

extern "C" void
prn_set_stream(int i)
{
  stream = i - 1;
}

//==============================================================================
//                               API FUNCTIONS
//==============================================================================

extern "C" int
openmc_set_seed(int64_t new_seed)
{
  seed = new_seed;
#pragma omp parallel
  for (int i = 0; i < N_STREAMS; i++) {
    prn_seed[i] = seed + i;
  }
  stream = STREAM_TRACKING;
#pragma end omp parallel
  return 0;
}
