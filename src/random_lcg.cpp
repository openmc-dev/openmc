#include "openmc/random_lcg.h"

#include <cmath>


namespace openmc {

// Starting seed
int64_t master_seed {1};

// LCG parameters
constexpr uint64_t prn_mult   {2806196910506780709LL};   // multiplication
                                                         //   factor, g
constexpr uint64_t prn_add    {1};                       // additive factor, c
constexpr uint64_t prn_mod    {0x8000000000000000};      // 2^63
constexpr uint64_t prn_mask   {0x7fffffffffffffff};      // 2^63 - 1
constexpr uint64_t prn_stride {152917LL};                // stride between
                                                         //   particles
constexpr double   prn_norm   {1.0 / prn_mod};           // 2^-63

//==============================================================================
// PRN
//==============================================================================

double prn(uint64_t* seed)
{
  // This algorithm uses bit-masking to find the next integer(8) value to be
  // used to calculate the random number.
  *seed = (prn_mult * (*seed) + prn_add) & prn_mask;

  // Once the integer is calculated, we just need to divide by 2**m,
  // represented here as multiplying by a pre-calculated factor
  return (*seed) * prn_norm;
}

//==============================================================================
// FUTURE_PRN
//==============================================================================

double future_prn(int64_t n, uint64_t seed)
{
  return future_seed(static_cast<uint64_t>(n), seed) * prn_norm;
}

//==============================================================================
// INIT_SEED
//==============================================================================

uint64_t init_seed(int64_t id, int offset)
{
  return future_seed(static_cast<uint64_t>(id) * prn_stride, master_seed + offset);
}

//==============================================================================
// INIT_PARTICLE_SEEDS
//==============================================================================

void init_particle_seeds(int64_t id, uint64_t* seeds)
{
  for (int i = 0; i < N_STREAMS; i++) {
    seeds[i] = future_seed(static_cast<uint64_t>(id) * prn_stride, master_seed + i);
  }
}

//==============================================================================
// ADVANCE_PRN_SEED
//==============================================================================

void advance_prn_seed(int64_t n, uint64_t* seed)
{
  *seed = future_seed(static_cast<uint64_t>(n), *seed);
}

//==============================================================================
// FUTURE_SEED
//==============================================================================

uint64_t future_seed(uint64_t n, uint64_t seed)
{
  // Make sure nskip is less than 2^M.
  n &= prn_mask;

  // The algorithm here to determine the parameters used to skip ahead is
  // described in F. Brown, "Random Number Generation with Arbitrary Stride,"
  // Trans. Am. Nucl. Soc. (Nov. 1994). This algorithm is able to skip ahead in
  // O(log2(N)) operations instead of O(N). Basically, it computes parameters G
  // and C which can then be used to find x_N = G*x_0 + C mod 2^M.

  // Initialize constants
  uint64_t g     {prn_mult};
  uint64_t c     {prn_add};
  uint64_t g_new {1};
  uint64_t c_new {0};

  while (n > 0) {
    // Check if the least significant bit is 1.
    if (n & 1) {
      g_new *= g;
      c_new = c_new * g + c;
    }
    c *= (g + 1);
    g *= g;

    // Move bits right, dropping least significant bit.
    n >>= 1;
  }

  // With G and C, we can now find the new seed.
  return (g_new * seed + c_new) & prn_mask;
}

//==============================================================================
//                               API FUNCTIONS
//==============================================================================

extern "C" int64_t openmc_get_seed() {return master_seed;}

extern "C" void openmc_set_seed(int64_t new_seed)
{
  master_seed = new_seed;
}

} // namespace openmc
