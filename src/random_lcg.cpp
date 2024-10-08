#include "openmc/random_lcg.h"

#include <cmath>

namespace openmc {

// Starting seed
int64_t master_seed {1};

// LCG parameters
constexpr uint64_t prn_mult {6364136223846793005ULL}; // multiplication
constexpr uint64_t prn_add {1442695040888963407ULL};  // additive factor, c
constexpr uint64_t prn_stride {152917LL}; // stride between particles

//==============================================================================
// PRN
//==============================================================================

// 64 bit implementation of the PCG-RXS-M-XS 64-bit state / 64-bit output
// geneator Adapted from: https://github.com/imneme/pcg-c, in particular
// https://github.com/imneme/pcg-c/blob/83252d9c23df9c82ecb42210afed61a7b42402d7/include/pcg_variants.h#L188-L192
// @techreport{oneill:pcg2014,
//    title = "PCG: A Family of Simple Fast Space-Efficient Statistically Good
//    Algorithms for Random Number Generation", author = "Melissa E. O'Neill",
//    institution = "Harvey Mudd College",
//    address = "Claremont, CA",
//    number = "HMC-CS-2014-0905",
//    year = "2014",
//    month = Sep,
//    xurl = "https://www.cs.hmc.edu/tr/hmc-cs-2014-0905.pdf",
//}
double prn(uint64_t* seed)
{
  // Advance the LCG
  *seed = (prn_mult * (*seed) + prn_add);

  // Permute the output
  uint64_t word =
    ((*seed >> ((*seed >> 59u) + 5u)) ^ *seed) * 12605985483714917081ull;
  uint64_t result = (word >> 43u) ^ word;

  // Convert output from unsigned integer to double
  return ldexp(result, -64);
}

//==============================================================================
// FUTURE_PRN
//==============================================================================

double future_prn(int64_t n, uint64_t seed)
{
  uint64_t fseed = future_seed(static_cast<uint64_t>(n), seed);
  return prn(&fseed);
}

//==============================================================================
// INIT_SEED
//==============================================================================

uint64_t init_seed(int64_t id, int offset)
{
  return future_seed(
    static_cast<uint64_t>(id) * prn_stride, master_seed + offset);
}

//==============================================================================
// INIT_PARTICLE_SEEDS
//==============================================================================

void init_particle_seeds(int64_t id, uint64_t* seeds)
{
  for (int i = 0; i < N_STREAMS; i++) {
    seeds[i] =
      future_seed(static_cast<uint64_t>(id) * prn_stride, master_seed + i);
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
  // The algorithm here to determine the parameters used to skip ahead is
  // described in F. Brown, "Random Number Generation with Arbitrary Stride,"
  // Trans. Am. Nucl. Soc. (Nov. 1994). This algorithm is able to skip ahead in
  // O(log2(N)) operations instead of O(N). Basically, it computes parameters G
  // and C which can then be used to find x_N = G*x_0 + C mod 2^M.

  // Initialize constants
  uint64_t g {prn_mult};
  uint64_t c {prn_add};
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
  return g_new * seed + c_new;
}

//==============================================================================
//                               API FUNCTIONS
//==============================================================================

extern "C" int64_t openmc_get_seed()
{
  return master_seed;
}

extern "C" void openmc_set_seed(int64_t new_seed)
{
  master_seed = new_seed;
}

} // namespace openmc
