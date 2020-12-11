#include "openmc/physics_common.h"

#include "openmc/settings.h"
#include "openmc/random_lcg.h"

namespace openmc {

//==============================================================================
// RUSSIAN_ROULETTE
//==============================================================================

HD void russian_roulette(Particle& p)
{
#ifdef __CUDA_ARCH__
  using gpu::weight_cutoff;
  using gpu::weight_survive;
#else
  using settings::weight_cutoff;
  using settings::weight_survive;
#endif
  if (p.wgt() < weight_cutoff) {
    if (prn(p.current_seed()) < p.wgt() / weight_survive) {
      p.wgt() = weight_survive;
      p.wgt_last() = p.wgt();
    } else {
      p.wgt() = 0.;
      p.wgt_last() = 0.;
      p.alive() = false;
    }
  }
}

} //namespace openmc
