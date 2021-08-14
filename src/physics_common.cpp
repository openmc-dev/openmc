#include "openmc/physics_common.h"

#include "openmc/random_lcg.h"
#include "openmc/settings.h"

namespace openmc {

//==============================================================================
// RUSSIAN_ROULETTE
//==============================================================================

void russian_roulette(Particle& p)
{
  if (p.wgt() < settings::weight_cutoff) {
    if (prn(p.current_seed()) < p.wgt() / settings::weight_survive) {
      p.wgt() = settings::weight_survive;
      p.wgt_last() = p.wgt();
    } else {
      p.wgt() = 0.;
      p.wgt_last() = 0.;
      p.alive() = false;
    }
  }
}

} // namespace openmc
