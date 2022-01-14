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
    if (settings::weight_survive * prn(p.current_seed()) < p.wgt()) {
      p.wgt() = settings::weight_survive;
    } else {
      p.wgt() = 0.;
      p.alive() = false;
    }
  }
}

} // namespace openmc
