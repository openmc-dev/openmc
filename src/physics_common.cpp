#include "openmc/physics_common.h"

#include "openmc/settings.h"
#include "openmc/random_lcg.h"

namespace openmc {

//==============================================================================
// RUSSIAN_ROULETTE
//==============================================================================

void russian_roulette(Particle& p, double weight_survive)
{
  if (prn(p.current_seed()) < p.wgt_ / weight_survive) {
    p.wgt_ = weight_survive;
  } else {
    p.wgt_ = 0.;
    p.alive_ = false;
  }
}

} //namespace openmc
