#include "openmc/physics_common.h"

#include "openmc/settings.h"
#include "openmc/random_lcg.h"

namespace openmc {

//==============================================================================
// RUSSIAN_ROULETTE
//==============================================================================

void russian_roulette(Particle* p)
{
  if (p->wgt_ < settings::weight_cutoff) {
    if (prn(p->current_seed()) < p->wgt_ / settings::weight_survive) {
      p->wgt_ = settings::weight_survive;
      p->wgt_last_ = p->wgt_;
    } else {
      p->wgt_ = 0.;
      p->wgt_last_ = 0.;
      p->alive_ = false;
    }
  }
}

} //namespace openmc
