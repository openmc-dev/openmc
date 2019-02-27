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
    if (prn() < p->wgt_ / settings::weight_survive) {
      p->wgt_ = settings::weight_survive;
      p->last_wgt_ = p->wgt_;
    } else {
      p->wgt_ = 0.;
      p->last_wgt_ = 0.;
      p->alive_ = false;
    }
  }
}

} //namespace openmc
