#include "openmc/physics_common.h"

#include "openmc/settings.h"
#include "openmc/random_lcg.h"

namespace openmc {

//==============================================================================
// RUSSIAN_ROULETTE
//==============================================================================

void russian_roulette(Particle* p)
{
  if (p->wgt < settings::weight_cutoff) {
    if (prn() < p->wgt / settings::weight_survive) {
      p->wgt = settings::weight_survive;
      p->last_wgt = p->wgt;
    } else {
      p->wgt = 0.;
      p->last_wgt = 0.;
      p->alive = false;
    }
  }
}

} //namespace openmc
