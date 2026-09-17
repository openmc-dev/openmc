#include "openmc/physics_common.h"

#include "openmc/random_lcg.h"
#include "openmc/settings.h"

namespace openmc {

//==============================================================================
// RUSSIAN_ROULETTE
//==============================================================================

void russian_roulette(Particle& p, double weight_survive)
{
  if (weight_survive * prn(p.current_seed()) < p.wgt()) {
    p.wgt() = weight_survive;
  } else {
    p.wgt() = 0.;
  }
}

void apply_russian_roulette(Particle& p)
{
  // Exit if survival biasing is turned off
  if (!settings::survival_biasing)
    return;

  // if survival normalization is on, use normalized weight cutoff and
  // normalized weight survive
  if (settings::survival_normalization) {
    if (p.wgt() < settings::weight_cutoff * p.wgt_born()) {
      russian_roulette(p, settings::weight_survive * p.wgt_born());
    }
  } else if (p.wgt() < settings::weight_cutoff) {
    russian_roulette(p, settings::weight_survive);
  }
}
} // namespace openmc
