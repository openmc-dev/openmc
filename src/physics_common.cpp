#include "openmc/physics_common.h"

#include "openmc/random_lcg.h"
#include "openmc/settings.h"

namespace openmc {

//==============================================================================
// RUSSIAN_ROULETTE
//==============================================================================

void russian_roulette(Particle& p, double weight_survive)
{
  // if survival normalization is applicable, use normalized weight survive
  if((settings::source_file || settings::surf_source_read)&&(settings::survival_normalization && settings::survival_biasing)){
    weight_survive = p.wgt_survive();
  }
  if (weight_survive * prn(p.current_seed()) < p.wgt()) {
    p.wgt() = weight_survive;
  } else {
    p.wgt() = 0.;
  } 
}

} // namespace openmc
