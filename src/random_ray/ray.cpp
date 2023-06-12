#include "openmc/random_ray/ray.h"

namespace openmc {

  Ray::Ray()
  {
    angular_flux_.resize(data::mg.num_energy_groups_);
  }
  
  void event_advance_ray(double distance_inactive, double distance_active)
  {
  }
  
  void attenuate_flux(double distance, bool is_active)
  {
  }

} // namespace openmc
