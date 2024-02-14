#ifndef OPENMC_RANDOM_RAY_RAY_H
#define OPENMC_RANDOM_RAY_RAY_H

#include "openmc/particle.h"
#include "openmc/random_ray/flat_source_domain.h"

namespace openmc {

/*
 * The RandomRay class encompasses data and methods for transporting random rays
 * through the model. It is a small extension of the Particle class.
 */

class RandomRay : public Particle {
public:
  //==========================================================================
  // Constructors
  RandomRay();
  RandomRay(uint64_t ray_id, int sampling_source, FlatSourceDomain* domain);

  //==========================================================================
  // Methods
  void event_advance_ray();
  void attenuate_flux(double distance, bool is_active);
  void initialize_ray(uint64_t ray_id, int sampling_source, FlatSourceDomain* domain);
  uint64_t transport_history_based_single_ray();

  //==========================================================================
  // Data

  std::vector<float> angular_flux_;
  std::vector<float> delta_psi_;
  double distance_travelled_ {0};
  int negroups_;
  bool is_active_ {false};
  bool is_alive_ {true};
  FlatSourceDomain* domain_ {nullptr};
}; // class RandomRay

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_RAY_H
