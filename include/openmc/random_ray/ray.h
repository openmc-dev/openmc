#ifndef OPENMC_RANDOM_RAY_RAY_H
#define OPENMC_RANDOM_RAY_RAY_H

#include "openmc/particle.h"

namespace openmc {

/*
 * The Ray class encompasses data and methods for transporting random rays
 * through the model. It is a small extension of the Particle class.
 */

class Ray : public Particle {
public:
  //==========================================================================
  // Constructors
  Ray();

  //==========================================================================
  // Methods
  void event_advance_ray(double distance_inactive, double distance_active);
  void attenuate_flux(double distance, bool is_active);
  void initialize_ray(uint64_t index_source, uint64_t nrays, int iter);
  uint64_t transport_history_based_single_ray(openmc::Particle& p, double distance_inactive, double distance_active)

  //==========================================================================
  // Data
  
  std::vector<float> angular_flux_;
  double distance_travelled_ {0};
  bool is_active_ {false};

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_RAY_H
