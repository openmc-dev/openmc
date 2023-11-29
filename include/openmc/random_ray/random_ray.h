#ifndef OPENMC_RANDOM_RAY_RAY_H
#define OPENMC_RANDOM_RAY_RAY_H

#include "openmc/particle.h"

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
  RandomRay(uint64_t index_source);

  //==========================================================================
  // Methods
  void event_advance_ray();
  void attenuate_flux(double distance, bool is_active);
  void initialize_ray(uint64_t index_source);
  uint64_t transport_history_based_single_ray();

  //==========================================================================
  // Data

  std::vector<float> angular_flux_;
  std::vector<float> delta_psi_;
  double distance_travelled_ {0};
  bool is_active_ {false};
  bool is_alive_ {true};
}; // class RandomRay

//==============================================================================
// Non-member functions
//==============================================================================
inline float cjosey_exponential(const float tau);

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_RAY_H
