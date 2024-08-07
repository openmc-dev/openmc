#ifndef OPENMC_RANDOM_RAY_H
#define OPENMC_RANDOM_RAY_H

#include "openmc/memory.h"
#include "openmc/particle.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/moment_matrix.h"
#include "openmc/source.h"

namespace openmc {

/*
 * The RandomRay class encompasses data and methods for transporting random rays
 * through the model. It is a small extension of the Particle class.
 */

// TODO: Inherit from GeometryState instead of Particle
class RandomRay : public Particle {
public:
  //----------------------------------------------------------------------------
  // Constructors
  RandomRay();
  RandomRay(uint64_t ray_id, FlatSourceDomain* domain, bool uncollided_ray);

  //----------------------------------------------------------------------------
  // Methods
  void event_advance_ray();
  void attenuate_flux(double distance, bool is_active);
  void attenuate_flux_flat_source(double distance, bool is_active);
  void attenuate_flux_linear_source(double distance, bool is_active);
  void event_advance_ray_first_collided();
  

  void initialize_ray(uint64_t ray_id, FlatSourceDomain* domain,bool uncollided_ray);
  uint64_t transport_history_based_single_ray();
  
  //----------------------------------------------------------------------------
  // Static data members
  static double distance_inactive_;          // Inactive (dead zone) ray length
  static double distance_active_;            // Active ray length
  static unique_ptr<Source> ray_source_;     // Starting source for ray sampling
  static RandomRaySourceShape source_shape_; // Flag for linear source

  static bool first_collided_source_;
  static int first_collided_rays_;
  static int first_collided_volume_rays_;

  static bool no_volume_calc; // Flag for FCS flux calculation 

  //----------------------------------------------------------------------------
  // Public data members
  vector<float> angular_flux_;
  vector<float> angular_flux_initial_;

private:
  //----------------------------------------------------------------------------
  // Private data members
  vector<float> delta_psi_;
  vector<MomentArray> delta_moments_;

  int negroups_;
  float ray_threshold {1e-20f};


  FlatSourceDomain* domain_ {nullptr}; // pointer to domain that has flat source
                                       // data needed for ray transport
  double distance_travelled_ {0};
  bool is_active_ {false};
  bool is_alive_ {true};
}; // class RandomRay

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_H
