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
  RandomRay(uint64_t ray_id, FlatSourceDomain* domain);

  //----------------------------------------------------------------------------
  // Methods
  void event_advance_ray();
  void attenuate_flux(double distance, bool is_active, double offset = 0.0);
  void attenuate_flux_inner(
    double distance, bool is_active, int64_t sr, int mesh_bin, Position r);
  void attenuate_flux_flat_source(
    SourceRegionHandle& srh, double distance, bool is_active, Position r);
  void attenuate_flux_flat_source_void(
    SourceRegionHandle& srh, double distance, bool is_active, Position r);
  void attenuate_flux_linear_source(
    SourceRegionHandle& srh, double distance, bool is_active, Position r);
  void attenuate_flux_linear_source_void(
    SourceRegionHandle& srh, double distance, bool is_active, Position r);

  void initialize_ray(uint64_t ray_id, FlatSourceDomain* domain);
  uint64_t transport_history_based_single_ray();
  SourceSite sample_prng();
  SourceSite sample_halton();

  //----------------------------------------------------------------------------
  // Static data members
  static double distance_inactive_;          // Inactive (dead zone) ray length
  static double distance_active_;            // Active ray length
  static unique_ptr<Source> ray_source_;     // Starting source for ray sampling
  static RandomRaySourceShape source_shape_; // Flag for linear source
  static RandomRaySampleMethod sample_method_; // Flag for sampling method

  static double avg_miss_rate_;     // Average ray miss rate per
                                    // iteration for reporting
  static int64_t n_source_regions_; // Total number of source regions
  static int64_t
    n_external_source_regions_; // Total number of source regions with
                                // non-zero external source terms
  static uint64_t total_geometric_intersections_; // Tracks the total number of
                                                  // geometric intersections by
                                                  // all rays for reporting

  // Kinetic simulation variables
  static int bd_order_; // Order of backwards difference approximation for
                        // time-derivatives
  static RandomRayTimeMethod
    time_method_; // Method for resolving angular flux time derivative term for
                  // flux propogation

  //----------------------------------------------------------------------------
  // Public data members
  vector<float> angular_flux_;
  bool ray_trace_only_ {false}; // If true, only perform geometry operations

  //---------------------------------------------------------------------------
  // Public data members for kinetic simulations
  vector<float> angular_flux_td_;
  vector<float> angular_flux_td_prime_;

private:
  //----------------------------------------------------------------------------
  // Private data members
  vector<float> delta_psi_;
  vector<MomentArray> delta_moments_;
  vector<int> mesh_bins_;
  vector<double> mesh_fractional_lengths_;

  int negroups_;
  FlatSourceDomain* domain_ {nullptr}; // pointer to domain that has flat source
                                       // data needed for ray transport
  double distance_travelled_ {0};
  bool is_active_ {false};
  bool is_alive_ {true};

  //---------------------------------------------------------------------------
  // Private data members for kinetic simulations
  vector<double> delta_psi_td_;
  vector<double> delta_psi_td_prime_;

}; // class RandomRay

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_H
