#ifndef OPENMC_RANDOM_RAY_H
#define OPENMC_RANDOM_RAY_H

#include "openmc/memory.h"
#include "openmc/particle.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/random_ray/moment_matrix.h"
#include "openmc/source.h"

#ifdef OPENMC_DAGMC_ENABLED
#include "DagMC.hpp"
#endif

namespace openmc {

// Container for MPI exchange
struct RayBufferContainer {
  Position position;
  Direction direction;
  double distance_travelled;
  vector<float> angular_flux;
  int surface;
  bool is_active;
  uint64_t ray_id; 
  int n_event;  // Number of events (surface crossings) the ray has undergone
  
  // GeometryState scalar fields
  int n_coord;
  int cell_instance;
  int n_coord_last;
  int material;
  int material_last;
  double sqrtkT;
  double sqrtkT_last;
  
  // GeometryState vector fields (sized to model::n_coord_levels at runtime)
  // LocalCoord is POD, so we can send it as contiguous bytes
  vector<LocalCoord> coord;
  
  // cell_last_ array
  vector<int> cell_last;
  
#ifdef OPENMC_DAGMC_ENABLED
  // DAGMC fields - fixed-size array to avoid variable-length vector
  Direction last_dir;
  moab::EntityHandle handles[MAX_N_HANDLES];
  int n_handles;  // Actual number of valid handles (may be less than MAX_N_HANDLES)
#endif
};

// Container for MPI exchange
struct RayExchangeData {
  Position position;
  Direction direction;
  double distance_travelled;
  int surface;
  bool is_active;
  uint64_t ray_id; 
  int n_event;  // Number of events (surface crossings) the ray has undergone
  
  // GeometryState scalar fields
  int n_coord;
  int cell_instance;
  int n_coord_last;
  int material;
  int material_last;
  double sqrtkT;
  double sqrtkT_last;
  
#ifdef OPENMC_DAGMC_ENABLED
  // DAGMC fields - fixed-size array to avoid variable-length vector
  Direction last_dir;
  moab::EntityHandle handles[MAX_N_HANDLES];
  int n_handles;  // Actual number of valid handles (may be less than MAX_N_HANDLES)
#endif
};

// Forward declare
class FlatSourceDomain;

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
  void restart_ray(FlatSourceDomain* domain, RayExchangeData& data, float* angular_flux,
                   LocalCoord* coord, int* cell_last_data);
  uint64_t transport_history_based_single_ray();
  SourceSite sample_prng();
  SourceSite sample_halton();

  bool has_left_subdomain();
  void pack_ray_for_buffer(double distance_buffer, Position position_buffer);

  //----------------------------------------------------------------------------
  // Static data members
  static double distance_inactive_;          // Inactive (dead zone) ray length
  static double distance_active_;            // Active ray length
  static unique_ptr<Source> ray_source_;     // Starting source for ray sampling
  static RandomRaySourceShape source_shape_; // Flag for linear source
  static RandomRayGeomDim geom_dim_;         // Flag for 2D vs 3D geometry
  static RandomRaySampleMethod sample_method_; // Flag for sampling method

  //----------------------------------------------------------------------------
  // Public data members
  vector<float> angular_flux_;
  RayBufferContainer exchange_data_;

  bool ray_trace_only_ {false}; // If true, only perform geometry operations
  int owner_rank_ {C_NONE}; // Rank that owns this ray based on its position

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
  bool is_local_ {true};
}; // class RandomRay

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_H
