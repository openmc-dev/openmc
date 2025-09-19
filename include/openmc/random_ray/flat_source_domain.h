#ifndef OPENMC_RANDOM_RAY_FLAT_SOURCE_DOMAIN_H
#define OPENMC_RANDOM_RAY_FLAT_SOURCE_DOMAIN_H

#include "openmc/constants.h"
#include "openmc/openmp_interface.h"
#include "openmc/position.h"
#include "openmc/random_ray/parallel_map.h"
#include "openmc/random_ray/source_region.h"
#include "openmc/source.h"
#include <unordered_map>
#include <unordered_set>

namespace openmc {

/*
 * The FlatSourceDomain class encompasses data and methods for storing
 * scalar flux and source region for all flat source regions in a
 * random ray simulation domain.
 */

class FlatSourceDomain {
public:
  //----------------------------------------------------------------------------
  // Constructors and Destructors
  FlatSourceDomain();
  virtual ~FlatSourceDomain() = default;

  //----------------------------------------------------------------------------
  // Methods
  virtual void update_single_neutron_source(SourceRegionHandle& srh);
  virtual void update_all_neutron_sources();
  void compute_k_eff();
  virtual void normalize_scalar_flux_and_volumes(
    double total_active_distance_per_iteration);

  int64_t add_source_to_scalar_flux();
  virtual void batch_reset();
  void convert_source_regions_to_tallies(int64_t start_sr_id);
  void reset_tally_volumes();
  void random_ray_tally();
  virtual void accumulate_iteration_flux();
  void output_to_vtk() const;
  void convert_external_sources();
  void count_external_source_regions();
  void set_adjoint_sources();
  void flux_swap();
  virtual double evaluate_flux_at_point(Position r, int64_t sr, int g) const;
  double compute_fixed_source_normalization_factor() const;
  void flatten_xs();
  void transpose_scattering_matrix();
  void serialize_final_fluxes(vector<double>& flux);
  void apply_meshes();
  void apply_mesh_to_cell_instances(int32_t i_cell, int32_t mesh_idx,
    int target_material_id, const vector<int32_t>& instances,
    bool is_target_void);
  void apply_mesh_to_cell_and_children(int32_t i_cell, int32_t mesh_idx,
    int32_t target_material_id, bool is_target_void);
  SourceRegionHandle get_subdivided_source_region_handle(
    SourceRegionKey sr_key, Position r, double dist, Direction u);
  void finalize_discovered_source_regions();
  void apply_transport_stabilization();
  int64_t n_source_regions() const
  {
    return source_regions_.n_source_regions();
  }
  int64_t n_source_elements() const
  {
    return source_regions_.n_source_regions() * negroups_;
  }
  int64_t lookup_base_source_region_idx(const GeometryState& p) const;
  SourceRegionKey lookup_source_region_key(const GeometryState& p) const;
  int64_t lookup_mesh_bin(int64_t sr, Position r) const;
  int lookup_mesh_idx(int64_t sr) const;

  //----------------------------------------------------------------------------
  // Static Data members
  static bool volume_normalized_flux_tallies_;
  static bool adjoint_; // If the user wants outputs based on the adjoint flux
  static double
    diagonal_stabilization_rho_; // Adjusts strength of diagonal stabilization
                                 // for transport corrected MGXS data

  // Static variables to store source region meshes and domains
  static std::unordered_map<int, vector<std::pair<Source::DomainType, int>>>
    mesh_domain_map_;

  //----------------------------------------------------------------------------
  // Static data members
  static RandomRayVolumeEstimator volume_estimator_;

  //----------------------------------------------------------------------------
  // Public Data members
  double k_eff_ {1.0};              // Eigenvalue
  bool mapped_all_tallies_ {false}; // If all source regions have been visited

  int64_t n_external_source_regions_ {0}; // Total number of source regions with
                                          // non-zero external source terms

  // 1D array representing source region starting offset for each OpenMC Cell
  // in model::cells
  vector<int64_t> source_region_offsets_;

  // 2D arrays stored in 1D representing values for all materials x energy
  // groups
  int n_materials_;
  vector<double> sigma_t_;
  vector<double> nu_sigma_f_;
  vector<double> sigma_f_;
  vector<double> chi_;

  // 3D arrays stored in 1D representing values for all materials x energy
  // groups x energy groups
  vector<double> sigma_s_;

  // The abstract container holding all source region-specific data
  SourceRegionContainer source_regions_;

  // Parallel hash map holding all source regions discovered during
  // a single iteration. This is a threadsafe data structure that is cleaned
  // out after each iteration and stored in the "source_regions_" container.
  // It is keyed with a SourceRegionKey, which combines the base source
  // region index and the mesh bin.
  ParallelMap<SourceRegionKey, SourceRegion, SourceRegionKey::HashFunctor>
    discovered_source_regions_;

  // Map that relates a SourceRegionKey to the index at which the source
  // region can be found in the "source_regions_" container.
  std::unordered_map<SourceRegionKey, int64_t, SourceRegionKey::HashFunctor>
    source_region_map_;

  // Map that relates a SourceRegionKey to the external source index. This map
  // is used to check if there are any point sources within a subdivided source
  // region at the time it is discovered.
  std::unordered_map<SourceRegionKey, vector<int>, SourceRegionKey::HashFunctor>
    external_point_source_map_;

  // Map that relates a base source region index to the external source index.
  // This map is used to check if there are any volumetric sources within a
  // subdivided source region at the time it is discovered.
  std::unordered_map<int64_t, vector<int>> external_volumetric_source_map_;

  // Map that relates a base source region index to a mesh index. This map
  // is used to check which subdivision mesh is present in a source region.
  std::unordered_map<int64_t, int> mesh_map_;

  // If transport corrected MGXS data is being used, there may be negative
  // in-group scattering cross sections that can result in instability in MOC
  // and random ray if used naively. This flag enables a stabilization
  // technique.
  bool is_transport_stabilization_needed_ {false};

protected:
  //----------------------------------------------------------------------------
  // Methods
  void apply_external_source_to_source_region(
    int src_idx, SourceRegionHandle& srh);
  void apply_external_source_to_cell_instances(int32_t i_cell, int src_idx,
    int target_material_id, const vector<int32_t>& instances);
  void apply_external_source_to_cell_and_children(
    int32_t i_cell, int src_idx, int32_t target_material_id);
  virtual void set_flux_to_flux_plus_source(int64_t sr, double volume, int g);
  void set_flux_to_source(int64_t sr, int g);
  virtual void set_flux_to_old_flux(int64_t sr, int g);

  //----------------------------------------------------------------------------
  // Private data members
  int negroups_; // Number of energy groups in simulation

  double
    simulation_volume_; // Total physical volume of the simulation domain, as
                        // defined by the 3D box of the random ray source

  // Volumes for each tally and bin/score combination. This intermediate data
  // structure is used when tallying quantities that must be normalized by
  // volume (i.e., flux). The vector is index by tally index, while the inner 2D
  // xtensor is indexed by bin index and score index in a similar manner to the
  // results tensor in the Tally class, though without the third dimension, as
  // SUM and SUM_SQ do not need to be tracked.
  vector<xt::xtensor<double, 2>> tally_volumes_;

}; // class FlatSourceDomain

//============================================================================
//! Non-member functions
//============================================================================

// Returns the inputted value in big endian byte ordering. If the system is
// little endian, the byte ordering is flipped. If the system is big endian,
// the inputted value is returned as is. This function is necessary as
// .vtk binary files use big endian byte ordering.
template<typename T>
T convert_to_big_endian(T in)
{
  // 4 byte integer
  uint32_t test = 1;

  // 1 byte pointer to first byte of test integer
  uint8_t* ptr = reinterpret_cast<uint8_t*>(&test);

  // If the first byte of test is 0, then the system is big endian. In this
  // case, we don't have to do anything as .vtk files are big endian
  if (*ptr == 0)
    return in;

  // Otherwise, the system is in little endian, so we need to flip the
  // endianness
  uint8_t* orig = reinterpret_cast<uint8_t*>(&in);
  uint8_t swapper[sizeof(T)];
  for (int i = 0; i < sizeof(T); i++) {
    swapper[i] = orig[sizeof(T) - i - 1];
  }
  T out = *reinterpret_cast<T*>(&swapper);
  return out;
}

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_FLAT_SOURCE_DOMAIN_H
