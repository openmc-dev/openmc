#ifndef OPENMC_RANDOM_RAY_FLAT_SOURCE_DOMAIN_H
#define OPENMC_RANDOM_RAY_FLAT_SOURCE_DOMAIN_H

#include "openmc/constants.h"
#include "openmc/openmp_interface.h"
#include "openmc/position.h"
#include "openmc/source.h"
#include "openmc/random_ray/source_region.h"
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
  virtual void update_neutron_source(double k_eff);
  double compute_k_eff(double k_eff_old) const;
  virtual void normalize_scalar_flux_and_volumes(
    double total_active_distance_per_iteration);

  int64_t add_source_to_scalar_flux();
  virtual void batch_reset();
  void convert_source_regions_to_tallies();
  void reset_tally_volumes();
  void random_ray_tally();
  virtual void accumulate_iteration_flux();
  void output_to_vtk() const;
  virtual void all_reduce_replicated_source_regions();
  void convert_external_sources();
  void count_external_source_regions();
  void set_adjoint_sources(const vector<double>& forward_flux);
  virtual void flux_swap();
  virtual double evaluate_flux_at_point(Position r, int64_t sr, int g) const;
  double compute_fixed_source_normalization_factor() const;
  void flatten_xs();
  void transpose_scattering_matrix();

  //----------------------------------------------------------------------------
  // Static Data members
  static bool volume_normalized_flux_tallies_;
  static bool adjoint_; // If the user wants outputs based on the adjoint flux

  //----------------------------------------------------------------------------
  // Static data members
  static RandomRayVolumeEstimator volume_estimator_;

  //----------------------------------------------------------------------------
  // Public Data members

  bool mapped_all_tallies_ {false}; // If all source regions have been visited

  int64_t n_source_regions_ {0}; // Total number of source regions in the model
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

  // A 1D array over all source regions, with each source region having a set of
  // volume tally tasks. This more complicated data structure is convenient for
  // ensuring that volumes are only tallied once per source region, regardless
  // of how many energy groups are used for tallying.
  vector<std::unordered_set<TallyTask, TallyTask::HashFunctor>> volume_task_;

  // A 1D array of source regions
  vector<SourceRegion> source_regions_;

protected:
  //----------------------------------------------------------------------------
  // Methods
  void apply_external_source_to_source_region(
    Discrete* discrete, double strength_factor, int64_t source_region);
  void apply_external_source_to_cell_instances(int32_t i_cell,
    Discrete* discrete, double strength_factor, int target_material_id,
    const vector<int32_t>& instances);
  void apply_external_source_to_cell_and_children(int32_t i_cell,
    Discrete* discrete, double strength_factor, int32_t target_material_id);
  virtual void set_flux_to_flux_plus_source(
    int sr, double volume, int g);
  void set_flux_to_source(int sr, int g);
  virtual void set_flux_to_old_flux(int sr, int g);

  //----------------------------------------------------------------------------
  // Private data members
  int negroups_;                  // Number of energy groups in simulation
  int64_t n_source_elements_ {0}; // Total number of source regions in the model
                                  // times the number of energy groups

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

template<typename T>
void parallel_fill(vector<T>& arr, T value)
{
#pragma omp parallel for schedule(static)
  for (int i = 0; i < arr.size(); i++) {
    arr[i] = value;
  }
}

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_FLAT_SOURCE_DOMAIN_H
