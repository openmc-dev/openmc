#ifndef OPENMC_RANDOM_RAY_FLAT_SOURCE_DOMAIN_H
#define OPENMC_RANDOM_RAY_FLAT_SOURCE_DOMAIN_H

#include "openmc/openmp_interface.h"
#include "openmc/position.h"
#include "openmc/source.h"

namespace openmc {

/*
 * The FlatSourceDomain class encompasses data and methods for storing
 * scalar flux and source region for all flat source regions in a
 * random ray simulation domain.
 */

class FlatSourceDomain {
public:
  //----------------------------------------------------------------------------
  // Helper Structs

  // A mapping object that is used to map between a specific random ray
  // source region and an OpenMC native tally bin that it should score to
  // every iteration.
  struct TallyTask {
    int tally_idx;
    int filter_idx;
    int score_idx;
    int score_type;
    TallyTask(int tally_idx, int filter_idx, int score_idx, int score_type)
      : tally_idx(tally_idx), filter_idx(filter_idx), score_idx(score_idx),
        score_type(score_type)
    {}
  };

  //----------------------------------------------------------------------------
  // Constructors
  FlatSourceDomain();

  //----------------------------------------------------------------------------
  // Methods
  void update_neutron_source(double k_eff);
  double compute_k_eff(double k_eff_old) const;
  void normalize_scalar_flux_and_volumes(
    double total_active_distance_per_iteration);
  int64_t add_source_to_scalar_flux();
  void batch_reset();
  void convert_source_regions_to_tallies();
  void random_ray_tally() const;
  void accumulate_iteration_flux();
  void output_to_vtk() const;
  void all_reduce_replicated_source_regions();
  void apply_fixed_source_to_source_region(
    Discrete* discrete, double strength_factor, int64_t source_region);
  void apply_fixed_source_to_cell_instances(int32_t i_cell, Discrete* discrete,
    double strength_factor, int target_material_id);
  void apply_fixed_source_to_cell_and_children(int32_t i_cell,
    Discrete* discrete, double strength_factor, int32_t target_material_id);
  void convert_fixed_sources();
  void count_fixed_source_regions();
  double calculate_total_volume_weighted_source_strength() const;

  //----------------------------------------------------------------------------
  // Public Data members

  bool mapped_all_tallies_ {false}; // If all source regions have been visited

  int64_t n_source_regions_ {0}; // Total number of source regions in the model
  int64_t n_fixed_source_regions_ {0}; // Total number of source regions with
                                       // non-zero fixed source terms

  // 1D array representing source region starting offset for each OpenMC Cell
  // in model::cells
  vector<int64_t> source_region_offsets_;

  // 1D arrays representing values for all source regions
  vector<OpenMPMutex> lock_;
  vector<int> was_hit_;
  vector<double> volume_;
  vector<int> position_recorded_;
  vector<Position> position_;

  // 2D arrays stored in 1D representing values for all source regions x energy
  // groups
  vector<float> scalar_flux_old_;
  vector<float> scalar_flux_new_;
  vector<float> source_;
  vector<float> fixed_source_;

  //----------------------------------------------------------------------------
  // Private data members
private:
  int negroups_;                  // Number of energy groups in simulation
  int64_t n_source_elements_ {0}; // Total number of source regions in the model
                                  // times the number of energy groups

  // 2D array representing values for all source regions x energy groups x tally
  // tasks
  vector<vector<TallyTask>> tally_task_;

  // 1D arrays representing values for all source regions
  vector<int> material_;
  vector<double> volume_t_;

  // 2D arrays stored in 1D representing values for all source regions x energy
  // groups
  vector<float> scalar_flux_final_;

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
