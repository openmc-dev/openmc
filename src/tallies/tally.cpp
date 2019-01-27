#include "openmc/tallies/tally.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_delayedgroup.h"
#include "openmc/tallies/filter_surface.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/message_passing.h"

#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp" // for empty_like
#include "xtensor/xview.hpp"

#include <array>
#include <cstddef>

namespace openmc {

//==============================================================================
// Global variable definitions
//==============================================================================

namespace model {
  std::vector<std::unique_ptr<Tally>> tallies;

  std::vector<int> active_tallies;
  std::vector<int> active_analog_tallies;
  std::vector<int> active_tracklength_tallies;
  std::vector<int> active_collision_tallies;
  //std::vector<int> active_meshsurf_tallies;
  //std::vector<int> active_surface_tallies;
}

double global_tally_absorption;
double global_tally_collision;
double global_tally_tracklength;
double global_tally_leakage;

//==============================================================================// Tally object implementation
//==============================================================================

void
Tally::set_filters(const int32_t filter_indices[], int n)
{
  // Clear old data.
  filters_.clear();
  strides_.clear();

  // Copy in the given filter indices.
  filters_.assign(filter_indices, filter_indices + n);

  for (int i = 0; i < n; ++i) {
    auto i_filt = filters_[i];
    //TODO: off-by-one
    if (i_filt < 1 || i_filt > model::tally_filters.size())
      throw std::out_of_range("Index in tally filter array out of bounds.");

    //TODO: off-by-one
    const auto* filt = model::tally_filters[i_filt-1].get();

    //TODO: off-by-one on each index
    if (dynamic_cast<const EnergyFilter*>(filt)) {
      if (dynamic_cast<const EnergyoutFilter*>(filt)) {
        energyout_filter_ = i + 1;
      } else {
        energyin_filter_ = i + 1;
      }
    } else if (dynamic_cast<const DelayedGroupFilter*>(filt)) {
      delayedgroup_filter_ = i + 1;
    } else if (dynamic_cast<const SurfaceFilter*>(filt)) {
      surface_filter_ = i + 1;
    } else if (dynamic_cast<const MeshFilter*>(filt)) {
      mesh_filter_ = i + 1;
    }
  }

  // Set the strides.  Filters are traversed in reverse so that the last filter
  // has the shortest stride in memory and the first filter has the longest
  // stride.
  strides_.resize(n, 0);
  int stride = 1;
  for (int i = n-1; i >= 0; --i) {
    strides_[i] = stride;
    //TODO: off-by-one
    stride *= model::tally_filters[filters_[i]-1]->n_bins_;
  }
  n_filter_bins_ = stride;
}

//==============================================================================
// Non-member functions
//==============================================================================

adaptor_type<2> global_tallies()
{
  // Get pointer to global tallies
  double* buffer;
  openmc_global_tallies(&buffer);

  // Adapt into xtensor
  std::array<size_t, 2> shape = {N_GLOBAL_TALLIES, 3};
  std::size_t size {3*N_GLOBAL_TALLIES};

  return xt::adapt(buffer, size, xt::no_ownership(), shape);
}

adaptor_type<3> tally_results(int idx)
{
  // Get pointer to tally results
  double* results;
  std::array<std::size_t, 3> shape;
  openmc_tally_results(idx, &results, shape.data());

  // Adapt array into xtensor with no ownership
  std::size_t size {shape[0] * shape[1] * shape[2]};
  return xt::adapt(results, size, xt::no_ownership(), shape);
}

#ifdef OPENMC_MPI
void reduce_tally_results()
{
  for (int i = 1; i <= n_tallies; ++i) {
    // Skip any tallies that are not active
    bool active;
    openmc_tally_get_active(i, &active);
    if (!active) continue;

    // Get view of accumulated tally values
    auto results = tally_results(i);
    auto values_view = xt::view(results, xt::all(), xt::all(), RESULT_VALUE);

    // Make copy of tally values in contiguous array
    xt::xtensor<double, 2> values = values_view;
    xt::xtensor<double, 2> values_reduced = xt::empty_like(values);

    // Reduce contiguous set of tally results
    MPI_Reduce(values.data(), values_reduced.data(), values.size(),
      MPI_DOUBLE, MPI_SUM, 0, mpi::intracomm);

    // Transfer values on master and reset on other ranks
    if (mpi::master) {
      values_view = values_reduced;
    } else {
      values_view = 0.0;
    }
  }

  // Get view of global tally values
  auto gt = global_tallies();
  auto gt_values_view = xt::view(gt, xt::all(), RESULT_VALUE);

  // Make copy of values in contiguous array
  xt::xtensor<double, 1> gt_values = gt_values_view;
  xt::xtensor<double, 1> gt_values_reduced = xt::empty_like(gt_values);

  // Reduce contiguous data
  MPI_Reduce(gt_values.data(), gt_values_reduced.data(), N_GLOBAL_TALLIES,
    MPI_DOUBLE, MPI_SUM, 0, mpi::intracomm);

  // Transfer values on master and reset on other ranks
  if (mpi::master) {
    gt_values_view = gt_values_reduced;
  } else {
    gt_values_view = 0.0;
  }

  // We also need to determine the total starting weight of particles from the
  // last realization
  double weight_reduced;
  MPI_Reduce(&total_weight, &weight_reduced, 1, MPI_DOUBLE, MPI_SUM,
    0, mpi::intracomm);
  if (mpi::master) total_weight = weight_reduced;
}
#endif

extern "C" void
setup_active_tallies_c()
{
  model::active_tallies.clear();
  model::active_analog_tallies.clear();
  model::active_tracklength_tallies.clear();
  model::active_collision_tallies.clear();
  //model::active_meshsurf_tallies.clear();
  //model::active_surface_tallies.clear();

  //TODO: off-by-one all through here
  for (auto i = 0; i < model::tallies.size(); ++i) {
    const auto& tally {*model::tallies[i]};

    if (tally.active_) {
      model::active_tallies.push_back(i + 1);
      switch (tally.type_) {

      case TALLY_VOLUME:
        switch (tally.estimator_) {
          case ESTIMATOR_ANALOG:
            model::active_analog_tallies.push_back(i + 1);
            break;
          case ESTIMATOR_TRACKLENGTH:
            model::active_tracklength_tallies.push_back(i + 1);
            break;
          case ESTIMATOR_COLLISION:
            model::active_collision_tallies.push_back(i + 1);
        }
        break;

      //case TALLY_MESH_SURFACE:
      //  model::active_meshsurf_tallies.push_back(i + 1);
      //  break;

      //case TALLY_SURFACE:
      //  model::active_surface_tallies.push_back(i + 1);
      }
    }
  }
}

extern "C" void
free_memory_tally_c()
{
  #pragma omp parallel
  {
    simulation::filter_matches.clear();
    model::tally_derivs.clear();
  }

  model::tally_filters.clear();

  model::tallies.clear();

  model::active_tallies.clear();
  model::active_analog_tallies.clear();
  model::active_tracklength_tallies.clear();
  model::active_collision_tallies.clear();
  //model::active_meshsurf_tallies.clear();
  //model::active_surface_tallies.clear();
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_tally_get_type(int32_t index, int32_t* type)
{
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  //TODO: off-by-one
  *type = model::tallies[index-1]->type_;
}

extern "C" int
openmc_tally_set_type(int32_t index, const char* type)
{
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  if (strcmp(type, "volume") == 0) {
    model::tallies[index-1]->type_ = TALLY_VOLUME;
  } else if (strcmp(type, "mesh-surface") == 0) {
    model::tallies[index-1]->type_ = TALLY_MESH_SURFACE;
  } else if (strcmp(type, "surface") == 0) {
    model::tallies[index-1]->type_ = TALLY_SURFACE;
  } else {
    std::stringstream errmsg;
    errmsg << "Unknown tally type: " << type;
    set_errmsg(errmsg);
    return OPENMC_E_INVALID_ARGUMENT;
  }
}

extern "C" int
openmc_tally_get_active(int32_t index, bool* active)
{
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  //TODO: off-by-one
  *active = model::tallies[index-1]->active_;
}

extern "C" int
openmc_tally_set_active(int32_t index, bool active)
{
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  //TODO: off-by-one
  model::tallies[index-1]->active_ = active;
}

extern "C" int
openmc_tally_get_filters(int32_t index, const int32_t** indices, int* n)
{
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  //TODO: off-by-one
  *indices = model::tallies[index-1]->filters().data();
  *n = model::tallies[index-1]->filters().size();
  return 0;
}

extern "C" int
openmc_tally_set_filters(int32_t index, int n, const int32_t* indices)
{
  // Make sure the index fits in the array bounds.
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  // Set the filters.
  try {
    //TODO: off-by-one
    model::tallies[index-1]->set_filters(indices, n);
  } catch (const std::out_of_range& ex) {
    set_errmsg(ex.what());
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  return 0;
}

extern "C" int
openmc_tally_get_nuclides(int32_t index, int** nuclides, int* n)
{
  // Make sure the index fits in the array bounds.
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  //TODO: off-by-one
  *n = model::tallies[index-1]->nuclides_.size();
  *nuclides = model::tallies[index-1]->nuclides_.data();
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Tally* tally_pointer(int indx) {return model::tallies[indx].get();}

  void
  extend_tallies_c(int n)
  {
    for (int i = 0; i < n; ++i)
      model::tallies.push_back(std::make_unique<Tally>());
  }

  int active_tallies_data(int i)
  {return model::active_tallies[i-1];}

  int active_tallies_size()
  {return model::active_tallies.size();}

  int active_analog_tallies_data(int i)
  {return model::active_analog_tallies[i-1];}

  int active_analog_tallies_size()
  {return model::active_analog_tallies.size();}

  int active_tracklength_tallies_data(int i)
  {return model::active_tracklength_tallies[i-1];}

  int active_tracklength_tallies_size()
  {return model::active_tracklength_tallies.size();}

  int active_collision_tallies_data(int i)
  {return model::active_collision_tallies[i-1];}

  int active_collision_tallies_size()
  {return model::active_collision_tallies.size();}

  int tally_get_type_c(Tally* tally) {return tally->type_;}

  void tally_set_type_c(Tally* tally, int type) {tally->type_ = type;}

  int tally_get_estimator_c(Tally* tally) {return tally->estimator_;}

  void tally_set_estimator_c(Tally* tally, int e) {tally->estimator_ = e;}

  void tally_set_filters_c(Tally* tally, int n, int32_t filter_indices[])
  {tally->set_filters(filter_indices, n);}

  int tally_get_n_filters_c(Tally* tally) {return tally->filters().size();}

  int32_t tally_get_filter_c(Tally* tally, int i) {return tally->filters(i);}

  int32_t tally_get_stride_c(Tally* tally, int i) {return tally->strides(i);}

  int32_t tally_get_n_filter_bins_c(Tally* tally)
  {return tally->n_filter_bins();}

  int tally_get_n_nuclide_bins_c(Tally* tally)
  {return tally->nuclides_.size();}

  int tally_get_nuclide_bins_c(Tally* tally, int i)
  {return tally->nuclides_[i-1];}

  void
  tally_set_nuclide_bins_c(Tally* tally, int n, int bins[], bool all_nuclides)
  {
    tally->nuclides_.clear();
    tally->nuclides_.assign(bins, bins + n);
    tally->all_nuclides_ = all_nuclides;
  }

  bool tally_get_all_nuclides_c(Tally* tally) {return tally->all_nuclides_;}

  int tally_get_energyin_filter_c(Tally* tally)
  {return tally->energyin_filter_;}

  int tally_get_energyout_filter_c(Tally* tally)
  {return tally->energyout_filter_;}

  int tally_get_delayedgroup_filter_c(Tally* tally)
  {return tally->delayedgroup_filter_;}

  int tally_get_surface_filter_c(Tally* tally)
  {return tally->surface_filter_;}

  int tally_get_mesh_filter_c(Tally* tally)
  {return tally->mesh_filter_;}

  int tally_get_deriv_c(Tally* tally) {return tally->deriv_;}

  int tally_set_deriv_c(Tally* tally, int deriv) {tally->deriv_ = deriv;}
}

} // namespace openmc
