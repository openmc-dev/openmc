#include "openmc/tallies/tally.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/tallies/filter.h"
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
}

double global_tally_absorption;
double global_tally_collision;
double global_tally_tracklength;
double global_tally_leakage;

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
free_memory_tally_c()
{
  #pragma omp parallel
  {
    simulation::filter_matches.clear();
  }

  model::tally_filters.clear();

  model::tallies.clear();
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_tally_get_filters(int32_t index, int32_t** indices, int* n)
{
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *indices = model::tallies[index-1]->filters_.data();
  *n = model::tallies[index-1]->filters_.size();
  return 0;
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

  void
  tally_set_filters_c(Tally* tally, int n, int32_t filter_indices[])
  {
    tally->filters_.clear();
    tally->filters_.assign(filter_indices, filter_indices + n);
  }

  int tally_get_n_filters_c(Tally* tally) {return tally->filters_.size();}

  int32_t tally_get_filter_c(Tally* tally, int i) {return tally->filters_[i];}
}

} // namespace openmc
