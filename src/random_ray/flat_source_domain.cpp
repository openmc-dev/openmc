#include "openmc/random_ray/flat_source_domain.h"

#include "openmc/cell.h"
#include "openmc/geometry.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/output.h"
#include "openmc/plot.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/simulation.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/timer.h"

#include <cstdio>

namespace openmc {

//==============================================================================
// FlatSourceDomain implementation
//==============================================================================

FlatSourceDomain::FlatSourceDomain() : negroups_(data::mg.num_energy_groups_)
{
  // Count the number of source regions, compute the cell offset
  // indices, and store the material type The reason for the offsets is that
  // some cell types may not have material fills, and therefore do not
  // produce FSRs. Thus, we cannot index into the global arrays directly
  for (const auto& c : model::cells) {
    if (c->type_ != Fill::MATERIAL) {
      source_region_offsets_.push_back(-1);
    } else {
      source_region_offsets_.push_back(n_source_regions_);
      n_source_regions_ += c->n_instances_;
      n_source_elements_ += c->n_instances_ * negroups_;
    }
  }

  // Initialize cell-wise arrays
  lock_.resize(n_source_regions_);
  material_.resize(n_source_regions_);
  position_recorded_.assign(n_source_regions_, 0);
  position_.resize(n_source_regions_);
  volume_.assign(n_source_regions_, 0.0);
  volume_t_.assign(n_source_regions_, 0.0);
  was_hit_.assign(n_source_regions_, 0);

  // Initialize element-wise arrays
  scalar_flux_new_.assign(n_source_elements_, 0.0);
  scalar_flux_old_.assign(n_source_elements_, 1.0);
  scalar_flux_final_.assign(n_source_elements_, 0.0);
  source_.resize(n_source_elements_);
  tally_task_.resize(n_source_elements_);

  // Initialize material array
  int64_t source_region_id = 0;
  for (int i = 0; i < model::cells.size(); i++) {
    Cell& cell = *model::cells[i];
    if (cell.type_ == Fill::MATERIAL) {
      for (int j = 0; j < cell.n_instances_; j++) {
        material_[source_region_id++] = cell.material(j);
      }
    }
  }

  // Sanity check
  if (source_region_id != n_source_regions_) {
    fatal_error("Unexpected number of source regions");
  }
}

void FlatSourceDomain::batch_reset()
{
  // Reset scalar fluxes, iteration volume tallies, and region hit flags to
  // zero
  parallel_fill<float>(scalar_flux_new_, 0.0f);
  parallel_fill<double>(volume_, 0.0);
  parallel_fill<int>(was_hit_, 0);
}

void FlatSourceDomain::accumulate_iteration_flux()
{
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements_; se++) {
    scalar_flux_final_[se] += scalar_flux_new_[se];
  }
}

// Compute new estimate of scattering + fission sources in each source region
// based on the flux estimate from the previous iteration.
void FlatSourceDomain::update_neutron_source(double k_eff)
{
  simulation::time_update_src.start();

  double inverse_k_eff = 1.0 / k_eff;

  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;

#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {
    int material = material_[sr];

    for (int e_out = 0; e_out < negroups_; e_out++) {
      float sigma_t = data::mg.macro_xs_[material].get_xs(
        MgxsType::TOTAL, e_out, nullptr, nullptr, nullptr, t, a);
      float scatter_source = 0.0f;
      float fission_source = 0.0f;

      for (int e_in = 0; e_in < negroups_; e_in++) {
        float scalar_flux = scalar_flux_old_[sr * negroups_ + e_in];
        float sigma_s = data::mg.macro_xs_[material].get_xs(
          MgxsType::NU_SCATTER, e_in, &e_out, nullptr, nullptr, t, a);
        float nu_sigma_f = data::mg.macro_xs_[material].get_xs(
          MgxsType::NU_FISSION, e_in, nullptr, nullptr, nullptr, t, a);
        float chi = data::mg.macro_xs_[material].get_xs(
          MgxsType::CHI_PROMPT, e_in, &e_out, nullptr, nullptr, t, a);
        scatter_source += sigma_s * scalar_flux;
        fission_source += nu_sigma_f * scalar_flux * chi;
      }

      fission_source *= inverse_k_eff;
      float new_isotropic_source = (scatter_source + fission_source) / sigma_t;
      source_[sr * negroups_ + e_out] = new_isotropic_source;
    }
  }

  simulation::time_update_src.stop();
}

// Normalizes flux and updates simulation-averaged volume estimate
void FlatSourceDomain::normalize_scalar_flux_and_volumes(
  double total_active_distance_per_iteration)
{
  float normalization_factor = 1.0 / total_active_distance_per_iteration;
  double volume_normalization_factor =
    1.0 / (total_active_distance_per_iteration * simulation::current_batch);

// Normalize scalar flux to total distance travelled by all rays this iteration
#pragma omp parallel for
  for (int64_t e = 0; e < scalar_flux_new_.size(); e++) {
    scalar_flux_new_[e] *= normalization_factor;
  }

// Accumulate cell-wise ray length tallies collected this iteration, then
// update the simulation-averaged cell-wise volume estimates
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions_; sr++) {
    volume_t_[sr] += volume_[sr];
    volume_[sr] = volume_t_[sr] * volume_normalization_factor;
  }
}

// Combine transport flux contributions and flat source contributions from the
// previous iteration to generate this iteration's estimate of scalar flux.
int64_t FlatSourceDomain::add_source_to_scalar_flux()
{
  int64_t n_hits = 0;

  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;

#pragma omp parallel for reduction(+ : n_hits)
  for (int sr = 0; sr < n_source_regions_; sr++) {

    // Check if this cell was hit this iteration
    int was_cell_hit = was_hit_[sr];
    if (was_cell_hit) {
      n_hits++;
    }

    double volume = volume_[sr];
    int material = material_[sr];
    for (int g = 0; g < negroups_; g++) {
      int64_t idx = (sr * negroups_) + g;

      // There are three scenarios we need to consider:
      if (was_cell_hit) {
        // 1. If the FSR was hit this iteration, then the new flux is equal to
        // the flat source from the previous iteration plus the contributions
        // from rays passing through the source region (computed during the
        // transport sweep)
        float sigma_t = data::mg.macro_xs_[material].get_xs(
          MgxsType::TOTAL, g, nullptr, nullptr, nullptr, t, a);
        scalar_flux_new_[idx] /= (sigma_t * volume);
        scalar_flux_new_[idx] += source_[idx];
      } else if (volume > 0.0) {
        // 2. If the FSR was not hit this iteration, but has been hit some
        // previous iteration, then we simply set the new scalar flux to be
        // equal to the contribution from the flat source alone.
        scalar_flux_new_[idx] = source_[idx];
      } else {
        // If the FSR was not hit this iteration, and it has never been hit in
        // any iteration (i.e., volume is zero), then we want to set this to 0
        // to avoid dividing anything by a zero volume.
        scalar_flux_new_[idx] = 0.0f;
      }
    }
  }

  // Return the number of source regions that were hit this iteration
  return n_hits;
}

// Generates new estimate of k_eff based on the differences between this
// iteration's estimate of the scalar flux and the last iteration's estimate.
double FlatSourceDomain::compute_k_eff(double k_eff_old) const
{
  double fission_rate_old = 0;
  double fission_rate_new = 0;

  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;

#pragma omp parallel for reduction(+ : fission_rate_old, fission_rate_new)
  for (int sr = 0; sr < n_source_regions_; sr++) {

    // If simulation averaged volume is zero, don't include this cell
    double volume = volume_[sr];
    if (volume == 0.0) {
      continue;
    }

    int material = material_[sr];

    double sr_fission_source_old = 0;
    double sr_fission_source_new = 0;

    for (int g = 0; g < negroups_; g++) {
      int64_t idx = (sr * negroups_) + g;
      double nu_sigma_f = data::mg.macro_xs_[material].get_xs(
        MgxsType::NU_FISSION, g, nullptr, nullptr, nullptr, t, a);
      sr_fission_source_old += nu_sigma_f * scalar_flux_old_[idx];
      sr_fission_source_new += nu_sigma_f * scalar_flux_new_[idx];
    }

    fission_rate_old += sr_fission_source_old * volume;
    fission_rate_new += sr_fission_source_new * volume;
  }

  double k_eff_new = k_eff_old * (fission_rate_new / fission_rate_old);

  return k_eff_new;
}

// This function is responsible for generating a mapping between random
// ray flat source regions (cell instances) and tally bins. The mapping
// takes the form of a "TallyTask" object, which accounts for one single
// score being applied to a single tally. Thus, a single source region
// may have anywhere from zero to many tally tasks associated with it ---
// meaning that the global "tally_task" data structure is in 2D. The outer
// dimension corresponds to the source element (i.e., each entry corresponds
// to a specific energy group within a specific source region), and the
// inner dimension corresponds to the tallying task itself. Mechanically,
// the mapping between FSRs and spatial filters is done by considering
// the location of a single known ray midpoint that passed through the
// FSR. I.e., during transport, the first ray to pass through a given FSR
// will write down its midpoint for use with this function. This is a cheap
// and easy way of mapping FSRs to spatial tally filters, but comes with
// the downside of adding the restriction that spatial tally filters must
// share boundaries with the physical geometry of the simulation (so as
// not to subdivide any FSR). It is acceptable for a spatial tally region
// to contain multiple FSRs, but not the other way around.

// TODO: In future work, it would be preferable to offer a more general
// (but perhaps slightly more expensive) option for handling arbitrary
// spatial tallies that would be allowed to subdivide FSRs.

// Besides generating the mapping structure, this function also keeps track
// of whether or not all flat source regions have been hit yet. This is
// required, as there is no guarantee that all flat source regions will
// be hit every iteration, such that in the first few iterations some FSRs
// may not have a known position within them yet to facilitate mapping to
// spatial tally filters. However, after several iterations, if all FSRs
// have been hit and have had a tally map generated, then this status will
// be passed back to the caller to alert them that this function doesn't
// need to be called for the remainder of the simulation.

void FlatSourceDomain::convert_source_regions_to_tallies()
{
  openmc::simulation::time_tallies.start();

  // Tracks if we've generated a mapping yet for all source regions.
  bool all_source_regions_mapped = true;

// Attempt to generate mapping for all source regions
#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {

    // If this source region has not been hit by a ray yet, then
    // we aren't going to be able to map it, so skip it.
    if (!position_recorded_[sr]) {
      all_source_regions_mapped = false;
      continue;
    }

    // A particle located at the recorded midpoint of a ray
    // crossing through this source region is used to estabilish
    // the spatial location of the source region
    Particle p;
    p.r() = position_[sr];
    p.r_last() = position_[sr];
    bool found = exhaustive_find_cell(p);

    // Loop over energy groups (so as to support energy filters)
    for (int g = 0; g < negroups_; g++) {

      // Set particle to the current energy
      p.g() = g;
      p.g_last() = g;
      p.E() = data::mg.energy_bin_avg_[p.g()];
      p.E_last() = p.E();

      int64_t source_element = sr * negroups_ + g;

      // If this task has already been populated, we don't need to do
      // it again.
      if (tally_task_[source_element].size() > 0) {
        continue;
      }

      // Loop over all active tallies. This logic is essentially identical
      // to what happens when scanning for applicable tallies during
      // MC transport.
      for (auto i_tally : model::active_tallies) {
        Tally& tally {*model::tallies[i_tally]};

        // Initialize an iterator over valid filter bin combinations.
        // If there are no valid combinations, use a continue statement
        // to ensure we skip the assume_separate break below.
        auto filter_iter = FilterBinIter(tally, p);
        auto end = FilterBinIter(tally, true, &p.filter_matches());
        if (filter_iter == end)
          continue;

        // Loop over filter bins.
        for (; filter_iter != end; ++filter_iter) {
          auto filter_index = filter_iter.index_;
          auto filter_weight = filter_iter.weight_;

          // Loop over scores
          for (auto score_index = 0; score_index < tally.scores_.size();
               score_index++) {
            auto score_bin = tally.scores_[score_index];
            // If a valid tally, filter, and score cobination has been found,
            // then add it to the list of tally tasks for this source element.
            tally_task_[source_element].emplace_back(
              i_tally, filter_index, score_index, score_bin);
          }
        }
      }
      // Reset all the filter matches for the next tally event.
      for (auto& match : p.filter_matches())
        match.bins_present_ = false;
    }
  }
  openmc::simulation::time_tallies.stop();

  mapped_all_tallies_ = all_source_regions_mapped;
}

// Tallying in random ray is not done directly during transport, rather,
// it is done only once after each power iteration. This is made possible
// by way of a mapping data structure that relates spatial source regions
// (FSRs) to tally/filter/score combinations. The mechanism by which the
// mapping is done (and the limitations incurred) is documented in the
// "convert_source_regions_to_tallies()" function comments above. The present
// tally function simply traverses the mapping data structure and executes
// the scoring operations to OpenMC's native tally result arrays.

void FlatSourceDomain::random_ray_tally() const
{
  openmc::simulation::time_tallies.start();

  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;

// We loop over all source regions and energy groups. For each
// element, we check if there are any scores needed and apply
// them.
#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {
    double volume = volume_[sr];
    double material = material_[sr];
    for (int g = 0; g < negroups_; g++) {
      int idx = sr * negroups_ + g;
      double flux = scalar_flux_new_[idx] * volume;
      for (auto& task : tally_task_[idx]) {
        double score;
        switch (task.score_type) {

        case SCORE_FLUX:
          score = flux;
          break;

        case SCORE_TOTAL:
          score = flux * data::mg.macro_xs_[material].get_xs(
                           MgxsType::TOTAL, g, NULL, NULL, NULL, t, a);
          break;

        case SCORE_FISSION:
          score = flux * data::mg.macro_xs_[material].get_xs(
                           MgxsType::FISSION, g, NULL, NULL, NULL, t, a);
          break;

        case SCORE_NU_FISSION:
          score = flux * data::mg.macro_xs_[material].get_xs(
                           MgxsType::NU_FISSION, g, NULL, NULL, NULL, t, a);
          break;

        case SCORE_EVENTS:
          score = 1.0;
          break;

        default:
          fatal_error("Invalid score specified in tallies.xml. Only flux, "
                      "total, fission, nu-fission, and events are supported in "
                      "random ray mode.");
          break;
        }
        Tally& tally {*model::tallies[task.tally_idx]};
#pragma omp atomic
        tally.results_(task.filter_idx, task.score_idx, TallyResult::VALUE) +=
          score;
      }
    }
  }
}

void FlatSourceDomain::all_reduce_replicated_source_regions()
{
#ifdef OPENMC_MPI

  // If we only have 1 MPI rank, no need
  // to reduce anything.
  if (mpi::n_procs <= 1)
    return;

  simulation::time_bank_sendrecv.start();

  // The "position_recorded" variable needs to be allreduced (and maxed),
  // as whether or not a cell was hit will affect some decisions in how the
  // source is calculated in the next iteration so as to avoid dividing
  // by zero. We take the max rather than the sum as the hit values are
  // expected to be zero or 1.
  MPI_Allreduce(MPI_IN_PLACE, position_recorded_.data(), n_source_regions_,
    MPI_INT, MPI_MAX, mpi::intracomm);

  // The position variable is more complicated to reduce than the others,
  // as we do not want the sum of all positions in each cell, rather, we
  // want to just pick any single valid position. Thus, we perform a gather
  // and then pick the first valid position we find for all source regions
  // that have had a position recorded. This operation does not need to
  // be broadcast back to other ranks, as this value is only used for the
  // tally conversion operation, which is only performed on the master rank.
  // While this is expensive, it only needs to be done for active batches,
  // and only if we have not mapped all the tallies yet. Once tallies are
  // fully mapped, then the position vector is fully populated, so this
  // operation can be skipped.

  // First, we broadcast the fully mapped tally status variable so that
  // all ranks are on the same page
  int mapped_all_tallies_i = static_cast<int>(mapped_all_tallies_);
  MPI_Bcast(&mapped_all_tallies_i, 1, MPI_INT, 0, mpi::intracomm);

  // Then, we perform the gather of position data, if needed
  if (simulation::current_batch > settings::n_inactive &&
      !mapped_all_tallies_i) {

    // Master rank will gather results and pick valid positions
    if (mpi::master) {
      // Initialize temporary vector for receiving positions
      vector<vector<Position>> all_position;
      all_position.resize(mpi::n_procs);
      for (int i = 0; i < mpi::n_procs; i++) {
        all_position[i].resize(n_source_regions_);
      }

      // Copy master rank data into gathered vector for convenience
      all_position[0] = position_;

      // Receive all data into gather vector
      for (int i = 1; i < mpi::n_procs; i++) {
        MPI_Recv(all_position[i].data(), n_source_regions_ * 3, MPI_DOUBLE, i,
          0, mpi::intracomm, MPI_STATUS_IGNORE);
      }

      // Scan through gathered data and pick first valid cell posiiton
      for (int sr = 0; sr < n_source_regions_; sr++) {
        if (position_recorded_[sr] == 1) {
          for (int i = 0; i < mpi::n_procs; i++) {
            if (all_position[i][sr].x != 0.0 || all_position[i][sr].y != 0.0 ||
                all_position[i][sr].z != 0.0) {
              position_[sr] = all_position[i][sr];
              break;
            }
          }
        }
      }
    } else {
      // Other ranks just send in their data
      MPI_Send(position_.data(), n_source_regions_ * 3, MPI_DOUBLE, 0, 0,
        mpi::intracomm);
    }
  }

  // For the rest of the source region data, we simply perform an all reduce,
  // as these values will be needed on all ranks for transport during the
  // next iteration.
  MPI_Allreduce(MPI_IN_PLACE, volume_.data(), n_source_regions_, MPI_DOUBLE,
    MPI_SUM, mpi::intracomm);

  MPI_Allreduce(MPI_IN_PLACE, was_hit_.data(), n_source_regions_, MPI_INT,
    MPI_SUM, mpi::intracomm);

  MPI_Allreduce(MPI_IN_PLACE, scalar_flux_new_.data(), n_source_elements_,
    MPI_FLOAT, MPI_SUM, mpi::intracomm);

  simulation::time_bank_sendrecv.stop();
#endif
}

// Outputs all basic material, FSR ID, multigroup flux, and
// fission source data to .vtk file that can be directly
// loaded and displayed by Paraview. Note that .vtk binary
// files require big endian byte ordering, so endianness
// is checked and flipped if necessary.
void FlatSourceDomain::output_to_vtk() const
{
  // Rename .h5 plot filename(s) to .vtk filenames
  for (int p = 0; p < model::plots.size(); p++) {
    PlottableInterface* plot = model::plots[p].get();
    plot->path_plot() =
      plot->path_plot().substr(0, plot->path_plot().find_last_of('.')) + ".vtk";
  }

  // Print header information
  print_plot();

  // Outer loop over plots
  for (int p = 0; p < model::plots.size(); p++) {

    // Get handle to OpenMC plot object and extract params
    Plot* openmc_plot = dynamic_cast<Plot*>(model::plots[p].get());

    // Random ray plots only support voxel plots
    if (!openmc_plot) {
      warning(fmt::format("Plot {} is invalid plot type -- only voxel plotting "
                          "is allowed in random ray mode.",
        p));
      continue;
    } else if (openmc_plot->type_ != Plot::PlotType::voxel) {
      warning(fmt::format("Plot {} is invalid plot type -- only voxel plotting "
                          "is allowed in random ray mode.",
        p));
      continue;
    }

    int Nx = openmc_plot->pixels_[0];
    int Ny = openmc_plot->pixels_[1];
    int Nz = openmc_plot->pixels_[2];
    Position origin = openmc_plot->origin_;
    Position width = openmc_plot->width_;
    Position ll = origin - width / 2.0;
    double x_delta = width.x / Nx;
    double y_delta = width.y / Ny;
    double z_delta = width.z / Nz;
    std::string filename = openmc_plot->path_plot();

    // Perform sanity checks on file size
    uint64_t bytes = Nx * Ny * Nz * (negroups_ + 1 + 1 + 1) * sizeof(float);
    write_message(5, "Processing plot {}: {}... (Estimated size is {} MB)",
      openmc_plot->id(), filename, bytes / 1.0e6);
    if (bytes / 1.0e9 > 1.0) {
      warning("Voxel plot specification is very large (>1 GB). Plotting may be "
              "slow.");
    } else if (bytes / 1.0e9 > 100.0) {
      fatal_error("Voxel plot specification is too large (>100 GB). Exiting.");
    }

    // Relate voxel spatial locations to random ray source regions
    vector<int> voxel_indices(Nx * Ny * Nz);

#pragma omp parallel for collapse(3)
    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
          Position sample;
          sample.z = ll.z + z_delta / 2.0 + z * z_delta;
          sample.y = ll.y + y_delta / 2.0 + y * y_delta;
          sample.x = ll.x + x_delta / 2.0 + x * x_delta;
          Particle p;
          p.r() = sample;
          bool found = exhaustive_find_cell(p);
          int i_cell = p.lowest_coord().cell;
          int64_t source_region_idx =
            source_region_offsets_[i_cell] + p.cell_instance();
          voxel_indices[z * Ny * Nx + y * Nx + x] = source_region_idx;
        }
      }
    }

    // Open file for writing
    std::FILE* plot = std::fopen(filename.c_str(), "wb");

    // Write vtk metadata
    std::fprintf(plot, "# vtk DataFile Version 2.0\n");
    std::fprintf(plot, "Dataset File\n");
    std::fprintf(plot, "BINARY\n");
    std::fprintf(plot, "DATASET STRUCTURED_POINTS\n");
    std::fprintf(plot, "DIMENSIONS %d %d %d\n", Nx, Ny, Nz);
    std::fprintf(plot, "ORIGIN 0 0 0\n");
    std::fprintf(plot, "SPACING %lf %lf %lf\n", x_delta, y_delta, z_delta);
    std::fprintf(plot, "POINT_DATA %d\n", Nx * Ny * Nz);

    // Plot multigroup flux data
    for (int g = 0; g < negroups_; g++) {
      std::fprintf(plot, "SCALARS flux_group_%d float\n", g);
      std::fprintf(plot, "LOOKUP_TABLE default\n");
      for (int fsr : voxel_indices) {
        int64_t source_element = fsr * negroups_ + g;
        float flux = scalar_flux_final_[source_element];
        flux /= (settings::n_batches - settings::n_inactive);
        flux = convert_to_big_endian<float>(flux);
        std::fwrite(&flux, sizeof(float), 1, plot);
      }
    }

    // Plot FSRs
    std::fprintf(plot, "SCALARS FSRs float\n");
    std::fprintf(plot, "LOOKUP_TABLE default\n");
    for (int fsr : voxel_indices) {
      float value = future_prn(10, fsr);
      value = convert_to_big_endian<float>(value);
      std::fwrite(&value, sizeof(float), 1, plot);
    }

    // Plot Materials
    std::fprintf(plot, "SCALARS Materials int\n");
    std::fprintf(plot, "LOOKUP_TABLE default\n");
    for (int fsr : voxel_indices) {
      int mat = material_[fsr];
      mat = convert_to_big_endian<int>(mat);
      std::fwrite(&mat, sizeof(int), 1, plot);
    }

    // Plot fission source
    std::fprintf(plot, "SCALARS total_fission_source float\n");
    std::fprintf(plot, "LOOKUP_TABLE default\n");
    for (int fsr : voxel_indices) {
      float total_fission = 0.0;
      int mat = material_[fsr];
      for (int g = 0; g < negroups_; g++) {
        int64_t source_element = fsr * negroups_ + g;
        float flux = scalar_flux_final_[source_element];
        flux /= (settings::n_batches - settings::n_inactive);
        float Sigma_f = data::mg.macro_xs_[mat].get_xs(
          MgxsType::FISSION, g, nullptr, nullptr, nullptr, 0, 0);
        total_fission += Sigma_f * flux;
      }
      total_fission = convert_to_big_endian<float>(total_fission);
      std::fwrite(&total_fission, sizeof(float), 1, plot);
    }

    std::fclose(plot);
  }
}

} // namespace openmc
