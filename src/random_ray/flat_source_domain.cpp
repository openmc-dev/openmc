#include "openmc/random_ray/flat_source_domain.h"

#include "openmc/cell.h"
#include "openmc/eigenvalue.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
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

// Static Variable Declarations
RandomRayVolumeEstimator FlatSourceDomain::volume_estimator_ {
  RandomRayVolumeEstimator::HYBRID};
bool FlatSourceDomain::volume_normalized_flux_tallies_ {false};

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
  volume_naive_.assign(n_source_regions_, 0.0);

  // Initialize element-wise arrays
  scalar_flux_new_.assign(n_source_elements_, 0.0);
  scalar_flux_final_.assign(n_source_elements_, 0.0);
  source_.resize(n_source_elements_);

  tally_task_.resize(n_source_elements_);
  volume_task_.resize(n_source_regions_);

  if (settings::run_mode == RunMode::EIGENVALUE) {
    // If in eigenvalue mode, set starting flux to guess of unity
    scalar_flux_old_.assign(n_source_elements_, 1.0);
  } else {
    // If in fixed source mode, set starting flux to guess of zero
    // and initialize external source arrays
    scalar_flux_old_.assign(n_source_elements_, 0.0);
    external_source_.assign(n_source_elements_, 0.0);
    external_source_present_.assign(n_source_regions_, false);
  }

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

  // Initialize tally volumes
  if (volume_normalized_flux_tallies_) {
    tally_volumes_.resize(model::tallies.size());
    for (int i = 0; i < model::tallies.size(); i++) {
      //  Get the shape of the 3D result tensor
      auto shape = model::tallies[i]->results().shape();

      // Create a new 2D tensor with the same size as the first
      // two dimensions of the 3D tensor
      tally_volumes_[i] =
        xt::xtensor<double, 2>::from_shape({shape[0], shape[1]});
    }
  }

  // Compute simulation domain volume based on ray source
  auto* is = dynamic_cast<IndependentSource*>(RandomRay::ray_source_.get());
  SpatialDistribution* space_dist = is->space();
  SpatialBox* sb = dynamic_cast<SpatialBox*>(space_dist);
  Position dims = sb->upper_right() - sb->lower_left();
  simulation_volume_ = dims.x * dims.y * dims.z;
}

void FlatSourceDomain::batch_reset()
{
  // Reset scalar fluxes, iteration volume tallies, and region hit flags to
  // zero
  parallel_fill<double>(scalar_flux_new_, 0.0);
  parallel_fill<double>(volume_, 0.0);
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
  // TODO: Currently assumes we are only using single temp/single angle data.
  const int t = 0;
  const int a = 0;

  // Add scattering source
#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {
    int material = material_[sr];

    for (int e_out = 0; e_out < negroups_; e_out++) {
      float sigma_t = data::mg.macro_xs_[material].get_xs(
        MgxsType::TOTAL, e_out, nullptr, nullptr, nullptr, t, a);
      float scatter_source = 0.0f;

      for (int e_in = 0; e_in < negroups_; e_in++) {
        float scalar_flux = scalar_flux_old_[sr * negroups_ + e_in];

        float sigma_s = data::mg.macro_xs_[material].get_xs(
          MgxsType::NU_SCATTER, e_in, &e_out, nullptr, nullptr, t, a);
        scatter_source += sigma_s * scalar_flux;
      }

      source_[sr * negroups_ + e_out] = scatter_source / sigma_t;
    }
  }

  if (settings::run_mode == RunMode::EIGENVALUE) {
    // Add fission source if in eigenvalue mode
#pragma omp parallel for
    for (int sr = 0; sr < n_source_regions_; sr++) {
      int material = material_[sr];

      for (int e_out = 0; e_out < negroups_; e_out++) {
        float sigma_t = data::mg.macro_xs_[material].get_xs(
          MgxsType::TOTAL, e_out, nullptr, nullptr, nullptr, t, a);
        float fission_source = 0.0f;

        for (int e_in = 0; e_in < negroups_; e_in++) {
          float scalar_flux = scalar_flux_old_[sr * negroups_ + e_in];
          float nu_sigma_f = data::mg.macro_xs_[material].get_xs(
            MgxsType::NU_FISSION, e_in, nullptr, nullptr, nullptr, t, a);
          float chi = data::mg.macro_xs_[material].get_xs(
            MgxsType::CHI_PROMPT, e_in, &e_out, nullptr, nullptr, t, a);
          fission_source += nu_sigma_f * scalar_flux * chi;
        }
        source_[sr * negroups_ + e_out] +=
          fission_source * inverse_k_eff / sigma_t;
      }
    }
  } else {
// Add external source if in fixed source mode
#pragma omp parallel for
    for (int se = 0; se < n_source_elements_; se++) {
      source_[se] += external_source_[se];
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
    volume_naive_[sr] = volume_[sr] * normalization_factor;
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

    double volume_simulation_avg = volume_[sr];
    double volume_iteration = volume_naive_[sr];

    // Increment the number of hits if cell was hit this iteration
    if (volume_iteration) {
      n_hits++;
    }

    // Check if an external source is present in this source region
    bool external_source_present =
      external_source_present_.size() && external_source_present_[sr];

    // The volume treatment depends on the volume estimator type
    // and whether or not an external source is present in the cell.
    double volume;
    switch (volume_estimator_) {
    case RandomRayVolumeEstimator::NAIVE:
      volume = volume_iteration;
      break;
    case RandomRayVolumeEstimator::SIMULATION_AVERAGED:
      volume = volume_simulation_avg;
      break;
    case RandomRayVolumeEstimator::HYBRID:
      if (external_source_present) {
        volume = volume_iteration;
      } else {
        volume = volume_simulation_avg;
      }
      break;
    default:
      fatal_error("Invalid volume estimator type");
    }

    int material = material_[sr];
    for (int g = 0; g < negroups_; g++) {
      int64_t idx = (sr * negroups_) + g;
      float flux;

      // There are three scenarios we need to consider:
      if (volume_iteration > 0.0) {
        // 1. If the FSR was hit this iteration, then the new flux is equal to
        // the flat source from the previous iteration plus the contributions
        // from rays passing through the source region (computed during the
        // transport sweep)
        float sigma_t = data::mg.macro_xs_[material].get_xs(
          MgxsType::TOTAL, g, nullptr, nullptr, nullptr, t, a);
        flux = scalar_flux_new_[idx];
        flux /= (sigma_t * volume);
        flux += source_[idx];
      } else if (volume_simulation_avg > 0.0) {
        // 2. If the FSR was not hit this iteration, but has been hit some
        // previous iteration, then we need to make a choice about what
        // to do. Naively we will usually want to set the flux to be equal
        // to the reduced source. However, in fixed source problems where
        // there is a strong external source present in the cell, and where
        // the cell has a very low cross section, this approximation will
        // cause a huge upward bias in the flux estimate of the cell (in these
        // conditions, the flux estimate can be orders of magnitude too large).
        // Thus, to avoid this bias, if any external source is present
        // in the cell we will use the previous iteration's flux estimate. This
        // injects a small degree of correlation into the simulation, but this
        // is going to be trivial when the miss rate is a few percent or less.
        if (external_source_present) {
          flux = scalar_flux_old_[idx];
        } else {
          flux = source_[idx];
        }
      } else {
        // If the FSR was not hit this iteration, and it has never been hit in
        // any iteration (i.e., volume is zero), then we want to set this to 0
        // to avoid dividing anything by a zero volume.
        flux = 0.0f;
      }

      // Write the new scalar flux to the array
      scalar_flux_new_[idx] = flux;
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

  // Vector for gathering fission source terms for Shannon entropy calculation
  vector<float> p(n_source_regions_, 0.0f);

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

    // Compute total fission rates in FSR
    sr_fission_source_old *= volume;
    sr_fission_source_new *= volume;

    // Accumulate totals
    fission_rate_old += sr_fission_source_old;
    fission_rate_new += sr_fission_source_new;

    // Store total fission rate in the FSR for Shannon calculation
    p[sr] = sr_fission_source_new;
  }

  double k_eff_new = k_eff_old * (fission_rate_new / fission_rate_old);

  double H = 0.0;
  // defining an inverse sum for better performance
  double inverse_sum = 1 / fission_rate_new;

#pragma omp parallel for reduction(+ : H)
  for (int sr = 0; sr < n_source_regions_; sr++) {
    // Only if FSR has non-negative and non-zero fission source
    if (p[sr] > 0.0f) {
      // Normalize to total weight of bank sites. p_i for better performance
      float p_i = p[sr] * inverse_sum;
      // Sum values to obtain Shannon entropy.
      H -= p_i * std::log2(p_i);
    }
  }

  // Adds entropy value to shared entropy vector in openmc namespace.
  simulation::entropy.push_back(H);

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
            // If a valid tally, filter, and score combination has been found,
            // then add it to the list of tally tasks for this source element.
            TallyTask task(i_tally, filter_index, score_index, score_bin);
            tally_task_[source_element].push_back(task);

            // Also add this task to the list of volume tasks for this source
            // region.
            volume_task_[sr].insert(task);
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

// Set the volume accumulators to zero for all tallies
void FlatSourceDomain::reset_tally_volumes()
{
  if (volume_normalized_flux_tallies_) {
#pragma omp parallel for
    for (int i = 0; i < tally_volumes_.size(); i++) {
      auto& tensor = tally_volumes_[i];
      tensor.fill(0.0); // Set all elements of the tensor to 0.0
    }
  }
}

// In fixed source mode, due to the way that volumetric fixed sources are
// converted and applied as volumetric sources in one or more source regions,
// we need to perform an additional normalization step to ensure that the
// reported scalar fluxes are in units per source neutron. This allows for
// direct comparison of reported tallies to Monte Carlo flux results.
// This factor needs to be computed at each iteration, as it is based on the
// volume estimate of each FSR, which improves over the course of the simulation
double FlatSourceDomain::compute_fixed_source_normalization_factor() const
{
  // If we are not in fixed source mode, then there are no external sources
  // so no normalization is needed.
  if (settings::run_mode != RunMode::FIXED_SOURCE) {
    return 1.0;
  }

  // Step 1 is to sum over all source regions and energy groups to get the
  // total external source strength in the simulation.
  double simulation_external_source_strength = 0.0;
#pragma omp parallel for reduction(+ : simulation_external_source_strength)
  for (int sr = 0; sr < n_source_regions_; sr++) {
    int material = material_[sr];
    double volume = volume_[sr] * simulation_volume_;
    for (int e = 0; e < negroups_; e++) {
      // Temperature and angle indices, if using multiple temperature
      // data sets and/or anisotropic data sets.
      // TODO: Currently assumes we are only using single temp/single
      // angle data.
      const int t = 0;
      const int a = 0;
      float sigma_t = data::mg.macro_xs_[material].get_xs(
        MgxsType::TOTAL, e, nullptr, nullptr, nullptr, t, a);
      simulation_external_source_strength +=
        external_source_[sr * negroups_ + e] * sigma_t * volume;
    }
  }

  // Step 2 is to determine the total user-specified external source strength
  double user_external_source_strength = 0.0;
  for (auto& ext_source : model::external_sources) {
    user_external_source_strength += ext_source->strength();
  }

  // The correction factor is the ratio of the user-specified external source
  // strength to the simulation external source strength.
  double source_normalization_factor =
    user_external_source_strength / simulation_external_source_strength;

  return source_normalization_factor;
}

// Tallying in random ray is not done directly during transport, rather,
// it is done only once after each power iteration. This is made possible
// by way of a mapping data structure that relates spatial source regions
// (FSRs) to tally/filter/score combinations. The mechanism by which the
// mapping is done (and the limitations incurred) is documented in the
// "convert_source_regions_to_tallies()" function comments above. The present
// tally function simply traverses the mapping data structure and executes
// the scoring operations to OpenMC's native tally result arrays.

void FlatSourceDomain::random_ray_tally()
{
  openmc::simulation::time_tallies.start();

  // Reset our tally volumes to zero
  reset_tally_volumes();

  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single
  // angle data.
  const int t = 0;
  const int a = 0;

  double source_normalization_factor =
    compute_fixed_source_normalization_factor();

// We loop over all source regions and energy groups. For each
// element, we check if there are any scores needed and apply
// them.
#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {
    // The fsr.volume_ is the unitless fractional simulation averaged volume
    // (i.e., it is the FSR's fraction of the overall simulation volume). The
    // simulation_volume_ is the total 3D physical volume in cm^3 of the entire
    // global simulation domain (as defined by the ray source box). Thus, the
    // FSR's true 3D spatial volume in cm^3 is found by multiplying its fraction
    // of the total volume by the total volume. Not important in eigenvalue
    // solves, but useful in fixed source solves for returning the flux shape
    // with a magnitude that makes sense relative to the fixed source strength.
    double volume = volume_[sr] * simulation_volume_;

    double material = material_[sr];
    for (int g = 0; g < negroups_; g++) {
      int idx = sr * negroups_ + g;
      double flux = scalar_flux_new_[idx] * source_normalization_factor;

      // Determine numerical score value
      for (auto& task : tally_task_[idx]) {
        double score;
        switch (task.score_type) {

        case SCORE_FLUX:
          score = flux * volume;
          break;

        case SCORE_TOTAL:
          score = flux * volume *
                  data::mg.macro_xs_[material].get_xs(
                    MgxsType::TOTAL, g, NULL, NULL, NULL, t, a);
          break;

        case SCORE_FISSION:
          score = flux * volume *
                  data::mg.macro_xs_[material].get_xs(
                    MgxsType::FISSION, g, NULL, NULL, NULL, t, a);
          break;

        case SCORE_NU_FISSION:
          score = flux * volume *
                  data::mg.macro_xs_[material].get_xs(
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

        // Apply score to the appropriate tally bin
        Tally& tally {*model::tallies[task.tally_idx]};
#pragma omp atomic
        tally.results_(task.filter_idx, task.score_idx, TallyResult::VALUE) +=
          score;
      }
    }

    // For flux tallies, the total volume of the spatial region is needed
    // for normalizing the flux. We store this volume in a separate tensor.
    // We only contribute to each volume tally bin once per FSR.
    if (volume_normalized_flux_tallies_) {
      for (const auto& task : volume_task_[sr]) {
        if (task.score_type == SCORE_FLUX) {
#pragma omp atomic
          tally_volumes_[task.tally_idx](task.filter_idx, task.score_idx) +=
            volume;
        }
      }
    }
  } // end FSR loop

  // Normalize any flux scores by the total volume of the FSRs scoring to that
  // bin. To do this, we loop over all tallies, and then all filter bins,
  // and then scores. For each score, we check the tally data structure to
  // see what index that score corresponds to. If that score is a flux score,
  // then we divide it by volume.
  if (volume_normalized_flux_tallies_) {
    for (int i = 0; i < model::tallies.size(); i++) {
      Tally& tally {*model::tallies[i]};
#pragma omp parallel for
      for (int bin = 0; bin < tally.n_filter_bins(); bin++) {
        for (int score_idx = 0; score_idx < tally.n_scores(); score_idx++) {
          auto score_type = tally.scores_[score_idx];
          if (score_type == SCORE_FLUX) {
            double vol = tally_volumes_[i](bin, score_idx);
            if (vol > 0.0) {
              tally.results_(bin, score_idx, TallyResult::VALUE) /= vol;
            }
          }
        }
      }
    }
  }

  openmc::simulation::time_tallies.stop();
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

  MPI_Allreduce(MPI_IN_PLACE, scalar_flux_new_.data(), n_source_elements_,
    MPI_DOUBLE, MPI_SUM, mpi::intracomm);

  simulation::time_bank_sendrecv.stop();
#endif
}

double FlatSourceDomain::evaluate_flux_at_point(
  Position r, int64_t sr, int g) const
{
  return scalar_flux_final_[sr * negroups_ + g] /
         (settings::n_batches - settings::n_inactive);
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
    vector<Position> voxel_positions(Nx * Ny * Nz);

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
          voxel_positions[z * Ny * Nx + y * Nx + x] = sample;
        }
      }
    }

    double source_normalization_factor =
      compute_fixed_source_normalization_factor();

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
      for (int i = 0; i < Nx * Ny * Nz; i++) {
        int64_t fsr = voxel_indices[i];
        int64_t source_element = fsr * negroups_ + g;
        float flux = evaluate_flux_at_point(voxel_positions[i], fsr, g);
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
    for (int i = 0; i < Nx * Ny * Nz; i++) {
      int64_t fsr = voxel_indices[i];

      float total_fission = 0.0;
      int mat = material_[fsr];
      for (int g = 0; g < negroups_; g++) {
        int64_t source_element = fsr * negroups_ + g;
        float flux = evaluate_flux_at_point(voxel_positions[i], fsr, g);
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

void FlatSourceDomain::apply_external_source_to_source_region(
  Discrete* discrete, double strength_factor, int64_t source_region)
{
  external_source_present_[source_region] = true;

  const auto& discrete_energies = discrete->x();
  const auto& discrete_probs = discrete->prob();

  for (int e = 0; e < discrete_energies.size(); e++) {
    int g = data::mg.get_group_index(discrete_energies[e]);
    external_source_[source_region * negroups_ + g] +=
      discrete_probs[e] * strength_factor;
  }
}

void FlatSourceDomain::apply_external_source_to_cell_instances(int32_t i_cell,
  Discrete* discrete, double strength_factor, int target_material_id,
  const vector<int32_t>& instances)
{
  Cell& cell = *model::cells[i_cell];

  if (cell.type_ != Fill::MATERIAL)
    return;

  for (int j : instances) {
    int cell_material_idx = cell.material(j);
    int cell_material_id = model::materials[cell_material_idx]->id();
    if (target_material_id == C_NONE ||
        cell_material_id == target_material_id) {
      int64_t source_region = source_region_offsets_[i_cell] + j;
      apply_external_source_to_source_region(
        discrete, strength_factor, source_region);
    }
  }
}

void FlatSourceDomain::apply_external_source_to_cell_and_children(
  int32_t i_cell, Discrete* discrete, double strength_factor,
  int32_t target_material_id)
{
  Cell& cell = *model::cells[i_cell];

  if (cell.type_ == Fill::MATERIAL) {
    vector<int> instances(cell.n_instances_);
    std::iota(instances.begin(), instances.end(), 0);
    apply_external_source_to_cell_instances(
      i_cell, discrete, strength_factor, target_material_id, instances);
  } else if (target_material_id == C_NONE) {
    std::unordered_map<int32_t, vector<int32_t>> cell_instance_list =
      cell.get_contained_cells(0, nullptr);
    for (const auto& pair : cell_instance_list) {
      int32_t i_child_cell = pair.first;
      apply_external_source_to_cell_instances(i_child_cell, discrete,
        strength_factor, target_material_id, pair.second);
    }
  }
}

void FlatSourceDomain::count_external_source_regions()
{
  n_external_source_regions_ = 0;
#pragma omp parallel for reduction(+ : n_external_source_regions_)
  for (int sr = 0; sr < n_source_regions_; sr++) {
    if (external_source_present_[sr]) {
      n_external_source_regions_++;
    }
  }
}

void FlatSourceDomain::convert_external_sources()
{
  // Loop over external sources
  for (int es = 0; es < model::external_sources.size(); es++) {
    Source* s = model::external_sources[es].get();
    IndependentSource* is = dynamic_cast<IndependentSource*>(s);
    Discrete* energy = dynamic_cast<Discrete*>(is->energy());
    const std::unordered_set<int32_t>& domain_ids = is->domain_ids();

    double strength_factor = is->strength();

    if (is->domain_type() == Source::DomainType::MATERIAL) {
      for (int32_t material_id : domain_ids) {
        for (int i_cell = 0; i_cell < model::cells.size(); i_cell++) {
          apply_external_source_to_cell_and_children(
            i_cell, energy, strength_factor, material_id);
        }
      }
    } else if (is->domain_type() == Source::DomainType::CELL) {
      for (int32_t cell_id : domain_ids) {
        int32_t i_cell = model::cell_map[cell_id];
        apply_external_source_to_cell_and_children(
          i_cell, energy, strength_factor, C_NONE);
      }
    } else if (is->domain_type() == Source::DomainType::UNIVERSE) {
      for (int32_t universe_id : domain_ids) {
        int32_t i_universe = model::universe_map[universe_id];
        Universe& universe = *model::universes[i_universe];
        for (int32_t i_cell : universe.cells_) {
          apply_external_source_to_cell_and_children(
            i_cell, energy, strength_factor, C_NONE);
        }
      }
    }
  } // End loop over external sources

  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single angle data.
  const int t = 0;
  const int a = 0;

// Divide the fixed source term by sigma t (to save time when applying each
// iteration)
#pragma omp parallel for
  for (int sr = 0; sr < n_source_regions_; sr++) {
    int material = material_[sr];
    for (int e = 0; e < negroups_; e++) {
      float sigma_t = data::mg.macro_xs_[material].get_xs(
        MgxsType::TOTAL, e, nullptr, nullptr, nullptr, t, a);
      external_source_[sr * negroups_ + e] /= sigma_t;
    }
  }
}
void FlatSourceDomain::flux_swap()
{
  scalar_flux_old_.swap(scalar_flux_new_);
}

} // namespace openmc
