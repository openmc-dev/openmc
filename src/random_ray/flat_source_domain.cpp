#include "openmc/random_ray/flat_source_domain.h"

#include "openmc/cell.h"
#include "openmc/constants.h"
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
#include "openmc/weight_windows.h"

#include <cstdio>

namespace openmc {

//==============================================================================
// FlatSourceDomain implementation
//==============================================================================

// Static Variable Declarations
RandomRayVolumeEstimator FlatSourceDomain::volume_estimator_ {
  RandomRayVolumeEstimator::HYBRID};
bool FlatSourceDomain::volume_normalized_flux_tallies_ {false};
bool FlatSourceDomain::adjoint_ {false};
double FlatSourceDomain::diagonal_stabilization_rho_ {1.0};
std::unordered_map<int, vector<std::pair<Source::DomainType, int>>>
  FlatSourceDomain::mesh_domain_map_;

FlatSourceDomain::FlatSourceDomain() : negroups_(data::mg.num_energy_groups_)
{
  // Count the number of source regions, compute the cell offset
  // indices, and store the material type The reason for the offsets is that
  // some cell types may not have material fills, and therefore do not
  // produce FSRs. Thus, we cannot index into the global arrays directly
  int base_source_regions = 0;
  for (const auto& c : model::cells) {
    if (c->type_ != Fill::MATERIAL) {
      source_region_offsets_.push_back(-1);
    } else {
      source_region_offsets_.push_back(base_source_regions);
      base_source_regions += c->n_instances();
    }
  }

  // Initialize source regions.
  bool is_linear = RandomRay::source_shape_ != RandomRaySourceShape::FLAT;
  source_regions_ = SourceRegionContainer(negroups_, is_linear);
  source_regions_.assign(
    base_source_regions, SourceRegion(negroups_, is_linear));

  // Initialize materials
  int64_t source_region_id = 0;
  for (int i = 0; i < model::cells.size(); i++) {
    Cell& cell = *model::cells[i];
    if (cell.type_ == Fill::MATERIAL) {
      for (int j = 0; j < cell.n_instances(); j++) {
        source_regions_.material(source_region_id++) = cell.material(j);
      }
    }
  }

  // Sanity check
  if (source_region_id != base_source_regions) {
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
// Reset scalar fluxes and iteration volume tallies to zero
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    source_regions_.volume(sr) = 0.0;
    source_regions_.volume_sq(sr) = 0.0;
  }

#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements(); se++) {
    source_regions_.scalar_flux_new(se) = 0.0;
  }
}

void FlatSourceDomain::accumulate_iteration_flux()
{
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements(); se++) {
    source_regions_.scalar_flux_final(se) +=
      source_regions_.scalar_flux_new(se);
  }
}

// Compute new estimate of scattering + fission sources in each source region
// based on the flux estimate from the previous iteration.
void FlatSourceDomain::update_neutron_source(double k_eff)
{
  simulation::time_update_src.start();

  double inverse_k_eff = 1.0 / k_eff;

// Reset all source regions to zero (important for void regions)
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements(); se++) {
    source_regions_.source(se) = 0.0;
  }

  // Add scattering + fission source
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    int material = source_regions_.material(sr);
    if (material == MATERIAL_VOID) {
      continue;
    }
    for (int g_out = 0; g_out < negroups_; g_out++) {
      double sigma_t = sigma_t_[material * negroups_ + g_out];
      double scatter_source = 0.0;
      double fission_source = 0.0;

      for (int g_in = 0; g_in < negroups_; g_in++) {
        double scalar_flux = source_regions_.scalar_flux_old(sr, g_in);
        double sigma_s =
          sigma_s_[material * negroups_ * negroups_ + g_out * negroups_ + g_in];
        double nu_sigma_f = nu_sigma_f_[material * negroups_ + g_in];
        double chi = chi_[material * negroups_ + g_out];

        scatter_source += sigma_s * scalar_flux;
        fission_source += nu_sigma_f * scalar_flux * chi;
      }
      source_regions_.source(sr, g_out) =
        (scatter_source + fission_source * inverse_k_eff) / sigma_t;
    }
  }

  // Add external source if in fixed source mode
  if (settings::run_mode == RunMode::FIXED_SOURCE) {
#pragma omp parallel for
    for (int64_t se = 0; se < n_source_elements(); se++) {
      source_regions_.source(se) += source_regions_.external_source(se);
    }
  }

  simulation::time_update_src.stop();
}

// Normalizes flux and updates simulation-averaged volume estimate
void FlatSourceDomain::normalize_scalar_flux_and_volumes(
  double total_active_distance_per_iteration)
{
  double normalization_factor = 1.0 / total_active_distance_per_iteration;
  double volume_normalization_factor =
    1.0 / (total_active_distance_per_iteration * simulation::current_batch);

// Normalize scalar flux to total distance travelled by all rays this
// iteration
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements(); se++) {
    source_regions_.scalar_flux_new(se) *= normalization_factor;
  }

// Accumulate cell-wise ray length tallies collected this iteration, then
// update the simulation-averaged cell-wise volume estimates
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    source_regions_.volume_t(sr) += source_regions_.volume(sr);
    source_regions_.volume_sq_t(sr) += source_regions_.volume_sq(sr);
    source_regions_.volume_naive(sr) =
      source_regions_.volume(sr) * normalization_factor;
    source_regions_.volume_sq(sr) =
      source_regions_.volume_sq_t(sr) / source_regions_.volume_t(sr);
    source_regions_.volume(sr) =
      source_regions_.volume_t(sr) * volume_normalization_factor;
  }
}

void FlatSourceDomain::set_flux_to_flux_plus_source(
  int64_t sr, double volume, int g)
{
  int material = source_regions_.material(sr);
  if (material == MATERIAL_VOID) {
    source_regions_.scalar_flux_new(sr, g) /= volume;
    if (settings::run_mode == RunMode::FIXED_SOURCE) {
      source_regions_.scalar_flux_new(sr, g) +=
        0.5f * source_regions_.external_source(sr, g) *
        source_regions_.volume_sq(sr);
    }
  } else {
    double sigma_t = sigma_t_[source_regions_.material(sr) * negroups_ + g];
    source_regions_.scalar_flux_new(sr, g) /= (sigma_t * volume);
    source_regions_.scalar_flux_new(sr, g) += source_regions_.source(sr, g);
  }
}

void FlatSourceDomain::set_flux_to_old_flux(int64_t sr, int g)
{
  source_regions_.scalar_flux_new(sr, g) =
    source_regions_.scalar_flux_old(sr, g);
}

void FlatSourceDomain::set_flux_to_source(int64_t sr, int g)
{
  source_regions_.scalar_flux_new(sr, g) = source_regions_.source(sr, g);
}

// Combine transport flux contributions and flat source contributions from the
// previous iteration to generate this iteration's estimate of scalar flux.
int64_t FlatSourceDomain::add_source_to_scalar_flux()
{
  int64_t n_hits = 0;
  double inverse_batch = 1.0 / simulation::current_batch;

#pragma omp parallel for reduction(+ : n_hits)
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {

    double volume_simulation_avg = source_regions_.volume(sr);
    double volume_iteration = source_regions_.volume_naive(sr);

    // Increment the number of hits if cell was hit this iteration
    if (volume_iteration) {
      n_hits++;
    }

    // Set the SR to small status if its expected number of hits
    // per iteration is less than 1.5
    if (source_regions_.n_hits(sr) * inverse_batch < MIN_HITS_PER_BATCH) {
      source_regions_.is_small(sr) = 1;
    } else {
      source_regions_.is_small(sr) = 0;
    }

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
      if (source_regions_.external_source_present(sr) ||
          source_regions_.is_small(sr)) {
        volume = volume_iteration;
      } else {
        volume = volume_simulation_avg;
      }
      break;
    default:
      fatal_error("Invalid volume estimator type");
    }

    for (int g = 0; g < negroups_; g++) {
      // There are three scenarios we need to consider:
      if (volume_iteration > 0.0) {
        // 1. If the FSR was hit this iteration, then the new flux is equal to
        // the flat source from the previous iteration plus the contributions
        // from rays passing through the source region (computed during the
        // transport sweep)
        set_flux_to_flux_plus_source(sr, volume, g);
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
        if (source_regions_.external_source_present(sr)) {
          set_flux_to_old_flux(sr, g);
        } else {
          set_flux_to_source(sr, g);
        }
      }
      // Halt if NaN implosion is detected
      if (!std::isfinite(source_regions_.scalar_flux_new(sr, g))) {
        fatal_error("A source region scalar flux is not finite. "
                    "This indicates a numerical instability in the "
                    "simulation. Consider increasing ray density or adjusting "
                    "the source region mesh.");
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

  // Vector for gathering fission source terms for Shannon entropy calculation
  vector<float> p(n_source_regions(), 0.0f);

#pragma omp parallel for reduction(+ : fission_rate_old, fission_rate_new)
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {

    // If simulation averaged volume is zero, don't include this cell
    double volume = source_regions_.volume(sr);
    if (volume == 0.0) {
      continue;
    }

    int material = source_regions_.material(sr);
    if (material == MATERIAL_VOID) {
      continue;
    }

    double sr_fission_source_old = 0;
    double sr_fission_source_new = 0;

    for (int g = 0; g < negroups_; g++) {
      double nu_sigma_f = nu_sigma_f_[material * negroups_ + g];
      sr_fission_source_old +=
        nu_sigma_f * source_regions_.scalar_flux_old(sr, g);
      sr_fission_source_new +=
        nu_sigma_f * source_regions_.scalar_flux_new(sr, g);
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
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
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

// It takes as an argument the starting index in the source region array,
// and it will operate from that index until the end of the array. This
// is useful as it can be called for both explicit user source regions or
// when a source region mesh is overlaid.

void FlatSourceDomain::convert_source_regions_to_tallies(int64_t start_sr_id)
{
  openmc::simulation::time_tallies.start();

  // Tracks if we've generated a mapping yet for all source regions.
  bool all_source_regions_mapped = true;

// Attempt to generate mapping for all source regions
#pragma omp parallel for
  for (int64_t sr = start_sr_id; sr < n_source_regions(); sr++) {

    // If this source region has not been hit by a ray yet, then
    // we aren't going to be able to map it, so skip it.
    if (!source_regions_.position_recorded(sr)) {
      all_source_regions_mapped = false;
      continue;
    }

    // A particle located at the recorded midpoint of a ray
    // crossing through this source region is used to estabilish
    // the spatial location of the source region
    Particle p;
    p.r() = source_regions_.position(sr);
    p.r_last() = source_regions_.position(sr);
    p.u() = {1.0, 0.0, 0.0};
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
      if (source_regions_.tally_task(sr, g).size() > 0) {
        continue;
      }

      // Loop over all active tallies. This logic is essentially identical
      // to what happens when scanning for applicable tallies during
      // MC transport.
      for (int i_tally = 0; i_tally < model::tallies.size(); i_tally++) {
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
          for (int score = 0; score < tally.scores_.size(); score++) {
            auto score_bin = tally.scores_[score];
            // If a valid tally, filter, and score combination has been found,
            // then add it to the list of tally tasks for this source element.
            TallyTask task(i_tally, filter_index, score, score_bin);
            source_regions_.tally_task(sr, g).push_back(task);

            // Also add this task to the list of volume tasks for this source
            // region.
            source_regions_.volume_task(sr).insert(task);
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
// volume estimate of each FSR, which improves over the course of the
// simulation
double FlatSourceDomain::compute_fixed_source_normalization_factor() const
{
  // If we are not in fixed source mode, then there are no external sources
  // so no normalization is needed.
  if (settings::run_mode != RunMode::FIXED_SOURCE || adjoint_) {
    return 1.0;
  }

  // Step 1 is to sum over all source regions and energy groups to get the
  // total external source strength in the simulation.
  double simulation_external_source_strength = 0.0;
#pragma omp parallel for reduction(+ : simulation_external_source_strength)
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    int material = source_regions_.material(sr);
    double volume = source_regions_.volume(sr) * simulation_volume_;
    for (int g = 0; g < negroups_; g++) {
      // For non-void regions, we store the external source pre-divided by
      // sigma_t. We need to multiply non-void regions back up by sigma_t
      // to get the total source strength in the expected units.
      double sigma_t = 1.0;
      if (material != MATERIAL_VOID) {
        sigma_t = sigma_t_[material * negroups_ + g];
      }
      simulation_external_source_strength +=
        source_regions_.external_source(sr, g) * sigma_t * volume;
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

  double source_normalization_factor =
    compute_fixed_source_normalization_factor();

// We loop over all source regions and energy groups. For each
// element, we check if there are any scores needed and apply
// them.
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    // The fsr.volume_ is the unitless fractional simulation averaged volume
    // (i.e., it is the FSR's fraction of the overall simulation volume). The
    // simulation_volume_ is the total 3D physical volume in cm^3 of the
    // entire global simulation domain (as defined by the ray source box).
    // Thus, the FSR's true 3D spatial volume in cm^3 is found by multiplying
    // its fraction of the total volume by the total volume. Not important in
    // eigenvalue solves, but useful in fixed source solves for returning the
    // flux shape with a magnitude that makes sense relative to the fixed
    // source strength.
    double volume = source_regions_.volume(sr) * simulation_volume_;

    double material = source_regions_.material(sr);
    for (int g = 0; g < negroups_; g++) {
      double flux =
        source_regions_.scalar_flux_new(sr, g) * source_normalization_factor;

      // Determine numerical score value
      for (auto& task : source_regions_.tally_task(sr, g)) {
        double score = 0.0;
        switch (task.score_type) {

        case SCORE_FLUX:
          score = flux * volume;
          break;

        case SCORE_TOTAL:
          if (material != MATERIAL_VOID) {
            score = flux * volume * sigma_t_[material * negroups_ + g];
          }
          break;

        case SCORE_FISSION:
          if (material != MATERIAL_VOID) {
            score = flux * volume * sigma_f_[material * negroups_ + g];
          }
          break;

        case SCORE_NU_FISSION:
          if (material != MATERIAL_VOID) {
            score = flux * volume * nu_sigma_f_[material * negroups_ + g];
          }
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
      for (const auto& task : source_regions_.volume_task(sr)) {
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

double FlatSourceDomain::evaluate_flux_at_point(
  Position r, int64_t sr, int g) const
{
  return source_regions_.scalar_flux_final(sr, g) /
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
    vector<double> weight_windows(Nx * Ny * Nz);
    float min_weight = 1e20;
#pragma omp parallel for collapse(3) reduction(min : min_weight)
    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
          Position sample;
          sample.z = ll.z + z_delta / 2.0 + z * z_delta;
          sample.y = ll.y + y_delta / 2.0 + y * y_delta;
          sample.x = ll.x + x_delta / 2.0 + x * x_delta;
          Particle p;
          p.r() = sample;
          p.r_last() = sample;
          p.E() = 1.0;
          p.E_last() = 1.0;
          p.u() = {1.0, 0.0, 0.0};

          bool found = exhaustive_find_cell(p);
          if (!found) {
            voxel_indices[z * Ny * Nx + y * Nx + x] = -1;
            voxel_positions[z * Ny * Nx + y * Nx + x] = sample;
            weight_windows[z * Ny * Nx + y * Nx + x] = 0.0;
            continue;
          }

          int i_cell = p.lowest_coord().cell;
          int64_t sr = source_region_offsets_[i_cell] + p.cell_instance();
          if (RandomRay::mesh_subdivision_enabled_) {
            int mesh_idx = base_source_regions_.mesh(sr);
            int mesh_bin;
            if (mesh_idx == C_NONE) {
              mesh_bin = 0;
            } else {
              mesh_bin = model::meshes[mesh_idx]->get_bin(p.r());
            }
            SourceRegionKey sr_key {sr, mesh_bin};
            auto it = source_region_map_.find(sr_key);
            if (it != source_region_map_.end()) {
              sr = it->second;
            } else {
              sr = -1;
            }
          }

          voxel_indices[z * Ny * Nx + y * Nx + x] = sr;
          voxel_positions[z * Ny * Nx + y * Nx + x] = sample;

          if (variance_reduction::weight_windows.size() == 1) {
            WeightWindow ww =
              variance_reduction::weight_windows[0]->get_weight_window(p);
            float weight = ww.lower_weight;
            weight_windows[z * Ny * Nx + y * Nx + x] = weight;
            if (weight < min_weight)
              min_weight = weight;
          }
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
    std::fprintf(plot, "ORIGIN %lf %lf %lf\n", ll.x, ll.y, ll.z);
    std::fprintf(plot, "SPACING %lf %lf %lf\n", x_delta, y_delta, z_delta);
    std::fprintf(plot, "POINT_DATA %d\n", Nx * Ny * Nz);

    int64_t num_neg = 0;
    int64_t num_samples = 0;
    float min_flux = 0.0;
    float max_flux = -1.0e20;
    // Plot multigroup flux data
    for (int g = 0; g < negroups_; g++) {
      std::fprintf(plot, "SCALARS flux_group_%d float\n", g);
      std::fprintf(plot, "LOOKUP_TABLE default\n");
      for (int i = 0; i < Nx * Ny * Nz; i++) {
        int64_t fsr = voxel_indices[i];
        int64_t source_element = fsr * negroups_ + g;
        float flux = 0;
        if (fsr >= 0) {
          flux = evaluate_flux_at_point(voxel_positions[i], fsr, g);
          if (flux < 0.0)
            flux = FlatSourceDomain::evaluate_flux_at_point(
              voxel_positions[i], fsr, g);
        }
        if (flux < 0.0) {
          num_neg++;
          if (flux < min_flux) {
            min_flux = flux;
          }
        }
        if (flux > max_flux)
          max_flux = flux;
        num_samples++;
        flux = convert_to_big_endian<float>(flux);
        std::fwrite(&flux, sizeof(float), 1, plot);
      }
    }

    // Slightly negative fluxes can be normal when sampling corners of linear
    // source regions. However, very common and high magnitude negative fluxes
    // may indicate numerical instability.
    if (num_neg > 0) {
      warning(fmt::format("{} plot samples ({:.4f}%) contained negative fluxes "
                          "(minumum found = {:.2e} maximum_found = {:.2e})",
        num_neg, (100.0 * num_neg) / num_samples, min_flux, max_flux));
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
      int mat = -1;
      if (fsr >= 0)
        mat = source_regions_.material(fsr);
      mat = convert_to_big_endian<int>(mat);
      std::fwrite(&mat, sizeof(int), 1, plot);
    }

    // Plot fission source
    if (settings::run_mode == RunMode::EIGENVALUE) {
      std::fprintf(plot, "SCALARS total_fission_source float\n");
      std::fprintf(plot, "LOOKUP_TABLE default\n");
      for (int i = 0; i < Nx * Ny * Nz; i++) {
        int64_t fsr = voxel_indices[i];
        float total_fission = 0.0;
        if (fsr >= 0) {
          int mat = source_regions_.material(fsr);
          if (mat != MATERIAL_VOID) {
            for (int g = 0; g < negroups_; g++) {
              int64_t source_element = fsr * negroups_ + g;
              float flux = evaluate_flux_at_point(voxel_positions[i], fsr, g);
              double sigma_f = sigma_f_[mat * negroups_ + g];
              total_fission += sigma_f * flux;
            }
          }
        }
        total_fission = convert_to_big_endian<float>(total_fission);
        std::fwrite(&total_fission, sizeof(float), 1, plot);
      }
    } else {
      std::fprintf(plot, "SCALARS external_source float\n");
      std::fprintf(plot, "LOOKUP_TABLE default\n");
      for (int i = 0; i < Nx * Ny * Nz; i++) {
        int64_t fsr = voxel_indices[i];
        int mat = source_regions_.material(fsr);
        float total_external = 0.0f;
        if (fsr >= 0) {
          for (int g = 0; g < negroups_; g++) {
            // External sources are already divided by sigma_t, so we need to
            // multiply it back to get the true external source.
            double sigma_t = 1.0;
            if (mat != MATERIAL_VOID) {
              sigma_t = sigma_t_[mat * negroups_ + g];
            }
            total_external += source_regions_.external_source(fsr, g) * sigma_t;
          }
        }
        total_external = convert_to_big_endian<float>(total_external);
        std::fwrite(&total_external, sizeof(float), 1, plot);
      }
    }

    // Plot weight window data
    if (variance_reduction::weight_windows.size() == 1) {
      std::fprintf(plot, "SCALARS weight_window_lower float\n");
      std::fprintf(plot, "LOOKUP_TABLE default\n");
      for (int i = 0; i < Nx * Ny * Nz; i++) {
        float weight = weight_windows[i];
        if (weight == 0.0)
          weight = min_weight;
        weight = convert_to_big_endian<float>(weight);
        std::fwrite(&weight, sizeof(float), 1, plot);
      }
    }

    std::fclose(plot);
  }
}

void FlatSourceDomain::apply_external_source_to_source_region(
  Discrete* discrete, double strength_factor, SourceRegionHandle& srh)
{
  srh.external_source_present() = 1;

  const auto& discrete_energies = discrete->x();
  const auto& discrete_probs = discrete->prob();

  for (int i = 0; i < discrete_energies.size(); i++) {
    int g = data::mg.get_group_index(discrete_energies[i]);
    srh.external_source(g) += discrete_probs[i] * strength_factor;
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
    int cell_material_id;
    if (cell_material_idx == MATERIAL_VOID) {
      cell_material_id = MATERIAL_VOID;
    } else {
      cell_material_id = model::materials[cell_material_idx]->id();
    }
    if (target_material_id == C_NONE ||
        cell_material_id == target_material_id) {
      int64_t source_region = source_region_offsets_[i_cell] + j;
      SourceRegionHandle srh =
        source_regions_.get_source_region_handle(source_region);
      apply_external_source_to_source_region(discrete, strength_factor, srh);
    }
  }
}

void FlatSourceDomain::apply_external_source_to_cell_and_children(
  int32_t i_cell, Discrete* discrete, double strength_factor,
  int32_t target_material_id)
{
  Cell& cell = *model::cells[i_cell];

  if (cell.type_ == Fill::MATERIAL) {
    vector<int> instances(cell.n_instances());
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
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    if (source_regions_.external_source_present(sr)) {
      n_external_source_regions_++;
    }
  }
}

void FlatSourceDomain::convert_external_sources()
{
  // Loop over external sources
  for (int es = 0; es < model::external_sources.size(); es++) {

    // Extract source information
    Source* s = model::external_sources[es].get();
    IndependentSource* is = dynamic_cast<IndependentSource*>(s);
    Discrete* energy = dynamic_cast<Discrete*>(is->energy());
    const std::unordered_set<int32_t>& domain_ids = is->domain_ids();
    double strength_factor = is->strength();

    // If there is no domain constraint specified, then this must be a point
    // source. In this case, we need to find the source region that contains the
    // point source and apply or relate it to the external source.
    if (is->domain_ids().size() == 0) {

      // Extract the point source coordinate and find the base source region at
      // that point
      auto sp = dynamic_cast<SpatialPoint*>(is->space());
      GeometryState gs;
      gs.r() = sp->r();
      gs.r_last() = sp->r();
      gs.u() = {1.0, 0.0, 0.0};
      bool found = exhaustive_find_cell(gs);
      if (!found) {
        fatal_error(fmt::format("Could not find cell containing external "
                                "point source at {}",
          sp->r()));
      }
      int i_cell = gs.lowest_coord().cell;
      int64_t sr = source_region_offsets_[i_cell] + gs.cell_instance();

      if (RandomRay::mesh_subdivision_enabled_) {
        // If mesh subdivision is enabled, we need to determine which subdivided
        // mesh bin the point source coordinate is in as well
        int mesh_idx = source_regions_.mesh(sr);
        int mesh_bin;
        if (mesh_idx == C_NONE) {
          mesh_bin = 0;
        } else {
          mesh_bin = model::meshes[mesh_idx]->get_bin(gs.r());
        }
        // With the source region and mesh bin known, we can use the
        // accompanying SourceRegionKey as a key into a map that stores the
        // corresponding external source index for the point source. Notably, we
        // do not actually apply the external source to any source regions here,
        // as if mesh subdivision is enabled, they haven't actually been
        // discovered & initilized yet. When discovered, they will read from the
        // point_source_map to determine if there are any point source terms
        // that should be applied.
        SourceRegionKey key {sr, mesh_bin};
        point_source_map_[key] = es;
      } else {
        // If we are not using mesh subdivision, we can apply the external
        // source directly to the source region as we do for volumetric domain
        // constraint sources.
        SourceRegionHandle srh = source_regions_.get_source_region_handle(sr);
        apply_external_source_to_source_region(energy, strength_factor, srh);
      }

    } else {
      // If not a point source, then use the volumetric domain constraints to
      // determine which source regions to apply the external source to.
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
    }
  } // End loop over external sources

// Divide the fixed source term by sigma t (to save time when applying each
// iteration)
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    int material = source_regions_.material(sr);
    if (material == MATERIAL_VOID) {
      continue;
    }
    for (int g = 0; g < negroups_; g++) {
      double sigma_t = sigma_t_[material * negroups_ + g];
      source_regions_.external_source(sr, g) /= sigma_t;
    }
  }
}

void FlatSourceDomain::flux_swap()
{
  source_regions_.flux_swap();
}

void FlatSourceDomain::flatten_xs()
{
  // Temperature and angle indices, if using multiple temperature
  // data sets and/or anisotropic data sets.
  // TODO: Currently assumes we are only using single temp/single angle data.
  const int t = 0;
  const int a = 0;

  n_materials_ = data::mg.macro_xs_.size();
  for (auto& m : data::mg.macro_xs_) {
    for (int g_out = 0; g_out < negroups_; g_out++) {
      if (m.exists_in_model) {
        double sigma_t =
          m.get_xs(MgxsType::TOTAL, g_out, NULL, NULL, NULL, t, a);
        sigma_t_.push_back(sigma_t);

        double nu_sigma_f =
          m.get_xs(MgxsType::NU_FISSION, g_out, NULL, NULL, NULL, t, a);
        nu_sigma_f_.push_back(nu_sigma_f);

        double sigma_f =
          m.get_xs(MgxsType::FISSION, g_out, NULL, NULL, NULL, t, a);
        sigma_f_.push_back(sigma_f);

        double chi =
          m.get_xs(MgxsType::CHI_PROMPT, g_out, &g_out, NULL, NULL, t, a);
        if (!std::isfinite(chi)) {
          // MGXS interface may return NaN in some cases, such as when material
          // is fissionable but has very small sigma_f.
          chi = 0.0;
        }
        chi_.push_back(chi);

        for (int g_in = 0; g_in < negroups_; g_in++) {
          double sigma_s =
            m.get_xs(MgxsType::NU_SCATTER, g_in, &g_out, NULL, NULL, t, a);
          sigma_s_.push_back(sigma_s);
          // For transport corrected XS data, diagonal elements may be negative.
          // In this case, set a flag to enable transport stabilization for the
          // simulation.
          if (g_out == g_in && sigma_s < 0.0)
            is_transport_stabilization_needed_ = true;
        }
      } else {
        sigma_t_.push_back(0);
        nu_sigma_f_.push_back(0);
        sigma_f_.push_back(0);
        chi_.push_back(0);
        for (int g_in = 0; g_in < negroups_; g_in++) {
          sigma_s_.push_back(0);
        }
      }
    }
  }
}

void FlatSourceDomain::set_adjoint_sources(const vector<double>& forward_flux)
{
  // Set the adjoint external source to 1/forward_flux. If the forward flux is
  // negative, zero, or extremely close to zero, set the adjoint source to zero,
  // as this is likely a very small source region that we don't need to bother
  // trying to vector particles towards. In the case of flux "being extremely
  // close to zero", we define this as being a fixed fraction of the maximum
  // forward flux, below which we assume the flux would be physically
  // undetectable.

  // First, find the maximum forward flux value
  double max_flux = 0.0;
#pragma omp parallel for reduction(max : max_flux)
  for (int64_t se = 0; se < n_source_elements(); se++) {
    double flux = forward_flux[se];
    if (flux > max_flux) {
      max_flux = flux;
    }
  }

  // Then, compute the adjoint source for each source region
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    for (int g = 0; g < negroups_; g++) {
      double flux = forward_flux[sr * negroups_ + g];
      if (flux <= ZERO_FLUX_CUTOFF * max_flux) {
        source_regions_.external_source(sr, g) = 0.0;
      } else {
        source_regions_.external_source(sr, g) = 1.0 / flux;
      }
      if (flux > 0.0) {
        source_regions_.external_source_present(sr) = 1;
      }
    }
  }

  // Divide the fixed source term by sigma t (to save time when applying each
  // iteration)
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    int material = source_regions_.material(sr);
    if (material == MATERIAL_VOID) {
      continue;
    }
    for (int g = 0; g < negroups_; g++) {
      double sigma_t = sigma_t_[material * negroups_ + g];
      source_regions_.external_source(sr, g) /= sigma_t;
    }
  }
}

void FlatSourceDomain::transpose_scattering_matrix()
{
  // Transpose the inner two dimensions for each material
  for (int m = 0; m < n_materials_; ++m) {
    int material_offset = m * negroups_ * negroups_;
    for (int i = 0; i < negroups_; ++i) {
      for (int j = i + 1; j < negroups_; ++j) {
        // Calculate indices of the elements to swap
        int idx1 = material_offset + i * negroups_ + j;
        int idx2 = material_offset + j * negroups_ + i;

        // Swap the elements to transpose the matrix
        std::swap(sigma_s_[idx1], sigma_s_[idx2]);
      }
    }
  }
}

void FlatSourceDomain::serialize_final_fluxes(vector<double>& flux)
{
  // Ensure array is correct size
  flux.resize(n_source_regions() * negroups_);
// Serialize the final fluxes for output
#pragma omp parallel for
  for (int64_t se = 0; se < n_source_elements(); se++) {
    flux[se] = source_regions_.scalar_flux_final(se);
  }
}

void FlatSourceDomain::apply_mesh_to_cell_instances(int32_t i_cell,
  int32_t mesh_idx, int target_material_id, const vector<int32_t>& instances,
  bool is_target_void)
{
  Cell& cell = *model::cells[i_cell];
  if (cell.type_ != Fill::MATERIAL)
    return;
  for (int32_t j : instances) {
    int cell_material_idx = cell.material(j);
    int cell_material_id = (cell_material_idx == C_NONE)
                             ? C_NONE
                             : model::materials[cell_material_idx]->id();

    if ((target_material_id == C_NONE && !is_target_void) ||
        cell_material_id == target_material_id) {
      int64_t sr = source_region_offsets_[i_cell] + j;
      if (source_regions_.mesh(sr) != C_NONE) {
        // print out the source region that is broken:
        fatal_error(fmt::format("Source region {} already has mesh idx {} "
                                "applied, but trying to apply mesh idx {}",
          sr, source_regions_.mesh(sr), mesh_idx));
      }
      source_regions_.mesh(sr) = mesh_idx;
    }
  }
}

void FlatSourceDomain::apply_mesh_to_cell_and_children(int32_t i_cell,
  int32_t mesh_idx, int32_t target_material_id, bool is_target_void)
{
  Cell& cell = *model::cells[i_cell];

  if (cell.type_ == Fill::MATERIAL) {
    vector<int> instances(cell.n_instances());
    std::iota(instances.begin(), instances.end(), 0);
    apply_mesh_to_cell_instances(
      i_cell, mesh_idx, target_material_id, instances, is_target_void);
  } else if (target_material_id == C_NONE && !is_target_void) {
    for (int j = 0; j < cell.n_instances(); j++) {
      std::unordered_map<int32_t, vector<int32_t>> cell_instance_list =
        cell.get_contained_cells(j, nullptr);
      for (const auto& pair : cell_instance_list) {
        int32_t i_child_cell = pair.first;
        apply_mesh_to_cell_instances(i_child_cell, mesh_idx, target_material_id,
          pair.second, is_target_void);
      }
    }
  }
}

void FlatSourceDomain::apply_meshes()
{
  // Skip if there are no mappings between mesh IDs and domains
  if (mesh_domain_map_.empty())
    return;

  // Loop over meshes
  for (int mesh_idx = 0; mesh_idx < model::meshes.size(); mesh_idx++) {
    Mesh* mesh = model::meshes[mesh_idx].get();
    int mesh_id = mesh->id();

    // Skip if mesh id is not present in the map
    if (mesh_domain_map_.find(mesh_id) == mesh_domain_map_.end())
      continue;

    // Loop over domains associated with the mesh
    for (auto& domain : mesh_domain_map_[mesh_id]) {
      Source::DomainType domain_type = domain.first;
      int domain_id = domain.second;

      if (domain_type == Source::DomainType::MATERIAL) {
        for (int i_cell = 0; i_cell < model::cells.size(); i_cell++) {
          if (domain_id == C_NONE) {
            apply_mesh_to_cell_and_children(i_cell, mesh_idx, domain_id, true);
          } else {
            apply_mesh_to_cell_and_children(i_cell, mesh_idx, domain_id, false);
          }
        }
      } else if (domain_type == Source::DomainType::CELL) {
        int32_t i_cell = model::cell_map[domain_id];
        apply_mesh_to_cell_and_children(i_cell, mesh_idx, C_NONE, false);
      } else if (domain_type == Source::DomainType::UNIVERSE) {
        int32_t i_universe = model::universe_map[domain_id];
        Universe& universe = *model::universes[i_universe];
        for (int32_t i_cell : universe.cells_) {
          apply_mesh_to_cell_and_children(i_cell, mesh_idx, C_NONE, false);
        }
      }
    }
  }
}

void FlatSourceDomain::prepare_base_source_regions()
{
  std::swap(source_regions_, base_source_regions_);
  source_regions_.negroups() = base_source_regions_.negroups();
  source_regions_.is_linear() = base_source_regions_.is_linear();
}

SourceRegionHandle FlatSourceDomain::get_subdivided_source_region_handle(
  int64_t sr, int mesh_bin, Position r, double dist, Direction u)
{
  SourceRegionKey sr_key {sr, mesh_bin};

  // Case 1: Check if the source region key is already present in the permanent
  // map. This is the most common condition, as any source region visited in a
  // previous power iteration will already be present in the permanent map. If
  // the source region key is found, we translate the key into a specific 1D
  // source region index and return a handle its position in the
  // source_regions_ vector.
  auto it = source_region_map_.find(sr_key);
  if (it != source_region_map_.end()) {
    int64_t sr = it->second;
    return source_regions_.get_source_region_handle(sr);
  }

  // Case 2: Check if the source region key is present in the temporary (thread
  // safe) map. This is a common occurrence in the first power iteration when
  // the source region has already been visited already by some other ray. We
  // begin by locking the temporary map before any operations are performed. The
  // lock is not global over the full data structure -- it will be dependent on
  // which key is used.
  discovered_source_regions_.lock(sr_key);

  // If the key is found in the temporary map, then we return a handle to the
  // source region that is stored in the temporary map.
  if (discovered_source_regions_.contains(sr_key)) {
    SourceRegionHandle handle {discovered_source_regions_[sr_key]};
    discovered_source_regions_.unlock(sr_key);
    return handle;
  }

  // Case 3: The source region key is not present anywhere, but it is only due
  // to floating point artifacts. These artifacts occur when the overlaid mesh
  // overlaps with actual geometry surfaces. In these cases, roundoff error may
  // result in the ray tracer detecting an additional (very short) segment
  // though a mesh bin that is actually past the physical source region
  // boundary. This is a result of the the multi-level ray tracing treatment in
  // OpenMC, which depending on the number of universes in the hierarchy etc can
  // result in the wrong surface being selected as the nearest. This can happen
  // in a lattice when there are two directions that both are very close in
  // distance, within the tolerance of FP_REL_PRECISION, and the are thus
  // treated as being equivalent so alternative logic is used. However, when we
  // go and ray trace on this with the mesh tracer we may go past the surface
  // bounding the current source region.
  //
  // To filter out this case, before we create the new source region, we double
  // check that the actual starting point of this segment (r) is still in the
  // same geometry source region that we started in. If an artifact is detected,
  // we discard the segment (and attenuation through it) as it is not really a
  // valid source region and will have only an infinitessimally small cell
  // combined with the mesh bin. Thankfully, this is a fairly rare condition,
  // and only triggers for very short ray lengths. It can be fixed by decreasing
  // the value of FP_REL_PRECISION in constants.h, but this may have unknown
  // consequences for the general ray tracer, so for now we do the below sanity
  // checks before generating phantom source regions. A significant extra cost
  // is incurred in instantiating the GeometryState object and doing a cell
  // lookup, but again, this is going to be an extremely rare thing to check
  // after the first power iteration has completed.

  // Sanity check on source region id
  GeometryState gs;
  gs.r() = r + TINY_BIT * u;
  gs.u() = {1.0, 0.0, 0.0};
  exhaustive_find_cell(gs);
  int gs_i_cell = gs.lowest_coord().cell;
  int64_t sr_found = source_region_offsets_[gs_i_cell] + gs.cell_instance();
  if (sr_found != sr) {
    discovered_source_regions_.unlock(sr_key);
    SourceRegionHandle handle;
    handle.is_numerical_fp_artifact_ = true;
    return handle;
  }

  // Sanity check on mesh bin
  int mesh_idx = base_source_regions_.mesh(sr);
  if (mesh_idx == C_NONE) {
    if (mesh_bin != 0) {
      discovered_source_regions_.unlock(sr_key);
      SourceRegionHandle handle;
      handle.is_numerical_fp_artifact_ = true;
      return handle;
    }
  } else {
    Mesh* mesh = model::meshes[mesh_idx].get();
    int bin_found = mesh->get_bin(r + TINY_BIT * u);
    if (bin_found != mesh_bin) {
      discovered_source_regions_.unlock(sr_key);
      SourceRegionHandle handle;
      handle.is_numerical_fp_artifact_ = true;
      return handle;
    }
  }

  // Case 4: The source region key is valid, but is not present anywhere. This
  // condition only occurs the first time the source region is discovered
  // (typically in the first power iteration). In this case, we need to handle
  // creation of the new source region and its storage into the parallel map.
  // The new source region is created by copying the base source region, so as
  // to inherit material, external source, and some flux properties etc. We
  // also pass the base source region id to allow the new source region to
  // know which base source region it is derived from.
  SourceRegion* sr_ptr = discovered_source_regions_.emplace(
    sr_key, {base_source_regions_.get_source_region_handle(sr), sr});
  discovered_source_regions_.unlock(sr_key);
  SourceRegionHandle handle {*sr_ptr};

  // Check if the new source region contains a point source and apply it if so
  auto it2 = point_source_map_.find(sr_key);
  if (it2 != point_source_map_.end()) {
    int es = it2->second;
    auto s = model::external_sources[es].get();
    auto is = dynamic_cast<IndependentSource*>(s);
    auto energy = dynamic_cast<Discrete*>(is->energy());
    double strength_factor = is->strength();
    apply_external_source_to_source_region(energy, strength_factor, handle);
    int material = handle.material();
    if (material != MATERIAL_VOID) {
      for (int g = 0; g < negroups_; g++) {
        double sigma_t = sigma_t_[material * negroups_ + g];
        handle.external_source(g) /= sigma_t;
      }
    }
  }

  return handle;
}

void FlatSourceDomain::finalize_discovered_source_regions()
{
  // Extract keys for entries with a valid volume.
  vector<SourceRegionKey> keys;
  for (const auto& pair : discovered_source_regions_) {
    if (pair.second.volume_ > 0.0) {
      keys.push_back(pair.first);
    }
  }

  if (!keys.empty()) {
    // Sort the keys, so as to ensure reproducible ordering given that source
    // regions may have been added to discovered_source_regions_ in an arbitrary
    // order due to shared memory threading.
    std::sort(keys.begin(), keys.end());

    // Remember the index of the first new source region
    int64_t start_sr_id = source_regions_.n_source_regions();

    // Append the source regions in the sorted key order.
    for (const auto& key : keys) {
      const SourceRegion& sr = discovered_source_regions_[key];
      source_region_map_[key] = source_regions_.n_source_regions();
      source_regions_.push_back(sr);
    }

    // Map all new source regions to tallies
    convert_source_regions_to_tallies(start_sr_id);
  }

  discovered_source_regions_.clear();
}

// This is the "diagonal stabilization" technique developed by Gunow et al. in:
//
// Geoffrey Gunow, Benoit Forget, Kord Smith, Stabilization of multi-group
// neutron transport with transport-corrected cross-sections, Annals of Nuclear
// Energy, Volume 126, 2019, Pages 211-219, ISSN 0306-4549,
// https://doi.org/10.1016/j.anucene.2018.10.036.
void FlatSourceDomain::apply_transport_stabilization()
{
  // Don't do anything if all in-group scattering
  // cross sections are positive
  if (!is_transport_stabilization_needed_) {
    return;
  }

  // Apply the stabilization factor to all source elements
#pragma omp parallel for
  for (int64_t sr = 0; sr < n_source_regions(); sr++) {
    int material = source_regions_.material(sr);
    if (material == MATERIAL_VOID) {
      continue;
    }
    for (int g = 0; g < negroups_; g++) {
      // Only apply stabilization if the diagonal (in-group) scattering XS is
      // negative
      double sigma_s =
        sigma_s_[material * negroups_ * negroups_ + g * negroups_ + g];
      if (sigma_s < 0.0) {
        double sigma_t = sigma_t_[material * negroups_ + g];
        double phi_new = source_regions_.scalar_flux_new(sr, g);
        double phi_old = source_regions_.scalar_flux_old(sr, g);

        // Equation 18 in the above Gunow et al. 2019 paper. For a default
        // rho of 1.0, this ensures there are no negative diagonal elements
        // in the iteration matrix. A lesser rho could be used (or exposed
        // as a user input parameter) to reduce the negative impact on
        // convergence rate though would need to be experimentally tested to see
        // if it doesn't become unstable. rho = 1.0 is good as it gives the
        // highest assurance of stability, and the impacts on convergence rate
        // are pretty mild.
        double D = diagonal_stabilization_rho_ * sigma_s / sigma_t;

        // Equation 16 in the above Gunow et al. 2019 paper
        source_regions_.scalar_flux_new(sr, g) =
          (phi_new - D * phi_old) / (1.0 - D);
      }
    }
  }
}

} // namespace openmc
