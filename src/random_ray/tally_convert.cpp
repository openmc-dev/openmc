#include "openmc/output.h"
#include "openmc/geometry.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/timer.h"
#include "openmc/mgxs_interface.h"
#include "openmc/particle.h"
#include "openmc/random_ray/tally_convert.h"
#include "openmc/random_ray/source_region.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace random_ray {

std::vector<std::vector<TallyTask>> tally_task;

}

//==============================================================================
// Non-method functions
//==============================================================================

bool convert_source_regions_to_tallies()
{
  int negroups = data::mg.num_energy_groups_;


  bool all_source_regions_mapped = true;

  #pragma omp parallel for
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {

    // If this source region has not been hit by a ray yet, then
    // we aren't going to be able to map it, so skip it.
    if (!random_ray::position_recorded[sr]) {
      all_source_regions_mapped = false;
      continue;
    }

    Particle p;
    p.r() = random_ray::position[sr];  

    bool found = exhaustive_find_cell(p);
    for (int e = 0; e < negroups; e++) {
      p.g() = e;
      p.g_last() = e;
      p.E() = data::mg.energy_bin_avg_[p.g()];
      p.E_last() = p.E();

      int64_t source_element = sr * negroups + e;
      
      // If this task has already been populated, we can move on
      if (random_ray::tally_task[source_element].size() > 0) {
        continue;
      }

      for (auto i_tally : model::active_tallies) {
        Tally& tally {*model::tallies[i_tally]};

        // Initialize an iterator over valid filter bin combinations.  If there are
        // no valid combinations, use a continue statement to ensure we skip the
        // assume_separate break below.
        auto filter_iter = FilterBinIter(tally, p);
        auto end = FilterBinIter(tally, true, &p.filter_matches());
        if (filter_iter == end)
          continue;

        // Loop over filter bins.
        for (; filter_iter != end; ++filter_iter) {
          auto filter_index = filter_iter.index_;
          auto filter_weight = filter_iter.weight_;

          // Loop over scores
          for (auto score_index = 0; score_index < tally.scores_.size(); score_index++) {
            auto score_bin = tally.scores_[score_index];
            random_ray::tally_task[source_element].emplace_back(i_tally, filter_index, score_index, score_bin);
          }
        }
      }
      // Reset all the filter matches for the next tally event.
      for (auto& match : p.filter_matches())
        match.bins_present_ = false;
    }
  }
  return all_source_regions_mapped;
}

void random_ray_tally()
{
  openmc::simulation::time_tallies.start();

  int negroups = data::mg.num_energy_groups_;
  
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
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {
    double volume = random_ray::volume[sr];
    double material = random_ray::material[sr];
    for (int e = 0; e < negroups; e++) {
      int idx = sr * negroups + e;
      double flux  = random_ray::scalar_flux_new[idx] * volume;
      for (auto& task : random_ray::tally_task[idx]) {
        double score;
        switch (task.score_type) {

          case SCORE_FLUX:
            score = flux;
            break;
          
          case SCORE_TOTAL:
            score = flux * data::mg.macro_xs_[material].get_xs(MgxsType::TOTAL, e, NULL, NULL, NULL, t, a);
            break;

          case SCORE_FISSION:
            score = flux * data::mg.macro_xs_[material].get_xs(MgxsType::FISSION, e, NULL, NULL, NULL, t, a);
            break;
          
          case SCORE_NU_FISSION:
            score = flux * data::mg.macro_xs_[material].get_xs(MgxsType::NU_FISSION, e, NULL, NULL, NULL, t, a);
            break;
          
          case SCORE_EVENTS:
            score = 1.0;
            break;
    
          default:
            fatal_error("Invalid score specified in tallies.xml. Only flux, total, fission, nu-fission, and events are supported in random ray mode.");
            break;
        }
        Tally& tally {*model::tallies[task.tally_idx]};
        #pragma omp atomic
        tally.results_(task.filter_idx, task.score_idx, TallyResult::VALUE) += score;
      }
    }
  }
  openmc::simulation::time_tallies.stop();
}

} // namespace openmc
