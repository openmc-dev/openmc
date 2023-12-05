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

// This function is responsible for generating a mapping between random
// ray flat source regions (cell instances) and tally bins. The mapping
// takes the form of a "TallyTask" object, which accounts for one single
// score being applied to a single tally. Thus, a single source region
// may have anywhere from zero to many tally tasks associated with it --
// meaning that the global "tally_task" data structure is in 2D. The outer
// dimension corresponds to the source element (i.e., each entry corresponds
// to a specific energy group within a specific source region), and the
// inner dimension corresponds to the tallying task itself. Mechanically,
// the mapping between FSRs and spatial filters is done by considering
// the location of a single known ray midpoint that passed through the
// FSR. I.e., during transport, the first ray to pass through a given FSR
// will write down its midpoint for use with this function. This is a cheap
// and easy way of mapping FSrs to spatial tally filters, but comes with
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

bool convert_source_regions_to_tallies()
{
  int negroups = data::mg.num_energy_groups_;
  
  openmc::simulation::time_tallies.start();

  // Tracks if we've generated a mapping yet for all source regions.
  bool all_source_regions_mapped = true;

  // Attempt to generate mapping for all source regions
  #pragma omp parallel for
  for (int sr = 0; sr < random_ray::n_source_regions; sr++) {

    // If this source region has not been hit by a ray yet, then
    // we aren't going to be able to map it, so skip it.
    if (!random_ray::position_recorded[sr]) {
      all_source_regions_mapped = false;
      continue;
    }

    // A particle located at the recorded midpoint of a ray
    // crossing through this source region is used to estabilish
    // the spatial location of the source region
    Particle p;
    p.r() = random_ray::position[sr];  
    bool found = exhaustive_find_cell(p);

    // Loop over energy groups (so as to support energy filters)
    for (int e = 0; e < negroups; e++) {

      // Set particle to the current energy
      p.g() = e;
      p.g_last() = e;
      p.E() = data::mg.energy_bin_avg_[p.g()];
      p.E_last() = p.E();

      int64_t source_element = sr * negroups + e;
      
      // If this task has already been populated, we don't need to do
      // it again.
      if (random_ray::tally_task[source_element].size() > 0) {
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
          for (auto score_index = 0; score_index < tally.scores_.size(); score_index++) {
            auto score_bin = tally.scores_[score_index];
            // If a valid tally, filter, and score cobination has been found,
            // then add it to the list of tally tasks for this source element.
            random_ray::tally_task[source_element].emplace_back(i_tally, filter_index, score_index, score_bin);
          }
        }
      }
      // Reset all the filter matches for the next tally event.
      for (auto& match : p.filter_matches())
        match.bins_present_ = false;
    }
  }
  openmc::simulation::time_tallies.stop();

  return all_source_regions_mapped;
}

// Tallying in random ray is not done directly during transport, rather,
// it is done only once after each power iteration. This is made possible
// by way of a mapping data structure that relates spatial source regions
// (FSRs) to tally/filter/score combinations. The mechanism by which the
// mapping is done (and the limitations incurred) is documented in the
// "convert_source_regions_to_tallies()" function comments above. The present
// tally function simply traverses the mapping data structure and executes
// the scoring operations to OpenMC's native tally result arrays.

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
}

} // namespace openmc
