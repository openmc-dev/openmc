#include "openmc/tallies/tally.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_delayedgroup.h"
#include "openmc/tallies/filter_surface.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/xml_interface.h"

#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp" // for empty_like
#include "xtensor/xview.hpp"

#include <array>
#include <cstddef>
#include <sstream>
#include <string>

namespace openmc {

//==============================================================================
// Functions defined in Fortran
//==============================================================================

extern "C" void
score_general_ce(Particle* p, int i_tally, int start_index, int filter_index,
  int i_nuclide, double atom_density, double flux);

extern "C" void
score_general_mg(Particle* p, int i_tally, int start_index, int filter_index,
  int i_nuclide, double atom_density, double flux);

extern "C" void
score_all_nuclides(Particle* p, int i_tally, double flux, int filter_index);

extern "C" int material_nuclide_index(int i_material, int i_nuclide);

extern "C" double material_atom_density(int i_material, int i);

extern "C" int tally_get_n_score_bins(int i_tally);

//==============================================================================
// Global variable definitions
//==============================================================================

namespace model {
  std::vector<std::unique_ptr<Tally>> tallies;

  std::vector<int> active_tallies;
  std::vector<int> active_analog_tallies;
  std::vector<int> active_tracklength_tallies;
  std::vector<int> active_collision_tallies;
  std::vector<int> active_meshsurf_tallies;
  std::vector<int> active_surface_tallies;
}

double global_tally_absorption;
double global_tally_collision;
double global_tally_tracklength;
double global_tally_leakage;

//==============================================================================
// Tally object implementation
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
    if (dynamic_cast<const EnergyoutFilter*>(filt)) {
      energyout_filter_ = i + 1;
    } else if (dynamic_cast<const DelayedGroupFilter*>(filt)) {
      delayedgroup_filter_ = i + 1;
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

void
Tally::init_triggers(pugi::xml_node node, int i_tally)
{
  //TODO: use id_ attribute when it's available in C++
  int id_;
  auto err = openmc_tally_get_id(i_tally, &id_);

  for (auto trigger_node: node.children("trigger")) {
    // Read the trigger type.
    TriggerMetric metric;
    if (check_for_node(trigger_node, "type")) {
      auto type_str = get_node_value(trigger_node, "type");
      if (type_str == "std_dev") {
        metric = TriggerMetric::standard_deviation;
      } else if (type_str == "variance") {
        metric = TriggerMetric::variance;
      } else if (type_str == "rel_err") {
        metric = TriggerMetric::relative_error;
      } else {
        std::stringstream msg;
        msg << "Unknown trigger type \"" << type_str << "\" in tally "  << id_;
        fatal_error(msg);
      }
    } else {
      std::stringstream msg;
      msg << "Must specify trigger type for tally " << id_
          << " in tally XML file";
      fatal_error(msg);
    }

    // Read the trigger threshold.
    double threshold;
    if (check_for_node(trigger_node, "threshold")) {
      threshold = std::stod(get_node_value(trigger_node, "threshold"));
    } else {
      std::stringstream msg;
      msg << "Must specify trigger threshold for tally " << id_
          << " in tally XML file";
      fatal_error(msg);
    }

    // Read the trigger scores.
    std::vector<std::string> trigger_scores;
    if (check_for_node(trigger_node, "scores")) {
      trigger_scores = get_node_array<std::string>(trigger_node, "scores");
    } else {
      trigger_scores.push_back("all");
    }

    //TODO: change this when tally scores are moved to C++
    // Get access to the tally's scores.
    int* tally_scores;
    int n_tally_scores;
    auto err = openmc_tally_get_scores(i_tally, &tally_scores, &n_tally_scores);

    // Parse the trigger scores and populate the triggers_ vector.
    for (auto score_str : trigger_scores) {
      if (score_str == "all") {
        triggers_.reserve(triggers_.size() + n_tally_scores);
        for (auto i_score = 0; i_score < n_tally_scores; ++i_score) {
          triggers_.push_back({metric, threshold, i_score});
        }
      } else {
        int i_score = 0;
        for (; i_score < n_tally_scores; ++i_score) {
          if (reaction_name(tally_scores[i_score]) == score_str) break;
        }
        if (i_score == n_tally_scores) {
          std::stringstream msg;
          msg << "Could not find the score \"" << score_str << "\" in tally "
              << id_ << " but it was listed in a trigger on that tally";
          fatal_error(msg);
        }
        triggers_.push_back({metric, threshold, i_score});
      }
    }
  }
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

//! Score tallies based on a simple count of events (for continuous energy).
//
//! Analog tallies ar etriggered at every collision, not every event.

extern "C" void
score_analog_tally_ce(Particle* p)
{
  for (auto i_tally : model::active_analog_tallies) {
    //TODO: off-by-one
    const Tally& tally {*model::tallies[i_tally-1]};
    auto n_score_bins = tally_get_n_score_bins(i_tally);

    //--------------------------------------------------------------------------
    // Loop through all relevant filters and find the filter bins and weights
    // for this event.

    // Find all valid bins in each filter if they have not already been found
    // for a previous tally.
    for (auto i_filt : tally.filters()) {
      //TODO: off-by-one
      auto& match {simulation::filter_matches[i_filt-1]};
      if (!match.bins_present_) {
        match.bins_.clear();
        match.weights_.clear();
        //TODO: off-by-one
        model::tally_filters[i_filt-1]
          ->get_all_bins(p, tally.estimator_, match);
        match.bins_present_ = true;
      }

      // If there are no valid bins for this filter, then there is nothing to
      // score so we can move on to the next tally.
      if (match.bins_.size() == 0) goto next_tally;

      // Set the index of the bin used in the first filter combination
      match.i_bin_ = 1;
    }

    //--------------------------------------------------------------------------
    // Loop over filter bins and nuclide bins

    for (bool filter_loop_done = false; !filter_loop_done; ) {

      //------------------------------------------------------------------------
      // Filter logic

      // Determine scoring index and weight for the current filter combination
      int filter_index = 1;
      double filter_weight = 1.;
      for (auto i = 0; i < tally.filters().size(); ++i) {
        auto i_filt = tally.filters(i);
        //TODO: off-by-one
        auto& match {simulation::filter_matches[i_filt-1]};
        auto i_bin = match.i_bin_;
        //TODO: off-by-one
        filter_index += (match.bins_[i_bin-1] - 1) * tally.strides(i);
        filter_weight *= match.weights_[i_bin-1];
      }

      //------------------------------------------------------------------------
      // Nuclide logic

      if (!tally.all_nuclides_) {
        for (auto i = 0; i < tally.nuclides_.size(); ++i) {
          auto i_nuclide = tally.nuclides_[i];

          // Tally this event in the present nuclide bin if that bin represents
          // the event nuclide or the total material.  Note that the i_nuclide
          // and flux arguments for score_general are not used for analog
          // tallies.
          if (i_nuclide == p->event_nuclide || i_nuclide == -1)
            score_general_ce(p, i_tally, i*n_score_bins, filter_index,
              -1, -1., filter_weight);
        }
      } else {
        // In the case that the user has requested to tally all nuclides, we
        // can take advantage of the fact that we know exactly how nuclide
        // bins correspond to nuclide indices.  First, tally the nuclide.
        auto i = p->event_nuclide;
        score_general_ce(p, i_tally, i*n_score_bins, filter_index,
          -1, -1., filter_weight);

        // Now tally the total material.
        i = tally.nuclides_.size();
        score_general_ce(p, i_tally, i*n_score_bins, filter_index,
          -1, -1., filter_weight);
      }

      //------------------------------------------------------------------------
      // Further filter logic

      // Increment the filter bins, starting with the last filter to find the
      // next valid bin combination
      filter_loop_done = true;
      for (int i = tally.filters().size()-1; i >= 0; --i) {
        auto i_filt = tally.filters(i);
        //TODO: off-by-one
        auto& match {simulation::filter_matches[i_filt-1]};
        if (match.i_bin_ < match.bins_.size()) {
          ++match.i_bin_;
          filter_loop_done = false;
          break;
        } else {
          match.i_bin_ = 1;
        }
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;

next_tally:
    ;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : simulation::filter_matches)
    match.bins_present_ = false;
}

//! Score tallies based on a simple count of events (for multigroup).
//
//! Analog tallies ar etriggered at every collision, not every event.

extern "C" void
score_analog_tally_mg(Particle* p)
{
  for (auto i_tally : model::active_analog_tallies) {
    //TODO: off-by-one
    const Tally& tally {*model::tallies[i_tally-1]};
    auto n_score_bins = tally_get_n_score_bins(i_tally);

    //--------------------------------------------------------------------------
    // Loop through all relevant filters and find the filter bins and weights
    // for this event.

    // Find all valid bins in each filter if they have not already been found
    // for a previous tally.
    for (auto i_filt : tally.filters()) {
      //TODO: off-by-one
      auto& match {simulation::filter_matches[i_filt-1]};
      if (!match.bins_present_) {
        match.bins_.clear();
        match.weights_.clear();
        //TODO: off-by-one
        model::tally_filters[i_filt-1]
          ->get_all_bins(p, tally.estimator_, match);
        match.bins_present_ = true;
      }

      // If there are no valid bins for this filter, then there is nothing to
      // score so we can move on to the next tally.
      if (match.bins_.size() == 0) goto next_tally;

      // Set the index of the bin used in the first filter combination
      match.i_bin_ = 1;
    }

    //--------------------------------------------------------------------------
    // Loop over filter bins and nuclide bins

    for (bool filter_loop_done = false; !filter_loop_done; ) {

      //------------------------------------------------------------------------
      // Filter logic

      // Determine scoring index and weight for the current filter combination
      int filter_index = 1;
      double filter_weight = 1.;
      for (auto i = 0; i < tally.filters().size(); ++i) {
        auto i_filt = tally.filters(i);
        //TODO: off-by-one
        auto& match {simulation::filter_matches[i_filt-1]};
        auto i_bin = match.i_bin_;
        //TODO: off-by-one
        filter_index += (match.bins_[i_bin-1] - 1) * tally.strides(i);
        filter_weight *= match.weights_[i_bin-1];
      }

      //------------------------------------------------------------------------
      // Nuclide logic

      for (auto i = 0; i < tally.nuclides_.size(); ++i) {
        auto i_nuclide = tally.nuclides_[i];

        double atom_density = 0.;
        if (i_nuclide > 0) {
          auto j = material_nuclide_index(p->material, i_nuclide);
          if (j == 0) continue;
          atom_density = material_atom_density(p->material, j);
        }

        score_general_mg(p, i_tally, i*n_score_bins, filter_index,
          i_nuclide, atom_density, filter_weight);
      }

      //------------------------------------------------------------------------
      // Further filter logic

      // Increment the filter bins, starting with the last filter to find the
      // next valid bin combination
      filter_loop_done = true;
      for (int i = tally.filters().size()-1; i >= 0; --i) {
        auto i_filt = tally.filters(i);
        //TODO: off-by-one
        auto& match {simulation::filter_matches[i_filt-1]};
        if (match.i_bin_ < match.bins_.size()) {
          ++match.i_bin_;
          filter_loop_done = false;
          break;
        } else {
          match.i_bin_ = 1;
        }
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;

next_tally:
    ;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : simulation::filter_matches)
    match.bins_present_ = false;
}

//! Score tallies using a tracklength estimate of the flux.
//
//! This is triggered at every event (surface crossing, lattice crossing, or
//! collision) and thus cannot be done for tallies that require post-collision
//! information.

extern "C" void
score_tracklength_tally(Particle* p, double distance)
{
  // Determine the tracklength estimate of the flux
  double flux = p->wgt * distance;

  for (auto i_tally : model::active_tracklength_tallies) {
    //TODO: off-by-one
    const Tally& tally {*model::tallies[i_tally-1]};
    auto n_score_bins = tally_get_n_score_bins(i_tally);

    //--------------------------------------------------------------------------
    // Loop through all relevant filters and find the filter bins and weights
    // for this event.

    // Find all valid bins in each filter if they have not already been found
    // for a previous tally.
    for (auto i_filt : tally.filters()) {
      //TODO: off-by-one
      auto& match {simulation::filter_matches[i_filt-1]};
      if (!match.bins_present_) {
        match.bins_.clear();
        match.weights_.clear();
        //TODO: off-by-one
        model::tally_filters[i_filt-1]
          ->get_all_bins(p, tally.estimator_, match);
        match.bins_present_ = true;
      }

      // If there are no valid bins for this filter, then there is nothing to
      // score so we can move on to the next tally.
      if (match.bins_.size() == 0) goto next_tally;

      // Set the index of the bin used in the first filter combination
      match.i_bin_ = 1;
    }

    //--------------------------------------------------------------------------
    // Loop over filter bins and nuclide bins

    for (bool filter_loop_done = false; !filter_loop_done; ) {

      //------------------------------------------------------------------------
      // Filter logic

      // Determine scoring index and weight for the current filter combination
      int filter_index = 1;
      double filter_weight = 1.;
      for (auto i = 0; i < tally.filters().size(); ++i) {
        auto i_filt = tally.filters(i);
        //TODO: off-by-one
        auto& match {simulation::filter_matches[i_filt-1]};
        auto i_bin = match.i_bin_;
        //TODO: off-by-one
        filter_index += (match.bins_[i_bin-1] - 1) * tally.strides(i);
        filter_weight *= match.weights_[i_bin-1];
      }

      //------------------------------------------------------------------------
      // Nuclide logic

      if (tally.all_nuclides_) {
        if (p->material != MATERIAL_VOID)
          score_all_nuclides(p, i_tally, flux*filter_weight, filter_index);
      } else {
        for (auto i = 0; i < tally.nuclides_.size(); ++i) {
          auto i_nuclide = tally.nuclides_[i];

          double atom_density = 0.;
          if (i_nuclide > 0) {
            if (p->material != MATERIAL_VOID) {
              auto j = material_nuclide_index(p->material, i_nuclide);
              if (j == 0) continue;
              atom_density = material_atom_density(p->material, j);
            }
          }

          //TODO: consider replacing this "if" with pointers or templates
          if (settings::run_CE) {
            score_general_ce(p, i_tally, i*n_score_bins, filter_index,
              i_nuclide, atom_density, flux*filter_weight);
          } else {
            score_general_mg(p, i_tally, i*n_score_bins, filter_index,
              i_nuclide, atom_density, flux*filter_weight);
          }
        }
      }

      //------------------------------------------------------------------------
      // Further filter logic

      // Increment the filter bins, starting with the last filter to find the
      // next valid bin combination
      filter_loop_done = true;
      for (int i = tally.filters().size()-1; i >= 0; --i) {
        auto i_filt = tally.filters(i);
        //TODO: off-by-one
        auto& match {simulation::filter_matches[i_filt-1]};
        if (match.i_bin_ < match.bins_.size()) {
          ++match.i_bin_;
          filter_loop_done = false;
          break;
        } else {
          match.i_bin_ = 1;
        }
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;

next_tally:
    ;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : simulation::filter_matches)
    match.bins_present_ = false;
}

//! Score tallies using a 1 / Sigma_t estimate of the flux.
//
//! This is triggered after every collision.  It is invalid for tallies that
//! require post-collison information because it can score reactions that didn't
//! actually occur, and we don't a priori know what the outcome will be for
//! reactions that we didn't sample.  It is assumed the material is not void
//! since collisions do not occur in voids.

extern "C" void
score_collision_tally(Particle* p)
{
  // Determine the collision estimate of the flux
  double flux;
  if (!settings::survival_biasing) {
    flux = p->last_wgt / simulation::material_xs.total;
  } else {
    flux = (p->last_wgt + p->absorb_wgt) / simulation::material_xs.total;
  }

  for (auto i_tally : model::active_collision_tallies) {
    //TODO: off-by-one
    const Tally& tally {*model::tallies[i_tally-1]};
    auto n_score_bins = tally_get_n_score_bins(i_tally);

    //--------------------------------------------------------------------------
    // Loop through all relevant filters and find the filter bins and weights
    // for this event.

    // Find all valid bins in each filter if they have not already been found
    // for a previous tally.
    for (auto i_filt : tally.filters()) {
      //TODO: off-by-one
      auto& match {simulation::filter_matches[i_filt-1]};
      if (!match.bins_present_) {
        match.bins_.clear();
        match.weights_.clear();
        //TODO: off-by-one
        model::tally_filters[i_filt-1]
          ->get_all_bins(p, tally.estimator_, match);
        match.bins_present_ = true;
      }

      // If there are no valid bins for this filter, then there is nothing to
      // score so we can move on to the next tally.
      if (match.bins_.size() == 0) goto next_tally;

      // Set the index of the bin used in the first filter combination
      match.i_bin_ = 1;
    }

    //--------------------------------------------------------------------------
    // Loop over filter bins and nuclide bins

    for (bool filter_loop_done = false; !filter_loop_done; ) {

      //------------------------------------------------------------------------
      // Filter logic

      // Determine scoring index and weight for the current filter combination
      int filter_index = 1;
      double filter_weight = 1.;
      for (auto i = 0; i < tally.filters().size(); ++i) {
        auto i_filt = tally.filters(i);
        //TODO: off-by-one
        auto& match {simulation::filter_matches[i_filt-1]};
        auto i_bin = match.i_bin_;
        //TODO: off-by-one
        filter_index += (match.bins_[i_bin-1] - 1) * tally.strides(i);
        filter_weight *= match.weights_[i_bin-1];
      }

      //------------------------------------------------------------------------
      // Nuclide logic

      if (tally.all_nuclides_) {
        score_all_nuclides(p, i_tally, flux*filter_weight, filter_index);
      } else {
        for (auto i = 0; i < tally.nuclides_.size(); ++i) {
          auto i_nuclide = tally.nuclides_[i];

          double atom_density = 0.;
          if (i_nuclide > 0) {
            auto j = material_nuclide_index(p->material, i_nuclide);
            if (j == 0) continue;
            atom_density = material_atom_density(p->material, j);
          }

          //TODO: consider replacing this "if" with pointers or templates
          if (settings::run_CE) {
            score_general_ce(p, i_tally, i*n_score_bins, filter_index,
              i_nuclide, atom_density, flux*filter_weight);
          } else {
            score_general_mg(p, i_tally, i*n_score_bins, filter_index,
              i_nuclide, atom_density, flux*filter_weight);
          }
        }
      }

      //------------------------------------------------------------------------
      // Further filter logic

      // Increment the filter bins, starting with the last filter to find the
      // next valid bin combination
      filter_loop_done = true;
      for (int i = tally.filters().size()-1; i >= 0; --i) {
        auto i_filt = tally.filters(i);
        //TODO: off-by-one
        auto& match {simulation::filter_matches[i_filt-1]};
        if (match.i_bin_ < match.bins_.size()) {
          ++match.i_bin_;
          filter_loop_done = false;
          break;
        } else {
          match.i_bin_ = 1;
        }
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;

next_tally:
    ;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : simulation::filter_matches)
    match.bins_present_ = false;
}

//! Score surface or mesh-surface tallies for particle currents.

static void
score_surface_tally_inner(Particle* p, const std::vector<int>& tallies)
{
  // No collision, so no weight change when survival biasing
  double flux = p->wgt;

  for (auto i_tally : tallies) {
    //TODO: off-by-one
    const Tally& tally {*model::tallies[i_tally-1]};
    auto n_score_bins = tally_get_n_score_bins(i_tally);
    auto results = tally_results(i_tally);

    //--------------------------------------------------------------------------
    // Loop through all relevant filters and find the filter bins and weights
    // for this event.

    // Find all valid bins in each filter if they have not already been found
    // for a previous tally.
    for (auto i_filt : tally.filters()) {
      //TODO: off-by-one
      auto& match {simulation::filter_matches[i_filt-1]};
      if (!match.bins_present_) {
        match.bins_.clear();
        match.weights_.clear();
        //TODO: off-by-one
        model::tally_filters[i_filt-1]
          ->get_all_bins(p, tally.estimator_, match);
        match.bins_present_ = true;
      }

      // If there are no valid bins for this filter, then there is nothing to
      // score so we can move on to the next tally.
      if (match.bins_.size() == 0) goto next_tally;

      // Set the index of the bin used in the first filter combination
      match.i_bin_ = 1;
    }

    //--------------------------------------------------------------------------
    // Loop over filter bins and nuclide bins

    for (bool filter_loop_done = false; !filter_loop_done; ) {

      //------------------------------------------------------------------------
      // Filter logic

      // Determine scoring index and weight for the current filter combination
      int filter_index = 1;
      double filter_weight = 1.;
      for (auto i = 0; i < tally.filters().size(); ++i) {
        auto i_filt = tally.filters(i);
        //TODO: off-by-one
        auto& match {simulation::filter_matches[i_filt-1]};
        auto i_bin = match.i_bin_;
        //TODO: off-by-one
        filter_index += (match.bins_[i_bin-1] - 1) * tally.strides(i);
        filter_weight *= match.weights_[i_bin-1];
      }

      //------------------------------------------------------------------------
      // Scoring loop

      // There is only one score type for current tallies so there is no need
      // for a further scoring function.
      double score = flux * filter_weight;
      for (auto score_index = 1; score_index < n_score_bins+1; ++score_index) {
        //TODO: off-by-one
        #pragma omp atomic
        results(filter_index-1, score_index-1, RESULT_VALUE) += score;
      }

      //------------------------------------------------------------------------
      // Further filter logic

      // Increment the filter bins, starting with the last filter to find the
      // next valid bin combination
      filter_loop_done = true;
      for (int i = tally.filters().size()-1; i >= 0; --i) {
        auto i_filt = tally.filters(i);
        //TODO: off-by-one
        auto& match {simulation::filter_matches[i_filt-1]};
        if (match.i_bin_ < match.bins_.size()) {
          ++match.i_bin_;
          filter_loop_done = false;
          break;
        } else {
          match.i_bin_ = 1;
        }
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;

next_tally:
    ;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : simulation::filter_matches)
    match.bins_present_ = false;
}

//! Score mesh-surface tallies for particle currents.

extern "C" void
score_meshsurface_tally(Particle* p)
{
  score_surface_tally_inner(p, model::active_meshsurf_tallies);
}

//! Score surface tallies for particle currents.

extern "C" void
score_surface_tally(Particle* p)
{
  score_surface_tally_inner(p, model::active_surface_tallies);
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
  model::active_meshsurf_tallies.clear();
  model::active_surface_tallies.clear();

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

      case TALLY_MESH_SURFACE:
        model::active_meshsurf_tallies.push_back(i + 1);
        break;

      case TALLY_SURFACE:
        model::active_surface_tallies.push_back(i + 1);
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

  return 0;
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

  return 0;
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

  return 0;
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

  return 0;
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

  int active_tallies_data(int i)
  {return model::active_tallies[i-1];}

  int active_tallies_size()
  {return model::active_tallies.size();}

  int active_analog_tallies_size()
  {return model::active_analog_tallies.size();}

  int active_tracklength_tallies_size()
  {return model::active_tracklength_tallies.size();}

  int active_collision_tallies_size()
  {return model::active_collision_tallies.size();}

  int active_meshsurf_tallies_size()
  {return model::active_meshsurf_tallies.size();}

  int active_surface_tallies_size()
  {return model::active_surface_tallies.size();}

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

  int tally_get_energyout_filter_c(Tally* tally)
  {return tally->energyout_filter_;}

  int tally_get_delayedgroup_filter_c(Tally* tally)
  {return tally->delayedgroup_filter_;}

  void tally_init_triggers(Tally* tally, int i_tally, pugi::xml_node* node)
  {tally->init_triggers(*node, i_tally);}

  int tally_get_deriv_c(Tally* tally) {return tally->deriv_;}

  int tally_set_deriv_c(Tally* tally, int deriv) {tally->deriv_ = deriv;}
}

} // namespace openmc
