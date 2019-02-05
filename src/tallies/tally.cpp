#include "openmc/tallies/tally.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/filter_cellfrom.h"
#include "openmc/tallies/filter_delayedgroup.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_legendre.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/filter_meshsurface.h"
#include "openmc/tallies/filter_surface.h"
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
//! An iterator over all combinations of a tally's matching filter bins.
//
//! This iterator handles two distinct tasks.  First, it maps the N-dimensional
//! space created by the indices of N filters onto a 1D sequence.  In other
//! words, it provides a single number that uniquely identifies a combination of
//! bins for many filters.  Second, it handles the task of finding each all
//! valid combinations of filter bins given that each filter can have 1 or 2 or
//! many bins that are valid for the current tally event.
//==============================================================================

class FilterBinIter
{
public:
  FilterBinIter(const Tally& tally, Particle* p, bool end)
    : tally_{tally}
  {
    // Handle the special case for an iterator that points to the end.
    if (end) {
      index_ = -1;
      return;
    }

    // Find all valid bins in each relevant filter if they have not already been
    // found for this event.
    for (auto i_filt : tally_.filters()) {
      //TODO: off-by-one
      auto& match {simulation::filter_matches[i_filt-1]};
      if (!match.bins_present_) {
        match.bins_.clear();
        match.weights_.clear();
        //TODO: off-by-one
        model::tally_filters[i_filt-1]
          ->get_all_bins(p, tally_.estimator_, match);
        match.bins_present_ = true;
      }

      // If there are no valid bins for this filter, then there are no valid
      // filter bin combinations so all iterators are end iterators.
      if (match.bins_.size() == 0) {
        index_ = -1;
        return;
      }

      // Set the index of the bin used in the first filter combination
      match.i_bin_ = 1;
    }

    // Compute the initial index and weight.
    compute_index_weight();
  }

  bool
  operator==(const FilterBinIter& other)
  {
    return index_ == other.index_;
  }

  bool
  operator!=(const FilterBinIter& other)
  {
    return !(*this == other);
  }

  FilterBinIter&
  operator++()
  {
    // Find the next valid combination of filter bins.  To do this, we search
    // backwards through the filters until we find the first filter whose bins
    // can be incremented.
    bool done_looping = true;
    for (int i = tally_.filters().size()-1; i >= 0; --i) {
      auto i_filt = tally_.filters(i);
      //TODO: off-by-one
      auto& match {simulation::filter_matches[i_filt-1]};
      if (match.i_bin_< match.bins_.size()) {
        // The bin for this filter can be incremented.  Increment it and do not
        // touch any of the remaining filters.
        ++match.i_bin_;
        done_looping = false;
        break;
      } else {
        // This bin cannot be incremented so reset it and continue to the next
        // filter.
        match.i_bin_ = 1;
      }
    }

    if (done_looping) {
      // We have visited every valid combination.  All done!
      index_ = -1;
    } else {
      // The loop found a new valid combination.  Compute the corresponding
      // index and weight.
      compute_index_weight();
    }

    return *this;
  }

  int index_ {1};
  double weight_ {1.};

private:
  void
  compute_index_weight()
  {
    index_ = 1;
    weight_ = 1.;
    for (auto i = 0; i < tally_.filters().size(); ++i) {
      auto i_filt = tally_.filters(i);
      //TODO: off-by-one
      auto& match {simulation::filter_matches[i_filt-1]};
      auto i_bin = match.i_bin_;
      //TODO: off-by-one
      index_ += (match.bins_[i_bin-1] - 1) * tally_.strides(i);
      weight_ *= match.weights_[i_bin-1];
    }
  }

  const Tally& tally_;
};

int
score_str_to_int(std::string score_str)
{
  if (score_str == "flux")
    return SCORE_FLUX;

  if (score_str == "total" || score_str == "(n,total)")
    return SCORE_TOTAL;

  if (score_str == "scatter")
    return SCORE_SCATTER;

  if (score_str == "nu-scatter")
    return SCORE_NU_SCATTER;

  if (score_str == "absorption")
    return SCORE_ABSORPTION;

  if (score_str == "fission" || score_str == "18")
    return SCORE_FISSION;

  if (score_str == "nu-fission")
    return SCORE_NU_FISSION;

  if (score_str == "decay-rate")
    return SCORE_DECAY_RATE;

  if (score_str == "delayed-nu-fission")
    return SCORE_DELAYED_NU_FISSION;

  if (score_str == "prompt-nu-fission")
    return SCORE_PROMPT_NU_FISSION;

  if (score_str == "kappa-fission")
    return SCORE_KAPPA_FISSION;

  if (score_str == "inverse-velocity")
    return SCORE_INVERSE_VELOCITY;

  if (score_str == "fission-q-prompt")
    return SCORE_FISS_Q_PROMPT;

  if (score_str == "fission-q-recoverable")
    return SCORE_FISS_Q_RECOV;

  if (score_str == "current")
    return SCORE_CURRENT;

  if (score_str == "events")
    return SCORE_EVENTS;

  if (score_str == "elastic" || score_str == "(n,elastic)")
    return ELASTIC;

  if (score_str == "n2n" || score_str == "(n,2n)")
    return N_2N;

  if (score_str == "n3n" || score_str == "(n,3n)")
    return N_3N;

  if (score_str == "n4n" || score_str == "(n,4n)")
    return N_4N;

  if (score_str == "(n,2nd)")
    return N_2ND;
  if (score_str == "(n,na)")
    return N_2NA;
  if (score_str == "(n,n3a)")
    return N_N3A;
  if (score_str == "(n,2na)")
    return N_2NA;
  if (score_str == "(n,3na)")
    return N_3NA;
  if (score_str == "(n,np)")
    return N_NP;
  if (score_str == "(n,n2a)")
    return N_N2A;
  if (score_str == "(n,2n2a)")
    return N_2N2A;
  if (score_str == "(n,nd)")
    return N_ND;
  if (score_str == "(n,nt)")
    return N_NT;
  if (score_str == "(n,nHe-3)")
    return N_N3HE;
  if (score_str == "(n,nd2a)")
    return N_ND2A;
  if (score_str == "(n,nt2a)")
    return N_NT2A;
  if (score_str == "(n,3nf)")
    return N_3NF;
  if (score_str == "(n,2np)")
    return N_2NP;
  if (score_str == "(n,3np)")
    return N_3NP;
  if (score_str == "(n,n2p)")
    return N_N2P;
  if (score_str == "(n,npa)")
    return N_NPA;
  if (score_str == "(n,n1)")
    return N_N1;
  if (score_str == "(n,nc)")
    return N_NC;
  if (score_str == "(n,gamma)")
    return N_GAMMA;
  if (score_str == "(n,p)")
    return N_P;
  if (score_str == "(n,d)")
    return N_D;
  if (score_str == "(n,t)")
    return N_T;
  if (score_str == "(n,3He)")
    return N_3HE;
  if (score_str == "(n,a)")
    return N_A;
  if (score_str == "(n,2a)")
    return N_2A;
  if (score_str == "(n,3a)")
    return N_3A;
  if (score_str == "(n,2p)")
    return N_2P;
  if (score_str == "(n,pa)")
    return N_PA;
  if (score_str == "(n,t2a)")
    return N_T2A;
  if (score_str == "(n,d2a)")
    return N_D2A;
  if (score_str == "(n,pd)")
    return N_PD;
  if (score_str == "(n,pt)")
    return N_PT;
  if (score_str == "(n,da)")
    return N_DA;

  // So far we have not identified this score string.  Check to see if it is a
  // deprecated score.
  if (score_str.rfind("scatter-", 0) == 0
      || score_str.rfind("nu-scatter-", 0) == 0
      || score_str.rfind("total-y", 0) == 0
      || score_str.rfind("flux-y", 0) == 0)
    fatal_error(score_str + " is no longer an available score");


  // Assume the given string is a reaction MT number.  Make sure it's a natural
  // number then return.
  int MT;
  try {
    MT = std::stoi(score_str);
  } catch (const std::invalid_argument& ex) {
    throw std::invalid_argument("Invalid tally score \"" + score_str + "\"");
  }
  if (MT < 1)
    throw std::invalid_argument("Invalid tally score \"" + score_str + "\"");
  return MT;
}

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
Tally::set_scores(pugi::xml_node node)
{
  if (!check_for_node(node, "scores"))
    fatal_error("No scores specified on tally " + std::to_string(id_));

  auto scores = get_node_array<std::string>(node, "scores");
  set_scores(scores);
}

void
Tally::set_scores(std::vector<std::string> scores)
{
  // Reset state and prepare for the new scores.
  scores_.clear();
  depletion_rx_ = false;
  scores_.reserve(scores.size());

  // Check for the presence of certain restrictive filters.
  bool energyout_present = energyout_filter_ != C_NONE;
  //bool delayedgroup_present = false;
  bool legendre_present = false;
  bool cell_present = false;
  bool cellfrom_present = false;
  bool surface_present = false;
  bool meshsurface_present = false;
  for (auto i_filt : filters_) {
    //TODO: off-by-one
    const auto* filt {model::tally_filters[i_filt-1].get()};
    if (dynamic_cast<const LegendreFilter*>(filt)) {
      legendre_present = true;
    } else if (dynamic_cast<const CellFromFilter*>(filt)) {
      cellfrom_present = true;
    } else if (dynamic_cast<const CellFilter*>(filt)) {
      cell_present = true;
    } else if (dynamic_cast<const SurfaceFilter*>(filt)) {
      surface_present = true;
    } else if (dynamic_cast<const MeshSurfaceFilter*>(filt)) {
      meshsurface_present = true;
    }
  }

  // Iterate over the given scores.
  for (auto score_str : scores) {
    // Make sure a delayed group filter wasn't used with an incompatible score.
    bool has_delayedgroup = delayedgroup_filter_ != C_NONE;
    if (delayedgroup_filter_ != C_NONE) {
      if (score_str != "delayed-nu-fission" && score_str != "decay-rate")
        fatal_error("Cannot tally " + score_str + "with a delayedgroup filter");
    }

    auto score = score_str_to_int(score_str);

    switch (score) {
    case SCORE_FLUX:
      if (!nuclides_.empty())
        if (!(nuclides_.size() == 1 && nuclides_[0] == -1))
          fatal_error("Cannot tally flux for an individual nuclide.");
      if (energyout_present)
        fatal_error("Cannot tally flux with an outgoing energy filter.");
      break;

    case SCORE_TOTAL:
    case SCORE_ABSORPTION:
    case SCORE_FISSION:
      if (energyout_present)
        fatal_error("Cannot tally " + score_str + " reaction rate with an "
          "outgoing energy filter");
      break;

    case SCORE_SCATTER:
      if (legendre_present)
        estimator_ = ESTIMATOR_ANALOG;
    case SCORE_NU_FISSION:
    case SCORE_DELAYED_NU_FISSION:
    case SCORE_PROMPT_NU_FISSION:
      if (energyout_present)
        estimator_ = ESTIMATOR_ANALOG;
      break;

    case SCORE_NU_SCATTER:
      if (settings::run_CE) {
        estimator_ = ESTIMATOR_ANALOG;
      } else {
        if (energyout_present || legendre_present)
          estimator_ = ESTIMATOR_ANALOG;
      }
      break;

    case N_2N:
    case N_3N:
    case N_4N:
    case N_GAMMA:
    case N_P:
    case N_A:
      depletion_rx_ = true;
      break;

    case SCORE_CURRENT:
      // Check which type of current is desired: mesh or surface currents.
      if (surface_present || cell_present || cellfrom_present) {
        if (meshsurface_present)
          fatal_error("Cannot tally mesh surface currents in the same tally as "
            "normal surface currents");
        type_ = TALLY_SURFACE;
      } else if (meshsurface_present) {
        type_ = TALLY_MESH_SURFACE;
      } else {
        fatal_error("Cannot tally currents without surface type filters");
      }
      break;
    }

    scores_.push_back(score);
  }

  // Make sure that no duplicate scores exist.
  for (auto it1 = scores_.begin(); it1 != scores_.end(); ++it1) {
    for (auto it2 = it1 + 1; it2 != scores_.end(); ++it2) {
      if (*it1 == *it2)
        fatal_error("Duplicate score of type \"" + reaction_name(*it1)
          + "\" found in tally " + std::to_string(id_));
    }
  }

  // Make sure all scores are compatible with multigroup mode.
  if (!settings::run_CE) {
    for (auto sc : scores_)
      if (sc > 0)
        fatal_error("Cannot tally " + reaction_name(sc) + " reaction rate "
          "in multi-group mode");
  }

  // Make sure current scores are not mixed in with volumetric scores.
  if (type_ == TALLY_SURFACE || type_ == TALLY_MESH_SURFACE) {
    if (scores_.size() != 1)
      fatal_error("Cannot tally other scores in the same tally as surface "
        "currents");
  }
}

void
Tally::init_triggers(pugi::xml_node node, int i_tally)
{
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

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p, false);
    auto end = FilterBinIter(tally, nullptr, true);
    if (filter_iter == end) continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
      if (!tally.all_nuclides_) {
        for (auto i = 0; i < tally.nuclides_.size(); ++i) {
          auto i_nuclide = tally.nuclides_[i];

          // Tally this event in the present nuclide bin if that bin represents
          // the event nuclide or the total material.  Note that the i_nuclide
          // and flux arguments for score_general are not used for analog
          // tallies.
          if (i_nuclide == p->event_nuclide || i_nuclide == -1)
            score_general_ce(p, i_tally, i*tally.scores_.size(), filter_index,
              -1, -1., filter_weight);
        }

      } else {
        // In the case that the user has requested to tally all nuclides, we
        // can take advantage of the fact that we know exactly how nuclide
        // bins correspond to nuclide indices.  First, tally the nuclide.
        auto i = p->event_nuclide;
        score_general_ce(p, i_tally, i*tally.scores_.size(), filter_index,
          -1, -1., filter_weight);

        // Now tally the total material.
        i = tally.nuclides_.size();
        score_general_ce(p, i_tally, i*tally.scores_.size(), filter_index,
          -1, -1., filter_weight);
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;
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

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p, false);
    auto end = FilterBinIter(tally, nullptr, true);
    if (filter_iter == end) continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
      for (auto i = 0; i < tally.nuclides_.size(); ++i) {
        auto i_nuclide = tally.nuclides_[i];

        double atom_density = 0.;
        if (i_nuclide > 0) {
          auto j = material_nuclide_index(p->material, i_nuclide);
          if (j == 0) continue;
          atom_density = material_atom_density(p->material, j);
        }

        score_general_mg(p, i_tally, i*tally.scores_.size(), filter_index,
          i_nuclide, atom_density, filter_weight);
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;
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

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p, false);
    auto end = FilterBinIter(tally, nullptr, true);
    if (filter_iter == end) continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
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
            score_general_ce(p, i_tally, i*tally.scores_.size(), filter_index,
              i_nuclide, atom_density, flux*filter_weight);
          } else {
            score_general_mg(p, i_tally, i*tally.scores_.size(), filter_index,
              i_nuclide, atom_density, flux*filter_weight);
          }
        }
      }

    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;
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

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p, false);
    auto end = FilterBinIter(tally, nullptr, true);
    if (filter_iter == end) continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
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
            score_general_ce(p, i_tally, i*tally.scores_.size(), filter_index,
              i_nuclide, atom_density, flux*filter_weight);
          } else {
            score_general_mg(p, i_tally, i*tally.scores_.size(), filter_index,
              i_nuclide, atom_density, flux*filter_weight);
          }
        }
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;
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
    auto results = tally_results(i_tally);

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p, false);
    auto end = FilterBinIter(tally, nullptr, true);
    if (filter_iter == end) continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over scores.
      // There is only one score type for current tallies so there is no need
      // for a further scoring function.
      double score = flux * filter_weight;
      for (auto score_index = 1; score_index < tally.scores_.size()+1;
           ++score_index) {
        //TODO: off-by-one
        #pragma omp atomic
        results(filter_index-1, score_index-1, RESULT_VALUE) += score;
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;
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
  model::active_meshsurf_tallies.clear();
  model::active_surface_tallies.clear();
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
openmc_tally_get_scores(int32_t index, int** scores, int* n)
{
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  //TODO: off-by-one
  *scores = model::tallies[index-1]->scores_.data();
  *n = model::tallies[index-1]->scores_.size();
  return 0;
}

extern "C" int
openmc_tally_set_scores(int32_t index, int n, const char** scores)
{
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  std::vector<std::string> scores_str(scores, scores+n);
  try {
    //TODO: off-by-one
    model::tallies[index-1]->set_scores(scores_str);
  } catch (const std::invalid_argument& ex) {
    set_errmsg(ex.what());
    return OPENMC_E_INVALID_ARGUMENT;
  }

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

  int tally_get_id_c(Tally* tally) {return tally->id_;}

  void tally_set_id_c(Tally* tally, int id) {tally->id_ = id;}

  int tally_get_type_c(Tally* tally) {return tally->type_;}

  void tally_set_type_c(Tally* tally, int type) {tally->type_ = type;}

  int tally_get_estimator_c(Tally* tally) {return tally->estimator_;}

  void tally_set_estimator_c(Tally* tally, int e) {tally->estimator_ = e;}

  bool tally_get_depletion_rx_c(Tally* tally) {return tally->depletion_rx_;}

  int tally_get_n_scores_c(Tally* tally) {return tally->scores_.size();}

  int tally_get_score_c(Tally* tally, int i) {return tally->scores_[i];}

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

  void tally_init_scores(Tally* tally, pugi::xml_node* node)
  {tally->set_scores(*node);}

  void tally_init_triggers(Tally* tally, int i_tally, pugi::xml_node* node)
  {tally->init_triggers(*node, i_tally);}

  int tally_get_deriv_c(Tally* tally) {return tally->deriv_;}

  int tally_set_deriv_c(Tally* tally, int deriv) {tally->deriv_ = deriv;}
}

} // namespace openmc
