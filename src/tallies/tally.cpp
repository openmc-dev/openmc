#include "openmc/tallies/tally.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/reaction_product.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
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

namespace simulation {
  xt::xtensor_fixed<double, xt::xshape<N_GLOBAL_TALLIES, 3>> global_tallies;
  int32_t n_realizations {0};
}

double global_tally_absorption;
double global_tally_collision;
double global_tally_tracklength;
double global_tally_leakage;

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
Tally::init_from_xml(pugi::xml_node node)
{
  if (check_for_node(node, "name")) name_ = get_node_value(node, "name");
}

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
    if (i_filt < 0 || i_filt >= model::tally_filters.size())
      throw std::out_of_range("Index in tally filter array out of bounds.");

    // Keep track of indices for special filters.
    const auto* filt = model::tally_filters[i_filt].get();
    if (dynamic_cast<const EnergyoutFilter*>(filt)) {
      energyout_filter_ = i;
    } else if (dynamic_cast<const DelayedGroupFilter*>(filt)) {
      delayedgroup_filter_ = i;
    }
  }

  // Set the strides.  Filters are traversed in reverse so that the last filter
  // has the shortest stride in memory and the first filter has the longest
  // stride.
  strides_.resize(n, 0);
  int stride = 1;
  for (int i = n-1; i >= 0; --i) {
    strides_[i] = stride;
    stride *= model::tally_filters[filters_[i]]->n_bins_;
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
  bool legendre_present = false;
  bool cell_present = false;
  bool cellfrom_present = false;
  bool surface_present = false;
  bool meshsurface_present = false;
  for (auto i_filt : filters_) {
    const auto* filt {model::tally_filters[i_filt].get()};
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
Tally::set_nuclides(pugi::xml_node node)
{
  nuclides_.clear();

  // By default, we tally just the total material rates.
  if (!check_for_node(node, "nuclides")) {
    nuclides_.push_back(-1);
    return;
  }

  if (get_node_value(node, "nuclides") == "all") {
    // This tally should bin every nuclide in the problem.  It should also bin
    // the total material rates.  To achieve this, set the nuclides_ vector to
    // 0, 1, 2, ..., -1.
    nuclides_.reserve(data::nuclides.size() + 1);
    for (auto i = 0; i < data::nuclides.size(); ++i)
      nuclides_.push_back(i);
    nuclides_.push_back(-1);
    all_nuclides_ = true;

  } else {
    // The user provided specifics nuclides.  Parse it as an array with either
    // "total" or a nuclide name like "U-235" in each position.
    auto words = get_node_array<std::string>(node, "nuclides");
    for (auto word : words) {
      if (word == "total") {
        nuclides_.push_back(-1);
      } else {
        auto search = data::nuclide_map.find(word);
        if (search == data::nuclide_map.end())
          fatal_error("Could not find the nuclide " + word
            + " specified in tally " + std::to_string(id_)
            + " in any material");
        nuclides_.push_back(search->second);
      }
    }
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

    // Parse the trigger scores and populate the triggers_ vector.
    const auto& tally {*model::tallies[i_tally]};
    for (auto score_str : trigger_scores) {
      if (score_str == "all") {
        triggers_.reserve(triggers_.size() + tally.scores_.size());
        for (auto i_score = 0; i_score < tally.scores_.size(); ++i_score) {
          triggers_.push_back({metric, threshold, i_score});
        }
      } else {
        int i_score = 0;
        for (; i_score < tally.scores_.size(); ++i_score) {
          if (reaction_name(tally.scores_[i_score]) == score_str) break;
        }
        if (i_score == tally.scores_.size()) {
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

void Tally::init_results()
{
  int n_scores = scores_.size() * nuclides_.size();
  results_ = xt::empty<double>({n_filter_bins_, n_scores, 3});
}

void Tally::accumulate()
{
  // Increment number of realizations
  n_realizations_ += settings::reduce_tallies ? 1 : mpi::n_procs;

  if (mpi::master || !settings::reduce_tallies) {
    // Calculate total source strength for normalization
    double total_source = 0.0;
    if (settings::run_mode == RUN_MODE_FIXEDSOURCE) {
      for (const auto& s : model::external_sources) {
        total_source += s.strength();
      }
    } else {
      total_source = 1.0;
    }

    // Accumulate each result
    for (int i = 0; i < results_.shape()[0]; ++i) {
      for (int j = 0; j < results_.shape()[1]; ++j) {
        double val = results_(i, j, RESULT_VALUE) /
          simulation::total_weight * total_source;
        results_(i, j, RESULT_VALUE) = 0.0;
        results_(i, j, RESULT_SUM) += val;
        results_(i, j, RESULT_SUM_SQ) += val*val;
      }
    }
  }
}

//==============================================================================
// Non-member functions
//==============================================================================

#ifdef OPENMC_MPI
void reduce_tally_results()
{
  for (int i_tally : model::active_tallies) {
    // Skip any tallies that are not active
    auto& tally {model::tallies[i_tally]};

    // Get view of accumulated tally values
    auto values_view = xt::view(tally->results_, xt::all(), xt::all(), RESULT_VALUE);

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
  auto& gt = simulation::global_tallies;
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
  MPI_Reduce(&simulation::total_weight, &weight_reduced, 1, MPI_DOUBLE, MPI_SUM,
    0, mpi::intracomm);
  if (mpi::master) simulation::total_weight = weight_reduced;
}
#endif

void
accumulate_tallies()
{
#ifdef OPENMC_MPI
  // Combine tally results onto master process
  if (settings::reduce_tallies) reduce_tally_results();
#endif

  // Increase number of realizations (only used for global tallies)
  simulation::n_realizations += settings::reduce_tallies ? 1 : mpi::n_procs;

  // Accumulate on master only unless run is not reduced then do it on all
  if (mpi::master || !settings::reduce_tallies) {
    auto& gt = simulation::global_tallies;

    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      if (simulation::current_batch > settings::n_inactive) {
        // Accumulate products of different estimators of k
        double k_col = gt(K_COLLISION, RESULT_VALUE) / simulation::total_weight;
        double k_abs = gt(K_ABSORPTION, RESULT_VALUE) / simulation::total_weight;
        double k_tra = gt(K_TRACKLENGTH, RESULT_VALUE) / simulation::total_weight;
        simulation::k_col_abs += k_col * k_abs;
        simulation::k_col_tra += k_col * k_tra;
        simulation::k_abs_tra += k_abs * k_tra;
      }
    }

    // Accumulate results for global tallies
    for (int i = 0; i < N_GLOBAL_TALLIES; ++i) {
      double val = gt(i, RESULT_VALUE)/simulation::total_weight;
      gt(i, RESULT_VALUE) = 0.0;
      gt(i, RESULT_SUM) += val;
      gt(i, RESULT_SUM_SQ) += val*val;
    }
  }

  // Accumulate results for each tally
  for (int i_tally : model::active_tallies) {
    auto& tally {model::tallies[i_tally]};
    tally->accumulate();
  }
}

void
setup_active_tallies()
{
  model::active_tallies.clear();
  model::active_analog_tallies.clear();
  model::active_tracklength_tallies.clear();
  model::active_collision_tallies.clear();
  model::active_meshsurf_tallies.clear();
  model::active_surface_tallies.clear();

  for (auto i = 0; i < model::tallies.size(); ++i) {
    const auto& tally {*model::tallies[i]};

    if (tally.active_) {
      model::active_tallies.push_back(i);
      switch (tally.type_) {

      case TALLY_VOLUME:
        switch (tally.estimator_) {
          case ESTIMATOR_ANALOG:
            model::active_analog_tallies.push_back(i);
            break;
          case ESTIMATOR_TRACKLENGTH:
            model::active_tracklength_tallies.push_back(i);
            break;
          case ESTIMATOR_COLLISION:
            model::active_collision_tallies.push_back(i);
        }
        break;

      case TALLY_MESH_SURFACE:
        model::active_meshsurf_tallies.push_back(i);
        break;

      case TALLY_SURFACE:
        model::active_surface_tallies.push_back(i);
      }

      // Check if tally contains depletion reactions and if so, set flag
      if (tally.depletion_rx_) simulation::need_depletion_rx = true;
    }
  }
}

extern "C" void
free_memory_tally_c()
{
  #pragma omp parallel
  {
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

extern "C" int
openmc_tally_set_nuclides(int32_t index, int n, const char** nuclides)
{
  // Make sure the index fits in the array bounds.
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  std::vector<std::string> words(nuclides, nuclides+n);
  std::vector<int> nucs;
  for (auto word : words){
    if (word == "total") {
      nucs.push_back(-1);
    } else {
      auto search = data::nuclide_map.find(word);
      if (search == data::nuclide_map.end()) {
        set_errmsg("Nuclide \"" + word  + "\" has not been loaded yet");
        return OPENMC_E_DATA;
      }
      nucs.push_back(search->second);
    }
  }

  //TODO: off-by-one
  model::tallies[index-1]->nuclides_ = nucs;

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

//! Reset tally results and number of realizations
extern "C" int
openmc_tally_reset(int32_t index)
{
  // Make sure the index fits in the array bounds.
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  // TODO: off-by-one
  auto& t {model::tallies[index-1]};
  t->n_realizations_ = 0;
  // TODO: Change to zero when xtensor is updated
  if (t->results_.size() != 1) {
    xt::view(t->results_, xt::all()) = 0.0;
  }
  return 0;
}

extern "C" int
openmc_tally_get_n_realizations(int32_t index, int32_t* n)
{
  // Make sure the index fits in the array bounds.
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  // TODO: off-by-one
  *n = model::tallies[index - 1]->n_realizations_;
  return 0;
}

//! \brief Returns a pointer to a tally results array along with its shape. This
//! allows a user to obtain in-memory tally results from Python directly.
extern "C" int
openmc_tally_results(int32_t index, double** results, size_t* shape)
{
  // Make sure the index fits in the array bounds.
  if (index < 1 || index > model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  // TODO: off-by-one
  const auto& t {model::tallies[index - 1]};
  // TODO: Change to zero when xtensor is updated
  if (t->results_.size() == 1) {
    set_errmsg("Tally results have not been allocated yet.");
    return OPENMC_E_ALLOCATE;
  }

  // Set pointer to results and copy shape
  *results = t->results_.data();
  auto s = t->results_.shape();
  shape[0] = s[0];
  shape[1] = s[1];
  shape[2] = s[2];
  return 0;
}

extern "C" int
openmc_global_tallies(double** ptr)
{
  *ptr = simulation::global_tallies.data();
  return 0;
}

extern "C" size_t tallies_size() { return model::tallies.size(); }

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

  void tally_init_from_xml(Tally* tally, pugi::xml_node* node)
  {tally->init_from_xml(*node);}

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

  int32_t tally_get_n_filter_bins_c(Tally* tally)
  {return tally->n_filter_bins();}

  int tally_get_n_nuclide_bins_c(Tally* tally)
  {return tally->nuclides_.size();}

  int tally_get_nuclide_bins_c(Tally* tally, int i)
  {return tally->nuclides_[i-1];}

  int tally_get_energyout_filter_c(Tally* tally)
  {return tally->energyout_filter_;}

  void tally_set_scores(Tally* tally, pugi::xml_node* node)
  {tally->set_scores(*node);}

  void tally_set_nuclides(Tally* tally, pugi::xml_node* node)
  {tally->set_nuclides(*node);}

  void tally_init_triggers(Tally* tally, int i_tally, pugi::xml_node* node)
  {tally->init_triggers(*node, i_tally);}

  int tally_get_deriv_c(Tally* tally) {return tally->deriv_;}

  int tally_set_deriv_c(Tally* tally, int deriv) {tally->deriv_ = deriv;}
}

} // namespace openmc
