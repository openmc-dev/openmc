#include "openmc/tallies/tally.h"

#include "openmc/array.h"
#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/reaction.h"
#include "openmc/reaction_product.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/filter_cellborn.h"
#include "openmc/tallies/filter_cellfrom.h"
#include "openmc/tallies/filter_collision.h"
#include "openmc/tallies/filter_delayedgroup.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_legendre.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/filter_meshborn.h"
#include "openmc/tallies/filter_meshsurface.h"
#include "openmc/tallies/filter_particle.h"
#include "openmc/tallies/filter_sph_harm.h"
#include "openmc/tallies/filter_surface.h"
#include "openmc/xml_interface.h"

#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp" // for empty_like
#include "xtensor/xview.hpp"
#include <fmt/core.h>

#include <algorithm> // for max
#include <cassert>
#include <cstddef> // for size_t
#include <string>

namespace openmc {

//==============================================================================
// Global variable definitions
//==============================================================================

namespace model {
//! a mapping of tally ID to index in the tallies vector
std::unordered_map<int, int> tally_map;
vector<unique_ptr<Tally>> tallies;
vector<int> active_tallies;
vector<int> active_analog_tallies;
vector<int> active_tracklength_tallies;
vector<int> active_collision_tallies;
vector<int> active_meshsurf_tallies;
vector<int> active_surface_tallies;
vector<int> active_pulse_height_tallies;
vector<int> pulse_height_cells;
vector<double> time_grid;
} // namespace model

namespace simulation {
xt::xtensor_fixed<double, xt::xshape<N_GLOBAL_TALLIES, 3>> global_tallies;
int32_t n_realizations {0};
} // namespace simulation

double global_tally_absorption;
double global_tally_collision;
double global_tally_tracklength;
double global_tally_leakage;

//==============================================================================
// Tally object implementation
//==============================================================================

Tally::Tally(int32_t id)
{
  index_ = model::tallies.size(); // Avoids warning about narrowing
  this->set_id(id);
  this->set_filters({});
}

Tally::Tally(pugi::xml_node node)
{
  index_ = model::tallies.size(); // Avoids warning about narrowing

  // Copy and set tally id
  if (!check_for_node(node, "id")) {
    throw std::runtime_error {"Must specify id for tally in tally XML file."};
  }
  int32_t id = std::stoi(get_node_value(node, "id"));
  this->set_id(id);

  if (check_for_node(node, "name"))
    name_ = get_node_value(node, "name");

  if (check_for_node(node, "multiply_density")) {
    multiply_density_ = get_node_value_bool(node, "multiply_density");
  }

  // =======================================================================
  // READ DATA FOR FILTERS

  // Check if user is using old XML format and throw an error if so
  if (check_for_node(node, "filter")) {
    throw std::runtime_error {
      "Tally filters must be specified independently of "
      "tallies in a <filter> element. The <tally> element itself should "
      "have a list of filters that apply, e.g., <filters>1 2</filters> "
      "where 1 and 2 are the IDs of filters specified outside of "
      "<tally>."};
  }

  // Determine number of filters
  vector<int> filter_ids;
  if (check_for_node(node, "filters")) {
    filter_ids = get_node_array<int>(node, "filters");
  }

  // Allocate and store filter user ids
  vector<Filter*> filters;
  for (int filter_id : filter_ids) {
    // Determine if filter ID is valid
    auto it = model::filter_map.find(filter_id);
    if (it == model::filter_map.end()) {
      throw std::runtime_error {fmt::format(
        "Could not find filter {} specified on tally {}", filter_id, id_)};
    }

    // Store the index of the filter
    filters.push_back(model::tally_filters[it->second].get());
  }

  // Set the filters
  this->set_filters(filters);

  // Check for the presence of certain filter types
  bool has_energyout = energyout_filter_ >= 0;
  int particle_filter_index = C_NONE;
  for (int64_t j = 0; j < filters_.size(); ++j) {
    int i_filter = filters_[j];
    const auto& f = model::tally_filters[i_filter].get();

    auto pf = dynamic_cast<ParticleFilter*>(f);
    if (pf)
      particle_filter_index = i_filter;

    // Change the tally estimator if a filter demands it
    FilterType filt_type = f->type();
    if (filt_type == FilterType::ENERGY_OUT ||
        filt_type == FilterType::LEGENDRE) {
      estimator_ = TallyEstimator::ANALOG;
    } else if (filt_type == FilterType::SPHERICAL_HARMONICS) {
      auto sf = dynamic_cast<SphericalHarmonicsFilter*>(f);
      if (sf->cosine() == SphericalHarmonicsCosine::scatter) {
        estimator_ = TallyEstimator::ANALOG;
      }
    } else if (filt_type == FilterType::SPATIAL_LEGENDRE ||
               filt_type == FilterType::ZERNIKE ||
               filt_type == FilterType::ZERNIKE_RADIAL) {
      estimator_ = TallyEstimator::COLLISION;
    }
  }

  // =======================================================================
  // READ DATA FOR NUCLIDES

  this->set_nuclides(node);

  // =======================================================================
  // READ DATA FOR SCORES

  this->set_scores(node);

  if (!check_for_node(node, "scores")) {
    fatal_error(fmt::format("No scores specified on tally {}.", id_));
  }

  // Set IFP if needed
  if (!settings::ifp_on) {
    // Determine if this tally has an IFP score
    bool has_ifp_score = false;
    for (int score : scores_) {
      if (score == SCORE_IFP_TIME_NUM || score == SCORE_IFP_BETA_NUM ||
          score == SCORE_IFP_DENOM) {
        has_ifp_score = true;
        break;
      }
    }

    // Check for errors
    if (has_ifp_score) {
      if (settings::run_mode == RunMode::EIGENVALUE) {
        if (settings::ifp_n_generation < 0) {
          settings::ifp_n_generation = DEFAULT_IFP_N_GENERATION;
          warning(fmt::format(
            "{} generations will be used for IFP (default value). It can be "
            "changed using the 'ifp_n_generation' settings.",
            settings::ifp_n_generation));
        }
        if (settings::ifp_n_generation > settings::n_inactive) {
          fatal_error("'ifp_n_generation' must be lower than or equal to the "
                      "number of inactive cycles.");
        }
        settings::ifp_on = true;
      } else {
        fatal_error(
          "Iterated Fission Probability can only be used in an eigenvalue "
          "calculation.");
      }
    }
  }

  // Set IFP parameters if needed
  if (settings::ifp_on) {
    for (int score : scores_) {
      switch (score) {
      case SCORE_IFP_TIME_NUM:
        if (settings::ifp_parameter == IFPParameter::None) {
          settings::ifp_parameter = IFPParameter::GenerationTime;
        } else if (settings::ifp_parameter == IFPParameter::BetaEffective) {
          settings::ifp_parameter = IFPParameter::Both;
        }
        break;
      case SCORE_IFP_BETA_NUM:
      case SCORE_IFP_DENOM:
        if (settings::ifp_parameter == IFPParameter::None) {
          settings::ifp_parameter = IFPParameter::BetaEffective;
        } else if (settings::ifp_parameter == IFPParameter::GenerationTime) {
          settings::ifp_parameter = IFPParameter::Both;
        }
        break;
      }
    }
  }

  // Check if tally is compatible with particle type
  if (!settings::photon_transport) {
    for (int score : scores_) {
      switch (score) {
      case SCORE_PULSE_HEIGHT:
        fatal_error(
          "For pulse-height tallies, photon transport needs to be activated.");
        break;
      }
    }
  }
  if (settings::photon_transport) {
    if (particle_filter_index == C_NONE) {
      for (int score : scores_) {
        switch (score) {
        case SCORE_INVERSE_VELOCITY:
          fatal_error("Particle filter must be used with photon "
                      "transport on and inverse velocity score");
          break;
        case SCORE_FLUX:
        case SCORE_TOTAL:
        case SCORE_SCATTER:
        case SCORE_NU_SCATTER:
        case SCORE_ABSORPTION:
        case SCORE_FISSION:
        case SCORE_NU_FISSION:
        case SCORE_CURRENT:
        case SCORE_EVENTS:
        case SCORE_DELAYED_NU_FISSION:
        case SCORE_PROMPT_NU_FISSION:
        case SCORE_DECAY_RATE:
          warning("You are tallying the '" + reaction_name(score) +
                  "' score and haven't used a particle filter. This score will "
                  "include contributions from all particles.");
          break;
        }
      }
    }
  } else {
    if (particle_filter_index >= 0) {
      const auto& f = model::tally_filters[particle_filter_index].get();
      auto pf = dynamic_cast<ParticleFilter*>(f);
      for (auto p : pf->particles()) {
        if (p != ParticleType::neutron) {
          warning(fmt::format(
            "Particle filter other than NEUTRON used with "
            "photon transport turned off. All tallies for particle type {}"
            " will have no scores",
            static_cast<int>(p)));
        }
      }
    }
  }

  // Check for a tally derivative.
  if (check_for_node(node, "derivative")) {
    int deriv_id = std::stoi(get_node_value(node, "derivative"));

    // Find the derivative with the given id, and store it's index.
    auto it = model::tally_deriv_map.find(deriv_id);
    if (it == model::tally_deriv_map.end()) {
      fatal_error(fmt::format(
        "Could not find derivative {} specified on tally {}", deriv_id, id_));
    }

    deriv_ = it->second;

    // Only analog or collision estimators are supported for differential
    // tallies.
    if (estimator_ == TallyEstimator::TRACKLENGTH) {
      estimator_ = TallyEstimator::COLLISION;
    }

    const auto& deriv = model::tally_derivs[deriv_];
    if (deriv.variable == DerivativeVariable::NUCLIDE_DENSITY ||
        deriv.variable == DerivativeVariable::TEMPERATURE) {
      for (int i_nuc : nuclides_) {
        if (has_energyout && i_nuc == -1) {
          fatal_error(fmt::format(
            "Error on tally {}: Cannot use a "
            "'nuclide_density' or 'temperature' derivative on a tally with an "
            "outgoing energy filter and 'total' nuclide rate. Instead, tally "
            "each nuclide in the material individually.",
            id_));
          // Note that diff tallies with these characteristics would work
          // correctly if no tally events occur in the perturbed material
          // (e.g. pertrubing moderator but only tallying fuel), but this
          // case would be hard to check for by only reading inputs.
        }
      }
    }
  }

  // If settings.xml trigger is turned on, create tally triggers
  if (settings::trigger_on) {
    this->init_triggers(node);
  }

  // =======================================================================
  // SET TALLY ESTIMATOR

  // Check if user specified estimator
  if (check_for_node(node, "estimator")) {
    std::string est = get_node_value(node, "estimator");
    if (est == "analog") {
      estimator_ = TallyEstimator::ANALOG;
    } else if (est == "tracklength" || est == "track-length" ||
               est == "pathlength" || est == "path-length") {
      // If the estimator was set to an analog estimator, this means the
      // tally needs post-collision information
      if (estimator_ == TallyEstimator::ANALOG ||
          estimator_ == TallyEstimator::COLLISION) {
        throw std::runtime_error {fmt::format("Cannot use track-length "
                                              "estimator for tally {}",
          id_)};
      }

      // Set estimator to track-length estimator
      estimator_ = TallyEstimator::TRACKLENGTH;

    } else if (est == "collision") {
      // If the estimator was set to an analog estimator, this means the
      // tally needs post-collision information
      if (estimator_ == TallyEstimator::ANALOG) {
        throw std::runtime_error {fmt::format("Cannot use collision estimator "
                                              "for tally ",
          id_)};
      }

      // Set estimator to collision estimator
      estimator_ = TallyEstimator::COLLISION;

    } else {
      throw std::runtime_error {
        fmt::format("Invalid estimator '{}' on tally {}", est, id_)};
    }
  }

#ifdef OPENMC_LIBMESH_ENABLED
  // ensure a tracklength tally isn't used with a libMesh filter
  for (auto i : this->filters_) {
    auto df = dynamic_cast<MeshFilter*>(model::tally_filters[i].get());
    if (df) {
      auto lm = dynamic_cast<LibMesh*>(model::meshes[df->mesh()].get());
      if (lm && estimator_ == TallyEstimator::TRACKLENGTH) {
        fatal_error("A tracklength estimator cannot be used with "
                    "an unstructured LibMesh tally.");
      }
    }
  }
#endif
}

Tally::~Tally()
{
  model::tally_map.erase(id_);
}

Tally* Tally::create(int32_t id)
{
  model::tallies.push_back(make_unique<Tally>(id));
  return model::tallies.back().get();
}

void Tally::set_id(int32_t id)
{
  assert(id >= 0 || id == C_NONE);

  // Clear entry in tally map if an ID was already assigned before
  if (id_ != C_NONE) {
    model::tally_map.erase(id_);
    id_ = C_NONE;
  }

  // Make sure no other tally has the same ID
  if (model::tally_map.find(id) != model::tally_map.end()) {
    throw std::runtime_error {
      fmt::format("Two tallies have the same ID: {}", id)};
  }

  // If no ID specified, auto-assign next ID in sequence
  if (id == C_NONE) {
    id = 0;
    for (const auto& t : model::tallies) {
      id = std::max(id, t->id_);
    }
    ++id;
  }

  // Update ID and entry in tally map
  id_ = id;
  model::tally_map[id] = index_;
}

std::vector<FilterType> Tally::filter_types() const
{
  std::vector<FilterType> filter_types;
  for (auto idx : this->filters())
    filter_types.push_back(model::tally_filters[idx]->type());
  return filter_types;
}

std::unordered_map<FilterType, int32_t> Tally::filter_indices() const
{
  std::unordered_map<FilterType, int32_t> filter_indices;
  for (int i = 0; i < this->filters().size(); i++) {
    const auto& f = model::tally_filters[this->filters(i)];

    filter_indices[f->type()] = i;
  }
  return filter_indices;
}

bool Tally::has_filter(FilterType filter_type) const
{
  for (auto idx : this->filters()) {
    if (model::tally_filters[idx]->type() == filter_type)
      return true;
  }
  return false;
}

void Tally::set_filters(span<Filter*> filters)
{
  // Clear old data.
  filters_.clear();
  strides_.clear();

  // Copy in the given filter indices.
  auto n = filters.size();
  filters_.reserve(n);

  for (auto* filter : filters) {
    add_filter(filter);
  }
}

void Tally::add_filter(Filter* filter)
{
  int32_t filter_idx = model::filter_map.at(filter->id());
  // if this filter is already present, do nothing and return
  if (std::find(filters_.begin(), filters_.end(), filter_idx) != filters_.end())
    return;

  // Keep track of indices for special filters
  if (filter->type() == FilterType::ENERGY_OUT) {
    energyout_filter_ = filters_.size();
  } else if (filter->type() == FilterType::DELAYED_GROUP) {
    delayedgroup_filter_ = filters_.size();
  } else if (filter->type() == FilterType::CELL) {
    cell_filter_ = filters_.size();
  } else if (filter->type() == FilterType::ENERGY) {
    energy_filter_ = filters_.size();
  }
  filters_.push_back(filter_idx);
}

void Tally::set_strides()
{
  // Set the strides.  Filters are traversed in reverse so that the last filter
  // has the shortest stride in memory and the first filter has the longest
  // stride.
  auto n = filters_.size();
  strides_.resize(n, 0);
  int stride = 1;
  for (int i = n - 1; i >= 0; --i) {
    strides_[i] = stride;
    stride *= model::tally_filters[filters_[i]]->n_bins();
  }
  n_filter_bins_ = stride;
}

void Tally::set_scores(pugi::xml_node node)
{
  if (!check_for_node(node, "scores"))
    fatal_error(fmt::format("No scores specified on tally {}", id_));

  auto scores = get_node_array<std::string>(node, "scores");
  set_scores(scores);
}

void Tally::set_scores(const vector<std::string>& scores)
{
  // Reset state and prepare for the new scores.
  scores_.clear();
  scores_.reserve(scores.size());

  // Check for the presence of certain restrictive filters.
  bool energyout_present = energyout_filter_ != C_NONE;
  bool legendre_present = false;
  bool cell_present = false;
  bool cellfrom_present = false;
  bool surface_present = false;
  bool meshsurface_present = false;
  bool non_cell_energy_present = false;
  for (auto i_filt : filters_) {
    const auto* filt {model::tally_filters[i_filt].get()};
    // Checking for only cell and energy filters for pulse-height tally
    if (!(filt->type() == FilterType::CELL ||
          filt->type() == FilterType::ENERGY)) {
      non_cell_energy_present = true;
    }
    if (filt->type() == FilterType::LEGENDRE) {
      legendre_present = true;
    } else if (filt->type() == FilterType::CELLFROM) {
      cellfrom_present = true;
    } else if (filt->type() == FilterType::CELL) {
      cell_present = true;
    } else if (filt->type() == FilterType::SURFACE) {
      surface_present = true;
    } else if (filt->type() == FilterType::MESH_SURFACE) {
      meshsurface_present = true;
    }
  }

  // Iterate over the given scores.
  for (auto score_str : scores) {
    // Make sure a delayed group filter wasn't used with an incompatible score.
    if (delayedgroup_filter_ != C_NONE) {
      if (score_str != "delayed-nu-fission" && score_str != "decay-rate")
        fatal_error("Cannot tally " + score_str + "with a delayedgroup filter");
    }

    // Determine integer code for score
    int score = reaction_type(score_str);

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
        fatal_error("Cannot tally " + score_str +
                    " reaction rate with an "
                    "outgoing energy filter");
      break;

    case SCORE_SCATTER:
      if (legendre_present)
        estimator_ = TallyEstimator::ANALOG;
    case SCORE_NU_FISSION:
    case SCORE_DELAYED_NU_FISSION:
    case SCORE_PROMPT_NU_FISSION:
      if (energyout_present)
        estimator_ = TallyEstimator::ANALOG;
      break;

    case SCORE_NU_SCATTER:
      if (settings::run_CE) {
        estimator_ = TallyEstimator::ANALOG;
      } else {
        if (energyout_present || legendre_present)
          estimator_ = TallyEstimator::ANALOG;
      }
      break;

    case SCORE_CURRENT:
      // Check which type of current is desired: mesh or surface currents.
      if (surface_present || cell_present || cellfrom_present) {
        if (meshsurface_present)
          fatal_error("Cannot tally mesh surface currents in the same tally as "
                      "normal surface currents");
        type_ = TallyType::SURFACE;
        estimator_ = TallyEstimator::ANALOG;
      } else if (meshsurface_present) {
        type_ = TallyType::MESH_SURFACE;
      } else {
        fatal_error("Cannot tally currents without surface type filters");
      }
      break;

    case HEATING:
      if (settings::photon_transport)
        estimator_ = TallyEstimator::COLLISION;
      break;

    case SCORE_PULSE_HEIGHT:
      if (non_cell_energy_present) {
        fatal_error("Pulse-height tallies are not compatible with filters "
                    "other than CellFilter and EnergyFilter");
      }
      type_ = TallyType::PULSE_HEIGHT;

      // Collecting indices of all cells covered by the filters in the pulse
      // height tally in global variable pulse_height_cells
      for (const auto& i_filt : filters_) {
        auto cell_filter =
          dynamic_cast<CellFilter*>(model::tally_filters[i_filt].get());
        if (cell_filter) {
          const auto& cells = cell_filter->cells();
          for (int i = 0; i < cell_filter->n_bins(); i++) {
            int cell_index = cells[i];
            if (!contains(model::pulse_height_cells, cell_index)) {
              model::pulse_height_cells.push_back(cell_index);
            }
          }
        }
      }
      break;

    case SCORE_IFP_TIME_NUM:
    case SCORE_IFP_BETA_NUM:
    case SCORE_IFP_DENOM:
      estimator_ = TallyEstimator::COLLISION;
      break;
    }

    scores_.push_back(score);
  }

  // Make sure that no duplicate scores exist.
  for (auto it1 = scores_.begin(); it1 != scores_.end(); ++it1) {
    for (auto it2 = it1 + 1; it2 != scores_.end(); ++it2) {
      if (*it1 == *it2)
        fatal_error(
          fmt::format("Duplicate score of type \"{}\" found in tally {}",
            reaction_name(*it1), id_));
    }
  }

  // Make sure all scores are compatible with multigroup mode.
  if (!settings::run_CE) {
    for (auto sc : scores_)
      if (sc > 0)
        fatal_error("Cannot tally " + reaction_name(sc) +
                    " reaction rate "
                    "in multi-group mode");
  }

  // Make sure current scores are not mixed in with volumetric scores.
  if (type_ == TallyType::SURFACE || type_ == TallyType::MESH_SURFACE) {
    if (scores_.size() != 1)
      fatal_error("Cannot tally other scores in the same tally as surface "
                  "currents.");
  }
  if ((surface_present || meshsurface_present) && scores_[0] != SCORE_CURRENT)
    fatal_error("Cannot tally score other than 'current' when using a surface "
                "or mesh-surface filter.");
}

void Tally::set_nuclides(pugi::xml_node node)
{
  nuclides_.clear();

  // By default, we tally just the total material rates.
  if (!check_for_node(node, "nuclides")) {
    nuclides_.push_back(-1);
    return;
  }

  // The user provided specifics nuclides.  Parse it as an array with either
  // "total" or a nuclide name like "U235" in each position.
  auto words = get_node_array<std::string>(node, "nuclides");
  this->set_nuclides(words);
}

void Tally::set_nuclides(const vector<std::string>& nuclides)
{
  nuclides_.clear();

  for (const auto& nuc : nuclides) {
    if (nuc == "total") {
      nuclides_.push_back(-1);
    } else {
      auto search = data::nuclide_map.find(nuc);
      if (search == data::nuclide_map.end()) {
        int err = openmc_load_nuclide(nuc.c_str(), nullptr, 0);
        if (err < 0)
          throw std::runtime_error {openmc_err_msg};
      }
      nuclides_.push_back(data::nuclide_map.at(nuc));
    }
  }
}

void Tally::init_triggers(pugi::xml_node node)
{
  for (auto trigger_node : node.children("trigger")) {
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
        fatal_error(fmt::format(
          "Unknown trigger type \"{}\" in tally {}", type_str, id_));
      }
    } else {
      fatal_error(fmt::format(
        "Must specify trigger type for tally {} in tally XML file", id_));
    }

    // Read the trigger threshold.
    double threshold;
    if (check_for_node(trigger_node, "threshold")) {
      threshold = std::stod(get_node_value(trigger_node, "threshold"));
      if (threshold <= 0) {
        fatal_error("Tally trigger threshold must be positive");
      }
    } else {
      fatal_error(fmt::format(
        "Must specify trigger threshold for tally {} in tally XML file", id_));
    }

    // Read whether to allow zero-tally bins to be ignored.
    bool ignore_zeros = false;
    if (check_for_node(trigger_node, "ignore_zeros")) {
      ignore_zeros = get_node_value_bool(trigger_node, "ignore_zeros");
    }

    // Read the trigger scores.
    vector<std::string> trigger_scores;
    if (check_for_node(trigger_node, "scores")) {
      trigger_scores = get_node_array<std::string>(trigger_node, "scores");
    } else {
      trigger_scores.push_back("all");
    }

    // Parse the trigger scores and populate the triggers_ vector.
    for (auto score_str : trigger_scores) {
      if (score_str == "all") {
        triggers_.reserve(triggers_.size() + this->scores_.size());
        for (auto i_score = 0; i_score < this->scores_.size(); ++i_score) {
          triggers_.push_back({metric, threshold, ignore_zeros, i_score});
        }
      } else {
        int i_score = 0;
        for (; i_score < this->scores_.size(); ++i_score) {
          if (this->scores_[i_score] == reaction_type(score_str))
            break;
        }
        if (i_score == this->scores_.size()) {
          fatal_error(
            fmt::format("Could not find the score \"{}\" in tally "
                        "{} but it was listed in a trigger on that tally",
              score_str, id_));
        }
        triggers_.push_back({metric, threshold, ignore_zeros, i_score});
      }
    }
  }
}

void Tally::init_results()
{
  int n_scores = scores_.size() * nuclides_.size();
  results_ = xt::empty<double>({n_filter_bins_, n_scores, 3});
}

void Tally::reset()
{
  n_realizations_ = 0;
  if (results_.size() != 0) {
    xt::view(results_, xt::all()) = 0.0;
  }
}

void Tally::accumulate()
{
  // Increment number of realizations
  n_realizations_ += settings::reduce_tallies ? 1 : mpi::n_procs;

  if (mpi::master || !settings::reduce_tallies) {
    // Calculate total source strength for normalization
    double total_source = 0.0;
    if (settings::run_mode == RunMode::FIXED_SOURCE) {
      total_source = model::external_sources_probability.integral();
    } else {
      total_source = 1.0;
    }

    // Account for number of source particles in normalization
    double norm =
      total_source / (settings::n_particles * settings::gen_per_batch);

    if (settings::solver_type == SolverType::RANDOM_RAY) {
      norm = 1.0;
    }

// Accumulate each result
#pragma omp parallel for
    for (int i = 0; i < results_.shape()[0]; ++i) {
      for (int j = 0; j < results_.shape()[1]; ++j) {
        double val = results_(i, j, TallyResult::VALUE) * norm;
        results_(i, j, TallyResult::VALUE) = 0.0;
        results_(i, j, TallyResult::SUM) += val;
        results_(i, j, TallyResult::SUM_SQ) += val * val;
      }
    }
  }
}

int Tally::score_index(const std::string& score) const
{
  for (int i = 0; i < scores_.size(); i++) {
    if (this->score_name(i) == score)
      return i;
  }
  return -1;
}

xt::xarray<double> Tally::get_reshaped_data() const
{
  std::vector<uint64_t> shape;
  for (auto f : filters()) {
    shape.push_back(model::tally_filters[f]->n_bins());
  }

  // add number of scores and nuclides to tally
  shape.push_back(results_.shape()[1]);
  shape.push_back(results_.shape()[2]);

  xt::xarray<double> reshaped_results = results_;
  reshaped_results.reshape(shape);
  return reshaped_results;
}

std::string Tally::score_name(int score_idx) const
{
  if (score_idx < 0 || score_idx >= scores_.size()) {
    fatal_error("Index in scores array is out of bounds.");
  }
  return reaction_name(scores_[score_idx]);
}

std::vector<std::string> Tally::scores() const
{
  std::vector<std::string> score_names;
  for (int score : scores_)
    score_names.push_back(reaction_name(score));
  return score_names;
}

std::string Tally::nuclide_name(int nuclide_idx) const
{
  if (nuclide_idx < 0 || nuclide_idx >= nuclides_.size()) {
    fatal_error("Index in nuclides array is out of bounds");
  }

  int nuclide = nuclides_.at(nuclide_idx);
  if (nuclide == -1) {
    return "total";
  }
  return data::nuclides.at(nuclide)->name_;
}

//==============================================================================
// Non-member functions
//==============================================================================

void read_tallies_xml()
{
  // Check if tallies.xml exists. If not, just return since it is optional
  std::string filename = settings::path_input + "tallies.xml";
  if (!file_exists(filename))
    return;

  write_message("Reading tallies XML file...", 5);

  // Parse tallies.xml file
  pugi::xml_document doc;
  doc.load_file(filename.c_str());
  pugi::xml_node root = doc.document_element();

  read_tallies_xml(root);
}

void read_tallies_xml(pugi::xml_node root)
{
  // Check for <assume_separate> setting
  if (check_for_node(root, "assume_separate")) {
    settings::assume_separate = get_node_value_bool(root, "assume_separate");
  }

  // Check for user meshes and allocate
  read_meshes(root);

  // We only need the mesh info for plotting
  if (settings::run_mode == RunMode::PLOTTING)
    return;

  // Read data for tally derivatives
  read_tally_derivatives(root);

  // ==========================================================================
  // READ FILTER DATA

  // Check for user filters and allocate
  for (auto node_filt : root.children("filter")) {
    auto f = Filter::create(node_filt);
  }

  // ==========================================================================
  // READ TALLY DATA

  // Check for user tallies
  int n = 0;
  for (auto node : root.children("tally"))
    ++n;
  if (n == 0 && mpi::master) {
    warning("No tallies present in tallies.xml file.");
  }

  for (auto node_tal : root.children("tally")) {
    model::tallies.push_back(make_unique<Tally>(node_tal));
  }
}

#ifdef OPENMC_MPI
void reduce_tally_results()
{
  // Don't reduce tally is no_reduce option is on
  if (settings::reduce_tallies) {
    for (int i_tally : model::active_tallies) {
      // Skip any tallies that are not active
      auto& tally {model::tallies[i_tally]};

      // Get view of accumulated tally values
      auto values_view = xt::view(tally->results_, xt::all(), xt::all(),
        static_cast<int>(TallyResult::VALUE));

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
  }

  // Note that global tallies are *always* reduced even when no_reduce option is
  // on.

  // Get view of global tally values
  auto& gt = simulation::global_tallies;
  auto gt_values_view =
    xt::view(gt, xt::all(), static_cast<int>(TallyResult::VALUE));

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
  if (mpi::master)
    simulation::total_weight = weight_reduced;
}
#endif

void accumulate_tallies()
{
#ifdef OPENMC_MPI
  // Combine tally results onto master process
  if (mpi::n_procs > 1 && settings::solver_type == SolverType::MONTE_CARLO) {
    reduce_tally_results();
  }
#endif

  // Increase number of realizations (only used for global tallies)
  simulation::n_realizations += 1;

  // Accumulate on master only unless run is not reduced then do it on all
  if (mpi::master || !settings::reduce_tallies) {
    auto& gt = simulation::global_tallies;

    if (settings::run_mode == RunMode::EIGENVALUE) {
      if (simulation::current_batch > settings::n_inactive) {
        // Accumulate products of different estimators of k
        double k_col = gt(GlobalTally::K_COLLISION, TallyResult::VALUE) /
                       simulation::total_weight;
        double k_abs = gt(GlobalTally::K_ABSORPTION, TallyResult::VALUE) /
                       simulation::total_weight;
        double k_tra = gt(GlobalTally::K_TRACKLENGTH, TallyResult::VALUE) /
                       simulation::total_weight;
        simulation::k_col_abs += k_col * k_abs;
        simulation::k_col_tra += k_col * k_tra;
        simulation::k_abs_tra += k_abs * k_tra;
      }
    }

    // Accumulate results for global tallies
    for (int i = 0; i < N_GLOBAL_TALLIES; ++i) {
      double val = gt(i, TallyResult::VALUE) / simulation::total_weight;
      gt(i, TallyResult::VALUE) = 0.0;
      gt(i, TallyResult::SUM) += val;
      gt(i, TallyResult::SUM_SQ) += val * val;
    }
  }

  // Accumulate results for each tally
  for (int i_tally : model::active_tallies) {
    auto& tally {model::tallies[i_tally]};
    tally->accumulate();
  }
}

void setup_active_tallies()
{
  model::active_tallies.clear();
  model::active_analog_tallies.clear();
  model::active_tracklength_tallies.clear();
  model::active_collision_tallies.clear();
  model::active_meshsurf_tallies.clear();
  model::active_surface_tallies.clear();
  model::active_pulse_height_tallies.clear();
  model::time_grid.clear();

  for (auto i = 0; i < model::tallies.size(); ++i) {
    const auto& tally {*model::tallies[i]};

    if (tally.active_) {
      model::active_tallies.push_back(i);
      switch (tally.type_) {

      case TallyType::VOLUME:
        switch (tally.estimator_) {
        case TallyEstimator::ANALOG:
          model::active_analog_tallies.push_back(i);
          break;
        case TallyEstimator::TRACKLENGTH:
          model::active_tracklength_tallies.push_back(i);
          if (auto time_filter = tally->get_filter<TimeFilter>()) {
            model::add_to_time_grid(time_filter.bins_);
          }
          break;
        case TallyEstimator::COLLISION:
          model::active_collision_tallies.push_back(i);
        }
        break;

      case TallyType::MESH_SURFACE:
        model::active_meshsurf_tallies.push_back(i);
        break;

      case TallyType::SURFACE:
        model::active_surface_tallies.push_back(i);
        break;

      case TallyType::PULSE_HEIGHT:
        model::active_pulse_height_tallies.push_back(i);
        break;
      }
    }
  }
}

void free_memory_tally()
{
  model::tally_derivs.clear();
  model::tally_deriv_map.clear();

  model::tally_filters.clear();
  model::filter_map.clear();

  model::tallies.clear();

  model::active_tallies.clear();
  model::active_analog_tallies.clear();
  model::active_tracklength_tallies.clear();
  model::active_collision_tallies.clear();
  model::active_meshsurf_tallies.clear();
  model::active_surface_tallies.clear();
  model::active_pulse_height_tallies.clear();

  model::tally_map.clear();
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int openmc_extend_tallies(
  int32_t n, int32_t* index_start, int32_t* index_end)
{
  if (index_start)
    *index_start = model::tallies.size();
  if (index_end)
    *index_end = model::tallies.size() + n - 1;
  for (int i = 0; i < n; ++i) {
    model::tallies.push_back(make_unique<Tally>(-1));
  }
  return 0;
}

extern "C" int openmc_get_tally_index(int32_t id, int32_t* index)
{
  auto it = model::tally_map.find(id);
  if (it == model::tally_map.end()) {
    set_errmsg(fmt::format("No tally exists with ID={}.", id));
    return OPENMC_E_INVALID_ID;
  }

  *index = it->second;
  return 0;
}

extern "C" void openmc_get_tally_next_id(int32_t* id)
{
  int32_t largest_tally_id = 0;
  for (const auto& t : model::tallies) {
    largest_tally_id = std::max(largest_tally_id, t->id_);
  }
  *id = largest_tally_id + 1;
}

extern "C" int openmc_tally_get_estimator(int32_t index, int* estimator)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *estimator = static_cast<int>(model::tallies[index]->estimator_);
  return 0;
}

extern "C" int openmc_tally_set_estimator(int32_t index, const char* estimator)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  auto& t {model::tallies[index]};

  std::string est = estimator;
  if (est == "analog") {
    t->estimator_ = TallyEstimator::ANALOG;
  } else if (est == "collision") {
    t->estimator_ = TallyEstimator::COLLISION;
  } else if (est == "tracklength") {
    t->estimator_ = TallyEstimator::TRACKLENGTH;
  } else {
    set_errmsg("Unknown tally estimator: " + est);
    return OPENMC_E_INVALID_ARGUMENT;
  }
  return 0;
}

extern "C" int openmc_tally_get_id(int32_t index, int32_t* id)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *id = model::tallies[index]->id_;
  return 0;
}

extern "C" int openmc_tally_set_id(int32_t index, int32_t id)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  model::tallies[index]->set_id(id);
  return 0;
}

extern "C" int openmc_tally_get_type(int32_t index, int32_t* type)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  *type = static_cast<int>(model::tallies[index]->type_);

  return 0;
}

extern "C" int openmc_tally_set_type(int32_t index, const char* type)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  if (strcmp(type, "volume") == 0) {
    model::tallies[index]->type_ = TallyType::VOLUME;
  } else if (strcmp(type, "mesh-surface") == 0) {
    model::tallies[index]->type_ = TallyType::MESH_SURFACE;
  } else if (strcmp(type, "surface") == 0) {
    model::tallies[index]->type_ = TallyType::SURFACE;
  } else if (strcmp(type, "pulse-height") == 0) {
    model::tallies[index]->type_ = TallyType::PULSE_HEIGHT;
  } else {
    set_errmsg(fmt::format("Unknown tally type: {}", type));
    return OPENMC_E_INVALID_ARGUMENT;
  }

  return 0;
}

extern "C" int openmc_tally_get_active(int32_t index, bool* active)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  *active = model::tallies[index]->active_;

  return 0;
}

extern "C" int openmc_tally_set_active(int32_t index, bool active)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  model::tallies[index]->active_ = active;

  return 0;
}

extern "C" int openmc_tally_get_writable(int32_t index, bool* writable)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  *writable = model::tallies[index]->writable();

  return 0;
}

extern "C" int openmc_tally_set_writable(int32_t index, bool writable)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  model::tallies[index]->set_writable(writable);

  return 0;
}

extern "C" int openmc_tally_get_multiply_density(int32_t index, bool* value)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  *value = model::tallies[index]->multiply_density();

  return 0;
}

extern "C" int openmc_tally_set_multiply_density(int32_t index, bool value)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  model::tallies[index]->set_multiply_density(value);

  return 0;
}

extern "C" int openmc_tally_get_scores(int32_t index, int** scores, int* n)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *scores = model::tallies[index]->scores_.data();
  *n = model::tallies[index]->scores_.size();
  return 0;
}

extern "C" int openmc_tally_set_scores(
  int32_t index, int n, const char** scores)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  vector<std::string> scores_str(scores, scores + n);
  try {
    model::tallies[index]->set_scores(scores_str);
  } catch (const std::invalid_argument& ex) {
    set_errmsg(ex.what());
    return OPENMC_E_INVALID_ARGUMENT;
  }

  return 0;
}

extern "C" int openmc_tally_get_nuclides(int32_t index, int** nuclides, int* n)
{
  // Make sure the index fits in the array bounds.
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *n = model::tallies[index]->nuclides_.size();
  *nuclides = model::tallies[index]->nuclides_.data();

  return 0;
}

extern "C" int openmc_tally_set_nuclides(
  int32_t index, int n, const char** nuclides)
{
  // Make sure the index fits in the array bounds.
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  vector<std::string> words(nuclides, nuclides + n);
  vector<int> nucs;
  for (auto word : words) {
    if (word == "total") {
      nucs.push_back(-1);
    } else {
      auto search = data::nuclide_map.find(word);
      if (search == data::nuclide_map.end()) {
        int err = openmc_load_nuclide(word.c_str(), nullptr, 0);
        if (err < 0) {
          set_errmsg(openmc_err_msg);
          return OPENMC_E_DATA;
        }
      }
      nucs.push_back(data::nuclide_map.at(word));
    }
  }

  model::tallies[index]->nuclides_ = nucs;

  return 0;
}

extern "C" int openmc_tally_get_filters(
  int32_t index, const int32_t** indices, size_t* n)
{
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *indices = model::tallies[index]->filters().data();
  *n = model::tallies[index]->filters().size();
  return 0;
}

extern "C" int openmc_tally_set_filters(
  int32_t index, size_t n, const int32_t* indices)
{
  // Make sure the index fits in the array bounds.
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  // Set the filters.
  try {
    // Convert indices to filter pointers
    vector<Filter*> filters;
    for (int64_t i = 0; i < n; ++i) {
      int32_t i_filt = indices[i];
      filters.push_back(model::tally_filters.at(i_filt).get());
    }
    model::tallies[index]->set_filters(filters);
  } catch (const std::out_of_range& ex) {
    set_errmsg("Index in tally filter array out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  return 0;
}

//! Reset tally results and number of realizations
extern "C" int openmc_tally_reset(int32_t index)
{
  // Make sure the index fits in the array bounds.
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  model::tallies[index]->reset();
  return 0;
}

extern "C" int openmc_tally_get_n_realizations(int32_t index, int32_t* n)
{
  // Make sure the index fits in the array bounds.
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *n = model::tallies[index]->n_realizations_;
  return 0;
}

//! \brief Returns a pointer to a tally results array along with its shape. This
//! allows a user to obtain in-memory tally results from Python directly.
extern "C" int openmc_tally_results(
  int32_t index, double** results, size_t* shape)
{
  // Make sure the index fits in the array bounds.
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  const auto& t {model::tallies[index]};
  if (t->results_.size() == 0) {
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

extern "C" int openmc_global_tallies(double** ptr)
{
  *ptr = simulation::global_tallies.data();
  return 0;
}

extern "C" size_t tallies_size()
{
  return model::tallies.size();
}

// given a tally ID, remove it from the tallies vector
extern "C" int openmc_remove_tally(int32_t index)
{
  // check that id is in the map
  if (index < 0 || index >= model::tallies.size()) {
    set_errmsg("Index in tallies array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  // delete the tally via iterator pointing to correct position
  // this calls the Tally destructor, removing the tally from the map as well
  model::tallies.erase(model::tallies.begin() + index);

  return 0;
}

} // namespace openmc
