#include "openmc/weight_windows.h"

#include "xtensor/xview.hpp"

#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/physics_common.h"
#include "openmc/search.h"
#include "openmc/tallies/tally.h"
#include "openmc/xml_interface.h"

#include <fmt/core.h>
#include <gsl/gsl-lite.hpp>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace variance_reduction {

std::unordered_map<int32_t, int32_t> ww_map;
openmc::vector<unique_ptr<WeightWindows>> weight_windows;

} // namespace variance_reduction

//==============================================================================
// Non-member functions
//==============================================================================

void apply_weight_windows(Particle& p)
{
  // skip dead or no energy
  if (p.E() <= 0 || !p.alive())
    return;

  bool in_domain = false;
  // TODO: this is a linear search - should do something more clever
  WeightWindow weight_window;
  for (const auto& ww : variance_reduction::weight_windows) {
    weight_window = ww->get_weight_window(p);
    if (weight_window.is_valid())
      break;
  }
  // particle is not in any of the ww domains, do nothing
  if (!weight_window.is_valid())
    return;

  // get the paramters
  double weight = p.wgt();

  // first check to see if particle should be killed for weight cutoff
  if (p.wgt() < weight_window.weight_cutoff) {
    p.wgt() = 0.0;
    return;
  }

  // check if particle is far above current weight window
  // only do this if the factor is not already set on the particle and a
  // maximum lower bound ratio is specified
  if (p.ww_factor() == 0.0 && weight_window.max_lb_ratio > 1.0 &&
      p.wgt() > weight_window.lower_weight * weight_window.max_lb_ratio) {
    p.ww_factor() =
      p.wgt() / (weight_window.lower_weight * weight_window.max_lb_ratio);
  }

  // move weight window closer to the particle weight if needed
  if (p.ww_factor() > 1.0)
    weight_window.scale(p.ww_factor());

  // if particle's weight is above the weight window split until they are within
  // the window
  if (weight > weight_window.upper_weight) {
    // do not further split the particle if above the limit
    if (p.n_split() >= settings::max_splits)
      return;

    double n_split = std::ceil(weight / weight_window.upper_weight);
    double max_split = weight_window.max_split;
    n_split = std::min(n_split, max_split);

    p.n_split() += n_split;

    // Create secondaries and divide weight among all particles
    int i_split = std::round(n_split);
    for (int l = 0; l < i_split - 1; l++) {
      p.create_secondary(weight / n_split, p.u(), p.E(), p.type());
    }
    // remaining weight is applied to current particle
    p.wgt() = weight / n_split;

  } else if (weight <= weight_window.lower_weight) {
    // if the particle weight is below the window, play Russian roulette
    double weight_survive =
      std::min(weight * weight_window.max_split, weight_window.survival_weight);
    russian_roulette(p, weight_survive);
  } // else particle is in the window, continue as normal
}

void free_memory_weight_windows()
{
  variance_reduction::ww_map.clear();
  variance_reduction::weight_windows.clear();
}

//==============================================================================
// WeightWindowSettings implementation
//==============================================================================

WeightWindows::WeightWindows(pugi::xml_node node)
{
  // Make sure required elements are present
  const vector<std::string> required_elems {"id", "particle_type",
    "energy_bounds", "lower_ww_bounds", "upper_ww_bounds"};
  for (const auto& elem : required_elems) {
    if (!check_for_node(node, elem.c_str())) {
      fatal_error(fmt::format("Must specify <{}> for weight windows.", elem));
    }
  }

  // Get weight windows ID
  int32_t id = std::stoi(get_node_value(node, "id"));
  this->set_id(id);

  // get the particle type
  auto particle_type_str = std::string(get_node_value(node, "particle_type"));
  particle_type_ = openmc::str_to_particle_type(particle_type_str);

  // Determine associated mesh
  int32_t mesh_id = std::stoi(get_node_value(node, "mesh"));
  mesh_idx_ = model::mesh_map.at(mesh_id);

  // energy bounds
  energy_bounds_ = get_node_array<double>(node, "energy_bounds");

  // get the survival value - optional
  if (check_for_node(node, "survival_ratio")) {
    survival_ratio_ = std::stod(get_node_value(node, "survival_ratio"));
    if (survival_ratio_ <= 1)
      fatal_error("Survival to lower weight window ratio must bigger than 1 "
                  "and less than the upper to lower weight window ratio.");
  }

  // get the max lower bound ratio - optional
  if (check_for_node(node, "max_lower_bound_ratio")) {
    max_lb_ratio_ = std::stod(get_node_value(node, "max_lower_bound_ratio"));
    if (max_lb_ratio_ < 1.0) {
      fatal_error("Maximum lower bound ratio must be larger than 1");
    }
  }

  // get the max split - optional
  if (check_for_node(node, "max_split")) {
    max_split_ = std::stod(get_node_value(node, "max_split"));
    if (max_split_ <= 1)
      fatal_error("max split must be larger than 1");
  }

  // weight cutoff - optional
  if (check_for_node(node, "weight_cutoff")) {
    weight_cutoff_ = std::stod(get_node_value(node, "weight_cutoff"));
    if (weight_cutoff_ <= 0)
      fatal_error("weight_cutoff must be larger than 0");
    if (weight_cutoff_ > 1)
      fatal_error("weight_cutoff must be less than 1");
  }

  // read the lower/upper weight bounds
  this->set_weight_windows(get_node_array<double>(node, "lower_ww_bounds"),
    get_node_array<double>(node, "upper_ww_bounds"));
}

void WeightWindows::set_id(int32_t id)
{
  Expects(id >= 0 || id == C_NONE);

  // Clear entry in mesh map in case one was already assigned
  if (id_ != C_NONE) {
    variance_reduction::ww_map.erase(id_);
    id_ = C_NONE;
  }

  // Ensure no other mesh has the same ID
  if (variance_reduction::ww_map.find(id) != variance_reduction::ww_map.end()) {
    throw std::runtime_error {
      fmt::format("Two weight windows have the same ID: {}", id)};
  }

  // If no ID is specified, auto-assign the next ID in the sequence
  if (id == C_NONE) {
    id = 0;
    for (const auto& m : variance_reduction::weight_windows) {
      id = std::max(id, m->id_);
    }
    ++id;
  }

  // Update ID and entry in the mesh map
  id_ = id;
  variance_reduction::ww_map[id] =
    variance_reduction::weight_windows.size() - 1;
}

WeightWindow WeightWindows::get_weight_window(const Particle& p) const
{
  // check for particle type
  if (particle_type_ != p.type()) {
    return {};
  }

  // Get mesh index for particle's position
  const auto& mesh = this->mesh();
  int ww_index = mesh.get_bin(p.r());

  // particle is outside the weight window mesh
  if (ww_index < 0)
    return {};

  // particle energy
  double E = p.E();

  // check to make sure energy is in range, expects sorted energy values
  if (E < energy_bounds_.front() || E > energy_bounds_.back())
    return {};

  // get the mesh bin in energy group
  int energy_bin =
    lower_bound_index(energy_bounds_.begin(), energy_bounds_.end(), E);

  // indices now points to the correct weight for the given energy
  ww_index += energy_bin * mesh.n_bins();

  // Create individual weight window
  WeightWindow ww;
  ww.lower_weight = lower_ww_[ww_index];
  ww.upper_weight = upper_ww_[ww_index];
  ww.survival_weight = ww.lower_weight * survival_ratio_;
  ww.max_lb_ratio = max_lb_ratio_;
  ww.max_split = max_split_;
  ww.weight_cutoff = weight_cutoff_;
  return ww;
}

double WeightWindows::bounds_size() const
{
  int num_spatial_bins = this->mesh().n_bins();
  int num_energy_bins = energy_bounds_.size() - 1;
  return num_spatial_bins * num_energy_bins;
}

void WeightWindows::set_weight_windows(
  gsl::span<const double> lower_bounds, gsl::span<const double> upper_bounds)
{
  // clear out old memory
  lower_ww_.clear();
  upper_ww_.clear();

  // set new weight window values
  lower_ww_.insert(lower_ww_.begin(), lower_bounds.begin(), lower_bounds.end());
  upper_ww_.insert(upper_ww_.begin(), upper_bounds.begin(), upper_bounds.end());

  // make sure that the upper and lower bounds have the same size
  if (upper_ww_.size() != lower_ww_.size()) {
    fatal_error("The upper and lower weight window lengths do not match.");
  }

  // check that the number of weight window entries is correct
  if (lower_ww_.size() != this->bounds_size()) {
    int num_energy_bins = energy_bounds_.size() - 1;
    int num_spatial_bins = this->mesh().n_bins();
    auto err_msg =
      fmt::format("In weight window domain {} the number of spatial "
                  "energy/spatial bins ({}) does not match the number "
                  "of weight bins ({})",
        id_, num_energy_bins * num_spatial_bins, lower_ww_.size());
    fatal_error(err_msg);
  }
}

void WeightWindows::set_weight_windows(
  gsl::span<const double> lower_bounds, double bounds_ratio)
{

  this->set_weight_windows(lower_bounds, lower_bounds);

  for (auto& e : upper_ww_) {
    e *= bounds_ratio;
  }
}

void WeightWindows::to_hdf5(hid_t group) const
{
  hid_t ww_group = create_group(group, fmt::format("weight_windows {}", id_));

  write_dataset(
    ww_group, "particle_type", openmc::particle_type_to_str(particle_type_));
  write_dataset(ww_group, "energy_bounds", energy_bounds_);
  write_dataset(ww_group, "lower_ww_bounds", lower_ww_);
  write_dataset(ww_group, "upper_ww_bounds", upper_ww_);
  write_dataset(ww_group, "survival_ratio", survival_ratio_);
  write_dataset(ww_group, "max_lower_bound_ratio", max_lb_ratio_);
  write_dataset(ww_group, "max_split", max_split_);
  write_dataset(ww_group, "weight_cutoff", weight_cutoff_);
  write_dataset(ww_group, "mesh", this->mesh().id_);

  close_group(ww_group);
}

//==============================================================================
// C API
//==============================================================================

extern "C" int openmc_set_weight_windows(
  int ww_id, size_t n, const double* lower_bounds, const double* upper_bounds)
{

  // look up the weight windows object
  const auto& wws =
    variance_reduction::weight_windows.at(variance_reduction::ww_map.at(ww_id));

  // check length of arrays
  if (n != wws->bounds_size()) {
    set_errmsg(fmt::format(
      "Incorrect size for weight window bounds for domain {}", wws->id()));
    return OPENMC_E_INVALID_ARGUMENT;
  }

  // set bounds
  wws->set_weight_windows({lower_bounds, n}, {upper_bounds, n});

  return 0;
}

extern "C" int openmc_update_weight_windows(int tally_idx, int ww_idx,
  const char* score, const char* value, const char* method)
{

  // get the requested tally
  const auto& tally = model::tallies.at(tally_idx);

  // check for a
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);

  // check the tally filters. We currently require that the tally can only have
  // a MeshFilter and (optionally) an EnergyFilter
  if (tally->filters().size() > 2) {
    set_errmsg(
      fmt::format("More than 2 filters found on tally {}", tally->id()));
    return OPENMC_E_INVALID_SIZE;
  }

  // ensure the provided score is valid
  auto tally_scores = tally->scores();
  auto score_it = std::find(tally_scores.begin(), tally_scores.end(), score);
  if (score_it == tally->scores().end()) {
    set_errmsg(fmt::format(
      "The score '{}' could not be found on tally {}", score, tally->id()));
    return OPENMC_E_INVALID_ARGUMENT;
  }

  // gather information from the tally (assume one group for now)
  int score_idx = score_it - tally_scores.begin();
  // sum results over all nuclides and select results related to a single score
  auto score_view =
    xt::view(xt::sum(tally->results(), {2}), score_idx, TallyResult::SUM);

  // now generate new weight window values
  double score_max = *std::max(score_view.begin(), score_view.end());

  // normalize the score by the max value
  xt::xarray<double> lower_bounds(score_view / score_max);

  wws->set_weight_windows(lower_bounds, 5.0);

  return 0;
}

} // namespace openmc
