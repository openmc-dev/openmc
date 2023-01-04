#include "openmc/weight_windows.h"

#include <set>

#include "xtensor/xview.hpp"
#include "xtensor/xstrided_view.hpp"

#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/physics_common.h"
#include "openmc/search.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/filter_particle.h"
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

WeightWindows::WeightWindows(int32_t id)
{
  index_ = variance_reduction::weight_windows.size();
  set_id(id);
  set_defaults();
}

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
  if (check_for_node(node, "energy_bounds"))
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

  set_defaults();
}

WeightWindows::~WeightWindows()
{
  variance_reduction::ww_map.erase(id());
}

WeightWindows* WeightWindows::create(int32_t id)
{
  variance_reduction::weight_windows.push_back(make_unique<WeightWindows>());
  return variance_reduction::weight_windows.back().get();
}

void WeightWindows::set_defaults()
{

  // ensure default values are set
  if (energy_bounds_.size() == 0) {
    int p_type = static_cast<int>(particle_type_);
    energy_bounds_.push_back(data::energy_min[p_type]);
    energy_bounds_.push_back(data::energy_max[p_type]);
  }
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
  variance_reduction::ww_map[id] = index_;
}


void WeightWindows::set_energy_bounds(gsl::span<double const> bounds) {
  energy_bounds_.clear();
  energy_bounds_.insert(energy_bounds_.begin(), bounds.begin(), bounds.end());
  // TODO: check that sizes still make sense
}

void WeightWindows::set_mesh(int32_t mesh_idx)
{
  if (mesh_idx < 0 || mesh_idx > model::meshes.size())
    fatal_error(fmt::format("Could not find a mesh for index {}", mesh_idx));

  mesh_idx_ = mesh_idx;
}

void WeightWindows::set_mesh(const std::unique_ptr<Mesh>& mesh)
{
  set_mesh(mesh.get());
}

void WeightWindows::set_mesh(const Mesh* mesh)
{
  set_mesh(model::mesh_map[mesh->id_]);
}

WeightWindow WeightWindows::get_weight_window(const Particle& p) const
{
  // check for particle type
  if (particle_type_ != p.type()) {
    return {};
  }

  // Get mesh index for particle's position
  const auto& mesh = this->mesh();
  int ww_index = mesh->get_bin(p.r());

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
  ww_index += energy_bin * mesh->n_bins();

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
  int num_spatial_bins = this->mesh()->n_bins();
  int num_energy_bins =
    energy_bounds_.size() > 0 ? energy_bounds_.size() - 1 : 1;
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
    int num_energy_bins =
      energy_bounds_.size() > 0 ? energy_bounds_.size() - 1 : 1;
    int num_spatial_bins = this->mesh()->n_bins();
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

void WeightWindows::update_weight_windows(const std::unique_ptr<Tally>& tally,
                                          const std::string& score,
                                          const std::string& value,
                                          const std::string& method) {
  // set of allowed filters
  const std::set<FilterType> allowed_filters = {FilterType::MESH, FilterType::ENERGY, FilterType::PARTICLE};

  const auto filters = tally->filters();
  const auto filter_types = tally->filter_types();
  auto filter_indices = tally->filter_indices();

  // make sure that all filters are allowed
  for (auto f_type : filter_types) {
    if (allowed_filters.find(f_type) == allowed_filters.end()) {
      fatal_error(fmt::format("Invalid filter type '{}' found on tally "
                              "used for weight window generation.", filter_type_strings[f_type]));
    }
  }

  // gather information from the tally (assume one group for now)
  int score_index = tally->score_index(score);

  if (score_index == C_NONE) {
    fatal_error(fmt::format("Score '{}' specified for weight window generation is not present on tally {}.",
    score, tally->id()));
  }

  // check for a particle filter
  // TODO: Add method checking for filter type (returns index?)
  int particle_idx = 0;
  if (tally->has_filter(FilterType::PARTICLE)) {
    const auto& particle_filter = model::tally_filters[tally->filters(filter_indices[FilterType::PARTICLE])];

    // make sure there is data for the WW particle type available
    const auto& particles = dynamic_cast<ParticleFilter*>(particle_filter.get())->particles();
    auto particle_it = std::find(particles.begin(), particles.end(), this->particle_type_);
    if (particle_it == particles.end()) {
      warning(fmt::format("Particle type '{}' not present on Filter {} for Tally {} used to update WeightWindows {}",
      particle_type_to_str(this->particle_type_), particle_filter->id(), tally->id(), this->id()));
      return;
    }

    // filter out other particle data
    int particle_idx = particle_it - particles.begin();
    int particle_dim = filter_indices[FilterType::PARTICLE];
  } else {
    // create an entry in filter_indices for the particle dimension
    // this will default to the last dimension
    filter_indices[FilterType::PARTICLE] = 2;
  }

  ///////////////////////////////////
  // TODO: Move into Tally class
  // re-shape results into a new xtensor
  std::array<uint64_t, 3> shape = {1, 1, 1};
  for (int i = 0; i < filters.size(); i++) {
    const auto& f = model::tally_filters[tally->filters(i)];
    shape[i] = f->n_bins();
  }
  ///////////////////////////////////

  // sum results over all nuclides and select results related to a single score
  auto score_view =
    xt::view(tally->results(), xt::all(), score_index, TallyResult::SUM);

  // apply the score view results to this xtensor
  // TODO: this is a big copy. Might want to figure out a way to
  // reduce memory usage here.

  xt::xarray<double> lower_bounds(score_view);
  // xt::xtensor<double, 3> lower_bounds(shape);
  lower_bounds.reshape(shape);
  std::copy(score_view.begin(), score_view.end(), lower_bounds.begin());
  // lower_bounds = xt::zeros_like(score_view);
  //lower_bounds = score_view;

  // down-select particle data
  // auto p_view = xt::view(xt::roll(lower_bounds, filter_indices[FilterType::PARTICLE]), particle_idx);

  // // adjust index of other filters if needed
  // for (auto filter_type : allowed_filters) {
  //   if (filter_indices[filter_type] > filter_indices[FilterType::PARTICLE])
  //     filter_indices[filter_type] -= 1;
  // }

  // move energy axis to the front
  auto e_view = xt::roll(lower_bounds, filter_indices[FilterType::ENERGY]);

  int e_shape = shape[filter_indices[FilterType::ENERGY]];

  // std::cout << "Here" << std::endl;

  for (int e = 0; e < e_shape; e++) {
    // select all
    auto group_view = xt::view(e_view, e, xt::all(), xt::all());

    double group_max = *std::max_element(group_view.begin(), group_view.end());

    // normalize values in this energy group by the maximum value in the group
    group_view /= group_max;
  }

  this->set_weight_windows(lower_bounds, 5.0);
}

void WeightWindows::export_to_hdf5(const std::string& filename) const
{
  hid_t file_id = file_open(filename, 'w');

  // Write file type
  write_attribute(file_id, "filetype", "weight_windows");

  // Write revisiion number for state point file
  write_attribute(file_id, "version", VERSION_WEIGHT_WINDOWS);

  hid_t weight_windows_group =
    create_group(file_id, fmt::format("weight_windows", id_));

  // write instance information to file
  this->to_hdf5(weight_windows_group);

  close_group(weight_windows_group);

  file_close(file_id);
}

void WeightWindows::to_hdf5(hid_t group) const
{

  hid_t ww_group = create_group(group, fmt::format("weight_windows {}", id()));

  write_dataset(
    ww_group, "particle_type", openmc::particle_type_to_str(particle_type_));
  write_dataset(ww_group, "energy_bounds", energy_bounds_);
  write_dataset(ww_group, "lower_ww_bounds", lower_ww_);
  write_dataset(ww_group, "upper_ww_bounds", upper_ww_);
  write_dataset(ww_group, "survival_ratio", survival_ratio_);
  write_dataset(ww_group, "max_lower_bound_ratio", max_lb_ratio_);
  write_dataset(ww_group, "max_split", max_split_);
  write_dataset(ww_group, "weight_cutoff", weight_cutoff_);
  write_dataset(ww_group, "mesh", this->mesh()->id());

  close_group(ww_group);
}

//==============================================================================
// C API
//==============================================================================

extern "C" int openmc_get_weight_windows_index(int32_t id, int32_t* idx)
{
  auto it = variance_reduction::ww_map.find(id);
  if (it == variance_reduction::ww_map.end()) {
    set_errmsg(fmt::format("No weight windows exist with ID={}", id));
    return OPENMC_E_INVALID_ID;
  }

  *idx = it->second;
  return 0;
}

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

extern "C" int openmc_update_weight_windows(int32_t tally_idx, int32_t ww_idx,
  const char* score, const char* value, const char* method)
{
  // get the requested tally
  const auto& tally = model::tallies.at(tally_idx);

  // get the WeightWindows object
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);

  // check the tally filters. We currently require that the tally can only have
  // a MeshFilter and (optionally) an EnergyFilter
  if (tally->filters().size() > 2) {
    set_errmsg(
      fmt::format("More than 2 filters found on tally {}", tally->id()));
    return OPENMC_E_INVALID_SIZE;
  }

  // ensure the provided score is valid
  int score_idx = tally->score_index(score);
  if (score_idx == -1) {
    set_errmsg(fmt::format(
      "The score '{}' could not be found on tally {}", score, tally->id()));
    return OPENMC_E_INVALID_ARGUMENT;
  }

  wws->update_weight_windows(tally, score, value, method);

  return 0;
}

extern "C" int openmc_weight_windows_set_mesh(int32_t ww_idx, int32_t mesh_idx)
{
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);
  wws->set_mesh(mesh_idx);
  return 0;
}

extern "C" int openmc_weight_windows_get_mesh(int32_t ww_idx, int32_t* mesh_idx)
{
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);
  *mesh_idx = model::mesh_map.at(wws->mesh()->id());
  return 0;
}

extern "C" int openmc_weight_windows_set_energy_bounds(
  int32_t ww_idx, double* e_bounds, size_t e_bounds_size)
{
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);
  wws->set_energy_bounds({e_bounds, e_bounds_size});
  return 0;
}

extern "C" int openmc_weight_windows_get_energy_bounds(
  int32_t ww_idx, const double** e_bounds, size_t* e_bounds_size)
{

  const auto& wws = variance_reduction::weight_windows[ww_idx].get();
  *e_bounds = wws->energy_bounds().data();
  *e_bounds_size = wws->energy_bounds().size();
  return 0;
}

extern "C" int openmc_extend_weight_windows(
  int32_t n, int32_t* index_start, int32_t* index_end)
{
  if (index_start)
    *index_start = variance_reduction::weight_windows.size();
  if (index_end)
    *index_end = variance_reduction::weight_windows.size() + n - 1;
  for (int i = 0; i < n; ++i)
    variance_reduction::weight_windows.push_back(
      make_unique<WeightWindows>(-1));
  return 0;
}

extern "C" size_t openmc_weight_windows_size()
{
  return variance_reduction::weight_windows.size();
}

extern "C" void openmc_weight_windows_export(const char* filename)
{
  hid_t ww_file = file_open(filename, 'w');

  hid_t weight_windows_group = create_group(ww_file, "weight_windows");

  for (const auto& ww : variance_reduction::weight_windows)
    ww->to_hdf5(weight_windows_group);

  close_group(weight_windows_group);

  file_close(ww_file);
}

} // namespace openmc
