#include "openmc/weight_windows.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <string>

#include "openmc/tensor.h"

#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/physics_common.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/filter_meshangular.h"
#include "openmc/tallies/filter_particle.h"
#include "openmc/tallies/tally.h"
#include "openmc/xml_interface.h"

#include <fmt/core.h>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace variance_reduction {

std::unordered_map<int32_t, int32_t> ww_map;
openmc::vector<unique_ptr<WeightWindows>> weight_windows;
openmc::vector<unique_ptr<WeightWindowsGenerator>> weight_windows_generators;
openmc::vector<unique_ptr<SourceBias>> source_biases;

} // namespace variance_reduction

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
  const vector<std::string> required_elems {
    "id", "particle_type", "lower_ww_bounds", "upper_ww_bounds"};
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
  particle_type_ = ParticleType {particle_type_str};

  // Determine associated mesh
  int32_t mesh_id = std::stoi(get_node_value(node, "mesh"));
  set_mesh(model::mesh_map.at(mesh_id));

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
  this->set_bounds(get_node_array<double>(node, "lower_ww_bounds"),
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
  auto wws = variance_reduction::weight_windows.back().get();
  variance_reduction::ww_map[wws->id()] =
    variance_reduction::weight_windows.size() - 1;
  return wws;
}

WeightWindows* WeightWindows::from_hdf5(
  hid_t wws_group, const std::string& group_name)
{
  // collect ID from the name of this group
  hid_t ww_group = open_group(wws_group, group_name);

  auto wws = WeightWindows::create();

  std::string particle_type;
  read_dataset(ww_group, "particle_type", particle_type);
  wws->particle_type_ = ParticleType {particle_type};

  read_dataset<double>(ww_group, "energy_bounds", wws->energy_bounds_);

  int32_t mesh_id;
  read_dataset(ww_group, "mesh", mesh_id);

  if (model::mesh_map.count(mesh_id) == 0) {
    fatal_error(
      fmt::format("Mesh {} used in weight windows does not exist.", mesh_id));
  }
  wws->set_mesh(model::mesh_map[mesh_id]);

  wws->lower_ww_ =
    tensor::Tensor<double>({static_cast<size_t>(wws->bounds_size()[0]),
      static_cast<size_t>(wws->bounds_size()[1])});
  wws->upper_ww_ =
    tensor::Tensor<double>({static_cast<size_t>(wws->bounds_size()[0]),
      static_cast<size_t>(wws->bounds_size()[1])});

  read_dataset<double>(ww_group, "lower_ww_bounds", wws->lower_ww_);
  read_dataset<double>(ww_group, "upper_ww_bounds", wws->upper_ww_);
  read_dataset(ww_group, "survival_ratio", wws->survival_ratio_);
  read_dataset(ww_group, "max_lower_bound_ratio", wws->max_lb_ratio_);
  read_dataset(ww_group, "max_split", wws->max_split_);
  read_dataset(ww_group, "weight_cutoff", wws->weight_cutoff_);

  close_group(ww_group);

  return wws;
}

void WeightWindows::set_defaults()
{
  // set energy bounds to the min/max energy supported by the data
  if (energy_bounds_.size() == 0) {
    int p_type = particle_type_.transport_index();
    if (p_type == C_NONE) {
      fatal_error("Weight windows particle is not supported for transport.");
    }
    energy_bounds_.push_back(data::energy_min[p_type]);
    energy_bounds_.push_back(data::energy_max[p_type]);
  }
}

void WeightWindows::allocate_ww_bounds()
{
  auto shape = bounds_size();
  if (shape[0] * shape[1] == 0) {
    auto msg = fmt::format(
      "Size of weight window bounds is zero for WeightWindows {}", id());
    warning(msg);
  }
  lower_ww_ = tensor::Tensor<double>(
    {static_cast<size_t>(shape[0]), static_cast<size_t>(shape[1])});
  lower_ww_.fill(-1);
  upper_ww_ = tensor::Tensor<double>(
    {static_cast<size_t>(shape[0]), static_cast<size_t>(shape[1])});
  upper_ww_.fill(-1);
}

void WeightWindows::set_id(int32_t id)
{
  assert(id >= 0 || id == C_NONE);

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

void WeightWindows::set_energy_bounds(span<const double> bounds)
{
  energy_bounds_.clear();
  energy_bounds_.insert(energy_bounds_.begin(), bounds.begin(), bounds.end());
  // if the mesh is set, allocate space for weight window bounds
  if (mesh_idx_ != C_NONE)
    allocate_ww_bounds();
}

void WeightWindows::set_particle_type(ParticleType p_type)
{
  if (!p_type.is_neutron() && !p_type.is_photon())
    fatal_error(fmt::format(
      "Particle type '{}' cannot be applied to weight windows.", p_type.str()));
  particle_type_ = p_type;
}

void WeightWindows::set_mesh(int32_t mesh_idx)
{
  if (mesh_idx < 0 || mesh_idx >= model::meshes.size())
    fatal_error(fmt::format("Could not find a mesh for index {}", mesh_idx));

  mesh_idx_ = mesh_idx;
  model::meshes[mesh_idx_]->prepare_for_point_location();
  allocate_ww_bounds();
}

void WeightWindows::set_mesh(const std::unique_ptr<Mesh>& mesh)
{
  set_mesh(mesh.get());
}

void WeightWindows::set_mesh(const Mesh* mesh)
{
  set_mesh(model::mesh_map[mesh->id_]);
}

std::pair<bool, WeightWindow> WeightWindows::get_weight_window(
  const Particle& p) const
{
  // check for particle type
  if (particle_type_ != p.type()) {
    return {false, {}};
  }

  // particle energy
  double E = p.E();

  // check to make sure energy is in range, expects sorted energy values
  if (E < energy_bounds_.front() || E > energy_bounds_.back())
    return {false, {}};

  // Get mesh index for particle's position
  const auto& mesh = this->mesh();
  int mesh_bin = mesh->get_bin(p.r());

  // particle is outside the weight window mesh
  if (mesh_bin < 0)
    return {false, {}};

  // get the mesh bin in energy group
  int energy_bin =
    lower_bound_index(energy_bounds_.begin(), energy_bounds_.end(), E);

  // mesh_bin += energy_bin * mesh->n_bins();
  // Create individual weight window
  WeightWindow ww;
  ww.lower_weight = lower_ww_(energy_bin, mesh_bin);
  ww.upper_weight = upper_ww_(energy_bin, mesh_bin);
  ww.survival_weight = ww.lower_weight * survival_ratio_;
  ww.max_lb_ratio = max_lb_ratio_;
  ww.max_split = max_split_;
  ww.weight_cutoff = weight_cutoff_;
  return {true, ww};
}

std::array<int, 2> WeightWindows::bounds_size() const
{
  int num_spatial_bins = this->mesh()->n_bins();
  int num_energy_bins =
    energy_bounds_.size() > 0 ? energy_bounds_.size() - 1 : 1;
  return {num_energy_bins, num_spatial_bins};
}

template<class T>
void WeightWindows::check_bounds(const T& lower, const T& upper) const
{
  // make sure that the upper and lower bounds have the same size
  if (lower.size() != upper.size()) {
    auto msg = fmt::format("The upper and lower weight window lengths do not "
                           "match.\n Lower size: {}\n Upper size: {}",
      lower.size(), upper.size());
    fatal_error(msg);
  }
  this->check_bounds(lower);
}

template<class T>
void WeightWindows::check_bounds(const T& bounds) const
{
  // check that the number of weight window entries is correct
  auto dims = this->bounds_size();
  if (bounds.size() != dims[0] * dims[1]) {
    auto err_msg =
      fmt::format("In weight window domain {} the number of spatial "
                  "energy/spatial bins ({}) does not match the number "
                  "of weight bins ({})",
        id_, dims, bounds.size());
    fatal_error(err_msg);
  }
}

void WeightWindows::set_bounds(const tensor::Tensor<double>& lower_bounds,
  const tensor::Tensor<double>& upper_bounds)
{

  this->check_bounds(lower_bounds, upper_bounds);

  // set new weight window values
  lower_ww_ = lower_bounds;
  upper_ww_ = upper_bounds;
}

void WeightWindows::set_bounds(
  const tensor::Tensor<double>& lower_bounds, double ratio)
{
  this->check_bounds(lower_bounds);

  // set new weight window values
  lower_ww_ = lower_bounds;
  upper_ww_ = lower_bounds;
  upper_ww_ *= ratio;
}

void WeightWindows::set_bounds(
  span<const double> lower_bounds, span<const double> upper_bounds)
{
  check_bounds(lower_bounds, upper_bounds);
  auto shape = this->bounds_size();
  lower_ww_ = tensor::Tensor<double>(
    {static_cast<size_t>(shape[0]), static_cast<size_t>(shape[1])});
  upper_ww_ = tensor::Tensor<double>(
    {static_cast<size_t>(shape[0]), static_cast<size_t>(shape[1])});

  // Copy weight window values from input spans into the tensors
  std::copy(lower_bounds.data(), lower_bounds.data() + lower_ww_.size(),
    lower_ww_.data());
  std::copy(upper_bounds.data(), upper_bounds.data() + upper_ww_.size(),
    upper_ww_.data());
}

void WeightWindows::set_bounds(span<const double> lower_bounds, double ratio)
{
  this->check_bounds(lower_bounds);

  auto shape = this->bounds_size();
  lower_ww_ = tensor::Tensor<double>(
    {static_cast<size_t>(shape[0]), static_cast<size_t>(shape[1])});
  upper_ww_ = tensor::Tensor<double>(
    {static_cast<size_t>(shape[0]), static_cast<size_t>(shape[1])});

  // Copy lower bounds into both arrays, then scale upper by ratio
  std::copy(lower_bounds.data(), lower_bounds.data() + lower_ww_.size(),
    lower_ww_.data());
  std::copy(lower_bounds.data(), lower_bounds.data() + upper_ww_.size(),
    upper_ww_.data());
  upper_ww_ *= ratio;
}

void WeightWindows::update_weights(const Tally* tally, const std::string& value,
  double threshold, double ratio, WeightWindowUpdateMethod method)
{
  ///////////////////////////
  // Setup and checks
  ///////////////////////////
  this->check_tally_update_compatibility(tally);

  // Dimensions of weight window arrays
  int e_bins = lower_ww_.shape(0);
  int64_t mesh_bins = lower_ww_.shape(1);

  // Initialize weight window arrays to -1.0 by default
#pragma omp parallel for collapse(2) schedule(static)
  for (int e = 0; e < e_bins; e++) {
    for (int64_t m = 0; m < mesh_bins; m++) {
      lower_ww_(e, m) = -1.0;
      upper_ww_(e, m) = -1.0;
    }
  }

  // determine which value to use
  const std::set<std::string> allowed_values = {"mean", "rel_err"};
  if (allowed_values.count(value) == 0) {
    fatal_error(fmt::format("Invalid value '{}' specified for weight window "
                            "generation. Must be one of: 'mean' or 'rel_err'",
      value));
  }

  // determine the index of the specified score
  int score_index = tally->score_index("flux");
  if (score_index == C_NONE) {
    fatal_error(
      fmt::format("A 'flux' score required for weight window generation "
                  "is not present on tally {}.",
        tally->id()));
  }

  ///////////////////////////
  // Extract tally data
  //
  // At the end of this section, mean and rel_err are
  // 2D tensors of tally data (n_e_groups, n_mesh_bins)
  //
  ///////////////////////////

  // build a shape for the tally results, this will always be
  // dimension 5 (3 filter dimensions, 1 score dimension, 1 results dimension)
  // Look for the size of the last dimension of the results tensor
  const auto& results = tally->results();
  const int results_dim = static_cast<int>(results.shape(2));
  std::array<int, 5> shape = {1, 1, 1, tally->n_scores(), results_dim};

  // set the shape for the filters applied on the tally
  for (int i = 0; i < tally->filters().size(); i++) {
    const auto& filter = model::tally_filters[tally->filters(i)];
    shape[i] = filter->n_bins();
  }

  // build the transpose information to re-order data according to filter type
  std::array<int, 5> transpose = {0, 1, 2, 3, 4};

  // track our filter types and where we've added new ones
  std::vector<FilterType> filter_types = tally->filter_types();

  // assign other filter types to dummy positions if needed
  if (!tally->has_filter(FilterType::PARTICLE))
    filter_types.push_back(FilterType::PARTICLE);

  if (!tally->has_filter(FilterType::ENERGY))
    filter_types.push_back(FilterType::ENERGY);

  // particle axis mapping
  transpose[0] =
    std::find(filter_types.begin(), filter_types.end(), FilterType::PARTICLE) -
    filter_types.begin();

  // energy axis mapping
  transpose[1] =
    std::find(filter_types.begin(), filter_types.end(), FilterType::ENERGY) -
    filter_types.begin();

  // mesh axis mapping
  transpose[2] =
    std::find(filter_types.begin(), filter_types.end(), FilterType::MESH) -
    filter_types.begin();

  // determine the index of the particle within its filter
  int particle_idx = 0;
  if (tally->has_filter(FilterType::PARTICLE)) {
    auto pf = tally->get_filter<ParticleFilter>();
    const auto& particles = pf->particles();

    auto p_it =
      std::find(particles.begin(), particles.end(), this->particle_type_);
    if (p_it == particles.end()) {
      auto msg = fmt::format("Particle type '{}' not present on Filter {} for "
                             "Tally {} used to update WeightWindows {}",
        this->particle_type_.str(), pf->id(), tally->id(), this->id());
      fatal_error(msg);
    }

    particle_idx = p_it - particles.begin();
  }

  // The tally results array is 3D: (n_filter_combos, n_scores, n_result_types).
  // The first dimension is a row-major flattening of up to 3 filter dimensions
  // (particle, energy, mesh) whose storage order depends on which filters the
  // tally has. We need to map our desired indices (particle, energy, mesh)
  // into the correct flat filter combination index.
  //
  // transpose[i] tells us which storage position holds dimension i:
  //   i=0 -> particle, i=1 -> energy, i=2 -> mesh
  // shape[j] gives the number of bins for filter storage position j.

  // Row-major strides for the 3 filter dimensions
  const int stride0 = shape[1] * shape[2];
  const int stride1 = shape[2];

  tensor::Tensor<double> sum(
    {static_cast<size_t>(e_bins), static_cast<size_t>(mesh_bins)});
  tensor::Tensor<double> sum_sq(
    {static_cast<size_t>(e_bins), static_cast<size_t>(mesh_bins)});

  const int i_sum = static_cast<int>(TallyResult::SUM);
  const int i_sum_sq = static_cast<int>(TallyResult::SUM_SQ);

  for (int e = 0; e < e_bins; e++) {
    for (int64_t m = 0; m < mesh_bins; m++) {
      // Place particle, energy, and mesh indices into their storage positions
      std::array<int, 3> idx = {0, 0, 0};
      idx[transpose[0]] = particle_idx;
      idx[transpose[1]] = e;
      idx[transpose[2]] = static_cast<int>(m);

      // Compute flat filter combination index (row-major over filter dims)
      int flat = idx[0] * stride0 + idx[1] * stride1 + idx[2];

      sum(e, m) = results(flat, score_index, i_sum);
      sum_sq(e, m) = results(flat, score_index, i_sum_sq);
    }
  }
  int n = tally->n_realizations_;

  //////////////////////////////////////////////
  //
  // Assign new weight windows
  //
  // Use references to the existing weight window data
  // to store and update the values
  //
  //////////////////////////////////////////////

  // up to this point the data arrays are views into the tally results (no
  // computation has been performed) now we'll switch references to the tally's
  // bounds to avoid allocating additional memory
  auto& new_bounds = this->lower_ww_;
  auto& rel_err = this->upper_ww_;

  // get mesh volumes
  auto mesh_vols = this->mesh()->volumes();

  // Calculate mean (new_bounds) and relative error
#pragma omp parallel for collapse(2) schedule(static)
  for (int e = 0; e < e_bins; e++) {
    for (int64_t m = 0; m < mesh_bins; m++) {
      // Calculate mean
      new_bounds(e, m) = sum(e, m) / n;
      // Calculate relative error
      if (sum(e, m) > 0.0) {
        double mean_val = new_bounds(e, m);
        double variance = (sum_sq(e, m) / n - mean_val * mean_val) / (n - 1);
        rel_err(e, m) = std::sqrt(variance) / mean_val;
      } else {
        rel_err(e, m) = INFTY;
      }
      if (value == "rel_err") {
        new_bounds(e, m) = 1.0 / rel_err(e, m);
      }
    }
  }

  // Divide by volume of mesh elements
#pragma omp parallel for collapse(2) schedule(static)
  for (int e = 0; e < e_bins; e++) {
    for (int64_t m = 0; m < mesh_bins; m++) {
      new_bounds(e, m) /= mesh_vols[m];
    }
  }

  if (method == WeightWindowUpdateMethod::MAGIC) {
    // For MAGIC, weight windows are proportional to the forward fluxes.
    // We normalize weight windows independently for each energy group.

    // Find group maximum and normalize (per energy group)
    for (int e = 0; e < e_bins; e++) {
      double group_max = 0.0;

      // Find maximum value across all elements in this energy group
#pragma omp parallel for schedule(static) reduction(max : group_max)
      for (int64_t m = 0; m < mesh_bins; m++) {
        if (new_bounds(e, m) > group_max) {
          group_max = new_bounds(e, m);
        }
      }

      // Normalize values in this energy group by the maximum value
      if (group_max > 0.0) {
        double norm_factor = 1.0 / (2.0 * group_max);
#pragma omp parallel for schedule(static)
        for (int64_t m = 0; m < mesh_bins; m++) {
          new_bounds(e, m) *= norm_factor;
        }
      }
    }
  } else {
    // For (FW-)CADIS, weight windows are inversely proportional to the adjoint
    // fluxes. We normalize the weight windows across all energy groups.
#pragma omp parallel for collapse(2) schedule(static)
    for (int e = 0; e < e_bins; e++) {
      for (int64_t m = 0; m < mesh_bins; m++) {
        // Take the inverse, but are careful not to divide by zero
        if (new_bounds(e, m) != 0.0) {
          new_bounds(e, m) = 1.0 / new_bounds(e, m);
        } else {
          new_bounds(e, m) = 0.0;
        }
      }
    }

    // Find the maximum value across all elements
    double max_val = 0.0;
#pragma omp parallel for collapse(2) schedule(static) reduction(max : max_val)
    for (int e = 0; e < e_bins; e++) {
      for (int64_t m = 0; m < mesh_bins; m++) {
        if (new_bounds(e, m) > max_val) {
          max_val = new_bounds(e, m);
        }
      }
    }

    // Parallel normalization
    if (max_val > 0.0) {
      double norm_factor = 1.0 / (2.0 * max_val);
#pragma omp parallel for collapse(2) schedule(static)
      for (int e = 0; e < e_bins; e++) {
        for (int64_t m = 0; m < mesh_bins; m++) {
          new_bounds(e, m) *= norm_factor;
        }
      }
    }
  }

  // Final processing
#pragma omp parallel for collapse(2) schedule(static)
  for (int e = 0; e < e_bins; e++) {
    for (int64_t m = 0; m < mesh_bins; m++) {
      // Values where the mean is zero should be ignored
      if (sum(e, m) <= 0.0) {
        new_bounds(e, m) = -1.0;
      }
      // Values where the relative error is higher than the threshold should be
      // ignored
      else if (rel_err(e, m) > threshold) {
        new_bounds(e, m) = -1.0;
      }
      // Set the upper bounds
      upper_ww_(e, m) = ratio * lower_ww_(e, m);
    }
  }
}

void WeightWindows::check_tally_update_compatibility(const Tally* tally)
{
  // define the set of allowed filters for the tally
  const std::set<FilterType> allowed_filters = {
    FilterType::MESH, FilterType::ENERGY, FilterType::PARTICLE};

  // retrieve a mapping of filter type to filter index for the tally
  auto filter_indices = tally->filter_indices();

  // a mesh filter is required for a tally used to update weight windows
  if (!filter_indices.count(FilterType::MESH)) {
    fatal_error(
      "A mesh filter is required for a tally to update weight window bounds");
  }

  // ensure the mesh filter is using the same mesh as this weight window object
  auto mesh_filter = tally->get_filter<MeshFilter>();

  // make sure that all of the filters present on the tally are allowed
  for (auto filter_pair : filter_indices) {
    if (allowed_filters.find(filter_pair.first) == allowed_filters.end()) {
      fatal_error(fmt::format("Invalid filter type '{}' found on tally "
                              "used for weight window generation.",
        model::tally_filters[tally->filters(filter_pair.second)]->type_str()));
    }
  }

  if (mesh_filter->mesh() != mesh_idx_) {
    int32_t mesh_filter_id = model::meshes[mesh_filter->mesh()]->id();
    int32_t ww_mesh_id = model::meshes[this->mesh_idx_]->id();
    fatal_error(fmt::format("Mesh filter {} uses a different mesh ({}) than "
                            "weight window {} mesh ({})",
      mesh_filter->id(), mesh_filter_id, id_, ww_mesh_id));
  }

  // if an energy filter exists, make sure the energy grid matches that of this
  // weight window object
  if (auto energy_filter = tally->get_filter<EnergyFilter>()) {
    std::vector<double> filter_bins = energy_filter->bins();
    std::set<double> filter_e_bounds(
      energy_filter->bins().begin(), energy_filter->bins().end());
    if (filter_e_bounds.size() != energy_bounds().size()) {
      fatal_error(
        fmt::format("Energy filter {} does not have the same number of energy "
                    "bounds ({}) as weight window object {} ({})",
          energy_filter->id(), filter_e_bounds.size(), id_,
          energy_bounds().size()));
    }

    for (auto e : energy_bounds()) {
      if (filter_e_bounds.count(e) == 0) {
        fatal_error(fmt::format(
          "Energy bounds of filter {} and weight windows {} do not match",
          energy_filter->id(), id_));
      }
    }
  }
}

void WeightWindows::to_hdf5(hid_t group) const
{
  hid_t ww_group = create_group(group, fmt::format("weight_windows_{}", id()));

  write_dataset(ww_group, "mesh", this->mesh()->id());
  write_dataset(ww_group, "particle_type", particle_type_.str());
  write_dataset(ww_group, "energy_bounds", energy_bounds_);
  write_dataset(ww_group, "lower_ww_bounds", lower_ww_);
  write_dataset(ww_group, "upper_ww_bounds", upper_ww_);
  write_dataset(ww_group, "survival_ratio", survival_ratio_);
  write_dataset(ww_group, "max_lower_bound_ratio", max_lb_ratio_);
  write_dataset(ww_group, "max_split", max_split_);
  write_dataset(ww_group, "weight_cutoff", weight_cutoff_);

  close_group(ww_group);
}

//==============================================================================
// SourceBias implementation
//==============================================================================

SourceBias::SourceBias(int32_t spatial_mesh_idx, int32_t angle_mesh_idx,
  vector<double> energy_bounds, int32_t ww_id, std::string forward_source_mesh)
  : spatial_mesh_idx_(spatial_mesh_idx), angle_mesh_idx_(angle_mesh_idx),
    energy_bounds_(std::move(energy_bounds)), ww_id_(ww_id),
    forward_source_mesh_(forward_source_mesh)
{
  int64_t spatial_bins = model::meshes[spatial_mesh_idx_]->n_bins();
  int64_t angle_bins =
    angle_mesh_idx_ == C_NONE ? 1 : model::meshes[angle_mesh_idx_]->n_bins();
  int64_t energy_bins = energy_bounds_.empty() ? 1 : energy_bounds_.size() - 1;

  flux_ = tensor::Tensor<double>({static_cast<size_t>(spatial_bins),
    static_cast<size_t>(angle_bins), static_cast<size_t>(energy_bins)});
}

void SourceBias::update(const Tally* tally)
{
  int score_index = tally->score_index("flux");
  if (score_index == C_NONE) {
    fatal_error(fmt::format("A 'flux' score is required on tally {} used "
                            "for source biasing.",
      tally->id()));
  }

  const int64_t spatial_bins = flux_.shape(0);
  const int64_t angle_bins = flux_.shape(1);
  const int64_t energy_bins = flux_.shape(2);

  const auto& results = tally->results();

  const int n_filters = static_cast<int>(tally->filters().size());
  std::vector<int64_t> filt_shape(n_filters);
  std::vector<int64_t> filt_stride(n_filters, 1);
  for (int i = 0; i < n_filters; ++i) {
    filt_shape[i] = model::tally_filters[tally->filters(i)]->n_bins();
  }
  for (int i = n_filters - 2; i >= 0; --i) {
    filt_stride[i] = filt_stride[i + 1] * filt_shape[i + 1];
  }

  std::vector<FilterType> filter_types = tally->filter_types();
  auto position_of = [&filter_types](FilterType type) -> int {
    auto it = std::find(filter_types.begin(), filter_types.end(), type);
    return it == filter_types.end()
             ? -1
             : static_cast<int>(it - filter_types.begin());
  };

  const int pos_mesh = position_of(FilterType::MESH);
  const int pos_energy = position_of(FilterType::ENERGY);
  const int pos_angle = position_of(FilterType::MESH_ANGULAR);

  if (pos_mesh == -1) {
    fatal_error(fmt::format(
      "Tally {} used for source biasing is missing a spatial mesh filter.",
      tally->id()));
  }
  if (angle_bins > 1 && pos_angle == -1) {
    fatal_error(
      fmt::format("Tally {} used for source biasing is missing the expected "
                  "mesh-angular filter.",
        tally->id()));
  }

  const int i_sum = static_cast<int>(TallyResult::SUM);
  const int64_t n = tally->n_realizations_;

#pragma omp parallel for collapse(3) schedule(static)
  for (int64_t m = 0; m < spatial_bins; ++m) {
    for (int64_t a = 0; a < angle_bins; ++a) {
      for (int64_t e = 0; e < energy_bins; ++e) {
        int64_t flat = m * filt_stride[pos_mesh];
        if (pos_energy != -1)
          flat += e * filt_stride[pos_energy];
        if (pos_angle != -1)
          flat += a * filt_stride[pos_angle];

        double raw_sum = results(flat, score_index, i_sum);
        flux_(m, a, e) = (n > 0 && raw_sum > 0.0) ? raw_sum / n : 0.0;
      }
    }
  }
}

void SourceBias::to_hdf5(hid_t group) const
{
  hid_t sb_group = create_group(group, fmt::format("source_bias_{}", ww_id_));

  write_dataset(
    sb_group, "spatial_mesh", model::meshes[spatial_mesh_idx_]->id());
  write_dataset(sb_group, "angle_mesh",
    angle_mesh_idx_ == C_NONE ? C_NONE : model::meshes[angle_mesh_idx_]->id());
  write_dataset(sb_group, "energy_bounds", energy_bounds_);
  // Get biased source strength B(r,Omega,g) = flux(r,Omega,g) *
  // unbiased_strength(r,Omega,g) and the per-voxel weight 1/flux(r,Omega,g)
  const int64_t spatial_bins = flux_.shape(0);
  const int64_t angle_bins = flux_.shape(1);
  const int64_t energy_bins = flux_.shape(2);

  tensor::Tensor<double> biased_strength({static_cast<size_t>(spatial_bins),
    static_cast<size_t>(angle_bins), static_cast<size_t>(energy_bins)});
  tensor::Tensor<double> weights({static_cast<size_t>(spatial_bins),
    static_cast<size_t>(angle_bins), static_cast<size_t>(energy_bins)});

  // Calculate biasing parameters and simultaneously check to see if any
  // voxels containing probability mass associated with the unbiased source
  // have not been sampled.
  double total_strength = 0.0;
  double total_missed_S = 0.0;
  for (int64_t m = 0; m < spatial_bins; ++m) {
    for (int64_t a = 0; a < angle_bins; ++a) {
      for (int64_t e = 0; e < energy_bins; ++e) {
        double psi = flux_(m, a, e);
        double b = psi * unbiased_strength_(m, a, e);
        biased_strength(m, a, e) = b;
        total_strength += b;
        if (psi <= 0.0) {
          total_missed_S += unbiased_strength_(m, a, e);
        }
      }
    }
  }

  if (total_missed_S > TINY_BIT) {
    fatal_error(fmt::format("Biased source does not sample whole support of "
                            "the unbiased source. Probability mass lost in "
                            "regions with nonpositive adjoint flux: {}.",
      total_missed_S));
  }

  // Rescale biased strengths so that they add up to probability 1.0.
  // This would be done anyway by the DiscreteIndex inside of CorrelatedSource
  // when sampling a phase space bin, but the normalization should be done
  // beforehand so that source's .strength() attribute still appears as 1.0.
  // Also apply sample weight of (total_strength / psi) so that total
  // probability mass still equals 1.0.
  for (int64_t m = 0; m < spatial_bins; ++m) {
    for (int64_t a = 0; a < angle_bins; ++a) {
      for (int64_t e = 0; e < energy_bins; ++e) {
        double psi = flux_(m, a, e);
        weights(m, a, e) = psi > 0.0 ? total_strength / psi : 0.0;
        biased_strength(m, a, e) /= total_strength;
      }
    }
  }

  write_dataset(sb_group, "biased_source_strength", biased_strength);
  write_dataset(sb_group, "weights", weights);

  close_group(sb_group);
}

void SourceBias::compute_unbiased_strength(int64_t n_samples_per_source)
{
  const int64_t spatial_bins = flux_.shape(0);
  const int64_t angle_bins = flux_.shape(1);
  const int64_t energy_bins = flux_.shape(2);

  unbiased_strength_ =
    tensor::Tensor<double>({static_cast<size_t>(spatial_bins),
      static_cast<size_t>(angle_bins), static_cast<size_t>(energy_bins)});
  for (int64_t m = 0; m < spatial_bins; ++m) {
    for (int64_t a = 0; a < angle_bins; ++a) {
      for (int64_t e = 0; e < energy_bins; ++e) {
        unbiased_strength_(m, a, e) = 0.0;
      }
    }
  }

  // A precomputed file takes priority over sampling model::external_sources
  if (file_exists(forward_source_mesh_)) {
    load_forward_source_mesh(forward_source_mesh_);
    return;
  }

  Mesh* spatial_mesh = model::meshes[spatial_mesh_idx_].get();
  Mesh* angle_mesh =
    angle_mesh_idx_ == C_NONE ? nullptr : model::meshes[angle_mesh_idx_].get();

  // Normalize by total strength across all forward external sources
  double total_strength = 0.0;
  for (const auto& source_ptr : model::external_sources) {
    auto* is = dynamic_cast<IndependentSource*>(source_ptr.get());
    if (!is) {
      fatal_error("Computing an unbiased source strength distribution "
                  "requires all external sources to be independent "
                  "sources.");
    }
    total_strength += is->strength();
  }
  if (total_strength <= 0.0) {
    fatal_error("Total external source strength must be positive to "
                "compute an unbiased source strength distribution.");
  }

  // Fixed local seed, not part of the sequence used for transport
  uint64_t seed = 1;

  tensor::Tensor<double> spatial_angle_counts(
    {static_cast<size_t>(spatial_bins), static_cast<size_t>(angle_bins)});

  for (const auto& source_ptr : model::external_sources) {
    Source* s = source_ptr.get();
    auto* is = dynamic_cast<IndependentSource*>(s);
    double relative_strength = is->strength() / total_strength;

    auto* energy_dist = dynamic_cast<Discrete*>(is->energy());
    if (!energy_dist) {
      fatal_error(
        "Source biasing requires all external sources to use a Discrete "
        "(multigroup) energy distribution, matching the requirement for "
        "random ray fixed source problems.");
    }
    vector<double> group_prob(energy_bins, 0.0);
    const auto& e_vals = energy_dist->x();
    const auto& e_probs = energy_dist->prob();
    for (std::size_t i = 0; i < e_vals.size(); ++i) {
      int g = data::mg.get_group_index(e_vals[i]);
      if (g >= 0 && g < energy_bins)
        group_prob[g] += e_probs[i];
    }

    // Monte Carlo estimate of the joint (spatial, angular) distribution.
    // Energy is not sampled when solver type is Random Ray.
    for (int64_t m = 0; m < spatial_bins; ++m) {
      for (int64_t a = 0; a < angle_bins; ++a) {
        spatial_angle_counts(m, a) = 0.0;
      }
    }

    for (int64_t i = 0; i < n_samples_per_source; ++i) {
      SourceSite site = is->sample(&seed);

      int32_t m = spatial_mesh->get_bin(site.r);
      if (m < 0)
        continue;

      int32_t a = angle_mesh ? angle_mesh->get_bin(site.u) : 0;
      if (angle_mesh && a < 0)
        continue;

      spatial_angle_counts(m, a) += 1.0;
    }

    for (int64_t m = 0; m < spatial_bins; ++m) {
      for (int64_t a = 0; a < angle_bins; ++a) {
        double frac = spatial_angle_counts(m, a) / n_samples_per_source;
        if (frac == 0.0)
          continue;
        for (int64_t e = 0; e < energy_bins; ++e) {
          unbiased_strength_(m, a, e) +=
            relative_strength * frac * group_prob[e];
        }
      }
    }
  }
}

void SourceBias::load_forward_source_mesh(const std::string& path)
{
  hid_t file = file_open(path, 'r');

  std::string filetype;
  read_attribute(file, "filetype", filetype);
  if (filetype != "forward_source") {
    file_close(file);
    fatal_error(fmt::format("File '{}' is not a forward source file.", path));
  }

  std::array<int, 2> file_version;
  read_attribute(file, "version", file_version);
  if (file_version[0] != VERSION_SOURCE_BIAS[0]) {
    file_close(file);
    fatal_error(fmt::format("File '{}' has version {} which is incompatible "
                            "with the expected version ({}).",
      path, file_version, VERSION_SOURCE_BIAS));
  }

  // Check that the file's spatial mesh matches the one this SourceBias was
  // built with.
  int32_t spatial_mesh_id;
  read_dataset(file, "spatial_mesh", spatial_mesh_id);
  int32_t expected_spatial_id = model::meshes[spatial_mesh_idx_]->id();
  if (spatial_mesh_id != expected_spatial_id) {
    file_close(file);
    fatal_error(fmt::format(
      "Spatial mesh in '{}' (id {}) does not match the spatial mesh used "
      "for source biasing (id {}).",
      path, spatial_mesh_id, expected_spatial_id));
  }

  int32_t angle_mesh_id;
  read_dataset(file, "angle_mesh", angle_mesh_id);
  int32_t expected_angle_id =
    angle_mesh_idx_ == C_NONE ? C_NONE : model::meshes[angle_mesh_idx_]->id();
  if (angle_mesh_id != expected_angle_id) {
    file_close(file);
    fatal_error(fmt::format(
      "Angular mesh in '{}' (id {}) does not match the angular mesh used "
      "for source biasing (id {}, or none expected).",
      path, angle_mesh_id, expected_angle_id));
  }

  vector<double> file_energy_bounds;
  read_dataset<double>(file, "energy_bounds", file_energy_bounds);
  if (file_energy_bounds.size() != energy_bounds_.size()) {
    file_close(file);
    fatal_error(fmt::format(
      "Energy group structure in '{}' ({} boundaries) does not match the "
      "structure used for source biasing ({} boundaries).",
      path, file_energy_bounds.size(), energy_bounds_.size()));
  }
  for (std::size_t i = 0; i < energy_bounds_.size(); ++i) {
    double tol = 1e-6 * std::max(1.0, std::abs(energy_bounds_[i]));
    if (std::abs(file_energy_bounds[i] - energy_bounds_[i]) > tol) {
      file_close(file);
      fatal_error(fmt::format(
        "Energy group boundary {} in '{}' ({}) does not match the value "
        "used for source biasing ({}).",
        i, path, file_energy_bounds[i], energy_bounds_[i]));
    }
  }

  read_dataset<double>(file, "unbiased_source_strength", unbiased_strength_);

  file_close(file);
}

WeightWindowsGenerator::WeightWindowsGenerator(pugi::xml_node node)
{
  // read information from the XML node
  int32_t mesh_id = std::stoi(get_node_value(node, "mesh"));
  int32_t mesh_idx = model::mesh_map[mesh_id];
  max_realizations_ = std::stoi(get_node_value(node, "max_realizations"));

  int32_t active_batches = settings::n_batches - settings::n_inactive;
  if (max_realizations_ > active_batches) {
    auto msg =
      fmt::format("The maximum number of specified tally realizations ({}) is "
                  "greater than the number of active batches ({}).",
        max_realizations_, active_batches);
    warning(msg);
  }
  auto tmp_str = get_node_value(node, "particle_type", false, true);
  auto particle_type = ParticleType {tmp_str};

  update_interval_ = std::stoi(get_node_value(node, "update_interval"));
  on_the_fly_ = get_node_value_bool(node, "on_the_fly");

  std::vector<double> e_bounds;
  if (check_for_node(node, "energy_bounds")) {
    e_bounds = get_node_array<double>(node, "energy_bounds");
  } else {
    int p_type = particle_type.transport_index();
    if (p_type == C_NONE) {
      fatal_error("Weight windows particle is not supported for transport.");
    }
    e_bounds.push_back(data::energy_min[p_type]);
    e_bounds.push_back(data::energy_max[p_type]);
  }

  // set method
  std::string method_string = get_node_value(node, "method");
  if (method_string == "magic") {
    method_ = WeightWindowUpdateMethod::MAGIC;
    if (settings::solver_type == SolverType::RANDOM_RAY &&
        FlatSourceDomain::adjoint_requested_) {
      fatal_error("Random ray weight window generation with MAGIC cannot be "
                  "done in adjoint mode.");
    }
  } else if (method_string == "fw_cadis") {
    method_ = WeightWindowUpdateMethod::FW_CADIS;
    if (settings::solver_type != SolverType::RANDOM_RAY) {
      fatal_error("FW-CADIS can only be run in random ray solver mode.");
    }
    FlatSourceDomain::adjoint_requested_ = true;
    if (check_for_node(node, "targets")) {
      FlatSourceDomain::fw_cadis_local_ = true;
      targets_ = get_node_array<size_t>(node, "targets");
      FlatSourceDomain::fw_cadis_local_targets_.insert(
        std::end(FlatSourceDomain::fw_cadis_local_targets_),
        std::begin(targets_), std::end(targets_));
    }
    if (check_for_node(node, "source_biasing")) {
      source_biasing_ = get_node_value_bool(node, "source_biasing");
      if (source_biasing_ &&
          check_for_node(node, "angular_biasing_quadrature")) {
        int32_t angle_mesh_id =
          std::stoi(get_node_value(node, "angular_biasing_quadrature"));
        angle_mesh_idx_ = model::mesh_map[angle_mesh_id];
      }
      if (source_biasing_ && check_for_node(node, "forward_source_mesh")) {
        forward_source_mesh_ = get_node_value(node, "forward_source_mesh");
      }
    }
  } else {
    fatal_error(fmt::format(
      "Unknown weight window update method '{}' specified", method_string));
  }

  // parse non-default update parameters if specified
  if (check_for_node(node, "update_parameters")) {
    pugi::xml_node params_node = node.child("update_parameters");
    if (check_for_node(params_node, "value"))
      tally_value_ = get_node_value(params_node, "value");
    if (check_for_node(params_node, "threshold"))
      threshold_ = std::stod(get_node_value(params_node, "threshold"));
    if (check_for_node(params_node, "ratio")) {
      ratio_ = std::stod(get_node_value(params_node, "ratio"));
    }
  }

  // check update parameter values
  if (tally_value_ != "mean" && tally_value_ != "rel_err") {
    fatal_error(fmt::format("Unsupported tally value '{}' specified for "
                            "weight window generation.",
      tally_value_));
  }
  if (threshold_ <= 0.0)
    fatal_error(fmt::format("Invalid relative error threshold '{}' (<= 0.0) "
                            "specified for weight window generation",
      ratio_));
  if (ratio_ <= 1.0)
    fatal_error(fmt::format("Invalid weight window ratio '{}' (<= 1.0) "
                            "specified for weight window generation",
      ratio_));

  // create a matching weight windows object
  auto wws = WeightWindows::create();
  ww_idx_ = wws->index();
  wws->set_mesh(mesh_idx);
  if (e_bounds.size() > 0)
    wws->set_energy_bounds(e_bounds);
  wws->set_particle_type(particle_type);
  wws->set_defaults();
}

void WeightWindowsGenerator::create_tally()
{
  const auto& wws = variance_reduction::weight_windows[ww_idx_];

  // create a tally based on the WWG information
  Tally* ww_tally = Tally::create();
  tally_idx_ = model::tally_map[ww_tally->id()];
  ww_tally->set_scores({"flux"});

  int32_t mesh_id = wws->mesh()->id();
  int32_t mesh_idx = model::mesh_map.at(mesh_id);
  // see if there's already a mesh filter using this mesh
  bool found_mesh_filter = false;
  for (const auto& f : model::tally_filters) {
    if (f->type() == FilterType::MESH) {
      const auto* mesh_filter = dynamic_cast<MeshFilter*>(f.get());
      if (mesh_filter->mesh() == mesh_idx && !mesh_filter->translated() &&
          !mesh_filter->rotated()) {
        ww_tally->add_filter(f.get());
        found_mesh_filter = true;
        break;
      }
    }
  }

  if (!found_mesh_filter) {
    auto mesh_filter = Filter::create("mesh");
    openmc_mesh_filter_set_mesh(mesh_filter->index(), model::mesh_map[mesh_id]);
    ww_tally->add_filter(mesh_filter);
  }

  const auto& e_bounds = wws->energy_bounds();
  if (e_bounds.size() > 0) {
    auto energy_filter = Filter::create("energy");
    openmc_energy_filter_set_bins(
      energy_filter->index(), e_bounds.size(), e_bounds.data());
    ww_tally->add_filter(energy_filter);
  }

  // add a particle filter
  auto particle_type = wws->particle_type();
  auto particle_filter = Filter::create("particle");
  auto pf = dynamic_cast<ParticleFilter*>(particle_filter);
  pf->set_particles({&particle_type, 1});
  ww_tally->add_filter(particle_filter);

  if (!source_biasing_ || method_ != WeightWindowUpdateMethod::FW_CADIS) {
    return;
  }

  // Now add a tally for source biasing
  Tally* sb_tally = Tally::create();
  sb_tally_idx_ = model::tally_map[sb_tally->id()];
  sb_tally->set_scores({"flux"});

  // Add same particle, energy, and spatial mesh filters
  for (int i = 0; i < ww_tally->filters().size(); ++i) {
    sb_tally->add_filter(model::tally_filters[ww_tally->filters(i)].get());
  }

  // Add angular dependency if requested
  if (angle_mesh_idx_ != C_NONE) {
    auto meshangle_filter = Filter::create("meshangular");
    auto maf = dynamic_cast<MeshAngularFilter*>(meshangle_filter);
    maf->set_mesh(angle_mesh_idx_);
    sb_tally->add_filter(meshangle_filter);
  }

  // Create the object that will hold the accumulated flux data used to
  // build a biased forward source once the simulation finishes
  auto sb = std::make_unique<SourceBias>(
    mesh_idx, angle_mesh_idx_, e_bounds, wws->id(), forward_source_mesh_);
  sb_idx_ = static_cast<int32_t>(variance_reduction::source_biases.size());
  variance_reduction::source_biases.push_back(std::move(sb));

  // Pre-compute unbiased source strength distribution
  variance_reduction::source_biases[sb_idx_]->compute_unbiased_strength();
}

void WeightWindowsGenerator::update() const
{
  const auto& wws = variance_reduction::weight_windows[ww_idx_];

  Tally* tally = model::tallies[tally_idx_].get();

  // If in random ray mode, only update on the last batch
  if (settings::solver_type == SolverType::RANDOM_RAY) {
    if (simulation::current_batch != settings::n_batches) {
      return;
    }
    // If in Monte Carlo mode and beyond the number of max realizations or
    // not at the correct update interval, skip the update
  } else if (max_realizations_ < tally->n_realizations_ ||
             tally->n_realizations_ % update_interval_ != 0) {
    return;
  }

  wws->update_weights(tally, tally_value_, threshold_, ratio_, method_);

  // if we're not doing on the fly generation, reset the tally results once
  // we're done with the update
  if (!on_the_fly_)
    tally->reset();

  // Update the source bias flux data the same way, from its own tally
  if (source_biasing_ && method_ == WeightWindowUpdateMethod::FW_CADIS) {
    Tally* sb_tally = model::tallies[sb_tally_idx_].get();
    variance_reduction::source_biases[sb_idx_]->update(sb_tally);

    if (!on_the_fly_)
      sb_tally->reset();
  }

  // TODO: deactivate or remove tally once weight window generation is
  // complete
}

//==============================================================================
// Non-member functions
//==============================================================================

std::pair<bool, WeightWindow> search_weight_window(const Particle& p)
{
  // TODO: this is a linear search - should do something more clever
  for (const auto& ww : variance_reduction::weight_windows) {
    auto [ww_found, weight_window] = ww->get_weight_window(p);
    if (ww_found)
      return {true, weight_window};
  }
  return {false, {}};
}

void apply_weight_windows(Particle& p)
{
  if (!settings::weight_windows_on)
    return;

  // Random ray rays are not Monte Carlo particles and must not be biased by
  // weight windows; the solver generates weight windows but never applies them
  if (settings::solver_type == SolverType::RANDOM_RAY)
    return;

  // WW on photon and neutron only
  if (!p.type().is_neutron() && !p.type().is_photon())
    return;

  // skip dead or no energy
  if (p.E() <= 0 || !p.alive())
    return;

  auto [ww_found, ww] = search_weight_window(p);
  if (ww_found && ww.is_valid()) {
    apply_weight_window(p, ww);
  } else {
    if (p.wgt_ww_born() == -1.0)
      p.wgt_ww_born() = 1.0;
  }
}

void apply_weight_window(Particle& p, WeightWindow weight_window)
{
  if (!weight_window.is_valid())
    return;

  // skip dead or no energy
  if (p.E() <= 0 || !p.alive())
    return;

  // If particle has not yet had its birth weight window value set, set it to
  // the current weight window.
  if (p.wgt_ww_born() == -1.0)
    p.wgt_ww_born() =
      (weight_window.lower_weight + weight_window.upper_weight) / 2;

  // Normalize weight windows based on particle's starting weight
  // and the value of the weight window the particle was born in.
  weight_window.scale(p.wgt_born() / p.wgt_ww_born());

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

  // If the particle's weight is above the weight window, split it until the
  // resulting particles are within the window. The comparisons use a relative
  // dead band so that the branch taken is insensitive to bit-level differences
  // in the window bounds (see WEIGHT_WINDOW_REL_TOL).
  if (weight > weight_window.upper_weight * (1.0 + WEIGHT_WINDOW_REL_TOL)) {
    // do not further split the particle if above the limit
    if (p.n_split() >= settings::max_history_splits)
      return;

    // Dividing by the same dead-banded bound used in the branch condition
    // keeps the number of splits stable when the weight-to-bound ratio sits
    // within rounding of an exact integer, which the weight window arithmetic
    // itself can produce (e.g., a roulette survivor assigned weight *
    // max_split, later split against an upper bound that is an exact multiple
    // of the same lower bound). Ratios within the dead band of an integer
    // consistently round down, and the branch condition guarantees the ratio
    // exceeds one; the lower clamp of 2 makes the always-splits invariant
    // explicit.
    double n_split = std::max(2.0,
      std::ceil(
        weight / ((1.0 + WEIGHT_WINDOW_REL_TOL) * weight_window.upper_weight)));
    double max_split = weight_window.max_split;
    n_split = std::min(n_split, max_split);

    p.n_split() += n_split;

    // Create secondaries and divide weight among all particles
    int i_split = std::round(n_split);
    for (int l = 0; l < i_split - 1; l++) {
      p.split(weight / n_split);
    }
    // remaining weight is applied to current particle
    p.wgt() = weight / n_split;

  } else if (weight <
             weight_window.lower_weight * (1.0 - WEIGHT_WINDOW_REL_TOL)) {
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

void finalize_variance_reduction()
{
  for (const auto& wwg : variance_reduction::weight_windows_generators) {
    wwg->create_tally();
  }
}

//==============================================================================
// C API
//==============================================================================

int verify_ww_index(int32_t index)
{
  if (index < 0 || index >= variance_reduction::weight_windows.size()) {
    set_errmsg(fmt::format("Index '{}' for weight windows is invalid", index));
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

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

extern "C" int openmc_weight_windows_get_id(int32_t index, int32_t* id)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows.at(index);
  *id = wws->id();
  return 0;
}

extern "C" int openmc_weight_windows_set_id(int32_t index, int32_t id)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows.at(index);
  wws->set_id(id);
  return 0;
}

extern "C" int openmc_weight_windows_update_magic(int32_t ww_idx,
  int32_t tally_idx, const char* value, double threshold, double ratio)
{
  if (int err = verify_ww_index(ww_idx))
    return err;

  if (tally_idx < 0 || tally_idx >= model::tallies.size()) {
    set_errmsg(fmt::format("Index '{}' for tally is invalid", tally_idx));
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  // get the requested tally
  const Tally* tally = model::tallies.at(tally_idx).get();

  // get the WeightWindows object
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);

  wws->update_weights(tally, value, threshold, ratio);

  return 0;
}

extern "C" int openmc_weight_windows_set_mesh(int32_t ww_idx, int32_t mesh_idx)
{
  if (int err = verify_ww_index(ww_idx))
    return err;
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);
  wws->set_mesh(mesh_idx);
  return 0;
}

extern "C" int openmc_weight_windows_get_mesh(int32_t ww_idx, int32_t* mesh_idx)
{
  if (int err = verify_ww_index(ww_idx))
    return err;
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);
  *mesh_idx = model::mesh_map.at(wws->mesh()->id());
  return 0;
}

extern "C" int openmc_weight_windows_set_energy_bounds(
  int32_t ww_idx, double* e_bounds, size_t e_bounds_size)
{
  if (int err = verify_ww_index(ww_idx))
    return err;
  const auto& wws = variance_reduction::weight_windows.at(ww_idx);
  wws->set_energy_bounds({e_bounds, e_bounds_size});
  return 0;
}

extern "C" int openmc_weight_windows_get_energy_bounds(
  int32_t ww_idx, const double** e_bounds, size_t* e_bounds_size)
{
  if (int err = verify_ww_index(ww_idx))
    return err;
  const auto& wws = variance_reduction::weight_windows[ww_idx].get();
  *e_bounds = wws->energy_bounds().data();
  *e_bounds_size = wws->energy_bounds().size();
  return 0;
}

extern "C" int openmc_weight_windows_set_particle(
  int32_t index, int32_t particle)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows.at(index);
  wws->set_particle_type(ParticleType {particle});
  return 0;
}

extern "C" int openmc_weight_windows_get_particle(
  int32_t index, int32_t* particle)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows.at(index);
  *particle = wws->particle_type().pdg_number();
  return 0;
}

extern "C" int openmc_weight_windows_get_bounds(int32_t index,
  const double** lower_bounds, const double** upper_bounds, size_t* size)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows[index];
  *size = wws->lower_ww_bounds().size();
  *lower_bounds = wws->lower_ww_bounds().data();
  *upper_bounds = wws->upper_ww_bounds().data();
  return 0;
}

extern "C" int openmc_weight_windows_set_bounds(int32_t index,
  const double* lower_bounds, const double* upper_bounds, size_t size)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows[index];
  wws->set_bounds(span<const double>(lower_bounds, size),
    span<const double>(upper_bounds, size));
  return 0;
}

extern "C" int openmc_weight_windows_get_survival_ratio(
  int32_t index, double* ratio)
{
  if (int err = verify_ww_index(index))
    return err;
  const auto& wws = variance_reduction::weight_windows[index];
  *ratio = wws->survival_ratio();
  return 0;
}

extern "C" int openmc_weight_windows_set_survival_ratio(
  int32_t index, double ratio)
{
  if (int err = verify_ww_index(index))
    return err;
  const auto& wws = variance_reduction::weight_windows[index];
  wws->survival_ratio() = ratio;
  std::cout << "Survival ratio: " << wws->survival_ratio() << std::endl;
  return 0;
}

extern "C" int openmc_weight_windows_get_max_lower_bound_ratio(
  int32_t index, double* lb_ratio)
{
  if (int err = verify_ww_index(index))
    return err;
  const auto& wws = variance_reduction::weight_windows[index];
  *lb_ratio = wws->max_lower_bound_ratio();
  return 0;
}

extern "C" int openmc_weight_windows_set_max_lower_bound_ratio(
  int32_t index, double lb_ratio)
{
  if (int err = verify_ww_index(index))
    return err;
  const auto& wws = variance_reduction::weight_windows[index];
  wws->max_lower_bound_ratio() = lb_ratio;
  return 0;
}

extern "C" int openmc_weight_windows_get_weight_cutoff(
  int32_t index, double* cutoff)
{
  if (int err = verify_ww_index(index))
    return err;
  const auto& wws = variance_reduction::weight_windows[index];
  *cutoff = wws->weight_cutoff();
  return 0;
}

extern "C" int openmc_weight_windows_set_weight_cutoff(
  int32_t index, double cutoff)
{
  if (int err = verify_ww_index(index))
    return err;
  const auto& wws = variance_reduction::weight_windows[index];
  wws->weight_cutoff() = cutoff;
  return 0;
}

extern "C" int openmc_weight_windows_get_max_split(
  int32_t index, int* max_split)
{
  if (int err = verify_ww_index(index))
    return err;
  const auto& wws = variance_reduction::weight_windows[index];
  *max_split = wws->max_split();
  return 0;
}

extern "C" int openmc_weight_windows_set_max_split(int32_t index, int max_split)
{
  if (int err = verify_ww_index(index))
    return err;
  const auto& wws = variance_reduction::weight_windows[index];
  wws->max_split() = max_split;
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
    variance_reduction::weight_windows.push_back(make_unique<WeightWindows>());
  return 0;
}

extern "C" size_t openmc_weight_windows_size()
{
  return variance_reduction::weight_windows.size();
}

extern "C" int openmc_weight_windows_export(const char* filename)
{

  if (!mpi::master)
    return 0;

  std::string name = filename ? filename : "weight_windows.h5";

  write_message(fmt::format("Exporting weight windows to {}...", name), 5);

  hid_t ww_file = file_open(name, 'w');

  // Write file type
  write_attribute(ww_file, "filetype", "weight_windows");

  // Write revisiion number for state point file
  write_attribute(ww_file, "version", VERSION_WEIGHT_WINDOWS);

  hid_t weight_windows_group = create_group(ww_file, "weight_windows");

  hid_t mesh_group = create_group(ww_file, "meshes");

  std::vector<int32_t> mesh_ids;
  std::vector<int32_t> ww_ids;
  for (const auto& ww : variance_reduction::weight_windows) {

    ww->to_hdf5(weight_windows_group);
    ww_ids.push_back(ww->id());

    // if the mesh has already been written, move on
    int32_t mesh_id = ww->mesh()->id();
    if (std::find(mesh_ids.begin(), mesh_ids.end(), mesh_id) != mesh_ids.end())
      continue;

    mesh_ids.push_back(mesh_id);
    ww->mesh()->to_hdf5(mesh_group);
  }

  write_attribute(mesh_group, "n_meshes", mesh_ids.size());
  write_attribute(mesh_group, "ids", mesh_ids);
  close_group(mesh_group);

  write_attribute(weight_windows_group, "n_weight_windows", ww_ids.size());
  write_attribute(weight_windows_group, "ids", ww_ids);
  close_group(weight_windows_group);

  file_close(ww_file);

  // If any source biasing data has been generated, export it too
  if (!variance_reduction::source_biases.empty()) {
    std::string sb_name = "source_bias.h5";

    write_message(
      fmt::format("Exporting source bias data to {}...", sb_name), 5);

    hid_t sb_file = file_open(sb_name, 'w');

    write_attribute(sb_file, "filetype", "source_bias");
    write_attribute(sb_file, "version", VERSION_SOURCE_BIAS);

    hid_t sb_mesh_group = create_group(sb_file, "meshes");
    std::vector<int32_t> sb_mesh_ids;

    auto write_mesh_once = [&](int32_t mesh_idx) {
      if (mesh_idx == C_NONE)
        return;
      int32_t mesh_id = model::meshes[mesh_idx]->id();
      if (std::find(sb_mesh_ids.begin(), sb_mesh_ids.end(), mesh_id) !=
          sb_mesh_ids.end())
        return;
      sb_mesh_ids.push_back(mesh_id);
      model::meshes[mesh_idx]->to_hdf5(sb_mesh_group);
    };

    for (const auto& sb : variance_reduction::source_biases) {
      sb->to_hdf5(sb_file);
      write_mesh_once(sb->spatial_mesh_idx());
      write_mesh_once(sb->angle_mesh_idx());
    }

    write_attribute(sb_mesh_group, "n_meshes", sb_mesh_ids.size());
    write_attribute(sb_mesh_group, "ids", sb_mesh_ids);
    close_group(sb_mesh_group);

    file_close(sb_file);
  }

  return 0;
}

extern "C" int openmc_weight_windows_import(const char* filename)
{
  std::string name = filename ? filename : "weight_windows.h5";

  if (mpi::master)
    write_message(fmt::format("Importing weight windows from {}...", name), 5);

  if (!file_exists(name)) {
    set_errmsg(fmt::format("File '{}' does not exist", name));
  }

  hid_t ww_file = file_open(name, 'r');

  // Check that filetype is correct
  std::string filetype;
  read_attribute(ww_file, "filetype", filetype);
  if (filetype != "weight_windows") {
    file_close(ww_file);
    set_errmsg(fmt::format("File '{}' is not a weight windows file.", name));
    return OPENMC_E_INVALID_ARGUMENT;
  }

  // Check that the file version is compatible
  std::array<int, 2> file_version;
  read_attribute(ww_file, "version", file_version);
  if (file_version[0] != VERSION_WEIGHT_WINDOWS[0]) {
    std::string err_msg =
      fmt::format("File '{}' has version {} which is incompatible with the "
                  "expected version ({}).",
        name, file_version, VERSION_WEIGHT_WINDOWS);
    set_errmsg(err_msg);
    return OPENMC_E_INVALID_ARGUMENT;
  }

  hid_t weight_windows_group = open_group(ww_file, "weight_windows");

  hid_t mesh_group = open_group(ww_file, "meshes");

  read_meshes(mesh_group);

  std::vector<std::string> names = group_names(weight_windows_group);

  for (const auto& name : names) {
    WeightWindows::from_hdf5(weight_windows_group, name);
  }

  close_group(weight_windows_group);

  file_close(ww_file);

  return 0;
}

} // namespace openmc
