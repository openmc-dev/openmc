#include "openmc/weight_windows.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <set>
#include <string>

#include "xtensor/xindex_view.hpp"
#include "xtensor/xio.hpp"
#include "xtensor/xmasked_view.hpp"
#include "xtensor/xnoalias.hpp"
#include "xtensor/xstrided_view.hpp"
#include "xtensor/xview.hpp"

#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/physics_common.h"
#include "openmc/random_ray/flat_source_domain.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_mesh.h"
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

} // namespace variance_reduction

//==============================================================================
// Non-member functions
//==============================================================================

void apply_weight_windows(Particle& p)
{
  if (!settings::weight_windows_on)
    return;

  // WW on photon and neutron only
  if (p.type() != ParticleType::neutron && p.type() != ParticleType::photon)
    return;

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

  // If particle has not yet had its birth weight window value set, set it to
  // the current weight window (or 1.0 if not born in a weight window).
  if (p.wgt_ww_born() == -1.0) {
    if (weight_window.is_valid()) {
      p.wgt_ww_born() =
        (weight_window.lower_weight + weight_window.upper_weight) / 2;
    } else {
      p.wgt_ww_born() = 1.0;
    }
  }

  // particle is not in any of the ww domains, do nothing
  if (!weight_window.is_valid())
    return;

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

  // if particle's weight is above the weight window split until they are within
  // the window
  if (weight > weight_window.upper_weight) {
    // do not further split the particle if above the limit
    if (p.n_split() >= settings::max_history_splits)
      return;

    double n_split = std::ceil(weight / weight_window.upper_weight);
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
  particle_type_ = openmc::str_to_particle_type(particle_type_str);

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
  wws->particle_type_ = openmc::str_to_particle_type(particle_type);

  read_dataset<double>(ww_group, "energy_bounds", wws->energy_bounds_);

  int32_t mesh_id;
  read_dataset(ww_group, "mesh", mesh_id);

  if (model::mesh_map.count(mesh_id) == 0) {
    fatal_error(
      fmt::format("Mesh {} used in weight windows does not exist.", mesh_id));
  }
  wws->set_mesh(model::mesh_map[mesh_id]);

  wws->lower_ww_ = xt::empty<double>(wws->bounds_size());
  wws->upper_ww_ = xt::empty<double>(wws->bounds_size());

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
    int p_type = static_cast<int>(particle_type_);
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
  lower_ww_ = xt::empty<double>(shape);
  lower_ww_.fill(-1);
  upper_ww_ = xt::empty<double>(shape);
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
  if (p_type != ParticleType::neutron && p_type != ParticleType::photon)
    fatal_error(
      fmt::format("Particle type '{}' cannot be applied to weight windows.",
        particle_type_to_str(p_type)));
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

WeightWindow WeightWindows::get_weight_window(const Particle& p) const
{
  // check for particle type
  if (particle_type_ != p.type()) {
    return {};
  }

  // Get mesh index for particle's position
  const auto& mesh = this->mesh();
  int mesh_bin = mesh->get_bin(p.r());

  // particle is outside the weight window mesh
  if (mesh_bin < 0)
    return {};

  // particle energy
  double E = p.E();

  // check to make sure energy is in range, expects sorted energy values
  if (E < energy_bounds_.front() || E > energy_bounds_.back())
    return {};

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
  return ww;
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

void WeightWindows::set_bounds(const xt::xtensor<double, 2>& lower_bounds,
  const xt::xtensor<double, 2>& upper_bounds)
{

  this->check_bounds(lower_bounds, upper_bounds);

  // set new weight window values
  lower_ww_ = lower_bounds;
  upper_ww_ = upper_bounds;
}

void WeightWindows::set_bounds(
  const xt::xtensor<double, 2>& lower_bounds, double ratio)
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
  lower_ww_ = xt::empty<double>(shape);
  upper_ww_ = xt::empty<double>(shape);

  // set new weight window values
  xt::view(lower_ww_, xt::all()) =
    xt::adapt(lower_bounds.data(), lower_ww_.shape());
  xt::view(upper_ww_, xt::all()) =
    xt::adapt(upper_bounds.data(), upper_ww_.shape());
}

void WeightWindows::set_bounds(span<const double> lower_bounds, double ratio)
{
  this->check_bounds(lower_bounds);

  auto shape = this->bounds_size();
  lower_ww_ = xt::empty<double>(shape);
  upper_ww_ = xt::empty<double>(shape);

  // set new weight window values
  xt::view(lower_ww_, xt::all()) =
    xt::adapt(lower_bounds.data(), lower_ww_.shape());
  xt::view(upper_ww_, xt::all()) =
    xt::adapt(lower_bounds.data(), upper_ww_.shape());
  upper_ww_ *= ratio;
}

void WeightWindows::update_weights(const Tally* tally, const std::string& value,
  double threshold, double ratio, WeightWindowUpdateMethod method)
{
  ///////////////////////////
  // Setup and checks
  ///////////////////////////
  this->check_tally_update_compatibility(tally);

  lower_ww_.fill(-1);
  upper_ww_.fill(-1);

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
  // At the end of this section, the mean and rel_err array
  // is a 2D view of tally data (n_e_groups, n_mesh_bins)
  //
  ///////////////////////////

  // build a shape for a view of the tally results, this will always be
  // dimension 5 (3 filter dimensions, 1 score dimension, 1 results dimension)
  std::array<int, 5> shape = {
    1, 1, 1, tally->n_scores(), static_cast<int>(TallyResult::SIZE)};

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

  // get a fully reshaped view of the tally according to tally ordering of
  // filters
  auto tally_values = xt::reshape_view(tally->results(), shape);

  // get a that is (particle, energy, mesh, scores, values)
  auto transposed_view = xt::transpose(tally_values, transpose);

  // determine the dimension and index of the particle data
  int particle_idx = 0;
  if (tally->has_filter(FilterType::PARTICLE)) {
    // get the particle filter
    auto pf = tally->get_filter<ParticleFilter>();
    const auto& particles = pf->particles();

    // find the index of the particle that matches these weight windows
    auto p_it =
      std::find(particles.begin(), particles.end(), this->particle_type_);
    // if the particle filter doesn't have particle data for the particle
    // used on this weight windows instance, report an error
    if (p_it == particles.end()) {
      auto msg = fmt::format("Particle type '{}' not present on Filter {} for "
                             "Tally {} used to update WeightWindows {}",
        particle_type_to_str(this->particle_type_), pf->id(), tally->id(),
        this->id());
      fatal_error(msg);
    }

    // use the index of the particle in the filter to down-select data later
    particle_idx = p_it - particles.begin();
  }

  // down-select data based on particle and score
  auto sum = xt::view(transposed_view, particle_idx, xt::all(), xt::all(),
    score_index, static_cast<int>(TallyResult::SUM));
  auto sum_sq = xt::view(transposed_view, particle_idx, xt::all(), xt::all(),
    score_index, static_cast<int>(TallyResult::SUM_SQ));
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

  // get mesh volumes and dimensions
  auto mesh_vols = this->mesh()->volumes();
  int e_bins = new_bounds.shape()[0];
  int mesh_bins = new_bounds.shape()[1];

  // Parallelize ALL calculations - eliminate ALL xtensor element-wise operations
  // First evaluate the xtensor views into concrete arrays for parallel access
  auto sum_evaluated = xt::eval(sum);
  auto sum_sq_evaluated = xt::eval(sum_sq);
  
  // Calculate mean (new_bounds) and relative error in parallel
  // Also handle rel_err inversion in the same loop to avoid separate parallel region
#pragma omp parallel for collapse(2) schedule(runtime)
  for (int e = 0; e < e_bins; e++) {
    for (int m = 0; m < mesh_bins; m++) {
      // Calculate mean
      new_bounds(e, m) = sum_evaluated(e, m) / n;
      
      // Calculate relative error
      if (sum_evaluated(e, m) > 0.0) {
        double mean_val = new_bounds(e, m);
        double variance = (sum_sq_evaluated(e, m) / n - mean_val * mean_val) / (n - 1);
        rel_err(e, m) = std::sqrt(variance) / mean_val;
      } else {
        rel_err(e, m) = INFTY;
      }
      
      // Handle rel_err inversion in the same loop if needed
      if (value == "rel_err") {
        new_bounds(e, m) = 1.0 / rel_err(e, m);
      }
    }
  }

  if (method == WeightWindowUpdateMethod::MAGIC) {
    // If we are computing weight windows with forward fluxes derived from a
    // Monte Carlo or forward random ray solve, we use the MAGIC algorithm.
    
    // First pass: divide by volume of mesh elements in parallel
#pragma omp parallel for collapse(2) schedule(runtime)
    for (int e = 0; e < e_bins; e++) {
      for (int i = 0; i < mesh_bins; i++) {
        new_bounds(e, i) /= mesh_vols[i];
      }
    }

    // Second pass: find group maximum and normalize (per energy group)
    for (int e = 0; e < e_bins; e++) {
      double group_max = 0.0;
      
      // Find maximum in this energy group using parallel reduction
#pragma omp parallel for schedule(runtime) reduction(max:group_max)
      for (int i = 0; i < mesh_bins; i++) {
        if (new_bounds(e, i) > group_max) {
          group_max = new_bounds(e, i);
        }
      }
      
      // Normalize values in this energy group by the maximum value
      if (group_max > 0.0) {
        double norm_factor = 2.0 * group_max;
#pragma omp parallel for schedule(runtime)
        for (int i = 0; i < mesh_bins; i++) {
          new_bounds(e, i) /= norm_factor;
        }
      }
    }
  } else {
    // If we are computing weight windows with adjoint fluxes derived from an
    // adjoint random ray solve, we use the FW-CADIS algorithm.
    
    // Combine volume division and inverse calculation into one parallel loop
#pragma omp parallel for collapse(2) schedule(runtime)
    for (int e = 0; e < e_bins; e++) {
      for (int i = 0; i < mesh_bins; i++) {
        // Divide by volume of mesh elements
        new_bounds(e, i) /= mesh_vols[i];
        
        // Take the inverse, but are careful not to divide by zero
        if (new_bounds(e, i) != 0.0) {
          new_bounds(e, i) = 1.0 / new_bounds(e, i);
        } else {
          new_bounds(e, i) = 0.0;
        }
      }
    }
    
    // Find the maximum value across all elements using parallel reduction
    double max_val = 0.0;
#pragma omp parallel for collapse(2) schedule(runtime) reduction(max:max_val)
    for (int e = 0; e < e_bins; e++) {
      for (int i = 0; i < mesh_bins; i++) {
        if (new_bounds(e, i) > max_val) {
          max_val = new_bounds(e, i);
        }
      }
    }
    
    // Combine normalization and minimum finding for efficiency
    if (max_val > 0.0) {
      double norm_factor = 2.0 * max_val;
      double min_val = INFTY;
      
      // Normalize and find minimum in one pass
#pragma omp parallel for collapse(2) schedule(runtime) reduction(min:min_val)
      for (int e = 0; e < e_bins; e++) {
        for (int i = 0; i < mesh_bins; i++) {
          new_bounds(e, i) /= norm_factor;
          if (new_bounds(e, i) != 0.0 && new_bounds(e, i) < min_val) {
            min_val = new_bounds(e, i);
          }
        }
      }
      
      // Apply minimum value to missed bins
#pragma omp parallel for collapse(2) schedule(runtime)
      for (int e = 0; e < e_bins; e++) {
        for (int i = 0; i < mesh_bins; i++) {
          if (new_bounds(e, i) == 0.0) {
            new_bounds(e, i) = min_val;
          }
        }
      }
    }
  }

  // Combine all final processing into a single parallel loop to reduce overhead
#pragma omp parallel for collapse(2) schedule(runtime)
  for (int e = 0; e < e_bins; e++) {
    for (int i = 0; i < mesh_bins; i++) {
      // Check if values where the mean is zero should be ignored
      if (sum_evaluated(e, i) <= 0.0) {
        new_bounds(e, i) = -1.0;
      }
      // Check if relative error is higher than threshold
      else if (rel_err(e, i) > threshold) {
        new_bounds(e, i) = -1.0;
      }
      
      // Update the upper bounds (always do this regardless of bounds value)
      upper_ww_(e, i) = ratio * lower_ww_(e, i);
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
  write_dataset(
    ww_group, "particle_type", openmc::particle_type_to_str(particle_type_));
  write_dataset(ww_group, "energy_bounds", energy_bounds_);
  write_dataset(ww_group, "lower_ww_bounds", lower_ww_);
  write_dataset(ww_group, "upper_ww_bounds", upper_ww_);
  write_dataset(ww_group, "survival_ratio", survival_ratio_);
  write_dataset(ww_group, "max_lower_bound_ratio", max_lb_ratio_);
  write_dataset(ww_group, "max_split", max_split_);
  write_dataset(ww_group, "weight_cutoff", weight_cutoff_);

  close_group(ww_group);
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
  auto tmp_str = get_node_value(node, "particle_type", true, true);
  auto particle_type = str_to_particle_type(tmp_str);

  update_interval_ = std::stoi(get_node_value(node, "update_interval"));
  on_the_fly_ = get_node_value_bool(node, "on_the_fly");

  std::vector<double> e_bounds;
  if (check_for_node(node, "energy_bounds")) {
    e_bounds = get_node_array<double>(node, "energy_bounds");
  } else {
    int p_type = static_cast<int>(particle_type);
    e_bounds.push_back(data::energy_min[p_type]);
    e_bounds.push_back(data::energy_max[p_type]);
  }

  // set method
  std::string method_string = get_node_value(node, "method");
  if (method_string == "magic") {
    method_ = WeightWindowUpdateMethod::MAGIC;
    if (settings::solver_type == SolverType::RANDOM_RAY &&
        FlatSourceDomain::adjoint_) {
      fatal_error("Random ray weight window generation with MAGIC cannot be "
                  "done in adjoint mode.");
    }
  } else if (method_string == "fw_cadis") {
    method_ = WeightWindowUpdateMethod::FW_CADIS;
    if (settings::solver_type != SolverType::RANDOM_RAY) {
      fatal_error("FW-CADIS can only be run in random ray solver mode.");
    }
    FlatSourceDomain::adjoint_ = true;
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
                            "specified for weight window generation"));

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
      if (mesh_filter->mesh() == mesh_idx && !mesh_filter->translated()) {
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
}

void WeightWindowsGenerator::update() const
{
  const auto& wws = variance_reduction::weight_windows[ww_idx_];

  Tally* tally = model::tallies[tally_idx_].get();

  // if we're beyond the number of max realizations or not at the corrrect
  // update interval, skip the update
  if (max_realizations_ < tally->n_realizations_ ||
      tally->n_realizations_ % update_interval_ != 0)
    return;

  wws->update_weights(tally, tally_value_, threshold_, ratio_, method_);

  // if we're not doing on the fly generation, reset the tally results once
  // we're done with the update
  if (!on_the_fly_)
    tally->reset();

  // TODO: deactivate or remove tally once weight window generation is
  // complete
}

//==============================================================================
// Non-member functions
//==============================================================================

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

extern "C" int openmc_weight_windows_set_particle(int32_t index, int particle)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows.at(index);
  wws->set_particle_type(static_cast<ParticleType>(particle));
  return 0;
}

extern "C" int openmc_weight_windows_get_particle(int32_t index, int* particle)
{
  if (int err = verify_ww_index(index))
    return err;

  const auto& wws = variance_reduction::weight_windows.at(index);
  *particle = static_cast<int>(wws->particle_type());
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
  wws->set_bounds({lower_bounds, size}, {upper_bounds, size});
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

  std::vector<std::string> names = group_names(weight_windows_group);

  for (const auto& name : names) {
    WeightWindows::from_hdf5(weight_windows_group, name);
  }

  close_group(weight_windows_group);

  file_close(ww_file);

  return 0;
}

} // namespace openmc
