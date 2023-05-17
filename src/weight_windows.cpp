#include "openmc/weight_windows.h"

#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/physics_common.h"
#include "openmc/search.h"
#include "openmc/xml_interface.h"

#include <fmt/core.h>
#include <gsl/gsl-lite.hpp>
// new 
#include "openmc/tallies/tally.h"
#include "openmc/math_functions.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace variance_reduction {

std::unordered_map<int32_t, int32_t> ww_map;
openmc::vector<unique_ptr<WeightWindows>> weight_windows;
  
// new in GVR 
bool global_on {false};
vector<double> iteration;
int32_t tally_idx {-1};
int source_space {1};
int n_statistics {10000};
int64_t nps_backup {-1};

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
  using namespace pugi;
  // check if GVR is used 
  if (check_for_node(node, "global")) {
    // turn on the flag
    variance_reduction::global_on = true;
    std::cout<<"Using GVR"<<std::endl;

    xml_node node_global = node.child("global");

    // Make sure required elements are present
    const vector<std::string> required_elems {"iteration", "tally",
      "source_space"};
    for (const auto& elem : required_elems) {
      if (!check_for_node(node_global, elem.c_str())) {
        fatal_error(fmt::format("Must specify <{}> for GVR.", elem));
      }
    }

    // iteration during WW generation
    variance_reduction::iteration = get_node_array<double>(node_global, "iteration"); 
    for (int i = 0; i<variance_reduction::iteration.size(); i++) std::cout<<"  "<<variance_reduction::iteration[i]<<"  ";
      std::cout<<"  "<<std::endl; 

    // tally for WW generation
    variance_reduction::tally_idx = std::stoi(get_node_value(node_global, "tally"));
      std::cout<<"  tally "<<variance_reduction::tally_idx<<"  "<<std::endl;

    // S value in the PS-GVR method
    variance_reduction::source_space = std::stod(get_node_value(node_global, "source_space"));
      std::cout<<"  S "<<variance_reduction::source_space<<"  "<<std::endl;     

    // n_statistics value in the PS-GVR method
    if (check_for_node(node_global, "n_statistics")) {
      variance_reduction::n_statistics = std::stod(get_node_value(node_global, "n_statistics"));
      std::cout<<"  n_statistics "<<variance_reduction::n_statistics<<"  "<<std::endl;
    } else std::cout<<"  default n_statistics "<<variance_reduction::n_statistics<<"  "<<std::endl;

  }
  
  if (variance_reduction::global_on) {
    // Make sure required elements are present for GVR or LVR
    const vector<std::string> required_elems {"id", "particle_type",
      "energy_bounds"};
    for (const auto& elem : required_elems) {
      if (!check_for_node(node, elem.c_str())) {
        fatal_error(fmt::format("Must specify <{}> for weight windows.", elem));
      }
    }
  } else {
    // Make sure required elements are present for WWM
    const vector<std::string> required_elems {"id", "particle_type",
      "energy_bounds", "lower_ww_bounds", "upper_ww_bounds"};
    for (const auto& elem : required_elems) {
      if (!check_for_node(node, elem.c_str())) {
        fatal_error(fmt::format("Must specify <{}> for weight windows.", elem));
      }
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

  if (variance_reduction::global_on) {
    // set the lower/upper weight bounds for GVR
    int mesh_size = this->mesh().n_bins();
    int energy_size = energy_bounds_.size()-1;
    for (int i = 0; i<mesh_size*energy_size; i++) {
      lower_ww_.push_back(-1);
      upper_ww_.push_back(-5);
    }
  } else {
    // read the lower/upper weight bounds for WWM
    lower_ww_ = get_node_array<double>(node, "lower_ww_bounds");
    upper_ww_ = get_node_array<double>(node, "upper_ww_bounds");
  }

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

  // make sure that the upper and lower bounds have the same size
  if (upper_ww_.size() != lower_ww_.size()) {
    fatal_error("The upper and lower weight window lengths do not match.");
  }

  // num spatial*energy bins must match num weight bins
  int num_spatial_bins = this->mesh().n_bins();
  int num_energy_bins = energy_bounds_.size() - 1;
  int num_weight_bins = lower_ww_.size();
  if (num_weight_bins != num_spatial_bins * num_energy_bins) {
    auto err_msg =
      fmt::format("In weight window domain {} the number of "
                  "energy/spatial bins ({}) does not match the number "
                  "of weight bins provided ({})",
        id_, num_energy_bins * num_spatial_bins, num_weight_bins);
    fatal_error(err_msg);
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
  
// new in GVR
//==============================================================================

std::pair<double, double> mean_stdev_GVR(const double* x, int n)
{
  double mean = x[static_cast<int>(TallyResult::SUM)] / n;
  double stdev =
    n > 1 ? std::sqrt(std::max(0.0,
              (x[static_cast<int>(TallyResult::SUM_SQ)] / n - mean * mean) /
                (n - 1)))
          : 0.0;
  return {mean, stdev};
}
  
void WeightWindows::calculate_WW() 
{
  // parameters for WW calculation
  double max_flux = -1;
  double min_flux = INFTY; 
  int mesh_size = this->mesh().n_bins();
  int energy_size = energy_bounds_.size()-1;
  std::cout<<"mesh_size "<<mesh_size<<", energy_size "<<energy_size<<std::endl;
  vector<double> flux_data(mesh_size*energy_size, 0);
  vector<double> max_flux_data(energy_size, 0);
  vector<double> min_flux_data(energy_size, 0);
  bool tally_found = false;

  // get the neutron flux from tally
  for (auto i_tally = 0; i_tally < model::tallies.size(); ++i_tally) {
    const auto& tally {*model::tallies[i_tally]}; 

    if (variance_reduction::tally_idx == tally.id_) {
      tally_found = true; 
  
      // Calculate t-value for confidence intervals
      double t_value = 1;
      if (openmc::settings::confidence_intervals) {
        auto alpha = 1 - CONFIDENCE_LEVEL;
        t_value = t_percentile(1 - alpha*0.5, tally.n_realizations_ - 1);
      }
  
      // find the max & min flux for each energy group
      for (int ee = 0; ee < energy_size; ++ee) {
        // celar all data before each energy group
        max_flux = -1;
        min_flux = INFTY;
        // mesh
        for (int mm = 0; mm < mesh_size; ++mm) {    
          double mean, stdev;
          std::tie(mean, stdev) = 
            mean_stdev_GVR(&tally.results_(mm+mesh_size*ee,0,0),
              tally.n_realizations_);
          if (mean > 0) {
            flux_data[mm+ee*mesh_size] = mean;
            // find the max flux
            if (mean > max_flux) max_flux = mean;
            // find the min flux
            if (mean < min_flux) min_flux = mean;
          }
        }
    
        // save the max & min flux for each energy group
        max_flux_data[ee] = max_flux;
        min_flux_data[ee] = min_flux;
        std::cout<<"max_flux = "<<max_flux<<", min_flux = "<<min_flux<<" "<<std::endl;
      }
    }
    if (tally_found) break;
  }
    if (!tally_found) fatal_error("The tally for GVR is not found. ");

  // calculate the WW based on PS-GVR
  double near_source = double(variance_reduction::n_statistics*variance_reduction::source_space)/double(variance_reduction::nps_backup);
  if (near_source > 0.2) near_source = 0.2;
  std::cout<<"near source "<<near_source<<std::endl;

  // for each energy group
  for (int ee = 0; ee < energy_size; ++ee) {
    // celar all data before each energy group
    double PS_k = 0.4/log(1/near_source);
    double PS_b = 0.5 - log(max_flux_data[ee]/min_flux_data[ee])*PS_k;

    // mesh
    for (int mm = 0; mm < mesh_size; ++mm) {    
      lower_ww_[mm+ee*mesh_size] = -1;
      upper_ww_[mm+ee*mesh_size] = -5;
      if (flux_data[mm+ee*mesh_size] <= 0) continue;
      if (flux_data[mm+ee*mesh_size] >= near_source*max_flux_data[ee]) lower_ww_[mm+ee*mesh_size] = log(flux_data[mm+ee*mesh_size]/min_flux_data[ee]) * PS_k + PS_b;
      else lower_ww_[mm+ee*mesh_size] = 0.1*flux_data[mm+ee*mesh_size]/(near_source*max_flux_data[ee]);
      upper_ww_[mm+ee*mesh_size] = 5*lower_ww_[mm+ee*mesh_size];
      //if ( std::abs(mm-34460) <= 20) std::cout<<"mm, "<<mm<<", ee "<<ee<<", flux "<<flux_data[mm+ee*mesh_size]<<", lower "<<lower_ww_[mm+ee*mesh_size]<<", upper "<<upper_ww_[mm+ee*mesh_size]<<std::endl;
    }
  }  

}

} // namespace openmc
