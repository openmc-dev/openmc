#include "openmc/error.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/search.h"
#include "openmc/weight_windows.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// WeightWindowMesh implementation
//==============================================================================

WeightWindowMesh::WeightWindowMesh() 
{
  // default parameters
  n_ww = false;                  // flag for neutron use weight window
  p_ww = false;                  // flag for photon use weight window
  user_defined_biasing = false;  // flag for weight biasing in energy
  
  // WWP
  // neutron
  n_upper_ratio    = 5.;         // upper weight window = upper_ratio * lower weight window
  n_survival_ratio = 3.;         // survival weight = survival_ratio * lower weight window
  n_max_split      = 5;          // max number of split particles
  n_multiplier     = 1.;         // multiplier for weight window lower bounds

  // photon
  p_upper_ratio    = 5.;         // upper weight window = upper_ratio * lower weight window
  p_survival_ratio = 3.;         // survival weight = survival_ratio * lower weight window
  p_max_split      = 5;          // max number of split particles
  p_multiplier     = 1.;         // multiplier for weight window lower bounds
}

WeightWindowMesh::WeightWindowMesh(pugi::xml_node node)
{
  using namespace pugi;

  WeightWindowMesh(); // call the default constructor
   
  // construct the mesh
  mesh_ = make_unique<RectilinearMesh>(node);

  // Energy group
  if (check_for_node(node, "energy")) {
    xml_node weightwindow_energy = node.child("energy");
    // energy group for neutron
    if (check_for_node(weightwindow_energy, "neutron")) {
      n_ww = true; // turn on the flag
      n_energy_group = get_node_array<double>(weightwindow_energy, "neutron");
      n_energy_group.insert(n_energy_group.begin(), 0);
    }
          
    // energy group for photon
    if (check_for_node(weightwindow_energy, "photon")) {
      if (!settings::photon_transport) { 
        fatal_error("Photon transport is not on but weight window for photon is used"); 
      }  
      p_ww = true; // turn on the flag
      p_energy_group = get_node_array<double>(weightwindow_energy, "photon");
      p_energy_group.insert(p_energy_group.begin(), 0);
    }   
  } else { 
    fatal_error("Must assign energy group for weight window"); 
  }  

  // read wwinp file
  if (check_for_node(node, "lower_ww")) {
    auto value = get_node_xarray<double>(node, "lower_ww");
    int groups = 0;  // total energy groups for neutron and photon
    if (n_ww) groups += n_energy_group.size() - 1;
    if (p_ww) groups += p_energy_group.size() - 1;
    if ( value.size() != ( mesh_->n_bins()* groups ) ) { 
        fatal_error("The number of lower weight window bounds is not the same as energy-space mesh cell numbers");  
        }
    // neutron lower weight window bound first     
    if (n_ww) {
      for (int i = 0; i < mesh_->n_bins()*(n_energy_group.size()-1); i++) {   
        n_ww_lower.push_back(value.at(i));
      }  
    } 
        // photon lower weight window bound later
    if (p_ww) {
      auto offset = mesh_->n_bins()*(n_energy_group.size()-1);
      for (int i = 0; i < mesh_->n_bins()*(p_energy_group.size()-1); i++) {   
        p_ww_lower.push_back(value.at(i + offset));
      }   
    } 
  } else { 
    fatal_error("Must assign lower weight window bound for weight window"); 
  }  

    
  // WWP-- weight window parameters
  // neutron
  if (check_for_node(node, "neutron_parameters")) {
    xml_node neutron_wwp = node.child("neutron_parameters");
        
    // upper weight window
    if (check_for_node(neutron_wwp,"upper")) {
      n_upper_ratio = std::stod(get_node_value(neutron_wwp, "upper"));
      if (n_upper_ratio < 2 ) fatal_error("Upper to lower weight window ratio must bigger than 2.");
    }
      
    // survival weight window
    if (check_for_node(neutron_wwp, "survival")) {
      n_survival_ratio = std::stod(get_node_value(neutron_wwp, "survival"));
      if (n_survival_ratio <= 1 || n_survival_ratio >= n_upper_ratio )  
        fatal_error("Survival to lower weight window ratio must bigger than 1 "
                    "and less than the upper to lower weight window ratio.");
    }
      
    // max split
    if (check_for_node(neutron_wwp, "max_split")) {
      n_max_split = std::stoi(get_node_value(neutron_wwp, "max_split"));
      if (n_max_split <= 1 ) fatal_error("Max split number must larger than 1.");
    }
      
    // multiplier for weight window lower bounds
    if (check_for_node(neutron_wwp, "multiplier")) {
      n_multiplier = std::stod(get_node_value(neutron_wwp, "multiplier"));
      if (n_multiplier <= 0 ) fatal_error("Multiplier for lower weight window must larger than 0.");
    }      
  }
  // neutron      
        
  // photon
  if (check_for_node(node, "photon_parameters")) {
    xml_node photon_wwp = node.child("photon_parameters");
        
    // upper weight window
    if (check_for_node(photon_wwp, "upper")) {
      p_upper_ratio = std::stod(get_node_value(photon_wwp, "upper"));
      if (p_upper_ratio < 2 ) fatal_error("Upper to lower weight window ratio must larger than 2.");
    }
      
    // survival weight window
    if (check_for_node(photon_wwp, "survival")) {
      p_survival_ratio = std::stod(get_node_value(photon_wwp, "survival"));
      if (p_survival_ratio <= 1 || p_survival_ratio >= p_upper_ratio )  
        fatal_error("Survival to lower weight window ratio must bigger than 1 "
                    "and less than the upper to lower weight window ratio.");
    }
      
    // max split
    if (check_for_node(photon_wwp, "max_split")) {
      p_max_split = std::stoi(get_node_value(photon_wwp, "max_split"));
      if (p_max_split <= 1 ) fatal_error("Max split number must larger than 1.");
    }
      
    // multiplier for weight window lower bounds
    if (check_for_node(photon_wwp, "multiplier")) {
      p_multiplier = std::stod(get_node_value(photon_wwp, "multiplier"));
      if (p_multiplier <= 0 ) fatal_error("Multiplier for lower weight window must larger than 0.");
    }     
  }
  // photon
    
  // user defined source weight biasing in energy
  if (check_for_node(node, "user_defined_biasing")) {
    user_defined_biasing = true;
    xml_node node_user_defined_biasing = node.child("user_defined_biasing");
    
    // energy group for source weight biasing
    if (check_for_node(node_user_defined_biasing, "biasing_energy")) {
      biasing_energy = get_node_array<double>(node_user_defined_biasing, "biasing_energy");
      biasing_energy.insert(biasing_energy.begin(), 0.);
    } else {
      fatal_error("Must provide energy groups for biasing.");
    }

    // origin probability for each energy group
    if (check_for_node(node_user_defined_biasing, "origin_probability")) {
      auto value = get_node_xarray<double>(node_user_defined_biasing, "origin_probability");
      if (value.size() != biasing_energy.size() - 1) fatal_error("Origin probabilities and biasing energies must have the same number of values.");
      origin_probability.push_back(0.);
      cumulative_probability.push_back(0.);
      for (int i = 0; i < value.size(); i++) {
        origin_probability.push_back(value.at(i));
        cumulative_probability.push_back(0.);
      }
      // normalization
      double total_probability = std::accumulate(origin_probability.begin(), origin_probability.end(), 0.0);
      for (int i = 1; i < origin_probability.size(); i++) {   
        origin_probability.at(i) = origin_probability.at(i)/total_probability;
        cumulative_probability.at(i) = cumulative_probability.at(i-1) + origin_probability.at(i);
      }
    } else {
      fatal_error("Must provide origin_probability for each group.");
    }
      
    // biasing weight for each energy group
    if (check_for_node(node_user_defined_biasing, "biasing")) {
      auto value = get_node_xarray<double>(node_user_defined_biasing, "biasing");
      if (value.size() != biasing_energy.size() - 1) {
        fatal_error("Biasing and biasing energies must have the same number of values.");
      }
      biasing.push_back(0);
      cumulative_biasing.push_back(0);
      for (double v : value) {
        biasing.push_back(v);
        cumulative_biasing.push_back(0);
      }
      // normalization
      double total_probability = std::accumulate(biasing.begin(), biasing.end(), 0.0);
      for (int i = 1; i < biasing.size(); i++) {
        biasing[i] = biasing[i]/total_probability;
        cumulative_biasing[i] = cumulative_biasing[i-1] + biasing[i];
      }
    } else {
      fatal_error("Must provide biasing for each energy group.");
    }
  }
}
  
//! source weight biasing in energy
void WeightWindowMesh::weight_biasing(SourceSite& site, uint64_t* seed) 
{
  int i = 0;
  double r = prn(seed);
  for (i = 0; i < cumulative_biasing.size() - 1; i++) 
    if (cumulative_biasing[i] <= r && r < cumulative_biasing[i+1])  break;
  site.E = biasing_energy[i] + ( biasing_energy[i+1] - biasing_energy[i] ) * prn(seed);
  site.wgt *= origin_probability[i+1] / biasing[i+1];
}
  
//! Get weight windows parameters given particle
WeightWindowMesh::WWParams WeightWindowMesh::get_params(Particle& p) const
{	
  // Particle's position and energy
  Position pos  = p.r();
  double E = p.E();
  
  std::vector<double> energy_group;
  std::vector<double> ww_lower;
  WWParams params;
	
  // Determine which set of weight window values to be used based on particle type
  if (p.type() == ParticleType::neutron) {
    energy_group = n_energy_group;
    ww_lower = n_ww_lower;
    params.lower_weight = n_multiplier;       
    params.upper_weight = n_upper_ratio;
    params.survival_weight = n_survival_ratio;
    params.max_split = n_max_split;
  } else if (p.type() == ParticleType::photon) {
    energy_group = p_energy_group;
    ww_lower = p_ww_lower;
    params.lower_weight = p_multiplier;       
    params.upper_weight = p_upper_ratio;
    params.survival_weight = p_survival_ratio;
    params.max_split = p_max_split;
  }

  // get the mesh bin in energy group
  int energy_bin = lower_bound_index(energy_group.begin(), energy_group.end(), E);

  // indices in weight window vector
  int indices = mesh_->get_bin(pos);
  indices += energy_bin*mesh_->n_bins(); 

  params.lower_weight = params.lower_weight*ww_lower[indices];  // equal to multiplier * lower weight window bound (from input file)
  params.upper_weight = params.lower_weight*params.upper_weight;           // equal to multiplied lower weight window bound * upper/lower ratio
  params.survival_weight = params.lower_weight*params.survival_weight;     // equal to multiplied lower weight window bound * survival/lower ratio
  
  return params;
}

} // namespace openmc