#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/search.h"
#include "openmc/weight_windows.h"
#include "openmc/xml_interface.h"

namespace openmc {

// read the weight window file
void read_weight_window_xml() {
  using namespace weight_window;
  using namespace settings;
  using namespace pugi;

  // Check if weight_window.xml exists                                                         
  std::string filename = path_input + "weight_windows.xml";
  
  // Parse settings.xml file                                                              
  xml_document doc;                                                                       
  if (file_exists(filename)) {
    auto result = doc.load_file(filename.c_str());                                          
    if (!result) {                                                                          
      fatal_error("Error processing weight_windows.xml file.");                                   
    }
    // Get root element                                                                     
    xml_node root = doc.document_element();     
    std::cout << " Reading weight window XML file ..." << std::endl;
    // try and read a weight window
    settings::ww_settings = std::make_shared<weight_window::WeightWindow>(root);
  } else {
    settings::ww_settings = std::make_shared<weight_window::WeightWindow>();
  }                                                                       
}

//==============================================================================
// WeightWindowMesh implementation
//==============================================================================

using namespace weight_window;

WeightWindow::WeightWindow() 
{
  // default parameters
  n_ww = false;                  // flag for neutron use weight window
  p_ww = false;                  // flag for photon use weight window
  user_defined_biasing = false;  // flag for weight biasing in energy

  // upper weight window = upper_ratio * lower weight window
  weight_params[ParticleType::neutron].upper_ratio = 5.; 
  // survival weight = survival_ratio * lower weight window  
  weight_params[ParticleType::neutron].survival_ratio = 3.; 
  // max number of split particles
  weight_params[ParticleType::neutron].max_split = 5; 
  // multiplier for weight window lower bounds
  weight_params[ParticleType::neutron].multiplier = 1.0; 
  
  // upper weight window = upper_ratio * lower weight window
  weight_params[ParticleType::photon].upper_ratio = 5.; 
  // survival weight = survival_ratio * lower weight window  
  weight_params[ParticleType::photon].survival_ratio = 3.; 
  // max number of split particles
  weight_params[ParticleType::photon].max_split = 5; 
  // multiplier for weight window lower bounds
  weight_params[ParticleType::photon].multiplier = 1.0; 

}

WeightWindow::WeightWindow(pugi::xml_node node)
{
  using namespace pugi;

  WeightWindow(); // call the default constructor
  
  // check for a mesh
  if (check_for_node(node, "mesh")) { 
    xml_node mesh = node.child("mesh");
    mesh_ = make_unique<RectilinearMesh>(mesh);
  }

  // make sure we have a particle section
  if (check_for_node(node, "particle")) {
    // read the particle section
    std::string particle_name = "neutron";
    openmc::ParticleType particle_type;
    if ( particle_name == "neutron") {
      particle_type = openmc::ParticleType::neutron;
      n_ww = true;
    }
    if ( particle_name == "photon") { 
      particle_type = openmc::ParticleType::photon;
      p_ww = true;
    }
    xml_node particle = node.child("particle");
    weight_params[particle_type] = read_particle_settings(particle);
    
  } else {
    fatal_error("Weight window file is missing a particle section");   
  }

  // user defined source weight biasing in energy
  if (check_for_node(node, "user_defined_biasing")) {
    read_user_biasing(node);  
  }
}

// read the particle specific settings 
WWParams WeightWindow::read_particle_settings(pugi::xml_node node) 
{ 
  using namespace pugi;

  WWParams settings;
  
  // 
  if (!check_for_node(node, "parameters")) {
    fatal_error("<parameters> section is missing from the" 
      "weight_windows.xml file.");
  } else {
    pugi::xml_node parameters = node.child("parameters");
  
    // get the upper bound - optional - could be defined as a mesh parameter
    settings.upper_ratio = std::stod(get_node_value(parameters, "upper_bound"));
    if(settings.upper_ratio < 2) fatal_error("upper bound must be larger than 2");

    // get the survival value
    settings.survival_ratio = std::stod(get_node_value(parameters, "survival"));
    if(settings.survival_ratio <= 1 
      || settings.survival_ratio >= settings.upper_ratio) 
        fatal_error("Survival to lower weight window ratio must bigger than 1"
                  "and less than the upper to lower weight window ratio.");   
     
    // get the max split
    settings.max_split = std::stod(get_node_value(parameters, "max_split"));
    if(settings.max_split <= 1) fatal_error("max split must be larger than 1");
     
    // multiplier
    settings.multiplier = std::stod(get_node_value(parameters, "multiplier"));
    if(settings.multiplier <= 0) fatal_error("multiplier must be larger than 0");
  }  

  // 
  if (!check_for_node(node, "energy")) {
    fatal_error("<energy> section is missing from the " 
      "weight_windows.xml file.");
  } else {
    settings.energy_bounds = get_node_array<double>(node, "energy");   
  }

  // read the lower weight bounds
  if (!check_for_node(node, "lower_ww")) {
    fatal_error("<lower_ww> section is missing from the " 
      "weight_windows.xml file.");
  } else {
    settings.lower_ww = get_node_array<double>(node, "lower_ww");    
    if (settings.lower_ww.size() != 
      mesh_->n_bins()*(settings.energy_bounds.size()-1)) {
         fatal_error("not enough entries in the lower_ww section");
      }
  }

  // read the upper weight bounds
  if (check_for_node(node, "upper_ww")) {
    settings.upper_ww = get_node_array<double>(node, "upper_ww");  
    // ensure there are enough entries
    if (settings.upper_ww.size() != 
      mesh_->n_bins()*(settings.energy_bounds.size()-1)) {
         fatal_error("not enough entries in the upper_ww section");
    }
  }

  // maybe multply out a vector for the upper bounds for consistency

  return settings;
}

// read any source biasing settings
void WeightWindow::read_user_biasing(pugi::xml_node node) {
  using namespace pugi;
  // set user defined biasing
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
  
//! source weight biasing in energy
void WeightWindow::weight_biasing(SourceSite& site, uint64_t* seed) 
{
  int i = 0;
  double r = prn(seed);
  for (i = 0; i < cumulative_biasing.size() - 1; i++) 
    if (cumulative_biasing[i] <= r && r < cumulative_biasing[i+1])  break;
  site.E = biasing_energy[i] + ( biasing_energy[i+1] - biasing_energy[i] ) * prn(seed);
  site.wgt *= origin_probability[i+1] / biasing[i+1];
}
  
//! Get weight windows parameters given particle - essentially
// given a location tell the particle the right bounds
ParticleWeightParams WeightWindow::get_params(Particle& p) const
{	
  // Particle's position and energy
  Position pos  = p.r();
  double E = p.E();

  // get the settings for the weight window
  ParticleType type = p.type();
  WWParams ww_settings = weight_params.at(type);
	
  // get the mesh bin in energy group
  int energy_bin = lower_bound_index(ww_settings.energy_bounds.begin(), 
                     ww_settings.energy_bounds.end(), E);

  // indices in weight window vector
  int indices = mesh_->get_bin(pos);            
  // indices points to the correct weight given 
  // an energy
  indices += energy_bin*mesh_->n_bins(); 

  ParticleWeightParams params;
  // set the weight for the current location
  params.lower_weight = ww_settings.multiplier*
                         ww_settings.lower_ww[indices];  
  // set the upper weight bound
  params.upper_weight = ww_settings.upper_ratio*
                         params.lower_weight;
  // set the survival weight
  params.survival_weight = params.lower_weight*
                          ww_settings.survival_ratio;
  // set the max split
  params.max_split = ww_settings.max_split;
  
  return params;
}

} // namespace openmc
