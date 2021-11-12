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
  std::string filename = path_input + "variance_reduction.xml";
  
  // Parse settings.xml file                                                              
  xml_document doc;                                                                       
  if (file_exists(filename)) {
    auto result = doc.load_file(filename.c_str());                                          
    if (!result) {                                                                          
      fatal_error("Error processing variance_reduction.xml file.");                                   
    }
    // Get root element                                                                     
    xml_node root = doc.document_element();

      // Display output message
    write_message("Reading weight window XML file...", 5);

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
  // upper weight window = upper_ratio * lower weight window
  weight_params[ParticleType::neutron].upper_ratio = 5.; 
  // survival weight = survival_ratio * lower weight window  
  weight_params[ParticleType::neutron].survival_ratio = 3.; 
  // max number of split particles
  weight_params[ParticleType::neutron].max_split = 5; 
  // multiplier for weight window lower bounds
  weight_params[ParticleType::neutron].multiplier = 1.0; 
  // weight cutoff for the problem
  weight_params[ParticleType::neutron].weight_cutoff = WEIGHT_CUTOFF; 
  
  // upper weight window = upper_ratio * lower weight window
  weight_params[ParticleType::photon].upper_ratio = 5.; 
  // survival weight = survival_ratio * lower weight window  
  weight_params[ParticleType::photon].survival_ratio = 3.; 
  // max number of split particles
  weight_params[ParticleType::photon].max_split = 5; 
  // multiplier for weight window lower bounds
  weight_params[ParticleType::photon].multiplier = 1.0;
  // weight cutoff for the problem
  weight_params[ParticleType::photon].weight_cutoff = WEIGHT_CUTOFF; 
}

WeightWindow::WeightWindow(pugi::xml_node node) : WeightWindow()
{
  using namespace pugi;
  
  // check for a mesh
  if (check_for_node(node, "mesh")) { 
    xml_node mesh = node.child("mesh");
    mesh_ = make_unique<RectilinearMesh>(mesh);
  }

  // make sure we have a particle section
  if (check_for_node(node, "particle")) {
    // read the particle specific settings
    xml_node particle = node.child("particle");

    // get the particle type
    std::string particle_type_str = std::string(get_node_value(particle, "type"));
    openmc::ParticleType particle_type = openmc::str_to_particle_type(particle_type_str);

    weight_params[particle_type] = read_particle_settings(particle);
    if ( particle_type_str == "neutron" )
      this->n_ww = true;
    if ( particle_type_str == "photon" )
      this->p_ww = true;
    
  } else {
    fatal_error("Weight window file is missing a particle section");   
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
    if (check_for_node(parameters, "upper_bound")) {
      settings.upper_ratio = std::stod(get_node_value(parameters, "upper_bound"));
      if(settings.upper_ratio < 2) fatal_error("upper bound must be larger than 2");
    }

    // get the survival value - optional
    if (check_for_node(parameters, "survival")) { 
      settings.survival_ratio = std::stod(get_node_value(parameters, "survival"));
      if(settings.survival_ratio <= 1 
        || settings.survival_ratio >= settings.upper_ratio) 
          fatal_error("Survival to lower weight window ratio must bigger than 1"
                    "and less than the upper to lower weight window ratio.");   
    }

    // get the max split - optional
    if (check_for_node(parameters, "max_split")) {
      settings.max_split = std::stod(get_node_value(parameters, "max_split"));
      if(settings.max_split <= 1) fatal_error("max split must be larger than 1");
    }

    // multiplier - optional 
    if (check_for_node(parameters, "multiplier")) {
      settings.multiplier = std::stod(get_node_value(parameters, "multiplier"));
      if(settings.multiplier <= 0) fatal_error("multiplier must be larger than 0");
    }

    // weight cutoff - optional
    if (check_for_node(parameters, "weight_cutoff")) {
      settings.weight_cutoff = std::stod(get_node_value(parameters, "weight_cutoff"));
      if(settings.weight_cutoff <= 0) fatal_error("weight_cutoff must be larger than 0");
      if(settings.weight_cutoff > 1) fatal_error("weight_cutoff must be less than 1");
    }
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
  
//! Get weight windows parameters given particle - essentially
// given a location tell the particle the right bounds
ParticleWeightParams WeightWindow::get_params(Particle& p) const
{	
  // Particle's position and energy
  Position pos  = p.r();
  double E = p.E();

  // get the settings for the weight window
  ParticleType type = p.type();

  // 
  WWParams ww_settings;
  auto it = weight_params.find(p.type());
  // return if

  // if we dont find the particle in the map return
  // in the window - no spliting will happen
  if( it != weight_params.end()) {
    ww_settings = it->second;
  } else {
    // no ww settings found - return in
    return ParticleWeightParams();
  }
	
  // indices in weight window vector
  // todo access vector of meshes
  int indices = mesh_->get_bin(pos);

  // no ww settings found - return in
  if ( indices < 0 ) return ParticleWeightParams();

  // find the min and max energy values
  const auto [min_e,max_e] = std::minmax_element(begin(energy_bounds),end(energy_bounds));

  // check to make sure energy is in range
  if ( E < min_e || E > max_e ) return ParticleWeightParams();
  
  // get the mesh bin in energy group
  int energy_bin = lower_bound_index(ww_settings.energy_bounds.begin(), 
                     ww_settings.energy_bounds.end(), E);

  // indices points to the correct weight given 
  // an energy
  indices += energy_bin*mesh_->n_bins(); 
  
  return ParticleWeightParams(ww_settings,indices);
}

} // namespace openmc
