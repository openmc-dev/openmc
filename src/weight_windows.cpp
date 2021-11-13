#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/search.h"
#include "openmc/weight_windows.h"
#include "openmc/xml_interface.h"

namespace openmc {

// read the weight window file
void read_variance_reduction_xml() {
  using namespace weight_window;
  using namespace settings;
  using namespace pugi;
      
  std::string filename = path_input + "variance_reduction.xml";
  
  // Parse variance_reduction.xml file                                
  xml_document doc;

  // try the read vr file
  if (file_exists(filename)) {
    auto result = doc.load_file(filename.c_str());
    
    if (!result) {
      fatal_error("Error processing variance_reduction.xml file.");             
    }
    // Get root element                                                                     
    xml_node root = doc.document_element();

    // Display output message
    write_message("Reading Variance Reduction XML file...", 5);

    // check for variance_reduction
    if (check_for_node(root, "variance_reduction")) {
      // check for weight window section
      xml_node variance_reduction = root.child("variance_reduction");
   
      if (check_for_node(variance_reduction,"weight_windows"))
	xml_node weight_windows = variance_reduction.child("weight_windows");
      else
	// Display output message
	warning("variance_reduction file has been read, but no variance reduction has been enabled.");
	
    else
      fatal_error("variance_reduction element is missing from the variance_reduction.xml file");             

    // try and read a weight window
    read_weight_windows(weight_windows);
  }                                                                        
}

// read the weight window section of a vr file
void read_weight_windows(pugi::xml_node node) 
{
  using namespace pugi;

  // read all the meshes in the VR file
  read_meshes(node);

  // make sure we have at least one settings section
  if (check_for_node(node, "settings")) {
    for (auto settings : node.children("settings")) {
      // do i need to check for id?
      int id = get_node_value(node, "id"));

      WWParams paramters = read_settings(settings);
      ww_map[id] = ww_params.size() - 1; // assign the id of the settings object
      ww_params.emplace_back(ww_params); // put the window in the map
    }
  } else {
    fatal_error("No settings element provided in the variance_reduction.xml file.");  
  }

  // make sure we have at least one domain entry
  if (check_for_node(node, "domain")) {
    for (auto domains : node.children("domain")) {
      // do i need to check for id?
      int id = get_node_value(node, "id"));

      WWDomain domain = read_domains(settings);
      ww_domain_map[id] = ww_domains.size() - 1; // assign the id of the settings object
      ww_domains.emplace_back(domain); // put the window in the map
    }
  } else {
    fatal_error("No domain element provided in the variance_reduction.xml file.");  
  }

  // check domains for consistency
  for ( auto domain : ww_domains ) {
    // num spatial*energy bins must match num weight bins
    int num_spatial_bins = model::meshes[domain->_ww_mesh_idx]->n_bins();
    int num_energy_bins = ww_params[domain->_ww_param_idx]->_energy_bounds.size() - 1;
    int num_weight_bins = ww_params[domain->_ww_param_idx]->lower_ww.size();

    if ( num_weight_bins != num_spatial_bins*num_energy_bins ) {
      std::string err_msg = "In domain " + std::string(domain->_id);
      err_msg += " the number of energy/spatial bins (" + std::string(num_energy_bins);
      err_msg += " /" +std::string(num_spatial_bins);
      err_msg += ") does not match the number of weight bins " + std::string(num_weight_bins);
      fatal_error(err_msg);
    }
  }
  settings::weightwindow_on = true;
}
  
//==============================================================================
// WeightWindowMesh implementation
//==============================================================================

using namespace weight_window;

WeightWindowParamters::WeightWindowParameters(pugi::xml_node node) {
  // get the id - we've already checked for existence up stream
  int _id = get_node_value(node, "id"));  
  
  // get the particle type
  if (check_for_node(node,"particle")) {
    std::string particle_type_str = std::string(get_node_value(particle, "type"));
    _particle_type = openmc::str_to_particle_type(particle_type_str);
  } else {
    fatal_error("No particle type specifier applied within the weight window setttings".);
  }

  // get the survival value - optional
  if (check_for_node(parameters, "survival")) { 
    _survival_ratio = std::stod(get_node_value(parameters, "survival"));
    if(_survival_ratio <= 1 ) 
      fatal_error("Survival to lower weight window ratio must bigger than 1"
		  "and less than the upper to lower weight window ratio.");   
    }
  
  // get the max split - optional
  if (check_for_node(parameters, "max_split")) {
    _max_split = std::stod(get_node_value(parameters, "max_split"));
    if(_max_split <= 1) fatal_error("max split must be larger than 1");
  }
  
  // weight cutoff - optional - but default is 1e-38
  if (check_for_node(parameters, "weight_cutoff")) {
    _weight_cutoff = std::stod(get_node_value(parameters, "weight_cutoff"));
    if(_weight_cutoff <= 0) fatal_error("weight_cutoff must be larger than 0");
    if(_weight_cutoff > 1) fatal_error("weight_cutoff must be less than 1");
  }

  // energy bounds 
  if (!check_for_node(node, "energy")) {
    fatal_error("<energy> section is missing from the " 
      "weight_windows.xml file.");
  } else {
    _energy_bounds = get_node_array<double>(node, "energy");   
  }

  // read the lower weight bounds
  if (!check_for_node(node, "lower_ww")) {
    fatal_error("<lower_ww> section is missing from the " 
      "variance_reduction.xml file.");
  } else {
    _lower_ww = get_node_array<double>(node, "lower_ww");    
  }

  // read the upper weight bounds
  if (!check_for_node(node, "upper_ww")) {
    fatal_error("<upper_ww> section is missing from the " 
      "variance_reduction.xml file.");
  } else {   
    _upper_ww = get_node_array<double>(node, "upper_ww");  
  }

  // make sure that the upper and lower bounds have the right size
  if ( _upper_ww.size() != _lower_ww.size() ) {
    fatal_error("The upper and lower weight window lengths do not match");
  }
}

// construct the weight window domain
WeightWindowDomain::WeightWindowDomain(const int32_t domain_id,
				       const int32_t mesh_idx,
				       const int32_t param_idx) {
  _ww_domain_id = domain_id;
  _ww_mesh_idx = mesh_idx;
  _ww_param_ixd = param_idx;
}

// read the specific weight window settings
WeightWindowDomain::WeightWindowDomain(pugi::xml_node node)
{
  WeightWindowDomain domain;

  // get the mesh id
  if (check_for_node(node,"mesh")) {
    int32_t mesh_id = get_node_value(particle, "mesh");
  } else {
    fatal_error("No mesh specifier in the domain.");
  }

  // get the settings id
  if (check_for_node(node,"settings")) {
    int32_t settings_id = get_node_value(particle, "settings");
  } else {
    fatal_error("No settings specifier in the domain .");
  }

  // set the indices to save time later
  int32_t mesh_idx = mesh_map[mesh_id];
  int32_t settings_idx = ww_map[settings_id];
    
  // construct it
  WeightWindowDomain(mesh_idx,settings_idx);
}

//! Get weight windows parameters given particle - essentially
// given a location tell the particle the right bounds
ParticleWeightParams WeightWindowDomain:get_params(Particle& p, bool& in_domain) const
{
  in_domain = false;

  // check for particle flavour
  ParticleType type = p.type();
  if ( ww_params[_ww_param_idx]->_particle_type != type ) {
    return ParticleWeightParams();
  }
  
  // Particle's position 
  Position pos  = p.r();
  // todo access vector of meshes
  int indices = model::meshes[_ww_mesh_idx]->get_bin(pos);
  
  // no ww settings found - return in
  if ( indices < 0 ) return ParticleWeightParams();

  // particle energy
  double E = p.E();
  double min_e, max_e;
    
  const vector<double> energy_bounds = ww_params[_ww_param_idx]->_energy_bounds;
    
  // find the min and max energy values
  const auto [min_e,max_e] = std::minmax_element(begin(energy_bounds),
						 end(energy_bounds));

  // check to make sure energy is in range
  if ( E < min_e || E > max_e ) return ParticleWeightParams();
  
  // get the mesh bin in energy group
  int energy_bin = lower_bound_index(energy_bounds.begin(), 
				     energy_bounds.end(), E);

  // indices now points to the correct weight given 
  // an energy
  indices += energy_bin*mesh_->n_bins();  

  in_domain = true;
  return ParticleWeightParams(ww_params[_ww_params_idx],indices);
}

} // namespace openmc
