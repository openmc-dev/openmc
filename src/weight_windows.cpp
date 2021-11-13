#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/search.h"
#include "openmc/weight_windows.h"
#include "openmc/xml_interface.h"

namespace openmc {

namespace variance_reduction {

std::unordered_map<int32_t, int32_t> ww_domain_map;
openmc::vector<unique_ptr<WeightWindowDomain>> ww_domains;

std::unordered_map<int32_t, int32_t> ww_map;
openmc::vector<unique_ptr<WeightWindowParameters>> ww_params;

} // namespace variance_reduction

// read the weight window file
void read_variance_reduction_xml() {

  std::string filename = settings::path_input + "variance_reduction.xml";

  // Parse variance_reduction.xml file
  pugi::xml_document doc;

  // try the read vr file
  if (file_exists(filename)) {
    auto result = doc.load_file(filename.c_str());

    if (!result) {
      fatal_error("Error processing variance_reduction.xml file.");
    }
    // Get root element
    pugi::xml_node root = doc.document_element();

    // Display output message
    write_message("Reading Variance Reduction XML file...", 5);

    // check for variance_reduction
    if (!check_for_node(root, "variance_reduction")) {
      fatal_error("variance_reduction element is missing from the "
                  "variance_reduction.xml file");
    }
    pugi::xml_node vr_node = root.child("variance_reduction");

    // check for weight window section
    if (check_for_node(vr_node, "weight_windows")) {
      pugi::xml_node weight_windows = vr_node.child("weight_windows");
      read_weight_windows(weight_windows);
    } else {
      // Display output message
      warning("variance_reduction file has been read, but no variance "
              "reduction has been enabled.");
    }
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
      variance_reduction::ww_params.emplace_back(
        std::make_unique<WeightWindowParameters>(settings));
      variance_reduction::ww_map[variance_reduction::ww_params.back()->id()] =
        variance_reduction::ww_params.size() - 1;
    }
  } else {
    warning("No settings element provided in the variance_reduction.xml file.");
  }

  // make sure we have at least one domain entry
  if (check_for_node(node, "domain")) {
    for (auto domain : node.children("domain")) {
      variance_reduction::ww_domains.emplace_back(
        std::make_unique<WeightWindowDomain>(domain));
      variance_reduction::ww_domain_map[variance_reduction::ww_domains.back()
                                          ->id()] =
        variance_reduction::ww_domains.size() - 1;
    }
  } else {
    warning("No domain element provided in the variance_reduction.xml file.");
  }

  // check domains for consistency
  for (const auto& domain : variance_reduction::ww_domains) {
    // num spatial*energy bins must match num weight bins
    int num_spatial_bins = model::meshes[domain->ww_mesh_idx()]->n_bins();
    int num_energy_bins = variance_reduction::ww_params[domain->ww_param_idx()]
                            ->energy_bounds()
                            .size() -
                          1;
    int num_weight_bins =
      variance_reduction::ww_params[domain->ww_param_idx()]->lower_ww().size();

    if ( num_weight_bins != num_spatial_bins*num_energy_bins ) {
      auto err_msg =
        fmt::format("In weight window domain {} the number of spatial "
                    "energy/spatial bins ({}) does not match the number "
                    "of weight bins ({})",
          domain->id(), num_energy_bins, num_weight_bins);
      fatal_error(err_msg);
    }
  }
  settings::weightwindow_on = true;
}

//==============================================================================
// WeightWindowMesh implementation
//==============================================================================

WeightWindowParameters::WeightWindowParameters(pugi::xml_node node)
{

  if (check_for_node(node, "id")) {
    int id_ = std::stoi(get_node_value(node, "id"));

    if (variance_reduction::ww_map.find(id_) !=
        variance_reduction::ww_map.end()) {
      auto msg = fmt::format(
        "Two or more weight window parameters use the same unique ID: {}", id_);
      fatal_error(msg);
    }
  }

  // get the particle type
  if (check_for_node(node,"particle")) {
    std::string particle_type_str = std::string(get_node_value(node, "type"));
    particle_type() = openmc::str_to_particle_type(particle_type_str);
  } else {
    fatal_error(
      "No particle type specifier applied within the weight window setttings.");
  }

  // get the survival value - optional
  if (check_for_node(node, "survival")) {
    survival_ratio() = std::stod(get_node_value(node, "survival"));
    if (survival_ratio() <= 1)
      fatal_error("Survival to lower weight window ratio must bigger than 1"
                  "and less than the upper to lower weight window ratio.");
  }

  // get the max split - optional
  if (check_for_node(node, "max_split")) {
    max_split() = std::stod(get_node_value(node, "max_split"));
    if (max_split() <= 1)
      fatal_error("max split must be larger than 1");
  }

  // weight cutoff - optional - but default is 1e-38
  if (check_for_node(node, "weight_cutoff")) {
    weight_cutoff() = std::stod(get_node_value(node, "weight_cutoff"));
    if (weight_cutoff() <= 0)
      fatal_error("weight_cutoff must be larger than 0");
    if (weight_cutoff() > 1)
      fatal_error("weight_cutoff must be less than 1");
  }

  // energy bounds
  if (!check_for_node(node, "energy")) {
    fatal_error("<energy> section is missing from the "
                "weight_windows.xml file.");
  } else {
    energy_bounds() = get_node_array<double>(node, "energy");
  }

  // read the lower weight bounds
  if (!check_for_node(node, "lower_ww")) {
    fatal_error("<lower_ww> section is missing from the "
                "variance_reduction.xml file.");
  } else {
    lower_ww() = get_node_array<double>(node, "lower_ww");
  }

  // read the upper weight bounds
  if (!check_for_node(node, "upper_ww")) {
    fatal_error("<upper_ww> section is missing from the "
                "variance_reduction.xml file.");
  } else {
    upper_ww() = get_node_array<double>(node, "upper_ww");
  }

  // make sure that the upper and lower bounds have the right size
  if (upper_ww().size() != lower_ww().size()) {
    fatal_error("The upper and lower weight window lengths do not match");
  }
}

// read the specific weight window settings
WeightWindowDomain::WeightWindowDomain(pugi::xml_node node)
{

  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));

    // Check to make sure this ID hasn't been used
    if (variance_reduction::ww_domain_map.find(id_) !=
        variance_reduction::ww_domain_map.end()) {
      auto msg = fmt::format(
        "Two or more weight window domains use the same unique ID: {}", id_);
      fatal_error(msg);
    }
  }

  int32_t mesh_id, settings_id;

  // get the mesh id
  if (check_for_node(node,"mesh")) {
    mesh_id = std::stoi(get_node_value(node, "mesh"));
  } else {
    fatal_error("No mesh specifier in the domain.");
  }

  // get the settings id
  if (check_for_node(node,"settings")) {
    settings_id = std::stoi(get_node_value(node, "settings"));
  } else {
    fatal_error("No settings specifier in the domain .");
  }

  // set the indices to save time later
  ww_mesh_idx_ = model::mesh_map[mesh_id];
  ww_param_idx_ = variance_reduction::ww_map[settings_id];
}

//! Get weight windows parameters given particle - essentially
// given a location tell the particle the right bounds
ParticleWeightParams WeightWindowDomain::get_params(
  Particle& p, bool& in_domain) const
{
  in_domain = false;

  // check for particle flavour
  ParticleType type = p.type();
  if (variance_reduction::ww_params[ww_param_idx()]->particle_type() != type) {
    return ParticleWeightParams();
  }

  // Particle's position
  Position pos  = p.r();
  // todo access vector of meshes
  int indices = model::meshes[ww_mesh_idx()]->get_bin(pos);

  // no ww settings found - return in
  if ( indices < 0 ) return ParticleWeightParams();

  // particle energy
  double E = p.E();

  const vector<double> energy_bounds =
    variance_reduction::ww_params[ww_param_idx()]->energy_bounds();

  // find the min and max energy values
  auto e_minmax =
    std::minmax_element(energy_bounds.begin(), energy_bounds.end());

  // check to make sure energy is in range
  if (E < *(e_minmax.first) || E > *(e_minmax.second))
    return ParticleWeightParams();

  // get the mesh bin in energy group
  int energy_bin =
    lower_bound_index(energy_bounds.begin(), energy_bounds.end(), E);

  // indices now points to the correct weight given
  // an energy
  indices += energy_bin * model::meshes[ww_mesh_idx()]->n_bins();

  in_domain = true;
  return ParticleWeightParams(
    variance_reduction::ww_params[ww_param_idx()], indices);
}

} // namespace openmc
