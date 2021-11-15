#include "openmc/weight_windows.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/search.h"
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

    // read all the meshes in the VR file
    read_meshes(root);

    // Display output message
    write_message("Reading variance reduction XML file...", 5);

    // check for weight window section
    if (check_for_node(root, "weight_windows")) {
      pugi::xml_node weight_windows = root.child("weight_windows");
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

  // make sure we have at least one settings section
  if (check_for_node(node, "settings")) {
    for (auto settings : node.children("settings")) {
      variance_reduction::ww_params.emplace_back(
        std::make_unique<WeightWindowParameters>(settings));
      variance_reduction::ww_map[variance_reduction::ww_params.back()->id()] =
        variance_reduction::ww_params.size() - 1;
    }
  } else {
    warning("No weight_window_settings element provided in the "
            "variance_reduction.xml file.");
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
    warning("No weight_window_domain element provided in the "
            "variance_reduction.xml file.");
  }

  // check domains for consistency
  for (const auto& domain : variance_reduction::ww_domains) {
    // num spatial*energy bins must match num weight bins
    int num_spatial_bins = model::meshes[domain->mesh_idx()]->n_bins();
    int num_energy_bins = variance_reduction::ww_params[domain->param_idx()]
                            ->energy_bounds()
                            .size() -
                          1;
    int num_weight_bins =
      variance_reduction::ww_params[domain->param_idx()]->lower_ww().size();

    if ( num_weight_bins != num_spatial_bins*num_energy_bins ) {
      auto err_msg =
        fmt::format("In weight window domain {} the number of spatial "
                    "energy/spatial bins ({}) does not match the number "
                    "of weight bins ({})",
          domain->id(), num_energy_bins, num_weight_bins);
      fatal_error(err_msg);
    }
  }
  settings::weight_windows_present = true;
}

//==============================================================================
// WeightWindowMesh implementation
//==============================================================================

WeightWindowParameters::WeightWindowParameters(pugi::xml_node node)
{
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));

    if (variance_reduction::ww_map.find(id_) !=
        variance_reduction::ww_map.end()) {
      auto msg = fmt::format(
        "Two or more weight window parameters use the same unique ID: {}", id_);
      fatal_error(msg);
    }
  }

  // get the particle type
  if (check_for_node(node, "particle")) {
    std::string particle_type_str =
      std::string(get_node_value(node, "particle"));
    particle_type() = openmc::str_to_particle_type(particle_type_str);
  } else {
    fatal_error(
      "No particle type specifier applied within the weight window setttings.");
  }

  // get the survival value - optional
  if (check_for_node(node, "survival_ratio")) {
    survival_ratio() = std::stod(get_node_value(node, "survival_ratio"));
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
  if (check_for_node(node, "energy_bins")) {
    energy_bounds() = get_node_array<double>(node, "energy_bins");
  } else {
    fatal_error("<energy_bins> is missing from the "
                "weight_windows.xml file.");
  }

  // read the lower weight bounds
  if (check_for_node(node, "lower_ww_bounds")) {
    lower_ww() = get_node_array<double>(node, "lower_ww_bounds");
  } else {
    fatal_error("<lower_ww_bounds> section is missing from the "
                "variance_reduction.xml file.");
  }

  // read the upper weight bounds
  if (check_for_node(node, "upper_ww_bounds")) {
    upper_ww() = get_node_array<double>(node, "upper_ww_bounds");
  } else {
    fatal_error("<upper_ww_bounds> node is missing from the "
                "variance_reduction.xml file.");
  }

  // make sure that the upper and lower bounds have the same size
  if (upper_ww().size() != lower_ww().size()) {
    fatal_error("The upper and lower weight window lengths do not match.");
  }
}

void WeightWindowParameters::to_statepoint(hid_t group) const
{
  hid_t params_group =
    create_group(group, "weight_window_parameters " + std::to_string(id_));

  write_dataset(params_group, "particle_type",
    openmc::particle_type_to_str(particle_type_));
  write_dataset(params_group, "energy_bounds", energy_bounds_);
  write_dataset(params_group, "lower_ww_bounds", lower_ww_);
  write_dataset(params_group, "upper_ww_bounds", upper_ww_);
  write_dataset(params_group, "survival_ratio", survival_ratio_);
  write_dataset(params_group, "max_split", max_split_);
  write_dataset(params_group, "weight_cutoff", weight_cutoff_);

  close_group(params_group);
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
  mesh_idx_ = model::mesh_map[mesh_id];
  param_idx_ = variance_reduction::ww_map[settings_id];
}

//! Get weight windows parameters given particle - essentially
// given a location tell the particle the right bounds
bool WeightWindowDomain::find_params(
  Particle& p, ParticleWeightParams& params) const
{

  // check for particle flavour
  ParticleType type = p.type();
  if (variance_reduction::ww_params[param_idx()]->particle_type() != type) {
    return false;
  }

  // todo access vector of meshes
  int indices = model::meshes[mesh_idx()]->get_bin(p.r());

  // no ww settings found - return in
  if (indices < 0)
    return false;

  // particle energy
  double E = p.E();

  const vector<double> energy_bounds =
    variance_reduction::ww_params[param_idx()]->energy_bounds();

  // find the min and max energy values
  auto e_minmax =
    std::minmax_element(energy_bounds.begin(), energy_bounds.end());

  // check to make sure energy is in range
  if (E < *(e_minmax.first) || E > *(e_minmax.second))
    return false;

  // get the mesh bin in energy group
  int energy_bin =
    lower_bound_index(energy_bounds.begin(), energy_bounds.end(), E);

  // indices now points to the correct weight given
  // an energy
  indices += energy_bin * model::meshes[mesh_idx()]->n_bins();

  params =
    ParticleWeightParams(variance_reduction::ww_params[param_idx()], indices);

  return true;
}

void WeightWindowDomain::to_statepoint(hid_t group) const
{
  hid_t domain_group =
    create_group(group, "weight_window_domain " + std::to_string(id_));

  write_dataset(domain_group, "mesh", model::meshes[mesh_idx_]->id_);
  write_dataset(
    domain_group, "settings", variance_reduction::ww_params[param_idx_]->id());

  close_group(domain_group);
}

} // namespace openmc
