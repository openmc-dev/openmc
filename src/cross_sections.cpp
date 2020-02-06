#include "openmc/cross_sections.h"

#include "openmc/constants.h"
#include "openmc/container_util.h"
#ifdef DAGMC
#include "openmc/dagmc.h"
#endif
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/string_utils.h"
#include "openmc/thermal.h"
#include "openmc/xml_interface.h"
#include "openmc/wmp.h"

#include "pugixml.hpp"

#include <cstdlib> // for getenv
#include <unordered_set>

namespace openmc {

//==============================================================================
// Global variable declarations
//==============================================================================

namespace data {

std::vector<Library> libraries;
std::map<LibraryKey, std::size_t> library_map;

}

//==============================================================================
// Library methods
//==============================================================================

Library::Library(pugi::xml_node node, const std::string& directory)
{
  // Get type of library
  if (check_for_node(node, "type")) {
    auto type = get_node_value(node, "type");
    if (type == "neutron") {
      type_ = Type::neutron;
    } else if (type == "thermal") {
      type_ = Type::thermal;
    } else if (type == "photon") {
      type_ = Type::photon;
    } else if (type == "wmp") {
      type_ = Type::wmp;
    } else {
      fatal_error("Unrecognized library type: " + type);
    }
  } else {
    fatal_error("Missing library type");
  }

  // Get list of materials
  if (check_for_node(node, "materials")) {
    materials_ = get_node_array<std::string>(node, "materials");
  }

  // determine path of cross section table
  if (!check_for_node(node, "path")) {
    fatal_error("Missing library path");
  }
  std::string path = get_node_value(node, "path");

  if (starts_with(path, "/")) {
    path_ = path;
  } else if (ends_with(directory, "/")) {
    path_ = directory + path;
  } else {
    path_ = directory + "/" + path;
  }

  if (!file_exists(path_)) {
    warning("Cross section library " + path_ + " does not exist.");
  }
}

//==============================================================================
// Non-member functions
//==============================================================================

void read_cross_sections_xml()
{
  pugi::xml_document doc;
  std::string filename = settings::path_input + "materials.xml";
#ifdef DAGMC
  std::string s;
  bool found_uwuw_mats = false;
  if (settings::dagmc) {
    found_uwuw_mats = get_uwuw_materials_xml(s);
  }

  if (found_uwuw_mats) {
    // if we found uwuw materials, load those
    doc.load_file(s.c_str());
  } else {
#endif
  // Check if materials.xml exists
  if (!file_exists(filename)) {
    fatal_error("Material XML file '" + filename + "' does not exist.");
  }
  // Parse materials.xml file
  doc.load_file(filename.c_str());
#ifdef DAGMC
  }
#endif

  auto root = doc.document_element();

  // Find cross_sections.xml file -- the first place to look is the
  // materials.xml file. If no file is found there, then we check the
  // OPENMC_CROSS_SECTIONS environment variable
  if (!check_for_node(root, "cross_sections")) {
    // No cross_sections.xml file specified in settings.xml, check
    // environment variable
    if (settings::run_CE) {
      char* envvar = std::getenv("OPENMC_CROSS_SECTIONS");
      if (!envvar) {
        fatal_error("No cross_sections.xml file was specified in "
          "materials.xml or in the OPENMC_CROSS_SECTIONS"
          " environment variable. OpenMC needs such a file to identify "
          "where to find data libraries. Please consult the"
          " user's guide at https://docs.openmc.org/ for "
          "information on how to set up data libraries.");
      }
      settings::path_cross_sections = envvar;
    } else {
      char* envvar = std::getenv("OPENMC_MG_CROSS_SECTIONS");
      if (!envvar) {
        fatal_error("No mgxs.h5 file was specified in "
              "materials.xml or in the OPENMC_MG_CROSS_SECTIONS environment "
              "variable. OpenMC needs such a file to identify where to "
              "find MG cross section libraries. Please consult the user's "
              "guide at http://openmc.readthedocs.io for information on "
              "how to set up MG cross section libraries.");
      }
      settings::path_cross_sections = envvar;
    }
  } else {
    settings::path_cross_sections = get_node_value(root, "cross_sections");
  }

  // Now that the cross_sections.xml or mgxs.h5 has been located, read it in
  if (settings::run_CE) {
    read_ce_cross_sections_xml();
  } else {
    data::mg.read_header(settings::path_cross_sections);
    put_mgxs_header_data_to_globals();
  }

  // Establish mapping between (type, material) and index in libraries
  int i = 0;
  for (const auto& lib : data::libraries) {
    for (const auto& name : lib.materials_) {
      LibraryKey key {lib.type_, name};
      data::library_map.insert({key, i});
    }
    ++i;
  }

  // Check that 0K nuclides are listed in the cross_sections.xml file
  for (const auto& name : settings::res_scat_nuclides) {
    LibraryKey key {Library::Type::neutron, name};
    if (data::library_map.find(key) == data::library_map.end()) {
      fatal_error("Could not find resonant scatterer " +
        name + " in cross_sections.xml file!");
    }
  }
}

void
read_ce_cross_sections(const std::vector<std::vector<double>>& nuc_temps,
  const std::vector<std::vector<double>>& thermal_temps)
{
  std::unordered_set<std::string> already_read;

  // Construct a vector of nuclide names because we haven't loaded nuclide data
  // yet, but we need to know the name of the i-th nuclide
  std::vector<std::string> nuclide_names(data::nuclide_map.size());
  std::vector<std::string> thermal_names(data::thermal_scatt_map.size());
  for (const auto& kv : data::nuclide_map) {
    nuclide_names[kv.second] = kv.first;
  }
  for (const auto& kv : data::thermal_scatt_map) {
    thermal_names[kv.second] = kv.first;
  }

  // Read cross sections
  for (const auto& mat : model::materials) {
    for (int i_nuc : mat->nuclide_) {
      // Find name of corresponding nuclide. Because we haven't actually loaded
      // data, we don't have the name available, so instead we search through
      // all key/value pairs in nuclide_map
      std::string& name = nuclide_names[i_nuc];

      // If we've already read this nuclide, skip it
      if (already_read.find(name) != already_read.end()) continue;

      LibraryKey key {Library::Type::neutron, name};
      int idx = data::library_map[key];
      std::string& filename = data::libraries[idx].path_;

      write_message("Reading " + name + " from " + filename, 6);

      // Open file and make sure version is sufficient
      hid_t file_id = file_open(filename, 'r');
      check_data_version(file_id);

      // Read nuclide data from HDF5
      hid_t group = open_group(file_id, name.c_str());
      int i_nuclide = data::nuclides.size();
      data::nuclides.push_back(std::make_unique<Nuclide>(
        group, nuc_temps[i_nuc], i_nuclide));

      close_group(group);
      file_close(file_id);

      // Determine if minimum/maximum energy for this nuclide is greater/less
      // than the previous
      if (data::nuclides[i_nuclide]->grid_.size() >= 1) {
        int neutron = static_cast<int>(Particle::Type::neutron);
        data::energy_min[neutron] = std::max(data::energy_min[neutron],
          data::nuclides[i_nuclide]->grid_[0].energy.front());
        data::energy_max[neutron] = std::min(data::energy_max[neutron],
          data::nuclides[i_nuclide]->grid_[0].energy.back());
      }

      // Add name and alias to dictionary
      already_read.insert(name);

      // Check if elemental data has been read, if needed
      std::string element = to_element(name);
      if (settings::photon_transport) {
        if (already_read.find(element) == already_read.end()) {
          // Read photon interaction data from HDF5 photon library
          LibraryKey key {Library::Type::photon, element};
          int idx = data::library_map[key];
          std::string& filename = data::libraries[idx].path_;
          write_message("Reading " + element + " from " + filename, 6);

          // Open file and make sure version is sufficient
          hid_t file_id = file_open(filename, 'r');
          check_data_version(file_id);

          // Read element data from HDF5
          hid_t group = open_group(file_id, element.c_str());
          data::elements.emplace_back(group, data::elements.size());

          // Determine if minimum/maximum energy for this element is greater/less than
          // the previous
          const auto& elem {data::elements.back()};
          if (elem.energy_.size() >= 1) {
            int photon = static_cast<int>(Particle::Type::photon);
            int n = elem.energy_.size();
            data::energy_min[photon] = std::max(data::energy_min[photon],
              std::exp(elem.energy_(1)));
            data::energy_max[photon] = std::min(data::energy_max[photon],
              std::exp(elem.energy_(n - 1)));
          }

          close_group(group);
          file_close(file_id);

          // Add element to set
          already_read.insert(element);
        }
      }

      // Read multipole file into the appropriate entry on the nuclides array
      if (settings::temperature_multipole) read_multipole_data(i_nuclide);
    }
  }

  for (auto& mat : model::materials) {
    for (const auto& table : mat->thermal_tables_) {
      // Get name of S(a,b) table
      int i_table = table.index_table;
      std::string& name = thermal_names[i_table];

      if (already_read.find(name) == already_read.end()) {
        LibraryKey key {Library::Type::thermal, name};
        int idx = data::library_map[key];
        std::string& filename = data::libraries[idx].path_;

        write_message("Reading " + name + " from " + filename, 6);

        // Open file and make sure version matches
        hid_t file_id = file_open(filename, 'r');
        check_data_version(file_id);

        // Read thermal scattering data from HDF5
        hid_t group = open_group(file_id, name.c_str());
        data::thermal_scatt.push_back(std::make_unique<ThermalScattering>(
          group, thermal_temps[i_table]));
        close_group(group);
        file_close(file_id);

        // Add name to dictionary
        already_read.insert(name);
      }
    } // thermal_tables_

    // Finish setting up materials (normalizing densities, etc.)
    mat->finalize();
  } // materials


  // Set up logarithmic grid for nuclides
  for (auto& nuc : data::nuclides) {
    nuc->init_grid();
  }
  int neutron = static_cast<int>(Particle::Type::neutron);
  simulation::log_spacing = std::log(data::energy_max[neutron] /
    data::energy_min[neutron]) / settings::n_log_bins;

  if (settings::photon_transport && settings::electron_treatment == ElectronTreatment::TTB) {
    // Determine if minimum/maximum energy for bremsstrahlung is greater/less
    // than the current minimum/maximum
    if (data::ttb_e_grid.size() >= 1) {
      int photon = static_cast<int>(Particle::Type::photon);
      int n_e = data::ttb_e_grid.size();
      data::energy_min[photon] = std::max(data::energy_min[photon], data::ttb_e_grid(1));
      data::energy_max[photon] = std::min(data::energy_max[photon], data::ttb_e_grid(n_e - 1));
    }

    // Take logarithm of energies since they are log-log interpolated
    data::ttb_e_grid = xt::log(data::ttb_e_grid);
  }

  // Show which nuclide results in lowest energy for neutron transport
  for (const auto& nuc : data::nuclides) {
    // If a nuclide is present in a material that's not used in the model, its
    // grid has not been allocated
    if (nuc->grid_.size() > 0) {
      double max_E = nuc->grid_[0].energy.back();
      int neutron = static_cast<int>(Particle::Type::neutron);
      if (max_E == data::energy_max[neutron]) {
        write_message("Maximum neutron transport energy: " +
          std::to_string(data::energy_max[neutron]) + " eV for " +
          nuc->name_, 7);
        if (mpi::master && data::energy_max[neutron] < 20.0e6) {
          warning("Maximum neutron energy is below 20 MeV. This may bias "
            " the results.");
        }
        break;
      }
    }
  }

  // Show minimum/maximum temperature
  write_message("Minimum neutron data temperature: " +
    std::to_string(data::temperature_min) + " K", 4);
  write_message("Maximum neutron data temperature: " +
    std::to_string(data::temperature_max) + " K", 4);

  // If the user wants multipole, make sure we found a multipole library.
  if (settings::temperature_multipole) {
    bool mp_found = false;
    for (const auto& nuc : data::nuclides) {
      if (nuc->multipole_) {
        mp_found = true;
        break;
      }
    }
    if (mpi::master && !mp_found) {
      warning("Windowed multipole functionality is turned on, but no multipole "
        "libraries were found. Make sure that windowed multipole data is "
        "present in your cross_sections.xml file.");
    }
  }
}

void read_ce_cross_sections_xml()
{
  // Check if cross_sections.xml exists
  const auto& filename = settings::path_cross_sections;
  if (!file_exists(filename)) {
    // Could not find cross_sections.xml file
    fatal_error("Cross sections XML file '" + filename +
      "' does not exist.");
  }

  write_message("Reading cross sections XML file...", 5);

  // Parse cross_sections.xml file
  pugi::xml_document doc;
  auto result = doc.load_file(filename.c_str());
  if (!result) {
    fatal_error("Error processing cross_sections.xml file.");
  }
  auto root = doc.document_element();

  std::string directory;
  if (check_for_node(root, "directory")) {
    // Copy directory information if present
    directory = get_node_value(root, "directory");
  } else {
    // If no directory is listed in cross_sections.xml, by default select the
    // directory in which the cross_sections.xml file resides
    auto pos = filename.rfind("/");
    if (pos == std::string::npos) {
      // no '/' found, probably a Windows directory
      pos = filename.rfind("\\");
    }
    directory = filename.substr(0, pos);
  }

  for (const auto& node_library : root.children("library")) {
    data::libraries.emplace_back(node_library, directory);
  }

  // Make sure file was not empty
  if (data::libraries.empty()) {
    fatal_error("No cross section libraries present in cross_sections.xml file.");
  }
}

void library_clear() {
  data::libraries.clear();
  data::library_map.clear();
}

} // namespace openmc
