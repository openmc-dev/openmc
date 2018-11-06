#include "openmc/cross_sections.h"

#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/settings.h"
#include "openmc/string_functions.h"
#include "openmc/string_utils.h"
#include "openmc/xml_interface.h"

#include "pugixml.hpp"

#include <cstdlib> // for getenv

namespace openmc {

//==============================================================================
// Global variable declarations
//==============================================================================

std::vector<Library> libraries;
std::map<LibraryKey, std::size_t> library_dict;

extern "C" void read_mg_cross_sections_header();

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
  // Check if materials.xml exists
  std::string filename = settings::path_input + "materials.xml";
  if (!file_exists(filename)) {
    fatal_error("Material XML file '" + filename + "' does not exist.");
  }

  // Parse materials.xml file
  pugi::xml_document doc;
  doc.load_file(filename.c_str());
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
          " user's guide at https://openmc.readthedocs.io for "
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

  // Find the windowed multipole library
  if (settings::run_mode != RUN_MODE_PLOTTING) {
    if (!check_for_node(root, "multipole_library")) {
      // No library location specified in materials.xml, check
      // environment variable
      char* envvar = std::getenv("OPENMC_MULTIPOLE_LIBRARY");
      if (envvar) settings::path_multipole = envvar;
    } else {
      settings::path_multipole = get_node_value(root, "multipole_library");
    }
    if (!ends_with(settings::path_multipole, "/")) {
      settings::path_multipole += "/";
    }
  }

  // Now that the cross_sections.xml or mgxs.h5 has been located, read it in
  if (settings::run_CE) {
    read_ce_cross_sections_xml();
  } else {
    read_mg_cross_sections_header();
  }

  // Establish mapping between (type, material) and index in libraries
  int i = 0;
  for (const auto& lib : libraries) {
    for (const auto& name : lib.materials_) {
      std::string lower_name = name;
      to_lower(lower_name);
      LibraryKey key {lib.type_, lower_name};
      library_dict.insert({key, i});
    }
    ++i;
  }

  // Check that 0K nuclides are listed in the cross_sections.xml file
  for (const auto& name : settings::res_scat_nuclides) {
    std::string lower_name = name;
    to_lower(lower_name);
    LibraryKey key {Library::Type::neutron, lower_name};
    if (library_dict.find(key) == library_dict.end()) {
      fatal_error("Could not find resonant scatterer " +
        name + " in cross_sections.xml file!");
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
    directory = filename.substr(0, pos);
  }

  for (const auto& node_library : root.children("library")) {
    libraries.emplace_back(node_library, directory);
  }

  // Make sure file was not empty
  if (libraries.empty()) {
    fatal_error("No cross section libraries present in cross_sections.xml file.");
  }
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" void library_clear() {
  libraries.clear();
  library_dict.clear();
}

extern "C" const char* library_path(int type, const char* name) {
  auto lib_type = static_cast<Library::Type>(type);
  LibraryKey key {lib_type, name};
  if (library_dict.find(key) == library_dict.end()) {
    return nullptr;
  } else {
    auto idx = library_dict[key];
    return libraries[idx].path_.c_str();
  }
}

extern "C" bool library_present(int type, const char* name) {
  auto lib_type = static_cast<Library::Type>(type);
  LibraryKey key {lib_type, name};
  return library_dict.find(key) != library_dict.end();
}

} // namespace openmc
