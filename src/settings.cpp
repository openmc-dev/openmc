#include "settings.h"

#include "error.h"
#include "openmc.h"
#include "string_utils.h"
#include "xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

char* openmc_path_input;
char* openmc_path_statepoint;
char* openmc_path_sourcepoint;
char* openmc_path_particle_restart;
std::string path_cross_sections;
std::string path_multipole;
std::string path_output;
std::string path_source;

//==============================================================================
// Functions
//==============================================================================

void read_settings(pugi::xml_node* root)
{
  // Look for deprecated cross_sections.xml file in settings.xml
  if (check_for_node(*root, "cross_sections")) {
    warning("Setting cross_sections in settings.xml has been deprecated."
        " The cross_sections are now set in materials.xml and the "
        "cross_sections input to materials.xml and the OPENMC_CROSS_SECTIONS"
        " environment variable will take precendent over setting "
        "cross_sections in settings.xml.");
    path_cross_sections = get_node_value(*root, "cross_sections");
  }

  // Look for deprecated windowed_multipole file in settings.xml
  if (openmc_run_mode != RUN_MODE_PLOTTING) {
    if (check_for_node(*root, "multipole_library")) {
      warning("Setting multipole_library in settings.xml has been "
          "deprecated. The multipole_library is now set in materials.xml and"
          " the multipole_library input to materials.xml and the "
          "OPENMC_MULTIPOLE_LIBRARY environment variable will take "
          "precendent over setting multipole_library in settings.xml.");
      path_multipole = get_node_value(*root, "multipole_library");
    }
    if (!ends_with(path_multipole, "/")) {
      path_multipole += "/";
    }
  }

  // Check for output options
  if (check_for_node(*root, "output")) {

    // Get pointer to output node
    pugi::xml_node node_output = root->child("output");

    // Set output directory if a path has been specified
    if (check_for_node(node_output, "path")) {
      path_output = get_node_value(node_output, "path");
      if (!ends_with(path_output, "/")) {
        path_output += "/";
      }
    }
  }
}

} // namespace openmc