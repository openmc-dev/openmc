#include "openmc/settings.h"

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/string_utils.h"
#include "openmc/xml_interface.h"

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

int temperature_method {TEMPERATURE_NEAREST};
bool temperature_multipole {false};
double temperature_tolerance {10.0};
double temperature_default {293.6};
std::array<double, 2> temperature_range {0.0, 0.0};

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

  // Get temperature settings
  if (check_for_node(*root, "temperature_default")) {
    temperature_default = std::stod(get_node_value(*root, "temperature_default"));
  }
  if (check_for_node(*root, "temperature_method")) {
    auto temp_str = get_node_value(*root, "temperature_method", true, true);
    if (temp_str == "nearest") {
      temperature_method = TEMPERATURE_NEAREST;
    } else if (temp_str == "interpolation") {
      temperature_method = TEMPERATURE_INTERPOLATION;
    } else {
      fatal_error("Unknown temperature method: " + temp_str);
    }
  }
  if (check_for_node(*root, "temperature_tolerance")) {
    temperature_tolerance = std::stod(get_node_value(*root, "temperature_tolerance"));
  }
  if (check_for_node(*root, "temperature_multipole")) {
    temperature_multipole = get_node_value_bool(*root, "temperature_multipole");
  }
  if (check_for_node(*root, "temperature_range")) {
    auto range = get_node_array<double>(*root, "temperature_range");
    temperature_range[0] = range[0];
    temperature_range[1] = range[1];
  }
}

} // namespace openmc
