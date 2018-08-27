#ifndef OPENMC_SETTINGS_H
#define OPENMC_SETTINGS_H

//! \file settings.h
//! \brief Settings for OpenMC

#include <array>
#include <string>

#include "pugixml.hpp"

namespace openmc {

//==============================================================================
// Global variable declarations
//==============================================================================

// Defined on Fortran side
extern "C" bool openmc_check_overlaps;
extern "C" bool openmc_particle_restart_run;
extern "C" bool openmc_restart_run;
extern "C" bool openmc_write_all_tracks;

// Defined in .cpp
// TODO: Make strings instead of char* once Fortran is gone
extern "C" char* openmc_path_input;
extern "C" char* openmc_path_statepoint;
extern "C" char* openmc_path_sourcepoint;
extern "C" char* openmc_path_particle_restart;
extern std::string path_cross_sections;
extern std::string path_multipole;
extern std::string path_output;
extern std::string path_source;

extern int temperature_method;
extern bool temperature_multipole;
extern double temperature_tolerance;
extern double temperature_default;
extern std::array<double, 2> temperature_range;

//==============================================================================
//! Read settings from XML file
//! \param[in] root XML node for <settings>
//==============================================================================

extern "C" void read_settings(pugi::xml_node* root);

} // namespace openmc

#endif // OPENMC_SETTINGS_H
