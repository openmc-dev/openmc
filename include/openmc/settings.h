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

namespace settings {

// Boolean flags
extern "C" bool check_overlaps;        //!< check overlaps in geometry?
extern "C" bool particle_restart_run;  //!< particle restart run?
extern "C" bool photon_transport;      //!< photon transport turned on?
extern "C" bool restart_run;           //!< restart run?
extern "C" bool run_CE;                //!< run with continuous-energy data?
extern "C" bool write_all_tracks;      //!< write track files for every particle?
extern "C" bool write_initial_source;  //!< write out initial source file?

// Paths to various files
// TODO: Make strings instead of char* once Fortran is gone
extern "C" char* path_input;
extern "C" char* path_statepoint;
extern "C" char* path_sourcepoint;
extern "C" char* path_particle_restart;
extern std::string path_cross_sections;
extern std::string path_multipole;
extern std::string path_output;
extern std::string path_source;

// Temperature settings
extern int temperature_method;
extern bool temperature_multipole;
extern double temperature_tolerance;
extern double temperature_default;
extern std::array<double, 2> temperature_range;

} // namespace settings

//==============================================================================
//! Read settings from XML file
//! \param[in] root XML node for <settings>
//==============================================================================

extern "C" void read_settings(pugi::xml_node* root);

} // namespace openmc

#endif // OPENMC_SETTINGS_H
