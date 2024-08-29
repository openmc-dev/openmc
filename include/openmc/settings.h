#ifndef OPENMC_SETTINGS_H
#define OPENMC_SETTINGS_H

//! \file settings.h
//! \brief Settings for OpenMC

#include <cstdint>
#include <string>
#include <unordered_set>

#include "pugixml.hpp"

#include "openmc/array.h"
#include "openmc/constants.h"
#include "openmc/vector.h"

namespace openmc {

// Type of surface source write
enum class SSWCellType {
  None,
  Both,
  From,
  To,
};

//==============================================================================
// Global variable declarations
//==============================================================================

namespace settings {

// Boolean flags
extern bool assume_separate;      //!< assume tallies are spatially separate?
extern bool check_overlaps;       //!< check overlaps in geometry?
extern bool confidence_intervals; //!< use confidence intervals for results?
extern bool
  create_fission_neutrons; //!< create fission neutrons (fixed source)?
extern bool create_delayed_neutrons; //!< create delayed fission neutrons?
extern "C" bool cmfd_run;            //!< is a CMFD run?
extern bool
  delayed_photon_scaling;   //!< Scale fission photon yield to include delayed
extern "C" bool entropy_on; //!< calculate Shannon entropy?
extern "C" bool
  event_based; //!< use event-based mode (instead of history-based)
extern bool legendre_to_tabular; //!< convert Legendre distributions to tabular?
extern bool material_cell_offsets; //!< create material cells offsets?
extern "C" bool output_summary;    //!< write summary.h5?
extern bool output_tallies;        //!< write tallies.out?
extern bool particle_restart_run;  //!< particle restart run?
extern "C" bool photon_transport;  //!< photon transport turned on?
extern "C" bool reduce_tallies;    //!< reduce tallies at end of batch?
extern bool res_scat_on;           //!< use resonance upscattering method?
extern "C" bool restart_run;       //!< restart run?
extern "C" bool run_CE;            //!< run with continuous-energy data?
extern bool source_latest;         //!< write latest source at each batch?
extern bool source_separate;       //!< write source to separate file?
extern bool source_write;          //!< write source in HDF5 files?
extern bool source_mcpl_write;     //!< write source in mcpl files?
extern bool surf_source_write;     //!< write surface source file?
extern bool surf_mcpl_write;       //!< write surface mcpl file?
extern bool surf_source_read;      //!< read surface source file?
extern bool survival_biasing;      //!< use survival biasing?
extern bool temperature_multipole; //!< use multipole data?
extern "C" bool trigger_on;        //!< tally triggers enabled?
extern bool trigger_predict;       //!< predict batches for triggers?
extern bool ufs_on;                //!< uniform fission site method on?
extern bool urr_ptables_on;        //!< use unresolved resonance prob. tables?
extern "C" bool weight_windows_on; //!< are weight windows are enabled?
extern bool weight_window_checkpoint_surface;   //!< enable weight window check
                                                //!< upon surface crossing?
extern bool weight_window_checkpoint_collision; //!< enable weight window check
                                                //!< upon collision?
extern bool write_all_tracks;     //!< write track files for every particle?
extern bool write_initial_source; //!< write out initial source file?

// Paths to various files
extern std::string path_cross_sections; //!< path to cross_sections.xml
extern std::string path_input;  //!< directory where main .xml files resides
extern std::string path_output; //!< directory where output files are written
extern std::string path_particle_restart; //!< path to a particle restart file
extern std::string path_sourcepoint;      //!< path to a source file
extern std::string path_statepoint;       //!< path to a statepoint file
extern std::string weight_windows_file;   //!< Location of weight window file to
                                          //!< load on simulation initialization

// This is required because the c_str() may not be the first thing in
// std::string. Sometimes it is, but it seems libc++ may not be like that
// on some computers, like the intel Mac.
extern "C" const char* path_statepoint_c; //!< C pointer to statepoint file name

extern "C" int32_t n_inactive;         //!< number of inactive batches
extern "C" int32_t max_lost_particles; //!< maximum number of lost particles
extern double
  rel_max_lost_particles; //!< maximum number of lost particles, relative to the
                          //!< total number of particles
extern "C" int32_t
  max_write_lost_particles;       //!< maximum number of lost particles
                                  //!< to be written to files
extern "C" int32_t gen_per_batch; //!< number of generations per batch
extern "C" int64_t n_particles;   //!< number of particles per generation

extern int64_t
  max_particles_in_flight;      //!< Max num. event-based particles in flight
extern int max_particle_events; //!< Maximum number of particle events
extern ElectronTreatment
  electron_treatment; //!< how to treat secondary electrons
extern array<double, 4>
  energy_cutoff; //!< Energy cutoff in [eV] for each particle type
extern array<double, 4>
  time_cutoff; //!< Time cutoff in [s] for each particle type
extern int
  legendre_to_tabular_points; //!< number of points to convert Legendres
extern int max_order;         //!< Maximum Legendre order for multigroup data
extern int n_log_bins;        //!< number of bins for logarithmic energy grid
extern int n_batches;         //!< number of (inactive+active) batches
extern int n_max_batches;     //!< Maximum number of batches
extern int max_tracks; //!< Maximum number of particle tracks written to file
extern ResScatMethod res_scat_method; //!< resonance upscattering method
extern double res_scat_energy_min; //!< Min energy in [eV] for res. upscattering
extern double res_scat_energy_max; //!< Max energy in [eV] for res. upscattering
extern vector<std::string>
  res_scat_nuclides;           //!< Nuclides using res. upscattering treatment
extern RunMode run_mode;       //!< Run mode (eigenvalue, fixed src, etc.)
extern SolverType solver_type; //!< Solver Type (Monte Carlo or Random Ray)
extern std::unordered_set<int>
  sourcepoint_batch; //!< Batches when source should be written
extern std::unordered_set<int>
  statepoint_batch; //!< Batches when state should be written
extern std::unordered_set<int>
  source_write_surf_id; //!< Surface ids where sources will be written
extern int
  max_history_splits; //!< maximum number of particle splits for weight windows
extern int64_t max_surface_particles; //!< maximum number of particles to be
                                      //!< banked on surfaces per process
extern int64_t ssw_cell_id;           //!< Cell id for the surface source
                                      //!< write setting
extern SSWCellType ssw_cell_type;     //!< Type of option for the cell
                                      //!< argument of surface source write
extern TemperatureMethod
  temperature_method; //!< method for choosing temperatures
extern double
  temperature_tolerance; //!< Tolerance in [K] on choosing temperatures
extern double temperature_default; //!< Default T in [K]
extern array<double, 2>
  temperature_range;           //!< Min/max T in [K] over which to load xs
extern int trace_batch;        //!< Batch to trace particle on
extern int trace_gen;          //!< Generation to trace particle on
extern int64_t trace_particle; //!< Particle ID to enable trace on
extern vector<array<int, 3>>
  track_identifiers;               //!< Particle numbers for writing tracks
extern int trigger_batch_interval; //!< Batch interval for triggers
extern "C" int verbosity;          //!< How verbose to make output
extern double weight_cutoff;       //!< Weight cutoff for Russian roulette
extern double weight_survive;      //!< Survival weight after Russian roulette

} // namespace settings

//==============================================================================
// Functions
//==============================================================================

//! Read settings from XML file
void read_settings_xml();

//! Read settings from XML node
//! \param[in] root XML node for <settings>
void read_settings_xml(pugi::xml_node root);

//! Select temperatures to read based on what is needed and available
void select_temperatures(std::string& object_name,
  const vector<double>& available_temperatures,
  const vector<double>& requested_temperatures, vector<int>& temps_to_read);

void free_memory_settings();

} // namespace openmc

#endif // OPENMC_SETTINGS_H
