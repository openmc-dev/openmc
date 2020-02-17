#ifndef OPENMC_SETTINGS_H
#define OPENMC_SETTINGS_H

//! \file settings.h
//! \brief Settings for OpenMC

#include <array>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

#include "pugixml.hpp"

#include "openmc/constants.h"

namespace openmc {

//==============================================================================
// Global variable declarations
//==============================================================================

namespace settings {

// Boolean flags
extern bool assume_separate;          //!< assume tallies are spatially separate?
extern bool check_overlaps;           //!< check overlaps in geometry?
extern bool confidence_intervals;     //!< use confidence intervals for results?
extern bool create_fission_neutrons;  //!< create fission neutrons (fixed source)?
extern "C" bool cmfd_run;             //!< is a CMFD run?
extern "C" bool dagmc;                //!< indicator of DAGMC geometry
extern bool delayed_photon_scaling;   //!< Scale fission photon yield to include delayed
extern "C" bool entropy_on;           //!< calculate Shannon entropy?
extern bool event_based;              //!< use event-based mode (instead of history-based)
extern bool legendre_to_tabular;      //!< convert Legendre distributions to tabular?
extern bool material_cell_offsets;    //!< create material cells offsets?
extern "C" bool output_summary;       //!< write summary.h5?
extern bool output_tallies;           //!< write tallies.out?
extern bool particle_restart_run;     //!< particle restart run?
extern "C" bool photon_transport;     //!< photon transport turned on?
extern "C" bool reduce_tallies;       //!< reduce tallies at end of batch?
extern bool res_scat_on;              //!< use resonance upscattering method?
extern "C" bool restart_run;          //!< restart run?
extern "C" bool run_CE;               //!< run with continuous-energy data?
extern bool source_latest;            //!< write latest source at each batch?
extern bool source_separate;          //!< write source to separate file?
extern bool source_write;             //!< write source in HDF5 files?
extern bool survival_biasing;         //!< use survival biasing?
extern bool temperature_multipole;    //!< use multipole data?
extern "C" bool trigger_on;           //!< tally triggers enabled?
extern bool trigger_predict;          //!< predict batches for triggers?
extern bool ufs_on;                   //!< uniform fission site method on?
extern bool urr_ptables_on;           //!< use unresolved resonance prob. tables?
extern bool write_all_tracks;         //!< write track files for every particle?
extern bool write_initial_source;     //!< write out initial source file?

// Paths to various files
extern std::string path_cross_sections;   //!< path to cross_sections.xml
extern std::string path_input;            //!< directory where main .xml files resides
extern std::string path_output;           //!< directory where output files are written
extern std::string path_particle_restart; //!< path to a particle restart file
extern std::string path_source;
extern std::string path_source_library;   //!< path to the source shared object
extern std::string path_sourcepoint;      //!< path to a source file
extern "C" std::string path_statepoint;   //!< path to a statepoint file

extern "C" int32_t n_batches;                //!< number of (inactive+active) batches
extern "C" int32_t n_inactive;               //!< number of inactive batches
extern "C" int32_t max_lost_particles;     //!< maximum number of lost particles
extern double rel_max_lost_particles;   //!< maximum number of lost particles, relative to the total number of particles
extern "C" int32_t gen_per_batch;            //!< number of generations per batch
extern "C" int64_t n_particles;              //!< number of particles per generation


extern int64_t max_particles_in_flight; //!< Max num. event-based particles in flight

extern ElectronTreatment electron_treatment;       //!< how to treat secondary electrons
extern std::array<double, 4> energy_cutoff;  //!< Energy cutoff in [eV] for each particle type
extern int legendre_to_tabular_points; //!< number of points to convert Legendres
extern int max_order;                //!< Maximum Legendre order for multigroup data
extern int n_log_bins;               //!< number of bins for logarithmic energy grid
extern int n_max_batches;            //!< Maximum number of batches
extern ResScatMethod res_scat_method; //!< resonance upscattering method
extern double res_scat_energy_min;   //!< Min energy in [eV] for res. upscattering
extern double res_scat_energy_max;   //!< Max energy in [eV] for res. upscattering
extern std::vector<std::string> res_scat_nuclides;  //!< Nuclides using res. upscattering treatment
extern RunMode run_mode;                 //!< Run mode (eigenvalue, fixed src, etc.)
extern std::unordered_set<int> sourcepoint_batch; //!< Batches when source should be written
extern std::unordered_set<int> statepoint_batch; //!< Batches when state should be written
extern TemperatureMethod temperature_method;           //!< method for choosing temperatures
extern double temperature_tolerance;     //!< Tolerance in [K] on choosing temperatures
extern double temperature_default;       //!< Default T in [K]
extern std::array<double, 2> temperature_range;  //!< Min/max T in [K] over which to load xs
extern int trace_batch;                  //!< Batch to trace particle on
extern int trace_gen;                    //!< Generation to trace particle on
extern int64_t trace_particle;           //!< Particle ID to enable trace on
extern std::vector<std::array<int, 3>> track_identifiers; //!< Particle numbers for writing tracks
extern int trigger_batch_interval;   //!< Batch interval for triggers
extern "C" int verbosity;                //!< How verbose to make output
extern double weight_cutoff;         //!< Weight cutoff for Russian roulette
extern double weight_survive;        //!< Survival weight after Russian roulette
} // namespace settings

//==============================================================================
// Functions
//==============================================================================

//! Read settings from XML file
//! \param[in] root XML node for <settings>
void read_settings_xml();

void free_memory_settings();

} // namespace openmc

#endif // OPENMC_SETTINGS_H
