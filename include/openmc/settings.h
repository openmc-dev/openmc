#ifndef OPENMC_SETTINGS_H
#define OPENMC_SETTINGS_H

//! \file settings.h
//! \brief Settings for OpenMC

#include <array>
#include <cstdint>
#include <string>

#include "pugixml.hpp"

namespace openmc {

//==============================================================================
// Global variable declarations
//==============================================================================

namespace settings {

// Boolean flags
extern "C" bool assume_separate;         //!< assume tallies are spatially separate?
extern "C" bool check_overlaps;          //!< check overlaps in geometry?
extern "C" bool cmfd_run;                //!< use CMFD?
extern "C" bool confidence_intervals;    //!< use confidence intervals for results?
extern "C" bool create_fission_neutrons; //!< create fission neutrons (fixed source)?
extern "C" bool entropy_on;              //!< calculate Shannon entropy?
extern "C" bool legendre_to_tabular;     //!< convert Legendre distributions to tabular?
extern "C" bool output_summary;          //!< write summary.h5?
extern "C" bool output_tallies;          //!< write tallies.out?
extern "C" bool particle_restart_run;    //!< particle restart run?
extern "C" bool photon_transport;        //!< photon transport turned on?
extern "C" bool reduce_tallies;          //!< reduce tallies at end of batch?
extern "C" bool res_scat_on;             //!< use resonance upscattering method?
extern "C" bool restart_run;             //!< restart run?
extern "C" bool run_CE;                  //!< run with continuous-energy data?
extern "C" bool source_latest;           //!< write latest source at each batch?
extern "C" bool source_separate;         //!< write source to separate file?
extern "C" bool source_write;            //!< write source in HDF5 files?
extern "C" bool survival_biasing;        //!< use survival biasing?
extern "C" bool temperature_multipole;   //!< use multipole data?
extern "C" bool trigger_on;              //!< tally triggers enabled?
extern "C" bool trigger_predict;         //!< predict batches for triggers?
extern "C" bool ufs_on;                  //!< uniform fission site method on?
extern "C" bool urr_ptables_on;          //!< use unresolved resonance prob. tables?
extern "C" bool write_all_tracks;        //!< write track files for every particle?
extern "C" bool write_initial_source;    //!< write out initial source file?
extern "C" bool dagmc;                   //!< indicator of DAGMC geometry

// Paths to various files
extern std::string path_cross_sections;   //!< path to cross_sections.xml
extern std::string path_input;            //!< directory where main .xml files resides
extern std::string path_multipole;        //!< directory containing multipole files
extern std::string path_output;           //!< directory where output files are written
extern std::string path_particle_restart; //!< path to a particle restart file
extern std::string path_source;
extern std::string path_sourcepoint;      //!< path to a source file
extern std::string path_statepoint;       //!< path to a statepoint file

extern "C" int32_t index_entropy_mesh;  //!< Index of entropy mesh in global mesh array
extern "C" int32_t index_ufs_mesh;      //!< Index of UFS mesh in global mesh array

extern "C" int32_t n_batches;      //!< number of (inactive+active) batches
extern "C" int32_t n_inactive;     //!< number of inactive batches
extern "C" int32_t gen_per_batch;  //!< number of generations per batch
extern "C" int64_t n_particles;    //!< number of particles per generation

extern "C" int electron_treatment;       //!< how to treat secondary electrons
extern "C" double energy_cutoff[4];      //!< Energy cutoff in [eV] for each particle type
extern "C" int legendre_to_tabular_points; //!< number of points to convert Legendres
extern "C" int max_order;                //!< Maximum Legendre order for multigroup data
extern "C" int n_log_bins;               //!< number of bins for logarithmic energy grid
extern "C" int n_max_batches;            //!< Maximum number of batches

extern "C" int res_scat_method;          //!< resonance upscattering method
extern "C" double res_scat_energy_min;   //!< Min energy in [eV] for res. upscattering
extern "C" double res_scat_energy_max;   //!< Max energy in [eV] for res. upscattering
extern "C" int run_mode;                 //!< Run mode (eigenvalue, fixed src, etc.)
extern "C" int temperature_method;       //!< method for choosing temperatures
extern "C" double temperature_tolerance; //!< Tolerance in [K] on choosing temperatures
extern "C" double temperature_default;   //!< Default T in [K]
extern "C" double temperature_range[2];  //!< Min/max T in [K] over which to load xs
extern "C" int trace_batch;              //!< Batch to trace particle on
extern "C" int trace_gen;                //!< Generation to trace particle on
extern "C" int64_t trace_particle;       //!< Particle ID to enable trace on
extern "C" int trigger_batch_interval;   //!< Batch interval for triggers
extern "C" int verbosity;                //!< How verbose to make output
extern "C" double weight_cutoff;         //!< Weight cutoff for Russian roulette
extern "C" double weight_survive;        //!< Survival weight after Russian roulette
} // namespace settings

//! Read settings from XML file
//! \param[in] root XML node for <settings>
extern "C" void read_settings_xml();

extern "C" void read_settings_xml_f(pugi::xml_node_struct* root_ptr);

} // namespace openmc

#endif // OPENMC_SETTINGS_H
