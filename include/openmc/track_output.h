#ifndef OPENMC_TRACK_OUTPUT_H
#define OPENMC_TRACK_OUTPUT_H

#include "openmc/particle.h"

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

//! Open HDF5 track file for writing and create track datatype
void open_track_file();

//! Close HDF5 resources for track file
void close_track_file();

//! Determine whether a given particle should collect/write track information
//
//! \param[in] p  Current particle
//! \return Whether to collect/write track information
bool check_track_criteria(const Particle& p);

//! Create a new track state history for a primary/secondary particle
//
//! \param[in] p  Current particle
void add_particle_track(Particle& p);

//! Store particle's current state
//
//! \param[in] p  Current particle
void write_particle_track(Particle& p);

//! Write full particle state history to HDF5 track file
//
//! \param[in] p  Current particle
void finalize_particle_track(Particle& p);

} // namespace openmc

#endif // OPENMC_TRACK_OUTPUT_H
