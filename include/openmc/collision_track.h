#ifndef OPENMC_COLLISION_TRACK_H
#define OPENMC_COLLISION_TRACK_H

#include <string>

namespace openmc {

class Particle;

//! Determine whether a collision satisfies the configured recording criteria.
//!
//! \param id_cell      ID of the cell where the collision occurred
//! \param mt_event     ENDF MT number of the reaction
//! \param nuclide      Name of the nuclide involved in the collision
//! \param id_universe  ID of the universe containing the collision
//! \param id_material  ID of the material where the collision occurred
//! \param energy_loss  Energy lost in the collision (E_before - E_after) [eV]
bool should_record_event(int id_cell, int mt_event, const std::string& nuclide,
  int id_universe, int id_material, double energy_loss);

//! Reserve space in the collision track bank according to user settings.
void reserve_bank_capacity();

//! Write collision track data to disk when the bank is full or the batch ends.
//!
//! \param last_batch Whether the current batch is the final batch of the run
void flush_bank(bool last_batch);

//! Record the current particle as a collision-track entry when applicable.
//!
//! \param particle Particle whose collision should be recorded if eligible
void collision_track_record(Particle& particle);

} // namespace openmc

#endif // OPENMC_COLLISION_TRACK_H
