//! \file geometry_aux.h
//! Auxilary functions for geometry initialization and general data handling.

#ifndef OPENMC_GEOMETRY_AUX_H
#define OPENMC_GEOMETRY_AUX_H

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

namespace openmc {

namespace model {
  extern std::unordered_map<int32_t, std::unordered_map<int32_t, int32_t>> universe_cell_counts;
  extern std::unordered_map<int32_t, int32_t> universe_level_counts;
} // namespace model

void read_geometry_xml();

//==============================================================================
//! Replace Universe, Lattice, and Material IDs with indices.
//==============================================================================

void adjust_indices();

//==============================================================================
//! Assign defaults to cells with undefined temperatures.
//==============================================================================

void assign_temperatures();

//==============================================================================
//! \brief Obtain a list of temperatures that each nuclide/thermal scattering
//! table appears at in the model. Later, this list is used to determine the
//! actual temperatures to read (which may be different if interpolation is
//! used)
//!
//! \param[out] nuc_temps  Vector of temperatures for each nuclide
//! \param[out] thermal_temps Vector of tempratures for each thermal scattering
//!   table
//==============================================================================

void get_temperatures(std::vector<std::vector<double>>& nuc_temps,
  std::vector<std::vector<double>>& thermal_temps);

//==============================================================================
//! \brief Perform final setup for geometry
//!
//! \param[out] nuc_temps  Vector of temperatures for each nuclide
//! \param[out] thermal_temps Vector of tempratures for each thermal scattering
//!   table
//==============================================================================

void finalize_geometry(std::vector<std::vector<double>>& nuc_temps,
  std::vector<std::vector<double>>& thermal_temps);

//==============================================================================
//! Figure out which Universe is the root universe.
//!
//! This function looks for a universe that is not listed in a Cell::fill or in
//! a Lattice.
//! \return The index of the root universe.
//==============================================================================

int32_t find_root_universe();

//==============================================================================
//! Populate all data structures needed for distribcells.
//==============================================================================

void prepare_distribcell();

//==============================================================================
//! Recursively search through the geometry and count cell instances.
//!
//! This function will update the Cell::n_instances value for each cell in the
//! geometry.
//! \param univ_indx The index of the universe to begin searching from (probably
//!   the root universe).
//==============================================================================

void count_cell_instances(int32_t univ_indx);

//==============================================================================
//! Recursively search through universes and count universe instances.
//! \param search_univ The index of the universe to begin searching from.
//! \param target_univ_id The ID of the universe to be counted.
//! \param univ_count_memo Memoized counts that make this function faster for
//!   large systems.  The first call to this function for each target_univ_id
//!   should start with an empty memo.
//! \return The number of instances of target_univ_id in the geometry tree under
//!   search_univ.
//==============================================================================

int count_universe_instances(int32_t search_univ, int32_t target_univ_id,
  std::unordered_map<int32_t, int32_t>& univ_count_memo);

//==============================================================================
//! Build a character array representing the path to a distribcell instance.
//! \param target_cell The index of the Cell in the global Cell array.
//! \param map The index of the distribcell mapping corresponding to the target
//!   cell.
//! \param target_offset An instance number for a distributed cell.
//! \return The unique traversal through the geometry tree that leads to the
//!   desired instance of the target cell.
//==============================================================================

std::string
distribcell_path(int32_t target_cell, int32_t map, int32_t target_offset);

//==============================================================================
//! Determine the maximum number of nested coordinate levels in the geometry.
//! \param univ The index of the universe to begin seraching from (probably the
//!   root universe).
//! \return The number of coordinate levels.
//==============================================================================

int maximum_levels(int32_t univ);

//==============================================================================
//! Deallocates global vectors and maps for cells, universes, and lattices.
//==============================================================================

void free_memory_geometry();

} // namespace openmc
#endif // OPENMC_GEOMETRY_AUX_H
