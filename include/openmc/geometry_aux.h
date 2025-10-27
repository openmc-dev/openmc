//! \file geometry_aux.h
//! Auxilary functions for geometry initialization and general data handling.

#ifndef OPENMC_GEOMETRY_AUX_H
#define OPENMC_GEOMETRY_AUX_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "openmc/vector.h"
#include "openmc/xml_interface.h"

namespace openmc {

namespace model {
extern std::unordered_map<int32_t, int32_t> universe_level_counts;
} // namespace model

//! Read geometry from XML file
void read_geometry_xml();

//! Read geometry from XML node
//! \param[in] root node of geometry XML element
void read_geometry_xml(pugi::xml_node root);

//==============================================================================
//! Replace Universe, Lattice, and Material IDs with indices.
//==============================================================================

void adjust_indices();

//==============================================================================
//! Assign defaults to cells with undefined temperatures.
//==============================================================================

void assign_temperatures();

//==============================================================================
//! Finalize densities (compute density multipliers).
//==============================================================================

void finalize_cell_densities();

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

void get_temperatures(
  vector<vector<double>>& nuc_temps, vector<vector<double>>& thermal_temps);

//==============================================================================
//! \brief Perform final setup for geometry
//==============================================================================

void finalize_geometry();

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
//! \param user_distribcells A set of cell indices to create distribcell data
//!   structures for regardless of whether or not they are part of a tally
//!   filter.
//==============================================================================

void prepare_distribcell(
  const std::vector<int32_t>* user_distribcells = nullptr);

//==============================================================================
//! Recursively search through the geometry and count universe instances.
//!
//! This function will update Universe.n_instances_ for each
//! universe in the geometry.
//==============================================================================

void count_universe_instances();

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

std::string distribcell_path(
  int32_t target_cell, int32_t map, int32_t target_offset);

//==============================================================================
//! Determine the maximum number of nested coordinate levels in the geometry.
//! \param univ The index of the universe to begin seraching from (probably the
//!   root universe).
//! \return The number of coordinate levels.
//==============================================================================

int maximum_levels(int32_t univ);

//==============================================================================
//! Check whether or not a universe is the root universe using its ID.
//! \param univ_id The ID of the universe to check.
//! \return Whether or not it is the root universe.
//==============================================================================

bool is_root_universe(int32_t univ_id);

//==============================================================================
//! Deallocates global vectors and maps for cells, universes, and lattices.
//==============================================================================

void free_memory_geometry();

} // namespace openmc
#endif // OPENMC_GEOMETRY_AUX_H
