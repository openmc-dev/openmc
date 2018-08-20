//! \file geometry_aux.h
//! Auxilary functions for geometry initialization and general data handling.

#ifndef OPENMC_GEOMETRY_AUX_H
#define OPENMC_GEOMETRY_AUX_H

#include <cstdint>


namespace openmc {

//==============================================================================
//! Replace Universe, Lattice, and Material IDs with indices.
//==============================================================================

extern "C" void adjust_indices();

//==============================================================================
//! Figure out which Universe is the root universe.
//!
//! This function looks for a universe that is not listed in a Cell::fill or in
//! a Lattice.
//! @return The index of the root universe.
//==============================================================================

extern "C" int32_t find_root_universe();

//==============================================================================
//! Allocate storage in Lattice and Cell objects for distribcell offset tables.
//==============================================================================

extern "C" void allocate_offset_tables(int n_maps);

//==============================================================================
//! Recursively search through the geometry and count cell instances.
//!
//! This function will update the Cell::n_instances value for each cell in the
//! geometry.
//! @param univ_indx The index of the universe to begin searching from (probably
//!   the root universe).
//==============================================================================

extern "C" void count_cell_instances(int32_t univ_indx);

//==============================================================================
//! Recursively search through universes and count universe instances.
//! @param search_univ The index of the universe to begin searching from.
//! @param target_univ_id The ID of the universe to be counted.
//! @return The number of instances of target_univ_id in the geometry tree under
//!   search_univ.
//==============================================================================

extern "C" int
count_universe_instances(int32_t search_univ, int32_t target_univ_id);

//==============================================================================
//! Populate Cell and Lattice distribcell offset tables.
//! @param target_univ_id The ID of the universe to be counted.
//! @param map The index of the distribcell map that defines the offsets for the
//!   target universe.
//==============================================================================

extern "C" void fill_offset_tables(int32_t target_univ_id, int map);

//==============================================================================
//! Find the length necessary for a string to contain a distribcell path.
//! @param target_cell The index of the Cell in the global Cell array.
//! @param map The index of the distribcell mapping corresponding to the target
//!   cell.
//! @param target_offset An instance number for a distributed cell.
//! @param root_univ The index of the root Universe in the global Universe
//!   array.
//! @return The size of a character array needed to fit the distribcell path.
//==============================================================================

extern "C" int
distribcell_path_len(int32_t target_cell, int32_t map, int32_t target_offset,
                     int32_t root_univ);

//==============================================================================
//! Build a character array representing the path to a distribcell instance.
//! @param target_cell The index of the Cell in the global Cell array.
//! @param map The index of the distribcell mapping corresponding to the target
//!   cell.
//! @param target_offset An instance number for a distributed cell.
//! @param root_univ The index of the root Universe in the global Universe
//!   array.
//! @param[out] path The unique traversal through the geometry tree that leads
//!   to the desired instance of the target cell.
//==============================================================================

extern "C" void
distribcell_path(int32_t target_cell, int32_t map, int32_t target_offset,
                 int32_t root_univ, char *path);

//==============================================================================
//! Determine the maximum number of nested coordinate levels in the geometry.
//! @param univ The index of the universe to begin seraching from (probably the
//!   root universe).
//! @return The number of coordinate levels.
//==============================================================================

extern "C" int maximum_levels(int32_t univ);

//==============================================================================
//! Deallocates global vectors and maps for cells, universes, and lattices.
//==============================================================================

extern "C" void free_memory_geometry_c();

} // namespace openmc
#endif // OPENMC_GEOMETRY_AUX_H
