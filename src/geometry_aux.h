//! \file geometry_aux.h
//! Auxilary functions for geometry initialization and general data handling.

#ifndef GEOMETRY_AUX_H
#define GEOMETRY_AUX_H

#include <cstdint>


namespace openmc {

//==============================================================================
//! Replace Universe, Lattice, and Material IDs with indices.
//==============================================================================

extern "C" void adjust_indices_c();

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

} // namespace openmc
#endif // GEOMETRY_AUX_H
