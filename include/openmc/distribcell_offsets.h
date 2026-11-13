#ifndef OPENMC_DISTRIBCELL_OFFSETS_H
#define OPENMC_DISTRIBCELL_OFFSETS_H

#include <cstdint>

#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Distributed cell offsets of a fill cell or lattice tile.
//
//! A map numbers the instances of a universe containing distributed cells. The
//! offset for a map is the number of instances of the map's universe before
//! the fill cell in its universe, or before the tile in its lattice. Offsets
//! are only read for maps whose universe is in the cell's fill or the tile, so
//! only those are stored, in the order of the sorted maps of the fill or tile.
//==============================================================================

class DistribcellOffsets {
public:
  DistribcellOffsets() = default;

  //! \param maps Sorted maps in the fill or tile
  //! \param values Offset for each map
  //! Both must outlive this object.
  DistribcellOffsets(const vector<int32_t>& maps, const int32_t* values);

  //! \return Offset for a map, or 0 if the map's universe isn't contained
  int32_t operator[](int32_t map) const
  {
    // Maps are numbered so that the maps in a universe are mostly consecutive,
    // so first look in the consecutive maps at the start of the list
    auto d = static_cast<uint32_t>(map - first_map_);
    if (d < static_cast<uint32_t>(n_first_))
      return values_[d];
    return find(map);
  }

  //! \return Whether the map's universe is contained
  bool contains(int32_t map) const;

private:
  int32_t find(int32_t map) const;

  const int32_t* values_ {nullptr};
  const vector<int32_t>* maps_ {nullptr};
  int32_t first_map_ {0}; //!< First map
  int32_t n_first_ {0};   //!< Number of consecutive maps from the first
};

} // namespace openmc

#endif // OPENMC_DISTRIBCELL_OFFSETS_H
