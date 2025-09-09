#ifndef OPENMC_UNIVERSE_H
#define OPENMC_UNIVERSE_H

#include "openmc/bounding_box.h"
#include "openmc/cell.h"

namespace openmc {

#ifdef OPENMC_DAGMC_ENABLED
class DAGUniverse;
#endif

class GeometryState;
class Universe;
class UniversePartitioner;

namespace model {

extern std::unordered_map<int32_t, int32_t> universe_map;
extern vector<unique_ptr<Universe>> universes;

} // namespace model

//==============================================================================
//! A geometry primitive that fills all space and contains cells.
//==============================================================================

class Universe {
public:
  int32_t id_;            //!< Unique ID
  vector<int32_t> cells_; //!< Cells within this universe
  int filled_with_triso_base_ = -1; //!< ID of cell filled with virtual lattice
  int32_t n_instances_;   //!< Number of instances of this universe

  //! \brief Write universe information to an HDF5 group.
  //! \param group_id An HDF5 group id.
  virtual void to_hdf5(hid_t group_id) const;

  virtual bool find_cell(GeometryState& p) const;

  virtual bool find_cell_in_virtual_lattice(GeometryState& p) const;

  BoundingBox bounding_box() const;

  /* By default, universes are CSG universes. The DAGMC
   * universe overrides standard behaviors, and in the future,
   * other things might too.
   */
  virtual GeometryType geom_type() const { return GeometryType::CSG; }

  unique_ptr<UniversePartitioner> partitioner_;
};

//==============================================================================
//! Speeds up geometry searches by grouping cells in a search tree.
//
//! Currently this object only works with universes that are divided up by a
//! bunch of z-planes.  It could be generalized to other planes, cylinders,
//! and spheres.
//==============================================================================

class UniversePartitioner {
public:
  explicit UniversePartitioner(const Universe& univ);

  //! Return the list of cells that could contain the given coordinates.
  const vector<int32_t>& get_cells(Position r, Direction u) const;

private:
  //! A sorted vector of indices to surfaces that partition the universe
  vector<int32_t> surfs_;

  //! Vectors listing the indices of the cells that lie within each partition
  //
  //! There are n+1 partitions with n surfaces.  `partitions_.front()` gives the
  //! cells that lie on the negative side of `surfs_.front()`.
  //! `partitions_.back()` gives the cells that lie on the positive side of
  //! `surfs_.back()`.  Otherwise, `partitions_[i]` gives cells sandwiched
  //! between `surfs_[i-1]` and `surfs_[i]`.
  vector<vector<int32_t>> partitions_;
};

} // namespace openmc
#endif // OPENMC_UNIVERSE_H
