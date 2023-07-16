#ifndef OPENMC_UNIVERSE_H
#define OPENMC_UNIVERSE_H

#include "openmc/aabb.h"
#include "openmc/cell.h"

namespace openmc {

#ifdef DAGMC
class DAGUniverse;
#endif

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

  //! \brief Write universe information to an HDF5 group.
  //! \param group_id An HDF5 group id.
  virtual void to_hdf5(hid_t group_id) const;

  virtual bool find_cell(Particle& p) const;

  virtual bool find_cell_in_list(
    const std::vector<int>& cells_to_search, Particle& p) const;
  virtual int find_cell_for_point(
    const std::vector<int>& cells_to_search, const Position& p) const;

  // Get bounding box of the entire universe
  BoundingBox bounding_box() const;
  // Get the bounding box for the partitioner set by the user (not yet
  // supported), or generate an optimal one and return it
  AABB partitioner_bounding_box() const;

  const GeometryType& geom_type() const { return geom_type_; }
  GeometryType& geom_type() { return geom_type_; }

  unique_ptr<UniversePartitioner> partitioner_;

private:
  GeometryType geom_type_ = GeometryType::CSG;
};

//==============================================================================
//! Speeds up geometry searches by grouping cells in a search tree.
//
//! The UniversePartitioner is an abstract class. Thus, the actual partitioning 
//! algorithm is implemented in subclasses like ZPlanePartitioner.
//==============================================================================

class UniversePartitioner {
public:
  //! Although the constructor is left to be defined by the derived class,
  //! it should ideally take a Universe as an argument 
  
  //! Dummy destructor 
  virtual ~UniversePartitioner();

  //! Return the list of cells that could contain the given coordinates.
  virtual const vector<int32_t>& get_cells(Position r, Direction u) const = 0;
  virtual const vector<int32_t>& get_cells_fallback(Position r, Direction u) const;

  virtual void export_to_file(const std::string& path) const;
};

} // namespace openmc
#endif // OPENMC_UNIVERSE_H
