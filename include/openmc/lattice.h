#ifndef OPENMC_LATTICE_H
#define OPENMC_LATTICE_H

#include <cstdint>
#include <unordered_map>

#include "hdf5.h"
#include "pugixml.hpp"

#include "openmc/array.h"
#include "openmc/constants.h"
#include "openmc/memory.h"
#include "openmc/pair.h"
#include "openmc/position.h"
#include "openmc/string.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Module constants
//==============================================================================

constexpr int32_t NO_OUTER_UNIVERSE{-1};

enum class LatticeType {
  rect, hex
};

//==============================================================================
// Global variables
//==============================================================================

class Lattice;

namespace model {
  extern std::unordered_map<int32_t, int32_t> lattice_map;
  extern vector<unique_ptr<Lattice>> lattices;
} // namespace model

namespace gpu {
extern __constant__ unique_ptr<Lattice>* lattices;
}

//==============================================================================
//! \class Lattice
//! \brief Abstract type for ordered array of universes.
//==============================================================================

class LatticeIter;
class ReverseLatticeIter;

class Lattice
{
public:
  int32_t id_;                         //!< Universe ID number
  string name_;                        //!< User-defined name
  LatticeType type_;
  vector<int32_t> universes_;          //!< Universes filling each lattice tile
  int32_t outer_ {NO_OUTER_UNIVERSE};  //!< Universe tiled outside the lattice
  vector<int32_t> offsets_;            //!< Distribcell offset table

  explicit Lattice(pugi::xml_node lat_node);
  Lattice() = delete;
  Lattice(Lattice const&) = delete;
  Lattice(Lattice&&) = default;

  virtual __host__ ~Lattice() {}

  HD virtual int32_t const& operator[](array<int, 3> const& i_xyz) = 0;

  virtual LatticeIter begin();
  LatticeIter end();

  virtual ReverseLatticeIter rbegin();
  ReverseLatticeIter rend();

  //! Convert internal universe values from IDs to indices using universe_map.
  void adjust_indices();

  //! Allocate offset table for distribcell.
  void allocate_offset_table(int n_maps)
  {offsets_.resize(n_maps * universes_.size(), C_NONE);}

  //! Populate the distribcell offset tables.
  int32_t fill_offset_table(int32_t offset, int32_t target_univ_id, int map,
    std::unordered_map<int32_t, int32_t>& univ_count_memo);

  //! \brief Check lattice indices.
  //! \param i_xyz[3] The indices for a lattice tile.
  //! \return true if the given indices fit within the lattice bounds.  False
  //!   otherwise.
  virtual HD bool are_valid_indices(array<int, 3> const& i_xyz) const = 0;

  //! \brief Find the next lattice surface crossing
  //! \param r A 3D Cartesian coordinate.
  //! \param u A 3D Cartesian direction.
  //! \param i_xyz The indices for a lattice tile.
  //! \return The distance to the next crossing and an array indicating how the
  //!   lattice indices would change after crossing that boundary.
  HD virtual pair<double, array<int, 3>> distance(
    Position r, Direction u, const array<int, 3>& i_xyz) const = 0;

  //! \brief Find the lattice tile indices for a given point.
  //! \param r A 3D Cartesian coordinate.
  //! \return An array containing the indices of a lattice tile.
  virtual HD void get_indices(Position r, Direction u, array<int, 3>& result) const = 0;

  //! \brief Get coordinates local to a lattice tile.
  //! \param r A 3D Cartesian coordinate.
  //! \param i_xyz The indices for a lattice tile.
  //! \return Local 3D Cartesian coordinates.
  virtual HD Position get_local_position(
    Position r, const array<int, 3>& i_xyz) const = 0;

  //! \brief Check flattened lattice index.
  //! \param indx The index for a lattice tile.
  //! \return true if the given index fit within the lattice bounds.  False
  //!   otherwise.
  virtual bool is_valid_index(int indx) const
  {return (indx >= 0) && (indx < universes_.size());}

  //! \brief Get the distribcell offset for a lattice tile.
  //! \param The map index for the target cell.
  //! \param i_xyz[3] The indices for a lattice tile.
  //! \return Distribcell offset i.e. the largest instance number for the target
  //!  cell found in the geometry tree under this lattice tile.
  virtual HD int32_t& offset(int map, array<int, 3> const& i_xyz) = 0;

  //! \brief Get the distribcell offset for a lattice tile.
  //! \param The map index for the target cell.
  //! \param indx The index for a lattice tile.
  //! \return Distribcell offset i.e. the largest instance number for the target
  //!  cell found in the geometry tree for this lattice index.
  virtual HD int32_t offset(int map, int indx) const = 0;

  //! \brief Convert an array index to a useful human-readable string.
  //! \param indx The index for a lattice tile.
  //! \return A string representing the lattice tile.
  virtual std::string index_to_string(int indx) const = 0;

  //! \brief Write lattice information to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

protected:
  bool is_3d_;  //!< Has divisions along the z-axis?

  virtual void to_hdf5_inner(hid_t group_id) const = 0;
};

//==============================================================================
//! An iterator over lattice universes.
//==============================================================================

class LatticeIter
{
public:
  int indx_;  //!< An index to a Lattice universes or offsets array.

  LatticeIter(Lattice &lat, int indx)
    : indx_(indx), lat_(lat)
  {}

  bool operator==(const LatticeIter &rhs) {return (indx_ == rhs.indx_);}

  bool operator!=(const LatticeIter &rhs) {return !(*this == rhs);}

  int32_t& operator*() {return lat_.universes_[indx_];}

  LatticeIter& operator++()
  {
    while (indx_ < lat_.universes_.size()) {
      ++indx_;
      if (lat_.is_valid_index(indx_)) return *this;
    }
    indx_ = lat_.universes_.size();
    return *this;
  }

protected:
  Lattice& lat_;
};

//==============================================================================
//! A reverse iterator over lattice universes.
//==============================================================================

class ReverseLatticeIter : public LatticeIter
{
public:
  ReverseLatticeIter(Lattice &lat, int indx)
    : LatticeIter {lat, indx}
  {}

  ReverseLatticeIter& operator++()
  {
    while (indx_ > -1) {
      --indx_;
      if (lat_.is_valid_index(indx_)) return *this;
    }
    indx_ = -1;
    return *this;
  }
};

//==============================================================================

class RectLattice : public Lattice
{
public:
  explicit RectLattice(pugi::xml_node lat_node);

  HD int32_t const& operator[](array<int, 3> const& i_xyz);

  HD bool are_valid_indices(array<int, 3> const& i_xyz) const;

  HD pair<double, array<int, 3>> distance(
    Position r, Direction u, const array<int, 3>& i_xyz) const;

  HD void get_indices(Position r, Direction u, array<int, 3>& result) const;

  HD Position get_local_position(Position r, const array<int, 3>& i_xyz) const;

  HD int32_t& offset(int map, array<int, 3> const& i_xyz);

  HD int32_t offset(int map, int indx) const;

  std::string index_to_string(int indx) const;

  void to_hdf5_inner(hid_t group_id) const;

private:
  array<int, 3> n_cells_;         //!< Number of cells along each axis
  Position lower_left_;           //!< Global lower-left corner of the lattice
  Position pitch_;                //!< Lattice tile width along each axis
};

//==============================================================================

class HexLattice : public Lattice
{
public:
  explicit HexLattice(pugi::xml_node lat_node);

  HD int32_t const& operator[](array<int, 3> const& i_xyz);

  LatticeIter begin();

  ReverseLatticeIter rbegin();

  HD bool are_valid_indices(array<int, 3> const& i_xyz) const;

  HD pair<double, array<int, 3>> distance(
    Position r, Direction u, const array<int, 3>& i_xyz) const;

  HD void get_indices(Position r, Direction u, array<int, 3>& result) const;

  HD Position get_local_position(Position r, const array<int, 3>& i_xyz) const;

  bool is_valid_index(int indx) const;

  HD int32_t& offset(int map, array<int, 3> const& i_xyz);

  HD int32_t offset(int map, int indx) const;

  std::string index_to_string(int indx) const;

  void to_hdf5_inner(hid_t group_id) const;

private:
  enum class Orientation {
      y, //!< Flat side of lattice parallel to y-axis
      x  //!< Flat side of lattice parallel to x-axis
  };

  //! Fill universes_ vector for 'y' orientation
  void fill_lattice_y(const vector<string>& univ_words);

  //! Fill universes_ vector for 'x' orientation
  void fill_lattice_x(const vector<string>& univ_words);

  int n_rings_;                   //!< Number of radial tile positions
  int n_axial_;                   //!< Number of axial tile positions
  Orientation orientation_;       //!< Orientation of lattice
  Position center_;               //!< Global center of lattice
  array<double, 2> pitch_;        //!< Lattice tile width and height
};

//==============================================================================
// Non-member functions
//==============================================================================

void read_lattices(pugi::xml_node node);

} //  namespace openmc
#endif // OPENMC_LATTICE_H
