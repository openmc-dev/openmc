#ifndef OPENMC_LATTICE_H
#define OPENMC_LATTICE_H

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "hdf5.h"
#include "pugixml.hpp"

#include "openmc/constants.h"
#include "openmc/position.h"


namespace openmc {

//==============================================================================
// Module constants
//==============================================================================

constexpr int32_t NO_OUTER_UNIVERSE{-1};

//==============================================================================
// Global variables
//==============================================================================

class Lattice;
extern std::vector<Lattice*> lattices_c;

extern std::unordered_map<int32_t, int32_t> lattice_map;

//==============================================================================
//! \class Lattice
//! \brief Abstract type for ordered array of universes.
//==============================================================================

class LatticeIter;
class ReverseLatticeIter;

class Lattice
{
public:
  int32_t id;                         //!< Universe ID number
  std::string name;                   //!< User-defined name
  std::vector<int32_t> universes;     //!< Universes filling each lattice tile
  int32_t outer {NO_OUTER_UNIVERSE};  //!< Universe tiled outside the lattice
  std::vector<int32_t> offsets;       //!< Distribcell offset table

  explicit Lattice(pugi::xml_node lat_node);

  virtual ~Lattice() {}

  virtual int32_t& operator[](const int i_xyz[3]) = 0;

  virtual LatticeIter begin();
  LatticeIter end();

  virtual ReverseLatticeIter rbegin();
  ReverseLatticeIter rend();

  //! Convert internal universe values from IDs to indices using universe_map.
  void adjust_indices();

  //! Allocate offset table for distribcell.
  void allocate_offset_table(int n_maps)
  {offsets.resize(n_maps * universes.size(), C_NONE);}

  //! Populate the distribcell offset tables.
  int32_t fill_offset_table(int32_t offset, int32_t target_univ_id, int map);

  //! \brief Check lattice indices.
  //! \param i_xyz[3] The indices for a lattice tile.
  //! \return true if the given indices fit within the lattice bounds.  False
  //!   otherwise.
  virtual bool are_valid_indices(const int i_xyz[3]) const = 0;

  //! \brief Find the next lattice surface crossing
  //! \param r A 3D Cartesian coordinate.
  //! \param u A 3D Cartesian direction.
  //! \param i_xyz[3] The indices for a lattice tile.
  //! \return The distance to the next crossing and an array indicating how the
  //!   lattice indices would change after crossing that boundary.
  virtual std::pair<double, std::array<int, 3>>
  distance(Position r, Direction u, const int i_xyz[3]) const
  = 0;

  //! \brief Find the lattice tile indices for a given point.
  //! \param r A 3D Cartesian coordinate.
  //! \return An array containing the indices of a lattice tile.
  virtual std::array<int, 3> get_indices(Position r) const = 0;

  //! \brief Get coordinates local to a lattice tile.
  //! \param r A 3D Cartesian coordinate.
  //! \param i_xyz[3] The indices for a lattice tile.
  //! \return Local 3D Cartesian coordinates.
  virtual Position
  get_local_position(Position r, const int i_xyz[3]) const = 0;

  //! \brief Check flattened lattice index.
  //! \param indx The index for a lattice tile.
  //! \return true if the given index fit within the lattice bounds.  False
  //!   otherwise.
  virtual bool is_valid_index(int indx) const
  {return (indx >= 0) && (indx < universes.size());}

  //! \brief Get the distribcell offset for a lattice tile.
  //! \param The map index for the target cell.
  //! \param i_xyz[3] The indices for a lattice tile.
  //! \return Distribcell offset i.e. the largest instance number for the target
  //!  cell found in the geometry tree under this lattice tile.
  virtual int32_t& offset(int map, const int i_xyz[3]) = 0;

  //! \brief Convert an array index to a useful human-readable string.
  //! \param indx The index for a lattice tile.
  //! \return A string representing the lattice tile.
  virtual std::string index_to_string(int indx) const = 0;

  //! \brief Write lattice information to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

protected:
  bool is_3d;  //!< Has divisions along the z-axis?

  virtual void to_hdf5_inner(hid_t group_id) const = 0;
};

//==============================================================================
//! An iterator over lattice universes.
//==============================================================================

class LatticeIter
{
public:
  int indx;  //!< An index to a Lattice universes or offsets array.

  LatticeIter(Lattice &lat_, int indx_)
    : lat(lat_),
      indx(indx_)
  {}

  bool operator==(const LatticeIter &rhs) {return (indx == rhs.indx);}

  bool operator!=(const LatticeIter &rhs) {return !(*this == rhs);}

  int32_t& operator*() {return lat.universes[indx];}

  LatticeIter& operator++()
  {
    while (indx < lat.universes.size()) {
      ++indx;
      if (lat.is_valid_index(indx)) return *this;
    }
    indx = lat.universes.size();
    return *this;
  }

protected:
  Lattice &lat;
};

//==============================================================================
//! A reverse iterator over lattice universes.
//==============================================================================

class ReverseLatticeIter : public LatticeIter
{
public:
  ReverseLatticeIter(Lattice &lat_, int indx_)
    : LatticeIter {lat_, indx_}
  {}

  ReverseLatticeIter& operator++()
  {
    while (indx > -1) {
      --indx;
      if (lat.is_valid_index(indx)) return *this;
    }
    indx = -1;
    return *this;
  }
};

//==============================================================================

class RectLattice : public Lattice
{
public:
  explicit RectLattice(pugi::xml_node lat_node);

  int32_t& operator[](const int i_xyz[3]);

  bool are_valid_indices(const int i_xyz[3]) const;

  std::pair<double, std::array<int, 3>>
  distance(Position r, Direction u, const int i_xyz[3]) const;

  std::array<int, 3> get_indices(Position r) const;

  Position
  get_local_position(Position r, const int i_xyz[3]) const;

  int32_t& offset(int map, const int i_xyz[3]);

  std::string index_to_string(int indx) const;

  void to_hdf5_inner(hid_t group_id) const;

private:
  std::array<int, 3> n_cells;    //!< Number of cells along each axis
  Position lower_left;           //!< Global lower-left corner of the lattice
  Position pitch;                //!< Lattice tile width along each axis

  // Convenience aliases
  int &nx {n_cells[0]};
  int &ny {n_cells[1]};
  int &nz {n_cells[2]};
};

//==============================================================================

class HexLattice : public Lattice
{
public:
  explicit HexLattice(pugi::xml_node lat_node);

  int32_t& operator[](const int i_xyz[3]);

  LatticeIter begin();

  ReverseLatticeIter rbegin();

  bool are_valid_indices(const int i_xyz[3]) const;

  std::pair<double, std::array<int, 3>>
  distance(Position r, Direction u, const int i_xyz[3]) const;

  std::array<int, 3> get_indices(Position r) const;

  Position
  get_local_position(Position r, const int i_xyz[3]) const;

  bool is_valid_index(int indx) const;

  int32_t& offset(int map, const int i_xyz[3]);

  std::string index_to_string(int indx) const;

  void to_hdf5_inner(hid_t group_id) const;

private:
  int n_rings;                   //!< Number of radial tile positions
  int n_axial;                   //!< Number of axial tile positions
  Position center;               //!< Global center of lattice
  std::array<double, 2> pitch;   //!< Lattice tile width and height
};

} //  namespace openmc
#endif // OPENMC_LATTICE_H
