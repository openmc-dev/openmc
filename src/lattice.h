#ifndef LATTICE_H
#define LATTICE_H

#include <array>
#include <cstdint>
#include <limits>  // For numeric_limits
#include <map>
#include <string>
#include <vector>

#include "hdf5.h"
#include "pugixml/pugixml.hpp"


namespace openmc {

//==============================================================================
// Module constants
//==============================================================================

constexpr int32_t NO_OUTER_UNIVERSE{-1};

//==============================================================================
// Global variables
//==============================================================================

//extern "C" int32_t n_lattice;

class Lattice;
//extern Lattice **lattices_c;
extern std::vector<Lattice*> lattices_c;

extern std::map<int32_t, int32_t> lattice_dict;

//==============================================================================
//! Abstract type for ordered array of universes
//==============================================================================

class Lattice
{
public:
  int32_t id;                        //! Universe ID number
  std::string name;                  //! User-defined name
  //std::vector<double> pitch;         //! Pitch along each basis
  std::vector<int32_t> universes;    //! Universes filling each lattice tile
  int32_t outer{NO_OUTER_UNIVERSE};  //! Universe tiled outside the lattice
  //std::vector<int32_t> offset;       //! Distribcell offsets

  explicit Lattice(pugi::xml_node lat_node);

  virtual ~Lattice() {}

  virtual int32_t& operator[](const int i_xyz[3]) = 0;

  //! Convert internal universe values from IDs to indices using universe_dict.
  virtual void adjust_indices() = 0;

  //! Check lattice indices.
  //! @param i_xyz[3] The indices for a lattice tile.
  //! @return true if the given indices fit within the lattice bounds.  False
  //!   otherwise.
  virtual bool are_valid_indices(const int i_xyz[3]) const = 0;

  //! Find the next lattice surface crossing
  //! @param xyz[3] A 3D Cartesian coordinate.
  //! @param uvw[3] A 3D Cartesian direction.
  //! @param i_xyz[3] The indices for a lattice tile.
  //! @return The distance to the next crossing and an array indicating how the
  //!   lattice indices would change after crossing that boundary.
  virtual std::pair<double, std::array<int, 3>>
  distance(const double xyz[3], const double uvw[3], const int i_xyz[3]) const
  = 0;

  //! Find the lattice tile indices for a given point.
  //! @param xyz[3] A 3D Cartesian coordinate.
  //! @return An array containing the indices of a lattice tile.
  virtual std::array<int, 3> get_indices(const double xyz[3]) const = 0;

  //! Get coordinates local to a lattice tile.
  //! @param global_xyz[3] A 3D Cartesian coordinate.
  //! @param i_xyz[3] The indices for a lattice tile.
  //! @return Local 3D Cartesian coordinates.
  virtual std::array<double, 3>
  get_local_xyz(const double global_xyz[3], const int i_xyz[3]) const = 0;

  //! Write all information needed to reconstruct the lattice to an HDF5 group.
  //! @param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

protected:
  bool is_3d;  //! Has divisions along the z-axis

  virtual void to_hdf5_inner(hid_t group_id) const = 0;
};

//==============================================================================
//==============================================================================

class RectLattice : public Lattice
{
public:
  explicit RectLattice(pugi::xml_node lat_node);

  virtual ~RectLattice() {}

  int32_t& operator[](const int i_xyz[3]);

  void adjust_indices();

  bool are_valid_indices(const int i_xyz[3]) const;

  std::pair<double, std::array<int, 3>>
  distance(const double xyz[3], const double uvw[3], const int i_xyz[3]) const;

  std::array<int, 3> get_indices(const double xyz[3]) const;

  std::array<double, 3>
  get_local_xyz(const double global_xyz[3], const int i_xyz[3]) const;

  void to_hdf5_inner(hid_t group_id) const;

protected:
  std::array<int, 3> n_cells;        //! Number of cells along each axis
  std::array<double, 3> lower_left;  //! Global lower-left corner of the lattice
  std::array<double, 3> pitch;       //! Lattice tile width along each axis
};

class HexLattice : public Lattice
{
public:
  explicit HexLattice(pugi::xml_node lat_node);

  virtual ~HexLattice() {}

  int32_t& operator[](const int i_xyz[3]);

  void adjust_indices();

  bool are_valid_indices(const int i_xyz[3]) const;

  std::pair<double, std::array<int, 3>>
  distance(const double xyz[3], const double uvw[3], const int i_xyz[3]) const;

  std::array<int, 3> get_indices(const double xyz[3]) const;

  std::array<double, 3>
  get_local_xyz(const double global_xyz[3], const int i_xyz[3]) const;

  void to_hdf5_inner(hid_t group_id) const;

protected:
  int n_rings;                   //! Number of radial tile positions
  int n_axial;                   //! Number of axial tile positions
  std::array<double, 3> center;  //! Global center of lattice
  std::array<double, 2> pitch;   //! Lattice tile width and height
};

} //  namespace openmc
#endif // LATTICE_H
