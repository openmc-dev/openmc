#ifndef LATTICE_H
#define LATTICE_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "hdf5.h"
#include "pugixml/pugixml.hpp"


namespace openmc {

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
  int32_t id;                     //! Universe ID number
  std::string name;               //! User-defined name
  //std::vector<double> pitch;      //! Pitch along each basis
  //std::vector<int32_t> universes; //! Universes filling each lattice tile
  //int32_t outer;                  //! Universe tiled outside the lattice
  //std::vector<int32_t> offset;    //! Distribcell offsets

  explicit Lattice(pugi::xml_node lat_node);

  virtual ~Lattice() {}

  //virtual bool are_valid_indices(const int i_xyz[3]) const = 0;

  //virtual void get_indices(const double global_xyz[3], int i_xyz[3]) const = 0;

  //virtual void get_local_xyz(const double global_xyz[3], const int i_xyz[3],
  //                           double local_xyz[3]) const = 0;

  //! Write all information needed to reconstruct the lattice to an HDF5 group.
  //! @param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;
};

//==============================================================================
//==============================================================================

class RectLattice : public Lattice
{
public:
  explicit RectLattice(pugi::xml_node lat_node);

  virtual ~RectLattice() {}
};

class HexLattice : public Lattice
{
public:
  explicit HexLattice(pugi::xml_node lat_node);

  virtual ~HexLattice() {}
};

} //  namespace openmc
#endif // LATTICE_H
