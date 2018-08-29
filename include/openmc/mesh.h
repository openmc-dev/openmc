#ifndef OPENMC_MESH_H
#define OPENMC_MESH_H

#include <cstdint> // for size_t
#include <memory> // for unique_ptr
#include <vector>
#include <unordered_map>

#include "hdf5.h"
#include "pugixml.hpp"
#include "xtensor/xarray.hpp"

#include "openmc/position.h"

namespace openmc {

//==============================================================================
//! Tessellation of n-dimensional Euclidean space by congruent squares or cubes
//==============================================================================

class RegularMesh {
public:
  // Constructors
  RegularMesh(pugi::xml_node node);

  // Methods
  int get_bin(Position r);
  int get_bin_from_indices(const int* ijk);
  void get_indices(Position r, int* ijk, bool* in_mesh);
  void get_indices_from_bin(int bin, int* ijk);
  bool intersects(Position r0, Position r1);
  void to_hdf5(hid_t group);

  int id_ {-1};  //!< User-specified ID
  int n_dimension_;
  double volume_frac_;
  xt::xarray<int> shape_;
  xt::xarray<double> lower_left_;
  xt::xarray<double> upper_right_;
  xt::xarray<double> width_;

private:
  bool intersects_1d(Position r0, Position r1);
  bool intersects_2d(Position r0, Position r1);
  bool intersects_3d(Position r0, Position r1);
};

//==============================================================================
// Global variables
//==============================================================================

extern std::vector<std::unique_ptr<RegularMesh>> meshes;

extern std::unordered_map<int32_t, int32_t> mesh_map;


} // namespace openmc

#endif // OPENMC_MESH_H
