//! \file hex_mesh.h
//! \brief Mesh types used for tallies, Shannon entropy, CMFD, etc.

#ifndef OPENMC_HEX_MESH_H
#define OPENMC_HEX_MESH_H

#include <cmath>
#include <unordered_map>

#include "hdf5.h"
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"
#include <gsl/gsl-lite.hpp>

#include "openmc/error.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/mesh.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/vector.h"
#include "openmc/xml_interface.h"

#ifdef DAGMC
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Core.hpp"
#include "moab/GeomUtil.hpp"
#include "moab/Matrix3.hpp"
#endif

#ifdef LIBMESH
#include "libmesh/bounding_box.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/explicit_system.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/point.h"
#endif

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

//==============================================================================
// Global variables
//==============================================================================

extern "C" const bool LIBMESH_ENABLED;

class HexagonalMesh : public PeriodicStructuredMesh {
public:
  // Constructors
  HexagonalMesh() = default;
  HexagonalMesh(pugi::xml_node node);

  using HexMeshIndex = std::array<int, 4>;

  using SpiralHexMeshIndex = std::array<int, 3>;

  struct HexMeshDistance {
    HexMeshDistance() = default;
    HexMeshDistance(int _index, bool _max_surface, double _distance)
      : next_index {_index}, max_surface {_max_surface}, distance {_distance}
    {}
    HexMeshIndex next_index {0, 0, 0, 0};
    bool max_surface {true};
    double distance {INFTY};
    bool operator<(const HexMeshDistance& o) const
    {
      return distance < o.distance;
    }
  };

  HexMeshIndex get_hexindices_from_bin(const int32_t) const;

  int32_t get_bin_from_hexindices(const HexMeshIndex& ijkl) const;

  int n_bins() const override;

  int n_surface_bins() const override;

  xt::xtensor<int, 1> get_x_shape() const;

  std::string bin_label(int bin) const;

  int32_t offset_in_ring(const HexMeshIndex& ijkl, int32_t r) const;

  HexMeshIndex rotate_hexindex(const HexMeshIndex& ijkl) const;

  int get_hexindex_in_direction(const Position& r, int i) const;

  int get_index_in_direction(double r, int i) const;

  virtual std::string get_mesh_type() const override;

  static const std::string mesh_type;

  Position sample_element(int32_t bin, uint64_t* seed) const override
  {
    return sample_element(get_hexindices_from_bin(bin), seed);
  };

  Position sample_element(const HexMeshIndex& ijkl, uint64_t* seed) const;

  HexMeshDistance distance_to_grid_boundary(const HexMeshIndex& ijk, int i,
    const Position& r0, const Direction& u, double l) const;

  StructuredMesh::MeshDistance distance_to_grid_boundary(const MeshIndex& ijk,
    int i, const Position& r0, const Direction& u, double l) const override;

  std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const override;

  void to_hdf5_inner(hid_t group) const override;

  Position get_position_from_hexindex(HexMeshIndex ijkl) const;

  HexMeshIndex get_hexindices(Position r, bool& in_mesh) const;

  bool in_hexmesh(HexMeshIndex& ijkl) const;

  double volume(const StructuredMesh::MeshIndex& ijk) const;

  // Data members
  double size_ {0};
  int hex_radius_ {0};
  int hex_count_ {0};
  double r_encl_ {0};

  double element_volume_ {0};
  double volume_frac_ {0};

  xt::xtensor<double, 1> width_;

  template<class T>
  void raytrace_mesh(
    Position r0, Position r1, const Direction& u, T tally) const;

  double find_surface_crossing(const Position& r, const Direction& u) const;

  int find_surface_crossing_index(
    const Position& r, const Direction& u, double l) const;

  double find_r_crossing(const Position& r, const Direction& u, double l) const;

  HexMeshDistance distance_to_hex_boundary(const HexMeshIndex& ijkl, int i,
    const Position& r0, const Direction& u, double l) const;

  Position get_hex_center(HexMeshIndex ijkl) const;

  double volume(const HexMeshIndex& ijkl) const;

  void bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins, vector<double>& lengths) const override;

  void surface_bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins) const override;

private:
  int init_plane_normals();

  int scale_grid_vectors(double s);
  // needed data members:
  // size_ length of basis vectors
  // hex-basis vectors r_,q_,s_
  // surface normal vectors n0_, n1_, n2_
  // max-index of hex-rings meaning indices run from -max,...,0,...,max
  // - this means:
  // the width of a hexagon is sqrt(3)*size_

  Position r_ {0, -1, 0};
  Position q_ {sqrt(3.0) * 0.5, 0.5, 0};
  Position s_ {-sqrt(3.0) * 0.5, 0.5, 0};

  Position n0_ {0, 0, 0};
  Position n1_ {0, 0, 0};
  Position n2_ {0, 0, 0};

  int hex_radius(const HexMeshIndex& ijkl) const;

  int hex_distance(const HexMeshIndex& ijkl0, const HexMeshIndex& ijkl1) const;
};

} // namespace openmc
#endif // OPENMC_HEX_MESH_H
