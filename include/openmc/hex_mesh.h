//! \file hex_mesh.h
//! \brief Mesh types used for tallies, Shannon entropy, CMFD, etc.

#ifndef OPENMC_HEX_MESH_H
#define OPENMC_HEX_MESH_H

#include <unordered_map>
#include <cmath>

#include "hdf5.h"
#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"
#include <gsl/gsl-lite.hpp>

#include "openmc/error.h"
#include "openmc/memory.h" // for unique_ptr
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

/*
class CylindricalMesh : public PeriodicStructuredMesh {
public:
  // Constructors
  CylindricalMesh() = default;
  CylindricalMesh(pugi::xml_node node);

  // Overridden methods
  virtual MeshIndex get_indices(Position r, bool& in_mesh) const override;

  int get_index_in_direction(double r, int i) const override;

  virtual std::string get_mesh_type() const override;

  static const std::string mesh_type;

  Position sample_element(const MeshIndex& ijk, uint64_t* seed) const override;

  MeshDistance distance_to_grid_boundary(const MeshIndex& ijk, int i,
    const Position& r0, const Direction& u, double l) const override;

  std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const override;

  void to_hdf5_inner(hid_t group) const override;

  double volume(const MeshIndex& ijk) const override;

  // grid accessors
  double r(int i) const { return grid_[0][i]; }
  double phi(int i) const { return grid_[1][i]; }
  double z(int i) const { return grid_[2][i]; }

  int set_grid();

  // Data members
  array<vector<double>, 3> grid_;

private:
  double find_r_crossing(
    const Position& r, const Direction& u, double l, int shell) const;
  double find_phi_crossing(
    const Position& r, const Direction& u, double l, int shell) const;
  StructuredMesh::MeshDistance find_z_crossing(
    const Position& r, const Direction& u, double l, int shell) const;

  bool full_phi_ {false};

  inline int sanitize_angular_index(int idx, bool full, int N) const
  {
    if ((idx > 0) and (idx <= N)) {
      return idx;
    } else if (full) {
      return (idx + N - 1) % N + 1;
    } else {
      return 0;
    }
  }

  inline int sanitize_phi(int idx) const
  {
    return sanitize_angular_index(idx, full_phi_, shape_[1]);
  }
};
*/

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

enum class ElementType { UNSUPPORTED = -1, LINEAR_TET, LINEAR_HEX };

//==============================================================================
// Global variables
//==============================================================================

extern "C" const bool LIBMESH_ENABLED;

class HexagonalMesh : public PeriodicStructuredMesh {
public:
  //Constructors
  HegagonalMesh() = default;
  HexagonalMesh(pugi::xml_node node);

  //Overridden methods
  //directly copied in from cylmesh from hereon in - go through
  int get_index_in_direction(double r, int i) const override;

  virtual std::string get_mesh_type() const override;

  static const std::string mesh_type;

  Position sample_element(const MeshIndex& ijk, uint64_t* seed) const override;

  MeshDistance distance_to_grid_boundary(const MeshIndex& ijk, int i,
    const Position& r0, const Direction& u, double l) const override;

  std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const override;

  void to_hdf5(hid_t group) const override;

  double volume(const MeshIndex& ijk) const override;

  // grid accessors
  double a(int i) const { return grid_[0][i]; }
  double b(int i) const { return grid_[1][i]; }
  double z(int i) const { return grid_[2][i]; }

  int set_grid();

  // Data members
  array<vector<double>, 3> grid_;

private:
  double find surface_crossing(
     const Position& r, const Direction& u) const;

  int find_surface_crossing_index(const Position& r, const Direction& u) const;

  double abs_a_ {0};
};
