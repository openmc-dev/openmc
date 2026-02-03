#ifndef OPENMC_XDG_H
#define OPENMC_XDG_H

namespace openmc {
extern "C" const bool XDG_ENABLED;
}

// always include the XML interface header
#include "openmc/xml_interface.h"

#ifdef OPENMC_XDG_ENABLED

#include "xdg/xdg.h"

#include "openmc/mesh.h"
#include "openmc/position.h"
namespace openmc {

class XDGMesh : public UnstructuredMesh{

public:
  // Constructors
  XDGMesh() = default;
  XDGMesh(pugi::xml_node node);
  XDGMesh(hid_t group);
  XDGMesh(const std::string& filename, double length_multiplier = 1.0);
  XDGMesh(std::shared_ptr<xdg::XDG> external_xdg);

  static const std::string mesh_lib_type;

  const std::shared_ptr<xdg::XDG>& xdg_instance() const { return xdg_; }

  // Overridden Methods

  //! Perform any preparation needed to support use in mesh filters
  void prepare_for_point_location() override;

  Position sample_element(int32_t bin, uint64_t* seed) const override;

  void bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins, vector<double>& lengths) const override;

  int get_bin(Position r) const override;

  bool bin_is_valid(int bin) const
  {
    return bin >= 0 && bin < n_bins();
  }

  xdg::MeshID bin_to_mesh_id(int bin) const;

  int mesh_id_to_bin(xdg::MeshID id) const;

  int n_bins() const override;

  int n_surface_bins() const override;

  std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const override;

  std::string library() const override;

  std::string mesh_library() const;

  //! Add a score to the mesh instance
  void add_score(const std::string& score) override {};

  //! Remove all scores from the mesh instance
  void remove_scores() override {};

  //! Set data for a score
  void set_score_data(const std::string& score, const vector<double>& values,
    const vector<double>& std_dev) override {};

  //! Write the mesh with any current tally data
  void write(const std::string& base_filename) const override;

  Position centroid(int bin) const override;

  int n_vertices() const override;

  Position vertex(int id) const override;

  std::vector<int> connectivity(int id) const override;

  //! Get the volume of a mesh bin
  //
  //! \param[in] bin Bin to return the volume for
  //! \return Volume of the bin
  double volume(int bin) const override;

  //! Get the distance to the nearest boundary for a given position and direction
  //! \param[in] g GeometryState object containing position and direction
  //! \return NextMeshCell struct containing distance, face index, and next indices
  NextMeshCell distance_to_bin_boundary(GeometryState& g) const;

  //! Get the distance to the nearest boundary for a given position and direction
  //! \param[in] r Position to check
  //! \param[in] u Direction to check
  //! \return Distance to the nearest boundary
  NextMeshCell distance_to_bin_boundary(int bin, const Position& r, const Direction& u) const;

private:
  void initialize() override;

  std::shared_ptr<xdg::XDG> xdg_; //!< XDG instance
  xdg::MeshLibrary mesh_library_ {xdg::MeshLibrary::LIBMESH}; //!< Mesh library type
};

} // namespace openmc

#endif // OPENMC_XDG_ENABLED

#endif // OPENMC_XDG_H