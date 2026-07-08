#ifndef OPENMC_XDG_H
#define OPENMC_XDG_H

namespace openmc {
extern "C" const bool XDG_ENABLED;
}

#ifdef OPENMC_XDG_ENABLED

#include "xdg/xdg.h"

#include "openmc/mesh.h"
#include "openmc/position.h"
#include "openmc/xml_interface.h"

namespace openmc {

class XDGMesh : public UnstructuredMesh {

public:
  //----------------------------------------------------------------------------
  // Constructors
  XDGMesh() = default;
  XDGMesh(pugi::xml_node node);
  XDGMesh(hid_t group);
  XDGMesh(const std::string& filename, double length_multiplier = 1.0);
  XDGMesh(std::shared_ptr<xdg::XDG> external_xdg);

  static const std::string mesh_type;

  //----------------------------------------------------------------------------
  // Methods

  //! Get the underlying XDG instance
  //!
  //! \return Shared pointer to the XDG instance
  const std::shared_ptr<xdg::XDG>& xdg_instance() const { return xdg_; }

  //! Check whether a bin index is valid
  //!
  //! \param[in] bin Bin index to check
  //! \return True if the bin index is in [0, n_bins())
  bool bin_is_valid(int bin) const { return bin >= 0 && bin < n_bins(); }

  std::string get_mesh_type() const override { return mesh_type; }

  //! Convert a mesh bin index to an XDG MeshID
  //!
  //! \param[in] bin Bin index to convert
  //! \return XDG MeshID corresponding to the bin
  xdg::MeshID bin_to_mesh_id(int bin) const;

  //! Convert an XDG MeshID to a mesh bin index
  //!
  //! \param[in] id XDG MeshID to convert
  //! \return Bin index corresponding to the MeshID, or -1 if invalid
  int mesh_id_to_bin(xdg::MeshID id) const;

  //! Get the name of the underlying mesh library
  //!
  //! \return Name of the mesh library being used in XDG
  std::string library() const override;

  //----------------------------------------------------------------------------
  // Overridden Methods

  //! Perform any preparation needed to support use in mesh filters
  void prepare_for_point_location() override;

  Position sample_element(int32_t bin, uint64_t* seed) const override;

  void bins_crossed(Position r0, Position r1, const Direction& u,
    vector<int>& bins, vector<double>& lengths) const override;

  int get_bin(Position r) const override;

  int n_bins() const override;

  int n_surface_bins() const override;

  std::pair<vector<double>, vector<double>> plot(
    Position plot_ll, Position plot_ur) const override;

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

private:
  void initialize() override;

  //----------------------------------------------------------------------------
  // Private data members
  std::shared_ptr<xdg::XDG> xdg_; //!< XDG instance
  xdg::MeshLibrary mesh_library_ {
    xdg::MeshLibrary::LIBMESH}; //!< Mesh library type
};

} // namespace openmc

#endif // OPENMC_XDG_ENABLED

#endif // OPENMC_XDG_H
