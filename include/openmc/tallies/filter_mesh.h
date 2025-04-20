#ifndef OPENMC_TALLIES_FILTER_MESH_H
#define OPENMC_TALLIES_FILTER_MESH_H

#include <cstdint>

#include "openmc/position.h"
#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Indexes the location of particle events to a mesh.  For tracklength tallies,
//! it will produce multiple valid bins and the bin weight will correspond to
//! the fraction of the track length that lies in that bin.
//==============================================================================

class MeshFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~MeshFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "mesh"; }
  FilterType type() const override { return FilterType::MESH; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  virtual int32_t mesh() const { return mesh_; }

  virtual void set_mesh(int32_t mesh);

  virtual void set_translation(const Position& translation);

  virtual void set_translation(const double translation[3]);

  virtual const Position& translation() const { return translation_; }

  virtual bool translated() const { return translated_; }

protected:
  //----------------------------------------------------------------------------
  // Data members

  int32_t mesh_;            //!< Index of the mesh
  bool translated_ {false}; //!< Whether or not the filter is translated
  Position translation_ {0.0, 0.0, 0.0}; //!< Filter translation
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MESH_H
