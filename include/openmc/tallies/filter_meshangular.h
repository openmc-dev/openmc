#ifndef OPENMC_TALLIES_FILTER_MESHSURFACE_H
#define OPENMC_TALLIES_FILTER_MESHSURFACE_H

#include "openmc/tallies/filter_mesh.h"

namespace openmc {

//==============================================================================
//! Indexes the direction of particle events to a mesh.
//==============================================================================

class MeshAngularFilter : public MeshFilter {
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "meshangular"; }
  FilterType type() const override { return FilterType::MESH_ANGULAR; }

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  //----------------------------------------------------------------------------
  // Accessors

  void set_translation(const Position& translation) const override
  {
    fatal_error("Angular mesh filters do not permit translation.");
  }

  void set_translation(const double translation[3]) const override
  {
    fatal_error("Angular mesh filters do not permit translation.");
  }
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MESHSURFACE_H
