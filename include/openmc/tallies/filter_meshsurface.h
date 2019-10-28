#ifndef OPENMC_TALLIES_FILTER_MESHSURFACE_H
#define OPENMC_TALLIES_FILTER_MESHSURFACE_H

#include "openmc/tallies/filter_mesh.h"

namespace openmc {

class MeshSurfaceFilter : public MeshFilter
{
public:
  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "meshsurface";}

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  void set_mesh(int32_t mesh) override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MESHSURFACE_H
