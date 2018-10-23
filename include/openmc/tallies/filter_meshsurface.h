#ifndef OPENMC_TALLIES_FILTER_MESHSURFACE_H
#define OPENMC_TALLIES_FILTER_MESHSURFACE_H

#include "openmc/tallies/filter_mesh.h"

namespace openmc {

class MeshSurfaceFilter : public MeshFilter
{
public:
  std::string type() const override {return "meshsurface";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(Particle* p, int estimator, FilterMatch& match)
  const override;

  std::string text_label(int bin) const override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MESHSURFACE_H
