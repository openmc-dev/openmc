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

  void get_all_bins(const Particle* p, TallyEstimator estimator, FilterMatch& match)
  const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  void set_mesh(int32_t mesh) override;

  enum class MeshDir {
    OUT_LEFT,  // x min
    IN_LEFT,  // x min
    OUT_RIGHT,  // x max
    IN_RIGHT,  // x max
    OUT_BACK,  // y min
    IN_BACK,  // y min
    OUT_FRONT,  // y max
    IN_FRONT,  // y max
    OUT_BOTTOM,  // z min
    IN_BOTTOM, // z min
    OUT_TOP, // z max
    IN_TOP // z max
  };
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MESHSURFACE_H
