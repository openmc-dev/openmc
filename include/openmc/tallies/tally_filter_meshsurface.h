#ifndef OPENMC_TALLY_FILTER_MESHSURFACE_H
#define OPENMC_TALLY_FILTER_MESHSURFACE_H

#include <cstdint>
#include <sstream>

#include "openmc/error.h"
#include "openmc/mesh.h"
#include "openmc/tallies/tally_filter.h"


namespace openmc {

class MeshSurfaceFilter : public TallyFilter
{
public:
  virtual ~MeshSurfaceFilter() override = default;

  virtual void
  from_xml(pugi::xml_node node) override
  {
    auto bins_ = get_node_array<int32_t>(node, "bins");
    if (bins_.size() != 1)
      fatal_error("Only one mesh can be specified per mesh filter.");

    auto id = bins_[0];
    auto search = mesh_map.find(id);
    if (search != mesh_map.end()) {
      mesh_ = search->second;
    } else{
      std::stringstream err_msg;
      err_msg << "Could not find cell " << id << " specified on tally filter.";
      fatal_error(err_msg);
    }
  }

  virtual void
  get_all_bins(Particle* p, int estimator, TallyFilterMatch& match)
  const override
  {
    meshes[mesh_]->surface_bins_crossed(p, match.bins);
    for (auto b : match.bins) match.weights.push_back(1.0);
  }

protected:
  int32_t mesh_;
};

} // namespace openmc
#endif // OPENMC_TALLY_FILTER_MESHSURFACE_H
