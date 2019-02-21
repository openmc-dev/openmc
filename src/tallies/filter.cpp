#include "openmc/tallies/filter.h"

#include <algorithm> // for max
#include <cstring> // for strcpy
#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"  // for MAX_LINE_LEN;
#include "openmc/error.h"
#include "openmc/tallies/filter_azimuthal.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/filter_cellborn.h"
#include "openmc/tallies/filter_cellfrom.h"
#include "openmc/tallies/filter_delayedgroup.h"
#include "openmc/tallies/filter_distribcell.h"
#include "openmc/tallies/filter_energyfunc.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_legendre.h"
#include "openmc/tallies/filter_material.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/filter_meshsurface.h"
#include "openmc/tallies/filter_mu.h"
#include "openmc/tallies/filter_particle.h"
#include "openmc/tallies/filter_polar.h"
#include "openmc/tallies/filter_sph_harm.h"
#include "openmc/tallies/filter_sptl_legendre.h"
#include "openmc/tallies/filter_surface.h"
#include "openmc/tallies/filter_universe.h"
#include "openmc/tallies/filter_zernike.h"

// explicit template instantiation definition
template class std::vector<openmc::FilterMatch>;

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {
  std::vector<FilterMatch> filter_matches;
}

namespace model {
  std::vector<std::unique_ptr<Filter>> tally_filters;
  std::unordered_map<int, int> filter_map;
}

//==============================================================================
// Non-member functions
//==============================================================================

Filter*
allocate_filter(const std::string& type)
{
  if (type == "azimuthal") {
    model::tally_filters.push_back(std::make_unique<AzimuthalFilter>());
  } else if (type == "cell") {
    model::tally_filters.push_back(std::make_unique<CellFilter>());
  } else if (type == "cellborn") {
    model::tally_filters.push_back(std::make_unique<CellbornFilter>());
  } else if (type == "cellfrom") {
    model::tally_filters.push_back(std::make_unique<CellFromFilter>());
  } else if (type == "distribcell") {
    model::tally_filters.push_back(std::make_unique<DistribcellFilter>());
  } else if (type == "delayedgroup") {
    model::tally_filters.push_back(std::make_unique<DelayedGroupFilter>());
  } else if (type == "energyfunction") {
    model::tally_filters.push_back(std::make_unique<EnergyFunctionFilter>());
  } else if (type == "energy") {
    model::tally_filters.push_back(std::make_unique<EnergyFilter>());
  } else if (type == "energyout") {
    model::tally_filters.push_back(std::make_unique<EnergyoutFilter>());
  } else if (type == "legendre") {
    model::tally_filters.push_back(std::make_unique<LegendreFilter>());
  } else if (type == "material") {
    model::tally_filters.push_back(std::make_unique<MaterialFilter>());
  } else if (type == "mesh") {
    model::tally_filters.push_back(std::make_unique<MeshFilter>());
  } else if (type == "meshsurface") {
    model::tally_filters.push_back(std::make_unique<MeshSurfaceFilter>());
  } else if (type == "mu") {
    model::tally_filters.push_back(std::make_unique<MuFilter>());
  } else if (type == "particle") {
    model::tally_filters.push_back(std::make_unique<ParticleFilter>());
  } else if (type == "polar") {
    model::tally_filters.push_back(std::make_unique<PolarFilter>());
  } else if (type == "surface") {
    model::tally_filters.push_back(std::make_unique<SurfaceFilter>());
  } else if (type == "spatiallegendre") {
    model::tally_filters.push_back(std::make_unique<SpatialLegendreFilter>());
  } else if (type == "sphericalharmonics") {
    model::tally_filters.push_back(std::make_unique<SphericalHarmonicsFilter>());
  } else if (type == "universe") {
    model::tally_filters.push_back(std::make_unique<UniverseFilter>());
  } else if (type == "zernike") {
    model::tally_filters.push_back(std::make_unique<ZernikeFilter>());
  } else if (type == "zernikeradial") {
    model::tally_filters.push_back(std::make_unique<ZernikeRadialFilter>());
  } else {
    throw std::runtime_error{"Unknown filter type: " + type};
  }
  return model::tally_filters.back().get();
}

extern "C" size_t tally_filters_size()
{
  return model::tally_filters.size();
}

//==============================================================================
// C API functions
//==============================================================================

int verify_filter(int32_t index)
{
  if (index < 1 || index > model::tally_filters.size()) {
    set_errmsg("Filter index is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

extern "C" int
openmc_filter_get_id(int32_t index, int32_t* id)
{
  int err = verify_filter(index);
  if (err) return err;

  // TODO: off-by-one
  *id = model::tally_filters[index-1]->id_;
  return 0;
}

extern "C" int
openmc_filter_set_id(int32_t index, int32_t id)
{
  int err = verify_filter(index);
  if (err) return err;

  if (model::filter_map.find(id) != model::filter_map.end()) {
    set_errmsg("Two filters have the same ID: " + std::to_string(id));
    return OPENMC_E_INVALID_ID;
  }

  // TODO: off-by-one
  model::tally_filters[index-1]->id_ = id;
  model::filter_map[id] = index - 1;
  return 0;
}

extern "C" int
openmc_filter_get_type(int32_t index, char* type)
{
  int err = verify_filter(index);
  if (err) return err;

  // TODO: off-by-one
  std::strcpy(type, model::tally_filters[index-1]->type().c_str());
  return 0;
}

extern "C" int
openmc_get_filter_index(int32_t id, int32_t* index)
{
  auto it = model::filter_map.find(id);
  if (it == model::filter_map.end()) {
    set_errmsg("No filter exists with ID=" + std::to_string(id) + ".");
    return OPENMC_E_INVALID_ID;
  }

  // TODO: off-by-one
  *index = it->second + 1;
  return 0;
}

extern "C" void
openmc_get_filter_next_id(int32_t* id)
{
  int32_t largest_filter_id = 0;
  for (const auto& t : model::tally_filters) {
    largest_filter_id = std::max(largest_filter_id, t->id_);
  }
  *id = largest_filter_id + 1;
}

extern "C" int
openmc_new_filter(const char* type, int32_t* index)
{
  allocate_filter(type);
  // TODO: off-by-one
  *index = model::tally_filters.size();
  return 0;
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  // filter_match_point moved to simulation.cpp


  int32_t filter_get_id(Filter* filt) {return filt->id_;}

  void filter_set_id(Filter* filt, int32_t id) {filt->id_ = id;}

  void filter_from_xml(Filter* filt, pugi::xml_node* node)
  {filt->from_xml(*node);}

  void
  filter_get_all_bins(Filter* filt, Particle* p, int estimator,
                      FilterMatch* match)
  {
    filt->get_all_bins(p, estimator, *match);
  }

  void filter_to_statepoint(Filter* filt, hid_t group)
  {filt->to_statepoint(group);}

  void filter_text_label(Filter* filt, int bin, char* label)
  {
    std::string label_str = filt->text_label(bin);
    int i = 0;
    for (; i < label_str.size() && i < MAX_LINE_LEN; i++)
      label[i] = label_str[i];
    label[i] = '\0';
  }

  void filter_initialize(Filter* filt) {filt->initialize();}

  int filter_n_bins(Filter* filt) {return filt->n_bins_;}

  int mesh_filter_get_mesh(MeshFilter* filt) {return filt->mesh();}

  int sphharm_filter_get_cosine(SphericalHarmonicsFilter* filt)
  {return static_cast<int>(filt->cosine_);}
}

} // namespace openmc
