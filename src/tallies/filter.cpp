#include "openmc/tallies/filter.h"

#include <algorithm> // for max
#include <cstring>   // for strcpy
#include <string>

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/constants.h" // for MAX_LINE_LEN;
#include "openmc/error.h"
#include "openmc/tallies/filter_azimuthal.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/filter_cell_instance.h"
#include "openmc/tallies/filter_cellborn.h"
#include "openmc/tallies/filter_cellfrom.h"
#include "openmc/tallies/filter_collision.h"
#include "openmc/tallies/filter_delayedgroup.h"
#include "openmc/tallies/filter_distribcell.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_energyfunc.h"
#include "openmc/tallies/filter_legendre.h"
#include "openmc/tallies/filter_material.h"
#include "openmc/tallies/filter_materialfrom.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/filter_meshsurface.h"
#include "openmc/tallies/filter_mu.h"
#include "openmc/tallies/filter_particle.h"
#include "openmc/tallies/filter_polar.h"
#include "openmc/tallies/filter_sph_harm.h"
#include "openmc/tallies/filter_sptl_legendre.h"
#include "openmc/tallies/filter_surface.h"
#include "openmc/tallies/filter_time.h"
#include "openmc/tallies/filter_universe.h"
#include "openmc/tallies/filter_zernike.h"
#include "openmc/xml_interface.h"

// explicit template instantiation definition
template class openmc::vector<openmc::FilterMatch>;

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
std::unordered_map<int, int> filter_map;
vector<unique_ptr<Filter>> tally_filters;
} // namespace model

//==============================================================================
// Non-member functions
//==============================================================================

extern "C" size_t tally_filters_size()
{
  return model::tally_filters.size();
}

//==============================================================================
// Filter implementation
//==============================================================================

Filter::Filter()
{
  index_ = model::tally_filters.size(); // Avoids warning about narrowing
}

Filter::~Filter()
{
  model::filter_map.erase(id_);
}

Filter* Filter::create(pugi::xml_node node)
{
  // Copy filter id
  if (!check_for_node(node, "id")) {
    fatal_error("Must specify id for filter in tally XML file.");
  }
  int filter_id = std::stoi(get_node_value(node, "id"));

  // Convert filter type to lower case
  std::string s;
  if (check_for_node(node, "type")) {
    s = get_node_value(node, "type", true);
  }

  // Allocate according to the filter type
  auto f = Filter::create(s, filter_id);

  // Read filter data from XML
  f->from_xml(node);
  return f;
}

Filter* Filter::create(const std::string& type, int32_t id)
{
  if (type == "azimuthal") {
    return Filter::create<AzimuthalFilter>(id);
  } else if (type == "cell") {
    return Filter::create<CellFilter>(id);
  } else if (type == "cellborn") {
    return Filter::create<CellBornFilter>(id);
  } else if (type == "cellfrom") {
    return Filter::create<CellFromFilter>(id);
  } else if (type == "cellinstance") {
    return Filter::create<CellInstanceFilter>(id);
  } else if (type == "distribcell") {
    return Filter::create<DistribcellFilter>(id);
  } else if (type == "delayedgroup") {
    return Filter::create<DelayedGroupFilter>(id);
  } else if (type == "energyfunction") {
    return Filter::create<EnergyFunctionFilter>(id);
  } else if (type == "energy") {
    return Filter::create<EnergyFilter>(id);
  } else if (type == "collision") {
    return Filter::create<CollisionFilter>(id);
  } else if (type == "energyout") {
    return Filter::create<EnergyoutFilter>(id);
  } else if (type == "legendre") {
    return Filter::create<LegendreFilter>(id);
  } else if (type == "material") {
    return Filter::create<MaterialFilter>(id);
  } else if (type == "materialfrom") {
    return Filter::create<MaterialFromFilter>(id);
  } else if (type == "mesh") {
    return Filter::create<MeshFilter>(id);
  } else if (type == "meshsurface") {
    return Filter::create<MeshSurfaceFilter>(id);
  } else if (type == "mu") {
    return Filter::create<MuFilter>(id);
  } else if (type == "particle") {
    return Filter::create<ParticleFilter>(id);
  } else if (type == "polar") {
    return Filter::create<PolarFilter>(id);
  } else if (type == "surface") {
    return Filter::create<SurfaceFilter>(id);
  } else if (type == "spatiallegendre") {
    return Filter::create<SpatialLegendreFilter>(id);
  } else if (type == "sphericalharmonics") {
    return Filter::create<SphericalHarmonicsFilter>(id);
  } else if (type == "time") {
    return Filter::create<TimeFilter>(id);
  } else if (type == "universe") {
    return Filter::create<UniverseFilter>(id);
  } else if (type == "zernike") {
    return Filter::create<ZernikeFilter>(id);
  } else if (type == "zernikeradial") {
    return Filter::create<ZernikeRadialFilter>(id);
  } else {
    throw std::runtime_error {fmt::format("Unknown filter type: {}", type)};
  }
  return nullptr;
}

void Filter::set_id(int32_t id)
{
  Expects(id >= 0 || id == C_NONE);

  // Clear entry in filter map if an ID was already assigned before
  if (id_ != C_NONE) {
    model::filter_map.erase(id_);
    id_ = C_NONE;
  }

  // Make sure no other filter has same ID
  if (model::filter_map.find(id) != model::filter_map.end()) {
    throw std::runtime_error {
      "Two filters have the same ID: " + std::to_string(id)};
  }

  // If no ID specified, auto-assign next ID in sequence
  if (id == C_NONE) {
    id = 0;
    for (const auto& f : model::tally_filters) {
      id = std::max(id, f->id_);
    }
    ++id;
  }

  // Update ID and entry in filter map
  id_ = id;
  model::filter_map[id] = index_;
}

//==============================================================================
// C API functions
//==============================================================================

int verify_filter(int32_t index)
{
  if (index < 0 || index >= model::tally_filters.size()) {
    set_errmsg("Filter index is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

extern "C" int openmc_filter_get_id(int32_t index, int32_t* id)
{
  if (int err = verify_filter(index))
    return err;

  *id = model::tally_filters[index]->id();
  return 0;
}

extern "C" int openmc_filter_set_id(int32_t index, int32_t id)
{
  if (int err = verify_filter(index))
    return err;

  model::tally_filters[index]->set_id(id);
  return 0;
}

extern "C" int openmc_filter_get_type(int32_t index, char* type)
{
  if (int err = verify_filter(index))
    return err;

  std::strcpy(type, model::tally_filters[index]->type_str().c_str());
  return 0;
}

extern "C" int openmc_filter_get_num_bins(int32_t index, int* n_bins)
{
  if (int err = verify_filter(index))
    return err;

  *n_bins = model::tally_filters[index]->n_bins();
  return 0;
}

extern "C" int openmc_get_filter_index(int32_t id, int32_t* index)
{
  auto it = model::filter_map.find(id);
  if (it == model::filter_map.end()) {
    set_errmsg("No filter exists with ID=" + std::to_string(id) + ".");
    return OPENMC_E_INVALID_ID;
  }

  *index = it->second;
  return 0;
}

extern "C" void openmc_get_filter_next_id(int32_t* id)
{
  int32_t largest_filter_id = 0;
  for (const auto& t : model::tally_filters) {
    largest_filter_id = std::max(largest_filter_id, t->id());
  }
  *id = largest_filter_id + 1;
}

extern "C" int openmc_new_filter(const char* type, int32_t* index)
{
  *index = model::tally_filters.size();
  Filter::create(type);
  return 0;
}

} // namespace openmc
