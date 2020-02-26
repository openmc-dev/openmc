#include "openmc/tallies/filter.h"

#include <algorithm> // for max
#include <cstring> // for strcpy
#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"  // for MAX_LINE_LEN;
#include "openmc/error.h"
#include "openmc/xml_interface.h"
#include "openmc/tallies/filter_azimuthal.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/filter_cellborn.h"
#include "openmc/tallies/filter_cellfrom.h"
#include "openmc/tallies/filter_cell_instance.h"
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

namespace model {
  std::unordered_map<int, int> filter_map;
  std::vector<std::unique_ptr<Filter>> tally_filters;
}

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
    model::tally_filters.push_back(std::make_unique<AzimuthalFilter>());
  } else if (type == "cell") {
    model::tally_filters.push_back(std::make_unique<CellFilter>());
  } else if (type == "cellborn") {
    model::tally_filters.push_back(std::make_unique<CellbornFilter>());
  } else if (type == "cellfrom") {
    model::tally_filters.push_back(std::make_unique<CellFromFilter>());
  } else if (type == "cellinstance") {
    model::tally_filters.push_back(std::make_unique<CellInstanceFilter>());
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

  // Assign ID
  model::tally_filters.back()->set_id(id);

  return model::tally_filters.back().get();
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
    throw std::runtime_error{"Two filters have the same ID: " + std::to_string(id)};
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

extern "C" int
openmc_filter_get_id(int32_t index, int32_t* id)
{
  if (int err = verify_filter(index)) return err;

  *id = model::tally_filters[index]->id();
  return 0;
}

extern "C" int
openmc_filter_set_id(int32_t index, int32_t id)
{
  if (int err = verify_filter(index)) return err;

  model::tally_filters[index]->set_id(id);
  return 0;
}

extern "C" int
openmc_filter_get_type(int32_t index, char* type)
{
  if (int err = verify_filter(index)) return err;

  std::strcpy(type, model::tally_filters[index]->type().c_str());
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

  *index = it->second;
  return 0;
}

extern "C" void
openmc_get_filter_next_id(int32_t* id)
{
  int32_t largest_filter_id = 0;
  for (const auto& t : model::tally_filters) {
    largest_filter_id = std::max(largest_filter_id, t->id());
  }
  *id = largest_filter_id + 1;
}

extern "C" int
openmc_new_filter(const char* type, int32_t* index)
{
  *index = model::tally_filters.size();
  Filter::create(type);
  return 0;
}

} // namespace openmc
