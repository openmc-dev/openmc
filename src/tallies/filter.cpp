#include "openmc/tallies/filter.h"

#include <algorithm> // for max
#include <cstring> // for strcpy
#include <string>

#include <fmt/core.h>

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
  Filter* tally_filters;
  int32_t n_tally_filters {0};
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

FilterType get_filter_type(const std::string& type)
{
  if (type == "azimuthal") {
    return FilterType::AzimuthalFilter;
  } else if (type == "cell") {
    return FilterType::CellFilter;
  } else if (type == "cellborn") {
    return FilterType::CellbornFilter;
  } else if (type == "cellfrom") {
    return FilterType::CellFromFilter;
  } else if (type == "cellinstance") {
    return FilterType::CellInstanceFilter;
  } else if (type == "distribcell") {
    return FilterType::DistribcellFilter;
  } else if (type == "delayedgroup") {
    return FilterType::DelayedGroupFilter;
  } else if (type == "energyfunction") {
    return FilterType::EnergyFunctionFilter;
  } else if (type == "energy") {
    return FilterType::EnergyFilter;
  } else if (type == "energyout") {
    return FilterType::EnergyoutFilter;
  } else if (type == "legendre") {
    return FilterType::LegendreFilter;
  } else if (type == "material") {
    return FilterType::MaterialFilter;
  } else if (type == "mesh") {
    return FilterType::MeshFilter;
  } else if (type == "meshsurface") {
    return FilterType::MeshSurfaceFilter;
  } else if (type == "mu") {
    return FilterType::MuFilter;
  } else if (type == "particle") {
    return FilterType::ParticleFilter;
  } else if (type == "polar") {
    return FilterType::PolarFilter;
  } else if (type == "surface") {
    return FilterType::SurfaceFilter;
  } else if (type == "spatiallegendre") {
    return FilterType::SpatialLegendreFilter;
  } else if (type == "sphericalharmonics") {
    return FilterType::SphericalHarmonicsFilter;
  } else if (type == "universe") {
    return FilterType::UniverseFilter;
  } else if (type == "zernike") {
    return FilterType::ZernikeFilter;
  } else if (type == "zernikeradial") {
    return FilterType::ZernikeRadialFilter;
  } else {
    throw std::runtime_error{fmt::format("Unknown filter type: {}", type)};
    return FilterType::AzimuthalFilter;
  }
}

Filter::Filter(pugi::xml_node node)
{
  index_ = model::tally_filters.size(); // Avoids warning about narrowing

  // Copy filter id
  if (!check_for_node(node, "id")) {
    fatal_error("Must specify id for filter in tally XML file.");
  }
  id_ = std::stoi(get_node_value(node, "id"));
  
  // Convert filter type to lower case
  std::string s;
  if (check_for_node(node, "type")) {
    s = get_node_value(node, "type", true);
  }

  type_ = get_filter_type(s);

  this->from_xml(node);
}

void Filter_from_xml(pugi::xml_node node)
{
  switch(type_){
    case FilterType::AzimuthalFilter          : return AzimuthalFilter_from_xml(node); break;
    case FilterType::CellFilter               : return CellFilter_from_xml(node); break;
    case FilterType::CellInstanceFilter       : return CellInstanceFilter_from_xml(node); break;
    case FilterType::CellbornFilter           : return CellbornFilter_from_xml(node); break;
    case FilterType::CellFromFilter           : return CellFromFilter_from_xml(node); break;
    case FilterType::DelayedGroupFilter       : return DelayedGroupFilter_from_xml(node); break;
    case FilterType::DistribcellFilter        : return DistribcellFilter_from_xml(node); break;
    case FilterType::EnergyFilter             : return EnergyFilter_from_xml(node); break;
    case FilterType::EnergyoutFilter          : return EnergyoutFilter_from_xml(node); break;
    case FilterType::EnergyFunctionFilter     : return EnergyFunctionFilter_from_xml(node); break;
    case FilterType::LegendreFilter           : return LegendreFilter_from_xml(node); break;
    case FilterType::MaterialFilter           : return MaterialFilter_from_xml(node); break;
    case FilterType::MeshFilter               : return MeshFilter_from_xml(node); break;
    case FilterType::MeshSurfaceFilter        : return MeshSurfaceFilter_from_xml(node); break;
    case FilterType::MuFilter                 : return MuFilter_from_xml(node); break;
    case FilterType::ParticleFilter           : return ParticleFilter_from_xml(node); break;
    case FilterType::PolarFilter              : return PolarFilter_from_xml(node); break;
    case FilterType::SphericalHarmonicsFilter : return SphericalHarmonicsFilter_from_xml(node); break;
    case FilterType::SpatialLegendreFilter    : return SpatialLegendreFilter_from_xml(node); break;
    case FilterType::SurfaceFilter            : return SurfaceFilter_from_xml(node); break;
    case FilterType::UniverseFilter           : return UniverseFilter_from_xml(node); break;
    case FilterType::ZernikeFilter            : return ZernikeFilter_from_xml(node); break;
    case FilterType::ZernikeRadialFilter      : return ZernikeRadialFilter_from_xml(node); break;
  }
}

Filter::~Filter()
{
  model::filter_map.erase(id_);
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
  
std::string Filter::type() const 
{
  switch(type_){
    case FilterType::AzimuthalFilter          : return "azimuthal";
    case FilterType::CellFilter               : return "cell";
    case FilterType::CellInstanceFilter       : return "cellborn";
    case FilterType::CellbornFilter           : return "cellfrom";
    case FilterType::CellFromFilter           : return "cellinstance";
    case FilterType::DelayedGroupFilter       : return "distribcell";
    case FilterType::DistribcellFilter        : return "delayedgroup";
    case FilterType::EnergyFilter             : return "energyfunction";
    case FilterType::EnergyoutFilter          : return "energy";
    case FilterType::EnergyFunctionFilter     : return "energyout";
    case FilterType::LegendreFilter           : return "legendre";
    case FilterType::MaterialFilter           : return "material";
    case FilterType::MeshFilter               : return "mesh";
    case FilterType::MeshSurfaceFilter        : return "meshsurface";
    case FilterType::MuFilter                 : return "mu";
    case FilterType::ParticleFilter           : return "particle";
    case FilterType::PolarFilter              : return "polar";
    case FilterType::SphericalHarmonicsFilter : return "surface";
    case FilterType::SpatialLegendreFilter    : return "spatiallegendre";
    case FilterType::SurfaceFilter            : return "sphericalharmonics";
    case FilterType::UniverseFilter           : return "universe";
    case FilterType::ZernikeFilter            : return "zernike";
    case FilterType::ZernikeRadialFilter      : return "zernikeradial";
  }
  return "undefined filter type";
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
