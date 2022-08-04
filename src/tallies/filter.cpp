#include "openmc/tallies/filter.h"

#include <algorithm> // for max
#include <cstring> // for strcpy
#include <string>

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/constants.h"  // for MAX_LINE_LEN;
#include "openmc/error.h"
#include "openmc/xml_interface.h"

// explicit template instantiation definition
template class std::vector<openmc::FilterMatch>;
template class openmc::vector<openmc::FilterMatch>;

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
  return model::n_tally_filters;
}

//==============================================================================
// Filter implementation
//==============================================================================

void Filter::copy_to_device()
{
  map_.copy_to_device();
  imap_.copy_to_device();
  bins_.copy_to_device();
  cells_.copy_to_device();
  cell_instances_.copy_to_device();
  groups_.copy_to_device();
  energy_.copy_to_device();
  y_.copy_to_device();
  materials_.copy_to_device();
  particles_.copy_to_device();
  surfaces_.copy_to_device();
  universes_.copy_to_device();
}

Filter::FilterType Filter::get_filter_type(const std::string& type)
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
  }
}

Filter::Filter(pugi::xml_node node, int32_t index) : index_(index)
{
  // Copy filter id
  if (!check_for_node(node, "id")) {
    fatal_error("Must specify id for filter in tally XML file.");
  }
  int32_t id = std::stoi(get_node_value(node, "id"));

  // Convert filter type to lower case
  std::string s;
  if (check_for_node(node, "type")) {
    s = get_node_value(node, "type", true);
  }

  type_ = get_filter_type(s);

  this->set_id(id);

  this->from_xml(node);
}

void Filter::from_xml(pugi::xml_node node)
{
  switch(type_){
    case FilterType::AzimuthalFilter          : AzimuthalFilter_from_xml(node); break;
    case FilterType::CellFilter               : CellFilter_from_xml(node); break;
    case FilterType::CellInstanceFilter       : CellInstanceFilter_from_xml(node); break;
    case FilterType::CellbornFilter           : CellFilter_from_xml(node); break; // Note - type uses parent
    case FilterType::CellFromFilter           : CellFilter_from_xml(node); break; // Note - type uses parent
    case FilterType::DelayedGroupFilter       : DelayedGroupFilter_from_xml(node); break;
    case FilterType::DistribcellFilter        : DistribcellFilter_from_xml(node); break;
    case FilterType::EnergyFilter             : EnergyFilter_from_xml(node); break;
    case FilterType::EnergyoutFilter          : EnergyFilter_from_xml(node); break; // Note - type uses parent
    case FilterType::EnergyFunctionFilter     : EnergyFunctionFilter_from_xml(node); break;
    case FilterType::LegendreFilter           : LegendreFilter_from_xml(node); break;
    case FilterType::MaterialFilter           : MaterialFilter_from_xml(node); break;
    case FilterType::MeshFilter               : MeshFilter_from_xml(node); break;
    case FilterType::MeshSurfaceFilter        : MeshFilter_from_xml(node); break; // Note -- type uses parent
    case FilterType::MuFilter                 : MuFilter_from_xml(node); break;
    case FilterType::ParticleFilter           : ParticleFilter_from_xml(node); break;
    case FilterType::PolarFilter              : PolarFilter_from_xml(node); break;
    case FilterType::SphericalHarmonicsFilter : SphericalHarmonicsFilter_from_xml(node); break;
    case FilterType::SpatialLegendreFilter    : SpatialLegendreFilter_from_xml(node); break;
    case FilterType::SurfaceFilter            : SurfaceFilter_from_xml(node); break;
    case FilterType::UniverseFilter           : UniverseFilter_from_xml(node); break;
    case FilterType::ZernikeFilter            : ZernikeFilter_from_xml(node); break;
    case FilterType::ZernikeRadialFilter      : ZernikeFilter_from_xml(node); break; // Note -- type uses parent
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
    for (int i = 0; i < model::n_tally_filters; i++) {
      Filter& f = model::tally_filters[i];
      id = std::max(id, f.id_);
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

std::string Filter::text_label(int bin) const
{
  switch(type_){
    case FilterType::AzimuthalFilter          : return AzimuthalFilter_text_label(bin); break;
    case FilterType::CellFilter               : return CellFilter_text_label(bin); break;
    case FilterType::CellInstanceFilter       : return CellInstanceFilter_text_label(bin); break;
    case FilterType::CellbornFilter           : return CellbornFilter_text_label(bin); break;
    case FilterType::CellFromFilter           : return CellFromFilter_text_label(bin); break;
    case FilterType::DelayedGroupFilter       : return DelayedGroupFilter_text_label(bin); break;
    case FilterType::DistribcellFilter        : return DistribcellFilter_text_label(bin); break;
    case FilterType::EnergyFilter             : return EnergyFilter_text_label(bin); break;
    case FilterType::EnergyoutFilter          : return EnergyoutFilter_text_label(bin); break;
    case FilterType::EnergyFunctionFilter     : return EnergyFunctionFilter_text_label(bin); break;
    case FilterType::LegendreFilter           : return LegendreFilter_text_label(bin); break;
    case FilterType::MaterialFilter           : return MaterialFilter_text_label(bin); break;
    case FilterType::MeshFilter               : return MeshFilter_text_label(bin); break;
    case FilterType::MeshSurfaceFilter        : return MeshSurfaceFilter_text_label(bin); break;
    case FilterType::MuFilter                 : return MuFilter_text_label(bin); break;
    case FilterType::ParticleFilter           : return ParticleFilter_text_label(bin); break;
    case FilterType::PolarFilter              : return PolarFilter_text_label(bin); break;
    case FilterType::SphericalHarmonicsFilter : return SphericalHarmonicsFilter_text_label(bin); break;
    case FilterType::SpatialLegendreFilter    : return SpatialLegendreFilter_text_label(bin); break;
    case FilterType::SurfaceFilter            : return SurfaceFilter_text_label(bin); break;
    case FilterType::UniverseFilter           : return UniverseFilter_text_label(bin); break;
    case FilterType::ZernikeFilter            : return ZernikeFilter_text_label(bin); break;
    case FilterType::ZernikeRadialFilter      : return ZernikeRadialFilter_text_label(bin); break;
  }
  return "undefined filter type";
}

void Filter::to_statepoint(hid_t filter_group) const
{
  write_dataset(filter_group, "type", type());
  write_dataset(filter_group, "n_bins", n_bins_);
  switch(type_){
    case FilterType::AzimuthalFilter          : AzimuthalFilter_to_statepoint(filter_group); break;
    case FilterType::CellFilter               : CellFilter_to_statepoint(filter_group); break;
    case FilterType::CellInstanceFilter       : CellInstanceFilter_to_statepoint(filter_group); break;
    case FilterType::CellbornFilter           : CellFilter_to_statepoint(filter_group); break; // Note - uses parent type
    case FilterType::CellFromFilter           : CellFilter_to_statepoint(filter_group); break; // Note - uses parent type
    case FilterType::DelayedGroupFilter       : DelayedGroupFilter_to_statepoint(filter_group); break;
    case FilterType::DistribcellFilter        : DistribcellFilter_to_statepoint(filter_group); break;
    case FilterType::EnergyFilter             : EnergyFilter_to_statepoint(filter_group); break;
    case FilterType::EnergyoutFilter          : EnergyFilter_to_statepoint(filter_group); break; // Note - uses parent type
    case FilterType::EnergyFunctionFilter     : EnergyFunctionFilter_to_statepoint(filter_group); break;
    case FilterType::LegendreFilter           : LegendreFilter_to_statepoint(filter_group); break;
    case FilterType::MaterialFilter           : MaterialFilter_to_statepoint(filter_group); break;
    case FilterType::MeshFilter               : MeshFilter_to_statepoint(filter_group); break;
    case FilterType::MeshSurfaceFilter        : MeshFilter_to_statepoint(filter_group); break; // Note - uses parent type
    case FilterType::MuFilter                 : MuFilter_to_statepoint(filter_group); break;
    case FilterType::ParticleFilter           : ParticleFilter_to_statepoint(filter_group); break;
    case FilterType::PolarFilter              : PolarFilter_to_statepoint(filter_group); break;
    case FilterType::SphericalHarmonicsFilter : SphericalHarmonicsFilter_to_statepoint(filter_group); break;
    case FilterType::SpatialLegendreFilter    : SpatialLegendreFilter_to_statepoint(filter_group); break;
    case FilterType::SurfaceFilter            : SurfaceFilter_to_statepoint(filter_group); break;
    case FilterType::UniverseFilter           : UniverseFilter_to_statepoint(filter_group); break;
    case FilterType::ZernikeFilter            : ZernikeFilter_to_statepoint(filter_group); break;
    case FilterType::ZernikeRadialFilter      : ZernikeFilter_to_statepoint(filter_group); break; // Note - uses parent type
  }
}

void Filter::get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  switch(type_){
    case FilterType::AzimuthalFilter          : AzimuthalFilter_get_all_bins(p, estimator, match); break;
    case FilterType::CellFilter               : CellFilter_get_all_bins(p, estimator, match); break;
    case FilterType::CellInstanceFilter       : CellInstanceFilter_get_all_bins(p, estimator, match); break;
    case FilterType::CellbornFilter           : CellbornFilter_get_all_bins(p, estimator, match); break;
    case FilterType::CellFromFilter           : CellFromFilter_get_all_bins(p, estimator, match); break;
    case FilterType::DelayedGroupFilter       : DelayedGroupFilter_get_all_bins(p, estimator, match); break;
    case FilterType::DistribcellFilter        : DistribcellFilter_get_all_bins(p, estimator, match); break;
    case FilterType::EnergyFilter             : EnergyFilter_get_all_bins(p, estimator, match); break;
    case FilterType::EnergyoutFilter          : EnergyoutFilter_get_all_bins(p, estimator, match); break;
    case FilterType::EnergyFunctionFilter     : EnergyFunctionFilter_get_all_bins(p, estimator, match); break;
    case FilterType::LegendreFilter           : LegendreFilter_get_all_bins(p, estimator, match); break;
    case FilterType::MaterialFilter           : MaterialFilter_get_all_bins(p, estimator, match); break;
    case FilterType::MeshFilter               : MeshFilter_get_all_bins(p, estimator, match); break;
    case FilterType::MeshSurfaceFilter        : MeshSurfaceFilter_get_all_bins(p, estimator, match); break;
    case FilterType::MuFilter                 : MuFilter_get_all_bins(p, estimator, match); break;
    case FilterType::ParticleFilter           : ParticleFilter_get_all_bins(p, estimator, match); break;
    case FilterType::PolarFilter              : PolarFilter_get_all_bins(p, estimator, match); break;
    case FilterType::SphericalHarmonicsFilter : SphericalHarmonicsFilter_get_all_bins(p, estimator, match); break;
    case FilterType::SpatialLegendreFilter    : SpatialLegendreFilter_get_all_bins(p, estimator, match); break;
    case FilterType::SurfaceFilter            : SurfaceFilter_get_all_bins(p, estimator, match); break;
    case FilterType::UniverseFilter           : UniverseFilter_get_all_bins(p, estimator, match); break;
    // TODO: ZernikeFilter relies on calc_zn, which uses a std::vector
    //case FilterType::ZernikeFilter            : ZernikeFilter_get_all_bins(p, estimator, match); break;
    case FilterType::ZernikeRadialFilter      : ZernikeRadialFilter_get_all_bins(p, estimator, match); break;
  }
}

void Filter::set_order(int order)
{
  switch(type_){
    case FilterType::LegendreFilter           : LegendreFilter_set_order(order); break;
    case FilterType::SphericalHarmonicsFilter : SphericalHarmonicsFilter_set_order(order); break;
    case FilterType::SpatialLegendreFilter    : SpatialLegendreFilter_set_order(order); break;
    case FilterType::ZernikeFilter            : ZernikeFilter_set_order(order); break;
    case FilterType::ZernikeRadialFilter      : ZernikeRadialFilter_set_order(order); break;
  }
}

void Filter::set_bins(gsl::span<const double> bins)
{
  switch(type_){
    case FilterType::AzimuthalFilter          :   AzimuthalFilter_set_bins(bins); break;
    case FilterType::EnergyFilter             :   EnergyFilter_set_bins(bins); break;
    case FilterType::EnergyoutFilter          :   EnergyFilter_set_bins(bins); break; // Note we map derived type to parent type
    case FilterType::MuFilter                 :   MuFilter_set_bins(bins); break;
    case FilterType::PolarFilter              :   PolarFilter_set_bins(bins); break;
  }
}

//==============================================================================
// C API functions
//==============================================================================

int verify_filter(int32_t index)
{
  if (index < 0 || index >= model::n_tally_filters) {
    set_errmsg("Filter index is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

extern "C" int
openmc_filter_get_id(int32_t index, int32_t* id)
{
  if (int err = verify_filter(index)) return err;

  *id = model::tally_filters[index].id();
  return 0;
}

extern "C" int
openmc_filter_set_id(int32_t index, int32_t id)
{
  if (int err = verify_filter(index)) return err;

  model::tally_filters[index].set_id(id);
  return 0;
}

extern "C" int
openmc_filter_get_type(int32_t index, char* type)
{
  if (int err = verify_filter(index)) return err;

  std::strcpy(type, model::tally_filters[index].type().c_str());
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
  for (int i = 0; i < model::n_tally_filters; i++) {
    Filter& t = model::tally_filters[i];
    largest_filter_id = std::max(largest_filter_id, t.id());
  }
  *id = largest_filter_id + 1;
}

extern "C" int
openmc_new_filter(const char* type, int32_t* index)
{
  printf("Error: creation of filters from C API not currently supported!\n");
  exit(1);
  /*
  *index = model::tally_filters.size();
  Filter::create(type);
  */
  return 0;
}

} // namespace openmc
