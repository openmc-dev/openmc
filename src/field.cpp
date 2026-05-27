#include "openmc/field.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/mesh.h"
#include "openmc/simulation.h"
#include "openmc/vector.h"

namespace openmc {

// -----------------------------------------------------------
// TemperatureField implementation
// -----------------------------------------------------------

TemperatureField::TemperatureField(
  Mesh* mesh_ptr, std::vector<double> values, std::string mapping)
{
  set_mesh(mesh_ptr);
  set_mapping(mapping);

  std::unique_ptr<FieldData<double>> data =
    std::make_unique<FieldData<double>>(values);
  set_data(std::move(data));
}

double TemperatureField::get_temperature(int bin)
{
  if (bin >= 0 && bin < values().size()) {
    return this->value(bin);
  }
  return -1.0;
}

double TemperatureField::get_sqrtkT(int bin)
{
  double temperature = get_temperature(bin);
  if (temperature >= 0) {
    return sqrt(temperature * K_BOLTZMANN);
  }
  return -1.0;
}

// -----------------------------------------------------------
// VelocityField implementation
// -----------------------------------------------------------

VelocityField::VelocityField(
  Mesh* mesh_ptr, std::vector<Direction> values, std::string mapping)
{
  set_mesh(mesh_ptr);
  set_mapping(mapping);

  std::unique_ptr<FieldData<Direction>> data =
    std::make_unique<FieldData<Direction>>(values);
  set_data(std::move(data));
}

int VelocityField::get_next_bin(const Position& r0, const Position& r1,
  int bin0, BCType& crossed_boundary, Position& intersection)
{
  // TODO - implement early mesh departure detection

  // Initialization
  double total_length = 0.0;
  vector<int> outward_surface_ids;
  vector<int> inward_surface_ids;
  vector<int> bins;
  vector<double> segment_lengths;

  // Raytracing
  mesh_ptr()->full_raytracing(
    r0, r1, outward_surface_ids, inward_surface_ids, bins, segment_lengths);

  // Consistency check
  if (segment_lengths.size() == 0) {
    fatal_error("Inconsistency in raytrace results.");
  }

  // If next point outside the mesh
  if (outward_surface_ids.size() != inward_surface_ids.size()) {

    // Retrieve ID from the last surface
    int surface_id = outward_surface_ids.back();

    // Calculate total length travelled from r0 to intersection
    for (auto l : segment_lengths) {
      total_length += l;
    }

    // Calculate intersection
    intersection = r0 + (r1 - r0) * total_length;

    // Translate surface ID in physical group
    int physical_group = mesh_ptr()->get_physical_group(surface_id);

    // Determine crossed_boundary type
    crossed_boundary = get_boundary_condition(physical_group);

    // Return bin default value when outside the mesh
    return C_NONE;
  }

  // Next point inside the mesh
  intersection = r1;
  crossed_boundary = BCType::NONE;
  return bins.back();
}

void VelocityField::randomly_place_on_inlet(
  Position& p, int& bin, uint64_t* seed)
{
  mesh_ptr()->randomly_place_on_physical_group(
    p, bin, seed, bc_map_[BCType::INLET]);
}

BCType VelocityField::get_boundary_condition(int physical_group)
{
  if (!bc_map_.empty()) {
    for (auto& pair : bc_map_) {
      std::vector<int>& vec = pair.second;
      if (std::find(vec.begin(), vec.end(), physical_group) != vec.end()) {
        return pair.first;
      }
    }
    fatal_error("Not found!");
  } else {
    fatal_error("Empty map for boundary conditions in the velocity field!");
  }
}

//==============================================================================
// C API
//==============================================================================

extern "C" int openmc_temperature_field_set_temperature(
  int32_t index, double temperature)
{
  if (index < 0 || index >= simulation::temperature_field.values().size()) {
    set_errmsg("Index in temperature field is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  simulation::temperature_field.assign(index, temperature);
  return 0;
}

} // namespace openmc
