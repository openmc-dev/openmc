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
    std::make_unique<SimpleFieldData<double>>(values);
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

int TemperatureField::get_bin(const Position& r)
{
  int bin = mesh_ptr()->get_bin(r);
  if (bin >= 0 && bin < values().size()) {
    return bin;
  } else {
    return C_NONE;
  }
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
    std::make_unique<SimpleFieldData<Direction>>(values);
  set_data(std::move(data));
}

Position VelocityField::find_departure_from_mesh(
  Position pa, Position pb, BCType& crossed_boundary)
{
  int physical_group;
  Position intersection = mesh_ptr()->departure_from_mesh(
    pb, pa, (pa - pb) / (pa.distance(pb)), physical_group);
  crossed_boundary = get_boundary_condition(physical_group);
  return intersection;
}

void VelocityField::randomly_place_on_inlet(
  Position& pa, int& cell, uint64_t* seed)
{
  mesh_ptr()->randomly_place_on_physical_group(
    pa, cell, seed, bc_map_[BCType::INLET]);
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
