#include "openmc/field.h"
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/mesh.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/vector.h"

#include <cmath>
#include <stdexcept>

namespace openmc {

ScalarField::ScalarField(
  Mesh* mesh_ptr, vector<double> values, const std::string& field_type)
{
  this->field_type_ = field_type;

  if (mesh_ptr != nullptr) {
    this->mesh_ptr_ = mesh_ptr;
  } else {
    fatal_error(fmt::format("No mesh found for {}!", field_type));
  }

  if (this->mesh_ptr_->n_bins() != values.size()) {
    fatal_error(fmt::format(
      "The number of elements in the mesh is not consistent with the "
      "number of values declared in {}!",
      field_type));
  }

  for (double v : values) {
    this->values_.push_back(v);
  }
}

MeshCrossing ScalarField::next_mesh_crossing(
  int current_bin, const Position& r, const Direction& u)
{
  return this->mesh_ptr()->next_mesh_crossing(current_bin, r, u);
}

void TemperatureField::validate_temperature(double temperature)
{
  if (!std::isfinite(temperature) || temperature < 0.0) {
    throw std::runtime_error {fmt::format(
      "Temperature of {} K is not a finite, non-negative value.", temperature)};
  }

  // data::temperature_min/max only bound anything once cross sections have
  // been read. Until then they hold their sentinel values of INFTY and 0.0,
  // and there is no range to check against.
  if (data::temperature_min > data::temperature_max)
    return;

  if (temperature < data::temperature_min - settings::temperature_tolerance) {
    throw std::runtime_error {
      fmt::format("Temperature of {} K is below minimum temperature at "
                  "which data is available of {} K.",
        temperature, data::temperature_min)};
  }
  if (temperature > data::temperature_max + settings::temperature_tolerance) {
    throw std::runtime_error {
      fmt::format("Temperature of {} K is above maximum temperature at "
                  "which data is available of {} K.",
        temperature, data::temperature_max)};
  }
}

TemperatureField::TemperatureField(Mesh* mesh_ptr, vector<double> values)
  : ScalarField(mesh_ptr, values, "TemperatureField")
{
  for (int i = 0; i < this->values().size(); ++i) {
    try {
      validate_temperature(this->value(i));
    } catch (const std::runtime_error& e) {
      throw std::runtime_error {
        fmt::format("Element {} of the temperature field: {}", i, e.what())};
    }
  }
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

  try {
    TemperatureField::validate_temperature(temperature);
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_INVALID_ARGUMENT;
  }

  simulation::temperature_field.value(index) = temperature;
  return 0;
}

extern "C" size_t openmc_temperature_field_size()
{
  return simulation::temperature_field.values().size();
}

extern "C" int openmc_temperature_field_get_value(
  int32_t index, double* temperature)
{
  if (index < 0 || index >= simulation::temperature_field.values().size()) {
    set_errmsg("Index in temperature field is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *temperature = simulation::temperature_field.value(index);
  return 0;
}

} // namespace openmc
