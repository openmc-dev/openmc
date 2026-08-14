#include "openmc/field.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/mesh.h"
#include "openmc/simulation.h"
#include "openmc/vector.h"

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
    fatal_error(
      fmt::format("The number of elements in the mesh is not consistent with the "
                  "number of values declared in {}!",
        field_type));
  }

  for (double v : values) {
    this->values_.push_back(v);
  }
}

double ScalarField::distance_to_next_boundary(
  int current_bin, const Position& r, const Direction& u, int& bin_next)
{
  return this->mesh_ptr()->distance_to_next_boundary(
    current_bin, r, u, bin_next);
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
