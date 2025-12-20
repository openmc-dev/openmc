#include "openmc/field.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/mesh.h"
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
      fmt::format("The number of bins in the mesh is not consistent with the "
                  "number of values declared in {}!",
        field_type));
  }

  for (double v : values) {
    this->values_.push_back(v);
  }
}

double ScalarField::distance_to_next_boundary(
  const Position& r, const Direction& u)
{
  return this->mesh_ptr()->distance_to_next_boundary(r, u);
}

double TemperatureField::get_temperature(const Position& r)
{
  // Get bin from position
  int i = this->mesh_ptr()->get_bin(r);

  // If we have a bin, we use it to locate the value
  if (i >= 0 && i < this->values().size()) {
    return this->value(i);
  }

  // No values were found (outside the mesh)
  return -1.0;
}

double TemperatureField::get_sqrtkT(const Position& r)
{
  double temperature = this->get_temperature(r);
  if (temperature >= 0) {
    return sqrt(temperature * K_BOLTZMANN);
  }
  return -1.0;
}

void TemperatureField::update_particle_temperature(Particle& p)
{
  // Save current temperature
  p.sqrtkT_last() = p.sqrtkT();

  // Determine the temperature based on the temperature field
  double field_sqrtkT = this->get_sqrtkT(p.r() + p.u() * TINY_BIT);

  // If particle inside the mesh, we use the temperature field
  if (field_sqrtkT >= 0.) {
    p.sqrtkT() = field_sqrtkT;

    // If particle outside the mesh, go back to the cell instance temperature
  } else {
    Cell& c {*model::cells[p.lowest_coord().cell()]};
    p.sqrtkT() = c.sqrtkT(p.cell_instance());
  }
}

} // namespace openmc
