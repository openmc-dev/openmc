#include "openmc/field.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/mesh.h"
#include "openmc/vector.h"

namespace openmc {

TemperatureField::TemperatureField(Mesh* mesh_ptr, vector<double> values)
{
  this->mesh_ptr_ = mesh_ptr;
  for (double v : values) {
    this->values_.push_back(v);
  }
}

double TemperatureField::distance_to_next_boundary(
  const Position& r, const Direction& u)
{
  if (this->mesh_ptr_ != nullptr) {
    return this->mesh_ptr_->distance_to_next_boundary(r, u);
  }
  fatal_error("TODO");
}

double TemperatureField::get_temperature(const Position& r)
{
  if (this->mesh_ptr_ != nullptr) {

    // Get bin from position
    int i = this->mesh_ptr_->get_bin(r);

    // If we have a bin then use it to locate the value
    if (i >= 0 && i < this->values_.size()) {
      // Return the temperature
      return this->values_[i];
    }

    // Return -1.0 if no values were found (probably outside of the mesh)
    return -1.0;

  } else {
    fatal_error("TODO");
  }
}

double TemperatureField::get_sqrtkT(const Position& r)
{
  double temperature = this->get_temperature(r);
  if (temperature >= 0) {
    return sqrt(temperature * K_BOLTZMANN);
  } else {
    // TODO
    fatal_error("");
    return temperature;
  }
}

void TemperatureField::update_particle_temperature(Particle& p) {

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
