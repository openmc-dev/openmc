#include "openmc/field.h"
#include "openmc/mesh.h"
#include "openmc/vector.h"

namespace openmc {

TemperatureField::TemperatureField(Mesh* mesh_ptr, vector<double> values) 
{
  this->mesh_ptr = mesh_ptr;
  for (double v: values) {
    this->values.push_back(v);
  }
}

double TemperatureField::distance_to_next_cell(Position r, Direction u) 
{
  if (this->mesh_ptr != nullptr) {
    return this->mesh_ptr->distance_to_next_boundary(r, u);
  }
  fatal_error("TODO");
}

double TemperatureField::get_temperature(Position r)
{
  // Get bin from position
  int i = this->mesh_ptr->get_bin(r);

  // If we have a bin then use it to locate the value
  if (i >= 0 && i < this->values.size()) {
    // Return the temperature
    return this->values[i];
  }
  
  // Return -1.0 if no values were found (probably outside of the mesh)  
  return -1.0;
}

double TemperatureField::get_sqrtkT(Position r)
{
  double temperature = this->get_temperature(r);
  if (temperature >= 0) {
    return sqrt(temperature * K_BOLTZMANN);
  } else {
    //TODO
    fatal_error("");
    return temperature;
  }
}

} // namespace openmc
