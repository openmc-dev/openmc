#ifndef OPENMC_FIELD_H
#define OPENMC_FIELD_H

#include "openmc/mesh.h"
#include "openmc/vector.h"

namespace openmc {

class ScalarField {
public:
  //----------------------------------------------------------------------------
  // Constructors
  ScalarField() = default;
  ScalarField(Mesh* mesh_ptr, vector<double> values, const std::string& field_name);
  ScalarField(Mesh* mesh_ptr, vector<double> values) : ScalarField(mesh_ptr, values, "ScalarField") {};

  //----------------------------------------------------------------------------
  // Methods

  //! Returns the distance to the next mesh boundary gien a particle position
  //! and direction. If the particle is initially outside, the distance will
  //! correspond to the nearest distance to the outer boundaries of the mesh.
  //
  //! \param[in] r Position of the particle
  //! \param[in] u Direction of the particle
  //! \return The distance in cm to the next mesh boundary
  double distance_to_next_boundary(const Position& r, const Direction& d);

  //----------------------------------------------------------------------------
  // Accessors

  // Field type
  const std::string& field_type() const { return this->field_type_; }

  // Mesh pointer
  Mesh* mesh_ptr() const {
    if (this->mesh_ptr_ == nullptr) {
      fatal_error(fmt::format("No mesh found for {}!", this->field_type_));
    } else {
      return this->mesh_ptr_;
    }
  }

  // Values
  double& value(int i) { return values_[i]; }
  const double& value(int i) const { return values_[i]; }
  const vector<double>& values() const { return values_; }

private:
  //----------------------------------------------------------------------------
  // Data members
  std::string field_type_; //! Name of field type
  Mesh* mesh_ptr_;         //!< Pointer to the geometric mesh
  vector<double> values_;  //!< Values associated with each mesh cell
};

class TemperatureField : public ScalarField {
public:
  //----------------------------------------------------------------------------
  // Constructors
  TemperatureField() = default;
  TemperatureField(Mesh* mesh_ptr, vector<double> values)
    : ScalarField(mesh_ptr, values, "TemperatureField") {};

  //----------------------------------------------------------------------------
  // Methods

  //! Returns the temperature in Kelvin corresponding to the given position.
  //
  //! \param[in] r Position of the particle
  //! \return Temperature in Kelvin
  double get_temperature(const Position& r);

  //! Returns the square root of the temperature multiplied by the Boltzmann
  //! constant in eV for the given position.
  //
  //! \param[in] r Position of the particle
  //! \return Sqrt(k_Boltzmann * temperature) in eV
  double get_sqrtkT(const Position& r);

  //! Update the temperature of a particle based on its position and direction.
  //! If the particle is inside the temperature field, its temperature is
  //! updated. If outside, the particle takes the temperature value
  //! associated with the current cell instance.
  //
  //! \param[inout] p Particle
  void update_particle_temperature(Particle& p);
};

} //  namespace openmc
#endif // OPENMC_FIELD_H
