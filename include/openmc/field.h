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
  ScalarField(
    Mesh* mesh_ptr, vector<double> values, const std::string& field_name);
  ScalarField(Mesh* mesh_ptr, vector<double> values)
    : ScalarField(mesh_ptr, values, "ScalarField") {};

  //----------------------------------------------------------------------------
  // Methods

  //! Find the next mesh crossing given a particle position and direction. If
  //! the particle is initially outside, the crossing will be at the nearest
  //! outer boundary of the mesh along its direction of travel.
  //
  //! \param[in] current_bin Current bin number
  //! \param[in] r Position of the particle
  //! \param[in] u Direction of the particle
  //! \return Distance to the crossing and the next bin number
  MeshCrossing next_mesh_crossing(
    int current_bin, const Position& r, const Direction& u);

  //----------------------------------------------------------------------------
  // Accessors

  // Field type
  const std::string& field_type() const { return this->field_type_; }

  // Mesh pointer
  Mesh* mesh_ptr() const
  {
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
  TemperatureField(Mesh* mesh_ptr, vector<double> values);

  //----------------------------------------------------------------------------
  // Methods

  //! Check that a temperature can be used with the currently loaded data
  //
  //! Rejects values that are not finite or that are negative. Once cross
  //! section data have been read, also rejects values outside the temperature
  //! range those data cover. Before that point no range is known and only the
  //! finiteness and sign of the value are checked.
  //
  //! \param[in] temperature Temperature in Kelvin
  //! \throws std::runtime_error if the temperature cannot be used
  static void validate_temperature(double temperature);

  //! Returns the temperature in Kelvin corresponding to a given bin number
  //! relative to the mesh.
  //
  //! \param[in] bin Bin number
  //! \return Temperature in Kelvin
  double get_temperature(int bin);

  //! Returns the square root of the temperature multiplied by the Boltzmann
  //! constant in eV for a given bin number relative to the mesh.
  //
  //! \param[in] bin Bin number
  //! \return Sqrt(k_Boltzmann * temperature) in eV
  double get_sqrtkT(int bin);

  //! Returns the bin number corresponding to the location of the particle.
  //
  //! \param[in] r Position of the particle
  //! \return Corresponding bin number or -1 if outside the mesh
  int get_bin(const Position& r);
};

} //  namespace openmc
#endif // OPENMC_FIELD_H
