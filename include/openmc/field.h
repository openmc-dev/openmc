#ifndef OPENMC_FIELD_H
#define OPENMC_FIELD_H

#include <any>
#include <unordered_map>

#include "openmc/mesh.h"
#include "openmc/vector.h"

namespace openmc {

template<typename T>
class Field;

template<typename T>
class FieldData {
public:
  virtual T evaluate(int bin, Position r) const = 0;
  virtual int size() = 0;
  virtual std::any values() const = 0;
};

template<typename T>
class SimpleFieldData : public FieldData<T> {
public:
  SimpleFieldData(std::vector<T> values) { values_ = values; }

  T evaluate(int bin, Position r) const { return values_[bin]; }

  int size() { return values_.size(); }

  std::any values() const override { return values_; }

private:
  std::vector<T> values_;
};

template<typename T>
class NestedFieldData : public FieldData<T> {
public:
  NestedFieldData(std::vector<Field<T>> values) { values_ = values; }

  T evaluate(int bin, Position r) const { return values_[bin].evaluate(r); }

  int size() { return values_.size(); }

  std::any values() const override { return values_; }

private:
  std::vector<Field<T>> values_;
};

enum class FieldMapping {
  NODAL, // Nodal representation (value defined for each vertex)
  CELL   // Cell-based representation (value defined for each cell)
};

enum class NodalEvaluation {
  UNDEFINED,     // Undefined evaluation
  INTERPOLATION, // Determine value using interpolation
  CLOSEST        // Determine value using closest point
};

template<typename T>
class Field {
public:
  // Constructors
  Field() = default;

  // Set mesh pointer
  void set_mesh(Mesh* value)
  {
    if (value != nullptr) {
      mesh_ = value;
    } else {
      fatal_error("No mesh found for this field!");
    }
  }

  // Set mapping
  void set_mapping(std::string value)
  {
    if (value == "nodal") {
      mapping_ = FieldMapping::NODAL;
    } else if (value == "cell") {
      mapping_ = FieldMapping::CELL;
    } else {
      fatal_error(fmt::format("Unrecognized mapping type: {}", value));
    }
  }

  // Set nodal evaluation method
  void set_nodal_evaluation(std::string value)
  {
    if (value == "interpolation") {
      nodal_evaluation_ = NodalEvaluation::INTERPOLATION;
    } else if (value == "closest") {
      nodal_evaluation_ = NodalEvaluation::CLOSEST;
    } else {
      fatal_error(fmt::format("Unrecognized nodal evaluation type: {}", value));
    }
  }

  // Set data
  void set_data(std::unique_ptr<FieldData<T>> data)
  {
    // Values/mesh consistency check
    // Check for nodal representation
    if (mapping() == FieldMapping::NODAL) {
      // TODO - check consistency compared to total number of expected unique
      // vertices
    }
    // Check for cell-based representation
    if (mapping() == FieldMapping::CELL) {
      if (mesh_ptr()->n_bins() != data->size()) {
        fatal_error("The number of bins in the mesh is not consistent with the "
                    "number of values declared for this field!");
      }
    }

    // Store data
    data_ = std::move(data);
  }

  // Evaluate a field value inside a mesh cell (nodal representation)
  T evaluate_inside_cell(int bin, Position r)
  {
    if (nodal_evaluation() == NodalEvaluation::INTERPOLATION) {
      return interpolate(bin, r);
    } else if (nodal_evaluation() == NodalEvaluation::CLOSEST) {
      // return closest_value(bin, r);
      //  TODO - remove the nodal evaluation option
    } else {
      fatal_error("Not implemented!");
    }
  }

  // Interpolate the field at a given position
  T interpolate(int bin, Position r)
  {

    Position n_r = mesh_ptr()->normalize_position(r);
    std::vector<int> v = mesh_ptr()->connectivity(bin);

    // Interpolate along x
    T c00 = value(v[0], r) * (1 - n_r[0]) + value(v[1], r) * n_r[0];
    T c01 = value(v[4], r) * (1 - n_r[0]) + value(v[5], r) * n_r[0];
    T c10 = value(v[2], r) * (1 - n_r[0]) + value(v[3], r) * n_r[0];
    T c11 = value(v[6], r) * (1 - n_r[0]) + value(v[7], r) * n_r[0];

    // Interpolate along y
    T c0 = c00 * (1 - n_r[1]) + c10 * n_r[1];
    T c1 = c01 * (1 - n_r[1]) + c11 * n_r[1];

    // Interpolate along z
    T result = c0 * (1 - n_r[2]) + c1 * n_r[2];
    return result;
  }

  //! Return the mesh bin associated with the given position
  //
  //! \param[in] r Position
  //! \return Mesh bin
  int get_mesh_bin(Position r) { return mesh_ptr()->get_bin(r); }

  //! Evaluate the field at a given position
  //
  //! \param[in] r Position where the field is evaluated
  //! \return Value of the field corresponding to the position
  T evaluate(Position r)
  {
    int bin = get_mesh_bin(r);
    if (mapping() == FieldMapping::NODAL) {
      return evaluate_inside_cell(bin, r);
    } else if (mapping() == FieldMapping::CELL) {
      return value(bin, r);
    } else {
      fatal_error("Not implemented!");
    }
  }

  //! Returns the distance to the next mesh boundary given a particle position
  //! and direction. If the particle is initially outside, the distance will
  //! correspond to the nearest distance to the outer boundaries of the mesh.
  //
  //! \param[in] current_bin Current bin number
  //! \param[in] r Position of the particle
  //! \param[in] u Direction of the particle
  //! \param[out] bin_next Next bin number
  //! \return The distance in cm to the next mesh boundary
  double distance_to_next_boundary(
    int current_bin, const Position& r, const Direction& d, int& bin_next);

  // Mapping accessor
  FieldMapping mapping() { return mapping_; }
  const FieldMapping mapping() const { return mapping_; }

  // Nodal evaluation
  NodalEvaluation nodal_evaluation() { return nodal_evaluation_; }
  const NodalEvaluation nodal_evaluation() const { return nodal_evaluation_; }

  // Mesh pointer accessor
  Mesh* mesh_ptr() const
  {
    if (mesh_ == nullptr) {
      fatal_error("No mesh found for this field!");
    } else {
      return mesh_;
    }
  }

  const FieldData<T>& data() const
  {
    if (data_ == nullptr) {
      fatal_error("No data found for this field!");
    } else {
      return *data_;
    }
  }

  T value(int bin, Position r) const { return data().evaluate(bin, r); }

private:
  FieldMapping mapping_; //!< Relationship between values and mesh
  NodalEvaluation nodal_evaluation_ =
    NodalEvaluation::UNDEFINED;        //!< Nodal evaluation method
  Mesh* mesh_;                         //!< Pointer to the geometric mesh
  std::unique_ptr<FieldData<T>> data_; //!< Data associated with the mesh
};

class TemperatureField : public Field<double> {
public:
  // Constructors
  TemperatureField() : Field<double>() {};
  TemperatureField(Mesh* mesh_ptr, std::vector<double> values,
    std::string mapping, std::string nodal_evaluation)
  {
    set_mesh(mesh_ptr);
    set_mapping(mapping);
    set_nodal_evaluation(nodal_evaluation);

    std::unique_ptr<FieldData<double>> data =
      std::make_unique<SimpleFieldData<double>>(values);
    set_data(std::move(data));
  }

  //! Returns the temperature in Kelvin corresponding to a given bin number
  //! relative to the mesh.
  //
  //! \param[in] r Position of the particle
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

  //! \return Sqrt(k_Boltzmann * temperature) in eV
  double get_sqrtkT(const Position& r);

  // Runtime cast to go through the values
  // Should be an extremely limited use case
  // TODO - all_values() method instead that would give a flatten
  // vector of all possible values so that we do not need std::any
  const vector<double> values() const
  {
    return std::any_cast<std::vector<double>>(data().values());
  }
};

// Boundary conditions type
enum class BCType { INLET, OUTLET, WALL };

// Boundary conditions map type
using BCMap = std::unordered_map<BCType, std::vector<int>>;

class VelocityField : public Field<Direction> {
public:
  // Constructors
  VelocityField() : Field<Direction>() {};
  VelocityField(Mesh* mesh_ptr, std::vector<Direction> values,
    std::string mapping, std::string nodal_evaluation);

  Position find_departure_from_mesh(
    Position pa, Position pb, BCType& crossed_boundary);

  void randomly_place_on_inlet(Position& pa, int& cell, uint64_t* seed);

  BCType get_boundary_condition(int physical_group);

  // Boundary conditions map accessors
  BCMap& bc_map() { return bc_map_; }
  const BCMap& bc_map() const { return bc_map_; }

private:
  BCMap bc_map_; //!< Boundary conditions map linking a boundary
                 //!< condition type to physical group numbers
};

} //  namespace openmc

#endif // OPENMC_FIELD_H
