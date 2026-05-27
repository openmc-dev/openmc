#ifndef OPENMC_FIELD_H
#define OPENMC_FIELD_H

#include <any>
#include <unordered_map>

#include "openmc/mesh.h"
#include "openmc/vector.h"

namespace openmc {

// -----------------------------------------------------------
// FieldData
// -----------------------------------------------------------

template<typename T>
class Field;

template<typename T>
class FieldData {
public:
  FieldData(std::vector<T> values) { values_ = values; }

  T evaluate(int bin) const { return values_[bin]; }

  void assign(int bin, T value) { values_[bin] = value; }

  int size() { return values_.size(); }

  std::vector<T> values() const { return values_; }

private:
  std::vector<T> values_;
};

// -----------------------------------------------------------
// Field
// -----------------------------------------------------------

enum class FieldMapping {
  NODAL, // Nodal representation (value defined for each vertex)
  CELL   // Cell-based representation (value defined for each cell)
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

  //! Set data field.
  //! Returns an error if the size of the data field is not consistent with its
  //! mapping type.
  //
  //! \param[in] data Data field
  void set_data(std::unique_ptr<FieldData<T>> data)
  {
    // Values/mesh size consistency check
    if (mapping() == FieldMapping::NODAL) {
      if (mesh_ptr()->n_vertices() != data->size()) {
        fatal_error("The number of bins in the mesh is not consistent with the "
                    "number of values declared for this field!");
      }
    } else if (mapping() == FieldMapping::CELL) {
      if (mesh_ptr()->n_bins() != data->size()) {
        fatal_error("The number of bins in the mesh is not consistent with the "
                    "number of values declared for this field!");
      }
    }

    data_ = std::move(data);
  }



  //! Return the mesh bin associated with the given position
  //
  //! \param[in] r Position
  //! \return Mesh bin
  int get_mesh_bin(Position r) { return mesh_ptr()->get_bin(r); }

  //! Evaluate the field at a given position with prior knowledge of the current
  //! bin.
  //
  //! \param[in] r Position
  //! \param[in] bin Bin number corresponding to the current position
  //! \return Value corresponding to the position
  T evaluate(const Position& r, int bin)
  {
    if (bin != C_NONE) {
      switch (mapping()) {
      case FieldMapping::NODAL:
        // TODO: implement other interpolation techniques
        return trilinear_interpolation(r, bin);
        break;
      case FieldMapping::CELL:
        return value(bin);
        break;
      default:
        fatal_error("Not implemented for this mapping type!");
        break;
      }
    } else {
      fatal_error("Bin outside the mesh.");
    }
  }

  //! Evaluate the field at a given position with prior knowledge of the bin and
  //! of the previous position. We assume that the transport between the two
  //! points is in a straight line, justifying the use of ray-tracing. If the
  //! next point is outside the mesh, we use the last bin that was located
  //! inside the mesh.
  //
  //! \param[in] r0 Previous position
  //! \param[in] r1 Current position
  //! \param[in] bin Bin corresponding to the previous position
  //! \return Value of the field corresponding to the current position
  T evaluate(Position r0, Position r1, int bin)
  {
    int next_bin = mesh_ptr()->get_last_bin_inside_mesh(r0, r1, bin);
    return evaluate(r1, next_bin);
  }

  //! Interpolate data field at a given position using trilinear interpolation.
  //! To avoid extrapolation, any normalized coordinate not in [0.0, 1.0] will
  //! trigger an error.
  //
  //! \param[in] r Position
  //! \param[in] bin Bin number corresponding to the position
  //! \return Interpolated value
  T trilinear_interpolation(const Position& r, int bin)
  {
    // Normalize coordinates
    Position n_r = mesh_ptr()->normalize_position(r);

    // Check that normalized coordinates are in [0.0, 1.0] to avoid
    // extrapolation
    for (int i = 0; i < 3; i++) {
      if ((n_r[i] < 0.0) || (n_r[i] > 1.0)) {
        fatal_error(
          "Normalized coordinates must be in [0.0, 1.0] for interpolation!");
      }
    }

    // Retrieve vertices
    vector<int> v = mesh_ptr()->connectivity(bin);

    // Interpolate along x
    T c00 = value(v[0]) * (1 - n_r[0]) + value(v[1]) * n_r[0];
    T c01 = value(v[4]) * (1 - n_r[0]) + value(v[5]) * n_r[0];
    T c10 = value(v[2]) * (1 - n_r[0]) + value(v[3]) * n_r[0];
    T c11 = value(v[6]) * (1 - n_r[0]) + value(v[7]) * n_r[0];

    // Interpolate along y
    T c0 = c00 * (1 - n_r[1]) + c10 * n_r[1];
    T c1 = c01 * (1 - n_r[1]) + c11 * n_r[1];

    // Interpolate along z
    T result = c0 * (1 - n_r[2]) + c1 * n_r[2];
    return result;
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
    int current_bin, const Position& r, const Direction& u, int& bin_next)
  {
    return this->mesh_ptr()->distance_to_next_boundary(
      current_bin, r, u, bin_next);
  }

  // Mapping accessor
  FieldMapping mapping() { return mapping_; }
  const FieldMapping mapping() const { return mapping_; }

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

  FieldData<T>& data()
  {
    if (data_ == nullptr) {
      fatal_error("No data found for this field!");
    } else {
      return *data_;
    }
  }

  T value(int bin) const { return data().evaluate(bin); }

  void assign(int bin, T value) { data().assign(bin, value); }

  //! Returns the list of values corresponding to the data stored in the field.
  //
  //! \return List of values
  const vector<T> values() const { return data().values(); }

private:
  FieldMapping mapping_;               //!< Relationship between values and mesh
  Mesh* mesh_;                         //!< Pointer to the geometric mesh
  std::unique_ptr<FieldData<T>> data_; //!< Data associated with the mesh
};

// -----------------------------------------------------------
// TemperatureField
// -----------------------------------------------------------

class TemperatureField : public Field<double> {
public:
  // Constructors
  TemperatureField() : Field<double>() {};
  TemperatureField(
    Mesh* mesh_ptr, std::vector<double> values, std::string mapping = "cell");

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
};

// -----------------------------------------------------------
// VelocityField
// -----------------------------------------------------------

// Boundary conditions type
enum class BCType { NONE, INLET, OUTLET, WALL };

// Boundary conditions map type
using BCMap = std::unordered_map<BCType, std::vector<int>>;

class VelocityField : public Field<Direction> {
public:
  // Constructors
  VelocityField() : Field<Direction>() {};
  VelocityField(
    Mesh* mesh_ptr, std::vector<Direction> values, std::string mapping);

  int get_next_bin(Position r0, Position r1, int current_bin,
    BCType& crossed_boundary, Position& intersection);

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
