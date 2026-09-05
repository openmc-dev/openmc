#ifndef OPENMC_FIELD_H
#define OPENMC_FIELD_H

#include <memory>
#include <string>
#include <unordered_map>

#include "fmt/format.h"
#include "pugixml.hpp"

#include "openmc/mesh.h"
#include "openmc/position.h"
#include "openmc/vector.h"

namespace openmc {

class Field;

namespace model {

//! Global container of all field instances, indexed by position.
extern vector<unique_ptr<Field>> fields;

//! Map from user-facing field ID to index in fields
extern std::unordered_map<int, int> field_map;

} // namespace model

// -----------------------------------------------------------
// FieldData
// -----------------------------------------------------------

//! Container for field data values defined in a mesh.
//!
//! Wraps a flat vector of values and provides element-wise access.
//! \tparam T Value type (e.g., double for scalar fields, Direction for vector
//!           fields)
template<typename T>
class FieldData {
public:
  // Constructor
  FieldData(vector<T> values) : values_(std::move(values)) {}

  //! Get the value for a given index.
  //
  //! \param[in] idx Index
  //! \return Associated value
  T evaluate(int idx) const { return values_[idx]; }

  //! Assign a value at a given index.
  //
  //! \param[in] idx Index
  //! \param[in] value Value to store
  void assign(int idx, T value) { values_[idx] = value; }

  //! Return the size of the data field.
  //
  //! \return Number of values
  int size() const { return values_.size(); }

  // Values accessors
  vector<T>& values() { return values_; }
  const vector<T>& values() const { return values_; }

private:
  vector<T> values_; //!< Stored data
};

// -----------------------------------------------------------
// Field
// -----------------------------------------------------------

//! Describes how field values are mapped onto a mesh
enum class FieldMapping {
  NODAL, // Nodal representation (values defined for each vertex)
  CELL   // Cell-based representation (values defined for each cell)
};

//! Abstract base class for all fields defined on a geometric mesh.
//
//! Provides common metadata (ID, type string, name) and a non-owning
//! pointer to the associated mesh. Concrete field behaviour is
//! implemented by MappedField and its subclasses.
class Field {
public:
  virtual ~Field() = default;

  Mesh* mesh_ptr() const
  {
    if (mesh_ == nullptr) {
      fatal_error("No mesh found for this field!");
    } else {
      return mesh_;
    }
  }

  FieldMapping mapping() const { return mapping_; }

protected:
  Mesh* mesh_;           //!< Non-owning pointer to the geometric mesh
  FieldMapping mapping_; //!< Relationship between values and mesh entities
                         //!< (nodes vs. cells)
};

// -----------------------------------------------------------
// MappedField
// -----------------------------------------------------------

//! A field that stores type data on a mesh with a specified mapping.
//!
//! Provides evaluation (including trilinear interpolation for nodal fields),
//! assignment, and mesh-boundary queries.
//!
//! \tparam T Value type stored per mesh entity (e.g., double for scalar fields,
//!           Direction for vector fields).
template<typename T>
class MappedField : public Field {
public:
  MappedField() = default;

  //! Construct a fully initialized MappedField.
  //
  //! \param[in] mesh_ptr Non-owning pointer to the mesh
  //! \param[in] values Field values.
  //! \param[in] mapping Mapping type: 'nodal' or 'cell'
  MappedField(Mesh* mesh_ptr, vector<T> values, std::string mapping)
  {
    set_mesh(mesh_ptr);
    set_mapping(mapping);

    std::unique_ptr<FieldData<T>> data =
      std::make_unique<FieldData<T>>(std::move(values));
    set_data(std::move(data));
  }

  //! Set the mesh pointer.
  //! Returns an error if the mesh pointer is not valid.
  //
  //! \param[in] value Mesh pointer
  void set_mesh(Mesh* value)
  {
    if (value != nullptr) {
      mesh_ = value;
    } else {
      fatal_error("No mesh found for this field!");
    }
  }

  //! Set the mapping type.
  //! Returns an error if the mapping type is not valid.
  //
  //! \param[in] value Mapping type
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

  //! Assign a value to a bin.
  //
  //! \param[in] bin Bin number
  //! \param[in] value Value to store
  void assign(int bin, T value) { data().assign(bin, value); }

  //! Return the mesh bin associated with a given position.
  //
  //! \param[in] r Position
  //! \return Bin number
  int get_bin(const Position& r) { return mesh_ptr()->get_bin(r); }

  //! Evaluate field at a given position inside the mesh knowing the current
  //! bin.
  //
  //! \param[in] r Position
  //! \param[in] bin Bin number corresponding to r
  //! \return Value corresponding to r
  T evaluate_in_mesh(const Position& r, int bin)
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

  //! Evaluate field at a given position knowing the previous position and bin.
  //
  //! The bin used to evaluate the field at r1 can either be:
  //! - the bin corresponding to r1, if r1 is inside the mesh,
  //! - the bin corresponding to the last crossed bin inside the mesh during ray
  //!   traversal analysis between r0 and r1, if r1 is outside the mesh,
  //! - bin0, if r1 is outside the mesh and no bins were found during ray
  //!   traversal.
  //
  //! \param[in] r0 Previous position (inside the mesh)
  //! \param[in] r1 Current position
  //! \param[in] bin0 Bin corresponding to r0
  //! \return Value corresponding to r1 relative to a clamped bin
  T evaluate_clamped(const Position& r0, const Position& r1, int bin0)
  {
    int next_bin = mesh_ptr()->get_bin_clamped(r0, r1, bin0);
    return evaluate_in_mesh(r1, next_bin);
  }

  //! Interpolate data field at a given position using trilinear interpolation.
  //! To avoid extrapolation, any normalized coordinate not in [0.0, 1.0] will
  //! trigger an error.
  //
  //! The position r and the bin number does not have to be related.
  //
  //! \param[in] r Position
  //! \param[in] bin Bin number
  //! \return Interpolated value
  T trilinear_interpolation(const Position& r, int bin)
  {
    // Normalize coordinates
    Position n_r = mesh_ptr()->normalize_coordinates(r, bin);

    // Protect from extrapolation
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
  //! \return Distance to the next mesh boundary
  double distance_to_next_boundary(
    int current_bin, const Position& r, const Direction& u, int& bin_next)
  {
    return this->mesh_ptr()->distance_to_next_boundary(
      current_bin, r, u, bin_next);
  }

  // Data field accessor
  FieldData<T>& data() const
  {
    if (data_ == nullptr) {
      fatal_error("No data found for this field!");
    } else {
      return *data_;
    }
  }

  // Data field value accessors
  const T value(int i) const { return data().evaluate(i); }
  const vector<T> values() const { return data().values(); }

private:
  std::unique_ptr<FieldData<T>> data_; //!< Data associated with the mesh
};

// -----------------------------------------------------------
// TemperatureField
// -----------------------------------------------------------

class TemperatureField : public MappedField<double> {
public:
  // Constructors
  TemperatureField() : MappedField<double>() {};
  TemperatureField(
    Mesh* mesh_ptr, vector<double> values, std::string mapping = "cell")
    : MappedField<double>(mesh_ptr, values, mapping) {};

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
};

// -----------------------------------------------------------
// VelocityField
// -----------------------------------------------------------

// Boundary conditions type
enum class BCType { NONE, INLET, OUTLET, WALL };

// Boundary conditions map type
using BCMap = std::unordered_map<BCType, vector<int>>;

class VelocityField : public MappedField<Direction> {
public:
  // Constructors
  VelocityField() : MappedField<Direction>() {};
  VelocityField(Mesh* mesh_ptr, vector<Direction> values, std::string mapping)
    : MappedField<Direction>(mesh_ptr, values, mapping) {};

  //! Find next bin associated with a given position (r1) knowing the previous
  //! position (r0) and the previous bin (bin0). The next bin is evaluated using
  //! raytracing. If r1 is outside the mesh, crossed_boundary will indicate
  //! the boundary condition associated with the last surface crossed before
  //! leaving the mesh. Also, intersection will be the last intersection with
  //! the mesh. If r1 is still inside the mesh, crossed_boundary will be NONE
  //! and the intersection will coincide with r1.
  //!
  //! We currently do not check whether the mesh was left during the travel,
  //! meaning that a point can leave and reenter the mesh. It is not problematic
  //! for convex geometries like regular meshes, but it can be for unstuctured
  //! meshes.
  //
  //! \param[in] r0 First position
  //! \param[in] r1 Second position
  //! \param[in] bin0 Bin number corresponding to r0
  //! \param[out] crossed_boundary Boundary type of the last crossed surface
  //! \param[out] intersection Last intersection with the mesh (or r1)
  //! \return Bin number corresponding to r1
  int get_next_bin(const Position& r0, const Position& r1, int bin0,
    BCType& crossed_boundary, Position& intersection);

  //! Update the position and the bin by moving the point to a random location
  //! on the inlet physical group.
  //
  //! \param[inout] p Position
  //! \param[inout] bin Bin number corresponding to p
  //! \param[in] seed Random number generator seed
  void randomly_place_on_inlet(Position& p, int& bin, uint64_t* seed);

  //! Retrieve the boundary condition associated with a physical group.
  //
  //! \param[in] physical_group Physical group
  //! \return Boundary condition
  BCType get_boundary_condition(int physical_group);

  // Boundary conditions map accessors
  BCMap& bc_map() { return bc_map_; }
  const BCMap& bc_map() const { return bc_map_; }

private:
  BCMap bc_map_; //!< Boundary conditions map linking a boundary condition type
                 //!< to physical group numbers
};

template<typename T>
T* get_field(int field_id)
{
  auto it = model::field_map.find(field_id);
  if (it == model::field_map.end()) {
    fatal_error(
      fmt::format("No <field> element with id={} was found.", field_id));
  }

  Field* base = model::fields[it->second].get();
  auto* field = dynamic_cast<T*>(base);
  if (!field) {
    fatal_error(fmt::format("Field id={} is of incompatible type.", field_id));
  }

  return field;
}

void read_fields(const pugi::xml_node& root);

} //  namespace openmc

#endif // OPENMC_FIELD_H
