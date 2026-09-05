#include "openmc/field.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/mesh.h"
#include "openmc/simulation.h"
#include "openmc/vector.h"

namespace openmc {

namespace model {

vector<unique_ptr<Field>> fields;
std::unordered_map<int32_t, int32_t> field_map;

} // namespace model

// -----------------------------------------------------------
// TemperatureField implementation
// -----------------------------------------------------------

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

// -----------------------------------------------------------
// VelocityField implementation
// -----------------------------------------------------------

int VelocityField::get_next_bin(const Position& r0, const Position& r1,
  int bin0, BCType& crossed_boundary, Position& intersection)
{
  // TODO - implement early mesh departure detection

  // The algorithm used to determine if r1 is outside the mesh assumes that
  // r0 is inside the mesh
  if (bin0 < 0 || bin0 >= mesh_ptr()->n_bins()) {
    fatal_error("r0 must be inside the mesh!");
  }

  double total_length = 0.0;
  vector<int> outward_surface_ids;
  vector<int> inward_surface_ids;
  vector<int> bins;
  vector<double> segment_lengths;

  mesh_ptr()->bins_and_surface_bins_crossed(
    r0, r1, outward_surface_ids, inward_surface_ids, bins, segment_lengths);

  if (segment_lengths.size() == 0) {

    // Try to move the DNP back a little in case it is right on a mesh boundary
    Direction u = r1 - r0;
    Position r0_nudged = r0 - (u / u.norm() * TINY_BIT);

    // Reset vectors
    outward_surface_ids.clear();
    inward_surface_ids.clear();
    bins.clear();
    segment_lengths.clear();

    mesh_ptr()->bins_and_surface_bins_crossed(r0_nudged, r1,
      outward_surface_ids, inward_surface_ids, bins, segment_lengths);

    if (segment_lengths.size() == 0) {
      fatal_error("Inconsistency in raytrace results.");
    }

    // Adjust total length
    total_length -= TINY_BIT;
  }

  // If next point outside the mesh
  if (outward_surface_ids.size() != inward_surface_ids.size()) {

    // Retrieve ID from the last surface
    int surface_id = outward_surface_ids.back();

    // Calculate total length travelled from r0 to intersection
    for (auto l : segment_lengths) {
      total_length += l;
    }

    // Calculate intersection
    intersection = r0 + (r1 - r0) * total_length;

    // Translate surface ID in physical group
    int physical_group = mesh_ptr()->get_physical_group(surface_id);

    // Determine crossed_boundary type
    crossed_boundary = get_boundary_condition(physical_group);

    // Return bin default value when outside the mesh
    return C_NONE;
  }

  // Next point inside the mesh
  intersection = r1;
  crossed_boundary = BCType::NONE;
  return bins.back();
}

void VelocityField::randomly_place_on_inlet(
  Position& p, int& bin, uint64_t* seed)
{
  mesh_ptr()->sample_on_physical_groups(p, bin, seed, bc_map_[BCType::INLET]);
}

BCType VelocityField::get_boundary_condition(int physical_group)
{
  if (!bc_map_.empty()) {
    for (auto& pair : bc_map_) {
      std::vector<int>& vec = pair.second;
      if (std::find(vec.begin(), vec.end(), physical_group) != vec.end()) {
        return pair.first;
      }
    }
    fatal_error("Not found!");
  } else {
    fatal_error("Empty map for boundary conditions in the velocity field!");
  }
}

void read_fields(const pugi::xml_node& root)
{
  for (auto node = root.child("field"); node;
       node = node.next_sibling("field")) {

    // Type
    if (!node.attribute("type")) {
      throw std::runtime_error(
        "A <field> element is missing the 'type' attribute");
    }
    std::string type = node.attribute("type").value();

    // ID
    if (!node.attribute("id")) {
      throw std::runtime_error(fmt::format(
        "A <field> element of type '{}' is missing the 'id' attribute.", type));
    }
    int id = std::stoi(node.attribute("id").value());

    // Check for duplicate
    if (model::field_map.find(id) != model::field_map.end()) {
      throw std::runtime_error(fmt::format("Duplicate field id={}.", id));
    }

    // Mesh
    if (!check_for_node(node, "mesh")) {
      throw std::runtime_error(
        fmt::format("Field id={} must specify a mesh.", id));
    }
    int mesh_id = std::stoi(get_node_value(node, "mesh"));
    if (model::mesh_map.find(mesh_id) == model::mesh_map.end()) {
      throw std::runtime_error(fmt::format(
        "Mesh {} specified in field id={} does not exist.", mesh_id, id));
    }
    Mesh* mesh_ptr = model::meshes[model::mesh_map.at(mesh_id)].get();

    // Values
    if (!check_for_node(node, "values")) {
      throw std::runtime_error(
        fmt::format("Field id={} must specify values.", id));
    }
    auto raw_values = get_node_array<double>(node, "values");

    // Mapping (optional, default "cell")
    std::string mapping = "cell";
    if (check_for_node(node, "mapping")) {
      mapping = get_node_value(node, "mapping");
    }
    if (mapping != "cell" && mapping != "nodal") {
      throw std::runtime_error(fmt::format(
        "Field id={}: unknown mapping '{}'. Must be 'cell' or 'nodal'.", id,
        mapping));
    }

    // Create a field
    std::unique_ptr<Field> field;
    std::string entity_name = (mapping == "cell") ? "element" : "node";
    size_t n_entities =
      (mapping == "cell") ? mesh_ptr->n_bins() : mesh_ptr->n_vertices();

    // Temperature field
    if (type == "temperature") {
      if (raw_values.size() != n_entities) {
        throw std::runtime_error(fmt::format(
          "Field id={}: temperature field has {} values but mesh has {} {}s.",
          id, raw_values.size(), n_entities, entity_name));
      }
      field = std::make_unique<TemperatureField>(mesh_ptr, raw_values, mapping);

      // Velocity field
    } else if (type == "velocity") {
      if (raw_values.size() != 3 * n_entities) {
        throw std::runtime_error(
          fmt::format("Field id={}: velocity field has {} values but expected "
                      "{} (3 x {} {}s).",
            id, raw_values.size(), 3 * n_entities, n_entities, entity_name));
      }
      vector<Position> vel_values;
      vel_values.reserve(n_entities);
      for (size_t i = 0; i < raw_values.size(); i += 3) {
        vel_values.emplace_back(
          raw_values[i], raw_values[i + 1], raw_values[i + 2]);
      }
      field = std::make_unique<VelocityField>(mesh_ptr, vel_values, mapping);

    } else {
      throw std::runtime_error(
        fmt::format("Unknown field type '{}' for field id={}.", type, id));
    }

    // Store
    model::field_map[id] = model::fields.size();
    model::fields.push_back(std::move(field));
  }
}

//==============================================================================
// C API
//==============================================================================

extern "C" int openmc_temperature_field_set_temperature(
  int32_t index, double temperature)
{
  if (index < 0 || index >= simulation::temperature_field->values().size()) {
    set_errmsg("Index in temperature field is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  simulation::temperature_field->assign(index, temperature);
  return 0;
}

extern "C" size_t openmc_temperature_field_size()
{
  return simulation::temperature_field->values().size();
}

extern "C" int openmc_temperature_field_get_value(
  int32_t index, double* temperature)
{
  if (index < 0 || index >= simulation::temperature_field->values().size()) {
    set_errmsg("Index in temperature field is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *temperature = simulation::temperature_field->value(index);
  return 0;
}

} // namespace openmc
