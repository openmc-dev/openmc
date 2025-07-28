
#include "openmc/cell.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iterator>
#include <set>
#include <sstream>
#include <string>

#include <fmt/core.h>
#include <gsl/gsl-lite.hpp>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/dagmc.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
std::unordered_map<int32_t, int32_t> cell_map;
vector<unique_ptr<Cell>> cells;

std::unordered_map<int32_t, int32_t> universe_map;
vector<unique_ptr<Universe>> universes;
} // namespace model

//==============================================================================
//! Convert region specification string to integer tokens.
//!
//! The characters (, ), |, and ~ count as separate tokens since they represent
//! operators.
//==============================================================================

vector<int32_t> tokenize(const std::string region_spec)
{
  // Check for an empty region_spec first.
  vector<int32_t> tokens;
  if (region_spec.empty()) {
    return tokens;
  }

  // Parse all halfspaces and operators except for intersection (whitespace).
  for (int i = 0; i < region_spec.size();) {
    if (region_spec[i] == '(') {
      tokens.push_back(OP_LEFT_PAREN);
      i++;

    } else if (region_spec[i] == ')') {
      tokens.push_back(OP_RIGHT_PAREN);
      i++;

    } else if (region_spec[i] == '|') {
      tokens.push_back(OP_UNION);
      i++;

    } else if (region_spec[i] == '~') {
      tokens.push_back(OP_COMPLEMENT);
      i++;

    } else if (region_spec[i] == '-' || region_spec[i] == '+' ||
               std::isdigit(region_spec[i])) {
      // This is the start of a halfspace specification.  Iterate j until we
      // find the end, then push-back everything between i and j.
      int j = i + 1;
      while (j < region_spec.size() && std::isdigit(region_spec[j])) {
        j++;
      }
      tokens.push_back(std::stoi(region_spec.substr(i, j - i)));
      i = j;

    } else if (std::isspace(region_spec[i])) {
      i++;

    } else {
      auto err_msg =
        fmt::format("Region specification contains invalid character, \"{}\"",
          region_spec[i]);
      fatal_error(err_msg);
    }
  }

  // Add in intersection operators where a missing operator is needed.
  int i = 0;
  while (i < tokens.size() - 1) {
    bool left_compat {(tokens[i] < OP_UNION) || (tokens[i] == OP_RIGHT_PAREN)};
    bool right_compat {(tokens[i + 1] < OP_UNION) ||
                       (tokens[i + 1] == OP_LEFT_PAREN) ||
                       (tokens[i + 1] == OP_COMPLEMENT)};
    if (left_compat && right_compat) {
      tokens.insert(tokens.begin() + i + 1, OP_INTERSECTION);
    }
    i++;
  }

  return tokens;
}

vector<vector<int32_t>> \
  generate_triso_distribution(\
    vector<int> lattice_shape,\
    vector<double> lattice_pitch,\
    vector<double> lattice_lower_left,\
    vector<std::int32_t> cell_rpn, int id)
{
  vector<vector<int32_t>> triso_distribution(lattice_shape[0]*lattice_shape[1]*lattice_shape[2]);
  vector<double> mesh_center(3);
  vector<int> mesh_ind(3);
  /*
  for (int i=0; i<lattice_shape[0]; i++) {
    for (int j=0; j<lattice_shape[1]; j++) {
      for (int k=0; k<lattice_shape[2]; k++) {
        mesh_center[0]=(i+0.5)*lattice_pitch[0]+lattice_lower_left[0];
        mesh_center[1]=(j+0.5)*lattice_pitch[1]+lattice_lower_left[1];
        mesh_center[2]=(k+0.5)*lattice_pitch[2]+lattice_lower_left[2];
        for (int32_t token : cell_rpn) {
          if (token >= OP_UNION) continue;
          if (model::surfaces[abs(token) - 1]->triso_in_mesh(mesh_center, lattice_pitch)) {
            triso_distribution[i+j*lattice_shape[0]+k*lattice_shape[0]*lattice_shape[1]].push_back(token);
            model::surfaces[abs(token) - 1]->connect_to_triso_base(id, "base");
          }
        }
      }
    }
  }*/
  for (int32_t token : cell_rpn) {
    if (token >= OP_UNION) continue;
    vector<double> triso_center=model::surfaces[abs(token) - 1]->get_center();
    for (int i=0; i<3; i++) {
      mesh_ind[i]=floor((triso_center[i]-lattice_lower_left[i])/lattice_pitch[i]);
    }
    for (int i=mesh_ind[0]-1; i<=mesh_ind[0]+1; i++) {
      for (int j=mesh_ind[1]-1; j<=mesh_ind[1]+1; j++) {
        for (int k=mesh_ind[2]-1; k<=mesh_ind[2]+1; k++) {
          if (i < 0 || i >= lattice_shape[0] ||\
              j < 0 || j >= lattice_shape[1] ||\
              k < 0 || k >= lattice_shape[2]) continue;
          mesh_center[0]=(i+0.5)*lattice_pitch[0]+lattice_lower_left[0];
          mesh_center[1]=(j+0.5)*lattice_pitch[1]+lattice_lower_left[1];
          mesh_center[2]=(k+0.5)*lattice_pitch[2]+lattice_lower_left[2];
          if (model::surfaces[abs(token) - 1]->triso_in_mesh(mesh_center, lattice_pitch)) {
            triso_distribution[i+j*lattice_shape[0]+k*lattice_shape[0]*lattice_shape[1]].push_back(token);
            model::surfaces[abs(token) - 1]->connect_to_triso_base(id, "base");
          }
        }
      }
    }
  }


  return triso_distribution;
}


//==============================================================================
//! Convert infix region specification to Reverse Polish Notation (RPN)
//!
//! This function uses the shunting-yard algorithm.
//==============================================================================

vector<int32_t> generate_rpn(int32_t cell_id, vector<int32_t> infix)
{
  vector<int32_t> rpn;
  vector<int32_t> stack;

  for (int32_t token : infix) {
    if (token < OP_UNION) {
      // If token is not an operator, add it to output
      rpn.push_back(token);
    } else if (token < OP_RIGHT_PAREN) {
      // Regular operators union, intersection, complement
      while (stack.size() > 0) {
        int32_t op = stack.back();

        if (op < OP_RIGHT_PAREN && ((token == OP_COMPLEMENT && token < op) ||
                                     (token != OP_COMPLEMENT && token <= op))) {
          // While there is an operator, op, on top of the stack, if the token
          // is left-associative and its precedence is less than or equal to
          // that of op or if the token is right-associative and its precedence
          // is less than that of op, move op to the output queue and push the
          // token on to the stack. Note that only complement is
          // right-associative.
          rpn.push_back(op);
          stack.pop_back();
        } else {
          break;
        }
      }

      stack.push_back(token);

    } else if (token == OP_LEFT_PAREN) {
      // If the token is a left parenthesis, push it onto the stack
      stack.push_back(token);

    } else {
      // If the token is a right parenthesis, move operators from the stack to
      // the output queue until reaching the left parenthesis.
      for (auto it = stack.rbegin(); *it != OP_LEFT_PAREN; it++) {
        // If we run out of operators without finding a left parenthesis, it
        // means there are mismatched parentheses.
        if (it == stack.rend()) {
          fatal_error(fmt::format(
            "Mismatched parentheses in region specification for cell {}",
            cell_id));
        }
        rpn.push_back(stack.back());
        stack.pop_back();
      }

      // Pop the left parenthesis.
      stack.pop_back();
    }
  }

  while (stack.size() > 0) {
    int32_t op = stack.back();

    // If the operator is a parenthesis it is mismatched.
    if (op >= OP_RIGHT_PAREN) {
      fatal_error(fmt::format(
        "Mismatched parentheses in region specification for cell {}", cell_id));
    }

    rpn.push_back(stack.back());
    stack.pop_back();
  }

  return rpn;
}

//==============================================================================
// Universe implementation
//==============================================================================

void Universe::to_hdf5(hid_t universes_group) const
{
  // Create a group for this universe.
  auto group = create_group(universes_group, fmt::format("universe {}", id_));

  // Write the geometry representation type.
  write_string(group, "geom_type", "csg", false);

  // Write the contained cells.
  if (cells_.size() > 0) {
    vector<int32_t> cell_ids;
    for (auto i_cell : cells_)
      cell_ids.push_back(model::cells[i_cell]->id_);
    write_dataset(group, "cells", cell_ids);
  }

  close_group(group);
}

bool Universe::find_cell(Particle& p) const
{
  if (filled_with_triso_base_ != -1) {
    Cell& c {*model::cells[model::cell_map[filled_with_triso_base_]]};
    vector<int> lat_ind(3);
    Position r {p.r_local()};
    lat_ind[0]=floor((r.x-c.vl_lower_left_[0])/c.vl_pitch_[0]);
    lat_ind[1]=floor((r.y-c.vl_lower_left_[1])/c.vl_pitch_[1]);
    lat_ind[2]=floor((r.z-c.vl_lower_left_[2])/c.vl_pitch_[2]);
    int32_t i_univ = p.coord(p.n_coord() - 1).universe;
    for (int token : c.vl_triso_distribution_[lat_ind[0]+lat_ind[1]*c.vl_shape_[0]+lat_ind[2]*c.vl_shape_[0]*c.vl_shape_[1]]) {
      vector<double> triso_center=model::surfaces[abs(token) - 1]->get_center();
      double triso_radius=model::surfaces[abs(token) - 1]->get_radius();
      if (model::cells[model::cell_map[model::surfaces[abs(token) - 1]->triso_base_index_]]->universe_!= i_univ) continue;
      if (abs(token)==abs(p.surface())) {
        if (p.surface() < 0) {
          p.coord(p.n_coord() - 1).cell = model::cell_map[model::surfaces[abs(token) - 1]->triso_particle_index_];
          return true;
        } else {
          p.coord(p.n_coord() - 1).cell = model::cell_map[filled_with_triso_base_];
          return true;
        }
      }
      if (pow(r.x-triso_center[0],2)+pow(r.y-triso_center[1],2)+pow(r.z-triso_center[2],2) < pow(triso_radius,2)) {
        p.coord(p.n_coord() - 1).cell = model::cell_map[model::surfaces[abs(token) - 1]->triso_particle_index_];
        return true;
      }
    }
    if (model::cells[model::cell_map[filled_with_triso_base_]]->universe_== i_univ) {
      p.coord(p.n_coord() - 1).cell = model::cell_map[filled_with_triso_base_];
      return true;
    }
  }
  const auto& cells {
    !partitioner_ ? cells_ : partitioner_->get_cells(p.r_local(), p.u_local())};

  for (auto it = cells.begin(); it != cells.end(); it++) {
    int32_t i_cell = *it;
    int32_t i_univ = p.coord(p.n_coord() - 1).universe;
    if (model::cells[i_cell]->universe_ != i_univ)
      continue;

    // Check if this cell contains the particle;
    Position r {p.r_local()};
    Direction u {p.u_local()};
    auto surf = p.surface();
    if (model::cells[i_cell]->contains(r, u, surf)) {
      p.coord(p.n_coord() - 1).cell = i_cell;
      return true;
    }
  }
  return false;
}

BoundingBox Universe::bounding_box() const
{
  BoundingBox bbox = {INFTY, -INFTY, INFTY, -INFTY, INFTY, -INFTY};
  if (cells_.size() == 0) {
    return {};
  } else {
    for (const auto& cell : cells_) {
      auto& c = model::cells[cell];
      bbox |= c->bounding_box();
    }
  }
  return bbox;
}

//==============================================================================
// Cell implementation
//==============================================================================

void Cell::set_rotation(const vector<double>& rot)
{
  if (fill_ == C_NONE) {
    fatal_error(fmt::format("Cannot apply a rotation to cell {}"
                            " because it is not filled with another universe",
      id_));
  }

  if (rot.size() != 3 && rot.size() != 9) {
    fatal_error(fmt::format("Non-3D rotation vector applied to cell {}", id_));
  }

  // Compute and store the rotation matrix.
  rotation_.clear();
  rotation_.reserve(rot.size() == 9 ? 9 : 12);
  if (rot.size() == 3) {
    double phi = -rot[0] * PI / 180.0;
    double theta = -rot[1] * PI / 180.0;
    double psi = -rot[2] * PI / 180.0;
    rotation_.push_back(std::cos(theta) * std::cos(psi));
    rotation_.push_back(-std::cos(phi) * std::sin(psi) +
                        std::sin(phi) * std::sin(theta) * std::cos(psi));
    rotation_.push_back(std::sin(phi) * std::sin(psi) +
                        std::cos(phi) * std::sin(theta) * std::cos(psi));
    rotation_.push_back(std::cos(theta) * std::sin(psi));
    rotation_.push_back(std::cos(phi) * std::cos(psi) +
                        std::sin(phi) * std::sin(theta) * std::sin(psi));
    rotation_.push_back(-std::sin(phi) * std::cos(psi) +
                        std::cos(phi) * std::sin(theta) * std::sin(psi));
    rotation_.push_back(-std::sin(theta));
    rotation_.push_back(std::sin(phi) * std::cos(theta));
    rotation_.push_back(std::cos(phi) * std::cos(theta));

    // When user specifies angles, write them at end of vector
    rotation_.push_back(rot[0]);
    rotation_.push_back(rot[1]);
    rotation_.push_back(rot[2]);
  } else {
    std::copy(rot.begin(), rot.end(), std::back_inserter(rotation_));
  }
}

double Cell::temperature(int32_t instance) const
{
  if (sqrtkT_.size() < 1) {
    throw std::runtime_error {"Cell temperature has not yet been set."};
  }

  if (instance >= 0) {
    double sqrtkT = sqrtkT_.size() == 1 ? sqrtkT_.at(0) : sqrtkT_.at(instance);
    return sqrtkT * sqrtkT / K_BOLTZMANN;
  } else {
    return sqrtkT_[0] * sqrtkT_[0] / K_BOLTZMANN;
  }
}

void Cell::set_temperature(double T, int32_t instance, bool set_contained)
{
  if (settings::temperature_method == TemperatureMethod::INTERPOLATION) {
    if (T < data::temperature_min) {
      throw std::runtime_error {"Temperature is below minimum temperature at "
                                "which data is available."};
    } else if (T > data::temperature_max) {
      throw std::runtime_error {"Temperature is above maximum temperature at "
                                "which data is available."};
    }
  }

  if (type_ == Fill::MATERIAL) {
    if (instance >= 0) {
      // If temperature vector is not big enough, resize it first
      if (sqrtkT_.size() != n_instances_)
        sqrtkT_.resize(n_instances_, sqrtkT_[0]);

      // Set temperature for the corresponding instance
      sqrtkT_.at(instance) = std::sqrt(K_BOLTZMANN * T);
    } else {
      // Set temperature for all instances
      for (auto& T_ : sqrtkT_) {
        T_ = std::sqrt(K_BOLTZMANN * T);
      }
    }
  } else {
    if (!set_contained) {
      throw std::runtime_error {
        fmt::format("Attempted to set the temperature of cell {} "
                    "which is not filled by a material.",
          id_)};
    }

    auto contained_cells = this->get_contained_cells(instance);
    for (const auto& entry : contained_cells) {
      auto& cell = model::cells[entry.first];
      Expects(cell->type_ == Fill::MATERIAL);
      auto& instances = entry.second;
      for (auto instance : instances) {
        cell->set_temperature(T, instance);
      }
    }
  }
}

void Cell::export_properties_hdf5(hid_t group) const
{
  // Create a group for this cell.
  auto cell_group = create_group(group, fmt::format("cell {}", id_));

  // Write temperature in [K] for one or more cell instances
  vector<double> temps;
  for (auto sqrtkT_val : sqrtkT_)
    temps.push_back(sqrtkT_val * sqrtkT_val / K_BOLTZMANN);
  write_dataset(cell_group, "temperature", temps);

  close_group(cell_group);
}

void Cell::import_properties_hdf5(hid_t group)
{
  auto cell_group = open_group(group, fmt::format("cell {}", id_));

  // Read temperatures from file
  vector<double> temps;
  read_dataset(cell_group, "temperature", temps);

  // Ensure number of temperatures makes sense
  auto n_temps = temps.size();
  if (n_temps > 1 && n_temps != n_instances_) {
    throw std::runtime_error(fmt::format(
      "Number of temperatures for cell {} doesn't match number of instances",
      id_));
  }

  // Modify temperatures for the cell
  sqrtkT_.clear();
  sqrtkT_.resize(temps.size());
  for (gsl::index i = 0; i < temps.size(); ++i) {
    this->set_temperature(temps[i], i);
  }

  close_group(cell_group);
}

void Cell::to_hdf5(hid_t cell_group) const
{

  // Create a group for this cell.
  auto group = create_group(cell_group, fmt::format("cell {}", id_));

  if (!name_.empty()) {
    write_string(group, "name", name_, false);
  }

  write_dataset(group, "universe", model::universes[universe_]->id_);

  to_hdf5_inner(group);

  // Write fill information.
  if (type_ == Fill::MATERIAL) {
    write_dataset(group, "fill_type", "material");
    std::vector<int32_t> mat_ids;
    for (auto i_mat : material_) {
      if (i_mat != MATERIAL_VOID) {
        mat_ids.push_back(model::materials[i_mat]->id_);
      } else {
        mat_ids.push_back(MATERIAL_VOID);
      }
    }
    if (mat_ids.size() == 1) {
      write_dataset(group, "material", mat_ids[0]);
    } else {
      write_dataset(group, "material", mat_ids);
    }

    std::vector<double> temps;
    for (auto sqrtkT_val : sqrtkT_)
      temps.push_back(sqrtkT_val * sqrtkT_val / K_BOLTZMANN);
    write_dataset(group, "temperature", temps);

  } else if (type_ == Fill::UNIVERSE) {
    write_dataset(group, "fill_type", "universe");
    write_dataset(group, "fill", model::universes[fill_]->id_);
    if (translation_ != Position(0, 0, 0)) {
      write_dataset(group, "translation", translation_);
    }
    if (!rotation_.empty()) {
      if (rotation_.size() == 12) {
        std::array<double, 3> rot {rotation_[9], rotation_[10], rotation_[11]};
        write_dataset(group, "rotation", rot);
      } else {
        write_dataset(group, "rotation", rotation_);
      }
    }

  } else if (type_ == Fill::LATTICE) {
    write_dataset(group, "fill_type", "lattice");
    write_dataset(group, "lattice", model::lattices[fill_]->id_);
  }

  close_group(group);
}

//==============================================================================
// CSGCell implementation
//==============================================================================

// default constructor
CSGCell::CSGCell()
{
  geom_type_ = GeometryType::CSG;
}

CSGCell::CSGCell(pugi::xml_node cell_node)
{
  geom_type_ = GeometryType::CSG;

  if (check_for_node(cell_node, "id")) {
    id_ = std::stoi(get_node_value(cell_node, "id"));
  } else {
    fatal_error("Must specify id of cell in geometry XML file.");
  }

  if (check_for_node(cell_node, "name")) {
    name_ = get_node_value(cell_node, "name");
  }

  if (check_for_node(cell_node, "universe")) {
    universe_ = std::stoi(get_node_value(cell_node, "universe"));
  } else {
    universe_ = 0;
  }

  // Check if the cell is the base of a virtual triso lattice
  bool virtual_lattice_present = check_for_node(cell_node, "virtual_lattice");
  if (virtual_lattice_present) {
    virtual_lattice_ = get_node_value_bool(cell_node, "virtual_lattice");
    if (virtual_lattice_) {
      if (check_for_node(cell_node, "lower_left") && check_for_node(cell_node, "pitch") && check_for_node(cell_node, "shape")) {
        vl_lower_left_ = get_node_array<double>(cell_node, "lower_left");
        vl_pitch_ = get_node_array<double>(cell_node, "pitch");
        vl_shape_ = get_node_array<int>(cell_node, "shape");
      } else {
        fatal_error(fmt::format("Lower_left, pitch and shape of the virtual lattice must be specified for cell {}", id_));
      }
    }
  } else {
    virtual_lattice_ = false;
  }
  
  if (check_for_node(cell_node, "triso_particle")) {
    triso_particle_ = get_node_value_bool(cell_node, "triso_particle");
  } else {
    triso_particle_ = false;
  }

  // Make sure that either material or fill was specified, but not both.
  bool fill_present = check_for_node(cell_node, "fill");
  bool material_present = check_for_node(cell_node, "material");
  if (!(fill_present || material_present)) {
    fatal_error(
      fmt::format("Neither material nor fill was specified for cell {}", id_));
  }
  if (fill_present && material_present) {
    fatal_error(fmt::format("Cell {} has both a material and a fill specified; "
                            "only one can be specified per cell",
      id_));
  }

  if (fill_present) {
    fill_ = std::stoi(get_node_value(cell_node, "fill"));
    if (fill_ == universe_) {
      fatal_error(fmt::format("Cell {} is filled with the same universe that"
                              "it is contained in.",
        id_));
    }
  } else {
    fill_ = C_NONE;
  }

  // Read the material element.  There can be zero materials (filled with a
  // universe), more than one material (distribmats), and some materials may
  // be "void".
  if (material_present) {
    vector<std::string> mats {
      get_node_array<std::string>(cell_node, "material", true)};
    if (mats.size() > 0) {
      material_.reserve(mats.size());
      for (std::string mat : mats) {
        if (mat.compare("void") == 0) {
          material_.push_back(MATERIAL_VOID);
        } else {
          material_.push_back(std::stoi(mat));
        }
      }
    } else {
      fatal_error(fmt::format(
        "An empty material element was specified for cell {}", id_));
    }
  }

  // Read the temperature element which may be distributed like materials.
  if (check_for_node(cell_node, "temperature")) {
    sqrtkT_ = get_node_array<double>(cell_node, "temperature");
    sqrtkT_.shrink_to_fit();

    // Make sure this is a material-filled cell.
    if (material_.size() == 0) {
      fatal_error(fmt::format(
        "Cell {} was specified with a temperature but no material. Temperature"
        "specification is only valid for cells filled with a material.",
        id_));
    }

    // Make sure all temperatures are non-negative.
    for (auto T : sqrtkT_) {
      if (T < 0) {
        fatal_error(fmt::format(
          "Cell {} was specified with a negative temperature", id_));
      }
    }

    // Convert to sqrt(k*T).
    for (auto& T : sqrtkT_) {
      T = std::sqrt(K_BOLTZMANN * T);
    }
  }

  // Read the region specification.
  std::string region_spec;
  if (check_for_node(cell_node, "region")) {
    region_spec = get_node_value(cell_node, "region");
  }

  // Get a tokenized representation of the region specification.
  region_ = tokenize(region_spec);
  region_.shrink_to_fit();

  // Convert user IDs to surface indices.
  for (auto& r : region_) {
    if (r < OP_UNION) {
      const auto& it {model::surface_map.find(abs(r))};
      if (it == model::surface_map.end()) {
        throw std::runtime_error {
          "Invalid surface ID " + std::to_string(abs(r)) +
          " specified in region for cell " + std::to_string(id_) + "."};
      }
      r = (r > 0) ? it->second + 1 : -(it->second + 1);
    }
  }

  // Convert the infix region spec to RPN.
  rpn_ = generate_rpn(id_, region_);

  // Check if this is a simple cell.
  simple_ = true;
  for (int32_t token : rpn_) {
    if ((token == OP_COMPLEMENT) || (token == OP_UNION)) {
      simple_ = false;
      break;
    }
  }

  // If this cell is simple, remove all the superfluous operator tokens.
  if (simple_) {
    size_t i0 = 0;
    size_t i1 = 0;
    while (i1 < rpn_.size()) {
      if (rpn_[i1] < OP_UNION) {
        rpn_[i0] = rpn_[i1];
        ++i0;
      }
      ++i1;
    }
    rpn_.resize(i0);
  }
  rpn_.shrink_to_fit();

  if (virtual_lattice_) {
    vl_triso_distribution_ = generate_triso_distribution(vl_shape_, vl_pitch_, vl_lower_left_, rpn_, id_);
    write_message("successfully generated virtual lattice", 5);
    for (int tri_sur_id : vl_triso_distribution_[floor(vl_shape_[0]/2)+floor(vl_shape_[1]/2)*vl_shape_[0]+\
    floor(vl_shape_[2]/2)*vl_shape_[0]*vl_shape_[1]])
    {
      std::cout<<model::surfaces[abs(tri_sur_id)-1]->id_<<std::endl;

    }
  }

  if (triso_particle_) {
    if (rpn_.size() != 1) {
      fatal_error(fmt::format("Wrong surface definition of triso particle cell {}", id_));
    } else {
      model::surfaces[abs(rpn_[0]) - 1]->connect_to_triso_base(id_, "particle");
    }
  }

  // Read the translation vector.
  if (check_for_node(cell_node, "translation")) {
    if (fill_ == C_NONE) {
      fatal_error(fmt::format("Cannot apply a translation to cell {}"
                              " because it is not filled with another universe",
        id_));
    }

    auto xyz {get_node_array<double>(cell_node, "translation")};
    if (xyz.size() != 3) {
      fatal_error(
        fmt::format("Non-3D translation vector applied to cell {}", id_));
    }
    translation_ = xyz;
  }

  // Read the rotation transform.
  if (check_for_node(cell_node, "rotation")) {
    auto rot {get_node_array<double>(cell_node, "rotation")};
    set_rotation(rot);
  }
}

//==============================================================================

bool CSGCell::contains(Position r, Direction u, int32_t on_surface) const
{
  if (simple_) {
    return contains_simple(r, u, on_surface);
  } else {
    return contains_complex(r, u, on_surface);
  }
}

//==============================================================================

std::pair<double, int32_t> CSGCell::distance(
  Position r, Direction u, int32_t on_surface, Particle* p) const
{
  double min_dist {INFTY};
  int32_t i_surf {std::numeric_limits<int32_t>::max()};
  double min_dis_vl;
  int32_t i_surf_vl;
  if (virtual_lattice_) {
    /*
    double tol_dis = p.collision_distance();
    double u_value = sqrt(pow(u.x,2)+pow(u.y,2)+pow(u.z,2)); //don't know if u has been normalized
    vector<double> norm_u = {u.x/u_value, u.y/u_value, u.z/u_value};
    vector<double> r_end = {r.x+tol_dis*norm_u[0], r.y+tol_dis*norm_u[1], r.z+tol_dis*norm_u[2]};
    double temp_pos_x, temp_pos_y, temp_pos_z;
    vector<vector<int>> passed_lattice={{floor((r.x-vl_lower_left_[0])/vl_pitch_[0]),\
                                          floor((r.y-vl_lower_left_[1])/vl_pitch_[1]),\
                                          floor((r.z-vl_lower_left_[2])/vl_pitch_[2]), 0}};
    int index_start;
    int index_end;
    for (int i = 0; i < 3; i++){
      if (passed_lattice[0][i] == vl_shape_[i]) {
        passed_lattice[0][i] = vl_shape_[i]-1;
      }
    }

    if (u.x > 0) {
      index_start = ceil((r.x-vl_lower_left_[0])/vl_pitch_[0]);
      index_end = floor((r_end[0]-vl_lower_left_[0])/vl_pitch_[0]);
      if (index_start <= index_end && index_start < vl_shape_[0]) {
        for (int i = index_start; i <= index_end; i++) {
          if (i >= vl_shape_[0]) break;
          temp_pos_x = i*vl_pitch_[0]+vl_lower_left_[0];
          temp_pos_y = (temp_pos_x-r.x)*norm_u[1]/norm_u[0]+r.y;
          temp_pos_z = (temp_pos_x-r.x)*norm_u[2]/norm_u[0]+r.z;
          passed_lattice.push_back({i, floor((temp_pos_y-vl_lower_left_[1])/vl_pitch_[1]),\
                                    floor((temp_pos_z-vl_lower_left_[2])/vl_pitch_[2]),\
                                    (temp_pos_x-r.x)/norm_u[0]})
        }
      }
    } 
    */

    double max_dis = p->collision_distance();
    double tol_dis = 0;
    vector<double> dis_to_bou(3), dis_to_bou_max(3);
    double u_value = sqrt(pow(u.x,2)+pow(u.y,2)+pow(u.z,2)); //don't know if u has been normalized
    vector<double> norm_u = {u.x/u_value, u.y/u_value, u.z/u_value};
    vector<int> lat_ind(3);
    vector<double> temp_pos ={r.x,r.y,r.z};
    int loop_time;
    for (int i=0; i<3; i++) {
      lat_ind[i]=floor((temp_pos[i]-vl_lower_left_[i])/vl_pitch_[i]);
    }

    dis_to_bou = {INFTY,INFTY,INFTY};
    for (int i=0; i<3; i++) {
      if (norm_u[i] > 0) {
        dis_to_bou[i]=std::abs(((lat_ind[i]+1)*vl_pitch_[i]+vl_lower_left_[i]-temp_pos[i])/norm_u[i]);
        dis_to_bou_max[i]=vl_pitch_[i]/norm_u[i];
      } else if (norm_u[i] < 0){
        dis_to_bou[i]=std::abs((lat_ind[i]*vl_pitch_[i]+vl_lower_left_[i]-temp_pos[i])/norm_u[i]);
        dis_to_bou_max[i]=-vl_pitch_[i]/norm_u[i];
      }
    }

    while(true) {
      
      

      if (lat_ind[0] < 0 || lat_ind[0] >= vl_shape_[0] ||\
          lat_ind[1] < 0 || lat_ind[1] >= vl_shape_[1] ||\
          lat_ind[2] < 0 || lat_ind[2] >= vl_shape_[2]) break;

      
      for (int token : vl_triso_distribution_[lat_ind[0]+lat_ind[1]*vl_shape_[0]+lat_ind[2]*vl_shape_[0]*vl_shape_[1]]) {
        bool coincident {std::abs(token) == std::abs(on_surface)};
        double d {model::surfaces[abs(token) - 1]->distance(r, u, coincident)};
        if (d < min_dist) {
          if (min_dist - d >= FP_PRECISION * min_dist) {
            min_dist = d;
            i_surf = -token;
          }
        }
      }
      

      int mes_bou_crossed=0;
      if  (dis_to_bou[1] < dis_to_bou[0]) {
        mes_bou_crossed=1;
      }
      if  (dis_to_bou[2] < dis_to_bou[mes_bou_crossed]) {
        mes_bou_crossed=2;
      }

      tol_dis = dis_to_bou[mes_bou_crossed];

      if (min_dist < tol_dis) {
        break;
      }

      if (norm_u[mes_bou_crossed] > 0) {
        lat_ind[mes_bou_crossed] += 1;
      } else {
        lat_ind[mes_bou_crossed] += -1;
      }

      dis_to_bou[mes_bou_crossed] += dis_to_bou_max[mes_bou_crossed];
      
      if (tol_dis > max_dis) {
        break;
      }
    }

    
  } else {
    for (int32_t token : rpn_) {
      // Ignore this token if it corresponds to an operator rather than a region.
      if (token >= OP_UNION)
        continue;

      // Calculate the distance to this surface.
      // Note the off-by-one indexing
      bool coincident {std::abs(token) == std::abs(on_surface)};
      double d {model::surfaces[abs(token) - 1]->distance(r, u, coincident)};

      // Check if this distance is the new minimum.
      if (d < min_dist) {
        if (min_dist - d >= FP_PRECISION * min_dist) {
          min_dist = d;
          i_surf = -token;
        }
      }
    }
  }
  
  return {min_dist, i_surf};
}

//==============================================================================

void CSGCell::to_hdf5_inner(hid_t group_id) const
{

  write_string(group_id, "geom_type", "csg", false);

  // Write the region specification.
  if (!region_.empty()) {
    std::stringstream region_spec {};
    for (int32_t token : region_) {
      if (token == OP_LEFT_PAREN) {
        region_spec << " (";
      } else if (token == OP_RIGHT_PAREN) {
        region_spec << " )";
      } else if (token == OP_COMPLEMENT) {
        region_spec << " ~";
      } else if (token == OP_INTERSECTION) {
      } else if (token == OP_UNION) {
        region_spec << " |";
      } else {
        // Note the off-by-one indexing
        auto surf_id = model::surfaces[abs(token) - 1]->id_;
        region_spec << " " << ((token > 0) ? surf_id : -surf_id);
      }
    }
    write_string(group_id, "region", region_spec.str(), false);
  }
}

BoundingBox CSGCell::bounding_box_simple() const
{
  BoundingBox bbox;
  for (int32_t token : rpn_) {
    bbox &= model::surfaces[abs(token) - 1]->bounding_box(token > 0);
  }
  return bbox;
}

void CSGCell::apply_demorgan(
  vector<int32_t>::iterator start, vector<int32_t>::iterator stop)
{
  while (start < stop) {
    if (*start < OP_UNION) {
      *start *= -1;
    } else if (*start == OP_UNION) {
      *start = OP_INTERSECTION;
    } else if (*start == OP_INTERSECTION) {
      *start = OP_UNION;
    }
    start++;
  }
}

vector<int32_t>::iterator CSGCell::find_left_parenthesis(
  vector<int32_t>::iterator start, const vector<int32_t>& rpn)
{
  // start search at zero
  int parenthesis_level = 0;
  auto it = start;
  while (it != rpn.begin()) {
    // look at two tokens at a time
    int32_t one = *it;
    int32_t two = *(it - 1);

    // decrement parenthesis level if there are two adjacent surfaces
    if (one < OP_UNION && two < OP_UNION) {
      parenthesis_level--;
      // increment if there are two adjacent operators
    } else if (one >= OP_UNION && two >= OP_UNION) {
      parenthesis_level++;
    }

    // if the level gets to zero, return the position
    if (parenthesis_level == 0) {
      // move the iterator back one before leaving the loop
      // so that all tokens in the parenthesis block are included
      it--;
      break;
    }

    // continue loop, one token at a time
    it--;
  }
  return it;
}

void CSGCell::remove_complement_ops(vector<int32_t>& rpn)
{
  auto it = std::find(rpn.begin(), rpn.end(), OP_COMPLEMENT);
  while (it != rpn.end()) {
    // find the opening parenthesis (if any)
    auto left = find_left_parenthesis(it, rpn);
    vector<int32_t> tmp(left, it + 1);

    // apply DeMorgan's law to any surfaces/operators between these
    // positions in the RPN
    apply_demorgan(left, it);
    // remove complement operator
    rpn.erase(it);
    // update iterator position
    it = std::find(rpn.begin(), rpn.end(), OP_COMPLEMENT);
  }
}

BoundingBox CSGCell::bounding_box_complex(vector<int32_t> rpn)
{
  // remove complements by adjusting surface signs and operators
  remove_complement_ops(rpn);

  vector<BoundingBox> stack(rpn.size());
  int i_stack = -1;

  for (auto& token : rpn) {
    if (token == OP_UNION) {
      stack[i_stack - 1] = stack[i_stack - 1] | stack[i_stack];
      i_stack--;
    } else if (token == OP_INTERSECTION) {
      stack[i_stack - 1] = stack[i_stack - 1] & stack[i_stack];
      i_stack--;
    } else {
      i_stack++;
      stack[i_stack] = model::surfaces[abs(token) - 1]->bounding_box(token > 0);
    }
  }

  Ensures(i_stack == 0);
  return stack.front();
}

BoundingBox CSGCell::bounding_box() const
{
  return simple_ ? bounding_box_simple() : bounding_box_complex(rpn_);
}

//==============================================================================

bool CSGCell::contains_simple(Position r, Direction u, int32_t on_surface) const
{
  for (int32_t token : rpn_) {
    // Assume that no tokens are operators. Evaluate the sense of particle with
    // respect to the surface and see if the token matches the sense. If the
    // particle's surface attribute is set and matches the token, that
    // overrides the determination based on sense().
    if (token == on_surface) {
    } else if (-token == on_surface) {
      return false;
    } else {
      // Note the off-by-one indexing
      bool sense = model::surfaces[abs(token) - 1]->sense(r, u);
      if (sense != (token > 0)) {
        return false;
      }
    }
  }
  return true;
}

//==============================================================================

bool CSGCell::contains_complex(
  Position r, Direction u, int32_t on_surface) const
{
  // Make a stack of booleans.  We don't know how big it needs to be, but we do
  // know that rpn.size() is an upper-bound.
  vector<bool> stack(rpn_.size());
  int i_stack = -1;

  for (int32_t token : rpn_) {
    // If the token is a binary operator (intersection/union), apply it to
    // the last two items on the stack. If the token is a unary operator
    // (complement), apply it to the last item on the stack.
    if (token == OP_UNION) {
      stack[i_stack - 1] = stack[i_stack - 1] || stack[i_stack];
      i_stack--;
    } else if (token == OP_INTERSECTION) {
      stack[i_stack - 1] = stack[i_stack - 1] && stack[i_stack];
      i_stack--;
    } else if (token == OP_COMPLEMENT) {
      stack[i_stack] = !stack[i_stack];
    } else {
      // If the token is not an operator, evaluate the sense of particle with
      // respect to the surface and see if the token matches the sense. If the
      // particle's surface attribute is set and matches the token, that
      // overrides the determination based on sense().
      i_stack++;
      if (token == on_surface) {
        stack[i_stack] = true;
      } else if (-token == on_surface) {
        stack[i_stack] = false;
      } else {
        // Note the off-by-one indexing
        bool sense = model::surfaces[abs(token) - 1]->sense(r, u);
        stack[i_stack] = (sense == (token > 0));
      }
    }
  }

  if (i_stack == 0) {
    // The one remaining bool on the stack indicates whether the particle is
    // in the cell.
    return stack[i_stack];
  } else {
    // This case occurs if there is no region specification since i_stack will
    // still be -1.
    return true;
  }
}

//==============================================================================
// UniversePartitioner implementation
//==============================================================================

UniversePartitioner::UniversePartitioner(const Universe& univ)
{
  // Define an ordered set of surface indices that point to z-planes.  Use a
  // functor to to order the set by the z0_ values of the corresponding planes.
  struct compare_surfs {
    bool operator()(const int32_t& i_surf, const int32_t& j_surf) const
    {
      const auto* surf = model::surfaces[i_surf].get();
      const auto* zplane = dynamic_cast<const SurfaceZPlane*>(surf);
      double zi = zplane->z0_;
      surf = model::surfaces[j_surf].get();
      zplane = dynamic_cast<const SurfaceZPlane*>(surf);
      double zj = zplane->z0_;
      return zi < zj;
    }
  };
  std::set<int32_t, compare_surfs> surf_set;

  // Find all of the z-planes in this universe.  A set is used here for the
  // O(log(n)) insertions that will ensure entries are not repeated.
  for (auto i_cell : univ.cells_) {
    for (auto token : model::cells[i_cell]->rpn_) {
      if (token < OP_UNION) {
        auto i_surf = std::abs(token) - 1;
        const auto* surf = model::surfaces[i_surf].get();
        if (const auto* zplane = dynamic_cast<const SurfaceZPlane*>(surf))
          surf_set.insert(i_surf);
      }
    }
  }

  // Populate the surfs_ vector from the ordered set.
  surfs_.insert(surfs_.begin(), surf_set.begin(), surf_set.end());

  // Populate the partition lists.
  partitions_.resize(surfs_.size() + 1);
  for (auto i_cell : univ.cells_) {
    // It is difficult to determine the bounds of a complex cell, so add complex
    // cells to all partitions.
    if (!model::cells[i_cell]->simple_) {
      for (auto& p : partitions_)
        p.push_back(i_cell);
      continue;
    }

    // Find the tokens for bounding z-planes.
    int32_t lower_token = 0, upper_token = 0;
    double min_z, max_z;
    for (auto token : model::cells[i_cell]->rpn_) {
      if (token < OP_UNION) {
        const auto* surf = model::surfaces[std::abs(token) - 1].get();
        if (const auto* zplane = dynamic_cast<const SurfaceZPlane*>(surf)) {
          if (lower_token == 0 || zplane->z0_ < min_z) {
            lower_token = token;
            min_z = zplane->z0_;
          }
          if (upper_token == 0 || zplane->z0_ > max_z) {
            upper_token = token;
            max_z = zplane->z0_;
          }
        }
      }
    }

    // If there are no bounding z-planes, add this cell to all partitions.
    if (lower_token == 0) {
      for (auto& p : partitions_)
        p.push_back(i_cell);
      continue;
    }

    // Find the first partition this cell lies in.  If the lower_token indicates
    // a negative halfspace, then the cell is unbounded in the lower direction
    // and it lies in the first partition onward.  Otherwise, it is bounded by
    // the positive halfspace given by the lower_token.
    int first_partition = 0;
    if (lower_token > 0) {
      for (int i = 0; i < surfs_.size(); ++i) {
        if (lower_token == surfs_[i] + 1) {
          first_partition = i + 1;
          break;
        }
      }
    }

    // Find the last partition this cell lies in.  The logic is analogous to the
    // logic for first_partition.
    int last_partition = surfs_.size();
    if (upper_token < 0) {
      for (int i = first_partition; i < surfs_.size(); ++i) {
        if (upper_token == -(surfs_[i] + 1)) {
          last_partition = i;
          break;
        }
      }
    }

    // Add the cell to all relevant partitions.
    for (int i = first_partition; i <= last_partition; ++i) {
      partitions_[i].push_back(i_cell);
    }
  }
}

const vector<int32_t>& UniversePartitioner::get_cells(
  Position r, Direction u) const
{
  // Perform a binary search for the partition containing the given coordinates.
  int left = 0;
  int middle = (surfs_.size() - 1) / 2;
  int right = surfs_.size() - 1;
  while (true) {
    // Check the sense of the coordinates for the current surface.
    const auto& surf = *model::surfaces[surfs_[middle]];
    if (surf.sense(r, u)) {
      // The coordinates lie in the positive halfspace.  Recurse if there are
      // more surfaces to check.  Otherwise, return the cells on the positive
      // side of this surface.
      int right_leaf = right - (right - middle) / 2;
      if (right_leaf != middle) {
        left = middle + 1;
        middle = right_leaf;
      } else {
        return partitions_[middle + 1];
      }

    } else {
      // The coordinates lie in the negative halfspace.  Recurse if there are
      // more surfaces to check.  Otherwise, return the cells on the negative
      // side of this surface.
      int left_leaf = left + (middle - left) / 2;
      if (left_leaf != middle) {
        right = middle - 1;
        middle = left_leaf;
      } else {
        return partitions_[middle];
      }
    }
  }
}

//==============================================================================
// Non-method functions
//==============================================================================

void read_cells(pugi::xml_node node)
{
  // Count the number of cells.
  int n_cells = 0;
  for (pugi::xml_node cell_node : node.children("cell")) {
    n_cells++;
  }

  // Loop over XML cell elements and populate the array.
  model::cells.reserve(n_cells);
  for (pugi::xml_node cell_node : node.children("cell")) {
    model::cells.push_back(make_unique<CSGCell>(cell_node));
  }

  // Fill the cell map.
  for (int i = 0; i < model::cells.size(); i++) {
    int32_t id = model::cells[i]->id_;
    auto search = model::cell_map.find(id);
    if (search == model::cell_map.end()) {
      model::cell_map[id] = i;
    } else {
      fatal_error(
        fmt::format("Two or more cells use the same unique ID: {}", id));
    }
  }

  read_dagmc_universes(node);

  // Populate the Universe vector and map.
  for (int i = 0; i < model::cells.size(); i++) {
    int32_t uid = model::cells[i]->universe_;
    auto it = model::universe_map.find(uid);
    if (it == model::universe_map.end()) {
      model::universes.push_back(make_unique<Universe>());
      model::universes.back()->id_ = uid;
      model::universes.back()->cells_.push_back(i);
      model::universe_map[uid] = model::universes.size() - 1;
    } else {
      model::universes[it->second]->cells_.push_back(i);
    }
    if (model::cells[i]->virtual_lattice_) {
      model::universes[it->second]->filled_with_triso_base_=model::cells[i]->id_;
    }
  }
  model::universes.shrink_to_fit();

  // Allocate the cell overlap count if necessary.
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }

  if (model::cells.size() == 0) {
    fatal_error("No cells were found in the geometry.xml file");
  }
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int openmc_cell_get_fill(
  int32_t index, int* type, int32_t** indices, int32_t* n)
{
  if (index >= 0 && index < model::cells.size()) {
    Cell& c {*model::cells[index]};
    *type = static_cast<int>(c.type_);
    if (c.type_ == Fill::MATERIAL) {
      *indices = c.material_.data();
      *n = c.material_.size();
    } else {
      *indices = &c.fill_;
      *n = 1;
    }
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

extern "C" int openmc_cell_set_fill(
  int32_t index, int type, int32_t n, const int32_t* indices)
{
  Fill filltype = static_cast<Fill>(type);
  if (index >= 0 && index < model::cells.size()) {
    Cell& c {*model::cells[index]};
    if (filltype == Fill::MATERIAL) {
      c.type_ = Fill::MATERIAL;
      c.material_.clear();
      for (int i = 0; i < n; i++) {
        int i_mat = indices[i];
        if (i_mat == MATERIAL_VOID) {
          c.material_.push_back(MATERIAL_VOID);
        } else if (i_mat >= 0 && i_mat < model::materials.size()) {
          c.material_.push_back(i_mat);
        } else {
          set_errmsg("Index in materials array is out of bounds.");
          return OPENMC_E_OUT_OF_BOUNDS;
        }
      }
      c.material_.shrink_to_fit();
    } else if (filltype == Fill::UNIVERSE) {
      c.type_ = Fill::UNIVERSE;
    } else {
      c.type_ = Fill::LATTICE;
    }
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

extern "C" int openmc_cell_set_temperature(
  int32_t index, double T, const int32_t* instance, bool set_contained)
{
  if (index < 0 || index >= model::cells.size()) {
    strcpy(openmc_err_msg, "Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  int32_t instance_index = instance ? *instance : -1;
  try {
    model::cells[index]->set_temperature(T, instance_index, set_contained);
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_UNASSIGNED;
  }
  return 0;
}

extern "C" int openmc_cell_get_temperature(
  int32_t index, const int32_t* instance, double* T)
{
  if (index < 0 || index >= model::cells.size()) {
    strcpy(openmc_err_msg, "Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  int32_t instance_index = instance ? *instance : -1;
  try {
    *T = model::cells[index]->temperature(instance_index);
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_UNASSIGNED;
  }
  return 0;
}

//! Get the bounding box of a cell
extern "C" int openmc_cell_bounding_box(
  const int32_t index, double* llc, double* urc)
{

  BoundingBox bbox;

  const auto& c = model::cells[index];
  bbox = c->bounding_box();

  // set lower left corner values
  llc[0] = bbox.xmin;
  llc[1] = bbox.ymin;
  llc[2] = bbox.zmin;

  // set upper right corner values
  urc[0] = bbox.xmax;
  urc[1] = bbox.ymax;
  urc[2] = bbox.zmax;

  return 0;
}

//! Get the name of a cell
extern "C" int openmc_cell_get_name(int32_t index, const char** name)
{
  if (index < 0 || index >= model::cells.size()) {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *name = model::cells[index]->name().data();

  return 0;
}

//! Set the name of a cell
extern "C" int openmc_cell_set_name(int32_t index, const char* name)
{
  if (index < 0 || index >= model::cells.size()) {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  model::cells[index]->set_name(name);

  return 0;
}

//==============================================================================
//! Define a containing (parent) cell
//==============================================================================

//! Used to locate a universe fill in the geometry
struct ParentCell {
  bool operator==(const ParentCell& other) const
  {
    return cell_index == other.cell_index &&
           lattice_index == other.lattice_index;
  }

  bool operator<(const ParentCell& other) const
  {
    return cell_index < other.cell_index ||
           (cell_index == other.cell_index &&
             lattice_index < other.lattice_index);
  }

  gsl::index cell_index;
  gsl::index lattice_index;
};

//! Structure used to insert ParentCell into hashed STL data structures
struct ParentCellHash {
  std::size_t operator()(const ParentCell& p) const
  {
    return 4096 * p.cell_index + p.lattice_index;
  }
};

//! Used to manage a traversal stack when locating parent cells of a cell
//! instance in the model
struct ParentCellStack {

  //! push method that adds to the parent_cells visited cells for this search
  //! universe
  void push(int32_t search_universe, const ParentCell& pc)
  {
    parent_cells_.push_back(pc);
    // add parent cell to the set of cells we've visited for this search
    // universe
    visited_cells_[search_universe].insert(pc);
  }

  //! removes the last parent_cell and clears the visited cells for the popped
  //! cell's universe
  void pop()
  {
    visited_cells_[this->current_univ()].clear();
    parent_cells_.pop_back();
  }

  //! checks whether or not the parent cell has been visited already for this
  //! search universe
  bool visited(int32_t search_universe, const ParentCell& parent_cell)
  {
    return visited_cells_[search_universe].count(parent_cell) != 0;
  }

  //! return the next universe to search for a parent cell
  int32_t current_univ() const
  {
    return model::cells[parent_cells_.back().cell_index]->universe_;
  }

  //! indicates whether nor not parent cells are present on the stack
  bool empty() const { return parent_cells_.empty(); }

  //! compute an instance for the provided distribcell index
  int32_t compute_instance(int32_t distribcell_index) const
  {
    int32_t instance = 0;
    for (const auto& parent_cell : this->parent_cells_) {
      auto& cell = model::cells[parent_cell.cell_index];
      if (cell->type_ == Fill::UNIVERSE) {
        instance += cell->offset_[distribcell_index];
      } else if (cell->type_ == Fill::LATTICE) {
        auto& lattice = model::lattices[cell->fill_];
        instance +=
          lattice->offset(distribcell_index, parent_cell.lattice_index);
      }
    }
    return instance;
  }

  // Accessors
  vector<ParentCell>& parent_cells() { return parent_cells_; }
  const vector<ParentCell>& parent_cells() const { return parent_cells_; }

  // Data Members
  vector<ParentCell> parent_cells_;
  std::unordered_map<int32_t, std::unordered_set<ParentCell, ParentCellHash>>
    visited_cells_;
};

vector<ParentCell> Cell::find_parent_cells(
  int32_t instance, const Position& r) const {

  // create a temporary particle
  Particle dummy_particle {};
  dummy_particle.r() = r;
  dummy_particle.u() = {0., 0., 1.};

  return find_parent_cells(instance, dummy_particle);
}

vector<ParentCell> Cell::find_parent_cells(
  int32_t instance, Particle& p) const {
  // look up the particle's location
  exhaustive_find_cell(p);
  const auto& coords = p.coord();

  // build a parent cell stack from the particle coordinates
  ParentCellStack stack;
  bool cell_found = false;
  for (auto it = coords.begin(); it != coords.end(); it++) {
    const auto& coord = *it;
    const auto& cell = model::cells[coord.cell];
    // if the cell at this level matches the current cell, stop adding to the stack
    if (coord.cell == model::cell_map[this->id_]) {
      cell_found = true;
      break;
    }

    // if filled with a lattice, get the lattice index from the next
    // level in the coordinates to push to the stack
    int lattice_idx = C_NONE;
    if (cell->type_ == Fill::LATTICE) {
      const auto& next_coord = *(it + 1);
      lattice_idx = model::lattices[next_coord.lattice]->get_flat_index(next_coord.lattice_i);
    }
    stack.push(coord.universe, {coord.cell, lattice_idx});
  }

  // if this loop finished because the cell was found and
  // the instance matches the one requested in the call
  // we have the correct path and can return the stack
  if (cell_found && stack.compute_instance(this->distribcell_index_) == instance) {
    return stack.parent_cells();
  }

  // fall back on an exhaustive search for the cell's parents
  return exhaustive_find_parent_cells(instance);
}


vector<ParentCell> Cell::exhaustive_find_parent_cells(
  int32_t instance) const
{
  ParentCellStack stack;
  // start with this cell's universe
  int32_t prev_univ_idx;
  int32_t univ_idx = this->universe_;

  while (true) {
    const auto& univ = model::universes[univ_idx];
    prev_univ_idx = univ_idx;

    // search for a cell that is filled w/ this universe
    for (const auto& cell : model::cells) {
      // if this is a material-filled cell, move on
      if (cell->type_ == Fill::MATERIAL)
        continue;

      if (cell->type_ == Fill::UNIVERSE) {
        // if this is in the set of cells previously visited for this universe,
        // move on
        if (stack.visited(univ_idx, {model::cell_map[cell->id_], C_NONE}))
          continue;

        // if this cell contains the universe we're searching for, add it to the
        // stack
        if (cell->fill_ == univ_idx) {
          stack.push(univ_idx, {model::cell_map[cell->id_], C_NONE});
          univ_idx = cell->universe_;
        }
      } else if (cell->type_ == Fill::LATTICE) {
        // retrieve the lattice and lattice universes
        const auto& lattice = model::lattices[cell->fill_];
        const auto& lattice_univs = lattice->universes_;

        // start search for universe
        auto lat_it = lattice_univs.begin();
        while (true) {
          // find the next lattice cell with this universe
          lat_it = std::find(lat_it, lattice_univs.end(), univ_idx);
          if (lat_it == lattice_univs.end())
            break;

          int lattice_idx = lat_it - lattice_univs.begin();

          // move iterator forward one to avoid finding the same entry
          lat_it++;
          if (stack.visited(
                univ_idx, {model::cell_map[cell->id_], lattice_idx}))
            continue;

          // add this cell and lattice index to the stack and exit loop
          stack.push(univ_idx, {model::cell_map[cell->id_], lattice_idx});
          univ_idx = cell->universe_;
          break;
        }
      }
      // if we've updated the universe, break
      if (prev_univ_idx != univ_idx)
        break;
    } // end cell loop search for universe

    // if we're at the top of the geometry and the instance matches, we're done
    if (univ_idx == model::root_universe &&
        stack.compute_instance(this->distribcell_index_) == instance)
      break;

    // if there is no match on the original cell's universe, report an error
    if (univ_idx == this->universe_) {
      fatal_error(
        fmt::format("Could not find the parent cells for cell {}, instance {}.",
          this->id_, instance));
    }

    // if we don't find a suitable update, adjust the stack and continue
    if (univ_idx == model::root_universe || univ_idx == prev_univ_idx) {
      stack.pop();
      univ_idx = stack.empty() ? this->universe_ : stack.current_univ();
    }

  } // end while

  // reverse the stack so the highest cell comes first
  std::reverse(stack.parent_cells().begin(), stack.parent_cells().end());
  return stack.parent_cells();
}

std::unordered_map<int32_t, vector<int32_t>> Cell::get_contained_cells(
  int32_t instance, Position* hint) const
{
  std::unordered_map<int32_t, vector<int32_t>> contained_cells;

  // if this is a material-filled cell it has no contained cells
  if (this->type_ == Fill::MATERIAL)
    return contained_cells;

  // find the pathway through the geometry to this cell
  vector<ParentCell> parent_cells;

  // if a positional hint is provided, attempt to do a fast lookup
  // of the parent cells
  parent_cells = hint ? find_parent_cells(instance, *hint)
                      : exhaustive_find_parent_cells(instance);

  // if this cell is filled w/ a material, it contains no other cells
  if (type_ != Fill::MATERIAL) {
    this->get_contained_cells_inner(contained_cells, parent_cells);
  }

  return contained_cells;
}

//! Get all cells within this cell
void Cell::get_contained_cells_inner(
  std::unordered_map<int32_t, vector<int32_t>>& contained_cells,
  vector<ParentCell>& parent_cells) const
{

  // filled by material, determine instance based on parent cells
  if (type_ == Fill::MATERIAL) {
    int instance = 0;
    if (this->distribcell_index_ >= 0) {
      for (auto& parent_cell : parent_cells) {
        auto& cell = model::cells[parent_cell.cell_index];
        if (cell->type_ == Fill::UNIVERSE) {
          instance += cell->offset_[distribcell_index_];
        } else if (cell->type_ == Fill::LATTICE) {
          auto& lattice = model::lattices[cell->fill_];
          instance += lattice->offset(
            this->distribcell_index_, parent_cell.lattice_index);
        }
      }
    }
    // add entry to contained cells
    contained_cells[model::cell_map[id_]].push_back(instance);
    // filled with universe, add the containing cell to the parent cells
    // and recurse
  } else if (type_ == Fill::UNIVERSE) {
    parent_cells.push_back({model::cell_map[id_], -1});
    auto& univ = model::universes[fill_];
    for (auto cell_index : univ->cells_) {
      auto& cell = model::cells[cell_index];
      cell->get_contained_cells_inner(contained_cells, parent_cells);
    }
    parent_cells.pop_back();
    // filled with a lattice, visit each universe in the lattice
    // with a recursive call to collect the cell instances
  } else if (type_ == Fill::LATTICE) {
    auto& lattice = model::lattices[fill_];
    for (auto i = lattice->begin(); i != lattice->end(); ++i) {
      auto& univ = model::universes[*i];
      parent_cells.push_back({model::cell_map[id_], i.indx_});
      for (auto cell_index : univ->cells_) {
        auto& cell = model::cells[cell_index];
        cell->get_contained_cells_inner(contained_cells, parent_cells);
      }
      parent_cells.pop_back();
    }
  }
}

//! Return the index in the cells array of a cell with a given ID
extern "C" int openmc_get_cell_index(int32_t id, int32_t* index)
{
  auto it = model::cell_map.find(id);
  if (it != model::cell_map.end()) {
    *index = it->second;
    return 0;
  } else {
    set_errmsg("No cell exists with ID=" + std::to_string(id) + ".");
    return OPENMC_E_INVALID_ID;
  }
}

//! Return the ID of a cell
extern "C" int openmc_cell_get_id(int32_t index, int32_t* id)
{
  if (index >= 0 && index < model::cells.size()) {
    *id = model::cells[index]->id_;
    return 0;
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

//! Set the ID of a cell
extern "C" int openmc_cell_set_id(int32_t index, int32_t id)
{
  if (index >= 0 && index < model::cells.size()) {
    model::cells[index]->id_ = id;
    model::cell_map[id] = index;
    return 0;
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

//! Return the translation vector of a cell
extern "C" int openmc_cell_get_translation(int32_t index, double xyz[])
{
  if (index >= 0 && index < model::cells.size()) {
    auto& cell = model::cells[index];
    xyz[0] = cell->translation_.x;
    xyz[1] = cell->translation_.y;
    xyz[2] = cell->translation_.z;
    return 0;
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

//! Set the translation vector of a cell
extern "C" int openmc_cell_set_translation(int32_t index, const double xyz[])
{
  if (index >= 0 && index < model::cells.size()) {
    if (model::cells[index]->fill_ == C_NONE) {
      set_errmsg(fmt::format("Cannot apply a translation to cell {}"
                             " because it is not filled with another universe",
        index));
      return OPENMC_E_GEOMETRY;
    }
    model::cells[index]->translation_ = Position(xyz);
    return 0;
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

//! Return the rotation matrix of a cell
extern "C" int openmc_cell_get_rotation(int32_t index, double rot[], size_t* n)
{
  if (index >= 0 && index < model::cells.size()) {
    auto& cell = model::cells[index];
    *n = cell->rotation_.size();
    std::memcpy(rot, cell->rotation_.data(), *n * sizeof(cell->rotation_[0]));
    return 0;
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

//! Set the flattened rotation matrix of a cell
extern "C" int openmc_cell_set_rotation(
  int32_t index, const double rot[], size_t rot_len)
{
  if (index >= 0 && index < model::cells.size()) {
    if (model::cells[index]->fill_ == C_NONE) {
      set_errmsg(fmt::format("Cannot apply a rotation to cell {}"
                             " because it is not filled with another universe",
        index));
      return OPENMC_E_GEOMETRY;
    }
    std::vector<double> vec_rot(rot, rot + rot_len);
    model::cells[index]->set_rotation(vec_rot);
    return 0;
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

//! Get the number of instances of the requested cell
extern "C" int openmc_cell_get_num_instances(
  int32_t index, int32_t* num_instances)
{
  if (index < 0 || index >= model::cells.size()) {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  *num_instances = model::cells[index]->n_instances_;
  return 0;
}

//! Extend the cells array by n elements
extern "C" int openmc_extend_cells(
  int32_t n, int32_t* index_start, int32_t* index_end)
{
  if (index_start)
    *index_start = model::cells.size();
  if (index_end)
    *index_end = model::cells.size() + n - 1;
  for (int32_t i = 0; i < n; i++) {
    model::cells.push_back(make_unique<CSGCell>());
  }
  return 0;
}

extern "C" int cells_size()
{
  return model::cells.size();
}

} // namespace openmc
