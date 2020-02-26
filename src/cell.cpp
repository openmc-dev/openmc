
#include "openmc/cell.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iterator>
#include <sstream>
#include <set>
#include <string>

#include <fmt/core.h>
#include <gsl/gsl>

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
#include "openmc/surface.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  std::vector<std::unique_ptr<Cell>> cells;
  std::unordered_map<int32_t, int32_t> cell_map;

  std::vector<std::unique_ptr<Universe>> universes;
  std::unordered_map<int32_t, int32_t> universe_map;
} // namespace model

//==============================================================================
//! Convert region specification string to integer tokens.
//!
//! The characters (, ), |, and ~ count as separate tokens since they represent
//! operators.
//==============================================================================

std::vector<int32_t>
tokenize(const std::string region_spec) {
  // Check for an empty region_spec first.
  std::vector<int32_t> tokens;
  if (region_spec.empty()) {
    return tokens;
  }

  // Parse all halfspaces and operators except for intersection (whitespace).
  for (int i = 0; i < region_spec.size(); ) {
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

    } else if (region_spec[i] == '-' || region_spec[i] == '+'
               || std::isdigit(region_spec[i])) {
      // This is the start of a halfspace specification.  Iterate j until we
      // find the end, then push-back everything between i and j.
      int j = i + 1;
      while (j < region_spec.size() && std::isdigit(region_spec[j])) {j++;}
      tokens.push_back(std::stoi(region_spec.substr(i, j-i)));
      i = j;

    } else if (std::isspace(region_spec[i])) {
      i++;

    } else {
      auto err_msg = fmt::format(
        "Region specification contains invalid character, \"{}\"", region_spec[i]);
      fatal_error(err_msg);
    }
  }

  // Add in intersection operators where a missing operator is needed.
  int i = 0;
  while (i < tokens.size()-1) {
    bool left_compat {(tokens[i] < OP_UNION) || (tokens[i] == OP_RIGHT_PAREN)};
    bool right_compat {(tokens[i+1] < OP_UNION)
                       || (tokens[i+1] == OP_LEFT_PAREN)
                       || (tokens[i+1] == OP_COMPLEMENT)};
    if (left_compat && right_compat) {
      tokens.insert(tokens.begin()+i+1, OP_INTERSECTION);
    }
    i++;
  }

  return tokens;
}

//==============================================================================
//! Convert infix region specification to Reverse Polish Notation (RPN)
//!
//! This function uses the shunting-yard algorithm.
//==============================================================================

std::vector<int32_t>
generate_rpn(int32_t cell_id, std::vector<int32_t> infix)
{
  std::vector<int32_t> rpn;
  std::vector<int32_t> stack;

  for (int32_t token : infix) {
    if (token < OP_UNION) {
      // If token is not an operator, add it to output
      rpn.push_back(token);
    } else if (token < OP_RIGHT_PAREN) {
      // Regular operators union, intersection, complement
      while (stack.size() > 0) {
        int32_t op = stack.back();

        if (op < OP_RIGHT_PAREN &&
             ((token == OP_COMPLEMENT && token < op) ||
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
            "Mismatched parentheses in region specification for cell {}", cell_id));
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

void
Universe::to_hdf5(hid_t universes_group) const
{
  // Create a group for this universe.
  auto group = create_group(universes_group, fmt::format("universe {}", id_));

  // Write the contained cells.
  if (cells_.size() > 0) {
    std::vector<int32_t> cell_ids;
    for (auto i_cell : cells_) cell_ids.push_back(model::cells[i_cell]->id_);
    write_dataset(group, "cells", cell_ids);
  }

  close_group(group);
}

BoundingBox Universe::bounding_box() const {
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

double
Cell::temperature(int32_t instance) const
{
  if (sqrtkT_.size() < 1) {
    throw std::runtime_error{"Cell temperature has not yet been set."};
  }

  if (instance >= 0) {
    double sqrtkT = sqrtkT_.size() == 1 ?
      sqrtkT_.at(0) :
      sqrtkT_.at(instance);
    return sqrtkT * sqrtkT / K_BOLTZMANN;
  } else {
    return sqrtkT_[0] * sqrtkT_[0] / K_BOLTZMANN;
  }
}

void
Cell::set_temperature(double T, int32_t instance)
{
  if (settings::temperature_method == TemperatureMethod::INTERPOLATION) {
    if (T < data::temperature_min) {
      throw std::runtime_error{"Temperature is below minimum temperature at "
        "which data is available."};
    } else if (T > data::temperature_max) {
      throw std::runtime_error{"Temperature is above maximum temperature at "
        "which data is available."};
    }
  }

  if (instance >= 0) {
    // If temperature vector is not big enough, resize it first
    if (sqrtkT_.size() != n_instances_) sqrtkT_.resize(n_instances_, sqrtkT_[0]);

    // Set temperature for the corresponding instance
    sqrtkT_.at(instance) = std::sqrt(K_BOLTZMANN * T);
  } else {
    // Set temperature for all instances
    for (auto& T_ : sqrtkT_) {
      T_ = std::sqrt(K_BOLTZMANN * T);
    }
  }
}

//==============================================================================
// CSGCell implementation
//==============================================================================

CSGCell::CSGCell() {} // empty constructor

CSGCell::CSGCell(pugi::xml_node cell_node)
{
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

  // Make sure that either material or fill was specified, but not both.
  bool fill_present = check_for_node(cell_node, "fill");
  bool material_present = check_for_node(cell_node, "material");
  if (!(fill_present || material_present)) {
    fatal_error(fmt::format(
      "Neither material nor fill was specified for cell {}", id_));
  }
  if (fill_present && material_present) {
    fatal_error(fmt::format("Cell {} has both a material and a fill specified; "
      "only one can be specified per cell", id_));
  }

  if (fill_present) {
    fill_ = std::stoi(get_node_value(cell_node, "fill"));
    if (fill_ == universe_) {
      fatal_error(fmt::format("Cell {} is filled with the same universe that"
        "it is contained in.", id_));
    }
  } else {
    fill_ = C_NONE;
  }

  // Read the material element.  There can be zero materials (filled with a
  // universe), more than one material (distribmats), and some materials may
  // be "void".
  if (material_present) {
    std::vector<std::string> mats
         {get_node_array<std::string>(cell_node, "material", true)};
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
      fatal_error(fmt::format("An empty material element was specified for cell {}",
        id_));
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
        "specification is only valid for cells filled with a material.", id_));
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
        throw std::runtime_error{"Invalid surface ID " + std::to_string(abs(r))
          + " specified in region for cell " + std::to_string(id_) + "."};
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

  // Read the translation vector.
  if (check_for_node(cell_node, "translation")) {
    if (fill_ == C_NONE) {
      fatal_error(fmt::format("Cannot apply a translation to cell {}"
        " because it is not filled with another universe", id_));
    }

    auto xyz {get_node_array<double>(cell_node, "translation")};
    if (xyz.size() != 3) {
      fatal_error(fmt::format(
        "Non-3D translation vector applied to cell {}", id_));
    }
    translation_ = xyz;
  }

  // Read the rotation transform.
  if (check_for_node(cell_node, "rotation")) {
    if (fill_ == C_NONE) {
      fatal_error(fmt::format("Cannot apply a rotation to cell {}"
        " because it is not filled with another universe", id_));
    }

    auto rot {get_node_array<double>(cell_node, "rotation")};
    if (rot.size() != 3 && rot.size() != 9) {
      fatal_error(fmt::format(
        "Non-3D rotation vector applied to cell {}", id_));
    }

    // Compute and store the rotation matrix.
    rotation_.reserve(rot.size() == 9 ? 9 : 12);
    if (rot.size() == 3) {
      double phi = -rot[0] * PI / 180.0;
      double theta = -rot[1] * PI / 180.0;
      double psi = -rot[2] * PI / 180.0;
      rotation_.push_back(std::cos(theta) * std::cos(psi));
      rotation_.push_back(-std::cos(phi) * std::sin(psi)
                          + std::sin(phi) * std::sin(theta) * std::cos(psi));
      rotation_.push_back(std::sin(phi) * std::sin(psi)
                          + std::cos(phi) * std::sin(theta) * std::cos(psi));
      rotation_.push_back(std::cos(theta) * std::sin(psi));
      rotation_.push_back(std::cos(phi) * std::cos(psi)
                          + std::sin(phi) * std::sin(theta) * std::sin(psi));
      rotation_.push_back(-std::sin(phi) * std::cos(psi)
                          + std::cos(phi) * std::sin(theta) * std::sin(psi));
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
}

//==============================================================================

bool
CSGCell::contains(Position r, Direction u, int32_t on_surface) const
{
  if (simple_) {
    return contains_simple(r, u, on_surface);
  } else {
    return contains_complex(r, u, on_surface);
  }
}

//==============================================================================

std::pair<double, int32_t>
CSGCell::distance(Position r, Direction u, int32_t on_surface, Particle* p) const
{
  double min_dist {INFTY};
  int32_t i_surf {std::numeric_limits<int32_t>::max()};

  for (int32_t token : rpn_) {
    // Ignore this token if it corresponds to an operator rather than a region.
    if (token >= OP_UNION) continue;

    // Calculate the distance to this surface.
    // Note the off-by-one indexing
    bool coincident {std::abs(token) == std::abs(on_surface)};
    double d {model::surfaces[abs(token)-1]->distance(r, u, coincident)};

    // Check if this distance is the new minimum.
    if (d < min_dist) {
      if (std::abs(d - min_dist) / min_dist >= FP_PRECISION) {
        min_dist = d;
        i_surf = -token;
      }
    }
  }

  return {min_dist, i_surf};
}

//==============================================================================

void
CSGCell::to_hdf5(hid_t cell_group) const
{
  // Create a group for this cell.
  auto group = create_group(cell_group, fmt::format("cell {}", id_));

  if (!name_.empty()) {
    write_string(group, "name", name_, false);
  }

  write_dataset(group, "universe", model::universes[universe_]->id_);

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
        auto surf_id = model::surfaces[abs(token)-1]->id_;
        region_spec << " " << ((token > 0) ? surf_id : -surf_id);
      }
    }
    write_string(group, "region", region_spec.str(), false);
  }

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

BoundingBox CSGCell::bounding_box_simple() const {
  BoundingBox bbox;
  for (int32_t token : rpn_) {
    bbox &= model::surfaces[abs(token)-1]->bounding_box(token > 0);
  }
  return bbox;
}

void CSGCell::apply_demorgan(std::vector<int32_t>::iterator start,
                             std::vector<int32_t>::iterator stop)
{
  while (start < stop) {
    if (*start < OP_UNION) { *start *= -1; }
    else if (*start == OP_UNION) { *start = OP_INTERSECTION; }
    else if (*start == OP_INTERSECTION) { *start = OP_UNION; }
    start++;
  }
}

std::vector<int32_t>::iterator
CSGCell::find_left_parenthesis(std::vector<int32_t>::iterator start,
                               const std::vector<int32_t>& rpn) {
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

void CSGCell::remove_complement_ops(std::vector<int32_t>& rpn) {
  auto it = std::find(rpn.begin(), rpn.end(), OP_COMPLEMENT);
  while (it != rpn.end()) {
    // find the opening parenthesis (if any)
    auto left = find_left_parenthesis(it, rpn);
    std::vector<int32_t> tmp(left, it+1);

    // apply DeMorgan's law to any surfaces/operators between these
    // positions in the RPN
    apply_demorgan(left, it);
    // remove complement operator
    rpn.erase(it);
    // update iterator position
    it = std::find(rpn.begin(), rpn.end(), OP_COMPLEMENT);
  }
}

BoundingBox CSGCell::bounding_box_complex(std::vector<int32_t> rpn) {
  // remove complements by adjusting surface signs and operators
  remove_complement_ops(rpn);

  std::vector<BoundingBox> stack(rpn.size());
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

BoundingBox CSGCell::bounding_box() const {
  return simple_ ? bounding_box_simple() : bounding_box_complex(rpn_);
}

//==============================================================================

bool
CSGCell::contains_simple(Position r, Direction u, int32_t on_surface) const
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
      bool sense = model::surfaces[abs(token)-1]->sense(r, u);
      if (sense != (token > 0)) {return false;}
    }
  }
  return true;
}

//==============================================================================

bool
CSGCell::contains_complex(Position r, Direction u, int32_t on_surface) const
{
  // Make a stack of booleans.  We don't know how big it needs to be, but we do
  // know that rpn.size() is an upper-bound.
  std::vector<bool> stack(rpn_.size());
  int i_stack = -1;

  for (int32_t token : rpn_) {
    // If the token is a binary operator (intersection/union), apply it to
    // the last two items on the stack. If the token is a unary operator
    // (complement), apply it to the last item on the stack.
    if (token == OP_UNION) {
      stack[i_stack-1] = stack[i_stack-1] || stack[i_stack];
      i_stack --;
    } else if (token == OP_INTERSECTION) {
      stack[i_stack-1] = stack[i_stack-1] && stack[i_stack];
      i_stack --;
    } else if (token == OP_COMPLEMENT) {
      stack[i_stack] = !stack[i_stack];
    } else {
      // If the token is not an operator, evaluate the sense of particle with
      // respect to the surface and see if the token matches the sense. If the
      // particle's surface attribute is set and matches the token, that
      // overrides the determination based on sense().
      i_stack ++;
      if (token == on_surface) {
        stack[i_stack] = true;
      } else if (-token == on_surface) {
        stack[i_stack] = false;
      } else {
        // Note the off-by-one indexing
        bool sense = model::surfaces[abs(token)-1]->sense(r, u);
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
// DAGMC Cell implementation
//==============================================================================
#ifdef DAGMC
DAGCell::DAGCell() : Cell{} {};

std::pair<double, int32_t>
DAGCell::distance(Position r, Direction u, int32_t on_surface, Particle* p) const
{
  // if we've changed direction or we're not on a surface,
  // reset the history and update last direction
  if (u != p->last_dir_ || on_surface == 0) {
    p->history_.reset();
    p->last_dir_ = u;
  }

  moab::ErrorCode rval;
  moab::EntityHandle vol = dagmc_ptr_->entity_by_index(3, dag_index_);
  moab::EntityHandle hit_surf;
  double dist;
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  rval = dagmc_ptr_->ray_fire(vol, pnt, dir, hit_surf, dist, &p->history_);
  MB_CHK_ERR_CONT(rval);
  int surf_idx;
  if (hit_surf != 0) {
    surf_idx = dagmc_ptr_->index_by_handle(hit_surf);
  } else {
    // indicate that particle is lost
    surf_idx = -1;
    dist = INFINITY;
  }

  return {dist, surf_idx};
}

bool DAGCell::contains(Position r, Direction u, int32_t on_surface) const
{
  moab::ErrorCode rval;
  moab::EntityHandle vol = dagmc_ptr_->entity_by_index(3, dag_index_);

  int result = 0;
  double pnt[3] = {r.x, r.y, r.z};
  double dir[3] = {u.x, u.y, u.z};
  rval = dagmc_ptr_->point_in_volume(vol, pnt, result, dir);
  MB_CHK_ERR_CONT(rval);
  return result;
}

void DAGCell::to_hdf5(hid_t group_id) const { return; }

BoundingBox DAGCell::bounding_box() const
{
  moab::ErrorCode rval;
  moab::EntityHandle vol = dagmc_ptr_->entity_by_index(3, dag_index_);
  double min[3], max[3];
  rval = dagmc_ptr_->getobb(vol, min, max);
  MB_CHK_ERR_CONT(rval);
  return {min[0], max[0], min[1], max[1], min[2], max[2]};
}

#endif

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
      for (auto& p : partitions_) p.push_back(i_cell);
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
      for (auto& p : partitions_) p.push_back(i_cell);
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

const std::vector<int32_t>&
UniversePartitioner::get_cells(Position r, Direction u) const
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
        return partitions_[middle+1];
      }

    } else {
      // The coordinates lie in the negative halfspace.  Recurse if there are
      // more surfaces to check.  Otherwise, return the cells on the negative
      // side of this surface.
      int left_leaf = left + (middle - left) / 2;
      if (left_leaf != middle) {
        right = middle-1;
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
  for (pugi::xml_node cell_node: node.children("cell")) {n_cells++;}
  if (n_cells == 0) {
    fatal_error("No cells found in geometry.xml!");
  }

  // Loop over XML cell elements and populate the array.
  model::cells.reserve(n_cells);
  for (pugi::xml_node cell_node : node.children("cell")) {
    model::cells.push_back(std::make_unique<CSGCell>(cell_node));
  }

  // Fill the cell map.
  for (int i = 0; i < model::cells.size(); i++) {
    int32_t id = model::cells[i]->id_;
    auto search = model::cell_map.find(id);
    if (search == model::cell_map.end()) {
      model::cell_map[id] = i;
    } else {
      fatal_error(fmt::format("Two or more cells use the same unique ID: {}", id));
    }
  }

  // Populate the Universe vector and map.
  for (int i = 0; i < model::cells.size(); i++) {
    int32_t uid = model::cells[i]->universe_;
    auto it = model::universe_map.find(uid);
    if (it == model::universe_map.end()) {
      model::universes.push_back(std::make_unique<Universe>());
      model::universes.back()->id_ = uid;
      model::universes.back()->cells_.push_back(i);
      model::universe_map[uid] = model::universes.size() - 1;
    } else {
      model::universes[it->second]->cells_.push_back(i);
    }
  }
  model::universes.shrink_to_fit();

  // Allocate the cell overlap count if necessary.
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_cell_get_fill(int32_t index, int* type, int32_t** indices, int32_t* n)
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

extern "C" int
openmc_cell_set_fill(int32_t index, int type, int32_t n,
                     const int32_t* indices)
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

extern "C" int
openmc_cell_set_temperature(int32_t index, double T, const int32_t* instance)
{
  if (index < 0 || index >= model::cells.size()) {
    strcpy(openmc_err_msg, "Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  int32_t instance_index = instance ? *instance : -1;
  try {
    model::cells[index]->set_temperature(T, instance_index);
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_UNASSIGNED;
  }
  return 0;
}

extern "C" int
openmc_cell_get_temperature(int32_t index, const int32_t* instance, double* T)
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
extern "C" int
openmc_cell_bounding_box(const int32_t index, double* llc, double* urc) {

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
extern "C" int
openmc_cell_get_name(int32_t index, const char** name) {
  if (index < 0 || index >= model::cells.size()) {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *name = model::cells[index]->name().data();

  return 0;
}

//! Set the name of a cell
extern "C" int
openmc_cell_set_name(int32_t index, const char* name) {
  if (index < 0 || index >= model::cells.size()) {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  model::cells[index]->set_name(name);

  return 0;
}


//! Return the index in the cells array of a cell with a given ID
extern "C" int
openmc_get_cell_index(int32_t id, int32_t* index)
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
extern "C" int
openmc_cell_get_id(int32_t index, int32_t* id)
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
extern "C" int
openmc_cell_set_id(int32_t index, int32_t id)
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

//! Extend the cells array by n elements
extern "C" int
openmc_extend_cells(int32_t n, int32_t* index_start, int32_t* index_end)
{
  if (index_start) *index_start = model::cells.size();
  if (index_end) *index_end = model::cells.size() + n - 1;
  for (int32_t i = 0; i < n; i++) {
    model::cells.push_back(std::make_unique<CSGCell>());
  }
  return 0;
}

#ifdef DAGMC
int32_t next_cell(DAGCell* cur_cell, DAGSurface* surf_xed)
{
  moab::EntityHandle surf =
    surf_xed->dagmc_ptr_->entity_by_index(2, surf_xed->dag_index_);
  moab::EntityHandle vol =
    cur_cell->dagmc_ptr_->entity_by_index(3, cur_cell->dag_index_);

  moab::EntityHandle new_vol;
  cur_cell->dagmc_ptr_->next_vol(surf, vol, new_vol);

  return cur_cell->dagmc_ptr_->index_by_handle(new_vol);
}
#endif

extern "C" int cells_size() { return model::cells.size(); }

} // namespace openmc
