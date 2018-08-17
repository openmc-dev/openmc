#include "openmc/cell.h"

#include <cmath>
#include <sstream>
#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/settings.h"
#include "openmc/surface.h"
#include "openmc/xml_interface.h"

//TODO: remove this include
#include <iostream>

//extern "C" {int32_t n_cells {0};}

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

int32_t n_cells {0};

std::vector<Cell*> cells;
std::unordered_map<int32_t, int32_t> cell_map;

std::vector<Universe*> universes;
std::unordered_map<int32_t, int32_t> universe_map;

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
      std::stringstream err_msg;
      err_msg << "Region specification contains invalid character, \""
              << region_spec[i] << "\"";
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
          std::stringstream err_msg;
          err_msg << "Mismatched parentheses in region specification for cell "
                  << cell_id;
          fatal_error(err_msg);
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
      std::stringstream err_msg;
      err_msg << "Mismatched parentheses in region specification for cell "
              << cell_id;
      fatal_error(err_msg);
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
  std::stringstream group_name;
  group_name << "universe " << id_;
  auto group = create_group(universes_group, group_name);

  // Write the contained cells.
  if (cells_.size() > 0) {
    std::vector<int32_t> cell_ids;
    for (auto i_cell : cells_) cell_ids.push_back(cells[i_cell]->id_);
    write_dataset(group, "cells", cell_ids);
  }

  close_group(group);
}

//==============================================================================
// Cell implementation
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
    std::stringstream err_msg;
    err_msg << "Neither material nor fill was specified for cell " << id_;
    fatal_error(err_msg);
  }
  if (fill_present && material_present) {
    std::stringstream err_msg;
    err_msg << "Cell " << id_ << " has both a material and a fill specified; "
            << "only one can be specified per cell";
    fatal_error(err_msg);
  }

  if (fill_present) {
    fill_ = std::stoi(get_node_value(cell_node, "fill"));
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
      std::stringstream err_msg;
      err_msg << "An empty material element was specified for cell " << id_;
      fatal_error(err_msg);
    }
  }

  // Read the temperature element which may be distributed like materials.
  if (check_for_node(cell_node, "temperature")) {
    sqrtkT_ = get_node_array<double>(cell_node, "temperature");
    sqrtkT_.shrink_to_fit();

    // Make sure this is a material-filled cell.
    if (material_.size() == 0) {
      std::stringstream err_msg;
      err_msg << "Cell " << id_ << " was specified with a temperature but "
           "no material. Temperature specification is only valid for cells "
           "filled with a material.";
      fatal_error(err_msg);
    }

    // Make sure all temperatures are non-negative.
    for (auto T : sqrtkT_) {
      if (T < 0) {
        std::stringstream err_msg;
        err_msg << "Cell " << id_
                << " was specified with a negative temperature";
        fatal_error(err_msg);
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
      r = copysign(surface_map[abs(r)] + 1, r);
    }
  }

  // Convert the infix region spec to RPN.
  rpn_ = generate_rpn(id_, region_);
  rpn_.shrink_to_fit();

  // Check if this is a simple cell.
  simple_ = true;
  for (int32_t token : rpn_) {
    if ((token == OP_COMPLEMENT) || (token == OP_UNION)) {
      simple_ = false;
      break;
    }
  }

  // Read the translation vector.
  if (check_for_node(cell_node, "translation")) {
    if (fill_ == C_NONE) {
      std::stringstream err_msg;
      err_msg << "Cannot apply a translation to cell " << id_
              << " because it is not filled with another universe";
      fatal_error(err_msg);
    }

    auto xyz {get_node_array<double>(cell_node, "translation")};
    if (xyz.size() != 3) {
      std::stringstream err_msg;
      err_msg << "Non-3D translation vector applied to cell " << id_;
      fatal_error(err_msg);
    }
    translation_ = xyz;
  }

  // Read the rotation transform.
  if (check_for_node(cell_node, "rotation")) {
    if (fill_ == C_NONE) {
      std::stringstream err_msg;
      err_msg << "Cannot apply a rotation to cell " << id_
              << " because it is not filled with another universe";
      fatal_error(err_msg);
    }

    auto rot {get_node_array<double>(cell_node, "rotation")};
    if (rot.size() != 3) {
      std::stringstream err_msg;
      err_msg << "Non-3D rotation vector applied to cell " << id_;
      fatal_error(err_msg);
    }

    // Store the rotation angles.
    rotation_.reserve(12);
    rotation_.push_back(rot[0]);
    rotation_.push_back(rot[1]);
    rotation_.push_back(rot[2]);

    // Compute and store the rotation matrix.
    auto phi = -rot[0] * PI / 180.0;
    auto theta = -rot[1] * PI / 180.0;
    auto psi = -rot[2] * PI / 180.0;
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
CSGCell::distance(Position r, Direction u, int32_t on_surface) const
{
  double min_dist {INFTY};
  int32_t i_surf {std::numeric_limits<int32_t>::max()};

  for (int32_t token : rpn_) {
    // Ignore this token if it corresponds to an operator rather than a region.
    if (token >= OP_UNION) continue;

    // Calculate the distance to this surface.
    // Note the off-by-one indexing
    bool coincident {token == on_surface};
    double d {surfaces[abs(token)-1]->distance(r, u, coincident)};

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
  std::stringstream group_name;
  group_name << "cell " << id_;
  auto group = create_group(cells_group, group_name);

  if (!name_.empty()) {
    write_string(group, "name", name_, false);
  }

  write_dataset(group, "universe", universes[universe_]->id_);

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
        region_spec << " "
             << copysign(surfaces[abs(token)-1]->id_, token);
      }
    }
    write_string(group, "region", region_spec.str(), false);
  }

  // Write fill information.
  if (type_ == FILL_MATERIAL) {
    write_dataset(group, "fill_type", "material");
    std::vector<int32_t> mat_ids;
    for (auto i_mat : material_) {
      if (i_mat != MATERIAL_VOID) {
        mat_ids.push_back(materials[i_mat]->id);
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

  } else if (type_ == FILL_UNIVERSE) {
    write_dataset(group, "fill_type", "universe");
    write_dataset(group, "fill", universes[fill_]->id_);
    if (translation_ != Position(0, 0, 0)) {
      write_dataset(group, "translation", translation_);
    }
    if (!rotation_.empty()) {
      std::array<double, 3> rot {rotation_[0], rotation_[1], rotation_[2]};
      write_dataset(group, "rotation", rot);
    }

  } else if (type_ == FILL_LATTICE) {
    write_dataset(group, "fill_type", "lattice");
    write_dataset(group, "lattice", lattices[fill_]->id_);
  }

  close_group(group);
}

//==============================================================================

bool
CSGCell::contains_simple(Position r, Direction u, int32_t on_surface) const
{
  for (int32_t token : rpn_) {
    if (token < OP_UNION) {
      // If the token is not an operator, evaluate the sense of particle with
      // respect to the surface and see if the token matches the sense. If the
      // particle's surface attribute is set and matches the token, that
      // overrides the determination based on sense().
      if (token == on_surface) {
      } else if (-token == on_surface) {
        return false;
      } else {
        // Note the off-by-one indexing
        bool sense = surfaces[abs(token)-1]->sense(r, u);
        if (sense != (token > 0)) {return false;}
      }
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
  bool stack[rpn_.size()];
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
        bool sense = surfaces[abs(token)-1]->sense(r, u);
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
// CAD Cell implementation
//==============================================================================
#ifdef CAD
CADCell::CADCell() : Cell{} {};

std::pair<double, int32_t> CADCell::distance(Position p, Direction u, int32_t on_surface) const {

  moab::EntityHandle vol = dagmc_ptr->entity_by_id(3, id);
  moab::EntityHandle hit_surf;
  double dist;
  dagmc_ptr->ray_fire(vol, p.xyz, u.xyz, hit_surf, dist);

  int surf_idx;
  if(hit_surf != 0) {
    surf_idx = dagmc_ptr->index_by_handle(hit_surf);
  }
  // indicate that particle is lost
  else {
    surf_idx = -1;
  }
  
  std::pair<double, int> result(dist, surf_idx);
  
  return result;
}
  
bool CADCell::contains(Position p, Direction u, int32_t on_surface) const {

  moab::EntityHandle vol = dagmc_ptr->entity_by_id(3, id);

  int result = 0;
  dagmc_ptr->point_in_volume(vol, p.xyz, result, u.xyz);
  
  return bool(result);
}

void CADCell::to_hdf5(hid_t group_id) const { return; }
#endif

//==============================================================================
// Non-method functions
extern "C" void
read_cells(pugi::xml_node* node)
{
  // Count the number of cells.
  for (pugi::xml_node cell_node: node->children("cell")) {n_cells++;}
  if (n_cells == 0) {
    fatal_error("No cells found in geometry.xml!");
  }

  // Loop over XML cell elements and populate the array.
  cells.reserve(n_cells);
  for (pugi::xml_node cell_node: node->children("cell")) {
    global_cells.push_back(new CSGCell(cell_node));
  }

  // Populate the Universe vector and map.
  for (int i = 0; i < cells.size(); i++) {
    int32_t uid = cells[i]->universe_;
    auto it = universe_map.find(uid);
    if (it == universe_map.end()) {
      universes.push_back(new Universe());
      universes.back()->id_ = uid;
      universes.back()->cells_.push_back(i);
      universe_map[uid] = universes.size() - 1;
    } else {
      universes[it->second]->cells_.push_back(i);
    }
  }
  universes.shrink_to_fit();

  // Allocate the cell overlap count if necessary.
  if (settings::check_overlaps) overlap_check_count.resize(n_cells, 0);
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_cell_get_fill(int32_t index, int* type, int32_t** indices, int32_t* n)
{
  if (index >= 1 && index <= cells.size()) {
    //TODO: off-by-one
    Cell& c {*cells[index - 1]};
    *type = c.type_;
    if (c.type_ == FILL_MATERIAL) {
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
  if (index >= 1 && index <= cells.size()) {
    //TODO: off-by-one
    Cell& c {*cells[index - 1]};
    if (type == FILL_MATERIAL) {
      c.type_ = FILL_MATERIAL;
      c.material_.clear();
      for (int i = 0; i < n; i++) {
        int i_mat = indices[i];
        if (i_mat == MATERIAL_VOID) {
          c.material_.push_back(MATERIAL_VOID);
        } else if (i_mat >= 1 && i_mat <= materials.size()) {
          //TODO: off-by-one
          c.material_.push_back(i_mat - 1);
        } else {
          set_errmsg("Index in materials array is out of bounds.");
          return OPENMC_E_OUT_OF_BOUNDS;
        }
      }
      c.material_.shrink_to_fit();
    } else if (type == FILL_UNIVERSE) {
      c.type_ = FILL_UNIVERSE;
    } else {
      c.type_ = FILL_LATTICE;
    }
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

//TODO: make sure data is loaded for this temperature
extern "C" int
openmc_cell_set_temperature(int32_t index, double T, const int32_t* instance)
{
  if (index >= 1 && index <= cells.size()) {
    //TODO: off-by-one
    Cell& c {*cells[index - 1]};

    if (instance) {
      if (*instance >= 0 && *instance < c.sqrtkT_.size()) {
        c.sqrtkT_[*instance] = std::sqrt(K_BOLTZMANN * T);
      } else {
        strcpy(openmc_err_msg, "Distribcell instance is out of bounds.");
        return OPENMC_E_OUT_OF_BOUNDS;
      }
    } else {
      for (auto& T_ : c.sqrtkT_) {
        T_ = std::sqrt(K_BOLTZMANN * T);
      }
    }

  } else {
    strcpy(openmc_err_msg, "Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  return 0;
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Cell* cell_pointer(int32_t cell_ind) {return cells[cell_ind];}

  int32_t cell_id(Cell* c) {return c->id_;}

  void cell_set_id(Cell* c, int32_t id) {c->id_ = id;}

  int cell_type(Cell* c) {return c->type_;}

  int cell_material(Cell* c, int32_t mat_ind) {return c->material[mat_ind-1];}
  
#ifdef CAD

  int32_t next_cell(CADCell* cur_cell, CADSurface *surf_xed ) {
    moab::EntityHandle surf = surf_xed->dagmc_ptr->entity_by_id(2,surf_xed->id);
    moab::EntityHandle vol = cur_cell->dagmc_ptr->entity_by_id(3,cur_cell->id);

    moab::EntityHandle new_vol;
    cur_cell->dagmc_ptr->next_vol(surf, vol, new_vol);

    return cur_cell->dagmc_ptr->index_by_handle(new_vol);
  }

#endif

  int32_t cell_universe(Cell* c) {return c->universe_;}

  int32_t cell_fill(Cell* c) {return c->fill_;}

  int32_t cell_n_instances(Cell* c) {return c->n_instances_;}

  int cell_distribcell_index(Cell* c) {return c->distribcell_index_;}

  int cell_material_size(Cell* c) {return c->material_.size();}

  //TODO: off-by-one
  int32_t cell_material(Cell* c, int i)
  {
    int32_t mat = c->material_[i-1];
    if (mat == MATERIAL_VOID) return MATERIAL_VOID;
    return mat + 1;
  }

  int cell_sqrtkT_size(Cell* c) {return c->sqrtkT_.size();}

  double cell_sqrtkT(Cell* c, int i) {return c->sqrtkT_[i];}

  int32_t cell_offset(Cell* c, int map) {return c->offset_[map];}

  void extend_cells_c(int32_t n)
  {
    cells.reserve(cells.size() + n);
    for (int32_t i = 0; i < n; i++) {
      global_cells.push_back(new CSGCell());
    }
    n_cells = cells.size();
  }

  int32_t universe_id(int i_univ) {return universes[i_univ]->id_;}

  void universes_to_hdf5(hid_t universes_group)
  {for (Universe* u : universes) u->to_hdf5(universes_group);}
}


} // namespace openmc
