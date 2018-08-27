#include "openmc/cell.h"

#include <cmath>
#include <sstream>
#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/surface.h"
#include "openmc/xml_interface.h"


namespace openmc {

//==============================================================================
// Constants
//==============================================================================

// TODO: Convert to enum
constexpr int32_t OP_LEFT_PAREN   {std::numeric_limits<int32_t>::max()};
constexpr int32_t OP_RIGHT_PAREN  {std::numeric_limits<int32_t>::max() - 1};
constexpr int32_t OP_COMPLEMENT   {std::numeric_limits<int32_t>::max() - 2};
constexpr int32_t OP_INTERSECTION {std::numeric_limits<int32_t>::max() - 3};
constexpr int32_t OP_UNION        {std::numeric_limits<int32_t>::max() - 4};

//==============================================================================
// Global variables
//==============================================================================

int32_t n_cells {0};

std::vector<Cell*> global_cells;
std::unordered_map<int32_t, int32_t> cell_map;

std::vector<Universe*> global_universes;
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
// Cell implementation
//==============================================================================

Cell::Cell(pugi::xml_node cell_node)
{
  if (check_for_node(cell_node, "id")) {
    id = std::stoi(get_node_value(cell_node, "id"));
  } else {
    fatal_error("Must specify id of cell in geometry XML file.");
  }

  if (check_for_node(cell_node, "name")) {
    name = get_node_value(cell_node, "name");
  }

  if (check_for_node(cell_node, "universe")) {
    universe = std::stoi(get_node_value(cell_node, "universe"));
  } else {
    universe = 0;
  }

  // Make sure that either material or fill was specified, but not both.
  bool fill_present = check_for_node(cell_node, "fill");
  bool material_present = check_for_node(cell_node, "material");
  if (!(fill_present || material_present)) {
    std::stringstream err_msg;
    err_msg << "Neither material nor fill was specified for cell " << id;
    fatal_error(err_msg);
  }
  if (fill_present && material_present) {
    std::stringstream err_msg;
    err_msg << "Cell " << id << " has both a material and a fill specified; "
            << "only one can be specified per cell";
    fatal_error(err_msg);
  }

  if (fill_present) {
    fill = std::stoi(get_node_value(cell_node, "fill"));
  } else {
    fill = C_NONE;
  }

  // Read the material element.  There can be zero materials (filled with a
  // universe), more than one material (distribmats), and some materials may
  // be "void".
  if (material_present) {
    std::vector<std::string> mats
         {get_node_array<std::string>(cell_node, "material", true)};
    if (mats.size() > 0) {
      material.reserve(mats.size());
      for (std::string mat : mats) {
        if (mat.compare("void") == 0) {
          material.push_back(MATERIAL_VOID);
        } else {
          material.push_back(std::stoi(mat));
        }
      }
    } else {
      std::stringstream err_msg;
      err_msg << "An empty material element was specified for cell " << id;
      fatal_error(err_msg);
    }
  }

  // Read the region specification.
  std::string region_spec;
  if (check_for_node(cell_node, "region")) {
    region_spec = get_node_value(cell_node, "region");
  }

  // Get a tokenized representation of the region specification.
  region = tokenize(region_spec);
  region.shrink_to_fit();

  // Convert user IDs to surface indices.
  for (auto &r : region) {
    if (r < OP_UNION) {
      r = copysign(surface_map[abs(r)] + 1, r);
    }
  }

  // Convert the infix region spec to RPN.
  rpn = generate_rpn(id, region);
  rpn.shrink_to_fit();

  // Check if this is a simple cell.
  simple = true;
  for (int32_t token : rpn) {
    if ((token == OP_COMPLEMENT) || (token == OP_UNION)) {
      simple = false;
      break;
    }
  }
}

//==============================================================================

bool
Cell::contains(Position r, Direction u, int32_t on_surface) const
{
  if (simple) {
    return contains_simple(r, u, on_surface);
  } else {
    return contains_complex(r, u, on_surface);
  }
}

//==============================================================================

std::pair<double, int32_t>
Cell::distance(Position r, Direction u, int32_t on_surface) const
{
  double min_dist {INFTY};
  int32_t i_surf {std::numeric_limits<int32_t>::max()};

  for (int32_t token : rpn) {
    // Ignore this token if it corresponds to an operator rather than a region.
    if (token >= OP_UNION) continue;

    // Calculate the distance to this surface.
    // Note the off-by-one indexing
    bool coincident {token == on_surface};
    double d {global_surfaces[abs(token)-1]->distance(r, u, coincident)};

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
Cell::to_hdf5(hid_t cell_group) const
{
  if (!name.empty()) {
    write_string(cell_group, "name", name, false);
  }

  //TODO: Fix the off-by-one indexing.
  write_dataset(cell_group, "universe", global_universes[universe-1]->id);

  // Write the region specification.
  if (!region.empty()) {
    std::stringstream region_spec {};
    for (int32_t token : region) {
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
             << copysign(global_surfaces[abs(token)-1]->id, token);
      }
    }
    write_string(cell_group, "region", region_spec.str(), false);
  }
}

//==============================================================================

bool
Cell::contains_simple(Position r, Direction u, int32_t on_surface) const
{
  for (int32_t token : rpn) {
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
        bool sense = global_surfaces[abs(token)-1]->sense(r, u);
        if (sense != (token > 0)) {return false;}
      }
    }
  }
  return true;
}

//==============================================================================

bool
Cell::contains_complex(Position r, Direction u, int32_t on_surface) const
{
  // Make a stack of booleans.  We don't know how big it needs to be, but we do
  // know that rpn.size() is an upper-bound.
  bool stack[rpn.size()];
  int i_stack = -1;

  for (int32_t token : rpn) {
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
        bool sense = global_surfaces[abs(token)-1]->sense(r, u);
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
// Non-method functions
//==============================================================================

extern "C" void
read_cells(pugi::xml_node* node)
{
  // Count the number of cells.
  for (pugi::xml_node cell_node: node->children("cell")) {n_cells++;}
  if (n_cells == 0) {
    fatal_error("No cells found in geometry.xml!");
  }

  // Loop over XML cell elements and populate the array.
  global_cells.reserve(n_cells);
  for (pugi::xml_node cell_node: node->children("cell")) {
    global_cells.push_back(new Cell(cell_node));
  }

  // Populate the Universe vector and map.
  for (int i = 0; i < global_cells.size(); i++) {
    int32_t uid = global_cells[i]->universe;
    auto it = universe_map.find(uid);
    if (it == universe_map.end()) {
      global_universes.push_back(new Universe());
      global_universes.back()->id = uid;
      global_universes.back()->cells.push_back(i);
      universe_map[uid] = global_universes.size() - 1;
    } else {
      global_universes[it->second]->cells.push_back(i);
    }
  }
  global_universes.shrink_to_fit();
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_cell_get_fill(int32_t index, int* type, int32_t** indices, int32_t* n)
{
  if (index >= 1 && index <= global_cells.size()) {
    //TODO: off-by-one
    Cell& c {*global_cells[index - 1]};
    *type = c.type;
    if (c.type == FILL_MATERIAL) {
      *indices = c.material.data();
      *n = c.material.size();
    } else {
      *indices = &c.fill;
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
  if (index >= 1 && index <= global_cells.size()) {
    //TODO: off-by-one
    Cell& c {*global_cells[index - 1]};
    if (type == FILL_MATERIAL) {
      c.type = FILL_MATERIAL;
      c.material.clear();
      for (int i = 0; i < n; i++) {
        int i_mat = indices[i];
        if (i_mat == MATERIAL_VOID) {
          c.material.push_back(MATERIAL_VOID);
        } else if (i_mat >= 1 && i_mat <= global_materials.size()) {
          //TODO: off-by-one
          c.material.push_back(i_mat - 1);
        } else {
          set_errmsg("Index in materials array is out of bounds.");
          return OPENMC_E_OUT_OF_BOUNDS;
        }
      }
      c.material.shrink_to_fit();
    } else if (type == FILL_UNIVERSE) {
      c.type = FILL_UNIVERSE;
    } else {
      c.type = FILL_LATTICE;
    }
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Cell* cell_pointer(int32_t cell_ind) {return global_cells[cell_ind];}

  int32_t cell_id(Cell* c) {return c->id;}

  void cell_set_id(Cell* c, int32_t id) {c->id = id;}

  int cell_type(Cell* c) {return c->type;}

  int32_t cell_universe(Cell* c) {return c->universe;}

  int32_t cell_fill(Cell* c) {return c->fill;}

  int32_t cell_n_instances(Cell* c) {return c->n_instances;}

  int cell_material_size(Cell* c) {return c->material.size();}

  //TODO: off-by-one
  int32_t cell_material(Cell* c, int i)
  {
    int32_t mat = c->material[i-1];
    if (mat == MATERIAL_VOID) return MATERIAL_VOID;
    return mat + 1;
  }

  bool cell_simple(Cell* c) {return c->simple;}

  bool cell_contains(Cell* c, double xyz[3], double uvw[3], int32_t on_surface)
  {
    Position r {xyz};
    Direction u {uvw};
    return c->contains(r, u, on_surface);
  }

  void cell_distance(Cell* c, double xyz[3], double uvw[3], int32_t on_surface,
                     double* min_dist, int32_t* i_surf)
  {
    Position r {xyz};
    Direction u {uvw};
    std::pair<double, int32_t> out = c->distance(r, u, on_surface);
    *min_dist = out.first;
    *i_surf = out.second;
  }

  int32_t cell_offset(Cell* c, int map) {return c->offset[map];}

  void cell_to_hdf5(Cell* c, hid_t group) {c->to_hdf5(group);}

  void extend_cells_c(int32_t n)
  {
    global_cells.reserve(global_cells.size() + n);
    for (int32_t i = 0; i < n; i++) {
      global_cells.push_back(new Cell());
    }
    n_cells = global_cells.size();
  }
}


} // namespace openmc
