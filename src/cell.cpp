#include "cell.h"

#include <cmath>
#include <limits>
#include <sstream>
#include <string>

#include "constants.h"
#include "error.h"
#include "hdf5_interface.h"
#include "lattice.h"
#include "surface.h"
#include "xml_interface.h"

//TODO: remove this include
#include <iostream>


namespace openmc {

//==============================================================================
// Constants
//==============================================================================

constexpr int32_t OP_LEFT_PAREN   {std::numeric_limits<int32_t>::max()};
constexpr int32_t OP_RIGHT_PAREN  {std::numeric_limits<int32_t>::max() - 1};
constexpr int32_t OP_COMPLEMENT   {std::numeric_limits<int32_t>::max() - 2};
constexpr int32_t OP_INTERSECTION {std::numeric_limits<int32_t>::max() - 3};
constexpr int32_t OP_UNION        {std::numeric_limits<int32_t>::max() - 4};

extern "C" double FP_PRECISION;

//==============================================================================
// Global variables
//==============================================================================

// Braces force n_cells to be defined here, not just declared.
extern "C" {int32_t n_cells {0};}

std::vector<Cell*> cells_c;
std::map<int32_t, int32_t> cell_dict;

std::vector<Universe*> universes_c;
std::map<int32_t, int32_t> universe_dict;

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
    id = stoi(get_node_value(cell_node, "id"));
  } else {
    fatal_error("Must specify id of cell in geometry XML file.");
  }

  //TODO: don't automatically lowercase cell and surface names
  if (check_for_node(cell_node, "name")) {
    name = get_node_value(cell_node, "name");
  }

  if (check_for_node(cell_node, "universe")) {
    universe = stoi(get_node_value(cell_node, "universe"));
  } else {
    universe = 0;
  }

  if (check_for_node(cell_node, "fill")) {
    fill = stoi(get_node_value(cell_node, "fill"));
  } else {
    fill = C_NONE;
  }

  std::string region_spec {""};
  if (check_for_node(cell_node, "region")) {
    region_spec = get_node_value(cell_node, "region");
  }

  // Get a tokenized representation of the region specification.
  region = tokenize(region_spec);
  region.shrink_to_fit();

  // Convert user IDs to surface indices.
  // Note that the index has 1 added to it in order to preserve the sign.
  for (auto it = region.begin(); it != region.end(); it++) {
    if (*it < OP_UNION) {
      *it = copysign(surface_dict[abs(*it)]+1, *it);
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

bool
Cell::contains(const double xyz[3], const double uvw[3],
               int32_t on_surface) const
{
  if (simple) {
    return contains_simple(xyz, uvw, on_surface);
  } else {
    return contains_complex(xyz, uvw, on_surface);
  }
}

std::pair<double, int32_t>
Cell::distance(const double xyz[3], const double uvw[3],
               int32_t on_surface) const
{
  double min_dist {INFTY};
  int32_t i_surf {std::numeric_limits<int32_t>::max()};

  for (int32_t token : rpn) {
    // Ignore this token if it corresponds to an operator rather than a region.
    if (token >= OP_UNION) {continue;}

    // Calculate the distance to this surface.
    // Note the off-by-one indexing
    bool coincident {token == on_surface};
    double d {surfaces_c[abs(token)-1]->distance(xyz, uvw, coincident)};
    //std::cout << token << " " << on_surface << " " << coincident << std::endl;

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

void
Cell::to_hdf5(hid_t cell_group) const
{
//  std::string group_name {"surface "};
//  group_name += std::to_string(id);
//
//  hid_t surf_group = create_group(group_id, group_name);

  if (!name.empty()) {
    write_string(cell_group, "name", name);
  }

  //TODO: Fix the off-by-one indexing.
  write_int(cell_group, "universe", universes_c[universe-1]->id);

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
        region_spec << " " << copysign(surfaces_c[abs(token)-1]->id, token);
      }
    }
    write_string(cell_group, "region", region_spec.str());
  }

//  close_group(cell_group);
}

bool
Cell::contains_simple(const double xyz[3], const double uvw[3],
                      int32_t on_surface) const
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
        bool sense = surfaces_c[abs(token)-1]->sense(xyz, uvw);
        if (sense != (token > 0)) {return false;}
      }
    }
  }
  return true;
}

bool
Cell::contains_complex(const double xyz[3], const double uvw[3],
                       int32_t on_surface) const
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
        bool sense = surfaces_c[abs(token)-1]->sense(xyz, uvw);;
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

extern "C" void
read_cells(pugi::xml_node *node)
{
  // Count the number of cells.
  for (pugi::xml_node cell_node: node->children("cell")) {n_cells++;}
  if (n_cells == 0) {
    fatal_error("No cells found in geometry.xml!");
  }

  // Allocate the vector of Cells.
  cells_c.reserve(n_cells);

  // Loop over XML cell elements and populate the array.
  for (pugi::xml_node cell_node: node->children("cell")) {
    cells_c.push_back(new Cell(cell_node));
  }

  // Populate the Universe vector and dictionary.
  for (Cell *c : cells_c) {
    int32_t uid = c->universe;
    auto it = universe_dict.find(uid);
    if (it == universe_dict.end()) {
      universes_c.push_back(new Universe());
      universes_c.back()->id = uid;
      universes_c.back()->cells.push_back(c);
      universe_dict[uid] = universes_c.size() - 1;
    } else {
      universes_c[it->second]->cells.push_back(c);
    }
  }
}

extern "C" void
adjust_indices_c()
{
  // Change cell.universe values from IDs to indices.
  for (Cell *c : cells_c) {
    auto it = universe_dict.find(c->universe);
    if (it != universe_dict.end()) {
      //TODO: Remove this off-by-one indexing.
      c->universe = it->second + 1;
    } else {
      std::stringstream err_msg;
      err_msg << "Could not find universe " << c->universe
              << " specified on cell " << c->id;
      fatal_error(err_msg);
    }
  }

  // Change all lattice universe values from IDs to indices.
  for (Lattice *l : lattices_c) {
    l->adjust_indices();
  }
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Cell* cell_pointer(int32_t cell_ind) {return cells_c[cell_ind];}

  int32_t cell_id(Cell *c) {return c->id;}

  void cell_set_id(Cell *c, int32_t id) {c->id = id;}

  int32_t cell_universe(Cell *c) {return c->universe;}

  void cell_set_universe(Cell *c, int32_t universe) {c->universe = universe;}

  bool cell_simple(Cell *c) {return c->simple;}

  bool cell_contains(Cell *c, double xyz[3], double uvw[3], int32_t on_surface)
  {return c->contains(xyz, uvw, on_surface);}

  void cell_distance(Cell *c, double xyz[3], double uvw[3], int32_t on_surface,
                     double &min_dist, int32_t &i_surf)
  {
    std::pair<double, int32_t> out = c->distance(xyz, uvw, on_surface);
    min_dist = out.first;
    i_surf = out.second;
  }

  void cell_to_hdf5(Cell *c, hid_t group) {c->to_hdf5(group);}

  void extend_cells_c(int32_t n)
  {
    cells_c.reserve(cells_c.size() + n);
    for (int32_t i = 0; i < n; i++) {
      cells_c.push_back(new Cell());
    }
    n_cells = cells_c.size();
  }
}

//extern "C" void free_memory_cells_c()
//{
//  delete cells_c;
//  cells_c = nullptr;
//  n_cells = 0;
//  cell_dict.clear();
//}

} // namespace openmc
