
#include "openmc/cell.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <iterator>
#include <set>
#include <sstream>
#include <string>

#include <fmt/core.h>

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

} // namespace model

//==============================================================================
// Cell implementation
//==============================================================================

int32_t Cell::n_instances() const
{
  return model::universes[universe_]->n_instances_;
}

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

double Cell::density_mult(int32_t instance) const
{
  if (instance >= 0) {
    return density_mult_.size() == 1 ? density_mult_.at(0)
                                     : density_mult_.at(instance);
  } else {
    return density_mult_[0];
  }
}

double Cell::density(int32_t instance) const
{
  const int32_t mat_index = material(instance);
  if (mat_index == MATERIAL_VOID)
    return 0.0;

  return density_mult(instance) * model::materials[mat_index]->density_gpcc();
}

void Cell::set_temperature(double T, int32_t instance, bool set_contained)
{
  if (settings::temperature_method == TemperatureMethod::INTERPOLATION) {
    if (T < (data::temperature_min - settings::temperature_tolerance)) {
      throw std::runtime_error {
        fmt::format("Temperature of {} K is below minimum temperature at "
                    "which data is available of {} K.",
          T, data::temperature_min)};
    } else if (T > (data::temperature_max + settings::temperature_tolerance)) {
      throw std::runtime_error {
        fmt::format("Temperature of {} K is above maximum temperature at "
                    "which data is available of {} K.",
          T, data::temperature_max)};
    }
  }

  if (type_ == Fill::MATERIAL) {
    if (instance >= 0) {
      // If temperature vector is not big enough, resize it first
      if (sqrtkT_.size() != n_instances())
        sqrtkT_.resize(n_instances(), sqrtkT_[0]);

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
      assert(cell->type_ == Fill::MATERIAL);
      auto& instances = entry.second;
      for (auto instance : instances) {
        cell->set_temperature(T, instance);
      }
    }
  }
}

void Cell::set_density(double density, int32_t instance, bool set_contained)
{
  if (type_ != Fill::MATERIAL && !set_contained) {
    fatal_error(
      fmt::format("Attempted to set the density multiplier of cell {} "
                  "which is not filled by a material.",
        id_));
  }

  if (type_ == Fill::MATERIAL) {
    const int32_t mat_index = material(instance);
    if (mat_index == MATERIAL_VOID)
      return;

    if (instance >= 0) {
      // If density multiplier vector is not big enough, resize it first
      if (density_mult_.size() != n_instances())
        density_mult_.resize(n_instances(), density_mult_[0]);

      // Set density multiplier for the corresponding instance
      density_mult_.at(instance) =
        density / model::materials[mat_index]->density_gpcc();
    } else {
      // Set density multiplier for all instances
      for (auto& x : density_mult_) {
        x = density / model::materials[mat_index]->density_gpcc();
      }
    }
  } else {
    auto contained_cells = this->get_contained_cells(instance);
    for (const auto& entry : contained_cells) {
      auto& cell = model::cells[entry.first];
      assert(cell->type_ == Fill::MATERIAL);
      auto& instances = entry.second;
      for (auto instance : instances) {
        cell->set_density(density, instance);
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

  // Write density for one or more cell instances
  if (type_ == Fill::MATERIAL && material_.size() > 0) {
    vector<double> density;
    for (int32_t i = 0; i < density_mult_.size(); ++i)
      density.push_back(this->density(i));

    write_dataset(cell_group, "density", density);
  }

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
  if (n_temps > 1 && n_temps != n_instances()) {
    fatal_error(fmt::format(
      "Number of temperatures for cell {} doesn't match number of instances",
      id_));
  }

  // Modify temperatures for the cell
  sqrtkT_.clear();
  sqrtkT_.resize(temps.size());
  for (int64_t i = 0; i < temps.size(); ++i) {
    this->set_temperature(temps[i], i);
  }

  // Read densities
  if (object_exists(cell_group, "density")) {
    vector<double> density;
    read_dataset(cell_group, "density", density);

    // Ensure number of densities makes sense
    auto n_density = density.size();
    if (n_density > 1 && n_density != n_instances()) {
      fatal_error(fmt::format("Number of densities for cell {} "
                              "doesn't match number of instances",
        id_));
    }

    // Set densities.
    for (int32_t i = 0; i < n_density; ++i) {
      this->set_density(density[i], i);
    }
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

    write_dataset(group, "density_mult", density_mult_);

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
      fatal_error(fmt::format("Cell {} is filled with the same universe that "
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

  // Read the density element which can be distributed similar to temperature.
  // These get assigned to the density multiplier, requiring a division by
  // the material density.
  // Note: calculating the actual density multiplier is deferred until materials
  // are finalized. density_mult_ contains the true density in the meantime.
  if (check_for_node(cell_node, "density")) {
    density_mult_ = get_node_array<double>(cell_node, "density");
    density_mult_.shrink_to_fit();
    xml_set_density_ = true;

    // Make sure this is a material-filled cell.
    if (material_.size() == 0) {
      fatal_error(fmt::format(
        "Cell {} was specified with a density but no material. Density"
        "specification is only valid for cells filled with a material.",
        id_));
    }

    // Make sure this is a non-void material.
    for (auto mat_id : material_) {
      if (mat_id == MATERIAL_VOID) {
        fatal_error(fmt::format(
          "Cell {} was specified with a density, but contains a void "
          "material. Density specification is only valid for cells "
          "filled with a non-void material.",
          id_));
      }
    }

    // Make sure all densities are non-negative and greater than zero.
    for (auto rho : density_mult_) {
      if (rho <= 0) {
        fatal_error(fmt::format(
          "Cell {} was specified with a density less than or equal to zero",
          id_));
      }
    }
  }

  // Read the region specification.
  std::string region_spec;
  if (check_for_node(cell_node, "region")) {
    region_spec = get_node_value(cell_node, "region");
  }

  // Get a tokenized representation of the region specification and apply De
  // Morgans law
  Region region(region_spec, id_);
  region_ = region;

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

void CSGCell::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "geom_type", "csg", false);
  write_string(group_id, "region", region_.str(), false);
}

//==============================================================================

vector<int32_t>::iterator CSGCell::find_left_parenthesis(
  vector<int32_t>::iterator start, const vector<int32_t>& infix)
{
  // start search at zero
  int parenthesis_level = 0;
  auto it = start;
  while (it != infix.begin()) {
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

//==============================================================================
// Region implementation
//==============================================================================

Region::Region(std::string region_spec, int32_t cell_id)
{
  // Check if region_spec is not empty.
  if (!region_spec.empty()) {
    // Parse all halfspaces and operators except for intersection (whitespace).
    for (int i = 0; i < region_spec.size();) {
      if (region_spec[i] == '(') {
        expression_.push_back(OP_LEFT_PAREN);
        i++;

      } else if (region_spec[i] == ')') {
        expression_.push_back(OP_RIGHT_PAREN);
        i++;

      } else if (region_spec[i] == '|') {
        expression_.push_back(OP_UNION);
        i++;

      } else if (region_spec[i] == '~') {
        expression_.push_back(OP_COMPLEMENT);
        i++;

      } else if (region_spec[i] == '-' || region_spec[i] == '+' ||
                 std::isdigit(region_spec[i])) {
        // This is the start of a halfspace specification.  Iterate j until we
        // find the end, then push-back everything between i and j.
        int j = i + 1;
        while (j < region_spec.size() && std::isdigit(region_spec[j])) {
          j++;
        }
        expression_.push_back(std::stoi(region_spec.substr(i, j - i)));
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
    while (i < expression_.size() - 1) {
      bool left_compat {
        (expression_[i] < OP_UNION) || (expression_[i] == OP_RIGHT_PAREN)};
      bool right_compat {(expression_[i + 1] < OP_UNION) ||
                         (expression_[i + 1] == OP_LEFT_PAREN) ||
                         (expression_[i + 1] == OP_COMPLEMENT)};
      if (left_compat && right_compat) {
        expression_.insert(expression_.begin() + i + 1, OP_INTERSECTION);
      }
      i++;
    }

    // Remove complement operators using DeMorgan's laws
    auto it = std::find(expression_.begin(), expression_.end(), OP_COMPLEMENT);
    while (it != expression_.end()) {
      // Erase complement
      expression_.erase(it);

      // Define stop given left parenthesis or not
      auto stop = it;
      if (*it == OP_LEFT_PAREN) {
        int depth = 1;
        do {
          stop++;
          if (*stop > OP_COMPLEMENT) {
            if (*stop == OP_RIGHT_PAREN) {
              depth--;
            } else {
              depth++;
            }
          }
        } while (depth > 0);
        it++;
      }

      // apply DeMorgan's law to any surfaces/operators between these
      // positions in the RPN
      apply_demorgan(it, stop);
      // update iterator position
      it = std::find(expression_.begin(), expression_.end(), OP_COMPLEMENT);
    }

    // Convert user IDs to surface indices.
    for (auto& r : expression_) {
      if (r < OP_UNION) {
        const auto& it {model::surface_map.find(abs(r))};
        if (it == model::surface_map.end()) {
          throw std::runtime_error {
            "Invalid surface ID " + std::to_string(abs(r)) +
            " specified in region for cell " + std::to_string(cell_id) + "."};
        }
        r = (r > 0) ? it->second + 1 : -(it->second + 1);
      }
    }

    // Check if this is a simple cell.
    simple_ = true;
    for (int32_t token : expression_) {
      if (token == OP_UNION) {
        simple_ = false;
        // Ensure intersections have precedence over unions
        add_precedence();
        break;
      }
    }

    // If this cell is simple, remove all the superfluous operator tokens.
    if (simple_) {
      for (auto it = expression_.begin(); it != expression_.end(); it++) {
        if (*it == OP_INTERSECTION || *it > OP_COMPLEMENT) {
          expression_.erase(it);
          it--;
        }
      }
    }
    expression_.shrink_to_fit();

  } else {
    simple_ = true;
  }
}

//==============================================================================

void Region::apply_demorgan(
  vector<int32_t>::iterator start, vector<int32_t>::iterator stop)
{
  do {
    if (*start < OP_UNION) {
      *start *= -1;
    } else if (*start == OP_UNION) {
      *start = OP_INTERSECTION;
    } else if (*start == OP_INTERSECTION) {
      *start = OP_UNION;
    }
    start++;
  } while (start < stop);
}

//==============================================================================
//! Add precedence for infix regions so intersections have higher
//! precedence than unions using parentheses.
//==============================================================================

int64_t Region::add_parentheses(int64_t start)
{
  int32_t start_token = expression_[start];
  // Add left parenthesis and set new position to be after parenthesis
  if (start_token == OP_UNION) {
    start += 2;
  }
  expression_.insert(expression_.begin() + start - 1, OP_LEFT_PAREN);

  // Keep track of return iterator distance. If we don't encounter a left
  // parenthesis, we return an iterator corresponding to wherever the right
  // parenthesis is inserted. If a left parenthesis is encountered, an iterator
  // corresponding to the left parenthesis is returned. Also note that we keep
  // track of a *distance* instead of an iterator because the underlying memory
  // allocation may change.
  std::size_t return_it_dist = 0;

  // Add right parenthesis
  // While the start iterator is within the bounds of infix
  while (start + 1 < expression_.size()) {
    start++;

    // If the current token is an operator and is different than the start token
    if (expression_[start] >= OP_UNION && expression_[start] != start_token) {
      // Skip wrapped regions but save iterator position to check precedence and
      // add right parenthesis, right parenthesis position depends on the
      // operator, when the operator is a union then do not include the operator
      // in the region, when the operator is an intersection then include the
      // operator and next surface
      if (expression_[start] == OP_LEFT_PAREN) {
        return_it_dist = start;
        int depth = 1;
        do {
          start++;
          if (expression_[start] > OP_COMPLEMENT) {
            if (expression_[start] == OP_RIGHT_PAREN) {
              depth--;
            } else {
              depth++;
            }
          }
        } while (depth > 0);
      } else {
        if (start_token == OP_UNION) {
          --start;
        }
        expression_.insert(expression_.begin() + start, OP_RIGHT_PAREN);
        if (return_it_dist > 0) {
          return return_it_dist;
        } else {
          return start - 1;
        }
      }
    }
  }
  // If we get here a right parenthesis hasn't been placed,
  // return iterator
  expression_.push_back(OP_RIGHT_PAREN);
  if (return_it_dist > 0) {
    return return_it_dist;
  } else {
    return start - 1;
  }
}

//==============================================================================

void Region::add_precedence()
{
  int32_t current_op = 0;
  std::size_t current_dist = 0;

  for (int64_t i = 0; i < expression_.size(); i++) {
    int32_t token = expression_[i];

    if (token == OP_UNION || token == OP_INTERSECTION) {
      if (current_op == 0) {
        // Set the current operator if is hasn't been set
        current_op = token;
        current_dist = i;
      } else if (token != current_op) {
        // If the current operator doesn't match the token, add parenthesis to
        // assert precedence
        if (current_op == OP_INTERSECTION) {
          i = add_parentheses(current_dist);
        } else {
          i = add_parentheses(i);
        }
        current_op = 0;
        current_dist = 0;
      }
    } else if (token > OP_COMPLEMENT) {
      // If the token is a parenthesis reset the current operator
      current_op = 0;
      current_dist = 0;
    }
  }
}

//==============================================================================
//! Convert infix region specification to Reverse Polish Notation (RPN)
//!
//! This function uses the shunting-yard algorithm.
//==============================================================================

vector<int32_t> Region::generate_postfix(int32_t cell_id) const
{
  vector<int32_t> rpn;
  vector<int32_t> stack;

  for (int32_t token : expression_) {
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

std::string Region::str() const
{
  std::stringstream region_spec {};
  if (!expression_.empty()) {
    for (int32_t token : expression_) {
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
  }
  return region_spec.str();
}

//==============================================================================

std::pair<double, int32_t> Region::distance(
  Position r, Direction u, int32_t on_surface) const
{
  double min_dist {INFTY};
  int32_t i_surf {std::numeric_limits<int32_t>::max()};

  for (int32_t token : expression_) {
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

  return {min_dist, i_surf};
}

//==============================================================================

bool Region::contains(Position r, Direction u, int32_t on_surface) const
{
  if (simple_) {
    return contains_simple(r, u, on_surface);
  } else {
    return contains_complex(r, u, on_surface);
  }
}

//==============================================================================

bool Region::contains_simple(Position r, Direction u, int32_t on_surface) const
{
  for (int32_t token : expression_) {
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

bool Region::contains_complex(Position r, Direction u, int32_t on_surface) const
{
  bool in_cell = true;
  int total_depth = 0;

  // For each token
  for (auto it = expression_.begin(); it != expression_.end(); it++) {
    int32_t token = *it;

    // If the token is a surface evaluate the sense
    // If the token is a union or intersection check to
    // short circuit
    if (token < OP_UNION) {
      if (token == on_surface) {
        in_cell = true;
      } else if (-token == on_surface) {
        in_cell = false;
      } else {
        // Note the off-by-one indexing
        bool sense = model::surfaces[abs(token) - 1]->sense(r, u);
        in_cell = (sense == (token > 0));
      }
    } else if ((token == OP_UNION && in_cell == true) ||
               (token == OP_INTERSECTION && in_cell == false)) {
      // If the total depth is zero return
      if (total_depth == 0) {
        return in_cell;
      }

      total_depth--;

      // While the iterator is within the bounds of the vector
      int depth = 1;
      do {
        // Get next token
        it++;
        int32_t next_token = *it;

        // If the token is an a parenthesis
        if (next_token > OP_COMPLEMENT) {
          // Adjust depth accordingly
          if (next_token == OP_RIGHT_PAREN) {
            depth--;
          } else {
            depth++;
          }
        }
      } while (depth > 0);
    } else if (token == OP_LEFT_PAREN) {
      total_depth++;
    } else if (token == OP_RIGHT_PAREN) {
      total_depth--;
    }
  }
  return in_cell;
}

//==============================================================================

BoundingBox Region::bounding_box(int32_t cell_id) const
{
  if (simple_) {
    return bounding_box_simple();
  } else {
    auto postfix = generate_postfix(cell_id);
    return bounding_box_complex(postfix);
  }
}

//==============================================================================

BoundingBox Region::bounding_box_simple() const
{
  BoundingBox bbox;
  for (int32_t token : expression_) {
    bbox &= model::surfaces[abs(token) - 1]->bounding_box(token > 0);
  }
  return bbox;
}

//==============================================================================

BoundingBox Region::bounding_box_complex(vector<int32_t> postfix) const
{
  vector<BoundingBox> stack(postfix.size());
  int i_stack = -1;

  for (auto& token : postfix) {
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

  assert(i_stack == 0);
  return stack.front();
}

//==============================================================================

vector<int32_t> Region::surfaces() const
{
  if (simple_) {
    return expression_;
  }

  vector<int32_t> surfaces = expression_;

  auto it = std::find_if(surfaces.begin(), surfaces.end(),
    [&](const auto& value) { return value >= OP_UNION; });

  while (it != surfaces.end()) {
    surfaces.erase(it);

    it = std::find_if(surfaces.begin(), surfaces.end(),
      [&](const auto& value) { return value >= OP_UNION; });
  }

  return surfaces;
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

  populate_universes();

  // Allocate the cell overlap count if necessary.
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }

  if (model::cells.size() == 0) {
    fatal_error("No cells were found in the geometry.xml file");
  }
}

void populate_universes()
{
  // Used to map universe index to the index of an implicit complement cell for
  // DAGMC universes
  std::unordered_map<int, int> implicit_comp_cells;

  // Populate the Universe vector and map.
  for (int index_cell = 0; index_cell < model::cells.size(); index_cell++) {
    int32_t uid = model::cells[index_cell]->universe_;
    auto it = model::universe_map.find(uid);
    if (it == model::universe_map.end()) {
      model::universes.push_back(make_unique<Universe>());
      model::universes.back()->id_ = uid;
      model::universes.back()->cells_.push_back(index_cell);
      model::universe_map[uid] = model::universes.size() - 1;
    } else {
#ifdef OPENMC_DAGMC_ENABLED
      // Skip implicit complement cells for now
      Universe* univ = model::universes[it->second].get();
      DAGUniverse* dag_univ = dynamic_cast<DAGUniverse*>(univ);
      if (dag_univ && (dag_univ->implicit_complement_idx() == index_cell)) {
        implicit_comp_cells[it->second] = index_cell;
        continue;
      }
#endif

      model::universes[it->second]->cells_.push_back(index_cell);
    }
  }

  // Add DAGUniverse implicit complement cells last
  for (const auto& it : implicit_comp_cells) {
    int index_univ = it.first;
    int index_cell = it.second;
    model::universes[index_univ]->cells_.push_back(index_cell);
  }

  model::universes.shrink_to_fit();
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

extern "C" int openmc_cell_set_density(
  int32_t index, double density, const int32_t* instance, bool set_contained)
{
  if (index < 0 || index >= model::cells.size()) {
    strcpy(openmc_err_msg, "Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  int32_t instance_index = instance ? *instance : -1;
  try {
    model::cells[index]->set_density(density, instance_index, set_contained);
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

extern "C" int openmc_cell_get_density(
  int32_t index, const int32_t* instance, double* density)
{
  if (index < 0 || index >= model::cells.size()) {
    strcpy(openmc_err_msg, "Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  int32_t instance_index = instance ? *instance : -1;
  try {
    if (model::cells[index]->type_ != Fill::MATERIAL) {
      fatal_error(
        fmt::format("Cell {}, instance {} is not filled with a material.",
          model::cells[index]->id_, instance_index));
    }

    int32_t mat_index = model::cells[index]->material(instance_index);
    if (mat_index == MATERIAL_VOID) {
      *density = 0.0;
    } else {
      *density = model::cells[index]->density_mult(instance_index) *
                 model::materials[mat_index]->density_gpcc();
    }
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

  int64_t cell_index;
  int64_t lattice_index;
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
    if (distribcell_index == C_NONE)
      return 0;

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
  int32_t instance, const Position& r) const
{

  // create a temporary particle
  GeometryState dummy_particle {};
  dummy_particle.r() = r;
  dummy_particle.u() = {0., 0., 1.};

  return find_parent_cells(instance, dummy_particle);
}

vector<ParentCell> Cell::find_parent_cells(
  int32_t instance, GeometryState& p) const
{
  // look up the particle's location
  exhaustive_find_cell(p);
  const auto& coords = p.coord();

  // build a parent cell stack from the particle coordinates
  ParentCellStack stack;
  bool cell_found = false;
  for (auto it = coords.begin(); it != coords.end(); it++) {
    const auto& coord = *it;
    const auto& cell = model::cells[coord.cell()];
    // if the cell at this level matches the current cell, stop adding to the
    // stack
    if (coord.cell() == model::cell_map[this->id_]) {
      cell_found = true;
      break;
    }

    // if filled with a lattice, get the lattice index from the next
    // level in the coordinates to push to the stack
    int lattice_idx = C_NONE;
    if (cell->type_ == Fill::LATTICE) {
      const auto& next_coord = *(it + 1);
      lattice_idx = model::lattices[next_coord.lattice()]->get_flat_index(
        next_coord.lattice_index());
    }
    stack.push(coord.universe(), {coord.cell(), lattice_idx});
  }

  // if this loop finished because the cell was found and
  // the instance matches the one requested in the call
  // we have the correct path and can return the stack
  if (cell_found &&
      stack.compute_instance(this->distribcell_index_) == instance) {
    return stack.parent_cells();
  }

  // fall back on an exhaustive search for the cell's parents
  return exhaustive_find_parent_cells(instance);
}

vector<ParentCell> Cell::exhaustive_find_parent_cells(int32_t instance) const
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
  *num_instances = model::cells[index]->n_instances();
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
