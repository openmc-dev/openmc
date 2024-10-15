#include "openmc/xml_interface.h"

#include <fmt/core.h>

#include "openmc/error.h"
#include "openmc/string_utils.h"
#include "openmc/vector.h"

namespace openmc {

std::string get_node_value(
  pugi::xml_node node, const char* name, bool lowercase, bool strip)
{
  // Search for either an attribute or child tag and get the data as a char*.
  const pugi::char_t* value_char;
  if (node.attribute(name)) {
    value_char = node.attribute(name).value();
  } else if (node.child(name)) {
    value_char = node.child_value(name);
  } else {
    fatal_error(fmt::format(
      "Node \"{}\" is not a member of the \"{}\" XML node", name, node.name()));
  }
  std::string value {value_char};

  // Convert to lower-case if needed
  if (lowercase)
    to_lower(value);

  // Strip leading/trailing whitespace if needed
  if (strip) {
    value.erase(0, value.find_first_not_of(" \t\r\n"));
    value.erase(value.find_last_not_of(" \t\r\n") + 1);
  }

  return value;
}

bool get_node_value_bool(pugi::xml_node node, const char* name)
{
  if (node.attribute(name)) {
    return node.attribute(name).as_bool();
  } else if (node.child(name)) {
    return node.child(name).text().as_bool();
  } else {
    fatal_error(fmt::format(
      "Node \"{}\" is not a member of the \"{}\" XML node", name, node.name()));
  }
  return false;
}

vector<Position>
get_node_position_array(pugi::xml_node node, const char* name, bool lowercase)
{
  vector<double> coords = get_node_array<double>(node, name, lowercase);
  if (coords.size() % 3 != 0) {
    fatal_error(fmt::format(
      "Incorect number of coordinates in Position array ({}) for \"{}\"", coords.size(), name));
  }
  vector<Position> positions;
  positions.resize(coords.size() / 3);
  auto it = coords.begin();
  for (int i = 0; i < positions.size(); i++) {
    positions[i] = {*it++, *it++, *it++};
  }
  return positions;
}

Position get_node_position(
  pugi::xml_node node, const char* name, bool lowercase)
{
  vector<double> arr = get_node_array<double>(node, name, lowercase);
  return Position(arr);
}

} // namespace openmc
