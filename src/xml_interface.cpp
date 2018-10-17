#include "openmc/xml_interface.h"

#include <algorithm>  // for transform
#include <sstream>

#include "openmc/error.h"


namespace openmc {

std::string
get_node_value(pugi::xml_node node, const char* name, bool lowercase,
               bool strip)
{
  // Search for either an attribute or child tag and get the data as a char*.
  const pugi::char_t* value_char;
  if (node.attribute(name)) {
    value_char = node.attribute(name).value();
  } else if (node.child(name)) {
    value_char = node.child_value(name);
  } else {
    std::stringstream err_msg;
    err_msg << "Node \"" << name << "\" is not a member of the \""
            << node.name() << "\" XML node";
    fatal_error(err_msg);
  }
  std::string value {value_char};

  // Convert to lower-case if needed
  if (lowercase) {
    std::transform(value.begin(), value.end(), value.begin(), ::tolower);
  }

  // Strip leading/trailing whitespace if needed
  if (strip) {
    value.erase(0, value.find_first_not_of(" \t\r\n"));
    value.erase(value.find_last_not_of(" \t\r\n") + 1);
  }

  return value;
}

bool
get_node_value_bool(pugi::xml_node node, const char* name)
{
  if (node.attribute(name)) {
    return node.attribute(name).as_bool();
  } else if (node.child(name)) {
    return node.child(name).text().as_bool();
  } else {
    std::stringstream err_msg;
    err_msg << "Node \"" << name << "\" is not a member of the \""
            << node.name() << "\" XML node";
    fatal_error(err_msg);
  }
  return false;
}

} // namespace openmc
