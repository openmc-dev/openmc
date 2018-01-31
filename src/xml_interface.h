#ifndef XML_INTERFACE_H
#define XML_INTERFACE_H

#include <algorithm>  // for std::transform
#include <sstream>
#include <string>

#include "pugixml/pugixml.hpp"


namespace openmc {

bool
check_for_node(pugi::xml_node node, const char *name)
{
  return node.attribute(name) || node.child(name);
}


std::string
get_node_value(pugi::xml_node node, const char *name)
{
  // Search for either an attribute or child tag and get the data as a char*.
  const pugi::char_t *value_char;
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

  // Convert to lowercase string.
  std::string value(value_char);
  std::transform(value.begin(), value.end(), value.begin(), ::tolower);

  // Remove whitespace.
  value.erase(0, value.find_first_not_of(" \t\r\n"));
  value.erase(value.find_last_not_of(" \t\r\n") + 1);

  return value;
}

} // namespace openmc
#endif // XML_INTERFACE_H
