#ifndef XML_INTERFACE_H
#define XML_INTERFACE_H

#include <algorithm>  // for std::transform
#include <string>

#include "pugixml/pugixml.hpp"


bool
check_for_node(const pugi::xml_node &node, const char *name)
{
  if (node.attribute(name)) {
    return true;
  } else if (node.child(name)) {
    return true;
  } else {
    return false;
  }
}


std::string
get_node_value(const pugi::xml_node &node, const char *name)
{
  // Search for either an attribute or child tag and get the data as a char*.
  const pugi::char_t *value_char;
  if (node.attribute(name)) {
    value_char = node.attribute(name).value();
  } else if (node.child(name)) {
    value_char = node.child_value(name);
  } else {
    std::string err_msg("Node \"");
    err_msg += name;
    err_msg += "\" is not a memeber of the \"";
    err_msg += node.name();
    err_msg += "\" XML node";
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

#endif // XML_INTERFACE_H
