#ifndef XML_INTERFACE_H
#define XML_INTERFACE_H

#include <algorithm>  // for std::transform
#include <sstream>
#include <string>
#include <vector>

#include "pugixml/pugixml.hpp"


namespace openmc {


inline std::vector<std::string>
split(const std::string in)
{
  std::vector<std::string> out;

  for (int i = 0; i < in.size(); ) {
    // Increment i until we find a non-whitespace character.
    if (std::isspace(in[i])) {
      i++;

    } else {
      // Find the next whitespace character at j.
      int j = i + 1;
      while (j < in.size() && std::isspace(in[j]) == 0) {j++;}

      // Push-back everything between i and j.
      out.push_back(in.substr(i, j-i));
      i = j + 1; // j is whitespace so leapfrog to j+1
    }
  }

  return out;
}


inline bool
check_for_node(pugi::xml_node node, const char *name)
{
  return node.attribute(name) || node.child(name);
}


inline std::string
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
