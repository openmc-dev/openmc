#ifndef XML_INTERFACE_H
#define XML_INTERFACE_H

#include <string>
#include <vector>

#include "pugixml/pugixml.hpp"


namespace openmc {

inline bool
check_for_node(pugi::xml_node node, const char *name)
{
  return node.attribute(name) || node.child(name);
}

std::string get_node_value(pugi::xml_node node, const char *name);

} // namespace openmc
#endif // XML_INTERFACE_H
