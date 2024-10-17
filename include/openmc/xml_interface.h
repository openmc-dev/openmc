#ifndef OPENMC_XML_INTERFACE_H
#define OPENMC_XML_INTERFACE_H

#include <cstddef> // for size_t
#include <sstream> // for stringstream
#include <string>

#include "pugixml.hpp"
#include "xtensor/xadapt.hpp"
#include "xtensor/xarray.hpp"

#include "openmc/position.h"
#include "openmc/vector.h"

namespace openmc {

inline bool check_for_node(pugi::xml_node node, const char* name)
{
  return node.attribute(name) || node.child(name);
}

std::string get_node_value(pugi::xml_node node, const char* name,
  bool lowercase = false, bool strip = false);

bool get_node_value_bool(pugi::xml_node node, const char* name);

template<typename T>
vector<T> get_node_array(
  pugi::xml_node node, const char* name, bool lowercase = false)
{
  // Get value of node attribute/child
  std::string s {get_node_value(node, name, lowercase)};

  // Read values one by one into vector
  std::stringstream iss {s};
  T value;
  vector<T> values;
  while (iss >> value)
    values.push_back(value);

  return values;
}

template<typename T>
xt::xarray<T> get_node_xarray(
  pugi::xml_node node, const char* name, bool lowercase = false)
{
  vector<T> v = get_node_array<T>(node, name, lowercase);
  vector<std::size_t> shape = {v.size()};
  return xt::adapt(v, shape);
}

std::vector<Position> get_node_position_array(
  pugi::xml_node node, const char* name, bool lowercase = false);

Position get_node_position(
  pugi::xml_node node, const char* name, bool lowercase = false);

} // namespace openmc
#endif // OPENMC_XML_INTERFACE_H
