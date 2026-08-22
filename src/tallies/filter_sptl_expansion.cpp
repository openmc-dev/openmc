#include "openmc/tallies/filter_sptl_expansion.h"

#include "openmc/xml_interface.h"

namespace openmc {

void SpatialExpansionFilter::from_xml(pugi::xml_node node)
{
  this->set_order(std::stoi(get_node_value(node, "order")));

  auto axis = get_node_value(node, "axis");
  switch (axis[0]) {
  case 'x':
    this->set_axis(Axis::x);
    break;
  case 'y':
    this->set_axis(Axis::y);
    break;
  case 'z':
    this->set_axis(Axis::z);
    break;
  default:
    throw std::runtime_error {
      "Axis for spatial expansion filters must be 'x', 'y', or 'z'"};
  }

  double min = std::stod(get_node_value(node, "min"));
  double max = std::stod(get_node_value(node, "max"));
  this->set_minmax(min, max);
}

void SpatialExpansionFilter::set_axis(Axis axis)
{
  axis_ = axis;
}

void SpatialExpansionFilter::set_minmax(double min, double max)
{
  if (max <= min) {
    throw std::invalid_argument {
      "Maximum value must be greater than minimum value"};
  }
  min_ = min;
  max_ = max;
}

double SpatialExpansionFilter::position(const Particle& p) const
{
  if (axis_ == Axis::x) {
    return p.r().x;
  } else if (axis_ == Axis::y) {
    return p.r().y;
  } else {
    return p.r().z;
  }
}

const char* SpatialExpansionFilter::axis_label() const
{
  if (axis_ == Axis::x) {
    return "x";
  } else if (axis_ == Axis::y) {
    return "y";
  } else {
    return "z";
  }
}

void SpatialExpansionFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "order", order_);
  write_dataset(filter_group, "axis", axis_label());
  write_dataset(filter_group, "min", min_);
  write_dataset(filter_group, "max", max_);
}

} // namespace openmc
