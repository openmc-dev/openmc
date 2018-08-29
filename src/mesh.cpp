#include "openmc/mesh.h"

#include <cmath>  // for ceil
#include <string>

#include "xtensor/xeval.hpp"
#include "xtensor/xmath.hpp"

#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<std::unique_ptr<RegularMesh>> meshes;

std::unordered_map<int32_t, int32_t> mesh_map;

//==============================================================================
// RegularMesh implementation
//==============================================================================

RegularMesh::RegularMesh(pugi::xml_node node)
{
  // Copy mesh id
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));

    // Check to make sure 'id' hasn't been used
    if (mesh_map.find(id_) != mesh_map.end()) {
      fatal_error("Two or more meshes use the same unique ID: " +
        std::to_string(id_));
    }
  }

  // Read mesh type
  if (check_for_node(node, "type")) {
    auto temp = get_node_value(node, "type", true, true);
    if (temp == "regular") {
      // TODO: move elsewhere
    } else {
      fatal_error("Invalid mesh type: " + temp);
    }
  }

  // Determine number of dimensions for mesh
  if (check_for_node(node, "dimension")) {
    shape_ = get_node_xarray<int>(node, "dimension");
    int n = n_dimension_ = shape_.size();
    if (n != 1 && n != 2 && n != 3) {
      fatal_error("Mesh must be one, two, or three dimensions.");
    }

    // Check that dimensions are all greater than zero
    if (xt::any(shape_ <= 0)) {
      fatal_error("All entries on the <dimension> element for a tally "
        "mesh must be positive.");
    }
  }

  // Check for lower-left coordinates
  if (check_for_node(node, "lower_left")) {
    // Read mesh lower-left corner location
    lower_left_ = get_node_xarray<double>(node, "lower_left");
  } else {
    fatal_error("Must specify <lower_left> on a mesh.");
  }

  if (check_for_node(node, "width")) {
    // Make sure both upper-right or width were specified
    if (check_for_node(node, "upper_right")) {
      fatal_error("Cannot specify both <upper_right> and <width> on a mesh.");
    }

    width_ = get_node_xarray<double>(node, "width");

    // Check to ensure width has same dimensions
    auto n = width_.size();
    if (n != lower_left_.size()) {
      fatal_error("Number of entries on <width> must be the same as "
        "the number of entries on <lower_left>.");
    }

    // Check for negative widths
    if (xt::any(width_ < 0.0)) {
      fatal_error("Cannot have a negative <width> on a tally mesh.");
    }

    // Set width and upper right coordinate
    upper_right_ = xt::eval(lower_left_ + shape_ * width_);

  } else if (check_for_node(node, "upper_right")) {
    upper_right_ = get_node_xarray<double>(node, "upper_right");

    // Check to ensure width has same dimensions
    auto n = upper_right_.size();
    if (n != lower_left_.size()) {
      fatal_error("Number of entries on <upper_right> must be the "
        "same as the number of entries on <lower_left>.");
    }

    // Check that upper-right is above lower-left
    if (xt::any(upper_right_ < lower_left_)) {
      fatal_error("The <upper_right> coordinates must be greater than "
        "the <lower_left> coordinates on a tally mesh.");
    }

    // Set width and upper right coordinate
    width_ = xt::eval((upper_right_ - lower_left_) / shape_);
  } else {
    fatal_error("Must specify either <upper_right> and <width> on a mesh.");
  }

  if (shape_.dimension() > 0) {
    if (shape_.size() != lower_left_.size()) {
      fatal_error("Number of entries on <lower_left> must be the same "
        "as the number of entries on <dimension>.");
    }

    // Set volume fraction
    volume_frac_ = 1.0/xt::prod(shape_)();
  }
}

int RegularMesh::get_bin(Position r)
{
  // Loop over the dimensions of the mesh
  for (int i = 0; i < n_dimension_; ++i) {
    // Check for cases where particle is outside of mesh
    if (r[i] < lower_left_[i]) {
      return -1;
    } else if (r[i] > upper_right_[i]) {
      return -1;
    }
  }

  // Determine indices
  int ijk[n_dimension_];
  bool in_mesh;
  get_indices(r, ijk, &in_mesh);
  if (!in_mesh) return -1;

  // Convert indices to bin
  return get_bin_from_indices(ijk);
}

int RegularMesh::get_bin_from_indices(const int* ijk)
{
  if (n_dimension_ == 1) {
    return ijk[0];
  } else if (n_dimension_ == 2) {
    return (ijk[1] - 1)*shape_[0] + ijk[0];
  } else if (n_dimension_ == 3) {
    return ((ijk[2] - 1)*shape_[1] + (ijk[1] - 1))*shape_[0] + ijk[0];
  }
}

void RegularMesh::get_indices(Position r, int* ijk, bool* in_mesh)
{
  // Find particle in mesh
  *in_mesh = true;
  for (int i = 0; i < n_dimension_; ++i) {
    ijk[i] = std::ceil((r[i] - lower_left_[i]) / width_[i]);

    // Check if indices are within bounds
    if (ijk[i] < 1 || ijk[i] > shape_[i]) *in_mesh = false;
  }
}

void RegularMesh::get_indices_from_bin(int bin, int* ijk)
{
  if (n_dimension_ == 1) {
    ijk[0] = bin;
  } else if (n_dimension_ == 2) {
    ijk[0] = (bin - 1) % shape_[0] + 1;
    ijk[1] = (bin - 1) / shape_[0] + 1;
  } else if (n_dimension_ == 3) {
    ijk[0] = (bin - 1) % shape_[0] + 1;
    ijk[1] = (bin - 1) % (shape_[0] * shape_[1]) / shape_[0] + 1;
    ijk[2] = (bin - 1) / (shape_[0] * shape_[1]) + 1;
  }
}

bool RegularMesh::intersects(Position r0, Position r1)
{
  switch(n_dimension_) {
    case 1:
      return intersects_1d(r0, r1);
    case 2:
      return intersects_2d(r0, r1);
    case 3:
      return intersects_3d(r0, r1);
  }
}

bool RegularMesh::intersects_1d(Position r0, Position r1)
{
  // Copy coordinates of mesh lower_left and upper_right
  double left = lower_left_[0];
  double right = upper_right_[0];

  // Check if line intersects either left or right surface
  if (r0.x < left) {
    return r1.x > left;
  } else if (r0.x < right) {
    return r1.x < left || r1.x > right;
  } else {
    return r1.x < right;
  }
}

bool RegularMesh::intersects_2d(Position r0, Position r1)
{
  // Copy coordinates of starting point
  double x0 = r0.x;
  double y0 = r0.y;

  // Copy coordinates of ending point
  double x1 = r1.x;
  double y1 = r1.y;

  // Copy coordinates of mesh lower_left
  double xm0 = lower_left_[0];
  double ym0 = lower_left_[1];

  // Copy coordinates of mesh upper_right
  double xm1 = upper_right_[0];
  double ym1 = upper_right_[1];

  // Check if line intersects left surface -- calculate the intersection point y
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1) {
      return true;
    }
  }

  // Check if line intersects back surface -- calculate the intersection point
  // x
  if ((y0 < ym0 && y1 > ym0) || (y0 > ym0 && y1 < ym0)) {
    double xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1) {
      return true;
    }
  }

  // Check if line intersects right surface -- calculate the intersection
  // point y
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1) {
      return true;
    }
  }

  // Check if line intersects front surface -- calculate the intersection point
  // x
  if ((y0 < ym1 && y1 > ym1) || (y0 > ym1 && y1 < ym1)) {
    double xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1) {
      return true;
    }
  }
  return false;
}

bool RegularMesh::intersects_3d(Position r0, Position r1)
{
  // Copy coordinates of starting point
  double x0 = r0.x;
  double y0 = r0.y;
  double z0 = r0.z;

  // Copy coordinates of ending point
  double x1 = r1.x;
  double y1 = r1.y;
  double z1 = r1.z;

  // Copy coordinates of mesh lower_left
  double xm0 = lower_left_[0];
  double ym0 = lower_left_[1];
  double zm0 = lower_left_[2];

  // Copy coordinates of mesh upper_right
  double xm1 = upper_right_[0];
  double ym1 = upper_right_[1];
  double zm1 = upper_right_[2];

  // Check if line intersects left surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm0 && x1 > xm0) || (x0 > xm0 && x1 < xm0)) {
    double yi = y0 + (xm0 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm0 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1 && zi >= zm0 && zi < zm1) {
      return true;
    }
  }

  // Check if line intersects back surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym0 && y1 > ym0) || (y0 > ym0 && y1 < ym0)) {
    double xi = x0 + (ym0 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym0 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1 && zi >= zm0 && zi < zm1) {
      return true;
    }
  }

  // Check if line intersects bottom surface -- calculate the intersection
  // point (x,y)
  if ((z0 < zm0 && z1 > zm0) || (z0 > zm0 && z1 < zm0)) {
    double xi = x0 + (zm0 - z0) * (x1 - x0) / (z1 - z0);
    double yi = y0 + (zm0 - z0) * (y1 - y0) / (z1 - z0);
    if (xi >= xm0 && xi < xm1 && yi >= ym0 && yi < ym1) {
      return true;
    }
  }

  // Check if line intersects right surface -- calculate the intersection point
  // (y,z)
  if ((x0 < xm1 && x1 > xm1) || (x0 > xm1 && x1 < xm1)) {
    double yi = y0 + (xm1 - x0) * (y1 - y0) / (x1 - x0);
    double zi = z0 + (xm1 - x0) * (z1 - z0) / (x1 - x0);
    if (yi >= ym0 && yi < ym1 && zi >= zm0 && zi < zm1) {
      return true;
    }
  }

  // Check if line intersects front surface -- calculate the intersection point
  // (x,z)
  if ((y0 < ym1 && y1 > ym1) || (y0 > ym1 && y1 < ym1)) {
    double xi = x0 + (ym1 - y0) * (x1 - x0) / (y1 - y0);
    double zi = z0 + (ym1 - y0) * (z1 - z0) / (y1 - y0);
    if (xi >= xm0 && xi < xm1 && zi >= zm0 && zi < zm1) {
      return true;
    }
  }

  // Check if line intersects top surface -- calculate the intersection point
  // (x,y)
  if ((z0 < zm1 && z1 > zm1) || (z0 > zm1 && z1 < zm1)) {
    double xi = x0 + (zm1 - z0) * (x1 - x0) / (z1 - z0);
    double yi = y0 + (zm1 - z0) * (y1 - y0) / (z1 - z0);
    if (xi >= xm0 && xi < xm1 && yi >= ym0 && yi < ym1) {
      return true;
    }
  }
  return false;
}

void RegularMesh::to_hdf5(hid_t group)
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  write_dataset(mesh_group, "type", "regular");
  write_dataset(mesh_group, "dimension", shape_);
  write_dataset(mesh_group, "lower_left", lower_left_);
  write_dataset(mesh_group, "upper_right", upper_right_);
  write_dataset(mesh_group, "width", width_);

  close_group(mesh_group);
}

} // namespace openmc
