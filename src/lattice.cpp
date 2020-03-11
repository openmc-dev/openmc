#include "openmc/lattice.h"

#include <cmath>
#include <string>
#include <vector>

#include <fmt/core.h>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/hdf5_interface.h"
#include "openmc/string_utils.h"
#include "openmc/xml_interface.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  std::vector<std::unique_ptr<Lattice>> lattices;
  std::unordered_map<int32_t, int32_t> lattice_map;
}

//==============================================================================
// Lattice implementation
//==============================================================================

Lattice::Lattice(pugi::xml_node lat_node)
{
  if (check_for_node(lat_node, "id")) {
    id_ = std::stoi(get_node_value(lat_node, "id"));
  } else {
    fatal_error("Must specify id of lattice in geometry XML file.");
  }

  if (check_for_node(lat_node, "name")) {
    name_ = get_node_value(lat_node, "name");
  }

  if (check_for_node(lat_node, "outer")) {
    outer_ = std::stoi(get_node_value(lat_node, "outer"));
  }
}

//==============================================================================

LatticeIter Lattice::begin()
{return LatticeIter(*this, 0);}

LatticeIter Lattice::end()
{return LatticeIter(*this, universes_.size());}

ReverseLatticeIter Lattice::rbegin()
{return ReverseLatticeIter(*this, universes_.size()-1);}

ReverseLatticeIter Lattice::rend()
{return ReverseLatticeIter(*this, -1);}

//==============================================================================

void
Lattice::adjust_indices()
{
  // Adjust the indices for the universes array.
  for (LatticeIter it = begin(); it != end(); ++it) {
    int uid = *it;
    auto search = model::universe_map.find(uid);
    if (search != model::universe_map.end()) {
      *it = search->second;
    } else {
      fatal_error(fmt::format(
        "Invalid universe number {} specified on lattice {}", uid, id_));
    }
  }

  // Adjust the index for the outer universe.
  if (outer_ != NO_OUTER_UNIVERSE) {
    auto search = model::universe_map.find(outer_);
    if (search != model::universe_map.end()) {
      outer_ = search->second;
    } else {
      fatal_error(fmt::format(
        "Invalid universe number {} specified on lattice {}", outer_, id_));
    }
  }
}

//==============================================================================

int32_t
Lattice::fill_offset_table(int32_t offset, int32_t target_univ_id, int map,
  std::unordered_map<int32_t, int32_t>& univ_count_memo)
{
  for (LatticeIter it = begin(); it != end(); ++it) {
    offsets_[map * universes_.size() + it.indx_] = offset;
    offset += count_universe_instances(*it, target_univ_id, univ_count_memo);
  }
  return offset;
}

//==============================================================================

void
Lattice::to_hdf5(hid_t lattices_group) const
{
  // Make a group for the lattice.
  std::string group_name {"lattice "};
  group_name += std::to_string(id_);
  hid_t lat_group = create_group(lattices_group, group_name);

  // Write the name and outer universe.
  if (!name_.empty()) {
    write_string(lat_group, "name", name_, false);
  }

  if (outer_ != NO_OUTER_UNIVERSE) {
    int32_t outer_id = model::universes[outer_]->id_;
    write_dataset(lat_group, "outer", outer_id);
  } else {
    write_dataset(lat_group, "outer", outer_);
  }

  // Call subclass-overriden function to fill in other details.
  to_hdf5_inner(lat_group);

  close_group(lat_group);
}

//==============================================================================
// RectLattice implementation
//==============================================================================

RectLattice::RectLattice(pugi::xml_node lat_node)
  : Lattice {lat_node}
{
  type_ = LatticeType::rect;

  // Read the number of lattice cells in each dimension.
  std::string dimension_str {get_node_value(lat_node, "dimension")};
  std::vector<std::string> dimension_words {split(dimension_str)};
  if (dimension_words.size() == 2) {
    n_cells_[0] = std::stoi(dimension_words[0]);
    n_cells_[1] = std::stoi(dimension_words[1]);
    n_cells_[2] = 1;
    is_3d_ = false;
  } else if (dimension_words.size() == 3) {
    n_cells_[0] = std::stoi(dimension_words[0]);
    n_cells_[1] = std::stoi(dimension_words[1]);
    n_cells_[2] = std::stoi(dimension_words[2]);
    is_3d_ = true;
  } else {
    fatal_error("Rectangular lattice must be two or three dimensions.");
  }

  // Read the lattice lower-left location.
  std::string ll_str {get_node_value(lat_node, "lower_left")};
  std::vector<std::string> ll_words {split(ll_str)};
  if (ll_words.size() != dimension_words.size()) {
    fatal_error("Number of entries on <lower_left> must be the same as the "
                "number of entries on <dimension>.");
  }
  lower_left_[0] = stod(ll_words[0]);
  lower_left_[1] = stod(ll_words[1]);
  if (is_3d_) {lower_left_[2] = stod(ll_words[2]);}

  // Read the lattice pitches.
  std::string pitch_str {get_node_value(lat_node, "pitch")};
  std::vector<std::string> pitch_words {split(pitch_str)};
  if (pitch_words.size() != dimension_words.size()) {
    fatal_error("Number of entries on <pitch> must be the same as the "
                "number of entries on <dimension>.");
  }
  pitch_[0] = stod(pitch_words[0]);
  pitch_[1] = stod(pitch_words[1]);
  if (is_3d_) {pitch_[2] = stod(pitch_words[2]);}

  // Read the universes and make sure the correct number was specified.
  std::string univ_str {get_node_value(lat_node, "universes")};
  std::vector<std::string> univ_words {split(univ_str)};
  if (univ_words.size() != nx*ny*nz) {
    fatal_error(fmt::format(
      "Expected {} universes for a rectangular lattice of size {}x{]x{} but {} "
      "were specified.", nx*ny*nz,  nx, ny, nz, univ_words.size()));
  }

  // Parse the universes.
  universes_.resize(nx*ny*nz, C_NONE);
  for (int iz = 0; iz < nz; iz++) {
    for (int iy = ny-1; iy > -1; iy--) {
      for (int ix = 0; ix < nx; ix++) {
        int indx1 = nx*ny*iz + nx*(ny-iy-1) + ix;
        int indx2 = nx*ny*iz + nx*iy + ix;
        universes_[indx1] = std::stoi(univ_words[indx2]);
      }
    }
  }
}

//==============================================================================

int32_t&
RectLattice::operator[](std::array<int, 3> i_xyz)
{
  int indx = nx*ny*i_xyz[2] + nx*i_xyz[1] + i_xyz[0];
  return universes_[indx];
}

//==============================================================================

bool
RectLattice::are_valid_indices(const int i_xyz[3]) const
{
  return (   (i_xyz[0] >= 0) && (i_xyz[0] < n_cells_[0])
          && (i_xyz[1] >= 0) && (i_xyz[1] < n_cells_[1])
          && (i_xyz[2] >= 0) && (i_xyz[2] < n_cells_[2]));
}

//==============================================================================

std::pair<double, std::array<int, 3>>
RectLattice::distance(Position r, Direction u, const std::array<int, 3>& i_xyz)
const
{
  // Get short aliases to the coordinates.
  double x = r.x;
  double y = r.y;
  double z = r.z;

  // Determine the oncoming edge.
  double x0 {copysign(0.5 * pitch_[0], u.x)};
  double y0 {copysign(0.5 * pitch_[1], u.y)};

  // Left and right sides
  double d {INFTY};
  std::array<int, 3> lattice_trans;
  if ((std::abs(x - x0) > FP_PRECISION) && u.x != 0) {
    d = (x0 - x) / u.x;
    if (u.x > 0) {
      lattice_trans = {1, 0, 0};
    } else {
      lattice_trans = {-1, 0, 0};
    }
  }

  // Front and back sides
  if ((std::abs(y - y0) > FP_PRECISION) && u.y != 0) {
    double this_d = (y0 - y) / u.y;
    if (this_d < d) {
      d = this_d;
      if (u.y > 0) {
        lattice_trans = {0, 1, 0};
      } else {
        lattice_trans = {0, -1, 0};
      }
    }
  }

  // Top and bottom sides
  if (is_3d_) {
    double z0 {copysign(0.5 * pitch_[2], u.z)};
    if ((std::abs(z - z0) > FP_PRECISION) && u.z != 0) {
      double this_d = (z0 - z) / u.z;
      if (this_d < d) {
        d = this_d;
        if (u.z > 0) {
          lattice_trans = {0, 0, 1};
        } else {
          lattice_trans = {0, 0, -1};
        }
      }
    }
  }

  return {d, lattice_trans};
}

//==============================================================================

std::array<int, 3>
RectLattice::get_indices(Position r, Direction u) const
{
  // Determine x index, accounting for coincidence
  double ix_ {(r.x - lower_left_.x) / pitch_.x};
  long ix_close {std::lround(ix_)};
  int ix;
  if (coincident(ix_, ix_close)) {
    ix = (u.x > 0) ? ix_close : ix_close - 1;
  } else {
    ix = std::floor(ix_);
  }

  // Determine y index, accounting for coincidence
  double iy_ {(r.y - lower_left_.y) / pitch_.y};
  long iy_close {std::lround(iy_)};
  int iy;
  if (coincident(iy_, iy_close)) {
    iy = (u.y > 0) ? iy_close : iy_close - 1;
  } else {
    iy = std::floor(iy_);
  }

  // Determine z index, accounting for coincidence
  int iz = 0;
  if (is_3d_) {
    double iz_ {(r.z - lower_left_.z) / pitch_.z};
    long iz_close {std::lround(iz_)};
    if (coincident(iz_, iz_close)) {
      iz = (u.z > 0) ? iz_close : iz_close - 1;
    } else {
      iz = std::floor(iz_);
    }
  }
  return {ix, iy, iz};
}

//==============================================================================

Position
RectLattice::get_local_position(Position r, const std::array<int, 3> i_xyz)
const
{
  r.x -= (lower_left_.x + (i_xyz[0] + 0.5)*pitch_.x);
  r.y -= (lower_left_.y + (i_xyz[1] + 0.5)*pitch_.y);
  if (is_3d_) {
    r.z -= (lower_left_.z + (i_xyz[2] + 0.5)*pitch_.z);
  }
  return r;
}

//==============================================================================

int32_t&
RectLattice::offset(int map, const int i_xyz[3])
{
  return offsets_[nx*ny*nz*map + nx*ny*i_xyz[2] + nx*i_xyz[1] + i_xyz[0]];
}

//==============================================================================

std::string
RectLattice::index_to_string(int indx) const
{
  int iz {indx / (nx * ny)};
  int iy {(indx - nx * ny * iz) / nx};
  int ix {indx - nx * ny * iz - nx * iy};
  std::string out {std::to_string(ix)};
  out += ',';
  out += std::to_string(iy);
  if (is_3d_) {
    out += ',';
    out += std::to_string(iz);
  }
  return out;
}

//==============================================================================

void
RectLattice::to_hdf5_inner(hid_t lat_group) const
{
  // Write basic lattice information.
  write_string(lat_group, "type", "rectangular", false);
  if (is_3d_) {
    write_dataset(lat_group, "pitch", pitch_);
    write_dataset(lat_group, "lower_left", lower_left_);
    write_dataset(lat_group, "dimension", n_cells_);
  } else {
    std::array<double, 2> pitch_short {{pitch_[0], pitch_[1]}};
    write_dataset(lat_group, "pitch", pitch_short);
    std::array<double, 2> ll_short {{lower_left_[0], lower_left_[1]}};
    write_dataset(lat_group, "lower_left", ll_short);
    std::array<int, 2> nc_short {{n_cells_[0], n_cells_[1]}};
    write_dataset(lat_group, "dimension", nc_short);
  }

  // Write the universe ids.  The convention here is to switch the ordering on
  // the y-axis to match the way universes are input in a text file.
  if (is_3d_) {
    hsize_t nx {static_cast<hsize_t>(n_cells_[0])};
    hsize_t ny {static_cast<hsize_t>(n_cells_[1])};
    hsize_t nz {static_cast<hsize_t>(n_cells_[2])};
    std::vector<int> out(nx*ny*nz);

    for (int m = 0; m < nz; m++) {
      for (int k = 0; k < ny; k++) {
        for (int j = 0; j < nx; j++) {
          int indx1 = nx*ny*m + nx*k + j;
          int indx2 = nx*ny*m + nx*(ny-k-1) + j;
          out[indx2] = model::universes[universes_[indx1]]->id_;
        }
      }
    }

    hsize_t dims[3] {nz, ny, nx};
    write_int(lat_group, 3, dims, "universes", out.data(), false);

  } else {
    hsize_t nx {static_cast<hsize_t>(n_cells_[0])};
    hsize_t ny {static_cast<hsize_t>(n_cells_[1])};
    std::vector<int> out(nx*ny);

    for (int k = 0; k < ny; k++) {
      for (int j = 0; j < nx; j++) {
        int indx1 = nx*k + j;
        int indx2 = nx*(ny-k-1) + j;
        out[indx2] = model::universes[universes_[indx1]]->id_;
      }
    }

    hsize_t dims[3] {1, ny, nx};
    write_int(lat_group, 3, dims, "universes", out.data(), false);
  }
}

//==============================================================================
// HexLattice implementation
//==============================================================================

HexLattice::HexLattice(pugi::xml_node lat_node)
  : Lattice {lat_node}
{
  type_ = LatticeType::hex;

  // Read the number of lattice cells in each dimension.
  n_rings_ = std::stoi(get_node_value(lat_node, "n_rings"));
  if (check_for_node(lat_node, "n_axial")) {
    n_axial_ = std::stoi(get_node_value(lat_node, "n_axial"));
    is_3d_ = true;
  } else {
    n_axial_ = 1;
    is_3d_ = false;
  }

  // Read the lattice orientation.  Default to 'y'.
  if (check_for_node(lat_node, "orientation")) {
    std::string orientation = get_node_value(lat_node, "orientation");
    if (orientation == "y") {
      orientation_ = Orientation::y;
    } else if (orientation == "x") {
      orientation_ = Orientation::x;
    } else {
      fatal_error("Unrecognized orientation '" + orientation
                  + "' for lattice " + std::to_string(id_));
    }
  } else {
    orientation_ = Orientation::y;
  }

  // Read the lattice center.
  std::string center_str {get_node_value(lat_node, "center")};
  std::vector<std::string> center_words {split(center_str)};
  if (is_3d_ && (center_words.size() != 3)) {
    fatal_error("A hexagonal lattice with <n_axial> must have <center> "
                "specified by 3 numbers.");
  } else if (!is_3d_ && center_words.size() != 2) {
    fatal_error("A hexagonal lattice without <n_axial> must have <center> "
                "specified by 2 numbers.");
  }
  center_[0] = stod(center_words[0]);
  center_[1] = stod(center_words[1]);
  if (is_3d_) {center_[2] = stod(center_words[2]);}

  // Read the lattice pitches.
  std::string pitch_str {get_node_value(lat_node, "pitch")};
  std::vector<std::string> pitch_words {split(pitch_str)};
  if (is_3d_ && (pitch_words.size() != 2)) {
    fatal_error("A hexagonal lattice with <n_axial> must have <pitch> "
                "specified by 2 numbers.");
  } else if (!is_3d_ && (pitch_words.size() != 1)) {
    fatal_error("A hexagonal lattice without <n_axial> must have <pitch> "
                "specified by 1 number.");
  }
  pitch_[0] = stod(pitch_words[0]);
  if (is_3d_) {pitch_[1] = stod(pitch_words[1]);}

  // Read the universes and make sure the correct number was specified.
  int n_univ = (3*n_rings_*n_rings_ - 3*n_rings_ + 1) * n_axial_;
  std::string univ_str {get_node_value(lat_node, "universes")};
  std::vector<std::string> univ_words {split(univ_str)};
  if (univ_words.size() != n_univ) {
    fatal_error(fmt::format(
      "Expected {} universes for a hexagonal lattice with {} rings and {} "
      "axial levels but {} were specified.", n_univ, n_rings_, n_axial_,
      univ_words.size()));
  }

  // Parse the universes.
  // Universes in hexagonal lattices are stored in a manner that represents
  // a skewed coordinate system: (x, alpha) in case of 'y' orientation
  // and (alpha,y) in 'x' one rather than (x, y).  There is
  // no obvious, direct relationship between the order of universes in the
  // input and the order that they will be stored in the skewed array so
  // the following code walks a set of index values across the skewed array
  // in a manner that matches the input order.  Note that i_x = 0, i_a = 0
  // or i_a = 0, i_y = 0 corresponds to the center of the hexagonal lattice.
  universes_.resize((2*n_rings_-1) * (2*n_rings_-1) * n_axial_, C_NONE);
  if (orientation_ == Orientation::y) {
    fill_lattice_y(univ_words);
  } else {
    fill_lattice_x(univ_words);
  }
}

//==============================================================================

void
HexLattice::fill_lattice_x(const std::vector<std::string>& univ_words)
{
  int input_index = 0;
  for (int m = 0; m < n_axial_; m++) {
    // Initialize lattice indecies.
    int i_a = -(n_rings_ - 1);
    int i_y = n_rings_ - 1;

    // Map upper region of hexagonal lattice which is found in the
    // first n_rings-1 rows of the input.
    for (int k = 0; k < n_rings_-1; k++) {

      // Iterate over the input columns.
      for (int j = 0; j < k+n_rings_; j++) {
        int indx = (2*n_rings_-1)*(2*n_rings_-1) * m
                    + (2*n_rings_-1) * (i_y+n_rings_-1)
                    + (i_a+n_rings_-1);
        universes_[indx] = std::stoi(univ_words[input_index]);
        input_index++;
        // Move to the next right neighbour cell
        i_a += 1;
      }

      // Return the lattice index to the start of the current row.
      i_a = -(n_rings_ - 1);
      i_y -= 1;
    }

    // Map the lower region from the centerline of cart to down side
    for (int k = 0; k < n_rings_; k++) {
      // Walk the index to the lower-right neighbor of the last row start.
      i_a = -(n_rings_ - 1) + k;

      // Iterate over the input columns.
      for (int j = 0; j < 2*n_rings_-k-1; j++) {
        int indx = (2*n_rings_-1)*(2*n_rings_-1) * m
                    + (2*n_rings_-1) * (i_y+n_rings_-1)
                    + (i_a+n_rings_-1);
        universes_[indx] = std::stoi(univ_words[input_index]);
        input_index++;
        // Move to the next right neighbour cell
        i_a += 1;
      }

      // Return lattice index to start of current row.
      i_y -= 1;
    }
  }
}

//==============================================================================

void
HexLattice::fill_lattice_y(const std::vector<std::string>& univ_words)
{
  int input_index = 0;
  for (int m = 0; m < n_axial_; m++) {
    // Initialize lattice indecies.
    int i_x = 1;
    int i_a = n_rings_ - 1;

    // Map upper triangular region of hexagonal lattice which is found in the
    // first n_rings-1 rows of the input.
    for (int k = 0; k < n_rings_-1; k++) {
      // Walk the index to lower-left neighbor of last row start.
      i_x -= 1;

      // Iterate over the input columns.
      for (int j = 0; j < k+1; j++) {
        int indx = (2*n_rings_-1)*(2*n_rings_-1) * m
                    + (2*n_rings_-1) * (i_a+n_rings_-1)
                    + (i_x+n_rings_-1);
        universes_[indx] = std::stoi(univ_words[input_index]);
        input_index++;
        // Walk the index to the right neighbor (which is not adjacent).
        i_x += 2;
        i_a -= 1;
      }

      // Return the lattice index to the start of the current row.
      i_x -= 2 * (k+1);
      i_a += (k+1);
    }

    // Map the middle square region of the hexagonal lattice which is found in
    // the next 2*n_rings-1 rows of the input.
    for (int k = 0; k < 2*n_rings_-1; k++) {
      if ((k % 2) == 0) {
        // Walk the index to the lower-left neighbor of the last row start.
        i_x -= 1;
      } else {
        // Walk the index to the lower-right neighbor of the last row start.
        i_x += 1;
        i_a -= 1;
      }

      // Iterate over the input columns.
      for (int j = 0; j < n_rings_ - (k % 2); j++) {
        int indx = (2*n_rings_-1)*(2*n_rings_-1) * m
                    + (2*n_rings_-1) * (i_a+n_rings_-1)
                    + (i_x+n_rings_-1);
        universes_[indx] = std::stoi(univ_words[input_index]);
        input_index++;
        // Walk the index to the right neighbor (which is not adjacent).
        i_x += 2;
        i_a -= 1;
      }

      // Return the lattice index to the start of the current row.
      i_x -= 2*(n_rings_ - (k % 2));
      i_a += n_rings_ - (k % 2);
    }

    // Map the lower triangular region of the hexagonal lattice.
    for (int k = 0; k < n_rings_-1; k++) {
      // Walk the index to the lower-right neighbor of the last row start.
      i_x += 1;
      i_a -= 1;

      // Iterate over the input columns.
      for (int j = 0; j < n_rings_-k-1; j++) {
        int indx = (2*n_rings_-1)*(2*n_rings_-1) * m
                    + (2*n_rings_-1) * (i_a+n_rings_-1)
                    + (i_x+n_rings_-1);
        universes_[indx] = std::stoi(univ_words[input_index]);
        input_index++;
        // Walk the index to the right neighbor (which is not adjacent).
        i_x += 2;
        i_a -= 1;
      }

      // Return lattice index to start of current row.
      i_x -= 2*(n_rings_ - k - 1);
      i_a += n_rings_ - k - 1;
    }
  }
}

//==============================================================================

int32_t&
HexLattice::operator[](std::array<int, 3> i_xyz)
{
  int indx = (2*n_rings_-1)*(2*n_rings_-1) * i_xyz[2]
              + (2*n_rings_-1) * i_xyz[1]
              + i_xyz[0];
  return universes_[indx];
}

//==============================================================================

LatticeIter HexLattice::begin()
{return LatticeIter(*this, n_rings_-1);}

ReverseLatticeIter HexLattice::rbegin()
{return ReverseLatticeIter(*this, universes_.size()-n_rings_);}

//==============================================================================

bool
HexLattice::are_valid_indices(const int i_xyz[3]) const
{
  return ((i_xyz[0] >= 0) && (i_xyz[1] >= 0) && (i_xyz[2] >= 0)
          && (i_xyz[0] < 2*n_rings_-1) && (i_xyz[1] < 2*n_rings_-1)
          && (i_xyz[0] + i_xyz[1] > n_rings_-2)
          && (i_xyz[0] + i_xyz[1] < 3*n_rings_-2)
          && (i_xyz[2] < n_axial_));
}

//==============================================================================

std::pair<double, std::array<int, 3>>
HexLattice::distance(Position r, Direction u, const std::array<int, 3>& i_xyz)
const
{
  // Short description of the direction vectors used here.  The beta, gamma, and
  // delta vectors point towards the flat sides of each hexagonal tile.
  // Y - orientation:
  //   basis0 = (1, 0)
  //   basis1 = (-1/sqrt(3), 1)   = +120 degrees from basis0
  //   beta   = (sqrt(3)/2, 1/2)  = +30 degrees from basis0
  //   gamma  = (sqrt(3)/2, -1/2) = -60 degrees from beta
  //   delta  = (0, 1)            = +60 degrees from beta
  // X - orientation:
  //   basis0 = (1/sqrt(3), -1)
  //   basis1 = (0, 1)            = +120 degrees from basis0
  //   beta   = (1, 0)            = +30 degrees from basis0
  //   gamma  = (1/2, -sqrt(3)/2) = -60 degrees from beta
  //   delta  = (1/2, sqrt(3)/2)  = +60 degrees from beta
  // The z-axis is considered separately.
  double beta_dir;
  double gamma_dir;
  double delta_dir;
  if (orientation_ == Orientation::y) {
    beta_dir = u.x * std::sqrt(3.0) / 2.0  + u.y / 2.0;
    gamma_dir = u.x * std::sqrt(3.0) / 2.0  - u.y / 2.0;
    delta_dir = u.y;
  } else {
    beta_dir = u.x;
    gamma_dir = u.x / 2.0  - u.y * std::sqrt(3.0) / 2.0;
    delta_dir = u.x / 2.0  + u.y * std::sqrt(3.0) / 2.0;
  }

  // Note that hexagonal lattice distance calculations are performed
  // using the particle's coordinates relative to the neighbor lattice
  // cells, not relative to the particle's current cell.  This is done
  // because there is significant disagreement between neighboring cells
  // on where the lattice boundary is due to finite precision issues.

  // beta direction
  double d {INFTY};
  std::array<int, 3> lattice_trans;
  double edge = -copysign(0.5*pitch_[0], beta_dir);  // Oncoming edge
  Position r_t;
  if (beta_dir > 0) {
    const std::array<int, 3> i_xyz_t {i_xyz[0]+1, i_xyz[1], i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  } else {
    const std::array<int, 3> i_xyz_t {i_xyz[0]-1, i_xyz[1], i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  }
  double beta;
  if (orientation_ == Orientation::y) {
    beta = r_t.x * std::sqrt(3.0) / 2.0 + r_t.y / 2.0;
  } else {
    beta = r_t.x;
  }
  if ((std::abs(beta - edge) > FP_PRECISION) && beta_dir != 0) {
    d = (edge - beta) / beta_dir;
    if (beta_dir > 0) {
      lattice_trans = {1, 0, 0};
    } else {
      lattice_trans = {-1, 0, 0};
    }
  }

  // gamma direction
  edge = -copysign(0.5*pitch_[0], gamma_dir);
  if (gamma_dir > 0) {
    const std::array<int, 3> i_xyz_t {i_xyz[0]+1, i_xyz[1]-1, i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  } else {
    const std::array<int, 3> i_xyz_t {i_xyz[0]-1, i_xyz[1]+1, i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  }
  double gamma;
  if (orientation_ == Orientation::y) {
    gamma = r_t.x * std::sqrt(3.0) / 2.0 - r_t.y / 2.0;
  } else {
    gamma = r_t.x  / 2.0 - r_t.y * std::sqrt(3.0) / 2.0;
  }
  if ((std::abs(gamma - edge) > FP_PRECISION) && gamma_dir != 0) {
    double this_d = (edge - gamma) / gamma_dir;
    if (this_d < d) {
      if (gamma_dir > 0) {
        lattice_trans = {1, -1, 0};
      } else {
        lattice_trans = {-1, 1, 0};
      }
      d = this_d;
    }
  }

  // delta direction
  edge = -copysign(0.5*pitch_[0], delta_dir);
  if (delta_dir > 0) {
    const std::array<int, 3> i_xyz_t {i_xyz[0], i_xyz[1]+1, i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  } else {
    const std::array<int, 3> i_xyz_t {i_xyz[0], i_xyz[1]-1, i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  }
  double delta;
  if (orientation_ == Orientation::y) {
    delta =  r_t.y;
  } else {
    delta = r_t.x  / 2.0 + r_t.y * std::sqrt(3.0) / 2.0;
  }
  if ((std::abs(delta - edge) > FP_PRECISION) && delta_dir != 0) {
    double this_d = (edge - delta) / delta_dir;
    if (this_d < d) {
      if (delta_dir > 0) {
        lattice_trans = {0, 1, 0};
      } else {
        lattice_trans = {0, -1, 0};
      }
      d = this_d;
    }
  }

  // Top and bottom sides
  if (is_3d_) {
    double z = r.z;
    double z0 {copysign(0.5 * pitch_[1], u.z)};
    if ((std::abs(z - z0) > FP_PRECISION) && u.z != 0) {
      double this_d = (z0 - z) / u.z;
      if (this_d < d) {
        d = this_d;
        if (u.z > 0) {
          lattice_trans = {0, 0, 1};
        } else {
          lattice_trans = {0, 0, -1};
        }
        d = this_d;
      }
    }
  }

  return {d, lattice_trans};
}

//==============================================================================

std::array<int, 3>
HexLattice::get_indices(Position r, Direction u) const
{
  // Offset the xyz by the lattice center.
  Position r_o {r.x - center_.x, r.y - center_.y, r.z};
  if (is_3d_) {r_o.z -= center_.z;}

  // Index the z direction, accounting for coincidence
  int iz = 0;
  if (is_3d_) {
    double iz_ {r_o.z / pitch_[1] + 0.5 * n_axial_};
    long iz_close {std::lround(iz_)};
    if (coincident(iz_, iz_close)) {
      iz = (u.z > 0) ? iz_close : iz_close - 1;
    } else {
      iz = std::floor(iz_);
    }
  }

  int i1, i2;
  if (orientation_ == Orientation::y) {
    // Convert coordinates into skewed bases.  The (x, alpha) basis is used to
    // find the index of the global coordinates to within 4 cells.
    double alpha = r_o.y - r_o.x / std::sqrt(3.0);
    i1 = std::floor(r_o.x / (0.5*std::sqrt(3.0) * pitch_[0]));
    i2 = std::floor(alpha / pitch_[0]);
  } else {
    // Convert coordinates into skewed bases.  The (alpha, y) basis is used to
    // find the index of the global coordinates to within 4 cells.
    double alpha = r_o.y - r_o.x * std::sqrt(3.0);
    i1 = std::floor(-alpha / (std::sqrt(3.0) * pitch_[0]));
    i2 = std::floor(r_o.y / (0.5*std::sqrt(3.0) * pitch_[0]));
  }

  // Add offset to indices (the center cell is (i1, i2) = (0, 0) but
  // the array is offset so that the indices never go below 0).
  i1 += n_rings_-1;
  i2 += n_rings_-1;

  // Calculate the (squared) distance between the particle and the centers of
  // the four possible cells.  Regular hexagonal tiles form a Voronoi
  // tessellation so the xyz should be in the hexagonal cell that it is closest
  // to the center of.  This method is used over a method that uses the
  // remainders of the floor divisions above because it provides better finite
  // precision performance.  Squared distances are used because they are more
  // computationally efficient than normal distances.

  // COINCIDENCE CHECK
  // if a distance to center, d, is within the coincidence tolerance of the
  // current minimum distance, d_min, the particle is on an edge or vertex.
  // In this case, the dot product of the position vector and direction vector
  // for the current indices, dp, and the dot product for the currently selected
  // indices, dp_min, are compared. The cell which the particle is moving into
  // is kept (i.e. the cell with the lowest dot product as the vectors will be
  // completely opposed if the particle is moving directly toward the center of
  // the cell).
  int i1_chg {};
  int i2_chg {};
  double d_min {INFTY};
  double dp_min {INFTY};
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      // get local coordinates
      const std::array<int, 3> i_xyz {i1 + j, i2 + i, 0};
      Position r_t = get_local_position(r, i_xyz);
      // calculate distance
      double d = r_t.x*r_t.x + r_t.y*r_t.y;
      // check for coincidence. Because the numerical error incurred
      // in hex geometry is higher than other geometries, the relative
      // coincidence is checked here so that coincidence is successfully
      // detected on large hex lattice with particles far from the origin
      // which have rounding errors larger than the FP_COINCIDENT thresdhold.
      bool on_boundary = coincident(1.0, d_min/d);
      if (d < d_min || on_boundary) {
        // normalize r_t and find dot product
        r_t /= std::sqrt(d);
        double dp = u.x * r_t.x + u.y * r_t.y;
        // do not update values if particle is on a
        // boundary and not moving into this cell
        if (on_boundary && dp > dp_min) continue;
        // update values
        d_min = d;
        i1_chg = j;
        i2_chg = i;
        dp_min = dp;
      }
    }
  }

  // update outgoing indices
  i1 += i1_chg;
  i2 += i2_chg;

  return {i1, i2, iz};
}

//==============================================================================

Position
HexLattice::get_local_position(Position r, const std::array<int, 3> i_xyz)
const
{
  if (orientation_ == Orientation::y) {
    // x_l = x_g - (center + pitch_x*cos(30)*index_x)
    r.x -= center_.x
           + std::sqrt(3.0)/2.0 * (i_xyz[0] - n_rings_ + 1) * pitch_[0];
    // y_l = y_g - (center + pitch_x*index_x + pitch_y*sin(30)*index_y)
    r.y -= (center_.y + (i_xyz[1] - n_rings_ + 1) * pitch_[0]
            + (i_xyz[0] - n_rings_ + 1) * pitch_[0] / 2.0);
  } else {
    // x_l = x_g - (center + pitch_x*index_a + pitch_y*sin(30)*index_y)
    r.x -= (center_.x + (i_xyz[0] - n_rings_ + 1) * pitch_[0]
            + (i_xyz[1] - n_rings_ + 1) * pitch_[0] / 2.0);
    // y_l = y_g - (center + pitch_y*cos(30)*index_y)
    r.y -= center_.y
           + std::sqrt(3.0)/2.0 * (i_xyz[1] - n_rings_ + 1) * pitch_[0];
  }

  if (is_3d_) {
      r.z -= center_.z - (0.5 * n_axial_ - i_xyz[2] - 0.5) * pitch_[1];
    }

  return r;
}

//==============================================================================

bool
HexLattice::is_valid_index(int indx) const
{
  int nx {2*n_rings_ - 1};
  int ny {2*n_rings_ - 1};
  int iz = indx / (nx * ny);
  int iy = (indx - nx*ny*iz) / nx;
  int ix = indx - nx*ny*iz - nx*iy;
  int i_xyz[3] {ix, iy, iz};
  return are_valid_indices(i_xyz);
}

//==============================================================================

int32_t&
HexLattice::offset(int map, const int i_xyz[3])
{
  int nx {2*n_rings_ - 1};
  int ny {2*n_rings_ - 1};
  int nz {n_axial_};
  return offsets_[nx*ny*nz*map + nx*ny*i_xyz[2] + nx*i_xyz[1] + i_xyz[0]];
}

//==============================================================================

std::string
HexLattice::index_to_string(int indx) const
{
  int nx {2*n_rings_ - 1};
  int ny {2*n_rings_ - 1};
  int iz {indx / (nx * ny)};
  int iy {(indx - nx * ny * iz) / nx};
  int ix {indx - nx * ny * iz - nx * iy};
  std::string out {std::to_string(ix - n_rings_ + 1)};
  out += ',';
  out += std::to_string(iy - n_rings_ + 1);
  if (is_3d_) {
    out += ',';
    out += std::to_string(iz);
  }
  return out;
}

//==============================================================================

void
HexLattice::to_hdf5_inner(hid_t lat_group) const
{
  // Write basic lattice information.
  write_string(lat_group, "type", "hexagonal", false);
  write_dataset(lat_group, "n_rings", n_rings_);
  write_dataset(lat_group, "n_axial", n_axial_);
  if (orientation_ == Orientation::y) {
    write_string(lat_group, "orientation", "y", false);
  } else {
    write_string(lat_group, "orientation", "x", false);
  }
  if (is_3d_) {
    write_dataset(lat_group, "pitch", pitch_);
    write_dataset(lat_group, "center", center_);
  } else {
    std::array<double, 1> pitch_short {{pitch_[0]}};
    write_dataset(lat_group, "pitch", pitch_short);
    std::array<double, 2> center_short {{center_[0], center_[1]}};
    write_dataset(lat_group, "center", center_short);
  }

  // Write the universe ids.
  hsize_t nx {static_cast<hsize_t>(2*n_rings_ - 1)};
  hsize_t ny {static_cast<hsize_t>(2*n_rings_ - 1)};
  hsize_t nz {static_cast<hsize_t>(n_axial_)};
  std::vector<int> out(nx*ny*nz);

  for (int m = 0; m < nz; m++) {
    for (int k = 0; k < ny; k++) {
      for (int j = 0; j < nx; j++) {
        int indx = nx*ny*m + nx*k + j;
        if (j + k < n_rings_ - 1) {
          // This array position is never used; put a -1 to indicate this.
          out[indx] = -1;
        } else if (j + k > 3*n_rings_ - 3) {
          // This array position is never used; put a -1 to indicate this.
          out[indx] = -1;
        } else {
          out[indx] = model::universes[universes_[indx]]->id_;
        }
      }
    }
  }

  hsize_t dims[3] {nz, ny, nx};
  write_int(lat_group, 3, dims, "universes", out.data(), false);
}

//==============================================================================
// Non-method functions
//==============================================================================

void read_lattices(pugi::xml_node node)
{
  for (pugi::xml_node lat_node : node.children("lattice")) {
    model::lattices.push_back(std::make_unique<RectLattice>(lat_node));
  }
  for (pugi::xml_node lat_node : node.children("hex_lattice")) {
    model::lattices.push_back(std::make_unique<HexLattice>(lat_node));
  }

  // Fill the lattice map.
  for (int i_lat = 0; i_lat < model::lattices.size(); i_lat++) {
    int id = model::lattices[i_lat]->id_;
    auto in_map = model::lattice_map.find(id);
    if (in_map == model::lattice_map.end()) {
      model::lattice_map[id] = i_lat;
    } else {
      fatal_error(fmt::format(
        "Two or more lattices use the same unique ID: {}", id));
    }
  }
}

} // namespace openmc
