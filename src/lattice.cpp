#include "openmc/lattice.h"

#include <cmath>
#include <sstream>
#include <vector>

#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry_aux.h"
#include "openmc/hdf5_interface.h"
#include "openmc/string_utils.h"
#include "openmc/xml_interface.h"


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<Lattice*> lattices_c;

std::unordered_map<int32_t, int32_t> lattice_map;

//==============================================================================
// Lattice implementation
//==============================================================================

Lattice::Lattice(pugi::xml_node lat_node)
{
  if (check_for_node(lat_node, "id")) {
    id = std::stoi(get_node_value(lat_node, "id"));
  } else {
    fatal_error("Must specify id of lattice in geometry XML file.");
  }

  if (check_for_node(lat_node, "name")) {
    name = get_node_value(lat_node, "name");
  }

  if (check_for_node(lat_node, "outer")) {
    outer = std::stoi(get_node_value(lat_node, "outer"));
  }
}

//==============================================================================

LatticeIter Lattice::begin()
{return LatticeIter(*this, 0);}

LatticeIter Lattice::end()
{return LatticeIter(*this, universes.size());}

ReverseLatticeIter Lattice::rbegin()
{return ReverseLatticeIter(*this, universes.size()-1);}

ReverseLatticeIter Lattice::rend()
{return ReverseLatticeIter(*this, -1);}

//==============================================================================

void
Lattice::adjust_indices()
{
  // Adjust the indices for the universes array.
  for (LatticeIter it = begin(); it != end(); ++it) {
    int uid = *it;
    auto search = universe_map.find(uid);
    if (search != universe_map.end()) {
      *it = search->second;
    } else {
      std::stringstream err_msg;
      err_msg << "Invalid universe number " << uid << " specified on "
           "lattice " << id;
      fatal_error(err_msg);
    }
  }

  // Adjust the index for the outer universe.
  if (outer != NO_OUTER_UNIVERSE) {
    auto search = universe_map.find(outer);
    if (search != universe_map.end()) {
      outer = search->second;
    } else {
      std::stringstream err_msg;
      err_msg << "Invalid universe number " << outer << " specified on "
           "lattice " << id;
      fatal_error(err_msg);
    }
  }
}

//==============================================================================

int32_t
Lattice::fill_offset_table(int32_t offset, int32_t target_univ_id, int map)
{
  for (LatticeIter it = begin(); it != end(); ++it) {
    offsets[map * universes.size() + it.indx] = offset;
    offset += count_universe_instances(*it, target_univ_id);
  }
  return offset;
}

//==============================================================================

void
Lattice::to_hdf5(hid_t lattices_group) const
{
  // Make a group for the lattice.
  std::string group_name {"lattice "};
  group_name += std::to_string(id);
  hid_t lat_group = create_group(lattices_group, group_name);

  // Write the name and outer universe.
  if (!name.empty()) {
    write_string(lat_group, "name", name, false);
  }

  if (outer != NO_OUTER_UNIVERSE) {
    int32_t outer_id = global_universes[outer]->id;
    write_dataset(lat_group, "outer", outer_id);
  } else {
    write_dataset(lat_group, "outer", outer);
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
  // Read the number of lattice cells in each dimension.
  std::string dimension_str {get_node_value(lat_node, "dimension")};
  std::vector<std::string> dimension_words {split(dimension_str)};
  if (dimension_words.size() == 2) {
    n_cells[0] = std::stoi(dimension_words[0]);
    n_cells[1] = std::stoi(dimension_words[1]);
    n_cells[2] = 1;
    is_3d = false;
  } else if (dimension_words.size() == 3) {
    n_cells[0] = std::stoi(dimension_words[0]);
    n_cells[1] = std::stoi(dimension_words[1]);
    n_cells[2] = std::stoi(dimension_words[2]);
    is_3d = true;
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
  lower_left[0] = stod(ll_words[0]);
  lower_left[1] = stod(ll_words[1]);
  if (is_3d) {lower_left[2] = stod(ll_words[2]);}

  // Read the lattice pitches.
  std::string pitch_str {get_node_value(lat_node, "pitch")};
  std::vector<std::string> pitch_words {split(pitch_str)};
  if (pitch_words.size() != dimension_words.size()) {
    fatal_error("Number of entries on <pitch> must be the same as the "
                "number of entries on <dimension>.");
  }
  pitch[0] = stod(pitch_words[0]);
  pitch[1] = stod(pitch_words[1]);
  if (is_3d) {pitch[2] = stod(pitch_words[2]);}

  // Read the universes and make sure the correct number was specified.
  std::string univ_str {get_node_value(lat_node, "universes")};
  std::vector<std::string> univ_words {split(univ_str)};
  if (univ_words.size() != nx*ny*nz) {
    std::stringstream err_msg;
    err_msg << "Expected " << nx*ny*nz
            << " universes for a rectangular lattice of size "
            << nx << "x" << ny << "x" << nz << " but " << univ_words.size()
            << " were specified.";
    fatal_error(err_msg);
  }

  // Parse the universes.
  universes.resize(nx*ny*nz, C_NONE);
  for (int iz = 0; iz < nz; iz++) {
    for (int iy = ny-1; iy > -1; iy--) {
      for (int ix = 0; ix < nx; ix++) {
        int indx1 = nx*ny*iz + nx*(ny-iy-1) + ix;
        int indx2 = nx*ny*iz + nx*iy + ix;
        universes[indx1] = std::stoi(univ_words[indx2]);
      }
    }
  }
}

//==============================================================================

int32_t&
RectLattice::operator[](const int i_xyz[3])
{
  int indx = nx*ny*i_xyz[2] + nx*i_xyz[1] + i_xyz[0];
  return universes[indx];
}

//==============================================================================

bool
RectLattice::are_valid_indices(const int i_xyz[3]) const
{
  return (   (i_xyz[0] >= 0) && (i_xyz[0] < n_cells[0])
          && (i_xyz[1] >= 0) && (i_xyz[1] < n_cells[1])
          && (i_xyz[2] >= 0) && (i_xyz[2] < n_cells[2]));
}

//==============================================================================

std::pair<double, std::array<int, 3>>
RectLattice::distance(Position r, Direction u, const int i_xyz[3]) const
{
  // Get short aliases to the coordinates.
  double x = r.x;
  double y = r.y;
  double z = r.z;

  // Determine the oncoming edge.
  double x0 {copysign(0.5 * pitch[0], u.x)};
  double y0 {copysign(0.5 * pitch[1], u.y)};

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
  if (is_3d) {
    double z0 {copysign(0.5 * pitch[2], u.z)};
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
RectLattice::get_indices(Position r) const
{
  int ix {static_cast<int>(std::ceil((r.x - lower_left.x) / pitch.x))-1};
  int iy {static_cast<int>(std::ceil((r.y - lower_left.y) / pitch.y))-1};
  int iz;
  if (is_3d) {
    iz = static_cast<int>(std::ceil((r.z - lower_left.z) / pitch.z))-1;
  } else {
    iz = 0;
  }
  return {ix, iy, iz};
}

//==============================================================================

Position
RectLattice::get_local_position(Position r, const int i_xyz[3]) const
{
  r.x -= (lower_left.x + (i_xyz[0] + 0.5)*pitch.x);
  r.y -= (lower_left.y + (i_xyz[1] + 0.5)*pitch.y);
  if (is_3d) {
    r.z -= (lower_left.z + (i_xyz[2] + 0.5)*pitch.z);
  }
  return r;
}

//==============================================================================

int32_t&
RectLattice::offset(int map, const int i_xyz[3])
{
  return offsets[nx*ny*nz*map + nx*ny*i_xyz[2] + nx*i_xyz[1] + i_xyz[0]];
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
  if (is_3d) {
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
  if (is_3d) {
    write_dataset(lat_group, "pitch", pitch);
    write_dataset(lat_group, "lower_left", lower_left);
    write_dataset(lat_group, "dimension", n_cells);
  } else {
    std::array<double, 2> pitch_short {{pitch[0], pitch[1]}};
    write_dataset(lat_group, "pitch", pitch_short);
    std::array<double, 2> ll_short {{lower_left[0], lower_left[1]}};
    write_dataset(lat_group, "lower_left", ll_short);
    std::array<int, 2> nc_short {{n_cells[0], n_cells[1]}};
    write_dataset(lat_group, "dimension", nc_short);
  }

  // Write the universe ids.  The convention here is to switch the ordering on
  // the y-axis to match the way universes are input in a text file.
  if (is_3d) {
    hsize_t nx {static_cast<hsize_t>(n_cells[0])};
    hsize_t ny {static_cast<hsize_t>(n_cells[1])};
    hsize_t nz {static_cast<hsize_t>(n_cells[2])};
    int out[nx*ny*nz];

    for (int m = 0; m < nz; m++) {
      for (int k = 0; k < ny; k++) {
        for (int j = 0; j < nx; j++) {
          int indx1 = nx*ny*m + nx*k + j;
          int indx2 = nx*ny*m + nx*(ny-k-1) + j;
          out[indx2] = global_universes[universes[indx1]]->id;
        }
      }
    }

    hsize_t dims[3] {nz, ny, nx};
    write_int(lat_group, 3, dims, "universes", out, false);

  } else {
    hsize_t nx {static_cast<hsize_t>(n_cells[0])};
    hsize_t ny {static_cast<hsize_t>(n_cells[1])};
    int out[nx*ny];

    for (int k = 0; k < ny; k++) {
      for (int j = 0; j < nx; j++) {
        int indx1 = nx*k + j;
        int indx2 = nx*(ny-k-1) + j;
        out[indx2] = global_universes[universes[indx1]]->id;
      }
    }

    hsize_t dims[3] {1, ny, nx};
    write_int(lat_group, 3, dims, "universes", out, false);
  }
}

//==============================================================================
// HexLattice implementation
//==============================================================================

HexLattice::HexLattice(pugi::xml_node lat_node)
  : Lattice {lat_node}
{
  // Read the number of lattice cells in each dimension.
  n_rings = std::stoi(get_node_value(lat_node, "n_rings"));
  if (check_for_node(lat_node, "n_axial")) {
    n_axial = std::stoi(get_node_value(lat_node, "n_axial"));
    is_3d = true;
  } else {
    n_axial = 1;
    is_3d = false;
  }

  // Read the lattice center.
  std::string center_str {get_node_value(lat_node, "center")};
  std::vector<std::string> center_words {split(center_str)};
  if (is_3d && (center_words.size() != 3)) {
    fatal_error("A hexagonal lattice with <n_axial> must have <center> "
                "specified by 3 numbers.");
  } else if (!is_3d && center_words.size() != 2) {
    fatal_error("A hexagonal lattice without <n_axial> must have <center> "
                "specified by 2 numbers.");
  }
  center[0] = stod(center_words[0]);
  center[1] = stod(center_words[1]);
  if (is_3d) {center[2] = stod(center_words[2]);}

  // Read the lattice pitches.
  std::string pitch_str {get_node_value(lat_node, "pitch")};
  std::vector<std::string> pitch_words {split(pitch_str)};
  if (is_3d && (pitch_words.size() != 2)) {
    fatal_error("A hexagonal lattice with <n_axial> must have <pitch> "
                "specified by 2 numbers.");
  } else if (!is_3d && (pitch_words.size() != 1)) {
    fatal_error("A hexagonal lattice without <n_axial> must have <pitch> "
                "specified by 1 number.");
  }
  pitch[0] = stod(pitch_words[0]);
  if (is_3d) {pitch[1] = stod(pitch_words[1]);}

  // Read the universes and make sure the correct number was specified.
  int n_univ = (3*n_rings*n_rings - 3*n_rings + 1) * n_axial;
  std::string univ_str {get_node_value(lat_node, "universes")};
  std::vector<std::string> univ_words {split(univ_str)};
  if (univ_words.size() != n_univ) {
    std::stringstream err_msg;
    err_msg << "Expected " << n_univ
            << " universes for a hexagonal lattice with " << n_rings
            << " rings and " << n_axial << " axial levels" << " but "
            << univ_words.size() << " were specified.";
    fatal_error(err_msg);
  }

  // Parse the universes.
  // Universes in hexagonal lattices are stored in a manner that represents
  // a skewed coordinate system: (x, alpha) rather than (x, y).  There is
  // no obvious, direct relationship between the order of universes in the
  // input and the order that they will be stored in the skewed array so
  // the following code walks a set of index values across the skewed array
  // in a manner that matches the input order.  Note that i_x = 0, i_a = 0
  // corresponds to the center of the hexagonal lattice.

  universes.resize((2*n_rings-1) * (2*n_rings-1) * n_axial, C_NONE);
  int input_index = 0;
  for (int m = 0; m < n_axial; m++) {
    // Initialize lattice indecies.
    int i_x = 1;
    int i_a = n_rings - 1;

    // Map upper triangular region of hexagonal lattice which is found in the
    // first n_rings-1 rows of the input.
    for (int k = 0; k < n_rings-1; k++) {
      // Walk the index to lower-left neighbor of last row start.
      i_x -= 1;

      // Iterate over the input columns.
      for (int j = 0; j < k+1; j++) {
        int indx = (2*n_rings-1)*(2*n_rings-1) * m
                    + (2*n_rings-1) * (i_a+n_rings-1)
                    + (i_x+n_rings-1);
        universes[indx] = std::stoi(univ_words[input_index]);
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
    for (int k = 0; k < 2*n_rings-1; k++) {
      if ((k % 2) == 0) {
        // Walk the index to the lower-left neighbor of the last row start.
        i_x -= 1;
      } else {
        // Walk the index to the lower-right neighbor of the last row start.
        i_x += 1;
        i_a -= 1;
      }

      // Iterate over the input columns.
      for (int j = 0; j < n_rings - (k % 2); j++) {
        int indx = (2*n_rings-1)*(2*n_rings-1) * m
                    + (2*n_rings-1) * (i_a+n_rings-1)
                    + (i_x+n_rings-1);
        universes[indx] = std::stoi(univ_words[input_index]);
        input_index++;
        // Walk the index to the right neighbor (which is not adjacent).
        i_x += 2;
        i_a -= 1;
      }

      // Return the lattice index to the start of the current row.
      i_x -= 2*(n_rings - (k % 2));
      i_a += n_rings - (k % 2);
    }

    // Map the lower triangular region of the hexagonal lattice.
    for (int k = 0; k < n_rings-1; k++) {
      // Walk the index to the lower-right neighbor of the last row start.
      i_x += 1;
      i_a -= 1;

      // Iterate over the input columns.
      for (int j = 0; j < n_rings-k-1; j++) {
        int indx = (2*n_rings-1)*(2*n_rings-1) * m
                    + (2*n_rings-1) * (i_a+n_rings-1)
                    + (i_x+n_rings-1);
        universes[indx] = std::stoi(univ_words[input_index]);
        input_index++;
        // Walk the index to the right neighbor (which is not adjacent).
        i_x += 2;
        i_a -= 1;
      }

      // Return lattice index to start of current row.
      i_x -= 2*(n_rings - k - 1);
      i_a += n_rings - k - 1;
    }
  }
}

//==============================================================================

int32_t&
HexLattice::operator[](const int i_xyz[3])
{
  int indx = (2*n_rings-1)*(2*n_rings-1) * i_xyz[2]
              + (2*n_rings-1) * i_xyz[1]
              + i_xyz[0];
  return universes[indx];
}

//==============================================================================

LatticeIter HexLattice::begin()
{return LatticeIter(*this, n_rings-1);}

ReverseLatticeIter HexLattice::rbegin()
{return ReverseLatticeIter(*this, universes.size()-n_rings);}

//==============================================================================

bool
HexLattice::are_valid_indices(const int i_xyz[3]) const
{
  return ((i_xyz[0] >= 0) && (i_xyz[1] >= 0) && (i_xyz[2] >= 0)
          && (i_xyz[0] < 2*n_rings-1) && (i_xyz[1] < 2*n_rings-1)
          && (i_xyz[0] + i_xyz[1] > n_rings-2)
          && (i_xyz[0] + i_xyz[1] < 3*n_rings-2)
          && (i_xyz[2] < n_axial));
}

//==============================================================================

std::pair<double, std::array<int, 3>>
HexLattice::distance(Position r, Direction u, const int i_xyz[3]) const
{
  // Compute the direction on the hexagonal basis.
  double beta_dir = u.x * std::sqrt(3.0) / 2.0 + u.y / 2.0;
  double gamma_dir = u.x * std::sqrt(3.0) / 2.0 - u.y / 2.0;

  // Note that hexagonal lattice distance calculations are performed
  // using the particle's coordinates relative to the neighbor lattice
  // cells, not relative to the particle's current cell.  This is done
  // because there is significant disagreement between neighboring cells
  // on where the lattice boundary is due to finite precision issues.

  // Upper-right and lower-left sides.
  double d {INFTY};
  std::array<int, 3> lattice_trans;
  double edge = -copysign(0.5*pitch[0], beta_dir);  // Oncoming edge
  Position r_t;
  if (beta_dir > 0) {
    const int i_xyz_t[3] {i_xyz[0]+1, i_xyz[1], i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  } else {
    const int i_xyz_t[3] {i_xyz[0]-1, i_xyz[1], i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  }
  double beta = r_t.x * std::sqrt(3.0) / 2.0 + r_t.y / 2.0;
  if ((std::abs(beta - edge) > FP_PRECISION) && beta_dir != 0) {
    d = (edge - beta) / beta_dir;
    if (beta_dir > 0) {
      lattice_trans = {1, 0, 0};
    } else {
      lattice_trans = {-1, 0, 0};
    }
  }

  // Lower-right and upper-left sides.
  edge = -copysign(0.5*pitch[0], gamma_dir);
  if (gamma_dir > 0) {
    const int i_xyz_t[3] {i_xyz[0]+1, i_xyz[1]-1, i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  } else {
    const int i_xyz_t[3] {i_xyz[0]-1, i_xyz[1]+1, i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  }
  double gamma = r_t.x * std::sqrt(3.0) / 2.0 - r_t.y / 2.0;
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

  // Upper and lower sides.
  edge = -copysign(0.5*pitch[0], u.y);
  if (u.y > 0) {
    const int i_xyz_t[3] {i_xyz[0], i_xyz[1]+1, i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  } else {
    const int i_xyz_t[3] {i_xyz[0], i_xyz[1]-1, i_xyz[2]};
    r_t = get_local_position(r, i_xyz_t);
  }
  if ((std::abs(r_t.y - edge) > FP_PRECISION) && u.y != 0) {
    double this_d = (edge - r_t.y) / u.y;
    if (this_d < d) {
      if (u.y > 0) {
        lattice_trans = {0, 1, 0};
      } else {
        lattice_trans = {0, -1, 0};
      }
      d = this_d;
    }
  }

  // Top and bottom sides
  if (is_3d) {
    double z = r.z;
    double z0 {copysign(0.5 * pitch[1], u.z)};
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
HexLattice::get_indices(Position r) const
{
  // Offset the xyz by the lattice center.
  Position r_o {r.x - center.x, r.y - center.y, r.z};
  if (is_3d) {r_o.z -= center.z;}

  // Index the z direction.
  std::array<int, 3> out;
  if (is_3d) {
    out[2] = static_cast<int>(std::ceil(r_o.z / pitch[1] + 0.5 * n_axial))-1;
  } else {
    out[2] = 0;
  }

  // Convert coordinates into skewed bases.  The (x, alpha) basis is used to
  // find the index of the global coordinates to within 4 cells.
  double alpha = r_o.y - r_o.x / std::sqrt(3.0);
  out[0] = static_cast<int>(std::floor(r_o.x
                                       / (0.5*std::sqrt(3.0) * pitch[0])));
  out[1] = static_cast<int>(std::floor(alpha / pitch[0]));

  // Add offset to indices (the center cell is (i_x, i_alpha) = (0, 0) but
  // the array is offset so that the indices never go below 0).
  out[0] += n_rings-1;
  out[1] += n_rings-1;

  // Calculate the (squared) distance between the particle and the centers of
  // the four possible cells.  Regular hexagonal tiles form a Voronoi
  // tessellation so the xyz should be in the hexagonal cell that it is closest
  // to the center of.  This method is used over a method that uses the
  // remainders of the floor divisions above because it provides better finite
  // precision performance.  Squared distances are used becasue they are more
  // computationally efficient than normal distances.
  int k {1};
  int k_min {1};
  double d_min {INFTY};
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      int i_xyz[3] {out[0] + j, out[1] + i, 0};
      Position r_t = get_local_position(r, i_xyz);
      double d = r_t.x*r_t.x + r_t.y*r_t.y;
      if (d < d_min) {
        d_min = d;
        k_min = k;
      }
      k++;
    }
  }

  // Select the minimum squared distance which corresponds to the cell the
  // coordinates are in.
  if (k_min == 2) {
    out[0] += 1;
  } else if (k_min == 3) {
    out[1] += 1;
  } else if (k_min == 4) {
    out[0] += 1;
    out[1] += 1;
  }

  return out;
}

//==============================================================================

Position
HexLattice::get_local_position(Position r, const int i_xyz[3]) const
{
  // x_l = x_g - (center + pitch_x*cos(30)*index_x)
  r.x -= (center.x + std::sqrt(3.0)/2.0 * (i_xyz[0] - n_rings + 1) * pitch[0]);
  // y_l = y_g - (center + pitch_x*index_x + pitch_y*sin(30)*index_y)
  r.y -= (center.y + (i_xyz[1] - n_rings + 1) * pitch[0]
                   + (i_xyz[0] - n_rings + 1) * pitch[0] / 2.0);
  if (is_3d) {
    r.z -= center.z - (0.5 * n_axial - i_xyz[2] - 0.5) * pitch[1];
  }

  return r;
}

//==============================================================================

bool
HexLattice::is_valid_index(int indx) const
{
  int nx {2*n_rings - 1};
  int ny {2*n_rings - 1};
  int nz {n_axial};
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
  int nx {2*n_rings - 1};
  int ny {2*n_rings - 1};
  int nz {n_axial};
  return offsets[nx*ny*nz*map + nx*ny*i_xyz[2] + nx*i_xyz[1] + i_xyz[0]];
}

//==============================================================================

std::string
HexLattice::index_to_string(int indx) const
{
  int nx {2*n_rings - 1};
  int ny {2*n_rings - 1};
  int iz {indx / (nx * ny)};
  int iy {(indx - nx * ny * iz) / nx};
  int ix {indx - nx * ny * iz - nx * iy};
  std::string out {std::to_string(ix - n_rings + 1)};
  out += ',';
  out += std::to_string(iy - n_rings + 1);
  if (is_3d) {
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
  write_dataset(lat_group, "n_rings", n_rings);
  write_dataset(lat_group, "n_axial", n_axial);
  if (is_3d) {
    write_dataset(lat_group, "pitch", pitch);
    write_dataset(lat_group, "center", center);
  } else {
    std::array<double, 1> pitch_short {{pitch[0]}};
    write_dataset(lat_group, "pitch", pitch_short);
    std::array<double, 2> center_short {{center[0], center[1]}};
    write_dataset(lat_group, "center", center_short);
  }

  // Write the universe ids.
  hsize_t nx {static_cast<hsize_t>(2*n_rings - 1)};
  hsize_t ny {static_cast<hsize_t>(2*n_rings - 1)};
  hsize_t nz {static_cast<hsize_t>(n_axial)};
  int out[nx*ny*nz];

  for (int m = 0; m < nz; m++) {
    for (int k = 0; k < ny; k++) {
      for (int j = 0; j < nx; j++) {
        int indx = nx*ny*m + nx*k + j;
        if (j + k < n_rings - 1) {
          // This array position is never used; put a -1 to indicate this.
          out[indx] = -1;
        } else if (j + k > 3*n_rings - 3) {
          // This array position is never used; put a -1 to indicate this.
          out[indx] = -1;
        } else {
          out[indx] = global_universes[universes[indx]]->id;
        }
      }
    }
  }

  hsize_t dims[3] {nz, ny, nx};
  write_int(lat_group, 3, dims, "universes", out, false);
}

//==============================================================================
// Non-method functions
//==============================================================================

extern "C" void
read_lattices(pugi::xml_node *node)
{
  for (pugi::xml_node lat_node : node->children("lattice")) {
    lattices_c.push_back(new RectLattice(lat_node));
  }
  for (pugi::xml_node lat_node : node->children("hex_lattice")) {
    lattices_c.push_back(new HexLattice(lat_node));
  }

  // Fill the lattice map.
  for (int i_lat = 0; i_lat < lattices_c.size(); i_lat++) {
    int id = lattices_c[i_lat]->id;
    auto in_map = lattice_map.find(id);
    if (in_map == lattice_map.end()) {
      lattice_map[id] = i_lat;
    } else {
      std::stringstream err_msg;
      err_msg << "Two or more lattices use the same unique ID: " << id;
      fatal_error(err_msg);
    }
  }
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Lattice* lattice_pointer(int lat_ind) {return lattices_c[lat_ind];}

  int32_t lattice_id(Lattice *lat) {return lat->id;}

  bool lattice_are_valid_indices(Lattice *lat, const int i_xyz[3])
  {return lat->are_valid_indices(i_xyz);}

  void lattice_distance(Lattice *lat, const double xyz[3], const double uvw[3],
                        const int i_xyz[3], double *d, int lattice_trans[3])
  {
    Position r {xyz};
    Direction u {uvw};
    std::pair<double, std::array<int, 3>> ld {lat->distance(r, u, i_xyz)};
    *d = ld.first;
    lattice_trans[0] = ld.second[0];
    lattice_trans[1] = ld.second[1];
    lattice_trans[2] = ld.second[2];
  }

  void lattice_get_indices(Lattice *lat, const double xyz[3], int i_xyz[3])
  {
    Position r {xyz};
    std::array<int, 3> inds = lat->get_indices(r);
    i_xyz[0] = inds[0];
    i_xyz[1] = inds[1];
    i_xyz[2] = inds[2];
  }

  void lattice_get_local_xyz(Lattice *lat, const double global_xyz[3],
                             const int i_xyz[3], double local_xyz[3])
  {
    Position global {global_xyz};
    Position local = lat->get_local_position(global, i_xyz);
    local_xyz[0] = local.x;
    local_xyz[1] = local.y;
    local_xyz[2] = local.z;
  }

  int32_t lattice_offset(Lattice *lat, int map, const int i_xyz[3])
  {return lat->offset(map, i_xyz);}

  int32_t lattice_outer(Lattice *lat) {return lat->outer;}

  void lattice_to_hdf5(Lattice *lat, hid_t group) {lat->to_hdf5(group);}

  int32_t lattice_universe(Lattice *lat, const int i_xyz[3])
  {return (*lat)[i_xyz];}
}

} // namespace openmc
