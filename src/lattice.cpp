#include "lattice.h"

#include <sstream>
#include <vector>

#include "error.h"
#include "hdf5_interface.h"
#include "xml_interface.h"

//TODO: this is only inlcuded for constants that should be moved elsewhere
//#include "surface.h"

//TODO: remove this include
#include <iostream>


namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

// Braces force n_lattices to be defined here, not just declared.
//extern "C" {int32_t n_lattices {0};}

//Lattice **lattices_c;
std::vector<Lattice*> lattices_c;

std::map<int32_t, int32_t> lattice_dict;

//==============================================================================
// Lattice implementation
//==============================================================================

Lattice::Lattice(pugi::xml_node lat_node)
{
  if (check_for_node(lat_node, "id")) {
    id = stoi(get_node_value(lat_node, "id"));
  } else {
    fatal_error("Must specify id of lattice in geometry XML file.");
  }

  if (check_for_node(lat_node, "name")) {
    name = get_node_value(lat_node, "name");
  }

  if (check_for_node(lat_node, "outer")) {
    outer = stoi(get_node_value(lat_node, "outer"));
  }
}

void
Lattice::to_hdf5(hid_t lat_group) const
{
  if (!name.empty()) {
    write_string(lat_group, "name", name);
  }
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
    n_cells[0] = stoi(dimension_words[0]);
    n_cells[1] = stoi(dimension_words[1]);
    n_cells[2] = 1;
    is_3d = false;
  } else if (dimension_words.size() == 3) {
    n_cells[0] = stoi(dimension_words[0]);
    n_cells[1] = stoi(dimension_words[1]);
    n_cells[2] = stoi(dimension_words[2]);
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
  int nx = n_cells[0];
  int ny = n_cells[1];
  int nz = n_cells[2];
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
  universes.resize(nx*ny*nz, -1);
  for (int iz = 0; iz < nz; iz++) {
    for (int iy = ny-1; iy > -1; iy--) {
      for (int ix = 0; ix < nx; ix++) {
        int indx1 = nx*ny*iz + nx*(ny-iy-1) + ix;
        int indx2 = nx*ny*iz + nx*iy + ix;
        universes[indx1] = stoi(univ_words[indx2]);
      }
    }
  }
}

std::pair<double, std::array<int, 3>>
RectLattice::distance(const double xyz[3], const double uvw[3]) const
{
  // Get short aliases to the coordinates.
  double x {xyz[0]};
  double y {xyz[1]};
  double z {xyz[2]};
  double u {uvw[0]};
  double v {uvw[1]};

  // Determine the oncoming edge.
  double x0 {copysign(0.5 * pitch[0], u)};
  double y0 {copysign(0.5 * pitch[1], v)};

  // Left and right sides
  double d {INFTY};
  std::array<int, 3> lattice_trans;
  if ((std::abs(x - x0) > FP_PRECISION) && u != 0) {
    d = (x0 - x) / u;
    if (u > 0) {
      lattice_trans = {1, 0, 0};
    } else {
      lattice_trans = {-1, 0, 0};
    }
  }

  // Front and back sides
  if ((std::abs(y - y0) > FP_PRECISION) && v != 0) {
    double this_d = (y0 - y) / v;
    if (this_d < d) {
      d = this_d;
      if (v > 0) {
        lattice_trans = {0, 1, 0};
      } else {
        lattice_trans = {0, -1, 0};
      }
    }
  }

  // Top and bottom sides
  if (is_3d) {
    double w {uvw[2]};
    double z0 {copysign(0.5 * pitch[2], w)};
    if ((std::abs(z - z0) > FP_PRECISION) && w != 0) {
      double this_d = (z0 - z) / w;
      if (this_d < d) {
        d = this_d;
        if (w > 0) {
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
// HexLattice implementation
//==============================================================================

HexLattice::HexLattice(pugi::xml_node lat_node)
  : Lattice {lat_node}
{
  // Read the number of lattice cells in each dimension.
  n_rings = stoi(get_node_value(lat_node, "n_rings"));
  if (check_for_node(lat_node, "n_axial")) {
    n_axial = stoi(get_node_value(lat_node, "n_axial"));
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
    fatal_error("A hexagonal lattice without <n_axial> must have <center> "
                "specified by 1 number.");
  }
  pitch[0] = stod(pitch_words[0]);
  if (is_3d) {pitch[1] = stod(pitch_words[1]);}

  // Read the universes and make sure the correct number was specified.
  //int n_univ = (2*n_rings - 1) * (2*n_rings - 1) * n_axial;
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

  universes.resize((2*n_rings-1) * (2*n_rings-1) * n_axial, -1);
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
        universes[indx] = stoi(univ_words[input_index]);
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
        universes[indx] = stoi(univ_words[input_index]);
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
        universes[indx] = stoi(univ_words[input_index]);
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

std::pair<double, std::array<int, 3>>
HexLattice::distance(const double xyz[3], const double uvw[3]) const
{
}

//==============================================================================

extern "C" void
read_lattices(pugi::xml_node *node)
{
  for (pugi::xml_node lat_node: node->children("lattice")) {
    lattices_c.push_back(new RectLattice(lat_node));
  }
  for (pugi::xml_node lat_node: node->children("hex_lattice")) {
    lattices_c.push_back(new HexLattice(lat_node));
  }

  // Fill the lattice dictionary.
  for (int i_lat = 0; i_lat < lattices_c.size(); i_lat++) {
    int id = lattices_c[i_lat]->id;
    auto in_dict = lattice_dict.find(id);
    if (in_dict == lattice_dict.end()) {
      lattice_dict[id] = i_lat;
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

  void lattice_distance(Lattice *lat, const double xyz[3], const double uvw[3],
                        double *d, int lattice_trans[3])
  {
    std::pair<double, std::array<int, 3>> ld {lat->distance(xyz, uvw)};
    *d = ld.first;
    lattice_trans[0] = ld.second[0];
    lattice_trans[1] = ld.second[1];
    lattice_trans[2] = ld.second[2];
  }

  void lattice_to_hdf5(Lattice *lat, hid_t group) {lat->to_hdf5(group);}
}

} // namespace openmc
