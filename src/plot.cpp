#include <fstream>

#include "openmc/plot.h"
#include "openmc/constants.h"
#include "openmc/settings.h"
#include "openmc/error.h"
#include "openmc/particle.h"
#include "openmc/geometry.h"
#include "openmc/cell.h"
#include "openmc/material.h"
#include "openmc/string_functions.h"

namespace openmc {

const int RED   = 1;
const int GREEN = 2;
const int BLUE  = 3;

const int WHITE[3] = {255, 255, 255};
const int NULLRGB[3] = {0, 0, 0};

//===============================================================================
// RUN_PLOT controls the logic for making one or many plots
//===============================================================================

int openmc_plot_geometry() {
  int err;

  for (auto i : n_plots) {
    ObjectPlot* pl = plots[i];

    std::stringstream ss;
    ss << "Processing plot " << pl->id << ": "
       << pl->path_plot << "...";
      write_message(ss.str(), 5);

      if (PLOT_TYPE::SLICE == pl->type) {
        // create 2D image
        // create_ppm(pl);
        continue;
      } else if (PLOT_TYPE::VOXEL == pl->type) {
        // create voxel file for 3D viewing
        // create_voxel(pl);
        continue;
      }

  }

  return 0;
}

//===============================================================================
// CREATE_PPM creates an image based on user input from a plots.xml <plot>
// specification in the portable pixmap format (PPM)
//===============================================================================

void create_ppm(ObjectPlot* pl) {

  int width = pl->pixels[0];
  int height = pl->pixels[1];

  double in_pixel = (pl->width[0])/double(width);
  double out_pixel = (pl->width[1])/double(height);

  std::vector< std::vector< std::vector<int> > > data;

  data.resize(width);
  for (auto & i : data) {
    i.resize(height);
    for (auto & j : i) {
      j.resize(3);
    }
  }

  int in_i, out_i;
  double xyz[3];
  if (PLOT_BASIS::XY == pl->basis) {
    in_i = 0;
    out_i = 1;
    xyz[0] = pl->origin[0] - pl->width[0] / TWO;
    xyz[1] = pl->origin[1] + pl->width[1] / TWO;
    xyz[2] = pl->origin[2];
  } else if (PLOT_BASIS::XZ == pl->basis) {
    in_i = 0;
    out_i = 2;
    xyz[0] = pl->origin[0] - pl->width[0] / TWO;
    xyz[1] = pl->origin[1];
    xyz[2] = pl->origin[2] + pl->width[1] / TWO;
  } else if (PLOT_BASIS::YZ == pl->basis) {
    in_i = 1;
    out_i = 2;
    xyz[0] = pl->origin[0];
    xyz[1] = pl->origin[1] - pl->width[0] / TWO;
    xyz[2] = pl->origin[2] + pl->width[1] / TWO;
  }

  double dir[3] = {HALF, HALF, HALF};
  Particle *p = new Particle();
  p->initialize();
  std::copy(xyz, xyz+3, p->coord[0].xyz);
  std::copy(dir, dir+3, p->coord[0].uvw);
  p->coord[0].universe = openmc_root_universe;

  // local variables
  int rgb[3];
  int id;
  for (int y = 0; y < height; y++) {
    p->coord[0].xyz[out_i] = xyz[out_i] - out_pixel*(y);
    for (int x = 0; x < width; x++) {
      p->coord[0].xyz[in_i] = xyz[in_i] + in_pixel*(x);
      position_rgb(p, pl, rgb, id);
      data[x][y][0] = rgb[0];
      data[x][y][1] = rgb[1];
      data[x][y][2] = rgb[2];
    }
  }

  output_ppm(pl, data);

}

//===============================================================================
// POSITION_RGB computes the red/green/blue values for a given plot with the
// current particle's position
//===============================================================================

void position_rgb(Particle* p, ObjectPlot* pl, int rgb[3], int &id) {
  bool found_cell;

  p->n_coord = 1;

  found_cell = find_cell(p, 0);

  int j = p->n_coord - 1;

  if (settings::check_overlaps) { check_cell_overlap(p); }

  // Set coordinate level if specified
  if (pl->level >= 0) {j = pl->level + 1;}

  Cell* c;

  if (!found_cell) {
    // If no cell, revert to default color
    std::copy(pl->not_found.rgb,
              pl->not_found.rgb + 3,
              rgb);
    id = -1;
  } else {
    if (PLOT_COLOR_BY::MATS == pl->color_by) {
      // Assign color based on material
      c = cells[p->coord[j].cell];
      if (c->type_ == FILL_UNIVERSE) {
        // If we stopped on a middle universe level, treat as if not found
        std::copy(pl->not_found.rgb,
                  pl->not_found.rgb + 3,
                  rgb);
        id = -1;
      } else if (p->material == MATERIAL_VOID) {
        // By default, color void cells white
        std::copy(WHITE, WHITE+3, rgb);
        id = -1;
      } else {
        std::copy(pl->colors[p->material - 1].rgb,
                  pl->colors[p->material - 1].rgb + 3,
                  rgb);
        id = materials[p->material - 1]->id;
      }

    } else if (PLOT_COLOR_BY::CELLS == pl->color_by) {
      // Assign color based on cell
      std::copy(pl->colors[p->coord[j].cell].rgb,
                pl->colors[p->coord[j].cell].rgb + 3,
                rgb);
      id = cells[p->coord[j].cell]->id_;
    } else {
      std::copy(NULLRGB, NULLRGB+3, rgb);
      id = -1;
    }

  } // endif found_cell
}

//===============================================================================
// OUTPUT_PPM writes out a previously generated image to a PPM file
//===============================================================================

  void output_ppm(ObjectPlot* pl,
                  std::vector< std::vector< std::vector<int> > > data)
{

  // Open PPM file for writing
  std::string fname = std::string(pl->path_plot);
  fname = strtrim(fname);
  std::ofstream of;

  of.open(fname);
  
  // Write header
  of << "P6" << std::endl;
  of << pl->pixels[0] << " " << pl->pixels[1] << std::endl;
  of << "255" << std::endl;
  of.close();

  of.open(fname, std::ios::binary | std::ios::app);
  // Write color for each pixel
  for (int y = 0; y < pl->pixels[1]; y++) {
    for (int x = 0; x < pl->pixels[0]; x++) {
      std::vector<int> rgb = data[x][y];
      of.write((char*)&rgb[0], 1);
      of.write((char*)&rgb[1], 1);
      of.write((char*)&rgb[2], 1);
    }
  }
  // Close file
  of.close();
}
  
void
voxel_init(hid_t file_id, const hsize_t* dims, hid_t* dspace, hid_t* dset,
           hid_t* memspace)
{
  // Create dataspace/dataset for voxel data
  *dspace = H5Screate_simple(3, dims, nullptr);
  *dset = H5Dcreate(file_id, "data", H5T_NATIVE_INT, *dspace, H5P_DEFAULT,
                    H5P_DEFAULT, H5P_DEFAULT);

  // Create dataspace for a slice of the voxel
  hsize_t dims_slice[2] {dims[1], dims[2]};
  *memspace = H5Screate_simple(2, dims_slice, nullptr);

  // Select hyperslab in dataspace
  hsize_t start[3] {0, 0, 0};
  hsize_t count[3] {1, dims[1], dims[2]};
  H5Sselect_hyperslab(*dspace, H5S_SELECT_SET, start, nullptr, count, nullptr);
}


void
voxel_write_slice(int x, hid_t dspace, hid_t dset, hid_t memspace, void* buf)
{
  hssize_t offset[3] {x - 1, 0, 0};
  H5Soffset_simple(dspace, offset);
  H5Dwrite(dset, H5T_NATIVE_INT, memspace, dspace, H5P_DEFAULT, buf);
}


void
voxel_finalize(hid_t dspace, hid_t dset, hid_t memspace)
{
  H5Dclose(dset);
  H5Sclose(dspace);
  H5Sclose(memspace);
}

} // namespace openmc
