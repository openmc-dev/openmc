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
#include "openmc/mesh.h"

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

  if (pl->index_meshlines_mesh >= 0) { draw_mesh_lines(pl, data); }
  
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
                  const std::vector< std::vector< std::vector<int> > > &data)
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
  // THIS IS HERE TO MATCH FORTRAN VERSION, NOT NECESSARY  
  of << std::endl; 
  of.close();
}

//===============================================================================
// DRAW_MESH_LINES draws mesh line boundaries on an image
//===============================================================================
  
void draw_mesh_lines(ObjectPlot *pl,
                     std::vector< std::vector< std::vector<int> > > &data) {

  std::vector<int> rgb; rgb.resize(3);
  rgb[0] = pl->meshlines_color.rgb[0];
  rgb[1] = pl->meshlines_color.rgb[1];
  rgb[2] = pl->meshlines_color.rgb[2];

  int outer, inner;
  switch(pl->basis){
  case PLOT_BASIS::XY :
    outer = 0;
    inner = 1;
    break;
  case PLOT_BASIS::XZ :
    outer = 0;
    inner = 2;
    break;
  case PLOT_BASIS::YZ :
    outer = 1;
    inner = 2;
    break;
  }

  double xyz_ll_plot[3], xyz_ur_plot[3];
  std::copy((double*)&pl->origin, (double*)&pl->origin + 3, xyz_ll_plot);
  std::copy((double*)&pl->origin, (double*)&pl->origin + 3, xyz_ur_plot);

  xyz_ll_plot[outer] = pl->origin[outer] - pl->width[0] / TWO;
  xyz_ll_plot[inner] = pl->origin[inner] - pl->width[1] / TWO;
  xyz_ur_plot[outer] = pl->origin[outer] + pl->width[0] / TWO;
  xyz_ur_plot[inner] = pl->origin[inner] + pl->width[1] / TWO;

  int width[3];
  width[0] = xyz_ur_plot[0] - xyz_ll_plot[0];
  width[1] = xyz_ur_plot[1] - xyz_ll_plot[1];
  width[2] = xyz_ur_plot[2] - xyz_ll_plot[2];

  auto &m = meshes[pl->index_meshlines_mesh];

  int ijk_ll[3], ijk_ur[3];
  bool in_mesh;
  m->get_indices(Position(xyz_ll_plot), &(ijk_ll[0]), &in_mesh);
  m->get_indices(Position(xyz_ur_plot), &(ijk_ur[0]), &in_mesh);

  // Fortran/C++ index correction
  ijk_ur[0]++; ijk_ur[1]++; ijk_ur[2]++;
  
  double frac;
  int outrange[3], inrange[3];
  double xyz_ll[3], xyz_ur[3];
  // sweep through all meshbins on this plane and draw borders  
  for (int i = ijk_ll[outer]; i <= ijk_ur[outer]; i++) {
    for (int j = ijk_ll[inner]; j <= ijk_ur[inner]; j++) {
      // check if we're in the mesh for this ijk
      if (i > 0 && i <= m->shape_[outer] && j >0 && j <= m->shape_[inner] ) {
      
        // get xyz's of lower left and upper right of this mesh cell
        xyz_ll[outer] = m->lower_left_[outer] + m->width_[outer] * (i - 1);
        xyz_ll[inner] = m->lower_left_[inner] + m->width_[inner] * (j - 1);
        xyz_ur[outer] = m->lower_left_[outer] + m->width_[outer] * i;
        xyz_ur[inner] = m->lower_left_[inner] + m->width_[inner] * j;

        // map the xyz ranges to pixel ranges
        frac = (xyz_ll[outer] - xyz_ll_plot[outer]) / width[outer];
        outrange[0] = int(frac * double(pl->pixels[0]));
        frac = (xyz_ur[outer] - xyz_ll_plot[outer]) / width[outer];
        outrange[1] = int(frac * double(pl->pixels[0]));

        frac = (xyz_ur[inner] - xyz_ll_plot[inner]) / width[inner];
        inrange[0] = int((ONE - frac) * (double)pl->pixels[1]);
        frac = (xyz_ll[inner] - xyz_ll_plot[inner]) / width[inner];
        inrange[1] = int((ONE - frac) * (double)pl->pixels[1]);

        // draw lines
        for (int out_ = outrange[0]; out_ <= outrange[1]; out_++) {
          for (int plus = 0; plus <= pl->meshlines_width; plus++) {
            data[out_][inrange[0] + plus] = rgb;
            data[out_][inrange[1] + plus] = rgb;
            data[out_][inrange[0] - plus] = rgb;
            data[out_][inrange[1] - plus] = rgb;
          }
        }

        for (int in_ = inrange[0]; in_ <= inrange[1]; in_++) {
          for (int plus = 0; plus <= pl->meshlines_width; plus++) {
            data[outrange[0] + plus][in_] = rgb;
            data[outrange[1] + plus][in_] = rgb;
            data[outrange[0] - plus][in_] = rgb;
            data[outrange[1] - plus][in_] = rgb;
          }
        }

      } // end if(in mesh)
    }
  } // end outer loops
        
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
