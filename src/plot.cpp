#include "openmc/plot.h"
#include "openmc/constants.h"
#include "openmc/settings.h"
#include "openmc/error.h"
#include "openmc/particle.h"
#include "openmc/geometry.h"

namespace openmc {

const int RED   = 1;
const int GREEN = 2;
const int BLUE  = 3;

//===============================================================================
// RUN_PLOT controls the logic for making one or many plots
//===============================================================================
  
int openmc_plot_geometry() {
  int err;

  for(auto i : n_plots) {
    ObjectPlot* pl = plots[i];
    
    std::stringstream ss;
    ss << "Processing plot " << pl->id_ << ": "
       << pl->path_plot_ << "...";
      write_message(ss.str(), 5);
      
      if (pl->type_ == PLOT_TYPE::SLICE) {
        // create 2D image
        // create_ppm(pl);
        continue;
      } else if (pl->type_ == PLOT_TYPE::VOXEL) {
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

  int width = pl->pixels_[0];
  int height = pl->pixels_[1];

  double in_pixel = (pl->width_[0])/double(width);
  double out_pixel = (pl->width_[1])/double(height);

  std::vector< std::vector< std::vector<int>>> data;

  data.resize(width);
  for (auto i : data) {
    i.resize(height);
    for (auto j : i) { j.resize(3); }
  }

  int in_i, out_i;
  double xyz[3];
  if (pl->basis_ == PLOT_BASIS::XY) {
    in_i = 0;
    out_i = 1;
    xyz[0] = pl->origin_[0] - pl->width_[0] / TWO;
    xyz[1] = pl->origin_[1] - pl->width_[1] / TWO;
    xyz[2] = pl->origin_[2];
  } else if (pl->basis_ == PLOT_BASIS::XZ) {
    in_i = 0;
    out_i = 2;
    xyz[0] = pl->origin_[0] - pl->width_[0] / TWO;
    xyz[1] = pl->origin_[1];
    xyz[2] = pl->origin_[2] - pl->width_[1] / TWO;
  } else if (pl->basis_ == PLOT_BASIS::YZ) {
    in_i = 1;
    out_i = 2;
    xyz[0] = pl->origin_[0];
    xyz[1] = pl->origin_[1] - pl->width_[0] / TWO;
    xyz[2] = pl->origin_[2] - pl->width_[1] / TWO;
  }

  double dir[3];
  Particle *p = new Particle();
  p->initialize();
  std::copy(xyz, xyz+3, p->coord[0].xyz);
  std::copy(dir, dir+3, p->coord[0].uvw);
  p->coord[0].universe = openmc_root_universe;

  // local variables
  int rgb[3];
  int id;
  for(int y = 0; y < height; y++) {
    p->coord[0].xyz[out_i] = xyz[out_i] - out_pixel*(y);
    for(int x = 0; x < width; x++) {
      p->coord[0].xyz[in_i] = xyz[in_i] + in_pixel*(x);
      // position_rgb(p, pl, rgb, id);

      std::copy(rgb, rgb+3, &(data[x][y][0]));
    }
  }

  //output_ppm(pl, data);
  
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
