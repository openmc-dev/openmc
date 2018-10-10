#ifndef OPENMC_PLOT_H
#define OPENMC_PLOT_H

#include <map>
#include <sstream>

#include "hdf5.h"
#include "position.h"
#include "openmc/constants.h"
#include "openmc/particle.h"

namespace openmc {

  typedef std::vector< std::vector< std::vector<int> > > ImageData;

  class ObjectPlot;

  int PLOT_LEVEL_LOWEST = -1;

  std::map<int, int> plot_dict;

  std::vector<int32_t> n_plots;

  std::vector<ObjectPlot*> plots;

  enum PLOT_TYPE {
    SLICE = 1,
    VOXEL = 2
  };

  enum PLOT_BASIS {
    XY = 1,
    XZ = 2,
    YZ = 3
  };

  enum PLOT_COLOR_BY {
    CELLS = 1,
    MATS = 2
  };

//===============================================================================
// ObjectColor holds color information for plotted objects
//===============================================================================

  struct ObjectColor {
    int rgb[3];
  };

//===============================================================================
// ObjectPlot holds plot information
//===============================================================================

  struct ObjectPlot {
    int id;
    int type;
    int color_by;
    Position origin;
    Position width;
    int basis;
    int pixels[3];
    int meshlines_width;
    int level;
    int index_meshlines_mesh;
    ObjectColor meshlines_color;
    ObjectColor not_found;
    ObjectColor colors[MAX_COORD];
    char path_plot[MAX_WORD_LEN];

  };


extern "C" void create_ppm(ObjectPlot* pl);

extern "C" void create_voxel(ObjectPlot *pl);

void draw_mesh_lines(ObjectPlot* pl,
                     std::vector< std::vector< std::vector<int> > > &data);

void output_ppm(ObjectPlot* pl,
                const std::vector< std::vector< std::vector<int> > > &data);

void position_rgb(Particle* p, ObjectPlot* pl, int rgb[3], int &id);

void voxel_init(hid_t file_id, const hsize_t* dims, hid_t* dspace,
                           hid_t* dset, hid_t* memspace);
void voxel_write_slice(int x, hid_t dspace, hid_t dset,
                                  hid_t memspace, void* buf);
void voxel_finalize(hid_t dspace, hid_t dset, hid_t memspace);

} // namespace openmc
#endif // OPENMC_PLOT_H
