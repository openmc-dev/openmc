#ifndef OPENMC_PLOT_H
#define OPENMC_PLOT_H

#include "hdf5.h"
#include "position.h"
#include <map>
#include <sstream>

namespace openmc {

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

  class ObjectColor {
  public:
    int rgb[3];
  };

//===============================================================================
// ObjectPlot holds plot information
//===============================================================================
  
  class ObjectPlot {
  public:
    int id;
    std::string path_plot;
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
    std::vector<ObjectColor> colors;
    
  };
  
extern "C" void voxel_init(hid_t file_id, const hsize_t* dims, hid_t* dspace,
                           hid_t* dset, hid_t* memspace);
extern "C" void voxel_write_slice(int x, hid_t dspace, hid_t dset,
                                  hid_t memspace, void* buf);
extern "C" void voxel_finalize(hid_t dspace, hid_t dset, hid_t memspace);
 
} // namespace openmc
#endif // OPENMC_PLOT_H
