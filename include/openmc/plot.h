#ifndef OPENMC_PLOT_H
#define OPENMC_PLOT_H

#include "hdf5.h"
#include "position.h"
#include <map>

namespace openmc {

  int PLOT_LEVEL_LOWEST = -1;

  std::map<int, int> plot_dict;
  
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
    int rgb_[3];
  };

//===============================================================================
// ObjectPlot holds plot information
//===============================================================================
  
  class ObjectPlot {

    int id_;
    std::string path_plot_;
    int type_;
    int color_by_;
    Position origin_;
    Position width_;
    int basis_;
    int pixels_[3];
    int meshlines_width_;
    int level_;
    int index_meshlines_mesh_;
    ObjectColor meshlines_color_;
    ObjectColor not_found_;
    std::vector<ObjectColor> colors_;
    
  };
  
extern "C" void voxel_init(hid_t file_id, const hsize_t* dims, hid_t* dspace,
                           hid_t* dset, hid_t* memspace);
extern "C" void voxel_write_slice(int x, hid_t dspace, hid_t dset,
                                  hid_t memspace, void* buf);
extern "C" void voxel_finalize(hid_t dspace, hid_t dset, hid_t memspace);
 
} // namespace openmc
#endif // OPENMC_PLOT_H
