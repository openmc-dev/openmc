#ifndef OPENMC_PLOT_H
#define OPENMC_PLOT_H

#include <map>
#include <sstream>

#include "hdf5.h"
#include "position.h"
#include "openmc/constants.h"
#include "openmc/particle.h"
#include "openmc/xml_interface.h"

namespace openmc {

//===============================================================================
// Global variables
//===============================================================================

  extern int PLOT_LEVEL_LOWEST; //!< lower bound on plot universe level

  extern std::map<int, int> plot_dict; //!< map of plot ids to index

  extern int n_plots; //!< number of plots in openmc run

  class ObjectPlot;
  extern std::vector<ObjectPlot*> plots; //!< Plot instance container

  typedef std::vector< std::vector< std::vector<int> > > ImageData;

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
    int rgb[3]; //!< RGB color values
  };

//===============================================================================
// ObjectPlot holds plot information
//===============================================================================

class ObjectPlot
{

public:
  // Constructor
  ObjectPlot(pugi::xml_node plot);

  // Methods
private:
  void set_id();
  void set_type();
  void set_output_path();
  void set_bg_color();
  void set_basis();
  void set_origin();
  void set_width();
  void set_universe();
  void set_default_colors();
  void set_user_colors();
  void set_meshlines();
  void set_mask();


  // Members
public:
  int id; //!< Plot ID
  int type; //!< Plot type (Slice/Voxel)
  int color_by; //!< Plot coloring (cell/material)
  Position origin; //!< Plot origin in geometry
  Position width; //!< Plot width in geometry
  int basis; //!< Plot basis (XY/XZ/YZ)
  int pixels[3]; //!< Plot size in pixels
  int meshlines_width; //!< Width of lines added to the plot
  int level; //!< Plot universe level
  int index_meshlines_mesh; //!< Index of the mesh to draw on the plot
  ObjectColor meshlines_color; //!< Color of meshlines on the plot
  ObjectColor not_found; //!< Plot background color
  std::vector<ObjectColor*> colors; //!< Plot colors
  std::string path_plot; //!< Plot output filename
  pugi::xml_node _plot_node;
};

//===============================================================================
// Non-member functions
//===============================================================================

//! Add mesh lines to image data of a plot object
//! \param[in] plot object
//! \param[out] image data associated with the plot object
void draw_mesh_lines(ObjectPlot* pl,
                     std::vector< std::vector< std::vector<int> > > &data);

//! Write a ppm image to file using a plot object's image data
//! \param[in] plot object
//! \param[out] image data associated with the plot object
void output_ppm(ObjectPlot* pl,
                const std::vector< std::vector< std::vector<int> > > &data);

//! Get the rgb color for a given particle position in a plot
//! \param[in] particle with position for current pixel
//! \param[in] plot object
//! \param[out] rgb color
//! \param[out] cell or material id for particle position
void position_rgb(Particle* p, ObjectPlot* pl, int rgb[3], int &id);


//! Initialize a voxel file
//! \param[in] id of an open hdf5 file
//! \param[in] dimensions of the voxel file (dx, dy, dz)
//! \param[out] dataspace pointer to voxel data
//! \param[out] dataset pointer to voxesl data
//! \param[out] pointer to memory space of voxel data
void voxel_init(hid_t file_id, const hsize_t* dims, hid_t* dspace,
                           hid_t* dset, hid_t* memspace);

//! Write a section of the voxel data to hdf5
//! \param[in] voxel slice
//! \param[out] dataspace pointer to voxel data
//! \param[out] dataset pointer to voxesl data
//! \param[out] pointer to data to write
void voxel_write_slice(int x, hid_t dspace, hid_t dset,
                                  hid_t memspace, void* buf);
//! Close voxel file entities
//! \param[in] data space to close
//! \param[in] dataset to close
//! \param[in] memory space to close
void voxel_finalize(hid_t dspace, hid_t dset, hid_t memspace);

//===============================================================================
// External functions
//===============================================================================

//! Read plots from plots.xml node
//! \param[in] plot node of plots.xml
extern "C" void read_plots(pugi::xml_node* plot_node);

//! Create a ppm image for a plot object
//! \param[in] plot object
extern "C" void create_ppm(ObjectPlot* pl);

//! Create an hdf5 voxel file for a plot object
//! \param[in] plot object
extern "C" void create_voxel(ObjectPlot *pl);


} // namespace openmc
#endif // OPENMC_PLOT_H
