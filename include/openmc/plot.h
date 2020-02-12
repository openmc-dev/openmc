#ifndef OPENMC_PLOT_H
#define OPENMC_PLOT_H

#include <unordered_map>
#include <sstream>

#include "pugixml.hpp"
#include "xtensor/xarray.hpp"

#include "hdf5.h"
#include "openmc/position.h"
#include "openmc/constants.h"
#include "openmc/cell.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/particle.h"
#include "openmc/xml_interface.h"
#include "openmc/random_lcg.h"

namespace openmc {

//===============================================================================
// Global variables
//===============================================================================

class Plot;

namespace model {

extern std::vector<Plot> plots; //!< Plot instance container
extern std::unordered_map<int, int> plot_map; //!< map of plot ids to index

extern uint64_t plotter_prn_seeds[N_STREAMS]; // Random number seeds used for plotter
extern int plotter_stream; // Stream index used by the plotter

} // namespace model

//===============================================================================
// RGBColor holds color information for plotted objects
//===============================================================================

struct RGBColor {
  //Constructors
  RGBColor() : red(0), green(0), blue(0) { };
  RGBColor(const int v[3]) : red(v[0]), green(v[1]), blue(v[2]) { };
  RGBColor(int r, int g, int b) : red(r), green(g), blue(b) { };

  RGBColor(const std::vector<int> &v) {
    if (v.size() != 3) {
      throw std::out_of_range("Incorrect vector size for RGBColor.");
    }
    red = v[0];
    green = v[1];
    blue = v[2];
  }

  bool operator ==(const RGBColor& other) {
    return red == other.red && green == other.green && blue == other.blue;
  }

  // Members
  uint8_t red, green, blue;
};

// some default colors
const RGBColor WHITE {255, 255, 255};
const RGBColor RED {255,   0,   0};


typedef xt::xtensor<RGBColor, 2> ImageData;

struct IdData {
  // Constructor
  IdData(size_t h_res, size_t v_res);

  // Methods
  void set_value(size_t y, size_t x, const Particle& p, int level);
  void set_overlap(size_t y, size_t x);

  // Members
  xt::xtensor<int32_t, 3> data_; //!< 2D array of cell & material ids
};

struct PropertyData {
  // Constructor
  PropertyData(size_t h_res, size_t v_res);

  // Methods
  void set_value(size_t y, size_t x, const Particle& p, int level);
  void set_overlap(size_t y, size_t x);

  // Members
  xt::xtensor<double, 3> data_; //!< 2D array of temperature & density data
};

enum class PlotType {
  slice = 1,
  voxel = 2
};

enum class PlotBasis {
  xy = 1,
  xz = 2,
  yz = 3
};

enum class PlotColorBy {
  cells = 0,
  mats = 1
};

//===============================================================================
// Plot class
//===============================================================================
class PlotBase {
public:
  template<class T> T get_map() const;

  // Members
public:
  Position origin_; //!< Plot origin in geometry
  Position width_; //!< Plot width in geometry
  PlotBasis basis_; //!< Plot basis (XY/XZ/YZ)
  std::array<size_t, 3> pixels_; //!< Plot size in pixels
  bool color_overlaps_; //!< Show overlapping cells?
  int level_; //!< Plot universe level
};

template<class T>
T PlotBase::get_map() const {

  size_t width = pixels_[0];
  size_t height = pixels_[1];

  // get pixel size
  double in_pixel = (width_[0])/static_cast<double>(width);
  double out_pixel = (width_[1])/static_cast<double>(height);

  // size data array
  T data(width, height);

  // setup basis indices and initial position centered on pixel
  int in_i, out_i;
  Position xyz = origin_;
  switch(basis_) {
  case PlotBasis::xy :
    in_i = 0;
    out_i = 1;
    break;
  case PlotBasis::xz :
    in_i = 0;
    out_i = 2;
    break;
  case PlotBasis::yz :
    in_i = 1;
    out_i = 2;
    break;
  default:
    UNREACHABLE();
  }

  // set initial position
  xyz[in_i] = origin_[in_i] - width_[0] / 2. + in_pixel / 2.;
  xyz[out_i] = origin_[out_i] + width_[1] / 2. - out_pixel / 2.;

  // arbitrary direction
  Direction dir = {0.7071, 0.7071, 0.0};

  #pragma omp parallel
  {
    Particle p;
    p.r() = xyz;
    p.u() = dir;
    p.coord_[0].universe = model::root_universe;
    int level = level_;
    int j{};

    #pragma omp for
    for (int y = 0; y < height; y++) {
      p.r()[out_i] =  xyz[out_i] - out_pixel * y;
      for (int x = 0; x < width; x++) {
        p.r()[in_i] = xyz[in_i] + in_pixel * x;
        p.n_coord_ = 1;
        // local variables
        bool found_cell = find_cell(&p, 0);
        j = p.n_coord_ - 1;
        if (level >=0) {j = level + 1;}
        if (found_cell) {
          data.set_value(y, x, p, j);
        }
        if (color_overlaps_ && check_cell_overlap(&p, false)) {
          data.set_overlap(y, x);
        }
      } // inner for
    } // outer for
  } // omp parallel

  return data;
}

class Plot : public PlotBase {

public:
  // Constructor
  Plot(pugi::xml_node plot);

  // Methods
private:
  void set_id(pugi::xml_node plot_node);
  void set_type(pugi::xml_node plot_node);
  void set_output_path(pugi::xml_node plot_node);
  void set_bg_color(pugi::xml_node plot_node);
  void set_basis(pugi::xml_node plot_node);
  void set_origin(pugi::xml_node plot_node);
  void set_width(pugi::xml_node plot_node);
  void set_universe(pugi::xml_node plot_node);
  void set_default_colors(pugi::xml_node plot_node);
  void set_user_colors(pugi::xml_node plot_node);
  void set_meshlines(pugi::xml_node plot_node);
  void set_mask(pugi::xml_node plot_node);
  void set_overlap_color(pugi::xml_node plot_node);

// Members
public:
  int id_; //!< Plot ID
  PlotType type_; //!< Plot type (Slice/Voxel)
  PlotColorBy color_by_; //!< Plot coloring (cell/material)
  int meshlines_width_; //!< Width of lines added to the plot
  int index_meshlines_mesh_ {-1}; //!< Index of the mesh to draw on the plot
  RGBColor meshlines_color_; //!< Color of meshlines on the plot
  RGBColor not_found_ {WHITE}; //!< Plot background color
  RGBColor overlap_color_ {RED}; //!< Plot overlap color
  std::vector<RGBColor> colors_; //!< Plot colors
  std::string path_plot_; //!< Plot output filename
};

//===============================================================================
// Non-member functions
//===============================================================================

//! Add mesh lines to image data of a plot object
//! \param[in] plot object
//! \param[out] image data associated with the plot object
void draw_mesh_lines(Plot pl, ImageData& data);

//! Write a ppm image to file using a plot object's image data
//! \param[in] plot object
//! \param[out] image data associated with the plot object
void output_ppm(Plot pl, const ImageData& data);

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

//! Read plot specifications from a plots.xml file
void read_plots_xml();

//! Create a ppm image for a plot object
//! \param[in] plot object
void create_ppm(Plot pl);

//! Create an hdf5 voxel file for a plot object
//! \param[in] plot object
void create_voxel(Plot pl);

//! Create a randomly generated RGB color
//! \return RGBColor with random value
RGBColor random_color();


} // namespace openmc
#endif // OPENMC_PLOT_H
