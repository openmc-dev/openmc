#ifndef OPENMC_PLOT_H
#define OPENMC_PLOT_H

#include <unordered_map>
#include <sstream>

#include "pugixml.hpp"
#include "xtensor/xarray.hpp"

#include "hdf5.h"
#include "openmc/position.h"
#include "openmc/constants.h"
#include "openmc/particle.h"
#include "openmc/xml_interface.h"

namespace openmc {

//===============================================================================
// Global variables
//===============================================================================

class Plot;

namespace model {

extern std::vector<Plot> plots; //!< Plot instance container
extern std::unordered_map<int, int> plot_map; //!< map of plot ids to index

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

  // Members
  uint8_t red, green, blue;
};

typedef xt::xtensor<RGBColor, 2> ImageData;
typedef xt::xtensor<int32_t, 3>  IdData;
typedef xt::xtensor<double, 3>  PropertyData;

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
  cells = 1,
  mats = 2
};

//===============================================================================
// Plot class
//===============================================================================
class PlotBase {
 public:
  IdData get_id_map() const;
  PropertyData get_property_map() const;

 private:
  template<class D, typename setter>
    D generate_data() const;

  // Members
 public:
  Position origin_; //!< Plot origin in geometry
  Position width_; //!< Plot width in geometry
  PlotBasis basis_; //!< Plot basis (XY/XZ/YZ)
  std::array<int, 3> pixels_; //!< Plot size in pixels
  int level_; //!< Plot universe level
};

class Plot : public PlotBase
{

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

// Members
public:
  int id_; //!< Plot ID
  PlotType type_; //!< Plot type (Slice/Voxel)
  PlotColorBy color_by_; //!< Plot coloring (cell/material)
  int meshlines_width_; //!< Width of lines added to the plot
  int index_meshlines_mesh_; //!< Index of the mesh to draw on the plot
  RGBColor meshlines_color_; //!< Color of meshlines on the plot
  RGBColor not_found_; //!< Plot background color
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

//! Get the rgb color for a given particle position in a plot
//! \param[in] particle with position for current pixel
//! \param[in] plot object
//! \param[out] rgb color
//! \param[out] cell or material id for particle position
void position_rgb(Particle p, Plot pl, RGBColor& rgb, int& id);

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
