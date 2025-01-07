#ifndef OPENMC_PLOT_H
#define OPENMC_PLOT_H

#include <cmath>
#include <set>
#include <sstream>
#include <unordered_map>

#include "pugixml.hpp"
#include "xtensor/xarray.hpp"

#include "hdf5.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"
#include "openmc/xml_interface.h"

namespace openmc {

//===============================================================================
// Global variables
//===============================================================================

class PlottableInterface;

namespace model {

extern std::unordered_map<int, int> plot_map; //!< map of plot ids to index
extern vector<std::unique_ptr<PlottableInterface>>
  plots; //!< Plot instance container

extern uint64_t plotter_seed; // Stream index used by the plotter

} // namespace model

//===============================================================================
// RGBColor holds color information for plotted objects
//===============================================================================

struct RGBColor {
  // Constructors
  RGBColor() : red(0), green(0), blue(0) {};
  RGBColor(const int v[3]) : red(v[0]), green(v[1]), blue(v[2]) {};
  RGBColor(int r, int g, int b) : red(r), green(g), blue(b) {};

  RGBColor(const vector<int>& v)
  {
    if (v.size() != 3) {
      throw std::out_of_range("Incorrect vector size for RGBColor.");
    }
    red = v[0];
    green = v[1];
    blue = v[2];
  }

  bool operator==(const RGBColor& other)
  {
    return red == other.red && green == other.green && blue == other.blue;
  }

  RGBColor& operator*=(const double x)
  {
    red *= x;
    green *= x;
    blue *= x;
    return *this;
  }

  // Members
  uint8_t red, green, blue;
};

// some default colors
const RGBColor WHITE {255, 255, 255};
const RGBColor RED {255, 0, 0};
const RGBColor BLACK {0, 0, 0};

/**
 * @class PlottableInterface
 * @brief Interface for plottable objects.
 *
 * PlottableInterface classes must have a unique ID in the plots.xml file.
 * They guarantee the ability to create output in some form. This interface
 * is designed to be implemented by classes that produce plot-relevant data
 * which can be visualized.
 */
class PlottableInterface {
private:
  void set_id(pugi::xml_node plot_node);
  int id_; // unique plot ID

  void set_bg_color(pugi::xml_node plot_node);
  void set_universe(pugi::xml_node plot_node);
  void set_default_colors(pugi::xml_node plot_node);
  void set_user_colors(pugi::xml_node plot_node);
  void set_overlap_color(pugi::xml_node plot_node);
  void set_mask(pugi::xml_node plot_node);

protected:
  // Plot output filename, derived classes have logic to set it
  std::string path_plot_;

public:
  enum class PlotColorBy { cells = 0, mats = 1 };

  // Creates the output image named path_plot_
  virtual void create_output() const = 0;

  // Print useful info to the terminal
  virtual void print_info() const = 0;

  const std::string& path_plot() const { return path_plot_; }
  std::string& path_plot() { return path_plot_; }
  int id() const { return id_; }
  int level() const { return level_; }

  // Public color-related data
  PlottableInterface(pugi::xml_node plot_node);
  virtual ~PlottableInterface() = default;
  int level_;                    // Universe level to plot
  bool color_overlaps_;          // Show overlapping cells?
  PlotColorBy color_by_;         // Plot coloring (cell/material)
  RGBColor not_found_ {WHITE};   // Plot background color
  RGBColor overlap_color_ {RED}; // Plot overlap color
  vector<RGBColor> colors_;      // Plot colors
};

typedef xt::xtensor<RGBColor, 2> ImageData;

struct IdData {
  // Constructor
  IdData(size_t h_res, size_t v_res);

  // Methods
  void set_value(size_t y, size_t x, const GeometryState& p, int level);
  void set_overlap(size_t y, size_t x);

  // Members
  xt::xtensor<int32_t, 3> data_; //!< 2D array of cell & material ids
};

struct PropertyData {
  // Constructor
  PropertyData(size_t h_res, size_t v_res);

  // Methods
  void set_value(size_t y, size_t x, const GeometryState& p, int level);
  void set_overlap(size_t y, size_t x);

  // Members
  xt::xtensor<double, 3> data_; //!< 2D array of temperature & density data
};

//===============================================================================
// Plot class
//===============================================================================

class SlicePlotBase {
public:
  template<class T>
  T get_map() const;

  enum class PlotBasis { xy = 1, xz = 2, yz = 3 };

  // Members
public:
  Position origin_;           //!< Plot origin in geometry
  Position width_;            //!< Plot width in geometry
  PlotBasis basis_;           //!< Plot basis (XY/XZ/YZ)
  array<size_t, 3> pixels_;   //!< Plot size in pixels
  bool slice_color_overlaps_; //!< Show overlapping cells?
  int slice_level_ {-1};      //!< Plot universe level
private:
};

template<class T>
T SlicePlotBase::get_map() const
{

  size_t width = pixels_[0];
  size_t height = pixels_[1];

  // get pixel size
  double in_pixel = (width_[0]) / static_cast<double>(width);
  double out_pixel = (width_[1]) / static_cast<double>(height);

  // size data array
  T data(width, height);

  // setup basis indices and initial position centered on pixel
  int in_i, out_i;
  Position xyz = origin_;
  switch (basis_) {
  case PlotBasis::xy:
    in_i = 0;
    out_i = 1;
    break;
  case PlotBasis::xz:
    in_i = 0;
    out_i = 2;
    break;
  case PlotBasis::yz:
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
  Direction dir = {1. / std::sqrt(2.), 1. / std::sqrt(2.), 0.0};

#pragma omp parallel
  {
    GeometryState p;
    p.r() = xyz;
    p.u() = dir;
    p.coord(0).universe = model::root_universe;
    int level = slice_level_;
    int j {};

#pragma omp for
    for (int y = 0; y < height; y++) {
      p.r()[out_i] = xyz[out_i] - out_pixel * y;
      for (int x = 0; x < width; x++) {
        p.r()[in_i] = xyz[in_i] + in_pixel * x;
        p.n_coord() = 1;
        // local variables
        bool found_cell = exhaustive_find_cell(p);
        j = p.n_coord() - 1;
        if (level >= 0) {
          j = level;
        }
        if (found_cell) {
          data.set_value(y, x, p, j);
        }
        if (slice_color_overlaps_ && check_cell_overlap(p, false)) {
          data.set_overlap(y, x);
        }
      } // inner for
    }
  }

  return data;
}

// Represents either a voxel or pixel plot
class Plot : public PlottableInterface, public SlicePlotBase {

public:
  enum class PlotType { slice = 1, voxel = 2 };

  Plot(pugi::xml_node plot, PlotType type);

private:
  void set_output_path(pugi::xml_node plot_node);
  void set_basis(pugi::xml_node plot_node);
  void set_origin(pugi::xml_node plot_node);
  void set_width(pugi::xml_node plot_node);
  void set_meshlines(pugi::xml_node plot_node);

public:
  // Add mesh lines to ImageData
  void draw_mesh_lines(ImageData& data) const;
  void create_image() const;
  void create_voxel() const;

  virtual void create_output() const;
  virtual void print_info() const;

  PlotType type_;                 //!< Plot type (Slice/Voxel)
  int meshlines_width_;           //!< Width of lines added to the plot
  int index_meshlines_mesh_ {-1}; //!< Index of the mesh to draw on the plot
  RGBColor meshlines_color_;      //!< Color of meshlines on the plot
};

/**
 * @class RaytracePlot
 * @brief Base class for plots that generate images through ray tracing.
 *
 * This class serves as a base for plots that create their visuals by tracing
 * rays from a camera through the problem geometry. It inherits from
 * PlottableInterface, ensuring that it provides an implementation for
 * generating output specific to ray-traced visualization. ProjectionPlot
 * and PhongPlot provide concrete implementations of this class.
 */
class RayTracePlot : public PlottableInterface {
public:
  RayTracePlot(pugi::xml_node plot);

  // Standard getters. No setting since it's done from XML.
  const Position& camera_position() const { return camera_position_; }
  const Position& look_at() const { return look_at_; }
  const double& horizontal_field_of_view() const
  {
    return horizontal_field_of_view_;
  }

  virtual void print_info() const;

  /* If starting the particle from outside the geometry, we have to
   * find a distance to the boundary in a non-standard surface intersection
   * check. It's an exhaustive search over surfaces in the top-level universe.
   */
  static int advance_to_boundary_from_void(GeometryState& p);

protected:
  Direction camera_x_axis() const
  {
    return {camera_to_model_[0], camera_to_model_[3], camera_to_model_[6]};
  }

  Direction camera_y_axis() const
  {
    return {camera_to_model_[1], camera_to_model_[4], camera_to_model_[7]};
  }

  Direction camera_z_axis() const
  {
    return {camera_to_model_[2], camera_to_model_[5], camera_to_model_[8]};
  }

  void set_output_path(pugi::xml_node plot_node);

  /*
   * Gets the starting position and direction for the pixel corresponding
   * to this horizontal and vertical position.
   */
  std::pair<Position, Direction> get_pixel_ray(int horiz, int vert) const;

  std::array<int, 2> pixels_; // pixel dimension of resulting image

private:
  void set_look_at(pugi::xml_node node);
  void set_camera_position(pugi::xml_node node);
  void set_field_of_view(pugi::xml_node node);
  void set_pixels(pugi::xml_node node);
  void set_orthographic_width(pugi::xml_node node);

  double horizontal_field_of_view_ {70.0}; // horiz. f.o.v. in degrees
  Position camera_position_;               // where camera is
  Position look_at_; // point camera is centered looking at

  Direction up_ {0.0, 0.0, 1.0}; // which way is up

  /* The horizontal thickness, if using an orthographic projection.
   * If set to zero, we assume using a perspective projection.
   */
  double orthographic_width_ {0.0};

  /*
   * Cached camera-to-model matrix with column vectors of axes. The x-axis is
   * the vector between the camera_position_ and look_at_; the y-axis is the
   * cross product of the x-axis with the up_ vector, and the z-axis is the
   * cross product of the x and y axes.
   */
  std::array<double, 9> camera_to_model_;
};

class ProjectionRay;

/**
 * @class ProjectionPlot
 * @brief Creates plots that are like colorful x-ray imaging
 *
 * ProjectionPlot is a specialized form of RayTracePlot designed for creating
 * projection plots. This involves tracing rays from a camera through the
 * problem geometry and rendering the results based on depth of penetration
 * through materials or cells and their colors.
 */
class ProjectionPlot : public RayTracePlot {

  friend class ProjectionRay;

public:
  ProjectionPlot(pugi::xml_node plot);

  virtual void create_output() const;
  virtual void print_info() const;

private:
  void set_opacities(pugi::xml_node node);
  void set_wireframe_thickness(pugi::xml_node node);
  void set_wireframe_ids(pugi::xml_node node);
  void set_wireframe_color(pugi::xml_node node);

  /* Checks if a vector of two TrackSegments is equivalent. We define this
   * to mean not having matching intersection lengths, but rather having
   * a matching sequence of surface/cell/material intersections.
   */
  struct TrackSegment;
  bool trackstack_equivalent(const vector<TrackSegment>& track1,
    const vector<TrackSegment>& track2) const;

  /* Used for drawing wireframe and colors. We record the list of
   * surface/cell/material intersections and the corresponding lengths as a ray
   * traverses the geometry, then color by iterating in reverse.
   */
  struct TrackSegment {
    int id;        // material or cell ID (which is being colored)
    double length; // length of this track intersection

    /* Recording this allows us to draw edges on the wireframe. For instance
     * if two surfaces bound a single cell, it allows drawing that sharp edge
     * where the surfaces intersect.
     */
    int surface; // last surface ID intersected in this segment
    TrackSegment(int id_a, double length_a, int surface_a)
      : id(id_a), length(length_a), surface(surface_a)
    {}
  };

  // which color IDs should be wireframed. If empty, all cells are wireframed.
  vector<int> wireframe_ids_;

  // Thickness of the wireframe lines. Can set to zero for no wireframe.
  int wireframe_thickness_ {1};

  RGBColor wireframe_color_ {BLACK}; // wireframe color
  vector<double> xs_; // macro cross section values for cell volume rendering
};

/**
 * @class PhongPlot
 * @brief Plots 3D objects as the eye might see them.
 *
 * Plots a geometry with single-scattered Phong lighting
 * plus a diffuse lighting contribution. The result is a
 * physically reasonable, aesthetic 3D view of a geometry.
 */
class PhongPlot : public RayTracePlot {
  friend class PhongRay;

public:
  PhongPlot(pugi::xml_node plot);

  virtual void create_output() const;
  virtual void print_info() const;

private:
  void set_opaque_ids(pugi::xml_node node);
  void set_light_position(pugi::xml_node node);
  void set_diffuse_fraction(pugi::xml_node node);

  std::set<int> opaque_ids_;

  double diffuse_fraction_ {0.1};

  // By default, the light is at the camera unless otherwise specified.
  Position light_location_;
};

// Base class that implements ray tracing logic, not necessarily through
// defined regions of the geometry but also outside of it.
class Ray : public GeometryState {

public:
  Ray(Position r, Direction u) { init_from_r_u(r, u); }

  // Called at every surface intersection within the model
  virtual void on_intersection() = 0;

  /*
   * Traces the ray through the geometry, calling on_intersection
   * at every surface boundary.
   */
  void trace();

  /*
   * Implementation code has read-only access to the distance
   * between the ray and the next cell it's intersecting
   */
  const BoundaryInfo& dist() { return dist_; }

  const int& first_surface() { return first_surface_; }

  const int& i_surface() { return i_surface_; }

  // Stops the ray and exits tracing when called from on_intersection
  void stop() { stop_ = true; }

  // Sets the dist_ variable
  void compute_distance();

protected:
  // Records how far the ray has traveled
  double traversal_distance_ {0.0};

private:
  // Max intersections before we assume ray tracing is caught in an infinite
  // loop:
  static const int MAX_INTERSECTIONS = 1000000;

  bool hit_something_ {false};
  bool intersection_found_ {true};
  bool stop_ {false};

  unsigned event_counter_ {0};

  // Records the first intersected surface on the model
  int first_surface_ {-1};
  int i_surface_ {-1};

  BoundaryInfo dist_;
};

class ProjectionRay : public Ray {
public:
  ProjectionRay(Position r, Direction u, const ProjectionPlot& plot,
    vector<ProjectionPlot::TrackSegment>& line_segments)
    : Ray(r, u), plot_(plot), line_segments_(line_segments)
  {}

  virtual void on_intersection() override;

private:
  /* Store a reference to the plot object which is running this ray, in order
   * to access some of the plot settings which influence the behavior where
   * intersections are.
   */
  const ProjectionPlot& plot_;

  /* The ray runs through the geometry, and records the lengths of ray segments
   * and cells they lie in along the way.
   */
  vector<ProjectionPlot::TrackSegment>& line_segments_;
};

class PhongRay : public Ray {
public:
  PhongRay(Position r, Direction u, const PhongPlot& plot)
    : Ray(r, u), plot_(plot)
  {
    result_color_ = plot_.not_found_;
  }

  virtual void on_intersection() override;

  const RGBColor& result_color() { return result_color_; }

private:
  const PhongPlot& plot_;

  /* After the ray is reflected, it is moving towards the
   * camera. It does that in order to see if the exposed surface
   * is shadowed by something else.
   */
  bool reflected_ {false};

  // Have to record the first hit ID, so that if the region
  // does get shadowed, we recall what its color should be
  // when tracing from the surface to the light.
  int orig_hit_id_ {-1};

  RGBColor result_color_;
};

//===============================================================================
// Non-member functions
//===============================================================================

/* Write a PPM image
 * filename - name of output file
 * data - image data to write
 */
void output_ppm(const std::string& filename, const ImageData& data);

#ifdef USE_LIBPNG
/* Write a PNG image
 * filename - name of output file
 * data - image data to write
 */
void output_png(const std::string& filename, const ImageData& data);
#endif

//! Initialize a voxel file
//! \param[in] id of an open hdf5 file
//! \param[in] dimensions of the voxel file (dx, dy, dz)
//! \param[out] dataspace pointer to voxel data
//! \param[out] dataset pointer to voxesl data
//! \param[out] pointer to memory space of voxel data
void voxel_init(hid_t file_id, const hsize_t* dims, hid_t* dspace, hid_t* dset,
  hid_t* memspace);

//! Write a section of the voxel data to hdf5
//! \param[in] voxel slice
//! \param[out] dataspace pointer to voxel data
//! \param[out] dataset pointer to voxesl data
//! \param[out] pointer to data to write
void voxel_write_slice(
  int x, hid_t dspace, hid_t dset, hid_t memspace, void* buf);

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

//! Read plot specifications from an XML Node
//! \param[in] XML node containing plot info
void read_plots_xml(pugi::xml_node root);

//! Clear memory
void free_memory_plot();

//! Create a randomly generated RGB color
//! \return RGBColor with random value
RGBColor random_color();

} // namespace openmc
#endif // OPENMC_PLOT_H
