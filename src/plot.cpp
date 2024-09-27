#include "openmc/plot.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>

#include "xtensor/xmanipulation.hpp"
#include "xtensor/xview.hpp"
#include <fmt/core.h>
#include <fmt/ostream.h>
#ifdef USE_LIBPNG
#include <png.h>
#endif

#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/dagmc.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/openmp_interface.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/progress_bar.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/string_utils.h"

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

constexpr int PLOT_LEVEL_LOWEST {-1}; //!< lower bound on plot universe level
constexpr int32_t NOT_FOUND {-2};
constexpr int32_t OVERLAP {-3};

IdData::IdData(size_t h_res, size_t v_res) : data_({v_res, h_res, 3}, NOT_FOUND)
{}

void IdData::set_value(size_t y, size_t x, const GeometryState& p, int level)
{
  // set cell data
  if (p.n_coord() <= level) {
    data_(y, x, 0) = NOT_FOUND;
    data_(y, x, 1) = NOT_FOUND;
  } else {
    data_(y, x, 0) = model::cells.at(p.coord(level).cell)->id_;
    data_(y, x, 1) = level == p.n_coord() - 1
                       ? p.cell_instance()
                       : cell_instance_at_level(p, level);
  }

  // set material data
  Cell* c = model::cells.at(p.lowest_coord().cell).get();
  if (p.material() == MATERIAL_VOID) {
    data_(y, x, 2) = MATERIAL_VOID;
    return;
  } else if (c->type_ == Fill::MATERIAL) {
    Material* m = model::materials.at(p.material()).get();
    data_(y, x, 2) = m->id_;
  }
}

void IdData::set_overlap(size_t y, size_t x)
{
  xt::view(data_, y, x, xt::all()) = OVERLAP;
}

PropertyData::PropertyData(size_t h_res, size_t v_res)
  : data_({v_res, h_res, 2}, NOT_FOUND)
{}

void PropertyData::set_value(
  size_t y, size_t x, const GeometryState& p, int level)
{
  Cell* c = model::cells.at(p.lowest_coord().cell).get();
  data_(y, x, 0) = (p.sqrtkT() * p.sqrtkT()) / K_BOLTZMANN;
  if (c->type_ != Fill::UNIVERSE && p.material() != MATERIAL_VOID) {
    Material* m = model::materials.at(p.material()).get();
    data_(y, x, 1) = m->density_gpcc_;
  }
}

void PropertyData::set_overlap(size_t y, size_t x)
{
  data_(y, x) = OVERLAP;
}

//==============================================================================
// Global variables
//==============================================================================

namespace model {

std::unordered_map<int, int> plot_map;
vector<std::unique_ptr<PlottableInterface>> plots;
uint64_t plotter_seed = 1;

} // namespace model

//==============================================================================
// RUN_PLOT controls the logic for making one or many plots
//==============================================================================

extern "C" int openmc_plot_geometry()
{

  for (auto& pl : model::plots) {
    write_message(5, "Processing plot {}: {}...", pl->id(), pl->path_plot());
    pl->create_output();
  }

  return 0;
}

void Plot::create_output() const
{
  if (PlotType::slice == type_) {
    // create 2D image
    create_image();
  } else if (PlotType::voxel == type_) {
    // create voxel file for 3D viewing
    create_voxel();
  }
}

void Plot::print_info() const
{
  // Plot type
  if (PlotType::slice == type_) {
    fmt::print("Plot Type: Slice\n");
  } else if (PlotType::voxel == type_) {
    fmt::print("Plot Type: Voxel\n");
  }

  // Plot parameters
  fmt::print("Origin: {} {} {}\n", origin_[0], origin_[1], origin_[2]);

  if (PlotType::slice == type_) {
    fmt::print("Width: {:4} {:4}\n", width_[0], width_[1]);
  } else if (PlotType::voxel == type_) {
    fmt::print("Width: {:4} {:4} {:4}\n", width_[0], width_[1], width_[2]);
  }

  if (PlotColorBy::cells == color_by_) {
    fmt::print("Coloring: Cells\n");
  } else if (PlotColorBy::mats == color_by_) {
    fmt::print("Coloring: Materials\n");
  }

  if (PlotType::slice == type_) {
    switch (basis_) {
    case PlotBasis::xy:
      fmt::print("Basis: XY\n");
      break;
    case PlotBasis::xz:
      fmt::print("Basis: XZ\n");
      break;
    case PlotBasis::yz:
      fmt::print("Basis: YZ\n");
      break;
    }
    fmt::print("Pixels: {} {}\n", pixels_[0], pixels_[1]);
  } else if (PlotType::voxel == type_) {
    fmt::print("Voxels: {} {} {}\n", pixels_[0], pixels_[1], pixels_[2]);
  }
}

void read_plots_xml()
{
  // Check if plots.xml exists; this is only necessary when the plot runmode is
  // initiated. Otherwise, we want to read plots.xml because it may be called
  // later via the API. In that case, its ok for a plots.xml to not exist
  std::string filename = settings::path_input + "plots.xml";
  if (!file_exists(filename) && settings::run_mode == RunMode::PLOTTING) {
    fatal_error(fmt::format("Plots XML file '{}' does not exist!", filename));
  }

  write_message("Reading plot XML file...", 5);

  // Parse plots.xml file
  pugi::xml_document doc;
  doc.load_file(filename.c_str());

  pugi::xml_node root = doc.document_element();

  read_plots_xml(root);
}

void read_plots_xml(pugi::xml_node root)
{
  for (auto node : root.children("plot")) {
    std::string id_string = get_node_value(node, "id", true);
    int id = std::stoi(id_string);
    if (check_for_node(node, "type")) {
      std::string type_str = get_node_value(node, "type", true);
      if (type_str == "slice")
        model::plots.emplace_back(
          std::make_unique<Plot>(node, Plot::PlotType::slice));
      else if (type_str == "voxel")
        model::plots.emplace_back(
          std::make_unique<Plot>(node, Plot::PlotType::voxel));
      else if (type_str == "projection")
        model::plots.emplace_back(std::make_unique<ProjectionPlot>(node));
      else
        fatal_error(
          fmt::format("Unsupported plot type '{}' in plot {}", type_str, id));

      model::plot_map[model::plots.back()->id()] = model::plots.size() - 1;
    } else {
      fatal_error(fmt::format("Must specify plot type in plot {}", id));
    }
  }
}

void free_memory_plot()
{
  model::plots.clear();
  model::plot_map.clear();
}

// creates an image based on user input from a plots.xml <plot>
// specification in the PNG/PPM format
void Plot::create_image() const
{

  size_t width = pixels_[0];
  size_t height = pixels_[1];

  ImageData data({width, height}, not_found_);

  // generate ids for the plot
  auto ids = get_map<IdData>();

  // assign colors
  for (size_t y = 0; y < height; y++) {
    for (size_t x = 0; x < width; x++) {
      int idx = color_by_ == PlotColorBy::cells ? 0 : 2;
      auto id = ids.data_(y, x, idx);
      // no setting needed if not found
      if (id == NOT_FOUND) {
        continue;
      }
      if (id == OVERLAP) {
        data(x, y) = overlap_color_;
        continue;
      }
      if (PlotColorBy::cells == color_by_) {
        data(x, y) = colors_[model::cell_map[id]];
      } else if (PlotColorBy::mats == color_by_) {
        if (id == MATERIAL_VOID) {
          data(x, y) = WHITE;
          continue;
        }
        data(x, y) = colors_[model::material_map[id]];
      } // color_by if-else
    }   // x for loop
  }     // y for loop

  // draw mesh lines if present
  if (index_meshlines_mesh_ >= 0) {
    draw_mesh_lines(data);
  }

// create image file
#ifdef USE_LIBPNG
  output_png(path_plot(), data);
#else
  output_ppm(path_plot(), data);
#endif
}

void PlottableInterface::set_id(pugi::xml_node plot_node)
{
  // Copy data into plots
  if (check_for_node(plot_node, "id")) {
    id_ = std::stoi(get_node_value(plot_node, "id"));
  } else {
    fatal_error("Must specify plot id in plots XML file.");
  }

  // Check to make sure 'id' hasn't been used
  if (model::plot_map.find(id_) != model::plot_map.end()) {
    fatal_error(
      fmt::format("Two or more plots use the same unique ID: {}", id_));
  }
}

// Checks if png or ppm is already present
bool file_extension_present(
  const std::string& filename, const std::string& extension)
{
  std::string file_extension_if_present =
    filename.substr(filename.find_last_of(".") + 1);
  if (file_extension_if_present == extension)
    return true;
  return false;
}

void Plot::set_output_path(pugi::xml_node plot_node)
{
  // Set output file path
  std::string filename;

  if (check_for_node(plot_node, "filename")) {
    filename = get_node_value(plot_node, "filename");
  } else {
    filename = fmt::format("plot_{}", id());
  }
  const std::string dir_if_present =
    filename.substr(0, filename.find_last_of("/") + 1);
  if (dir_if_present.size() > 0 && !dir_exists(dir_if_present)) {
    fatal_error(fmt::format("Directory '{}' does not exist!", dir_if_present));
  }
  // add appropriate file extension to name
  switch (type_) {
  case PlotType::slice:
#ifdef USE_LIBPNG
    if (!file_extension_present(filename, "png"))
      filename.append(".png");
#else
    if (!file_extension_present(filename, "ppm"))
      filename.append(".ppm");
#endif
    break;
  case PlotType::voxel:
    if (!file_extension_present(filename, "h5"))
      filename.append(".h5");
    break;
  }

  path_plot_ = filename;

  // Copy plot pixel size
  vector<int> pxls = get_node_array<int>(plot_node, "pixels");
  if (PlotType::slice == type_) {
    if (pxls.size() == 2) {
      pixels_[0] = pxls[0];
      pixels_[1] = pxls[1];
    } else {
      fatal_error(
        fmt::format("<pixels> must be length 2 in slice plot {}", id()));
    }
  } else if (PlotType::voxel == type_) {
    if (pxls.size() == 3) {
      pixels_[0] = pxls[0];
      pixels_[1] = pxls[1];
      pixels_[2] = pxls[2];
    } else {
      fatal_error(
        fmt::format("<pixels> must be length 3 in voxel plot {}", id()));
    }
  }
}

void PlottableInterface::set_bg_color(pugi::xml_node plot_node)
{
  // Copy plot background color
  if (check_for_node(plot_node, "background")) {
    vector<int> bg_rgb = get_node_array<int>(plot_node, "background");
    if (bg_rgb.size() == 3) {
      not_found_ = bg_rgb;
    } else {
      fatal_error(fmt::format("Bad background RGB in plot {}", id()));
    }
  }
}

void Plot::set_basis(pugi::xml_node plot_node)
{
  // Copy plot basis
  if (PlotType::slice == type_) {
    std::string pl_basis = "xy";
    if (check_for_node(plot_node, "basis")) {
      pl_basis = get_node_value(plot_node, "basis", true);
    }
    if ("xy" == pl_basis) {
      basis_ = PlotBasis::xy;
    } else if ("xz" == pl_basis) {
      basis_ = PlotBasis::xz;
    } else if ("yz" == pl_basis) {
      basis_ = PlotBasis::yz;
    } else {
      fatal_error(
        fmt::format("Unsupported plot basis '{}' in plot {}", pl_basis, id()));
    }
  }
}

void Plot::set_origin(pugi::xml_node plot_node)
{
  // Copy plotting origin
  auto pl_origin = get_node_array<double>(plot_node, "origin");
  if (pl_origin.size() == 3) {
    origin_ = pl_origin;
  } else {
    fatal_error(fmt::format("Origin must be length 3 in plot {}", id()));
  }
}

void Plot::set_width(pugi::xml_node plot_node)
{
  // Copy plotting width
  vector<double> pl_width = get_node_array<double>(plot_node, "width");
  if (PlotType::slice == type_) {
    if (pl_width.size() == 2) {
      width_.x = pl_width[0];
      width_.y = pl_width[1];
    } else {
      fatal_error(
        fmt::format("<width> must be length 2 in slice plot {}", id()));
    }
  } else if (PlotType::voxel == type_) {
    if (pl_width.size() == 3) {
      pl_width = get_node_array<double>(plot_node, "width");
      width_ = pl_width;
    } else {
      fatal_error(
        fmt::format("<width> must be length 3 in voxel plot {}", id()));
    }
  }
}

void PlottableInterface::set_universe(pugi::xml_node plot_node)
{
  // Copy plot universe level
  if (check_for_node(plot_node, "level")) {
    level_ = std::stoi(get_node_value(plot_node, "level"));
    if (level_ < 0) {
      fatal_error(fmt::format("Bad universe level in plot {}", id()));
    }
  } else {
    level_ = PLOT_LEVEL_LOWEST;
  }
}

void PlottableInterface::set_default_colors(pugi::xml_node plot_node)
{
  // Copy plot color type and initialize all colors randomly
  std::string pl_color_by = "cell";
  if (check_for_node(plot_node, "color_by")) {
    pl_color_by = get_node_value(plot_node, "color_by", true);
  }
  if ("cell" == pl_color_by) {
    color_by_ = PlotColorBy::cells;
    colors_.resize(model::cells.size());
  } else if ("material" == pl_color_by) {
    color_by_ = PlotColorBy::mats;
    colors_.resize(model::materials.size());
  } else {
    fatal_error(fmt::format(
      "Unsupported plot color type '{}' in plot {}", pl_color_by, id()));
  }

  for (auto& c : colors_) {
    c = random_color();
    // make sure we don't interfere with some default colors
    while (c == RED || c == WHITE) {
      c = random_color();
    }
  }
}

void PlottableInterface::set_user_colors(pugi::xml_node plot_node)
{
  for (auto cn : plot_node.children("color")) {
    // Make sure 3 values are specified for RGB
    vector<int> user_rgb = get_node_array<int>(cn, "rgb");
    if (user_rgb.size() != 3) {
      fatal_error(fmt::format("Bad RGB in plot {}", id()));
    }
    // Ensure that there is an id for this color specification
    int col_id;
    if (check_for_node(cn, "id")) {
      col_id = std::stoi(get_node_value(cn, "id"));
    } else {
      fatal_error(fmt::format(
        "Must specify id for color specification in plot {}", id()));
    }
    // Add RGB
    if (PlotColorBy::cells == color_by_) {
      if (model::cell_map.find(col_id) != model::cell_map.end()) {
        col_id = model::cell_map[col_id];
        colors_[col_id] = user_rgb;
      } else {
        warning(fmt::format(
          "Could not find cell {} specified in plot {}", col_id, id()));
      }
    } else if (PlotColorBy::mats == color_by_) {
      if (model::material_map.find(col_id) != model::material_map.end()) {
        col_id = model::material_map[col_id];
        colors_[col_id] = user_rgb;
      } else {
        warning(fmt::format(
          "Could not find material {} specified in plot {}", col_id, id()));
      }
    }
  } // color node loop
}

void Plot::set_meshlines(pugi::xml_node plot_node)
{
  // Deal with meshlines
  pugi::xpath_node_set mesh_line_nodes = plot_node.select_nodes("meshlines");

  if (!mesh_line_nodes.empty()) {
    if (PlotType::voxel == type_) {
      warning(fmt::format("Meshlines ignored in voxel plot {}", id()));
    }

    if (mesh_line_nodes.size() == 1) {
      // Get first meshline node
      pugi::xml_node meshlines_node = mesh_line_nodes[0].node();

      // Check mesh type
      std::string meshtype;
      if (check_for_node(meshlines_node, "meshtype")) {
        meshtype = get_node_value(meshlines_node, "meshtype");
      } else {
        fatal_error(fmt::format(
          "Must specify a meshtype for meshlines specification in plot {}",
          id()));
      }

      // Ensure that there is a linewidth for this meshlines specification
      std::string meshline_width;
      if (check_for_node(meshlines_node, "linewidth")) {
        meshline_width = get_node_value(meshlines_node, "linewidth");
        meshlines_width_ = std::stoi(meshline_width);
      } else {
        fatal_error(fmt::format(
          "Must specify a linewidth for meshlines specification in plot {}",
          id()));
      }

      // Check for color
      if (check_for_node(meshlines_node, "color")) {
        // Check and make sure 3 values are specified for RGB
        vector<int> ml_rgb = get_node_array<int>(meshlines_node, "color");
        if (ml_rgb.size() != 3) {
          fatal_error(
            fmt::format("Bad RGB for meshlines color in plot {}", id()));
        }
        meshlines_color_ = ml_rgb;
      }

      // Set mesh based on type
      if ("ufs" == meshtype) {
        if (!simulation::ufs_mesh) {
          fatal_error(
            fmt::format("No UFS mesh for meshlines on plot {}", id()));
        } else {
          for (int i = 0; i < model::meshes.size(); ++i) {
            if (const auto* m =
                  dynamic_cast<const RegularMesh*>(model::meshes[i].get())) {
              if (m == simulation::ufs_mesh) {
                index_meshlines_mesh_ = i;
              }
            }
          }
          if (index_meshlines_mesh_ == -1)
            fatal_error("Could not find the UFS mesh for meshlines plot");
        }
      } else if ("entropy" == meshtype) {
        if (!simulation::entropy_mesh) {
          fatal_error(
            fmt::format("No entropy mesh for meshlines on plot {}", id()));
        } else {
          for (int i = 0; i < model::meshes.size(); ++i) {
            if (const auto* m =
                  dynamic_cast<const RegularMesh*>(model::meshes[i].get())) {
              if (m == simulation::entropy_mesh) {
                index_meshlines_mesh_ = i;
              }
            }
          }
          if (index_meshlines_mesh_ == -1)
            fatal_error("Could not find the entropy mesh for meshlines plot");
        }
      } else if ("tally" == meshtype) {
        // Ensure that there is a mesh id if the type is tally
        int tally_mesh_id;
        if (check_for_node(meshlines_node, "id")) {
          tally_mesh_id = std::stoi(get_node_value(meshlines_node, "id"));
        } else {
          std::stringstream err_msg;
          fatal_error(fmt::format("Must specify a mesh id for meshlines tally "
                                  "mesh specification in plot {}",
            id()));
        }
        // find the tally index
        int idx;
        int err = openmc_get_mesh_index(tally_mesh_id, &idx);
        if (err != 0) {
          fatal_error(fmt::format("Could not find mesh {} specified in "
                                  "meshlines for plot {}",
            tally_mesh_id, id()));
        }
        index_meshlines_mesh_ = idx;
      } else {
        fatal_error(fmt::format("Invalid type for meshlines on plot {}", id()));
      }
    } else {
      fatal_error(fmt::format("Mutliple meshlines specified in plot {}", id()));
    }
  }
}

void PlottableInterface::set_mask(pugi::xml_node plot_node)
{
  // Deal with masks
  pugi::xpath_node_set mask_nodes = plot_node.select_nodes("mask");

  if (!mask_nodes.empty()) {
    if (mask_nodes.size() == 1) {
      // Get pointer to mask
      pugi::xml_node mask_node = mask_nodes[0].node();

      // Determine how many components there are and allocate
      vector<int> iarray = get_node_array<int>(mask_node, "components");
      if (iarray.size() == 0) {
        fatal_error(
          fmt::format("Missing <components> in mask of plot {}", id()));
      }

      // First we need to change the user-specified identifiers to indices
      // in the cell and material arrays
      for (auto& col_id : iarray) {
        if (PlotColorBy::cells == color_by_) {
          if (model::cell_map.find(col_id) != model::cell_map.end()) {
            col_id = model::cell_map[col_id];
          } else {
            fatal_error(fmt::format("Could not find cell {} specified in the "
                                    "mask in plot {}",
              col_id, id()));
          }
        } else if (PlotColorBy::mats == color_by_) {
          if (model::material_map.find(col_id) != model::material_map.end()) {
            col_id = model::material_map[col_id];
          } else {
            fatal_error(fmt::format("Could not find material {} specified in "
                                    "the mask in plot {}",
              col_id, id()));
          }
        }
      }

      // Alter colors based on mask information
      for (int j = 0; j < colors_.size(); j++) {
        if (contains(iarray, j)) {
          if (check_for_node(mask_node, "background")) {
            vector<int> bg_rgb = get_node_array<int>(mask_node, "background");
            colors_[j] = bg_rgb;
          } else {
            colors_[j] = WHITE;
          }
        }
      }

    } else {
      fatal_error(fmt::format("Mutliple masks specified in plot {}", id()));
    }
  }
}

void PlottableInterface::set_overlap_color(pugi::xml_node plot_node)
{
  color_overlaps_ = false;
  if (check_for_node(plot_node, "show_overlaps")) {
    color_overlaps_ = get_node_value_bool(plot_node, "show_overlaps");
    // check for custom overlap color
    if (check_for_node(plot_node, "overlap_color")) {
      if (!color_overlaps_) {
        warning(fmt::format(
          "Overlap color specified in plot {} but overlaps won't be shown.",
          id()));
      }
      vector<int> olap_clr = get_node_array<int>(plot_node, "overlap_color");
      if (olap_clr.size() == 3) {
        overlap_color_ = olap_clr;
      } else {
        fatal_error(fmt::format("Bad overlap RGB in plot {}", id()));
      }
    }
  }

  // make sure we allocate the vector for counting overlap checks if
  // they're going to be plotted
  if (color_overlaps_ && settings::run_mode == RunMode::PLOTTING) {
    settings::check_overlaps = true;
    model::overlap_check_count.resize(model::cells.size(), 0);
  }
}

PlottableInterface::PlottableInterface(pugi::xml_node plot_node)
{
  set_id(plot_node);
  set_bg_color(plot_node);
  set_universe(plot_node);
  set_default_colors(plot_node);
  set_user_colors(plot_node);
  set_mask(plot_node);
  set_overlap_color(plot_node);
}

Plot::Plot(pugi::xml_node plot_node, PlotType type)
  : PlottableInterface(plot_node), type_(type), index_meshlines_mesh_ {-1}
{
  set_output_path(plot_node);
  set_basis(plot_node);
  set_origin(plot_node);
  set_width(plot_node);
  set_meshlines(plot_node);
  slice_level_ = level_; // Copy level employed in SlicePlotBase::get_map
  slice_color_overlaps_ = color_overlaps_;
}

//==============================================================================
// OUTPUT_PPM writes out a previously generated image to a PPM file
//==============================================================================

void output_ppm(const std::string& filename, const ImageData& data)
{
  // Open PPM file for writing
  std::string fname = filename;
  fname = strtrim(fname);
  std::ofstream of;

  of.open(fname);

  // Write header
  of << "P6\n";
  of << data.shape()[0] << " " << data.shape()[1] << "\n";
  of << "255\n";
  of.close();

  of.open(fname, std::ios::binary | std::ios::app);
  // Write color for each pixel
  for (int y = 0; y < data.shape()[1]; y++) {
    for (int x = 0; x < data.shape()[0]; x++) {
      RGBColor rgb = data(x, y);
      of << rgb.red << rgb.green << rgb.blue;
    }
  }
  of << "\n";
}

//==============================================================================
// OUTPUT_PNG writes out a previously generated image to a PNG file
//==============================================================================

#ifdef USE_LIBPNG
void output_png(const std::string& filename, const ImageData& data)
{
  // Open PNG file for writing
  std::string fname = filename;
  fname = strtrim(fname);
  auto fp = std::fopen(fname.c_str(), "wb");

  // Initialize write and info structures
  auto png_ptr =
    png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
  auto info_ptr = png_create_info_struct(png_ptr);

  // Setup exception handling
  if (setjmp(png_jmpbuf(png_ptr)))
    fatal_error("Error during png creation");

  png_init_io(png_ptr, fp);

  // Write header (8 bit colour depth)
  int width = data.shape()[0];
  int height = data.shape()[1];
  png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB,
    PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
  png_write_info(png_ptr, info_ptr);

  // Allocate memory for one row (3 bytes per pixel - RGB)
  std::vector<png_byte> row(3 * width);

  // Write color for each pixel
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      RGBColor rgb = data(x, y);
      row[3 * x] = rgb.red;
      row[3 * x + 1] = rgb.green;
      row[3 * x + 2] = rgb.blue;
    }
    png_write_row(png_ptr, row.data());
  }

  // End write
  png_write_end(png_ptr, nullptr);

  // Clean up data structures
  std::fclose(fp);
  png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
  png_destroy_write_struct(&png_ptr, &info_ptr);
}
#endif

//==============================================================================
// DRAW_MESH_LINES draws mesh line boundaries on an image
//==============================================================================

void Plot::draw_mesh_lines(ImageData& data) const
{
  RGBColor rgb;
  rgb = meshlines_color_;

  int ax1, ax2;
  switch (basis_) {
  case PlotBasis::xy:
    ax1 = 0;
    ax2 = 1;
    break;
  case PlotBasis::xz:
    ax1 = 0;
    ax2 = 2;
    break;
  case PlotBasis::yz:
    ax1 = 1;
    ax2 = 2;
    break;
  default:
    UNREACHABLE();
  }

  Position ll_plot {origin_};
  Position ur_plot {origin_};

  ll_plot[ax1] -= width_[0] / 2.;
  ll_plot[ax2] -= width_[1] / 2.;
  ur_plot[ax1] += width_[0] / 2.;
  ur_plot[ax2] += width_[1] / 2.;

  Position width = ur_plot - ll_plot;

  // Find the (axis-aligned) lines of the mesh that intersect this plot.
  auto axis_lines =
    model::meshes[index_meshlines_mesh_]->plot(ll_plot, ur_plot);

  // Find the bounds along the second axis (accounting for low-D meshes).
  int ax2_min, ax2_max;
  if (axis_lines.second.size() > 0) {
    double frac = (axis_lines.second.back() - ll_plot[ax2]) / width[ax2];
    ax2_min = (1.0 - frac) * pixels_[1];
    if (ax2_min < 0)
      ax2_min = 0;
    frac = (axis_lines.second.front() - ll_plot[ax2]) / width[ax2];
    ax2_max = (1.0 - frac) * pixels_[1];
    if (ax2_max > pixels_[1])
      ax2_max = pixels_[1];
  } else {
    ax2_min = 0;
    ax2_max = pixels_[1];
  }

  // Iterate across the first axis and draw lines.
  for (auto ax1_val : axis_lines.first) {
    double frac = (ax1_val - ll_plot[ax1]) / width[ax1];
    int ax1_ind = frac * pixels_[0];
    for (int ax2_ind = ax2_min; ax2_ind < ax2_max; ++ax2_ind) {
      for (int plus = 0; plus <= meshlines_width_; plus++) {
        if (ax1_ind + plus >= 0 && ax1_ind + plus < pixels_[0])
          data(ax1_ind + plus, ax2_ind) = rgb;
        if (ax1_ind - plus >= 0 && ax1_ind - plus < pixels_[0])
          data(ax1_ind - plus, ax2_ind) = rgb;
      }
    }
  }

  // Find the bounds along the first axis.
  int ax1_min, ax1_max;
  if (axis_lines.first.size() > 0) {
    double frac = (axis_lines.first.front() - ll_plot[ax1]) / width[ax1];
    ax1_min = frac * pixels_[0];
    if (ax1_min < 0)
      ax1_min = 0;
    frac = (axis_lines.first.back() - ll_plot[ax1]) / width[ax1];
    ax1_max = frac * pixels_[0];
    if (ax1_max > pixels_[0])
      ax1_max = pixels_[0];
  } else {
    ax1_min = 0;
    ax1_max = pixels_[0];
  }

  // Iterate across the second axis and draw lines.
  for (auto ax2_val : axis_lines.second) {
    double frac = (ax2_val - ll_plot[ax2]) / width[ax2];
    int ax2_ind = (1.0 - frac) * pixels_[1];
    for (int ax1_ind = ax1_min; ax1_ind < ax1_max; ++ax1_ind) {
      for (int plus = 0; plus <= meshlines_width_; plus++) {
        if (ax2_ind + plus >= 0 && ax2_ind + plus < pixels_[1])
          data(ax1_ind, ax2_ind + plus) = rgb;
        if (ax2_ind - plus >= 0 && ax2_ind - plus < pixels_[1])
          data(ax1_ind, ax2_ind - plus) = rgb;
      }
    }
  }
}

/* outputs a binary file that can be input into silomesh for 3D geometry
 * visualization.  It works the same way as create_image by dragging a particle
 * across the geometry for the specified number of voxels. The first 3 int's in
 * the binary are the number of x, y, and z voxels.  The next 3 double's are
 * the widths of the voxels in the x, y, and z directions. The next 3 double's
 * are the x, y, and z coordinates of the lower left point. Finally the binary
 * is filled with entries of four int's each. Each 'row' in the binary contains
 * four int's: 3 for x,y,z position and 1 for cell or material id.  For 1
 * million voxels this produces a file of approximately 15MB.
 */
void Plot::create_voxel() const
{
  // compute voxel widths in each direction
  array<double, 3> vox;
  vox[0] = width_[0] / static_cast<double>(pixels_[0]);
  vox[1] = width_[1] / static_cast<double>(pixels_[1]);
  vox[2] = width_[2] / static_cast<double>(pixels_[2]);

  // initial particle position
  Position ll = origin_ - width_ / 2.;

  // Open binary plot file for writing
  std::ofstream of;
  std::string fname = std::string(path_plot_);
  fname = strtrim(fname);
  hid_t file_id = file_open(fname, 'w');

  // write header info
  write_attribute(file_id, "filetype", "voxel");
  write_attribute(file_id, "version", VERSION_VOXEL);
  write_attribute(file_id, "openmc_version", VERSION);

#ifdef GIT_SHA1
  write_attribute(file_id, "git_sha1", GIT_SHA1);
#endif

  // Write current date and time
  write_attribute(file_id, "date_and_time", time_stamp().c_str());
  array<int, 3> pixels;
  std::copy(pixels_.begin(), pixels_.end(), pixels.begin());
  write_attribute(file_id, "num_voxels", pixels);
  write_attribute(file_id, "voxel_width", vox);
  write_attribute(file_id, "lower_left", ll);

  // Create dataset for voxel data -- note that the dimensions are reversed
  // since we want the order in the file to be z, y, x
  hsize_t dims[3];
  dims[0] = pixels_[2];
  dims[1] = pixels_[1];
  dims[2] = pixels_[0];
  hid_t dspace, dset, memspace;
  voxel_init(file_id, &(dims[0]), &dspace, &dset, &memspace);

  SlicePlotBase pltbase;
  pltbase.width_ = width_;
  pltbase.origin_ = origin_;
  pltbase.basis_ = PlotBasis::xy;
  pltbase.pixels_ = pixels_;
  pltbase.slice_color_overlaps_ = color_overlaps_;

  ProgressBar pb;
  for (int z = 0; z < pixels_[2]; z++) {
    // update z coordinate
    pltbase.origin_.z = ll.z + z * vox[2];

    // generate ids using plotbase
    IdData ids = pltbase.get_map<IdData>();

    // select only cell/material ID data and flip the y-axis
    int idx = color_by_ == PlotColorBy::cells ? 0 : 2;
    xt::xtensor<int32_t, 2> data_slice =
      xt::view(ids.data_, xt::all(), xt::all(), idx);
    xt::xtensor<int32_t, 2> data_flipped = xt::flip(data_slice, 0);

    // Write to HDF5 dataset
    voxel_write_slice(z, dspace, dset, memspace, data_flipped.data());

    // update progress bar
    pb.set_value(
      100. * static_cast<double>(z + 1) / static_cast<double>((pixels_[2])));
  }

  voxel_finalize(dspace, dset, memspace);
  file_close(file_id);
}

void voxel_init(hid_t file_id, const hsize_t* dims, hid_t* dspace, hid_t* dset,
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

void voxel_write_slice(
  int x, hid_t dspace, hid_t dset, hid_t memspace, void* buf)
{
  hssize_t offset[3] {x, 0, 0};
  H5Soffset_simple(dspace, offset);
  H5Dwrite(dset, H5T_NATIVE_INT, memspace, dspace, H5P_DEFAULT, buf);
}

void voxel_finalize(hid_t dspace, hid_t dset, hid_t memspace)
{
  H5Dclose(dset);
  H5Sclose(dspace);
  H5Sclose(memspace);
}

RGBColor random_color(void)
{
  return {int(prn(&model::plotter_seed) * 255),
    int(prn(&model::plotter_seed) * 255), int(prn(&model::plotter_seed) * 255)};
}

ProjectionPlot::ProjectionPlot(pugi::xml_node node) : PlottableInterface(node)
{
  set_output_path(node);
  set_look_at(node);
  set_camera_position(node);
  set_field_of_view(node);
  set_pixels(node);
  set_opacities(node);
  set_orthographic_width(node);
  set_wireframe_thickness(node);
  set_wireframe_ids(node);
  set_wireframe_color(node);

  if (check_for_node(node, "orthographic_width") &&
      check_for_node(node, "field_of_view"))
    fatal_error("orthographic_width and field_of_view are mutually exclusive "
                "parameters.");
}

void ProjectionPlot::set_wireframe_color(pugi::xml_node plot_node)
{
  // Copy plot background color
  if (check_for_node(plot_node, "wireframe_color")) {
    vector<int> w_rgb = get_node_array<int>(plot_node, "wireframe_color");
    if (w_rgb.size() == 3) {
      wireframe_color_ = w_rgb;
    } else {
      fatal_error(fmt::format("Bad wireframe RGB in plot {}", id()));
    }
  }
}

void ProjectionPlot::set_output_path(pugi::xml_node node)
{
  // Set output file path
  std::string filename;

  if (check_for_node(node, "filename")) {
    filename = get_node_value(node, "filename");
  } else {
    filename = fmt::format("plot_{}", id());
  }

#ifdef USE_LIBPNG
  if (!file_extension_present(filename, "png"))
    filename.append(".png");
#else
  if (!file_extension_present(filename, "ppm"))
    filename.append(".ppm");
#endif
  path_plot_ = filename;
}

// Advances to the next boundary from outside the geometry
// Returns -1 if no intersection found, and the surface index
// if an intersection was found.
int ProjectionPlot::advance_to_boundary_from_void(GeometryState& p)
{
  constexpr double scoot = 1e-5;
  double min_dist = {INFINITY};
  auto coord = p.coord(0);
  Universe* uni = model::universes[model::root_universe].get();
  int intersected_surface = -1;
  for (auto c_i : uni->cells_) {
    auto dist = model::cells.at(c_i)->distance(coord.r, coord.u, 0, &p);
    if (dist.first < min_dist) {
      min_dist = dist.first;
      intersected_surface = dist.second;
    }
  }
  if (min_dist > 1e300)
    return -1;
  else { // advance the particle
    for (int j = 0; j < p.n_coord(); ++j)
      p.coord(j).r += (min_dist + scoot) * p.coord(j).u;
    return std::abs(intersected_surface);
  }
}

bool ProjectionPlot::trackstack_equivalent(
  const std::vector<TrackSegment>& track1,
  const std::vector<TrackSegment>& track2) const
{
  if (wireframe_ids_.empty()) {
    // Draw wireframe for all surfaces/cells/materials
    if (track1.size() != track2.size())
      return false;
    for (int i = 0; i < track1.size(); ++i) {
      if (track1[i].id != track2[i].id ||
          track1[i].surface != track2[i].surface) {
        return false;
      }
    }
    return true;
  } else {
    // This runs in O(nm) where n is the intersection stack size
    // and m is the number of IDs we are wireframing. A simpler
    // algorithm can likely be found.
    for (const int id : wireframe_ids_) {
      int t1_i = 0;
      int t2_i = 0;

      // Advance to first instance of the ID
      while (t1_i < track1.size() && t2_i < track2.size()) {
        while (t1_i < track1.size() && track1[t1_i].id != id)
          t1_i++;
        while (t2_i < track2.size() && track2[t2_i].id != id)
          t2_i++;

        // This one is really important!
        if ((t1_i == track1.size() && t2_i != track2.size()) ||
            (t1_i != track1.size() && t2_i == track2.size()))
          return false;
        if (t1_i == track1.size() && t2_i == track2.size())
          break;
        // Check if surface different
        if (track1[t1_i].surface != track2[t2_i].surface)
          return false;

        // Pretty sure this should not be used:
        // if (t2_i != track2.size() - 1 &&
        //     t1_i != track1.size() - 1 &&
        //     track1[t1_i+1].id != track2[t2_i+1].id) return false;
        if (t2_i != 0 && t1_i != 0 &&
            track1[t1_i - 1].surface != track2[t2_i - 1].surface)
          return false;

        // Check if neighboring cells are different
        // if (track1[t1_i ? t1_i - 1 : 0].id != track2[t2_i ? t2_i - 1 : 0].id)
        // return false; if (track1[t1_i < track1.size() - 1 ? t1_i + 1 : t1_i
        // ].id !=
        //    track2[t2_i < track2.size() - 1 ? t2_i + 1 : t2_i].id) return
        //    false;
        t1_i++, t2_i++;
      }
    }
    return true;
  }
}

void ProjectionPlot::create_output() const
{
  // Get centerline vector for camera-to-model. We create vectors around this
  // that form a pixel array, and then trace rays along that.
  auto up = up_ / up_.norm();
  Direction looking_direction = look_at_ - camera_position_;
  looking_direction /= looking_direction.norm();
  if (std::abs(std::abs(looking_direction.dot(up)) - 1.0) < 1e-9)
    fatal_error("Up vector cannot align with vector between camera position "
                "and look_at!");
  Direction cam_yaxis = looking_direction.cross(up);
  cam_yaxis /= cam_yaxis.norm();
  Direction cam_zaxis = cam_yaxis.cross(looking_direction);
  cam_zaxis /= cam_zaxis.norm();

  // Transformation matrix for directions
  std::vector<double> camera_to_model = {looking_direction.x, cam_yaxis.x,
    cam_zaxis.x, looking_direction.y, cam_yaxis.y, cam_zaxis.y,
    looking_direction.z, cam_yaxis.z, cam_zaxis.z};

  // Now we convert to the polar coordinate system with the polar angle
  // measuring the angle from the vector up_. Phi is the rotation about up_. For
  // now, up_ is hard-coded to be +z.
  constexpr double DEGREE_TO_RADIAN = M_PI / 180.0;
  double horiz_fov_radians = horizontal_field_of_view_ * DEGREE_TO_RADIAN;
  double p0 = static_cast<double>(pixels_[0]);
  double p1 = static_cast<double>(pixels_[1]);
  double vert_fov_radians = horiz_fov_radians * p1 / p0;
  double dphi = horiz_fov_radians / p0;
  double dmu = vert_fov_radians / p1;

  size_t width = pixels_[0];
  size_t height = pixels_[1];
  ImageData data({width, height}, not_found_);

  // This array marks where the initial wireframe was drawn.
  // We convolve it with a filter that gets adjusted with the
  // wireframe thickness in order to thicken the lines.
  xt::xtensor<int, 2> wireframe_initial({width, height}, 0);

  /* Holds all of the track segments for the current rendered line of pixels.
   * old_segments holds a copy of this_line_segments from the previous line.
   * By holding both we can check if the cell/material intersection stack
   * differs from the left or upper neighbor. This allows a robustly drawn
   * wireframe. If only checking the left pixel (which requires substantially
   * less memory), the wireframe tends to be spotty and be disconnected for
   * surface edges oriented horizontally in the rendering.
   *
   * Note that a vector of vectors is required rather than a 2-tensor,
   * since the stack size varies within each column.
   */
  const int n_threads = num_threads();
  std::vector<std::vector<std::vector<TrackSegment>>> this_line_segments(
    n_threads);
  for (int t = 0; t < n_threads; ++t) {
    this_line_segments[t].resize(pixels_[0]);
  }

  // The last thread writes to this, and the first thread reads from it.
  std::vector<std::vector<TrackSegment>> old_segments(pixels_[0]);

#pragma omp parallel
  {
    const int n_threads = num_threads();
    const int tid = thread_num();

    GeometryState p;
    p.u() = {1.0, 0.0, 0.0};

    int vert = tid;
    for (int iter = 0; iter <= pixels_[1] / n_threads; iter++) {

      // Save bottom line of current work chunk to compare against later
      // I used to have this inside the below if block, but it causes a
      // spurious line to be drawn at the bottom of the image. Not sure
      // why, but moving it here fixes things.
      if (tid == n_threads - 1)
        old_segments = this_line_segments[n_threads - 1];

      if (vert < pixels_[1]) {

        for (int horiz = 0; horiz < pixels_[0]; ++horiz) {

          // Projection mode below decides ray starting conditions
          Position init_r;
          Direction init_u;

          // Generate the starting position/direction of the ray
          if (orthographic_width_ == 0.0) { // perspective projection
            double this_phi =
              -horiz_fov_radians / 2.0 + dphi * horiz + 0.5 * dphi;
            double this_mu =
              -vert_fov_radians / 2.0 + dmu * vert + M_PI / 2.0 + 0.5 * dmu;
            Direction camera_local_vec;
            camera_local_vec.x = std::cos(this_phi) * std::sin(this_mu);
            camera_local_vec.y = std::sin(this_phi) * std::sin(this_mu);
            camera_local_vec.z = std::cos(this_mu);
            init_u = camera_local_vec.rotate(camera_to_model);
            init_r = camera_position_;
          } else { // orthographic projection
            init_u = looking_direction;

            double x_pix_coord = (static_cast<double>(horiz) - p0 / 2.0) / p0;
            double y_pix_coord = (static_cast<double>(vert) - p1 / 2.0) / p0;

            init_r = camera_position_;
            init_r += cam_yaxis * x_pix_coord * orthographic_width_;
            init_r += cam_zaxis * y_pix_coord * orthographic_width_;
          }

          // Resets internal geometry state of particle
          p.init_from_r_u(init_r, init_u);

          bool hitsomething = false;
          bool intersection_found = true;
          int loop_counter = 0;

          this_line_segments[tid][horiz].clear();

          int first_surface =
            -1; // surface first passed when entering the model
          bool first_inside_model = true; // false after entering the model
          while (intersection_found) {
            bool inside_cell = false;

            int32_t i_surface = std::abs(p.surface()) - 1;
            if (i_surface > 0 &&
                model::surfaces[i_surface]->geom_type_ == GeometryType::DAG) {
#ifdef DAGMC
              int32_t i_cell = next_cell(i_surface,
                p.cell_last(p.n_coord() - 1), p.lowest_coord().universe);
              inside_cell = i_cell >= 0;
#else
              fatal_error(
                "Not compiled for DAGMC, but somehow you have a DAGCell!");
#endif
            } else {
              inside_cell = exhaustive_find_cell(p);
            }

            if (inside_cell) {

              // This allows drawing wireframes with surface intersection
              // edges on the model boundary for the same cell.
              if (first_inside_model) {
                this_line_segments[tid][horiz].emplace_back(
                  color_by_ == PlotColorBy::mats ? p.material()
                                                 : p.lowest_coord().cell,
                  0.0, first_surface);
                first_inside_model = false;
              }

              hitsomething = true;
              intersection_found = true;
              auto dist = distance_to_boundary(p);
              this_line_segments[tid][horiz].emplace_back(
                color_by_ == PlotColorBy::mats ? p.material()
                                               : p.lowest_coord().cell,
                dist.distance, std::abs(dist.surface_index));

              // Advance particle
              for (int lev = 0; lev < p.n_coord(); ++lev) {
                p.coord(lev).r += dist.distance * p.coord(lev).u;
              }
              p.surface() = dist.surface_index;
              p.n_coord_last() = p.n_coord();
              p.n_coord() = dist.coord_level;
              if (dist.lattice_translation[0] != 0 ||
                  dist.lattice_translation[1] != 0 ||
                  dist.lattice_translation[2] != 0) {
                cross_lattice(p, dist);
              }

            } else {
              first_surface = advance_to_boundary_from_void(p);
              intersection_found =
                first_surface != -1; // -1 if no surface found
            }
            loop_counter++;
            if (loop_counter > MAX_INTERSECTIONS)
              fatal_error("Infinite loop in projection plot");
          }

          // Now color the pixel based on what we have intersected...
          // Loops backwards over intersections.
          Position current_color(
            not_found_.red, not_found_.green, not_found_.blue);
          const auto& segments = this_line_segments[tid][horiz];
          for (unsigned i = segments.size(); i-- > 0;) {
            int colormap_idx = segments[i].id;
            RGBColor seg_color = colors_[colormap_idx];
            Position seg_color_vec(
              seg_color.red, seg_color.green, seg_color.blue);
            double mixing = std::exp(-xs_[colormap_idx] * segments[i].length);
            current_color =
              current_color * mixing + (1.0 - mixing) * seg_color_vec;
            RGBColor result;
            result.red = static_cast<uint8_t>(current_color.x);
            result.green = static_cast<uint8_t>(current_color.y);
            result.blue = static_cast<uint8_t>(current_color.z);
            data(horiz, vert) = result;
          }

          // Check to draw wireframe in horizontal direction. No inter-thread
          // comm.
          if (horiz > 0) {
            if (!trackstack_equivalent(this_line_segments[tid][horiz],
                  this_line_segments[tid][horiz - 1])) {
              wireframe_initial(horiz, vert) = 1;
            }
          }
        }
      } // end "if" vert in correct range

      // We require a barrier before comparing vertical neighbors' intersection
      // stacks. i.e. all threads must be done with their line.
#pragma omp barrier

      // Now that the horizontal line has finished rendering, we can fill in
      // wireframe entries that require comparison among all the threads. Hence
      // the omp barrier being used. It has to be OUTSIDE any if blocks!
      if (vert < pixels_[1]) {
        // Loop over horizontal pixels, checking intersection stack of upper
        // neighbor

        const std::vector<std::vector<TrackSegment>>* top_cmp = nullptr;
        if (tid == 0)
          top_cmp = &old_segments;
        else
          top_cmp = &this_line_segments[tid - 1];

        for (int horiz = 0; horiz < pixels_[0]; ++horiz) {
          if (!trackstack_equivalent(
                this_line_segments[tid][horiz], (*top_cmp)[horiz])) {
            wireframe_initial(horiz, vert) = 1;
          }
        }
      }

      // We need another barrier to ensure threads don't proceed to modify their
      // intersection stacks on that horizontal line while others are
      // potentially still working on the above.
#pragma omp barrier
      vert += n_threads;
    }
  } // end omp parallel

  // Now thicken the wireframe lines and apply them to our image
  for (int vert = 0; vert < pixels_[1]; ++vert) {
    for (int horiz = 0; horiz < pixels_[0]; ++horiz) {
      if (wireframe_initial(horiz, vert)) {
        if (wireframe_thickness_ == 1)
          data(horiz, vert) = wireframe_color_;
        for (int i = -wireframe_thickness_ / 2; i < wireframe_thickness_ / 2;
             ++i)
          for (int j = -wireframe_thickness_ / 2; j < wireframe_thickness_ / 2;
               ++j)
            if (i * i + j * j < wireframe_thickness_ * wireframe_thickness_) {

              // Check if wireframe pixel is out of bounds
              int w_i = std::max(std::min(horiz + i, pixels_[0] - 1), 0);
              int w_j = std::max(std::min(vert + j, pixels_[1] - 1), 0);
              data(w_i, w_j) = wireframe_color_;
            }
      }
    }
  }

#ifdef USE_LIBPNG
  output_png(path_plot(), data);
#else
  output_ppm(path_plot(), data);
#endif
}

void ProjectionPlot::print_info() const
{
  fmt::print("Plot Type: Projection\n");
  fmt::print("Camera position: {} {} {}\n", camera_position_.x,
    camera_position_.y, camera_position_.z);
  fmt::print("Look at: {} {} {}\n", look_at_.x, look_at_.y, look_at_.z);
  fmt::print(
    "Horizontal field of view: {} degrees\n", horizontal_field_of_view_);
  fmt::print("Pixels: {} {}\n", pixels_[0], pixels_[1]);
}

void ProjectionPlot::set_opacities(pugi::xml_node node)
{
  xs_.resize(colors_.size(), 1e6); // set to large value for opaque by default

  for (auto cn : node.children("color")) {
    // Make sure 3 values are specified for RGB
    double user_xs = std::stod(get_node_value(cn, "xs"));
    int col_id = std::stoi(get_node_value(cn, "id"));

    // Add RGB
    if (PlotColorBy::cells == color_by_) {
      if (model::cell_map.find(col_id) != model::cell_map.end()) {
        col_id = model::cell_map[col_id];
        xs_[col_id] = user_xs;
      } else {
        warning(fmt::format(
          "Could not find cell {} specified in plot {}", col_id, id()));
      }
    } else if (PlotColorBy::mats == color_by_) {
      if (model::material_map.find(col_id) != model::material_map.end()) {
        col_id = model::material_map[col_id];
        xs_[col_id] = user_xs;
      } else {
        warning(fmt::format(
          "Could not find material {} specified in plot {}", col_id, id()));
      }
    }
  }
}

void ProjectionPlot::set_orthographic_width(pugi::xml_node node)
{
  if (check_for_node(node, "orthographic_width")) {
    double orthographic_width =
      std::stod(get_node_value(node, "orthographic_width", true));
    if (orthographic_width < 0.0)
      fatal_error("Requires positive orthographic_width");
    orthographic_width_ = orthographic_width;
  }
}

void ProjectionPlot::set_wireframe_thickness(pugi::xml_node node)
{
  if (check_for_node(node, "wireframe_thickness")) {
    int wireframe_thickness =
      std::stoi(get_node_value(node, "wireframe_thickness", true));
    if (wireframe_thickness < 0)
      fatal_error("Requires non-negative wireframe thickness");
    wireframe_thickness_ = wireframe_thickness;
  }
}

void ProjectionPlot::set_wireframe_ids(pugi::xml_node node)
{
  if (check_for_node(node, "wireframe_ids")) {
    wireframe_ids_ = get_node_array<int>(node, "wireframe_ids");
    // It is read in as actual ID values, but we have to convert to indices in
    // mat/cell array
    for (auto& x : wireframe_ids_)
      x = color_by_ == PlotColorBy::mats ? model::material_map[x]
                                         : model::cell_map[x];
  }
  // We make sure the list is sorted in order to later use
  // std::binary_search.
  std::sort(wireframe_ids_.begin(), wireframe_ids_.end());
}

void ProjectionPlot::set_pixels(pugi::xml_node node)
{
  vector<int> pxls = get_node_array<int>(node, "pixels");
  if (pxls.size() != 2)
    fatal_error(
      fmt::format("<pixels> must be length 2 in projection plot {}", id()));
  pixels_[0] = pxls[0];
  pixels_[1] = pxls[1];
}

void ProjectionPlot::set_camera_position(pugi::xml_node node)
{
  vector<double> camera_pos = get_node_array<double>(node, "camera_position");
  if (camera_pos.size() != 3) {
    fatal_error(
      fmt::format("look_at element must have three floating point values"));
  }
  camera_position_.x = camera_pos[0];
  camera_position_.y = camera_pos[1];
  camera_position_.z = camera_pos[2];
}

void ProjectionPlot::set_look_at(pugi::xml_node node)
{
  vector<double> look_at = get_node_array<double>(node, "look_at");
  if (look_at.size() != 3) {
    fatal_error("look_at element must have three floating point values");
  }
  look_at_.x = look_at[0];
  look_at_.y = look_at[1];
  look_at_.z = look_at[2];
}

void ProjectionPlot::set_field_of_view(pugi::xml_node node)
{
  // Defaults to 70 degree horizontal field of view (see .h file)
  if (check_for_node(node, "field_of_view")) {
    double fov = std::stod(get_node_value(node, "field_of_view", true));
    if (fov < 180.0 && fov > 0.0) {
      horizontal_field_of_view_ = fov;
    } else {
      fatal_error(fmt::format(
        "Field of view for plot {} out-of-range. Must be in (0, 180).", id()));
    }
  }
}

extern "C" int openmc_id_map(const void* plot, int32_t* data_out)
{

  auto plt = reinterpret_cast<const SlicePlotBase*>(plot);
  if (!plt) {
    set_errmsg("Invalid slice pointer passed to openmc_id_map");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  if (plt->slice_color_overlaps_ && model::overlap_check_count.size() == 0) {
    model::overlap_check_count.resize(model::cells.size());
  }

  auto ids = plt->get_map<IdData>();

  // write id data to array
  std::copy(ids.data_.begin(), ids.data_.end(), data_out);

  return 0;
}

extern "C" int openmc_property_map(const void* plot, double* data_out)
{

  auto plt = reinterpret_cast<const SlicePlotBase*>(plot);
  if (!plt) {
    set_errmsg("Invalid slice pointer passed to openmc_id_map");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  if (plt->slice_color_overlaps_ && model::overlap_check_count.size() == 0) {
    model::overlap_check_count.resize(model::cells.size());
  }

  auto props = plt->get_map<PropertyData>();

  // write id data to array
  std::copy(props.data_.begin(), props.data_.end(), data_out);

  return 0;
}

} // namespace openmc
