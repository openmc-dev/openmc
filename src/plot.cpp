#include "openmc/plot.h"

#include <algorithm>
#include <fstream>
#include <sstream>

#include <fmt/core.h>
#include <fmt/ostream.h>
#include "xtensor/xview.hpp"

#include "openmc/constants.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
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

IdData::IdData(size_t h_res, size_t v_res)
  : data_({v_res, h_res, 2}, NOT_FOUND)
{ }

void
IdData::set_value(size_t y, size_t x, const Particle& p, int level) {
  Cell* c = model::cells[p.coord_[level].cell].get();
  data_(y,x,0) = c->id_;
  if (p.material_ == MATERIAL_VOID) {
    data_(y,x,1) = MATERIAL_VOID;
    return;
  } else if (c->type_ != Fill::UNIVERSE) {
    Material* m = model::materials[p.material_].get();
    data_(y,x,1) = m->id_;
  }
}

void IdData::set_overlap(size_t y, size_t x) {
  xt::view(data_, y, x, xt::all()) = OVERLAP;
}

PropertyData::PropertyData(size_t h_res, size_t v_res)
  : data_({v_res, h_res, 2}, NOT_FOUND)
{ }

void
PropertyData::set_value(size_t y, size_t x, const Particle& p, int level) {
  Cell* c = model::cells[p.coord_[level].cell].get();
  data_(y,x,0) = (p.sqrtkT_ * p.sqrtkT_) / K_BOLTZMANN;
  if (c->type_ != Fill::UNIVERSE && p.material_ != MATERIAL_VOID) {
    Material* m = model::materials[p.material_].get();
    data_(y,x,1) = m->density_gpcc_;
  }
}

void PropertyData::set_overlap(size_t y, size_t x) {
  data_(y, x) = OVERLAP;
}

//==============================================================================
// Global variables
//==============================================================================

namespace model {

std::vector<Plot> plots;
std::unordered_map<int, int> plot_map;
uint64_t plotter_seed = 1;

} // namespace model

//==============================================================================
// RUN_PLOT controls the logic for making one or many plots
//==============================================================================

extern "C"
int openmc_plot_geometry()
{
  for (auto pl : model::plots) {
    write_message(fmt::format("Processing plot {}: {}...",
      pl.id_, pl.path_plot_), 5);

    if (PlotType::slice == pl.type_) {
      // create 2D image
      create_ppm(pl);
    } else if (PlotType::voxel == pl.type_) {
      // create voxel file for 3D viewing
      create_voxel(pl);
    }
  }
  return 0;
}


void read_plots_xml()
{
  // Check if plots.xml exists
  std::string filename = settings::path_input + "plots.xml";
  if (!file_exists(filename)) {
    fatal_error("Plots XML file '" + filename + "' does not exist!");
  }

  write_message("Reading plot XML file...", 5);

  // Parse plots.xml file
  pugi::xml_document doc;
  doc.load_file(filename.c_str());

  pugi::xml_node root = doc.document_element();
  for (auto node : root.children("plot")) {
    Plot pl(node);
    model::plots.push_back(pl);
    model::plot_map[pl.id_] = model::plots.size() - 1;
  }
}

//==============================================================================
// CREATE_PPM creates an image based on user input from a plots.xml <plot>
// specification in the portable pixmap format (PPM)
//==============================================================================

void create_ppm(Plot pl)
{

  size_t width = pl.pixels_[0];
  size_t height = pl.pixels_[1];

  ImageData data({width, height}, pl.not_found_);

  // generate ids for the plot
  auto ids = pl.get_map<IdData>();

  // assign colors
  for (size_t y = 0; y < height; y++) {
    for (size_t x = 0; x < width; x++) {
      auto id = ids.data_(y, x, pl.color_by_);
      // no setting needed if not found
      if (id == NOT_FOUND) { continue; }
      if (id == OVERLAP) {
        data(x,y) = pl.overlap_color_;
        continue;
      }
      if (PlotColorBy::cells == pl.color_by_) {
        data(x,y) = pl.colors_[model::cell_map[id]];
      } else if (PlotColorBy::mats == pl.color_by_) {
        if (id == MATERIAL_VOID) {
          data(x,y) = WHITE;
          continue;
        }
        data(x,y) = pl.colors_[model::material_map[id]];
      } // color_by if-else
    } // x for loop
  } // y for loop

  // draw mesh lines if present
  if (pl.index_meshlines_mesh_ >= 0) {draw_mesh_lines(pl, data);}

  // write ppm data to file
  output_ppm(pl, data);
}

void
Plot::set_id(pugi::xml_node plot_node)
{
  // Copy data into plots
  if (check_for_node(plot_node, "id")) {
    id_ = std::stoi(get_node_value(plot_node, "id"));
  } else {
    fatal_error("Must specify plot id in plots XML file.");
  }

  // Check to make sure 'id' hasn't been used
  if (model::plot_map.find(id_) != model::plot_map.end()) {
    fatal_error(fmt::format("Two or more plots use the same unique ID: {}", id_));
  }
}

void
Plot::set_type(pugi::xml_node plot_node)
{
  // Copy plot type
  // Default is slice
  type_ = PlotType::slice;
  // check type specified on plot node
  if (check_for_node(plot_node, "type")) {
    std::string type_str = get_node_value(plot_node, "type", true);
    // set type using node value
    if (type_str == "slice") {
      type_ = PlotType::slice;
    }
    else if (type_str == "voxel") {
      type_ = PlotType::voxel;
    } else {
      // if we're here, something is wrong
      fatal_error(fmt::format("Unsupported plot type '{}' in plot {}",
        type_str, id_));
    }
  }
}

void
Plot::set_output_path(pugi::xml_node plot_node)
{
  // Set output file path
  std::string filename;

  if (check_for_node(plot_node, "filename")) {
    filename = get_node_value(plot_node, "filename");
  } else {
    filename = fmt::format("plot_{}", id_);
  }
  // add appropriate file extension to name
  switch(type_) {
  case PlotType::slice:
    filename.append(".ppm");
    break;
  case PlotType::voxel:
    filename.append(".h5");
    break;
  }

  path_plot_ = filename;

  // Copy plot pixel size
  std::vector<int> pxls = get_node_array<int>(plot_node, "pixels");
  if (PlotType::slice == type_) {
    if (pxls.size() == 2) {
      pixels_[0] = pxls[0];
      pixels_[1] = pxls[1];
    } else {
      fatal_error(fmt::format("<pixels> must be length 2 in slice plot {}", id_));
    }
  } else if (PlotType::voxel == type_) {
    if (pxls.size() == 3) {
      pixels_[0] = pxls[0];
      pixels_[1] = pxls[1];
      pixels_[2] = pxls[2];
    } else {
      fatal_error(fmt::format("<pixels> must be length 3 in voxel plot {}", id_));
    }
  }
}

void
Plot::set_bg_color(pugi::xml_node plot_node)
{
  // Copy plot background color
  if (check_for_node(plot_node, "background")) {
    std::vector<int> bg_rgb = get_node_array<int>(plot_node, "background");
    if (PlotType::voxel == type_) {
      if (mpi::master) {
        warning(fmt::format("Background color ignored in voxel plot {}", id_));
      }
    }
    if (bg_rgb.size() == 3) {
      not_found_ = bg_rgb;
    } else {
      fatal_error(fmt::format("Bad background RGB in plot {}", id_));
    }
  }
}

void
Plot::set_basis(pugi::xml_node plot_node)
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
      fatal_error(fmt::format("Unsupported plot basis '{}' in plot {}",
        pl_basis, id_));
    }
  }
}

void
Plot::set_origin(pugi::xml_node plot_node)
{
  // Copy plotting origin
  auto pl_origin = get_node_array<double>(plot_node, "origin");
  if (pl_origin.size() == 3) {
    origin_ = pl_origin;
  } else {
    fatal_error(fmt::format("Origin must be length 3 in plot {}", id_));
  }
}

void
Plot::set_width(pugi::xml_node plot_node)
{
  // Copy plotting width
  std::vector<double> pl_width = get_node_array<double>(plot_node, "width");
  if (PlotType::slice == type_) {
    if (pl_width.size() == 2) {
      width_.x = pl_width[0];
      width_.y = pl_width[1];
    } else {
      fatal_error(fmt::format("<width> must be length 2 in slice plot {}", id_));
    }
  } else if (PlotType::voxel == type_) {
    if (pl_width.size() == 3) {
      pl_width = get_node_array<double>(plot_node, "width");
      width_ = pl_width;
    } else {
      fatal_error(fmt::format("<width> must be length 3 in voxel plot {}", id_));
    }
  }
}

void
Plot::set_universe(pugi::xml_node plot_node)
{
  // Copy plot universe level
  if (check_for_node(plot_node, "level")) {
    level_ = std::stoi(get_node_value(plot_node, "level"));
    if (level_ < 0) {
      fatal_error(fmt::format("Bad universe level in plot {}", id_));
    }
  } else {
    level_ = PLOT_LEVEL_LOWEST;
  }
}

void
Plot::set_default_colors(pugi::xml_node plot_node)
{
  // Copy plot color type and initialize all colors randomly
  std::string pl_color_by = "cell";
  if (check_for_node(plot_node, "color_by")) {
    pl_color_by = get_node_value(plot_node, "color_by", true);
  }
  if ("cell" == pl_color_by) {
    color_by_ = PlotColorBy::cells;
    colors_.resize(model::cells.size());
  } else if("material" == pl_color_by) {
    color_by_ = PlotColorBy::mats;
    colors_.resize(model::materials.size());
  } else {
    fatal_error(fmt::format("Unsupported plot color type '{}' in plot {}",
      pl_color_by, id_));
  }

  for (auto& c : colors_) {
    c = random_color();
    // make sure we don't interfere with some default colors
    while (c == RED || c == WHITE) {
      c = random_color();
    }
  }
}

void
Plot::set_user_colors(pugi::xml_node plot_node)
{
  if (!plot_node.select_nodes("color").empty() && PlotType::voxel == type_) {
    if (mpi::master) {
      warning(fmt::format("Color specifications ignored in voxel plot {}", id_));
    }
  }

  for (auto cn : plot_node.children("color")) {
    // Make sure 3 values are specified for RGB
    std::vector<int> user_rgb = get_node_array<int>(cn, "rgb");
    if (user_rgb.size() != 3) {
      fatal_error(fmt::format("Bad RGB in plot {}", id_));
    }
    // Ensure that there is an id for this color specification
    int col_id;
    if (check_for_node(cn, "id")) {
      col_id = std::stoi(get_node_value(cn, "id"));
    } else {
      fatal_error(fmt::format(
        "Must specify id for color specification in plot {}", id_));
    }
    // Add RGB
    if (PlotColorBy::cells == color_by_) {
      if (model::cell_map.find(col_id) != model::cell_map.end()) {
        col_id = model::cell_map[col_id];
        colors_[col_id] = user_rgb;
      } else {
        fatal_error(fmt::format("Could not find cell {} specified in plot {}",
          col_id, id_));
      }
    } else if (PlotColorBy::mats == color_by_) {
      if (model::material_map.find(col_id) != model::material_map.end()) {
        col_id = model::material_map[col_id];
        colors_[col_id] = user_rgb;
      } else {
        fatal_error(fmt::format(
          "Could not find material {} specified in plot {}", col_id, id_));
      }
    }
  } // color node loop
}

void
Plot::set_meshlines(pugi::xml_node plot_node)
{
  // Deal with meshlines
  pugi::xpath_node_set mesh_line_nodes = plot_node.select_nodes("meshlines");

  if (!mesh_line_nodes.empty()) {
    if (PlotType::voxel == type_) {
      warning(fmt::format("Meshlines ignored in voxel plot {}", id_));
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
          "Must specify a meshtype for meshlines specification in plot {}", id_));
      }

      // Ensure that there is a linewidth for this meshlines specification
      std::string meshline_width;
      if (check_for_node(meshlines_node, "linewidth")) {
        meshline_width = get_node_value(meshlines_node, "linewidth");
        meshlines_width_ = std::stoi(meshline_width);
      } else {
        fatal_error(fmt::format(
          "Must specify a linewidth for meshlines specification in plot {}", id_));
      }

      // Check for color
      if (check_for_node(meshlines_node, "color")) {
        // Check and make sure 3 values are specified for RGB
        std::vector<int> ml_rgb = get_node_array<int>(meshlines_node, "color");
        if (ml_rgb.size() != 3) {
          fatal_error(fmt::format("Bad RGB for meshlines color in plot {}", id_));
        }
        meshlines_color_ = ml_rgb;
      }

      // Set mesh based on type
      if ("ufs" == meshtype) {
        if (!simulation::ufs_mesh) {
          fatal_error(fmt::format("No UFS mesh for meshlines on plot {}", id_));
        } else {
          for (int i = 0; i < model::meshes.size(); ++i) {
            if (const auto* m
                = dynamic_cast<const RegularMesh*>(model::meshes[i].get())) {
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
          fatal_error(fmt::format("No entropy mesh for meshlines on plot {}", id_));
        } else {
          for (int i = 0; i < model::meshes.size(); ++i) {
            if (const auto* m
                = dynamic_cast<const RegularMesh*>(model::meshes[i].get())) {
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
            "mesh specification in plot {}", id_));
        }
        // find the tally index
        int idx;
        int err = openmc_get_mesh_index(tally_mesh_id, &idx);
        if (err != 0) {
          fatal_error(fmt::format("Could not find mesh {} specified in "
            "meshlines for plot {}", tally_mesh_id, id_));
        }
        index_meshlines_mesh_ = idx;
      } else {
        fatal_error(fmt::format("Invalid type for meshlines on plot {}", id_ ));
      }
    } else {
      fatal_error(fmt::format("Mutliple meshlines specified in plot {}", id_));
    }
  }
}

void
Plot::set_mask(pugi::xml_node plot_node)
{
  // Deal with masks
  pugi::xpath_node_set mask_nodes = plot_node.select_nodes("mask");

  if (!mask_nodes.empty()) {
    if (PlotType::voxel == type_) {
      if (mpi::master) {
        warning(fmt::format("Mask ignored in voxel plot {}", id_));
      }
    }

    if (mask_nodes.size() == 1) {
      // Get pointer to mask
      pugi::xml_node mask_node = mask_nodes[0].node();

      // Determine how many components there are and allocate
      std::vector<int> iarray = get_node_array<int>(mask_node, "components");
      if (iarray.size() == 0) {
        fatal_error(fmt::format("Missing <components> in mask of plot {}", id_));
      }

      // First we need to change the user-specified identifiers to indices
      // in the cell and material arrays
      for (auto& col_id : iarray) {
        if (PlotColorBy::cells == color_by_) {
          if (model::cell_map.find(col_id) != model::cell_map.end()) {
            col_id  = model::cell_map[col_id];
          }
          else {
            fatal_error(fmt::format("Could not find cell {} specified in the "
              "mask in plot {}", col_id, id_));
          }
        } else if (PlotColorBy::mats == color_by_) {
          if (model::material_map.find(col_id) != model::material_map.end()) {
            col_id = model::material_map[col_id];
          }
          else {
            fatal_error(fmt::format("Could not find material {} specified in "
              "the mask in plot {}", col_id, id_));
          }
        }
      }

      // Alter colors based on mask information
      for (int j = 0; j < colors_.size(); j++) {
        if (std::find(iarray.begin(), iarray.end(), j) == iarray.end()) {
          if (check_for_node(mask_node, "background")) {
            std::vector<int> bg_rgb = get_node_array<int>(mask_node, "background");
            colors_[j] = bg_rgb;
          } else {
            colors_[j] = WHITE;
          }
        }
      }

    } else {
      fatal_error(fmt::format("Mutliple masks specified in plot {}", id_));
    }
  }
}

void Plot::set_overlap_color(pugi::xml_node plot_node) {
  color_overlaps_ = false;
  if (check_for_node(plot_node, "show_overlaps")) {
    color_overlaps_ = get_node_value_bool(plot_node, "show_overlaps");
    // check for custom overlap color
    if (check_for_node(plot_node, "overlap_color")) {
      if (!color_overlaps_) {
        warning(fmt::format(
          "Overlap color specified in plot {} but overlaps won't be shown.", id_));
      }
      std::vector<int> olap_clr = get_node_array<int>(plot_node, "overlap_color");
      if (olap_clr.size() == 3) {
        overlap_color_ = olap_clr;
      } else {
        fatal_error(fmt::format("Bad overlap RGB in plot {}", id_));
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

Plot::Plot(pugi::xml_node plot_node)
  : index_meshlines_mesh_{-1}, overlap_color_{RED}
{
  set_id(plot_node);
  set_type(plot_node);
  set_output_path(plot_node);
  set_bg_color(plot_node);
  set_basis(plot_node);
  set_origin(plot_node);
  set_width(plot_node);
  set_universe(plot_node);
  set_default_colors(plot_node);
  set_user_colors(plot_node);
  set_meshlines(plot_node);
  set_mask(plot_node);
  set_overlap_color(plot_node);
} // End Plot constructor

//==============================================================================
// OUTPUT_PPM writes out a previously generated image to a PPM file
//==============================================================================

void output_ppm(Plot pl, const ImageData& data)
{
  // Open PPM file for writing
  std::string fname = pl.path_plot_;
  fname = strtrim(fname);
  std::ofstream of;

  of.open(fname);

  // Write header
  of << "P6\n";
  of << pl.pixels_[0] << " " << pl.pixels_[1] << "\n";
  of << "255\n";
  of.close();

  of.open(fname, std::ios::binary | std::ios::app);
  // Write color for each pixel
  for (int y = 0; y < pl.pixels_[1]; y++) {
    for (int x = 0; x < pl.pixels_[0]; x++) {
      RGBColor rgb = data(x,y);
      of << rgb.red << rgb.green << rgb.blue;
    }
  }
  of << "\n";
}

//==============================================================================
// DRAW_MESH_LINES draws mesh line boundaries on an image
//==============================================================================

void draw_mesh_lines(Plot pl, ImageData& data)
{
  RGBColor rgb;
  rgb = pl.meshlines_color_;

  int ax1, ax2;
  switch(pl.basis_) {
  case PlotBasis::xy :
    ax1 = 0;
    ax2 = 1;
    break;
  case PlotBasis::xz :
    ax1 = 0;
    ax2 = 2;
    break;
  case PlotBasis::yz :
    ax1 = 1;
    ax2 = 2;
    break;
  default:
    UNREACHABLE();
  }

  Position ll_plot {pl.origin_};
  Position ur_plot {pl.origin_};

  ll_plot[ax1] -= pl.width_[0] / 2.;
  ll_plot[ax2] -= pl.width_[1] / 2.;
  ur_plot[ax1] += pl.width_[0] / 2.;
  ur_plot[ax2] += pl.width_[1] / 2.;

  Position width = ur_plot - ll_plot;

  // Find the (axis-aligned) lines of the mesh that intersect this plot.
  auto axis_lines = model::meshes[pl.index_meshlines_mesh_]
    ->plot(ll_plot, ur_plot);

  // Find the bounds along the second axis (accounting for low-D meshes).
  int ax2_min, ax2_max;
  if (axis_lines.second.size() > 0) {
    double frac = (axis_lines.second.back() - ll_plot[ax2]) / width[ax2];
    ax2_min = (1.0 - frac) * pl.pixels_[1];
    if (ax2_min < 0) ax2_min = 0;
    frac = (axis_lines.second.front() - ll_plot[ax2]) / width[ax2];
    ax2_max = (1.0 - frac) * pl.pixels_[1];
    if (ax2_max > pl.pixels_[1]) ax2_max = pl.pixels_[1];
  } else {
    ax2_min = 0;
    ax2_max = pl.pixels_[1];
  }

  // Iterate across the first axis and draw lines.
  for (auto ax1_val : axis_lines.first) {
    double frac = (ax1_val - ll_plot[ax1]) / width[ax1];
    int ax1_ind = frac * pl.pixels_[0];
    for (int ax2_ind = ax2_min; ax2_ind < ax2_max; ++ax2_ind) {
      for (int plus = 0; plus <= pl.meshlines_width_; plus++) {
        if (ax1_ind+plus >= 0 && ax1_ind+plus < pl.pixels_[0])
          data(ax1_ind+plus, ax2_ind) = rgb;
        if (ax1_ind-plus >= 0 && ax1_ind-plus < pl.pixels_[0])
          data(ax1_ind-plus, ax2_ind) = rgb;
      }
    }
  }

  // Find the bounds along the first axis.
  int ax1_min, ax1_max;
  if (axis_lines.first.size() > 0) {
    double frac = (axis_lines.first.front() - ll_plot[ax1]) / width[ax1];
    ax1_min = frac * pl.pixels_[0];
    if (ax1_min < 0) ax1_min = 0;
    frac = (axis_lines.first.back() - ll_plot[ax1]) / width[ax1];
    ax1_max = frac * pl.pixels_[0];
    if (ax1_max > pl.pixels_[0]) ax1_max = pl.pixels_[0];
  } else {
    ax1_min = 0;
    ax1_max = pl.pixels_[0];
  }

  // Iterate across the second axis and draw lines.
  for (auto ax2_val : axis_lines.second) {
    double frac = (ax2_val - ll_plot[ax2]) / width[ax2];
    int ax2_ind = (1.0 - frac) * pl.pixels_[1];
    for (int ax1_ind = ax1_min; ax1_ind < ax1_max; ++ax1_ind) {
      for (int plus = 0; plus <= pl.meshlines_width_; plus++) {
        if (ax2_ind+plus >= 0 && ax2_ind+plus < pl.pixels_[1])
          data(ax1_ind, ax2_ind+plus) = rgb;
        if (ax2_ind-plus >= 0 && ax2_ind-plus < pl.pixels_[1])
          data(ax1_ind, ax2_ind-plus) = rgb;
      }
    }
  }
}

//==============================================================================
// CREATE_VOXEL outputs a binary file that can be input into silomesh for 3D
// geometry visualization.  It works the same way as create_ppm by dragging a
// particle across the geometry for the specified number of voxels. The first 3
// int's in the binary are the number of x, y, and z voxels.  The next 3
// double's are the widths of the voxels in the x, y, and z directions. The
// next 3 double's are the x, y, and z coordinates of the lower left
// point. Finally the binary is filled with entries of four int's each. Each
// 'row' in the binary contains four int's: 3 for x,y,z position and 1 for
// cell or material id.  For 1 million voxels this produces a file of
// approximately 15MB.
// =============================================================================

void create_voxel(Plot pl)
{
  // compute voxel widths in each direction
  std::array<double, 3> vox;
  vox[0] = pl.width_[0]/(double)pl.pixels_[0];
  vox[1] = pl.width_[1]/(double)pl.pixels_[1];
  vox[2] = pl.width_[2]/(double)pl.pixels_[2];

  // initial particle position
  Position ll = pl.origin_ - pl.width_ / 2.;

  // Open binary plot file for writing
  std::ofstream of;
  std::string fname = std::string(pl.path_plot_);
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
  std::array<int, 3> pixels;
  std::copy(pl.pixels_.begin(), pl.pixels_.end(), pixels.begin());
  write_attribute(file_id, "num_voxels", pixels);
  write_attribute(file_id, "voxel_width", vox);
  write_attribute(file_id, "lower_left", ll);

  // Create dataset for voxel data -- note that the dimensions are reversed
  // since we want the order in the file to be z, y, x
  hsize_t dims[3];
  dims[0] = pl.pixels_[2];
  dims[1] = pl.pixels_[1];
  dims[2] = pl.pixels_[0];
  hid_t dspace, dset, memspace;
  voxel_init(file_id, &(dims[0]), &dspace, &dset, &memspace);

  PlotBase pltbase;
  pltbase.width_ = pl.width_;
  pltbase.origin_ = pl.origin_;
  pltbase.basis_ = PlotBasis::xy;
  pltbase.pixels_ = pl.pixels_;
  pltbase.level_ = -1; // all universes for voxel files
  pltbase.color_overlaps_ = pl.color_overlaps_;

  ProgressBar pb;
  for (int z = 0; z < pl.pixels_[2]; z++) {
    // update progress bar
    pb.set_value(100.*(double)z/(double)(pl.pixels_[2]-1));

    // update z coordinate
    pltbase.origin_.z = ll.z + z * vox[2];

    // generate ids using plotbase
    IdData ids = pltbase.get_map<IdData>();

    // select only cell/material ID data and flip the y-axis
    int idx = pl.color_by_ == PlotColorBy::cells ? 0 : 1;
    xt::xtensor<int32_t, 2> data_slice = xt::view(ids.data_, xt::all(), xt::all(), idx);
    xt::xtensor<int32_t, 2> data_flipped = xt::flip(data_slice, 0);

    // Write to HDF5 dataset
    voxel_write_slice(z, dspace, dset, memspace, data_flipped.data());
  }

  voxel_finalize(dspace, dset, memspace);
  file_close(file_id);
}

void
voxel_init(hid_t file_id, const hsize_t* dims,
           hid_t* dspace, hid_t* dset, hid_t* memspace)
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
  hssize_t offset[3] {x, 0, 0};
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

RGBColor random_color(void) {
  return {int(prn(&model::plotter_seed)*255),
          int(prn(&model::plotter_seed)*255),
          int(prn(&model::plotter_seed)*255)};
}

extern "C" int openmc_id_map(const void* plot, int32_t* data_out)
{

  auto plt = reinterpret_cast<const PlotBase*>(plot);
  if (!plt) {
    set_errmsg("Invalid slice pointer passed to openmc_id_map");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  if (plt->color_overlaps_ && model::overlap_check_count.size() == 0) {
    model::overlap_check_count.resize(model::cells.size());
  }

  auto ids = plt->get_map<IdData>();

  // write id data to array
  std::copy(ids.data_.begin(), ids.data_.end(), data_out);

  return 0;
}

extern "C" int openmc_property_map(const void* plot, double* data_out) {

  auto plt = reinterpret_cast<const PlotBase*>(plot);
  if (!plt) {
    set_errmsg("Invalid slice pointer passed to openmc_id_map");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  if (plt->color_overlaps_ && model::overlap_check_count.size() == 0) {
    model::overlap_check_count.resize(model::cells.size());
  }

  auto props = plt->get_map<PropertyData>();

  // write id data to array
  std::copy(props.data_.begin(), props.data_.end(), data_out);

  return 0;
}


} // namespace openmc
