#include <fstream>
#include <sstream>

#include "openmc/plot.h"
#include "openmc/constants.h"
#include "openmc/settings.h"
#include "openmc/error.h"
#include "openmc/particle.h"
#include "openmc/geometry.h"
#include "openmc/cell.h"
#include "openmc/material.h"
#include "openmc/string_functions.h"
#include "openmc/mesh.h"
#include "openmc/output.h"
#include "openmc/hdf5_interface.h"
#include "openmc/random_lcg.h"
#include "openmc/output.h"
#include "openmc/progress_bar.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================
  
int PLOT_LEVEL_LOWEST = -1;

std::unordered_map<int, int> plot_map;

int n_plots;

std::vector<Plot> plots;

const RGBColor WHITE = {255, 255, 255};

//==============================================================================
// RUN_PLOT controls the logic for making one or many plots
//==============================================================================

extern "C"
int openmc_plot_geometry()
{
  int err;

  for (auto pl : plots) {
    std::stringstream ss;
    ss << "Processing plot " << pl.id_ << ": "
       << pl.path_plot_ << "...";
    write_message(ss.str(), 5);

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


void
read_plots(pugi::xml_node* plots_node)
{
  n_plots = 0;
  for (auto node : plots_node->children("plot")) {
    Plot pl(node);
    plots.push_back(pl);
    plot_map[pl.id_] = n_plots++;
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

  double in_pixel = (pl.width_[0])/static_cast<double>(width);
  double out_pixel = (pl.width_[1])/static_cast<double>(height);

  ImageData data;
  data.resize({width, height});

  int in_i, out_i;
  double xyz[3];
  switch(pl.basis_) {
  case PlotBasis::xy :
    in_i = 0;
    out_i = 1;
    xyz[0] = pl.origin_[0] - pl.width_[0] / 2.;
    xyz[1] = pl.origin_[1] + pl.width_[1] / 2.;
    xyz[2] = pl.origin_[2];
    break;
  case PlotBasis::xz :
    in_i = 0;
    out_i = 2;
    xyz[0] = pl.origin_[0] - pl.width_[0] / 2.;
    xyz[1] = pl.origin_[1];
    xyz[2] = pl.origin_[2] + pl.width_[1] / 2.;
    break;
  case PlotBasis::yz :
    in_i = 1;
    out_i = 2;
    xyz[0] = pl.origin_[0];
    xyz[1] = pl.origin_[1] - pl.width_[0] / 2.;
    xyz[2] = pl.origin_[2] + pl.width_[1] / 2.;
    break;
  }

  double dir[3] = {0.5, 0.5, 0.5};

#pragma omp parallel
{
  Particle p;
  p.initialize();
  std::copy(xyz, xyz+3, p.coord[0].xyz);
  std::copy(dir, dir+3, p.coord[0].uvw);
  p.coord[0].universe = openmc_root_universe;

#pragma omp for
  for (int y = 0; y < height; y++) {
    p.coord[0].xyz[out_i] = xyz[out_i] - out_pixel * y;
    for (int x = 0; x < width; x++) {
      // local variables
      RGBColor rgb;
      int id;
      p.coord[0].xyz[in_i] = xyz[in_i] + in_pixel * x;
      position_rgb(p, pl, rgb, id);
      data(x,y) = rgb;
    }
  }
}
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
  if (plot_map.find(id_) != plot_map.end()) {
    std::stringstream err_msg;
    err_msg << "Two or more plots use the same unique ID: " << id_;
    fatal_error(err_msg.str());
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
    std::stringstream err_msg;
    err_msg << "Unsupported plot type '" << type_str
            << "' in plot " << id_;
    fatal_error(err_msg.str());
    }
  }
}

void
Plot::set_output_path(pugi::xml_node plot_node)
{
  // Set output file path
  std::stringstream filename;
  filename << "plot_" << id_;

  if (check_for_node(plot_node, "filename")) {
    filename << get_node_value(plot_node, "filename");
  }
  // add appropriate file extension to name
  switch(type_) {
  case PlotType::slice:
    filename << ".ppm";
    break;
  case PlotType::voxel:
    filename << ".h5";
    break;
  }

  path_plot_ = filename.str();

  // Copy plot pixel size
  std::vector<int> pxls = get_node_array<int>(plot_node, "pixels");
  if (PlotType::slice == type_) {
    if (pxls.size() == 2) {
      pixels_[0] = pxls[0];
      pixels_[1] = pxls[1];
    } else {
      std::stringstream err_msg;
      err_msg << "<pixels> must be length 2 in slice plot "
              << id_;
      fatal_error(err_msg.str());
    }
  } else if (PlotType::voxel == type_) {
    if (pxls.size() == 3) {
      pixels_[0] = pxls[0];
      pixels_[1] = pxls[1];
      pixels_[2] = pxls[2];
    } else {
      std::stringstream err_msg;
      err_msg << "<pixels> must be length 3 in voxel plot "
              << id_;
      fatal_error(err_msg.str());
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
      if (openmc_master) {
        std::stringstream err_msg;
        err_msg << "Background color ignored in voxel plot "
                << id_;
        warning(err_msg.str());
      }
    }
    if (bg_rgb.size() == 3) {
      not_found_ = bg_rgb;
    } else {
      std::stringstream err_msg;
      err_msg << "Bad background RGB in plot "
              << id_;
      fatal_error(err_msg);
    }
  } else {
    // default to a white background
    not_found_ = WHITE;
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
      std::stringstream err_msg;
      err_msg << "Unsupported plot basis '" << pl_basis
              << "' in plot " << id_;
      fatal_error(err_msg);
    }
  }
}

void
Plot::set_origin(pugi::xml_node plot_node)
{
  // Copy plotting origin
  std::vector<double> pl_origin = get_node_array<double>(plot_node, "origin");
  if (pl_origin.size() == 3) {
    origin_[0] = pl_origin[0];
    origin_[1] = pl_origin[1];
    origin_[2] = pl_origin[2];
  } else {
    std::stringstream err_msg;
    err_msg << "Origin must be length 3 in plot "
            << id_;
    fatal_error(err_msg);
  }
}

void
Plot::set_width(pugi::xml_node plot_node)
{
  // Copy plotting width
  std::vector<double> pl_width = get_node_array<double>(plot_node, "width");
  if (PlotType::slice == type_) {
    if (pl_width.size() == 2) {
      width_[0] = pl_width[0];
      width_[1] = pl_width[1];
    } else {
      std::stringstream err_msg;
      err_msg << "<width> must be length 2 in slice plot "
              << id_;
      fatal_error(err_msg);
    }
  } else if (PlotType::voxel == type_) {
    if (pl_width.size() == 3) {
      pl_width = get_node_array<double>(plot_node, "width");
      width_[0] = pl_width[0];
      width_[1] = pl_width[1];
      width_[2] = pl_width[2];
    } else {
      std::stringstream err_msg;
      err_msg << "<width> must be length 3 in voxel plot "
              << id_;
      fatal_error(err_msg);
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
      std::stringstream err_msg;
      err_msg << "Bad universe level in plot " << id_;
      fatal_error(err_msg);
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
    colors_.resize(n_cells);
    for (int i = 0; i < n_cells; i++) {
      colors_[i] = random_color();
    }

  } else if("material" == pl_color_by) {
    color_by_ = PlotColorBy::mats;
    colors_.resize(n_materials);
    for (int i = 0; i < materials.size(); i++) {
      colors_[i] = random_color();
    }
  } else {
    std::stringstream err_msg;
    err_msg << "Unsupported plot color type '" << pl_color_by
            << "' in plot " << id_;
    fatal_error(err_msg);
  }
}

void
Plot::set_user_colors(pugi::xml_node plot_node)
{
  if (!plot_node.select_nodes("color").empty() && PlotType::voxel == type_) {
    if (openmc_master) {
      std::stringstream err_msg;
      err_msg << "Color specifications ignored in voxel plot "
              << id_;
      warning(err_msg);
    }
  }
  
  for (auto cn : plot_node.children("color")) {
    // Make sure 3 values are specified for RGB
    std::vector<int> user_rgb = get_node_array<int>(cn, "rgb");
    if (user_rgb.size() != 3) {
      std::stringstream err_msg;
      err_msg << "Bad RGB in plot " << id_;
      fatal_error(err_msg);
    }
    // Ensure that there is an id for this color specification
    int col_id;
    if (check_for_node(cn, "id")) {
      col_id = std::stoi(get_node_value(cn, "id"));
    } else {
      std::stringstream err_msg;
      err_msg << "Must specify id for color specification in plot "
              << id_;
      fatal_error(err_msg);
    }
    // Add RGB
    if (PlotColorBy::cells == color_by_) {
      if (cell_map.find(col_id) != cell_map.end()) {
        col_id = cell_map[col_id];
        colors_[col_id] = user_rgb;
      } else {
        std::stringstream err_msg;
        err_msg << "Could not find cell " << col_id
                << " specified in plot " << id_;
        fatal_error(err_msg);
      }
    } else if (PlotColorBy::mats == color_by_) {
      if (material_map.find(col_id) != material_map.end()) {
        col_id = material_map[col_id];
        colors_[col_id] = user_rgb;
      } else {
        std::stringstream err_msg;
        err_msg << "Could not find material " << col_id
                << " specified in plot " << id_;
        fatal_error(err_msg);
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
      std::stringstream msg;
      msg << "Meshlines ignored in voxel plot " << id_;
      warning(msg);
    }

    if (mesh_line_nodes.size() == 1) {
      // Get first meshline node
      pugi::xml_node meshlines_node = mesh_line_nodes[0].node();

      // Check mesh type
      std::string meshtype;
      if (check_for_node(meshlines_node, "meshtype")) {
        meshtype = get_node_value(meshlines_node, "meshtype");
      } else {
        std::stringstream err_msg;
        err_msg << "Must specify a meshtype for meshlines specification in plot " << id_;
        fatal_error(err_msg);
      }

      // Ensure that there is a linewidth for this meshlines specification
      std::string meshline_width;
      if (check_for_node(meshlines_node, "linewidth")) {
        meshline_width = get_node_value(meshlines_node, "linewidth");
        meshlines_width_ = std::stoi(meshline_width);
      } else {
        std::stringstream err_msg;
        err_msg << "Must specify a linewidth for meshlines specification in plot " << id_;
        fatal_error(err_msg);
      }

      // Check for color
      if (check_for_node(meshlines_node, "color")) {
        // Check and make sure 3 values are specified for RGB
        std::vector<int> ml_rgb = get_node_array<int>(meshlines_node, "color");
        if (ml_rgb.size() != 3) {
          std::stringstream err_msg;
          err_msg << "Bad RGB for meshlines color in plot " << id_;
          fatal_error(err_msg);
        }
        meshlines_color_ = ml_rgb;
      }

      // Set mesh based on type
      if ("ufs" == meshtype) {
        if (settings::index_ufs_mesh < 0) {
          std::stringstream err_msg;
          err_msg << "No UFS mesh for meshlines on plot " << id_;
          fatal_error(err_msg);
        } else {
          index_meshlines_mesh_ = settings::index_ufs_mesh;
        }
      } else if ("cmfd" == meshtype) {
        if (!settings::cmfd_run) {
          std::stringstream err_msg;
          err_msg << "Need CMFD run to plot CMFD mesh for meshlines on plot " << id_;
          fatal_error(err_msg);
        } else {
          index_meshlines_mesh_ = settings::index_cmfd_mesh;
        }
      } else if ("entropy" == meshtype) {
        if (settings::index_entropy_mesh < 0) {
          std::stringstream err_msg;
          err_msg <<"No entropy mesh for meshlines on plot " << id_;
          fatal_error(err_msg);
        } else {
          index_meshlines_mesh_ = settings::index_entropy_mesh;
        }
      } else if ("tally" == meshtype) {
        // Ensure that there is a mesh id if the type is tally
        int tally_mesh_id;
        if (check_for_node(meshlines_node, "id")) {
          tally_mesh_id = std::stoi(get_node_value(meshlines_node, "id"));
        } else {
          std::stringstream err_msg;
          err_msg << "Must specify a mesh id for meshlines tally "
                  << "mesh specification in plot " << id_;
          fatal_error(err_msg);
        }
        // find the tally index
        int idx;
        int err = openmc_get_mesh_index(tally_mesh_id, &idx);
        if (err != 0) {
          std::stringstream err_msg;
          err_msg << "Could not find mesh " << tally_mesh_id
                  << " specified in meshlines for plot " << id_;
          fatal_error(err_msg);
        }
        index_meshlines_mesh_ = idx;
      } else {
        std::stringstream err_msg;
        err_msg << "Invalid type for meshlines on plot " << id_ ;
        fatal_error(err_msg);
      }
    } else {
      std::stringstream err_msg;
      err_msg << "Mutliple meshlines specified in plot " << id_;
      fatal_error(err_msg);
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
      if (openmc_master) {
        std::stringstream wrn_msg;
        wrn_msg << "Mask ignored in voxel plot " << id_;
        warning(wrn_msg);
      }
    }

    if (mask_nodes.size() == 1) {
      // Get pointer to mask
      pugi::xml_node mask_node = mask_nodes[0].node();

      // Determine how many components there are and allocate
      std::vector<int> iarray = get_node_array<int>(mask_node, "components");      
      if (iarray.size() == 0) {
        std::stringstream err_msg;
        err_msg << "Missing <components> in mask of plot " << id_;
        fatal_error(err_msg);
      }

      // First we need to change the user-specified identifiers to indices
      // in the cell and material arrays
      for (auto& col_id : iarray) {
        if (PlotColorBy::cells == color_by_) {
          if (cell_map.find(col_id) != cell_map.end()) {
            col_id  = cell_map[col_id];
          }
          else {
            std::stringstream err_msg;
            err_msg << "Could not find cell " << col_id
                    << " specified in the mask in plot " << id_;
            fatal_error(err_msg);
          }
        } else if (PlotColorBy::mats == color_by_) {
          if (material_map.find(col_id) != material_map.end()) {
            col_id = material_map[col_id];
          }
          else {
            std::stringstream err_msg;
            err_msg << "Could not find material " << col_id
                    << " specified in the mask in plot " << id_;
            fatal_error(err_msg);
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
      std::stringstream err_msg;
      err_msg << "Mutliple masks specified in plot " << id_;
      fatal_error(err_msg);
    }
  }
}

Plot::Plot(pugi::xml_node plot_node):
index_meshlines_mesh_(-1)
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
} // End Plot constructor

//==============================================================================
// POSITION_RGB computes the red/green/blue values for a given plot with the
// current particle's position
//==============================================================================


void position_rgb(Particle p, Plot pl, RGBColor& rgb, int& id)
{
  p.n_coord = 1;

  bool found_cell = find_cell(&p, 0);

  int j = p.n_coord - 1;

  if (settings::check_overlaps) {check_cell_overlap(&p);}

  // Set coordinate level if specified
  if (pl.level_ >= 0) {j = pl.level_ + 1;}

  if (!found_cell) {
    // If no cell, revert to default color
    rgb = pl.not_found_;
    id = -1;
  } else {
    if (PlotColorBy::mats == pl.color_by_) {
      // Assign color based on material
      Cell* c = cells[p.coord[j].cell];
      if (c->type_ == FILL_UNIVERSE) {
        // If we stopped on a middle universe level, treat as if not found
        rgb = pl.not_found_;
        id = -1;
      } else if (p.material == MATERIAL_VOID) {
        // By default, color void cells white
        rgb = WHITE;
        id = -1;
      } else {
        rgb = pl.colors_[p.material - 1];
        id = materials[p.material - 1]->id_;
      }
    } else if (PlotColorBy::cells == pl.color_by_) {
      // Assign color based on cell
      rgb = pl.colors_[p.coord[j].cell];
      id = cells[p.coord[j].cell]->id_;
    }
  } // endif found_cell
}

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
  of << "P6" << "\n";
  of << pl.pixels_[0] << " " << pl.pixels_[1] << "\n";
  of << "255" << "\n";
  of.close();

  of.open(fname, std::ios::binary | std::ios::app);
  // Write color for each pixel
  for (int y = 0; y < pl.pixels_[1]; y++) {
    for (int x = 0; x < pl.pixels_[0]; x++) {
      RGBColor rgb = data(x,y);
      of << rgb.red << rgb.green << rgb.blue;
    }
  }

  // Close file
  // THIS IS HERE TO MATCH FORTRAN VERSION, NOT TECHNICALLY NECESSARY
  of << "\n";
  of.close();
}

//==============================================================================
// DRAW_MESH_LINES draws mesh line boundaries on an image
//==============================================================================

void draw_mesh_lines(Plot pl, ImageData& data)
{
  RGBColor rgb;
  rgb = pl.meshlines_color_;

  int outer, inner;
  switch(pl.basis_) {
  case PlotBasis::xy :
    outer = 0;
    inner = 1;
    break;
  case PlotBasis::xz :
    outer = 0;
    inner = 2;
    break;
  case PlotBasis::yz :
    outer = 1;
    inner = 2;
    break;
  }

  double xyz_ll_plot[3], xyz_ur_plot[3];
  xyz_ll_plot[0] = pl.origin_[0];
  xyz_ll_plot[1] = pl.origin_[1];
  xyz_ll_plot[2] = pl.origin_[2];

  xyz_ur_plot[0] = pl.origin_[0];
  xyz_ur_plot[1] = pl.origin_[1];
  xyz_ur_plot[2] = pl.origin_[2];
  
  xyz_ll_plot[outer] = pl.origin_[outer] - pl.width_[0] / 2.;
  xyz_ll_plot[inner] = pl.origin_[inner] - pl.width_[1] / 2.;
  xyz_ur_plot[outer] = pl.origin_[outer] + pl.width_[0] / 2.;
  xyz_ur_plot[inner] = pl.origin_[inner] + pl.width_[1] / 2.;

  int width[3];
  width[0] = xyz_ur_plot[0] - xyz_ll_plot[0];
  width[1] = xyz_ur_plot[1] - xyz_ll_plot[1];
  width[2] = xyz_ur_plot[2] - xyz_ll_plot[2];

  auto& m = meshes[pl.index_meshlines_mesh_];

  int ijk_ll[3], ijk_ur[3];
  bool in_mesh;
  m->get_indices(Position(xyz_ll_plot), &(ijk_ll[0]), &in_mesh);
  m->get_indices(Position(xyz_ur_plot), &(ijk_ur[0]), &in_mesh);

  // Fortran/C++ index correction
  ijk_ur[0]++; ijk_ur[1]++; ijk_ur[2]++;

  double xyz_ll[3], xyz_ur[3];
  // sweep through all meshbins on this plane and draw borders
  for (int i = ijk_ll[outer]; i <= ijk_ur[outer]; i++) {
    for (int j = ijk_ll[inner]; j <= ijk_ur[inner]; j++) {
      // check if we're in the mesh for this ijk
      if (i > 0 && i <= m->shape_[outer] && j >0 && j <= m->shape_[inner] ) {
        int outrange[3], inrange[3];
        // get xyz's of lower left and upper right of this mesh cell
        xyz_ll[outer] = m->lower_left_[outer] + m->width_[outer] * (i - 1);
        xyz_ll[inner] = m->lower_left_[inner] + m->width_[inner] * (j - 1);
        xyz_ur[outer] = m->lower_left_[outer] + m->width_[outer] * i;
        xyz_ur[inner] = m->lower_left_[inner] + m->width_[inner] * j;

        // map the xyz ranges to pixel ranges
        double frac = (xyz_ll[outer] - xyz_ll_plot[outer]) / width[outer];
        outrange[0] = int(frac * double(pl.pixels_[0]));
        frac = (xyz_ur[outer] - xyz_ll_plot[outer]) / width[outer];
        outrange[1] = int(frac * double(pl.pixels_[0]));

        frac = (xyz_ur[inner] - xyz_ll_plot[inner]) / width[inner];
        inrange[0] = int((1. - frac) * (double)pl.pixels_[1]);
        frac = (xyz_ll[inner] - xyz_ll_plot[inner]) / width[inner];
        inrange[1] = int((1. - frac) * (double)pl.pixels_[1]);

        // draw lines
        for (int out_ = outrange[0]; out_ <= outrange[1]; out_++) {
          for (int plus = 0; plus <= pl.meshlines_width_; plus++) {
            data(out_, inrange[0] + plus) = rgb;
            data(out_, inrange[1] + plus) = rgb;
            data(out_, inrange[0] - plus) = rgb;
            data(out_, inrange[1] - plus) = rgb;
          }
        }

        for (int in_ = inrange[0]; in_ <= inrange[1]; in_++) {
          for (int plus = 0; plus <= pl.meshlines_width_; plus++) {
            data(outrange[0] + plus, in_) = rgb;
            data(outrange[1] + plus, in_) = rgb;
            data(outrange[0] - plus, in_) = rgb;
            data(outrange[1] - plus, in_) = rgb;
          }
        }

      } // end if(in mesh)
    }
  } // end outer loops
} // end draw_mesh_lines

//==============================================================================
// CREATE_VOXEL outputs a binary file that can be input into silomesh for 3D
// geometry visualization.  It works the same way as create_ppm by dragging a
// particle across the geometry for the specified number of voxels. The first 3
// int(4)'s in the binary are the number of x, y, and z voxels.  The next 3
// real(8)'s are the widths of the voxels in the x, y, and z directions. The
// next 3 real(8)'s are the x, y, and z coordinates of the lower left
// point. Finally the binary is filled with entries of four int(4)'s each. Each
// 'row' in the binary contains four int(4)'s: 3 for x,y,z position and 1 for
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
  std::array<double, 3> ll;
  ll[0] = pl.origin_[0] - pl.width_[0] / 2.;
  ll[1] = pl.origin_[1] - pl.width_[1] / 2.;
  ll[2] = pl.origin_[2] - pl.width_[2] / 2.;

  // allocate and initialize particle
  double dir[3] = {0.5, 0.5, 0.5};
  Particle p;
  p.initialize();
  std::copy(ll.begin(), ll.begin()+ll.size(), p.coord[0].xyz);
  std::copy(dir, dir+3, p.coord[0].uvw);
  p.coord[0].universe = openmc_root_universe;

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
  hsize_t three = 3;
  write_attribute(file_id, "num_voxels", pl.pixels_);
  write_attribute(file_id, "voxel_width", vox);
  write_attribute(file_id, "lower_left", ll);

  // Create dataset for voxel data -- note that the dimensions are reversed
  // since we want the order in the file to be z, y, x
  hsize_t dims[3];
  dims[0] = pl.pixels_[0];
  dims[1] = pl.pixels_[1];
  dims[2] = pl.pixels_[2];
  hid_t dspace, dset, memspace;
  voxel_init(file_id, &(dims[0]), &dspace, &dset, &memspace);

  // move to center of voxels
  ll[0] = ll[0] + vox[0] / 2.;
  ll[1] = ll[1] + vox[1] / 2.;
  ll[2] = ll[2] + vox[2] / 2.;

  int data[pl.pixels_[1]][pl.pixels_[2]];

  ProgressBar pb;

  RGBColor rgb;
  int id;
  for (int x = 0; x < pl.pixels_[0]; x++) {
    pb.set_value(100.*(double)x/(double)(pl.pixels_[0]-1));
    for (int y = 0; y < pl.pixels_[1]; y++) {
      for (int z = 0; z < pl.pixels_[2]; z++) {
        // get voxel color
        position_rgb(p, pl, rgb, id);
        // write to plot data
        data[y][z] = id;
        // advance particle in z direction
        p.coord[0].xyz[2] = p.coord[0].xyz[2] + vox[2];
      }
      // advance particle in y direction
      p.coord[0].xyz[1] = p.coord[0].xyz[1] + vox[1];
      p.coord[0].xyz[2] = ll[2];
    }
    // advance particle in x direction
    p.coord[0].xyz[0] = p.coord[0].xyz[0] + vox[0];
    p.coord[0].xyz[1] = ll[1];
    p.coord[0].xyz[2] = ll[2];
    // Write to HDF5 dataset
    voxel_write_slice(x, dspace, dset, memspace, &(data[0]));
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

RGBColor random_color() {
  RGBColor rgb;
  rgb.red = int(prn()*255);
  rgb.green = int(prn()*255);
  rgb.blue = int(prn()*255);
  return rgb;
  //  return {int(prn()*255), int(prn()*255), int(prn()*255)};
}

} // namespace openmc
