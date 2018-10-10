#include <fstream>

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

namespace openmc {

const int RED   = 0;
const int GREEN = 1;
const int BLUE  = 2;

const int WHITE[3] = {255, 255, 255};
const int NULLRGB[3] = {0, 0, 0};

//===============================================================================
// RUN_PLOT controls the logic for making one or many plots
//===============================================================================

int openmc_plot_geometry() {
  int err;

  for(int i = 0; i < n_plots; i++) {
    ObjectPlot* pl = plots[i];

    std::stringstream ss;
    ss << "Processing plot " << pl->id << ": "
       << pl->path_plot << "...";
      write_message(ss.str(), 5);

      if (PLOT_TYPE::SLICE == pl->type) {
        // create 2D image
        create_ppm(pl);
        continue;
      } else if (PLOT_TYPE::VOXEL == pl->type) {
        // create voxel file for 3D viewing
        create_voxel(pl);
        continue;
      }

  }
  return 0;
}

void read_plots(pugi::xml_node plots_node) {

  std::vector<pugi::xml_node> plot_nodes;
  plot_nodes = get_child_nodes(plots_node, "plot");

  n_plots = plot_nodes.size();

  for(auto plot : plot_nodes) {
    //    ObjectPlot* pl = new ObjectPlot(plot);
  }
}

//===============================================================================
// CREATE_PPM creates an image based on user input from a plots.xml <plot>
// specification in the portable pixmap format (PPM)
//===============================================================================

void create_ppm(ObjectPlot* pl) {

  int width = pl->pixels[0];
  int height = pl->pixels[1];

  double in_pixel = (pl->width[0])/double(width);
  double out_pixel = (pl->width[1])/double(height);

  ImageData data;
  data.resize(width);
  for (auto & i : data) {
    i.resize(height);
    for (auto & j : i) {
      j.resize(3);
    }
  }

  int in_i, out_i;
  double xyz[3];
  switch(pl->basis) {
  case PLOT_BASIS::XY :
    in_i = 0;
    out_i = 1;
    xyz[0] = pl->origin[0] - pl->width[0] / TWO;
    xyz[1] = pl->origin[1] + pl->width[1] / TWO;
    xyz[2] = pl->origin[2];
    break;
  case PLOT_BASIS::XZ :
    in_i = 0;
    out_i = 2;
    xyz[0] = pl->origin[0] - pl->width[0] / TWO;
    xyz[1] = pl->origin[1];
    xyz[2] = pl->origin[2] + pl->width[1] / TWO;
    break;
  case PLOT_BASIS::YZ :
    in_i = 1;
    out_i = 2;
    xyz[0] = pl->origin[0];
    xyz[1] = pl->origin[1] - pl->width[0] / TWO;
    xyz[2] = pl->origin[2] + pl->width[1] / TWO;
    break;
  }

  double dir[3] = {HALF, HALF, HALF};
  Particle *p = new Particle();
  p->initialize();
  std::copy(xyz, xyz+3, p->coord[0].xyz);
  std::copy(dir, dir+3, p->coord[0].uvw);
  p->coord[0].universe = openmc_root_universe;

  // local variables
  int rgb[3];
  int id;
  for (int y = 0; y < height; y++) {
    p->coord[0].xyz[out_i] = xyz[out_i] - out_pixel*(y);
    for (int x = 0; x < width; x++) {
      p->coord[0].xyz[in_i] = xyz[in_i] + in_pixel*(x);
      position_rgb(p, pl, rgb, id);
      data[x][y][RED] = rgb[RED];
      data[x][y][GREEN] = rgb[GREEN];
      data[x][y][BLUE] = rgb[BLUE];
    }
  }

  if (pl->index_meshlines_mesh >= 0) { draw_mesh_lines(pl, data); }

  output_ppm(pl, data);

}

ObjectPlot::ObjectPlot(pugi::xml_node plot_node) {

  // Copy data into plots  
  if (check_for_node(plot_node, "id")) {
    id = std::stoi(get_node_value(plot_node, "id"));
  } else {
    fatal_error("Must specify plot id in plots XML file.");
  }

  // Check to make sure 'id' hasn't been used
  if (plot_dict.find(id) != plot_dict.end()) {
    std::stringstream err_msg;
    err_msg << "Two or more plots use the same unique ID: ";
    err_msg << id;
    fatal_error(err_msg.str());
  }

  // Copy plot type
  // Default is slice
  std::string type_str = "slice"; 
  if (check_for_node(plot_node, "type")) {
    type_str = get_node_value(plot_node, "type");
    to_lower(type_str);
    if (type_str == "slice") {
      type = PLOT_TYPE::SLICE;
    }
    else if (type_str == "voxel") {
      type = PLOT_TYPE::VOXEL;
    } else {
      std::stringstream err_msg;
      err_msg << "Unsupported plot type '" << type_str
              << "' in plot " << id;
      fatal_error(err_msg.str());
    }
  }

  // Set output file path
  std::stringstream filename;
  filename << "plot_" << id;

  if (check_for_node(plot_node, "filename")) {
    switch(type) {
    case PLOT_TYPE::SLICE:
      filename << ".ppm";
      break;
    case PLOT_TYPE::VOXEL:
      filename << ".h5";
      break;
    }
  }

  // Copy plot pixel size
  std::vector<int> pxls;
  if (PLOT_TYPE::SLICE == type) {
    if (node_word_count(plot_node, "pixels") == 2) {
      pxls = get_node_array<int>(plot_node, "pixels");
      pixels[0] = pxls[0];
      pixels[1] = pxls[1];
    } else {
      std::stringstream err_msg;
      err_msg << "<pixels> must be length 2 in slice plot "
              << id;
      fatal_error(err_msg.str());
    }
  } else if (PLOT_TYPE::VOXEL == type) {
    if (node_word_count(plot_node, "pixels") == 3) {
      pxls = get_node_array<int>(plot_node, "pixels");
      pixels[0] = pxls[0];
      pixels[1] = pxls[1];
      pixels[2] = pxls[2];      
    } else {
      std::stringstream err_msg;
      err_msg << "<pixels> must be length 3 in voxel plot "
              << id;
      fatal_error(err_msg.str());
    }
  }

  // Copy plot background color
  std::vector<int> bg_rgb;
  if (check_for_node(plot_node, "background")) {
    if (PLOT_TYPE::VOXEL == type) {
      if (openmc_master) {
        std::stringstream err_msg;
        err_msg << "Background color ignored in voxel plot "
                << id;
        warning(err_msg.str());
      }
    }
    if (node_word_count(plot_node, "background") == 3) {
      bg_rgb = get_node_array<int>(plot_node, "background");
      not_found.rgb[0] = bg_rgb[0];
      not_found.rgb[1] = bg_rgb[1];
      not_found.rgb[2] = bg_rgb[2];
    } else {
      std::stringstream err_msg;
      err_msg << "Bad background RGB in plot "
              << id;
      fatal_error(err_msg);
    }
  } else {
    // default to a white background
    not_found.rgb[0] = 255;
    not_found.rgb[1] = 255;
    not_found.rgb[2] = 255;
  }

  // Copy plot basis
  if (PLOT_TYPE::SLICE == type) {
    std::string pl_basis = "xy";
    if (check_for_node(plot_node, "basis")) {
      pl_basis = get_node_value(plot_node, "basis");
    }
    to_lower(pl_basis);
    if ("xy" == pl_basis) {
      basis = PLOT_BASIS::XY;
    } else if ("xz" == pl_basis) {
      basis = PLOT_BASIS::XZ;
    } else if ("yz" == pl_basis) {
      basis = PLOT_BASIS::YZ;
    } else {
      std::stringstream err_msg;
      err_msg << "Unsupported plot basis '" << pl_basis
              << "' in plot " << id;
      fatal_error(err_msg);
    }
  }

  // Copy plotting origin
  std::vector<double> pl_origin;
  if (node_word_count(plot_node, "origin") == 3) {
    pl_origin = get_node_array<double>(plot_node, "origin");
    origin[0] = pl_origin[0];
    origin[1] = pl_origin[1];
    origin[2] = pl_origin[2];
  } else {
    std::stringstream err_msg;
    err_msg << "Origin must be length 3 in plot "
            << id;
    fatal_error(err_msg);
  }

  // Copy plotting width
  std::vector<int> pl_width;
  if (PLOT_TYPE::SLICE == type) {
    if (node_word_count(plot_node, "width") == 2) {
      pl_width = get_node_array<int>(plot_node, "width");
      width[0] = pl_width[0];
      width[1] = pl_width[1];
    } else {
      std::stringstream err_msg;
      err_msg << "<width> must be length 2 in slice plot "
              << id;
      fatal_error(err_msg);
    }
  } else if (PLOT_TYPE::VOXEL == type) {
    if (node_word_count(plot_node, "width") == 3) {
      pl_width = get_node_array<int>(plot_node, "width");
      width[0] = pl_width[0];
      width[1] = pl_width[1];
      width[2] = pl_width[2];
    } else {
      std::stringstream err_msg;
      err_msg << "<width> must be length 3 in voxel plot "
              << id;
      fatal_error(err_msg);
    }
  }

  // Copy plot universe level
  if (check_for_node(plot_node, "level")) {
    level = std::stoi(get_node_value(plot_node, "level"));
    if (level < 0) {
      std::stringstream err_msg;
      err_msg << "Bad universe level in plot " << id ;
      fatal_error(err_msg);
    }
  } else {
    level = PLOT_LEVEL_LOWEST;
  }

  // Copy plot color type and initialize all colors randomly
  std::string pl_color_by = "cell";
  if (check_for_node(plot_node, "color_by")) {
    pl_color_by = get_node_value(plot_node, "color_by");
    to_lower(pl_color_by);
  }
  if ("cell" == pl_color_by) {
    color_by = PLOT_COLOR_BY::CELLS;

    //    colors.resize(n_cells);
    for(int i = 0; i < n_cells; i++) {
      colors[i].rgb[0] = int(prn()*255);
      colors[i].rgb[1] = int(prn()*255);
      colors[i].rgb[2] = int(prn()*255);
    }

  } else if("material" == pl_color_by) {

    color_by = PLOT_COLOR_BY::MATS;
      
    //    colors.resize(materials.size());
    for(int i = 0; i < n_cells; i++) {
      colors[i].rgb[0] = int(prn()*255);
      colors[i].rgb[1] = int(prn()*255);
      colors[i].rgb[2] = int(prn()*255);
    }
  } else {
    std::stringstream err_msg;
    err_msg << "Unsupported plot color type '" << pl_color_by
            << "' in plot " << id;
    fatal_error(err_msg);
  }

  // Get the number of <color> nodes and get a list of them
  std::vector<pugi::xml_node> color_nodes;
  color_nodes = get_child_nodes(plot_node, "color");
    
  // Copy user-specified colors
  if (0 != color_nodes.size()) {

    if (PLOT_TYPE::VOXEL == type) {
      if (openmc_master) {
        std::stringstream err_msg;
        err_msg << "Color specifications ignored in voxel plot "
                << id;
        warning(err_msg);
      }
    }

    for(auto cn : color_nodes) {

      // Check and make sure 3 values are specified for RGB
      if (node_word_count(cn, "rgb") != 3) {
        std::stringstream err_msg;
        err_msg << "Bad RGB in plot " << id;
        fatal_error(err_msg);
      }

      // Ensure that there is an id for this color specification
      // Ensure that there is an id for this color specification
      int col_id;
      if (check_for_node(cn, "id")) {
        col_id = std::stoi(get_node_value(cn, "id"));
      } else {
        std::stringstream err_msg;
        err_msg << "Must specify id for color specification in plot "
                << id;
        fatal_error(err_msg);
      }
      // Add RGB
      if (PLOT_COLOR_BY::CELLS == color_by) {
        std::vector<int> cell_rgb;
        if (cell_map.find(col_id) != cell_map.end()) {
          col_id = cell_map[col_id];
          cell_rgb = get_node_array<int>(cn, "rgb");
          colors[col_id].rgb[0] = cell_rgb[0];
          colors[col_id].rgb[1] = cell_rgb[1];
          colors[col_id].rgb[2] = cell_rgb[2];          
        } else {
          std::stringstream err_msg;
          err_msg << "Could not find cell " << col_id
                  << " specified in plot " << id;
          fatal_error(err_msg);
        }
      } else if (PLOT_COLOR_BY::MATS == color_by) {
        std::vector<int> mat_rgb;
        if (material_map.find(col_id) != material_map.end()) {
          col_id = material_map[col_id];
          mat_rgb = get_node_array<int>(cn, "rgb");
          colors[col_id].rgb[0] = mat_rgb[0];
          colors[col_id].rgb[1] = mat_rgb[1];
          colors[col_id].rgb[2] = mat_rgb[2];          
        } else {
          std::stringstream err_msg;
          err_msg << "Could not find material " << col_id
                  << " specified in plot " << id;            
          fatal_error(err_msg);
        }
      }
    } // color node loop
  }
      

  // Deal with meshlines
  std::vector<pugi::xml_node> mesh_line_nodes;
  mesh_line_nodes = get_child_nodes(plot_node, "meshlines");

  int n_meshlines = mesh_line_nodes.size();

  if (n_meshlines != 0) {

    if (PLOT_TYPE::VOXEL == type) {
      std::stringstream msg;
      msg << "Meshlines ignored in voxel plot " << id;
      warning(msg);
    }

    if (0 == n_meshlines) {
      // Skip if no meshlines are specified
    } else if (1 == n_meshlines) {
      // Get first meshline node
      pugi::xml_node meshlines_node = mesh_line_nodes[0];
      
      // Check mesh type
      std::string meshline_type;
      if (check_for_node(meshlines_node, "meshtype")) {
        meshline_type = get_node_value(meshlines_node, "meshtype");
      } else {
        std::stringstream err_msg;
        err_msg << "Must specify a meshtype for meshlines specification in plot " << id;
        fatal_error(err_msg);
      }

      // Check mesh type
      std::string meshline_width;
      if (check_for_node(meshlines_node, "linewidth")) {
        meshline_width = get_node_value(meshlines_node, "linewidth");
      } else {
        std::stringstream err_msg;
        err_msg << "Must specify a linewidth for meshlines specification in plot " << id;
        fatal_error(err_msg);
      }

      // Check for color
      std::vector<int> ml_rgb;
      if (check_for_node(meshlines_node, "color")) {
        // Check and make sure 3 values are specified for RGB
        if (node_word_count(meshlines_node, "color") != 3) {
          std::stringstream err_msg;
          err_msg << "Bad RGB for meshlines color in plot " << id;
          fatal_error(err_msg);
        }
        ml_rgb = get_node_array<int>(meshlines_node, "color");
        meshlines_color.rgb[0] = ml_rgb[0];
        meshlines_color.rgb[1] = ml_rgb[1];
        meshlines_color.rgb[2] = ml_rgb[2];
      } else {
        meshlines_color.rgb[0] = 0;
        meshlines_color.rgb[1] = 0;
        meshlines_color.rgb[2] = 0;
      }
      
    } 
  }
} // End ObjectPlot constructor
  
//===============================================================================
// POSITION_RGB computes the red/green/blue values for a given plot with the
// current particle's position
//===============================================================================

void position_rgb(Particle* p, ObjectPlot* pl, int rgb[3], int &id) {
  bool found_cell;

  p->n_coord = 1;

  found_cell = find_cell(p, 0);

  int j = p->n_coord - 1;

  if (settings::check_overlaps) { check_cell_overlap(p); }

  // Set coordinate level if specified
  if (pl->level >= 0) {j = pl->level + 1;}

  Cell* c;

  if (!found_cell) {
    // If no cell, revert to default color
    std::copy(pl->not_found.rgb,
              pl->not_found.rgb + 3,
              rgb);
    id = -1;
  } else {
    if (PLOT_COLOR_BY::MATS == pl->color_by) {
      // Assign color based on material
      c = cells[p->coord[j].cell];
      if (c->type_ == FILL_UNIVERSE) {
        // If we stopped on a middle universe level, treat as if not found
        std::copy(pl->not_found.rgb,
                  pl->not_found.rgb + 3,
                  rgb);
        id = -1;
      } else if (p->material == MATERIAL_VOID) {
        // By default, color void cells white
        std::copy(WHITE, WHITE+3, rgb);
        id = -1;
      } else {
        std::copy(pl->colors[p->material - 1].rgb,
                  pl->colors[p->material - 1].rgb + 3,
                  rgb);
        id = materials[p->material - 1]->id;
      }

    } else if (PLOT_COLOR_BY::CELLS == pl->color_by) {
      // Assign color based on cell
      std::copy(pl->colors[p->coord[j].cell].rgb,
                pl->colors[p->coord[j].cell].rgb + 3,
                rgb);
      id = cells[p->coord[j].cell]->id_;
    } else {
      std::copy(NULLRGB, NULLRGB+3, rgb);
      id = -1;
    }

  } // endif found_cell
}

//===============================================================================
// OUTPUT_PPM writes out a previously generated image to a PPM file
//===============================================================================

void output_ppm(ObjectPlot* pl, const ImageData &data)
{
  // Open PPM file for writing
  std::string fname = std::string(pl->path_plot);
  fname = strtrim(fname);
  std::ofstream of;

  of.open(fname);

  // Write header
  of << "P6" << std::endl;
  of << pl->pixels[0] << " " << pl->pixels[1] << std::endl;
  of << "255" << std::endl;
  of.close();

  of.open(fname, std::ios::binary | std::ios::app);
  // Write color for each pixel
  for (int y = 0; y < pl->pixels[1]; y++) {
    for (int x = 0; x < pl->pixels[0]; x++) {
      std::vector<int> rgb = data[x][y];
      of.write((char*)&rgb[RED], 1);
      of.write((char*)&rgb[GREEN], 1);
      of.write((char*)&rgb[BLUE], 1);
    }
  }
  // Close file
  // THIS IS HERE TO MATCH FORTRAN VERSION, NOT NECESSARY
  of << std::endl;
  of.close();
}

//===============================================================================
// DRAW_MESH_LINES draws mesh line boundaries on an image
//===============================================================================

void draw_mesh_lines(ObjectPlot *pl, ImageData &data)
{
  std::vector<int> rgb; rgb.resize(3);
  rgb[RED] = pl->meshlines_color.rgb[RED];
  rgb[GREEN] = pl->meshlines_color.rgb[GREEN];
  rgb[BLUE] = pl->meshlines_color.rgb[BLUE];

  int outer, inner;
  switch(pl->basis){
  case PLOT_BASIS::XY :
    outer = 0;
    inner = 1;
    break;
  case PLOT_BASIS::XZ :
    outer = 0;
    inner = 2;
    break;
  case PLOT_BASIS::YZ :
    outer = 1;
    inner = 2;
    break;
  }

  double xyz_ll_plot[3], xyz_ur_plot[3];
  std::copy((double*)&pl->origin, (double*)&pl->origin + 3, xyz_ll_plot);
  std::copy((double*)&pl->origin, (double*)&pl->origin + 3, xyz_ur_plot);

  xyz_ll_plot[outer] = pl->origin[outer] - pl->width[0] / TWO;
  xyz_ll_plot[inner] = pl->origin[inner] - pl->width[1] / TWO;
  xyz_ur_plot[outer] = pl->origin[outer] + pl->width[0] / TWO;
  xyz_ur_plot[inner] = pl->origin[inner] + pl->width[1] / TWO;

  int width[3];
  width[0] = xyz_ur_plot[0] - xyz_ll_plot[0];
  width[1] = xyz_ur_plot[1] - xyz_ll_plot[1];
  width[2] = xyz_ur_plot[2] - xyz_ll_plot[2];

  auto &m = meshes[pl->index_meshlines_mesh];

  int ijk_ll[3], ijk_ur[3];
  bool in_mesh;
  m->get_indices(Position(xyz_ll_plot), &(ijk_ll[0]), &in_mesh);
  m->get_indices(Position(xyz_ur_plot), &(ijk_ur[0]), &in_mesh);

  // Fortran/C++ index correction
  ijk_ur[0]++; ijk_ur[1]++; ijk_ur[2]++;

  double frac;
  int outrange[3], inrange[3];
  double xyz_ll[3], xyz_ur[3];
  // sweep through all meshbins on this plane and draw borders
  for (int i = ijk_ll[outer]; i <= ijk_ur[outer]; i++) {
    for (int j = ijk_ll[inner]; j <= ijk_ur[inner]; j++) {
      // check if we're in the mesh for this ijk
      if (i > 0 && i <= m->shape_[outer] && j >0 && j <= m->shape_[inner] ) {

        // get xyz's of lower left and upper right of this mesh cell
        xyz_ll[outer] = m->lower_left_[outer] + m->width_[outer] * (i - 1);
        xyz_ll[inner] = m->lower_left_[inner] + m->width_[inner] * (j - 1);
        xyz_ur[outer] = m->lower_left_[outer] + m->width_[outer] * i;
        xyz_ur[inner] = m->lower_left_[inner] + m->width_[inner] * j;

        // map the xyz ranges to pixel ranges
        frac = (xyz_ll[outer] - xyz_ll_plot[outer]) / width[outer];
        outrange[0] = int(frac * double(pl->pixels[0]));
        frac = (xyz_ur[outer] - xyz_ll_plot[outer]) / width[outer];
        outrange[1] = int(frac * double(pl->pixels[0]));

        frac = (xyz_ur[inner] - xyz_ll_plot[inner]) / width[inner];
        inrange[0] = int((ONE - frac) * (double)pl->pixels[1]);
        frac = (xyz_ll[inner] - xyz_ll_plot[inner]) / width[inner];
        inrange[1] = int((ONE - frac) * (double)pl->pixels[1]);

        // draw lines
        for (int out_ = outrange[0]; out_ <= outrange[1]; out_++) {
          for (int plus = 0; plus <= pl->meshlines_width; plus++) {
            data[out_][inrange[0] + plus] = rgb;
            data[out_][inrange[1] + plus] = rgb;
            data[out_][inrange[0] - plus] = rgb;
            data[out_][inrange[1] - plus] = rgb;
          }
        }

        for (int in_ = inrange[0]; in_ <= inrange[1]; in_++) {
          for (int plus = 0; plus <= pl->meshlines_width; plus++) {
            data[outrange[0] + plus][in_] = rgb;
            data[outrange[1] + plus][in_] = rgb;
            data[outrange[0] - plus][in_] = rgb;
            data[outrange[1] - plus][in_] = rgb;
          }
        }

      } // end if(in mesh)
    }
  } // end outer loops

}

//===============================================================================
// CREATE_VOXEL outputs a binary file that can be input into silomesh for 3D
// geometry visualization.  It works the same way as create_ppm by dragging a
// particle across the geometry for the specified number of voxels. The first
// 3 int(4)'s in the binary are the number of x, y, and z voxels.  The next 3
// real(8)'s are the widths of the voxels in the x, y, and z directions. The next
// 3 real(8)'s are the x, y, and z coordinates of the lower left point. Finally
// the binary is filled with entries of four int(4)'s each. Each 'row' in the
// binary contains four int(4)'s: 3 for x,y,z position and 1 for cell or material
// id.  For 1 million voxels this produces a file of approximately 15MB.
//===============================================================================

void create_voxel(ObjectPlot *pl) {

  // compute voxel widths in each direction
  double vox[3];
  vox[0] = pl->width[0]/(double)pl->pixels[0];
  vox[1] = pl->width[1]/(double)pl->pixels[1];
  vox[2] = pl->width[2]/(double)pl->pixels[2];

  // initial particle position
  double ll[3];
  ll[0] = pl->origin[0] - pl->width[0] / TWO;
  ll[1] = pl->origin[1] - pl->width[1] / TWO;
  ll[2] = pl->origin[2] - pl->width[2] / TWO;

  // allocate and initialize particle
  double dir[3] = {HALF, HALF, HALF};
  Particle *p = new Particle();
  p->initialize();
  std::copy(ll, ll + 3, p->coord[0].xyz);
  std::copy(dir, dir + 3, p->coord[0].uvw);
  p->coord[0].universe = openmc_root_universe;

  // Open binary plot file for writing
  std::ofstream of;
  std::string fname = std::string(pl->path_plot);
  fname = strtrim(fname);
  hid_t file_id = file_open(fname, 'w');

  // write header info
  write_attribute(file_id, "filetype", 'voxel');
  write_attribute(file_id, "version", VERSION_VOXEL);
  write_attribute(file_id, "openmc_version", VERSION);

#ifdef GIT_SHA1
  write_attribute(file_id, "git_sha1", GIT_SHA1);
#endif

  // Write current date and time
  //  write_attribute(file_id, "date_and_time", time_stamp()); <-- RE-ADD!!!
  hsize_t three = 3;
  write_attr_int(file_id, 1, &three, "num_voxels", pl->pixels);
  write_attr_double(file_id, 1, &three, "voxel_width", vox);
  write_attr_double(file_id, 1, &three, "lower_left", ll);

  // Create dataset for voxel data -- note that the dimensions are reversed
  // since we want the order in the file to be z, y, x
  hsize_t dims[3];
  dims[0] = pl->pixels[0];
  dims[1] = pl->pixels[1];
  dims[2] = pl->pixels[2];
  hid_t dspace, dset, memspace;
  voxel_init(file_id, &(dims[0]), &dspace, &dset, &memspace);

  // move to center of voxels
  ll[0] = ll[0] + vox[0] / TWO;
  ll[1] = ll[1] + vox[1] / TWO;
  ll[2] = ll[2] + vox[2] / TWO;

  int data[pl->pixels[1]][pl->pixels[2]];

  int rgb[3], id;
  for (int x = 0; x < pl->pixels[0]; x++) {
    // progress bar here
    for (int y = 0; y < pl->pixels[1]; y++) {
      for(int z = 0; z < pl->pixels[2]; z++) {
        // get voxel color
        position_rgb(p, pl, rgb, id);
        // write to plot data
        data[y][z] = id;
        // advance particle in z direction
        p->coord[0].xyz[2] = p->coord[0].xyz[2] + vox[2];
      }
      // advance particle in y direction
      p->coord[0].xyz[1] = p->coord[0].xyz[1] + vox[1];
      p->coord[0].xyz[2] = ll[2];
    }
    // advance particle in x direction
    p->coord[0].xyz[0] = p->coord[0].xyz[0] + vox[0];
    p->coord[0].xyz[1] = ll[1];
    p->coord[0].xyz[2] = ll[2];
    // Write to HDF5 dataset
    voxel_write_slice(x, dspace, dset, memspace, &(data[0]));
  }

  voxel_finalize(dspace, dset, memspace);
  file_close(file_id);

}


void
voxel_init(hid_t file_id, const hsize_t* dims, hid_t* dspace, hid_t* dset,
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

} // namespace openmc
