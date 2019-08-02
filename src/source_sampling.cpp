#ifndef PYNE_IS_AMALGAMATED
#include "source_sampling.h"
#endif

// Global sampler instance
static pyne::Sampler* sampler = NULL;
// Global variable for mode range
const int SUBVOXEL_START = 3;

// Fortran API
void pyne::sampling_setup_(int* mode, int* cell_list_size) {
  if (sampler == NULL) {
    std::string filename ("source.h5m");
    std::string src_tag_name ("source_density");
    std::string e_bounds_file ("e_bounds");
    std::vector<double> e_bounds = read_e_bounds(e_bounds_file);
    std::map<std::string, std::string> tag_names;
    tag_names.insert(std::pair<std::string, std::string> ("src_tag_name",
          "source_density"));
    tag_names.insert(std::pair<std::string, std::string> ("bias_tag_name",
          "biased_source_density"));
    tag_names.insert(std::pair<std::string, std::string> ("cell_number_tag_name",
          "cell_number"));
    tag_names.insert(std::pair<std::string, std::string> ("cell_fracs_tag_name",
          "cell_fracs"));
    sampler = new pyne::Sampler(filename, tag_names, e_bounds, *mode);
    *cell_list_size = sampler->get_cell_list_size();
  }
}

void pyne::particle_birth_(double* rands,
                           double* x,
                           double* y,
                           double* z,
                           double* e,
                           double* w,
                           int* cell_list) {
    std::vector<double> rands2(rands, rands + 6);
    pyne::SourceParticle src = sampler->particle_birth(rands2);
    *x = src.get_x();
    *y = src.get_y();
    *z = src.get_z();
    *e = src.get_e();
    *w = src.get_w();
    std::vector<int> c_list = src.get_cell_list();
    int cell_list_size = sampler->get_cell_list_size();
    std::copy(c_list.begin(), c_list.end(), cell_list);
}

std::vector<double> pyne::read_e_bounds(std::string e_bounds_file){
  std::vector<double> e_bounds;
  std::ifstream inputFile(e_bounds_file.c_str());
  double value;
  if (inputFile) {
    while (inputFile >> value)
      e_bounds.push_back(value);
  }
  else {
    throw std::runtime_error("File " + e_bounds_file + " not found or no read permission");
  }
  return e_bounds;
}


// C++ API
pyne::Sampler::Sampler(std::string filename, 
                 std::string src_tag_name, 
                 std::vector<double> e_bounds, 
                 bool uniform)
  : filename(filename), src_tag_name(src_tag_name), e_bounds(e_bounds) {
  bias_mode = (uniform) ? UNIFORM : ANALOG;
  setup();
}

pyne::Sampler::Sampler(std::string filename,  
                 std::string src_tag_name, 
                 std::vector<double> e_bounds, 
                 std::string bias_tag_name)
  : filename(filename), 
    src_tag_name(src_tag_name), 
    e_bounds(e_bounds), 
    bias_tag_name(bias_tag_name) {
  bias_mode = USER;
  setup();
}

pyne::Sampler::Sampler(std::string filename,
                 std::map<std::string, std::string> tag_names,
                 std::vector<double> e_bounds, 
                 int mode)
  : filename(filename),
    tag_names(tag_names),
    e_bounds(e_bounds),
    mode(mode) {
  // determine the bias_mode
  if (mode == 0){
    bias_mode = ANALOG; 
  } else if (mode == 1) {
    bias_mode = UNIFORM;
  } else if (mode == 2) {
    bias_mode = USER;
  } else if (mode == 3) {
    bias_mode = ANALOG;
  } else if (mode == 4) {
    bias_mode = UNIFORM;
  } else if (mode == 5) {
    bias_mode = USER;
  }

  // find out the src_tag_name and bias_tag_name
  if (tag_names.find("src_tag_name") == tag_names.end()) {
    // src_tag_name not found
    throw std::invalid_argument("src_tag_name not found");
  } else {
    // found src_tag_name
    src_tag_name = tag_names["src_tag_name"];
  }
  if (bias_mode == USER) {
    // bias_tag_name required
    if (tag_names.find("bias_tag_name") == tag_names.end()) {
      // bias_tag_name not found
      throw std::invalid_argument("bias_tag_name not found");
    } else {
      // found bias_tag_name
      bias_tag_name = tag_names["bias_tag_name"];
    }
  }
  setup();
}

pyne::SourceParticle pyne::Sampler::particle_birth(std::vector<double> rands) {
  // select mesh volume and energy group
  // For Unstructured mesh, p_src_num_cells and max_num_cells are set to 1
  // For Cartisian mesh, max_num_cells is obtained via cell_fracs tag
  int pdf_idx =at->sample_pdf(rands[0], rands[1]);
  int ve_idx = pdf_idx/p_src_num_cells/num_e_groups;
  int c_idx = (pdf_idx/num_e_groups)%p_src_num_cells;
  int e_idx = pdf_idx % num_e_groups;

  // Sample uniformly within the selected mesh volume element and energy
  // group.
  std::vector<double> xyz_rands;
  xyz_rands.push_back(rands[2]);
  xyz_rands.push_back(rands[3]);
  xyz_rands.push_back(rands[4]);
  moab::CartVect pos = sample_xyz(ve_idx, xyz_rands);
  cell_list.resize(0);
  if (ve_type == moab::MBHEX) {
     if (mesh_mode == SUBVOXEL) {
         cell_list.emplace_back(cell_number[ve_idx*max_num_cells + c_idx]);
     }
     else { // Voxel
        for (int c=0; c<max_num_cells; c++) {
            cell_list.emplace_back(cell_number[ve_idx*max_num_cells + c]);
        }
     }
  }
  pyne::SourceParticle src = SourceParticle(pos[0], pos[1], pos[2],
      sample_e(e_idx, rands[5]), sample_w(pdf_idx), cell_list);
  return src;
}


void pyne::Sampler::setup() {
  moab::ErrorCode rval;
  moab::EntityHandle loaded_file_set;
  // Create MOAB instance
  mesh = new moab::Core();
  rval = mesh->create_meshset(moab::MESHSET_SET, loaded_file_set);
  rval = mesh->load_file(filename.c_str(), &loaded_file_set);
  if (rval != moab::MB_SUCCESS)
    throw std::invalid_argument("Could not load mesh file.");

  // Get mesh volume elements
  moab::Range ves;
  rval = mesh->get_entities_by_dimension(loaded_file_set, 3, ves);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem entities of dimension 3");
  num_ves = ves.size();
  int num_hex, num_tet;
  rval = mesh->get_number_entities_by_type(loaded_file_set, moab::MBHEX, num_hex);
  rval = mesh->get_number_entities_by_type(loaded_file_set, moab::MBTET, num_tet);
  if (num_hex == num_ves) {
    ve_type = moab::MBHEX;
    verts_per_ve = 8;
  } else if (num_tet == num_ves) {
    ve_type = moab::MBTET;
    verts_per_ve = 4;
  }
  else throw std::invalid_argument("Mesh file must contain only tets or hexes.");

  // Assign MeshMode: VOXEL, SUBVOXEL, TET
  // Accept mode: 0, 1, 2, 3, 4, 5
  if (ve_type == moab::MBHEX){
     if (mode < SUBVOXEL_START){
        mesh_mode = VOXEL;
     } else {
        mesh_mode = SUBVOXEL;
     }
  } else {
     mesh_mode = TET;
  }

  if (ve_type == moab::MBHEX){
      // cell_number_tag
      if (tag_names.find("cell_number_tag_name") == tag_names.end()) {
        // cell_number_tag_name not found
        throw std::invalid_argument("cell_number_tag_name not found");
      } else {
        // found cell_number_tag_name
        cell_number_tag_name = tag_names["cell_number_tag_name"];
      }
      // cell_fracs_tag
      if (tag_names.find("cell_fracs_tag_name") == tag_names.end()) {
        // cell_fracs_tag_name not found
        throw std::invalid_argument("cell_fracs_tag_name not found");
      } else {
        // found cell_fracs_tag_name
        cell_fracs_tag_name = tag_names["cell_fracs_tag_name"];
      }
  }

  // Process all the spatial and tag data and create an alias table.
  std::vector<double> volumes(num_ves);
  mesh_geom_data(ves, volumes);
  mesh_tag_data(ves, volumes);
}

void pyne::Sampler::mesh_geom_data(moab::Range ves, std::vector<double> &volumes) {
  // Get connectivity.
  moab::ErrorCode rval;
  std::vector<moab::EntityHandle> connect;
  rval = mesh->get_connectivity_by_type(ve_type, connect);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem getting mesh connectivity.");

  // Grab the coordinates that define 4 connected points within a mesh volume
  // element and setup a data structure to allow uniform sampling with each 
  // mesh volume element.
  double coords[verts_per_ve*3];
  int v;
  for (v=0; v<num_ves; ++v) {
    rval = mesh->get_coords(&connect[verts_per_ve*v], verts_per_ve, &coords[0]);
    if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem vertex coordinates.");
    volumes[v] = measure(ve_type, verts_per_ve, &coords[0]);
    if (ve_type == moab::MBHEX) {
      moab::CartVect o(coords[0], coords[1], coords[2]);
      moab::CartVect x(coords[3], coords[4], coords[5]);
      moab::CartVect y(coords[9], coords[10], coords[11]);
      moab::CartVect z(coords[12], coords[13], coords[14]);
      edge_points ep = {o, x-o, y-o, z-o};
      all_edge_points.push_back(ep);
    } else if (ve_type == moab::MBTET) {
      moab::CartVect o(coords[0], coords[1], coords[2]);
      moab::CartVect x(coords[3], coords[4], coords[5]);
      moab::CartVect y(coords[6], coords[7], coords[8]);
      moab::CartVect z(coords[9], coords[10], coords[11]);
      edge_points ep = {o, x-o, y-o, z-o};
      all_edge_points.push_back(ep);
    }
  }
}

void pyne::Sampler::mesh_tag_data(moab::Range ves, 
                                  const std::vector<double> volumes) {
  moab::ErrorCode rval;
  moab::Tag src_tag;
  moab::Tag cell_number_tag;
  moab::Tag cell_fracs_tag;
  rval = mesh->tag_get_handle(src_tag_name.c_str(),
                              moab::MB_TAG_VARLEN, 
                              moab::MB_TYPE_DOUBLE, 
                              src_tag);
  // THIS rval FAILS because we do not know number of energy groups a priori.
  // That's okay. That's what the next line is all about:
  num_e_groups = num_groups(src_tag);

  // Set the default value of max_num_cells to 1, so that the structured mesh
  // and unstructured mesh r2s can use the same form of pdf size description.
  max_num_cells = 1;
  p_src_num_cells = 1;
  int tag_size;
  // set the default value of cell_number to -1, cell_fracs to 1.0.
  cell_number.resize(num_ves, -1);
  cell_fracs.resize(num_ves, 1.0);
  if (ve_type == moab::MBHEX) {
      // Read the cell_number tag and cell_fracs tag
      rval = mesh->tag_get_handle(cell_number_tag_name.c_str(),
                                  cell_number_tag);
      rval = mesh->tag_get_handle(cell_fracs_tag_name.c_str(),
                                  cell_fracs_tag);
      has_cell_fracs = check_cell_fracs(cell_fracs_tag);
      max_num_cells = get_max_num_cells(cell_fracs_tag);
      if (mesh_mode == SUBVOXEL) {
         p_src_num_cells = max_num_cells;
         num_e_groups /= p_src_num_cells;
	 // cell_fracs must exist in SUBVOXEL mode
         if (has_cell_fracs == false) {
             throw std::runtime_error("No cell_fracs tag found in sub-voxel R2S. Wrong source file used.");
         }
         cell_fracs.resize(num_ves*p_src_num_cells);
         rval = mesh->tag_get_data(cell_fracs_tag, ves, &cell_fracs[0]);
      }
      // check and set cell_number
      if (has_cell_fracs) {
          cell_number.resize(num_ves*max_num_cells);
          rval = mesh->tag_get_data(cell_number_tag, ves, &cell_number[0]);
      }
  }
  std::vector<double> pdf(num_ves*num_e_groups*p_src_num_cells);
  rval = mesh->tag_get_data(src_tag, ves, &pdf[0]);
  if (rval != moab::MB_SUCCESS)
    throw std::runtime_error("Problem getting source tag data.");

  // Multiply the source densities by the VE volumes
  int v, c, e;
  for (v=0; v<num_ves; ++v) {
      for (c=0; c<p_src_num_cells; ++c) {
          for (e=0; e<num_e_groups; ++e) {
              pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups + e] *=
                  volumes[v]*cell_fracs[v*p_src_num_cells + c];
          }
      }
  }
  normalize_pdf(pdf);

  // Setup alias table based off PDF or biased PDF
  if (bias_mode == ANALOG) {
    at = new AliasTable(pdf);
  } else {
    std::vector<double> bias_pdf = read_bias_pdf(ves, volumes, pdf);
    normalize_pdf(bias_pdf);
    //  Create alias table based off biased pdf and calculate birth weights.
    biased_weights.resize(num_ves*p_src_num_cells*num_e_groups);
    for (int i=0; i<biased_weights.size(); ++i) {
      biased_weights[i] = pdf[i]/bias_pdf[i];
    }
    at = new AliasTable(bias_pdf);
  }
}

std::vector<double> pyne::Sampler::read_bias_pdf(moab::Range ves, 
                                                 std::vector<double> volumes,
                                                 std::vector<double> pdf) {
    std::vector<double> bias_pdf(num_ves*p_src_num_cells*num_e_groups);
    int v, c, e;
    moab::ErrorCode rval;
    if (bias_mode == UNIFORM) {
      // Sub-voxel Uniform sampling: uniform in space, analog in energy. Biased PDF is
      // found by normalizing the total photon emission density to 1 in each
      // mesh volume element and multiplying by the volume of the element.
      double q_in_group;
      for (v=0; v<num_ves; ++v) {
        for (c=0; c<p_src_num_cells; ++c) {
            q_in_group = 0.0;
            for (e=0; e<num_e_groups; ++e) {
                q_in_group += pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups + e];
            }

            if (q_in_group > 0) {
                for (e=0; e<num_e_groups; ++e) {
                    bias_pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups + e] =
                        volumes[v]*cell_fracs[v*p_src_num_cells + c]*
                        pdf[v*p_src_num_cells*num_e_groups +
                        c*num_e_groups + e]/q_in_group;
                }
            } else {
                for (e=0; e<num_e_groups; ++e) {
                  bias_pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups + e] = 0.0;
                }
            }

        }
      }
    } else if (bias_mode == USER) {
      // Get the biased PDF from the mesh
      moab::Tag bias_tag;
      rval = mesh->tag_get_handle(bias_tag_name.c_str(),
                                  moab::MB_TAG_VARLEN,
                                  moab::MB_TYPE_DOUBLE,
                                  bias_tag);
      num_bias_groups = num_groups(bias_tag);
      if (num_bias_groups == num_e_groups * p_src_num_cells) {
        // Spatial, cell and energy biasing. The supplied bias PDF values are
        // applied to each specific energy group and sub-voxels in a mesh
        // volume element.
        rval = mesh->tag_get_data(bias_tag, ves, &bias_pdf[0]);
        if (rval != moab::MB_SUCCESS)
          throw std::runtime_error("Problem getting bias tag data.");
        for (v=0; v<num_ves; ++v) {
            for (c=0; c<p_src_num_cells; c++) {
                for (e=0; e<num_e_groups; ++e)
                    bias_pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups + e] *=
                       volumes[v]*cell_fracs[v*p_src_num_cells + c];
            }
        }
      } else if (num_bias_groups == 1) {
        // Spatial biasing only: the supplied bias PDF values are applied
        // to all energy groups within a mesh volume element, which are
        // sampled in analog.
        std::vector<double> spatial_pdf(num_ves);
        rval = mesh->tag_get_data(bias_tag, ves, &spatial_pdf[0]);
        if (rval != moab::MB_SUCCESS)
          throw std::runtime_error("Problem getting bias tag data.");
        double q_in_group;
        for (v=0; v<num_ves; ++v) {
          q_in_group = 0;
          for (c=0; c<p_src_num_cells; ++c){
              for (e=0; e<num_e_groups; ++e){
                q_in_group += pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups + e];
              }
          }
          if (q_in_group > 0){
            for (c=0; c<p_src_num_cells; ++c){
                for (e=0; e<num_e_groups; ++e){
                    bias_pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups + e] =
                        spatial_pdf[v]*volumes[v]*cell_fracs[v*p_src_num_cells + c]*
                        pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups + e]/
                        q_in_group;
                }
            }
          } else {
            for (c=0; c<p_src_num_cells; ++c)
                for (e=0; e<num_e_groups; ++e){
                    bias_pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups + e] =  0;
                }
          }
        }
      } else if (num_bias_groups == num_e_groups) {
        // Voxel and energy biasing. Apply the energy bias to all the sub-voxel in the voxel
        std::vector<double> spa_erg_pdf(num_ves*num_e_groups);
        rval = mesh->tag_get_data(bias_tag, ves, &spa_erg_pdf[0]);
        if (rval != moab::MB_SUCCESS)
          throw std::runtime_error("Problem getting bias tag data.");
        double q_in_group;
        for (v=0; v<num_ves; ++v) {
            for (e=0; e<num_e_groups; ++e) {
                q_in_group = 0.0;
                for (c=0; c<p_src_num_cells; ++c) {
                    q_in_group += pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups +e];
                }
                if (q_in_group >0) {
                    for (c=0; c<p_src_num_cells; ++c) {
                        bias_pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups +e] =
                            spa_erg_pdf[v*num_e_groups+e]*volumes[v]*cell_fracs[v*p_src_num_cells + c]*
                            pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups +e]/q_in_group;
                    }
                } else {
                    for (c=0; c<p_src_num_cells; ++c) {
                        bias_pdf[v*p_src_num_cells*num_e_groups + c*num_e_groups + e] = 0.0;
                    }
                }
            }
        }
      } else {
        throw std::length_error("Length of bias tag must equal length of the"
                                "  p_src_num_cells*num_e_group, num_e_groups, or 1.");
      }
     }
    double q_in_all = 0.0;
    for (int i=0; i<bias_pdf.size(); i++)
      q_in_all += bias_pdf[i];
    if (q_in_all <= 0.0) {
      throw std::runtime_error("Bias data are ALL ZERO!");
    }
 return bias_pdf;
}


moab::CartVect pyne::Sampler::sample_xyz(int ve_idx, std::vector<double> rands) {
  double s = rands[0];
  double t = rands[1];
  double u = rands[2];

  // Transform s, t, u to uniformly sample a tetrahedron. See:
  // C. Rocchini and P. Cignoni, “Generating Random Points in a Tetrahedron,” 
  //  Journal of Graphics Tools, 5, 200–202 (2001).
  if (ve_type == moab::MBTET) {
    if (s + t > 1) {
      s = 1.0 - s;
      t = 1.0 - t;
    }
    if (s + t + u > 1) {
      if (t + u > 1) {
        double old_t = t;
        t = 1.0 - u;
        u = 1.0 - s - old_t;
      }else if (t + u <= 1) {
        double old_s = s;
        s = 1.0 - t - u;
        u = old_s + t + u - 1;
      }
    }
  }

 return s*all_edge_points[ve_idx].x_vec + \
        t*all_edge_points[ve_idx].y_vec + \
        u*all_edge_points[ve_idx].z_vec + \
          all_edge_points[ve_idx].o_point;
}

double pyne::Sampler::sample_e(int e_idx, double rand) {
   double e_min = e_bounds[e_idx];
   double e_max = e_bounds[e_idx + 1];
   return rand * (e_max - e_min) + e_min;
}

double pyne::Sampler::sample_w(int pdf_idx) {
  return (bias_mode == ANALOG) ? 1.0 : biased_weights[pdf_idx];
}

void pyne::Sampler::normalize_pdf(std::vector<double> & pdf) {
  double sum = 0;
  for (int i=0; i<pdf.size(); ++i)
    sum += pdf[i];
  for (int i=0; i<pdf.size(); ++i)
    pdf[i] /= sum;
}

int pyne::Sampler::num_groups(moab::Tag tag) {
  moab::ErrorCode rval;
  int tag_size;
  rval = mesh->tag_get_bytes(tag, *(&tag_size));
  if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem getting tag size.");
  return tag_size/sizeof(double);
}

bool pyne::Sampler::check_cell_fracs(moab::Tag cell_fracs_tag) {
  moab::ErrorCode rval;
  int tag_size;
  rval = mesh->tag_get_bytes(cell_fracs_tag, *(&tag_size));
  if (rval != moab::MB_SUCCESS) {
      std::cout<<"Warning: Old version source file used. No cell_number and cell_fracs tag found!"<<std::endl;
      std::cout<<"Default cell_number [-1] and cell_fracs [1.0] will be used."<<std::endl;
      has_cell_fracs = false;
  } else {
    has_cell_fracs = true;
  }
  return has_cell_fracs;
}

int pyne::Sampler::get_max_num_cells(moab::Tag cell_fracs_tag) {
  moab::ErrorCode rval;
  int tag_size;
  if (has_cell_fracs) {
      rval = mesh->tag_get_bytes(cell_fracs_tag, *(&tag_size));
      return tag_size/sizeof(double);
  } else {
      return 1;
  }
}

int pyne::Sampler::get_cell_list_size() {
   // cell_list_size should be:
   // 0: for unstructured mesh
   // 1: for sub-voxel R2S
   // max_num_cells: for voxel R2S
   if (ve_type == moab::MBTET or has_cell_fracs == false) {
      return 0;
   } else {
      if (mesh_mode == VOXEL) {
         return max_num_cells;
      } else {
         return 1;
      }
   }
}

// Random-number sampling using the Walker-Vose alias method,
// Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2013)
// M. D. Vose, IEEE T. Software Eng. 17, 972 (1991)
// A. J. Walker, Electronics Letters 10, 127 (1974); ACM TOMS 3, 253 (1977)

pyne::AliasTable::AliasTable(std::vector<double> p) {
  n = p.size();
  prob.resize(n);
  alias.resize(n);
  std::vector<double> small(n);
  std::vector<double> large(n);
  int i, a, g;

  for (i=0; i<n; ++i) 
    p[i] *= n;

  // Set separate index lists for small and large probabilities:
  int n_s = 0;
  int n_l = 0;
  for (i=n-1; i>=0; --i) {
    // at variance from Schwarz, we revert the index order
    if (p[i] < 1)
      small[n_s++] = i;
    else
      large[n_l++] = i;
  }

  // Work through index lists
  while(n_s && n_l) {
    a = small[--n_s]; // Schwarz's l
    g = large[--n_l]; // Schwarz's g
    prob[a] = p[a];
    alias[a] = g;
    p[g] = p[g] + p[a] - 1;
    if (p[g] < 1)
      small[n_s++] = g;
    else
      large[n_l++] = g;
  }

  while(n_l)
    prob[large[--n_l]] = 1;

  while(n_s)
    // can only happen through numeric instability
    prob[small[--n_s] ] = 1;
}

int pyne::AliasTable::sample_pdf(double rand1, double rand2) {
  int i = (int) n * rand1;
  return rand2 < prob[i] ? i : alias[i];
}

pyne::SourceParticle::SourceParticle() {
    x = -1.0;
    y = -1.0;
    z = -1.0;
    e = -1.0;
    w = -1.0;
}

pyne::SourceParticle::SourceParticle(double _x, double _y, double _z,
                                     double _e, double _w, std::vector<int> _cell_list)
    : x(_x), y(_y), z(_z),
      e(_e), w(_w), cell_list(_cell_list) {}

pyne::SourceParticle::~SourceParticle() {};
