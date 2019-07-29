/// \file source_sampling.h
/// \author Elliott Biondo (biondo\@wisc.edu)
///
/// \brief Mesh-based Monte Carlo source sampling.
/// 
/// The Sampler class is used for Monte Carlo source sampling from mesh-based
/// sources.  The source density distribution and optional biased source density
/// distribution are defined on a MOAB mesh. Upon instantiation, a Sampler  
/// object reads this mesh and creates an alias table for randomly sampling
/// particle birth parameters. The particle_birth member function is supplied 
/// with 6 pseudo-random numbers and returns the position, energy, and weight 
/// of a particle upon birth. 
/// There are three sampling modes: analog, uniform, and user-speficied
/// In analog sampling, no source biasing is used and birth weights
/// are all 1. In uniform sampling, the position of the particle (but not the 
/// energy) is sampled uniformly and weights are adjusted accordingly. In 
/// user-speficied mode, a supplied biased source density distribution is used 
/// for sampling and particle weights are adjusted accordingly. The biased 
/// source density distribution must have the same number of energy groups as 
/// the unbiased distribution. Alternatively, it may have exactly 1 energy
/// group, in which case only spatial biasing is done, and energies are sampled
/// in analog.
 
#ifndef PYNE_6OR6BJURKJHHTOFWXO2VMQM5EY
#define PYNE_6OR6BJURKJHHTOFWXO2VMQM5EY

#include <assert.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <stdexcept> 
#include <sstream>
#include <string>
#include <map>

#include "moab/Range.hpp"
#include "moab/Core.hpp"
#ifndef PYNE_IS_AMALGAMATED
#include "measure.h"
#endif
#include "moab/CartVect.hpp"

#ifdef __cplusplus
extern "C" {
#endif

namespace pyne {

  /// MCNP interface for source sampling setup
  /// \param mode The sampling mode: 
  /// Voxel(DEFAULT) R2S: 0 = analog, 1 = uniform, 2 = user-specified
  /// SubVoxel(SUBVOXEL) R2S: 3 = analog, 4 = uniform, 5 = user-specified
  /// cell_list_size if the variable used for create the cell_list
  /// cell_list_size could be
  ///      0: for unstructured mesh
  ///      1: for sub-voxel mode,
  ///      max_num_cells: for voxel mode
  void sampling_setup_(int* mode, int* cell_list_size);
  /// MCNP interface to sample particle birth parameters after sampling setup
  /// \param rands Six pseudo-random numbers supplied from the Fortran side.
  /// \param x The sampled x position returned by this function
  /// \param y The sampled y position returned by this function
  /// \param z The sampled z position returned by this function
  /// \param e The sampled energy returned by this function
  /// \param w The sampled y statistical weight returned by this function
  /// \param c The sampled cell number
  void particle_birth_(double* rands,
                            double* x,
                            double* y,
                            double* z,
                            double* e,
                            double* w,
                            int* cell_list);
  /// Helper function for MCNP interface that reads energy boudaries from a file
  /// \param e_bounds_file A file containing the energy group boundaries.
  std::vector<double> read_e_bounds(std::string e_bounds_file);

  /// Stores 4 connected points in a mesh volume element
  struct edge_points {
    moab::CartVect o_point;
    moab::CartVect x_vec;
    moab::CartVect y_vec;
    moab::CartVect z_vec;
  };
  
  /// A data structure for O(1) source sampling
  class AliasTable {
  public:
    /// Constructor
    /// \param p A normalized probability distribution function
    AliasTable(std::vector<double> p);
    /// Samples the alias table
    /// \param rand1 A random number in range [0, 1].
    /// \param rand2 A random number in range [0, 1].
    int sample_pdf(double rand1, double rand2);
    ~AliasTable(){};
    int n; /// Number of bins in the PDF.
    std::vector<double> prob; /// Probabilities.
    std::vector<int> alias; /// Alias probabilities.
  };

  // class Source particle
  class SourceParticle {
    public:
    SourceParticle();
    /// Constructor for source particle
    /// \param x The x coordinate of the source particle
    /// \param y The y coordinate of the source particle
    /// \param z The z coordinate of the source particle
    /// \param e The energy of the source particle
    /// \param w The weight of the source particle
    /// \param c The cell number of the source particle
    SourceParticle(double x,
                   double y,
                   double z,
                   double e,
                   double w,
                   std::vector<int> cell_list);
    ~SourceParticle();

    double get_x() {return x;};
    double get_y() {return y;};
    double get_z() {return z;};
    double get_e() {return e;};
    double get_w() {return w;};
    std::vector<int> get_cell_list() {return cell_list;};
    private:
    double x; // x coordinate
    double y; // y coordinate
    double z; // z coordinate
    double e; // energy
    double w; // weight
    /// cell list. For sub-voxel mode, the size of cell_list is 1.
    std::vector<int> cell_list;
  };
  
  /// Problem modes
  enum BiasMode {USER, ANALOG, UNIFORM};
  enum MeshMode {VOXEL, SUBVOXEL, TET};
  
  /// Mesh based Monte Carlo source sampling.
  class Sampler {
  public:
    /// Constuctor for analog and uniform sampling
    /// \param filename The path to the MOAB mesh (.h5m) file
    /// \param src_tag_name The name of the tag that describes the unbiased 
    ///                     source density distribution.
    /// \param e_bounds The energy boundaries, note there are N + 1 energy
    ///                 bounds for N energy groups
    /// \param uniform If false, analog sampling is used. If true, uniform
    ///                sampling is used.
    Sampler(std::string filename, 
            std::string src_tag_name, 
            std::vector<double> e_bounds, 
            bool uniform);
    /// Constuctor for analog and uniform sampling
    /// \param filename The path to the MOAB mesh (.h5m) file
    /// \param src_tag_name The name of the tag with the unbiased source density
    ///                     distribution.
    /// \param e_bounds The energy boundaries, note there are N + 1 energy
    ///                 bounds for N energy groups
    /// \param bias_tag_name The name of the tag describing the biased source
    ///                       density distribution. The tag must have the same
    ///                       number of energy groups as <src_tag_name> or 1.
    ///                       If 1 (i.e. spatial biasing only), all energy groups
    ///                       within a mesh volume element are sampled equally.
    Sampler(std::string filename, 
            std::string src_tag_name, 
            std::vector<double> e_bounds, 
            std::string bias_tag_name);
    /// Constuctor for overall sampler
    /// \param filename The filename of the h5m file
    /// \param tag_names The map of src_tag_name and bias_tag_name
    /// \param e_bounds The energy boundaries, note there are N + 1 energy
    ///                 bounds for N energy groups
    /// \param mode The mode number, 0, 1, 2, 3, 4 or 5
    Sampler(std::string filename,
            std::map<std::string, std::string> tag_names,
            std::vector<double> e_bounds,
            int mode);

    /// Samples particle birth parameters
    /// \param rands Six pseudo-random numbers in range [0, 1].
    /// \return A SourceParticle object containing the x position, y, position,
    ///         z, position, e, energy and w, weight of a particle.
    pyne::SourceParticle particle_birth(std::vector<double> rands);

    /// Return cell_list_size
    int get_cell_list_size();

    ~Sampler() {
      delete mesh;
      delete at;
    };
  
  // member variables
  private:
    // problem parameters
    std::string filename; ///< MOAB mesh file path
    std::string src_tag_name; ///< Unbiased source density distribution
    std::string bias_tag_name; ///< Biased source density distribution
    std::string cell_number_tag_name; ///< Cell number tag
    std::string cell_fracs_tag_name; ///< Cell volume fraction tag
    std::map<std::string, std::string> tag_names; /// < tag names
    std::vector<double> e_bounds;  ///< Energy boundaries
    int num_e_groups; ///< Number of groups in tag \a _src_tag_name
    int num_bias_groups; ///< Number of groups tag \a _bias_tag_name
    int max_num_cells; /// Max number of cells in voxels
    // Equivalent cell number in photon source.
    // For voxel R2S, p_src_num_cells = 1
    // For sub-voxel R2S, p_src_num_cells = max_num_cells
    // For unstructured R2S, p_src_num_cells = 1.
    int p_src_num_cells;
    int cell_list_size;
    bool has_cell_fracs;
    BiasMode bias_mode; ///< Bias mode: ANALOG, UNIFORM, USER
    MeshMode mesh_mode; ///< Mesh mode: VOXEL, SUBVOXEL, TET
    int mode; ///< Sampler mode, currently support 0, 1, 2, 3, 4, 5
    // mesh
    moab::Interface* mesh; ///< MOAB mesh
    int num_ves; ///< Number of mesh volume elements on \a mesh.
    moab::EntityType ve_type; ///< Type of mesh volume: moab::TET or moab::HEX
    int verts_per_ve; ///< Number of verticles per mesh volume element
    // sampling
    std::vector<edge_points> all_edge_points; ///< Four connected points on a VE.
    std::vector<double> biased_weights; ///< Birth weights for biased sampling.
    std::vector<int> cell_number; ///< Tag cell_number
    std::vector<int> cell_list; ///< Cell list passed to Fortran
    std::vector<double> cell_fracs; ///< Tag cell_fracs
    AliasTable* at; ///< Alias table used for sampling.
  
  // member functions
  private:
    // instantiation
    void setup();
    void mesh_geom_data(moab::Range ves, std::vector<double> &volumes);
    void mesh_tag_data(moab::Range ves, const std::vector<double> volumes);
    // select birth parameters
    moab::CartVect sample_xyz(int ve_idx, std::vector<double> rands);
    double sample_e(int e_idx, double rand);
    double sample_w(int pdf_idx);
    // helper functions
    void normalize_pdf(std::vector<double> & pdf);
    int num_groups(moab::Tag tag);
    std::vector<double> read_bias_pdf(moab::Range ves, std::vector<double> volumes, 
                                      std::vector<double> pdf);
    // Get max_num_cells
    int get_max_num_cells(moab::Tag cell_fracs_tag);
    // get has_cell_fracs
    bool check_cell_fracs(moab::Tag cell_fracs_tag);
  };
} //end namespace pyne

#ifdef __cplusplus
} // extern "C"
#endif

#endif // PYNE_6OR6BJURKJHHTOFWXO2VMQM5EY
