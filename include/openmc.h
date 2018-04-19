#ifndef OPENMC_H
#define OPENMC_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

  struct Bank {
    double wgt;
    double xyz[3];
    double uvw[3];
    double E;
    int delayed_group;
  };

  void openmc_calculate_volumes();
  int openmc_cell_get_fill(int32_t index, int* type, int32_t** indices, int32_t* n);
  int openmc_cell_get_id(int32_t index, int32_t* id);
  int openmc_cell_set_fill(int32_t index, int type, int32_t n, const int32_t* indices);
  int openmc_cell_set_id(int32_t index, int32_t id);
  int openmc_cell_set_temperature(int32_t index, double T, const int32_t* instance);
  int openmc_energy_filter_get_bins(int32_t index, double** energies, int32_t* n);
  int openmc_energy_filter_set_bins(int32_t index, int32_t n, const double* energies);
  int openmc_extend_cells(int32_t n, int32_t* index_start, int32_t* index_end);
  int openmc_extend_filters(int32_t n, int32_t* index_start, int32_t* index_end);
  int openmc_extend_materials(int32_t n, int32_t* index_start, int32_t* index_end);
  int openmc_extend_meshes(int32_t n, int32_t* index_start, int32_t* index_end);
  int openmc_extend_sources(int32_t n, int32_t* index_start, int32_t* index_end);
  int openmc_extend_tallies(int32_t n, int32_t* index_start, int32_t* index_end);
  int openmc_filter_get_id(int32_t index, int32_t* id);
  int openmc_filter_get_type(int32_t index, const char** type);
  int openmc_filter_set_id(int32_t index, int32_t id);
  int openmc_filter_set_type(int32_t index, const char* type);
  void openmc_finalize();
  int openmc_find(double* xyz, int rtype, int32_t* id, int32_t* instance);
  int openmc_get_cell_index(int32_t id, int32_t* index);
  int openmc_get_filter_index(int32_t id, int32_t* index);
  void openmc_get_filter_next_id(int32_t* id);
  int openmc_get_keff(double k_combined[]);
  int openmc_get_material_index(int32_t id, int32_t* index);
  int openmc_get_mesh_index(int32_t id, int32_t* index);
  int openmc_get_nuclide_index(const char name[], int* index);
  int64_t openmc_get_seed();
  int openmc_get_tally_index(int32_t id, int32_t* index);
  void openmc_hard_reset();
  void openmc_init(const int* intracomm);
  int openmc_legendre_filter_get_order(int32_t index, int* order);
  int openmc_legendre_filter_set_order(int32_t index, int order);
  int openmc_load_nuclide(char name[]);
  int openmc_material_add_nuclide(int32_t index, const char name[], double density);
  int openmc_material_get_densities(int32_t index, int** nuclides, double** densities, int* n);
  int openmc_material_get_id(int32_t index, int32_t* id);
  int openmc_material_set_density(int32_t index, double density);
  int openmc_material_set_densities(int32_t index, int n, const char** name, const double* density);
  int openmc_material_set_id(int32_t index, int32_t id);
  int openmc_material_filter_get_bins(int32_t index, int32_t** bins, int32_t* n);
  int openmc_material_filter_set_bins(int32_t index, int32_t n, const int32_t* bins);
  int openmc_mesh_filter_get_mesh(int32_t index, int32_t* index_mesh);
  int openmc_mesh_filter_set_mesh(int32_t index, int32_t index_mesh);
  int openmc_mesh_get_id(int32_t index, int32_t* id);
  int openmc_mesh_get_dimension(int32_t index, int** id, int* n);
  int openmc_mesh_get_params(int32_t index, double** ll, double** ur, double** width, int* n);
  int openmc_mesh_set_id(int32_t index, int32_t id);
  int openmc_mesh_set_dimension(int32_t index, int n, const int* dims);
  int openmc_mesh_set_params(int32_t index, const double* ll, const double* ur, const double* width, int n);
  int openmc_meshsurface_filter_get_mesh(int32_t index, int32_t* index_mesh);
  int openmc_meshsurface_filter_set_mesh(int32_t index, int32_t index_mesh);
  int openmc_next_batch(int* status);
  int openmc_nuclide_name(int index, char** name);
  void openmc_plot_geometry();
  void openmc_reset();
  int openmc_run();
  void openmc_set_seed(int64_t new_seed);
  void openmc_simulation_finalize();
  void openmc_simulation_init();
  int openmc_source_bank(struct Bank** ptr, int64_t* n);
  int openmc_source_set_strength(int32_t index, double strength);
  int openmc_spatial_legendre_filter_get_order(int32_t index, int* order);
  int openmc_spatial_legendre_filter_get_params(int32_t index, int* axis, double* min, double* max);
  int openmc_spatial_legendre_filter_set_order(int32_t index, int order);
  int openmc_spatial_legendre_filter_set_params(int32_t index, const int* axis,
                                                const double* min, const double* max);
  int openmc_sphharm_filter_get_order(int32_t index, int* order);
  int openmc_sphharm_filter_get_cosine(int32_t index, char cosine[]);
  int openmc_sphharm_filter_set_order(int32_t index, int order);
  int openmc_sphharm_filter_set_cosine(int32_t index, const char cosine[]);
  void openmc_statepoint_write(const char filename[]);
  int openmc_tally_get_active(int32_t index, bool* active);
  int openmc_tally_get_id(int32_t index, int32_t* id);
  int openmc_tally_get_filters(int32_t index, int32_t** indices, int* n);
  int openmc_tally_get_n_realizations(int32_t index, int32_t* n);
  int openmc_tally_get_nuclides(int32_t index, int** nuclides, int* n);
  int openmc_tally_get_scores(int32_t index, int** scores, int* n);
  int openmc_tally_results(int32_t index, double** ptr, int shape_[3]);
  int openmc_tally_set_active(int32_t index, bool active);
  int openmc_tally_set_filters(int32_t index, int n, const int32_t* indices);
  int openmc_tally_set_id(int32_t index, int32_t id);
  int openmc_tally_set_nuclides(int32_t index, int n, const char** nuclides);
  int openmc_tally_set_scores(int32_t index, int n, const char** scores);
  int openmc_tally_set_type(int32_t index, const char* type);
  int openmc_zernike_filter_get_order(int32_t index, int* order);
  int openmc_zernike_filter_get_params(int32_t index, double* x, double* y, double* r);
  int openmc_zernike_filter_set_order(int32_t index, int order);
  int openmc_zernike_filter_set_params(int32_t index, const double* x,
                                       const double* y, const double* r);

  // Error codes
  extern int OPENMC_E_UNASSIGNED;
  extern int OPENMC_E_ALLOCATE;
  extern int OPENMC_E_OUT_OF_BOUNDS;
  extern int OPENMC_E_INVALID_SIZE;
  extern int OPENMC_E_INVALID_ARGUMENT;
  extern int OPENMC_E_INVALID_TYPE;
  extern int OPENMC_E_INVALID_ID;
  extern int OPENMC_E_GEOMETRY;
  extern int OPENMC_E_DATA;
  extern int OPENMC_E_PHYSICS;
  extern int OPENMC_E_WARNING;

  // Global variables
  extern char openmc_err_msg[256];
  extern double keff;
  extern double keff_std;
  extern int32_t n_batches;
  extern int32_t n_cells;
  extern int32_t n_filters;
  extern int32_t n_inactive;
  extern int32_t n_lattices;
  extern int32_t n_materials;
  extern int32_t n_meshes;
  extern int64_t n_particles;
  extern int32_t n_plots;
  extern int32_t n_realizations;
  extern int32_t n_sab_tables;
  extern int32_t n_sources;
  extern int32_t n_surfaces;
  extern int32_t n_tallies;
  extern int32_t n_universes;
  extern int run_mode;
  extern bool simulation_initialized;
  extern int verbosity;

#ifdef __cplusplus
}
#endif

#endif // OPENMC_H
