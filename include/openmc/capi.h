#ifndef OPENMC_CAPI_H
#define OPENMC_CAPI_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int openmc_calculate_volumes();
int openmc_cell_filter_get_bins(
  int32_t index, const int32_t** cells, int32_t* n);
int openmc_cell_get_fill(
  int32_t index, int* type, int32_t** indices, int32_t* n);
int openmc_cell_get_id(int32_t index, int32_t* id);
int openmc_cell_get_temperature(
  int32_t index, const int32_t* instance, double* T);
int openmc_cell_get_translation(int32_t index, double xyz[]);
int openmc_cell_get_rotation(int32_t index, double rot[], size_t* n);
int openmc_cell_get_name(int32_t index, const char** name);
int openmc_cell_get_num_instances(int32_t index, int32_t* num_instances);
int openmc_cell_set_name(int32_t index, const char* name);
int openmc_cell_set_fill(
  int32_t index, int type, int32_t n, const int32_t* indices);
int openmc_cell_set_id(int32_t index, int32_t id);
int openmc_cell_set_temperature(
  int32_t index, double T, const int32_t* instance, bool set_contained = false);
int openmc_cell_set_translation(int32_t index, const double xyz[]);
int openmc_cell_set_rotation(int32_t index, const double rot[], size_t rot_len);
int openmc_energy_filter_get_bins(
  int32_t index, const double** energies, size_t* n);
int openmc_energy_filter_set_bins(
  int32_t index, size_t n, const double* energies);
int openmc_energyfunc_filter_get_energy(
  int32_t index, size_t* n, const double** energy);
int openmc_energyfunc_filter_get_y(int32_t index, size_t* n, const double** y);
int openmc_energyfunc_filter_set_data(
  int32_t index, size_t n, const double* energies, const double* y);
int openmc_extend_cells(int32_t n, int32_t* index_start, int32_t* index_end);
int openmc_extend_filters(int32_t n, int32_t* index_start, int32_t* index_end);
int openmc_extend_materials(
  int32_t n, int32_t* index_start, int32_t* index_end);
int openmc_extend_meshes(
  int32_t n, const char* type, int32_t* index_start, int32_t* index_end);
int openmc_extend_tallies(int32_t n, int32_t* index_start, int32_t* index_end);
int openmc_filter_get_id(int32_t index, int32_t* id);
int openmc_filter_get_type(int32_t index, char* type);
int openmc_filter_set_id(int32_t index, int32_t id);
int openmc_finalize();
int openmc_find_cell(const double* xyz, int32_t* index, int32_t* instance);
int openmc_cell_bounding_box(const int32_t index, double* llc, double* urc);
int openmc_global_bounding_box(double* llc, double* urc);
int openmc_fission_bank(void** ptr, int64_t* n);
int openmc_get_cell_index(int32_t id, int32_t* index);
int openmc_get_filter_index(int32_t id, int32_t* index);
void openmc_get_filter_next_id(int32_t* id);
int openmc_get_keff(double k_combined[]);
int openmc_get_material_index(int32_t id, int32_t* index);
int openmc_get_mesh_index(int32_t id, int32_t* index);
int openmc_get_n_batches(int* n_batches, bool get_max_batches);
int openmc_get_nuclide_index(const char name[], int* index);
int openmc_add_unstructured_mesh(
  const char filename[], const char library[], int* id);
int64_t openmc_get_seed();
int openmc_get_tally_index(int32_t id, int32_t* index);
void openmc_get_tally_next_id(int32_t* id);
int openmc_global_tallies(double** ptr);
int openmc_hard_reset();
int openmc_init(int argc, char* argv[], const void* intracomm);
bool openmc_is_statepoint_batch();
int openmc_legendre_filter_get_order(int32_t index, int* order);
int openmc_legendre_filter_set_order(int32_t index, int order);
int openmc_load_nuclide(const char* name, const double* temps, int n);
int openmc_material_add_nuclide(
  int32_t index, const char name[], double density);
int openmc_material_get_densities(
  int32_t index, const int** nuclides, const double** densities, int* n);
int openmc_material_get_id(int32_t index, int32_t* id);
int openmc_material_get_fissionable(int32_t index, bool* fissionable);
int openmc_material_get_density(int32_t index, double* density);
int openmc_material_get_volume(int32_t index, double* volume);
int openmc_material_set_density(
  int32_t index, double density, const char* units);
int openmc_material_set_densities(
  int32_t index, int n, const char** name, const double* density);
int openmc_material_set_id(int32_t index, int32_t id);
int openmc_material_get_name(int32_t index, const char** name);
int openmc_material_set_name(int32_t index, const char* name);
int openmc_material_set_volume(int32_t index, double volume);
int openmc_material_filter_get_bins(
  int32_t index, const int32_t** bins, size_t* n);
int openmc_material_filter_set_bins(
  int32_t index, size_t n, const int32_t* bins);
int openmc_mesh_filter_get_mesh(int32_t index, int32_t* index_mesh);
int openmc_mesh_filter_set_mesh(int32_t index, int32_t index_mesh);
int openmc_mesh_filter_get_translation(int32_t index, double translation[3]);
int openmc_mesh_filter_set_translation(int32_t index, double translation[3]);
int openmc_mesh_get_id(int32_t index, int32_t* id);
int openmc_mesh_set_id(int32_t index, int32_t id);
int openmc_meshsurface_filter_get_mesh(int32_t index, int32_t* index_mesh);
int openmc_meshsurface_filter_set_mesh(int32_t index, int32_t index_mesh);
int openmc_new_filter(const char* type, int32_t* index);
int openmc_next_batch(int* status);
int openmc_nuclide_name(int index, const char** name);
int openmc_plot_geometry();
int openmc_id_map(const void* slice, int32_t* data_out);
int openmc_property_map(const void* slice, double* data_out);
int openmc_rectilinear_mesh_get_grid(int32_t index, double** grid_x, int* nx,
  double** grid_y, int* ny, double** grid_z, int* nz);
int openmc_rectilinear_mesh_set_grid(int32_t index, const double* grid_x,
  const int nx, const double* grid_y, const int ny, const double* grid_z,
  const int nz);
int openmc_regular_mesh_get_dimension(int32_t index, int** id, int* n);
int openmc_regular_mesh_get_params(
  int32_t index, double** ll, double** ur, double** width, int* n);
int openmc_regular_mesh_set_dimension(int32_t index, int n, const int* dims);
int openmc_regular_mesh_set_params(int32_t index, int n, const double* ll,
  const double* ur, const double* width);
int openmc_reset();
int openmc_reset_timers();
int openmc_run();
void openmc_set_seed(int64_t new_seed);
int openmc_set_n_batches(
  int32_t n_batches, bool set_max_batches, bool add_statepoint_batch);
int openmc_simulation_finalize();
int openmc_simulation_init();
int openmc_source_bank(void** ptr, int64_t* n);
int openmc_spatial_legendre_filter_get_order(int32_t index, int* order);
int openmc_spatial_legendre_filter_get_params(
  int32_t index, int* axis, double* min, double* max);
int openmc_spatial_legendre_filter_set_order(int32_t index, int order);
int openmc_spatial_legendre_filter_set_params(
  int32_t index, const int* axis, const double* min, const double* max);
int openmc_sphharm_filter_get_order(int32_t index, int* order);
int openmc_sphharm_filter_get_cosine(int32_t index, char cosine[]);
int openmc_sphharm_filter_set_order(int32_t index, int order);
int openmc_sphharm_filter_set_cosine(int32_t index, const char cosine[]);
int openmc_statepoint_write(const char* filename, bool* write_source);
int openmc_tally_allocate(int32_t index, const char* type);
int openmc_tally_get_active(int32_t index, bool* active);
int openmc_tally_get_estimator(int32_t index, int* estimator);
int openmc_tally_get_id(int32_t index, int32_t* id);
int openmc_tally_get_filters(int32_t index, const int32_t** indices, size_t* n);
int openmc_tally_get_n_realizations(int32_t index, int32_t* n);
int openmc_tally_get_nuclides(int32_t index, int** nuclides, int* n);
int openmc_tally_get_scores(int32_t index, int** scores, int* n);
int openmc_tally_get_type(int32_t index, int32_t* type);
int openmc_tally_get_writable(int32_t index, bool* writable);
int openmc_tally_reset(int32_t index);
int openmc_tally_results(int32_t index, double** ptr, size_t shape_[3]);
int openmc_tally_set_active(int32_t index, bool active);
int openmc_tally_set_estimator(int32_t index, const char* estimator);
int openmc_tally_set_filters(int32_t index, size_t n, const int32_t* indices);
int openmc_tally_set_id(int32_t index, int32_t id);
int openmc_tally_set_nuclides(int32_t index, int n, const char** nuclides);
int openmc_tally_set_scores(int32_t index, int n, const char** scores);
int openmc_tally_set_type(int32_t index, const char* type);
int openmc_tally_set_writable(int32_t index, bool writable);
int openmc_zernike_filter_get_order(int32_t index, int* order);
int openmc_zernike_filter_get_params(
  int32_t index, double* x, double* y, double* r);
int openmc_zernike_filter_set_order(int32_t index, int order);
int openmc_zernike_filter_set_params(
  int32_t index, const double* x, const double* y, const double* r);

//! Sets the mesh and energy grid for CMFD reweight
//! \param[in] meshtyally_id id of CMFD Mesh Tally
//! \param[in] cmfd_indices indices storing spatial and energy dimensions of
//! CMFD problem \param[in] norm CMFD normalization factor
void openmc_initialize_mesh_egrid(
  const int meshtally_id, const int* cmfd_indices, const double norm);

//! Sets the mesh and energy grid for CMFD reweight
//! \param[in] feedback whether or not to run CMFD feedback
//! \param[in] cmfd_src computed CMFD source
void openmc_cmfd_reweight(const bool feedback, const double* cmfd_src);

//! Sets the fixed variables that are used for CMFD linear solver
//! \param[in] indptr CSR format index pointer array of loss matrix
//! \param[in] len_indptr length of indptr
//! \param[in] indices CSR format index array of loss matrix
//! \param[in] n_elements number of non-zero elements in CMFD loss matrix
//! \param[in] dim dimension n of nxn CMFD loss matrix
//! \param[in] spectral spectral radius of CMFD matrices and tolerances
//! \param[in] map coremap for problem, storing accelerated regions
//! \param[in] use_all_threads whether to use all threads when running CMFD
//! solver
void openmc_initialize_linsolver(const int* indptr, int len_indptr,
  const int* indices, int n_elements, int dim, double spectral, const int* map,
  bool use_all_threads);

//! Runs a Gauss Seidel linear solver to solve CMFD matrix equations
//! linear solver
//! \param[in] A_data CSR format data array of coefficient matrix
//! \param[in] b right hand side vector
//! \param[out] x unknown vector
//! \param[in] tol tolerance on final error
//! \return number of inner iterations required to reach convergence
int openmc_run_linsolver(
  const double* A_data, const double* b, double* x, double tol);

//! Export physical properties for model
//! \param[in] filename Filename to write to
//! \return Error code
int openmc_properties_export(const char* filename);

//! Import physical properties for model
//! \param[in] filename Filename to read from
// \return Error code
int openmc_properties_import(const char* filename);

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

#ifdef __cplusplus
}
#endif

#endif // OPENMC_CAPI_H
