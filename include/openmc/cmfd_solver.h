#ifndef OPENMC_CMFD_SOLVER_H
#define OPENMC_CMFD_SOLVER_H

#include <cmath>

#include "xtensor/xtensor.hpp"

namespace openmc {

//===============================================================================
// Global variables
//===============================================================================

// CSR format index pointer array of loss matrix
extern std::vector<int> indptr;

// CSR format index array of loss matrix
extern std::vector<int> indices;

// Dimension n of nxn CMFD loss matrix
extern int dim;

// Spectral radius of CMFD matrices and tolerances
extern double spectral;

// Maximum dimension in x, y, and z directions
extern int nx;
extern int ny;
extern int nz;

// Number of energy groups
extern int ng;

// Indexmap storing all x, y, z positions of accelerated regions
extern xt::xtensor<int, 2> indexmap;

//===============================================================================
// Non-member functions
//===============================================================================

//! returns the index in CSR index array corresponding to the diagonal element
//! of a specified row
//! \param[in] row of interest
//! \return index in CSR index array corresponding to diagonal element
int get_diagonal_index(int row);

//! sets the elements of indexmap based on input coremap
//! \param[in] user-defined coremap
void set_indexmap(int* coremap);

//! solves a one group CMFD linear system
//! \param[in] CSR format data array of coefficient matrix
//! \param[in] right hand side vector
//! \param[out] unknown vector
//! \param[in] tolerance on final error
//! \return number of inner iterations required to reach convergence
int cmfd_linsolver_1g(double* A_data, double* b, double* x, double tol);

//! solves a two group CMFD linear system
//! \param[in] CSR format data array of coefficient matrix
//! \param[in] right hand side vector
//! \param[out] unknown vector
//! \param[in] tolerance on final error
//! \return number of inner iterations required to reach convergence
int cmfd_linsolver_2g(double* A_data, double* b, double* x, double tol);

//! solves a general CMFD linear system
//! \param[in] CSR format data array of coefficient matrix
//! \param[in] right hand side vector
//! \param[out] unknown vector
//! \param[in] tolerance on final error
//! \return number of inner iterations required to reach convergence
int cmfd_linsolver_ng(double* A_data, double* b, double* x, double tol);

//! converts a matrix index to spatial and group indices
//! \param[in] iteration counter over row
//! \param[out] iteration counter for groups
//! \param[out] iteration counter for x
//! \param[out] iteration counter for y
//! \param[out] iteration counter for z
void matrix_to_indices(int irow, int& g, int& i, int& j, int& k);

//===============================================================================
// External functions
//===============================================================================

//! sets the fixed variables that are used for the linear solver
//! \param[in] CSR format index pointer array of loss matrix
//! \param[in] length of indptr
//! \param[in] CSR format index array of loss matrix
//! \param[in] number of non-zero elements in CMFD loss matrix
//! \param[in] dimension n of nxn CMFD loss matrix
//! \param[in] spectral radius of CMFD matrices and tolerances
//! \param[in] indices storing spatial and energy dimensions of CMFD problem
//! \param[in] coremap for problem, storing accelerated regions
extern "C" void openmc_initialize_linsolver(int* indptr, int len_indptr,
                                            int* indices, int n_elements,
                                            int dim, double spectral,
                                            int* cmfd_indices, int* map);

//! runs a Gauss Seidel linear solver to solve CMFD matrix equations
//! linear solver
//! \param[in] CSR format data array of coefficient matrix
//! \param[in] right hand side vector
//! \param[out] unknown vector
//! \param[in] tolerance on final error
//! \return number of inner iterations required to reach convergence
extern "C" int openmc_run_linsolver(double* A_data, double* b, double* x,
                                    double tol);


} // namespace openmc
#endif // OPENMC_CMFD_SOLVER_H