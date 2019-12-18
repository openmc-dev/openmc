#include "openmc/cmfd_solver.h"

#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "xtensor/xtensor.hpp"

#include "openmc/error.h"
#include "openmc/constants.h"
#include "openmc/capi.h"

namespace openmc {

namespace cmfd {

//==============================================================================
// Global variables
//==============================================================================

std::vector<int> indptr;

std::vector<int> indices;

int dim;

double spectral;

int nx, ny, nz, ng;

xt::xtensor<int, 2> indexmap;

int use_all_threads;

} // namespace cmfd

//==============================================================================
// MATRIX_TO_INDICES converts a matrix index to spatial and group
// indices
//==============================================================================

void matrix_to_indices(int irow, int& g, int& i, int& j, int& k)
{
  g = irow % cmfd::ng;
  i = cmfd::indexmap(irow/cmfd::ng, 0);
  j = cmfd::indexmap(irow/cmfd::ng, 1);
  k = cmfd::indexmap(irow/cmfd::ng, 2);
}

//==============================================================================
// GET_DIAGONAL_INDEX returns the index in CSR index array corresponding to
// the diagonal element of a specified row
//==============================================================================

int get_diagonal_index(int row)
{
  for (int j = cmfd::indptr[row]; j < cmfd::indptr[row+1]; j++) {
    if (cmfd::indices[j] == row)
      return j;
  }

  // Return -1 if not found
  return -1;
}

//==============================================================================
// SET_INDEXMAP sets the elements of indexmap based on input coremap
//==============================================================================

void set_indexmap(const int* coremap)
{
  for (int z = 0; z < cmfd::nz; z++) {
    for (int y = 0; y < cmfd::ny; y++) {
      for (int x = 0; x < cmfd::nx; x++) {
        if (coremap[(z*cmfd::ny*cmfd::nx) + (y*cmfd::nx) + x] != CMFD_NOACCEL) {
          int counter = coremap[(z*cmfd::ny*cmfd::nx) + (y*cmfd::nx) + x];
          cmfd::indexmap(counter, 0) = x;
          cmfd::indexmap(counter, 1) = y;
          cmfd::indexmap(counter, 2) = z;
        }
      }
    }
  }
}

//==============================================================================
// CMFD_LINSOLVER_1G solves a one group CMFD linear system
//==============================================================================

int cmfd_linsolver_1g(const double* A_data, const double* b, double* x,
                      double tol)
{
  // Set overrelaxation parameter
  double w = 1.0;

  // Perform Gauss-Seidel iterations
  for (int igs = 1; igs <= 10000; igs++) {
    double err = 0.0;

    // Copy over x vector
    std::vector<double> tmpx {x, x+cmfd::dim};

    // Perform red/black Gauss-Seidel iterations
    for (int irb = 0; irb < 2; irb++) {

      // Loop around matrix rows
      #pragma omp parallel for reduction (+:err) if(cmfd::use_all_threads)
      for (int irow = 0; irow < cmfd::dim; irow++) {
        int g, i, j, k;
        matrix_to_indices(irow, g, i, j, k);

        // Filter out black cells
        if ((i+j+k) % 2 != irb) continue;

        // Get index of diagonal for current row
        int didx = get_diagonal_index(irow);

        // Perform temporary sums, first do left of diag, then right of diag
        double tmp1 = 0.0;
        for (int icol = cmfd::indptr[irow]; icol < didx; icol++)
          tmp1 += A_data[icol] * x[cmfd::indices[icol]];
        for (int icol = didx + 1; icol < cmfd::indptr[irow + 1]; icol++)
          tmp1 += A_data[icol] * x[cmfd::indices[icol]];

        // Solve for new x
        double x1 = (b[irow] - tmp1) / A_data[didx];

        // Perform overrelaxation
        x[irow] = (1.0 - w) * x[irow] + w * x1;

        // Compute residual and update error
        double res = (tmpx[irow] - x[irow]) / tmpx[irow];
        err += res * res;
      }
    }

    // Check convergence
    err = std::sqrt(err / cmfd::dim);
    if (err < tol)
      return igs;

    // Calculate new overrelaxation parameter
    w = 1.0/(1.0 - 0.25 * cmfd::spectral * w);
  }

  // Throw error, as max iterations met
  fatal_error("Maximum Gauss-Seidel iterations encountered.");

  // Return -1 by default, although error thrown before reaching this point
  return -1;
}

//==============================================================================
// CMFD_LINSOLVER_2G solves a two group CMFD linear system
//==============================================================================

int cmfd_linsolver_2g(const double* A_data, const double* b, double* x,
                      double tol)
{
  // Set overrelaxation parameter
  double w = 1.0;

  // Perform Gauss-Seidel iterations
  for (int igs = 1; igs <= 10000; igs++) {
    double err = 0.0;

    // Copy over x vector
    std::vector<double> tmpx {x, x+cmfd::dim};

    // Perform red/black Gauss-Seidel iterations
    for (int irb = 0; irb < 2; irb++) {

      // Loop around matrix rows
      #pragma omp parallel for reduction (+:err) if(cmfd::use_all_threads) 
      for (int irow = 0; irow < cmfd::dim; irow+=2) {
        int g, i, j, k;
        matrix_to_indices(irow, g, i, j, k);

        // Filter out black cells
        if ((i+j+k) % 2 != irb) continue;

        // Get index of diagonals for current row and next row
        int d1idx = get_diagonal_index(irow);
        int d2idx = get_diagonal_index(irow+1);

        // Get block diagonal
        double m11 = A_data[d1idx];     // group 1 diagonal
        double m12 = A_data[d1idx + 1]; // group 1 right of diagonal (sorted by col)
        double m21 = A_data[d2idx - 1]; // group 2 left of diagonal (sorted by col)
        double m22 = A_data[d2idx];     // group 2 diagonal

        // Analytically invert the diagonal
        double dm = m11*m22 - m12*m21;
        double d11 = m22/dm;
        double d12 = -m12/dm;
        double d21 = -m21/dm;
        double d22 = m11/dm;

        // Perform temporary sums, first do left of diag, then right of diag
        double tmp1 = 0.0;
        double tmp2 = 0.0;
        for (int icol = cmfd::indptr[irow]; icol < d1idx; icol++)
          tmp1 += A_data[icol] * x[cmfd::indices[icol]];
        for (int icol = cmfd::indptr[irow+1]; icol < d2idx-1; icol++)
          tmp2 += A_data[icol] * x[cmfd::indices[icol]];
        for (int icol = d1idx + 2; icol < cmfd::indptr[irow + 1]; icol++)
          tmp1 += A_data[icol] * x[cmfd::indices[icol]];
        for (int icol = d2idx + 1; icol < cmfd::indptr[irow + 2]; icol++)
          tmp2 += A_data[icol] * x[cmfd::indices[icol]];

        // Adjust with RHS vector
        tmp1 = b[irow] - tmp1;
        tmp2 = b[irow + 1] - tmp2;

        // Solve for new x
        double x1 = d11*tmp1 + d12*tmp2;
        double x2 = d21*tmp1 + d22*tmp2;

        // Perform overrelaxation
        x[irow] = (1.0 - w) * x[irow] + w * x1;
        x[irow + 1] = (1.0 - w) * x[irow + 1] + w * x2;

        // Compute residual and update error
        double res = (tmpx[irow] - x[irow]) / tmpx[irow];
        err += res * res;
      }
    }

    // Check convergence
    err = std::sqrt(err / cmfd::dim);
    if (err < tol)
      return igs;

    // Calculate new overrelaxation parameter
    w = 1.0/(1.0 - 0.25 * cmfd::spectral * w);
  }

  // Throw error, as max iterations met
  fatal_error("Maximum Gauss-Seidel iterations encountered.");

  // Return -1 by default, although error thrown before reaching this point
  return -1;
}

//==============================================================================
// CMFD_LINSOLVER_NG solves a general CMFD linear system
//==============================================================================

int cmfd_linsolver_ng(const double* A_data, const double* b, double* x,
                      double tol)
{
  // Set overrelaxation parameter
  double w = 1.0;

  // Perform Gauss-Seidel iterations
  for (int igs = 1; igs <= 10000; igs++) {
    double err = 0.0;

    // Copy over x vector
    std::vector<double> tmpx {x, x+cmfd::dim};

    // Loop around matrix rows
    for (int irow = 0; irow < cmfd::dim; irow++) {
      // Get index of diagonal for current row
      int didx = get_diagonal_index(irow);

      // Perform temporary sums, first do left of diag, then right of diag
      double tmp1 = 0.0;
      for (int icol = cmfd::indptr[irow]; icol < didx; icol++)
        tmp1 += A_data[icol] * x[cmfd::indices[icol]];
      for (int icol = didx + 1; icol < cmfd::indptr[irow + 1]; icol++)
        tmp1 += A_data[icol] * x[cmfd::indices[icol]];

      // Solve for new x
      double x1 = (b[irow] - tmp1) / A_data[didx];

      // Perform overrelaxation
      x[irow] = (1.0 - w) * x[irow] + w * x1;

      // Compute residual and update error
      double res = (tmpx[irow] - x[irow]) / tmpx[irow];
      err += res * res;
    }

    // Check convergence
    err = std::sqrt(err / cmfd::dim);
    if (err < tol)
      return igs;

    // Calculate new overrelaxation parameter
    w = 1.0/(1.0 - 0.25 * cmfd::spectral * w);
  }

  // Throw error, as max iterations met
  fatal_error("Maximum Gauss-Seidel iterations encountered.");

  // Return -1 by default, although error thrown before reaching this point
  return -1;
}

//==============================================================================
// OPENMC_INITIALIZE_LINSOLVER sets the fixed variables that are used for the
// linear solver
//==============================================================================

extern "C"
void openmc_initialize_linsolver(const int* indptr, int len_indptr,
                                 const int* indices, int n_elements, int dim,
                                 double spectral, const int* cmfd_indices,
                                 const int* map, bool use_all_threads)
{
  // Store elements of indptr
  for (int i = 0; i < len_indptr; i++)
    cmfd::indptr.push_back(indptr[i]);

  // Store elements of indices
  for (int i = 0; i < n_elements; i++)
    cmfd::indices.push_back(indices[i]);

  // Set dimenion of CMFD problem and specral radius
  cmfd::dim = dim;
  cmfd::spectral = spectral;

  // Set number of groups
  cmfd::ng = cmfd_indices[3];

  // Set problem dimensions and indexmap if 1 or 2 group problem
  if (cmfd::ng == 1 || cmfd::ng == 2) {
    cmfd::nx = cmfd_indices[0];
    cmfd::ny = cmfd_indices[1];
    cmfd::nz = cmfd_indices[2];

    // Resize indexmap and set its elements
    cmfd::indexmap.resize({static_cast<size_t>(dim), 3});
    set_indexmap(map);
  }

  // Use all threads allocated to OpenMC simulation to run CMFD solver
  cmfd::use_all_threads = use_all_threads;
}

//==============================================================================
// OPENMC_RUN_LINSOLVER runs a Gauss Seidel linear solver to solve CMFD matrix
// equations
//==============================================================================

extern "C"
int openmc_run_linsolver(const double* A_data, const double* b, double* x,
                         double tol)
{
  switch (cmfd::ng) {
  case 1:
    return cmfd_linsolver_1g(A_data, b, x, tol);
  case 2:
    return cmfd_linsolver_2g(A_data, b, x, tol);
  default:
    return cmfd_linsolver_ng(A_data, b, x, tol);
  }
}

void free_memory_cmfd()
{
  cmfd::indptr.clear();
  cmfd::indices.clear();
  // Resize indexmap to be an empty array
  cmfd::indexmap.resize({0});
}

} // namespace openmc
