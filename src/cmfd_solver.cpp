#include <vector>

#include "openmc/cmfd_solver.h"
#include "openmc/error.h"
#include "openmc/constants.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

std::vector<int> indptr;

std::vector<int> indices;

int dim;

double spectral;

int nx, ny, nz, ng;

xt::xtensor<int, 2> indexmap;

//==============================================================================
// GET_DIAGONAL_INDEX returns the index in CSR index array corresponding to
// the diagonal element of a specified row
//==============================================================================

int get_diagonal_index(int row) {
  for (int j = indptr[row]; j < indptr[row+1]; j++) {
    if (indices[j] == row)
      return j;
  }

  // Return -1 if not found
  return -1;
}

//==============================================================================
// SET_INDEXMAP sets the elements of indexmap based on input coremap
//==============================================================================

void set_indexmap(int* coremap) {
  for (int z = 0; z < nz; z++) {
    for (int y = 0; y < ny; y++) {
      for (int x = 0; x < nx; x++) {
        if (coremap[(z*ny*nx) + (y*nx) + x] != CMFD_NOACCEL) {
          int counter = coremap[(z*ny*nx) + (y*nx) + x];
          indexmap(counter, 0) = x;
          indexmap(counter, 1) = y;
          indexmap(counter, 2) = z;
        }
      }
    }
  }
}

//==============================================================================
// CMFD_LINSOLVER_1G solves a one group CMFD linear system
//==============================================================================

int cmfd_linsolver_1g(double* A_data, double* b, double* x, double tol) {
  // Set overrelaxation parameter
  double w = 1.0;

  // Perform Gauss-Seidel iterations
  for (int igs = 1; igs <= 10000; igs++) {
    double tmpx[dim];
    double err = 0.0;

    // Copy over x vector
    std::copy(x, x+dim, tmpx);

    // Perform red/black Gauss-Seidel iterations
    for (int irb = 0; irb < 2; irb++) {

      // Loop around matrix rows
      for (int irow = 0; irow < dim; irow++) {
        int g, i, j, k;
        matrix_to_indices(irow, g, i, j, k);

        // Filter out black cells
        if ((i+j+k) % 2 != irb) continue;

        // Get index of diagonal for current row
        int didx = get_diagonal_index(irow);

        // Perform temporary sums, first do left of diag, then right of diag
        double tmp1 = 0.0;
        for (int icol = indptr[irow]; icol < didx; icol++)
          tmp1 += A_data[icol] * x[indices[icol]];
        for (int icol = didx + 1; icol < indptr[irow + 1]; icol++)
          tmp1 += A_data[icol] * x[indices[icol]];

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
    err = std::sqrt(err / dim);
    if (err < tol)
      return igs;

    // Calculate new overrelaxation parameter
    w = 1.0/(1.0 - 0.25 * spectral * w);
  }

  // Throw error, as max iterations met
  fatal_error("Maximum Gauss-Seidel iterations encountered.");

  // Return -1 by default, although error thrown before reaching this point
  return -1;
}

//==============================================================================
// CMFD_LINSOLVER_2G solves a two group CMFD linear system
//==============================================================================

int cmfd_linsolver_2g(double* A_data, double* b, double* x, double tol) {
  // Set overrelaxation parameter
  double w = 1.0;

  // Perform Gauss-Seidel iterations
  for (int igs = 1; igs <= 10000; igs++) {
    double tmpx[dim];
    double err = 0.0;

    // Copy over x vector
    std::copy(x, x+dim, tmpx);

    // Perform red/black Gauss-Seidel iterations
    for (int irb = 0; irb < 2; irb++) {

      // Loop around matrix rows
      for (int irow = 0; irow < dim; irow+=2) {
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
        for (int icol = indptr[irow]; icol < d1idx; icol++)
          tmp1 += A_data[icol] * x[indices[icol]];
        for (int icol = indptr[irow+1]; icol < d2idx-1; icol++)
          tmp2 += A_data[icol] * x[indices[icol]];
        for (int icol = d1idx + 2; icol < indptr[irow + 1]; icol++)
          tmp1 += A_data[icol] * x[indices[icol]];
        for (int icol = d2idx + 1; icol < indptr[irow + 2]; icol++)
          tmp2 += A_data[icol] * x[indices[icol]];

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
    err = std::sqrt(err / dim);
    if (err < tol)
      return igs;

    // Calculate new overrelaxation parameter
    w = 1.0/(1.0 - 0.25 * spectral * w);
  }

  // Throw error, as max iterations met
  fatal_error("Maximum Gauss-Seidel iterations encountered.");

  // Return -1 by default, although error thrown before reaching this point
  return -1;
}

//==============================================================================
// CMFD_LINSOLVER_NG solves a general CMFD linear system
//==============================================================================

int cmfd_linsolver_ng(double* A_data, double* b, double* x, double tol) {
  // Set overrelaxation parameter
  double w = 1.0;

  // Perform Gauss-Seidel iterations
  for (int igs = 1; igs <= 10000; igs++) {
    double tmpx[dim];
    double err = 0.0;

    // Copy over x vector
    std::copy(x, x+dim, tmpx);

    // Loop around matrix rows
    for (int irow = 0; irow < dim; irow++) {
      // Get index of diagonal for current row
      int didx = get_diagonal_index(irow);

      // Perform temporary sums, first do left of diag, then right of diag
      double tmp1 = 0.0;
      for (int icol = indptr[irow]; icol < didx; icol++)
        tmp1 += A_data[icol] * x[indices[icol]];
      for (int icol = didx + 1; icol < indptr[irow + 1]; icol++)
        tmp1 += A_data[icol] * x[indices[icol]];

      // Solve for new x
      double x1 = (b[irow] - tmp1) / A_data[didx];

      // Perform overrelaxation
      x[irow] = (1.0 - w) * x[irow] + w * x1;

      // Compute residual and update error
      double res = (tmpx[irow] - x[irow]) / tmpx[irow];
      err += res * res;
    }

    // Check convergence
    err = std::sqrt(err / dim);
    if (err < tol)
      return igs;

    // Calculate new overrelaxation parameter
    w = 1.0/(1.0 - 0.25 * spectral * w);
  }

  // Throw error, as max iterations met
  fatal_error("Maximum Gauss-Seidel iterations encountered.");

  // Return -1 by default, although error thrown before reaching this point
  return -1;
}

//==============================================================================
// MATRIX_TO_INDICES converts a matrix index to spatial and group
// indices
//==============================================================================

void matrix_to_indices(int irow, int& g, int& i, int& j, int& k) {
  g = irow % ng;
  i = indexmap(irow/ng, 0);
  j = indexmap(irow/ng, 1);
  k = indexmap(irow/ng, 2);
}

//==============================================================================
// OPENMC_INITIALIZE_LINSOLVER sets the fixed variables that are used for the
// linear solver
//==============================================================================

extern "C"
void openmc_initialize_linsolver(int* indptr, int len_indptr, int* indices,
                                 int n_elements, int dim, double spectral,
                                 int* cmfd_indices, int* map) {
  // Store elements of indptr
  for (int i = 0; i < len_indptr; i++)
    openmc::indptr.push_back(indptr[i]);

  // Store elements of indices
  for (int i = 0; i < n_elements; i++)
    openmc::indices.push_back(indices[i]);

  // Set dimenion of CMFD problem and specral radius
  openmc::dim = dim;
  openmc::spectral = spectral;

  // Set number of groups
  openmc::ng = cmfd_indices[3];

  // Set problem dimensions and indexmap if 1 or 2 group problem
  if (openmc::ng == 1 || openmc::ng == 2) {
    openmc::nx = cmfd_indices[0];
    openmc::ny = cmfd_indices[1];
    openmc::nz = cmfd_indices[2];

    // Resize indexmap and set its elements
    openmc::indexmap.resize({static_cast<size_t>(dim), 3});
    set_indexmap(map);
  }
}

//==============================================================================
// OPENMC_RUN_LINSOLVER runs a Gauss Seidel linear solver to solve CMFD matrix
// equations
//==============================================================================

extern "C"
int openmc_run_linsolver(double* A_data, double* b, double* x, double tol) {
  switch (ng) {
  case 1:
    return cmfd_linsolver_1g(A_data, b, x, tol);
  case 2:
    return cmfd_linsolver_2g(A_data, b, x, tol);
  default:
    return cmfd_linsolver_ng(A_data, b, x, tol);
  }
}


} // namespace openmc