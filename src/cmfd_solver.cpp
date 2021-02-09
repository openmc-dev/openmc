#include "openmc/cmfd_solver.h"

#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "xtensor/xtensor.hpp"

#include "openmc/bank.h"
#include "openmc/error.h"
#include "openmc/constants.h"
#include "openmc/capi.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/tallies/filter_energy.h"
#include "openmc/tallies/filter_mesh.h"
#include "openmc/tallies/tally.h"

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

StructuredMesh* mesh;

std::vector<double> egrid;

double norm;

} // namespace cmfd

//==============================================================================
// GET_CMFD_ENERGY_BIN returns the energy bin for a source site energy
//==============================================================================

int get_cmfd_energy_bin(const double E)
{
  // Check if energy is out of grid bounds
  if (E < cmfd::egrid[0]) {
    // throw warning message
    warning("Detected source point below energy grid");
    return 0;
  } else if (E >= cmfd::egrid[cmfd::ng]) {
    // throw warning message
    warning("Detected source point above energy grid");
    return cmfd::ng - 1;
  } else {
    // Iterate through energy grid to find matching bin
    for (int g = 0; g < cmfd::ng; g++) {
      if (E >= cmfd::egrid[g] && E < cmfd::egrid[g+1]) {
        return g;
      }
    }
  }
  // Return -1 by default
  return -1;
}

//==============================================================================
// COUNT_BANK_SITES bins fission sites according to CMFD mesh and energy
//==============================================================================

xt::xtensor<double, 1> count_bank_sites(xt::xtensor<int, 1>& bins, bool* outside)
{
  // Determine shape of array for counts
  std::size_t cnt_size = cmfd::nx * cmfd::ny * cmfd::nz * cmfd::ng;
  std::vector<std::size_t> cnt_shape = {cnt_size};

  // Create array of zeros
  xt::xarray<double> cnt {cnt_shape, 0.0};
  bool outside_ = false;

  auto bank_size = simulation::source_bank.size();
  for (int i = 0; i < bank_size; i++) {
    const auto& site = simulation::source_bank[i];

    // determine scoring bin for CMFD mesh
    int mesh_bin = cmfd::mesh->get_bin(site.r);

    // if outside mesh, skip particle
    if (mesh_bin < 0) {
      outside_ = true;
      continue;
    }

    // determine scoring bin for CMFD energy
    int energy_bin = get_cmfd_energy_bin(site.E);

    // add to appropriate bin
    cnt(mesh_bin*cmfd::ng+energy_bin) += site.wgt;

    // store bin index which is used again when updating weights
    bins[i] = mesh_bin*cmfd::ng+energy_bin;
  }

  // Create copy of count data. Since ownership will be acquired by xtensor,
  // std::allocator must be used to avoid Valgrind mismatched free() / delete
  // warnings.
  int total = cnt.size();
  double* cnt_reduced = std::allocator<double>{}.allocate(total);

#ifdef OPENMC_MPI
  // collect values from all processors
  MPI_Reduce(cnt.data(), cnt_reduced, total, MPI_DOUBLE, MPI_SUM, 0,
    mpi::intracomm);

  // Check if there were sites outside the mesh for any processor
  MPI_Reduce(&outside_, outside, 1, MPI_C_BOOL, MPI_LOR, 0, mpi::intracomm);

#else
  std::copy(cnt.data(), cnt.data() + total, cnt_reduced);
  *outside = outside_;
#endif

  // Adapt reduced values in array back into an xarray
  auto arr = xt::adapt(cnt_reduced, total, xt::acquire_ownership(), cnt_shape);
  xt::xarray<double> counts = arr;

  return counts;
}

//==============================================================================
// OPENMC_CMFD_REWEIGHT performs reweighting of particles in source bank
//==============================================================================

extern "C"
void openmc_cmfd_reweight(const bool feedback, const double* cmfd_src)
{
  // Get size of source bank and cmfd_src
  auto bank_size = simulation::source_bank.size();
  std::size_t src_size = cmfd::nx * cmfd::ny * cmfd::nz * cmfd::ng;

  // count bank sites for CMFD mesh, store bins in bank_bins for reweighting
  xt::xtensor<int, 1> bank_bins({bank_size}, 0);
  bool sites_outside;
  xt::xtensor<double, 1> sourcecounts = count_bank_sites(bank_bins,
                                                         &sites_outside);

  // Compute CMFD weightfactors
  xt::xtensor<double, 1> weightfactors = xt::xtensor<double, 1>({src_size}, 1.);
  if (mpi::master) {
    if (sites_outside) {
      fatal_error("Source sites outside of the CMFD mesh");
    }

    double norm = xt::sum(sourcecounts)()/cmfd::norm;
    for (int i = 0; i < src_size; i++) {
      if (sourcecounts[i] > 0 && cmfd_src[i] > 0) {
        weightfactors[i] = cmfd_src[i] * norm / sourcecounts[i];
      }
    }
  }

  if (!feedback) return;

#ifdef OPENMC_MPI
  // Send weightfactors to all processors
  MPI_Bcast(weightfactors.data(), src_size, MPI_DOUBLE, 0, mpi::intracomm);
#endif

  // Iterate through fission bank and update particle weights
  for (int64_t i = 0; i < bank_size; i++) {
    auto& site = simulation::source_bank[i];
    site.wgt *= weightfactors(bank_bins(i));
  }
}

//==============================================================================
// OPENMC_INITIALIZE_MESH_EGRID sets the mesh and energy grid for CMFD reweight
//==============================================================================

extern "C"
void openmc_initialize_mesh_egrid(const int meshtally_id, const int* cmfd_indices,
                                  const double norm)
{
  // Make sure all CMFD memory is freed
  free_memory_cmfd();

  // Set CMFD indices
  cmfd::nx = cmfd_indices[0];
  cmfd::ny = cmfd_indices[1];
  cmfd::nz = cmfd_indices[2];
  cmfd::ng = cmfd_indices[3];

  // Set CMFD reweight properties
  cmfd::norm = norm;

  // Find index corresponding to tally id
  int32_t tally_index;
  openmc_get_tally_index(meshtally_id, &tally_index);

  // Get filters assocaited with tally
  const auto& tally_filters = model::tallies[tally_index]->filters();

  // Get mesh filter index
  auto meshfilter_index = tally_filters[0];

  // Store energy filter index if defined, otherwise set to -1
  auto energy_index = (tally_filters.size() == 2) ? tally_filters[1] : -1;

  // Get mesh index from mesh filter index
  int32_t mesh_index;
  openmc_mesh_filter_get_mesh(meshfilter_index, &mesh_index);

  // Get mesh from mesh index
  cmfd::mesh = dynamic_cast<StructuredMesh*>(model::meshes[mesh_index].get());

  // Get energy bins from energy index, otherwise use default
  if (energy_index != -1) {
    auto efilt_base = model::tally_filters[energy_index].get();
    auto* efilt = dynamic_cast<EnergyFilter*>(efilt_base);
    cmfd::egrid = efilt->bins();
  } else {
    cmfd::egrid = {0.0, INFTY};
  }
}

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
        int idx = (z*cmfd::ny*cmfd::nx) + (y*cmfd::nx) + x;
        if (coremap[idx] != CMFD_NOACCEL) {
          int counter = coremap[idx];
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
                                 double spectral, const int* map,
                                 bool use_all_threads)
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

  // Set indexmap if 1 or 2 group problem
  if (cmfd::ng == 1 || cmfd::ng == 2) {
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
  // Clear std::vectors
  cmfd::indptr.clear();
  cmfd::indices.clear();
  cmfd::egrid.clear();

  // Resize xtensors to be empty
  cmfd::indexmap.resize({0});

  // Set pointers to null
  cmfd::mesh = nullptr;
}

} // namespace openmc
