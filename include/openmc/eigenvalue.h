//! \file eigenvalue.h
//! \brief Data/functions related to k-eigenvalue calculations

#ifndef OPENMC_EIGENVALUE_H
#define OPENMC_EIGENVALUE_H

#include <cstdint> // for int64_t

#include "xtensor/xtensor.hpp"
#include <hdf5.h>

#include "openmc/array.h"
#include "openmc/particle.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace simulation {

extern double keff_generation;     //!<  Single-generation k on each processor
extern array<double, 2> k_sum;     //!< Used to reduce sum and sum_sq
extern vector<double> entropy;     //!< Shannon entropy at each generation
extern xt::xtensor<double, 1> source_frac; //!< Source fraction for UFS

//! For alpha-eigenvalue simulation
extern array<double, 2> alpha_sum;       //!< The alpha eigenvalue
extern array<double, 2> k_alpha_sum;     //!< Multiplication factor
extern array<double, 2> rho_sum;         //!< Reactivity
extern array<double, 2> beta_sum;        //!< Delayed fission fraction
extern array<double, 2> tr_sum;          //!< Removal time (or mean life time)
extern array<double, 2> alpha_other_sum; //!< The other alpha in in-hour Eq.
//! Note: In alpha-eigenvalue mode, keff_generation and mean of k_sum converge 
//!       to one, while mean of k_alpha_sum converge to the "actual" keff that 
//!       is calculated based on the alpha-eigenfunction flux.
//!       If searching for the fundamental mode, the "other alpha" is an
//!       estimate of the left-most mode. If searching for the left-most mode,
//!       the "other alpha" is an estimate of the fundamental mode.

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

//! Collect/normalize the tracklength keff from each process
//! 
//! In alpha_mode, also collect/normalize the tracklength neutron density (Cn), 
//! prompt production (Cp), and delayed production (Cd)
void calculate_generation_keff();

//! Calculate mean/standard deviation of keff during active generations
//!
//! This function sets the global variables keff and keff_std which represent
//! the mean and standard deviation of the mean of k-effective over active
//! generations. It also broadcasts the value from the master process.
//!
//! In alpha_mode, also calculate the mean/standard deviation of alpha_eff
void calculate_average_keff();

//! Calculates a minimum variance estimate of k-effective
//!
//! The minimum variance estimate is based on a linear combination of the
//! collision, absorption, and tracklength estimates. The theory behind this can
//! be found in M. Halperin, "Almost linearly-optimum combination of unbiased
//! estimates," J. Am. Stat. Assoc., 56, 36-43 (1961),
//! doi:10.1080/01621459.1961.10482088. The implementation here follows that
//! described in T. Urbatsch et al., "Estimation and interpretation of keff
//! confidence intervals in MCNP," Nucl. Technol., 111, 169-182 (1995).
//!
//! \param[out] k_combined Estimate of k-effective and its standard deviation
//! \return Error status
extern "C" int openmc_get_keff(double* k_combined);

//! Sample/redistribute source sites from accumulated fission sites
void synchronize_bank();

//! Calculates the Shannon entropy of the fission source distribution to assess
//! source convergence
void shannon_entropy();

//! Determines the source fraction in each UFS mesh cell and reweights the
//! source bank so that the sum of the weights is equal to n_particles. The
//! 'source_frac' variable is used later to bias the production of fission sites
void ufs_count_sites();

//! Get UFS weight corresponding to particle's location
double ufs_get_weight(const Particle& p);

//! Write data related to k-eigenvalue to statepoint
//! \param[in] group HDF5 group
void write_eigenvalue_hdf5(hid_t group);

//! Read data related to k-eigenvalue from statepoint
//! \param[in] group HDF5 group
void read_eigenvalue_hdf5(hid_t group);

} // namespace openmc

#endif // OPENMC_EIGENVALUE_H
