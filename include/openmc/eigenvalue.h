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

extern double keff_generation; //!<  Single-generation k on each processor
extern array<double, 2> k_sum; //!< Used to reduce sum and sum_sq
extern vector<double> entropy; //!< Shannon entropy at each generation
extern xt::xtensor<double, 1> source_frac; //!< Source fraction for UFS

// Delayed neutron kinetics parameters
extern double keff_prompt_generation; //!< Single-generation k_prompt
extern vector<double> k_prompt;       //!< k_prompt for each generation
extern double keff_prompt;            //!< Mean k_prompt over active generations
extern double keff_prompt_std;        //!< Standard deviation of k_prompt
extern double beta_eff;               //!< Effective delayed neutron fraction
extern double beta_eff_std;           //!< Standard deviation of beta_eff
extern double alpha_k_based;          //!< Alpha eigenvalue (k-based method)
extern double alpha_k_based_std;      //!< Standard deviation of alpha (k-based)
extern double alpha_static;           //!< Alpha eigenvalue (COG Static method)
extern double alpha_static_std;    //!< Standard deviation of alpha (COG Static)
extern double prompt_gen_time;     //!< Prompt neutron generation time
extern double prompt_gen_time_std; //!< Standard deviation of gen time
extern int kinetics_tally_index;   //!< Index of internal kinetics tally

// Alpha eigenvalue calculation (COG Static method) - iteration state
extern double alpha_previous;          //!< Previous iteration's alpha value
extern double pseudo_absorption_sigma; //!< Pseudo-absorption cross section
extern int alpha_iteration;            //!< Current alpha iteration number
extern bool alpha_converged;           //!< Alpha convergence flag

} // namespace simulation

//==============================================================================
// Non-member functions
//==============================================================================

//! Collect/normalize the tracklength keff from each process
void calculate_generation_keff();

//! Collect/normalize the tracklength k_prompt from each process
void calculate_generation_prompt_keff();

//! Calculate mean/standard deviation of keff during active generations
//!
//! This function sets the global variables keff and keff_std which represent
//! the mean and standard deviation of the mean of k-effective over active
//! generations. It also broadcasts the value from the master process.
void calculate_average_keff();

//! Calculate delayed neutron kinetics parameters
//!
//! This function calculates k_prompt, beta_eff, and alpha eigenvalues
//! over active generations. Results are stored in simulation namespace.
void calculate_kinetics_parameters();

//! Setup internal tallies for alpha eigenvalue calculations
//!
//! Creates tallies with prompt chain scores needed for alpha calculations
void setup_kinetics_tallies();

//! Run alpha eigenvalue calculation (COG Static method)
//!
//! This function runs additional batches after normal eigenvalue calculation
//! to determine alpha through iterative refinement with pseudo-absorption.
//! The method seeks alpha such that K'(alpha) = 1.0, based on the alpha
//! static method implemented in the COG Monte Carlo code.
void run_alpha_iterations();

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
