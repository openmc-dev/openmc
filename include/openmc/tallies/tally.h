#ifndef OPENMC_TALLIES_TALLY_H
#define OPENMC_TALLIES_TALLY_H

#include "openmc/constants.h"

#include "xtensor/xtensor.hpp"

namespace openmc {

//==============================================================================
// Global variable declarations
//==============================================================================

extern "C" double total_weight;

// Threadprivate variables

extern "C" double global_tally_absorption;
#pragma omp threadprivate(global_tally_absorption)

//==============================================================================
// Non-member functions
//==============================================================================

// Alias for the type returned by xt::adapt(...). N is the dimension of the
// multidimensional array
template <std::size_t N>
using adaptor_type = xt::xtensor_adaptor<xt::xbuffer_adaptor<double*&, xt::no_ownership>, N>;

//! Get the global tallies as a multidimensional array
//! \return Global tallies array
adaptor_type<2> global_tallies();

//! Get tally results as a multidimensional array
//! \param idx Index in tallies array
//! \return Tally results array
adaptor_type<3> tally_results(int idx);

#ifdef OPENMC_MPI
//! Collect all tally results onto master process
extern "C" void reduce_tally_results();
#endif

} // namespace openmc

#endif // OPENMC_TALLIES_TALLY_H
