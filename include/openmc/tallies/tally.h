#ifndef OPENMC_TALLIES_TALLY_H
#define OPENMC_TALLIES_TALLY_H

#include "openmc/constants.h"

#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include <memory> // for unique_ptr
#include <string>
#include <vector>

namespace openmc {

//==============================================================================
//! A user-specified flux-weighted (or current) measurement.
//==============================================================================

class Tally {
public:
  Tally() {}

  std::vector<int32_t> filters_;
};

//==============================================================================
// Global variable declarations
//==============================================================================

extern "C" double total_weight;

namespace model {
  extern std::vector<std::unique_ptr<Tally>> tallies;
}

// Threadprivate variables
extern "C" double global_tally_absorption;
extern "C" double global_tally_collision;
extern "C" double global_tally_tracklength;
extern "C" double global_tally_leakage;
#pragma omp threadprivate(global_tally_absorption, global_tally_collision, \
  global_tally_tracklength, global_tally_leakage)

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

extern "C" void free_memory_tally_c();

} // namespace openmc

#endif // OPENMC_TALLIES_TALLY_H
