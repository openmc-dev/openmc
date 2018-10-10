#ifndef OPENMC_TALLIES_TALLY_H
#define OPENMC_TALLIES_TALLY_H

#include "openmc/constants.h"

#include "xtensor/xtensor.hpp"

namespace openmc {

extern "C" double total_weight;

xt::xtensor<double, 3> tally_results(int idx);

xt::xtensor<double, 2> global_tallies();

#ifdef OPENMC_MPI
//! Collect all tally results onto master process
extern "C" void reduce_tally_results();
#endif

} // namespace openmc

#endif // OPENMC_TALLIES_TALLY_H
