#ifndef OPENMC_TALLIES_TALLY_H
#define OPENMC_TALLIES_TALLY_H

#include "openmc/constants.h"

#include "xtensor/xtensor.hpp"

namespace openmc {

extern "C" double total_weight;

template <std::size_t N>
using adaptor_type = xt::xtensor_adaptor<xt::xbuffer_adaptor<double*&, xt::no_ownership>, N>;

adaptor_type<2> global_tallies();
adaptor_type<3> tally_results(int idx);

#ifdef OPENMC_MPI
//! Collect all tally results onto master process
extern "C" void reduce_tally_results();
#endif

} // namespace openmc

#endif // OPENMC_TALLIES_TALLY_H
