//! \file subcritical.h
//! \brief Data/functions related to subcritical (fixed source) multiplication calculations

#ifndef OPENMC_SUBCRITICAL_H
#define OPENMC_SUBCRITICAL_H

#include <utility>

namespace openmc {

void convert_to_subcritical_k(double& k, double& k_std);

double convert_to_subcritical_k(double k);

} // namespace openmc
#endif // OPENMC_SUBCRITICAL_H