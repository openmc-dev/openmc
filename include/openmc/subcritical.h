//! \file subcritical.h
//! \brief Data/functions related to subcritical (fixed source) multiplication calculations

#ifndef OPENMC_SUBCRITICAL_H
#define OPENMC_SUBCRITICAL_H

#include "hdf5.h"
#include <utility>

namespace openmc {

void convert_to_subcritical_k(double& k, double& k_std);

double convert_to_subcritical_k(double k);

double calculate_ks(double k, double kq);

double calculate_sigma_ks(double k, double k_std, double kq, double kq_std);

extern "C" int openmc_get_subcritical_kq(double* k_combined);

void write_subcritical_hdf5(hid_t group);

} // namespace openmc
#endif // OPENMC_SUBCRITICAL_H