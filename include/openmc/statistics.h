//! \file statistics.h
//! A collection of helper functions for calculating sample statistics.

#ifndef OPENMC_STATISTICS_H
#define OPENMC_STATISTICS_H

#include <vector>

namespace openmc {

double get_mean(const std::vector<double>& samples);
double get_median(const std::vector<double>& samples);
double get_variance(const std::vector<double>& samples);
double get_skewness(const std::vector<double>& samples);
double get_kurtosis(const std::vector<double>& samples);

} // namespace openmc
#endif // OPENMC_STATISTICS_H
