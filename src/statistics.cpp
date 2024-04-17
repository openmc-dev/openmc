#include "openmc/statistics.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace openmc {

double get_mean(const std::vector<double>& samples)
{
  const int n = samples.size();
  double sum = 0.0;
  for (double val : samples) {
    sum += val;
  }

  return sum / n;
}

double get_median(const std::vector<double>& samples)
{
  double median;
  const int n = samples.size();

  // Sort samples
  std::vector<double> samples_sorted = samples;
  std::sort(samples_sorted.begin(), samples_sorted.begin() + n);

  if (n % 2 == 0) {
    median = (samples_sorted[n / 2 - 1] + samples_sorted[n / 2]) / 2.0;
  } else {
    median = samples_sorted[n / 2];
  }
  return median;
}

double get_variance(const std::vector<double>& samples)
{
  const int n = samples.size();
  double sum = 0.0;
  double sum_sq = 0.0;
  for (double val : samples) {
    sum += val;
    sum_sq += val * val;
  }
  const double mean = sum / n;

  return (sum_sq / n - mean * mean) / (n - 1);
}

double get_skewness(const std::vector<double>& samples)
{
  // Based on G1 in https://www.jstor.org/stable/2988433

  const double n = samples.size();

  // Mean
  double mean = 0.0;
  for (double val : samples) {
    mean += val;
  }
  mean /= n;

  // Sample moments
  double m2 = 0.0, m3 = 0.0;
  for (double x : samples) {
    const double val = x - mean;
    const double val_sq = val * val;
    m2 += val_sq;
    m3 += val_sq * val;
  }
  m2 /= n;
  m3 /= n;

  // Skewness
  const double g1 = m3 / std::pow(m2, 1.5);
  return std::sqrt(n * (n - 1)) / (n - 2) * g1;
}

double get_kurtosis(const std::vector<double>& samples)
{
  // Based on G2 in https://www.jstor.org/stable/2988433
  const double n = samples.size();

  // Mean
  double mean = 0.0;
  for (double val : samples) {
    mean += val;
  }
  mean /= n;

  // Sample moments
  double m2 = 0.0, m4 = 0.0;
  for (double x : samples) {
    const double val = x - mean;
    const double val_sq = val * val;
    m2 += val_sq;
    m4 += val_sq * val_sq;
  }
  m2 /= n;
  m4 /= n;

  // Kurtosis
  const double g2 = m4 / (m2 * m2) - 3.0;
  return (n - 1) / ((n - 2) * (n - 3)) * ((n + 1) * g2 + 6);
}

} // namespace openmc
