#include "openmc/random_dist.h"

#include <cmath>

#include "openmc/constants.h"
#include "openmc/random_lcg.h"

namespace openmc {

double uniform_distribution(double a, double b, uint64_t* seed)
{
  return a + (b - a) * prn(seed);
}

double maxwell_spectrum(double T, uint64_t* seed)
{
  // Set the random numbers
  double r1 = prn(seed);
  double r2 = prn(seed);
  double r3 = prn(seed);

  // determine cosine of pi/2*r
  double c = std::cos(PI / 2. * r3);

  // Determine outgoing energy
  return -T * (std::log(r1) + std::log(r2) * c * c);
}

double watt_spectrum(double a, double b, uint64_t* seed)
{
  double w = maxwell_spectrum(a, seed);
  return w + 0.25 * a * a * b +
         uniform_distribution(-1., 1., seed) * std::sqrt(a * a * b * w);
}

double normal_variate(double mean, double standard_deviation, uint64_t* seed)
{
  // Sample a normal variate using Marsaglia's polar method
  double x, y, r2;
  do {
    x = uniform_distribution(-1., 1., seed);
    y = uniform_distribution(-1., 1., seed);
    r2 = x * x + y * y;
  } while (r2 > 1 || r2 == 0);
  double z = std::sqrt(-2.0 * std::log(r2) / r2);
  return mean + standard_deviation * z * x;
}

double muir_spectrum(double e0, double m_rat, double kt, uint64_t* seed)
{
  // https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-05411-MS
  double sigma = std::sqrt(4. * e0 * kt / m_rat);
  return normal_variate(e0, sigma, seed);
}

} // namespace openmc
