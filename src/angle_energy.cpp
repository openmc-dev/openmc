#include "openmc/angle_energy.h"

#include <cmath> // for sqrt

#include "openmc/random_lcg.h"

namespace openmc {

double get_jac_and_transform(
  double E_in, double& mu, double& E_out, uint64_t* seed, double awr)
{
  double E_cm = E_out;
  double mu_lab = mu;
  double D = mu_lab * mu_lab - 1.0 + (awr + 1.0) * (awr + 1.0) * E_cm / E_in;
  if (D < 0.0)
    return 0.0;
  D = std::sqrt(D);

  if ((mu_lab + D) <= 0.0)
    return 0.0;

  double E_out1 =
    E_in * ((mu_lab + D) / (awr + 1.0)) * ((mu_lab + D) / (awr + 1.0));
  double mu_cm1 =
    mu_lab * std::sqrt(E_out1 / E_cm) - std::sqrt(E_in / E_cm) / (awr + 1.0);

  double mult = 1.0;
  if (mu_lab - D <= 0.0) {
    mu = mu_cm1;
    E_out = E_out1;
  } else {
    double E_out2 =
      E_in * ((mu_lab - D) / (awr + 1.0)) * ((mu_lab - D) / (awr + 1.0));
    double mu_cm2 =
      mu_lab * std::sqrt(E_out2 / E_cm) - std::sqrt(E_in / E_cm) / (awr + 1.0);
    mult = 2.0;
    if (prn(seed) < 0.5) {
      mu = mu_cm1;
      E_out = E_out1;
    } else {
      mu = mu_cm2;
      E_out = E_out2;
    }
  }
  return mult * E_out * (awr + 1.0) / (D * std::sqrt(E_cm * E_in));
}

double get_jac_and_transform(double E_in, double& mu, double& E_out,
  uint64_t* seed, double awr, double E_com, double E_t)
{
  double E_cm = E_out;
  double mu_lab = mu;
  double D = mu_lab * mu_lab - 1.0 + E_cm / E_com;
  if (D < 0.0)
    return 0.0;
  D = std::sqrt(D);

  if ((mu_lab + D) <= 0.0)
    return 0.0;

  double E_out1 = E_com * (mu_lab + D) * (mu_lab + D);
  double mu_cm1 = (awr + 1.0) / (2 * awr * std::sqrt(E_t * E_cm)) *
                  ((awr / (awr + 1.0)) * (awr / (awr + 1.0)) * E_t + E_cm +
                    2.0 * mu_lab * std::sqrt(E_in * E_out1) - E_in - E_out1);

  double mult = 1.0;
  if (mu_lab - D <= 0.0) {
    mu = mu_cm1;
    E_out = E_out1;
  } else {
    double E_out2 = E_com * (mu_lab - D) * (mu_lab - D);
    double mu_cm2 = (awr + 1.0) / (2 * awr * std::sqrt(E_t * E_cm)) *
                    ((awr / (awr + 1.0)) * (awr / (awr + 1.0)) * E_t + E_cm +
                      2.0 * mu_lab * std::sqrt(E_in * E_out2) - E_in - E_out2);
    mult = 2.0;
    if (prn(seed) < 0.5) {
      mu = mu_cm1;
      E_out = E_out1;
    } else {
      mu = mu_cm2;
      E_out = E_out2;
    }
  }
  return mult * E_out / (D * std::sqrt(E_cm * E_com));
}

} // namespace openmc
