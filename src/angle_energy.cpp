#include "openmc/angle_energy.h"

#include <cmath> // for sqrt

#include "openmc/random_lcg.h"

namespace openmc {

double get_jac_and_transform(
  double E_in, double& mu, double& E_out, uint64_t* seed, double awr)
{
  double E_com = E_in / ((awr + 1.0) * (awr + 1.0));
  return get_jac_and_transform_impl(E_in, mu, E_out, seed, awr, E_com);
}

double get_jac_and_transform_impl(
  double E_com, double& mu, double& E_out, uint64_t* seed, double awr)
{
  double E_cm = E_out;
  double mu_lab = mu;
  double D = mu_lab * mu_lab - 1.0 + E_cm / E_com;
  if (D <= 0.0)
    return 0.0;
  D = std::sqrt(D);

  if (mu_lab <= 0.0 && E_cm <= E_cm)
    return 0.0;
  double E_out1 = E_com * (mu_lab + D) * (mu_lab + D);
  double mult;
  if (E_cm > E_com) {
    mult = 1.0;
    E_out = E_out1;
  } else {
    mult = 2.0;
    if (prn(seed) < 0.5) {
      E_out = E_out1;
    } else {
      E_out = E_com * (mu_lab - D) * (mu_lab - D);
    }
  }
  mu = mu_lab * std::sqrt(E_out / E_cm) - std::sqrt(E_com / E_cm);
  return mult * E_out / (D * std::sqrt(E_cm * E_com));
}

} // namespace openmc
