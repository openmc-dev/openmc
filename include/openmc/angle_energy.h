#ifndef OPENMC_ANGLE_ENERGY_H
#define OPENMC_ANGLE_ENERGY_H

#include <cmath> // for sqrt
#include <cstdint>

#include "openmc/random_lcg.h"

namespace openmc {

//==============================================================================
//! Abstract type that defines a correlated or uncorrelated angle-energy
//! distribution that is a function of incoming energy. Each derived type must
//! implement a sample() method that returns an outgoing energy and
//! scattering cosine given an incoming energy.
//==============================================================================

class AngleEnergy {
public:
  //! Sample an outgoing energy and scattering cosine
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  virtual void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const = 0;

  //! Sample an outgoing energy and evaluate the angular PDF
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[in] mu Scattering cosine with respect to current direction
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[inout] seed Pseudorandom seed pointer
  //! \param[in] is_com Is scattering cosine is given in center of mass
  //! coordinates \param[in] awr Weight of nucleus in neutron masses \return
  //! Probability density for the scattering cosine
  virtual double sample_energy_and_pdf(double E_in, double mu, double& E_out,
    uint64_t* seed, bool is_com, double awr) const = 0;
  virtual ~AngleEnergy() = default;
};

inline double get_jac_and_transform(
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

inline double get_jac_and_transform(double E_in, double& mu, double& E_out,
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

#endif // OPENMC_ANGLE_ENERGY_H
