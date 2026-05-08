#ifndef OPENMC_ANGLE_ENERGY_H
#define OPENMC_ANGLE_ENERGY_H

#include <cstdint>

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

double get_jac_and_transform(
  double E_in, double& mu, double& E_out, uint64_t* seed, double awr);

double get_jac_and_transform(double E_in, double& mu, double& E_out,
  uint64_t* seed, double awr, double E_com, double E_t);

} // namespace openmc

#endif // OPENMC_ANGLE_ENERGY_H
