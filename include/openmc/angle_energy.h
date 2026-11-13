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

//! Transform a lab-frame scattering cosine and center-of-mass outgoing energy
//! for a stationary target, and return the Jacobian of the transformation
//!
//! Computes the energy of the center-of-mass motion from E_in and awr and
//! calls get_jac_and_transform_impl.
//!
//! \param[in] E_in Incoming energy in the lab frame in [eV]
//! \param[in,out] mu On input, the scattering cosine in the lab frame. On
//!   output, the corresponding scattering cosine in the center-of-mass frame.
//! \param[in,out] E_out On input, the outgoing energy in the center-of-mass
//!   frame in [eV]. On output, the outgoing energy in the lab frame in [eV].
//! \param[inout] seed Pseudorandom seed pointer. A random number is consumed
//!   only when two lab-frame outgoing energies are possible.
//! \param[in] awr Weight of nucleus in neutron masses
//! \return Jacobian that converts a center-of-mass angular PDF into a lab-frame
//!   one (including a factor of 2 when two lab energies are possible), or 0 if
//!   the lab cosine cannot be reached from the given center-of-mass energy, in
//!   which case mu and E_out are left unchanged
double get_jac_and_transform(
  double E_in, double& mu, double& E_out, uint64_t* seed, double awr);

//! Transform a lab-frame scattering cosine and center-of-mass outgoing energy
//! given the energy of the center-of-mass motion, and return the Jacobian of
//! the transformation
//!
//! \param[in] E_com Kinetic energy in the lab frame of a neutron moving with
//!   the center-of-mass velocity in [eV]
//! \param[in,out] mu On input, the scattering cosine in the lab frame. On
//!   output, the corresponding scattering cosine in the center-of-mass frame.
//! \param[in,out] E_out On input, the outgoing energy in the center-of-mass
//!   frame in [eV]. On output, the outgoing energy in the lab frame in [eV].
//! \param[inout] seed Pseudorandom seed pointer. A random number is consumed
//!   only when two lab-frame outgoing energies are possible.
//! \param[in] awr Weight of nucleus in neutron masses
//! \return Jacobian that converts a center-of-mass angular PDF into a lab-frame
//!   one (including a factor of 2 when two lab energies are possible), or 0 if
//!   the lab cosine cannot be reached from the given center-of-mass energy, in
//!   which case mu and E_out are left unchanged
double get_jac_and_transform_impl(
  double E_com, double& mu, double& E_out, uint64_t* seed, double awr);

} // namespace openmc

#endif // OPENMC_ANGLE_ENERGY_H
