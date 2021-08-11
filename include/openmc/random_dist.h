#ifndef OPENMC_RANDOM_DIST_H
#define OPENMC_RANDOM_DIST_H

#include <cstdint> // for uint64_t

namespace openmc {

//==============================================================================
//! Sample a uniform distribution [a, b)
//
//! \param a Lower bound of uniform distribution
//! \param b Upper bound of uniform distribtion
//! \param seed A pointer to the pseudorandom seed
//! \return Sampled variate
//==============================================================================

double uniform_distribution(double a, double b, uint64_t* seed);

//==============================================================================
//! Samples an energy from the Maxwell fission distribution based on a direct
//! sampling scheme.
//!
//! The probability distribution function for a Maxwellian is given as
//! p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T). This PDF can be sampled using
//! rule C64 in the Monte Carlo Sampler LA-9721-MS.
//!
//! \param T The tabulated function of the incoming energy
//! \param seed A pointer to the pseudorandom seed
//! \return The sampled outgoing energy
//==============================================================================

extern "C" double maxwell_spectrum(double T, uint64_t* seed);

//==============================================================================
//! Samples an energy from a Watt energy-dependent fission distribution.
//!
//! Although fitted parameters exist for many nuclides, generally the
//! continuous tabular distributions (LAW 4) should be used in lieu of the Watt
//! spectrum. This direct sampling scheme is an unpublished scheme based on the
//! original Watt spectrum derivation (See F. Brown's MC lectures).
//!
//! \param a Watt parameter a
//! \param b Watt parameter b
//! \param seed A pointer to the pseudorandom seed
//! \return The sampled outgoing energy
//==============================================================================

extern "C" double watt_spectrum(double a, double b, uint64_t* seed);

//==============================================================================
//! Samples an energy from the Gaussian energy-dependent fission distribution.
//!
//! Samples from a Normal distribution with a given mean and standard deviation
//! The PDF is defined as s(x) = (1/2*sigma*sqrt(2) * e-((mu-x)/2*sigma)^2
//! Its sampled according to
//! http://www-pdg.lbl.gov/2009/reviews/rpp2009-rev-monte-carlo-techniques.pdf
//! section 33.4.4
//!
//! \param mean mean of the Gaussian distribution
//! \param std_dev standard deviation of the Gaussian distribution
//! \param seed A pointer to the pseudorandom seed
//! \result The sampled outgoing energy
//==============================================================================

extern "C" double normal_variate(double mean, double std_dev, uint64_t* seed);

//==============================================================================
//! Samples an energy from the Muir (Gaussian) energy-dependent distribution.
//!
//! This is another form of the Gaussian distribution but with more easily
//! modifiable parameters
//! https://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-05411-MS
//!
//! \param e0 peak neutron energy [eV]
//! \param m_rat ratio of the fusion reactants to AMU
//! \param kt the ion temperature of the reactants [eV]
//! \param seed A pointer to the pseudorandom seed
//! \result The sampled outgoing energy
//==============================================================================

extern "C" double muir_spectrum(
  double e0, double m_rat, double kt, uint64_t* seed);

} // namespace openmc

#endif // OPENMC_RANDOM_DIST_H
