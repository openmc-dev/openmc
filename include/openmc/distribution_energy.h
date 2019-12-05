//! \file distribution_energy.h
//! Energy distributions that depend on incident particle energy

#ifndef OPENMC_DISTRIBUTION_ENERGY_H
#define OPENMC_DISTRIBUTION_ENERGY_H

#include <vector>

#include "xtensor/xtensor.hpp"
#include "hdf5.h"

#include "openmc/constants.h"
#include "openmc/endf.h"

namespace openmc {

//===============================================================================
//! Abstract class defining an energy distribution that is a function of the
//! incident energy of a projectile. Each derived type must implement a sample()
//! function that returns a sampled outgoing energy given an incoming energy
//===============================================================================

class EnergyDistribution {
public:
  virtual double sample(double E, uint64_t* seed) const = 0;
  virtual ~EnergyDistribution() = default;
};

//===============================================================================
//! Discrete photon energy distribution
//===============================================================================

class DiscretePhoton : public EnergyDistribution {
public:
  explicit DiscretePhoton(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;
private:
  int primary_flag_; //!< Indicator of whether the photon is a primary or
                     //!< non-primary photon.
  double energy_; //!< Photon energy or binding energy
  double A_; //!< Atomic weight ratio of the target nuclide
};

//===============================================================================
//! Level inelastic scattering distribution
//===============================================================================

class LevelInelastic : public EnergyDistribution {
public:
  explicit LevelInelastic(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;
private:
  double threshold_; //!< Energy threshold in lab, (A + 1)/A * |Q|
  double mass_ratio_; //!< (A/(A+1))^2
};

//===============================================================================
//! An energy distribution represented as a tabular distribution with histogram
//! or linear-linear interpolation. This corresponds to ACE law 4, which NJOY
//! produces for a number of ENDF energy distributions.
//===============================================================================

class ContinuousTabular : public EnergyDistribution {
public:
  explicit ContinuousTabular(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;
private:
  //! Outgoing energy for a single incoming energy
  struct CTTable {
    Interpolation interpolation; //!< Interpolation law
    int n_discrete; //!< Number of of discrete energies
    xt::xtensor<double, 1> e_out; //!< Outgoing energies in [eV]
    xt::xtensor<double, 1> p; //!< Probability density
    xt::xtensor<double, 1> c; //!< Cumulative distribution
  };

  int n_region_; //!< Number of inteprolation regions
  std::vector<int> breakpoints_; //!< Breakpoints between regions
  std::vector<Interpolation> interpolation_; //!< Interpolation laws
  std::vector<double> energy_; //!< Incident energy in [eV]
  std::vector<CTTable> distribution_; //!< Distributions for each incident energy
};

//===============================================================================
//! Evaporation spectrum corresponding to ACE law 9 and ENDF File 5, LF=9.
//===============================================================================

class Evaporation : public EnergyDistribution {
public:
  explicit Evaporation(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;
private:
  Tabulated1D theta_; //!< Incoming energy dependent parameter
  double u_; //!< Restriction energy
};

//===============================================================================
//! Energy distribution of neutrons emitted from a Maxwell fission spectrum.
//! This corresponds to ACE law 7 and ENDF File 5, LF=7.
//===============================================================================

class MaxwellEnergy : public EnergyDistribution {
public:
  explicit MaxwellEnergy(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;
private:
  Tabulated1D theta_; //!< Incoming energy dependent parameter
  double u_; //!< Restriction energy
};

//===============================================================================
//! Energy distribution of neutrons emitted from a Watt fission spectrum. This
//! corresponds to ACE law 11 and ENDF File 5, LF=11.
//===============================================================================

class WattEnergy : public EnergyDistribution {
public:
  explicit WattEnergy(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const;
private:
  Tabulated1D a_; //!< Energy-dependent 'a' parameter
  Tabulated1D b_; //!< Energy-dependent 'b' parameter
  double u_; //!< Restriction energy
};

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_ENERGY_H
