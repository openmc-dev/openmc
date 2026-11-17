//! \file distribution_energy.h
//! Energy distributions that depend on incident particle energy

#ifndef OPENMC_DISTRIBUTION_ENERGY_H
#define OPENMC_DISTRIBUTION_ENERGY_H

#include "hdf5.h"
#include "openmc/tensor.h"

#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/vector.h"

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
  double sample(double E, uint64_t* seed) const override;

private:
  friend class SecondaryDataReader;

  int primary_flag_; //!< Indicator of whether the photon is a primary or
                     //!< non-primary photon.
  double energy_;    //!< Photon energy or binding energy
  double A_;         //!< Atomic weight ratio of the target nuclide
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
  double sample(double E, uint64_t* seed) const override;

private:
  friend class SecondaryDataReader;

  double threshold_;  //!< Energy threshold in lab, (A + 1)/A * |Q|
  double mass_ratio_; //!< (A/(A+1))^2
};

//===============================================================================
//! An energy distribution represented as a tabular distribution with histogram
//! or linear-linear interpolation. This corresponds to ACE law 4, which NJOY
//! produces for a number of ENDF energy distributions.
//===============================================================================

class ContinuousTabular : public EnergyDistribution {
public:
  //! Outgoing-energy distribution at one incident energy.
  struct CTTable {
    Interpolation interpolation;
    int n_discrete;
    tensor::Tensor<double> e_out;
    tensor::Tensor<double> p;
    tensor::Tensor<double> c;
  };

  explicit ContinuousTabular(hid_t group);

  //! Sample energy distribution
  //! \param[in] E Incident particle energy in [eV]
  //! \param[inout] seed Pseudorandom number seed pointer
  //! \return Sampled energy in [eV]
  double sample(double E, uint64_t* seed) const override;

  const vector<double>& energy() const { return energy_; }
  const vector<CTTable>& distribution() const { return distribution_; }
  //! Validate and replace all tables atomically, preserving the incident grid.
  //! Discrete-only and continuous tables cannot be mixed when incident-energy
  //! interpolation uses unit-base interpolation.
  void set_spectrum(const vector<double>& energy,
    const vector<int>& interpolation, const vector<int>& n_discrete,
    const vector<vector<double>>& e_out, const vector<vector<double>>& p);

private:
  friend class SecondaryDataReader;

  static CTTable make_table(Interpolation interpolation, int n_discrete,
    const vector<double>& e_out, const vector<double>& p,
    bool native_grid = false);

  int n_region_;                        //!< Number of inteprolation regions
  vector<int> breakpoints_;             //!< Breakpoints between regions
  vector<Interpolation> interpolation_; //!< Interpolation laws
  vector<double> energy_;               //!< Incident energy in [eV]
  vector<CTTable> distribution_; //!< Distributions for each incident energy
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
  double sample(double E, uint64_t* seed) const override;

private:
  friend class SecondaryDataReader;

  Tabulated1D theta_; //!< Incoming energy dependent parameter
  double u_;          //!< Restriction energy
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
  double sample(double E, uint64_t* seed) const override;

private:
  friend class SecondaryDataReader;

  Tabulated1D theta_; //!< Incoming energy dependent parameter
  double u_;          //!< Restriction energy
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
  double sample(double E, uint64_t* seed) const override;

private:
  friend class SecondaryDataReader;

  Tabulated1D a_; //!< Energy-dependent 'a' parameter
  Tabulated1D b_; //!< Energy-dependent 'b' parameter
  double u_;      //!< Restriction energy
};

} // namespace openmc

#endif // OPENMC_DISTRIBUTION_ENERGY_H
