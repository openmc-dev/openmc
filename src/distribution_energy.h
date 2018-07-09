#ifndef OPENMC_DISTRIBUTION_ENERGY_H
#define OPENMC_DISTRIBUTION_ENERGY_H

#include <vector>

#include "xtensor/xtensor.hpp"
#include "constants.h"
#include "endf.h"
#include "hdf5.h"

namespace openmc {

//===============================================================================
//! Abstract class defining an energy distribution that is a function of the
//! incident energy of a projectile. Each derived type must implement a sample()
//! function that returns a sampled outgoing energy given an incoming energy
//===============================================================================

class EnergyDistribution {
public:
  virtual double sample(double E) const = 0;
  virtual ~EnergyDistribution() = default;
};

class DiscretePhoton : public EnergyDistribution {
public:
  explicit DiscretePhoton(hid_t group);
  double sample(double E) const;
private:
  int primary_flag_;
  double energy_;
  double A_;
};

class LevelInelastic : public EnergyDistribution {
public:
  explicit LevelInelastic(hid_t group);
  double sample(double E) const;
private:
  double threshold_;
  double mass_ratio_;
};

class ContinuousTabular : public EnergyDistribution {
public:
  explicit ContinuousTabular(hid_t group);
  double sample(double E) const;
private:
  struct CTTable {
    Interpolation interpolation;
    int n_discrete;
    xt::xtensor<double, 1> e_out;
    xt::xtensor<double, 1> p;
    xt::xtensor<double, 1> c;
  };

  int n_region_;
  std::vector<int> breakpoints_;
  std::vector<Interpolation> interpolation_;
  std::vector<double> energy_;
  std::vector<CTTable> distribution_;
};

class Evaporation : public EnergyDistribution {
public:
  explicit Evaporation(hid_t group);
  double sample(double E) const;
private:
  Tabulated1D theta_;
  double u_;
};

class MaxwellEnergy : public EnergyDistribution {
public:
  explicit MaxwellEnergy(hid_t group);
  double sample(double E) const;
private:
  Tabulated1D theta_;
  double u_;
};

class WattEnergy : public EnergyDistribution {
public:
  explicit WattEnergy(hid_t group);
  double sample(double E) const;
private:
  Tabulated1D a_;
  Tabulated1D b_;
  double u_;
};

}

#endif // OPENMC_DISTRIBUTION_ENERGY_H
