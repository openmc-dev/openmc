//! \file secondary_kalbach.h
//! Kalbach-Mann angle-energy distribution

#ifndef OPENMC_SECONDARY_KALBACH_H
#define OPENMC_SECONDARY_KALBACH_H

#include "hdf5.h"
#include "xtensor/xtensor.hpp"

#include "openmc/angle_energy.h"
#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Correlated angle-energy distribution with the angular distribution
//! represented using Kalbach-Mann systematics. This corresponds to ACE law 44
//! and ENDF File 6, LAW=1, LANG=2.
//==============================================================================

class KalbachMann : public AngleEnergy {
public:
  explicit KalbachMann(hid_t group);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const override;

private:
  //! Outgoing energy/angle at a single incoming energy
  struct KMTable {
    int n_discrete;               //!< Number of discrete lines
    Interpolation interpolation;  //!< Interpolation law
    xt::xtensor<double, 1> e_out; //!< Outgoing energies [eV]
    xt::xtensor<double, 1> p;     //!< Probability density
    xt::xtensor<double, 1> c;     //!< Cumulative distribution
    xt::xtensor<double, 1> r;     //!< Pre-compound fraction
    xt::xtensor<double, 1> a;     //!< Parameterized function
  };

  bool is_photon_;                      //!< Whether the projectile is a photon
  int n_region_;                        //!< Number of interpolation regions
  vector<int> breakpoints_;             //!< Breakpoints between regions
  vector<Interpolation> interpolation_; //!< Interpolation laws
  vector<double> energy_;        //!< Energies [eV] at which distributions
                                 //!< are tabulated
  vector<KMTable> distribution_; //!< Distribution at each energy
};

} // namespace openmc

#endif // OPENMC_SECONDARY_KALBACH_H
