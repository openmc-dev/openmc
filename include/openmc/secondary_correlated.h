//! \file secondary_correlated.h
//! Correlated angle-energy distribution

#ifndef OPENMC_SECONDARY_CORRELATED_H
#define OPENMC_SECONDARY_CORRELATED_H

#include "hdf5.h"
#include "xtensor/xtensor.hpp"

#include "openmc/angle_energy.h"
#include "openmc/distribution.h"
#include "openmc/endf.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Correlated angle-energy distribution corresponding to ACE law 61 and ENDF
//! File 6, LAW=1, LANG!=2.
//==============================================================================

class CorrelatedAngleEnergy : public AngleEnergy {
public:
  //! Outgoing energy/angle at a single incoming energy
  struct CorrTable {
    int n_discrete;                    //!< Number of discrete lines
    Interpolation interpolation;       //!< Interpolation law
    xt::xtensor<double, 1> e_out;      //!< Outgoing energies [eV]
    xt::xtensor<double, 1> p;          //!< Probability density
    xt::xtensor<double, 1> c;          //!< Cumulative distribution
    vector<unique_ptr<Tabular>> angle; //!< Angle distribution
  };

  explicit CorrelatedAngleEnergy(hid_t group);

  //! Sample distribution for an angle and energy
  //! \param[in] E_in Incoming energy in [eV]
  //! \param[out] E_out Outgoing energy in [eV]
  //! \param[out] mu Outgoing cosine with respect to current direction
  //! \param[inout] seed Pseudorandom seed pointer
  void sample(
    double E_in, double& E_out, double& mu, uint64_t* seed) const override;
  Distribution& sample_dist(double E_in, double& E_out, uint64_t* seed) const;
  double sample_energy_and_pdf(
    double E_in, double mu, double& E_out, uint64_t* seed) const override;

  // energy property
  vector<double>& energy() { return energy_; }
  const vector<double>& energy() const { return energy_; }

  // distribution property
  vector<CorrTable>& distribution() { return distribution_; }
  const vector<CorrTable>& distribution() const { return distribution_; }

private:
  int n_region_;                        //!< Number of interpolation regions
  vector<int> breakpoints_;             //!< Breakpoints between regions
  vector<Interpolation> interpolation_; //!< Interpolation laws
  vector<double> energy_;          //!< Energies [eV] at which distributions
                                   //!< are tabulated
  vector<CorrTable> distribution_; //!< Distribution at each energy
};

} // namespace openmc

#endif // OPENMC_SECONDARY_CORRELATED_H
