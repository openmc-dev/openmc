#include "openmc/bremsstrahlung.h"

#include "openmc/constants.h"
#include "openmc/material.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/settings.h"

#include "xtensor/xmath.hpp"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {

xt::xtensor<double, 1> ttb_e_grid;
xt::xtensor<double, 1> ttb_k_grid;
std::vector<Bremsstrahlung> ttb;

} // namespace data

//==============================================================================
// Non-member functions
//==============================================================================

void thick_target_bremsstrahlung(Particle& p, double* E_lost)
{
  if (p.material_ == MATERIAL_VOID) return;

  int photon = static_cast<int>(Particle::Type::photon);
  if (p.E_ < settings::energy_cutoff[photon]) return;

  // Get bremsstrahlung data for this material and particle type
  BremsstrahlungData* mat;
  if (p.type_ == Particle::Type::positron) {
    mat = &model::materials[p.material_]->ttb_->positron;
  } else {
    mat = &model::materials[p.material_]->ttb_->electron;
  }

  double e = std::log(p.E_);
  auto n_e = data::ttb_e_grid.size();

  // Find the lower bounding index of the incident electron energy
  size_t j = lower_bound_index(data::ttb_e_grid.cbegin(),
    data::ttb_e_grid.cend(), e);
  if (j == n_e - 1) --j;

  // Get the interpolation bounds
  double e_l = data::ttb_e_grid(j);
  double e_r = data::ttb_e_grid(j+1);
  double y_l = mat->yield(j);
  double y_r = mat->yield(j+1);

  // Calculate the interpolation weight w_j+1 of the bremsstrahlung energy PDF
  // interpolated in log energy, which can be interpreted as the probability
  // of index j+1
  double f = (e - e_l)/(e_r - e_l);

  // Get the photon number yield for the given energy using linear
  // interpolation on a log-log scale
  double y = std::exp(y_l + (y_r - y_l)*f);

  // Sample number of secondary bremsstrahlung photons
  int n = y + prn(p.current_seed());

  *E_lost = 0.0;
  if (n == 0) return;

  // Sample index of the tabulated PDF in the energy grid, j or j+1
  double c_max;
  int i_e;
  if (prn(p.current_seed()) <= f || j == 0) {
    i_e = j + 1;

    // Interpolate the maximum value of the CDF at the incoming particle
    // energy on a log-log scale
    double p_l = mat->pdf(i_e, i_e - 1);
    double p_r = mat->pdf(i_e, i_e);
    double c_l = mat->cdf(i_e, i_e - 1);
    double a = std::log(p_r/p_l)/(e_r - e_l) + 1.0;
    c_max = c_l + std::exp(e_l)*p_l/a*(std::exp(a*(e - e_l)) - 1.0);
  } else {
    i_e = j;

    // Maximum value of the CDF
    c_max = mat->cdf(i_e, i_e);
  }

  // Sample the energies of the emitted photons
  for (int i = 0; i < n; ++i) {
    // Generate a random number r and determine the index i for which
    // cdf(i) <= r*cdf,max <= cdf(i+1)
    double c = prn(p.current_seed())*c_max;
    int i_w = lower_bound_index(&mat->cdf(i_e, 0), &mat->cdf(i_e, 0) + i_e, c);

    // Sample the photon energy
    double w_l = data::ttb_e_grid(i_w);
    double w_r = data::ttb_e_grid(i_w + 1);
    double p_l = mat->pdf(i_e, i_w);
    double p_r = mat->pdf(i_e, i_w + 1);
    double c_l = mat->cdf(i_e, i_w);
    double a = std::log(p_r/p_l)/(w_r - w_l) + 1.0;
    double w = std::exp(w_l)*std::pow(a*(c - c_l)/(std::exp(w_l)*p_l) + 1.0, 1.0/a);

    if (w > settings::energy_cutoff[photon]) {
      // Create secondary photon
      p.create_secondary(p.u(), w, Particle::Type::photon);
      *E_lost += w;
    }
  }
}

} // namespace openmc
