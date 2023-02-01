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

double* device_ttb_e_grid {nullptr};
size_t ttb_e_grid_size {0};

} // namespace data

//==============================================================================
// Bremsstrahlung implementation
//==============================================================================

void BremsstrahlungData::copy_to_device()
{
  device_pdf_ = pdf_.data();
  device_cdf_ = cdf_.data();
  device_yield_ = yield_.data();
  #pragma omp target enter data map(to: device_pdf_[:pdf_.size()])
  #pragma omp target enter data map(to: device_cdf_[:cdf_.size()])
  #pragma omp target enter data map(to: device_yield_[:yield_.size()])
}

void BremsstrahlungData::release_from_device()
{
  #pragma omp target exit data map(release: device_pdf_[:pdf_.size()])
  #pragma omp target exit data map(release: device_cdf_[:cdf_.size()])
  #pragma omp target exit data map(release: device_yield_[:yield_.size()])
}

double BremsstrahlungData::pdf(gsl::index i, gsl::index j) const
{
  return *(device_pdf_ + i * data::ttb_e_grid_size + j);
}

double BremsstrahlungData::cdf(gsl::index i, gsl::index j) const
{
  return *(device_cdf_ + i * data::ttb_e_grid_size + j);
}

void Bremsstrahlung::copy_to_device()
{
  electron.copy_to_device();
  positron.copy_to_device();
}

void Bremsstrahlung::release_from_device()
{
  electron.release_from_device();
  positron.release_from_device();
}

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
    mat = &model::materials[p.material_].ttb_.positron;
  } else {
    mat = &model::materials[p.material_].ttb_.electron;
  }

  double e = std::log(p.E_);
  auto n_e = data::ttb_e_grid_size;

  // Find the lower bounding index of the incident electron energy
  size_t j = lower_bound_index(data::device_ttb_e_grid,
    data::device_ttb_e_grid + data::ttb_e_grid_size, e);
  if (j == n_e - 1) --j;

  // Get the interpolation bounds
  double e_l = data::device_ttb_e_grid[j];
  double e_r = data::device_ttb_e_grid[j+1];
  double y_l = mat->device_yield_[j];
  double y_r = mat->device_yield_[j+1];

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
    auto start = mat->device_cdf_ + i_e * data::ttb_e_grid_size;
    int i_w = lower_bound_index(start, start + i_e, c);

    // Sample the photon energy
    double w_l = data::device_ttb_e_grid[i_w];
    double w_r = data::device_ttb_e_grid[i_w + 1];
    double p_l = mat->pdf(i_e, i_w);
    double p_r = mat->pdf(i_e, i_w + 1);
    double c_l = mat->cdf(i_e, i_w);
    double a = std::log(p_r/p_l)/(w_r - w_l) + 1.0;
    double w = std::exp(w_l)*std::pow(a*(c - c_l)/(std::exp(w_l)*p_l) + 1.0, 1.0/a);

    if (w > settings::energy_cutoff[photon]) {
      // Create secondary photon
      p.create_secondary(p.wgt_, p.u(), w, Particle::Type::photon);
      *E_lost += w;
    }
  }
}

} // namespace openmc
