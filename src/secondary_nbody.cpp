#include "openmc/secondary_nbody.h"

#include <cmath> // for log

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/nuclide.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/tallies/tally_scoring.h"

namespace openmc {

//==============================================================================
// NBodyPhaseSpace implementation
//==============================================================================

NBodyPhaseSpace::NBodyPhaseSpace(hid_t group)
{
  read_attribute(group, "n_particles", n_bodies_);
  read_attribute(group, "total_mass", mass_ratio_);
  read_attribute(group, "atomic_weight_ratio", A_);
  read_attribute(group, "q_value", Q_);
}

void NBodyPhaseSpace::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // By definition, the distribution of the angle is isotropic for an N-body
  // phase space distribution
  mu = uniform_distribution(-1., 1., seed);

  // Determine E_max parameter
  double Ap = mass_ratio_;
  double E_max = (Ap - 1.0) / Ap * (A_ / (A_ + 1.0) * E_in + Q_);

  // x is essentially a Maxwellian distribution
  double x = maxwell_spectrum(1.0, seed);

  double y;
  double r1, r2, r3, r4, r5, r6;
  switch (n_bodies_) {
  case 3:
    y = maxwell_spectrum(1.0, seed);
    break;
  case 4:
    r1 = prn(seed);
    r2 = prn(seed);
    r3 = prn(seed);
    y = -std::log(r1 * r2 * r3);
    break;
  case 5:
    r1 = prn(seed);
    r2 = prn(seed);
    r3 = prn(seed);
    r4 = prn(seed);
    r5 = prn(seed);
    r6 = prn(seed);
    y = -std::log(r1 * r2 * r3 * r4) -
        std::log(r5) * std::pow(std::cos(PI / 2.0 * r6), 2);
    break;
  default:
    throw std::runtime_error {"N-body phase space with >5 bodies."};
  }

  // Now determine v and E_out
  double v = x / (x + y);
  E_out = E_max * v;
}

void NBodyPhaseSpace::get_pdf(double det_pos[4], double E_in, double& E_out,
  uint64_t* seed, Particle& p, std::vector<double>& mu_cm,
  std::vector<double>& Js, std::vector<Particle>& ghost_particles,
  std::vector<double>& pdfs_lab) const
{
  // By definition, the distribution of the angle is isotropic for an N-body
  // phase space distribution

  // Determine E_max parameter
  double Ap = mass_ratio_;
  double E_max = (Ap - 1.0) / Ap * (A_ / (A_ + 1.0) * E_in + Q_);

  // x is essentially a Maxwellian distribution
  double x = maxwell_spectrum(1.0, seed);

  double y;
  double r1, r2, r3, r4, r5, r6;
  switch (n_bodies_) {
  case 3:
    y = maxwell_spectrum(1.0, seed);
    break;
  case 4:
    r1 = prn(seed);
    r2 = prn(seed);
    r3 = prn(seed);
    y = -std::log(r1 * r2 * r3);
    break;
  case 5:
    r1 = prn(seed);
    r2 = prn(seed);
    r3 = prn(seed);
    r4 = prn(seed);
    r5 = prn(seed);
    r6 = prn(seed);
    y = -std::log(r1 * r2 * r3 * r4) -
        std::log(r5) * std::pow(std::cos(PI / 2.0 * r6), 2);
    break;
  default:
    throw std::runtime_error {"N-body phase space with >5 bodies."};
  }

  // Now determine v and E_out
  double v = x / (x + y);
  E_out = E_max * v;

  get_pdf_to_point_elastic(det_pos, p, mu_cm, Js, ghost_particles, E_out / 1e6);
  for (std::size_t i = 0; i < mu_cm.size(); ++i) {
    // Assuming Js.size() is the same as mu_cm.size()
    double mu_c = mu_cm[i];
    double derivative = Js[i];
    double pdf_cm = 0.5;
    pdfs_lab.push_back(pdf_cm / std::abs(derivative));
  }

  const auto& nuc {data::nuclides[p.event_nuclide()]};
  const auto& rx {nuc->reactions_[p.event_index_mt()]};
  if (!rx->scatter_in_cm_) {
    fatal_error("didn't implement lab");
  }
}

} // namespace openmc
