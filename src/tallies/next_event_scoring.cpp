#include "openmc/tallies/next_event_scoring.h"

#include <algorithm> // for clamp

#include "openmc/secondary_uncorrelated.h"
#include "openmc/source.h"

namespace openmc {

void score_point_tally_elastic(
  Particle& p, int i_nuclide, const Reaction& rx, int i_product, Direction v_t)
{

  const auto& nuc {data::nuclides[i_nuclide]};
  double awr = nuc->awr_;

  // Neutron velocity in LAB
  Direction v_n = std::sqrt(p.E()) * p.u();
  auto u_n = v_n / v_n.norm();

  // Velocity of center-of-mass
  Direction v_cm = (v_n + awr * v_t) / (awr + 1.0);
  auto u_cm = v_cm / v_cm.norm();

  double E_in = p.E();
  double E_com = v_cm.dot(v_cm);
  double E_out = (v_n - v_cm).dot(v_n - v_cm);

  auto& d = rx.products_[i_product].distribution_[0];
  auto d_ = dynamic_cast<UncorrelatedAngleEnergy*>(d.get());

  auto pdf = [&](Direction u, double& E) {
    double mu = u.dot(u_cm);
    double mu_l = u.dot(u_n);
    E = E_out;
    double jac =
      get_jac_and_transform_impl(E_com, mu, E, p.current_seed(), awr);
    double mu_cm =
      1.0 + mu_l * std::sqrt(E_in * E) / E_out - (E_in + E) / (2.0 * E_out);
    mu_cm = std::clamp(mu_cm, -1.0, 1.0);
    if (!d_->angle().empty()) {
      return jac * d_->angle().evaluate(p.E(), mu_cm) / (2.0 * PI);
    } else {
      return jac * 0.5 / (2.0 * PI);
    }
  };
  score_point_tally_impl(p.r(), p.type(), p.time(), pdf);
}

void score_point_tally_inelastic(
  Particle& p, int i_nuclide, const Reaction& rx, int i_product, double yield)
{
  const auto& nuc {data::nuclides[i_nuclide]};
  double awr = nuc->awr_;
  auto u_n = p.u();
  auto E_n = p.E();
  auto is_com = rx.scatter_in_cm_;

  auto pdf = [&](Direction u, double& E) {
    double mu = u.dot(u_n);
    return rx.products_[i_product].sample_energy_and_pdf(
             E_n, mu, E, p.current_seed(), is_com, awr) /
           (2.0 * PI) * yield;
  };
  score_point_tally_impl(p.r(), p.type(), p.time(), pdf);
}

void score_point_tally_fission(
  Particle& p, int i_nuclide, const Reaction& rx, int i_product)
{
  const auto& nuc {data::nuclides[i_nuclide]};
  double awr = nuc->awr_;
  auto u_n = p.u();
  auto E_n = p.E();
  auto is_com = rx.scatter_in_cm_;

  auto pdf = [&](Direction u, double& E) {
    double mu = u.dot(u_n);
    return rx.products_[i_product].sample_energy_and_pdf(
             E_n, mu, E, p.current_seed(), is_com, awr) /
           (2.0 * PI);
  };
  score_point_tally_impl(p.r(), p.type(), p.time(), pdf);
}

void score_point_tally_sab(Particle& p, int i_nuclide, const ThermalData& sab,
  const NuclideMicroXS& micro)
{
  const auto& nuc {data::nuclides[i_nuclide]};
  double awr = nuc->awr_;
  auto u_n = p.u();
  auto E_n = p.E();
  auto pdf = [&](Direction u, double& E) {
    double mu = u.dot(u_n);
    return sab.sample_energy_and_pdf(
             micro, E_n, mu, E, p.current_seed(), false, awr) /
           (2.0 * PI);
  };
  score_point_tally_impl(p.r(), p.type(), p.time(), pdf);
}

void score_point_tally_source(SourceSite& site, int source_index)
{
  auto src_ = model::external_sources[source_index].get();
  auto src = dynamic_cast<IndependentSource*>(src_);
  if (!src)
    fatal_error("Only independent source is valid for point detectors.");
  auto pdf = [&](Direction u, double& E) {
    E = site.E;
    return src->angle()->evaluate(u);
  };
  score_point_tally_impl(site.r, site.particle, site.time, pdf);
}

} // namespace openmc
