#include "openmc/physics_mg.h"

#include <stdexcept>

#include "xtensor/xarray.hpp"
#include <fmt/core.h>

#include "openmc/bank.h"
#include "openmc/constants.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/math_functions.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/particle.h"
#include "openmc/physics_common.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/tallies/tally.h"
#include "openmc/weight_windows.h"

namespace openmc {

void collision_mg(Particle& p)
{
  // Add to the collision counter for the particle
  p.n_collision()++;

  if (settings::weight_window_checkpoint_collision)
    apply_weight_windows(p);

  // Sample the reaction type
  sample_reaction(p);

  // Display information about collision
  if ((settings::verbosity >= 10) || p.trace()) {
    write_message(fmt::format("    Energy Group = {}", p.g()), 1);
  }
}

void sample_reaction(Particle& p)
{
  // Create fission bank sites. Note that while a fission reaction is sampled,
  // it never actually "happens", i.e. the weight of the particle does not
  // change when sampling fission sites. The following block handles all
  // absorption (including fission)

  if (model::materials[p.material()]->fissionable()) {
    if (settings::run_mode == RunMode::EIGENVALUE ||
        (settings::run_mode == RunMode::FIXED_SOURCE &&
          settings::create_fission_neutrons)) {
      create_fission_sites(p);
    }
  }

  // If survival biasing is being used, the following subroutine adjusts the
  // weight of the particle. Otherwise, it checks to see if absorption occurs.
  if (p.macro_xs().absorption > 0.) {
    absorption(p);
  }
  if (!p.alive())
    return;

  // Sample a scattering event to determine the energy of the exiting neutron
  scatter(p);

  // Play Russian roulette if survival biasing is turned on
  if (settings::survival_biasing) {
    if (p.wgt() < settings::weight_cutoff) {
      russian_roulette(p, settings::weight_survive);
    }
  }
}

void scatter(Particle& p)
{
  data::mg.macro_xs_[p.material()].sample_scatter(p.g_last(), p.g(), p.mu(),
    p.wgt(), p.current_seed(), p.mg_xs_cache().t, p.mg_xs_cache().a);

  // Rotate the angle
  p.u() = rotate_angle(p.u(), p.mu(), nullptr, p.current_seed());

  // Update energy value for downstream compatability (in tallying)
  p.E() = data::mg.energy_bin_avg_[p.g()];

  // Set event component
  p.event() = TallyEvent::SCATTER;
}

void create_fission_sites(Particle& p)
{
  // If uniform fission source weighting is turned on, we increase or decrease
  // the expected number of fission sites produced
  double weight = settings::ufs_on ? ufs_get_weight(p) : 1.0;

  // Determine the expected number of neutrons produced
  double nu_t = p.wgt() / simulation::keff * weight * p.macro_xs().nu_fission /
                p.macro_xs().total;

  // Sample the number of neutrons produced
  int nu = static_cast<int>(nu_t);
  if (prn(p.current_seed()) <= (nu_t - int(nu_t))) {
    nu++;
  }

  // If no neutrons were produced then don't continue
  if (nu == 0)
    return;

  // Initialize the counter of delayed neutrons encountered for each delayed
  // group.
  double nu_d[MAX_DELAYED_GROUPS] = {0.};

  // Clear out particle's nu fission bank
  p.nu_bank().clear();

  p.fission() = true;

  // Determine whether to place fission sites into the shared fission bank
  // or the secondary particle bank.
  bool use_fission_bank = (settings::run_mode == RunMode::EIGENVALUE);

  // Counter for the number of fission sites successfully stored to the shared
  // fission bank or the secondary particle bank
  int n_sites_stored;

  for (n_sites_stored = 0; n_sites_stored < nu; n_sites_stored++) {
    // Initialize fission site object with particle data
    SourceSite site;
    site.r = p.r();
    site.particle = ParticleType::neutron;
    site.time = p.time();
    site.wgt = 1. / weight;

    // Sample the cosine of the angle, assuming fission neutrons are emitted
    // isotropically
    double mu = 2. * prn(p.current_seed()) - 1.;

    // Sample the azimuthal angle uniformly in [0, 2.pi)
    double phi = 2. * PI * prn(p.current_seed());
    site.u.x = mu;
    site.u.y = std::sqrt(1. - mu * mu) * std::cos(phi);
    site.u.z = std::sqrt(1. - mu * mu) * std::sin(phi);

    // Sample secondary energy distribution for the fission reaction
    int dg;
    int gout;
    data::mg.macro_xs_[p.material()].sample_fission_energy(
      p.g(), dg, gout, p.current_seed(), p.mg_xs_cache().t, p.mg_xs_cache().a);

    // Store the energy and delayed groups on the fission bank
    site.E = gout;

    // We add 1 to the delayed_group bc in MG, -1 is prompt, but in the rest
    // of the code, 0 is prompt.
    site.delayed_group = dg + 1;

    // If delayed product production, sample time of emission
    if (dg != -1) {
      auto& macro_xs = data::mg.macro_xs_[p.material()];
      double decay_rate =
        macro_xs.get_xs(MgxsType::DECAY_RATE, 0, nullptr, nullptr, &dg, 0, 0);
      site.time -= std::log(prn(p.current_seed())) / decay_rate;

      // Reject site if it exceeds time cutoff
      double t_cutoff = settings::time_cutoff[static_cast<int>(site.particle)];
      if (site.time > t_cutoff) {
        continue;
      }
    }

    // Store fission site in bank
    if (use_fission_bank) {
      int64_t idx = simulation::fission_bank.thread_safe_append(site);
      if (idx == -1) {
        warning(
          "The shared fission bank is full. Additional fission sites created "
          "in this generation will not be banked. Results may be "
          "non-deterministic.");

        // Break out of loop as no more sites can be added to fission bank
        break;
      }
    } else {
      p.secondary_bank().push_back(site);
    }

    // Set parent and progeny ID
    site.parent_id = p.id();
    site.progeny_id = p.n_progeny()++;

    // Set the delayed group on the particle as well
    p.delayed_group() = dg + 1;

    // Increment the number of neutrons born delayed
    if (p.delayed_group() > 0) {
      nu_d[dg]++;
    }

    // Write fission particles to nuBank
    p.nu_bank().emplace_back();
    NuBank* nu_bank_entry = &p.nu_bank().back();
    nu_bank_entry->wgt = site.wgt;
    nu_bank_entry->E = site.E;
    nu_bank_entry->delayed_group = site.delayed_group;
  }

  // If shared fission bank was full, and no fissions could be added,
  // set the particle fission flag to false.
  if (n_sites_stored == 0) {
    p.fission() = false;
    return;
  }

  // Set nu to the number of fission sites successfully stored. If the fission
  // bank was not found to be full then these values are already equivalent.
  nu = n_sites_stored;

  // Store the total weight banked for analog fission tallies
  p.n_bank() = nu;
  p.wgt_bank() = nu / weight;
  for (size_t d = 0; d < MAX_DELAYED_GROUPS; d++) {
    p.n_delayed_bank(d) = nu_d[d];
  }
}

void absorption(Particle& p)
{
  if (settings::survival_biasing) {
    // Determine weight absorbed in survival biasing
    double wgt_absorb = p.wgt() * p.macro_xs().absorption / p.macro_xs().total;

    // Adjust weight of particle by the probability of absorption
    p.wgt() -= wgt_absorb;

    // Score implicit absorpion estimate of keff
    p.keff_tally_absorption() +=
      wgt_absorb * p.macro_xs().nu_fission / p.macro_xs().absorption;
  } else {
    if (p.macro_xs().absorption > prn(p.current_seed()) * p.macro_xs().total) {
      p.keff_tally_absorption() +=
        p.wgt() * p.macro_xs().nu_fission / p.macro_xs().absorption;
      p.wgt() = 0.0;
      p.event() = TallyEvent::ABSORB;
    }
  }
}

} // namespace openmc
