#include "openmc/physics_mg.h"

#include <stdexcept>

#include <fmt/core.h>
#include "xtensor/xarray.hpp"

#include "openmc/bank.h"
#include "openmc/constants.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/math_functions.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/physics_common.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/tallies/tally.h"

namespace openmc {

void
collision_mg(Particle* p)
{
  // Add to the collision counter for the particle
  p->n_collision_++;

  // Sample the reaction type
  sample_reaction(p);

  // Display information about collision
  if ((settings::verbosity >= 10) || p->trace_) {
    write_message(fmt::format("    Energy Group = {}", p->g_), 1);
  }
}

void
sample_reaction(Particle* p)
{
  // Create fission bank sites. Note that while a fission reaction is sampled,
  // it never actually "happens", i.e. the weight of the particle does not
  // change when sampling fission sites. The following block handles all
  // absorption (including fission)

  if (model::materials[p->material_]->fissionable_) {
    if (settings::run_mode == RunMode::EIGENVALUE  ||
       (settings::run_mode == RunMode::FIXED_SOURCE &&
        settings::create_fission_neutrons)) {
      create_fission_sites(p);
    }
  }

  // If survival biasing is being used, the following subroutine adjusts the
  // weight of the particle. Otherwise, it checks to see if absorption occurs.
  if (p->macro_xs_.absorption > 0.) {
    absorption(p);
  } else {
    p->wgt_absorb_ = 0.;
  }
  if (!p->alive_) return;

  // Sample a scattering event to determine the energy of the exiting neutron
  scatter(p);

  // Play Russian roulette if survival biasing is turned on
  if (settings::survival_biasing) {
    russian_roulette(p);
    if (!p->alive_) return;
  }
}

void
scatter(Particle* p)
{
  data::mg.macro_xs_[p->material_].sample_scatter(p->g_last_, p->g_, p->mu_,
                                                  p->wgt_, p->current_seed());

  // Rotate the angle
  p->u() = rotate_angle(p->u(), p->mu_, nullptr, p->current_seed());

  // Update energy value for downstream compatability (in tallying)
  p->E_ = data::mg.energy_bin_avg_[p->g_];

  // Set event component
  p->event_ = TallyEvent::SCATTER;
}

void
create_fission_sites(Particle* p)
{
  // If uniform fission source weighting is turned on, we increase or decrease
  // the expected number of fission sites produced
  double weight = settings::ufs_on ? ufs_get_weight(p) : 1.0;

  // Determine the expected number of neutrons produced
  double nu_t = p->wgt_ / simulation::keff * weight *
       p->macro_xs_.nu_fission / p->macro_xs_.total;

  // Sample the number of neutrons produced
  int nu = static_cast<int>(nu_t);
  if (prn(p->current_seed()) <= (nu_t - int(nu_t))) {
    nu++;
  }

  // Begin banking the source neutrons
  // First, if our bank is full then don't continue
  if (nu == 0) return;

  // Initialize the counter of delayed neutrons encountered for each delayed
  // group.
  double nu_d[MAX_DELAYED_GROUPS] = {0.};

  // Clear out particle's nu fission bank
  p->nu_bank_.clear();

  p->fission_ = true;
  int skipped = 0;

  // Determine whether to place fission sites into the shared fission bank
  // or the secondary particle bank.
  bool use_fission_bank = (settings::run_mode == RunMode::EIGENVALUE);

  for (int i = 0; i < nu; ++i) {
    // Initialize fission site object with particle data
    Particle::Bank site;
    site.r = p->r();
    site.particle = Particle::Type::neutron;
    site.wgt = 1. / weight;
    site.parent_id = p->id_;
    site.progeny_id = p->n_progeny_++;

    // Sample the cosine of the angle, assuming fission neutrons are emitted
    // isotropically
    double mu = 2.*prn(p->current_seed()) - 1.;

    // Sample the azimuthal angle uniformly in [0, 2.pi)
    double phi = 2. * PI * prn(p->current_seed() );
    site.u.x = mu;
    site.u.y = std::sqrt(1. - mu * mu) * std::cos(phi);
    site.u.z = std::sqrt(1. - mu * mu) * std::sin(phi);

    // Sample secondary energy distribution for the fission reaction
    int dg;
    int gout;
    data::mg.macro_xs_[p->material_].sample_fission_energy(p->g_, dg, gout,
      p->current_seed());

    // Store the energy and delayed groups on the fission bank
    site.E = gout;

    // We add 1 to the delayed_group bc in MG, -1 is prompt, but in the rest
    // of the code, 0 is prompt.
    site.delayed_group = dg + 1;

    // Store fission site in bank
    if (use_fission_bank) {
      int64_t idx = simulation::fission_bank.thread_safe_append(site);
      if (idx == -1) {
        warning("The shared fission bank is full. Additional fission sites created "
            "in this generation will not be banked.");
        skipped++;
        break;
      }
    } else {
      p->secondary_bank_.push_back(site);
    }

    // Set the delayed group on the particle as well
    p->delayed_group_ = dg + 1;

    // Increment the number of neutrons born delayed
    if (p->delayed_group_ > 0) {
      nu_d[dg]++;
    }

    // Write fission particles to nuBank
    if (use_fission_bank) {
      p->nu_bank_.emplace_back();
      Particle::NuBank* nu_bank_entry = &p->nu_bank_.back();
      nu_bank_entry->wgt              = site.wgt;
      nu_bank_entry->E                = site.E;
      nu_bank_entry->delayed_group    = site.delayed_group;
    }
  }

  // If shared fission bank was full, and no fissions could be added,
  // set the particle fission flag to false.
  if (nu == skipped) {
    p->fission_ = false;
    return;
  }

  // If shared fission bank was full, but some fissions could be added,
  // reduce nu accordingly
  nu -= skipped;

  // Store the total weight banked for analog fission tallies
  p->n_bank_ = nu;
  p->wgt_bank_ = nu / weight;
  for (size_t d = 0; d < MAX_DELAYED_GROUPS; d++) {
    p->n_delayed_bank_[d] = nu_d[d];
  }
}

void
absorption(Particle* p)
{
  if (settings::survival_biasing) {
    // Determine weight absorbed in survival biasing
    p->wgt_absorb_ = p->wgt_ * p->macro_xs_.absorption / p->macro_xs_.total;

    // Adjust weight of particle by the probability of absorption
    p->wgt_ -= p->wgt_absorb_;
    p->wgt_last_ = p->wgt_;

    // Score implicit absorpion estimate of keff
    p->keff_tally_absorption_ += p->wgt_absorb_ * p->macro_xs_.nu_fission /
        p->macro_xs_.absorption;
  } else {
    if (p->macro_xs_.absorption > prn(p->current_seed()) * p->macro_xs_.total) {
      p->keff_tally_absorption_ += p->wgt_ * p->macro_xs_.nu_fission /
           p->macro_xs_.absorption;
      p->alive_ = false;
      p->event_ = TallyEvent::ABSORB;
    }

  }
}

} //namespace openmc
