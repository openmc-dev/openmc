#include "openmc/physics_mg.h"

#include <stdexcept>
#include <sstream>

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
  if ((settings::verbosity >= 10) || (simulation::trace)) {
    std::stringstream msg;
    msg << "    Energy Group = " << p->g_;
    write_message(msg, 1);
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
    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      create_fission_sites(
        p, simulation::fission_bank.data(), &simulation::n_bank,
        simulation::fission_bank.size());
    } else if ((settings::run_mode == RUN_MODE_FIXEDSOURCE) &&
               (settings::create_fission_neutrons)) {
      create_fission_sites(p, p->secondary_bank_, &(p->n_secondary_),
                           MAX_SECONDARY);
    }
  }

  // If survival biasing is being used, the following subroutine adjusts the
  // weight of the particle. Otherwise, it checks to see if absorption occurs.
  if (simulation::material_xs.absorption > 0.) {
    absorption(p);
  } else {
    p->absorb_wgt_ = 0.;
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
  // Adjust indices for Fortran to C++ indexing
  // TODO: Remove when no longer needed
  int gin = p->last_g_ - 1;
  int gout = p->g_ - 1;
  int i_mat = p->material_;
  data::macro_xs[i_mat].sample_scatter(gin, gout, p->mu_, p->wgt_);

  // Adjust return value for fortran indexing
  // TODO: Remove when no longer needed
  p->g_ = gout + 1;

  // Rotate the angle
  p->u() = rotate_angle(p->u(), p->mu_, nullptr);

  // Update energy value for downstream compatability (in tallying)
  p->E_ = data::energy_bin_avg[gout];

  // Set event component
  p->event_ = EVENT_SCATTER;
}

void
create_fission_sites(Particle* p, Particle::Bank* bank_array, int64_t* size_bank,
     int64_t bank_array_size)
{
  // TODO: Heat generation from fission

  // If uniform fission source weighting is turned on, we increase or decrease
  // the expected number of fission sites produced
  double weight = settings::ufs_on ? ufs_get_weight(p) : 1.0;

  // Determine the expected number of neutrons produced
  double nu_t = p->wgt_ / simulation::keff * weight *
       simulation::material_xs.nu_fission / simulation::material_xs.total;

  // Sample the number of neutrons produced
  int nu = static_cast<int>(nu_t);
  if (prn() <= (nu_t - int(nu_t))) {
    nu++;
  }

  // Check for the bank size getting hit. For fixed source calculations, this
  // is a fatal error; for eigenvalue calculations, it just means that k-eff
  // was too high for a single batch.
  if (*size_bank + nu > bank_array_size) {
    if (settings::run_mode == RUN_MODE_FIXEDSOURCE) {
      throw std::runtime_error{"Secondary particle bank size limit reached."
           " If you are running a subcritical multiplication problem,"
           " k-effective may be too close to one."};
    } else {
      if (mpi::master) {
        std::stringstream msg;
        msg << "Maximum number of sites in fission bank reached. This can"
             " result in irreproducible results using different numbers of"
             " processes/threads.";
        warning(msg);
      }
    }
  }

  // Begin banking the source neutrons
  // First, if our bank is full then don't continue
  if ((nu == 0) || (*size_bank == bank_array_size)) return;

  // Initialize the counter of delayed neutrons encountered for each delayed
  // group.
  double nu_d[MAX_DELAYED_GROUPS] = {0.};

  p->fission_ = true;
  for (size_t i = static_cast<size_t>(*size_bank);
       i < static_cast<size_t>(std::min(*size_bank + nu, bank_array_size)); i++) {
    // Bank source neutrons by copying the particle data
    bank_array[i].r = p->r();

    // Set that the bank particle is a neutron
    bank_array[i].particle = Particle::Type::neutron;

    // Set the weight of the fission bank site
    bank_array[i].wgt = 1. / weight;

    // Sample the cosine of the angle, assuming fission neutrons are emitted
    // isotropically
    double mu = 2. * prn() - 1.;

    // Sample the azimuthal angle uniformly in [0, 2.pi)
    double phi = 2. * PI * prn();
    bank_array[i].u.x = mu;
    bank_array[i].u.y = std::sqrt(1. - mu * mu) * std::cos(phi);
    bank_array[i].u.z = std::sqrt(1. - mu * mu) * std::sin(phi);

    // Sample secondary energy distribution for the fission reaction and set
    // the energy in the fission bank
    int dg;
    int gout;
    data::macro_xs[p->material_].sample_fission_energy(p->g_ - 1, dg, gout);
    bank_array[i].E = gout + 1;
    bank_array[i].delayed_group = dg + 1;

    // Set the delayed group on the particle as well
    p->delayed_group_ = dg + 1;

    // Increment the number of neutrons born delayed
    if (p->delayed_group_ > 0) {
      nu_d[dg]++;
    }
  }

  // Increment number of bank sites
  *size_bank = std::min(*size_bank + nu, bank_array_size);

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
    p->absorb_wgt_ = p->wgt_ *
         simulation::material_xs.absorption / simulation::material_xs.total;

    // Adjust weight of particle by the probability of absorption
    p->wgt_ -= p->absorb_wgt_;
    p->last_wgt_ = p->wgt_;

    // Score implicit absorpion estimate of keff
#pragma omp atomic
    global_tally_absorption += p->absorb_wgt_ *
         simulation::material_xs.nu_fission /
         simulation::material_xs.absorption;
  } else {
    if (simulation::material_xs.absorption >
        prn() * simulation::material_xs.total) {
#pragma omp atomic
      global_tally_absorption += p->wgt_ * simulation::material_xs.nu_fission /
           simulation::material_xs.absorption;
      p->alive_ = false;
      p->event_ = EVENT_ABSORB;
    }

  }
}

} //namespace openmc
