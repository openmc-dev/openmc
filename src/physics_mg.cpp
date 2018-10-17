#include "openmc/physics_mg.h"

#include <stdexcept>
#include <sstream>

#include "xtensor/xarray.hpp"

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
collision_mg(Particle* p, const double* energy_bin_avg,
             const MaterialMacroXS* material_xs)
{
  // Add to the collision counter for the particle
  p->n_collision++;

  // Sample the reaction type
  sample_reaction(p, energy_bin_avg, material_xs);

  // Display information about collision
  if ((settings::verbosity >= 10) || (simulation::trace)) {
    std::stringstream msg;
    msg << "    Energy Group = " << p->g;
    write_message(msg, 1);
  }
}

void
sample_reaction(Particle* p, const double* energy_bin_avg,
                const MaterialMacroXS* material_xs)
{
  // Create fission bank sites. Note that while a fission reaction is sampled,
  // it never actually "happens", i.e. the weight of the particle does not
  // change when sampling fission sites. The following block handles all
  // absorption (including fission)

  if (materials[p->material - 1]->fissionable) {
    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      Bank* result_bank;
      int64_t result_bank_size;
      // Get pointer to fission bank from Fortran side
      openmc_fission_bank(&result_bank, &result_bank_size);
      create_fission_sites(p, result_bank, &n_bank, result_bank_size,
                           material_xs);
    } else if ((settings::run_mode == RUN_MODE_FIXEDSOURCE) &&
               (settings::create_fission_neutrons)) {
      create_fission_sites(p, p->secondary_bank, &(p->n_secondary),
                           MAX_SECONDARY, material_xs);
    }
  }

  // If survival biasing is being used, the following subroutine adjusts the
  // weight of the particle. Otherwise, it checks to see if absorption occurs.
  if (material_xs->absorption > 0.) {
    absorption(p, material_xs);
  } else {
    p->absorb_wgt = 0.;
  }
  if (!p->alive) return;

  // Sample a scattering event to determine the energy of the exiting neutron
  scatter(p, energy_bin_avg);

  // Play Russian roulette if survival biasing is turned on
  if (settings::survival_biasing) {
    russian_roulette(p);
    if (!p->alive) return;
  }
}

void
scatter(Particle* p, const double* energy_bin_avg)
{
  // Adjust indices for Fortran to C++ indexing
  // TODO: Remove when no longer needed
  int gin = p->last_g - 1;
  int gout = p->g - 1;
  int i_mat = p->material - 1;
  macro_xs[i_mat].sample_scatter(gin, gout, p->mu, p->wgt);

  // Adjust return value for fortran indexing
  // TODO: Remove when no longer needed
  p->g = gout + 1;

  // Rotate the angle
  rotate_angle_c(p->coord[0].uvw, p->mu, nullptr);

  // Update energy value for downstream compatability (in tallying)
  p->E = energy_bin_avg[gout];

  // Set event component
  p->event = EVENT_SCATTER;
}

void
create_fission_sites(Particle* p, Bank* bank_array, int64_t* size_bank,
     int64_t bank_array_size, const MaterialMacroXS* material_xs)
{
  // TODO: Heat generation from fission

  // If uniform fission source weighting is turned on, we increase or decrease
  // the expected number of fission sites produced
  double weight = settings::ufs_on ? ufs_get_weight(p) : 1.0;

  // Determine the expected number of neutrons produced
  double nu_t = p->wgt / simulation::keff * weight * material_xs->nu_fission /
       material_xs->total;

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

  p->fission = true;
  for (size_t i = static_cast<size_t>(*size_bank);
       i < static_cast<size_t>(std::min(*size_bank + nu, bank_array_size)); i++) {
    // Bank source neutrons by copying the particle data
    bank_array[i].xyz[0] = p->coord[0].xyz[0];
    bank_array[i].xyz[1] = p->coord[0].xyz[1];
    bank_array[i].xyz[2] = p->coord[0].xyz[2];

    // Set that the bank particle is a neutron
    bank_array[i].particle = static_cast<int>(ParticleType::neutron);

    // Set the weight of the fission bank site
    bank_array[i].wgt = 1. / weight;

    // Sample the cosine of the angle, assuming fission neutrons are emitted
    // isotropically
    double mu = 2. * prn() - 1.;

    // Sample the azimuthal angle uniformly in [0, 2.pi)
    double phi = 2. * PI * prn();
    bank_array[i].uvw[0] = mu;
    bank_array[i].uvw[1] = std::sqrt(1. - mu * mu) * std::cos(phi);
    bank_array[i].uvw[2] = std::sqrt(1. - mu * mu) * std::sin(phi);

    // Sample secondary energy distribution for the fission reaction and set
    // the energy in the fission bank
    int dg;
    int gout;
    macro_xs[p->material - 1].sample_fission_energy(p->g - 1, dg, gout);
    bank_array[i].E = static_cast<double>(gout + 1);
    bank_array[i].delayed_group = dg + 1;

    // Set the delayed group on the particle as well
    p->delayed_group = dg + 1;

    // Increment the number of neutrons born delayed
    if (p->delayed_group > 0) {
      nu_d[dg]++;
    }
  }

  // Increment number of bank sites
  *size_bank = std::min(*size_bank + nu, bank_array_size);

  // Store the total weight banked for analog fission tallies
  p->n_bank = nu;
  p->wgt_bank = nu / weight;
  for (size_t d = 0; d < MAX_DELAYED_GROUPS; d++) {
    p->n_delayed_bank[d] = nu_d[d];
  }
}

void
absorption(Particle* p, const MaterialMacroXS* material_xs)
{
  if (settings::survival_biasing) {
    // Determine weight absorbed in survival biasing
    p->absorb_wgt = p->wgt * material_xs->absorption / material_xs->total;

    // Adjust weight of particle by the probability of absorption
    p->wgt -= p->absorb_wgt;
    p->last_wgt = p->wgt;

    // Score implicit absorpion estimate of keff
#pragma omp atomic
    global_tally_absorption += p->absorb_wgt * material_xs->nu_fission /
         material_xs->absorption;
  } else {
    if (material_xs->absorption > prn() * material_xs->total) {
#pragma omp atomic
      global_tally_absorption += p->wgt * material_xs->nu_fission /
           material_xs->absorption;
      p->alive = false;
      p->event = EVENT_ABSORB;
    }

  }
}

} //namespace openmc
