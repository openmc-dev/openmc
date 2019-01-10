#include "openmc/physics.h"

#include "openmc/bank.h"
#include "openmc/constants.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/math_functions.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/physics_common.h"
#include "openmc/random_lcg.h"
#include "openmc/reaction.h"
#include "openmc/secondary_uncorrelated.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/thermal.h"
#include "openmc/tallies/tally.h"

#include <algorithm> // for max, min, max_element
#include <cmath> // for sqrt, exp, log, abs, copysign
#include <sstream>

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

void collision(Particle* p)
{
  // Add to collision counter for particle
  ++(p->n_collision);

  // Sample reaction for the material the particle is in
  switch (static_cast<ParticleType>(p->type)) {
  case ParticleType::neutron:
    sample_neutron_reaction(p);
    break;
  case ParticleType::photon:
    sample_photon_reaction(p);
    break;
  case ParticleType::electron:
    sample_electron_reaction(p);
    break;
  case ParticleType::positron:
    sample_positron_reaction(p);
    break;
  }

  // Kill particle if energy falls below cutoff
  if (p->E < settings::energy_cutoff[p->type - 1]) {
    p->alive = false;
    p->wgt = 0.0;
    p->last_wgt = 0.0;
  }

  // Display information about collision
  if (settings::verbosity >= 10 || simulation::trace) {
    std::stringstream msg;
    if (static_cast<ParticleType>(p->type) == ParticleType::neutron) {
      msg << "    " << reaction_name(p->event_MT) << " with " <<
        data::nuclides[p->event_nuclide-1]->name_ << ". Energy = " << p->E << " eV.";
    } else {
      msg << "    " << reaction_name(p->event_MT) << " with " <<
        data::elements[p->event_nuclide-1].name_ << ". Energy = " << p->E << " eV.";
    }
    write_message(msg, 1);
  }
}

void sample_neutron_reaction(Particle* p)
{
  int i_nuclide;
  int i_nuc_mat;
  sample_nuclide(p, SCORE_TOTAL, &i_nuclide, &i_nuc_mat);

  // Save which nuclide particle had collision with
  // TODO: off-by-one
  p->event_nuclide = i_nuclide + 1;

  // Create fission bank sites. Note that while a fission reaction is sampled,
  // it never actually "happens", i.e. the weight of the particle does not
  // change when sampling fission sites. The following block handles all
  // absorption (including fission)

  const auto& nuc {data::nuclides[i_nuclide]};

  if (nuc->fissionable_) {
    Reaction* rx = sample_fission(i_nuclide, p->E);
    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      create_fission_sites(p, i_nuclide, rx, simulation::fission_bank.data(),
        &simulation::n_bank, simulation::fission_bank.size());
    } else if (settings::run_mode == RUN_MODE_FIXEDSOURCE &&
      settings::create_fission_neutrons) {
      create_fission_sites(p, i_nuclide, rx, p->secondary_bank,
        &p->n_secondary, MAX_SECONDARY);
    }
  }

  // Create secondary photons
  if (settings::photon_transport) {
    prn_set_stream(STREAM_PHOTON);
    sample_secondary_photons(p, i_nuclide);
    prn_set_stream(STREAM_TRACKING);
  }

  // If survival biasing is being used, the following subroutine adjusts the
  // weight of the particle. Otherwise, it checks to see if absorption occurs

  if (simulation::micro_xs[i_nuclide].absorption > 0.0) {
    absorption(p, i_nuclide);
  } else {
    p->absorb_wgt = 0.0;
  }
  if (!p->alive) return;

  // Sample a scattering reaction and determine the secondary energy of the
  // exiting neutron
  scatter(p, i_nuclide, i_nuc_mat);

  // Advance URR seed stream 'N' times after energy changes
  if (p->E != p->last_E) {
    prn_set_stream(STREAM_URR_PTABLE);
    advance_prn_seed(data::nuclides.size());
    prn_set_stream(STREAM_TRACKING);
  }

  // Play russian roulette if survival biasing is turned on
  if (settings::survival_biasing) {
    russian_roulette(p);
    if (!p->alive) return;
  }
}

void
create_fission_sites(Particle* p, int i_nuclide, const Reaction* rx, Bank* bank_array,
  int64_t* size_bank, int64_t bank_capacity)
{
  // TODO: Heat generation from fission

  // If uniform fission source weighting is turned on, we increase or decrease
  // the expected number of fission sites produced
  double weight = settings::ufs_on ? ufs_get_weight(p) : 1.0;

  // Determine the expected number of neutrons produced
  double nu_t = p->wgt / simulation::keff * weight * simulation::micro_xs[
    i_nuclide].nu_fission / simulation::micro_xs[i_nuclide].total;

  // Sample the number of neutrons produced
  int nu = static_cast<int>(nu_t);
  if (prn() <= (nu_t - nu)) ++nu;

  // Check for the bank size getting hit. For fixed source calculations, this
  // is a fatal error; for eigenvalue calculations, it just means that k-eff
  // was too high for a single batch.
  if (*size_bank + nu > bank_capacity) {
    if (settings::run_mode == RUN_MODE_FIXEDSOURCE) {
      throw std::runtime_error{"Secondary particle bank size limit reached."
           " If you are running a subcritical multiplication problem,"
           " k-effective may be too close to one."};
    } else {
      if (mpi::master) {
        warning("Maximum number of sites in fission bank reached. This can"
          " result in irreproducible results using different numbers of"
          " processes/threads.");
      }
    }
  }

  // Begin banking the source neutrons
  // First, if our bank is full then don't continue
  if (nu == 0 || *size_bank == bank_capacity) return;

  // Initialize the counter of delayed neutrons encountered for each delayed
  // group.
  double nu_d[MAX_DELAYED_GROUPS] = {0.};

  p->fission = true;
  for (size_t i = *size_bank; i < std::min(*size_bank + nu, bank_capacity); ++i) {
    // Bank source neutrons by copying the particle data
    bank_array[i].xyz[0] = p->coord[0].xyz[0];
    bank_array[i].xyz[1] = p->coord[0].xyz[1];
    bank_array[i].xyz[2] = p->coord[0].xyz[2];

    // Set that the bank particle is a neutron
    bank_array[i].particle = static_cast<int>(ParticleType::neutron);

    // Set the weight of the fission bank site
    bank_array[i].wgt = 1. / weight;

    // Sample delayed group and angle/energy for fission reaction
    sample_fission_neutron(i_nuclide, rx, p->E, &bank_array[i]);

    // Set the delayed group on the particle as well
    p->delayed_group = bank_array[i].delayed_group;

    // Increment the number of neutrons born delayed
    if (p->delayed_group > 0) {
      nu_d[p->delayed_group-1]++;
    }
  }

  // Increment number of bank sites
  *size_bank = std::min(*size_bank + nu, bank_capacity);

  // Store the total weight banked for analog fission tallies
  p->n_bank = nu;
  p->wgt_bank = nu / weight;
  for (size_t d = 0; d < MAX_DELAYED_GROUPS; d++) {
    p->n_delayed_bank[d] = nu_d[d];
  }
}

// TODO: Finish converting photon physics functions

// void sample_photon_reaction(Particle* p)
// {
//   // Kill photon if below energy cutoff -- an extra check is made here because
//   // photons with energy below the cutoff may have been produced by neutrons
//   // reactions or atomic relaxation
//   if (p->E < energy_cutoff(PHOTON)) {
//     p->E = 0.0
//     p->alive = false
//     return
//   }

//   // Sample element within material
//   i_element = sample_element(p)
//   p->event_nuclide = i_element

//   // Calculate photon energy over electron rest mass equivalent
//   alpha = p->E/MASS_ELECTRON_EV

//   // For tallying purposes, this routine might be called directly. In that
//   // case, we need to sample a reaction via the cutoff variable
//   prob = 0.0
//   cutoff = prn() * micro_photon_xs(i_element) % total

//   associate (elm => elements(i_element))
//     // Coherent (Rayleigh) scattering
//     prob = prob + micro_photon_xs(i_element) % coherent
//     if (prob > cutoff) {
//       rayleigh_scatter(elm, alpha, mu)
//       p->coord(1) % uvw = rotate_angle(p->coord(1) % uvw, mu)
//       p->event_MT = COHERENT
//       return
//     }

//     // Incoherent (Compton) scattering
//     prob = prob + micro_photon_xs(i_element) % incoherent
//     if (prob > cutoff) {
//       compton_scatter(elm, alpha, alpha_out, mu, i_shell, true)

//       // Determine binding energy of shell. The binding energy is 0.0 if
//       // doppler broadening is not used.
//       if (i_shell == 0) {
//         e_b = 0.0
//       } else {
//         e_b = elm % binding_energy(i_shell)
//       }

//       // Create Compton electron
//       E_electron = (alpha - alpha_out)*MASS_ELECTRON_EV - e_b
//       mu_electron = (alpha - alpha_out*mu) &
//             / std::sqrt(alpha**2 + alpha_out**2 - 2.0*alpha*alpha_out*mu)
//       phi = 2.0*PI*prn()
//       uvw = rotate_angle(p->coord(1) % uvw, mu_electron, phi)
//       particle_create_secondary(p, uvw, E_electron, ELECTRON, true)

//       // TODO: Compton subshell data does not match atomic relaxation data
//       // Allow electrons to fill orbital and produce auger electrons
//       // and fluorescent photons
//       if (i_shell > 0) {
//         atomic_relaxation(p, elm, i_shell)
//       }

//       phi = phi + PI
//       p->E = alpha_out*MASS_ELECTRON_EV
//       p->coord(1) % uvw = rotate_angle(p->coord(1) % uvw, mu, phi)
//       p->event_MT = INCOHERENT
//       return
//     }

//     // Photoelectric effect
//     prob_after = prob + micro_photon_xs(i_element) % photoelectric
//     if (prob_after > cutoff) {
//       do i_shell = 1, size(elm % shells)
//         // Get grid index and interpolation factor
//         i_grid = micro_photon_xs(i_element) % index_grid
//         f      = micro_photon_xs(i_element) % interp_factor

//         // Check threshold of reaction
//         i_start = elm % shells(i_shell) % threshold
//         if (i_grid <= i_start) cycle

//         // Evaluation subshell photoionization cross section
//         xs = std::exp(elm % shells(i_shell) % cross_section(i_grid - i_start) + &
//               f*(elm % shells(i_shell) % cross_section(i_grid + 1 - i_start) - &
//               elm % shells(i_shell) % cross_section(i_grid - i_start)))

//         prob = prob + xs
//         if (prob > cutoff) {
//           E_electron = p->E - elm % shells(i_shell) % binding_energy

//           // Sample mu using non-relativistic Sauter distribution.
//           // See Eqns 3.19 and 3.20 in "Implementing a photon physics
//           // model in Serpent 2" by Toni Kaltiaisenaho
//           SAMPLE_MU: do
//             r = prn()
//             if (FOUR * (1.0 - r) * r >= prn()) {
//               rel_vel = std::sqrt(E_electron * (E_electron + 2.0 * MASS_ELECTRON_EV))&
//                     / (E_electron + MASS_ELECTRON_EV)
//               mu = (2.0 * r + rel_vel - 1.0) / &
//                     (2.0 * rel_vel * r - rel_vel + 1.0)
//               exit SAMPLE_MU
//             }
//           end do SAMPLE_MU

//           phi = 2.0*PI*prn()
//           uvw(1) = mu
//           uvw(2) = std::sqrt(1.0 - mu*mu)*std::cos(phi)
//           uvw(3) = std::sqrt(1.0 - mu*mu)*std::sin(phi)

//           // Create secondary electron
//           particle_create_secondary(p, uvw, E_electron, ELECTRON, &
//                 run_CE=true)

//           // Allow electrons to fill orbital and produce auger electrons
//           // and fluorescent photons
//           atomic_relaxation(p, elm, i_shell)
//           p->event_MT = 533 + elm % shells(i_shell) % index_subshell
//           p->alive = false
//           p->E = 0.0

//           return
//         }
//       end do
//     }
//     prob = prob_after

//     // Pair production
//     prob = prob + micro_photon_xs(i_element) % pair_production
//     if (prob > cutoff) {
//       pair_production(elm, alpha, E_electron, E_positron, mu_electron, &
//             mu_positron)

//       // Create secondary electron
//       uvw = rotate_angle(p->coord(1) % uvw, mu_electron)
//       particle_create_secondary(p, uvw, E_electron, ELECTRON, true)

//       // Create secondary positron
//       uvw = rotate_angle(p->coord(1) % uvw, mu_positron)
//       particle_create_secondary(p, uvw, E_positron, POSITRON, true)

//       p->event_MT = PAIR_PROD
//       p->alive = false
//       p->E = 0.0
//     }

//   end associate
// }

void sample_electron_reaction(Particle* p)
{
  // TODO: create reaction types

  if (settings::electron_treatment == ELECTRON_TTB) {
    double E_lost;
    thick_target_bremsstrahlung(p, &E_lost);
  }

  p->E = 0.0;
  p->alive = false;
}

void sample_positron_reaction(Particle* p)
{
  // TODO: create reaction types

  if (settings::electron_treatment == ELECTRON_TTB) {
    double E_lost;
    thick_target_bremsstrahlung(p, &E_lost);
  }

  // Sample angle isotropically
  double mu = 2.0*prn() - 1.0;
  double phi = 2.0*PI*prn();
  std::array<double, 3> uvw;
  uvw[0] = mu;
  uvw[1] = std::sqrt(1.0 - mu*mu)*std::cos(phi);
  uvw[2] = std::sqrt(1.0 - mu*mu)*std::sin(phi);

  // Create annihilation photon pair traveling in opposite directions
  int photon = static_cast<int>(ParticleType::photon);
  p->create_secondary(uvw.data(), MASS_ELECTRON_EV, photon, true);

  uvw[0] = -uvw[0];
  uvw[1] = -uvw[1];
  uvw[2] = -uvw[2];
  p->create_secondary(uvw.data(), MASS_ELECTRON_EV, photon, true);

  p->E = 0.0;
  p->alive = false;
}

void sample_nuclide(const Particle* p, int mt, int* i_nuclide, int* i_nuc_mat)
{
  // Sample cumulative distribution function
  double cutoff;
  switch (mt) {
  case SCORE_TOTAL:
    cutoff = prn() * simulation::material_xs.total;
    break;
  case SCORE_SCATTER:
    cutoff = prn() * (simulation::material_xs.total -
      simulation::material_xs.absorption);
    break;
  case SCORE_FISSION:
    cutoff = prn() * simulation::material_xs.fission;
    break;
  }

  // Get pointers to nuclide/density arrays
  int* nuclides;
  double* densities;
  int n;
  openmc_material_get_densities(p->material, &nuclides, &densities, &n);

  *i_nuc_mat = 0;
  double prob = 0.0;
  while (prob < cutoff) {
    // Check to make sure that a nuclide was sampled
    if (*i_nuc_mat > n) {
      p->write_restart();
      fatal_error("Did not sample any nuclide during collision.");
    }

    // Find atom density
    // TODO: off-by-one
    *i_nuclide = nuclides[*i_nuc_mat] - 1;
    double atom_density = densities[*i_nuc_mat];

    // Determine microscopic cross section
    double sigma;
    switch (mt) {
    case SCORE_TOTAL:
      sigma = atom_density * simulation::micro_xs[*i_nuclide].total;
      break;
    case SCORE_SCATTER:
      sigma = atom_density * (simulation::micro_xs[*i_nuclide].total -
        simulation::micro_xs[*i_nuclide].absorption);
        break;
    case SCORE_FISSION:
      sigma = atom_density * simulation::micro_xs[*i_nuclide].fission;
      break;
    }

    // Increment probability to compare to cutoff
    prob += sigma;

    ++(*i_nuc_mat);
  }
}

// TODO: Finish converting photon physics functions

// void sample_element(Particle* p)
// {
//   associate (mat => materials(p->material))
//     // Sample cumulative distribution function
//     cutoff = prn() * simulation::material_xs.total

//     i = 0
//     prob = 0.0
//     do while (prob < cutoff)
//       i = i + 1

//       // Check to make sure that a nuclide was sampled
//       if (i > mat % n_nuclides) {
//         particle_write_restart(p)
//         fatal_error("Did not sample any element during collision.")
//       }

//       // Find atom density
//       i_element    = mat % element(i)
//       atom_density = mat % atom_density(i)

//       // Determine microscopic cross section
//       sigma = atom_density * micro_photon_xs(i_element) % total

//       // Increment probability to compare to cutoff
//       prob = prob + sigma
//     end do
//   end associate
// }

Reaction* sample_fission(int i_nuclide, double E)
{
  // Get pointer to nuclide
  const auto& nuc {data::nuclides[i_nuclide]};

  // If we're in the URR, by default use the first fission reaction. We also
  // default to the first reaction if we know that there are no partial fission
  // reactions
  if (simulation::micro_xs[i_nuclide].use_ptable || !nuc->has_partial_fission_) {
    return nuc->fission_rx_[0];
  }

  // Check to see if we are in a windowed multipole range.  WMP only supports
  // the first fission reaction.
  if (nuc->multipole_) {
    if (E >= nuc->multipole_->E_min_ && E <= nuc->multipole_->E_max_) {
      return nuc->fission_rx_[0];
    }
  }

  // Get grid index and interpolatoin factor and sample fission cdf
  int i_temp = simulation::micro_xs[i_nuclide].index_temp;
  int i_grid = simulation::micro_xs[i_nuclide].index_grid;
  double f = simulation::micro_xs[i_nuclide].interp_factor;
  double cutoff = prn() * simulation::micro_xs[i_nuclide].fission;
  double prob = 0.0;

  // Loop through each partial fission reaction type
  for (auto& rx : nuc->fission_rx_) {
    // if energy is below threshold for this reaction, skip it
    int threshold = rx->xs_[i_temp].threshold;
    if (i_grid < threshold) continue;

    // add to cumulative probability
    prob += (1.0 - f) * rx->xs_[i_temp].value[i_grid - threshold]
            + f*rx->xs_[i_temp].value[i_grid - threshold + 1];

    // Create fission bank sites if fission occurs
    if (prob > cutoff) return rx;
  }
}

void sample_photon_product(int i_nuclide, double E, int* i_rx, int* i_product)
{
  // Get grid index and interpolation factor and sample photon production cdf
  int i_temp = simulation::micro_xs[i_nuclide].index_temp;
  int i_grid = simulation::micro_xs[i_nuclide].index_grid;
  double f = simulation::micro_xs[i_nuclide].interp_factor;
  double cutoff = prn() * simulation::micro_xs[i_nuclide].photon_prod;
  double prob = 0.0;

  // Loop through each reaction type
  const auto& nuc {data::nuclides[i_nuclide]};
  for (int i = 0; i < nuc->reactions_.size(); ++i) {
    const auto& rx = nuc->reactions_[i];
    int threshold = rx->xs_[i_temp].threshold;

    // if energy is below threshold for this reaction, skip it
    if (i_grid < threshold) continue;

    // Evaluate neutron cross section
    double xs = ((1.0 - f) * rx->xs_[i_temp].value[i_grid - threshold]
      + f*(rx->xs_[i_temp].value[i_grid - threshold + 1]));

    for (int j = 0; j < rx->products_.size(); ++j) {
      if (rx->products_[j].particle_ == ParticleType::photon) {
        // add to cumulative probability
        prob += (*rx->products_[j].yield_)(E) * xs;

        *i_rx = i;
        *i_product = j;
        if (prob > cutoff) return;
      }
    }
  }
}

void absorption(Particle* p, int i_nuclide)
{
  if (settings::survival_biasing) {
    // Determine weight absorbed in survival biasing
    p->absorb_wgt = p->wgt * simulation::micro_xs[i_nuclide].absorption /
          simulation::micro_xs[i_nuclide].total;

    // Adjust weight of particle by probability of absorption
    p->wgt -= p->absorb_wgt;
    p->last_wgt = p->wgt;

    // Score implicit absorption estimate of keff
    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      global_tally_absorption += p->absorb_wgt * simulation::micro_xs[
        i_nuclide].nu_fission / simulation::micro_xs[i_nuclide].absorption;
    }
  } else {
    // See if disappearance reaction happens
    if (simulation::micro_xs[i_nuclide].absorption >
        prn() * simulation::micro_xs[i_nuclide].total) {
      // Score absorption estimate of keff
      if (settings::run_mode == RUN_MODE_EIGENVALUE) {
        global_tally_absorption += p->wgt * simulation::micro_xs[
          i_nuclide].nu_fission / simulation::micro_xs[i_nuclide].absorption;
      }

      p->alive = false;
      p->event = EVENT_ABSORB;
      p->event_MT = N_DISAPPEAR;
    }
  }
}

void scatter(Particle* p, int i_nuclide, int i_nuc_mat)
{
  // copy incoming direction
  Direction u_old {p->coord[0].uvw};

  // Get pointer to nuclide and grid index/interpolation factor
  const auto& nuc {data::nuclides[i_nuclide]};
  const auto& micro {simulation::micro_xs[i_nuclide]};
  int i_temp =  micro.index_temp;
  int i_grid =  micro.index_grid - 1;
  double f = micro.interp_factor;

  // For tallying purposes, this routine might be called directly. In that
  // case, we need to sample a reaction via the cutoff variable
  double cutoff = prn() * (micro.total - micro.absorption);
  bool sampled = false;

  // Calculate elastic cross section if it wasn't precalculated
  if (micro.elastic == CACHE_INVALID) {
    nuc->calculate_elastic_xs();
  }

  double prob = micro.elastic - micro.thermal;
  if (prob > cutoff) {
    // =======================================================================
    // NON-S(A,B) ELASTIC SCATTERING

    // Determine temperature
    double kT = nuc->multipole_ ? p->sqrtkT*p->sqrtkT : nuc->kTs_[i_temp];

    // Perform collision physics for elastic scattering
    elastic_scatter(i_nuclide, nuc->reactions_[0].get(), kT,
      &p->E, p->coord[0].uvw, &p->mu, &p->wgt);

    p->event_MT = ELASTIC;
    sampled = true;
  }

  prob = micro.elastic;
  if (prob > cutoff && !sampled) {
    // =======================================================================
    // S(A,B) SCATTERING

    sab_scatter(i_nuclide, micro.index_sab, &p->E, p->coord[0].uvw, &p->mu);

    p->event_MT = ELASTIC;
    sampled = true;
  }

  if (!sampled) {
    // =======================================================================
    // INELASTIC SCATTERING

    int j = 0;
    int i;
    while (prob < cutoff) {
      i = nuc->index_inelastic_scatter_[j];
      ++j;

      // Check to make sure inelastic scattering reaction sampled
      if (i >= nuc->reactions_.size()) {
        p->write_restart();
        fatal_error("Did not sample any reaction for nuclide " + nuc->name_);
      }

      // if energy is below threshold for this reaction, skip it
      const auto& xs {nuc->reactions_[i]->xs_[i_temp]};
      int threshold = xs.threshold - 1;
      if (i_grid < threshold) continue;

      // add to cumulative probability
      prob += (1.0 - f)*xs.value[i_grid - threshold] +
        f*xs.value[i_grid - threshold + 1];
    }

    // Perform collision physics for inelastic scattering
    const auto& rx {nuc->reactions_[i]};
    inelastic_scatter(nuc.get(), rx.get(), p);
    p->event_MT = rx->mt_;
  }

  // Set event component
  p->event = EVENT_SCATTER;

  // Sample new outgoing angle for isotropic-in-lab scattering
  if (material_isotropic(p->material, i_nuc_mat)) {
    // Sample isotropic-in-lab outgoing direction
    double mu = 2.0*prn() - 1.0;
    double phi = 2.0*PI*prn();
    Direction u_new;
    u_new.x = mu;
    u_new.y = std::sqrt(1.0 - mu*mu)*std::cos(phi);
    u_new.z = std::sqrt(1.0 - mu*mu)*std::sin(phi);

    p->mu = u_old.dot(u_new);

    // Change direction of particle
    p->coord[0].uvw[0] = u_new.x;
    p->coord[0].uvw[1] = u_new.y;
    p->coord[0].uvw[2] = u_new.z;
  }
}

void elastic_scatter(int i_nuclide, const Reaction* rx, double kT, double* E,
  double* uvw, double* mu_lab, double* wgt)
{
  // get pointer to nuclide
  const auto& nuc {data::nuclides[i_nuclide]};

  double vel = std::sqrt(*E);
  double awr = nuc->awr_;

  // Neutron velocity in LAB
  Direction u {uvw};
  Direction v_n = vel*u;

  // Sample velocity of target nucleus
  Direction v_t {};
  if (!simulation::micro_xs[i_nuclide].use_ptable) {
    v_t = sample_target_velocity(nuc.get(), *E, u, v_n,
      simulation::micro_xs[i_nuclide].elastic, kT, wgt);
  }

  // Velocity of center-of-mass
  Direction v_cm = (v_n + awr*v_t)/(awr + 1.0);

  // Transform to CM frame
  v_n -= v_cm;

  // Find speed of neutron in CM
  vel = v_n.norm();

  // Sample scattering angle, checking if it is an ncorrelated angle-energy
  // distribution
  double mu_cm;
  auto& d = rx->products_[0].distribution_[0];
  auto d_ = dynamic_cast<UncorrelatedAngleEnergy*>(d.get());
  if (d_) {
    mu_cm = d_->angle().sample(*E);
  } else {
    mu_cm = 2.0*prn() - 1.0;
  }

  // Determine direction cosines in CM
  Direction u_cm = v_n/vel;

  // Rotate neutron velocity vector to new angle -- note that the speed of the
  // neutron in CM does not change in elastic scattering. However, the speed
  // will change when we convert back to LAB
  v_n = vel * rotate_angle(u_cm, mu_cm, nullptr);

  // Transform back to LAB frame
  v_n += v_cm;

  *E = v_n.dot(v_n);
  vel = std::sqrt(*E);

  // compute cosine of scattering angle in LAB frame by taking dot product of
  // neutron's pre- and post-collision angle
  *mu_lab = u.dot(v_n) / vel;

  // Set energy and direction of particle in LAB frame
  u = v_n / vel;
  uvw[0] = u.x;
  uvw[1] = u.y;
  uvw[2] = u.z;

  // Because of floating-point roundoff, it may be possible for mu_lab to be
  // outside of the range [-1,1). In these cases, we just set mu_lab to exactly
  // -1 or 1
  if (std::abs(*mu_lab) > 1.0) *mu_lab = std::copysign(1.0, *mu_lab);
}

void sab_scatter(int i_nuclide, int i_sab, double* E, double* uvw, double* mu)
{
  // Determine temperature index
  const auto& micro {simulation::micro_xs[i_nuclide]};
  int i_temp = micro.index_temp_sab;

  // Sample energy and angle
  double E_out;
  data::thermal_scatt[i_sab]->data_[i_temp].sample(micro, *E, &E_out, mu);

  // Set energy to outgoing, change direction of particle
  *E = E_out;
  rotate_angle_c(uvw, *mu, nullptr);
}

Direction sample_target_velocity(const Nuclide* nuc, double E, Direction u,
  Direction v_neut, double xs_eff, double kT, double* wgt)
{
  // check if nuclide is a resonant scatterer
  ResScatMethod sampling_method;
  if (nuc->resonant_) {

    // sampling method to use
    sampling_method = settings::res_scat_method;

    // upper resonance scattering energy bound (target is at rest above this E)
    if (E > settings::res_scat_energy_max) {
      return {};

    // lower resonance scattering energy bound (should be no resonances below)
    } else if (E < settings::res_scat_energy_min) {
      sampling_method = ResScatMethod::cxs;
    }

  // otherwise, use free gas model
  } else {
    if (E >= FREE_GAS_THRESHOLD * kT && nuc->awr_ > 1.0) {
      return {};
    } else {
      sampling_method = ResScatMethod::cxs;
    }
  }

  // use appropriate target velocity sampling method
  switch (sampling_method) {
  case ResScatMethod::cxs:

    // sample target velocity with the constant cross section (cxs) approx.
    return sample_cxs_target_velocity(nuc->awr_, E, u, kT);

  case ResScatMethod::dbrc:
  case ResScatMethod::rvs: {
    double E_red = std::sqrt(nuc->awr_ * E / kT);
    double E_low = std::pow(std::max(0.0, E_red - 4.0), 2) * kT / nuc->awr_;
    double E_up = (E_red + 4.0)*(E_red + 4.0) * kT / nuc->awr_;

    // find lower and upper energy bound indices
    // lower index
    int i_E_low;
    if (E_low < nuc->energy_0K_.front()) {
      i_E_low = 0;
    } else if (E_low > nuc->energy_0K_.back()) {
      i_E_low = nuc->energy_0K_.size() - 2;
    } else {
      i_E_low = lower_bound_index(nuc->energy_0K_.begin(),
        nuc->energy_0K_.end(), E_low);
    }

    // upper index
    int i_E_up;
    if (E_up < nuc->energy_0K_.front()) {
      i_E_up = 0;
    } else if (E_up > nuc->energy_0K_.back()) {
      i_E_up = nuc->energy_0K_.size() - 2;
    } else {
      i_E_up = lower_bound_index(nuc->energy_0K_.begin(),
        nuc->energy_0K_.end(), E_up);
    }

    if (i_E_up == i_E_low) {
      // Handle degenerate case -- if the upper/lower bounds occur for the same
      // index, then using cxs is probably a good approximation
      return sample_cxs_target_velocity(nuc->awr_, E, u, kT);
    }

    if (sampling_method == ResScatMethod::dbrc) {
      // interpolate xs since we're not exactly at the energy indices
      double xs_low = nuc->elastic_0K_[i_E_low];
      double m = (nuc->elastic_0K_[i_E_low + 1] - xs_low)
        / (nuc->energy_0K_[i_E_low + 1] - nuc->energy_0K_[i_E_low]);
      xs_low += m * (E_low - nuc->energy_0K_[i_E_low]);
      double xs_up = nuc->elastic_0K_[i_E_up];
      m = (nuc->elastic_0K_[i_E_up + 1] - xs_up)
        / (nuc->energy_0K_[i_E_up + 1] - nuc->energy_0K_[i_E_up]);
      xs_up += m * (E_up - nuc->energy_0K_[i_E_up]);

      // get max 0K xs value over range of practical relative energies
      double xs_max = *std::max_element(&nuc->elastic_0K_[i_E_low + 1],
        &nuc->elastic_0K_[i_E_up + 1]);
      xs_max = std::max({xs_low, xs_max, xs_up});

      while (true) {
        double E_rel;
        Direction v_target;
        while (true) {
          // sample target velocity with the constant cross section (cxs) approx.
          v_target = sample_cxs_target_velocity(nuc->awr_, E, u, kT);
          Direction v_rel = v_neut - v_target;
          E_rel = v_rel.dot(v_rel);
          if (E_rel < E_up) break;
        }

        // perform Doppler broadening rejection correction (dbrc)
        double xs_0K = nuc->elastic_xs_0K(E_rel);
        double R = xs_0K / xs_max;
        if (prn() < R) return v_target;
      }

    } else if (sampling_method == ResScatMethod::rvs) {
      // interpolate xs CDF since we're not exactly at the energy indices
      // cdf value at lower bound attainable energy
      double m = (nuc->xs_cdf_[i_E_low] - nuc->xs_cdf_[i_E_low - 1])
        / (nuc->energy_0K_[i_E_low + 1] - nuc->energy_0K_[i_E_low]);
      double cdf_low = nuc->xs_cdf_[i_E_low - 1]
            + m * (E_low - nuc->energy_0K_[i_E_low]);
      if (E_low <= nuc->energy_0K_.front()) cdf_low = 0.0;

      // cdf value at upper bound attainable energy
      m = (nuc->xs_cdf_[i_E_up] - nuc->xs_cdf_[i_E_up - 1])
        / (nuc->energy_0K_[i_E_up + 1] - nuc->energy_0K_[i_E_up]);
      double cdf_up = nuc->xs_cdf_[i_E_up - 1]
        + m*(E_up - nuc->energy_0K_[i_E_up]);

      while (true) {
        // directly sample Maxwellian
        double E_t = -kT * std::log(prn());

        // sample a relative energy using the xs cdf
        double cdf_rel = cdf_low + prn()*(cdf_up - cdf_low);
        int i_E_rel = lower_bound_index(&nuc->xs_cdf_[i_E_low-1],
          &nuc->xs_cdf_[i_E_up+1], cdf_rel);
        double E_rel = nuc->energy_0K_[i_E_low + i_E_rel];
        double m = (nuc->xs_cdf_[i_E_low + i_E_rel]
              - nuc->xs_cdf_[i_E_low + i_E_rel - 1])
              / (nuc->energy_0K_[i_E_low + i_E_rel + 1]
              -  nuc->energy_0K_[i_E_low + i_E_rel]);
        E_rel += (cdf_rel - nuc->xs_cdf_[i_E_low + i_E_rel - 1]) / m;

        // perform rejection sampling on cosine between
        // neutron and target velocities
        double mu = (E_t + nuc->awr_ * (E - E_rel)) /
          (2.0 * std::sqrt(nuc->awr_ * E * E_t));

        if (std::abs(mu) < 1.0) {
          // set and accept target velocity
          E_t /= nuc->awr_;
          return std::sqrt(E_t) * rotate_angle(u, mu, nullptr);
        }
      }
    }
  } // case RVS, DBRC
  } // switch (sampling_method)
}

Direction
sample_cxs_target_velocity(double awr, double E, Direction u, double kT)
{
  double beta_vn = std::sqrt(awr * E / kT);
  double alpha = 1.0/(1.0 + std::sqrt(PI)*beta_vn/2.0);

  double beta_vt_sq;
  double mu;
  while (true) {
    // Sample two random numbers
    double r1 = prn();
    double r2 = prn();

    if (prn() < alpha) {
      // With probability alpha, we sample the distribution p(y) =
      // y*e^(-y). This can be done with sampling scheme C45 frmo the Monte
      // Carlo sampler

      beta_vt_sq = -std::log(r1*r2);

    } else {
      // With probability 1-alpha, we sample the distribution p(y) = y^2 *
      // e^(-y^2). This can be done with sampling scheme C61 from the Monte
      // Carlo sampler

      double c = std::cos(PI/2.0 * prn());
      beta_vt_sq = -std::log(r1) - std::log(r2)*c*c;
    }

    // Determine beta * vt
    double beta_vt = std::sqrt(beta_vt_sq);

    // Sample cosine of angle between neutron and target velocity
    mu = 2.0*prn() - 1.0;

    // Determine rejection probability
    double accept_prob = std::sqrt(beta_vn*beta_vn + beta_vt_sq -
      2*beta_vn*beta_vt*mu) / (beta_vn + beta_vt);

    // Perform rejection sampling on vt and mu
    if (prn() < accept_prob) break;
  }

  // Determine speed of target nucleus
  double vt = std::sqrt(beta_vt_sq*kT/awr);

  // Determine velocity vector of target nucleus based on neutron's velocity
  // and the sampled angle between them
  return vt * rotate_angle(u, mu, nullptr);
}

void sample_fission_neutron(int i_nuclide, const Reaction* rx, double E_in, Bank* site)
{
  // Sample cosine of angle -- fission neutrons are always emitted
  // isotropically. Sometimes in ACE data, fission reactions actually have
  // an angular distribution listed, but for those that do, it's simply just
  // a uniform distribution in mu
  double mu = 2.0 * prn() - 1.0;

  // Sample azimuthal angle uniformly in [0,2*pi)
  double phi = 2.0*PI*prn();
  site->uvw[0] = mu;
  site->uvw[1] = std::sqrt(1.0 - mu*mu) * std::cos(phi);
  site->uvw[2] = std::sqrt(1.0 - mu*mu) * std::sin(phi);

  // Determine total nu, delayed nu, and delayed neutron fraction
  const auto& nuc {data::nuclides[i_nuclide]};
  double nu_t = nuc->nu(E_in, Nuclide::EmissionMode::total);
  double nu_d = nuc->nu(E_in, Nuclide::EmissionMode::delayed);
  double beta = nu_d / nu_t;

  if (prn() < beta) {
    // ====================================================================
    // DELAYED NEUTRON SAMPLED

    // sampled delayed precursor group
    double xi = prn()*nu_d;
    double prob = 0.0;
    int group;
    for (group = 1; group < nuc->n_precursor_; ++group) {
      // determine delayed neutron precursor yield for group j
      double yield = (*rx->products_[group].yield_)(E_in);

      // Check if this group is sampled
      prob += yield;
      if (xi < prob) break;
    }

    // if the sum of the probabilities is slightly less than one and the
    // random number is greater, j will be greater than nuc %
    // n_precursor -- check for this condition
    group = std::min(group, nuc->n_precursor_);

    // set the delayed group for the particle born from fission
    site->delayed_group = group;

    int n_sample = 0;
    while (true) {
      // sample from energy/angle distribution -- note that mu has already been
      // sampled above and doesn't need to be resampled
      rx->products_[group].sample(E_in, site->E, mu);

      // resample if energy is greater than maximum neutron energy
      // TODO: off-by-one
      constexpr int neutron = static_cast<int>(ParticleType::neutron) - 1;
      if (site->E < data::energy_max[neutron]) break;

      // check for large number of resamples
      ++n_sample;
      if (n_sample == MAX_SAMPLE) {
        // particle_write_restart(p)
        fatal_error("Resampled energy distribution maximum number of times "
          "for nuclide " + nuc->name_);
      }
    }

  } else {
    // ====================================================================
    // PROMPT NEUTRON SAMPLED

    // set the delayed group for the particle born from fission to 0
    site->delayed_group = 0;

    // sample from prompt neutron energy distribution
    int n_sample = 0;
    while (true) {
      rx->products_[0].sample(E_in, site->E, mu);

      // resample if energy is greater than maximum neutron energy
      // TODO: off-by-one
      constexpr int neutron = static_cast<int>(ParticleType::neutron) - 1;
      if (site->E < data::energy_max[neutron]) break;

      // check for large number of resamples
      ++n_sample;
      if (n_sample == MAX_SAMPLE) {
        // particle_write_restart(p)
        fatal_error("Resampled energy distribution maximum number of times "
          "for nuclide " + nuc->name_);
      }
    }
  }
}

void inelastic_scatter(const Nuclide* nuc, const Reaction* rx, Particle* p)
{
  // copy energy of neutron
  double E_in = p->E;

  // sample outgoing energy and scattering cosine
  double E;
  double mu;
  rx->products_[0].sample(E_in, E, mu);

  // if scattering system is in center-of-mass, transfer cosine of scattering
  // angle and outgoing energy from CM to LAB
  if (rx->scatter_in_cm_) {
    double E_cm = E;

    // determine outgoing energy in lab
    double A = nuc->awr_;
    E = E_cm + (E_in + 2.0*mu*(A + 1.0) * std::sqrt(E_in*E_cm))
          / ((A + 1.0)*(A + 1.0));

    // determine outgoing angle in lab
    mu = mu*std::sqrt(E_cm/E) + 1.0/(A+1.0) * std::sqrt(E_in/E);
  }

  // Because of floating-point roundoff, it may be possible for mu to be
  // outside of the range [-1,1). In these cases, we just set mu to exactly -1
  // or 1
  if (std::abs(mu) > 1.0) mu = std::copysign(1.0, mu);

  // Set outgoing energy and scattering angle
  p->E = E;
  p->mu = mu;

  // change direction of particle
  rotate_angle_c(p->coord[0].uvw, mu, nullptr);

  // evaluate yield
  double yield = (*rx->products_[0].yield_)(E_in);
  if (std::floor(yield) == yield) {
    // If yield is integral, create exactly that many secondary particles
    for (int i = 0; i < static_cast<int>(std::round(yield)) - 1; ++i) {
      int neutron = static_cast<int>(ParticleType::neutron);
      p->create_secondary(p->coord[0].uvw, p->E, neutron, true);
    }
  } else {
    // Otherwise, change weight of particle based on yield
    p->wgt *= yield;
  }
}

void sample_secondary_photons(Particle* p, int i_nuclide)
{
  // Sample the number of photons produced
  double y_t = p->wgt * simulation::micro_xs[i_nuclide].photon_prod /
    simulation::micro_xs[i_nuclide].total;
  int y = static_cast<int>(y_t);
  if (prn() <= y_t - y) ++y;

  // Sample each secondary photon
  for (int i = 0; i < y; ++i) {
    // Sample the reaction and product
    int i_rx;
    int i_product;
    sample_photon_product(i_nuclide, p->E, &i_rx, &i_product);

    // Sample the outgoing energy and angle
    auto& rx = data::nuclides[i_nuclide]->reactions_[i_rx];
    double E;
    double mu;
    rx->products_[i_product].sample(p->E, E, mu);

    // Sample the new direction
    double uvw[3];
    std::copy(p->coord[0].uvw, p->coord[0].uvw + 3, uvw);
    rotate_angle_c(uvw, mu, nullptr);

    // Create the secondary photon
    int photon = static_cast<int>(ParticleType::photon);
    p->create_secondary(uvw, E, photon, true);
  }
}

} // namespace openmc
