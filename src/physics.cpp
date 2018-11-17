#include "openmc/physics.h"

#include "openmc/bank.h"
#include "openmc/constants.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/physics_common.h"
#include "openmc/random_lcg.h"
#include "openmc/reaction.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

#include <algorithm> // for max, min
#include <cmath> // for sqrt, exp, log
#include <sstream>

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

void collision(Particle* p)
{
  // Add to collision counter for particle
  ++p->n_collision;

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
  p->event_nuclide = i_nuclide;

  // Create fission bank sites. Note that while a fission reaction is sampled,
  // it never actually "happens", i.e. the weight of the particle does not
  // change when sampling fission sites. The following block handles all
  // absorption (including fission)

  const auto& nuc {data::nuclides[i_nuclide-1]};

  if (nuc->fissionable_) {
    int i_rx = sample_fission(i_nuclide, p->E);
    if (settings::run_mode == RUN_MODE_EIGENVALUE) {
      create_fission_sites(p, i_nuclide, i_rx, simulation::fission_bank.data(),
        &simulation::n_bank, simulation::fission_bank.size());
    } else if (settings::run_mode == RUN_MODE_FIXEDSOURCE &&
      settings::create_fission_neutrons) {
      create_fission_sites(p, i_nuclide, i_rx, p->secondary_bank,
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

  if (simulation::micro_xs[i_nuclide-1].absorption > 0.0) {
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
create_fission_sites(Particle* p, int i_nuclide, int i_rx, Bank* bank_array,
  int64_t* size_bank, int64_t bank_capacity)
{
  // TODO: Heat generation from fission

  // If uniform fission source weighting is turned on, we increase or decrease
  // the expected number of fission sites produced
  double weight = settings::ufs_on ? ufs_get_weight(p) : 1.0;

  // Determine the expected number of neutrons produced
  double nu_t = p->wgt / simulation::keff * weight * simulation::micro_xs[
    i_nuclide-1].nu_fission / simulation::micro_xs[i_nuclide-1].total;

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
    sample_fission_neutron(i_nuclide, i_rx, p->E, &bank_array[i]);

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
//             / std::sqrt(alpha**2 + alpha_out**2 - TWO*alpha*alpha_out*mu)
//       phi = TWO*PI*prn()
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
//               rel_vel = std::sqrt(E_electron * (E_electron + TWO * MASS_ELECTRON_EV))&
//                     / (E_electron + MASS_ELECTRON_EV)
//               mu = (TWO * r + rel_vel - 1.0) / &
//                     (TWO * rel_vel * r - rel_vel + 1.0)
//               exit SAMPLE_MU
//             }
//           end do SAMPLE_MU

//           phi = TWO*PI*prn()
//           uvw(1) = mu
//           uvw(2) = std::sqrt(1.0 - mu*mu)*cos(phi)
//           uvw(3) = std::sqrt(1.0 - mu*mu)*sin(phi)

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

// void sample_electron_reaction(Particle* p)
// {
//   // TODO: create reaction types

//   if (electron_treatment == ELECTRON_TTB) {
//     thick_target_bremsstrahlung(p, E_lost)
//   }

//   p->E = 0.0
//   p->alive = false
// }

// void sample_positron_reaction(Particle* p)
// {
//   // TODO: create reaction types

//   if (electron_treatment == ELECTRON_TTB) {
//     thick_target_bremsstrahlung(p, E_lost)
//   }

//   // Sample angle isotropically
//   mu = TWO*prn() - 1.0
//   phi = TWO*PI*prn()
//   uvw(1) = mu
//   uvw(2) = std::sqrt(1.0 - mu*mu)*cos(phi)
//   uvw(3) = std::sqrt(1.0 - mu*mu)*sin(phi)

//   // Create annihilation photon pair traveling in opposite directions
//   particle_create_secondary(p, uvw, MASS_ELECTRON_EV, PHOTON, true)
//   particle_create_secondary(p, -uvw, MASS_ELECTRON_EV, PHOTON, true)

//   p->E = 0.0
//   p->alive = false
// }

// void sample_nuclide(Particle* p, int mt, int i_nuclide, int i_nuc_mat)
// {
//   // Get pointer to current material
//   mat => materials(p->material)

//   // Sample cumulative distribution function
//   select case (base)
//   case ('total')
//     cutoff = prn() * material_xs % total
//   case ('scatter')
//     cutoff = prn() * (material_xs % total - material_xs % absorption)
//   case ('fission')
//     cutoff = prn() * material_xs % fission
//   end select

//   i_nuc_mat = 0
//   prob = 0.0
//   do while (prob < cutoff)
//     i_nuc_mat = i_nuc_mat + 1

//     // Check to make sure that a nuclide was sampled
//     if (i_nuc_mat > mat % n_nuclides) {
//       particle_write_restart(p)
//       fatal_error("Did not sample any nuclide during collision.")
//     }

//     // Find atom density
//     i_nuclide    = mat % nuclide(i_nuc_mat)
//     atom_density = mat % atom_density(i_nuc_mat)

//     // Determine microscopic cross section
//     select case (base)
//     case ('total')
//       sigma = atom_density * micro_xs(i_nuclide) % total
//     case ('scatter')
//       sigma = atom_density * (micro_xs(i_nuclide) % total - &
//             micro_xs(i_nuclide) % absorption)
//     case ('fission')
//       sigma = atom_density * micro_xs(i_nuclide) % fission
//     end select

//     // Increment probability to compare to cutoff
//     prob = prob + sigma
//   end do
// }

// void sample_element(Particle* p)
// {
//   associate (mat => materials(p->material))
//     // Sample cumulative distribution function
//     cutoff = prn() * material_xs % total

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

// int sample_fission(int i_nuclide, double E)
// {
//   // Get pointer to nuclide
//   nuc => nuclides(i_nuclide)

//   // If we're in the URR, by default use the first fission reaction. We also
//   // default to the first reaction if we know that there are no partial fission
//   // reactions
//   if (micro_xs(i_nuclide) % use_ptable || &
//         !nuc % has_partial_fission) {
//     i_reaction = nuc % index_fission(1)
//     return
//   }

//   // Check to see if we are in a windowed multipole range.  WMP only supports
//   // the first fission reaction.
//   if (nuc % mp_present) {
//     if (E >= nuc % multipole % E_min && &
//           E <= nuc % multipole % E_max) {
//       i_reaction = nuc % index_fission(1)
//       return
//     }
//   }

//   // Get grid index and interpolatoin factor and sample fission cdf
//   i_temp = micro_xs(i_nuclide) % index_temp
//   i_grid = micro_xs(i_nuclide) % index_grid
//   f      = micro_xs(i_nuclide) % interp_factor
//   cutoff = prn() * micro_xs(i_nuclide) % fission
//   prob   = 0.0

//   // Loop through each partial fission reaction type

//   FISSION_REACTION_LOOP: do i = 1, nuc % n_fission
//     i_reaction = nuc % index_fission(i)

//     associate (rx => nuc % reactions(i_reaction))
//       // if energy is below threshold for this reaction, skip it
//       threshold = rx % xs_threshold(i_temp)
//       if (i_grid < threshold) cycle

//       // add to cumulative probability
//       prob = prob + ((1.0 - f) * rx % xs(i_temp, i_grid - threshold + 1) &
//             + f*(rx % xs(i_temp, i_grid - threshold + 2)))
//     end associate

//     // Create fission bank sites if fission occurs
//     if (prob > cutoff) exit FISSION_REACTION_LOOP
//   end do FISSION_REACTION_LOOP

// end subroutine sample_fission
// }

// void sample_photon_product(int i_nuclide, double E, int* i_rx, int* i_product);
// {
//   // Get pointer to nuclide
//   associate (nuc => nuclides(i_nuclide))

//     // Get grid index and interpolation factor and sample photon production cdf
//     i_temp = micro_xs(i_nuclide) % index_temp
//     i_grid = micro_xs(i_nuclide) % index_grid
//     f      = micro_xs(i_nuclide) % interp_factor
//     cutoff = prn() * micro_xs(i_nuclide) % photon_prod
//     prob   = 0.0

//     // Loop through each reaction type
//     REACTION_LOOP: do i_reaction = 1, size(nuc % reactions)
//       associate (rx => nuc % reactions(i_reaction))
//         threshold = rx % xs_threshold(i_temp)

//         // if energy is below threshold for this reaction, skip it
//         if (i_grid < threshold) cycle

//         do i_product = 1, rx % products_size()
//           if (rx % product_particle(i_product) == PHOTON) {
//             // add to cumulative probability
//             yield = rx % product_yield(i_product, E)
//             prob = prob + ((1.0 - f) * rx % xs(i_temp, i_grid - threshold + 1) &
//                   + f*(rx % xs(i_temp, i_grid - threshold + 2))) * yield

//             if (prob > cutoff) return
//             last_valid_reaction = i_reaction
//             last_valid_product = i_product
//           }
//         end do
//       end associate
//     end do REACTION_LOOP
//   end associate

//   i_reaction = last_valid_reaction
//   i_product = last_valid_product
// }

// void absorption(Particle* p, int i_nuclide)
// {
//   if (survival_biasing) {
//     // Determine weight absorbed in survival biasing
//     p->absorb_wgt = p->wgt * micro_xs(i_nuclide) % absorption / &
//           micro_xs(i_nuclide) % total

//     // Adjust weight of particle by probability of absorption
//     p->wgt = p->wgt - p->absorb_wgt
//     p->last_wgt = p->wgt

//     // Score implicit absorption estimate of keff
//     if (run_mode == MODE_EIGENVALUE) {
//       global_tally_absorption = global_tally_absorption + p->absorb_wgt * &
//             micro_xs(i_nuclide) % nu_fission / micro_xs(i_nuclide) % absorption
//     }
//   } else {
//     // See if disappearance reaction happens
//     if (micro_xs(i_nuclide) % absorption > &
//           prn() * micro_xs(i_nuclide) % total) {
//       // Score absorption estimate of keff
//       if (run_mode == MODE_EIGENVALUE) {
//         global_tally_absorption = global_tally_absorption + p->wgt * &
//               micro_xs(i_nuclide) % nu_fission / micro_xs(i_nuclide) % absorption
//       }

//       p->alive = false
//       p->event = EVENT_ABSORB
//       p->event_MT = N_DISAPPEAR
//     }
//   }
// }

// void scatter(Particle*, int i_nuclide, int i_nuc_mat)
// {
//   // copy incoming direction
//   uvw_old(:) = p->coord(1) % uvw

//   // Get pointer to nuclide and grid index/interpolation factor
//   nuc    => nuclides(i_nuclide)
//   i_temp =  micro_xs(i_nuclide) % index_temp
//   i_grid =  micro_xs(i_nuclide) % index_grid
//   f      =  micro_xs(i_nuclide) % interp_factor

//   // For tallying purposes, this routine might be called directly. In that
//   // case, we need to sample a reaction via the cutoff variable
//   cutoff = prn() * (micro_xs(i_nuclide) % total - &
//         micro_xs(i_nuclide) % absorption)
//   sampled = false

//   // Calculate elastic cross section if it wasn't precalculated
//   if (micro_xs(i_nuclide) % elastic == CACHE_INVALID) {
//     nuc % calculate_elastic_xs(micro_xs(i_nuclide))
//   }

//   prob = micro_xs(i_nuclide) % elastic - micro_xs(i_nuclide) % thermal
//   if (prob > cutoff) {
//     // =======================================================================
//     // NON-S(A,B) ELASTIC SCATTERING

//     // Determine temperature
//     if (nuc % mp_present) {
//       kT = p->sqrtkT**2
//     } else {
//       kT = nuc % kTs(micro_xs(i_nuclide) % index_temp)
//     }

//     // Perform collision physics for elastic scattering
//     elastic_scatter(i_nuclide, nuc % reactions(1), kT, p->E, &
//                           p->coord(1) % uvw, p->mu, p->wgt)

//     p->event_MT = ELASTIC
//     sampled = true
//   }

//   prob = micro_xs(i_nuclide) % elastic
//   if (prob > cutoff && !sampled) {
//     // =======================================================================
//     // S(A,B) SCATTERING

//     sab_scatter(i_nuclide, micro_xs(i_nuclide) % index_sab, p->E, &
//                       p->coord(1) % uvw, p->mu)

//     p->event_MT = ELASTIC
//     sampled = true
//   }

//   if (!sampled) {
//     // =======================================================================
//     // INELASTIC SCATTERING

//     j = 0
//     do while (prob < cutoff)
//       j = j + 1
//       i = nuc % index_inelastic_scatter(j)

//       // Check to make sure inelastic scattering reaction sampled
//       if (i > size(nuc % reactions)) {
//         particle_write_restart(p)
//         fatal_error("Did not sample any reaction for nuclide " &
//               &// trim(nuc % name))
//       }

//       associate (rx => nuc % reactions(i))
//         // if energy is below threshold for this reaction, skip it
//         threshold = rx % xs_threshold(i_temp)
//         if (i_grid < threshold) cycle

//         // add to cumulative probability
//         prob = prob + ((1.0 - f)*rx % xs(i_temp, i_grid - threshold + 1) &
//               + f*(rx % xs(i_temp, i_grid - threshold + 2)))
//       end associate
//     end do

//     // Perform collision physics for inelastic scattering
//     inelastic_scatter(nuc, nuc%reactions(i), p)
//     p->event_MT = nuc % reactions(i) % MT

//   }

//   // Set event component
//   p->event = EVENT_SCATTER

//   // Sample new outgoing angle for isotropic-in-lab scattering
//   associate (mat => materials(p->material))
//     if (mat % has_isotropic_nuclides) {
//       if (materials(p->material) % p0(i_nuc_mat)) {
//         // Sample isotropic-in-lab outgoing direction
//         uvw_new(1) = TWO * prn() - 1.0
//         phi = TWO * PI * prn()
//         uvw_new(2) = cos(phi) * std::sqrt(1.0 - uvw_new(1)*uvw_new(1))
//         uvw_new(3) = sin(phi) * std::sqrt(1.0 - uvw_new(1)*uvw_new(1))
//         p->mu = dot_product(uvw_old, uvw_new)

//         // Change direction of particle
//         p->coord(1) % uvw = uvw_new
//       }
//     }
//   end associate
// }

// void elastic_scatter(int i_nuclide, const Reaction& rx, double kT, double* E,
//   Direction& u, double* mu_lab, double* wgt)
// {
//   // get pointer to nuclide
//   nuc => nuclides(i_nuclide)

//   vel = std::sqrt(E)
//   awr = nuc % awr

//   // Neutron velocity in LAB
//   v_n = vel * uvw

//   // Sample velocity of target nucleus
//   if (!micro_xs(i_nuclide) % use_ptable) {
//     sample_target_velocity(nuc, v_t, E, uvw, v_n, wgt, &
//           micro_xs(i_nuclide) % elastic, kT)
//   } else {
//     v_t = 0.0
//   }

//   // Velocity of center-of-mass
//   v_cm = (v_n + awr*v_t)/(awr + 1.0)

//   // Transform to CM frame
//   v_n = v_n - v_cm

//   // Find speed of neutron in CM
//   vel = std::sqrt(dot_product(v_n, v_n))

//   // Sample scattering angle
//   mu_cm = rxn % sample_elastic_mu(E)

//   // Determine direction cosines in CM
//   uvw_cm = v_n/vel

//   // Rotate neutron velocity vector to new angle -- note that the speed of the
//   // neutron in CM does not change in elastic scattering. However, the speed
//   // will change when we convert back to LAB
//   v_n = vel * rotate_angle(uvw_cm, mu_cm)

//   // Transform back to LAB frame
//   v_n = v_n + v_cm

//   E = dot_product(v_n, v_n)
//   vel = std::sqrt(E)

//   // compute cosine of scattering angle in LAB frame by taking dot product of
//   // neutron's pre- and post-collision angle
//   mu_lab = dot_product(uvw, v_n) / vel

//   // Set energy and direction of particle in LAB frame
//   uvw = v_n / vel

//   // Because of floating-point roundoff, it may be possible for mu_lab to be
//   // outside of the range [-1,1). In these cases, we just set mu_lab to exactly
//   // -1 or 1

//   if (abs(mu_lab) > 1.0) mu_lab = sign(1.0,mu_lab)
// }

// void sab_scatter(int i_nuclide, int i_sab, double* E, Direction* u, double* mu)
// {
//   // Sample from C++ side
//   ptr = C_LOC(micro_xs(i_nuclide))
//   sab_tables(i_sab) % sample(ptr, E, E_out, mu)

//   // Set energy to outgoing, change direction of particle
//   E = E_out
//   uvw = rotate_angle(uvw, mu)
// }

// void sample_target_velocity(int i_nuclide, Direction* v_target, double E, Direction u,
//   Direction v_neut, double* wgt, double xs_eff, double kT)
// {
//   awr = nuc % awr

//   // check if nuclide is a resonant scatterer
//   if (nuc % resonant) {

//     // sampling method to use
//     sampling_method = res_scat_method

//     // upper resonance scattering energy bound (target is at rest above this E)
//     if (E > res_scat_energy_max) {
//       v_target = 0.0
//       return

//     // lower resonance scattering energy bound (should be no resonances below)
//     } else if (E < res_scat_energy_min) {
//       sampling_method = RES_SCAT_CXS
//     }

//   // otherwise, use free gas model
//   } else {
//     if (E >= FREE_GAS_THRESHOLD * kT && awr > 1.0) {
//       v_target = 0.0
//       return
//     } else {
//       sampling_method = RES_SCAT_CXS
//     }
//   }

//   // use appropriate target velocity sampling method
//   select case (sampling_method)
//   case (RES_SCAT_CXS)

//     // sample target velocity with the constant cross section (cxs) approx.
//     sample_cxs_target_velocity(nuc, v_target, E, uvw, kT)

//   case (RES_SCAT_WCM)

//     // sample target velocity with the constant cross section (cxs) approx.
//     sample_cxs_target_velocity(nuc, v_target, E, uvw, kT)

//     // adjust weight as prescribed by the weight correction method (wcm)
//     E_rel = dot_product((v_neut - v_target), (v_neut - v_target))
//     xs_0K = elastic_xs_0K(E_rel, nuc)
//     wcf = xs_0K / xs_eff
//     wgt = wcf * wgt

//   case (RES_SCAT_DBRC, RES_SCAT_ARES)
//     E_red = std::sqrt(awr * E / kT)
//     E_low = std::max(0.0, E_red - FOUR)**2 * kT / awr
//     E_up  = (E_red + FOUR)**2 * kT / awr

//     // find lower and upper energy bound indices
//     // lower index
//     n_grid = size(nuc % energy_0K)
//     if (E_low < nuc % energy_0K(1)) {
//       i_E_low = 1
//     } else if (E_low > nuc % energy_0K(n_grid)) {
//       i_E_low = n_grid - 1
//     } else {
//       i_E_low = binary_search(nuc % energy_0K, n_grid, E_low)
//     }

//     // upper index
//     if (E_up < nuc % energy_0K(1)) {
//       i_E_up = 1
//     } else if (E_up > nuc % energy_0K(n_grid)) {
//       i_E_up = n_grid - 1
//     } else {
//       i_E_up = binary_search(nuc % energy_0K, n_grid, E_up)
//     }

//     if (i_E_up == i_E_low) {
//       // Handle degenerate case -- if the upper/lower bounds occur for the same
//       // index, then using cxs is probably a good approximation
//       sample_cxs_target_velocity(nuc, v_target, E, uvw, kT)

//     } else {
//       if (sampling_method == RES_SCAT_DBRC) {
//         // interpolate xs since we're not exactly at the energy indices
//         xs_low = nuc % elastic_0K(i_E_low)
//         m = (nuc % elastic_0K(i_E_low + 1) - xs_low) &
//               / (nuc % energy_0K(i_E_low + 1) - nuc % energy_0K(i_E_low))
//         xs_low = xs_low + m * (E_low - nuc % energy_0K(i_E_low))
//         xs_up = nuc % elastic_0K(i_E_up)
//         m = (nuc % elastic_0K(i_E_up + 1) - xs_up) &
//               / (nuc % energy_0K(i_E_up + 1) - nuc % energy_0K(i_E_up))
//         xs_up = xs_up + m * (E_up - nuc % energy_0K(i_E_up))

//         // get max 0K xs value over range of practical relative energies
//         xs_max = std::max(xs_low, &
//               maxval(nuc % elastic_0K(i_E_low + 1 : i_E_up)), xs_up)

//         DBRC_REJECT_LOOP: do
//           TARGET_ENERGY_LOOP: do
//             // sample target velocity with the constant cross section (cxs) approx.
//             sample_cxs_target_velocity(nuc, v_target, E, uvw, kT)
//             E_rel = dot_product((v_neut - v_target), (v_neut - v_target))
//             if (E_rel < E_up) exit TARGET_ENERGY_LOOP
//           end do TARGET_ENERGY_LOOP

//           // perform Doppler broadening rejection correction (dbrc)
//           xs_0K = elastic_xs_0K(E_rel, nuc)
//           R = xs_0K / xs_max
//           if (prn() < R) exit DBRC_REJECT_LOOP
//         end do DBRC_REJECT_LOOP

//       } else if (sampling_method == RES_SCAT_ARES) {
//         // interpolate xs CDF since we're not exactly at the energy indices
//         // cdf value at lower bound attainable energy
//         m = (nuc % xs_cdf(i_E_low) - nuc % xs_cdf(i_E_low - 1)) &
//               / (nuc % energy_0K(i_E_low + 1) - nuc % energy_0K(i_E_low))
//         cdf_low = nuc % xs_cdf(i_E_low - 1) &
//               + m * (E_low - nuc % energy_0K(i_E_low))
//         if (E_low <= nuc % energy_0K(1)) cdf_low = 0.0

//         // cdf value at upper bound attainable energy
//         m = (nuc % xs_cdf(i_E_up) - nuc % xs_cdf(i_E_up - 1)) &
//               / (nuc % energy_0K(i_E_up + 1) - nuc % energy_0K(i_E_up))
//         cdf_up = nuc % xs_cdf(i_E_up - 1) &
//               + m * (E_up - nuc % energy_0K(i_E_up))

//         ARES_REJECT_LOOP: do

//           // directly sample Maxwellian
//           E_t = -kT * std::log(prn())

//           // sample a relative energy using the xs cdf
//           cdf_rel = cdf_low + prn() * (cdf_up - cdf_low)
//           i_E_rel = binary_search(nuc % xs_cdf(i_E_low-1:i_E_up), &
//                 i_E_up - i_E_low + 2, cdf_rel)
//           E_rel = nuc % energy_0K(i_E_low + i_E_rel - 1)
//           m = (nuc % xs_cdf(i_E_low + i_E_rel - 1) &
//                 - nuc % xs_cdf(i_E_low + i_E_rel - 2)) &
//                 / (nuc % energy_0K(i_E_low + i_E_rel) &
//                 -  nuc % energy_0K(i_E_low + i_E_rel - 1))
//           E_rel = E_rel + (cdf_rel - nuc % xs_cdf(i_E_low + i_E_rel - 2)) / m

//           // perform rejection sampling on cosine between
//           // neutron and target velocities
//           mu = (E_t + awr * (E - E_rel)) / (TWO * std::sqrt(awr * E * E_t))

//           if (abs(mu) < 1.0) {
//             // set and accept target velocity
//             E_t = E_t / awr
//             v_target = std::sqrt(E_t) * rotate_angle(uvw, mu)
//             exit ARES_REJECT_LOOP
//           }
//         end do ARES_REJECT_LOOP
//       }
//     }
//   end select
// }

// void sample_cxs_target_velocity(int i_nuclide, Direction* v_target, double E, Direction u,
//   double kT)
// {
//   awr = nuc % awr

//   beta_vn = std::sqrt(awr * E / kT)
//   alpha = 1.0/(1.0 + std::sqrt(pi)*beta_vn/TWO)

//   do
//     // Sample two random numbers
//     r1 = prn()
//     r2 = prn()

//     if (prn() < alpha) {
//       // With probability alpha, we sample the distribution p(y) =
//       // y*e^(-y). This can be done with sampling scheme C45 frmo the Monte
//       // Carlo sampler

//       beta_vt_sq = -std::log(r1*r2)

//     } else {
//       // With probability 1-alpha, we sample the distribution p(y) = y^2 *
//       // e^(-y^2). This can be done with sampling scheme C61 from the Monte
//       // Carlo sampler

//       c = cos(PI/TWO * prn())
//       beta_vt_sq = -std::log(r1) - std::log(r2)*c*c
//     }

//     // Determine beta * vt
//     beta_vt = std::sqrt(beta_vt_sq)

//     // Sample cosine of angle between neutron and target velocity
//     mu = TWO*prn() - 1.0

//     // Determine rejection probability
//     accept_prob = std::sqrt(beta_vn*beta_vn + beta_vt_sq - 2*beta_vn*beta_vt*mu) &
//           /(beta_vn + beta_vt)

//     // Perform rejection sampling on vt and mu
//     if (prn() < accept_prob) exit
//   end do

//   // Determine speed of target nucleus
//   vt = std::sqrt(beta_vt_sq*kT/awr)

//   // Determine velocity vector of target nucleus based on neutron's velocity
//   // and the sampled angle between them
//   v_target = vt * rotate_angle(uvw, mu)
// }

// void sample_fission_neutron(int i_nuclide, const Reaction& rx, double E_in, Bank* site)
// {
//   // Sample cosine of angle -- fission neutrons are always emitted
//   // isotropically. Sometimes in ACE data, fission reactions actually have
//   // an angular distribution listed, but for those that do, it's simply just
//   // a uniform distribution in mu
//   mu = TWO * prn() - 1.0

//   // Sample azimuthal angle uniformly in [0,2*pi)
//   phi = TWO*PI*prn()
//   site % uvw(1) = mu
//   site % uvw(2) = std::sqrt(1.0 - mu*mu) * cos(phi)
//   site % uvw(3) = std::sqrt(1.0 - mu*mu) * sin(phi)

//   // Determine total nu, delayed nu, and delayed neutron fraction
//   nu_t = nuc % nu(E_in, EMISSION_TOTAL)
//   nu_d = nuc % nu(E_in, EMISSION_DELAYED)
//   beta = nu_d / nu_t

//   if (prn() < beta) {
//     // ====================================================================
//     // DELAYED NEUTRON SAMPLED

//     // sampled delayed precursor group
//     xi = prn()*nu_d
//     prob = 0.0
//     do group = 1, nuc % n_precursor

//       // determine delayed neutron precursor yield for group j
//       yield = rxn % product_yield(1 + group, E_in)

//       // Check if this group is sampled
//       prob = prob + yield
//       if (xi < prob) exit
//     end do

//     // if the sum of the probabilities is slightly less than one and the
//     // random number is greater, j will be greater than nuc %
//     // n_precursor -- check for this condition
//     group = min(group, nuc % n_precursor)

//     // set the delayed group for the particle born from fission
//     site % delayed_group = group

//     n_sample = 0
//     do
//       // sample from energy/angle distribution -- note that mu has already been
//       // sampled above and doesn't need to be resampled
//       rxn % product_sample(1 + group, E_in, site % E, mu)

//       // resample if energy is greater than maximum neutron energy
//       if (site % E < energy_max(NEUTRON)) exit

//       // check for large number of resamples
//       n_sample = n_sample + 1
//       if (n_sample == MAX_SAMPLE) {
//         // particle_write_restart(p)
//         fatal_error("Resampled energy distribution maximum number of " &
//               // "times for nuclide " // nuc % name)
//       }
//     end do

//   } else {
//     // ====================================================================
//     // PROMPT NEUTRON SAMPLED

//     // set the delayed group for the particle born from fission to 0
//     site % delayed_group = 0

//     // sample from prompt neutron energy distribution
//     n_sample = 0
//     do
//       rxn % product_sample(1, E_in, site % E, mu)

//       // resample if energy is greater than maximum neutron energy
//       if (site % E < energy_max(NEUTRON)) exit

//       // check for large number of resamples
//       n_sample = n_sample + 1
//       if (n_sample == MAX_SAMPLE) {
//         // particle_write_restart(p)
//         fatal_error("Resampled energy distribution maximum number of " &
//               // "times for nuclide " // nuc % name)
//       }
//     end do
//   }
// }

// void inelastic_scatter(int i_nuclide, const Reaction& rx, Particle* p)
// {
//   // copy energy of neutron
//   E_in = p->E

//   // sample outgoing energy and scattering cosine
//   rxn % product_sample(1, E_in, E, mu)

//   // if scattering system is in center-of-mass, transfer cosine of scattering
//   // angle and outgoing energy from CM to LAB
//   if (rxn % scatter_in_cm) {
//     E_cm = E

//     // determine outgoing energy in lab
//     A = nuc%awr
//     E = E_cm + (E_in + TWO * mu * (A+1.0) * std::sqrt(E_in * E_cm)) &
//           / ((A+1.0)*(A+1.0))

//     // determine outgoing angle in lab
//     mu = mu * std::sqrt(E_cm/E) + 1.0/(A+1.0) * std::sqrt(E_in/E)
//   }

//   // Because of floating-point roundoff, it may be possible for mu to be
//   // outside of the range [-1,1). In these cases, we just set mu to exactly -1
//   // or 1
//   if (abs(mu) > 1.0) mu = sign(1.0,mu)

//   // Set outgoing energy and scattering angle
//   p->E = E
//   p->mu = mu

//   // change direction of particle
//   p->coord(1) % uvw = rotate_angle(p->coord(1) % uvw, mu)

//   // evaluate yield
//   yield = rxn % product_yield(1, E_in)
//   if (mod(yield, 1.0) == 0.0) {
//     // If yield is integral, create exactly that many secondary particles
//     do i = 1, nint(yield) - 1
//       particle_create_secondary(p, p->coord(1) % uvw, p->E, &
//             NEUTRON, run_CE=true)
//     end do
//   } else {
//     // Otherwise, change weight of particle based on yield
//     p->wgt = yield * p->wgt
//   }
// }

// void sample_secondary_photons(Particle* p, int i_nuclide)
// {
//   // Sample the number of photons produced
//   nu_t =  p->wgt * micro_xs(i_nuclide) % photon_prod / &
//         micro_xs(i_nuclide) % total
//   if (prn() > nu_t - int(nu_t)) {
//     nu = int(nu_t)
//   } else {
//     nu = int(nu_t) + 1
//   }

//   // Sample each secondary photon
//   do i = 1, nu

//     // Sample the reaction and product
//     sample_photon_product(i_nuclide, p->E, i_reaction, i_product)

//     // Sample the outgoing energy and angle
//     nuclides(i_nuclide) % reactions(i_reaction) % &
//           product_sample(i_product, p->E, E, mu)

//     // Sample the new direction
//     uvw = rotate_angle(p->coord(1) % uvw, mu)

//     // Create the secondary photon
//     particle_create_secondary(p, uvw, E, PHOTON, run_CE=true)
//   end do
// }

} // namespace openmc
