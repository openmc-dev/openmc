#include "openmc/physics.h"

#include "openmc/bank.h"
#include "openmc/bremsstrahlung.h"
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
#include "openmc/string_utils.h"
#include "openmc/thermal.h"
#include "openmc/tallies/tally.h"

#include <fmt/core.h>

#include <algorithm> // for max, min, max_element
#include <cmath> // for sqrt, exp, log, abs, copysign

namespace openmc {

//==============================================================================
// Non-member functions
//==============================================================================

void collision(Particle* p)
{
  // Add to collision counter for particle
  ++(p->n_collision_);

  // Sample reaction for the material the particle is in
  switch (static_cast<Particle::Type>(p->type_)) {
  case Particle::Type::neutron:
    sample_neutron_reaction(p);
    break;
  case Particle::Type::photon:
    sample_photon_reaction(p);
    break;
  case Particle::Type::electron:
    sample_electron_reaction(p);
    break;
  case Particle::Type::positron:
    sample_positron_reaction(p);
    break;
  }

  // Kill particle if energy falls below cutoff
  int type = static_cast<int>(p->type_);
  if (p->E_ < settings::energy_cutoff[type]) {
    p->alive_ = false;
    p->wgt_ = 0.0;
  }

  // Display information about collision
  if (settings::verbosity >= 10 || p->trace_) {
    std::string msg;
    if (p->event_ == TallyEvent::KILL) {
      msg = fmt::format("    Killed. Energy = {} eV.", p->E_);
    } else if (p->type_ == Particle::Type::neutron) {
      msg = fmt::format("    {} with {}. Energy = {} eV.",
        reaction_name(p->event_mt_), data::nuclides[p->event_nuclide_]->name_,
        p->E_);
    } else if (p->type_ == Particle::Type::photon) {
      msg = fmt::format("    {} with {}. Energy = {} eV.",
        reaction_name(p->event_mt_),
        to_element(data::nuclides[p->event_nuclide_]->name_), p->E_);
    } else {
      msg = fmt::format("    Disappeared. Energy = {} eV.", p->E_);
    }
    write_message(msg, 1);
  }
}

void sample_neutron_reaction(Particle* p)
{
  // Sample a nuclide within the material
  int i_nuclide = sample_nuclide(p);

  // Save which nuclide particle had collision with
  p->event_nuclide_ = i_nuclide;

  // Create fission bank sites. Note that while a fission reaction is sampled,
  // it never actually "happens", i.e. the weight of the particle does not
  // change when sampling fission sites. The following block handles all
  // absorption (including fission)

  const auto& nuc {data::nuclides[i_nuclide]};

  if (nuc->fissionable_) {
    Reaction* rx = sample_fission(i_nuclide, p);
    if (settings::run_mode == RunMode::EIGENVALUE) {
      create_fission_sites(p, i_nuclide, rx);
    } else if (settings::run_mode == RunMode::FIXED_SOURCE &&
      settings::create_fission_neutrons) {
      create_fission_sites(p, i_nuclide, rx);

      // Make sure particle population doesn't grow out of control for
      // subcritical multiplication problems.
      if (p->secondary_bank_.size() >= 10000) {
        fatal_error("The secondary particle bank appears to be growing without "
        "bound. You are likely running a subcritical multiplication problem "
        "with k-effective close to or greater than one.");
      }
    }
  }

  // Create secondary photons
  if (settings::photon_transport) {
    p->stream_ = STREAM_PHOTON;
    sample_secondary_photons(p, i_nuclide);
	p->stream_ = STREAM_TRACKING;
  }

  // If survival biasing is being used, the following subroutine adjusts the
  // weight of the particle. Otherwise, it checks to see if absorption occurs

  if (p->neutron_xs_[i_nuclide].absorption > 0.0) {
    absorption(p, i_nuclide);
  } else {
    p->wgt_absorb_ = 0.0;
  }
  if (!p->alive_) return;

  // Sample a scattering reaction and determine the secondary energy of the
  // exiting neutron
  scatter(p, i_nuclide);

  // Advance URR seed stream 'N' times after energy changes
  if (p->E_ != p->E_last_) {
    p->stream_ = STREAM_URR_PTABLE;
    advance_prn_seed(data::nuclides.size(), p->current_seed());
    p->stream_ = STREAM_TRACKING;
  }

  // Play russian roulette if survival biasing is turned on
  if (settings::survival_biasing) {
    russian_roulette(p);
    if (!p->alive_) return;
  }
}

void
create_fission_sites(Particle* p, int i_nuclide, const Reaction* rx)
{
  // If uniform fission source weighting is turned on, we increase or decrease
  // the expected number of fission sites produced
  double weight = settings::ufs_on ? ufs_get_weight(p) : 1.0;

  // Determine the expected number of neutrons produced
  double nu_t = p->wgt_ / simulation::keff * weight * p->neutron_xs_[
    i_nuclide].nu_fission / p->neutron_xs_[i_nuclide].total;

  // Sample the number of neutrons produced
  int nu = static_cast<int>(nu_t);
  if (prn(p->current_seed()) <= (nu_t - nu)) ++nu;

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

    // Sample delayed group and angle/energy for fission reaction
    sample_fission_neutron(i_nuclide, rx, p->E_, &site, p->current_seed());

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
    p->delayed_group_ = site.delayed_group;

    // Increment the number of neutrons born delayed
    if (p->delayed_group_ > 0) {
      nu_d[p->delayed_group_-1]++;
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

void sample_photon_reaction(Particle* p)
{
  // Kill photon if below energy cutoff -- an extra check is made here because
  // photons with energy below the cutoff may have been produced by neutrons
  // reactions or atomic relaxation
  int photon = static_cast<int>(Particle::Type::photon);
  if (p->E_ < settings::energy_cutoff[photon]) {
    p->E_ = 0.0;
    p->alive_ = false;
    return;
  }

  // Sample element within material
  int i_element = sample_element(p);
  const auto& micro {p->photon_xs_[i_element]};
  const auto& element {data::elements[i_element]};

  // Calculate photon energy over electron rest mass equivalent
  double alpha = p->E_/MASS_ELECTRON_EV;

  // For tallying purposes, this routine might be called directly. In that
  // case, we need to sample a reaction via the cutoff variable
  double prob = 0.0;
  double cutoff = prn(p->current_seed()) * micro.total;

  // Coherent (Rayleigh) scattering
  prob += micro.coherent;
  if (prob > cutoff) {
    double mu = element.rayleigh_scatter(alpha, p->current_seed());
    p->u() = rotate_angle(p->u(), mu, nullptr, p->current_seed());
    p->event_ = TallyEvent::SCATTER;
    p->event_mt_ = COHERENT;
    return;
  }

  // Incoherent (Compton) scattering
  prob += micro.incoherent;
  if (prob > cutoff) {
    double alpha_out, mu;
    int i_shell;
    element.compton_scatter(alpha, true, &alpha_out, &mu, &i_shell, p->current_seed());

    // Determine binding energy of shell. The binding energy is 0.0 if
    // doppler broadening is not used.
    double e_b;
    if (i_shell == -1) {
      e_b = 0.0;
    } else {
      e_b = element.binding_energy_[i_shell];
    }

    // Create Compton electron
    double phi = 2.0*PI*prn(p->current_seed());
    double E_electron = (alpha - alpha_out)*MASS_ELECTRON_EV - e_b;
    int electron = static_cast<int>(Particle::Type::electron);
    if (E_electron >= settings::energy_cutoff[electron]) {
      double mu_electron = (alpha - alpha_out*mu)
        / std::sqrt(alpha*alpha + alpha_out*alpha_out - 2.0*alpha*alpha_out*mu);
      Direction u = rotate_angle(p->u(), mu_electron, &phi, p->current_seed());
      p->create_secondary(u, E_electron, Particle::Type::electron);
    }

    // TODO: Compton subshell data does not match atomic relaxation data
    // Allow electrons to fill orbital and produce auger electrons
    // and fluorescent photons
    if (i_shell >= 0) {
      const auto& shell = element.shells_[i_shell];
      element.atomic_relaxation(shell, *p);
    }

    phi += PI;
    p->E_ = alpha_out*MASS_ELECTRON_EV;
    p->u() = rotate_angle(p->u(), mu, &phi, p->current_seed());
    p->event_ = TallyEvent::SCATTER;
    p->event_mt_ = INCOHERENT;
    return;
  }

  // Photoelectric effect
  double prob_after = prob + micro.photoelectric;
  if (prob_after > cutoff) {
    for (const auto& shell : element.shells_) {
      // Get grid index and interpolation factor
      int i_grid = micro.index_grid;
      double f = micro.interp_factor;

      // Check threshold of reaction
      int i_start = shell.threshold;
      if (i_grid < i_start) continue;

      // Evaluation subshell photoionization cross section
      double xs = std::exp(shell.cross_section(i_grid - i_start) +
        f*(shell.cross_section(i_grid + 1 - i_start) -
        shell.cross_section(i_grid - i_start)));

      prob += xs;
      if (prob > cutoff) {
        double E_electron = p->E_ - shell.binding_energy;

        // Sample mu using non-relativistic Sauter distribution.
        // See Eqns 3.19 and 3.20 in "Implementing a photon physics
        // model in Serpent 2" by Toni Kaltiaisenaho
        double mu;
        while (true) {
          double r = prn(p->current_seed());
          if (4.0*(1.0 - r)*r >= prn(p->current_seed())) {
            double rel_vel = std::sqrt(E_electron * (E_electron +
              2.0*MASS_ELECTRON_EV)) / (E_electron + MASS_ELECTRON_EV);
            mu = (2.0*r + rel_vel - 1.0) / (2.0*rel_vel*r - rel_vel + 1.0);
            break;
          }
        }

        double phi = 2.0*PI*prn(p->current_seed());
        Direction u;
        u.x = mu;
        u.y = std::sqrt(1.0 - mu*mu)*std::cos(phi);
        u.z = std::sqrt(1.0 - mu*mu)*std::sin(phi);

        // Create secondary electron
        p->create_secondary(u, E_electron, Particle::Type::electron);

        // Allow electrons to fill orbital and produce auger electrons
        // and fluorescent photons
        element.atomic_relaxation(shell, *p);
        p->event_ = TallyEvent::ABSORB;
        p->event_mt_ = 533 + shell.index_subshell;
        p->alive_ = false;
        p->E_ = 0.0;
        return;
      }
    }
  }
  prob = prob_after;

  // Pair production
  prob += micro.pair_production;
  if (prob > cutoff) {
    double E_electron, E_positron;
    double mu_electron, mu_positron;
    element.pair_production(alpha, &E_electron, &E_positron,
      &mu_electron, &mu_positron, p->current_seed());

    // Create secondary electron
    Direction u = rotate_angle(p->u(), mu_electron, nullptr, p->current_seed());
    p->create_secondary(u, E_electron, Particle::Type::electron);

    // Create secondary positron
    u = rotate_angle(p->u(), mu_positron, nullptr, p->current_seed());
    p->create_secondary(u, E_positron, Particle::Type::positron);

    p->event_ = TallyEvent::ABSORB;
    p->event_mt_ = PAIR_PROD;
    p->alive_ = false;
    p->E_ = 0.0;
  }
}

void sample_electron_reaction(Particle* p)
{
  // TODO: create reaction types

  if (settings::electron_treatment == ElectronTreatment::TTB) {
    double E_lost;
    thick_target_bremsstrahlung(*p, &E_lost);
  }

  p->E_ = 0.0;
  p->alive_ = false;
  p->event_ = TallyEvent::ABSORB;
}

void sample_positron_reaction(Particle* p)
{
  // TODO: create reaction types

  if (settings::electron_treatment == ElectronTreatment::TTB) {
    double E_lost;
    thick_target_bremsstrahlung(*p, &E_lost);
  }

  // Sample angle isotropically
  double mu = 2.0*prn(p->current_seed()) - 1.0;
  double phi = 2.0*PI*prn(p->current_seed());
  Direction u;
  u.x = mu;
  u.y = std::sqrt(1.0 - mu*mu)*std::cos(phi);
  u.z = std::sqrt(1.0 - mu*mu)*std::sin(phi);

  // Create annihilation photon pair traveling in opposite directions
  p->create_secondary(u, MASS_ELECTRON_EV, Particle::Type::photon);
  p->create_secondary(-u, MASS_ELECTRON_EV, Particle::Type::photon);

  p->E_ = 0.0;
  p->alive_ = false;
  p->event_ = TallyEvent::ABSORB;
}

int sample_nuclide(Particle* p)
{
  // Sample cumulative distribution function
  double cutoff = prn(p->current_seed()) * p->macro_xs_.total;

  // Get pointers to nuclide/density arrays
  const auto& mat {model::materials[p->material_]};
  int n = mat->nuclide_.size();

  double prob = 0.0;
  for (int i = 0; i < n; ++i) {
    // Get atom density
    int i_nuclide = mat->nuclide_[i];
    double atom_density = mat->atom_density_[i];

    // Increment probability to compare to cutoff
    prob += atom_density * p->neutron_xs_[i_nuclide].total;
    if (prob >= cutoff) return i_nuclide;
  }

  // If we reach here, no nuclide was sampled
  p->write_restart();
  throw std::runtime_error{"Did not sample any nuclide during collision."};
}

int sample_element(Particle* p)
{
  // Sample cumulative distribution function
  double cutoff = prn(p->current_seed()) * p->macro_xs_.total;

  // Get pointers to elements, densities
  const auto& mat {model::materials[p->material_]};

  double prob = 0.0;
  for (int i = 0; i < mat->element_.size(); ++i) {
    // Find atom density
    int i_element = mat->element_[i];
    double atom_density = mat->atom_density_[i];

    // Determine microscopic cross section
    double sigma = atom_density * p->photon_xs_[i_element].total;

    // Increment probability to compare to cutoff
    prob += sigma;
    if (prob > cutoff) {
      // Save which nuclide particle had collision with for tally purpose
      p->event_nuclide_ = mat->nuclide_[i];

      return i_element;
    }
  }

  // If we made it here, no element was sampled
  p->write_restart();
  fatal_error("Did not sample any element during collision.");
}

Reaction* sample_fission(int i_nuclide, Particle* p)
{
  // Get pointer to nuclide
  const auto& nuc {data::nuclides[i_nuclide]};

  // If we're in the URR, by default use the first fission reaction. We also
  // default to the first reaction if we know that there are no partial fission
  // reactions
  if (p->neutron_xs_[i_nuclide].use_ptable || !nuc->has_partial_fission_) {
    return nuc->fission_rx_[0];
  }

  // Check to see if we are in a windowed multipole range.  WMP only supports
  // the first fission reaction.
  if (nuc->multipole_) {
    if (p->E_ >= nuc->multipole_->E_min_ && p->E_ <= nuc->multipole_->E_max_) {
      return nuc->fission_rx_[0];
    }
  }

  // Get grid index and interpolatoin factor and sample fission cdf
  int i_temp = p->neutron_xs_[i_nuclide].index_temp;
  int i_grid = p->neutron_xs_[i_nuclide].index_grid;
  double f = p->neutron_xs_[i_nuclide].interp_factor;
  double cutoff = prn(p->current_seed()) * p->neutron_xs_[i_nuclide].fission;
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

  // If we reached here, no reaction was sampled
  throw std::runtime_error{"No fission reaction was sampled for " + nuc->name_};
}

void sample_photon_product(int i_nuclide, Particle* p, int* i_rx, int* i_product)
{
  // Get grid index and interpolation factor and sample photon production cdf
  int i_temp = p->neutron_xs_[i_nuclide].index_temp;
  int i_grid = p->neutron_xs_[i_nuclide].index_grid;
  double f = p->neutron_xs_[i_nuclide].interp_factor;
  double cutoff = prn(p->current_seed()) * p->neutron_xs_[i_nuclide].photon_prod;
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
      if (rx->products_[j].particle_ == Particle::Type::photon) {
        // add to cumulative probability
        prob += (*rx->products_[j].yield_)(p->E_) * xs;

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
    p->wgt_absorb_ = p->wgt_ * p->neutron_xs_[i_nuclide].absorption /
          p->neutron_xs_[i_nuclide].total;

    // Adjust weight of particle by probability of absorption
    p->wgt_ -= p->wgt_absorb_;
    p->wgt_last_ = p->wgt_;

    // Score implicit absorption estimate of keff
    if (settings::run_mode == RunMode::EIGENVALUE) {
      p->keff_tally_absorption_ += p->wgt_absorb_ * p->neutron_xs_[
        i_nuclide].nu_fission / p->neutron_xs_[i_nuclide].absorption;
    }
  } else {
    // See if disappearance reaction happens
    if (p->neutron_xs_[i_nuclide].absorption >
        prn(p->current_seed()) * p->neutron_xs_[i_nuclide].total) {
      // Score absorption estimate of keff
      if (settings::run_mode == RunMode::EIGENVALUE) {
        p->keff_tally_absorption_ += p->wgt_ * p->neutron_xs_[
          i_nuclide].nu_fission / p->neutron_xs_[i_nuclide].absorption;
      }

      p->alive_ = false;
      p->event_ = TallyEvent::ABSORB;
      p->event_mt_ = N_DISAPPEAR;
    }
  }
}

void scatter(Particle* p, int i_nuclide)
{
  // copy incoming direction
  Direction u_old {p->u()};

  // Get pointer to nuclide and grid index/interpolation factor
  const auto& nuc {data::nuclides[i_nuclide]};
  const auto& micro {p->neutron_xs_[i_nuclide]};
  int i_temp =  micro.index_temp;
  int i_grid =  micro.index_grid;
  double f = micro.interp_factor;

  // For tallying purposes, this routine might be called directly. In that
  // case, we need to sample a reaction via the cutoff variable
  double cutoff = prn(p->current_seed()) * (micro.total - micro.absorption);
  bool sampled = false;

  // Calculate elastic cross section if it wasn't precalculated
  if (micro.elastic == CACHE_INVALID) {
    nuc->calculate_elastic_xs(*p);
  }

  double prob = micro.elastic - micro.thermal;
  if (prob > cutoff) {
    // =======================================================================
    // NON-S(A,B) ELASTIC SCATTERING

    // Determine temperature
    double kT = nuc->multipole_ ? p->sqrtkT_*p->sqrtkT_ : nuc->kTs_[i_temp];

    // Perform collision physics for elastic scattering
    elastic_scatter(i_nuclide, *nuc->reactions_[0], kT, p);

    p->event_mt_ = ELASTIC;
    sampled = true;
  }

  prob = micro.elastic;
  if (prob > cutoff && !sampled) {
    // =======================================================================
    // S(A,B) SCATTERING

    sab_scatter(i_nuclide, micro.index_sab, p);

    p->event_mt_ = ELASTIC;
    sampled = true;
  }

  if (!sampled) {
    // =======================================================================
    // INELASTIC SCATTERING

    int j = 0;
    int i = 0;
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
      if (i_grid < xs.threshold) continue;

      // add to cumulative probability
      prob += (1.0 - f)*xs.value[i_grid - xs.threshold] +
        f*xs.value[i_grid - xs.threshold + 1];
    }

    // Perform collision physics for inelastic scattering
    const auto& rx {nuc->reactions_[i]};
    inelastic_scatter(nuc.get(), rx.get(), p);
    p->event_mt_ = rx->mt_;
  }

  // Set event component
  p->event_ = TallyEvent::SCATTER;

  // Sample new outgoing angle for isotropic-in-lab scattering
  const auto& mat {model::materials[p->material_]};
  if (!mat->p0_.empty()) {
    int i_nuc_mat = mat->mat_nuclide_index_[i_nuclide];
    if (mat->p0_[i_nuc_mat]) {
      // Sample isotropic-in-lab outgoing direction
      double mu = 2.0*prn(p->current_seed()) - 1.0;
      double phi = 2.0*PI*prn(p->current_seed());

      // Change direction of particle
      p->u().x = mu;
      p->u().y = std::sqrt(1.0 - mu*mu)*std::cos(phi);
      p->u().z = std::sqrt(1.0 - mu*mu)*std::sin(phi);
      p->mu_ = u_old.dot(p->u());
    }
  }
}

void elastic_scatter(int i_nuclide, const Reaction& rx, double kT,
  Particle* p)
{
  // get pointer to nuclide
  const auto& nuc {data::nuclides[i_nuclide]};

  double vel = std::sqrt(p->E_);
  double awr = nuc->awr_;

  // Neutron velocity in LAB
  Direction v_n = vel*p->u();

  // Sample velocity of target nucleus
  Direction v_t {};
  if (!p->neutron_xs_[i_nuclide].use_ptable) {
    v_t = sample_target_velocity(nuc.get(), p->E_, p->u(), v_n,
      p->neutron_xs_[i_nuclide].elastic, kT, p->current_seed());
  }

  // Velocity of center-of-mass
  Direction v_cm = (v_n + awr*v_t)/(awr + 1.0);

  // Transform to CM frame
  v_n -= v_cm;

  // Find speed of neutron in CM
  vel = v_n.norm();

  // Sample scattering angle, checking if angle distribution is present (assume
  // isotropic otherwise)
  double mu_cm;
  auto& d = rx.products_[0].distribution_[0];
  auto d_ = dynamic_cast<UncorrelatedAngleEnergy*>(d.get());
  if (!d_->angle().empty()) {
    mu_cm = d_->angle().sample(p->E_, p->current_seed());
  } else {
    mu_cm = 2.0*prn(p->current_seed()) - 1.0;
  }

  // Determine direction cosines in CM
  Direction u_cm = v_n/vel;

  // Rotate neutron velocity vector to new angle -- note that the speed of the
  // neutron in CM does not change in elastic scattering. However, the speed
  // will change when we convert back to LAB
  v_n = vel * rotate_angle(u_cm, mu_cm, nullptr, p->current_seed());

  // Transform back to LAB frame
  v_n += v_cm;

  p->E_ = v_n.dot(v_n);
  vel = std::sqrt(p->E_);

  // compute cosine of scattering angle in LAB frame by taking dot product of
  // neutron's pre- and post-collision angle
  p->mu_ = p->u().dot(v_n) / vel;

  // Set energy and direction of particle in LAB frame
  p->u() = v_n / vel;

  // Because of floating-point roundoff, it may be possible for mu_lab to be
  // outside of the range [-1,1). In these cases, we just set mu_lab to exactly
  // -1 or 1
  if (std::abs(p->mu_) > 1.0) p->mu_ = std::copysign(1.0, p->mu_);
}

void sab_scatter(int i_nuclide, int i_sab, Particle* p)
{
  // Determine temperature index
  const auto& micro {p->neutron_xs_[i_nuclide]};
  int i_temp = micro.index_temp_sab;

  // Sample energy and angle
  double E_out;
  data::thermal_scatt[i_sab]->data_[i_temp].sample(micro, p->E_, &E_out, &p->mu_, p->current_seed());

  // Set energy to outgoing, change direction of particle
  p->E_ = E_out;
  p->u() = rotate_angle(p->u(), p->mu_, nullptr, p->current_seed());
}

Direction sample_target_velocity(const Nuclide* nuc, double E, Direction u,
  Direction v_neut, double xs_eff, double kT, uint64_t* seed)
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
    return sample_cxs_target_velocity(nuc->awr_, E, u, kT, seed);

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
      return sample_cxs_target_velocity(nuc->awr_, E, u, kT, seed);
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
          v_target = sample_cxs_target_velocity(nuc->awr_, E, u, kT, seed);
          Direction v_rel = v_neut - v_target;
          E_rel = v_rel.dot(v_rel);
          if (E_rel < E_up) break;
        }

        // perform Doppler broadening rejection correction (dbrc)
        double xs_0K = nuc->elastic_xs_0K(E_rel);
        double R = xs_0K / xs_max;
        if (prn(seed) < R) return v_target;
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
        double E_t = -kT * std::log(prn(seed));

        // sample a relative energy using the xs cdf
        double cdf_rel = cdf_low + prn(seed)*(cdf_up - cdf_low);
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
          return std::sqrt(E_t) * rotate_angle(u, mu, nullptr, seed);
        }
      }
    }
    } // case RVS, DBRC
  } // switch (sampling_method)

  UNREACHABLE();
}

Direction
sample_cxs_target_velocity(double awr, double E, Direction u, double kT, uint64_t* seed)
{
  double beta_vn = std::sqrt(awr * E / kT);
  double alpha = 1.0/(1.0 + std::sqrt(PI)*beta_vn/2.0);

  double beta_vt_sq;
  double mu;
  while (true) {
    // Sample two random numbers
    double r1 = prn(seed);
    double r2 = prn(seed);

    if (prn(seed) < alpha) {
      // With probability alpha, we sample the distribution p(y) =
      // y*e^(-y). This can be done with sampling scheme C45 frmo the Monte
      // Carlo sampler

      beta_vt_sq = -std::log(r1*r2);

    } else {
      // With probability 1-alpha, we sample the distribution p(y) = y^2 *
      // e^(-y^2). This can be done with sampling scheme C61 from the Monte
      // Carlo sampler

      double c = std::cos(PI/2.0 * prn(seed));
      beta_vt_sq = -std::log(r1) - std::log(r2)*c*c;
    }

    // Determine beta * vt
    double beta_vt = std::sqrt(beta_vt_sq);

    // Sample cosine of angle between neutron and target velocity
    mu = 2.0*prn(seed) - 1.0;

    // Determine rejection probability
    double accept_prob = std::sqrt(beta_vn*beta_vn + beta_vt_sq -
      2*beta_vn*beta_vt*mu) / (beta_vn + beta_vt);

    // Perform rejection sampling on vt and mu
    if (prn(seed) < accept_prob) break;
  }

  // Determine speed of target nucleus
  double vt = std::sqrt(beta_vt_sq*kT/awr);

  // Determine velocity vector of target nucleus based on neutron's velocity
  // and the sampled angle between them
  return vt * rotate_angle(u, mu, nullptr, seed);
}

void sample_fission_neutron(int i_nuclide, const Reaction* rx, double E_in, Particle::Bank* site, uint64_t* seed)
{
  // Sample cosine of angle -- fission neutrons are always emitted
  // isotropically. Sometimes in ACE data, fission reactions actually have
  // an angular distribution listed, but for those that do, it's simply just
  // a uniform distribution in mu
  double mu = 2.0 * prn(seed) - 1.0;

  // Sample azimuthal angle uniformly in [0,2*pi)
  double phi = 2.0*PI*prn(seed);
  site->u.x = mu;
  site->u.y = std::sqrt(1.0 - mu*mu) * std::cos(phi);
  site->u.z = std::sqrt(1.0 - mu*mu) * std::sin(phi);

  // Determine total nu, delayed nu, and delayed neutron fraction
  const auto& nuc {data::nuclides[i_nuclide]};
  double nu_t = nuc->nu(E_in, Nuclide::EmissionMode::total);
  double nu_d = nuc->nu(E_in, Nuclide::EmissionMode::delayed);
  double beta = nu_d / nu_t;

  if (prn(seed) < beta) {
    // ====================================================================
    // DELAYED NEUTRON SAMPLED

    // sampled delayed precursor group
    double xi = prn(seed)*nu_d;
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
      rx->products_[group].sample(E_in, site->E, mu, seed);

      // resample if energy is greater than maximum neutron energy
      constexpr int neutron = static_cast<int>(Particle::Type::neutron);
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
      rx->products_[0].sample(E_in, site->E, mu, seed);

      // resample if energy is greater than maximum neutron energy
      constexpr int neutron = static_cast<int>(Particle::Type::neutron);
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
  double E_in = p->E_;

  // sample outgoing energy and scattering cosine
  double E;
  double mu;
  rx->products_[0].sample(E_in, E, mu, p->current_seed());

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
  p->E_ = E;
  p->mu_ = mu;

  // change direction of particle
  p->u() = rotate_angle(p->u(), mu, nullptr, p->current_seed());

  // evaluate yield
  double yield = (*rx->products_[0].yield_)(E_in);
  if (std::floor(yield) == yield) {
    // If yield is integral, create exactly that many secondary particles
    for (int i = 0; i < static_cast<int>(std::round(yield)) - 1; ++i) {
      p->create_secondary(p->u(), p->E_, Particle::Type::neutron);
    }
  } else {
    // Otherwise, change weight of particle based on yield
    p->wgt_ *= yield;
  }
}

void sample_secondary_photons(Particle* p, int i_nuclide)
{
  // Sample the number of photons produced
  double y_t = p->wgt_ * p->neutron_xs_[i_nuclide].photon_prod /
    p->neutron_xs_[i_nuclide].total;
  int y = static_cast<int>(y_t);
  if (prn(p->current_seed()) <= y_t - y) ++y;

  // Sample each secondary photon
  for (int i = 0; i < y; ++i) {
    // Sample the reaction and product
    int i_rx;
    int i_product;
    sample_photon_product(i_nuclide, p, &i_rx, &i_product);

    // Sample the outgoing energy and angle
    auto& rx = data::nuclides[i_nuclide]->reactions_[i_rx];
    double E;
    double mu;
    rx->products_[i_product].sample(p->E_, E, mu, p->current_seed());

    // Sample the new direction
    Direction u = rotate_angle(p->u(), mu, nullptr, p->current_seed());

    // Create the secondary photon
    p->create_secondary(u, E, Particle::Type::photon);
  }
}

} // namespace openmc
