#include "openmc/tallies/tally_scoring.h"
#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/distribution_multi.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/ifp.h"
#include "openmc/material.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/reaction_product.h"
#include "openmc/search.h"
#include "openmc/secondary_correlated.h"
#include "openmc/secondary_uncorrelated.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/string_utils.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/filter_delayedgroup.h"
#include "openmc/tallies/filter_energy.h"
#include <csignal>
#include <fstream>
#include <iostream>
#include <string>
#include <typeinfo>
int ghost_counter = 0;
int col_counter = 0;
namespace openmc {

//==============================================================================
// FilterBinIter implementation
//==============================================================================

FilterBinIter::FilterBinIter(const Tally& tally, Particle& p)
  : filter_matches_ {p.filter_matches()}, tally_ {tally}
{
  // Find all valid bins in each relevant filter if they have not already been
  // found for this event.
  for (auto i_filt : tally_.filters()) {
    auto& match {filter_matches_[i_filt]};
    if (!match.bins_present_) {
      match.bins_.clear();
      match.weights_.clear();
      model::tally_filters[i_filt]->get_all_bins(p, tally_.estimator_, match);
      match.bins_present_ = true;
    }

    // If there are no valid bins for this filter, then there are no valid
    // filter bin combinations so all iterators are end iterators.
    if (match.bins_.size() == 0) {
      index_ = -1;
      return;
    }

    // Set the index of the bin used in the first filter combination
    match.i_bin_ = 0;
  }

  // Compute the initial index and weight.
  this->compute_index_weight();
}

FilterBinIter::FilterBinIter(
  const Tally& tally, bool end, vector<FilterMatch>* particle_filter_matches)
  : filter_matches_ {*particle_filter_matches}, tally_ {tally}
{
  // Handle the special case for an iterator that points to the end.
  if (end) {
    index_ = -1;
    return;
  }

  for (auto i_filt : tally_.filters()) {
    auto& match {filter_matches_[i_filt]};
    if (!match.bins_present_) {
      match.bins_.clear();
      match.weights_.clear();
      for (auto i = 0; i < model::tally_filters[i_filt]->n_bins(); ++i) {
        match.bins_.push_back(i);
        match.weights_.push_back(1.0);
      }
      match.bins_present_ = true;
    }

    if (match.bins_.size() == 0) {
      index_ = -1;
      return;
    }

    match.i_bin_ = 0;
  }

  // Compute the initial index and weight.
  this->compute_index_weight();
}

FilterBinIter& FilterBinIter::operator++()
{
  // Find the next valid combination of filter bins.  To do this, we search
  // backwards through the filters until we find the first filter whose bins
  // can be incremented.
  bool visited_all_combinations = true;
  for (int i = tally_.filters().size() - 1; i >= 0; --i) {
    auto i_filt = tally_.filters(i);
    auto& match {filter_matches_[i_filt]};
    if (match.i_bin_ < match.bins_.size() - 1) {
      // The bin for this filter can be incremented.  Increment it and do not
      // touch any of the remaining filters.
      ++match.i_bin_;
      visited_all_combinations = false;
      break;
    } else {
      // This bin cannot be incremented so reset it and continue to the next
      // filter.
      match.i_bin_ = 0;
    }
  }

  if (visited_all_combinations) {
    // We have visited every valid combination.  All done!
    index_ = -1;
  } else {
    // The loop found a new valid combination.  Compute the corresponding
    // index and weight.
    compute_index_weight();
  }

  return *this;
}

void FilterBinIter::compute_index_weight()
{
  index_ = 0;
  weight_ = 1.;
  for (auto i = 0; i < tally_.filters().size(); ++i) {
    auto i_filt = tally_.filters(i);
    auto& match {filter_matches_[i_filt]};
    auto i_bin = match.i_bin_;
    index_ += match.bins_[i_bin] * tally_.strides(i);
    weight_ *= match.weights_[i_bin];
  }
}

//==============================================================================
// Non-member functions
//==============================================================================

//! Helper function used to increment tallies with a delayed group filter.

void score_fission_delayed_dg(int i_tally, int d_bin, double score,
  int score_index, vector<FilterMatch>& filter_matches)
{
  // Save the original delayed group bin
  auto& tally {*model::tallies[i_tally]};
  auto i_filt = tally.filters(tally.delayedgroup_filter_);
  auto& dg_match {filter_matches[i_filt]};
  auto i_bin = dg_match.i_bin_;
  auto original_bin = dg_match.bins_[i_bin];
  dg_match.bins_[i_bin] = d_bin;

  // Determine the filter scoring index
  auto filter_index = 0;
  double filter_weight = 1.;
  for (auto i = 0; i < tally.filters().size(); ++i) {
    auto i_filt = tally.filters(i);
    auto& match {filter_matches[i_filt]};
    auto i_bin = match.i_bin_;
    filter_index += match.bins_[i_bin] * tally.strides(i);
    filter_weight *= match.weights_[i_bin];
  }

// Update the tally result
#pragma omp atomic
  tally.results_(filter_index, score_index, TallyResult::VALUE) +=
    score * filter_weight;

  // Reset the original delayed group bin
  dg_match.bins_[i_bin] = original_bin;
}

//! Helper function to retrieve fission q value from a nuclide

double get_nuc_fission_q(const Nuclide& nuc, const Particle& p, int score_bin)
{
  if (score_bin == SCORE_FISS_Q_PROMPT) {
    if (nuc.fission_q_prompt_) {
      return (*nuc.fission_q_prompt_)(p.E_last());
    }
  } else if (score_bin == SCORE_FISS_Q_RECOV) {
    if (nuc.fission_q_recov_) {
      return (*nuc.fission_q_recov_)(p.E_last());
    }
  }
  return 0.0;
}

//! Helper function to score fission energy
//
//! Pulled out to support both the fission_q scores and energy deposition
//! score

double score_fission_q(const Particle& p, int score_bin, const Tally& tally,
  double flux, int i_nuclide, double atom_density)
{
  if (tally.estimator_ == TallyEstimator::ANALOG) {
    const Nuclide& nuc {*data::nuclides[p.event_nuclide()]};
    if (settings::survival_biasing) {
      // No fission events occur if survival biasing is on -- need to
      // calculate fraction of absorptions that would have resulted in
      // fission scaled by the Q-value
      if (p.neutron_xs(p.event_nuclide()).total > 0) {
        return p.wgt_last() * get_nuc_fission_q(nuc, p, score_bin) *
               p.neutron_xs(p.event_nuclide()).fission * flux /
               p.neutron_xs(p.event_nuclide()).total;
      }
    } else {
      // Skip any non-absorption events
      if (p.event() == TallyEvent::SCATTER)
        return 0.0;
      // All fission events will contribute, so again we can use particle's
      // weight entering the collision as the estimate for the fission
      // reaction rate
      if (p.neutron_xs(p.event_nuclide()).absorption > 0) {
        return p.wgt_last() * get_nuc_fission_q(nuc, p, score_bin) *
               p.neutron_xs(p.event_nuclide()).fission * flux /
               p.neutron_xs(p.event_nuclide()).absorption;
      }
    }
  } else {
    if (i_nuclide >= 0) {
      const Nuclide& nuc {*data::nuclides[i_nuclide]};
      return get_nuc_fission_q(nuc, p, score_bin) * atom_density * flux *
             p.neutron_xs(i_nuclide).fission;
    } else if (p.material() != MATERIAL_VOID) {
      const Material& material {*model::materials[p.material()]};
      double score {0.0};
      for (auto i = 0; i < material.nuclide_.size(); ++i) {
        auto j_nuclide = material.nuclide_[i];
        auto atom_density = material.atom_density_(i);
        const Nuclide& nuc {*data::nuclides[j_nuclide]};
        score += get_nuc_fission_q(nuc, p, score_bin) * atom_density *
                 p.neutron_xs(j_nuclide).fission;
      }
      return score * flux;
    }
  }
  return 0.0;
}

//! Helper function to obtain the kerma coefficient for a given nuclide

double get_nuclide_neutron_heating(
  const Particle& p, const Nuclide& nuc, int rxn_index, int i_nuclide)
{
  size_t mt = nuc.reaction_index_[rxn_index];
  if (mt == C_NONE)
    return 0.0;

  const auto& micro = p.neutron_xs(i_nuclide);
  auto i_temp = micro.index_temp;
  if (i_temp < 0)
    return 0.0; // Can be true due to multipole

  // Determine total kerma
  const auto& rx {*nuc.reactions_[mt]};
  double kerma = rx.xs(micro);
  if (kerma == 0.0)
    return 0.0;

  if (settings::run_mode == RunMode::EIGENVALUE) {
    // Determine kerma for fission as (EFR + EB)*sigma_f
    double kerma_fission =
      nuc.fragments_
        ? ((*nuc.fragments_)(p.E_last()) + (*nuc.betas_)(p.E_last())) *
            p.neutron_xs(i_nuclide).fission
        : 0.0;

    // Determine non-fission kerma as difference
    double kerma_non_fission = kerma - kerma_fission;

    // Re-weight non-fission kerma by keff to properly balance energy release
    // and deposition. See D. P. Griesheimer, S. J. Douglass, and M. H. Stedry,
    // "Self-consistent energy normalization for quasistatic reactor
    // calculations", Proc. PHYSOR, Cambridge, UK, Mar 29-Apr 2, 2020.
    kerma = simulation::keff * kerma_non_fission + kerma_fission;
  }
  return kerma;
}

//! Helper function to obtain neutron heating [eV]

double score_neutron_heating(const Particle& p, const Tally& tally, double flux,
  int rxn_bin, int i_nuclide, double atom_density)
{
  // Get heating macroscopic "cross section"
  double heating_xs;
  if (i_nuclide >= 0) {
    const Nuclide& nuc {*data::nuclides[i_nuclide]};
    heating_xs = get_nuclide_neutron_heating(p, nuc, rxn_bin, i_nuclide);
    if (tally.estimator_ == TallyEstimator::ANALOG) {
      heating_xs /= p.neutron_xs(i_nuclide).total;
    } else {
      heating_xs *= atom_density;
    }
  } else {
    if (p.material() != MATERIAL_VOID) {
      heating_xs = 0.0;
      const Material& material {*model::materials[p.material()]};
      for (auto i = 0; i < material.nuclide_.size(); ++i) {
        int j_nuclide = material.nuclide_[i];
        double atom_density {material.atom_density_(i)};
        const Nuclide& nuc {*data::nuclides[j_nuclide]};
        heating_xs += atom_density *
                      get_nuclide_neutron_heating(p, nuc, rxn_bin, j_nuclide);
      }
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        heating_xs /= p.macro_xs().total;
      }
    }
  }
  double score = heating_xs * flux;
  if (tally.estimator_ == TallyEstimator::ANALOG) {
    // All events score to a heating tally bin. We actually use a
    // collision estimator in place of an analog one since there is no
    // reaction-wise heating cross section
    score *= p.wgt_last();
  }
  return score;
}

//! Helper function to obtain reaction Q value for photons and charged particles
double get_reaction_q_value(const Particle& p)
{
  if (p.type() == ParticleType::photon && p.event_mt() == PAIR_PROD) {
    // pair production
    return -2 * MASS_ELECTRON_EV;
  } else if (p.type() == ParticleType::positron) {
    // positron annihilation
    return 2 * MASS_ELECTRON_EV;
  } else {
    return 0.0;
  }
}

//! Helper function to obtain particle heating [eV]

double score_particle_heating(const Particle& p, const Tally& tally,
  double flux, int rxn_bin, int i_nuclide, double atom_density)
{
  if (p.type() == ParticleType::neutron)
    return score_neutron_heating(
      p, tally, flux, rxn_bin, i_nuclide, atom_density);
  if (i_nuclide == -1 || i_nuclide == p.event_nuclide() ||
      p.event_nuclide() == -1) {
    // For pair production and positron annihilation, we need to account for the
    // reaction Q value
    double Q = get_reaction_q_value(p);

    // Get the pre-collision energy of the particle.
    auto E = p.E_last();

    // The energy deposited is the sum of the incident energy and the reaction
    // Q-value less the energy of any outgoing particles
    double score = E + Q - p.E() - p.bank_second_E();

    score *= p.wgt_last();

    // if no event_nuclide (charged particle) scale energy deposition by
    // fractional charge density
    if (i_nuclide != -1 && p.event_nuclide() == -1) {
      const auto& mat {model::materials[p.material()]};
      int z = data::nuclides[i_nuclide]->Z_;
      auto i = mat->mat_nuclide_index_[i_nuclide];
      score *= (z * mat->atom_density_[i] / mat->charge_density());
    }
    return score;
  }
  return 0.0;
}

//! Helper function for nu-fission tallies with energyout filters.
//
//! In this case, we may need to score to multiple bins if there were multiple
//! neutrons produced with different energies.

void score_fission_eout(Particle& p, int i_tally, int i_score, int score_bin)
{
  auto& tally {*model::tallies[i_tally]};
  auto i_eout_filt = tally.filters()[tally.energyout_filter_];
  auto i_bin = p.filter_matches(i_eout_filt).i_bin_;
  auto bin_energyout = p.filter_matches(i_eout_filt).bins_[i_bin];

  const EnergyoutFilter& eo_filt {
    *dynamic_cast<EnergyoutFilter*>(model::tally_filters[i_eout_filt].get())};

  // Note that the score below is weighted by keff. Since the creation of
  // fission sites is weighted such that it is expected to create n_particles
  // sites, we need to multiply the score by keff to get the true nu-fission
  // rate. Otherwise, the sum of all nu-fission rates would be ~1.0.

  // loop over number of particles banked
  for (auto i = 0; i < p.n_bank(); ++i) {
    const auto& bank = p.nu_bank(i);

    // get the delayed group
    auto g = bank.delayed_group;

    // determine score based on bank site weight and keff
    double score = simulation::keff * bank.wgt;

    // Add derivative information for differential tallies.  Note that the
    // i_nuclide and atom_density arguments do not matter since this is an
    // analog estimator.
    if (tally.deriv_ != C_NONE)
      apply_derivative_to_score(p, i_tally, 0, 0., SCORE_NU_FISSION, score);

    if (!settings::run_CE && eo_filt.matches_transport_groups()) {

      // determine outgoing energy group from fission bank
      auto g_out = static_cast<int>(bank.E);

      // modify the value so that g_out = 0 corresponds to the highest energy
      // bin
      g_out = eo_filt.n_bins() - g_out - 1;

      // change outgoing energy bin
      p.filter_matches(i_eout_filt).bins_[i_bin] = g_out;

    } else {

      double E_out;
      if (settings::run_CE) {
        E_out = bank.E;
      } else {
        E_out = data::mg.energy_bin_avg_[static_cast<int>(bank.E)];
      }

      // Set EnergyoutFilter bin index
      if (E_out < eo_filt.bins().front() || E_out > eo_filt.bins().back()) {
        continue;
      } else {
        auto i_match = lower_bound_index(
          eo_filt.bins().begin(), eo_filt.bins().end(), E_out);
        p.filter_matches(i_eout_filt).bins_[i_bin] = i_match;
      }
    }

    // Case for tallying prompt neutrons
    if (score_bin == SCORE_NU_FISSION ||
        (score_bin == SCORE_PROMPT_NU_FISSION && g == 0)) {

      // Find the filter scoring index for this filter combination
      int filter_index = 0;
      double filter_weight = 1.0;
      for (auto j = 0; j < tally.filters().size(); ++j) {
        auto i_filt = tally.filters(j);
        auto& match {p.filter_matches(i_filt)};
        auto i_bin = match.i_bin_;
        filter_index += match.bins_[i_bin] * tally.strides(j);
        filter_weight *= match.weights_[i_bin];
      }

// Update tally results
#pragma omp atomic
      tally.results_(filter_index, i_score, TallyResult::VALUE) +=
        score * filter_weight;

    } else if (score_bin == SCORE_DELAYED_NU_FISSION && g != 0) {

      // If the delayed group filter is present, tally to corresponding delayed
      // group bin if it exists
      if (tally.delayedgroup_filter_ >= 0) {

        // Get the index of the delayed group filter
        auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];

        const DelayedGroupFilter& dg_filt {*dynamic_cast<DelayedGroupFilter*>(
          model::tally_filters[i_dg_filt].get())};

        // Loop over delayed group bins until the corresponding bin is found
        for (auto d_bin = 0; d_bin < dg_filt.n_bins(); ++d_bin) {
          if (dg_filt.groups()[d_bin] == g) {
            // Find the filter index and weight for this filter combination
            double filter_weight = 1.;
            for (auto j = 0; j < tally.filters().size(); ++j) {
              auto i_filt = tally.filters(j);
              auto& match {p.filter_matches(i_filt)};
              auto i_bin = match.i_bin_;
              filter_weight *= match.weights_[i_bin];
            }

            score_fission_delayed_dg(i_tally, d_bin, score * filter_weight,
              i_score, p.filter_matches());
          }
        }

        // If the delayed group filter is not present, add score to tally
      } else {

        // Find the filter index and weight for this filter combination
        int filter_index = 0;
        double filter_weight = 1.;
        for (auto j = 0; j < tally.filters().size(); ++j) {
          auto i_filt = tally.filters(j);
          auto& match {p.filter_matches(i_filt)};
          auto i_bin = match.i_bin_;
          filter_index += match.bins_[i_bin] * tally.strides(j);
          filter_weight *= match.weights_[i_bin];
        }

// Update tally results
#pragma omp atomic
        tally.results_(filter_index, i_score, TallyResult::VALUE) +=
          score * filter_weight;
      }
    }
  }

  // Reset outgoing energy bin and score index
  p.filter_matches(i_eout_filt).bins_[i_bin] = bin_energyout;
}

double get_nuclide_xs(const Particle& p, int i_nuclide, int score_bin)
{
  const auto& nuc {*data::nuclides[i_nuclide]};

  // Get reaction object, or return 0 if reaction is not present
  auto m = nuc.reaction_index_[score_bin];
  if (m == C_NONE)
    return 0.0;
  const auto& rx {*nuc.reactions_[m]};
  const auto& micro {p.neutron_xs(i_nuclide)};

  // In the URR, the (n,gamma) cross section is sampled randomly from
  // probability tables. Make sure we use the sampled value (which is equal to
  // absorption - fission) rather than the dilute average value
  if (micro.use_ptable && score_bin == N_GAMMA) {
    return micro.absorption - micro.fission;
  }

  auto i_temp = micro.index_temp;
  if (i_temp >= 0) { // Can be false due to multipole
    // Get index on energy grid and interpolation factor
    auto i_grid = micro.index_grid;
    auto f = micro.interp_factor;

    // Calculate interpolated cross section
    double xs = rx.xs(micro);

    if (settings::run_mode == RunMode::EIGENVALUE &&
        score_bin == HEATING_LOCAL) {
      // Determine kerma for fission as (EFR + EGP + EGD + EB)*sigma_f
      double kerma_fission =
        nuc.fragments_
          ? ((*nuc.fragments_)(p.E_last()) + (*nuc.betas_)(p.E_last()) +
              (*nuc.prompt_photons_)(p.E_last()) +
              (*nuc.delayed_photons_)(p.E_last())) *
              micro.fission
          : 0.0;

      // Determine non-fission kerma as difference
      double kerma_non_fission = xs - kerma_fission;

      // Re-weight non-fission kerma by keff to properly balance energy release
      // and deposition. See D. P. Griesheimer, S. J. Douglass, and M. H.
      // Stedry, "Self-consistent energy normalization for quasistatic reactor
      // calculations", Proc. PHYSOR, Cambridge, UK, Mar 29-Apr 2, 2020.
      xs = simulation::keff * kerma_non_fission + kerma_fission;
    }
    return xs;
  } else {
    // For multipole, calculate (n,gamma) from other reactions
    return rx.mt_ == N_GAMMA ? micro.absorption - micro.fission : 0.0;
  }
  return 0.0;
}

//! Update tally results for continuous-energy tallies with a tracklength or
//! collision estimator.

void score_general_ce_nonanalog(Particle& p, int i_tally, int start_index,
  int filter_index, double filter_weight, int i_nuclide, double atom_density,
  double flux)
{
  Tally& tally {*model::tallies[i_tally]};

  // Get the pre-collision energy of the particle.
  auto E = p.E_last();

  using Type = ParticleType;

  for (auto i = 0; i < tally.scores_.size(); ++i) {
    auto score_bin = tally.scores_[i];
    auto score_index = start_index + i;
    double score = 0.0;

    switch (score_bin) {
    case SCORE_FLUX:
      score = flux;
      break;

    case SCORE_TOTAL:
      if (i_nuclide >= 0) {
        if (p.type() == Type::neutron) {
          score = p.neutron_xs(i_nuclide).total * atom_density * flux;
        } else if (p.type() == Type::photon) {
          score = p.photon_xs(i_nuclide).total * atom_density * flux;
        }
      } else {
        score = p.macro_xs().total * flux;
      }
      break;

    case SCORE_INVERSE_VELOCITY:
      if (p.type() != Type::neutron)
        continue;

      // Score inverse velocity in units of s/cm.
      score = flux / p.speed();
      break;

    case SCORE_SCATTER:
      if (p.type() != Type::neutron && p.type() != Type::photon)
        continue;

      if (i_nuclide >= 0) {
        if (p.type() == Type::neutron) {
          const auto& micro = p.neutron_xs(i_nuclide);
          score = (micro.total - micro.absorption) * atom_density * flux;
        } else {
          const auto& micro = p.photon_xs(i_nuclide);
          score = (micro.coherent + micro.incoherent) * atom_density * flux;
        }
      } else {
        if (p.type() == Type::neutron) {
          score = (p.macro_xs().total - p.macro_xs().absorption) * flux;
        } else {
          score = (p.macro_xs().coherent + p.macro_xs().incoherent) * flux;
        }
      }
      break;

    case SCORE_ABSORPTION:
      if (p.type() != Type::neutron && p.type() != Type::photon)
        continue;

      if (i_nuclide >= 0) {
        if (p.type() == Type::neutron) {
          score = p.neutron_xs(i_nuclide).absorption * atom_density * flux;
        } else {
          const auto& xs = p.photon_xs(i_nuclide);
          score =
            (xs.total - xs.coherent - xs.incoherent) * atom_density * flux;
        }
      } else {
        if (p.type() == Type::neutron) {
          score = p.macro_xs().absorption * flux;
        } else {
          score =
            (p.macro_xs().photoelectric + p.macro_xs().pair_production) * flux;
        }
      }
      break;

    case SCORE_FISSION:
      if (p.macro_xs().fission == 0)
        continue;

      if (i_nuclide >= 0) {
        score = p.neutron_xs(i_nuclide).fission * atom_density * flux;
      } else {
        score = p.macro_xs().fission * flux;
      }
      break;

    case SCORE_NU_FISSION:
      if (p.macro_xs().fission == 0)
        continue;

      if (i_nuclide >= 0) {
        score = p.neutron_xs(i_nuclide).nu_fission * atom_density * flux;
      } else {
        score = p.macro_xs().nu_fission * flux;
      }
      break;

    case SCORE_PROMPT_NU_FISSION:
      if (p.macro_xs().fission == 0)
        continue;
      if (i_nuclide >= 0) {
        score = p.neutron_xs(i_nuclide).fission *
                data::nuclides[i_nuclide]->nu(
                  E, ReactionProduct::EmissionMode::prompt) *
                atom_density * flux;
      } else {
        score = 0.;
        // Add up contributions from each nuclide in the material.
        if (p.material() != MATERIAL_VOID) {
          const Material& material {*model::materials[p.material()]};
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto j_nuclide = material.nuclide_[i];
            auto atom_density = material.atom_density_(i);
            score += p.neutron_xs(j_nuclide).fission *
                     data::nuclides[j_nuclide]->nu(
                       E, ReactionProduct::EmissionMode::prompt) *
                     atom_density * flux;
          }
        }
      }
      break;

    case SCORE_DELAYED_NU_FISSION:
      if (p.macro_xs().fission == 0)
        continue;
      if (i_nuclide >= 0) {
        if (tally.delayedgroup_filter_ != C_NONE) {
          auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
          const DelayedGroupFilter& filt {*dynamic_cast<DelayedGroupFilter*>(
            model::tally_filters[i_dg_filt].get())};
          // Tally each delayed group bin individually
          for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
            auto d = filt.groups()[d_bin];
            auto yield = data::nuclides[i_nuclide]->nu(
              E, ReactionProduct::EmissionMode::delayed, d);
            score =
              p.neutron_xs(i_nuclide).fission * yield * atom_density * flux;
            score_fission_delayed_dg(
              i_tally, d_bin, score, score_index, p.filter_matches());
          }
          continue;
        } else {
          // If the delayed group filter is not present, compute the score
          // by multiplying the delayed-nu-fission macro xs by the flux
          score = p.neutron_xs(i_nuclide).fission *
                  data::nuclides[i_nuclide]->nu(
                    E, ReactionProduct::EmissionMode::delayed) *
                  atom_density * flux;
        }
      } else {
        // Need to add up contributions for each nuclide
        if (tally.delayedgroup_filter_ != C_NONE) {
          auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
          const DelayedGroupFilter& filt {*dynamic_cast<DelayedGroupFilter*>(
            model::tally_filters[i_dg_filt].get())};
          if (p.material() != MATERIAL_VOID) {
            const Material& material {*model::materials[p.material()]};
            for (auto i = 0; i < material.nuclide_.size(); ++i) {
              auto j_nuclide = material.nuclide_[i];
              auto atom_density = material.atom_density_(i);
              // Tally each delayed group bin individually
              for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                auto d = filt.groups()[d_bin];
                auto yield = data::nuclides[j_nuclide]->nu(
                  E, ReactionProduct::EmissionMode::delayed, d);
                score =
                  p.neutron_xs(j_nuclide).fission * yield * atom_density * flux;
                score_fission_delayed_dg(
                  i_tally, d_bin, score, score_index, p.filter_matches());
              }
            }
          }
          continue;
        } else {
          score = 0.;
          if (p.material() != MATERIAL_VOID) {
            const Material& material {*model::materials[p.material()]};
            for (auto i = 0; i < material.nuclide_.size(); ++i) {
              auto j_nuclide = material.nuclide_[i];
              auto atom_density = material.atom_density_(i);
              score += p.neutron_xs(j_nuclide).fission *
                       data::nuclides[j_nuclide]->nu(
                         E, ReactionProduct::EmissionMode::delayed) *
                       atom_density * flux;
            }
          }
        }
      }
      break;

    case SCORE_DECAY_RATE:
      if (p.macro_xs().fission == 0)
        continue;
      if (i_nuclide >= 0) {
        const auto& nuc {*data::nuclides[i_nuclide]};
        if (!nuc.fissionable_)
          continue;
        const auto& rxn {*nuc.fission_rx_[0]};
        if (tally.delayedgroup_filter_ != C_NONE) {
          auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
          const DelayedGroupFilter& filt {*dynamic_cast<DelayedGroupFilter*>(
            model::tally_filters[i_dg_filt].get())};
          // Tally each delayed group bin individually
          for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
            auto d = filt.groups()[d_bin];
            auto yield = nuc.nu(E, ReactionProduct::EmissionMode::delayed, d);
            auto rate = rxn.products_[d].decay_rate_;
            score = p.neutron_xs(i_nuclide).fission * yield * flux *
                    atom_density * rate;
            score_fission_delayed_dg(
              i_tally, d_bin, score, score_index, p.filter_matches());
          }
          continue;
        } else {
          score = 0.;
          // We need to be careful not to overshoot the number of
          // delayed groups since this could cause the range of the
          // rxn.products_ array to be exceeded. Hence, we use the size
          // of this array and not the MAX_DELAYED_GROUPS constant for
          // this loop.
          for (auto d = 1; d < rxn.products_.size(); ++d) {
            const auto& product = rxn.products_[d];
            if (product.particle_ != Type::neutron)
              continue;

            auto yield = nuc.nu(E, ReactionProduct::EmissionMode::delayed, d);
            auto rate = product.decay_rate_;
            score += p.neutron_xs(i_nuclide).fission * flux * yield *
                     atom_density * rate;
          }
        }
      } else {
        if (tally.delayedgroup_filter_ != C_NONE) {
          auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
          const DelayedGroupFilter& filt {*dynamic_cast<DelayedGroupFilter*>(
            model::tally_filters[i_dg_filt].get())};
          if (p.material() != MATERIAL_VOID) {
            const Material& material {*model::materials[p.material()]};
            for (auto i = 0; i < material.nuclide_.size(); ++i) {
              auto j_nuclide = material.nuclide_[i];
              auto atom_density = material.atom_density_(i);
              const auto& nuc {*data::nuclides[j_nuclide]};
              if (nuc.fissionable_) {
                const auto& rxn {*nuc.fission_rx_[0]};
                // Tally each delayed group bin individually
                for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                  auto d = filt.groups()[d_bin];
                  auto yield =
                    nuc.nu(E, ReactionProduct::EmissionMode::delayed, d);
                  auto rate = rxn.products_[d].decay_rate_;
                  score = p.neutron_xs(j_nuclide).fission * yield * flux *
                          atom_density * rate;
                  score_fission_delayed_dg(
                    i_tally, d_bin, score, score_index, p.filter_matches());
                }
              }
            }
          }
          continue;
        } else {
          score = 0.;
          if (p.material() != MATERIAL_VOID) {
            const Material& material {*model::materials[p.material()]};
            for (auto i = 0; i < material.nuclide_.size(); ++i) {
              auto j_nuclide = material.nuclide_[i];
              auto atom_density = material.atom_density_(i);
              const auto& nuc {*data::nuclides[j_nuclide]};
              if (nuc.fissionable_) {
                const auto& rxn {*nuc.fission_rx_[0]};
                // We need to be careful not to overshoot the number of
                // delayed groups since this could cause the range of the
                // rxn.products_ array to be exceeded. Hence, we use the size
                // of this array and not the MAX_DELAYED_GROUPS constant for
                // this loop.
                for (auto d = 1; d < rxn.products_.size(); ++d) {
                  const auto& product = rxn.products_[d];
                  if (product.particle_ != Type::neutron)
                    continue;

                  auto yield =
                    nuc.nu(E, ReactionProduct::EmissionMode::delayed, d);
                  auto rate = product.decay_rate_;
                  score += p.neutron_xs(j_nuclide).fission * yield *
                           atom_density * flux * rate;
                }
              }
            }
          }
        }
      }
      break;

    case SCORE_KAPPA_FISSION:
      if (p.macro_xs().fission == 0.)
        continue;
      score = 0.;
      // Kappa-fission values are determined from the Q-value listed for the
      // fission cross section.
      if (i_nuclide >= 0) {
        const auto& nuc {*data::nuclides[i_nuclide]};
        if (nuc.fissionable_) {
          const auto& rxn {*nuc.fission_rx_[0]};
          score = rxn.q_value_ * p.neutron_xs(i_nuclide).fission *
                  atom_density * flux;
        }
      } else if (p.material() != MATERIAL_VOID) {
        const Material& material {*model::materials[p.material()]};
        for (auto i = 0; i < material.nuclide_.size(); ++i) {
          auto j_nuclide = material.nuclide_[i];
          auto atom_density = material.atom_density_(i);
          const auto& nuc {*data::nuclides[j_nuclide]};
          if (nuc.fissionable_) {
            const auto& rxn {*nuc.fission_rx_[0]};
            score += rxn.q_value_ * p.neutron_xs(j_nuclide).fission *
                     atom_density * flux;
          }
        }
      }
      break;

    case SCORE_EVENTS:
// Simply count the number of scoring events
#pragma omp atomic
      tally.results_(filter_index, score_index, TallyResult::VALUE) += 1.0;
      continue;

    case ELASTIC:
      if (p.type() != Type::neutron)
        continue;

      if (i_nuclide >= 0) {
        if (p.neutron_xs(i_nuclide).elastic == CACHE_INVALID)
          data::nuclides[i_nuclide]->calculate_elastic_xs(p);
        score = p.neutron_xs(i_nuclide).elastic * atom_density * flux;
      } else {
        score = 0.;
        if (p.material() != MATERIAL_VOID) {
          const Material& material {*model::materials[p.material()]};
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto j_nuclide = material.nuclide_[i];
            auto atom_density = material.atom_density_(i);
            if (p.neutron_xs(j_nuclide).elastic == CACHE_INVALID)
              data::nuclides[j_nuclide]->calculate_elastic_xs(p);
            score += p.neutron_xs(j_nuclide).elastic * atom_density * flux;
          }
        }
      }
      break;

    case SCORE_FISS_Q_PROMPT:
    case SCORE_FISS_Q_RECOV:
      if (p.macro_xs().fission == 0.)
        continue;
      score =
        score_fission_q(p, score_bin, tally, flux, i_nuclide, atom_density);
      break;

    case SCORE_IFP_TIME_NUM:
      if (settings::ifp_on) {
        if ((p.type() == Type::neutron) && (p.fission())) {
          if (is_generation_time_or_both()) {
            const auto& lifetimes =
              simulation::ifp_source_lifetime_bank[p.current_work() - 1];
            if (lifetimes.size() == settings::ifp_n_generation) {
              score = lifetimes[0] * p.wgt_last();
            }
          }
        }
      }
      break;

    case SCORE_IFP_BETA_NUM:
      if (settings::ifp_on) {
        if ((p.type() == Type::neutron) && (p.fission())) {
          if (is_beta_effective_or_both()) {
            const auto& delayed_groups =
              simulation::ifp_source_delayed_group_bank[p.current_work() - 1];
            if (delayed_groups.size() == settings::ifp_n_generation) {
              if (delayed_groups[0] > 0) {
                score = p.wgt_last();
              }
            }
          }
        }
      }
      break;

    case SCORE_IFP_DENOM:
      if (settings::ifp_on) {
        if ((p.type() == Type::neutron) && (p.fission())) {
          int ifp_data_size;
          if (is_beta_effective_or_both()) {
            ifp_data_size = static_cast<int>(
              simulation::ifp_source_delayed_group_bank[p.current_work() - 1]
                .size());
          } else {
            ifp_data_size = static_cast<int>(
              simulation::ifp_source_lifetime_bank[p.current_work() - 1]
                .size());
          }
          if (ifp_data_size == settings::ifp_n_generation) {
            score = p.wgt_last();
          }
        }
      }
      break;

    case N_2N:
    case N_3N:
    case N_4N:
    case N_GAMMA:
    case N_P:
    case N_A:
      // This case block only works if cross sections for these reactions have
      // been precalculated. When they are not, we revert to the default case,
      // which looks up cross sections
      if (!simulation::need_depletion_rx)
        goto default_case;

      if (p.type() != Type::neutron)
        continue;

      int m;
      switch (score_bin) {
        // clang-format off
      case N_GAMMA: m = 0; break;
      case N_P:     m = 1; break;
      case N_A:     m = 2; break;
      case N_2N:    m = 3; break;
      case N_3N:    m = 4; break;
      case N_4N:    m = 5; break;
        // clang-format on
      }
      if (i_nuclide >= 0) {
        score = p.neutron_xs(i_nuclide).reaction[m] * atom_density * flux;
      } else {
        score = 0.;
        if (p.material() != MATERIAL_VOID) {
          const Material& material {*model::materials[p.material()]};
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto j_nuclide = material.nuclide_[i];
            auto atom_density = material.atom_density_(i);
            score += p.neutron_xs(j_nuclide).reaction[m] * atom_density * flux;
          }
        }
      }
      break;

    case COHERENT:
    case INCOHERENT:
    case PHOTOELECTRIC:
    case PAIR_PROD:
      if (p.type() != Type::photon)
        continue;

      if (i_nuclide >= 0) {
        const auto& micro = p.photon_xs(i_nuclide);
        double xs = (score_bin == COHERENT)        ? micro.coherent
                    : (score_bin == INCOHERENT)    ? micro.incoherent
                    : (score_bin == PHOTOELECTRIC) ? micro.photoelectric
                                                   : micro.pair_production;
        score = xs * atom_density * flux;
      } else {
        double xs = (score_bin == COHERENT)     ? p.macro_xs().coherent
                    : (score_bin == INCOHERENT) ? p.macro_xs().incoherent
                    : (score_bin == PHOTOELECTRIC)
                      ? p.macro_xs().photoelectric
                      : p.macro_xs().pair_production;
        score = xs * flux;
      }
      break;

    case HEATING:
      score = score_particle_heating(
        p, tally, flux, HEATING, i_nuclide, atom_density);
      break;

    default:
    default_case:

      // The default block is really only meant for redundant neutron reactions
      // (e.g. 444, 901)
      if (p.type() != Type::neutron)
        continue;

      // Any other cross section has to be calculated on-the-fly
      if (score_bin < 2)
        fatal_error("Invalid score type on tally " + std::to_string(tally.id_));
      score = 0.;
      if (i_nuclide >= 0) {
        score = get_nuclide_xs(p, i_nuclide, score_bin) * atom_density * flux;
      } else if (p.material() != MATERIAL_VOID) {
        const Material& material {*model::materials[p.material()]};
        for (auto i = 0; i < material.nuclide_.size(); ++i) {
          auto j_nuclide = material.nuclide_[i];
          auto atom_density = material.atom_density_(i);
          score +=
            get_nuclide_xs(p, j_nuclide, score_bin) * atom_density * flux;
        }
      }
    }

    // Add derivative information on score for differential tallies.
    if (tally.deriv_ != C_NONE)
      apply_derivative_to_score(
        p, i_tally, i_nuclide, atom_density, score_bin, score);

// Update tally results
#pragma omp atomic
    tally.results_(filter_index, score_index, TallyResult::VALUE) +=
      score * filter_weight;
  }
}

//! Update tally results for continuous-energy tallies with an analog estimator.
//
//! For analog tallies, the flux estimate depends on the score type so the flux
//! argument is really just used for filter weights. The atom_density argument
//! is not used for analog tallies.

void score_general_ce_analog(Particle& p, int i_tally, int start_index,
  int filter_index, double filter_weight, int i_nuclide, double atom_density,
  double flux)
{
  Tally& tally {*model::tallies[i_tally]};

  // Get the pre-collision energy of the particle.
  auto E = p.E_last();

  // Determine how much weight was absorbed due to survival biasing
  double wgt_absorb = settings::survival_biasing
                        ? p.wgt_last() *
                            p.neutron_xs(p.event_nuclide()).absorption /
                            p.neutron_xs(p.event_nuclide()).total
                        : 0.0;

  using Type = ParticleType;

  for (auto i = 0; i < tally.scores_.size(); ++i) {
    auto score_bin = tally.scores_[i];
    auto score_index = start_index + i;
    double score = 0.0;

    switch (score_bin) {
    case SCORE_FLUX:
      // All events score to a flux bin. We actually use a collision estimator
      // in place of an analog one since there is no way to count 'events'
      // exactly for the flux
      if (p.type() == Type::neutron || p.type() == Type::photon) {
        score = flux * p.wgt_last() / p.macro_xs().total;
      } else {
        score = 0.;
      }
      break;

    case SCORE_TOTAL:
      // All events will score to the total reaction rate. We can just use
      // use the weight of the particle entering the collision as the score
      score = p.wgt_last() * flux;
      break;

    case SCORE_INVERSE_VELOCITY:
      if (p.type() != Type::neutron)
        continue;

      // All events score to an inverse velocity bin. We actually use a
      // collision estimator in place of an analog one since there is no way
      // to count 'events' exactly for the inverse velocity
      score = flux * p.wgt_last() / (p.macro_xs().total * p.speed());
      break;

    case SCORE_SCATTER:
      if (p.type() != Type::neutron && p.type() != Type::photon)
        continue;

      // Skip any event where the particle didn't scatter
      if (p.event() != TallyEvent::SCATTER)
        continue;
      // Since only scattering events make it here, again we can use the
      // weight entering the collision as the estimator for the reaction rate
      score = (p.wgt_last() - wgt_absorb) * flux;
      break;

    case SCORE_NU_SCATTER:
      if (p.type() != Type::neutron)
        continue;

      // Only analog estimators are available.
      // Skip any event where the particle didn't scatter
      if (p.event() != TallyEvent::SCATTER)
        continue;
      // For scattering production, we need to use the pre-collision weight
      // times the yield as the estimate for the number of neutrons exiting a
      // reaction with neutrons in the exit channel
      score = (p.wgt_last() - wgt_absorb) * flux;

      // Don't waste time on very common reactions we know have multiplicities
      // of one.
      if (p.event_mt() != ELASTIC && p.event_mt() != N_LEVEL &&
          !(p.event_mt() >= N_N1 && p.event_mt() <= N_NC)) {
        // Get yield and apply to score
        auto m =
          data::nuclides[p.event_nuclide()]->reaction_index_[p.event_mt()];
        const auto& rxn {*data::nuclides[p.event_nuclide()]->reactions_[m]};
        score *= (*rxn.products_[0].yield_)(E);
      }
      break;

    case SCORE_ABSORPTION:
      if (p.type() != Type::neutron && p.type() != Type::photon)
        continue;

      if (settings::survival_biasing) {
        // No absorption events actually occur if survival biasing is on --
        // just use weight absorbed in survival biasing
        score = wgt_absorb * flux;
      } else {
        // Skip any event where the particle wasn't absorbed
        if (p.event() == TallyEvent::SCATTER)
          continue;
        // All fission and absorption events will contribute here, so we
        // can just use the particle's weight entering the collision
        score = p.wgt_last() * flux;
      }
      break;

    case SCORE_FISSION:
      if (p.macro_xs().fission == 0)
        continue;
      if (settings::survival_biasing) {
        // No fission events occur if survival biasing is on -- use collision
        // estimator instead
        if (p.neutron_xs(p.event_nuclide()).total > 0) {
          score = p.wgt_last() * p.neutron_xs(p.event_nuclide()).fission /
                  p.neutron_xs(p.event_nuclide()).total * flux;
        } else {
          score = 0.;
        }
      } else {
        // Skip any non-absorption events
        if (p.event() == TallyEvent::SCATTER)
          continue;
        // All fission events will contribute, so again we can use particle's
        // weight entering the collision as the estimate for the fission
        // reaction rate
        score = p.wgt_last() * p.neutron_xs(p.event_nuclide()).fission /
                p.neutron_xs(p.event_nuclide()).absorption * flux;
      }
      break;

    case SCORE_NU_FISSION:
      if (p.macro_xs().fission == 0)
        continue;
      if (settings::survival_biasing || p.fission()) {
        if (tally.energyout_filter_ != C_NONE) {
          // Fission has multiple outgoing neutrons so this helper function
          // is used to handle scoring the multiple filter bins.
          score_fission_eout(p, i_tally, score_index, score_bin);
          continue;
        }
      }
      if (settings::survival_biasing) {
        // No fission events occur if survival biasing is on -- use collision
        // estimator instead
        if (p.neutron_xs(p.event_nuclide()).total > 0) {
          score = p.wgt_last() * p.neutron_xs(p.event_nuclide()).nu_fission /
                  p.neutron_xs(p.event_nuclide()).total * flux;
        } else {
          score = 0.;
        }
      } else {
        // Skip any non-fission events
        if (!p.fission())
          continue;
        // If there is no outgoing energy filter, than we only need to score
        // to one bin. For the score to be 'analog', we need to score the
        // number of particles that were banked in the fission bank. Since
        // this was weighted by 1/keff, we multiply by keff to get the proper
        // score.
        score = simulation::keff * p.wgt_bank() * flux;
      }
      break;

    case SCORE_PROMPT_NU_FISSION:
      if (p.macro_xs().fission == 0)
        continue;
      if (settings::survival_biasing || p.fission()) {
        if (tally.energyout_filter_ != C_NONE) {
          // Fission has multiple outgoing neutrons so this helper function
          // is used to handle scoring the multiple filter bins.
          score_fission_eout(p, i_tally, score_index, score_bin);
          continue;
        }
      }
      if (settings::survival_biasing) {
        // No fission events occur if survival biasing is on -- need to
        // calculate fraction of absorptions that would have resulted in
        // prompt-nu-fission
        if (p.neutron_xs(p.event_nuclide()).total > 0) {
          score = p.wgt_last() * p.neutron_xs(p.event_nuclide()).fission *
                  data::nuclides[p.event_nuclide()]->nu(
                    E, ReactionProduct::EmissionMode::prompt) /
                  p.neutron_xs(p.event_nuclide()).total * flux;
        } else {
          score = 0.;
        }
      } else {
        // Skip any non-fission events
        if (!p.fission())
          continue;
        // If there is no outgoing energy filter, than we only need to score
        // to one bin. For the score to be 'analog', we need to score the
        // number of particles that were banked in the fission bank. Since
        // this was weighted by 1/keff, we multiply by keff to get the proper
        // score.
        auto n_delayed = std::accumulate(
          p.n_delayed_bank(), p.n_delayed_bank() + MAX_DELAYED_GROUPS, 0);
        auto prompt_frac = 1. - n_delayed / static_cast<double>(p.n_bank());
        score = simulation::keff * p.wgt_bank() * prompt_frac * flux;
      }
      break;

    case SCORE_DELAYED_NU_FISSION:
      if (p.macro_xs().fission == 0)
        continue;
      if (settings::survival_biasing || p.fission()) {
        if (tally.energyout_filter_ != C_NONE) {
          // Fission has multiple outgoing neutrons so this helper function
          // is used to handle scoring the multiple filter bins.
          score_fission_eout(p, i_tally, score_index, score_bin);
          continue;
        }
      }
      if (settings::survival_biasing) {
        // No fission events occur if survival biasing is on -- need to
        // calculate fraction of absorptions that would have resulted in
        // delayed-nu-fission
        if (p.neutron_xs(p.event_nuclide()).total > 0 &&
            data::nuclides[p.event_nuclide()]->fissionable_) {
          if (tally.delayedgroup_filter_ != C_NONE) {
            auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
            const DelayedGroupFilter& filt {*dynamic_cast<DelayedGroupFilter*>(
              model::tally_filters[i_dg_filt].get())};
            // Tally each delayed group bin individually
            for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
              auto dg = filt.groups()[d_bin];
              auto yield = data::nuclides[p.event_nuclide()]->nu(
                E, ReactionProduct::EmissionMode::delayed, dg);
              score = p.wgt_last() * yield *
                      p.neutron_xs(p.event_nuclide()).fission /
                      p.neutron_xs(p.event_nuclide()).total * flux;
              score_fission_delayed_dg(
                i_tally, d_bin, score, score_index, p.filter_matches());
            }
            continue;
          } else {
            // If the delayed group filter is not present, compute the score
            // by multiplying the absorbed weight by the fraction of the
            // delayed-nu-fission xs to the absorption xs
            score = p.wgt_last() * p.neutron_xs(p.event_nuclide()).fission *
                    data::nuclides[p.event_nuclide()]->nu(
                      E, ReactionProduct::EmissionMode::delayed) /
                    p.neutron_xs(p.event_nuclide()).total * flux;
          }
        }
      } else {
        // Skip any non-fission events
        if (!p.fission())
          continue;
        // If there is no outgoing energy filter, than we only need to score
        // to one bin. For the score to be 'analog', we need to score the
        // number of particles that were banked in the fission bank. Since
        // this was weighted by 1/keff, we multiply by keff to get the proper
        // score. Loop over the neutrons produced from fission and check which
        // ones are delayed. If a delayed neutron is encountered, add its
        // contribution to the fission bank to the score.
        if (tally.delayedgroup_filter_ != C_NONE) {
          auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
          const DelayedGroupFilter& filt {*dynamic_cast<DelayedGroupFilter*>(
            model::tally_filters[i_dg_filt].get())};
          // Tally each delayed group bin individually
          for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
            auto d = filt.groups()[d_bin];
            score = simulation::keff * p.wgt_bank() / p.n_bank() *
                    p.n_delayed_bank(d - 1) * flux;
            score_fission_delayed_dg(
              i_tally, d_bin, score, score_index, p.filter_matches());
          }
          continue;
        } else {
          // Add the contribution from all delayed groups
          auto n_delayed = std::accumulate(
            p.n_delayed_bank(), p.n_delayed_bank() + MAX_DELAYED_GROUPS, 0);
          score =
            simulation::keff * p.wgt_bank() / p.n_bank() * n_delayed * flux;
        }
      }
      break;

    case SCORE_DECAY_RATE:
      if (p.macro_xs().fission == 0)
        continue;
      if (settings::survival_biasing) {
        // No fission events occur if survival biasing is on -- need to
        // calculate fraction of absorptions that would have resulted in
        // delayed-nu-fission
        const auto& nuc {*data::nuclides[p.event_nuclide()]};
        if (p.neutron_xs(p.event_nuclide()).total > 0 && nuc.fissionable_) {
          const auto& rxn {*nuc.fission_rx_[0]};
          if (tally.delayedgroup_filter_ != C_NONE) {
            auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
            const DelayedGroupFilter& filt {*dynamic_cast<DelayedGroupFilter*>(
              model::tally_filters[i_dg_filt].get())};
            // Tally each delayed group bin individually
            for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
              auto d = filt.groups()[d_bin];
              auto yield = nuc.nu(E, ReactionProduct::EmissionMode::delayed, d);
              auto rate = rxn.products_[d].decay_rate_;
              score = p.wgt_last() * yield *
                      p.neutron_xs(p.event_nuclide()).fission /
                      p.neutron_xs(p.event_nuclide()).total * rate * flux;
              score_fission_delayed_dg(
                i_tally, d_bin, score, score_index, p.filter_matches());
            }
            continue;
          } else {
            // If the delayed group filter is not present, compute the score
            // by multiplying the absorbed weight by the fraction of the
            // delayed-nu-fission xs to the absorption xs for all delayed
            // groups
            score = 0.;
            // We need to be careful not to overshoot the number of
            // delayed groups since this could cause the range of the
            // rxn.products_ array to be exceeded. Hence, we use the size
            // of this array and not the MAX_DELAYED_GROUPS constant for
            // this loop.
            for (auto d = 1; d < rxn.products_.size(); ++d) {
              const auto& product = rxn.products_[d];
              if (product.particle_ != Type::neutron)
                continue;

              auto yield = nuc.nu(E, ReactionProduct::EmissionMode::delayed, d);
              auto rate = product.decay_rate_;
              score += rate * p.wgt_last() *
                       p.neutron_xs(p.event_nuclide()).fission * yield /
                       p.neutron_xs(p.event_nuclide()).total * flux;
            }
          }
        }
      } else {
        // Skip any non-fission events
        if (!p.fission())
          continue;
        // If there is no outgoing energy filter, than we only need to score
        // to one bin. For the score to be 'analog', we need to score the
        // number of particles that were banked in the fission bank. Since
        // this was weighted by 1/keff, we multiply by keff to get the proper
        // score. Loop over the neutrons produced from fission and check which
        // ones are delayed. If a delayed neutron is encountered, add its
        // contribution to the fission bank to the score.
        score = 0.;
        for (auto i = 0; i < p.n_bank(); ++i) {
          const auto& bank = p.nu_bank(i);
          auto g = bank.delayed_group;
          if (g != 0) {
            const auto& nuc {*data::nuclides[p.event_nuclide()]};
            const auto& rxn {*nuc.fission_rx_[0]};
            auto rate = rxn.products_[g].decay_rate_;
            score += simulation::keff * bank.wgt * rate * flux;
            if (tally.delayedgroup_filter_ != C_NONE) {
              auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
              const DelayedGroupFilter& filt {
                *dynamic_cast<DelayedGroupFilter*>(
                  model::tally_filters[i_dg_filt].get())};
              // Find the corresponding filter bin and then score
              for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                auto d = filt.groups()[d_bin];
                if (d == g)
                  score_fission_delayed_dg(
                    i_tally, d_bin, score, score_index, p.filter_matches());
              }
              score = 0.;
            }
          }
        }
      }
      break;

    case SCORE_KAPPA_FISSION:
      if (p.macro_xs().fission == 0.)
        continue;
      score = 0.;
      // Kappa-fission values are determined from the Q-value listed for the
      // fission cross section.
      if (settings::survival_biasing) {
        // No fission events occur if survival biasing is on -- need to
        // calculate fraction of absorptions that would have resulted in
        // fission scaled by the Q-value
        const auto& nuc {*data::nuclides[p.event_nuclide()]};
        if (p.neutron_xs(p.event_nuclide()).total > 0 && nuc.fissionable_) {
          const auto& rxn {*nuc.fission_rx_[0]};
          score = p.wgt_last() * rxn.q_value_ *
                  p.neutron_xs(p.event_nuclide()).fission /
                  p.neutron_xs(p.event_nuclide()).total * flux;
        }
      } else {
        // Skip any non-absorption events
        if (p.event() == TallyEvent::SCATTER)
          continue;
        // All fission events will contribute, so again we can use particle's
        // weight entering the collision as the estimate for the fission
        // reaction rate
        const auto& nuc {*data::nuclides[p.event_nuclide()]};
        if (p.neutron_xs(p.event_nuclide()).absorption > 0 &&
            nuc.fissionable_) {
          const auto& rxn {*nuc.fission_rx_[0]};
          score = p.wgt_last() * rxn.q_value_ *
                  p.neutron_xs(p.event_nuclide()).fission /
                  p.neutron_xs(p.event_nuclide()).absorption * flux;
        }
      }
      break;

    case SCORE_EVENTS:
// Simply count the number of scoring events
#pragma omp atomic
      tally.results_(filter_index, score_index, TallyResult::VALUE) += 1.0;
      continue;

    case ELASTIC:
      if (p.type() != Type::neutron)
        continue;

      // Check if event MT matches
      if (p.event_mt() != ELASTIC)
        continue;
      score = (p.wgt_last() - wgt_absorb) * flux;
      break;

    case SCORE_FISS_Q_PROMPT:
    case SCORE_FISS_Q_RECOV:
      if (p.macro_xs().fission == 0.)
        continue;
      score =
        score_fission_q(p, score_bin, tally, flux, i_nuclide, atom_density);
      break;

    case N_2N:
    case N_3N:
    case N_4N:
    case N_GAMMA:
    case N_P:
    case N_A:
      // This case block only works if cross sections for these reactions have
      // been precalculated. When they are not, we revert to the default case,
      // which looks up cross sections
      if (!simulation::need_depletion_rx)
        goto default_case;

      if (p.type() != Type::neutron)
        continue;

      // Check if the event MT matches
      if (p.event_mt() != score_bin)
        continue;
      score = (p.wgt_last() - wgt_absorb) * flux;
      break;

    case COHERENT:
    case INCOHERENT:
    case PHOTOELECTRIC:
    case PAIR_PROD:
      if (p.type() != Type::photon)
        continue;

      if (score_bin == PHOTOELECTRIC) {
        // Photoelectric events are assigned an MT value corresponding to the
        // shell cross section. Also, photons below the energy cutoff are
        // assumed to have been absorbed via photoelectric absorption
        if ((p.event_mt() < 534 || p.event_mt() > 572) &&
            p.event_mt() != REACTION_NONE)
          continue;
      } else {
        if (p.event_mt() != score_bin)
          continue;
      }
      score = p.wgt_last() * flux;
      break;

    case HEATING:
      score = score_particle_heating(
        p, tally, flux, HEATING, i_nuclide, atom_density);
      break;

    default:
    default_case:

      // The default block is really only meant for redundant neutron reactions
      // (e.g. 444, 901)
      if (p.type() != Type::neutron)
        continue;

      // Any other score is assumed to be a MT number. Thus, we just need
      // to check if it matches the MT number of the event
      if (p.event_mt() != score_bin)
        continue;
      score = (p.wgt_last() - wgt_absorb) * flux;
    }

    // Add derivative information on score for differential tallies.
    if (tally.deriv_ != C_NONE)
      apply_derivative_to_score(
        p, i_tally, i_nuclide, atom_density, score_bin, score);

// Update tally results
#pragma omp atomic
    tally.results_(filter_index, score_index, TallyResult::VALUE) +=
      score * filter_weight;
  }
}

//! Update tally results for multigroup tallies with any estimator.
//
//! For analog tallies, the flux estimate depends on the score type so the flux
//! argument is really just used for filter weights.

void score_general_mg(Particle& p, int i_tally, int start_index,
  int filter_index, double filter_weight, int i_nuclide, double atom_density,
  double flux)
{

  auto& tally {*model::tallies[i_tally]};

  // Set the direction and group to use with get_xs
  Direction p_u;
  int p_g;
  double wgt_absorb = 0.0;
  if (tally.estimator_ == TallyEstimator::ANALOG ||
      tally.estimator_ == TallyEstimator::COLLISION) {
    if (settings::survival_biasing) {
      // Determine weight that was absorbed
      wgt_absorb = p.wgt_last() * p.neutron_xs(p.event_nuclide()).absorption /
                   p.neutron_xs(p.event_nuclide()).total;

      // Then we either are alive and had a scatter (and so g changed),
      // or are dead and g did not change
      if (p.alive()) {
        p_u = p.u_last();
        p_g = p.g_last();
      } else {
        p_u = p.u_local();
        p_g = p.g();
      }
    } else if (p.event() == TallyEvent::SCATTER) {

      // Then the energy group has been changed by the scattering routine
      // meaning gin is now in p % last_g
      p_u = p.u_last();
      p_g = p.g_last();
    } else {

      // No scatter, no change in g.
      p_u = p.u_local();
      p_g = p.g();
    }
  } else {

    // No actual collision so g has not changed.
    p_u = p.u_local();
    p_g = p.g();
  }

  // For shorthand, assign pointers to the material and nuclide xs set
  auto& nuc_xs = (i_nuclide >= 0) ? data::mg.nuclides_[i_nuclide]
                                  : data::mg.macro_xs_[p.material()];
  auto& macro_xs = data::mg.macro_xs_[p.material()];

  // Find the temperature and angle indices of interest
  int macro_t = p.mg_xs_cache().t;
  int macro_a = macro_xs.get_angle_index(p_u);
  int nuc_t = 0;
  int nuc_a = 0;
  if (i_nuclide >= 0) {
    nuc_t = nuc_xs.get_temperature_index(p.sqrtkT());
    nuc_a = nuc_xs.get_angle_index(p_u);
  }

  for (auto i = 0; i < tally.scores_.size(); ++i) {
    auto score_bin = tally.scores_[i];
    auto score_index = start_index + i;

    double score;

    switch (score_bin) {

    case SCORE_FLUX:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // All events score to a flux bin. We actually use a collision estimator
        // in place of an analog one since there is no way to count 'events'
        // exactly for the flux
        score = flux * p.wgt_last() / p.macro_xs().total;
      } else {
        score = flux;
      }
      break;

    case SCORE_TOTAL:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // All events will score to the total reaction rate. We can just use
        // use the weight of the particle entering the collision as the score
        score = flux * p.wgt_last();
        if (i_nuclide >= 0) {
          score *= atom_density *
                   nuc_xs.get_xs(MgxsType::TOTAL, p_g, nuc_t, nuc_a) /
                   macro_xs.get_xs(MgxsType::TOTAL, p_g, macro_t, macro_a);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux *
                  nuc_xs.get_xs(MgxsType::TOTAL, p_g, nuc_t, nuc_a);
        } else {
          score = p.macro_xs().total * flux;
        }
      }
      break;

    case SCORE_INVERSE_VELOCITY:
      if (tally.estimator_ == TallyEstimator::ANALOG ||
          tally.estimator_ == TallyEstimator::COLLISION) {
        // All events score to an inverse velocity bin. We actually use a
        // collision estimator in place of an analog one since there is no way
        // to count 'events' exactly for the inverse velocity
        score = flux * p.wgt_last();
        if (i_nuclide >= 0) {
          score *=
            nuc_xs.get_xs(MgxsType::INVERSE_VELOCITY, p_g, nuc_t, nuc_a) /
            macro_xs.get_xs(MgxsType::TOTAL, p_g, macro_t, macro_a);
        } else {
          score *=
            macro_xs.get_xs(MgxsType::INVERSE_VELOCITY, p_g, macro_t, macro_a) /
            macro_xs.get_xs(MgxsType::TOTAL, p_g, macro_t, macro_a);
        }
      } else {
        if (i_nuclide >= 0) {
          score =
            flux * nuc_xs.get_xs(MgxsType::INVERSE_VELOCITY, p_g, nuc_t, nuc_a);
        } else {
          score = flux * macro_xs.get_xs(
                           MgxsType::INVERSE_VELOCITY, p_g, macro_t, macro_a);
        }
      }
      break;

    case SCORE_SCATTER:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // Skip any event where the particle didn't scatter
        if (p.event() != TallyEvent::SCATTER)
          continue;
        // Since only scattering events make it here, again we can use the
        // weight entering the collision as the estimator for the reaction rate
        score = (p.wgt_last() - wgt_absorb) * flux;
        if (i_nuclide >= 0) {
          score *= atom_density *
                   nuc_xs.get_xs(MgxsType::SCATTER_FMU, p.g_last(), &p.g(),
                     &p.mu(), nullptr, nuc_t, nuc_a) /
                   macro_xs.get_xs(MgxsType::SCATTER_FMU, p.g_last(), &p.g(),
                     &p.mu(), nullptr, macro_t, macro_a);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux *
                  nuc_xs.get_xs(MgxsType::SCATTER, p_g, nullptr, &p.mu(),
                    nullptr, nuc_t, nuc_a);
        } else {
          score = flux * macro_xs.get_xs(MgxsType::SCATTER, p_g, nullptr,
                           &p.mu(), nullptr, macro_t, macro_a);
        }
      }
      break;

    case SCORE_NU_SCATTER:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // Skip any event where the particle didn't scatter
        if (p.event() != TallyEvent::SCATTER)
          continue;
        // For scattering production, we need to use the pre-collision weight
        // times the multiplicity as the estimate for the number of neutrons
        // exiting a reaction with neutrons in the exit channel
        score = (p.wgt_last() - wgt_absorb) * flux;
        // Since we transport based on material data, the angle selected
        // was not selected from the f(mu) for the nuclide.  Therefore
        // adjust the score by the actual probability for that nuclide.
        if (i_nuclide >= 0) {
          score *= atom_density *
                   nuc_xs.get_xs(MgxsType::NU_SCATTER_FMU, p.g_last(), &p.g(),
                     &p.mu(), nullptr, nuc_t, nuc_a) /
                   macro_xs.get_xs(MgxsType::NU_SCATTER_FMU, p.g_last(), &p.g(),
                     &p.mu(), nullptr, macro_t, macro_a);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux *
                  nuc_xs.get_xs(MgxsType::NU_SCATTER, p_g, nuc_t, nuc_a);
        } else {
          score =
            flux * macro_xs.get_xs(MgxsType::NU_SCATTER, p_g, macro_t, macro_a);
        }
      }
      break;

    case SCORE_ABSORPTION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No absorption events actually occur if survival biasing is on --
          // just use weight absorbed in survival biasing
          score = wgt_absorb * flux;
        } else {
          // Skip any event where the particle wasn't absorbed
          if (p.event() == TallyEvent::SCATTER)
            continue;
          // All fission and absorption events will contribute here, so we
          // can just use the particle's weight entering the collision
          score = p.wgt_last() * flux;
        }
        if (i_nuclide >= 0) {
          score *= atom_density *
                   nuc_xs.get_xs(MgxsType::ABSORPTION, p_g, nuc_t, nuc_a) /
                   macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux *
                  nuc_xs.get_xs(MgxsType::ABSORPTION, p_g, nuc_t, nuc_a);
        } else {
          score = p.macro_xs().absorption * flux;
        }
      }
      break;

    case SCORE_FISSION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // fission
          score = wgt_absorb * flux;
        } else {
          // Skip any non-absorption events
          if (p.event() == TallyEvent::SCATTER)
            continue;
          // All fission events will contribute, so again we can use particle's
          // weight entering the collision as the estimate for the fission
          // reaction rate
          score = p.wgt_last() * flux;
        }
        if (i_nuclide >= 0) {
          score *= atom_density *
                   nuc_xs.get_xs(MgxsType::FISSION, p_g, nuc_t, nuc_a) /
                   macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
        } else {
          score *= macro_xs.get_xs(MgxsType::FISSION, p_g, macro_t, macro_a) /
                   macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux *
                  nuc_xs.get_xs(MgxsType::FISSION, p_g, nuc_t, nuc_a);
        } else {
          score =
            flux * macro_xs.get_xs(MgxsType::FISSION, p_g, macro_t, macro_a);
        }
      }
      break;

    case SCORE_NU_FISSION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing || p.fission()) {
          if (tally.energyout_filter_ != C_NONE) {
            // Fission has multiple outgoing neutrons so this helper function
            // is used to handle scoring the multiple filter bins.
            score_fission_eout(p, i_tally, score_index, score_bin);
            continue;
          }
        }
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // nu-fission
          score = wgt_absorb * flux;
          if (i_nuclide >= 0) {
            score *=
              atom_density *
              nuc_xs.get_xs(MgxsType::NU_FISSION, p_g, nuc_t, nuc_a) /
              macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
          } else {
            score *=
              macro_xs.get_xs(MgxsType::NU_FISSION, p_g, macro_t, macro_a) /
              macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
          }
        } else {
          // Skip any non-fission events
          if (!p.fission())
            continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score.
          score = simulation::keff * p.wgt_bank() * flux;
          if (i_nuclide >= 0) {
            score *= atom_density *
                     nuc_xs.get_xs(MgxsType::FISSION, p_g, nuc_t, nuc_a) /
                     macro_xs.get_xs(MgxsType::FISSION, p_g, macro_t, macro_a);
          }
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux *
                  nuc_xs.get_xs(MgxsType::NU_FISSION, p_g, nuc_t, nuc_a);
        } else {
          score =
            flux * macro_xs.get_xs(MgxsType::NU_FISSION, p_g, macro_t, macro_a);
        }
      }
      break;

    case SCORE_PROMPT_NU_FISSION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing || p.fission()) {
          if (tally.energyout_filter_ != C_NONE) {
            // Fission has multiple outgoing neutrons so this helper function
            // is used to handle scoring the multiple filter bins.
            score_fission_eout(p, i_tally, score_index, score_bin);
            continue;
          }
        }
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // prompt-nu-fission
          score = wgt_absorb * flux;
          if (i_nuclide >= 0) {
            score *=
              atom_density *
              nuc_xs.get_xs(MgxsType::PROMPT_NU_FISSION, p_g, nuc_t, nuc_a) /
              macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
          } else {
            score *=
              macro_xs.get_xs(
                MgxsType::PROMPT_NU_FISSION, p_g, macro_t, macro_a) /
              macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
          }
        } else {
          // Skip any non-fission events
          if (!p.fission())
            continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score.
          auto n_delayed = std::accumulate(
            p.n_delayed_bank(), p.n_delayed_bank() + MAX_DELAYED_GROUPS, 0);
          auto prompt_frac = 1. - n_delayed / static_cast<double>(p.n_bank());
          score = simulation::keff * p.wgt_bank() * prompt_frac * flux;
          if (i_nuclide >= 0) {
            score *= atom_density *
                     nuc_xs.get_xs(MgxsType::FISSION, p_g, nuc_t, nuc_a) /
                     macro_xs.get_xs(MgxsType::FISSION, p_g, macro_t, macro_a);
          }
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux *
                  nuc_xs.get_xs(MgxsType::PROMPT_NU_FISSION, p_g, nuc_t, nuc_a);
        } else {
          score = flux * macro_xs.get_xs(
                           MgxsType::PROMPT_NU_FISSION, p_g, macro_t, macro_a);
        }
      }
      break;

    case SCORE_DELAYED_NU_FISSION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing || p.fission()) {
          if (tally.energyout_filter_ != C_NONE) {
            // Fission has multiple outgoing neutrons so this helper function
            // is used to handle scoring the multiple filter bins.
            score_fission_eout(p, i_tally, score_index, score_bin);
            continue;
          }
        }
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // delayed-nu-fission
          double abs_xs =
            macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
          if (abs_xs > 0.) {
            if (tally.delayedgroup_filter_ != C_NONE) {
              auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
              const DelayedGroupFilter& filt {
                *dynamic_cast<DelayedGroupFilter*>(
                  model::tally_filters[i_dg_filt].get())};
              // Tally each delayed group bin individually
              for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                auto d = filt.groups()[d_bin] - 1;
                score = wgt_absorb * flux;
                if (i_nuclide >= 0) {
                  score *= nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d, nuc_t, nuc_a) /
                           abs_xs;
                } else {
                  score *= macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d, macro_t, macro_a) /
                           abs_xs;
                }
                score_fission_delayed_dg(
                  i_tally, d_bin, score, score_index, p.filter_matches());
              }
              continue;
            } else {
              // If the delayed group filter is not present, compute the score
              // by multiplying the absorbed weight by the fraction of the
              // delayed-nu-fission xs to the absorption xs
              score = wgt_absorb * flux;
              if (i_nuclide >= 0) {
                score *= nuc_xs.get_xs(
                           MgxsType::DELAYED_NU_FISSION, p_g, nuc_t, nuc_a) /
                         abs_xs;
              } else {
                score *= macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                           macro_t, macro_a) /
                         abs_xs;
              }
            }
          }
        } else {
          // Skip any non-fission events
          if (!p.fission())
            continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score. Loop over the neutrons produced from fission and check which
          // ones are delayed. If a delayed neutron is encountered, add its
          // contribution to the fission bank to the score.
          if (tally.delayedgroup_filter_ != C_NONE) {
            auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
            const DelayedGroupFilter& filt {*dynamic_cast<DelayedGroupFilter*>(
              model::tally_filters[i_dg_filt].get())};
            // Tally each delayed group bin individually
            for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
              auto d = filt.groups()[d_bin];
              score = simulation::keff * p.wgt_bank() / p.n_bank() *
                      p.n_delayed_bank(d - 1) * flux;
              if (i_nuclide >= 0) {
                score *=
                  atom_density *
                  nuc_xs.get_xs(MgxsType::FISSION, p_g, nuc_t, nuc_a) /
                  macro_xs.get_xs(MgxsType::FISSION, p_g, macro_t, macro_a);
              }
              score_fission_delayed_dg(
                i_tally, d_bin, score, score_index, p.filter_matches());
            }
            continue;
          } else {
            // Add the contribution from all delayed groups
            auto n_delayed = std::accumulate(
              p.n_delayed_bank(), p.n_delayed_bank() + MAX_DELAYED_GROUPS, 0);
            score =
              simulation::keff * p.wgt_bank() / p.n_bank() * n_delayed * flux;
            if (i_nuclide >= 0) {
              score *=
                atom_density *
                nuc_xs.get_xs(MgxsType::FISSION, p_g, nuc_t, nuc_a) /
                macro_xs.get_xs(MgxsType::FISSION, p_g, macro_t, macro_a);
            }
          }
        }
      } else {
        if (tally.delayedgroup_filter_ != C_NONE) {
          auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
          const DelayedGroupFilter& filt {*dynamic_cast<DelayedGroupFilter*>(
            model::tally_filters[i_dg_filt].get())};
          // Tally each delayed group bin individually
          for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
            auto d = filt.groups()[d_bin] - 1;
            if (i_nuclide >= 0) {
              score = flux * atom_density *
                      nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g, nullptr,
                        nullptr, &d, nuc_t, nuc_a);
            } else {
              score = flux * macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                               nullptr, nullptr, &d, macro_t, macro_a);
            }
            score_fission_delayed_dg(
              i_tally, d_bin, score, score_index, p.filter_matches());
          }
          continue;
        } else {
          if (i_nuclide >= 0) {
            score =
              flux * atom_density *
              nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g, nuc_t, nuc_a);
          } else {
            score = flux * macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             macro_t, macro_a);
          }
        }
      }
      break;

    case SCORE_DECAY_RATE:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // delayed-nu-fission
          double abs_xs =
            macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
          if (abs_xs > 0) {
            if (tally.delayedgroup_filter_ != C_NONE) {
              auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
              const DelayedGroupFilter& filt {
                *dynamic_cast<DelayedGroupFilter*>(
                  model::tally_filters[i_dg_filt].get())};
              // Tally each delayed group bin individually
              for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                auto d = filt.groups()[d_bin] - 1;
                score = wgt_absorb * flux;
                if (i_nuclide >= 0) {
                  score *= nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                             nullptr, &d, nuc_t, nuc_a) *
                           nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d, nuc_t, nuc_a) /
                           abs_xs;
                } else {
                  score *= macro_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                             nullptr, &d, macro_t, macro_a) *
                           macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d, macro_t, macro_a) /
                           abs_xs;
                }
                score_fission_delayed_dg(
                  i_tally, d_bin, score, score_index, p.filter_matches());
              }
              continue;
            } else {
              // If the delayed group filter is not present, compute the score
              // by multiplying the absorbed weight by the fraction of the
              // delayed-nu-fission xs to the absorption xs for all delayed
              // groups
              score = 0.;
              for (auto d = 0; d < data::mg.num_delayed_groups_; ++d) {
                if (i_nuclide >= 0) {
                  score += wgt_absorb * flux *
                           nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                             nullptr, &d, nuc_t, nuc_a) *
                           nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d, nuc_t, nuc_a) /
                           abs_xs;
                } else {
                  score += wgt_absorb * flux *
                           macro_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                             nullptr, &d, macro_t, macro_a) *
                           macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d, macro_t, macro_a) /
                           abs_xs;
                }
              }
            }
          }
        } else {
          // Skip any non-fission events
          if (!p.fission())
            continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score. Loop over the neutrons produced from fission and check which
          // ones are delayed. If a delayed neutron is encountered, add its
          // contribution to the fission bank to the score.
          score = 0.;
          for (auto i = 0; i < p.n_bank(); ++i) {
            const auto& bank = p.nu_bank(i);
            auto d = bank.delayed_group - 1;
            if (d != -1) {
              if (i_nuclide >= 0) {
                score +=
                  simulation::keff * atom_density * bank.wgt * flux *
                  nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d,
                    nuc_t, nuc_a) *
                  nuc_xs.get_xs(MgxsType::FISSION, p_g, nuc_t, nuc_a) /
                  macro_xs.get_xs(MgxsType::FISSION, p_g, macro_t, macro_a);
              } else {
                score += simulation::keff * bank.wgt * flux *
                         macro_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                           nullptr, &d, macro_t, macro_a);
              }
              if (tally.delayedgroup_filter_ != C_NONE) {
                auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
                const DelayedGroupFilter& filt {
                  *dynamic_cast<DelayedGroupFilter*>(
                    model::tally_filters[i_dg_filt].get())};
                // Find the corresponding filter bin and then score
                for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                  auto dg = filt.groups()[d_bin];
                  if (dg == d + 1)
                    score_fission_delayed_dg(
                      i_tally, d_bin, score, score_index, p.filter_matches());
                }
                score = 0.;
              }
            }
          }
          if (tally.delayedgroup_filter_ != C_NONE)
            continue;
        }
      } else {
        if (tally.delayedgroup_filter_ != C_NONE) {
          auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
          const DelayedGroupFilter& filt {*dynamic_cast<DelayedGroupFilter*>(
            model::tally_filters[i_dg_filt].get())};
          // Tally each delayed group bin individually
          for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
            auto d = filt.groups()[d_bin] - 1;
            if (i_nuclide >= 0) {
              score = atom_density * flux *
                      nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr, nullptr,
                        &d, nuc_t, nuc_a) *
                      nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g, nullptr,
                        nullptr, &d, nuc_t, nuc_a);
            } else {
              score = flux *
                      macro_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                        nullptr, &d, macro_t, macro_a) *
                      macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                        nullptr, nullptr, &d, macro_t, macro_a);
            }
            score_fission_delayed_dg(
              i_tally, d_bin, score, score_index, p.filter_matches());
          }
          continue;
        } else {
          score = 0.;
          for (auto d = 0; d < data::mg.num_delayed_groups_; ++d) {
            if (i_nuclide >= 0) {
              score += atom_density * flux *
                       nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                         nullptr, &d, nuc_t, nuc_a) *
                       nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g, nullptr,
                         nullptr, &d, nuc_t, nuc_a);
            } else {
              score += flux *
                       macro_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                         nullptr, &d, macro_t, macro_a) *
                       macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                         nullptr, nullptr, &d, macro_t, macro_a);
            }
          }
        }
      }
      break;

    case SCORE_KAPPA_FISSION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // fission scaled by the Q-value
          score = wgt_absorb * flux;
        } else {
          // Skip any non-absorption events
          if (p.event() == TallyEvent::SCATTER)
            continue;
          // All fission events will contribute, so again we can use particle's
          // weight entering the collision as the estimate for the fission
          // reaction rate
          score = p.wgt_last() * flux;
        }
        if (i_nuclide >= 0) {
          score *= atom_density *
                   nuc_xs.get_xs(MgxsType::KAPPA_FISSION, p_g, nuc_t, nuc_a) /
                   macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
        } else {
          score *=
            macro_xs.get_xs(MgxsType::KAPPA_FISSION, p_g, macro_t, macro_a) /
            macro_xs.get_xs(MgxsType::ABSORPTION, p_g, macro_t, macro_a);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux *
                  nuc_xs.get_xs(MgxsType::KAPPA_FISSION, p_g, nuc_t, nuc_a);
        } else {
          score = flux * macro_xs.get_xs(
                           MgxsType::KAPPA_FISSION, p_g, macro_t, macro_a);
        }
      }
      break;

    case SCORE_EVENTS:
// Simply count the number of scoring events
#pragma omp atomic
      tally.results_(filter_index, score_index, TallyResult::VALUE) += 1.0;
      continue;

    default:
      continue;
    }

// Update tally results
#pragma omp atomic
    tally.results_(filter_index, score_index, TallyResult::VALUE) +=
      score * filter_weight;
  }
}

void score_analog_tally_ce(Particle& p)
{
  // Since electrons/positrons are not transported, we assign a flux of zero.
  // Note that the heating score does NOT use the flux and will be non-zero for
  // electrons/positrons.
  double flux =
    (p.type() == ParticleType::neutron || p.type() == ParticleType::photon)
      ? 1.0
      : 0.0;

  for (auto i_tally : model::active_analog_tallies) {
    const Tally& tally {*model::tallies[i_tally]};

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p);
    auto end = FilterBinIter(tally, true, &p.filter_matches());
    if (filter_iter == end)
      continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
      for (auto i = 0; i < tally.nuclides_.size(); ++i) {
        auto i_nuclide = tally.nuclides_[i];

        // Tally this event in the present nuclide bin if that bin represents
        // the event nuclide or the total material.  Note that the atomic
        // density argument for score_general is not used for analog tallies.
        if (i_nuclide == p.event_nuclide() || i_nuclide == -1)
          score_general_ce_analog(p, i_tally, i * tally.scores_.size(),
            filter_index, filter_weight, i_nuclide, -1.0, flux);
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate)
      break;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : p.filter_matches())
    match.bins_present_ = false;
}

void score_analog_tally_mg(Particle& p)
{
  for (auto i_tally : model::active_analog_tallies) {
    const Tally& tally {*model::tallies[i_tally]};

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p);
    auto end = FilterBinIter(tally, true, &p.filter_matches());
    if (filter_iter == end)
      continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
      for (auto i = 0; i < tally.nuclides_.size(); ++i) {
        auto i_nuclide = tally.nuclides_[i];

        double atom_density = 0.;
        if (i_nuclide >= 0) {
          auto j =
            model::materials[p.material()]->mat_nuclide_index_[i_nuclide];
          if (j == C_NONE)
            continue;
          atom_density = model::materials[p.material()]->atom_density_(j);
        }

        score_general_mg(p, i_tally, i * tally.scores_.size(), filter_index,
          filter_weight, i_nuclide, atom_density, 1.0);
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate)
      break;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : p.filter_matches())
    match.bins_present_ = false;
}

void score_tracklength_tally(Particle& p, double distance)
{
  // Determine the tracklength estimate of the flux
  double flux = p.wgt() * distance;

  // Set 'none' value for log union grid index
  int i_log_union = C_NONE;

  for (auto i_tally : model::active_tracklength_tallies) {
    const Tally& tally {*model::tallies[i_tally]};

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p);
    auto end = FilterBinIter(tally, true, &p.filter_matches());
    if (filter_iter == end)
      continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
      for (auto i = 0; i < tally.nuclides_.size(); ++i) {
        auto i_nuclide = tally.nuclides_[i];

        double atom_density = 0.;
        if (i_nuclide >= 0) {
          if (p.material() != MATERIAL_VOID) {
            const auto& mat = model::materials[p.material()];
            auto j = mat->mat_nuclide_index_[i_nuclide];
            if (j == C_NONE) {
              // Determine log union grid index
              if (i_log_union == C_NONE) {
                int neutron = static_cast<int>(ParticleType::neutron);
                i_log_union = std::log(p.E() / data::energy_min[neutron]) /
                              simulation::log_spacing;
              }

              // Update micro xs cache
              if (!tally.multiply_density()) {
                p.update_neutron_xs(i_nuclide, i_log_union);
                atom_density = 1.0;
              }
            } else {
              atom_density =
                tally.multiply_density() ? mat->atom_density_(j) : 1.0;
            }
          }
        }

        // TODO: consider replacing this "if" with pointers or templates
        if (settings::run_CE) {
          score_general_ce_nonanalog(p, i_tally, i * tally.scores_.size(),
            filter_index, filter_weight, i_nuclide, atom_density, flux);
        } else {
          score_general_mg(p, i_tally, i * tally.scores_.size(), filter_index,
            filter_weight, i_nuclide, atom_density, flux);
        }
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate)
      break;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : p.filter_matches())
    match.bins_present_ = false;
}

void score_collision_tally(Particle& p)
{
  // Determine the collision estimate of the flux
  double flux = 0.0;
  if (p.type() == ParticleType::neutron || p.type() == ParticleType::photon) {
    flux = p.wgt_last() / p.macro_xs().total;
  }

  // Set 'none value for log union grid index
  int i_log_union = C_NONE;

  for (auto i_tally : model::active_collision_tallies) {
    const Tally& tally {*model::tallies[i_tally]};

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p);
    auto end = FilterBinIter(tally, true, &p.filter_matches());
    if (filter_iter == end)
      continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
      for (auto i = 0; i < tally.nuclides_.size(); ++i) {
        auto i_nuclide = tally.nuclides_[i];

        double atom_density = 0.;
        if (i_nuclide >= 0) {
          const auto& mat = model::materials[p.material()];
          auto j = mat->mat_nuclide_index_[i_nuclide];
          if (j == C_NONE) {
            // Determine log union grid index
            if (i_log_union == C_NONE) {
              int neutron = static_cast<int>(ParticleType::neutron);
              i_log_union = std::log(p.E() / data::energy_min[neutron]) /
                            simulation::log_spacing;
            }

            // Update micro xs cache
            if (!tally.multiply_density()) {
              p.update_neutron_xs(i_nuclide, i_log_union);
              atom_density = 1.0;
            }
          } else {
            atom_density =
              tally.multiply_density() ? mat->atom_density_(j) : 1.0;
          }
        }

        // TODO: consider replacing this "if" with pointers or templates
        if (settings::run_CE) {
          score_general_ce_nonanalog(p, i_tally, i * tally.scores_.size(),
            filter_index, filter_weight, i_nuclide, atom_density, flux);
        } else {
          score_general_mg(p, i_tally, i * tally.scores_.size(), filter_index,
            filter_weight, i_nuclide, atom_density, flux);
        }
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate)
      break;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : p.filter_matches())
    match.bins_present_ = false;
}

void score_point_tally(Particle& p)
{
  if ((p.event_mt() == 101) || (p.event_mt() == 18) ||
      (p.event_index_mt() == -1234))
    return; // absorption or fission or s(a,b)

  for (auto i_tally : model::active_point_tallies) {

    double yield = 1;
    col_counter++;
    const auto& nuc {data::nuclides[p.event_nuclide()]};
    const auto& rx {nuc->reactions_[p.event_index_mt()]};
    auto& d = rx->products_[0].distribution_[0];
    auto d_ = dynamic_cast<UncorrelatedAngleEnergy*>(d.get());
    // Initialize
    std::vector<double> mu_cm;
    std::vector<double> Js;
    std::vector<Particle> ghost_particles;
    std::vector<double> pdfs_lab;
    // std::vector<double> fluxes;
    if (std::isnan(p.r().x)) {
      return; // sometimes the get_pdf function turns p.r() into nan
    }

    // Get position (x,y,z) of detector
    double det_pos[4];
    get_det_pos(det_pos, i_tally);

    if (p.event_mt() == 2 && p.event_index_mt() != -1234) {
      get_pdf_to_point_elastic(det_pos, p, mu_cm, Js, ghost_particles);
      for (std::size_t i = 0; i < mu_cm.size(); ++i) {
        // Assuming Js.size() is the same as mu_cm.size()
        double mu_c = mu_cm[i];
        double derivative = Js[i];
        double pdf_cm = d_->angle().get_pdf(p.E_last(), mu_c, p.current_seed());
        pdfs_lab.push_back(pdf_cm / std::abs(derivative));
      }
    }
    if (p.event_mt() != 2) { // Inelastic
      // make sure v_t is 0
      // copy energy of neutron
      double E_in = p.E_last();

      double E_out;
      rx->products_[0].get_pdf(i_tally, E_in, E_out, p.current_seed(), p, mu_cm,
        Js, ghost_particles, pdfs_lab);

      yield = (*rx->products_[0].yield_)(p.E_last());
      if (std::floor(yield) != yield && yield > 0) {
        yield = 1;
      }
      for (auto& p_ghost : ghost_particles) {
        p_ghost.wgt() = p_ghost.wgt() * yield;
        ;
      }
    }

    for (size_t index = 0; index < ghost_particles.size(); ++index) {
      auto& ghost_p = ghost_particles[index];
      double pdf_lab = pdfs_lab[index];
      score_ghost_particle(ghost_p, pdf_lab, i_tally);
      // calculate shielding
    } // for loop on ghost particles
  }
}

void get_det_pos(double (&det_pos)[4], int i_tally)
{
  const Tally& tally {*model::tallies[i_tally]};
  if (tally.positions_.size() == 4) {
    for (auto i = 0; i < 3; ++i) {
      auto pos_coord = tally.positions_[i];
      det_pos[i] = std::stod(pos_coord);
    }
    auto R0 = tally.positions_[3];
    det_pos[4] = std::stod(R0);
  } else {
    fatal_error("user must use 3 positions and R0");
  }
}

void score_point_tally_from_source(const SourceSite* src)
{
  if (!src->ext) {
    return;
  }
  double flux = 0;
  for (auto i_tally : model::active_point_tallies) {
    double det_pos[4]; // Get position (x,y,z) of detector
    get_det_pos(det_pos, i_tally);
    Direction u_lab {det_pos[0] - src->r.x, // towards the detector
      det_pos[1] - src->r.y, det_pos[2] - src->r.z};
    Direction u_lab_unit = u_lab / u_lab.norm(); // normalize
    // Get the angle distribution
    double pdf_mu_det = 0;
    Source& mysource = *(model::external_sources[src->source_index]);
    if (auto* independentSource = dynamic_cast<IndependentSource*>(&mysource))
    // Check if the actual type is IndependentSource
    {
      UnitSphereDistribution* angleDistribution = independentSource->angle();
      if (angleDistribution) {
        Direction my_u_ref = angleDistribution->u_ref_;
        Direction my_u_ref_unit = my_u_ref / my_u_ref.norm();
        if (typeid(*angleDistribution) == typeid(PolarAzimuthal))
        // polar and assuming phi is isotropic
        {
          double my_det_mu = my_u_ref_unit.dot(u_lab_unit);
          if (std::abs(my_det_mu) > 1.0) {
            my_det_mu = std::copysign(1.0, my_det_mu);
          }
          auto* polarAzimuthalDistribution =
            dynamic_cast<PolarAzimuthal*>(angleDistribution);
          Distribution* muDistribution = polarAzimuthalDistribution->mu();
          Distribution* phiDistribution = polarAzimuthalDistribution->phi();
          pdf_mu_det = muDistribution->get_pdf(my_det_mu);
        } else if (typeid(*angleDistribution) == typeid(Isotropic))
        // Isotropic
        {
          pdf_mu_det = 0.5;
        } else if (typeid(*angleDistribution) == typeid(Monodirectional))
        // pencil
        {
          if (my_u_ref_unit.dot(u_lab_unit) == 1) {
            double total_distance = u_lab.norm();
            pdf_mu_det =
              (2 * PI * total_distance * total_distance); // to cancel the Rs
          }    // pdf is not defined for discrete case
        } else // unexpected type
        {
          fatal_error("unexpected type");
        }
      }

      Particle ghost_particle = Particle();
      ghost_particle.initialize_ghost_particle_from_source(src, u_lab_unit);

      score_ghost_particle(ghost_particle, pdf_mu_det, i_tally);
    }
  }
}

void score_surface_tally(Particle& p, const vector<int>& tallies)
{
  double current = p.wgt_last();

  for (auto i_tally : tallies) {
    auto& tally {*model::tallies[i_tally]};

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p);
    auto end = FilterBinIter(tally, true, &p.filter_matches());
    if (filter_iter == end)
      continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over scores.
      // There is only one score type for current tallies so there is no need
      // for a further scoring function.
      double score = current * filter_weight;
      for (auto score_index = 0; score_index < tally.scores_.size();
           ++score_index) {
#pragma omp atomic
        tally.results_(filter_index, score_index, TallyResult::VALUE) += score;
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate)
      break;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : p.filter_matches())
    match.bins_present_ = false;
}

void boostf(double A[4], double B[4], double X[4])
{
  //
  //     boosts B(labfram) to A rest frame and gives output in X
  //
  double W;
  int j;

  if ((A[0] * A[0] - A[1] * A[1] - A[2] * A[2] - A[3] * A[3]) <= 0) {
  }

  W = sqrt(A[0] * A[0] - A[1] * A[1] - A[2] * A[2] - A[3] * A[3]);

  if (W == 0 || W == (-A[0]))

    X[0] = (B[0] * A[0] - B[3] * A[3] - B[2] * A[2] - B[1] * A[1]) / W;
  for (j = 1; j <= 3; j++) {
    X[j] = B[j] - A[j] * (B[0] + X[0]) / (A[0] + W);
  }

  return;
}

double get_MFP(Particle& ghost_particle, double total_distance)
{
  // calculate shilding
  double remaining_distance = total_distance;
  double total_MFP = 0;
  ghost_particle.event_calculate_xs();
  ghost_particle.boundary() = distance_to_boundary(ghost_particle);
  double advance_distance = ghost_particle.boundary().distance();

  while (advance_distance < remaining_distance) // Advance to next boundary
  {
    total_MFP += advance_distance * ghost_particle.macro_xs().total;
    // Advance particle in space and time
    for (int j = 0; j < ghost_particle.n_coord(); ++j) {
      ghost_particle.coord(j).r() +=
        advance_distance * ghost_particle.coord(j).u();
    }
    remaining_distance -= advance_distance;
    ghost_particle.time() += advance_distance / ghost_particle.speed();
    ghost_particle.event_cross_surface();
    ghost_particle.event_calculate_xs();
    ghost_particle.boundary() = distance_to_boundary(ghost_particle);
    advance_distance = ghost_particle.boundary().distance();
  }
  total_MFP += remaining_distance *
               ghost_particle.macro_xs().total; // advance to next boundary
  for (int j = 0; j < ghost_particle.n_coord(); ++j) {
    ghost_particle.coord(j).r() +=
      remaining_distance * ghost_particle.coord(j).u();
  }
  ghost_particle.time() += remaining_distance / ghost_particle.speed();
  return total_MFP;
}

void get_pdf_to_point_elastic(double det_pos[4], Particle& p,
  std::vector<double>& mu_cm, std::vector<double>& Js,
  std::vector<Particle>& ghost_particles, double E3k_cm_given)
{
  Direction u_lab {det_pos[0] - p.r().x, // towards the detector
    det_pos[1] - p.r().y, det_pos[2] - p.r().z};
  Direction u_lab_unit = u_lab / u_lab.norm(); // normalize

  double m1 = p.getMass() / 1e6; // mass of incoming particle in MeV
  const auto& nuc {data::nuclides[p.event_nuclide()]};
  double awr = nuc->awr_;
  double m2 = m1 * awr; // mass of target
  double m3 = m1;       // mass of outgoing particle to detector
  double m4 = m2;       // mass of recoil target  system

  double E1_tot =
    p.E_last() / 1e6 + m1; // total Energy of incoming particle in MeV
  double p1_tot = std::sqrt(
    E1_tot * E1_tot - m1 * m1); // total momenta of incoming particle in MeV
  // without this the get_pdf function turns p.r() into nan
  Direction p1 = p1_tot * p.u_last();    // 3 momentum of incoming particle
  Direction p2 = p.v_t() * m2 / C_LIGHT; // 3 momentum of target in lab
  double E2_tot = std::sqrt(p2.norm() * p2.norm() + m2 * m2);
  double E_cm = E1_tot + E2_tot;
  Direction p_cm = p1 + p2;
  double p_tot_cm = p_cm.norm();

  double cos_lab = u_lab_unit.dot(p_cm) / (p_tot_cm); // between cm and p3
  if (std::abs(cos_lab) > 1.0) {
    cos_lab = std::copysign(1.0, cos_lab);
  }

  double theta = std::acos(cos_lab);
  double sin_lab_sq = 1 - cos_lab * cos_lab;

  double M_cm = std::sqrt(
    E_cm * E_cm -
    p_tot_cm * p_tot_cm); // mass of the center of mass (incoming and target)
  double gamma = E_cm / M_cm;
  double p1_cm[4];
  double A[4] = {E_cm, p_cm.x, p_cm.y, p_cm.z};
  // double invA[4] = {E_cm, -p_cm.x , -p_cm.y , -p_cm.z};
  // boostf( invA ,p1_cm,  maybe_p1_lab); boost back to lab
  double B[4] = {E1_tot, p1.x, p1.y, p1.z};
  boostf(A, B, p1_cm);
  double p1_tot_cm =
    std::sqrt(p1_cm[1] * p1_cm[1] + p1_cm[2] * p1_cm[2] + p1_cm[3] * p1_cm[3]);
  double E3_cm = (M_cm * M_cm + m3 * m3 - m4 * m4) / (2 * M_cm);
  if (E3k_cm_given >= 0.0) {
    E3_cm = E3k_cm_given + m3;
    m4 = std::sqrt(M_cm * M_cm + m3 * m3 - 2 * M_cm * E3_cm);
  }
  double p3_tot_cm = std::sqrt(E3_cm * E3_cm - m3 * m3);
  double cond = (M_cm / p_tot_cm) * (p3_tot_cm / m3);
  double insq = (M_cm * M_cm * p3_tot_cm * p3_tot_cm -
                 m3 * m3 * p_tot_cm * p_tot_cm * sin_lab_sq);
  double p3_tot_1 = 0;
  double p3_tot_2 = 0;
  double E3k_1 = 0;
  double E3k_2 = 0;
  Direction p3_1 = {0, 0, 0};
  Direction p3_2 = {0, 0, 0};
  double Fp3cm_1[4];
  double Fp3cm_2[4];
  const auto& rx {nuc->reactions_[0]};
  // auto& d = rx->products_[0].distribution_[0];
  // auto d_ = dynamic_cast<UncorrelatedAngleEnergy*>(d.get());

  double q = (p_tot_cm / E_cm) * (E3_cm / p3_tot_cm);
  double approx_tol = 0.0001;

  if (insq >= 0) //( (cond > 1) || ( (cond < 1) && (theta < std::asin(cond)) ) )
  {
    // first solution

    p3_tot_1 = ((M_cm * M_cm + m3 * m3 - m4 * m4) * p_tot_cm * cos_lab +
                 2 * E_cm * std::sqrt(insq)) /
               2 / (M_cm * M_cm + p_tot_cm * p_tot_cm * sin_lab_sq);
    if (p3_tot_1 <= 0)
      return;
    p3_1 = u_lab_unit * p3_tot_1;
    double E3_tot_1 = std::sqrt(p3_tot_1 * p3_tot_1 + m3 * m3);
    E3k_1 = (E3_tot_1 - m3) * 1e6; // back to eV
    double B1[4] = {E3_tot_1, p3_1.x, p3_1.y, p3_1.z};
    boostf(A, B1, Fp3cm_1);

    double p3cm_tot_1 =
      std::sqrt(Fp3cm_1[1] * Fp3cm_1[1] + Fp3cm_1[2] * Fp3cm_1[2] +
                Fp3cm_1[3] * Fp3cm_1[3]);
    double mucm_1 =
      (Fp3cm_1[1] * p1_cm[1] + Fp3cm_1[2] * p1_cm[2] + Fp3cm_1[3] * p1_cm[3]) /
      (p1_tot_cm * p3cm_tot_1); // good until here
    if (std::abs(mucm_1) > 1.0) {
      mucm_1 = std::copysign(1.0, mucm_1);
    }
    // double pdf1cm = d_->angle().get_pdf(p.E_last(),mucm_1,p.current_seed());
    // pdfs_cm.push_back(pdf1cm);
    mu_cm.push_back(mucm_1);

    double mucm03_1 =
      (Fp3cm_1[1] * p_cm.x + Fp3cm_1[2] * p_cm.y + Fp3cm_1[3] * p_cm.z) /
      (p_tot_cm * p3cm_tot_1);

    if (std::abs(mucm03_1) > 1.0) {
      mucm03_1 = std::copysign(1.0, mucm03_1);
    }
    double sincm1 = std::sqrt(
      1 - mucm03_1 * mucm03_1); // if this is zero derivative is inf so pdf is 0
    double sin_ratio1 = std::sqrt(sin_lab_sq) / sincm1;
    double derivative1 =
      gamma * (1 + q * mucm03_1) * (sin_ratio1 * sin_ratio1 * sin_ratio1);
    if (sin_lab_sq < approx_tol) {
      derivative1 = ((cos_lab) / (gamma * mucm03_1 * (1 + q * mucm03_1))) *
                    ((cos_lab) / (gamma * mucm03_1 * (1 + q * mucm03_1)));
    }
    Js.push_back(derivative1);

    Particle ghost_particle = Particle();
    ghost_particle.initialize_ghost_particle(p, u_lab_unit, E3k_1);
    ghost_particles.push_back(ghost_particle);

    if (true) //((cond < 1) && (theta < std::asin(cond)))
    {
      // second solution

      p3_tot_2 = ((M_cm * M_cm + m3 * m3 - m4 * m4) * p_tot_cm * cos_lab -
                   2 * E_cm * std::sqrt(insq)) /
                 2 / (M_cm * M_cm + p_tot_cm * p_tot_cm * sin_lab_sq);
      if (p3_tot_2 < 0)
        return;
      p3_2 = u_lab_unit * p3_tot_2;
      double E3_tot_2 = std::sqrt(p3_tot_2 * p3_tot_2 + m3 * m3);
      E3k_2 = (E3_tot_2 - m3) * 1e6;
      if (p3_tot_2 < 0 || E3k_2 < 0)
        return;
      double B2[4] = {E3_tot_2, p3_2.x, p3_2.y, p3_2.z};
      boostf(A, B2, Fp3cm_2);
      double p3cm_tot_2 =
        std::sqrt(Fp3cm_2[1] * Fp3cm_2[1] + Fp3cm_2[2] * Fp3cm_2[2] +
                  Fp3cm_2[3] * Fp3cm_2[3]);
      double mucm_2 = (Fp3cm_2[1] * p1_cm[1] + Fp3cm_2[2] * p1_cm[2] +
                        Fp3cm_2[3] * p1_cm[3]) /
                      (p1_tot_cm * p3cm_tot_2);
      if (std::abs(mucm_2) > 1) {
        mucm_2 = std::copysign(1.0, mucm_2);
      }
      // double pdf2cm =
      // d_->angle().get_pdf(p.E_last(),mucm_2,p.current_seed());
      // pdfs_cm.push_back(pdf2cm);
      mu_cm.push_back(mucm_2);

      double mucm03_2 =
        (Fp3cm_2[1] * p_cm.x + Fp3cm_2[2] * p_cm.y + Fp3cm_2[3] * p_cm.z) /
        (p_tot_cm * p3cm_tot_1);
      if (std::abs(mucm03_2) > 1.0) {
        mucm03_2 = std::copysign(1.0, mucm03_2);
      }
      double sincm2 = std::sqrt(1 - mucm03_2 * mucm03_2);
      double sin_ratio2 = std::sqrt(sin_lab_sq) / sincm2;
      double derivative2 =
        gamma * (1 + q * mucm03_2) * (sin_ratio2 * sin_ratio2 * sin_ratio2);
      if (sin_lab_sq < approx_tol) {
        derivative2 = ((cos_lab) / (gamma * mucm03_2 * (1 + q * mucm03_2))) *
                      ((cos_lab) / (gamma * mucm03_2 * (1 + q * mucm03_2)));
      }
      Js.push_back(derivative2);

      Particle ghost_particle = Particle();
      ghost_particle.initialize_ghost_particle(p, u_lab_unit, E3k_2);
      ghost_particles.push_back(ghost_particle);
    }
  }
}
void score_ghost_particle(Particle& ghost_p, double pdf_lab, int i_tally)
{

  double myflux;
  double det_pos[4];
  get_det_pos(det_pos, i_tally);
  Direction u_lab {det_pos[0] - ghost_p.r().x, // towards the detector
    det_pos[1] - ghost_p.r().y, det_pos[2] - ghost_p.r().z};
  Direction u_lab_unit = u_lab / u_lab.norm(); // normalize
  double total_distance = u_lab.norm();
  double R0 = det_pos[3];
  double total_MFP1 = get_MFP(ghost_p, total_distance);

  if (total_distance < R0) {
    if (ghost_p.macro_xs().total == 0) {
      myflux = (ghost_p.wgt() * pdf_lab) / (2 / 3 * PI * R0 * R0);

    } else {
      // mutliplying in 1=10000/10000 to avoid float point error
      myflux = (10000 * ghost_p.wgt() * pdf_lab *
                 (1 - exp(-ghost_p.macro_xs().total * R0))) /
               (10000 * 2 / 3 * PI * R0 * R0 * R0 * ghost_p.macro_xs().total);
    }
  } else {
    myflux = (ghost_p.wgt()) * exp(-total_MFP1) /
             (2 * PI * total_distance * total_distance) * pdf_lab;
  }

  if (std::isnan(myflux) || std::abs(myflux) > 1e10) {

    // Cause a segmentation fault to crash the program
    int* ptr = nullptr;
    *ptr = 42; // This will trigger a segmentation fault
  }

  if (myflux < 0) {
    fatal_error("negetive flux");
  }
  if (ghost_p.type() != ParticleType::neutron) {
    myflux = 0;
  }
  const Tally& tally {*model::tallies[i_tally]};
  // Initialize an iterator over valid filter bin combinations.  If there are
  // no valid combinations, use a continue statement to ensure we skip the
  // assume_separate break below.
  auto filter_iter = FilterBinIter(tally, ghost_p);
  auto end = FilterBinIter(tally, true, &ghost_p.filter_matches());
  if (filter_iter == end)
    return;

  // Loop over filter bins.

  for (; filter_iter != end; ++filter_iter) {
    auto filter_index = filter_iter.index_;
    auto filter_weight = filter_iter.weight_;

    // Loop over nuclide bins.
    for (auto i = 0; i < tally.nuclides_.size(); ++i) {
      auto i_nuclide = tally.nuclides_[i];

      double atom_density = 0.;
      if (i_nuclide >= 0) {
        auto j =
          model::materials[ghost_p.material()]->mat_nuclide_index_[i_nuclide];
        if (j == C_NONE)
          continue;
        atom_density = model::materials[ghost_p.material()]->atom_density_(j);
      }
      // TODO: consider replacing this "if" with pointers or templates
      if (settings::run_CE) {
        score_general_ce_nonanalog(ghost_p, i_tally, i * tally.scores_.size(),
          filter_index, filter_weight, i_nuclide, atom_density, myflux);
      } else {
        fatal_error("multi group not implemnted for point tally");
      }
    }
  }

  // If the user has specified that we can assume all tallies are spatially
  // separate, this implies that once a tally has been scored to, we needn't
  // check the others. This cuts down on overhead when there are many
  // tallies specified
  if (settings::assume_separate)
    return;

  // Reset all the filter matches for the next tally event.
  for (auto& match : ghost_p.filter_matches())
    match.bins_present_ = false;
}

void score_pulse_height_tally(Particle& p, const vector<int>& tallies)
{
  // The pulse height tally in OpenMC hijacks the logic of CellFilter and
  // EnergyFilter to score specific quantities related to particle pulse height.
  // This is achieved by setting the pulse-height cell of the tally to the cell
  // of the particle being scored, and the energy to the particle's last
  // recorded energy (E_last()). After the tally is scored, the values are reset
  // to ensure proper accounting and avoid interference with subsequent
  // calculations or tallies.

  // Save original cell/energy information
  int orig_n_coord = p.n_coord();
  int orig_cell = p.coord(0).cell();
  double orig_E_last = p.E_last();

  for (auto i_tally : tallies) {
    auto& tally {*model::tallies[i_tally]};

    // Determine all CellFilter in the tally
    for (const auto& filter : tally.filters()) {
      auto cell_filter =
        dynamic_cast<CellFilter*>(model::tally_filters[filter].get());
      if (cell_filter != nullptr) {

        const auto& cells = cell_filter->cells();
        // Loop over all cells in the CellFilter
        for (auto cell_index = 0; cell_index < cells.size(); ++cell_index) {
          int cell_id = cells[cell_index];

          // Temporarily change cell of particle
          p.n_coord() = 1;
          p.coord(0).cell() = cell_id;

          // Determine index of cell in model::pulse_height_cells
          auto it = std::find(model::pulse_height_cells.begin(),
            model::pulse_height_cells.end(), cell_id);
          int index = std::distance(model::pulse_height_cells.begin(), it);

          // Temporarily change energy of particle to pulse-height value
          p.E_last() = p.pht_storage()[index];

          // Initialize an iterator over valid filter bin combinations. If
          // there are no valid combinations, use a continue statement to ensure
          // we skip the assume_separate break below.
          auto filter_iter = FilterBinIter(tally, p);
          auto end = FilterBinIter(tally, true, &p.filter_matches());
          if (filter_iter == end)
            continue;

          // Loop over filter bins.
          for (; filter_iter != end; ++filter_iter) {
            auto filter_index = filter_iter.index_;
            auto filter_weight = filter_iter.weight_;

            // Loop over scores.
            for (auto score_index = 0; score_index < tally.scores_.size();
                 ++score_index) {
#pragma omp atomic
              tally.results_(filter_index, score_index, TallyResult::VALUE) +=
                filter_weight;
            }
          }

          // Reset all the filter matches for the next tally event.
          for (auto& match : p.filter_matches())
            match.bins_present_ = false;
        }
      }
    }
    // Restore cell/energy
    p.n_coord() = orig_n_coord;
    p.coord(0).cell() = orig_cell;
    p.E_last() = orig_E_last;
  }
}
//
} // namespace openmc
