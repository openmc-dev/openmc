#include "openmc/tallies/tally_scoring.h"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/ifp.h"
#include "openmc/material.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/reaction_product.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/string_utils.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/sensitivity.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell.h"
#include "openmc/tallies/filter_delayedgroup.h"
#include "openmc/tallies/filter_energy.h"

#include <string>

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

void score_tracklength_tally_general(
  Particle& p, double flux, const vector<int>& tallies)
{
  // Set 'none' value for log union grid index
  int i_log_union = C_NONE;

  for (auto i_tally : tallies) {
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

void score_timed_tracklength_tally(Particle& p, double total_distance)
{
  double speed = p.speed();
  double total_dt = total_distance / speed;

  // save particle last state
  auto time_last = p.time_last();
  auto r_last = p.r_last();

  // move particle back
  p.move_distance(-total_distance);
  p.time() -= total_dt;
  p.lifetime() -= total_dt;

  double distance_traveled = 0.0;
  while (distance_traveled < total_distance) {

    double distance = std::min(distance_to_time_boundary(p.time(), speed),
      total_distance - distance_traveled);
    double dt = distance / speed;

    // Save particle last state for tracklength tallies
    p.time_last() = p.time();
    p.r_last() = p.r();

    // Advance particle in space and time
    p.move_distance(distance);
    p.time() += dt;
    p.lifetime() += dt;

    // Determine the tracklength estimate of the flux
    double flux = p.wgt() * distance;

    score_tracklength_tally_general(
      p, flux, model::active_timed_tracklength_tallies);
    distance_traveled += distance;
  }

  p.time_last() = time_last;
  p.r_last() = r_last;
}

void score_tracklength_tally(Particle& p, double distance)
{

  // Determine the tracklength estimate of the flux
  double flux = p.wgt() * distance;

  score_tracklength_tally_general(p, flux, model::active_tracklength_tallies);
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

        // Add a check to see if tally is a sensitivity tally, if so use special function
        // check if sens_ != C_NONE
        if (tally.sens_ != C_NONE) {
          score_collision_sensitivity_tally(p, i_tally, i*tally.scores_.size(), filter_index,
            filter_weight, i_nuclide, atom_density, flux);
          continue;
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

void score_collision_sensitivity_tally(Particle& p, int i_tally, int start_index, int filter_index,
  double filter_weight, int i_nuclide, double atom_density, double flux)
{
  SensitivityTally& tally = dynamic_cast<SensitivityTally&>(*model::tallies[i_tally]);
  
  // Add check to see if filter is an importance filter, if not return

  // Get the pre-collision energy of the particle.
  auto E = p.E_last();
  
  // Determine how much weight was absorbed due to survival biasing 
  double wgt_absorb = 0.0;  
  if (settings::survival_biasing) {           
      wgt_absorb = p.wgt_last() * p.neutron_xs(i_nuclide).absorption /
                              p.neutron_xs(i_nuclide).total;
      flux = (p.wgt_last() - wgt_absorb) / p.macro_xs().total;                             
  }
  
  for (auto i = 0; i < tally.scores_.size(); ++i) {
    auto score_bin = tally.scores_[i];
    auto score_index = start_index + i;
    double score = 0.0;
        
    switch (score_bin) {
    case SCORE_FLUX:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // All events score to a flux bin. We actually use a collision estimator
        // in place of an analog one since there is no way to count 'events'
        // exactly for the flux
        if (settings::survival_biasing) {
          // We need to account for the fact that some weight was already
          // absorbed
          score = p.wgt_last() - wgt_absorb;
        } else {
          score = p.wgt_last();
        }

        if (p.type() == ParticleType::neutron ||
          p.type() == ParticleType::photon) {
          score *= flux / p.macro_xs().total;
        } else {
          score = 0.;
        }
      } else {
        score = flux;
      }
      break;


    case SCORE_TOTAL:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // All events will score to the total reaction rate. We can just use
        // use the weight of the particle entering the collision as the score
        if (settings::survival_biasing) {
          // We need to account for the fact that some weight was already
          // absorbed
          score = (p.wgt_last() - wgt_absorb) * flux;
        } else {
          score = p.wgt_last() * flux;
        }

      } else {
        if (i_nuclide >= 0) {
          if (p.type() == ParticleType::neutron) {
            score = p.neutron_xs(i_nuclide).total * atom_density * flux;
          } else if (p.type() == ParticleType::photon) {
            score = p.photon_xs(i_nuclide).total * atom_density * flux;
          }
        } else {
          score = p.macro_xs().total * flux;
        }
      }
      break;


    case SCORE_INVERSE_VELOCITY:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // All events score to an inverse velocity bin. We actually use a
        // collision estimator in place of an analog one since there is no way
        // to count 'events' exactly for the inverse velocity
        if (settings::survival_biasing) {
          // We need to account for the fact that some weight was already
          // absorbed
          score = p.wgt_last() - wgt_absorb;
        } else {
          score = p.wgt_last();
        }
        score *= flux / p.macro_xs().total;
      } else {
        score = flux;
      }
      // Score inverse velocity in units of s/cm.
      score /= std::sqrt(2. * E / MASS_NEUTRON_EV) * C_LIGHT * 100.;
      break;


    case SCORE_SCATTER:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // Skip any event where the particle didn't scatter
        if (p.event() != TallyEvent::SCATTER) continue;
        // Since only scattering events make it here, again we can use the
        // weight entering the collision as the estimator for the reaction rate
        score = p.wgt_last() * flux;
      } else {
        if (i_nuclide >= 0) {
          score = (p.neutron_xs(i_nuclide).total
            - p.neutron_xs(i_nuclide).absorption) * atom_density * flux;
        } else {
          score = (p.macro_xs().total
            - p.macro_xs().absorption) * flux;
        }
      }
      break;


    case SCORE_NU_SCATTER:
      // Only analog estimators are available.
      // Skip any event where the particle didn't scatter
      if (p.event() != TallyEvent::SCATTER) continue;
      // For scattering production, we need to use the pre-collision weight
      // times the yield as the estimate for the number of neutrons exiting a
      // reaction with neutrons in the exit channel
      if (p.event_mt() == ELASTIC || p.event_mt() == N_LEVEL ||
        (p.event_mt() >= N_N1 && p.event_mt() <= N_NC)) {
        // Don't waste time on very common reactions we know have
        // multiplicities of one.
        score = p.wgt_last() * flux;
      } else {
        // Get yield and apply to score
        auto m = data::nuclides[p.event_nuclide()]->reaction_index_[p.event_mt()];
        const auto& rxn {*data::nuclides[p.event_nuclide()]->reactions_[m]};
        score = p.wgt_last() * flux * (*rxn.products_[0].yield_)(E);
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
          if (p.event() == TallyEvent::SCATTER) continue;
          // All fission and absorption events will contribute here, so we
          // can just use the particle's weight entering the collision
          score = p.wgt_last() * flux;
        }
      } else {
        if (i_nuclide >= 0) {
          score = p.neutron_xs(i_nuclide).absorption * atom_density
            * flux;
        } else {
          score = p.macro_xs().absorption * flux;
        }
      }
      break;


    case SCORE_FISSION:
      if (p.macro_xs().absorption == 0) continue;
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // fission
          if (p.neutron_xs(p.event_nuclide()).absorption > 0) {
            score = wgt_absorb
              * p.neutron_xs(p.event_nuclide()).fission
              / p.neutron_xs(p.event_nuclide()).absorption * flux;
          } else {
            score = 0.;
          }
        } else {
          // Skip any non-absorption events
          if (p.event() == TallyEvent::SCATTER) continue;
          // All fission events will contribute, so again we can use particle's
          // weight entering the collision as the estimate for the fission
          // reaction rate
          score = p.wgt_last()
            * p.neutron_xs(p.event_nuclide()).fission
            / p.neutron_xs(p.event_nuclide()).absorption * flux;
        }
      } else {
        if (i_nuclide >= 0) {
          score = p.neutron_xs(i_nuclide).fission * atom_density * flux;
        } else {
          score = p.macro_xs().fission * flux;
        }
      }
      break;


    case SCORE_NU_FISSION:
      if (p.macro_xs().absorption == 0) continue;
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
          if (p.neutron_xs(p.event_nuclide()).absorption > 0) {
            score = wgt_absorb
              * p.neutron_xs(p.event_nuclide()).nu_fission
              / p.neutron_xs(p.event_nuclide()).absorption * flux;
          } else {
            score = 0.;
          }
        } else {
          // Skip any non-fission events
          if (!p.fission()) continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score.
          score = simulation::keff * p.wgt_bank() * flux;
        }
      } else {
        if (i_nuclide >= 0) {
          score = p.neutron_xs(i_nuclide).nu_fission * atom_density
            * flux;
        } else {
          score = p.macro_xs().nu_fission * flux;
        }
      }
      break;


    case SCORE_PROMPT_NU_FISSION:
      if (p.macro_xs().absorption == 0) continue;
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
          if (p.neutron_xs(p.event_nuclide()).absorption > 0) {
            score = wgt_absorb
              * p.neutron_xs(p.event_nuclide()).fission
              * data::nuclides[p.event_nuclide()]
              ->nu(E, ReactionProduct::EmissionMode::prompt)
              / p.neutron_xs(p.event_nuclide()).absorption * flux;
          } else {
            score = 0.;
          }
        } else {
          // Skip any non-fission events
          if (!p.fission()) continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score.
          auto n_delayed = std::accumulate(p.n_delayed_bank(),
            p.n_delayed_bank()+MAX_DELAYED_GROUPS, 0);
          auto prompt_frac = 1. - n_delayed / static_cast<double>(p.n_bank());
          score = simulation::keff * p.wgt_bank() * prompt_frac * flux;
        }
      } else {
        if (i_nuclide >= 0) {
          score = p.neutron_xs(i_nuclide).fission
            * data::nuclides[i_nuclide]
            ->nu(E, ReactionProduct::EmissionMode::prompt)
            * atom_density * flux;
        } else {
          score = 0.;
          // Add up contributions from each nuclide in the material.
          if (p.material() != MATERIAL_VOID) {
            const Material& material {*model::materials[p.material()]};
            for (auto i = 0; i < material.nuclide_.size(); ++i) {
              auto j_nuclide = material.nuclide_[i];
              auto atom_density = material.atom_density_(i);
              score += p.neutron_xs(j_nuclide).fission
                * data::nuclides[j_nuclide]
                ->nu(E, ReactionProduct::EmissionMode::prompt)
                * atom_density * flux;
            }
          }
        }
      }
      break;


    case SCORE_DELAYED_NU_FISSION:
      if (p.macro_xs().absorption == 0) continue;
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
          if (p.neutron_xs(p.event_nuclide()).absorption > 0
            && data::nuclides[p.event_nuclide()]->fissionable_) {
            if (tally.delayedgroup_filter_ != C_NONE) {
              auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
              const DelayedGroupFilter& filt
                {*dynamic_cast<DelayedGroupFilter*>(
                model::tally_filters[i_dg_filt].get())};
              // Tally each delayed group bin individually
              for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                auto dg = filt.groups()[d_bin];
                auto yield = data::nuclides[p.event_nuclide()]
                  ->nu(E, ReactionProduct::EmissionMode::delayed, dg);
                score = wgt_absorb * yield
                  * p.neutron_xs(p.event_nuclide()).fission
                  / p.neutron_xs(p.event_nuclide()).absorption * flux;
                score_fission_delayed_dg(i_tally, d_bin, score,
                  score_index, p.filter_matches());
              }
              continue;
            } else {
              // If the delayed group filter is not present, compute the score
              // by multiplying the absorbed weight by the fraction of the
              // delayed-nu-fission xs to the absorption xs
              score = wgt_absorb
                * p.neutron_xs(p.event_nuclide()).fission
                * data::nuclides[p.event_nuclide()]
                ->nu(E, ReactionProduct::EmissionMode::delayed)
                / p.neutron_xs(p.event_nuclide()).absorption *flux;
            }
          }
        } else {
          // Skip any non-fission events
          if (!p.fission()) continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score. Loop over the neutrons produced from fission and check which
          // ones are delayed. If a delayed neutron is encountered, add its
          // contribution to the fission bank to the score.
          if (tally.delayedgroup_filter_ != C_NONE) {
            auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
            const DelayedGroupFilter& filt
              {*dynamic_cast<DelayedGroupFilter*>(
              model::tally_filters[i_dg_filt].get())};
            // Tally each delayed group bin individually
            for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
              auto d = filt.groups()[d_bin];
              score = simulation::keff * p.wgt_bank() / p.n_bank()
                * p.n_delayed_bank(d-1) * flux;
              score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches());
            }
            continue;
          } else {
            // Add the contribution from all delayed groups
            auto n_delayed = std::accumulate(p.n_delayed_bank(),
              p.n_delayed_bank()+MAX_DELAYED_GROUPS, 0);
            score = simulation::keff * p.wgt_bank() / p.n_bank() * n_delayed
              * flux;
          }
        }
      } else {
        if (i_nuclide >= 0) {
          if (tally.delayedgroup_filter_ != C_NONE) {
            auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
            const DelayedGroupFilter& filt
              {*dynamic_cast<DelayedGroupFilter*>(
              model::tally_filters[i_dg_filt].get())};
            // Tally each delayed group bin individually
            for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
              auto d = filt.groups()[d_bin];
              auto yield = data::nuclides[i_nuclide]
                ->nu(E, ReactionProduct::EmissionMode::delayed, d);
              score = p.neutron_xs(i_nuclide).fission * yield
                * atom_density * flux;
              score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches());
            }
            continue;
          } else {
            // If the delayed group filter is not present, compute the score
            // by multiplying the delayed-nu-fission macro xs by the flux
            score = p.neutron_xs(i_nuclide).fission
              * data::nuclides[i_nuclide]
              ->nu(E, ReactionProduct::EmissionMode::delayed)
              * atom_density * flux;
          }
        } else {
          if (tally.delayedgroup_filter_ != C_NONE) {
            auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
            const DelayedGroupFilter& filt
              {*dynamic_cast<DelayedGroupFilter*>(
              model::tally_filters[i_dg_filt].get())};
            if (p.material() != MATERIAL_VOID) {
              const Material& material {*model::materials[p.material()]};
              for (auto i = 0; i < material.nuclide_.size(); ++i) {
                auto j_nuclide = material.nuclide_[i];
                auto atom_density = material.atom_density_(i);
                // Tally each delayed group bin individually
                for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                  auto d = filt.groups()[d_bin];
                  auto yield = data::nuclides[j_nuclide]
                    ->nu(E, ReactionProduct::EmissionMode::delayed, d);
                  score = p.neutron_xs(j_nuclide).fission * yield
                    * atom_density * flux;
                  score_fission_delayed_dg(i_tally, d_bin, score,
                    score_index, p.filter_matches());
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
                score += p.neutron_xs(j_nuclide).fission
                  * data::nuclides[j_nuclide]
                  ->nu(E, ReactionProduct::EmissionMode::delayed)
                  * atom_density * flux;
              }
            }
          }
        }
      }
      break;


    case SCORE_DECAY_RATE:
      if (p.macro_xs().absorption == 0) continue;
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // delayed-nu-fission
          const auto& nuc {*data::nuclides[p.event_nuclide()]};
          if (p.neutron_xs(p.event_nuclide()).absorption > 0
            && nuc.fissionable_) {
            const auto& rxn {*nuc.fission_rx_[0]};
            if (tally.delayedgroup_filter_ != C_NONE) {
              auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
              const DelayedGroupFilter& filt
                {*dynamic_cast<DelayedGroupFilter*>(
                model::tally_filters[i_dg_filt].get())};
              // Tally each delayed group bin individually
              for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                auto d = filt.groups()[d_bin];
                auto yield
                  = nuc.nu(E, ReactionProduct::EmissionMode::delayed, d);
                auto rate = rxn.products_[d].decay_rate_;
                score = wgt_absorb * yield
                  * p.neutron_xs(p.event_nuclide()).fission
                  / p.neutron_xs(p.event_nuclide()).absorption
                  * rate * flux;
                score_fission_delayed_dg(i_tally, d_bin, score,
                  score_index, p.filter_matches());
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
              for (auto d = 0; d < rxn.products_.size() - 2; ++d) {
                auto yield
                  = nuc.nu(E, ReactionProduct::EmissionMode::delayed, d+1);
                auto rate = rxn.products_[d+1].decay_rate_;
                score += rate * wgt_absorb
                  * p.neutron_xs(p.event_nuclide()).fission * yield
                  / p.neutron_xs(p.event_nuclide()).absorption * flux;
              }
            }
          }
        } else {
          // Skip any non-fission events
          if (!p.fission()) continue;
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
                const DelayedGroupFilter& filt
                  {*dynamic_cast<DelayedGroupFilter*>(
                  model::tally_filters[i_dg_filt].get())};
                // Find the corresponding filter bin and then score
                for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                  auto d = filt.groups()[d_bin];
                  if (d == g)
                    score_fission_delayed_dg(i_tally, d_bin, score,
                      score_index, p.filter_matches());
                }
                score = 0.;
              }
            }
          }
        }
      } else {
        if (i_nuclide >= 0) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          if (!nuc.fissionable_) continue;
          const auto& rxn {*nuc.fission_rx_[0]};
          if (tally.delayedgroup_filter_ != C_NONE) {
            auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
            const DelayedGroupFilter& filt
              {*dynamic_cast<DelayedGroupFilter*>(
              model::tally_filters[i_dg_filt].get())};
            // Tally each delayed group bin individually
            for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
              auto d = filt.groups()[d_bin];
              auto yield
                = nuc.nu(E, ReactionProduct::EmissionMode::delayed, d);
              auto rate = rxn.products_[d].decay_rate_;
              score = p.neutron_xs(i_nuclide).fission * yield * flux
                * atom_density * rate;
              score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches());
            }
            continue;
          } else {
            score = 0.;
            // We need to be careful not to overshoot the number of
            // delayed groups since this could cause the range of the
            // rxn.products_ array to be exceeded. Hence, we use the size
            // of this array and not the MAX_DELAYED_GROUPS constant for
            // this loop.
            for (auto d = 0; d < rxn.products_.size() - 2; ++d) {
              auto yield
                = nuc.nu(E, ReactionProduct::EmissionMode::delayed, d+1);
              auto rate = rxn.products_[d+1].decay_rate_;
              score += p.neutron_xs(i_nuclide).fission * flux
                * yield * atom_density * rate;
            }
          }
        } else {
          if (tally.delayedgroup_filter_ != C_NONE) {
            auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
            const DelayedGroupFilter& filt
              {*dynamic_cast<DelayedGroupFilter*>(
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
                    auto yield
                      = nuc.nu(E, ReactionProduct::EmissionMode::delayed, d);
                    auto rate = rxn.products_[d].decay_rate_;
                    score = p.neutron_xs(j_nuclide).fission * yield
                      * flux * atom_density * rate;
                    score_fission_delayed_dg(i_tally, d_bin, score,
                      score_index, p.filter_matches());
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
                  for (auto d = 0; d < rxn.products_.size() - 2; ++d) {
                    auto yield
                      = nuc.nu(E, ReactionProduct::EmissionMode::delayed, d+1);
                    auto rate = rxn.products_[d+1].decay_rate_;
                    score += p.neutron_xs(j_nuclide).fission
                      * yield * atom_density * flux * rate;
                  }
                }
              }
            }
          }
        }
      }
      break;


    case SCORE_KAPPA_FISSION:
      if (p.macro_xs().absorption == 0.) continue;
      score = 0.;
      // Kappa-fission values are determined from the Q-value listed for the
      // fission cross section.
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // fission scaled by the Q-value
          const auto& nuc {*data::nuclides[p.event_nuclide()]};
          if (p.neutron_xs(p.event_nuclide()).absorption > 0
            && nuc.fissionable_) {
            const auto& rxn {*nuc.fission_rx_[0]};
            score = wgt_absorb * rxn.q_value_
              * p.neutron_xs(p.event_nuclide()).fission
              / p.neutron_xs(p.event_nuclide()).absorption * flux;
          }
        } else {
          // Skip any non-absorption events
          if (p.event() == TallyEvent::SCATTER) continue;
          // All fission events will contribute, so again we can use particle's
          // weight entering the collision as the estimate for the fission
          // reaction rate
          const auto& nuc {*data::nuclides[p.event_nuclide()]};
          if (p.neutron_xs(p.event_nuclide()).absorption > 0
            && nuc.fissionable_) {
            const auto& rxn {*nuc.fission_rx_[0]};
            score = p.wgt_last() * rxn.q_value_
              * p.neutron_xs(p.event_nuclide()).fission
              / p.neutron_xs(p.event_nuclide()).absorption * flux;
          }
        }
      } else {
        if (i_nuclide >= 0) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          if (nuc.fissionable_) {
            const auto& rxn {*nuc.fission_rx_[0]};
            score = rxn.q_value_ * p.neutron_xs(i_nuclide).fission
              * atom_density * flux;
          }
        } else if (p.material() != MATERIAL_VOID) {
          const Material& material {*model::materials[p.material()]};
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto j_nuclide = material.nuclide_[i];
            auto atom_density = material.atom_density_(i);
            const auto& nuc {*data::nuclides[j_nuclide]};
            if (nuc.fissionable_) {
              const auto& rxn {*nuc.fission_rx_[0]};
              score += rxn.q_value_ * p.neutron_xs(j_nuclide).fission
                * atom_density * flux;
            }
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
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // Check if event MT matches
        if (p.event_mt() != ELASTIC) continue;
        score = p.wgt_last() * flux;
      } else {
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
              score += p.neutron_xs(j_nuclide).elastic * atom_density
                * flux;
            }
          }
        }
      }
      break;

    case SCORE_FISS_Q_PROMPT:
    case SCORE_FISS_Q_RECOV:
      if (p.macro_xs().absorption == 0.) continue;
      score = score_fission_q(p, score_bin, tally, flux, i_nuclide, atom_density);
      break;


    case N_2N:
    case N_3N:
    case N_4N:
    case N_GAMMA:
    case N_P:
    case N_A:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // Check if the event MT matches
        if (p.event_mt() != score_bin) continue;
        score = p.wgt_last() * flux;
      } else {
        int m;
        switch (score_bin) {
        case N_GAMMA: m = 0; break;
        case N_P: m = 1; break;
        case N_A: m = 2; break;
        case N_2N: m = 3; break;
        case N_3N: m = 4; break;
        case N_4N: m = 5; break;
        }
        if (i_nuclide >= 0) {
          score = p.neutron_xs(i_nuclide).reaction[m] * atom_density
            * flux;
        } else {
          score = 0.;
          if (p.material() != MATERIAL_VOID) {
            const Material& material {*model::materials[p.material()]};
            for (auto i = 0; i < material.nuclide_.size(); ++i) {
              auto j_nuclide = material.nuclide_[i];
              auto atom_density = material.atom_density_(i);
              score += p.neutron_xs(j_nuclide).reaction[m]
                * atom_density * flux;
            }
          }
        }
      }
      break;


    case HEATING:
      if (p.type() == ParticleType::neutron) {
        score = score_neutron_heating(p, tally, flux, HEATING,
            i_nuclide, atom_density);
      } else {
        // The energy deposited is the difference between the pre-collision and
        // post-collision energy...
        score = E - p.E();

        // ...less the energy of any secondary particles since they will be
        // transported individually later
        const auto& bank = p.secondary_bank();
        for (auto it = bank.end() - p.bank_second_E(); it < bank.end(); ++it) {
          score -= it->E;
        }

        score *= p.wgt_last();
      }
      break;

    default:

      // The default block is really only meant for redundant neutron reactions
      // (e.g. 444, 901)
      if (p.type() != ParticleType::neutron) continue;

      if (tally.estimator_ == TallyEstimator::ANALOG) {

        // Any other score is assumed to be a MT number. Thus, we just need
        // to check if it matches the MT number of the event
        if (p.event_mt() != score_bin) continue;
        score = p.wgt_last()*flux;
      } else {
        // Any other cross section has to be calculated on-the-fly
        if (score_bin < 2) fatal_error("Invalid score type on tally "
          + std::to_string(tally.id_));
        score = 0.;
        if (i_nuclide >= 0) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          auto m = nuc.reaction_index_[score_bin];
          if (m == C_NONE) continue;
          const auto& rxn {*nuc.reactions_[m]};
          auto i_temp = p.neutron_xs(i_nuclide).index_temp;
          if (i_temp >= 0) { // Can be false due to multipole
            auto i_grid = p.neutron_xs(i_nuclide).index_grid;
            auto f = p.neutron_xs(i_nuclide).interp_factor;
            const auto& xs {rxn.xs_[i_temp]};
            if (i_grid >= xs.threshold) {
              score = ((1.0 - f) * xs.value[i_grid-xs.threshold]
                + f * xs.value[i_grid-xs.threshold+1]) * atom_density * flux;
            }
          }
        } else if (p.material() != MATERIAL_VOID) {
          const Material& material {*model::materials[p.material()]};
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto j_nuclide = material.nuclide_[i];
            auto atom_density = material.atom_density_(i);
            const auto& nuc {*data::nuclides[j_nuclide]};
            auto m = nuc.reaction_index_[score_bin];
            if (m == C_NONE) continue;
            const auto& rxn {*nuc.reactions_[m]};
            auto i_temp = p.neutron_xs(j_nuclide).index_temp;
            if (i_temp >= 0) { // Can be false due to multipole
              auto i_grid = p.neutron_xs(j_nuclide).index_grid;
              auto f = p.neutron_xs(j_nuclide).interp_factor;
              const auto& xs {rxn.xs_[i_temp]};
              if (i_grid >= xs.threshold) {
                score += ((1.0 - f) * xs.value[i_grid-xs.threshold]
                  + f * xs.value[i_grid-xs.threshold+1]) * atom_density
                  * flux;
              }
            }
          }
        }
      }
    }

    // Add sensitivity information on score for sensitivity tallies.
    // Retrieve particle's cumulative sensitivity
    
    if (score == 0.0) return;
    
    const auto& sens {model::tally_sens[tally.sens_]};
    const auto cumulative_sensitivities = p.cumulative_sensitivities(tally.sens_);
    
    if (settings::run_mode == RunMode::EIGENVALUE) {
      
      #pragma omp atomic
      tally.denominator_ += score*filter_weight;  //add an if statement if we want denominator? ;
          
      // Update tally results
      for (auto idx = 0; idx < cumulative_sensitivities.size(); idx++){
        #pragma omp atomic
        tally.results_(idx, score_index, SensitivityTallyResult::VALUE) += cumulative_sensitivities[idx]*score*filter_weight;
      }
      
      if (sens.sens_nuclide == p.fission_nuclide()){
        switch (sens.variable) {
        
        case SensitivityVariable::CROSS_SECTION:
        {
          if (sens.sens_reaction == SCORE_FISSION ){
            // Get the energy of the parent particle.
            auto E = p.E_parent();
            // Bin the energy.
            if (E >= sens.energy_bins_.front() && E <= sens.energy_bins_.back()) {
              auto bin = lower_bound_index(sens.energy_bins_.begin(), sens.energy_bins_.end(), E);
              #pragma omp atomic
              tally.previous_results_(bin, score_index, SensitivityTallyResult::VALUE) += score*filter_weight;
            }
          }
        }
          break;
        case SensitivityVariable::MULTIPOLE:
        {
          // check if in resonance range
          const auto& nuc {*data::nuclides[sens.sens_nuclide]};
          if (multipole_in_range(nuc, p.E_parent())){
            // Calculate derivative of the fission cross section at p->E_parent()
            double sig_s, sig_a, sig_f;
            std::tie(sig_s, sig_a, sig_f)
              = nuc.multipole_->evaluate(p.E_parent(), p.sqrtkT());
            auto derivative = nuc.multipole_->evaluate_pole_deriv_fission(p.E_parent(), p.sqrtkT()); //This actually needs to be parent kT
            // sum/bin 1/micro_sigma_scatter * derivative
            int start = derivative.first;
            int size  = derivative.second.size();
            double sen_score = score*filter_weight/sig_f;
            for (int deriv_idx = start; deriv_idx < start + size ; deriv_idx++){
              #pragma omp atomic
              tally.previous_results_(deriv_idx, score_index, SensitivityTallyResult::VALUE) += sen_score*derivative.second[deriv_idx - start];
            }
          }
          // multiply score by derivative of fission cross section wrt multipole parameters / fission cross section
          // at parent energy. so I need to also evaluate the fission derivative..        
        }
          break;
  
        case SensitivityVariable::CURVE_FIT:
        {
          // check if in resonance range
          const auto& nuc {*data::nuclides[sens.sens_nuclide]};
          if (multipole_in_range(nuc, p.E_parent())){
            // Calculate derivative of the fission cross section at p->E_parent_
            double sig_s, sig_a, sig_f;
            std::tie(sig_s, sig_a, sig_f)
              = nuc.multipole_->evaluate(p.E_parent(), p.sqrtkT());
            auto derivative = nuc.multipole_->evaluate_fit_deriv_fission(p.E_parent(), p.sqrtkT()); //This actually needs to be parent kT
            // sum/bin 1/micro_sigma_scatter * derivative
            int start = derivative.first;
            int size  = derivative.second.size();
            double sen_score = score*filter_weight/sig_f;
            for (int deriv_idx = start; deriv_idx < start + size ; deriv_idx++){
              #pragma omp atomic
              tally.previous_results_(deriv_idx, score_index, SensitivityTallyResult::VALUE) += sen_score*derivative.second[deriv_idx - start];
            }
          }
          // multiply score by derivative of fission cross section wrt multipole parameters / fission cross section
          // at parent energy. so I need to also evaluate the fission derivative..        
        }
          break;
        }
      }
    } else {      
            
      #pragma omp atomic
      tally.denominator_ += score*filter_weight;
      
      for (int idx = 0; idx < cumulative_sensitivities.size(); idx++){
        #pragma omp atomic
        tally.results_(idx, score_index, SensitivityTallyResult::VALUE) += cumulative_sensitivities[idx]*score*filter_weight;      
      } 
      
      // // direct effect???      
      // double atom_density = 0.;
      // if (sens.sens_nuclide >= 0) {
      //   auto j = model::materials[p.material()]->mat_nuclide_index_[sens.sens_nuclide];
      //   if (j == C_NONE) continue;
      //   atom_density = model::materials[p.material()]->atom_density_(j);
      // }
      // 
      // switch (sens.variable) {
      // 
      // case SensitivityVariable::CROSS_SECTION:
      // {
      //   // Calculate the sensitivity with respect to the cross section
      //   // at this energy
      // 
      //   // Get the post-collision energy of the particle.
      //   auto E = p.E();
      // 
      //   // Get the correct cross section
      //   double macro_xs;
      //   switch (sens.sens_reaction) {
      //   case SCORE_TOTAL:
      //     if (sens.sens_nuclide >=0){
      //         macro_xs = p.neutron_xs(sens.sens_nuclide).total * atom_density;
      //     } else {
      //         macro_xs = p.macro_xs().total;
      //     }
      //     break;
      //   case SCORE_SCATTER:
      //     if (sens.sens_nuclide >=0){
      //         macro_xs = (p.neutron_xs(sens.sens_nuclide).total 
      //         - p.neutron_xs(sens.sens_nuclide).absorption) * atom_density;
      //     } else {
      //         macro_xs = p.macro_xs().total - p.macro_xs().absorption;
      //     }
      //     break;
      //   case ELASTIC:
      //     if (sens.sens_nuclide >= 0) {
      //         if (p.neutron_xs(sens.sens_nuclide).elastic == CACHE_INVALID)
      //           data::nuclides[sens.sens_nuclide]->calculate_elastic_xs(p);
      //         macro_xs = p.neutron_xs(sens.sens_nuclide).elastic * atom_density;
      //       } 
      //     break;
      //   case SCORE_ABSORPTION: 
      //     if (sens.sens_nuclide >=0){
      //         macro_xs = p.neutron_xs(sens.sens_nuclide).absorption * atom_density;
      //     } else {
      //         macro_xs = p.macro_xs().absorption;
      //     }
      //     break;
      //   case SCORE_FISSION:
      //     if (p.macro_xs().absorption == 0) continue;
      // 
      //     if (sens.sens_nuclide >= 0) {
      //       macro_xs = p.neutron_xs(sens.sens_nuclide).fission * atom_density;
      //     } else {
      //       macro_xs = p.macro_xs().fission;
      //     }
      //     break;      
      //   default:
      //     if (sens.sens_nuclide >= 0) {
      //       macro_xs = get_nuclide_xs(p, sens.sens_nuclide, sens.sens_reaction) * atom_density;
      //     }
      //     break;
      //   }  
      //   // Bin the contribution.
      //   if (E >= sens.energy_bins_.front() && E <= sens.energy_bins_.back()) {
      //     auto bin = lower_bound_index(sens.energy_bins_.begin(), sens.energy_bins_.end(), E);
      //     tally.previous_results_(bin, score_index, SensitivityTallyResult::VALUE) += flux * macro_xs;
      //   }                
      // }
      // break;     
      // }
      //}      
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
} // namespace openmc
