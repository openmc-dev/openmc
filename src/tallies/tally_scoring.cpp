#include "openmc/tallies/tally_scoring.h"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
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
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_delayedgroup.h"
#include "openmc/tallies/filter_energy.h"

#include <string>

namespace openmc {

//==============================================================================
// FilterBinIter implementation
//==============================================================================

FilterBinIter::FilterBinIter(const Tally& tally, Particle& p)
  : filter_matches_{p.filter_matches_}, tally_{tally}
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

FilterBinIter::FilterBinIter(const Tally& tally, bool end,
    std::vector<FilterMatch>* particle_filter_matches)
  : filter_matches_{*particle_filter_matches}, tally_{tally}
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

FilterBinIter&
FilterBinIter::operator++()
{
  // Find the next valid combination of filter bins.  To do this, we search
  // backwards through the filters until we find the first filter whose bins
  // can be incremented.
  bool visited_all_combinations = true;
  for (int i = tally_.filters().size()-1; i >= 0; --i) {
    auto i_filt = tally_.filters(i);
    auto& match {filter_matches_[i_filt]};
    if (match.i_bin_ < match.bins_.size()-1) {
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

void
FilterBinIter::compute_index_weight()
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

void
score_fission_delayed_dg(int i_tally, int d_bin, double score, int score_index,
    std::vector<FilterMatch>& filter_matches)
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
  tally.results_(filter_index, score_index, TallyResult::VALUE) += score*filter_weight;

  // Reset the original delayed group bin
  dg_match.bins_[i_bin] = original_bin;
}

//! Helper function to retrieve fission q value from a nuclide

double get_nuc_fission_q(const Nuclide& nuc, const Particle& p, int score_bin)
{
  if (score_bin == SCORE_FISS_Q_PROMPT) {
    if (nuc.fission_q_prompt_) {
      return (*nuc.fission_q_prompt_)(p.E_last_);
    }
  } else if (score_bin == SCORE_FISS_Q_RECOV) {
    if (nuc.fission_q_recov_) {
      return (*nuc.fission_q_recov_)(p.E_last_);
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
    const Nuclide& nuc {*data::nuclides[p.event_nuclide_]};
    if (settings::survival_biasing) {
      // No fission events occur if survival biasing is on -- need to
      // calculate fraction of absorptions that would have resulted in
      // fission scaled by the Q-value
      if (p.neutron_xs_[p.event_nuclide_].absorption > 0) {
        return p.wgt_absorb_ * get_nuc_fission_q(nuc, p, score_bin)
               * p.neutron_xs_[p.event_nuclide_].fission * flux
               / p.neutron_xs_[p.event_nuclide_].absorption;
      }
    } else {
      // Skip any non-absorption events
      if (p.event_ == TallyEvent::SCATTER) return 0.0;
      // All fission events will contribute, so again we can use particle's
      // weight entering the collision as the estimate for the fission
      // reaction rate
      if (p.neutron_xs_[p.event_nuclide_].absorption > 0) {
        return p.wgt_last_ * get_nuc_fission_q(nuc, p, score_bin)
               * p.neutron_xs_[p.event_nuclide_].fission * flux
               / p.neutron_xs_[p.event_nuclide_].absorption;
      }
    }
  } else {
    if (i_nuclide >= 0) {
      const Nuclide& nuc {*data::nuclides[i_nuclide]};
      return get_nuc_fission_q(nuc, p, score_bin) * atom_density * flux
             * p.neutron_xs_[i_nuclide].fission;
    } else if (p.material_ != MATERIAL_VOID) {
      const Material& material {*model::materials[p.material_]};
      double score {0.0};
      for (auto i = 0; i < material.nuclide_.size(); ++i) {
        auto j_nuclide = material.nuclide_[i];
        auto atom_density = material.atom_density_(i);
        const Nuclide& nuc {*data::nuclides[j_nuclide]};
        score += get_nuc_fission_q(nuc, p, score_bin) * atom_density
          * p.neutron_xs_[j_nuclide].fission;
      }
      return score * flux;
    }
  }
  return 0.0;
}

//! Helper function to obtain the kerma coefficient for a given nuclide

double get_nuclide_neutron_heating(const Particle& p, const Nuclide& nuc,
    int rxn_index, int i_nuclide)
{
  size_t mt = nuc.reaction_index_[rxn_index];
  if (mt == C_NONE) return 0.0;

  auto i_temp = p.neutron_xs_[i_nuclide].index_temp;
  if (i_temp < 0) return 0.0; // Can be true due to multipole

  const auto& rxn {*nuc.reactions_[mt]};
  const auto& xs {rxn.xs_[i_temp]};
  auto i_grid = p.neutron_xs_[i_nuclide].index_grid;
  if (i_grid < xs.threshold) return 0.0;

  // Determine total kerma
  auto f = p.neutron_xs_[i_nuclide].interp_factor;
  double kerma = (1.0 - f) * xs.value[i_grid-xs.threshold]
    + f * xs.value[i_grid-xs.threshold+1];

  if (settings::run_mode == RunMode::EIGENVALUE) {
    // Determine kerma for fission as (EFR + EB)*sigma_f
    double kerma_fission = nuc.fragments_ ?
      ((*nuc.fragments_)(p.E_last_) + (*nuc.betas_)(p.E_last_))
      *p.neutron_xs_[i_nuclide].fission : 0.0;

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
  double score;
  // Get heating macroscopic "cross section"
  double heating_xs;
  if (i_nuclide >= 0) {
    const Nuclide& nuc {*data::nuclides[i_nuclide]};
    heating_xs = get_nuclide_neutron_heating(p, nuc, rxn_bin, i_nuclide);
    if (tally.estimator_ == TallyEstimator::ANALOG) {
      heating_xs /= p.neutron_xs_[i_nuclide].total;
    } else {
      heating_xs *= atom_density;
    }
  } else {
    if (p.material_ != MATERIAL_VOID) {
      heating_xs = 0.0;
      const Material& material {*model::materials[p.material_]};
      for (auto i = 0; i< material.nuclide_.size(); ++i) {
        int j_nuclide = material.nuclide_[i];
        double atom_density {material.atom_density_(i)};
        const Nuclide& nuc {*data::nuclides[j_nuclide]};
        heating_xs += atom_density * get_nuclide_neutron_heating(p, nuc, rxn_bin, j_nuclide);
      }
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        heating_xs /= p.macro_xs_.total;
      }
    }
  }
  score = heating_xs * flux;
  if (tally.estimator_ == TallyEstimator::ANALOG) {
    // All events score to a heating tally bin. We actually use a
    // collision estimator in place of an analog one since there is no
    // reaction-wise heating cross section
    if (settings::survival_biasing) {
      // Account for the fact that some weight has been absorbed
      score *= p.wgt_last_ + p.wgt_absorb_;
    } else {
      score *= p.wgt_last_;
    }
  }
  return score;
}

//! Helper function for nu-fission tallies with energyout filters.
//
//! In this case, we may need to score to multiple bins if there were multiple
//! neutrons produced with different energies.

void
score_fission_eout(Particle& p, int i_tally, int i_score, int score_bin)
{
  auto& tally {*model::tallies[i_tally]};
  auto i_eout_filt = tally.filters()[tally.energyout_filter_];
  auto i_bin = p.filter_matches_[i_eout_filt].i_bin_;
  auto bin_energyout = p.filter_matches_[i_eout_filt].bins_[i_bin];

  const EnergyoutFilter& eo_filt
    {*dynamic_cast<EnergyoutFilter*>(model::tally_filters[i_eout_filt].get())};

  // Note that the score below is weighted by keff. Since the creation of
  // fission sites is weighted such that it is expected to create n_particles
  // sites, we need to multiply the score by keff to get the true nu-fission
  // rate. Otherwise, the sum of all nu-fission rates would be ~1.0.

  // loop over number of particles banked
  for (auto i = 0; i < p.n_bank_; ++i) {
    const auto& bank = p.nu_bank_[i];

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
      p.filter_matches_[i_eout_filt].bins_[i_bin] = g_out;

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
        auto i_match = lower_bound_index(eo_filt.bins().begin(),
          eo_filt.bins().end(), E_out);
        p.filter_matches_[i_eout_filt].bins_[i_bin] = i_match;
      }

    }

    // Case for tallying prompt neutrons
    if (score_bin == SCORE_NU_FISSION
      || (score_bin == SCORE_PROMPT_NU_FISSION && g == 0)) {

      // Find the filter scoring index for this filter combination
      int filter_index = 0;
      double filter_weight = 1.0;
      for (auto j = 0; j < tally.filters().size(); ++j) {
        auto i_filt = tally.filters(j);
        auto& match {p.filter_matches_[i_filt]};
        auto i_bin = match.i_bin_;
        filter_index += match.bins_[i_bin] * tally.strides(j);
        filter_weight *= match.weights_[i_bin];
      }

      // Update tally results
      #pragma omp atomic
      tally.results_(filter_index, i_score, TallyResult::VALUE) += score*filter_weight;

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
              auto& match {p.filter_matches_[i_filt]};
              auto i_bin = match.i_bin_;
              filter_weight *= match.weights_[i_bin];
            }

            score_fission_delayed_dg(i_tally, d_bin, score*filter_weight,
              i_score, p.filter_matches_);
          }
        }

      // If the delayed group filter is not present, add score to tally
      } else {

        // Find the filter index and weight for this filter combination
        int filter_index = 0;
        double filter_weight = 1.;
        for (auto j = 0; j < tally.filters().size(); ++j) {
          auto i_filt = tally.filters(j);
          auto& match {p.filter_matches_[i_filt]};
          auto i_bin = match.i_bin_;
          filter_index += match.bins_[i_bin] * tally.strides(j);
          filter_weight *= match.weights_[i_bin];
        }

        // Update tally results
        #pragma omp atomic
        tally.results_(filter_index, i_score, TallyResult::VALUE) += score*filter_weight;
      }
    }
  }

  // Reset outgoing energy bin and score index
  p.filter_matches_[i_eout_filt].bins_[i_bin] = bin_energyout;
}

double get_nuclide_xs(const Particle& p, int i_nuclide, int score_bin) {
  const auto& nuc {*data::nuclides[i_nuclide]};

  // Get reaction object, or return 0 if reaction is not present
  auto m = nuc.reaction_index_[score_bin];
  if (m == C_NONE) return 0.0;
  const auto& rxn {*nuc.reactions_[m]};
  const auto& micro {p.neutron_xs_[i_nuclide]};

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
    const auto& xs {rxn.xs_[i_temp]};
    double value;
    if (i_grid >= xs.threshold) {
      value = ((1.0 - f) * xs.value[i_grid-xs.threshold]
        + f * xs.value[i_grid-xs.threshold+1]);
    } else {
      value = 0.0;
    }

    if (settings::run_mode == RunMode::EIGENVALUE && score_bin == HEATING_LOCAL) {
      // Determine kerma for fission as (EFR + EGP + EGD + EB)*sigma_f
      double kerma_fission = nuc.fragments_ ?
        ((*nuc.fragments_)(p.E_last_) + (*nuc.betas_)(p.E_last_) +
         (*nuc.prompt_photons_)(p.E_last_) + (*nuc.delayed_photons_)(p.E_last_))
        * micro.fission : 0.0;

      // Determine non-fission kerma as difference
      double kerma_non_fission = value - kerma_fission;

      // Re-weight non-fission kerma by keff to properly balance energy release
      // and deposition. See D. P. Griesheimer, S. J. Douglass, and M. H. Stedry,
      // "Self-consistent energy normalization for quasistatic reactor
      // calculations", Proc. PHYSOR, Cambridge, UK, Mar 29-Apr 2, 2020.
      value = simulation::keff * kerma_non_fission + kerma_fission;
    }
    return value;
  } else {
    // For multipole, calculate (n,gamma) from other reactions
    return rxn.mt_ == N_GAMMA ? micro.absorption - micro.fission : 0.0;
  }
  return 0.0;
}

//! Update tally results for continuous-energy tallies with any estimator.
//
//! For analog tallies, the flux estimate depends on the score type so the flux
//! argument is really just used for filter weights.  The atom_density argument
//! is not used for analog tallies.

void
score_general_ce(Particle& p, int i_tally, int start_index, int filter_index,
  double filter_weight, int i_nuclide, double atom_density, double flux)
{
  Tally& tally {*model::tallies[i_tally]};

  // Get the pre-collision energy of the particle.
  auto E = p.E_last_;

  using Type = Particle::Type;

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
          score = p.wgt_last_ + p.wgt_absorb_;
        } else {
          score = p.wgt_last_;
        }

        if (p.type_ == Type::neutron || p.type_ == Type::photon) {
          score *= flux / p.macro_xs_.total;
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
          score = (p.wgt_last_ + p.wgt_absorb_) * flux;
        } else {
          score = p.wgt_last_ * flux;
        }

      } else {
        if (i_nuclide >= 0) {
          if (p.type_ == Type::neutron) {
            score = p.neutron_xs_[i_nuclide].total * atom_density * flux;
          } else if (p.type_ == Type::photon) {
            score = p.photon_xs_[i_nuclide].total * atom_density * flux;
          }
        } else {
          score = p.macro_xs_.total * flux;
        }
      }
      break;


    case SCORE_INVERSE_VELOCITY:
      if (p.type_ != Type::neutron) continue;

      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // All events score to an inverse velocity bin. We actually use a
        // collision estimator in place of an analog one since there is no way
        // to count 'events' exactly for the inverse velocity
        if (settings::survival_biasing) {
          // We need to account for the fact that some weight was already
          // absorbed
          score = p.wgt_last_ + p.wgt_absorb_;
        } else {
          score = p.wgt_last_;
        }
        score *= flux / p.macro_xs_.total;
      } else {
        score = flux;
      }
      // Score inverse velocity in units of s/cm.
      score /= std::sqrt(2. * E / MASS_NEUTRON_EV) * C_LIGHT * 100.;
      break;


    case SCORE_SCATTER:
      if (p.type_ != Type::neutron && p.type_ != Type::photon) continue;

      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // Skip any event where the particle didn't scatter
        if (p.event_ != TallyEvent::SCATTER) continue;
        // Since only scattering events make it here, again we can use the
        // weight entering the collision as the estimator for the reaction rate
        score = p.wgt_last_ * flux;
      } else {
        if (i_nuclide >= 0) {
          if (p.type_ == Type::neutron) {
            const auto& micro = p.neutron_xs_[i_nuclide];
            score = (micro.total - micro.absorption) * atom_density * flux;
          } else {
            const auto& micro = p.photon_xs_[i_nuclide];
            score = (micro.coherent + micro.incoherent) * atom_density * flux;
          }
        } else {
          if (p.type_ == Type::neutron) {
            score = (p.macro_xs_.total - p.macro_xs_.absorption) * flux;
          } else {
            score = (p.macro_xs_.coherent + p.macro_xs_.incoherent) * flux;
          }
        }
      }
      break;


    case SCORE_NU_SCATTER:
      if (p.type_ != Type::neutron) continue;

      // Only analog estimators are available.
      // Skip any event where the particle didn't scatter
      if (p.event_ != TallyEvent::SCATTER) continue;
      // For scattering production, we need to use the pre-collision weight
      // times the yield as the estimate for the number of neutrons exiting a
      // reaction with neutrons in the exit channel
      if (p.event_mt_ == ELASTIC || p.event_mt_ == N_LEVEL ||
        (p.event_mt_ >= N_N1 && p.event_mt_ <= N_NC)) {
        // Don't waste time on very common reactions we know have
        // multiplicities of one.
        score = p.wgt_last_ * flux;
      } else {
        // Get yield and apply to score
        auto m = data::nuclides[p.event_nuclide_]->reaction_index_[p.event_mt_];
        const auto& rxn {*data::nuclides[p.event_nuclide_]->reactions_[m]};
        score = p.wgt_last_ * flux * (*rxn.products_[0].yield_)(E);
      }
      break;


    case SCORE_ABSORPTION:
      if (p.type_ != Type::neutron && p.type_ != Type::photon) continue;

      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No absorption events actually occur if survival biasing is on --
          // just use weight absorbed in survival biasing
          score = p.wgt_absorb_ * flux;
        } else {
          // Skip any event where the particle wasn't absorbed
          if (p.event_ == TallyEvent::SCATTER) continue;
          // All fission and absorption events will contribute here, so we
          // can just use the particle's weight entering the collision
          score = p.wgt_last_ * flux;
        }
      } else {
        if (i_nuclide >= 0) {
          if (p.type_ == Type::neutron) {
            score = p.neutron_xs_[i_nuclide].absorption * atom_density * flux;
          } else {
            const auto& xs = p.photon_xs_[i_nuclide];
            score = (xs.total - xs.coherent - xs.incoherent) * atom_density * flux;
          }
        } else {
          if (p.type_ == Type::neutron) {
            score = p.macro_xs_.absorption * flux;
          } else {
            score = (p.macro_xs_.photoelectric + p.macro_xs_.pair_production) * flux;
          }
        }
      }
      break;


    case SCORE_FISSION:
      if (p.macro_xs_.absorption == 0) continue;
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // fission
          if (p.neutron_xs_[p.event_nuclide_].absorption > 0) {
            score = p.wgt_absorb_
              * p.neutron_xs_[p.event_nuclide_].fission
              / p.neutron_xs_[p.event_nuclide_].absorption * flux;
          } else {
            score = 0.;
          }
        } else {
          // Skip any non-absorption events
          if (p.event_ == TallyEvent::SCATTER) continue;
          // All fission events will contribute, so again we can use particle's
          // weight entering the collision as the estimate for the fission
          // reaction rate
          score = p.wgt_last_
            * p.neutron_xs_[p.event_nuclide_].fission
            / p.neutron_xs_[p.event_nuclide_].absorption * flux;
        }
      } else {
        if (i_nuclide >= 0) {
          score = p.neutron_xs_[i_nuclide].fission * atom_density * flux;
        } else {
          score = p.macro_xs_.fission * flux;
        }
      }
      break;


    case SCORE_NU_FISSION:
      if (p.macro_xs_.absorption == 0) continue;
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing || p.fission_) {
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
          if (p.neutron_xs_[p.event_nuclide_].absorption > 0) {
            score = p.wgt_absorb_
              * p.neutron_xs_[p.event_nuclide_].nu_fission
              / p.neutron_xs_[p.event_nuclide_].absorption * flux;
          } else {
            score = 0.;
          }
        } else {
          // Skip any non-fission events
          if (!p.fission_) continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score.
          score = simulation::keff * p.wgt_bank_ * flux;
        }
      } else {
        if (i_nuclide >= 0) {
          score = p.neutron_xs_[i_nuclide].nu_fission * atom_density
            * flux;
        } else {
          score = p.macro_xs_.nu_fission * flux;
        }
      }
      break;


    case SCORE_PROMPT_NU_FISSION:
      if (p.macro_xs_.absorption == 0) continue;
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing || p.fission_) {
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
          if (p.neutron_xs_[p.event_nuclide_].absorption > 0) {
            score = p.wgt_absorb_
              * p.neutron_xs_[p.event_nuclide_].fission
              * data::nuclides[p.event_nuclide_]
              ->nu(E, ReactionProduct::EmissionMode::prompt)
              / p.neutron_xs_[p.event_nuclide_].absorption * flux;
          } else {
            score = 0.;
          }
        } else {
          // Skip any non-fission events
          if (!p.fission_) continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score.
          auto n_delayed = std::accumulate(p.n_delayed_bank_,
            p.n_delayed_bank_+MAX_DELAYED_GROUPS, 0);
          auto prompt_frac = 1. - n_delayed / static_cast<double>(p.n_bank_);
          score = simulation::keff * p.wgt_bank_ * prompt_frac * flux;
        }
      } else {
        if (i_nuclide >= 0) {
          score = p.neutron_xs_[i_nuclide].fission
            * data::nuclides[i_nuclide]
            ->nu(E, ReactionProduct::EmissionMode::prompt)
            * atom_density * flux;
        } else {
          score = 0.;
          // Add up contributions from each nuclide in the material.
          if (p.material_ != MATERIAL_VOID) {
            const Material& material {*model::materials[p.material_]};
            for (auto i = 0; i < material.nuclide_.size(); ++i) {
              auto j_nuclide = material.nuclide_[i];
              auto atom_density = material.atom_density_(i);
              score += p.neutron_xs_[j_nuclide].fission
                * data::nuclides[j_nuclide]
                ->nu(E, ReactionProduct::EmissionMode::prompt)
                * atom_density * flux;
            }
          }
        }
      }
      break;


    case SCORE_DELAYED_NU_FISSION:
      if (p.macro_xs_.absorption == 0) continue;
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing || p.fission_) {
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
          if (p.neutron_xs_[p.event_nuclide_].absorption > 0
            && data::nuclides[p.event_nuclide_]->fissionable_) {
            if (tally.delayedgroup_filter_ != C_NONE) {
              auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
              const DelayedGroupFilter& filt
                {*dynamic_cast<DelayedGroupFilter*>(
                model::tally_filters[i_dg_filt].get())};
              // Tally each delayed group bin individually
              for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                auto dg = filt.groups()[d_bin];
                auto yield = data::nuclides[p.event_nuclide_]
                  ->nu(E, ReactionProduct::EmissionMode::delayed, dg);
                score = p.wgt_absorb_ * yield
                  * p.neutron_xs_[p.event_nuclide_].fission
                  / p.neutron_xs_[p.event_nuclide_].absorption * flux;
                score_fission_delayed_dg(i_tally, d_bin, score,
                  score_index, p.filter_matches_);
              }
              continue;
            } else {
              // If the delayed group filter is not present, compute the score
              // by multiplying the absorbed weight by the fraction of the
              // delayed-nu-fission xs to the absorption xs
              score = p.wgt_absorb_
                * p.neutron_xs_[p.event_nuclide_].fission
                * data::nuclides[p.event_nuclide_]
                ->nu(E, ReactionProduct::EmissionMode::delayed)
                / p.neutron_xs_[p.event_nuclide_].absorption *flux;
            }
          }
        } else {
          // Skip any non-fission events
          if (!p.fission_) continue;
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
              score = simulation::keff * p.wgt_bank_ / p.n_bank_
                * p.n_delayed_bank_[d-1] * flux;
              score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches_);
            }
            continue;
          } else {
            // Add the contribution from all delayed groups
            auto n_delayed = std::accumulate(p.n_delayed_bank_,
              p.n_delayed_bank_+MAX_DELAYED_GROUPS, 0);
            score = simulation::keff * p.wgt_bank_ / p.n_bank_ * n_delayed
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
              score = p.neutron_xs_[i_nuclide].fission * yield
                * atom_density * flux;
              score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches_);
            }
            continue;
          } else {
            // If the delayed group filter is not present, compute the score
            // by multiplying the delayed-nu-fission macro xs by the flux
            score = p.neutron_xs_[i_nuclide].fission
              * data::nuclides[i_nuclide]
              ->nu(E, ReactionProduct::EmissionMode::delayed)
              * atom_density * flux;
          }
        } else {
          // Need to add up contributions for each nuclide
          if (tally.delayedgroup_filter_ != C_NONE) {
            auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
            const DelayedGroupFilter& filt
              {*dynamic_cast<DelayedGroupFilter*>(
              model::tally_filters[i_dg_filt].get())};
            if (p.material_ != MATERIAL_VOID) {
              const Material& material {*model::materials[p.material_]};
              for (auto i = 0; i < material.nuclide_.size(); ++i) {
                auto j_nuclide = material.nuclide_[i];
                auto atom_density = material.atom_density_(i);
                // Tally each delayed group bin individually
                for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                  auto d = filt.groups()[d_bin];
                  auto yield = data::nuclides[j_nuclide]
                    ->nu(E, ReactionProduct::EmissionMode::delayed, d);
                  score = p.neutron_xs_[j_nuclide].fission * yield
                    * atom_density * flux;
                  score_fission_delayed_dg(i_tally, d_bin, score,
                    score_index, p.filter_matches_);
                }
              }
            }
            continue;
          } else {
            score = 0.;
            if (p.material_ != MATERIAL_VOID) {
              const Material& material {*model::materials[p.material_]};
              for (auto i = 0; i < material.nuclide_.size(); ++i) {
                auto j_nuclide = material.nuclide_[i];
                auto atom_density = material.atom_density_(i);
                score += p.neutron_xs_[j_nuclide].fission
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
      if (p.macro_xs_.absorption == 0) continue;
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // delayed-nu-fission
          const auto& nuc {*data::nuclides[p.event_nuclide_]};
          if (p.neutron_xs_[p.event_nuclide_].absorption > 0
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
                score = p.wgt_absorb_ * yield
                  * p.neutron_xs_[p.event_nuclide_].fission
                  / p.neutron_xs_[p.event_nuclide_].absorption
                  * rate * flux;
                score_fission_delayed_dg(i_tally, d_bin, score,
                  score_index, p.filter_matches_);
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
                score += rate * p.wgt_absorb_
                  * p.neutron_xs_[p.event_nuclide_].fission * yield
                  / p.neutron_xs_[p.event_nuclide_].absorption * flux;
              }
            }
          }
        } else {
          // Skip any non-fission events
          if (!p.fission_) continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score. Loop over the neutrons produced from fission and check which
          // ones are delayed. If a delayed neutron is encountered, add its
          // contribution to the fission bank to the score.
          score = 0.;
          for (auto i = 0; i < p.n_bank_; ++i) {
            const auto& bank = p.nu_bank_[i];
            auto g = bank.delayed_group;
            if (g != 0) {
              const auto& nuc {*data::nuclides[p.event_nuclide_]};
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
                      score_index, p.filter_matches_);
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
              score = p.neutron_xs_[i_nuclide].fission * yield * flux
                * atom_density * rate;
              score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches_);
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
              score += p.neutron_xs_[i_nuclide].fission * flux
                * yield * atom_density * rate;
            }
          }
        } else {
          if (tally.delayedgroup_filter_ != C_NONE) {
            auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
            const DelayedGroupFilter& filt
              {*dynamic_cast<DelayedGroupFilter*>(
              model::tally_filters[i_dg_filt].get())};
            if (p.material_ != MATERIAL_VOID) {
              const Material& material {*model::materials[p.material_]};
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
                    score = p.neutron_xs_[j_nuclide].fission * yield
                      * flux * atom_density * rate;
                    score_fission_delayed_dg(i_tally, d_bin, score,
                      score_index, p.filter_matches_);
                  }
                }
              }
            }
            continue;
          } else {
            score = 0.;
            if (p.material_ != MATERIAL_VOID) {
              const Material& material {*model::materials[p.material_]};
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
                    score += p.neutron_xs_[j_nuclide].fission
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
      if (p.macro_xs_.absorption == 0.) continue;
      score = 0.;
      // Kappa-fission values are determined from the Q-value listed for the
      // fission cross section.
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // fission scaled by the Q-value
          const auto& nuc {*data::nuclides[p.event_nuclide_]};
          if (p.neutron_xs_[p.event_nuclide_].absorption > 0
            && nuc.fissionable_) {
            const auto& rxn {*nuc.fission_rx_[0]};
            score = p.wgt_absorb_ * rxn.q_value_
              * p.neutron_xs_[p.event_nuclide_].fission
              / p.neutron_xs_[p.event_nuclide_].absorption * flux;
          }
        } else {
          // Skip any non-absorption events
          if (p.event_ == TallyEvent::SCATTER) continue;
          // All fission events will contribute, so again we can use particle's
          // weight entering the collision as the estimate for the fission
          // reaction rate
          const auto& nuc {*data::nuclides[p.event_nuclide_]};
          if (p.neutron_xs_[p.event_nuclide_].absorption > 0
            && nuc.fissionable_) {
            const auto& rxn {*nuc.fission_rx_[0]};
            score = p.wgt_last_ * rxn.q_value_
              * p.neutron_xs_[p.event_nuclide_].fission
              / p.neutron_xs_[p.event_nuclide_].absorption * flux;
          }
        }
      } else {
        if (i_nuclide >= 0) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          if (nuc.fissionable_) {
            const auto& rxn {*nuc.fission_rx_[0]};
            score = rxn.q_value_ * p.neutron_xs_[i_nuclide].fission
              * atom_density * flux;
          }
        } else if (p.material_ != MATERIAL_VOID) {
          const Material& material {*model::materials[p.material_]};
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto j_nuclide = material.nuclide_[i];
            auto atom_density = material.atom_density_(i);
            const auto& nuc {*data::nuclides[j_nuclide]};
            if (nuc.fissionable_) {
              const auto& rxn {*nuc.fission_rx_[0]};
              score += rxn.q_value_ * p.neutron_xs_[j_nuclide].fission
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
      if (p.type_ != Type::neutron) continue;

      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // Check if event MT matches
        if (p.event_mt_ != ELASTIC) continue;
        score = p.wgt_last_ * flux;
      } else {
        if (i_nuclide >= 0) {
          if (p.neutron_xs_[i_nuclide].elastic == CACHE_INVALID)
            data::nuclides[i_nuclide]->calculate_elastic_xs(p);
          score = p.neutron_xs_[i_nuclide].elastic * atom_density * flux;
        } else {
          score = 0.;
          if (p.material_ != MATERIAL_VOID) {
            const Material& material {*model::materials[p.material_]};
            for (auto i = 0; i < material.nuclide_.size(); ++i) {
              auto j_nuclide = material.nuclide_[i];
              auto atom_density = material.atom_density_(i);
              if (p.neutron_xs_[j_nuclide].elastic == CACHE_INVALID)
                data::nuclides[j_nuclide]->calculate_elastic_xs(p);
              score += p.neutron_xs_[j_nuclide].elastic * atom_density
                * flux;
            }
          }
        }
      }
      break;

    case SCORE_FISS_Q_PROMPT:
    case SCORE_FISS_Q_RECOV:
      if (p.macro_xs_.absorption == 0.) continue;
      score = score_fission_q(p, score_bin, tally, flux, i_nuclide, atom_density);
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
      if (!simulation::need_depletion_rx) goto default_case;

      if (p.type_ != Type::neutron) continue;

      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // Check if the event MT matches
        if (p.event_mt_ != score_bin) continue;
        score = p.wgt_last_ * flux;
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
          score = p.neutron_xs_[i_nuclide].reaction[m] * atom_density
            * flux;
        } else {
          score = 0.;
          if (p.material_ != MATERIAL_VOID) {
            const Material& material {*model::materials[p.material_]};
            for (auto i = 0; i < material.nuclide_.size(); ++i) {
              auto j_nuclide = material.nuclide_[i];
              auto atom_density = material.atom_density_(i);
              score += p.neutron_xs_[j_nuclide].reaction[m]
                * atom_density * flux;
            }
          }
        }
      }
      break;


    case COHERENT:
    case INCOHERENT:
    case PHOTOELECTRIC:
    case PAIR_PROD:
      if (p.type_ != Type::photon) continue;

      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (score_bin == PHOTOELECTRIC) {
          // Photoelectric events are assigned an MT value corresponding to the
          // shell cross section. Also, photons below the energy cutoff are
          // assumed to have been absorbed via photoelectric absorption
          if ((p.event_mt_ < 534 || p.event_mt_ > 572) &&
              p.event_mt_ != REACTION_NONE) continue;
        } else {
          if (p.event_mt_ != score_bin) continue;
        }
        score = p.wgt_last_ * flux;
      } else {
        if (i_nuclide >= 0) {
          const auto& micro = p.photon_xs_[i_nuclide];
          double xs =
            (score_bin == COHERENT) ? micro.coherent :
            (score_bin == INCOHERENT) ? micro.incoherent :
            (score_bin == PHOTOELECTRIC) ? micro.photoelectric :
            micro.pair_production;
          score = xs * atom_density * flux;
        } else {
          double xs =
            (score_bin == COHERENT) ? p.macro_xs_.coherent :
            (score_bin == INCOHERENT) ? p.macro_xs_.incoherent :
            (score_bin == PHOTOELECTRIC) ? p.macro_xs_.photoelectric :
            p.macro_xs_.pair_production;
          score = xs * flux;
        }
      }
      break;


    case HEATING:
      if (p.type_ == Type::neutron) {
        score = score_neutron_heating(p, tally, flux, HEATING,
            i_nuclide, atom_density);
      } else {
        // The energy deposited is the difference between the pre-collision and
        // post-collision energy...
        score = E - p.E_;

        // ...less the energy of any secondary particles since they will be
        // transported individually later
        const auto& bank = p.secondary_bank_;
        for (auto it = bank.end() - p.n_bank_second_; it < bank.end(); ++it) {
          score -= it->E;
        }

        score *= p.wgt_last_;
      }
      break;

    default:
    default_case:

      // The default block is really only meant for redundant neutron reactions
      // (e.g. 444, 901)
      if (p.type_ != Type::neutron) continue;

      if (tally.estimator_ == TallyEstimator::ANALOG) {

        // Any other score is assumed to be a MT number. Thus, we just need
        // to check if it matches the MT number of the event
        if (p.event_mt_ != score_bin) continue;
        score = p.wgt_last_*flux;
      } else {
        // Any other cross section has to be calculated on-the-fly
        if (score_bin < 2) fatal_error("Invalid score type on tally "
          + std::to_string(tally.id_));
        score = 0.;
        if (i_nuclide >= 0) {
          score = get_nuclide_xs(p, i_nuclide, score_bin) * atom_density * flux;
        } else if (p.material_ != MATERIAL_VOID) {
          const Material& material {*model::materials[p.material_]};
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto j_nuclide = material.nuclide_[i];
            auto atom_density = material.atom_density_(i);
            score += get_nuclide_xs(p, j_nuclide, score_bin) * atom_density * flux;
          }
        }
      }
    }

    // Add derivative information on score for differential tallies.
    if (tally.deriv_ != C_NONE)
      apply_derivative_to_score(p, i_tally, i_nuclide, atom_density, score_bin,
        score);

    // Update tally results
    #pragma omp atomic
    tally.results_(filter_index, score_index, TallyResult::VALUE) += score*filter_weight;
  }
}

//! Update tally results for multigroup tallies with any estimator.
//
//! For analog tallies, the flux estimate depends on the score type so the flux
//! argument is really just used for filter weights.

void
score_general_mg(Particle& p, int i_tally, int start_index, int filter_index,
  double filter_weight, int i_nuclide, double atom_density, double flux)
{
  auto& tally {*model::tallies[i_tally]};

  // Set the direction and group to use with get_xs
  Direction p_u;
  int p_g;
  if (tally.estimator_ == TallyEstimator::ANALOG
    || tally.estimator_ == TallyEstimator::COLLISION) {

    if (settings::survival_biasing) {

      // Then we either are alive and had a scatter (and so g changed),
      // or are dead and g did not change
      if (p.alive_) {
        p_u = p.u_last_;
        p_g = p.g_last_;
      } else {
        p_u = p.u_local();
        p_g = p.g_;
      }
    } else if (p.event_ == TallyEvent::SCATTER) {

      // Then the energy group has been changed by the scattering routine
      // meaning gin is now in p % last_g
      p_u = p.u_last_;
      p_g = p.g_last_;
    } else {

      // No scatter, no change in g.
      p_u = p.u_local();
      p_g = p.g_;
    }
  } else {

    // No actual collision so g has not changed.
    p_u = p.u_local();
    p_g = p.g_;
  }

  // For shorthand, assign pointers to the material and nuclide xs set
  auto& nuc_xs = data::mg.nuclides_[i_nuclide];
  auto& macro_xs = data::mg.macro_xs_[p.material_];

  // Find the temperature and angle indices of interest
  macro_xs.set_angle_index(p_u);
  if (i_nuclide >= 0) {
    nuc_xs.set_temperature_index(p.sqrtkT_);
    nuc_xs.set_angle_index(p_u);
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
        if (settings::survival_biasing) {
          // We need to account for the fact that some weight was already
          // absorbed
          score = p.wgt_last_ + p.wgt_absorb_;
        } else {
          score = p.wgt_last_;
        }
        score *= flux / p.macro_xs_.total;
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
          score = p.wgt_last_ + p.wgt_absorb_;
        } else {
          score = p.wgt_last_;
        }
        //TODO: should flux be multiplied in above instead of below?
        if (i_nuclide >= 0) {
          score *= flux * atom_density * nuc_xs.get_xs(MgxsType::TOTAL, p_g)
            / macro_xs.get_xs(MgxsType::TOTAL, p_g);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux * nuc_xs.get_xs(MgxsType::TOTAL, p_g);
        } else {
          score = p.macro_xs_.total * flux;
        }
      }
      break;


    case SCORE_INVERSE_VELOCITY:
      if (tally.estimator_ == TallyEstimator::ANALOG
        || tally.estimator_ == TallyEstimator::COLLISION) {
        // All events score to an inverse velocity bin. We actually use a
        // collision estimator in place of an analog one since there is no way
        // to count 'events' exactly for the inverse velocity
        if (settings::survival_biasing) {
          // We need to account for the fact that some weight was already
          // absorbed
          score = p.wgt_last_ + p.wgt_absorb_;
        } else {
          score = p.wgt_last_;
        }
        if (i_nuclide >= 0) {
          score *= flux * nuc_xs.get_xs(MgxsType::INVERSE_VELOCITY, p_g)
            / macro_xs.get_xs(MgxsType::TOTAL, p_g);
        } else {
          score *= flux * macro_xs.get_xs(MgxsType::INVERSE_VELOCITY, p_g)
            / macro_xs.get_xs(MgxsType::TOTAL, p_g);
        }
      } else {
        if (i_nuclide >= 0) {
          score = flux * nuc_xs.get_xs(MgxsType::INVERSE_VELOCITY, p_g);
        } else {
          score = flux * macro_xs.get_xs(MgxsType::INVERSE_VELOCITY, p_g);
        }
      }
      break;


    case SCORE_SCATTER:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // Skip any event where the particle didn't scatter
        if (p.event_ != TallyEvent::SCATTER) continue;
        // Since only scattering events make it here, again we can use the
        // weight entering the collision as the estimator for the reaction rate
        score = p.wgt_last_ * flux;
        if (i_nuclide >= 0) {
          score *= atom_density * nuc_xs.get_xs(
              MgxsType::SCATTER_FMU, p.g_last_, &p.g_, &p.mu_, nullptr)
            / macro_xs.get_xs(
              MgxsType::SCATTER_FMU, p.g_last_, &p.g_, &p.mu_, nullptr);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux * nuc_xs.get_xs(
            MgxsType::SCATTER, p_g, nullptr, &p.mu_, nullptr);
        } else {
          score = flux * macro_xs.get_xs(
            MgxsType::SCATTER, p_g, nullptr, &p.mu_, nullptr);
        }
      }
      break;


    case SCORE_NU_SCATTER:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        // Skip any event where the particle didn't scatter
        if (p.event_ != TallyEvent::SCATTER) continue;
        // For scattering production, we need to use the pre-collision weight
        // times the multiplicity as the estimate for the number of neutrons
        // exiting a reaction with neutrons in the exit channel
        score = p.wgt_ * flux;
        // Since we transport based on material data, the angle selected
        // was not selected from the f(mu) for the nuclide.  Therefore
        // adjust the score by the actual probability for that nuclide.
        if (i_nuclide >= 0) {
          score *= atom_density
            * nuc_xs.get_xs(MgxsType::NU_SCATTER_FMU, p.g_last_, &p.g_,
                            &p.mu_, nullptr)
            / macro_xs.get_xs(MgxsType::NU_SCATTER_FMU, p.g_last_, &p.g_,
                              &p.mu_, nullptr);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux * nuc_xs.get_xs(MgxsType::NU_SCATTER, p_g);
        } else {
          score = flux * macro_xs.get_xs(MgxsType::NU_SCATTER, p_g);
        }
      }
      break;


    case SCORE_ABSORPTION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No absorption events actually occur if survival biasing is on --
          // just use weight absorbed in survival biasing
          score = p.wgt_absorb_ * flux;
        } else {
          // Skip any event where the particle wasn't absorbed
          if (p.event_ == TallyEvent::SCATTER) continue;
          // All fission and absorption events will contribute here, so we
          // can just use the particle's weight entering the collision
          score = p.wgt_last_ * flux;
        }
        if (i_nuclide >= 0) {
          score *= atom_density * nuc_xs.get_xs(MgxsType::ABSORPTION, p_g)
            / macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux
            * nuc_xs.get_xs(MgxsType::ABSORPTION, p_g);
        } else {
          score = p.macro_xs_.absorption * flux;
        }
      }
      break;


    case SCORE_FISSION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing) {
          // No fission events occur if survival biasing is on -- need to
          // calculate fraction of absorptions that would have resulted in
          // fission
          score = p.wgt_absorb_ * flux;
        } else {
          // Skip any non-absorption events
          if (p.event_ == TallyEvent::SCATTER) continue;
          // All fission events will contribute, so again we can use particle's
          // weight entering the collision as the estimate for the fission
          // reaction rate
          score = p.wgt_last_ * flux;
        }
        if (i_nuclide >= 0) {
          score *= atom_density * nuc_xs.get_xs(MgxsType::FISSION, p_g)
            / macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
        } else {
          score *= macro_xs.get_xs(MgxsType::FISSION, p_g)
            / macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux * nuc_xs.get_xs(MgxsType::FISSION, p_g);
        } else {
          score = flux * macro_xs.get_xs(MgxsType::FISSION, p_g);
        }
      }
      break;


    case SCORE_NU_FISSION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing || p.fission_) {
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
          score = p.wgt_absorb_ * flux;
          if (i_nuclide >= 0) {
            score *= atom_density * nuc_xs.get_xs(MgxsType::NU_FISSION, p_g)
              / macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
          } else {
            score *= macro_xs.get_xs(MgxsType::NU_FISSION, p_g)
              / macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
          }
        } else {
          // Skip any non-fission events
          if (!p.fission_) continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score.
          score = simulation::keff * p.wgt_bank_ * flux;
          if (i_nuclide >= 0) {
            score *= atom_density * nuc_xs.get_xs(MgxsType::FISSION, p_g)
              / macro_xs.get_xs(MgxsType::FISSION, p_g);
          }
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux
            * nuc_xs.get_xs(MgxsType::NU_FISSION, p_g);
        } else {
          score = flux * macro_xs.get_xs(MgxsType::NU_FISSION, p_g);
        }
      }
      break;


    case SCORE_PROMPT_NU_FISSION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing || p.fission_) {
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
          score = p.wgt_absorb_ * flux;
          if (i_nuclide >= 0) {
            score *= atom_density *
              nuc_xs.get_xs(MgxsType::PROMPT_NU_FISSION, p_g)
              / macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
          } else {
            score *= macro_xs.get_xs(MgxsType::PROMPT_NU_FISSION, p_g)
              / macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
          }
        } else {
          // Skip any non-fission events
          if (!p.fission_) continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score.
          auto n_delayed = std::accumulate(p.n_delayed_bank_,
            p.n_delayed_bank_+MAX_DELAYED_GROUPS, 0);
          auto prompt_frac = 1. - n_delayed / static_cast<double>(p.n_bank_);
          score = simulation::keff * p.wgt_bank_ * prompt_frac * flux;
          if (i_nuclide >= 0) {
            score *= atom_density * nuc_xs.get_xs(MgxsType::FISSION, p_g)
              / macro_xs.get_xs(MgxsType::FISSION, p_g);
          }
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux *
            nuc_xs.get_xs(MgxsType::PROMPT_NU_FISSION, p_g);
        } else {
          score = flux * macro_xs.get_xs(MgxsType::PROMPT_NU_FISSION, p_g);
        }
      }
      break;


    case SCORE_DELAYED_NU_FISSION:
      if (tally.estimator_ == TallyEstimator::ANALOG) {
        if (settings::survival_biasing || p.fission_) {
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
          double abs_xs = macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
          if (abs_xs > 0.) {
            if (tally.delayedgroup_filter_ != C_NONE) {
              auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
              const DelayedGroupFilter& filt
                {*dynamic_cast<DelayedGroupFilter*>(
                model::tally_filters[i_dg_filt].get())};
              // Tally each delayed group bin individually
              for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                auto d = filt.groups()[d_bin] - 1;
                score = p.wgt_absorb_ * flux;
                if (i_nuclide >= 0) {
                  score *= nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                                         nullptr, nullptr, &d)
                    / abs_xs;
                } else {
                  score *= macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                                           nullptr, nullptr, &d)
                    / abs_xs;
                }
                score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches_);
              }
              continue;
            } else {
              // If the delayed group filter is not present, compute the score
              // by multiplying the absorbed weight by the fraction of the
              // delayed-nu-fission xs to the absorption xs
              score = p.wgt_absorb_ * flux;
              if (i_nuclide >= 0) {
                score *= nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g)
                  / abs_xs;
              } else {
                score *= macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g)
                  / abs_xs;
              }
            }
          }
        } else {
          // Skip any non-fission events
          if (!p.fission_) continue;
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
              score = simulation::keff * p.wgt_bank_ / p.n_bank_
                * p.n_delayed_bank_[d-1] * flux;
              if (i_nuclide >= 0) {
                score *= atom_density * nuc_xs.get_xs(MgxsType::FISSION, p_g)
                  / macro_xs.get_xs(MgxsType::FISSION, p_g);
              }
              score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches_);
            }
            continue;
          } else {
            // Add the contribution from all delayed groups
            auto n_delayed = std::accumulate(p.n_delayed_bank_,
              p.n_delayed_bank_+MAX_DELAYED_GROUPS, 0);
            score = simulation::keff * p.wgt_bank_ / p.n_bank_ * n_delayed
              * flux;
            if (i_nuclide >= 0) {
              score *= atom_density * nuc_xs.get_xs(MgxsType::FISSION, p_g)
                / macro_xs.get_xs(MgxsType::FISSION, p_g);
            }
          }
        }
      } else {
        if (tally.delayedgroup_filter_ != C_NONE) {
          auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
          const DelayedGroupFilter& filt
            {*dynamic_cast<DelayedGroupFilter*>(
            model::tally_filters[i_dg_filt].get())};
          // Tally each delayed group bin individually
          for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
            auto d = filt.groups()[d_bin] - 1;
            if (i_nuclide >= 0) {
              score = flux * atom_density * nuc_xs.get_xs(
                MgxsType::DELAYED_NU_FISSION, p_g, nullptr, nullptr, &d);
            } else {
              score = flux * macro_xs.get_xs(
                MgxsType::DELAYED_NU_FISSION, p_g, nullptr, nullptr, &d);
            }
            score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches_);
          }
          continue;
        } else {
          if (i_nuclide >= 0) {
            score = flux * atom_density *
              nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g);
          } else {
            score = flux * macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g);
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
          double abs_xs = macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
          if (abs_xs > 0) {
            if (tally.delayedgroup_filter_ != C_NONE) {
              auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
              const DelayedGroupFilter& filt
                {*dynamic_cast<DelayedGroupFilter*>(
                model::tally_filters[i_dg_filt].get())};
              // Tally each delayed group bin individually
              for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                auto d = filt.groups()[d_bin] - 1;
                score = p.wgt_absorb_ * flux;
                if (i_nuclide >= 0) {
                  score *= nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g,
                                         nullptr, nullptr, &d)
                    * nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                                    nullptr, nullptr, &d) / abs_xs;
                } else {
                  score *= macro_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                                           nullptr, &d)
                    * macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                                      nullptr, nullptr, &d) / abs_xs;
                }
                score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches_);
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
                  score += p.wgt_absorb_ * flux
                    * nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                                    nullptr, &d)
                    * nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                                    nullptr, nullptr, &d) / abs_xs;
                } else {
                  score += p.wgt_absorb_ * flux
                    * macro_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                                      nullptr, &d)
                    * macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                                      nullptr, nullptr, &d) / abs_xs;
                }
              }
            }
          }
        } else {
          // Skip any non-fission events
          if (!p.fission_) continue;
          // If there is no outgoing energy filter, than we only need to score
          // to one bin. For the score to be 'analog', we need to score the
          // number of particles that were banked in the fission bank. Since
          // this was weighted by 1/keff, we multiply by keff to get the proper
          // score. Loop over the neutrons produced from fission and check which
          // ones are delayed. If a delayed neutron is encountered, add its
          // contribution to the fission bank to the score.
          score = 0.;
          for (auto i = 0; i < p.n_bank_; ++i) {
            const auto& bank = p.nu_bank_[i];
            auto d = bank.delayed_group - 1;
            if (d != -1) {
              if (i_nuclide >= 0) {
                score += simulation::keff * atom_density * bank.wgt * flux
                  * nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d)
                  * nuc_xs.get_xs(MgxsType::FISSION, p_g)
                  / macro_xs.get_xs(MgxsType::FISSION, p_g);
              } else {
                score += simulation::keff * bank.wgt * flux
                  * macro_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d);
              }
              if (tally.delayedgroup_filter_ != C_NONE) {
                auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
                const DelayedGroupFilter& filt
                  {*dynamic_cast<DelayedGroupFilter*>(
                  model::tally_filters[i_dg_filt].get())};
                // Find the corresponding filter bin and then score
                for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                  auto dg = filt.groups()[d_bin];
                  if (dg == d + 1)
                    score_fission_delayed_dg(i_tally, d_bin, score,
                                             score_index, p.filter_matches_);
                }
                score = 0.;
              }
            }
          }
          if (tally.delayedgroup_filter_ != C_NONE) continue;
        }
      } else {
        if (tally.delayedgroup_filter_ != C_NONE) {
          auto i_dg_filt = tally.filters()[tally.delayedgroup_filter_];
          const DelayedGroupFilter& filt
            {*dynamic_cast<DelayedGroupFilter*>(
            model::tally_filters[i_dg_filt].get())};
          // Tally each delayed group bin individually
          for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
            auto d = filt.groups()[d_bin] - 1;
            if (i_nuclide >= 0) {
              score = atom_density * flux
                * nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                                nullptr, &d)
                * nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g, nullptr,
                                nullptr, &d);
            } else {
              score = flux
                * macro_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr,
                                  nullptr, &d)
                * macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g, nullptr,
                                  nullptr, &d);
            }
            score_fission_delayed_dg(i_tally, d_bin, score, score_index, p.filter_matches_);
          }
          continue;
        } else {
          score = 0.;
          for (auto d = 0; d < data::mg.num_delayed_groups_; ++d) {
            if (i_nuclide >= 0) {
              score += atom_density * flux
                * nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d)
                * nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g, nullptr, nullptr, &d);
            } else {
              score += flux
                * macro_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d)
                * macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g, nullptr, nullptr, &d);
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
          score = p.wgt_absorb_ * flux;
        } else {
          // Skip any non-absorption events
          if (p.event_ == TallyEvent::SCATTER) continue;
          // All fission events will contribute, so again we can use particle's
          // weight entering the collision as the estimate for the fission
          // reaction rate
          score = p.wgt_last_ * flux;
        }
        if (i_nuclide >= 0) {
          score *= atom_density
            * nuc_xs.get_xs(MgxsType::KAPPA_FISSION, p_g)
            / macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
        } else {
          score *=
            macro_xs.get_xs(MgxsType::KAPPA_FISSION, p_g)
            / macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux * nuc_xs.get_xs(MgxsType::KAPPA_FISSION, p_g);
        } else {
          score = flux * macro_xs.get_xs(MgxsType::KAPPA_FISSION, p_g);
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
    tally.results_(filter_index, score_index, TallyResult::VALUE) += score*filter_weight;
  }
}

//! Tally rates for when the user requests a tally on all nuclides.

void
score_all_nuclides(Particle& p, int i_tally, double flux,
  int filter_index, double filter_weight)
{
  const Tally& tally {*model::tallies[i_tally]};
  const Material& material {*model::materials[p.material_]};

  // Score all individual nuclide reaction rates.
  for (auto i = 0; i < material.nuclide_.size(); ++i) {
    auto i_nuclide = material.nuclide_[i];
    auto atom_density = material.atom_density_(i);

    //TODO: consider replacing this "if" with pointers or templates
    if (settings::run_CE) {
      score_general_ce(p, i_tally, i_nuclide*tally.scores_.size(), filter_index,
        filter_weight, i_nuclide, atom_density, flux);
    } else {
      score_general_mg(p, i_tally, i_nuclide*tally.scores_.size(), filter_index,
        filter_weight, i_nuclide, atom_density, flux);
    }
  }

  // Score total material reaction rates.
  int i_nuclide = -1;
  double atom_density = 0.;
  auto n_nuclides = data::nuclides.size();
  //TODO: consider replacing this "if" with pointers or templates
  if (settings::run_CE) {
    score_general_ce(p, i_tally, n_nuclides*tally.scores_.size(), filter_index,
      filter_weight, i_nuclide, atom_density, flux);
  } else {
    score_general_mg(p, i_tally, n_nuclides*tally.scores_.size(), filter_index,
      filter_weight, i_nuclide, atom_density, flux);
  }
}

void score_analog_tally_ce(Particle& p)
{
  // Since electrons/positrons are not transported, we assign a flux of zero.
  // Note that the heating score does NOT use the flux and will be non-zero for
  // electrons/positrons.
  double flux =
    (p.type_ == Particle::Type::neutron || p.type_ == Particle::Type::photon) ?
    1.0 : 0.0;

  for (auto i_tally : model::active_analog_tallies) {
    const Tally& tally {*model::tallies[i_tally]};

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p);
    auto end = FilterBinIter(tally, true, &p.filter_matches_);
    if (filter_iter == end) continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
      if (!tally.all_nuclides_) {
        for (auto i = 0; i < tally.nuclides_.size(); ++i) {
          auto i_nuclide = tally.nuclides_[i];

          // Tally this event in the present nuclide bin if that bin represents
          // the event nuclide or the total material.  Note that the atomic
          // density argument for score_general is not used for analog tallies.
          if (i_nuclide == p.event_nuclide_ || i_nuclide == -1)
            score_general_ce(p, i_tally, i*tally.scores_.size(), filter_index,
              filter_weight, i_nuclide, -1.0, flux);
        }

      } else {
        // In the case that the user has requested to tally all nuclides, we
        // can take advantage of the fact that we know exactly how nuclide
        // bins correspond to nuclide indices.  First, tally the nuclide.
        auto i = p.event_nuclide_;
        score_general_ce(p, i_tally, i*tally.scores_.size(), filter_index,
          filter_weight, -1, -1.0, flux);

        // Now tally the total material.
        i = tally.nuclides_.size();
        score_general_ce(p, i_tally, i*tally.scores_.size(), filter_index,
          filter_weight, -1, -1.0, flux);
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : p.filter_matches_)
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
    auto end = FilterBinIter(tally, true, &p.filter_matches_);
    if (filter_iter == end) continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
      for (auto i = 0; i < tally.nuclides_.size(); ++i) {
        auto i_nuclide = tally.nuclides_[i];

        double atom_density = 0.;
        if (i_nuclide >= 0) {
          auto j = model::materials[p.material_]->mat_nuclide_index_[i_nuclide];
          if (j == C_NONE) continue;
          atom_density = model::materials[p.material_]->atom_density_(j);
        }

        score_general_mg(p, i_tally, i*tally.scores_.size(), filter_index,
          filter_weight, i_nuclide, atom_density, 1.0);
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : p.filter_matches_)
    match.bins_present_ = false;
}

void
score_tracklength_tally(Particle& p, double distance)
{
  // Determine the tracklength estimate of the flux
  double flux = p.wgt_ * distance;

  for (auto i_tally : model::active_tracklength_tallies) {
    const Tally& tally {*model::tallies[i_tally]};

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p);
    auto end = FilterBinIter(tally, true, &p.filter_matches_);
    if (filter_iter == end) continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
      if (tally.all_nuclides_) {
        if (p.material_ != MATERIAL_VOID)
          score_all_nuclides(p, i_tally, flux*filter_weight, filter_index,
            filter_weight);

      } else {
        for (auto i = 0; i < tally.nuclides_.size(); ++i) {
          auto i_nuclide = tally.nuclides_[i];

          double atom_density = 0.;
          if (i_nuclide >= 0) {
            if (p.material_ != MATERIAL_VOID) {
              auto j = model::materials[p.material_]->mat_nuclide_index_[i_nuclide];
              if (j == C_NONE) continue;
              atom_density = model::materials[p.material_]->atom_density_(j);
            }
          }

          //TODO: consider replacing this "if" with pointers or templates
          if (settings::run_CE) {
            score_general_ce(p, i_tally, i*tally.scores_.size(), filter_index,
              filter_weight, i_nuclide, atom_density, flux);
          } else {
            score_general_mg(p, i_tally, i*tally.scores_.size(), filter_index,
              filter_weight, i_nuclide, atom_density, flux);
          }
        }
      }

    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : p.filter_matches_)
    match.bins_present_ = false;
}

void score_collision_tally(Particle& p)
{
  // Determine the collision estimate of the flux
  double flux = 0.0;
  if (p.type_ == Particle::Type::neutron || p.type_ == Particle::Type::photon) {
    if (!settings::survival_biasing) {
      flux = p.wgt_last_ / p.macro_xs_.total;
    } else {
      flux = (p.wgt_last_ + p.wgt_absorb_) / p.macro_xs_.total;
    }
  }

  for (auto i_tally : model::active_collision_tallies) {
    const Tally& tally {*model::tallies[i_tally]};

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p);
    auto end = FilterBinIter(tally, true, &p.filter_matches_);
    if (filter_iter == end) continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over nuclide bins.
      if (tally.all_nuclides_) {
        score_all_nuclides(p, i_tally, flux*filter_weight, filter_index,
          filter_weight);

      } else {
        for (auto i = 0; i < tally.nuclides_.size(); ++i) {
          auto i_nuclide = tally.nuclides_[i];

          double atom_density = 0.;
          if (i_nuclide >= 0) {
            auto j = model::materials[p.material_]->mat_nuclide_index_[i_nuclide];
            if (j == C_NONE) continue;
            atom_density = model::materials[p.material_]->atom_density_(j);
          }

          //TODO: consider replacing this "if" with pointers or templates
          if (settings::run_CE) {
            score_general_ce(p, i_tally, i*tally.scores_.size(), filter_index,
              filter_weight, i_nuclide, atom_density, flux);
          } else {
            score_general_mg(p, i_tally, i*tally.scores_.size(), filter_index,
              filter_weight, i_nuclide, atom_density, flux);
          }
        }
      }
    }

    // If the user has specified that we can assume all tallies are spatially
    // separate, this implies that once a tally has been scored to, we needn't
    // check the others. This cuts down on overhead when there are many
    // tallies specified
    if (settings::assume_separate) break;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : p.filter_matches_)
    match.bins_present_ = false;
}

void
score_surface_tally(Particle& p, const std::vector<int>& tallies)
{
  // No collision, so no weight change when survival biasing
  double flux = p.wgt_;

  for (auto i_tally : tallies) {
    auto& tally {*model::tallies[i_tally]};

    // Initialize an iterator over valid filter bin combinations.  If there are
    // no valid combinations, use a continue statement to ensure we skip the
    // assume_separate break below.
    auto filter_iter = FilterBinIter(tally, p);
    auto end = FilterBinIter(tally, true, &p.filter_matches_);
    if (filter_iter == end) continue;

    // Loop over filter bins.
    for (; filter_iter != end; ++filter_iter) {
      auto filter_index = filter_iter.index_;
      auto filter_weight = filter_iter.weight_;

      // Loop over scores.
      // There is only one score type for current tallies so there is no need
      // for a further scoring function.
      double score = flux * filter_weight;
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
    if (settings::assume_separate) break;
  }

  // Reset all the filter matches for the next tally event.
  for (auto& match : p.filter_matches_)
    match.bins_present_ = false;
}

} // namespace openmc
