#include "openmc/tallies/tally_scoring.h"
#include <csignal>
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
#include "openmc/distribution_multi.h"
#include "openmc/secondary_uncorrelated.h"
#include "openmc/geometry.h"
//#include "/home/open_mc/openmc/src/tallies/MyCalcs.cpp"
#include <typeinfo>
#include <string>
int ghost_counter=0;
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
        double xs = (score_bin == COHERENT)
                      ? micro.coherent
                      : (score_bin == INCOHERENT) ? micro.incoherent
                                                  : (score_bin == PHOTOELECTRIC)
                                                      ? micro.photoelectric
                                                      : micro.pair_production;
        score = xs * atom_density * flux;
      } else {
        double xs = (score_bin == COHERENT)
                      ? p.macro_xs().coherent
                      : (score_bin == INCOHERENT)
                          ? p.macro_xs().incoherent
                          : (score_bin == PHOTOELECTRIC)
                              ? p.macro_xs().photoelectric
                              : p.macro_xs().pair_production;
        score = xs * flux;
      }
      break;

    case HEATING:
      if (p.type() == Type::neutron) {
        score = score_neutron_heating(
          p, tally, flux, HEATING, i_nuclide, atom_density);
      } else {
        // The energy deposited is the difference between the pre-collision and
        // post-collision energy...
        score = E - p.E();

        // ...less the energy of any secondary particles since they will be
        // transported individually later
        const auto& bank = p.secondary_bank();
        for (auto it = bank.end() - p.n_bank_second(); it < bank.end(); ++it) {
          score -= it->E;
        }

        score *= p.wgt_last();
      }
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
      if (p.type() == Type::neutron) {
        score = score_neutron_heating(
          p, tally, flux, HEATING, i_nuclide, atom_density);
      } else {
        // The energy deposited is the difference between the pre-collision and
        // post-collision energy...
        score = E - p.E();

        // ...less the energy of any secondary particles since they will be
        // transported individually later
        const auto& bank = p.secondary_bank();
        for (auto it = bank.end() - p.n_bank_second(); it < bank.end(); ++it) {
          score -= it->E;
        }

        score *= p.wgt_last();
      }
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
  macro_xs.set_angle_index(p_u);
  if (i_nuclide >= 0) {
    nuc_xs.set_temperature_index(p.sqrtkT());
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
          score *= atom_density * nuc_xs.get_xs(MgxsType::TOTAL, p_g) /
                   macro_xs.get_xs(MgxsType::TOTAL, p_g);
        }
      } else {
        if (i_nuclide >= 0) {
          score = atom_density * flux * nuc_xs.get_xs(MgxsType::TOTAL, p_g);
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
          score *= nuc_xs.get_xs(MgxsType::INVERSE_VELOCITY, p_g) /
                   macro_xs.get_xs(MgxsType::TOTAL, p_g);
        } else {
          score *= macro_xs.get_xs(MgxsType::INVERSE_VELOCITY, p_g) /
                   macro_xs.get_xs(MgxsType::TOTAL, p_g);
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
        if (p.event() != TallyEvent::SCATTER)
          continue;
        // Since only scattering events make it here, again we can use the
        // weight entering the collision as the estimator for the reaction rate
        score = (p.wgt_last() - wgt_absorb) * flux;
        if (i_nuclide >= 0) {
          score *= atom_density *
                   nuc_xs.get_xs(MgxsType::SCATTER_FMU, p.g_last(), &p.g(),
                     &p.mu(), nullptr) /
                   macro_xs.get_xs(MgxsType::SCATTER_FMU, p.g_last(), &p.g(),
                     &p.mu(), nullptr);
        }
      } else {
        if (i_nuclide >= 0) {
          score =
            atom_density * flux *
            nuc_xs.get_xs(MgxsType::SCATTER, p_g, nullptr, &p.mu(), nullptr);
        } else {
          score = flux * macro_xs.get_xs(
                           MgxsType::SCATTER, p_g, nullptr, &p.mu(), nullptr);
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
                     &p.mu(), nullptr) /
                   macro_xs.get_xs(MgxsType::NU_SCATTER_FMU, p.g_last(), &p.g(),
                     &p.mu(), nullptr);
        }
      } else {
        if (i_nuclide >= 0) {
          score =
            atom_density * flux * nuc_xs.get_xs(MgxsType::NU_SCATTER, p_g);
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
          score *= atom_density * nuc_xs.get_xs(MgxsType::ABSORPTION, p_g) /
                   macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
        }
      } else {
        if (i_nuclide >= 0) {
          score =
            atom_density * flux * nuc_xs.get_xs(MgxsType::ABSORPTION, p_g);
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
          score *= atom_density * nuc_xs.get_xs(MgxsType::FISSION, p_g) /
                   macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
        } else {
          score *= macro_xs.get_xs(MgxsType::FISSION, p_g) /
                   macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
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
            score *= atom_density * nuc_xs.get_xs(MgxsType::NU_FISSION, p_g) /
                     macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
          } else {
            score *= macro_xs.get_xs(MgxsType::NU_FISSION, p_g) /
                     macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
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
            score *= atom_density * nuc_xs.get_xs(MgxsType::FISSION, p_g) /
                     macro_xs.get_xs(MgxsType::FISSION, p_g);
          }
        }
      } else {
        if (i_nuclide >= 0) {
          score =
            atom_density * flux * nuc_xs.get_xs(MgxsType::NU_FISSION, p_g);
        } else {
          score = flux * macro_xs.get_xs(MgxsType::NU_FISSION, p_g);
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
            score *= atom_density *
                     nuc_xs.get_xs(MgxsType::PROMPT_NU_FISSION, p_g) /
                     macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
          } else {
            score *= macro_xs.get_xs(MgxsType::PROMPT_NU_FISSION, p_g) /
                     macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
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
            score *= atom_density * nuc_xs.get_xs(MgxsType::FISSION, p_g) /
                     macro_xs.get_xs(MgxsType::FISSION, p_g);
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
          double abs_xs = macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
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
                             nullptr, nullptr, &d) /
                           abs_xs;
                } else {
                  score *= macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d) /
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
                score *=
                  nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g) / abs_xs;
              } else {
                score *=
                  macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g) / abs_xs;
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
                score *= atom_density * nuc_xs.get_xs(MgxsType::FISSION, p_g) /
                         macro_xs.get_xs(MgxsType::FISSION, p_g);
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
              score *= atom_density * nuc_xs.get_xs(MgxsType::FISSION, p_g) /
                       macro_xs.get_xs(MgxsType::FISSION, p_g);
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
                        nullptr, &d);
            } else {
              score = flux * macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                               nullptr, nullptr, &d);
            }
            score_fission_delayed_dg(
              i_tally, d_bin, score, score_index, p.filter_matches());
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
              const DelayedGroupFilter& filt {
                *dynamic_cast<DelayedGroupFilter*>(
                  model::tally_filters[i_dg_filt].get())};
              // Tally each delayed group bin individually
              for (auto d_bin = 0; d_bin < filt.n_bins(); ++d_bin) {
                auto d = filt.groups()[d_bin] - 1;
                score = wgt_absorb * flux;
                if (i_nuclide >= 0) {
                  score *= nuc_xs.get_xs(
                             MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d) *
                           nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d) /
                           abs_xs;
                } else {
                  score *= macro_xs.get_xs(
                             MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d) *
                           macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d) /
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
                           nuc_xs.get_xs(
                             MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d) *
                           nuc_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d) /
                           abs_xs;
                } else {
                  score += wgt_absorb * flux *
                           macro_xs.get_xs(
                             MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d) *
                           macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                             nullptr, nullptr, &d) /
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
                score += simulation::keff * atom_density * bank.wgt * flux *
                         nuc_xs.get_xs(
                           MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d) *
                         nuc_xs.get_xs(MgxsType::FISSION, p_g) /
                         macro_xs.get_xs(MgxsType::FISSION, p_g);
              } else {
                score += simulation::keff * bank.wgt * flux *
                         macro_xs.get_xs(
                           MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d);
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
              score =
                atom_density * flux *
                nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d) *
                nuc_xs.get_xs(
                  MgxsType::DELAYED_NU_FISSION, p_g, nullptr, nullptr, &d);
            } else {
              score = flux *
                      macro_xs.get_xs(
                        MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d) *
                      macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                        nullptr, nullptr, &d);
            }
            score_fission_delayed_dg(
              i_tally, d_bin, score, score_index, p.filter_matches());
          }
          continue;
        } else {
          score = 0.;
          for (auto d = 0; d < data::mg.num_delayed_groups_; ++d) {
            if (i_nuclide >= 0) {
              score +=
                atom_density * flux *
                nuc_xs.get_xs(MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d) *
                nuc_xs.get_xs(
                  MgxsType::DELAYED_NU_FISSION, p_g, nullptr, nullptr, &d);
            } else {
              score += flux *
                       macro_xs.get_xs(
                         MgxsType::DECAY_RATE, p_g, nullptr, nullptr, &d) *
                       macro_xs.get_xs(MgxsType::DELAYED_NU_FISSION, p_g,
                         nullptr, nullptr, &d);
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
          score *= atom_density * nuc_xs.get_xs(MgxsType::KAPPA_FISSION, p_g) /
                   macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
        } else {
          score *= macro_xs.get_xs(MgxsType::KAPPA_FISSION, p_g) /
                   macro_xs.get_xs(MgxsType::ABSORPTION, p_g);
        }
      } else {
        if (i_nuclide >= 0) {
          score =
            atom_density * flux * nuc_xs.get_xs(MgxsType::KAPPA_FISSION, p_g);
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
            auto j =
              model::materials[p.material()]->mat_nuclide_index_[i_nuclide];
            if (j == C_NONE)
              continue;
            atom_density = model::materials[p.material()]->atom_density_(j);
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
  fmt::print("collision tally called\n");
  // Determine the collision estimate of the flux
  double flux = 0.0;
  if (p.type() == ParticleType::neutron || p.type() == ParticleType::photon) {
    flux = p.wgt_last() / p.macro_xs().total;
  }

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
          auto j =
            model::materials[p.material()]->mat_nuclide_index_[i_nuclide];
          if (j == C_NONE)
            continue;
          atom_density = model::materials[p.material()]->atom_density_(j);
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
Position GetRotVector(double phi ,Position u_lab ,Position k )
  {
    return u_lab*std::cos(phi) + k.crossProduct(u_lab) * std::sin(phi) + k * k.dot(u_lab) * (1 - std::cos(phi));
  }

void score_point_tally(Particle& p)
{
  col_counter ++;
  fmt::print("------------------------collison happened------------------------\n");
  fmt::print("col counter = {}\n",col_counter);
  //std::cout << "mass in ev  " << p.getMass() << std::endl ;
  // Determine the collision estimate of the flux
  bool verbose=false;//true;
  double ReturnArray[4]= {std::nan(""),std::nan(""),std::nan(""),std::nan("")};
  double flux = 0.0;
  double flux1 = 0.0;
  double flux2 = 0.0;
  const auto& nuc {data::nuclides[p.event_nuclide()]};
  double awr = nuc->awr_;
  double dl = 0;
  getMu_COM(0,0,0,p,awr,ReturnArray , 0, dl);
  const auto& rx {nuc->reactions_[0]};
  auto& d = rx->products_[0].distribution_[0];
  auto d_ = dynamic_cast<UncorrelatedAngleEnergy*>(d.get());
  Direction u_lab {0-p.r().x,0-p.r().y,0-p.r().z};
  u_lab = u_lab/u_lab.norm();
  double cos_lab = u_lab.dot(p.u_last());
  double theta_lab = std::acos(cos_lab);
  double sin_lab = std::sin(theta_lab);
  double mu_COM_or_itay = (std::sqrt(awr*awr - sin_lab*sin_lab)*cos_lab - sin_lab*sin_lab)/awr;
  double mu_COM1 = ReturnArray[0];
  double E1 = ReturnArray[2];
  double mu_COM2 = ReturnArray[1];
  double E2 = ReturnArray[3];
  double ReturnArrayPlus[4] = {std::nan(""),std::nan(""),std::nan(""),std::nan("")};
  double ReturnArrayMinus[4] = {std::nan(""),std::nan(""),std::nan(""),std::nan("")};
// calculate new detector place
  // axis of rotation
  // cross product of incoming and outgoing secondary - for the plane of the collison
  Position k = ( p.u_last() ).crossProduct(u_lab);
  
  // santity for cross product 
  
 
   
 // normalize k
  k = k/k.norm();
  double dphi = 0.00001;
  Position u_lab_plus = GetRotVector(dphi,u_lab ,k);
  Position u_lab_minus = GetRotVector(-dphi,u_lab ,k);
 // now we rotate ulab - which is the vector to the dectctor in the plane of the collision by angle phi
  /*
 fmt::print("santity for cross product{0} \n",1);
  fmt::print("u last= {0} , {1} , {2}\n",p.u_last().x,p.u_last().y,p.u_last().z);
  fmt::print("u lab= {0} , {1} , {2}\n",u_lab.x,u_lab.y,u_lab.z);
  fmt::print("k= {0} , {1} , {2}\n",k.x,k.y,k.z);
  fmt::print("santity for Rotation{0} \n",1);
  fmt::print("u last= {0} , {1} , {2}\n",p.u_last().x,p.u_last().y,p.u_last().z);
  fmt::print("u lab original = {0} , {1} , {2}\n",u_lab.x,u_lab.y,u_lab.z);
  fmt::print("u lab plus= {0} , {1} , {2}\n",u_lab_plus.x,u_lab_plus.y,u_lab_plus.z);
  fmt::print("u lab minus= {0} , {1} , {2}\n",u_lab_minus.x,u_lab_minus.y,u_lab_minus.z);
  fmt::print("det plus = {0} , {1} , {2}\n",(u_lab_plus + p.r()).x,(u_lab_plus + p.r()).y,(u_lab_plus + p.r()).z);
  fmt::print("det minus = {0} , {1} , {2}\n",(u_lab_minus + p.r()).x,(u_lab_minus + p.r()).y,(u_lab_minus + p.r()).z);
  fmt::print("theta original , plus , minus = {0} , {1} , {2}\n",theta_lab, std::acos(u_lab_plus.dot(p.u_last())) ,std::acos(u_lab_minus.dot(p.u_last())));
*/
  
// Rodrigues' Rotation Formula

  getMu_COM((u_lab_plus + p.r()).x,(u_lab_plus + p.r()).y,(u_lab_plus + p.r()).z,p,awr,ReturnArrayPlus ,100 , dl);
  getMu_COM((u_lab_minus + p.r()).x,(u_lab_minus + p.r()).y,(u_lab_minus + p.r()).z,p,awr,ReturnArrayMinus ,-100 , dl );
  double MuPlus1 = ReturnArrayPlus[0]; // not sure about changing mu lab correctly
  double MuPlus2 = ReturnArrayPlus[1];
  double MuMinus1 = ReturnArrayMinus[0];
  double MuMinus2 = ReturnArrayMinus[1];
  double theta_pdf1 = d_->angle().get_pdf_value(p.E_last(),mu_COM1,p.current_seed());
  double theta_pdf2 = d_->angle().get_pdf_value(p.E_last(),mu_COM2,p.current_seed());
  double derivative1 = std::abs((MuPlus1-MuMinus1)/(2*dphi)/sin_lab);
  double derivative2 = std::abs((MuPlus2-MuMinus2)/(2*dphi)/sin_lab);   // divide by zero can cause nan!
  /*
  // one sided derivative
  double derivative1 = std::abs(MuPlus1-mu_COM1)/(dphi);
  double derivative2 = std::abs(MuPlus2-mu_COM2)/(dphi);
  */ 

  //double derivative =1;
  double theta_pdf_lab1 = theta_pdf1 * derivative1;
  double theta_pdf_lab2 = theta_pdf2 * derivative2;
  //double E_ghost = p.E_last()*(1+awr*awr+2*awr*mu_COM)/(1+awr)/(1+awr);
  
  //double m_incoming =MASS_NEUTRON_EV;
  double E_ghost1 = E1;
  double E_ghost2 = E2;


// first solution
  if(!std::isnan(mu_COM1))
  {
  Particle ghost_particle=Particle();
  ghost_particle.initilze_ghost_particle(p,u_lab,E_ghost1);
  ghost_counter++;
  //fmt::print("---------------------ghost particle created {}---------------------\n",ghost_counter);



  //calculate shilding
  double total_distance = std::sqrt(p.r().x*p.r().x+p.r().y*p.r().y+p.r().z*p.r().z);
  double remaining_distance = std::sqrt(p.r().x*p.r().x+p.r().y*p.r().y+p.r().z*p.r().z);
  double total_MFP = 0;
  ghost_particle.event_calculate_xs();
  ghost_particle.boundary() = distance_to_boundary(ghost_particle);
  double advance_distance = ghost_particle.boundary().distance;
  while(advance_distance<remaining_distance)
  {
    total_MFP += advance_distance*ghost_particle.macro_xs().total;
      //Advance particle in space and time
    for (int j = 0; j < ghost_particle.n_coord(); ++j) {
      ghost_particle.coord(j).r += advance_distance * ghost_particle.coord(j).u;
    }
      remaining_distance-=advance_distance;
      //fmt::print("advane distance 1 ={} \n",advance_distance);
      //fmt::print("advane XS 1= {}\n",ghost_particle.macro_xs().total);
      //fmt::print("ghost particle speed= {}\n",ghost_particle.speed());
      ghost_particle.time() += advance_distance / ghost_particle.speed();
      ghost_particle.event_cross_surface();
      //fmt::print("pos = {0} , {1} , {2}\n",ghost_particle.r().x,ghost_particle.r().y,ghost_particle.r().z);
      ghost_particle.event_calculate_xs();
      ghost_particle.boundary() = distance_to_boundary(ghost_particle);
      advance_distance = ghost_particle.boundary().distance;
      //fmt::print("advane distance 2 ={} \n",advance_distance);
      //fmt::print("advane XS 2= {}\n",ghost_particle.macro_xs().total);

  }
  
  //// kill patricle
  //ghost::ghost_particles[0].event_death();
  //ghost::ghost_particles.resize(0);
  if(!std::isnan(flux1))
  {
  flux1 = ghost_particle.wgt()*exp(-total_MFP)/(2*3.14*total_distance*total_distance)*theta_pdf_lab1;
  }
  if(std::isnan(flux1))
  {
    flux1 = 0;
  }
  if(verbose)
  {
    fmt::print("------------------------flux contribution------------------------\n");
    fmt::print("parent particle pos = {0} , {1} , {2}\n",p.r().x,p.r().y,p.r().z);
    fmt::print("parent particle u = {0} , {1} , {2}\n",p.u_last().x,p.u_last().y,p.u_last().z);
    fmt::print("parent particle E = {}\n",p.E_last());
    fmt::print("ghost particle E = {}\n",ghost_particle.E());
    fmt::print("ghost particle u = {0} , {1} , {2}\n",ghost_particle.u().x,ghost_particle.u().y,ghost_particle.u().z);
    fmt::print("ghost particle Mu COM 1 Arik= {}\n",mu_COM1);
    fmt::print("ghost particle Mu COM 2 Arik= {}\n",mu_COM2);
    fmt::print("ghost particle Mu COM Or and Itay = {}\n",mu_COM_or_itay);
    fmt::print("flux1 = {}\n",flux1);
    fmt::print("theta1 cm pdf ={} \n",theta_pdf1);
    fmt::print("theta2 cm pdf ={} \n",theta_pdf2);
    fmt::print("MuPlus for sol 0 ={} \n",MuPlus1);
    fmt::print("MuPlus for sol 1 ={}  \n",MuPlus2);
    fmt::print("MuMinus for sol 0 ={} \n",MuMinus1);
    fmt::print("MuMinus for sol 1 ={}  \n",MuMinus2);
    fmt::print("derivative1 ={} \n",derivative1);
    fmt::print("derivative2 ={} \n",derivative2);
    fmt::print("theta1 lab pdf ={} \n",theta_pdf_lab1);
    fmt::print("theta2 lab pdf ={} \n",theta_pdf_lab2);
    fmt::print("santity for cross product{0} \n",1);
  fmt::print("u last= {0} , {1} , {2}\n",p.u_last().x,p.u_last().y,p.u_last().z);
  fmt::print("u lab= {0} , {1} , {2}\n",u_lab.x,u_lab.y,u_lab.z);
  fmt::print("k= {0} , {1} , {2}\n",k.x,k.y,k.z);
  fmt::print("santity for Rotation{0} \n",1);
  fmt::print("u last= {0} , {1} , {2}\n",p.u_last().x,p.u_last().y,p.u_last().z);
  fmt::print("u lab original = {0} , {1} , {2}\n",u_lab.x,u_lab.y,u_lab.z);
  fmt::print("u lab plus= {0} , {1} , {2}\n",u_lab_plus.x,u_lab_plus.y,u_lab_plus.z);
  fmt::print("u lab minus= {0} , {1} , {2}\n",u_lab_minus.x,u_lab_minus.y,u_lab_minus.z);
  fmt::print("det plus = {0} , {1} , {2}\n",(u_lab_plus + p.r()).x,(u_lab_plus + p.r()).y,(u_lab_plus + p.r()).z);
  fmt::print("det minus = {0} , {1} , {2}\n",(u_lab_minus + p.r()).x,(u_lab_minus + p.r()).y,(u_lab_minus + p.r()).z);
  fmt::print("theta original , plus , minus = {0} , {1} , {2}\n",theta_lab, std::acos(u_lab_plus.dot(p.u_last())) ,std::acos(u_lab_minus.dot(p.u_last())));


    //fmt::print("Mu - ={} \n",MuMinus);
    //fmt::print("Mu + ={} \n",MuPlus);
  }
  }

  // second soultion
   if(!std::isnan(mu_COM2))
  {
  Particle ghost_particle=Particle();
  ghost_particle.initilze_ghost_particle(p,u_lab,E_ghost2);
  ghost_counter++;
  //fmt::print("---------------------ghost particle created {}---------------------\n",ghost_counter);



  //calculate shilding
  double total_distance = std::sqrt(p.r().x*p.r().x+p.r().y*p.r().y+p.r().z*p.r().z);
  double remaining_distance = std::sqrt(p.r().x*p.r().x+p.r().y*p.r().y+p.r().z*p.r().z);
  double total_MFP = 0;
  ghost_particle.event_calculate_xs();
  ghost_particle.boundary() = distance_to_boundary(ghost_particle);
  double advance_distance = ghost_particle.boundary().distance;
  while(advance_distance<remaining_distance)
  {
    total_MFP += advance_distance*ghost_particle.macro_xs().total;
      //Advance particle in space and time
    for (int j = 0; j < ghost_particle.n_coord(); ++j) {
      ghost_particle.coord(j).r += advance_distance * ghost_particle.coord(j).u;
    }
      remaining_distance-=advance_distance;
      //fmt::print("advane distance 1 ={} \n",advance_distance);
      //fmt::print("advane XS 1= {}\n",ghost_particle.macro_xs().total);
      //fmt::print("ghost particle speed= {}\n",ghost_particle.speed());
      ghost_particle.time() += advance_distance / ghost_particle.speed();
      ghost_particle.event_cross_surface();
      //fmt::print("pos = {0} , {1} , {2}\n",ghost_particle.r().x,ghost_particle.r().y,ghost_particle.r().z);
      ghost_particle.event_calculate_xs();
      ghost_particle.boundary() = distance_to_boundary(ghost_particle);
      advance_distance = ghost_particle.boundary().distance;
      //fmt::print("advane distance 2 ={} \n",advance_distance);
      //fmt::print("advane XS 2= {}\n",ghost_particle.macro_xs().total);

  }
  
  //// kill patricle
  //ghost::ghost_particles[0].event_death();
  //ghost::ghost_particles.resize(0);
if(!std::isnan(flux2))
  {
  flux2 = ghost_particle.wgt()*exp(-total_MFP)/(2*3.14*total_distance*total_distance)*theta_pdf_lab2;
  }
  if(std::isnan(flux2))
  {
    flux2 = 0;
  }
  if(verbose)
  {
    fmt::print("flux params = {0} {1} {2} {3}\n",ghost_particle.wgt() , exp(-total_MFP) , (2*3.14*total_distance*total_distance) ,theta_pdf_lab2 );

    fmt::print("------------------------flux contribution------------------------\n");
    fmt::print("parent particle pos = {0} , {1} , {2}\n",p.r().x,p.r().y,p.r().z);
    fmt::print("parent particle u = {0} , {1} , {2}\n",p.u_last().x,p.u_last().y,p.u_last().z);
    fmt::print("parent particle E = {}\n",p.E_last());
    fmt::print("ghost particle E = {}\n",ghost_particle.E());
    fmt::print("ghost particle u = {0} , {1} , {2}\n",ghost_particle.u().x,ghost_particle.u().y,ghost_particle.u().z);
    fmt::print("ghost particle Mu COM 1 Arik= {}\n",mu_COM1);
    fmt::print("ghost particle Mu COM 2 Arik= {}\n",mu_COM2);
    fmt::print("ghost particle Mu COM Or and Itay = {}\n",mu_COM_or_itay);
    fmt::print("flux2 = {}\n",flux2);
    fmt::print("theta1 cm pdf ={} \n",theta_pdf1);
    fmt::print("theta2 cm pdf ={} \n",theta_pdf2);
    fmt::print("derivative1 ={} \n",derivative1);
    fmt::print("derivative2 ={} \n",derivative2);
    fmt::print("theta1 lab pdf ={} \n",theta_pdf_lab1);
    fmt::print("theta2 lab pdf ={} \n",theta_pdf_lab2);

     fmt::print("MuPlus for sol 0 ={} \n",MuPlus1);
    fmt::print("MuPlus for sol 1 ={}  \n",MuPlus2);

    fmt::print("MuMinus for sol 0 ={} \n",MuMinus1);
    fmt::print("MuMinus for sol 1 ={}  \n",MuMinus2);


       fmt::print("santity for cross product{0} \n",1);
  fmt::print("u last= {0} , {1} , {2}\n",p.u_last().x,p.u_last().y,p.u_last().z);
  fmt::print("u lab= {0} , {1} , {2}\n",u_lab.x,u_lab.y,u_lab.z);
  fmt::print("k= {0} , {1} , {2}\n",k.x,k.y,k.z);
  fmt::print("santity for Rotation{0} \n",1);
  fmt::print("u last= {0} , {1} , {2}\n",p.u_last().x,p.u_last().y,p.u_last().z);
  fmt::print("u lab original = {0} , {1} , {2}\n",u_lab.x,u_lab.y,u_lab.z);
  fmt::print("u lab plus= {0} , {1} , {2}\n",u_lab_plus.x,u_lab_plus.y,u_lab_plus.z);
  fmt::print("u lab minus= {0} , {1} , {2}\n",u_lab_minus.x,u_lab_minus.y,u_lab_minus.z);
  fmt::print("det plus = {0} , {1} , {2}\n",(u_lab_plus + p.r()).x,(u_lab_plus + p.r()).y,(u_lab_plus + p.r()).z);
  fmt::print("det minus = {0} , {1} , {2}\n",(u_lab_minus + p.r()).x,(u_lab_minus + p.r()).y,(u_lab_minus + p.r()).z);
  fmt::print("theta original , plus , minus = {0} , {1} , {2}\n",theta_lab, std::acos(u_lab_plus.dot(p.u_last())) ,std::acos(u_lab_minus.dot(p.u_last())));
  fmt::print("p weight = {0} last ={1}\n",p.wgt(),p.wgt_last());
//fmt::print("Mu - ={} \n",MuMinus);
    //fmt::print("Mu + ={} \n",MuPlus);
  }
  }


  flux = flux1 + flux2;

  if (p.type() != ParticleType::neutron) {
    if(p.event_mt() != 2)
    {
      flux = 0;
    }
  }


  for (auto i_tally : model::active_point_tallies) {
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


/*


void getMu_COM(double x_det , double y_det , double z_det ,Particle p_col , double awr , double incoming_mass ,double ReturnArray[],int diff_mode,double dl )
{
  //double inv_gamma = p_col.getMass()/ (p_col.E_last() + p_col.getMass());
  //double p_momentum =  (1/inv_gamma)*(p_col.getMass() / 1e9 )*std::sqrt(1 - inv_gamma * inv_gamma); // GeV/c
  //struct retVals { int i1, i2; }; // Declare a local structure
  double m1= p_col.getMass()/1e6; // mass of incoming particle in MeV
  double E1_tot = p_col.E_last()/1e6 + m1; // total Energy of incoming particle in MeV
  double p1_tot = std::sqrt(E1_tot*E1_tot  - m1*m1);
  //return retVals {10, 20};
  // std::cout<<"  p col "<<p_col<<std::endl;
  // double r1[4];

  
  //r1[0] = 0;
  //r1[1] = x_det;
  //r1[2] = y_det;
  //r1[3] = z_det;

  //for (int i = 0; i < 4; i++)
   // {std::cout << r1[i];}
  double r1[4]  = {0, x_det, y_det, z_det};  // detector position lab {ignor, x, y, z}
  double r2[4]= {0, p_col.r().x, p_col.r().y, p_col.r().z}; // collision position lab {ignor, x, y, z} 
  double r3[4]; // r1-r2 vector from collision to detector
  double m2= m1*awr; // mass of target matirial
  std::cout <<"m2 = "<<m2<<std::endl;
  double m3= m1; // mass of outgoing particle to detector //
  double m4= m2; // mass of recoil target  system
  double p1[3]={p1_tot*p_col.u_last().x,p1_tot* p_col.u_last().y,p1_tot* p_col.u_last().z}; // 3 momentum of incoming particle
  double p2[3]={0, 0, 0}; //3 momentum of target in lab // need seed from openmc
  for (int i = 0; i < 4; i++)
    {std::cout <<"r1 "<<i<<" "<< r1[i] <<"  ";}
  std::cout <<std::endl;
  //std::raise(SIGTRAP);
  for (int i = 0; i < 4; i++)
    {std::cout <<"r2 "<<i<<" "<< r2[i] <<"  ";}
  std::cout <<std::endl;


  
    // calculate
  double Fp1[4]; //four momentum of incoming particle  
  double Fp2[4]; //four momentum of target
  double Fp3[4]; //four momentum of particle going to the detector
  double UFp3[4]; //unit 3  momentum of particle going to the detector
  double CM[4]; //four momentum of center of mass frame in lab
  double LAB[4]; //four momentum of lab in center of mass frame 
  double mCM; // mass of center of mass system
  double pCM; // momentum of center of mass system
  double p3LAB; // momentum of out particle
  
  double CMFp1[4]; //four momentum of incoming particle in CM  
  double CMFp3[4]; //four momentum of out particle in CM  
  double CMFp4[4]; //four momentum of out target in CM  
  double Fp4[4]; //four momentum of out target in LAB  
  double CMUp1[4]; //unit three vector of incoming particle in CM  
  double CMUp3[4]; //unit three vector of out particle in CM  
  double cosCM; // cosine of out going particle to detector in CM frame
  double cosLAB; // cosine of out going particle to detector in LAB frame
  double CME3; // Energy of out particle in CM
  double CMp3; // momentum of out particle in CM
  double Ur3[4], UCM[4];
  double aa, bb, cc;
  //double ReturnArray[2];
  
  

  Fp1[0]=0;
  Fp2[0]=0;
  CM[0]=0;
  LAB[0]=0;
  r3[0]=0;


  for(int i=0; i<3; i++){
    Fp1[i+1]=p1[i];
    Fp2[i+1]=p2[i];
    CM[i+1]=Fp1[i+1]+Fp2[i+1];
    LAB[i+1]=-CM[i+1];
    r3[i+1]=r1[i+1]-r2[i+1];
  }
 
  Fp1[0]=sqrt(Fp1[1]*Fp1[1]+Fp1[2]*Fp1[2]+Fp1[3]*Fp1[3]+m1*m1);
  Fp2[0]=sqrt(Fp2[1]*Fp2[1]+Fp2[2]*Fp2[2]+Fp2[3]*Fp2[3]+m2*m2);
  CM[0]=Fp1[0]+Fp2[0];
  LAB[0]=CM[0];
  r3[0]=0;

  mCM=sqrt(CM[0]*CM[0]-CM[1]*CM[1]-CM[2]*CM[2]-CM[3]*CM[3]);
  pCM=sqrt(CM[1]*CM[1]+CM[2]*CM[2]+CM[3]*CM[3]);
  CME3=(mCM*mCM-m4*m4+m3*m3)/(2*mCM); // energy of out going particle  in CM
  CMp3=sqrt(CME3*CME3-m3*m3);
  
  std::cout<<mCM<<"  "<<pCM<<"  "<<CME3<<std::endl;

  boostf(CM,Fp1,CMFp1);
  double diff = 1;
  Vunit(CM,UCM); 
  Vunit(r3,Ur3); 
  cosLAB=Vdot(UCM,Ur3);
  if (diff_mode == 1 && cosLAB+dl<=1)
     {cosLAB = cosLAB+dl;}
  if (diff_mode == -1 && cosLAB-dl>=-1)
     {cosLAB = cosLAB-dl;} 
  std::cout<<"cosLAB=  "<<cosLAB<<std::endl;
  Vunit(CMFp1,CMUp1);


  aa=pCM*pCM*cosLAB*cosLAB-CM[0]*CM[0];
  bb=2*pCM*cosLAB*CME3*mCM;
  cc=CME3*mCM*CME3*mCM-m3*m3*CM[0]*CM[0];

  p3LAB=0;
  p3LAB=(-bb+sqrt(bb*bb-4*aa*cc))/2.0/aa;
  if(p3LAB<=0) p3LAB=(-bb-sqrt(bb*bb-4*aa*cc))/2/aa;
  if(p3LAB<=0) {
    std::cout<<" detector out of range" <<std::endl;
    //return -1;
  }


  std::cout<<"p3LAB= "<<p3LAB<<std::endl;

  Fp3[0]=sqrt(p3LAB*p3LAB+m3*m3);
  for(int i=0; i<3; i++){
    Fp3[i+1]=p3LAB*Ur3[i+1];
  }
  
  boostf(CM,Fp3,CMFp3);
  Vunit(CMFp3,CMUp3); 
  cosCM=Vdot(UCM,CMUp3);

  std::cout<<"cosCM= "<<cosCM<<std::endl;
  ReturnArray[0] = cosCM;
  ReturnArray[1] = (Fp3[0]-m3)*1e6;  //retrun Energy in eV
  

  std::cout<<" -----------  sanity tests ----------------------------------"<<std::endl;
 
 Vunit(Fp3,UFp3); 

  std::cout<<Vdot(UCM,UFp3)<<" = cosLAB =  "<<cosLAB<<std::endl;

  std::cout<<CMFp3[0]<<" m3 energy in CM  "<<CME3<<std::endl;

  CMFp4[0]=0;
  for(int i=0; i<3; i++){
    CMFp4[i+1]=-CMFp3[i+1];
  }
  CMFp4[0]=sqrt(m4*m4+CMFp4[1]*CMFp4[1]+CMFp4[2]*CMFp4[2]+CMFp4[3]*CMFp4[3]);
  boostf(LAB,CMFp4,Fp4);

  for(int i=0; i<4; i++){
    std::cout<<CM[i]<<" = energy momentum conserv in LAB = "<<Fp4[i]+Fp3[i]<<std::endl;
  }
 
  double test[4];
  boostf(LAB,CMFp1,test);
  for (int i=0; i<4; i++){
    std::cout<<test[i]<<std::endl;
    std::cout<<Fp1[i]<<std::endl;
  }
  //return retVals {cosCM, CME3};

}


void boostf( double A[4], double B[4], double X[4])
{
  //
  //     boosts B(labfram) to A rest frame and gives output in X
  //
  double W;
  int j;
  
  if ((A[0]*A[0]-A[1]*A[1]-A[2]*A[2]-A[3]*A[3])<=0) { 
    std::cout <<"negative sqrt in boostf"<<A[0]<<A[1]<<A[2]<<A[3]<<std::endl;}
      
  W=sqrt(A[0]*A[0]-A[1]*A[1]-A[2]*A[2]-A[3]*A[3]);

  if(W==0 || W==(-A[0])) std::cout <<"divid by 0 in boostf"<<std::endl;

  X[0]=(B[0]*A[0]-B[3]*A[3]-B[2]*A[2]-B[1]*A[1])/W;
    for(j=1; j<=3; j++) {
      X[j]=B[j]-A[j]*(B[0]+X[0])/(A[0]+W);
    } 
  return;
}




double Vdot(double A[4],double B[4])
{
  int j;

  double dot = 0;

  for(j=1; j<=3; j++) {
    dot = dot + A[j]*B[j];
  }

       
  return dot;   
}


void Vunit(double A[4] ,double B[4])
{
  double fff;
  int j;

  fff = 0;
  
  for(j=1; j<=3; j++) {
    fff = fff + A[j]*A[j];
  }
  
  if (fff==0) {
    std::cout <<"in vunit divid by zero" << std::endl;
    return;
  }
  
  for(j=1; j<=3; j++) {
    B[j] = A[j]/sqrt(fff);
  }
  B[0] = 0;
  
  return;   
}


void Vcros(double A[4],double B[4],double C[4])
{
  C[1] = A[2]*B[3]-A[3]*B[2];
  C[2] = A[3]*B[1]-A[1]*B[3];
  C[3] = A[1]*B[2]-A[2]*B[1];
  C[0] = 0;

  if (C[1]==0 && C[2]==0 && C[3]==0.) { 
    std::cout << "vcross zero" << std::endl;
  }
  
  return;   
}
*/


      //Author: Arik Kreisel


void getMu_COM(double x_det , double y_det , double z_det ,Particle p_col , double awr ,double ReturnArray[],int diff_mode,double dl )
{

  double cs=0;
  int m=0;
  
  double m1= p_col.getMass()/1e6; // mass of incoming particle in MeV
  double E1_tot = p_col.E_last()/1e6 + m1; // total Energy of incoming particle in MeV
  double p1_tot = std::sqrt(E1_tot*E1_tot  - m1*m1);
  double r1[4]  = {0, x_det, y_det, z_det};  // detector position lab {ignor, x, y, z}
  double r2[4]= {0, p_col.r().x, p_col.r().y, p_col.r().z}; // collision position lab {ignor, x, y, z} 
  double r3[4]; // r1-r2 vector from collision to detector
  double m2= m1*awr; // mass of target matirial
  double m3= m1; // mass of outgoing particle to detector  (rest mass?)
  double m4= m2; // mass of recoil target  system
  double p1[3]={p1_tot*p_col.u_last().x,p1_tot* p_col.u_last().y,p1_tot* p_col.u_last().z}; // 3 momentum of incoming particle
  double p2[3]={0, 0, 0}; //3 momentum of target in lab  
 
  // calculate
  double Fp1[4]; //four momentum of incoming particle  in LAB 
  double Fp2[4]; //four momentum of target in LAB
  double Fp3[4]; //four momentum of particle going to the detector in LAB
  double UFp3[4]; //unit 3  momentum of particle going to the detector
  double CM[4]; //four momentum of center of mass frame in lab
  double LAB[4]; //four momentum of lab in center of mass frame 
  double mCM; // mass of center of mass system
  double pCM; // momentum of center of mass system
  double p3LAB[2]; // momentum of out particle
  
  double CMFp1[4]; //four momentum of incoming particle in CM  
  double CMFp2[4]; //four momentum of incoming target in CM  
  double CMFp3[4]; //four momentum of out particle in CM  
  double CMFp4[4]; //four momentum of out target in CM  
  double Fp4[4]; //four momentum of out target in LAB  
  double CMUp1[4]; //unit three vector of incoming particle in CM  
  double CMUp3[4]; //unit three vector of out particle in CM  
  double cosCM; // cosine of out going particle to detector in CM frame
  double cosLAB; // cosine of out going particle to detector in LAB frame
  double CME3; // Energy of out particle in CM
  double CMp3; // momentum of out particle in CM
  double  Ur3[4], UCM[4], ULAB[4];
  double aa, bb, cc;


  Fp1[0]=0;
  Fp2[0]=0;
  CM[0]=0;
  LAB[0]=0;
  r3[0]=0;


  for(int i=0; i<3; i++){
    Fp1[i+1]=p1[i];
    Fp2[i+1]=p2[i];
    CM[i+1]=Fp1[i+1]+Fp2[i+1];
    LAB[i+1]=-CM[i+1];
    r3[i+1]=r1[i+1]-r2[i+1];
  }
 
  Fp1[0]=sqrt(Fp1[1]*Fp1[1]+Fp1[2]*Fp1[2]+Fp1[3]*Fp1[3]+m1*m1);
  Fp2[0]=sqrt(Fp2[1]*Fp2[1]+Fp2[2]*Fp2[2]+Fp2[3]*Fp2[3]+m2*m2);
  CM[0]=Fp1[0]+Fp2[0];
  LAB[0]=CM[0];
  r3[0]=0;

  mCM=sqrt(CM[0]*CM[0]-CM[1]*CM[1]-CM[2]*CM[2]-CM[3]*CM[3]);
  pCM=sqrt(CM[1]*CM[1]+CM[2]*CM[2]+CM[3]*CM[3]);
  CME3=(mCM*mCM-m4*m4+m3*m3)/(2*mCM); // energy of out going particle  in CM
  CMp3=sqrt(CME3*CME3-m3*m3);

  //std::cout<<" mCM= "<<mCM<<" pCM= "<<pCM<<" CME3= "<<CME3<<std::endl;

  boostf(CM,Fp1,CMFp1);
  boostf(CM,Fp2,CMFp2);

  Vunit(CM,UCM); 
  Vunit(LAB,ULAB); 
  Vunit(r3,Ur3); 
  cosLAB=Vdot(UCM,Ur3);
  if (diff_mode == 1 && cosLAB+dl<=1)
     {cosLAB = cosLAB+dl;}
  if (diff_mode == -1 && cosLAB-dl>=-1)
     {cosLAB = cosLAB-dl;} 

  //std::cout<<"cosLAB=  "<<cosLAB<<std::endl;
  Vunit(CMFp1,CMUp1);


  aa=pCM*pCM*cosLAB*cosLAB-CM[0]*CM[0];
  bb=2*pCM*cosLAB*CME3*mCM;
  cc=CME3*mCM*CME3*mCM-m3*m3*CM[0]*CM[0];

   int j=1;

   if(bb*bb-4*aa*cc<0) {
    //std::cout<<" detector out of range" <<std::endl;
    return; //continue;
  }
   p3LAB[0]=(-bb+sqrt(bb*bb-4*aa*cc))/2.0/aa;
   if(p3LAB[0]<=0) {p3LAB[0]=(-bb-sqrt(bb*bb-4*aa*cc))/2/aa;
     if(p3LAB[0]<=0) {
       std::cout<<" detector out of range" <<std::endl;
       return; //continue;
     }
   }
   else if(p3LAB[0]>0){p3LAB[1]=(-bb-sqrt(bb*bb-4*aa*cc))/2/aa;
     if(p3LAB[1]>0) j=j+1;
   }

  //std::cout<<" p3LAB1= "<<(-bb+sqrt(bb*bb-4*aa*cc))/2.0/aa<<" p3LAB2= "<<(-bb-sqrt(bb*bb-4*aa*cc))/2.0/aa<<std::endl;

  for (int l=0;l<j;l++){

    //std::cout<<"l= "<<l<<std::endl;

    Fp3[0]=sqrt(p3LAB[l]*p3LAB[l]+m3*m3);
    for(int i=0; i<3; i++){
      Fp3[i+1]=p3LAB[l]*Ur3[i+1];
    }
  
    boostf(CM,Fp3,CMFp3);
    Vunit(CMFp3,CMUp3); 
    cosCM=Vdot(UCM,CMUp3); // input to openMC Cross section calculation
    ReturnArray[l] = cosCM;
    ReturnArray[l+2] = (Fp3[0]-m3)*1e6;  //retrun Energy in eV
    /*
    
    std::cout<<" -----------  running diff mode: "<< diff_mode << " with dl = "<< dl <<" -----------------------------------"<<std::endl;
    std::cout<<" -----------  sanity tests for solution "<< l <<" -----------------------------------"<<std::endl;
    std::cout<<"cosCM= "<<cosCM<<std::endl;
    Vunit(Fp3,UFp3); 
    
    std::cout<<Vdot(UCM,UFp3)<<" = cosLAB =  "<<cosLAB<<std::endl;

    std::cout<<CMFp3[0]<<" m3 energy in CM  "<<CME3<<std::endl;
    
    CMFp4[0]=0;
    for(int i=0; i<3; i++){
      CMFp4[i+1]=-CMFp3[i+1];
    }
    CMFp4[0]=sqrt(m4*m4+CMFp4[1]*CMFp4[1]+CMFp4[2]*CMFp4[2]+CMFp4[3]*CMFp4[3]);
    boostf(LAB,CMFp4,Fp4);
  
    for(int i=0; i<4; i++){
      std::cout<<i<<"  "<<CM[i]<<" = energy momentum conserv in LAB = "<<Fp4[i]+Fp3[i]<<std::endl;
      std::cout<<i<<"  "<<CMFp1[i]+CMFp2[i]<<" = energy momentum conserv in CM = "<<CMFp4[i]+CMFp3[i]<<std::endl;
    }


 
    double test[4];
    boostf(LAB,CMFp3,test);
    for (int i=0; i<4; i++){
      std::cout<<i<<"  "<<test[i]<<std::endl;
      std::cout<<i<<"  "<<Fp3[i]<<std::endl;
    }
    
    std::cout <<" m3 enerrgy in lab ="<<Fp3[0]<<std::endl;
    std::cout <<" m4 enerrgy in lab ="<<Fp4[0]<<std::endl;
    std::cout <<" sqrt(E3^2-P3^2) in lab ="<<sqrt(Fp3[0]*Fp3[0]-Fp3[1]*Fp3[1]-Fp3[2]*Fp3[2]-Fp3[3]*Fp3[3])<<std::endl;
    std::cout <<" sqrt(E4^2-P4^2) in lab ="<<sqrt(Fp4[0]*Fp4[0]-Fp4[1]*Fp4[1]-Fp4[2]*Fp4[2]-Fp4[3]*Fp4[3])<<std::endl;
    std::cout <<" m1 enerrgy in lab ="<<Fp1[0]<<std::endl;

    double tanLab=sqrt(1-cosCM*cosCM)/(CM[0]/mCM*(pCM*CME3/CM[0]/CMp3+cosCM));
    double fCOSlab=cos(atan(tanLab));
    

    std::cout<<"cos( atan( tanLab= "<<cos(atan(tanLab))<<std::endl;
    */
  }
  
 }


void boostf( double A[4], double B[4], double X[4])
{
  //
  //     boosts B(labfram) to A rest frame and gives output in X
  //
  double W;
  int j;
  
  if ((A[0]*A[0]-A[1]*A[1]-A[2]*A[2]-A[3]*A[3])<=0) { 
    std::cout <<"negative sqrt in boostf"<<A[0]<<A[1]<<A[2]<<A[3]<<std::endl;}
      
  W=sqrt(A[0]*A[0]-A[1]*A[1]-A[2]*A[2]-A[3]*A[3]);

  if(W==0 || W==(-A[0])) std::cout <<"divid by 0 in boostf"<<std::endl;

  X[0]=(B[0]*A[0]-B[3]*A[3]-B[2]*A[2]-B[1]*A[1])/W;
    for(j=1; j<=3; j++) {
      X[j]=B[j]-A[j]*(B[0]+X[0])/(A[0]+W);
    } 

  return;
}




double Vdot(double A[4],double B[4])
{
  int j;

  double dot = 0;

  for(j=1; j<=3; j++) {
    dot = dot + A[j]*B[j];
  }

       
  return dot;   
}


void Vunit(double A[4] ,double B[4])
{
  double fff;
  int j;

  fff = 0;
  
  for(j=1; j<=3; j++) {
    fff = fff + A[j]*A[j];
  }
  
  if (fff==0) {
    std::cout <<"in vunit divid by zero" << std::endl;
    return;
  }
  
  for(j=1; j<=3; j++) {
    B[j] = A[j]/sqrt(fff);
  }
  B[0] = 0;
  
  return;   
}


void Vcros(double A[4],double B[4],double C[4])
{
  C[1] = A[2]*B[3]-A[3]*B[2];
  C[2] = A[3]*B[1]-A[1]*B[3];
  C[3] = A[1]*B[2]-A[2]*B[1];
  C[0] = 0;

  if (C[1]==0 && C[2]==0 && C[3]==0.) { 
    std::cout << "vcross zero" << std::endl;
  }
  
  return;   
}




} // namespace openmc
