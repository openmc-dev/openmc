#include "openmc/tallies/trigger.h"

#include <cmath>
#include <sstream>
#include <utility> // for std::pair

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/reaction.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/tallies/tally.h"

namespace openmc {

static std::pair<double, double>
get_tally_uncertainty(int i_tally, int score_index, int filter_index)
{
  int n;
  int err = openmc_tally_get_n_realizations(i_tally, &n);

  auto results = tally_results(i_tally);
  //TODO: off-by-one
  auto sum = results(filter_index-1, score_index-1, RESULT_SUM);
  auto sum_sq = results(filter_index-1, score_index-1, RESULT_SUM_SQ);

  auto mean = sum / n;
  double std_dev = std::sqrt((sum_sq/n - mean*mean) / (n-1));
  double rel_err = (mean != 0.) ? std_dev / std::abs(mean) : 0.;

  return {std_dev, rel_err};
}

// Functions defined in F90.
extern "C" int n_tally_triggers(int i_tally);
extern "C" Trigger* get_tally_trigger(int i_tally, int i_trig);

//! Finds the limiting limiting tally trigger.
//
//! param[out] ratio The uncertainty/threshold ratio for the most limiting
//!   tally trigger
//! param[out] tally_id The ID number of the most limiting tally
//! param[out] score The most limiting tally score bin

void
check_tally_triggers(double& ratio, int& tally_id, int& score)
{
  ratio = 0.;
  //TODO: off-by-one
  for (auto i_tally = 1; i_tally < model::tallies.size()+1; ++i_tally) {
    const Tally& t {*model::tallies[i_tally-1]};

    // Ignore tallies with less than two realizations.
    int n_reals;
    int err = openmc_tally_get_n_realizations(i_tally, &n_reals);
    if (n_reals < 2) continue;

    //TODO: off-by-one
    for (auto i_trig = 1; i_trig < n_tally_triggers(i_tally)+1; ++i_trig) {
      auto& trigger {*get_tally_trigger(i_tally, i_trig)};

      trigger.std_dev = 0.;
      trigger.rel_err = 0.;
      trigger.variance = 0.;

      const auto& results = tally_results(i_tally);
      //TODO: off-by-one
      for (auto filter_index = 1; filter_index < results.shape()[0]+1;
           ++filter_index) {
        //TODO: off-by-one
        for (auto score_index = 1; score_index < results.shape()[1]+1;
              ++score_index) {
          auto uncert_pair = get_tally_uncertainty(i_tally, score_index,
            filter_index);
          double std_dev = uncert_pair.first;
          double rel_err = uncert_pair.second;
          if (trigger.std_dev < std_dev) {
            trigger.std_dev = std_dev;
            trigger.variance = std_dev * std_dev;
          }
          if (trigger.rel_err < rel_err) {
            trigger.rel_err = rel_err;
          }

          double uncertainty;
          switch (trigger.type) {
            case VARIANCE:
              uncertainty = trigger.variance;
              break;
            case STANDARD_DEVIATION:
              uncertainty = trigger.std_dev;
              break;
            case RELATIVE_ERROR:
              uncertainty = trigger.rel_err;
          }

          double this_ratio = uncertainty / trigger.threshold;
          if (trigger.type == VARIANCE) {
            this_ratio = std::sqrt(ratio);
          }

          if (this_ratio > ratio) {
            ratio = this_ratio;
            int* scores;
            int junk;
            err = openmc_tally_get_scores(i_tally, &scores, &junk);
            //TODO: off-by-one
            score = scores[trigger.score_index-1];
            err = openmc_tally_get_id(i_tally, &tally_id);
          }
        }
      }
    }
  }
}

//! Computes the uncertainty/threshold ratio for the eigenvalue trigger.

double
check_keff_trigger()
{
  if (settings::run_mode != RUN_MODE_EIGENVALUE) return 0.;
  if (settings::keff_trigger.type == 0) return 0.;

  double k_combined[2];
  int err = openmc_get_keff(k_combined);

  double uncertainty = 0.;
  switch (settings::keff_trigger.type) {
  case VARIANCE:
    uncertainty = k_combined[1] * k_combined[1];
    break;
  case STANDARD_DEVIATION:
    uncertainty = k_combined[1];
    break;
  case RELATIVE_ERROR:
    uncertainty = k_combined[1] / k_combined[0];
  }

  double ratio = uncertainty / settings::keff_trigger.threshold;
  if (settings::keff_trigger.type == VARIANCE)
    ratio = std::sqrt(ratio);
  return ratio;
}

//! See if tally and eigenvalue uncertainties are under trigger thresholds.

void
check_triggers()
{
  // Make some aliases.
  const auto current_batch {simulation::current_batch};
  const auto n_batches {settings::n_batches};
  const auto interval {settings::trigger_batch_interval};

  // See if the current batch is one for which the triggers must be checked.
  if (!settings::trigger_on) return;
  if (current_batch < n_batches) return;
  if (((current_batch - n_batches) % interval) != 0) return;

  // Check the eigenvalue and tally triggers.
  double keff_ratio = check_keff_trigger();
  double tally_ratio;
  int tally_id, score;
  check_tally_triggers(tally_ratio, tally_id, score);

  // If all the triggers are satisfied, alert the user and return.
  if (std::max(keff_ratio, tally_ratio) <= 1.) {
    simulation::satisfy_triggers = true;
    write_message("Triggers satisfied for batch "
      + std::to_string(current_batch), 7);
    return;
  }

  // At least one trigger is unsatisfied.  Let the user know which one.
  simulation::satisfy_triggers = false;
  std::stringstream msg;
  msg << "Triggers unsatisfied, max unc./thresh. is ";
  if (keff_ratio >= tally_ratio) {
    msg << keff_ratio << " for eigenvalue";
  } else {
    msg << tally_ratio << " for " << reaction_name(score) << " in tally "
        << tally_id;
  }
  write_message(msg, 7);

  // Estimate batches til triggers are satisfied.
  if (settings::trigger_predict) {
    // This calculation assumes tally variance is proportional to 1/N where N is
    // the number of batches.
    auto max_ratio = std::max(keff_ratio, tally_ratio);
    auto n_active = current_batch - settings::n_inactive;
    auto n_pred_batches = static_cast<int>(n_active * max_ratio * max_ratio)
      + settings::n_inactive + 1;

    std::stringstream msg;
    msg << "The estimated number of batches is " << n_pred_batches;
    if (n_pred_batches > settings::n_max_batches) {
      msg << " --- greater than max batches";
      warning(msg);
    } else {
      write_message(msg, 7);
    }
  }
}

} // namespace openmc
