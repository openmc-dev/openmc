#include "openmc/tallies/derivative.h"

#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/tallies/tally.h"
#include "openmc/xml_interface.h"

#include <fmt/core.h>

template class std::vector<openmc::TallyDerivative>;

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  std::vector<TallyDerivative> tally_derivs;
  std::unordered_map<int, int> tally_deriv_map;
}

//==============================================================================
// TallyDerivative implementation
//==============================================================================

TallyDerivative::TallyDerivative(pugi::xml_node node)
{
  if (check_for_node(node, "id")) {
    id = std::stoi(get_node_value(node, "id"));
  } else {
    fatal_error("Must specify an ID for <derivative> elements in the tally "
                "XML file");
  }

  if (id <= 0)
    fatal_error("<derivative> IDs must be an integer greater than zero");

  std::string variable_str = get_node_value(node, "variable");

  if (variable_str == "density") {
    variable = DerivativeVariable::DENSITY;

  } else if (variable_str == "nuclide_density") {
    variable = DerivativeVariable::NUCLIDE_DENSITY;

    std::string nuclide_name = get_node_value(node, "nuclide");
    bool found = false;
    for (auto i = 0; i < data::nuclides.size(); ++i) {
      if (data::nuclides[i]->name_ == nuclide_name) {
        found = true;
        diff_nuclide = i;
      }
    }
    if (!found) {
      fatal_error(fmt::format("Could not find the nuclide \"{}\" specified in "
        "derivative {} in any material.", nuclide_name, id));
    }

  } else if (variable_str == "temperature") {
    variable = DerivativeVariable::TEMPERATURE;

  } else {
    fatal_error(fmt::format("Unrecognized variable \"{}\" on derivative {}",
      variable_str, id));
  }

  diff_material = std::stoi(get_node_value(node, "material"));
}

//==============================================================================
// Non-method functions
//==============================================================================

void
read_tally_derivatives(pugi::xml_node node)
{
  // Populate the derivatives array.
  for (auto deriv_node : node.children("derivative"))
    model::tally_derivs.emplace_back(deriv_node);

  // Fill the derivative map.
  for (auto i = 0; i < model::tally_derivs.size(); ++i) {
    auto id = model::tally_derivs[i].id;
    auto search = model::tally_deriv_map.find(id);
    if (search == model::tally_deriv_map.end()) {
      model::tally_deriv_map[id] = i;
    } else {
      fatal_error("Two or more derivatives use the same unique ID: "
                  + std::to_string(id));
    }
  }

  // Make sure derivatives were not requested for an MG run.
  if (!settings::run_CE && !model::tally_derivs.empty())
    fatal_error("Differential tallies not supported in multi-group mode");
}

void
apply_derivative_to_score(const Particle* p, int i_tally, int i_nuclide,
  double atom_density, int score_bin, double& score)
{
  const Tally& tally {*model::tallies[i_tally]};

  if (score == 0.0) return;

  // If our score was previously c then the new score is
  // c * (1/f * d_f/d_p + 1/c * d_c/d_p)
  // where (1/f * d_f/d_p) is the (logarithmic) flux derivative and p is the
  // perturbated variable.

  const auto& deriv {model::tally_derivs[tally.deriv_]};
  const auto flux_deriv = p->flux_derivs_[tally.deriv_];

  // Handle special cases where we know that d_c/d_p must be zero.
  if (score_bin == SCORE_FLUX) {
    score *= flux_deriv;
    return;
  } else if (p->material_ == MATERIAL_VOID) {
    score *= flux_deriv;
    return;
  }
  const Material& material {*model::materials[p->material_]};
  if (material.id_ != deriv.diff_material) {
    score *= flux_deriv;
    return;
  }

  switch (deriv.variable) {

  //============================================================================
  // Density derivative:
  // c = Sigma_MT
  // c = sigma_MT * N
  // c = sigma_MT * rho * const
  // d_c / d_rho = sigma_MT * const
  // (1 / c) * (d_c / d_rho) = 1 / rho

  case DerivativeVariable::DENSITY:
    switch (tally.estimator_) {

    case TallyEstimator::ANALOG:
    case TallyEstimator::COLLISION:
      switch (score_bin) {

      case SCORE_TOTAL:
      case SCORE_SCATTER:
      case SCORE_ABSORPTION:
      case SCORE_FISSION:
      case SCORE_NU_FISSION:
        score *= flux_deriv + 1. / material.density_gpcc_;
        break;

      default:
        fatal_error("Tally derivative not defined for a score on tally "
          + std::to_string(tally.id_));
      }
      break;

    default:
      fatal_error("Differential tallies are only implemented for analog and "
        "collision estimators.");
    }
    break;

  //============================================================================
  // Nuclide density derivative:
  // If we are scoring a reaction rate for a single nuclide then
  // c = Sigma_MT_i
  // c = sigma_MT_i * N_i
  // d_c / d_N_i = sigma_MT_i
  // (1 / c) * (d_c / d_N_i) = 1 / N_i
  // If the score is for the total material (i_nuclide = -1)
  // c = Sum_i(Sigma_MT_i)
  // d_c / d_N_i = sigma_MT_i
  // (1 / c) * (d_c / d_N) = sigma_MT_i / Sigma_MT
  // where i is the perturbed nuclide.

  case DerivativeVariable::NUCLIDE_DENSITY:
    switch (tally.estimator_) {

    case TallyEstimator::ANALOG:
      if (p->event_nuclide_ != deriv.diff_nuclide) {
        score *= flux_deriv;
        return;
      }

      switch (score_bin) {

      case SCORE_TOTAL:
      case SCORE_SCATTER:
      case SCORE_ABSORPTION:
      case SCORE_FISSION:
      case SCORE_NU_FISSION:
        {
          // Find the index of the perturbed nuclide.
          int i;
          for (i = 0; i < material.nuclide_.size(); ++i)
            if (material.nuclide_[i] == deriv.diff_nuclide) break;
          score *= flux_deriv + 1. / material.atom_density_(i);
        }
        break;

      default:
        fatal_error("Tally derivative not defined for a score on tally "
          + std::to_string(tally.id_));
      }
      break;

    case TallyEstimator::COLLISION:
      switch (score_bin) {

      case SCORE_TOTAL:
        if (i_nuclide == -1 && p->macro_xs_.total > 0.0) {
          score *= flux_deriv
            + p->neutron_xs_[deriv.diff_nuclide].total
            / p->macro_xs_.total;
        } else if (i_nuclide == deriv.diff_nuclide
                   && p->neutron_xs_[i_nuclide].total) {
          score *= flux_deriv + 1. / atom_density;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_SCATTER:
        if (i_nuclide == -1 && (p->macro_xs_.total
                                - p->macro_xs_.absorption) > 0.0) {
          score *= flux_deriv
            + (p->neutron_xs_[deriv.diff_nuclide].total
            - p->neutron_xs_[deriv.diff_nuclide].absorption)
            / (p->macro_xs_.total
            - p->macro_xs_.absorption);
        } else if (i_nuclide == deriv.diff_nuclide) {
          score *= flux_deriv + 1. / atom_density;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_ABSORPTION:
        if (i_nuclide == -1 && p->macro_xs_.absorption > 0.0) {
          score *= flux_deriv
            + p->neutron_xs_[deriv.diff_nuclide].absorption
            / p->macro_xs_.absorption;
        } else if (i_nuclide == deriv.diff_nuclide
                   && p->neutron_xs_[i_nuclide].absorption) {
          score *= flux_deriv + 1. / atom_density;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_FISSION:
        if (i_nuclide == -1 && p->macro_xs_.fission > 0.0) {
          score *= flux_deriv
            + p->neutron_xs_[deriv.diff_nuclide].fission
            / p->macro_xs_.fission;
        } else if (i_nuclide == deriv.diff_nuclide
                   && p->neutron_xs_[i_nuclide].fission) {
          score *= flux_deriv + 1. / atom_density;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_NU_FISSION:
        if (i_nuclide == -1 && p->macro_xs_.nu_fission > 0.0) {
          score *= flux_deriv
            + p->neutron_xs_[deriv.diff_nuclide].nu_fission
            / p->macro_xs_.nu_fission;
        } else if (i_nuclide == deriv.diff_nuclide
                   && p->neutron_xs_[i_nuclide].nu_fission) {
          score *= flux_deriv + 1. / atom_density;
        } else {
          score *= flux_deriv;
        }
        break;

      default:
        fatal_error("Tally derivative not defined for a score on tally "
          + std::to_string(tally.id_));
      }
      break;

    default:
      fatal_error("Differential tallies are only implemented for analog and "
        "collision estimators.");
    }
    break;

  //============================================================================
  // Temperature derivative:
  // If we are scoring a reaction rate for a single nuclide then
  // c = Sigma_MT_i
  // c = sigma_MT_i * N_i
  // d_c / d_T = (d_sigma_Mt_i / d_T) * N_i
  // (1 / c) * (d_c / d_T) = (d_sigma_MT_i / d_T) / sigma_MT_i
  // If the score is for the total material (i_nuclide = -1)
  // (1 / c) * (d_c / d_T) = Sum_i((d_sigma_MT_i / d_T) * N_i) / Sigma_MT_i
  // where i is the perturbed nuclide.  The d_sigma_MT_i / d_T term is
  // computed by multipole_deriv_eval.  It only works for the resolved
  // resonance range and requires multipole data.

  case DerivativeVariable::TEMPERATURE:
    switch (tally.estimator_) {

    case TallyEstimator::ANALOG:
      {
        // Find the index of the event nuclide.
        int i;
        for (i = 0; i < material.nuclide_.size(); ++i)
          if (material.nuclide_[i] == p->event_nuclide_) break;

        const auto& nuc {*data::nuclides[p->event_nuclide_]};
        if (!multipole_in_range(&nuc, p->E_last_)) {
          score *= flux_deriv;
          break;
        }

        switch (score_bin) {

        case SCORE_TOTAL:
          if (p->neutron_xs_[p->event_nuclide_].total) {
            double dsig_s, dsig_a, dsig_f;
            std::tie(dsig_s, dsig_a, dsig_f)
              = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
            score *= flux_deriv + (dsig_s + dsig_a) * material.atom_density_(i)
              / p->macro_xs_.total;
          } else {
            score *= flux_deriv;
          }
          break;

        case SCORE_SCATTER:
          if (p->neutron_xs_[p->event_nuclide_].total
              - p->neutron_xs_[p->event_nuclide_].absorption) {
            double dsig_s, dsig_a, dsig_f;
            std::tie(dsig_s, dsig_a, dsig_f)
              = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
            score *= flux_deriv + dsig_s * material.atom_density_(i)
              / (p->macro_xs_.total
              - p->macro_xs_.absorption);
          } else {
            score *= flux_deriv;
          }
          break;

        case SCORE_ABSORPTION:
          if (p->neutron_xs_[p->event_nuclide_].absorption) {
            double dsig_s, dsig_a, dsig_f;
            std::tie(dsig_s, dsig_a, dsig_f)
              = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
            score *= flux_deriv + dsig_a * material.atom_density_(i)
              / p->macro_xs_.absorption;
          } else {
            score *= flux_deriv;
          }
          break;

        case SCORE_FISSION:
          if (p->neutron_xs_[p->event_nuclide_].fission) {
            double dsig_s, dsig_a, dsig_f;
            std::tie(dsig_s, dsig_a, dsig_f)
              = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
            score *= flux_deriv + dsig_f * material.atom_density_(i)
              / p->macro_xs_.fission;
          } else {
            score *= flux_deriv;
          }
          break;

        case SCORE_NU_FISSION:
          if (p->neutron_xs_[p->event_nuclide_].fission) {
            double nu = p->neutron_xs_[p->event_nuclide_].nu_fission
              / p->neutron_xs_[p->event_nuclide_].fission;
            double dsig_s, dsig_a, dsig_f;
            std::tie(dsig_s, dsig_a, dsig_f)
              = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
            score *= flux_deriv + nu * dsig_f * material.atom_density_(i)
              / p->macro_xs_.nu_fission;
          } else {
            score *= flux_deriv;
          }
          break;

        default:
          fatal_error("Tally derivative not defined for a score on tally "
            + std::to_string(tally.id_));
        }
      }
      break;

    case TallyEstimator::COLLISION:
      if (i_nuclide != -1) {
        const auto& nuc {data::nuclides[i_nuclide]};
        if (!multipole_in_range(nuc.get(), p->E_last_)) {
          score *= flux_deriv;
          return;
        }
      }

      switch (score_bin) {

      case SCORE_TOTAL:
        if (i_nuclide == -1 && p->macro_xs_.total > 0.0) {
          double cum_dsig = 0;
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto i_nuc = material.nuclide_[i];
            const auto& nuc {*data::nuclides[i_nuc]};
            if (multipole_in_range(&nuc, p->E_last_)
                && p->neutron_xs_[i_nuc].total) {
              double dsig_s, dsig_a, dsig_f;
              std::tie(dsig_s, dsig_a, dsig_f)
                = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
              cum_dsig += (dsig_s + dsig_a) * material.atom_density_(i);
            }
          }
          score *= flux_deriv + cum_dsig / p->macro_xs_.total;
        } else if (p->neutron_xs_[i_nuclide].total) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          score *= flux_deriv
            + (dsig_s + dsig_a) / p->neutron_xs_[i_nuclide].total;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_SCATTER:
        if (i_nuclide == -1 && (p->macro_xs_.total
            - p->macro_xs_.absorption)) {
          double cum_dsig = 0;
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto i_nuc = material.nuclide_[i];
            const auto& nuc {*data::nuclides[i_nuc]};
            if (multipole_in_range(&nuc, p->E_last_)
                && (p->neutron_xs_[i_nuc].total
                - p->neutron_xs_[i_nuc].absorption)) {
              double dsig_s, dsig_a, dsig_f;
              std::tie(dsig_s, dsig_a, dsig_f)
                = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
              cum_dsig += dsig_s * material.atom_density_(i);
            }
          }
          score *= flux_deriv + cum_dsig / (p->macro_xs_.total
            - p->macro_xs_.absorption);
        } else if (p->neutron_xs_[i_nuclide].total
                   - p->neutron_xs_[i_nuclide].absorption) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          score *= flux_deriv + dsig_s / (p->neutron_xs_[i_nuclide].total
            - p->neutron_xs_[i_nuclide].absorption);
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_ABSORPTION:
        if (i_nuclide == -1 && p->macro_xs_.absorption > 0.0) {
          double cum_dsig = 0;
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto i_nuc = material.nuclide_[i];
            const auto& nuc {*data::nuclides[i_nuc]};
            if (multipole_in_range(&nuc, p->E_last_)
                && p->neutron_xs_[i_nuc].absorption) {
              double dsig_s, dsig_a, dsig_f;
              std::tie(dsig_s, dsig_a, dsig_f)
                = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
              cum_dsig += dsig_a * material.atom_density_(i);
            }
          }
          score *= flux_deriv + cum_dsig / p->macro_xs_.absorption;
        } else if (p->neutron_xs_[i_nuclide].absorption) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          score *= flux_deriv
            + dsig_a / p->neutron_xs_[i_nuclide].absorption;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_FISSION:
        if (i_nuclide == -1 && p->macro_xs_.fission > 0.0) {
          double cum_dsig = 0;
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto i_nuc = material.nuclide_[i];
            const auto& nuc {*data::nuclides[i_nuc]};
            if (multipole_in_range(&nuc, p->E_last_)
                && p->neutron_xs_[i_nuc].fission) {
              double dsig_s, dsig_a, dsig_f;
              std::tie(dsig_s, dsig_a, dsig_f)
                = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
              cum_dsig += dsig_f * material.atom_density_(i);
            }
          }
          score *= flux_deriv + cum_dsig / p->macro_xs_.fission;
        } else if (p->neutron_xs_[i_nuclide].fission) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          score *= flux_deriv
            + dsig_f / p->neutron_xs_[i_nuclide].fission;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_NU_FISSION:
        if (i_nuclide == -1 && p->macro_xs_.nu_fission > 0.0) {
          double cum_dsig = 0;
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto i_nuc = material.nuclide_[i];
            const auto& nuc {*data::nuclides[i_nuc]};
            if (multipole_in_range(&nuc, p->E_last_)
                && p->neutron_xs_[i_nuc].fission) {
              double nu = p->neutron_xs_[i_nuc].nu_fission
                / p->neutron_xs_[i_nuc].fission;
              double dsig_s, dsig_a, dsig_f;
              std::tie(dsig_s, dsig_a, dsig_f)
                = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
              cum_dsig += nu * dsig_f * material.atom_density_(i);
            }
          }
          score *= flux_deriv + cum_dsig / p->macro_xs_.nu_fission;
        } else if (p->neutron_xs_[i_nuclide].fission) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          score *= flux_deriv
            + dsig_f / p->neutron_xs_[i_nuclide].fission;
        } else {
          score *= flux_deriv;
        }
        break;

      default:
        break;
      }
      break;

    default:
      fatal_error("Differential tallies are only implemented for analog and "
        "collision estimators.");
    }
    break;
  }
}

void
score_track_derivative(Particle* p, double distance)
{
  // A void material cannot be perturbed so it will not affect flux derivatives.
  if (p->material_ == MATERIAL_VOID) return;
  const Material& material {*model::materials[p->material_]};

  for (auto idx = 0; idx < model::tally_derivs.size(); idx++) {
    const auto& deriv = model::tally_derivs[idx];
    auto& flux_deriv = p->flux_derivs_[idx];
    if (deriv.diff_material != material.id_) continue;

    switch (deriv.variable) {

    case DerivativeVariable::DENSITY:
      // phi is proportional to e^(-Sigma_tot * dist)
      // (1 / phi) * (d_phi / d_rho) = - (d_Sigma_tot / d_rho) * dist
      // (1 / phi) * (d_phi / d_rho) = - Sigma_tot / rho * dist
      flux_deriv -= distance * p->macro_xs_.total
        / material.density_gpcc_;
      break;

    case DerivativeVariable::NUCLIDE_DENSITY:
      // phi is proportional to e^(-Sigma_tot * dist)
      // (1 / phi) * (d_phi / d_N) = - (d_Sigma_tot / d_N) * dist
      // (1 / phi) * (d_phi / d_N) = - sigma_tot * dist
      flux_deriv -= distance
        * p->neutron_xs_[deriv.diff_nuclide].total;
      break;

    case DerivativeVariable::TEMPERATURE:
      for (auto i = 0; i < material.nuclide_.size(); ++i) {
        const auto& nuc {*data::nuclides[material.nuclide_[i]]};
        if (multipole_in_range(&nuc, p->E_last_)) {
          // phi is proportional to e^(-Sigma_tot * dist)
          // (1 / phi) * (d_phi / d_T) = - (d_Sigma_tot / d_T) * dist
          // (1 / phi) * (d_phi / d_T) = - N (d_sigma_tot / d_T) * dist
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_, p->sqrtkT_);
          flux_deriv -= distance * (dsig_s + dsig_a)
            * material.atom_density_(i);
        }
      }
      break;
    }
  }
}

void score_collision_derivative(Particle* p)
{
  // A void material cannot be perturbed so it will not affect flux derivatives.
  if (p->material_ == MATERIAL_VOID) return;

  const Material& material {*model::materials[p->material_]};

  for (auto idx = 0; idx < model::tally_derivs.size(); idx++) {
    const auto& deriv = model::tally_derivs[idx];
    auto& flux_deriv = p->flux_derivs_[idx];

    if (deriv.diff_material != material.id_) continue;

    switch (deriv.variable) {

    case DerivativeVariable::DENSITY:
      // phi is proportional to Sigma_s
      // (1 / phi) * (d_phi / d_rho) = (d_Sigma_s / d_rho) / Sigma_s
      // (1 / phi) * (d_phi / d_rho) = 1 / rho
      flux_deriv += 1. / material.density_gpcc_;
      break;

    case DerivativeVariable::NUCLIDE_DENSITY:
      if (p->event_nuclide_ != deriv.diff_nuclide) continue;
      // Find the index in this material for the diff_nuclide.
      int i;
      for (i = 0; i < material.nuclide_.size(); ++i)
        if (material.nuclide_[i] == deriv.diff_nuclide) break;
      // Make sure we found the nuclide.
      if (material.nuclide_[i] != deriv.diff_nuclide) {
        fatal_error(fmt::format(
          "Could not find nuclide {} in material {} for tally derivative {}",
          data::nuclides[deriv.diff_nuclide]->name_, material.id_, deriv.id));
      }
      // phi is proportional to Sigma_s
      // (1 / phi) * (d_phi / d_N) = (d_Sigma_s / d_N) / Sigma_s
      // (1 / phi) * (d_phi / d_N) = sigma_s / Sigma_s
      // (1 / phi) * (d_phi / d_N) = 1 / N
      flux_deriv += 1. / material.atom_density_(i);
      break;

    case DerivativeVariable::TEMPERATURE:
      // Loop over the material's nuclides until we find the event nuclide.
      for (auto i_nuc : material.nuclide_) {
        const auto& nuc {*data::nuclides[i_nuc]};
        if (i_nuc == p->event_nuclide_ && multipole_in_range(&nuc, p->E_last_)) {
          // phi is proportional to Sigma_s
          // (1 / phi) * (d_phi / d_T) = (d_Sigma_s / d_T) / Sigma_s
          // (1 / phi) * (d_phi / d_T) = (d_sigma_s / d_T) / sigma_s
          const auto& micro_xs {p->neutron_xs_[i_nuc]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          flux_deriv += dsig_s / (micro_xs.total - micro_xs.absorption);
          // Note that this is an approximation!  The real scattering cross
          // section is
          // Sigma_s(E'->E, u'->u) = Sigma_s(E') * P(E'->E, u'->u).
          // We are assuming that d_P(E'->E, u'->u) / d_T = 0 and only
          // computing d_S(E') / d_T.  Using this approximation in the vicinity
          // of low-energy resonances causes errors (~2-5% for PWR pincell
          // eigenvalue derivatives).
        }
      }
      break;
    }
  }
}

}// namespace openmc
