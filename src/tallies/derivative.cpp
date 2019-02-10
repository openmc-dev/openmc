#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/tallies/derivative.h"
#include "openmc/xml_interface.h"

#include <sstream>

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
    variable = DIFF_DENSITY;

  } else if (variable_str == "nuclide_density") {
    variable = DIFF_NUCLIDE_DENSITY;

    std::string nuclide_name = get_node_value(node, "nuclide");
    bool found = false;
    for (auto i = 0; i < data::nuclides.size(); ++i) {
      if (data::nuclides[i]->name_ == nuclide_name) {
        found = true;
        //TODO: off-by-one
        diff_nuclide = i + 1;
      }
    }
    if (!found) {
      std::stringstream out;
      out << "Could not find the nuclide \"" << nuclide_name
          << "\" specified in derivative " << id << " in any material.";
      fatal_error(out);
    }

  } else if (variable_str == "temperature") {
    variable = DIFF_TEMPERATURE;

  } else {
    std::stringstream out;
    out << "Unrecognized variable \"" << variable_str
        << "\" on derivative " << id;
    fatal_error(out);
  }

  diff_material = std::stoi(get_node_value(node, "material"));
}

//==============================================================================
// Non-method functions
//==============================================================================

extern "C" void
read_tally_derivatives(pugi::xml_node* node)
{
  // Populate the derivatives array.  This must be done in parallel because
  // the derivatives are threadprivate.
  #pragma omp parallel
  {
    for (auto deriv_node : node->children("derivative"))
      model::tally_derivs.push_back(deriv_node);
  }

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

//! Adjust diff tally flux derivatives for a particle tracking event.

extern "C" void
score_track_derivative(Particle* p, double distance)
{
  // A void material cannot be perturbed so it will not affect flux derivatives.
  if (p->material == MATERIAL_VOID) return;
  //TODO: off-by-one
  const Material& material {*model::materials[p->material-1]};

  for (auto& deriv : model::tally_derivs) {
    if (deriv.diff_material != material.id_) continue;

    switch (deriv.variable) {

    case DIFF_DENSITY:
      // phi is proportional to e^(-Sigma_tot * dist)
      // (1 / phi) * (d_phi / d_rho) = - (d_Sigma_tot / d_rho) * dist
      // (1 / phi) * (d_phi / d_rho) = - Sigma_tot / rho * dist
      deriv.flux_deriv -= distance * simulation::material_xs.total
        / material.density_gpcc_;
      break;

    case DIFF_NUCLIDE_DENSITY:
      // phi is proportional to e^(-Sigma_tot * dist)
      // (1 / phi) * (d_phi / d_N) = - (d_Sigma_tot / d_N) * dist
      // (1 / phi) * (d_phi / d_N) = - sigma_tot * dist
      //TODO: off-by-one
      deriv.flux_deriv -= distance
        * simulation::micro_xs[deriv.diff_nuclide-1].total;
      break;

    case DIFF_TEMPERATURE:
      for (auto i = 0; i < material.nuclide_.size(); ++i) {
        const auto& nuc {*data::nuclides[material.nuclide_[i]]};
        //TODO: off-by-one
        if (multipole_in_range(&nuc, p->last_E)) {
          // phi is proportional to e^(-Sigma_tot * dist)
          // (1 / phi) * (d_phi / d_T) = - (d_Sigma_tot / d_T) * dist
          // (1 / phi) * (d_phi / d_T) = - N (d_sigma_tot / d_T) * dist
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E, p->sqrtkT);
          deriv.flux_deriv -= distance * (dsig_s + dsig_a)
            * material.atom_density_(i);
        }
      }
      break;
    }
  }
}

//! Adjust diff tally flux derivatives for a particle scattering event.
//
//! Note that this subroutine will be called after absorption events in
//! addition to scattering events, but any flux derivatives scored after an
//! absorption will never be tallied.  The paricle will be killed before any
//! further tallies are scored.

extern "C" void
score_collision_derivative(Particle* p)
{
  // A void material cannot be perturbed so it will not affect flux derivatives.
  if (p->material == MATERIAL_VOID) return;
  //TODO: off-by-one
  const Material& material {*model::materials[p->material-1]};

  for (auto& deriv : model::tally_derivs) {
    if (deriv.diff_material != material.id_) continue;

    switch (deriv.variable) {

    case DIFF_DENSITY:
      // phi is proportional to Sigma_s
      // (1 / phi) * (d_phi / d_rho) = (d_Sigma_s / d_rho) / Sigma_s
      // (1 / phi) * (d_phi / d_rho) = 1 / rho
      deriv.flux_deriv += 1. / material.density_gpcc_;
      break;

    case DIFF_NUCLIDE_DENSITY:
      //TODO: off-by-one throughout on diff_nuclide
      if (p->event_nuclide != deriv.diff_nuclide) continue;
      // Find the index in this material for the diff_nuclide.
      int i;
      for (i = 0; i < material.nuclide_.size(); ++i)
        if (material.nuclide_[i] == deriv.diff_nuclide - 1) break;
      // Make sure we found the nuclide.
      if (material.nuclide_[i] != deriv.diff_nuclide - 1) {
        std::stringstream err_msg;
        err_msg << "Could not find nuclide "
          << data::nuclides[deriv.diff_nuclide-1]->name_ << " in material "
          << material.id_ << " for tally derivative " << deriv.id;
        fatal_error(err_msg);
      }
      // phi is proportional to Sigma_s
      // (1 / phi) * (d_phi / d_N) = (d_Sigma_s / d_N) / Sigma_s
      // (1 / phi) * (d_phi / d_N) = sigma_s / Sigma_s
      // (1 / phi) * (d_phi / d_N) = 1 / N
      deriv.flux_deriv += 1. / material.atom_density_(i);
      break;

    case DIFF_TEMPERATURE:
      // Loop over the material's nuclides until we find the event nuclide.
      for (auto i_nuc : material.nuclide_) {
        const auto& nuc {*data::nuclides[i_nuc]};
        //TODO: off-by-one
        if (i_nuc == p->event_nuclide - 1
            && multipole_in_range(&nuc, p->last_E)) {
          // phi is proportional to Sigma_s
          // (1 / phi) * (d_phi / d_T) = (d_Sigma_s / d_T) / Sigma_s
          // (1 / phi) * (d_phi / d_T) = (d_sigma_s / d_T) / sigma_s
          const auto& micro_xs {simulation::micro_xs[i_nuc]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->last_E, p->sqrtkT);
          deriv.flux_deriv += dsig_s / (micro_xs.total - micro_xs.absorption);
          // Note that this is an approximation!  The real scattering cross
          // section is
          // Sigma_s(E'->E, uvw'->uvw) = Sigma_s(E') * P(E'->E, uvw'->uvw). 
          // We are assuming that d_P(E'->E, uvw'->uvw) / d_T = 0 and only
          // computing d_S(E') / d_T.  Using this approximation in the vicinity
          // of low-energy resonances causes errors (~2-5% for PWR pincell
          // eigenvalue derivatives).
        }
      }
      break;
    }
  }
}

//! Set the flux derivatives on differential tallies to zero.
extern "C" void
zero_flux_derivs()
{
  for (auto& deriv : model::tally_derivs) deriv.flux_deriv = 0.;
}

//==============================================================================
// Fortran interop
//==============================================================================

extern "C" int n_tally_derivs() {return model::tally_derivs.size();}

extern "C" TallyDerivative*
tally_deriv_c(int i)
{
  return &(model::tally_derivs[i]);
}

}// namespace openmc
