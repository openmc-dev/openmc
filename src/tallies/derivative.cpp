#include "openmc/error.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/tallies/derivative.h"
#include "openmc/xml_interface.h"

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
