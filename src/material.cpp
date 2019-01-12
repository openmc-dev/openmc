#include "openmc/material.h"

#include <cmath>
#include <string>
#include <sstream>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xoperation.hpp"
#include "xtensor/xview.hpp"

#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/search.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {

std::vector<Material*> materials;
std::unordered_map<int32_t, int32_t> material_map;

} // namespace model

//==============================================================================
// Material implementation
//==============================================================================

Material::Material(pugi::xml_node node)
{
  if (check_for_node(node, "id")) {
    id_ = std::stoi(get_node_value(node, "id"));
  } else {
    fatal_error("Must specify id of material in materials XML file.");
  }

  if (check_for_node(node, "temperature")) {
    temperature_ = std::stod(get_node_value(node, "temperature"));
  }

  if (check_for_node(node, "volume")) {
    volume_ = std::stod(get_node_value(node, "volume"));
  }
}

void Material::init_bremsstrahlung()
{
  // Create new object
  ttb_ = std::make_unique<Bremsstrahlung>();

  // Get the size of the energy grids
  auto n_k = data::ttb_k_grid.size();
  auto n_e = data::ttb_e_grid.size();

  // Get pointers to nuclides, elements, densities
  int32_t index;
  openmc_get_material_index(id_, &index);
  int* nuclide_;
  double* atom_density_;
  int n;
  openmc_material_get_densities(index, &nuclide_, &atom_density_, &n);
  int* element_ = material_element(index);

  for (int particle = 0; particle < 2; ++particle) {
    // Loop over logic twice, once for electron, once for positron
    BremsstrahlungData* ttb = (particle == 0) ? &ttb_->electron : &ttb_->positron;
    bool positron = (particle == 1);

    // Allocate arrays for TTB data
    ttb->pdf = xt::zeros<double>({n_e, n_e});
    ttb->cdf = xt::zeros<double>({n_e, n_e});
    ttb->yield = xt::empty<double>({n_e});

    // Allocate temporary arrays
    xt::xtensor<double, 1> stopping_power_collision({n_e}, 0.0);
    xt::xtensor<double, 1> stopping_power_radiative({n_e}, 0.0);
    xt::xtensor<double, 2> dcs({n_e, n_k}, 0.0);

    double Z_eq_sq = 0.0;
    double sum_density = 0.0;

    // Calculate the molecular DCS and the molecular total stopping power using
    // Bragg's additivity rule.
    // TODO: The collision stopping power cannot be accurately calculated using
    // Bragg's additivity rule since the mean excitation energies and the
    // density effect corrections cannot simply be summed together. Bragg's
    // additivity rule fails especially when a higher-density compound is
    // composed of elements that are in lower-density form at normal temperature
    // and pressure (at which the NIST stopping powers are given). It will be
    // used to approximate the collision stopping powers for now, but should be
    // fixed in the future.
    for (int i = 0; i < n; ++i) {
      // Get pointer to current element
      // TODO: off-by-one
      const auto& elm = data::elements[element_[i] - 1];
      // TODO: off-by-one
      double awr = data::nuclides[nuclide_[i] - 1]->awr_;

      // Get atomic density and mass density of nuclide given atom/weight percent
      double atom_density = (atom_density_[0] > 0.0) ?
        atom_density_[i] : -atom_density_[i] / awr;
      double mass_density = atom_density * awr;

      // Calculate the "equivalent" atomic number Zeq of the material
      Z_eq_sq += atom_density * elm.Z_ * elm.Z_;
      sum_density += atom_density;

      // Accumulate material DCS
      dcs += (atom_density * elm.Z_ * elm.Z_) * elm.dcs_;

      // Accumulate material collision stopping power
      stopping_power_collision += (mass_density * MASS_NEUTRON / N_AVOGADRO)
        * elm.stopping_power_collision_;

      // Accumulate material radiative stopping power
      stopping_power_radiative += (mass_density * MASS_NEUTRON / N_AVOGADRO)
        * elm.stopping_power_radiative_;
    }
    Z_eq_sq /= sum_density;

    // Calculate the positron DCS and radiative stopping power. These are
    // obtained by multiplying the electron DCS and radiative stopping powers by
    // a factor r, which is a numerical approximation of the ratio of the
    // radiative stopping powers for positrons and electrons. Source: F. Salvat,
    // J. M. Fernández-Varea, and J. Sempau, "PENELOPE-2011: A Code System for
    // Monte Carlo Simulation of Electron and Photon Transport," OECD-NEA,
    // Issy-les-Moulineaux, France (2011).
    if (positron) {
      for (int i = 0; i < n_e; ++i) {
        double t = std::log(1.0 + 1.0e6*data::ttb_e_grid(i)/(Z_eq_sq*MASS_ELECTRON_EV));
        double r = 1.0 - std::exp(-1.2359e-1*t + 6.1274e-2*std::pow(t, 2)
          - 3.1516e-2*std::pow(t, 3) + 7.7446e-3*std::pow(t, 4)
          - 1.0595e-3*std::pow(t, 5) + 7.0568e-5*std::pow(t, 6)
          - 1.808e-6*std::pow(t, 7));
        stopping_power_radiative(i) *= r;
        auto dcs_i = xt::view(dcs, i, xt::all());
        dcs_i *= r;
      }
    }

    // Total material stopping power
    xt::xtensor<double, 1> stopping_power = stopping_power_collision +
      stopping_power_radiative;

    // Loop over photon energies
    xt::xtensor<double, 1> f({n_e}, 0.0);
    xt::xtensor<double, 1> z({n_e}, 0.0);
    for (int i = 0; i < n_e - 1; ++i) {
      double w = data::ttb_e_grid(i);

      // Loop over incident particle energies
      for (int j = i; j < n_e; ++j) {
        double e = data::ttb_e_grid(j);

        // Reduced photon energy
        double k = w / e;

        // Find the lower bounding index of the reduced photon energy
        int i_k = lower_bound_index(data::ttb_k_grid.cbegin(),
          data::ttb_k_grid.cend(), k);

        // Get the interpolation bounds
        double k_l = data::ttb_k_grid(i_k);
        double k_r = data::ttb_k_grid(i_k + 1);
        double x_l = dcs(j, i_k);
        double x_r = dcs(j, i_k + 1);

        // Find the value of the DCS using linear interpolation in reduced
        // photon energy k
        double x = x_l + (k - k_l)*(x_r - x_l)/(k_r - k_l);

        // Ratio of the velocity of the charged particle to the speed of light
        double beta = std::sqrt(e*(e + 2.0*MASS_ELECTRON_EV)) /
          (e + MASS_ELECTRON_EV);

        // Compute the integrand of the PDF
        f(j) = x / (beta*beta * stopping_power(j) * w);
      }

      // Number of points to integrate
      int n = n_e - i;

      // Integrate the PDF using cubic spline integration over the incident
      // particle energy
      if (n > 2) {
        spline_c(n, &data::ttb_e_grid(i), &f(i), &z(i));

        double c = 0.0;
        for (int j = i; j < n_e - 1; ++j) {
          c += spline_integrate_c(n, &data::ttb_e_grid(i), &f(i), &z(i),
            data::ttb_e_grid(j), data::ttb_e_grid(j+1));

          ttb->pdf(j+1,i) = c;
        }

      // Integrate the last two points using trapezoidal rule in log-log space
      } else {
        double e_l = std::log(data::ttb_e_grid(i));
        double e_r = std::log(data::ttb_e_grid(i+1));
        double x_l = std::log(f(i));
        double x_r = std::log(f(i+1));

        ttb->pdf(i+1,i) = 0.5*(e_r - e_l)*(std::exp(e_l + x_l) + std::exp(e_r + x_r));
      }
    }

    // Loop over incident particle energies
    for (int j = 1; j < n_e; ++j) {
      // Set last element of PDF to small non-zero value to enable log-log
      // interpolation
      ttb->pdf(j,j) = std::exp(-500.0);

      // Loop over photon energies
      double c = 0.0;
      for (int i = 0; i < j; ++i) {
        // Integrate the CDF from the PDF using the trapezoidal rule in log-log
        // space
        double w_l = std::log(data::ttb_e_grid(i));
        double w_r = std::log(data::ttb_e_grid(i+1));
        double x_l = std::log(ttb->pdf(j,i));
        double x_r = std::log(ttb->pdf(j,i+1));

        c += 0.5*(w_r - w_l)*(std::exp(w_l + x_l) + std::exp(w_r + x_r));
        ttb->cdf(j,i+1) = c;
      }

      // Set photon number yield
      ttb->yield(j) = c;
    }

    // Use logarithm of number yield since it is log-log interpolated
    ttb->yield = xt::where(ttb->yield > 0.0, xt::log(ttb->yield), -500.0);
  }
}

//==============================================================================
// Non-method functions
//==============================================================================

extern "C" void
read_materials(pugi::xml_node* node)
{
  // Loop over XML material elements and populate the array.
  for (pugi::xml_node material_node : node->children("material")) {
    model::materials.push_back(new Material(material_node));
  }
  model::materials.shrink_to_fit();

  // Populate the material map.
  for (int i = 0; i < model::materials.size(); i++) {
    int32_t mid = model::materials[i]->id_;
    auto search = model::material_map.find(mid);
    if (search == model::material_map.end()) {
      model::material_map[mid] = i;
    } else {
      std::stringstream err_msg;
      err_msg << "Two or more materials use the same unique ID: " << mid;
      fatal_error(err_msg);
    }
  }
}

//==============================================================================
// C API
//==============================================================================

extern "C" int
openmc_material_get_volume(int32_t index, double* volume)
{
  if (index >= 1 && index <= model::materials.size()) {
    Material* m = model::materials[index - 1];
    if (m->volume_ >= 0.0) {
      *volume = m->volume_;
      return 0;
    } else {
      std::stringstream msg;
      msg << "Volume for material with ID=" << m->id_ << " not set.";
      set_errmsg(msg);
      return OPENMC_E_UNASSIGNED;
    }
  } else {
    set_errmsg("Index in materials array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

extern "C" int
openmc_material_set_volume(int32_t index, double volume)
{
  if (index >= 1 && index <= model::materials.size()) {
    Material* m = model::materials[index - 1];
    if (volume >= 0.0) {
      m->volume_ = volume;
      return 0;
    } else {
      set_errmsg("Volume must be non-negative");
      return OPENMC_E_INVALID_ARGUMENT;
    }
  } else {
    set_errmsg("Index in materials array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

//==============================================================================
// Fortran compatibility functions
//==============================================================================

extern "C" {
  Material* material_pointer(int32_t indx) {return model::materials[indx];}

  int32_t material_id(Material* mat) {return mat->id_;}

  void material_set_id(Material* mat, int32_t id, int32_t index)
  {
    mat->id_ = id;
    //TODO: off-by-one
    model::material_map[id] = index - 1;
  }

  bool material_fissionable(Material* mat) {return mat->fissionable;}

  void material_set_fissionable(Material* mat, bool fissionable)
  {
    mat->fissionable = fissionable;
  }

  void material_init_bremsstrahlung(Material* mat)
  {
    mat->init_bremsstrahlung();
  }

  void extend_materials_c(int32_t n)
  {
    model::materials.reserve(model::materials.size() + n);
    for (int32_t i = 0; i < n; i++) {
      model::materials.push_back(new Material());
    }
  }

  void free_memory_material_c()
  {
    for (Material *mat : model::materials) {delete mat;}
    model::materials.clear();
    model::material_map.clear();
  }
}

} // namespace openmc
