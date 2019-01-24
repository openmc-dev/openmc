#include "openmc/material.h"

#include <algorithm> // for min, max, sort
#include <cmath>
#include <iterator>
#include <string>
#include <sstream>

#include "xtensor/xbuilder.hpp"
#include "xtensor/xoperation.hpp"
#include "xtensor/xview.hpp"

#include "openmc/capi.h"
#include "openmc/cross_sections.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/math_functions.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/string_utils.h"
#include "openmc/thermal.h"
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

  if (check_for_node(node, "name")) {
    name_ = get_node_value(node, "name");
  }

  if (check_for_node(node, "depletable")) {
    depletable_ = get_node_value_bool(node, "depletable");
  }

  bool sum_density {false};
  pugi::xml_node density_node = node.child("density");
  std::string units;
  if (density_node) {
    units = get_node_value(density_node, "units");
    if (units == "sum") {
      sum_density = true;
    } else if (units == "macro") {
      if (check_for_node(density_node, "value")) {
        density_ = std::stod(get_node_value(density_node, "value"));
      } else {
        density_ = 1.0;
      }
    } else {
      double val = std::stod(get_node_value(density_node, "value"));
      if (val <= 0.0) {
        fatal_error("Need to specify a positive density on material "
          + std::to_string(id_) + ".");
      }

      if (units == "g/cc" || units == "g/cm3") {
        density_ = -val;
      } else if (units == "kg/m3") {
        density_ = -1.0e-3 * val;
      } else if (units == "atom/b-cm") {
        density_ = val;
      } else if (units == "atom/cc" || units == "atom/cm3") {
        density_ = 1.0e-24 * val;
      } else {
        fatal_error("Unknown units '" + units + "' specified on material "
          + std::to_string(id_) + ".");
      }
    }
  } else {
    fatal_error("Must specify <density> element in material "
      + std::to_string(id_) + ".");
  }

  if (node.child("element")) {
    fatal_error("Unable to add an element to material " + std::to_string(id_) +
      " since the element option has been removed from the xml input. "
      "Elements can only be added via the Python API, which will expand "
      "elements into their natural nuclides.");
  }

  // =======================================================================
  // READ AND PARSE <nuclide> TAGS

  // Check to ensure material has at least one nuclide
  if (!check_for_node(node, "nuclide") && !check_for_node(node, "macroscopic")) {
    fatal_error("No macroscopic data or nuclides specified on material "
      + std::to_string(id_));
  }

  // Create list of macroscopic x/s based on those specified, just treat
  // them as nuclides. This is all really a facade so the user thinks they
  // are entering in macroscopic data but the code treats them the same
  // as nuclides internally.
  // Get pointer list of XML <macroscopic>
  auto node_macros = node.children("macroscopic");
  int num_macros = std::distance(node_macros.begin(), node_macros.end());

  std::vector<std::string> names;
  std::vector<double> densities;
  if (settings::run_CE && num_macros > 0) {
    fatal_error("Macroscopic can not be used in continuous-energy mode.");
  } else if (num_macros > 1) {
    fatal_error("Only one macroscopic object permitted per material, "
      + std::to_string(id_));
  } else if (num_macros == 1) {
    pugi::xml_node node_nuc = *node_macros.begin();

    // Check for empty name on nuclide
    if (!check_for_node(node_nuc, "name")) {
      fatal_error("No name specified on macroscopic data in material "
        + std::to_string(id_));
    }

    // store nuclide name
    std::string name = get_node_value(node_nuc, "name", false, true);
    names.push_back(name);

    // Set density for macroscopic data
    if (units == "macro") {
      densities.push_back(1.0);
    } else {
      fatal_error("Units can only be macro for macroscopic data " + name);
    }
  } else {
    // Create list of nuclides based on those specified
    for (auto node_nuc : node.children("nuclide")) {
      // Check for empty name on nuclide
      if (!check_for_node(node_nuc, "name")) {
        fatal_error("No name specified on nuclide in material "
          + std::to_string(id_));
      }

      // store nuclide name
      std::string name = get_node_value(node_nuc, "name", false, true);
      names.push_back(name);

      // Check if no atom/weight percents were specified or if both atom and
      // weight percents were specified
      if (units == "macro") {
        densities.push_back(1.0);
      } else {
        bool has_ao = check_for_node(node_nuc, "ao");
        bool has_wo = check_for_node(node_nuc, "wo");

        if (!has_ao && !has_wo) {
          fatal_error("No atom or weight percent specified for nuclide: " + name);
        } else if (has_ao && has_wo) {
          fatal_error("Cannot specify both atom and weight percents for a "
            "nuclide: " + name);
        }

        // Copy atom/weight percents
        if (has_ao) {
          densities.push_back(std::stod(get_node_value(node_nuc, "ao")));
        } else {
          densities.push_back(-std::stod(get_node_value(node_nuc, "wo")));
        }
      }
    }
  }

  // =======================================================================
  // READ AND PARSE <isotropic> element

  std::vector<std::string> iso_lab;
  if (check_for_node(node, "isotropic")) {
    iso_lab = get_node_array<std::string>(node, "isotropic");
  }

  // ========================================================================
  // COPY NUCLIDES TO ARRAYS IN MATERIAL

  // allocate arrays in Material object
  auto n = names.size();
  nuclide_.reserve(n);
  atom_density_ = xt::empty<double>({n});
  if (settings::photon_transport) element_.reserve(n);

  for (int i = 0; i < names.size(); ++i) {
    const auto& name {names[i]};

    // Check that this nuclide is listed in the cross_sections.xml file
    LibraryKey key {Library::Type::neutron, name};
    if (data::library_map.find(key) == data::library_map.end()) {
      fatal_error("Could not find nuclide " + name + " in cross_sections.xml.");
    }

    // If this nuclide hasn't been encountered yet, we need to add its name
    // and alias to the nuclide_dict
    if (data::nuclide_map.find(name) == data::nuclide_map.end()) {
      int index = data::nuclide_map.size();
      data::nuclide_map[name] = index;
      nuclide_.push_back(index);
    } else {
      nuclide_.push_back(data::nuclide_map[name]);
    }

    // If the corresponding element hasn't been encountered yet and photon
    // transport will be used, we need to add its symbol to the element_dict
    if (settings::photon_transport) {
      int pos = name.find_first_of("0123456789");
      std::string element = name.substr(pos);

      // Make sure photon cross section data is available
      LibraryKey key {Library::Type::photon, element};
      if (data::library_map.find(key) == data::library_map.end()) {
        fatal_error("Could not find element " + element
          + " in cross_sections.xml.");
      }

      if (data::element_map.find(element) == data::element_map.end()) {
        int index = data::element_map.size();
        data::element_map[element] = index;
        element_.push_back(index);
      } else {
        element_.push_back(data::element_map[element]);
      }
    }

    // Copy name and atom/weight percent
    //mat % names(j) = name
    atom_density_(i) = densities[i];
  }

  if (settings::run_CE) {
    // By default, isotropic-in-lab is not used
    if (iso_lab.size() > 0) {
      has_isotropic_nuclides_ = true;
      p0_.resize(n);

      // Apply isotropic-in-lab treatment to specified nuclides
      for (const auto& nuc : iso_lab) {
        for (int j = 0; j < names.size(); ++j) {
          if (names[j] == nuc) p0_[j] = true;
        }
      }
    }
  }

  // Check to make sure either all atom percents or all weight percents are
  // given
  if (!(xt::all(atom_density_ >= 0.0) || xt::all(atom_density_ <= 0.0))) {
    fatal_error("Cannot mix atom and weight percents in material "
      + std::to_string(id_));
  }

  // Determine density if it is a sum value
  if (sum_density) density_ = xt::sum(atom_density_)();


  if (check_for_node(node, "temperature")) {
    temperature_ = std::stod(get_node_value(node, "temperature"));
  }

  if (check_for_node(node, "volume")) {
    volume_ = std::stod(get_node_value(node, "volume"));
  }

  // =======================================================================
  // READ AND PARSE <sab> TAG FOR THERMAL SCATTERING DATA
  if (settings::run_CE) {
    // Loop over <sab> elements

    std::vector<std::string> sab_names;
    for (auto node_sab : node.children("sab")) {
      // Determine name of thermal scattering table
      if (!check_for_node(node_sab, "name")) {
        fatal_error("Need to specify <name> for thermal scattering table.");
      }
      std::string name = get_node_value(node_sab, "name");
      sab_names.push_back(name);

      // Read the fraction of nuclei affected by this thermal scattering table
      double fraction = 1.0;
      if (check_for_node(node_sab, "fraction")) {
        fraction = std::stod(get_node_value(node_sab, "fraction"));
      }

      // Check that the thermal scattering table is listed in the
      // cross_sections.xml file
      LibraryKey key {Library::Type::thermal, name};
      if (data::library_map.find(key) == data::library_map.end()) {
        fatal_error("Could not find thermal scattering data " + name +
          " in cross_sections.xml file.");
      }

      // Determine index of thermal scattering data in global
      // data::thermal_scatt array
      int index_table;
      if (data::thermal_scatt_map.find(name) == data::thermal_scatt_map.end()) {
        index_table = data::thermal_scatt_map.size();
        data::thermal_scatt_map[name] = index_table;
      } else {
        index_table = data::thermal_scatt_map[name];
      }

      // Add entry to thermal tables vector. For now, we put the nuclide index
      // as zero since we don't know which nuclides the table is being applied
      // to yet (this is assigned in init_thermal)
      thermal_tables_.push_back({index_table, 0, fraction});
    }
  }
}

void Material::finalize()
{
  // Set fissionable if any nuclide is fissionable
  for (const auto& i_nuc : nuclide_) {
    if (data::nuclides[i_nuc]->fissionable_) {
      fissionable_ = true;
      break;
    }
  }

  // Generate material bremsstrahlung data for electrons and positrons
  if (settings::photon_transport && settings::electron_treatment == ELECTRON_TTB) {
    this->init_bremsstrahlung();
  }

  // Assign thermal scattering tables
  this->init_thermal();

  // Normalize density
  this->normalize_density();
}

void Material::normalize_density()
{
  bool percent_in_atom = (atom_density_(1) > 0.0);
  bool density_in_atom = (density_ > 0.0);

  double sum_percent = 0.0;
  for (int i = 0; i < nuclide_.size(); ++i) {
    // determine atomic weight ratio
    int i_nuc = nuclide_[i];
    double awr = settings::run_CE ?
      data::nuclides[i_nuc]->awr_ : data::nuclides_MG[i_nuc].awr;

    // if given weight percent, convert all values so that they are divided
    // by awr. thus, when a sum is done over the values, it's actually
    // sum(w/awr)
    if (!percent_in_atom) atom_density_(i) = -atom_density_(i) / awr;
  }

  // determine normalized atom percents. if given atom percents, this is
  // straightforward. if given weight percents, the value is w/awr and is
  // divided by sum(w/awr)
  atom_density_ /= xt::sum(atom_density_)();

  // Change density in g/cm^3 to atom/b-cm. Since all values are now in
  // atom percent, the sum needs to be re-evaluated as 1/sum(x*awr)
  if (!density_in_atom) {
    sum_percent = 0.0;
    for (int i = 0; i < nuclide_.size(); ++i) {
      int i_nuc = nuclide_[i];
      double awr = settings::run_CE ?
        data::nuclides[i_nuc]->awr_ : data::nuclides_MG[i_nuc].awr;
      sum_percent += atom_density_(i)*awr;
    }
    sum_percent = 1.0 / sum_percent;
    density_ *= -N_AVOGADRO / MASS_NEUTRON * sum_percent;
  }

  // Calculate nuclide atom densities
  atom_density_ *= density_;

  // Calculate density in g/cm^3.
  double density_gpcc_ = 0.0;
  for (int i = 0; i < nuclide_.size(); ++i) {
    int i_nuc = nuclide_[i];
    double awr = settings::run_CE ? data::nuclides[i_nuc]->awr_ : 1.0;
    density_gpcc_ += atom_density_(i) * awr * MASS_NEUTRON / N_AVOGADRO;
  }
}

void Material::init_thermal()
{
  std::vector<ThermalTable> tables;

  for (const auto& table : thermal_tables_) {
    // In order to know which nuclide the S(a,b) table applies to, we need
    // to search through the list of nuclides for one which has a matching
    // name
    bool found = false;
    for (int j = 0; j < nuclide_.size(); ++j) {
      const auto& name {data::nuclides[nuclide_[j]]->name_};
      if (contains(data::thermal_scatt[table.index_table]->nuclides_, name)) {
        tables.push_back({table.index_table, j, table.fraction});
        found = true;
      }
    }

    // Check to make sure thermal scattering table matched a nuclide
    if (!found) {
      fatal_error("Thermal scattering table " + data::thermal_scatt[
        table.index_table]->name_  + " did not match any nuclide on material "
        + std::to_string(id_));
    }
  }

  // Make sure each nuclide only appears in one table.
  for (int j = 0; j < tables.size(); ++j) {
    for (int k = j+1; k < tables.size(); ++k) {
      if (tables[j].index_nuclide == tables[k].index_nuclide) {
        int index = nuclide_[tables[j].index_nuclide];
        auto name = data::nuclides[index]->name_;
        fatal_error(name + " in material " + std::to_string(id_) + " was found "
          "in multiple thermal scattering tables. Each nuclide can appear in "
          "only one table per material.");
      }
    }
  }

  // If there are multiple S(a,b) tables, we need to make sure that the
  // entries in i_sab_nuclides are sorted or else they won't be applied
  // correctly in the cross_section module.
  std::sort(tables.begin(), tables.end(), [](ThermalTable a, ThermalTable b) {
    return a.index_nuclide < b.index_nuclide;
  });

  // Update the list of thermal tables
  thermal_tables_ = tables;
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
        spline(n, &data::ttb_e_grid(i), &f(i), &z(i));

        double c = 0.0;
        for (int j = i; j < n_e - 1; ++j) {
          c += spline_integrate(n, &data::ttb_e_grid(i), &f(i), &z(i),
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

void Material::calculate_xs(const Particle& p) const
{
  // Set all material macroscopic cross sections to zero
  simulation::material_xs.total = 0.0;
  simulation::material_xs.absorption = 0.0;
  simulation::material_xs.fission = 0.0;
  simulation::material_xs.nu_fission = 0.0;

  if (p.type == static_cast<int>(ParticleType::neutron)) {
    this->calculate_neutron_xs(p);
  } else if (p.type == static_cast<int>(ParticleType::photon)) {
    this->calculate_photon_xs(p);
  }
}

void Material::calculate_neutron_xs(const Particle& p) const
{
  // TODO: off-by-one
  int neutron = static_cast<int>(ParticleType::neutron) - 1;

  // Find energy index on energy grid
  int i_grid = std::log(p.E/data::energy_min[neutron])/simulation::log_spacing;

  // Determine if this material has S(a,b) tables
  bool check_sab = (thermal_tables_.size() > 0);

  // Initialize position in i_sab_nuclides
  int j = 0;

  // Add contribution from each nuclide in material
  for (int i = 0; i < nuclide_.size(); ++i) {
    // ======================================================================
    // CHECK FOR S(A,B) TABLE

    int i_sab = C_NONE;
    double sab_frac = 0.0;

    // Check if this nuclide matches one of the S(a,b) tables specified.
    // This relies on thermal_tables_ being sorted by .index_nuclide
    if (check_sab) {
      const auto& sab {thermal_tables_[j]};
      if (i == sab.index_nuclide) {
        // Get index in sab_tables
        i_sab = sab.index_table;
        sab_frac = sab.fraction;

        // If particle energy is greater than the highest energy for the
        // S(a,b) table, then don't use the S(a,b) table
        if (p.E > data::thermal_scatt[i_sab]->threshold()) i_sab = C_NONE;

        // Increment position in i_sab_nuclides
        ++j;

        // Don't check for S(a,b) tables if there are no more left
        if (j == thermal_tables_.size()) check_sab = false;
      }
    }

    // ======================================================================
    // CALCULATE MICROSCOPIC CROSS SECTION

    // Determine microscopic cross sections for this nuclide
    int i_nuclide = nuclide_[i];

    // Calculate microscopic cross section for this nuclide
    const auto& micro {simulation::micro_xs[i_nuclide]};
    if (p.E != micro.last_E
        || p.sqrtkT != micro.last_sqrtkT
        || i_sab != micro.index_sab
        || sab_frac != micro.sab_frac) {
      data::nuclides[i_nuclide]->calculate_xs(i_sab, p.E, i_grid,
        p.sqrtkT, sab_frac);
    }

    // ======================================================================
    // ADD TO MACROSCOPIC CROSS SECTION

    // Copy atom density of nuclide in material
    double atom_density = atom_density_(i);

    // Add contributions to cross sections
    const auto& macro {simulation::material_xs};
    simulation::material_xs.total += atom_density * micro.total;
    simulation::material_xs.absorption += atom_density * micro.absorption;
    simulation::material_xs.fission += atom_density * micro.fission;
    simulation::material_xs.nu_fission += atom_density * micro.nu_fission;
  }
}

void Material::calculate_photon_xs(const Particle& p) const
{
  simulation::material_xs.coherent = 0.0;
  simulation::material_xs.incoherent = 0.0;
  simulation::material_xs.photoelectric = 0.0;
  simulation::material_xs.pair_production = 0.0;

  // Add contribution from each nuclide in material
  for (int i = 0; i < nuclide_.size(); ++i) {
    // ========================================================================
    // CALCULATE MICROSCOPIC CROSS SECTION

    // Determine microscopic cross sections for this nuclide
    int i_element = element_[i];

    // Calculate microscopic cross section for this nuclide
    const auto& micro {simulation::micro_photon_xs[i_element]};
    if (p.E != micro.last_E) {
      data::elements[i_element].calculate_xs(p.E);
    }

    // ========================================================================
    // ADD TO MACROSCOPIC CROSS SECTION

    // Copy atom density of nuclide in material
    double atom_density = atom_density_(i);

    // Add contributions to material macroscopic cross sections
    simulation::material_xs.total += atom_density * micro.total;
    simulation::material_xs.coherent += atom_density * micro.coherent;
    simulation::material_xs.incoherent += atom_density * micro.incoherent;
    simulation::material_xs.photoelectric += atom_density * micro.photoelectric;
    simulation::material_xs.pair_production += atom_density * micro.pair_production;
  }
}

int Material::set_density(double density, std::string units)
{
  if (nuclide_.empty()) {
    set_errmsg("No nuclides exist in material yet.");
    return OPENMC_E_ALLOCATE;
  }

  if (units == "atom/b-cm") {
    // Set total density based on value provided
    density_ = density;

    // Determine normalized atom percents
    double sum_percent = xt::sum(atom_density_)();
    atom_density_ /= sum_percent;

    // Recalculate nuclide atom densities based on given density
    atom_density_ *= density;

    // Calculate density in g/cm^3.
    density_gpcc_ = 0.0;
    for (int i = 0; i < nuclide_.size(); ++i) {
      int i_nuc = nuclide_[i];
      double awr = data::nuclides[i_nuc]->awr_;
      density_gpcc_ += atom_density_(i) * awr * MASS_NEUTRON / N_AVOGADRO;
    }
  } else if (units == "g/cm3" || units == "g/cc") {
    // Determine factor by which to change densities
    double previous_density_gpcc = density_gpcc_;
    double f = density / previous_density_gpcc;

    // Update densities
    density_gpcc_ = density;
    density_ *= f;
    atom_density_ *= f;
  } else {
    set_errmsg("Invalid units '" + units + "' specified.");
    return OPENMC_E_INVALID_ARGUMENT;
  }
  return 0;
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
openmc_material_add_nuclide(int32_t index, const char* name, double density)
{
  int err = 0;
  if (index >= 1 && index <= model::materials.size()) {
    Material* m = model::materials[index - 1];

    // Check if nuclide is already in material
    for (int i = 0; i < m->nuclide_.size(); ++i) {
      int i_nuc = m->nuclide_[i];
      if (data::nuclides[i_nuc]->name_ == name) {
        double awr = data::nuclides[i_nuc]->awr_;
        m->density_ += density - m->atom_density_(i);
        m->density_gpcc_ += (density - m->atom_density_(i))
          * awr * MASS_NEUTRON / N_AVOGADRO;
        m->atom_density_(i) = density;
        return 0;
      }
    }

    // If nuclide wasn't found, extend nuclide/density arrays
    err = openmc_load_nuclide(name);

    if (err == 0) {
      // Append new nuclide/density
      int i_nuc = data::nuclide_map[name];
      m->nuclide_.push_back(i_nuc);

      auto n = m->nuclide_.size();

      // Create copy of atom_density_ array with one extra entry
      xt::xtensor<double, 1> atom_density = xt::zeros<double>({n});
      xt::view(atom_density, xt::range(0, n-1)) = m->atom_density_;
      atom_density(n) = density;
      m->atom_density_ = atom_density;

      m->density_ += density;
      m->density_gpcc_ += density * data::nuclides[i_nuc]->awr_
        * MASS_NEUTRON / N_AVOGADRO;
    }
  } else {
    set_errmsg("Index in materials array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return err;
}

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
openmc_material_set_density(int32_t index, double density, const char* units)
{
  if (index >= 1 && index <= model::materials.size()) {
    return model::materials[index - 1]->set_density(density, units);
  } else {
    set_errmsg("Index in materials array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

extern "C" int
openmc_material_set_densities(int32_t index, int n, const char** name, const double* density)
{
  if (index >= 1 && index <= model::materials.size()) {
    // TODO: off-by-one
    auto& mat {model::materials[index - 1]};
    if (n != mat->nuclide_.size()) {
      mat->nuclide_.resize(n);
      mat->atom_density_ = xt::zeros<double>({n});
    }

    double sum_density = 0.0;
    for (int i = 0; i < n; ++i) {
      std::string nuc {name[i]};
      if (data::nuclide_map.find(nuc) == data::nuclide_map.end()) {
        int err = openmc_load_nuclide(nuc.c_str());
        if (err < 0) return err;
      }

      mat->nuclide_[i] = data::nuclide_map[nuc];
      mat->atom_density_(i) = density[i];
      sum_density += density[i];
    }

    // Set total density to the sum of the vector

    int err = mat->set_density(sum_density, "atom/b-cm");

    // Assign S(a,b) tables
    mat->init_thermal();
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

  bool material_fissionable(Material* mat) {return mat->fissionable_;}

  void material_set_fissionable(Material* mat, bool fissionable)
  {
    mat->fissionable_ = fissionable;
  }

  void material_calculate_xs_c(Material* mat, const Particle* p)
  {
    mat->calculate_xs(*p);
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
