#include <fstream>

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/geometry.h"
#include "openmc/majorant.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/search.h"
#include "openmc/simulation.h"
#include "openmc/settings.h"
#include "openmc/thermal.h"
#include "openmc/universe.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {
std::unique_ptr<NeutronMajorant> n_majorant;
std::unique_ptr<PhotonMajorant> p_majorant;
std::string n_majorant_file;
std::string p_majorant_file;

} // namespace data

//==============================================================================
// Majorant implementation
//==============================================================================

Majorant::Majorant(int i_universe)
  : maj_universe_(i_universe)
{
  // Find all materials contained in the majorant's universe. This also obtains
  // the maximum density multiplier applied to that material.
  std::set<int> unique_materials;
  if (maj_universe_ == C_NONE || maj_universe_ >= model::universes.size()) {
    fatal_error( fmt::format("Invalid majorant universe: {}", maj_universe_));
  }

  const auto & maj_uni = model::universes[maj_universe_];
  for (int i_cell : maj_uni->cells_) {
    const auto & cell = model::cells[i_cell];

    // If the cell is filled with a material, it won't have any sub-cells.
    if (cell->type_ == Fill::MATERIAL) {
      // Loop over instances. TODO: confirm if this is unecessary and use 0 instead?
      for (int instance = 0; instance < cell->n_instances(); ++instance) {
        int i_material = cell->material(instance);

        // Check to see if we've found the contained material yet. If not, add to the set
        // of materials discovered and add to the map of density multipliers.
        if (unique_materials.count(i_material) == 0) {
          unique_materials.emplace(i_material);
          max_density_mult_[i_material] = cell->density_mult(instance);
        } else {
          // We've found this material already. Need to take the maximum density multiplier.
          max_density_mult_.at(i_material) =
            std::max(max_density_mult_.at(i_material),
                     cell->density_mult(instance));
        }
      }
    } else {
      // This cell is filled with a universe or lattice. Need to get the list of cells and
      // cell instances.
      const auto contained_cells = cell->get_contained_cells();
      for (const auto & [i_con_cell, contained_instances] : contained_cells) {
        const auto & contained_cell = model::cells[i_con_cell];

        // Loop over contained cell instances.
        for (auto instance : contained_instances) {
          // Check to see if we've found the contained material instance yet. If not, add
          // to the set of materials discovered and add to the map of density multipliers.
          int i_material = contained_cell->material(instance);
          if (unique_materials.count(i_material) == 0) {
            unique_materials.emplace(i_material);
            max_density_mult_[i_material] = cell->density_mult(instance);
          } else {
            // We've found this material already. Need to take the maximum density multiplier
            // for the contained instance.
            max_density_mult_.at(i_material) =
              std::max(max_density_mult_.at(i_material), cell->density_mult(instance));
          }
        }
      }
    }
  }

  // Clear the contained materials vector and insert the elements from the set.
  contained_materials_.clear();
  for (auto i_mat : unique_materials) {
    contained_materials_.push_back(i_mat);
  }
}

void Majorant::compute_majorant()
{
  // Fill with zeros.
  xs_.resize(grid_.energy.size(), 0.0);

  std::vector<double> material_maj_xs;
  material_maj_xs.resize(grid_.energy.size(), 0.0);
  for (int i_material : contained_materials_) {
    // Populate the per-material majorant cross section.
    fill_material_maj_xs(*model::materials[i_material], max_density_mult_.at(i_material),
      grid_.energy, material_maj_xs);

    // Compute the full majorant by taking the max over each material cross section.
    for (int i_energy = 0; i_energy < xs_.size(); ++i_energy) {
      xs_[i_energy] = std::max(xs_[i_energy], material_maj_xs[i_energy]);
    }
    std::fill(material_maj_xs.begin(), material_maj_xs.end(), 0.0);
  }
}

void Majorant::write_ascii(const std::string& filename) const
{
  std::ofstream of(filename);
  for (int i = 0; i < xs_.size(); i++) {
    of << grid_.energy[i] << "," << xs_[i] << "\n";
  }
  of.close();
}

//==============================================================================
// NeutronMajorant implementation
//==============================================================================

NeutronMajorant::NeutronMajorant(int i_universe)
  : Majorant(i_universe)
{ }

double NeutronMajorant::calculate_neutron_xs(double energy) const
{
  const int i_grid = get_i_grid(energy, grid_);
  return interpolate_lin_1D(grid_.energy[i_grid], grid_.energy[i_grid + 1], xs_[i_grid], xs_[i_grid + 1], energy);
}

void NeutronMajorant::compute_unionized_grid()
{
  // In the event the majorant needs to be re-generated (e.g. in-memory for
  // multiphysics), we need to reset the unionized grid.
  grid_.energy.clear();

  // This function generates a unionized cross section grid between smooth cross
  // sections and URR probability table grids.
  std::set<int> processed_nuclides;
  for (int i_mat : contained_materials_) {
    const auto & mat = model::materials[i_mat];
    for (auto i_nuclide : mat->nuclide_) {
      // Only unionize nuclides we haven't checked yet.
      if (processed_nuclides.count(i_nuclide) > 0) {
        continue;
      }

      const auto & nuclide = data::nuclides[i_nuclide];
      // ======================================================================
      // Unionizing the URR temperature grid.
      if (nuclide->urr_present_ && settings::urr_ptables_on) {
        for (const auto & nuc_urr : nuclide->urr_data_) {
          grid_.energy.insert(grid_.energy.end(), nuc_urr.energy_.begin(), nuc_urr.energy_.end());
        }
      }

      // ======================================================================
      // Unionize the smooth cross section grid.
      for (const auto & n_grid : nuclide->grid_) {
        grid_.energy.insert(grid_.energy.end(), n_grid.energy.begin(), n_grid.energy.end());
      }

      processed_nuclides.insert(i_nuclide);
    }
  }
  std::sort(grid_.energy.begin(), grid_.energy.end());
  auto unique_end = std::unique(grid_.energy.begin(), grid_.energy.end());
  grid_.energy.resize(std::distance(grid_.energy.begin(), unique_end));

  // remove all values below the minimum neutron energy
  auto min_it = grid_.energy.begin();
  while (*min_it < data::energy_min[i_neutron_]) { min_it++; }
  grid_.energy.erase(grid_.energy.begin(), min_it + 1);
  // insert the minimum neutron energy at the beginning
  grid_.energy.insert(grid_.energy.begin(), data::energy_min[i_neutron_]);

  // remove all values above the maximum neutron energy
  auto max_it = --grid_.energy.end();
  while (*max_it > data::energy_max[i_neutron_]) { max_it--; }
  grid_.energy.erase(max_it - 1, grid_.energy.end());
  // insert the maximum neutron energy at the end
  grid_.energy.insert(grid_.energy.end(), data::energy_max[i_neutron_]);

  // Initialize the grid for fast lookups.
  grid_.init();
}

void NeutronMajorant::fill_material_maj_xs(const Material & mat, double max_density_mult,
  const std::vector<double> & to_grid, std::vector<double> & mat_maj)
{
  for (int i_energy = 0; i_energy < to_grid.size(); ++i_energy) {
    mat_maj[i_energy] = 0.0;
    const double union_energy = to_grid[i_energy];

    int mat_sab_table_idx = 0;
    bool check_sab = (mat.thermal_tables_.size() > 0);

    for (int i = 0; i < mat.nuclide_.size(); ++i) {
      const auto & nuclide = data::nuclides[mat.nuclide_[i]];
      // ======================================================================
      // CHECK FOR S(A,B) TABLE
      int i_sab = C_NONE;
      double sab_frac = 0.0;

      // Check if this nuclide matches one of the S(a,b) tables specified.
      // This relies on thermal_tables_ being sorted by .index_nuclide
      if (check_sab) {
        const auto& sab {mat.thermal_tables_[mat_sab_table_idx]};
        if (i == sab.index_nuclide) {
          // Get index in sab_tables
          i_sab = sab.index_table;
          sab_frac = sab.fraction;

          // If particle energy is greater than the highest energy for the
          // S(a,b) table, then don't use the S(a,b) table
          if ((union_energy - 1e-6) > data::thermal_scatt[i_sab]->energy_max_) {
            i_sab = C_NONE;
          }

          // Increment position in thermal_tables_
          ++mat_sab_table_idx;

          // Don't check for S(a,b) tables if there are no more left
          if (mat_sab_table_idx == mat.thermal_tables_.size()) {
            check_sab = false;
          }
        }
      }

      // ======================================================================
      // Compute the maximum smooth total cross section. This is either the
      // free gas cross section at energies larger than the Bragg edge, or
      // the bound cross section in the thermal scattering region.
      double micro_smooth_tot_xs = 0.0;
      if (i_sab >= 0) {
        // Thermal scattering cross sections using S(a,b) tables.
        micro_smooth_tot_xs = calculate_max_sab_tot_xs(union_energy, i_sab, sab_frac, *nuclide);
      } else {
        // Free gas smooth cross section
        micro_smooth_tot_xs = calculate_max_smooth_xs(union_energy, *nuclide);
      }

      // ======================================================================
      // Compute the URR cross section. This shouldn't intersect with the
      // S(a,b) cross section.
      double micro_urr_xs = calculate_max_urr_xs(union_energy, *nuclide, micro_smooth_tot_xs);

      // ======================================================================
      // Accumulate the macroscopic cross section.
      // TODO: density multipliers for per-cell material densities.
      mat_maj[i_energy] += std::max(micro_smooth_tot_xs, micro_urr_xs) * mat.atom_density(i, max_density_mult);
    }
  }
}

double NeutronMajorant::calculate_max_smooth_xs(double energy, const Nuclide & nuc) const
{
  double max_smooth_tot_xs = 0.0;
  for (int i_temp = 0; i_temp < nuc.kTs_.size(); ++i_temp) {
    const auto & nuc_grid = nuc.grid_[i_temp];
    int i_grid = get_i_grid(energy, nuc_grid);
    auto total = nuc.xs_[i_temp].slice(openmc::tensor::all, 0);
    double xs = interpolate_lin_1D(nuc_grid.energy[i_grid], nuc_grid.energy[i_grid + 1], total[i_grid], total[i_grid + 1], energy);
    max_smooth_tot_xs = std::max(max_smooth_tot_xs, xs);
  }

  return max_smooth_tot_xs;
}

double NeutronMajorant::calculate_max_urr_xs(double energy, const Nuclide & nuc, double smooth_xs) const
{
  if (!nuc.urr_present_) {
    return 0.0;
  }

  double max_urr_xs = 0.0;
  for (const auto & urr : nuc.urr_data_) {
    if (!(urr.energy_in_bounds(energy - 1e-6) || urr.energy_in_bounds(energy + 1e-6))) {
      continue;
    }

    int i_energy = lower_bound_index(&urr.energy_.front(), &urr.energy_.back(), energy);

    // Find the maximum URR cross sections for the two bounding energy points.
    double max_urr_xs_E0 = 0.0;
    double max_urr_xs_E1 = 0.0;
    for (int i_cdf = 0; i_cdf < urr.n_cdf(); ++i_cdf) {
      max_urr_xs_E0 = std::max(max_urr_xs_E0, urr.xs_values_(i_energy, i_cdf).total);
      max_urr_xs_E1 = std::max(max_urr_xs_E1, urr.xs_values_(i_energy + 1, i_cdf).total);
    }

    // Interpolate the bounding energy points.
    double interp_urr_xs = 0.0;
    if (urr.interp_ == Interpolation::lin_lin) {
      interp_urr_xs =
        interpolate_lin_1D(urr.energy_[i_energy], urr.energy_[i_energy + 1], max_urr_xs_E0, max_urr_xs_E1, energy);
    } else if (urr.interp_ == Interpolation::log_log) {
      interp_urr_xs =
        interpolate_log_1D(urr.energy_[i_energy], urr.energy_[i_energy + 1], max_urr_xs_E0, max_urr_xs_E1, energy);
    }

    // Multiply by the smooth cross section (after interpolation) if required.
    if (urr.multiply_smooth_) {
      interp_urr_xs *= smooth_xs;
    }

    max_urr_xs = std::max(max_urr_xs, interp_urr_xs);
  }

  return max_urr_xs;
}

double NeutronMajorant::calculate_max_sab_tot_xs(double energy, int i_sab, double sab_frac, const Nuclide & nuc) const
{
  const auto & thermal = data::thermal_scatt[i_sab];

  // Loop over the nuclide's temperature grid to ensure we're consistent.
  double max_sab_total = 0.0;
  for (int i_nuc_temp = 0; i_nuc_temp < nuc.kTs_.size(); ++i_nuc_temp) {
    double nuc_kT = nuc.kTs_[i_nuc_temp] * nuc.kTs_[i_nuc_temp];

    // Compute the elastic and inelastic scattering cross sections. The S(a,b)
    // cross sections are interpolated to match the nuclide temperature point.
    double thermal_elastic;
    double thermal_inelastic;
    const auto & tkTs = thermal->kTs_;
    if (tkTs.size() > 1) {
      if (nuc_kT < tkTs.front()) {
        thermal->data_.front().calculate_xs(energy, &thermal_elastic, &thermal_inelastic);
      } else if (nuc_kT > tkTs.back()) {
        thermal->data_.back().calculate_xs(energy, &thermal_elastic, &thermal_inelastic);
      } else {
        // Find temperatures that bound the actual temperature
        int i_sab_temp = 0;
        while (tkTs[i_sab_temp + 1] < nuc_kT && i_sab_temp + 1 < tkTs.size() - 1) {
          ++i_sab_temp;
        }
        // Interpolate the scattering cross sections to the nuclide temperature grid point.
        double T0_elastic, T1_elastic, T0_inelastic, T1_inelastic;
        thermal->data_[i_sab_temp].calculate_xs(energy, &T0_elastic, &T0_inelastic);
        thermal->data_[i_sab_temp + 1].calculate_xs(energy, &T1_elastic, &T1_inelastic);
        thermal_elastic = interpolate_lin_1D(tkTs[i_sab_temp], tkTs[i_sab_temp + 1], T0_elastic, T1_elastic, nuc_kT);
        thermal_inelastic = interpolate_lin_1D(tkTs[i_sab_temp], tkTs[i_sab_temp + 1], T0_inelastic, T1_inelastic, nuc_kT);
      }
    } else {
      thermal->data_[0].calculate_xs(energy, &thermal_elastic, &thermal_inelastic);
    }

    // Compute the free gas total and elastic cross sections interpolated on the majorant grid.
    const auto & nuc_grid = nuc.grid_[i_nuc_temp];
    int i_grid = get_i_grid(energy, nuc_grid);
    const auto & free_tot = nuc.xs_[i_nuc_temp].slice(openmc::tensor::all, 0);
    const auto & free_ela = nuc.reactions_[0]->xs_[i_nuc_temp].value;
    double tot_xs = interpolate_lin_1D(nuc_grid.energy[i_grid], nuc_grid.energy[i_grid + 1], free_tot[i_grid], free_tot[i_grid + 1], energy);
    double ela_xs = interpolate_lin_1D(nuc_grid.energy[i_grid], nuc_grid.energy[i_grid + 1], free_ela[i_grid], free_ela[i_grid + 1], energy);

    double thermal_xs = sab_frac * (thermal_elastic + thermal_inelastic);
    double sab_corrected_total = tot_xs + thermal_xs - sab_frac * ela_xs;
    max_sab_total = std::max(sab_corrected_total, max_sab_total);
  }

  return max_sab_total;
}

int NeutronMajorant::get_i_grid(double energy, const Nuclide::EnergyGrid & grid) const
{
  // Find energy index on energy grid
  int i_log_union = std::log(energy * data::energy_min_rcp[i_neutron_]) * simulation::log_spacing_rcp;

  int i_grid;
  if (i_log_union < 0) {
    i_grid = 0;
  } else if (i_log_union >= (grid.grid_index.size() - 2)) {
    i_grid = grid.energy.size() - 2;
  } else {
    // Determine bounding indices based on which equal log-spaced
    // interval the energy is in
    int i_low  = grid.grid_index[i_log_union];
    int i_high = grid.grid_index[i_log_union + 1] + 1;

    // Perform binary search over reduced range
    i_grid = i_low + lower_bound_index(&grid.energy[i_low], &grid.energy[i_high], energy);
  }

  // check for rare case where two energy points are the same
  if (grid.energy[i_grid] == grid.energy[i_grid + 1]) ++i_grid;

  return i_grid;
}

//==============================================================================
// PhotonMajorant implementation
//==============================================================================

PhotonMajorant::PhotonMajorant(int i_universe)
  : Majorant(i_universe)
{ }

void PhotonMajorant::compute_unionized_grid()
{
  // In the event the majorant needs to be re-generated (e.g. in-memory for
  // multiphysics), we need to reset the unionized grid.
  grid_.energy.clear();

  // This function generates a unionized cross section grid for all elements.
  std::set<int> processed_elements;
  for (int i_mat : contained_materials_) {
    const auto & mat = model::materials[i_mat];
    for (int i = 0; i < mat->nuclide_.size(); ++i) {
      // Only unionize elements we haven't checked yet.
      if (processed_elements.count(mat->element_[i]) > 0) {
        continue;
      }

      const auto & element  = data::elements[mat->element_[i]];
      grid_.energy.insert(grid_.energy.end(), element->energy_.begin(), element->energy_.end());

      processed_elements.insert(mat->element_[i]);
    }
  }
  std::sort(grid_.energy.begin(), grid_.energy.end());
  auto unique_end = std::unique(grid_.energy.begin(), grid_.energy.end());
  grid_.energy.resize(std::distance(grid_.energy.begin(), unique_end));

  // remove all values below the minimum photon energy
  auto min_it = grid_.energy.begin();
  while (*min_it < std::log(data::energy_min[i_photon_])) { min_it++; }
  grid_.energy.erase(grid_.energy.begin(), min_it + 1);
  // insert the minimum photon energy at the beginning
  grid_.energy.insert(grid_.energy.begin(), std::log(data::energy_min[i_photon_]));

  // remove all values above the maximum photon energy
  auto max_it = --grid_.energy.end();
  while (*max_it > std::log(data::energy_max[i_photon_])) { max_it--; }
  grid_.energy.erase(max_it - 1, grid_.energy.end());
  // insert the maximum photon energy at the end
  grid_.energy.insert(grid_.energy.end(), std::log(data::energy_max[i_photon_]));
}

double PhotonMajorant::calculate_photon_xs(double energy) const
{
  double log_energy = std::log(energy);
  int i_grid = get_i_grid(log_energy, grid_.energy);

  // calculate interpolation factor
  double f =
    (log_energy - grid_.energy[i_grid]) / (grid_.energy[i_grid + 1] - grid_.energy[i_grid]);

  // interpolate the total cross section
  return std::exp(xs_[i_grid] + f * (xs_[i_grid + 1] - xs_[i_grid]));
}

void PhotonMajorant::fill_material_maj_xs(const Material & mat, double max_density_mult,
  const std::vector<double> & to_grid, std::vector<double> & mat_maj)
{
  for (int i_energy = 0; i_energy < to_grid.size(); ++i_energy) {
    mat_maj[i_energy] = 0.0;
    const double union_log_energy = to_grid[i_energy];

    for (int i = 0; i < mat.nuclide_.size(); ++i) {
      const int i_element = mat.element_[i];

      mat_maj[i_energy] += calculate_elem_tot_xs(union_log_energy, *data::elements[i_element]) * mat.atom_density(i, max_density_mult);
    }
    mat_maj[i_energy] = std::log(mat_maj[i_energy]);
  }
}

double PhotonMajorant::calculate_elem_tot_xs(double log_energy, const PhotonInteraction & elem) const
{
  int i_grid = get_i_grid(log_energy, elem.energy_);

  // calculate interpolation factor
  double f =
    (log_energy - elem.energy_(i_grid)) / (elem.energy_(i_grid + 1) - elem.energy_(i_grid));

  // Calculate microscopic coherent cross section
  double coherent = std::exp(
    elem.coherent_(i_grid) + f * (elem.coherent_(i_grid + 1) - elem.coherent_(i_grid)));

  // Calculate microscopic incoherent cross section
  double incoherent = std::exp(
    elem.incoherent_(i_grid) + f * (elem.incoherent_(i_grid + 1) - elem.incoherent_(i_grid)));

  // Calculate microscopic photoelectric cross section
  double photoelectric = 0.0;
  tensor::View<const double> xs_lower = elem.cross_sections_.slice(i_grid);
  tensor::View<const double> xs_upper = elem.cross_sections_.slice(i_grid + 1);

  for (int i = 0; i < xs_upper.size(); ++i)
    if (xs_lower(i) != 0)
      photoelectric +=
        std::exp(xs_lower(i) + f * (xs_upper(i) - xs_lower(i)));

  // Calculate microscopic pair production cross section
  double pair_production = std::exp(
    elem.pair_production_total_(i_grid) +
    f * (elem.pair_production_total_(i_grid + 1) - elem.pair_production_total_(i_grid)));

  // Calculate microscopic total cross section
  double total =
    coherent + incoherent + photoelectric + pair_production;
  return total;
}

int PhotonMajorant::get_i_grid(double log_energy, const std::vector<double> & energy_grid) const
{
  int n_grid = energy_grid.size();
  int i_grid;
  if (log_energy <= energy_grid[0]) {
    i_grid = 0;
  } else if (log_energy > energy_grid[n_grid - 1]) {
    i_grid = n_grid - 2;
  } else {
    // We use upper_bound_index here because sometimes photons are created with
    // energies that exactly match a grid point
    i_grid = upper_bound_index(energy_grid.cbegin(), energy_grid.cend(), log_energy);
  }

  // check for case where two energy points are the same
  if (energy_grid[i_grid] == energy_grid[i_grid + 1])
    ++i_grid;

  return i_grid;
}

int PhotonMajorant::get_i_grid(double log_energy, const tensor::Tensor<double> & energy_grid) const
{
  int n_grid = energy_grid.size();
  int i_grid;
  if (log_energy <= energy_grid[0]) {
    i_grid = 0;
  } else if (log_energy > energy_grid(n_grid - 1)) {
    i_grid = n_grid - 2;
  } else {
    // We use upper_bound_index here because sometimes photons are created with
    // energies that exactly match a grid point
    i_grid = upper_bound_index(energy_grid.cbegin(), energy_grid.cend(), log_energy);
  }

  // check for case where two energy points are the same
  if (energy_grid(i_grid) == energy_grid(i_grid + 1))
    ++i_grid;

  return i_grid;
}

//! Create/load a majorant cross section for photons or neutrons. Errors if they
//  exist already.
void create_majorants()
{
  write_message("Creating a neutron majorant cross section");
  data::n_majorant = std::make_unique<NeutronMajorant>(model::root_universe);
  data::n_majorant->compute_unionized_grid();
  data::n_majorant->compute_majorant();
  data::n_majorant->write_ascii("neutron_majorant.txt");

  if (settings::photon_transport) {
    write_message("Creating a photon majorant cross section");
    data::p_majorant = std::make_unique<PhotonMajorant>(model::root_universe);
    data::p_majorant->compute_unionized_grid();
    data::p_majorant->compute_majorant();
    data::p_majorant->write_ascii("photon_majorant.txt");
  }
}

//! Reset the photon and neutron majorant cross sections.
void reset_majorants()
{
  openmc::data::n_majorant.reset(nullptr);
  openmc::data::p_majorant.reset(nullptr);
}
} // namespace openmc
