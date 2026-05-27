#include <fstream>

#include <fmt/core.h>

#include "openmc/constants.h"
#include "openmc/majorant.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/search.h"
#include "openmc/simulation.h"
#include "openmc/thermal.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace data {
std::unique_ptr<NeutronMajorant> n_majorant;
std::string majorant_file;

} // namespace data

//==============================================================================
// Majorant implementation
//==============================================================================

Majorant::Majorant(const std::string & majorant_file, double min_E_transport, double max_E_transport)
{
  std::ifstream majorant_data(majorant_file);

  std::string line;
  while (std::getline(majorant_data, line)) {
    auto delim_pos = line.find(",");
    auto energy = std::stod(line.substr(0, delim_pos));
    if (energy < min_E_transport) {
      continue;
    }
    if (energy > max_E_transport) {
      break;
    }
    grid_.energy.push_back(energy);
    xs_.push_back(std::stod(line.substr(delim_pos + 1)));
  }

  grid_.init();
}

void
Majorant::write_ascii(const std::string& filename) const
{
  std::ofstream of(filename);
  for (int i = 0; i < xs_.size(); i++) {
    of << grid_.energy[i] << "\t" << xs_[i] << "\n";
  }
  of.close();
}

double
Majorant::interpolate_lin_1D(double x_0, double x_1, double y_0, double y_1, double x) const
{
  double f = (x - x_0) / (x_1 - x_0);
  return (1.0 - f) * y_0 + f * y_1;
}

double
Majorant::interpolate_log_1D(double x_0, double x_1, double y_0, double y_1, double x) const
{
  double f = std::log(x / x_0) / std::log(x_1 / x_0);
  return std::exp((1.0 - f) * std::log(y_0) + f * std::log(y_1));
}

//==============================================================================
// NeutronMajorant implementation
//==============================================================================

NeutronMajorant::NeutronMajorant(const std::string & majorant_file)
  : Majorant(majorant_file,
             data::energy_min[ParticleType::neutron().transport_index()],
             data::energy_max[ParticleType::neutron().transport_index()])
{ }

double
NeutronMajorant::calculate_neutron_xs(double energy) const
{
  int i_grid = get_i_grid(energy, grid_);

  // calculate interpolation factor
  double f = (energy - grid_.energy[i_grid]) /
              (grid_.energy[i_grid + 1]- grid_.energy[i_grid]);

  double xs = (1.0 - f) * xs_[i_grid] + f * xs_[i_grid + 1];

  return xs;
}

void
NeutronMajorant::init()
{
  // Unionize the grid.
  compute_unionized_grid();

  // Setup the majorant given the new grid.
  setup_majorant();
}

void
NeutronMajorant::compute_unionized_grid()
{
  write_message("Unionizing nuclide cross section grids.");

  // This function generates a unionized cross section grid between smooth cross
  // sections and URR probability table grids.
  for (const auto & mat : model::materials) {
    for (auto nuclide_idx : mat->nuclide_) {
      const auto & nuclide = data::nuclides[nuclide_idx];

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
        grid_.energy.insert(grid_.energy.end(), n_grid.energy.begin(),n_grid.energy.end());
      }
    }
  }
  std::sort(grid_.energy.begin(), grid_.energy.end());
  std::unique(grid_.energy.begin(), grid_.energy.end());

  int neutron = ParticleType::neutron().transport_index();
  // remove all values below the minimum neutron energy
  auto min_it = grid_.energy.begin();
  while (*min_it < data::energy_min[neutron]) { min_it++; }
  grid_.energy.erase(grid_.energy.begin(), min_it + 1);
  // insert the minimum neutron energy at the beginning
  grid_.energy.insert(grid_.energy.begin(), data::energy_min[neutron]);

  // remove all values above the maximum neutron energy
  auto max_it = --grid_.energy.end();
  while (*max_it > data::energy_max[neutron]) { max_it--; }
  grid_.energy.erase(max_it - 1, grid_.energy.end());
  // insert the maximum neutron energy at the end
  grid_.energy.insert(grid_.energy.end(), data::energy_max[neutron]);

  grid_.init();
}

void
NeutronMajorant::setup_majorant()
{
  // Fill with zeros.
  xs_.resize(grid_.energy.size(), 0.0);

  std::vector<double> material_maj_xs;
  material_maj_xs.resize(grid_.energy.size(), 0.0);
  for (const auto & mat : model::materials) {
    write_message("Computing majorant total cross section for " + mat->name() + ".");

    // Populate the per-material majorant cross section.
    fill_material_maj_xs((*mat.get()), grid_.energy, material_maj_xs);

    // Compute the full majorant by taking the max over each material cross section.
    for (int i_energy = 0; i_energy < xs_.size(); ++i_energy) {
      xs_[i_energy] = std::max(xs_[i_energy], material_maj_xs[i_energy]);
    }
    std::fill(material_maj_xs.begin(), material_maj_xs.end(), 0.0);
  }
}

void
NeutronMajorant::fill_material_maj_xs(const Material & mat, const std::vector<double> & to_grid, std::vector<double> & mat_maj)
{
  for (int i_energy = 0; i_energy < to_grid.size(); ++i_energy) {
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
        micro_smooth_tot_xs = calculate_max_sab_tot_xs(union_energy, i_sab, sab_frac, (*nuclide.get()));
      } else {
        // Free gas smooth cross section
        micro_smooth_tot_xs = calculate_max_smooth_xs(union_energy, (*nuclide.get()));
      }

      // ======================================================================
      // Compute the URR cross section. This shouldn't intersect with the
      // S(a,b) cross section.
      double micro_urr_xs = calculate_max_urr_xs(union_energy, (*nuclide.get()), micro_smooth_tot_xs);

      // ======================================================================
      // Accumulate the macroscopic cross section.
      // TODO: density multipliers for per-cell material densities.
      mat_maj[i_energy] += std::max(micro_smooth_tot_xs, micro_urr_xs) * mat.atom_density(i);
    }
  }
}

double
NeutronMajorant::calculate_max_smooth_xs(double energy, const Nuclide & nuc) const
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

double
NeutronMajorant::calculate_max_urr_xs(double energy, const Nuclide & nuc, double smooth_xs) const
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

double
NeutronMajorant::calculate_max_sab_tot_xs(double energy, int i_sab, double sab_frac, const Nuclide & nuc) const
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

int
NeutronMajorant::get_i_grid(double energy, const Nuclide::EnergyGrid & grid) const
{
  // Find energy index on energy grid
  int neutron = ParticleType::neutron().transport_index();
  int i_log_union = std::log(energy * data::energy_min_rcp[neutron]) * simulation::log_spacing_rcp;

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
// Static functions
//==============================================================================

void create_majorants()
{
  if (data::majorant_file != "") {
    write_message("Loading majorant from " + data::majorant_file);
    // We can load the majorant from a file instead.
    data::n_majorant = std::make_unique<NeutronMajorant>(data::majorant_file);
    return;
  }

  write_message("Creating majorant cross section...");
  data::n_majorant = std::make_unique<NeutronMajorant>();
  data::n_majorant->init();
  data::n_majorant->write_ascii("macro_majorant.txt");
}

} // namespace openmc
