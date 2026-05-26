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

namespace data {
std::vector<std::unique_ptr<Majorant>> nuclide_majorants;
std::unique_ptr<Majorant> n_majorant;
std::string majorant_file;
}

void create_majorant()
{
  write_message("Creating majorant cross section...");
  if (data::majorant_file != "") {
    write_message("Loading majorant from " + data::majorant_file);
    // We can load the majorant from a file instead.
    data::n_majorant = std::make_unique<Majorant>(data::majorant_file);
    return;
  }

  data::n_majorant = std::make_unique<Majorant>();
  fatal_error("The rest has not been implemented!");

  // data::n_majorant->grid_.init();
  // data::n_majorant->write_ascii("macro_majorant.txt");
}

Majorant::Majorant()
{
  // Unionize the grid.
  compute_unionized_grid();

  // Setup the majorant given the new grid.
  setup_majorant();
}

Majorant::Majorant(const std::vector<double>& energy,
                   const std::vector<double>& xs) : xs_(xs)
{
  grid_.energy = energy;
  grid_.init();
}

Majorant::Majorant(const std::string & majorant_file)
{
  const auto neutron_idx = ParticleType::neutron().transport_index();
  std::ifstream majorant_data(majorant_file);

  std::string line;
  while (std::getline(majorant_data, line)) {
    auto delim_pos = line.find(",");
    auto energy = std::stod(line.substr(0, delim_pos));
    if (energy < data::energy_min[neutron_idx]) {
      continue;
    }
    if (energy > data::energy_max[neutron_idx]) {
      break;
    }
    grid_.energy.push_back(energy);
    xs_.push_back(std::stod(line.substr(delim_pos + 1)));
  }

  grid_.init();
}

double
Majorant::calculate_xs(double energy) const
{
  // Find energy index on energy grid
  int neutron = ParticleType::neutron().transport_index();
  int i_log_union = std::log(energy * data::energy_min_rcp[neutron]) * simulation::log_spacing_rcp;

  int i_grid;
  if (i_log_union < 0) {
    i_grid = 0;
  } else if (i_log_union >= (grid_.grid_index.size() - 2)) {
    i_grid = grid_.energy.size() - 2;
  } else {
    // Determine bounding indices based on which equal log-spaced
    // interval the energy is in
    int i_low  = grid_.grid_index[i_log_union];
    int i_high = grid_.grid_index[i_log_union + 1] + 1;

    // Perform binary search over reduced range
    i_grid = i_low + lower_bound_index(&grid_.energy[i_low], &grid_.energy[i_high], energy);
  }

  // check for rare case where two energy points are the same
  if (grid_.energy[i_grid] == grid_.energy[i_grid + 1]) ++i_grid;

  // calculate interpolation factor
  double f = (energy - grid_.energy[i_grid]) /
              (grid_.energy[i_grid + 1]- grid_.energy[i_grid]);

  double xs = (1.0 - f) * xs_[i_grid] + f * xs_[i_grid + 1];

  return xs;
}

void Majorant::write_ascii(const std::string& filename) const
{
  std::ofstream of(filename);
  for (int i = 0; i < xs_.size(); i++) {
    of << grid_.energy[i] << "\t" << xs_[i] << "\n";
  }
  of.close();
}

void
Majorant::compute_unionized_grid()
{
  write_message("Unionizing nuclide cross section grids.");

  // This function generates a unionized cross section grid between smooth cross
  // sections and URR probability table grids.
  std::vector<double> grid_copy;
  for (const auto & mat : model::materials) {
    if (mat->ncrystal_mat_) {
      fatal_error("Delta tracking is not supported when using NCrystal!");
    }

    for (auto nuclide_idx : mat->nuclide_) {
      const auto & nuclide = data::nuclides[nuclide_idx];

      // ======================================================================
      // Unionizing the URR temperature grid. Loop over temperature points.
      if (nuclide->urr_present_ && settings::urr_ptables_on) {
        for (const auto & nuc_urr : nuclide->urr_data_) {
          grid_copy = grid_.energy;
          vector_union_1D(grid_copy, nuc_urr.energy_, grid_.energy);
        }
      }

      // ======================================================================
      // Unionize the smooth cross section grid. Loop over temperature points.
      for (const auto & n_grid : nuclide->grid_) {
        grid_copy = grid_.energy;
        vector_union_1D(grid_copy, n_grid.energy, grid_.energy);
      }
    }
  }
  grid_.init();
}

void
Majorant::setup_majorant()
{
  // Fill with zeros.
  xs_.resize(grid_.energy.size(), 0.0);

  std::vector<double> material_maj;
  for (const auto & mat : model::materials) {
    write_message("Computing majorant total cross section for " + mat->name() + ".");
    material_maj.clear();
    material_maj.resize(grid_.energy.size(), 0.0);

    // Populate the per-material majorant cross section.
    fill_material_maj_xs((*mat.get()), grid_.energy, material_maj);
  }
}

void
Majorant::fill_material_maj_xs(const Material & mat, const std::vector<double> & to_grid, std::vector<double> & mat_maj)
{
  bool check_sab = (mat.thermal_tables_.size() > 0);

  for (int e_idx = 0; e_idx < to_grid.size(); ++e_idx) {
    for (auto nuclide_idx : mat.nuclide_) {
      const auto & nuclide = data::nuclides[nuclide_idx];
      // ======================================================================
      // Start with the smooth cross section.
      double micro_smooth_xs = calculate_max_smooth_xs(to_grid[e_idx], (*nuclide.get()));

      // ======================================================================
      // Compute the URR cross section.
      double urr_xs = calculate_max_urr_xs(to_grid[e_idx], (*nuclide.get()), micro_smooth_xs);

      // ======================================================================
      // Compute the S(a,b) cross section.
    }
  }
}

double
Majorant::calculate_max_smooth_xs(double energy, const Nuclide & nuc)
{
  double max_smooth_xs = 0.0;
  for (int i = 0; i < nuc.kTs_.size(); ++i) {
    int i_grid = get_i_grid(energy, nuc.grid_[i]);
    double xs = interpolate_lin_1D(grid_.energy[i_grid], grid_.energy[i_grid + 1], xs_[i_grid], xs_[i_grid + 1], energy);
    max_smooth_xs = std::max(max_smooth_xs, xs);
  }

  return max_smooth_xs;
}

double
Majorant::calculate_max_urr_xs(double energy, const Nuclide & nuc, double smooth)
{
  if (!nuc.urr_present_) {
    return 0.0;
  }

  double max_urr_xs = 0.0;
  for (const auto & urr : nuc.urr_data_) {
    if (!urr.energy_in_bounds(energy)) {
      continue;
    }
    // Linear search to find the left and right interpolation points.
    // TODO: if this is a problem we can perform binary search in energy
    int i_energy;
    for (int i = 0; i < urr.energy_.size() - 1; ++i) {
      if (urr.energy_[i] <= energy && energy < urr.energy_[i + 1]) {
        i_energy = i;
        break;
      }
    }

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
      interp_urr_xs *= smooth;
    }

    // Take the maximum over temperature.
    max_urr_xs = std::max(max_urr_xs, interp_urr_xs);
  }

  return max_urr_xs;
}

int
Majorant::get_i_grid(double energy, const Nuclide::EnergyGrid & grid)
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

double
Majorant::interpolate_lin_1D(double E_0, double E_1, double xs_0, double xs_1, double E)
{
  double f = (E - E_0) / (E_1 - E_0);
  return (1.0 - f) * xs_0 + f * xs_1;
}

double interpolate_log_1D(double E_0, double E_1, double xs_0, double xs_1, double E)
{
  double f = std::log(E / E_0) / std::log(E_1 / E_0);
  return std::exp((1.0 - f) * std::log(xs_0) + f * std::log(xs_1));
}

void
Majorant::vector_union_1D(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & result)
{
  int i = 0;
  int j = 0;

  while (i < a.size() && j < b.size()) {
    if (a[i] < b[j]) {
      if (result.empty() || result.back() != a[i]) {
        result.push_back(a[i]);
      }
      i++;
    } else if (b[j] < a[i]) {
      if (result.empty() || result.back() != b[j]) {
        result.push_back(b[j]);
      }
      j++;
    } else {
      if (result.empty() || result.back() != a[i]) {
        result.push_back(a[i]);
      }
      i++;
      j++;
    }
  }

  while (i < a.size()) {
    if (result.empty() || result.back() != a[i]) {
      result.push_back(a[i]);
    }
    i++;
  }

  while (j < b.size()) {
    if (result.empty() || result.back() != b[j]) {
      result.push_back(b[j]);
    }
    j++;
  }
}
} // namespace openmc
