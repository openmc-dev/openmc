//! \file majorant.h
//! \brief Majorant cross section type

#ifndef OPENMC_MAJORANT_H
#define OPENMC_MAJORANT_H

#include <vector>

#include "openmc/settings.h"
#include "openmc/nuclide.h"
#include "openmc/material.h"

namespace openmc {

class Majorant;

namespace data {
  extern std::unique_ptr<Majorant> n_majorant;
  extern std::string majorant_file;
}

class Majorant {

public:
  Majorant();
  Majorant(const std::vector<double>& energy, const std::vector<double>& xs);
  Majorant(const std::string & majorant_file);

  void write_ascii(const std::string& filename) const;

  //! \brief Calculate the microscopic cross section at a given energy
  double calculate_xs(double energy) const;

  // data members
  std::vector<int> nuclides; // index of nuclides applied
  std::vector<double> xs_; // cross section values
  Nuclide::EnergyGrid grid_;
  constexpr static double safety_factor {1.01};

private:
  //! \brief Unionize the smooth and URR cross section grids for all nuclides in the problem.
  void compute_unionized_grid();

  //! \brief Compute the majorant cross section.
  void setup_majorant();

  //! \brief Compute a per-material majorant cross section.
  void fill_material_maj_xs(const Material & mat, const std::vector<double> & to_grid, std::vector<double> & mat_maj);

  //! \brief Compute the maximum smooth cross section for a given energy point.
  double calculate_max_smooth_xs(double energy, const Nuclide & nuc) const;

  //! \brief Compute the maximum URR cross section for a given energy point.
  double calculate_max_urr_xs(double energy, const Nuclide & nuc, double smooth) const;

  //! \brief Compute the maximum correction factor for the S(a,b) total cross section.
  double calculate_max_sab_tot_xs(double energy, int i_sab, double sab_frac, const Nuclide & nuc) const;

  //! \brief Get the grid index for energy interpolation.
  int get_i_grid(double energy, const Nuclide::EnergyGrid & grid) const;

  //! \brief Helper function to perform linear-linear interpolation.
  double interpolate_lin_1D(double x_0, double x_1, double y_0, double y_1, double x) const;
}; // class Majorant

  void create_majorant();
}

#endif // OPENMC_MAJORANT_H
