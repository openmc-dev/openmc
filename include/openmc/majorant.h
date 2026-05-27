//! \file majorant.h
//! \brief Majorant cross section type

#ifndef OPENMC_MAJORANT_H
#define OPENMC_MAJORANT_H

#include <vector>

#include "openmc/settings.h"
#include "openmc/nuclide.h"
#include "openmc/material.h"

namespace openmc {

class NeutronMajorant;

//==============================================================================
// Global variables
//==============================================================================

namespace data {
  extern std::unique_ptr<NeutronMajorant> n_majorant;
  extern std::string majorant_file;
} // namespace data

//==============================================================================

class Majorant {
public:
  //----------------------------------------------------------------------------
  // Constructors

  Majorant() = default;
  Majorant(const std::string & majorant_file, double min_E_transport, double max_E_transport);

  //----------------------------------------------------------------------------
  // Methods

  //! \brief Initialize the majorant cross section.
  virtual void init() = 0;

  //! \brief Write the majorant cross section to a CSV file for visualization
  // TODO: remove this when done prototyping.
  //
  //! \param[in] filename The path/name for the majorant file
  void write_ascii(const std::string& filename) const;

  //----------------------------------------------------------------------------
  // Data members

  constexpr static double safety_factor {1.01}; //!< A dilation factor to ensure floating-point round
                                                //!< off and inexact majorant URR and S(a,b) cross sections
                                                //!< don't bias results

protected:
  //----------------------------------------------------------------------------
  // Protected Methods

  //! \brief Helper function to perform linear interpolation.
  //
  //! \param[in] x_0 The first x coordinate
  //! \param[in] x_1 The second x coordinate
  //! \param[in] y_0 The y coordinate associated with x_0
  //! \param[in] y_1 The y coordinate associated with x_1
  //! \param[in] x A x point between x_0 and x_1 to find a y value at
  double interpolate_lin_1D(double x_0, double x_1, double y_0, double y_1, double x) const;

  //! \brief Helper function to perform log interpolation.
  double interpolate_log_1D(double x_0, double x_1, double y_0, double y_1, double x) const;

  //----------------------------------------------------------------------------
  // Protected data members

  Nuclide::EnergyGrid grid_; //!< The unionized energy grid
  std::vector<double> xs_;   //!< Macroscopic majorant cross sections at each grid point in grid_
}; // class Majorant

//==============================================================================

class NeutronMajorant : public Majorant {

public:
  //----------------------------------------------------------------------------
  // Constructors

  NeutronMajorant() = default;
  NeutronMajorant(const std::string & majorant_file);

  //----------------------------------------------------------------------------
  // Public Methods

  virtual void init() override final;

  //! \brief Calculate the macroscopic majorant cross section in units of [cm^-1]
  //
  //! \param[in] energy The energy to compute the cross section at in [eV]
  double calculate_neutron_xs(double energy) const;

private:
  //----------------------------------------------------------------------------
  // Private Methods

  //! \brief A function to unionize the smooth and URR cross section grids for all nuclides in the problem.
  void compute_unionized_grid();

  //! \brief Compute the majorant cross section.
  void setup_majorant();

  //! \brief Compute a per-material macroscopic majorant cross section in units of [barn]
  //
  //! \param[in] mat The material to compute the majorant cross section of
  //! \param[in] to_grid The grid points to evaluate the majorant at in [eV]
  //! \param[out] mat_maj The array to write the macroscopic majorant to. The resulting cross section has units of [cm^-1]
  void fill_material_maj_xs(const Material & mat, const std::vector<double> & to_grid, std::vector<double> & mat_maj);

  //! \brief Compute the maximum smooth microscopic total cross section in units of [barn].
  //
  //! \param[in] energy The energy to evaluate the cross section at in [eV]
  //! \param[in] nuc The nuclide to compute the microscopic total cross section of
  double calculate_max_smooth_xs(double energy, const Nuclide & nuc) const;

  //! \brief Compute the maximum microscopic total URR cross section in units of [barn].
  //
  //! \param[in] energy The energy to evaluate the cross section at in [eV]
  //! \param[in] nuc The nuclide to compute the microscopic total cross section of
  //! \param[in] smooth_xs The smooth total cross section in units of [eV] to use (if needed by the ptable)
  double calculate_max_urr_xs(double energy, const Nuclide & nuc, double smooth_xs) const;

  //! \brief Compute the maximum microscopic S(a,b) total cross section in units of [barn]
  //
  //! \param[in] energy The energy to evaluate the cross section at in [eV]
  //! \param[in] i_sab The index into the thermal scattering table array for this nuclide
  //! \param[in] sab_frac The fraction of the bound cross section to use vs the free gas cross section
  //! \param[in] nuc The nuclide to compute the microscopic total cross section of
  double calculate_max_sab_tot_xs(double energy, int i_sab, double sab_frac, const Nuclide & nuc) const;

  //! \brief Get the grid index for energy interpolation
  //
  //! \param[in] energy The energy to evaluate the cross section at in [eV]
  //! \param[in] grid The energy grid to search for an energy grid index
  int get_i_grid(double energy, const Nuclide::EnergyGrid & grid) const;
}; // class NeutronMajorant

//==============================================================================

//==============================================================================
// Non-member functions
//==============================================================================

//! \brief A function to create the global majorant cross section(s).
void create_majorants();
}

#endif // OPENMC_MAJORANT_H
