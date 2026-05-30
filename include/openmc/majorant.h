//! \file majorant.h
//! \brief Majorant cross section type

#ifndef OPENMC_MAJORANT_H
#define OPENMC_MAJORANT_H

#include <vector>

#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/photon.h"

namespace openmc {

class NeutronMajorant;
class PhotonMajorant;

//==============================================================================
// Global variables
//==============================================================================

namespace data {
  extern std::unique_ptr<NeutronMajorant> n_majorant;
  extern std::unique_ptr<PhotonMajorant> p_majorant;
  extern std::string n_majorant_file;
  extern std::string p_majorant_file;
} // namespace data

//==============================================================================

class Majorant {
public:
  //----------------------------------------------------------------------------
  // Constructors

  Majorant(int i_universe);
  Majorant(const std::string & majorant_file, int p_transport_indx);

  //----------------------------------------------------------------------------
  // Methods

  //! \brief A function to unionize particle energy grids.
  virtual void compute_unionized_grid() = 0;

  //! \brief Populate the majorant cross section.
  void compute_majorant();

  //! \brief Write the majorant cross section to a CSV file for visualization
  // TODO: remove this when done prototyping.
  //
  //! \param[in] filename The path/name for the majorant file
  void write_ascii(const std::string& filename) const;

protected:
  //----------------------------------------------------------------------------
  // Protected Methods

  //! \brief Find materials in maj_universe_ by traversing the geometry tree.
  //
  void discover_contained_materials();

  //! \brief Compute a per-material macroscopic majorant cross section in units of [cm^-1]
  //
  //! \param[in] mat The material to compute the majorant cross section of
  //! \param[in] to_grid The grid points to evaluate the majorant at in [eV]
  //! \param[out] mat_maj The array to write the macroscopic majorant to. The resulting cross section has units of [cm^-1]
  virtual void fill_material_maj_xs(const Material & mat, const std::vector<double> & to_grid, std::vector<double> & mat_maj) = 0;

  //! \brief Helper function to perform linear interpolation.
  //
  //! \param[in] x_0 The first x coordinate
  //! \param[in] x_1 The second x coordinate
  //! \param[in] y_0 The y coordinate associated with x_0
  //! \param[in] y_1 The y coordinate associated with x_1
  //! \param[in] x The point between x_0 and x_1 to find a y value at
  inline double interpolate_lin_1D(double x_0, double x_1, double y_0, double y_1, double x) const
  {
    const double f = (x - x_0) / (x_1 - x_0);
    return (1.0 - f) * y_0 + f * y_1;
  }

  //! \brief Helper function to perform log interpolation.
  //
  //! \param[in] x_0 The first x coordinate
  //! \param[in] x_1 The second x coordinate
  //! \param[in] y_0 The y coordinate associated with x_0
  //! \param[in] y_1 The y coordinate associated with x_1
  //! \param[in] x The point between x_0 and x_1 to find a y value at
  inline double interpolate_log_1D(double x_0, double x_1, double y_0, double y_1, double x) const
  {
    const double f = std::log(x / x_0) / std::log(x_1 / x_0);
    return std::exp((1.0 - f) * std::log(y_0) + f * std::log(y_1));
  }

  //----------------------------------------------------------------------------
  // Protected data members

  int maj_universe_ = C_NONE;            //!< Index into the universe array for the universe which this
                                         //   majorant uses to fetch material properties.
  std::vector<int> contained_materials_; //!< A vector of materials contained in maj_universe_
  Nuclide::EnergyGrid grid_;             //!< The unionized energy grid
  std::vector<double> xs_;               //!< Macroscopic majorant cross sections at each grid point in grid_
}; // class Majorant

//==============================================================================

class NeutronMajorant : public Majorant {

public:
  //----------------------------------------------------------------------------
  // Constructors

  NeutronMajorant(int i_universe);
  NeutronMajorant(const std::string & majorant_file);

  //----------------------------------------------------------------------------
  // Public Methods

  virtual void compute_unionized_grid() override final;

  //! \brief Calculate the macroscopic majorant cross section in units of [cm^-1]
  //
  //! \param[in] energy The energy to compute the cross section at in [eV]
  double calculate_neutron_xs(double energy) const;

  //----------------------------------------------------------------------------
  // Data members

  constexpr static double safety_factor_ {1.01}; //!< A dilation factor to ensure floating-point round
                                                 //!< off and inexact majorant URR and S(a,b) cross sections
                                                 //!< don't bias results

protected:
  //----------------------------------------------------------------------------
  // Protected Methods

  //! \brief Compute a per-material macroscopic majorant cross section in units of [cm^-1]
  //
  //! \param[in] mat The material to compute the majorant cross section of
  //! \param[in] to_grid The grid points to evaluate the majorant at in [eV]
  //! \param[out] mat_maj The array to write the macroscopic majorant to. The resulting cross section has units of [cm^-1]
  virtual void fill_material_maj_xs(const Material & mat, const std::vector<double> & to_grid, std::vector<double> & mat_maj) override;

private:
  //----------------------------------------------------------------------------
  // Private Methods

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

  //----------------------------------------------------------------------------
  // Private data members

  static constexpr int i_neutron_ = ParticleType::neutron().transport_index();
}; // class NeutronMajorant

//==============================================================================

class PhotonMajorant : public Majorant {
public:
  //----------------------------------------------------------------------------
  // Constructors

  PhotonMajorant(int i_universe);
  PhotonMajorant(const std::string & majorant_file);

  //----------------------------------------------------------------------------
  // Public Methods

  virtual void compute_unionized_grid() override final;

  //! \brief Calculate the macroscopic majorant cross section in units of [cm^-1]
  //
  //! \param[in] energy The energy to compute the cross section at in [eV]
  double calculate_photon_xs(double energy) const;

  //----------------------------------------------------------------------------
  // Data members

  constexpr static double safety_factor_ {1.01}; //!< A dilation factor to catch interpolation error

protected:
  //----------------------------------------------------------------------------
  // Protected Methods

  //! \brief Compute a per-material macroscopic majorant cross section in units of [cm^-1]
  //
  //! \param[in] mat The material to compute the majorant cross section of
  //! \param[in] to_grid The grid points to evaluate the majorant at in [eV]
  //! \param[out] mat_maj The array to write the macroscopic majorant to. The resulting cross section has units of [cm^-1]
  virtual void fill_material_maj_xs(const Material & mat, const std::vector<double> & to_grid, std::vector<double> & mat_maj) override;

private:
  //----------------------------------------------------------------------------
  // Private Methods

  //! \brief Compute the maximum smooth microscopic total cross section in units of [barn].
  //
  //! \param[in] energy The energy to evaluate the cross section at in [eV]
  //! \param[in] nuc The nuclide to compute the microscopic total cross section of
  double calculate_elem_tot_xs(double log_energy, const PhotonInteraction & elem) const;

  //! \brief Get the grid index for energy interpolation
  //
  //! \param[in] energy The energy to evaluate the cross section at in [eV]
  //! \param[in] grid The energy grid to search for an energy grid index
  int get_i_grid(double log_energy, const std::vector<double> & energy_grid) const;
  int get_i_grid(double log_energy, const tensor::Tensor<double> & energy_grid) const;

  //----------------------------------------------------------------------------
  // Private data members

  static constexpr int i_photon_ = ParticleType::photon().transport_index();
}; // class PhotonMajorant

//==============================================================================
// Static functions
//==============================================================================

//! \brief A function to create majorant cross sections.
void create_majorants();

//! \brief A function to reset majorant cross sectiosn.
void reset_majorants();
} // namespace openmc

#endif // OPENMC_MAJORANT_H
