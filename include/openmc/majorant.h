//! \file majorant.h
//! Majorant cross section type

#ifndef OPENMC_MAJORANT_H
#define OPENMC_MAJORANT_H

#include <vector>

#include "openmc/nuclide.h"

namespace openmc {

class NeutronMajorant;
class PhotonMajorant;

//==============================================================================
// Global variables
//==============================================================================

namespace data {
extern std::unique_ptr<NeutronMajorant> n_majorant;
extern std::unique_ptr<PhotonMajorant> p_majorant;
} // namespace data

//==============================================================================

class Majorant {
public:
  //----------------------------------------------------------------------------
  // Constructor

  Majorant(int i_universe);

  //----------------------------------------------------------------------------
  // Methods

  //! A function to unionize particle energy grids.
  virtual void compute_unionized_grid() = 0;

  //! A function to populate the majorant cross section.
  void compute_majorant();

protected:
  //----------------------------------------------------------------------------
  // Protected Methods

  //! Compute a per-material macroscopic majorant cross section in units
  //! of [cm^-1]
  //!
  //! \param[in] mat The material to compute the majorant cross section of
  //! \param[in] to_grid The grid points to evaluate the majorant at in [eV]
  //! \param[out] mat_maj The array to write the macroscopic majorant to.
  //!   The resulting cross section has units of [cm^-1]
  virtual void fill_material_maj_xs(int i_material, double max_density_mult,
    const std::vector<double>& to_grid, std::vector<double>& mat_maj) const = 0;

  //! Post-processes the energy grid by calling std::sort(), std::unique().
  //! This also removes all energy values below the transport minimum and
  //! above the transport maximum for a given particle type.
  //!
  //! \param[in] particle_type The particle type transport index to use when
  //!   fetching transport minimum / maximum energies
  //! \param[out] grid The energy grid to post-process. This is performed
  //!   in-place.
  void post_process_grid(int particle_type, Nuclide::EnergyGrid& grid) const;

  //----------------------------------------------------------------------------
  // Protected data members

  //! Index into the universe array for the universe which this majorant uses
  //! to fetch material properties.
  int maj_universe_ = C_NONE;
  //!< A vector of materials contained in maj_universe_.
  std::vector<int> contained_materials_;
  //! A map of each material index and the corresponding maximum density
  //! multiplier applied to that material by a cell.
  std::unordered_map<int, double> max_density_mult_;
  //! The unionized energy grid.
  Nuclide::EnergyGrid grid_;
  //! Macroscopic majorant cross sections at each energy point in grid_.
  std::vector<double> xs_;
}; // class Majorant

//==============================================================================

class NeutronMajorant : public Majorant {

public:
  //----------------------------------------------------------------------------
  // Constructors

  NeutronMajorant(int i_universe);

  //----------------------------------------------------------------------------
  // Public Methods

  virtual void compute_unionized_grid() override final;

  //! Calculate the macroscopic majorant cross section.
  //!
  //! \param[in] energy The energy to compute the cross section at in [eV]
  //! \return The neutron majorant cross section at the given energy in [cm^-1]
  double calculate_neutron_xs(double energy) const;

  //----------------------------------------------------------------------------
  // Data members

  //! A dilation factor to ensure floating-point round off and inexact majorant
  //! URR and S(a,b) cross sections don't bias results.
  constexpr static double safety_factor_ {1.01};

protected:
  //----------------------------------------------------------------------------
  // Protected Methods

  //! Compute a per-material macroscopic majorant cross section.
  //!
  //! \param[in] i_material Index into the materials array.
  //! \param[in] to_grid The grid points to evaluate the majorant at in [eV].
  //! \param[out] mat_maj The array to write the macroscopic majorant to.
  //!   The resulting cross section has units of [cm^-1].
  virtual void fill_material_maj_xs(int i_material, double max_density_mult,
    const std::vector<double>& to_grid,
    std::vector<double>& mat_maj) const override;

private:
  //----------------------------------------------------------------------------
  // Private Methods

  //! Compute the maximum smooth microscopic total cross section.
  //!
  //! \param[in] i_nuclide Index into the nuclides array.
  //! \param[in] energy The energy to evaluate the cross section at in [eV].
  //! \return The maximum smooth cross section in [barn].
  double calculate_max_smooth_xs(int i_nuclide, double energy) const;

  //! Compute the maximum microscopic total URR cross section.
  //!
  //! \param[in] i_nuclide Index into the nuclides array.
  //! \param[in] energy The energy to evaluate the cross section at in [eV].
  //! \param[in] smooth_xs The smooth total cross section in units of [eV] to
  //!   use (if needed by the ptable).
  //! \return The maximum URR total cross section in [barn].
  double calculate_max_urr_xs(
    int i_nuclide, double energy, double smooth_xs) const;

  //! Compute the maximum microscopic S(a,b) total cross section.
  //!
  //! \param[in] energy The energy to evaluate the cross section at in [eV].
  //! \param[in] i_sab The index into the thermal scattering table array for
  //!   this nuclide.
  //! \param[in] sab_frac The fraction of the bound cross section to use vs the
  //!   free gas cross section.
  //! \param[in] nuc The nuclide to compute the microscopic total cross section
  //!   of.
  //! \return The maximum S(a,b) total cross section in [barn].
  double calculate_max_sab_tot_xs(
    int i_nuclide, int i_sab, double sab_frac, double energy) const;

  //! Get the grid index for energy interpolation.
  //!
  //! \param[in] energy The energy to evaluate the cross section at in [eV].
  //! \param[in] grid The energy grid to search for an energy grid index.
  //! \return The grid index.
  int get_i_grid(double energy, const Nuclide::EnergyGrid& grid) const;

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

  //----------------------------------------------------------------------------
  // Public Methods

  virtual void compute_unionized_grid() override final;

  //! Calculate the macroscopic majorant cross section.
  //!
  //! \param[in] energy The energy to compute the cross section at in [eV]
  //! \return The photon majorant cross section at the given energy in [cm^-1]
  double calculate_photon_xs(double energy) const;

  //----------------------------------------------------------------------------
  // Data members

  //! A dilation factor to catch interpolation error.
  constexpr static double safety_factor_ {1.01};

protected:
  //----------------------------------------------------------------------------
  // Protected Methods

  //! Compute a per-material macroscopic majorant cross section.
  //!
  //! \param[in] i_material Index into the materials array
  //! \param[in] to_grid The grid points to evaluate the majorant at in [eV]
  //! \param[out] mat_maj The array to write the macroscopic majorant to. The
  //!   resulting cross section has units of [cm^-1]
  virtual void fill_material_maj_xs(int i_material, double max_density_mult,
    const std::vector<double>& to_grid,
    std::vector<double>& mat_maj) const override;

private:
  //----------------------------------------------------------------------------
  // Private Methods

  //! Compute the maximum smooth microscopic total cross section.
  //!
  //! \param[in] i_element Index in the elements array.
  //! \param[in] energy The energy to evaluate the cross section at in [eV].
  //! \return The maximum microscopic total cross section in [barn].
  double calculate_elem_tot_xs(int i_element, double log_energy) const;

  //! Get the grid index for energy interpolation.
  //!
  //! \param[in] energy The energy to evaluate the cross section at in [eV].
  //! \param[in] grid The energy grid to search for an energy grid index.
  //! \return The grid index.
  int get_i_grid(
    double log_energy, const std::vector<double>& energy_grid) const;
  int get_i_grid(
    double log_energy, const tensor::Tensor<double>& energy_grid) const;

  //----------------------------------------------------------------------------
  // Private data members

  static constexpr int i_photon_ = ParticleType::photon().transport_index();
}; // class PhotonMajorant

//==============================================================================
// Static functions
//==============================================================================

//! A function to create majorant cross sections.
void create_majorants();

//! A function to reset majorant cross sections.
void reset_majorants();
} // namespace openmc

#endif // OPENMC_MAJORANT_H
