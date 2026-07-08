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

  //! Virtual function for the minimum transport energy in the majorant.
  //!
  //! \return The minimum transport energy associated with the majorant
  virtual double min_transport_energy() const = 0;

  //! Virtual function for the maximum transport energy in the majorant.
  //!
  //! \return The maximum transport energy associated with the majorant
  virtual double max_transport_energy() const = 0;

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
  //! \param[out] grid The energy grid to post-process. This is performed
  //!   in-place.
  void post_process_grid(Nuclide::EnergyGrid& grid) const;

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

  //! Minimum neutron transport energy.
  //
  //! \return The minimum transport energy associated with the majorant [eV]
  virtual double min_transport_energy() const override { return data::energy_min[i_neutron_]; }

  //! Maximum neutron transport energy.
  //
  //! \return The maximum transport energy associated with the majorant [eV]
  virtual double max_transport_energy() const override { return data::energy_max[i_neutron_]; }

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

  //! Minimum neutron transport energy.
  //
  //! \return The minimum transport energy associated with the majorant [eV]
  virtual double min_transport_energy() const override { return std::log(data::energy_min[i_photon_]); }

  //! Maximum neutron transport energy.
  //
  //! \return The maximum transport energy associated with the majorant [eV]
  virtual double max_transport_energy() const override { return std::log(data::energy_max[i_photon_]); }

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

  //! Get the grid index for energy interpolation. This is templated due to the
  //! use of tensor::Tensor<double> in PhotonInteraction.
  //!
  //! \param[in] energy The energy to evaluate the cross section at in [eV].
  //! \param[in] grid The energy grid to search for an energy grid index.
  //! \return The grid index.
  template <typename T>
  int get_i_grid(double log_energy, const T& energy_grid) const {
    int n_grid = energy_grid.size();
    int i_grid;
    if (log_energy <= energy_grid[0]) {
      i_grid = 0;
    } else if (log_energy > energy_grid[n_grid - 1]) {
      i_grid = n_grid - 2;
    } else {
      // We use upper_bound_index here because sometimes photons are created with
      // energies that exactly match a grid point
      i_grid =
        upper_bound_index(energy_grid.cbegin(), energy_grid.cend(), log_energy);
    }

    // check for case where two energy points are the same
    if (energy_grid[i_grid] == energy_grid[i_grid + 1])
      ++i_grid;

    return i_grid;
  }

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
