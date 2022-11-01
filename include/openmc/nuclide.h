//! \file nuclide.h
//! \brief Nuclide type and other associated types/data

#ifndef OPENMC_NUCLIDE_H
#define OPENMC_NUCLIDE_H

#include <array>
#include <memory> // for unique_ptr
#include <unordered_map>
#include <utility> // for pair
#include <vector>

#include <gsl/gsl>
#include <hdf5.h>

#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/particle.h"
#include "openmc/reaction.h"
#include "openmc/reaction_product.h"
#include "openmc/urr.h"
#include "openmc/vector.h"
#include "openmc/wmp.h"

namespace openmc {

//==============================================================================
// Data for a nuclide
//==============================================================================

class Nuclide {
public:
  // Types, aliases
  using EmissionMode = ReactionProduct::EmissionMode;
  struct EnergyGrid {
    std::vector<int> grid_index;
    std::vector<double> energy;
  };

  // Constructors/destructors
  Nuclide(hid_t group, const std::vector<double>& temperature);
  ~Nuclide();

  //! Initialize logarithmic grid for energy searches
  void init_grid();

  #pragma omp declare target
  //! \brief Calculate a microscopic cross section
  //
  //! \param[in] i_log_union Index on the logarithmic unionized grid
  //! \param[in] p particle Particle state used to determine cross section
  //! \param[out] micro Pointer to a struct to return main cross section data in. May be left as nullptr.
  //! \param[out] reaction Array for returning depletion reaction data to. May be left as nullptr.
  //! \param[out] cache Pointer to a full micro XS cache entry for the nuclide. May be left as nullptr
  void calculate_xs(int i_log_union, Particle& p, MicroXS* micro, double* reaction, NuclideMicroXS* cache);

  // Methods
  double nu(double E, EmissionMode mode, int group=0) const;
  double calculate_elastic_xs(int i_temp, int i_grid, double interp_factor ) const;

  //! Determines the microscopic 0K elastic cross section at a trial relative
  //! energy used in resonance scattering
  double elastic_xs_0K(double E) const;
  #pragma omp end declare target

  //! \brief Calculate reaction rate based on group-wise flux distribution
  //
  //! \param[in] MT ENDF MT value for desired reaction
  //! \param[in] temperature Temperature in [K]
  //! \param[in] energy Energy group boundaries in [eV]
  //! \param[in] flux Flux in each energy group (not normalized per eV)
  //! \return Reaction rate
  double collapse_rate(int MT, double temperature, gsl::span<const double> energy,
    gsl::span<const double> flux) const;

  void flatten_xs_data();

  void flatten_wmp_data();

  void copy_to_device();
  void release_from_device();

  // Data members
  std::string name_; //!< Name of nuclide, e.g. "U235"
  int Z_; //!< Atomic number
  int A_; //!< Mass number
  int metastable_; //!< Metastable state
  double awr_; //!< Atomic weight ratio
  gsl::index index_; //!< Index in the nuclides array

  // Temperature dependent cross section data
  vector<double> kTs_; //!< temperatures in eV (k*T)
  std::vector<EnergyGrid> grid_; //!< Energy grid at each temperature
  std::vector<xt::xtensor<double, 2>> xs_; //!< Cross sections at each temperature

  // Flattened 1D temperature dependent cross section data
  int total_energy_gridpoints_;
  int total_index_gridpoints_;
  int* flat_temp_offsets_ {nullptr};
  int* flat_grid_index_;
  double* flat_grid_energy_;
  double* flat_xs_;

  // Multipole data
  std::unique_ptr<WindowedMultipole> multipole_;
  WindowedMultipole* device_multipole_;
  const WindowedMultipole* multipole() const { return device_multipole_; }

  // Fission data
  bool fissionable_ {false}; //!< Whether nuclide is fissionable
  bool has_partial_fission_ {false}; //!< has partial fission reactions?
  std::vector<ReactionFlatContainer*> fission_rx_; //!< Fission reactions
  int n_precursor_ {0}; //!< Number of delayed neutron precursors
  std::unique_ptr<Function1DFlatContainer> total_nu_; //!< Total neutron yield
  std::unique_ptr<Function1DFlatContainer> fission_q_prompt_; //!< Prompt fission energy release
  std::unique_ptr<Function1DFlatContainer> fission_q_recov_; //!< Recoverable fission energy release
  std::unique_ptr<Function1DFlatContainer> prompt_photons_; //!< Prompt photon energy release
  std::unique_ptr<Function1DFlatContainer> delayed_photons_; //!< Delayed photon energy release
  std::unique_ptr<Function1DFlatContainer> fragments_; //!< Fission fragment energy release
  std::unique_ptr<Function1DFlatContainer> betas_; //!< Delayed beta energy release
  Function1DFlatContainer* device_total_nu_;

  // Resonance scattering information
  bool resonant_ {false};
  vector<double> energy_0K_;
  vector<double> elastic_0K_;
  vector<double> xs_cdf_;

  // Unresolved resonance range information
  bool urr_present_ {false};
  int urr_inelastic_ {C_NONE};
  vector<UrrData> urr_data_;

  vector<ReactionFlatContainer> reactions_; //!< Reactions
  std::array<size_t, 902> reaction_index_; //!< Index of each reaction
  vector<int> index_inelastic_scatter_;

  ReactionFlatContainer** device_fission_rx_;
private:
  void create_derived(const Function1DFlatContainer* prompt_photons, const Function1DFlatContainer* delayed_photons);

  //! Determine temperature index and interpolation factor
  //
  //! \param[in] T Temperature in [K]
  //! \return Temperature index and interpolation factor
  std::pair<gsl::index, double> find_temperature(double T) const;

  #pragma omp declare target
  static int XS_TOTAL;
  static int XS_ABSORPTION;
  static int XS_FISSION;
  static int XS_NU_FISSION;
  static int XS_PHOTON_PROD;
  #pragma omp end declare target
};

//==============================================================================
// Non-member functions
//==============================================================================

//! Checks for the right version of nuclear data within HDF5 files
void check_data_version(hid_t file_id);

bool multipole_in_range(const Nuclide& nuc, double E);

//==============================================================================
// Global variables
//==============================================================================

namespace data {

// Minimum/maximum transport energy for each particle type. Order corresponds to
// that of the ParticleType enum
#pragma omp declare target
extern std::array<double, 2> energy_min;
extern std::array<double, 2> energy_max;
#pragma omp end declare target

//! Minimum temperature in [K] that nuclide data is available at
extern double temperature_min;

//! Maximum temperature in [K] that nuclide data is available at
extern double temperature_max;

extern std::unordered_map<std::string, int> nuclide_map;
#pragma omp declare target
extern Nuclide* nuclides;
extern size_t nuclides_size;
#pragma omp end declare target
extern size_t nuclides_capacity;

} // namespace data

//==============================================================================
// Non-member functions
//==============================================================================

void nuclides_clear();

} // namespace openmc

#endif // OPENMC_NUCLIDE_H
