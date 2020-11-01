//! \file nuclide.h
//! \brief Nuclide type and other associated types/data

#ifndef OPENMC_NUCLIDE_H
#define OPENMC_NUCLIDE_H

#include <unordered_map>
#include <utility> // for pair

#include <gsl/gsl>
#include <hdf5.h>

#include "openmc/array.h"
#include "openmc/constants.h"
#include "openmc/endf.h"
#include "openmc/memory.h"
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
    // See vector.h if you're trying to understand the CUDA version. Otherwise,
    // the things called replicated_vector are just normal std::vectors.
    replicated_vector<int> grid_index;
    replicated_vector<xsfloat> energy;
  };

  struct CrossSectionSet {
    xsfloat total;
    xsfloat absorption;
    xsfloat fission;
    xsfloat nu_fission;
    xsfloat photon_production;

    // Construction gives zero XS
    CrossSectionSet()
      : total(0), absorption(0), fission(0), nu_fission(0), photon_production(0)
    {}
  };

  // Constructors/destructors
  Nuclide(hid_t group, const vector<xsfloat>& temperature);
  ~Nuclide();

  //! Initialize logarithmic grid for energy searches
  void init_grid();

  void calculate_xs(int i_sab, int i_log_union, double sab_frac, Particle& p);

  void calculate_sab_xs(int i_sab, double sab_frac, Particle& p);

  // Methods
  xsfloat nu(xsfloat E, EmissionMode mode, int group=0) const;
  void calculate_elastic_xs(Particle& p) const;

  //! Determines the microscopic 0K elastic cross section at a trial relative
  //! energy used in resonance scattering
  xsfloat elastic_xs_0K(xsfloat E) const;

  //! \brief Determines cross sections in the unresolved resonance range
  //! from probability tables.
  void calculate_urr_xs(int i_temp, Particle& p) const;

  //! \brief Calculate reaction rate based on group-wise flux distribution
  //
  //! \param[in] MT ENDF MT value for desired reaction
  //! \param[in] temperature Temperature in [K]
  //! \param[in] energy Energy group boundaries in [eV]
  //! \param[in] flux Flux in each energy group (not normalized per eV)
  //! \return Reaction rate
  xsfloat collapse_rate(int MT, xsfloat temperature, gsl::span<const xsfloat> energy,
    gsl::span<const double> flux) const;

  // Data members
  std::string name_; //!< Name of nuclide, e.g. "U235"
  int Z_; //!< Atomic number
  int A_; //!< Mass number
  int metastable_; //!< Metastable state
  double awr_; //!< Atomic weight ratio
  gsl::index index_; //!< Index in the nuclides array

  // Temperature dependent cross section data
  vector<xsfloat> kTs_;                //!< temperatures in eV (k*T)
  replicated_vector<EnergyGrid> grid_; //!< Energy grid at each temperature
  replicated_vector<replicated_vector<CrossSectionSet>>
    xs_; //!< Cross sections at each temperature

  // Multipole data
  unique_ptr<WindowedMultipole> multipole_;

  // Fission data
  bool fissionable_ {false}; //!< Whether nuclide is fissionable
  bool has_partial_fission_ {false}; //!< has partial fission reactions?
  vector<Reaction*> fission_rx_;     //!< Fission reactions
  int n_precursor_ {0}; //!< Number of delayed neutron precursors
  unique_ptr<Function1D> total_nu_;         //!< Total neutron yield
  unique_ptr<Function1D> fission_q_prompt_; //!< Prompt fission energy release
  unique_ptr<Function1D>
    fission_q_recov_; //!< Recoverable fission energy release
  unique_ptr<Function1D> prompt_photons_;  //!< Prompt photon energy release
  unique_ptr<Function1D> delayed_photons_; //!< Delayed photon energy release
  unique_ptr<Function1D> fragments_;       //!< Fission fragment energy release
  unique_ptr<Function1D> betas_;           //!< Delayed beta energy release

  // Resonance scattering information
  bool resonant_ {false};
  vector<xsfloat> energy_0K_;
  vector<xsfloat> elastic_0K_;
  vector<xsfloat> xs_cdf_;

  // Unresolved resonance range information
  bool urr_present_ {false};
  int urr_inelastic_ {C_NONE};
  vector<UrrData> urr_data_;

  vector<unique_ptr<Reaction>> reactions_; //!< Reactions
  array<size_t, 902> reaction_index_;      //!< Index of each reaction
  vector<int> index_inelastic_scatter_;

private:
  void create_derived(const Function1D* prompt_photons, const Function1D* delayed_photons);

  //! Determine temperature index and interpolation factor
  //
  //! \param[in] T Temperature in [K]
  //! \return Temperature index and interpolation factor
  std::pair<gsl::index, double> find_temperature(double T) const;
};

//==============================================================================
// Non-member functions
//==============================================================================

//! Checks for the right version of nuclear data within HDF5 files
void check_data_version(hid_t file_id);

bool multipole_in_range(const Nuclide& nuc, xsfloat E);

//==============================================================================
// Global variables
//==============================================================================

namespace data {

// Minimum/maximum transport energy for each particle type. Order corresponds to
// that of the ParticleType enum
extern array<xsfloat, 2> energy_min;
extern array<xsfloat, 2> energy_max;

//! Minimum temperature in [K] that nuclide data is available at
extern xsfloat temperature_min;

//! Maximum temperature in [K] that nuclide data is available at
extern xsfloat temperature_max;

extern std::unordered_map<std::string, int> nuclide_map;
extern vector<unique_ptr<Nuclide>> nuclides;

} // namespace data

//==============================================================================
// Non-member functions
//==============================================================================

void nuclides_clear();

} // namespace openmc

#endif // OPENMC_NUCLIDE_H
