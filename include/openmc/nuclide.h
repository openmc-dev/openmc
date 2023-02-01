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
#include "openmc/thermal.h"
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
  // Template function for performing a microscopic cross section lookup.
  // This function is templated due to emprical results showing high performance sensitivity
  // to the amount of data being returned by the function. We have also experimented with
  // passing in several pointer arguments to provide space for what gets optionally returned
  // to the caller, but the LLVM compiler appears to generate much faster code if the
  // function returns an object directly rather than writing to a reference passed in.
  // The template type is used at the end of the function as part of the return call as
  // the correct type is generated via a constructor that knows what to do with all the possible
  // parameters that might need to be passed back to the caller.
  template <typename T>
  T calculate_xs(int i_log_union, Particle& p, bool need_depletion_rx, double E, double sqrtkT)
  {
    double reaction[DEPLETION_RX_SIZE] = {};

    // ======================================================================
    // CHECK FOR SAB TABLE BEGIN
    // ======================================================================

    int i_sab = C_NONE;
    double sab_frac = 0.0;

    const auto& mat = model::materials[p.material_];
    for (int s = 0; s < mat.thermal_tables_.size(); s++) {
      const auto& sab {mat.thermal_tables(s)};
      if (index_ == sab.index_nuclide) {
        // Get index in sab_tables
        i_sab = sab.index_table;
        sab_frac = sab.fraction;

        // If particle energy is greater than the highest energy for the
        // S(a,b) table, then don't use the S(a,b) table
        if (E > data::device_thermal_scatt[i_sab].energy_max_) i_sab = C_NONE;
      }
    }

    // ======================================================================
    // CHECK FOR SAB TABLE END
    // ======================================================================

    #ifndef NO_MICRO_XS_CACHE
    {
      // Check if a microscopic XS lookup is even required. If all state variables are the same,
      // we don't need to perform a lookup at all.
      if (     E      == cache.last_E
          && sqrtkT == cache.last_sqrtkT
          && i_sab     == cache.index_sab
          && sab_frac  == cache.sab_frac
         )
      {
        // If the cache is still valid, then we can pass back any needed values directly
        // from the cache
        if (micro) {
          micro->total      = cache.total;
          micro->absorption = cache.absorption;
          micro->fission    = cache.fission;
          micro->nu_fission = cache.nu_fission;
        }
        if (need_depletion_rx) {
          for ( int r = 0; r < DEPLETION_RX_SIZE; r++) {
            reaction[r] = cache.reaction[r];
          }
        }
        return;
      }
    }
    #endif

    double total;
    double elastic = CACHE_INVALID;
    double absorption;
    double fission;
    double nu_fission;
    double photon_prod = 0.0;

    bool use_mp = false;
    // Check to see if there is multipole data present at this energy
    if (multipole()) {
      use_mp = (E >= multipole()->E_min_ && E <= multipole()->E_max_);
    }

    int i_temp = -1;
    int i_grid;
    double f;
    // Evaluate multipole or interpolate
    if (use_mp) {
      // ======================================================================
      // MULTIPOLE LOOKUP BEGIN
      // ======================================================================

      // Call multipole kernel
      double sig_s, sig_a, sig_f;
      std::tie(sig_s, sig_a, sig_f) = multipole()->evaluate(E, sqrtkT);

      total = sig_s + sig_a;
      elastic = sig_s;
      absorption = sig_a;
      fission = sig_f;
      nu_fission = fissionable_ ?
        sig_f * this->nu(E, EmissionMode::total) : 0.0;

      if (need_depletion_rx) {
        // Only non-zero reaction is (n,gamma)
        reaction[0] = sig_a - sig_f;
      }

      // Ensure these values are set
      // Note, the only time either is used is in one of 4 places:
      // 1. physics.cpp - scatter - For inelastic scatter.
      // 2. physics.cpp - sample_fission - For partial fissions.
      // 3. tally.F90 - score_general - For tallying on MTxxx reactions.
      // 4. nuclide.cpp - calculate_urr_xs - For unresolved purposes.
      // It is worth noting that none of these occur in the resolved
      // resonance range, so the value here does not matter.  index_temp is
      // set to -1 to force a segfault in case a developer messes up and tries
      // to use it with multipole.
      i_grid = -1;
      f = 0.0;

      // ======================================================================
      // MULTIPOLE LOOKUP END
      // ======================================================================
    } else {
      // ======================================================================
      // POINTWISE LOOKUP BEGIN
      // ======================================================================

      // Find the appropriate temperature index.
      double kT = sqrtkT*sqrtkT;
      switch (settings::temperature_method) {
        case TemperatureMethod::NEAREST:
          {
            double max_diff = INFTY;
            for (int t = 0; t < kTs_.size(); ++t) {
              double diff = std::abs(kTs_[t] - kT);
              if (diff < max_diff) {
                i_temp = t;
                max_diff = diff;
              }
            }
          }
          break;

        case TemperatureMethod::INTERPOLATION:
          // Find temperatures that bound the actual temperature
          for (i_temp = 0; i_temp < kTs_.size() - 1; ++i_temp) {
            if (kTs_[i_temp] <= kT && kT < kTs_[i_temp + 1]) break;
          }

          // Randomly sample between temperature i and i+1
          f = (kT - kTs_[i_temp]) / (kTs_[i_temp + 1] - kTs_[i_temp]);
          if (f > prn(p.current_seed())) ++i_temp;
          break;
      }

      // Offset index grid
      int index_offset = i_temp * (settings::n_log_bins + 1);
      int* grid_index = &flat_grid_index_[index_offset];

      // Offset energy grid
      int energy_offset = flat_temp_offsets_[i_temp];
      double* energy = &flat_grid_energy_[energy_offset];

      // Offset xs
      int xs_offset = flat_temp_offsets_[i_temp] * 5;
      double* xs = &flat_xs_[xs_offset];

      // Determine # of gridpoints for this temperature
      int num_gridpoints;
      if (i_temp < kTs_.size() - 1) {
        num_gridpoints = flat_temp_offsets_[i_temp + 1] - energy_offset;
      } else {
        num_gridpoints = total_energy_gridpoints_ - energy_offset;
      }

      // Determine the energy grid index using a logarithmic mapping to
      // reduce the energy range over which a binary search needs to be
      // performed

      if (E < energy[0]) {
        i_grid = 0;
      } else if (E > energy[num_gridpoints-1]) {
        i_grid = num_gridpoints - 2;
      } else {
        // Determine bounding indices based on which equal log-spaced
        // interval the energy is in
        int i_low  = grid_index[i_log_union];
        int i_high = grid_index[i_log_union + 1] + 1;

        // Perform binary search over reduced range
        // Note the STL-based binary search seems to work on llvm/V100 but not elsewhere
        //i_grid = i_low + lower_bound_index(&energy[i_low], &energy[i_high], E);

        // Iterative linear search (a few percent faster on device due to reduced branching)
        for (; i_low < i_high - 1; i_low++) {
          if (E < energy[i_low + 1])
            break;
        }
        i_grid = i_low;
      }

      // check for rare case where two energy points are the same
      if (energy[i_grid] == energy[i_grid + 1]) ++i_grid;

      // 1D indexing conversion for lower XS value
      int i_grid1D = i_grid * 5;
      int i_next1D = (i_grid + 1) * 5;

      // Execute All Lookups (in CUDA, it would be nice to try __ldg())
      double total_low            = xs[i_grid1D + XS_TOTAL];
      double absorption_low       = xs[i_grid1D + XS_ABSORPTION];
      double fission_low          = xs[i_grid1D + XS_FISSION];
      double nu_fission_low       = xs[i_grid1D + XS_NU_FISSION];
      double photon_prod_low      = xs[i_grid1D + XS_PHOTON_PROD];

      double total_next       = xs[i_next1D + XS_TOTAL];
      double absorption_next  = xs[i_next1D + XS_ABSORPTION];
      double fission_next     = xs[i_next1D + XS_FISSION];
      double nu_fission_next  = xs[i_next1D + XS_NU_FISSION];
      double photon_prod_next = xs[i_next1D + XS_PHOTON_PROD];

      // Calculate interpolation factor and complement
      f = (E - energy[i_grid]) /
        (energy[i_grid + 1]- energy[i_grid]);
      double f_comp = 1.0 - f;

      // Perform Interpolation
      total       = f_comp * total_low       + f * total_next;
      absorption  = f_comp * absorption_low  + f * absorption_next;
      fission     = f_comp * fission_low     + f * fission_next;
      nu_fission  = f_comp * nu_fission_low  + f * nu_fission_next;
      photon_prod = f_comp * photon_prod_low + f * photon_prod_next;

      // Depletion-related reactions
      if (need_depletion_rx) {
        for (int j = 0; j < DEPLETION_RX.size(); ++j) {
          // If reaction is present and energy is greater than threshold, set the
          // reaction xs appropriately
          int i_rx = reaction_index_[DEPLETION_RX[j]];
          if (i_rx >= 0) {
            const auto& rx = reactions_[i_rx].obj();

            // Physics says that (n,gamma) is not a threshold reaction, so we don't
            // need to specifically check its threshold index
            if (j == 0) {
              reaction[0] = rx.xs(i_temp, i_grid, f);
              continue;
            }

            int threshold = rx.xs_threshold(i_temp);
            if (i_grid >= threshold) {
              reaction[j] = rx.xs(i_temp, i_grid, f);
            } else if (j >= 3) {
              // One can show that the the threshold for (n,(x+1)n) is always
              // higher than the threshold for (n,xn). Thus, if we are below
              // the threshold for, e.g., (n,2n), there is no reason to check
              // the threshold for (n,3n) and (n,4n).
              break;
            }
          }
        }
      } // end depletion RX conditional

      // ======================================================================
      // POINTWISE LOOKUP END
      // ======================================================================
    }

    int index_sab = C_NONE;
    double thermal = 0.0;
    double thermal_elastic = 0.0;
    int index_temp_sab;

    // ======================================================================
    // SAB BEGIN
    // ======================================================================

    // If there is S(a,b) data for this nuclide, we need to set the sab_scatter
    // and sab_elastic cross sections and correct the total and elastic cross
    // sections.
    if (i_sab >= 0)
    {
      // Set flag that S(a,b) treatment should be used for scattering
      index_sab = i_sab;

      // Calculate the S(a,b) cross section
      int sab_i_temp;
      double sab_elastic;
      double sab_inelastic;
      data::device_thermal_scatt[i_sab].calculate_xs(E, sqrtkT, &sab_i_temp, &sab_elastic, &sab_inelastic, p.current_seed());

      // Store the S(a,b) cross sections.
      thermal = sab_frac * (sab_elastic + sab_inelastic);
      thermal_elastic = sab_frac * sab_elastic;

      // Calculate free atom elastic cross section
      if (elastic == CACHE_INVALID) {
        elastic = this->calculate_elastic_xs(i_temp, i_grid, f);
      }

      // Correct total and elastic cross sections
      total = total + thermal - sab_frac * elastic;
      elastic = thermal + (1.0 - sab_frac) * elastic;

      // Save temperature index and thermal fraction
      index_temp_sab = sab_i_temp;
    } else {
      sab_frac = 0.0;
    }
    // ======================================================================
    // SAB END
    // ======================================================================

    bool use_ptable = false;

    // ======================================================================
    // URR BEGIN
    // ======================================================================

    // If the particle is in the unresolved resonance range and there are
    // probability tables, we need to determine cross sections from the table
    if (settings::urr_ptables_on && urr_present_ && !use_mp) {
      int n = urr_data_[i_temp].n_energy_;
      if ((E > urr_data_[i_temp].device_energy_[0]) &&
          (E < urr_data_[i_temp].device_energy_[n-1]))
      {
        use_ptable = true;

        // Create a shorthand for the URR data
        const auto& urr = urr_data_[i_temp];

        // Determine the energy table
        int i_energy = 0;
        while (E >= urr.device_energy_[i_energy + 1]) {++i_energy;};

        // Sample the probability table using the cumulative distribution

        // Random nmbers for the xs calculation are sampled from a separate stream.
        // This guarantees the randomness and, at the same time, makes sure we
        // reuse random numbers for the same nuclide at different temperatures,
        // therefore preserving correlation of temperature in probability tables.
        p.stream_ = STREAM_URR_PTABLE;
        //TODO: to maintain the same random number stream as the Fortran code this
        //replaces, the seed is set with index_ + 1 instead of index_
        double r = future_prn(static_cast<int64_t>(index_ + 1), *p.current_seed());
        p.stream_ = STREAM_TRACKING;

        int i_low = 0;
        while (urr.prob(i_energy, URRTableParam::CUM_PROB, i_low) <= r) {++i_low;};

        int i_up = 0;
        while (urr.prob(i_energy + 1, URRTableParam::CUM_PROB, i_up) <= r) {++i_up;};

        // Determine elastic, fission, and capture cross sections from the
        // probability table
        double p_elastic = 0.;
        double p_fission = 0.;
        double p_capture = 0.;
        double p_f;
        if (urr.interp_ == Interpolation::lin_lin) {
          // Determine the interpolation factor on the table
          p_f = (E - urr.device_energy_[i_energy]) /
            (urr.device_energy_[i_energy + 1] - urr.device_energy_[i_energy]);

          p_elastic = (1. - p_f) * urr.prob(i_energy, URRTableParam::ELASTIC, i_low) +
            p_f * urr.prob(i_energy + 1, URRTableParam::ELASTIC, i_up);
          p_fission = (1. - p_f) * urr.prob(i_energy, URRTableParam::FISSION, i_low) +
            p_f * urr.prob(i_energy + 1, URRTableParam::FISSION, i_up);
          p_capture = (1. - p_f) * urr.prob(i_energy, URRTableParam::N_GAMMA, i_low) +
            p_f * urr.prob(i_energy + 1, URRTableParam::N_GAMMA, i_up);
        } else if (urr.interp_ == Interpolation::log_log) {
          // Determine interpolation factor on the table
          p_f = std::log(E / urr.device_energy_[i_energy]) /
            std::log(urr.device_energy_[i_energy + 1] / urr.device_energy_[i_energy]);

          // Calculate the elastic cross section/factor
          if ((urr.prob(i_energy, URRTableParam::ELASTIC, i_low) > 0.) &&
              (urr.prob(i_energy + 1, URRTableParam::ELASTIC, i_up) > 0.)) {
            p_elastic =
              std::exp((1. - p_f) *
                  std::log(urr.prob(i_energy, URRTableParam::ELASTIC, i_low)) +
                  p_f * std::log(urr.prob(i_energy + 1, URRTableParam::ELASTIC, i_up)));
          } else {
            p_elastic = 0.;
          }

          // Calculate the fission cross section/factor
          if ((urr.prob(i_energy, URRTableParam::FISSION, i_low) > 0.) &&
              (urr.prob(i_energy + 1, URRTableParam::FISSION, i_up) > 0.)) {
            p_fission =
              std::exp((1. - p_f) *
                  std::log(urr.prob(i_energy, URRTableParam::FISSION, i_low)) +
                  p_f * std::log(urr.prob(i_energy + 1, URRTableParam::FISSION, i_up)));
          } else {
            p_fission = 0.;
          }

          // Calculate the capture cross section/factor
          if ((urr.prob(i_energy, URRTableParam::N_GAMMA, i_low) > 0.) &&
              (urr.prob(i_energy + 1, URRTableParam::N_GAMMA, i_up) > 0.)) {
            p_capture =
              std::exp((1. - p_f) *
                  std::log(urr.prob(i_energy, URRTableParam::N_GAMMA, i_low)) +
                  p_f * std::log(urr.prob(i_energy + 1, URRTableParam::N_GAMMA, i_up)));
          } else {
            p_capture = 0.;
          }
        }

        // Determine the treatment of inelastic scattering
        double p_inelastic = 0.;
        if (urr.inelastic_flag_ != C_NONE) {
          // Determine inelastic scattering cross section
          auto rx = reactions_[urr_inelastic_].obj();
          p_inelastic = rx.xs(i_temp, i_grid,f);
        }

        // Multiply by smooth cross-section if needed
        if (urr.multiply_smooth_) {
          p_elastic *= this->calculate_elastic_xs(i_temp, i_grid, f);
          p_capture *= (absorption - fission);
          p_fission *= fission;
        }

        // Check for negative values
        if (p_elastic < 0.) {p_elastic = 0.;}
        if (p_fission < 0.) {p_fission = 0.;}
        if (p_capture < 0.) {p_capture = 0.;}

        // Set elastic, absorption, fission, total, and capture x/s. Note that the
        // total x/s is calculated as a sum of partials instead of the table-provided
        // value
        elastic = p_elastic;
        absorption = p_capture + p_fission;
        fission = p_fission;
        total = p_elastic + p_inelastic + p_capture + p_fission;

        if (need_depletion_rx) {
          reaction[0] = p_capture;
        }

        // Determine nu-fission cross-section
        if (fissionable_) {
          nu_fission = nu(E, EmissionMode::total) * fission;
        }
      }
    }
    // ======================================================================
    // URR END
    // ======================================================================

    // ======================================================================
    // Return Struct (depends on template type)
    // ======================================================================

    return T(
        total,
        absorption,
        fission,
        nu_fission,
        elastic,
        thermal,
        thermal_elastic,
        photon_prod,
        reaction,
        i_grid,
        i_temp,
        f,
        index_sab,
        index_temp_sab,
        sab_frac,
        use_ptable,
        E,
        sqrtkT);
  }

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
  WindowedMultipole* device_multipole_ {nullptr};
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
  Function1DFlatContainer* device_total_nu_ {nullptr};

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

  ReactionFlatContainer** device_fission_rx_ {nullptr};
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
