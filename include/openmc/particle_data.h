#ifndef OPENMC_PARTICLE_DATA_H
#define OPENMC_PARTICLE_DATA_H

#include "openmc/array.h"
#include "openmc/constants.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"
#include "openmc/tallies/filter_match.h"
#include "openmc/vector.h"

#ifdef OPENMC_DAGMC_ENABLED
#include "DagMC.hpp"
#endif

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

// Since cross section libraries come with different numbers of delayed groups
// (e.g. ENDF/B-VII.1 has 6 and JEFF 3.1.1 has 8 delayed groups) and we don't
// yet know what cross section library is being used when the tallies.xml file
// is read in, we want to have an upper bound on the size of the array we
// use to store the bins for delayed group tallies.
constexpr int MAX_DELAYED_GROUPS {8};

constexpr double CACHE_INVALID {-1.0};

//==========================================================================
// Aliases and type definitions

//! Particle types
enum class ParticleType { neutron, photon, electron, positron };

//! Saved ("banked") state of a particle
//! NOTE: This structure's MPI type is built in initialize_mpi() of
//! initialize.cpp. Any changes made to the struct here must also be
//! made when building the Bank MPI type in initialize_mpi().
//! NOTE: This structure is also used on the python side, and is defined
//! in lib/core.py. Changes made to the type here must also be made to the
//! python defintion.
struct SourceSite {
  Position r;
  Direction u;
  double E;
  double E_parent;
  double time {0.0};
  double wgt {1.0};
  int delayed_group {0};
  int fission_nuclide;
  int surf_id {SURFACE_NONE};
  ParticleType particle;

  // Extra attributes that don't show up in source written to file
  int parent_nuclide {-1};
  int64_t parent_id;
  int64_t progeny_id;
};

//! State of a particle used for particle track files
struct TrackState {
  Position r;           //!< Position in [cm]
  Direction u;          //!< Direction
  double E;             //!< Energy in [eV]
  double time {0.0};    //!< Time in [s]
  double wgt {1.0};     //!< Weight
  int cell_id;          //!< Cell ID
  int cell_instance;    //!< Cell instance
  int material_id {-1}; //!< Material ID (default value indicates void)
};

//! Full history of a single particle's track states
struct TrackStateHistory {
  ParticleType particle;
  std::vector<TrackState> states;
};

//! Saved ("banked") state of a particle, for nu-fission tallying
struct NuBank {
  double E;          //!< particle energy
  double wgt;        //!< particle weight
  int delayed_group; //!< particle delayed group
};

class LocalCoord {
public:
  void rotate(const vector<double>& rotation);

  //! clear data from a single coordinate level
  void reset();

  // accessors
  Position& r() { return r_; }
  const Position& r() const { return r_; }

  Direction& u() { return u_; }
  const Direction& u() const { return u_; }

  int& cell() { return cell_; }
  const int& cell() const { return cell_; }

  int& universe() { return universe_; }
  const int& universe() const { return universe_; }

  int& lattice() { return lattice_; }
  int lattice() const { return lattice_; }

  array<int, 3>& lattice_index() { return lattice_index_; }
  const array<int, 3>& lattice_index() const { return lattice_index_; }

  bool& rotated() { return rotated_; }
  const bool& rotated() const { return rotated_; }

private:
  // Data members
  Position r_;  //!< particle position
  Direction u_; //!< particle direction
  int cell_ {-1};
  int universe_ {-1};
  int lattice_ {-1};
  array<int, 3> lattice_index_ {{-1, -1, -1}};
  bool rotated_ {false}; //!< Is the level rotated?
};

//==============================================================================
//! Cached microscopic cross sections for a particular nuclide at the current
//! energy
//==============================================================================

struct NuclideMicroXS {
  // Microscopic cross sections in barns
  double total;      //!< total cross section
  double absorption; //!< absorption (disappearance)
  double fission;    //!< fission
  double nu_fission; //!< neutron production from fission

  double elastic;         //!< If sab_frac is not 1 or 0, then this value is
                          //!<   averaged over bound and non-bound nuclei
  double thermal;         //!< Bound thermal elastic & inelastic scattering
  double thermal_elastic; //!< Bound thermal elastic scattering
  double photon_prod;     //!< microscopic photon production xs

  // Cross sections for depletion reactions (note that these are not stored in
  // macroscopic cache)
  double reaction[DEPLETION_RX.size()];

  // Indicies and factors needed to compute cross sections from the data tables
  int index_grid;       //!< Index on nuclide energy grid
  int index_temp;       //!< Temperature index for nuclide
  double interp_factor; //!< Interpolation factor on nuc. energy grid
  int index_sab {-1};   //!< Index in sab_tables
  int index_temp_sab;   //!< Temperature index for sab_tables
  double sab_frac;      //!< Fraction of atoms affected by S(a,b)
  bool use_ptable;      //!< In URR range with probability tables?

  // Energy and temperature last used to evaluate these cross sections.  If
  // these values have changed, then the cross sections must be re-evaluated.
  double last_E {0.0};       //!< Last evaluated energy
  double last_sqrtkT {0.0};  //!< Last temperature in sqrt(Boltzmann constant
                             //!< * temperature (eV))
  double ncrystal_xs {-1.0}; //!< NCrystal cross section
};

//==============================================================================
//! Cached microscopic photon cross sections for a particular element at the
//! current energy
//==============================================================================

struct ElementMicroXS {
  int index_grid;         //!< index on element energy grid
  double last_E {0.0};    //!< last evaluated energy in [eV]
  double interp_factor;   //!< interpolation factor on energy grid
  double total;           //!< microscopic total photon xs
  double coherent;        //!< microscopic coherent xs
  double incoherent;      //!< microscopic incoherent xs
  double photoelectric;   //!< microscopic photoelectric xs
  double pair_production; //!< microscopic pair production xs
};

//==============================================================================
// MacroXS contains cached macroscopic cross sections for the material a
// particle is traveling through
//==============================================================================

struct MacroXS {
  double total;       //!< macroscopic total xs
  double absorption;  //!< macroscopic absorption xs
  double fission;     //!< macroscopic fission xs
  double nu_fission;  //!< macroscopic production xs
  double photon_prod; //!< macroscopic photon production xs

  // Photon cross sections
  double coherent;        //!< macroscopic coherent xs
  double incoherent;      //!< macroscopic incoherent xs
  double photoelectric;   //!< macroscopic photoelectric xs
  double pair_production; //!< macroscopic pair production xs
};

//==============================================================================
// Cache contains the cached data for an MGXS object
//==============================================================================

struct CacheDataMG {
  int material {-1}; //!< material index
  double sqrtkT;     //!< last temperature corresponding to t
  int t {0};         //!< temperature index
  int a {0};         //!< angle index
  Direction u;       //!< angle that corresponds to a
};

//==============================================================================
// Information about nearest boundary crossing
//==============================================================================

class BoundaryInfo {
public:
  void reset()
  {
    distance_ = INFINITY;
    surface_ = SURFACE_NONE;
    coord_level_ = 0;
    lattice_translation_ = {0, 0, 0};
  }
  double& distance() { return distance_; }
  const double& distance() const { return distance_; }

  int& surface() { return surface_; }
  const int& surface() const { return surface_; }

  int coord_level() const { return coord_level_; }
  int& coord_level() { return coord_level_; }

  array<int, 3>& lattice_translation() { return lattice_translation_; }
  const array<int, 3>& lattice_translation() const
  {
    return lattice_translation_;
  }

  // TODO: off-by-one
  int surface_index() const { return std::abs(surface()) - 1; }

private:
  // Data members
  double distance_ {INFINITY}; //!< distance to nearest boundary
  int surface_ {
    SURFACE_NONE};      //!< surface token, non-zero if boundary is surface
  int coord_level_ {0}; //!< coordinate level after crossing boundary
  array<int, 3> lattice_translation_ {
    0, 0, 0}; //!< which way lattice indices will change
};

/*
 * Contains all geometry state information for a particle.
 */
class GeometryState {
public:
  GeometryState();

  /*
   * GeometryState does not store any ID info, so give some reasonable behavior
   * here. The Particle class redefines this. This is only here for the error
   * reporting behavior that occurs in geometry.cpp. The explanation for
   * mark_as_lost is the same.
   */
  virtual void mark_as_lost(const char* message);
  void mark_as_lost(const std::string& message);
  void mark_as_lost(const std::stringstream& message);

  // resets all coordinate levels for the particle
  void clear()
  {
    for (auto& level : coord_) {
      level.reset();
    }
    n_coord_ = 1;

    for (auto& cell : cell_last_) {
      cell = C_NONE;
    }
    n_coord_last_ = 1;
  }

  //! moves the particle by the specified distance to its next location
  //! \param distance the distance the particle is moved
  void move_distance(double distance);

  void advance_to_boundary_from_void();

  // Initialize all internal state from position and direction
  void init_from_r_u(Position r_a, Direction u_a)
  {
    clear();
    surface() = SURFACE_NONE;
    material() = C_NONE;
    r() = r_a;
    u() = u_a;
    r_last_current() = r_a;
    r_last() = r_a;
    u_last() = u_a;
  }

  // Unique ID. This is not geometric info, but the
  // error reporting in geometry.cpp requires this.
  // We could save this to implement it in Particle,
  // but that would require virtuals.
  int64_t& id() { return id_; }
  const int64_t& id() const { return id_; }

  // Number of current coordinate levels
  int& n_coord() { return n_coord_; }
  const int& n_coord() const { return n_coord_; }

  // Offset for distributed properties
  int& cell_instance() { return cell_instance_; }
  const int& cell_instance() const { return cell_instance_; }

  // Coordinates for all nesting levels
  LocalCoord& coord(int i) { return coord_[i]; }
  const LocalCoord& coord(int i) const { return coord_[i]; }
  const vector<LocalCoord>& coord() const { return coord_; }

  // Innermost universe nesting coordinates
  LocalCoord& lowest_coord() { return coord_[n_coord_ - 1]; }
  const LocalCoord& lowest_coord() const { return coord_[n_coord_ - 1]; }

  // Last coordinates on all nesting levels, before crossing a surface
  int& n_coord_last() { return n_coord_last_; }
  const int& n_coord_last() const { return n_coord_last_; }
  int& cell_last(int i) { return cell_last_[i]; }
  const int& cell_last(int i) const { return cell_last_[i]; }

  // Coordinates at birth
  Position& r_born() { return r_born_; }
  const Position& r_born() const { return r_born_; }

  // Coordinates of last collision or reflective/periodic surface
  // crossing for current tallies
  Position& r_last_current() { return r_last_current_; }
  const Position& r_last_current() const { return r_last_current_; }

  // Previous direction and spatial coordinates before a collision
  Position& r_last() { return r_last_; }
  const Position& r_last() const { return r_last_; }
  Position& u_last() { return u_last_; }
  const Position& u_last() const { return u_last_; }

  // Accessors for position in global coordinates
  Position& r() { return coord_[0].r(); }
  const Position& r() const { return coord_[0].r(); }

  // Accessors for position in local coordinates
  Position& r_local() { return coord_[n_coord_ - 1].r(); }
  const Position& r_local() const { return coord_[n_coord_ - 1].r(); }

  // Accessors for direction in global coordinates
  Direction& u() { return coord_[0].u(); }
  const Direction& u() const { return coord_[0].u(); }

  // Accessors for direction in local coordinates
  Direction& u_local() { return coord_[n_coord_ - 1].u(); }
  const Direction& u_local() const { return coord_[n_coord_ - 1].u(); }

  // Surface token for the surface that the particle is currently on
  int& surface() { return surface_; }
  const int& surface() const { return surface_; }

  // Surface index based on the current value of the surface_ attribute
  int surface_index() const
  {
    // TODO: off-by-one
    return std::abs(surface_) - 1;
  }

  // Boundary information
  BoundaryInfo& boundary() { return boundary_; }

#ifdef OPENMC_DAGMC_ENABLED
  // DagMC state variables
  moab::DagMC::RayHistory& history() { return history_; }
  Direction& last_dir() { return last_dir_; }
#endif

  // material of current and last cell
  int& material() { return material_; }
  const int& material() const { return material_; }
  int& material_last() { return material_last_; }
  const int& material_last() const { return material_last_; }

  // temperature of current and last cell
  double& sqrtkT() { return sqrtkT_; }
  const double& sqrtkT() const { return sqrtkT_; }
  double& sqrtkT_last() { return sqrtkT_last_; }

  // density multiplier of the current and last cell
  double& density_mult() { return density_mult_; }
  const double& density_mult() const { return density_mult_; }
  double& density_mult_last() { return density_mult_last_; }

private:
  int64_t id_ {-1}; //!< Unique ID

  int n_coord_ {1};          //!< number of current coordinate levels
  int cell_instance_;        //!< offset for distributed properties
  vector<LocalCoord> coord_; //!< coordinates for all levels

  int n_coord_last_ {1};  //!< number of current coordinates
  vector<int> cell_last_; //!< coordinates for all levels

  Position r_born_;         //!< coordinates at birth
  Position r_last_current_; //!< coordinates of the last collision or
                            //!< reflective/periodic surface crossing for
                            //!< current tallies
  Position r_last_;         //!< previous coordinates
  Direction u_last_;        //!< previous direction coordinates

  int surface_ {
    SURFACE_NONE}; //!< surface token for surface the particle is currently on

  BoundaryInfo boundary_; //!< Info about the next intersection

  int material_ {-1};      //!< index for current material
  int material_last_ {-1}; //!< index for last material

  double sqrtkT_ {-1.0};     //!< sqrt(k_Boltzmann * temperature) in eV
  double sqrtkT_last_ {0.0}; //!< last temperature

  double density_mult_ {1.0};      //!< density multiplier
  double density_mult_last_ {1.0}; //!< last density multiplier

#ifdef OPENMC_DAGMC_ENABLED
  moab::DagMC::RayHistory history_;
  Direction last_dir_;
#endif
};

//============================================================================
//! Defines how particle data is laid out in memory
//============================================================================

/*
 * This class was added in order to separate the layout and access of particle
 * data from particle physics operations during a development effort to get
 * OpenMC running on GPUs. In the event-based Monte Carlo method, one creates
 * an array of particles on which actions like cross section lookup and surface
 * crossing are done en masse, which works best on vector computers of yore and
 * modern GPUs. It has been shown in the below publication [1] that arranging
 * particle data into a structure of arrays rather than an array of structures
 * enhances performance on GPUs. For instance, rather than having an
 * std::vector<Particle> where consecutive particle energies would be separated
 * by about 400 bytes, one would create a structure which has a single
 * std::vector<double> of energies.  The motivation here is that more coalesced
 * memory accesses occur, in the parlance of GPU programming.
 *
 * So, this class enables switching between the array-of-structures and
 * structure- of-array data layout at compile time. In GPU branches of the
 * code, our Particle class inherits from a class that provides an array of
 * particle energies, and can access them using the E() method (defined below).
 * In the CPU code, we inherit from this class which gives the conventional
 * layout of particle data, useful for history-based tracking.
 *
 * As a result, we always use the E(), r_last(), etc. methods to access
 * particle data in order to keep a unified interface between
 * structure-of-array and array-of-structure code on either CPU or GPU code
 * while sharing the same physics code on each codebase.
 *
 * [1] Hamilton, Steven P., Stuart R. Slattery, and Thomas M. Evans.
 *   “Multigroup Monte Carlo on GPUs: Comparison of History- and Event-Based
 *   Algorithms.” Annals of Nuclear Energy 113 (March 2018): 506–18.
 *   https://doi.org/10.1016/j.anucene.2017.11.032.
 */
class ParticleData : public GeometryState {
private:
  //==========================================================================
  // Data members -- see public: below for descriptions

  vector<NuclideMicroXS> neutron_xs_;
  vector<ElementMicroXS> photon_xs_;
  MacroXS macro_xs_;
  CacheDataMG mg_xs_cache_;

  ParticleType type_ {ParticleType::neutron};

  double E_;
  double E_last_;
  double E_parent_; //!< energy of parent neutron in eV
  int g_ {0};
  int g_last_;
  
  // Other birth data
  int fission_nuclide_;  //!< this particle was born as a result of this nuclide fissioning
  // a double for fission cross section at birth? if so, I need to also add it to the bank...

  double wgt_ {1.0};
  double wgt_born_ {1.0};
  double wgt_ww_born_ {-1.0};
  double mu_;
  double time_ {0.0};
  double time_last_ {0.0};
  double wgt_last_ {1.0};

  bool fission_ {false};
  TallyEvent event_;
  int event_nuclide_;
  int event_mt_;
  int delayed_group_ {0};
  int parent_nuclide_ {-1};

  int n_bank_ {0};
  double bank_second_E_ {0.0};
  double wgt_bank_ {0.0};
  int n_delayed_bank_[MAX_DELAYED_GROUPS];

  int cell_born_ {-1};

  // Iterated Fission Probability
  double lifetime_ {0.0}; //!< neutron lifetime [s]

  int n_collision_ {0};

  bool write_track_ {false};

  uint64_t seeds_[N_STREAMS];
  int stream_;

  vector<SourceSite> secondary_bank_;

  int64_t current_work_;

  vector<double> flux_derivs_;
  
  std::vector<std::vector<double>> cumulative_sensitivities_;  // for sensitivities for this particle

  vector<FilterMatch> filter_matches_;

  vector<TrackStateHistory> tracks_;

  vector<NuBank> nu_bank_;

  vector<double> pht_storage_;

  double keff_tally_absorption_ {0.0};
  double keff_tally_collision_ {0.0};
  double keff_tally_tracklength_ {0.0};
  double keff_tally_leakage_ {0.0};

  bool trace_ {false};

  double collision_distance_;

  int n_event_ {0};

  int n_split_ {0};
  double ww_factor_ {0.0};

  int64_t n_progeny_ {0};

public:
  //----------------------------------------------------------------------------
  // Constructors
  ParticleData();

  //==========================================================================
  // Methods and accessors

  // Cross section caches
  NuclideMicroXS& neutron_xs(int i)
  {
    return neutron_xs_[i];
  } // Microscopic neutron cross sections
  const NuclideMicroXS& neutron_xs(int i) const { return neutron_xs_[i]; }

  // Microscopic photon cross sections
  ElementMicroXS& photon_xs(int i) { return photon_xs_[i]; }

  // Macroscopic cross sections
  MacroXS& macro_xs() { return macro_xs_; }
  const MacroXS& macro_xs() const { return macro_xs_; }

  // Multigroup macroscopic cross sections
  CacheDataMG& mg_xs_cache() { return mg_xs_cache_; }
  const CacheDataMG& mg_xs_cache() const { return mg_xs_cache_; }

  // Particle type (n, p, e, gamma, etc)
  ParticleType& type() { return type_; }
  const ParticleType& type() const { return type_; }

  // Current particle energy, energy before collision,
  // and corresponding multigroup group indices. Energy
  // units are eV.
  double& E() { return E_; }
  const double& E() const { return E_; }
  double& E_last() { return E_last_; }
  const double& E_last() const { return E_last_; }
  double& E_parent() { return E_parent_; }              // for sensitivity analysis
  const double& E_parent() const { return E_parent_; }  // for SA
  int& g() { return g_; }
  const int& g() const { return g_; }
  int& g_last() { return g_last_; }
  const int& g_last() const { return g_last_; }

  // Statistic weight of particle. Setting to zero indicates that the particle
  // is dead.
  double& wgt() { return wgt_; }
  double wgt() const { return wgt_; }

  // Statistic weight of particle at birth
  double& wgt_born() { return wgt_born_; }
  double wgt_born() const { return wgt_born_; }

  // Weight window value at birth
  double& wgt_ww_born() { return wgt_ww_born_; }
  const double& wgt_ww_born() const { return wgt_ww_born_; }

  // Statistic weight of particle at last collision
  double& wgt_last() { return wgt_last_; }
  const double& wgt_last() const { return wgt_last_; }

  // Whether particle is alive
  bool alive() const { return wgt_ != 0.0; }

  // Polar scattering angle after a collision
  double& mu() { return mu_; }
  const double& mu() const { return mu_; }

  // Tracks the time of a particle as it traverses the problem.
  // Units are seconds.
  double& time() { return time_; }
  const double& time() const { return time_; }
  double& time_last() { return time_last_; }
  const double& time_last() const { return time_last_; }

  // Particle lifetime
  double& lifetime() { return lifetime_; }
  const double& lifetime() const { return lifetime_; }

  // What event took place, described in greater detail below
  TallyEvent& event() { return event_; }
  const TallyEvent& event() const { return event_; }
  bool& fission() { return fission_; }            // true if implicit fission
  int& event_nuclide() { return event_nuclide_; } // index of collision nuclide
  const int& event_nuclide() const { return event_nuclide_; }
  int& event_mt() { return event_mt_; } // MT number of collision
  const int& event_mt() const { return event_mt_; }
  int& delayed_group() { return delayed_group_; } // delayed group
  const int& delayed_group() const { return delayed_group_; }
  const int& parent_nuclide() const { return parent_nuclide_; }
  int& parent_nuclide() { return parent_nuclide_; } // Parent nuclide

  // Post-collision data
  double& bank_second_E()
  {
    return bank_second_E_;
  } // energy of last reaction secondaries
  const double& bank_second_E() const { return bank_second_E_; }

  int& n_bank() { return n_bank_; }        // number of banked fission sites
  double& wgt_bank() { return wgt_bank_; } // weight of banked fission sites
  int* n_delayed_bank()
  {
    return n_delayed_bank_;
  } // number of delayed fission sites
  int& n_delayed_bank(int i)
  {
    return n_delayed_bank_[i];
  } // number of delayed fission sites

  // Index of cell particle is born in
  int& cell_born() { return cell_born_; }
  const int& cell_born() const { return cell_born_; }

  // Total number of collisions suffered by particle
  int& n_collision() { return n_collision_; }
  const int& n_collision() const { return n_collision_; }

  // whether this track is to be written
  bool& write_track() { return write_track_; }

  // RNG state
  uint64_t& seeds(int i) { return seeds_[i]; }
  uint64_t* seeds() { return seeds_; }
  int& stream() { return stream_; }

  // secondary particle bank
  SourceSite& secondary_bank(int i) { return secondary_bank_[i]; }
  decltype(secondary_bank_)& secondary_bank() { return secondary_bank_; }

  // Current simulation work index
  int64_t& current_work() { return current_work_; }
  const int64_t& current_work() const { return current_work_; }

  // Used in tally derivatives
  double& flux_derivs(int i) { return flux_derivs_[i]; }
  const double& flux_derivs(int i) const { return flux_derivs_[i]; }
  
  // Used in sensitivity analysis
  std::vector<double>& cumulative_sensitivities(int i) { return cumulative_sensitivities_[i]; }
  const std::vector<double>& cumulative_sensitivities(int i) const { return cumulative_sensitivities_[i]; }
  int& fission_nuclide() { return fission_nuclide_; }
  const int& fission_nuclide() const { return fission_nuclide_; }

  // Matches of tallies
  decltype(filter_matches_)& filter_matches() { return filter_matches_; }
  FilterMatch& filter_matches(int i) { return filter_matches_[i]; }

  // Tracks to output to file
  decltype(tracks_)& tracks() { return tracks_; }

  // Bank of recently fissioned particles
  decltype(nu_bank_)& nu_bank() { return nu_bank_; }
  NuBank& nu_bank(int i) { return nu_bank_[i]; }

  // Interim pulse height tally storage
  vector<double>& pht_storage() { return pht_storage_; }

  // Global tally accumulators
  double& keff_tally_absorption() { return keff_tally_absorption_; }
  double& keff_tally_collision() { return keff_tally_collision_; }
  double& keff_tally_tracklength() { return keff_tally_tracklength_; }
  double& keff_tally_leakage() { return keff_tally_leakage_; }

  // Shows debug info
  bool& trace() { return trace_; }

  // Distance to the next collision
  double& collision_distance() { return collision_distance_; }

  // Number of events particle has undergone
  int& n_event() { return n_event_; }

  // Number of times variance reduction has caused a particle split
  int n_split() const { return n_split_; }
  int& n_split() { return n_split_; }

  // Particle-specific factor for on-the-fly weight window adjustment
  double ww_factor() const { return ww_factor_; }
  double& ww_factor() { return ww_factor_; }

  // Number of progeny produced by this particle
  int64_t& n_progeny() { return n_progeny_; }

  //! Gets the pointer to the particle's current PRN seed
  uint64_t* current_seed() { return seeds_ + stream_; }
  const uint64_t* current_seed() const { return seeds_ + stream_; }

  //! Force recalculation of neutron xs by setting last energy to zero
  void invalidate_neutron_xs()
  {
    for (auto& micro : neutron_xs_)
      micro.last_E = 0.0;
  }

  //! Get track information based on particle's current state
  TrackState get_track_state() const;

  void zero_delayed_bank()
  {
    for (int& n : n_delayed_bank_) {
      n = 0;
    }
  }

  void zero_flux_derivs()
  {
    for (double& d : flux_derivs_) {
      d = 0;
    }
  }
  
  void resize_flux_derivs(int newSize)
  {
    flux_derivs_.resize(newSize, 0.0);
  }
    
  // Resize and initialize sensitivity vectors
  void resize_init_cumulative_sensitivities(int newSize)
  {
    cumulative_sensitivities_.resize(newSize, {0.0});
  }
  
  void resize_init_cumulative_sensitivities_vec(int indx, int newSize)
  {
    cumulative_sensitivities_[indx].resize(newSize, 0.0);
  }
      
  void resize_alloc_filter_matches(int newSize)
  {
    filter_matches_.resize(newSize);
    for (auto &match: filter_matches_){
      match.bins_.resize(1, 0.0);
      match.weights_.resize(1, 0.0);
      match.i_bin_ = 0;
      // bins_present_ default is false from header
    }
  }
  
  void initialize_cumulative_sensitivities()
  {
    for (auto& it : cumulative_sensitivities_){
      std::fill(it.begin(), it.end(), 0.0);
    }
  }
};

} // namespace openmc

#endif // OPENMC_PARTICLE_DATA_H
