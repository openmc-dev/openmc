#ifndef OPENMC_PARTICLE_H
#define OPENMC_PARTICLE_H

//! \file particle.h
//! \brief Particle type

#include <array>
#include <cstdint>
#include <memory> // for unique_ptr
#include <sstream>
#include <string>
#include <vector>

#include "openmc/constants.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/tallies/filter_match.h"

#ifdef DAGMC
#include "DagMC.hpp"
#endif

/*
//#define NEUTRON_XS_SIZE 300 // Depleted SMR uses 296
#define NEUTRON_XS_SIZE 35 // Depleted SMR uses 296
#define PHOTON_XS_SIZE 9 // Pincell example uses 9
#define COORD_SIZE 6 // Depleted SMR uses 6
//#define SECONDARY_BANK_SIZE 200 // 100 not enough to pass regression tests, but 200 works. TODO: narrow this down.
#define SECONDARY_BANK_SIZE 10 // 100 not enough to pass regression tests, but 200 works. TODO: narrow this down.
#define FLUX_DERIVS_SIZE 5 // This is the min required to pass regression tests (diff_tally is limiter)
#define FILTER_MATCHES_SIZE 14 // tallies regression test is the limiter here. More realistic tests only need 2. This can be set at runtime init though.
//#define FILTER_MATCHES_SIZE 140 // tallies regression test is the limiter here. More realistic tests only need 2. This can be set at runtime init though.
#define NU_BANK_SIZE 16 // infinite_cell regression test
*/
// Minimal for HM-SMall
//#define NEUTRON_XS_SIZE 69 // HM-Small
//#define NEUTRON_XS_SIZE 272 // HM-Large
//#define NEUTRON_XS_SIZE 92 // SMR
#ifdef NO_MICRO_XS_CACHE
#define NEUTRON_XS_SIZE 1
#else
#define NEUTRON_XS_SIZE 296 // Depleted SMR uses 296
#endif
//#define NEUTRON_XS_SIZE 35 // Fresh Pincell
#define PHOTON_XS_SIZE 9 // Pincell example uses 9
#define COORD_SIZE 6 // Depleted SMR uses 6
//#define COORD_SIZE 3 // Depleted SMR uses 6
//#define SECONDARY_BANK_SIZE 200 // 100 not enough to pass regression tests, but 200 works. TODO: narrow this down.
#define SECONDARY_BANK_SIZE 5 // 100 not enough to pass regression tests, but 200 works. TODO: narrow this down.
#define FLUX_DERIVS_SIZE 1 // This is the min required to pass regression tests (diff_tally is limiter)
#define NU_BANK_SIZE 16 // infinite_cell regression test

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

//==============================================================================
// Class declarations
//==============================================================================

// Forward declare the Surface class for use in function arguments.
class Surface;

class LocalCoord {
public:
  //void rotate(const std::vector<double>& rotation);
  #pragma omp declare target
  void rotate(const double* rotation);
  #pragma omp end declare target

  //! clear data from a single coordinate level
  #pragma omp declare target
  void reset();
  #pragma omp end declare target

  Position r; //!< particle position
  Direction u; //!< particle direction
  int cell {-1};
  int universe {-1};
  int lattice {-1};
  int lattice_x {-1};
  int lattice_y {-1};
  int lattice_z {-1};
  bool rotated {false};  //!< Is the level rotated?
};


//==============================================================================
//! Cached microscopic cross sections for a particular nuclide at the current
//! energy
//==============================================================================
struct NuclideMicroXS {
  // Microscopic cross sections in barns
  double total;            //!< total cross section
  double absorption;       //!< absorption (disappearance)
  double fission;          //!< fission
  double nu_fission;       //!< neutron production from fission

  double elastic;          //!< If sab_frac is not 1 or 0, then this value is
                           //!<   averaged over bound and non-bound nuclei
  double thermal;          //!< Bound thermal elastic & inelastic scattering
  double thermal_elastic;  //!< Bound thermal elastic scattering
  double photon_prod;      //!< microscopic photon production xs

  // Cross sections for depletion reactions (note that these are not stored in
  // macroscopic cache)
  double reaction[DEPLETION_RX_SIZE];

  // Indicies and factors needed to compute cross sections from the data tables
  int index_grid;        //!< Index on nuclide energy grid
  int index_temp;        //!< Temperature index for nuclide
  double interp_factor;  //!< Interpolation factor on nuc. energy grid
  int index_sab {-1};    //!< Index in sab_tables
  int index_temp_sab;    //!< Temperature index for sab_tables
  double sab_frac;       //!< Fraction of atoms affected by S(a,b)
  bool use_ptable;       //!< In URR range with probability tables?

  // Energy and temperature last used to evaluate these cross sections.  If
  // these values have changed, then the cross sections must be re-evaluated.
  double last_E {0.0};      //!< Last evaluated energy
  double last_sqrtkT {0.0}; //!< Last temperature in sqrt(Boltzmann constant
                            //!< * temperature (eV))
                            
  NuclideMicroXS() = default;

  NuclideMicroXS(
      double total,
      double absorption,
      double fission,
      double nu_fission,
      double elastic,
      double thermal,
      double thermal_elastic,
      double photon_prod,
      double reaction_in[DEPLETION_RX_SIZE],
      int index_grid,
      int index_temp,
      double interp_factor,
      int index_sab,
      int index_temp_sab,
      double sab_frac,
      bool use_ptable,
      double last_E,
      double las_sqrtkT
      ) :
    total(total),
    absorption(absorption),
    fission(fission),
    nu_fission(nu_fission),
    elastic( elastic),
    thermal(thermal),
    thermal_elastic(thermal_elastic),
    photon_prod(photon_prod),
    index_grid(index_grid),
    index_temp(index_temp),
    interp_factor(interp_factor),
    index_sab(index_sab),
    index_temp_sab(index_temp_sab),
    sab_frac(sab_frac),
    use_ptable(use_ptable),
    last_E(last_E),
    last_sqrtkT(las_sqrtkT)
  {
    for (int r = 0; r < DEPLETION_RX_SIZE; r++) {
      reaction[r] = reaction_in[r];
    }
  }
};

class NuclideMicroXSCache {
  public:
  NuclideMicroXS neutron_xs_[NEUTRON_XS_SIZE]; //!< Microscopic neutron cross sections
  #pragma omp declare target
  NuclideMicroXS  operator [](int64_t i) const
  {
    #ifdef NO_MICRO_XS_CACHE
    return neutron_xs_[0];
    #else
    return neutron_xs_[i];
    #endif
  }
  NuclideMicroXS& operator [](int64_t i)
  {
    #ifdef NO_MICRO_XS_CACHE
    return neutron_xs_[0];
    #else
    return neutron_xs_[i];
    #endif
  }
  #pragma omp end declare target
  void clear()
  {
    if (settings::run_CE) {
      for (auto& micro : neutron_xs_) micro.last_E = 0.0;
    }
  }
};

//==============================================================================
//! Cached microscopic photon cross sections for a particular element at the
//! current energy
//==============================================================================

struct ElementMicroXS {
  int index_grid; //!< index on element energy grid
  double last_E {0.0}; //!< last evaluated energy in [eV]
  double interp_factor; //!< interpolation factor on energy grid
  double total; //!< microscopic total photon xs
  double coherent; //!< microscopic coherent xs
  double incoherent; //!< microscopic incoherent xs
  double photoelectric; //!< microscopic photoelectric xs
  double pair_production; //!< microscopic pair production xs
};

//==============================================================================
// MACROXS contains cached macroscopic cross sections for the material a
// particle is traveling through
//==============================================================================

struct MacroXS {
  double total;         //!< macroscopic total xs
  double absorption;    //!< macroscopic absorption xs
  double fission;       //!< macroscopic fission xs
  double nu_fission;    //!< macroscopic production xs
  double photon_prod;   //!< macroscopic photon production xs
  
  // Cross sections for depletion reactions
  double reaction[DEPLETION_RX_SIZE];

  // Photon cross sections
  double coherent;        //!< macroscopic coherent xs
  double incoherent;      //!< macroscopic incoherent xs
  double photoelectric;   //!< macroscopic photoelectric xs
  double pair_production; //!< macroscopic pair production xs
};

//==============================================================================
// MICROXS contains microscopic neutron cross sections for a nuclide
//==============================================================================

struct MicroXS {
  double total;         //!< macroscopic total xs
  double absorption;    //!< macroscopic absorption xs
  double fission;       //!< macroscopic fission xs
  double nu_fission;    //!< macroscopic production xs
  double reaction[DEPLETION_RX_SIZE];

  MicroXS() = default;

  MicroXS(
      double total,
      double absorption,
      double fission,
      double nu_fission,
      double elastic,
      double thermal,
      double thermal_elastic,
      double photon_prod,
      double reaction_in[DEPLETION_RX_SIZE],
      int index_grid,
      int index_temp,
      double interp_factor,
      int index_sab,
      int index_temp_sab,
      double sab_frac,
      bool use_ptable,
      double last_E,
      double last_sqrtkT
      ) :
    total(total),
    absorption(absorption),
    fission(fission),
    nu_fission(nu_fission)
  {
    for (int r = 0; r < DEPLETION_RX_SIZE; r++) {
      reaction[r] = reaction_in[r];
    }
  }

};

//==============================================================================
// Information about nearest boundary crossing
//==============================================================================

struct BoundaryInfo {
  double distance {INFINITY};   //!< distance to nearest boundary
  int surface_index {0}; //!< if boundary is surface, index in surfaces vector
  int coord_level;   //!< coordinate level after crossing boundary
  std::array<int, 3> lattice_translation {}; //!< which way lattice indices will change
};

//============================================================================
//! State of a particle being transported through geometry
//============================================================================

class Particle {
public:
  //==========================================================================
  // Aliases and type definitions

  //! Particle types
  enum class Type {
    neutron, photon, electron, positron
  };

  //! Saved ("banked") state of a particle
  //! NOTE: This structure's MPI type is built in initialize_mpi() of
  //! initialize.cpp. Any changes made to the struct here must also be
  //! made when building the Bank MPI type in initialize_mpi().
  //! NOTE: This structure is also used on the python side, and is defined
  //! in lib/core.py. Changes made to the type here must also be made to the
  //! python defintion.
  struct Bank {
    Position r;
    Direction u;
    double E;
    double wgt;
    int delayed_group;
    int surf_id;
    Type particle;
    int64_t parent_id;
    int64_t progeny_id;
  };

  //! Saved ("banked") state of a particle, for nu-fission tallying
  struct NuBank {
    double E;  //!< particle energy
    double wgt; //!< particle weight
    int delayed_group; //!< particle delayed group
  };

  //==========================================================================
  // Constructors

  #pragma omp declare target
  Particle();
  #pragma omp end declare target

  //==========================================================================
  // Methods and accessors

  // Accessors for position in global coordinates
  #pragma omp declare target
  Position& r() { return coord_[0].r; }
  const Position& r() const { return coord_[0].r; }

  // Accessors for position in local coordinates
  Position& r_local() { return coord_[n_coord_ - 1].r; }
  const Position& r_local() const { return coord_[n_coord_ - 1].r; }

  // Accessors for direction in global coordinates
  Direction& u() { return coord_[0].u; }
  const Direction& u() const { return coord_[0].u; }

  // Accessors for direction in local coordinates
  Direction& u_local() { return coord_[n_coord_ - 1].u; }
  const Direction& u_local() const { return coord_[n_coord_ - 1].u; }

  //! Whether particle is alive
  bool alive() const { return wgt_ != 0.0; }

  //! resets all coordinate levels for the particle
  void clear();
  #pragma omp end declare target

  //! create a secondary particle
  //
  //! stores the current phase space attributes of the particle in the
  //! secondary bank and increments the number of sites in the secondary bank.
  //! \param wgt Weight of the secondary particle
  //! \param u Direction of the secondary particle
  //! \param E Energy of the secondary particle in [eV]
  //! \param type Particle type
  #pragma omp declare target
  void create_secondary(double wgt, Direction u, double E, Type type);
  #pragma omp end declare target

  //! initialize from a source site
  //
  //! initializes a particle from data stored in a source site. The source
  //! site may have been produced from an external source, from fission, or
  //! simply as a secondary particle.
  //! \param src Source site data
  #pragma omp declare target
  void from_source(const Bank& src);
  #pragma omp end declare target

  // Coarse-grained particle events
  #pragma omp declare target
  void event_tracklength_tally(bool need_depletion_rx);
  void event_calculate_xs(bool need_depletion_rx);
  bool event_calculate_xs_dispatch();
  void event_calculate_xs_execute(bool need_depletion_rx);
  void event_collide();
  void event_advance();
  void event_cross_surface();
  void event_revive_from_secondary();
  void event_death();
  void accumulate_keff_tallies_global();
  void accumulate_keff_tallies_local(double& absorption, double& collision, double& tracklength, double& leakage);
  #pragma omp end declare target

  //! Cross a surface and handle boundary conditions
  #pragma omp declare target
  void cross_surface();
  #pragma omp end declare target

  //! Cross a vacuum boundary condition.
  //
  //! \param surf The surface (with the vacuum boundary condition) that the
  //!   particle struck.
  #pragma omp declare target
  void cross_vacuum_bc(const Surface& surf);
  #pragma omp end declare target

  //! Cross a reflective boundary condition.
  //
  //! \param surf The surface (with the reflective boundary condition) that the
  //!   particle struck.
  //! \param new_u The direction of the particle after reflection.
  #pragma omp declare target
  void cross_reflective_bc(const Surface& surf, Direction new_u);
  #pragma omp end declare target


  //! Cross a periodic boundary condition.
  //
  //! \param surf The surface (with the periodic boundary condition) that the
  //!   particle struck.
  //! \param new_r The position of the particle after translation/rotation.
  //! \param new_u The direction of the particle after translation/rotation.
  //! \param new_surface The signed index of the surface that the particle will
  //!   reside on after translation/rotation.
  void cross_periodic_bc(const Surface& surf, Position new_r, Direction new_u,
                         int new_surface);

  //! mark a particle as lost and create a particle restart file
  //! \param message A warning message to display
  void mark_as_lost(const char* message);
  #pragma omp declare target
  void mark_as_lost_short();
  #pragma omp end declare target

  void mark_as_lost(const std::string& message)
  {mark_as_lost(message.c_str());}

  void mark_as_lost(const std::stringstream& message)
  {mark_as_lost(message.str());}

  //! create a particle restart HDF5 file
  void write_restart() const;

  //! Gets the pointer to the particle's current PRN seed
  #pragma omp declare target
  uint64_t* current_seed() {return seeds_ + stream_;}
  const uint64_t* current_seed() const {return seeds_ + stream_;}
  #pragma omp end declare target

  //==========================================================================
  // Data members

  // Cross section caches

  // TODO: neutron_xs_ can eventually be converted to an allocated array, with size fixed at runtime
  //std::vector<NuclideMicroXS> neutron_xs_; //!< Microscopic neutron cross sections
  //NuclideMicroXS neutron_xs_[NEUTRON_XS_SIZE]; //!< Microscopic neutron cross sections
  NuclideMicroXSCache neutron_xs_; //!< Microscopic neutron cross sections

  // TODO: photon_xs_ can eventually be converted to an allocated array, with size fixed at runtime
  //std::vector<ElementMicroXS> photon_xs_; //!< Microscopic photon cross sections
  ElementMicroXS photon_xs_[PHOTON_XS_SIZE]; //!< Microscopic photon cross sections

  MacroXS macro_xs_; //!< Macroscopic cross sections

  int64_t id_;  //!< Unique ID
  Type type_ {Type::neutron};   //!< Particle type (n, p, e, etc.)

  int n_coord_ {1};              //!< number of current coordinate levels
  int cell_instance_;            //!< offset for distributed properties

  // TODO: coord_ can eventually be converted to an allocated array, with size fixed at runtime
  //std::vector<LocalCoord> coord_; //!< coordinates for all levels
  LocalCoord coord_[COORD_SIZE]; //!< coordinates for all levels

  // Particle coordinates before crossing a surface
  int n_coord_last_ {1};      //!< number of current coordinates
  // TODO: cell_last__ can eventually be converted to an allocated array, with size fixed at runtime
  //std::vector<int> cell_last_;  //!< coordinates for all levels
  int cell_last_[COORD_SIZE];  //!< coordinates for all levels

  // Energy data
  double E_;       //!< post-collision energy in eV
  double E_last_;  //!< pre-collision energy in eV
  int g_ {0};      //!< post-collision energy group (MG only)
  int g_last_;     //!< pre-collision energy group (MG only)

  // Other physical data
  double wgt_ {1.0};     //!< particle weight
  double mu_;      //!< angle of scatter

  // Other physical data
  Position r_last_current_; //!< coordinates of the last collision or
                            //!< reflective/periodic surface crossing for
                            //!< current tallies
  Position r_last_;   //!< previous coordinates
  Direction u_last_;  //!< previous direction coordinates
  double wgt_last_ {1.0};   //!< pre-collision particle weight

  // What event took place
  bool fission_ {false}; //!< did particle cause implicit fission
  TallyEvent event_;          //!< scatter, absorption
  int event_nuclide_;  //!< index in nuclides array
  int event_mt_;       //!< reaction MT
  int delayed_group_ {0};  //!< delayed group

  // Post-collision physical data
  int n_bank_ {0};        //!< number of fission sites banked
  int n_bank_second_ {0}; //!< number of secondary particles banked
  double wgt_bank_ {0.0}; //!< weight of fission sites banked
  int n_delayed_bank_[MAX_DELAYED_GROUPS];  //!< number of delayed fission
                                            //!< sites banked

  // Indices for various arrays
  int surface_ {0};         //!< index for surface particle is on
  int cell_born_ {-1};      //!< index for cell particle was born in
  int material_ {-1};       //!< index for current material
  int material_last_ {-1};  //!< index for last material

  // Boundary information
  BoundaryInfo boundary_;

  // Temperature of current cell
  double sqrtkT_ {-1.0};      //!< sqrt(k_Boltzmann * temperature) in eV
  double sqrtkT_last_ {0.0};  //!< last temperature

  // Statistical data
  int n_collision_ {0};  //!< number of collisions

  // Track output
  bool write_track_ {false};

  // Current PRNG state
  uint64_t seeds_[N_STREAMS]; // current seeds
  int      stream_;           // current RNG stream

  // Secondary particle bank
  // This one is necessarilly going to be large/wasteful, unless we add in shared
  // secondary bank among threads or something along those lines
  //std::vector<Particle::Bank> secondary_bank_;
  Particle::Bank secondary_bank_[SECONDARY_BANK_SIZE];
  uint64_t secondary_bank_length_ {0};

  int64_t current_work_; // current work index

  // TODO: flux_derivs can eventually be converted to an allocated array, with size fixed at runtime
  //std::vector<double> flux_derivs_;  // for derivatives for this particle
  double flux_derivs_[FLUX_DERIVS_SIZE];  // for derivatives for this particle

  // TODO: filter_matches_ can eventually be converted to an allocated array, with size fixed at runtime
  //std::vector<FilterMatch> filter_matches_; // tally filter matches
  FilterMatch* filter_matches_; // tally filter matches

  std::vector<std::vector<Position>> tracks_; // tracks for outputting to file

  //std::vector<NuBank> nu_bank_; // bank of most recently fissioned particles
  NuBank nu_bank_[NU_BANK_SIZE]; // bank of most recently fissioned particles

  // Global tally accumulators
  double keff_tally_absorption_ {0.0};
  double keff_tally_collision_ {0.0};
  double keff_tally_tracklength_ {0.0};
  double keff_tally_leakage_ {0.0};

  bool trace_ {false};     //!< flag to show debug information

  double collision_distance_; // distance to particle's next closest collision

  double advance_distance_; // distance the particle actually advanced this event

  int n_event_ {0}; // number of events executed in this particle's history

  // DagMC state variables
  #ifdef DAGMC
  moab::DagMC::RayHistory history_;
  Direction last_dir_;
  #endif

  int64_t n_progeny_ {0}; // Number of progeny produced by this particle

  bool operator<(const Particle& rhs) const
  {
    if( alive() && !rhs.alive() )
      return true;
    if( !alive() && rhs.alive() )
      return false;
    if( !alive() && !rhs.alive() )
      return false;
    // case: both are alive, so we sort by id
    return id_ < rhs.id_;
  }
};

//============================================================================
//! Functions
//============================================================================

std::string particle_type_to_str(Particle::Type type);

Particle::Type str_to_particle_type(std::string str);

} // namespace openmc

#endif // OPENMC_PARTICLE_H
