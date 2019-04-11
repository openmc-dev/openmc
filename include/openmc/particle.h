#ifndef OPENMC_PARTICLE_H
#define OPENMC_PARTICLE_H

//! \file particle.h
//! \brief Particle type

#include <array>
#include <cstdint>
#include <memory> // for unique_ptr
#include <sstream>
#include <string>

#include "openmc/constants.h"
#include "openmc/position.h"

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

// Maximum number of lost particles
constexpr int MAX_LOST_PARTICLES {10};

// Maximum number of lost particles, relative to the total number of particles
constexpr double REL_MAX_LOST_PARTICLES {1.0e-6};

constexpr double CACHE_INVALID {-1.0};

//==============================================================================
// Class declarations
//==============================================================================

struct LocalCoord {
  Position r; //!< particle position
  Direction u; //!< particle direction
  int cell {-1};
  int universe {-1};
  int lattice {-1};
  int lattice_x {-1};
  int lattice_y {-1};
  int lattice_z {-1};
  bool rotated {false};  //!< Is the level rotated?

  //! clear data from a single coordinate level
  void reset();
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
  double reaction[DEPLETION_RX.size()];

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

  // Photon cross sections
  double coherent;        //!< macroscopic coherent xs
  double incoherent;      //!< macroscopic incoherent xs
  double photoelectric;   //!< macroscopic photoelectric xs
  double pair_production; //!< macroscopic pair production xs
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
  struct Bank {
    Position r;
    Direction u;
    double E;
    double wgt;
    int delayed_group;
    Type particle;
  };

  //==========================================================================
  // Constructors

  Particle();

  //==========================================================================
  // Methods and accessors

  // Accessors for position in global coordinates
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

  //! resets all coordinate levels for the particle
  void clear();

  //! create a secondary particle
  //
  //! stores the current phase space attributes of the particle in the
  //! secondary bank and increments the number of sites in the secondary bank.
  //! \param u Direction of the secondary particle
  //! \param E Energy of the secondary particle in [eV]
  //! \param type Particle type
  void create_secondary(Direction u, double E, Type type) const;

  //! initialize from a source site
  //
  //! initializes a particle from data stored in a source site. The source
  //! site may have been produced from an external source, from fission, or
  //! simply as a secondary particle.
  //! \param src Source site data
  void from_source(const Bank* src);

  //! Transport a particle from birth to death
  void transport();

  //! Cross a surface and handle boundary conditions
  void cross_surface();

  //! mark a particle as lost and create a particle restart file
  //! \param message A warning message to display
  void mark_as_lost(const char* message);

  void mark_as_lost(const std::string& message)
  {mark_as_lost(message.c_str());}

  void mark_as_lost(const std::stringstream& message)
  {mark_as_lost(message.str());}

  //! create a particle restart HDF5 file
  void write_restart() const;

  //==========================================================================
  // Data members

  // Cross section caches
  std::vector<NuclideMicroXS> neutron_xs_; //!< Microscopic neutron cross sections
  std::vector<ElementMicroXS> photon_xs_; //!< Microscopic photon cross sections
  MacroXS macro_xs_; //!< Macroscopic cross sections

  int64_t id_;  //!< Unique ID
  Type type_ {Type::neutron};   //!< Particle type (n, p, e, etc.)

  int n_coord_ {1};              //!< number of current coordinate levels
  int cell_instance_;            //!< offset for distributed properties
  std::vector<LocalCoord> coord_; //!< coordinates for all levels

  // Particle coordinates before crossing a surface
  int n_coord_last_ {1};      //!< number of current coordinates
  std::vector<int> cell_last_;  //!< coordinates for all levels

  // Energy data
  double E_;       //!< post-collision energy in eV
  double E_last_;  //!< pre-collision energy in eV
  int g_ {0};      //!< post-collision energy group (MG only)
  int g_last_;     //!< pre-collision energy group (MG only)

  // Other physical data
  double wgt_ {1.0};     //!< particle weight
  double mu_;      //!< angle of scatter
  bool alive_ {true};     //!< is particle alive?

  // Other physical data
  Position r_last_current_; //!< coordinates of the last collision or
                            //!< reflective/periodic surface crossing for
                            //!< current tallies
  Position r_last_;   //!< previous coordinates
  Direction u_last_;  //!< previous direction coordinates
  double wgt_last_ {1.0};   //!< pre-collision particle weight
  double wgt_absorb_ {0.0}; //!< weight absorbed for survival biasing

  // What event took place
  bool fission_ {false}; //!< did particle cause implicit fission
  int event_;          //!< scatter, absorption
  int event_nuclide_;  //!< index in nuclides array
  int event_mt_;       //!< reaction MT
  int delayed_group_ {0};  //!< delayed group

  // Post-collision physical data
  int n_bank_ {0};        //!< number of fission sites banked
  double wgt_bank_ {0.0}; //!< weight of fission sites banked
  int n_delayed_bank_[MAX_DELAYED_GROUPS];  //!< number of delayed fission
                                            //!< sites banked

  // Indices for various arrays
  int surface_ {0};             //!< index for surface particle is on
  int cell_born_ {-1};      //!< index for cell particle was born in
  int material_ {-1};       //!< index for current material
  int material_last_ {-1};  //!< index for last material

  // Temperature of current cell
  double sqrtkT_ {-1.0};      //!< sqrt(k_Boltzmann * temperature) in eV
  double sqrtkT_last_ {0.0};  //!< last temperature

  // Statistical data
  int n_collision_ {0};  //!< number of collisions

  // Track output
  bool write_track_ {false};
};

} // namespace openmc

#endif // OPENMC_PARTICLE_H
