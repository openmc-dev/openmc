#ifndef OPENMC_PARTICLE_H
#define OPENMC_PARTICLE_H

//! \file particle.h
//! \brief Particle type

#include <cstdint>
#include <array>

#include "openmc/capi.h"

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

// Maximum number of secondary particles created
constexpr int MAX_SECONDARY {1000};

// Maximum number of lost particles
constexpr int MAX_LOST_PARTICLES {10};

// Maximum number of lost particles, relative to the total number of particles
constexpr double REL_MAX_LOST_PARTICLES {1.0e-6};

//! Particle types
enum class ParticleType {
  neutron, photon, electron, positron
};

extern "C" {

  struct LocalCoord {
    int cell {-1};
    int universe {-1};
    int lattice {-1};
    int lattice_x {-1};
    int lattice_y {-1};
    int lattice_z {-1};
    double xyz[3]; //!< particle position
    double uvw[3]; //!< particle direction
    bool rotated {false};  //!< Is the level rotated?

    //! clear data from a single coordinate level
    void reset();
  };

  //============================================================================
  //! State of a particle being transported through geometry
  //============================================================================

  struct Particle {
    int64_t id;  //!< Unique ID
    int type;    //!< Particle type (n, p, e, etc.)

    int n_coord;                  //!< number of current coordinate levels
    int cell_instance;            //!< offset for distributed properties
    LocalCoord coord[MAX_COORD];  //!< coordinates for all levels

    // Particle coordinates before crossing a surface
    int last_n_coord;          //!< number of current coordinates
    int last_cell[MAX_COORD];  //!< coordinates for all levels

    // Energy data
    double E;       //!< post-collision energy in eV
    double last_E;  //!< pre-collision energy in eV
    int g;          //!< post-collision energy group (MG only)
    int last_g;     //!< pre-collision energy group (MG only)

    // Other physical data
    double wgt;     //!< particle weight
    double mu;      //!< angle of scatter
    bool alive;     //!< is particle alive?

    // Other physical data
    double last_xyz_current[3];  //!< coordinates of the last collision or
                                 //!< reflective/periodic surface crossing for
                                 //!< current tallies
    double last_xyz[3];          //!< previous coordinates
    double last_uvw[3];          //!< previous direction coordinates
    double last_wgt;             //!< pre-collision particle weight
    double absorb_wgt;           //!< weight absorbed for survival biasing

    // What event took place
    bool fission;       //!< did particle cause implicit fission
    int event;          //!< scatter, absorption
    int event_nuclide;  //!< index in nuclides array
    int event_MT;       //!< reaction MT
    int delayed_group;  //!< delayed group

    // Post-collision physical data
    int n_bank;        //!< number of fission sites banked
    double wgt_bank;   //!< weight of fission sites banked
    int n_delayed_bank[MAX_DELAYED_GROUPS];  //!< number of delayed fission
                                             //!< sites banked

    // Indices for various arrays
    int surface;        //!< index for surface particle is on
    int cell_born;      //!< index for cell particle was born in
    int material;       //!< index for current material
    int last_material;  //!< index for last material

    // Temperature of current cell
    double sqrtkT;       //!< sqrt(k_Boltzmann * temperature) in eV
    double last_sqrtkT;  //!< last temperature

    // Statistical data
    int n_collision;  //!< number of collisions

    // Track output
    bool write_track {false};

    // Secondary particles created
    int64_t n_secondary {};
    Bank secondary_bank[MAX_SECONDARY];

    //! resets all coordinate levels for the particle
    void clear();

    //! create a secondary particle
    //
    //! stores the current phase space attributes of the particle in the
    //! secondary bank and increments the number of sites in the secondary bank.
    //! \param uvw Direction of the secondary particle
    //! \param E Energy of the secondary particle in [eV]
    //! \param type Particle type
    //! \param run_CE Whether continuous-energy data is being used
    void create_secondary(const double* uvw, double E, int type, bool run_CE);

    //! sets default attributes for a particle
    void initialize();

    //! initialize from a source site
    //
    //! initializes a particle from data stored in a source site. The source
    //! site may have been produced from an external source, from fission, or
    //! simply as a secondary particle.
    //! \param src Source site data
    //! \param run_CE Whether continuous-energy data is being used
    //! \param energy_bin_avg An array of energy group bin averages
    void from_source(const Bank* src, bool run_CE, const double* energy_bin_avg);

    //! mark a particle as lost and create a particle restart file
    //! \param message A warning message to display
    void mark_as_lost(const char* message);

    //! create a particle restart HDF5 file
    void write_restart();
  };


  //============================================================================
  // Fortran compatibility functions
  //============================================================================

  void reset_coord(LocalCoord* c);
  void particle_clear(Particle* p);
  void particle_create_secondary(Particle* p, const double* uvw, double E,
                                 int type, bool run_CE);
  void particle_initialize(Particle* p);
  void particle_from_source(Particle* p, const Bank* src, bool run_CE,
                            const double* energy_bin_avg);
  void particle_mark_as_lost(Particle* p, const char* message);
  void particle_write_restart(Particle* p);

} // extern "C"

} // namespace openmc

#endif // OPENMC_PARTICLE_H
