#ifndef OPENMC_PARTICLE_H
#define OPENMC_PARTICLE_H

//! \file particle.h
//! \brief Particle type

#include <cstdint>
#include <string>

#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle_data.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"
#include "openmc/tallies/filter_match.h"
#include "openmc/vector.h"

namespace openmc {

// Forward declare the Surface class for use in Particle::cross_vacuum_bc, etc.
class Surface;

/*
 * The Particle class encompasses data and methods for transporting particles
 * through their lifecycle. Its base class defines particle data layout in
 * memory. A more detailed description of the rationale behind this approach
 * can be found in particle_data.h.
 */

class Particle : public ParticleData {
public:
  //==========================================================================
  // Constructors

  Particle() = default;

  //==========================================================================
  // Methods

  double speed() const;

  //! moves the particle by the distance length to its next location
  //! \param length the distance the particle is moved
  void move_distance(double length);

  //! create a secondary particle
  //
  //! stores the current phase space attributes of the particle in the
  //! secondary bank and increments the number of sites in the secondary bank.
  //! \param wgt Weight of the secondary particle
  //! \param u Direction of the secondary particle
  //! \param E Energy of the secondary particle in [eV]
  //! \param type Particle type
  //! \return Whether a secondary particle was created
  bool create_secondary(double wgt, Direction u, double E, ParticleType type);

  //! initialize from a source site
  //
  //! initializes a particle from data stored in a source site. The source
  //! site may have been produced from an external source, from fission, or
  //! simply as a secondary particle.
  //! \param src Source site data
  void from_source(const SourceSite* src);

  // Coarse-grained particle events
  void event_calculate_xs();
  void event_advance();
  void event_cross_surface();
  void event_collide();
  void event_revive_from_secondary();
  void event_death();

  //! pulse-height recording
  void pht_collision_energy();
  void pht_secondary_particles();

  //! Cross a surface and handle boundary conditions
  void cross_surface(const Surface& surf);

  //! Cross a vacuum boundary condition.
  //
  //! \param surf The surface (with the vacuum boundary condition) that the
  //!   particle struck.
  void cross_vacuum_bc(const Surface& surf);

  //! Cross a reflective boundary condition.
  //
  //! \param surf The surface (with the reflective boundary condition) that the
  //!   particle struck.
  //! \param new_u The direction of the particle after reflection.
  void cross_reflective_bc(const Surface& surf, Direction new_u);

  //! Cross a periodic boundary condition.
  //
  //! \param surf The surface (with the periodic boundary condition) that the
  //!   particle struck.
  //! \param new_r The position of the particle after translation/rotation.
  //! \param new_u The direction of the particle after translation/rotation.
  //! \param new_surface The signed index of the surface that the particle will
  //!   reside on after translation/rotation.
  void cross_periodic_bc(
    const Surface& surf, Position new_r, Direction new_u, int new_surface);

  //! mark a particle as lost and create a particle restart file
  //! \param message A warning message to display
  virtual void mark_as_lost(const char* message) override;
  using GeometryState::mark_as_lost;

  //! create a particle restart HDF5 file
  void write_restart() const;

  //! Update microscopic cross section cache
  //
  //! \param[in] i_nuclide Index in data::nuclides
  //! \param[in] i_grid Index on log union grid
  //! \param[in] i_sab Index in data::thermal_scatt
  //! \param[in] sab_frac  S(a,b) table fraction
  //! \param[in] ncrystal_xs Thermal scattering xs from NCrystal
  void update_neutron_xs(int i_nuclide, int i_grid = C_NONE, int i_sab = C_NONE,
    double sab_frac = 0.0, double ncrystal_xs = -1.0);
};

//============================================================================
//! Functions
//============================================================================

std::string particle_type_to_str(ParticleType type);

ParticleType str_to_particle_type(std::string str);

void add_surf_source_to_bank(Particle& p, const Surface& surf);

} // namespace openmc

#endif // OPENMC_PARTICLE_H
