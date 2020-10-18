#include "openmc/boundary_condition.h"
#include "openmc/surface.h"

namespace openmc {

void
VacuumBC::handle_particle(Particle& p, const Surface& surf) const
{
  p.cross_vacuum_bc(surf);
}

void
ReflectiveBC::handle_particle(Particle& p, const Surface& surf) const
{
  Direction u = surf.reflect(p.r(), p.u(), &p);
  u /= u.norm();

  p.cross_reflective_bc(surf, u);
}

void
WhiteBC::handle_particle(Particle& p, const Surface& surf) const
{
  Direction u = surf.diffuse_reflect(p.r(), p.u(), p.current_seed());
  u /= u.norm();

  p.cross_reflective_bc(surf, u);
}

void
PeriodicBC::handle_particle(Particle& p, const Surface& surf) const
{
  // Get a pointer to the partner periodic surface
  auto surf_p = dynamic_cast<const PeriodicSurface*>(&surf);
  auto other = dynamic_cast<const PeriodicSurface*>(
    model::surfaces[surf_p->i_periodic_].get());

  // Compute the new particle location and direction
  Position r {p.r()};
  Direction u {p.u()};
  bool rotational = other->periodic_translate(surf_p, r, u);

  // Pick the particle's new surface
  // TODO: off-by-one
  int new_surface = rotational ?
    surf_p->i_periodic_ + 1 :
    ((p.surface_ > 0) ? surf_p->i_periodic_ + 1 : -(surf_p->i_periodic_ + 1));

  p.cross_periodic_bc(surf, r, u, new_surface);
}

} // namespace openmc
