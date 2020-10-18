#include "openmc/boundary_condition.h"
#include "openmc/surface.h"

namespace openmc {

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

//void
//PeriodicBC::handle_particle(Particle& p, const Surface& surf) const
//{
//}

} // namespace openmc
