#include "openmc/physics_mg.h"

#include "openmc/constants.h"
#include "openmc/math_functions.h"
#include "openmc/mgxs_interface.h"

namespace openmc {

void scatter(Particle* p, const double* energy_bin_avg)
{
  // Adjust indices for Fortran to C++ indexing
  // TODO: Remove when no longer needed
  int gin = p->last_g - 1;
  int gout = p->g - 1;
  int i_mat = p->material - 1;
  macro_xs[i_mat].sample_scatter(gin, gout, p->mu, p->wgt);

  // Adjust return value for fortran indexing
  // TODO: Remove when no longer needed
  p->g = gout + 1;

  // Rotate the angle
  rotate_angle_c(p->coord[0].uvw, p->mu, nullptr);

  // Update energy value for downstream compatability (in tallying)
  p->E = energy_bin_avg[gout];

  // Set event component
  p->event = EVENT_SCATTER;
}

} //namespace openmc
