#include "openmc/streamline_integrator.h"
#include "openmc/field.h"
#include "openmc/particle_data.h"

namespace openmc {

void RK4StreamlineIntegrator::next_step(
  double& tn, Position& yn, int cell_n, VelocityField* field)
{
  // Intermediate velocity estimates
  Direction k1 = field->evaluate_in_mesh(yn, cell_n);
  Direction k2 = field->evaluate_clamped(yn, yn + dt() / 2. * k1, cell_n);
  Direction k3 = field->evaluate_clamped(yn, yn + dt() / 2. * k2, cell_n);
  Direction k4 = field->evaluate_clamped(yn, yn + dt() * k3, cell_n);

  // Step forward
  yn += dt() / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
  tn += dt();
}

} // namespace openmc
