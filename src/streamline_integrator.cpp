#include "openmc/streamline_integrator.h"
#include "openmc/field.h"

namespace openmc {

void RK4StreamlineIntegrator::next_step(
  double& tn, Position& yn, VelocityField& field)
{
  int datacell_id;

  // Calculate k1
  Direction k1 = field.evaluate(yn) * dt();

  // Calculate k2
  Position p2 = yn + k1 / 2.;
  Direction k2 = field.evaluate(p2) * dt();

  // Calculate k3
  Position p3 = yn + k2 / 2.;
  Direction k3 = field.evaluate(p3) * dt();

  // Calculate k4
  Position p4 = yn + k3;
  Direction k4 = field.evaluate(p4) * dt();

  // Step forward
  yn += (k1 + k2 * 2 + k3 * 2 + k4) / 6.;
  tn += dt();
}

} // namespace openmc
