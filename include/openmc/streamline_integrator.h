#ifndef OPENMC_STREAMLINE_INTEGRATOR_H
#define OPENMC_STREAMLINE_INTEGRATOR_H

#include "openmc/error.h"
#include "openmc/field.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"

#include <string>

namespace openmc {

//! Abstract integrator
class StreamlineIntegrator {
public:
  //! Advance to the next integration step
  //!
  //! \param [inout] dt Time [s]
  //! \param [inout] y_n Position
  //! \param [in] cell_n Cell containing y_n
  //! \param [in] field Pointer to the velocity field
  virtual void next_step(
    double& t_n, Position& y_n, int cell_n, VelocityField& field) = 0;

  //! Time step accessors
  double& dt() { return dt_; }
  const double& dt() const { return dt_; }

private:
  double dt_; //!< Time step [s]
};

//! Runge Kutta 4 integrator
class RK4StreamlineIntegrator : public StreamlineIntegrator {
public:
  //! Constructor
  RK4StreamlineIntegrator(double time_step)
  {
    dt() = time_step;
  }

  void next_step(
    double& t_n, Position& y_n, int cell_n, VelocityField& field) override;

private:
  std::string method_ = "Runge Kutta 4"; //!< Method name
};

} // namespace openmc

#endif // OPENMC_STREAMLINE_INTEGRATOR_H
