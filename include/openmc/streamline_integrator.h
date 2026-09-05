#ifndef OPENMC_STREAMLINE_INTEGRATOR_H
#define OPENMC_STREAMLINE_INTEGRATOR_H

#include "openmc/field.h"
#include "openmc/particle_data.h"

#include <string>

namespace openmc {

//! Streamline integrator
class StreamlineIntegrator {
public:
  virtual ~StreamlineIntegrator() = default;

  //! Advance to the next integration step.
  //!
  //! \param [inout] dt Time step
  //! \param [inout] y_n Position
  //! \param [in] cell_n Cell containing y_n
  //! \param [in] field Pointer to the velocity field
  virtual void next_step(
    double& t_n, Position& y_n, int cell_n, VelocityField* field) = 0;

  //! Time step accessors
  double& dt() { return dt_; }
  const double& dt() const { return dt_; }

private:
  double dt_; //!< Time step
};

//! Runge Kutta 4 integrator
class RK4StreamlineIntegrator : public StreamlineIntegrator {
public:
  RK4StreamlineIntegrator(double time_step) { dt() = time_step; }

  ~RK4StreamlineIntegrator() override = default;

  void next_step(
    double& t_n, Position& y_n, int cell_n, VelocityField* field) override;

private:
  std::string method_ = "Runge Kutta 4"; //!< Method name
};

} // namespace openmc

#endif // OPENMC_STREAMLINE_INTEGRATOR_H
