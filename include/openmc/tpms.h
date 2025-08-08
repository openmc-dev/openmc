#ifndef OPENMC_TPMS_H
#define OPENMC_TPMS_H

#include <algorithm>
#include <numeric>

// #include "boost/math/tools/roots.hpp"
#include "openmc/bisect.h"
#include "openmc/position.h"
#include "openmc/vector.h"

#include "tpms_functions.h"

namespace openmc {

//==============================================================================
//! Generic Triply Periodic Minimal Surface (TPMS) class to solve that contains
//! the ray-tracing algorithm.
//==============================================================================

class TPMS {
public:
  struct rootFinding {
    bool isRoot;
    double xa;
    double xb;
    int status;
  };

public:
  TPMS(double _x0, double _y0, double _z0, double _a, double _b, double _c,
    double _d, double _e, double _f, double _g, double _h, double _i);
  virtual double fk(double k, Position r, Direction u) const { return 0.; };
  virtual double fpk(double k, Position r, Direction u) const { return 0.; };
  virtual double fppk(double k, Position r, Direction u) const { return 0.; };
  virtual double sampling_frequency(Direction u) const { return 1.; };
  rootFinding root_in_interval(double L0, double L1, Position r, Direction u);
  double ray_tracing(Position r, Direction u, double max_range);

public:
  double x0, y0, z0;
  double a, b, c, d, e, f, g, h, i;
};

//==============================================================================
//! TPMS class to define TPMS with constant isovalue and pitch. Simpler and
//! Faster has no interpolation of pitch and isovalue are needed.
//==============================================================================

class TPMSClassic : public TPMS {
public:
  TPMSClassic(double _isovalue, double _pitch, double _x0, double _y0,
    double _z0, double _a, double _b, double _c, double _d, double _e,
    double _f, double _g, double _h, double _i);

public:
  double isovalue, pitch;
};

//==============================================================================
//! TPMS class to define TPMS with isovalue and pitch defined as a function of
//! x,y and z. Useful to design heterogeneous lattice.
//==============================================================================

class TPMSFunction : public TPMS {
public:
  TPMSFunction(const FunctionForTPMS& _fIsovalue,
    const FunctionForTPMS& _fPitch, double _x0, double _y0, double _z0,
    double _a, double _b, double _c, double _d, double _e, double _f, double _g,
    double _h, double _i);
  double get_pitch(Position r) const;
  double get_isovalue(Position r) const;
  std::vector<double> get_pitch_first_partial_derivatives(Position r) const;
  std::vector<double> get_isovalue_first_partial_derivatives(Position r) const;
  std::vector<std::vector<double>> get_pitch_second_partial_derivatives(
    Position r) const;
  std::vector<std::vector<double>> get_isovalue_second_partial_derivatives(
    Position r) const;

public:
  const FunctionForTPMS& fIsovalue;
  const FunctionForTPMS& fPitch;
};

//==============================================================================
//! TPMS class with Schwarz-P equations implemented with constant pitch and
//! isovalue.
//==============================================================================

class SchwarzP : public TPMSClassic {
  using TPMSClassic::TPMSClassic;

public:
  double evaluate(Position r) const;
  Direction normal(Position r) const;
  double fk(double k, Position r, Direction u) const;
  double fpk(double k, Position r, Direction u) const;
  double fppk(double k, Position r, Direction u) const;
  double sampling_frequency(Direction u) const;
};

//==============================================================================
//! TPMS class with Gyroid equations implemented with constant pitch and
//! isovalue.
//==============================================================================

class Gyroid : public TPMSClassic {
  using TPMSClassic::TPMSClassic;

public:
  double evaluate(Position r) const;
  Direction normal(Position r) const;
  double fk(double k, Position r, Direction u) const;
  double fpk(double k, Position r, Direction u) const;
  double fppk(double k, Position r, Direction u) const;
  double sampling_frequency(Direction u) const;
};

//==============================================================================
//! TPMS class with Diamond equations implemented with constant pitch and
//! isovalue.
//==============================================================================

class Diamond : public TPMSClassic {
  using TPMSClassic::TPMSClassic;

public:
  double evaluate(Position r) const;
  Direction normal(Position r) const;
  double fk(double k, Position r, Direction u) const;
  double fpk(double k, Position r, Direction u) const;
  double fppk(double k, Position r, Direction u) const;
  double sampling_frequency(Direction u) const;
};

//==============================================================================
//! TPMS class with Schwarz-P equations implemented with functionalized pitch
//! and isovalue.
//==============================================================================

class FunctionSchwarzP : public TPMSFunction {
  using TPMSFunction::TPMSFunction;

public:
  double evaluate(Position r) const;
  Direction normal(Position r) const;
  double fk(double k, Position r, Direction u) const;
  double fpk(double k, Position r, Direction u) const;
  double fppk(double k, Position r, Direction u) const;
  double sampling_frequency(Direction u) const;
};

//==============================================================================
//! TPMS class with Gyroid equations implemented with functionalized pitch
//! and isovalue.
//==============================================================================

class FunctionGyroid : public TPMSFunction {
  using TPMSFunction::TPMSFunction;

public:
  double evaluate(Position r) const;
  Direction normal(Position r) const;
  double fk(double k, Position r, Direction u) const;
  double fpk(double k, Position r, Direction u) const;
  double fppk(double k, Position r, Direction u) const;
  double sampling_frequency(Direction u) const;
};

//==============================================================================
//! TPMS class with Diamond equations implemented with functionalized pitch
//! and isovalue.
//==============================================================================

class FunctionDiamond : public TPMSFunction {
  using TPMSFunction::TPMSFunction;

public:
  double evaluate(Position r) const;
  Direction normal(Position r) const;
  double fk(double k, Position r, Direction u) const;
  double fpk(double k, Position r, Direction u) const;
  double fppk(double k, Position r, Direction u) const;
  double sampling_frequency(Direction u) const;
};

} // namespace openmc
#endif