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

class TPMSClassic : public TPMS {
public:
  TPMSClassic(double _isovalue, double _pitch, double _x0, double _y0,
    double _z0, double _a, double _b, double _c, double _d, double _e,
    double _f, double _g, double _h, double _i);

public:
  double isovalue, pitch;
};

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