#ifndef OPENMC_TPMS_H
#define OPENMC_TPMS_H

// #include "boost/math/tools/roots.hpp"
#include "openmc/bisect.h"

#include "openmc/vector.h"
#include "openmc/position.h"

namespace openmc {

class TPMS
{
public:
    struct rootFinding {bool isRoot; double xa; double xb; int status;};
public:
    TPMS(double _cst, double _pitch, double _x0, double _y0, double _z0, double _a, double _b, double _c, double _d, double _e, double _f, double _g, double _h, double _i);
    virtual double fk(double k, Position r, Direction u) const {return 0.;};
    virtual double fpk(double k, Position r, Direction u) const {return 0.;};
    virtual double fppk(double k, Position r, Direction u) const {return 0.;};
    virtual double sampling_frequency(Direction u) const {return 1.;};
    rootFinding root_in_interval(double L0, double L1, Position r, Direction u);
    double ray_tracing(Position r, Direction u, double max_range);

public:
    double cst, pitch;
    double x0, y0, z0;
    double a, b, c, d, e, f, g, h, i;

public:
    double XLIM = 1.0e6;
};

class SchwarzP : public TPMS
{
using TPMS::TPMS;
public:
    double evaluate(Position r) const;
    Direction normal(Position r) const;
    double fk(double k, Position r, Direction u) const;
    double fpk(double k, Position r, Direction u) const;
    double fppk(double k, Position r, Direction u) const;
    double sampling_frequency(Direction u) const;
};

class Gyroid : public TPMS
{
using TPMS::TPMS;
public:
    double evaluate(Position r) const;
    Direction normal(Position r) const;
    double fk(double k, Position r, Direction u) const;
    double fpk(double k, Position r, Direction u) const;
    double fppk(double k, Position r, Direction u) const;
    double sampling_frequency(Direction u) const;
};

class Diamond : public TPMS
{
using TPMS::TPMS;
public:
    double evaluate(Position r) const;
    Direction normal(Position r) const;
    double fk(double k, Position r, Direction u) const;
    double fpk(double k, Position r, Direction u) const;
    double fppk(double k, Position r, Direction u) const;
    double sampling_frequency(Direction u) const;
};

}
#endif