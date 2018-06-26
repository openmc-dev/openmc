#ifndef GEOMETRY_H
#define GEOMETRY_H

namespace openmc {

extern "C" int openmc_root_universe;

struct Position {
  double x = 0.;
  double y = 0.;
  double z = 0.;

  Position() = default;
  Position(double x_, double y_, double z_) : x{x_}, y{y_}, z{z_} { };
  Position(const double xyz[]) : x{xyz[0]}, y{xyz[1]}, z{xyz[2]} { };

  Position& operator+=(Position);
  Position& operator+=(double);
  Position& operator-=(Position);
  Position& operator-=(double);
  Position& operator*=(Position);
  Position& operator*=(double);
  const double& operator[](int i) const {
    switch (i) {
      case 0: return x;
      case 1: return y;
      case 2: return z;
    }
  }
  double& operator[](int i) {
    switch (i) {
      case 0: return x;
      case 1: return y;
      case 2: return z;
    }
  }

  inline double dot(Position other) {
    return x*other.x + y*other.y + z*other.z;
  }
};

inline Position operator+(Position a, Position b) { return a += b; }
inline Position operator+(Position a, double b)   { return a += b; }
inline Position operator+(double a, Position b)   { return b += a; }

inline Position operator-(Position a, Position b) { return a -= b; }
inline Position operator-(Position a, double b)   { return a -= b; }
inline Position operator-(double a, Position b)   { return b -= a; }

inline Position operator*(Position a, Position b) { return a *= b; }
inline Position operator*(Position a, double b)   { return a *= b; }
inline Position operator*(double a, Position b)   { return b *= a; }


struct Angle : Position {
  double& u() { return x; }
  double& v() { return y; }
  double& w() { return z; }
  Angle() = default;
  Angle(double u, double v, double w) : Position{u, v, w} { };
  Angle(const double uvw[]) : Position{uvw} { };
  Angle(Position r) : Position{r} { };
};

} // namespace openmc

#endif // GEOMETRY_H
