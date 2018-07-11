#ifndef OPENMC_ENDF_H
#define OPENMC_ENDF_H

#include <vector>

#include "constants.h"
#include "hdf5.h"

namespace openmc {

Interpolation int2interp(int i);
bool is_fission(int MT);

class Function1D {
public:
  virtual double operator()(double x) const = 0;
};

class Polynomial : public Function1D {
public:
  explicit Polynomial(hid_t dset);
  double operator()(double x) const;
private:
  std::vector<double> coef_;
};

class Tabulated1D : public Function1D {
public:
  Tabulated1D() = default;
  explicit Tabulated1D(hid_t dset);
  double operator()(double x) const;
private:
  std::size_t n_regions_ {0};
  std::vector<int> nbt_;
  std::vector<Interpolation> int_;
  std::size_t n_pairs_;
  std::vector<double> x_;
  std::vector<double> y_;
};

} // namespace openmc

#endif // OPENMC_ENDF_H
