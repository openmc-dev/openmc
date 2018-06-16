//! \file scattdata.h
//! A collection of multi-group scattering data classes

#ifndef SCATTDATA_H
#define SCATTDATA_H

#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>

#include "constants.h"
#include "math_functions.h"
#include "random_lcg.h"
#include "error.h"

namespace openmc {

// temporary declaations so we can name our friend functions
class ScattDataLegendre;
class ScattDataTabular;

//==============================================================================
// SCATTDATA contains all the data needed to describe the scattering energy and
// angular distribution data
//==============================================================================

class ScattData {
  protected:
    void generic_init(int order, int_1dvec& in_gmin, int_1dvec& in_gmax,
                      double_2dvec& in_energy, double_2dvec& in_mult);
    void generic_combine(const int max_order,
                         const std::vector<ScattData*>& those_scatts,
                         const double_1dvec& scalars, int_1dvec& in_gmin,
                         int_1dvec& in_gmax, double_2dvec& sparse_mult,
                         double_3dvec& sparse_scatter);
  public:
    double_2dvec energy; // Normalized p0 matrix for sampling Eout
    double_2dvec mult;   // nu-scatter multiplication (nu-scatt/scatt)
    double_3dvec dist;   // Angular distribution
    int_1dvec gmin;      // minimum outgoing group
    int_1dvec gmax;      // maximum outgoing group
    double_1dvec scattxs; // Isotropic Sigma_{s,g_{in}}
    virtual double calc_f(const int gin, const int gout, const double mu) = 0;
    virtual void sample(const int gin, int& gout, double& mu, double& wgt) = 0;
    virtual void init(int_1dvec& in_gmin, int_1dvec& in_gmax,
                      double_2dvec& in_mult, double_3dvec& coeffs) = 0;
    void sample_energy(int gin, int& gout, int& i_gout);
    double get_xs(const int xstype, const int gin, const int* gout,
                  const double* mu);
    virtual void combine(const std::vector<ScattData*>& those_scatts,
                         const double_1dvec& scalars) = 0;
    virtual int get_order() = 0;
    virtual double_3dvec get_matrix(const int max_order) = 0;
};

//==============================================================================
// ScattDataLegendre represents the angular distributions as Legendre kernels
//==============================================================================

class ScattDataLegendre: public ScattData {
  protected:
    // Maximal value for rejection sampling from a rectangle
    double_2dvec max_val;
    friend void convert_legendre_to_tabular(ScattDataLegendre& leg,
                                            ScattDataTabular& tab, int n_mu);
  public:
    void init(int_1dvec& in_gmin, int_1dvec& in_gmax, double_2dvec& in_mult,
              double_3dvec& coeffs);
    void update_max_val();
    double calc_f(const int gin, const int gout, const double mu);
    void sample(const int gin, int& gout, double& mu, double& wgt);
    void combine(const std::vector<ScattData*>& those_scatts,
                 const double_1dvec& scalars);
    int get_order() {return dist[0][0].size() - 1;};
    double_3dvec get_matrix(const int max_order);
};

//==============================================================================
// ScattDataHistogram represents the angular distributions as a histogram, as it
// would be if it came from a "mu" tally in OpenMC
//==============================================================================

class ScattDataHistogram: public ScattData {
  protected:
    double_1dvec mu;
    double dmu;
    double_3dvec fmu;
  public:
    void init(int_1dvec& in_gmin, int_1dvec& in_gmax, double_2dvec& in_mult,
              double_3dvec& coeffs);
    double calc_f(const int gin, const int gout, const double mu);
    void sample(const int gin, int& gout, double& mu, double& wgt);
    void combine(const std::vector<ScattData*>& those_scatts,
                 const double_1dvec& scalars);
    int get_order() {return dist[0][0].size();};
    double_3dvec get_matrix(const int max_order);
};

//==============================================================================
// ScattDataTabular represents the angular distributions as a table of mu and
// f(mu)
//==============================================================================

class ScattDataTabular: public ScattData {
  protected:
    double_1dvec mu;
    double dmu;
    double_3dvec fmu;
    friend void convert_legendre_to_tabular(ScattDataLegendre& leg,
                                            ScattDataTabular& tab, int n_mu);
  public:
    void init(int_1dvec& in_gmin, int_1dvec& in_gmax, double_2dvec& in_mult,
              double_3dvec& coeffs);
    double calc_f(const int gin, const int gout, const double mu);
    void sample(const int gin, int& gout, double& mu, double& wgt);
    void combine(const std::vector<ScattData*>& those_scatts,
                 const double_1dvec& scalars);
    int get_order() {return dist[0][0].size();};
    double_3dvec get_matrix(const int max_order);
};

//==============================================================================
// Function to convert Legendre functions to tabular
//==============================================================================

void
convert_legendre_to_tabular(ScattDataLegendre& leg, ScattDataTabular& tab,
                            int n_mu);

} // namespace openmc
#endif // SCATTDATA_H