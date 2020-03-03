//! \file scattdata.h
//! A collection of multi-group scattering data classes

#ifndef OPENMC_SCATTDATA_H
#define OPENMC_SCATTDATA_H

#include <vector>

#include "xtensor/xtensor.hpp"

#include "openmc/constants.h"

namespace openmc {

// forward declarations so we can name our friend functions
class ScattDataLegendre;
class ScattDataTabular;

//==============================================================================
// SCATTDATA contains all the data needed to describe the scattering energy and
// angular distribution data
//==============================================================================

class ScattData {
  public:
    virtual ~ScattData() = default;
  protected:
    //! \brief Initializes the attributes of the base class.
    void
    base_init(int order, const xt::xtensor<int, 1>& in_gmin,
         const xt::xtensor<int, 1>& in_gmax, const double_2dvec& in_energy,
         const double_2dvec& in_mult);

    //! \brief Combines microscopic ScattDatas into a macroscopic one.
    void
    base_combine(size_t max_order, const std::vector<ScattData*>& those_scatts,
         const std::vector<double>& scalars, xt::xtensor<int, 1>& in_gmin,
         xt::xtensor<int, 1>& in_gmax, double_2dvec& sparse_mult,
         double_3dvec& sparse_scatter);

  public:

    double_2dvec energy;            // Normalized p0 matrix for sampling Eout
    double_2dvec mult;              // nu-scatter multiplication (nu-scatt/scatt)
    double_3dvec dist;              // Angular distribution
    xt::xtensor<int, 1> gmin;    // minimum outgoing group
    xt::xtensor<int, 1> gmax;    // maximum outgoing group
    xt::xtensor<double, 1> scattxs; // Isotropic Sigma_{s,g_{in}}

    //! \brief Calculates the value of normalized f(mu).
    //!
    //! The value of f(mu) is normalized as in the integral of f(mu)dmu across
    //! [-1,1] is 1.
    //!
    //! @param gin Incoming energy group of interest.
    //! @param gout Outgoing energy group of interest.
    //! @param mu Cosine of the change-in-angle of interest.
    //! @return The value of f(mu).
    virtual double
    calc_f(int gin, int gout, double mu) = 0;

    //! \brief Samples the outgoing energy and angle from the ScattData info.
    //!
    //! @param gin Incoming energy group.
    //! @param gout Sampled outgoing energy group.
    //! @param mu Sampled cosine of the change-in-angle.
    //! @param wgt Weight of the particle to be adjusted.
    //! @param seed Pseudorandom number seed pointer
    virtual void
    sample(int gin, int& gout, double& mu, double& wgt, uint64_t* seed) = 0;

    //! \brief Initializes the ScattData object from a given scatter and
    //!   multiplicity matrix.
    //!
    //! @param in_gmin List of minimum outgoing groups for every incoming group
    //! @param in_gmax List of maximum outgoing groups for every incoming group
    //! @param in_mult Input sparse multiplicity matrix
    //! @param coeffs Input sparse scattering matrix
    virtual void
    init(const xt::xtensor<int, 1>& in_gmin, const xt::xtensor<int, 1>& in_gmax,
         const double_2dvec& in_mult, const double_3dvec& coeffs) = 0;

    //! \brief Combines the microscopic data.
    //!
    //! @param those_scatts Microscopic objects to combine.
    //! @param scalars Scalars to multiply the microscopic data by.
    virtual void
    combine(const std::vector<ScattData*>& those_scatts,
         const std::vector<double>& scalars) = 0;

    //! \brief Getter for the dimensionality of the scattering order.
    //!
    //! If Legendre this is the "n" in "Pn"; for Tabular, this is the number
    //! of points, and for Histogram this is the number of bins.
    //!
    //! @return The order.
    virtual size_t
    get_order() = 0;

    //! \brief Builds a dense scattering matrix from the constituent parts
    //!
    //! @param max_order If Legendre this is the maximum value of "n" in "Pn"
    //!   requested; ignored otherwise.
    //! @return The dense scattering matrix.
    virtual xt::xtensor<double, 3>
    get_matrix(size_t max_order) = 0;

    //! \brief Samples the outgoing energy from the ScattData info.
    //!
    //! @param gin Incoming energy group.
    //! @param gout Sampled outgoing energy group.
    //! @param i_gout Sampled outgoing energy group index.
    //! @param seed Pseudorandom number seed pointer
    void
    sample_energy(int gin, int& gout, int& i_gout, uint64_t* seed);

    //! \brief Provides a cross section value given certain parameters
    //!
    //! @param xstype Type of cross section requested, according to the
    //!   enumerated constants.
    //! @param gin Incoming energy group.
    //! @param gout Outgoing energy group; use nullptr if irrelevant, or if a
    //!   sum is requested.
    //! @param mu Cosine of the change-in-angle, for scattering quantities;
    //!   use nullptr if irrelevant.
    //! @return Requested cross section value.
    double
    get_xs(MgxsType xstype, int gin, const int* gout, const double* mu);
};

//==============================================================================
// ScattDataLegendre represents the angular distributions as Legendre kernels
//==============================================================================

class ScattDataLegendre: public ScattData {

  protected:

    // Maximal value for rejection sampling from a rectangle
    double_2dvec max_val;

    // Friend convert_legendre_to_tabular so it has access to protected
    // parameters
    friend void
    convert_legendre_to_tabular(ScattDataLegendre& leg, ScattDataTabular& tab);

  public:

    void
    init(const xt::xtensor<int, 1>& in_gmin, const xt::xtensor<int, 1>& in_gmax,
         const double_2dvec& in_mult, const double_3dvec& coeffs);

    void
    combine(const std::vector<ScattData*>& those_scatts,
            const std::vector<double>& scalars);

    //! \brief Find the maximal value of the angular distribution to use as a
    // bounding box with rejection sampling.
    void
    update_max_val();

    double
    calc_f(int gin, int gout, double mu);

    void
    sample(int gin, int& gout, double& mu, double& wgt, uint64_t* seed);

    size_t
    get_order() {return dist[0][0].size() - 1;};

    xt::xtensor<double, 3>
    get_matrix(size_t max_order);
};

//==============================================================================
// ScattDataHistogram represents the angular distributions as a histogram, as it
// would be if it came from a "mu" tally in OpenMC
//==============================================================================

class ScattDataHistogram: public ScattData {

  protected:

    xt::xtensor<double, 1> mu; // Angle distribution mu bin boundaries
    double dmu;                // Quick storage of the mu spacing
    double_3dvec fmu;          // The angular distribution histogram

  public:

    void
    init(const xt::xtensor<int, 1>& in_gmin, const xt::xtensor<int, 1>& in_gmax,
         const double_2dvec& in_mult, const double_3dvec& coeffs);

    void
    combine(const std::vector<ScattData*>& those_scatts,
            const std::vector<double>& scalars);

    double
    calc_f(int gin, int gout, double mu);

    void
    sample(int gin, int& gout, double& mu, double& wgt, uint64_t* seed);

    size_t
    get_order() {return dist[0][0].size();};

    xt::xtensor<double, 3>
    get_matrix(size_t max_order);
};

//==============================================================================
// ScattDataTabular represents the angular distributions as a table of mu and
// f(mu)
//==============================================================================

class ScattDataTabular: public ScattData {

  protected:

    xt::xtensor<double, 1> mu; // Angle distribution mu grid points
    double dmu;                // Quick storage of the mu spacing
    double_3dvec fmu;          // The angular distribution function

    // Friend convert_legendre_to_tabular so it has access to protected
    // parameters
    friend void
    convert_legendre_to_tabular(ScattDataLegendre& leg, ScattDataTabular& tab);

  public:

    void
    init(const xt::xtensor<int, 1>& in_gmin, const xt::xtensor<int, 1>& in_gmax,
         const double_2dvec& in_mult, const double_3dvec& coeffs);

    void
    combine(const std::vector<ScattData*>& those_scatts,
            const std::vector<double>& scalars);

    double
    calc_f(int gin, int gout, double mu);

    void
    sample(int gin, int& gout, double& mu, double& wgt, uint64_t* seed);

    size_t
    get_order() {return dist[0][0].size();};

    xt::xtensor<double, 3>
    get_matrix(size_t max_order);
};

//==============================================================================
// Function to convert Legendre functions to tabular
//==============================================================================

//! \brief Converts a ScattDatalegendre to a ScattDataHistogram
//!
//! @param leg The initial ScattDataLegendre object.
//! @param leg The resultant ScattDataTabular object.
//! @param n_mu The number of mu points to use when building the
//!   ScattDataTabular object.
void
convert_legendre_to_tabular(ScattDataLegendre& leg, ScattDataTabular& tab,
                            int n_mu);

} // namespace openmc
#endif // OPENMC_SCATTDATA_H
