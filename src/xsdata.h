//! \file xsdata.h
//! A collection of classes for containing the Multi-Group Cross Section data

#ifndef XSDATA_H
#define XSDATA_H

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <string>
#include <valarray>
#include <vector>
#include <iostream>

#include "constants.h"
#include "hdf5_interface.h"
#include "math_functions.h"
#include "random_lcg.h"
#include "scattdata.h"


namespace openmc {

//==============================================================================
// XSDATA contains the temperature-independent cross section data for an MGXS
//==============================================================================

class XsData {
  private:
    void _scatter_from_hdf5(hid_t xsdata_grp, int n_pol, int n_azi,
         int energy_groups, int scatter_format, int final_scatter_format,
         int order_data, int max_order, int legendre_to_tabular_points);
    void _fissionable_from_hdf5(hid_t xsdata_grp, int n_pol, int n_azi,
         int energy_groups, int delayed_groups, bool is_isotropic);
  public:
    // The following quantities have the following dimensions:
    // [phi][theta][incoming group]
    double_3dvec total;
    double_3dvec absorption;
    double_3dvec nu_fission;
    double_3dvec prompt_nu_fission;
    double_3dvec kappa_fission;
    double_3dvec fission;
    double_3dvec inverse_velocity;
    // decay_rate has the following dimensions:
    // [phi][theta][delayed group]
    double_3dvec decay_rate;
    // delayed_nu_fission has the following dimensions:
    // [phi][theta][incoming group][delayed group]
    double_4dvec delayed_nu_fission;
    // chi_prompt has the following dimensions:
    // [phi][theta][incoming group][outgoing group]
    double_4dvec chi_prompt;
    // chi_delayed has the following dimensions:
    // [phi][theta][incoming group][outgoing group][delayed group]
    double_5dvec chi_delayed;
    // scatter has the following dimensions: [phi][theta]
    std::vector<std::vector<ScattData*> > scatter;

    XsData() = default;
    XsData(int num_groups, int num_delayed_groups, bool fissionable,
           int scatter_format, int n_pol, int n_azi);
    void from_hdf5(hid_t xsdata_grp, bool fissionable, int scatter_format,
                   int final_scatter_format, int order_data, int max_order,
                   int legendre_to_tabular_points, bool is_isotropic);
    void combine(std::vector<XsData*> those_xs, double_1dvec& scalars);
    bool equiv(const XsData& that);
};


} //namespace openmc
#endif // XSDATA_H