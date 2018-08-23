//! \file xsdata.h
//! A collection of classes for containing the Multi-Group Cross Section data

#ifndef OPENMC_XSDATA_H
#define OPENMC_XSDATA_H

#include <memory>
#include <vector>

#include "openmc/hdf5_interface.h"
#include "openmc/scattdata.h"


namespace openmc {

//==============================================================================
// XSDATA contains the temperature-independent cross section data for an MGXS
//==============================================================================

class XsData {

  private:
    //! \brief Reads scattering data from the HDF5 file
    void
    scatter_from_hdf5(hid_t xsdata_grp, int n_pol, int n_azi, int energy_groups,
         int scatter_format, int final_scatter_format, int order_data,
         int max_order, int legendre_to_tabular_points);

    //! \brief Reads fission data from the HDF5 file
    void
    fission_from_hdf5(hid_t xsdata_grp, int n_pol, int n_azi, int energy_groups,
         int delayed_groups, bool is_isotropic);

  public:

    // The following quantities have the following dimensions:
    // [angle][incoming group]
    double_2dvec total;
    double_2dvec absorption;
    double_2dvec nu_fission;
    double_2dvec prompt_nu_fission;
    double_2dvec kappa_fission;
    double_2dvec fission;
    double_2dvec inverse_velocity;

    // decay_rate has the following dimensions:
    // [angle][delayed group]
    double_2dvec decay_rate;
    // delayed_nu_fission has the following dimensions:
    // [angle][incoming group][delayed group]
    double_3dvec delayed_nu_fission;
    // chi_prompt has the following dimensions:
    // [angle][incoming group][outgoing group]
    double_3dvec chi_prompt;
    // chi_delayed has the following dimensions:
    // [angle][incoming group][outgoing group][delayed group]
    double_4dvec chi_delayed;
    // scatter has the following dimensions: [angle]
    std::vector<std::shared_ptr<ScattData> > scatter;

    XsData() = default;

    //! \brief Constructs the XsData object metadata.
    //!
    //! @param num_groups Number of energy groups.
    //! @param num_delayed_groups Number of delayed groups.
    //! @param fissionable Is this a fissionable data set or not.
    //! @param scatter_format The scattering representation of the file.
    //! @param n_pol Number of polar angles.
    //! @param n_azi Number of azimuthal angles.
    XsData(int num_groups, int num_delayed_groups, bool fissionable,
           int scatter_format, int n_pol, int n_azi);

    //! \brief Loads the XsData object from the HDF5 file
    //!
    //! @param xs_id HDF5 group id for the cross section data.
    //! @param fissionable Is this a fissionable data set or not.
    //! @param scatter_format The scattering representation of the file.
    //! @param final_scatter_format The scattering representation after reading;
    //!   this is different from scatter_format if converting a Legendre to
    //!   a tabular representation.
    //! @param order_data The dimensionality of the scattering data in the file.
    //! @param max_order Maximum order requested by the user;
    //!   this is only used for Legendre scattering.
    //! @param legendre_to_tabular Flag to denote if any Legendre provided
    //!   should be converted to a Tabular representation.
    //! @param legendre_to_tabular_points If a conversion is requested, this
    //!   provides the number of points to use in the tabular representation.
    //! @param is_isotropic Is this an isotropic or angular with respect to
    //!   the incoming particle.
    //! @param n_pol Number of polar angles.
    //! @param n_azi Number of azimuthal angles.
    void
    from_hdf5(hid_t xsdata_grp, bool fissionable, int scatter_format,
         int final_scatter_format, int order_data, int max_order,
         int legendre_to_tabular_points, bool is_isotropic, int n_pol,
         int n_azi);

    //! \brief Combines the microscopic data to a macroscopic object.
    //!
    //! @param micros Microscopic objects to combine.
    //! @param scalars Scalars to multiply the microscopic data by.
    void
    combine(const std::vector<XsData*>& those_xs, const double_1dvec& scalars);

    //! \brief Checks to see if this and that are able to be combined
    //!
    //! This comparison is used when building macroscopic cross sections
    //! from microscopic cross sections.
    //! @param that The other XsData to compare to this one.
    //! @return True if they can be combined.
    bool
    equiv(const XsData& that);
};


} //namespace openmc
#endif // OPENMC_XSDATA_H
