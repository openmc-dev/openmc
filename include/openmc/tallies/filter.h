#ifndef OPENMC_TALLIES_FILTER_H
#define OPENMC_TALLIES_FILTER_H

#include <cstdint>
#include <string>
#include <unordered_map>

#include "pugixml.hpp"
#include <gsl/gsl-lite.hpp>

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/memory.h"
#include "openmc/particle.h"
#include "openmc/tallies/filter_match.h"
#include "openmc/vector.h"

namespace openmc {

enum class FilterType {
  AZIMUTHAL,
  CELLBORN,
  CELLFROM,
  CELL,
  CELL_INSTANCE,
  COLLISION,
  DELAYED_GROUP,
  DISTRIBCELL,
  ENERGY_FUNCTION,
  ENERGY,
  ENERGY_OUT,
  LEGENDRE,
  MATERIAL,
  MATERIALFROM,
  MESH,
  MESHBORN,
  MESH_SURFACE,
  MU,
  PARTICLE,
  POLAR,
  SPHERICAL_HARMONICS,
  SPATIAL_LEGENDRE,
  SURFACE,
  TIME,
  UNIVERSE,
  ZERNIKE,
  ZERNIKE_RADIAL
};

//==============================================================================
//! Modifies tally score events.
//==============================================================================

class Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions

  Filter();
  virtual ~Filter();

  //! Create a new tally filter
  //
  //! \tparam T Type of the filter
  //! \param[in] id  Unique ID for the filter. If none is passed, an ID is
  //!    automatically assigned
  //! \return Pointer to the new filter object
  template<typename T>
  static T* create(int32_t id = -1);

  //! Create a new tally filter
  //
  //! \param[in] type  Type of the filter
  //! \param[in] id  Unique ID for the filter. If none is passed, an ID is
  //!    automatically assigned
  //! \return Pointer to the new filter object
  static Filter* create(const std::string& type, int32_t id = -1);

  //! Create a new tally filter from an XML node
  //
  //! \param[in] node XML node
  //! \return Pointer to the new filter object
  static Filter* create(pugi::xml_node node);

  //! Uses an XML input to fill the filter's data fields.
  virtual void from_xml(pugi::xml_node node) = 0;

  //----------------------------------------------------------------------------
  // Methods

  virtual std::string type_str() const = 0;
  virtual FilterType type() const = 0;

  //! Matches a tally event to a set of filter bins and weights.
  //!
  //! \param[in] p Particle being tracked
  //! \param[in] estimator Tally estimator being used
  //! \param[out] match will contain the matching bins and corresponding
  //!   weights; note that there may be zero matching bins
  virtual void get_all_bins(
    const Particle& p, TallyEstimator estimator, FilterMatch& match) const = 0;

  //! Writes data describing this filter to an HDF5 statepoint group.
  virtual void to_statepoint(hid_t filter_group) const
  {
    write_dataset(filter_group, "type", type_str());
    write_dataset(filter_group, "n_bins", n_bins_);
  }

  //! Return a string describing a filter bin for the tallies.out file.
  //
  //! For example, an `EnergyFilter` might return the string
  //! "Incoming Energy [0.625E-6, 20.0)".
  virtual std::string text_label(int bin) const = 0;

  //----------------------------------------------------------------------------
  // Accessors

  //! Get unique ID of filter
  //! \return Unique ID
  int32_t id() const { return id_; }

  //! Assign a unique ID to the filter
  //! \param[in]  Unique ID to assign. A value of -1 indicates that an ID should
  //!   be automatically assigned
  void set_id(int32_t id);

  //! Get number of bins
  //! \return Number of bins
  int n_bins() const { return n_bins_; }

  gsl::index index() const { return index_; }

  //----------------------------------------------------------------------------
  // Data members

protected:
  int n_bins_;

private:
  int32_t id_ {C_NONE};
  gsl::index index_;
};

//==============================================================================
// Global variables
//==============================================================================

namespace model {
extern "C" int32_t n_filters;
extern std::unordered_map<int, int> OPENMC_API filter_map;
extern vector<unique_ptr<Filter>> tally_filters;
} // namespace model

//==============================================================================
// Non-member functions
//==============================================================================

//! Make sure index corresponds to a valid filter
int verify_filter(int32_t index);

//==============================================================================
// Filter implementation
//==============================================================================

template<typename T>
T* Filter::create(int32_t id)
{
  static_assert(std::is_base_of<Filter, T>::value,
    "Type specified is not derived from openmc::Filter");
  // Create filter and add to filters vector
  auto filter = make_unique<T>();
  auto ptr_out = filter.get();
  model::tally_filters.emplace_back(std::move(filter));
  // Assign ID
  model::tally_filters.back()->set_id(id);

  return ptr_out;
}

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_H
