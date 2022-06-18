#ifndef OPENMC_TALLIES_FILTER_H
#define OPENMC_TALLIES_FILTER_H

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <gsl/gsl>

#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/particle.h"
#include "openmc/tallies/filter_match.h"
#include "pugixml.hpp"


namespace openmc {

//==============================================================================
//! Modifies tally score events.
//==============================================================================

enum class LegendreAxis {
  x, y, z
};

enum class SphericalHarmonicsCosine {
  scatter, particle
};

class Filter
{
public:
  // Types of Filter
  enum class FilterType {
    AzimuthalFilter,
    CellFilter,
    CellInstanceFilter,
    CellbornFilter,
    CellFromFilter,
    DelayedGroupFilter,
    DistribcellFilter,
    EnergyFilter,
    EnergyoutFilter,
    EnergyFunctionFilter,
    LegendreFilter,
    MaterialFilter,
    MeshFilter,
    MeshSurfaceFilter,
    MuFilter,
    ParticleFilter,
    PolarFilter,
    SphericalHarmonicsFilter,
    SpatialLegendreFilter,
    SurfaceFilter,
    UniverseFilter,
    ZernikeFilter,
    ZernikeRadialFilter
  };

    enum class MeshDir {
    OUT_LEFT,  // x min
    IN_LEFT,  // x min
    OUT_RIGHT,  // x max
    IN_RIGHT,  // x max
    OUT_BACK,  // y min
    IN_BACK,  // y min
    OUT_FRONT,  // y max
    IN_FRONT,  // y max
    OUT_BOTTOM,  // z min
    IN_BOTTOM, // z min
    OUT_TOP, // z max
    IN_TOP // z max
  };

  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions
  
  ~Filter();

  Filter();

  //! Uses an XML input to fill the filter's data fields.
  void from_xml(pugi::xml_node node);

  void AzimuthalFilter_from_xml(pugi::xml_node node);
  void CellFilter_from_xml(pugi::xml_node node);
  void CellInstanceFilter_from_xml(pugi::xml_node node);
  void CellbornFilter_from_xml(pugi::xml_node node);
  void CellFromFilter_from_xml(pugi::xml_node node);
  void DelayedGroupFilter_from_xml(pugi::xml_node node);
  void DistribcellFilter_from_xml(pugi::xml_node node);
  void EnergyFilter_from_xml(pugi::xml_node node);
  void EnergyoutFilter_from_xml(pugi::xml_node node);
  void EnergyFunctionFilter_from_xml(pugi::xml_node node);
  void LegendreFilter_from_xml(pugi::xml_node node);
  void MaterialFilter_from_xml(pugi::xml_node node);
  void MeshFilter_from_xml(pugi::xml_node node);
  void MeshSurfaceFilter_from_xml(pugi::xml_node node);
  void MuFilter_from_xml(pugi::xml_node node);
  void ParticleFilter_from_xml(pugi::xml_node node);
  void PolarFilter_from_xml(pugi::xml_node node);
  void SphericalHarmonicsFilter_from_xml(pugi::xml_node node);
  void SpatialLegendreFilter_from_xml(pugi::xml_node node);
  void SurfaceFilter_from_xml(pugi::xml_node node);
  void UniverseFilter_from_xml(pugi::xml_node node);
  void ZernikeFilter_from_xml(pugi::xml_node node);
  void ZernikeRadialFilter_from_xml(pugi::xml_node node);

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const;

  //! Matches a tally event to a set of filter bins and weights.
  //!
  //! \param[in] p Particle being tracked
  //! \param[in] estimator Tally estimator being used
  //! \param[out] match will contain the matching bins and corresponding
  //!   weights; note that there may be zero matching bins
  virtual void
  get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const = 0;

  //! Writes data describing this filter to an HDF5 statepoint group.
  void to_statepoint(hid_t filter_group) const;

  void AzimuthalFilter_to_statepoint(hid_t filter_group) const;
  void CellFilter_to_statepoint(hid_t filter_group) const;
  void CellInstanceFilter_to_statepoint(hid_t filter_group) const;
  void CellbornFilter_to_statepoint(hid_t filter_group) const;
  void CellFromFilter_to_statepoint(hid_t filter_group) const;
  void DelayedGroupFilter_to_statepoint(hid_t filter_group) const;
  void DistribcellFilter_to_statepoint(hid_t filter_group) const;
  void EnergyFilter_to_statepoint(hid_t filter_group) const;
  void EnergyoutFilter_to_statepoint(hid_t filter_group) const;
  void EnergyFunctionFilter_to_statepoint(hid_t filter_group) const;
  void LegendreFilter_to_statepoint(hid_t filter_group) const;
  void MaterialFilter_to_statepoint(hid_t filter_group) const;
  void MeshFilter_to_statepoint(hid_t filter_group) const;
  void MeshSurfaceFilter_to_statepoint(hid_t filter_group) const;
  void MuFilter_to_statepoint(hid_t filter_group) const;
  void ParticleFilter_to_statepoint(hid_t filter_group) const;
  void PolarFilter_to_statepoint(hid_t filter_group) const;
  void SphericalHarmonicsFilter_to_statepoint(hid_t filter_group) const;
  void SpatialLegendreFilter_to_statepoint(hid_t filter_group) const;
  void SurfaceFilter_to_statepoint(hid_t filter_group) const;
  void UniverseFilter_to_statepoint(hid_t filter_group) const;
  void ZernikeFilter_to_statepoint(hid_t filter_group) const;
  void ZernikeRadialFilter_to_statepoint(hid_t filter_group) const;

  //! Return a string describing a filter bin for the tallies.out file.
  //
  //! For example, an `EnergyFilter` might return the string
  //! "Incoming Energy [0.625E-6, 20.0)".
  std::string text_label(int bin) const;

  std::string AzimuthalFilter_text_label(int bin) const;
  std::string CellFilter_text_label(int bin) const;
  std::string CellInstanceFilter_text_label(int bin) const;
  std::string CellbornFilter_text_label(int bin) const;
  std::string CellFromFilter_text_label(int bin) const;
  std::string DelayedGroupFilter_text_label(int bin) const;
  std::string DistribcellFilter_text_label(int bin) const;
  std::string EnergyFilter_text_label(int bin) const;
  std::string EnergyoutFilter_text_label(int bin) const;
  std::string EnergyFunctionFilter_text_label(int bin) const;
  std::string LegendreFilter_text_label(int bin) const;
  std::string MaterialFilter_text_label(int bin) const;
  std::string MeshFilter_text_label(int bin) const;
  std::string MeshSurfaceFilter_text_label(int bin) const;
  std::string MuFilter_text_label(int bin) const;
  std::string ParticleFilter_text_label(int bin) const;
  std::string PolarFilter_text_label(int bin) const;
  std::string SphericalHarmonicsFilter_text_label(int bin) const;
  std::string SpatialLegendreFilter_text_label(int bin) const;
  std::string SurfaceFilter_text_label(int bin) const;
  std::string UniverseFilter_text_label(int bin) const;
  std::string ZernikeFilter_text_label(int bin) const;
  std::string ZernikeRadialFilter_text_label(int bin) const;

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

private:
  FilterType type_;
  int n_bins_;
  int32_t id_ {C_NONE};
  gsl::index index_;
  std::vector<double> bins_;
  std::vector<int32_t> cells_;
  std::unordered_map<int32_t, int> map_;
  std::vector<CellInstance> cell_instances_;
  std::unordered_map<CellInstance, gsl::index, CellInstanceHash> map_;
  std::vector<int> groups_;
  int32_t cell_;
  bool matches_transport_groups_ {false};
  std::vector<double> energy_;
  double x_;
  std::vector<double> y_; //TODO: There is a collision here. ZernikeFilter has it as a scalar double, EnergyFunctionFilter has it as a vector
  double r_;
  int order_;
  std::vector<int32_t> materials_;
  int32_t mesh_;
  std::vector<Particle::Type> particles_;
  SphericalHarmonicsCosine cosine_ {SphericalHarmonicsCosine::particle};
  LegendreAxis axis_;
  double min_;
  double max_;
  std::vector<int32_t> surfaces_;
  std::vector<int32_t> universes_;
};

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  extern "C" int32_t n_filters;
  extern std::unordered_map<int, int> filter_map;
  extern Filter* tally_filters;
  extern int32_t n_tally_filters;
}

//==============================================================================
// Non-member functions
//==============================================================================

//! Make sure index corresponds to a valid filter
int verify_filter(int32_t index);

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_H
