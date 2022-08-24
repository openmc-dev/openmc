#ifndef OPENMC_TALLIES_FILTER_H
#define OPENMC_TALLIES_FILTER_H

#include <cmath>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <gsl/gsl>

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mesh.h"
#include "openmc/particle.h"
#include "openmc/static_map.h"
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

  Filter(pugi::xml_node node, int32_t index);

  //! Uses an XML input to fill the filter's data fields.
  void from_xml(pugi::xml_node node);

  void AzimuthalFilter_from_xml(pugi::xml_node node);
  void CellFilter_from_xml(pugi::xml_node node);
  void CellInstanceFilter_from_xml(pugi::xml_node node);
  //void CellbornFilter_from_xml(pugi::xml_node node);
  //void CellFromFilter_from_xml(pugi::xml_node node);
  void DelayedGroupFilter_from_xml(pugi::xml_node node);
  void DistribcellFilter_from_xml(pugi::xml_node node);
  void EnergyFilter_from_xml(pugi::xml_node node);
  //void EnergyoutFilter_from_xml(pugi::xml_node node);
  void EnergyFunctionFilter_from_xml(pugi::xml_node node);
  void LegendreFilter_from_xml(pugi::xml_node node);
  void MaterialFilter_from_xml(pugi::xml_node node);
  void MeshFilter_from_xml(pugi::xml_node node);
  //void MeshSurfaceFilter_from_xml(pugi::xml_node node);
  void MuFilter_from_xml(pugi::xml_node node);
  void ParticleFilter_from_xml(pugi::xml_node node);
  void PolarFilter_from_xml(pugi::xml_node node);
  void SphericalHarmonicsFilter_from_xml(pugi::xml_node node);
  void SpatialLegendreFilter_from_xml(pugi::xml_node node);
  void SurfaceFilter_from_xml(pugi::xml_node node);
  void UniverseFilter_from_xml(pugi::xml_node node);
  void ZernikeFilter_from_xml(pugi::xml_node node);
  //void ZernikeRadialFilter_from_xml(pugi::xml_node node);

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const;
  FilterType get_type() const {return type_;}
  FilterType get_filter_type(const std::string& type);

  //! Matches a tally event to a set of filter bins and weights.
  //!
  //! \param[in] p Particle being tracked
  //! \param[in] estimator Tally estimator being used
  //! \param[out] match will contain the matching bins and corresponding
  //!   weights; note that there may be zero matching bins
  void get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;

  void AzimuthalFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void CellFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void CellInstanceFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void CellbornFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void CellFromFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void DelayedGroupFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void DistribcellFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void EnergyFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void EnergyoutFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void EnergyFunctionFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void LegendreFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void MaterialFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void MeshFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void MeshSurfaceFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void MuFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void ParticleFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void PolarFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void SphericalHarmonicsFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void SpatialLegendreFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void SurfaceFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void UniverseFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void ZernikeFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;
  void ZernikeRadialFilter_get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match) const;


  //! Writes data describing this filter to an HDF5 statepoint group.
  void to_statepoint(hid_t filter_group) const;

  void AzimuthalFilter_to_statepoint(hid_t filter_group) const;
  void CellFilter_to_statepoint(hid_t filter_group) const;
  void CellInstanceFilter_to_statepoint(hid_t filter_group) const;
  //void CellbornFilter_to_statepoint(hid_t filter_group) const;
  //void CellFromFilter_to_statepoint(hid_t filter_group) const;
  void DelayedGroupFilter_to_statepoint(hid_t filter_group) const;
  void DistribcellFilter_to_statepoint(hid_t filter_group) const;
  void EnergyFilter_to_statepoint(hid_t filter_group) const;
  //void EnergyoutFilter_to_statepoint(hid_t filter_group) const;
  void EnergyFunctionFilter_to_statepoint(hid_t filter_group) const;
  void LegendreFilter_to_statepoint(hid_t filter_group) const;
  void MaterialFilter_to_statepoint(hid_t filter_group) const;
  void MeshFilter_to_statepoint(hid_t filter_group) const;
  //void MeshSurfaceFilter_to_statepoint(hid_t filter_group) const;
  void MuFilter_to_statepoint(hid_t filter_group) const;
  void ParticleFilter_to_statepoint(hid_t filter_group) const;
  void PolarFilter_to_statepoint(hid_t filter_group) const;
  void SphericalHarmonicsFilter_to_statepoint(hid_t filter_group) const;
  void SpatialLegendreFilter_to_statepoint(hid_t filter_group) const;
  void SurfaceFilter_to_statepoint(hid_t filter_group) const;
  void UniverseFilter_to_statepoint(hid_t filter_group) const;
  void ZernikeFilter_to_statepoint(hid_t filter_group) const;
  //void ZernikeRadialFilter_to_statepoint(hid_t filter_group) const;

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

  // Superset of type specific accessors

  // Defined in header

  const vector<int32_t>& cells() const { return cells_; }
  const vector<CellInstance>& cell_instances() const { return cell_instances_; }
  const vector<int>& groups() const { return groups_; }
  int32_t cell() const { return cell_; }
  const vector<double>& energy() const { return energy_; }
  const vector<double>& y() const { return y_; }
  const vector<double>& bins() const { return bins_; }
  bool matches_transport_groups() const { return matches_transport_groups_; }
  int order() const { return order_; }
  vector<int32_t>& materials() { return materials_; }
  const vector<int32_t>& materials() const { return materials_; }
  int32_t mesh() const {return mesh_;}
  void set_mesh(int32_t mesh){
    mesh_ = mesh;
    if (type_ == FilterType::MeshFilter)
      n_bins_ = model::meshes[mesh_].n_bins();
    if (type_ == FilterType::MeshSurfaceFilter)
      n_bins_ = model::meshes[mesh_].n_surface_bins();
  }
  const vector<Particle::Type>& particles() const { return particles_; }
  SphericalHarmonicsCosine cosine() const { return cosine_; }
  LegendreAxis axis() const { return axis_; }
  double min() const { return min_; }
  double max() const { return max_; }
  double x() const { return x_; }
  void set_x(double x) { x_ = x; }
  double yy() const { return yy_; }
  void set_y(double y) { yy_ = y; }
  double r() const { return r_; }
  void set_r(double r) { r_ = r; }

  // Defined once

  void set_cells(gsl::span<int32_t> cells);
  void set_cell_instances(gsl::span<CellInstance> instances);
  void set_cell(int32_t cell);
  void set_groups(gsl::span<int> groups);
  void set_data(gsl::span<const double> energy, gsl::span<const double> y);
  void set_materials(gsl::span<const int32_t> materials);
  void set_particles(gsl::span<Particle::Type> particles);
  void set_cosine(gsl::cstring_span cosine);
  void set_axis(LegendreAxis axis);
  void set_minmax(double min, double max);
  void set_surfaces(gsl::span<int32_t> surfaces);
  void set_universes(gsl::span<int32_t> universes);

  // Defined by class (need dispatching)

  void set_bins(gsl::span<const double> bins);

  void AzimuthalFilter_set_bins(gsl::span<const double> bins);
  void EnergyFilter_set_bins(gsl::span<const double> bins);
  void MuFilter_set_bins(gsl::span<const double> bins);
  void PolarFilter_set_bins(gsl::span<const double> bins);

  void set_order(int order);

  void LegendreFilter_set_order(int order);
  void SphericalHarmonicsFilter_set_order(int order);
  void SpatialLegendreFilter_set_order(int order);
  void ZernikeFilter_set_order(int order);
  void ZernikeRadialFilter_set_order(int order);

  void copy_to_device();

  //----------------------------------------------------------------------------
  // Data members

private:
  FilterType type_;
  int n_bins_;
  int32_t id_ {C_NONE};
  gsl::index index_;
  vector<double> bins_;
  vector<int32_t> cells_;
  static_map<int32_t, int> map_;
  vector<CellInstance> cell_instances_;
  static_map<CellInstance, gsl::index, CellInstanceHash> imap_;
  vector<int> groups_;
  int32_t cell_;
  bool matches_transport_groups_ {false};
  vector<double> energy_;
  double x_;
  // Note: there was a collision between EnergyFunction (vector y_) and Zernike (scalar y_), so Zernike uses yy_
  vector<double> y_;
  double yy_;
  double r_;
  int order_;
  vector<int32_t> materials_;
  int32_t mesh_;
  vector<Particle::Type> particles_;
  SphericalHarmonicsCosine cosine_ {SphericalHarmonicsCosine::particle};
  LegendreAxis axis_;
  double min_;
  double max_;
  vector<int32_t> surfaces_;
  vector<int32_t> universes_;
};

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  extern "C" int32_t n_filters;
  // TODO: Need to declare the filter_map as target
  extern std::unordered_map<int, int> filter_map;
  #pragma omp declare target
  extern Filter* tally_filters;
  extern int32_t n_tally_filters;
  #pragma omp end declare target
}

//==============================================================================
// Non-member functions
//==============================================================================

//! Make sure index corresponds to a valid filter
int verify_filter(int32_t index);

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_H
