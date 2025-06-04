#ifndef OPENMC_TALLIES_FILTER_MESHMATERIAL_H
#define OPENMC_TALLIES_FILTER_MESHMATERIAL_H

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "openmc/position.h"
#include "openmc/random_ray/source_region.h"
#include "openmc/span.h"
#include "openmc/tallies/filter.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
//! Helper structs that define a combination of a mesh element index and a
//! material index and a functor for hashing to place in an unordered_map
//==============================================================================

struct ElementMat {
  //! Check for equality
  bool operator==(const ElementMat& other) const
  {
    return index_element == other.index_element && index_mat == other.index_mat;
  }

  int32_t index_element;
  int32_t index_mat;
};

struct ElementMatHash {
  std::size_t operator()(const ElementMat& k) const
  {
    size_t seed = 0;
    hash_combine(seed, k.index_element);
    hash_combine(seed, k.index_mat);
    return seed;
  }
};

//==============================================================================
//! Indexes the location of particle events to combinations of mesh element
//! index and material.  For tracklength tallies, it will produce multiple valid
//! bins and the bin weight will correspond to the fraction of the track length
//! that lies in that bin.
//==============================================================================

class MeshMaterialFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~MeshMaterialFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "meshmaterial"; }
  FilterType type() const override { return FilterType::MESH_MATERIAL; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors

  int32_t mesh() const { return mesh_; }

  void set_mesh(int32_t mesh);

  void set_bins(vector<ElementMat>&& bins);

  virtual void set_translation(const Position& translation);

  virtual void set_translation(const double translation[3]);

  virtual const Position& translation() const { return translation_; }

  virtual bool translated() const { return translated_; }

private:
  //----------------------------------------------------------------------------
  // Data members

  int32_t mesh_;            //!< Index of the mesh
  bool translated_ {false}; //!< Whether or not the filter is translated
  Position translation_ {0.0, 0.0, 0.0}; //!< Filter translation

  //! The indices of the mesh element-material combinations binned by this
  //! filter.
  vector<ElementMat> bins_;

  //! The set of materials used in this filter
  std::unordered_set<int32_t> materials_;

  //! A map from mesh element-material indices to filter bin indices.
  std::unordered_map<ElementMat, int32_t, ElementMatHash> map_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MESHMATERIAL_H
