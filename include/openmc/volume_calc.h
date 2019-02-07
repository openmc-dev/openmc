#ifndef VOLUME_CALC_H
#define VOLUME_CALC_H

#include "openmc/position.h"

#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include <string>
#include <vector>

namespace openmc {

//==============================================================================
// Volume calculation class
//==============================================================================

class VolumeCalculation {
public:
  // Aliases, types
  using int_2dvec = std::vector<std::vector<int>>;

  struct Result {
    std::array<double, 2> volume; //!< Mean/standard deviation of volume
    std::vector<int> nuclides; //!< Index of nuclides
    std::vector<double> atoms; //!< Number of atoms for each nuclide
    std::vector<double> uncertainty; //!< Uncertainty on number of atoms
  }; // Results for a single domain

  // Constructors
  VolumeCalculation() = default;
  VolumeCalculation(pugi::xml_node node);

  // Methods
  std::vector<Result> execute() const;
  void write_volume(const std::string& filename, const std::vector<Result>& results) const;

  // Data members
  int domain_type_; //!< Type of domain (cell, material, etc.)
  int n_samples_; //!< Number of samples to use
  Position lower_left_; //!< Lower-left position of bounding box
  Position upper_right_; //!< Upper-right position of bounding box
  std::vector<int> domain_ids_; //!< IDs of domains to find volumes of

private:
  void check_hit(int i_domain, int i_material, int_2dvec& indices,
    int_2dvec& hits) const;

};

//==============================================================================
// Global variables
//==============================================================================

namespace model {
extern std::vector<VolumeCalculation> volume_calcs;
}

} // namespace openmc

#endif // VOLUME_CALC_H
