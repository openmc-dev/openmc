#ifndef OPENMC_VOLUME_CALC_H
#define OPENMC_VOLUME_CALC_H

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
  struct Result {
    std::array<double, 2> volume; //!< Mean/standard deviation of volume
    std::vector<int> nuclides; //!< Index of nuclides
    std::vector<double> atoms; //!< Number of atoms for each nuclide
    std::vector<double> uncertainty; //!< Uncertainty on number of atoms
  }; // Results for a single domain

  // Constructors
  VolumeCalculation(pugi::xml_node node);

  // Methods

  //! \brief Stochastically determine the volume of a set of domains along with the
  //!   average number densities of nuclides within the domain
  //
  //! \return Vector of results for each user-specified domain
  std::vector<Result> execute() const;

  //! \brief Write volume calculation results to HDF5 file
  //
  //! \param[in] filename  Path to HDF5 file to write
  //! \param[in] results   Vector of results for each domain
  void to_hdf5(const std::string& filename, const std::vector<Result>& results) const;

  // Data members
  int domain_type_; //!< Type of domain (cell, material, etc.)
  int n_samples_; //!< Number of samples to use
  Position lower_left_; //!< Lower-left position of bounding box
  Position upper_right_; //!< Upper-right position of bounding box
  std::vector<int> domain_ids_; //!< IDs of domains to find volumes of

private:
  //! \brief Check whether a material has already been hit for a given domain.
  //! If not, add new entries to the vectors
  //
  //! \param[in] i_material Index in global materials vector
  //! \param[in,out] indices Vector of material indices
  //! \param[in,out] hits Number of hits corresponding to each material
  void check_hit(int i_material, std::vector<int>& indices,
    std::vector<int>& hits) const;

};

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  extern std::vector<VolumeCalculation> volume_calcs;
}

//==============================================================================
// Non-member functions
//==============================================================================

void free_memory_volume();

} // namespace openmc

#endif // OPENMC_VOLUME_CALC_H
