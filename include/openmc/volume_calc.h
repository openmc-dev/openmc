#ifndef OPENMC_VOLUME_CALC_H
#define OPENMC_VOLUME_CALC_H

#include "openmc/position.h"

#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

#include <string>
#include <vector>
#include <gsl/gsl>

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
    size_t num_samples;

    Result& operator +=( const Result& other) {
      Expects(volume.size() == other.volume.size());
      Expects(atoms.size() == atoms.size());

      auto& a_samples = num_samples;
      auto& b_samples = other.num_samples;

      size_t total_samples = num_samples + other.num_samples;

      for (int i = 0; i < volume.size(); i++) {
        // calculate weighted average of volume results
        auto& a_vol = volume[0];
        auto& b_vol = other.volume[0];
        volume[0] = (a_samples * a_vol + b_samples * b_vol) / total_samples;

        // propagate error
        auto& a_err = volume[1];
        auto& b_err = other.volume[1];
        volume[1] = std::sqrt(a_samples * a_err * a_err + b_samples * b_err * b_err) / total_samples;
      }

      for (int i = 0; i < atoms.size(); i++) {
        // calculate weighted average of atom results
        auto& a_atoms = atoms[i];
        auto& b_atoms = other.atoms[i];
        atoms[i] = (a_samples * a_atoms + b_samples * b_atoms) / total_samples;

        // propagate error
        auto& a_err = uncertainty[i];
        auto& b_err = other.uncertainty[i];
        uncertainty[i] = std::sqrt(a_samples * a_err * a_err + b_samples * b_err * b_err) / total_samples;
      }

      // update number of samples on the returned set of results;
      num_samples = total_samples;

      return *this;
    }
    
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
  size_t n_samples_; //!< Number of samples to use
  double trigger_ {-1.0}; //!< Error threshold for domain volumes
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

  //! \brief Perform calculation for domain volumes and average nuclide density 
  //!   using n_samples_
  //
  //! \param[in] seed_offset Seed offset used for independent calculations
  //! \return Vector of results for each user-specified domain
  std::vector<Result> _execute(size_t seed_offset = 0) const;


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
