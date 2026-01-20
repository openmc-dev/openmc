#ifndef OPENMC_VOLUME_CALC_H
#define OPENMC_VOLUME_CALC_H

#include <algorithm> // for find
#include <cstdint>
#include <string>
#include <vector>

#include "openmc/array.h"
#include "openmc/cell.h"
#include "openmc/openmp_interface.h"
#include "openmc/position.h"
#include "openmc/tallies/trigger.h"
#include "openmc/timer.h"
#include "openmc/vector.h"

#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"

namespace openmc {

//==============================================================================
// Volume calculation class
//==============================================================================

class VolumeCalculation {
private:
  // Aliases, types
  struct NuclResult {
    vector<int> nuclides;       //!< Index of nuclides
    vector<double> atoms;       //!< Number of atoms for each nuclide
    vector<double> uncertainty; //!< Uncertainty on number of atoms
  }; // Results of nuclides calculation for a single domain

public:
  //! \brief Tally corresponding to a single material (AoS)
  struct VolTally {
    double score {0.0};            //!< Current batch scores accumulator
    array<double, 2> score_acc {}; //!< Scores and squared scores accumulator
    int32_t index {0};             //!< Material ID

    VolTally() = default;
    inline VolTally(int i, double s = 0.0)
    {
      index = i;
      score = s;
    }

    //! \brief Add batch scores means to a tally
    //! \param[in] batch_size_1 Inversed batch size
    inline void finalize_batch(const double batch_size_1);

    //! \brief Pass counters data from a given tally to this
    //! \param[in] vol_tally Data source
    inline void assign_tally(const VolTally& vol_tally);

    //! \brief Add counters data from a given tally to this
    //! \param[in] vol_tally Data source
    inline void append_tally(const VolTally& vol_tally);

    //! \brief Determines given trigger condition satisfaction for this tally
    //
    //! \param[in] trigger_type Type of trigger condition
    //! \param[in] threshold    Value for trigger condition (either volume
    //! fraction variance or squared rel. err. dependent on the trigger type)
    //! \param[in] n_samples    Statistics size
    //! \return                 True if the trigger condition is satisfied
    inline bool trigger_state(const TriggerMetric trigger_type,
      const double threshold, const size_t& n_samples) const;
  }; // public just because it is used externally for the MPI struct definition

  //! \brief Online results of calculation specific for each thread
  struct CalcResults {
    uint64_t n_samples; //!< Number of samples
    int iterations;     //!< Number of iterations needed to obtain the results
    double cost; //!< Product of spent time and number of threads/processes
    Timer sampling_time; // Timer for measurment of the simulation
    vector<vector<VolumeCalculation::VolTally>>
      vol_tallies; //!< Volume tallies for each domain
    vector<VolumeCalculation::NuclResult>
      nuc_results; //!< Nuclides of each domain

    CalcResults(const VolumeCalculation& vol_calc);
    CalcResults() = default;

    //! \brief Reset all counters
    void reset();

    //! \brief Append another counters to this
    void append(const CalcResults& other);

#ifdef OPENMC_MPI
    //! \brief Collects results from all MPI processes to this
    void collect_MPI();
#endif
  };

  // Constructors
  VolumeCalculation(pugi::xml_node node);
  VolumeCalculation() = default;

  // Methods

  //! \brief Stochastically determine the volume of a set of domains along with
  //!        the average number densities of nuclides within the domain
  //
  //! \param[in,out] results Results of calculation entity for filling
  //! \return Vector of results for each user-specified domain
  void execute(CalcResults& results) const;

  //! \brief Print volume calculation results
  //
  //! \param[in] results   Full volume calculation results
  void show_results(const CalcResults& results) const;

  //! \brief Write volume calculation results to HDF5 file
  //
  //! \param[in] filename  Path to HDF5 file to write
  //! \param[in] results   Results entity
  void to_hdf5(const std::string& filename, const CalcResults& results) const;

private:
  //! \brief Rejection estimator
  inline void score_hit(const Particle& p, CalcResults& results) const;

  //! \brief Check whether a material has already been hit for a given domain.
  //! If not, add new entries to the vectors
  //
  //! \param[in] i_material      Index in global materials vector
  //! \param[in] contrib         Scoring value
  //! \param[in,out] vol_tallies Vector of tallies corresponding to each
  //! material
  void check_hit(const int32_t i_material, const double contrib,
    vector<VolTally>& vol_tallies) const;

  //! \brief Reduce vector of volumetric tallies from each thread to a single
  //! copy
  //
  //! \param[in] local_results Results specific to each thread
  //! \param[out] results      Reduced results
  void reduce_results(
    const CalcResults& local_results, CalcResults& results) const;

  //! \brief Prints a statistics parameter
  //
  //! \param[in] label   Name of parameter
  //! \param[in] units   Units of measure
  //! \param[in] value   Value of parameter
  void show_vol_stat(
    const std::string label, const std::string units, const double value) const;

  //! \brief Prints volume result for a domain
  //
  //! \param[in] domain_type  Either material or cell, ect.
  //! \param[in] domain_id    Number of this domain
  //! \param[in] region_name  Domain description
  //! \param[in] mean         Center of confidence interval
  //! \param[in] stddev       Half-width of confidence interval
  void show_volume(const std::string domain_type, const int domain_id,
    const std::string region_name, const double mean,
    const double stddev) const;

  // Tally filter and map types
  enum class TallyDomain { UNIVERSE, MATERIAL, CELL };

  // Data members
  TallyDomain domain_type_;      //!< Type of domain (cell, material, etc.)
  size_t n_samples_;             //!< Number of samples to use
  double volume_sample_;         //!< Volume of bounding primitive
  double threshold_ {-1.0};      //!< Error threshold for domain volumes
  double threshold_cnd_;         //!< Pre-computed value for trigger condition
  int max_iterations_ {INT_MAX}; //!< Limit of iterations number (necessary
                                 //!< maximum value of data type by default)
  TriggerMetric trigger_type_ {
    TriggerMetric::not_active}; //!< Trigger metric for the volume calculation
  Position lower_left_;         //!< Lower-left position of bounding box
  Position upper_right_;        //!< Upper-right position of bounding box
  vector<int> domain_ids_;      //!< IDs of domains to find volumes of

  constexpr static int _INDEX_TOTAL =
    -999; //!< Index of zero-element tally for entire domain totals should be
          //!< out of material ID range

  //! \brief Computes estimated mean and std.dev. for a tally
  //
  //! \param[in] n_samples    Statistic's size
  //! \param[in] coeff_norm   Normalization coefficient to multiply
  //! \param[in] vol_tally    Tally
  //! \return                 Array of mean and stddev
  array<double, 2> get_tally_results(const size_t& n_samples,
    const double coeff_norm, const VolTally& vol_tally) const;
};

//==============================================================================
// Global variables
//==============================================================================

namespace model {
extern vector<VolumeCalculation> volume_calcs;
}

//==============================================================================
// Non-member functions
//==============================================================================

void free_memory_volume();

} // namespace openmc

#endif // OPENMC_VOLUME_CALC_H
