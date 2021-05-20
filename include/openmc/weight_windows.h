#ifndef OPENMC_WEIGHT_WINDOWS_H
#define OPENMC_WEIGHT_WINDOWS_H 1

#include "openmc/mesh.h"
#include "openmc/particle.h"

namespace openmc {

// Weight Window Mesh class 
class WeightWindowMesh : public RectilinearMesh {
public:
  // Constructors - default
  WeightWindowMesh();

  // Constructors
  WeightWindowMesh(pugi::xml_node node);

  //! source weight biasing in energy
  //
  //! \param[in] site, particle in bank, modify the weight based on the input energy biasing
  //! \param[in] seed, random number seed
  void weight_biasing(SourceSite& site, uint64_t* seed);
   
  struct WWParams {
    double lower_weight;
    double upper_weight;
    double survival_weight;
    int max_split;
  };
  
  WWParams params;  //!< weight windows parameters
  
  // Get weight windows parameters given particle
  WeightWindowMesh::WWParams get_params(Particle& p) const;

  // weight window mesh and energy group

  int ww_type;   //!< weight window input file type
  std::vector<double> n_energy_group; //!< energy group for neutron
  std::vector<double> p_energy_group; //!< energy group for photon
  std::vector<double> n_ww_lower;  //!< lower weight window for mesh for neutron
  std::vector<double> p_ww_lower;  //!< lower weight window for mesh for photon
  bool n_ww;   //!< flag for neutron use weight window
  bool p_ww;   //!< flag for photon use weight window
  
  // WWP
  // neutron
  double n_upper_ratio;  //!< upper weight window = upper_ratio * lower weight window
  double n_survival_ratio;  //!< survival weight = survival_ratio * lower weight window
  int n_max_split;   //!< max number of split particles
  double n_multiplier;   //!< multiplier for weight window lower bounds

  // photon
  double p_upper_ratio;   //!< upper weight window = upper_ratio * lower weight window
  double p_survival_ratio;  //!< survival weight = survival_ratio * lower weight window
  int p_max_split;   //!< max number of split particles
  double p_multiplier;  //!< multiplier for weight window lower bounds

  // source weight biasing in energy
  bool user_defined_biasing;      //!< use user difined weight or not
  std::vector<double> biasing_energy;   //!< energy group for weight biasing
  std::vector<double> origin_probability; //!< probability for each group
  std::vector<double> cumulative_probability; //!< cumulative probability for each group
  std::vector<double> biasing;   //!< biasing for each energy group
  std::vector<double> cumulative_biasing;   //!< cumulative probability for biasing for each energy group
 
  // ptr for the structure mesh
  std::unique_ptr<Mesh> mesh_;

};
// Weight Window Mesh class

} // namespace openmc
#endif //OPENMC_WEIGHT_WINDOWS_H