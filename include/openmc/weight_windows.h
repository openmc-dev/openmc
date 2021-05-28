#ifndef OPENMC_WEIGHT_WINDOWS_H
#define OPENMC_WEIGHT_WINDOWS_H 1

#include "openmc/mesh.h"
#include "openmc/particle.h"

namespace openmc {

void read_weight_window_xml(); //!< read the weight window file

namespace weight_window {

// struct to pass specific location data 
struct ParticleWeightParams {
  double lower_weight;
  double upper_weight;
  double survival_weight;
  int max_split;
  double weight_cutoff;
};

// generic struct to contain the ww settings
struct WWParams {
  double upper_ratio;
  double multiplier;
  double survival_ratio;
  int max_split;
  double weight_cutoff;
  std::vector<double> energy_bounds;
  std::vector<double> lower_ww;
  std::vector<double> upper_ww;
};

// Weight Window Mesh class 
class WeightWindow {
public:
  // Constructors - default
  WeightWindow();

  // Constructors
  WeightWindow(pugi::xml_node node);

  //! source weight biasing in energy
  //
  //! \param[in] site, particle in bank, modify the weight based on the input energy biasing
  //! \param[in] seed, random number seed
  void weight_biasing(SourceSite& site, uint64_t* seed);

  WWParams read_particle_settings(pugi::xml_node node);

  void read_user_biasing(pugi::xml_node node);
  
  // Get weight windows parameters given particle
  ParticleWeightParams get_params(Particle& p) const;

  bool n_ww;   //!< flag for neutron use weight window
  bool p_ww;   //!< flag for photon use weight window
  
  // source weight biasing in energy
  bool user_defined_biasing;      //!< use user difined weight or not
  std::vector<double> biasing_energy;   //!< energy group for weight biasing
  std::vector<double> origin_probability; //!< probability for each group
  std::vector<double> cumulative_probability; //!< cumulative probability for each group
  std::vector<double> biasing;   //!< biasing for each energy group
  std::vector<double> cumulative_biasing;   //!< cumulative probability for biasing for each energy group
 
  // ptr for the structure mesh
  std::shared_ptr<RectilinearMesh> mesh_; //!< The mesh for the problem
  std::map<ParticleType,WWParams> weight_params;  //!< weight windows parameters
};
// Weight Window class
} // namespace weight_window
} // namespace openmc
#endif //OPENMC_WEIGHT_WINDOWS_H
