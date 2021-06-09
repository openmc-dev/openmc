#ifndef OPENMC_WEIGHT_WINDOWS_H
#define OPENMC_WEIGHT_WINDOWS_H 1

#include "openmc/mesh.h"
#include "openmc/particle.h"

#define WEIGHT_CUTOFF 1.0E-38 // default low weight cutoff

namespace openmc {

void read_weight_window_xml(); //!< read the weight window file

namespace weight_window {

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

// struct to pass specific location data 
struct ParticleWeightParams {
  // 0,1,0.5,1,1e-12
  ParticleWeightParams() : lower_weight(0), upper_weight(1), survival_weight(0.5), max_split(1), weight_cutoff(WEIGHT_CUTOFF) {}
  // constructor
  ParticleWeightParams(const WWParams ww_settings, const int indices) {
    // set the weight for the current location
    lower_weight = ww_settings.multiplier*
                   ww_settings.lower_ww[indices];  
    // set the upper weight bound
    upper_weight = ww_settings.upper_ratio*
                    lower_weight;
    // set the survival weight
    survival_weight = lower_weight*ww_settings.survival_ratio;
    // set the max split
    max_split = ww_settings.max_split;

    // set the weight cutoff
    weight_cutoff = ww_settings.weight_cutoff;  
  }

  double lower_weight;
  double upper_weight;
  double survival_weight;
  int max_split;
  double weight_cutoff;
};

// Weight Window Mesh class 
class WeightWindow {
public:
  // Constructors - default
  WeightWindow();

  // Constructors
  WeightWindow(pugi::xml_node node);

  // Get weight windows parameters given particle
  ParticleWeightParams get_params(Particle& p) const;

  bool n_ww {false};   //!< flag for neutron use weight window
  bool p_ww {false};   //!< flag for photon use weight window

private:
  
  WWParams read_particle_settings(pugi::xml_node node);
  
  // ptr for the structure mesh
  std::shared_ptr<RectilinearMesh> mesh_; //!< The mesh for the problem
  std::map<ParticleType,WWParams> weight_params;  //!< weight windows parameters
};
// Weight Window class
} // namespace weight_window
} // namespace openmc
#endif //OPENMC_WEIGHT_WINDOWS_H
