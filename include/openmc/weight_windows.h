#ifndef OPENMC_WEIGHT_WINDOWS_H
#define OPENMC_WEIGHT_WINDOWS_H 1

#include "openmc/mesh.h"
#include "openmc/particle.h"

#define WEIGHT_CUTOFF 1.0E-38 // default low weight cutoff

namespace openmc {

void read_variance_reduction_xml(); //!< read the weight window file

namespace weight_window {

// domain map provides
extern std::unordered_map<int32_t, int32_t> ww_domain_map;
extern openmc::vector<unique_ptr<>> ww_domains;

// ww_map provides
extern std::unordered_map<int32_t, int32_t> ww_map;
extern openmc::vector<unique_ptr<WeightWindowParameters>> ww_params;
  
// struct to pass specific location data 
struct ParticleWeightParams {
  // 0,1,0.5,1,1e-12
  ParticleWeightParams() : lower_weight(0), upper_weight(1), survival_weight(0.5), max_split(1), weight_cutoff(WEIGHT_CUTOFF) {}
  // constructor
  ParticleWeightParams(const unique_ptr<WeightWindowParameters> params,
		       const int indices) {
    
    // set the weight for the current location
    lower_weight = params._lower_ww[indices];  
    // set the upper weight bound
    upper_weight = ww_settings._upper_ww[indices];
    // set the survival weight
    survival_weight = lower_weight*params._survival_ratio;
    // set the max split
    max_split = params._max_split;
    // set the weight cutoff
    weight_cutoff = params._weight_cutoff;  
  }

  double lower_weight;
  double upper_weight;
  double survival_weight;
  int max_split;
  double weight_cutoff;
};

class WeightWindowParameters {
public:
  // Constructors - default
  WeightWindowParameters();

  WeightWindowParameters(pugi::xml_node node);

private:
  
  openmc::ParticleType _particle_type; //!< 
  double _survival_ratio;
  int _max_split;
  double _weight_cutoff;
  std::vector<double> _energy_bounds;
  std::vector<double> _lower_ww;
  std::vector<double> _upper_ww;
  int32_t _id;
}
  
class WeightWindowDomain {
public:
  // Constructrors - default
  WeightWindowDomain();

  // Constructors
  WeightWindowDomain(pugi::xml_node node);
  
  // Constructors
  WeightWindowDomain(int32_t mesh_idx, int32_t param_idx);

private:
  int32_t _ww_domain_id; //!< the id of the weight window domain
  int32_t _ww_mesh_idx;  //!< the idx of the mesh this domain uses
  int32_t _ww_param_idx; //!< the idx of the ww params this domain uses
  
}
  
// Weight Window class 
class WeightWindow {
public:
  // Constructors - default
  WeightWindow();

  // Constructors
  WeightWindow(pugi::xml_node node);

  // Get weight windows parameters given particle
  ParticleWeightParams get_params(Particle& p) const;
    
};
  
// Weight Window class
} // namespace weight_window
} // namespace openmc
#endif //OPENMC_WEIGHT_WINDOWS_H
