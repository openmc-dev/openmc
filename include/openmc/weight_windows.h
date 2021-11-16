#ifndef OPENMC_WEIGHT_WINDOWS_H
#define OPENMC_WEIGHT_WINDOWS_H 1

#include "openmc/mesh.h"
#include "openmc/particle.h"

#define WEIGHT_CUTOFF 1.0E-38 // default low weight cutoff

namespace openmc {

void read_weight_windows(pugi::xml_node node);
void read_variance_reduction_xml();

class WeightWindowDomain;
class WeightWindowParameters;

namespace variance_reduction {

extern std::unordered_map<int32_t, int32_t> ww_domain_map;
extern openmc::vector<unique_ptr<WeightWindowDomain>> ww_domains;

extern std::unordered_map<int32_t, int32_t> ww_map;
extern openmc::vector<unique_ptr<WeightWindowParameters>> ww_params;

} // namespace variance_reduction

class WeightWindowParameters {
public:
  // Constructors - default
  WeightWindowParameters();

  WeightWindowParameters(pugi::xml_node node);

  // Methods
  void to_statepoint(hid_t group) const;

  // Accessors
  int32_t id() { return id_; }

  const ParticleType& particle_type() const { return particle_type_; }
  ParticleType& particle_type() { return particle_type_; }

  const std::vector<double>& energy_bins() const { return energy_bins_; }
  std::vector<double>& energy_bins() { return energy_bins_; }

  const std::vector<double>& lower_ww() const { return lower_ww_; }
  std::vector<double>& lower_ww() { return lower_ww_; }

  const std::vector<double>& upper_ww() const { return upper_ww_; }
  std::vector<double>& upper_ww() { return upper_ww_; }

  int max_split() const { return max_split_; }
  int& max_split() { return max_split_; }

  double survival_ratio() const { return survival_ratio_; }
  double& survival_ratio() { return survival_ratio_; }

  double weight_cutoff() const { return weight_cutoff_; }
  double& weight_cutoff() { return weight_cutoff_; }

private:
  int32_t id_;
  ParticleType particle_type_;
  openmc::vector<double> energy_bins_;
  openmc::vector<double> lower_ww_;
  openmc::vector<double> upper_ww_;
  double survival_ratio_;
  int max_split_;
  double weight_cutoff_;
};

// struct to pass specific location data
struct ParticleWeightParams {
  // 0,1,0.5,1,1e-12
  ParticleWeightParams() : lower_weight(0), upper_weight(1), survival_weight(0.5), max_split(1), weight_cutoff(WEIGHT_CUTOFF) {}
  // constructor
  ParticleWeightParams(
    const unique_ptr<WeightWindowParameters>& params, const int indices)
  {
    // set the weight for the current location
    lower_weight = params->lower_ww()[indices];
    // set the upper weight bound
    upper_weight = params->upper_ww()[indices];
    // set the survival weight
    survival_weight = lower_weight * params->survival_ratio();
    // set the max split
    max_split = params->max_split();
    // set the weight cutoff
    weight_cutoff = params->weight_cutoff();
  }

  double lower_weight;
  double upper_weight;
  double survival_weight;
  int max_split;
  double weight_cutoff;
};


class WeightWindowDomain {
public:
  // Constructrors - default
  WeightWindowDomain();

  // Constructors
  WeightWindowDomain(pugi::xml_node node);

  // Methods
  bool find_params(Particle& p, ParticleWeightParams& params) const;

  void to_statepoint(hid_t group) const;

  int32_t id() const { return id_; }

  int32_t mesh_idx() const { return mesh_idx_; }

  int32_t param_idx() const { return param_idx_; }

private:
  int32_t id_;        //!< the id of the weight window domain
  int32_t mesh_idx_;  //!< the idx of the mesh this domain uses
  int32_t param_idx_; //!< the idx of the ww params this domain uses
};

} // namespace openmc
#endif //OPENMC_WEIGHT_WINDOWS_H
