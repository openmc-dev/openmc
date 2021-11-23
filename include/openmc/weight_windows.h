#ifndef OPENMC_WEIGHT_WINDOWS_H
#define OPENMC_WEIGHT_WINDOWS_H 1

#include "openmc/mesh.h"
#include "openmc/particle.h"

namespace openmc {

void read_weight_windows(pugi::xml_node node);
void read_variance_reduction_xml();

class WeightWindowDomain;
class WeightWindows;

namespace variance_reduction {

extern std::unordered_map<int32_t, int32_t> ww_domain_map;
extern openmc::vector<unique_ptr<WeightWindowDomain>> weight_window_domains;

extern std::unordered_map<int32_t, int32_t> ww_map;
extern openmc::vector<unique_ptr<WeightWindows>> weight_windows;

} // namespace variance_reduction

class WeightWindows {
public:
  // Constructors
  WeightWindows();

  WeightWindows(pugi::xml_node node);

  // Methods
  void to_statepoint(hid_t group) const;

  //! Return a weight window at the specified index
  WeightWindow weight_window(const int index) const;

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
  // Data members
  int32_t id_;
  ParticleType particle_type_;
  openmc::vector<double> energy_bins_;
  openmc::vector<double> lower_ww_;
  openmc::vector<double> upper_ww_;
  double survival_ratio_;
  int max_split_;
  double weight_cutoff_;
};

class WeightWindowDomain {
public:
  // Constructrors - default
  WeightWindowDomain();

  // Constructors
  WeightWindowDomain(pugi::xml_node node);

  // Methods
  bool get_weight_window(Particle& p) const;

  void to_statepoint(hid_t group) const;

  // Accessors
  int32_t id() const { return id_; }
  int32_t mesh_idx() const { return mesh_idx_; }
  int32_t weight_windows_idx() const { return weight_windows_idx_; }

private:
  int32_t id_;        //!< the id of the weight window domain
  int32_t mesh_idx_;  //!< the idx of the mesh this domain uses
  int32_t
    weight_windows_idx_; //!< the idx of the weight windows this domain uses
};

} // namespace openmc
#endif //OPENMC_WEIGHT_WINDOWS_H
