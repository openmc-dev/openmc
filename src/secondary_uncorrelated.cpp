#include "openmc/secondary_uncorrelated.h"

#include <string> // for string

#include <fmt/core.h>

#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/nuclide.h"
#include "openmc/particle.h"
#include "openmc/random_dist.h"
#include "openmc/tallies/tally_scoring.h"

namespace openmc {

//==============================================================================
// UncorrelatedAngleEnergy implementation
//==============================================================================

UncorrelatedAngleEnergy::UncorrelatedAngleEnergy(hid_t group)
{
  // Check if angle group is present & read
  if (object_exists(group, "angle")) {
    hid_t angle_group = open_group(group, "angle");
    angle_ = AngleDistribution {angle_group};
    close_group(angle_group);
  }

  // Check if energy group is present & read
  if (object_exists(group, "energy")) {
    hid_t energy_group = open_group(group, "energy");

    std::string type;
    read_attribute(energy_group, "type", type);
    using UPtrEDist = unique_ptr<EnergyDistribution>;
    if (type == "discrete_photon") {
      energy_ = UPtrEDist {new DiscretePhoton {energy_group}};
    } else if (type == "level") {
      energy_ = UPtrEDist {new LevelInelastic {energy_group}};
    } else if (type == "continuous") {
      energy_ = UPtrEDist {new ContinuousTabular {energy_group}};
    } else if (type == "maxwell") {
      energy_ = UPtrEDist {new MaxwellEnergy {energy_group}};
    } else if (type == "evaporation") {
      energy_ = UPtrEDist {new Evaporation {energy_group}};
    } else if (type == "watt") {
      energy_ = UPtrEDist {new WattEnergy {energy_group}};
    } else {
      warning(
        fmt::format("Energy distribution type '{}' not implemented.", type));
    }
    close_group(energy_group);
  }
}

void UncorrelatedAngleEnergy::sample(
  double E_in, double& E_out, double& mu, uint64_t* seed) const
{
  // Sample cosine of scattering angle
  if (!angle_.empty()) {
    mu = angle_.sample(E_in, seed);
  } else {
    // no angle distribution given => assume isotropic for all energies
    mu = uniform_distribution(-1., 1., seed);
  }

  // Sample outgoing energy
  E_out = energy_->sample(E_in, seed);
}

void UncorrelatedAngleEnergy::get_pdf(double det_pos[4], double E_in,
  double& E_out, uint64_t* seed, Particle& p, std::vector<double>& mu_cm,
  std::vector<double>& Js, std::vector<Particle>& ghost_particles,
  std::vector<double>& pdfs_lab) const
{
  bool COM = false;
  const auto& nuc {data::nuclides[p.event_nuclide()]};
  if (p.event_index_mt() != -999) {
    const auto& rx {nuc->reactions_[p.event_index_mt()]};
    COM = rx->scatter_in_cm_;
  }

  // Sample outgoing energy
  E_out = energy_->sample(E_in, seed);

  if (COM) {

    get_pdf_to_point_elastic(
      det_pos, p, mu_cm, Js, ghost_particles, E_out / 1e6);
    for (std::size_t i = 0; i < mu_cm.size(); ++i) {
      // Assuming Js.size() is the same as mu_cm.size()
      double mu_c = mu_cm[i];
      double derivative = Js[i];
      double pdf_cm;
      if (!angle_.empty()) {
        pdf_cm = angle_.get_pdf(E_in, mu_c, seed);
      } else {
        // no angle distribution given => assume isotropic for all energies
        pdf_cm = 0.5;
      }
      pdfs_lab.push_back(pdf_cm / std::abs(derivative));
    }
  }

  if (!COM) {
    // finding mu_cm, E_out is in lab
    // fatal_error("Didnt implemt lab");
    Direction u_lab {det_pos[0] - p.r().x, // towards the detector
      det_pos[1] - p.r().y, det_pos[2] - p.r().z};
    Direction u_lab_unit = u_lab / u_lab.norm(); // normalize
    double E_lab = E_out;
    Particle ghost_particle = Particle();
    ghost_particle.initialize_ghost_particle(p, u_lab_unit, E_lab);
    ghost_particles.push_back(ghost_particle);
    double pdf_mu_lab;
    if (!angle_.empty()) {
      pdf_mu_lab = angle_.get_pdf(E_in, u_lab_unit.dot(p.u_last()), seed);
    } else {
      // no angle distribution given => assume isotropic for all energies
      pdf_mu_lab = 0.5;
    }

    pdfs_lab.push_back(pdf_mu_lab);
  }
}

} // namespace openmc
