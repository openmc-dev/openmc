#include <cmath> // for M_PI
#include <unordered_map>

#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/particle.h"

class SerializedSource {
  protected:
    double radius_;
    double energy_;

    // Protect the constructor so that the class can only be created by serialisation.
    SerializedSource(double radius, double energy) {
      radius_ = radius;
      energy_ = energy;
    }

  public:
    // Getters for the values that we want to use in sampling.
    double radius() { return radius_; }
    double energy() { return energy_; }

    static SerializedSource from_string(const char* parameters) {
      std::unordered_map<std::string, std::string> parameter_mapping;

      std::stringstream ss(parameters);
      std::string parameter;
      while (std::getline(ss, parameter, ',')) {
        parameter.erase(0, parameter.find_first_not_of(' '));
        std::string key = parameter.substr(0, parameter.find_first_of('='));
        std::string value = parameter.substr(parameter.find_first_of('=') + 1, parameter.length());
        parameter_mapping[key] = value;
      }

      return SerializedSource(std::stod(parameter_mapping["radius"]), std::stod(parameter_mapping["energy"]));
    }
};

// you must have external C linkage here otherwise
// dlopen will not find the file
extern "C" openmc::Particle::Bank sample_source(uint64_t* seed, const char* parameters) {
  SerializedSource source = SerializedSource::from_string(parameters);

  openmc::Particle::Bank particle;
  // wgt
  particle.particle = openmc::Particle::Type::neutron;
  particle.wgt = 1.0;
  // position
  double angle = 2. * M_PI * openmc::prn(seed);
  double radius = source.radius();
  particle.r.x = radius * std::cos(angle);
  particle.r.y = radius * std::sin(angle);
  particle.r.z = 0.0;
  // angle
  particle.u = {1.0, 0.0, 0.0};
  particle.E = source.energy();
  particle.delayed_group = 0;

  return particle;
}
