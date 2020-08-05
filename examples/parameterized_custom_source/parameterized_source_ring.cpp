#include <cmath> // for M_PI
#include <unordered_map>

#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/particle.h"

class ParameterizedSource : public openmc::CustomSource {
  protected:
    double radius_;
    double energy_;

    // Protect the constructor so that the class can only be created by serialisation.
    ParameterizedSource(double radius, double energy) {
      radius_ = radius;
      energy_ = energy;
    }

  public:
    // Getters for the values that we want to use in sampling.
    double radius() { return radius_; }
    double energy() { return energy_; }

    // Defines a function that can create a pointer to a new instance of this class
    // by deserializing from the provided string.
    static ParameterizedSource* from_string(const char* parameters) {
      std::unordered_map<std::string, std::string> parameter_mapping;

      std::stringstream ss(parameters);
      std::string parameter;
      while (std::getline(ss, parameter, ',')) {
        parameter.erase(0, parameter.find_first_not_of(' '));
        std::string key = parameter.substr(0, parameter.find_first_of('='));
        std::string value = parameter.substr(parameter.find_first_of('=') + 1, parameter.length());
        parameter_mapping[key] = value;
      }

      return new ParameterizedSource(std::stod(parameter_mapping["radius"]), std::stod(parameter_mapping["energy"]));
    }

    // Samples from an instance of this class.
    openmc::Particle::Bank sample_source(uint64_t* seed) {
      openmc::Particle::Bank particle;
      // wgt
      particle.particle = openmc::Particle::Type::neutron;
      particle.wgt = 1.0;
      // position
      double angle = 2. * M_PI * openmc::prn(seed);
      double radius = this->radius();
      particle.r.x = radius * std::cos(angle);
      particle.r.y = radius * std::sin(angle);
      particle.r.z = 0.0;
      // angle
      particle.u = {1.0, 0.0, 0.0};
      particle.E = this->energy();
      particle.delayed_group = 0;

      return particle;
    }
};

// A function to create a pointer to an instance of this class when generated
// via a plugin call using dlopen/dlsym.
// You must have external C linkage here otherwise dlopen will not find the file
extern "C" ParameterizedSource* create(const char* parameters) {
  return ParameterizedSource::from_string(parameters);
}
