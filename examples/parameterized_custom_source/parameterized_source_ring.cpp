#include <cmath> // for M_PI
#include <memory> // for unique_ptr
#include <unordered_map>

#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/particle.h"

class Source : public openmc::CustomSource
{
  protected:
    double radius_;
    double energy_;

    // Protect the constructor as we only want the class to be created by the from_string method.
    Source(double radius, double energy)
    {
      radius_ = radius;
      energy_ = energy;
    }

  public:
    // Getters for the values that we want to use in sampling.
    double radius() { return radius_; }
    double energy() { return energy_; }

    // Defines a function that can create a unique pointer to a new instance of this class
    // by extracting the parameters from the provided string.
    static std::unique_ptr<Source> from_string(const char* parameters)
    {
      std::unordered_map<std::string, std::string> parameter_mapping;

      std::stringstream ss(parameters);
      std::string parameter;
      while (std::getline(ss, parameter, ',')) {
        parameter.erase(0, parameter.find_first_not_of(' '));
        std::string key = parameter.substr(0, parameter.find_first_of('='));
        std::string value = parameter.substr(parameter.find_first_of('=') + 1, parameter.length());
        parameter_mapping[key] = value;
      }

      return std::unique_ptr<Source> (
        new Source(std::stod(parameter_mapping["radius"]), std::stod(parameter_mapping["energy"]))
      );
    }

    // Samples from an instance of this class.
    openmc::Particle::Bank sample_source(uint64_t* seed)
    {
      openmc::Particle::Bank particle;
      // wgt
      particle.particle = openmc::Particle::Type::neutron;
      particle.wgt = 1.0;
      // position
      double angle = 2.0 * M_PI * openmc::prn(seed);
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

// A function to create a unique pointer to an instance of this class when generated
// via a plugin call using dlopen/dlsym.
// You must have external C linkage here otherwise dlopen will not find the file
extern "C" std::unique_ptr<Source> openmc_create_source(const char* parameters)
{
  return Source::from_string(parameters);
}
