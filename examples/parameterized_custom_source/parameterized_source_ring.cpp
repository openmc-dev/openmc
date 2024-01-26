#include <cmath> // for M_PI
#include <memory> // for unique_ptr
#include <unordered_map>

#include "openmc/random_lcg.h"
#include "openmc/source.h"
#include "openmc/particle.h"

class RingSource : public openmc::Source {
  public:
    RingSource(double radius, double energy) : radius_(radius), energy_(energy) { }

    // Defines a function that can create a unique pointer to a new instance of this class
    // by extracting the parameters from the provided string.
    static std::unique_ptr<RingSource> from_string(std::string parameters)
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

      double radius = std::stod(parameter_mapping["radius"]);
      double energy = std::stod(parameter_mapping["energy"]);
      return std::make_unique<RingSource>(radius, energy);
    }

    // Samples from an instance of this class.
    openmc::SourceSite sample(uint64_t* seed) const
    {
      openmc::SourceSite particle;
      // particle type
      particle.particle = openmc::ParticleType::neutron;
      // position
      double angle = 2.0 * M_PI * openmc::prn(seed);
      double radius = this->radius_;
      particle.r.x = radius * std::cos(angle);
      particle.r.y = radius * std::sin(angle);
      particle.r.z = 0.0;
      // angle
      particle.u = {1.0, 0.0, 0.0};
      particle.E = this->energy_;

      return particle;
    }

  private:
    double radius_;
    double energy_;
};

// A function to create a unique pointer to an instance of this class when generated
// via a plugin call using dlopen/dlsym.
// You must have external C linkage here otherwise dlopen will not find the file
extern "C" std::unique_ptr<RingSource> openmc_create_source(std::string parameters)
{
  return RingSource::from_string(parameters);
}
