#include "openmc/source.h"

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#define HAS_DYNAMIC_LINKING
#endif

#include <algorithm> // for move

#ifdef HAS_DYNAMIC_LINKING
#include <dlfcn.h> // for dlopen, dlsym, dlclose, dlerror
#endif

#include "xtensor/xadapt.hpp"
#include <fmt/core.h>

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/mcpl_interface.h"
#include "openmc/memory.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/state_point.h"
#include "openmc/string_utils.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {

vector<unique_ptr<Source>> external_sources;
}

//==============================================================================
// Source create implementation
//==============================================================================

unique_ptr<Source> Source::create(pugi::xml_node node)
{
  // if the source type is present, use it to determine the type
  // of object to create
  if (check_for_node(node, "type")) {
    std::string source_type = get_node_value(node, "type");
    if (source_type == "independent") {
      return make_unique<IndependentSource>(node);
    } else if (source_type == "file") {
      return make_unique<FileSource>(node);
    } else if (source_type == "compiled") {
      return make_unique<CompiledSourceWrapper>(node);
    } else if (source_type == "mesh") {
      return make_unique<MeshSource>(node);
    } else {
      fatal_error(fmt::format("Invalid source type '{}' found.", source_type));
    }
  } else {
    // support legacy source format
    if (check_for_node(node, "file")) {
      return make_unique<FileSource>(node);
    } else if (check_for_node(node, "library")) {
      return make_unique<CompiledSourceWrapper>(node);
    } else {
      return make_unique<IndependentSource>(node);
    }
  }
}

//==============================================================================
// IndependentSource implementation
//==============================================================================

IndependentSource::IndependentSource(
  UPtrSpace space, UPtrAngle angle, UPtrDist energy, UPtrDist time)
  : space_ {std::move(space)}, angle_ {std::move(angle)},
    energy_ {std::move(energy)}, time_ {std::move(time)}
{}

IndependentSource::IndependentSource(pugi::xml_node node)
{
  // Check for particle type
  if (check_for_node(node, "particle")) {
    auto temp_str = get_node_value(node, "particle", true, true);
    if (temp_str == "neutron") {
      particle_ = ParticleType::neutron;
    } else if (temp_str == "photon") {
      particle_ = ParticleType::photon;
      settings::photon_transport = true;
    } else {
      fatal_error(std::string("Unknown source particle type: ") + temp_str);
    }
  }

  // Check for source strength
  if (check_for_node(node, "strength")) {
    strength_ = std::stod(get_node_value(node, "strength"));
  }

  // Check for external source file
  if (check_for_node(node, "file")) {

  } else {

    // Spatial distribution for external source
    if (check_for_node(node, "space")) {
      space_ = SpatialDistribution::create(node.child("space"));
    } else {
      // If no spatial distribution specified, make it a point source
      space_ = UPtrSpace {new SpatialPoint()};
    }

    // Determine external source angular distribution
    if (check_for_node(node, "angle")) {
      angle_ = UnitSphereDistribution::create(node.child("angle"));
    } else {
      angle_ = UPtrAngle {new Isotropic()};
    }

    // Determine external source energy distribution
    if (check_for_node(node, "energy")) {
      pugi::xml_node node_dist = node.child("energy");
      energy_ = distribution_from_xml(node_dist);
    } else {
      // Default to a Watt spectrum with parameters 0.988 MeV and 2.249 MeV^-1
      energy_ = UPtrDist {new Watt(0.988e6, 2.249e-6)};
    }

    // Determine external source time distribution
    if (check_for_node(node, "time")) {
      pugi::xml_node node_dist = node.child("time");
      time_ = distribution_from_xml(node_dist);
    } else {
      // Default to a Constant time T=0
      double T[] {0.0};
      double p[] {1.0};
      time_ = UPtrDist {new Discrete {T, p, 1}};
    }

    // Check for domains to reject from
    if (check_for_node(node, "domain_type")) {
      std::string domain_type = get_node_value(node, "domain_type");
      if (domain_type == "cell") {
        domain_type_ = DomainType::CELL;
      } else if (domain_type == "material") {
        domain_type_ = DomainType::MATERIAL;
      } else if (domain_type == "universe") {
        domain_type_ = DomainType::UNIVERSE;
      } else {
        fatal_error(std::string(
          "Unrecognized domain type for source rejection: " + domain_type));
      }

      auto ids = get_node_array<int>(node, "domain_ids");
      domain_ids_.insert(ids.begin(), ids.end());
    }
  }
}

SourceSite IndependentSource::sample(uint64_t* seed) const
{
  SourceSite site;
  site.particle = particle_;

  // Repeat sampling source location until a good site has been found
  bool found = false;
  int n_reject = 0;
  static int n_accept = 0;

  while (!found) {
    // Set particle type
    Particle p;
    p.type() = particle_;
    p.u() = {0.0, 0.0, 1.0};

    // Sample spatial distribution
    p.r() = space_->sample(seed);

    // Now search to see if location exists in geometry
    found = exhaustive_find_cell(p);

    // Check if spatial site is in fissionable material
    if (found) {
      auto space_box = dynamic_cast<SpatialBox*>(space_.get());
      if (space_box) {
        if (space_box->only_fissionable()) {
          // Determine material
          auto mat_index = p.material();
          if (mat_index == MATERIAL_VOID) {
            found = false;
          } else {
            found = model::materials[mat_index]->fissionable();
          }
        }
      }

      // Rejection based on cells/materials/universes
      if (!domain_ids_.empty()) {
        found = false;
        if (domain_type_ == DomainType::MATERIAL) {
          auto mat_index = p.material();
          if (mat_index != MATERIAL_VOID) {
            found = contains(domain_ids_, model::materials[mat_index]->id());
          }
        } else {
          for (int i = 0; i < p.n_coord(); i++) {
            auto id = (domain_type_ == DomainType::CELL)
                        ? model::cells[p.coord(i).cell]->id_
                        : model::universes[p.coord(i).universe]->id_;
            if ((found = contains(domain_ids_, id)))
              break;
          }
        }
      }
    }

    // Check for rejection
    if (!found) {
      ++n_reject;
      if (n_reject >= EXTSRC_REJECT_THRESHOLD &&
          static_cast<double>(n_accept) / n_reject <= EXTSRC_REJECT_FRACTION) {
        fatal_error("More than 95% of external source sites sampled were "
                    "rejected. Please check your external source's spatial "
                    "definition.");
      }
    }

    site.r = p.r();
  }

  // Sample angle
  site.u = angle_->sample(seed);

  // Check for monoenergetic source above maximum particle energy
  auto p = static_cast<int>(particle_);
  auto energy_ptr = dynamic_cast<Discrete*>(energy_.get());
  if (energy_ptr) {
    auto energies = xt::adapt(energy_ptr->x());
    if (xt::any(energies > data::energy_max[p])) {
      fatal_error("Source energy above range of energies of at least "
                  "one cross section table");
    }
  }

  while (true) {
    // Sample energy spectrum
    site.E = energy_->sample(seed);

    // Resample if energy falls above maximum particle energy
    if (site.E < data::energy_max[p])
      break;

    n_reject++;
    if (n_reject >= EXTSRC_REJECT_THRESHOLD &&
        static_cast<double>(n_accept) / n_reject <= EXTSRC_REJECT_FRACTION) {
      fatal_error("More than 95% of external source sites sampled were "
                  "rejected. Please check your external source energy spectrum "
                  "definition.");
    }
  }

  // Sample particle creation time
  site.time = time_->sample(seed);

  // Increment number of accepted samples
  ++n_accept;

  return site;
}

//==============================================================================
// FileSource implementation
//==============================================================================
FileSource::FileSource(pugi::xml_node node)
{
  auto path = get_node_value(node, "file", false, true);
  if (ends_with(path, ".mcpl") || ends_with(path, ".mcpl.gz")) {
    sites_ = mcpl_source_sites(path);
  } else {
    this->load_sites_from_file(path);
  }
}

FileSource::FileSource(const std::string& path)
{
  load_sites_from_file(path);
}

void FileSource::load_sites_from_file(const std::string& path)
{
  // Check if source file exists
  if (!file_exists(path)) {
    fatal_error(fmt::format("Source file '{}' does not exist.", path));
  }

  // Read the source from a binary file instead of sampling from some
  // assumed source distribution
  write_message(6, "Reading source file from {}...", path);

  // Open the binary file
  hid_t file_id = file_open(path, 'r', true);

  // Check to make sure this is a source file
  std::string filetype;
  read_attribute(file_id, "filetype", filetype);
  if (filetype != "source" && filetype != "statepoint") {
    fatal_error("Specified starting source file not a source file type.");
  }

  // Read in the source particles
  read_source_bank(file_id, sites_, false);

  // Close file
  file_close(file_id);
}

SourceSite FileSource::sample(uint64_t* seed) const
{
  size_t i_site = sites_.size() * prn(seed);
  return sites_[i_site];
}

//==============================================================================
// CompiledSourceWrapper implementation
//==============================================================================
CompiledSourceWrapper::CompiledSourceWrapper(pugi::xml_node node)
{
  // Get shared library path and parameters
  auto path = get_node_value(node, "library", false, true);
  std::string parameters;
  if (check_for_node(node, "parameters")) {
    parameters = get_node_value(node, "parameters", false, true);
  }
  setup(path, parameters);
}

void CompiledSourceWrapper::setup(
  const std::string& path, const std::string& parameters)
{
#ifdef HAS_DYNAMIC_LINKING
  // Open the library
  shared_library_ = dlopen(path.c_str(), RTLD_LAZY);
  if (!shared_library_) {
    fatal_error("Couldn't open source library " + path);
  }

  // reset errors
  dlerror();

  // get the function to create the custom source from the library
  auto create_compiled_source = reinterpret_cast<create_compiled_source_t*>(
    dlsym(shared_library_, "openmc_create_source"));

  // check for any dlsym errors
  auto dlsym_error = dlerror();
  if (dlsym_error) {
    std::string error_msg = fmt::format(
      "Couldn't open the openmc_create_source symbol: {}", dlsym_error);
    dlclose(shared_library_);
    fatal_error(error_msg);
  }

  // create a pointer to an instance of the custom source
  compiled_source_ = create_compiled_source(parameters);

#else
  fatal_error("Custom source libraries have not yet been implemented for "
              "non-POSIX systems");
#endif
}

CompiledSourceWrapper::~CompiledSourceWrapper()
{
  // Make sure custom source is cleared before closing shared library
  if (compiled_source_.get())
    compiled_source_.reset();

#ifdef HAS_DYNAMIC_LINKING
  dlclose(shared_library_);
#else
  fatal_error("Custom source libraries have not yet been implemented for "
              "non-POSIX systems");
#endif
}

//==============================================================================
// MeshSource implementation
//==============================================================================

MeshSource::MeshSource(pugi::xml_node node)
{
  int32_t mesh_id = stoi(get_node_value(node, "mesh"));
  int32_t mesh_idx = model::mesh_map.at(mesh_id);
  const auto& mesh = model::meshes[mesh_idx];

  std::vector<double> strengths;
  // read all source distributions and populate strengths vector for MeshSpatial
  // object
  for (auto source_node : node.children("source")) {
    sources_.emplace_back(Source::create(source_node));
    strengths.push_back(sources_.back()->strength());
  }

  // the number of source distributions should either be one or equal to the
  // number of mesh elements
  if (sources_.size() > 1 && sources_.size() != mesh->n_bins()) {
    fatal_error(fmt::format("Incorrect number of source distributions ({}) for "
                            "mesh source with {} elements.",
      sources_.size(), mesh->n_bins()));
  }

  space_ = std::make_unique<MeshSpatial>(mesh_idx, strengths);
}

SourceSite MeshSource::sample(uint64_t* seed) const
{
  // sample location and element from mesh
  auto mesh_location = space_->sample_mesh(seed);

  // Sample source for the chosen element
  int32_t element = mesh_location.first;
  SourceSite site = source(element)->sample(seed);

  // Replace spatial position with the one already sampled
  site.r = mesh_location.second;

  return site;
}

//==============================================================================
// Non-member functions
//==============================================================================

void initialize_source()
{
  write_message("Initializing source particles...", 5);

// Generation source sites from specified distribution in user input
#pragma omp parallel for
  for (int64_t i = 0; i < simulation::work_per_rank; ++i) {
    // initialize random number seed
    int64_t id = simulation::total_gen * settings::n_particles +
                 simulation::work_index[mpi::rank] + i + 1;
    uint64_t seed = init_seed(id, STREAM_SOURCE);

    // sample external source distribution
    simulation::source_bank[i] = sample_external_source(&seed);
  }

  // Write out initial source
  if (settings::write_initial_source) {
    write_message("Writing out initial source...", 5);
    std::string filename = settings::path_output + "initial_source.h5";
    hid_t file_id = file_open(filename, 'w', true);
    write_source_bank(file_id, simulation::source_bank, simulation::work_index);
    file_close(file_id);
  }
}

SourceSite sample_external_source(uint64_t* seed)
{
  // Determine total source strength
  double total_strength = 0.0;
  for (auto& s : model::external_sources)
    total_strength += s->strength();

  // Sample from among multiple source distributions
  int i = 0;
  if (model::external_sources.size() > 1) {
    double xi = prn(seed) * total_strength;
    double c = 0.0;
    for (; i < model::external_sources.size(); ++i) {
      c += model::external_sources[i]->strength();
      if (xi < c)
        break;
    }
  }

  // Sample source site from i-th source distribution
  SourceSite site {model::external_sources[i]->sample(seed)};

  // If running in MG, convert site.E to group
  if (!settings::run_CE) {
    site.E = lower_bound_index(data::mg.rev_energy_bins_.begin(),
      data::mg.rev_energy_bins_.end(), site.E);
    site.E = data::mg.num_energy_groups_ - site.E - 1.;
  }

  return site;
}

void free_memory_source()
{
  model::external_sources.clear();
}

//==============================================================================
// C API
//==============================================================================

extern "C" int openmc_sample_external_source(
  size_t n, uint64_t* seed, void* sites)
{
  if (!sites || !seed) {
    set_errmsg("Received null pointer.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  if (model::external_sources.empty()) {
    set_errmsg("No external sources have been defined.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  auto sites_array = static_cast<SourceSite*>(sites);
  for (size_t i = 0; i < n; ++i) {
    sites_array[i] = sample_external_source(seed);
  }
  return 0;
}

} // namespace openmc
