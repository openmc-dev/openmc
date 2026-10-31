#include "openmc/source.h"

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#define HAS_DYNAMIC_LINKING
#endif

#include <algorithm> // for max
#include <cmath>     // for sin, cos, abs
#include <numeric>   // for accumulate
#include <utility>   // for move

#ifdef HAS_DYNAMIC_LINKING
#include <dlfcn.h> // for dlopen, dlsym, dlclose, dlerror
#endif

#include "openmc/tensor.h"
#include <fmt/core.h>

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/math.h"
#include "openmc/math_functions.h"
#include "openmc/mcpl_interface.h"
#include "openmc/memory.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/state_point.h"
#include "openmc/string_utils.h"
#include "openmc/surface.h"
#include "openmc/xml_interface.h"

namespace openmc {

std::atomic<int64_t> source_n_accept {0};
std::atomic<int64_t> source_n_reject {0};

namespace {

void validate_particle_type(ParticleType type, const std::string& context)
{
  if (type.is_transportable())
    return;

  fatal_error(
    fmt::format("Unsupported source particle type '{}' (PDG {}) in {}.",
      type.str(), type.pdg_number(), context));
}

} // namespace

//==============================================================================
// Global variables
//==============================================================================

namespace model {

vector<unique_ptr<Source>> external_sources;

vector<unique_ptr<Source>> adjoint_sources;

DiscreteIndex external_sources_probability;

} // namespace model

//==============================================================================
// Source implementation
//==============================================================================

Source::Source(pugi::xml_node node)
{
  // Check for source strength
  if (check_for_node(node, "strength")) {
    strength_ = std::stod(get_node_value(node, "strength"));
    if (strength_ < 0.0) {
      fatal_error("Source strength is negative.");
    }
  }

  // Check for additional defined constraints
  read_constraints(node);
}

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
    } else if (source_type == "tokamak") {
      return make_unique<TokamakSource>(node);
    } else if (source_type == "stellarator") {
      return make_unique<StellaratorSource>(node);
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

void Source::read_constraints(pugi::xml_node node)
{
  // Check for constraints node. For backwards compatibility, if no constraints
  // node is given, still try searching for domain constraints from top-level
  // node.
  pugi::xml_node constraints_node = node.child("constraints");
  if (constraints_node) {
    node = constraints_node;
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
      fatal_error(
        std::string("Unrecognized domain type for constraint: " + domain_type));
    }

    auto ids = get_node_array<int>(node, "domain_ids");
    domain_ids_.insert(ids.begin(), ids.end());
  }

  if (check_for_node(node, "time_bounds")) {
    auto ids = get_node_array<double>(node, "time_bounds");
    if (ids.size() != 2) {
      fatal_error("Time bounds must be represented by two numbers.");
    }
    time_bounds_ = std::make_pair(ids[0], ids[1]);
  }
  if (check_for_node(node, "energy_bounds")) {
    auto ids = get_node_array<double>(node, "energy_bounds");
    if (ids.size() != 2) {
      fatal_error("Energy bounds must be represented by two numbers.");
    }
    energy_bounds_ = std::make_pair(ids[0], ids[1]);
  }

  if (check_for_node(node, "fissionable")) {
    only_fissionable_ = get_node_value_bool(node, "fissionable");
  }

  // Check for how to handle rejected particles
  if (check_for_node(node, "rejection_strategy")) {
    std::string rejection_strategy = get_node_value(node, "rejection_strategy");
    if (rejection_strategy == "kill") {
      rejection_strategy_ = RejectionStrategy::KILL;
    } else if (rejection_strategy == "resample") {
      rejection_strategy_ = RejectionStrategy::RESAMPLE;
    } else {
      fatal_error(std::string(
        "Unrecognized strategy source rejection: " + rejection_strategy));
    }
  }
}

void check_rejection_fraction(int64_t n_reject, int64_t n_accept)
{
  // Don't check unless we've hit a minimum number of total sites rejected
  if (n_reject < EXTSRC_REJECT_THRESHOLD)
    return;

  // Compute fraction of accepted sites and compare against minimum
  double fraction = static_cast<double>(n_accept) / n_reject;
  if (fraction <= settings::source_rejection_fraction) {
    fatal_error(fmt::format(
      "Too few source sites satisfied the constraints (minimum source "
      "rejection fraction = {}). Please check your source definition or "
      "set a lower value of Settings.source_rejection_fraction.",
      settings::source_rejection_fraction));
  }
}

SourceSite Source::sample_with_constraints(uint64_t* seed) const
{
  bool accepted = false;
  int64_t n_local_reject = 0;
  SourceSite site {};

  while (!accepted) {
    // Sample a source site without considering constraints yet
    site = this->sample(seed);

    if (constraints_applied()) {
      accepted = true;
    } else {
      // Check whether sampled site satisfies constraints
      accepted = satisfies_spatial_constraints(site.r) &&
                 satisfies_energy_constraints(site.E) &&
                 satisfies_time_constraints(site.time);
      if (!accepted) {
        ++n_local_reject;

        // Check per-particle rejection limit
        if (n_local_reject >= MAX_SOURCE_REJECTIONS_PER_SAMPLE) {
          fatal_error("Exceeded maximum number of source rejections per "
                      "sample. Please check your source definition.");
        }

        // For the "kill" strategy, accept particle but set weight to 0 so that
        // it is terminated immediately
        if (rejection_strategy_ == RejectionStrategy::KILL) {
          accepted = true;
          site.wgt = 0.0;
        }
      }
    }
  }

  // Flush local rejection count, update accept counter, and check overall
  // rejection fraction
  if (n_local_reject > 0) {
    source_n_reject += n_local_reject;
  }
  ++source_n_accept;
  check_rejection_fraction(source_n_reject, source_n_accept);

  return site;
}

bool Source::satisfies_energy_constraints(double E) const
{
  return E > energy_bounds_.first && E < energy_bounds_.second;
}

bool Source::satisfies_time_constraints(double time) const
{
  return time > time_bounds_.first && time < time_bounds_.second;
}

bool Source::satisfies_spatial_constraints(Position r) const
{
  GeometryState geom_state;
  geom_state.r() = r;
  geom_state.u() = {0.0, 0.0, 1.0};

  // Reject particle if it's not in the geometry at all
  bool found = exhaustive_find_cell(geom_state);
  if (!found)
    return false;

  // Check the geometry state against specified domains
  bool accepted = true;
  if (!domain_ids_.empty()) {
    if (domain_type_ == DomainType::MATERIAL) {
      auto mat_index = geom_state.material();
      if (mat_index == MATERIAL_VOID) {
        accepted = false;
      } else {
        accepted = contains(domain_ids_, model::materials[mat_index]->id());
      }
    } else {
      for (int i = 0; i < geom_state.n_coord(); i++) {
        auto id =
          (domain_type_ == DomainType::CELL)
            ? model::cells[geom_state.coord(i).cell()].get()->id_
            : model::universes[geom_state.coord(i).universe()].get()->id_;
        if ((accepted = contains(domain_ids_, id)))
          break;
      }
    }
  }

  // Check if spatial site is in fissionable material
  if (accepted && only_fissionable_) {
    // Determine material
    auto mat_index = geom_state.material();
    if (mat_index == MATERIAL_VOID) {
      accepted = false;
    } else {
      accepted = model::materials[mat_index]->fissionable();
    }
  }

  return accepted;
}

//==============================================================================
// IndependentSource implementation
//==============================================================================

IndependentSource::IndependentSource(
  UPtrSpace space, UPtrAngle angle, UPtrDist energy, UPtrDist time)
  : space_ {std::move(space)}, angle_ {std::move(angle)},
    energy_ {std::move(energy)}, time_ {std::move(time)}
{}

IndependentSource::IndependentSource(pugi::xml_node node) : Source(node)
{
  // Check for particle type
  if (check_for_node(node, "particle")) {
    auto temp_str = get_node_value(node, "particle", false, true);
    particle_ = ParticleType(temp_str);
    if (particle_ == ParticleType::photon() ||
        particle_ == ParticleType::electron() ||
        particle_ == ParticleType::positron()) {
      settings::photon_transport = true;
    }
  }
  validate_particle_type(particle_, "IndependentSource");

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

    // For backwards compatibility, check for only fissionable setting on box
    // source
    auto space_box = dynamic_cast<SpatialBox*>(space_.get());
    if (space_box) {
      if (!only_fissionable_) {
        only_fissionable_ = space_box->only_fissionable();
      }
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

      // For decay photon sources, use the absolute photon emission rate in
      // [photons/s] as the source strength
      if (dynamic_cast<DecaySpectrum*>(energy_.get())) {
        if (strength_ != 1.0) {
          warning(fmt::format(
            "Source strength of {} is ignored because the source uses a "
            "DecaySpectrum energy distribution. The source strength will be "
            "set from the DecaySpectrum emission rate.",
            strength_));
        }
        strength_ = energy_->integral();
      }
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
  }
}

SourceSite IndependentSource::sample(uint64_t* seed) const
{
  SourceSite site {};
  site.particle = particle_;
  double r_wgt = 1.0;
  double E_wgt = 1.0;

  // Repeat sampling source location until a good site has been accepted
  bool accepted = false;
  int64_t n_local_reject = 0;

  while (!accepted) {

    // Sample spatial distribution
    auto [r, r_wgt_temp] = space_->sample(seed);
    site.r = r;
    r_wgt = r_wgt_temp;

    // Check if sampled position satisfies spatial constraints
    accepted = satisfies_spatial_constraints(site.r);

    // Check for rejection
    if (!accepted) {
      ++n_local_reject;
      if (n_local_reject >= MAX_SOURCE_REJECTIONS_PER_SAMPLE) {
        fatal_error("Exceeded maximum number of source rejections per "
                    "sample. Please check your source definition.");
      }
    }
  }

  // Sample angle
  auto [u, u_wgt] = angle_->sample(seed);
  site.u = u;

  site.wgt = r_wgt * u_wgt;

  // Sample energy and time for neutron and photon sources
  if (settings::solver_type != SolverType::RANDOM_RAY) {
    // Check for monoenergetic source above maximum particle energy
    auto p = particle_.transport_index();
    auto energy_ptr = dynamic_cast<Discrete*>(energy_.get());
    auto decay_spectrum = dynamic_cast<DecaySpectrum*>(energy_.get());
    if (energy_ptr) {
      auto energies =
        tensor::Tensor<double>(energy_ptr->x().data(), energy_ptr->x().size());
      if ((energies > data::energy_max[p]).any()) {
        fatal_error("Source energy above range of energies of at least "
                    "one cross section table");
      }
    }

    while (true) {
      // Sample energy spectrum. For decay photon sources, also get the parent
      // nuclide index to store in the source site for tallying purposes.
      if (decay_spectrum) {
        auto sample = decay_spectrum->sample_with_parent(seed);
        site.E = sample.energy;
        E_wgt = sample.weight;
        site.parent_nuclide = sample.parent_nuclide;
      } else {
        auto [E, E_wgt_temp] = energy_->sample(seed);
        site.E = E;
        E_wgt = E_wgt_temp;
      }

      // Resample if energy falls above maximum particle energy
      if (site.E < data::energy_max[p] &&
          (satisfies_energy_constraints(site.E)))
        break;

      ++n_local_reject;
      if (n_local_reject >= MAX_SOURCE_REJECTIONS_PER_SAMPLE) {
        fatal_error("Exceeded maximum number of source rejections per "
                    "sample. Please check your source definition.");
      }
    }

    // Sample particle creation time
    auto [time, time_wgt] = time_->sample(seed);
    site.time = time;

    site.wgt *= (E_wgt * time_wgt);
  }

  // Flush local rejection count into global counter
  if (n_local_reject > 0) {
    source_n_reject += n_local_reject;
  }

  return site;
}

//==============================================================================
// FileSource implementation
//==============================================================================

FileSource::FileSource(pugi::xml_node node) : Source(node)
{
  auto path = get_node_value(node, "file", false, true);
  load_sites_from_file(path);
}

FileSource::FileSource(const std::string& path)
{
  load_sites_from_file(path);
}

void FileSource::load_sites_from_file(const std::string& path)
{
  // If MCPL file, use the dedicated file reader
  if (ends_with(path, ".mcpl") || ends_with(path, ".mcpl.gz")) {
    sites_ = mcpl_source_sites(path);
  } else {
    // Check if source file exists
    if (!file_exists(path)) {
      fatal_error(fmt::format("Source file '{}' does not exist.", path));
    }

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

  // Make sure particles in source file have valid types. If any particle is a
  // photon, electron, or positron, enable photon transport so that the
  // appropriate cross sections are loaded.
  for (const auto& site : this->sites_) {
    validate_particle_type(site.particle, "FileSource");
    if (site.particle == ParticleType::photon() ||
        site.particle == ParticleType::electron() ||
        site.particle == ParticleType::positron()) {
      settings::photon_transport = true;
    }
  }
}

SourceSite FileSource::sample(uint64_t* seed) const
{
  // Sample a particle randomly from list
  size_t i_site = sites_.size() * prn(seed);
  SourceSite site = sites_[i_site];

  // Surface source files store unsigned surface IDs. If the ID refers to a CSG
  // surface containing the source site, determine the signed half-space from
  // the particle direction. Otherwise, ignore the surface ID and allow the
  // normal cell search to locate the particle.
  if (site.surf_id != SURFACE_NONE) {
    auto it = model::surface_map.find(std::abs(site.surf_id));
    if (it != model::surface_map.end()) {
      const auto& surf = *model::surfaces[it->second];
      if (surf.geom_type() == GeometryType::CSG &&
          std::abs(surf.evaluate(site.r)) < FP_COINCIDENT) {
        int surf_id = std::abs(site.surf_id);
        site.surf_id =
          (site.u.dot(surf.normal(site.r)) > 0.0) ? surf_id : -surf_id;
        return site;
      }
    }
    site.surf_id = SURFACE_NONE;
  }

  return site;
}

//==============================================================================
// CompiledSourceWrapper implementation
//==============================================================================

CompiledSourceWrapper::CompiledSourceWrapper(pugi::xml_node node) : Source(node)
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
// MeshElementSpatial implementation
//==============================================================================

std::pair<Position, double> MeshElementSpatial::sample(uint64_t* seed) const
{
  return {model::meshes[mesh_index_]->sample_element(elem_index_, seed), 1.0};
}

//==============================================================================
// MeshSource implementation
//==============================================================================

MeshSource::MeshSource(pugi::xml_node node) : Source(node)
{
  int32_t mesh_id = stoi(get_node_value(node, "mesh"));
  int32_t mesh_idx = model::mesh_map.at(mesh_id);
  const auto& mesh = model::meshes[mesh_idx];

  std::vector<double> strengths;
  // read all source distributions and populate strengths vector for MeshSpatial
  // object
  for (auto source_node : node.children("source")) {
    auto src = Source::create(source_node);
    if (auto ptr = dynamic_cast<IndependentSource*>(src.get())) {
      src.release();
      sources_.emplace_back(ptr);
    } else {
      fatal_error(
        "The source assigned to each element must be an IndependentSource.");
    }
    strengths.push_back(sources_.back()->strength());
  }

  // Set spatial distributions for each mesh element
  for (int elem_index = 0; elem_index < sources_.size(); ++elem_index) {
    sources_[elem_index]->set_space(
      std::make_unique<MeshElementSpatial>(mesh_idx, elem_index));
  }

  // Make sure sources use valid particle types
  for (const auto& src : sources_) {
    validate_particle_type(src->particle_type(), "MeshSource");
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
  // Sample a mesh element based on the relative strengths
  int32_t element = space_->sample_element_index(seed);

  // Sample the distribution for the specific mesh element; note that the
  // spatial distribution has been set for each element using MeshElementSpatial
  return source(element)->sample_with_constraints(seed);
}

//==============================================================================
// TokamakSource implementation
//==============================================================================

TokamakSource::TokamakSource(pugi::xml_node node) : Source(node)
{
  // Read geometry parameters
  major_radius_ = std::stod(get_node_value(node, "major_radius"));
  minor_radius_ = std::stod(get_node_value(node, "minor_radius"));
  elongation_ = std::stod(get_node_value(node, "elongation"));
  triangularity_ = std::stod(get_node_value(node, "triangularity"));
  shafranov_shift_ = std::stod(get_node_value(node, "shafranov_shift"));

  // Read optional vertical shift
  if (check_for_node(node, "vertical_shift")) {
    vertical_shift_ = std::stod(get_node_value(node, "vertical_shift"));
  } else {
    vertical_shift_ = 0.0;
  }

  // Read optional toroidal angle bounds
  if (check_for_node(node, "phi_start")) {
    phi_start_ = std::stod(get_node_value(node, "phi_start"));
  } else {
    phi_start_ = 0.0;
  }
  if (check_for_node(node, "phi_extent")) {
    phi_extent_ = std::stod(get_node_value(node, "phi_extent"));
  } else {
    phi_extent_ = 2.0 * PI;
  }
  if (check_for_node(node, "n_alpha")) {
    n_alpha_ = std::stoi(get_node_value(node, "n_alpha"));
  } else {
    n_alpha_ = 101; // Default
  }

  // Read emission profile
  r_over_a_ = get_node_array<double>(node, "r_over_a");
  emission_density_ = get_node_array<double>(node, "emission_density");

  // Read energy distribution(s)
  for (auto energy_node : node.children("energy")) {
    energy_dists_.push_back(distribution_from_xml(energy_node));
  }

  // Read optional time distribution; default to a delta distribution at t=0
  // for the same behavior as IndependentSource
  if (check_for_node(node, "time")) {
    time_ = distribution_from_xml(node.child("time"));
  } else {
    double T[] {0.0};
    double p[] {1.0};
    time_ = UPtrDist {new Discrete {T, p, 1}};
  }

  // Validate inputs
  if (emission_density_.size() != r_over_a_.size()) {
    fatal_error("TokamakSource: emission_density and r_over_a must have the "
                "same length.");
  }
  if (r_over_a_.size() < 2) {
    fatal_error(
      "TokamakSource: At least 2 radial points are required for profiles.");
  }
  if (r_over_a_.front() != 0.0) {
    fatal_error("TokamakSource: r_over_a must start at 0.");
  }
  if (r_over_a_.back() != 1.0) {
    fatal_error("TokamakSource: r_over_a must end at 1.");
  }
  for (size_t i = 1; i < r_over_a_.size(); ++i) {
    if (r_over_a_[i] <= r_over_a_[i - 1]) {
      fatal_error("TokamakSource: r_over_a must be strictly increasing.");
    }
  }
  for (size_t i = 0; i < emission_density_.size(); ++i) {
    if (emission_density_[i] < 0.0) {
      fatal_error("TokamakSource: emission_density values cannot be negative.");
    }
  }
  if (major_radius_ <= 0.0) {
    fatal_error("TokamakSource: major_radius must be > 0.");
  }
  if (minor_radius_ <= 0.0) {
    fatal_error("TokamakSource: minor_radius must be > 0.");
  }
  if (minor_radius_ >= major_radius_) {
    fatal_error("TokamakSource: minor_radius must be less than major_radius.");
  }
  if (elongation_ <= 0.0) {
    fatal_error("TokamakSource: elongation must be > 0.");
  }
  if (triangularity_ < -1.0 || triangularity_ > 1.0) {
    fatal_error("TokamakSource: triangularity must be in the range [-1, 1].");
  }
  if (shafranov_shift_ < 0.0) {
    fatal_error("TokamakSource: shafranov_shift must be >= 0.");
  }
  if (shafranov_shift_ >= 0.5 * minor_radius_) {
    fatal_error("TokamakSource: shafranov_shift must be less than half the "
                "minor radius.");
  }
  if (phi_extent_ <= 0.0 || phi_extent_ > 2.0 * PI) {
    fatal_error("TokamakSource: phi_extent must be > 0 and <= 2*pi.");
  }
  if (n_alpha_ <= 2) {
    fatal_error("TokamakSource: n_alpha must be > 2.");
  }
  if (n_alpha_ < 51) {
    warning("TokamakSource: n_alpha values below 51 may introduce noticeable "
            "discretization bias in source sampling.");
  }
  if (energy_dists_.empty()) {
    fatal_error("TokamakSource: At least one energy distribution is required.");
  }
  if (energy_dists_.size() != 1 && energy_dists_.size() != r_over_a_.size()) {
    fatal_error("TokamakSource: energy distributions must be either 1 (for all "
                "r) or match the number of r_over_a points.");
  }

  // Compute normalized geometry parameters
  epsilon_ = minor_radius_ / major_radius_;
  delta_tilde_ = shafranov_shift_ / minor_radius_;

  // Initialize isotropic angular distribution
  angle_ = UPtrAngle {new Isotropic()};

  precompute_sampling_distributions();
}

void TokamakSource::precompute_sampling_distributions()
{
  // Use precomputed normalized geometry parameters
  double eps = epsilon_;    // Inverse aspect ratio (a/R0)
  double Dt = delta_tilde_; // Normalized Shafranov shift (Delta/a)
  double delta = triangularity_;

  //==========================================================================
  // RADIAL CDF (computed first since it's simpler and sampled first)
  //==========================================================================
  // The marginal radial PDF is obtained by analytically integrating the joint
  // distribution f(r_tilde, alpha) over alpha. The result is:
  //
  //   p(r_tilde) ~ S(r_tilde) * [(1 + eps*Dt)*r_tilde
  //                              - (3/8)*c1*eps*r_tilde^2
  //                              - 2*eps*Dt*r_tilde^3]
  //
  // where the Bessel function coefficients are:
  //   c0 = J_0(delta) + J_2(delta)
  //   c1 = (J_1(2*delta) + J_3(2*delta)) / c0
  //
  // For delta -> 0, c0 -> 1 and c1 -> 0, giving the circular cross-section
  // limit.

  // Compute Bessel function coefficients. openmc::cyl_bessel_j handles
  // negative arguments (negative triangularity) via the parity relation
  // J_n(-x) = (-1)^n * J_n(x).
  double J0_d = cyl_bessel_j(0, delta);
  double J2_d = cyl_bessel_j(2, delta);
  double J1_2d = cyl_bessel_j(1, 2.0 * delta);
  double J3_2d = cyl_bessel_j(3, 2.0 * delta);
  double c0 = J0_d + J2_d;
  double c1 = (J1_2d + J3_2d) / c0;

  // Coefficients for the radial polynomial: A*r - B*r^2 - C*r^3
  radial_poly_a_ = 1.0 + eps * Dt;
  radial_poly_b_ = 0.375 * c1 * eps; // 3/8 * c1 * eps
  radial_poly_c_ = 2.0 * eps * Dt;

  // Build a refined radial grid that retains the user-specified grid points.
  // The emission density is interpreted as linear-linear between those points.
  constexpr int MIN_SUBINTERVALS = 8;
  constexpr double MAX_GRID_SPACING = 1.0e-3;
  vector<double> radial_grid {r_over_a_.front()};
  vector<double> radial_emission {emission_density_.front()};
  for (size_t i = 1; i < r_over_a_.size(); ++i) {
    double r_lo = r_over_a_[i - 1];
    double r_hi = r_over_a_[i];
    double s_lo = emission_density_[i - 1];
    double s_hi = emission_density_[i];
    int n_subintervals = std::max(MIN_SUBINTERVALS,
      static_cast<int>(std::ceil((r_hi - r_lo) / MAX_GRID_SPACING)));
    for (int j = 1; j <= n_subintervals; ++j) {
      double t = static_cast<double>(j) / n_subintervals;
      radial_grid.push_back(r_lo + t * (r_hi - r_lo));
      radial_emission.push_back(s_lo + t * (s_hi - s_lo));
    }
  }

  vector<double> radial_pdf(radial_grid.size());
  for (size_t i = 0; i < radial_grid.size(); ++i) {
    double r = radial_grid[i];
    // p(r) ~ S(r) * [A*r - B*r^2 - C*r^3]
    double geometric_factor =
      radial_poly_a_ * r - radial_poly_b_ * r * r - radial_poly_c_ * r * r * r;
    radial_pdf[i] = radial_emission[i] * std::max(0.0, geometric_factor);
  }

  // Check that the refined profile contains positive probability mass before
  // constructing the normalized tabular distribution.
  double total = 0.0;
  for (size_t i = 1; i < radial_grid.size(); ++i) {
    total += 0.5 * (radial_pdf[i - 1] + radial_pdf[i]) *
             (radial_grid[i] - radial_grid[i - 1]);
  }
  if (total <= 0.0) {
    fatal_error(
      "TokamakSource: Integrated emission density is zero or negative. "
      "Check emission_density profile.");
  }
  radial_dist_ = make_unique<Tabular>(radial_grid.data(), radial_pdf.data(),
    radial_grid.size(), Interpolation::lin_lin);

  //==========================================================================
  // POLOIDAL CDFs (for conditional sampling of alpha given r)
  //==========================================================================
  // The conditional distribution P(alpha | r) is a mixture:
  //   P(alpha | r) ~ sum_k w_k(r) * I_hat_k * p_k(alpha)
  // where:
  //   - w_k(r) are the "dynamic" Bernstein weight functions (depend on r)
  //   - I_hat_k are the "static" normalized integrals (precomputed constants)
  //   - p_k(alpha) are the normalized basis distributions (precomputed CDFs)
  //
  // The static weights I_hat_k = I_k / (2*pi*c0) are:
  //   I_hat_0 = 1 + eps*Dt
  //   I_hat_1 = 1 + eps*Dt - (3/16)*c1*eps
  //   I_hat_2 = 1 - (3/8)*c1*eps
  //   I_hat_3 = 1 + eps*Dt
  //   I_hat_4 = 1 + (1/2)*eps*Dt - (3/16)*c1*eps
  //   I_hat_5 = 1 - eps*Dt - (3/8)*c1*eps

  // Compute static weights analytically
  poloidal_integrals_[0] = 1.0 + eps * Dt;
  poloidal_integrals_[1] = 1.0 + eps * Dt - 0.1875 * c1 * eps; // 3/16 = 0.1875
  poloidal_integrals_[2] = 1.0 - 0.375 * c1 * eps;             // 3/8 = 0.375
  poloidal_integrals_[3] = 1.0 + eps * Dt;
  poloidal_integrals_[4] = 1.0 + 0.5 * eps * Dt - 0.1875 * c1 * eps;
  poloidal_integrals_[5] = 1.0 - eps * Dt - 0.375 * c1 * eps;

  // Build the alpha grid on [0, pi] (half domain due to up-down symmetry)
  int n_alpha = n_alpha_;
  vector<double> alpha_grid(n_alpha);
  double dalpha = PI / (n_alpha - 1);
  for (int i = 0; i < n_alpha; ++i) {
    alpha_grid[i] = i * dalpha;
  }

  // Compute basis function values g_k(alpha) for tabular distributions
  // Using Bernstein form:
  //   R_tilde = b0*(1-r)^2 + 2*b1*r*(1-r) + b2*r^2
  //   J_tilde = b3*(1-r) + b4*r
  // with:
  //   b0(alpha) = 1 + eps*Dt
  //   b1(alpha) = b0 + (eps/2)*cos(psi),  psi = alpha + delta*sin(alpha)
  //   b2(alpha) = 1 + eps*cos(psi)
  //   b3(alpha) = cos(delta*sin(alpha))
  //               + (delta/4)*(cos(alpha - delta*sin(alpha))
  //                          - cos(3*alpha + delta*sin(alpha)))
  //   b4(alpha) = b3(alpha) - 2*Dt*cos(alpha)

  array<vector<double>, N_POLOIDAL_BASIS> basis;
  for (int k = 0; k < N_POLOIDAL_BASIS; ++k) {
    basis[k].resize(n_alpha);
  }

  for (int i = 0; i < n_alpha; ++i) {
    double alpha = alpha_grid[i];
    double sin_alpha = openmc::sin(alpha);
    double cos_alpha = openmc::cos(alpha);
    double delta_sin_alpha = delta * sin_alpha;
    double psi = alpha + delta_sin_alpha;
    double cos_psi = openmc::cos(psi);

    // Bernstein coefficients b0-b4
    double b0 = 1.0 + eps * Dt;
    double b1 = b0 + 0.5 * eps * cos_psi;
    double b2 = 1.0 + eps * cos_psi;
    double b3 = openmc::cos(delta_sin_alpha) +
                0.25 * delta *
                  (openmc::cos(alpha - delta_sin_alpha) -
                    openmc::cos(3.0 * alpha + delta_sin_alpha));
    double b4 = b3 - 2.0 * Dt * cos_alpha;

    // 6 basis functions g_k(alpha) = b_i * b_j
    basis[0][i] = b0 * b3; // w0 = (1-r)^3
    basis[1][i] = b1 * b3; // w1 = 2*r*(1-r)^2
    basis[2][i] = b2 * b3; // w2 = r^2*(1-r)
    basis[3][i] = b0 * b4; // w3 = r*(1-r)^2
    basis[4][i] = b1 * b4; // w4 = 2*r^2*(1-r)
    basis[5][i] = b2 * b4; // w5 = r^3
  }

  // Build a linear-linear distribution for each basis function p_k(alpha)
  for (int k = 0; k < N_POLOIDAL_BASIS; ++k) {
    poloidal_dists_[k] = make_unique<Tabular>(
      alpha_grid.data(), basis[k].data(), n_alpha, Interpolation::lin_lin);
  }
}

double TokamakSource::sample_r_over_a(uint64_t* seed) const
{
  return radial_dist_->sample(seed).first;
}

double TokamakSource::mixture_weight(int k, double r) const
{
  double s = 1.0 - r;
  switch (k) {
  case 0:
    return s * s * s * poloidal_integrals_[0];
  case 1:
    return 2.0 * r * s * s * poloidal_integrals_[1];
  case 2:
    return r * r * s * poloidal_integrals_[2];
  case 3:
    return r * s * s * poloidal_integrals_[3];
  case 4:
    return 2.0 * r * r * s * poloidal_integrals_[4];
  case 5:
    return r * r * r * poloidal_integrals_[5];
  default:
    UNREACHABLE();
  }
}

double TokamakSource::sample_poloidal_angle(double r_norm, uint64_t* seed) const
{
  // Sample from the conditional distribution P(alpha | r_tilde) using
  // mixture sampling with 6 precomputed basis distributions.
  //
  // The conditional is: P(alpha | r) ~ sum_k w_k(r) * I_hat_k * p_k(alpha)
  // where:
  //   - w_k(r) are the "dynamic" Bernstein weight functions
  //   - I_hat_k are the "static" normalized integrals (precomputed in
  //   poloidal_integrals_)
  //   - p_k(alpha) are the normalized, precomputed basis distributions
  //
  // The normalization sum_k w_k(r) * I_hat_k equals the radial geometric
  // polynomial evaluated at r, which is known analytically.
  //
  // Algorithm:
  // 1. Compute total from analytical normalization
  // 2. Lazily evaluate mixture weights with early exit to select component k
  // 3. Sample alpha from the selected basis distribution

  // Analytical normalization: sum_k w_k(r) * I_hat_k
  double total =
    radial_poly_a_ - radial_poly_b_ * r_norm - radial_poly_c_ * r_norm * r_norm;
  double xi = prn(seed) * total;

  // Sample component via lazy evaluation with early exit
  // Order optimized for peaked emission profiles: 0, 1, 4, 5, 3, 2
  constexpr int order[] = {0, 1, 4, 5, 3, 2};
  double cumsum = 0.0;
  int component = order[N_POLOIDAL_BASIS - 1];
  for (int i = 0; i < N_POLOIDAL_BASIS; ++i) {
    cumsum += mixture_weight(order[i], r_norm);
    if (xi < cumsum) {
      component = order[i];
      break;
    }
  }

  // Sample alpha from [0, pi]
  double alpha = poloidal_dists_[component]->sample(seed).first;

  // Exploit up-down symmetry: randomly flip to [pi, 2*pi] with 50% probability
  // This is equivalent to flipping the sign of Z in the final position
  if (prn(seed) >= 0.5) {
    alpha = 2.0 * PI - alpha;
  }
  return alpha;
}

std::pair<double, double> TokamakSource::sample_energy(
  double r_norm, uint64_t* seed) const
{
  if (energy_dists_.size() == 1) {
    // Single distribution for all r
    return energy_dists_[0]->sample(seed);
  }

  // Multiple distributions: stochastic selection between bracketing r points
  // Find the interval containing r_norm
  size_t i = lower_bound_index(r_over_a_.begin(), r_over_a_.end(), r_norm);

  // Handle boundary cases
  if (i >= energy_dists_.size() - 1) {
    return energy_dists_.back()->sample(seed);
  }

  // Stochastic interpolation: randomly select one of the two bracketing
  // distributions based on distance to each
  double t = (r_norm - r_over_a_[i]) / (r_over_a_[i + 1] - r_over_a_[i]);
  size_t idx = (prn(seed) < t) ? i + 1 : i;
  return energy_dists_[idx]->sample(seed);
}

Position TokamakSource::flux_to_cartesian(
  double r, double alpha, double phi) const
{
  // Flux surface parameterization:
  // R = R0 + r*cos(alpha + delta*sin(alpha)) + Delta*(1 - (r/a)^2)
  // Z = kappa * r * sin(alpha)
  // x = R * cos(phi)
  // y = R * sin(phi)
  // z = Z

  double psi = alpha + triangularity_ * openmc::sin(alpha);
  double r_over_a_sq = (r * r) / (minor_radius_ * minor_radius_);

  double R = major_radius_ + r * openmc::cos(psi) +
             shafranov_shift_ * (1.0 - r_over_a_sq);
  double Z = elongation_ * r * openmc::sin(alpha);

  double x = R * openmc::cos(phi);
  double y = R * openmc::sin(phi);
  double z = Z;

  return {x, y, z};
}

SourceSite TokamakSource::sample(uint64_t* seed) const
{
  SourceSite site;
  site.particle = ParticleType::neutron();
  site.wgt = 1.0;
  site.delayed_group = 0;

  // 1. Sample r/a from radial CDF
  double r_norm = sample_r_over_a(seed);
  double r = r_norm * minor_radius_;

  // 2. Sample poloidal angle from conditional distribution P(alpha|r)
  double alpha = sample_poloidal_angle(r_norm, seed);

  // 3. Sample toroidal angle uniformly in [phi_start, phi_start + phi_extent]
  double phi = phi_start_ + phi_extent_ * prn(seed);

  // 4. Convert to Cartesian coordinates
  site.r = flux_to_cartesian(r, alpha, phi);

  // 4a. Apply vertical shift if non-zero
  if (vertical_shift_ != 0.0) {
    site.r.z += vertical_shift_;
  }

  // 5. Sample isotropic direction
  site.u = angle_->sample(seed).first;

  // 6. Sample energy from distribution(s), applying the importance weight so
  // that biased distributions are handled correctly
  auto [E, E_wgt] = sample_energy(r_norm, seed);
  site.E = E;

  // 7. Sample particle creation time
  auto [time, time_wgt] = time_->sample(seed);
  site.time = time;

  site.wgt *= E_wgt * time_wgt;

  return site;
}

//==============================================================================
// StellaratorSource implementation
//==============================================================================

namespace {

//! Maximum over t in [0, 1] of the quadratic through q(0)=f0, q(1/2)=fh,
//! q(1)=f1. Exact, since the sampled density is quadratic in the radial
//! interpolation parameter within a bin.
double quad_max01(double f0, double fh, double f1)
{
  double c = 2.0 * (f0 + f1 - 2.0 * fh);
  double b = f1 - f0 - c;
  double m = std::max(f0, f1);
  if (c < 0.0) {
    double t_star = -b / (2.0 * c);
    if (t_star > 0.0 && t_star < 1.0) {
      m = std::max(m, f0 + t_star * (b + t_star * c));
    }
  }
  return m;
}

// Bernstein basis functions are nonnegative and sum to one on [0,1], so the
// largest Bernstein coefficient bounds the cubic for rejection sampling.
double cubic_rejection_bound(const array<double, 4>& c)
{
  double bound = std::max({c[0], c[0] + c[1] / 3,
    c[0] + 2 * c[1] / 3 + c[2] / 3, c[0] + c[1] + c[2] + c[3]});
  // Add a small relative margin for rounding in the bound and PDF evaluation.
  return bound * (1.0 + 1.e-12);
}

} // namespace

StellaratorSource::StellaratorSource(pugi::xml_node node) : Source(node)
{
  // Read number of field periods
  if (check_for_node(node, "num_field_periods")) {
    num_field_periods_ = std::stoi(get_node_value(node, "num_field_periods"));
  } else {
    num_field_periods_ = 1;
  }

  // Read radial grid and emission profile
  rho_ = get_node_array<double>(node, "rho");
  emission_density_ = get_node_array<double>(node, "emission_density");

  // Read Fourier mode numbers and coefficient tables
  mode_m_ = get_node_array<int>(node, "mode_m");
  mode_n_ = get_node_array<int>(node, "mode_n");
  rmnc_ = get_node_array<double>(node, "rmnc");
  zmns_ = get_node_array<double>(node, "zmns");
  if (check_for_node(node, "rmns") || check_for_node(node, "zmnc")) {
    if (!check_for_node(node, "rmns") || !check_for_node(node, "zmnc")) {
      fatal_error("StellaratorSource: rmns and zmnc must both be given for "
                  "non-stellarator-symmetric equilibria.");
    }
    rmns_ = get_node_array<double>(node, "rmns");
    zmnc_ = get_node_array<double>(node, "zmnc");
    asym_ = true;
  }

  // Read energy distribution(s)
  for (auto energy_node : node.children("energy")) {
    energy_dists_.push_back(distribution_from_xml(energy_node));
  }

  // Read optional time distribution; default to a delta distribution at t=0
  if (check_for_node(node, "time")) {
    time_ = distribution_from_xml(node.child("time"));
  } else {
    double T[] {0.0};
    double p[] {1.0};
    time_ = UPtrDist {new Discrete {T, p, 1}};
  }

  // Validate inputs
  size_t n_rho = rho_.size();
  size_t n_modes = mode_m_.size();
  if (num_field_periods_ < 1) {
    fatal_error("StellaratorSource: num_field_periods must be >= 1.");
  }
  if (n_rho < 2) {
    fatal_error(
      "StellaratorSource: At least 2 radial points are required for profiles.");
  }
  if (rho_.front() != 0.0) {
    fatal_error("StellaratorSource: rho must start at 0.");
  }
  if (rho_.back() != 1.0) {
    fatal_error("StellaratorSource: rho must end at 1.");
  }
  for (size_t i = 1; i < n_rho; ++i) {
    if (rho_[i] <= rho_[i - 1]) {
      fatal_error("StellaratorSource: rho must be strictly increasing.");
    }
  }
  if (emission_density_.size() != n_rho) {
    fatal_error("StellaratorSource: emission_density and rho must have the "
                "same length.");
  }
  for (double s : emission_density_) {
    if (s < 0.0) {
      fatal_error(
        "StellaratorSource: emission_density values cannot be negative.");
    }
  }
  if (std::none_of(emission_density_.begin(), emission_density_.end(),
        [](double s) { return s > 0.0; })) {
    fatal_error("StellaratorSource: emission_density must contain at least "
                "one positive value.");
  }
  if (n_modes == 0) {
    fatal_error("StellaratorSource: At least one Fourier mode is required.");
  }
  if (mode_n_.size() != n_modes) {
    fatal_error(
      "StellaratorSource: mode_m and mode_n must have the same length.");
  }
  for (int m : mode_m_) {
    if (m < 0) {
      fatal_error("StellaratorSource: poloidal mode numbers must be >= 0 in "
                  "the combined Fourier form cos/sin(m*theta - n*Nfp*zeta).");
    }
  }
  if (rmnc_.size() != n_rho * n_modes || zmns_.size() != n_rho * n_modes) {
    fatal_error("StellaratorSource: rmnc and zmns must have length "
                "len(rho)*len(mode_m).");
  }
  if (asym_ &&
      (rmns_.size() != n_rho * n_modes || zmnc_.size() != n_rho * n_modes)) {
    fatal_error("StellaratorSource: rmns and zmnc must have length "
                "len(rho)*len(mode_m).");
  }
  if (energy_dists_.empty()) {
    fatal_error(
      "StellaratorSource: At least one energy distribution is required.");
  }
  if (energy_dists_.size() != 1 && energy_dists_.size() != n_rho) {
    fatal_error("StellaratorSource: energy distributions must be either 1 "
                "(for all rho) or match the number of rho points.");
  }

  // Initialize isotropic angular distribution
  angle_ = UPtrAngle {new Isotropic()};

  precompute_sampling_distributions();
}

void StellaratorSource::precompute_sampling_distributions()
{
  size_t n_rho = rho_.size();
  size_t n_modes = mode_m_.size();
  size_t n_bins = n_rho - 1;

  // Angular quadrature grid. The density R*tau is a trigonometric polynomial
  // of poloidal degree <= 3*m_max and toroidal degree <= 3*n_max per field
  // period, so a uniform trapezoidal product rule with more points than the
  // bandwidth integrates it exactly (discrete Fourier orthogonality). The
  // denser grid below also estimates the angular maximum for rejection
  // sampling; a heuristic safety factor allows for peaks between grid points.
  int m_max = 0;
  int n_max = 0;
  for (size_t k = 0; k < n_modes; ++k) {
    m_max = std::max(m_max, mode_m_[k]);
    n_max = std::max(n_max, std::abs(mode_n_[k]));
  }
  int n_theta = std::max(64, 8 * m_max);
  int n_zeta = std::max(64, 8 * n_max);
  size_t n_grid = static_cast<size_t>(n_theta) * n_zeta;

  // Per-mode trigonometric tables on the angle grid (zeta spans one field
  // period; the geometry is periodic with period 2*pi/Nfp)
  double dtheta = 2.0 * PI / n_theta;
  double dzeta = 2.0 * PI / (num_field_periods_ * n_zeta);
  vector<double> cos_mt(n_modes * n_theta), sin_mt(n_modes * n_theta);
  vector<double> cos_nz(n_modes * n_zeta), sin_nz(n_modes * n_zeta);
  for (size_t k = 0; k < n_modes; ++k) {
    for (int i = 0; i < n_theta; ++i) {
      double a = mode_m_[k] * i * dtheta;
      cos_mt[k * n_theta + i] = std::cos(a);
      sin_mt[k * n_theta + i] = std::sin(a);
    }
    for (int j = 0; j < n_zeta; ++j) {
      double a = mode_n_[k] * num_field_periods_ * j * dzeta;
      cos_nz[k * n_zeta + j] = std::cos(a);
      sin_nz[k * n_zeta + j] = std::sin(a);
    }
  }

  // Keep only adjacent surfaces, reusing the older surface's storage after
  // processing each interval. Each surface is evaluated once.
  vector<double> R_s(2 * n_grid), Z_s(2 * n_grid);
  vector<double> Rt_s(2 * n_grid), Zt_s(2 * n_grid);

  // Evaluate positions and poloidal derivatives at each angular grid point.
  auto evaluate_surface = [&](size_t i) {
    const double* rc = &rmnc_[i * n_modes];
    const double* zs = &zmns_[i * n_modes];
    const double* rs = asym_ ? &rmns_[i * n_modes] : nullptr;
    const double* zc = asym_ ? &zmnc_[i * n_modes] : nullptr;
    for (int it = 0; it < n_theta; ++it) {
      for (int jz = 0; jz < n_zeta; ++jz) {
        double R = 0.0, Z = 0.0, Rt = 0.0, Zt = 0.0;
        for (size_t k = 0; k < n_modes; ++k) {
          // cos/sin(m*theta - n*Nfp*zeta) via angle-addition from tables
          double cm = cos_mt[k * n_theta + it];
          double sm = sin_mt[k * n_theta + it];
          double cn = cos_nz[k * n_zeta + jz];
          double sn = sin_nz[k * n_zeta + jz];
          double c = cm * cn + sm * sn;
          double s = sm * cn - cm * sn;
          double m = mode_m_[k];
          R += rc[k] * c;
          Rt -= m * rc[k] * s;
          Z += zs[k] * s;
          Zt += m * zs[k] * c;
          if (asym_) {
            R += rs[k] * s;
            Rt += m * rs[k] * c;
            Z += zc[k] * c;
            Zt -= m * zc[k] * s;
          }
        }
        if (R <= 0.0) {
          fatal_error(
            "StellaratorSource: flux surfaces reach R <= 0; check the "
            "Fourier coefficients (units should be cm).");
        }
        size_t p = (i % 2) * n_grid + it * n_zeta + jz;
        R_s[p] = R;
        Z_s[p] = Z;
        Rt_s[p] = Rt;
        Zt_s[p] = Zt;
      }
    }
  };

  // Accumulate signed integrals and bounds for both orientations. The total
  // volume determines the Jacobian sign after all intervals are processed.
  double quad_wgt = 4.0 * PI * PI / n_grid; // full-torus angular weight
  double v_tot = 0.0;
  vector<double> vp_lo(n_bins), vp_mid(n_bins), vp_hi(n_bins);
  vector<double> negative_envelope(n_bins);
  envelope_.assign(n_bins, 0.0);
  double f_peak = 0.0;
  double f_min = 0.0;
  evaluate_surface(0);
  for (size_t b = 0; b < n_bins; ++b) {
    evaluate_surface(b + 1);
    double inv_drho = 1.0 / (rho_[b + 1] - rho_[b]);
    size_t lo = (b % 2) * n_grid;
    size_t hi = ((b + 1) % 2) * n_grid;
    const double* R0 = &R_s[lo];
    const double* R1 = &R_s[hi];
    const double* Z0 = &Z_s[lo];
    const double* Z1 = &Z_s[hi];
    const double* Rt0 = &Rt_s[lo];
    const double* Rt1 = &Rt_s[hi];
    const double* Zt0 = &Zt_s[lo];
    const double* Zt1 = &Zt_s[hi];
    double s0 = 0.0, sh = 0.0, s1 = 0.0;
    double env = 0.0;
    double negative_env = 0.0;
    for (size_t p = 0; p < n_grid; ++p) {
      // Linear interpolation makes radial derivatives constant within a bin.
      // R and the poloidal derivatives are linear, so R*tau is quadratic.
      double Rr = (R1[p] - R0[p]) * inv_drho;
      double Zr = (Z1[p] - Z0[p]) * inv_drho;
      double f0 = R0[p] * (Rr * Zt0[p] - Rt0[p] * Zr);
      double f1 = R1[p] * (Rr * Zt1[p] - Rt1[p] * Zr);
      double Rh = 0.5 * (R0[p] + R1[p]);
      double Rth = 0.5 * (Rt0[p] + Rt1[p]);
      double Zth = 0.5 * (Zt0[p] + Zt1[p]);
      double fh = Rh * (Rr * Zth - Rth * Zr);
      s0 += f0;
      sh += fh;
      s1 += f1;

      // Endpoints and midpoint determine the radial quadratic at this angle.
      // Bound both signs until the global Jacobian orientation is known.
      env = std::max(env, quad_max01(f0, fh, f1));
      negative_env = std::max(negative_env, quad_max01(-f0, -fh, -f1));
      f_peak = std::max({f_peak, f0, fh, f1});
      f_min = std::min({f_min, f0, fh, f1});
    }
    vp_lo[b] = s0 * quad_wgt;
    vp_mid[b] = sh * quad_wgt;
    vp_hi[b] = s1 * quad_wgt;
    v_tot += sh * quad_wgt * (rho_[b + 1] - rho_[b]);
    // Heuristic allowance for angular peaks between quadrature points
    envelope_[b] = 1.05 * env;
    negative_envelope[b] = 1.05 * negative_env;
  }
  jacobian_sign_ = (v_tot >= 0.0) ? 1.0 : -1.0;
  if (jacobian_sign_ < 0.0) {
    // Reverse the integrals and extrema to use a positive volume element.
    envelope_.swap(negative_envelope);
    double old_peak = f_peak;
    f_peak = -f_min;
    f_min = -old_peak;
    for (size_t b = 0; b < n_bins; ++b) {
      vp_lo[b] = -vp_lo[b];
      vp_mid[b] = -vp_mid[b];
      vp_hi[b] = -vp_hi[b];
    }
  }
  if (f_min < -1.0e-6 * f_peak) {
    fatal_error(
      "StellaratorSource: the flux-coordinate Jacobian changes sign within "
      "the plasma volume; the flux surfaces self-intersect or the Fourier "
      "coefficient tables are inconsistent.");
  }

  // Integrate S(rho)*V'(rho) independently on each radial interval. Scaling
  // the emission profile by its maximum preserves the normalized distribution.
  double s_max =
    *std::max_element(emission_density_.begin(), emission_density_.end());
  vector<double> masses(n_bins);
  radial_pdf_.resize(n_bins);
  for (size_t b = 0; b < n_bins; ++b) {
    double h = rho_[b + 1] - rho_[b];
    double v0 = vp_lo[b];
    double vh = vp_mid[b];
    double v1 = vp_hi[b];

    // Recover V'(t) = a[0] + a[1]*t + a[2]*t^2 for t in [0, 1]. Each
    // interval retains its own one-sided radial derivatives at the endpoints.
    array<double, 3> a {v0, 4 * vh - 3 * v0 - v1, 2 * (v0 + v1 - 2 * vh)};

    // A convex quadratic can have its minimum inside the interval.
    double minimum = std::min(v0, v1);
    if (a[2] > 0) {
      double t = -a[1] / (2 * a[2]);
      if (t > 0 && t < 1)
        minimum = std::min(minimum, a[0] + t * (a[1] + t * a[2]));
    }
    if (!std::isfinite(v0 + vh + v1) || minimum < 0)
      fatal_error("StellaratorSource: invalid differential volume polynomial.");
    double s0 = emission_density_[b] / s_max;
    double s1 = emission_density_[b + 1] / s_max;
    double ds = s1 - s0;
    // Simpson's rule exactly integrates the cubic S(rho)*V'(rho).
    masses[b] = h * ((s0 / 6) * v0 + ((s0 + s1) / 3) * vh + (s1 / 6) * v1);

    // Multiply V'(t) by the linear emission profile to form the radial PDF.
    auto& c = radial_pdf_[b];
    c = {s0 * a[0], s0 * a[1] + ds * a[0], s0 * a[2] + ds * a[1], ds * a[2]};
    double bound = cubic_rejection_bound(c);
    if (bound > 0) {
      if (!std::isfinite(bound) || !std::isfinite(masses[b]) || masses[b] <= 0)
        fatal_error("StellaratorSource: invalid radial interval mass.");
      // Store the acceptance probability for uniform proposals within the bin.
      for (double& v : c)
        v /= bound;
    }
  }
  // Select intervals by integrated probability, including their radial widths.
  radial_bins_.assign({masses.data(), masses.size()});
}

double StellaratorSource::eval_density(int bin, double t, double theta,
  double zeta, double* R_out, double* Z_out) const
{
  size_t n_modes = mode_m_.size();
  // Coefficients at the two surfaces bracketing the sampled radial position
  const double* rc0 = &rmnc_[bin * n_modes];
  const double* rc1 = &rmnc_[(bin + 1) * n_modes];
  const double* zs0 = &zmns_[bin * n_modes];
  const double* zs1 = &zmns_[(bin + 1) * n_modes];
  const double* rs0 = asym_ ? &rmns_[bin * n_modes] : nullptr;
  const double* rs1 = asym_ ? &rmns_[(bin + 1) * n_modes] : nullptr;
  const double* zc0 = asym_ ? &zmnc_[bin * n_modes] : nullptr;
  const double* zc1 = asym_ ? &zmnc_[(bin + 1) * n_modes] : nullptr;
  double inv_drho = 1.0 / (rho_[bin + 1] - rho_[bin]);

  // Accumulate positions and derivatives; suffixes r and t denote rho and
  // theta.
  double R = 0.0, Rr = 0.0, Rt = 0.0;
  double Z = 0.0, Zr = 0.0, Zt = 0.0;
  for (size_t k = 0; k < n_modes; ++k) {
    double a = mode_m_[k] * theta - mode_n_[k] * num_field_periods_ * zeta;
    double c = std::cos(a);
    double s = std::sin(a);
    double m = mode_m_[k];

    // Interpolate coefficients in t; their rho derivatives are the surface
    // differences divided by the bin width. Angular derivatives act on the
    // Fourier basis and introduce a factor of m.
    double drc = rc1[k] - rc0[k];
    double rc = rc0[k] + t * drc;
    double dzs = zs1[k] - zs0[k];
    double zs = zs0[k] + t * dzs;
    R += rc * c;
    Rr += drc * inv_drho * c;
    Rt -= m * rc * s;
    Z += zs * s;
    Zr += dzs * inv_drho * s;
    Zt += m * zs * c;
    if (asym_) {
      // Add the sine-R and cosine-Z terms for non-stellarator-symmetric shapes.
      double drs = rs1[k] - rs0[k];
      double rs = rs0[k] + t * drs;
      double dzc = zc1[k] - zc0[k];
      double zc = zc0[k] + t * dzc;
      R += rs * s;
      Rr += drs * inv_drho * s;
      Rt += m * rs * c;
      Z += zc * c;
      Zr += dzc * inv_drho * c;
      Zt -= m * zc * s;
    }
  }
  // Reuse the evaluated position when an angular proposal is accepted.
  if (R_out)
    *R_out = R;
  if (Z_out)
    *Z_out = Z;

  // Since zeta is the cylindrical azimuth, the volume Jacobian is R*tau,
  // with tau = R_rho*Z_theta - R_theta*Z_rho. S(rho) cancels from the
  // conditional angular density; the global sign accounts for handedness.
  return jacobian_sign_ * R * (Rr * Zt - Rt * Zr);
}

std::pair<double, double> StellaratorSource::sample_energy(
  double rho, uint64_t* seed) const
{
  if (energy_dists_.size() == 1) {
    // Single distribution for all rho
    return energy_dists_[0]->sample(seed);
  }

  // Multiple distributions: stochastic selection between bracketing points
  size_t i = lower_bound_index(rho_.begin(), rho_.end(), rho);
  if (i >= energy_dists_.size() - 1) {
    return energy_dists_.back()->sample(seed);
  }
  double t = (rho - rho_[i]) / (rho_[i + 1] - rho_[i]);
  size_t idx = (prn(seed) < t) ? i + 1 : i;
  return energy_dists_[idx]->sample(seed);
}

SourceSite StellaratorSource::sample(uint64_t* seed) const
{
  SourceSite site;

  // 1. Select a radial bin by its integrated mass
  int bin = radial_bins_.sample(seed);
  const auto& c = radial_pdf_[bin];

  // 2. Sample uniformly within the radial bin and reject based on the PDF
  double t = prn(seed);
  int64_t n_reject = 0;
  while (true) {
    double pdf = c[0] + t * (c[1] + t * (c[2] + t * c[3]));
    if (prn(seed) < pdf)
      break;
    if (++n_reject > MAX_SOURCE_REJECTIONS_PER_SAMPLE) {
      fatal_error("StellaratorSource: exceeded the maximum number of "
                  "rejections while sampling within a radial interval.");
    }
    // Retain the selected bin so rejection does not change its probability.
    t = prn(seed);
  }
  double rho = rho_[bin] + t * (rho_[bin + 1] - rho_[bin]);

  // 3. Rejection-sample (theta, zeta) from p(theta, zeta | rho) ~ R*tau
  double env = envelope_[bin];
  double R, Z, zeta;
  n_reject = 0;
  while (true) {
    double theta = uniform_distribution(0.0, 2.0 * PI, seed);
    zeta = uniform_distribution(0.0, 2.0 * PI, seed);
    double f = eval_density(bin, t, theta, zeta, &R, &Z);
    if (prn(seed) * env < f)
      break;
    if (++n_reject > MAX_SOURCE_REJECTIONS_PER_SAMPLE) {
      fatal_error("StellaratorSource: exceeded the maximum number of "
                  "rejections while sampling the poloidal/toroidal angles.");
    }
  }

  // 4. Convert to Cartesian coordinates (zeta is the cylindrical angle)
  site.r = {R * std::cos(zeta), R * std::sin(zeta), Z};

  // 5. Sample isotropic direction
  site.u = angle_->sample(seed).first;

  // 6. Sample energy from distribution(s), applying the importance weight
  auto [E, E_wgt] = sample_energy(rho, seed);
  site.E = E;

  // 7. Sample particle creation time
  auto [time, time_wgt] = time_->sample(seed);
  site.time = time;

  site.wgt *= E_wgt * time_wgt;

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
  // Sample from among multiple source distributions
  int i = 0;
  int n_sources = model::external_sources.size();
  if (n_sources > 1) {
    if (settings::uniform_source_sampling) {
      i = prn(seed) * n_sources;
    } else {
      i = model::external_sources_probability.sample(seed);
    }
  }

  // Sample source site from i-th source distribution
  SourceSite site {model::external_sources[i]->sample_with_constraints(seed)};

  // For uniform source sampling, multiply the weight by the ratio of the actual
  // probability of sampling source i to the biased probability of sampling
  // source i, which is (strength_i / total_strength) / (1 / n)
  if (n_sources > 1 && settings::uniform_source_sampling) {
    double total_strength = model::external_sources_probability.integral();
    site.wgt *=
      model::external_sources[i]->strength() * n_sources / total_strength;
  }

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
  model::adjoint_sources.clear();
  reset_source_rejection_counters();
}

void reset_source_rejection_counters()
{
  source_n_accept = 0;
  source_n_reject = 0;
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

  // Derive independent per-particle seeds from the base seed so that
  // each iteration has its own RNG state for thread-safe parallel sampling.
  uint64_t base_seed = *seed;

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; ++i) {
    uint64_t particle_seed = init_seed(base_seed + i, STREAM_SOURCE);
    sites_array[i] = sample_external_source(&particle_seed);
  }
  return 0;
}

} // namespace openmc
