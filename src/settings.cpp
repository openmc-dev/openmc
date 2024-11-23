#include "openmc/settings.h"

#include <cmath>  // for ceil, pow
#include <limits> // for numeric_limits
#include <string>

#include <fmt/core.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/distribution.h"
#include "openmc/distribution_multi.h"
#include "openmc/distribution_spatial.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/mcpl_interface.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/output.h"
#include "openmc/plot.h"
#include "openmc/random_lcg.h"
#include "openmc/random_ray/random_ray.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/string_utils.h"
#include "openmc/tallies/trigger.h"
#include "openmc/volume_calc.h"
#include "openmc/weight_windows.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace settings {

// Default values for boolean flags
bool assume_separate {false};
bool check_overlaps {false};
bool cmfd_run {false};
bool confidence_intervals {false};
bool create_delayed_neutrons {true};
bool create_fission_neutrons {true};
bool delayed_photon_scaling {true};
bool entropy_on {false};
bool event_based {false};
bool legendre_to_tabular {true};
bool material_cell_offsets {true};
bool output_summary {true};
bool output_tallies {true};
bool particle_restart_run {false};
bool photon_transport {false};
bool reduce_tallies {true};
bool res_scat_on {false};
bool restart_run {false};
bool run_CE {true};
bool source_latest {false};
bool source_separate {false};
bool source_write {true};
bool source_mcpl_write {false};
bool surf_source_write {false};
bool surf_mcpl_write {false};
bool surf_source_read {false};
bool survival_biasing {false};
bool temperature_multipole {false};
bool trigger_on {false};
bool trigger_predict {false};
bool uniform_source_sampling {false};
bool ufs_on {false};
bool urr_ptables_on {true};
bool weight_windows_on {false};
bool weight_window_checkpoint_surface {false};
bool weight_window_checkpoint_collision {true};
bool write_all_tracks {false};
bool write_initial_source {false};

std::string path_cross_sections;
std::string path_input;
std::string path_output;
std::string path_particle_restart;
std::string path_sourcepoint;
std::string path_statepoint;
const char* path_statepoint_c {path_statepoint.c_str()};
std::string weight_windows_file;

int32_t n_inactive {0};
int32_t max_lost_particles {10};
double rel_max_lost_particles {1.0e-6};
int32_t max_write_lost_particles {-1};
int32_t gen_per_batch {1};
int64_t n_particles {-1};

int64_t max_particles_in_flight {100000};
int max_particle_events {1000000};

ElectronTreatment electron_treatment {ElectronTreatment::TTB};
array<double, 4> energy_cutoff {0.0, 1000.0, 0.0, 0.0};
array<double, 4> time_cutoff {INFTY, INFTY, INFTY, INFTY};
int legendre_to_tabular_points {C_NONE};
int max_order {0};
int n_log_bins {8000};
int n_batches;
int n_max_batches;
int max_history_splits {10'000'000};
int max_tracks {1000};
ResScatMethod res_scat_method {ResScatMethod::rvs};
double res_scat_energy_min {0.01};
double res_scat_energy_max {1000.0};
vector<std::string> res_scat_nuclides;
RunMode run_mode {RunMode::UNSET};
SolverType solver_type {SolverType::MONTE_CARLO};
std::unordered_set<int> sourcepoint_batch;
std::unordered_set<int> statepoint_batch;
std::unordered_set<int> source_write_surf_id;
int64_t ssw_max_particles;
int64_t ssw_max_files;
int64_t ssw_cell_id {C_NONE};
SSWCellType ssw_cell_type {SSWCellType::None};
TemperatureMethod temperature_method {TemperatureMethod::NEAREST};
double temperature_tolerance {10.0};
double temperature_default {293.6};
array<double, 2> temperature_range {0.0, 0.0};
int trace_batch;
int trace_gen;
int64_t trace_particle;
vector<array<int, 3>> track_identifiers;
int trigger_batch_interval {1};
int verbosity {7};
double weight_cutoff {0.25};
double weight_survive {1.0};

} // namespace settings

//==============================================================================
// Functions
//==============================================================================

void get_run_parameters(pugi::xml_node node_base)
{
  using namespace settings;
  using namespace pugi;

  // Check number of particles
  if (!check_for_node(node_base, "particles")) {
    fatal_error("Need to specify number of particles.");
  }

  // Get number of particles if it wasn't specified as a command-line argument
  if (n_particles == -1) {
    n_particles = std::stoll(get_node_value(node_base, "particles"));
  }

  // Get maximum number of in flight particles for event-based mode
  if (check_for_node(node_base, "max_particles_in_flight")) {
    max_particles_in_flight =
      std::stoll(get_node_value(node_base, "max_particles_in_flight"));
  }

  // Get maximum number of events allowed per particle
  if (check_for_node(node_base, "max_particle_events")) {
    max_particle_events =
      std::stoll(get_node_value(node_base, "max_particle_events"));
  }

  // Get number of basic batches
  if (check_for_node(node_base, "batches")) {
    n_batches = std::stoi(get_node_value(node_base, "batches"));
  }
  if (!trigger_on)
    n_max_batches = n_batches;

  // Get max number of lost particles
  if (check_for_node(node_base, "max_lost_particles")) {
    max_lost_particles =
      std::stoi(get_node_value(node_base, "max_lost_particles"));
  }

  // Get relative number of lost particles
  if (check_for_node(node_base, "rel_max_lost_particles")) {
    rel_max_lost_particles =
      std::stod(get_node_value(node_base, "rel_max_lost_particles"));
  }

  // Get relative number of lost particles
  if (check_for_node(node_base, "max_write_lost_particles")) {
    max_write_lost_particles =
      std::stoi(get_node_value(node_base, "max_write_lost_particles"));
  }

  // Get number of inactive batches
  if (run_mode == RunMode::EIGENVALUE ||
      solver_type == SolverType::RANDOM_RAY) {
    if (check_for_node(node_base, "inactive")) {
      n_inactive = std::stoi(get_node_value(node_base, "inactive"));
    }
    if (check_for_node(node_base, "generations_per_batch")) {
      gen_per_batch =
        std::stoi(get_node_value(node_base, "generations_per_batch"));
    }

    // Preallocate space for keff and entropy by generation
    int m = settings::n_max_batches * settings::gen_per_batch;
    simulation::k_generation.reserve(m);
    simulation::entropy.reserve(m);

    // Get the trigger information for keff
    if (check_for_node(node_base, "keff_trigger")) {
      xml_node node_keff_trigger = node_base.child("keff_trigger");

      if (check_for_node(node_keff_trigger, "type")) {
        auto temp = get_node_value(node_keff_trigger, "type", true, true);
        if (temp == "std_dev") {
          keff_trigger.metric = TriggerMetric::standard_deviation;
        } else if (temp == "variance") {
          keff_trigger.metric = TriggerMetric::variance;
        } else if (temp == "rel_err") {
          keff_trigger.metric = TriggerMetric::relative_error;
        } else {
          fatal_error("Unrecognized keff trigger type " + temp);
        }
      } else {
        fatal_error("Specify keff trigger type in settings XML");
      }

      if (check_for_node(node_keff_trigger, "threshold")) {
        keff_trigger.threshold =
          std::stod(get_node_value(node_keff_trigger, "threshold"));
        if (keff_trigger.threshold <= 0) {
          fatal_error("keff trigger threshold must be positive");
        }
      } else {
        fatal_error("Specify keff trigger threshold in settings XML");
      }
    }
  }

  // Random ray variables
  if (solver_type == SolverType::RANDOM_RAY) {
    xml_node random_ray_node = node_base.child("random_ray");
    if (check_for_node(random_ray_node, "distance_active")) {
      RandomRay::distance_active_ =
        std::stod(get_node_value(random_ray_node, "distance_active"));
      if (RandomRay::distance_active_ <= 0.0) {
        fatal_error("Random ray active distance must be greater than 0");
      }
    } else {
      fatal_error("Specify random ray active distance in settings XML");
    }
    if (check_for_node(random_ray_node, "distance_inactive")) {
      RandomRay::distance_inactive_ =
        std::stod(get_node_value(random_ray_node, "distance_inactive"));
      if (RandomRay::distance_inactive_ < 0) {
        fatal_error(
          "Random ray inactive distance must be greater than or equal to 0");
      }
    } else {
      fatal_error("Specify random ray inactive distance in settings XML");
    }
    if (check_for_node(random_ray_node, "source")) {
      xml_node source_node = random_ray_node.child("source");
      // Get point to list of <source> elements and make sure there is at least
      // one
      RandomRay::ray_source_ = Source::create(source_node);
    } else {
      fatal_error("Specify random ray source in settings XML");
    }
    if (check_for_node(random_ray_node, "volume_estimator")) {
      std::string temp_str =
        get_node_value(random_ray_node, "volume_estimator", true, true);
      if (temp_str == "simulation_averaged") {
        FlatSourceDomain::volume_estimator_ =
          RandomRayVolumeEstimator::SIMULATION_AVERAGED;
      } else if (temp_str == "naive") {
        FlatSourceDomain::volume_estimator_ = RandomRayVolumeEstimator::NAIVE;
      } else if (temp_str == "hybrid") {
        FlatSourceDomain::volume_estimator_ = RandomRayVolumeEstimator::HYBRID;
      } else {
        fatal_error("Unrecognized volume estimator: " + temp_str);
      }
    }
    if (check_for_node(random_ray_node, "source_shape")) {
      std::string temp_str =
        get_node_value(random_ray_node, "source_shape", true, true);
      if (temp_str == "flat") {
        RandomRay::source_shape_ = RandomRaySourceShape::FLAT;
      } else if (temp_str == "linear") {
        RandomRay::source_shape_ = RandomRaySourceShape::LINEAR;
      } else if (temp_str == "linear_xy") {
        RandomRay::source_shape_ = RandomRaySourceShape::LINEAR_XY;
      } else {
        fatal_error("Unrecognized source shape: " + temp_str);
      }
    }
    if (check_for_node(random_ray_node, "volume_normalized_flux_tallies")) {
      FlatSourceDomain::volume_normalized_flux_tallies_ =
        get_node_value_bool(random_ray_node, "volume_normalized_flux_tallies");
    }
    if (check_for_node(random_ray_node, "adjoint")) {
      FlatSourceDomain::adjoint_ =
        get_node_value_bool(random_ray_node, "adjoint");
    }
  }
}

void read_settings_xml()
{
  using namespace settings;
  using namespace pugi;
  // Check if settings.xml exists
  std::string filename = settings::path_input + "settings.xml";
  if (!file_exists(filename)) {
    if (run_mode != RunMode::PLOTTING) {
      fatal_error("Could not find any XML input files! In order to run OpenMC, "
                  "you first need a set of input files; at a minimum, this "
                  "includes settings.xml, geometry.xml, and materials.xml or a "
                  "single model XML file. Please consult the user's guide at "
                  "https://docs.openmc.org for further information.");
    } else {
      // The settings.xml file is optional if we just want to make a plot.
      return;
    }
  }

  // Parse settings.xml file
  xml_document doc;
  auto result = doc.load_file(filename.c_str());
  if (!result) {
    fatal_error("Error processing settings.xml file.");
  }

  // Get root element
  xml_node root = doc.document_element();

  // Verbosity
  if (check_for_node(root, "verbosity")) {
    verbosity = std::stoi(get_node_value(root, "verbosity"));
  }

  // To this point, we haven't displayed any output since we didn't know what
  // the verbosity is. Now that we checked for it, show the title if necessary
  if (mpi::master) {
    if (verbosity >= 2)
      title();
  }

  write_message("Reading settings XML file...", 5);

  read_settings_xml(root);
}

void read_settings_xml(pugi::xml_node root)
{
  using namespace settings;
  using namespace pugi;

  // Find if a multi-group or continuous-energy simulation is desired
  if (check_for_node(root, "energy_mode")) {
    std::string temp_str = get_node_value(root, "energy_mode", true, true);
    if (temp_str == "mg" || temp_str == "multi-group") {
      run_CE = false;
    } else if (temp_str == "ce" || temp_str == "continuous-energy") {
      run_CE = true;
    }
  }

  // Check for user meshes and allocate
  read_meshes(root);

  // Look for deprecated cross_sections.xml file in settings.xml
  if (check_for_node(root, "cross_sections")) {
    warning(
      "Setting cross_sections in settings.xml has been deprecated."
      " The cross_sections are now set in materials.xml and the "
      "cross_sections input to materials.xml and the OPENMC_CROSS_SECTIONS"
      " environment variable will take precendent over setting "
      "cross_sections in settings.xml.");
    path_cross_sections = get_node_value(root, "cross_sections");
  }

  if (!run_CE) {
    // Scattering Treatments
    if (check_for_node(root, "max_order")) {
      max_order = std::stoi(get_node_value(root, "max_order"));
    } else {
      // Set to default of largest int - 1, which means to use whatever is
      // contained in library. This is largest int - 1 because for legendre
      // scattering, a value of 1 is added to the order; adding 1 to the largest
      // int gets you the largest negative integer, which is not what we want.
      max_order = std::numeric_limits<int>::max() - 1;
    }
  }

  // Check for a trigger node and get trigger information
  if (check_for_node(root, "trigger")) {
    xml_node node_trigger = root.child("trigger");

    // Check if trigger(s) are to be turned on
    trigger_on = get_node_value_bool(node_trigger, "active");

    if (trigger_on) {
      if (check_for_node(node_trigger, "max_batches")) {
        n_max_batches = std::stoi(get_node_value(node_trigger, "max_batches"));
      } else {
        fatal_error("<max_batches> must be specified with triggers");
      }

      // Get the batch interval to check triggers
      if (!check_for_node(node_trigger, "batch_interval")) {
        trigger_predict = true;
      } else {
        trigger_batch_interval =
          std::stoi(get_node_value(node_trigger, "batch_interval"));
        if (trigger_batch_interval <= 0) {
          fatal_error("Trigger batch interval must be greater than zero");
        }
      }
    }
  }

  // Check run mode if it hasn't been set from the command line
  xml_node node_mode;
  if (run_mode == RunMode::UNSET) {
    if (check_for_node(root, "run_mode")) {
      std::string temp_str = get_node_value(root, "run_mode", true, true);
      if (temp_str == "eigenvalue") {
        run_mode = RunMode::EIGENVALUE;
      } else if (temp_str == "fixed source") {
        run_mode = RunMode::FIXED_SOURCE;
      } else if (temp_str == "plot") {
        run_mode = RunMode::PLOTTING;
      } else if (temp_str == "particle restart") {
        run_mode = RunMode::PARTICLE;
      } else if (temp_str == "volume") {
        run_mode = RunMode::VOLUME;
      } else {
        fatal_error("Unrecognized run mode: " + temp_str);
      }

      // Assume XML specifies <particles>, <batches>, etc. directly
      node_mode = root;
    } else {
      warning("<run_mode> should be specified.");

      // Make sure that either eigenvalue or fixed source was specified
      node_mode = root.child("eigenvalue");
      if (node_mode) {
        run_mode = RunMode::EIGENVALUE;
      } else {
        node_mode = root.child("fixed_source");
        if (node_mode) {
          run_mode = RunMode::FIXED_SOURCE;
        } else {
          fatal_error("<eigenvalue> or <fixed_source> not specified.");
        }
      }
    }
  }

  // Check solver type
  if (check_for_node(root, "random_ray")) {
    solver_type = SolverType::RANDOM_RAY;
    if (run_CE)
      fatal_error("multi-group energy mode must be specified in settings XML "
                  "when using the random ray solver.");
  }

  if (run_mode == RunMode::EIGENVALUE || run_mode == RunMode::FIXED_SOURCE) {
    // Read run parameters
    get_run_parameters(node_mode);

    // Check number of active batches, inactive batches, max lost particles and
    // particles
    if (n_batches <= n_inactive) {
      fatal_error("Number of active batches must be greater than zero.");
    } else if (n_inactive < 0) {
      fatal_error("Number of inactive batches must be non-negative.");
    } else if (n_particles <= 0) {
      fatal_error("Number of particles must be greater than zero.");
    } else if (max_lost_particles <= 0) {
      fatal_error("Number of max lost particles must be greater than zero.");
    } else if (rel_max_lost_particles <= 0.0 || rel_max_lost_particles >= 1.0) {
      fatal_error("Relative max lost particles must be between zero and one.");
    }
  }

  // Copy plotting random number seed if specified
  if (check_for_node(root, "plot_seed")) {
    auto seed = std::stoll(get_node_value(root, "plot_seed"));
    model::plotter_seed = seed;
  }

  // Copy random number seed if specified
  if (check_for_node(root, "seed")) {
    auto seed = std::stoll(get_node_value(root, "seed"));
    openmc_set_seed(seed);
  }

  // Check for electron treatment
  if (check_for_node(root, "electron_treatment")) {
    auto temp_str = get_node_value(root, "electron_treatment", true, true);
    if (temp_str == "led") {
      electron_treatment = ElectronTreatment::LED;
    } else if (temp_str == "ttb") {
      electron_treatment = ElectronTreatment::TTB;
    } else {
      fatal_error("Unrecognized electron treatment: " + temp_str + ".");
    }
  }

  // Check for photon transport
  if (check_for_node(root, "photon_transport")) {
    photon_transport = get_node_value_bool(root, "photon_transport");

    if (!run_CE && photon_transport) {
      fatal_error("Photon transport is not currently supported in "
                  "multigroup mode");
    }
  }

  // Number of bins for logarithmic grid
  if (check_for_node(root, "log_grid_bins")) {
    n_log_bins = std::stoi(get_node_value(root, "log_grid_bins"));
    if (n_log_bins < 1) {
      fatal_error("Number of bins for logarithmic grid must be greater "
                  "than zero.");
    }
  }

  // Number of OpenMP threads
  if (check_for_node(root, "threads")) {
    if (mpi::master)
      warning("The <threads> element has been deprecated. Use "
              "the OMP_NUM_THREADS environment variable to set the number of "
              "threads.");
  }

  // ==========================================================================
  // EXTERNAL SOURCE

  // Get point to list of <source> elements and make sure there is at least one
  for (pugi::xml_node node : root.children("source")) {
    model::external_sources.push_back(Source::create(node));
  }

  // Check if the user has specified to read surface source
  if (check_for_node(root, "surf_source_read")) {
    surf_source_read = true;
    // Get surface source read node
    xml_node node_ssr = root.child("surf_source_read");

    std::string path = "surface_source.h5";
    // Check if the user has specified different file for surface source reading
    if (check_for_node(node_ssr, "path")) {
      path = get_node_value(node_ssr, "path", false, true);
    }
    model::external_sources.push_back(make_unique<FileSource>(path));
  }

  vector<double> source_strengths;
  for (auto& s : model::external_sources) {
    source_strengths.push_back(s->strength());
  }
  model::external_sources_probability.assign(source_strengths);

  // If no source specified, default to isotropic point source at origin with
  // Watt spectrum. No default source is needed in random ray mode.
  if (model::external_sources.empty() &&
      settings::solver_type != SolverType::RANDOM_RAY) {
    double T[] {0.0};
    double p[] {1.0};
    model::external_sources.push_back(make_unique<IndependentSource>(
      UPtrSpace {new SpatialPoint({0.0, 0.0, 0.0})},
      UPtrAngle {new Isotropic()}, UPtrDist {new Watt(0.988e6, 2.249e-6)},
      UPtrDist {new Discrete(T, p, 1)}));
  }

  // Check if we want to write out source
  if (check_for_node(root, "write_initial_source")) {
    write_initial_source = get_node_value_bool(root, "write_initial_source");
  }

  // Survival biasing
  if (check_for_node(root, "survival_biasing")) {
    survival_biasing = get_node_value_bool(root, "survival_biasing");
  }

  // Probability tables
  if (check_for_node(root, "ptables")) {
    urr_ptables_on = get_node_value_bool(root, "ptables");
  }

  // Cutoffs
  if (check_for_node(root, "cutoff")) {
    xml_node node_cutoff = root.child("cutoff");
    if (check_for_node(node_cutoff, "weight")) {
      weight_cutoff = std::stod(get_node_value(node_cutoff, "weight"));
    }
    if (check_for_node(node_cutoff, "weight_avg")) {
      weight_survive = std::stod(get_node_value(node_cutoff, "weight_avg"));
    }
    if (check_for_node(node_cutoff, "energy_neutron")) {
      energy_cutoff[0] =
        std::stod(get_node_value(node_cutoff, "energy_neutron"));
    } else if (check_for_node(node_cutoff, "energy")) {
      warning("The use of an <energy> cutoff is deprecated and should "
              "be replaced by <energy_neutron>.");
      energy_cutoff[0] = std::stod(get_node_value(node_cutoff, "energy"));
    }
    if (check_for_node(node_cutoff, "energy_photon")) {
      energy_cutoff[1] =
        std::stod(get_node_value(node_cutoff, "energy_photon"));
    }
    if (check_for_node(node_cutoff, "energy_electron")) {
      energy_cutoff[2] =
        std::stof(get_node_value(node_cutoff, "energy_electron"));
    }
    if (check_for_node(node_cutoff, "energy_positron")) {
      energy_cutoff[3] =
        std::stod(get_node_value(node_cutoff, "energy_positron"));
    }
    if (check_for_node(node_cutoff, "time_neutron")) {
      time_cutoff[0] = std::stod(get_node_value(node_cutoff, "time_neutron"));
    }
    if (check_for_node(node_cutoff, "time_photon")) {
      time_cutoff[1] = std::stod(get_node_value(node_cutoff, "time_photon"));
    }
    if (check_for_node(node_cutoff, "time_electron")) {
      time_cutoff[2] = std::stod(get_node_value(node_cutoff, "time_electron"));
    }
    if (check_for_node(node_cutoff, "time_positron")) {
      time_cutoff[3] = std::stod(get_node_value(node_cutoff, "time_positron"));
    }
  }

  // Particle trace
  if (check_for_node(root, "trace")) {
    auto temp = get_node_array<int64_t>(root, "trace");
    if (temp.size() != 3) {
      fatal_error("Must provide 3 integers for <trace> that specify the "
                  "batch, generation, and particle number.");
    }
    trace_batch = temp.at(0);
    trace_gen = temp.at(1);
    trace_particle = temp.at(2);
  }

  // Particle tracks
  if (check_for_node(root, "track")) {
    // Get values and make sure there are three per particle
    auto temp = get_node_array<int>(root, "track");
    if (temp.size() % 3 != 0) {
      fatal_error(
        "Number of integers specified in 'track' is not "
        "divisible by 3.  Please provide 3 integers per particle to be "
        "tracked.");
    }

    // Reshape into track_identifiers
    int n_tracks = temp.size() / 3;
    for (int i = 0; i < n_tracks; ++i) {
      track_identifiers.push_back(
        {temp[3 * i], temp[3 * i + 1], temp[3 * i + 2]});
    }
  }

  // Shannon entropy
  if (solver_type == SolverType::RANDOM_RAY) {
    if (check_for_node(root, "entropy_mesh")) {
      fatal_error("Random ray uses FSRs to compute the Shannon entropy. "
                  "No user-defined entropy mesh is supported.");
    }
    entropy_on = true;
  } else if (solver_type == SolverType::MONTE_CARLO) {
    if (check_for_node(root, "entropy_mesh")) {
      int temp = std::stoi(get_node_value(root, "entropy_mesh"));
      if (model::mesh_map.find(temp) == model::mesh_map.end()) {
        fatal_error(fmt::format(
          "Mesh {} specified for Shannon entropy does not exist.", temp));
      }

      auto* m = dynamic_cast<RegularMesh*>(
        model::meshes[model::mesh_map.at(temp)].get());
      if (!m)
        fatal_error("Only regular meshes can be used as an entropy mesh");
      simulation::entropy_mesh = m;

      // Turn on Shannon entropy calculation
      entropy_on = true;

    } else if (check_for_node(root, "entropy")) {
      fatal_error(
        "Specifying a Shannon entropy mesh via the <entropy> element "
        "is deprecated. Please create a mesh using <mesh> and then reference "
        "it by specifying its ID in an <entropy_mesh> element.");
    }
  }
  // Uniform fission source weighting mesh
  if (check_for_node(root, "ufs_mesh")) {
    auto temp = std::stoi(get_node_value(root, "ufs_mesh"));
    if (model::mesh_map.find(temp) == model::mesh_map.end()) {
      fatal_error(fmt::format("Mesh {} specified for uniform fission site "
                              "method does not exist.",
        temp));
    }

    auto* m =
      dynamic_cast<RegularMesh*>(model::meshes[model::mesh_map.at(temp)].get());
    if (!m)
      fatal_error("Only regular meshes can be used as a UFS mesh");
    simulation::ufs_mesh = m;

    // Turn on uniform fission source weighting
    ufs_on = true;

  } else if (check_for_node(root, "uniform_fs")) {
    fatal_error(
      "Specifying a UFS mesh via the <uniform_fs> element "
      "is deprecated. Please create a mesh using <mesh> and then reference "
      "it by specifying its ID in a <ufs_mesh> element.");
  }

  // Check if the user has specified to write state points
  if (check_for_node(root, "state_point")) {

    // Get pointer to state_point node
    auto node_sp = root.child("state_point");

    // Determine number of batches at which to store state points
    if (check_for_node(node_sp, "batches")) {
      // User gave specific batches to write state points
      auto temp = get_node_array<int>(node_sp, "batches");
      for (const auto& b : temp) {
        statepoint_batch.insert(b);
      }
    } else {
      // If neither were specified, write state point at last batch
      statepoint_batch.insert(n_batches);
    }
  } else {
    // If no <state_point> tag was present, by default write state point at
    // last batch only
    statepoint_batch.insert(n_batches);
  }

  // Check if the user has specified to write source points
  if (check_for_node(root, "source_point")) {
    // Get source_point node
    xml_node node_sp = root.child("source_point");

    // Determine batches at which to store source points
    if (check_for_node(node_sp, "batches")) {
      // User gave specific batches to write source points
      auto temp = get_node_array<int>(node_sp, "batches");
      for (const auto& b : temp) {
        sourcepoint_batch.insert(b);
      }
    } else {
      // If neither were specified, write source points with state points
      sourcepoint_batch = statepoint_batch;
    }

    // Check if the user has specified to write binary source file
    if (check_for_node(node_sp, "separate")) {
      source_separate = get_node_value_bool(node_sp, "separate");
    }
    if (check_for_node(node_sp, "write")) {
      source_write = get_node_value_bool(node_sp, "write");
    }
    if (check_for_node(node_sp, "mcpl")) {
      source_mcpl_write = get_node_value_bool(node_sp, "mcpl");

      // Make sure MCPL support is enabled
      if (source_mcpl_write && !MCPL_ENABLED) {
        fatal_error(
          "Your build of OpenMC does not support writing MCPL source files.");
      }
    }
    if (check_for_node(node_sp, "overwrite_latest")) {
      source_latest = get_node_value_bool(node_sp, "overwrite_latest");
      source_separate = source_latest;
    }
  } else {
    // If no <source_point> tag was present, by default we keep source bank in
    // statepoint file and write it out at statepoints intervals
    source_separate = false;
    sourcepoint_batch = statepoint_batch;
  }

  // Check is the user specified to convert strength to statistical weight
  if (check_for_node(root, "uniform_source_sampling")) {
    uniform_source_sampling =
      get_node_value_bool(root, "uniform_source_sampling");
  }

  // Check if the user has specified to write surface source
  if (check_for_node(root, "surf_source_write")) {
    surf_source_write = true;
    // Get surface source write node
    xml_node node_ssw = root.child("surf_source_write");

    // Determine surface ids at which crossing particles are to be banked.
    // If no surfaces are specified, all surfaces in the model will be used
    // to bank source points.
    if (check_for_node(node_ssw, "surface_ids")) {
      auto temp = get_node_array<int>(node_ssw, "surface_ids");
      for (const auto& b : temp) {
        source_write_surf_id.insert(b);
      }
    }

    // Get maximum number of particles to be banked per surface
    if (check_for_node(node_ssw, "max_particles")) {
      ssw_max_particles = std::stoll(get_node_value(node_ssw, "max_particles"));
    } else {
      fatal_error("A maximum number of particles needs to be specified "
                  "using the 'max_particles' parameter to store surface "
                  "source points.");
    }

    // Get maximum number of surface source files to be created
    if (check_for_node(node_ssw, "max_source_files")) {
      ssw_max_files = std::stoll(get_node_value(node_ssw, "max_source_files"));
    } else {
      ssw_max_files = 1;
    }

    if (check_for_node(node_ssw, "mcpl")) {
      surf_mcpl_write = get_node_value_bool(node_ssw, "mcpl");

      // Make sure MCPL support is enabled
      if (surf_mcpl_write && !MCPL_ENABLED) {
        fatal_error("Your build of OpenMC does not support writing MCPL "
                    "surface source files.");
      }
    }
    // Get cell information
    if (check_for_node(node_ssw, "cell")) {
      ssw_cell_id = std::stoll(get_node_value(node_ssw, "cell"));
      ssw_cell_type = SSWCellType::Both;
    }
    if (check_for_node(node_ssw, "cellfrom")) {
      if (ssw_cell_id != C_NONE) {
        fatal_error(
          "'cell', 'cellfrom' and 'cellto' cannot be used at the same time.");
      }
      ssw_cell_id = std::stoll(get_node_value(node_ssw, "cellfrom"));
      ssw_cell_type = SSWCellType::From;
    }
    if (check_for_node(node_ssw, "cellto")) {
      if (ssw_cell_id != C_NONE) {
        fatal_error(
          "'cell', 'cellfrom' and 'cellto' cannot be used at the same time.");
      }
      ssw_cell_id = std::stoll(get_node_value(node_ssw, "cellto"));
      ssw_cell_type = SSWCellType::To;
    }
  }

  // If source is not separate and is to be written out in the statepoint file,
  // make sure that the sourcepoint batch numbers are contained in the
  // statepoint list
  if (!source_separate) {
    for (const auto& b : sourcepoint_batch) {
      if (!contains(statepoint_batch, b)) {
        fatal_error(
          "Sourcepoint batches are not a subset of statepoint batches.");
      }
    }
  }

  // Check if the user has specified to not reduce tallies at the end of every
  // batch
  if (check_for_node(root, "no_reduce")) {
    reduce_tallies = !get_node_value_bool(root, "no_reduce");
  }

  // Check if the user has specified to use confidence intervals for
  // uncertainties rather than standard deviations
  if (check_for_node(root, "confidence_intervals")) {
    confidence_intervals = get_node_value_bool(root, "confidence_intervals");
  }

  // Check for output options
  if (check_for_node(root, "output")) {
    // Get pointer to output node
    pugi::xml_node node_output = root.child("output");

    // Check for summary option
    if (check_for_node(node_output, "summary")) {
      output_summary = get_node_value_bool(node_output, "summary");
    }

    // Check for ASCII tallies output option
    if (check_for_node(node_output, "tallies")) {
      output_tallies = get_node_value_bool(node_output, "tallies");
    }

    // Set output directory if a path has been specified
    if (check_for_node(node_output, "path")) {
      path_output = get_node_value(node_output, "path");
      if (!ends_with(path_output, "/")) {
        path_output += "/";
      }
    }
  }

  // Resonance scattering parameters
  if (check_for_node(root, "resonance_scattering")) {
    xml_node node_res_scat = root.child("resonance_scattering");

    // See if resonance scattering is enabled
    if (check_for_node(node_res_scat, "enable")) {
      res_scat_on = get_node_value_bool(node_res_scat, "enable");
    } else {
      res_scat_on = true;
    }

    // Determine what method is used
    if (check_for_node(node_res_scat, "method")) {
      auto temp = get_node_value(node_res_scat, "method", true, true);
      if (temp == "rvs") {
        res_scat_method = ResScatMethod::rvs;
      } else if (temp == "dbrc") {
        res_scat_method = ResScatMethod::dbrc;
      } else {
        fatal_error(
          "Unrecognized resonance elastic scattering method: " + temp + ".");
      }
    }

    // Minimum energy for resonance scattering
    if (check_for_node(node_res_scat, "energy_min")) {
      res_scat_energy_min =
        std::stod(get_node_value(node_res_scat, "energy_min"));
    }
    if (res_scat_energy_min < 0.0) {
      fatal_error("Lower resonance scattering energy bound is negative");
    }

    // Maximum energy for resonance scattering
    if (check_for_node(node_res_scat, "energy_max")) {
      res_scat_energy_max =
        std::stod(get_node_value(node_res_scat, "energy_max"));
    }
    if (res_scat_energy_max < res_scat_energy_min) {
      fatal_error("Upper resonance scattering energy bound is below the "
                  "lower resonance scattering energy bound.");
    }

    // Get resonance scattering nuclides
    if (check_for_node(node_res_scat, "nuclides")) {
      res_scat_nuclides =
        get_node_array<std::string>(node_res_scat, "nuclides");
    }
  }

  // Get volume calculations
  for (pugi::xml_node node_vol : root.children("volume_calc")) {
    model::volume_calcs.emplace_back(node_vol);
  }

  // Get temperature settings
  if (check_for_node(root, "temperature_default")) {
    temperature_default =
      std::stod(get_node_value(root, "temperature_default"));
  }
  if (check_for_node(root, "temperature_method")) {
    auto temp = get_node_value(root, "temperature_method", true, true);
    if (temp == "nearest") {
      temperature_method = TemperatureMethod::NEAREST;
    } else if (temp == "interpolation") {
      temperature_method = TemperatureMethod::INTERPOLATION;
    } else {
      fatal_error("Unknown temperature method: " + temp);
    }
  }
  if (check_for_node(root, "temperature_tolerance")) {
    temperature_tolerance =
      std::stod(get_node_value(root, "temperature_tolerance"));
  }
  if (check_for_node(root, "temperature_multipole")) {
    temperature_multipole = get_node_value_bool(root, "temperature_multipole");

    // Multipole currently doesn't work with photon transport
    if (temperature_multipole && photon_transport) {
      fatal_error("Multipole data cannot currently be used in conjunction with "
                  "photon transport.");
    }
  }
  if (check_for_node(root, "temperature_range")) {
    auto range = get_node_array<double>(root, "temperature_range");
    temperature_range[0] = range.at(0);
    temperature_range[1] = range.at(1);
  }

  // Check for tabular_legendre options
  if (check_for_node(root, "tabular_legendre")) {
    // Get pointer to tabular_legendre node
    xml_node node_tab_leg = root.child("tabular_legendre");

    // Check for enable option
    if (check_for_node(node_tab_leg, "enable")) {
      legendre_to_tabular = get_node_value_bool(node_tab_leg, "enable");
    }

    // Check for the number of points
    if (check_for_node(node_tab_leg, "num_points")) {
      legendre_to_tabular_points =
        std::stoi(get_node_value(node_tab_leg, "num_points"));
      if (legendre_to_tabular_points <= 1 && !run_CE) {
        fatal_error(
          "The 'num_points' subelement/attribute of the "
          "<tabular_legendre> element must contain a value greater than 1");
      }
    }
  }

  // Check whether create delayed neutrons in fission
  if (check_for_node(root, "create_delayed_neutrons")) {
    create_delayed_neutrons =
      get_node_value_bool(root, "create_delayed_neutrons");
  }

  // Check whether create fission sites
  if (run_mode == RunMode::FIXED_SOURCE) {
    if (check_for_node(root, "create_fission_neutrons")) {
      create_fission_neutrons =
        get_node_value_bool(root, "create_fission_neutrons");
    }
  }

  // Check whether to scale fission photon yields
  if (check_for_node(root, "delayed_photon_scaling")) {
    delayed_photon_scaling =
      get_node_value_bool(root, "delayed_photon_scaling");
  }

  // Check whether to use event-based parallelism
  if (check_for_node(root, "event_based")) {
    event_based = get_node_value_bool(root, "event_based");
  }

  // Check whether material cell offsets should be generated
  if (check_for_node(root, "material_cell_offsets")) {
    material_cell_offsets = get_node_value_bool(root, "material_cell_offsets");
  }

  // Weight window information
  for (pugi::xml_node node_ww : root.children("weight_windows")) {
    variance_reduction::weight_windows.emplace_back(
      std::make_unique<WeightWindows>(node_ww));
  }

  // Enable weight windows by default if one or more are present
  if (variance_reduction::weight_windows.size() > 0)
    settings::weight_windows_on = true;

  // read weight windows from file
  if (check_for_node(root, "weight_windows_file")) {
    weight_windows_file = get_node_value(root, "weight_windows_file");
  }

  // read settings for weight windows value, this will override
  // the automatic setting even if weight windows are present
  if (check_for_node(root, "weight_windows_on")) {
    weight_windows_on = get_node_value_bool(root, "weight_windows_on");
  }

  if (check_for_node(root, "max_history_splits")) {
    settings::max_history_splits =
      std::stoi(get_node_value(root, "max_history_splits"));
  }

  if (check_for_node(root, "max_tracks")) {
    settings::max_tracks = std::stoi(get_node_value(root, "max_tracks"));
  }

  // Create weight window generator objects
  if (check_for_node(root, "weight_window_generators")) {
    auto wwgs_node = root.child("weight_window_generators");
    for (pugi::xml_node node_wwg :
      wwgs_node.children("weight_windows_generator")) {
      variance_reduction::weight_windows_generators.emplace_back(
        std::make_unique<WeightWindowsGenerator>(node_wwg));
    }
    // if any of the weight windows are intended to be generated otf, make sure
    // they're applied
    for (const auto& wwg : variance_reduction::weight_windows_generators) {
      if (wwg->on_the_fly_) {
        settings::weight_windows_on = true;
        break;
      }
    }
  }

  // Set up weight window checkpoints
  if (check_for_node(root, "weight_window_checkpoints")) {
    xml_node ww_checkpoints = root.child("weight_window_checkpoints");
    if (check_for_node(ww_checkpoints, "collision")) {
      weight_window_checkpoint_collision =
        get_node_value_bool(ww_checkpoints, "collision");
    }
    if (check_for_node(ww_checkpoints, "surface")) {
      weight_window_checkpoint_surface =
        get_node_value_bool(ww_checkpoints, "surface");
    }
  }
}

void free_memory_settings()
{
  settings::statepoint_batch.clear();
  settings::sourcepoint_batch.clear();
  settings::source_write_surf_id.clear();
  settings::res_scat_nuclides.clear();
}

//==============================================================================
// C API functions
//==============================================================================

extern "C" int openmc_set_n_batches(
  int32_t n_batches, bool set_max_batches, bool add_statepoint_batch)
{
  if (settings::n_inactive >= n_batches) {
    set_errmsg("Number of active batches must be greater than zero.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  if (simulation::current_batch >= n_batches) {
    set_errmsg("Number of batches must be greater than current batch.");
    return OPENMC_E_INVALID_ARGUMENT;
  }

  if (!settings::trigger_on) {
    // Set n_batches and n_max_batches to same value
    settings::n_batches = n_batches;
    settings::n_max_batches = n_batches;
  } else {
    // Set n_batches and n_max_batches based on value of set_max_batches
    if (set_max_batches) {
      settings::n_max_batches = n_batches;
    } else {
      settings::n_batches = n_batches;
    }
  }

  // Update size of k_generation and entropy
  int m = settings::n_max_batches * settings::gen_per_batch;
  simulation::k_generation.reserve(m);
  simulation::entropy.reserve(m);

  // Add value of n_batches to statepoint_batch
  if (add_statepoint_batch &&
      !(contains(settings::statepoint_batch, n_batches)))
    settings::statepoint_batch.insert(n_batches);

  return 0;
}

extern "C" int openmc_get_n_batches(int* n_batches, bool get_max_batches)
{
  *n_batches = get_max_batches ? settings::n_max_batches : settings::n_batches;

  return 0;
}

} // namespace openmc
