#include "openmc/volume_calc.h"

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/distribution_multi.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/nuclide.h"
#include "openmc/openmp_interface.h"
#include "openmc/output.h"
#include "openmc/plot.h"
#include "openmc/random_dist.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/timer.h"
#include "openmc/xml_interface.h"

#include "xtensor/xadapt.hpp"
#include "xtensor/xview.hpp"
#include <fmt/core.h>

#include <algorithm> // for copy
#include <cmath>     // for pow, sqrt
#include <unordered_set>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
vector<VolumeCalculation> volume_calcs;
} // namespace model

#ifdef OPENMC_MPI
static MPI_Datatype mpi_vol_results;  //!< MPI struct for CalcResults
static MPI_Datatype mpi_volume_tally; //!< MPI struct for VolTally
#endif

//==============================================================================
// VolumeCalculation implementation
//==============================================================================

VolumeCalculation::VolumeCalculation(pugi::xml_node node)
{
  // Read domain type (cell, material or universe)
  std::string domain_type = get_node_value(node, "domain_type");
  if (domain_type == "cell") {
    domain_type_ = TallyDomain::CELL;
  } else if (domain_type == "material") {
    domain_type_ = TallyDomain::MATERIAL;
  } else if (domain_type == "universe") {
    domain_type_ = TallyDomain::UNIVERSE;
  } else {
    fatal_error(std::string("Unrecognized domain type for stochastic "
                            "volume calculation: " +
                            domain_type));
  }

  // Read domain IDs, bounding corodinates and number of samples
  domain_ids_ = get_node_array<int>(node, "domain_ids");
  lower_left_ = get_node_array<double>(node, "lower_left");
  upper_right_ = get_node_array<double>(node, "upper_right");
  n_samples_ = std::stoull(get_node_value(node, "samples"));

  // Determine volume of bounding box
  const Position d {upper_right_ - lower_left_};
  volume_sample_ = d.x * d.y * d.z;

  if (check_for_node(node, "threshold")) {
    pugi::xml_node threshold_node = node.child("threshold");

    threshold_ = std::stod(get_node_value(threshold_node, "threshold"));
    if (threshold_ <= 0.0) {
      fatal_error(fmt::format("Invalid error threshold {} provided for a "
                              "volume calculation.",
        threshold_));
    }

    // Prior calculation of the trigger volume fractions-related values t'
    // deriviated from the given volumes-related values t to prevent excessive
    // repeated computations during Monte Carlo execution. The values of t' are
    // computed via the sample mean \bar{x} and the adjusted sample variance
    // s^2.
    std::string tmp = get_node_value(threshold_node, "type");
    if (tmp == "variance") {
      trigger_type_ = TriggerMetric::variance;
      // Condition: s^2 < t' / N
      // t' = s^2 = t / (V_b)^2, in sq. volume fraction units
      threshold_cnd_ = threshold_ / std::pow(volume_sample_, 2);
      ;
    } else if (tmp == "std_dev") {
      trigger_type_ = TriggerMetric::standard_deviation;
      // Condition: s^2 < t' / N
      // t' = s^2 = (t / V_b)^2, in sq. volume fraction units
      threshold_cnd_ = std::pow(threshold_ / volume_sample_, 2);
    } else if (tmp == "rel_err") {
      trigger_type_ = TriggerMetric::relative_error;
      // Condition: s^2 / \bar{x}^2 < t' / N
      // t' = s^2 / \bar{x}^2 = t^2, in relative units
      threshold_cnd_ = threshold_ * threshold_;
    } else {
      fatal_error(fmt::format(
        "Invalid volume calculation trigger type '{}' provided.", tmp));
    }

    if (check_for_node(threshold_node, "max_iterations")) {
      max_iterations_ =
        std::stoi(get_node_value(threshold_node, "max_iterations"));
      if (max_iterations_ <= 0) {
        fatal_error(fmt::format(
          "Invalid error max_iterations {} provided.", max_iterations_));
      }
    }
  }

  // Ensure there are no duplicates by copying elements to a set and then
  // comparing the length with the original vector
  std::unordered_set<int> unique_ids(domain_ids_.cbegin(), domain_ids_.cend());
  if (unique_ids.size() != domain_ids_.size()) {
    throw std::runtime_error {"Domain IDs for a volume calculation "
                              "must be unique."};
  }

  // Type of volume estimator
  std::string tmp = get_node_value(node, "estimator_type");
  if (tmp == "hit") {
    mode_ = EstMode::REJECTION;
  } else if (tmp == "ray") {
    mode_ = EstMode::RAYTRACE;
  } else {
    fatal_error(
      fmt::format("Invalid volume calculation mode '{}' provided.", tmp));
  }
}

void VolumeCalculation::execute(CalcResults& master_results) const
{
#ifdef OPENMC_MPI
  // MPI types are commited in the beginning of calculation and freed on the
  // return, remaining a MPI type after the return produces a MPICH warning for
  // memory leak
  this->initialize_MPI_struct();
#endif

  // Check to make sure domain IDs are valid
  for (auto uid : domain_ids_) {
    switch (domain_type_) {
    case TallyDomain::CELL:
      if (model::cell_map.find(uid) == model::cell_map.end()) {
        throw std::runtime_error {fmt::format(
          "Cell {} in volume calculation does not exist in geometry.", uid)};
      }
      break;
    case TallyDomain::MATERIAL:
      if (model::material_map.find(uid) == model::material_map.end()) {
        throw std::runtime_error {fmt::format(
          "Material {} in volume calculation does not exist in geometry.",
          uid)};
      }
      break;
    case TallyDomain::UNIVERSE:
      if (model::universe_map.find(uid) == model::universe_map.end()) {
        throw std::runtime_error {fmt::format(
          "Universe {} in volume calculation does not exist in geometry.",
          uid)};
      }
    }
  }

  // Shared data that is collected from all threads
  const int n_domains = domain_ids_.size();

  // Divide work over MPI processes
  uint64_t min_samples = n_samples_ / mpi::n_procs;
  uint64_t remainder = n_samples_ % mpi::n_procs;
  uint64_t i_start, i_end;
  if (mpi::rank < remainder) {
    i_start = (min_samples + 1) * mpi::rank;
    i_end = i_start + min_samples + 1;
  } else {
    i_start =
      (min_samples + 1) * remainder + (mpi::rank - remainder) * min_samples;
    i_end = i_start + min_samples;
  }

  while (true) {

#pragma omp parallel
    {
      // Temporary variables that are private to each thread
      CalcResults results(*this);
      results.sampling_time.start();

      // Sample locations and count scores
#pragma omp for
      for (size_t i = i_start; i < i_end; i++) {
        uint64_t id = master_results.iterations * n_samples_ + i;
        uint64_t seed = init_seed(id, STREAM_VOLUME);

        Position r {uniform_distribution(lower_left_.x, upper_right_.x, &seed),
          uniform_distribution(lower_left_.y, upper_right_.y, &seed),
          uniform_distribution(lower_left_.z, upper_right_.z, &seed)};

        results.n_samples++;

        switch (mode_) {
        case EstMode::REJECTION: {
          constexpr double flt3 = 1. / std::sqrt(3.);
          Direction u {flt3, flt3, flt3};

          // Create zero-length ray, it is a bit excessive due to undemanded
          // internal ray variables initialization
          VolEstRay ray_hit(r, u, 0., 1., *this, results);
          ray_hit.score_hit();
        } break;
        case EstMode::RAYTRACE: {
          Direction u = isotropic_direction(&seed);

          // Compute lengths of the bounding box's chord segments
          const auto chord_len = get_box_chord(r, u);
          const double ch_len_tot = -chord_len.first + chord_len.second;

          r += chord_len.first * u;
          VolEstRay ray(r, u, ch_len_tot, 1. / ch_len_tot, *this, results);
          ray.trace(); // Trace from a boundary to another
        }
        }

        // This passing across all tallies after each sample can be
        // inefficient for the case of large number of domains and small
        // size of batch, but the batch size is currently assigned to 1
        // for keep the input format being unchanged
        const double batch_size_1 = 1. / static_cast<double>(1);
        for (auto& vt : results.vol_tallies) {
          for (auto& vol_tally : vt) {
            vol_tally.finalize_batch(batch_size_1);
          }
        }
      } // sample/batch loop

      results.sampling_time.stop();
      results.cost = results.sampling_time.elapsed();

      // At this point, each thread has its own volume tallies lists and we
      // now need to reduce them. OpenMP is not nearly smart enough to do this
      // on its own, so we have to manually reduce them
      reduce_results(results, master_results);
    } // omp parallel

    // bump iteration counter
    master_results.iterations++;

#ifdef OPENMC_MPI
    master_results.collect_MPI(); // collect results to master process
#endif
    // Process volume estimation results in master process for the trigger state
    // determination
    bool stop_calc =
      mpi::master && (trigger_type_ == TriggerMetric::not_active ||
                       master_results.iterations == max_iterations_);

    if (!stop_calc) {
      // Compute current trigger state across totals (0th elements) only
      for (const auto& vt : master_results.vol_tallies) {
        stop_calc = vt[0].trigger_state(
          trigger_type_, threshold_cnd_, master_results.n_samples);
        if (!stop_calc)
          break;
      }
    }

#ifdef OPENMC_MPI
    // Send the state of calculation continuation just obtained in master
    // process to all processes
    MPI_Bcast(&stop_calc, 1, MPI_CXX_BOOL, 0, mpi::intracomm);
#endif

    if (!stop_calc)
      continue; // while loop
    // No trigger is applied or the trigger condition is satisfied, we're
    // done

    if (mpi::master) {
      // Normalize all results on the bounding primitive volume and compute
      // stddev
      for (auto& vt : master_results.vol_tallies) {
        for (auto& vol_tally : vt) {
          vol_tally.score_acc = get_tally_results(
            master_results.n_samples, volume_sample_, vol_tally);
        }
      }

      // Compute nuclides
      for (int i_domain = 0; i_domain < n_domains; ++i_domain) {
        // Create 2D array to store atoms/uncertainty for each nuclide. Later,
        // this is compressed into vectors storing only those nuclides that are
        // non-zero
        auto n_nuc =
          settings::run_CE ? data::nuclides.size() : data::mg.nuclides_.size();
        xt::xtensor<double, 2> atoms({n_nuc, 2}, 0.0);

        for (int j = 0; j < master_results.vol_tallies[i_domain].size(); ++j) {
          const int i_material = master_results.vol_tallies[i_domain][j].index;
          if (i_material == MATERIAL_VOID || i_material == _INDEX_TOTAL)
            continue;

          const auto& mat = model::materials[i_material];
          for (int k = 0; k < mat->nuclide_.size(); ++k) {
            auto& volume = master_results.vol_tallies[i_domain][j].score_acc;
            // Collect calculated nuclide amounts N [atoms] and stddev as
            // N = V [cm^3] * \ro [atoms/b-cm] * 1.e24 [b-cm/cm^3]
            atoms(mat->nuclide_[k], 0) +=
              volume[0] * mat->atom_density(k) * 1.0e24;
            atoms(mat->nuclide_[k], 1) +=
              volume[1] * mat->atom_density(k) * 1.0e24;
          }
        }

        // Get reference to result for this domain
        auto& result {master_results.nuc_results[i_domain]};
        // Convert full arrays to vectors
        for (int j = 0; j < n_nuc; ++j) {
          if (atoms(j, 0) > 0.0) {
            result.nuclides.push_back(j);
            result.atoms.push_back(atoms(j, 0));
            result.uncertainty.push_back(atoms(j, 1));
          }
        }
      } // end domains loop
    }

#ifdef OPENMC_MPI
    // MPI types are commited in the beginning of calculation and freed on
    // return, lefting MPI types after return produces MPICH warnings about
    // memory leak
    this->delete_MPI_struct();
#endif

    return;

  } // end while
}

void VolumeCalculation::show_vol_stat(
  const std::string label, const std::string units, const double value) const
{
  fmt::print("{0:<20} = {2:10.4e} {1:<}\n", label, units, value);
}

void VolumeCalculation::show_volume(const std::string domain_type,
  const int domain_id, const std::string region_name, const double mean,
  const double stddev) const
{
  fmt::print(" {0:>9}{1:>6}: {2:10.4e} +/- {3:10.4e} cm^3", domain_type,
    domain_id, mean, stddev);
  if (!region_name.empty()) {
    fmt::print("   //{:<}", region_name);
  }
  fmt::print("\n");
}

void VolumeCalculation::show_results(const CalcResults& results) const
{
  // Show tracing statistics
  write_message(5, " ");
  show_vol_stat(
    "Total sample size", "", static_cast<double>(results.n_samples));
  show_vol_stat("Running cost", "thread-sec", results.cost);

  switch (mode_) {
  case EstMode::REJECTION:
    show_vol_stat("Cost of hitting", "thread-sec/hit",
      static_cast<double>(results.cost) /
        static_cast<double>(results.n_samples));
    break;
  case EstMode::RAYTRACE:
    show_vol_stat("Rays traced", "", results.n_rays);
    show_vol_stat("Average ray length", "segments/ray",
      static_cast<double>(results.n_segs) /
        static_cast<double>(results.n_rays));
    show_vol_stat("Cost of tracing", "thread-sec/segment",
      static_cast<double>(results.cost) / static_cast<double>(results.n_segs));
    if (results.n_errors != 0)
      show_vol_stat("Error rate", "errors/ray",
        static_cast<double>(results.n_errors) /
          static_cast<double>(results.n_rays));
  }
  write_message(5, " ");

  std::string domain_type;
  if (domain_type_ == TallyDomain::CELL) {
    domain_type = "Cell";
  } else if (domain_type_ == TallyDomain::MATERIAL) {
    domain_type = "Material";
  } else {
    domain_type = "Universe";
  }

  // Display domain volumes
  for (int j = 0; j < domain_ids_.size(); j++) {
    std::string region_name {""};
    if (domain_type_ == TallyDomain::CELL) {
      int cell_idx = model::cell_map[domain_ids_[j]];
      region_name = model::cells[cell_idx]->name();
    } else if (domain_type_ == TallyDomain::MATERIAL) {
      int mat_idx = model::material_map[domain_ids_[j]];
      region_name = model::materials[mat_idx]->name();
    }
    if (region_name.size())
      region_name.insert(0, " "); // prepend space for formatting

    show_volume(domain_type, domain_ids_[j], region_name,
      results.vol_tallies[j][0].score_acc[0],
      results.vol_tallies[j][0].score_acc[1]);
  }
  write_message(4, " "); // Blank line afer results printed
}

void VolumeCalculation::to_hdf5(
  const std::string& filename, const CalcResults& results) const
{
  // Create HDF5 file
  hid_t file_id = file_open(filename, 'w');

  // Write header info
  write_attribute(file_id, "filetype", "volume");
  write_attribute(file_id, "version", VERSION_VOLUME);
  write_attribute(file_id, "openmc_version", VERSION);
#ifdef GIT_SHA1
  write_attribute(file_id, "git_sha1", GIT_SHA1);
#endif

  // Write current date and time
  write_attribute(file_id, "date_and_time", time_stamp());

  // Write basic metadata
  write_attribute(file_id, "samples", n_samples_);
  write_attribute(file_id, "lower_left", lower_left_);
  write_attribute(file_id, "upper_right", upper_right_);
  // Write estimator info
  std::string estimator_str;
  switch (mode_) {
  case EstMode::REJECTION:
    estimator_str = "hit";
    break;
  case EstMode::RAYTRACE:
    estimator_str = "ray";
    break;
  }
  write_attribute(file_id, "estimator_type", estimator_str);
  // Write trigger info
  if (trigger_type_ != TriggerMetric::not_active) {
    write_attribute(file_id, "iterations", results.iterations);
    write_attribute(file_id, "threshold", threshold_);
    std::string trigger_str;
    switch (trigger_type_) {
    case TriggerMetric::variance:
      trigger_str = "variance";
      break;
    case TriggerMetric::standard_deviation:
      trigger_str = "std_dev";
      break;
    case TriggerMetric::relative_error:
      trigger_str = "rel_err";
      break;
    default:
      break;
    }
    write_attribute(file_id, "trigger_type", trigger_str);
    // check max_iterations on default value
    if (max_iterations_ !=
        std::numeric_limits<decltype(max_iterations_)>::max())
      write_attribute(file_id, "max_iterations", max_iterations_);
  } else {
    write_attribute(file_id, "iterations", 1);
  }

  if (domain_type_ == TallyDomain::CELL) {
    write_attribute(file_id, "domain_type", "cell");
  } else if (domain_type_ == TallyDomain::MATERIAL) {
    write_attribute(file_id, "domain_type", "material");
  } else if (domain_type_ == TallyDomain::UNIVERSE) {
    write_attribute(file_id, "domain_type", "universe");
  }

  for (int i = 0; i < domain_ids_.size(); ++i) {
    hid_t group_id =
      create_group(file_id, fmt::format("domain_{}", domain_ids_[i]));

    // Write volume for domain
    const auto& result {results.nuc_results[i]};
    write_dataset(group_id, "volume", results.vol_tallies[i][0].score_acc);

    // Create array of nuclide names from the vector
    auto n_nuc = result.nuclides.size();

    vector<std::string> nucnames;
    for (int i_nuc : result.nuclides) {
      nucnames.push_back(settings::run_CE ? data::nuclides[i_nuc]->name_
                                          : data::mg.nuclides_[i_nuc].name);
    }

    // Create array of total # of atoms with uncertainty for each nuclide
    xt::xtensor<double, 2> atom_data({n_nuc, 2});
    xt::view(atom_data, xt::all(), 0) = xt::adapt(result.atoms);
    xt::view(atom_data, xt::all(), 1) = xt::adapt(result.uncertainty);

    // Write results
    write_dataset(group_id, "nuclides", nucnames);
    write_dataset(group_id, "atoms", atom_data);

    close_group(group_id);
  }

  file_close(file_id);
}

void VolumeCalculation::check_hit(const int32_t i_material,
  const double contrib, vector<VolTally>& vol_tallies) const
{
  // Contribute to entire domain result tally
  vol_tallies[0].score += contrib;

  // Check if this material was previously hit and if so, contribute score
  for (int j = 1; j < vol_tallies.size(); j++) {
    if (vol_tallies[j].index == i_material) {
      vol_tallies[j].score += contrib;
      return;
    }
  }

  // The material was not previously hit, append an entry to the material
  // indices and hits lists
  vol_tallies.push_back(VolTally(i_material, contrib));
}

void VolumeCalculation::reduce_results(
  const CalcResults& local_results, CalcResults& results) const
{
  auto n_threads = num_threads();

  // Collect scalar variables
#pragma omp for ordered schedule(static, 1)
  for (int i = 0; i < n_threads; ++i) {
#pragma omp ordered
    {
      results.n_samples += local_results.n_samples;
      results.n_rays += local_results.n_rays;
      results.n_segs += local_results.n_segs;
      results.n_errors += local_results.n_errors;
      results.cost += local_results.cost;
    }
  }

  // Collect vectored domain-wise results
  for (int i_domain = 0; i_domain < domain_ids_.size(); ++i_domain) {
    const vector<VolTally>& local_vol_tall =
      local_results.vol_tallies[i_domain];
    vector<VolTally>& vol_tall = results.vol_tallies[i_domain];

#pragma omp for ordered schedule(static, 1)
    for (int i = 0; i < n_threads; ++i) {
#pragma omp ordered
      for (int j = 0; j < local_vol_tall.size(); ++j) {
        // Check if this material has been added to the master list and if
        // so, accumulate scores
        const auto ind {local_vol_tall[j].index};
        const auto it = std::find_if(vol_tall.begin(), vol_tall.end(),
          [ind](const VolTally& vt) { return vt.index == ind; });
        if (it == vol_tall.end()) {
          vol_tall.push_back(local_vol_tall[j]);
        } else {
          vol_tall[it - vol_tall.begin()].append_tally(local_vol_tall[j]);
        }
      }
    }
  }
}

#ifdef OPENMC_MPI
void VolumeCalculation::initialize_MPI_struct() const
{
  // This code is a slightly modified replica of initialize_mpi() from
  // initialize.cpp. It works under GCC in the Release configuration, but not
  // sure that the adress offsets of structure's memebrs should be necessary the
  // same everywhere as using an optimizing compiler.
  CalcResults cr(*this);
  MPI_Aint cr_disp[5], cr_d;
  MPI_Get_address(&cr, &cr_d);
  MPI_Get_address(&cr.n_errors, &cr_disp[0]);
  MPI_Get_address(&cr.n_rays, &cr_disp[1]);
  MPI_Get_address(&cr.n_segs, &cr_disp[2]);
  MPI_Get_address(&cr.n_samples, &cr_disp[3]);
  MPI_Get_address(&cr.cost, &cr_disp[4]);
  for (int i = 0; i < 5; i++) {
    cr_disp[i] -= cr_d;
  }

  int cr_blocks[] {1, 1, 1, 1, 1};
  MPI_Datatype cr_types[] {
    MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE};
  MPI_Type_create_struct(5, cr_blocks, cr_disp, cr_types, &mpi_vol_results);
  MPI_Type_commit(&mpi_vol_results);

  VolTally vt;
  MPI_Aint vt_disp[3], vt_d;
  MPI_Get_address(&vt, &vt_d);
  MPI_Get_address(&vt.score, &vt_disp[0]);
  MPI_Get_address(&vt.score_acc, &vt_disp[1]);
  MPI_Get_address(&vt.index, &vt_disp[2]);
  for (int i = 0; i < 3; i++) {
    vt_disp[i] -= vt_d;
  }

  int vt_blocks[] {1, 2, 1};
  MPI_Datatype vt_types[] {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  MPI_Type_create_struct(3, vt_blocks, vt_disp, vt_types, &mpi_volume_tally);
  MPI_Type_commit(&mpi_volume_tally);
}

void VolumeCalculation::delete_MPI_struct() const
{
  MPI_Type_free(&mpi_vol_results);
  MPI_Type_free(&mpi_volume_tally);
}

void VolumeCalculation::CalcResults::collect_MPI()
{
  vector<int> domain_sizes(vol_tallies.size());

  if (!mpi::master) {
    // n_domain + 2 MPI messages will be send in total to node mpi::master

    // To determination of an unique tag for each MPI message (as below)
    int mpi_offset = mpi::rank * (vol_tallies.size() + 2) + 2;

    // Pass root data of the struct
    MPI_Send(
      (void*)this, 1, mpi_vol_results, 0, mpi_offset - 2, mpi::intracomm);

    // Pass sizes of domain-wise data
    for (int i_domain = 0; i_domain < vol_tallies.size(); i_domain++) {
      domain_sizes[i_domain] = vol_tallies[i_domain].size();
    }

    MPI_Send(domain_sizes.data(), domain_sizes.size(), MPI_INT, 0,
      mpi_offset - 1, mpi::intracomm);

    // Pass domain-wise data of struct
    for (int i_domain = 0; i_domain < vol_tallies.size(); i_domain++) {
      MPI_Send(vol_tallies[i_domain].data(), domain_sizes[i_domain],
        mpi_volume_tally, 0, mpi_offset + i_domain, mpi::intracomm);
    }

    this->reset(); // Delete passed to main process data

  } else {

    // n_domain + 2 MPI messages will be recieved in total on node mpi::master

    for (int i_proc = 1; i_proc < mpi::n_procs; i_proc++) {

      CalcResults res_buff(*this); // temporary storage for recived data

      // To determination of an unique tag for each MPI message (as above)
      int mpi_offset = i_proc * (vol_tallies.size() + 2) + 2;

      // Pass root data of struct
      MPI_Recv(&res_buff, 1, mpi_vol_results, i_proc, mpi_offset - 2,
        mpi::intracomm, MPI_STATUS_IGNORE);

      // Pass sizes of domain-wise data
      MPI_Recv(domain_sizes.data(), domain_sizes.size(), MPI_INT, i_proc,
        mpi_offset - 1, mpi::intracomm, MPI_STATUS_IGNORE);

      // Pass domain-wise data of struct
      for (int i_domain = 0; i_domain < vol_tallies.size(); i_domain++) {
        res_buff.vol_tallies[i_domain].resize(domain_sizes[i_domain]);
        MPI_Recv(res_buff.vol_tallies[i_domain].data(), domain_sizes[i_domain],
          mpi_volume_tally, i_proc, mpi_offset + i_domain, mpi::intracomm,
          MPI_STATUS_IGNORE);
      }

      this->append(res_buff);
    }
  }
}
#endif

VolumeCalculation::CalcResults::CalcResults(const VolumeCalculation& vol_calc)
{
  n_samples = 0;
  n_rays = 0;
  n_segs = 0;
  n_errors = 0;
  iterations = 0;
  cost = 0.;
  for (int i = 0; i < vol_calc.domain_ids_.size(); i++) {
    vector<VolTally> vt_vect; // Tally group for a domain
    vt_vect.push_back(
      VolTally(_INDEX_TOTAL)); // Zero-element tally for entire domain totals
    vol_tallies.push_back(vt_vect);
    nuc_results.push_back(NuclResult());
  }
}

void VolumeCalculation::CalcResults::reset()
{
  n_samples = 0;
  n_rays = 0;
  n_segs = 0;
  n_errors = 0;
  cost = 0.;
  for (auto& vt : vol_tallies) {
    std::fill(vt.begin(), vt.end(), VolTally());
  }
}

void VolumeCalculation::CalcResults::append(const CalcResults& other)
{
  n_samples += other.n_samples;
  n_rays += other.n_rays;
  n_segs += other.n_segs;
  n_errors += other.n_errors;
  cost += other.cost;

  // The domain-wise vectors this.vol_tallies and other.vol_tallies are
  // conformed each to other by definition
  for (auto id = 0; id < vol_tallies.size(); id++) {
    // Merging current domain vector from other.vol_tallies into this via
    // pair-wise comparisons
    for (const auto& vt_other : other.vol_tallies[id]) {
      bool already_appended = false;
      for (auto& vt : vol_tallies[id]) {
        if (vt.index == vt_other.index) {
          vt.append_tally(vt_other);
          already_appended = true;
          break;
        }
      }
      if (!already_appended)
        vol_tallies[id].push_back(vt_other);
    }
  }
}

inline VolumeCalculation::VolTally::VolTally(const int i_material,
  const double contrib, const double score_acc_, const double score2_acc_)
{
  score = contrib;
  score_acc[0] = score_acc_;
  score_acc[1] = score2_acc_;
  index = i_material;
}

inline void VolumeCalculation::VolTally::finalize_batch(
  const double batch_size_1)
{
  if (score != 0.) {
    score *= batch_size_1;
    score_acc[0] += score;
    score_acc[1] += score * score;
    score = 0.;
  }
}

inline void VolumeCalculation::VolTally::assign_tally(const VolTally& vol_tally)
{
  score = vol_tally.score;
  score_acc = vol_tally.score_acc;
  index = vol_tally.index;
}

inline void VolumeCalculation::VolTally::append_tally(const VolTally& vol_tally)
{
  score += vol_tally.score;
  score_acc[0] += vol_tally.score_acc[0];
  score_acc[1] += vol_tally.score_acc[1];
}

inline bool VolumeCalculation::VolTally::trigger_state(
  const TriggerMetric trigger_type, const double threshold,
  const size_t& n_samples) const
{
  // For sample contribution to volume fraction limited by 1, the maximal
  // allowed n_samples value is around 1.e102, but this is still much larger
  // than the size_t limit equal to ~1.8e19
  const double ns1 = static_cast<double>(n_samples - 1);
  const double ns = static_cast<double>(n_samples);
  const double mean_xi_sq = score_acc[0] * score_acc[0];

  // Adjusted sample variance: s^2 = (\bar{x^2} - \bar{x}^2) / (N-1)
  // \bar{x}^2 = mean_xi_sq / N^2, \bar{\x^2} = score_acc[1] / N, N = ns
  switch (trigger_type) {
  case TriggerMetric::variance:
  case TriggerMetric::standard_deviation:
    // Condition: s^2 / N < t'
    // Equivalent implementation:
    // N^2 * (\bar{x^2} - \bar{x}^2) < t' * (N-1) * N^2
    return score_acc[1] * ns - mean_xi_sq < threshold * ns1 * ns * ns;
  case TriggerMetric::relative_error:
    // Condition: (s^2 / \mu^2) / N < t'
    // Equivalent implementation:
    // N^2 * (\bar{x^2} - \bar{x}^2) < t' * (N-1) * (N * \bar{x})^2
    return score_acc[1] * ns - mean_xi_sq < threshold * ns1 * mean_xi_sq;
  default:
    return true;
  }
}

array<double, 2> VolumeCalculation::get_tally_results(const size_t& n_samples,
  const double coeff_norm, const VolTally& vol_tally) const
{
  array<double, 2> volume;
  const double ns_1 = 1. / static_cast<double>(n_samples);
  volume[0] = vol_tally.score_acc[0] * ns_1;
  volume[1] = vol_tally.score_acc[1] * ns_1;
  volume[1] = std::sqrt(
    (volume[1] - volume[0] * volume[0]) / static_cast<double>(n_samples - 1));
  volume[0] *= coeff_norm;
  volume[1] *= coeff_norm;
  return volume;
}

std::pair<double, double> VolumeCalculation::get_box_chord(
  const Position& r, const Direction& u) const
{
  // Compute distanses to each box plane orthogonal to an axis
  Direction u_1 = {1., 1., 1.};
  u_1 = u_1 / u;
  Position xi = (lower_left_ - r) * u_1;
  const array<double, 3> dist1 = {xi.x, xi.y, xi.z};
  xi = (upper_right_ - r) * u_1;
  const array<double, 3> dist2 = {xi.x, xi.y, xi.z};

  // Find the minimal forward (positive values) and backward (negative values)
  // distances across the computed ones (probably there is some STL
  // alternatives)
  std::pair<double, double> chord_lengths {std::minmax(dist1[0], dist2[0])};
  for (int i = 1; i < 3; i++) {
    const std::pair<double, double> dist_mm = std::minmax(dist1[i], dist2[i]);
    chord_lengths.first = std::max(chord_lengths.first, dist_mm.first);
    chord_lengths.second = std::min(chord_lengths.second, dist_mm.second);
  }
  return chord_lengths;
}

void VolEstRay::on_intersection()
{
  if (traversal_distance_ == 0.) {
    return; // No crossing model
  }

  results_.n_segs++;
  if (traversal_distance_last_ == 0.) {
    results_.n_rays++; // First segment of new ray
  }

  // At this point, current GeometryState parameters represent the cell behind
  // crossed surface, but the segment length is known for the previous
  // cell only, therefore we use below the last cell keept in GeometryState
  const auto score = (std::min(traversal_distance_, traversal_distance_max_) -
                       traversal_distance_last_) *
                     coeff_mult_;

  //----------------------------------------------------------------------------
  // Tracing error diagnostic
  // TODO: those can be implemented here clever diagnostic and guidance for
  // user to fix input mistakes
  //----------------------------------------------------------------------------

  if (traversal_distance_ >= traversal_distance_max_) {
    stop();
  } else {
    traversal_distance_last_ = traversal_distance_;
  }

  // In a case of single-segment ray (leakage after 1st surface crossing),
  // current segment material ID seems to be contained in material(), but in
  // other cases it is in material_last(). Due to this unclear behavior, it is
  // used below a more fundamental way for material determination -- via the
  // stable last cell record.
  vol_scoring(VolumeCalculation::EstMode::RAYTRACE, score,
    (current_cell(VolumeCalculation::EstMode::RAYTRACE, n_coord_last() - 1)
        ->material(n_coord_last() - 1)));
}

void VolEstRay::score_hit()
{
  if (exhaustive_find_cell(*this))
    vol_scoring(VolumeCalculation::EstMode::REJECTION, 1.,
      material()); // One hit score
}

void VolEstRay::vol_scoring(
  const VolumeCalculation::EstMode mode, const double score, const int id_mat)
{
  const auto n_domains = vol_calc_.domain_ids_.size();

  switch (vol_calc_.domain_type_) {
  case VolumeCalculation::TallyDomain::MATERIAL:
    if (id_mat != MATERIAL_VOID) {
      for (auto i_domain = 0; i_domain < n_domains; i_domain++) {
        if (model::materials[id_mat]->id_ == vol_calc_.domain_ids_[i_domain]) {
          vol_calc_.check_hit(id_mat, score, results_.vol_tallies[i_domain]);
          break;
        }
      }
    }
    break;
  case VolumeCalculation::TallyDomain::CELL:
    for (auto level = 0; level < n_coord_last(); ++level) {
      for (auto i_domain = 0; i_domain < n_domains; i_domain++) {
        if (current_cell(mode, level)->id_ == vol_calc_.domain_ids_[i_domain]) {
          vol_calc_.check_hit(id_mat, score, results_.vol_tallies[i_domain]);
          break;
        }
      }
    }
    break;
  case VolumeCalculation::TallyDomain::UNIVERSE:
    for (auto level = 0; level < n_coord_last(); ++level) {
      for (auto i_domain = 0; i_domain < n_domains; ++i_domain) {
        if (model::universes[current_cell(mode, level)->universe_]->id_ ==
            vol_calc_.domain_ids_[i_domain]) {
          vol_calc_.check_hit(id_mat, score, results_.vol_tallies[i_domain]);
          break;
        }
      }
    }
  }
}

inline std::unique_ptr<Cell>& VolEstRay::current_cell(
  VolumeCalculation::EstMode mode, int level)
{
  switch (mode) {
  case VolumeCalculation::EstMode::REJECTION:
    return model::cells[coord(level).cell()]; // Current position
  case VolumeCalculation::EstMode::RAYTRACE:
    return model::cells[cell_last(level)]; // Previous segment
  }
}

void free_memory_volume()
{
  openmc::model::volume_calcs.clear();
}

} // namespace openmc

//==============================================================================
// OPENMC_CALCULATE_VOLUMES runs each of the stochastic volume calculations
// that the user has specified and writes results to HDF5 files
//==============================================================================

int openmc_calculate_volumes()
{
  using namespace openmc;

  if (mpi::master) {
    header("STOCHASTIC VOLUME CALCULATION", 3);
  }
  Timer time_volume;
  time_volume.start();

  for (int i = 0; i < model::volume_calcs.size(); ++i) {
    write_message(4, "Running volume calculation {}", i + 1);

    // Run volume calculation
    const auto& vol_calc {model::volume_calcs[i]};
    VolumeCalculation::CalcResults results(vol_calc);
    try {
      vol_calc.execute(results);
    } catch (const std::exception& e) {
      set_errmsg(e.what());
      return OPENMC_E_UNASSIGNED;
    }

    if (mpi::master) {

      // Output volume calculation results and statistics
      vol_calc.show_results(results);

      // Write volumes to HDF5 file
      std::string filename =
        fmt::format("{}volume_{}.h5", settings::path_output, i + 1);
      vol_calc.to_hdf5(filename, results);
    }
  }

  // Show elapsed time
  time_volume.stop();
  write_message(6, "Elapsed time: {} sec", time_volume.elapsed());

  return 0;
}
