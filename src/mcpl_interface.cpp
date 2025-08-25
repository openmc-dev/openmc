#include "openmc/mcpl_interface.h"

#include "openmc/bank.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/message_passing.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/state_point.h"
#include "openmc/vector.h"

#include <fmt/core.h>

#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// WARNING: These declarations MUST EXACTLY MATCH the structure and function
// signatures of the libmcpl being loaded at runtime. Any discrepancy will
// likely lead to crashes or incorrect behavior. This is a maintenance risk.
// MCPL 2.2.0

#pragma pack(push, 1)
struct mcpl_particle_repr_t {
  double ekin;
  double polarisation[3];
  double position[3];
  double direction[3];
  double time;
  double weight;
  int32_t pdgcode;
  uint32_t userflags;
};
#pragma pack(pop)

// Opaque struct definitions replicating the MCPL C-API to ensure ABI
// compatibility without including mcpl.h. These must be kept in sync.
struct mcpl_file_t {
  void* internal;
};
struct mcpl_outfile_t {
  void* internal;
};

// Function pointer types for the dynamically loaded MCPL library
using mcpl_open_file_fpt = mcpl_file_t* (*)(const char* filename);
using mcpl_hdr_nparticles_fpt = uint64_t (*)(mcpl_file_t* file_handle);
using mcpl_read_fpt = const mcpl_particle_repr_t* (*)(mcpl_file_t* file_handle);
using mcpl_close_file_fpt = void (*)(mcpl_file_t* file_handle);

using mcpl_hdr_add_data_fpt = void (*)(mcpl_outfile_t* file_handle,
  const char* key, int32_t ldata, const char* data);
using mcpl_create_outfile_fpt = mcpl_outfile_t* (*)(const char* filename);
using mcpl_hdr_set_srcname_fpt = void (*)(
  mcpl_outfile_t* outfile_handle, const char* srcname);
using mcpl_add_particle_fpt = void (*)(
  mcpl_outfile_t* outfile_handle, const mcpl_particle_repr_t* particle);
using mcpl_close_outfile_fpt = void (*)(mcpl_outfile_t* outfile_handle);
using mcpl_hdr_add_stat_sum_fpt = void (*)(
  mcpl_outfile_t* outfile_handle, const char* key, double value);

namespace openmc {

#ifdef _WIN32
using LibraryHandleType = HMODULE;
#else
using LibraryHandleType = void*;
#endif

std::string get_last_library_error()
{
#ifdef _WIN32
  DWORD error_code = GetLastError();
  if (error_code == 0)
    return "No error reported by system."; // More accurate than "No error."
  LPSTR message_buffer = nullptr;
  size_t size =
    FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                     FORMAT_MESSAGE_IGNORE_INSERTS,
      NULL, error_code, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
      (LPSTR)&message_buffer, 0, NULL);
  std::string message(message_buffer, size);
  LocalFree(message_buffer);
  while (
    !message.empty() && (message.back() == '\n' || message.back() == '\r')) {
    message.pop_back();
  }
  return message;
#else
  const char* err = dlerror();
  return err ? std::string(err) : "No error reported by dlerror.";
#endif
}

struct McplApi {
  mcpl_open_file_fpt open_file;
  mcpl_hdr_nparticles_fpt hdr_nparticles;
  mcpl_read_fpt read;
  mcpl_close_file_fpt close_file;
  mcpl_create_outfile_fpt create_outfile;
  mcpl_hdr_set_srcname_fpt hdr_set_srcname;
  mcpl_hdr_add_data_fpt hdr_add_data;
  mcpl_add_particle_fpt add_particle;
  mcpl_close_outfile_fpt close_outfile;
  mcpl_hdr_add_stat_sum_fpt hdr_add_stat_sum;

  explicit McplApi(LibraryHandleType lib_handle)
  {
    if (!lib_handle)
      throw std::runtime_error(
        "MCPL library handle is null during API binding.");

    auto load_symbol_platform = [lib_handle](const char* name) {
      void* sym = nullptr;
#ifdef _WIN32
      sym = (void*)GetProcAddress(lib_handle, name);
#else
      sym = dlsym(lib_handle, name);
#endif
      if (!sym) {
        throw std::runtime_error(
          fmt::format("Failed to load MCPL symbol '{}': {}", name,
            get_last_library_error()));
      }
      return sym;
    };

    open_file = reinterpret_cast<mcpl_open_file_fpt>(
      load_symbol_platform("mcpl_open_file"));
    hdr_nparticles = reinterpret_cast<mcpl_hdr_nparticles_fpt>(
      load_symbol_platform("mcpl_hdr_nparticles"));
    read = reinterpret_cast<mcpl_read_fpt>(load_symbol_platform("mcpl_read"));
    close_file = reinterpret_cast<mcpl_close_file_fpt>(
      load_symbol_platform("mcpl_close_file"));
    create_outfile = reinterpret_cast<mcpl_create_outfile_fpt>(
      load_symbol_platform("mcpl_create_outfile"));
    hdr_set_srcname = reinterpret_cast<mcpl_hdr_set_srcname_fpt>(
      load_symbol_platform("mcpl_hdr_set_srcname"));
    hdr_add_data = reinterpret_cast<mcpl_hdr_add_data_fpt>(
      load_symbol_platform("mcpl_hdr_add_data"));
    add_particle = reinterpret_cast<mcpl_add_particle_fpt>(
      load_symbol_platform("mcpl_add_particle"));
    close_outfile = reinterpret_cast<mcpl_close_outfile_fpt>(
      load_symbol_platform("mcpl_close_outfile"));

    // Try to load mcpl_hdr_add_stat_sum (available in MCPL >= 2.1.0)
    // Set to nullptr if not available for graceful fallback
    try {
      hdr_add_stat_sum = reinterpret_cast<mcpl_hdr_add_stat_sum_fpt>(
        load_symbol_platform("mcpl_hdr_add_stat_sum"));
    } catch (const std::runtime_error&) {
      hdr_add_stat_sum = nullptr;
    }
  }
};

static LibraryHandleType g_mcpl_lib_handle = nullptr;
static std::unique_ptr<McplApi> g_mcpl_api;
static bool g_mcpl_init_attempted = false;
static bool g_mcpl_successfully_loaded = false;
static std::string g_mcpl_load_error_msg;
static std::once_flag g_mcpl_init_flag;

void append_error(std::string& existing_msg, const std::string& new_error)
{
  if (!existing_msg.empty()) {
    existing_msg += "; ";
  }
  existing_msg += new_error;
}

void initialize_mcpl_interface_impl()
{
  g_mcpl_init_attempted = true;
  g_mcpl_load_error_msg.clear();

  // Try mcpl-config
  if (!g_mcpl_lib_handle) {
    FILE* pipe = nullptr;
#ifdef _WIN32
    pipe = _popen("mcpl-config --show libpath", "r");
#else
    pipe = popen("mcpl-config --show libpath 2>/dev/null", "r");
#endif
    if (pipe) {
      char buffer[512];
      if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        std::string shlibpath = buffer;
        // Remove trailing whitespace
        while (!shlibpath.empty() &&
               std::isspace(static_cast<unsigned char>(shlibpath.back()))) {
          shlibpath.pop_back();
        }

        if (!shlibpath.empty()) {
#ifdef _WIN32
          g_mcpl_lib_handle = LoadLibraryA(shlibpath.c_str());
#else
          g_mcpl_lib_handle = dlopen(shlibpath.c_str(), RTLD_LAZY);
#endif
          if (!g_mcpl_lib_handle) {
            append_error(
              g_mcpl_load_error_msg, fmt::format("From mcpl-config ({}): {}",
                                       shlibpath, get_last_library_error()));
          }
        }
      }
#ifdef _WIN32
      _pclose(pipe);
#else
      pclose(pipe);
#endif
    } else { // pipe failed to open
      append_error(g_mcpl_load_error_msg,
        "mcpl-config command not found or failed to execute");
    }
  }

  // Try standard library names
  if (!g_mcpl_lib_handle) {
#ifdef _WIN32
    const char* standard_names[] = {"mcpl.dll", "libmcpl.dll"};
#else
    const char* standard_names[] = {"libmcpl.so", "libmcpl.dylib"};
#endif
    for (const char* name : standard_names) {
#ifdef _WIN32
      g_mcpl_lib_handle = LoadLibraryA(name);
#else
      g_mcpl_lib_handle = dlopen(name, RTLD_LAZY);
#endif
      if (g_mcpl_lib_handle)
        break;
    }
    if (!g_mcpl_lib_handle) {
      append_error(
        g_mcpl_load_error_msg, fmt::format("Using standard names (e.g. {}): {}",
                                 standard_names[0], get_last_library_error()));
    }
  }

  if (!g_mcpl_lib_handle) {
    if (mpi::master) {
      warning(fmt::format("MCPL library could not be loaded. MCPL-dependent "
                          "features will be unavailable. Load attempts: {}",
        g_mcpl_load_error_msg.empty()
          ? "No specific error during load attempts."
          : g_mcpl_load_error_msg));
    }
    g_mcpl_successfully_loaded = false;
    return;
  }

  try {
    g_mcpl_api = std::make_unique<McplApi>(g_mcpl_lib_handle);
    g_mcpl_successfully_loaded = true;
    // Do not call dlclose/FreeLibrary at exit. Leaking the handle is safer
    // and standard practice for libraries used for the application's lifetime.
  } catch (const std::runtime_error& e) {
    append_error(g_mcpl_load_error_msg,
      fmt::format(
        "MCPL library loaded, but failed to bind symbols: {}", e.what()));
    if (mpi::master) {
      warning(g_mcpl_load_error_msg);
    }
#ifdef _WIN32
    FreeLibrary(g_mcpl_lib_handle);
#else
    dlclose(g_mcpl_lib_handle);
#endif
    g_mcpl_lib_handle = nullptr;
    g_mcpl_successfully_loaded = false;
  }
}

void initialize_mcpl_interface_if_needed()
{
  std::call_once(g_mcpl_init_flag, initialize_mcpl_interface_impl);
}

bool is_mcpl_interface_available()
{
  initialize_mcpl_interface_if_needed();
  return g_mcpl_successfully_loaded;
}

inline void ensure_mcpl_ready_or_fatal()
{
  initialize_mcpl_interface_if_needed();
  if (!g_mcpl_successfully_loaded) {
    fatal_error("MCPL functionality is required, but the MCPL library is not "
                "available or failed to initialize. Please ensure MCPL is "
                "installed and its library can be found (e.g., via PATH on "
                "Windows, LD_LIBRARY_PATH on Linux, or DYLD_LIBRARY_PATH on "
                "macOS). You can often install MCPL with 'pip install mcpl' or "
                "'conda install mcpl'.");
  }
}

SourceSite mcpl_particle_to_site(const mcpl_particle_repr_t* particle_repr)
{
  SourceSite site;
  switch (particle_repr->pdgcode) {
  case 2112:
    site.particle = ParticleType::neutron;
    break;
  case 22:
    site.particle = ParticleType::photon;
    break;
  case 11:
    site.particle = ParticleType::electron;
    break;
  case -11:
    site.particle = ParticleType::positron;
    break;
  default:
    fatal_error(fmt::format(
      "MCPL: Encountered unexpected PDG code {} when converting to SourceSite.",
      particle_repr->pdgcode));
    break;
  }

  // Copy position and direction
  site.r.x = particle_repr->position[0];
  site.r.y = particle_repr->position[1];
  site.r.z = particle_repr->position[2];
  site.u.x = particle_repr->direction[0];
  site.u.y = particle_repr->direction[1];
  site.u.z = particle_repr->direction[2];
  // MCPL stores kinetic energy in [MeV], time in [ms]
  site.E = particle_repr->ekin * 1e6;
  site.time = particle_repr->time * 1e-3;
  site.wgt = particle_repr->weight;
  return site;
}

vector<SourceSite> mcpl_source_sites(std::string path)
{
  ensure_mcpl_ready_or_fatal();
  vector<SourceSite> sites;

  mcpl_file_t* mcpl_file = g_mcpl_api->open_file(path.c_str());
  if (!mcpl_file) {
    fatal_error(fmt::format("MCPL: Could not open file '{}'. It might be "
                            "missing, inaccessible, or not a valid MCPL file.",
      path));
  }

  size_t n_particles_in_file = g_mcpl_api->hdr_nparticles(mcpl_file);
  size_t n_skipped = 0;
  if (n_particles_in_file > 0) {
    sites.reserve(n_particles_in_file);
  }

  for (size_t i = 0; i < n_particles_in_file; ++i) {
    const mcpl_particle_repr_t* p_repr = g_mcpl_api->read(mcpl_file);
    if (!p_repr) {
      warning(fmt::format("MCPL: Read error or unexpected end of file '{}' "
                          "after reading {} of {} expected particles.",
        path, sites.size(), n_particles_in_file));
      break;
    }
    if (p_repr->pdgcode == 2112 || p_repr->pdgcode == 22 ||
        p_repr->pdgcode == 11 || p_repr->pdgcode == -11) {
      sites.push_back(mcpl_particle_to_site(p_repr));
    } else {
      n_skipped++;
    }
  }

  g_mcpl_api->close_file(mcpl_file);

  if (n_skipped > 0 && n_particles_in_file > 0) {
    double percent_skipped =
      100.0 * static_cast<double>(n_skipped) / n_particles_in_file;
    warning(fmt::format(
      "MCPL: Skipped {} of {} total particles ({:.1f}%) in file '{}' because "
      "their type is not supported by OpenMC.",
      n_skipped, n_particles_in_file, percent_skipped, path));
  }

  if (sites.empty()) {
    if (n_particles_in_file > 0) {
      fatal_error(fmt::format(
        "MCPL file '{}' contained {} particles, but none were of the supported "
        "types (neutron, photon, electron, positron). OpenMC cannot proceed "
        "without source particles.",
        path, n_particles_in_file));
    } else {
      fatal_error(fmt::format(
        "MCPL file '{}' is empty or contains no particle data.", path));
    }
  }
  return sites;
}

void write_mcpl_source_bank_internal(mcpl_outfile_t* file_id,
  span<SourceSite> local_source_bank,
  const vector<int64_t>& bank_index_all_ranks)
{
  if (mpi::master) {
    if (!file_id) {
      fatal_error("MCPL: Internal error - master rank called "
                  "write_mcpl_source_bank_internal with null file_id.");
    }
    vector<SourceSite> receive_buffer;

    for (int rank_idx = 0; rank_idx < mpi::n_procs; ++rank_idx) {
      size_t num_sites_on_rank = static_cast<size_t>(
        bank_index_all_ranks[rank_idx + 1] - bank_index_all_ranks[rank_idx]);
      if (num_sites_on_rank == 0)
        continue;

      span<const SourceSite> sites_to_write;
#ifdef OPENMC_MPI
      if (rank_idx == mpi::rank) {
        sites_to_write = openmc::span<const SourceSite>(
          local_source_bank.data(), num_sites_on_rank);
      } else {
        if (receive_buffer.size() < num_sites_on_rank) {
          receive_buffer.resize(num_sites_on_rank);
        }
        MPI_Recv(receive_buffer.data(), num_sites_on_rank, mpi::source_site,
          rank_idx, rank_idx, mpi::intracomm, MPI_STATUS_IGNORE);
        sites_to_write = openmc::span<const SourceSite>(
          receive_buffer.data(), num_sites_on_rank);
      }
#else
      sites_to_write = openmc::span<const SourceSite>(
        local_source_bank.data(), num_sites_on_rank);
#endif
      for (const auto& site : sites_to_write) {
        mcpl_particle_repr_t p_repr {};
        p_repr.position[0] = site.r.x;
        p_repr.position[1] = site.r.y;
        p_repr.position[2] = site.r.z;
        p_repr.direction[0] = site.u.x;
        p_repr.direction[1] = site.u.y;
        p_repr.direction[2] = site.u.z;
        p_repr.ekin = site.E * 1e-6;
        p_repr.time = site.time * 1e3;
        p_repr.weight = site.wgt;
        switch (site.particle) {
        case ParticleType::neutron:
          p_repr.pdgcode = 2112;
          break;
        case ParticleType::photon:
          p_repr.pdgcode = 22;
          break;
        case ParticleType::electron:
          p_repr.pdgcode = 11;
          break;
        case ParticleType::positron:
          p_repr.pdgcode = -11;
          break;
        default:
          continue;
        }
        g_mcpl_api->add_particle(file_id, &p_repr);
      }
    }
  } else {
#ifdef OPENMC_MPI
    if (!local_source_bank.empty()) {
      MPI_Send(local_source_bank.data(), local_source_bank.size(),
        mpi::source_site, 0, mpi::rank, mpi::intracomm);
    }
#endif
  }
}

void write_mcpl_source_point(const char* filename, span<SourceSite> source_bank,
  const vector<int64_t>& bank_index)
{
  ensure_mcpl_ready_or_fatal();

  std::string filename_(filename);
  const auto extension = get_file_extension(filename_);
  if (extension.empty()) {
    filename_.append(".mcpl");
  } else if (extension != "mcpl") {
    warning(fmt::format("Specified filename '{}' has an extension '.{}', but "
                        "an MCPL file (.mcpl) will be written using this name.",
      filename, extension));
  }

  mcpl_outfile_t* file_id = nullptr;

  if (mpi::master) {
    file_id = g_mcpl_api->create_outfile(filename_.c_str());
    if (!file_id) {
      fatal_error(fmt::format(
        "MCPL: Failed to create output file '{}'. Check permissions and path.",
        filename_));
    }
    std::string src_line;
    if (VERSION_DEV) {
      src_line = fmt::format("OpenMC {}.{}.{}-dev{}", VERSION_MAJOR,
        VERSION_MINOR, VERSION_RELEASE, VERSION_COMMIT_COUNT);
    } else {
      src_line = fmt::format(
        "OpenMC {}.{}.{}", VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE);
    }
    g_mcpl_api->hdr_set_srcname(file_id, src_line.c_str());

    // Initialize stat:sum with -1 to indicate incomplete file (issue #3514)
    // This follows MCPL >= 2.1.0 convention for tracking simulation statistics
    // The -1 value indicates "not available" if file creation is interrupted
    if (g_mcpl_api->hdr_add_stat_sum) {
      // Using key "openmc_np1" following tkittel's recommendation
      // Initial value of -1 prevents misleading values in case of crashes
      g_mcpl_api->hdr_add_stat_sum(file_id, "openmc_np1", -1.0);
    }
  }

  write_mcpl_source_bank_internal(file_id, source_bank, bank_index);

  if (mpi::master) {
    if (file_id) {
      // Update stat:sum with actual particle count before closing (issue #3514)
      // This represents the original number of source particles in the
      // simulation (not the number of particles in the file)
      if (g_mcpl_api->hdr_add_stat_sum) {
        // Calculate total source particles from active batches
        // Per issue #3514: this should be the original number of source
        // particles, not the number written to the file
        int64_t total_source_particles =
          static_cast<int64_t>(settings::n_batches - settings::n_inactive) *
          settings::gen_per_batch * settings::n_particles;
        // Update with actual count - this overwrites the initial -1 value
        g_mcpl_api->hdr_add_stat_sum(
          file_id, "openmc_np1", static_cast<double>(total_source_particles));
      }

      g_mcpl_api->close_outfile(file_id);
    }
  }
}

// Collision track feature with MCPL
void write_mcpl_collision_track_internal(mcpl_outfile_t* file_id,
  span<CollisionTrackSite> collision_track_bank,
  const vector<int64_t>& bank_index_all_ranks)
{
  if (mpi::master) {
    if (!file_id) {
      fatal_error("MCPL: Internal error - master rank called "
                  "write_mcpl_source_bank_internal with null file_id.");
    }
    vector<CollisionTrackSite> receive_buffer;

    for (int rank_idx = 0; rank_idx < mpi::n_procs; ++rank_idx) {
      size_t num_sites_on_rank = static_cast<size_t>(
        bank_index_all_ranks[rank_idx + 1] - bank_index_all_ranks[rank_idx]);
      if (num_sites_on_rank == 0)
        continue;

      span<const CollisionTrackSite> sites_to_write;
#ifdef OPENMC_MPI
      if (rank_idx == mpi::rank) {
        sites_to_write = openmc::span<const CollisionTrackSite>(
          collision_track_bank.data(), num_sites_on_rank);
      } else {
        if (receive_buffer.size() < num_sites_on_rank) {
          receive_buffer.resize(num_sites_on_rank);
        }
        MPI_Recv(receive_buffer.data(), num_sites_on_rank,
          mpi::collision_track_site, rank_idx, rank_idx, mpi::intracomm,
          MPI_STATUS_IGNORE);
        sites_to_write = openmc::span<const CollisionTrackSite>(
          receive_buffer.data(), num_sites_on_rank);
      }
#else
      sites_to_write = openmc::span<const CollisionTrackSite>(
        collision_track_bank.data(), num_sites_on_rank);
#endif
      int index = 0;
      for (const auto& site : collision_track_bank) {
        // Binary blob should be added before the mcpl_add_particle function
        std::ostringstream custom_data_stream;
        custom_data_stream << " dE : " << site.dE
                           << " ; event_mt : " << site.event_mt
                           << " ; delayed_group : " << site.delayed_group
                           << " ; cell_id : " << site.cell_id
                           << " ; nuclide_id : " << site.nuclide_id
                           << " ; material_id : " << site.material_id
                           << " ; universe_id : " << site.universe_id
                           << " ; n_collision : " << site.n_collision
                           << " ; parent_id : " << site.parent_id
                           << " ; progeny_id : " << site.progeny_id;

        std::string custom_data_str = custom_data_stream.str();
        std::ostringstream custom_key;
        custom_key << "blob_" << index;
        std::string CK = custom_key.str();
        g_mcpl_api->hdr_add_data(
          file_id, CK.c_str(), custom_data_str.size(), custom_data_str.c_str());
        index++;
      }

      for (const auto& site : sites_to_write) {
        mcpl_particle_repr_t p_repr {};
        p_repr.position[0] = site.r.x;
        p_repr.position[1] = site.r.y;
        p_repr.position[2] = site.r.z;
        p_repr.direction[0] = site.u.x;
        p_repr.direction[1] = site.u.y;
        p_repr.direction[2] = site.u.z;
        p_repr.ekin = site.E * 1e-6;
        p_repr.time = site.time * 1e3;
        p_repr.weight = site.wgt;
        switch (site.particle) {
        case ParticleType::neutron:
          p_repr.pdgcode = 2112;
          break;
        case ParticleType::photon:
          p_repr.pdgcode = 22;
          break;
        case ParticleType::electron:
          p_repr.pdgcode = 11;
          break;
        case ParticleType::positron:
          p_repr.pdgcode = -11;
          break;
        default:
          continue;
        }
        g_mcpl_api->add_particle(file_id, &p_repr);
      }
    }
  } else {
#ifdef OPENMC_MPI
    if (!collision_track_bank.empty()) {
      MPI_Send(collision_track_bank.data(), collision_track_bank.size(),
        mpi::collision_track_site, 0, mpi::rank, mpi::intracomm);
    }
#endif
  }
}

void write_mcpl_collision_track(const char* filename,
  span<CollisionTrackSite> collision_track_bank,
  const vector<int64_t>& bank_index)
{
  ensure_mcpl_ready_or_fatal();

  std::string filename_(filename);
  const auto extension = get_file_extension(filename_);
  if (extension.empty()) {
    filename_.append(".mcpl");
  } else if (extension != "mcpl") {
    warning(fmt::format("Specified filename '{}' has an extension '.{}', but "
                        "an MCPL file (.mcpl) will be written using this name.",
      filename, extension));
  }

  mcpl_outfile_t* file_id = nullptr;

  if (mpi::master) {
    file_id = g_mcpl_api->create_outfile(filename_.c_str());
    if (!file_id) {
      fatal_error(fmt::format(
        "MCPL: Failed to create output file '{}'. Check permissions and path.",
        filename_));
    }
    std::string src_line;
    if (VERSION_DEV) {
      src_line = fmt::format("OpenMC {}.{}.{}-dev{}", VERSION_MAJOR,
        VERSION_MINOR, VERSION_RELEASE, VERSION_COMMIT_COUNT);
    } else {
      src_line = fmt::format(
        "OpenMC {}.{}.{}", VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE);
    }

    g_mcpl_api->hdr_set_srcname(file_id, src_line.c_str());
  }
  write_mcpl_collision_track_internal(
    file_id, collision_track_bank, bank_index);

  if (mpi::master) {
    if (file_id) {
      g_mcpl_api->close_outfile(file_id);
    }
  }
}

} // namespace openmc
