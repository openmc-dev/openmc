#include "openmc/soa_particle.h"
#include "openmc/settings.h"

#include "openmc/geometry.h" // model::n_coord_levels
#include "openmc/nuclide.h"  // data::nuclides
#include "openmc/photon.h"   // data::elements
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"

namespace openmc {

void allocate_soa_data()
{
  using namespace soa;

  // Pick the number of particles in flight
  int particles_in_flight;
  if (settings::event_based) {
    particles_in_flight = settings::max_particles_in_flight;
  } else {
#ifdef _OPENMP
    particles_in_flight = omp_get_max_threads();
#else
    particles_in_flight = 1;
#endif
  }

  // Allocate sufficient room in all arrays
  neutron_xs.resize(data::nuclides.size() * particles_in_flight);
  photon_xs.resize(data::elements.size() * particles_in_flight);
  macro_xs.resize(particles_in_flight);
  id.resize(particles_in_flight);
  type.resize(particles_in_flight);
  n_coord.resize(particles_in_flight);
  cell_instance.resize(particles_in_flight);
  coord.resize(particles_in_flight);
  n_coord_last.resize(particles_in_flight);
  cell_last.resize(particles_in_flight);
  E.resize(particles_in_flight);
  E_last.resize(particles_in_flight);
  g.resize(particles_in_flight);
  g_last.resize(particles_in_flight);
  wgt.resize(particles_in_flight);
  mu.resize(particles_in_flight);
  alive.resize(particles_in_flight);
  r_last_current.resize(particles_in_flight);
  r_last.resize(particles_in_flight);
  u_last.resize(particles_in_flight);
  wgt_last.resize(particles_in_flight);
  wgt_absorb.resize(particles_in_flight);
  fission.resize(particles_in_flight);
  event.resize(particles_in_flight);
  event_nuclide.resize(particles_in_flight);
  event_mt.resize(particles_in_flight);
  delayed_group.resize(particles_in_flight);
  n_bank.resize(particles_in_flight);
  n_bank_second.resize(particles_in_flight);
  wgt_bank.resize(particles_in_flight);
  n_delayed_bank.resize(particles_in_flight * MAX_DELAYED_GROUPS);
  surface.resize(particles_in_flight);
  cell_born.resize(particles_in_flight);
  material.resize(particles_in_flight);
  material_last.resize(particles_in_flight);
  boundary.resize(particles_in_flight);
  sqrtkT.resize(particles_in_flight);
  sqrtkT_last.resize(particles_in_flight);
  n_collision.resize(particles_in_flight);
  write_track.resize(particles_in_flight);
  seeds.resize(N_STREAMS * particles_in_flight);
  stream.resize(particles_in_flight);
  secondary_bank.resize(max_secondary_particles * particles_in_flight);
  secondary_bank_current_indx.resize(particles_in_flight);
  current_work.resize(particles_in_flight);
  flux_derivs.resize(model::tally_derivs.size() * particles_in_flight);
  filter_matches.resize(model::tally_filters.size() * particles_in_flight);
  nu_bank.resize(max_fissions * particles_in_flight);
  nu_bank_current_indx.resize(particles_in_flight);
  keff_tally_absorption.resize(particles_in_flight);
  keff_tally_collision.resize(particles_in_flight);
  keff_tally_tracklength.resize(particles_in_flight);
  keff_tally_leakage.resize(particles_in_flight);
  trace.resize(particles_in_flight);
  collision_distance.resize(particles_in_flight);
  n_event.resize(particles_in_flight);
#ifdef DAGMC
  history.resize(particles_in_flight);
  last_dir.resize(particles_in_flight);
#endif
  n_progeny.resize(particles_in_flight);

  // Cache a few frequently accessed variables.
  n_nuclides = data::nuclides.size();
  n_elements = data::elements.size();
  n_coord_levels = model::n_coord_levels;
  n_tally_derivs = model::tally_derivs.size();
  n_tally_filters = model::tally_filters.size();

#ifdef __CUDACC__
  // Now that data is laid out in memory, we can copy shared memory pointers to
  // device constant memory for fast access from that side. Similarly, there
  // are a few constants to put in shared memory e.g. n_nuclides.
  cudaMemcpyToSymbol(gpu::n_nuclides, n_nuclides, sizeof(int));
  cudaMemcpyToSymbol(gpu::n_elements, n_elements, sizeof(int));
  cudaMemcpyToSymbol(gpu::n_coord_levels, n_coord_levels, sizeof(int));
  cudaMemcpyToSymbol(gpu::n_tally_derivs, n_tally_derivs, sizeof(int));
  cudaMemcpyToSymbol(gpu::n_tally_filters, n_tally_filters, sizeof(int));

  void* tmp;
  tmp = neutron_xs.data();
  cudaMemcpyToSymbol(gpu::neutron_xs, &tmp, sizeof(void*));
  tmp = photon_xs.data();
  cudaMemcpyToSymbol(gpu::photon_xs, &tmp, sizeof(void*));
  tmp = macro_xs.data();
  cudaMemcpyToSymbol(gpu::macro_xs, &tmp, sizeof(void*));
  tmp = id.data();
  cudaMemcpyToSymbol(gpu::id, &tmp, sizeof(void*));
  tmp = type.data();
  cudaMemcpyToSymbol(gpu::type, &tmp, sizeof(void*));
  tmp = n_coord.data();
  cudaMemcpyToSymbol(gpu::n_coord, &tmp, sizeof(void*));
  tmp = coord.data();
  cudaMemcpyToSymbol(gpu::cell_instance, &tmp, sizeof(void*));
  tmp = coord.data();
  cudaMemcpyToSymbol(gpu::coord, &tmp, sizeof(void*));
  tmp = n_coord_last.data();
  cudaMemcpyToSymbol(gpu::n_coord_last, &tmp, sizeof(void*));
  tmp = cell_last.data();
  cudaMemcpyToSymbol(gpu::cell_last, &tmp, sizeof(void*));
  tmp = E.data();
  cudaMemcpyToSymbol(gpu::E, &tmp, sizeof(void*));
  tmp = E_last.data();
  cudaMemcpyToSymbol(gpu::E_last, &tmp, sizeof(void*));
  tmp = g.data();
  cudaMemcpyToSymbol(gpu::g, &tmp, sizeof(void*));
  tmp = g_last.data();
  cudaMemcpyToSymbol(gpu::g_last, &tmp, sizeof(void*));
  tmp = wgt.data();
  cudaMemcpyToSymbol(gpu::wgt, &tmp, sizeof(void*));
  tmp = mu.data();
  cudaMemcpyToSymbol(gpu::mu, &tmp, sizeof(void*));
  tmp = alive.data();
  cudaMemcpyToSymbol(gpu::alive, &tmp, sizeof(void*));
  tmp = r_last_current.data();
  cudaMemcpyToSymbol(gpu::r_last_current, &tmp, sizeof(void*));
  tmp = u_last.data();
  cudaMemcpyToSymbol(gpu::r_last, &tmp, sizeof(void*));
  tmp = u_last.data();
  cudaMemcpyToSymbol(gpu::u_last, &tmp, sizeof(void*));
  tmp = wgt_last.data();
  cudaMemcpyToSymbol(gpu::wgt_last, &tmp, sizeof(void*));
  tmp = wgt_absorb.data();
  cudaMemcpyToSymbol(gpu::wgt_absorb, &tmp, sizeof(void*));
  tmp = fission.data();
  cudaMemcpyToSymbol(gpu::fission, &tmp, sizeof(void*));
  tmp = event.data();
  cudaMemcpyToSymbol(gpu::event, &tmp, sizeof(void*));
  tmp = event_nuclide.data();
  cudaMemcpyToSymbol(gpu::event_nuclide, &tmp, sizeof(void*));
  tmp = event_mt.data();
  cudaMemcpyToSymbol(gpu::event_mt, &tmp, sizeof(void*));
  tmp = delayed_group.data();
  cudaMemcpyToSymbol(gpu::delayed_group, &tmp, sizeof(void*));
  tmp = n_bank.data();
  cudaMemcpyToSymbol(gpu::n_bank, &tmp, sizeof(void*));
  tmp = n_bank_second.data();
  cudaMemcpyToSymbol(gpu::n_bank_second, &tmp, sizeof(void*));
  tmp = wgt_bank.data();
  cudaMemcpyToSymbol(gpu::wgt_bank, &tmp, sizeof(void*));
  tmp = n_delayed_bank.data();
  cudaMemcpyToSymbol(gpu::n_delayed_bank, &tmp, sizeof(void*));
  tmp = surface.data();
  cudaMemcpyToSymbol(gpu::surface, &tmp, sizeof(void*));
  tmp = cell_born.data();
  cudaMemcpyToSymbol(gpu::cell_born, &tmp, sizeof(void*));
  tmp = material.data();
  cudaMemcpyToSymbol(gpu::material, &tmp, sizeof(void*));
  tmp = material_last.data();
  cudaMemcpyToSymbol(gpu::material_last, &tmp, sizeof(void*));
  tmp = boundary.data();
  cudaMemcpyToSymbol(gpu::boundary, &tmp, sizeof(void*));
  tmp = sqrtkT.data();
  cudaMemcpyToSymbol(gpu::sqrtkT, &tmp, sizeof(void*));
  tmp = sqrtkT_last.data();
  cudaMemcpyToSymbol(gpu::sqrtkT_last, &tmp, sizeof(void*));
  tmp = n_collision.data();
  cudaMemcpyToSymbol(gpu::n_collision, &tmp, sizeof(void*));
  tmp = write_track.data();
  cudaMemcpyToSymbol(gpu::write_track, &tmp, sizeof(void*));
  tmp = seeds.data();
  cudaMemcpyToSymbol(gpu::seeds, &tmp, sizeof(void*));
  tmp = stream.data();
  cudaMemcpyToSymbol(gpu::stream, &tmp, sizeof(void*));
  tmp = secondary_bank_current_indx.data();
  cudaMemcpyToSymbol(gpu::secondary_bank_current_indx, &tmp, sizeof(void*));
  tmp = secondary_bank.data();
  cudaMemcpyToSymbol(gpu::secondary_bank, &tmp, sizeof(void*));
  tmp = current_work.data();
  cudaMemcpyToSymbol(gpu::current_work, &tmp, sizeof(void*));
  tmp = flux_derivs.data();
  cudaMemcpyToSymbol(gpu::flux_derivs, &tmp, sizeof(void*));
  tmp = filter_matches.data();
  cudaMemcpyToSymbol(gpu::filter_matches, &tmp, sizeof(void*));
  tmp = nu_bank_current_indx.data();
  cudaMemcpyToSymbol(gpu::nu_bank_current_indx, &tmp, sizeof(void*));
  tmp = nu_bank.data();
  cudaMemcpyToSymbol(gpu::nu_bank, &tmp, sizeof(void*));
  tmp = keff_tally_absorption.data();
  cudaMemcpyToSymbol(gpu::keff_tally_absorption, &tmp, sizeof(void*));
  tmp = keff_tally_collision.data();
  cudaMemcpyToSymbol(gpu::keff_tally_collision, &tmp, sizeof(void*));
  tmp = keff_tally_tracklength.data();
  cudaMemcpyToSymbol(gpu::keff_tally_tracklength, &tmp, sizeof(void*));
  tmp = keff_tally_leakage.data();
  cudaMemcpyToSymbol(gpu::keff_tally_leakage, &tmp, sizeof(void*));
  tmp = trace.data();
  cudaMemcpyToSymbol(gpu::trace, &tmp, sizeof(void*));
  tmp = collision_distance.data();
  cudaMemcpyToSymbol(gpu::collision_distance, &tmp, sizeof(void*));
  tmp = n_event.data();
  cudaMemcpyToSymbol(gpu::n_event, &tmp, sizeof(void*));
  tmp = n_progeny.data();
  cudaMemcpyToSymbol(gpu::n_progeny, &tmp, sizeof(void*));
  tmp = neutron_xs.data();
  cudaMemcpyToSymbol(gpu::neutron_xs, &tmp, sizeof(void*));
  tmp = photon_xs.data();
  cudaMemcpyToSymbol(gpu::photon_xs, &tmp, sizeof(void*));
  tmp = macro_xs.data();
  cudaMemcpyToSymbol(gpu::macro_xs, &tmp, sizeof(void*));
  tmp = id.data();
  cudaMemcpyToSymbol(gpu::id, &tmp, sizeof(void*));
  tmp = type.data();
  cudaMemcpyToSymbol(gpu::type, &tmp, sizeof(void*));
  tmp = n_coord.data();
  cudaMemcpyToSymbol(gpu::n_coord, &tmp, sizeof(void*));
  tmp = cell_instance.data();
  cudaMemcpyToSymbol(gpu::cell_instance, &tmp, sizeof(void*));
  tmp = coord.data();
  cudaMemcpyToSymbol(gpu::coord, &tmp, sizeof(void*));
  tmp = n_coord_last.data();
  cudaMemcpyToSymbol(gpu::n_coord_last, &tmp, sizeof(void*));
  tmp = cell_last.data();
  cudaMemcpyToSymbol(gpu::cell_last, &tmp, sizeof(void*));
  tmp = E.data();
  cudaMemcpyToSymbol(gpu::E, &tmp, sizeof(void*));
  tmp = E_last.data();
  cudaMemcpyToSymbol(gpu::E_last, &tmp, sizeof(void*));
  tmp = g.data();
  cudaMemcpyToSymbol(gpu::g, &tmp, sizeof(void*));
  tmp = g_last.data();
  cudaMemcpyToSymbol(gpu::g_last, &tmp, sizeof(void*));
  tmp = wgt.data();
  cudaMemcpyToSymbol(gpu::wgt, &tmp, sizeof(void*));
  tmp = mu.data();
  cudaMemcpyToSymbol(gpu::mu, &tmp, sizeof(void*));
  tmp = alive.data();
  cudaMemcpyToSymbol(gpu::alive, &tmp, sizeof(void*));
  tmp = r_last_current.data();
  cudaMemcpyToSymbol(gpu::r_last_current, &tmp, sizeof(void*));
  tmp = r_last.data();
  cudaMemcpyToSymbol(gpu::r_last, &tmp, sizeof(void*));
  tmp = u_last.data();
  cudaMemcpyToSymbol(gpu::u_last, &tmp, sizeof(void*));
  tmp = wgt_last.data();
  cudaMemcpyToSymbol(gpu::wgt_last, &tmp, sizeof(void*));
  tmp = wgt_absorb.data();
  cudaMemcpyToSymbol(gpu::wgt_absorb, &tmp, sizeof(void*));
  tmp = fission.data();
  cudaMemcpyToSymbol(gpu::fission, &tmp, sizeof(void*));
  tmp = event.data();
  cudaMemcpyToSymbol(gpu::event, &tmp, sizeof(void*));
  tmp = event_nuclide.data();
  cudaMemcpyToSymbol(gpu::event_nuclide, &tmp, sizeof(void*));
  tmp = event_mt.data();
  cudaMemcpyToSymbol(gpu::event_mt, &tmp, sizeof(void*));
  tmp = delayed_group.data();
  cudaMemcpyToSymbol(gpu::delayed_group, &tmp, sizeof(void*));
  tmp = n_bank.data();
  cudaMemcpyToSymbol(gpu::n_bank, &tmp, sizeof(void*));
  tmp = n_bank_second.data();
  cudaMemcpyToSymbol(gpu::n_bank_second, &tmp, sizeof(void*));
  tmp = wgt_bank.data();
  cudaMemcpyToSymbol(gpu::wgt_bank, &tmp, sizeof(void*));
  tmp = n_delayed_bank.data();
  cudaMemcpyToSymbol(gpu::n_delayed_bank, &tmp, sizeof(void*));
  tmp = surface.data();
  cudaMemcpyToSymbol(gpu::surface, &tmp, sizeof(void*));
  tmp = cell_born.data();
  cudaMemcpyToSymbol(gpu::cell_born, &tmp, sizeof(void*));
  tmp = material.data();
  cudaMemcpyToSymbol(gpu::material, &tmp, sizeof(void*));
  tmp = material_last.data();
  cudaMemcpyToSymbol(gpu::material_last, &tmp, sizeof(void*));
  tmp = boundary.data();
  cudaMemcpyToSymbol(gpu::boundary, &tmp, sizeof(void*));
  tmp = sqrtkT.data();
  cudaMemcpyToSymbol(gpu::sqrtkT, &tmp, sizeof(void*));
  tmp = sqrtkT_last.data();
  cudaMemcpyToSymbol(gpu::sqrtkT_last, &tmp, sizeof(void*));
  tmp = n_collision.data();
  cudaMemcpyToSymbol(gpu::n_collision, &tmp, sizeof(void*));
  tmp = write_track.data();
  cudaMemcpyToSymbol(gpu::write_track, &tmp, sizeof(void*));
  tmp = seeds.data();
  cudaMemcpyToSymbol(gpu::seeds, &tmp, sizeof(void*));
  tmp = stream.data();
  cudaMemcpyToSymbol(gpu::stream, &tmp, sizeof(void*));
  tmp = secondary_bank_current_indx.data();
  cudaMemcpyToSymbol(gpu::secondary_bank_current_indx, &tmp, sizeof(void*));
  tmp = secondary_bank.data();
  cudaMemcpyToSymbol(gpu::secondary_bank, &tmp, sizeof(void*));
  tmp = current_work.data();
  cudaMemcpyToSymbol(gpu::current_work, &tmp, sizeof(void*));
  tmp = flux_derivs.data();
  cudaMemcpyToSymbol(gpu::flux_derivs, &tmp, sizeof(void*));
  tmp = filter_matches.data();
  cudaMemcpyToSymbol(gpu::filter_matches, &tmp, sizeof(void*));
  tmp = nu_bank_current_indx.data();
  cudaMemcpyToSymbol(gpu::nu_bank_current_indx, &tmp, sizeof(void*));
  tmp = nu_bank.data();
  cudaMemcpyToSymbol(gpu::nu_bank, &tmp, sizeof(void*));
  tmp = keff_tally_absorption.data();
  cudaMemcpyToSymbol(gpu::keff_tally_absorption, &tmp, sizeof(void*));
  tmp = keff_tally_collision.data();
  cudaMemcpyToSymbol(gpu::keff_tally_collision, &tmp, sizeof(void*));
  tmp = keff_tally_tracklength.data();
  cudaMemcpyToSymbol(gpu::keff_tally_tracklength, &tmp, sizeof(void*));
  tmp = keff_tally_leakage.data();
  cudaMemcpyToSymbol(gpu::keff_tally_leakage, &tmp, sizeof(void*));
  tmp = trace.data();
  cudaMemcpyToSymbol(gpu::trace, &tmp, sizeof(void*));
  tmp = collision_distance.data();
  cudaMemcpyToSymbol(gpu::collision_distance, &tmp, sizeof(void*));
  tmp = n_event.data();
  cudaMemcpyToSymbol(gpu::n_event, &tmp, sizeof(void*));
  tmp = n_progeny.data();
  cudaMemcpyToSymbol(gpu::n_progeny, &tmp, sizeof(void*));
#endif
}

namespace soa {
int n_nuclides;
int n_elements;
int n_coord_levels;
int n_tally_derivs;
int n_tally_filters;

namespace gpu {
__constant__ int n_nuclides;
__constant__ int n_elements;
__constant__ int n_coord_levels;
__constant__ int n_tally_derivs;
__constant__ int n_tally_filters;
} // namespace gpu

vector<NuclideMicroXS> neutron_xs;
vector<ElementMicroXS> photon_xs;
vector<MacroXS> macro_xs;
vector<int64_t> id;
vector<ParticleType> type;
vector<int> n_coord;
vector<int> cell_instance;
vector<LocalCoord> coord;
vector<int> n_coord_last;
vector<int> cell_last;
vector<double> E;
vector<double> E_last;
vector<int> g;
vector<int> g_last;
vector<double> wgt;
vector<double> mu;
vector<char> alive;
vector<Position> r_last_current;
vector<Position> r_last;
vector<Direction> u_last;
vector<double> wgt_last;
vector<double> wgt_absorb;
vector<char> fission;
vector<TallyEvent> event;
vector<int> event_nuclide;
vector<int> event_mt;
vector<int> delayed_group;
vector<int> n_bank;
vector<int> n_bank_second;
vector<double> wgt_bank;
vector<int> n_delayed_bank; // MAX_DELAYED_GROUPS pitch
vector<int> surface;
vector<int> cell_born;
vector<int> material;
vector<int> material_last;
vector<BoundaryInfo> boundary;
vector<double> sqrtkT;
vector<double> sqrtkT_last;
vector<int> n_collision;
vector<char> write_track;
vector<uint64_t> seeds; // N_STREAMS pitch
vector<int> stream;
vector<ParticleBank> secondary_bank; // max_secondary_particles pitch
vector<int> secondary_bank_current_indx;
vector<int64_t> current_work;
vector<double> flux_derivs;         // tally_derivs.size() pitch
vector<FilterMatch> filter_matches; // tally_filters.size() pitch
vector<NuBank> nu_bank;             // max_fissions pitch
vector<int> nu_bank_current_indx;
vector<double> keff_tally_absorption;
vector<double> keff_tally_collision;
vector<double> keff_tally_tracklength;
vector<double> keff_tally_leakage;
vector<char> trace;
vector<double> collision_distance;
vector<int> n_event;
vector<int64_t> n_progeny;

namespace gpu {
__constant__ NuclideMicroXS* neutron_xs;
__constant__ ElementMicroXS* photon_xs;
__constant__ MacroXS* macro_xs;
__constant__ int64_t* id;
__constant__ ParticleType* type;
__constant__ int* n_coord;
__constant__ int* cell_instance;
__constant__ LocalCoord* coord;
__constant__ int* n_coord_last;
__constant__ int* cell_last;
__constant__ double* E;
__constant__ double* E_last;
__constant__ int* g;
__constant__ int* g_last;
__constant__ double* wgt;
__constant__ double* mu;
__constant__ char* alive; // vector<bool* can't return references
__constant__ Position* r_last_current;
__constant__ Position* r_last;
__constant__ Direction* u_last;
__constant__ double* wgt_last;
__constant__ double* wgt_absorb;
__constant__ char* fission;
__constant__ TallyEvent* event;
__constant__ int* event_nuclide;
__constant__ int* event_mt;
__constant__ int* delayed_group;
__constant__ int* n_bank;
__constant__ int* n_bank_second;
__constant__ double* wgt_bank;
__constant__ int* n_delayed_bank; // MAX_DELAYED_GROUPS pitch
__constant__ int* surface;
__constant__ int* cell_born;
__constant__ int* material;
__constant__ int* material_last;
__constant__ BoundaryInfo* boundary;
__constant__ double* sqrtkT;
__constant__ double* sqrtkT_last;
__constant__ int* n_collision;
__constant__ char* write_track;
__constant__ uint64_t* seeds; // N_STREAMS pitch
__constant__ int* stream;
__constant__ int* secondary_bank_current_indx;
__constant__ ParticleBank* secondary_bank; // max_secondary_particles pitch
__constant__ int64_t* current_work;
__constant__ double* flux_derivs;         // tally_derivs.size() pitch
__constant__ FilterMatch* filter_matches; // tally_filters.size() pitch
__constant__ int* nu_bank_current_indx;
__constant__ NuBank* nu_bank;
__constant__ double* keff_tally_absorption;
__constant__ double* keff_tally_collision;
__constant__ double* keff_tally_tracklength;
__constant__ double* keff_tally_leakage;
__constant__ char* trace;
__constant__ double* collision_distance;
__constant__ int* n_event;
__constant__ int64_t* n_progeny;
__constant__ NuclideMicroXS* neutron_xs;
__constant__ ElementMicroXS* photon_xs;
__constant__ MacroXS* macro_xs;
__constant__ int64_t* id;
__constant__ ParticleType* type;
__constant__ int* n_coord;
__constant__ int* cell_instance;
__constant__ LocalCoord* coord;
__constant__ int* n_coord_last;
__constant__ int* cell_last;
__constant__ double* E;
__constant__ double* E_last;
__constant__ int* g;
__constant__ int* g_last;
__constant__ double* wgt;
__constant__ double* mu;
__constant__ char* alive;
__constant__ Position* r_last_current;
__constant__ Position* r_last;
__constant__ Direction* u_last;
__constant__ double* wgt_last;
__constant__ double* wgt_absorb;
__constant__ char* fission;
__constant__ TallyEvent* event;
__constant__ int* event_nuclide;
__constant__ int* event_mt;
__constant__ int* delayed_group;
__constant__ int* n_bank;
__constant__ int* n_bank_second;
__constant__ double* wgt_bank;
__constant__ int* n_delayed_bank;
__constant__ int* surface;
__constant__ int* cell_born;
__constant__ int* material;
__constant__ int* material_last;
__constant__ BoundaryInfo* boundary;
__constant__ double* sqrtkT;
__constant__ double* sqrtkT_last;
__constant__ int* n_collision;
__constant__ char* write_track;
__constant__ uint64_t* seeds;
__constant__ int* stream;
__constant__ int* secondary_bank_current_indx;
__constant__ ParticleBank* secondary_bank;
__constant__ int64_t* current_work;
__constant__ double* flux_derivs;
__constant__ FilterMatch* filter_matches;
__constant__ int* nu_bank_current_indx;
__constant__ NuBank* nu_bank;
__constant__ double* keff_tally_absorption;
__constant__ double* keff_tally_collision;
__constant__ double* keff_tally_tracklength;
__constant__ double* keff_tally_leakage;
__constant__ char* trace;
__constant__ double* collision_distance;
__constant__ int* n_event;
__constant__ int64_t* n_progeny;
} // namespace gpu
} // namespace soa

} // namespace openmc
