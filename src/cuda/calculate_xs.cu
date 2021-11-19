#include "openmc/cuda/calculate_xs.h"
#include "openmc/geometry.h"         // find_cell
#include "openmc/reaction_product.h" // EmissionMode
#include "openmc/search.h"
#include "openmc/settings.h" // BLOCKSIZE

namespace openmc {
namespace gpu {

__constant__ unique_ptr<Material>* materials;
__constant__ unique_ptr<Nuclide>* nuclides;
__constant__ Particle* particles;
__constant__ xsfloat energy_min_neutron;
__constant__ xsfloat energy_max_neutron;
__constant__ xsfloat log_spacing;
__constant__ unsigned number_nuclides;
__constant__ bool need_depletion_rx;

__managed__ unsigned managed_calculate_fuel_queue_index;
__managed__ unsigned managed_calculate_nonfuel_queue_index;

__global__ void __launch_bounds__(BLOCKSIZE) process_calculate_xs_events_device(
  EventQueueItem* __restrict__ queue, unsigned queue_size)
{
  using EmissionMode = ReactionProduct::EmissionMode;

  unsigned tid = threadIdx.x + blockDim.x * blockIdx.x;
  if (tid >= queue_size)
    return;
  Particle p(queue[tid].idx);
  auto const E = __ldg(&queue[tid].E);
  auto const mat_idx = __ldg(&queue[tid].material);

  // Store pre-collision particle properties
  p.wgt_last() = p.wgt();
  p.E_last() = E;
  p.u_last() = p.u();
  p.r_last() = p.r();

  // Reset event variables
  p.event() = TallyEvent::KILL;
  p.event_nuclide() = NUCLIDE_NONE;
  p.event_mt() = REACTION_NONE;

  p.macro_xs().total = 0.0;
  p.macro_xs().neutron.absorption = 0.0;
  p.macro_xs().neutron.fission = 0.0;
  p.macro_xs().neutron.nu_fission = 0.0;

  // Skip void material
  if (mat_idx == -1)
    return;

  Material const& m = *materials[mat_idx];

  unsigned i_log_union = std::log(E / energy_min_neutron) / log_spacing;

  // Add contribution from each nuclide in material
  auto const n_nuclides = m.nuclide_.size();
  for (int i = 0; i < n_nuclides; ++i) {
    auto const& i_nuclide = m.nuclide_[i];
    auto* __restrict__ micro {&p.neutron_xs(i_nuclide)};

    if (E != micro->last_E || p.sqrtkT() != micro->last_sqrtkT) {
      auto const& nuclide = *nuclides[i_nuclide];
      micro->elastic = CACHE_INVALID;
      micro->thermal = 0.0;
      micro->thermal_elastic = 0.0;
      micro->use_ptable = false;

      // Find the appropriate temperature index. why would someone use
      // nearest?
      xsfloat kT = p.sqrtkT() * p.sqrtkT();
      xsfloat f;
      int i_temp;

      switch (gpu::temperature_method) {
      case TemperatureMethod::NEAREST: {
        double max_diff = INFTY;
        for (int t = 0; t < nuclide.kTs_.size(); ++t) {
          double diff = std::abs(nuclide.kTs_[t] - kT);
          if (diff < max_diff) {
            i_temp = t;
            max_diff = diff;
          }
        }
      } break;

      case TemperatureMethod::INTERPOLATION:
        // Find temperatures that bound the actual temperature
        for (i_temp = 0; i_temp < nuclide.kTs_.size() - 1; ++i_temp) {
          if (nuclide.kTs_[i_temp] <= kT && kT < nuclide.kTs_[i_temp + 1])
            break;
        }

        // Randomly sample between temperature i and i+1
        f = (kT - nuclide.kTs_[i_temp]) /
            (nuclide.kTs_[i_temp + 1] - nuclide.kTs_[i_temp]);
        if (f > prn(p.current_seed()))
          ++i_temp;
        break;
      }

      const auto& grid {nuclide.grid_[i_temp]};
      int i_grid;
      if (E < grid.energy.front()) {
        i_grid = 0;
      } else if (E > grid.energy.back()) {
        i_grid = grid.energy.size() - 2;
      } else {
        // Determine bounding indices based on which equal log-spaced
        // interval the energy is in
        int i_low = __ldg(&grid.grid_index[i_log_union]);
        int i_high = __ldg(&grid.grid_index[i_log_union + 1]) + 1;

        // Perform binary search over reduced range
        i_grid = i_low + lower_bound_index_linear(
                           &grid.energy[i_low], &grid.energy[i_high], E);
      }
      const auto xs_left {nuclide.xs_[i_temp][i_grid]};
      const auto xs_right {nuclide.xs_[i_temp][i_grid + 1]};
      // check for rare case where two energy points are the same
      if (grid.energy[i_grid] == grid.energy[i_grid + 1])
        ++i_grid;

      // calculate interpolation factor
      f = (E - grid.energy[i_grid]) /
          (grid.energy[i_grid + 1] - grid.energy[i_grid]);

      micro->index_temp = i_temp;
      micro->index_grid = i_grid;
      micro->interp_factor = f;

      // Calculate all microscopic cross sections
      micro->total = (1.0 - f) * xs_left.total + f * xs_right.total;
      micro->absorption =
        (1.0 - f) * xs_left.absorption + f * xs_right.absorption;

      if (nuclide.fissionable_) {
        // Calculate microscopic nuclide total cross section
        micro->fission = (1.0 - f) * xs_left.fission + f * xs_right.fission;

        // Calculate microscopic nuclide nu-fission cross section
        micro->nu_fission =
          (1.0 - f) * xs_left.nu_fission + f * xs_right.nu_fission;
      } else {
        micro->fission = 0.0;
        micro->nu_fission = 0.0;
      }

      // Calculate microscopic nuclide photon production cross section
      micro->photon_prod =
        (1.0 - f) * xs_left.photon_production + f * xs_right.photon_production;

      micro->index_sab = C_NONE;
      micro->sab_frac = 0.0;
      micro->last_E = E;
      micro->last_sqrtkT = p.sqrtkT();

      // Calculate URR cross sections if needed
      if (gpu::urr_ptables_on && nuclide.urr_present_) {
        if (nuclide.urr_data_[i_temp].energy_in_bounds(E)) {
          micro->use_ptable = true;
          // TODO check storing by value
          const auto& urr = nuclide.urr_data_[i_temp];

          int i_energy = 0;
          while (E >= urr.energy_[i_energy + 1]) {
            ++i_energy;
          };

          p.stream() = STREAM_URR_PTABLE;
          double r =
            future_prn(static_cast<int64_t>(nuclide.index_), *p.current_seed());
          p.stream() = STREAM_TRACKING;

          int i_low = 0;
          while (urr.cdf_values_(i_energy, i_low) <= r) {
            ++i_low;
          };

          int i_up = 0;
          while (urr.cdf_values_(i_energy + 1, i_up) <= r) {
            ++i_up;
          };

          // Determine elastic, fission, and capture cross sections from the
          // probability table
          xsfloat elastic = 0.;
          xsfloat fission = 0.;
          xsfloat capture = 0.;
          xsfloat f;
          if (urr.interp_ == Interpolation::lin_lin) {
            // Determine the interpolation factor on the table
            f = (E - urr.energy_[i_energy]) /
                (urr.energy_[i_energy + 1] - urr.energy_[i_energy]);

            elastic = (1. - f) * urr.xs_values_(i_energy, i_low).elastic +
                      f * urr.xs_values_(i_energy + 1, i_up).elastic;
            fission = (1. - f) * urr.xs_values_(i_energy, i_low).fission +
                      f * urr.xs_values_(i_energy + 1, i_up).fission;
            capture = (1. - f) * urr.xs_values_(i_energy, i_low).n_gamma +
                      f * urr.xs_values_(i_energy + 1, i_up).n_gamma;
          } else if (urr.interp_ == Interpolation::log_log) {
            // Determine interpolation factor on the table
            f = std::log(E / urr.energy_[i_energy]) /
                std::log(urr.energy_[i_energy + 1] / urr.energy_[i_energy]);

            // Calculate the elastic cross section/factor
            if ((urr.xs_values_(i_energy, i_low).elastic > 0.) &&
                (urr.xs_values_(i_energy + 1, i_up).elastic > 0.)) {
              elastic = std::exp(
                (1. - f) * std::log(urr.xs_values_(i_energy, i_low).elastic) +
                f * std::log(urr.xs_values_(i_energy + 1, i_up).elastic));
            } else {
              elastic = 0.;
            }

            // Calculate the fission cross section/factor
            if ((urr.xs_values_(i_energy, i_low).fission > 0.) &&
                (urr.xs_values_(i_energy + 1, i_up).fission > 0.)) {
              fission = std::exp(
                (1. - f) * std::log(urr.xs_values_(i_energy, i_low).fission) +
                f * std::log(urr.xs_values_(i_energy + 1, i_up).fission));
            } else {
              fission = 0.;
            }

            // Calculate the capture cross section/factor
            if ((urr.xs_values_(i_energy, i_low).n_gamma > 0.) &&
                (urr.xs_values_(i_energy + 1, i_up).n_gamma > 0.)) {
              capture = std::exp(
                (1. - f) * std::log(urr.xs_values_(i_energy, i_low).n_gamma) +
                f * std::log(urr.xs_values_(i_energy + 1, i_up).n_gamma));
            } else {
              capture = 0.;
            }
          }

          // Determine the treatment of inelastic scattering
          xsfloat inelastic = 0.;
          if (urr.inelastic_flag_ != C_NONE) {
            // get interpolation factor
            f = micro->interp_factor;

            // Determine inelastic scattering cross section
            Reaction* rx = nuclide.reactions_[nuclide.urr_inelastic_].get();
            int xs_index = micro->index_grid - rx->xs_[i_temp].threshold;
            if (xs_index >= 0) {
              inelastic = (1. - f) * rx->xs_[i_temp].value[xs_index] +
                          f * rx->xs_[i_temp].value[xs_index + 1];
            }
          }

          // Multiply by smooth cross-section if needed
          if (urr.multiply_smooth_) {
            const auto& xs = nuclide.reactions_[0]->xs_[i_temp].value;
            f = micro->interp_factor;
            micro->elastic = (1.0 - f) * xs[i_grid] + f * xs[i_grid + 1];
            elastic *= micro->elastic;
            capture *= (micro->absorption - micro->fission);
            fission *= micro->fission;
          }

          // Check for negative values
          if (elastic < 0.) {
            elastic = 0.;
          }
          if (fission < 0.) {
            fission = 0.;
          }
          if (capture < 0.) {
            capture = 0.;
          }

          // Set elastic, absorption, fission, total, and capture x/s. Note that
          // the total x/s is calculated as a sum of partials instead of the
          // table-provided value
          micro->elastic = elastic;
          micro->absorption = capture + fission;
          micro->fission = fission;
          micro->total = elastic + inelastic + capture + fission;

          // if (simulation::need_depletion_rx) {
          //   micro.reaction[0] = capture;
          // }

          // Determine nu-fission cross-section
          if (nuclide.fissionable_) {
            micro->nu_fission =
              nuclide.nu(E, EmissionMode::total) * micro->fission;
          }
        }
      }
    }

    double const& atom_density = m.atom_density_[i];
    p.macro_xs().total += atom_density * micro->total;
    p.macro_xs().neutron.absorption += atom_density * micro->absorption;
    p.macro_xs().neutron.fission += atom_density * micro->fission;
    p.macro_xs().neutron.nu_fission += atom_density * micro->nu_fission;
  }
}

} // namespace gpu
} // namespace openmc
