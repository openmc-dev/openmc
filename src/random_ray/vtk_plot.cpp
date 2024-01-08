#include "openmc/random_ray/vtk_plot.h"
#include "openmc/mgxs_interface.h"
#include "openmc/output.h"
#include "openmc/plot.h"

namespace openmc {

// Outputs all basic material, FSR ID, multigroup flux, and
// fission source data to .vtk file that can be directly
// loaded and displayed by Paraview.
void plot_3D_vtk()
{
  // Print header information
  print_plot();

  // Outer loop over plots
  for (int p = 0; p < model::plots.size(); p++) {

    // Get handle to OpenMC plot object and extract params
    Plot* openmc_plot = dynamic_cast<Plot*>(model::plots[p].get());

    // Random ray plots only support voxel plots
    if (openmc_plot == nullptr) {
      warning(fmt::format("Plot {} is invalid plot type -- only voxel plotting "
                          "is allowed in random ray mode.",
        p));
      continue;
    } else if (openmc_plot->type_ != Plot::PlotType::voxel) {
      warning(fmt::format("Plot {} is invalid plot type -- only voxel plotting "
                          "is allowed in random ray mode.",
        p));
      continue;
    }

    int Nx = openmc_plot->pixels_[0];
    int Ny = openmc_plot->pixels_[1];
    int Nz = openmc_plot->pixels_[2];
    Position origin = openmc_plot->origin_;
    Position width = openmc_plot->width_;
    Position ll = origin - width / 2.0;
    double x_delta = width.x / Nx;
    double y_delta = width.y / Ny;
    double z_delta = width.z / Nz;
    std::string filename = openmc_plot->path_plot();
    filename.replace(filename.size() - 3, 3, ".vtk");

    int negroups = data::mg.num_energy_groups_;

    // Perform sanity checks on file size
    uint64_t bytes = Nx * Ny * Nz * (negroups + 1 + 1 + 1) * sizeof(float);
    write_message(5, "Processing plot {}: {}... (Estimated size is {} MB)",
      openmc_plot->id(), filename, bytes / 1.0e6);
    if (bytes / 1.0e9 > 1.0) {
      warning("Voxel plot specification is very large (>1 GB). Plotting may be "
              "slow.");
    } else if (bytes / 1.0e9 > 100.0) {
      fatal_error("Voxel plot specification is too large (>100 GB). Exiting.");
    }

    // Relate voxel spatial locations to random ray source regions
    std::vector<int> voxel_indices(Nx * Ny * Nz);

#pragma omp parallel for collapse(3)
    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
          Position sample;
          sample.z = ll.z + z_delta / 2.0 + z * z_delta;
          sample.y = ll.y + y_delta / 2.0 + y * y_delta;
          sample.x = ll.x + x_delta / 2.0 + x * x_delta;
          Particle p;
          p.r() = sample;
          bool found = exhaustive_find_cell(p);
          int i_cell = p.lowest_coord().cell;
          int64_t source_region_idx =
            random_ray::source_region_offsets[i_cell] + p.cell_instance();
          voxel_indices[z * Ny * Nx + y * Nx + x] = source_region_idx;
        }
      }
    }

    // Open file for writing
    FILE* plot = fopen(filename.c_str(), "w");

    // Write vtk metadata
    fprintf(plot, "# vtk DataFile Version 2.0\n");
    fprintf(plot, "Dataset File\n");
    fprintf(plot, "BINARY\n");
    fprintf(plot, "DATASET STRUCTURED_POINTS\n");
    fprintf(plot, "DIMENSIONS %d %d %d\n", Nx, Ny, Nz);
    fprintf(plot, "ORIGIN 0 0 0\n");
    fprintf(plot, "SPACING %lf %lf %lf\n", x_delta, y_delta, z_delta);
    fprintf(plot, "POINT_DATA %d\n", Nx * Ny * Nz);

    // Plot multigroup flux data
    for (int g = 0; g < negroups; g++) {
      fprintf(plot, "SCALARS flux_group_%d float\n", g);
      fprintf(plot, "LOOKUP_TABLE default\n");
      for (int fsr : voxel_indices) {
        int64_t source_element = fsr * negroups + g;
        float flux = random_ray::scalar_flux_final[source_element];
        flux /= (settings::n_batches - settings::n_inactive);
        flux = flip_endianness<float>(flux);
        fwrite(&flux, sizeof(float), 1, plot);
      }
    }

    // Plot FSRs
    fprintf(plot, "SCALARS FSRs float\n");
    fprintf(plot, "LOOKUP_TABLE default\n");
    for (int fsr : voxel_indices) {
      float value = future_prn(10, fsr);
      value = flip_endianness<float>(value);
      fwrite(&value, sizeof(float), 1, plot);
    }

    // Plot Materials
    fprintf(plot, "SCALARS Materials int\n");
    fprintf(plot, "LOOKUP_TABLE default\n");
    for (int fsr : voxel_indices) {
      int mat = random_ray::material[fsr];
      mat = flip_endianness<int>(mat);
      fwrite(&mat, sizeof(int), 1, plot);
    }

    // Plot fission source
    fprintf(plot, "SCALARS total_fission_source float\n");
    fprintf(plot, "LOOKUP_TABLE default\n");
    for (int fsr : voxel_indices) {
      float total_fission = 0.0;
      int mat = random_ray::material[fsr];
      for (int g = 0; g < negroups; g++) {
        int64_t source_element = fsr * negroups + g;
        float flux = random_ray::scalar_flux_final[source_element];
        flux /= (settings::n_batches - settings::n_inactive);
        float Sigma_f = data::mg.macro_xs_[mat].get_xs(
          MgxsType::FISSION, g, nullptr, nullptr, nullptr, 0, 0);
        total_fission += Sigma_f * flux;
      }
      total_fission = flip_endianness<float>(total_fission);
      fwrite(&total_fission, sizeof(float), 1, plot);
    }

    fclose(plot);
  }
}

} // namespace openmc
