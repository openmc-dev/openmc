#include <stdio.h>

#include "openmc/random_ray/iteration.h"
#include "openmc/random_ray/tally_convert.h"
#include "openmc/random_ray/vtk_plot.h"
#include "openmc/output.h"
#include "openmc/geometry.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally_scoring.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/eigenvalue.h"
#include "openmc/timer.h"
#include "openmc/mgxs_interface.h"
#include "openmc/message_passing.h"
#include "openmc/plot.h"

namespace openmc {

float eswap_float( float f )
{    
	char * ptr = (char *) &f;
	char orig[4];
	orig[0] = ptr[0];
	orig[1] = ptr[1];
	orig[2] = ptr[2];
	orig[3] = ptr[3];
	char swapper[4];
	swapper[0] = orig[3];
	swapper[1] = orig[2];
	swapper[2] = orig[1];
	swapper[3] = orig[0];
	float f2 = *((float *)swapper);
	return f2;
}

int eswap_int( int f )
{    
	char * ptr = (char *) &f;
	char orig[4];
	orig[0] = ptr[0];
	orig[1] = ptr[1];
	orig[2] = ptr[2];
	orig[3] = ptr[3];
	char swapper[4];
	swapper[0] = orig[3];
	swapper[1] = orig[2];
	swapper[2] = orig[1];
	swapper[3] = orig[0];
	int f2 = *((int *)swapper);
	return f2;
}

void plot_3D_vtk()
{
  //char* fname = "plots.vtk";
	//plot = fopen(fname, "w");
  
  // Get handle to plot
  Plot& openmc_plot = *dynamic_cast<Plot*>(model::plots[0].get());
  int Nx = openmc_plot.pixels_[0];
  int Ny = openmc_plot.pixels_[1];
  int Nz = openmc_plot.pixels_[2];
  Position origin = openmc_plot.origin_;           //!< Plot origin in geometry
  Position width = openmc_plot.width_;            //!< Plot width in geometry
  
  std::string filename = openmc_plot.path_plot();
  filename.replace(filename.size() - 3, 3, ".vtk");
	FILE* plot = fopen(filename.c_str(), "w");
  
  // Use box source for plotting bounds
  //Source* s = model::external_sources[0].get();
  //IndependentSource* is = dynamic_cast<IndependentSource*>(s);
  //SpatialDistribution* space_dist = is->space();
  //SpatialBox* sb = dynamic_cast<SpatialBox*>(space_dist);

  Position ll = origin - width/2.0;

	double x_delta = width.x / Nx;
	double y_delta = width.y / Ny;
	double z_delta = width.z / Nz;

	fprintf(plot,"# vtk DataFile Version 2.0\n");
	fprintf(plot, "Dataset File\n");
	fprintf(plot, "BINARY\n");
	fprintf(plot, "DATASET STRUCTURED_POINTS\n");
	fprintf(plot, "DIMENSIONS %d %d %d\n", Nx, Ny, Nz);
	fprintf(plot, "ORIGIN 0 0 0\n");
	fprintf(plot, "SPACING %lf %lf %lf\n", x_delta, y_delta, z_delta);
	fprintf(plot, "POINT_DATA %d\n", Nx*Ny*Nz);
  
  std::vector<int> voxel_indices(Nx*Ny*Nz);

  #pragma omp parallel for collapse(3)
  for (int z = 0; z < Nz; z++) {
    for (int y = 0; y < Ny; y++) {
      for (int x = 0; x < Nx; x++) {
        Position sample;
        sample.z = ll.z + z_delta/2.0 + z * z_delta;
        sample.y = ll.y + y_delta/2.0 + y * y_delta;
        sample.x = ll.x + x_delta/2.0 + x * x_delta;
        Particle p;
        p.r() = sample;  
        bool found = exhaustive_find_cell(p);
        int i_cell = p.lowest_coord().cell;
        int64_t source_region_idx = random_ray::source_region_offsets[i_cell] + p.cell_instance();
        voxel_indices[z*Ny*Nx + y*Nx + x] = source_region_idx;
      }
    }
  }
  
  int negroups = data::mg.num_energy_groups_;
  
  // Plot multigroup flux data
  for( int g = 0; g < negroups; g++ ) {
    fprintf(plot, "SCALARS flux_group_%d float\n", g);
    fprintf(plot, "LOOKUP_TABLE default\n");
    for (int fsr : voxel_indices) {
      int64_t source_element = fsr * negroups + g;
      float flux = random_ray::scalar_flux_old[source_element];
      flux = eswap_float(flux);
      fwrite(&flux, sizeof(float), 1, plot);
    }
  }
  
  // Plot FSRs
  fprintf(plot, "SCALARS FSRs float\n");
  fprintf(plot, "LOOKUP_TABLE default\n");
  for (int fsr : voxel_indices) {
    float value = future_prn(10, fsr);
    value = eswap_float(value);
    fwrite(&value, sizeof(float), 1, plot);
  }

  // Plot Materials
  fprintf(plot, "SCALARS Materials int\n");
  fprintf(plot, "LOOKUP_TABLE default\n");
  for (int fsr : voxel_indices) {
    int mat = random_ray::material[fsr];
    mat = eswap_int(mat);
    fwrite(&mat, sizeof(int), 1, plot);
  }
  
  // Plot fission source
  fprintf(plot, "SCALARS total_fission_source float\n");
  fprintf(plot, "LOOKUP_TABLE default\n");
  for (int fsr : voxel_indices) {
    float total_fission = 0.0;
    int mat = random_ray::material[fsr];
    for( int g = 0; g < negroups; g++ ) {
      int64_t source_element = fsr * negroups + g;
      float flux = random_ray::scalar_flux_old[source_element];
      float Sigma_f = data::mg.macro_xs_[mat].get_xs(MgxsType::FISSION, g, nullptr, nullptr, nullptr, 0, 0);
      total_fission += Sigma_f * flux;
    }
    total_fission = eswap_float(total_fission);
    fwrite(&total_fission, sizeof(float), 1, plot);
  }

	fclose(plot);
	printf("Finished plotting!\n");
}

} // namespace openmc
