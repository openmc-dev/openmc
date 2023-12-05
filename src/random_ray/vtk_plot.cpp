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
	FILE * fast;
  char* fname = "plots.vtk";
	fast = fopen(fname, "w");
  
  Source* s = model::external_sources[0].get();

  // Use box source for plotting bounds
  IndependentSource* is = dynamic_cast<IndependentSource*>(s);
  SpatialDistribution* space_dist = is->space();
  SpatialBox* sb = dynamic_cast<SpatialBox*>(space_dist);
  Position ll = sb->lower_left();
  Position ur = sb->upper_right();

  int Nx = 1000;
  int Ny = 1000;
  int Nz = 1;

	double x_delta = (ur.x - ll.x) / Nx;
	double y_delta = (ur.y - ll.y) / Ny;
	double z_delta = (ur.z - ll.z) / Nz;

	fprintf(fast,"# vtk DataFile Version 2.0\n");
	fprintf(fast, "Dataset File\n");
	fprintf(fast, "BINARY\n");
	fprintf(fast, "DATASET STRUCTURED_POINTS\n");
	fprintf(fast, "DIMENSIONS %d %d %d\n", Nx, Ny, Nz);
	fprintf(fast, "ORIGIN 0 0 0\n");
	fprintf(fast, "SPACING %lf %lf %lf\n", x_delta, y_delta, z_delta);
	fprintf(fast, "POINT_DATA %d\n", Nx*Ny*Nz);
  
  Position sample;
  std::vector<int> voxel_indices;

  for( int z = 0; z < Nz; z++ )
  {
    sample.z = ll.z + z_delta/2.0 + z * z_delta;
    for( int y = 0; y < Ny; y++ )
    {
      sample.y = ll.y + y_delta/2.0 + y * y_delta;
      for( int x = 0; x < Nx; x++ )
      {
        sample.x = ll.x + x_delta/2.0 + x * x_delta;
        Particle p;
        p.r() = sample;  
        bool found = exhaustive_find_cell(p);
        int i_cell = p.lowest_coord().cell;
        int64_t source_region_idx = random_ray::source_region_offsets[i_cell] + p.cell_instance();
        voxel_indices.push_back(source_region_idx);
      }
    }
  }
  
  int negroups = data::mg.num_energy_groups_;
  
  // Plot multigroup flux data
  for( int g = 0; g < negroups; g++ ) {
    fprintf(fast, "SCALARS flux_group_%d float\n", g);
    fprintf(fast, "LOOKUP_TABLE default\n");
    for (int fsr : voxel_indices) {
      int64_t source_element = fsr * negroups + g;
      float flux = random_ray::scalar_flux_old[source_element];
      flux = eswap_float(flux);
      fwrite(&flux, sizeof(float), 1, fast);
    }
  }
  
  // Plot FSRs
  fprintf(fast, "SCALARS FSRs float\n");
  fprintf(fast, "LOOKUP_TABLE default\n");
  for (int fsr : voxel_indices) {
    float value = future_prn(10, fsr);
    value = eswap_float(value);
    fwrite(&value, sizeof(float), 1, fast);
  }

  // Plot Materials
  fprintf(fast, "SCALARS Materials int\n");
  fprintf(fast, "LOOKUP_TABLE default\n");
  for (int fsr : voxel_indices) {
    int mat = random_ray::material[fsr];
    mat = eswap_int(mat);
    fwrite(&mat, sizeof(int), 1, fast);
  }

	fclose(fast);
	printf("Finished plotting!\n");
}

} // namespace openmc
