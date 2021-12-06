# This Makefile will build a statically linked OpenMC executable
# in CUDA mode. It seems that CMake is doing something which screws
# up this process (either has to do with creating libopenmc.so, or
# that it's not using nvlink as the linker and I'm not sure how to
# make it do this). Moreover, it's impossible to have CMake use -dlto
# along with its method of specifying GPU architectures. This
# file sidesteps all that.
#
# Notably, because this doesn't create a functioning libopenmc.so,
# you can't call C API functions through the python interface if
# you build this way. If you really need that, you'll have to build
# a non-LTO GPU OpenMC.

NVCC=/opt/cuda/bin/nvcc
ARCH=sm_61
COMPILE_OPTIONS=-forward-unknown-to-host-compiler -DGIT_SHA1=\"definitelyIncorrectSHA\" -DGSL_THROW_ON_CONTRACT_VIOLATION -DHD="__host__ __device__" -Dgsl_CONFIG_ALLOWS_NONSTRICT_SPAN_COMPARISON=1 -Dgsl_CONFIG_DEFAULTS_VERSION=1 -Dlibopenmc_EXPORTS -I/home/gavin/Code/openmc/include -I/home/gavin/Code/openmc/vendor/pugixml/src -I/home/gavin/Code/openmc/vendor/faddeeva -I/home/gavin/Code/openmc/vendor/xtensor/include -I/home/gavin/Code/openmc/vendor/xtl/include -I/home/gavin/Code/openmc/vendor/gsl-lite/include -I/home/gavin/Code/openmc/vendor/fmt/include -arch=$(ARCH) -dlto -DH5Oget_info_by_idx_vers=1 -DH5O_info_t_vers=1 -fopenmp -O2 -fno-omit-frame-pointer -Xptxas=-dlcm=ca --expt-relaxed-constexpr -std=c++14 -x cu -Dxsfloat=double -DHOST=__host__ -DHD="__host__ __device__" -DMANAGED=__managed__
SOURCEFILES=bank.o \
  boundary_condition.o \
  bremsstrahlung.o \
  cell.o \
  cmfd_solver.o \
  cross_sections.o \
  dagmc.o \
  distribution.o \
  distribution_angle.o \
  distribution_energy.o \
  distribution_multi.o \
  distribution_spatial.o \
  eigenvalue.o \
  endf.o \
  error.o \
  event.o \
  finalize.o \
  geometry.o \
  geometry_aux.o \
  hdf5_interface.o \
  initialize.o \
  lattice.o \
  material.o \
  math_functions.o \
  mesh.o \
  message_passing.o \
  mgxs.o \
  mgxs_interface.o \
  nuclide.o \
  output.o \
  particle.o \
  particle_data.o \
  particle_restart.o \
  photon.o \
  physics.o \
  physics_common.o \
  physics_mg.o \
  plot.o \
  position.o \
  progress_bar.o \
  random_lcg.o \
  random_dist.o \
  reaction.o \
  reaction_product.o \
  scattdata.o \
  secondary_correlated.o \
  secondary_kalbach.o \
  secondary_nbody.o \
  secondary_thermal.o \
  secondary_uncorrelated.o \
  settings.o \
  simulation.o \
  soa_particle.o \
  source.o \
  state_point.o \
  string_utils.o \
  summary.o \
  surface.o \
  thermal.o \
  timer.o \
  track_output.o \
  urr.o \
  volume_calc.o \
  wmp.o \
  xml_interface.o \
  xsdata.o \
  tallies/derivative.o \
  tallies/filter.o \
  tallies/filter_azimuthal.o \
  tallies/filter_cell.o \
  tallies/filter_cell_instance.o \
  tallies/filter_cellborn.o \
  tallies/filter_cellfrom.o \
  tallies/filter_collision.o \
  tallies/filter_delayedgroup.o \
  tallies/filter_distribcell.o \
  tallies/filter_energy.o \
  tallies/filter_energyfunc.o \
  tallies/filter_legendre.o \
  tallies/filter_material.o \
  tallies/filter_mesh.o \
  tallies/filter_meshsurface.o \
  tallies/filter_mu.o \
  tallies/filter_particle.o \
  tallies/filter_polar.o \
  tallies/filter_sph_harm.o \
  tallies/filter_sptl_legendre.o \
  tallies/filter_surface.o \
  tallies/filter_universe.o \
  tallies/filter_zernike.o \
  tallies/tally.o \
  tallies/tally_scoring.o \
  tallies/trigger.o \
  cuda/advance.o \
  cuda/calculate_xs.o \
  cuda/refill.o \
  main.o

openmc: objectfiles
	nvcc -dlto -arch=$(ARCH) $(SOURCEFILES) ../build/lib/libfaddeeva.a ../build/lib/libfmt.a ../build/lib/libpugixml.a -lgomp -lhdf5_hl -lhdf5 -o openmc
.PHONY: objectfiles clean
objectfiles: $(SOURCEFILES)
%.o: %.cpp
	$(NVCC) $(COMPILE_OPTIONS) -dc $< -o $@
%.o: %.cu
	$(NVCC) $(COMPILE_OPTIONS) -c $< -o $@
clean:
	rm *.o 2> /dev/null
	rm tallies/*.o 2> /dev/null
	rm cuda/*.o 2> /dev/null
