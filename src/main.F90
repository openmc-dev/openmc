program main

  use constants
  use message_passing
  use openmc_api, only: openmc_init, openmc_finalize, openmc_run, &
                        openmc_plot_geometry, openmc_calculate_volumes
  use particle_restart, only: run_particle_restart
  use settings, only: run_mode

  implicit none

#ifdef MPI
  integer :: mpi_err ! MPI error code
#endif

  ! Initialize run -- when run with MPI, pass communicator
#ifdef MPI
#ifdef MPIF08
  call openmc_init(MPI_COMM_WORLD % MPI_VAL)
#else
  call openmc_init(MPI_COMM_WORLD)
#endif
#else
  call openmc_init()
#endif

  ! start problem based on mode
  select case (run_mode)
  case (MODE_FIXEDSOURCE, MODE_EIGENVALUE)
    call openmc_run()
  case (MODE_PLOTTING)
    call openmc_plot_geometry()
  case (MODE_PARTICLE)
    if (master) call run_particle_restart()
  case (MODE_VOLUME)
    call openmc_calculate_volumes()
  end select

  ! finalize run
  call openmc_finalize()

#ifdef MPI
  ! If MPI is in use and enabled, terminate it
  call MPI_FINALIZE(mpi_err)
#endif


end program main
