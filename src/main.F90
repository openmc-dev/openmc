program main

  use constants
  use finalize,          only: openmc_finalize
  use global
  use initialize,        only: openmc_init
  use message_passing
  use particle_restart,  only: run_particle_restart
  use plot,              only: run_plot
  use simulation,        only: run_simulation
  use volume_calc,       only: run_volume_calculations

  implicit none

  ! Initialize run -- when run with MPI, pass communicator
#ifdef MPI
  call openmc_init(MPI_COMM_WORLD)
#else
  call openmc_init()
#endif

  ! start problem based on mode
  select case (run_mode)
  case (MODE_FIXEDSOURCE, MODE_EIGENVALUE)
    call run_simulation()
  case (MODE_PLOTTING)
    call run_plot()
  case (MODE_PARTICLE)
    if (master) call run_particle_restart()
  case (MODE_VOLUME)
    call run_volume_calculations()
  end select

  ! finalize run
  call openmc_finalize()

end program main
