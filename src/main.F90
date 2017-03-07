program main

  use constants
  use eigenvalue,        only: run_eigenvalue
  use finalize,          only: finalize_run
  use fixed_source,      only: run_fixedsource
  use global
  use initialize,        only: initialize_run
  use particle_restart,  only: run_particle_restart
  use plot,              only: run_plot

  implicit none

  ! set up problem
  call initialize_run()

  ! start problem based on mode
  select case (run_mode)
  case (MODE_FIXEDSOURCE)
    call run_fixedsource()
  case (MODE_EIGENVALUE)
    call run_eigenvalue()
  case (MODE_PLOTTING)
    call run_plot()
  case (MODE_PARTICLE)
    if (master) call run_particle_restart()
  end select

  ! finalize run
  call finalize_run()

end program main
