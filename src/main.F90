program main

  use constants
  use criticality,  only: run_criticality
  use finalize,     only: finalize_run
  use fixed_source, only: run_fixedsource
  use global
  use initialize,   only: initialize_run
  use plotter,      only: run_plot

  implicit none

  ! set up problem
  call initialize_run()

  ! start problem based on mode
  select case (run_mode)
  case (MODE_FIXEDSOURCE)
     call run_fixedsource()
  case (MODE_CRITICALITY)
     call run_criticality()
  case (MODE_PLOTTING)
     call run_plot()
  case (MODE_TALLIES)
     ! For tallies-only mode, we just skip straight to finalize_run to write out
     ! the tally results
     n_realizations = restart_batch - n_inactive
  end select
     
  ! finalize run
  call finalize_run()

end program main 
