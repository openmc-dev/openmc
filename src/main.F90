program main

  use constants,    only: MODE_CRITICALITY, MODE_PLOTTING, MODE_FIXEDSOURCE
  use criticality,  only: run_criticality
  use finalize,     only: finalize_run
  use fixed_source, only: run_fixedsource
  use global,       only: run_mode
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
  end select
     
  ! finalize run
  call finalize_run()
  
end program main
