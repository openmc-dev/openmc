module finalize

# ifdef PETSC
  use cmfd_output,    only: finalize_cmfd
# endif
  use global
  use output,         only: print_runtime
  use tally,          only: write_tallies, tally_statistics
  use timing,         only: timer_start, timer_stop

#ifdef MPI
  use mpi
#endif

#ifdef HDF5
  use hdf5_interface, only: hdf5_finalize
#endif

  implicit none

contains

!===============================================================================
! FINALIZE_RUN does all post-simulation tasks such as calculating tally
! statistics and writing out tallies
!===============================================================================

  subroutine finalize_run()

    ! Start finalization timer
    call timer_start(time_finalize)

    if (run_mode /= MODE_PLOTTING) then
       if (output_tallies) then
          ! Calculate statistics for tallies and write to tallies.out
          call tally_statistics()
          if (master) call write_tallies()
       end if
    end if

# ifdef PETSC
    ! finalize cmfd
    if (cmfd_run) call finalize_cmfd()
# endif

    ! stop timers and show timing statistics
    call timer_stop(time_finalize)
    call timer_stop(time_total)
    if (master .and. (run_mode /= MODE_PLOTTING .and. &
         run_mode /= MODE_TALLIES)) call print_runtime()

    ! deallocate arrays
    call free_memory()

#ifdef HDF5
    ! Close HDF5 interface and release memory
    call hdf5_finalize()
#endif
    
#ifdef MPI
    ! If MPI is in use and enabled, terminate it
    call MPI_FINALIZE(mpi_err)
#endif

  end subroutine finalize_run

end module finalize
