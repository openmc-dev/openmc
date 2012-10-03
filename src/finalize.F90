module finalize

  use cmfd_output,    only: finalize_cmfd
  use global
  use output,         only: print_runtime
  use source,         only: write_source_binary
  use tally,          only: write_tallies, tally_statistics
  use timing,         only: timer_start, timer_stop

#ifdef HDF5
  use hdf5_interface, only: hdf5_write_results, hdf5_close_output
#endif

  implicit none

contains

!===============================================================================
! FINALIZE_RUN does all post-simulation tasks such as calculating tally
! statistics, writing out tallies, and writing hdf5 output.
!===============================================================================

  subroutine finalize_run()

    ! Start finalization timer
    call timer_start(time_finalize)

    if (run_mode /= MODE_PLOTTING) then
       ! Calculate statistics for tallies and write to tallies.out
       call tally_statistics()
       if (master) call write_tallies()

       ! Write out binary source
       if (write_source) call write_source_binary()
    end if

    ! finalize cmfd
    if (cmfd_run) call finalize_cmfd()

    ! stop timers and show timing statistics
    call timer_stop(time_finalize)
    call timer_stop(time_total)
    if (master .and. (run_mode /= MODE_PLOTTING)) call print_runtime()

#ifdef HDF5
    ! Write time statistics to HDF5 output 
    if (master) then
       call hdf5_write_results()
       call hdf5_close_output()
    end if
#endif

    ! deallocate arrays
    call free_memory()
    
#ifdef MPI
    ! If MPI is in use and enabled, terminate it
    call MPI_FINALIZE(mpi_err)
#endif

  end subroutine finalize_run

end module finalize
