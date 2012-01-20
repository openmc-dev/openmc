module finalize

  use global
  use output,         only: print_runtime
  use tally,          only: write_tallies, tally_statistics
  use timing,         only: timer_stop

#ifdef HDF5
  use hdf5_interface, only: hdf5_write_timing, hdf5_close_output
#endif

  implicit none

contains

!===============================================================================
! FINALIZE_RUN
!===============================================================================

  subroutine finalize_run()

    ! Calculate statistics for tallies and write to tallies.out
    if (.not. plotting) then
       call tally_statistics()
       if (master) call write_tallies()
    end if

    ! show timing statistics
    call timer_stop(time_total)
    if (master .and. (.not. plotting)) call print_runtime()

#ifdef HDF5
    if (master) then
       call hdf5_write_timing()
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
