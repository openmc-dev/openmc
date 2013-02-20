module finalize

# ifdef PETSC
  use cmfd_output,    only: finalize_cmfd
# endif
  use global
  use output,         only: print_runtime, print_results, write_tallies
  use tally,          only: tally_statistics

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
    call time_finalize % start()

    if (run_mode /= MODE_PLOTTING) then
      ! Calculate statistics for tallies and write to tallies.out
      if (master) call tally_statistics()
      if (output_tallies) then
        if (master) call write_tallies()
      end if
    end if

#ifdef PETSC
    ! finalize cmfd
    if (cmfd_run) call finalize_cmfd()
#endif

    ! stop timers and show timing statistics
    call time_finalize % stop()
    call time_total % stop()
    if (master .and. (run_mode /= MODE_PLOTTING .and. &
         run_mode /= MODE_TALLIES)) then
      call print_runtime()
      call print_results()
    end if

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
