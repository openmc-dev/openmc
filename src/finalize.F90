module finalize

  use global
  use output,         only: print_runtime, print_results, &
                            print_overlap_check, write_tallies
  use tally,          only: tally_statistics

#ifdef MPI
  use message_passing
#endif

  use hdf5_interface, only: hdf5_bank_t
  use hdf5, only: h5tclose_f, h5close_f

  implicit none

contains

!===============================================================================
! FINALIZE_RUN does all post-simulation tasks such as calculating tally
! statistics and writing out tallies
!===============================================================================

  subroutine finalize_run()

    integer :: hdf5_err

    ! Start finalization timer
    call time_finalize%start()

    if (run_mode /= MODE_PLOTTING .and. run_mode /= MODE_PARTICLE) then
      ! Calculate statistics for tallies and write to tallies.out
      if (master) then
        if (n_realizations > 1) call tally_statistics()
      end if
      if (output_tallies) then
        if (master) call write_tallies()
      end if
      if (check_overlaps) call reduce_overlap_count()
    end if

    ! Stop timers and show timing statistics
    call time_finalize%stop()
    call time_total%stop()
    if (master .and. (run_mode /= MODE_PLOTTING .and. &
         run_mode /= MODE_PARTICLE)) then
      call print_runtime()
      call print_results()
      if (check_overlaps) call print_overlap_check()
    end if

    ! Deallocate arrays
    call free_memory()

    ! Release compound datatypes
    call h5tclose_f(hdf5_bank_t, hdf5_err)

    ! Close FORTRAN interface.
    call h5close_f(hdf5_err)

#ifdef MPI
    ! Free all MPI types
    call MPI_TYPE_FREE(MPI_BANK, mpi_err)

    ! If MPI is in use and enabled, terminate it
    call MPI_FINALIZE(mpi_err)
#endif

  end subroutine finalize_run

!===============================================================================
! REDUCE_OVERLAP_COUNT accumulates cell overlap check counts to master
!===============================================================================

  subroutine reduce_overlap_count()

#ifdef MPI
      if (master) then
        call MPI_REDUCE(MPI_IN_PLACE, overlap_check_cnt, n_cells, &
             MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
      else
        call MPI_REDUCE(overlap_check_cnt, overlap_check_cnt, n_cells, &
             MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
      end if
#endif

  end subroutine reduce_overlap_count

end module finalize
