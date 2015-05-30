module finalize

  use constants
  use global
  use output,         only: header, print_runtime, print_results, &
                            print_overlap_check, write_tallies, print_testing, &
                            print_domain_interactions
  use state_point,    only: write_distribmat_comps
  use tally,          only: tally_statistics

#ifdef MPI
  use mpi
#endif

#ifdef HDF5
  use hdf5_interface,  only: h5tclose_f, h5close_f, hdf5_err, &
                             hdf5_tallyresult_t, hdf5_bank_t, h5fclose_f
#endif

  implicit none

contains

!===============================================================================
! FINALIZE_RUN does all post-simulation tasks such as calculating tally
! statistics and writing out tallies
!===============================================================================

  subroutine finalize_run()

    integer :: i

    if (run_mode == MODE_TESTING) then
      if (master) call print_testing()
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif
      return
    end if

    ! Start finalization timer
    call time_finalize % start()

    if (run_mode /= MODE_PLOTTING .and. run_mode /= MODE_PARTICLE &
        .and. run_mode /= MODE_DISTRIBUTION) then
      ! Calculate statistics for tallies and write to tallies.out
      if (master .or. (dd_run .and. domain_decomp % local_master)) then
        if (n_realizations > 1) call tally_statistics()
        if (output_tallies) call write_tallies()
      end if
      if (check_overlaps) call reduce_overlap_count()
      if (output_distribmats) call write_distribmat_comps(OUTPUT_MATFILE)
    end if

    ! Stop timers and show timing statistics
    call time_finalize % stop()
    call time_total % stop()
    if (master .and. (run_mode /= MODE_PLOTTING .and. &
         run_mode /= MODE_PARTICLE)) then
      call print_runtime()
      call print_results()
      if (check_overlaps) call print_overlap_check()
    end if

    ! Print number of domain interactions
    if (dd_run .and. domain_decomp % count_interactions) then
#ifdef MPI
      call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
#endif
      if (domain_decomp % local_master) call print_domain_interactions()
    end if

    ! Clearn up OTF file reading
    do i = 1, n_materials
      if (materials(i) % otf_compositions) then
        call materials(i) % comp_file % close()
      end if
    end do
    ! Close the OTF file
    if (otf_matfile_open) then
      do i = 1, n_materials
        if (materials(i) % otf_compositions) then
          call h5fclose_f(materials(i) % comp_file % file_id, hdf5_err)
          otf_matfile_open = .false.
          exit
        end if
      end do
    end if

    ! Deallocate arrays
    call free_memory()

#ifdef HDF5
    ! Release compound datatypes
    call h5tclose_f(hdf5_tallyresult_t, hdf5_err)
    call h5tclose_f(hdf5_bank_t, hdf5_err)

    ! Close FORTRAN interface.
    call h5close_f(hdf5_err)
#endif

#ifdef MPI
    ! Free all MPI types
    call MPI_TYPE_FREE(MPI_BANK, mpi_err)
    call MPI_TYPE_FREE(MPI_TALLYRESULT, mpi_err)

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
