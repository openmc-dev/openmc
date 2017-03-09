module finalize

  use hdf5, only: h5tclose_f, h5close_f
  use global
  use hdf5_interface, only: hdf5_bank_t
  use message_passing

#ifdef PURXS
  use purxs_api, only:&
       URR_isotopes,&
       URR_num_isotopes
#endif

  implicit none

contains

!===============================================================================
! FINALIZE_RUN does all post-simulation tasks such as calculating tally
! statistics and writing out tallies
!===============================================================================

  subroutine openmc_finalize()

    integer :: hdf5_err

    ! Deallocate arrays
    call free_memory()
#ifdef PURXS
    ! deallocate URR isotopes
    do i = 1, URR_num_isotopes
      call URR_isotopes(i) % dealloc()
    end do
#endif
    if (run_mode /= MODE_PURXS) then

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
    end if

  end subroutine openmc_finalize

end module finalize
