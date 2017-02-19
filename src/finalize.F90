module finalize

  use hdf5, only: h5tclose_f, h5close_f

  use global
  use hdf5_interface, only: hdf5_bank_t
  use message_passing

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

  end subroutine openmc_finalize

end module finalize
