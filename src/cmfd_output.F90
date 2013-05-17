module cmfd_output

#ifdef PETSC

! This modules cleans up cmfd objects and echos the results

  implicit none
  private
  public :: finalize_cmfd

contains

!===============================================================================
! FINALIZE_CMFD
!===============================================================================

  subroutine finalize_cmfd() 

    use global,      only: cmfd, cmfd_write_balance, cmfd_write_hdf5, &
                           master, mpi_err
    use cmfd_header, only: deallocate_cmfd

    ! finalize petsc
    call PetscFinalize(mpi_err)

    ! write out final neutron balance
    if (master .and. cmfd_write_balance) call write_neutron_balance()

    ! deallocate cmfd object
    call deallocate_cmfd(cmfd)

  end subroutine finalize_cmfd 

!===============================================================================
! WRITE_NEUTRON_BALANCE 
!===============================================================================

  subroutine write_neutron_balance()

    use cmfd_data,    only: neutron_balance
    use constants,    only: CMFD_BALANCE, MAX_FILE_LEN
    use global,       only: path_output

    character(MAX_FILE_LEN) :: filename

    filename = trim(path_output) // 'neutron_balance.out'

    ! open file for output
    open(UNIT=CMFD_BALANCE, FILE=filename, ACTION='write')

    ! write out the tally
    call neutron_balance(CMFD_BALANCE) 

    ! close file
    close(CMFD_BALANCE)

  end subroutine write_neutron_balance

#endif

end module cmfd_output
