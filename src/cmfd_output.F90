module cmfd_output

! This modules cleans up cmfd objects and outputs neutron balance info

  implicit none
  private
  public :: finalize_cmfd

contains

!===============================================================================
! FINALIZE_CMFD finalizes PETSc and frees memory associated with CMFD
!===============================================================================

  subroutine finalize_cmfd() 

    use global,      only: cmfd, cmfd_write_balance, &
                           master, mpi_err
    use cmfd_header, only: deallocate_cmfd

    ! Finalize PETSc 
#ifdef PETSC
    call PetscFinalize(mpi_err)
# endif

    ! Write out final neutron balance
    if (master .and. cmfd_write_balance) call write_neutron_balance()

    ! Deallocate cmfd object
    call deallocate_cmfd(cmfd)

  end subroutine finalize_cmfd 

!===============================================================================
! WRITE_NEUTRON_BALANCE writes an ASCII neutron balance file out
!===============================================================================

  subroutine write_neutron_balance()

    use cmfd_data,    only: neutron_balance
    use constants,    only: CMFD_BALANCE, MAX_FILE_LEN
    use global,       only: path_output

    character(MAX_FILE_LEN) :: filename

    filename = trim(path_output) // 'neutron_balance.out'

    ! Open file for output
    open(UNIT=CMFD_BALANCE, FILE=filename, ACTION='write')

    ! Write out the tally
    call neutron_balance(CMFD_BALANCE) 

    ! Close file
    close(CMFD_BALANCE)

  end subroutine write_neutron_balance

end module cmfd_output
