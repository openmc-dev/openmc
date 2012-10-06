module cmfd_output

# ifdef PETSC

! This modules cleans up cmfd objects and echos the results

  implicit none
  private
  public :: finalize_cmfd, cmfd_output_summary, write_hdf5_output 

contains

!===============================================================================
! FINALIZE_CMFD
!===============================================================================

  subroutine finalize_cmfd() 

    use global,                only: cmfd,                                     &
                                     cmfd_write_balance, cmfd_write_hdf5,      &
                                     master, mpi_err
    use cmfd_header,           only: deallocate_cmfd

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
    use constants,    only: CMFD_BALANCE

    ! open file for output
    open(unit=CMFD_BALANCE,file="neutron_balance.out")

    ! write out the tally
    call neutron_balance(CMFD_BALANCE) 

    ! close file
    close(CMFD_BALANCE)

  end subroutine write_neutron_balance

!===============================================================================
! CMFD_OUTPUT_SUMMARY 
!===============================================================================

  subroutine cmfd_output_summary()

    use global,       only: cmfd, time_cmfd, cmfd_run_adjoint
    use output,       only: header 
    use timing,       only: timer_get_value 

    real(8) :: elapsed

    ! get elapsed time
    elapsed = timer_get_value(time_cmfd)

    ! print out heading
!    call header('CMFD Execution Summary',1)

    ! write out information
!    call printStdOutCard(OFF=0, DSC='Total time elapsed', RV=elapsed, PS=' s')
!    call printStdOutCard(OFF=0, DSC='Final k-effective', RV=cmfd%keff)
!    if (cmfd_run_adjoint)                                                      &
!      call printStdOutCard(OFF=0, DSC='Adjoint eigenvalue', RV=cmfd%adj_keff)

    end subroutine cmfd_output_summary

!===============================================================================
! WRITE_HDF5_OUTPUT 
!===============================================================================

  subroutine write_hdf5_output()

# ifdef HDF5

    use global,       only: cmfd, current_batch
    use hdf5
    use string,       only: to_str

    character(100) :: grpname 
    character(100) :: path
    integer(HID_T) :: group_id
    integer(HID_T) :: dataspace_id         ! Data space identifier
    integer(HID_T) :: dataset_id           ! Dataset identifier
    integer(HSIZE_T), dimension(1) :: dim1 ! scalar for hdf5 dimensions
    integer(HSIZE_T), dimension(4) :: dim4 ! vector for hdf5 dimensions
    integer(Fortran_Integer) :: hdf5_err   ! error
    integer :: nx         ! number of mesh cells in x direction
    integer :: ny         ! number of mesh cells in y direction
    integer :: nz         ! number of mesh cells in z direction
    integer :: ng         ! number of energy groups
    integer(hid_t) :: dataType
    integer :: r8 = 8

    ! get hdf5 data type
    dataType = h5kind_to_type(r8,H5_REAL_KIND)

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! set group name
    grpname = 'batch'//trim(to_str(current_batch))
    ! create a new group for this active batch
!   call h5gcreate_f(cmfd%file_id,grpname,group_id,hdf5_err)

    ! write out openmc source
    path = trim(grpname)//'/openmc_src'
    dim4 = (/ng,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,hdf5_err)
!   call h5dcreate_f(cmfd%file_id,trim(path),dataType,dataspace_id,   &
!                    dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,dataType,cmfd%openmc_src,dim4,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write out cmfd source
    path = trim(grpname)//'/cmfd_src'
    dim4 = (/ng,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,hdf5_err)
!   call h5dcreate_f(cmfd%file_id,trim(path),dataType,dataspace_id,   &
!                    dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,dataType,cmfd%cmfd_src,dim4,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write out keff
    path = trim(grpname)//'/cmfd_keff'
    dim1 = (/1/)
    call h5screate_simple_f(1,dim1,dataspace_id,hdf5_err)
!   call h5dcreate_f(cmfd%file_id,trim(path),dataType,dataspace_id,   &
!                    dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,dataType,cmfd%keff,dim1,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! close the group
    call h5gclose_f(group_id,hdf5_err)

# endif

  end subroutine write_hdf5_output

# endif

end module cmfd_output
