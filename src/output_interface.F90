module output_interface

  use constants
  use global
  use tally_header,  only: TallyResult

#ifdef HDF5
  use hdf5_interface
#elif MPI
  use mpi_interface
#endif

  implicit none

#if MPI
  integer :: fh_state_point ! mpi file for state point
  integer :: fh_source      ! mpi file for source
#endif

  interface write_data
    module procedure write_double
    module procedure write_double_1Darray
    module procedure write_integer
    module procedure write_integer_1Darray
    module procedure write_long
    module procedure write_string
  end interface write_data


contains

!===============================================================================
! FILE_CREATE
!===============================================================================

  subroutine file_create(filename, fh_str)

    character(MAX_FILE_LEN) :: filename
    character(MAX_WORD_LEN) :: fh_str

#ifdef HDF5
# ifdef MPI
    if (trim(fh_str) == 'state_point') then
      if(master) call hdf5_file_create(trim(trim(filename) // '.h5'), &
                      hdf5_fh)
   else
     print *,'Opening source'
     call hdf5_parallel_file_create(trim(filename) // '.h5', &
          hdf5_fh)
    endif
# else
   if (trim(fh_str) == 'state_point') then
     call hdf5_file_create(trim(filename) // '.h5', &
          hdf5_fh)
   else
     call hdf5_file_create(trim(filename) // '.h5', &
          hdf5_fh)
   endif
# endif
#elif MPI
   if (trim(fh_str) == 'state_point') then
     call mpi_file_create(trim(filename) // '.binary/' , fh_state_point) 
   else
     call mpi_file_create(trim(filename) // '.binary', fh_source)
   endif
#else
   if (trim(fh_str) == 'state_point') then
     open(UNIT=UNIT_STATE, FILE=trim(filename) // '.binary', STATUS='replace', &
          ACCESS='stream')
   else
     open(UNIT=UNIT_SOURCE, FILE=trim(filename) // '.binary', STATUS='replace', &
          ACCESS='stream')
   endif
#endif

  end subroutine file_create

!===============================================================================
! FILE_CLOSE
!===============================================================================

  subroutine file_close(fh_str)

    character(MAX_WORD_LEN) :: fh_str

#ifdef HDF5
# ifdef MPI
    if (trim(fh_str) == 'state_point') then
     if(master) call hdf5_file_close(hdf5_fh)
   else
     call hdf5_file_close(hdf5_fh)
    endif
# else
     call hdf5_file_close(hdf5_fh)
   endif
# endif
#elif MPI
   if (trim(fh_str) == 'state_point') then
     call mpi_close_file(fh_state_point) 
   else
     call mpi_close_file(fh_source)
   endif
#else
   if (trim(fh_str) == 'state_point') then
     close(UNIT=UNIT_STATE)
   else
     close(UNIT=UNIT_SOURCE)
   endif
#endif

  end subroutine file_close

!===============================================================================
! WRITE_INTEGER
!===============================================================================

  subroutine write_double(buffer, name, group)

    real(8),        intent(in)           :: buffer
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_double(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_double

!===============================================================================
! WRITE_INTEGER_1Darray
!===============================================================================

  subroutine write_double_1Darray(buffer, name, group, length)

    integer,        intent(in)           :: length
    real(8),        intent(in)           :: buffer(:)
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_double_1Darray(temp_group, name, buffer, length)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_double_1Darray

!===============================================================================
! WRITE_INTEGER
!===============================================================================

  subroutine write_integer(buffer, name, group)

    integer,        intent(in)           :: buffer
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_integer(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_integer

!===============================================================================
! WRITE_INTEGER_1Darray
!===============================================================================

  subroutine write_integer_1Darray(buffer, name, group, length)

    integer,        intent(in)           :: length
    integer,        intent(in)           :: buffer(:)
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_integer_1Darray(temp_group, name, buffer, length)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_integer_1Darray

!===============================================================================
! WRITE_LONG
!===============================================================================

  subroutine write_long(buffer, name, group)

    integer(8),     intent(in)           :: buffer
    character(*),   intent(in)           :: name
    character(*),   intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_long(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_long


!===============================================================================
! WRITE_STRING
!===============================================================================

  subroutine write_string(buffer, name, group)

    character(*), intent(in)           :: buffer
    character(*), intent(in)           :: name
    character(*), intent(in), optional :: group

#ifdef HDF5
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    endif
    call hdf5_write_string(temp_group, name, buffer)
    if (present(group)) call hdf5_close_group()
#elif MPI

#else

#endif

  end subroutine write_string

!===============================================================================
! WRITE_GLOBAL_TALLIES
!===============================================================================

  subroutine write_tally_result(buffer, name, group, length)

    character(*),      intent(in), optional :: group
    character(*),      intent(in)           :: name
    integer,           intent(in)           :: length
    type(TallyResult), intent(in)           :: buffer(length)

    integer          :: hdf5_err
    integer(HSIZE_T) :: dims(1)
    integer(HID_T)   :: dspace
    integer(HID_T)   :: dset
    type(c_ptr)      :: f_ptr

#ifdef HDF5

    ! Open up sub-group if present
    if (present(group)) then
      call hdf5_open_group(group)
    else
      temp_group = hdf5_fh
    end if

    ! Set overall size of vector to write
    dims(1) = length

    ! Create up a dataspace for size
    call h5screate_simple_f(1, dims, dspace, hdf5_err)

    ! Create the dataset
    call h5dcreate_f(temp_group, name, hdf5_tallyresult_t, dspace, dset, &
         hdf5_err)

    ! Set pointer to first value and write
    f_ptr = c_loc(buffer(1))
    call h5dwrite_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)

    ! Close ids
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    if (present(group)) then
      call hdf5_close_group()
    end if

#elif MPI

#else

#endif 
   
  end subroutine write_tally_result

!===============================================================================
! WRITE_SOURCE_BANK
!===============================================================================

  subroutine write_source_bank()

#ifdef HDF5
    integer(HSIZE_T) :: dims(1)
    integer(HID_T)   :: dset
    integer(HID_T)   :: dspace
    integer(HID_T)   :: memspace
    integer(HID_T)   :: plist
    integer          :: rank
    integer(8)       :: offset(1)
    type(c_ptr)      :: f_ptr
#endif

#ifdef HDF5
# ifdef MPI

    ! Set size of total dataspace for all procs and rank
    dims(1) = n_particles
    rank = 1

    ! Create that dataspace
    call h5screate_simple_f(rank, dims, dspace, hdf5_err)

    ! Create the dataset for that dataspace
    call h5dcreate_f(hdf5_fh, "source", hdf5_bank_t, dspace, dset, hdf5_err)

    ! Close the dataspace
    call h5sclose_f(dspace, hdf5_err)

    ! Create another data space but for each proc individually
    dims(1) = work
    call h5screate_simple_f(rank, dims, memspace, hdf5_err)

    ! Get the individual local proc dataspace
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Select hyperslab for this dataspace
    offset(1) = bank_first - 1_8
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims, hdf5_err)

    ! Set up the property list for parallel writing
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
    call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank(1))

    ! Write data to file in parallel
    call h5dwrite_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
         file_space_id = dspace, mem_space_id = memspace, &
         xfer_prp = plist)

    ! Close all ids
    call h5sclose_f(dspace, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

# else

# endif 

#elif MPI

#else

#endif

  end subroutine write_source_bank

end module output_interface
