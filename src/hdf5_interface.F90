module hdf5_interface

#ifdef HDF5

  use hdf5
  use h5lt
  use, intrinsic :: ISO_C_BINDING
  
#ifdef MPI
   use mpi, only: MPI_COMM_WORLD, MPI_INFO_NULL
#endif

  implicit none

  integer          :: hdf5_err   ! HDF5 error code
  integer          :: hdf5_rank  ! rank of data
  integer(HID_T)   :: dset       ! data set handle
  integer(HID_T)   :: dspace     ! data or file space handle
  integer(HID_T)   :: memspace   ! data space handle for individual procs
  integer(HID_T)   :: plist      ! property list handle
  integer(HSIZE_T) :: dims1(1)   ! dims type for 1-D array
  integer(HSIZE_T) :: dims2(2)   ! dims type for 2-D array
  integer(HSIZE_T) :: dims3(3)   ! dims type for 3-D array
  integer(HSIZE_T) :: dims4(4)   ! dims type for 4-D array
  type(c_ptr)      :: f_ptr      ! pointer to data

  ! Generic HDF5 write procedure interface
  interface hdf5_write_data
    module procedure hdf5_write_double
    module procedure hdf5_write_double_1Darray
    module procedure hdf5_write_double_2Darray
    module procedure hdf5_write_double_3Darray
    module procedure hdf5_write_double_4Darray
    module procedure hdf5_write_integer
    module procedure hdf5_write_integer_1Darray
    module procedure hdf5_write_integer_2Darray
    module procedure hdf5_write_integer_3Darray
    module procedure hdf5_write_integer_4Darray
    module procedure hdf5_write_long
    module procedure hdf5_write_string
#ifdef MPI
    module procedure hdf5_write_double_parallel
    module procedure hdf5_write_double_1Darray_parallel
    module procedure hdf5_write_double_2Darray_parallel
    module procedure hdf5_write_double_3Darray_parallel
    module procedure hdf5_write_double_4Darray_parallel
    module procedure hdf5_write_integer_parallel
    module procedure hdf5_write_integer_1Darray_parallel
    module procedure hdf5_write_integer_2Darray_parallel
    module procedure hdf5_write_integer_3Darray_parallel
    module procedure hdf5_write_integer_4Darray_parallel
    module procedure hdf5_write_long_parallel
    module procedure hdf5_write_string_parallel
#endif
  end interface hdf5_write_data

  ! Generic HDF5 read procedure interface
  interface hdf5_read_data
    module procedure hdf5_read_double
    module procedure hdf5_read_double_1Darray
    module procedure hdf5_read_double_2Darray
    module procedure hdf5_read_double_3Darray
    module procedure hdf5_read_double_4Darray
    module procedure hdf5_read_integer
    module procedure hdf5_read_integer_1Darray
    module procedure hdf5_read_integer_2Darray
    module procedure hdf5_read_integer_3Darray
    module procedure hdf5_read_integer_4Darray
    module procedure hdf5_read_long
    module procedure hdf5_read_string
#ifdef MPI
    module procedure hdf5_read_double_parallel
    module procedure hdf5_read_double_1Darray_parallel
    module procedure hdf5_read_double_2Darray_parallel
    module procedure hdf5_read_double_3Darray_parallel
    module procedure hdf5_read_double_4Darray_parallel
    module procedure hdf5_read_integer_parallel
    module procedure hdf5_read_integer_1Darray_parallel
    module procedure hdf5_read_integer_2Darray_parallel
    module procedure hdf5_read_integer_3Darray_parallel
    module procedure hdf5_read_integer_4Darray_parallel
    module procedure hdf5_read_long_parallel
    module procedure hdf5_read_string_parallel
#endif
  end interface hdf5_read_data

contains

!===============================================================================
! HDF5_FILE_CREATE creates HDF5 file
!===============================================================================

  subroutine hdf5_file_create(filename, file_id)

    character(*),   intent(in)    :: filename ! name of file
    integer(HID_T), intent(inout) :: file_id  ! file handle

    ! Create the file
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err)

  end subroutine hdf5_file_create

!===============================================================================
! HDF5_FILE_OPEN opens HDF5 file
!===============================================================================

  subroutine hdf5_file_open(filename, file_id, mode)

    character(*),  intent(in)      :: filename ! name of file
    character(*),  intent(in)      :: mode     ! access mode to file
    integer(HID_T), intent(inout)  :: file_id  ! file handle

    integer :: open_mode ! HDF5 open mode

    ! Determine access type
    open_mode = H5F_ACC_RDONLY_F
    if (trim(mode) == 'w') then
      open_mode = H5F_ACC_RDWR_F
    end if

    ! Open file
    call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err)

  end subroutine hdf5_file_open

!===============================================================================
! HDF5_FILE_CLOSE closes HDF5 file
!===============================================================================

  subroutine hdf5_file_close(file_id)

    integer(HID_T), intent(inout) :: file_id ! file handle

    ! Close the file
    call h5fclose_f(file_id, hdf5_err)

  end subroutine hdf5_file_close

#ifdef MPI

!===============================================================================
! HDF5_FILE_CREATE_PARALLEL creates HDF5 file with parallel I/O
!===============================================================================

  subroutine hdf5_file_create_parallel(filename, file_id)

    character(*),   intent(in)    :: filename ! name of file
    integer(HID_T), intent(inout) :: file_id  ! file handle

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdf5_err)
    call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_err)

    ! Create the file collectively
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err, &
                     access_prp = plist)

    ! Close the property list
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_file_create_parallel

!===============================================================================
! HDF5_FILE_OPEN_PARALLEL opens HDF5 file with parallel I/O
!===============================================================================

  subroutine hdf5_file_open_parallel(filename, file_id, mode)

    character(*),  intent(in)     :: filename ! name of file
    character(*),  intent(in)     :: mode     ! access mode
    integer(HID_T), intent(inout) :: file_id  ! file handle

    integer        :: open_mode ! HDF5 access mode

    ! Setup file access property list with parallel I/O access
    call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdf5_err)
    call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_err)

    ! Determine access type
    open_mode = H5F_ACC_RDONLY_F
    if (trim(mode) == 'w') then
      open_mode = H5F_ACC_RDWR_F
    end if

    ! Create the file collectively
    call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err, &
                   access_prp = plist)

    ! Close the property list
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_file_open_parallel

#endif

!===============================================================================
! HDF5_OPEN_GROUP creates/opens HDF5 group to temp_group
!===============================================================================

  subroutine hdf5_open_group(hdf5_fh, group, hdf5_grp)

    character(*),   intent(in)    :: group    ! name of group
    integer(HID_T), intent(in)    :: hdf5_fh  ! file handle of main output file
    integer(HID_T), intent(inout) :: hdf5_grp ! handle for group

    logical :: status ! does the group exist

    ! Check if group exists
    call h5ltpath_valid_f(hdf5_fh, trim(group), .true., status, hdf5_err) 

    ! Either create or open group
    if (status) then
      call h5gopen_f(hdf5_fh, trim(group), hdf5_grp, hdf5_err)
    else
      call h5gcreate_f(hdf5_fh, trim(group), hdf5_grp, hdf5_err)
    end if

  end subroutine hdf5_open_group

!===============================================================================
! HDF5_CLOSE_GROUP closes HDF5 temp_group
!===============================================================================

  subroutine hdf5_close_group(hdf5_grp)

    integer(HID_T), intent(inout) :: hdf5_grp

    ! Close the group
    call h5gclose_f(hdf5_grp, hdf5_err)

  end subroutine hdf5_close_group

!===============================================================================
! HDF5_WRITE_INTEGER writes integer scalar data
!===============================================================================

  subroutine hdf5_write_integer(group, name, buffer)

    integer(HID_T), intent(in) :: group  ! name of group
    character(*),   intent(in) :: name   ! name of data
    integer,        intent(in) :: buffer ! data to write

    ! Set rank and dimensions
    hdf5_rank = 1
    dims1(1) = 1

    call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims1, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_write_integer

!===============================================================================
! HDF5_READ_INTEGER reads integer scalar data
!===============================================================================

  subroutine hdf5_read_integer(group, name, buffer)

    integer(HID_T), intent(in)    :: group  ! name of group
    character(*),   intent(in)    :: name   ! name of data
    integer,        intent(inout) :: buffer ! read data to here 

    integer :: buffer_copy(1) ! need an array for read

    ! Set up dimensions
    dims1(1) = 1

    ! Read data
    call h5ltread_dataset_int_f(group, name, buffer_copy, dims1, hdf5_err)
    buffer = buffer_copy(1)

  end subroutine hdf5_read_integer

!===============================================================================
! HDF5_WRITE_INTEGER_1DARRAY writes integer 1-D array
!===============================================================================

  subroutine hdf5_write_integer_1Darray(group, name, buffer, len)

    integer,        intent(in) :: len       ! length of array to write
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    integer,        intent(in) :: buffer(:) ! data to write

    ! Set rank and dimensions of data
    hdf5_rank = 1
    dims1(1) = len

    ! Write data
    call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims1, &
         buffer, hdf5_err)

  end subroutine hdf5_write_integer_1Darray

!===============================================================================
! HDF5_READ_INTEGER_1DARRAY reads integer 1-D array
!===============================================================================

  subroutine hdf5_read_integer_1Darray(group, name, buffer, length)

    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    integer,        intent(inout) :: buffer(:) ! read data to here
    integer,        intent(in)    :: length    ! length of array

    ! Set dimensions
    dims1(1) = length

    ! Read data
    call h5ltread_dataset_int_f(group, name, buffer, dims1, hdf5_err)

  end subroutine hdf5_read_integer_1Darray

!===============================================================================
! HDF5_WRITE_INTEGER_2DARRAY writes integer 2-D array
!===============================================================================

  subroutine hdf5_write_integer_2Darray(group, name, buffer, length)

    integer,        intent(in) :: length(2) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    integer,        intent(in) :: buffer(length(1),length(2)) ! data to write

    ! Set rank and dimensions
    hdf5_rank = 2
    dims2 = length

    ! Write data
    call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims2, &
         buffer, hdf5_err)

  end subroutine hdf5_write_integer_2Darray

!===============================================================================
! HDF5_READ_INTEGER_2DARRAY reads integer 2-D array
!===============================================================================

  subroutine hdf5_read_integer_2Darray(group, name, buffer, length)

    integer,        intent(in)    :: length(2) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    integer,        intent(inout) :: buffer(length(1),length(2)) ! data to read

    ! Set rank and dimensions
    dims2 = length

    ! Write data
    call h5ltread_dataset_int_f(group, name, buffer, dims2, hdf5_err)

  end subroutine hdf5_read_integer_2Darray

!===============================================================================
! HDF5_WRITE_INTEGER_3DARRAY writes integer 3-D array
!===============================================================================

  subroutine hdf5_write_integer_3Darray(group, name, buffer, length)

    integer,        intent(in) :: length(3) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    integer,        intent(in) :: buffer(length(1),length(2), &
                                         length(3)) ! data to write

    ! Set rank and dimensions
    hdf5_rank = 3
    dims3 = length

    ! Write data
    call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims3, &
         buffer, hdf5_err)

  end subroutine hdf5_write_integer_3Darray

!===============================================================================
! HDF5_READ_INTEGER_3DARRAY reads integer 3-D array
!===============================================================================

  subroutine hdf5_read_integer_3Darray(group, name, buffer, length)

    integer,        intent(in)    :: length(3) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    integer,        intent(inout) :: buffer(length(1),length(2), &
                                            length(3)) ! data to read

    ! Set rank and dimensions
    dims3 = length

    ! Write data
    call h5ltread_dataset_int_f(group, name, buffer, dims3, hdf5_err)

  end subroutine hdf5_read_integer_3Darray

!===============================================================================
! HDF5_WRITE_INTEGER_4DARRAY writes integer 4-D array
!===============================================================================

  subroutine hdf5_write_integer_4Darray(group, name, buffer, length)

    integer,        intent(in)    :: length(4) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    integer,        intent(in)    :: buffer(length(1),length(2), &
                                            length(3),length(4)) ! data to write

    ! Set rank and dimensions
    hdf5_rank = 4
    dims4 = length

    ! Write data
    call h5ltmake_dataset_int_f(group, name, hdf5_rank, dims4, &
         buffer, hdf5_err)

  end subroutine hdf5_write_integer_4Darray

!===============================================================================
! HDF5_READ_INTEGER_4DARRAY reads integer 4-D array
!===============================================================================

  subroutine hdf5_read_integer_4Darray(group, name, buffer, length)

    integer,        intent(in)    :: length(4) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    integer,        intent(inout) :: buffer(length(1),length(2), &
                                            length(3),length(4)) ! data to read

    ! Set rank and dimensions
    dims4 = length

    ! Write data
    call h5ltread_dataset_int_f(group, name, buffer, dims4, hdf5_err)

  end subroutine hdf5_read_integer_4Darray

!===============================================================================
! HDF5_WRITE_DOUBLE writes integer scalar data
!===============================================================================

  subroutine hdf5_write_double(group, name, buffer)

    integer(HID_T), intent(in) :: group  ! name of group
    character(*),   intent(in) :: name   ! name of data
    real(8),        intent(in) :: buffer ! data to write

    ! Set rank and dimensions
    hdf5_rank = 1
    dims1(1) = 1

    call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims1, &
         (/ buffer /), hdf5_err)

  end subroutine hdf5_write_double

!===============================================================================
! HDF5_READ_DOUBLE reads double scalar data
!===============================================================================

  subroutine hdf5_read_double(group, name, buffer)

    integer(HID_T), intent(in)    :: group  ! name of group
    character(*),   intent(in)    :: name   ! name of data
    real(8),        intent(inout) :: buffer ! read data to here 

    real(8) :: buffer_copy(1) ! need an array for read

    ! Set up dimensions
    dims1(1) = 1

    ! Read data
    call h5ltread_dataset_double_f(group, name, buffer_copy, dims1, hdf5_err)
    buffer = buffer_copy(1)

  end subroutine hdf5_read_double

!===============================================================================
! HDF5_WRITE_DOUBLE_1DARRAY writes double 1-D array
!===============================================================================

  subroutine hdf5_write_double_1Darray(group, name, buffer, length)

    integer,        intent(in) :: length    ! length of array to write
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    real(8),        intent(in) :: buffer(:) ! data to write

    ! Set rank and dimensions of data
    hdf5_rank = 1
    dims1(1) = length

    ! Write data
    call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims1, &
         buffer, hdf5_err)

  end subroutine hdf5_write_double_1Darray

!===============================================================================
! HDF5_READ_DOUBLE_1DARRAY reads double 1-D array
!===============================================================================

  subroutine hdf5_read_double_1Darray(group, name, buffer, length)

    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    real(8),        intent(inout) :: buffer(:) ! read data to here
    integer,        intent(in)    :: length    ! length of array

    ! Set dimensions
    dims1(1) = length

    ! Read data
    call h5ltread_dataset_double_f(group, name, buffer, dims1, hdf5_err)

  end subroutine hdf5_read_double_1Darray

!===============================================================================
! HDF5_WRITE_DOUBLE_2DARRAY writes double 2-D array
!===============================================================================

  subroutine hdf5_write_double_2Darray(group, name, buffer, length)

    integer,        intent(in) :: length(2) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    real(8),        intent(in) :: buffer(length(1),length(2)) ! data to write

    ! Set rank and dimensions
    hdf5_rank = 2
    dims2 = length

    ! Write data
    call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims2, &
         buffer, hdf5_err)

  end subroutine hdf5_write_double_2Darray

!===============================================================================
! HDF5_READ_DOUBLE_2DARRAY reads double 2-D array
!===============================================================================

  subroutine hdf5_read_double_2Darray(group, name, buffer, length)

    integer,        intent(in)    :: length(2) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    real(8),        intent(inout) :: buffer(length(1),length(2)) ! data to read

    ! Set rank and dimensions
    dims2 = length

    ! Write data
    call h5ltread_dataset_double_f(group, name, buffer, dims2, hdf5_err)

  end subroutine hdf5_read_double_2Darray

!===============================================================================
! HDF5_WRITE_DOUBLE_3DARRAY writes double 3-D array
!===============================================================================

  subroutine hdf5_write_double_3Darray(group, name, buffer, length)

    integer,        intent(in) :: length(3) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    real(8),        intent(in) :: buffer(length(1),length(2), &
                                         length(3)) ! data to write

    ! Set rank and dimensions
    hdf5_rank = 3
    dims3 = length

    ! Write data
    call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims3, &
         buffer, hdf5_err)

  end subroutine hdf5_write_double_3Darray

!===============================================================================
! HDF5_READ_DOUBLE_3DARRAY reads double 3-D array
!===============================================================================

  subroutine hdf5_read_double_3Darray(group, name, buffer, length)

    integer,        intent(in)    :: length(3) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    real(8),        intent(inout) :: buffer(length(1),length(2), &
                                            length(3)) ! data to read

    ! Set rank and dimensions
    dims3 = length

    ! Write data
    call h5ltread_dataset_double_f(group, name, buffer, dims3, hdf5_err)

  end subroutine hdf5_read_double_3Darray

!===============================================================================
! HDF5_WRITE_DOUBLE_4DARRAY writes double 4-D array
!===============================================================================

  subroutine hdf5_write_double_4Darray(group, name, buffer, length)

    integer,        intent(in)    :: length(4) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    real(8),        intent(in)    :: buffer(length(1),length(2), &
                                            length(3),length(4)) ! data to write

    ! Set rank and dimensions
    hdf5_rank = 4
    dims4 = length

    ! Write data
    call h5ltmake_dataset_double_f(group, name, hdf5_rank, dims4, &
         buffer, hdf5_err)

  end subroutine hdf5_write_double_4Darray

!===============================================================================
! HDF5_READ_DOUBLE_4DARRAY reads double 4-D array
!===============================================================================

  subroutine hdf5_read_double_4Darray(group, name, buffer, length)

    integer,        intent(in)    :: length(4) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    real(8),        intent(inout) :: buffer(length(1),length(2), &
                                            length(3),length(4)) ! data to read

    ! Set rank and dimensions
    dims4 = length

    ! Write data
    call h5ltread_dataset_double_f(group, name, buffer, dims4, hdf5_err)

  end subroutine hdf5_read_double_4Darray

!===============================================================================
! HDF5_WRITE_LONG writes long integer scalar data
!===============================================================================

  subroutine hdf5_write_long(group, name, buffer, long_type)

    integer(HID_T),     intent(in) :: group     ! name of group
    character(*),       intent(in) :: name      ! name of data
    integer(8), target, intent(in) :: buffer    ! data to write
    integer(HID_T),     intent(in) :: long_type ! HDF5 long type

    ! Set up rank and dimensions
    hdf5_rank = 1
    dims1(1) = 1

    ! Create dataspace and dataset
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)
    call h5dcreate_f(group, name, long_type, dspace, dset, hdf5_err)

    ! Write eight-byte integer
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, long_type, f_ptr, hdf5_err)

    ! Close dataspace and dataset for long integer
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

  end subroutine hdf5_write_long

!===============================================================================
! HDF5_READ_LONG read long integer scalar data
!===============================================================================

  subroutine hdf5_read_long(group, name, buffer, long_type)

    integer(HID_T),     intent(in)  :: group     ! name of group
    character(*),       intent(in)  :: name      ! name of data
    integer(8), target, intent(out) :: buffer    ! read data to here
    integer(HID_T),     intent(in)  :: long_type ! long integer type

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Get pointer to buffer
    f_ptr = c_loc(buffer)

    ! Read data from dataset
    call h5dread_f(dset, long_type, f_ptr, hdf5_err)

    ! Close dataset
    call h5dclose_f(dset, hdf5_err)

  end subroutine hdf5_read_long

!===============================================================================
! HDF5_WRITE_STRING writes string data
!===============================================================================

  subroutine hdf5_write_string(group, name, buffer, length)

    integer(HID_T), intent(in)    :: group  ! name of group
    character(*),   intent(in)    :: name   ! name of data
    character(*),   intent(in)    :: buffer ! data to write
    integer,        intent(in)    :: length

    character(len=length), dimension(1) :: str_tmp

!   Fortran 2003 implementation not compatible with IBM compiler Feb 2013
!   type(c_ptr), dimension(1), target :: wdata
!   character(len=length, kind=c_char), dimension(1), target :: c_str
!   dims1(1) = 1
!   call h5screate_simple_f(1, dims1, dspace, hdf5_err)
!   call h5dcreate_f(group, name, H5T_STRING, dspace, dset, hdf5_err)
!   c_str(1) = buffer
!   wdata(1) = c_loc(c_str(1))
!   f_ptr = c_loc(wdata(1))

    ! Number of strings to write
    dims1(1) = 1

    ! Insert null character at end of string when writing
    call h5tset_strpad_f(H5T_STRING, H5T_STR_NULLPAD_F, hdf5_err)

    ! Create the dataspace and dataset
    call h5screate_simple_f(1, dims1, dspace, hdf5_err)
    call h5dcreate_f(group, name, H5T_STRING, dspace, dset, hdf5_err)

    ! Set up dimesnions of string to write
    dims2 = (/length, 1/) ! full array of strings to write 
    dims1(1) = length     ! length of string

    ! Copy over string buffer to a rank 1 array
    str_tmp(1) = buffer

    ! Write the variable dataset
    call h5dwrite_vl_f(dset, H5T_STRING, str_tmp, dims2, dims1, hdf5_err, &
         mem_space_id=dspace)

    ! Close all
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

  end subroutine hdf5_write_string

!===============================================================================
! HDF5_READ_STRING reads string data
!===============================================================================

  subroutine hdf5_read_string(group, name, buffer, length)

    integer(HID_T), intent(in)    :: group  ! name of group
    character(*),   intent(in)    :: name   ! name of data
    character(*),   intent(inout) :: buffer ! read data to here
    integer,        intent(in)    :: length ! length of string to read

    character(len=length), dimension(1) :: str_tmp

    ! Fortran 2003 implementation not compatible with IBM Feb 2013 compiler
!    type(c_ptr), dimension(1), target :: buf_ptr
!    character(len=length, kind=c_char), pointer :: chr_ptr
!    f_ptr = c_loc(buf_ptr(1))
!    call h5dread_f(dset, H5T_STRING, f_ptr, hdf5_err, xfer_prp=plist)
!    call c_f_pointer(buf_ptr(1), chr_ptr) 
!    buffer = chr_ptr
!    nullify(chr_ptr)

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Get dataspace to read
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Set dimensions
    dims2 = (/length, 1/)
    dims1(1) = length

    ! Read in the data
    call h5dread_vl_f(dset, H5T_STRING, str_tmp, dims2, dims1, hdf5_err, &
         mem_space_id=dspace, xfer_prp = plist)

    ! Copy over buffer
    buffer = str_tmp(1)

    ! Close dataset
    call h5dclose_f(dset, hdf5_err)

  end subroutine hdf5_read_string

!===============================================================================
! HDF5_WRITE_ATTRIBUTE_STRING writes a string attribute to a variables
!===============================================================================

  subroutine hdf5_write_attribute_string(group, var, attr_type, attr_str)

    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: var       ! name of varaible to set attr
    character(*),   intent(in) :: attr_type ! the attr type id
    character(*),   intent(in) :: attr_str  ! attribute sting

    call h5ltset_attribute_string_f(group, var, attr_type, attr_str, hdf5_err)

  end subroutine hdf5_write_attribute_string

# ifdef MPI

!===============================================================================
! HDF5_WRITE_INTEGER_PARALLEL writes integer scalar data in parallel
!===============================================================================

  subroutine hdf5_write_integer_parallel(group, name, buffer, collect)

    integer(HID_T), intent(in) :: group   ! name of group
    character(*),   intent(in) :: name    ! name of data
    integer,target, intent(in) :: buffer  ! data to write
    logical,        intent(in) :: collect ! collect I/O

    ! Set rank and dimensions
    hdf5_rank = 1
    dims1(1) = 1

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, H5T_NATIVE_INTEGER, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_integer_parallel

!===============================================================================
! HDF5_READ_INTEGER_PARALLEL reads integer scalar data
!===============================================================================

  subroutine hdf5_read_integer_parallel(group, name, buffer, collect)

    integer(HID_T),  intent(in)    :: group  ! name of group
    character(*),    intent(in)    :: name   ! name of data
    integer, target, intent(inout) :: buffer ! read data to here
    logical,         intent(in)    :: collect ! collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_integer_parallel

!===============================================================================
! HDF5_WRITE_INTEGER_1DARRAY_PARALLEL writes integer 1-D array in parallel
!===============================================================================

  subroutine hdf5_write_integer_1Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in) :: length    ! length of array to write
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    integer,target, intent(in) :: buffer(length) ! data to write
    logical,        intent(in) :: collect   ! collect I/O

    ! Set rank and dimensions of data
    hdf5_rank = 1
    dims1(1) = length

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, H5T_NATIVE_INTEGER, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_integer_1Darray_parallel

!===============================================================================
! HDF5_WRITE_INTEGER_1DARRAY_PARALLEL reads integer 1-D array in parallel
!===============================================================================

  subroutine hdf5_read_integer_1Darray_parallel(group, name, buffer, length, &
             collect)

    integer,         intent(in)    :: length     ! length of array
    integer(HID_T),  intent(in)    :: group      ! name of group
    character(*),    intent(in)    :: name       ! name of data
    integer, target, intent(inout) :: buffer(length)  ! read data to here
    logical,         intent(in)    :: collect    ! collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_integer_1Darray_parallel

!===============================================================================
! HDF5_WRITE_INTEGER_2DARRAY_PARALLEL writes integer 2-D array in parallel
!===============================================================================

  subroutine hdf5_write_integer_2Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in) :: length(2) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    integer,target, intent(in) :: buffer(length(1),length(2)) ! data to write
    logical,        intent(in) :: collect ! collective I/O

    ! Set rank and dimensions
    hdf5_rank = 2
    dims2 = length

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims2, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, H5T_NATIVE_INTEGER, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_integer_2Darray_parallel

!===============================================================================
! HDF5_READ_INTEGER_2DARRAY_PARALLEL reads integer 2-D array in parallel
!===============================================================================

  subroutine hdf5_read_integer_2Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in)    :: length(2) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    integer,target, intent(inout) :: buffer(length(1),length(2)) ! data to read
    logical,        intent(in)    :: collect ! collect I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_integer_2Darray_parallel

!===============================================================================
! HDF5_WRITE_INTEGER_3DARRAY_PARALLEL writes integer 3-D array in parallel
!===============================================================================

  subroutine hdf5_write_integer_3Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in) :: length(3) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    integer,target, intent(in) :: buffer(length(1),length(2), &
                                         length(3)) ! data to write
    logical,        intent(in) :: collect ! collective I/O

    ! Set rank and dimensions
    hdf5_rank = 3
    dims3 = length

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims3, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, H5T_NATIVE_INTEGER, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_integer_3Darray_parallel

!===============================================================================
! HDF5_READ_INTEGER_3DARRAY_PARALLEL reads integer 3-D array in parallel
!===============================================================================

  subroutine hdf5_read_integer_3Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in)    :: length(3) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    integer,target, intent(inout) :: buffer(length(1),length(2), &
                                            length(3)) ! data to read
    logical,        intent(in)    :: collect ! collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_integer_3Darray_parallel

!===============================================================================
! HDF5_WRITE_INTEGER_4DARRAY_PARALLEL writes integer 4-D array in parallel
!===============================================================================

  subroutine hdf5_write_integer_4Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in)    :: length(4) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    integer,target, intent(in)    :: buffer(length(1),length(2), &
                                            length(3),length(4)) ! data to write
    logical,        intent(in)    :: collect ! collective I/O

    ! Set rank and dimensions
    hdf5_rank = 4
    dims4 = length

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims4, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, H5T_NATIVE_INTEGER, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_integer_4Darray_parallel

!===============================================================================
! HDF5_READ_INTEGER_4DARRAY_PARALLEL reads integer 4-D array in parallel
!===============================================================================

  subroutine hdf5_read_integer_4Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in)    :: length(4) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    integer,target, intent(inout) :: buffer(length(1),length(2), &
                                            length(3),length(4)) ! data to read
    logical,        intent(in)    :: collect   ! collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_integer_4Darray_parallel

!===============================================================================
! HDF5_WRITE_DOUBLE_PARALLEL writes double scalar data in parallel
!===============================================================================

  subroutine hdf5_write_double_parallel(group, name, buffer, collect)

    integer(HID_T), intent(in) :: group   ! name of group
    character(*),   intent(in) :: name    ! name of data
    real(8),target, intent(in) :: buffer  ! data to write
    logical,        intent(in) :: collect ! collect I/O

    ! Set rank and dimensions
    hdf5_rank = 1
    dims1(1) = 1

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, H5T_NATIVE_DOUBLE, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_double_parallel

!===============================================================================
! HDF5_READ_DOUBLE_PARALLEL reads double scalar data
!===============================================================================

  subroutine hdf5_read_double_parallel(group, name, buffer, collect)

    integer(HID_T),  intent(in)    :: group  ! name of group
    character(*),    intent(in)    :: name   ! name of data
    real(8), target, intent(inout) :: buffer ! read data to here
    logical,         intent(in)    :: collect ! collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_double_parallel

!===============================================================================
! HDF5_WRITE_DOUBLE_1DARRAY_PARALLEL writes double 1-D array in parallel
!===============================================================================

  subroutine hdf5_write_double_1Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in) :: length    ! length of array to write
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    real(8),target, intent(in) :: buffer(length) ! data to write
    logical,        intent(in) :: collect   ! collect I/O

    ! Set rank and dimensions of data
    hdf5_rank = 1
    dims1(1) = length

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, H5T_NATIVE_DOUBLE, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_double_1Darray_parallel

!===============================================================================
! HDF5_WRITE_DOUBLE_1DARRAY_PARALLEL reads double 1-D array in parallel
!===============================================================================

  subroutine hdf5_read_double_1Darray_parallel(group, name, buffer, length, &
             collect)

    integer,         intent(in)    :: length    ! length of array
    integer(HID_T),  intent(in)    :: group     ! name of group
    character(*),    intent(in)    :: name      ! name of data
    real(8),target,  intent(inout) :: buffer(length) ! read data to here
    logical,         intent(in)    :: collect   ! collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_double_1Darray_parallel

!===============================================================================
! HDF5_WRITE_DOUBLE_2DARRAY_PARALLEL writes double 2-D array in parallel
!===============================================================================

  subroutine hdf5_write_double_2Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in) :: length(2) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    real(8),target, intent(in) :: buffer(length(1),length(2)) ! data to write
    logical,        intent(in) :: collect ! collective I/O

    ! Set rank and dimensions
    hdf5_rank = 2
    dims2 = length

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims2, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, H5T_NATIVE_DOUBLE, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer(1,1))
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_double_2Darray_parallel

!===============================================================================
! HDF5_READ_DOUBLE_2DARRAY_PARALLEL reads double 2-D array in parallel
!===============================================================================

  subroutine hdf5_read_double_2Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in)    :: length(2) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    real(8),target, intent(inout) :: buffer(length(1),length(2)) ! data to read
    logical,        intent(in)    :: collect ! collect I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_double_2Darray_parallel

!===============================================================================
! HDF5_WRITE_DOUBLE_3DARRAY_PARALLEL writes double 3-D array in parallel
!===============================================================================

  subroutine hdf5_write_double_3Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in) :: length(3) ! length of array dimensions
    integer(HID_T), intent(in) :: group     ! name of group
    character(*),   intent(in) :: name      ! name of data
    real(8),target, intent(in) :: buffer(length(1),length(2), &
                                         length(3)) ! data to write
    logical,        intent(in) :: collect ! collective I/O

    ! Set rank and dimensions
    hdf5_rank = 3
    dims3 = length

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims3, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, H5T_NATIVE_DOUBLE, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_double_3Darray_parallel

!===============================================================================
! HDF5_READ_DOUBLE_3DARRAY_PARALLEL reads double 3-D array in parallel
!===============================================================================

  subroutine hdf5_read_double_3Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in)    :: length(3) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    real(8),target, intent(inout) :: buffer(length(1),length(2), &
                                            length(3)) ! data to read
    logical,        intent(in)    :: collect ! collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_double_3Darray_parallel

!===============================================================================
! HDF5_WRITE_DOUBLE_4DARRAY_PARALLEL writes double 4-D array in parallel
!===============================================================================

  subroutine hdf5_write_double_4Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in)    :: length(4) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    real(8),target, intent(in)    :: buffer(length(1),length(2), &
                                            length(3),length(4)) ! data to write
    logical,        intent(in)    :: collect ! collective I/O

    ! Set rank and dimensions
    hdf5_rank = 4
    dims4 = length

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims4, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, H5T_NATIVE_DOUBLE, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_double_4Darray_parallel

!===============================================================================
! HDF5_READ_DOUBLE_4DARRAY_PARALLEL reads double 4-D array in parallel
!===============================================================================

  subroutine hdf5_read_double_4Darray_parallel(group, name, buffer, length, &
             collect)

    integer,        intent(in)    :: length(4) ! length of array dimensions
    integer(HID_T), intent(in)    :: group     ! name of group
    character(*),   intent(in)    :: name      ! name of data
    real(8),target, intent(inout) :: buffer(length(1),length(2), &
                                            length(3),length(4)) ! data to read
    logical,        intent(in)    :: collect   ! collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_double_4Darray_parallel

!===============================================================================
! HDF5_WRITE_LONG_PARALLEL writes long integer scalar data in parallel
!===============================================================================

  subroutine hdf5_write_long_parallel(group, name, buffer, long_type, collect)

    integer(HID_T),     intent(in) :: group     ! name of group
    character(*),       intent(in) :: name      ! name of data
    integer(8), target, intent(in) :: buffer    ! data to write
    integer(HID_T),     intent(in) :: long_type ! HDF5 long type
    logical,            intent(in) :: collect   ! collective I/O

    ! Set up rank and dimensions
    hdf5_rank = 1
    dims1(1) = 1

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Create dataspace
    call h5screate_simple_f(hdf5_rank, dims1, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group, name, long_type, dspace, dset, hdf5_err)

    ! Write data
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, long_type, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close all 
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_long_parallel

!===============================================================================
! HDF5_READ_LONG_PARALLEL read long integer scalar data in parallel
!===============================================================================

  subroutine hdf5_read_long_parallel(group, name, buffer, long_type, collect)

    integer(HID_T),     intent(in)  :: group     ! name of group
    character(*),       intent(in)  :: name      ! name of data
    integer(8), target, intent(out) :: buffer    ! read data to here
    integer(HID_T),     intent(in)  :: long_type ! long integer type
    logical,            intent(in)  :: collect   ! collective I/O

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Read data
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, long_type, f_ptr, hdf5_err, xfer_prp=plist)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_long_parallel

!===============================================================================
! HDF5_WRITE_STRING_PARALLEL writes string data in parallel
!===============================================================================

  subroutine hdf5_write_string_parallel(group, name, buffer, length, collect)

    integer(HID_T), intent(in)    :: group   ! name of group
    character(*),   intent(in)    :: name    ! name of data
    character(*),   intent(in)    :: buffer  ! data to write
    integer,        intent(in)    :: length  ! length of string
    logical,        intent(in)    :: collect ! collective I/O

    character(len=length), dimension(1) :: str_tmp

!   Fortran 2003 implementation not compatible with IBM compiler Feb 2013
!   type(c_ptr), dimension(1), target :: wdata
!   character(len=length, kind=c_char), dimension(1), target :: c_str
!   dims1(1) = 1
!   call h5screate_simple_f(1, dims1, dspace, hdf5_err)
!   call h5dcreate_f(group, name, H5T_STRING, dspace, dset, hdf5_err)
!   c_str(1) = buffer
!   wdata(1) = c_loc(c_str(1))
!   f_ptr = c_loc(wdata(1))

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Number of strings to write
    dims1(1) = 1

    ! Insert null character at end of string when writing
    call h5tset_strpad_f(H5T_STRING, H5T_STR_NULLPAD_F, hdf5_err)

    ! Create the dataspace and dataset
    call h5screate_simple_f(1, dims1, dspace, hdf5_err)
    call h5dcreate_f(group, name, H5T_STRING, dspace, dset, hdf5_err)

    ! Set up dimesnions of string to write
    dims2 = (/length, 1/) ! full array of strings to write 
    dims1(1) = length     ! length of string

    ! Copy over string buffer to a rank 1 array
    str_tmp(1) = buffer

    ! Write the variable dataset
    call h5dwrite_vl_f(dset, H5T_STRING, str_tmp, dims2, dims1, hdf5_err, &
         mem_space_id=dspace, xfer_prp=plist)

    ! Close all
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_write_string_parallel

!===============================================================================
! HDF5_READ_STRING_PARALLEL reads string data in parallel
!===============================================================================

  subroutine hdf5_read_string_parallel(group, name, buffer, length, collect)

    integer(HID_T), intent(in)    :: group   ! name of group
    character(*),   intent(in)    :: name    ! name of data
    character(*),   intent(inout) :: buffer  ! read data to here
    integer,        intent(in)    :: length  ! length of string
    logical,        intent(in)    :: collect ! collective I/O

    character(len=length), dimension(1) :: str_tmp

    ! Fortran 2003 implementation not compatible with IBM Feb 2013 compiler
!    type(c_ptr), dimension(1), target :: buf_ptr
!    character(len=length, kind=c_char), pointer :: chr_ptr
!    f_ptr = c_loc(buf_ptr(1))
!    call h5dread_f(dset, H5T_STRING, f_ptr, hdf5_err, xfer_prp=plist)
!    call c_f_pointer(buf_ptr(1), chr_ptr) 
!    buffer = chr_ptr
!    nullify(chr_ptr)

    ! Create property list for independent or collective read
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)

    ! Set independent or collective option
    if (collect) then
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)
    else
      call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_INDEPENDENT_F, hdf5_err)
    end if

    ! Open dataset
    call h5dopen_f(group, name, dset, hdf5_err)

    ! Get dataspace to read
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Set dimensions
    dims2 = (/length, 1/)
    dims1(1) = length

    ! Read in the data
    call h5dread_vl_f(dset, H5T_STRING, str_tmp, dims2, dims1, hdf5_err, &
         mem_space_id=dspace, xfer_prp = plist)

    ! Copy over buffer
    buffer = str_tmp(1)

    ! Close dataset and property list
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

  end subroutine hdf5_read_string_parallel

# endif

#endif

end module hdf5_interface
