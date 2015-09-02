module hdf5_interface

  ! This module provides the high-level procedures which greatly simplify
  ! writing/reading different types of data to HDF5 files. In order to get it to
  ! work with gfotran 4.6, all the write_<type>_ND subroutines had to be split
  ! into two procedures, one accepting an assumed-shape array and another one
  ! with an explicit-shape array since in gfortran 4.6 C_LOC does not work with
  ! an assumed-shape array. When we move to gfortran 4.9+, these procedures can
  ! be combined into one simply accepting an assumed-shape array.

  use tally_header, only: TallyResult

  use hdf5
  use h5lt
  use, intrinsic :: ISO_C_BINDING

#ifdef PHDF5
  use message_passing, only: MPI_COMM_WORLD, MPI_INFO_NULL
#endif

  implicit none
  private

  integer(HID_T), public :: hdf5_tallyresult_t ! Compound type for TallyResult
  integer(HID_T), public :: hdf5_bank_t        ! Compound type for Bank
  integer(HID_T), public :: hdf5_integer8_t    ! type for integer(8)

  interface write_dataset
    module procedure write_double
    module procedure write_double_1D
    module procedure write_double_2D
    module procedure write_double_3D
    module procedure write_double_4D
    module procedure write_integer
    module procedure write_integer_1D
    module procedure write_integer_2D
    module procedure write_integer_3D
    module procedure write_integer_4D
    module procedure write_long
    module procedure write_string
    module procedure write_tally_result_1D
    module procedure write_tally_result_2D
  end interface write_dataset

  interface read_dataset
    module procedure read_double
    module procedure read_double_1D
    module procedure read_double_2D
    module procedure read_double_3D
    module procedure read_double_4D
    module procedure read_integer
    module procedure read_integer_1D
    module procedure read_integer_2D
    module procedure read_integer_3D
    module procedure read_integer_4D
    module procedure read_long
    module procedure read_string
    module procedure read_tally_result_1D
    module procedure read_tally_result_2D
  end interface read_dataset

  public :: write_dataset
  public :: read_dataset
  public :: file_create
  public :: file_open
  public :: file_close
  public :: open_group
  public :: close_group
  public :: write_source_bank
  public :: read_source_bank
  public :: write_attribute_string

contains

!===============================================================================
! FILE_CREATE creates HDF5 file
!===============================================================================

  function file_create(filename, parallel) result(file_id)
    character(*),      intent(in)    :: filename ! name of file
    logical, optional, intent(in)    :: parallel ! whether to write in serial
    integer(HID_T) :: file_id

    integer(HID_T) :: plist      ! property list handle
    integer        :: hdf5_err   ! HDF5 error code
    logical :: parallel_

    ! Check for serial option
    if (present(parallel)) then
      parallel_ = parallel
    else
      parallel_ = .false.
    end if

    if (parallel_) then
      ! Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdf5_err)
#ifdef PHDF5
#ifdef MPIF08
      call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD%MPI_VAL, &
           MPI_INFO_NULL%MPI_VAL, hdf5_err)
#else
      call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_err)
#endif
#endif

      ! Create the file collectively
      call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err, &
                       access_prp = plist)

      ! Close the property list
      call h5pclose_f(plist, hdf5_err)
    else
      ! Create the file
      call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdf5_err)
    end if

  end function file_create

!===============================================================================
! FILE_OPEN opens HDF5 file
!===============================================================================

  function file_open(filename, mode, parallel) result(file_id)
    character(*),      intent(in)    :: filename ! name of file
    character(*),      intent(in)    :: mode     ! access mode to file
    logical, optional, intent(in)    :: parallel ! whether to write in serial
    integer(HID_T) :: file_id

    logical :: parallel_
    integer(HID_T) :: plist     ! property list handle
    integer        :: hdf5_err  ! HDF5 error code
    integer        :: open_mode ! HDF5 open mode

    ! Check for serial option
    if (present(parallel)) then
      parallel_ = parallel
    else
      parallel_ = .false.
    end if

    ! Determine access type
    open_mode = H5F_ACC_RDONLY_F
    if (trim(mode) == 'w') open_mode = H5F_ACC_RDWR_F

    if (parallel_) then
      ! Setup file access property list with parallel I/O access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist, hdf5_err)
#ifdef PHDF5
#ifdef MPIF08
      call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD%MPI_VAL, &
           MPI_INFO_NULL%MPI_VAL, hdf5_err)
#else
      call h5pset_fapl_mpio_f(plist, MPI_COMM_WORLD, MPI_INFO_NULL, hdf5_err)
#endif
#endif

      ! Open the file collectively
      call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err, &
                     access_prp = plist)

      ! Close the property list
      call h5pclose_f(plist, hdf5_err)
    else
      ! Open file
      call h5fopen_f(trim(filename), open_mode, file_id, hdf5_err)
    end if

  end function file_open

!===============================================================================
! FILE_CLOSE closes HDF5 file
!===============================================================================

  subroutine file_close(file_id)
    integer(HID_T), intent(in) :: file_id

    integer :: hdf5_err

    call h5fclose_f(file_id, hdf5_err)
  end subroutine file_close

!===============================================================================
! OPEN_GROUP opens an existing HDF5 group
!===============================================================================

  function open_group(group_id, name) result(newgroup_id)
    integer(HID_T), intent(in) :: group_id
    character(*),   intent(in) :: name ! name of group
    integer(HID_T) :: newgroup_id

    logical :: exists   ! does the group exist
    integer :: hdf5_err ! HDF5 error code

    ! Check if group exists
    call h5ltpath_valid_f(group_id, trim(name), .true., exists, hdf5_err)

    ! Either create or open group
    if (exists) call h5gopen_f(group_id, trim(name), newgroup_id, hdf5_err)
  end function open_group

!===============================================================================
! CREATE_GROUP creates a new HDF5 group
!===============================================================================

  function create_group(group_id, name) result(newgroup_id)
    integer(HID_T), intent(in) :: group_id
    character(*),   intent(in) :: name ! name of group
    integer(HID_T) :: newgroup_id

    integer :: hdf5_err ! HDF5 error code
    logical :: exists   ! does the group exist

    ! Check if group exists
    call h5ltpath_valid_f(group_id, trim(name), .true., exists, hdf5_err)

    ! create group
    if (.not. exists) &
         call h5gcreate_f(group_id, trim(name), newgroup_id, hdf5_err)
  end function create_group

!===============================================================================
! CLOSE_GROUP closes HDF5 temp_group
!===============================================================================

  subroutine close_group(group_id)
    integer(HID_T), intent(inout) :: group_id

    integer :: hdf5_err ! HDF5 error code

    call h5gclose_f(group_id, hdf5_err)
  end subroutine close_group

!===============================================================================
! WRITE_DOUBLE writes double precision scalar data
!===============================================================================

  subroutine write_double(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    real(8),      intent(in), target   :: buffer  ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up independentive vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Create dataspace and dataset
    call h5screate_f(H5S_SCALAR_F, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_DOUBLE, &
                     dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_double

!===============================================================================
! READ_DOUBLE reads double precision scalar data
!===============================================================================

  subroutine read_double(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    real(8),      intent(inout), target :: buffer  ! read data to here
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_double

!===============================================================================
! WRITE_DOUBLE_1DARRAY writes double precision 1-D array data
!===============================================================================

  subroutine write_double_1D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_double_1D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_double_1D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_double_1D

  subroutine write_double_1D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(dims(1)) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(1, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_DOUBLE, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_double_1D_explicit

!===============================================================================
! READ_DOUBLE_1DARRAY reads double precision 1-D array data
!===============================================================================

  subroutine read_double_1D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(inout), target :: buffer(:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call read_double_1D_explicit(group_id, dims, name, buffer, indep)
    else
      call read_double_1D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine read_double_1D

  subroutine read_double_1D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(inout), target :: buffer(dims(1)) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_double_1D_explicit

!===============================================================================
! WRITE_DOUBLE_2DARRAY writes double precision 2-D array data
!===============================================================================

  subroutine write_double_2D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(2)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_double_2D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_double_2D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_double_2D

  subroutine write_double_2D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in)       :: dims(2)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(dims(1),dims(2))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(2, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_DOUBLE, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_double_2D_explicit

!===============================================================================
! READ_DOUBLE_2DARRAY reads double precision 2-D array data
!===============================================================================

  subroutine read_double_2D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(inout), target :: buffer(:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(2)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call read_double_2D_explicit(group_id, dims, name, buffer, indep)
    else
      call read_double_2D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine read_double_2D

  subroutine read_double_2D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(2)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(inout), target :: buffer(dims(1),dims(2))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_double_2D_explicit

!===============================================================================
! WRITE_DOUBLE_3DARRAY writes double precision 3-D array data
!===============================================================================

  subroutine write_double_3D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(3)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_double_3D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_double_3D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_double_3D

  subroutine write_double_3D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(3)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(dims(1),dims(2),dims(3))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(3, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_DOUBLE, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_double_3D_explicit

!===============================================================================
! READ_DOUBLE_3DARRAY reads double precision 3-D array data
!===============================================================================

  subroutine read_double_3D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(inout), target :: buffer(:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(3)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call read_double_3D_explicit(group_id, dims, name, buffer, indep)
    else
      call read_double_3D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine read_double_3D

  subroutine read_double_3D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(3)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(inout), target :: buffer(dims(1),dims(2),dims(3))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_double_3D_explicit

!===============================================================================
! WRITE_DOUBLE_4DARRAY writes double precision 4-D array data
!===============================================================================

  subroutine write_double_4D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:,:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(4)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_double_4D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_double_4D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_double_4D

  subroutine write_double_4D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(4)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(dims(1),dims(2),dims(3),dims(4))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(4, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_DOUBLE, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_double_4D_explicit

!===============================================================================
! READ_DOUBLE_4DARRAY reads double precision 4-D array data
!===============================================================================

  subroutine read_double_4D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(inout), target :: buffer(:,:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(4)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call read_double_4D_explicit(group_id, dims, name, buffer, indep)
    else
      call read_double_4D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine read_double_4D

  subroutine read_double_4D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(4)
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(inout), target :: buffer(dims(1),dims(2),dims(3),dims(4))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_double_4D_explicit

!===============================================================================
! WRITE_INTEGER writes integer precision scalar data
!===============================================================================

  subroutine write_integer(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    integer,      intent(in), target   :: buffer  ! data to write
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Create dataspace and dataset
    call h5screate_f(H5S_SCALAR_F, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, &
                     dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_integer

!===============================================================================
! READ_INTEGER reads integer precision scalar data
!===============================================================================

  subroutine read_integer(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    integer,      intent(inout), target :: buffer  ! read data to here
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_integer

!===============================================================================
! WRITE_INTEGER_1DARRAY writes integer precision 1-D array data
!===============================================================================

  subroutine write_integer_1D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_integer_1D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_integer_1D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_integer_1D

  subroutine write_integer_1D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(dims(1)) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(1, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_integer_1D_explicit

!===============================================================================
! READ_INTEGER_1DARRAY reads integer precision 1-D array data
!===============================================================================

  subroutine read_integer_1D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(inout), target :: buffer(:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call read_integer_1D_explicit(group_id, dims, name, buffer, indep)
    else
      call read_integer_1D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine read_integer_1D

  subroutine read_integer_1D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(inout), target :: buffer(dims(1)) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_integer_1D_explicit

!===============================================================================
! WRITE_INTEGER_2DARRAY writes integer precision 2-D array data
!===============================================================================

  subroutine write_integer_2D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(2)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_integer_2D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_integer_2D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_integer_2D

  subroutine write_integer_2D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(2)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(dims(1),dims(2))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(2, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_integer_2D_explicit

!===============================================================================
! READ_INTEGER_2DARRAY reads integer precision 2-D array data
!===============================================================================

  subroutine read_integer_2D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(inout), target :: buffer(:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(2)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call read_integer_2D_explicit(group_id, dims, name, buffer, indep)
    else
      call read_integer_2D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine read_integer_2D

  subroutine read_integer_2D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(2)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(inout), target :: buffer(dims(1),dims(2))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_integer_2D_explicit

!===============================================================================
! WRITE_INTEGER_3DARRAY writes integer precision 3-D array data
!===============================================================================

  subroutine write_integer_3D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(3)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_integer_3D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_integer_3D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_integer_3D

  subroutine write_integer_3D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(3)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(dims(1),dims(2),dims(3))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(3, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_integer_3D_explicit

!===============================================================================
! READ_INTEGER_3DARRAY reads integer precision 3-D array data
!===============================================================================

  subroutine read_integer_3D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(inout), target :: buffer(:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(3)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call read_integer_3D_explicit(group_id, dims, name, buffer, indep)
    else
      call read_integer_3D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine read_integer_3D

  subroutine read_integer_3D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(3)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(inout), target :: buffer(dims(1),dims(2),dims(3))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_integer_3D_explicit

!===============================================================================
! WRITE_INTEGER_4DARRAY writes integer precision 4-D array data
!===============================================================================

  subroutine write_integer_4D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:,:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(4)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call write_integer_4D_explicit(group_id, dims, name, buffer, indep)
    else
      call write_integer_4D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine write_integer_4D

  subroutine write_integer_4D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(4)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(dims(1),dims(2),dims(3),dims(4))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5screate_simple_f(4, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_NATIVE_INTEGER, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_integer_4D_explicit

!===============================================================================
! READ_INTEGER_4DARRAY reads integer precision 4-D array data
!===============================================================================

  subroutine read_integer_4D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(inout), target :: buffer(:,:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(4)

    dims(:) = shape(buffer)
    if (present(indep)) then
      call read_integer_4D_explicit(group_id, dims, name, buffer, indep)
    else
      call read_integer_4D_explicit(group_id, dims, name, buffer)
    end if
  end subroutine read_integer_4D

  subroutine read_integer_4D_explicit(group_id, dims, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(4)
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(inout), target :: buffer(dims(1),dims(2),dims(3),dims(4))
    logical,      intent(in), optional :: indep   ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, H5T_NATIVE_INTEGER, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_integer_4D_explicit

!===============================================================================
! WRITE_LONG writes long integer scalar data
!===============================================================================

  subroutine write_long(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    integer(8),   intent(in), target   :: buffer  ! data to write
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Create dataspace and dataset
    call h5screate_f(H5S_SCALAR_F, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), hdf5_integer8_t, &
                     dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_f(dset, hdf5_integer8_t, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_f(dset, hdf5_integer8_t, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_long

!===============================================================================
! READ_LONG reads long integer scalar data
!===============================================================================

  subroutine read_long(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    integer(8),   intent(inout), target :: buffer  ! read data to here
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    type(c_ptr) :: f_ptr

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_f(dset, hdf5_integer8_t, f_ptr, hdf5_err, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_f(dset, hdf5_integer8_t, f_ptr, hdf5_err)
    end if

    call h5dclose_f(dset, hdf5_err)
  end subroutine read_long

!===============================================================================
! WRITE_STRING writes string data
!===============================================================================

  subroutine write_string(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in)           :: buffer  ! read data to here
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: n
    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    integer(HSIZE_T) :: dims1(1)
    integer(HSIZE_T) :: dims2(2)
    type(c_ptr) :: f_ptr
    character(len=len_trim(buffer)), dimension(1) :: str_tmp

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Insert null character at end of string when writing
    call h5tset_strpad_f(H5T_STRING, H5T_STR_NULLPAD_F, hdf5_err)

    ! Create the dataspace and dataset
    dims1(1) = 1
    call h5screate_simple_f(1, dims1, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), H5T_STRING, dspace, dset, hdf5_err)

    ! Set up dimesnions of string to write
    n = len_trim(buffer)
    dims2(:) = [n, 1] ! full array of strings to write
    dims1(1) = n      ! length of string

    ! Copy over string buffer to a rank 1 array
    str_tmp(1) = buffer

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dwrite_vl_f(dset, H5T_STRING, str_tmp, dims2, dims1, hdf5_err, &
           mem_space_id=dspace, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dwrite_vl_f(dset, H5T_STRING, str_tmp, dims2, dims1, hdf5_err, &
           mem_space_id=dspace)
    end if

    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_string

!===============================================================================
! READ_STRING reads string data
!===============================================================================

  subroutine read_string(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(inout)        :: buffer  ! read data to here
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: n
    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist   ! property list
    integer(HID_T) :: dset    ! data set handle
    integer(HID_T) :: dspace  ! data or file space handle
    integer(HSIZE_T) :: dims1(1)
    integer(HSIZE_T) :: dims2(2)
    type(c_ptr) :: f_ptr
    character(len=len_trim(buffer)), dimension(1) :: str_tmp

    ! Set up collective vs. independent I/O
    data_xfer_mode = H5FD_MPIO_COLLECTIVE_F
    if (present(indep)) then
      if (indep) data_xfer_mode = H5FD_MPIO_INDEPENDENT_F
    end if

    ! Set up dimesnions of string to write
    n = len_trim(buffer)
    dims2(:) = [n, 1] ! full array of strings to write
    dims1(1) = n      ! length of string

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    call h5dget_space_f(dset, dspace, hdf5_err)

    if (using_mpio_device(group_id)) then
#ifdef PHDF5
      call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
      call h5pset_dxpl_mpio_f(plist, data_xfer_mode, hdf5_err)
      call h5dread_vl_f(dset, H5T_STRING, str_tmp, dims2, dims1, hdf5_err, &
           mem_space_id=dspace, xfer_prp=plist)
      call h5pclose_f(plist, hdf5_err)
#endif
    else
      call h5dread_vl_f(dset, H5T_STRING, str_tmp, dims2, dims1, hdf5_err, &
           mem_space_id=dspace)
    end if

    ! Copy over buffer
    buffer = str_tmp(1)

    ! Close dataset
    call h5dclose_f(dset, hdf5_err)
  end subroutine read_string

!===============================================================================
! WRITE_ATTRIBUTE_STRING
!===============================================================================

  subroutine write_attribute_string(group_id, var, attr_type, attr_str)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)   :: var       ! variable name for attr
    character(*), intent(in)   :: attr_type ! attr identifier type
    character(*), intent(in)   :: attr_str  ! string for attr id type

    integer :: hdf5_err

    call h5ltset_attribute_string_f(group_id, var, attr_type, attr_str, hdf5_err)
  end subroutine write_attribute_string

!===============================================================================
! WRITE_TALLY_RESULT writes an OpenMC TallyResult type
!===============================================================================

  subroutine write_tally_result_1D(group_id, name, buffer)
    integer(HID_T), intent(in)      :: group_id
    character(*),      intent(in)         :: name      ! name of data
    type(TallyResult), intent(in), target :: buffer(:) ! data to write

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    call write_tally_result_1D_explicit(group_id, dims, name, buffer)
  end subroutine write_tally_result_1D

  subroutine write_tally_result_1D_explicit(group_id, dims, name, buffer)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*),      intent(in)         :: name      ! name of data
    type(TallyResult), intent(in), target :: buffer(dims(1))

    integer :: hdf5_err
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    call h5screate_simple_f(1, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), hdf5_tallyresult_t, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_tally_result_1D_explicit

  subroutine write_tally_result_2D(group_id, name, buffer)
    integer(HID_T), intent(in)      :: group_id
    character(*),      intent(in)         :: name      ! name of data
    type(TallyResult), intent(in), target :: buffer(:,:) ! data to write

    integer(HSIZE_T) :: dims(2)

    dims(:) = shape(buffer)
    call write_tally_result_2D_explicit(group_id, dims, name, buffer)
  end subroutine write_tally_result_2D

  subroutine write_tally_result_2D_explicit(group_id, dims, name, buffer)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(2)
    character(*),      intent(in)         :: name        ! name of data
    type(TallyResult), intent(in), target :: buffer(dims(1),dims(2))

    integer :: hdf5_err
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    type(c_ptr) :: f_ptr

    call h5screate_simple_f(2, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, trim(name), hdf5_tallyresult_t, &
         dspace, dset, hdf5_err)
    f_ptr = c_loc(buffer)
    call h5dwrite_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)
  end subroutine write_tally_result_2D_explicit

!===============================================================================
! READ_TALLY_RESULT reads OpenMC TallyResult data
!===============================================================================

  subroutine read_tally_result_1D(group_id, name, buffer)
    integer(HID_T),    intent(in)            :: group_id
    character(*),      intent(in)            :: name      ! name of data
    type(TallyResult), intent(inout), target :: buffer(:) ! read data here

    integer(HSIZE_T) :: dims(1)

    dims(:) = shape(buffer)
    call read_tally_result_1D_explicit(group_id, dims, name, buffer)
  end subroutine read_tally_result_1D

  subroutine read_tally_result_1D_explicit(group_id, dims, name, buffer)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(1)
    character(*), intent(in) :: name      ! name of data
    type(TallyResult), intent(inout), target :: buffer(dims(1))

    integer :: hdf5_err
    integer(HID_T) :: dset ! data set handle
    type(c_ptr) :: f_ptr

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
  end subroutine read_tally_result_1D_explicit

  subroutine read_tally_result_2D(group_id, name, buffer)
    integer(HID_T), intent(in) :: group_id
    character(*),      intent(in)            :: name      ! name of data
    type(TallyResult), intent(inout), target :: buffer(:,:)

    integer(HSIZE_T) :: dims(2)

    dims(:) = shape(buffer)
    call read_tally_result_2D_explicit(group_id, dims, name, buffer)
  end subroutine read_tally_result_2D

  subroutine read_tally_result_2D_explicit(group_id, dims, name, buffer)
    integer(HID_T), intent(in) :: group_id
    integer(HSIZE_T), intent(in) :: dims(2)
    character(*),      intent(in)            :: name        ! name of data
    type(TallyResult), intent(inout), target :: buffer(dims(1),dims(2))

    integer :: hdf5_err
    integer(HID_T) :: dset ! data set handle
    type(c_ptr) :: f_ptr

    call h5dopen_f(group_id, trim(name), dset, hdf5_err)
    f_ptr = c_loc(buffer)
    call h5dread_f(dset, hdf5_tallyresult_t, f_ptr, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
  end subroutine read_tally_result_2D_explicit

!===============================================================================
! WRITE_SOURCE_BANK writes OpenMC source_bank data
!===============================================================================

  subroutine write_source_bank(group_id)
    use bank_header, only: Bank
    use global, only: n_particles, work, source_bank

    integer(HID_T), intent(in) :: group_id

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data or file space handle
    integer(HID_T) :: memspace ! memory space handle
    integer(HSIZE_T) :: dims(1)
    type(c_ptr) :: f_ptr
#ifdef PHDF5
    integer(HSIZE_T) :: offset(1) ! source data offset
#endif

#ifdef PHDF5
    ! Set size of total dataspace for all procs and rank
    dims(1) = n_particles
    call h5screate_simple_f(1, dims, dspace, hdf5_err)
    call h5dcreate_f(group_id, "source_bank", hdf5_bank_t, dspace, dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

    ! Create another data space but for each proc individually
    dims(1) = work
    call h5screate_simple_f(rank, dims, memspace, hdf5_err)

    ! Get the individual local proc dataspace
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Select hyperslab for this dataspace
    offset(1) = work_index(rank)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims, hdf5_err)

    ! Set up the property list for parallel writing
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
    call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank)

    ! Write data to file in parallel
    call h5dwrite_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
         file_space_id=dspace, mem_space_id=memspace, &
         xfer_prp=plist)

    ! Close all ids
    call h5sclose_f(dspace, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

#else

    ! Set size
    dims(1) = work

    ! Create dataspace
    call h5screate_simple_f(1, dims, dspace, hdf5_err)

    ! Create dataset
    call h5dcreate_f(group_id, "source_bank", hdf5_bank_t, &
         dspace, dset, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank)

    ! Write dataset to file
    call h5dwrite_f(dset, hdf5_bank_t, f_ptr, hdf5_err)

    ! Close all ids
    call h5dclose_f(dset, hdf5_err)
    call h5sclose_f(dspace, hdf5_err)

#endif

  end subroutine write_source_bank

!===============================================================================
! READ_SOURCE_BANK reads OpenMC source_bank data
!===============================================================================

  subroutine read_source_bank(group_id)
    use bank_header, only: Bank
    use global, only: work, source_bank

    integer(HID_T), intent(in) :: group_id

    integer :: hdf5_err
    integer :: data_xfer_mode
    integer(HID_T) :: plist    ! property list
    integer(HID_T) :: dset     ! data set handle
    integer(HID_T) :: dspace   ! data space handle
    integer(HID_T) :: memspace ! memory space handle
    integer(HSIZE_T) :: dims(1)
    type(c_ptr) :: f_ptr
#ifdef PHDF5
    integer(HSIZE_T) :: offset(1) ! offset of data
#endif

#ifdef PHDF5

    ! Open the dataset
    call h5dopen_f(group_id, "source_bank", dset, hdf5_err)

    ! Create another data space but for each proc individually
    dims(1) = work
    call h5screate_simple_f(1, dims, memspace, hdf5_err)

    ! Get the individual local proc dataspace
    call h5dget_space_f(dset, dspace, hdf5_err)

    ! Select hyperslab for this dataspace
    offset(1) = work_index(rank)
    call h5sselect_hyperslab_f(dspace, H5S_SELECT_SET_F, offset, dims, hdf5_err)

    ! Set up the property list for parallel writing
    call h5pcreate_f(H5P_DATASET_XFER_F, plist, hdf5_err)
    call h5pset_dxpl_mpio_f(plist, H5FD_MPIO_COLLECTIVE_F, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank)

    ! Read data from file in parallel
    call h5dread_f(dset, hdf5_bank_t, f_ptr, hdf5_err, &
         file_space_id=dspace, mem_space_id=memspace, &
         xfer_prp=plist)

    ! Close all ids
    call h5sclose_f(dspace, hdf5_err)
    call h5sclose_f(memspace, hdf5_err)
    call h5dclose_f(dset, hdf5_err)
    call h5pclose_f(plist, hdf5_err)

#else

    ! Open dataset
    call h5dopen_f(group_id, "source_bank", dset, hdf5_err)

    ! Set up pointer to data
    f_ptr = c_loc(source_bank)

    ! Read dataset from file
    call h5dread_f(dset, hdf5_bank_t, f_ptr, hdf5_err)

    ! Close all ids
    call h5dclose_f(dset, hdf5_err)

#endif

  end subroutine read_source_bank

  function using_mpio_device(obj_id) result(mpio)
    integer(HID_T), intent(in) :: obj_id
    logical :: mpio

    integer :: hdf5_err
    integer :: driver
    integer(HID_T) :: file_id
    integer(HID_T) :: fapl_id

    ! Determine file that this object is part of
    call h5iget_file_id_f(obj_id, file_id, hdf5_err)

    ! Get file access property list
    call h5fget_access_plist_f(file_id, fapl_id, hdf5_err)

    ! Get low-level driver identifier
    call h5pget_driver_f(fapl_id, driver, hdf5_err)

    ! Close file access property list access
    call h5pclose_f(fapl_id, hdf5_err)

    ! Close file access -- note that this only decreases the reference count so
    ! that the file is not actually closed
    call h5fclose_f(file_id, hdf5_err)

    mpio = (driver == H5FD_MPIO_F)
  end function using_mpio_device

end module hdf5_interface
