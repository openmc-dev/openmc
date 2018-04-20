module hdf5_interface

!==============================================================================
! HDF5_INTERFACE -- This module provides the high-level procedures which greatly
! simplify writing/reading different types of data to HDF5 files. In order to
! get it to work with gfotran 4.6, all the write_<type>_ND subroutines had to be
! split into two procedures, one accepting an assumed-shape array and another
! one with an explicit-shape array since in gfortran 4.6 C_LOC does not work
! with an assumed-shape array. When we move to gfortran 4.9+, these procedures
! can be combined into one simply accepting an assumed-shape array.
!==============================================================================

  use, intrinsic :: ISO_C_BINDING

  use hdf5
  use h5lt

  use error, only: fatal_error
#ifdef PHDF5
  use message_passing, only: mpi_intracomm, MPI_INFO_NULL
#endif
  use string, only: to_c_string

  implicit none
  private

  interface write_dataset
    module procedure write_double_0D
    module procedure write_double_1D
    module procedure write_double_2D
    module procedure write_double_3D
    module procedure write_double_4D
    module procedure write_integer_0D
    module procedure write_integer_1D
    module procedure write_integer_2D
    module procedure write_integer_3D
    module procedure write_integer_4D
    module procedure write_long
    module procedure write_string
    module procedure write_string_1D
  end interface write_dataset

  interface read_dataset
    module procedure read_double_0D
    module procedure read_double_1D
    module procedure read_double_2D
    module procedure read_double_3D
    module procedure read_double_4D
    module procedure read_integer_0D
    module procedure read_integer_1D
    module procedure read_integer_2D
    module procedure read_integer_3D
    module procedure read_integer_4D
    module procedure read_long
    module procedure read_string
    module procedure read_string_1D
    module procedure read_complex_2D
  end interface read_dataset

  interface read_attribute
    module procedure read_attribute_double
    module procedure read_attribute_double_1D
    module procedure read_attribute_double_2D
    module procedure read_attribute_integer
    module procedure read_attribute_integer_1D
    module procedure read_attribute_integer_2D
    module procedure read_attribute_string
    module procedure read_attribute_string_1D
    module procedure read_attribute_logical
  end interface read_attribute

  interface write_attribute
    module procedure write_attribute_double
    module procedure write_attribute_double_1D
    module procedure write_attribute_integer
    module procedure write_attribute_integer_1D
    module procedure write_attribute_string
  end interface write_attribute

  public :: write_dataset
  public :: read_dataset
  public :: attribute_exists
  public :: write_attribute
  public :: read_attribute
  public :: file_open
  public :: file_close
  public :: create_group
  public :: object_exists
  public :: open_group
  public :: close_group
  public :: open_dataset
  public :: close_dataset
  public :: get_shape
  public :: get_ndims
  public :: get_groups
  public :: get_datasets
  public :: get_name
  public :: HID_T, HSIZE_T, SIZE_T

  interface
    function attribute_typesize(obj_id, name) result(sz) bind(C)
      import HID_T, C_CHAR, C_SIZE_T
      integer(HID_T), value :: obj_id
      character(kind=C_CHAR), intent(in) :: name(*)
      integer(C_SIZE_T) :: sz
    end function attribute_typesize

    function dataset_typesize(dset) result(sz) bind(C)
      import HID_T, C_SIZE_T
      integer(HID_T), value :: dset
      integer(C_SIZE_T) :: sz
    end function dataset_typesize

    subroutine file_close(file_id) bind(C)
      import HID_T
      integer(HID_T), value :: file_id
    end subroutine file_close

    subroutine get_shape_c(obj_id, dims) bind(C, name='get_shape')
      import HID_T, HSIZE_T
      integer(HID_T), value :: obj_id
      integer(HSIZE_T), intent(out) :: dims(*)
    end subroutine get_shape_c

    subroutine get_shape_attr(obj_id, name, dims) bind(C)
      import HID_T, HSIZE_T, C_CHAR
      integer(HID_T), value :: obj_id
      character(kind=C_CHAR), intent(in) :: name(*)
      integer(HSIZE_T), intent(out) :: dims(*)
    end subroutine get_shape_attr

    subroutine read_double_c(obj_id, name, buffer, indep) &
         bind(C, name='read_double')
      import HID_T, C_DOUBLE, C_BOOL, C_PTR
      integer(HID_T), value :: obj_id
      type(C_PTR), value :: name
      real(C_DOUBLE), intent(out) :: buffer(*)
      logical(C_BOOL), value :: indep
    end subroutine read_double_c

    subroutine read_attr_int_c(obj_id, name, buffer) &
         bind(C, name='read_attr_int')
      import HID_T, C_CHAR, C_INT
      integer(HID_T), value :: obj_id
      character(kind=C_CHAR), intent(in) :: name(*)
      integer(C_INT), intent(out) :: buffer(*)
    end subroutine read_attr_int_c

    subroutine read_attr_double_c(obj_id, name, buffer) &
         bind(C, name='read_attr_double')
      import HID_T, C_CHAR, C_DOUBLE
      integer(HID_T), value :: obj_id
      character(kind=C_CHAR), intent(in) :: name(*)
      real(C_DOUBLE), intent(out) :: buffer(*)
    end subroutine read_attr_double_c

    subroutine read_attr_string_c(obj_id, name, slen, buffer) &
         bind(C, name='read_attr_string')
      import HID_T, C_CHAR, C_SIZE_T
      integer(HID_T), value :: obj_id
      character(kind=C_CHAR), intent(in) :: name(*)
      integer(C_SIZE_T), value :: slen
      character(kind=C_CHAR), intent(out) :: buffer(*)
    end subroutine read_attr_string_c

    subroutine read_int_c(obj_id, name, buffer, indep) &
         bind(C, name='read_int')
      import HID_T, C_INT, C_BOOL, C_PTR
      integer(HID_T), value :: obj_id
      type(C_PTR), value :: name
      integer(C_INT), intent(out) :: buffer(*)
      logical(C_BOOL), value :: indep
    end subroutine read_int_c

    subroutine read_llong_c(obj_id, name, buffer, indep) &
         bind(C, name='read_llong')
      import HID_T, C_INT, C_BOOL, C_PTR, C_LONG_LONG
      integer(HID_T), value :: obj_id
      type(C_PTR), value :: name
      integer(C_LONG_LONG), intent(out) :: buffer(*)
      logical(C_BOOL), value :: indep
    end subroutine read_llong_c

    subroutine read_string_c(obj_id, name, slen, buffer, indep) &
         bind(C, name='read_string')
      import HID_T, C_PTR, C_SIZE_T, C_CHAR, C_BOOL
      integer(HID_T), value :: obj_id
      type(C_PTR), value :: name
      integer(C_SIZE_T), value :: slen
      character(kind=C_CHAR), intent(out) :: buffer(*)
      logical(C_BOOL), value :: indep
    end subroutine read_string_c

    subroutine read_complex_c(obj_id, name, buffer, indep) &
         bind(C, name='read_complex')
      import HID_T, C_PTR, C_DOUBLE_COMPLEX, C_BOOL
      integer(HID_T), value :: obj_id
      type(C_PTR), value :: name
      complex(C_DOUBLE_COMPLEX), intent(out) :: buffer(*)
      logical(C_BOOL), value :: indep
    end subroutine read_complex_c

    subroutine write_attr_double_c(obj_id, ndim, dims, name, buffer) &
         bind(C, name='write_attr_double')
      import HID_T, HSIZE_T, C_INT, C_DOUBLE, C_CHAR
      integer(HID_T), value :: obj_id
      integer(C_INT), value :: ndim
      integer(HSIZE_T), intent(in) :: dims(*)
      character(kind=C_CHAR), intent(in) :: name(*)
      real(C_DOUBLE), intent(in) :: buffer(*)
    end subroutine write_attr_double_c

    subroutine write_attr_int_c(obj_id, ndim, dims, name, buffer) &
         bind(C, name='write_attr_int')
      import HID_T, HSIZE_T, C_INT, C_CHAR
      integer(HID_T), value :: obj_id
      integer(C_INT), value :: ndim
      integer(HSIZE_T), intent(in) :: dims(*)
      character(kind=C_CHAR), intent(in) :: name(*)
      integer(C_INT), intent(in) :: buffer(*)
    end subroutine write_attr_int_c

    subroutine write_attr_string_c(obj_id, name, buffer) &
         bind(C, name='write_attr_string')
      import HID_T, C_CHAR
      integer(HID_T), value :: obj_id
      character(kind=C_CHAR), intent(in) :: name(*)
      character(kind=C_CHAR), intent(in) :: buffer(*)
    end subroutine write_attr_string_c

    subroutine write_double_c(group_id, ndim, dims, name, buffer, indep) &
         bind(C, name='write_double')
      import HID_T, HSIZE_T, C_INT, C_DOUBLE, C_CHAR, C_BOOL
      integer(HID_T), value :: group_id
      integer(C_INT), value :: ndim
      integer(HSIZE_T), intent(in) :: dims(*)
      character(kind=C_CHAR), intent(in) :: name(*)
      real(C_DOUBLE), intent(in) :: buffer(*)
      logical(C_BOOL), value :: indep
    end subroutine write_double_c

    subroutine write_int_c(group_id, ndim, dims, name, buffer, indep) &
         bind(C, name='write_int')
      import HID_T, HSIZE_T, C_INT, C_CHAR, C_BOOL
      integer(HID_T), value :: group_id
      integer(C_INT), value :: ndim
      integer(HSIZE_T), intent(in) :: dims(*)
      character(kind=C_CHAR), intent(in) :: name(*)
      integer(C_INT), intent(in) :: buffer(*)
      logical(C_BOOL), value :: indep
    end subroutine write_int_c

    subroutine write_llong_c(group_id, ndim, dims, name, buffer, indep) &
         bind(C, name='write_llong')
      import HID_T, HSIZE_T, C_INT, C_CHAR, C_BOOL, C_LONG_LONG
      integer(HID_T), value :: group_id
      integer(C_INT), value :: ndim
      integer(HSIZE_T), intent(in) :: dims(*)
      character(kind=C_CHAR), intent(in) :: name(*)
      integer(C_LONG_LONG), intent(in) :: buffer(*)
      logical(C_BOOL), value :: indep
    end subroutine write_llong_c

    subroutine write_string_c(group_id, ndim, dims, slen, name, buffer, indep) &
         bind(C, name='write_string')
      import HID_T, HSIZE_T, C_INT, C_CHAR, C_BOOL, C_SIZE_T
      integer(HID_T), value :: group_id
      integer(C_INT), value :: ndim
      integer(HSIZE_T), intent(in) :: dims(*)
      integer(C_SIZE_T), value :: slen
      character(kind=C_CHAR), intent(in) :: name(*)
      character(kind=C_CHAR), intent(in) :: buffer(*)
      logical(C_BOOL), value :: indep
    end subroutine write_string_c
  end interface

contains

!===============================================================================
! FILE_OPEN opens HDF5 file
!===============================================================================

  function file_open(filename, mode, parallel) result(file_id)
    character(*),      intent(in) :: filename ! name of file
    character,         value :: mode     ! access mode to file
    logical, optional, intent(in) :: parallel ! whether to write in serial
    integer(HID_T) :: file_id

    character(kind=C_CHAR) :: mode_
    logical(C_BOOL) :: parallel_

    interface
      function file_open_c(name, mode, parallel) bind(C, name='file_open') result(file_id)
        import HID_T, C_CHAR, C_BOOL, C_INT
        character(kind=C_CHAR) :: name(*)
        character(kind=C_CHAR), value :: mode
        logical(C_BOOL), value :: parallel
        integer(HID_T) :: file_id
      end function file_open_c
    end interface

    mode_ = mode
    parallel_ = .false.
#ifdef PHDF5
    if (present(parallel)) parallel_ = parallel
#endif

    file_id = file_open_c(to_c_string(filename), mode, parallel_)
  end function file_open

!===============================================================================
! GET_GROUPS Gets a list of all the groups in a given location.
!===============================================================================

  subroutine get_groups(object_id, names)
    integer(HID_T), intent(in)  :: object_id
    character(len=150), allocatable, intent(out) :: names(:)

    integer :: n_members, i, group_count, type
    integer :: hdf5_err
    character(len=150) :: name

    ! Get number of members in this location
    call h5gn_members_f(object_id, './', n_members, hdf5_err)

    ! Get the number of groups
    group_count = 0
    do i = 0, n_members - 1
      call h5gget_obj_info_idx_f(object_id, "./", i, name, type, hdf5_err)
      if (type == H5G_GROUP_F) then
        group_count = group_count + 1
      end if
    end do

    ! Now we can allocate the storage for the ids
    allocate(names(group_count))
    group_count = 0
    do i = 0, n_members - 1
      call h5gget_obj_info_idx_f(object_id, "./", i, name, type, hdf5_err)
      if (type == H5G_GROUP_F) then
        group_count = group_count + 1
        names(group_count) = trim(name)
      end if
    end do

  end subroutine get_groups

!===============================================================================
! ATTRIBUTE_EXISTS checks to see if an attribute exists in the object
!===============================================================================

  function attribute_exists(object_id, name) result(exists)
    integer(HID_T), intent(in) :: object_id
    character(*),   intent(in) :: name ! name of group
    logical :: exists

    interface
      function attribute_exists_c(obj_id, name) result(exists) &
           bind(C, name='attribute_exists')
        import HID_T, C_CHAR, C_BOOL
        integer(HID_T), value :: obj_id
        character(kind=C_CHAR), intent(in) :: name(*)
        logical(C_BOOL) :: exists
      end function attribute_exists_c
    end interface

    exists = attribute_exists_c(object_id, to_c_string(name))
  end function attribute_exists

!===============================================================================
! CHECK_GROUP Checks to see if a group exists in the object
!===============================================================================

  function object_exists(object_id, name) result(exists)
    integer(HID_T), intent(in) :: object_id
    character(*),   intent(in) :: name ! name of group
    logical :: exists

    interface
      function object_exists_c(obj_id, name) result(exists) &
           bind(C, name='object_exists')
        import HID_T, C_CHAR, C_BOOL
        integer(HID_T), value :: obj_id
        character(kind=C_CHAR), intent(in) :: name(*)
        logical(C_BOOL) :: exists
      end function object_exists_c
    end interface

    ! Check if group exists
    exists = object_exists_c(object_id, to_c_string(name))
  end function object_exists

!===============================================================================
! GET_DATASETS Gets a list of all the datasets in a given location.
!===============================================================================

  subroutine get_datasets(object_id, names)
    integer(HID_T), intent(in)  :: object_id
    character(len=150), allocatable, intent(out) :: names(:)

    integer :: n_members, i, dset_count, type
    integer :: hdf5_err
    character(len=150) :: name


    ! Get number of members in this location
    call h5gn_members_f(object_id, './', n_members, hdf5_err)

    ! Get the number of datasets
    dset_count = 0
    do i = 0, n_members - 1
      call h5gget_obj_info_idx_f(object_id, "./", i, name, type, hdf5_err)
      if (type == H5G_DATASET_F ) then
        dset_count = dset_count + 1
      end if
    end do

    ! Now we can allocate the storage for the ids
    allocate(names(dset_count))
    dset_count = 0
    do i = 0, n_members - 1
      call h5gget_obj_info_idx_f(object_id, "./", i, name, type, hdf5_err)
      if (type == H5G_DATASET_F ) then
        dset_count = dset_count + 1
        names(dset_count) = trim(name)
      end if
    end do

  end subroutine get_datasets

!===============================================================================
! GET_NAME Obtains the name of the current group in group_id
!===============================================================================

  function get_name(group_id, name_len_) result(name)
    integer(HID_T),            intent(in) :: group_id
    integer(SIZE_T), optional, intent(in) :: name_len_

    character(len=150) :: name ! name of group
    integer(SIZE_T) :: name_len, name_file_len
    integer :: hdf5_err ! HDF5 error code

    if (present(name_len_)) then
      name_len = name_len_
    else
      name_len = 150
    end if

    call h5iget_name_f(group_id, name, name_len, name_file_len, hdf5_err)
  end function get_name

!===============================================================================
! OPEN_GROUP opens an existing HDF5 group
!===============================================================================

  function open_group(group_id, name) result(newgroup_id)
    integer(HID_T), intent(in) :: group_id
    character(*),   intent(in) :: name ! name of group
    integer(HID_T) :: newgroup_id

    interface
      function open_group_c(group_id, name) result(newgroup_id) &
           bind(C, name='open_group')
        import HID_T, C_CHAR
        integer(HID_T), value :: group_id
        character(kind=C_CHAR), intent(in) :: name(*)
        integer(HID_T) :: newgroup_id
      end function open_group_c
    end interface

    newgroup_id = open_group_c(group_id, to_c_string(name))
  end function open_group

!===============================================================================
! CREATE_GROUP creates a new HDF5 group
!===============================================================================

  function create_group(group_id, name) result(newgroup_id)
    integer(HID_T), intent(in) :: group_id
    character(*),   intent(in) :: name ! name of group
    integer(HID_T) :: newgroup_id

    interface
      function create_group_c(group_id, name) result(newgroup_id) &
           bind(C, name='create_group')
        import HID_T, C_CHAR
        integer(HID_T), value :: group_id
        character(kind=C_CHAR), intent(in) :: name(*)
        integer(HID_T) :: newgroup_id
      end function create_group_c
    end interface

    newgroup_id = create_group_c(group_id, to_c_string(name))
  end function create_group

!===============================================================================
! CLOSE_GROUP closes HDF5 temp_group
!===============================================================================

  subroutine close_group(group_id)
    integer(HID_T), intent(inout) :: group_id
    interface
      subroutine close_group_c(group_id) bind(C, name='close_group')
        import HID_T
        integer(HID_T), value :: group_id
      end subroutine close_group_c
    end interface

    call close_group_c(group_id)
  end subroutine close_group

!===============================================================================
! OPEN_DATASET opens an existing HDF5 dataset
!===============================================================================

  function open_dataset(group_id, name) result(dataset_id)
    integer(HID_T), intent(in) :: group_id
    character(*),   intent(in) :: name ! name of dataset
    integer(HID_T)             :: dataset_id

    interface
      function open_dataset_c(group_id, name) result(dset_id) &
           bind(C, name='open_dataset')
        import HID_T, C_CHAR
        integer(HID_T), value :: group_id
        character(kind=C_CHAR), intent(in) :: name(*)
        integer(HID_T) :: dset_id
      end function open_dataset_c
    end interface

    dataset_id = open_dataset_c(group_id, to_c_string(name))
  end function open_dataset

!===============================================================================
! CLOSE_DATASET closes an HDF5 dataset
!===============================================================================

  subroutine close_dataset(dataset_id)
    integer(HID_T), intent(inout) :: dataset_id

    interface
      subroutine close_dataset_c(dset_id) bind(C, name='close_dataset')
        import HID_T
        integer(HID_T), value :: dset_id
      end subroutine close_dataset_c
    end interface

    call close_dataset_c(dataset_id)
  end subroutine close_dataset

!===============================================================================
! WRITE_DOUBLE_ND writes double precision N-D array data
!===============================================================================

  subroutine write_double_0D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    real(8),      intent(in), target   :: buffer  ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(0)
    logical(C_BOOL) :: indep_
    real(C_DOUBLE) :: buffer_(1)

    indep_ = .false.
    if (present(indep)) indep_ = indep
    buffer_(1) = buffer

    call write_double_c(group_id, 0, dims, to_c_string(name), buffer_, indep_)
  end subroutine write_double_0D

  subroutine write_double_1D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(1)
    logical(C_BOOL) :: indep_

    dims(1) = size(buffer)
    indep_ = .false.
    if (present(indep)) indep_ = indep

    call write_double_c(group_id, 1, dims, to_c_string(name), buffer, indep_)
  end subroutine write_double_1D

  subroutine write_double_2D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(2)
    logical(C_BOOL) :: indep_

    ! Reverse shape of array since it will be written from C
    dims(1) = size(buffer, 2)
    dims(2) = size(buffer, 1)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    call write_double_c(group_id, 2, dims, to_c_string(name), buffer, indep_)
  end subroutine write_double_2D

  subroutine write_double_3D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(3)
    logical(C_BOOL) :: indep_

    ! Reverse shape of array since it will be written from C
    dims(1) = size(buffer, 3)
    dims(2) = size(buffer, 2)
    dims(3) = size(buffer, 1)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    call write_double_c(group_id, 3, dims, to_c_string(name), buffer, indep_)
  end subroutine write_double_3D

  subroutine write_double_4D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    real(8),      intent(in), target   :: buffer(:,:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(4)
    logical(C_BOOL) :: indep_

    ! Reverse shape of array since it will be written from C
    dims(1) = size(buffer, 4)
    dims(2) = size(buffer, 3)
    dims(3) = size(buffer, 2)
    dims(4) = size(buffer, 1)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    call write_double_c(group_id, 4, dims, to_c_string(name), buffer, indep_)
  end subroutine write_double_4D

!===============================================================================
! READ_DOUBLE_ND reads double precision N-D array data
!===============================================================================

  subroutine read_double_0D(buffer, obj_id, name, indep)
    real(8),      target,   intent(inout) :: buffer
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    real(C_DOUBLE) :: buffer_(1)
    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_double_c(obj_id, c_loc(name_), buffer_, indep_)
    else
      call read_double_c(obj_id, C_NULL_PTR, buffer_, indep_)
    end if
    buffer = buffer_(1)
  end subroutine read_double_0D

  subroutine read_double_1D(buffer, obj_id, name, indep)
    real(8), target,        intent(inout) :: buffer(:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_double_c(obj_id, c_loc(name_), buffer, indep_)
    else
      call read_double_c(obj_id, C_NULL_PTR, buffer, indep_)
    end if
  end subroutine read_double_1D

  subroutine read_double_2D(buffer, obj_id, name, indep)
    real(8), target,        intent(inout) :: buffer(:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_double_c(obj_id, c_loc(name_), buffer, indep_)
    else
      call read_double_c(obj_id, C_NULL_PTR, buffer, indep_)
    end if
  end subroutine read_double_2D

  subroutine read_double_3D(buffer, obj_id, name, indep)
    real(8), target,        intent(inout) :: buffer(:,:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_double_c(obj_id, c_loc(name_), buffer, indep_)
    else
      call read_double_c(obj_id, C_NULL_PTR, buffer, indep_)
    end if
  end subroutine read_double_3D

  subroutine read_double_4D(buffer, obj_id, name, indep)
    real(8), target,        intent(inout) :: buffer(:,:,:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_double_c(obj_id, c_loc(name_), buffer, indep_)
    else
      call read_double_c(obj_id, C_NULL_PTR, buffer, indep_)
    end if
  end subroutine read_double_4D

!===============================================================================
! WRITE_INTEGER_ND writes integer precision N-D array data
!===============================================================================

  subroutine write_integer_0D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    integer,      intent(in), target   :: buffer  ! data to write
    logical,      intent(in), optional :: indep ! independent I/O

    integer(HSIZE_T) :: dims(0)
    logical(C_BOOL) :: indep_
    integer(C_INT) :: buffer_(1)

    indep_ = .false.
    if (present(indep)) indep_ = indep
    buffer_(1) = buffer

    call write_int_c(group_id, 0, dims, to_c_string(name), buffer_, indep_)
  end subroutine write_integer_0D

  subroutine write_integer_1D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(1)
    logical(C_BOOL) :: indep_

    dims(1) = size(buffer)
    indep_ = .false.
    if (present(indep)) indep_ = indep

    call write_int_c(group_id, 1, dims, to_c_string(name), buffer, indep_)
  end subroutine write_integer_1D

  subroutine write_integer_2D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(2)
    logical(C_BOOL) :: indep_

    ! Reverse shape of array since it will be written from C
    dims(1) = size(buffer, 2)
    dims(2) = size(buffer, 1)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    call write_int_c(group_id, 2, dims, to_c_string(name), buffer, indep_)
  end subroutine write_integer_2D

  subroutine write_integer_3D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(3)
    logical(C_BOOL) :: indep_

    ! Reverse shape of array since it will be written from C
    dims(1) = size(buffer, 3)
    dims(2) = size(buffer, 2)
    dims(3) = size(buffer, 1)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    call write_int_c(group_id, 3, dims, to_c_string(name), buffer, indep_)
  end subroutine write_integer_3D

  subroutine write_integer_4D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name      ! name of data
    integer,      intent(in), target   :: buffer(:,:,:,:) ! data to write
    logical,      intent(in), optional :: indep   ! independent I/O

    integer(HSIZE_T) :: dims(4)
    logical(C_BOOL) :: indep_

    ! Reverse shape of array since it will be written from C
    dims(1) = size(buffer, 4)
    dims(2) = size(buffer, 3)
    dims(3) = size(buffer, 2)
    dims(4) = size(buffer, 1)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    call write_int_c(group_id, 3, dims, to_c_string(name), buffer, indep_)
  end subroutine write_integer_4D

!===============================================================================
! READ_INTEGER_ND reads integer precision N-D array data
!===============================================================================

  subroutine read_integer_0D(buffer, obj_id, name, indep)
    integer,      target,   intent(inout) :: buffer
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer(C_INT) :: buffer_(1)
    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_int_c(obj_id, c_loc(name_), buffer_, indep_)
    else
      call read_int_c(obj_id, C_NULL_PTR, buffer_, indep_)
    end if
    buffer = buffer_(1)
  end subroutine read_integer_0D

  subroutine read_integer_1D(buffer, obj_id, name, indep)
    integer, target,        intent(inout) :: buffer(:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_int_c(obj_id, c_loc(name_), buffer, indep_)
    else
      call read_int_c(obj_id, C_NULL_PTR, buffer, indep_)
    end if
  end subroutine read_integer_1D

  subroutine read_integer_2D(buffer, obj_id, name, indep)
    integer, target,        intent(inout) :: buffer(:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_int_c(obj_id, c_loc(name_), buffer, indep_)
    else
      call read_int_c(obj_id, C_NULL_PTR, buffer, indep_)
    end if
  end subroutine read_integer_2D

  subroutine read_integer_3D(buffer, obj_id, name, indep)
    integer, target,        intent(inout) :: buffer(:,:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_int_c(obj_id, c_loc(name_), buffer, indep_)
    else
      call read_int_c(obj_id, C_NULL_PTR, buffer, indep_)
    end if
  end subroutine read_integer_3D

  subroutine read_integer_4D(buffer, obj_id, name, indep)
    integer, target,        intent(inout) :: buffer(:,:,:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_int_c(obj_id, c_loc(name_), buffer, indep_)
    else
      call read_int_c(obj_id, C_NULL_PTR, buffer, indep_)
    end if
  end subroutine read_integer_4D

!===============================================================================
! WRITE_LONG writes long integer scalar data
!===============================================================================

  subroutine write_long(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    integer(8),   intent(in), target   :: buffer  ! data to write
    logical,      intent(in), optional :: indep ! independent I/O

    integer(HSIZE_T) :: dims(0)
    logical(C_BOOL) :: indep_
    integer(C_LONG_LONG) :: buffer_(1)

    indep_ = .false.
    if (present(indep)) indep_ = indep
    buffer_(1) = buffer

    call write_llong_c(group_id, 0, dims, to_c_string(name), buffer_, indep_)
  end subroutine write_long

!===============================================================================
! READ_LONG reads long integer scalar data
!===============================================================================

  subroutine read_long(buffer, obj_id, name, indep)
    integer(8),   target,   intent(inout) :: buffer
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    integer(C_LONG_LONG) :: buffer_(1)
    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_llong_c(obj_id, c_loc(name_), buffer_, indep_)
    else
      call read_llong_c(obj_id, C_NULL_PTR, buffer_, indep_)
    end if
    buffer = buffer_(1)
  end subroutine read_long

!===============================================================================
! WRITE_STRING writes string data
!===============================================================================

  subroutine write_string(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), target   :: buffer  ! read data to here
    logical,      intent(in), optional :: indep ! independent I/O

    integer(HSIZE_T) :: dims(0)
    integer(C_SIZE_T) :: slen
    logical(C_BOOL) :: indep_

    indep_ = .false.
    if (present(indep)) indep_ = indep
    slen = len_trim(buffer)

    call write_string_c(group_id, 0, dims, slen, to_c_string(name), &
         to_c_string(buffer), indep_)
  end subroutine write_string

!===============================================================================
! READ_STRING reads string data
!===============================================================================

  subroutine read_string(buffer, obj_id, name, indep)
    character(*), target,   intent(inout) :: buffer
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep ! independent I/O

    integer(HID_T) :: dset_id
    integer(C_SIZE_T) :: i, n
    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), allocatable, target :: buffer_(:)

    if (present(name)) then
      dset_id = open_dataset(obj_id, name)
    else
      dset_id = obj_id
    end if

    ! Allocate a C char array to get string
    n = dataset_typesize(dset_id)
    allocate(buffer_(n))

    indep_ = .false.
    if (present(indep)) indep_ = indep
    call read_string_c(dset_id, C_NULL_PTR, n, buffer_, indep_)

    buffer = ''
    do i = 1, n
      if (buffer_(i) == C_NULL_CHAR) cycle
      buffer(i:i) = buffer_(i)
    end do

    call close_dataset(dset_id)
  end subroutine read_string

!===============================================================================
! WRITE_STRING_1D writes string 1-D array data
!===============================================================================

  subroutine write_string_1D(group_id, name, buffer, indep)
    integer(HID_T), intent(in) :: group_id
    character(*), intent(in)           :: name    ! name for data
    character(*), intent(in), target   :: buffer(:)  ! read data to here
    logical,      intent(in), optional :: indep ! independent I/O

    integer :: i
    integer(HSIZE_T) :: dims(1)
    integer(C_SIZE_T) :: m
    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), allocatable :: buffer_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! Copy array of characters into an array of C chars with the right memory
    ! layout
    dims(1) = size(buffer)
    m = maxval(len_trim(buffer)) + 1
    allocate(buffer_(dims(1)*m))
    do i = 0, dims(1) - 1
      buffer_(i*m+1 : (i+1)*m) = to_c_string(buffer(i+1))
    end do

    call write_string_c(group_id, 1, dims, m, to_c_string(name), &
         buffer_, indep_)
  end subroutine write_string_1D

!===============================================================================
! READ_STRING_1D reads string 1-D array data
!===============================================================================

  subroutine read_string_1D(buffer, obj_id, name, indep)
    character(*), target,   intent(inout) :: buffer(:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep ! independent I/O

    integer(HID_T) :: dset_id
    integer(C_SIZE_T) :: i, j, k, n, m
    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), allocatable, target :: buffer_(:)

    if (present(name)) then
      dset_id = open_dataset(obj_id, name)
    else
      dset_id = obj_id
    end if

    ! Allocate a C char array to get strings
    n = dataset_typesize(dset_id)
    m = size(buffer)
    allocate(buffer_(n*m))

    indep_ = .false.
    if (present(indep)) indep_ = indep
    call read_string_c(dset_id, C_NULL_PTR, n, buffer_, indep_)

    ! Convert null-terminated C strings into Fortran strings
    do i = 1, m
      buffer(i) = ''
      do j = 1, n
        k = (i-1)*n + j
        if (buffer_(k) == C_NULL_CHAR) exit
        buffer(i)(j:j) = buffer_(k)
      end do
    end do

    call close_dataset(dset_id)
  end subroutine read_string_1D

!===============================================================================
! WRITE_ATTRIBUTE_STRING
!===============================================================================

  subroutine write_attribute_string(obj_id, name, buffer)
    integer(HID_T), intent(in)       :: obj_id    ! object to write attribute to
    character(*), intent(in)         :: name      ! name of attribute
    character(*), intent(in), target :: buffer    ! string to write


    call write_attr_string_c(obj_id, to_c_string(name), to_c_string(buffer))
  end subroutine write_attribute_string

  subroutine read_attribute_double(buffer, obj_id, name)
    real(8),        intent(inout), target :: buffer
    integer(HID_T), intent(in)            :: obj_id
    character(*),   intent(in)            :: name

    real(C_DOUBLE) :: buffer_(1)

    call read_attr_double_c(obj_id, to_c_string(name), buffer_)
    buffer = buffer_(1)
  end subroutine read_attribute_double

  subroutine write_attribute_double(obj_id, name, buffer)
    integer(HID_T), intent(in)  :: obj_id
    character(*),   intent(in)  :: name
    real(8), intent(in), target :: buffer

    integer(HSIZE_T) :: dims(0)
    real(C_DOUBLE) :: buffer_(1)

    buffer_(1) = buffer
    call write_attr_double_c(obj_id, 0, dims, to_c_string(name), buffer_)
  end subroutine write_attribute_double

  subroutine read_attribute_double_1D(buffer, obj_id, name)
    real(8), target, allocatable, intent(inout) :: buffer(:)
    integer(HID_T),  intent(in)    :: obj_id
    character(*),    intent(in)    :: name

    integer(HSIZE_T) :: dims(1)

    if (.not. allocated(buffer)) then
      call get_shape_attr(obj_id, to_c_string(name), dims)
      allocate(buffer(dims(1)))
    end if

    call read_attr_double_c(obj_id, to_c_string(name), buffer)
  end subroutine read_attribute_double_1D

  subroutine write_attribute_double_1D(obj_id, name, buffer)
    integer(HID_T),  intent(in) :: obj_id
    character(*),    intent(in) :: name
    real(8), target, intent(in) :: buffer(:)

    integer(HSIZE_T) :: dims(1)

    dims(1) = size(buffer)
    call write_attr_double_c(obj_id, 1, dims, to_c_string(name), buffer)
  end subroutine write_attribute_double_1D

  subroutine read_attribute_double_2D(buffer, obj_id, name)
    real(8), target, allocatable, intent(inout) :: buffer(:,:)
    integer(HID_T),  intent(in)    :: obj_id
    character(*),    intent(in)    :: name

    integer(HSIZE_T) :: dims(2)

    if (.not. allocated(buffer)) then
      call get_shape_attr(obj_id, to_c_string(name), dims)
      allocate(buffer(dims(2), dims(1)))
    end if

    call read_attr_double_c(obj_id, to_c_string(name), buffer)
  end subroutine read_attribute_double_2D

  subroutine read_attribute_integer(buffer, obj_id, name)
    integer,        intent(inout), target :: buffer
    integer(HID_T), intent(in)            :: obj_id
    character(*),   intent(in)            :: name

    integer(C_INT) :: buffer_(1)

    call read_attr_int_c(obj_id, to_c_string(name), buffer_)
    buffer = buffer_(1)
  end subroutine read_attribute_integer

  subroutine write_attribute_integer(obj_id, name, buffer)
    integer(HID_T), intent(in)  :: obj_id
    character(*),   intent(in)  :: name
    integer, intent(in), target :: buffer

    integer(HSIZE_T) :: dims(0)
    integer(C_INT) :: buffer_(1)

    buffer_(1) = buffer
    call write_attr_int_c(obj_id, 0, dims, to_c_string(name), buffer_)
  end subroutine write_attribute_integer

  subroutine read_attribute_integer_1D(buffer, obj_id, name)
    integer, target, allocatable, intent(inout) :: buffer(:)
    integer(HID_T),  intent(in)    :: obj_id
    character(*),    intent(in)    :: name

    integer(HSIZE_T) :: dims(1)

    if (.not. allocated(buffer)) then
      call get_shape_attr(obj_id, to_c_string(name), dims)
      allocate(buffer(dims(1)))
    end if

    call read_attr_int_c(obj_id, to_c_string(name), buffer)
  end subroutine read_attribute_integer_1D

  subroutine write_attribute_integer_1D(obj_id, name, buffer)
    integer(HID_T),  intent(in) :: obj_id
    character(*),    intent(in) :: name
    integer, target, intent(in) :: buffer(:)

    integer(HSIZE_T) :: dims(1)

    dims(1) = size(buffer)
    call write_attr_int_c(obj_id, 1, dims, to_c_string(name), buffer)
  end subroutine write_attribute_integer_1D

  subroutine read_attribute_integer_2D(buffer, obj_id, name)
    integer, target, allocatable, intent(inout) :: buffer(:,:)
    integer(HID_T),  intent(in)    :: obj_id
    character(*),    intent(in)    :: name

    integer(HSIZE_T) :: dims(2)

    if (.not. allocated(buffer)) then
      call get_shape_attr(obj_id, to_c_string(name), dims)
      allocate(buffer(dims(2), dims(1)))
    end if

    call read_attr_int_c(obj_id, to_c_string(name), buffer)
  end subroutine read_attribute_integer_2D

  subroutine read_attribute_string(buffer, obj_id, name)
    character(*),   intent(inout) :: buffer  ! read data to here
    integer(HID_T), intent(in)    :: obj_id
    character(*),   intent(in)    :: name    ! name for data

    integer(C_SIZE_T) :: i, n
    character(kind=C_CHAR), allocatable, target :: buffer_(:)

    ! Allocate a C char array to get string
    n = attribute_typesize(obj_id, to_c_string(name))
    allocate(buffer_(n))

    ! Read attribute
    call read_attr_string_c(obj_id, to_c_string(name), n, buffer_)

    ! Copy back to Fortran string
    buffer = ''
    do i = 1, n
      if (buffer_(i) == C_NULL_CHAR) cycle
      buffer(i:i) = buffer_(i)
    end do
  end subroutine read_attribute_string

  subroutine read_attribute_string_1D(buffer, obj_id, name)
    character(*), target, allocatable, intent(inout) :: buffer(:)
    integer(HID_T), intent(in) :: obj_id
    character(*),   intent(in) :: name

    integer(C_SIZE_T) :: i, j, k, n, m
    integer(HSIZE_T) :: dims(1)
    character(kind=C_CHAR), allocatable, target :: buffer_(:)

    if (.not. allocated(buffer)) then
      call get_shape_attr(obj_id, to_c_string(name), dims)
      allocate(buffer(dims(1)))
    end if

    ! Allocate a C char array to get strings
    n = attribute_typesize(obj_id, to_c_string(name))
    m = size(buffer)
    allocate(buffer_((n+1)*m))

    ! Read attribute
    call read_attr_string_c(obj_id, to_c_string(name), n, buffer_)

    ! Convert null-terminated C strings into Fortran strings
    do i = 1, m
      buffer(i) = ''
      do j = 1, n
        k = (i-1)*(n+1) + j
        if (buffer_(k) == C_NULL_CHAR) exit
        buffer(i)(j:j) = buffer_(k)
      end do
    end do
  end subroutine read_attribute_string_1D

  subroutine read_attribute_logical(buffer, obj_id, name)
    logical,        intent(inout), target :: buffer
    integer(HID_T), intent(in)            :: obj_id
    character(*),   intent(in)            :: name

    integer :: tmp

    call read_attribute_integer(tmp, obj_id, name)
    buffer = (tmp /= 0)
  end subroutine read_attribute_logical

  subroutine get_shape(obj_id, dims)
    integer(HID_T),   intent(in)  :: obj_id
    integer(HSIZE_T), intent(out) :: dims(:)

    integer :: i
    integer(HSIZE_T) :: dims_c(size(dims))

    call get_shape_c(obj_id, dims_c)
    do i = 1, size(dims)
      dims(i) = dims_c(size(dims) - i + 1)
    end do
  end subroutine get_shape

  subroutine get_ndims(obj_id, ndims)
    integer(HID_T), intent(in)  :: obj_id
    integer,        intent(out) :: ndims

    interface
      function dataset_ndims(dset) result(ndims) bind(C)
        import HID_T, C_INT
        integer(HID_T), value :: dset
        integer(C_INT) :: ndims
      end function dataset_ndims
    end interface

    ndims = dataset_ndims(obj_id)
  end subroutine get_ndims

!===============================================================================
! READ_COMPLEX_2D reads double precision complex 2-D array data as output by
! the h5py HDF5 python module.
!===============================================================================

  subroutine read_complex_2D(buffer, obj_id, name, indep)
    complex(C_DOUBLE_COMPLEX), target, intent(inout) :: buffer(:,:)
    integer(HID_T),         intent(in)    :: obj_id
    character(*), optional, intent(in)    :: name
    logical,      optional, intent(in)    :: indep  ! independent I/O

    logical(C_BOOL) :: indep_
    character(kind=C_CHAR), target, allocatable :: name_(:)

    indep_ = .false.
    if (present(indep)) indep_ = indep

    ! If 'name' argument is passed, obj_id is interpreted to be a group and
    ! 'name' is the name of the dataset we should read from
    if (present(name)) then
      name_ = to_c_string(name)
      call read_complex_c(obj_id, c_loc(name_), buffer, indep_)
    else
      call read_complex_c(obj_id, C_NULL_PTR, buffer, indep_)
    end if
  end subroutine read_complex_2D

end module hdf5_interface
