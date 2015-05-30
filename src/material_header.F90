module material_header

  use constants,        only: OTF_HEADROOM, MAX_LINE_LEN
  use dict_header,      only: DictIntInt
  use output_interface, only: BinaryOutput

#ifdef HDF5
  use hdf5
#endif

  implicit none

#ifdef HDF5
    integer :: hdf5_err
#endif

!===============================================================================
! COMPOSITION describes a material by its constituent nuclide atom fractions
!===============================================================================

  type Composition
    real(8), allocatable  :: atom_density(:)
  end type Composition

!===============================================================================
! COMPOSITIONFILE is a container for a composition file for distributed mats
!===============================================================================

  type CompositionFile
    character(MAX_LINE_LEN) :: path         ! path to the file
    character(MAX_LINE_LEN) :: group        ! group in HDF5 file
    type(BinaryOutput)      :: fh           ! file handle
    integer                 :: n_nuclides   ! number of comps per row
    integer                 :: n_instances  ! number of comp rows
    logical                 :: initialized = .false.
#ifdef HDF5
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer(HID_T) :: dset       ! data set handle
    integer(HID_T) :: dspace     ! data or file space handle
    integer(HID_T) :: memspace   ! data space handle for individual procs
    integer(HID_T) :: plist      ! property list handle
    integer(HSIZE_T) :: block1(1)
    integer(HSIZE_T) :: dims1(1)
#endif
  contains
    procedure :: init => composition_file_init
    procedure :: load => composition_file_load
    procedure :: close => composition_file_close
  end type CompositionFile

  integer, parameter :: COMPFILE_TYPE_BINARY = 1
  integer, parameter :: COMPFILE_TYPE_HDF5 = 2

!===============================================================================
! DENSITY describes a material by its constituent nuclide densities
!===============================================================================

  type Density
    real(8), allocatable :: density(:)       ! nuclide densities in X/Y
    integer              :: num              ! size of density
  end type Density

!===============================================================================
! MATERIAL describes a material by its constituent nuclides
!===============================================================================

  type Material
    integer                        :: id         ! unique identifier
    character(len=52)              :: name = ""  ! User-defined name
    integer                        :: n_nuclides ! number of nuclides
    integer                        :: n_comp     ! number of compositions
    integer                        :: cell       ! assigned cell 
                                                 ! only for distributed material
    integer                        :: map        ! map number for this material
    integer, allocatable           :: nuclide(:) ! index in nuclides array
    type(Density)                  :: density    ! material density in atom/b-cm
    type(Composition), allocatable :: comp(:)    ! atom densities in atom/b-cm

    ! Energy grid information
    integer              :: n_grid    ! # of union material grid points
    real(8), allocatable :: e_grid(:) ! union material grid energies

    ! Unionized energy grid information
    integer, allocatable :: nuclide_grid_index(:,:) ! nuclide e_grid pointers

    ! S(a,b) data references
    integer              :: n_sab = 0         ! number of S(a,b) tables
    integer, allocatable :: i_sab_nuclides(:) ! index of corresponding nuclide
    integer, allocatable :: i_sab_tables(:)   ! index in sab_tables

    ! Temporary names read during initialization
    character(12), allocatable :: names(:)     ! isotope names
    character(12), allocatable :: sab_names(:) ! name of S(a,b) table

    ! Distribution variables
    logical                        :: distrib_dens ! distributed densities
    logical                        :: distrib_comp ! distributed compositions

    ! On-the-fly allocation controls
    real(8), allocatable           :: otf_comp(:, :)
    logical :: otf_compositions = .false.
    type(CompositionFile) :: comp_file          ! compositions file
    integer :: size_comp_array                  ! Size of composition array
                                                ! when variable
    integer :: next_comp_idx = 1                ! Next index in the composition
                                                ! array when variable
    type(DictIntInt) :: comp_index_map          ! real_inst --> local_inst
    type(DictIntInt) :: reverse_comp_index_map  ! local_inst --> real_inst

    ! Does this material contain fissionable nuclides?
    logical :: fissionable = .false.

  contains
    procedure :: get_density
    procedure :: otf_comp_index
    procedure :: grow_composition_array

  end type Material
  
contains

!===============================================================================
! GET_DENSITY returns the atom density for a given instance of this material,
! which could be different inside many different fill-cell instances
!===============================================================================

  function get_density(this, i, j) result(density)

    class(Material), intent(inout) :: this
    integer,         intent(in)    :: i       ! i_th composition
    integer,         intent(in)    :: j       ! j_th nuclide in the composition
    real(8)                        :: density ! density to be returned

    if (this % otf_compositions) then
      density = this % otf_comp(j, this % otf_comp_index(i))
    else
      density = this % comp(i) % atom_density(j)
    end if

  end function get_density

!===============================================================================
! OTF_FILTER_INDEX returns the filter index in the results array when OTF tally
!===============================================================================

    function otf_comp_index(this, real_inst) result(idx)
      
      class(Material), intent(inout) :: this 
      integer,         intent(in)    :: real_inst
      
      integer :: idx

      ! If we're doing on-the-fly memory allocation, we use the map
      if (this % comp_index_map % has_key(real_inst)) then

        idx = this % comp_index_map % get_key(real_inst)

      else

        ! This is the first time this composition index has been needed

        ! Grow the composition array if we've used it all
        if (this % next_comp_idx > this % size_comp_array) &
            call this % grow_composition_array()

        ! Update the map
        call this % comp_index_map % add_key(real_inst, this % next_comp_idx)
        call this % reverse_comp_index_map % add_key( &
            this % next_comp_idx, real_inst)

        ! Set the return index
        idx = this % next_comp_idx
        
        ! Load the compositions for this instance
        this % otf_comp(:, idx) = this % comp_file % load(real_inst)

        ! Increment the next index
        this % next_comp_idx = this % next_comp_idx + 1

      end if

    end function otf_comp_index

!===============================================================================
! GROW_COMPOSITION_ARRAY
!===============================================================================

    subroutine grow_composition_array(this)

      class(Material), intent(inout) :: this 

      integer :: newsize
      real(8), allocatable :: temp(:, :)

      newsize = ceiling(real(this % size_comp_array, 8) * OTF_HEADROOM)

      ! Allocate results array with increased size
      allocate(temp(this % n_nuclides, newsize))

      ! Copy original results to temporary array
      temp(:, 1:this % size_comp_array) = &
          this % otf_comp(:, 1:this % size_comp_array)

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=this % otf_comp)

      ! Update size
      this % size_comp_array = newsize

    end subroutine grow_composition_array

!===============================================================================
! COMPOSITION_FILE_INIT
!===============================================================================

    subroutine composition_file_init(this, length, file_id, plist)

      class(CompositionFile), intent(inout) :: this 
      integer, intent(in) :: length
#ifndef HDF5
      integer, intent(in) :: file_id
      integer, intent(in) :: plist
#else
      integer(HID_T), intent(in) :: file_id
      integer(HID_T), intent(in) :: plist

      this % file_id = file_id
      this % plist = plist
      this % block1 = length
      this % dims1 = length

      ! Open the group
      call h5gopen_f(this % file_id, trim(this % group), this % group_id, hdf5_err)
  
      ! Open the dataset
      call h5dopen_f(this % group_id, 'comps', this % dset, hdf5_err)

      ! Open the dataspace and memory space
      call h5dget_space_f(this % dset, this % dspace, hdf5_err)
      call h5screate_simple_f(1, this % block1, this % memspace, hdf5_err)

      this % initialized = .true.

#endif

    end subroutine composition_file_init

!===============================================================================
! COMPOSITION_FILE_LOAD
!===============================================================================

    function composition_file_load(this, real_inst) result(comp)
      
      class(CompositionFile), intent(inout) :: this 
      integer,                intent(in)    :: real_inst
 
      real(8), allocatable :: comp(:)
     
#ifdef HDF5
      integer(HSIZE_T) :: start1(1)  ! start type for 1-D array
      integer(HSIZE_T) :: count1(1)  ! count type for 1-D array

      if (.not. this % initialized) then
        print *, 'Trying to load OTF materials with an uninitialized file!'
        stop
      end if

      allocate(comp(this % block1(1)))

      start1 = (real_inst - 1) * this % block1
      count1 = 1

      ! Select the hyperslab
      call h5sselect_hyperslab_f(this % dspace, H5S_SELECT_SET_F, start1, &
          count1, hdf5_err, block = this % block1)

      ! Read data from file into memory
      call H5dread_f(this % dset, H5T_NATIVE_DOUBLE, comp, &
          this % dims1, hdf5_err, xfer_prp = this % plist, &
          mem_space_id = this % memspace, file_space_id = this % dspace)

#endif

    end function composition_file_load

!===============================================================================
! COMPOSITION_FILE_CLOSE
!===============================================================================

    subroutine composition_file_close(this)

      class(CompositionFile), intent(inout) :: this

#ifdef HDF5
      ! Close the dataspace and memory space
      call h5sclose_f(this % dspace, hdf5_err)
      call h5sclose_f(this % memspace, hdf5_err)

      ! Close dataset and group
      call h5dclose_f(this % dset, hdf5_err)
      call h5gclose_f(this % group_id, hdf5_err)
#endif

    end subroutine composition_file_close

end module material_header
