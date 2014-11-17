module material_header

  use constants,        only: OTF_HEADROOM, MAX_LINE_LEN
  use dict_header,      only: DictIntInt
  use output_interface, only: BinaryOutput

  implicit none

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
    integer                 :: n_nuclides   ! number of comps per row
    integer                 :: n_instances  ! number of comp rows
  contains
    procedure :: load => composition_file_load
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
    integer                        :: n_nuclides ! number of nuclides
    integer                        :: n_comp     ! number of compositions
    integer                        :: cell       ! assigned cell 
                                                 ! only for distributed material
    integer                        :: map        ! map number for this material
    integer, allocatable           :: nuclide(:) ! index in nuclides array
    type(Density)                  :: density    ! material density in atom/b-cm
    type(Composition), allocatable :: comp(:)    ! atom densities in atom/b-cm

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
    logical :: otf_compositions = .false.
    type(CompositionFile) :: comp_file         ! compositions file
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
      density = this % comp(this % otf_comp_index(i)) % atom_density(j)
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

      ! If we're doing on-the-fly memory allocation, we must use the map
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
        this % comp(idx) = this % comp_file % load(real_inst)

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
      type(Composition), allocatable :: temp(:)

      newsize = ceiling(real(this % size_comp_array, 8) * OTF_HEADROOM)

      ! Allocate results array with increased size
      allocate(temp(newsize))

      ! Copy original results to temporary array
      temp(1:this % size_comp_array) = this % comp(1:this % size_comp_array)

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=this % comp)

      ! Update size
      this % size_comp_array = newsize

    end subroutine grow_composition_array

!===============================================================================
! COMPOSITION_FILE_LOAD
!===============================================================================

    function composition_file_load(this, real_inst) result(comp)
      
      class(CompositionFile), intent(inout) :: this 
      integer,                intent(in)    :: real_inst
      
      type(Composition) :: comp
      type(BinaryOutput) :: fh

      allocate(comp % atom_density(this % n_nuclides))

      call fh % file_open(this % path, 'r', serial = .false., &
          direct_access = .true., record_len = 8 * this % n_nuclides)
      call fh % read_data(comp % atom_density, 'comps', &
          length = this % n_nuclides, record = real_inst, offset = 16)
      call fh % file_close()

    end function composition_file_load

end module material_header
