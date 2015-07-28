module material_header

  use list_header, only: ListChar

  implicit none

!===============================================================================
! MATERIAL describes a material by its constituent nuclides
!===============================================================================

  type Material
    integer              :: id              ! unique identifier
    integer              :: n_nuclides      ! number of nuclides
    integer, allocatable :: nuclide(:)      ! index in nuclides array
    real(8)              :: density         ! total atom density in atom/b-cm
    real(8), allocatable :: atom_density(:) ! nuclide atom density in atom/b-cm

    ! S(a,b) data references
    integer              :: n_sab = 0         ! number of S(a,b) tables
    integer, allocatable :: i_sab_nuclides(:) ! index of corresponding nuclide
    integer, allocatable :: i_sab_tables(:)   ! index in sab_tables

    ! Temporary names read during initialization
    character(12), allocatable :: names(:)     ! isotope names
    character(12), allocatable :: sab_names(:) ! name of S(a,b) table

    ! Does this material contain fissionable nuclides?
    logical :: fissionable = .false.

    ! Are any nuclides in this material given in terms of natural elements
    logical :: nat_elements = .false.

    ! String array of energy-xs data that need to be written
    type(ListChar), allocatable :: write_xs(:)

    ! String array of MT's for sec. angular dist's that need to be written
    type(ListChar), allocatable :: write_angle(:)

  end type Material

end module material_header
