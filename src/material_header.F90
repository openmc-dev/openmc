module material_header

  implicit none

!===============================================================================
! MATERIAL describes a material by its constituent nuclides
!===============================================================================

  type Material
    integer              :: id              ! unique identifier
    character(len=104) :: name = ""         ! User-defined name
    integer              :: n_nuclides      ! number of nuclides
    integer, allocatable :: nuclide(:)      ! index in nuclides array
    real(8)              :: density         ! total atom density in atom/b-cm
    real(8), allocatable :: atom_density(:) ! nuclide atom density in atom/b-cm

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

    ! Does this material contain fissionable nuclides?
    logical :: fissionable = .false.

    ! enforce isotropic scattering in lab
    logical, allocatable :: p0(:)

  end type Material

end module material_header
