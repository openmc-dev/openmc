module material_header

  implicit none

!===============================================================================
! COMPOSITION describes a material by its constituent nuclide atom fractions
!===============================================================================

  type Composition
    real(8),allocatable  :: atom_density(:)

  end type Composition


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
    integer, allocatable           :: nuclide(:) ! index in nuclides array
    type(Density)                  :: density    ! atom densities in atom/b-cm
    type(Composition), allocatable :: comp(:)    ! atom fractions via compositions

    ! S(a,b) data references
    integer              :: n_sab = 0         ! number of S(a,b) tables
    integer, allocatable :: i_sab_nuclides(:) ! index of corresponding nuclide
    integer, allocatable :: i_sab_tables(:)   ! index in sab_tables

    ! Temporary names read during initialization
    character(12), allocatable :: names(:)     ! isotope names
    character(12), allocatable :: sab_names(:) ! name of S(a,b) table

    ! Distribution Variables
    logical                        :: distrib_dens ! distributed densities
    logical                        :: distrib_comp ! distributed compositions
    integer                        :: comp_map     ! composition map ID
    integer                        :: dens_map     ! density map ID
  contains
    procedure :: get_density => get_density
  end type Material



contains


  function get_density(this, i, j) result(density)

    class(Material)     :: this
    integer, intent(in) :: i                ! i_th composition
    integer, intent(in) :: j                ! j_th nuclide in the composition
    real(8)             :: density          ! density to be returned

    density = this % comp(i) % atom_density(j)

  end function get_density
    
end module material_header
