module material_header

  use constants
  use nuclide_header, only: Nuclide

  implicit none

!===============================================================================
! MATERIAL describes a material by its constituent nuclides
!===============================================================================

  type Material
    integer              :: id              ! unique identifier
    character(len=104)   :: name = ""       ! User-defined name
    integer              :: n_nuclides      ! number of nuclides
    integer, allocatable :: nuclide(:)      ! index in nuclides array
    integer, allocatable :: element(:)      ! index in elements array
    real(8)              :: density         ! total atom density in atom/b-cm
    real(8), allocatable :: atom_density(:) ! nuclide atom density in atom/b-cm
    real(8)              :: density_gpcc    ! total density in g/cm^3

    ! Energy grid information
    integer              :: n_grid    ! # of union material grid points
    real(8), allocatable :: e_grid(:) ! union material grid energies

    ! Unionized energy grid information
    integer, allocatable :: nuclide_grid_index(:,:) ! nuclide e_grid pointers

    ! S(a,b) data
    integer              :: n_sab = 0         ! number of S(a,b) tables
    integer, allocatable :: i_sab_nuclides(:) ! index of corresponding nuclide
    integer, allocatable :: i_sab_tables(:)   ! index in sab_tables
    real(8), allocatable :: sab_fracs(:)      ! how often to use S(a,b)

    ! Temporary names read during initialization
    character(20), allocatable :: names(:)     ! isotope names
    character(20), allocatable :: sab_names(:) ! name of S(a,b) table

    ! Does this material contain fissionable nuclides? Is it depletable?
    logical :: fissionable = .false.
    logical :: depletable = .false.

    ! enforce isotropic scattering in lab
    logical, allocatable :: p0(:)

  contains
    procedure :: set_density => material_set_density
  end type Material

contains

!===============================================================================
! MATERIAL_SET_DENSITY sets the total density of a material in atom/b-cm.
!===============================================================================

  function material_set_density(m, density, nuclides) result(err)
    class(Material), intent(inout) :: m
    real(8), intent(in) :: density
    type(Nuclide), intent(in) :: nuclides(:)
    integer :: err

    integer :: i
    real(8) :: sum_percent
    real(8) :: awr

    err = -1
    if (allocated(m % atom_density)) then
      ! Set total density based on value provided
      m % density = density

      ! Determine normalized atom percents
      sum_percent = sum(m % atom_density)
      m % atom_density(:) = m % atom_density / sum_percent

      ! Recalculate nuclide atom densities based on given density
      m % atom_density(:) = density * m % atom_density

      ! Calculate density in g/cm^3.
      m % density_gpcc = ZERO
      do i = 1, m % n_nuclides
        awr = nuclides(m % nuclide(i)) % awr
        m % density_gpcc = m % density_gpcc &
             + m % atom_density(i) * awr * MASS_NEUTRON / N_AVOGADRO
      end do
      err = 0
    end if
  end function material_set_density

end module material_header
