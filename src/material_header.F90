module material_header

  use constants
  use error, only: fatal_error
  use nuclide_header, only: Nuclide
  use sab_header, only: SAlphaBeta
  use stl_vector, only: VectorReal, VectorInt
  use string, only: to_str

  implicit none

!===============================================================================
! MATERIAL describes a material by its constituent nuclides
!===============================================================================

  type Material
    integer              :: id              ! unique identifier
    character(len=104)   :: name = ""       ! User-defined name
    integer              :: n_nuclides      ! number of nuclides
    integer, allocatable :: nuclide(:)      ! index in nuclides array
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
    procedure :: assign_sab_tables => material_assign_sab_tables
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

!===============================================================================
! ASSIGN_SAB_TABLES assigns S(alpha,beta) tables to specific nuclides within
! materials so the code knows when to apply bound thermal scattering data
!===============================================================================

  subroutine material_assign_sab_tables(this, nuclides, sab_tables)
    class(Material), intent(inout) :: this
    type(Nuclide), intent(in) :: nuclides(:)
    type(SAlphaBeta), intent(in) :: sab_tables(:)

    integer :: j            ! index over nuclides in material
    integer :: k            ! index over S(a,b) tables in material
    integer :: m            ! position for sorting
    integer :: temp_nuclide ! temporary value for sorting
    integer :: temp_table   ! temporary value for sorting
    real(8) :: temp_frac    ! temporary value for sorting
    logical :: found
    type(VectorInt)  :: i_sab_tables
    type(VectorInt)  :: i_sab_nuclides
    type(VectorReal) :: sab_fracs

    if (.not. allocated(this % i_sab_tables)) return

    ASSIGN_SAB: do k = 1, size(this % i_sab_tables)
      ! In order to know which nuclide the S(a,b) table applies to, we need
      ! to search through the list of nuclides for one which has a matching
      ! name
      found = .false.
      associate (sab => sab_tables(this % i_sab_tables(k)))
        FIND_NUCLIDE: do j = 1, size(this % nuclide)
          if (any(sab % nuclides == nuclides(this % nuclide(j)) % name)) then
            call i_sab_tables % push_back(this % i_sab_tables(k))
            call i_sab_nuclides % push_back(j)
            call sab_fracs % push_back(this % sab_fracs(k))
            found = .true.
          end if
        end do FIND_NUCLIDE
      end associate

      ! Check to make sure S(a,b) table matched a nuclide
      if (.not. found) then
        call fatal_error("S(a,b) table " // trim(this % &
             sab_names(k)) // " did not match any nuclide on material " &
             // trim(to_str(this % id)))
      end if
    end do ASSIGN_SAB

    ! Make sure each nuclide only appears in one table.
    do j = 1, i_sab_nuclides % size()
      do k = j+1, i_sab_nuclides % size()
        if (i_sab_nuclides % data(j) == i_sab_nuclides % data(k)) then
          call fatal_error(trim( &
               nuclides(this % nuclide(i_sab_nuclides % data(j))) % name) &
               // " in material " // trim(to_str(this % id)) // " was found &
               &in multiple S(a,b) tables. Each nuclide can only appear in &
               &one S(a,b) table per material.")
        end if
      end do
    end do

    ! Update i_sab_tables and i_sab_nuclides
    deallocate(this % i_sab_tables)
    deallocate(this % sab_fracs)
    if (allocated(this % i_sab_nuclides)) deallocate(this % i_sab_nuclides)
    m = i_sab_tables % size()
    allocate(this % i_sab_tables(m))
    allocate(this % i_sab_nuclides(m))
    allocate(this % sab_fracs(m))
    this % i_sab_tables(:) = i_sab_tables % data(1:m)
    this % i_sab_nuclides(:) = i_sab_nuclides % data(1:m)
    this % sab_fracs(:) = sab_fracs % data(1:m)

    ! Clear entries in vectors for next material
    call i_sab_tables % clear()
    call i_sab_nuclides % clear()
    call sab_fracs % clear()

    ! If there are multiple S(a,b) tables, we need to make sure that the
    ! entries in i_sab_nuclides are sorted or else they won't be applied
    ! correctly in the cross_section module. The algorithm here is a simple
    ! insertion sort -- don't need anything fancy!

    if (size(this % i_sab_tables) > 1) then
      SORT_SAB: do k = 2, size(this % i_sab_tables)
        ! Save value to move
        m = k
        temp_nuclide = this % i_sab_nuclides(k)
        temp_table   = this % i_sab_tables(k)
        temp_frac    = this % i_sab_tables(k)

        MOVE_OVER: do
          ! Check if insertion value is greater than (m-1)th value
          if (temp_nuclide >= this % i_sab_nuclides(m-1)) exit

          ! Move values over until hitting one that's not larger
          this % i_sab_nuclides(m) = this % i_sab_nuclides(m-1)
          this % i_sab_tables(m)   = this % i_sab_tables(m-1)
          this % sab_fracs(m)      = this % sab_fracs(m-1)
          m = m - 1

          ! Exit if we've reached the beginning of the list
          if (m == 1) exit
        end do MOVE_OVER

        ! Put the original value into its new position
        this % i_sab_nuclides(m) = temp_nuclide
        this % i_sab_tables(m)   = temp_table
        this % sab_fracs(m)      = temp_frac
      end do SORT_SAB
    end if

    ! Deallocate temporary arrays for names of nuclides and S(a,b) tables
    if (allocated(this % names)) deallocate(this % names)
  end subroutine material_assign_sab_tables

end module material_header
