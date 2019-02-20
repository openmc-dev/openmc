module nuclide_header

  use, intrinsic :: ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

  use constants

  implicit none

!===============================================================================
! NUCLIDEMICROXS contains cached microscopic cross sections for a particular
! nuclide at the current energy
!===============================================================================

  ! Arbitrary value to indicate invalid cache state for elastic scattering
  ! (NuclideMicroXS % elastic)
  real(8), parameter :: CACHE_INVALID = -1

  type, bind(C) :: NuclideMicroXS
    ! Microscopic cross sections in barns
    real(C_DOUBLE) :: total
    real(C_DOUBLE) :: absorption       ! absorption (disappearance)
    real(C_DOUBLE) :: fission          ! fission
    real(C_DOUBLE) :: nu_fission       ! neutron production from fission

    real(C_DOUBLE) :: elastic          ! If sab_frac is not 1 or 0, then this value is
                                !   averaged over bound and non-bound nuclei
    real(C_DOUBLE) :: thermal          ! Bound thermal elastic & inelastic scattering
    real(C_DOUBLE) :: thermal_elastic  ! Bound thermal elastic scattering
    real(C_DOUBLE) :: photon_prod      ! microscopic photon production xs

    ! Cross sections for depletion reactions (note that these are not stored in
    ! macroscopic cache)
    real(C_DOUBLE) :: reaction(6)

    ! Indicies and factors needed to compute cross sections from the data tables
    integer(C_INT) :: index_grid        ! Index on nuclide energy grid
    integer(C_INT) :: index_temp        ! Temperature index for nuclide
    real(C_DOUBLE) :: interp_factor     ! Interpolation factor on nuc. energy grid
    integer(C_INT) :: index_sab = NONE  ! Index in sab_tables
    integer(C_INT) :: index_temp_sab    ! Temperature index for sab_tables
    real(C_DOUBLE) :: sab_frac          ! Fraction of atoms affected by S(a,b)
    logical(C_BOOL) :: use_ptable        ! In URR range with probability tables?

    ! Energy and temperature last used to evaluate these cross sections.  If
    ! these values have changed, then the cross sections must be re-evaluated.
    real(C_DOUBLE) :: last_E = ZERO       ! Last evaluated energy
    real(C_DOUBLE) :: last_sqrtkT = ZERO  ! Last temperature in sqrt(Boltzmann
                                   !   constant * temperature (eV))
  end type NuclideMicroXS

!===============================================================================
! MATERIALMACROXS contains cached macroscopic cross sections for the material a
! particle is traveling through
!===============================================================================

  type, bind(C) :: MaterialMacroXS
    real(C_DOUBLE) :: total         ! macroscopic total xs
    real(C_DOUBLE) :: absorption    ! macroscopic absorption xs
    real(C_DOUBLE) :: fission       ! macroscopic fission xs
    real(C_DOUBLE) :: nu_fission    ! macroscopic production xs
    real(C_DOUBLE) :: photon_prod   ! macroscopic photon production xs

    ! Photon cross sections
    real(C_DOUBLE) :: coherent        ! macroscopic coherent xs
    real(C_DOUBLE) :: incoherent      ! macroscopic incoherent xs
    real(C_DOUBLE) :: photoelectric   ! macroscopic photoelectric xs
    real(C_DOUBLE) :: pair_production ! macroscopic pair production xs
  end type MaterialMacroXS

  ! Cross section caches
  type(NuclideMicroXS), allocatable, target :: micro_xs(:)  ! Cache for each nuclide
!$omp threadprivate(micro_xs)

  interface
    function nuclides_size() bind(C) result(n)
      import C_INT
      integer(C_INT) :: n
    end function
  end interface

contains

  function micro_xs_ptr() result(ptr) bind(C)
    type(C_PTR) :: ptr
    ptr = C_LOC(micro_xs(1))
  end function


!===============================================================================
! FREE_MEMORY_NUCLIDE deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_nuclide()
    interface
      subroutine library_clear() bind(C)
      end subroutine

      subroutine nuclides_clear() bind(C)
      end subroutine
    end interface

    call nuclides_clear()
    call library_clear()

  end subroutine free_memory_nuclide

end module nuclide_header
