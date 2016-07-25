module multipole_header

  implicit none

  !========================================================================
  ! Multipole related constants

  ! Formalisms
  integer, parameter :: FORM_MLBW = 2, &
                        FORM_RM   = 3, &
                        FORM_RML  = 7

  ! Constants that determine which value to access
  integer, parameter :: MP_EA = 1       ! Pole

  ! Reich-Moore indices
  integer, parameter :: RM_RT = 2, &    ! Residue total
                        RM_RA = 3, &    ! Residue absorption
                        RM_RF = 4       ! Residue fission

  ! Multi-level Breit Wigner indices
  integer, parameter :: MLBW_RT = 2, &  ! Residue total
                        MLBW_RX = 3, &  ! Residue compettitive
                        MLBW_RA = 4, &  ! Residue absorption
                        MLBW_RF = 5     ! Residue fission

  ! Polynomial fit indices
  integer, parameter :: FIT_T = 1, &    ! Total
                        FIT_A = 2, &    ! Absorption
                        FIT_F = 3       ! Fission

  ! Value of 'true' when checking if nuclide is fissionable
  integer, parameter :: MP_FISS = 1

!===============================================================================
! MULTIPOLE contains all the components needed for the windowed multipole
! temperature dependent cross section libraries for the resolved resonance
! region.
!===============================================================================

  type MultipoleArray

    !=========================================================================
    ! Isotope Properties
    logical                 :: fissionable = .false.  ! Is this isotope fissionable?
    integer                 :: length                 ! Number of poles
    integer, allocatable    :: l_value(:)             ! The l index of the pole
    real(8), allocatable    :: pseudo_k0RS(:)         ! The value (sqrt(2*mass neutron)/reduced planck constant)
                                                      ! * AWR/(AWR + 1) * scattering radius for each l
    complex(8), allocatable :: data(:,:)              ! Contains all of the pole-residue data
    real(8)                 :: sqrtAWR                ! Square root of the atomic weight ratio

    !=========================================================================
    ! Windows

    integer :: windows                      ! Number of windows
    integer :: fit_order                    ! Order of the fit. 1 linear, 2 quadratic, etc.
    real(8) :: start_E                      ! Start energy for the windows
    real(8) :: end_E                        ! End energy for the windows
    real(8) :: spacing                      ! The actual spacing in sqrt(E) space.
    ! spacing = sqrt(multipole_w%endE - multipole_w%startE)/multipole_w%windows
    integer, allocatable :: w_start(:)      ! Contains the index of the pole at the start of the window
    integer, allocatable :: w_end(:)        ! Contains the index of the pole at the end of the window
    real(8), allocatable :: curvefit(:,:,:) ! Contains the fitting function.  (reaction type, coeff index, window index)

    integer, allocatable :: broaden_poly(:) ! if 1, broaden, if 0, don't.

    !=========================================================================
    ! Storage Helpers
    integer :: num_l
    integer :: max_w

    integer :: formalism

    contains
      procedure :: allocate => multipole_allocate ! Allocates Multipole
  end type MultipoleArray

contains

!===============================================================================
! MULTIPOLE_ALLOCATE allocates necessary data for Multipole.
!===============================================================================

    subroutine multipole_allocate(multipole)
      class(MultipoleArray), intent(inout)        :: multipole   ! Multipole object to allocate.

      ! This function assumes length, numL, fissionable, windows, fitorder,
      ! and formalism are known

      ! Allocate the pole-residue storage.
      ! MLBW has one more pole than Reich-Moore, and fissionable nuclides
      ! have further one more.
      if (multipole % formalism == FORM_MLBW) then
        if (multipole % fissionable) then
          allocate(multipole % data(5, multipole % length))
        else
          allocate(multipole % data(4, multipole % length))
        end if
      else if (multipole % formalism == FORM_RM) then
        if (multipole % fissionable) then
          allocate(multipole % data(4, multipole % length))
        else
          allocate(multipole % data(3, multipole % length))
        end if
      end if

      ! Allocate the l value table for each pole-residue set.
      allocate(multipole % l_value(multipole % length))

      ! Allocate the table of pseudo_k0RS values at each l.
      allocate(multipole % pseudo_k0RS(multipole % num_l))

      ! Allocate window start, window end
      allocate(multipole % w_start(multipole % windows))
      allocate(multipole % w_end(multipole % windows))

      ! Allocate broaden_poly
      allocate(multipole % broaden_poly(multipole % windows))

      ! Allocate curvefit
      if(multipole % fissionable) then
        allocate(multipole % curvefit(FIT_F, multipole % fit_order+1, multipole % windows))
      else
        allocate(multipole % curvefit(FIT_A, multipole % fit_order+1, multipole % windows))
      end if
    end subroutine
end module multipole_header
