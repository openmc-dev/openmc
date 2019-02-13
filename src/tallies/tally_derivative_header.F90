module tally_derivative_header

  use, intrinsic :: ISO_C_BINDING

  implicit none

!===============================================================================
! TALLYDERIVATIVE describes a first-order derivative that can be applied to
! tallies.
!===============================================================================

  type, bind(C), public :: TallyDerivative
    integer(C_INT) :: id
    integer(C_INT) :: variable
    integer(C_INT) :: diff_material
    integer(C_INT) :: diff_nuclide
    real(C_DOUBLE) :: flux_deriv
  end type TallyDerivative

  interface
    function n_tally_derivs() result(n) bind(C)
      import C_INT
      integer(C_INT) :: n
    end function

    function tally_deriv_c(i) result(deriv) bind(C)
      import C_INT, TallyDerivative
      integer(C_INT), value :: i
      type(TallyDerivative), pointer :: deriv
    end function
  end interface

end module tally_derivative_header
