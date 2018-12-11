module eigenvalue

  use, intrinsic :: ISO_C_BINDING

  implicit none

  interface
    function openmc_get_keff(k_combined) result(err) bind(C)
      import C_INT, C_DOUBLE
      real(C_DOUBLE), intent(out) :: k_combined(2)
      integer(C_INT) :: err
    end function
  end interface

end module eigenvalue
