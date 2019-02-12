module math

  use, intrinsic :: ISO_C_BINDING

  implicit none
  private
  public :: t_percentile

  interface

    pure function t_percentile(p, df) bind(C) result(t)
      use ISO_C_BINDING
      implicit none
      real(C_DOUBLE), value, intent(in) :: p
      integer(C_INT), value, intent(in) :: df
      real(C_DOUBLE) :: t
    end function t_percentile
  end interface

end module math
