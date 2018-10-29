module tally_filter_sph_harm

  use, intrinsic :: ISO_C_BINDING

  use tally_filter_header

  implicit none
  private

  integer, public, parameter :: COSINE_SCATTER = 1
  integer, public, parameter :: COSINE_PARTICLE = 2

  type, public, extends(TallyFilter) :: SphericalHarmonicsFilter
  contains
    procedure :: cosine
  end type SphericalHarmonicsFilter

contains

  function cosine(this) result(val)
    class(SphericalHarmonicsFilter) :: this
    integer                         :: val
    interface
      function sphharm_filter_get_cosine(filt) result(val) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: filt
        integer(C_INT) :: val
      end function
    end interface
    val = sphharm_filter_get_cosine(this % ptr)
  end function cosine

end module tally_filter_sph_harm
