module tally_filter_particle

  use, intrinsic :: ISO_C_BINDING

  use tally_filter_header

  implicit none
  private

!===============================================================================
! PARTICLEFILTER specifies what particles can score to a tally
!===============================================================================

  type, public, extends(TallyFilter) :: ParticleFilter
  contains
    procedure :: particles
  end type ParticleFilter

contains

  function particles(this, i) result(ptype)
    class(ParticleFilter), intent(in) :: this
    integer,               intent(in) :: i
    integer                           :: ptype
    interface
      function particle_filter_particles(filt, i) result(ptype) bind(C)
        import C_PTR, C_INT
        type(C_PTR),    value :: filt
        integer(C_INT), value :: i
        integer(C_INT)        :: ptype
      end function
    end interface
    ptype = particle_filter_particles(this % ptr, i)
  end function particles

end module tally_filter_particle
