module tally_filter_header

  use particle_header, only: Particle

  use hdf5

  implicit none

  type, abstract :: TallyFilter
    integer :: n_bins = 0
    integer :: type
  contains
    procedure(get_next_bin_), deferred :: get_next_bin
    procedure(get_score_), deferred :: get_score
    procedure(to_statepoint_), deferred :: to_statepoint
    procedure(to_summary_), deferred :: to_summary
    procedure(initialize_), deferred :: initialize
  end type TallyFilter

  type TallyFilterContainer
    class(TallyFilter), allocatable :: obj
  end type TallyFilterContainer

  abstract interface

    function get_next_bin_(this, p, estimator, current_bin) result(next_bin)
      import TallyFilter
      import Particle
      class(TallyFilter), intent(in) :: this
      type(Particle),     intent(in) :: p
      integer,            intent(in) :: estimator
      integer,            intent(in) :: current_bin
      integer                        :: next_bin
    end function get_next_bin_

    function get_score_(this, bin) result(score)
      import TallyFilter
      class(TallyFilter), intent(in) :: this
      integer,            intent(in) :: bin
      real(8)                        :: score
    end function get_score_

    subroutine to_statepoint_(this, filter_group)
      import TallyFilter
      import HID_T
      class(TallyFilter), intent(in) :: this
      integer(HID_T),     intent(in) :: filter_group
    end subroutine to_statepoint_

    subroutine to_summary_(this, filter_group)
      import TallyFilter
      import HID_T
      class(TallyFilter), intent(in) :: this
      integer(HID_T),     intent(in) :: filter_group
    end subroutine to_summary_

    subroutine initialize_(this)
      import TallyFilter
      class(TallyFilter), intent(inout) :: this
    end subroutine initialize_

  end interface

end module tally_filter_header
