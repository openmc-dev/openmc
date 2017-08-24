module tally_filter_mu

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use algorithm,          only: binary_search
  use constants,          only: ONE, MAX_LINE_LEN, NO_BIN_FOUND
  use error,              only: fatal_error
  use hdf5_interface
  use particle_header,    only: Particle
  use string,             only: to_str
  use tally_filter_header

  implicit none
  private

!===============================================================================
! MUFILTER bins the incoming-outgoing direction cosine.  This is only used for
! scatter reactions.
!===============================================================================

  type, public, extends(TallyFilter) :: MuFilter
    real(8), allocatable :: bins(:)
  contains
    procedure :: get_all_bins => get_all_bins_mu
    procedure :: to_statepoint => to_statepoint_mu
    procedure :: text_label => text_label_mu
  end type MuFilter

contains

  subroutine get_all_bins_mu(this, p, estimator, match)
    class(MuFilter), intent(in)  :: this
    type(Particle),  intent(in)  :: p
    integer,         intent(in)  :: estimator
    type(TallyFilterMatch), intent(inout) :: match

    integer :: n
    integer :: bin

    n = this % n_bins

    ! Search to find incoming energy bin.
    bin = binary_search(this % bins, n + 1, p % mu)
    if (bin /= NO_BIN_FOUND) then
      call match % bins % push_back(bin)
      call match % weights % push_back(ONE)
    end if
  end subroutine get_all_bins_mu

  subroutine to_statepoint_mu(this, filter_group)
    class(MuFilter), intent(in) :: this
    integer(HID_T),  intent(in) :: filter_group

    call write_dataset(filter_group, "type", "mu")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % bins)
  end subroutine to_statepoint_mu

  function text_label_mu(this, bin) result(label)
    class(MuFilter), intent(in) :: this
    integer,         intent(in) :: bin
    character(MAX_LINE_LEN)     :: label

    real(8) :: E0, E1

    E0 = this % bins(bin)
    E1 = this % bins(bin + 1)
    label = "Change-in-Angle [" // trim(to_str(E0)) // ", " &
         // trim(to_str(E1)) // ")"
  end function text_label_mu

end module tally_filter_mu
