module tally_filter_particle

  use, intrinsic :: ISO_C_BINDING

  use constants,          only: ONE, MAX_LINE_LEN
  use hdf5_interface
  use particle_header,    only: Particle
  use string,             only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private

!===============================================================================
! PARTICLEFILTER specifies what particles can score to a tally
!===============================================================================

  type, public, extends(TallyFilter) :: ParticleFilter
    integer, allocatable :: particles(:)
  contains
    procedure :: from_xml
    procedure :: get_all_bins
    procedure :: to_statepoint
    procedure :: text_label
  end type ParticleFilter

contains

  subroutine from_xml(this, node)
    class(ParticleFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: n

    ! Determine how many bins were given
    n = node_word_count(node, "bins")

    ! Allocate and store bins
    this % n_bins = n
    allocate(this % particles(n))
    call get_node_array(node, "bins", this % particles)
  end subroutine from_xml

  subroutine get_all_bins(this, p, estimator, match)
    class(ParticleFilter), intent(in)  :: this
    type(Particle),    intent(in)  :: p
    integer,           intent(in)  :: estimator
    type(TallyFilterMatch), intent(inout) :: match

    integer :: i

    do i = 1, this % n_bins
      if (this % particles(i) == p % type) then
        call match % bins % push_back(i)
        call match % weights % push_back(ONE)
      end if
    end do
  end subroutine get_all_bins

  subroutine to_statepoint(this, filter_group)
    class(ParticleFilter), intent(in) :: this
    integer(HID_T),    intent(in) :: filter_group

    call write_dataset(filter_group, "type", "particle")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", this % particles)
  end subroutine to_statepoint

  function text_label(this, bin) result(label)
    class(ParticleFilter), intent(in) :: this
    integer,           intent(in) :: bin
    character(MAX_LINE_LEN)       :: label

    label = "Particle " // to_str(this % particles(bin))
  end function text_label

end module tally_filter_particle
