module volume_header

  use constants, only: FILTER_CELL, FILTER_MATERIAL, FILTER_UNIVERSE
  use error,     only: fatal_error
  use xml_interface

  implicit none

  type VolumeCalculation
    integer :: domain_type
    integer, allocatable :: domain_id(:)
    real(8) :: lower_left(3)
    real(8) :: upper_right(3)
    integer :: samples
  contains
    procedure :: from_xml => volume_from_xml
  end type VolumeCalculation

  type(VolumeCalculation), allocatable :: volume_calcs(:)

contains

  subroutine volume_from_xml(this, node_vol)
    class(VolumeCalculation), intent(out) :: this
    type(XMLNode), intent(in) :: node_vol

    integer :: num_domains
    character(10) :: temp_str

    ! Check domain type
    call get_node_value(node_vol, "domain_type", temp_str)
    select case (temp_str)
    case ('cell')
      this % domain_type = FILTER_CELL
    case ('material')
      this % domain_type = FILTER_MATERIAL
    case ('universe')
      this % domain_type = FILTER_UNIVERSE
    case default
      call fatal_error("Unrecognized domain type for stochastic volume &
           &calculation: " // trim(temp_str))
    end select

    ! Read cell IDs
    if (check_for_node(node_vol, "domain_ids")) then
      num_domains = node_word_count(node_vol, "domain_ids")
    else
      call fatal_error("Must specify at least one cell for a volume calculation")
    end if
    allocate(this % domain_id(num_domains))
    call get_node_array(node_vol, "domain_ids", this % domain_id)

    ! Read lower-left and upper-right bounding coordinates
    call get_node_array(node_vol, "lower_left", this % lower_left)
    call get_node_array(node_vol, "upper_right", this % upper_right)

    ! Read number of samples
    call get_node_value(node_vol, "samples", this % samples)
  end subroutine volume_from_xml

!===============================================================================
! FREE_MEMORY_VOLUME deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_volume()
    if (allocated(volume_calcs)) deallocate(volume_calcs)
  end subroutine free_memory_volume

end module volume_header
