module tally_derivative_header

  use constants
  use dict_header, only: DictIntInt
  use error, only: fatal_error
  use nuclide_header, only: nuclide_dict
  use string, only: to_str, to_lower
  use xml_interface

  implicit none
  private

!===============================================================================
! TALLYDERIVATIVE describes a first-order derivative that can be applied to
! tallies.
!===============================================================================

  type, public :: TallyDerivative
    integer :: id
    integer :: variable
    integer :: diff_material
    integer :: diff_nuclide
    real(8) :: flux_deriv
  contains
    procedure :: from_xml
  end type TallyDerivative

  type(TallyDerivative), public, allocatable :: tally_derivs(:)
!$omp threadprivate(tally_derivs)

  ! Dictionary that maps user IDs to indices in 'tally_derivs'
  type(DictIntInt), public :: tally_deriv_dict

contains

  subroutine from_xml(this, node)
    class(TallyDerivative), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    character(MAX_WORD_LEN) :: temp_str
    character(MAX_WORD_LEN) :: word

    ! Copy the derivative id.
    if (check_for_node(node, "id")) then
      call get_node_value(node, "id", this % id)
    else
      call fatal_error("Must specify an ID for <derivative> elements in the&
           & tally XML file")
    end if

    ! Make sure the id is > 0.
    if (this % id <= 0) then
      call fatal_error("<derivative> IDs must be an integer greater than &
           &zero")
    end if

    ! Make sure this id has not already been used.
    if (tally_deriv_dict % has_key(this % id)) then
      call fatal_error("Two or more <derivative>'s use the same unique &
           &ID: " // trim(to_str(this % id)))
    end if

    ! Read the independent variable name.
    call get_node_value(node, "variable", temp_str)
    temp_str = to_lower(temp_str)

    select case(temp_str)
    case("density")
      this % variable = DIFF_DENSITY

    case("nuclide_density")
      this % variable = DIFF_NUCLIDE_DENSITY

      call get_node_value(node, "nuclide", word)
      word = trim(to_lower(word))
      if (.not. nuclide_dict % has_key(word)) then
        call fatal_error("Could not find the nuclide " &
             // trim(word) // " specified in derivative " &
             // trim(to_str(this % id)) // " in any material.")
      end if
      this % diff_nuclide = nuclide_dict % get_key(word)

    case("temperature")
      this % variable = DIFF_TEMPERATURE
    end select

    call get_node_value(node, "material", this % diff_material)

  end subroutine from_xml

end module tally_derivative_header
