!===============================================================================
! EXPAND_NATURAL_ELEMENT converts natural elements specified using an <element>
! tag within a material into individual isotopes based on data from the NIST
! Atomic Weights and Isotopic Compositions database.  It was last updated July
! 2010.  In some cases, modifications have been made to work with ENDF/B-VII.1
! where evaluations of particular isotopes don't exist.  The case statements in
! this subroutine were automatically written by the 'write_atomic_weights.py'
! script in the /src/utils directory.
!===============================================================================

  subroutine expand_natural_element(name, xs, density, list_names, &
       list_density)

    character(*),   intent(in)    :: name
    character(*),   intent(in)    :: xs
    real(8),        intent(in)    :: density
    type(ListChar), intent(inout) :: list_names
    type(ListReal), intent(inout) :: list_density

    character(2) :: element_name

    element_name = name(1:2)
    call lower_case(element_name)

    select case (element_name)

    !### CASES GO HERE ###
    case default
      message = "Cannot expand element: " // name
      call fatal_error()

    end select

  end subroutine expand_natural_element

