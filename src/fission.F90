module fission

  use ace_header,    only: Nuclide
  use constants
  use error,         only: fatal_error
  use global,        only: message
  use interpolation, only: interpolate_tab1
  use search,        only: binary_search

  implicit none

contains

!===============================================================================
! NU_TOTAL calculates the total number of neutrons emitted per fission for a
! given nuclide and incoming neutron energy
!===============================================================================

  function nu_total(nuc, E) result(nu)

    type(Nuclide), pointer :: nuc ! nuclide from which to find nu
    real(8), intent(in)    :: E   ! energy of incoming neutron
    real(8)                :: nu  ! number of total neutrons emitted per fission

    integer :: i  ! loop index
    integer :: NC ! number of polynomial coefficients
    real(8) :: c  ! polynomial coefficient

    if (nuc % nu_t_type == NU_NONE) then
      message = "No neutron emission data for table: " // nuc % name
      call fatal_error()
    elseif (nuc % nu_t_type == NU_POLYNOMIAL) then
      ! determine number of coefficients
      NC = int(nuc % nu_t_data(1))

      ! sum up polynomial in energy
      nu = ZERO
      do i = 0, NC - 1
        c = nuc % nu_t_data(i+2)
        nu = nu + c * E**i
      end do
    elseif (nuc % nu_t_type == NU_TABULAR) then
      ! use ENDF interpolation laws to determine nu
      nu = interpolate_tab1(nuc % nu_t_data, E)
    end if

  end function nu_total

!===============================================================================
! NU_PROMPT calculates the total number of prompt neutrons emitted per fission
! for a given nuclide and incoming neutron energy
!===============================================================================

  function nu_prompt(nuc, E) result(nu)

    type(Nuclide), pointer :: nuc ! nuclide from which to find nu
    real(8), intent(in)    :: E   ! energy of incoming neutron
    real(8)                :: nu  ! number of prompt neutrons emitted per fission

    integer :: i  ! loop index
    integer :: NC ! number of polynomial coefficients
    real(8) :: c  ! polynomial coefficient

    if (nuc % nu_p_type == NU_NONE) then
      ! since no prompt or delayed data is present, this means all neutron
      ! emission is prompt -- WARNING: This currently returns zero. The calling
      ! routine needs to know this situation is occurring since we don't want
      ! to call nu_total unnecessarily if it's already been called
      nu = ZERO
    elseif (nuc % nu_p_type == NU_POLYNOMIAL) then
      ! determine number of coefficients
      NC = int(nuc % nu_p_data(1))

      ! sum up polynomial in energy
      nu = ZERO
      do i = 0, NC - 1
        c = nuc % nu_p_data(i+2)
        nu = nu + c * E**i
      end do
    elseif (nuc % nu_p_type == NU_TABULAR) then
      ! use ENDF interpolation laws to determine nu
      nu = interpolate_tab1(nuc % nu_p_data, E)
    end if

  end function nu_prompt

!===============================================================================
! NU_DELAYED calculates the total number of delayed neutrons emitted per fission
! for a given nuclide and incoming neutron energy
!===============================================================================

  function nu_delayed(nuc, E) result(nu)

    type(Nuclide), pointer :: nuc ! nuclide from which to find nu
    real(8), intent(in)    :: E   ! energy of incoming neutron
    real(8)                :: nu  ! number of delayed neutrons emitted per fission

    if (nuc % nu_d_type == NU_NONE) then
      nu = ZERO
    elseif (nuc % nu_d_type == NU_TABULAR) then
      ! use ENDF interpolation laws to determine nu
      nu = interpolate_tab1(nuc % nu_d_data, E)
    end if

  end function nu_delayed

end module fission
