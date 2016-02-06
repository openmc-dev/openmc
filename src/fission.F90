module fission

  use nuclide_header, only: NuclideCE
  use constants
  use error,          only: fatal_error
  use interpolation,  only: interpolate_tab1
  use search,         only: binary_search

  implicit none

contains

!===============================================================================
! NU_TOTAL calculates the total number of neutrons emitted per fission for a
! given nuclide and incoming neutron energy
!===============================================================================

  pure function nu_total(nuc, E) result(nu)
    type(NuclideCE), intent(in) :: nuc ! nuclide from which to find nu
    real(8),       intent(in) :: E   ! energy of incoming neutron
    real(8)                   :: nu  ! number of total neutrons emitted per fission

    integer :: i  ! loop index
    integer :: NC ! number of polynomial coefficients
    real(8) :: c  ! polynomial coefficient

    if (nuc % nu_t_type == NU_NONE) then
      nu = ERROR_REAL
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

  pure function nu_prompt(nuc, E) result(nu)
    type(NuclideCE), intent(in) :: nuc ! nuclide from which to find nu
    real(8),       intent(in) :: E   ! energy of incoming neutron
    real(8)                   :: nu  ! number of prompt neutrons emitted per fission

    integer :: i  ! loop index
    integer :: NC ! number of polynomial coefficients
    real(8) :: c  ! polynomial coefficient

    if (nuc % nu_p_type == NU_NONE) then
      ! since no prompt or delayed data is present, this means all neutron
      ! emission is prompt -- WARNING: This currently returns zero. The calling
      ! routine needs to know this situation is occurring since we don't want
      ! to call nu_total unnecessarily if it has already been called.
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

  pure function nu_delayed(nuc, E) result(nu)
    type(NuclideCE), intent(in) :: nuc ! nuclide from which to find nu
    real(8),       intent(in) :: E   ! energy of incoming neutron
    real(8)                   :: nu  ! number of delayed neutrons emitted per fission

    if (nuc % nu_d_type == NU_NONE) then
      ! since no prompt or delayed data is present, this means all neutron
      ! emission is prompt -- WARNING: This currently returns zero. The calling
      ! routine needs to know this situation is occurring since we don't want
      ! to call nu_delayed unnecessarily if it has already been called.
      nu = ZERO
    elseif (nuc % nu_d_type == NU_TABULAR) then
      ! use ENDF interpolation laws to determine nu
      nu = interpolate_tab1(nuc % nu_d_data, E)
    end if

  end function nu_delayed

!===============================================================================
! YIELD_DELAYED calculates the fractional yield of delayed neutrons emitted for
! a given nuclide and incoming neutron energy in a given delayed group.
!===============================================================================

  pure function yield_delayed(nuc, E, g) result(yield)
    type(NuclideCE), intent(in) :: nuc   ! nuclide from which to find nu
    real(8), intent(in)       :: E     ! energy of incoming neutron
    real(8)                   :: yield ! delayed neutron precursor yield
    integer, intent(in)       :: g     ! the delayed neutron precursor group
    integer                   :: d     ! precursor group
    integer                   :: lc    ! index before start of energies/nu values
    integer                   :: NR    ! number of interpolation regions
    integer                   :: NE    ! number of energies tabulated

    yield = ZERO

    if (g > nuc % n_precursor .or. g < 1) then
      ! if the precursor group is outside the range of precursor groups for
      ! the input nuclide, return ZERO.
      yield = ZERO
    else if (nuc % nu_d_type == NU_NONE) then
      ! since no prompt or delayed data is present, this means all neutron
      ! emission is prompt -- WARNING: This currently returns zero. The calling
      ! routine needs to know this situation is occurring since we don't want
      ! to call yield_delayed unnecessarily if it has already been called.
      yield = ZERO
    else if (nuc % nu_d_type == NU_TABULAR) then

      lc = 1

      ! loop over delayed groups and determine the yield for the desired group
      do d = 1, nuc % n_precursor

        ! determine number of interpolation regions and energies
        NR = int(nuc % nu_d_precursor_data(lc + 1))
        NE = int(nuc % nu_d_precursor_data(lc + 2 + 2*NR))

        ! check if this is the desired group
        if (d == g) then

          ! determine delayed neutron precursor yield for group g
          yield = interpolate_tab1(nuc % nu_d_precursor_data( &
               lc+1:lc+2+2*NR+2*NE), E)

          exit
        end if

        ! advance pointer
        lc = lc + 2 + 2*NR + 2*NE + 1
      end do
    end if

  end function yield_delayed

end module fission
