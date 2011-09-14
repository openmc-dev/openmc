module fission

  use constants
  use cross_section_header, only: Nuclide
  use error,                only: fatal_error
  use search,               only: binary_search

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

    integer :: i           ! loop index
    integer :: j           ! index on nu energy grid / precursor group
    integer :: loc         ! index before start of energies/nu values
    integer :: NC          ! number of polynomial coefficients
    integer :: NR          ! number of interpolation regions
    integer :: NE          ! number of energies tabulated
    real(8) :: c           ! polynomial coefficient
    real(8) :: f           ! interpolation factor
    character(MAX_LINE_LEN) :: msg ! error message

    if (nuc % nu_t_type == NU_NONE) then
       msg = "No neutron emission data for table: " // nuc % name
       call fatal_error(msg)
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
       ! determine number of interpolation regions -- as far as I can tell, no
       ! nu data has multiple interpolation regions. Furthermore, it seems all
       ! are lin-lin.
       NR = int(nuc % nu_t_data(1))
       if (NR /= 0) then
          msg = "Multiple interpolation regions not supported while &
               &attempting to determine total nu."
          call fatal_error(msg)
       end if

       ! determine number of energies
       loc = 2 + 2*NR
       NE = int(nuc % nu_t_data(loc))

       ! do binary search over tabuled energies to determine appropriate index
       ! and interpolation factor
       if (E < nuc % nu_t_data(loc+1)) then
          j = 1
          f = ZERO
       elseif (E > nuc % nu_t_data(loc+NE)) then
          j = NE - 1
          f = ONE
       else
          j = binary_search(nuc % nu_t_data(loc+1), NE, E)
          f = (E - nuc % nu_t_data(loc+j)) / &
               & (nuc % nu_t_data(loc+j+1) - nuc % nu_t_data(loc+j))
       end if

       ! determine nu total
       loc = loc + NE
       nu = nuc % nu_t_data(loc+j) + f * & 
            & (nuc % nu_t_data(loc+j+1) - nuc % nu_t_data(loc+j))
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

    integer :: i           ! loop index
    integer :: j           ! index on nu energy grid / precursor group
    integer :: loc         ! index before start of energies/nu values
    integer :: NC          ! number of polynomial coefficients
    integer :: NR          ! number of interpolation regions
    integer :: NE          ! number of energies tabulated
    real(8) :: c           ! polynomial coefficient
    real(8) :: f           ! interpolation factor
    character(MAX_LINE_LEN) :: msg ! error message

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
       ! determine number of interpolation regions
       NR = int(nuc % nu_p_data(1))
       if (NR /= 0) then
          msg = "Multiple interpolation regions not supported while & 
               &attempting to determine prompt nu."
          call fatal_error(msg)
       end if

       ! determine number of energies
       loc = 2 + 2*NR
       NE = int(nuc % nu_p_data(loc))

       ! do binary search over tabuled energies to determine appropriate index
       ! and interpolation factor
       if (E < nuc % nu_p_data(loc+1)) then
          j = 1
          f = ZERO
       elseif (E > nuc % nu_p_data(loc+NE)) then
          j = NE - 1
          f = ONE
       else
          j = binary_search(nuc % nu_p_data(loc+1), NE, E)
          f = (E - nuc % nu_p_data(loc+j)) / &
               & (nuc % nu_p_data(loc+j+1) - nuc % nu_p_data(loc+j))
       end if

       ! determine nu total
       loc = loc + NE
       nu = nuc % nu_p_data(loc+j) + f * & 
            & (nuc % nu_p_data(loc+j+1) - nuc % nu_p_data(loc+j))
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

    integer :: j           ! index on nu energy grid / precursor group
    integer :: loc         ! index before start of energies/nu values
    integer :: NR          ! number of interpolation regions
    integer :: NE          ! number of energies tabulated
    real(8) :: f           ! interpolation factor
    character(MAX_LINE_LEN) :: msg ! error message

    if (nuc % nu_d_type == NU_NONE) then
       nu = ZERO
    elseif (nuc % nu_d_type == NU_TABULAR) then
       ! determine number of interpolation regions
       NR = int(nuc % nu_d_data(1))
       if (NR /= 0) then
          msg = "Multiple interpolation regions not supported while & 
               &attempting to determine delayed nu."
          call fatal_error(msg)
       end if

       ! determine number of energies
       loc = 2 + 2*NR
       NE = int(nuc % nu_d_data(loc))

       ! do binary search over tabuled energies to determine appropriate index
       ! and interpolation factor
       if (E < nuc % nu_d_data(loc+1)) then
          j = 1
          f = ZERO
       elseif (E > nuc % nu_d_data(loc+NE)) then
          j = NE - 1
          f = ONE
       else
          j = binary_search(nuc % nu_d_data(loc+1), NE, E)
          f = (E - nuc % nu_d_data(loc+j)) / &
               & (nuc % nu_d_data(loc+j+1) - nuc % nu_d_data(loc+j))
       end if

       ! determine nu total
       loc = loc + NE
       nu = nuc % nu_d_data(loc+j) + f * & 
            & (nuc % nu_d_data(loc+j+1) - nuc % nu_d_data(loc+j))
    end if

  end function nu_delayed

end module fission
