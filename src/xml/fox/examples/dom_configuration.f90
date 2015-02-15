program dom_configuration
! Use FoX DOM API to determine 
! default values and settablility 
! of DOM Configuration options.

  use FoX_dom

  implicit none 
 
  integer, parameter :: max_name_len = 50
  integer, parameter :: max_configs = 50
  type(DOMConfiguration), pointer :: config
  integer :: num_configs
  integer :: config_name_len
  character(len=max_name_len), dimension(max_configs) :: config_names
  integer :: i
  integer :: defaults
  integer :: settable
! Uncomment these to help with changing settings
!  logical :: set
!  logical :: def
  
  settable = 0
  defaults = 0

  config => newDOMConfig()
  config_name_len = len(getParameterNames(config))
  num_configs = size(getParameterNames(config))
  config_names(1:num_configs) = getParameterNames(config)
  write (*,"(a)") '=========================================================================='
  write (*,"(a3, 3x, a7, 3x, a7, 3x, a)"), 'num', 'can set', 'default', 'name'
  write (*,"(a)") '=========================================================================='
  ! Loop over all configuration options and report settability and default value.
  do i = 1, num_configs
    write (*,"(i3, 6x, l1, 9x, l1, 6x, a)") i, & 
      &  canSetParameter(config, config_names(i), .true.), & 
      &  getParameter(config, config_names(i)), trim(config_names(i))
    if (canSetParameter(config, config_names(i), .true.)) then
      settable = ibset(settable, i)
    else
      settable = ibclr(settable, i)
    endif
    if (getParameter(config, config_names(i))) then
      defaults = ibset(defaults, i)
    else
      defaults = ibclr(defaults, i)
    endif
  end do
  write (*,"(a)") '=========================================================================='
  ! These numbers should be parameters inside m_dom_configuration.m4...
  print*, 'defaults = ', defaults
  print*, 'settable = ', settable
  call destroy(config)

! Uncomment these to help with changing settings in m_dom_configurations.m4
!  do
!    read(*,*) i, set, def
!    if (set) then 
!      settable = ibset(settable, i)
!    else 
!      settable = ibclr(settable, i)
!    endif
!    if (def) then
!      defaults = ibset(defaults, i)
!    else
!      defaults = ibclr(defaults, i)
!    endif
!    print*, 'defaults = ', defaults
!    print*, 'settable = ', settable
!  enddo

end program dom_configuration
