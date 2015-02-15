module m_handlers

  use FoX_common
  use FoX_sax

  implicit none
  private

    public :: error_handler
    public :: fatalError_handler

contains

  subroutine error_handler(msg)
    character(len=*), intent(in) :: msg
    write(*,'(a)') "Error encountered and caught"
    write(*,'(a)') msg
  end subroutine error_handler

  subroutine fatalError_handler(msg)
    character(len=*), intent(in) :: msg
    write(*,'(a)') "Fatal error encountered and caught"
    write(*,'(a)') msg
  end subroutine fatalError_handler

end module m_handlers

program sax_well_formed
  !
  ! Check for well-formedness - nothing happens on success;
  ! error message omitted and execution halt on any failure.
  !
  use FoX_sax
  use m_handlers
  implicit none

  integer :: iostat
  type(xml_t)  :: fxml

  call open_xml_file(fxml, "test.xml", iostat=iostat)
  if (iostat /= 0) then
    write(*,*) "Cannot open file."
    stop
  endif

  call parse(fxml,&
    error_handler=error_handler,                 &
    fatalError_handler=fatalError_handler, &
    namespaces=.true.,  &
    validate=.true.)

  call close_xml_t(fxml)

  print*, "Finished"

end program sax_well_formed
