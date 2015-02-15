module m_handlers_two

  use FoX_common
  use FoX_sax

  implicit none
  private

  ! This shows how you might use a SAX processor
  ! to extract character data from within certain nodes.

  ! It will print out the character content of all scalar
  ! nodes inside parameter nodes.

  ! It will also remember the value of the scalar within the 
  ! parameter node called DM.EnergyTolerance

  ! Note that when dealing with character data, we need to 
  ! concatenate all the character handling events in question:
  ! we do this using some of the character handling routines 
  ! provided as a convenience by FoX.

  public :: characters_handler
  public :: endElement_handler
  public :: startDocument_handler
  public :: startElement_handler

  logical :: inScalar
  logical :: inParameter
  logical :: etolFound
  real :: etol
  character(200) :: c

  public :: etol

contains

  subroutine characters_handler(chunk)
    character(len=*), intent(in) :: chunk

    if (inScalar.and.inParameter) then
      print*, "Found some scalar parameter data:"
      print*, trim(chunk)
      if (etolFound) then
        c = trim(c)//" "//chunk
      endif
    endif

  end subroutine characters_handler

  subroutine endElement_handler(URI, localname, name)
    character(len=*), intent(in)     :: URI
    character(len=*), intent(in)     :: localname
    character(len=*), intent(in)     :: name

    ! Note if we are leaving a scalar or parameter element.
    if (URI=="http://www.xml-cml.org/schema" &
      .and. localName=="scalar") then
      inScalar = .false.
      if (etolFound) then
        ! pull the data out of the concatenated string:
        call rts(c, etol)
      endif
    elseif (URI=="http://www.xml-cml.org/schema" &
      .and. localName=="parameter") then
      inParameter = .false.
      etolFound = .false.
    endif
  end subroutine endElement_handler
  
  subroutine startDocument_handler
    ! Initialize state variables
    inScalar = .false.
    inParameter = .false.
    ! Initalize other module variables
    c = ''
  end subroutine startDocument_handler
  
  subroutine startElement_handler(URI, localname, name, attributes)
    character(len=*), intent(in)   :: URI
    character(len=*), intent(in)   :: localname
    character(len=*), intent(in)   :: name
    type(dictionary_t), intent(in) :: attributes

    integer :: i

    ! Note if we are entering a scalar or parameter element.
    if (URI=="http://www.xml-cml.org/schema" &
      .and. localName=="scalar") then
      inScalar = .true.
    elseif (URI=="http://www.xml-cml.org/schema" &
      .and. localName=="parameter") then
      inParameter = .true.
      ! Loop over the attributes, looking for the name ...
      do i = 1, getLength(attributes)
        if (getQName(attributes, i)=="name") then
          ! And if the attribute is called name, check if the name
          ! is what we are interested in.
          print*, "name=", getValue(attributes, i)
          etolfound = getValue(attributes, i)=="DM.EnergyTolerance"
        endif
      enddo
    endif
  end subroutine startElement_handler

end module m_handlers_two

program sax_example

  use FoX_sax

  use m_handlers_two ! above

  implicit none

  integer :: iostat
  type(xml_t)  :: fxml

  call open_xml_file(fxml, "h2o.xml", iostat=iostat)
  if (iostat /= 0) then
    write(*,*) "Cannot open file."
    stop
  endif

  call parse(fxml,&
    characters_handler=characters_handler, &
    endElement_handler=endElement_handler, & 
    startDocument_handler=startDocument_handler, & 
    startElement_handler=startElement_handler)

  call close_xml_t(fxml)

  print*, "The energy tolerance requested was:", etol

end program sax_example
