module m_wkml_color

  use FoX_wxml, only : xmlf_t
#ifndef DUMMYLIB
  use FoX_wxml, only : xml_NewElement, xml_EndElement, xml_AddCharacters
  use FoX_common
  use m_common_error

  use m_wkml_color_def, only: colorArray
#endif

  implicit none
  private

  type color_t
    private
    character(len=8) :: hex
  end type color_t

#ifndef DUMMYLIB
  integer, parameter :: defaultAlpha = 249
#endif

  type(color_t), parameter :: defaultCI(5) = &
    (/color_t('eef00000'), color_t('eeb00030'), &
    color_t('ee700070'), color_t('ee3000b0'), color_t('ee0000f0')/)

#ifndef DUMMYLIB
  character(len=*), parameter :: hexdigits = '0123456789abcdefABCDEF'
#endif

  interface kmlGetCustomColor
    module procedure kmlGetColor_index
    module procedure kmlGetColor_byName
  end interface kmlGetCustomColor

#ifndef DUMMYLIB
  interface kmlAddColor
    module procedure kmlAddColor_c
    module procedure kmlAddColor_h
  end interface kmlAddColor

  interface kmlAddBgColor
    module procedure kmlAddBgColor_c
    module procedure kmlAddBgColor_h
  end interface kmlAddBgColor

  interface kmlAddTextColor
    module procedure kmlAddTextColor_c
    module procedure kmlAddTextColor_h
  end interface kmlAddTextColor
#endif

  public :: color_t
#ifndef DUMMYLIB
  public :: defaultCI
#endif

  public :: kmlGetCustomColor
  public :: kmlSetCustomColor

#ifndef DUMMYLIB
  public :: kmlGetColorHex

  public :: kmlAddColor
  public :: kmlAddBgColor
  public :: kmlAddTextColor

! this one is used for MCHSIM like coloring
  public :: kmlMakeColorMap
#endif

contains

#ifndef DUMMYLIB
  function checkColorHex(h) result(p)
    character(len=8) :: h
    logical :: p
    p = (verify(h, hexdigits)==0)
  end function checkColorHex

  subroutine kmlAddColor_c(xf, col)
    type(xmlf_t), intent(inout) :: xf
    type(color_t), intent(in) :: col

    call xml_NewElement(xf, 'color')
    call xml_AddCharacters(xf, col%hex)
    call xml_EndElement(xf, 'color')
  end subroutine kmlAddColor_c


  subroutine kmlAddColor_h(xf, col)
    type(xmlf_t), intent(inout) :: xf
    character(len=8), intent(in) :: col
#ifndef DUMMYLIB
    integer, pointer :: p
    
    if (.not.checkColorHex(col)) then
      print*, col
      print*, verify(col, hexdigits)
      print*, "Invalid color value"
      nullify(p)
      print*, p
    endif

    call xml_NewElement(xf, 'color')
    call xml_AddCharacters(xf, col)
    call xml_EndElement(xf, 'color')
#endif
  end subroutine kmlAddColor_h

  subroutine kmlAddBgColor_c(xf, col)
    type(xmlf_t), intent(inout) :: xf
    type(color_t), intent(in) :: col

    call xml_NewElement(xf, 'bgColor')
    call xml_AddCharacters(xf, col%hex)
    call xml_EndElement(xf, 'bgColor')
  end subroutine kmlAddBgColor_c


  subroutine kmlAddBgColor_h(xf, col)
    type(xmlf_t), intent(inout) :: xf
    character(len=8), intent(in) :: col

    if (.not.checkColorHex(col)) then
      call FoX_error("Invalid color value")
    endif

    call xml_NewElement(xf, 'bgColor')
    call xml_AddCharacters(xf, col)
    call xml_EndElement(xf, 'bgColor')
  end subroutine kmlAddBgColor_h

  subroutine kmlAddTextColor_c(xf, col)
    type(xmlf_t), intent(inout) :: xf
    type(color_t), intent(in) :: col

    call xml_NewElement(xf, 'textcolor')
    call xml_AddCharacters(xf, col%hex)
    call xml_EndElement(xf, 'textcolor')
  end subroutine kmlAddTextColor_c


  subroutine kmlAddTextColor_h(xf, col)
    type(xmlf_t), intent(inout) :: xf
    character(len=8), intent(in) :: col

    if (.not.checkColorHex(col)) then
      call FoX_error("Invalid color value")
    endif

    call xml_NewElement(xf, 'textcolor')
    call xml_AddCharacters(xf, col)
    call xml_EndElement(xf, 'textcolor')
  end subroutine kmlAddTextColor_h
#endif

  subroutine kmlSetCustomColor(myCI, colorhex)
    type(color_t), intent(out) :: myCI
    character(len=8), intent(in) :: colorhex
#ifndef DUMMYLIB
    if (.not.checkColorHex(colorhex)) then
      call FoX_error("Invalid color value")
    endif

    myCI%hex = colorhex
#endif
  end subroutine kmlSetCustomColor

#ifndef DUMMYLIB
  function kmlGetColorHex(col) result(h)
    type(color_t), intent(in) :: col
    character(len=8) :: h

    h = col%hex
  end function kmlGetColorHex
#endif

  function kmlGetColor_index(i) result(c)
    integer, intent(in) :: i
    type(color_t) :: c

#ifndef DUMMYLIB
    if (i>size(colorArray).or.(i<1)) then
      call FoX_error("Invalid index in kmlGetColor (by index)")
    endif
    
    ! Note - 'x2' argument is the format string
    c%hex = str(defaultAlpha, 'x2')//str(colorArray(i)%b, 'x2') &
      //str(colorArray(i)%g, 'x2')//str(colorArray(i)%r, 'x2')
#endif
  end function kmlGetColor_index

  function kmlGetColor_byName(name) result(c)
    character(len=*), intent(in) :: name
    type(color_t) :: c

#ifndef DUMMYLIB
    integer :: i
    logical :: found

    found = .false.
    do i = 1, size(colorArray)
      if (trim(colorArray(i)%name)==name) then
        c%hex = str(defaultAlpha, 'x2')//str(colorArray(i)%b, 'x2') &
          //str(colorArray(i)%g, 'x2')//str(colorArray(i)%r, 'x2')
        found = .true.
      endif
    enddo

    if (.not.found) then
      call FoX_error("Invalid name in kmlGetColor (by name)")
    endif
#endif
  end function kmlGetColor_byName

#ifndef DUMMYLIB
  function kmlMakeColorMap(numcolors, start, finish) result(colormap)
    integer, intent(in) :: numcolors
    type(color_t), intent(in), optional :: start, finish
    type(color_t), pointer :: colormap(:)

    integer :: i, red, blue

    if (present(start).neqv.present(finish)) then
      print*, 'Must specify both of start and finish in kmlMakeColorMap'
    endif

    ! FIXME do something with start & finish
    allocate(colormap(numcolors))
    do i = 1, numcolors
      red = (256.0/numcolors)*i-1
      blue = 256 - red
      colormap(i)%hex = "e0"//str(red, "x2")//str(40)//str(blue, "x2")
    enddo
  end function kmlMakeColorMap
#endif

end module m_wkml_color
