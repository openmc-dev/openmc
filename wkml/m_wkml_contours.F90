module m_wkml_contours
#ifndef DUMMYLIB
  use m_common_error
  use m_contours, only: contourObject, destroy_contourObject
  use m_contours, only: make_contours_on_simplest_grid
  use m_contours, only: make_contours_on_a_simple_grid
  use m_contours, only: make_contours_on_a_complex_grid
#endif

  use fox_m_fsys_realtypes, only : sp, dp
  use m_wkml_color, only: color_t
  use FoX_wxml, only: xmlf_t

#ifndef DUMMYLIB
  use FoX_common, only: str, operator(//)

  use m_wkml_color, only: kmlGetColorHex
  use m_wkml_lowlevel, only: kmlOpenFolder, kmlCloseFolder, kmlOpenPlacemark, kmlClosePlacemark
  use m_wkml_features, only: kmlCreateLine, kmlStartRegion, kmlAddInnerBoundary, kmlEndRegion
  use m_wkml_styling, only: kmlCreateLineStyle, kmlCreatePolygonStyle 
#endif
 
  implicit none
  private

  interface kmlCreateContours
    module procedure kmlCreateContours_sp
    module procedure kmlCreateContours_longlat_sp
    module procedure kmlCreateContours_longlat2_sp
!    module procedure kmlCreateContours_dp
!    module procedure kmlCreateContours_longlat_dp
  end interface kmlCreateContours

  public :: kmlCreateContours

contains

  subroutine kmlCreateContours_sp(xf, east, west, north, south, values, &
    mask, colormap, height, contour_values, num_levels, name, lines, regions)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: east, west, north, south
    real(sp), intent(in), optional :: values(:,:)
    real(sp), intent(in), optional :: mask
    type(color_t), intent(in), optional :: colormap(:)
    real(sp), intent(in), optional :: height
    real(sp), intent(in), optional :: contour_values(:)
    integer, intent(in), optional :: num_levels
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: lines, regions

#ifndef DUMMYLIB
    type(contourObject) :: o

    o = make_contours_on_simplest_grid(east, west, north, south, values, contour_values, num_levels, mask)

    call kmlOutputContours(xf, o, height, colormap, name, lines, regions)

    call destroy_contourObject(o)
#endif
  end subroutine kmlCreateContours_sp

  subroutine kmlCreateContours_longlat_sp(xf, longitude, latitude, values, &
    mask, colormap, height, contour_values, num_levels, name, lines, regions)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: longitude(:), latitude(:)
    real(sp), intent(in), optional :: values(:,:)
    real(sp), intent(in), optional :: mask
    type(color_t), intent(in), optional :: colormap(:)
    real(sp), intent(in), optional :: height
    real(sp), intent(in), optional :: contour_values(:)
    integer, intent(in), optional :: num_levels
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: lines, regions

#ifndef DUMMYLIB
    type(contourObject) :: o

    o = make_contours_on_a_simple_grid(longitude, latitude, values, contour_values, num_levels, mask)

    call kmlOutputContours(xf, o, height, colormap, name, lines, regions)

    call destroy_contourObject(o)
#endif
  end subroutine kmlCreateContours_longlat_sp

  subroutine kmlCreateContours_longlat2_sp(xf, longitude, latitude, values, &
    mask, colormap, height, contour_values, num_levels, name, lines, regions)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: longitude(:,:), latitude(:,:)
    real(sp), intent(in), optional :: values(:,:)
    real(sp), intent(in), optional :: mask
    type(color_t), intent(in), optional :: colormap(:)
    real(sp), intent(in), optional :: height
    real(sp), intent(in), optional :: contour_values(:)
    integer, intent(in), optional :: num_levels
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: lines, regions

#ifndef DUMMYLIB
    type(contourObject) :: o

    o = make_contours_on_a_complex_grid(longitude, latitude, values, contour_values, num_levels, mask)

    call kmlOutputContours(xf, o, height, colormap, name, lines, regions)

    call destroy_contourObject(o)
#endif
  end subroutine kmlCreateContours_longlat2_sp

#ifndef DUMMYLIB
  subroutine kmlOutputContours(xf, o, height, colormap, name, lines, regions)
    type(xmlf_t), intent(inout) :: xf
    type(contourObject), intent(in) :: o
    real(sp), intent(in), optional :: height
    type(color_t), intent(in), optional :: colormap(:)
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: lines, regions

    logical :: lines_, regions_

    if (present(lines)) then
      lines_ = lines
    else
      lines_ = .false.
    endif
    if (present(regions)) then
      regions_ = regions
    else
      regions_ = .false.
    endif

    call kmlOpenFolder(xf, name=name)
    if (lines_) call outputContourLines(xf, o, height=height, colormap=colormap)
    if (regions_) call outputContourRegions(xf, o, height=height, colormap=colormap)
    call kmlCloseFolder(xf)

  end subroutine kmlOutputContours


  subroutine kmlCreateContours_old(xf, z, x, y, x_complex, y_complex, &
    contour_values, ncv, ignore_gt, lines, regions, colormap, name)
    type(xmlf_t), intent(inout) :: xf
    real, intent(in) :: z(:,:)
    real, intent(in), optional :: x(:), y(:)
    real, intent(in), optional :: x_complex(:,:), y_complex(:,:)
    real, intent(in), optional :: contour_values(:)
    integer, intent(in), optional :: ncv
    real, intent(in), optional :: ignore_gt
    logical, intent(in), optional :: lines, regions
    type(color_t), optional :: colormap(:)
    character(len=*), intent(in), optional :: name

    logical :: lines_, regions_
    type(contourObject) :: o

    if (.not.(present(x).and.present(y).and. &
      .not.present(x_complex).and..not.present(y_complex)) &
      .and..not.(present(x_complex).and.present(y_complex).and. &
      .not.present(x).and..not.present(y))) then
      call FoX_error('Can only specify simple x & y OR complex x & y')
    endif

    if (present(lines)) then
      lines_ = lines
    else
      lines_ = .false.
    endif
    if (present(regions)) then
      regions_ = regions
    else
      regions_ = .false.
    endif
    
    if (present(x)) then
      o = make_contours_on_a_simple_grid(x, y, z, contour_values, ncv, ignore_gt)
    else
      o = make_contours_on_a_complex_grid(x_complex, y_complex, z, contour_values, ncv, ignore_gt)
    endif

    call kmlOpenFolder(xf, name=name)
    if (lines_) call outputContourLines(xf, o, colormap=colormap)
    if (regions_) call outputContourRegions(xf, o, colormap=colormap)
    call kmlCloseFolder(xf)

    call destroy_contourObject(o)

  end subroutine kmlCreateContours_old


  subroutine outputContourLines(xf, cp, height, colormap)
    type(xmlf_t), intent(inout) :: xf
    type(contourObject), intent(in) :: cp
    real, intent(in), optional :: height
    type(color_t), intent(in), optional :: colormap(:)

    integer :: colornum
    integer :: i, j

    if (present(colormap)) then
      if (size(colormap)<size(cp%contours)) then
        call FoX_error("More contours than colours in outputContourLines")
      endif
    endif

    call kmlOpenFolder(xf, name='Contour Lines')

    call kmlOpenPlacemark(xf, name='Bounding Box')
    call kmlCreateLineStyle(xf, width=5, colorhex='ffffffff')
    call kmlCreateLine(xf, cp%lines(0)%list(1)%x, cp%lines(0)%list(1)%y)
    call kmlClosePlacemark(xf)

    do i = 1, size(cp%contours)
      if (i==size(cp%contours).and.associated(cp%zmax)) then
        call kmlOpenFolder(xf, name='Border of undefined region')
      else
        call kmlOpenFolder(xf, name=str(cp%contours(i)))
      endif
      do j = 1, size(cp%lines(i)%list)
        call kmlOpenPlacemark(xf)
        if (present(colormap)) then
          call kmlCreateLineStyle(xf, width=3, color=colormap(i))
        else
          colornum = i*255/(size(cp%contours)+1)
          call kmlCreateLineStyle(xf, width=3, colorhex='ff'//str(255-colornum, 'x2')//'00'//str(colornum, 'x2'))
        endif
        call kmlCreateLine(xf, cp%lines(i)%list(j)%x, cp%lines(i)%list(j)%y)
        call kmlClosePlacemark(xf)
      enddo
      call kmlCloseFolder(xf)
    enddo

    call kmlCloseFolder(xf) ! Contour Lines


  end subroutine outputContourLines


  subroutine outputContourRegions(xf, cp, height, colormap)
    type(xmlf_t), intent(inout) :: xf
    type(contourObject), intent(in) :: cp
    real(sp), intent(in), optional :: height
    type(color_t), intent(in), optional :: colormap(:)

    integer :: colornum
    integer :: i, j, k
    character(len=8) :: colorhex

    if (present(colormap)) then
      if (size(colormap)<size(cp%contours+1)) then
        call FoX_error("More contours than colours in outputContourRegions")
      endif
    endif

    call kmlOpenFolder(xf, name='Contour Regions')

    do i = 0, size(cp%contours)
      if (i==0) then
        call kmlOpenFolder(xf, name="Less than "//str(cp%contours(1)))
      elseif (i<size(cp%contours)-1) then
        call kmlOpenFolder(xf, name="Between "//str(cp%contours(i))//" and "//str(cp%contours(i+1)))
      elseif (i==size(cp%contours)-1) then
        if (associated(cp%zmax)) then
          call kmlOpenFolder(xf, name="More than "//str(cp%contours(i)))
        else
          call kmlOpenFolder(xf, name="Between "//str(cp%contours(i))//" and "//str(cp%contours(i+1)))
        endif
      elseif (i==size(cp%contours)) then
        if (associated(cp%zmax)) then
          call kmlOpenFolder(xf, name="Region of undefined values")
        else
          call kmlOpenFolder(xf, name="More than "//str(cp%contours(i)))
        endif
      endif

      do j = 1, size(cp%polys(i)%list)
        if (associated(cp%zmax).and.i==size(cp%contours)) then
          colorhex="00000000"
        elseif (present(colormap)) then
          colorhex = kmlGetColorHex(colormap(i))
        else
          colornum = i*255/(size(cp%contours)+1)
          colorhex = "c0"//str(255-colornum, "x2")//"ff"//str(colornum, "x2")
        endif
        ! FIXME here we should use height to adjust height of the polygon, but currently we haven't stored
        ! that information in cp.
        call kmlStartRegion(xf, cp%polys(i)%list(j)%outerBoundary%x, cp%polys(i)%list(j)%outerBoundary%y, &
          name=str(i)//' '//str(j), fillcolorhex=colorhex)
        do k = 1, size(cp%polys(i)%list(j)%innerBoundaries%list)
          call kmlAddInnerBoundary(xf, cp%polys(i)%list(j)%outerBoundary%x, cp%polys(i)%list(j)%outerBoundary%y)
        enddo
        call kmlEndRegion(xf)
      enddo
      call kmlCloseFolder(xf)
    enddo
    call kmlCloseFolder(xf)


  end subroutine outputContourRegions
#endif
end module m_wkml_contours
