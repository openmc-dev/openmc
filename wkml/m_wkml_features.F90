module m_wkml_features

  use fox_m_fsys_realtypes, only : sp, dp
  use FoX_wxml, only: xmlf_t
#ifndef DUMMYLIB
  use m_common_error, only: FoX_error
  use FoX_common, only:  str, operator(//)
  use FoX_wxml, only : xml_AddAttribute, xml_NewElement, xml_EndElement, xml_AddCharacters, xmlf_OpenTag
  use m_wkml_lowlevel, only: kmlOpenFolder, kmlCloseFolder, kmlOpenPlacemark, kmlClosePlacemark, &
    kmlAddCoordinates, kmlOpenInnerBoundaryIs, kmlCloseInnerBoundaryIs, kmlOpenOuterBoundaryIs, kmlCloseOuterBoundaryIs, &
    kmlAddDescription, kmlAddStyleURL, kmlOpenPoint, kmlClosePoint, kmlAddAltitudeMode, kmlAddExtrude, &
    kmlOpenLineString, kmlCloseLineString, kmlOpenLinearRing, kmlCloseLinearRing, kmlOpenPolygon, kmlClosePolygon, &
    kmlAddwhen, kmlOpenTimeStamp, kmlCloseTimeStamp, kmlAddname,kmlAddtessellate
  use m_wkml_styling, only: kmlCreatePointStyle, kmlCreateLineStyle, kmlCreatePolygonStyle, &
    kmlOpenStyle, kmlCloseStyle
  use m_wkml_chart
#endif
  use m_wkml_Color, only: col => color_t

  implicit none
  private

  interface kmlAddPoint
    module procedure kmlAddPoint_dp
    module procedure kmlAddPoint_sp
  end interface kmlAddPoint

  interface kmlCreatePoints
    module procedure kmlCreatePoints_0d_sp
    module procedure kmlCreatePoints_1d_sp
    module procedure kmlCreatePoints_2d_sp
    module procedure kmlCreatePoints_0d_dp
    module procedure kmlCreatePoints_1d_dp
    module procedure kmlCreatePoints_2d_dp
  end interface kmlCreatePoints

  interface kmlCreateLine
    module procedure kmlCreateLine_1d_sp
    module procedure kmlCreateLine_2d_sp
    module procedure kmlCreateLine_1d_dp
    module procedure kmlCreateLine_2d_dp
  end interface kmlCreateLine

  interface kmlStartRegion
    module procedure kmlStartPolygon_1d_sp
    module procedure kmlStartPolygon_2d_sp
    module procedure kmlStartPolygon_1d_dp
    module procedure kmlStartPolygon_2d_dp
  end interface kmlStartRegion

  interface kmlAddInnerBoundary
    module procedure kmlAddInnerBoundary_1d_sp
    module procedure kmlAddInnerBoundary_2d_sp
    module procedure kmlAddInnerBoundary_1d_dp
    module procedure kmlAddInnerBoundary_2d_dp
  end interface kmlAddInnerBoundary

  public :: kmlAddPoint
  public :: kmlCreatePoints
  public :: kmlCreateLine
  public :: kmlStartRegion
  public :: kmlAddInnerBoundary
  public :: kmlEndRegion

contains

! First, routines for adding a series of individual points

! a) Add a single point (single and double precision)
  subroutine kmlAddPoint_sp(xf, longitude, latitude, z, extrude, altitudeMode)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: longitude,latitude
    real(sp), intent(in), optional :: z
    logical, intent(in), optional :: extrude
    character(len=*), intent(in), optional :: altitudeMode
! No need for tessellate - it has no effect on a Point according to KML documentation

#ifndef DUMMYLIB
    call kmlOpenPoint(xf)
    if (present(extrude)) call kmlAddExtrude(xf, extrude)
    if (present(altitudeMode)) call kmlAddAltitudeMode(xf, altitudeMode)
    call kmlAddCoordinates(xf, longitude, latitude, z)
    call kmlClosePoint(xf)
#endif
  end subroutine kmlAddPoint_sp

  subroutine kmlAddPoint_dp(xf, longitude, latitude, z, extrude, altitudeMode)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: longitude,latitude
    real(dp), intent(in), optional :: z
    logical, intent(in), optional :: extrude
    character(len=*), intent(in), optional :: altitudeMode
! No need for tessellate - it has no effect on a Point according to KML documentation

#ifndef DUMMYLIB
    call kmlOpenPoint(xf)
    if (present(extrude)) call kmlAddExtrude(xf, extrude)
    if (present(altitudeMode)) call kmlAddAltitudeMode(xf, altitudeMode)
    call kmlAddCoordinates(xf, longitude, latitude, z)
    call kmlClosePoint(xf)
#endif
  end subroutine kmlAddPoint_dp

! Now - interfaces for a series of points. We have 0d (a single point, but with all options,
! for symmetry), 1d (2 arrays, one for long, one for lat) and 2d (one array, 1st dim 2, holding
! long and lat).
! 0d/2d work by simply rearranging args and immediately calling 1d.

! We do this once for single precision, and again for double precision. This is
! done by cut & pasting but really we should autogenerate through m4.

! We do not also provide assumed-size-array versions of these, since doing this
! by hand involves too much room for error. After autogeneration it would be easy
! though.

! Interface 0: Single long/lat point, optional altitude
  subroutine kmlCreatePoints_0d_sp(xf, longitude, latitude, altitude, &
    extrude, altitudeMode, description,description_numbers,&
    name, color, colorname, colorhex, scale, styleURL,&
    charttype,chartscale,chartsize,chartdata,charttitle,&
    chartlabel,dataname,values)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: longitude
    real(sp), intent(in) :: latitude
    real(sp), intent(in), optional :: altitude
    logical, intent(in), optional :: extrude
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorname
    character(len=8), intent(in), optional :: colorhex
    real(sp), intent(in), optional :: scale ! FIXME what if this is an integer
    character(len=*), intent(in), optional :: description(:)
    real(sp), intent(in), optional :: description_numbers(:)
    character(len=*), intent(in), optional :: styleURL

! variables for kmlAddChart
    character(len=*), intent(in),optional :: charttype,chartscale,chartsize,charttitle,chartlabel
    real(sp), intent(in), optional :: chartdata(:)
! variable for extended data
    character(len=*), intent(in),optional :: dataname
    integer :: k
    real(sp), intent(in), optional :: values(:)

#ifndef DUMMYLIB
    ! Need to check presence explicitly since we mutate altitude with (//)
    if (present(altitude)) then
      call kmlCreatePoints(xf, (/longitude/), (/latitude/), (/altitude/), &
        extrude=extrude, altitudeMode=altitudeMode, &
        name=name, color=color, colorname=colorname, colorhex=colorhex, &
        scale=scale, styleURL=styleURL,description_numbers=description_numbers,&
        description=description,charttype=charttype,&
        chartscale=chartscale,chartsize=chartsize,chartdata=chartdata,&
        charttitle=charttitle,chartlabel=chartlabel,dataname=dataname,values=values)
    else
      call kmlCreatePoints(xf, (/longitude/), (/latitude/), &
        extrude=extrude, altitudeMode=altitudeMode, &
        name=name, color=color, colorname=colorname, colorhex=colorhex, &
        scale=scale,description=description,&
        description_numbers=description_numbers,styleURL=styleURL,&
        charttype=charttype,chartscale=chartscale,chartsize=chartsize,&
        charttitle=charttitle,chartdata=chartdata,chartlabel=chartlabel,&
        dataname=dataname,values=values)
    endif
#endif
  end subroutine kmlCreatePoints_0d_sp

! Interface 1: separate long & lat arrays, optional altitude
  subroutine kmlCreatePoints_1d_sp(xf, longitude, latitude, altitude, &
    extrude, altitudeMode, description,description_numbers,&
    name, color, colorname, colorhex, &
    scale, styleURL,time,&
    charttype,chartscale,chartsize,chartdata,charttitle,chartlabel,&
    dataname, values)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: longitude(:)
    real(sp), intent(in) :: latitude(:)
    real(sp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: extrude
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorname
    character(len=8), intent(in), optional :: colorhex
    real(sp), intent(in), optional :: scale ! FIXME what if this is an integer
    character(len=*), intent(in), optional :: description(:)
    real(sp), intent(in), optional :: description_numbers(:)
    character(len=*), intent(in), optional :: styleURL
    character(len=*), intent(in), optional :: time(:)
    integer :: i, n

! variables for kmlAddChart
    character(len=*), intent(in),optional :: charttype,chartscale,chartsize,charttitle
    character(len=*), intent(in),optional :: chartlabel
    real(sp), intent(in),optional :: chartdata(:)
! variable for extended data
    character(len=*), intent(in),optional :: dataname
    integer :: k
    real(sp), intent(in), optional :: values(:)

#ifndef DUMMYLIB
    n = size(longitude)
    if (n/=size(latitude)) then
      call FoX_error("Incommensurate sizes for longitude and latitude arrays in kmlCreatePoints")
    endif
    if (present(altitude)) then
      if (n/=size(altitude)) then
        call FoX_error("Incommensurate sizes for longitude and altitude arraysin kmlCreatePoints")
      endif
    endif

    if (present(styleURL)) then
      if (present(color).or.present(scale)) then
        call FoX_error("cannot specify styleURL as well as color or scale")
      endif
    endif

    call kmlOpenFolder(xf, name=name)

    do i = 1, n
      call kmlOpenPlacemark(xf)
      if (present(color).or.present(colorname).or.present(colorhex).or.present(scale)) then
        if (present(scale)) then
          call kmlCreatePointStyle(xf, scale=scale, color=color, colorname=colorname, colorhex=colorhex)
        else
          call kmlCreatePointStyle(xf, color=color, colorname=colorname, colorhex=colorhex)
        endif
      elseif (present(styleURL)) then
        call kmlAddStyleURL(xf, styleURL)
      endif

! remove description and will replace by extended data in the future 24042008
      if (present(description).and.present(description_numbers)) then
        call kmlAddDescription(xf, description(i)//achar(13)//str(description_numbers(i)))
      elseif (present(description)) then
        call kmlAddDescription(xf, description(i))
      elseif (present(description_numbers)) then
        call kmlAddDescription(xf, str(description_numbers(i)))
      endif

! add by GT 24042008 for adding chart functions
      if (present(charttype).and.present(chartsize).and.present(chartdata)&
       .and.present(chartscale).and.present(charttitle).and.present(chartlabel)) then
        call kmlAddChart(xf,charttype,chartsize,chartdata,chartscale,&
        charttitle,chartlabel)
      end if

      ! adding time funtion by GT 17102007
      if (present(time)) then
        call kmlOpenTimeStamp(xf)
         call kmlAddwhen(xf,time(i))
        call kmlCloseTimeStamp(xf)
      end if

      ! Extended Data
      if (present(dataname)) then
        call xml_NewElement(xf,'ExtendedData')
        do k=1,size(values)
          call xml_NewElement(xf,'Data')
           call xml_AddAttribute(xf,'name', dataname)
           call xml_NewElement(xf,'displayName')
            call xml_AddCharacters(xf,dataname)
           call xml_EndElement(xf,'displayName')
           call xml_NewElement(xf,'value')
            call xml_AddCharacters(xf,str(values(k),fmt="r5"))
           call xml_EndElement(xf,'value')
          call xml_EndElement(xf,'Data')
        end do
        call xml_EndElement(xf,'ExtendedData')
      end if
        
      if (present(altitude)) then
        call kmlAddPoint(xf, longitude(i), latitude(i), altitude(i), extrude=extrude, altitudeMode=altitudeMode)
      else
        call kmlAddPoint(xf, longitude(i), latitude(i), extrude=extrude, altitudeMode=altitudeMode)
      endif
      call kmlClosePlacemark(xf)
    enddo

    call kmlCloseFolder(xf)
#endif
  end subroutine kmlCreatePoints_1d_sp

! Interface 2: combined long & lat array, optional altitude
  subroutine kmlCreatePoints_2d_sp(xf, coords, altitude, &
    extrude, altitudeMode,description,description_numbers,&
    name, color, colorname, colorhex, scale, styleURL,&
    charttype,chartscale,chartsize,chartdata,charttitle,chartlabel,&
    dataname,values)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: coords(:,:)
    real(sp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: extrude
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorname
    character(len=8), intent(in), optional :: colorhex
    real(sp), intent(in), optional :: scale ! FIXME what if this is an integer
    character(len=*), intent(in), optional :: description(:)
    real(sp), intent(in), optional :: description_numbers(:)
    character(len=*), intent(in), optional :: styleURL

! variables for kmlAddChart
    character(len=*), intent(in),optional :: charttype,chartscale,chartsize,charttitle,chartlabel
    real(sp), intent(in), optional :: chartdata(:)
! variable for extended data
    character(len=*), intent(in),optional :: dataname
    integer :: k
    real(sp), intent(in), optional :: values(:)

#ifndef DUMMYLIB    
    if (size(coords,1)==2) then
      call kmlCreatePoints(xf, coords(1,:), coords(2,:), altitude, &
        extrude=extrude, altitudeMode=altitudeMode, &
        name=name, color=color, scale=scale, styleURL=styleURL,&
        description_numbers=description_numbers,description=description,&
        charttype=charttype,chartscale=chartscale,chartsize=chartsize,&
        chartdata=chartdata,chartlabel=chartlabel,dataname=dataname,values=values)
    elseif (size(coords,1)==3) then
      if (present(altitude)) then
        call FoX_error("Cannot specify 3-dimensional coords with separate altitude in kmlCreatePoints")
      endif
      call kmlCreatePoints(xf, coords(1,:), coords(2,:), coords(3,:), &
        extrude=extrude, altitudeMode=altitudeMode, &
        name=name, color=color, scale=scale,styleURL=styleURL,&
        description_numbers=description_numbers,description=description,&
        charttype=charttype,chartscale=chartscale,chartsize=chartsize,&
        charttitle=charttitle,chartdata=chartdata,chartlabel=chartlabel,&
        dataname=dataname,values=values)
    else
      call FoX_error("coords array first dimension is wrong in kmlCreatePoints - must be 2 or 3")
    endif
#endif
  end subroutine kmlCreatePoints_2d_sp

! Interface 0: Single long/lat point, optional altitude
  subroutine kmlCreatePoints_0d_dp(xf, longitude, latitude, altitude, &
    extrude, altitudeMode,description,description_numbers,&
    name, color, colorname, colorhex, scale, styleURL,&
    charttype,chartscale,chartsize,chartdata,charttitle,&
    chartlabel,dataname,values)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: longitude
    real(dp), intent(in) :: latitude
    real(dp), intent(in), optional :: altitude
    logical, intent(in), optional :: extrude
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorname
    character(len=8), intent(in), optional :: colorhex
    real(dp), intent(in), optional :: scale ! FIXME what if this is an integer
    character(len=*), intent(in), optional :: description(:)
    real(dp), intent(in), optional :: description_numbers(:)
    character(len=*), intent(in), optional :: styleURL

! variables for kmlAddChart
    character(len=*), intent(in),optional :: charttype,chartscale,chartsize,charttitle,chartlabel
    real(dp), intent(in), optional :: chartdata(:)
! variable for extended data
    character(len=*), intent(in),optional :: dataname
    integer :: k
    real(dp), intent(in), optional :: values(:)


#ifndef DUMMYLIB
    ! Need to check presence explicitly since we mutate altitude with (//)
    if (present(altitude)) then
      call kmlCreatePoints(xf, (/longitude/), (/latitude/), (/altitude/), &
        extrude=extrude, altitudeMode=altitudeMode, &
        name=name, color=color, colorname=colorname, colorhex=colorhex, &
        scale=scale, styleURL=styleURL,charttype=charttype,&
        description_numbers=description_numbers,description=description,&
        chartscale=chartscale,chartsize=chartsize,chartdata=chartdata,&
        charttitle=charttitle,chartlabel=chartlabel,dataname=dataname, values=values)
    else
      call kmlCreatePoints(xf, (/longitude/), (/latitude/), &
        extrude=extrude, altitudeMode=altitudeMode, &
        name=name, color=color, colorname=colorname, colorhex=colorhex, &
        scale=scale, styleURL=styleURL,charttype=charttype,&
        description_numbers=description_numbers,description=description,&
        chartscale=chartscale,chartsize=chartsize,chartdata=chartdata,&
        charttitle=charttitle,chartlabel=chartlabel,dataname=dataname, values=values)
    endif
#endif
  end subroutine kmlCreatePoints_0d_dp

! Interface 1: separate long & lat arrays, optional altitude
  subroutine kmlCreatePoints_1d_dp(xf, longitude, latitude, altitude, &
    extrude, altitudeMode, description_numbers,description,&
    name, color, colorname, colorhex, scale, styleURL, time,&
    charttype,chartscale,chartsize,chartdata,charttitle,chartlabel,&
    dataname, values)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: longitude(:)
    real(dp), intent(in) :: latitude(:)
    real(dp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: extrude
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorname
    character(len=8), intent(in), optional :: colorhex
    real(dp), intent(in), optional :: scale ! FIXME what if this is an integer
    character(len=*), intent(in), optional :: description(:) ! FIXME need to check length of array if present
    real(dp), intent(in), optional :: description_numbers(:)
    character(len=*), intent(in), optional :: styleURL
    character(len=*), intent(in), optional :: time(:)    

    integer :: i, n
    character, pointer :: tempStyleUrl(:)

! variables for kmlAddChart
    character(len=*), intent(in),optional :: charttype,chartscale,chartsize,charttitle,chartlabel
    real(dp), intent(in), optional :: chartdata(:)
! variable for extended data
    character(len=*), intent(in),optional :: dataname
    integer :: k
    real(dp), intent(in), optional :: values(:)

#ifndef DUMMYLIB
    n = size(longitude)
    if (n/=size(latitude)) then
      call FoX_error("Incommensurate sizes for longitude and latitude arrays in kmlCreatePoints")
    endif
    if (present(altitude)) then
      if (n/=size(altitude)) then
        call FoX_error("Incommensurate sizes for longitude and altitude arraysin kmlCreatePoints")
      endif
    endif

    if (present(styleURL)) then
      if (present(color).or.present(scale)) then
        call FoX_error("cannot specify styleURL as well as color or scale")
      endif
    endif

    call kmlOpenFolder(xf, name=name)

    do i = 1, n
      call kmlOpenPlacemark(xf)
      if (present(color).or.present(colorname).or.present(colorhex).or.present(scale)) then
        if (present(scale)) then
          call kmlCreatePointStyle(xf, scale=scale, color=color, colorname=colorname, colorhex=colorhex)
        else
          call kmlCreatePointStyle(xf, color=color, colorname=colorname, colorhex=colorhex)
        endif
      elseif (present(styleURL)) then
        call kmlAddStyleURL(xf, styleURL)
      endif
! remove description and use extended data instead 24042008
      if (present(description).and.present(description_numbers)) then
        call kmlAddDescription(xf, description(i)//achar(13)//str(description_numbers(i)))
      elseif (present(description)) then
        call kmlAddDescription(xf, description(i))
      elseif (present(description_numbers)) then
        call kmlAddDescription(xf, str(description_numbers(i)))
      endif

! add by GT 24042008 for adding chart functions
      if (present(charttype).and.present(chartsize).and.  & 
          present(chartdata).and.present(chartscale).and. & 
          present(charttitle) ) then
        call kmlAddChart(xf,charttype,chartsize,chartdata,chartscale,charttitle,chartlabel)
      end if

    ! adding time funtion by GT 17102007
      if (present(time)) then
        call kmlOpenTimeStamp(xf)
         call kmlAddwhen(xf,time(i))
        call kmlCloseTimeStamp(xf)
      end if

    ! add by GT for extended data 24/04/2008
      if (present(dataname)) then
        call xml_NewElement(xf,'ExtendedData')
        do k=1,size(values)
          call xml_NewElement(xf,'Data')
           call xml_AddAttribute(xf,'name', dataname)
           call xml_NewElement(xf,'displayName')
            call xml_AddCharacters(xf,dataname)
           call xml_EndElement(xf,'displayName')
           call xml_NewElement(xf,'value')
            call xml_AddCharacters(xf,str(values(k),fmt="r5"))
           call xml_EndElement(xf,'value')
          call xml_EndElement(xf,'Data')
        end do
        call xml_EndElement(xf,'ExtendedData')
      end if

      if (present(altitude)) then
        call kmlAddPoint(xf, longitude(i), latitude(i), altitude(i), extrude=extrude, altitudeMode=altitudeMode)
      else
        call kmlAddPoint(xf, longitude(i), latitude(i), extrude=extrude, altitudeMode=altitudeMode)
      endif
      call kmlClosePlacemark(xf)
    enddo

    call kmlCloseFolder(xf)
#endif
  end subroutine kmlCreatePoints_1d_dp

! Interface 2: combined long & lat array, optional altitude
  subroutine kmlCreatePoints_2d_dp(xf, coords, altitude, &
    extrude, altitudeMode, &
    name, color, colorname, colorhex, &
    description, description_numbers,scale, styleURL,&
    charttype,chartscale,chartsize,chartdata,charttitle,&
    chartlabel,dataname, values)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: extrude
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorname
    character(len=8), intent(in), optional :: colorhex
    real(dp), intent(in), optional :: scale ! FIXME what if this is an integer
    character(len=*), intent(in), optional :: description(:)
    real(dp), intent(in), optional :: description_numbers(:)
    character(len=*), intent(in), optional :: styleURL
!    character(len=*), intent(in), optional :: time    

! variables for kmlAddChart
    character(len=*), intent(in),optional :: charttype,chartscale,chartsize,charttitle,chartlabel
    real(dp), intent(in), optional :: chartdata(:)
! variable for extended data
    character(len=*), intent(in),optional :: dataname
    integer :: k
    real(dp), intent(in), optional :: values(:)

#ifndef DUMMYLIB
    if (size(coords,1)==2) then
      call kmlCreatePoints(xf, coords(1,:), coords(2,:), altitude, &
        extrude=extrude, altitudeMode=altitudeMode, &
        name=name, color=color, colorname=colorname, colorhex=colorhex, &
        scale=scale,styleURL=styleURL,charttype=charttype,&
        description_numbers=description_numbers,description=description,&
        chartscale=chartscale,chartsize=chartsize,chartdata=chartdata,&
        charttitle=charttitle,chartlabel=chartlabel,dataname=dataname, values=values)
    elseif (size(coords,1)==3) then
      if (present(altitude)) then
        call FoX_error("Cannot specify 3-dimensional coords with separate altitude in kmlCreatePoints")
      endif
      call kmlCreatePoints(xf, coords(1,:), coords(2,:), coords(3,:), &
        extrude=extrude, altitudeMode=altitudeMode, &
        name=name, color=color, colorname=colorname, colorhex=colorhex, &
        scale=scale,styleURL=styleURL,charttype=charttype, &
        description_numbers=description_numbers,description=description,&
        chartscale=chartscale,chartsize=chartsize,chartdata=chartdata,&
        charttitle=charttitle,chartlabel=chartlabel,&
        dataname=dataname, values=values)
    else
      call FoX_error("coords array first dimension is wrong in kmlCreatePoints - must be 2 or 3")
    endif
#endif
  end subroutine kmlCreatePoints_2d_dp


! Now commands to create lines - only one interface really (array (x), array(y), optional(z))
! but specified for both single and double precision, assumed-size and assumed-shape arrays



! all those kmlCreateLine subroutines assumes that long and lat has to be the same size
! this function is the easiest case, wchih just allows to create a line segment
! lat and long is a sigle value  07/03/2008 GT

  subroutine kmlCreateLine_seg_sh_dp(xf,xi,yi,xe,ye,zi,ze,tessellate,altmode,&
                                     name,linewidth,description_ch,styleURL)
    type(xmlf_t), intent(inout) ::xf
    real(dp), intent(in) :: xi,yi,xe,ye ! start x coor and end coord
    character(len=*),intent(in), optional :: name, styleURL,altmode,linewidth,description_ch
    logical, intent(in),optional :: tessellate
    real(dp), intent(in),optional :: zi,ze
#ifndef DUMMYLIB
    integer :: i, nodes
    character(len=1) :: palm

    palm="#"
    call kmlOpenPlacemark(xf)
    if(present(name)) then
      call kmlAddname(xf,name)
    else
      call kmlAddname(xf,'')
    end if
    if (present(description_ch)) then
      call kmlAdddescription(xf, description_ch)
    end if

    if (present(styleURL)) then
      call kmlAddstyleUrl(xf,palm//styleURL)
    else
      call kmlAddstyleUrl(xf,palm//styleURL)
    end if
    call kmlOpenLineString(xf)
    if(present(tessellate)) then
      call kmlAddtessellate(xf,tessellate)
    else
      call kmlAddtessellate(xf,.true.)
    end if
    if(present(altmode)) then
      call kmlAddaltitudeMode(xf,altmode)
    else
      call kmlAddaltitudeMode(xf,"relativeToGround")
    end if
    call xml_NewElement(xf,'coordinates')
     call xml_AddCharacters(xf,xi)
     call xml_AddCharacters(xf,',')
     call xml_AddCharacters(xf,yi)
     if (present(zi)) then
       call xml_AddCharacters(xf,',')
       call xml_AddCharacters(xf,zi)
     end if
     call xml_AddCharacters(xf,' ')
     call xml_AddCharacters(xf,xe)
     call xml_AddCharacters(xf,',')
     call xml_AddCharacters(xf,ye)
     if (present(ze)) then
       call xml_AddCharacters(xf,',')
       call xml_AddCharacters(xf,ze)
     end if
     call xml_AddCharacters(xf,' ')  !adding space
    call xml_EndElement(xf,'coordinates')
    call kmlCloseLineString(xf)
    call kmlClosePlacemark(xf)
#endif
  end subroutine kmlCreateLine_seg_sh_dp


! Interface 1: separate long & lat arrays, optional altitude
  subroutine kmlCreateLine_1d_sp(xf, longitude, latitude, altitude, &
    closed, extrude, tessellate, altitudeMode, &
    name, color, colorname, colorhex, &
    width, description, styleURL)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: longitude(:)
    real(sp), intent(in) :: latitude(:)
    real(sp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: closed
    logical, intent(in), optional :: extrude
    logical, intent(in), optional :: tessellate
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorname
    character(len=8), intent(in), optional :: colorhex
    integer, intent(in), optional :: width
    character(len=*), intent(in), optional :: description
    character(len=*), intent(in), optional :: styleURL

#ifndef DUMMYLIB    
    integer :: n
    logical :: closed_, needPlacemark

    if (present(closed)) then
      closed_ = closed
    else
      closed_ = .false.
    endif
    needPlacemark = (xmlf_OpenTag(xf)/="Placemark".and.xmlf_OpenTag(xf)/="outerBoundaryIs".and.xmlf_OpenTag(xf)/="innerBoundaryIs")

    if (needPlacemark) call kmlOpenPlacemark(xf, name=name, description=description, styleurl=styleurl)
    n = size(longitude)
    if (closed_.and.n<3) then
      call FoX_error("Not enough points on closed path")
    endif

    n = size(longitude)
    if (n/=size(latitude)) then
      call FoX_error("Incommensurate sizes for longitude and latitude arrays in kmlCreateLine")
    endif
    if (present(altitude)) then
      if (n/=size(altitude)) then
        call FoX_error("Incommensurate sizes for longitude and altitude arrays in kmlCreateLine")
      endif
    endif

    if (present(styleURL)) then
      if (present(color).or.present(width)) then
        call FoX_error("cannot specify styleURL as well as color or width")
      endif
    endif

    if (closed_) then
      call kmlOpenLinearRing(xf, altitudeMode=altitudeMode, tessellate=tessellate, extrude=extrude)
    else
      call kmlOpenLineString(xf, altitudeMode=altitudeMode, tessellate=tessellate, extrude=extrude)
    endif

    if (present(color).or.present(colorname).or.present(colorhex).or.present(width)) then
      call kmlCreateLineStyle(xf, color=color, colorname=colorname, colorhex=colorhex, width=width)
    elseif (present(styleURL)) then
      call kmlAddStyleURL(xf, styleURL)
    endif

    if (present(description)) call kmlAddDescription(xf, description)
    call kmlAddCoordinates(xf, longitude, latitude, altitude, repeat=closed_)

    if (closed_) then
      call kmlCloseLinearRing(xf)
    else
      call kmlCloseLineString(xf)
    endif
    if (needPlacemark) call kmlClosePlacemark(xf)
#endif
  end subroutine kmlCreateLine_1d_sp

! Interface 2: combined long & lat array, optional altitude
  subroutine kmlCreateLine_2d_sp(xf, coords, altitude, &
    closed, extrude, tessellate, altitudeMode, &
    name, color, colorname, colorhex, &
    width, description, styleURL)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: coords(:,:)
    real(sp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: closed
    logical, intent(in), optional :: extrude
    logical, intent(in), optional :: tessellate
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorname
    character(len=8), intent(in), optional :: colorhex
    integer, intent(in), optional :: width
    character(len=*), intent(in), optional :: description
    character(len=*), intent(in), optional :: styleURL

#ifndef DUMMYLIB    
    if (size(coords,1)==2) then
      call kmlCreateLine(xf, coords(1,:), coords(2,:), altitude, &
        closed=closed, extrude=extrude, tessellate=tessellate, altitudeMode=altitudeMode, &
        name=name, color=color, colorname=colorname, colorhex=colorhex, &
        width=width, description=description, styleURL=styleURL)
    elseif (size(coords,1)==3) then
      if (present(altitude)) then
        call FoX_error("Cannot specify 3-dimensional coords with separate altitude in kmlCreateLine")
      endif
      call kmlCreateLine(xf, coords(1,:), coords(2,:), coords(3,:), &
        closed=closed, extrude=extrude, tessellate=tessellate, altitudeMode=altitudeMode, &
        name=name, color=color, colorname=colorname, colorhex=colorhex, &
        width=width, description=description, styleURL=styleURL)
    else
      call FoX_error("coords array first dimension is wrong in kmlCreateLine - must be 2 or 3")
    endif
#endif
  end subroutine kmlCreateLine_2d_sp

! Interface 1: separate long & lat arrays, optional altitude
  subroutine kmlCreateLine_1d_dp(xf, longitude, latitude, altitude, &
    closed, extrude, tessellate, altitudeMode, &
    name, color, colorname, colorhex, &
    width, description, styleURL)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: longitude(:)
    real(dp), intent(in) :: latitude(:)
    real(dp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: closed
    logical, intent(in), optional :: extrude
    logical, intent(in), optional :: tessellate
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorname
    character(len=8), intent(in), optional :: colorhex
    integer, intent(in), optional :: width
    character(len=*), intent(in), optional :: description
    character(len=*), intent(in), optional :: styleURL

#ifndef DUMMYLIB    
    integer :: n
    logical :: closed_, needPlacemark

    if (present(closed)) then
      closed_ = closed
    else
      closed_ = .false.
    endif
    needPlacemark = (xmlf_OpenTag(xf)/="Placemark".and.xmlf_OpenTag(xf)/="outerBoundaryIs".and.xmlf_OpenTag(xf)/="innerBoundaryIs")

    if (needPlacemark) call kmlOpenPlacemark(xf, name=name, description=description, styleurl=styleurl)
    n = size(longitude)
    if (closed_.and.n<3) then
      call FoX_error("Not enough points on closed path")
    endif

    n = size(longitude)
    if (n/=size(latitude)) then
      call FoX_error("Incommensurate sizes for longitude and latitude arrays in kmlCreateLine")
    endif
    if (present(altitude)) then
      if (n/=size(altitude)) then
        call FoX_error("Incommensurate sizes for longitude and altitude arrays in kmlCreateLine")
      endif
    endif

    if (present(styleURL)) then
      if (present(color).or.present(width)) then
        call FoX_error("cannot specify styleURL as well as color or width")
      endif
    endif

    if (closed_) then
      call kmlOpenLinearRing(xf, altitudeMode=altitudeMode, tessellate=tessellate, extrude=extrude)
    else
      call kmlOpenLineString(xf, altitudeMode=altitudeMode, tessellate=tessellate, extrude=extrude)
    endif

    if (present(color).or.present(colorname).or.present(colorhex).or.present(width)) then
      call kmlCreateLineStyle(xf, color=color, colorname=colorname, colorhex=colorhex, width=width)
    elseif (present(styleURL)) then
      call kmlAddStyleURL(xf, styleURL)
    endif

    if (present(description)) call kmlAddDescription(xf, description)
    call kmlAddCoordinates(xf, longitude, latitude, altitude, repeat=closed_)

    if (closed_) then
      call kmlCloseLinearRing(xf)
    else
      call kmlCloseLineString(xf)
    endif
    if (needPlacemark) call kmlClosePlacemark(xf)
#endif
  end subroutine kmlCreateLine_1d_dp

! Interface 2: combined long & lat array, optional altitude
  subroutine kmlCreateLine_2d_dp(xf, coords, altitude, &
    closed, extrude, tessellate, altitudeMode, &
    name, color, colorname, colorhex, &
    width, description, styleURL)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: closed
    logical, intent(in), optional :: extrude
    logical, intent(in), optional :: tessellate
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorname
    character(len=8), intent(in), optional :: colorhex
    integer, intent(in), optional :: width
    character(len=*), intent(in), optional :: description
    character(len=*), intent(in), optional :: styleURL

#ifndef DUMMYLIB    
    if (size(coords,1)==2) then
      call kmlCreateLine(xf, coords(1,:), coords(2,:), altitude, &
        closed=closed, extrude=extrude, tessellate=tessellate, altitudeMode=altitudeMode, &
        name=name, color=color, colorname=colorname, colorhex=colorhex, &
        width=width, description=description, styleURL=styleURL)
    elseif (size(coords,1)==3) then
      if (present(altitude)) then
        call FoX_error("Cannot specify 3-dimensional coords with separate altitude in kmlCreateLine")
      endif
      call kmlCreateLine(xf, coords(1,:), coords(2,:), coords(3,:), &
        closed=closed, extrude=extrude, tessellate=tessellate, altitudeMode=altitudeMode, &
        name=name, color=color, colorname=colorname, colorhex=colorhex, &
        width=width, description=description, styleURL=styleURL)
    else
      call FoX_error("coords array first dimension is wrong in kmlCreateLine - must be 2 or 3")
    endif
#endif
  end subroutine kmlCreateLine_2d_dp


  subroutine kmlStartPolygon_1d_sp(xf, longitude, latitude, altitude, &
    extrude, tessellate, altitudeMode, &
    name, fillcolor, fillcolorname, fillcolorhex, &
    linecolor, linecolorname, linecolorhex, linewidth, description, styleURL)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: longitude(:)
    real(sp), intent(in) :: latitude(:)
    real(sp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: extrude
    logical, intent(in), optional :: tessellate
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: fillcolor
    character(len=*), intent(in), optional :: fillcolorname
    character(len=8), intent(in), optional :: fillcolorhex
    type(col), intent(in), optional :: linecolor
    character(len=*), intent(in), optional :: linecolorname
    character(len=8), intent(in), optional :: linecolorhex
    integer, intent(in), optional :: linewidth
    character(len=*), intent(in), optional :: description
    character(len=*), intent(in), optional :: styleURL

#ifndef DUMMYLIB    
    integer :: n
    logical :: needPlacemark, outline, fill

    n = size(longitude)
    outline = present(linecolor).or.present(linecolorname).or.present(linecolorhex).or.present(linewidth)
    fill = present(fillcolor).or.present(fillcolorname).or.present(fillcolorhex)

    call kmlOpenPlacemark(xf, name=name, description=description)
    if (fill.or.outline) then
      call kmlOpenStyle(xf)
      if (outline) then
        call kmlCreateLineStyle(xf, color=linecolor, colorname=linecolorname, colorhex=linecolorhex, &
          width=linewidth)
      endif
      if (fill) then
        call kmlCreatePolygonStyle(xf, color=fillcolor, colorname=fillcolorname, colorhex=fillcolorhex, &
          fill=fill, outline=outline)
      endif
      call kmlCloseStyle(xf)
    elseif (present(styleUrl)) then
      call kmlAddStyleUrl(xf, styleUrl)
    endif

    call kmlOpenPolygon(xf, extrude=extrude, tessellate=tessellate, altitudeMode=altitudeMode)
    call kmlOpenOuterBoundaryIs(xf)
    call kmlCreateLine(xf, longitude, latitude, altitude, closed=.true.)
    call kmlCloseOuterBoundaryIs(xf)

    ! Leave polygon unclosed, we might add innerboundaries ...
#endif
  end subroutine kmlStartPolygon_1d_sp

  subroutine kmlStartPolygon_2d_sp(xf, coords, altitude, &
    extrude, tessellate, altitudeMode, &
    name, fillcolor, fillcolorname, fillcolorhex, &
    linecolor, linecolorname, linecolorhex, linewidth, description, styleURL)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: coords(:,:)
    real(sp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: extrude
    logical, intent(in), optional :: tessellate
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: fillcolor
    character(len=*), intent(in), optional :: fillcolorname
    character(len=8), intent(in), optional :: fillcolorhex
    type(col), intent(in), optional :: linecolor
    character(len=*), intent(in), optional :: linecolorname
    character(len=8), intent(in), optional :: linecolorhex
    integer, intent(in), optional :: linewidth
    character(len=*), intent(in), optional :: description
    character(len=*), intent(in), optional :: styleURL

#ifndef DUMMYLIB    
    if (size(coords,1)==2) then
      call kmlStartRegion(xf, coords(1,:), coords(2,:), altitude, &
        extrude=extrude, tessellate=tessellate, altitudeMode=altitudeMode, &
        name=name, fillcolor=fillcolor, fillcolorname=fillcolorname, fillcolorhex=fillcolorhex, &
        linecolor=linecolor, linecolorname=linecolorname, linecolorhex=linecolorhex, &
        linewidth=linewidth, description=description, styleURL=styleURL)
    elseif (size(coords,1)==3) then
      if (present(altitude)) then
        call FoX_error("Cannot specify 3-dimensional coords with separate altitude in kmlStartPolygon")
      endif
      call kmlStartRegion(xf, coords(1,:), coords(2,:), coords(3,:), &
        extrude=extrude, tessellate=tessellate, altitudeMode=altitudeMode, &
        name=name, fillcolor=fillcolor, fillcolorname=fillcolorname, fillcolorhex=fillcolorhex, &
        linecolor=linecolor, linecolorname=linecolorname, linecolorhex=linecolorhex, &
        linewidth=linewidth, description=description, styleURL=styleURL)
    else
      call FoX_error("coords array first dimension is wrong in kmlStartPolygon - must be 2 or 3")
    endif
#endif
  end subroutine kmlStartPolygon_2d_sp

  subroutine kmlStartPolygon_1d_dp(xf, longitude, latitude, altitude, &
    extrude, tessellate, altitudeMode, &
    name, fillcolor, fillcolorname, fillcolorhex, &
    linecolor, linecolorname, linecolorhex, linewidth, description, styleURL)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: longitude(:)
    real(dp), intent(in) :: latitude(:)
    real(dp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: extrude
    logical, intent(in), optional :: tessellate
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: fillcolor
    character(len=*), intent(in), optional :: fillcolorname
    character(len=8), intent(in), optional :: fillcolorhex
    type(col), intent(in), optional :: linecolor
    character(len=*), intent(in), optional :: linecolorname
    character(len=8), intent(in), optional :: linecolorhex
    integer, intent(in), optional :: linewidth
    character(len=*), intent(in), optional :: description
    character(len=*), intent(in), optional :: styleURL

#ifndef DUMMYLIB    
    integer :: n
    logical :: needPlacemark, outline, fill

    n = size(longitude)
    outline = present(linecolor).or.present(linecolorname).or.present(linecolorhex).or.present(linewidth)
    fill = present(fillcolor).or.present(fillcolorname).or.present(fillcolorhex)

    call kmlOpenPlacemark(xf, name=name, description=description)
    if (fill.or.outline) then
      call kmlOpenStyle(xf)
      if (outline) then
        call kmlCreateLineStyle(xf, color=linecolor, colorname=linecolorname, colorhex=linecolorhex, &
          width=linewidth)
      endif
      if (fill) then
        call kmlCreatePolygonStyle(xf, color=fillcolor, colorname=fillcolorname, colorhex=fillcolorhex, &
          fill=fill, outline=outline)
      endif
      call kmlCloseStyle(xf)
    elseif (present(styleUrl)) then
      call kmlAddStyleUrl(xf, styleUrl)
    endif

    call kmlOpenPolygon(xf, extrude=extrude, tessellate=tessellate, altitudeMode=altitudeMode)

    call kmlOpenOuterBoundaryIs(xf)
    call kmlCreateLine(xf, longitude, latitude, altitude, closed=.true.)
    call kmlCloseOuterBoundaryIs(xf)

    ! Leave polygon unclosed, we might add innerboundaries ...
#endif
  end subroutine kmlStartPolygon_1d_dp

  subroutine kmlStartPolygon_2d_dp(xf, coords, altitude, extrude, &
                                   tessellate, altitudeMode, name,&
                                   fillcolor, fillcolorname, fillcolorhex, &
                                   linecolor, linecolorname, linecolorhex, &
                                   linewidth, description, styleURL)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in), optional :: altitude(:)
    logical, intent(in), optional :: extrude
    logical, intent(in), optional :: tessellate
    character(len=*), intent(in), optional :: altitudeMode
    character(len=*), intent(in), optional :: name
    type(col), intent(in), optional :: fillcolor
    character(len=*), intent(in), optional :: fillcolorname
    character(len=8), intent(in), optional :: fillcolorhex
    type(col), intent(in), optional :: linecolor
    character(len=*), intent(in), optional :: linecolorname
    character(len=8), intent(in), optional :: linecolorhex
    integer, intent(in), optional :: linewidth
    character(len=*), intent(in), optional :: description
    character(len=*), intent(in), optional :: styleURL

#ifndef DUMMYLIB    
    if (size(coords,1)==2) then
      call kmlStartRegion(xf, coords(1,:), coords(2,:), altitude, &
        extrude=extrude, tessellate=tessellate, altitudeMode=altitudeMode, &
        name=name, fillcolor=fillcolor, fillcolorname=fillcolorname, fillcolorhex=fillcolorhex, &
        linecolor=linecolor, linecolorname=linecolorname, linecolorhex=linecolorhex, &
        linewidth=linewidth, description=description, styleURL=styleURL)
    elseif (size(coords,1)==3) then
      if (present(altitude)) then
        call FoX_error("Cannot specify 3-dimensional coords with separate altitude in kmlStartPolygon")
      endif
      call kmlStartRegion(xf, coords(1,:), coords(2,:), coords(3,:), &
        extrude=extrude, tessellate=tessellate, altitudeMode=altitudeMode, &
        name=name, fillcolor=fillcolor, fillcolorname=fillcolorname, fillcolorhex=fillcolorhex, &
        linecolor=linecolor, linecolorname=linecolorname, linecolorhex=linecolorhex, &
        linewidth=linewidth, description=description, styleURL=styleURL)
    else
      call FoX_error("coords array first dimension is wrong in kmlStartPolygon - must be 2 or 3")
    endif
#endif
  end subroutine kmlStartPolygon_2d_dp

  subroutine kmlAddInnerBoundary_1d_sp(xf, longitude, latitude, altitude)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: longitude(:)
    real(sp), intent(in) :: latitude(:)
    real(sp), intent(in), optional :: altitude(:)

#ifndef DUMMYLIB
    integer :: n

    n = size(longitude)

    if (xmlf_OpenTag(xf)/="Polygon") then
      call FoX_error("Can only add an inner boundary inside a polygon")
    endif

    call kmlOpenInnerBoundaryIs(xf)
    call kmlCreateLine(xf, longitude, latitude, altitude, closed=.true.)
    call kmlCloseInnerBoundaryIs(xf)
#endif
  end subroutine kmlAddInnerBoundary_1d_sp

  subroutine kmlAddInnerBoundary_2d_sp(xf, coords, altitude)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: coords(:,:)
    real(sp), intent(in), optional :: altitude(:)

#ifndef DUMMYLIB    
    if (size(coords,1)==2) then
      call kmlAddInnerBoundary(xf, coords(1,:), coords(2,:), altitude)
    elseif (size(coords,1)==3) then
      if (present(altitude)) then
        call FoX_error("Cannot specify 3-dimensional coords with separate altitude in kmlAddInnerBoundary")
      endif
      call kmlAddInnerBoundary(xf, coords(1,:), coords(2,:), coords(3,:))
    else
      call FoX_error("coords array first dimension is wrong in kmlAddInnerBoundary - must be 2 or 3")
    endif
#endif
  end subroutine kmlAddInnerBoundary_2d_sp

  subroutine kmlAddInnerBoundary_1d_dp(xf, longitude, latitude, altitude)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: longitude(:)
    real(dp), intent(in) :: latitude(:)
    real(dp), intent(in), optional :: altitude(:)

#ifndef DUMMYLIB    
    integer :: n

    n = size(longitude)

    if (xmlf_OpenTag(xf)/="Polygon") then
      call FoX_error("Can only add an inner boundary inside a polygon")
    endif

    call kmlOpenInnerBoundaryIs(xf)
    call kmlCreateLine(xf, longitude, latitude, altitude, closed=.true.)
    call kmlCloseInnerBoundaryIs(xf)
#endif
  end subroutine kmlAddInnerBoundary_1d_dp

  subroutine kmlAddInnerBoundary_2d_dp(xf, coords, altitude)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: coords(:,:)
    real(dp), intent(in), optional :: altitude(:)

#ifndef DUMMYLIB
    if (size(coords,1)==2) then
      call kmlAddInnerBoundary(xf, coords(1,:), coords(2,:), altitude)
    elseif (size(coords,1)==3) then
      if (present(altitude)) then
        call FoX_error("Cannot specify 3-dimensional coords with separate altitude in kmlAddInnerBoundary")
      endif
      call kmlAddInnerBoundary(xf, coords(1,:), coords(2,:), coords(3,:))
    else
      call FoX_error("coords array first dimension is wrong in kmlAddInnerBoundary - must be 2 or 3")
    endif
#endif
  end subroutine kmlAddInnerBoundary_2d_dp
 
  subroutine kmlEndRegion(xf)
    type(xmlf_t), intent(inout) :: xf

#ifndef DUMMYLIB
    call kmlClosePolygon(xf)
    call kmlClosePlacemark(xf)
#endif
  end subroutine kmlEndRegion

end module m_wkml_features
