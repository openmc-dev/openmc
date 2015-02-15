define(`TOHW_m4_wkml_createCells3',`dnl
  subroutine kmlcreateCells3_$1(xf,longitude,latitude,values,myCI,mask,outline,time,vizvalues,dataname)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! this subroutine is going to read X, Y ,Z, stylecolor
      ! each XYZ is a vector, this is used for testing glimmer or netcdf situation 01302007 GT

      type(xmlf_t),     intent(inout)        :: xf 
      real($1),         intent(in)           :: longitude(:)
      real($1),         intent(in)           :: latitude(:)
      real($1),         intent(in)           :: values(:,:)
      type(color_t),    intent(in)           :: myCI(:) 
      real($1),         intent(in), optional :: mask ! usually represent no data
      logical,          intent(in), optional :: outline
      character(len=*), intent(in), optional :: time
      real($1),         intent(in), optional :: vizvalues(:)
      character(len=*), intent(in), optional :: dataname
#ifndef DUMMYLIB

      integer :: i, j, k, x, y
      integer :: nx, ny, nnx, nny  ! numbers at X(long), numbers at Y(Lat)

      type(color_t), allocatable :: valuehex(:,:)

      integer,dimension(4) :: xp=(/0,1,1,0/)  ! id for coordintes
      integer,dimension(4) :: yp=(/0,0,1,1/)

      character(LEN=8) :: stylecolor, colorhextmp

      character(15) :: lonchar, latchar, elchar
      character(50) :: coords
      real($1) :: valueres

!      if (present(valuescale)) then
!      values=valuescale*values
!      end if

! get the size of x and y vector
      nx=size(longitude)
      ny=size(latitude)
!      nnx=size(values,1)
!      nny=size(values,2)


!        do i=1,nx-1
!         write(*,400), lon(i)
!        end do
! allocate the memory for x and y
!     allocate(lon(nx))
!     allocate(lat(ny))
! allocate the memory for values
     allocate(valuehex(nx,ny))
!    if (present(vizvalues)) then
      do k=1, size(vizvalues)-1
        do i=1,nx
         do j=1,ny
          if (values(i,j) .lt. vizvalues(1)) then
          valuehex(i,j)= myCI(1)
          end if
          if ((values(i,j) >= vizvalues(k)) .and. (values(i,j) < vizvalues(k+1))) then
          valuehex(i,j)= myCI(k+1)
          end if
          if (values(i,j) .GT. vizvalues(k+1)) then
          valuehex(i,j)= myCI(k+2)
          end if
         end do
        end do
      end do

!     else
!     dividing passed in values to how many colors scales
!       valueres=(MAXVAL(values,MASK = values .LT. mask)-MINVAL(values))/size(myCI)
!       do k=1, size(myCI)
!        do i=1,nx
!         do j=1,ny
!         if (values(i,j) >= MINVAL(values)+valueres*(k-1)) then
!             valuehex(i,j)= myCI(k) !sometime this line is not used
!         end if
!         end do
!        end do
!       end do

!     end if
! adding style function in 071307 GT
!       do i=1,size(myCI)
!         call kmlCreatePolygonStyle(xf,color=myCI(i),id=str(i))
!       end do


      do i=1,nx-1
           do j=1, ny-1
!          if(all(values(i:i+1,j:j+1)==mask)) cycle
          if (values(i,j) == mask) cycle
          call kmlOpenPlacemark(xf)
           call kmlAddname(xf,"srf_dep")

           ! adding time funtion by GT 10042008
           if (present(time)) then
            call kmlOpenTimeStamp(xf)
            call kmlAddwhen(xf,time)
            call kmlCloseTimeStamp(xf)
           end if
!          add by GT for extended data 21/04/2008
           if (present(dataname)) then
           call xml_NewElement(xf,"ExtendedData")
            call xml_NewElement(xf,"Data")
               call xml_AddAttribute(xf,"name", dataname)
               call xml_NewElement(xf,"displayName")
                 call xml_AddCharacters(xf,dataname)
               call xml_EndElement(xf,"displayName")
               call xml_NewElement(xf,"value")
                 call xml_AddCharacters(xf,str(values(i,j),fmt="r5"))
               call xml_EndElement(xf,"value")
            call xml_EndElement(xf,"Data")
            call xml_EndElement(xf,"ExtendedData")
            end if
!           call kmlAddstyleUrl(xf,"#"//stylecolor)
            if (present(outline)) then
            call kmlCreatePolygonStyle(xf,color=valuehex(i,j),outline=outline)
            else
            call kmlCreatePolygonStyle(xf,color=valuehex(i,j))
            end if
!            call kmlAddstyleUrl(xf,"#"//valuehex(i,j))

           call kmlOpenPolygon(xf)
             call kmlAddextrude(xf,.true.)
             call kmlAddaltitudeMode(xf,"clampToGround")
           call kmlOpenouterBoundaryIs(xf)
             call kmlOpenLinearRing(xf)
                call xml_NewElement(xf,name="coordinates")
                      coords=str(longitude(i))//","//str(latitude(j))//","//str(values(i,j))
                      call xml_AddCharacters(xf,coords)
                      call xml_AddNewLine(xf) ! this function is missing in FOX2.0.2 version
                      coords=str(longitude(i))//","//str(latitude(j+1))//","//str(values(i,j))
                      call xml_AddCharacters(xf,coords)
                      call xml_AddNewLine(xf) ! this function is missing in FOX2.0.2 version
                      coords=str(longitude(i+1))//","//str(latitude(j+1))//","//str(values(i,j))
                      call xml_AddCharacters(xf,coords)
                      call xml_AddNewLine(xf) ! this function is missing in FOX2.0.2 version
                      coords=str(longitude(i+1))//","//str(latitude(j))//","//str(values(i,j))
                      call xml_AddCharacters(xf,coords)
                      call xml_AddNewLine(xf) ! this function is missing in FOX2.0.2 version
               call xml_EndElement(xf,name="coordinates")
             call kmlCloseLinearRing(xf)
           call kmlCloseouterBoundaryIs(xf)
           call kmlClosePolygon(xf)
          call kmlClosePlacemark(xf)
          end do
       end do
!      deallocate(longitude)
!      deallocate(latitude)
!      deallocate(values)
       deallocate(valuehex)
#endif
  end subroutine kmlCreateCells3_$1')`'dnl
dnl
define(`TOHW_m4_wkml_createCells'`',`dnl
  subroutine $1_$2(xf, &
$3, values, &
    mask, colormap, height, contour_values, num_levels, name)
    type(xmlf_t) :: xf
    real($2), intent(in) :: $4
    real($2), intent(in) :: values(:,:)
    real($2), intent(in), optional :: mask
    type(color_t), target, optional :: colormap(:)
    real($2), intent(in), optional :: height
    real($2), intent(in), optional :: contour_values(:)
    integer, intent(in), optional :: num_levels
    character(len=*), intent(in), optional :: name
#ifndef DUMMYLIB
    integer  :: i, ic, j, k, m, n, numcolors
    real($2) :: square(3,4), lat, long, average
    real($2) :: minvalue, lat_inc, long_inc, valueres !resolution of input value
    character(len=15), allocatable :: styleURL(:) ! FIXME this ought to be dynamically sized,
                                     ! but this allows up to 10^9 separate IDs in one doc.
    type(color_t), pointer :: defaultMap(:), thisColor

    m = size(values, 1)
    n = size(values, 2)

    if (present(contour_values).and.present(num_levels)) then
      call FoX_error("Cannot specify both contour_values and num_levels in kmlCreateCells")
    elseif (present(contour_values)) then
      if (present(colormap)) then
        if (size(colormap)/=size(contour_values)+1) then
          call FoX_error("colormap must be one item longer than contour_values in kmlCreateCells")
        endif
      endif
      numcolors = size(contour_values)+1
    elseif (present(num_levels)) then
      if (present(colormap)) then
        if (size(colormap)/=num_levels+1) then
          call FoX_error("colormap must be one item longer than num_levels in kmlCreateCells")
        endif
      endif
      numcolors = num_levels+1
    else
      if (present(colormap)) then
        numcolors = size(colormap)
      else
        numcolors = 5
      endif
    endif

    if (.not.present(colormap)) defaultMap => kmlMakeColorMap(numcolors)
    
    minvalue = minval(values)
    if (present(mask)) then
      valueres = (maxval(values, mask=(values<mask))-minvalue)/(numcolors-1)
    else
      valueres = (maxval(values)-minvalue)/(numcolors-1)
    endif

    call kmlOpenFolder(xf, name=name)

! When we have a working style-handling policy, replace here.
!!$    allocate(styleURL(numcolors))
!!$    do i = 1, numcolors
!!$      styleURL(i) = xmlf_NewId(xf)
!!$      if (present(colormap)) then
!!$        call kmlCreatePolygonStyle(xf, color=colormap(i), id=trim(styleURL(i)))
!!$      else
!!$        call kmlCreatePolygonStyle(xf, color=defaultMap(i), id=trim(styleURL(i)))
!!$      endif
!!$    end do

$5

    call kmlCloseFolder(xf)

    if (.not.present(colormap)) deallocate(defaultMap)
#endif
  end subroutine $1_$2
')`'dnl
dnl
define(`TOHW_m4_wkml_spdp_createCells', `'`dnl
TOHW_m4_wkml_createCells(`$1', `sp', `$2', `$3', `$4')`'dnl
TOHW_m4_wkml_createCells(`$1', `dp', `$2', `$3', `$4')`'dnl
')`'dnl
dnl
module m_wkml_coverage

  use fox_m_fsys_realtypes, only: sp, dp
  use FoX_wxml, only : xmlf_t
  use m_wkml_color, only: color_t
#ifndef DUMMYLIB
  use FoX_wxml, only : xml_NewElement, xml_EndElement, xml_AddAttribute, xml_AddNewLine, xml_AddCharacters
  use m_common_error, only: FoX_error
  use FoX_common, only: str
  use m_common_error

  use m_wkml_lowlevel, only: kmlOpenFolder, kmlCloseFolder, kmlopenplacemark, & 
       kmlAddname, kmlAddstyleurl, kmlopenpolygon, kmladdextrude, kmladdaltitudemode, & 
       kmlopenouterboundaryis, kmlopenlinearring, kmlcloselinearring, & 
       kmlcloseouterboundaryis, kmlclosepolygon, kmlcloseplacemark, &
       kmlOpenTimeStamp, kmlCloseTimeStamp, kmlAddwhen
  use m_wkml_color, only: kmlSetCustomColor, kmlMakeColorMap
  use m_wkml_features, only: kmlStartRegion, kmlEndRegion
  use m_wkml_styling, only: kmlCreatePolygonStyle
  use m_wkml_chart
#endif

  implicit none
  private

  interface kmlCreateRGBCells
    module procedure kmlCreateRGBCells_sp
    module procedure kmlCreateRGBCells_dp
  end interface kmlCreateRGBCells

  interface kmlCreateCells
    module procedure kmlCreateCells_sp
    module procedure kmlCreateCells_dp
    module procedure kmlCreateCells_longlat_sp
    module procedure kmlCreateCells_longlat_dp
    module procedure kmlCreateCells_longlat2_sp
    module procedure kmlCreateCells_longlat2_dp
    module procedure kmlCreateCells3_dp
    module procedure kmlCreateCells3_sp
  end interface kmlCreateCells

  public :: kmlCreateRGBCells
  public :: kmlCreateCells

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !this subroutine is used for creating cells (pixels),this version is used for SIG and assing rgb color auto

  subroutine kmlCreateRGBCells_sp(xf, east, west, south, north, reflectance, rgb, numbit)
    type(xmlf_t) :: xf
    real(sp), intent(in)  :: east, west, south, north
    real(sp), intent(in), optional :: reflectance(:,:)
    integer, intent(in), optional :: numbit !! color interval
    character, intent(in), optional :: rgb
#ifndef DUMMYLIB
    call kmlCreateRGBCells(xf, real(east, dp), real(west, dp), &
      real(south, dp), real(north, dp), real(reflectance, dp), rgb, numbit)
#endif
  end subroutine kmlCreateRGBCells_sp

  subroutine kmlCreateRGBCells_dp(xf, east, west, south, north, reflectance, rgb, numbit)
    type(xmlf_t) :: xf
    real(dp), intent(in)  :: east, west, south, north
    real(dp), intent(in), optional :: reflectance(:,:)
    integer, intent(in), optional :: numbit !! color interval
    character, intent(in), optional :: rgb
#ifndef DUMMYLIB
    integer :: numbit_
    integer :: i, dn

    type(color_t), allocatable :: colormap(:)

    if (rgb/="r" .and. rgb/="g" .and. rgb/="b") then
      call FoX_error("Must use one of r, g, or b in CreateRGBCells")
    endif

    if (present(numbit)) then
      if (numbit<1.or.numbit>256) then
        call FoX_error("numbit out of range")
      elseif (mod(256, numbit)/=0) then
        call FoX_error("numbit must be a power of 2")
      endif
      numbit_ = numbit ! Check for sensible values
    else
      numbit_ = 256
    endif

    allocate(colormap(numbit_))
    do i = 1, numbit_
      dn = 256/i - 1
      if (rgb=="b") then
        call kmlSetCustomColor(colormap(i), "EE"//str(dn, "x")//"00"//"00")
      elseif (rgb=="g") then
        call kmlSetCustomColor(colormap(i), "EE"//"00"//str(dn, "x")//"00")
      elseif (rgb=="r") then
        call kmlSetCustomColor(colormap(i), "EE"//"00"//"00"//str(dn, "x"))
      endif
    enddo

    call kmlCreateCells(xf, east=east, west=west, south=south, north=north, &
      values=reflectance, mask=1.0d0, colormap=colormap)
#endif
  end subroutine kmlCreateRGBCells_dp

! createCells was called createCells2/createCells3
! Its function is to produce coloured cells, with the colour varying
! according to value.
! value must be passed in as a 2D array.
! x/y coords may be specified either by:
! 1: simple EWSN coords to specify the corners of the array
!     ie, a very simple rectangular spaced grid with equal spacing of points.
! 2: two 1-D arrays, one of longitude, one of latitude (of the same lengths as value(:,:))
!     ie a rectangular grid, but with irregular spacing
! 3: two 2-D arrays, one for longitude, one for latitude (both of the same size as value(:,:))
!     ie a topologically rectilinear, but otherwise irregular grid.

! FIXME in the simplest case (1: above) we should make sure and document where the EWNS are
! expected to be in relationship to the grid points. Currently we assume coincident.

TOHW_m4_wkml_spdp_createCells(`kmlCreateCells', `east, west, south, north', `east, west, south, north', `'`dnl
    lat_inc = (east-west)/m ! Increment in latitude
    long_inc = (north-south)/n ! Increment in longitude
    do i = 1, m
      long = west+long_inc*(i-0.5) ! Subtract 0.5 so that cells are *centred* on the long/lat point.
      do j = 1, n
        ! Plot from north to south in order to match data structure
        lat = north+lat_inc*(i-0.5)
        if (present(mask)) then
          if (values(i,j)>=mask) cycle
        endif
        square(1, :) = (/long, long, long+long_inc, long+long_inc/)       ! x-coords
        square(2, :) = (/lat, lat-lat_inc, lat-lat_inc, lat/)             ! y-coords
        if (present(height)) then                                           ! z-coords
          square(3,:) = height*((/values(i,j), values(i+1,j), values(i+1,j+1), values(i+1,j+1)/)-minValue)
        endif
        if (present(contour_values)) then
          ! New logic by GT: this version works, 07032009
          do ic=1, size(contour_values)-1
            if (values(i,j) .lt. contour_values(1)) then
              k = ic
            else if ((values(i,j) >= contour_values(ic)) &
                .and. (values(i,j) < contour_values(ic+1))) then
              k = ic+1
            end if
          end do
        else
          k = int(floor((values(i, j)-minvalue)/valueres))
        endif
        if (present(colormap)) then
          thisColor => colormap(k+1)
        else
          thisColor => defaultMap(k+1)
        endif
        if (present(colormap)) then
          if (present(height)) then
            call kmlStartRegion(xf, square, &
              extrude=.true., altitudeMode="relativeToGround", fillcolor=colorMap(k+1))
          else
            call kmlStartRegion(xf, square(:2,:), fillcolor=colorMap(k+1))
          endif
        else
          if (present(height)) then
            call kmlStartRegion(xf, square, &
              extrude=.true., altitudeMode="relativeToGround", fillcolor=defaultMap(k+1))
          else
            call kmlStartRegion(xf, square(:2,:), fillcolor=defaultMap(k+1))
          endif
        endif
        call kmlEndRegion(xf)
      end do
    end do
')`'


TOHW_m4_wkml_spdp_createCells(`kmlCreateCells_longlat', `longitude, latitude', `longitude(:), latitude(:)', `'`dnl
    ! To help with .m4 and make the autogenerated code like GTs code. Rewrite valueres here...
    ! but only for longlat arguments - not sure why!
    if (.not.present(mask)) then
       valueres = maxval(values)/(numcolors-1)
    endif
    do i = 1, m-1
      do j = 1, n-1
        if (present(mask)) then
          if (any(values(i:i+1, j:j+1)>=mask)) cycle ! Dont draw the cell if any of its vertices are masked out
        endif
        square(1, :) = (/longitude(i), longitude(i+1), longitude(i+1), longitude(i)/) ! x-coords
        square(2, :) = (/latitude(j), latitude(j), latitude(j+1), latitude(j+1)/)     ! y-coords
        if (present(height)) then                                           ! z-coords
          square(3,:) = height*((/values(i,j), values(i+1,j), values(i+1,j+1), values(i+1,j+1)/)-minValue)
        endif
        average = sum(values(i:i+1,j:j+1))/4.0d0
        ! Colour the cell according to the average of the 4 values defining the cell.
        if (present(contour_values)) then
           do ic=1, size(contour_values)-1
              if (average .lt. contour_values(1)) then  
                 k = 0
              else if ((average >= contour_values(ic)) .and. (average < contour_values(ic+1))) then
                 k = ic
              else if ((average >= contour_values(ic))) then
                 k = ic+1
              end if
           end do
        else
          k = int(floor((average-minvalue)/valueres))
        endif
        if (present(colormap)) then
          thisColor => colormap(k+1)
        else
          thisColor => defaultMap(k+1)
        endif
        if (present(colormap)) then
          if (present(height)) then
            call kmlStartRegion(xf, square, &
              extrude=.true., altitudeMode="relativeToGround", &
              fillcolor=colorMap(k+1), description=str(values(i,j),FMT="r5"))
          else
            call kmlStartRegion(xf, square(:2,:), &
              fillcolor=colorMap(k+1), description=str(values(i,j),FMT="r5"))
          endif
        else
          if (present(height)) then
            call kmlStartRegion(xf, square, &
              extrude=.true., altitudeMode="relativeToGround", &
              fillcolor=defaultMap(k+1), description=str(values(i,j),FMT="r5"))
          else
            call kmlStartRegion(xf, square(:2,:), &
              fillcolor=defaultMap(k+1), description=str(values(i,j),FMT="r5"))
          endif
        endif
        call kmlEndRegion(xf)
      end do
    end do
')`'

TOHW_m4_wkml_spdp_createCells(`kmlCreateCells_longlat2', `longitude, latitude', `longitude(:,:), latitude(:,:)', `'`dnl
    do i = 1, m-1
      do j = 1, n-1
        if (present(mask)) then
          if (any(values(i:i+1, j:j+1)>=mask)) cycle ! Dont draw the cell if any of its vertices are masked out
        endif
        square(1, :) = (/longitude(i,j), longitude(i+1,j), longitude(i+1,j+1), longitude(i,j+1)/) ! x-coords
        square(2, :) = (/latitude(i,j), latitude(i+1,j), latitude(i+1,j+1), latitude(i,j+1)/)     ! y-coords
        if (present(height)) then                                                                         ! z-coords
          square(3,:) = height*((/values(i,j), values(i+1,j), values(i+1,j+1), values(i+1,j+1)/)-minValue)
        endif
        average = sum(values(i:i+1,j:j+1))/4.0d0
        ! Colour the cell according to the average of the 4 values defining the cell.
        if (present(contour_values)) then
          do ic=1, size(contour_values)-1    
            if (average .lt. contour_values(1)) then 
              k = 0                              
            else if ((average >= contour_values(ic)) .and. (average < contour_values(ic+1))) then 
              k = ic+1 
            end if      
          end do                         
        else
          k = int(floor((average-minvalue)/valueres))
        endif
        if (present(colormap)) then
          if (present(height)) then
            call kmlStartRegion(xf, square, &
              extrude=.true., altitudeMode="relativeToGround", fillcolor=colorMap(k+1))
          else
            call kmlStartRegion(xf, square(:2,:), fillcolor=colorMap(k+1))
          endif
        else
          if (present(height)) then
            call kmlStartRegion(xf, square, &
              extrude=.true., altitudeMode="relativeToGround", fillcolor=defaultMap(k+1))
          else
            call kmlStartRegion(xf, square(:2,:), fillcolor=defaultMap(k+1))
          endif
        endif
        call kmlEndRegion(xf)
      end do
    end do
')`'

TOHW_m4_wkml_createCells3(`dp')

TOHW_m4_wkml_createCells3(`sp')

end module m_wkml_coverage
