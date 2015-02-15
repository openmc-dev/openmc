module m_wkml_styling

  use fox_m_fsys_realtypes, only : sp, dp
  use FoX_wxml, only: xmlf_t
  use m_wkml_color, only: col => color_t
#ifndef DUMMYLIB
  use m_common_error
  use FoX_wxml, only: xmlf_OpenTag, &
    xml_NewElement, xml_EndElement, xml_AddAttribute
  use m_wkml_color, only: kmlAddColor, kmlAddBgColor, kmlAddTextColor, kmlGetCustomColor
  use m_wkml_lowlevel, only: kmlAddFill, kmlAddWidth, kmlAddColorMode, kmlAddHref, &
    kmlAddkey, kmlAddOutline, kmlAddScale, kmlAddStyleUrl, kmlAddHeading, &
    kmlOpenIcon, kmlCloseIcon, kmlOpenPair, kmlClosePair, kmlAddwidth,kmlAddcolorMode, &
    kmlOpenPlacemark,kmlClosePlacemark,kmlAddname,kmlOpenPolygon,kmlClosePolygon
#endif
 
  implicit none
  private

  interface kmlCreatePointStyle
    module procedure kmlAddPointStyle_s_none_hd_none
    module procedure kmlAddPointStyle_s_int
    module procedure kmlAddPointStyle_s_sp
    module procedure kmlAddPointStyle_s_dp
    module procedure kmlAddPointStyle_s_none_hd_dp
    module procedure kmlAddPointStyle_s_none_hd_int
    module procedure kmlAddPointStyle_s_none_hd_sp
  end interface kmlCreatePointStyle

  public :: kmlCreatePointStyle
  public :: kmlCreateLineStyle
  public :: kmlCreatePolygonStyle
#ifndef DUMMYLIB
  public :: kmlOpenStyle, kmlCloseStyle

  ! add by GT 18/04/2008
  public :: kmlAddLegend
#endif
contains
#ifndef DUMMYLIB
  subroutine mostOfPointStyle(xf, color, colorhex, colorname, colormode, heading, iconhref, id)
    type(xmlf_t), intent(inout) :: xf
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode
    real(dp), intent(in), optional :: heading
    character(len=*), intent(in), optional :: iconhref
    character(len=*), intent(in), optional :: id

    logical :: needStyle

    needStyle = (xmlf_openTag(xf)/='Style')
    if (needStyle) then
      call kmlOpenStyle(xf, id)
      call kmlOpenIconStyle(xf)
    else
      call kmlOpenIconStyle(xf, id)
    endif

    if (count((/present(color), present(colorhex), present(colorname)/))>1) then
      print*, "Cannot specify more than one of color, colorhex, and colorname"
    endif
    if (present(color)) then
      call kmlAddColor(xf, color)
    elseif (present(colorhex)) then
      call kmlAddColor(xf, colorhex)
    elseif (present(colorname)) then
      call kmlAddColor(xf, kmlGetCustomColor(colorname))
    endif

    if (present(colormode)) call kmlAddColorMode(xf, colormode)
    if (present(heading)) call kmlAddHeading(xf, heading)
    if (present(iconhref)) then
      call kmlOpenIcon(xf)
      call kmlAddHref(xf, iconhref)
      call kmlCloseIcon(xf)
    endif
  end subroutine mostOfPointStyle
#endif

  subroutine kmlAddPointStyle_s_none_hd_none(xf, color, colorhex, colorname, colormode, iconhref, id)
    ! We temporarily ignore hotspot
    type(xmlf_t), intent(inout) :: xf
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode
    character(len=*), intent(in), optional :: iconhref
    character(len=*), intent(in), optional :: id
#ifndef DUMMYLIB
    logical :: needStyle
    needStyle = (xmlf_openTag(xf)/='Style')

    call mostOfPointStyle(xf, color, colorhex, colorname, colormode, iconhref=iconhref, id=id)
    call kmlCloseIconStyle(xf)
    if (needStyle) call kmlCloseStyle(xf)
#endif
  end subroutine kmlAddPointStyle_s_none_hd_none

  subroutine kmlAddPointStyle_s_none_hd_int(xf, color, colorhex, colorname, colormode, heading, iconhref, id)
    ! We temporarily ignore hotspot
    type(xmlf_t), intent(inout) :: xf
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode
    integer, intent(in) :: heading
    character(len=*), intent(in), optional :: iconhref
    character(len=*), intent(in), optional :: id
#ifndef DUMMYLIB
    logical :: needStyle
    needStyle = (xmlf_openTag(xf)/='Style')

    call mostOfPointStyle(xf, color, colorhex, colorname, colormode, real(heading,dp), iconhref, id)
    call kmlCloseIconStyle(xf)
    if (needStyle) call kmlCloseStyle(xf)
#endif
  end subroutine kmlAddPointStyle_s_none_hd_int

  subroutine kmlAddPointStyle_s_none_hd_sp(xf, color, colorhex, colorname, colormode, heading, iconhref, id)
    ! We temporarily ignore hotspot
    type(xmlf_t), intent(inout) :: xf
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode
    real(sp), intent(in) :: heading
    character(len=*), intent(in), optional :: iconhref
    character(len=*), intent(in), optional :: id
#ifndef DUMMYLIB
    logical :: needStyle
    needStyle = (xmlf_openTag(xf)/='Style')

    call mostOfPointStyle(xf, color, colorhex, colorname, colormode, real(heading,dp), iconhref, id)
    call kmlCloseIconStyle(xf)
    if (needStyle) call kmlCloseStyle(xf)
#endif
  end subroutine kmlAddPointStyle_s_none_hd_sp

  subroutine kmlAddPointStyle_s_none_hd_dp(xf, color, colorhex, colorname, colormode, heading, iconhref, id)
    ! We temporarily ignore hotspot
    type(xmlf_t), intent(inout) :: xf
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode
    real(dp), intent(in) :: heading
    character(len=*), intent(in), optional :: iconhref
    character(len=*), intent(in), optional :: id
#ifndef DUMMYLIB
    logical :: needStyle
    needStyle = (xmlf_openTag(xf)/='Style')

    call mostOfPointStyle(xf, color, colorhex, colorname, colormode, heading, iconhref, id)
    call kmlCloseIconStyle(xf)
    if (needStyle) call kmlCloseStyle(xf)
#endif
  end subroutine kmlAddPointStyle_s_none_hd_dp

  subroutine kmlAddPointStyle_s_int(xf, scale, color, colorhex, colorname, colormode, heading, iconhref, id)
    ! We temporarily ignore hotspot
    type(xmlf_t), intent(inout) :: xf
    integer, intent(in) :: scale
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode
    integer, intent(in), optional :: heading
    character(len=*), intent(in), optional :: iconhref
    character(len=*), intent(in), optional :: id
#ifndef DUMMYLIB
    logical :: needStyle
    needStyle = (xmlf_openTag(xf)/='Style')

    if (present(heading)) then
      call mostOfPointStyle(xf, color, colorhex, colorname, colormode, real(heading,dp), iconhref, id)
    else
      call mostOfPointStyle(xf, color, colorhex, colorname, colormode, iconhref=iconhref, id=id)
    endif
    call kmlAddScale(xf, scale)
    call kmlCloseIconStyle(xf)
    if (needStyle) call kmlCloseStyle(xf)
#endif
  end subroutine kmlAddPointStyle_s_int

  subroutine kmlAddPointStyle_s_sp(xf, scale, color, colorhex, colorname, colormode, heading, iconhref, id)
    ! We temporarily ignore hotspot
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: scale
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode
    real(sp), intent(in), optional :: heading
    character(len=*), intent(in), optional :: iconhref
    character(len=*), intent(in), optional :: id
#ifndef DUMMYLIB
    logical :: needStyle
    needStyle = (xmlf_openTag(xf)/='Style')

    if (present(heading)) then
      call mostOfPointStyle(xf, color, colorhex, colorname, colormode, real(heading,dp), iconhref, id)
    else
      call mostOfPointStyle(xf, color, colorhex, colorname, colormode, iconhref=iconhref, id=id)
    endif
    call kmlAddScale(xf, scale)
    call kmlCloseIconStyle(xf)
    if (needStyle) call kmlCloseStyle(xf)
#endif
  end subroutine kmlAddPointStyle_s_sp

  subroutine kmlAddPointStyle_s_dp(xf, scale, color, colorhex, colorname, colormode, heading, iconhref, id)
    ! We temporarily ignore hotspot
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: scale
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode
    real(dp), intent(in), optional :: heading
    character(len=*), intent(in), optional :: iconhref
    character(len=*), intent(in), optional :: id
#ifndef DUMMYLIB
    logical :: needStyle
    needStyle = (xmlf_openTag(xf)/='Style')

    call mostOfPointStyle(xf, color, colorhex, colorname, colormode, heading, iconhref, id)
    call kmlAddScale(xf, scale)
    call kmlCloseIconStyle(xf)
    if (needStyle) call kmlCloseStyle(xf)
#endif
  end subroutine kmlAddPointStyle_s_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef DUMMYLIB
  subroutine kmlCreatePointStyleOld(xf,normalurl,normalscale,stylemap,highlighturl,highlightscale,&
    normalcolor,highlightcolor,color,shape)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: normalurl, highlighturl,normalcolor
    character(len=*), intent(in), optional ::  highlightcolor
    character(len=*), intent(in), optional :: shape, stylemap  
    character(len=8), intent(in),optional :: color
    character(len=1) :: palm
    real, intent(in),optional :: normalscale,highlightscale

    palm="#"

    !         print*, 'call kmlCreatePointStyle'

    if (present(color)) then
      if (present(normalcolor)) then
        call kmlOpenStyle(xf,normalcolor)
      else
        call kmlOpenStyle(xf,"inactive-"//stylemap)
      end if
      call kmlOpenIconStyle(xf)

      !          print*, 'call kmlCreatePointStyle addcolor'
      call kmlAddcolor(xf,color)
      call kmlAddcolorMode(xf,"normal")

      call kmlOpenIcon(xf)
      if (present(normalurl)) then
        call kmlAddhref(xf,normalurl)
      else
        if (present(shape)) then
          ! print for debug
          !            print*, 'call kmlCreatePointStyle shape'
          call kmlAddhref(xf,"http://cete.niees.group.cam.ac.uk/googlemaptest/icons/"//shape//".png")
        else
          !              print*, 'call kmlCreatePointStyle circle'
          call kmlAddhref(xf,"http://cete.niees.group.cam.ac.uk/googlemaptest/icons/circle.png")             
        end if
      end if
      call kmlCloseIcon(xf)
      call kmlCloseIconStyle(xf)
      call kmlOpenLabelStyle(xf)
      !          print*, 'call OpenLabelStyle'  
      if (present(normalscale)) then
        call kmlAddscale(xf,normalscale)
      else
        call kmlAddscale(xf,0.2)
      end if
      call kmlCloseLabelStyle(xf)
      !          print*, 'call CloseLabelStyle'
      call kmlCloseStyle(xf)

      if (present(highlightcolor)) then

        call kmlOpenStyle(xf,highlightcolor)
        !         print*, 'call kmlOpenStyle highlightcolor'
      else
        call kmlOpenStyle(xf,"active-"//stylemap)
        !         print*, 'call kmlOpenStyle active'
      end if
      call kmlOpenIconStyle(xf)
      !            print*, 'call kmlOpenIconStyle'
      if (present(highlightcolor)) then
        call kmlAddcolor(xf,highlightcolor)
      else
        call kmlAddcolor(xf,"FF00FFFF")
      end if
      call kmlAddcolorMode(xf,"normal")

      call kmlOpenIcon(xf)
      if (present(highlighturl)) then
        call kmlAddhref(xf,highlighturl)
      else
        if (present(shape)) then
          !               print*, 'call if present shape'
          call kmlAddhref(xf,"http://cete.niees.group.cam.ac.uk/googlemaptest/icons/"//shape//".png")
        else
          call kmlAddhref(xf,"http://cete.niees.group.cam.ac.uk/googlemaptest/icons/circle.png")
        end if
      end if
      call kmlCloseIcon(xf)
      call kmlCloseIconStyle(xf)
      call kmlOpenLabelStyle(xf)
      if (present(highlightscale)) then
        call kmlAddscale(xf,highlightscale)
      else
        call kmlAddscale(xf,0.5)
      end if
      call kmlCloseLabelStyle(xf)
      call kmlCloseStyle(xf)

      if (present(stylemap)) then
        call kmlOpenStyleMap(xf,stylemap)
      else
        call kmlOpenStyleMap(xf,"point")
      end if
      call kmlOpenPair(xf)
      call kmlAddkey(xf,"normal")
      if (present(normalcolor)) then
        call kmlAddstyleUrl(xf,palm//normalcolor)
      else
        call kmlAddstyleUrl(xf,palm//"inactive-"//stylemap)
      end if
      call kmlClosePair(xf)
      call kmlOpenPair(xf)
      call kmlAddkey(xf,"highlight")
      if (present(highlightcolor)) then
        call kmlAddstyleUrl(xf,palm//highlightcolor)
      else
        call kmlAddstyleUrl(xf,palm//"active-"//stylemap)
      end if
      call kmlClosePair(xf)
      call kmlCloseStyleMap(xf)


    else !if user did not giving color, then using default
      call kmlOpenStyle(xf,"inactive-"//stylemap)
      !          print*, 'start kmlOpenstyle'
      call kmlOpenIconStyle(xf)
      call kmlAddcolor(xf,"FF0000FF") ! default is red
      call kmlAddcolorMode(xf,"normal")
      call kmlOpenIcon(xf)
      !               print*, 'start kmlOpenIcon'
      if (present(shape)) then
        call kmlAddhref(xf,"http://cete.niees.group.cam.ac.uk/googlemaptest/icons/"//shape//".png")
      else
        call kmlAddhref(xf,"http://cete.niees.group.cam.ac.uk/googlemaptest/icons/circle.png")
      end if
      call kmlCloseIcon(xf)
      call kmlCloseIconStyle(xf)
      call kmlOpenLabelStyle(xf)
      !          print*, 'start kmlOpenLabelStyle'
      if (present(normalscale)) then
        call kmlAddscale(xf,normalscale)
      else
        call kmlAddscale(xf,0.2)
      end if
      call kmlCloseLabelStyle(xf)
      call kmlCloseStyle(xf)
      call kmlOpenStyle(xf,"active-"//stylemap)
      call kmlOpenIconStyle(xf)
      call kmlAddcolor(xf,"FF00FFFF") ! default is yellow
      call kmlAddcolorMode(xf,"normal")
      call kmlOpenIcon(xf)
      if (present(shape)) then
        call kmlAddhref(xf,"http://cete.niees.group.cam.ac.uk/googlemaptest/icons/"//shape//".png")
      else
        call kmlAddhref(xf,"http://cete.niees.group.cam.ac.uk/googlemaptest/icons/circle.png")
      end if
      call kmlCloseIcon(xf)
      call kmlCloseIconStyle(xf)
      call kmlOpenLabelStyle(xf)
      if (present(highlightscale)) then
        call kmlAddscale(xf,highlightscale)
      else
        call kmlAddscale(xf,0.5)
      end if
      call kmlCloseLabelStyle(xf)
      call kmlCloseStyle(xf)

      if (present(stylemap)) then
        call kmlOpenStyleMap(xf,stylemap)
      else
        call kmlOpenStyleMap(xf,"point")
      end if
      call kmlOpenPair(xf)
      call kmlAddkey(xf,"normal")
      call kmlAddstyleUrl(xf,palm//"inactive-"//stylemap)
      call kmlClosePair(xf)
      call kmlOpenPair(xf)
      call kmlAddkey(xf,"highlight")
      call kmlAddstyleUrl(xf,palm//"active-"//stylemap)
      call kmlClosePair(xf)
      call kmlCloseStyleMap(xf)

    end if

  end subroutine kmlCreatePointStyleOld

  ! Here we have four almost identical implementations; we need them
  ! so that we can have width as an effectively optional argument that may be 
  ! integer or real. (we can't do that with an actual optional argument and
  ! preserve the keyword form). All of these are presented through a single
  ! interface above.

  subroutine mostOfLineStyle(xf, color, colorhex, colorname, colormode, id)
    type(xmlf_t), intent(inout) :: xf
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode, id

    logical :: needStyle

    needStyle = (xmlf_openTag(xf)/='Style')
    if (needStyle) then
      call kmlOpenStyle(xf, id)
      call kmlOpenLineStyle(xf)
    else
      call kmlOpenLineStyle(xf, id)
    endif

    if (count((/present(color), present(colorhex), present(colorname)/))>1) then
      print*, "Cannot specify more than one of color, colorhex, and colorname"
    endif
    if (present(color)) then
      call kmlAddColor(xf, color)
    elseif (present(colorhex)) then
      call kmlAddColor(xf, colorhex)
    elseif (present(colorname)) then
      call kmlAddColor(xf, kmlGetCustomColor(colorname))
    endif

    if (present(colormode)) call kmlAddColorMode(xf, colormode)

  end subroutine mostOfLineStyle
#endif

  subroutine kmlCreateLineStyle(xf, width, color, colorhex, colorname, colormode, id)
    type(xmlf_t), intent(inout) :: xf
    integer, intent(in), optional :: width
    type(col), intent(in) , optional:: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode, id
#ifndef DUMMYLIB
    logical :: needStyle
    needStyle = (xmlf_openTag(xf)/='Style')

    call mostOfLineStyle(xf, color, colorhex, colorname, colormode, id)
    if (present(width)) call kmlAddWidth(xf, width)
    call kmlCloseLineStyle(xf)
    if (needStyle) call kmlCloseStyle(xf)
#endif
  end subroutine kmlCreateLineStyle

#ifndef DUMMYLIB
!!!! add the old version by GT 08/03/2008
  subroutine kmlCreateLineStyle_old(xf,stylemap,color,width)
          type(xmlf_t), intent(inout) :: xf
         character(len=*), intent(in) :: stylemap
         character(len=*), intent(in), optional ::  color
         integer, intent(in),optional :: width
         character(len=1) :: palm
         !color has to be: blue, green,yellow,purple,red,orange

         palm="#"
         call kmlOpenStyle(xf,stylemap)
         if(present(color)) then
            if(present(width)) then
              call kmlAddLineStyle(xf, color=color,width=width)
            else
              call kmlAddLineStyle(xf, color=color)
            end if
         else
            if(present(width)) then
              call kmlAddLineStyle(xf, color="EE0000FF",width=width)
             else
              call kmlAddLineStyle(xf, color="EE0000FF")
             end if
          end if
          call kmlCloseStyle(xf)
          end subroutine kmlCreateLineStyle_old
!!! add the old version by GT 08/03/2008
         subroutine kmlAddLineStyle(xf, color, width,colormode,id)
         type(xmlf_t), intent(inout) :: xf
         character(len=*), intent(in) :: color
         character(len=*), intent(in), optional :: colormode, id
         integer, intent(in), optional :: width
         if (present(id)) call kmlOpenStyle(xf,id)
          call kmlOpenLineStyle(xf)
            call kmlAddcolor(xf,color)
            if (present(colormode)) then
            call kmlAddcolorMode(xf,colormode)
            else
            call kmlAddcolorMode(xf,"normal")
            end if
            if  (present(width)) then
            call kmlAddwidth(xf,width)
            else
            call kmlAddwidth(xf,2)
            end if
          call kmlCloseLineStyle(xf)
         end subroutine kmlAddLineStyle


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddLabelStyle_scale(xf, scale1,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    real, intent(in) :: scale1
    !         character(len=*), intent(in) :: scale1
    call xml_NewElement(xf,'LabelStyle')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    call kmlAddscale(xf,scale1)
    call xml_EndElement(xf,'LabelStyle')
  end subroutine kmlAddLabelStyle_scale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddLabelStyle_color(xf, scale1, color, colormode,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    real, intent(in) :: scale1
    character(len=*), intent(in) :: color,colormode
    call xml_NewElement(xf,'LabelStyle')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    call kmlAddcolor(xf,color)
    call kmlAddcolorMode(xf,colormode)
    call kmlAddscale(xf,scale1)
    call xml_EndElement(xf,'LabelStyle')
  end subroutine kmlAddLabelStyle_color


  subroutine kmlAddBalloonStyle(xf, bgcolor, textcolor,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: bgcolor, textcolor
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'BolloonStyle')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    call kmlAddbgcolor(xf, bgcolor)
    call kmlAddtextcolor(xf, textcolor)
    call xml_EndElement(xf, 'BalloonStyle')
  end subroutine kmlAddBalloonStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Specifies how a Feature is displayed in the list view. The list view is a hierarchy of containers and children; in Google Earth, this is the Places panel.
  subroutine kmlAddListStyle_bgcolor(xf, bgcolor)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: bgcolor
    call xml_NewElement(xf, 'ListStyle')
    call kmlAddbgcolor(xf, bgcolor)
    call xml_EndElement(xf, 'ListStyle')
  end subroutine kmlAddListStyle_bgcolor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenListStyle(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'ListStyle')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenListStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseListStyle(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'ListStyle')
  end subroutine kmlCloseListStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenLabelStyle(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'LabelStyle')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenLabelStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseLabelStyle(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'LabelStyle')
  end subroutine kmlCloseLabelStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenBalloonStyle(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'BalloonStyle')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenBalloonStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseBalloonStyle(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'BalloonStyle')
  end subroutine kmlCloseBalloonStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenColorStyle(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'ColorStyle')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenColorStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseColorStyle(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'ColorStyle')
  end subroutine kmlCloseColorStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenIconStyle(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'IconStyle')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenIconStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseIconStyle(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'IconStyle')
  end subroutine kmlCloseIconStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenLineStyle(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'LineStyle')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenLineStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseLineStyle(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'LineStyle')
  end subroutine kmlCloseLineStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenPolyStyle(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'PolyStyle')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenPolyStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlClosePolyStyle(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'PolyStyle')
  end subroutine kmlClosePolyStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenStyle(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Style')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseStyle(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Style')
  end subroutine kmlCloseStyle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenStyleMap(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'StyleMap')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenStyleMap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseStyleMap(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'StyleMap')
  end subroutine kmlCloseStyleMap
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<!-- abstract element; do not create -->
  !<!-- StyleSelector id="ID" -->                 <!-- Style,StyleMap -->
  !<!-- /StyleSelector -->

  subroutine kmlOpenStyleSelector(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'StyleSelector')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenStyleSelector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseStyleSelector(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'StyleSelector')
  end subroutine kmlCloseStyleSelector
#endif


  subroutine kmlCreatePolygonStyle(xf, color, colorhex, colorname, colormode, fill, outline, id)
    type(xmlf_t), intent(inout) :: xf
    logical, intent(in), optional :: fill
    logical, intent(in), optional :: outline
    type(col), intent(in), optional :: color
    character(len=*), intent(in), optional :: colorhex
    character(len=*), intent(in), optional :: colorname
    character(len=*), intent(in), optional :: colormode
    character(len=*), intent(in), optional :: id
#ifndef DUMMYLIB
    logical :: needStyle

    needStyle = (xmlf_OpenTag(xf)/='Style')

    if (needStyle) call kmlOpenStyle(xf,id)
    call kmlOpenPolyStyle(xf)

    if (count((/present(color), present(colorhex), present(colorname)/))>1) then
      print*, "Cannot specify more than one of color, colorhex, and colorname"
    endif
    if (present(color)) then
      call kmlAddColor(xf, color)
    elseif (present(colorhex)) then
      call kmlAddColor(xf, colorhex)
    elseif (present(colorname)) then
      call kmlAddColor(xf, kmlGetCustomColor(colorname))
    endif

    if (present(colormode)) call kmlAddColorMode(xf, colormode)

    if (present(fill)) call kmlAddFill(xf,fill)
    if (present(outline)) call kmlAddoutline(xf,outline)

    call kmlClosePolyStyle(xf)
    if (needStyle) call kmlCloseStyle(xf)
#endif
  end subroutine kmlCreatePolygonStyle
#ifndef DUMMYLIB
  subroutine kmlAddLegend(xf,myCI,vizvalues)
    ! add by GY 18/04/2008
    use m_wkml_color !required for derived type color
    use FoX_common   !required for function str

    type(xmlf_t), intent(inout) :: xf
    integer :: i
    double precision, intent(in) :: vizvalues(:)
    type(color_t),intent(in) :: myCI(:)
         do i=1,size(myCI)
          if (i == 1) then
          call kmlOpenPlacemark(xf)
           call kmlAddname(xf,"less than"//str(vizvalues(1),fmt="r3"))
            call kmlCreatePolygonStyle(xf,color=myCI(i))
            call kmlOpenPolygon(xf)
            call kmlClosePolygon(xf)
          call kmlClosePlacemark(xf)
          end if
          if (i /= 1 .and. i /= size(myCI)) then
          call kmlOpenPlacemark(xf)
           call kmlAddname(xf,str(vizvalues(i-1),fmt="r3")//"-"//str(vizvalues(i),fmt="r3"))
            call kmlCreatePolygonStyle(xf,color=myCI(i))
            call kmlOpenPolygon(xf)
            call kmlClosePolygon(xf)
          call kmlClosePlacemark(xf)
          end if
          if (i == size(myCI)) then 
          call kmlOpenPlacemark(xf)
           call kmlAddname(xf,"greater than"//str(vizvalues(i-1),fmt="r3"))
            call kmlCreatePolygonStyle(xf,color=myCI(i))
            call kmlOpenPolygon(xf)
            call kmlClosePolygon(xf)
          call kmlClosePlacemark(xf)
          end if
          end do
   end subroutine kmlAddLegend
#endif

end module m_wkml_styling
