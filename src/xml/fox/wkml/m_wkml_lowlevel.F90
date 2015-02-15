module m_wkml_lowlevel

#ifdef DUMMYLIB
  use FoX_wxml, only : xmlf_t
#else
  use m_common_error
  use fox_m_fsys_realtypes
  use FoX_wxml
  use FoX_common
  use FoX_utils, only: URI, parseURI, expressURI, destroyURI
#endif

  implicit none
  private
  
  ! Subroutines public to users of FoX_wkml (rexported in that module)
  public kmlOpenDocument
  public kmlOpenFolder
  public kmlCloseFolder
  public kmlCloseDocument

#ifndef DUMMYLIB
  ! Subroutines public inside wkml for kmlAdd
  public kmlAddname,kmlAddopen,kmlAddoutline,kmlAddextrude, &
    kmlAddfill,kmlAddvisibility,kmlAddtessellate,kmlAddrefreshVisibility,kmlAddflyToView,kmlAddaddress,kmlAddphoneNumber,  &
    kmlAdddescription,kmlAdddescription_ch,kmlAdddescription_sp,kmlAdddescription_dp,&
    kmlAddaltitudeMode,kmlAddrefreshMode,kmlAddviewrefreshMode,kmlAddviewRefreshTime,&
    kmlAddbegin,kmlAddend,kmlAddcolorMode,kmlAddmessage,kmlAddcookie,kmlAddlinkName, &
    kmlAddlinkDescription, &
    kmlAdddrawOrder,kmlAddsouth,kmlAddeast,kmlAddwest,kmlAddheading,kmlAddhref,kmlAddkey,kmlAddlongitude,kmlAddlatitude, &
    kmlAddaltitude,kmlAddminRefreshPeriod,kmlAddlistItemType,kmlAddmaxAltitude,kmlAddminAltitude,kmlAddminLodPixels, &
    kmlAddmaxLodPixels,kmlAddminFadeExtent,kmlAddmaxFadeExtent,kmlAddrange,kmlAddtilt,kmlAddroll,kmlAddrequest,kmlAddrotation, &
    kmlAddscale,kmlAddtargetHref,kmlAddviewBoundScale,kmlAddwhen,kmlAddwidth, & 
    kmlAddIcon_href,kmlAddIcon_refresh,kmlAddIcon_view, kmlAddStyleURL
  public :: kmlAddCoordinates

  ! Subroutines public inside wkml for kmlOpen and kmlClose
  public kmlOpencoordinates,kmlClosecoordinates,kmlOpenItemIcon,kmlCloseItemIcon,kmlAddstate, &
    kmlOpenAddressDetail,kmlCloseAddressDetail, kmlOpenChange,kmlCloseChange, &
    kmlOpenContainer,kmlCloseContainer,kmlOpenCreate,kmlCloseCreate,kmlOpenDelete,kmlCloseDelete, &
    kmlOpenFeature,kmlCloseFeature,kmlOpenGeometry,kmlCloseGeometry, &
    kmlOpenGeometryCollection,kmlCloseGeometryCollection,kmlOpenGroundOverlay,kmlCloseGroundOverlay,kmlOpenIcon,kmlCloseIcon, &
    kmlOpenLatLonAltBox,kmlCloseLatLonAltBox,kmlOpenLatLonBox,kmlCloseLatLonBox, &  
    kmlOpenLink,kmlCloseLink,kmlOpenLocation,kmlCloseLocation,kmlOpenLod,kmlCloseLod,kmlOpenLookAt,kmlCloseLookAt, &
    kmlOpenModel,kmlCloseModel,kmlOpenMultiGeometry,kmlCloseMultiGeometry,kmlOpenNetworkLink,kmlCloseNetworkLink, &
    kmlOpenNetworkLinkControl,kmlCloseNetworkLinkControl,kmlOpenObject,kmlCloseObject,kmlOpenObjArrayField,kmlCloseObjArrayField, &
    kmlOpenObjField,kmlCloseObjField,kmlOpenOrientation,kmlCloseOrientation,kmlOpenOverlay,kmlCloseOverlay, &
    kmlOpenPair,kmlClosePair,kmlOpenPlacemark,kmlClosePlacemark,kmlOpenPoint,kmlClosePoint, &
    kmlOpenRegion,kmlCloseRegion,kmlOpenResponse,kmlCloseResponse,kmlOpenScale,kmlCloseScale, &
    kmlOpenSchema,kmlCloseSchema,kmlOpenSchemaField,kmlCloseSchemaField,kmlOpenScreenOverlay,kmlCloseScreenOverlay, &  
    kmlOpenSimpleArrayField,kmlCloseSimpleArrayField,kmlOpenSimpleField,kmlCloseSimpleField,kmlOpenSnippet,kmlCloseSnippet, &
    kmlOpenStatus,kmlCloseStatus,kmlOpenTimePrimitive,kmlCloseTimePrimitive,kmlOpenTimeSpan,kmlCloseTimeSpan,kmlOpenTimeStamp, &
    kmlCloseTimeStamp,kmlOpenUpdate,kmlCloseUpdate,kmlOpenUrl,kmlCloseUrl

  ! kmlAdddescription support more data type
  interface kmlAdddescription
    module procedure kmlAdddescription_ch
    module procedure kmlAdddescription_sp
    module procedure kmlAdddescription_dp
  end interface

  interface kmlAddlongitude
    module procedure kmlAddlongitude_dp
    module procedure kmlAddlongitude_sp
  end interface

  interface kmlAddlatitude
    module procedure kmlAddlatitude_dp
    module procedure kmlAddlatitude_sp
  end interface

  interface kmlAddCoordinates
    module procedure kmlAddCoordinates_sp
    module procedure kmlAddCoordinates_dp
    module procedure kmlAddCoordinates_array_sp
    module procedure kmlAddCoordinates_array_dp
  end interface kmlAddCoordinates

  interface kmlAddScale
    module procedure kmlAddScale_int
    module procedure kmlAddScale_sp
    module procedure kmlAddScale_dp
  end interface kmlAddScale


  public :: kmlOpenInnerBoundaryIs, kmlCloseInnerBoundaryIs
  public :: kmlOpenOuterBoundaryIs, kmlCloseOuterBoundaryIs
  public :: kmlOpenLineString, kmlCloseLineString
  public :: kmlOpenLinearRing, kmlCloseLinearRing
  public :: kmlOpenPolygon, kmlClosePolygon
#endif
contains
#ifndef DUMMYLIB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddNamespace(xf, prefix, URI)
    type(xmlf_t), intent(inout) :: xf

    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: URI

    if (xmlf_OpenTag(xf) /= "") &
      call FoX_error("Cannot do kmlAddNamespace after document output")

    call xml_DeclareNamespace(xf, URI, prefix)
  end subroutine kmlAddNamespace


  subroutine kmlAddname(xf, name)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) ::  name
    call xml_NewElement(xf, 'name')
    call xml_AddCharacters(xf, name)
    call xml_EndElement(xf, 'name')
  end subroutine kmlAddname


  subroutine kmlAddopen(xf, folderopen)
    type(xmlf_t), intent(inout) :: xf
    logical, intent(in) :: folderopen
    character :: f
    f = merge('1', '0', folderopen)
    call xml_NewElement(xf, 'open')
    call xml_AddCharacters(xf, f)
    call xml_EndElement(xf, 'open')
  end subroutine kmlAddopen


  subroutine kmlAddoutline(xf, outline)
    type(xmlf_t), intent(inout) :: xf
    logical, intent(in) :: outline
    character :: o
    o = merge('1', '0', outline)
    call xml_NewElement(xf, 'outline')
    call xml_AddCharacters(xf, outline)
    call xml_EndElement(xf, 'outline')
  end subroutine kmlAddoutline


  subroutine kmlAddextrude(xf, extrude)
    type(xmlf_t), intent(inout) :: xf
    logical, intent(in) :: extrude
    character :: e
    e = merge('1', '0', extrude)
    call xml_NewElement(xf, 'extrude')
    call xml_AddCharacters(xf, e)
    call xml_EndElement(xf, 'extrude')
  end subroutine kmlAddextrude


  subroutine kmlAddfill(xf, fill)
    type(xmlf_t), intent(inout) :: xf
    logical, intent(in) :: fill
    character :: f
    f = merge('1', '0', fill)
    call xml_NewElement(xf, 'fill')
    call xml_AddCharacters(xf, f)
    call xml_EndElement(xf, 'fill')
  end subroutine kmlAddfill


  subroutine kmlAddvisibility(xf, visibility)
    type(xmlf_t), intent(inout) :: xf
    logical, intent(in) :: visibility
    character :: v
    v = merge('1', '0', visibility)
    call xml_NewElement(xf, 'visibility')
    call xml_AddCharacters(xf, v)
    call xml_EndElement(xf, 'visibility')
  end subroutine kmlAddvisibility


  subroutine kmlAddtessellate(xf, tessellate)
    type(xmlf_t), intent(inout) :: xf
    logical, intent(in) :: tessellate
    character :: t
    t = merge('1', '0', tessellate)
    call xml_NewElement(xf, 'tessellate')
    call xml_AddCharacters(xf, t)
    call xml_EndElement(xf, 'tessellate')
  end subroutine kmlAddtessellate


  subroutine kmlAddrefreshVisibility(xf, rvisibility)
    type(xmlf_t), intent(inout) :: xf
    character(len=1), intent(in) :: rvisibility
    call xml_NewElement(xf,'refreshVisibility')
    call xml_AddCharacters(xf,rvisibility)
    call xml_EndElement(xf,'refreshVisibility')
  end subroutine kmlAddrefreshVisibility
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddflyToView(xf, flytoview)
    type(xmlf_t), intent(inout) :: xf
    character(len=1), intent(in) :: flytoview
    call xml_NewElement(xf,'flyToView')
    call xml_AddCharacters(xf,flytoview)
    call xml_EndElement(xf,'flyToView')
  end subroutine kmlAddflyToView
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddaddress(xf, address)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: address
    call xml_NewElement(xf,'address')
    call xml_AddCharacters(xf,address)
    call xml_EndElement(xf,'address')
  end subroutine kmlAddaddress
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddphoneNumber(xf, phone)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: phone
    call xml_NewElement(xf,'phoneNumber')
    call xml_AddCharacters(xf,phone)
    call xml_EndElement(xf,'phoneNumber')
  end subroutine kmlAddphoneNumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !input as string
  subroutine kmlAdddescription_ch(xf, description)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: description
    call xml_NewElement(xf,'description')
    call xml_AddCharacters(xf,description)
    call xml_EndElement(xf,'description')
  end subroutine kmlAdddescription_ch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! input as real and convert to string using wxml_common function str
  subroutine kmlAdddescription_sp(xf, description2)
    use Fox_wxml
    use FoX_common
    type(xmlf_t), intent(inout) :: xf
    real, intent(in) :: description2
    call xml_NewElement(xf,'description')
    !           call xml_AddCharacters(xf,str(description2))
    call xml_AddCharacters(xf,description2,'r')
    call xml_EndElement(xf,'description')
  end subroutine kmlAdddescription_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! input as real and convert to string using wxml_common function str
  subroutine kmlAdddescription_dp(xf, description2)
    use Fox_wxml
    use FoX_common
    type(xmlf_t), intent(inout) :: xf
    double precision, intent(in) :: description2
    call xml_NewElement(xf,'description')
    !           call xml_AddCharacters(xf,str(description2))
    call xml_AddCharacters(xf,description2,'r')
    call xml_EndElement(xf,'description')
  end subroutine kmlAdddescription_dp


  subroutine kmlAddaltitudeMode(xf, altitudemode)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: altitudemode
    if (altitudemode/='clampToGround' &
      .and.altitudemode/='relativeToGround' &
      .and.altitudeMode/='absolute') then
      call FoX_error('unknown altitudeMode value')
    endif
    call xml_NewElement(xf, 'altitudeMode')
    call xml_AddCharacters(xf, altitudemode)
    call xml_EndElement(xf,'altitudeMode')
  end subroutine kmlAddaltitudeMode


  subroutine kmlAddrefreshMode(xf, refreshmode)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: refreshmode
    call xml_NewElement(xf,'refreshMode')
    call xml_AddCharacters(xf,refreshmode)
    call xml_EndElement(xf,'refreshMode')
  end subroutine kmlAddrefreshMode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddrefreshInterval(xf, refreshinterval)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: refreshinterval
    call xml_NewElement(xf,'refreshInterval')
    call xml_AddCharacters(xf,refreshinterval)
    call xml_EndElement(xf,'refreshInterval')
  end subroutine kmlAddrefreshInterval
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddviewrefreshMode(xf, viewrmode)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: viewrmode
    call xml_NewElement(xf,'viewRefreshMode')
    call xml_AddCharacters(xf,viewrmode)
    call xml_EndElement(xf,'viewRefreshMode')
  end subroutine kmlAddviewrefreshMode
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddviewRefreshTime(xf, viewrtime)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: viewrtime
    call xml_NewElement(xf,'viewRefreshTime')
    call xml_AddCharacters(xf,viewrtime)
    call xml_EndElement(xf,'viewRefreshTime')
  end subroutine kmlAddviewRefreshTime


  !<begin>1876-08-01</begin>
  subroutine kmlAddbegin(xf, begin)  ! 1984-08-11   date format
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: begin
    call xml_NewElement(xf,'begin')
    call xml_AddCharacters(xf,begin)
    call xml_EndElement(xf,'begin')
  end subroutine kmlAddbegin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<end>1876-08-01</end>
  subroutine kmlAddend(xf, enddate)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: enddate
    call xml_NewElement(xf,'end')
    call xml_AddCharacters(xf,enddate)
    call xml_EndElement(xf,'end')
  end subroutine kmlAddend


  !<colorMode>normal<colorMode>       <!-- kml:colorModeEnum:normal or random -->
  subroutine kmlAddcolorMode(xf, colormode)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: colormode
    if (colorMode/="normal".and.colormode/="random") then
      call FoX_error("Invalid value for colormode in kmlAddColorMode")
    endif
    call xml_NewElement(xf, 'colorMode')
    call xml_AddCharacters(xf, colormode)
    call xml_EndElement(xf,'colorMode')
  end subroutine kmlAddcolorMode


  ! <message>This is a pop-up message. You will only see this once</message>
  subroutine kmlAddmessage(xf, message)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: message
    call xml_NewElement(xf,'message')
    call xml_AddCharacters(xf,message)
    call xml_EndElement(xf,'message')
  end subroutine kmlAddmessage
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<cookie>cookie=sometext</cookie>
  subroutine kmlAddcookie(xf, cookie)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: cookie
    call xml_NewElement(xf,'cookie')
    call xml_AddCharacters(xf,cookie)
    call xml_EndElement(xf,'cookie')
  end subroutine kmlAddcookie
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<linkName>New KML features</linkName>
  subroutine kmlAddlinkName(xf, linkName)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: linkName
    call xml_NewElement(xf,'linkName')
    call xml_AddCharacters(xf,linkName)
    call xml_EndElement(xf,'linkName')
  end subroutine kmlAddlinkName
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<linkDescription><![CDATA[KML now has new features available!]]></linkDescription>
  subroutine kmlAddlinkDescription(xf, linkDescription)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: linkDescription
    call xml_NewElement(xf,'linkDescription')
    call xml_AddCharacters(xf,linkDescription)
    call xml_EndElement(xf,'linkDescription')
  end subroutine kmlAddlinkDescription
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! <drawOrder>0</drawOrder>                  <!-- int -->
  subroutine kmlAdddrawOrder(xf, drawOrder)
    type(xmlf_t), intent(inout) :: xf
    character(len=1), intent(in) :: drawOrder
    call xml_NewElement(xf,'drawOrder')
    call xml_AddCharacters(xf,drawOrder)
    call xml_EndElement(xf,'drawOrder')
  end subroutine kmlAdddrawOrder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<north>48.25475939255556</north>
  subroutine kmlAddnorth(xf, north)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: north
    call xml_NewElement(xf,'north')
    call xml_AddCharacters(xf,north)
    call xml_EndElement(xf,'north')
  end subroutine kmlAddnorth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<south>48.25207367852141</south>
  subroutine kmlAddsouth(xf, south)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: south
    call xml_NewElement(xf,'south')
    call xml_AddCharacters(xf,south)
    call xml_EndElement(xf,'south')
  end subroutine kmlAddsouth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<east>-90.86591508839973</east>
  subroutine kmlAddeast(xf, east)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: east
    call xml_NewElement(xf,'east')
    call xml_AddCharacters(xf,east)
    call xml_EndElement(xf,'east')
  end subroutine kmlAddeast
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<west>-90.8714285289695</west>
  subroutine kmlAddwest(xf, west)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: west
    call xml_NewElement(xf,'west')
    call xml_AddCharacters(xf,west)
    call xml_EndElement(xf,'west')
  end subroutine kmlAddwest


  !<heading>0</heading>  0-180
  subroutine kmlAddheading(xf, heading)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: heading
    ! Note that according to the schema the heading 
    ! is a kml:angle360Type which can be between 
    ! -360 and 360 degrees. This is not what the 
    ! google docs say (0 -> 360), but we follow the
    ! schema here. Also, as this derived from a schema 1.0 
    ! data type, we need to include the fmt argument to 
    ! avoid exponentail notation.
    !
    ! Schema: http://schemas.opengis.net/kml/2.2.0/ogckml22.xsd
    ! Google docs: http://code.google.com/apis/kml/documentation/kmlreference.html#iconstyle
    !
    if (heading<-360.0_dp.or.heading>360.0_dp) then
      call FoX_error("invalid value for heading")
    endif
    call xml_NewElement(xf, 'heading')
    call xml_AddCharacters(xf, heading, fmt='r10')
    call xml_EndElement(xf, 'heading')
  end subroutine kmlAddheading


  !<href>C:/Documents and Settings/All Users/Documents/My Pictures/Sample Pictures/Sunset.jpg</href>
  subroutine kmlAddhref(xf, url)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: url
    type(URI), pointer :: u
    call xml_NewElement(xf,'href')
    u => parseURI(url)
    if (.not.associated(u)) then
      call FoX_error("Invalid URI")
    endif
    call xml_AddCharacters(xf,expressURI(u))
    call destroyURI(u)
    call xml_EndElement(xf,'href')
  end subroutine kmlAddhref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<key>normal</key>   normal or highlight
  subroutine kmlAddkey(xf, key)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: key
    call xml_NewElement(xf,'key')
    call xml_AddCharacters(xf,key)
    call xml_EndElement(xf,'key')
  end subroutine kmlAddkey
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<longitude>39.55375305703105</longitude>
  subroutine kmlAddlongitude_dp(xf, longitude)
    type(xmlf_t), intent(inout) :: xf
    double precision, intent(in) :: longitude
    !        character(len=*), intent(in) :: longitude
    call xml_NewElement(xf,'longitude')
    call xml_AddCharacters(xf,longitude)
    call xml_EndElement(xf,'longitude')
  end subroutine kmlAddlongitude_dp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<longitude>39.55375305703105</longitude>
  subroutine kmlAddlongitude_sp(xf, longitude)
    type(xmlf_t), intent(inout) :: xf
    real, intent(in) :: longitude
    !        character(len=*), intent(in) :: longitude
    call xml_NewElement(xf,'longitude')
    call xml_AddCharacters(xf,longitude)
    call xml_EndElement(xf,'longitude')
  end subroutine kmlAddlongitude_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<latitude>-118.9813220168456</latitude>
  subroutine kmlAddlatitude_dp(xf, latitude)
    type(xmlf_t), intent(inout) :: xf
    double precision, intent(in) :: latitude
    !        character(len=*), intent(in) :: latitude
    call xml_NewElement(xf,'latitude')
    call xml_AddCharacters(xf,latitude)
    call xml_EndElement(xf,'latitude')
  end subroutine kmlAddlatitude_dp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<latitude>-118.9813220168456</latitude>
  subroutine kmlAddlatitude_sp(xf, latitude)
    type(xmlf_t), intent(inout) :: xf
    real, intent(in) :: latitude
    !        character(len=*), intent(in) :: latitude
    call xml_NewElement(xf,'latitude')
    call xml_AddCharacters(xf,latitude)
    call xml_EndElement(xf,'latitude')
  end subroutine kmlAddlatitude_sp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<altitude>1223</altitude>
  subroutine kmlAddaltitude(xf, altitude)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: altitude
    call xml_NewElement(xf,'altitude')
    call xml_AddCharacters(xf,altitude)
    call xml_EndElement(xf,'altitude')
  end subroutine kmlAddaltitude
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! seconds default=0
  subroutine kmlAddminRefreshPeriod(xf, minRefreshPeriod)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: minRefreshPeriod
    call xml_NewElement(xf,'minRefreshPeriod')
    call xml_AddCharacters(xf,minRefreshPeriod)
    call xml_EndElement(xf,'minRefreshPeriod')
  end subroutine kmlAddminRefreshPeriod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<listItemType>checkHideChildren</listItemType>
  subroutine kmlAddlistItemType(xf, listItemType)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: listItemType
    call xml_NewElement(xf,'listItemType')
    call xml_AddCharacters(xf,listItemType)
    call xml_EndElement(xf,'listItemType')
  end subroutine kmlAddlistItemType
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<maxAltitude>0</maxAltitude>
  !Defaults to 0; specified in meters above sea level (and is affected by the <altitudeMode> specification)
  subroutine kmlAddmaxAltitude(xf, maxAltitude)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: maxAltitude
    call xml_NewElement(xf,'maxAltitude')
    call xml_AddCharacters(xf,maxAltitude)
    call xml_EndElement(xf,'maxAltitude')
  end subroutine kmlAddmaxAltitude
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<minAltitude>0</minAltitude>
  subroutine kmlAddminAltitude(xf, minAltitude)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: minAltitude
    call xml_NewElement(xf,'minAltitude')
    call xml_AddCharacters(xf,minAltitude)
    call xml_EndElement(xf,'minAltitude')
  end subroutine kmlAddminAltitude
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<minLodPixels>128</minLodPixels>
  !<minLodPixels> (default = 0) Measurement in screen pixels that represents the minimum limit 
  !of the visibility range for a given Region. Google Earth calculates the size of the Region
  ! when projected onto screen space. Then it computes the square root of the Region's area
  ! (if, for example, the Region is square and the viewpoint is directly above the Region,
  ! and the Region is not tilted, this measurement is equal to the width of the projected
  ! Region). If this measurement falls within the limits defined by <minLodPixels> and 
  !<maxLodPixels> (and if the <LatLonAltBox> is in view), the Region is active. If
  ! this limit is not reached, the associated geometry is considered to be too far 
  !from the user's viewpoint to be drawn.
  subroutine kmlAddminLodPixels(xf, minLodPixels)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: minLodPixels
    call xml_NewElement(xf,'minLodPixels')
    call xml_AddCharacters(xf,minLodPixels)
    call xml_EndElement(xf,'minLodPixels')
  end subroutine kmlAddminLodPixels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<maxLodPixels>1024</maxLodPixels>
  !<maxLodPixels> (default = -1) Measurement in screen pixels that represents the maximum 
  !limit of the visibility range for a given Region.
  subroutine kmlAddmaxLodPixels(xf, maxLodPixels)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: maxLodPixels
    call xml_NewElement(xf,'maxLodPixels')
    call xml_AddCharacters(xf,maxLodPixels)
    call xml_EndElement(xf,'maxLodPixels')
  end subroutine kmlAddmaxLodPixels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<minFadeExtent>128</minFadeExtent>
  !minFadeExtent> (default = 0) Distance over which the geometry fades, from fully opaque to fully transparent. This ramp value, expressed in screen pixels, is applied at the minimum end of the LOD (visibility) limits.
  subroutine kmlAddminFadeExtent(xf, minFadeExtent)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: minFadeExtent
    call xml_NewElement(xf,'minFadeExtent')
    call xml_AddCharacters(xf,minFadeExtent)
    call xml_EndElement(xf,'minFadeExtent')
  end subroutine kmlAddminFadeExtent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<maxFadeExtent>128</maxFadeExtent>
  !<maxFadeExtent> (default = 0) Distance over which the geometry fades, from fully transparent to fully opaque. This ramp value, expressed in screen pixels, is applied at the maximum end of the LOD (visibility) limits.
  subroutine kmlAddmaxFadeExtent(xf, maxFadeExtent)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: maxFadeExtent
    call xml_NewElement(xf,'maxFadeExtent')
    call xml_AddCharacters(xf,maxFadeExtent)
    call xml_EndElement(xf,'maxFadeExtent')
  end subroutine kmlAddmaxFadeExtent
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<range>500</range>
  subroutine kmlAddrange(xf, range1)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: range1
    call xml_NewElement(xf,'range')
    call xml_AddCharacters(xf,range1)
    call xml_EndElement(xf,'range')
  end subroutine kmlAddrange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<tilt>45</tilt>
  subroutine kmlAddtilt(xf, tilt)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: tilt
    call xml_NewElement(xf,'tilt')
    call xml_AddCharacters(xf,tilt)
    call xml_EndElement(xf,'tilt')
  end subroutine kmlAddtilt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<roll>0</roll>          <!-- kml:angle360 -->
  subroutine kmlAddroll(xf, roll)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: roll
    call xml_NewElement(xf,'roll')
    call xml_AddCharacters(xf,roll)
    call xml_EndElement(xf,'roll')
  end subroutine kmlAddroll
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<request>geocode</request>  !currently, it is geocode
  subroutine kmlAddrequest(xf, geocode)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: geocode
    call xml_NewElement(xf,'request')
    call xml_AddCharacters(xf,geocode)
    call xml_EndElement(xf,'request')
  end subroutine kmlAddrequest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<rotation>39.37878630116985</rotation>
  !Indicates the angle of rotation of the parent object. A value of 0 means no rotation.
  ! The value is an angle in degrees counterclockwise starting from north. Use ±180 to
  ! indicate the rotation of the parent object from 0. The center of the <rotation>, 
  !if not (.5,.5), is specified in <rotationXY>.
  subroutine kmlAddrotation(xf, rotation)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: rotation
    call xml_NewElement(xf,'rotation')
    call xml_AddCharacters(xf,rotation)
    call xml_EndElement(xf,'rotation')
  end subroutine kmlAddrotation

  subroutine kmlAddScale_int(xf, scale)
    type(xmlf_t), intent(inout) :: xf
    integer, intent(in) :: scale

    call xml_NewElement(xf, 'scale')
    call xml_AddCharacters(xf, scale)
    call xml_EndElement(xf, 'scale')
  end subroutine kmlAddScale_int

  subroutine kmlAddScale_sp(xf, scale)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: scale

    call xml_NewElement(xf, 'scale')
    call xml_AddCharacters(xf, str(scale, "r"))
    call xml_EndElement(xf, 'scale')
  end subroutine kmlAddScale_sp

  subroutine kmlAddScale_dp(xf, scale)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: scale

    call xml_NewElement(xf, 'scale')
    call xml_AddCharacters(xf, str(scale, "r"))
    call xml_EndElement(xf, 'scale')
  end subroutine kmlAddScale_dp


  !<targetHref>http://www/~sam/January14Data/Point.kml</targetHref>
  subroutine kmlAddtargetHref(xf, targetHref)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: targetHref
    call xml_NewElement(xf,'targetHref')
    call xml_AddCharacters(xf,targetHref)
    call xml_EndElement(xf,'targetHref')
  end subroutine kmlAddtargetHref
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<viewBoundScale>
  !Scales the BBOX parameters before sending them to the server. Default = 1. A value < 1 specifies to use less than the full view (screen). A value >1 specifies to fetch an area that extends beyond the edges of the current view.
  subroutine kmlAddviewBoundScale(xf, viewBoundScale)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: viewBoundScale
    call xml_NewElement(xf,'viewBoundScale')
    call xml_AddCharacters(xf,viewBoundScale)
    call xml_EndElement(xf,'viewBoundScale')
  end subroutine kmlAddviewBoundScale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<when>1997</when>
  !<when>1997-07</when>
  !<when>1997-07-16</when>
  !<when>1997-07-16T07:30:15Z</when>
  !<when>1997-07-16T10:30:15+03:00</when>
  subroutine kmlAddwhen(xf, when)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: when
    call xml_NewElement(xf,'when')
    call xml_AddCharacters(xf,when)
    call xml_EndElement(xf,'when')
  end subroutine kmlAddwhen


  subroutine kmlAddwidth(xf, width)
    type(xmlf_t), intent(inout) :: xf
    integer, intent(in) :: width
    call xml_NewElement(xf, 'width')
    call xml_AddCharacters(xf, width)
    call xml_EndElement(xf, 'width')
  end subroutine kmlAddwidth


  !<coordinates>-90.86948943473118,48.25450093195546</coordinates>
  !A single tuple consisting of floating point values for longitude, 
  !latitude, and altitude (in that order). Longitude and latitude 
  !values are in degrees, where longitude ≥ -180 and <= 180 latitude 
  ! ≥ -90 and ≤ 90 altitude values (optional) are in meters above 
  !sea level Do not include spaces between the three values that describe a coordinate.
  subroutine kmlAddcoordinates_sp(xf, x, y, z)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: x,y
    real(sp), intent(in), optional :: z
    call xml_NewElement(xf, 'coordinates')
    call xml_AddCharacters(xf, x, 'r6')
    call xml_AddCharacters(xf, ',')
    call xml_AddCharacters(xf, y, 'r6')
    if (present(z)) then
      call xml_AddCharacters(xf, ',')
      call xml_AddCharacters(xf, z, 'r6')
    endif
    call xml_EndElement(xf, 'coordinates')
  end subroutine kmlAddcoordinates_sp

  subroutine kmlAddcoordinates_dp(xf, x, y, z)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: x,y
    real(dp), intent(in), optional :: z
    call xml_NewElement(xf, 'coordinates')
    call xml_AddCharacters(xf, x, 'r6')
    call xml_AddCharacters(xf, ',')
    call xml_AddCharacters(xf, y, 'r6')
    if (present(z)) then
      call xml_AddCharacters(xf, ',')
      call xml_AddCharacters(xf, z, 'r6')
    endif
    call xml_EndElement(xf,'coordinates')
  end subroutine kmlAddcoordinates_dp


  subroutine kmlAddIcon_href(xf, url,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in) :: url
    call xml_NewElement(xf,'Icon')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    call kmlAddhref(xf,url)
    call xml_EndElement(xf,'Icon')
  end subroutine kmlAddIcon_href
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddIcon_refresh(xf, url, rmode, rinterval,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in) :: url, rmode, rinterval
    call xml_NewElement(xf,'Icon')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    call kmlAddhref(xf,url)
    call kmlAddrefreshMode(xf,rmode)
    call kmlAddrefreshInterval(xf,rinterval)
    call xml_EndElement(xf,'Icon')
  end subroutine kmlAddIcon_refresh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlAddIcon_view(xf, url, rmode, rinterval,vrmode, vrtime, vbscale,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in) :: url, rmode, rinterval,vrmode, vrtime,vbscale
    call xml_NewElement(xf,'Icon')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    call kmlAddhref(xf,url)
    call kmlAddrefreshMode(xf,rmode)
    call kmlAddrefreshInterval(xf,rinterval)
    call kmlAddviewRefreshMode(xf,vrmode)
    call kmlAddviewRefreshTime(xf,vrtime)
    call kmlAddviewBoundScale(xf, vbscale)
    call xml_EndElement(xf,'Icon')
  end subroutine kmlAddIcon_view
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpencoordinates(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_NewElement(xf,'coordinates')
  end subroutine kmlOpencoordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlClosecoordinates(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'coordinates')
  end subroutine kmlClosecoordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenItemIcon(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_NewElement(xf,'ItemIcon')
  end subroutine kmlOpenItemIcon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseItemIcon(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'ItemIcon')
  end subroutine kmlCloseItemIcon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<state> - specifies the current state of the NetworkLink or Folder. Possible values are open, closed, error, fetching0, fetching1, and fetching2. These values can be combined by inserting a space between two values (no comma).
  !<!-- kml:itemIconModeEnum:open, closed, error, fetching0, fetching1, or fetching2 -->
  subroutine kmlAddstate(xf, state)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: state
    call xml_NewElement(xf,'state')
    call xml_AddCharacters(xf,state)
    call xml_EndElement(xf,'state')
  end subroutine kmlAddstate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenAddressDetail(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'AddressDetail')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenAddressDetail
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseAddressDetail(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'AddressDetail')
  end subroutine kmlCloseAddressDetail
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenChange(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Change')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenChange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseChange(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Change')
  end subroutine kmlCloseChange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenContainer(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Container')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenContainer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseContainer(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Container')
  end subroutine kmlCloseContainer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenCreate(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Create')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenCreate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseCreate(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Create')
  end subroutine kmlCloseCreate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenDelete(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Delete')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenDelete
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseDelete(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Delete')
  end subroutine kmlCloseDelete
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenDocument(xf,name,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in):: name
    character(len=*), intent(in), optional :: id
#ifndef DUMMYLIB
    call xml_NewElement(xf,'Document')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    call kmladdname(xf,name)
#endif
  end subroutine kmlOpenDocument
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseDocument(xf)
    type(xmlf_t), intent(inout) :: xf
#ifndef DUMMYLIB
    call xml_EndElement(xf,'Document')
#endif
  end subroutine kmlCloseDocument
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenFolder(xf,id,name)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional:: name
#ifndef DUMMYLIB
    call xml_NewElement(xf,'Folder')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    if (present(name)) then
      call kmladdname(xf,name)
    end if
#endif
  end subroutine kmlOpenFolder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseFolder(xf)
    type(xmlf_t), intent(inout) :: xf
#ifndef DUMMYLIB
    call xml_EndElement(xf,'Folder')
#endif
  end subroutine kmlCloseFolder
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef DUMMYLIB
  subroutine kmlOpenFeature(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Feature')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenFeature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseFeature(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Feature')
  end subroutine kmlCloseFeature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenGeometry(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Geometry')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenGeometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseGeometry(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Geometry')
  end subroutine kmlCloseGeometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenGeometryCollection(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'GeometryCollection')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenGeometryCollection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseGeometryCollection(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'GeometryCollection')
  end subroutine kmlCloseGeometryCollection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenGroundOverlay(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'GroundOverlay')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenGroundOverlay
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseGroundOverlay(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'GroundOverlay')
  end subroutine kmlCloseGroundOverlay
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenIcon(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Icon')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenIcon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseIcon(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Icon')
  end subroutine kmlCloseIcon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenLatLonAltBox(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'LatLonAltBox')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenLatLonAltBox
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseLatLonAltBox(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'LatLonAltBox')
  end subroutine kmlCloseLatLonAltBox
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenLatLonBox(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'LatLonBox')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenLatLonBox
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseLatLonBox(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'LatLonBox')
  end subroutine kmlCloseLatLonBox

  !<styleUrl>#myIconStyleID</styleUrl>
  !<styleUrl>http://someserver.com/somestylefile.xml#restaurant</styleUrl>
  subroutine kmlAddstyleUrl(xf, styleUrl)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: styleUrl
    call xml_NewElement(xf, 'styleUrl')
    call xml_AddCharacters(xf, styleUrl)
    call xml_EndElement(xf, 'styleUrl')
  end subroutine kmlAddstyleUrl

  subroutine kmlAddCoordinates_array_sp(xf, x, y, z, repeat)
    type(xmlf_t), intent(inout) :: xf
    real(sp), intent(in) :: x(:)
    real(sp), intent(in) :: y(:)
    real(sp), intent(in), optional :: z(:)
    logical, intent(in), optional :: repeat
    
    integer :: nodes, i
    logical :: repeat_
    if (present(repeat)) then
      repeat_ = repeat
    else
      repeat_ = .false.
    endif

    nodes = size(x)
    if (nodes/=size(y)) then
      call FoX_error("Inconsistent array lengths in kmlAddCoordinates")
    endif
    if (present(z)) then
      if (nodes/=size(z)) then
        call FoX_error("Inconsistent array lengths in kmlAddCoordinates")
      endif
    endif

    call xml_NewElement(xf, 'coordinates')
    call xml_AddNewline(xf)
    do i = 1, nodes
      call xml_AddCharacters(xf, str(x(i), 'r6')//','//str(y(i), 'r6'))
      if (present(z)) then
        call xml_AddCharacters(xf, ','//str(z(i), 'r6'))
      end if
      call xml_AddNewline(xf)
    end do
    if (repeat_) then
      call xml_AddCharacters(xf, str(x(1), 'r6')//','//str(y(1), 'r6'))
      if (present(z)) then
        call xml_AddCharacters(xf, ','//str(z(1), 'r6'))
      end if
      call xml_AddNewline(xf)
    endif
    call xml_EndElement(xf, 'coordinates')
  end subroutine kmlAddCoordinates_array_sp

  subroutine kmlAddCoordinates_array_dp(xf, x, y, z, repeat)
    type(xmlf_t), intent(inout) :: xf
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in), optional :: z(:)
    logical, intent(in), optional :: repeat
    
    integer :: nodes, i
    logical :: repeat_
    if (present(repeat)) then
      repeat_ = repeat
    else
      repeat_ = .false.
    endif

    nodes = size(x)
    if (nodes/=size(y)) then
      call FoX_error("Inconsistent array lengths in kmlAddCoordinates")
    endif
    if (present(z)) then
      if (nodes/=size(z)) then
        call FoX_error("Inconsistent array lengths in kmlAddCoordinates")
      endif
    endif

    call xml_NewElement(xf, 'coordinates')
    call xml_AddNewline(xf)
    do i = 1, nodes
      call xml_AddCharacters(xf, str(x(i), 'r6')//','//str(y(i), 'r6'))
      if (present(z)) then
        call xml_AddCharacters(xf, ','//str(z(i), 'r6'))
      end if
      call xml_AddNewline(xf)
    end do
    if (repeat_) then
      call xml_AddCharacters(xf, str(x(1), 'r6')//','//str(y(1), 'r6'))
      if (present(z)) then
        call xml_AddCharacters(xf, ','//str(z(1), 'r6'))
      end if
      call xml_AddNewline(xf)
    endif
    call xml_EndElement(xf, 'coordinates')
  end subroutine kmlAddCoordinates_array_dp


  subroutine kmlOpenLink(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Link')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenLink
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseLink(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Link')
  end subroutine kmlCloseLink
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenLocation(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Location')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenLocation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseLocation(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Location')
  end subroutine kmlCloseLocation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenLod(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Lod')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenLod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseLod(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Lod')
  end subroutine kmlCloseLod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenLookAt(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'LookAt')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenLookAt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseLookAt(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'LookAt')
  end subroutine kmlCloseLookAt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenModel(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Model')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenModel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseModel(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Model')
  end subroutine kmlCloseModel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenMultiGeometry(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'MultiGeometry')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenMultiGeometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseMultiGeometry(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'MultiGeometry')
  end subroutine kmlCloseMultiGeometry
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenNetworkLink(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'NetworkLink')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenNetworkLink
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseNetworkLink(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'NetworkLink')
  end subroutine kmlCloseNetworkLink
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenNetworkLinkControl(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'NetworkLinkControl')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenNetworkLinkControl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseNetworkLinkControl(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'NetworkLinkControl')
  end subroutine kmlCloseNetworkLinkControl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<!-- abstract element; do not create -->
  subroutine kmlOpenObject(xf,id,targetid)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id, targetid
    call xml_NewElement(xf,'Object')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    if (present(targetid)) call xml_AddAttribute(xf,'targetId', targetid)
  end subroutine kmlOpenObject
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseObject(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Object')
  end subroutine kmlCloseObject
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenObjArrayField(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'ObjArrayField')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenObjArrayField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseObjArrayField(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'ObjArrayField')
  end subroutine kmlCloseObjArrayField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenObjField(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'ObjField')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenObjField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseObjField(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'ObjField')
  end subroutine kmlCloseObjField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenOrientation(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Orientation')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenOrientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseOrientation(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Orientation')
  end subroutine kmlCloseOrientation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenOverlay(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Overlay')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenOverlay
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseOverlay(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Overlay')
  end subroutine kmlCloseOverlay
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenPair(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Pair')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenPair
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlClosePair(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Pair')
  end subroutine kmlClosePair


  subroutine kmlOpenPlacemark(xf, id, name, description, styleurl)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: name
    character(len=*), intent(in), optional :: description
    character(len=*), intent(in), optional :: styleurl
    call xml_NewElement(xf,'Placemark')
    if (present(id)) call xml_AddAttribute(xf, 'id', id)
    if (present(name)) call kmlAddName(xf, name)
    if (present(description)) call kmlAddDescription(xf, description)
    if (present(styleurl)) call kmlAddStyleURL(xf, styleurl)
  end subroutine kmlOpenPlacemark


  subroutine kmlClosePlacemark(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Placemark')
  end subroutine kmlClosePlacemark


  subroutine kmlOpenPoint(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Point')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenPoint
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlClosePoint(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Point')
  end subroutine kmlClosePoint


  subroutine kmlOpenRegion(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Region')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenRegion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseRegion(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Region')
  end subroutine kmlCloseRegion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !he accuracy attribute indicates how accurately the given address was able to be geocoded
  subroutine kmlOpenResponse(xf,accuracy)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: accuracy
    call xml_NewElement(xf,'Response')
    if (present(accuracy)) call xml_AddAttribute(xf,'accuracy', accuracy)
  end subroutine kmlOpenResponse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseResponse(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Response')
  end subroutine kmlCloseResponse
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenScale(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Scale')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenScale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseScale(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Scale')
  end subroutine kmlCloseScale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Specifies a custom KML schema that is typically used to extend and add metadata to KML objects. The "name" attribute is required. Currently, the only value for parent is "Placemark."

  subroutine kmlOpenSchema(xf,name1,parent)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name1
    character(len=*), intent(in), optional :: parent
    call xml_NewElement(xf,'Schema')
    call xml_AddAttribute(xf,'name', name1)
    if (present(parent)) call xml_AddAttribute(xf,'parent', parent)
  end subroutine kmlOpenSchema
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseSchema(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Schema')
  end subroutine kmlCloseSchema
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<!-- abstract element; do not create -->
  subroutine kmlOpenSchemaField(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'SchemaField')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenSchemaField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseSchemaField(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'SchemaField')
  end subroutine kmlCloseSchemaField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenScreenOverlay(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'ScreenOverlay')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenScreenOverlay
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseScreenOverlay(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'ScreenOverlay')
  end subroutine kmlCloseScreenOverlay
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenSimpleArrayField(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'SimpleArrayField')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenSimpleArrayField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseSimpleArrayField(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'SimpleArrayField')
  end subroutine kmlCloseSimpleArrayField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !<SimpleField name="AREA" type="double"></SimpleField>
  subroutine kmlOpenSimpleField(xf,name1,type1)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in) :: name1, type1
    call xml_NewElement(xf,'SimpleField')
    call xml_AddAttribute(xf,'name', name1)
    call xml_AddAttribute(xf,'type', type1)
  end subroutine kmlOpenSimpleField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseSimpleField(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'SimpleField')
  end subroutine kmlCloseSimpleField
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !A short description of the feature. In Google Earth, this description is displayed in the
  ! Places panel under the name of the feature. If a Snippet is not supplied, the first two 
  !lines of the <description> are used. In Google Earth, if a Placemark contains both a  
  !description and a Snippet, the <Snippet> appears beneath the Placemark in the Places 
  !panel, and the <description> appears in the Placemark's description balloon. This tag
  ! does not support HTML markup. <Snippet> has a maxLines attribute, an integer that
  ! specifies the maximum number of lines to display. Default for maxLines is 2.

  subroutine kmlOpenSnippet(xf,maxlines)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: maxlines
    call xml_NewElement(xf,'Snippet')
    if (present(maxlines)) call xml_AddAttribute(xf,'maxLines', maxlines)
  end subroutine kmlOpenSnippet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseSnippet(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Snippet')
  end subroutine kmlCloseSnippet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenStatus(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Status')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenStatus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseStatus(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Status')
  end subroutine kmlCloseStatus
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenTimePrimitive(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'TimePrimitive')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenTimePrimitive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseTimePrimitive(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'TimePrimitive')
  end subroutine kmlCloseTimePrimitive
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenTimeSpan(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'TimeSpan')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenTimeSpan
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseTimeSpan(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'TimeSpan')
  end subroutine kmlCloseTimeSpan
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlOpenTimeStamp(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'TimeStamp')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenTimeStamp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseTimeStamp(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'TimeStamp')
  end subroutine kmlCloseTimeStamp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  subroutine kmlOpenUpdate(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Update')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenUpdate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseUpdate(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Update')
  end subroutine kmlCloseUpdate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !This element is deprecated in KML Release 2.1 and has been replaced by <Link>, which provides the additional functionality of Regions.
  subroutine kmlOpenUrl(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf,'Url')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenUrl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine kmlCloseUrl(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Url')
  end subroutine kmlCloseUrl


  subroutine kmlOpenOuterBoundaryIs(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf, 'outerBoundaryIs')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenOuterBoundaryIs

  subroutine kmlCloseouterBoundaryIs(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf, 'outerBoundaryIs')
  end subroutine kmlCloseouterBoundaryIs

  subroutine kmlOpenInnerBoundaryIs(xf,id)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    call xml_NewElement(xf, 'innerBoundaryIs')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
  end subroutine kmlOpenInnerBoundaryIs

  subroutine kmlCloseInnerBoundaryIs(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf, 'innerBoundaryIs')
  end subroutine kmlCloseInnerBoundaryIs

  subroutine kmlOpenLinearRing(xf, id, altitudeMode, tessellate, extrude)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: altitudeMode
    logical, intent(in), optional :: tessellate
    logical, intent(in), optional :: extrude

    if (present(extrude).and.present(altitudeMode)) then
      if (extrude.and.altitudeMode=='clampToGround') then
        print*, "Inconsistent settings for extrude and altitudeMode"
      endif
    endif

    call xml_NewElement(xf, 'LinearRing')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    if (present(tessellate)) call kmlAddtessellate(xf, tessellate)
    if (present(altitudeMode)) call kmlAddAltitudeMode(xf, altitudeMode)
    if (present(extrude)) call kmlAddExtrude(xf, extrude)

  end subroutine kmlOpenLinearRing

  subroutine kmlCloseLinearRing(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'LinearRing')
  end subroutine kmlCloseLinearRing


  subroutine kmlOpenLineString(xf, id, altitudeMode, tessellate, extrude)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: altitudeMode
    logical, intent(in), optional :: tessellate
    logical, intent(in), optional :: extrude

    if (present(extrude).and.present(altitudeMode)) then
      if (extrude.and.altitudeMode=='clampToGround') then
        print*, "Inconsistent settings for extrude and altitudeMode"
      endif
    endif

    call xml_NewElement(xf, 'LineString')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    if (present(tessellate)) call kmlAddtessellate(xf, tessellate)
    if (present(altitudeMode)) call kmlAddAltitudeMode(xf, altitudeMode)
    if (present(extrude)) call kmlAddExtrude(xf, extrude)

  end subroutine kmlOpenLineString

  subroutine kmlCloseLineString(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'LineString')
  end subroutine kmlCloseLineString

  subroutine kmlOpenPolygon(xf, id, altitudeMode, tessellate, extrude)
    type(xmlf_t), intent(inout) :: xf
    character(len=*), intent(in), optional :: id
    character(len=*), intent(in), optional :: altitudeMode
    logical, intent(in), optional :: tessellate
    logical, intent(in), optional :: extrude

    if (present(extrude).and.present(altitudeMode)) then
      if (extrude.and.altitudeMode=='clampToGround') then
        print*, "Inconsistent settings for extrude and altitudeMode"
      endif
    endif

    call xml_NewElement(xf, 'Polygon')
    if (present(id)) call xml_AddAttribute(xf,'id', id)
    if (present(tessellate)) call kmlAddtessellate(xf, tessellate)
    if (present(altitudeMode)) call kmlAddAltitudeMode(xf, altitudeMode)
    if (present(extrude)) call kmlAddExtrude(xf, extrude)

  end subroutine kmlOpenPolygon

  subroutine kmlClosePolygon(xf)
    type(xmlf_t), intent(inout) :: xf
    call xml_EndElement(xf,'Polygon')
  end subroutine kmlClosePolygon  

#endif

end module m_wkml_lowlevel
