program wkml_example
  ! Example point plotting in KML using FoX_wkml

  use FoX_wkml

  implicit none

  type(xmlf_t) :: myfile

  ! Data on the location of car accidents in Cambridge
  real :: latitude(10) = (/52.178179, 52.210026, 52.204842, &
          52.176018, 52.173411, 52.209044, 52.213750, 52.226834, 52.214031, 52.207677/)
  real :: longitude(10) = (/0.143026, 0.111787, 0.135543, & 
          0.143683, 0.112742, 0.141804, 0.123205, 0.111618, 0.110557, 0.127060/)

  call kmlBeginFile(myfile, "wkml_example.kml", -1)
  call kmlCreatePoints(myfile, longitude, latitude)
  call kmlFinishFile(myfile)


end program wkml_example
