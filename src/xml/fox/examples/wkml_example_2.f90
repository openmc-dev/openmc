program path
   ! Draw a line in the sky from the UK to Florida
   use FoX_wkml

   type(xmlf_t) :: myfile
   character(len=25) :: filename
   real :: latitude(10) = (/52.0000, 50.0000, 48.0000, 45.0000, 40.0000, 35.0000, 30.0000, 28.0000, 25.0000, 23.0000/)
   real :: longitude(10) = (/0.0000, -5.0000, -10.0000, -20.0000, -30.0000, -40.0000, -50.0000, -60.0000, -70.0000, -80.0000/)
   real :: zed(10) = (/2000.00, 10000.00, 30000.00, 60000.00, 80000.00, 80000.00, 80000.00, 60000.00, 40000.00, 30000.00/)

   filename='wkml_example_2.kml'


   call kmlBeginFile(myfile, filename, -1)
    call kmlOpenDocument(myfile, "kmlCreateLine-path")
      ! create the line style and th ereference id is "path"
      call kmlCreateLineStyle(myfile, width=5,colorname="yellow", id="path")
      ! create line and reference to the style "path" using argument styleURL="#path"
      call kmlCreateLine(myfile, longitude, latitude, altitude=zed, &
      altitudeMode="absolute",styleURL="#path")
    call kmlCloseDocument(myfile)
   call kmlFinishFile(myfile)


end program path
