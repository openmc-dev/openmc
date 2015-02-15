program test

use FoX_wkml
use m_contours_test_data_sp, only : long, lat, cdata, setup_values

implicit none

type(xmlf_t) :: xf
real :: latpoints(61,69)
real :: longpoints(61,69)
integer :: i

call setup_values()

forall(i = 1:69) latpoints(:,i) = lat(i)
forall(i = 1:61) longpoints(i,:) = long(i)

call kmlBeginFile(xf,"test.xml", -1)

call kmlCreateContours(xf,longpoints,latpoints,values=cdata,name='test contours', & 
                       contour_values=(/0.0,0.001,0.002,0.003,0.004,   &
                                        0.005,0.006,0.007,0.01,0.02/), & 
                       regions=.true.)

call kmlFinishFile(xf)

end program test
