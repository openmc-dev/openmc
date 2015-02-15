program test

use FoX_wkml
use m_contours_test_data_sp, only : long, lat, cdata, setup_values

implicit none

type(xmlf_t) :: xf

call setup_values()

call kmlBeginFile(xf,"test.xml", -1)

call kmlCreateContours(xf,long,lat,values=cdata,name='test contours', & 
                       contour_values=(/0.0,0.001,0.002,0.003,0.004,   &
                                        0.005,0.006,0.007,0.01,0.02/), & 
                       regions=.true.)

call kmlFinishFile(xf)

end program test
