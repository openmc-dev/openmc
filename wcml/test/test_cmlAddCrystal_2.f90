program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  double precision, dimension(3) :: lens, angs

  lens = (/1.0d0, 1.0d0, 1.0d0/)
  angs = (/90.0d0, 90.0d0, 90.0d0/)

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  
  call cmlAddCrystal(xf, lens(1), lens(2), lens(3), angs(1), angs(2), angs(3), &
       lenfmt="r3", angfmt="r3")

  call cmlFinishFile(xf)

end program test
