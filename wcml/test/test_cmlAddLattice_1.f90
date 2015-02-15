program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  real, dimension(3,3) :: ucell

  ucell = reshape((/  1.0, 0.5, 0.5,   &
                      0.0, 1.0, 0.0,   &
                      0.0, 0.0, 1.0/), &
                   (/3,3/))

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  
  call cmlAddLattice(xf, cell=ucell, fmt="s3")

  call cmlFinishFile(xf)

end program test
