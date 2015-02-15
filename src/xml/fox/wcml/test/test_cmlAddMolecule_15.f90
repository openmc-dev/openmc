program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  double precision, dimension(3,4) :: coords
  character(len=2), dimension(4) :: elems

  coords = reshape((/  0.0d0,    0.0d0, 0.0d0,   &
                       0.0d0,  1.203d0, 0.0d0,   &
                    -0.934d0, -0.579d0, 0.0d0,   &
                     0.934d0, -0.579d0, 0.0d0/), &
                   (/3,4/))

  elems = (/"C ", "O ", "H ", "H "/)

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  
  call cmlAddMolecule(xf, coords=coords, elements=elems, title="Formaldehyde")

  call cmlFinishFile(xf)

end program test
