program test

  use FoX_wcml
  implicit none

  character(len=*), parameter :: filename = 'test.xml'
  type(xmlf_t) :: xf

  real, dimension(3,4) :: coords
  character(len=2), dimension(4) :: elems

  coords = reshape((/  0.0,    0.0, 0.0,   &
                       0.0,  1.203, 0.0,   &
                    -0.934, -0.579, 0.0,   &
                     0.934, -0.579, 0.0/), &
                   (/3,4/))

  elems = (/"C ", "O ", "H ", "H "/)

  call cmlBeginFile(xf, filename, unit=-1)
  call cmlStartCml(xf)
  
  call cmlStartMolecule(xf, title="Formaldehyde")
  call cmlAddAtoms(xf, coords=coords, elements=elems)
  call cmlEndMolecule(xf)

  call cmlFinishFile(xf)

end program test
