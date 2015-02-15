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
  
  call cmlAddMolecule(xf, coords=coords, elements=elems, title="Formaldehyde", &
    atomIds=(/"c1", "o1", "h1", "h2"/), bondAtom1Refs=(/"c1", "c1", "c1"/), &
    bondAtom2Refs=(/"o1", "h1", "h1"/), bondOrders=(/"D", "S", "S"/))

  call cmlFinishFile(xf)

end program test
