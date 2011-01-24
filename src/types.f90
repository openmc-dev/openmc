module types

  implicit none

  type :: Surface
     integer :: uid
     integer :: type
     real(8), allocatable :: coeffs(:)
     integer, allocatable :: neighbor_pos(:)
     integer, allocatable :: neighbor_neg(:)
     integer :: bc
  end type Surface

  type :: Cell
     integer :: uid
     integer :: n_items
     integer, allocatable :: boundary_list(:)
     integer :: material
  end type Cell

  type :: Neutron
    integer  :: uid     ! Unique ID
    real(8)  :: xyz(3)  ! location
    real(8)  :: uvw(3)  ! directional cosines
    integer  :: cell    ! current cell
    integer  :: surface ! current surface
    real(8)  :: wgt     ! particle weight
    logical  :: alive   ! is particle alive?
  end type Neutron

  type :: Bank
     integer :: uid    ! Unique ID
     real(8) :: xyz(3) ! Location of bank particle
  end type Bank

  type :: Isotope
     integer :: uid     ! unique identifier
     integer :: zaid    ! ZAID, e.g. 92235
     character(3) :: xs ! cross section identifier, e.g. 70c
     real(8) :: density ! density in atom/b-cm
  end type Isotope

  type :: Material
     integer :: uid
     integer :: n_isotopes
     type(Isotope), allocatable :: isotopes(:)
  end type Material

  type :: ExtSource
     integer :: type                    ! type of source, e.g. 'box' or 'cell'
     real(8), allocatable :: values(:)  ! values for particular source type
  end type ExtSource

end module types
