module types

  implicit none

  type :: Material
    real(8) :: Sigma_f          ! Macro fission xs
    real(8) :: Sigma_s          ! Macro scattering xs
    real(8) :: Sigma_a          ! Macro absorption xs
    real(8) :: Sigma_t          ! Macro total xs
    real(8) :: Sigma_nf         ! Macro nu-fission xs 
  end type Material
  
  type :: region
    real(8) :: xyz_min(3)       ! Minimum coordinates of region
    real(8) :: xyz_max(3)       ! Maximum coordinates of region
    integer  :: mat             ! Material in region
  end type region

  type :: bank
     real(8) :: xyz(3)          ! location
     integer  :: uid            ! Unique ID
  end type bank
  
  type :: particle
    real(8) :: xyz(3)           ! location
    real(8) :: uvw(3)           ! directional cosines
    real(8) :: wgt              ! particle weight
    integer  :: ijk(3)          ! interval in global mesh
    integer  :: local_ijk(3)    ! interval in local mesh
    integer  :: uid             ! Unique ID
  end type particle
  
  type :: fissionBankPtr
    type(bank), pointer :: ptr  ! Used for sorting fission bank
  end type fissionBankPtr
  
  ! For each basic type (Cell, Surface, and Material), we have
  ! corresponding types CellList, etc that have an extra pointer for a
  ! linked list. As the input file is being read, the list types are
  ! used since they can be dynamically allocated. Once the input file
  ! is complete, the linked list is converted to a normal array.

  type :: Surface
     integer :: id
     integer :: type
     real(8), allocatable :: coeffs(:)
     integer, allocatable :: neighbor_pos(:)
     integer, allocatable :: neighbor_neg(:)
     integer :: bc
  end type Surface

  type :: SurfaceList
     integer :: id
     integer :: type
     real(8), allocatable :: coeffs(:)
     integer, allocatable :: neighbor_pos(:)
     integer, allocatable :: neighbor_neg(:)
     integer :: bc
     type(SurfaceList), pointer :: next
  end type SurfaceList

  type :: Cell
     integer :: id
     integer :: n_items
     integer, allocatable :: boundary_list(:)
     integer :: material
  end type Cell

  type :: CellList
     integer :: id
     integer, allocatable :: boundary_list(:)
     integer :: material
     type(CellList), pointer :: next
  end type CellList

  type :: Neutron
    integer  :: uid    ! Unique ID
    real(8)  :: xyz(3) ! location
    real(8)  :: uvw(3) ! directional cosines
    real(8)  :: wgt    ! particle weight
  end type Neutron


end module types
