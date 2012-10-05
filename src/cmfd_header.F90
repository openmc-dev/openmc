module cmfd_header 

#ifdef HDF5
  use hdf5
#endif

implicit none
private
public :: allocate_cmfd, deallocate_cmfd

  type, public :: cmfd_type

    ! indices for problem
    integer :: indices(4)

    ! albedo boundary condition
    real(8) :: albedo(6)

    ! core overlay map
    integer, allocatable :: coremap(:,:,:)
    integer, allocatable :: indexmap(:,:)
    integer :: mat_dim = 9999

    ! energy grid
    real(8), allocatable :: egrid(:)

    ! cross sections
    real(8), allocatable :: totalxs(:,:,:,:)
    real(8), allocatable :: p1scattxs(:,:,:,:)
    real(8), allocatable :: scattxs(:,:,:,:,:)
    real(8), allocatable :: nfissxs(:,:,:,:,:)

    ! diffusion coefficient
    real(8), allocatable :: diffcof(:,:,:,:)
    real(8), allocatable :: diffusion(:,:,:,:)

    ! current
    real(8), allocatable :: current(:,:,:,:,:)

    ! flux
    real(8), allocatable :: flux(:,:,:,:)

    ! coupling coefficients and equivalence parameters
    real(8), allocatable :: dtilde(:,:,:,:,:)
    real(8), allocatable :: dhat(:,:,:,:,:)

    ! dimensions of mesh cells ([hu,hv,hw],xloc,yloc,zloc)
    real(8), allocatable :: hxyz(:,:,:,:)

    ! source distributions
    real(8), allocatable :: cmfd_src(:,:,:,:)
    real(8), allocatable :: openmc_src(:,:,:,:)

    ! source sites in each mesh box
    real(8), allocatable :: sourcecounts(:,:,:,:)

    ! weight adjustment factors 
    real(8), allocatable :: weightfactors(:,:,:,:)

    ! eigenvector/eigenvalue from cmfd run
    real(8), allocatable :: phi(:)
    real(8) :: keff = 0.0_8

    ! eigenvector/eigenvalue from adjoint run
    real(8), allocatable :: adj_phi(:)
    real(8) :: adj_keff = 0.0_8

    ! residual for neutron balance
    real(8), allocatable :: resnb(:,:,:,:)

    ! openmc source normalization factor
    real(8) :: norm = 1.0_8

    ! Shannon entropy from cmfd fission source
    real(8) :: entropy 

  end type cmfd_type

contains

!==============================================================================
! ALLOCATE_CMFD
!==============================================================================

  subroutine allocate_cmfd(this)

    type(cmfd_type) :: this 

    integer :: nx  ! number of mesh cells in x direction
    integer :: ny  ! number of mesh cells in y direction
    integer :: nz  ! number of mesh cells in z direction
    integer :: ng  ! number of energy groups

   ! extract spatial and energy indices from object
    nx = this % indices(1)
    ny = this % indices(2)
    nz = this % indices(3)
    ng = this % indices(4)

    ! allocate flux, cross sections and diffusion coefficient
    if (.not. allocated(this % flux))       allocate(this % flux(ng,nx,ny,nz))
    if (.not. allocated(this % totalxs))    allocate(this % totalxs(ng,nx,ny,nz))
    if (.not. allocated(this % p1scattxs))  allocate(this % p1scattxs(ng,nx,ny,nz))
    if (.not. allocated(this % scattxs))    allocate(this % scattxs(ng,ng,nx,ny,nz))
    if (.not. allocated(this % nfissxs))    allocate(this % nfissxs(ng,ng,nx,ny,nz))
    if (.not. allocated(this % diffcof))    allocate(this % diffcof(ng,nx,ny,nz))
    if (.not. allocated(this % diffusion))  allocate(this%diffusion(ng,nx,ny,nz))

    ! allocate dtilde and dhat 
    if (.not. allocated(this % dtilde))     allocate(this % dtilde(6,ng,nx,ny,nz))
    if (.not. allocated(this % dhat))       allocate(this % dhat(6,ng,nx,ny,nz))

    ! allocate dimensions for each box (here for general case)
    if (.not. allocated(this % hxyz))       allocate(this % hxyz(3,nx,ny,nz))

    ! allocate surface currents
    if (.not. allocated(this % current))    allocate(this % current(12,ng,nx,ny,nz))

    ! allocate source distributions
    if (.not. allocated(this % cmfd_src)) allocate(this % cmfd_src(ng,nx,ny,nz))
    if (.not. allocated(this % openmc_src)) allocate(this % openmc_src(ng,nx,ny,nz))

    ! allocate source weight modification vars
    if (.not. allocated(this % sourcecounts)) allocate(this % sourcecounts(ng,nx,ny,nz))
    if (.not. allocated(this % weightfactors)) allocate(this % weightfactors(ng,nx,ny,nz))

    ! set everthing to 0 except weight multiply factors if feedback isnt on
    this % flux          = 0.0_8
    this % totalxs       = 0.0_8
    this % p1scattxs     = 0.0_8
    this % scattxs       = 0.0_8
    this % nfissxs       = 0.0_8
    this % diffcof       = 0.0_8
    this % diffusion     = 0.0_8
    this % dtilde        = 0.0_8
    this % dhat          = 0.0_8
    this % hxyz          = 0.0_8
    this % current       = 0.0_8
    this % cmfd_src      = 0.0_8
    this % openmc_src    = 0.0_8
    this % sourcecounts  = 0.0_8
    this % weightfactors = 1.0_8

  end subroutine allocate_cmfd

!===============================================================================
! DEALLOCATE_CMFD 
!===============================================================================

  subroutine deallocate_cmfd(this)

    type(cmfd_type) :: this

    if (allocated(this % egrid))         deallocate(this % egrid)
    if (allocated(this % totalxs))       deallocate(this % totalxs)
    if (allocated(this % p1scattxs))     deallocate(this % p1scattxs)
    if (allocated(this % scattxs))       deallocate(this % scattxs)
    if (allocated(this % nfissxs))       deallocate(this % nfissxs)
    if (allocated(this % diffcof))       deallocate(this % diffcof)
    if (allocated(this % diffusion))     deallocate(this % diffusion)
    if (allocated(this % current))       deallocate(this % current)
    if (allocated(this % flux))          deallocate(this % flux)
    if (allocated(this % dtilde))        deallocate(this % dtilde)
    if (allocated(this % dhat))          deallocate(this % dhat)
    if (allocated(this % hxyz))          deallocate(this % hxyz)
    if (allocated(this % coremap))       deallocate(this % coremap)
    if (allocated(this % indexmap))      deallocate(this % indexmap)
    if (allocated(this % phi))           deallocate(this % phi)
    if (allocated(this % sourcecounts))  deallocate(this % sourcecounts)
    if (allocated(this % weightfactors)) deallocate(this % weightfactors)
    if (allocated(this % cmfd_src))      deallocate(this % cmfd_src)
    if (allocated(this % openmc_src))    deallocate(this % openmc_src)

  end subroutine deallocate_cmfd

end module cmfd_header
