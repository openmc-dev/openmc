module cmfd_header

  implicit none
  private
  public :: allocate_cmfd, deallocate_cmfd

!===============================================================================
! cmfd is used to store diffusion parameters and other information for CMFD
! analysis.  
!===============================================================================

  type, public ::  cmfd_obj

    ! array indices([1-x,2-y,3-z,4-g],upper bound)
    integer              :: indices(4)

    ! cross sections
    real(8), allocatable :: totalxs(:,:,:,:)
    real(8), allocatable :: p1scattxs(:,:,:,:)
    real(8), allocatable :: scattxs(:,:,:,:,:)
    real(8), allocatable :: nfissxs(:,:,:,:,:)

    ! diffusion coefficient
    real(8), allocatable :: diffcof(:,:,:,:)

    ! current 
    real(8), allocatable :: current(:,:,:,:,:)

    ! flux
    real(8), allocatable :: flux(:,:,:,:)

    ! coupling coefficients
    real(8), allocatable :: dtilde(:,:,:,:,:)
    real(8), allocatable :: dhat(:,:,:,:,:)

    ! core albedo boundary conditions
    real(8)              :: albedo(6)

    ! dimensions of mesh cells ([hu,hv,hw],xloc,yloc,zloc)
    real(8), allocatable :: hxyz(:,:,:,:)

    ! source probability distribution
    real(8), allocatable :: sourcepdf(:,:,:,:)

    ! source sites in each mesh box
    real(8), allocatable :: sourcecounts(:,:,:,:)

    ! weight adjustment factors
    real(8), allocatable :: weightfactors(:,:,:,:)

    ! core map for no reflector accel
    integer, allocatable :: coremap(:,:,:)
    integer, allocatable :: indexmap(:,:)
    integer :: mat_dim

    ! eigenvector/eigenvalue from cmfd run
    real(8), allocatable :: phi(:)
    real(8) :: keff = 0.0_8

    ! residual for neutron balance
    real(8), allocatable :: resnb(:,:,:,:)

  end type cmfd_obj

contains

!===============================================================================
! ALLOCATE_CMFD allocates all of the space for the cmfd object based on tallies
!===============================================================================

  subroutine allocate_cmfd(this)

    type(cmfd_obj) :: this

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
    if (.not. allocated(this % flux)) allocate(this % flux(ng,nx,ny,nz))
    if (.not. allocated(this % totalxs)) allocate(this % totalxs(ng,nx,ny,nz))
    if (.not. allocated(this % p1scattxs)) allocate(this % p1scattxs(ng,nx,ny,nz))
    if (.not. allocated(this % scattxs)) allocate(this % scattxs(ng,ng,nx,ny,nz))
    if (.not. allocated(this % nfissxs)) allocate(this % nfissxs(ng,ng,nx,ny,nz))
    if (.not. allocated(this % diffcof)) allocate(this % diffcof(ng,nx,ny,nz))

    ! allocate dtilde and dhat
    if (.not. allocated(this % dtilde)) allocate(this % dtilde(6,ng,nx,ny,nz))
    if (.not. allocated(this % dhat)) allocate(this % dhat(6,ng,nx,ny,nz))

    ! allocate dimensions for each box (here for general case)
    if (.not. allocated(this % hxyz)) allocate(this % hxyz(3,nx,ny,nz))

    ! allocate this fission source pdf
    !allocate( this % sourcepdf(ng,nx,ny,nz) )

    ! allocate surface currents
    if (.not. allocated(this % current)) allocate(this % current(12,ng,nx,ny,nz))

  end subroutine allocate_cmfd

!===============================================================================
! DEALLOCATE_CMFD 
!===============================================================================

  subroutine deallocate_cmfd(this)

    type(cmfd_obj) :: this

    ! deallocate cmfd
    if (allocated(this % totalxs))       deallocate(this % totalxs)
    if (allocated(this % p1scattxs))     deallocate(this % p1scattxs)
    if (allocated(this % scattxs))       deallocate(this % scattxs)
    if (allocated(this % nfissxs))       deallocate(this % nfissxs)
    if (allocated(this % diffcof))       deallocate(this % diffcof)
    if (allocated(this % current))       deallocate(this % current)
    if (allocated(this % flux))          deallocate(this % flux)
    if (allocated(this % dtilde))        deallocate(this % dtilde)
    if (allocated(this % dhat))          deallocate(this % dhat)
    if (allocated(this % hxyz))          deallocate(this % hxyz)
    if (allocated(this % coremap))       deallocate(this % coremap)
    if (allocated(this % indexmap))      deallocate(this % indexmap)
    if (allocated(this % phi))           deallocate(this % phi)
    if (allocated(this % sourcepdf))     deallocate(this % sourcepdf)
    if (allocated(this % sourcecounts))  deallocate(this % sourcecounts)
    if (allocated(this % weightfactors)) deallocate(this % weightfactors)

  end subroutine deallocate_cmfd

end module cmfd_header
