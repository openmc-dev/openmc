module cmfd_header

  use constants,  only: CMFD_NOACCEL, ZERO, ONE
  use mesh_header, only: RegularMesh
  use set_header, only: SetInt
  use tally_header, only: TallyContainer
  use timer_header, only: Timer

  implicit none
  private
  public :: allocate_cmfd, deallocate_cmfd

  type, public :: cmfd_type

    ! Indices for problem
    integer :: indices(4)

    ! Albedo boundary condition
    real(8) :: albedo(6)

    ! Core overlay map
    integer, allocatable :: coremap(:,:,:)
    integer, allocatable :: indexmap(:,:)
    integer :: mat_dim = CMFD_NOACCEL

    ! Energy grid
    real(8), allocatable :: egrid(:)

    ! Cross sections
    real(8), allocatable :: totalxs(:,:,:,:)
    real(8), allocatable :: p1scattxs(:,:,:,:)
    real(8), allocatable :: scattxs(:,:,:,:,:)
    real(8), allocatable :: nfissxs(:,:,:,:,:)

    ! Diffusion coefficient
    real(8), allocatable :: diffcof(:,:,:,:)

    ! Current
    real(8), allocatable :: current(:,:,:,:,:)

    ! Flux
    real(8), allocatable :: flux(:,:,:,:)

    ! Coupling coefficients and equivalence parameters
    real(8), allocatable :: dtilde(:,:,:,:,:)
    real(8), allocatable :: dhat(:,:,:,:,:)

    ! Dimensions of mesh cells ([hu,hv,hw],xloc,yloc,zloc)
    real(8), allocatable :: hxyz(:,:,:,:)

    ! Source distributions
    real(8), allocatable :: cmfd_src(:,:,:,:)
    real(8), allocatable :: openmc_src(:,:,:,:)

    ! Source sites in each mesh box
    real(8), allocatable :: sourcecounts(:,:)

    ! Weight adjustment factors
    real(8), allocatable :: weightfactors(:,:,:,:)

    ! Eigenvector/eigenvalue from cmfd run
    real(8), allocatable :: phi(:)
    real(8) :: keff = ZERO

    ! Eigenvector/eigenvalue from adjoint run
    real(8), allocatable :: adj_phi(:)
    real(8) :: adj_keff = ZERO

    ! Residual for neutron balance
    real(8), allocatable :: resnb(:,:,:,:)

    ! Openmc source normalization factor
    real(8) :: norm = ONE

    ! "Shannon entropy" from cmfd fission source
    real(8), allocatable :: entropy(:)

    ! RMS of neutron balance equations
    real(8), allocatable :: balance(:)

    ! RMS deviation of OpenMC and CMFD normalized source
    real(8), allocatable :: src_cmp(:)

    ! Dominance ratio
    real(8), allocatable :: dom(:)

    ! List of CMFD k
    real(8), allocatable :: k_cmfd(:)

    ! Balance keff
    real(8) :: keff_bal

  end type cmfd_type

  ! Main object
  type(cmfd_type), public :: cmfd

  type(RegularMesh), public, pointer :: cmfd_mesh => null()

  ! Pointers for different tallies
  type(TallyContainer), public, pointer :: cmfd_tallies(:) => null()

  ! Timing objects
  type(Timer), public :: time_cmfd      ! timer for whole cmfd calculation
  type(Timer), public :: time_cmfdbuild ! timer for matrix build
  type(Timer), public :: time_cmfdsolve ! timer for solver

  ! Flag for active core map
  logical, public :: cmfd_coremap = .false.

  ! Flag to reset dhats to zero
  logical, public :: dhat_reset = .false.

  ! Flag to activate neutronic feedback via source weights
  logical, public :: cmfd_feedback = .false.

  ! Adjoint method type
  character(len=10), public :: cmfd_adjoint_type = 'physical'

  ! Number of incomplete ilu factorization levels
  integer, public :: cmfd_ilu_levels = 1

  ! Batch to begin cmfd
  integer, public :: cmfd_begin = 1

  ! Tally reset list
  integer, public :: n_cmfd_resets
  type(SetInt), public :: cmfd_reset

  ! Compute effective downscatter cross section
  logical, public :: cmfd_downscatter = .false.

  ! Convergence monitoring
  logical, public :: cmfd_power_monitor = .false.

  ! Cmfd output
  logical, public :: cmfd_write_matrices = .false.

  ! Run an adjoint calculation (last batch only)
  logical, public :: cmfd_run_adjoint = .false.

  ! CMFD run logicals
  logical, public :: cmfd_on = .false.

  ! CMFD display info
  character(len=25), public :: cmfd_display = 'balance'

  ! Estimate of spectral radius of CMFD matrices and tolerances
  real(8), public :: cmfd_spectral = ZERO
  real(8), public :: cmfd_shift = 1.e6
  real(8), public :: cmfd_ktol = 1.e-8_8
  real(8), public :: cmfd_stol = 1.e-8_8
  real(8), public :: cmfd_atoli = 1.e-10_8
  real(8), public :: cmfd_rtoli = 1.e-5_8

contains

!==============================================================================
! ALLOCATE_CMFD allocates all data in of cmfd type
!==============================================================================

  subroutine allocate_cmfd(this, n_batches)

    integer, intent(in)            :: n_batches ! number of batches in calc
    type(cmfd_type), intent(inout) :: this      ! cmfd instance

    integer :: nx  ! number of mesh cells in x direction
    integer :: ny  ! number of mesh cells in y direction
    integer :: nz  ! number of mesh cells in z direction
    integer :: ng  ! number of energy groups

   ! Extract spatial and energy indices from object
    nx = this % indices(1)
    ny = this % indices(2)
    nz = this % indices(3)
    ng = this % indices(4)

    ! Allocate flux, cross sections and diffusion coefficient
    if (.not. allocated(this % flux))       allocate(this % flux(ng,nx,ny,nz))
    if (.not. allocated(this % totalxs))    allocate(this % totalxs(ng,nx,ny,nz))
    if (.not. allocated(this % p1scattxs))  allocate(this % p1scattxs(ng,nx,ny,nz))
    if (.not. allocated(this % scattxs))    allocate(this % scattxs(ng,ng,nx,ny,nz))
    if (.not. allocated(this % nfissxs))    allocate(this % nfissxs(ng,ng,nx,ny,nz))
    if (.not. allocated(this % diffcof))    allocate(this % diffcof(ng,nx,ny,nz))

    ! Allocate dtilde and dhat
    if (.not. allocated(this % dtilde))     allocate(this % dtilde(6,ng,nx,ny,nz))
    if (.not. allocated(this % dhat))       allocate(this % dhat(6,ng,nx,ny,nz))

    ! Allocate dimensions for each box (here for general case)
    if (.not. allocated(this % hxyz))       allocate(this % hxyz(3,nx,ny,nz))

    ! Allocate surface currents
    if (.not. allocated(this % current))    allocate(this % current(12,ng,nx,ny,nz))

    ! Allocate source distributions
    if (.not. allocated(this % cmfd_src)) allocate(this % cmfd_src(ng,nx,ny,nz))
    if (.not. allocated(this % openmc_src)) allocate(this % openmc_src(ng,nx,ny,nz))

    ! Allocate source weight modification vars
    if (.not. allocated(this % sourcecounts)) allocate(this % sourcecounts(ng,nx*ny*nz))
    if (.not. allocated(this % weightfactors)) allocate(this % weightfactors(ng,nx,ny,nz))

    ! Allocate batchwise parameters
    if (.not. allocated(this % entropy)) allocate(this % entropy(n_batches))
    if (.not. allocated(this % balance)) allocate(this % balance(n_batches))
    if (.not. allocated(this % src_cmp)) allocate(this % src_cmp(n_batches))
    if (.not. allocated(this % dom)) allocate(this % dom(n_batches))
    if (.not. allocated(this % k_cmfd)) allocate(this % k_cmfd(n_batches))

    ! Set everthing to 0 except weight multiply factors if feedback isnt on
    this % flux          = ZERO
    this % totalxs       = ZERO
    this % p1scattxs     = ZERO
    this % scattxs       = ZERO
    this % nfissxs       = ZERO
    this % diffcof       = ZERO
    this % dtilde        = ZERO
    this % dhat          = ZERO
    this % hxyz          = ZERO
    this % current       = ZERO
    this % cmfd_src      = ZERO
    this % openmc_src    = ZERO
    this % sourcecounts  = ZERO
    this % weightfactors = ONE
    this % balance       = ZERO
    this % src_cmp       = ZERO
    this % dom           = ZERO
    this % k_cmfd        = ZERO
    this % entropy       = ZERO

  end subroutine allocate_cmfd

!===============================================================================
! DEALLOCATE_CMFD frees all memory of cmfd type
!===============================================================================

  subroutine deallocate_cmfd(this)

    type(cmfd_type), intent(inout) :: this ! cmfd instance

    if (allocated(this % egrid))         deallocate(this % egrid)
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
    if (allocated(this % sourcecounts))  deallocate(this % sourcecounts)
    if (allocated(this % weightfactors)) deallocate(this % weightfactors)
    if (allocated(this % cmfd_src))      deallocate(this % cmfd_src)
    if (allocated(this % openmc_src))    deallocate(this % openmc_src)
    if (allocated(this % balance))       deallocate(this % balance)
    if (allocated(this % src_cmp))       deallocate(this % src_cmp)
    if (allocated(this % dom))           deallocate(this % dom)
    if (allocated(this % k_cmfd))        deallocate(this % k_cmfd)
    if (allocated(this % entropy))       deallocate(this % entropy)
    if (allocated(this % resnb))         deallocate(this % resnb)

  end subroutine deallocate_cmfd

end module cmfd_header
