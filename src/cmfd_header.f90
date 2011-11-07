module cmfd_header

  implicit none

!===============================================================================
! cmfd is used to store diffusion parameters and other information for CMFD
! analysis.  
!===============================================================================

  type cmfd_obj

    ! array indices([1-x,2-y,3-z,4-g],upper bound)
    integer              :: indices(4)

    ! cross sections
    real(8), allocatable :: totalxs(:,:,:,:)
    real(8), allocatable :: scattxs(:,:,:,:,:)
    real(8), allocatable :: nfissxs(:,:,:,:,:)

    ! diffusion coefficient
    real(8), allocatable :: diffcof(:,:,:,:)

    ! currents
    real(8), allocatable :: currentX(:,:,:,:)
    real(8), allocatable :: currentY(:,:,:,:)
    real(8), allocatable :: currentZ(:,:,:,:)

    ! flux
    real(8), allocatable :: flux(:,:,:,:)

    ! coupling coefficients
    real(8), allocatable :: dtilda(:,:,:,:,:)
    real(8), allocatable :: dhat(:,:,:,:,:)

    ! core albedo boundary conditions
    real(8)              :: albedo(6)

    ! dimensions of mesh cells (xloc,yloc,zloc,[hu,hv,hw])
    real(8), allocatable :: hxyz(:,:,:,:)

    ! source probability distribution
    real(8), allocatable :: sourcepdf(:,:,:,:)

    ! source sites in each mesh box
    real(8), allocatable :: sourcecounts(:,:,:,:)

    ! weight adjustment factors
    real(8), allocatable :: weightfactors(:,:,:,:)

    ! core map for xs association
    integer, allocatable :: coremap(:,:,:)

    ! we may need to add the mesh object
    ! add accumulation of important parameters

  end type cmfd_obj

end module cmfd_header
