module cmfd_slepc_solver

  use cmfd_loss_operator, only: loss_operator,init_M_operator,                 &
 &                        build_loss_matrix,destroy_M_operator
  use cmfd_prod_operator, only: prod_operator,init_F_operator,                 &
 &                        build_prod_matrix,destroy_F_operator

  implicit none

#include <finclude/petsc.h90>
#include <finclude/slepcsys.h>
#include <finclude/slepceps.h>

  type(loss_operator) :: loss
  type(prod_operator) :: prod

  Vec         :: phi               ! eigenvector
  EPS         :: eps               ! slepc eigenvalue object
  ST          :: st                ! slepc spectral trans object
  KSP         :: ksp               ! linear solver object
  PC          :: pc                ! preconditioner object
  integer     :: ierr              ! error flag
  real(8)     :: keff              ! the converged eigenvalue

contains

!===============================================================================
! CMFD_SLEPC_EXECUTE
!===============================================================================

  subroutine cmfd_slepc_execute()

    ! initialize data
    call init_data()

    ! initialize solver
    call init_solver()

    ! build operators
    call build_loss_matrix(loss)
    call build_prod_matrix(prod)

    ! set operators to EPS object
    call EPSSetOperators(eps,prod%F,loss%M,ierr)

    ! solve the system
    call EPSSolve(eps,ierr)

    ! extracts results to cmfd object
    call extract_results()

    ! deallocate all slepc data
    call finalize()

  end subroutine cmfd_slepc_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

    integer :: n ! problem size

    ! set up matrices
    call init_M_operator(loss)
    call init_F_operator(prod)

    ! get problem size
    n = loss%n

    ! set up eigenvector
    call VecCreate(PETSC_COMM_SELF,phi,ierr)
    call VecSetSizes(phi,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(phi,ierr)

  end subroutine init_data

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver()

    character(LEN=20) :: epstype,sttype,ksptype,pctype

    ! create EPS Object
    call EPSCreate(PETSC_COMM_SELF,eps,ierr)
    call EPSSetProblemType(eps,EPS_GNHEP,ierr)
    call EPSSetType(eps,EPSPOWER,ierr)
    call EPSSetFromOptions(eps,ierr)
    call EPSSetWhichEigenpairs(eps,EPS_LARGEST_MAGNITUDE,ierr)

    ! get ST, KSP and PC objects
    call EPSGetST(eps,st,ierr)
    call STGetKSP(st,ksp,ierr)
    call KSPGetPC(ksp,pc,ierr)

    ! set GMRES default
    call KSPSetType(ksp,KSPGMRES,ierr)

    ! set precursor type
    call PCSetType(pc,PCILU,ierr)
    call PCFactorSetLevels(pc,4,ierr)
    call PCSetFromOptions(pc,ierr)

    ! get all types and print
    call EPSGetType(eps,epstype,ierr)
    call STGetType(st,sttype,ierr)
    call KSPGetType(ksp,ksptype,ierr)
    call PCGetType(pc,pctype,ierr)

    ! display information to user
!   write(*,'(/,A)') 'SLEPC SOLVER OPTIONS:'
!   write(*,*) '---------------------'
!   write(*,*) 'EPS TYPE IS: ',epstype
!   write(*,*) 'ST TYPE IS: ',sttype
!   write(*,*) 'KSP TYPE IS: ',ksptype
!   write(*,*) 'PC TYPE IS: ',pctype

  end subroutine init_solver

!==============================================================================
! EXTRACT_RESULTS
!==============================================================================

  subroutine extract_results()

    use global, only: cmfd

    integer              :: n         ! problem size
    integer              :: i_eig = 0 ! eigenvalue to extract
    PetscViewer          :: viewer    ! petsc output object
    PetscScalar, pointer :: phi_v(:)  ! pointer to eigenvector info
 
    ! get problem size
    n = loss%n

    ! also allocate in cmfd object
    if (.not. allocated(cmfd%phi)) allocate(cmfd%phi(n))

    ! extract run information
    call EPSGetEigenpair(eps,i_eig,keff,PETSC_NULL,phi,PETSC_NULL_OBJECT,ierr)

    ! convert petsc phi_object to cmfd_obj
    call VecGetArrayF90(phi,phi_v,ierr)
    cmfd%phi = phi_v
    call VecRestoreArrayF90(phi,phi_v,ierr)

    ! save eigenvalue
    cmfd%keff = keff

    ! write out results
    call PetscViewerBinaryOpen(PETSC_COMM_SELF,'fluxvec.bin',FILE_MODE_WRITE,  &
   &                           viewer,ierr)
    call VecView(phi,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine extract_results

!==============================================================================
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize data objects
    call destroy_M_operator(loss)
    call destroy_F_operator(prod)
    call VecDestroy(phi,ierr)

    ! finalize solver objects
    call EPSDestroy(eps,ierr)

  end subroutine finalize

end module cmfd_slepc_solver
