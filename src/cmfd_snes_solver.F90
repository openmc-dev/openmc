module cmfd_snes_solver

  use cmfd_loss_operator, only: loss_operator,init_M_operator,                 &
 &                        build_loss_matrix,destroy_M_operator
  use cmfd_prod_operator, only: prod_operator,init_F_operator,                 &
 &                        build_prod_matrix,destroy_F_operator
  use cmfd_slepc_solver,  only: cmfd_slepc_execute

implicit none

#include <finclude/petsc.h90>


  type(loss_operator) :: loss
  type(prod_operator) :: prod

  Mat         :: jac         ! jacobian matrix
  Vec         :: resvec      ! residual vector
  Vec         :: xvec        ! results
  KSP         :: ksp         ! linear solver context
  PC          :: pc          ! preconditioner
  SNES        :: snes        ! nonlinear solver context
  integer     :: ierr        ! error flag

contains

!===============================================================================
! CMFD_SNES_EXECUTE
!===============================================================================

  subroutine cmfd_snes_execute()

    ! call slepc solver 
    call cmfd_slepc_execute()

    ! initialize data
    call init_data()

    ! initialize solver
    call init_solver()

    ! solve the system
    call SNESSolve(snes,PETSC_NULL,xvec,ierr)

    ! extracts results to cmfd object
    call extract_results()

    ! deallocate all slepc data
    call finalize()

  end subroutine cmfd_snes_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

    use global, only: cmfd

    integer :: k  ! implied do counter
    integer :: n  ! problem size
    real(8), pointer :: xptr(:) ! solution pointer

    ! set up matrices
    call init_M_operator(loss)
    call init_F_operator(prod)

    ! get problem size
    n = loss%n

    ! create PETSc vectors
    call VecCreate(PETSC_COMM_SELF,resvec,ierr)
    call VecSetSizes(resvec,PETSC_DECIDE,n+1,ierr)
    call VecSetFromOptions(resvec,ierr)
    call VecCreate(PETSC_COMM_SELF,xvec,ierr)
    call VecSetSizes(xvec,PETSC_DECIDE,n+1,ierr)
    call VecSetFromOptions(xvec,ierr)

    ! set flux in guess
    call VecSetValues(xvec,n,(/(k,k=0,n-1)/),cmfd%phi,INSERT_VALUES,ierr)
    call VecAssemblyBegin(xvec,ierr)
    call VecAssemblyEnd(xvec,ierr)

    ! set keff in guess
    call VecGetArrayF90(xvec,xptr,ierr)
    xptr(n+1) = 1.0_8/cmfd%keff
    call VecRestoreArrayF90(xvec,xptr,ierr)

  end subroutine init_data

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver()

    ! create SNES context
    call SNESCreate(PETSC_COMM_SELF,snes,ierr)

    ! set the residual function
    call SNESSetFunction(snes,resvec,compute_nonlinear_residual,PETSC_NULL,ierr)

    ! set GMRES solver
    call SNESGetKSP(snes,ksp,ierr)
    call KSPSetType(ksp,KSPGMRES,ierr)

    ! set preconditioner
    call KSPGetPC(ksp,pc,ierr)
    call PCSetType(pc,PCNONE,ierr)

    ! set matrix free finite difference
    call MatCreateSNESMF(snes,jac,ierr)
    call SNESSetJacobian(snes,jac,jac,MatMFFDComputeJacobian,PETSC_NULL,ierr)

    ! set SNES options
    call SNESSetFromOptions(snes,ierr)

  end subroutine init_solver

!===============================================================================
! COMPUTE_NONLINEAR_RESIDUAL
!===============================================================================

  subroutine compute_nonlinear_residual(snes,x,res,ierr)

    ! arguments
    SNES        :: snes          ! nonlinear solver context
    Vec         :: x             ! independent vector
    Vec         :: res           ! residual vector
    integer     :: ierr          ! error flag

    Vec         :: phi           ! flux vector
    Vec         :: rphi          ! flux part of residual
    Vec         :: phiM          ! M part of residual flux calc
    integer     :: n             ! problem size
    real(8)     :: lambda        ! eigenvalue
    real(8)     :: reslamb       ! residual for lambda

    real(8), pointer :: xptr(:)  ! pointer to solution vector
    real(8), pointer :: rptr(:)  ! pointer to residual vector

    ! get problem size
    n = loss%n

    ! get pointers to vectors
    call VecGetArrayF90(x,xptr,ierr)
    call VecGetArrayF90(res,rptr,ierr)

    ! create petsc vector for flux
    call VecCreate(PETSC_COMM_SELF,phi,ierr)
    call VecSetSizes(phi,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(phi,ierr)
    call VecCreate(PETSC_COMM_SELF,rphi,ierr)
    call VecSetSizes(rphi,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(rphi,ierr)

    ! extract flux and place in petsc vector 
    call VecPlaceArray(phi,xptr,ierr)
    call VecPlaceArray(rphi,rptr,ierr)

    ! extract eigenvalue
    lambda = xptr(n+1)

    ! create operators
    call build_loss_matrix(loss)
    call build_prod_matrix(prod)

    ! create new petsc vectors to perform math
    call VecCreate(PETSC_COMM_SELF,phiM,ierr)
    call VecSetSizes(phiM,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(phiM,ierr)

    ! calculate flux part of residual vector
    call MatMult(loss%M,phi,phiM,ierr)
    call MatMult(prod%F,phi,rphi,ierr)
    call VecAYPX(rphi,-1.0_8*lambda,phiM,ierr)

    ! set eigenvalue part of residual vector
    call VecDot(phi,phi,reslamb,ierr)

    ! map to ptr
    rptr(n+1) = 0.5_8 - 0.5_8*reslamb

    ! reset arrays that are not used
    call VecResetArray(phi,ierr)
    call VecResetArray(rphi,ierr)

    ! restore arrays for residual and solution
    call VecRestoreArrayF90(x,xptr,ierr)
    call VecRestoreArrayF90(res,rptr,ierr)

    ! destroy all temp vectors
    call VecDestroy(phi,ierr)
    call VecDestroy(phiM,ierr)
    call VecDestroy(rphi,ierr)

  end subroutine compute_nonlinear_residual

!==============================================================================
! EXTRACT_RESULTS
!==============================================================================

  subroutine extract_results()

    use global, only: cmfd

    integer :: n ! problem size
    PetscViewer :: viewer
    PetscScalar, pointer :: xptr(:) ! pointer to eigenvector info

    ! get problem size
    n = loss%n

    ! also allocate in cmfd object
    if (.not. allocated(cmfd%phi)) allocate(cmfd%phi(n))

    ! convert petsc phi_object to cmfd_obj
    call VecGetArrayF90(xvec,xptr,ierr)
    cmfd%phi = xptr(1:n)

    ! save eigenvalue
    cmfd%keff = 1.0_8 / xptr(n+1)
    call VecRestoreArrayF90(xvec,xptr,ierr)

     ! write out results
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'residual.bin',FILE_MODE_WRITE, &
                               viewer,ierr)
    call VecView(resvec,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine extract_results

!==============================================================================
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize data objects
    call destroy_M_operator(loss)
    call destroy_F_operator(prod)
    call VecDestroy(xvec,ierr)
    call VecDestroy(resvec,ierr)

    ! finalize solver objects
    call SNESDestroy(snes,ierr)

  end subroutine finalize

end module cmfd_snes_solver
