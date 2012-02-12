module cmfd_snes_solver

#ifdef PETSC
  use cmfd_loss_operator,     only: loss_operator,init_M_operator,             &
 &                                  build_loss_matrix,destroy_M_operator
  use cmfd_prod_operator,     only: prod_operator,init_F_operator,             &
 &                                  build_prod_matrix,destroy_F_operator
  use cmfd_jacobian_operator, only: jacobian_operator,init_J_operator,         &
 &                                  build_jacobian_matrix,destroy_J_operator,  &
 &                                  operators
  use cmfd_slepc_solver,      only: cmfd_slepc_execute

implicit none

#include <finclude/petsc.h90>

  type(jacobian_operator) :: jac_prec
  type(operators) :: ctx

  Mat         :: jac         ! jacobian matrix
! Mat         :: jac_prec    ! preconditioned jacobian matrix
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
print *,'executing slepc'
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ! call slepc solver 
    call cmfd_slepc_execute()
print *,'initing data'
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ! initialize data
    call init_data()
print *,'initing solver'
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    ! initialize solver
    call init_solver()
print *,'solving system'
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
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

    integer              :: k         ! implied do counter
    integer              :: n         ! problem size
    integer              :: row_start ! local row start
    integer              :: row_end   ! local row end
    integer, allocatable :: nnzv(:)   ! vector of number of nonzeros in jacobian
    real(8), pointer     :: xptr(:)   ! solution pointer

    ! set up operator matrices
    call init_M_operator(ctx%loss)
    call init_F_operator(ctx%prod)
    call init_J_operator(jac_prec,ctx)

    ! get problem size
    n = jac_prec%n - 1 

    ! create PETSc vectors
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n+1,resvec,ierr)
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n+1,xvec,ierr)

    ! get the local dimensions for each process
    call VecGetOwnershipRange(xvec,row_start,row_end,ierr)

    if (row_end == n + 1) row_end = n

    ! set flux in guess
    call VecSetValues(xvec,row_end-row_start-1,(/(k,k=row_start,row_end-1)/),  &
   &                  cmfd%phi(row_start+1:row_end),INSERT_VALUES,ierr)
    call VecAssemblyBegin(xvec,ierr)
    call VecAssemblyEnd(xvec,ierr)

    ! set keff in guess
    if (row_end == n) then
      call VecGetArrayF90(xvec,xptr,ierr)
      xptr(size(xptr)) = 1.0_8/cmfd%keff
      call VecRestoreArrayF90(xvec,xptr,ierr)
    end if

  end subroutine init_data

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver()

    character(LEN=20) :: snestype,ksptype,pctype

    ! turn on mf_operator option
    call PetscOptionsSetValue("-snes_mf_operator","TRUE",ierr)

    ! create SNES context
    call SNESCreate(PETSC_COMM_WORLD,snes,ierr)

    ! set the residual function
    call SNESSetFunction(snes,resvec,compute_nonlinear_residual,PETSC_NULL,ierr)

    ! set GMRES solver
    call SNESGetKSP(snes,ksp,ierr)
    call KSPSetType(ksp,KSPGMRES,ierr)

    ! set preconditioner
    call KSPGetPC(ksp,pc,ierr)
!   call PCSetType(pc,PCILU,ierr)
!   call PCFactorSetLevels(pc,25,ierr)
    call PCSetType(pc,PCBJACOBI,ierr)
    call KSPSetFromOptions(ksp,ierr)

    ! create matrix free jacobian
    call MatCreateSNESMF(snes,jac,ierr)

    ! set matrix free finite difference
    call SNESSetJacobian(snes,jac,jac_prec,build_jacobian_matrix,ctx,ierr)

    ! set lags
    call SNESSetLagJacobian(snes,-2,ierr)
    call SNESSetLagPreconditioner(snes,-1,ierr)

    ! set convergence
    call SNESSetTolerances(snes,1.0e-8_8,1.0e-8_8,1.0e-8_8,50,10000,ierr)

    ! set SNES options
    call SNESSetFromOptions(snes,ierr)

    ! turn off line searching
    call SNESLineSearchSet(snes,SNESLineSearchNo,PETSC_NULL,ierr)

    ! get all types and print
    call SNESGetType(snes,snestype,ierr)
    call KSPGetType(ksp,ksptype,ierr)
    call PCGetType(pc,pctype,ierr)

    ! display information to user
!   write(*,'(/,A)') 'SNES SOLVER OPTIONS:'
!   write(*,*) '---------------------'
!   write(*,*) 'SNES TYPE IS: ',snestype
!   write(*,*) 'KSP TYPE IS: ',ksptype
!   write(*,*) 'PC TYPE IS: ',pctype

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

    ! local variables
    Vec         :: phi           ! flux vector
    Vec         :: rphi          ! flux part of residual
    Vec         :: phiM          ! M part of residual flux calc
    integer     :: n             ! problem size
    integer     :: row_start     ! local row start
    integer     :: row_end       ! local row end
    integer     :: n_procs       ! number of processors
    real(8)     :: lambda        ! eigenvalue
    real(8)     :: reslamb       ! residual for lambda

    real(8), pointer :: xptr(:)  ! pointer to solution vector
    real(8), pointer :: rptr(:)  ! pointer to residual vector
print *,'In residual'

    ! get number of processors
    call MPI_COMM_SIZE(MPI_COMM_WORLD,n_procs,ierr)

    ! get problem size
    call VecGetSize(x,n,ierr)
    n = n - 1 ! subtract off last row 

    ! get the local dimensions for each process
    call VecGetOwnershipRange(x,row_start,row_end,ierr)

    ! get pointers to vectors
    call VecGetArrayF90(x,xptr,ierr)
    call VecGetArrayF90(res,rptr,ierr)

    ! create petsc vector for flux
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,phi,ierr)
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,rphi,ierr)

    ! extract flux and place in petsc vector 
    call VecPlaceArray(phi,xptr,ierr)
    call VecPlaceArray(rphi,rptr,ierr)

    ! extract eigenvalue and broadcast (going to want to make this more general in future)
    if (row_end == n+1) then
      lambda = xptr(size(xptr))
    end if
    call MPI_BCAST(lambda,1,MPI_REAL8,n_procs-1,MPI_COMM_WORLD,ierr)

    ! create operators
    call build_loss_matrix(ctx%loss)
    call build_prod_matrix(ctx%prod)

    ! create new petsc vectors to perform math
    call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,phiM,ierr)

    ! calculate flux part of residual vector
    call MatMult(ctx%loss%M,phi,phiM,ierr)
    call MatMult(ctx%prod%F,phi,rphi,ierr)
    call VecAYPX(rphi,-1.0_8*lambda,phiM,ierr)

    ! set eigenvalue part of residual vector
    call VecDot(phi,phi,reslamb,ierr)

    ! map to ptr
    if (row_end == n+1) rptr(size(rptr)) = 0.5_8 - 0.5_8*reslamb

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

!===============================================================================
! EXTRACT_RESULTS
!===============================================================================

  subroutine extract_results()

    use global, only: cmfd

    integer :: n ! problem size
    PetscScalar, pointer :: xptr(:) ! pointer to eigenvector info

    ! get problem size
    n = ctx%loss%n

    ! also allocate in cmfd object
    if (.not. allocated(cmfd%phi)) allocate(cmfd%phi(n))

    ! convert petsc phi_object to cmfd_obj
    call VecGetArrayF90(xvec,xptr,ierr)
    cmfd%phi = xptr(1:n)

    ! save eigenvalue
    cmfd%keff = 1.0_8 / xptr(n+1)
    call VecRestoreArrayF90(xvec,xptr,ierr)

  end subroutine extract_results

!==============================================================================
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize data objects
    call destroy_M_operator(ctx%loss)
    call destroy_F_operator(ctx%prod)
    call destroy_J_operator(jac_prec)
    call VecDestroy(xvec,ierr)
    call VecDestroy(resvec,ierr)

    ! finalize solver objects
    call SNESDestroy(snes,ierr)

  end subroutine finalize

#endif

end module cmfd_snes_solver
