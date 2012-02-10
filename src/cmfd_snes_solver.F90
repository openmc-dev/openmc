module cmfd_snes_solver

#ifdef PETSC
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
  Mat         :: jac_prec    ! preconditioned jacobian matrix
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
    integer, allocatable :: nnzv(:) ! vector of number of nonzeros in jacobian
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

   ! allocate and create number of nonzeros vector
   if (.not. allocated(nnzv)) allocate(nnzv(n+1))
   nnzv = loss%nnz + 1
   nnzv(n+1) = n+1

   ! create the preconditioner matrix
   call MatCreateSeqAIJ(PETSC_COMM_SELF,n+1,n+1,PETSC_NULL_INTEGER,          &
  &                     nnzv,jac_prec,ierr)
   call MatSetOption(jac_prec,MAT_NEW_NONZERO_LOCATIONS,PETSC_TRUE,ierr)
   call MatSetOption(jac_prec,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE,ierr)
   call MatSetOption(jac_prec,MAT_USE_HASH_TABLE,PETSC_TRUE,ierr)

  end subroutine init_data

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver()

    character(LEN=20) :: snestype,ksptype,pctype

    ! turn on mf_operator option
    call PetscOptionsSetValue("-snes_mf_operator","TRUE",ierr)

    ! create SNES context
    call SNESCreate(PETSC_COMM_SELF,snes,ierr)

    ! set the residual function
    call SNESSetFunction(snes,resvec,compute_nonlinear_residual,PETSC_NULL,ierr)

    ! set GMRES solver
    call SNESGetKSP(snes,ksp,ierr)
    call KSPSetType(ksp,KSPGMRES,ierr)

    ! set preconditioner
    call KSPGetPC(ksp,pc,ierr)
    call PCSetType(pc,PCILU,ierr)
    call PCFactorSetLevels(pc,25,ierr)
    call KSPSetFromOptions(ksp,ierr)

    ! create matrix free jacobian
    call MatCreateSNESMF(snes,jac,ierr)

    ! set matrix free finite difference
    call SNESSetJacobian(snes,jac,jac_prec,compute_jacobian,PETSC_NULL,ierr)

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

!===============================================================================
! COMPUTE_JACOBIAN 
!===============================================================================

  subroutine compute_jacobian(snes,x,jac,jac_prec,flag,user,ierr)

     ! formal variables
     SNES          :: snes      ! the snes context
     Vec           :: x         ! the solution vector
     Mat           :: jac       ! the jacobian matrix
     Mat           :: jac_prec  ! the jacobian preconditioner
     MatStructure  :: flag      ! not used
     integer       :: user(*)   ! not used
     integer       :: ierr      ! petsc error flag

     ! local varibles
     Vec                  :: phi      ! flux vector
     Vec                  :: source   ! source vector
     integer              :: n        ! problem size
     integer              :: k        ! implied do loop counter
     integer              :: ncols    ! number of nonzeros in cols
     integer              :: irow     ! row counter
     integer, allocatable :: nnzv(:)  ! vector of number of nonzeros for jac
     integer, allocatable :: cols(:)  ! vector of column numbers
     real(8)              :: lambda   ! eigenvalue
     real(8), pointer     :: xptr(:)  ! pointer to solution vector
     real(8), pointer     :: sptr(:)  ! pointer to source vector
     real(8), allocatable :: vals(:)  ! vector of row values

     ! create operators
     call build_loss_matrix(loss)
     call build_prod_matrix(prod)

     ! get problem size
     n = loss%n

     ! allocate cols and rows and initialize to zero
     if (.not. allocated(cols)) allocate(cols(loss%nnz))
     if (.not. allocated(vals)) allocate(vals(loss%nnz))
     cols = 0
     vals = 0.0_8

     ! get pointers to residual vector 
     call VecGetArrayF90(x,xptr,ierr)

     ! create petsc vector for flux
     call VecCreate(PETSC_COMM_SELF,phi,ierr)
     call VecSetSizes(phi,PETSC_DECIDE,n,ierr)
     call VecSetFromOptions(phi,ierr)

     ! extract flux and eigenvalue 
     call VecPlaceArray(phi,xptr,ierr)
     lambda = xptr(n+1)

     ! compute math (M-lambda*F) M is overwritten here
     call MatAXPY(loss%M,-1.0_8*lambda,prod%F,SUBSET_NONZERO_PATTERN,ierr)

     ! create tmp petsc vector for source
     call VecCreate(PETSC_COMM_SELF,source,ierr)
     call VecSetSizes(source,PETSC_DECIDE,n,ierr)
     call VecSetFromOptions(source,ierr)

     ! perform math (-F*phi --> source)
     call MatMult(prod%F,phi,source,ierr)
     call VecScale(source,-1.0_8,ierr)

     ! get pointer to source
     call VecGetArrayF90(source,sptr,ierr)

     ! begin loop to insert things into matrix
     do irow = 0,n-1

       ! get row of matrix
       call MatGetRow(loss%M,irow,ncols,cols,vals,ierr)

       ! set that row to Jacobian matrix
       call MatSetValues(jac_prec,1,irow,ncols,cols(1:ncols),vals,INSERT_VALUES,ierr)

       ! restore the row
       call MatRestoreRow(loss%M,irow,ncols,cols,vals,ierr)

       ! insert source value
       call MatSetValue(jac_prec,irow,n,sptr(irow+1),INSERT_VALUES,ierr)

     end do

     ! set values in last row of matrix
     call MatSetValues(jac_prec,1,n,n,(/(k,k=0,n-1)/),xptr(1:n),INSERT_VALUES,ierr)
     call MatSetValue(jac_prec,n,n,1.0_8,INSERT_VALUES,ierr)

     ! assemble matrix
     call MatAssemblyBegin(jac_prec,MAT_FINAL_ASSEMBLY,ierr)
     call MatAssemblyEnd(jac_prec,MAT_FINAL_ASSEMBLY,ierr)
     call MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY,ierr)
     call MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY,ierr)

     ! reset all vectors
     call VecResetArray(phi,ierr)

     ! restore all vectors
     call VecRestoreArrayF90(x,xptr,ierr)
     call VecRestoreArrayF90(source,sptr,ierr)

     ! destroy all temporary objects
     call VecDestroy(phi,ierr)
     call VecDestroy(source,ierr)

     ! deallocate all temporary space
     if (allocated(cols)) deallocate(cols)
     if (allocated(vals)) deallocate(vals)

  end subroutine compute_jacobian

!===============================================================================
! EXTRACT_RESULTS
!===============================================================================

  subroutine extract_results()

    use global, only: cmfd

    integer :: n ! problem size
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

#endif

end module cmfd_snes_solver
