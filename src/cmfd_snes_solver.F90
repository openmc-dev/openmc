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

    use global, only: cmfd,rank,n_procs_cmfd

    integer              :: k         ! implied do counter
    integer              :: n         ! problem size
    integer              :: row_start ! local row start
    integer              :: row_end   ! local row end
    real(8), pointer     :: xptr(:)   ! solution pointer

    ! set up operator matrices
    call init_M_operator(ctx%loss)
    call init_F_operator(ctx%prod)
    call init_J_operator(jac_prec,ctx)

    ! get problem size
    n = jac_prec%n

    ! create PETSc vectors
    call VecCreateMPI(PETSC_COMM_WORLD,jac_prec%localn,PETSC_DECIDE,resvec,ierr)
    call VecCreateMPI(PETSC_COMM_WORLD,jac_prec%localn,PETSC_DECIDE,xvec,ierr)

    ! get the local dimensions for each process
    call VecGetOwnershipRange(xvec,row_start,row_end,ierr)

    if (rank == n_procs_cmfd - 1) row_end = n

    ! set flux in guess
    call VecSetValues(xvec,row_end-row_start,(/(k,k=row_start,row_end-1)/),  &
   &                  cmfd%phi(row_start+1:row_end),INSERT_VALUES,ierr)
    call VecAssemblyBegin(xvec,ierr)
    call VecAssemblyEnd(xvec,ierr)

    ! set keff in guess
    if (rank == n_procs_cmfd - 1) then
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
    call PCSetType(pc,PCHYPRE,ierr)
    call PCSetFromOptions(pc,ierr)
    call KSPSetFromOptions(ksp,ierr)

    ! create matrix free jacobian
    call MatCreateSNESMF(snes,jac,ierr)

    ! set matrix free finite difference
    call SNESSetJacobian(snes,jac,jac_prec,build_jacobian_matrix,ctx,ierr)

    ! set lags
    call SNESSetLagJacobian(snes,-2,ierr)
    call SNESSetLagPreconditioner(snes,-1,ierr)

    ! set convergence
!   call SNESSetTolerances(snes,1.0e-8_8,1.0e-8_8,1.0e-8_8,50,10000,ierr)

    ! set SNES options
    call SNESSetFromOptions(snes,ierr)

    ! turn off line searching
!   call SNESLineSearchSet(snes,SNESLineSearchNo,PETSC_NULL,ierr)

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

    use global, only: rank,n_procs_cmfd,path_input

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
    real(8)     :: lambda        ! eigenvalue
    real(8)     :: reslamb       ! residual for lambda

    real(8), pointer :: xptr(:)  ! pointer to solution vector
    real(8), pointer :: rptr(:)  ! pointer to residual vector
PetscViewer :: viewer

    ! create operators
    call build_loss_matrix(ctx%loss)
    call build_prod_matrix(ctx%prod)

    ! get problem size
    n = ctx%loss%n 

    ! get pointers to vectors
    call VecGetArrayF90(x,xptr,ierr)
    call VecGetArrayF90(res,rptr,ierr)

    ! create petsc vector for flux
    call VecCreateMPI(PETSC_COMM_WORLD,ctx%loss%localn,PETSC_DECIDE,phi,ierr)
    call VecCreateMPI(PETSC_COMM_WORLD,ctx%loss%localn,PETSC_DECIDE,rphi,ierr)

    ! extract flux and place in petsc vector 
    call VecPlaceArray(phi,xptr,ierr)
    call VecPlaceArray(rphi,rptr,ierr)

    ! extract eigenvalue and broadcast (going to want to make this more general in future)
    if (rank == n_procs_cmfd - 1) lambda = xptr(size(xptr)) 
    call MPI_BCAST(lambda,1,MPI_REAL8,n_procs_cmfd-1,PETSC_COMM_WORLD,ierr)

    ! create new petsc vectors to perform math
    call VecCreateMPI(PETSC_COMM_WORLD,ctx%loss%localn,PETSC_DECIDE,phiM,ierr)

    ! calculate flux part of residual vector
    call MatMult(ctx%loss%M,phi,phiM,ierr)
    call MatMult(ctx%prod%F,phi,rphi,ierr)
    call VecAYPX(rphi,-1.0_8*lambda,phiM,ierr)

    ! set eigenvalue part of residual vector
    call VecDot(phi,phi,reslamb,ierr)

    ! map to ptr
    if (rank == n_procs_cmfd) rptr(size(rptr)) = 0.5_8 - 0.5_8*reslamb

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

    ! write out matrix in binary file (debugging)
    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,trim(path_input)//'residual.bin' &
   &     ,FILE_MODE_WRITE,viewer,ierr)
    call VecView(res,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD,trim(path_input)//'x.bin' &
   &     ,FILE_MODE_WRITE,viewer,ierr)
    call VecView(x,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine compute_nonlinear_residual

!===============================================================================
! EXTRACT_RESULTS
!===============================================================================

  subroutine extract_results()

    use global, only: cmfd,rank,n_procs_cmfd

    integer :: n ! problem size
    integer              :: row_start ! local row start
    integer              :: row_end   ! local row end
    real(8),allocatable  :: mybuf(:)  ! temp buffer
    PetscScalar, pointer :: xptr(:) ! pointer to eigenvector info

    ! get problem size
    n = ctx%loss%n

    ! also allocate in cmfd object
    if (.not. allocated(cmfd%phi)) allocate(cmfd%phi(n))
    if (.not. allocated(mybuf)) allocate(mybuf(n))

    ! get ownership range
    call VecGetOwnershipRange(xvec,row_start,row_end,ierr)

    ! resize the last proc
    if (rank == n_procs_cmfd - 1) row_end = row_end - 1

    ! convert petsc phi_object to cmfd_obj
    call VecGetArrayF90(xvec,xptr,ierr)
    cmfd%phi(row_start+1:row_end) = xptr(1:row_end - row_start) 

    ! reduce result to all 
    mybuf = 0.0_8
    call MPI_ALLREDUCE(cmfd%phi,mybuf,n,MPI_REAL8,MPI_SUM,PETSC_COMM_WORLD,ierr)

    ! move buffer to object and deallocate
    cmfd%phi = mybuf
    if(allocated(mybuf)) deallocate(mybuf)

    ! save eigenvalue
    if(rank == n_procs_cmfd - 1) cmfd%keff = 1.0_8 / xptr(size(xptr))
    call MPI_BCAST(cmfd%keff,1,MPI_REAL8,n_procs_cmfd-1,PETSC_COMM_WORLD,ierr)
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
