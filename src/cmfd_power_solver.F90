module cmfd_power_solver

#ifdef PETSC
  use cmfd_loss_operator, only: loss_operator,init_M_operator,                 &
 &                        build_loss_matrix,destroy_M_operator
  use cmfd_prod_operator, only: prod_operator,init_F_operator,                 &
 &                        build_prod_matrix,destroy_F_operator


  implicit none
  private
  public :: cmfd_power_execute

#include <finclude/petsc.h90>

  type(loss_operator) :: loss
  type(prod_operator) :: prod

  Vec         :: phi_n      ! new flux eigenvector
  Vec         :: phi_o      ! old flux eigenvector
  Vec         :: S_n        ! new source vector
  Vec         :: S_o        ! old source vector
  KSP         :: krylov     ! krylov solver
  PC          :: prec       ! preconditioner for krylov
  integer     :: ierr       ! error flag
  real(8)     :: k_n        ! new k-eigenvalue
  real(8)     :: k_o        ! old k-eigenvalue
  logical     :: iconv      ! did the problem converged

contains

!===============================================================================
! CMFD_POWER_EXECUTE
!===============================================================================

  subroutine cmfd_power_execute()

    ! initialize solver
    call init_solver()

    ! initialize matrices and vectors
    call init_data()

    ! set up M loss matrix
    call build_loss_matrix(loss) 

    ! set up F production matrix
    call build_prod_matrix(prod)

    ! set up krylov info
    call KSPSetOperators(krylov, loss%M, loss%M, SAME_NONZERO_PATTERN, ierr)
    call KSPSetUp(krylov,ierr)

    ! calculate preconditioner (ILU)
    call PCFactorGetMatrix(prec,loss%M,ierr)

    ! begin power iteration 
    call execute_power_iter()

    ! extract results
    call extract_results()

    ! deallocate petsc objects
    call finalize()

  end subroutine cmfd_power_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

    integer :: n      ! problem size
    real(8) :: guess  ! initial guess

    ! set up matrices
    call init_M_operator(loss)
    call init_F_operator(prod)

    ! get problem size
    n = loss%n

    ! set up flux vectors
    call VecCreate(PETSC_COMM_SELF,phi_n,ierr)
    call VecSetSizes(phi_n,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(phi_n,ierr)
    call VecCreate(PETSC_COMM_SELF,phi_o,ierr)
    call VecSetSizes(phi_o,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(phi_o,ierr)

    ! set up source vectors
    call VecCreate(PETSC_COMM_SELF,S_n,ierr)
    call VecSetSizes(S_n,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(S_n,ierr)
    call VecCreate(PETSC_COMM_SELF,S_o,ierr)
    call VecSetSizes(S_o,PETSC_DECIDE,n,ierr)
    call VecSetFromOptions(S_o,ierr)

    ! set initial guess
    guess = 1.0_8 
    call VecSet(phi_n,guess,ierr)
    call VecSet(phi_o,guess,ierr)
    k_n = guess
    k_o = guess

  end subroutine init_data

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver()

    real(8)     :: solvertol ! krylov tolerance

    ! set tolerance
    solvertol = 1.0e-7_8

    ! set up krylov solver
    call KSPCreate(PETSC_COMM_SELF,krylov,ierr)
    call KSPSetTolerances(krylov,solvertol,PETSC_DEFAULT_DOUBLE_PRECISION,     &
   &                      PETSC_DEFAULT_DOUBLE_PRECISION,                      &
   &                      PETSC_DEFAULT_INTEGER,ierr)
    call KSPSetType(krylov,KSPGMRES,ierr)
    call KSPSetInitialGuessNonzero(krylov,PETSC_TRUE,ierr)
    call KSPSetInitialGuessNonzero(krylov,PETSC_TRUE,ierr)
    call KSPGetPC(krylov,prec,ierr)
    call PCSetType(prec,PCILU,ierr)
    call PCFactorSetLevels(prec,5,ierr)
    call KSPSetFromOptions(krylov,ierr)

  end subroutine init_solver

!===============================================================================
! EXECUTE_POWER_ITER  in the main power iteration routine 
!                     for the cmfd calculation
!===============================================================================

  subroutine execute_power_iter()

    real(8)     :: num       ! numerator for eigenvalue update
    real(8)     :: den       ! denominator for eigenvalue update
    real(8)     :: one=1.0_8 ! one
    integer     :: i         ! iteration counter

    ! reset convergence flag
    iconv = .FALSE.

    ! begin power iteration
    do i = 1,10000

      ! compute source vector
      call MatMult(prod%F,phi_o,S_o,ierr)

      ! normalize source vector
      call VecScale(S_o,one/k_o,ierr)

      ! compute new flux vector
      call KSPSolve(krylov,S_o,phi_n,ierr)

      ! compute new source vector
      call MatMult(prod%F,phi_n,S_n,ierr)

      ! compute new k-eigenvalue
      call VecSum(S_n,num,ierr)
      call VecSum(S_o,den,ierr)
      k_n = num/den

      ! renormalize the old source
      call VecScale(S_o,k_o,ierr)

      ! check convergence
      call convergence()

      ! to break or not to break
      if (iconv) exit

      ! record old values
      call VecCopy(phi_n,phi_o,ierr)
      k_o = k_n

    end do

  end subroutine execute_power_iter 

!===============================================================================
! CONVERGENCE checks the convergence of eigenvalue, eigenvector and source
!===============================================================================

  subroutine convergence()

    real(8)     :: ktol = 1.e-6_8 ! tolerance on keff
    real(8)     :: stol = 1.e-5_8 ! tolerance on source
    real(8)     :: kerr           ! error in keff
    real(8)     :: serr           ! error in source
    real(8)     :: one = -1.0_8   ! one
    real(8)     :: norm_n         ! L2 norm of new source
    real(8)     :: norm_o         ! L2 norm of old source
    integer     :: floc           ! location of max error in flux
    integer     :: sloc           ! location of max error in source
    integer     :: ierr           ! petsc error code
    integer     :: n              ! vector size

    ! reset convergence flag
    iconv = .FALSE.

    ! calculate error in keff
    kerr = abs(k_o - k_n)/k_n

    ! calculate max error in source
    call VecNorm(S_n,NORM_2,norm_n,ierr)
    call VecNorm(S_o,NORM_2,norm_o,ierr)
    serr = abs(norm_n-norm_o)/norm_n

    ! check for convergence
    if(kerr < ktol .and. serr < stol) iconv = .TRUE.

    ! print out to user (TODO: make formatted)
    print *,k_n,kerr,serr

  end subroutine convergence

!==============================================================================
! EXTRACT_RESULTS
!==============================================================================

  subroutine extract_results()

    use global, only: cmfd,path_input

    PetscViewer :: viewer
    PetscScalar, pointer :: phi_v(:) ! pointer to eigenvector info
    integer :: n ! problem size

    ! get problem size
    n = loss%n

    ! also allocate in cmfd object
    if (.not. allocated(cmfd%phi)) allocate(cmfd%phi(n))

    ! convert petsc phi_object to cmfd_obj
    call VecGetArrayF90(phi_n,phi_v,ierr)
    cmfd%phi = phi_v

    call VecRestoreArrayF90(phi_n,phi_v,ierr)

    ! save eigenvalue
    cmfd%keff = k_n

    ! write out results
    call PetscViewerBinaryOpen(PETSC_COMM_SELF,trim(path_input)//'fluxvec.bin' &
   &     ,FILE_MODE_WRITE,viewer,ierr)
    call VecView(phi_n,viewer,ierr)
    call PetscViewerDestroy(viewer,ierr)

  end subroutine extract_results

!==============================================================================
! FINALIZE
!==============================================================================

  subroutine finalize()

    ! finalize data objects
    call destroy_M_operator(loss)
    call destroy_F_operator(prod)

    call VecDestroy(phi_n,ierr)
    call VecDestroy(phi_o,ierr)
    call VecDestroy(S_n,ierr)
    call VecDestroy(S_o,ierr)

    ! finalize solver objects
!   call KSPDestroy(krylov,ierr)

  end subroutine finalize

#endif

end module cmfd_power_solver
