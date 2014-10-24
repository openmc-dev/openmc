module solver_interface

  use error,          only: fatal_error
  use matrix_header,  only: Matrix
  use vector_header,  only: Vector

#ifdef PETSC
  use petscksp
  use petscsnes
#endif

  implicit none
  private

  ! GMRES solver type 
  type, public :: GMRESSolver 
#ifdef PETSC
    type(ksp) :: ksp_ ! Krylov linear solver instance
    type(pc)  :: pc_  ! Preconditioner instance
#endif
   contains
#ifdef PETSC
     procedure :: create       => petsc_gmres_create
     procedure :: set_oper     => petsc_gmres_set_oper
     procedure :: destroy      => petsc_gmres_destroy
     procedure :: solve        => petsc_gmres_solve
#endif
  end type GMRESSolver

  ! Derived type to contain list of data needed to be passed to procedures
  type, public :: Jfnk_ctx
    procedure (res_interface), pointer, nopass :: res_proc_ptr
    procedure (jac_interface), pointer, nopass :: jac_proc_ptr
  end type Jfnk_ctx

  ! JFNK solver type 
  type, public :: JFNKSolver 
#ifdef PETSC
    type(ksp)  :: ksp_   ! Krylov linear solver instance
    type(pc)   :: pc_    ! Preconditioner instance
    type(snes) :: snes_  ! Nonlinear solver instance
    type(mat)  :: jac_mf ! Matrix free jacobian instance
    integer    :: ls     ! Line search instance
#endif
   contains
#ifdef PETSC
     procedure :: create        => petsc_jfnk_create
     procedure :: destroy       => petsc_jfnk_destroy
     procedure :: set_functions => petsc_jfnk_set_functions
     procedure :: solve         => petsc_jfnk_solve
#endif
  end type JFNKSolver

  ! Abstract interface stating how jacobian and residual routines look
  abstract interface
    subroutine res_interface(x, res)
      import :: Vector
      type(Vector), intent(in)    :: x   ! solution vector
      type(Vector), intent(inout) :: res ! residual vector
    end subroutine res_interface

    subroutine jac_interface(x)
      import :: Vector
      type(Vector), intent(in) :: x ! solution vector
    end subroutine jac_interface
  end interface

#ifdef PETSC
  integer :: petsc_err ! petsc error code
#endif

contains

#ifdef PETSC
!===============================================================================
! PETSC_GMRES_CREATE sets up a PETSc GMRES solver
!===============================================================================

  subroutine petsc_gmres_create(self)

    class(GMRESSolver), intent(inout) :: self ! GMRES solver instance

    integer :: ilu_levels = 5
    real(8) :: rtol = 1.0e-10_8
    real(8) :: atol = 1.0e-10_8

    call KSPCreate(PETSC_COMM_WORLD, self % ksp_, petsc_err)
    call KSPSetTolerances(self % ksp_, rtol, atol, &
         PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER, petsc_err)
    call KSPSetType(self % ksp_, 'gmres', petsc_err)
    call KSPSetInitialGuessNonzero(self % ksp_, PETSC_TRUE, petsc_err)
    call KSPGetPC(self % ksp_, self % pc_, petsc_err)
    call PCFactorSetLevels(self % pc_, ilu_levels, petsc_err)
    call KSPSetFromOptions(self % ksp_, petsc_err)

  end subroutine petsc_gmres_create

!===============================================================================
! PETSC_GMRES_SET_OPER sets the matrix opetors for the GMRES solver
!===============================================================================

  subroutine petsc_gmres_set_oper(self, prec_mat, mat_in)

    class(GMRESSolver), intent(inout) :: self     ! GMRES solver instanace
    type(Matrix), intent(inout)       :: prec_mat ! preconditioner matrix
    type(Matrix), intent(inout)       :: mat_in   ! coefficient matrix

    call KSPSetOperators(self % ksp_, mat_in % petsc_mat, prec_mat % petsc_mat, &
         SAME_NONZERO_PATTERN, petsc_err)
    call KSPSetUp(self % ksp_, petsc_err) 

  end subroutine petsc_gmres_set_oper

!===============================================================================
! PETSC_GMRES_DESTROY frees all memory associated with the GMRES solver
!===============================================================================

  subroutine petsc_gmres_destroy(self)

   class(GMRESSolver), intent(inout) :: self ! GMRES solver instance

    call KSPDestroy(self % ksp_, petsc_err)

  end subroutine petsc_gmres_destroy

!===============================================================================
! PETSC_GMRES_SOLVE solves the linear system
!===============================================================================

  subroutine petsc_gmres_solve(self, b, x)

    class(GMRESSolver), intent(inout) :: self ! GMRES solver instance
    type(Vector), intent(inout)       :: b    ! right hand side vector
    type(Vector), intent(inout)       :: x    ! solution vector

    call KSPSolve(self % ksp_, b % petsc_vec, x % petsc_vec, petsc_err)

  end subroutine petsc_gmres_solve

!===============================================================================
! PETSC_JFNK_CREATE sets up a JFNK solver using PETSc SNES
!===============================================================================

  subroutine petsc_jfnk_create(self)

    class(JFNKSolver), intent(inout) :: self ! JFNK solver instance

    ! Turn on mf_operator option for matrix free jacobian
    call PetscOptionsSetValue("-snes_mf_operator", "TRUE", petsc_err)

    ! Create the SNES context
    call SNESCreate(PETSC_COMM_WORLD, self % snes_, petsc_err)

    ! Set up the GMRES solver
    call SNESGetKSP(self % snes_, self % ksp_, petsc_err)
    call KSPSetType(self % ksp_, 'gmres', petsc_err)    

    ! Apply options
    call SNESGetLineSearch(self % snes_, self % ls, petsc_err)
    call SNESLineSearchSetType(self % ls, 'basic', petsc_err)
    call SNESSetFromOptions(self % snes_, petsc_err)

    ! Set up preconditioner
    call KSPGetPC(self % ksp_, self % pc_, petsc_err)
    call PCSetType(self % pc_, 'ilu', petsc_err)
    call PCFactorSetLevels(self % pc_, 5, petsc_err)
    call PCSetFromOptions(self % pc_, petsc_err)
    call KSPSetFromOptions(self % ksp_, petsc_err)

  end subroutine petsc_jfnk_create

!===============================================================================
! PETSC_JFNK_DESTROY frees all memory associated with JFNK solver
!===============================================================================

  subroutine petsc_jfnk_destroy(self)

    class(JFNKSolver), intent(inout) :: self ! JFNK solver instance

    call MatDestroy(self % jac_mf, petsc_err)
    call SNESDestroy(self % snes_, petsc_err)

  end subroutine petsc_jfnk_destroy

!===============================================================================
! PETSC_JFNK_SET_FUNCTIONS sets user functions and matrix free objects
!===============================================================================

  subroutine petsc_jfnk_set_functions(self, ctx, res, jac_prec)

    class(JFNKSolver), intent(inout) :: self     ! JFNK solver instance
    type(Jfnk_ctx),    intent(inout) :: ctx      ! JFNK context instance
    type(Vector),      intent(inout) :: res      ! residual vector
    type(Matrix),      intent(inout) :: jac_prec ! preconditioner matrix

    ! Set residual procedure
    call SNESSetFunction(self % snes_, res % petsc_vec, &
         petsc_jfnk_compute_residual, ctx, petsc_err)

    ! Create the matrix free jacobian
    call MatCreateSNESMF(self % snes_, self % jac_mf, petsc_err)

    ! Set Jacobian procedure
    call SNESSetJacobian(self % snes_, self % jac_mf, jac_prec % petsc_mat, &
         petsc_jfnk_compute_jacobian, ctx, petsc_err)

    ! Set up Jacobian Lags
    call SNESSetLagJacobian(self % snes_, -2, petsc_err)
    call SNESSetLagPreconditioner(self % snes_, -1, petsc_err)

  end subroutine petsc_jfnk_set_functions

!===============================================================================
! PETSC_JFNK_SOLVE solves the nonlinear system
!===============================================================================

  subroutine petsc_jfnk_solve(self, xvec)

    class(JFNKSolver), intent(inout) :: self ! JFNK instance
    type(Vector), intent(inout)      :: xvec ! solution vector

    call SNESSolve(self % snes_, PETSC_NULL_DOUBLE, xvec % petsc_vec, petsc_err)

  end subroutine petsc_jfnk_solve

!===============================================================================
! PETSC_JFNK_COMPUTE_RESIDUAL buffer routine to user specifed residual routine
!===============================================================================

  subroutine petsc_jfnk_compute_residual(snes_, x, res, ctx, ierr)

    type(snes), intent(inout)     :: snes_ ! PETSc SNES object
    type(vec), intent(inout)      :: x     ! PETSc solution vector
    type(vec), intent(inout)      :: res   ! PETSc residual vector
    integer, intent(inout)        :: ierr  ! error code
    type(Jfnk_ctx), intent(inout) :: ctx ! JFNK context instance

    type(Vector) :: xvec   ! solution vector
    type(Vector) :: resvec ! residual vector

    ! We need to use the x and res that come from PETSc because the pointer
    ! location changes and therefore we can not use module variables
    ! for the residual and x vectors here

    ! Need to point an OpenMC vector to PETSc vector
    call VecGetArrayF90(x, xvec % val, ierr)
    call VecGetArrayF90(res, resvec % val, ierr)

    ! Call user residual routine
    call ctx % res_proc_ptr(xvec, resvec)

    ! Need to restore the PETSc vector
    call VecRestoreArrayF90(x, xvec % val, ierr)
    call VecRestoreArrayF90(res, resvec % val, ierr)

  end subroutine petsc_jfnk_compute_residual

!===============================================================================
! PETSC_JFNK_COMPUTE_JACOBIAN buffer routine to user specified jacobian routine
!===============================================================================

  subroutine petsc_jfnk_compute_jacobian(snes_, x, jac_mf, jac_prec, flag, &
             ctx, ierr)

    type(snes), intent(inout)     :: snes_     ! PETSc snes instance
    type(vec), intent(inout)      :: x         ! PETSc solution vector
    type(mat), intent(inout)      :: jac_mf    ! PETSc matrix free jacobian
    type(mat), intent(inout)      :: jac_prec  ! PETSc matrix jacobian precond.
    integer, intent(inout)        :: flag      ! unused madatory flag
    type(Jfnk_ctx), intent(inout) :: ctx       ! JFNK context instance
    integer, intent(inout)        :: ierr      ! error code

    type(Vector) :: xvec ! solution vector

    ! Again, we use the vector that comes from Petsc to build the Jacobian
    ! matrix

    ! Need to point OpenMC vector to PETSc Vector
    call VecGetArrayF90(x, xvec % val, ierr)

    ! Evaluate user jacobian routine
    call ctx % jac_proc_ptr(xvec)

    ! Restore the PETSc vector
    call VecRestoreArrayF90(x, xvec % val, ierr)

    call MatAssemblyBegin(jac_mf, MAT_FINAL_ASSEMBLY, petsc_err)
    call MatAssemblyEnd(jac_mf, MAT_FINAL_ASSEMBLY, petsc_err)

  end subroutine petsc_jfnk_compute_jacobian
#endif

end module solver_interface
