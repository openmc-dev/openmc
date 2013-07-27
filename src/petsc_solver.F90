module petsc_solver

  use matrix_header,  only: Matrix
  use vector_header,  only: Vector

  implicit none
  private

#ifdef PETSC
#  include <finclude/petsc.h90>
#endif

  ! PETSc GMRES solver type 
  type, public :: Petsc_gmres 
#ifdef PETSC
    KSP :: ksp
    PC  :: pc
#endif
   contains
     procedure :: create       => gmres_create
     procedure :: set_oper     => gmres_set_oper
     procedure :: destroy      => gmres_destroy
     procedure :: solve        => gmres_solve
  end type Petsc_gmres

  ! Derived type to contain list of data needed to be passed to procedures
  type, public :: Jfnk_ctx
    procedure (res_interface), pointer, nopass :: res_proc_ptr
    procedure (jac_interface), pointer, nopass :: jac_proc_ptr
  end type Jfnk_ctx

  ! Petsc SNES JFNK solver type 
  type, public :: Petsc_jfnk
#ifdef PETSC
    KSP  :: ksp
    PC   :: pc
    SNES :: snes
    Mat  :: jac_mf
    SNESLineSearch :: ls
#endif
   contains
     procedure :: create        => jfnk_create
     procedure :: destroy       => jfnk_destroy
     procedure :: set_functions => jfnk_set_functions
     procedure :: solve         => jfnk_solve
  end type Petsc_jfnk

  integer :: petsc_err

  ! Abstract interface stating how jacobian and residual routines look
  abstract interface
    subroutine res_interface(x, res)
      import :: Vector
      type(Vector) :: x
      type(Vector) :: res
    end subroutine res_interface

    subroutine jac_interface(x)
      import :: Vector
      type(Vector) :: x
    end subroutine jac_interface
  end interface

contains

!===============================================================================
! GMRES_CREATE sets up a PETSc GMRES solver
!===============================================================================

  subroutine gmres_create(self)

    class(Petsc_gmres) :: self

    integer :: ilu_levels = 5
    real(8) :: rtol = 1.0e-10_8
    real(8) :: atol = 1.0e-10_8

#ifdef PETSC
    call KSPCreate(PETSC_COMM_WORLD, self % ksp, petsc_err)
    call KSPSetTolerances(self % ksp, rtol, atol, &
         PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER, petsc_err)
    call KSPSetType(self % ksp, KSPGMRES, petsc_err)
    call KSPSetInitialGuessNonzero(self % ksp, PETSC_TRUE, petsc_err)
    call KSPGetPC(self % ksp, self % pc, petsc_err)
    call PCFactorSetLevels(self % pc, ilu_levels, petsc_err)
    call KSPSetFromOptions(self % ksp, petsc_err)
#endif

  end subroutine gmres_create

!===============================================================================
! GMRES_SET_OPER sets the matrix opetors for the GMRES solver
!===============================================================================

  subroutine gmres_set_oper(self, prec_mat, mat_in)

    class(Petsc_gmres) :: self
    type(Matrix)       :: prec_mat
    type(Matrix)       :: mat_in

#ifdef PETSC
    call KSPSetOperators(self % ksp, mat_in % petsc_mat, prec_mat % petsc_mat, &
         SAME_NONZERO_PATTERN, petsc_err)
    call KSPSetUp(self % ksp, petsc_err) 
#endif

  end subroutine gmres_set_oper

!===============================================================================
! GMRES_DESTROY frees all memory associated with the GMRES solver
!===============================================================================

  subroutine gmres_destroy(self)

   class(Petsc_gmres) :: self

#ifdef PETSC
    call KSPDestroy(self % ksp, petsc_err)
#endif

  end subroutine gmres_destroy

!===============================================================================
! GMRES_SOLVE solves the linear system
!===============================================================================

  subroutine gmres_solve(self, b, x)

    class(Petsc_gmres) :: self
    type(Vector)       :: b
    type(Vector)       :: x

#ifdef PETSC
    call KSPSolve(self % ksp, b % petsc_vec, x % petsc_vec, petsc_err)
#endif

  end subroutine gmres_solve

!===============================================================================
! JFNK_CREATE sets up a JFNK solver using PETSc SNES
!===============================================================================

  subroutine jfnk_create(self)

    class(Petsc_jfnk) :: self

#ifdef PETSC
    ! Turn on mf_operator option for matrix free jacobian
    call PetscOptionsSetValue("-snes_mf_operator", "TRUE", petsc_err)

    ! Create the SNES context
    call SNESCreate(PETSC_COMM_WORLD, self % snes, petsc_err)

    ! Set up the GMRES solver
    call SNESGetKSP(self % snes, self % ksp, petsc_err)
    call KSPSetType(self % ksp, KSPGMRES, petsc_err)    

    ! Apply options
    call SNESGetLineSearch(self % snes, self % ls, petsc_err)
    call SNESLineSearchSetType(self % ls, SNESLINESEARCHBASIC, petsc_err)
    call SNESSetFromOptions(self % snes, petsc_err)

    ! Set up preconditioner
    call KSPGetPC(self % ksp, self % pc, petsc_err)
    call PCSetType(self % pc, PCILU, petsc_err)
    call PCFactorSetLevels(self % pc, 5, petsc_err)
    call PCSetFromOptions(self % pc, petsc_err)
    call KSPSetFromOptions(self % ksp, petsc_err)
#endif

  end subroutine jfnk_create

!===============================================================================
! JFNK_DESTROY frees all memory associated with JFNK solver
!===============================================================================

  subroutine jfnk_destroy(self)

    class(Petsc_jfnk) :: self

#ifdef PETSC
    call MatDestroy(self % jac_mf, petsc_err)
    call SNESDestroy(self % snes, petsc_err)
#endif

  end subroutine jfnk_destroy

!===============================================================================
! JFNK_SET_FUNCTIONS sets user functions and matrix free objects
!===============================================================================

  subroutine jfnk_set_functions(self, ctx, res, jac_prec)

    class(Petsc_jfnk) :: self
    type(Jfnk_ctx)    :: ctx
    type(Vector)      :: res
    type(Matrix)      :: jac_prec

# ifdef PETSC
    ! Set residual procedure
    call SNESSetFunction(self % snes, res % petsc_vec, jfnk_compute_residual, &
         ctx, petsc_err)

    ! Create the matrix free jacobian
    call MatCreateSNESMF(self % snes, self % jac_mf, petsc_err)

    ! Set Jacobian procedure
    call SNESSetJacobian(self % snes, self % jac_mf, jac_prec % petsc_mat, &
         jfnk_compute_jacobian, ctx, petsc_err)

    ! Set up Jacobian Lags
    call SNESSetLagJacobian(self % snes, -2, petsc_err)
    call SNESSetLagPreconditioner(self % snes, -1, petsc_err)
#endif

  end subroutine jfnk_set_functions

!===============================================================================
! JFNK_SOLVE solves the nonlinear system
!===============================================================================

  subroutine jfnk_solve(self, xvec)

    class(Petsc_jfnk) :: self
    type(Vector)      :: xvec

#ifdef PETSC
    call SNESSolve(self % snes, PETSC_NULL_DOUBLE, xvec % petsc_vec, petsc_err)
#endif

  end subroutine jfnk_solve

!===============================================================================
! JFNK_COMPUTE_RESIDUAL buffer routine to user specifed residual routine
!===============================================================================
# ifdef PETSC
  subroutine jfnk_compute_residual(snes, x, res, ctx, ierr)

    SNES           :: snes
    Vec            :: x
    Vec            :: res
    integer        :: ierr
    type(Jfnk_ctx) :: ctx

    type(Vector) :: xvec
    type(Vector) :: resvec

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

  end subroutine jfnk_compute_residual
#endif
!===============================================================================
! JFNK_COMPUTE_JACOBIAN
!===============================================================================
#ifdef PETSC
  subroutine jfnk_compute_jacobian(snes, x, jac_mf, jac_prec, flag, ctx, &
             ierr)

    SNES           :: snes
    Vec            :: x
    Mat            :: jac_mf
    Mat            :: jac_prec
    MatStructure   :: flag 
    type(Jfnk_ctx) :: ctx
    integer        :: ierr

    type(Vector) :: xvec

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

  end subroutine jfnk_compute_jacobian
#endif
end module petsc_solver
