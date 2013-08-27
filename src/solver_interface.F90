module solver_interface

  use error,          only: fatal_error
  use global,         only: message
  use matrix_header,  only: Matrix
  use vector_header,  only: Vector

#ifndef PETSC
  use gmres_solver,   only: ilut, pgmres
#endif

  implicit none
  private

#ifdef PETSC
#  include <finclude/petsc.h90>
#endif

  ! GMRES solver type 
  type, public :: GMRESSolver 
#ifdef PETSC
    KSP :: ksp
    PC  :: pc
#else
    integer :: im      ! iterations before restart
    integer :: iout    ! unit number for printing gmres iters
    integer :: maxits  ! maximum iterations
    integer :: lfil    ! fill in parameter
    integer :: iwk     ! size of MSR storage
    integer, allocatable :: ju(:) ! row pointers for ILUT result
    integer, allocatable :: jlu(:) ! MSR storage of j
    integer, allocatable :: jr(:) ! temp work array
    integer, allocatable :: jwu(:) ! temp work array
    integer, allocatable :: jwl(:) ! temp work array
    real(8) :: tol     ! gmres tolerance
    real(8) :: droptol ! ILU drop tolerance
    real(8), allocatable :: alu(:) ! MSR storage of ILU A
    real(8), allocatable :: vv(:) ! GMRES work array
    real(8), allocatable :: wu(:) ! temp work array
    real(8), allocatable :: wl(:) ! temp work array
    real(8), allocatable :: rhs(:) ! copy of rhs array
    type(Matrix), pointer :: A ! coefficient matrix
#endif
   contains
#ifdef PETSC
     procedure :: create       => petsc_gmres_create
     procedure :: set_oper     => petsc_gmres_set_oper
     procedure :: destroy      => petsc_gmres_destroy
     procedure :: solve        => petsc_gmres_solve
#else
     procedure :: create       => gmres_create
     procedure :: set_oper     => gmres_set_oper
     procedure :: destroy      => gmres_destroy
     procedure :: solve        => gmres_solve
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
    KSP  :: ksp
    PC   :: pc
    SNES :: snes
    Mat  :: jac_mf
    SNESLineSearch :: ls
#endif
   contains
#ifdef PETSC
     procedure :: create        => petsc_jfnk_create
     procedure :: destroy       => petsc_jfnk_destroy
     procedure :: set_functions => petsc_jfnk_set_functions
     procedure :: solve         => petsc_jfnk_solve
#else
     procedure :: create        => jfnk_create
     procedure :: destroy       => jfnk_destroy
     procedure :: set_functions => jfnk_set_functions
     procedure :: solve         => jfnk_solve
#endif
  end type JFNKSolver

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

#ifdef PETSC
  integer :: petsc_err ! petsc error code
#endif

contains

#ifndef PETSC
!===============================================================================
! GMRES_CREATE sets up a PETSc GMRES solver
!===============================================================================

  subroutine gmres_create(self)

    class(GMRESSolver) :: self

    ! Set up default values
    self % im = 40
    self % iout = 0
    self % maxits = 1000
    self % lfil = 0
    self % tol = 1.e-10_8
    self % droptol = 0.0_8

  end subroutine gmres_create

!===============================================================================
! GMRES_SET_OPER sets the matrix opetors for the GMRES solver
!===============================================================================

  subroutine gmres_set_oper(self, prec_mat, mat_in)

    class(GMRESSolver) :: self
    type(Matrix)       :: prec_mat
    type(Matrix), target :: mat_in

    integer :: ierr

    ! Set pointer to matrix
    self % A => mat_in

    ! Set up data
    self % iwk = self % A % nnz + 1

    ! Set up arrays
    allocate(self % ju(self % A % n))
    allocate(self % jlu(self % iwk))
    allocate(self % alu(self % iwk))
    allocate(self % jr(self % A % n))
    allocate(self % jwu(self % A % n))
    allocate(self % jwl(self % A % n))
    allocate(self % vv(self % A % n * (self % im + 1)))
    allocate(self % wu(self % A % n + 1))
    allocate(self % wl(self % A % n))
    allocate(self % rhs(self % A % n))

    ! Perform ILU Preconditioning
    call ilut(self % A % n, self % A % val, self % A % col, self % A % row, &
              self % lfil, self % droptol, self % alu, self % jlu, self % ju, &
              self % iwk, self % wu, self % wl, self % jr, self % jwl, &
              self % jwu, ierr)

  end subroutine gmres_set_oper

!===============================================================================
! GMRES_DESTROY frees all memory associated with the GMRES solver
!===============================================================================

  subroutine gmres_destroy(self)

   class(GMRESSolver) :: self

   ! Destroy all arrays associated with GMRES solver
   deallocate(self % ju)
   deallocate(self % jlu)
   deallocate(self % alu)
   deallocate(self % jr)
   deallocate(self % jwu)
   deallocate(self % jwl)
   deallocate(self % vv)
   deallocate(self % wu)
   deallocate(self % wl)
   deallocate(self % rhs)

  end subroutine gmres_destroy

!===============================================================================
! GMRES_SOLVE solves the linear system
!===============================================================================

  subroutine gmres_solve(self, b, x)

    class(GMRESSolver) :: self
    type(Vector)       :: b
    type(Vector)       :: x

    integer :: ierr ! error code

    ! Copy rhs
    self % rhs = b % val

    ! Perform GMRES calculation
    call pgmres(self % A % n, self % im, self % rhs, x % val, self % vv, &
                self % tol, self % maxits, self % iout, self % A % val, &
                self % A % col, self % A % row, self % alu, self % jlu, &
                self % ju, ierr)

  end subroutine gmres_solve

!===============================================================================
! JFNK_CREATE sets up a JFNK solver using PETSc SNES
!===============================================================================

  subroutine jfnk_create(self)

    class(JFNKSolver) :: self

    message = 'JFNK cannot be used without PETSc.'
    call fatal_error()

  end subroutine jfnk_create

!===============================================================================
! JFNK_DESTROY frees all memory associated with JFNK solver
!===============================================================================

  subroutine jfnk_destroy(self)

    class(JFNKSolver) :: self

    message = 'JFNK cannot be used without PETSc.'
    call fatal_error()

  end subroutine jfnk_destroy

!===============================================================================
! JFNK_SET_FUNCTIONS sets user functions and matrix free objects
!===============================================================================

  subroutine jfnk_set_functions(self, ctx, res, jac_prec)

    class(JFNKSolver) :: self
    type(Jfnk_ctx)    :: ctx
    type(Vector)      :: res
    type(Matrix)      :: jac_prec

    message = 'JFNK cannot be used without PETSc.'
    call fatal_error()

  end subroutine jfnk_set_functions

!===============================================================================
! JFNK_SOLVE solves the nonlinear system
!===============================================================================

  subroutine jfnk_solve(self, xvec)

    class(JFNKSolver) :: self
    type(Vector)      :: xvec

    message = 'JFNK cannot be used without PETSc.'
    call fatal_error()

  end subroutine jfnk_solve

#else
!===============================================================================
! PETSC_GMRES_CREATE sets up a PETSc GMRES solver
!===============================================================================

  subroutine petsc_gmres_create(self)

    class(GMRESSolver) :: self

    integer :: ilu_levels = 5
    real(8) :: rtol = 1.0e-10_8
    real(8) :: atol = 1.0e-10_8

    call KSPCreate(PETSC_COMM_WORLD, self % ksp, petsc_err)
    call KSPSetTolerances(self % ksp, rtol, atol, &
         PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER, petsc_err)
    call KSPSetType(self % ksp, KSPGMRES, petsc_err)
    call KSPSetInitialGuessNonzero(self % ksp, PETSC_TRUE, petsc_err)
    call KSPGetPC(self % ksp, self % pc, petsc_err)
    call PCFactorSetLevels(self % pc, ilu_levels, petsc_err)
    call KSPSetFromOptions(self % ksp, petsc_err)

  end subroutine petsc_gmres_create

!===============================================================================
! PETSC_GMRES_SET_OPER sets the matrix opetors for the GMRES solver
!===============================================================================

  subroutine petsc_gmres_set_oper(self, prec_mat, mat_in)

    class(GMRESSolver) :: self
    type(Matrix)       :: prec_mat
    type(Matrix)       :: mat_in

    call KSPSetOperators(self % ksp, mat_in % petsc_mat, prec_mat % petsc_mat, &
         SAME_NONZERO_PATTERN, petsc_err)
    call KSPSetUp(self % ksp, petsc_err) 

  end subroutine petsc_gmres_set_oper

!===============================================================================
! PETSC_GMRES_DESTROY frees all memory associated with the GMRES solver
!===============================================================================

  subroutine petsc_gmres_destroy(self)

   class(GMRESSolver) :: self

    call KSPDestroy(self % ksp, petsc_err)

  end subroutine petsc_gmres_destroy

!===============================================================================
! PETSC_GMRES_SOLVE solves the linear system
!===============================================================================

  subroutine petsc_gmres_solve(self, b, x)

    class(GMRESSolver) :: self
    type(Vector)       :: b
    type(Vector)       :: x

    call KSPSolve(self % ksp, b % petsc_vec, x % petsc_vec, petsc_err)

  end subroutine petsc_gmres_solve

!===============================================================================
! PETSC_JFNK_CREATE sets up a JFNK solver using PETSc SNES
!===============================================================================

  subroutine petsc_jfnk_create(self)

    class(JFNKSolver) :: self

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

  end subroutine petsc_jfnk_create

!===============================================================================
! PETSC_JFNK_DESTROY frees all memory associated with JFNK solver
!===============================================================================

  subroutine petsc_jfnk_destroy(self)

    class(JFNKSolver) :: self

    call MatDestroy(self % jac_mf, petsc_err)
    call SNESDestroy(self % snes, petsc_err)

  end subroutine petsc_jfnk_destroy

!===============================================================================
! PETSC_JFNK_SET_FUNCTIONS sets user functions and matrix free objects
!===============================================================================

  subroutine petsc_jfnk_set_functions(self, ctx, res, jac_prec)

    class(JFNKSolver) :: self
    type(Jfnk_ctx)    :: ctx
    type(Vector)      :: res
    type(Matrix)      :: jac_prec

    ! Set residual procedure
    call SNESSetFunction(self % snes, res % petsc_vec, &
         petsc_jfnk_compute_residual, ctx, petsc_err)

    ! Create the matrix free jacobian
    call MatCreateSNESMF(self % snes, self % jac_mf, petsc_err)

    ! Set Jacobian procedure
    call SNESSetJacobian(self % snes, self % jac_mf, jac_prec % petsc_mat, &
         petsc_jfnk_compute_jacobian, ctx, petsc_err)

    ! Set up Jacobian Lags
    call SNESSetLagJacobian(self % snes, -2, petsc_err)
    call SNESSetLagPreconditioner(self % snes, -1, petsc_err)

  end subroutine petsc_jfnk_set_functions

!===============================================================================
! PETSC_JFNK_SOLVE solves the nonlinear system
!===============================================================================

  subroutine petsc_jfnk_solve(self, xvec)

    class(JFNKSolver) :: self
    type(Vector)      :: xvec

    call SNESSolve(self % snes, PETSC_NULL_DOUBLE, xvec % petsc_vec, petsc_err)

  end subroutine petsc_jfnk_solve

!===============================================================================
! PETSC_JFNK_COMPUTE_RESIDUAL buffer routine to user specifed residual routine
!===============================================================================

  subroutine petsc_jfnk_compute_residual(snes, x, res, ctx, ierr)

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

  end subroutine petsc_jfnk_compute_residual

!===============================================================================
! PETSC_JFNK_COMPUTE_JACOBIAN buffer routine to user specified jacobian routine
!===============================================================================

  subroutine petsc_jfnk_compute_jacobian(snes, x, jac_mf, jac_prec, flag, ctx, &
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

  end subroutine petsc_jfnk_compute_jacobian
#endif

end module solver_interface
