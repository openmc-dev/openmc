module petsc_solver

  use matrix_header,  only: Matrix
  use vector_header,  only: Vector

  implicit none
  private

#ifdef PETSC
#  include <finclude/petsc.h90>
#endif

  type, public :: Petsc_gmres 
#ifdef PETSC
    KSP :: ksp
    PC  :: pc
#endif
   contains
     procedure :: create       => gmres_create
     procedure :: precondition => gmres_precondition
     procedure :: set_oper     => gmres_set_oper
     procedure :: destroy      => gmres_destroy
     procedure :: solve        => gmres_solve
  end type Petsc_gmres

  integer :: Petsc_err

contains

!===============================================================================
! GMRES_CREATE
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
#endif

  end subroutine gmres_create

!===============================================================================
! GMRES_PRECONDITION
!===============================================================================

  subroutine gmres_precondition(self, mat)

    class(Petsc_gmres) :: self
    type(Matrix)       :: mat

#ifdef PETSC
    call KSPSetUp(self % ksp, petsc_err) 
    call PCFactorGetMatrix(self % pc, mat % petsc_mat, petsc_err) 
#endif

  end subroutine gmres_precondition

!===============================================================================
! GMRES_SET_OPER
!===============================================================================

  subroutine gmres_set_oper(self, prec_mat, mat)

    class(Petsc_gmres) :: self
    type(Matrix)       :: prec_mat
    type(Matrix)       :: mat

#ifdef PETSC
    call KSPSetOperators(self % ksp, mat % petsc_mat, prec_mat % petsc_mat, &
         SAME_NONZERO_PATTERN, petsc_err)
#endif

  end subroutine gmres_set_oper

!===============================================================================
! GMRES_DESTROY
!===============================================================================

  subroutine gmres_destroy(self)

   class(Petsc_gmres) :: self

#ifdef PETSC
    call KSPDestroy(self % ksp, petsc_err)
#endif

  end subroutine gmres_destroy

!===============================================================================
! GMRES_SOLVE
!===============================================================================

  subroutine gmres_solve(self, b, x)

    class(Petsc_gmres) :: self
    type(Vector)       :: b
    type(Vector)       :: x

    integer :: petsc_err

#ifdef PETSC
    call KSPSolve(self % ksp, b % petsc_vec, x % petsc_vec, petsc_err)
#endif

  end subroutine gmres_solve

end module petsc_solver
