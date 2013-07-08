module petsc_solver

  implicit none
  private

#  include <finclude/petsc.h90>

  type, public :: Petsc_gmres 
    KSP :: ksp
    PC  :: pc
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

    call KSPCreate(PETSC_COMM_WORLD, self % ksp, petsc_err)
    call KSPSetTolerances(self % ksp, rtol, atol, &
         PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER, petsc_err)
    call KSPSetType(self % ksp, KSPGMRES, petsc_err)
    call KSPSetInitialGuessNonzero(self % ksp, PETSC_TRUE, petsc_err)
    call KSPGetPC(self % ksp, self % pc, petsc_err)
    call PCFactorSetLevels(self % pc, ilu_levels, petsc_err)

  end subroutine gmres_create

!===============================================================================
! GMRES_PRECONDITION
!===============================================================================

  subroutine gmres_precondition(self, matrix)

    class(Petsc_gmres) :: self
    Mat          :: matrix

    call KSPSetUp(self % ksp, petsc_err) 
    call PCFactorGetMatrix(self % pc, matrix, petsc_err) 

  end subroutine gmres_precondition

!===============================================================================
! GMRES_SET_OPER
!===============================================================================

  subroutine gmres_set_oper(self, prec_matrix, matrix)

    class(Petsc_gmres) :: self
    Mat                :: prec_matrix
    Mat                :: matrix

    call KSPSetOperators(self % ksp, matrix, prec_matrix, &
         SAME_NONZERO_PATTERN, petsc_err)

  end subroutine gmres_set_oper

!===============================================================================
! GMRES_DESTROY
!===============================================================================

  subroutine gmres_destroy(self)

   class(Petsc_gmres) :: self

    call KSPDestroy(self % ksp, petsc_err)

  end subroutine gmres_destroy

!===============================================================================
! GMRES_SOLVE
!===============================================================================

  subroutine gmres_solve(self, b, x)

    class(Petsc_gmres) :: self
    Vec :: b
    Vec :: x

    integer :: petsc_err

    call KSPSolve(self % ksp, b, x, petsc_err)

  end subroutine gmres_solve

end module petsc_solver
