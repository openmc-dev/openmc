module petsc_solver

  implicit none
  private

# ifdef PETSC
#  include <finclude/petsc.h90>
# endif

  type, public :: gmres
#  ifdef PETSC
    KSP :: ksp
    PC  :: pc
#  endif
   contains
     procedure :: create       => gmres_create
     procedure :: precondition => gmres_precondition
     procedure :: set_oper     => gmres_set_oper
     procedure :: destroy      => gmres_destroy
  end type gmres

  integer :: petsc_err

contains

!===============================================================================
! GMRES_CREATE
!===============================================================================

  subroutine gmres_create(self)

    class(gmres) :: self

    integer :: ilu_levels = 5
    real(8) :: rtol = 1.0e-10_8
    real(8) :: atol = 1.0e-10_8

#  ifdef PETSC
    call KSPCreate(PETSC_COMM_WORLD, self % ksp, petsc_err)
    call KSPSetTolerances(self % ksp, rtol, atol, &
         PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER, petsc_err)
    call KSPSetType(self % ksp, KSPGMRES, petsc_err)
    call KSPSetInitialGuessNonzero(self % ksp, PETSC_TRUE, petsc_err)
    call KSPGetPC(self % ksp, self % pc, petsc_err)
    call PCFactorSetLevels(self % pc, ilu_levels, petsc_err)
    call KSPSetUp(self % ksp, petsc_err) 
#  endif

  end subroutine gmres_create

!===============================================================================
! GMRES_PRECONDITION
!===============================================================================

  subroutine gmres_precondition(self, matrix)

    class(gmres) :: self
    Mat          :: matrix

#  ifdef PETSC
    call PCFactorGetMatrix(self % pc, matrix, petsc_err) 
#  endif

  end subroutine gmres_precondition

!===============================================================================
! GMRES_SET_OPER
!===============================================================================

  subroutine gmres_set_oper(self, prec_matrix, matrix)

    class(gmres) :: self
    Mat          :: prec_matrix
    Mat          :: matrix

#  ifdef PETSC
    call KSPSetOperators(self % ksp, matrix, prec_matrix, &
         SAME_NONZERO_PATTERN, petsc_err)
#  endif

  end subroutine gmres_set_oper

!===============================================================================
! GMRES_DESTROY
!===============================================================================

  subroutine gmres_destroy(self)

   class(gmres) :: self

#  ifdef PETSC
    call KSPDestroy(self % ksp, petsc_err)
#  endif 

  end subroutine gmres_destroy

end module petsc_solver 
