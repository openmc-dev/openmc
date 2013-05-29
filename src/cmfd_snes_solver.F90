module cmfd_snes_solver

# ifdef PETSC

  use cmfd_loss_operator,     only: loss_operator,init_M_operator,           &
                                    build_loss_matrix,destroy_M_operator
  use cmfd_prod_operator,     only: prod_operator,init_F_operator,           &
                                    build_prod_matrix,destroy_F_operator
  use cmfd_jacobian_operator, only: jacobian_operator,init_J_operator,       &
                                    build_jacobian_matrix,destroy_J_operator,&
                                    operators
  use cmfd_power_solver,      only: cmfd_power_execute

  implicit none
  private
  public :: cmfd_snes_execute

# include <finclude/petsc.h90>

  type(jacobian_operator) :: jac_prec
  type(operators) :: ctx

  Mat         :: jac                    ! jacobian matrix
  Vec         :: resvec                 ! residual vector
  Vec         :: xvec                   ! results
  KSP         :: ksp                    ! linear solver context
  PC          :: pc                     ! preconditioner
  SNES        :: snes                   ! nonlinear solver context
  integer     :: ierr                   ! error flag

  logical     :: adjoint_calc = .false. ! adjoint calculation

contains

!===============================================================================
! CMFD_SNES_EXECUTE
!===============================================================================

  subroutine cmfd_snes_execute(adjoint)

    logical, optional :: adjoint ! adjoint calculation

    ! check for adjoint
    if (present(adjoint)) adjoint_calc = adjoint

    ! seed with power iteration
    call cmfd_power_execute(k_tol=1.E-3_8, s_tol=1.E-3_8, adjoint=adjoint_calc)

    ! initialize data
    call init_data()

    ! initialize solver
    call init_solver()

    ! precondition
    call precondition_matrix()

    ! solve the system
    call SNESSolve(snes, PETSC_NULL_DOUBLE, xvec, ierr)

    ! extracts results to cmfd object
    call extract_results()

    ! deallocate all slepc data
    call finalize()

  end subroutine cmfd_snes_execute

!===============================================================================
! INIT_DATA allocates matrices vectors for CMFD solution
!===============================================================================

  subroutine init_data()

    use constants, only: ONE
    use global,    only: cmfd, n_procs_cmfd, rank

    integer              :: k         ! implied do counter
    integer              :: n         ! problem size
    integer              :: row_start ! local row start
    integer              :: row_end   ! local row end
    real(8), pointer     :: xptr(:)   ! solution pointer

    ! set up operator matrices
    call init_M_operator(ctx%loss)
    call init_F_operator(ctx%prod)
    call init_J_operator(jac_prec, ctx)

    ! get problem size
    n = jac_prec%n

    ! create PETSc vectors
    call VecCreateMPI(PETSC_COMM_WORLD, jac_prec%localn, PETSC_DECIDE, &
         resvec, ierr)
    call VecCreateMPI(PETSC_COMM_WORLD, jac_prec%localn, PETSC_DECIDE, &
         xvec, ierr)

    ! get the local dimensions for each process
    call VecGetOwnershipRange(xvec, row_start, row_end, ierr)

    if (rank == n_procs_cmfd - 1) row_end = n

    ! set flux in guess
    if (adjoint_calc) then
      call VecSetValues(xvec, row_end-row_start, (/(k,k=row_start,row_end-1)/), &
           cmfd%adj_phi(row_start+1:row_end), INSERT_VALUES, ierr)
    else
      call VecSetValues(xvec, row_end-row_start, (/(k,k=row_start,row_end-1)/), &
           cmfd%phi(row_start+1:row_end), INSERT_VALUES, ierr)
    end if
    call VecAssemblyBegin(xvec, ierr)
    call VecAssemblyEnd(xvec, ierr)

    ! set keff in guess
    if (rank == n_procs_cmfd - 1) then
      call VecGetArrayF90(xvec, xptr, ierr)
      if (adjoint_calc) then
        xptr(size(xptr)) = ONE/cmfd%adj_keff
      else
        xptr(size(xptr)) = ONE/cmfd%keff
      end if
      call VecRestoreArrayF90(xvec, xptr, ierr)
    end if

    ! nullify pointers
    if (associated(xptr)) nullify(xptr)

  end subroutine init_data

!===============================================================================
! INIT_SOLVER
!===============================================================================

  subroutine init_solver()

    use global,  only: cmfd_snes_monitor, n_procs_cmfd

    ! turn on mf_operator option
    call PetscOptionsSetValue("-snes_mf_operator", "TRUE", ierr)

    ! create SNES context
    call SNESCreate(PETSC_COMM_WORLD, snes, ierr)

    ! set the residual function
    call SNESSetFunction(snes, resvec, compute_nonlinear_residual, &
         PETSC_NULL_DOUBLE, ierr)

    ! set GMRES solver
    call SNESGetKSP(snes, ksp, ierr)
    call KSPSetType(ksp, KSPGMRES, ierr)

    ! create matrix free jacobian
    call MatCreateSNESMF(snes, jac, ierr)

    ! set matrix free finite difference
    call SNESSetJacobian(snes, jac, jac_prec, build_jacobian_matrix, ctx, ierr)

    ! set lags
    call SNESSetLagJacobian(snes, -2, ierr)
    call SNESSetLagPreconditioner(snes, -1, ierr)

    ! set convergence
!   call SNESSetTolerances(snes, 1.0e-8_8, 1.0e-8_8, 1.0e-8_8, 50, 10000, ierr)

    ! set SNES options
    call SNESSetFromOptions(snes, ierr)

    ! turn off line searching
!   call SNESLineSearchSet(snes, SNESLineSearchNo, PETSC_NULL_OBJECT, ierr)

  end subroutine init_solver

!===============================================================================
! PRECONDITION_MATRIX 
!===============================================================================

  subroutine precondition_matrix()

    use global,                 only: cmfd_snes_monitor, cmfd_solver_type, &
                                      n_procs_cmfd, cmfd_ilu_levels, master
    use string,                 only: to_str
    use, intrinsic :: ISO_FORTRAN_ENV

    character(LEN=20) :: snestype, ksptype, pctype

    ! set up preconditioner
    call KSPGetPC(ksp, pc, ierr)
    if (n_procs_cmfd == 1) then
      call PCSetType(pc, PCILU, ierr)
      call PCFactorSetLevels(pc, cmfd_ilu_levels, ierr)
    else
      call PetscOptionsSetValue("-pc_type", "bjacobi", ierr)
      call PetscOptionsSetValue("-sub_pc_type", "ilu", ierr)
      call PetscOptionsSetValue("-sub_pc_factor_levels", &
           trim(to_str(cmfd_ilu_levels)), ierr)
    end if

    ! get options
    call PCSetFromOptions(pc, ierr)
    call KSPSetFromOptions(ksp, ierr)

    ! finalize ksp setup
!   call KSPSetUp(ksp, ierr)

    ! get all types and print
    call SNESGetType(snes, snestype, ierr)
    call KSPGetType(ksp, ksptype, ierr)
    call PCGetType(pc, pctype, ierr)

    ! display solver info to user
    if (cmfd_snes_monitor .and. master) then
      write(OUTPUT_UNIT,'(A)') ''
      write(OUTPUT_UNIT,'(A)') '########################################################'
      write(OUTPUT_UNIT,'(A)') '################ JFNK Nonlinear Solver  ################'
      write(OUTPUT_UNIT,'(A)') '########################################################'
      write(OUTPUT_UNIT,'(A)') ''
      write(OUTPUT_UNIT,101) 'NONLINEAR SOLVER: ', snestype
      write(OUTPUT_UNIT,101) 'LINEAR SOLVER:    ', ksptype
      write(OUTPUT_UNIT,101) 'PRECONDITIONER:   ', pctype
      write(OUTPUT_UNIT,100) 'ILU levels: ', cmfd_ilu_levels
      write(OUTPUT_UNIT,'(A)') ''
      write(OUTPUT_UNIT,'(A)') '---------------------------------------------'
      write(OUTPUT_UNIT,'(A)') ''
    end if

 100 FORMAT(A,1X,I0)
 101 FORMAT(A,1X,A)
 
  end subroutine precondition_matrix

!===============================================================================
! COMPUTE_NONLINEAR_RESIDUAL
!===============================================================================

  subroutine compute_nonlinear_residual(snes, x, res, ierr)

    use global,       only: n_procs_cmfd, cmfd_write_matrices, &
                            cmfd_adjoint_type, rank

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
    PetscViewer :: viewer

    logical :: physical_adjoint = .false.

    ! check for physical adjoint
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'physical') &
         physical_adjoint = .true.

    ! create operators
    call build_loss_matrix(ctx%loss, adjoint = physical_adjoint)
    call build_prod_matrix(ctx%prod, adjoint = physical_adjoint)

    ! check for adjoint
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'math') &
         call compute_adjoint()

    ! get problem size
    n = ctx%loss%n 

    ! get pointers to vectors
    call VecGetArrayF90(x, xptr, ierr)
    call VecGetArrayF90(res, rptr, ierr)

    ! create petsc vector for flux
    call VecCreateMPI(PETSC_COMM_WORLD, ctx%loss%localn, &
         PETSC_DECIDE, phi, ierr)
    call VecCreateMPI(PETSC_COMM_WORLD, ctx%loss%localn, &
         PETSC_DECIDE, rphi, ierr)

    ! extract flux and place in petsc vector 
    call VecPlaceArray(phi, xptr, ierr)
    call VecPlaceArray(rphi, rptr, ierr)

    ! extract eigenvalue and broadcast (going to want to make this more general in future)
    if (rank == n_procs_cmfd - 1) lambda = xptr(size(xptr)) 
    call MPI_BCAST(lambda, 1, MPI_REAL8, n_procs_cmfd-1, PETSC_COMM_WORLD, ierr)

    ! create new petsc vectors to perform math
    call VecCreateMPI(PETSC_COMM_WORLD, ctx%loss%localn, PETSC_DECIDE, &
         phiM, ierr)

    ! calculate flux part of residual vector
    call MatMult(ctx%loss%M, phi, phiM, ierr)
    call MatMult(ctx%prod%F, phi, rphi, ierr)
    call VecAYPX(rphi, -lambda, phiM, ierr)

    ! set eigenvalue part of residual vector
    call VecDot(phi, phi, reslamb, ierr)

    ! map to ptr
    if (rank == n_procs_cmfd - 1) rptr(size(rptr)) = 0.5_8 - 0.5_8*reslamb

    ! reset arrays that are not used
    call VecResetArray(phi, ierr)
    call VecResetArray(rphi, ierr)

    ! restore arrays for residual and solution
    call VecRestoreArrayF90(x, xptr, ierr)
    call VecRestoreArrayF90(res, rptr, ierr)

    ! destroy all temp vectors
    call VecDestroy(phi, ierr)
    call VecDestroy(phiM, ierr)
    call VecDestroy(rphi, ierr)

    ! nullify all pointers
    if (associated(xptr)) nullify(xptr)
    if (associated(rptr)) nullify(rptr)

    ! write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
      if (adjoint_calc) then
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'adj_res.bin', &
             FILE_MODE_WRITE, viewer, ierr)
      else
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'res.bin', &
             FILE_MODE_WRITE, viewer, ierr)
      end if
      call VecView(res, viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)

      if (adjoint_calc) then
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'adj_x.bin', &
             FILE_MODE_WRITE, viewer, ierr)
      else
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'x.bin', &
             FILE_MODE_WRITE, viewer, ierr)
      end if
       call VecView(x, viewer, ierr)
       call PetscViewerDestroy(viewer, ierr)

    end if

  end subroutine compute_nonlinear_residual

!===============================================================================
! COMPUTE_ADJOINT
!===============================================================================

  subroutine compute_adjoint()

    use global,  only: cmfd_write_matrices

    PetscViewer :: viewer

    ! transpose matrices
    call MatTranspose(ctx%loss%M, MAT_REUSE_MATRIX, ctx%loss%M, ierr)
    call MatTranspose(ctx%prod%F, MAT_REUSE_MATRIX, ctx%prod%F, ierr)

    ! write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
      call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'adj_lossmat.bin', &
           FILE_MODE_WRITE, viewer, ierr)
      call MatView(ctx%loss%M, viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)

      call PetscViewerBinaryOpen(PETSC_COMM_WORLD, 'adj_prodmat.bin', &
           FILE_MODE_WRITE, viewer, ierr)
      call MatView(ctx%prod%F, viewer, ierr)
      call PetscViewerDestroy(viewer, ierr)
    end if

  end subroutine compute_adjoint

!===============================================================================
! EXTRACT_RESULTS 
!===============================================================================

  subroutine extract_results()

    use constants, only: ZERO, ONE
    use global,    only: cmfd, n_procs_cmfd, rank

    integer              :: n         ! problem size
    integer              :: row_start ! local row start
    integer              :: row_end   ! local row end
    real(8)              :: keff      ! keff of problem
    real(8),allocatable  :: mybuf(:)  ! temp buffer
    PetscScalar, pointer :: xptr(:)   ! pointer to eigenvector info

    ! get problem size
    n = ctx%loss%n

    ! also allocate in cmfd object
    if (adjoint_calc) then
      if (.not. allocated(cmfd%adj_phi)) allocate(cmfd%adj_phi(n))
    else
      if (.not. allocated(cmfd%phi)) allocate(cmfd%phi(n))
    end if
    if (.not. allocated(mybuf)) allocate(mybuf(n))

    ! get ownership range
    call VecGetOwnershipRange(xvec, row_start, row_end, ierr)

    ! resize the last proc
    if (rank == n_procs_cmfd - 1) row_end = row_end - 1

    ! convert petsc phi_object to cmfd_obj
    call VecGetArrayF90(xvec, xptr, ierr)
    if (adjoint_calc) then
      cmfd%adj_phi(row_start+1:row_end) = xptr(1:row_end - row_start)
    else
      cmfd%phi(row_start+1:row_end) = xptr(1:row_end - row_start) 
    end if

    ! reduce result to all 
    mybuf = ZERO
    if (adjoint_calc) then
      call MPI_ALLREDUCE(cmfd%adj_phi, mybuf, n, MPI_REAL8, MPI_SUM, &
           PETSC_COMM_WORLD, ierr)
    else
      call MPI_ALLREDUCE(cmfd%phi, mybuf, n, MPI_REAL8, MPI_SUM, &
           PETSC_COMM_WORLD, ierr)
    end if

    ! move buffer to object and deallocate
    cmfd%phi = mybuf
    if(allocated(mybuf)) deallocate(mybuf)

    ! save eigenvalue
    if(rank == n_procs_cmfd - 1) keff = ONE / xptr(size(xptr))
    if (adjoint_calc) then
      cmfd%adj_keff = keff
      call MPI_BCAST(cmfd%adj_keff, 1, MPI_REAL8, n_procs_cmfd-1, &
           PETSC_COMM_WORLD, ierr)
    else
      cmfd%keff = keff
      call MPI_BCAST(cmfd%keff, 1, MPI_REAL8, n_procs_cmfd-1, &
           PETSC_COMM_WORLD, ierr)
    end if
    call VecRestoreArrayF90(xvec, xptr, ierr)

    ! nullify pointers and deallocate local variables
    if (associated(xptr)) nullify(xptr)
    if (allocated(mybuf)) deallocate(mybuf)

  end subroutine extract_results

!===============================================================================
! FINALIZE 
!===============================================================================

  subroutine finalize()

    ! finalize data objects
    call destroy_M_operator(ctx%loss)
    call destroy_F_operator(ctx%prod)
    call destroy_J_operator(jac_prec)
    call VecDestroy(xvec, ierr)
    call VecDestroy(resvec, ierr)
    call MatDestroy(jac, ierr)

    ! finalize solver objects
    call SNESDestroy(snes, ierr)

  end subroutine finalize

# endif

end module cmfd_snes_solver
