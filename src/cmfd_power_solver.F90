module cmfd_power_solver

! This module contains routines to execute the power iteration solver

  use cmfd_loss_operator, only: loss_operator,init_M_operator,                 &
 &                        build_loss_matrix,destroy_M_operator
  use cmfd_prod_operator, only: prod_operator,init_F_operator,                 &
 &                        build_prod_matrix,destroy_F_operator

  implicit none
  private
  public :: cmfd_power_execute 

# include <finclude/petsc.h90>

  type(loss_operator) :: loss   ! M loss matrix
  type(prod_operator) :: prod   ! F production matrix

  Vec         :: phi_n                ! new flux eigenvector
  Vec         :: phi_o                ! old flux eigenvector
  Vec         :: S_n                  ! new source vector
  Vec         :: S_o                  ! old source vector
  KSP         :: krylov               ! krylov solver
  KSP         :: sub_krylov           ! sub-ksp for bjacobi
  PC          :: prec                 ! preconditioner for krylov
  PC          :: sub_prec             ! sub-prec for bjacobi
  integer     :: ierr                 ! error flag
  integer     :: nlocal
  integer     :: first
  real(8)     :: k_n                  ! new k-eigenvalue
  real(8)     :: k_o                  ! old k-eigenvalue
  real(8)     :: ktol = 1.e-8_8       ! tolerance on keff
  real(8)     :: stol = 1.e-8_8       ! tolerance on source
  logical     :: iconv                ! did the problem converged

  logical     :: adjoint_calc = .FALSE. ! run an adjoint calculation

contains

!===============================================================================
! CMFD_POWER_EXECUTE
!===============================================================================

  subroutine cmfd_power_execute(k_tol,s_tol,adjoint)

    use global,  only: cmfd_adjoint_type, n_procs_cmfd

    real(8), optional :: k_tol    ! tolerance on keff
    real(8), optional :: s_tol    ! tolerance on source
    logical, optional :: adjoint  ! adjoint calc

    logical :: physical_adjoint = .FALSE.

    ! set tolerances if present
    if (present(k_tol)) ktol = k_tol
    if (present(s_tol)) stol = s_tol

    ! check for adjoint execution
    if (present(adjoint)) adjoint_calc = adjoint

    ! check for physical adjoint
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'physical')          &
        physical_adjoint = .TRUE.

    ! initialize solver
    call init_solver()

    ! initialize matrices and vectors
    call init_data()

    ! set up M loss matrix
    call build_loss_matrix(loss, adjoint=physical_adjoint) 

    ! set up F production matrix
    call build_prod_matrix(prod, adjoint=physical_adjoint)

    ! check for adjoint calculation
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'math')              &
        call compute_adjoint()

    ! set up krylov info
    call KSPSetOperators(krylov, loss%M, loss%M, SAME_NONZERO_PATTERN, ierr)

    ! precondition matrix
    call precondition_matrix()

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
    n = loss%localn

    ! set up flux vectors
    call VecCreateMPI(PETSC_COMM_WORLD,n,PETSC_DECIDE,phi_n,ierr)
    call VecCreateMPI(PETSC_COMM_WORLD,n,PETSC_DECIDE,phi_o,ierr)

    ! set up source vectors
    call VecCreateMPI(PETSC_COMM_WORLD,n,PETSC_DECIDE,S_n,ierr)
    call VecCreateMPI(PETSC_COMM_WORLD,n,PETSC_DECIDE,S_o,ierr)

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

    real(8)  :: rtol ! relative tolerance
    real(8)  :: atol ! absolute tolerance  

    ! set tolerance
    rtol = 1.0e-10_8
    atol = 1.0e-10_8

    ! set up krylov solver
    call KSPCreate(PETSC_COMM_WORLD,krylov,ierr)
    call KSPSetTolerances(krylov,rtol,atol,                                    &
   &                      PETSC_DEFAULT_DOUBLE_PRECISION,                      &
   &                      PETSC_DEFAULT_INTEGER,ierr)
    call KSPSetType(krylov,KSPGMRES,ierr)
    call KSPSetInitialGuessNonzero(krylov,PETSC_TRUE,ierr)

  end subroutine init_solver

!===============================================================================
! PRECONDITION_MATRIX 
!===============================================================================

  subroutine precondition_matrix()

    use global,                 only: cmfd_power_monitor, cmfd_solver_type,    &
                                      n_procs_cmfd, cmfd_ilu_levels, master
    use string,                 only: to_str

    character(len=20) :: ksptype,pctype
 
    ! set up preconditioner
    call KSPGetPC(krylov,prec,ierr)
    if (n_procs_cmfd == 1) then
      call PCSetType(prec,PCILU,ierr)
      call PCFactorSetLevels(prec,cmfd_ilu_levels,ierr)
      call KSPSetUp(krylov,ierr)
      call PCFactorGetMatrix(prec,loss%M,ierr)
    else
      call PetscOptionsSetValue("-pc_type","bjacobi",ierr)
      call PetscOptionsSetValue("-sub_pc_type","ilu",ierr)  
      call PetscOptionsSetValue("-sub_pc_factor_levels",trim(to_str(         &
                                cmfd_ilu_levels)),ierr)
      call PCSetFromOptions(prec,ierr)
      call KSPSetUp(krylov,ierr)
    end if

    ! get options
    if (trim(cmfd_solver_type) == 'power') call KSPSetFromOptions(krylov,ierr)

    ! get all types and print
    call KSPGetType(krylov,ksptype,ierr)
    call PCGetType(prec,pctype,ierr)

    ! print heading information
    if (cmfd_power_monitor .and. master) then
      write(*,*)
      write(*,*) '########################################################'
      write(*,*) '################ Power Iteration Solver ################'
      write(*,*) '########################################################'
      write(*,*)
      write(*,*) 'Eigenvalue Tolerance:',ktol
      write(*,*) 'Source Tolerance:    ',stol
      write(*,*)
      write(*,*) 'Linear Solver Type:  ',ksptype
      write(*,*) 'Preconditioner Type: ',pctype
      write(*,*) 'ILU Fill Levels:',cmfd_ilu_levels
      write(*,*)
      write(*,*) '---------------------------------------------'
      write(*,*)
    end if

  end subroutine precondition_matrix

!===============================================================================
! COMPUTE_ADJOINT 
!===============================================================================

  subroutine compute_adjoint()

    use global,  only: cmfd_write_matrices

    PetscViewer :: viewer

    ! transpose matrices
    call MatTranspose(loss%M,MAT_REUSE_MATRIX,loss%M,ierr)
    call MatTranspose(prod%F,MAT_REUSE_MATRIX,prod%F,ierr)

    ! write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
      call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'adj_lossmat.bin'            &
     &     ,FILE_MODE_WRITE,viewer,ierr)
      call MatView(loss%M,viewer,ierr)
      call PetscViewerDestroy(viewer,ierr)

      call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'adj_prodmat.bin'            &
     &     ,FILE_MODE_WRITE,viewer,ierr)
      call MatView(prod%F,viewer,ierr)
      call PetscViewerDestroy(viewer,ierr)
    end if

  end subroutine compute_adjoint

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
      call convergence(i)

      ! to break or not to break
      if (iconv) exit

      ! record old values
      call VecCopy(phi_n,phi_o,ierr)
      k_o = k_n

    end do

  end subroutine execute_power_iter 

!===============================================================================
! CONVERGENCE
!===============================================================================

  subroutine convergence(iter)

    use global,  only: cmfd_power_monitor, master

    integer     :: iter           ! iteration number

    real(8)     :: kerr           ! error in keff
    real(8)     :: serr           ! error in source
    real(8)     :: norm_n         ! L2 norm of new source
    real(8)     :: norm_o         ! L2 norm of old source
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

    ! print out to user
    if (cmfd_power_monitor .and. master) then
      write(*,FMT='(I0,":",T10,"k-eff: ",F0.8,T30,"k-error: ",1PE12.5,T55,     &
                    "src-error: ",1PE12.5)') iter,k_n,kerr,serr
    end if

  end subroutine convergence

!===============================================================================
! EXTRACT_RESULTS 
!===============================================================================

  subroutine extract_results()

    use global,           only: cmfd, n_procs_cmfd, cmfd_write_matrices

    integer              :: n         ! problem size
    integer              :: row_start ! local row start
    integer              :: row_end   ! local row end
    real(8),allocatable  :: mybuf(:)  ! temp buffer
    PetscScalar, pointer :: phi_v(:)  ! pointer to eigenvector info
    PetscViewer          :: viewer    ! petsc viewer for binary write

    ! get problem size
    n = loss%n

    ! also allocate in cmfd object
    if (adjoint_calc) then
      if (.not. allocated(cmfd%adj_phi)) allocate(cmfd%adj_phi(n))
    else
      if (.not. allocated(cmfd%phi)) allocate(cmfd%phi(n))
    end if
    if (.not. allocated(mybuf)) allocate(mybuf(n))

    ! get ownership range
    call VecGetOwnershipRange(phi_n,row_start,row_end,ierr)

    ! convert petsc phi_object to cmfd_obj
    call VecGetArrayF90(phi_n,phi_v,ierr)
    if (adjoint_calc) then
      cmfd%adj_phi(row_start+1:row_end) = phi_v
    else
      cmfd%phi(row_start+1:row_end) = phi_v 
    end if
    call VecRestoreArrayF90(phi_n,phi_v,ierr)

    ! save eigenvalue
    if(adjoint_calc) then
      cmfd%adj_keff = k_n
    else
      cmfd%keff = k_n
    end if

    ! reduce result to all 
    mybuf = 0.0_8
    if (adjoint_calc) then
      call MPI_ALLREDUCE(cmfd%adj_phi,mybuf,n,MPI_REAL8,MPI_SUM,PETSC_COMM_WORLD,ierr)
    else
      call MPI_ALLREDUCE(cmfd%phi,mybuf,n,MPI_REAL8,MPI_SUM,PETSC_COMM_WORLD,ierr)
    end if

    ! normalize phi to 1
    if (adjoint_calc) then
      cmfd%adj_phi = cmfd%adj_phi/sqrt(sum(cmfd%adj_phi*cmfd%adj_phi))
    else
      cmfd%phi = cmfd%phi/sqrt(sum(cmfd%phi*cmfd%phi))
    end if

    ! write out results
    if (cmfd_write_matrices) then
      if (adjoint_calc) then
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'adj_fluxvec.bin' &
       &     ,FILE_MODE_WRITE,viewer,ierr)
      else
        call PetscViewerBinaryOpen(PETSC_COMM_WORLD,'fluxvec.bin' &
       &     ,FILE_MODE_WRITE,viewer,ierr)
      end if
      call VecView(phi_n,viewer,ierr)
      call PetscViewerDestroy(viewer,ierr)
    end if

    ! nullify pointer and deallocate local vars
    if (associated(phi_v)) nullify(phi_v)
    if (allocated(mybuf)) deallocate(mybuf)

  end subroutine extract_results

!===============================================================================
! FINALIZE 
!===============================================================================

  subroutine finalize()

    use global,  only: n_procs_cmfd

    ! finalize solver objects
    call KSPDestroy(krylov,ierr)
    call KSPDestroy(sub_krylov,ierr)

    ! finalize data objects
    if (n_procs_cmfd > 1) call destroy_M_operator(loss) ! only destroy for jacobi
    call destroy_F_operator(prod)

    call VecDestroy(phi_n,ierr)
    call VecDestroy(phi_o,ierr)
    call VecDestroy(S_n,ierr)
    call VecDestroy(S_o,ierr)

  end subroutine finalize

end module cmfd_power_solver
