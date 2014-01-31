module cmfd_power_solver

! This module contains routines to execute the power iteration solver

  use cmfd_loss_operator, only: init_loss_matrix, build_loss_matrix
  use cmfd_prod_operator, only: init_prod_matrix, build_prod_matrix
  use matrix_header,      only: Matrix
  use solver_interface,   only: GMRESSolver
  use vector_header,      only: Vector

  implicit none
  private
  public :: cmfd_power_execute 

  logical :: iconv                ! did the problem converged
  real(8) :: k_n                  ! new k-eigenvalue
  real(8) :: k_o                  ! old k-eigenvalue
  real(8) :: ktol = 1.e-8_8       ! tolerance on keff
  real(8) :: stol = 1.e-8_8       ! tolerance on source
  real(8) :: norm_n               ! current norm of source vector
  real(8) :: norm_o               ! old norm of source vector
  real(8) :: kerr                 ! error in keff
  real(8) :: serr                 ! error in source
  logical :: adjoint_calc         ! run an adjoint calculation
  type(Matrix) :: loss            ! cmfd loss matrix
  type(Matrix) :: prod            ! cmfd prod matrix
  type(Vector) :: phi_n           ! new flux vector
  type(Vector) :: phi_o           ! old flux vector
  type(Vector) :: s_n             ! new source vector
  type(Vector) :: s_o             ! old flux vector
  type(Vector) :: serr_v          ! error in source
  type(GMRESSolver) :: gmres      ! gmres solver

contains

!===============================================================================
! CMFD_POWER_EXECUTE sets up and runs power iteration solver for CMFD
!===============================================================================

  subroutine cmfd_power_execute(k_tol, s_tol, adjoint)

    use global,  only: cmfd_adjoint_type, time_cmfdbuild, time_cmfdsolve

    real(8), intent(in), optional :: k_tol    ! tolerance on keff
    real(8), intent(in), optional :: s_tol    ! tolerance on source
    logical, intent(in), optional :: adjoint  ! adjoint calc

    logical :: physical_adjoint = .false. ! physical adjoint default false

    ! Set tolerances if present
    if (present(k_tol)) ktol = k_tol
    if (present(s_tol)) stol = s_tol

    ! Check for adjoint execution
    adjoint_calc = .false.
    if (present(adjoint)) adjoint_calc = adjoint

    ! Check for physical adjoint
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'physical') &
        physical_adjoint = .true.

    ! Start timer for build
    call time_cmfdbuild % start()

    ! Initialize solver
#ifdef PETSC
    call gmres % create()
#endif

    ! Initialize matrices and vectors
    call init_data(physical_adjoint)

    ! Check for mathematical adjoint calculation
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'math') &
        call compute_adjoint()

    ! Set up krylov info
#ifdef PETSC
    call gmres % set_oper(loss, loss)
#endif

    ! Stop timer for build
    call time_cmfdbuild % stop()

    ! Begin power iteration 
    call time_cmfdsolve % start()
    call execute_power_iter()
    call time_cmfdsolve % stop()

    ! Extract results
    call extract_results()

    ! Deallocate data 
    call finalize()

  end subroutine cmfd_power_execute

!===============================================================================
! INIT_DATA allocates matrices and vectors for CMFD solution
!===============================================================================

  subroutine init_data(adjoint)

    use constants, only: ONE, ZERO
#ifdef PETSC
    use global,    only: cmfd_write_matrices
#endif

    logical, intent(in) :: adjoint ! adjoint calcualtion

    integer :: n      ! problem size
    real(8) :: guess  ! initial guess

    ! Set up matrices
    call init_loss_matrix(loss)
    call init_prod_matrix(prod)

    ! Get problem size
    n = loss % n

    ! Set up flux vectors
    call phi_n % create(n)
    call phi_o % create(n)

    ! Set up source vectors
    call s_n % create(n)
    call s_o % create(n)
    call serr_v % create(n)

    ! Set initial guess
    guess = ONE
    phi_n % val = guess
    phi_o % val = guess
    k_n = guess
    k_o = guess

    ! Fill in loss matrix
    call build_loss_matrix(loss, adjoint=adjoint) 

    ! Fill in production matrix
    call build_prod_matrix(prod, adjoint=adjoint)

    ! Setup petsc for everything
    call loss % assemble()
    call prod % assemble()
#ifdef PETSC
    call loss % setup_petsc()
    call prod % setup_petsc()
    call phi_n % setup_petsc()
    call phi_o % setup_petsc()
    call s_o % setup_petsc()
    call s_n % setup_petsc()
    if (cmfd_write_matrices) call loss % write_petsc_binary('loss.bin')
    if (cmfd_write_matrices) call prod % write_petsc_binary('prod.bin')
#endif

    ! Set norms to 0
    norm_n = ZERO
    norm_o = ZERO

  end subroutine init_data

!===============================================================================
! COMPUTE_ADJOINT computes a mathematical adjoint of CMFD problem 
!===============================================================================

  subroutine compute_adjoint()

    use global,  only: cmfd_write_matrices

    ! Transpose matrices
    call loss % transpose()
    call prod % transpose()

    ! Write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
      call loss % write_petsc_binary('adj_lossmat.bin')
      call prod % write_petsc_binary('adj_prodmat.bin')
    end if

  end subroutine compute_adjoint

!===============================================================================
! EXECUTE_POWER_ITER  is the main power iteration routine
!                     for the cmfd calculation
!===============================================================================

  subroutine execute_power_iter()

    integer :: i ! iteration counter

    ! Reset convergence flag
    iconv = .false.

    ! Begin power iteration
    do i = 1, 10000

      ! Compute source vector
      call prod % vector_multiply(phi_o, s_o)

      ! Normalize source vector
      s_o % val = s_o % val / k_o

      ! Compute new flux vector
#ifdef PETSC
      call gmres % solve(s_o, phi_n)
#endif

      ! Compute new source vector
      call prod % vector_multiply(phi_n, s_n)

      ! Compute new k-eigenvalue
      k_n = sum(s_n % val) / sum(s_o % val)

      ! Renormalize the old source
      s_o % val = s_o % val * k_o

      ! Check convergence
      call convergence(i)

      ! Break loop if converged
      if (iconv) exit

      ! Record old values
      phi_o % val = phi_n % val
      k_o = k_n
      norm_o = norm_n

    end do

  end subroutine execute_power_iter 

!===============================================================================
! CONVERGENCE checks the convergence of the CMFD problem
!===============================================================================

  subroutine convergence(iter)

    use constants,  only: ONE, TINY_BIT
    use global,     only: cmfd_power_monitor, master
    use, intrinsic :: ISO_FORTRAN_ENV

    integer, intent(in) :: iter ! iteration number

    ! Reset convergence flag
    iconv = .false.

    ! Calculate error in keff
    kerr = abs(k_o - k_n)/k_n

    ! Calculate max error in source
    where (s_n % val > TINY_BIT)
      serr_v % val = ((s_n % val - s_o % val)/s_n % val)**2
    end where
    serr = sqrt(ONE/dble(s_n % n) * sum(serr_v % val))

    ! Check for convergence
    if(kerr < ktol .and. serr < stol) iconv = .true.

    ! Save the L2 norm of the source
    norm_n = serr

    ! Print out to user
    if (cmfd_power_monitor .and. master) then
      write(OUTPUT_UNIT,FMT='(I0,":",T10,"k-eff: ",F0.8,T30,"k-error: ", &
           &1PE12.5,T55, "src-error: ",1PE12.5)') iter, k_n, kerr, serr
    end if

  end subroutine convergence

!===============================================================================
! EXTRACT_RESULTS takes results and puts them in CMFD global data object
!===============================================================================

  subroutine extract_results()

    use global, only: cmfd, cmfd_write_matrices, current_batch

    character(len=25)    :: filename  ! name of file to write data 
    integer              :: n         ! problem size

    ! Get problem size
    n = loss % n

    ! Allocate in cmfd object if not already allocated
    if (adjoint_calc) then
      if (.not. allocated(cmfd%adj_phi)) allocate(cmfd%adj_phi(n))
    else
      if (.not. allocated(cmfd%phi)) allocate(cmfd%phi(n))
    end if

    ! Save values 
    if (adjoint_calc) then
      cmfd % adj_phi = phi_n % val
    else
      cmfd % phi = phi_n % val 
    end if

    ! Save eigenvalue
    if(adjoint_calc) then
      cmfd%adj_keff = k_n
    else
      cmfd%keff = k_n
    end if

    ! Normalize phi to 1
    if (adjoint_calc) then
      cmfd%adj_phi = cmfd%adj_phi/sqrt(sum(cmfd%adj_phi*cmfd%adj_phi))
    else
      cmfd%phi = cmfd%phi/sqrt(sum(cmfd%phi*cmfd%phi))
    end if

    ! Save dominance ratio
    cmfd % dom(current_batch) = norm_n/norm_o

    ! Write out results
    if (cmfd_write_matrices) then
      if (adjoint_calc) then
        filename = 'adj_fluxvec.bin'
      else
        filename = 'fluxvec.bin'
      end if
      call phi_n % write_petsc_binary(filename)
    end if

  end subroutine extract_results

!===============================================================================
! FINALIZE frees all memory associated with power iteration
!===============================================================================

  subroutine finalize()

    ! Destroy all objects 
#ifdef PETSC
    call gmres  % destroy()
#endif
    call loss   % destroy() 
    call prod   % destroy()
    call phi_n  % destroy()
    call phi_o  % destroy()
    call s_n    % destroy()
    call s_o    % destroy()
    call serr_v % destroy

  end subroutine finalize

end module cmfd_power_solver
