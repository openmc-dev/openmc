module cmfd_power_solver

! This module contains routines to execute the power iteration solver

  use cmfd_loss_operator, only: init_loss_matrix, build_loss_matrix
  use cmfd_prod_operator, only: init_prod_matrix, build_prod_matrix
  use matrix_header,      only: Matrix
  use petsc_solver,       only: Petsc_gmres
  use vector_header,      only: Vector

  implicit none
  private
  public :: cmfd_power_execute 

  logical :: iconv                ! did the problem converged
  real(8) :: k_n                  ! new k-eigenvalue
  real(8) :: k_o                  ! old k-eigenvalue
  real(8) :: ktol = 1.e-8_8       ! tolerance on keff
  real(8) :: stol = 1.e-8_8       ! tolerance on source
  logical :: adjoint_calc = .false. ! run an adjoint calculation
  type(Matrix) :: loss            ! cmfd loss matrix
  type(Matrix) :: prod            ! cmfd prod matrix
  type(Vector) :: phi_n           ! new flux vector
  type(Vector) :: phi_o           ! old flux vector
  type(Vector) :: S_n             ! new source vector
  type(Vector) :: S_o             ! old flux vector
  type(Petsc_gmres) :: gmres      ! gmres solver

contains

!===============================================================================
! CMFD_POWER_EXECUTE sets up and runs power iteration solver for CMFD
!===============================================================================

  subroutine cmfd_power_execute(k_tol, s_tol, adjoint)

    use global,  only: cmfd_adjoint_type, time_cmfdbuild, time_cmfdsolve

    real(8), optional :: k_tol    ! tolerance on keff
    real(8), optional :: s_tol    ! tolerance on source
    logical, optional :: adjoint  ! adjoint calc

    logical :: physical_adjoint = .false.

    ! Set tolerances if present
    if (present(k_tol)) ktol = k_tol
    if (present(s_tol)) stol = s_tol

    ! Check for adjoint execution
    if (present(adjoint)) adjoint_calc = adjoint

    ! Check for physical adjoint
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'physical') &
        physical_adjoint = .true.

    ! Start timer for build
    call time_cmfdbuild % start()

    ! Initialize solver
    call gmres % create()

    ! Initialize matrices and vectors
    call init_data(physical_adjoint)

    ! Check for adjoint calculation
    if (adjoint_calc .and. trim(cmfd_adjoint_type) == 'math') &
        call compute_adjoint()

    ! Set up krylov info
    call gmres % set_oper(loss, loss)

    ! Precondition matrix
    call gmres % precondition(loss)

    ! Stop timer for build
    call time_cmfdbuild % stop()

    ! Begin power iteration 
    call time_cmfdsolve % start()
    call execute_power_iter()
    call time_cmfdsolve % stop()

    ! Extract results
    call extract_results()

    ! Deallocate petsc objects
    call finalize()

  end subroutine cmfd_power_execute

!===============================================================================
! INIT_DATA allocates matrices and vectors for CMFD solution
!===============================================================================

  subroutine init_data(adjoint)

    use constants, only: ONE

    logical :: adjoint

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
    call S_n % create(n)
    call S_o % create(n)

    ! Set initial guess
    guess = ONE
    phi_n % val = guess
    phi_o % val = guess
    k_n = guess
    k_o = guess

    ! Set up loss matrix
    call build_loss_matrix(loss, adjoint=adjoint) 

    ! Set up production matrix
    call build_prod_matrix(prod, adjoint=adjoint)

    ! Setup petsc for everything
    call loss % assemble()
    call prod % assemble()
    call loss % setup_petsc()
    call prod % setup_petsc()
    call phi_n % setup_petsc()
    call phi_o % setup_petsc()
    call S_o % setup_petsc()
    call S_n % setup_petsc()

  end subroutine init_data

!===============================================================================
! COMPUTE_ADJOINT computes a mathematical adjoint of CMFD problem 
!===============================================================================

  subroutine compute_adjoint()

    use global,  only: cmfd_write_matrices

    ! transpose matrices
    call loss % transpose()
    call prod % transpose()

    ! write out matrix in binary file (debugging)
    if (cmfd_write_matrices) then
      call loss % write_petsc_binary('adj_lossmat.bin')
      call prod % write_petsc_binary('adj_prodmat.bin')
    end if

  end subroutine compute_adjoint

!===============================================================================
! EXECUTE_POWER_ITER  in the main power iteration routine
!                     for the cmfd calculation
!===============================================================================

  subroutine execute_power_iter()

    integer     :: i         ! iteration counter

    ! reset convergence flag
    iconv = .false.

    ! begin power iteration
    do i = 1, 10000

      ! compute source vector
      call prod % vector_multiply(phi_o, S_o)

      ! normalize source vector
      S_o % val = S_o % val / k_o

      ! compute new flux vector
      call gmres % solve(S_o, phi_n)

      ! compute new source vector
      call prod % vector_multiply(phi_n, S_n)

      ! compute new k-eigenvalue
      k_n = sum(S_n % val) / sum(S_o % val)

      ! renormalize the old source
      S_o % val = S_o % val * k_o

      ! check convergence
      call convergence(i)

      ! to break or not to break
      if (iconv) exit

      ! record old values
      phi_o % val = phi_n % val
      k_o = k_n

    end do

  end subroutine execute_power_iter 

!===============================================================================
! CONVERGENCE checks the convergence of the CMFD problem
!===============================================================================

  subroutine convergence(iter)

    use constants,  only: ONE
    use global,     only: cmfd_power_monitor, master
    use, intrinsic :: ISO_FORTRAN_ENV

    integer     :: iter           ! iteration number

    real(8)     :: kerr           ! error in keff
    real(8)     :: serr           ! error in source

    ! reset convergence flag
    iconv = .false.

    ! calculate error in keff
    kerr = abs(k_o - k_n)/k_n

    ! calculate max error in source
    serr = sqrt(ONE/dble(S_n % n) * sum(((S_n % val - S_o % val)/S_n % val)**2))

    ! check for convergence
    if(kerr < ktol .and. serr < stol) iconv = .true.

    ! print out to user
    if (cmfd_power_monitor .and. master) then
      write(OUTPUT_UNIT,FMT='(I0,":",T10,"k-eff: ",F0.8,T30,"k-error: ",1PE12.5,T55, &
           "src-error: ",1PE12.5)') iter, k_n, kerr, serr
    end if

  end subroutine convergence

!===============================================================================
! EXTRACT_RESULTS takes results and puts them in CMFD global data object
!===============================================================================

  subroutine extract_results()

    use global, only: cmfd, cmfd_write_matrices

    character(len=25)    :: filename  ! name of file to write data 
    integer              :: n         ! problem size

    ! get problem size
    n = loss % n

    ! also allocate in cmfd object
    if (adjoint_calc) then
      if (.not. allocated(cmfd%adj_phi)) allocate(cmfd%adj_phi(n))
    else
      if (.not. allocated(cmfd%phi)) allocate(cmfd%phi(n))
    end if

    ! save values 
    if (adjoint_calc) then
      cmfd % adj_phi = phi_n % val
    else
      cmfd % phi = phi_n % val 
    end if

    ! save eigenvalue
    if(adjoint_calc) then
      cmfd%adj_keff = k_n
    else
      cmfd%keff = k_n
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

    ! destroy all objects 
    call gmres % destroy() 
    call loss  % destroy() 
    call prod  % destroy()
    call phi_n % destroy()
    call phi_o % destroy()
    call S_n   % destroy()
    call S_o   % destroy()

  end subroutine finalize

end module cmfd_power_solver
