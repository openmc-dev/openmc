module URR_edits

  use URR_constants,      only: ISOTROPIC,&
       EQUIPROBABLE_32_BINS,&
       TABULAR
  use URR_error, only: exit_status,&
                       log_message,&
                       INFO
  use URR_openmc_wrapper, only: Material,&
       materials,&
       e_grid,&
       Nuclide,&
       nuclides,&
       to_str,&
       n_materials,&
       Reaction,&
       prn,&
       fatal_error
  use URR_isotope,        only: Isotope,&
                                isotopes
  use URR_resonance,      only: wigner_level_spacing

  implicit none
  private
  public :: write_angle,&
            write_xs

contains


!> Write abscissa-ordinate pairs to a text file
  subroutine write_coords(unit_num, filename, x, y)

    integer :: unit_num ! unit number for output file
    integer :: i        ! coordinate pair index
    character(*) :: filename ! name of output file
    real(8) :: x(:) ! vector of abscissae
    real(8) :: y(:) ! vector of ordinates

    open(unit = unit_num, file = trim(adjustl(filename)))

    if (size(x) == size(y)) then
      do i = 1, size(x)
        write(unit_num, '(ES24.16, ES24.16)') x(i), y(i)
      end do
    else
      call fatal_error('Mismatched vector lengths in write_coords')
    end if

    close(unit_num)

  end subroutine write_coords


!> Write out user-requested pointwise, continuous-energy cross sections
  subroutine write_xs()

    integer :: i      ! material index
    integer :: j      ! nuclide index
    integer :: k      ! loop index for list of xs grids to be written
    integer :: i_E    ! energy point index
    integer :: i_nuc  ! index in global nuclides array
    integer :: x_size ! length of vector of abscissae
    integer :: y_size ! length of vector of ordinates
    real(8), allocatable :: x_vals(:)        ! vector of abscissae
    real(8), allocatable :: y_vals(:)        ! vector of ordinates
    character(80) :: filename                ! name of xs grid output file
    type(Material), pointer :: mat  ! material pointer
    type(Nuclide),  pointer :: nuc  ! nuclide pointer
    type(Isotope),  pointer :: tope ! isotope pointer

    ! loop over all materials
    do i = 1, n_materials
      mat => materials(i)

      if (mat % nat_elements) then
        call log_message(INFO, 'cross_sections element(s) in material(s) with &
             &natural elements will be ignored')
        cycle
      end if

      ! loop over all nuclides in material
      do j = 1, mat % n_nuclides
        i_nuc = mat % nuclide(j)
        nuc => nuclides(i_nuc)
        if (nuc % i_isotope /= 0) then
          tope => isotopes(nuc % i_isotope)
        else
          nullify(tope)
        end if

        ! loop over all requested energy-xs grid outputs
        do k = 1, mat % write_xs(j) % size()

          ! determine which energy-xs grid is requested
          select case (trim(adjustl(mat % write_xs(j) % get_item(k))))

          ! write unionized energy grid values
          case ('unionized')

            ! unionized energy grid values
            x_size = size(e_grid)
            allocate(x_vals(x_size))

            ! if we're writing the unionized grid energy values, xs values don't
            ! matter, so give the energy values as both the abscissae and
            ! ordinates 
            x_vals = e_grid
            y_size = x_size
            allocate(y_vals(y_size))
            y_vals = x_vals
            filename = "unionized-energy-grid.dat"

          ! write a nuclide's energy-total xs grid values
          case ('total')

            x_size = size(nuc % energy)
            allocate(x_vals(x_size))
            x_vals = nuc % energy
            y_size = size(nuc % total)
            allocate(y_vals(y_size))
            y_vals = nuc % total
            filename = trim(adjustl(nuc % name)) // "-total.dat"

          ! write a nuclide's energy-elastic xs grid values
          case ('elastic')

            x_size = size(nuc % energy)
            allocate(x_vals(x_size))
            x_vals = nuc % energy
            y_size = size(nuc % elastic)
            allocate(y_vals(y_size))
            y_vals = nuc % elastic
            filename = trim(adjustl(nuc % name)) // "-elastic.dat"

          ! write a nuclide's energy-fission xs grid values
          case ('fission')

            x_size = size(nuc % energy)
            allocate(x_vals(x_size))
            x_vals = nuc % energy
            y_size = size(nuc % fission)
            allocate(y_vals(y_size))
            y_vals = nuc % fission
            filename = trim(adjustl(nuc % name)) // "-fission.dat"

          ! write a nuclide's energy-capture xs grid values
          case ('capture')

            x_size = size(nuc % energy)
            allocate(x_vals(x_size))
            x_vals = nuc % energy
            y_size = size(nuc % absorption)
            allocate(y_vals(y_size))
            y_vals = nuc % absorption - nuc % fission
            filename = trim(adjustl(nuc % name)) // "-capture.dat"

          ! write a nuclide's energy-capture xs grid values
          case ('inelastic')

            x_size = size(nuc % energy)
            allocate(x_vals(x_size))
            x_vals = nuc % energy
            y_size = size(nuc % total)
            allocate(y_vals(y_size))
            y_vals = nuc % total - nuc % absorption - nuc % elastic
            filename = trim(adjustl(nuc % name)) // "-inelastic.dat"

          ! write a nuclide's URR energy-total xs grid values
          case ('urr-total')

            x_size = size(tope % urr_E)
            allocate(x_vals(x_size))
            x_vals = tope % urr_E
            y_size = size(tope % point_xs)
            allocate(y_vals(y_size))
            do i_E = 1, y_size
              y_vals(i_E) = tope % point_xs(i_E) % t
            end do
            filename = trim(adjustl(to_str(tope % ZAI))) // "-urr-total.dat"

          ! write a nuclide's URR energy-elastic xs grid values
          case ('urr-elastic')

            x_size = size(tope % urr_E)
            allocate(x_vals(x_size))
            x_vals = tope % urr_E
            y_size = size(tope % point_xs)
            allocate(y_vals(y_size))
            do i_E = 1, y_size
              y_vals(i_E) = tope % point_xs(i_E) % n
            end do
            filename = trim(adjustl(to_str(tope % ZAI))) // "-urr-elastic.dat"

          ! write a nuclide's URR energy-fission xs grid values
          case ('urr-fission')

            x_size = size(tope % urr_E)
            allocate(x_vals(x_size))
            x_vals = tope % urr_E
            y_size = size(tope % point_xs)
            allocate(y_vals(y_size))
            do i_E = 1, y_size
              y_vals(i_E) = tope % point_xs(i_E) % f
            end do
            filename = trim(adjustl(to_str(tope % ZAI))) // "-urr-fission.dat"

          ! write a nuclide's URR energy-inelastic xs grid values
          case ('urr-inelastic')

            x_size = size(tope % urr_E)
            allocate(x_vals(x_size))
            x_vals = tope % urr_E
            y_size = size(tope % point_xs)
            allocate(y_vals(y_size))
            do i_E = 1, y_size
              y_vals(i_E) = tope % point_xs(i_E) % x
            end do
            filename = trim(adjustl(to_str(tope % ZAI))) // "-urr-inelastic.dat"

          ! write a nuclide's URR energy-capture xs grid values
          case ('urr-capture')

            x_size = size(tope % urr_E)
            allocate(x_vals(x_size))
            x_vals = tope % urr_E
            y_size = size(tope % point_xs)
            allocate(y_vals(y_size))
            do i_E = 1, y_size
              y_vals(i_E) = tope % point_xs(i_E) % g
            end do
            filename = trim(adjustl(to_str(tope % ZAI))) // "-urr-capture.dat"

          ! the requested xs is not recognized
          case default
            call fatal_error('Not an allowed energy-xs grid option')
          end select

          ! write the energy-xs value pairs to a file
          call write_coords(99, &
            & filename, &
            & x_vals, &
            & y_vals)
          deallocate(x_vals)
          if (allocated(y_vals)) deallocate(y_vals)
        end do
      end do
    end do

    nullify(mat)
    nullify(nuc)
    nullify(tope)

  end subroutine write_xs


!> Write out ACE format elastic secondary angular distribution data
  subroutine write_angle()

    type(Material), pointer :: mat => null() ! material pointer
    type(Nuclide), pointer  :: nuc => null() ! nuclide pointer
    type(Reaction), pointer :: rxn => null() ! reaction pointer
    integer :: i_mat  ! material index
    integer :: i_nuc  ! nuclide index
    integer :: i_rxn  ! reaction index
    integer :: i_MT   ! MT index
    integer :: i_E    ! energy index
    integer :: i_mu   ! cosine index
    integer :: interp ! type of interpolation
    integer :: type   ! angular distribution type
    integer :: n      ! number of incoming energy bins
    integer :: lc     ! location in data array
    integer :: NP     ! number of points in cos distribution

    ! loop over all materials
    MAT_LOOP: do i_mat = 1, n_materials
      mat => materials(i_mat)

      if (mat % nat_elements) then
        call log_message(INFO, 'angular_dist element(s) in material(s) with &
             &natural elements will be ignored')
        cycle
      end if

      ! loop over all nuclides in material
      NUC_LOOP: do i_nuc = 1, mat % n_nuclides
        nuc => nuclides(mat % nuclide(i_nuc))

        ! loop over all requested reaction MT numbers
        RXN_LOOP: do i_rxn = 1, mat % write_angle(i_nuc) % size()

          ! determine which reaction is requested
          select case(trim(adjustl(mat % write_angle(i_nuc) % get_item(i_rxn))))
          case('2')
            ! get elastic scattering reaction
            i_MT = 1
            do
              rxn => nuc % reactions(i_MT)
              if (rxn % MT == 2) exit
              i_MT = i_MT + 1
            end do

          case default
            call fatal_error('Writing secondary angular distribution for this&
                 & reaction is not supported.')

          end select

          ! check if reaction has angular distribution
          if (.not. rxn % has_angle_dist)&
               call fatal_error('No sec. angular dist. exists for writing')

          ! determine number of incoming energies
          n = rxn % adist % n_energy

          open(unit=90, file='sec-ang-dist-MT'//trim(adjustl(to_str(rxn%MT)))&
               //'-'//trim(adjustl(to_str(nuc%zaid)))//'.dat')
          
          ENERGY_LOOP: do i_E = 2, n
            write(90, '(ES24.16)') rxn % adist % energy(i_E)
            lc  = rxn % adist % location(i_E)
            type = rxn % adist % type(i_E)
            
            if (type == ISOTROPIC) then
              call fatal_error('isotropic secondary angular distribution')
            
            elseif (type == EQUIPROBABLE_32_BINS) then
              do i_mu = 1, 33
                write(90, '(ES24.16)') rxn % adist % energy(i_E),&
                                       rxn % adist % data(lc + i_mu)
              end do

            elseif (type == TABULAR) then
              interp = int(rxn % adist % data(lc + 1))
              NP     = int(rxn % adist % data(lc + 2))
              lc = lc + 2
              do i_mu = 1, NP
                write(90, '(ES24.16,ES24.16,ES24.16,ES24.16)')&
                     rxn % adist % energy(i_E),&
                     rxn % adist % data(lc + 0*NP + i_mu),&
                     rxn % adist % data(lc + 1*NP + i_mu),&
                     rxn % adist % data(lc + 2*NP + i_mu)
              
              end do
            end if
          end do ENERGY_LOOP

          close(90)
        
        end do RXN_LOOP
      end do NUC_LOOP
    end do MAT_LOOP

  end subroutine write_angle
  

!> Generate realizations of level spacings by sampling the Wigner's Surmise
!! distribution
  subroutine wigner_dist_samples(n)

    integer :: n ! number of samples
    integer :: i ! iteration index
    real(8) :: D_mean    ! average level spacing in eV
    real(8) :: D_vals(n) ! sampled level spacings

    D_mean = 20.0_8

    do i = 1, n
      D_vals(i) = wigner_level_spacing(D_mean, prn())
    end do

    call write_coords(99, 'wigner-dist-samples.dat',&
         dble([(i, i = 1, n)]), D_vals(:) / D_mean)

  end subroutine wigner_dist_samples


end module URR_edits
