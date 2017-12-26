module physics_mg
  ! This module contains the multi-group specific physics routines so as to not
  ! hinder performance of the CE versions with multiple if-thens.

  use bank_header
  use constants
  use error,                  only: fatal_error, warning, write_message
  use material_header,        only: Material, materials
  use math,                   only: rotate_angle
  use mesh_header,            only: meshes
  use mgxs_header
  use message_passing
  use nuclide_header,         only: material_xs
  use particle_header,        only: Particle
  use physics_common
  use random_lcg,             only: prn
  use scattdata_header
  use settings
  use simulation_header
  use string,                 only: to_str
  use tally_header

  implicit none

contains

!===============================================================================
! COLLISION_MG samples a nuclide and reaction and then calls the appropriate
! routine for that reaction
!===============================================================================

  subroutine collision_mg(p)

    type(Particle), intent(inout) :: p

    ! Add to collision counter for particle
    p % n_collision = p % n_collision + 1

    ! Sample nuclide/reaction for the material the particle is in
    call sample_reaction(p)

    ! Display information about collision
    if (verbosity >= 10 .or. trace) then
      call write_message("    " // "Energy Group = " // trim(to_str(p % g)))
    end if

  end subroutine collision_mg

!===============================================================================
! SAMPLE_REACTION samples a nuclide based on the macroscopic cross sections for
! each nuclide within a material and then samples a reaction for that nuclide
! and calls the appropriate routine to process the physics. Note that there is
! special logic when suvival biasing is turned on since fission and
! disappearance are treated implicitly.
!===============================================================================

  subroutine sample_reaction(p)

    type(Particle), intent(inout) :: p

    type(Material), pointer :: mat

    mat => materials(p % material)

    ! Create fission bank sites. Note that while a fission reaction is sampled,
    ! it never actually "happens", i.e. the weight of the particle does not
    ! change when sampling fission sites. The following block handles all
    ! absorption (including fission)

    if (mat % fissionable) then
      if (run_mode == MODE_EIGENVALUE) then
        call create_fission_sites(p, fission_bank, n_bank)
      elseif (run_mode == MODE_FIXEDSOURCE .and. create_fission_neutrons) then
        call create_fission_sites(p, p % secondary_bank, p % n_secondary)
      end if
    end if

    ! If survival biasing is being used, the following subroutine adjusts the
    ! weight of the particle. Otherwise, it checks to see if absorption occurs

    if (material_xs % absorption > ZERO) then
      call absorption(p)
    else
      p % absorb_wgt = ZERO
    end if
    if (.not. p % alive) return

    ! Sample a scattering reaction and determine the secondary energy of the
    ! exiting neutron
    call scatter(p)

    ! Play russian roulette if survival biasing is turned on
    if (survival_biasing) then
      call russian_roulette(p)
      if (.not. p % alive) return
    end if

  end subroutine sample_reaction

!===============================================================================
! ABSORPTION
!===============================================================================

  subroutine absorption(p)

    type(Particle), intent(inout) :: p

    if (survival_biasing) then
      ! Determine weight absorbed in survival biasing
      p % absorb_wgt = (p % wgt * &
                        material_xs % absorption / material_xs % total)

      ! Adjust weight of particle by probability of absorption
      p % wgt = p % wgt - p % absorb_wgt
      p % last_wgt = p % wgt

      ! Score implicit absorption estimate of keff
!$omp atomic
      global_tallies(RESULT_VALUE, K_ABSORPTION) = &
           global_tallies(RESULT_VALUE, K_ABSORPTION) + p % absorb_wgt * &
           material_xs % nu_fission / material_xs % absorption
    else
      ! See if disappearance reaction happens
      if (material_xs % absorption > prn() * material_xs % total) then
        ! Score absorption estimate of keff
!$omp atomic
        global_tallies(RESULT_VALUE, K_ABSORPTION) = &
             global_tallies(RESULT_VALUE, K_ABSORPTION) + p % wgt * &
             material_xs % nu_fission / material_xs % absorption

        p % alive = .false.
        p % event = EVENT_ABSORB
      end if
    end if

  end subroutine absorption

!===============================================================================
! SCATTER
!===============================================================================

  subroutine scatter(p)

    type(Particle), intent(inout)  :: p

    call macro_xs(p % material) % obj % sample_scatter(p % coord(1) % uvw, &
                                                       p % last_g, p % g, &
                                                       p % mu, p % wgt)

    ! Update energy value for downstream compatability (in tallying)
    p % E = energy_bin_avg(p % g)

    ! Convert change in angle (mu) to new direction
    p % coord(1) % uvw = rotate_angle(p % coord(1) % uvw, p % mu)

    ! Set event component
    p % event = EVENT_SCATTER

  end subroutine scatter

!===============================================================================
! CREATE_FISSION_SITES determines the average total, prompt, and delayed
! neutrons produced from fission and creates appropriate bank sites.
!===============================================================================

  subroutine create_fission_sites(p, bank_array, size_bank)
    type(Particle), intent(inout) :: p
    type(Bank),     intent(inout) :: bank_array(:)
    integer(8),     intent(inout) :: size_bank

    integer :: nu_d(MAX_DELAYED_GROUPS) ! number of delayed neutrons born
    integer :: i                        ! loop index
    integer :: dg                       ! delayed group
    integer :: gout                     ! group out
    integer :: nu                       ! actual number of neutrons produced
    integer :: mesh_bin                 ! mesh bin for source site
    real(8) :: nu_t                     ! total nu
    real(8) :: mu                       ! fission neutron angular cosine
    real(8) :: phi                      ! fission neutron azimuthal angle
    real(8) :: weight                   ! weight adjustment for ufs method
    class(Mgxs), pointer :: xs

    ! Get Pointers
    xs => macro_xs(p % material) % obj

    ! TODO: Heat generation from fission

    ! If uniform fission source weighting is turned on, we increase of decrease
    ! the expected number of fission sites produced

    if (ufs) then
      associate (m => meshes(index_ufs_mesh))
        ! Determine indices on ufs mesh for current location
        call m % get_bin(p % coord(1) % xyz, mesh_bin)

        if (mesh_bin == NO_BIN_FOUND) then
          call p % write_restart()
          call fatal_error("Source site outside UFS mesh!")
        end if

        if (source_frac(1, mesh_bin) /= ZERO) then
          weight = m % volume_frac / source_frac(1, mesh_bin)
        else
          weight = ONE
        end if
      end associate
    else
      weight = ONE
    end if

    ! Determine expected number of neutrons produced
    nu_t = p % wgt / keff * weight * &
         material_xs % nu_fission / material_xs % total

    ! Sample number of neutrons produced
    if (prn() > nu_t - int(nu_t)) then
      nu = int(nu_t)
    else
      nu = int(nu_t) + 1
    end if

    ! Check for bank size getting hit. For fixed source calculations, this is a
    ! fatal error. For eigenvalue calculations, it just means that k-effective
    ! was too high for a single batch.
    if (size_bank + nu > size(bank_array)) then
      if (run_mode == MODE_FIXEDSOURCE) then
        call fatal_error("Secondary particle bank size limit reached. If you &
             &are running a subcritical multiplication problem, k-effective &
             &may be too close to one.")
      else
        if (master) call warning("Maximum number of sites in fission bank &
             &reached. This can result in irreproducible results using different &
             &numbers of processes/threads.")
      end if
    end if

    ! Bank source neutrons
    if (nu == 0 .or. size_bank == size(bank_array)) return

    ! Initialize counter of delayed neutrons encountered for each delayed group
    ! to zero.
    nu_d(:) = 0

    p % fission = .true. ! Fission neutrons will be banked
    do i = int(size_bank,4) + 1, int(min(size_bank + nu, int(size(bank_array),8)),4)
      ! Bank source neutrons by copying particle data
      bank_array(i) % xyz = p % coord(1) % xyz

      ! Set weight of fission bank site
      bank_array(i) % wgt = ONE/weight

      ! Sample cosine of angle -- fission neutrons are treated as being emitted
      ! isotropically.
      mu = TWO * prn() - ONE

      ! Sample azimuthal angle uniformly in [0,2*pi)
      phi = TWO * PI * prn()
      bank_array(i) % uvw(1) = mu
      bank_array(i) % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
      bank_array(i) % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

      ! Sample secondary energy distribution for fission reaction and set energy
      ! in fission bank
      call xs % sample_fission_energy(p % g, bank_array(i) % uvw, dg, gout)

      bank_array(i) % E = real(gout, 8)
      bank_array(i) % delayed_group = dg

      ! Set delayed group on particle too
      p % delayed_group = dg

      ! Increment the number of neutrons born delayed
      if (p % delayed_group > 0) then
        nu_d(p % delayed_group) = nu_d(p % delayed_group) + 1
      end if

    end do

    ! increment number of bank sites
    size_bank = min(size_bank + nu, int(size(bank_array),8))

    ! Store total weight banked for analog fission tallies
    p % n_bank   = nu
    p % wgt_bank = nu/weight
    p % n_delayed_bank(:) = nu_d(:)

  end subroutine create_fission_sites

end module physics_mg
