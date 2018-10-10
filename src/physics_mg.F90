module physics_mg
  ! This module contains the multi-group specific physics routines so as to not
  ! hinder performance of the CE versions with multiple if-thens.

  use bank_header
  use constants
  use error,                  only: fatal_error, warning, write_message
  use material_header,        only: Material, materials
  use math,                   only: rotate_angle
  use mgxs_interface
  use message_passing
  use nuclide_header,         only: MaterialMacroXS, material_xs
  use particle_header
  use physics_common
  use random_lcg,             only: prn
  use settings
  use simulation_header
  use string,                 only: to_str
  use tally_header

  implicit none

  interface
    subroutine scatter(p, energy_bin_avg) bind(C)
      import Particle, C_DOUBLE
      type(Particle), intent(inout) :: p
      real(C_DOUBLE), intent(in)    :: energy_bin_avg(*)
    end subroutine scatter

    subroutine create_fission_sites(p, bank_array, size_bank, bank_array_size, &
         material_xs) bind(C)
      import Particle, Bank, C_INT64_T, MaterialMacroXS
      type(Particle),        intent(inout) :: p
      type(Bank),            intent(inout) :: bank_array(*)
      integer(C_INT64_T),    intent(inout) :: size_bank
      integer(C_INT64_T),    intent(in)    :: bank_array_size
      type(MaterialMacroXS), intent(in)    :: material_xs
    end subroutine create_fission_sites

    subroutine absorption(p, material_xs) bind(C)
      import Particle, MaterialMacroXS
      type(Particle),        intent(inout) :: p
      type(MaterialMacroXS), intent(in)    :: material_xs
    end subroutine absorption
  end interface

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
        call create_fission_sites(p, fission_bank, n_bank, &
                                  size(fission_bank, KIND=C_INT64_T), &
                                  material_xs)
      elseif (run_mode == MODE_FIXEDSOURCE .and. create_fission_neutrons) then
        call create_fission_sites(p, p % secondary_bank, p % n_secondary, &
                                  size(p % secondary_bank, KIND=C_INT64_T), &
                                  material_xs)
      end if
    end if

    ! If survival biasing is being used, the following subroutine adjusts the
    ! weight of the particle. Otherwise, it checks to see if absorption occurs

    if (material_xs % absorption > ZERO) then
      call absorption(p, material_xs)
    else
      p % absorb_wgt = ZERO
    end if
    if (.not. p % alive) return

    ! Sample a scattering reaction and determine the secondary energy of the
    ! exiting neutron
    call scatter(p, energy_bin_avg)

    ! Play russian roulette if survival biasing is turned on
    if (survival_biasing) then
      call russian_roulette(p)
      if (.not. p % alive) return
    end if

  end subroutine sample_reaction

end module physics_mg
