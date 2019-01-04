module physics

  use constants
  use error,                  only: fatal_error
  use material_header,        only: Material, materials
  use math
  use message_passing
  use nuclide_header
  use particle_header
  use photon_header
  use random_lcg,             only: prn
  use settings
  use simulation_header

  implicit none

  interface
    subroutine collision(p) bind(C)
      import Particle
      type(Particle), intent(inout) :: p
    end subroutine
  end interface

contains

!===============================================================================
! SAMPLE_ELEMENT
!===============================================================================

  function sample_element(p) result(i_element) bind(C)
    type(Particle), intent(in) :: p
    integer(C_INT) :: i_element

    integer :: i
    real(8) :: prob
    real(8) :: cutoff
    real(8) :: atom_density ! atom density of nuclide in atom/b-cm
    real(8) :: sigma        ! microscopic total xs for nuclide

    associate (mat => materials(p % material))
      ! Sample cumulative distribution function
      cutoff = prn() * material_xs % total

      i = 0
      prob = ZERO
      do while (prob < cutoff)
        i = i + 1

        ! Check to make sure that a nuclide was sampled
        if (i > mat % n_nuclides) then
          call particle_write_restart(p)
          call fatal_error("Did not sample any element during collision.")
        end if

        ! Find atom density
        i_element    = mat % element(i)
        atom_density = mat % atom_density(i)

        ! Determine microscopic cross section
        sigma = atom_density * micro_photon_xs(i_element) % total

        ! Increment probability to compare to cutoff
        prob = prob + sigma
      end do
    end associate

  end function sample_element

end module physics
