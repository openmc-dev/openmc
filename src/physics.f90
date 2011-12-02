module physics

  use constants
  use cross_section_header, only: Nuclide, Reaction, DistEnergy
  use endf,                 only: reaction_name, is_fission, is_scatter
  use error,                only: fatal_error, warning
  use fission,              only: nu_total, nu_prompt, nu_delayed
  use geometry,             only: find_cell, distance_to_boundary,              &
                                  cross_surface, cross_lattice
  use geometry_header,      only: Universe, BASE_UNIVERSE
  use global
  use interpolation,        only: interpolate_tab1
  use output,               only: write_message
  use particle_header,      only: Particle, LocalCoord
  use random_lcg,           only: prn
  use search,               only: binary_search
  use string,               only: int_to_str
  use tally,                only: score_tally, score_surface_current

  implicit none

contains

!===============================================================================
! TRANSPORT encompasses the main logic for moving a particle through geometry.
!===============================================================================

  subroutine transport(p)

    type(Particle), pointer :: p

    integer        :: surface_crossed ! surface which particle is on
    integer        :: last_cell       ! most recent cell particle was in
    integer        :: n_event         ! number of collisions/crossings
    real(8)        :: d_boundary      ! distance to nearest boundary
    real(8)        :: d_collision     ! sampled distance to collision
    real(8)        :: distance        ! distance particle travels
    logical        :: found_cell      ! found cell which particle is in?
    logical        :: lattice_crossed ! is surface crossing in lattice?
    type(LocalCoord), pointer :: coord

    if (p % coord % cell == NONE) then
       call find_cell(p, found_cell)
       ! Particle couldn't be located
       if (.not. found_cell) then
          write(message, '(A,3ES11.3)') & 
               "Could not locate cell for particle at: ", p % coord0 % xyz
          call fatal_error()
       end if

       ! set birth cell attribute
       p % cell_born = p % coord % cell
    end if

    if (verbosity >= 9) then
       message = "Simulating Particle " // trim(int_to_str(p % id))
       call write_message()
    end if

    if (verbosity >= 10) then
       message = "    Born in cell " // trim(int_to_str(&
            cells(p % coord % cell) % id))
       call write_message()
    end if

    ! Initialize number of events to zero
    n_event = 0

    ! find energy index, interpolation factor
    do while (p % alive)

       ! Calculate microscopic and macroscopic cross sections
       call calculate_xs(p)

       ! Find the distance to the nearest boundary
       call distance_to_boundary(p, d_boundary, surface_crossed, lattice_crossed)

       ! Sample a distance to collision
       d_collision = -log(prn()) / material_xs % total
       
       ! Select smaller of the two distances
       distance = min(d_boundary, d_collision)

       ! Advance particle
       coord => p % coord0
       do while (associated(coord))
          coord % xyz = coord % xyz + distance * coord % uvw
          coord => coord % next
       end do

       if (d_collision > d_boundary) then
          last_cell = p % coord % cell
          p % coord % cell = NONE
          if (lattice_crossed) then
             p % surface = NONE
             call cross_lattice(p)
          else
             p % surface = surface_crossed
             call cross_surface(p, last_cell)
          end if
       else
          ! collision
          p % surface = NONE
          call collision(p)

          ! Save coordinates at collision for tallying purposes
          p % last_xyz = p % coord0 % xyz

          ! Set all uvws to base level -- right now, after a collision, only the
          ! base level uvws are changed
          coord => p % coord0
          do while(associated(coord))
             coord % uvw = p % coord0 % uvw
             coord => coord % next
          end do
       end if

       ! If particle has too many events, display warning and kill it
       n_event = n_event + 1
       if (n_event == MAX_EVENTS) then
          message = "Particle " // trim(int_to_str(p%id)) // " underwent " & 
               // "maximum number of events."
          call warning()
          p % alive = .false.
       end if

    end do

  end subroutine transport

!===============================================================================
! CALCULATE_XS determines the macroscopic cross sections for the material the
! particle is currently traveling through.
!===============================================================================

  subroutine calculate_xs(p)

    type(Particle), pointer :: p

    integer                  :: i             ! loop index over nuclides
    integer                  :: index_nuclide ! index into nuclides array
    integer                  :: index_sab
    real(8)                  :: atom_density  ! atom density of a nuclide
    real(8)                  :: sab_threshold ! threshold for S(a,b) table
    type(Material),  pointer :: mat => null() ! current material

    ! If the material is the same as the last material and the energy of the
    ! particle hasn't changed, we don't need to lookup cross sections again.

    if (p % material == p % last_material) return

    ! Set all material macroscopic cross sections to zero
    material_xs % total      = ZERO
    material_xs % elastic    = ZERO
    material_xs % absorption = ZERO
    material_xs % fission    = ZERO
    material_xs % nu_fission = ZERO

    mat => materials(p % material)

    ! Find energy index on unionized grid
    call find_energy_index(p)

    ! Check if there's an S(a,b) table for this material
    if (mat % has_sab_table) then
       sab_threshold = sab_tables(mat % sab_table) % threshold_inelastic
    else
       sab_threshold = ZERO
    end if

    ! Add contribution from each nuclide in material
    do i = 1, mat % n_nuclides
       ! Determine microscopic cross sections for this nuclide
       index_nuclide = mat % nuclide(i)

       ! Determine whether to use S(a,b) based on energy of particle and whether
       ! the nuclide matches the S(a,b) table
       if (p % E < sab_threshold .and. i == mat % sab_nuclide) then
          index_sab = mat % sab_table
       else
          index_sab = 0
       end if

       ! Calculate microscopic cross section for this nuclide
       call calculate_nuclide_xs(p, index_nuclide, index_sab)

       ! Copy atom density of nuclide in material
       atom_density = mat % atom_density(i)

       ! Add contributions to material macroscopic total cross section
       material_xs % total = material_xs % total + &
            atom_density * micro_xs(index_nuclide) % total
       
       ! Add contributions to material macroscopic scattering cross section
       material_xs % elastic = material_xs % elastic + &
            atom_density * micro_xs(index_nuclide) % elastic
       
       ! Add contributions to material macroscopic absorption cross section
       material_xs % absorption = material_xs % absorption + & 
            atom_density * micro_xs(index_nuclide) % absorption
       
       ! Add contributions to material macroscopic fission cross section
       material_xs % fission = material_xs % fission + &
            atom_density * micro_xs(index_nuclide) % fission
       
       ! Add contributions to material macroscopic nu-fission cross section
       material_xs % nu_fission = material_xs % nu_fission + &
            atom_density * micro_xs(index_nuclide) % nu_fission
    end do

  end subroutine calculate_xs

!===============================================================================
! CALCULATE_NUCLIDE_XS determines microscopic cross sections for a nuclide of a
! given index in the nuclides array at the energy of the given particle
!===============================================================================

  subroutine calculate_nuclide_xs(p, index_nuclide, index_sab)

    type(Particle), pointer :: p
    integer, intent(in)     :: index_nuclide ! index into nuclides array
    integer, intent(in)     :: index_sab     ! index into sab_tables array

    integer                :: i         ! index into nuclides array
    integer                :: IE        ! index on nuclide energy grid
    integer                :: IE_sab    ! index on S(a,b) energy grid
    real(8)                :: f         ! interp factor on nuclide energy grid
    real(8)                :: f_sab     ! interp factor on S(a,b) energy grid
    real(8)                :: inelastic ! S(a,b) inelastic cross section
    real(8)                :: elastic   ! S(a,b) elastic cross section
    real(8)                :: nu        ! total # of neutrons emitted per fission
    type(Nuclide),   pointer :: nuc => null()
    type(SAB_Table), pointer :: sab => null()

    ! Copy index of nuclide
    i = index_nuclide

    ! Set pointer to nuclide
    nuc => nuclides(i)

    ! TODO: Check if last energy/temp combination is same as current. If so, we
    ! can return.

    ! TODO: If not using unionized energy grid, we need to find the index on the
    ! nuclide energy grid using lethargy mapping or whatever other technique

    ! search nuclide energy grid
    IE = nuc % grid_index(p % IE)
    f = (p%E - nuc%energy(IE))/(nuc%energy(IE+1) - nuc%energy(IE))

    micro_xs(i) % index_grid = IE
    micro_xs(i) % interp_factor = f

    ! Initialize sab treatment to false
    micro_xs(i) % use_sab     = .false.
    micro_xs(i) % elastic_sab = ZERO

    ! Initialize nuclide cross-sections to zero
    micro_xs(i) % fission    = ZERO
    micro_xs(i) % nu_fission = ZERO

    ! Calculate microscopic nuclide total cross section
    micro_xs(i) % total = &
         (ONE-f) * nuc % total(IE) + f * nuc % total(IE+1)

    ! Calculate microscopic nuclide total cross section
    micro_xs(i) % elastic = &
         (ONE-f) * nuc % elastic(IE) + f * nuc % elastic(IE+1)

    ! Calculate microscopic nuclide absorption cross section
    micro_xs(i) % absorption = &
         (ONE-f) * nuc % absorption(IE) + f * nuc % absorption(IE+1)

    if (nuc % fissionable) then
       ! Calculate microscopic nuclide total cross section
       micro_xs(i) % fission = &
            (ONE-f) * nuc % fission(IE) + f * nuc % fission(IE+1)

       ! Calculate microscopic nuclide nu-fission cross section
       nu = nu_total(nuc, p % E)
       micro_xs(i) % nu_fission = nu * micro_xs(i) % fission
    end if

    ! If there is S(a,b) data for this nuclide, we need to do a few
    ! things. Since the total cross section was based on non-S(a,b) data, we
    ! need to correct it by subtracting the non-S(a,b) elastic cross section and
    ! then add back in the calculated S(a,b) elastic+inelastic cross section.

    if (index_sab > 0) then
       micro_xs(i) % use_sab = .true.

       ! Get pointer to S(a,b) table
       sab => sab_tables(index_sab)

       ! Get index and interpolation factor for inelastic grid
       if (p%E < sab % inelastic_e_in(1)) then
          IE_sab = 1
          f_sab = ZERO
       else
          IE_sab = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, p%E)
          f_sab = (p%E - sab%inelastic_e_in(IE_sab)) / & 
               (sab%inelastic_e_in(IE_sab+1) - sab%inelastic_e_in(IE_sab))
       end if

       ! Calculate S(a,b) inelastic scattering cross section
       inelastic = (ONE-f_sab) * sab % inelastic_sigma(IE_sab) + f_sab * &
            sab % inelastic_sigma(IE_sab + 1)

       ! Check for elastic data
       if (p % E < sab % threshold_elastic) then
          ! Determine whether elastic scattering is given in the coherent or
          ! incoherent approximation. For coherent, the cross section is
          ! represented as P/E whereas for incoherent, it is simply P

          if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
             if (p % E < sab % elastic_e_in(1)) then
                ! If energy is below that of the lowest Bragg peak, the elastic
                ! cross section will be zero
                elastic = ZERO
             else
                IE_sab = binary_search(sab % elastic_e_in, sab % n_elastic_e_in, p%E)
                elastic = sab % elastic_P(IE_sab) / p % E
             end if
          else
             ! Determine index on elastic energy grid
             if (p % E < sab % elastic_e_in(1)) then
                IE_sab = 1
             else
                IE_sab = binary_search(sab % elastic_e_in, sab % n_elastic_e_in, p%E)
             end if

             ! Get interpolation factor for elastic grid
             f_sab = (p%E - sab%elastic_e_in(IE_sab))/(sab%elastic_e_in(IE_sab+1) - &
                  sab%elastic_e_in(IE_sab))

             ! Calculate S(a,b) elastic scattering cross section
             elastic = (ONE-f_sab) * sab % elastic_P(IE_sab) + f_sab * &
                  sab % elastic_P(IE_sab + 1)
          end if
       else
          ! No elastic data
          elastic = ZERO
       end if

       ! Correct total and elastic cross sections
       micro_xs(i) % total = micro_xs(i) % total - micro_xs(i) % elastic &
            + inelastic + elastic
       micro_xs(i) % elastic = inelastic + elastic

       ! Store S(a,b) elastic cross section for sampling later
       micro_xs(i) % elastic_sab = elastic
    end if

  end subroutine calculate_nuclide_xs

!===============================================================================
! FIND_ENERGY_INDEX determines the index on the union energy grid and the
! interpolation factor for a particle at a certain energy
!===============================================================================

  subroutine find_energy_index(p)

    type(Particle), pointer :: p

    integer :: IE     ! index on union energy grid
    real(8) :: E      ! energy of particle
    real(8) :: interp ! interpolation factor

    ! copy particle's energy
    E = p % E

    ! if particle's energy is outside of energy grid range, set to first or last
    ! index. Otherwise, do a binary search through the union energy grid.
    if (E < e_grid(1)) then
       IE = 1
    elseif (E > e_grid(n_grid)) then
       IE = n_grid - 1
    else
       IE = binary_search(e_grid, n_grid, E)
    end if
    
    ! calculate the interpolation factor -- note this will be outside of [0,1)
    ! for a particle outside the energy range of the union grid
    interp = (E - e_grid(IE))/(e_grid(IE+1) - e_grid(IE))

    ! set particle attributes
    p % IE     = IE
    p % interp = interp
    
  end subroutine find_energy_index

!===============================================================================
! COLLISION samples a nuclide and reaction and then calls the appropriate
! routine for that reaction
!===============================================================================

  subroutine collision(p)

    type(Particle), pointer :: p

    integer :: MT            ! ENDF reaction number
    logical :: scattered     ! was this a scattering reaction?
    logical :: fissioned     ! was this a fission reaction?

    ! Store pre-collision particle properties
    p % last_wgt = p % wgt
    p % last_E   = p % E

    ! Add to collision counter for particle
    p % n_collision = p % n_collision + 1

    ! Sample nuclide/reaction for the material the particle is in
    call sample_reaction(p, MT)

    ! check for very low energy
    if (p % E < 1.0e-100_8) then
       p % alive = .false.
       message = "Killing neutron with extremely low energy"
       call warning()
    end if

    ! Score collision estimator tallies for any macro tallies -- this is done
    ! after a collision has occurred rather than before because we need
    ! information on the outgoing energy for any tallies with an outgoing energy
    ! filter

    ! Check if particle scattered or fissioned
    if (survival_biasing) then
       fissioned = .false.
       scattered = .true.
    else
       fissioned = is_fission(MT)
       scattered = is_scatter(MT)
    end if

    if (tallies_on) then
       call score_tally(p, scattered, fissioned)
       call score_surface_current(p)
    end if

    ! Reset number of particles banked during collision
    p % n_bank = 0

    ! find energy index, interpolation factor
    call find_energy_index(p)

  end subroutine collision

!===============================================================================
! SAMPLE_REACTION samples a nuclide based on the macroscopic cross sections for
! each nuclide within a material and then samples a reaction for that nuclide
! and calls the appropriate routine to process the physics. Note that there is
! special logic when suvival biasing is turned on since fission and
! disappearance are treated implicitly.
!===============================================================================

  subroutine sample_reaction(p, MT)

    type(Particle), pointer :: p
    integer, intent(out)    :: MT ! ENDF MT number of reaction that occured

    integer :: i             ! index over nuclides in a material
    integer :: index_nuclide ! index in nuclides array
    integer :: IE            ! index on nuclide energy grid
    real(8) :: f             ! interpolation factor
    real(8) :: sigma         ! microscopic total xs for nuclide
    real(8) :: prob          ! cumulative probability
    real(8) :: cutoff        ! random number
    real(8) :: atom_density  ! atom density of nuclide in atom/b-cm
    type(Material), pointer :: mat => null()
    type(Nuclide),  pointer :: nuc => null()
    type(Reaction), pointer :: rxn => null()

    ! Get pointer to current material
    mat => materials(p % material)

    ! ==========================================================================
    ! SAMPLE NUCLIDE WITHIN THE MATERIAL

    cutoff = prn() * material_xs % total
    prob = ZERO

    i = 0
    do while (prob < cutoff)
       i = i + 1

       ! Check to make sure that a nuclide was sampled
       if (i > mat % n_nuclides) then
          message = "Did not sample any nuclide during collision."
          call fatal_error()
       end if

       ! Find atom density and microscopic total cross section
       index_nuclide = mat % nuclide(i)
       atom_density = mat % atom_density(i)
       sigma = atom_density * micro_xs(index_nuclide) % total

       ! Increment probability to compare to cutoff
       prob = prob + sigma
    end do

    ! Get pointer to table, nuclide grid index and interpolation factor
    nuc => nuclides(index_nuclide)
    IE  =  micro_xs(index_nuclide) % index_grid
    f   =  micro_xs(index_nuclide) % interp_factor

    ! ==========================================================================
    ! DISAPPEARANCE REACTIONS (ANALOG) OR IMPLICIT CAPTURE (SURVIVAL BIASING)

    if (survival_biasing) then
       ! adjust weight of particle by probability of absorption
       p % wgt = p % wgt * (ONE - micro_xs(index_nuclide) % absorption / &
            micro_xs(index_nuclide) % total)

    else
       ! set cutoff variable for analog cases
       cutoff = prn() * micro_xs(index_nuclide) % total
       prob = ZERO

       ! Add disappearance cross-section to prob
       prob = prob + micro_xs(index_nuclide) % absorption - &
            micro_xs(index_nuclide) % fission

       ! See if disappearance reaction happens
       if (prob > cutoff) then
          p % alive = .false.
          MT = N_DISAPPEAR
          return
       end if
    end if

    ! ==========================================================================
    ! FISSION EVENTS (ANALOG) OR BANK EXPECTED FISSION SITES (IMPLICIT)

    if (nuc % fissionable) then
       if (survival_biasing) then
          ! If survival biasing is turned on, then no fission events actually
          ! occur since absorption is treated implicitly. However, we still need
          ! to bank sites so we sample a fission reaction (if there are
          ! multiple) and bank the expected number of fission neutrons
          ! created.

          if (nuc % has_partial_fission) then
             cutoff = prn() * micro_xs(index_nuclide) % fission
             prob = ZERO

             i = 0
             do while (prob < cutoff)
                i = i + 1

                ! Check to make sure partial fission reaction sampled
                if (i > nuc % n_fission) then
                   message = "Did not sample any partial fission reaction " // &
                        "for survival biasing in " // trim(nuc % name)
                   call fatal_error()
                end if

                rxn => nuc % reactions(nuc % index_fission(i))

                ! if energy is below threshold for this reaction, skip it
                if (IE < rxn%IE) cycle

                ! add to cumulative probability
                prob = prob + ((ONE-f)*rxn%sigma(IE-rxn%IE+1) & 
                     & + f*(rxn%sigma(IE-rxn%IE+2)))
             end do
          else
             ! For nuclides with only total fission reaction, get a pointer to
             ! the fission reaction
             rxn => nuc % reactions(nuc % index_fission(1))
          end if

          ! Bank expected number of fission neutrons
          call create_fission_sites(p, index_nuclide, rxn)

       else
          ! If survival biasing is not turned on, then fission events can occur
          ! just like any other reaction. Here we loop through the fission
          ! reactions for the nuclide and see if any of them occur
          
          do i = 1, nuc % n_fission
             rxn => nuc % reactions(nuc % index_fission(i))

             ! if energy is below threshold for this reaction, skip it
             if (IE < rxn%IE) cycle

             ! add to cumulative probability
             if (nuc % has_partial_fission) then
                prob = prob + ((ONE-f)*rxn%sigma(IE-rxn%IE+1) & 
                     & + f*(rxn%sigma(IE-rxn%IE+2)))
             else
                prob = prob + micro_xs(index_nuclide) % fission
             end if

             ! Create fission bank sites if fission occus
             if (prob > cutoff) then
                call create_fission_sites(p, index_nuclide, rxn, .true.)
                p % alive = .false.
                MT = rxn % MT
                return
             end if
          end do
       end if

    end if

    ! ==========================================================================
    ! WEIGHT CUTOFF (SURVIVAL BIASING ONLY)

    if (survival_biasing) then
       if (p % wgt < weight_cutoff) then
          if (prn() < p % wgt / weight_survive) then
             p % wgt = weight_survive
          else
             p % wgt = ZERO
             p % alive = .false.
          end if
       end if

       ! At this point, we also need to set the cutoff variable for cases with
       ! survival biasing. The cutoff will be a random number times the
       ! scattering cross section

       cutoff = prn() * (micro_xs(index_nuclide) % total - &
            micro_xs(index_nuclide) % absorption)
       prob = ZERO
    end if

    ! ==========================================================================
    ! SCATTERING REACTIONS

    prob = prob + micro_xs(index_nuclide) % elastic
    if (prob > cutoff) then
       ! =======================================================================
       ! ELASTIC SCATTERING

       if (micro_xs(index_nuclide) % use_sab) then
          ! S(a,b) scattering
          call sab_scatter(p, index_nuclide, mat % sab_table)

       else
          ! get pointer to elastic scattering reaction
          rxn => nuc % reactions(1)

          ! Perform collision physics for elastic scattering
          call elastic_scatter(p, nuc, rxn)

       end if

       ! Set MT to be returned
       MT = 2
    else
       ! =======================================================================
       ! INELASTIC SCATTERING

       ! note that indexing starts from 2 since nuc % reactions(1) is elastic
       ! scattering
       i = 1
       do while (prob < cutoff)
          i = i + 1

          ! Check to make sure inelastic scattering reaction sampled
          if (i > nuc % n_reaction) then
             message = "Did not sample any reaction for nuclide " // &
                  trim(nuc % name) // " on material " // &
                  trim(int_to_str(mat % id))
             call fatal_error()
          end if

          rxn => nuc % reactions(i)

          ! Skip fission reactions
          if (rxn%MT == N_FISSION .or. rxn%MT == N_F .or. rxn%MT == N_NF &
               .or. rxn%MT == N_2NF .or. rxn%MT == N_3NF) cycle
             
          ! some materials have gas production cross sections with MT > 200 that
          ! are duplicates. Also MT=4 is total level inelastic scattering which
          ! should be skipped
          if (rxn%MT >= 200 .or. rxn%MT == N_LEVEL) cycle
          
          ! if energy is below threshold for this reaction, skip it
          if (IE < rxn%IE) cycle

          ! add to cumulative probability
          prob = prob + ((ONE-f)*rxn%sigma(IE-rxn%IE+1) & 
               & + f*(rxn%sigma(IE-rxn%IE+2)))
       end do

       ! Perform collision physics for inelastics scattering
       call inelastic_scatter(p, nuc, rxn)

       ! Set MT to be returned
       MT = rxn % MT
    end if

  end subroutine sample_reaction

!===============================================================================
! ELASTIC_SCATTER treats the elastic scattering of a neutron with a
! target. Currently this assumes target-at-rest kinematics -- obviously will
! need to be fixed
!===============================================================================

  subroutine elastic_scatter(p, nuc, rxn)

    type(Particle), pointer :: p
    type(Nuclide),  pointer :: nuc
    type(Reaction), pointer :: rxn

    real(8) :: awr     ! atomic weight ratio of target
    real(8) :: mu      ! cosine of polar angle
    real(8) :: vel     ! magnitude of velocity
    real(8) :: v_n(3)  ! velocity of neutron
    real(8) :: v_cm(3) ! velocity of center-of-mass
    real(8) :: v_t(3)  ! velocity of target nucleus
    real(8) :: u       ! x-direction
    real(8) :: v       ! y-direction
    real(8) :: w       ! z-direction
    real(8) :: E       ! energy

    vel = sqrt(p % E)
    awr = nuc % awr

    ! Neutron velocity in LAB
    v_n = vel * p % coord0 % uvw

    ! Sample velocity of target nucleus
    call sample_target_velocity(p, nuc, v_t)

    ! Velocity of center-of-mass
    v_cm = (v_n + awr*v_t)/(awr + ONE)

    ! Transform to CM frame
    v_n = v_n - v_cm

    ! Find speed of neutron in CM
    vel = sqrt(dot_product(v_n, v_n))

    ! Sample scattering angle
    mu = sample_angle(rxn, p % E)

    ! Determine direction cosines in CM
    u = v_n(1)/vel
    v = v_n(2)/vel
    w = v_n(3)/vel

    ! Change direction cosines according to mu
    call rotate_angle(u, v, w, mu)

    ! Rotate neutron velocity vector to new angle -- note that the speed of the
    ! neutron in CM does not change in elastic scattering. However, the speed
    ! will change when we convert back to LAB
    v_n = vel * (/ u, v, w /)

    ! Transform back to LAB frame
    v_n = v_n + v_cm

    E = dot_product(v_n, v_n)
    vel = sqrt(E)

    ! Set energy and direction of particle in LAB frame
    p % E = E
    p % coord0 % uvw = v_n / vel

    ! Copy scattering cosine for tallies
    p % mu = mu

  end subroutine elastic_scatter

!===============================================================================
! SAB_SCATTER performs thermal scattering of a particle with a bound scatterer
! according to a specified S(a,b) table.
!===============================================================================

  subroutine sab_scatter(p, index_nuclide, index_sab)

    type(Particle), pointer :: p
    integer, intent(in)     :: index_nuclide ! index in micro_xs
    integer, intent(in)     :: index_sab     ! index in sab_tables

    integer :: i            ! incoming energy bin
    integer :: j            ! outgoing energy bin
    integer :: k            ! outgoing cosine bin
    integer :: n_energy_out ! number of outgoing energy bins
    real(8) :: f            ! interpolation factor
    real(8) :: r            ! used for skewed sampling
    real(8) :: E            ! outgoing energy
    real(8) :: E_ij         ! outgoing energy j for E_in(i)
    real(8) :: E_i1j        ! outgoing energy j for E_in(i+1)
    real(8) :: mu           ! outgoing cosine
    real(8) :: mu_ijk       ! outgoing cosine k for E_in(i) and E_out(j)
    real(8) :: mu_i1jk      ! outgoing cosine k for E_in(i+1) and E_out(j)
    real(8) :: prob         ! probability for sampling Bragg edge
    real(8) :: u, v, w      ! directional cosines
    type(SAB_Table), pointer :: sab => null()

    ! Get pointer to S(a,b) table
    sab => sab_tables(index_sab)

    ! Determine whether inelastic or elastic scattering will occur
    if (prn() < micro_xs(index_nuclide) % elastic_sab / &
         micro_xs(index_nuclide) % elastic) then
       ! elastic scattering

       ! Get index and interpolation factor for elastic grid
       if (p%E < sab % elastic_e_in(1)) then
          i = 1
          f = ZERO
       else
          i = binary_search(sab % elastic_e_in, sab % n_elastic_e_in, p%E)
          f = (p%E - sab%elastic_e_in(i)) / & 
               (sab%elastic_e_in(i+1) - sab%elastic_e_in(i))
       end if

       ! Select treatment based on elastic mode
       if (sab % elastic_mode == SAB_ELASTIC_DISCRETE) then
          ! With this treatment, we interpolate between two discrete cosines
          ! corresponding to neighboring incoming energies. This is used for
          ! data derived in the incoherent approximation

          ! Sample outgoing cosine bin
          k = 1 + int(prn() * sab % n_elastic_mu)

          ! Determine outgoing cosine corresponding to E_in(i) and E_in(i+1)
          mu_ijk  = sab % elastic_mu(k,i)
          mu_i1jk = sab % elastic_mu(k,i+1)

          ! Cosine of angle between incoming and outgoing neutron
          mu = (1 - f)*mu_ijk + f*mu_i1jk

       elseif (sab % elastic_mode == SAB_ELASTIC_EXACT) then
          ! This treatment is used for data derived in the coherent
          ! approximation, i.e. for crystalline structures that have Bragg
          ! edges.
          
          ! Sample a Bragg edge between 1 and i
          prob = prn() * sab % elastic_P(i+1)
          k = binary_search(sab % elastic_P(1:i+1), i+1, prob)

          ! Characteristic scattering cosine for this Bragg egg
          mu = ONE - 2.0*sab % elastic_e_in(k) / p % E

       end if

       ! Outgoing energy is same as incoming energy
       E = p % E

    else
       ! Determine number of outgoing energy and angle bins
       n_energy_out = sab % n_inelastic_e_out
       
       ! Get index and interpolation factor for inelastic grid
       if (p%E < sab % inelastic_e_in(1)) then
          i = 1
          f = ZERO
       else
          i = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, p%E)
          f = (p%E - sab%inelastic_e_in(i)) / & 
               (sab%inelastic_e_in(i+1) - sab%inelastic_e_in(i))
       end if

       ! Now that we have an incoming energy bin, we need to determine the
       ! outgoing energy bin. This will depend on the "secondary energy
       ! mode". If the mode is 0, then the outgoing energy bin is chosen from a
       ! set of equally-likely bins. However, if the mode is 1, then the first
       ! two and last two bins are skewed to have lower probabilities than the
       ! other bins (0.1 for the first and last bins and 0.4 for the second and
       ! second to last bins, relative to a normal bin probability of 1)

       if (sab % secondary_mode == SAB_SECONDARY_EQUAL) then
          ! All bins equally likely
          j = 1 + int(prn() * n_energy_out)
       elseif (sab % secondary_mode == SAB_SECONDARY_SKEWED) then
          r = prn() * (n_energy_out - 3)
          if (r > ONE) then
             ! equally likely N-4 middle bins
             j = int(r) + 2
          elseif (r > 0.6) then
             ! second to last bin has relative probability of 0.4
             j = n_energy_out - 1
          elseif (r > 0.5) then
             ! last bin has relative probability of 0.1
             j = n_energy_out
          elseif (r > 0.1) then
             ! second bin has relative probability of 0.4
             j = 2
          else
             ! first bin has relative probability of 0.1
             j = 1
          end if
       else
          message = "Invalid secondary energy mode on S(a,b) table " // &
               trim(sab % name)
       end if

       ! Determine outgoing energy corresponding to E_in(i) and E_in(i+1)
       E_ij  = sab % inelastic_e_out(j,i)
       E_i1j = sab % inelastic_e_out(j,i+1)

       ! Outgoing energy
       E = (1 - f)*E_ij + f*E_i1j

       ! Sample outgoing cosine bin
       k = 1 + int(prn() * sab % n_inelastic_mu)

       ! Determine outgoing cosine corresponding to E_in(i) and E_in(i+1)
       mu_ijk  = sab % inelastic_mu(k,j,i)
       mu_i1jk = sab % inelastic_mu(k,j,i+1)

       ! Cosine of angle between incoming and outgoing neutron
       mu = (1 - f)*mu_ijk + f*mu_i1jk
    end if

    ! copy directional cosines
    u = p % coord0 % uvw(1)
    v = p % coord0 % uvw(2)
    w = p % coord0 % uvw(3)

    ! change direction of particle
    call rotate_angle(u, v, w, mu)
    p % coord0 % uvw = (/ u, v, w /)

    ! change energy of particle
    p % E = E

    ! Copy scattering cosine for tallies
    p % mu = mu

  end subroutine sab_scatter

!===============================================================================
! SAMPLE_TARGET_VELOCITY samples the target velocity based on the free gas
! scattering formulation used by most Monte Carlo codes. Excellent documentation
! for this method can be found in FRA-TM-123.
!===============================================================================

  subroutine sample_target_velocity(p, nuc, v_target)

    type(Particle), pointer :: p
    type(Nuclide),  pointer :: nuc
    real(8), intent(out)    :: v_target(3)

    real(8) :: u, v, w     ! direction of target 
    real(8) :: kT          ! equilibrium temperature of target in MeV
    real(8) :: alpha       ! probability of sampling f2 over f1
    real(8) :: mu          ! cosine of angle between neutron and target vel
    real(8) :: r1, r2      ! pseudo-random numbers
    real(8) :: c           ! cosine used in maxwell sampling
    real(8) :: accept_prob ! probability of accepting combination of vt and mu
    real(8) :: beta_vn     ! beta * speed of neutron
    real(8) :: beta_vt     ! beta * speed of target
    real(8) :: beta_vt_sq  ! (beta * speed of target)^2
    real(8) :: vt          ! speed of target

    ! Determine equilibrium temperature in MeV
    kT = nuc % kT

    ! Check if energy is above threshold
    if (p % E >= FREE_GAS_THRESHOLD * kT) then
       v_target = ZERO
       return
    end if

    ! calculate beta
    beta_vn = sqrt(nuc%awr * p%E / kT)

    alpha = ONE/(ONE + sqrt(pi)*beta_vn/2.0)

    do
       ! Sample two random numbers
       r1 = prn()
       r2 = prn()

       if (prn() < alpha) then
          ! With probability alpha, we sample the distribution p(y) =
          ! y*e^(-y). This can be done with sampling scheme C45 frmo the Monte
          ! Carlo sampler

          beta_vt_sq = -log(r1*r2)

       else
          ! With probability 1-alpha, we sample the distribution p(y) = y^2 *
          ! e^(-y^2). This can be done with sampling scheme C61 from the Monte
          ! Carlo sampler

          c = cos(PI/2.0 * prn())
          beta_vt_sq = -log(r1) - log(r2)*c*c
       end if

       ! Determine beta * vt
       beta_vt = sqrt(beta_vt_sq)

       ! Sample cosine of angle between neutron and target velocity
       mu = 2.0*prn() - ONE

       ! Determine rejection probability
       accept_prob = sqrt(beta_vn*beta_vn + beta_vt_sq - 2*beta_vn*beta_vt*mu) &
            /(beta_vn + beta_vt)

       ! Perform rejection sampling on vt and mu
       if (prn() < accept_prob) exit
    end do

    ! determine direction of target velocity based on the neutron's velocity
    ! vector and the sampled angle between them
    u = p % coord0 % uvw(1)
    v = p % coord0 % uvw(2)
    w = p % coord0 % uvw(3)
    call rotate_angle(u, v, w, mu)

    ! determine speed of target nucleus
    vt = sqrt(beta_vt_sq*kT/nuc % awr)

    ! determine velocity vector of target nucleus
    v_target(1) = u*vt
    v_target(2) = v*vt
    v_target(3) = W*vt

  end subroutine sample_target_velocity

!===============================================================================
! CREATE_FISSION_SITES determines the average total, prompt, and delayed
! neutrons produced from fission and creates appropriate bank sites.
!===============================================================================

  subroutine create_fission_sites(p, index_nuclide, rxn, event)

    type(Particle), pointer :: p
    integer, intent(in)     :: index_nuclide
    type(Reaction), pointer :: rxn
    logical, optional       :: event

    integer :: i            ! loop index
    integer :: j            ! index on nu energy grid / precursor group
    integer :: lc           ! index before start of energies/nu values
    integer :: NR           ! number of interpolation regions
    integer :: NE           ! number of energies tabulated
    integer :: nu           ! actual number of neutrons produced
    integer :: law          ! energy distribution law
    integer :: n_sample     ! number of times resampling
    real(8) :: E            ! incoming energy of neutron
    real(8) :: E_out        ! outgoing energy of fission neutron
    real(8) :: nu_t         ! total nu
    real(8) :: nu_p         ! prompt nu
    real(8) :: nu_d         ! delayed nu
    real(8) :: mu           ! fission neutron angular cosine
    real(8) :: phi          ! fission neutron azimuthal angle
    real(8) :: beta         ! delayed neutron fraction
    real(8) :: xi           ! random number
    real(8) :: yield        ! delayed neutron precursor yield
    real(8) :: prob         ! cumulative probability
    logical :: actual_event ! did fission actually occur? (no survival biasing)
    type(Nuclide),    pointer :: nuc
    type(DistEnergy), pointer :: edist => null()

    ! Get pointer to nuclide
    nuc => nuclides(index_nuclide)

    ! check whether actual fission event occurred for when survival biasing is
    ! turned off -- assume by default that no event occurs
    if (present(event)) then
       actual_event = event
    else
       actual_event = .false.
    end if

    ! copy energy of neutron
    E = p % E

    ! Determine total nu
    nu_t = nu_total(nuc, E)

    ! Determine prompt nu
    if (nuc % nu_p_type == NU_NONE) then
       nu_p = nu_t
    else
       nu_p = nu_prompt(nuc, E)
    end if
          
    ! Determine delayed nu
    nu_d = nu_delayed(nuc, E)

    ! Determine delayed neutron fraction
    beta = nu_d / nu_t

    ! TODO: Heat generation from fission

    ! Sample number of neutrons produced
    if (actual_event) then
       nu_t = p % wgt / keff * nu_t
    else
       nu_t = p % last_wgt * micro_xs(index_nuclide) % fission / (keff * &
            micro_xs(index_nuclide) % total) * nu_t
    end if
    if (prn() > nu_t - int(nu_t)) then
       nu = int(nu_t)
    else
       nu = int(nu_t) + 1
    end if

    ! Bank source neutrons
    if (nu == 0 .or. n_bank == 3*work) return
    do i = int(n_bank,4) + 1, int(min(n_bank + nu, 3*work),4)
       ! Bank source neutrons by copying particle data
       fission_bank(i) % id  = p % id
       fission_bank(i) % xyz = p % coord0 % xyz

       ! sample cosine of angle
       mu = sample_angle(rxn, E)

       ! sample between delayed and prompt neutrons
       if (prn() < beta) then
          ! ====================================================================
          ! DELAYED NEUTRON SAMPLED

          ! sampled delayed precursor group
          xi = prn()
          lc = 1
          prob = ZERO
          do j = 1, nuc % n_precursor
             ! determine number of interpolation regions and energies
             NR = int(nuc % nu_d_precursor_data(lc + 1))
             NE = int(nuc % nu_d_precursor_data(lc + 2 + 2*NR))

             ! determine delayed neutron precursor yield for group j
             yield = interpolate_tab1(nuc % nu_d_precursor_data( &
                  lc+1:lc+2+2*NR+2*NE), E)

             ! Check if this group is sampled
             prob = prob + yield
             if (xi < prob) exit

             ! advance pointer
             lc = lc + 2 + 2*NR + 2*NE + 1
          end do

          ! select energy distribution for group j
          law = nuc % nu_d_edist(j) % law
          edist => nuc % nu_d_edist(j)

          ! sample from energy distribution
          n_sample = 0
          do
             if (law == 44 .or. law == 61) then
                call sample_energy(edist, E, E_out, mu)
             else
                call sample_energy(edist, E, E_out)
             end if

             ! resample if energy is >= 20 MeV
             if (E_out < 20) exit

             ! check for large number of resamples
             n_sample = n_sample + 1
             if (n_sample == MAX_SAMPLE) then
                message = "Resampled energy distribution maximum number of " // &
                     "times for nuclide " // nuc % name
                call fatal_error()
             end if
          end do

       else
          ! ====================================================================
          ! PROMPT NEUTRON SAMPLED

          ! sample from prompt neutron energy distribution
          law = rxn % edist % law
          n_sample = 0
          do
             if (law == 44 .or. law == 61) then
                call sample_energy(rxn%edist, E, E_out, prob)
             else
                call sample_energy(rxn%edist, E, E_out)
             end if

             ! resample if energy is >= 20 MeV
             if (E_out < 20) exit

             ! check for large number of resamples
             n_sample = n_sample + 1
             if (n_sample == MAX_SAMPLE) then
                message = "Resampled energy distribution maximum number of " // &
                     "times for nuclide " // nuc % name
                call fatal_error()
             end if
          end do

       end if

       ! Sample azimuthal angle uniformly in [0,2*pi)
       phi = TWO*PI*prn()
       fission_bank(i) % uvw(1) = mu
       fission_bank(i) % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
       fission_bank(i) % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

       ! set energy of fission neutron
       fission_bank(i) % E = E_out
    end do

    ! increment number of bank sites
    n_bank = min(n_bank + nu, 3*work)
    p % n_bank = nu

  end subroutine create_fission_sites

!===============================================================================
! INELASTIC_SCATTER handles all reactions with a single secondary neutron (other
! than fission), i.e. level scattering, (n,np), (n,na), etc.
!===============================================================================

  subroutine inelastic_scatter(p, nuc, rxn)

    type(Particle), pointer :: p
    type(Nuclide),  pointer :: nuc
    type(Reaction), pointer :: rxn

    integer :: n_secondary ! number of secondary particles
    integer :: law         ! secondary energy distribution law
    real(8) :: A           ! atomic weight ratio of nuclide
    real(8) :: E_in        ! incoming energy
    real(8) :: mu          ! cosine of scattering angle
    real(8) :: E           ! outgoing energy in laboratory
    real(8) :: E_cm        ! outgoing energy in center-of-mass
    real(8) :: u,v,w       ! direction cosines
    real(8) :: Q           ! Q-value of reaction
    
    ! copy energy of neutron
    E_in = p % E

    ! determine A and Q
    A = nuc % awr
    Q = rxn % Q_value

    ! determine secondary energy distribution law
    law = rxn % edist % law

    ! sample scattering angle
    mu = sample_angle(rxn, E_in)

    ! sample outgoing energy
    if (law == 44 .or. law == 61) then
       call sample_energy(rxn%edist, E_in, E, mu)
    elseif (law == 66) then
       call sample_energy(rxn%edist, E_in, E, A=A, Q=Q)
    else
       call sample_energy(rxn%edist, E_in, E)
    end if

    ! if scattering system is in center-of-mass, transfer cosine of scattering
    ! angle and outgoing energy from CM to LAB
    if (rxn % TY < 0) then
       E_cm = E

       ! determine outgoing energy in lab
       E = E_cm + (E_in + TWO * mu * (A+ONE) * sqrt(E_in * E_cm)) & 
            & / ((A+ONE)*(A+ONE))

       ! determine outgoing angle in lab
       mu = mu * sqrt(E_cm/E) + ONE/(A+ONE) * sqrt(E_in/E)
    end if

    ! copy directional cosines
    u = p % coord0 % uvw(1)
    v = p % coord0 % uvw(2)
    w = p % coord0 % uvw(3)

    ! change direction of particle
    call rotate_angle(u, v, w, mu)
    p % coord0 % uvw = (/ u, v, w /)

    ! change energy of particle
    p % E = E

    ! Copy scattering cosine for tallies
    p % mu = mu

    ! change weight of particle based on multiplicity
    n_secondary = abs(rxn % TY)
    p % wgt = n_secondary * p % wgt

  end subroutine inelastic_scatter

!===============================================================================
! SAMPLE_ANGLE samples the cosine of the angle between incident and exiting
! particle directions either from 32 equiprobable bins or from a tabular
! distribution.
!===============================================================================

  function sample_angle(rxn, E) result(mu)

    type(Reaction), pointer    :: rxn ! reaction
    real(8),        intent(in) :: E   ! incoming energy

    real(8)        :: xi      ! random number on [0,1)
    integer        :: interp  ! type of interpolation
    integer        :: type    ! angular distribution type
    integer        :: i       ! incoming energy bin
    integer        :: n       ! number of incoming energy bins
    integer        :: lc      ! location in data array
    integer        :: NP      ! number of points in cos distribution
    integer        :: k       ! index on cosine grid
    real(8)        :: r       ! interpolation factor on incoming energy
    real(8)        :: frac    ! interpolation fraction on cosine
    real(8)        :: mu0     ! cosine in bin k
    real(8)        :: mu1     ! cosine in bin k+1
    real(8)        :: mu      ! final cosine sampled
    real(8)        :: c_k     ! cumulative frequency at k
    real(8)        :: c_k1    ! cumulative frequency at k+1
    real(8)        :: p0,p1   ! probability distribution

    ! check if reaction has angular distribution -- if not, sample outgoing
    ! angle isotropically
    if (.not. rxn % has_angle_dist) then
       mu = TWO * prn() - ONE
       return
    end if

    ! determine number of incoming energies
    n = rxn % adist % n_energy

    ! find energy bin and calculate interpolation factor -- if the energy is
    ! outside the range of the tabulated energies, choose the first or last bins
    if (E < rxn % adist % energy(1)) then
       i = 1
       r = ZERO
    elseif (E > rxn % adist % energy(n)) then
       i = n - 1
       r = ONE
    else
       i = binary_search(rxn % adist % energy, n, E)
       r = (E - rxn % adist % energy(i)) / & 
            & (rxn % adist % energy(i+1) - rxn % adist % energy(i))
    end if

    ! Sample between the ith and (i+1)th bin
    if (r > prn()) i = i + 1

    ! check whether this is a 32-equiprobable bin or a tabular distribution
    lc  = rxn % adist % location(i)
    type = rxn % adist % type(i)
    if (type == ANGLE_ISOTROPIC) then
       mu = TWO * prn() - ONE
    elseif (type == ANGLE_32_EQUI) then
       ! sample cosine bin
       xi = prn()
       k = 1 + int(32.0_8*xi)

       ! calculate cosine
       mu0 = rxn % adist % data(lc + k)
       mu1 = rxn % adist % data(lc + k+1)
       mu = mu0 + (32.0_8 * xi - k) * (mu1 - mu0)

    elseif (type == ANGLE_TABULAR) then
       interp = int(rxn % adist % data(lc + 1))
       NP     = int(rxn % adist % data(lc + 2))

       ! determine outgoing cosine bin
       xi = prn()
       lc = lc + 2
       c_k = rxn % adist % data(lc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = rxn % adist % data(lc + 2*NP + k+1)
          if (xi < c_k1) exit
          c_k = c_k1
       end do

       p0  = rxn % adist % data(lc + NP + k)
       mu0 = rxn % adist % data(lc + k)
       if (interp == HISTOGRAM) then
          ! Histogram interpolation
          mu = mu0 + (xi - c_k)/p0

       elseif (interp == LINEAR_LINEAR) then
          ! Linear-linear interpolation -- not sure how you come about the
          ! formula given in the MCNP manual
          p1  = rxn % adist % data(lc + NP + k+1)
          mu1 = rxn % adist % data(lc + k+1)

          frac = (p1 - p0)/(mu1 - mu0)
          if (frac == ZERO) then
             mu = mu0 + (xi - c_k)/p0
          else
             mu = mu0 + (sqrt(p0*p0 + 2*frac*(xi - c_k))-p0)/frac
          end if
       else
          message = "Unknown interpolation type: " // trim(int_to_str(interp))
          call fatal_error()
       end if

       if (abs(mu) > ONE) then
          message = "Sampled cosine of angle outside [-1, 1)."
          call warning()

          mu = sign(ONE,mu)
       end if
         
    else
       message = "Unknown angular distribution type: " // trim(int_to_str(type))
       call fatal_error()
    end if
    
  end function sample_angle

!===============================================================================
! ROTATE_ANGLE rotates direction cosines through a polar angle whose cosine is
! mu and through an azimuthal angle sampled uniformly. Note that this is done
! with direct sampling rather than rejection as is done in MCNP and SERPENT.
!===============================================================================

  subroutine rotate_angle(u, v, w, mu)

!    type(Particle), pointer :: p
    real(8), intent(inout) :: u
    real(8), intent(inout) :: v
    real(8), intent(inout) :: w
    real(8), intent(in)    :: mu ! cosine of angle in lab

    real(8) :: phi, sinphi, cosphi
    real(8) :: a,b
    real(8) :: u0, v0, w0

    ! Copy original directional cosines
    u0 = u
    v0 = v
    w0 = w

    ! Sample azimuthal angle in [0,2pi)
    phi = TWO * PI * prn()

    ! Precompute factors to save flops
    sinphi = sin(phi)
    cosphi = cos(phi)
    a = sqrt(ONE - mu*mu)
    b = sqrt(ONE - w0*w0)

    ! Need to treat special case where sqrt(1 - w**2) is close to zero by
    ! expanding about the v component rather than the w component
    if (b > 1e-10) then
       u = mu*u0 + a*(u0*w0*cosphi - v0*sinphi)/b
       v = mu*v0 + a*(v0*w0*cosphi + u0*sinphi)/b
       w = mu*w0 - a*b*cosphi
    else
       b = sqrt(ONE - v0*v0)
       u = mu*u0 + a*(u0*v0*cosphi + w0*sinphi)/b
       v = mu*v0 - a*b*cosphi
       w = mu*w0 + a*(v0*w0*cosphi - u0*sinphi)/b
    end if

  end subroutine rotate_angle
    
!===============================================================================
! SAMPLE_ENERGY
!===============================================================================

  recursive subroutine sample_energy(edist, E_in, E_out, mu_out, A, Q)

    type(DistEnergy),  pointer       :: edist
    real(8), intent(in)              :: E_in
    real(8), intent(out)             :: E_out
    real(8), intent(inout), optional :: mu_out
    real(8), intent(in),    optional :: A
    real(8), intent(in),    optional :: Q

    integer :: i           ! index on incoming energy grid
    integer :: k           ! sampled index on outgoing grid
    integer :: l           ! sampled index on incoming grid
    integer :: n_sample    ! number of rejections
    integer :: lc          ! dummy index
    integer :: NR          ! number of interpolation regions
    integer :: NE          ! number of energies
    integer :: NET         ! number of outgoing energies
    integer :: INTTp       ! combination of INTT and ND
    integer :: INTT        ! 1 = histogram, 2 = linear-linear
    integer :: JJ          ! 1 = histogram, 2 = linear-linear
    integer :: ND          ! number of discrete lines
    integer :: NP          ! number of points in distribution

    real(8) :: p_valid     ! probability of law validity

    real(8) :: E_i_1, E_i_K   ! endpoints on outgoing grid i
    real(8) :: E_i1_1, E_i1_K ! endpoints on outgoing grid i+1
    real(8) :: E_1, E_K       ! endpoints interpolated between i and i+1

    real(8) :: E_l_k, E_l_k1  ! adjacent E on outgoing grid l
    real(8) :: p_l_k, p_l_k1  ! adjacent p on outgoing grid l
    real(8) :: c_k, c_k1      ! cumulative probability

    real(8) :: KM_A           ! Kalbach-Mann parameter R
    real(8) :: KM_R           ! Kalbach-Mann parameter R
    real(8) :: A_k, A_k1      ! Kalbach-Mann A on outgoing grid l
    real(8) :: R_k, R_k1      ! Kalbach-Mann R on outgoing grid l

    real(8) :: Watt_a, Watt_b ! Watt spectrum parameters

    real(8) :: mu_k    ! angular cosine in bin k
    real(8) :: mu_k1   ! angular cosine in bin k+1
    real(8) :: p_k     ! angular pdf in bin k
    real(8) :: p_k1    ! angular pdf in bin k+1

    real(8) :: E_cm
    real(8) :: r           ! interpolation factor on incoming energy
    real(8) :: frac        ! interpolation factor on outgoing energy
    real(8) :: U           ! restriction energy
    real(8) :: T           ! nuclear temperature

    real(8) :: Ap          ! total mass ratio for n-body dist
    integer :: n_bodies    ! number of bodies for n-body dist
    real(8) :: E_max       ! parameter for n-body dist
    real(8) :: x, y, v     ! intermediate variables for n-body dist
    real(8) :: r1, r2, r3, r4, r5, r6

    ! ==========================================================================
    ! SAMPLE ENERGY DISTRIBUTION IF THERE ARE MULTIPLE

    if (associated(edist % next)) then
       if (edist % p_valid % n_regions > 0) then
          p_valid = interpolate_tab1(edist % p_valid, E_in)

          if (prn() > p_valid) then
             if (edist % law == 44 .or. edist % law == 61) then
                call sample_energy(edist%next, E_in, E_out, mu_out)
             elseif (edist % law == 66) then
                call sample_energy(edist%next, E_in, E_out, A=A, Q=Q)
             else
                call sample_energy(edist%next, E_in, E_out)
             end if
             return
          end if
       end if
    end if

    ! Determine which secondary energy distribution law to use
    select case (edist % law)
    case (1)
       ! =======================================================================
       ! TABULAR EQUIPROBABLE ENERGY BINS

       ! read number of interpolation regions, incoming energies, and outgoing
       ! energies
       NR  = int(edist % data(1))
       NE  = int(edist % data(2 + 2*NR))
       NET = int(edist % data(3 + 2*NR + NE))
       if (NR > 0) then
          message = "Multiple interpolation regions not supported while &
               &attempting to sample equiprobable energy bins."
          call fatal_error()
       end if

       ! determine index on incoming energy grid and interpolation factor
       lc = 2 + 2*NR
       i = binary_search(edist % data(lc+1), NE, E_in)
       r = (E_in - edist%data(lc+i)) / &
            & (edist%data(lc+i+1) - edist%data(lc+i))

       ! Sample outgoing energy bin
       r1 = prn()
       k = 1 + int(NET * r1)

       ! Determine E_1 and E_K
       lc    = 3 + 3*NR + NE + (i-1)*NET
       E_i_1 = edist % data(lc + 1)
       E_i_K = edist % data(lc + NET)

       lc     = 3 + 3*NR + NE + i*NET
       E_i1_1 = edist % data(lc + 1)
       E_i1_K = edist % data(lc + NET)

       E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
       E_K = E_i_K + r*(E_i1_K - E_i_K)

       ! Randomly select between the outgoing table for incoming energy E_i and
       ! E_(i+1)
       if (prn() < r) then
          l = i + 1
       else
          l = i
       end if

       ! Determine E_l_k and E_l_k+1
       lc     = 3 + 2*NR + NE + (l-1)*NET
       E_l_k  = edist % data(lc+k)
       E_l_k1 = edist % data(lc+k+1)

       ! Determine E' (denoted here as E_out)
       r2 = prn()
       E_out  = E_l_k + r2*(E_l_k1 - E_l_k)

       ! Now interpolate between incident energy bins i and i + 1
       if (l == i) then
          E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
       else
          E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
       end if

    case (3)
       ! =======================================================================
       ! INELASTIC LEVEL SCATTERING

       E_cm = edist%data(2) * (E_in - edist%data(1))
       
       E_out = E_cm

    case (4)
       ! =======================================================================
       ! CONTINUOUS TABULAR DISTRIBUTION

       ! read number of interpolation regions and incoming energies 
       NR  = int(edist % data(1))
       NE  = int(edist % data(2 + 2*NR))
       if (NR > 0) then
          message = "Multiple interpolation regions not supported while &
               &attempting to sample continuous tabular distribution."
          call fatal_error()
       end if

       ! find energy bin and calculate interpolation factor -- if the energy is
       ! outside the range of the tabulated energies, choose the first or last
       ! bins
       lc = 2 + 2*NR
       if (E_in < edist % data(lc+1)) then
          i = 1
          r = ZERO
       elseif (E_in > edist % data(lc+NE)) then
          i = NE - 1
          r = ONE
       else
          i = binary_search(edist % data(lc+1), NE, E_in)
          r = (E_in - edist%data(lc+i)) / & 
               & (edist%data(lc+i+1) - edist%data(lc+i))
       end if

       ! Sample between the ith and (i+1)th bin
       r2 = prn()
       if (r > r2) then
          l = i + 1
       else
          l = i
       end if

       ! interpolation for energy E1 and EK
       lc    = int(edist%data(2 + 2*NR + NE + i))
       NP    = int(edist%data(lc + 2))
       E_i_1 = edist%data(lc + 2 + 1)
       E_i_K = edist%data(lc + 2 + NP)

       lc     = int(edist%data(2 + 2*NR + NE + i + 1))
       NP     = int(edist%data(lc + 2))
       E_i1_1 = edist%data(lc + 2 + 1)
       E_i1_K = edist%data(lc + 2 + NP)

       E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
       E_K = E_i_K + r*(E_i1_K - E_i_K)

       ! determine location of outgoing energies, pdf, cdf for E(l)
       lc = int(edist % data(2 + 2*NR + NE + l))

       ! determine type of interpolation and number of discrete lines
       INTTp = int(edist % data(lc + 1))
       NP    = int(edist % data(lc + 2))
       if (INTTp > 10) then
          INTT = mod(INTTp,10)
          ND = (INTTp - INTT)/10
       else
          INTT = INTTp
          ND = 0
       end if

       if (ND > 0) then
          ! discrete lines present
          message = "Discrete lines in continuous tabular distributed not &
               &yet supported"
          call fatal_error()
       end if

       ! determine outgoing energy bin
       r1 = prn()
       lc = lc + 2 ! start of EOUT
       c_k = edist % data(lc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = edist % data(lc + 2*NP + k+1)
          if (r1 < c_k1) exit
          c_k = c_k1
       end do

       E_l_k = edist % data(lc+k)
       p_l_k = edist % data(lc+NP+k)
       if (INTT == HISTOGRAM) then
          ! Histogram interpolation
          E_out = E_l_k + (r1 - c_k)/p_l_k

       elseif (INTT == LINEAR_LINEAR) then
          ! Linear-linear interpolation -- not sure how you come about the
          ! formula given in the MCNP manual
          E_l_k1 = edist % data(lc+k+1)
          p_l_k1 = edist % data(lc+NP+k+1)

          frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
          if (frac == ZERO) then
             E_out = E_l_k + (r1 - c_k)/p_l_k
          else
             E_out = E_l_k + (sqrt(p_l_k*p_l_k + 2*frac*(r1 - c_k)) - & 
                  & p_l_k)/frac
          end if
       else
          message = "Unknown interpolation type: " // trim(int_to_str(INTT))
          call fatal_error()
       end if

       ! Now interpolate between incident energy bins i and i + 1
       if (l == i) then
          E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
       else
          E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
       end if

    case (5)
       ! =======================================================================
       ! GENERAL EVAPORATION SPECTRUM

    case (7)
       ! =======================================================================
       ! MAXWELL FISSION SPECTRUM

       ! read number of interpolation regions and incoming energies 
       NR = int(edist % data(1))
       NE = int(edist % data(2 + 2*NR))

       ! determine nuclear temperature from tabulated function
       T = interpolate_tab1(edist % data, E_in)

       ! determine restriction energy
       lc = 2 + 2*NR + 2*NE
       U = edist % data(lc + 1)
       
       n_sample = 0
       do
          ! sample maxwell fission spectrum
          E_out = maxwell_spectrum(T)

          ! accept energy based on restriction energy
          if (E_out <= E_in - U) exit

          ! check for large number of rejections
          n_sample = n_sample + 1
          if (n_sample == MAX_SAMPLE) then
             message = "Too many rejections on Maxwell fission spectrum."
             call fatal_error()
          end if
       end do
       
    case (9)
       ! =======================================================================
       ! EVAPORATION SPECTRUM

       ! read number of interpolation regions and incoming energies 
       NR = int(edist % data(1))
       NE = int(edist % data(2 + 2*NR))

       ! determine nuclear temperature from tabulated function
       T = interpolate_tab1(edist % data, E_in)

       ! determine restriction energy
       lc = 2 + 2*NR + 2*NE
       U = edist % data(lc + 1)

       ! sample outgoing energy based on evaporation spectrum probability
       ! density function
       n_sample = 0
       do
          r1 = prn()
          r2 = prn()
          E_out = -T * log(r1*r2)
          if (E_out <= E_in - U) exit

          ! check for large number of rejections
          n_sample = n_sample + 1
          if (n_sample == MAX_SAMPLE) then
             message = "Too many rejections on evaporation spectrum."
             call fatal_error()
          end if
       end do
       
    case (11)
       ! =======================================================================
       ! ENERGY-DEPENDENT WATT SPECTRUM

       ! read number of interpolation regions and incoming energies for
       ! parameter 'a'
       NR = int(edist % data(1))
       NE = int(edist % data(2 + 2*NR))

       ! determine Watt parameter 'a' from tabulated function
       Watt_a = interpolate_tab1(edist % data, E_in)

       ! determine Watt parameter 'b' from tabulated function
       lc = 2 + 2*(NR + NE)
       Watt_b = interpolate_tab1(edist % data, E_in, lc + 1)

       ! read number of interpolation regions and incoming energies for
       ! parameter 'a'
       NR = int(edist % data(lc + 1))
       NE = int(edist % data(lc + 2 + 2*NR))

       ! determine restriction energy
       lc = lc + 2 + 2*(NR + NE)
       U = edist % data(lc + 1)

       n_sample = 0
       do
          ! Sample energy-dependent Watt fission spectrum
          E_out = watt_spectrum(Watt_a, Watt_b)

          ! accept energy based on restriction energy
          if (E_out <= E_in - U) exit

          ! check for large number of rejections
          n_sample = n_sample + 1
          if (n_sample == MAX_SAMPLE) then
             message = "Too many rejections on Watt spectrum."
             call fatal_error()
          end if
       end do

    case (44)
       ! =======================================================================
       ! KALBACH-MANN CORRELATED SCATTERING

       if (.not. present(mu_out)) then
          message = "Law 44 called without giving mu_out as argument."
          call fatal_error()
       end if

       ! read number of interpolation regions and incoming energies 
       NR = int(edist % data(1))
       NE = int(edist % data(2 + 2*NR))
       if (NR > 0) then
          message = "Multiple interpolation regions not supported while &
               &attempting to sample Kalbach-Mann distribution."
          call fatal_error()
       end if

       ! find energy bin and calculate interpolation factor -- if the energy is
       ! outside the range of the tabulated energies, choose the first or last
       ! bins
       lc = 2 + 2*NR
       if (E_in < edist % data(lc+1)) then
          i = 1
          r = ZERO
       elseif (E_in > edist % data(lc+NE)) then
          i = NE - 1
          r = ONE
       else
          i = binary_search(edist % data(lc+1), NE, E_in)
          r = (E_in - edist%data(lc+i)) / & 
               & (edist%data(lc+i+1) - edist%data(lc+i))
       end if

       ! Sample between the ith and (i+1)th bin
       r2 = prn()
       if (r > r2) then
          l = i + 1
       else
          l = i
       end if

       ! determine endpoints on grid i
       lc    = int(edist%data(2+2*NR+NE + i)) ! start of LDAT for i
       NP    = int(edist%data(lc + 2))
       E_i_1 = edist%data(lc + 2 + 1)
       E_i_K = edist%data(lc + 2 + NP)

       ! determine endpoints on grid i+1
       lc     = int(edist%data(2+2*NR+NE + i+1)) ! start of LDAT for i+1
       NP     = int(edist%data(lc + 2))
       E_i1_1 = edist%data(lc + 2 + 1)
       E_i1_K = edist%data(lc + 2 + NP)

       E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
       E_K = E_i_K + r*(E_i1_K - E_i_K)

       ! determine location of outgoing energies, pdf, cdf for E(l)
       lc = int(edist % data(2 + 2*NR + NE + l))

       ! determine type of interpolation and number of discrete lines
       INTTp = int(edist % data(lc + 1))
       NP    = int(edist % data(lc + 2))
       if (INTTp > 10) then
          INTT = mod(INTTp,10)
          ND = (INTTp - INTT)/10
       else
          INTT = INTTp
          ND = 0
       end if

       if (ND > 0) then
          ! discrete lines present
          message = "Discrete lines in continuous tabular distributed not &
               &yet supported"
          call fatal_error()
       end if

       ! determine outgoing energy bin
       r1 = prn()
       lc = lc + 2 ! start of EOUT
       c_k = edist % data(lc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = edist % data(lc + 2*NP + k+1)
          if (r1 < c_k1) exit
          c_k = c_k1
       end do

       E_l_k = edist % data(lc+k)
       p_l_k = edist % data(lc+NP+k)
       if (INTT == HISTOGRAM) then
          ! Histogram interpolation
          E_out = E_l_k + (r1 - c_k)/p_l_k

          ! Determine Kalbach-Mann parameters
          KM_R = edist % data(lc + 3*NP + k)
          KM_A = edist % data(lc + 4*NP + k)

       elseif (INTT == LINEAR_LINEAR) then
          ! Linear-linear interpolation -- not sure how you come about the
          ! formula given in the MCNP manual
          E_l_k1 = edist % data(lc+k+1)
          p_l_k1 = edist % data(lc+NP+k+1)

          ! Find E prime
          frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
          if (frac == ZERO) then
             E_out = E_l_k + (r1 - c_k)/p_l_k
          else
             E_out = E_l_k + (sqrt(p_l_k*p_l_k + 2*frac*(r1 - c_k)) - & 
                  & p_l_k)/frac
          end if

          ! Determine Kalbach-Mann parameters
          R_k  = edist % data(lc + 3*NP + k)
          R_k1 = edist % data(lc + 3*NP + k+1)
          A_k  = edist % data(lc + 4*NP + k)
          A_k1 = edist % data(lc + 4*NP + k+1)
          
          KM_R = R_k + (R_k1 - R_k)*(E_out - E_l_k)/(E_l_k1 - E_l_k)
          KM_A = A_k + (A_k1 - A_k)*(E_out - E_l_k)/(E_l_k1 - E_l_k)
       else
          message = "Unknown interpolation type: " // trim(int_to_str(INTT))
          call fatal_error()
       end if

       ! Now interpolate between incident energy bins i and i + 1
       if (l == i) then
          E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
       else
          E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
       end if

       ! Sampled correlated angle from Kalbach-Mann parameters
       r3 = prn()
       r4 = prn()
       T = (TWO*r4 - ONE) * sinh(KM_A)
       if (r3 > KM_R) then
          mu_out = log(T + sqrt(T*T + ONE))/KM_A
       else
          mu_out = log(r4*exp(KM_A) + (ONE - r4)*exp(-KM_A))/KM_A
       end if

    case (61)
       ! =======================================================================
       ! CORRELATED ENERGY AND ANGLE DISTRIBUTION

       if (.not. present(mu_out)) then
          message = "Law 44 called without giving mu_out as argument."
          call fatal_error()
       end if

       ! read number of interpolation regions and incoming energies 
       NR = int(edist % data(1))
       NE = int(edist % data(2 + 2*NR))
       if (NR > 0) then
          message = "Multiple interpolation regions not supported while &
               &attempting to sample correlated energy-angle distribution."
          call fatal_error()
       end if

       ! find energy bin and calculate interpolation factor -- if the energy is
       ! outside the range of the tabulated energies, choose the first or last
       ! bins
       lc = 2 + 2*NR
       if (E_in < edist % data(lc+1)) then
          i = 1
          r = ZERO
       elseif (E_in > edist % data(lc+NE)) then
          i = NE - 1
          r = ONE
       else
          i = binary_search(edist % data(lc+1), NE, E_in)
          r = (E_in - edist%data(lc+i)) / & 
               & (edist%data(lc+i+1) - edist%data(lc+i))
       end if

       ! Sample between the ith and (i+1)th bin
       r2 = prn()
       if (r > r2) then
          l = i + 1
       else
          l = i
       end if

       ! determine endpoints on grid i
       lc    = int(edist%data(2+2*NR+NE + i)) ! start of LDAT for i
       NP    = int(edist%data(lc + 2))
       E_i_1 = edist%data(lc + 2 + 1)
       E_i_K = edist%data(lc + 2 + NP)

       ! determine endpoints on grid i+1
       lc     = int(edist%data(2+2*NR+NE + i+1)) ! start of LDAT for i+1
       NP     = int(edist%data(lc + 2))
       E_i1_1 = edist%data(lc + 2 + 1)
       E_i1_K = edist%data(lc + 2 + NP)

       E_1 = E_i_1 + r*(E_i1_1 - E_i_1)
       E_K = E_i_K + r*(E_i1_K - E_i_K)

       ! determine location of outgoing energies, pdf, cdf for E(l)
       lc = int(edist % data(2 + 2*NR + NE + l))

       ! determine type of interpolation and number of discrete lines
       INTTp = int(edist % data(lc + 1))
       NP    = int(edist % data(lc + 2))
       if (INTTp > 10) then
          INTT = mod(INTTp,10)
          ND = (INTTp - INTT)/10
       else
          INTT = INTTp
          ND = 0
       end if

       if (ND > 0) then
          ! discrete lines present
          message = "Discrete lines in continuous tabular distributed not &
               &yet supported"
          call fatal_error()
       end if

       ! determine outgoing energy bin
       r1 = prn()
       lc = lc + 2 ! start of EOUT
       c_k = edist % data(lc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = edist % data(lc + 2*NP + k+1)
          if (r1 < c_k1) exit
          c_k = c_k1
       end do

       E_l_k = edist % data(lc+k)
       p_l_k = edist % data(lc+NP+k)
       if (INTT == HISTOGRAM) then
          ! Histogram interpolation
          E_out = E_l_k + (r1 - c_k)/p_l_k

       elseif (INTT == LINEAR_LINEAR) then
          ! Linear-linear interpolation -- not sure how you come about the
          ! formula given in the MCNP manual
          E_l_k1 = edist % data(lc+k+1)
          p_l_k1 = edist % data(lc+NP+k+1)

          ! Find E prime
          frac = (p_l_k1 - p_l_k)/(E_l_k1 - E_l_k)
          if (frac == ZERO) then
             E_out = E_l_k + (r1 - c_k)/p_l_k
          else
             E_out = E_l_k + (sqrt(p_l_k*p_l_k + 2*frac*(r1 - c_k)) - & 
                  & p_l_k)/frac
          end if
       else
          message = "Unknown interpolation type: " // trim(int_to_str(INTT))
          call fatal_error()
       end if

       ! Now interpolate between incident energy bins i and i + 1
       if (l == i) then
          E_out = E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
       else
          E_out = E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
       end if

       ! Find location of correlated angular distribution
       lc = int(edist % data(lc+3*NP+k))

       ! Check if angular distribution is isotropic
       if (lc == 0) then
          mu_out = TWO * prn() - ONE
          return
       end if

       ! interpolation type and number of points in angular distribution
       JJ = int(edist % data(lc + 1))
       NP = int(edist % data(lc + 2))

       ! determine outgoing cosine bin
       r3 = prn()
       lc = lc + 2
       c_k = edist % data(lc + 2*NP + 1)
       do k = 1, NP-1
          c_k1 = edist % data(lc + 2*NP + k+1)
          if (r3 < c_k1) exit
          c_k = c_k1
       end do

       p_k  = edist % data(lc + NP + k)
       mu_k = edist % data(lc + k)
       if (JJ == HISTOGRAM) then
          ! Histogram interpolation
          mu_out = mu_k + (r3 - c_k)/p_k

       elseif (JJ == LINEAR_LINEAR) then
          ! Linear-linear interpolation -- not sure how you come about the
          ! formula given in the MCNP manual
          p_k1  = edist % data(lc + NP + k+1)
          mu_k1 = edist % data(lc + k+1)

          frac = (p_k1 - p_k)/(mu_k1 - mu_k)
          if (frac == ZERO) then
             mu_out = mu_k + (r3 - c_k)/p_k
          else
             mu_out = mu_k + (sqrt(p_k*p_k + 2*frac*(r3 - c_k))-p_k)/frac
          end if
       else
          message = "Unknown interpolation type: " // trim(int_to_str(JJ))
          call fatal_error()
       end if

    case (66)
       ! =======================================================================
       ! N-BODY PHASE SPACE DISTRIBUTION

       ! read number of bodies in phase space and total mass ratio
       n_bodies = int(edist % data(1))
       Ap       = edist % data(2)

       ! determine E_max parameter
       E_max = (Ap - ONE)/Ap * (A/(A+ONE)*E_in + Q)

       ! x is essentially a Maxwellian distribution
       x = maxwell_spectrum(ONE)

       select case (n_bodies)
       case (3)
          y = maxwell_spectrum(ONE)
       case (4)
          r1 = prn()
          r2 = prn()
          r3 = prn()
          y = -log(r1*r2*r3)
       case (5)
          r1 = prn()
          r2 = prn()
          r3 = prn()
          r4 = prn()
          r5 = prn()
          r6 = prn()
          y = -log(r1*r2*r3*r4) - log(r5) * cos(PI/2.*r6)**2
       end select

       ! now determine v and E_out
       v = x/(x+y)
       E_out = E_max * v

    case (67)
       ! =======================================================================
       ! LABORATORY ENERGY-ANGLE LAW

    end select
    
  end subroutine sample_energy

!===============================================================================
! MAXWELL_SPECTRUM samples an energy from the Maxwell fission distribution based
! on a direct sampling scheme. The probability distribution function for a
! Maxwellian is given as p(x) = 2/(T*sqrt(pi))*sqrt(x/T)*exp(-x/T). This PDF can
! be sampled using rule C64 in the Monte Carlo Sampler LA-9721-MS.
!===============================================================================

  function maxwell_spectrum(T) result(E_out)

    real(8), intent(in)  :: T     ! tabulated function of incoming E
    real(8)              :: E_out ! sampled energy

    real(8) :: r1, r2, r3  ! random numbers
    real(8) :: c           ! cosine of pi/2*r3

    r1 = prn()
    r2 = prn()
    r3 = prn()

    ! determine cosine of pi/2*r
    c = cos(PI/2.*r3)

    ! determine outgoing energy
    E_out = -T*(log(r1) + log(r2)*c*c)

  end function maxwell_spectrum

!===============================================================================
! WATT_SPECTRUM samples the outgoing energy from a Watt energy-dependent fission
! spectrum. Although fitted parameters exist for many nuclides, generally the
! continuous tabular distributions (LAW 4) should be used in lieu of the Watt
! spectrum. This direct sampling scheme is an unpublished scheme based on the
! original Watt spectrum derivation (See F. Brown's MC lectures).
!===============================================================================

  function watt_spectrum(a, b) result(E_out)

    real(8), intent(in) :: a     ! Watt parameter a
    real(8), intent(in) :: b     ! Watt parameter b
    real(8)             :: E_out ! energy of emitted neutron

    real(8) :: w ! sampled from Maxwellian

    w     = maxwell_spectrum(a)
    E_out = w + a*a*b/4. + (2.*prn() - ONE)*sqrt(a*a*b*w)

  end function watt_spectrum

!===============================================================================
! WIGNER samples a Wigner distribution of energy level spacings. Note that this
! scheme is C50 in the Monte Carlo Sampler from Los Alamos (LA-9721-MS).
!===============================================================================

  function wigner(D_avg) result (D)

    real(8), intent(in) :: D_avg ! average level spacing
    real(8)             :: D     ! sampled level spacing

    real(8) :: c

    c = -4.*D_avg*D_avg/PI * log(prn())
    D = sqrt(c)

  end function wigner

!===============================================================================
! CHI_SQUARED samples a chi-squared distribution with n degrees of freedom. The
! distribution of resonance widths in the unresolved region is given by a
! chi-squared distribution. For the special case of n=1, this is a Porter-Thomas
! distribution. For cases with n odd, rule C64 is used whereas for cases with n
! even, rule C45 is used.
!===============================================================================

  function chi_squared(n, G_avg) result(G)

    integer, intent(in)           :: n     ! number of degrees of freedom
    real(8), intent(in), optional :: G_avg ! average resonance width

    integer :: i       ! loop index
    real(8) :: G       ! sampled random variable (or resonance width)
    real(8) :: x, y, c ! dummy variables
    real(8) :: r1, r2  ! psuedorandom numbers

    select case (mod(n,2))
    case (0)
       ! Even number of degrees of freedom can be sampled via rule C45. We can
       ! sample x as -2/n*log(product(r_i, i = 1 to n/2))
       x = ONE
       do i = 1, n/2
          x = x * prn()
       end do
       x = -2./n * log(x)

    case (1)
       ! Odd number of degrees of freedom can be sampled via rule C64. We can
       ! sample x as -2/n*(log(r)*cos^2(pi/2*r) + log(product(r_i, i = 1 to
       ! floor(n/2)))

       ! Note that we take advantage of integer division on n/2
       y = ONE
       do i = 1, n/2
          y = y * prn()
       end do

       r1 = prn()
       r2 = prn()
       c = cos(PI/2.*r2)
       x = -2./n * (log(y) + log(r1)*c*c)
    end select

    ! If sampling a chi-squared distribution for a resonance width and the
    ! average resonance width has been given, return the sampled resonance
    ! width.
    if (present(G_avg)) then
       G = x * G_avg
    else
       G = x
    end if

  end function chi_squared

end module physics
