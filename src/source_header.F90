module source_header

  use, intrinsic :: ISO_C_BINDING

  use bank_header, only: Bank
  use constants
  use distribution_univariate
  use distribution_multivariate
  use error
  use geometry, only: find_cell
  use material_header, only: materials
  use nuclide_header, only: energy_min_neutron, energy_max_neutron
  use particle_header, only: Particle
  use string, only: to_lower
  use xml_interface

  implicit none
  private
  public :: free_memory_source
  public :: openmc_extend_sources
  public :: openmc_source_set_strength

  integer :: n_accept = 0  ! Number of samples accepted
  integer :: n_reject = 0  ! Number of samples rejected

!===============================================================================
! SOURCEDISTRIBUTION describes an external source of particles for a
! fixed-source problem or for the starting source in a k eigenvalue problem
!===============================================================================

  type, public :: SourceDistribution
    real(8) :: strength = ONE ! source strength
    class(SpatialDistribution),    allocatable :: space  ! spatial distribution
    class(UnitSphereDistribution), allocatable :: angle  ! angle distribution
    class(Distribution),           allocatable :: energy ! energy distribution
  contains
    procedure :: from_xml
    procedure :: sample
  end type SourceDistribution

  ! Number of external source distributions
  integer(C_INT32_T), public, bind(C) :: n_sources = 0

  ! External source distributions
  type(SourceDistribution), public, allocatable :: external_source(:)

contains

  subroutine from_xml(this, node, path_source)
    class(SourceDistribution), intent(inout) :: this
    type(XMLNode), intent(in) :: node
    character(MAX_FILE_LEN), intent(out) :: path_source

    integer :: n
    logical :: file_exists
    character(MAX_WORD_LEN) :: type
    type(XMLNode) :: node_space
    type(XMLNode) :: node_angle
    type(XMLNode) :: node_dist

    ! Check for source strength
    if (check_for_node(node, "strength")) then
      call get_node_value(node, "strength", this % strength)
    end if

    ! Check for external source file
    if (check_for_node(node, "file")) then
      ! Copy path of source file
      call get_node_value(node, "file", path_source)

      ! Check if source file exists
      inquire(FILE=path_source, EXIST=file_exists)
      if (.not. file_exists) then
        call fatal_error("Source file '" // trim(path_source) &
             // "' does not exist!")
      end if

    else

      ! Spatial distribution for external source
      if (check_for_node(node, "space")) then

        ! Get pointer to spatial distribution
        node_space = node % child("space")

        ! Check for type of spatial distribution
        type = ''
        if (check_for_node(node_space, "type")) &
             call get_node_value(node_space, "type", type)
        select case (to_lower(type))
        case ('cartesian')
          allocate(CartesianIndependent :: this % space)

        case ('box')
          allocate(SpatialBox :: this % space)

        case ('fission')
          allocate(SpatialBox :: this % space)
          select type(space => this % space)
          type is (SpatialBox)
            space % only_fissionable = .true.
          end select

        case ('point')
          allocate(SpatialPoint :: this % space)

        case default
          call fatal_error("Invalid spatial distribution for external source: "&
               // trim(type))
        end select

        ! Read spatial distribution from XML
        call this % space % from_xml(node_space)

      else
        ! If no spatial distribution specified, make it a point source
        allocate(SpatialPoint :: this % space)
        select type (space => this % space)
        type is (SpatialPoint)
          space % xyz(:) = [ZERO, ZERO, ZERO]
        end select
      end if

      ! Determine external source angular distribution
      if (check_for_node(node, "angle")) then

        ! Get pointer to angular distribution
        node_angle = node % child("angle")

        ! Check for type of angular distribution
        type = ''
        if (check_for_node(node_angle, "type")) &
             call get_node_value(node_angle, "type", type)
        select case (to_lower(type))
        case ('isotropic')
          allocate(Isotropic :: this % angle)

        case ('monodirectional')
          allocate(Monodirectional :: this % angle)

        case ('mu-phi')
          allocate(PolarAzimuthal :: this % angle)

        case default
          call fatal_error("Invalid angular distribution for external source: "&
               // trim(type))
        end select

        ! Read reference directional unit vector
        if (check_for_node(node_angle, "reference_uvw")) then
          n = node_word_count(node_angle, "reference_uvw")
          if (n /= 3) then
            call fatal_error('Angular distribution reference direction must have &
                 &three parameters specified.')
          end if
          call get_node_array(node_angle, "reference_uvw", &
               this % angle % reference_uvw)
        else
          ! By default, set reference unit vector to be positive z-direction
          this % angle % reference_uvw(:) = [ZERO, ZERO, ONE]
        end if

        ! Read parameters for angle distribution
        select type (angle => this % angle)
        type is (Monodirectional)
          call get_node_array(node_angle, "reference_uvw", &
               this % angle % reference_uvw)

        type is (PolarAzimuthal)
          if (check_for_node(node_angle, "mu")) then
            node_dist = node_angle % child("mu")
            call distribution_from_xml(angle % mu, node_dist)
          else
            allocate(Uniform :: angle%mu)
            select type (mu => angle%mu)
            type is (Uniform)
              mu % a = -ONE
              mu % b = ONE
            end select
          end if

          if (check_for_node(node_angle, "phi")) then
            node_dist = node_angle % child("phi")
            call distribution_from_xml(angle % phi, node_dist)
          else
            allocate(Uniform :: angle%phi)
            select type (phi => angle%phi)
            type is (Uniform)
              phi % a = ZERO
              phi % b = TWO*PI
            end select
          end if
        end select

      else
        ! Set default angular distribution isotropic
        allocate(Isotropic :: this % angle)
        this % angle % reference_uvw(:) = [ZERO, ZERO, ONE]
      end if

      ! Determine external source energy distribution
      if (check_for_node(node, "energy")) then
        node_dist = node % child("energy")
        call distribution_from_xml(this % energy, node_dist)
      else
        ! Default to a Watt spectrum with parameters 0.988 MeV and 2.249 MeV^-1
        allocate(Watt :: this % energy)
        select type(energy => this % energy)
        type is (Watt)
          energy % a = 0.988e6_8
          energy % b = 2.249e-6_8
        end select
      end if
    end if

  end subroutine from_xml

  function sample(this) result(site)
    class(SourceDistribution), intent(in) :: this
    type(Bank) :: site

    logical :: found      ! Does the source particle exist within geometry?
    type(Particle) :: p   ! Temporary particle for using find_cell

    ! Set weight to one by default
    site % wgt = ONE

    ! Repeat sampling source location until a good site has been found
    found = .false.
    do while (.not. found)
      ! Set particle defaults
      call p % initialize()

      ! Sample spatial distribution
      site % xyz(:) = this % space % sample()

      ! Fill p with needed data
      p % coord(1) % xyz(:) = site % xyz
      p % coord(1) % uvw(:) = [ ONE, ZERO, ZERO ]

      ! Now search to see if location exists in geometry
      call find_cell(p, found)

      ! Check if spatial site is in fissionable material
      if (found) then
        select type (space => this % space)
        type is (SpatialBox)
          if (space % only_fissionable) then
            if (p % material == MATERIAL_VOID) then
              found = .false.
            elseif (.not. materials(p % material) % fissionable) then
              found = .false.
            end if
          end if
        end select
      end if

      ! Check for rejection
      if (.not. found) then
        n_reject = n_reject + 1
        if (n_reject >= EXTSRC_REJECT_THRESHOLD .and. &
             real(n_accept, 8)/n_reject <= EXTSRC_REJECT_FRACTION) then
          call fatal_error("More than 95% of external source sites sampled &
               &were rejected. Please check your external source definition.")
        end if
      end if
    end do

    ! Increment number of accepted samples
    n_accept = n_accept + 1

    call p % clear()

    ! Sample angle
    site % uvw(:) = this % angle % sample()

    ! Check for monoenergetic source above maximum neutron energy
    select type (energy => this % energy)
    type is (Discrete)
      if (any(energy % x > energy_max_neutron)) then
        call fatal_error("Source energy above range of energies of at least &
             &one cross section table")
      else if (any(energy % x < energy_min_neutron)) then
        call fatal_error("Source energy below range of energies of at least &
             &one cross section table")
      end if
    end select

    do
      ! Sample energy spectrum
      site % E = this % energy % sample()

      ! Resample if energy falls outside minimum or maximum neutron energy
      if (site % E < energy_max_neutron .and. site % E > energy_min_neutron) exit
    end do

    ! Set delayed group
    site % delayed_group = 0

  end function sample

!===============================================================================
! FREE_MEMORY_SOURCE deallocates global arrays defined in this module
!===============================================================================

  subroutine free_memory_source()
    n_sources = 0
    if (allocated(external_source)) deallocate(external_source)
  end subroutine free_memory_source

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_extend_sources(n, index_start, index_end) result(err) bind(C)
    ! Extend the external_source array by n elements
    integer(C_INT32_T), value, intent(in) :: n
    integer(C_INT32_T), optional, intent(out) :: index_start
    integer(C_INT32_T), optional, intent(out) :: index_end
    integer(C_INT) :: err

    type(SourceDistribution), allocatable :: temp(:) ! temporary array

    if (n_sources == 0) then
      ! Allocate external_source array
      allocate(external_source(n))
    else
      ! Allocate external_source array with increased size
      allocate(temp(n_sources + n))

      ! Copy original source array to temporary array
      temp(1:n_sources) = external_source

      ! Move allocation from temporary array
      call move_alloc(FROM=temp, TO=external_source)
    end if

    ! Return indices in external_source array
    if (present(index_start)) index_start = n_sources + 1
    if (present(index_end)) index_end = n_sources + n
    n_sources = n_sources + n

    err = 0
  end function openmc_extend_sources


  function openmc_source_set_strength(index, strength) result(err) bind(C)
    integer(C_INT32_T), value, intent(in) :: index
    real(C_DOUBLE),     value, intent(in) :: strength
    integer(C_INT) :: err

    if (index >= 1 .and. index <= n_sources) then
      if (strength > ZERO) then
        external_source(index) % strength = strength
        err = 0
      else
        err = E_INVALID_ARGUMENT
        call set_errmsg("Source strength must be positive.")
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in external source array is out of bounds.")
    end if
  end function openmc_source_set_strength

end module source_header
