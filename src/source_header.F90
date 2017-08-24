module source_header

  use constants
  use distribution_univariate
  use distribution_multivariate
  use error, only: fatal_error
  use string, only: to_lower
  use xml_interface

  implicit none
  private

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
    procedure :: from_xml => source_from_xml
  end type SourceDistribution

  ! External source distributions
  type(SourceDistribution), public, allocatable :: external_source(:)

contains

  subroutine source_from_xml(this, node, path_source)
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

  end subroutine source_from_xml

end module source_header
