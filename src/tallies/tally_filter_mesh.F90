module tally_filter_mesh

  use, intrinsic :: ISO_C_BINDING

  use hdf5

  use constants
  use dict_header,         only: EMPTY
  use error
  use mesh_header,         only: RegularMesh, meshes, n_meshes, mesh_dict
  use hdf5_interface
  use particle_header,     only: Particle
  use string,              only: to_str
  use tally_filter_header
  use xml_interface

  implicit none
  private
  public :: openmc_mesh_filter_set_mesh

!===============================================================================
! MESHFILTER indexes the location of particle events to a regular mesh.  For
! tracklength tallies, it will produce multiple valid bins and the bin weight
! will correspond to the fraction of the track length that lies in that bin.
!===============================================================================

  type, public, extends(TallyFilter) :: MeshFilter
    integer :: mesh
  contains
    procedure :: from_xml
    procedure :: get_all_bins => get_all_bins_mesh
    procedure :: to_statepoint => to_statepoint_mesh
    procedure :: text_label => text_label_mesh
  end type MeshFilter

contains

  subroutine from_xml(this, node)
    class(MeshFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: i_mesh
    integer :: id
    integer :: n
    integer :: val

    n = node_word_count(node, "bins")

    if (n /= 1) call fatal_error("Only one mesh can be &
         &specified per mesh filter.")

    ! Determine id of mesh
    call get_node_value(node, "bins", id)

    ! Get pointer to mesh
    val = mesh_dict % get(id)
    if (val /= EMPTY) then
      i_mesh = val
    else
      call fatal_error("Could not find mesh " // trim(to_str(id)) &
           // " specified on filter.")
    end if

    ! Determine number of bins
    this % n_bins = product(meshes(i_mesh) % dimension)

    ! Store the index of the mesh
    this % mesh = i_mesh
  end subroutine from_xml

  subroutine get_all_bins_mesh(this, p, estimator, match)
    class(MeshFilter), intent(in)  :: this
    type(Particle),    intent(in)  :: p
    integer,           intent(in)  :: estimator
    type(TallyFilterMatch), intent(inout) :: match

    integer, parameter :: MAX_SEARCH_ITER = 100 ! Maximum number of times we can
                                                !  can loop while trying to find
                                                !  the first intersection.

    integer :: j                    ! loop index for direction
    integer :: n
    integer :: ijk0(3)              ! indices of starting coordinates
    integer :: ijk1(3)              ! indices of ending coordinates
    integer :: search_iter          ! loop count for intersection search
    integer :: bin
    real(8) :: weight               ! weight to be pushed back
    real(8) :: uvw(3)               ! cosine of angle of particle
    real(8) :: xyz0(3)              ! starting/intermediate coordinates
    real(8) :: xyz1(3)              ! ending coordinates of particle
    real(8) :: xyz_cross            ! coordinates of next boundary
    real(8) :: d(3)                 ! distance to each bounding surface
    real(8) :: total_distance       ! distance of entire particle track
    real(8) :: distance             ! distance traveled in mesh cell
    logical :: start_in_mesh        ! starting coordinates inside mesh?
    logical :: end_in_mesh          ! ending coordinates inside mesh?
    type(RegularMesh), pointer :: m

    weight = ERROR_REAL

    ! Get a pointer to the mesh.
    m => meshes(this % mesh)
    n = m % n_dimension

    if (estimator /= ESTIMATOR_TRACKLENGTH) then
      ! If this is an analog or collision tally, then there can only be one
      ! valid mesh bin.
      call m % get_bin(p % coord(1) % xyz, bin)
      if (bin /= NO_BIN_FOUND) then
        call match % bins % push_back(bin)
        call match % weights % push_back(ONE)
      end if
      return
    end if

    ! A track can span multiple mesh bins so we need to handle a lot of
    ! intersection logic for tracklength tallies.

    ! ========================================================================
    ! Determine if the track intersects the tally mesh.

    ! Copy the starting and ending coordinates of the particle.  Offset these
    ! just a bit for the purposes of determining if there was an intersection
    ! in case the mesh surfaces coincide with lattice/geometric surfaces which
    ! might produce finite-precision errors.
    xyz0 = p % last_xyz + TINY_BIT * p % coord(1) % uvw
    xyz1 = p % coord(1) % xyz - TINY_BIT * p % coord(1) % uvw

    ! Determine indices for starting and ending location.
    call m % get_indices(xyz0, ijk0(:n), start_in_mesh)
    call m % get_indices(xyz1, ijk1(:n), end_in_mesh)

    ! If this is the first iteration of the filter loop, check if the track
    ! intersects any part of the mesh.
    if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) then
      if (.not. m % intersects(xyz0, xyz1)) return
    end if

    ! ========================================================================
    ! Figure out which mesh cell to tally.

    ! Copy the un-modified coordinates the particle direction.
    xyz0 = p % last_xyz
    xyz1 = p % coord(1) % xyz
    uvw = p % coord(1) % uvw

    ! Compute the length of the entire track.
    total_distance = sqrt(sum((xyz1 - xyz0)**2))

    ! We are looking for the first valid mesh bin.  Check to see if the
    ! particle starts inside the mesh.
    if (any(ijk0(:n) < 1) .or. any(ijk0(:n) > m % dimension)) then
      ! The particle does not start in the mesh.  Note that we nudged the
      ! start and end coordinates by a TINY_BIT each so we will have
      ! difficulty resolving tracks that are less than 2*TINY_BIT in length.
      ! If the track is that short, it is also insignificant so we can
      ! safely ignore it in the tallies.
      if (total_distance < 2*TINY_BIT) return

      ! The particle does not start in the mesh so keep iterating the ijk0
      ! indices to cross the nearest mesh surface until we've found a valid
      ! bin.  MAX_SEARCH_ITER prevents an infinite loop.
      search_iter = 0
      do while (any(ijk0(:n) < 1) .or. any(ijk0(:n) > m % dimension))
        if (search_iter == MAX_SEARCH_ITER) then
          call warning("Failed to find a mesh intersection on a tally mesh &
               &filter.")
          return
        end if

        do j = 1, n
          if (abs(uvw(j)) < FP_PRECISION) then
            d(j) = INFINITY
          else if (uvw(j) > 0) then
            xyz_cross = m % lower_left(j) + ijk0(j) * m % width(j)
            d(j) = (xyz_cross - xyz0(j)) / uvw(j)
          else
            xyz_cross = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
            d(j) = (xyz_cross - xyz0(j)) / uvw(j)
          end if
        end do
        j = minloc(d(:n), 1)
        if (uvw(j) > ZERO) then
          ijk0(j) = ijk0(j) + 1
        else
          ijk0(j) = ijk0(j) - 1
        end if

        search_iter = search_iter + 1
      end do
      distance = d(j)
      xyz0 = xyz0 + distance * uvw
    end if

    do
      ! ========================================================================
      ! Compute the length of the track segment in the appropiate mesh cell and
      ! return.

      if (all(ijk0(:n) == ijk1(:n))) then
        ! The track ends in this cell.  Use the particle end location rather
        ! than the mesh surface.
        distance = sqrt(sum((xyz1 - xyz0)**2))
      else
        ! The track exits this cell.  Determine the distance to the closest mesh
        ! surface.
        do j = 1, n
          if (abs(uvw(j)) < FP_PRECISION) then
            d(j) = INFINITY
          else if (uvw(j) > 0) then
            xyz_cross = m % lower_left(j) + ijk0(j) * m % width(j)
            d(j) = (xyz_cross - xyz0(j)) / uvw(j)
          else
            xyz_cross = m % lower_left(j) + (ijk0(j) - 1) * m % width(j)
            d(j) = (xyz_cross - xyz0(j)) / uvw(j)
          end if
        end do
        j = minloc(d(:n), 1)
        distance = d(j)
      end if

      ! Assign the next tally bin and the score.
      bin = m % get_bin_from_indices(ijk0(:n))
      call match % bins % push_back(bin)
      call match % weights % push_back(distance / total_distance)

      ! Find the next mesh cell that the particle enters.

      ! If the particle track ends in that bin, then we are done.
      if (all(ijk0(:n) == ijk1(:n))) exit

      ! Translate the starting coordintes by the distance to that face. This
      ! should be the xyz that we computed the distance to in the last
      ! iteration of the filter loop.
      xyz0 = xyz0 + distance * uvw

      ! Increment the indices into the next mesh cell.
      if (uvw(j) > ZERO) then
        ijk0(j) = ijk0(j) + 1
      else
        ijk0(j) = ijk0(j) - 1
      end if

      ! If the next indices are invalid, then the track has left the mesh and
      ! we are done.
      if (any(ijk0(:n) < 1) .or. any(ijk0(:n) > m % dimension)) exit
    end do

  end subroutine get_all_bins_mesh

  subroutine to_statepoint_mesh(this, filter_group)
    class(MeshFilter), intent(in) :: this
    integer(HID_T),    intent(in) :: filter_group

    call write_dataset(filter_group, "type", "mesh")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", meshes(this % mesh) % id)
  end subroutine to_statepoint_mesh

  function text_label_mesh(this, bin) result(label)
    class(MeshFilter), intent(in) :: this
    integer,           intent(in) :: bin
    character(MAX_LINE_LEN)       :: label

    integer, allocatable       :: ijk(:)

    associate (m => meshes(this % mesh))
      allocate(ijk(m % n_dimension))
      call m % get_indices_from_bin(bin, ijk)
      if (m % n_dimension == 1) then
        label = "Mesh Index (" // trim(to_str(ijk(1))) // ")"
      elseif (m % n_dimension == 2) then
        label = "Mesh Index (" // trim(to_str(ijk(1))) // ", " // &
             trim(to_str(ijk(2))) // ")"
      elseif (m % n_dimension == 3) then
        label = "Mesh Index (" // trim(to_str(ijk(1))) // ", " // &
             trim(to_str(ijk(2))) // ", " // trim(to_str(ijk(3))) // ")"
      end if
    end associate
  end function text_label_mesh

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_mesh_filter_set_mesh(index, index_mesh) result(err) bind(C)
    ! Set the mesh for a mesh filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: index_mesh
    integer(C_INT) :: err

    err = 0
    if (index >= 1 .and. index <= n_filters) then
      if (allocated(filters(index) % obj)) then
        select type (f => filters(index) % obj)
        type is (MeshFilter)
          if (index_mesh >= 1 .and. index_mesh <= n_meshes) then
            f % mesh = index_mesh
            f % n_bins = product(meshes(index_mesh) % dimension)
          else
            err = E_OUT_OF_BOUNDS
            call set_errmsg("Index in 'meshes' array is out of bounds.")
          end if
        class default
          err = E_INVALID_TYPE
          call set_errmsg("Tried to set mesh on a non-mesh filter.")
        end select
      else
        err = E_ALLOCATE
        call set_errmsg("Filter type has not been set yet.")
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in filters array out of bounds.")
    end if
  end function openmc_mesh_filter_set_mesh

end module tally_filter_mesh
