module tally_filter_meshsurface

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
  public :: openmc_meshsurface_filter_get_mesh
  public :: openmc_meshsurface_filter_set_mesh

!===============================================================================
! MESHFILTER indexes the location of particle events to a regular mesh.  For
! tracklength tallies, it will produce multiple valid bins and the bin weight
! will correspond to the fraction of the track length that lies in that bin.
!===============================================================================

  type, public, extends(TallyFilter) :: MeshSurfaceFilter
    integer :: mesh
  contains
    procedure :: from_xml
    procedure :: get_all_bins
    procedure :: to_statepoint
    procedure :: text_label
  end type MeshSurfaceFilter

contains

  subroutine from_xml(this, node)
    class(MeshSurfaceFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: i_mesh
    integer :: id
    integer :: n
    integer :: n_dim
    integer :: val

    n = node_word_count(node, "bins")

    if (n /= 1) call fatal_error("Only one mesh can be &
         &specified per meshsurface filter.")

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
    n_dim = meshes(i_mesh) % n_dimension
    this % n_bins = 4*n_dim*product(meshes(i_mesh) % dimension)

    ! Store the index of the mesh
    this % mesh = i_mesh
  end subroutine from_xml

  subroutine get_all_bins(this, p, estimator, match)
    class(MeshSurfaceFilter), intent(in)  :: this
    type(Particle),           intent(in)  :: p
    integer,                  intent(in)  :: estimator
    type(TallyFilterMatch), intent(inout) :: match

    integer :: j                    ! loop indices
    integer :: n_dim                ! num dimensions of the mesh
    integer :: d1                   ! dimension index
    integer :: ijk0(3)              ! indices of starting coordinates
    integer :: ijk1(3)              ! indices of ending coordinates
    integer :: n_cross              ! number of surface crossings
    integer :: i_mesh               ! flattened mesh bin index
    integer :: i_surf               ! surface index (1--12)
    integer :: i_bin                ! actual index for filter
    real(8) :: uvw(3)               ! cosine of angle of particle
    real(8) :: xyz0(3)              ! starting/intermediate coordinates
    real(8) :: xyz1(3)              ! ending coordinates of particle
    real(8) :: xyz_cross(3)         ! coordinates of bounding surfaces
    real(8) :: d(3)                 ! distance to each bounding surface
    real(8) :: distance             ! actual distance traveled
    logical :: start_in_mesh        ! particle's starting xyz in mesh?
    logical :: end_in_mesh          ! particle's ending xyz in mesh?

    ! Copy starting and ending location of particle
    xyz0 = p % last_xyz_current
    xyz1 = p % coord(1) % xyz

    associate (m => meshes(this % mesh))
      n_dim = m % n_dimension

      ! Determine indices for starting and ending location
      call m % get_indices(xyz0, ijk0, start_in_mesh)
      call m % get_indices(xyz1, ijk1, end_in_mesh)

      ! Check to see if start or end is in mesh -- if not, check if track still
      ! intersects with mesh
      if ((.not. start_in_mesh) .and. (.not. end_in_mesh)) then
        if (.not. m % intersects(xyz0, xyz1)) return
      end if

      ! Calculate number of surface crossings
      n_cross = sum(abs(ijk1(:n_dim) - ijk0(:n_dim)))
      if (n_cross == 0) return

      ! Copy particle's direction
      uvw = p % coord(1) % uvw

      ! Bounding coordinates
      do d1 = 1, n_dim
        if (uvw(d1) > 0) then
          xyz_cross(d1) = m % lower_left(d1) + ijk0(d1) * m % width(d1)
        else
          xyz_cross(d1) = m % lower_left(d1) + (ijk0(d1) - 1) * m % width(d1)
        end if
      end do

      do j = 1, n_cross
        ! Set the distances to infinity
        d = INFINITY

        ! Calculate distance to each bounding surface. We need to treat
        ! special case where the cosine of the angle is zero since this would
        ! result in a divide-by-zero.
        do d1 = 1, n_dim
          if (uvw(d1) == 0) then
            d(d1) = INFINITY
          else
            d(d1) = (xyz_cross(d1) - xyz0(d1))/uvw(d1)
          end if
        end do

        ! Determine the closest bounding surface of the mesh cell by
        ! calculating the minimum distance. Then use the minimum distance and
        ! direction of the particle to determine which surface was crossed.
        distance = minval(d)

        ! Loop over the dimensions
        do d1 = 1, n_dim

          ! Check whether distance is the shortest distance
          if (distance == d(d1)) then

            ! Check whether particle is moving in positive d1 direction
            if (uvw(d1) > 0) then

              ! Outward current on d1 max surface
              if (all(ijk0(:n_dim) >= 1) .and. &
                   all(ijk0(:n_dim) <= m % dimension)) then
                i_surf = d1 * 4 - 1
                i_mesh = m % get_bin_from_indices(ijk0)
                i_bin = 4*n_dim*(i_mesh - 1) + i_surf

                call match % bins % push_back(i_bin)
                call match % weights % push_back(ONE)
              end if

              ! Advance position
              ijk0(d1) = ijk0(d1) + 1
              xyz_cross(d1) = xyz_cross(d1) + m % width(d1)

              ! If the particle crossed the surface, tally the inward current on
              ! d1 min surface
              if (all(ijk0(:n_dim) >= 1) .and. &
                   all(ijk0(:n_dim) <= m % dimension)) then
                i_surf = d1 * 4 - 2
                i_mesh = m % get_bin_from_indices(ijk0)
                i_bin = 4*n_dim*(i_mesh - 1) + i_surf

                call match % bins % push_back(i_bin)
                call match % weights % push_back(ONE)
              end if

            else
              ! The particle is moving in the negative d1 direction

              ! Outward current on d1 min surface
              if (all(ijk0(:n_dim) >= 1) .and. &
                   all(ijk0(:n_dim) <= m % dimension)) then
                i_surf = d1 * 4 - 3
                i_mesh = m % get_bin_from_indices(ijk0)
                i_bin = 4*n_dim*(i_mesh - 1) + i_surf

                call match % bins % push_back(i_bin)
                call match % weights % push_back(ONE)
              end if

              ! Advance position
              ijk0(d1) = ijk0(d1) - 1
              xyz_cross(d1) = xyz_cross(d1) - m % width(d1)

              ! If the particle crossed the surface, tally the inward current on
              ! d1 max surface
              if (all(ijk0(:n_dim) >= 1) .and. &
                   all(ijk0(:n_dim) <= m % dimension)) then
                i_surf = d1 * 4
                i_mesh = m % get_bin_from_indices(ijk0)
                i_bin = 4*n_dim*(i_mesh - 1) + i_surf

                call match % bins % push_back(i_bin)
                call match % weights % push_back(ONE)
              end if
            end if
          end if
        end do

        ! Calculate new coordinates
        xyz0 = xyz0 + distance * uvw
      end do
    end associate

  end subroutine get_all_bins

  subroutine to_statepoint(this, filter_group)
    class(MeshSurfaceFilter), intent(in) :: this
    integer(HID_T),           intent(in) :: filter_group

    call write_dataset(filter_group, "type", "meshsurface")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", meshes(this % mesh) % id)
  end subroutine to_statepoint

  function text_label(this, bin) result(label)
    class(MeshSurfaceFilter), intent(in) :: this
    integer,                  intent(in) :: bin
    character(MAX_LINE_LEN)              :: label

    integer :: i_mesh
    integer :: i_surf
    integer :: n_dim
    integer, allocatable :: ijk(:)

    associate (m => meshes(this % mesh))
      n_dim = m % n_dimension
      allocate(ijk(n_dim))

      ! Get flattend mesh index and surface index
      i_mesh = (bin - 1) / (4*n_dim) + 1
      i_surf = mod(bin - 1, 4*n_dim) + 1

      ! Get mesh index part of label
      call m % get_indices_from_bin(i_mesh, ijk)
      if (m % n_dimension == 1) then
        label = "Mesh Index (" // trim(to_str(ijk(1))) // ")"
      elseif (m % n_dimension == 2) then
        label = "Mesh Index (" // trim(to_str(ijk(1))) // ", " // &
             trim(to_str(ijk(2))) // ")"
      elseif (m % n_dimension == 3) then
        label = "Mesh Index (" // trim(to_str(ijk(1))) // ", " // &
             trim(to_str(ijk(2))) // ", " // trim(to_str(ijk(3))) // ")"
      end if

      ! Get surface part of label
      select case (i_surf)
      case (OUT_LEFT)
        label = trim(label) // " Outgoing, x-min"
      case (IN_LEFT)
        label = trim(label) // " Incoming, x-min"
      case (OUT_RIGHT)
        label = trim(label) // " Outgoing, x-max"
      case (IN_RIGHT)
        label = trim(label) // " Incoming, x-max"
      case (OUT_BACK)
        label = trim(label) // " Outgoing, y-min"
      case (IN_BACK)
        label = trim(label) // " Incoming, y-min"
      case (OUT_FRONT)
        label = trim(label) // " Outgoing, y-max"
      case (IN_FRONT)
        label = trim(label) // " Incoming, y-max"
      case (OUT_BOTTOM)
        label = trim(label) // " Outgoing, z-min"
      case (IN_BOTTOM)
        label = trim(label) // " Incoming, z-min"
      case (OUT_TOP)
        label = trim(label) // " Outgoing, z-max"
      case (IN_TOP)
        label = trim(label) // " Incoming, z-max"
      end select
    end associate
  end function text_label

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_meshsurface_filter_get_mesh(index, index_mesh) result(err) bind(C)
    ! Get the mesh for a mesh surface filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), intent(out)       :: index_mesh
    integer(C_INT) :: err

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MeshSurfaceFilter)
        index_mesh = f % mesh
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to set mesh on a non-mesh filter.")
      end select
    end if
  end function openmc_meshsurface_filter_get_mesh


  function openmc_meshsurface_filter_set_mesh(index, index_mesh) result(err) bind(C)
    ! Set the mesh for a mesh surface filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), value, intent(in) :: index_mesh
    integer(C_INT) :: err

    integer :: n_dim

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MeshSurfaceFilter)
        if (index_mesh >= 1 .and. index_mesh <= n_meshes) then
          f % mesh = index_mesh
          n_dim = meshes(index_mesh) % n_dimension
          f % n_bins = 4*n_dim*product(meshes(index_mesh) % dimension)
        else
          err = E_OUT_OF_BOUNDS
          call set_errmsg("Index in 'meshes' array is out of bounds.")
        end if
      class default
        err = E_INVALID_TYPE
        call set_errmsg("Tried to set mesh on a non-mesh filter.")
      end select
    end if
  end function openmc_meshsurface_filter_set_mesh

end module tally_filter_meshsurface
