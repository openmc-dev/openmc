module tally_filter_meshsurface

  use, intrinsic :: ISO_C_BINDING

  use constants
  use error
  use mesh_header
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

  type, public, extends(CppTallyFilter) :: MeshSurfaceFilter
    integer :: mesh
  contains
    procedure :: from_xml
    procedure :: to_statepoint
    procedure :: text_label
  end type MeshSurfaceFilter

contains

  subroutine from_xml(this, node)
    class(MeshSurfaceFilter), intent(inout) :: this
    type(XMLNode), intent(in) :: node

    integer :: id
    integer :: n
    integer(C_INT) :: err

    call this % from_xml_cpp_inner(node)

    n = node_word_count(node, "bins")

    if (n /= 1) call fatal_error("Only one mesh can be &
         &specified per meshsurface filter.")

    ! Determine id of mesh
    call get_node_value(node, "bins", id)

    ! Get pointer to mesh
    err = openmc_get_mesh_index(id, this % mesh)
    if (err /= 0) then
      call fatal_error("Could not find mesh " // trim(to_str(id)) &
           // " specified on filter.")
    end if

    ! Determine number of bins
    this % n_bins = this % n_bins_cpp()
  end subroutine from_xml

  subroutine to_statepoint(this, filter_group)
    class(MeshSurfaceFilter), intent(in) :: this
    integer(HID_T),           intent(in) :: filter_group

    type(RegularMesh) :: m

    m = meshes(this % mesh)
    call write_dataset(filter_group, "type", "meshsurface")
    call write_dataset(filter_group, "n_bins", this % n_bins)
    call write_dataset(filter_group, "bins", m % id())
  end subroutine to_statepoint

  function text_label(this, bin) result(label)
    class(MeshSurfaceFilter), intent(in) :: this
    integer,                  intent(in) :: bin
    character(MAX_LINE_LEN)              :: label

    integer :: i_mesh
    integer :: i_surf
    integer :: n_dim
    integer, allocatable :: ijk(:)
    type(RegularMesh) :: m

    m = meshes(this % mesh)
    n_dim = m % n_dimension()
    allocate(ijk(n_dim))

    ! Get flattend mesh index and surface index
    i_mesh = (bin - 1) / (4*n_dim) + 1
    i_surf = mod(bin - 1, 4*n_dim) + 1

    ! Get mesh index part of label
    call m % get_indices_from_bin(i_mesh, ijk)
    if (m % n_dimension() == 1) then
      label = "Mesh Index (" // trim(to_str(ijk(1))) // ")"
    elseif (m % n_dimension() == 2) then
      label = "Mesh Index (" // trim(to_str(ijk(1))) // ", " // &
            trim(to_str(ijk(2))) // ")"
    elseif (m % n_dimension() == 3) then
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
  end function text_label

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_meshsurface_filter_get_mesh(index, index_mesh) result(err) bind(C)
    ! Get the mesh for a mesh surface filter
    integer(C_INT32_T), value, intent(in) :: index
    integer(C_INT32_T), intent(out)       :: index_mesh
    integer(C_INT) :: err

    interface
      function meshsurface_filter_get_mesh(filt) result(index_mesh) bind(C)
        import C_PTR, C_INT
        type(C_PTR), value :: filt
        integer(C_INT)     :: index_mesh
      end function meshsurface_filter_get_mesh
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MeshSurfaceFilter)
        index_mesh = f % mesh
        index_mesh = meshsurface_filter_get_mesh(f % ptr)
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

    interface
      subroutine meshsurface_filter_set_mesh(filt, mesh) bind(C)
        import C_PTR, C_INT
        type(C_PTR),    value :: filt
        integer(C_INT), value :: mesh
      end subroutine meshsurface_filter_set_mesh
    end interface

    err = verify_filter(index)
    if (err == 0) then
      select type (f => filters(index) % obj)
      type is (MeshSurfaceFilter)
        if (index_mesh >= 0 .and. index_mesh < n_meshes()) then
          call meshsurface_filter_set_mesh(f % ptr, index_mesh)
          f % n_bins = f % n_bins_cpp()
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
