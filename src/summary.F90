module summary

  use hdf5

  use constants
  use endf,            only: reaction_name
  use error,           only: write_message
  use geometry_header
  use hdf5_interface
  use material_header, only: Material, n_materials
  use mesh_header,     only: RegularMesh
  use message_passing
  use mgxs_header,     only: nuclides_MG
  use nuclide_header
  use output,          only: time_stamp
  use settings,        only: run_CE
  use surface_header
  use string,          only: to_str
  use tally_header,    only: TallyObject
  use tally_filter_distribcell, only: find_offset

  implicit none
  private

  public :: write_summary

contains

!===============================================================================
! WRITE_SUMMARY
!===============================================================================

  subroutine write_summary()

    integer(HID_T) :: file_id

    ! Display output message
    call write_message("Writing summary.h5 file...", 5)

    ! Create a new file using default properties.
    file_id = file_create("summary.h5")

    call write_header(file_id)
    call write_nuclides(file_id)
    call write_geometry(file_id)
    call write_materials(file_id)

    ! Terminate access to the file.
    call file_close(file_id)

  end subroutine write_summary

!===============================================================================
! WRITE_HEADER
!===============================================================================

  subroutine write_header(file_id)
    integer(HID_T), intent(in) :: file_id

    ! Write filetype and version info
    call write_attribute(file_id, "filetype", "summary")
    call write_attribute(file_id, "version", VERSION_SUMMARY)
    call write_attribute(file_id, "openmc_version", VERSION)
#ifdef GIT_SHA1
    call write_attribute(file_id, "git_sha1", GIT_SHA1)
#endif

    ! Write current date and time
    call write_attribute(file_id, "date_and_time", time_stamp())

  end subroutine write_header

!===============================================================================
! WRITE_NUCLIDES
!===============================================================================

  subroutine write_nuclides(file_id)
    integer(HID_T), intent(in) :: file_id
    integer(HID_T) :: nuclide_group
    integer :: i
    character(12), allocatable :: nucnames(:)
    real(8), allocatable :: awrs(:)

    ! Write useful data from nuclide objects
    nuclide_group = create_group(file_id, "nuclides")
    call write_attribute(nuclide_group, "n_nuclides", n_nuclides)

    ! Build array of nuclide names and awrs
    allocate(nucnames(n_nuclides))
    allocate(awrs(n_nuclides))
    do i = 1, n_nuclides
      if (run_CE) then
        nucnames(i) = nuclides(i) % name
        awrs(i)     = nuclides(i) % awr
      else
        nucnames(i) = nuclides_MG(i) % obj % name
        awrs(i)     = nuclides_MG(i) % obj % awr
      end if
    end do

    ! Write nuclide names and awrs
    call write_dataset(nuclide_group, "names", nucnames)
    call write_dataset(nuclide_group, "awrs", awrs)

    call close_group(nuclide_group)

    deallocate(nucnames, awrs)

  end subroutine write_nuclides

!===============================================================================
! WRITE_GEOMETRY
!===============================================================================

  subroutine write_geometry(file_id)
    integer(HID_T), intent(in) :: file_id

    integer :: i, j, k, m
    integer, allocatable :: lattice_universes(:,:,:)
    integer, allocatable :: cell_materials(:)
    integer, allocatable :: cell_ids(:)
    real(8), allocatable :: cell_temperatures(:)
    integer(HID_T) :: geom_group
    integer(HID_T) :: cells_group, cell_group
    integer(HID_T) :: surfaces_group, surface_group
    integer(HID_T) :: universes_group, univ_group
    integer(HID_T) :: lattices_group, lattice_group
    real(8), allocatable :: coeffs(:)
    character(:), allocatable :: region_spec
    type(Cell),     pointer :: c
    class(Surface), pointer :: s
    class(Lattice), pointer :: lat

    ! Use H5LT interface to write number of geometry objects
    geom_group = create_group(file_id, "geometry")
    call write_attribute(geom_group, "n_cells", n_cells)
    call write_attribute(geom_group, "n_surfaces", n_surfaces)
    call write_attribute(geom_group, "n_universes", n_universes)
    call write_attribute(geom_group, "n_lattices", n_lattices)

    ! ==========================================================================
    ! WRITE INFORMATION ON CELLS

    ! Create a cell group (nothing directly written in this group) then close
    cells_group = create_group(geom_group, "cells")

    ! Write information on each cell
    CELL_LOOP: do i = 1, n_cells
      c => cells(i)
      cell_group = create_group(cells_group, "cell " // trim(to_str(c%id)))

      ! Write name for this cell
      call write_dataset(cell_group, "name", c%name)

      ! Write universe for this cell
      call write_dataset(cell_group, "universe", universes(c%universe)%id)

      ! Write information on what fills this cell
      select case (c%type)
      case (FILL_MATERIAL)
        call write_dataset(cell_group, "fill_type", "material")

        if (size(c % material) == 1) then
          if (c % material(1) == MATERIAL_VOID) then
            call write_dataset(cell_group, "material", MATERIAL_VOID)
          else
            call write_dataset(cell_group, "material", &
                 materials(c % material(1)) % id)
          end if
        else
          allocate(cell_materials(size(c % material)))
          do j = 1, size(c % material)
            if (c % material(j) == MATERIAL_VOID) then
              cell_materials(j) = MATERIAL_VOID
            else
              cell_materials(j) = materials(c % material(j)) % id
            end if
          end do
          call write_dataset(cell_group, "material", cell_materials)
          deallocate(cell_materials)
        end if

        allocate(cell_temperatures(size(c % sqrtkT)))
        cell_temperatures(:) = c % sqrtkT(:)
        cell_temperatures(:) = cell_temperatures(:)**2 / K_BOLTZMANN
        call write_dataset(cell_group, "temperature", cell_temperatures)
        deallocate(cell_temperatures)

      case (FILL_UNIVERSE)
        call write_dataset(cell_group, "fill_type", "universe")
        call write_dataset(cell_group, "fill", universes(c%fill)%id)
        if (allocated(c%offset)) then
          if (size(c%offset) > 0) then
            call write_dataset(cell_group, "offset", c%offset)
          end if
        end if

        if (allocated(c%translation)) then
          call write_dataset(cell_group, "translation", c%translation)
        end if
        if (allocated(c%rotation)) then
          call write_dataset(cell_group, "rotation", c%rotation)
        end if

      case (FILL_LATTICE)
        call write_dataset(cell_group, "fill_type", "lattice")
        call write_dataset(cell_group, "lattice", lattices(c%fill)%obj%id)
      end select

      ! Write list of bounding surfaces
      region_spec = ""
      do j = 1, size(c%region)
        k = c%region(j)
        select case(k)
        case (OP_LEFT_PAREN)
          region_spec = trim(region_spec) // " ("
        case (OP_RIGHT_PAREN)
          region_spec = trim(region_spec) // " )"
        case (OP_COMPLEMENT)
          region_spec = trim(region_spec) // " ~"
        case (OP_INTERSECTION)
        case (OP_UNION)
          region_spec = trim(region_spec) // " |"
        case default
          region_spec = trim(region_spec) // " " // to_str(&
               sign(surfaces(abs(k))%obj%id, k))
        end select
      end do
      call write_dataset(cell_group, "region", adjustl(region_spec))

      call close_group(cell_group)
    end do CELL_LOOP

    call close_group(cells_group)

    ! ==========================================================================
    ! WRITE INFORMATION ON SURFACES

    ! Create surfaces group
    surfaces_group = create_group(geom_group, "surfaces")

    ! Write information on each surface
    SURFACE_LOOP: do i = 1, n_surfaces
      s => surfaces(i)%obj
      surface_group = create_group(surfaces_group, "surface " // &
           trim(to_str(s%id)))

      ! Write name for this surface
      call write_dataset(surface_group, "name", s%name)

      ! Write surface type
      select type (s)
      type is (SurfaceXPlane)
        call write_dataset(surface_group, "type", "x-plane")
        allocate(coeffs(1))
        coeffs(1) = s%x0

      type is (SurfaceYPlane)
        call write_dataset(surface_group, "type", "y-plane")
        allocate(coeffs(1))
        coeffs(1) = s%y0

      type is (SurfaceZPlane)
        call write_dataset(surface_group, "type", "z-plane")
        allocate(coeffs(1))
        coeffs(1) = s%z0

      type is (SurfacePlane)
        call write_dataset(surface_group, "type", "plane")
        allocate(coeffs(4))
        coeffs(:) = [s%A, s%B, s%C, s%D]

      type is (SurfaceXCylinder)
        call write_dataset(surface_group, "type", "x-cylinder")
        allocate(coeffs(3))
        coeffs(:) = [s%y0, s%z0, s%r]

      type is (SurfaceYCylinder)
        call write_dataset(surface_group, "type", "y-cylinder")
        allocate(coeffs(3))
        coeffs(:) = [s%x0, s%z0, s%r]

      type is (SurfaceZCylinder)
        call write_dataset(surface_group, "type", "z-cylinder")
        allocate(coeffs(3))
        coeffs(:) = [s%x0, s%y0, s%r]

      type is (SurfaceSphere)
        call write_dataset(surface_group, "type", "sphere")
        allocate(coeffs(4))
        coeffs(:) = [s%x0, s%y0, s%z0, s%r]

      type is (SurfaceXCone)
        call write_dataset(surface_group, "type", "x-cone")
        allocate(coeffs(4))
        coeffs(:) = [s%x0, s%y0, s%z0, s%r2]

      type is (SurfaceYCone)
        call write_dataset(surface_group, "type", "y-cone")
        allocate(coeffs(4))
        coeffs(:) = [s%x0, s%y0, s%z0, s%r2]

      type is (SurfaceZCone)
        call write_dataset(surface_group, "type", "z-cone")
        allocate(coeffs(4))
        coeffs(:) = [s%x0, s%y0, s%z0, s%r2]

      type is (SurfaceQuadric)
        call write_dataset(surface_group, "type", "quadric")
        allocate(coeffs(10))
        coeffs(:) = [s%A, s%B, s%C, s%D, s%E, s%F, s%G, s%H, s%J, s%K]

      end select
      call write_dataset(surface_group, "coefficients", coeffs)
      deallocate(coeffs)

      ! Write boundary type
      select case (s%bc)
      case (BC_TRANSMIT)
        call write_dataset(surface_group, "boundary_type", "transmission")
      case (BC_VACUUM)
        call write_dataset(surface_group, "boundary_type", "vacuum")
      case (BC_REFLECT)
        call write_dataset(surface_group, "boundary_type", "reflective")
      case (BC_PERIODIC)
        call write_dataset(surface_group, "boundary_type", "periodic")
      end select

      call close_group(surface_group)
    end do SURFACE_LOOP

    call close_group(surfaces_group)

    ! ==========================================================================
    ! WRITE INFORMATION ON UNIVERSES

    ! Create universes group (nothing directly written here) then close
    universes_group = create_group(geom_group, "universes")

    ! Write information on each universe
    UNIVERSE_LOOP: do i = 1, n_universes
      associate (u => universes(i))
        univ_group = create_group(universes_group, "universe " // &
             trim(to_str(u%id)))

        ! Write list of cells in this universe
        if (size(u % cells) > 0) then
          allocate(cell_ids(size(u % cells)))
          do j = 1, size(u % cells)
            cell_ids(j) = cells(u % cells(j)) % id
          end do
          call write_dataset(univ_group, "cells", cell_ids)
          deallocate(cell_ids)
        end if

        call close_group(univ_group)
      end associate
    end do UNIVERSE_LOOP

    call close_group(universes_group)

    ! ==========================================================================
    ! WRITE INFORMATION ON LATTICES

    ! Create lattices group (nothing directly written here) then close
    lattices_group = create_group(geom_group, "lattices")

    ! Write information on each lattice
    LATTICE_LOOP: do i = 1, n_lattices
      lat => lattices(i)%obj
      lattice_group = create_group(lattices_group, "lattice " // trim(to_str(lat%id)))

      ! Write name, pitch, and outer universe
      call write_dataset(lattice_group, "name", lat%name)
      call write_dataset(lattice_group, "pitch", lat%pitch)
      if (lat % outer > 0) then
        call write_dataset(lattice_group, "outer", universes(lat % outer) % id)
      else
        call write_dataset(lattice_group, "outer", lat % outer)
      end if

      select type (lat)
      type is (RectLattice)
        ! Write lattice type.
        call write_dataset(lattice_group, "type", "rectangular")

        ! Write lattice dimensions, lower left corner, and pitch
        if (lat % is_3d) then
          call write_dataset(lattice_group, "dimension", lat % n_cells)
          call write_dataset(lattice_group, "lower_left", lat % lower_left)
        else
          call write_dataset(lattice_group, "dimension", lat % n_cells(1:2))
          call write_dataset(lattice_group, "lower_left", lat % lower_left)
        end if

        ! Write lattice universes.
        allocate(lattice_universes(lat%n_cells(1), lat%n_cells(2), &
             &lat%n_cells(3)))
        do j = 1, lat%n_cells(1)
          do k = 0, lat%n_cells(2) - 1
            do m = 1, lat%n_cells(3)
              lattice_universes(j, k+1, m) = &
                   universes(lat%universes(j, lat%n_cells(2) - k, m))%id
            end do
          end do
        end do

      type is (HexLattice)
        ! Write lattice type.
        call write_dataset(lattice_group, "type", "hexagonal")

        ! Write number of lattice cells.
        call write_dataset(lattice_group, "n_rings", lat%n_rings)
        call write_dataset(lattice_group, "n_axial", lat%n_axial)

        ! Write lattice center
        call write_dataset(lattice_group, "center", lat%center)

        ! Write lattice universes.
        allocate(lattice_universes(2*lat%n_rings - 1, 2*lat%n_rings - 1, &
             &lat%n_axial))
        do m = 1, lat%n_axial
          do k = 1, 2*lat%n_rings - 1
            do j = 1, 2*lat%n_rings - 1
              if (j + k < lat%n_rings + 1) then
                ! This array position is never used; put a -1 to indicate this
                lattice_universes(j,k,m) = -1
                cycle
              else if (j + k > 3*lat%n_rings - 1) then
                ! This array position is never used; put a -1 to indicate this
                lattice_universes(j,k,m) = -1
                cycle
              end if
              lattice_universes(j,k,m) = universes(lat%universes(j,k,m))%id
            end do
          end do
        end do
      end select

      ! Write lattice universes
      call write_dataset(lattice_group, "universes", lattice_universes)
      deallocate(lattice_universes)

      call close_group(lattice_group)
    end do LATTICE_LOOP

    call close_group(lattices_group)
    call close_group(geom_group)

  end subroutine write_geometry

!===============================================================================
! WRITE_MATERIALS
!===============================================================================

  subroutine write_materials(file_id)
    integer(HID_T), intent(in) :: file_id

    integer :: i
    integer :: j
    character(20), allocatable :: nucnames(:)
    integer(HID_T) :: materials_group
    integer(HID_T) :: material_group
    type(Material), pointer :: m

    materials_group = create_group(file_id, "materials")

    ! write number of materials
    call write_dataset(file_id, "n_materials", n_materials)

    ! Write information on each material
    do i = 1, n_materials
      m => materials(i)
      material_group = create_group(materials_group, "material " // &
           trim(to_str(m%id)))

      if (m % depletable) then
        call write_attribute(material_group, "depletable", 1)
      else
        call write_attribute(material_group, "depletable", 0)
      end if

      ! Write name for this material
      call write_dataset(material_group, "name", m % name)

      ! Write atom density with units
      call write_dataset(material_group, "atom_density", m % density)

      ! Copy ZAID for each nuclide to temporary array
      allocate(nucnames(m%n_nuclides))
      do j = 1, m%n_nuclides
        if (run_CE) then
          nucnames(j) = nuclides(m%nuclide(j))%name
        else
          nucnames(j) = nuclides_MG(m%nuclide(j))%obj%name
        end if
      end do

      ! Write temporary array to 'nuclides'
      call write_dataset(material_group, "nuclides", nucnames)

      ! Deallocate temporary array
      deallocate(nucnames)

      ! Write atom densities
      call write_dataset(material_group, "nuclide_densities", m%atom_density)

      if (m%n_sab > 0) then
        call write_dataset(material_group, "sab_names", m%sab_names)
      end if

      call close_group(material_group)
    end do

    call close_group(materials_group)

  end subroutine write_materials

end module summary
