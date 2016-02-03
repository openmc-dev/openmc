module summary

  use ace_header,      only: Reaction, UrrData
  use constants
  use endf,            only: reaction_name
  use geometry_header, only: Cell, Universe, Lattice, RectLattice, &
                             &HexLattice
  use global
  use hdf5_interface
  use material_header, only: Material
  use mesh_header,     only: RegularMesh
  use nuclide_header
  use output,          only: time_stamp
  use surface_header
  use string,          only: to_str
  use tally_header,    only: TallyObject

  use hdf5

  implicit none
  private

  public :: write_summary

contains

!===============================================================================
! WRITE_SUMMARY
!===============================================================================

  subroutine write_summary()

    integer(HID_T) :: file_id

    ! Create a new file using default properties.
    file_id = file_create("summary.h5")

    ! Write header information
    call write_header(file_id)

    if (run_CE) then
      call write_dataset(file_id, "run_CE", 1)
    else
      call write_dataset(file_id, "run_CE", 0)
    end if

    ! Write number of particles
    call write_dataset(file_id, "n_particles", n_particles)
    call write_dataset(file_id, "n_batches", n_batches)
    call write_attribute_string(file_id, "n_particles", &
         "description", "Number of particles per generation")
    call write_attribute_string(file_id, "n_batches", &
         "description", "Total number of batches")

    ! Write eigenvalue information
    if (run_mode == MODE_EIGENVALUE) then
      ! write number of inactive/active batches and generations/batch
      call write_dataset(file_id, "n_inactive", n_inactive)
      call write_dataset(file_id, "n_active", n_active)
      call write_dataset(file_id, "gen_per_batch", gen_per_batch)

      ! Add description of each variable
      call write_attribute_string(file_id, "n_inactive", &
           "description", "Number of inactive batches")
      call write_attribute_string(file_id, "n_active", &
           "description", "Number of active batches")
      call write_attribute_string(file_id, "gen_per_batch", &
           "description", "Number of generations per batch")
    end if

    call write_geometry(file_id)
    call write_materials(file_id)
    if (n_tallies > 0) then
      call write_tallies(file_id)
    end if

    ! Terminate access to the file.
    call file_close(file_id)

  end subroutine write_summary

!===============================================================================
! WRITE_HEADER
!===============================================================================

  subroutine write_header(file_id)
    integer(HID_T), intent(in) :: file_id

    ! Write filetype and revision
    call write_dataset(file_id, "filetype", "summary")
    call write_dataset(file_id, "revision", REVISION_SUMMARY)

    ! Write version information
    call write_dataset(file_id, "version_major", VERSION_MAJOR)
    call write_dataset(file_id, "version_minor", VERSION_MINOR)
    call write_dataset(file_id, "version_release", VERSION_RELEASE)

    ! Write current date and time
    call write_dataset(file_id, "date_and_time", time_stamp())

    ! Write MPI information
    call write_dataset(file_id, "n_procs", n_procs)
    call write_attribute_string(file_id, "n_procs", "description", &
         "Number of MPI processes")

  end subroutine write_header

!===============================================================================
! WRITE_GEOMETRY
!===============================================================================

  subroutine write_geometry(file_id)
    integer(HID_T), intent(in) :: file_id

    integer          :: i, j, k, m
    integer, allocatable :: lattice_universes(:,:,:)
    integer, allocatable :: cell_materials(:)
    integer(HID_T) :: geom_group
    integer(HID_T) :: cells_group, cell_group
    integer(HID_T) :: surfaces_group, surface_group
    integer(HID_T) :: universes_group, univ_group
    integer(HID_T) :: lattices_group, lattice_group
    real(8), allocatable :: coeffs(:)
    character(REGION_SPEC_LEN) :: region_spec
    type(Cell),     pointer :: c
    class(Surface), pointer :: s
    type(Universe), pointer :: u
    class(Lattice), pointer :: lat

    ! Use H5LT interface to write number of geometry objects
    geom_group = create_group(file_id, "geometry")
    call write_dataset(geom_group, "n_cells", n_cells)
    call write_dataset(geom_group, "n_surfaces", n_surfaces)
    call write_dataset(geom_group, "n_universes", n_universes)
    call write_dataset(geom_group, "n_lattices", n_lattices)

    ! ==========================================================================
    ! WRITE INFORMATION ON CELLS

    ! Create a cell group (nothing directly written in this group) then close
    cells_group = create_group(geom_group, "cells")

    ! Write information on each cell
    CELL_LOOP: do i = 1, n_cells
      c => cells(i)
      cell_group = create_group(cells_group, "cell " // trim(to_str(c%id)))

      ! Write internal OpenMC index for this cell
      call write_dataset(cell_group, "index", i)

      ! Write name for this cell
      call write_dataset(cell_group, "name", c%name)

      ! Write universe for this cell
      call write_dataset(cell_group, "universe", universes(c%universe)%id)

      ! Write information on what fills this cell
      select case (c%type)
      case (CELL_NORMAL)
        call write_dataset(cell_group, "fill_type", "normal")
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

      case (CELL_FILL)
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

      case (CELL_LATTICE)
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

      call write_dataset(cell_group, "distribcell_index", c % distribcell_index)

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

      ! Write internal OpenMC index for this surface
      call write_dataset(surface_group, "index", i)

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

      ! Write boundary condition
      select case (s%bc)
      case (BC_TRANSMIT)
        call write_dataset(surface_group, "boundary_condition", "transmission")
      case (BC_VACUUM)
        call write_dataset(surface_group, "boundary_condition", "vacuum")
      case (BC_REFLECT)
        call write_dataset(surface_group, "boundary_condition", "reflective")
      case (BC_PERIODIC)
        call write_dataset(surface_group, "boundary_condition", "periodic")
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
      u => universes(i)
      univ_group = create_group(universes_group, "universe " // &
           trim(to_str(u%id)))

      ! Write internal OpenMC index for this universe
      call write_dataset(univ_group, "index", i)

      ! Write list of cells in this universe
      if (u%n_cells > 0) call write_dataset(univ_group, "cells", u%cells)

      call close_group(univ_group)
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

      ! Write internal OpenMC index for this lattice
      call write_dataset(lattice_group, "index", i)

      ! Write name, pitch, and outer universe
      call write_dataset(lattice_group, "name", lat%name)
      call write_dataset(lattice_group, "pitch", lat%pitch)
      call write_dataset(lattice_group, "outer", lat%outer)

      ! Write distribcell offsets if present
      if (allocated(lat%offset)) then
        if (size(lat%offset) > 0) then
          call write_dataset(lattice_group, "offsets", lat%offset)
        end if
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
    integer :: i_list
    character(12), allocatable :: nucnames(:)
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

      ! Write internal OpenMC index for this material
      call write_dataset(material_group, "index", i)

      ! Write name for this material
      call write_dataset(material_group, "name", m%name)

      ! Write atom density with units
      call write_dataset(material_group, "atom_density", m%density)
      call write_attribute_string(material_group, "atom_density", "units", &
           "atom/b-cm")

      ! Copy ZAID for each nuclide to temporary array
      allocate(nucnames(m%n_nuclides))
      do j = 1, m%n_nuclides
        if (run_CE) then
          i_list = nuclides(m%nuclide(j))%listing
        else
          i_list = nuclides_MG(m%nuclide(j))%obj%listing
        end if
        nucnames(j) = xs_listings(i_list)%alias
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

!===============================================================================
! WRITE_TALLIES
!===============================================================================

  subroutine write_tallies(file_id)
    integer(HID_T), intent(in) :: file_id

    integer :: i, j, k
    integer :: i_list, i_xs
    integer :: n_order      ! loop index for moment orders
    integer :: nm_order     ! loop index for Ynm moment orders
    integer(HID_T) :: tallies_group
    integer(HID_T) :: mesh_group
    integer(HID_T) :: tally_group
    integer(HID_T) :: filter_group
    character(20), allocatable :: str_array(:)
    type(RegularMesh), pointer :: m
    type(TallyObject), pointer :: t

    tallies_group = create_group(file_id, "tallies")

    ! Write total number of meshes
    call write_dataset(tallies_group, "n_meshes", n_meshes)

    ! Write information for meshes
    MESH_LOOP: do i = 1, n_meshes
      m => meshes(i)
      mesh_group = create_group(tallies_group, "mesh " // trim(to_str(m%id)))

      ! Write internal OpenMC index for this mesh
      call write_dataset(mesh_group, "index", i)

      ! Write type and number of dimensions
      call write_dataset(mesh_group, "type", "regular")

      ! Write mesh information
      call write_dataset(mesh_group, "dimension", m%dimension)
      call write_dataset(mesh_group, "lower_left", m%lower_left)
      call write_dataset(mesh_group, "upper_right", m%upper_right)
      call write_dataset(mesh_group, "width", m%width)

      call close_group(mesh_group)
    end do MESH_LOOP

    ! Write number of tallies
    call write_dataset(tallies_group, "n_tallies", n_tallies)

    TALLY_METADATA: do i = 1, n_tallies
      ! Get pointer to tally
      t => tallies(i)
      tally_group = create_group(tallies_group, "tally " // trim(to_str(t%id)))

      ! Write internal OpenMC index for this tally
      call write_dataset(tally_group, "index", i)

      ! Write the name for this tally
      call write_dataset(tally_group, "name", t%name)

      ! Write number of filters
      call write_dataset(tally_group, "n_filters", t%n_filters)

      FILTER_LOOP: do j = 1, t%n_filters
        filter_group = create_group(tally_group, "filter " // trim(to_str(j)))

        ! Write number of bins for this filter
        call write_dataset(filter_group, "n_bins", t%filters(j)%n_bins)

        ! Write filter bins
        if (t%filters(j)%type == FILTER_ENERGYIN .or. &
             t%filters(j)%type == FILTER_ENERGYOUT .or. &
             t%filters(j)%type == FILTER_MU .or. &
             t%filters(j)%type == FILTER_POLAR .or. &
             t%filters(j)%type == FILTER_AZIMUTHAL) then
          call write_dataset(filter_group, "bins", t%filters(j)%real_bins)
        else
          call write_dataset(filter_group, "bins", t%filters(j)%int_bins)
        end if

        ! Write name of type
        select case (t%filters(j)%type)
        case(FILTER_UNIVERSE)
          call write_dataset(filter_group, "type", "universe")
        case(FILTER_MATERIAL)
          call write_dataset(filter_group, "type", "material")
        case(FILTER_CELL)
          call write_dataset(filter_group, "type", "cell")
        case(FILTER_CELLBORN)
          call write_dataset(filter_group, "type", "cellborn")
        case(FILTER_SURFACE)
          call write_dataset(filter_group, "type", "surface")
        case(FILTER_MESH)
          call write_dataset(filter_group, "type", "mesh")
        case(FILTER_ENERGYIN)
          call write_dataset(filter_group, "type", "energy")
        case(FILTER_ENERGYOUT)
          call write_dataset(filter_group, "type", "energyout")
        case(FILTER_DISTRIBCELL)
          call write_dataset(filter_group, "type", "distribcell")
        case(FILTER_MU)
          call write_dataset(filter_group, "type", "mu")
        case(FILTER_POLAR)
          call write_dataset(filter_group, "type", "polar")
        case(FILTER_AZIMUTHAL)
          call write_dataset(filter_group, "type", "azimuthal")
        case(FILTER_DELAYEDGROUP)
          call write_dataset(filter_group, "type", "delayedgroup")
        end select

        call close_group(filter_group)
      end do FILTER_LOOP

      ! Create temporary array for nuclide bins
      allocate(str_array(t%n_nuclide_bins))
      NUCLIDE_LOOP: do j = 1, t%n_nuclide_bins
        if (t%nuclide_bins(j) > 0) then
          i_list = nuclides(t%nuclide_bins(j))%listing
          i_xs = index(xs_listings(i_list)%alias, '.')
          if (i_xs > 0) then
            str_array(j) = xs_listings(i_list)%alias(1:i_xs - 1)
          else
            str_array(j) = xs_listings(i_list)%alias
          end if
        else
          str_array(j) = 'total'
        end if
      end do NUCLIDE_LOOP

      ! Write and deallocate nuclide bins
      call write_dataset(tally_group, "nuclides", str_array)
      deallocate(str_array)

      ! Write number of score bins
      call write_dataset(tally_group, "n_score_bins", t%n_score_bins)
      allocate(str_array(size(t%score_bins)))
      do j = 1, size(t%score_bins)
        str_array(j) = reaction_name(t%score_bins(j))
      end do
      call write_dataset(tally_group, "score_bins", str_array)

      deallocate(str_array)

      ! Write explicit moment order strings for each score bin
      k = 1
      allocate(str_array(t%n_score_bins))
      MOMENT_LOOP: do j = 1, t%n_user_score_bins
        select case(t%score_bins(k))
        case (SCORE_SCATTER_N, SCORE_NU_SCATTER_N)
          str_array(k) = 'P' // trim(to_str(t%moment_order(k)))
          k = k + 1
        case (SCORE_SCATTER_PN, SCORE_NU_SCATTER_PN)
          do n_order = 0, t%moment_order(k)
            str_array(k) = 'P' // trim(to_str(n_order))
            k = k + 1
          end do
        case (SCORE_SCATTER_YN, SCORE_NU_SCATTER_YN, SCORE_FLUX_YN, &
             SCORE_TOTAL_YN)
          do n_order = 0, t%moment_order(k)
            do nm_order = -n_order, n_order
              str_array(k) = 'Y' // trim(to_str(n_order)) // ',' // &
                   trim(to_str(nm_order))
              k = k + 1
            end do
          end do
        case default
          str_array(k) = ''
          k = k + 1
        end select
      end do MOMENT_LOOP

      call write_dataset(tally_group, "moment_orders", str_array)
      deallocate(str_array)

      call close_group(tally_group)
    end do TALLY_METADATA

    call close_group(tallies_group)

  end subroutine write_tallies

!===============================================================================
! WRITE_TIMING
!===============================================================================

  subroutine write_timing(file_id)
    integer(HID_T), intent(in) :: file_id

    integer(8)     :: total_particles
    integer(HID_T) :: time_group
    real(8)        :: speed

    time_group = create_group(file_id, "timing")

    ! Write timing data
    call write_dataset(time_group, "time_initialize", time_initialize%elapsed)
    call write_dataset(time_group, "time_read_xs", time_read_xs%elapsed)
    call write_dataset(time_group, "time_transport", time_transport%elapsed)
    call write_dataset(time_group, "time_bank", time_bank%elapsed)
    call write_dataset(time_group, "time_bank_sample", time_bank_sample%elapsed)
    call write_dataset(time_group, "time_bank_sendrecv", time_bank_sendrecv%elapsed)
    call write_dataset(time_group, "time_tallies", time_tallies%elapsed)
    call write_dataset(time_group, "time_inactive", time_inactive%elapsed)
    call write_dataset(time_group, "time_active", time_active%elapsed)
    call write_dataset(time_group, "time_finalize", time_finalize%elapsed)
    call write_dataset(time_group, "time_total", time_total%elapsed)

    ! Add descriptions to timing data
    call write_attribute_string(time_group, "time_initialize", "description", &
         "Total time elapsed for initialization (s)")
    call write_attribute_string(time_group, "time_read_xs", "description", &
         "Time reading cross-section libraries (s)")
    call write_attribute_string(time_group, "time_transport", "description", &
         "Time in transport only (s)")
    call write_attribute_string(time_group, "time_bank", "description", &
         "Total time synchronizing fission bank (s)")
    call write_attribute_string(time_group, "time_bank_sample", "description", &
         "Time between generations sampling source sites (s)")
    call write_attribute_string(time_group, "time_bank_sendrecv", "description", &
         "Time between generations SEND/RECVing source sites (s)")
    call write_attribute_string(time_group, "time_tallies", "description", &
         "Time between batches accumulating tallies (s)")
    call write_attribute_string(time_group, "time_inactive", "description", &
         "Total time in inactive batches (s)")
    call write_attribute_string(time_group, "time_active", "description", &
         "Total time in active batches (s)")
    call write_attribute_string(time_group, "time_finalize", "description", &
         "Total time for finalization (s)")
    call write_attribute_string(time_group, "time_total", "description", &
         "Total time elapsed (s)")

    ! Write calculation rate
    total_particles = n_particles * n_batches * gen_per_batch
    speed = real(total_particles) / (time_inactive%elapsed + &
         time_active%elapsed)
    call write_dataset(time_group, "neutrons_per_second", speed)

    call close_group(time_group)
  end subroutine write_timing

end module summary
