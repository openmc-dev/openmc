module summary

  use ace_header,      only: Reaction, UrrData, Nuclide
  use constants
  use endf,            only: reaction_name
  use geometry_header, only: Cell, Surface, Universe, Lattice, RectLattice, &
                             &HexLattice
  use global
  use hdf5_interface
  use material_header, only: Material
  use mesh_header,     only: StructuredMesh
  use output,          only: time_stamp
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
    call write_nuclides(file_id)
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
    integer(HID_T) :: geom_group
    integer(HID_T) :: cells_group, cell_group
    integer(HID_T) :: surfaces_group, surface_group
    integer(HID_T) :: universes_group, univ_group
    integer(HID_T) :: lattices_group, lattice_group
    type(Cell),     pointer :: c
    type(Surface),  pointer :: s
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
        if (c%material == MATERIAL_VOID) then
          call write_dataset(cell_group, "material", -1)
        else
          call write_dataset(cell_group, "material", materials(c%material)%id)
        end if

      case (CELL_FILL)
        call write_dataset(cell_group, "fill_type", "universe")
        call write_dataset(cell_group, "fill", universes(c%fill)%id)
        call write_dataset(cell_group, "maps", size(c%offset))
        if (size(c%offset) > 0) then
          call write_dataset(cell_group, "offset", c%offset)
        end if

        if (allocated(c%translation)) then
          call write_dataset(cell_group, "translated", 1)
          call write_dataset(cell_group, "translation", c%translation)
        else
          call write_dataset(cell_group, "translated", 0)
        end if

        if (allocated(c%rotation)) then
          call write_dataset(cell_group, "rotated", 1)
          call write_dataset(cell_group, "rotation", c%rotation)
        else
          call write_dataset(cell_group, "rotated", 0)
        end if

      case (CELL_LATTICE)
        call write_dataset(cell_group, "fill_type", "lattice")
        call write_dataset(cell_group, "lattice", lattices(c%fill)%obj%id)
      end select

      ! Write list of bounding surfaces
      if (c%n_surfaces > 0) then
        call write_dataset(cell_group, "surfaces", c%surfaces)
      end if

      call close_group(cell_group)
    end do CELL_LOOP

    call close_group(cells_group)

    ! ==========================================================================
    ! WRITE INFORMATION ON SURFACES

    ! Create surfaces group
    surfaces_group = create_group(geom_group, "surfaces")

    ! Write information on each surface
    SURFACE_LOOP: do i = 1, n_surfaces
      s => surfaces(i)
      surface_group = create_group(surfaces_group, "surface " // &
           trim(to_str(s%id)))

      ! Write internal OpenMC index for this surface
      call write_dataset(surface_group, "index", i)

      ! Write name for this surface
      call write_dataset(surface_group, "name", s%name)

      ! Write surface type
      select case (s%type)
      case (SURF_PX)
        call write_dataset(surface_group, "type", "X Plane")
      case (SURF_PY)
        call write_dataset(surface_group, "type", "Y Plane")
      case (SURF_PZ)
        call write_dataset(surface_group, "type", "Z Plane")
      case (SURF_PLANE)
        call write_dataset(surface_group, "type", "Plane")
      case (SURF_CYL_X)
        call write_dataset(surface_group, "type", "X Cylinder")
      case (SURF_CYL_Y)
        call write_dataset(surface_group, "type", "Y Cylinder")
      case (SURF_CYL_Z)
        call write_dataset(surface_group, "type", "Z Cylinder")
      case (SURF_SPHERE)
        call write_dataset(surface_group, "type", "Sphere")
      case (SURF_CONE_X)
        call write_dataset(surface_group, "type", "X Cone")
      case (SURF_CONE_Y)
        call write_dataset(surface_group, "type", "Y Cone")
      case (SURF_CONE_Z)
        call write_dataset(surface_group, "type", "Z Cone")
      end select

      ! Write coefficients for surface
      call write_dataset(surface_group, "coefficients", s%coeffs)

      ! Write positive neighbors
      if (allocated(s%neighbor_pos)) then
        call write_dataset(surface_group, "neighbors_positive", s%neighbor_pos)
      end if

      ! Write negative neighbors
      if (allocated(s%neighbor_neg)) then
        call write_dataset(surface_group, "neighbors_negative", s%neighbor_neg)
      end if

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

      ! Write name for this lattice
      call write_dataset(lattice_group, "name", lat%name)

      ! Write lattice type
      select type (lat)
      type is (RectLattice)
        ! Write lattice type.
        call write_dataset(lattice_group, "type", "rectangular")

        ! Write lattice dimensions, lower left corner, and pitch
        call write_dataset(lattice_group, "dimension", lat%n_cells)
        call write_dataset(lattice_group, "lower_left", lat%lower_left)
        call write_dataset(lattice_group, "pitch", lat%pitch)

        call write_dataset(lattice_group, "outer", lat%outer)
        call write_dataset(lattice_group, "offset_size", size(lat%offset))
        call write_dataset(lattice_group, "maps", size(lat%offset,1))

        if (size(lat%offset) > 0) then
          call write_dataset(lattice_group, "offsets", lat%offset)
        end if

        ! Write lattice universes.
        allocate(lattice_universes(lat%n_cells(1), lat%n_cells(2), &
             &lat%n_cells(3)))
        do j = 1, lat%n_cells(1)
          do k = 1, lat%n_cells(2)
            do m = 1, lat%n_cells(3)
              lattice_universes(j,k,m) = universes(lat%universes(j,k,m))%id
            end do
          end do
        end do
        call write_dataset(lattice_group, "universes", lattice_universes)
        deallocate(lattice_universes)

      type is (HexLattice)
        ! Write lattice type.
        call write_dataset(lattice_group, "type", "hexagonal")

        ! Write number of lattice cells.
        call write_dataset(lattice_group, "n_rings", lat%n_rings)
        call write_dataset(lattice_group, "n_axial", lat%n_axial)

        ! Write lattice center, pitch and outer universe.
        call write_dataset(lattice_group, "center", lat%center)
        call write_dataset(lattice_group, "pitch", lat%pitch)

        call write_dataset(lattice_group, "outer", lat%outer)
        call write_dataset(lattice_group, "offset_size", size(lat%offset))
        call write_dataset(lattice_group, "maps", size(lat%offset,1))

        if (size(lat%offset) > 0) then
          call write_dataset(lattice_group, "offsets", lat%offset)
        end if

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
        call write_dataset(lattice_group, "universes", lattice_universes)
        deallocate(lattice_universes)
      end select

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

    integer          :: i
    integer          :: j
    integer, allocatable :: zaids(:)
    integer(HID_T) :: materials_group
    integer(HID_T) :: material_group
    integer(HID_T) :: sab_group
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
      allocate(zaids(m%n_nuclides))
      do j = 1, m%n_nuclides
        zaids(j) = nuclides(m%nuclide(j))%zaid
      end do

      ! Write temporary array to 'nuclides'
      call write_dataset(material_group, "nuclides", zaids)

      ! Deallocate temporary array
      deallocate(zaids)

      ! Write atom densities
      call write_dataset(material_group, "nuclide_densities", m%atom_density)

      ! Write S(a,b) information if present
      call write_dataset(material_group, "n_sab", m%n_sab)

      if (m%n_sab > 0) then
        call write_dataset(material_group, "i_sab_nuclides", m%i_sab_nuclides)
        call write_dataset(material_group, "i_sab_tables", m%i_sab_tables)

        sab_group = create_group(material_group, "sab_tables")
        do j = 1, m%n_sab
          call write_dataset(sab_group, to_str(j), m%sab_names(j))
        end do
        call close_group(sab_group)
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

    integer           :: i, j
    integer, allocatable :: temp_array(:) ! nuclide bin array
    integer(HID_T) :: tallies_group
    integer(HID_T) :: mesh_group
    integer(HID_T) :: tally_group
    integer(HID_T) :: filter_group
    type(StructuredMesh), pointer :: m
    type(TallyObject), pointer :: t

    tallies_group = create_group(file_id, "tallies")

    ! Write total number of meshes
    call write_dataset(tallies_group, "n_meshes", n_meshes)

    ! Write information for meshes
    MESH_LOOP: do i = 1, n_meshes
      m => meshes(i)
      mesh_group = create_group(tallies_group, "mesh " // trim(to_str(m%id)))

      ! Write type and number of dimensions
      call write_dataset(mesh_group, "type", m%type)
      call write_dataset(mesh_group, "n_dimension", m%n_dimension)

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

      ! Write the name for this tally
      call write_dataset(tally_group, "name_size", len(t%name))
      if (len(t%name) > 0) then
        call write_dataset(tally_group, "name", t%name)
      endif

      ! Write size of each tally
      call write_dataset(tally_group, "total_score_bins", t%total_score_bins)
      call write_dataset(tally_group, "total_filter_bins", t%total_filter_bins)

      ! Write number of filters
      call write_dataset(tally_group, "n_filters", t%n_filters)

      FILTER_LOOP: do j = 1, t%n_filters
        filter_group = create_group(tally_group, "filter " // trim(to_str(j)))

        ! Write type of filter
        call write_dataset(filter_group, "type", t%filters(j)%type)

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
          call write_dataset(filter_group, "type_name", "universe")
        case(FILTER_MATERIAL)
          call write_dataset(filter_group, "type_name", "material")
        case(FILTER_CELL)
          call write_dataset(filter_group, "type_name", "cell")
        case(FILTER_CELLBORN)
          call write_dataset(filter_group, "type_name", "cellborn")
        case(FILTER_SURFACE)
          call write_dataset(filter_group, "type_name", "surface")
        case(FILTER_MESH)
          call write_dataset(filter_group, "type_name", "mesh")
        case(FILTER_ENERGYIN)
          call write_dataset(filter_group, "type_name", "energy")
        case(FILTER_ENERGYOUT)
          call write_dataset(filter_group, "type_name", "energyout")
        case(FILTER_MU)
          call write_dataset(filter_group, "type_name", "mu")
        case(FILTER_POLAR)
          call write_dataset(filter_group, "type_name", "polar")
        case(FILTER_AZIMUTHAL)
          call write_dataset(filter_group, "type_name", "azimuthal")
        end select

        call close_group(filter_group)
      end do FILTER_LOOP

      ! Write number of nuclide bins
      call write_dataset(tally_group, "n_nuclide_bins", t%n_nuclide_bins)

      ! Create temporary array for nuclide bins
      allocate(temp_array(t%n_nuclide_bins))
      NUCLIDE_LOOP: do j = 1, t%n_nuclide_bins
        if (t%nuclide_bins(j) > 0) then
          temp_array(j) = nuclides(t%nuclide_bins(j))%zaid
        else
          temp_array(j) = t%nuclide_bins(j)
        end if
      end do NUCLIDE_LOOP

      ! Write and deallocate nuclide bins
      call write_dataset(tally_group, "nuclide_bins", temp_array)
      deallocate(temp_array)

      ! Write number of score bins
      call write_dataset(tally_group, "n_score_bins", t%n_score_bins)
      call write_dataset(tally_group, "score_bins", t%score_bins)

      call close_group(tally_group)
    end do TALLY_METADATA

    call close_group(tallies_group)

  end subroutine write_tallies

!===============================================================================
! WRITE_NUCLIDES
!===============================================================================

  subroutine write_nuclides(file_id)
    integer(HID_T), intent(in) :: file_id

    integer        :: i, j
    integer        :: size_total
    integer        :: size_xs
    integer        :: size_angle
    integer        :: size_energy
    integer(HID_T) :: nuclides_group, nuclide_group
    integer(HID_T) :: reactions_group, rxn_group
    type(Nuclide),  pointer :: nuc
    type(Reaction), pointer :: rxn
    type(UrrData),  pointer :: urr

    nuclides_group = create_group(file_id, "nuclides")

    ! write number of nuclides
    call write_dataset(nuclides_group, "n_nuclides", n_nuclides_total)

    ! Write information on each nuclide
    NUCLIDE_LOOP: do i = 1, n_nuclides_total
      nuc => nuclides(i)
      nuclide_group = create_group(nuclides_group, nuc%name)

      ! Write internal OpenMC index for this nuclide
      call write_dataset(nuclide_group, "index", i)

      ! Determine size of cross-sections
      size_xs = (5 + nuc%n_reaction) * nuc%n_grid * 8
      size_total = size_xs

      ! Write some basic attributes
      call write_dataset(nuclide_group, "zaid", nuc%zaid)
      call write_dataset(nuclide_group, "alias", xs_listings(nuc%listing)%alias)
      call write_dataset(nuclide_group, "awr", nuc%awr)
      call write_dataset(nuclide_group, "kT", nuc%kT)
      call write_dataset(nuclide_group, "n_grid", nuc%n_grid)
      call write_dataset(nuclide_group, "n_reactions", nuc%n_reaction)
      call write_dataset(nuclide_group, "n_fission", nuc%n_fission)
      call write_dataset(nuclide_group, "size_xs", size_xs)

      ! =======================================================================
      ! WRITE INFORMATION ON EACH REACTION

      ! Create overall group for reactions and close it
      reactions_group = create_group(nuclide_group, "reactions")

      RXN_LOOP: do j = 1, nuc%n_reaction
        ! Information on each reaction
        rxn => nuc%reactions(j)
        rxn_group = create_group(reactions_group, trim(reaction_name(rxn%MT)))

        ! Determine size of angle distribution
        if (rxn%has_angle_dist) then
          size_angle = rxn%adist%n_energy * 16 + size(rxn%adist%data) * 8
        else
          size_angle = 0
        end if

        ! Determine size of energy distribution
        if (rxn%has_energy_dist) then
          size_energy = size(rxn%edist%data) * 8
        else
          size_energy = 0
        end if

        ! Write information on reaction
        call write_dataset(rxn_group, "Q_value", rxn%Q_value)
        call write_dataset(rxn_group, "multiplicity", rxn%multiplicity)
        call write_dataset(rxn_group, "threshold", rxn%threshold)
        call write_dataset(rxn_group, "size_angle", size_angle)
        call write_dataset(rxn_group, "size_energy", size_energy)

        ! Accumulate data size
        size_total = size_total + size_angle + size_energy

        call close_group(rxn_group)
      end do RXN_LOOP

      call close_group(reactions_group)

      ! =======================================================================
      ! WRITE INFORMATION ON URR PROBABILITY TABLES

      if (nuc%urr_present) then
        urr => nuc%urr_data
        call write_dataset(nuclide_group, "urr_n_energy", urr%n_energy)
        call write_dataset(nuclide_group, "urr_n_prob", urr%n_prob)
        call write_dataset(nuclide_group, "urr_interp", urr%interp)
        call write_dataset(nuclide_group, "urr_inelastic", urr%inelastic_flag)
        call write_dataset(nuclide_group, "urr_absorption", urr%absorption_flag)
        call write_dataset(nuclide_group, "urr_min_E", urr%energy(1))
        call write_dataset(nuclide_group, "urr_max_E", urr%energy(urr%n_energy))
      end if

      ! Write total memory used
      call write_dataset(nuclide_group, "size_total", size_total)

      call close_group(nuclide_group)
    end do NUCLIDE_LOOP

    call close_group(nuclides_group)

  end subroutine write_nuclides

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
