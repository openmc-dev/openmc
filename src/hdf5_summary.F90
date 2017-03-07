module hdf5_summary

#ifdef HDF5

  use ace_header,      only: Reaction, UrrData, Nuclide
  use constants
  use endf,            only: reaction_name
  use geometry_header, only: Cell, Surface, Universe, Lattice
  use global
  use material_header, only: Material
  use mesh_header,     only: StructuredMesh
  use output_interface
  use output,          only: time_stamp
  use string,          only: to_str
  use tally_header,    only: TallyObject

  implicit none

  type(BinaryOutput) :: su

contains

!===============================================================================
! HDF5_WRITE_SUMMARY
!===============================================================================

  subroutine hdf5_write_summary()

    character(MAX_FILE_LEN) :: filename = "summary.h5"

    ! Create a new file using default properties.
    call su % file_create(filename)

    ! Write header information
    call hdf5_write_header()

    ! Write eigenvalue information
    if (run_mode == MODE_EIGENVALUE) then

      ! Write number of particles
      call su % write_data(n_particles, "n_particles")

      ! Use H5LT interface to write n_batches, n_inactive, and n_active
      call su % write_data(n_batches, "n_batches")
      call su % write_data(n_inactive, "n_inactive")
      call su % write_data(n_active, "n_active")
      call su % write_data(gen_per_batch, "gen_per_batch")

      ! Add description of each variable
      call su % write_attribute_string("n_particles", &
           "description", "Number of particles per generation")
      call su % write_attribute_string("n_batches", &
           "description", "Total number of batches")
      call su % write_attribute_string("n_inactive", &
           "description", "Number of inactive batches")
      call su % write_attribute_string("n_active", &
           "description", "Number of active batches")
      call su % write_attribute_string("gen_per_batch", &
           "description", "Number of generations per batch")
    end if

    call hdf5_write_geometry()
    call hdf5_write_materials()
    call hdf5_write_nuclides()
    if (n_tallies > 0) then
      call hdf5_write_tallies()
    end if

    ! Terminate access to the file.
    call su % file_close() 

  end subroutine hdf5_write_summary

!===============================================================================
! HDF5_WRITE_HEADER
!===============================================================================

  subroutine hdf5_write_header()

    ! Write version information
    call su % write_data(VERSION_MAJOR, "version_major")
    call su % write_data(VERSION_MINOR, "version_minor")
    call su % write_data(VERSION_RELEASE, "version_release") 

    ! Write current date and time
    call su % write_data(time_stamp(), "date_and_time")

    ! Write MPI information
    call su % write_data(n_procs, "n_procs")
    call su % write_attribute_string("n_procs", "description", &
         "Number of MPI processes") 

  end subroutine hdf5_write_header

!===============================================================================
! HDF5_WRITE_GEOMETRY
!===============================================================================

  subroutine hdf5_write_geometry()

    integer          :: i, j, k, m
    integer          :: n_x, n_y, n_z
    integer          :: length(3)
    integer, allocatable :: lattice_universes(:,:,:)
    type(Cell),     pointer :: c => null()
    type(Surface),  pointer :: s => null()
    type(Universe), pointer :: u => null()
    type(Lattice),  pointer :: lat => null()

    ! Use H5LT interface to write number of geometry objects
    call su % write_data(n_cells, "n_cells", group="geometry")
    call su % write_data(n_surfaces, "n_surfaces", group="geometry")
    call su % write_data(n_universes, "n_universes", group="geometry")
    call su % write_data(n_lattices, "n_lattices", group="geometry")

    ! ==========================================================================
    ! WRITE INFORMATION ON CELLS

    ! Create a cell group (nothing directly written in this group) then close
    call su % open_group("geometry/cells")
    call su % close_group()

    ! Write information on each cell
    CELL_LOOP: do i = 1, n_cells
      c => cells(i)

      ! Write universe for this cell
      call su % write_data(universes(c % universe) % id, "universe", &
           group="geometry/cells/cell " // trim(to_str(c % id)))

      ! Write information on what fills this cell
      select case (c % type)
      case (CELL_NORMAL)
        call su % write_data("normal", "fill_type", &
             group="geometry/cells/cell " // trim(to_str(c % id)))
        if (c % material == MATERIAL_VOID) then
          call su % write_data(-1, "material", &
               group="geometry/cells/cell " // trim(to_str(c % id)))
        else
          call su % write_data(materials(c % material) % id, "material", &
               group="geometry/cells/cell " // trim(to_str(c % id)))
        end if
      case (CELL_FILL)
        call su % write_data("universe", "fill_type", &
             group="geometry/cells/cell " // trim(to_str(c % id)))
        call su % write_data(universes(c % fill) % id, "material", &
             group="geometry/cells/cell " // trim(to_str(c % id))) 
      case (CELL_LATTICE)
        call su % write_data("lattice", "fill_type", &
             group="geometry/cells/cell " // trim(to_str(c % id)))
        call su % write_data(lattices(c % fill) % id, "lattice", &
             group="geometry/cells/cell " // trim(to_str(c % id))) 
      end select

      ! Write list of bounding surfaces
      if (c % n_surfaces > 0) then
        call su % write_data(c % surfaces, "surfaces", length= c % n_surfaces, &
             group="geometry/cells/cell " // trim(to_str(c % id)))
      end if

    end do CELL_LOOP

    ! ==========================================================================
    ! WRITE INFORMATION ON SURFACES

    ! Create surfaces group (nothing directly written here) then close
    call su % open_group("geometry/surfaces")
    call su % close_group()

    ! Write information on each surface
    SURFACE_LOOP: do i = 1, n_surfaces
      s => surfaces(i)

      ! Write surface type
      select case (s % type)
      case (SURF_PX)
        call su % write_data("X Plane", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_PY)
        call su % write_data("Y Plane", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_PZ)
        call su % write_data("Z Plane", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_PLANE)
        call su % write_data("Plane", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CYL_X)
        call su % write_data("X Cylinder", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CYL_Y)
        call su % write_data("Y Cylinder", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CYL_Z)
        call su % write_data("Z Cylinder", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_SPHERE)
        call su % write_data("Sphere", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CONE_X)
        call su % write_data("X Cone", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CONE_Y)
        call su % write_data("Y Cone", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (SURF_CONE_Z)
        call su % write_data("Z Cone", "type", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      end select

      ! Write coefficients for surface
      call su % write_data(s % coeffs, "coefficients", length=size(s % coeffs), &
           group="geometry/surfaces/surface " // trim(to_str(s % id)))

      ! Write positive neighbors
      if (allocated(s % neighbor_pos)) then
        call su % write_data(s % neighbor_pos, "neighbors_positive", &
             length=size(s % neighbor_pos), &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      end if

      ! Write negative neighbors
      if (allocated(s % neighbor_neg)) then
        call su % write_data(s % neighbor_neg, "neighbors_negative", &
             length=size(s % neighbor_neg), &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      end if

      ! Write boundary condition
      select case (s % bc)
      case (BC_TRANSMIT)
        call su % write_data("transmission", "boundary_condition", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (BC_VACUUM)
        call su % write_data("vacuum", "boundary_condition", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (BC_REFLECT)
        call su % write_data("reflective", "boundary_condition", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      case (BC_PERIODIC)
        call su % write_data("periodic", "boundary_condition", &
             group="geometry/surfaces/surface " // trim(to_str(s % id)))
      end select

    end do SURFACE_LOOP

    ! ==========================================================================
    ! WRITE INFORMATION ON UNIVERSES

    ! Create universes group (nothing directly written here) then close
    call su % open_group("geometry/universes")
    call su % close_group()

    ! Write information on each universe
    UNIVERSE_LOOP: do i = 1, n_universes
      u => universes(i)

      ! Write list of cells in this universe
      if (u % n_cells > 0) then
        call su % write_data(u % cells, "cells", length=u % n_cells, &
             group="geometry/universes/universe " // trim(to_str(u % id)))
      end if

    end do UNIVERSE_LOOP

    ! ==========================================================================
    ! WRITE INFORMATION ON LATTICES

    ! Create lattices group (nothing directly written here) then close
    call su % open_group("geometry/lattices")
    call su % close_group()

    ! Write information on each lattice
    LATTICE_LOOP: do i = 1, n_lattices
      lat => lattices(i)

      ! Write lattice type
      select case(lat % type)
      case (LATTICE_RECT)
        call su % write_data("rectangular", "type", &
             group="geometry/lattices/lattice " // trim(to_str(lat % id)))
      case (LATTICE_HEX)
        call su % write_data("hexagonal", "type", &
             group="geometry/lattices/lattice " // trim(to_str(lat % id)))
      end select

      ! Write lattice dimensions, lower left corner, and width of element
      call su % write_data(lat % dimension, "dimension", &
           length=lat % n_dimension, &
           group="geometry/lattices/lattice " // trim(to_str(lat % id)))
      call su % write_data(lat % lower_left, "lower_left", &
           length=lat % n_dimension, &
           group="geometry/lattices/lattice " // trim(to_str(lat % id)))
      call su % write_data(lat % width, "width", &
           length=lat % n_dimension, &
           group="geometry/lattices/lattice " // trim(to_str(lat % id)))

      ! Determine dimensions of lattice
      n_x = lat % dimension(1)
      n_y = lat % dimension(2)
      if (lat % n_dimension == 3) then
        n_z = lat % dimension(3)
      else
        n_z = 1
      end if
        
      ! Write lattice universes
      allocate(lattice_universes(n_x, n_y, n_z))
      do j = 1, n_x
        do k = 1, n_y
          do m = 1, n_z
            lattice_universes(j,k,m) = universes(lat % universes(j,k,m)) % id
          end do
        end do
      end do
      length = [n_x, n_y, n_z]
      call su % write_data(lattice_universes, "universes", length=length, &
           group="geometry/lattices/lattice " // trim(to_str(lat % id)))
      deallocate(lattice_universes)

    end do LATTICE_LOOP

  end subroutine hdf5_write_geometry

!===============================================================================
! HDF5_WRITE_MATERIALS
!===============================================================================

  subroutine hdf5_write_materials()

    integer          :: i
    integer          :: j
    integer, allocatable :: zaids(:)
    type(Material), pointer :: m => null()

    ! Use H5LT interface to write number of materials
    call su % write_data(n_materials, "n_materials", group="materials")

    ! Write information on each material
    do i = 1, n_materials
      m => materials(i)

      ! Write atom density with units
      call su % write_data(m % density, "atom_density", &
           group="materials/material " // trim(to_str(m % id)))
      call su % write_attribute_string("atom_density", "units", "atom/b-cm", &
           group="materials/material " // trim(to_str(m % id)))

      ! Copy ZAID for each nuclide to temporary array
      allocate(zaids(m % n_nuclides))
      do j = 1, m % n_nuclides
        zaids(j) = nuclides(m % nuclide(j)) % zaid
      end do

      ! Write temporary array to 'nuclides'
      call su % write_data(zaids, "nuclides", length=m % n_nuclides, &
           group="materials/material " // trim(to_str(m % id)))

      ! Deallocate temporary array
      deallocate(zaids)

      ! Write atom densities
      call su % write_data(m % atom_density, "nuclide_densities", &
           length=m % n_nuclides, &
           group="materials/material " // trim(to_str(m % id)))

      ! Write S(a,b) information if present
      if (m % n_sab > 0) then
        call su % write_data(m % i_sab_nuclides, "i_sab_nuclides", &
             length=m % n_sab, &
             group="materials/material " // trim(to_str(m % id)))
        call su % write_data(m % i_sab_tables, "i_sab_tables", &
             length=m % n_sab, &
             group="materials/material " // trim(to_str(m % id))) 
      end if

    end do

  end subroutine hdf5_write_materials

!===============================================================================
! HDF5_WRITE_TALLIES
!===============================================================================

  subroutine hdf5_write_tallies()

    integer           :: i, j
    integer, allocatable :: temp_array(:) ! nuclide bin array
    type(StructuredMesh), pointer :: m => null()
    type(TallyObject), pointer :: t => null()

    ! Write total number of meshes
    call su % write_data(n_meshes, "n_meshes", group="tallies")

    ! Write information for meshes
    MESH_LOOP: do i = 1, n_meshes
      m => meshes(i)

      ! Write type and number of dimensions
      call su % write_data(m % type, "type", &
           group="tallies/mesh " // trim(to_str(m % id)))

      call su % write_data(m % n_dimension, "n_dimension", &
           group="tallies/mesh " // trim(to_str(m % id)))

      ! Write mesh information
      call su % write_data(m % dimension, "dimension", &
           length=m % n_dimension, &
           group="tallies/mesh " // trim(to_str(m % id)))
      call su % write_data(m % lower_left, "lower_left", &
           length=m % n_dimension, &
           group="tallies/mesh " // trim(to_str(m % id)))
      call su % write_data(m % upper_right, "upper_right", &
           length=m % n_dimension, &
           group="tallies/mesh " // trim(to_str(m % id)))
      call su % write_data(m % width, "width", &
           length=m % n_dimension, &
           group="tallies/mesh " // trim(to_str(m % id)))

    end do MESH_LOOP

    ! Write number of tallies
    call su % write_data(n_tallies, "n_tallies", group="tallies")

    TALLY_METADATA: do i = 1, n_tallies
      ! Get pointer to tally
      t => tallies(i)

      ! Write size of each tally
      call su % write_data(t % total_score_bins, "total_score_bins", &
           group="tallies/tally " // trim(to_str(t % id)))
      call su % write_data(t % total_filter_bins, "total_filter_bins", &
           group="tallies/tally " // trim(to_str(t % id)))

      ! Write number of filters
      call su % write_data(t % n_filters, "n_filters", &
           group="tallies/tally " // trim(to_str(t % id)))

      FILTER_LOOP: do j = 1, t % n_filters
        ! Write type of filter
        call su % write_data(t % filters(j) % type, "type", &
             group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))

        ! Write number of bins for this filter
        call su % write_data(t % filters(j) % n_bins, "n_bins", &
             group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))

        ! Write filter bins
        if (t % filters(j) % type == FILTER_ENERGYIN .or. &
             t % filters(j) % type == FILTER_ENERGYOUT) then
          call su % write_data(t % filters(j) % real_bins, "bins", &
               length=size(t % filters(j) % real_bins), &
               group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))
        else
          call su % write_data(t % filters(j) % int_bins, "bins", &
               length=size(t % filters(j) % int_bins), &
               group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))
        end if

        ! Write name of type
        select case (t % filters(j) % type)
        case(FILTER_UNIVERSE)
          call su % write_data("universe", "type_name", &
               group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))
        case(FILTER_MATERIAL)
          call su % write_data("material", "type_name", &
               group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))
        case(FILTER_CELL)
          call su % write_data("cell", "type_name", &
               group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))
        case(FILTER_CELLBORN)
          call su % write_data("cellborn", "type_name", &
               group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))
        case(FILTER_SURFACE)
          call su % write_data("surface", "type_name", &
               group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))
        case(FILTER_MESH)
          call su % write_data("mesh", "type_name", &
               group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))
        case(FILTER_ENERGYIN)
          call su % write_data("energy", "type_name", &
          group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))
        case(FILTER_ENERGYOUT)
          call su % write_data("energyout", "type_name", &
          group="tallies/tally " // trim(to_str(t % id)) &
                // "/filter " // trim(to_str(j)))
        end select

      end do FILTER_LOOP

      ! Write number of nuclide bins
      call su % write_data(t % n_nuclide_bins, "n_nuclide_bins", &
           group="tallies/tally " // trim(to_str(t % id)))

      ! Create temporary array for nuclide bins
      allocate(temp_array(t % n_nuclide_bins))
      NUCLIDE_LOOP: do j = 1, t % n_nuclide_bins
        if (t % nuclide_bins(j) > 0) then
          temp_array(j) = nuclides(t % nuclide_bins(j)) % zaid
        else
          temp_array(j) = t % nuclide_bins(j)
        end if
      end do NUCLIDE_LOOP

      ! Write and deallocate nuclide bins
      call su % write_data(temp_array, "nuclide_bins", length=t % n_nuclide_bins, &
           group="tallies/tally " // trim(to_str(t % id)))
      deallocate(temp_array)

      ! Write number of score bins
      call su % write_data(t % n_score_bins, "n_score_bins", &
           group="tallies/tally " // trim(to_str(t % id)))
      call su % write_data(t % score_bins, "score_bins", length=t % n_score_bins, &
           group="tallies/tally " // trim(to_str(t % id)))

    end do TALLY_METADATA

  end subroutine hdf5_write_tallies

!===============================================================================
! HDF5_WRITE_NUCLIDES
!===============================================================================

  subroutine hdf5_write_nuclides()

    integer        :: i, j
    integer        :: size_total
    integer        :: size_xs
    integer        :: size_angle
    integer        :: size_energy
    type(Nuclide),  pointer :: nuc => null()
    type(Reaction), pointer :: rxn => null()
    type(UrrData),  pointer :: urr => null()

    ! Use H5LT interface to write number of nuclides
    call su % write_data(n_nuclides_total, "n_nuclides", group="nuclides")

    ! Write information on each nuclide
    NUCLIDE_LOOP: do i = 1, n_nuclides_total
      nuc => nuclides(i)

      ! Determine size of cross-sections
      size_xs = (5 + nuc % n_reaction) * nuc % n_grid * 8
      size_total = size_xs

      ! Write some basic attributes
      call su % write_data(nuc % zaid, "zaid", &
           group="nuclides/" // trim(nuc % name))
      call su % write_data(nuc % awr, "awr", &
           group="nuclides/" // trim(nuc % name))
      call su % write_data(nuc % kT, "kT", &
           group="nuclides/" // trim(nuc % name))
      call su % write_data(nuc % n_grid, "n_grid", &
           group="nuclides/" // trim(nuc % name))
      call su % write_data(nuc % n_reaction, "n_reactions", &
           group="nuclides/" // trim(nuc % name))
      call su % write_data(nuc % n_fission, "n_fission", &
           group="nuclides/" // trim(nuc % name))
      call su % write_data(size_xs, "size_xs", &
           group="nuclides/" // trim(nuc % name))

      ! =======================================================================
      ! WRITE INFORMATION ON EACH REACTION

      ! Create overall group for reactions and close it
      call su % open_group("nuclides/" // trim(nuc % name) // "/reactions")
      call su % close_group()

      RXN_LOOP: do j = 1, nuc % n_reaction
        ! Information on each reaction
        rxn => nuc % reactions(j)

        ! Determine size of angle distribution
        if (rxn % has_angle_dist) then
          size_angle = rxn % adist % n_energy * 16 + size(rxn % adist % data) * 8
        else
          size_angle = 0
        end if

        ! Determine size of energy distribution
        if (rxn % has_energy_dist) then
          size_energy = size(rxn % edist % data) * 8
        else
          size_energy = 0
        end if

        ! Write information on reaction
        call su % write_data(rxn % Q_value, "Q_value", &
             group="nuclides/" // trim(nuc % name) // "/reactions/" // &
             trim(reaction_name(rxn % MT)))
        call su % write_data(rxn % multiplicity, "multiplicity", &
             group="nuclides/" // trim(nuc % name) // "/reactions/" // &
             trim(reaction_name(rxn % MT)))
        call su % write_data(rxn % threshold, "threshold", &
             group="nuclides/" // trim(nuc % name) // "/reactions/" // &
             trim(reaction_name(rxn % MT)))
        call su % write_data(size_angle, "size_angle", &
             group="nuclides/" // trim(nuc % name) // "/reactions/" // &
             trim(reaction_name(rxn % MT)))
        call su % write_data(size_energy, "size_energy", &
             group="nuclides/" // trim(nuc % name) // "/reactions/" // &
             trim(reaction_name(rxn % MT)))

        ! Accumulate data size
        size_total = size_total + size_angle + size_energy
      end do RXN_LOOP

      ! =======================================================================
      ! WRITE INFORMATION ON URR PROBABILITY TABLES

      if (nuc % urr_present) then
        urr => nuc % urr_data
        call su % write_data(urr % n_energy, "urr_n_energy", &
             group="nuclides/" // trim(nuc % name))
        call su % write_data(urr % n_prob, "urr_n_prob", &
             group="nuclides/" // trim(nuc % name))
        call su % write_data(urr % interp, "urr_interp", &
             group="nuclides/" // trim(nuc % name))
        call su % write_data(urr % inelastic_flag, "urr_inelastic", &
             group="nuclides/" // trim(nuc % name))
        call su % write_data(urr % absorption_flag, "urr_absorption", &
             group="nuclides/" // trim(nuc % name))
        call su % write_data(urr % energy(1), "urr_min_E", &
             group="nuclides/" // trim(nuc % name))
        call su % write_data(urr % energy(urr % n_energy), "urr_max_E", &
             group="nuclides/" // trim(nuc % name))
      end if

      ! Write total memory used
      call su % write_data(size_total, "size_total", &
           group="nuclides/" // trim(nuc % name))

    end do NUCLIDE_LOOP

  end subroutine hdf5_write_nuclides

!===============================================================================
! HDF5_WRITE_TIMING
!===============================================================================

  subroutine hdf5_write_timing()

    integer(8)       :: total_particles
    real(8)          :: speed

    ! Write timing data
    call su % write_data(time_initialize % elapsed, "time_initialize", &
         group="timing")
    call su % write_data(time_read_xs % elapsed, "time_read_xs", &
         group="timing")
    call su % write_data(time_unionize % elapsed, "time_unionize", &
         group="timing")
    call su % write_data(time_transport % elapsed, "time_transport", &
         group="timing")
    call su % write_data(time_bank % elapsed, "time_bank", &
         group="timing")
    call su % write_data(time_bank_sample % elapsed, "time_bank_sample", &
         group="timing")
    call su % write_data(time_bank_sendrecv % elapsed, "time_bank_sendrecv", &
         group="timing")
    call su % write_data(time_tallies % elapsed, "time_tallies", &
         group="timing")
    call su % write_data(time_inactive % elapsed, "time_inactive", &
         group="timing")
    call su % write_data(time_active % elapsed, "time_active", &
         group="timing")
    call su % write_data(time_finalize % elapsed, "time_finalize", &
         group="timing")
    call su % write_data(time_total % elapsed, "time_total", &
         group="timing")

    ! Add descriptions to timing data
    call su % write_attribute_string("time_initialize", "description", &
         "Total time elapsed for initialization (s)", group="timing")
    call su % write_attribute_string("time_read_xs", "description", &
         "Time reading cross-section libraries (s)", group="timing")
    call su % write_attribute_string("time_unionize", "description", &
         "Time unionizing energy grid (s)", group="timing")
    call su % write_attribute_string("time_transport", "description", &
         "Time in transport only (s)", group="timing")
    call su % write_attribute_string("time_bank", "description", &
         "Total time synchronizing fission bank (s)", group="timing")
    call su % write_attribute_string("time_bank_sample", "description", &
         "Time between generations sampling source sites (s)", group="timing")
    call su % write_attribute_string("time_bank_sendrecv", "description", &
         "Time between generations SEND/RECVing source sites (s)", &
         group="timing")
    call su % write_attribute_string("time_tallies", "description", &
         "Time between batches accumulating tallies (s)", group="timing")
    call su % write_attribute_string("time_inactive", "description", &
         "Total time in inactive batches (s)", group="timing")
    call su % write_attribute_string("time_active", "description", &
         "Total time in active batches (s)", group="timing")
    call su % write_attribute_string("time_finalize", "description", &
         "Total time for finalization (s)", group="timing")
    call su % write_attribute_string("time_total", "description", &
         "Total time elapsed (s)", group="timing")

    ! Write calculation rate
    total_particles = n_particles * n_batches * gen_per_batch
    speed = real(total_particles) / (time_inactive % elapsed + &
         time_active % elapsed)
    call su % write_data(speed, "neutrons_per_second", group="timing")

  end subroutine hdf5_write_timing

#endif

end module hdf5_summary
