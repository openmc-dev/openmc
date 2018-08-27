module summary

  use constants
  use endf,            only: reaction_name
  use error,           only: write_message
  use geometry_header
  use hdf5_interface
  use material_header, only: Material, n_materials, openmc_material_get_volume
  use mesh_header,     only: RegularMesh
  use message_passing
  use mgxs_interface
  use nuclide_header
  use output,          only: time_stamp
  use settings,        only: run_CE
  use surface_header
  use string,          only: to_str
  use tally_header,    only: TallyObject

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
    file_id = file_open("summary.h5", 'w')

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
    integer(HID_T) :: macro_group
    integer :: i
    character(12), allocatable :: nuc_names(:)
    character(12), allocatable :: macro_names(:)
    real(8), allocatable :: awrs(:)
    integer :: num_nuclides
    integer :: num_macros
    integer :: j
    integer :: k

    ! Find how many of these nuclides are macroscopic objects
    if (run_CE) then
      ! Then none are macroscopic
      num_nuclides = n_nuclides
      num_macros = 0
    else
      num_nuclides = 0
      num_macros = 0
      do i = 1, n_nuclides
        if (get_awr_c(i) /= MACROSCOPIC_AWR) then
          num_nuclides = num_nuclides + 1
        else
          num_macros = num_macros + 1
        end if
      end do
    end if

    ! Build array of nuclide names and awrs while only sorting nuclides from
    ! macroscopics
    if (num_nuclides > 0) then
      allocate(nuc_names(num_nuclides))
      allocate(awrs(num_nuclides))
    end if
    if (num_macros > 0) then
      allocate(macro_names(num_macros))
    end if

    j = 1
    k = 1
    do i = 1, n_nuclides
      if (run_CE) then
        nuc_names(i) = nuclides(i) % name
        awrs(i)     = nuclides(i) % awr
      else
        if (get_awr_c(i) /= MACROSCOPIC_AWR) then
          call get_name_c(i, len(nuc_names(j)), nuc_names(j))
          nuc_names(j) = trim(nuc_names(j))
          awrs(j)      = get_awr_c(i)
          j = j + 1
        else
          call get_name_c(i, len(macro_names(k)), macro_names(k))
          macro_names(k) = trim(macro_names(k))
          k = k + 1
        end if
      end if
    end do

    nuclide_group = create_group(file_id, "nuclides")
    call write_attribute(nuclide_group, "n_nuclides", num_nuclides)
    macro_group = create_group(file_id, "macroscopics")
    call write_attribute(macro_group, "n_macroscopics", num_macros)
    ! Write nuclide names and awrs
    if (num_nuclides > 0) then
      ! Write useful data from nuclide objects
      call write_dataset(nuclide_group, "names", nuc_names)
      call write_dataset(nuclide_group, "awrs", awrs)
    end if
    if (num_macros > 0) then
      ! Write useful data from macroscopic objects
      call write_dataset(macro_group, "names", macro_names)
    end if
    call close_group(nuclide_group)
    call close_group(macro_group)


    if (allocated(nuc_names)) deallocate(nuc_names, awrs)
    if (allocated(macro_names)) deallocate(macro_names)

  end subroutine write_nuclides

!===============================================================================
! WRITE_GEOMETRY
!===============================================================================

  subroutine write_geometry(file_id)
    integer(HID_T), intent(in) :: file_id

    integer :: i, j, cell_fill
    integer, allocatable :: cell_materials(:)
    integer, allocatable :: cell_ids(:)
    real(8), allocatable :: cell_temperatures(:)
    integer(HID_T) :: geom_group
    integer(HID_T) :: cells_group, cell_group
    integer(HID_T) :: surfaces_group
    integer(HID_T) :: universes_group, univ_group
    integer(HID_T) :: lattices_group
    type(Cell),     pointer :: c
    class(Lattice), pointer :: lat

    ! Use H5LT interface to write number of geometry objects
    geom_group = create_group(file_id, "geometry")
    call write_attribute(geom_group, "n_cells", n_cells)
    call write_attribute(geom_group, "n_surfaces", n_surfaces)
    call write_attribute(geom_group, "n_universes", n_universes)
    call write_attribute(geom_group, "n_lattices", size(lattices))

    ! ==========================================================================
    ! WRITE INFORMATION ON CELLS

    ! Create a cell group (nothing directly written in this group) then close
    cells_group = create_group(geom_group, "cells")

    ! Write information on each cell
    CELL_LOOP: do i = 1, n_cells
      c => cells(i)
      cell_group = create_group(cells_group, "cell " // trim(to_str(c%id())))

      call c % to_hdf5(cell_group)

      ! Write information on what fills this cell
      select case (c%type())
      case (FILL_MATERIAL)
        call write_dataset(cell_group, "fill_type", "material")

        if (c % material_size() == 1) then
          if (c % material(1) == MATERIAL_VOID) then
            call write_dataset(cell_group, "material", MATERIAL_VOID)
          else
            call write_dataset(cell_group, "material", &
                 materials(c % material(1)) % id())
          end if
        else
          allocate(cell_materials(c % material_size()))
          do j = 1, c % material_size()
            if (c % material(j) == MATERIAL_VOID) then
              cell_materials(j) = MATERIAL_VOID
            else
              cell_materials(j) = materials(c % material(j)) % id()
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
        call write_dataset(cell_group, "fill", universes(c%fill()+1)%id)

        if (allocated(c%translation)) then
          call write_dataset(cell_group, "translation", c%translation)
        end if
        if (allocated(c%rotation)) then
          call write_dataset(cell_group, "rotation", c%rotation)
        end if

      case (FILL_LATTICE)
        call write_dataset(cell_group, "fill_type", "lattice")
        ! Do not access the 'lattices' array with 'c % fill() + 1' directly; it
        ! causes a segfault in GCC 7.3.0
        cell_fill = c % fill() + 1
        call write_dataset(cell_group, "lattice", lattices(cell_fill)%obj%id())
      end select

      call close_group(cell_group)
    end do CELL_LOOP

    call close_group(cells_group)

    ! ==========================================================================
    ! WRITE INFORMATION ON SURFACES

    ! Create surfaces group
    surfaces_group = create_group(geom_group, "surfaces")

    ! Write information on each surface
    SURFACE_LOOP: do i = 1, n_surfaces
      call surfaces(i) % to_hdf5(surfaces_group)
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
            cell_ids(j) = cells(u % cells(j)) % id()
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

    lattices_group = create_group(geom_group, "lattices")

    do i = 1, size(lattices)
      lat => lattices(i)%obj
      call lat % to_hdf5(lattices_group)
    end do

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
    integer :: k
    integer :: n
    integer :: err
    character(20), allocatable :: nuc_names(:)
    character(20), allocatable :: macro_names(:)
    real(8) :: volume
    real(8), allocatable :: nuc_densities(:)
    integer :: num_nuclides
    integer :: num_macros
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
           trim(to_str(m%id())))

      if (m % depletable) then
        call write_attribute(material_group, "depletable", 1)
      else
        call write_attribute(material_group, "depletable", 0)
      end if

      err = openmc_material_get_volume(i, volume)
      if (err == 0 .and. volume > ZERO) then
        call write_attribute(material_group, "volume", volume)
      end if

      ! Write name for this material
      call write_dataset(material_group, "name", m % name)

      ! Write atom density with units
      call write_dataset(material_group, "atom_density", m % density)

      if (run_CE) then
        num_nuclides = m % n_nuclides
        num_macros = 0
      else
        ! Find the number of macroscopic and nuclide data in this material
        num_nuclides  = 0
        num_macros = 0
        do j = 1, m % n_nuclides
          if (get_awr_c(m % nuclide(j)) /= MACROSCOPIC_AWR) then
            num_nuclides = num_nuclides + 1
          else
            num_macros = num_macros + 1
          end if
        end do
      end if

      ! Copy ZAID or macro name for each nuclide to temporary array
      if (num_nuclides > 0) then
        allocate(nuc_names(num_nuclides))
        allocate(nuc_densities(num_nuclides))
      end if
      if (run_CE) then
        do j = 1, m % n_nuclides
          nuc_names(j) = nuclides(m%nuclide(j))%name
          nuc_densities(j) = m % atom_density(j)
        end do
      else
        if (num_macros > 0) then
          allocate(macro_names(num_macros))
        end if

        k = 1
        n = 1
        do j = 1, m % n_nuclides
          if (get_awr_c(m % nuclide(j)) /= MACROSCOPIC_AWR) then
            call get_name_c(m % nuclide(j), len(nuc_names(k)), nuc_names(k))
            nuc_names(k) = trim(nuc_names(k))
            nuc_densities(k) = m % atom_density(j)
            k = k + 1
          else
            call get_name_c(m % nuclide(j), len(macro_names(n)), macro_names(n))
            macro_names(n) = trim(macro_names(n))
            n = n + 1
          end if
        end do
      end if

      ! Write temporary array to 'nuclides'
      if (num_nuclides > 0) then
        call write_dataset(material_group, "nuclides", nuc_names)
        ! Deallocate temporary array
        deallocate(nuc_names)
        ! Write nuclide atom densities
        call write_dataset(material_group, "nuclide_densities", nuc_densities)
        deallocate(nuc_densities)
      end if

      ! Write temporary array to 'macroscopics'
      if (num_macros > 0) then
        call write_dataset(material_group, "macroscopics", macro_names)
        ! Deallocate temporary array
        deallocate(macro_names)
      end if

      if (m%n_sab > 0) then
        call write_dataset(material_group, "sab_names", m%sab_names)
      end if

      call close_group(material_group)
    end do

    call close_group(materials_group)

  end subroutine write_materials

end module summary
