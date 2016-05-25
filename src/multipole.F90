module multipole

  use constants
  use global
  use hdf5
  use hdf5_interface
  use multipole_header, only: MultipoleArray, FIT_T, FIT_A, FIT_F, &
                              MP_FISS, FORM_MLBW, FORM_RM
  use search,           only: binary_search

  implicit none

contains

!===============================================================================
! MULTIPOLE_READ Reads in a multipole HDF5 file with the original API
! specification.  Subject to change as the library format matures.
!===============================================================================

  subroutine multipole_read(filename, multipole, i_table)
    character(len=*), intent(in)              :: filename  ! Filename of the
                                                           !  multipole library
                                                           !  to load
    type(MultipoleArray), intent(out), target :: multipole ! The object to fill
    integer, intent(in)                       :: i_table   ! index in nuclides/
                                                           !  sab_tables
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id
    integer :: is_fissionable
    integer :: cut1, cut2                   ! Old indices just outside MP region
    integer :: new_n_grid                   ! Number of points in new E grid
    real(8), allocatable :: new_energy(:)   ! New energy grid
    real(8), allocatable :: new_xs(:)       ! New cross sections
    integer :: i
    integer :: IE                           ! Reaction threshold

    associate (nuc => nuclides(i_table))

      !=========================================================================
      ! Copy in data from the file.

      ! Open file for reading and move into the /isotope group.
      file_id = file_open(filename, 'r', parallel=.true.)
      group_id = open_group(file_id, "/nuclide")

      ! Load in all the array size scalars.
      call read_dataset(multipole % length, group_id, "length")
      call read_dataset(multipole % windows, group_id, "windows")
      call read_dataset(multipole % num_l, group_id, "num_l")
      call read_dataset(multipole % fit_order, group_id, "fit_order")
      call read_dataset(multipole % max_w, group_id, "max_w")
      call read_dataset(is_fissionable, group_id, "fissionable")
      if (is_fissionable == MP_FISS) then
        multipole % fissionable = .true.
      else
        multipole % fissionable = .false.
      end if
      call read_dataset(multipole % formalism, group_id, "formalism")

      call read_dataset(multipole % spacing, group_id, "spacing")
      call read_dataset(multipole % sqrtAWR, group_id, "sqrtAWR")
      call read_dataset(multipole % start_E, group_id, "start_E")
      call read_dataset(multipole % end_E, group_id, "end_E")

      ! Allocate the multipole array components.
      call multipole % allocate()

      ! Read in arrays.
      call read_dataset(multipole % data, group_id, "data")
      call read_dataset(multipole % pseudo_k0RS, group_id, "pseudo_K0RS")
      call read_dataset(multipole % l_value, group_id, "l_value")
      call read_dataset(multipole % w_start, group_id, "w_start")
      call read_dataset(multipole % w_end, group_id, "w_end")
      call read_dataset(multipole % broaden_poly, group_id, "broaden_poly")

      call read_dataset(multipole % curvefit, group_id, "curvefit")

      ! Close the file.
      call close_group(group_id)
      call file_close(file_id)

      !=========================================================================
      ! Remove the uneeded/inconsitent pointwise data.  This step enforces the
      ! assumption that no inelastic scattering reactions can occur in the
      ! multipole region.  The energy grid is replaced with one that removes all
      ! energies covered by multiple and adds four new points.  Two new points
      ! mark the edges of the multipole region and cross sections will be
      ! interpolated to these points.  The other two points are used to zero the
      ! cross sections inside the multipole region.

      ! Find the points just outside the multipole region.
      cut1 = binary_search(nuc % energy, nuc % n_grid, insert_pts(1)) + 2
      cut2 = binary_search(nuc % energy, nuc % n_grid, insert_pts(4)) - 1

      ! Generate the new energy grid.
      new_n_grid = nuc % n_grid - (cut2 - cut1 - 1)
      allocate(new_energy(new_n_grid))
      new_energy(1:cut1) = nuc % energy(1:cut1)
      new_energy(cut1+1:new_n_grid) = nuc % energy(cut2:nuc % n_grid)

      ! Adjust the total cross section.
      allocate(new_xs(new_n_grid))
      new_xs(1:cut1-1) = nuc % total(1:cut1-1)
      new_xs(cut1:cut1+1) = ZERO
      new_xs(cut1+2:new_n_grid) = nuc % total(cut2+1:nuc % n_grid)
      call move_alloc(new_xs, nuc % total)

      ! Adjust the elastic cross section.
      allocate(new_xs(new_n_grid))
      new_xs(1:cut1-1) = nuc % elastic(1:cut1-1)
      new_xs(cut1:cut1+1) = ZERO
      new_xs(cut1+2:new_n_grid) = nuc % elastic(cut2+1:nuc % n_grid)
      call move_alloc(new_xs, nuc % elastic)

      ! Adjust the fission cross section.
      allocate(new_xs(new_n_grid))
      new_xs(1:cut1-1) = nuc % fission(1:cut1-1)
      new_xs(cut1:cut1+1) = ZERO
      new_xs(cut1+2:new_n_grid) = nuc % fission(cut2+1:nuc % n_grid)
      call move_alloc(new_xs, nuc % fission)

      ! Adjust the nu-fission cross section.
      allocate(new_xs(new_n_grid))
      new_xs(1:cut1-1) = nuc % nu_fission(1:cut1-1)
      new_xs(cut1:cut1+1) = ZERO
      new_xs(cut1+2:new_n_grid) = nuc % nu_fission(cut2+1:nuc % n_grid)
      call move_alloc(new_xs, nuc % nu_fission)

      ! Adjust the absorption cross section.
      allocate(new_xs(new_n_grid))
      new_xs(1:cut1-1) = nuc % absorption(1:cut1-1)
      new_xs(cut1:cut1+1) = ZERO
      new_xs(cut1+2:new_n_grid) = nuc % absorption(cut2+1:nuc % n_grid)
      call move_alloc(new_xs, nuc % absorption)

      ! Adjust other cross sections.
      do i = 1, nuc % n_reaction
        associate (rxn => nuc % reactions(i))
          if (.not. allocated(rxn % sigma)) cycle  ! Skip unallocated reactions
          IE = rxn % threshold
          if (rxn % threshold >= cut2) then
            ! The threshold is above the multipole range.  All we need to do
            ! is adjust the threshold index to match the new grid.
            rxn % threshold = rxn % threshold - (cut2 - cut1 - 1)
          else if (rxn % threshold <= cut1) then
            ! The threhold is below the multipole range.  Remove the multipole
            ! region just like we did with the other reactions.
            ! The new grid removed (cut2 - cut1 - 1) points.
            allocate(new_xs(size(rxn % sigma) - (cut2 - cut1 - 1)))
            new_xs(1:cut1-IE) = rxn % sigma(1:cut1-IE)
            new_xs(cut1-IE+1:cut1-IE+2) = ZERO
            new_xs(cut1-IE+3:size(new_xs)) = &
                 rxn % sigma(cut2-IE+2:size(rxn % sigma))
            call move_alloc(new_xs, rxn % sigma)
          else
            ! The threshold lies within the multipole range.  Remove the first
            ! cut2-IE+1 points.
            allocate(new_xs(size(rxn % sigma) - (cut2-IE+1)))
            new_xs(1:size(new_xs)) = rxn % sigma(cut2-IE+2:size(rxn % sigma))
            call move_alloc(new_xs, rxn % sigma)
            rxn % threshold = cut2 + 1
          end if
        end associate
      end do

      ! Apply the new energy grid.
      nuc % n_grid = new_n_grid
      call move_alloc(new_energy, nuc % energy)

    end associate

  end subroutine multipole_read

end module multipole
