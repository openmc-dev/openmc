module multipole

  use hdf5

  use constants
  use error,            only: fatal_error
  use hdf5_interface
  use multipole_header, only: MultipoleArray, FIT_T, FIT_A, FIT_F, &
                              MP_FISS, FORM_MLBW, FORM_RM
  use nuclide_header,   only: nuclides

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

    ! Intermediate loading components
    integer :: is_fissionable
    character(len=10) :: version

    associate (nuc => nuclides(i_table))

      ! Open file for reading and move into the /isotope group
      file_id = file_open(filename, 'r', parallel=.true.)
      group_id = open_group(file_id, "/nuclide")

      ! Check the file version number.
      call read_dataset(version, file_id, "version")
      if (version /= VERSION_MULTIPOLE) call fatal_error("The current multipole&
           & format version is " // trim(VERSION_MULTIPOLE) // " but the file "&
           // trim(filename) // " uses version " // trim(version) // ".")

      ! Load in all the array size scalars
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

      ! Allocate the multipole array components
      call multipole % allocate()

      ! Read in arrays
      call read_dataset(multipole % data, group_id, "data")
      call read_dataset(multipole % pseudo_k0RS, group_id, "pseudo_K0RS")
      call read_dataset(multipole % l_value, group_id, "l_value")
      call read_dataset(multipole % w_start, group_id, "w_start")
      call read_dataset(multipole % w_end, group_id, "w_end")
      call read_dataset(multipole % broaden_poly, group_id, "broaden_poly")

      call read_dataset(multipole % curvefit, group_id, "curvefit")

      call close_group(group_id)

      ! Close file
      call file_close(file_id)

    end associate

  end subroutine multipole_read

end module multipole
