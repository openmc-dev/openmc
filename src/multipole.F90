module multipole

  use constants
  use global
  use hdf5
  use hdf5_interface
  use multipole_header, only: MultipoleArray, FIT_T, FIT_A, FIT_F, &
                              MP_FISS, FORM_MLBW, FORM_RM

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
    integer :: NMT
    integer :: i, j
    integer, allocatable :: MT(:)
    logical :: accumulated_fission
    character(len=24) :: MT_n ! Takes the form '/nuclide/reactions/MT???'
    integer :: is_fissionable

    associate (nuc => nuclides(i_table))

      ! Open file for reading and move into the /isotope group
      file_id = file_open(filename, 'r', parallel=.true.)
      group_id = open_group(file_id, "/nuclide")

      ! Load in all the array size scalars
      call read_dataset(group_id, "length", multipole % length)
      call read_dataset(group_id, "windows", multipole % windows)
      call read_dataset(group_id, "num_l", multipole % num_l)
      call read_dataset(group_id, "fit_order", multipole % fit_order)
      call read_dataset(group_id, "max_w", multipole % max_w)
      call read_dataset(group_id, "fissionable", is_fissionable)
      if (is_fissionable == MP_FISS) then
        multipole % fissionable = .true.
      else
        multipole % fissionable = .false.
      end if
      call read_dataset(group_id, "formalism", multipole % formalism)

      call read_dataset(group_id, "spacing", multipole % spacing)
      call read_dataset(group_id, "sqrtAWR", multipole % sqrtAWR)
      call read_dataset(group_id, "start_E", multipole % start_E)
      call read_dataset(group_id, "end_E", multipole % end_E)

      ! Allocate the multipole array components
      call multipole % allocate()

      ! Read in arrays
      call read_dataset(group_id, "data", multipole % data)
      call read_dataset(group_id, "pseudo_K0RS", multipole % pseudo_k0RS)
      call read_dataset(group_id, "l_value", multipole % l_value)
      call read_dataset(group_id, "w_start", multipole % w_start)
      call read_dataset(group_id, "w_end", multipole % w_end)
      call read_dataset(group_id, "broaden_poly", multipole % broaden_poly)

      call read_dataset(group_id, "curvefit", multipole % curvefit)

      ! Delete ACE pointwise data
      call read_dataset(group_id, "n_grid", nuc % n_grid)

      deallocate(nuc % energy)
      deallocate(nuc % total)
      deallocate(nuc % elastic)
      deallocate(nuc % fission)
      deallocate(nuc % nu_fission)
      deallocate(nuc % absorption)

      allocate(nuc % energy(nuc % n_grid))
      allocate(nuc % total(nuc % n_grid))
      allocate(nuc % elastic(nuc % n_grid))
      allocate(nuc % fission(nuc % n_grid))
      allocate(nuc % nu_fission(nuc % n_grid))
      allocate(nuc % absorption(nuc % n_grid))

      nuc % total(:) = ZERO
      nuc % absorption(:) = ZERO
      nuc % fission(:) = ZERO

      ! Read in new energy axis (converting eV to MeV)
      call read_dataset(group_id, "energy_points", nuc % energy)
      nuc % energy = nuc % energy / 1.0e6_8

      ! Get count and list of MT tables
      call read_dataset(group_id, "MT_count", NMT)
      allocate(MT(NMT))

      call read_dataset(group_id, "MT_list", MT)

      call close_group(group_id)

      accumulated_fission = .false.

      ! Loop over each MT entry and load it into a reaction.
      do i = 1, NMT
        write(MT_n, '(A, I3.3)') '/nuclide/reactions/MT', MT(i)

        group_id = open_group(file_id, MT_n)

        ! Each MT needs to be treated slightly differently.
        select case (MT(i))
          case(ELASTIC)
            call read_dataset(group_id, "MT_sigma", nuc % elastic)
            nuc % total(:) = nuc % total + nuc % elastic
          case(N_FISSION)
            call read_dataset(group_id, "MT_sigma", nuc % fission)
            nuc % total(:) = nuc % total + nuc % fission
            nuc % absorption(:) = nuc % absorption + nuc % fission
            accumulated_fission = .true.
          case default
            ! Search through all of our secondary reactions
            do j = 1, nuc % n_reaction
              if (nuc % reactions(j) % MT == MT(i)) then
                ! Match found

                ! Individual Fission components exist, so remove the combined
                ! fission cross section.
                if ( (MT(i) == N_F .or. MT(i) == N_NF .or. MT(i) == N_2NF &
                      .or. MT(i) == N_3NF) .and. accumulated_fission) then
                  nuc % total(:) = nuc % total - nuc % fission
                  nuc % absorption(:) = nuc % absorption - nuc % fission
                  nuc % fission(:) = ZERO
                  accumulated_fission = .false.
                end if

                deallocate(nuc % reactions(j) % sigma)
                allocate(nuc % reactions(j) % sigma(nuc % n_grid))

                call read_dataset(group_id, "MT_sigma", &
                                  nuc % reactions(j) % sigma)
                call read_dataset(group_id, "Q_value", &
                                  nuc % reactions(j) % Q_value)
                call read_dataset(group_id, "threshold", &
                                  nuc % reactions(j) % threshold)
                nuc % reactions(j) % threshold = 1 ! TODO: reconsider implications.
                nuc % reactions(j) % Q_value = nuc % reactions(j) % Q_value &
                                               / 1.0e6_8

                ! Accumulate total
                if (MT(i) /= N_LEVEL .and. MT(i) <= N_DA) then
                  nuc % total(:) = nuc % total + nuc % reactions(j) % sigma
                end if

                ! Accumulate absorption
                if (MT(i) >= N_GAMMA .and. MT(i) <= N_DA) then
                  nuc % absorption(:) = nuc % absorption &
                                        + nuc % reactions(j) % sigma
                end if

                ! Accumulate fission (if needed)
                if ( (MT(i) == N_F .or. MT(i) == N_NF .or. MT(i) == N_2NF &
                      .or. MT(i) == N_3NF) ) then
                  nuc % fission(:) = nuc % fission + nuc % reactions(j) % sigma
                  nuc % absorption(:) = nuc % absorption &
                                        + nuc % reactions(j) % sigma
                end if
              end if
            end do
        end select

        call close_group(group_id)
      end do

      ! Close file
      call file_close(file_id)

    end associate

  end subroutine multipole_read

end module multipole
