module mgxs_data

  use constants
  use algorithm,       only: find
  use error,           only: fatal_error
  use geometry_header, only: get_temperatures
  use global
  use hdf5_interface
  use material_header, only: Material
  use mgxs_header
  use output,          only: write_message
  use set_header,      only: SetChar
  use stl_vector,      only: VectorReal
  use string,          only: to_lower
  implicit none

contains

!===============================================================================
! READ_XS reads all the cross sections for the problem and stores them in
! nuclides and sab_tables arrays
!===============================================================================

  subroutine read_mgxs()
    integer :: i            ! index in materials array
    integer :: j            ! index over nuclides in material
    integer :: i_xsdata     ! index in xsdata_dict
    integer :: i_nuclide    ! index in nuclides array
    character(20)  :: name  ! name of library to load
    integer :: representation ! Data representation
    character(MAX_LINE_LEN) :: temp_str
    type(Material), pointer :: mat
    type(SetChar) :: already_read
    integer(HID_T) :: file_id
    integer(HID_T) :: xsdata_group
    logical :: file_exists
    logical :: get_kfiss, get_fiss
    integer :: l
    type(DictCharInt) :: xsdata_dict
    type(VectorReal), allocatable :: temps(:)

    ! Check if MGXS Library exists
    inquire(FILE=path_cross_sections, EXIST=file_exists)
    if (.not. file_exists) then
      ! Could not find MGXS Library file
      call fatal_error("Cross sections HDF5 file '" &
           &// trim(path_cross_sections) // "' does not exist!")
    end if

    call write_message("Loading Cross Section Data...", 5)

    ! Get temperatures
    call get_temperatures(cells, materials, material_dict, nuclide_dict, &
                          n_nuclides_total, temps)

    ! Open file for reading
    file_id = file_open(path_cross_sections, 'r', parallel=.true.)

    ! allocate arrays for MGXS storage and cross section cache
    allocate(nuclides_MG(n_nuclides_total))
!$omp parallel
    allocate(micro_xs(n_nuclides_total))
!$omp end parallel

    ! Find out if we need fission & kappa fission
    ! (i.e., are there any SCORE_FISSION or SCORE_KAPPA_FISSION tallies?)
    get_kfiss = .false.
    get_fiss  = .false.
    do i = 1, n_tallies
      do l = 1, tallies(i) % n_score_bins
        if (tallies(i) % score_bins(l) == SCORE_KAPPA_FISSION) then
          get_kfiss = .true.
        end if
        if (tallies(i) % score_bins(l) == SCORE_FISSION .or. &
             tallies(i) % score_bins(l) == SCORE_NU_FISSION) then
          get_fiss = .true.
        end if
      end do
      if (get_kfiss .and. get_fiss) exit
    end do

    ! ==========================================================================
    ! READ ALL MGXS CROSS SECTION TABLES

    ! Loop over all files
    MATERIAL_LOOP: do i = 1, n_materials
      mat => materials(i)

      NUCLIDE_LOOP: do j = 1, mat % n_nuclides
        name = mat % names(j)

        if (.not. already_read % contains(name)) then
          i_xsdata = xsdata_dict % get_key(to_lower(name))
          i_nuclide = mat % nuclide(j)

          call write_message("Loading " // trim(name) // " Data...", 5)

          ! Check to make sure cross section set exists in the library
          if (object_exists(file_id, trim(name))) then
            xsdata_group = open_group(file_id, trim(name))
          else
            call fatal_error("Data for '" // trim(name) // "' does not exist in "&
                 &// trim(path_cross_sections))
          end if

          ! First find out the data representation
          if (attribute_exists(xsdata_group, "representation")) then
            call read_attribute(temp_str, xsdata_group, "representation")
            if (trim(temp_str) == 'isotropic') then
              representation = MGXS_ISOTROPIC
            else if (trim(temp_str) == 'angle') then
              representation = MGXS_ANGLE
            else
              call fatal_error("Invalid Data Representation!")
            end if
          else
            ! Default to isotropic representation
            representation = MGXS_ISOTROPIC
          end if

          ! Now allocate accordingly
          select case(representation)
          case(MGXS_ISOTROPIC)
            allocate(MgxsIso :: nuclides_MG(i_nuclide) % obj)
          case(MGXS_ANGLE)
            allocate(MgxsAngle :: nuclides_MG(i_nuclide) % obj)
          end select

          ! Now read in the data specific to the type we just declared
          call nuclides_MG(i_nuclide) % obj % from_hdf5(xsdata_group, &
               energy_groups, temps(i_nuclide), temperature_method, &
               temperature_tolerance, get_kfiss, get_fiss, max_order, &
               legendre_to_tabular, legendre_to_tabular_points)

          ! Add name to dictionary
          call already_read % add(name)
        end if
      end do NUCLIDE_LOOP
    end do MATERIAL_LOOP

    ! Avoid some valgrind leak errors
    call already_read % clear()

    ! Loop around material
    MATERIAL_LOOP3: do i = 1, n_materials

      ! Get material
      mat => materials(i)

      ! Loop around nuclides in material
      NUCLIDE_LOOP2: do j = 1, mat % n_nuclides
        ! Is this fissionable?
        if (nuclides_MG(mat % nuclide(j)) % obj % fissionable) then
          mat % fissionable = .true.
        end if
        if (mat % fissionable) then
          exit NUCLIDE_LOOP2
        end if

      end do NUCLIDE_LOOP2
    end do MATERIAL_LOOP3

  end subroutine read_mgxs

!===============================================================================
! CREATE_MACRO_XS generates the macroscopic x/s from the microscopic input data
!===============================================================================

  subroutine create_macro_xs()
    integer :: i_mat ! index in materials array
    type(Material), pointer :: mat ! current material
    type(VectorReal), allocatable :: kTs(:)

    allocate(macro_xs(n_materials))

    ! Get temperatures to read for each material
    call get_mat_kTs(kTs)

    do i_mat = 1, n_materials
      mat => materials(i_mat)

      ! Check to see how our nuclides are represented
      ! Force all to be the same type
      ! Therefore type(nuclides(mat % nuclide(1)) % obj) dictates type(macroxs)
      select type(nuc => nuclides_MG(mat % nuclide(1)) % obj)
      type is (MgxsIso)
        allocate(MgxsIso :: macro_xs(i_mat) % obj)
      type is (MgxsAngle)
        allocate(MgxsAngle :: macro_xs(i_mat) % obj)
      end select
      ! Do not read materials which we do not actually use in the problem to
      ! save space
      if (allocated(kTs(i_mat) % data)) then
        call macro_xs(i_mat) % obj % combine(kTs(i_mat), mat, nuclides_MG, &
                                             energy_groups, max_order, &
                                             temperature_tolerance, &
                                             temperature_method)
      end if
    end do
  end subroutine create_macro_xs

!===============================================================================
! GET_MAT_kTs returns a list of temperatures (in eV) that each
! material appears at in the model.
!===============================================================================

  subroutine get_mat_kTs(kTs)
    type(VectorReal), allocatable, intent(out) :: kTs(:)

    integer :: i, j
    integer :: i_material  ! Index in materials array
    real(8) :: kT ! temperature in eV

    allocate(kTs(size(materials)))

    do i = 1, size(cells)
      do j = 1, size(cells(i) % material)
        ! Skip any non-material cells and void materials
        if (cells(i) % material(j) == NONE .or. &
             cells(i) % material(j) == MATERIAL_VOID) cycle

        ! Get temperature of cell (rounding to nearest integer)
        if (size(cells(i) % sqrtkT) > 1) then
          kT = cells(i) % sqrtkT(j)**2
        else
          kT = cells(i) % sqrtkT(1)**2
        end if

        i_material = material_dict % get_key(cells(i) % material(j))

        ! Add temperature if it hasn't already been added
        if (find(kTs(i_material), kT) == -1) then
          call kTs(i_material) % push_back(kT)
        end if

      end do
    end do

  end subroutine get_mat_kTs


end module mgxs_data
