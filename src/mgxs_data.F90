module mgxs_data

use constants
  use error,           only: fatal_error
  use global
  use material_header, only: Material
  use mgxs_header
  use output,          only: write_message
  use set_header,      only: SetChar
  use string,          only: to_lower
  use xml_interface

  implicit none

contains

!===============================================================================
! READ_XS reads all the cross sections for the problem and stores them in
! nuclides and sab_tables arrays
!===============================================================================

  subroutine read_mgxs()

    integer :: i            ! index in materials array
    integer :: j            ! index over nuclides in material
    integer :: i_xsdata     ! index in <xsdata> list
    integer :: i_nuclide    ! index in nuclides
    character(20)  :: name  ! name of isotope, e.g. 92235.03c
    integer :: representation ! Data representation
    type(Material),    pointer :: mat
    type(SetChar) :: already_read
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_xsdata
    type(NodeList), pointer :: node_xsdata_list => null()
    logical :: file_exists
    character(MAX_LINE_LEN) :: temp_str
    logical :: get_kfiss, get_fiss
    integer :: l
    type(DictCharInt) :: xsdata_dict

    ! Check if cross_sections.xml exists
    inquire(FILE=path_cross_sections, EXIST=file_exists)
    if (.not. file_exists) then
      ! Could not find cross_sections.xml file
      call fatal_error("Cross sections XML file '" &
           &// trim(path_cross_sections) // "' does not exist!")
    end if

    call write_message("Loading Cross Section Data...", 5)

    ! Parse cross_sections.xml file
    call open_xmldoc(doc, path_cross_sections)

    ! Get node list of all <xsdata>
    call get_node_list(doc, "xsdata", node_xsdata_list)

    ! Build dictionary mapping nuclide names to an index in the <xsdata> node
    ! list
    do i = 1, get_list_size(node_xsdata_list)
      ! Get pointer to xsdata table XML node
      call get_list_item(node_xsdata_list, i, node_xsdata)

      ! Get name and create pair (name, i)
      call get_node_value(node_xsdata, "name", name)
      call xsdata_dict % add_key(to_lower(name), i)
    end do

    ! allocate arrays for ACE table storage and cross section cache
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

          ! Get pointer to xsdata table XML node
          call get_list_item(node_xsdata_list, i_xsdata, node_xsdata)

          call write_message("Loading " // trim(name) // " Data...", 5)

          ! First find out the data representation
          if (check_for_node(node_xsdata, "representation")) then
            call get_node_value(node_xsdata, "representation", temp_str)
            temp_str = trim(to_lower(temp_str))
            if (temp_str == 'isotropic' .or. temp_str == 'iso') then
              representation = MGXS_ISOTROPIC
            else if (temp_str == 'angle') then
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
          call nuclides_MG(i_nuclide) % obj % init_file(node_xsdata, &
               energy_groups, get_kfiss, get_fiss, max_order)

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
    integer :: scatt_type

    allocate(macro_xs(n_materials))

    do i_mat = 1, n_materials
      mat => materials(i_mat)

      ! Check to see how our nuclides are represented
      ! Force all to be the same type
      ! Therefore type(nuclides(mat % nuclide(1)) % obj) dictates type(macroxs)
      ! At the same time, we will find the scattering type, as that will dictate
      ! how we allocate the scatter object within macroxs
      scatt_type = nuclides_MG(mat % nuclide(1)) % obj % scatt_type
      select type(nuc => nuclides_MG(mat % nuclide(1)) % obj)
      type is (MgxsIso)
        allocate(MgxsIso :: macro_xs(i_mat) % obj)
      type is (MgxsAngle)
        allocate(MgxsAngle :: macro_xs(i_mat) % obj)
      end select
      call macro_xs(i_mat) % obj % combine(mat, nuclides_MG, energy_groups, &
                                           max_order, scatt_type)
    end do
  end subroutine create_macro_xs

end module mgxs_data
