module mgxs_data

use constants
  use error,           only: fatal_error
  use global
  use macroxs_header
  use material_header, only: Material
  use nuclide_header
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
    integer :: i_listing    ! index in xs_listings array
    integer :: i_nuclide    ! index in nuclides
    character(12)  :: name  ! name of isotope, e.g. 92235.03c
    character(12)  :: alias ! alias of isotope, e.g. U-235.03c
    integer :: representation ! Data representation
    type(Material),    pointer :: mat
    type(SetChar) :: already_read
    type(Node), pointer :: doc => null()
    type(Node), pointer :: node_xsdata
    type(NodeList), pointer :: node_xsdata_list => null()
    logical :: file_exists
    integer :: error_code
    character(MAX_LINE_LEN) :: error_text, temp_str
    logical :: get_kfiss, get_fiss
    integer :: l

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
    n_listings = get_list_size(node_xsdata_list)

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
        if (tallies(i) % score_bins(l) == SCORE_FISSION) then
          get_fiss = .true.
        end if
      end do
      if (get_kfiss .and. get_fiss) exit
    end do

    ! ==========================================================================
    ! READ ALL ACE CROSS SECTION TABLES

    ! Loop over all files
    MATERIAL_LOOP: do i = 1, n_materials
      mat => materials(i)

      NUCLIDE_LOOP: do j = 1, mat % n_nuclides
        name = mat % names(j)

        if (.not. already_read % contains(name)) then
          i_listing = xs_listing_dict % get_key(to_lower(name))
          i_nuclide = mat % nuclide(j)
          name  = xs_listings(i_listing) % name
          alias = xs_listings(i_listing) % alias

          ! Get pointer to xsdata table XML node
          call get_list_item(node_xsdata_list, i_listing, node_xsdata)

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
            allocate(NuclideIso :: nuclides_MG(i_nuclide) % obj)
          case(MGXS_ANGLE)
            allocate(NuclideAngle :: nuclides_MG(i_nuclide) % obj)
          end select

          ! Now read in the data specific to the type we just declared
          call NuclideMG_init(nuclides_MG(i_nuclide) % obj, node_xsdata, &
                               energy_groups, get_kfiss, get_fiss, error_code, &
                               error_text)

          ! Keep track of what listing is associated with this nuclide
          nuclides_MG(i_nuclide) % obj % listing = i_listing

          ! Handle any errors
          if (error_code /= 0) then
            call fatal_error(trim(error_text))
          end if

          ! Add name and alias to dictionary
          call already_read % add(name)
          call already_read % add(alias)
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
! SAME_NUCLIDE_LIST creates a linked list for each nuclide containing the
! indices in the nuclides array of all other instances of that nuclide.  For
! example, the same nuclide may exist at multiple temperatures resulting
! in multiple entries in the nuclides array for a single zaid number.
!===============================================================================

  subroutine same_NuclideMG_list()

    integer :: i ! index in nuclides array
    integer :: j ! index in nuclides array

    do i = 1, n_nuclides_total
      do j = 1, n_nuclides_total
        if (nuclides_MG(i) % obj % zaid == nuclides_MG(j) % obj % zaid) then
          call nuclides_MG(i) % obj % nuc_list % push_back(j)
        end if
      end do
    end do

  end subroutine same_NuclideMG_list

!===============================================================================
! NUCLIDE_*_INIT reads in the data from the XML file, as already accessed
!===============================================================================

    subroutine NuclideMG_init(this, node_xsdata, groups, get_kfiss, get_fiss, &
                                  error_code, error_text)
      class(NuclideMG), intent(inout) :: this        ! Working Object
      type(Node), pointer, intent(in)  :: node_xsdata ! Data from data.xml
      integer, intent(in)              :: groups      ! Number of Energy groups
      logical, intent(in)              :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)              :: get_fiss ! Should we get fiss data?
      integer, intent(inout)           :: error_code  ! Code signifying error
      character(MAX_LINE_LEN), intent(inout) :: error_text ! Error message to print

      type(Node), pointer     :: node_legendre_mu
      character(MAX_LINE_LEN) :: temp_str
      logical                 :: enable_leg_mu

      ! Initialize error data
      error_code = 0
      error_text = ''

      ! Load the data
      call get_node_value(node_xsdata, "name", this % name)
      this % name = to_lower(this % name)
      if (check_for_node(node_xsdata, "kT")) then
        call get_node_value(node_xsdata, "kT", this % kT)
      else
        this % kT = ZERO
      end if
      if (check_for_node(node_xsdata, "zaid")) then
        call get_node_value(node_xsdata, "zaid", this % zaid)
      else
        this % zaid = -1
      end if
      if (check_for_node(node_xsdata, "scatt_type")) then
        call get_node_value(node_xsdata, "scatt_type", temp_str)
        temp_str = trim(to_lower(temp_str))
        if (temp_str == 'legendre') then
          this % scatt_type = ANGLE_LEGENDRE
        else if (temp_str == 'histogram') then
          this % scatt_type = ANGLE_HISTOGRAM
        else if (temp_str == 'tabular') then
          this % scatt_type = ANGLE_TABULAR
        else
          error_code = 1
          error_text = "Invalid Scatt Type Option!"
          return
        end if
      else
        this % scatt_type = ANGLE_LEGENDRE
      end if

      if (check_for_node(node_xsdata, "order")) then
        call get_node_value(node_xsdata, "order", this % order)
      else
        error_code = 1
        error_text = "Order Must Be Provided!"
        return
      end if

      ! Get scattering treatment
      if (check_for_node(node_xsdata, "tabular_legendre")) then
        call get_node_ptr(node_xsdata, "tabular_legendre", node_legendre_mu)
        if (check_for_node(node_legendre_mu, "enable")) then
          call get_node_value(node_legendre_mu, "enable", temp_str)
          temp_str = trim(to_lower(temp_str))
          if (temp_str == 'true' .or. temp_str == '1') then
            enable_leg_mu = .true.
          elseif (temp_str == 'false' .or. temp_str == '0') then
            enable_leg_mu = .false.
          else
            call fatal_error("Unrecognized tabular_legendre/enable: " // temp_str)
          end if
        else
          enable_leg_mu = .true.
          this % legendre_mu_points = 33
        end if
        if (enable_leg_mu .and. &
             check_for_node(node_legendre_mu, "num_points")) then
          call get_node_value(node_legendre_mu, "num_points", &
                              this % legendre_mu_points)
          if (this % legendre_mu_points <= 0) then
            call fatal_error("num_points element must be positive and non-zero!")
          end if
          this % legendre_mu_points = -1 * this % legendre_mu_points
        end if
      else
        this % legendre_mu_points = 1
      end if

      if (check_for_node(node_xsdata, "fissionable")) then
        call get_node_value(node_xsdata, "fissionable", temp_str)
        temp_str = to_lower(temp_str)
        if (trim(temp_str) == 'true' .or. trim(temp_str) == '1') then
          this % fissionable = .true.
        else
          this % fissionable = .false.
        end if
      else
        error_code = 1
        error_text = "Fissionable element must be set!"
        return
      end if

      select type(this)
      type is (NuclideIso)
        call NuclideIso_init(this, node_xsdata, groups, get_kfiss, get_fiss, &
                              error_code, error_text)
      type is (NuclideAngle)
        call NuclideAngle_init(this, node_xsdata, groups, get_kfiss, get_fiss, &
                                error_code, error_text)
      end select

    end subroutine NuclideMG_init

    subroutine NuclideIso_init(this, node_xsdata, groups, get_kfiss, get_fiss, &
                                error_code, error_text)
      class(NuclideIso), intent(inout)  :: this        ! Working Object
      type(Node), pointer, intent(in)    :: node_xsdata ! Data from data.xml
      integer, intent(in)                :: groups      ! Number of Energy groups
      logical, intent(in)                :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)                :: get_fiss    ! Need fiss data?
      integer, intent(inout)             :: error_code  ! Code signifying error
      character(MAX_LINE_LEN), intent(inout) :: error_text ! Error message to print

      real(8), allocatable :: temp_arr(:)
      integer :: arr_len
      integer :: order_dim

      ! Load the more specific data
      if (this % fissionable) then

        if (check_for_node(node_xsdata, "chi")) then
          ! Get chi
          allocate(this % chi(groups))
          call get_node_array(node_xsdata, "chi", this % chi)

          ! Get nu_fission (as a vector)
          if (check_for_node(node_xsdata, "nu_fission")) then
            allocate(temp_arr(groups * 1))
            call get_node_array(node_xsdata, "nu_fission", temp_arr)
            allocate(this % nu_fission(groups, 1))
            this % nu_fission = reshape(temp_arr, (/groups, 1/))
            deallocate(temp_arr)
          else
            error_code = 1
            error_text = "If fissionable, must provide nu_fission!"
            return
          end if

        else
          ! Get nu_fission (as a matrix)
          if (check_for_node(node_xsdata, "nu_fission")) then

            allocate(temp_arr(groups*groups))
            call get_node_array(node_xsdata, "nu_fission", temp_arr)
            allocate(this % nu_fission(groups, groups))
            this % nu_fission = reshape(temp_arr, (/groups, groups/))
            deallocate(temp_arr)
          else
            error_code = 1
            error_text = "If fissionable, must provide nu_fission!"
            return
          end if
        end if
        if (get_fiss) then
          allocate(this % fission(groups))
          if (check_for_node(node_xsdata, "fission")) then
            call get_node_array(node_xsdata, "fission", this % fission)
          else
            error_code = 1
            error_text = "Fission data missing, required due to fission&
                         & tallies in tallies.xml file!"
            return
          end if
        end if
        if (get_kfiss) then
          allocate(this % k_fission(groups))
          if (check_for_node(node_xsdata, "kappa_fission")) then
            call get_node_array(node_xsdata, "kappa_fission", this % k_fission)
          else
            error_code = 1
            error_text = "kappa_fission data missing, required due to &
                          &kappa-fission tallies in tallies.xml file!"
            return
          end if
        end if
      end if

      allocate(this % absorption(groups))
      if (check_for_node(node_xsdata, "absorption")) then
        call get_node_array(node_xsdata, "absorption", this % absorption)
      else
        error_code = 1
        error_text = "Must provide absorption!"
        return
      end if

      if (this % scatt_type == ANGLE_LEGENDRE) then
        order_dim = this % order + 1
      else if (this % scatt_type == ANGLE_HISTOGRAM) then
        order_dim = this % order
      else if (this % scatt_type == ANGLE_TABULAR) then
        order_dim = this % order
      end if

      allocate(this % scatter(groups, groups, order_dim))
      if (check_for_node(node_xsdata, "scatter")) then
        allocate(temp_arr(groups * groups * order_dim))
        call get_node_array(node_xsdata, "scatter", temp_arr)
        this % scatter = reshape(temp_arr, (/groups, groups, order_dim/))
        deallocate(temp_arr)
      else
        error_code = 1
        error_text = "Must provide scatter!"
        return
      end if


      allocate(this % total(groups))
      if (check_for_node(node_xsdata, "total")) then
        call get_node_array(node_xsdata, "total", this % total)
      else
        this % total = this % absorption + sum(this%scatter(:,:,1),dim=1)
      end if

      ! Get Mult Data
      allocate(this % mult(groups, groups))
      if (check_for_node(node_xsdata, "multiplicity")) then
        arr_len = get_arraysize_double(node_xsdata, "multiplicity")
        if (arr_len == groups * groups) then
          allocate(temp_arr(arr_len))
          call get_node_array(node_xsdata, "multiplicity", temp_arr)
          this % mult = reshape(temp_arr, (/groups, groups/))
          deallocate(temp_arr)
        else
          error_code = 1
          error_text = "Multiplicity Length Not Same as number of groups squared!"
          return
        end if
      else
        this % mult = ONE
      end if

    end subroutine NuclideIso_init

    subroutine NuclideAngle_init(this, node_xsdata, groups, get_kfiss, get_fiss, &
                                  error_code, error_text)
      class(NuclideAngle), intent(inout) :: this        ! Working Object
      type(Node), pointer, intent(in)     :: node_xsdata ! Data from data.xml
      integer, intent(in)                 :: groups      ! Number of Energy groups
      logical, intent(in)                 :: get_kfiss   ! Need Kappa-Fission?
      logical, intent(in)                 :: get_fiss    ! Should we get fiss data?
      integer, intent(inout)              :: error_code  ! Code signifying error
      character(MAX_LINE_LEN), intent(inout) :: error_text ! Error message to print

      real(8), allocatable :: temp_arr(:)
      integer :: arr_len
      real(8) :: dangle
      integer :: iangle
      integer :: order_dim

      if (this % scatt_type == ANGLE_LEGENDRE) then
        order_dim = this % order + 1
      else if (this % scatt_type == ANGLE_HISTOGRAM) then
        order_dim = this % order
      else if (this % scatt_type == ANGLE_TABULAR) then
        order_dim = this % order
      end if

      if (check_for_node(node_xsdata, "num_polar")) then
        call get_node_value(node_xsdata, "num_polar", this % Npol)
      else
        error_code = 1
        error_text = "num_polar Must Be Provided!"
        return
      end if

      if (check_for_node(node_xsdata, "num_azimuthal")) then
        call get_node_value(node_xsdata, "num_azimuthal", this % Nazi)
      else
        error_code = 1
        error_text = "num_azimuthal Must Be Provided!"
        return
      end if

      ! Load angle data, if present (else equally spaced)
      allocate(this % polar(this % Npol))
      allocate(this % azimuthal(this % Nazi))
      if (check_for_node(node_xsdata, "polar")) then
        error_code = 1
        error_text = "User-Specified polar angle bins not yet supported!"
        return
        call get_node_array(node_xsdata, "polar", this % polar)
      else
        dangle = PI / real(this % Npol,8)
        do iangle = 1, this % Npol
          this % polar(iangle) = (real(iangle,8) - 0.5_8) * dangle
        end do
      end if
      if (check_for_node(node_xsdata, "azimuthal")) then
        error_code = 1
        error_text = "User-Specified azimuthal angle bins not yet supported!"
        return
        call get_node_array(node_xsdata, "azimuthal", this % azimuthal)
      else
        dangle = TWO * PI / real(this % Nazi,8)
        do iangle = 1, this % Nazi
          this % azimuthal(iangle) = -PI + (real(iangle,8) - 0.5_8) * dangle
        end do
      end if

      ! Load the more specific data
      if (this % fissionable) then

        if (check_for_node(node_xsdata, "chi")) then
          ! Get chi
          allocate(temp_arr(groups * this % Nazi * this % Npol))
          call get_node_array(node_xsdata, "chi", temp_arr)
          allocate(this % chi(groups, this % Nazi, this % Npol))
          this % chi = reshape(temp_arr, (/groups, this % Nazi, this % Npol/))
          deallocate(temp_arr)

          ! Get nu_fission (as a vector)
          if (check_for_node(node_xsdata, "nu_fission")) then
            allocate(temp_arr(groups * this % Nazi * this % Npol))
            call get_node_array(node_xsdata, "nu_fission", temp_arr)
            allocate(this % nu_fission(groups, 1, this % Nazi, this % Npol))
            this % nu_fission = reshape(temp_arr, (/groups, 1, this % Nazi, &
                                                    this % Npol/))
            deallocate(temp_arr)
          else
            error_code = 1
            error_text = "If fissionable, must provide nu_fission!"
            return
          end if

        else
          ! Get nu_fission (as a matrix)
          if (check_for_node(node_xsdata, "nu_fission")) then

            allocate(temp_arr(groups * this % Nazi * this % Npol))
            call get_node_array(node_xsdata, "nu_fission", temp_arr)
            allocate(this % nu_fission(groups, groups, this % Nazi, this % Npol))
            this % nu_fission = reshape(temp_arr, (/groups, groups, &
                                                    this % Nazi, this % Npol/))
            deallocate(temp_arr)
          else
            error_code = 1
            error_text = "If fissionable, must provide nu_fission!"
            return
          end if
        end if
        if (get_fiss) then
          if (check_for_node(node_xsdata, "fission")) then
            allocate(temp_arr(groups * this % Nazi * this % Npol))
            call get_node_array(node_xsdata, "fission", temp_arr)
            allocate(this % fission(groups, this % Nazi, this % Npol))
            this % fission = reshape(temp_arr, (/groups, this % Nazi, this % Npol/))
            deallocate(temp_arr)
          else
            error_code = 1
            error_text = "Fission data missing, required due to fission&
                         & tallies in tallies.xml file!"
            return
          end if
        end if
        if (get_kfiss) then
          if (check_for_node(node_xsdata, "kappa_fission")) then
            allocate(temp_arr(groups * this % Nazi * this % Npol))
            call get_node_array(node_xsdata, "kappa_fission", temp_arr)
            allocate(this % k_fission(groups, this % Nazi, this % Npol))
            this % k_fission = reshape(temp_arr, (/groups, this % Nazi, this % Npol/))
            deallocate(temp_arr)
          else
            error_code = 1
            error_text = "kappa_fission data missing, required due to &
                          &kappa-fission tallies in tallies.xml file!"
            return
          end if
        end if
      end if

      if (check_for_node(node_xsdata, "absorption")) then
        allocate(temp_arr(groups * this % Nazi * this % Npol))
        call get_node_array(node_xsdata, "absorption", temp_arr)
        allocate(this % absorption(groups, this % Nazi, this % Npol))
        this % absorption = reshape(temp_arr, (/groups, this % Nazi, this % Npol/))
        deallocate(temp_arr)
      else
        error_code = 1
        error_text = "Must provide absorption!"
        return
      end if

      allocate(this % scatter(groups, groups, order_dim, this % Nazi, this % Npol))
      if (check_for_node(node_xsdata, "scatter")) then
        allocate(temp_arr(groups * groups * order_dim * this % Nazi * this%Npol))
        call get_node_array(node_xsdata, "scatter", temp_arr)
        this % scatter = reshape(temp_arr, (/groups, groups, order_dim, &
                                             this%Nazi,this%Npol/))
        deallocate(temp_arr)
      else
        error_code = 1
        error_text = "Must provide scatter!"
        return
      end if

      if (check_for_node(node_xsdata, "total")) then
        allocate(temp_arr(groups * this % Nazi * this % Npol))
        call get_node_array(node_xsdata, "total", temp_arr)
        allocate(this % total(groups, this % Nazi, this % Npol))
        this % total = reshape(temp_arr, (/groups, this % Nazi, this % Npol/))
        deallocate(temp_arr)
      else
        this % total = this % absorption + sum(this%scatter(:,:,1,:,:),dim=1)
      end if

      ! Get Mult Data
      allocate(this % mult(groups, groups, this % Nazi, this % Npol))
      if (check_for_node(node_xsdata, "multiplicity")) then
        arr_len = get_arraysize_double(node_xsdata, "multiplicity")
        if (arr_len == groups * groups * this % Nazi * this % Npol) then
          allocate(temp_arr(arr_len))
          call get_node_array(node_xsdata, "multiplicity", temp_arr)
          this % mult = reshape(temp_arr, (/groups, groups, this % Nazi, this % Npol/))
          deallocate(temp_arr)
        else
          error_code = 1
          error_text = "Multiplicity Length Does Not Match!"
          return
        end if
      else
        this % mult = ONE
      end if

    end subroutine NuclideAngle_init


!===============================================================================
! CREATE_MACRO_XS generates the macroscopic x/s from the microscopic input data
!===============================================================================

  subroutine create_macro_xs()
    integer :: i_mat ! index in materials array
    integer :: i             ! loop index over nuclides
    integer :: l             ! Loop over score bins
    type(Material), pointer :: mat ! current material
    logical :: get_kfiss, get_fiss
    integer :: error_code
    character(MAX_LINE_LEN) :: error_text
    integer :: representation
    integer :: scatt_type
    integer :: legendre_mu_points

    ! Find out if we need fission & kappa fission
    ! (i.e., are there any SCORE_FISSION or SCORE_KAPPA_FISSION tallies?)
    get_kfiss = .false.
    get_fiss  = .false.
    do i = 1, n_tallies
      do l = 1, tallies(i) % n_score_bins
        if (tallies(i) % score_bins(l) == SCORE_KAPPA_FISSION) then
          get_kfiss = .true.
        end if
        if (tallies(i) % score_bins(l) == SCORE_FISSION) then
          get_fiss = .true.
        end if
      end do
      if (get_kfiss .and. get_fiss) &
           exit
    end do

    allocate(macro_xs(n_materials))

    do i_mat = 1, n_materials
      mat => materials(i_mat)

      ! Check to see how our nuclides are represented
      ! Force all to be the same type
      ! Therefore type(nuclides(mat % nuclide(1)) % obj) dictates type(macroxs)
      ! At the same time, we will find the scattering type, as that will dictate
      ! how we allocate the scatter object within macroxs
      legendre_mu_points = nuclides_MG(mat % nuclide(1)) % obj % legendre_mu_points
      select type(nuc => nuclides_MG(mat % nuclide(1)) % obj)
      type is (NuclideIso)
        representation = MGXS_ISOTROPIC
      type is (NuclideAngle)
        representation = MGXS_ANGLE
      end select
      scatt_type = nuclides_MG(mat % nuclide(1)) % obj % scatt_type

      ! Now allocate accordingly
      select case(representation)
      case(MGXS_ISOTROPIC)
        allocate(MacroXSIso :: macro_xs(i_mat) % obj)
      case(MGXS_ANGLE)
        allocate(MacroXSAngle :: macro_xs(i_mat) % obj)
      end select

      call macro_xs(i_mat) % obj % init(mat, nuclides_MG, energy_groups, &
                                        get_kfiss, get_fiss, max_order, &
                                        scatt_type, legendre_mu_points, &
                                        error_code, error_text)
      ! Handle any errors
      if (error_code /= 0) then
        call fatal_error(trim(error_text))
      end if
    end do
  end subroutine create_macro_xs

end module mgxs_data