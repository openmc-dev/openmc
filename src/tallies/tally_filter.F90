module tally_filter

  use, intrinsic :: ISO_C_BINDING

  use hdf5, only: HID_T

  use error
  use string,              only: to_f_string
  use tally_filter_header

  ! Inherit other filters
  use tally_filter_azimuthal
  use tally_filter_cell
  use tally_filter_cellborn
  use tally_filter_cellfrom
  use tally_filter_delayedgroup
  use tally_filter_distribcell
  use tally_filter_energy
  use tally_filter_energyfunc
  use tally_filter_legendre
  use tally_filter_material
  use tally_filter_mesh
  use tally_filter_meshsurface
  use tally_filter_mu
  use tally_filter_polar
  use tally_filter_sph_harm
  use tally_filter_sptl_legendre
  use tally_filter_surface
  use tally_filter_universe
  use tally_filter_zernike

  implicit none

contains

!===============================================================================
!                               C API FUNCTIONS
!===============================================================================

  function openmc_filter_get_type(index, type) result(err) bind(C)
    ! Get the type of a filter
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR), intent(out) :: type(*)
    integer(C_INT) :: err

    integer :: i
    character(20) :: type_

    err = verify_filter(index)
    if (err == 0) then
      ! Get type as a Fortran string
      select type (f => filters(index) % obj)
      type is (AzimuthalFilter)
        type_ = 'azimuthal'
      type is (CellFilter)
        type_ = 'cell'
      type is (CellbornFilter)
        type_ = 'cellborn'
      type is (CellfromFilter)
        type_ = 'cellfrom'
      type is (DelayedGroupFilter)
        type_ = 'delayedgroup'
      type is (DistribcellFilter)
        type_ = 'distribcell'
      type is (EnergyFilter)
        type_ = 'energy'
      type is (EnergyoutFilter)
        type_ = 'energyout'
      type is (EnergyFunctionFilter)
        type_ = 'energyfunction'
      type is (LegendreFilter)
        type_ = 'legendre'
      type is (MaterialFilter)
        type_ = 'material'
      type is (MeshFilter)
        type_ = 'mesh'
      type is (MeshSurfaceFilter)
        type_ = 'meshsurface'
      type is (MuFilter)
        type_ = 'mu'
      type is (PolarFilter)
        type_ = 'polar'
      type is (SphericalHarmonicsFilter)
        type_ = 'sphericalharmonics'
      type is (SpatialLegendreFilter)
        type_ = 'spatiallegendre'
      type is (SurfaceFilter)
        type_ = 'surface'
      type is (UniverseFilter)
        type_ = 'universe'
      type is (ZernikeFilter)
        type_ = 'zernike'
      end select

      ! Convert Fortran string to null-terminated C string. We assume the
      ! caller has allocated a char array buffer
      do i = 1, len_trim(type_)
        type(i) = type_(i:i)
      end do
      type(len_trim(type_) + 1) = C_NULL_CHAR
    end if

  end function openmc_filter_get_type


  function openmc_filter_set_type(index, type) result(err) bind(C)
    ! Set the type of a filter
    integer(C_INT32_T), value, intent(in) :: index
    character(kind=C_CHAR), intent(in) :: type(*)
    integer(C_INT) :: err

    character(:), allocatable :: type_

    ! Convert C string to Fortran string
    type_ = to_f_string(type)

    err = 0
    if (index >= 1 .and. index <= n_filters) then
      if (allocated(filters(index) % obj)) then
        err = E_ALLOCATE
        call set_errmsg("Filter type has already been set.")
      else
        select case (type_)
        case ('azimuthal')
          allocate(AzimuthalFilter :: filters(index) % obj)
        case ('cell')
          allocate(CellFilter :: filters(index) % obj)
        case ('cellborn')
          allocate(CellbornFilter :: filters(index) % obj)
        case ('cellfrom')
          allocate(CellfromFilter :: filters(index) % obj)
        case ('delayedgroup')
          allocate(DelayedGroupFilter :: filters(index) % obj)
        case ('distribcell')
          allocate(DistribcellFilter :: filters(index) % obj)
        case ('energy')
          allocate(EnergyFilter :: filters(index) % obj)
        case ('energyout')
          allocate(EnergyoutFilter :: filters(index) % obj)
        case ('energyfunction')
          allocate(EnergyFunctionFilter :: filters(index) % obj)
        case ('legendre')
          allocate(LegendreFilter :: filters(index) % obj)
        case ('material')
          allocate(MaterialFilter :: filters(index) % obj)
        case ('mesh')
          allocate(MeshFilter :: filters(index) % obj)
        case ('meshsurface')
          allocate(MeshSurfaceFilter :: filters(index) % obj)
        case ('mu')
          allocate(MuFilter :: filters(index) % obj)
        case ('polar')
          allocate(PolarFilter :: filters(index) % obj)
        case ('sphericalharmonics')
          allocate(SphericalHarmonicsFilter :: filters(index) % obj)
        case ('spatiallegendre')
          allocate(SpatialLegendreFilter :: filters(index) % obj)
        case ('surface')
          allocate(SurfaceFilter :: filters(index) % obj)
        case ('universe')
          allocate(UniverseFilter :: filters(index) % obj)
        case ('zernike')
          allocate(ZernikeFilter :: filters(index) % obj)
        case default
          err = E_UNASSIGNED
          call set_errmsg("Unknown filter type: " // trim(type_))
        end select
      end if
    else
      err = E_OUT_OF_BOUNDS
      call set_errmsg("Index in filters array is out of bounds.")
    end if
  end function openmc_filter_set_type

end module tally_filter
