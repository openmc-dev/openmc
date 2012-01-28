module cmfd_output

  implicit none
  private
  public :: write_cmfd_hdf5

contains

!===============================================================================
! WRITE_CMFD_HDF5 writes an hdf5 output file with the cmfd object for restarts
!===============================================================================

  subroutine write_cmfd_hdf5()

    use global, only: cmfd

#ifdef HDF5
    use hdf5
    use global, only: hdf5_output_file,hdf5_err
#endif

! character(LEN=7), parameter :: filename = "cmfd.h5" ! File name
    character(LEN=4), parameter :: grpname = "cmfd" ! Group name

! integer(HID_T) :: file_id ! File identifier
    integer(HID_T) :: group_id ! Group identifier
    integer(HID_T) :: dataspace_id ! Data space identifier
    integer(HID_T) :: dataset_id ! Dataset identifier

    integer(HSIZE_T), dimension(1) :: dim1 ! vector for hdf5 dimensions
    integer(HSIZE_T), dimension(3) :: dim3 ! vector for hdf5 dimensions
    integer(HSIZE_T), dimension(4) :: dim4 ! vector for hdf5 dimensions
    integer(HSIZE_T), dimension(5) :: dim5 ! vector for hdf5 dimensions

    integer :: nx ! number of mesh cells in x direction
    integer :: ny ! number of mesh cells in y direction
    integer :: nz ! number of mesh cells in z direction
    integer :: ng ! number of energy groups

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! create the CMFD group
    call h5gcreate_f(hdf5_output_file, grpname, group_id, hdf5_err)

    ! write indices from cmfd object
    dim1 = (/4/)
    call h5screate_simple_f(1,dim1,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/indices",H5T_NATIVE_INTEGER,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%indices,dim1,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write totalxs from cmfd object
    dim4 = (/ng,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/totalxs",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%totalxs,dim4,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write p1scattxs from cmfd object
    dim4 = (/ng,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/p1scattxs",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%p1scattxs,dim4,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write scattxs from cmfd object
    dim5 = (/ng,ng,nx,ny,nz/)
    call h5screate_simple_f(5,dim5,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/scattxs",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%scattxs,dim5,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write nfissxs from cmfd object
    dim5 = (/ng,ng,nx,ny,nz/)
    call h5screate_simple_f(5,dim5,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/nfissxs",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%nfissxs,dim5,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write diffcof from cmfd object
    dim4 = (/ng,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/diffcof",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%diffcof,dim4,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write current from cmfd object
    dim5 = (/12,ng,nx,ny,nz/)
    call h5screate_simple_f(5,dim5,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/current",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%current,dim5,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

   ! write flux from cmfd object
    dim4 = (/ng,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/flux",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%flux,dim4,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write dtilde from cmfd object
    dim5 = (/6,ng,nx,ny,nz/)
    call h5screate_simple_f(5,dim5,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/dtilde",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%dtilde,dim5,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write dhat from cmfd object
    dim5 = (/6,ng,nx,ny,nz/)
    call h5screate_simple_f(5,dim5,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/dhat",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%dhat,dim5,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write albedo from cmfd object
    dim1 = (/6/)
    call h5screate_simple_f(1,dim1,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/albedo",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%albedo,dim1,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write hxyz from cmfd object
    dim4 = (/3,nx,ny,nz/)
    call h5screate_simple_f(4,dim4,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/hxyz",H5T_NATIVE_DOUBLE,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,cmfd%hxyz,dim4,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write coremap from cmfd object
    dim3 = (/nx,ny,nz/)
    call h5screate_simple_f(3,dim3,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/coremap",H5T_NATIVE_INTEGER,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%coremap,dim3,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! write mat_dim from cmfd object
    dim1 = (/1/)
    call h5screate_simple_f(1,dim1,dataspace_id,hdf5_err)
    call h5dcreate_f(hdf5_output_file,"cmfd/mat_dim",H5T_NATIVE_INTEGER,dataspace_id, &
   & dataset_id,hdf5_err)
    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,cmfd%mat_dim,dim1,hdf5_err)
    call h5sclose_f(dataspace_id,hdf5_err)
    call h5dclose_f(dataset_id,hdf5_err)

    ! close the CMFD group
    call h5gclose_f(group_id, hdf5_err)

  end subroutine write_cmfd_hdf5

end module cmfd_output
