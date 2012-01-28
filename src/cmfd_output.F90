module cmfd_output

  implicit none
  private
  public :: write_cmfd_hdf5,write_cmfd_vtk

  integer :: nx         ! number of mesh cells in x direction
  integer :: ny         ! number of mesh cells in y direction
  integer :: nz         ! number of mesh cells in z direction
  integer :: ng         ! number of energy groups

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

!===============================================================================
! WRITE_CMFD_VTK outputs mesh data in vtk file for viewing
!===============================================================================

  subroutine write_cmfd_vtk()

    use global,      only: cmfd,meshes
    use mesh_header, only: StructuredMesh
    use vtk_writer

    integer :: i          ! x loop counter
    integer :: j          ! y loop counter
    integer :: k          ! z loop counter
    integer :: g          ! group counter
    integer :: n_idx      ! index in eigenvector
    real(8) :: x_m        ! -x coordinate
    real(8) :: x_p        ! +x coordinate
    real(8) :: y_m        ! -y coordinate
    real(8) :: y_p        ! +y coordinate
    real(8) :: z_m        ! -z coordinate
    real(8) :: z_p        ! +z coordinate
    real(8) :: x_uns(1:8) ! array of x points
    real(8) :: y_uns(1:8) ! array of y points
    real(8) :: z_uns(1:8) ! array of z points

    type(StructuredMesh), pointer :: m => null() ! pointer to mesh

    ! vtk specific variables
    integer :: E_IO              ! error code
    integer :: nn                ! number of nodes
    integer :: nc                ! number of cells
    integer :: con(8)            ! connectivity vector
    integer :: off(1:1)          ! offset, number of nodes in cell
    integer(1) :: cell_id(1:1)   ! cell type
    real(8) :: real_buffer(1:1)  ! real data buffer 8-byte
    character(len=40) :: varname ! name of output variable
    character(len=3) :: str_g    ! string for energy group #

    ! extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! point to mesh object
    m => meshes(1)

    ! set up vtk file
    E_IO = VTK_INI_XML(output_format = 'ASCII', &
   & filename = 'cmfd_unst.vtu', &
   & mesh_topology = 'UnstructuredGrid')

    ! set vtk parameters
    nn = 8
    nc = 1
    con = (/0,1,2,3,4,5,6,7/)
    off = (/8/)
    cell_id = (/11/)

    ! begin loop to construct mesh
    ZLOOP: do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          ! check for non accelerated region
          if (allocated(cmfd%coremap)) then
            if (cmfd%coremap(i,j,k) == 99999) then
              cycle
            end if
          end if

          ! calculate all coordinates
          x_m = dble(i - 1)*m%width(1) + m%origin(1)
          x_p = dble(i)*m%width(1) + m%origin(1)
          y_m = dble(j - 1)*m%width(2) + m%origin(2)
          y_p = dble(j)*m%width(2) + m%origin(2)

          z_m = dble(k - 1)*m%width(3) + m%origin(3)
          z_p = dble(k)*m%width(3) + m%origin(3)

          ! set up points arrays
          x_uns = (/x_m,x_p,x_m,x_p,x_m,x_p,x_m,x_p/)
          y_uns = (/y_m,y_m,y_p,y_p,y_m,y_m,y_p,y_p/)
          z_uns = (/z_m,z_m,z_m,z_m,z_p,z_p,z_p,z_p/)

          ! set up geometry piece
          E_IO = VTK_GEO_XML(nn,nc,x_uns,y_uns,z_uns)

          ! open data block in vtk file
          E_IO = VTK_DAT_XML('cell','open')

          ! loop around energy
          GROUP: do g = 1,ng

            ! convert group int to str
            write(str_g,'(I3)') g

            ! write out flux
            n_idx = get_matrix_idx(g,i,j,k)
            real_buffer = (/cmfd%phi(n_idx)/)
            varname = 'flux_'//trim(adjustl(str_g))
            E_IO = VTK_VAR_XML(nc,varname,real_buffer)

          end do GROUP

          ! close data block in vtk file
          E_IO = VTK_DAT_XML('cell','close')
          
          ! write out connectivity
          E_IO = VTK_CON_XML(nc,con,off,cell_id)

          ! close geometry piece
          E_IO = VTK_GEO_XML()

        end do XLOOP

      end do YLOOP

    end do ZLOOP

    ! close vtk file
    E_IO = VTK_END_XML()

  end subroutine write_cmfd_vtk

!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix 
!===============================================================================

  function get_matrix_idx(g,i,j,k)

    use global, only: cmfd

    ! arguments
    integer :: get_matrix_idx  ! the index location in matrix
    integer :: i               ! current x index
    integer :: j               ! current y index
    integer :: k               ! current z index
    integer :: g               ! current group index

    ! local variables
    integer :: nidx            ! index in matrix

    ! check if coremap is used
    if (allocated(cmfd % coremap)) then

      ! get idx from core map
      nidx = ng*(cmfd % coremap(i,j,k)) - (ng - g)

    else

      ! compute index
      nidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if

    ! record value to function
    get_matrix_idx = nidx

  end function get_matrix_idx

end module cmfd_output
