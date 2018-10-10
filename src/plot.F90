module plot

  use, intrinsic :: ISO_C_BINDING

  use constants
  use error,           only: fatal_error, write_message
  use geometry,        only: find_cell, check_cell_overlap
  use geometry_header, only: Cell, root_universe, cells
  use hdf5_interface
  use output,          only: time_stamp
  use material_header, only: materials
  use mesh_header,     only: meshes, RegularMesh
  use particle_header
  use plot_header
  use progress_header, only: ProgressBar
  use settings,        only: check_overlaps
  use string,          only: to_str

  implicit none
  private

  public :: openmc_plot_geometry

  integer, parameter :: RED = 1
  integer, parameter :: GREEN = 2
  integer, parameter :: BLUE = 3

  interface
    subroutine position_rgb(p, pl, rgb, id) bind(C)
      import Particle, ObjectPlot, C_INT
      type(Particle), intent(inout) :: p
      type(ObjectPlot), intent(in) :: pl
      integer(C_INT),   intent(out) :: rgb(3)
      integer(C_INT), intent(out) :: id
    end subroutine position_rgb


    subroutine create_ppm(pl) bind(C)
      import ObjectPlot
      type(ObjectPlot), intent(in) :: pl
    end subroutine create_ppm

    subroutine create_voxel(pl) bind(C)
      import ObjectPlot
      type(ObjectPlot), intent(in) :: pl
    end subroutine create_voxel

  end interface
contains

!===============================================================================
! RUN_PLOT controls the logic for making one or many plots
!===============================================================================

  function openmc_plot_geometry() result(err) bind(C)
    integer(C_INT) :: err

    integer :: i ! loop index for plots

    do i = 1, n_plots
      associate (pl => plots(i))
        ! Display output message
        call write_message("Processing plot " // trim(to_str(pl % id)) &
             // ": " // trim(pl % path_plot) // " ...", 5)

        if (pl % type == PLOT_TYPE_SLICE) then
          ! create 2d image
          call create_ppm(pl)
        else if (pl % type == PLOT_TYPE_VOXEL) then
          ! create dump for 3D silomesh utility script
          call create_voxel(pl)
        end if
      end associate
    end do

    err = 0
  end function openmc_plot_geometry

end module plot
