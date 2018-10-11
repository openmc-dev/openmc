module plot

  use, intrinsic :: ISO_C_BINDING

  use error,           only: write_message
  use particle_header
  use plot_header
  use string,          only: to_str

  implicit none
  private

  public :: openmc_plot_geometry

  interface
    function openmc_plot_geometry_c() bind(C) result(err)
      import C_INT
      integer(C_INT) :: err
    end function openmc_plot_geometry_c
     
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

    err = openmc_plot_geometry_c()
    return

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
