module geometry

  use constants
  use error,                  only: fatal_error, warning, write_message
  use geometry_header
  use particle_header
  use simulation_header
  use settings
  use surface_header
  use stl_vector,             only: VectorInt
  use string,                 only: to_str

  use, intrinsic :: ISO_C_BINDING

  implicit none

  interface
    function count_universe_instances(search_univ, target_univ_id) bind(C) &
         result(count)
      import C_INT32_T, C_INT
      integer(C_INT32_T), intent(in), value :: search_univ
      integer(C_INT32_T), intent(in), value :: target_univ_id
      integer(C_INT)                        :: count
    end function count_universe_instances

    subroutine check_cell_overlap(p) bind(C)
      import Particle
      type(Particle), intent(in) :: p
    end subroutine check_cell_overlap

    function find_cell_c(p, search_surf) &
         bind(C, name="find_cell") result(found)
      import Particle, C_INT, C_BOOL
      type(Particle),  intent(inout)        :: p
      integer(C_INT),  intent(in), value    :: search_surf
      logical(C_BOOL)                       :: found
    end function find_cell_c

    subroutine cross_lattice(p, lattice_translation) &
         bind(C, name="cross_lattice")
      import Particle, C_INT
      type(Particle), intent(inout) :: p
      integer(C_INT), intent(in)    :: lattice_translation(3)
    end subroutine cross_lattice

    subroutine distance_to_boundary(p, dist, surface_crossed, &
         lattice_translation, next_level) bind(C)
      import Particle, C_DOUBLE, C_INT
      type(Particle), intent(inout) :: p
      real(C_DOUBLE), intent(out)   :: dist
      integer(C_INT), intent(out)   :: surface_crossed
      integer(C_INT), intent(out)   :: lattice_translation(3)
      integer(C_INT), intent(out)   :: next_level
    end subroutine distance_to_boundary

    subroutine neighbor_lists() bind(C)
    end subroutine neighbor_lists

#ifdef CAD
    
    function next_cell_c(current_cell, surface_crossed) &
      bind(C, name="next_cell") result(new_cell)
      import C_PTR, C_INT32_T
      type(C_PTR), intent(in), value :: current_cell
      type(C_PTR), intent(in), value :: surface_crossed
      integer(C_INT32_T)             :: new_cell
    end function next_cell_c
    
    function is_implicit_complement_C(cell) &
         bind(C, name="is_implicit_complement") result(res)
      import C_PTR, C_BOOL
      type(C_PTR), intent(in), value :: cell
      logical(C_BOOL)                :: res
    end function is_implicit_complement_C

#endif

 end interface
 
contains

  function cell_contains(c, p) result(in_cell)
    type(Cell), intent(in) :: c
    type(Particle), intent(in) :: p
    logical :: in_cell
    in_cell = cell_contains_c(c%ptr, p%coord(p%n_coord)%xyz, &
                              p%coord(p%n_coord)%uvw, p%surface)
  end function cell_contains

#ifdef CAD
  
  function next_cell(c, s) result(new_cell)
    type(Cell), intent(in) :: c
    type(Surface), intent(in) :: s
    integer :: new_cell
    new_cell = next_cell_c(c%ptr, s%ptr)
  end function next_cell

  function is_implicit_complement(c) result(res)
    type(Cell), intent(in) :: c
    logical:: res
    res = is_implicit_complement_c(c%ptr)
  end function is_implicit_complement

#endif
  
!===============================================================================
! FIND_CELL determines what cell a source particle is in within a particular
! universe. If the base universe is passed, the particle should be found as long
! as it's within the geometry
!===============================================================================

  subroutine find_cell(p, found, search_surf)
    type(Particle),    intent(inout) :: p
    logical,           intent(inout) :: found
    integer, optional, intent(in)    :: search_surf

    if (present(search_surf)) then
      found = find_cell_c(p, search_surf)
    else
      found = find_cell_c(p, 0)
    end if

  end subroutine find_cell

end module geometry
