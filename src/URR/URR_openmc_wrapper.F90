module URR_openmc_wrapper

  use ace_header,      only: Nuclide,&
       Reaction
  use error,           only: fatal_error,&
       warning
  use global,          only: master,&
       nuclides,&
       n_materials,&
       materials,&
       e_grid
  use list_header,     only: ListReal,&
                             ListInt
  use material_header, only: Material
  use output,          only: write_message
  use random_lcg,      only: prn
  use search,          only: binary_search
  use string,          only: to_str

  implicit none
  private
  public :: Nuclide,&
       fatal_error,&
       warning,&
       Reaction,&
       master,&
       nuclides,&
       n_materials,&
       materials,&
       e_grid,&
            ListReal,&
            ListInt,&
            Material,&
            write_message,&
            prn,&
            binary_search,&
            to_str

end module URR_openmc_wrapper
