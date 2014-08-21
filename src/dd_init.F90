module dd_init

  use constants
  use dd_header,  only: dd_type
  use error,      only: fatal_error, warning
  use global,     only: domain_decomp, n_procs, rank, message, mpi_err
  use mesh,       only: bin_to_mesh_indices, mesh_indices_to_bin
  use output,     only: write_message
  use search,     only: binary_search
  use string,     only: to_str
  
  use mpi
  
  implicit none
  private
  public :: initialize_domain_decomp

contains

!===============================================================================
! INITIALIZE_DOMAIN_DECOMP
!===============================================================================

  subroutine initialize_domain_decomp()
    
    integer :: d, nd            ! neighbor and 2nd-neighbor indices
    integer :: neighbor_meshbin
    type(dd_type), pointer :: dd => domain_decomp

    message = "Initializing domain decomposition parameters..."
    call write_message(6)

    if (n_procs < dd % n_domains) then
      message = "Not enough processors for domain decomposition. Must have" // &
                " at least one per domain." // &
                " Domains: " // trim(to_str(dd % n_domains)) // &
                " Procs: " // trim(to_str(n_procs))
      call fatal_error()
    end if

    ! Determine number of processors per domain
    call calculate_domain_n_procs()

    ! Determine master ranks local to each domain
    allocate(dd % domain_masters(dd % n_domains))
    dd % domain_masters(1) = 0
    do d = 2, dd % n_domains
      dd % domain_masters(d) = dd % domain_masters(d - 1) + &
          dd % domain_n_procs(d - 1)
    end do 

    ! Set domain for this processor
    if (rank >= dd % domain_masters(dd % n_domains)) then
      dd % meshbin = dd % n_domains
    else
      dd % meshbin = binary_search(dd % domain_masters, dd % n_domains, rank)
    end if
    
    ! Find domain meshbin ijk
    call bin_to_mesh_indices(dd % mesh, dd % meshbin, dd % ijk)
    
    ! Determine meshbins of domains in the local neighborhood
    call set_neighbor_meshbins()
 
    ! Find the maximum number of processors working on any domain   
    dd % max_domain_procs = maxval(dd % domain_n_procs)

    ! Set mapping dict of mesh bin indices --> relative neighbor indices
    ! This is the reverse of dd_neighbor_meshbins
    do d = 1, N_CARTESIAN_NEIGHBORS
      neighbor_meshbin = dd % neighbor_meshbins(d)
      if (neighbor_meshbin == NO_BIN_FOUND) cycle
      call dd % bins_dict % add_key(neighbor_meshbin, d)
      do nd = 1, N_CARTESIAN_NEIGHBORS
        neighbor_meshbin = &
            dd % neighbor_meshbins(d * N_CARTESIAN_NEIGHBORS + nd)
        if (neighbor_meshbin == NO_BIN_FOUND) cycle
        if (neighbor_meshbin /= dd % meshbin) then
          call dd % bins_dict % add_key(neighbor_meshbin, &
              d * N_CARTESIAN_NEIGHBORS + nd)
        end if
      end do
    end do
    
    ! Initialize different MPI communicators for each domain, for fission bank
    ! and tally synchronization
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, dd % meshbin, rank, dd % comm, mpi_err)
    call MPI_COMM_SIZE(dd % comm, dd % n_domain_procs, mpi_err)
    call MPI_COMM_RANK(dd % comm, dd % rank, mpi_err)
    
    ! Set local_master
    if (dd % rank == 0) then
      dd % local_master = .true.
    else
      dd % local_master = .false.
    end if
    
    ! Allocate scattering information arrays and buffers
    allocate(dd % n_scatters_neighborhood(N_DD_NEIGHBORS,N_CARTESIAN_NEIGHBORS))
    allocate(dd % n_scatters_domain(N_CARTESIAN_NEIGHBORS)) 
    allocate(dd % n_scatters_local(N_CARTESIAN_NEIGHBORS))
    allocate(dd % scatter_offest(N_CARTESIAN_NEIGHBORS))
    allocate(dd % send_rank_info(dd % max_domain_procs * N_CARTESIAN_NEIGHBORS))
    allocate(dd % recv_rank_info(dd % max_domain_procs * N_CARTESIAN_NEIGHBORS))
    allocate(dd % proc_finish(dd % max_domain_procs))

    dd % n_scatters_local = 0
    dd % n_scatters_domain = 0
  
  end subroutine initialize_domain_decomp

!===============================================================================
! CALCULATE_DOMAIN_N_PROCS determines how many processes should work on each
! domain from the user-specified load distribution, or evenly distributes
! processors accross all domains if that wasn't specified
!===============================================================================

  subroutine calculate_domain_n_procs()

    integer                :: d
    integer                :: max_
    logical                :: forward
    real(8), allocatable   :: frac_nodes(:)
    type(dd_type), pointer :: dd => domain_decomp

    allocate(dd % domain_n_procs(dd % n_domains))

    allocate(frac_nodes(dd % n_domains))

    if (.not. allocated(dd % domain_load_dist)) then

      if (.not. mod(n_procs, dd % n_domains) == 0) then
        message = "Number of processes not evenly divisible by number of &
                  &domains. No <nodemap> was specified, so processes will be &
                  &distributed un-evenly"
        call warning()
      end if

      dd % domain_n_procs = n_procs / dd % n_domains
      do d = 1, dd % n_domains
        if (d - 1 < mod(n_procs, dd % n_domains)) &
          dd % domain_n_procs(d) = dd % domain_n_procs(d) + 1
      end do
      
    else

      ! Normalize load distribution
      dd % domain_load_dist = dd % domain_load_dist / sum(dd % domain_load_dist)

      ! Determine exact fractional number of nodes required
      frac_nodes = dd % domain_load_dist * real(n_procs, 8)
      
      ! Assign at least one processes per domain
      ! TODO: add support for allowing no processes to waste time on a domain
      ! with a user-input load of exactly zero.  For those domains, a fatal
      ! error would need to be thrown if a particle happened travel there. This
      ! feature would make sense if a domain is completely unreachable, outside
      ! of the problem boundary conditions (e.g., corner nodes with BEAVRS when
      ! using a StructuredMesh domain mesh
      dd % domain_n_procs = max(1, ceiling(frac_nodes))
      
      ! Perform a simple symmetric peak-shaving to make it match
      ! TODO: different matching strategies need to be explored, and support for
      ! having one processes handle multiple domains needs to be added.
      forward = .true.
      do while (sum(dd % domain_n_procs) > n_procs)
        max_ = maxval(dd % domain_n_procs)
        
        if (forward) then

          do d = 1, dd % n_domains
            if (dd % domain_n_procs(d) == max_) then
              dd % domain_n_procs(d) = dd % domain_n_procs(d) - 1
              exit
            end if
          end do
          
          forward = .false.

        else

          do d = dd % n_domains, 1, -1
            if (dd % domain_n_procs(d) == max_) then
              dd % domain_n_procs(d) = dd % domain_n_procs(d) - 1
              exit
            end if
          end do

          forward = .true.

        end if
        
      end do

    end if

    deallocate(frac_nodes)
    
  end subroutine calculate_domain_n_procs

!===============================================================================
! SET_NEIGHBOR_BINS determines the domian meshbins for the neighbors of the
! current domain with 1 and 2 degrees of separation
!===============================================================================

  subroutine set_neighbor_meshbins()
  
    type(dd_type), pointer :: dd => domain_decomp
    
    dd % neighbor_meshbins = (/ &
      ! First-order neighbors
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/-1,  0,  0/)), &  !-x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 1,  0,  0/)), &  !+x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0, -1,  0/)), &  !-y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0,  1,  0/)), &  !+y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0,  0, -1/)), &  !-z
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0,  0,  1/)), &  !+z
      ! Second-order neighbors (neighbors of the first-order neighbors)
      ! -x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/-2,  0,  0/)), & !-x
      dd % meshbin                                           , & !+x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/-1, -1,  0/)), & !-y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/-1,  1,  0/)), & !+y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/-1,  0, -1/)), & !-z
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/-1,  0,  1/)), & !+z
      ! +x
      dd % meshbin                                           , & !-x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 2,  0,  0/)), & !+x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 1, -1,  0/)), & !-y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 1,  1,  0/)), & !+y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 1,  0, -1/)), & !-z
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 1,  0,  1/)), & !+z
      ! -y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/-1, -1,  0/)), & !-x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 1, -1,  0/)), & !+x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0, -2,  0/)), & !-y
      dd % meshbin                                           , & !+y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0, -1, -1/)), & !-z
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0, -1,  1/)), & !+z
      ! +y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/-1,  1,  0/)), & !-x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 1,  1,  0/)), & !+x
      dd % meshbin                                           , & !-y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0,  2,  0/)), & !+y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0,  1, -1/)), & !-z
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0,  1,  1/)), & !+z
      ! -z
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/-1,  0, -1/)), & !-x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 1,  0, -1/)), & !+x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0, -1, -1/)), & !-y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0,  1, -1/)), & !+y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0,  0, -2/)), & !-z
      dd % meshbin                                           , & !+z
      ! +z
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/-1,  0,  1/)), & !-x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 1,  0,  1/)), & !+x
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0, -1,  1/)), & !-y
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0,  1,  1/)), & !+y
      dd % meshbin                                           , & !-z
      mesh_indices_to_bin(dd % mesh, dd % ijk + (/ 0,  0,  2/))  & !+z
    /)

  end subroutine set_neighbor_meshbins

end module dd_init
