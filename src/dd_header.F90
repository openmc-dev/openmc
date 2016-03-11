module dd_header

  use dict_header,      only: DictIntInt
  use mesh_header,      only: RegularMesh
  use particle_header,  only: Particle, ParticleBuffer

  implicit none

  type, public :: DomainDecomType

    ! Domain mesh information
    type(RegularMesh), pointer :: mesh => null()
    integer :: n_domains
    integer :: meshbin               ! mesh domain of current processor
    integer :: ijk(3)                ! ijk corresponding to the meshbin
    logical :: allow_truncation      ! Whether or not we allow particles to
                                     ! leak out of the mesh (and stop tracking
                                     ! tracking them if they do)

    ! Info about the local 2nd-order neighborhood
    integer :: neighbor_meshbins(42) ! domain meshbins
    type(DictIntInt) :: bins_dict    ! reverse of neighbor_meshbins

    ! DD communication information
    integer :: comm   ! communicator for fission bank and tallies in this domain
    integer :: comm_domain_masters ! communicator for local domain masters
    integer :: rank   ! this processes's rank in domain sub-communicator
    logical :: local_master         ! master process in this domain?
    integer :: n_domain_procs       ! number of processors in domain sub-comm
    integer :: max_domain_procs     ! max number of procs in any other sub-comm

    ! User-specified ditribution of load
    real(8), allocatable :: domain_load_dist(:) ! index is domain meshbin

    ! Calculated number of procs per domain
    integer, allocatable :: domain_n_procs(:) ! index is domain meshbin

    ! Global ranks of the local masters in each domain, used for computing send
    ! ranks
    integer, allocatable :: domain_masters(:) ! index is domain meshbin

    ! For non-DD runs we re-ruse the same particle data structure for each
    ! history, but for DD runs we need to store particles if they would be
    ! transmitted to a neighboring domain.
    type(Particle), allocatable :: particle_buffer(:)
    integer(8) :: size_particle_buffer ! size of particle_buffer

    ! During simulation we keep track of where particles will scatter out to
    ! with p % outscatter_destination, and then send them later during
    ! syncronize_banks
    type(ParticleBuffer), allocatable :: send_buffer(:)
    type(ParticleBuffer), allocatable :: recv_buffer(:)
    integer(8) :: size_send_buffer ! size of buffer, resized if needed
    integer(8) :: size_recv_buffer ! size of buffer, resized if needed

    ! Scatter information
    integer              :: n_global_scatters ! # of scatters in last DD stage
    integer, allocatable :: send_rank_info(:) ! # of particles to send to bin
    integer, allocatable :: recv_rank_info(:) ! # of particles to recv from bin
    integer              :: n_inscatt         ! # of particles received
    integer, allocatable :: proc_finish(:)    ! finish idx in scatt array

    ! Before synchronizing scatters, all processors need to know how many
    ! particles are going everywhere
    integer, allocatable :: n_scatters_neighborhood(:,:) ! (from_bin, to_bin)
    integer, allocatable :: n_scatters_domain(:) ! (to_bin) just this domain
    integer, allocatable :: n_scatters_local(:) ! (to bin) just this processor
    integer, allocatable :: scatter_offest(:) ! (to bin) offset this processor

    ! To help with setting nodemaps, we can run in a mode that counts all
    ! particle interactions in a domain(a counter in the while(alive) particle
    ! loop.  This will be printed as output for the domain masters only
    logical    :: count_interactions = .false.
    integer(8) :: n_interaction = 0_8

  contains
    procedure :: deallocate => deallocate_dd

  end type DomainDecomType

contains

!===============================================================================
! DEALLOCATE_DD frees all memory of dd type
!===============================================================================

  subroutine deallocate_dd(this)

    class(DomainDecomType), intent(inout) :: this ! dd instance

    if (associated(this % mesh)) deallocate(this % mesh)
    if (allocated(this % domain_load_dist)) deallocate(this % domain_load_dist)
    if (allocated(this % domain_n_procs)) deallocate(this % domain_n_procs)
    if (allocated(this % domain_masters)) deallocate(this % domain_masters)
    if (allocated(this % particle_buffer)) deallocate(this % particle_buffer)
    if (allocated(this % send_buffer)) deallocate(this % send_buffer)
    if (allocated(this % recv_buffer)) deallocate(this % recv_buffer)
    if (allocated(this % send_rank_info)) deallocate(this % send_rank_info)
    if (allocated(this % recv_rank_info)) deallocate(this % recv_rank_info)
    if (allocated(this % proc_finish)) deallocate(this % proc_finish)
    if (allocated(this % n_scatters_neighborhood)) &
         deallocate(this % n_scatters_neighborhood)
    if (allocated(this % n_scatters_domain)) &
         deallocate(this % n_scatters_domain)
    if (allocated(this % n_scatters_local)) deallocate(this % n_scatters_local)
    if (allocated(this % scatter_offest)) deallocate(this % scatter_offest)
    call this % bins_dict % clear()

  end subroutine deallocate_dd

end module dd_header
