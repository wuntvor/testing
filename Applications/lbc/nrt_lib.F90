module nrt_lib
   implicit none
#ifndef USE_ADCL
   include 'mpif.h'
#endif

#ifdef USE_MPI
#ifdef USE_ADCL
   include 'ADCL.inc'
#endif
#endif
    integer(4), parameter :: MAXDIMS=3

    integer, parameter :: I2B = selected_int_kind(2)
    integer, parameter :: I4B = selected_int_kind(4)
    integer, parameter :: I8B = selected_int_kind(8)
    integer, parameter :: R4B = selected_real_kind(4)
    integer, parameter :: R8B = selected_real_kind(8)

   ! set in params File which problem to solve
   integer, parameter                  :: cylinder=1,cavity=2,shock=3
   integer, parameter                  :: gaussian=4,cavity_same=5,flute=6
   integer, parameter                  :: corotating_vortex=7,cylinder_impulse=8
   integer, parameter                  :: caa_cavity=9,plate_trail=10,taylor_vortex=11
   integer, parameter                  :: gauss_convect=12,planar_standing_wave=13 
   integer, parameter                  :: fluid=0,wall=4,inlet=1,outlet=2
   ! special bcs including lid, non reflecting bc
   integer, parameter                  :: lid=3,sponge_wall=5,reset_wall=6,nr_wall=7
   integer, parameter                  :: per_wall=8,nr_wall_in=9,nr_wall_out=10 

   ! Model constants
   real(R8B), parameter :: cs  = 0.577350318431854_R8B  !1./sqrt(3.) !LU dx=dy=dt=1
   real(R8B), parameter :: cs2 = 0.333333333333333_R8B  
   real(R8B), parameter :: cs2inv = 3.0_R8B  
   real(R8B), parameter :: pi = 3.14159265358979
   ! simulation parameters
   integer,target       :: gtstep_cur 

   type sim_parameter
      real(R8B) :: rho0
      integer :: problem
      integer :: sponge_size ! thickness of sponge layer 
      ! description which models etc.
      character(len=14)    :: comm_method 
      character*16         :: problem_name
      character(len=4)     :: layoutname
      character(len=16)     :: modelname
      ! position of obstacle 
      integer         :: obst_x,obst_y,obst_z,obst_r
      ! max iteration, interval of output
      integer         :: tMax,tOut,tRamp,tOutBeg,tInitial

      ! problem parameters
      real(R8B)            :: Re,omega,umax,umaxRamp,nu
      real(R8B)            :: rhoges,lchar,obst_l
      ! number of Obstacles. Will be assigned in read_params
      integer         :: nObs
      logical         :: setObstacles,save_output      
      logical         :: goend,initial 
      integer         :: gx(3)
      ! global dimensions of domain
#ifdef F2003_ALLOC_DERIVED_TYPE
      real(R8B), dimension(:),allocatable      :: g_x,g_y,g_z
#else
      real(R8B), dimension(:),pointer          :: g_x => NULL() ,g_y => NULL() ,g_z => NULL() 
#endif
   end type sim_parameter


   type measure
      ! wall time measurement
      real(R8B)            :: tStart,tEnd,duration,mlups,total_duration
      ! communication measurement
      real(R8B)            :: comm_duration,tSt_comm,tEnd_comm,sbnd_duration,cceq_duration
      real(R8B)            :: coll_duration,prop_duration,ccmv_duration,bncb_duration,init_duration
      real(R8B)            :: mpl_cp_buf, mpl_cp_buf_back, mpl_exch, mpl_comm_size 
      real(R8B)            :: mpl_sync 
      ! mean values for communication measurement
      real(R8B)            :: mean_total,max_total
      real(R8B)            :: mean_comm,mean_sbnd,mean_cceq
      real(R8B)            :: mean_coll,mean_prop,mean_ccmv,mean_bncb
      real(R8B)            :: mean_mpl_cp_buf, mean_mpl_cp_buf_back, mean_mpl_exch, mean_mpl_comm_size 
      real(R8B)            :: mean_mpl_sync 
       ! memory consumption, global, per thread, and additional costs for master
      real(R8B)            :: mem,mem_thread,mem_master
   end type measure

   integer, parameter :: glob_nnghb = 18 !number of neighbors, is copied into prc%nnghb !CHANGED_REDUCE_MPI was 18
   integer, parameter :: rdcd_ngh = 6 !number of reduced neighbors
   integer, parameter :: n_ghost = 1     !number of ghost layers. is copied into prc%n_ghost_nodes

   type lb_block
      ! dimensions of current block
      integer                              :: lx(3),num_in(3)
      ! densities fIn, equilibrium densities, macr. speed
      ! macr. density
      integer :: nobs                            ! number of obstacles
!      integer :: nobs_sponge                     ! number of obstacles in sponge layer
      integer :: nnr_wall                        ! number of non-reflecting boundary nodes
!      integer :: nobs_reset                      ! number of obstacles in sponge layer
      integer :: nper_wall                       ! number of periodic boundary nodes  
#ifdef F2003_ALLOC_DERIVED_TYPE
!      real(R8B),dimension(:,:,:,:), allocatable,target :: fOut,u,u0
      real(R8B),dimension(:,:,:,:), allocatable :: u,u0
      real(R8B),dimension(:,:,:,:), pointer :: fOut => NULL()
      real(R8B),dimension(:,:,:),   allocatable :: rho
      real(R8B),dimension(:,:,:),   allocatable :: omega,gomega
      integer,dimension(:,:,:),allocatable      :: state
      integer,dimension(:,:),allocatable        :: obs  ! obstacle vector
      integer,dimension(:,:),allocatable        :: obs_sponge  ! obstacle vector
      integer,dimension(:,:),allocatable        :: obs_nrwall  ! obstacle vector
      integer,dimension(:,:),allocatable        :: obs_per  ! periodic vector
!      real(R8B),dimension(:,:),allocatable        :: per_val      ! 0 values for nrbc
      real(R8B),dimension(:,:),allocatable        :: nrwall_0val  ! 0 values for nrbc
      integer,dimension(:,:),allocatable        :: obs_reset   ! obstacle vector
      ! global variables, only needed when INIT_WITH_ROOT is set
      real(R8B),dimension(:,:,:),   allocatable :: grho
      real(R8B),dimension(:,:,:,:), allocatable :: gfIn, gu
      integer,dimension(:,:,:),allocatable      :: gstate
      ! coordinates
      real(R8B), dimension(:),allocatable      :: x,y,z
      ! refinement levels 
!#ifdef REFINE
!      integer,dimension(:,:,:),allocatable :: rfn_state
!      integer :: level
!      type(lb_block),pointer :: children(:) => NULL() 
!      integer,dimension(:,:,:),allocatable :: interpolate_children
!      type(lb_block),pointer :: parent => NULL() 
!      integer,dimension(:,:,:),allocatable :: interpolate_parent
!#endif /* REFINE */
!      real(R8B),allocatable,target  :: fIn(:,:,:,:)
      real(R8B),pointer  :: fIn(:,:,:,:) => NULL()
#else /* F2003_ALLOC_DERIVED_TYPE */
      real(R8B),dimension(:,:,:,:), pointer :: fOut => NULL() , u => NULL() , u0 => NULL() 
      real(R8B),dimension(:,:,:),   pointer :: rho => NULL() 
      real(R8B),dimension(:,:,:),   pointer :: omega => NULL() ,gomega => NULL() 
      integer,dimension(:,:,:),pointer      :: state => NULL() 
      integer,dimension(:,:),pointer        :: obs => NULL()   ! obstacle vector
!      integer,dimension(:,:),pointer        :: obs_sponge => NULL()   ! obstacle vector
      integer,dimension(:,:),pointer        :: obs_nrwall => NULL()   ! obstacle vector
      integer,dimension(:,:),pointer        :: obs_per => NULL()   ! periodic vector
!      real(R8B),dimension(:,:),pointer      :: per_val => NULL()       ! 0 values for nrbc
      real(R8B),dimension(:,:),pointer      :: nrwall_0val => NULL()   ! 0 values for nrbc
!      integer,dimension(:,:),pointer        :: obs_reset => NULL()   ! obstacle vector
      ! global variables, only needed when INIT_WITH_ROOT is set
      real(R8B),dimension(:,:,:),   pointer :: grho => NULL() 
      real(R8B),dimension(:,:,:,:), pointer :: gfIn => NULL() , gu => NULL() 
      integer,dimension(:,:,:),pointer      :: gstate => NULL() 
      ! coordinates
      real(R8B), dimension(:),pointer       :: x => NULL() ,y => NULL() ,z => NULL() 
      real(R8B),pointer      :: fIn(:,:,:,:) => NULL() 
#endif /* F2003_ALLOC_DERIVED_TYPE */


   end type lb_block

   integer :: send_dat(rdcd_ngh),recv_dat(rdcd_ngh) ! why do handles inside prc create errors??


   type mpl_var
#ifdef USE_CAF
      integer,allocatable :: caf_cart_comm[:,:,:]
      integer,allocatable :: caf_comm[:]
#endif
      character(len=14)               :: comm_method
      character(len=14)               :: dens_layout 
      integer                    :: rk, size, comm, ierr, cart_comm
      integer                    :: root_th, root_crd(3),crd(3)
      integer                    :: ndir_pdf(18)! => NULL()
!     integer                    :: send_buf(rdcd_ngh),recv_buf(rdcd_ngh)  !FIXME Handles inside prc creates very strange errors
      integer,pointer            :: drecv_pdf(:,:) => NULL() 
      integer,pointer            :: dsend_pdf(:,:) => NULL() 
      ! number of neighboring subdomains
      integer                    :: nnghb
      integer                    :: nnghb_avail       
      integer                    :: n_ghost_nodes       
      ! neighboring subdomains (compare model%c)
      !     1  x-       7  x+y+     13 x-z-  
      !     2  x+       8  x+y-     14 x-z+
      !     3  y-       9  x-y-     15 y+z+
      !     4  y+       10 x-y+     16 y+z-
      !     5  z-       11 x+z+     17 y-z-
      !     6  z+       12 x+z-     18 y-z+

      integer                  :: nghb(glob_nnghb) 
      integer                  :: team(rdcd_ngh) 
      integer                  :: ng_crd(3,glob_nnghb)
      ! Communication request handles
#if defined USE_MPI
      integer                  :: ireq(2*glob_nnghb) 
      integer :: status(mpi_status_size,2*glob_nnghb)
#endif

      ! boundaries of subdomain 
      ! indices: 1:3 crd(1:3)
      !          4   x = 1, y = 2, z = 3
      !          5   beginning = 1, end = 2
      integer,pointer          :: bnd(:,:,:,:,:)    => NULL() 
      
#ifdef USE_MPI
      integer                  :: stat(MPI_STATUS_SIZE) 
#ifdef USE_ADCL
      integer                  :: adcl_vmap, adcl_vec_fIn, adcl_vec_fOut, adcl_topo, & 
                                  adcl_request_fIn, adcl_request_fOut, adcl_active_req
!      real(R8B),allocatable    :: adcl_data(:,:,:,:)
      real(R8B),pointer        :: adcl_data(:,:,:,:) => NULL() 
!#!endif
#endif /* USE_ADCL */
#endif /* USE_MPI */
      ! number of mpi-processes in each cartesian direction
      integer                  :: np(3)    ! 
      integer                  :: pos_sends(glob_nnghb,3)
      integer                  :: pos_send (glob_nnghb,3)
      integer                  :: pos_recvs(glob_nnghb,3)
      integer                  :: pos_recv (glob_nnghb,3)
      integer                  :: length(glob_nnghb)
      integer                  :: length_recv(glob_nnghb)
#ifdef USE_CAF
#endif
   end type mpl_var
   integer, dimension(18,3) :: mpl_dir 

   type mpl_buf
      integer          :: idx_len
#ifdef F2003_ALLOC_DERIVED_TYPE
      integer,allocatable          :: idx(:)
      integer,allocatable          :: idx_recv(:)
      real(R8B),allocatable :: send(:,:,:,:) 
      real(R8B),allocatable :: recv(:,:,:,:) 
#else
      integer,pointer              :: idx(:) => NULL() 
      integer,pointer              :: idx_recv(:) => NULL() 
      real(R8B),pointer     :: send(:,:,:,:)  => NULL() 
      real(R8B),pointer     :: recv(:,:,:,:)  => NULL() 
#endif
   end type



#ifdef USE_CAF
#ifdef CAF_DER_TYPE
   type(mpl_buf) :: mpl_buffer(glob_nnghb)[*] 
#else
   real(R8B),allocatable,target,dimension(:,:,:,:)[:] :: send1,recv1,send2,recv2,send3,recv3
   real(R8B),allocatable,target,dimension(:,:,:,:)[:] :: send4,recv4,send5,recv5,send6,recv6
!   real(R8B),pointer,dimension(:,:,:,:) :: send_pnti,recv_pnti,send_pntn,recv_pntn

#endif
#else /* NOT USE_CAF -> USE_MPI*/
!#ifdef BUFFER_EACH /* BUFFER_EACH*/
   type(mpl_buf),save :: mpl_buffer(glob_nnghb) 
!#else /* BUFFER_EACH*/
!   type(mpl_buf) :: mpl_buffer(rdcd_ngh) 
!#endif /* BUFFER_EACH*/
#endif /* USE_CAF */  



#ifdef D2Q9
   integer(4), parameter                    :: NDIM= 2 
#endif
#ifdef D3Q19
   integer(4), parameter                    :: NDIM= 3 
#endif

   data mpl_dir(:,1) /-1,+1, 0, 0, 0, 0, +1, +1, -1, -1, +1, +1, -1, -1, 0, 0, 0, 0  /
   data mpl_dir(:,2) / 0, 0,-1,+1, 0, 0, +1, -1, -1, +1,  0,  0,  0,  0,+1,+1,-1,-1  /
   data mpl_dir(:,3) / 0, 0, 0, 0,-1,+1,  0,  0,  0,  0, +1, -1, -1, +1,+1,-1,-1,+1  / 

end module nrt_lib
!------------------------------------------------------------------------


