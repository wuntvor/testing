module mpl_lib 
!------------------------------------------------------------------------
! 
! Message Passing Library 
!
! Supports MPI and Coarray Fortran. 
! MPI is activated via the preprocessing Variable USE_MPI
!         Then a routine has to be chosen
! CAF ... USE_CAF
!
!
!
!
!
!
! Author: Manuel Hasert     
!         hasert@hlrs.de;  manuel.hasert@gmail.com
!         m.hasert@grs-sim.de
!         2009
!


   use nrt_lib
   use timing 

   implicit none
   integer :: target_rank

#ifdef USE_MPI
   integer,parameter :: mpl_proc_null = mpi_proc_null
   integer,parameter :: tag_rho=2
!   integer,parameter :: left_req=1,right_req=2,upper_req=3,lower_req=4,top_req=5,bottom_req=6
!   integer,parameter :: diagxy_req=7,diagyx_req=8,diagyz_req=9,diagzy_req=10,diagxz_req=11,diagzx_req=12
!   integer,parameter :: diagmxy_req=13,diagmyx_req=14,diagmyz_req=15,diagmzy_req=16,diagmxz_req=17,diagmzx_req=18
#else 
   integer,parameter :: mpl_proc_null = -2
#endif
#ifdef USE_CAF
   integer :: cobnd(3) 
#endif
end module mpl_lib


!------------------------------------------------------------------------







module mpl_set
   use mpl_lib
   use LBMODEL
   use tools
   implicit none
#include "include/replace.h"
contains


   subroutine comm_res(result_state,prc) 
   !------------------------------------------------------------------------
   !
   ! initialize parallel communication
   !
      type(mpl_var)   :: prc
      integer                       :: result_state,tmp
      integer,dimension(prc%size)   :: result_array

      if(prc%size > 1) then
#ifdef USE_MPI
         call mpi_allreduce(result_state,tmp,1,MPI_INTEGER,mpi_min,prc%cart_comm,prc%ierr)
         result_state=tmp
#else
      write(*,*) "result must be comm'ed"
      stop
#endif
      endif

   end subroutine comm_res


   subroutine mpl_init(meas,prc) 
   !------------------------------------------------------------------------
   !
   ! initialize parallel communication
   !
      type(measure),intent(inout) :: meas
      type(mpl_var)   :: prc


      prc%ierr           =  0
      prc%root_th        = -1
      prc%root_crd(:)    = -1
      ! copy values from global to local variable
#ifdef BUFFER_EACH
      prc%nnghb         = glob_nnghb
#else
      prc%nnghb         = rdcd_ngh 
#endif

      prc%n_ghost_nodes = n_ghost
      ! reset measure variable
      meas%comm_duration   = 0._R8B
      meas%duration        = 0._R8B
      meas%init_duration   = 0._R8B
      meas%sbnd_duration   = 0._R8B
      meas%cceq_duration   = 0._R8B
      meas%coll_duration   = 0._R8B
      meas%prop_duration   = 0._R8B
      meas%ccmv_duration   = 0._R8B
      meas%bncb_duration   = 0._R8B
      meas%mpl_cp_buf      = 0._R8B
      meas%mpl_cp_buf_back = 0._R8B
      meas%mpl_exch        = 0._R8B
      meas%mpl_sync        = 0._R8B
      meas%mpl_comm_size   = 0._R8B


#ifdef USE_MPI
      call MPI_INIT(prc%ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,prc%rk,prc%ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,prc%size,prc%ierr)



#ifdef USE_ADCL
      call ADCL_Init(prc%ierr)
      call ADCL_Vmap_halo_allocate(prc%n_ghost_nodes,prc%adcl_vmap,prc%ierr)
#endif
#endif  
#ifdef USE_CAF
      prc%size          = num_images()
      prc%rk            = this_image()-1
#endif

#ifndef USE_MPI 
#ifndef USE_CAF
      prc%rk            = 0
      prc%size          = 1
      prc%root_th       = 0
      prc%root_crd      = (/0,0,0/)
      prc%crd           = (/0,0,0/)
#endif
#endif

if(prc%rk==0) then
write(*,*) 
write(*,*) 
write(*,*) "             oooo   .o8                       "
write(*,*) "              888   888                       "
write(*,*) "              888   888oooo.   .ooooo.        "
write(*,*) "              888   d88   88b d88    Y8       "
write(*,*) "              888   888   888 888             "
write(*,*) "              888   888   888 888   .o8       "
write(*,*) "             o888o  `Y8bod8P  `Y8bod8P       "
write(*,*) ""
write(*,*) "              Lattice Boltzmann Kernel       "
write(*,*) " " 
endif
#if defined USE_MPI && defined USE_CAF
      if(prc%rk==0) then
         write(*,*) "Caution: CAF and MPI are both active !!!"
         write(*,*)
      endif
#endif


end subroutine mpl_init 





subroutine mpl_init_cart_comm(prc)
!------------------------------------------------------------------------
!
! Initialize the domain for parallel calculation.
! - a cartesian topology is created
! - calculate boundaries of the subdomain
! - identify neighbors of each process
!
! this is the reduced communication routine!
! only cartesian directions are exchanged after another
!
   type(mpl_var)                    :: prc
   integer                          :: allocstat
   logical                          :: periodic(3),reorder

   allocstat=0
   allocate(prc%bnd(prc%np(1),prc%np(2),prc%np(3),3,2),stat=allocstat)

   ! check if number of processes correspond 
   if(prc%np(1)*prc%np(2)*prc%np(3) /= prc%size) then 
      if(prc%rk == 0 ) then
         write(*,*) "Number of given Procs in lbc.param do not correspond with applied Procs of mpirun -np #"
         write(*,*) "Given in lbcparams", prc%np(1)*prc%np(2)*prc%np(3)   
         write(*,*) "Given as arguments", prc%size
      end if
#ifdef USE_MPI
      call MPI_FINALIZE(prc%ierr)
#endif
      stop
   end if

   periodic(1:3) = .false.
   reorder       = .true. 

#ifdef USE_MPI

   ! initialize the cartesian mpi-world and get current rank and coordinates
   call MPI_CART_CREATE(MPI_COMM_WORLD,3,prc%np(1:NDIM),&
                     &  periodic(1:NDIM),reorder,prc%cart_comm,prc%ierr)

   call MPI_COMM_RANK(prc%cart_comm,prc%rk,prc%ierr)

   call MPI_CART_COORDS(prc%cart_comm,prc%rk,3,prc%crd(1:NDIM),prc%ierr) 

#endif
#ifdef USE_CAF
   allocate(prc%caf_cart_comm[prc%np(1),prc%np(2),*])
   prc%crd = this_image(prc%caf_cart_comm)-1

#endif
      ! Set root thread number and coordinates
      prc%root_th     = 0 
      prc%root_crd    = (/0,0,0/)

end subroutine mpl_init_cart_comm






subroutine mpl_init_domain(lb_dom,s_par,prc)
!------------------------------------------------------------------------
!
! Initialize the domain for parallel calculation.
! - a cartesian topology is created
! - calculate boundaries of the subdomain
! - identify neighbors of each process
!
! this is the reduced communication routine!
! only cartesian directions are exchanged after another
!
   type(lb_block),intent(inout)     :: lb_dom
   type(sim_parameter),intent(in)   :: s_par
   type(mpl_var)                    :: prc

   call mpl_init_cart_comm(prc)

   call mpl_divide_and_set_domain_size(lb_dom,s_par,prc)

   call mpl_assign_neighbors(lb_dom,prc)

   call mpl_create_buffers(lb_dom,prc)

   call mpl_init_persistent(prc)

end subroutine mpl_init_domain 
!------------------------------------------------------------------------





subroutine mpl_init_persistent(prc)
!------------------------------------------------------------------------
!
!
   type(mpl_var)                    :: prc
   integer                          :: neighbor,counter

#ifdef PERSISTENT2
      counter = 0

   do neighbor = 1,prc%nnghb
      if (prc%nghb(neighbor) /= mpl_proc_null .and. prc%nghb(neighbor) .ge. 0) then 
      counter = counter+1
#ifdef MPI_SUBARRAY
      call mpi_recv_init(fIn,1,recv_dat(counter),prc%nghb(counter),mpi_any_tag,      & 
                     prc%cart_comm,prc%ireq(counter),prc%ierr)

      call mpi_send_init(fIn,1,send_dat(counter),prc%nghb(counter),1,                &
                     prc%cart_comm,prc%ireq(prc%nnghb_avail+counter),prc%ierr)

#else /* MPI_SUBARRAY */
      call mpi_send_init(mpl_buffer(neighbor)%send,                &
                     prc%length(neighbor),                         & 
                     mpi_double_precision,                         &
                     prc%nghb(neighbor),  1,                       &
                     prc%cart_comm,prc%ireq(prc%nnghb_avail+counter),prc%ierr)

      call mpi_recv_init(mpl_buffer(neighbor)%recv,                &
                     prc%length(neighbor),                         & 
                     mpi_double_precision,                         &
                     prc%nghb(neighbor),mpi_any_tag,               &
                     prc%cart_comm,prc%ireq(counter),prc%ierr)
#endif /* MPI_SUBARRAY */
      endif
   enddo
#endif /* PERSISTENT */

end subroutine mpl_init_persistent






subroutine mpl_divide_and_set_domain_size(lb_dom,s_par,prc)
!------------------------------------------------------------------------
!
!
   type(lb_block)                   :: lb_dom
   type(sim_parameter),intent(in)   :: s_par
   type(mpl_var)                    :: prc
   integer                          :: i,j,k
      ! fill up the sending/receiving buffer sizes (from pos_sends to pos_send, same for recv)
      ! set all to zero first (there are more buffer entrys than needed for reduced comm)

#ifndef USE_MPI 
#ifndef USE_CAF
   ! if no parallelism is used, set global = local parameters
   prc%bnd(1,1,1,1,1) = 1
   prc%bnd(1,1,1,1,2) = s_par%gx(1) 
   prc%bnd(1,1,1,2,1) = 1
   prc%bnd(1,1,1,2,2) = s_par%gx(2)
   prc%bnd(1,1,1,3,1) = 1
   prc%bnd(1,1,1,3,2) = s_par%gx(3)
   lb_dom%lx(1)       = s_par%gx(1)
   lb_dom%lx(2)       = s_par%gx(2)
   lb_dom%lx(3)       = s_par%gx(3)
#endif
#endif

#ifdef USE_MPI
   if(prc%rk == prc%root_th) then
#endif
!FIXME Caution: this calculation is still errorous. Check lb_dom%lx(1)=17 and use 10 threads
!               and 10 subdivisions in x-direction -> last two are not correct!!
      if(real(s_par%gx(1)) / real(prc%np(1)) < 2. .and. prc%np(1) > 1 ) & 
                     & write(*,*) "caution: Do not use that many threads in x-Direction"
      if(real(s_par%gx(2)) / real(prc%np(2)) < 2. .and. prc%np(2) > 1 ) &
                     & write(*,*)  "caution: Do not use that many threads in y-Direction"
      if(real(s_par%gx(3)) / real(prc%np(3)) < 2. .and. prc%np(3) > 1 ) &
                     & write(*,*) "caution: Do not use that many threads in z-Direction"

         ! calculate the boundaries for each subdomain
         do k = 1,prc%np(3)
            do j = 1,prc%np(2)
               do i = 1,prc%np(1)
                  prc%bnd(i,j,k,1,1) = (i-1) * CEILING(real(s_par%gx(1))/real(prc%np(1)))+1
                  prc%bnd(i,j,k,2,1) = (j-1) * CEILING(real(s_par%gx(2))/real(prc%np(2)))+1
                  prc%bnd(i,j,k,3,1) = (k-1) * CEILING(real(s_par%gx(3))/real(prc%np(3)))+1
                  if(i == prc%np(1) ) then
                     prc%bnd(i,j,k,1,2) = s_par%gx(1) 
                  else 
                     prc%bnd(i,j,k,1,2) = i*CEILING(real(s_par%gx(1))/real(prc%np(1)))
                  end if
                  if(j == prc%np(2) ) then
                     prc%bnd(i,j,k,2,2) = s_par%gx(2) 
                  else 
                     prc%bnd(i,j,k,2,2) = j*CEILING(real(s_par%gx(2))/real(prc%np(2)))
                  end if
                  if(k == prc%np(3) ) then
                     prc%bnd(i,j,k,3,2) = s_par%gx(3) 
                  else 
                     prc%bnd(i,j,k,3,2) = k*CEILING(real(s_par%gx(3))/real(prc%np(3)))
                  end if
               end do
            end do
         end do
#ifdef USE_MPI
      end if
      call MPI_BCAST(prc%bnd,prc%np(1)*prc%np(2)*prc%np(3)*3*2,MPI_INTEGER,0,prc%cart_comm,prc%ierr)
      if(s_par%problem==cavity_same) then
         ! if the complete cavity is calculated on each process (for benchmarking reasons),
         ! each process uses the globally defined domain size.
         prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,:,1) = 1
         prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1,2) = s_par%gx(1)
         prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,2,2) = s_par%gx(2)
         prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,3,2) = s_par%gx(3)
      end if
#endif 
#ifdef USE_CAF
      ! the boundary of all coarrays have to be the same! as the first domain is the largest > 
      ! assign to all coarrays
      cobnd(1) = prc%bnd(1,1,1,1,2) &
           &   - prc%bnd(1,1,1,1,1)+1 
      cobnd(2) = prc%bnd(1,1,1,2,2) &
           &   - prc%bnd(1,1,1,2,1)+1
      cobnd(3) = prc%bnd(1,1,1,3,2) & 
           &   - prc%bnd(1,1,1,3,1)+1
#endif
      ! get the size of each subdomain
      lb_dom%lx(1) = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1,2) &
            &   - prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1,1)+1 
      lb_dom%lx(2) = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,2,2) &
            &   - prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,2,1)+1
      lb_dom%lx(3) = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,3,2) & 
            &   - prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,3,1)+1

#ifdef VERBOSE
      call mpl_barrier()

      write(*,'(a,i3,a,i4,a,i4,a,i4,a,3i3)') & 
      & " Rank ",prc%rk,"  Domain size:  ",lb_dom%lx(1)," ", & 
      & lb_dom%lx(2)," ",lb_dom%lx(3),"  mpi-crd ",prc%crd

      call mpl_barrier()
#else
      if(prc%rk==0) then


      endif
#endif
      if(prc%rk == 0) write(*,*)


end subroutine mpl_divide_and_set_domain_size





subroutine mpl_create_buffers(lb_dom,prc)
!------------------------------------------------------------------------
!
   implicit none
   type(lb_block)      :: lb_dom
   type(mpl_var)       :: prc
   integer             :: lx,ly,lz 
   integer             :: neighbor
   integer             :: subsizes(4)
#if defined MPI_SUBARRAY
   integer             :: sizes(4),starts(4),s(3),counter
#endif      

      lx = lb_dom%lx(1)
      ly = lb_dom%lx(2)
      lz = lb_dom%lx(3)


      ! Assign the neighbors to prc%nghb
      ! orthogonal neighbors
      ! get the buffer length for each communication buffer     
      prc%length(:) = 0

      do neighbor=1,prc%nnghb
         prc%length(neighbor) = (prc%pos_send(neighbor,1) - prc%pos_sends(neighbor,1) + 1) &
                  &    * (prc%pos_send(neighbor,2) - prc%pos_sends(neighbor,2) + 1) &
                  &    * (prc%pos_send(neighbor,3) - prc%pos_sends(neighbor,3) + 1)

         prc%length_recv(neighbor) = (prc%pos_recv(neighbor,1) - prc%pos_recvs(neighbor,1) + 1) &
                  &    * (prc%pos_recv(neighbor,2) - prc%pos_recvs(neighbor,2) + 1) &
                  &    * (prc%pos_recv(neighbor,3) - prc%pos_recvs(neighbor,3) + 1)
      end do

#if defined MPI_SUBARRAY
      ! create the subarray derived mpi type for usage of isend irecv
   counter=0
   do neighbor=1,prc%nnghb
      ! receive buffer
#ifdef PERSISTENT2
      if (prc%nghb(neighbor) /= mpl_proc_null .and. prc%nghb(neighbor) .ge. 0) then 
      counter = counter+1
#else
      counter=neighbor
#endif
      s(1) = prc%pos_recv(counter,1) - prc%pos_recvs(counter,1) + 1
      s(2) = prc%pos_recv(counter,2) - prc%pos_recvs(counter,2) + 1
      s(3) = prc%pos_recv(counter,3) - prc%pos_recvs(counter,3) + 1
      sizes    = (/NDX(nnod, lx+2, ly+2, lz+2) /)
      subsizes = (/NDX(nnod, s(1),s(2),s(3)) /) 
      starts   = (/NDX(0, prc%pos_recvs(counter,1), prc%pos_recvs(counter,2),prc%pos_recvs(counter,3))/)
      prc%length(counter) = subsizes(1)*subsizes(2)*subsizes(3)*subsizes(4)

      call mpi_type_create_subarray(4,sizes,subsizes,starts,      &
         &  mpi_order_fortran,mpi_double_precision,recv_dat(counter),prc%ierr)
      call mpi_type_commit(recv_dat(counter),prc%ierr)

      ! send buffer
      s(1) = prc%pos_send(counter,1) - prc%pos_sends(counter,1) + 1
      s(2) = prc%pos_send(counter,2) - prc%pos_sends(counter,2) + 1
      s(3) = prc%pos_send(counter,3) - prc%pos_sends(counter,3) + 1
      subsizes = (/NDX(nnod,s(1),s(2),s(3) ) /) 
      starts   = (/NDX(0, prc%pos_sends(counter,1), prc%pos_sends(counter,2),prc%pos_sends(counter,3))/)
      call mpi_type_create_subarray(4,sizes,subsizes,starts,&
         &  mpi_order_fortran,mpi_double_precision,send_dat(counter),prc%ierr)
      call mpi_type_commit(send_dat(counter),prc%ierr)
#ifdef PERSISTENT2
   endif
#endif

     end do

#else /* if not MPI_SUBARRAY #if defined MPI_SENDRECV_BUF || defined MPI_SENDRECV_ALL || defined MPI_ISEND_IRECV_BUF */
   ! create send and receive buffers
   ! This loop identifies, based on the bitmap test of  
   ! which densities to copy into comm buffer

   do neighbor = 1, prc%nnghb
         
         subsizes = (/ prc%ndir_pdf(neighbor),                                 & !nnod,&
                      (prc%pos_recv(neighbor,1) - prc%pos_recvs(neighbor,1) + 1),     &
                      (prc%pos_recv(neighbor,2) - prc%pos_recvs(neighbor,2) + 1),     &
                      (prc%pos_recv(neighbor,3) - prc%pos_recvs(neighbor,3) + 1) /) 
         prc%length(neighbor) = subsizes(1)*subsizes(2)*subsizes(3)*subsizes(4)

#if defined USE_CAF 
#if defined CAF_DER_TYPE
      allocate(mpl_buffer(neighbor)%send(subsizes(1),subsizes(2),subsizes(3),subsizes(4)))
      allocate(mpl_buffer(neighbor)%recv(subsizes(1),subsizes(2),subsizes(3),subsizes(4)))
#else
      if(neighbor == 1) then
      allocate(send1(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      allocate(recv1(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      elseif(neighbor == 2 ) then
      allocate(send2(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      allocate(recv2(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      elseif(neighbor == 3 ) then
      allocate(send3(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      allocate(recv3(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      elseif(neighbor == 4 ) then
      allocate(send4(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      allocate(recv4(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      elseif(neighbor == 5 ) then
      allocate(send5(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      allocate(recv5(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      elseif(neighbor == 6 ) then
      allocate(send6(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      allocate(recv6(subsizes(1),subsizes(2),subsizes(3),subsizes(4))[*])
      endif
#endif
#endif
#if defined USE_MPI 
      allocate(mpl_buffer(neighbor)%send(subsizes(1),subsizes(2),subsizes(3),subsizes(4)))
      allocate(mpl_buffer(neighbor)%recv(subsizes(1),subsizes(2),subsizes(3),subsizes(4)))
#endif

    end do
#endif



end subroutine mpl_create_buffers
!------------------------------------------------------------------------






!------------------------------------------------------------------------
subroutine mpl_assign_neighbors(lb_dom,prc)
!
!
   type(mpl_var)     :: prc
   integer           :: target_rank, ii,jj,lc(3),istart,iend,zero_pos
   type(lb_block)    :: lb_dom
   integer           :: ind, rel_pos,ghostn 
   integer           :: ll
   integer           :: xs(3),ys(3),zs(3) 
   integer           :: xr(3),yr(3),zr(3) 
   integer           :: counter,counterrecv
   logical           :: added_send,added_recv

      prc%nghb(:)     = mpl_proc_null
      prc%nnghb_avail = 0 
      ! fill up the sending/receiving buffer sizes (from pos_sends to pos_send, same for recv)
      ! set all to zero first (there are more buffer entries than needed for reduced comm)

!      allocate(prc%ndir_pdf( prc%nnghb))
      allocate(prc%dsend_pdf(prc%nnghb,nnod))
      allocate(prc%drecv_pdf(prc%nnghb,nnod))


      ! reduced mpi communication
      ! set all pdf directions to zero
      prc%drecv_pdf = 0
      prc%dsend_pdf = 0

      prc%pos_sends = 0
      prc%pos_recvs = 0 
      prc%pos_send  = 0
      prc%pos_recv  = 0
#ifdef BUFFER_EACH
   ! when each buffer is exchanged separately, ghost nodes don't have to be communicated
   ghostn = 0
#else /* ORTHOGONAL DIRECTIONS */
   ! orthogonal direction exchange needs ghost node values for writing to corners 
   ghostn = 1
#endif

!      if(prc%np(1)*prc%np(2)*prc%np(3) > 1) then
! where is the zero-vector?
zero_pos=0
do ii=1,nnod
if(cx(ii,1) ==0 .and. cx(ii,2)==0 .and. cx(ii,3)==0) zero_pos=ii
enddo

! this here should maybe performed via the cx() vectors
rel_pos=0
if(zero_pos==1) then
    istart=2; iend=prc%nnghb+1
    rel_pos=-1
else
    istart=1; iend=prc%nnghb
endif

!--------------------------------
! now do the loop over all neighbors
do ii=istart,iend
    ind=ii+rel_pos
    target_rank = -2

   ! get the coordinates of the current neighbor 
    lc(1) = prc%crd(1)+cx(ii,1)
    lc(2) = prc%crd(2)+cx(ii,2)
    lc(3) = prc%crd(3)+cx(ii,3)

   ! if it's a valid coordinate, get rank information
   if(lc(1) .ge. 0 .and. lc(1) .lt. prc%np(1)) then
   if(lc(2) .ge. 0 .and. lc(2) .lt. prc%np(2)) then
   if(lc(3) .ge. 0 .and. lc(3) .lt. prc%np(3)) then
#if defined USE_CAF
       target_rank = image_index(prc%caf_cart_comm,lc+1) 
    prc%ng_crd(:,ind) = lc+1
#else
#if defined USE_MPI
       call MPI_CART_RANK(prc%cart_comm,lc,target_rank,prc%ierr)
#endif
#endif
   prc%nnghb_avail = prc%nnghb_avail + 1
   ! assign
   prc%nghb(ind) = target_rank
   endif
   endif
   endif

! maybe we should right away put the buffer node information right here

   ! treat x direction
   if(cx(ii,1) ==0) then
      xs(1) = 1  - ghostn
      xs(2) = lb_dom%lx(1) + ghostn
      xr(1) = 1  - ghostn
      xr(2) = lb_dom%lx(1) + ghostn
   elseif(cx(ii,1)==1) then
      xs(1:2) = lb_dom%lx(1)
      xr(1:2) = lb_dom%lx(1)+1
   else ! -1
      xs(1:2) = 1
      xr(1:2) = 0
   endif

   ! treat y direction
   if(cx(ii,2) ==0) then
      ys(1) = 1  - ghostn
      ys(2) = lb_dom%lx(2) + ghostn
      yr(1) = 1  - ghostn
      yr(2) = lb_dom%lx(2) + ghostn
   elseif(cx(ii,2)==1) then
      ys(1:2) = lb_dom%lx(2)
      yr(1:2) = lb_dom%lx(2)+1
   else ! -1
      ys(1:2) = 1
      yr(1:2) = 0
   endif

   ! treat z direction
   if(cx(ii,3) ==0) then
      zs(1) = 1  - ghostn
      zs(2) = lb_dom%lx(3) + ghostn
      zr(1) = 1  - ghostn
      zr(2) = lb_dom%lx(3) + ghostn
   elseif(cx(ii,3)==1) then
      zs(1:2) = lb_dom%lx(3)
      zr(1:2) = lb_dom%lx(3)+1
   else ! -1
      zs(1:2) = 1
      zr(1:2) = 0
   endif


   prc%pos_sends(ind,1:3) = (/xs(1),  ys(1), zs(1)/)
   prc%pos_send (ind,1:3) = (/xs(2),  ys(2), zs(2)/)
   prc%pos_recvs(ind,1:3) = (/xr(1),  yr(1), zr(1)/)
   prc%pos_recv (ind,1:3) = (/xr(2),  yr(2), zr(2)/)

   ! which densities have to be exchanged??
   counter=0
   counterrecv=0
   do ll=1,nnod
#ifdef COMM_REDUCED
   added_send=.false.
   added_recv=.false.
      do jj=1,3
      if(cx(ll,jj) == cx(ii,jj) .and. abs(cx(ll,jj)) .gt. 0) then
         if(added_send .eqv. .false.) then
         counter=counter+1
         prc%dsend_pdf(ind,counter)=ll
         added_send = .true.
         endif
      endif 
      if(cx(ll,jj) == -cx(ii,jj) .and. abs(cx(ll,jj)) .gt. 0) then
         if(added_recv .eqv. .false.) then
         counterrecv=counterrecv+1
         prc%drecv_pdf(ind,counterrecv)=ll
         added_recv = .true.
         endif
      endif 
       enddo
#else
   counter=counter+1
   prc%dsend_pdf(ind,ll)   = ll
   prc%drecv_pdf(ind,ll)   = ll
#endif /* COMM_REDUCED */
   enddo
   prc%ndir_pdf(ind) = counter

#ifdef DEBUG_MPI
   write(66+prc%rk,*) "Index",      ind,"  dir",cx(ii,:)
   write(66+prc%rk,*) prc%pos_sends(ind,:) 
   write(66+prc%rk,*) prc%pos_send( ind,:) 
   write(66+prc%rk,*) prc%pos_recvs(ind,:) 
   write(66+prc%rk,*) prc%pos_recv( ind,:) 
   write(66+prc%rk,*) 
   write(66+prc%rk,*) "densities",counter
   write(66+prc%rk,*) "send",prc%dsend_pdf(ind,:) 
   write(66+prc%rk,*) 
   write(66+prc%rk,*) "recv",prc%drecv_pdf(ind,:) 
   write(66+prc%rk,*) 
   write(66+prc%rk,*) 
#endif /* DEBUG_MPI */



enddo

#ifdef DEBUG_MPI
   write(*,*) prc%rk,"neighbors",prc%nghb
#endif




#ifdef USE_CAF
!-----------------------------------
! build team (get neighbors) and include invoking image

prc%team=0
counter  = 0 
do i=1,prc%nnghb
   if(prc%nghb(i) > 0) then
      counter=counter+1
 !     prc%team(counter) = prc%nghb(i)
   endif
enddo
!allocate(prc%team(counter))
allocate(g_team(counter+1))
! invoking image is part of team
counter=1
g_team(counter) = prc%rk+1 
do i=1,prc%nnghb
   if(prc%nghb(i) > 0) then
      counter=counter+1
!      prc%team(counter) = prc%nghb(i)
      g_team(counter) = prc%nghb(i)
   endif
enddo
#endif



end subroutine mpl_assign_neighbors
!------------------------------------------------------------------------





#ifdef USE_ADCL
!------------------------------------------------------------------------
subroutine mpl_adcl_init(lb_dom,s_par,prc) 
   use LBMODEL
   type(lb_block),intent(inout)     :: lb_dom
   type(sim_parameter),intent(in)   :: s_par
   type(mpl_var)                    :: prc
   integer                     :: dims(4)
   integer                     :: lx,ly,lz
#ifdef LAYOUT_LIJK
   integer                     :: allocstat
#endif

   lx = lb_dom%lx(1)
   ly = lb_dom%lx(2)
   lz = lb_dom%lx(3)

      dims = (/lx+2,ly+2,lz+2,nnod/)
#ifdef LAYOUT_LIJK
      allocate(prc%adcl_data(0:lx+1,0:ly+1,0:lz+1,nnod),stat=allocstat)
      if(allocstat/=0) write(*,*) "Error allocating prc%adcl_data"
      call ADCL_Vector_register_generic(3,dims,nnod,prc%adcl_vmap,mpi_double_precision,&
                                       & prc%adcl_data,prc%adcl_vec,prc%ierr) 
#endif
#ifdef LAYOUT_IJKL
      call ADCL_Vector_register_generic(3,dims,nnod,prc%adcl_vmap,mpi_double_precision,&
                                       & lb_dom%fIn,prc%adcl_vec_fin,prc%ierr) 
      call ADCL_Vector_register_generic(3,dims,nnod,prc%adcl_vmap,mpi_double_precision,&
                                       & lb_dom%fOut,prc%adcl_vec_fout,prc%ierr) 
#endif
      call ADCL_Topology_create_extended(prc%cart_comm, prc%adcl_topo, prc%ierr )
      call ADCL_Request_create(prc%adcl_vec_fIn,prc%adcl_topo,ADCL_FNCTSET_NEIGHBORHOOD, prc%adcl_request_fIn, prc%ierr)
      call ADCL_Request_create(prc%adcl_vec_fOut,prc%adcl_topo,ADCL_FNCTSET_NEIGHBORHOOD, prc%adcl_request_fOut, prc%ierr)

end subroutine mpl_adcl_init
!------------------------------------------------------------------------
#endif /*USE_ADCL*/








subroutine mpl_finish(prc) 
!------------------------------------------------------------------------
! finish up parallel communication
! - ADCL and MPI finishing routines are called, 
! - files are closed and 
! - buffer arrays are deallocated

  type(mpl_var)   :: prc
  integer        :: ierr,i,neighbor


#ifdef USE_MPI
#ifdef USE_ADCL
   ! finalize ADCL
   call ADCL_Request_free ( prc%adcl_request_fIn, ierr )
   call ADCL_Request_free ( prc%adcl_request_fOut, ierr )
   call ADCL_Vector_deregister ( prc%adcl_vec_fIn, ierr )
   call ADCL_Vector_deregister ( prc%adcl_vec_fOut, ierr )
   call ADCL_Vmap_free ( prc%adcl_vmap, ierr )
   call ADCL_Topology_free ( prc%adcl_topo, ierr )
   call ADCL_finalize( ierr )
#endif

#if defined(PERSISTENT2)
   ! kill persistent mpi requests
   do neighbor = 1,prc%nnghb_avail
      call mpi_request_free(prc%ireq(neighbor  ),prc%ierr) 
      call mpi_request_free(prc%ireq(prc%nnghb_avail+neighbor  ),prc%ierr)
   enddo
#endif

  call MPI_COMM_FREE(prc%cart_comm,prc%ierr)

  call MPI_FINALIZE(prc%ierr)

#endif
   deallocate(prc%bnd,stat=ierr)
   deallocate(prc%dsend_pdf,stat=ierr)
   deallocate(prc%drecv_pdf,stat=ierr)

   do i = 1,prc%nnghb
#ifndef USE_CAF 
      deallocate(mpl_buffer(i)%send,stat=ierr)
      deallocate(mpl_buffer(i)%recv,stat=ierr)
#endif
#ifdef CAF_DER_TYPE
      deallocate(mpl_buffer(i)%send,stat=ierr)
      deallocate(mpl_buffer(i)%recv,stat=ierr)
#endif
#if defined MPI_NOT_ALL_DENS
      deallocate(mpl_buffer(i)%idx,stat=ierr)
      deallocate(mpl_buffer(i)%idx_recv,stat=ierr)
#endif
   end do
#if defined USE_CAF
#ifndef CAF_DER_TYPE
      deallocate(send1,stat=ierr)
      deallocate(recv1,stat=ierr)
      deallocate(send2,stat=ierr)
      deallocate(recv2,stat=ierr)
      deallocate(send3,stat=ierr)
      deallocate(recv3,stat=ierr)
      deallocate(send4,stat=ierr)
      deallocate(recv4,stat=ierr)
      deallocate(send5,stat=ierr)
      deallocate(recv5,stat=ierr)
      deallocate(send6,stat=ierr)
      deallocate(recv6,stat=ierr)
#endif
!   deallocate(prc%team,stat=ierr)
   deallocate(g_team,stat=ierr)
   deallocate(prc%caf_cart_comm,stat=ierr)
#endif 

#ifdef WRITE_TSTATUS
   close(42)
#endif


if(prc%rk==0) then
write(*,*) 
write(*,*) "     ... done!"
write(*,*) 
endif
end subroutine mpl_finish 
!------------------------------------------------------------------------




#ifdef USE_CAF


  !------------------------------------------------------------------------
  !  This is the CAF COmmunication routine
  !------------------------------------------------------------------------



  subroutine mpl_exchange(lb_dom,fIn,prc,meas)
  !------------------------------------------------------------------------
  !
  ! Data from the border cells are written to the ghost cells and vice versa
  !
  ! Note: This is the exchange routine, where diagonal communication is done 
  !       via two orthogonal communcations. This reduces the overall communcations to six.
  !       Diagonally related domains interchange the node values
  !       by first doing the communication in the first orthogonal direction, 
  !       then do a sync, and then exchange in the other orthogonal direction.
  !

      type(lb_block),intent(inout) :: lb_dom
      type(mpl_var)  :: prc
      type(measure)  :: meas
      integer        :: ii,jj,kk,i,j,k,l,lx,ly,lz
      integer        :: counter,i_req(2*prc%nnghb)
      real(R8B)      :: copy_start,copy_end,comm_start,comm_end
      integer        :: next
      real(R8B)      :: sync_start,sync_end
      real(R8B)      :: copyb_st,copyb_end,copied_amount 
#ifdef OLD_CAF
      real(R8B) :: fIn(NDX(nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))[prc%np(1),prc%np(2),*]
#else
      real(R8B) :: fIn(NDX(nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))
#endif


#ifdef CAF_BENCH

!-------------------------------------------
!  small CoArray Fortran bandwidth benchmark
!  with differently allocated arrays
!

   real(R8B),allocatable,dimension(:,:,:,:)    :: test1
   real(R8B),allocatable,dimension(:,:,:,:)[:] :: test2 
   integer :: s1,s2,s3,s4 

   s1=size(send1,1)
   s2=size(send1,2)
   s3=size(send1,3)
   s4=size(send1,4)
   if(prc%rk==0) then
      write(*,*) 
      write(*,*) "CAF BENCHMARK"
      write(*,*) 
   endif
   ! Receive Array
   allocate(test1(s1,s2,s3,s4))
#if CAF_BENCH>1
   if(prc%rk==0)       write(*,*) "CAF-Arrays: locally defined allocatable array"

   ! Send Array
   allocate(test2(s1,s2,s3,s4)[*])
   sync all
   test2=real(prc%rk+1)
#else
   if(prc%rk==0)       write(*,*) "CAF-Arrays: globally defined send buffer" 

   ! Send Array
   send1=real(prc%rk+1)
#endif

   sync all
   call cpu_time_measure(meas%tStart)


#if CAF_BENCH>1
   test1(:,:,:,:) = test2(:,:,:,:)[2]
#else
   if(prc%rk==0)    test1(:,:,:,:) = send1(:,:,:,:)[2]
#endif

   sync all 
   call cpu_time_measure(meas%tEnd)

   write(*,*) prc%rk+1,"Res         ",test1(1,1,1,1)
if(prc%rk==0) then
   write(*,*) "duration was",meas%tEnd - meas%tStart
   write(*,*) "Buf size    ",real(size(test1,1)*size(test1,2)*size(test1,3)*size(test1,4)*8)
   write(*,*) "BW          ",real(size(test1,1)*size(test1,2)*size(test1,3)*size(test1,4)*8)&
/(meas%tEnd - meas%tStart)/1000000.
endif
sync all


   deallocate(test1)
#if CAF_BENCH>1
   deallocate(test2)
#endif

stop

! End of CAF Benchmark
!-------------------------------------------

#endif /* CAF_BENCH*/



#ifdef SYNCHRONIZE_COMM 
   call mpl_barrier()
#endif

   call cpu_time_measure(meas%tSt_comm) 

   lx = lb_dom%lx(1)
   ly = lb_dom%lx(2)
   lz = lb_dom%lx(3)


   do i = 1, prc%nnghb,2
      next = i+1

!if(i==1) then
!   send_pnti => send1
!   send_pntn => send2
!   recv_pnti => recv1
!   recv_pntn => recv2
!elseif(i==3) then
!   send_pnti => send3
!   send_pntn => send4
!   recv_pnti => recv3
!   recv_pntn => recv4
!elseif(i==5) then
!   send_pnti => send5
!   send_pntn => send6
!   recv_pnti => recv5
!   recv_pntn => recv6
!endif

      ! do loop over the six cartesian neighbors. 
      ! diagonal speeds are communicated by recopying the buffer 
      ! after each communication.
      call cpu_time_measure(copy_start)

      ! is there a neighbor in x- / y- / z- direction?
      if(prc%nghb(i) /= mpl_proc_null .and. prc%nghb(i) .ge. 0) then

      do kk = 1, prc%pos_send(i,3)-prc%pos_sends(i,3)+1
         do jj = 1, prc%pos_send(i,2)-prc%pos_sends(i,2)+1
            do ii = 1, prc%pos_send(i,1)-prc%pos_sends(i,1)+1
               do l = 1,prc%ndir_pdf
#ifdef CAF_DER_TYPE
                  mpl_buffer(i)%send(l,ii,jj,kk)  =   & 
fIn(NDX(prc%dsend_pdf(i,l),ii+prc%pos_sends(i,1)-1, jj+prc%pos_sends(i,2)-1, kk+prc%pos_sends(i,3)-1))
#else
!send_pnti(l,ii,jj,kk) = &
!fIn(NDX(prc%dsend_pdf(i,l),ii+prc%pos_sends(i,1)-1, jj+prc%pos_sends(i,2)-1, kk+prc%pos_sends(i,3)-1))
if(i==1) then
                  send1(l,ii,jj,kk) = & 
fIn(NDX(prc%dsend_pdf(i,l),ii+prc%pos_sends(i,1)-1, jj+prc%pos_sends(i,2)-1, kk+prc%pos_sends(i,3)-1))
elseif(i==3) then
                  send3(l,ii,jj,kk) = & 
fIn(NDX(prc%dsend_pdf(i,l),ii+prc%pos_sends(i,1)-1, jj+prc%pos_sends(i,2)-1, kk+prc%pos_sends(i,3)-1))
elseif(i==5) then
                  send5(l,ii,jj,kk) = & 
fIn(NDX(prc%dsend_pdf(i,l),ii+prc%pos_sends(i,1)-1, jj+prc%pos_sends(i,2)-1, kk+prc%pos_sends(i,3)-1))
endif
#endif
                 enddo
            enddo
         enddo
      enddo
      endif

      if(prc%nghb(next) /= mpl_proc_null .and. prc%nghb(next) .ge. 0) then
      do kk = 1, prc%pos_send(next,3)-prc%pos_sends(next,3)+1
         do jj = 1, prc%pos_send(next,2)-prc%pos_sends(next,2)+1
            do ii = 1, prc%pos_send(next,1)-prc%pos_sends(next,1)+1
               do l = 1,prc%ndir_pdf
#ifdef CAF_DER_TYPE
                  mpl_buffer(next)%send(l,ii,jj,kk) =    & 
fIn(NDX(prc%dsend_pdf(next,l), ii+prc%pos_sends(next,1)-1, jj+prc%pos_sends(next,2)-1, kk+prc%pos_sends(next,3)-1))
#else
!send_pntn(l,ii,jj,kk) = &
!fIn(NDX(prc%dsend_pdf(next,l), ii+prc%pos_sends(next,1)-1, jj+prc%pos_sends(next,2)-1, kk+prc%pos_sends(next,3)-1))
if(i==1) then
                  send2(l,ii,jj,kk) = & 
fIn(NDX(prc%dsend_pdf(next,l), ii+prc%pos_sends(next,1)-1, jj+prc%pos_sends(next,2)-1, kk+prc%pos_sends(next,3)-1))
elseif(i==3) then
                  send4(l,ii,jj,kk) = & 
fIn(NDX(prc%dsend_pdf(next,l), ii+prc%pos_sends(next,1)-1, jj+prc%pos_sends(next,2)-1, kk+prc%pos_sends(next,3)-1))
elseif(i==5) then
                  send6(l,ii,jj,kk) = & 
fIn(NDX(prc%dsend_pdf(next,l), ii+prc%pos_sends(next,1)-1, jj+prc%pos_sends(next,2)-1, kk+prc%pos_sends(next,3)-1))
endif
#endif
                 enddo
            enddo
         enddo
      enddo
      endif
      call cpu_time_measure(copy_end)

     

!---------------------------------------
! Do  Exchange         

      call cpu_time_measure(sync_start)
!      call sync_images(prc%team)
      call sync_images(g_team)
      call cpu_time_measure(sync_end)
 
      call cpu_time_measure(comm_start)

      if(prc%nghb(i) /= mpl_proc_null .and. prc%nghb(i) .ge. 0) then
#ifdef CAF_DER_TYPE
          mpl_buffer(i)%recv(:,:,:,:) =       mpl_buffer(next)[prc%nghb(i)]%send(:,:,:,:)
#else
if(i==1) then
          recv1(:,:,:,:) = send2(:,:,:,:)[prc%nghb(i)]
elseif(i==3) then
          recv3(:,:,:,:) = send4(:,:,:,:)[prc%nghb(i)]
elseif(i==3) then
          recv5(:,:,:,:) = send6(:,:,:,:)[prc%nghb(i)]
endif
          call sync_memory()
#endif
      endif
      if(prc%nghb(next) /= mpl_proc_null .and. prc%nghb(next) .ge. 0) then
#ifdef CAF_DER_TYPE
          mpl_buffer(next)%recv(:,:,:,:) = mpl_buffer(i)[prc%nghb(next)]%send(:,:,:,:)
#else
if(i==1) then
          recv2(:,:,:,:) = send1(:,:,:,:)[prc%nghb(next)]
elseif(i==3) then
          recv4(:,:,:,:) = send3(:,:,:,:)[prc%nghb(next)]
elseif(i==5) then
          recv6(:,:,:,:) = send5(:,:,:,:)[prc%nghb(next)]
endif
          call sync_memory()
#endif
      endif
   call cpu_time_measure(comm_end)


      call cpu_time_measure(sync_start)
!      call sync_images(prc%team)
      call sync_images(g_team)
      call cpu_time_measure(sync_end)
      meas%mpl_sync = meas%mpl_sync +sync_end - sync_start


!-------------------------------------
! Copy back from buffer 
 
      call cpu_time_measure(copyb_st) 
      
      if(prc%nghb(i) /= mpl_proc_null .and. prc%nghb(i) .ge. 0) then
      meas%mpl_comm_size   = meas%mpl_comm_size   + real(prc%length(i)*R8B)
      do kk = 1, prc%pos_recv(i,3)-prc%pos_recvs(i,3)+1
         do jj = 1, prc%pos_recv(i,2)-prc%pos_recvs(i,2)+1
            do ii = 1, prc%pos_recv(i,1)-prc%pos_recvs(i,1)+1
               do l = 1,prc%ndir_pdf
#ifdef CAF_DER_TYPE
fIn(NDX(prc%drecv_pdf(i,l),ii+prc%pos_recvs(i,1)-1,jj+prc%pos_recvs(i,2)-1, kk+prc%pos_recvs(i,3)-1)) & 
            =  mpl_buffer(i)%recv(l,ii,jj,kk)  
#else
!fIn(NDX(prc%drecv_pdf(i,l),ii+prc%pos_recvs(i,1)-1,jj+prc%pos_recvs(i,2)-1, kk+prc%pos_recvs(i,3)-1)) & 
!               =  recv_pnti(l,ii,jj,kk) 
if(i==1) then
fIn(NDX(prc%drecv_pdf(i,l),ii+prc%pos_recvs(i,1)-1,jj+prc%pos_recvs(i,2)-1, kk+prc%pos_recvs(i,3)-1)) & 
               =  recv1(l,ii,jj,kk) 
elseif(i==3) then
fIn(NDX(prc%drecv_pdf(i,l),ii+prc%pos_recvs(i,1)-1,jj+prc%pos_recvs(i,2)-1, kk+prc%pos_recvs(i,3)-1)) & 
               =  recv3(l,ii,jj,kk) 
elseif(i==5) then
fIn(NDX(prc%drecv_pdf(i,l),ii+prc%pos_recvs(i,1)-1,jj+prc%pos_recvs(i,2)-1, kk+prc%pos_recvs(i,3)-1)) & 
               =  recv5(l,ii,jj,kk) 
endif
#endif
               enddo
            enddo
         enddo
      enddo
      endif

      if(prc%nghb(next) /= mpl_proc_null .and. prc%nghb(next) .ge. 0) then
      meas%mpl_comm_size   = meas%mpl_comm_size   + real(prc%length(next)*R8B)
      do kk = 1, prc%pos_recv(next,3)-prc%pos_recvs(next,3)+1
         do jj = 1, prc%pos_recv(next,2)-prc%pos_recvs(next,2)+1
            do ii = 1, prc%pos_recv(next,1)-prc%pos_recvs(next,1)+1
               do l = 1,prc%ndir_pdf
#ifdef CAF_DER_TYPE
fIn(NDX(prc%drecv_pdf(next,l),ii+prc%pos_recvs(next,1)-1,jj+prc%pos_recvs(next,2)-1, kk+prc%pos_recvs(next,3)-1))   & 
            =  mpl_buffer(next)%recv(l,ii,jj,kk)  
#else
!fIn(NDX(prc%drecv_pdf(next,l),ii+prc%pos_recvs(next,1)-1,jj+prc%pos_recvs(next,2)-1, kk+prc%pos_recvs(next,3)-1))   & 
!               =  recv_pntn(l,ii,jj,kk) 
if(i==1) then
fIn(NDX(prc%drecv_pdf(next,l),ii+prc%pos_recvs(next,1)-1,jj+prc%pos_recvs(next,2)-1, kk+prc%pos_recvs(next,3)-1))   & 
               =  recv2(l,ii,jj,kk) 
elseif(i==3) then
fIn(NDX(prc%drecv_pdf(next,l),ii+prc%pos_recvs(next,1)-1,jj+prc%pos_recvs(next,2)-1, kk+prc%pos_recvs(next,3)-1))   & 
               =  recv4(l,ii,jj,kk) 
elseif(i==5) then
fIn(NDX(prc%drecv_pdf(next,l),ii+prc%pos_recvs(next,1)-1,jj+prc%pos_recvs(next,2)-1, kk+prc%pos_recvs(next,3)-1))   & 
               =  recv6(l,ii,jj,kk) 
endif
#endif
               enddo
            enddo
         enddo
      enddo
      endif
      call cpu_time_measure(copyb_end) 

      meas%mpl_cp_buf      = meas%mpl_cp_buf      + copy_end - copy_start
      meas%mpl_cp_buf_back = meas%mpl_cp_buf_back + copyb_end - copyb_st
      meas%mpl_exch        = meas%mpl_exch        + comm_end - comm_start
      enddo

!---------------------------------------
! End of CAF Communication 




#if defined USE_CAF && defined OLD_CAF
!-----------------------------------------------------------------------
! Coarray communication    OLD 
!-----------------------------------------------------------------------
! 
! 
   call cpu_time_measure(comm_start)
   call sync_images(g_team)
   call cpu_time_measure(comm_end)
   meas%mpl_sync = meas%mpl_sync + comm_end - comm_start

!------------------------------
! x+/-

   if(prc%crd(1) < prc%np(1) - 1) then  ! x+
      call cpu_time_measure(comm_start) 
      fIn(NDX(1:nnod,lx+1,0:ly+1,0:lz+1))  & 
    & = fIn(NDX(1:nnod,1,   0:ly+1,0:lz+1)) [prc%crd(1)+1+1,prc%crd(2)+1,prc%crd(3)+1]
      call cpu_time_measure(comm_end) 
      meas%mpl_exch = meas%mpl_exch + comm_end - comm_start
   end if
   if(prc%crd(1) > 0 ) then             ! x-
      call cpu_time_measure(comm_start) 
      fIn(NDX(1:nnod,0   ,0:ly+1,0:lz+1))  & 
      = fIn(NDX(1:nnod,lx,  0:ly+1,0:lz+1))[prc%crd(1)-1 +1,prc%crd(2)+1,prc%crd(3)+1]
      call cpu_time_measure(comm_end) 
      meas%mpl_exch = meas%mpl_exch + comm_end - comm_start
   end if

   call cpu_time_measure(comm_start) 
   call sync_images(g_team)
   call cpu_time_measure(comm_end) 
   meas%mpl_sync = meas%mpl_sync + comm_end - comm_start


!------------------------------
! y+/-

  if(prc%crd(2) < prc%np(2) - 1) then  ! y+
      call cpu_time_measure(comm_start) 
      fIn(NDX(1:nnod,0:lx+1,ly+1,0:lz+1))  & 
      = fIn(NDX(1:nnod,0:lx+1,1   ,0:lz+1))[prc%crd(1)+1,prc%crd(2)+1+1,prc%crd(3)+1]
      call cpu_time_measure(comm_end) 
      meas%mpl_exch = meas%mpl_exch + comm_end - comm_start
   end if
   if(prc%crd(2) > 0 ) then            ! y-
      call cpu_time_measure(comm_start) 
      fIn(NDX(1:nnod,0:lx+1,0,0:lz+1))     & 
      = fIn(NDX(1:nnod,0:lx+1,ly  ,0:lz+1))[prc%crd(1)+1,prc%crd(2)+1-1,prc%crd(3)+1]
      call cpu_time_measure(comm_end) 
      meas%mpl_exch = meas%mpl_exch + comm_end - comm_start
   end if

   call cpu_time_measure(comm_start) 
   call sync_images(g_team)
   call cpu_time_measure(comm_end) 
   meas%mpl_sync = meas%mpl_sync + comm_end - comm_start

!------------------------------
! z+/-

   if(prc%crd(3) < prc%np(3) - 1) then  ! z+
      call cpu_time_measure(comm_start) 
      fIn(NDX(1:nnod,0:lx+1,0:ly+1,lz+1))  & 
      = fIn(NDX(1:nnod,0:lx+1,0:ly+1,1   ))[prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1+1]
      call cpu_time_measure(comm_end) 
      meas%mpl_exch = meas%mpl_exch + comm_end - comm_start
   end if
   if(prc%crd(3) > 0 ) then            ! z-
      call cpu_time_measure(comm_start) 
      fIn(NDX(1:nnod,0:lx+1,0:ly+1,0))  & 
     = fIn(NDX(1:nnod,1:lx,1:ly,lz  )) [prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1-1]
      call cpu_time_measure(comm_end) 
      meas%mpl_exch = meas%mpl_exch + comm_end - comm_start
   end if

!----------------------------------------------------------------------------

   call cpu_time_measure(comm_start) 
   call sync_images(g_team)
   call cpu_time_measure(comm_end) 
   meas%mpl_sync = meas%mpl_sync + comm_end - comm_start

   call cpu_time_measure(comm_end) 
#endif /* USE_CAF */


   call cpu_time_measure(meas%tEnd_comm)
   meas%comm_duration = meas%tEnd_comm - meas%tSt_comm + meas%comm_duration

  end subroutine mpl_exchange
!------------------------------------------------------------------------








#else /* not USE_CAF */





  !------------------------------------------------------------------------
  !  This is the MPI Communication routine
  !------------------------------------------------------------------------


  subroutine mpl_exchange(lb_dom,fIn,prc,meas)
  !------------------------------------------------------------------------
  !
  ! Data from the border cells are written to the ghost cells and vice versa
  !
  ! Note: This is the exchange routine, where diagonal communication is done 
  !       via two orthogonal communcations. This reduces the overall communcations to six.
  !       Diagonally related domains interchange the node values
  !       by first doing the communication in the first orthogonal direction, 
  !       then do a sync, and then exchange in the other orthogonal direction.
  !

      type(lb_block),intent(inout) :: lb_dom
      type(mpl_var)  :: prc
      type(measure)  :: meas
      integer        :: neighbor
      real(R8B) :: fIn(NDX(nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))

#ifdef SYNCHRONIZE_COMM 
   call mpl_barrier()
#endif

   call cpu_time_measure(meas%tSt_comm) 


#if defined USE_MPI

#ifdef USE_ADCL

!-----------------------------------------------------------------------
! ADCL Communication 
!-----------------------------------------------------------------------

#ifdef LAYOUT_LIJK
   write(*,*) "For ADCL, please use Layout IJKL for good performance. Stopping..."
   stop
#endif /* LAYOUT_LIJK */

!----------------------------------
! ADCL data exchange (ADCL Request) 

   call cpu_time_measure(comm_start) 
   if ( prc%adcl_active_req == 0 ) then 
      call ADCL_Request_start(prc%adcl_request_fIn, prc%ierr)
   else
      call ADCL_Request_start(prc%adcl_request_fOut, prc%ierr)
   end if
   call cpu_time_measure(comm_end) 
   meas%mpl_exch = meas%mpl_exch + comm_end - comm_start

   ! Measure the size of the exchange-buffers
   if(prc%nghb(1) /= mpl_proc_null) meas%mpl_comm_size = meas%mpl_comm_size   + real(prc%length(1)*R8B)
   if(prc%nghb(2) /= mpl_proc_null) meas%mpl_comm_size = meas%mpl_comm_size   + real(prc%length(2)*R8B)
   if(prc%nghb(3) /= mpl_proc_null) meas%mpl_comm_size = meas%mpl_comm_size   + real(prc%length(3)*R8B)
   if(prc%nghb(4) /= mpl_proc_null) meas%mpl_comm_size = meas%mpl_comm_size   + real(prc%length(4)*R8B)
   if(prc%nghb(5) /= mpl_proc_null) meas%mpl_comm_size = meas%mpl_comm_size   + real(prc%length(5)*R8B)
   if(prc%nghb(6) /= mpl_proc_null) meas%mpl_comm_size = meas%mpl_comm_size   + real(prc%length(6)*R8B)


#else /* if no USE_ADCL !! */


!------------------------------------------------------------------------------
! MPI Isend Irecv, Sendrecv and Persistent,... 
!-----------------------------------------------------------------------------


   do neighbor = 1, prc%nnghb

      ! Only start buffer copy and communication if neigh is existing
      if (prc%nghb(neighbor) /= mpl_proc_null .and. prc%nghb(neighbor) .ge. 0) then 

         call mpl_fill_buffer(neighbor,lb_dom,fIn,prc,meas)

         call mpl_communicate_buffer(neighbor, &
#ifdef MPI_SUBARRAY
                  lb_dom,fIn, &
#endif
                  prc,meas)

#ifndef PERSISTENT2
         ! Either, the buffer is copied back after each neighbor exchange
         call mpl_read_buffer(neighbor,lb_dom,fIn,prc,meas)
      else
         ! If no neighbor is existing in this direction, reset the handle
         ! Do NOT reset the handle for the persistent communication
         prc%ireq(neighbor)           = MPI_REQUEST_NULL 
         prc%ireq(neighbor+prc%nnghb) = MPI_REQUEST_NULL 
#endif
      endif 

   enddo

#if defined(PERSISTENT2)  || defined(BUFFER_EACH)
#if defined(PERSISTENT2)
      call mpi_startall(2*prc%nnghb_avail, prc%ireq, prc%ierr)
      call mpi_waitall( 2*prc%nnghb_avail, prc%ireq, prc%status, prc%ierr) 
#else
      call mpi_waitall(2*prc%nnghb,prc%ireq,prc%status,prc%ierr) 
#endif
   do neighbor = 1, prc%nnghb
      if (prc%nghb(neighbor) /= mpl_proc_null .and. prc%nghb(neighbor) .ge. 0) then 
         call mpl_read_buffer(neighbor,lb_dom,fIn,prc,meas)
      endif 
   enddo
#endif

#endif /* USE_ADCL */
#endif /* USE_MPI */



   call cpu_time_measure(meas%tEnd_comm)
   meas%comm_duration = meas%tEnd_comm - meas%tSt_comm + meas%comm_duration




#ifdef DEBUG_MPI
   call mpi_barrier(mpi_comm_world,prc%ierr)
stop
#endif


  end subroutine mpl_exchange
!------------------------------------------------------------------------




#endif /* USE_CAF */



  subroutine mpl_fill_buffer(neighbor,lb_dom,fIn,prc,meas)
  !------------------------------------------------------------------------
  !
  !

      type(lb_block),intent(inout) :: lb_dom
      type(mpl_var)  :: prc
      type(measure)  :: meas
      integer,intent(in)  :: neighbor
      integer        :: ii,jj,kk,l
      real(R8B)      :: copy_start,copy_end
      real(R8B) :: fIn(NDX(nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))

#ifndef MPI_SUBARRAY
#ifdef DEBUG_MPI
      write(50+prc%rk,*) "Neighbor index is",neighbor
      write(50+prc%rk,*) "sending to",prc%nghb(neighbor)
      write(50+prc%rk,*) "Limits are" 
      write(50+prc%rk,*) prc%pos_sends(neighbor,1),prc%pos_send(neighbor,1)
      write(50+prc%rk,*) prc%pos_sends(neighbor,2),prc%pos_send(neighbor,2)
      write(50+prc%rk,*) prc%pos_sends(neighbor,3),prc%pos_send(neighbor,3)
#endif


      call cpu_time_measure(copy_start)
      do kk = 1, prc%pos_send(neighbor,3)-prc%pos_sends(neighbor,3)+1
         do jj = 1, prc%pos_send(neighbor,2)-prc%pos_sends(neighbor,2)+1
            do ii = 1, prc%pos_send(neighbor,1)-prc%pos_sends(neighbor,1)+1
               do l = 1,prc%ndir_pdf(neighbor)
                  mpl_buffer(neighbor)%send(l,ii,jj,kk) =    & 

#ifndef DEBUG_MPI

fIn(NDX(prc%dsend_pdf(neighbor,l),ii+prc%pos_sends(neighbor,1)-1, jj+prc%pos_sends(neighbor,2)-1, kk+prc%pos_sends(neighbor,3)-1))


#else
      (prc%rk+1)*100000+prc%dsend_pdf(neighbor,l)+10000*(ii+prc%pos_sends(neighbor,1)-1)+&
      1000*(jj+prc%pos_sends(neighbor,2)-1)+100*(kk+prc%pos_sends(neighbor,3)-1)
      write(50+prc%rk,*) prc%dsend_pdf(neighbor,l),ii+prc%pos_sends(neighbor,1)-1,&
       jj+prc%pos_sends(neighbor,2)-1, kk+prc%pos_sends(neighbor,3)-1,mpl_buffer(neighbor)%send(l,ii,jj,kk)
#endif
               enddo
            enddo
         enddo
      enddo

      call cpu_time_measure(copy_end)
      meas%mpl_cp_buf = meas%mpl_cp_buf + copy_end - copy_start
#endif /* MPI_SUBARRAY */

  end subroutine mpl_fill_buffer




  subroutine mpl_read_buffer(neighbor,lb_dom,fIn,prc,meas)
  !------------------------------------------------------------------------
  !
  !   Fill up the communication buffer
  !

      type(lb_block),intent(inout) :: lb_dom
      type(mpl_var)  :: prc
      type(measure)  :: meas
      integer,intent(in)  :: neighbor
      integer        :: ii,jj,kk,l
      real(R8B)      :: copy_start,copy_end
      real(R8B) :: fIn(NDX(nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))
#ifndef MPI_SUBARRAY

#ifdef DEBUG_MPI
      write(41+prc%rk,*) "recving from",prc%nghb(neighbor)
      write(41+prc%rk,*) "Limits are" 
      write(41+prc%rk,*) prc%pos_recvs(neighbor,1),prc%pos_recv(neighbor,1)
      write(41+prc%rk,*) prc%pos_recvs(neighbor,2),prc%pos_recv(neighbor,2)
      write(41+prc%rk,*) prc%pos_recvs(neighbor,3),prc%pos_recv(neighbor,3)
#endif /* DEBUG_MPI */


      call cpu_time_measure(copy_start)

      do kk = 1, prc%pos_recv(neighbor,3)-prc%pos_recvs(neighbor,3)+1
         do jj = 1, prc%pos_recv(neighbor,2)-prc%pos_recvs(neighbor,2)+1
            do ii = 1, prc%pos_recv(neighbor,1)-prc%pos_recvs(neighbor,1)+1
               do l = 1,prc%ndir_pdf(neighbor)

fIn(NDX(prc%drecv_pdf(neighbor,l),ii+prc%pos_recvs(neighbor,1)-1,jj+prc%pos_recvs(neighbor,2)-1, kk+prc%pos_recvs(neighbor,3)-1)) & 
            =  mpl_buffer(neighbor)%recv(l,ii,jj,kk)  

#ifdef DEBUG_MPI
      write(41+prc%rk,*) prc%drecv_pdf(neighbor,l),ii+prc%pos_recvs(neighbor,1)-1,&
      jj+prc%pos_recvs(neighbor,2)-1, kk+prc%pos_recvs(neighbor,3)-1,mpl_buffer(neighbor)%recv(l,ii,jj,kk)
#endif
               enddo
            enddo
         enddo
      enddo



      call cpu_time_measure(copy_end)
      meas%mpl_cp_buf_back = meas%mpl_cp_buf_back + copy_end - copy_start
#endif /* MPI_SUBARRAY */

  end subroutine mpl_read_buffer



  subroutine mpl_communicate_buffer(neighbor, & 
#ifdef MPI_SUBARRAY
lb_dom,fIn, &
#endif
prc,meas)
  !------------------------------------------------------------------------
  !
  !   Exchange the communication buffers
  !   Implementation of different MPI communications
  ! 
  !
      type(mpl_var)  :: prc
      type(measure)  :: meas
      integer,intent(in)  :: neighbor
      real(R8B)      :: comm_start,comm_end
#ifdef MPI_SUBARRAY
      type(lb_block),intent(inout) :: lb_dom
      real(R8B) :: fIn(NDX(nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))
#endif


      call cpu_time_measure(comm_start)

#if defined SENDRECV
      call mpi_sendrecv(  mpl_buffer(neighbor)%send,                                     &
                          prc%length(neighbor),mpi_double_precision,prc%nghb(neighbor), 1,      &
                          mpl_buffer(neighbor)%recv,                                     & 
                          prc%length(neighbor),mpi_double_precision,prc%nghb(neighbor),1,       &
                          prc%cart_comm,prc%stat,prc%ierr)

#elif defined ISEND_IRECV
      call mpi_irecv(mpl_buffer(neighbor)%recv,                &
                     prc%length(neighbor),                     & 
                     mpi_double_precision,              &
                     prc%nghb(neighbor),mpi_any_tag,           &
                     prc%cart_comm,prc%ireq(neighbor),prc%ierr)

      call mpi_isend(mpl_buffer(neighbor)%send,                &
                     prc%length(neighbor),                     & 
                     mpi_double_precision,              &
                     prc%nghb(neighbor),  1,                   &
                     prc%cart_comm,prc%ireq(prc%nnghb+neighbor),prc%ierr)

#ifndef BUFFER_EACH
      call mpi_wait(prc%ireq(neighbor  ),prc%status(:,neighbor  ),prc%ierr) 
      call mpi_wait(prc%ireq(prc%nnghb+neighbor  ),prc%status(:,prc%nnghb+neighbor),prc%ierr)
#endif /* BUFFER_EACH */


!      call mpi_waitall(2*prc%nnghb,prc%ireq,prc%status,prc%ierr) 
#elif defined PERSISTENT 

      ! initialize communications

      call mpi_send_init(mpl_buffer(neighbor)%send,                &
                     prc%length(neighbor),                         & 
                     mpi_double_precision,                         &
                     prc%nghb(neighbor),  1,                       &
                     prc%cart_comm,prc%ireq(prc%nnghb+neighbor),prc%ierr)


      call mpi_recv_init(mpl_buffer(neighbor)%recv,                &
                     prc%length(neighbor),                         & 
                     mpi_double_precision,                         &
                     prc%nghb(neighbor),mpi_any_tag,               &
                     prc%cart_comm,prc%ireq(neighbor),prc%ierr)

      call mpi_start(prc%ireq(prc%nnghb+neighbor), prc%ierr)
      call mpi_start(prc%ireq(neighbor), prc%ierr)
      call mpi_wait(prc%ireq(neighbor  ),prc%status(:,neighbor  ),prc%ierr) 
      call mpi_wait(prc%ireq(prc%nnghb+neighbor  ),prc%status(:,prc%nnghb+neighbor),prc%ierr)
      call mpi_request_free(prc%ireq(neighbor  ),prc%ierr) 
      call mpi_request_free(prc%ireq(prc%nnghb+neighbor  ),prc%ierr)

!      call mpi_startall(count, array_of_requests, prc%ierr)
!      call mpi_waitall()

#elif defined PERSISTENT2 


!      call mpi_startall(2*prc%nnghb, prc%ireq, prc%ierr)
!      call mpi_start(prc%ireq(neighbor),           prc%ierr)
!      call mpi_start(prc%ireq(prc%nnghb+neighbor), prc%ierr)


#elif defined MPI_SUBARRAY
#ifndef PERSISTENT2
      call mpi_irecv(   fIn,       & 
                        1,                & 
                        recv_dat(neighbor),      & 
                        prc%nghb(neighbor),      & 
                        mpi_any_tag,      & 
                        prc%cart_comm,    & 
                        prc%ireq(neighbor),         & 
                        prc%ierr)
write(*,*) prc%rk,"irecv error    ",prc%ierr 
      if(prc%nghb(neighbor) /= mpl_proc_null) meas%mpl_comm_size   = meas%mpl_comm_size   + real(prc%length(neighbor)*R8B)

      call mpi_isend(  fIn,                                                 &
                          1,send_dat(neighbor),prc%nghb(neighbor),1,                   &
                          prc%cart_comm,prc%ireq(prc%nnghb + neighbor),prc%ierr)
write(*,*) prc%rk,"isend error    ",prc%ierr 

write(*,*) prc%rk,"request handels",prc%ireq 
write(*,*) prc%rk,"waiting"
      call mpi_wait(prc%ireq(neighbor  ),prc%status(:,neighbor  ),prc%ierr) 
write(*,*) prc%rk,"waiting"
      call mpi_wait(prc%ireq(prc%nnghb+neighbor),prc%status(:,prc%nnghb+neighbor),prc%ierr)
#endif
#endif








      call cpu_time_measure(comm_end)
      meas%mpl_exch      = meas%mpl_exch + comm_end - comm_start
      meas%mpl_comm_size = meas%mpl_comm_size   + real(prc%length(neighbor)*R8B)


  end subroutine mpl_communicate_buffer









#ifdef INIT_WITH_ROOT
  subroutine mpl_decompose_domain(lb_dom,s_par,prc)
  !------------------------------------------------------------------------
  !
  !  The domain is decomposed and assigned to the threads.
  !  A loop goes over the subdomains step by step.
  !  If the current subdomain is the domain of the root thread, then the arrays are simply copied
  !  For all other domains, the root thread sends the part of the array to the other threads
  !
  !
   implicit none
   type(lb_block),intent(inout)     :: lb_dom
   type(sim_parameter),intent(in)   :: s_par 
   type(mpl_var)   :: prc
   integer                :: crd_array(3)
   integer                :: i,j,k,length,lx,ly,lz
   integer                :: b_l(3),b_u(3),ul(3)
   integer :: ab(prc%np(1),prc%np(2),prc%np(3),3),ae(prc%np(1),prc%np(2),prc%np(3),3)
#ifdef USE_CAF
    real(R8B),allocatable       :: fIn_buf(:,:,:,:)[:,:,:],rho_buf(:,:,:)[:,:,:]
    real(R8B),allocatable       :: u_buf(:,:,:,:)[:,:,:]
    real(R8B),allocatable       :: omega_buf(:,:,:)[:,:,:]
    integer,allocatable    :: state_buf(:,:,:)[:,:,:]
 




    allocate(fIn_buf(NDX(nnod,cobnd(1),cobnd(2),cobnd(3)))[prc%np(1),prc%np(2),*])
    allocate(rho_buf(cobnd(1),cobnd(2),cobnd(3))[prc%np(1),prc%np(2),*])
    allocate(u_buf(NDX(NDIM,cobnd(1),cobnd(2),cobnd(3)))[prc%np(1),prc%np(2),*])
    allocate(state_buf(cobnd(1),cobnd(2),cobnd(3))[prc%np(1),prc%np(2),*])
    allocate(omega_buf(cobnd(1),cobnd(2),cobnd(3))[prc%np(1),prc%np(2),*])
#endif 

   lx =  lb_dom%lx(1)
   ly =  lb_dom%lx(2)
   lz =  lb_dom%lx(3)

   ab = prc%bnd(:,:,:,:,1)
   ae = prc%bnd(:,:,:,:,2)

    if(prc%rk == prc%root_th) then
       ! master thread distributes the subdomains
       ! go over all subdomains, one by one
       
       do k=1,prc%np(3)
          do j=1,prc%np(2)
             do i=1,prc%np(1)
                b_l(1:3) = prc%bnd(i,j,k,1:3,1) 
                b_u(1:3) = prc%bnd(i,j,k,1:3,2) 
                length = (b_u(1)-b_l(1) + 1) * (b_u(2)-b_l(2)+1) * (b_u(3)-b_l(3)+1)
                if(i == (prc%root_crd(1)+1) .and. & 
                   j == (prc%root_crd(2)+1) .and. &
                   k == (prc%root_crd(3)+1)) then 

                   lb_dom%rho(1:lx,1:ly,1:lz)&
                     & = lb_dom%grho(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                     &        prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                     &        prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2))

                   lb_dom%u(NDX(:,1:lx,1:ly,1:lz))&
                     & = lb_dom%gu( & 
NDX(:,ab(i,j,k,1):ae(i,j,k,1),ab(i,j,k,2):ae(i,j,k,2),ab(i,j,k,3):ae(i,j,k,3)))
!&NDX(1:NDIM,prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2),prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2),prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)))

                  lb_dom%fIn(NDX(:,1:lx,1:ly,1:lz))&
                     & = lb_dom%gfIn( & 
NDX(:,ab(i,j,k,1):ae(i,j,k,1),ab(i,j,k,2):ae(i,j,k,2),ab(i,j,k,3):ae(i,j,k,3)))
!&NDX(1:nnod,prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2),prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)))

                   lb_dom%state(1:lx,1:ly,1:lz)&
                     & = lb_dom%gstate(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                     &        prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                     &        prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2))

#ifdef SPONGE
                   lb_dom%omega(1:lx,1:ly,1:lz)&
                     & = lb_dom%gomega(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                     &        prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                     &        prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2))
#endif
                   lb_dom%x(1:lx) = s_par%g_x(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2)) 
                   lb_dom%y(1:ly) = s_par%g_y(prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2)) 
                   lb_dom%z(1:lz) = s_par%g_z(prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)) 

                else  
#ifdef USE_MPI
                crd_array=(/(i-1),(j-1),(k-1)/)
                call MPI_CART_RANK(prc%cart_comm,crd_array,target_rank,prc%ierr)
                call MPI_SEND(lb_dom%grho(b_l(1):b_u(1),b_l(2):b_u(2),b_l(3):b_u(3)),length, & 
                           &  mpi_double_precision,target_rank,tag_rho,prc%cart_comm,prc%ierr)
                call MPI_SEND(lb_dom%gu(NDX(1:NDIM,b_l(1):b_u(1),b_l(2):b_u(2),b_l(3):b_u(3))),length*NDIM,&
                          &  mpi_double_precision,target_rank,tag_rho,prc%cart_comm,prc%ierr)
                call MPI_SEND(lb_dom%gfIn(NDX(1:nnod,b_l(1):b_u(1),b_l(2):b_u(2),b_l(3):b_u(3))),length*nnod,&
                           &  mpi_double_precision,target_rank,tag_rho,prc%cart_comm,prc%ierr)
                call MPI_SEND(lb_dom%gstate(b_l(1):b_u(1),b_l(2):b_u(2),b_l(3):b_u(3)),length,MPI_INTEGER,&
                           &  target_rank,tag_rho,prc%cart_comm,prc%ierr)
#ifdef SPONGE
                call MPI_SEND(lb_dom%gomega(b_l(1):b_u(1),b_l(2):b_u(2),b_l(3):b_u(3)),length,MPI_double_precision,&
                           &  target_rank,tag_rho,prc%cart_comm,prc%ierr)
#endif
#endif /*USE_MPI*/
#ifdef USE_CAF
                   rho_buf(1:lx,1:ly,1:lz)[i,j,k]&
                     & = lb_dom%grho(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                     &        prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                     &        prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2))

                   fIn_buf(NDX(1:nnod,1:lx,1:ly,1:lz))[i,j,k]&
                     & = lb_dom%gfIn( & 
NDX(:,ab(i,j,k,1):ae(i,j,k,1),ab(i,j,k,2):ae(i,j,k,2),ab(i,j,k,3):ae(i,j,k,3)))
!& NDX(1:nnod,prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2),prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)))

                   u_buf(NDX(1:NDIM,1:lx,1:ly,1:lz))[i,j,k]&
                     & = lb_dom%gu( & 
NDX(:,ab(i,j,k,1):ae(i,j,k,1),ab(i,j,k,2):ae(i,j,k,2),ab(i,j,k,3):ae(i,j,k,3)))
!& NDX(1:NDIM,prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2),prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)))

                   state_buf(1:lx,1:ly,1:lz)[i,j,k]&
                     & = lb_dom%gstate(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                     &        prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                     &        prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2))
#ifdef SPONGE
                   omega_buf(1:lx,1:ly,1:lz)[i,j,k]&
                     & = lb_dom%gomega(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                     &        prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                     &        prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2))
#endif /* SPONGE */

#endif /*USE_CAF*/
                end if
             end do  
          end do  
       end do  
    else ! if current thread is not root thread 
#ifdef USE_MPI
       ! all other threads: 
       ! wait until mpi_recv is done
       b_l(1:3) = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1:3,1) 
       b_u(1:3) = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1:3,2) 
       length = (b_u(1)-b_l(1) + 1) * (b_u(2)-b_l(2)+1) * (b_u(3)-b_l(3)+1)
       ul(1:3) = b_u - b_l + 1
       ul(1:3) = (/lx,ly,lz/)
       
       call MPI_RECV(       lb_dom%rho(1:ul(1),1:ul(2),1:ul(3)),length,mpi_double_precision,&
                        &  prc%root_th,tag_rho,prc%cart_comm,prc%stat,prc%ierr)
       call MPI_RECV(  lb_dom%u(NDX(1:NDIM,1:ul(1),1:ul(2),1:ul(3))),length*NDIM,mpi_double_precision,&
                        &  prc%root_th,tag_rho,prc%cart_comm,prc%stat,prc%ierr)
       call MPI_RECV(lb_dom%fIn(NDX(1:nnod,1:ul(1),1:ul(2),1:ul(3))),length*nnod,mpi_double_precision,&
                        &  prc%root_th,tag_rho,prc%cart_comm,prc%stat,prc%ierr)
       call MPI_RECV(     lb_dom%state(1:ul(1),1:ul(2),1:ul(3)),length,MPI_INTEGER,&
                        &  prc%root_th,tag_rho,prc%cart_comm,prc%stat,prc%ierr)
#ifdef SPONGE
       call MPI_RECV(     lb_dom%omega(1:ul(1),1:ul(2),1:ul(3)),length,mpi_double_precision,&
                        &  prc%root_th,tag_rho,prc%cart_comm,prc%stat,prc%ierr)
#endif
#endif /* USE_MPI */
   end if
#ifdef USE_CAF
   call sync_images(g_team)
   if (prc%rk /= 0 ) then ! copy from coarrays to local arrays
      lb_dom%rho(1:lx,1:ly,1:lz)        = rho_buf(1:lx,1:ly,1:lz)
      lb_dom%u(NDX(1:NDIM,1:lx,1:ly,1:lz))   = u_buf(NDX(1:NDIM,1:lx,1:ly,1:lz))
      lb_dom%fIn(NDX(1:nnod,1:lx,1:ly,1:lz)) = fIn_buf(NDX(1:nnod,1:lx,1:ly,1:lz))
      lb_dom%state(1:lx,1:ly,1:lz)      = state_buf(1:lx,1:ly,1:lz)
#ifdef SPONGE
      lb_dom%omega(1:lx,1:ly,1:lz)      = omega_buf(1:lx,1:ly,1:lz)
#endif
   end if
   ! deallocate buffers
   deallocate(rho_buf)
   deallocate(fIn_buf)
   deallocate(u_buf)
   deallocate(state_buf)
#ifdef SPONGE
   deallocate(omega_buf)
#endif /* SPONGE  */

#endif /* USE_CAF */
 end subroutine mpl_decompose_domain
!------------------------------------------------------------------------
#endif /* INIT_WITH_ROOT */




  subroutine mpl_barrier()
  !------------------------------------------------------------------------
  !
  !  The domain is decomposed and assigned to the threads.
  !  A loop goes over the subdomains step by step.
  !  If the current subdomain is the domain of the root thread, then the arrays are simply copied
  !  For all other domains, the root thread sends the part of the array to the other threads
  !
  !
   implicit none
   integer         :: ierr

#ifdef USE_MPI
   call mpi_barrier(mpi_comm_world,ierr)
#endif
#ifdef USE_CAF
   call sync_all()
#endif

   end subroutine mpl_barrier







#ifdef USE_OLD_MPI

subroutine mpl_set_nodes_to_send(lb_dom,prc)
!------------------------------------------------------------------------
!
! OLD ROUTINE!! 
!
!#!define USE_NEW
   type(mpl_var)                    :: prc
   type(lb_block)                   :: lb_dom
   integer                     :: lx,ly,lz,cnt,ind,l,ll
   integer                     :: xs(3),ys(3),zs(3) 
   integer                     :: xr(3),yr(3),zr(3) 
   integer                     :: zeropos,start_it,end_it,relative,ind_rel,ghostn
#ifdef USE_NEW 
   integer                     :: sig(3),dimi  
#endif
      ! fill up the sending/receiving buffer sizes (from pos_sends to pos_send, same for recv)
      ! set all to zero first (there are more buffer entries than needed for reduced comm)

      lx =  lb_dom%lx(1)
      ly =  lb_dom%lx(2)
      lz =  lb_dom%lx(3)

      ! First count, how many densities have to be exchanged in one direction. Choose x-direction
      ! This could only be a problem, if the stencil is not symmetric!! this should always be the case.
#ifdef COMM_REDUCED
      prc%ndir_pdf = 0
      do l = 1,nnod
         if(cx(l,1) == 1) prc%ndir_pdf=prc%ndir_pdf+1  
      enddo
#else
      prc%ndir_pdf = nnod
#endif

      allocate(prc%dsend_pdf(prc%nnghb,prc%ndir_pdf))
      allocate(prc%drecv_pdf(prc%nnghb,prc%ndir_pdf))


      ! reduced mpi communication
      ! set all pdf directions to zero
      prc%drecv_pdf = 0
      prc%dsend_pdf = 0

      prc%pos_sends = 0
      prc%pos_recvs = 0 
      prc%pos_send  = 0
      prc%pos_recv  = 0

!#ifdef BUFFER_EACH
   ! fill up the sending/receiving buffer sizes (from pos_sends to pos_send, same for recv)
! Identify zero
zeropos=-1
relative=0
do ind=1,nnod
   if(cx(ind,1)==0 .and. cx(ind,2) == 0 .and. cx(ind,3)==0) zeropos=ind
enddo

if(zeropos==1) then
start_it=2
relative=-1
end_it=prc%nnghb+1
elseif(zeropos==nnod) then
start_it=1
end_it=prc%nnghb
endif



do ind=start_it,end_it
   ind_rel = ind + relative
#ifdef BUFFER_EACH
   ! when each buffer is exchanged separately, ghost nodes don't have to be communicated
   ghostn = 0
#else /* ORTHOGONAL DIRECTIONS */
   ! orthogonal direction exchange needs ghost node values for writing to corners 
   ghostn = 1
#endif
   ! treat x direction
   if(cx(ind,1) ==0) then
      xs(1) = 1  - ghostn
      xs(2) = lx + ghostn
      xr(1) = 1  - ghostn
      xr(2) = lx + ghostn
   elseif(cx(ind,1)==1) then
      xs(1:2) = lx
      xr(1:2) = lx+1
   else ! -1
      xs(1:2) = 1
      xr(1:2) = 0
   endif

   ! treat y direction
   if(cx(ind,2) ==0) then
      ys(1) = 1  - ghostn
      ys(2) = ly + ghostn
      yr(1) = 1  - ghostn
      yr(2) = ly + ghostn
   elseif(cx(ind,2)==1) then
      ys(1:2) = ly
      yr(1:2) = ly+1
   else ! -1
      ys(1:2) = 1
      yr(1:2) = 0
   endif

   ! treat z direction
   if(cx(ind,3) ==0) then
      zs(1) = 1  - ghostn
      zs(2) = lz + ghostn
      zr(1) = 1  - ghostn
      zr(2) = lz + ghostn
   elseif(cx(ind,3)==1) then
      zs(1:2) = lz
      zr(1:2) = lz+1
   else ! -1
      zs(1:2) = 1
      zr(1:2) = 0
   endif




   prc%pos_sends(ind_rel,1:3) = (/xs(1),  ys(1), zs(1)/)
   prc%pos_send (ind_rel,1:3) = (/xs(2),  ys(2), zs(2)/)
   prc%pos_recvs(ind_rel,1:3) = (/xr(1),  yr(1), zr(1)/)
   prc%pos_recv (ind_rel,1:3) = (/xr(2),  yr(2), zr(2)/)

   do ll=1,nnod
   prc%dsend_pdf(ind_rel,l)   = ll
   prc%drecv_pdf(ind_rel,l)   = ll
   enddo

#ifdef DEBUG_MPI
   write(66+prc%rk,*) "Index",ind_rel,"  dir",cx(ind,:)
   write(66+prc%rk,*) prc%pos_sends(ind_rel,:) 
   write(66+prc%rk,*) prc%pos_send(ind_rel,:) 
   write(66+prc%rk,*) prc%pos_recvs(ind_rel,:) 
   write(66+prc%rk,*) prc%pos_recv(ind_rel,:) 
#endif DEBUG_MPI
enddo







!#else /* BUFFER_EACH */
#ifdef USE_NEW /* NOT WORKING YET!!!*/ 
      ! this will fill up the density locations 1 to 6
      do ind  = 1,5,2

      if(ind == 1) then  ! index 1 and to
      ! first direction : x-
      prc%pos_sends(ind,1:3) = (/1,    0,    0    /) 
      prc%pos_send (ind,1:3) = (/1,    ly+1, lz+1 /) 
      prc%pos_recvs(ind,1:3) = (/0,    0,    0    /) 
      prc%pos_recv (ind,1:3) = (/0,    ly+1, lz+1 /) 
      ! first direction : x+
      prc%pos_sends(ind+1,1:3) = (/lx, 0,    0    /) 
      prc%pos_send (ind+1,1:3) = (/lx, ly+1, lz+1 /) 
      prc%pos_recvs(ind+1,1:3) = (/lx+1, 0,    0    /) 
      prc%pos_recv (ind+1,1:3) = (/lx+1, ly+1, lz+1 /) 
      dimi = 1 
      sig  = (/-1,0,0 /)

      elseif( ind==3) then
      ! first direction : y-
      prc%pos_sends(ind,1:3) = (/0,    1,    0    /) 
      prc%pos_send (ind,1:3) = (/lx+1, 1,    lz+1 /) 
      prc%pos_recvs(ind,1:3) = (/0,    0,    0    /) 
      prc%pos_recv (ind,1:3) = (/lx+1, 0,    lz+1 /) 

      prc%pos_sends(ind+1,1:3) = (/0,    ly,       0    /) 
      prc%pos_send (ind+1,1:3) = (/lx+1, ly,       lz+1 /) 
      prc%pos_recvs(ind+1,1:3) = (/0,    ly+1, 0    /) 
      prc%pos_recv (ind+1,1:3) = (/lx+1, ly+1, lz+1 /) 
      dimi = 2 
      sig  = (/0,-1,0 /)

      elseif( ind==5) then
      ! first direction : z-
      prc%pos_sends(ind,1:3) = (/0,    0,    1    /) 
      prc%pos_send (ind,1:3) = (/lx+1, ly+1, 1    /) 
      prc%pos_recvs(ind,1:3) = (/0,    0,    0    /) 
      prc%pos_recv (ind,1:3) = (/lx+1, ly+1, 0    /) 
      prc%pos_sends(ind+1,1:3) = (/0,    0,    lz       /)
      prc%pos_send (ind+1,1:3) = (/lx+1, ly+1, lz       /)
      prc%pos_recvs(ind+1,1:3) = (/0,    0,    lz+1 /) 
      prc%pos_recv (ind+1,1:3) = (/lx+1, ly+1, lz+1 /) 
      dimi = 3 
      sig  = (/0,0,-1 /)
      endif

      cnt = 0
      do l = 1,nnod
#ifdef COMM_REDUCED
         if(cx(l,dimi) == -1) then
         cnt = cnt+1
         prc%dsend_pdf(ind,cnt)   = l 
         prc%drecv_pdf(ind+1,cnt) = l 
         do ll=1,nnod
            if(cx(ll,1) == sig(1)*cx(l,1) .and. &
               cx(ll,2) == sig(2)*cx(l,2) .and. &
               cx(ll,3) == sig(3)*cx(l,3)) then
               prc%drecv_pdf(ind,cnt)   = ll
               prc%dsend_pdf(ind+1,cnt) = ll
            endif
         enddo
         endif 
#else
         ! if all densities should be communicated, simply fill up array
         prc%dsend_pdf(ind,l)   = l
         prc%drecv_pdf(ind,l)   = l
         prc%dsend_pdf(ind+1,l) = l
         prc%drecv_pdf(ind+1,l) = l
#endif
      enddo
      enddo

#else /* USE OLD ?*/

#ifdef ORIGINAL

!-------------------------
! old implementation

      ind = 1
      ! first direction : x-
      prc%pos_sends(ind,1:3) = (/1,    0,    0    /) 
      prc%pos_send (ind,1:3) = (/1,    ly+1, lz+1 /) 
      prc%pos_recvs(ind,1:3) = (/0,    0,    0    /) 
      prc%pos_recv (ind,1:3) = (/0,    ly+1, lz+1 /) 
      ! first direction : x+
      prc%pos_sends(ind+1,1:3) = (/lx,       0,    0    /) 
      prc%pos_send (ind+1,1:3) = (/lx,       ly+1, lz+1 /) 
      prc%pos_recvs(ind+1,1:3) = (/lx+1, 0,    0    /) 
      prc%pos_recv (ind+1,1:3) = (/lx+1, ly+1, lz+1 /) 

      cnt = 0
      do l = 1,nnod
#ifdef COMM_REDUCED
         if(cx(l,1) == -1) then
         cnt = cnt+1
         prc%dsend_pdf(ind,cnt) = l 
         prc%drecv_pdf(ind+1,cnt) = l 
         do ll=1,nnod
            if(cx(ll,1) == -cx(l,1) .and. cx(ll,2) == cx(l,2) .and. cx(ll,3) == cx(l,3)) then
               prc%drecv_pdf(ind,cnt)   = ll
               prc%dsend_pdf(ind+1,cnt) = ll
            endif
         enddo
         endif 
#else
         prc%dsend_pdf(ind,l) = l
         prc%drecv_pdf(ind,l) = l
         prc%dsend_pdf(ind+1,l) = l
         prc%drecv_pdf(ind+1,l) = l
#endif
      enddo


      ind = 3
      ! first direction : y-
      prc%pos_sends(ind,1:3) = (/0,    1,    0    /) 
      prc%pos_send (ind,1:3) = (/lx+1, 1,    lz+1 /) 
      prc%pos_recvs(ind,1:3) = (/0,    0,    0    /) 
      prc%pos_recv (ind,1:3) = (/lx+1, 0,    lz+1 /) 

      prc%pos_sends(ind+1,1:3) = (/0,    ly,       0    /) 
      prc%pos_send (ind+1,1:3) = (/lx+1, ly,       lz+1 /) 
      prc%pos_recvs(ind+1,1:3) = (/0,    ly+1, 0    /) 
      prc%pos_recv (ind+1,1:3) = (/lx+1, ly+1, lz+1 /) 

      cnt = 0
      do l = 1,nnod
#ifdef COMM_REDUCED
         if(cx(l,2) == -1) then
         cnt = cnt+1
         prc%dsend_pdf(ind,cnt) = l 
         prc%drecv_pdf(ind+1,cnt) = l 
         do ll=1,nnod
            if(cx(ll,1) == cx(l,1) .and. cx(ll,2) == -cx(l,2) .and. cx(ll,3) == cx(l,3)) then
               prc%drecv_pdf(ind,cnt) = ll
               prc%dsend_pdf(ind+1,cnt) = ll
            endif
         enddo
         endif 
#else
         prc%dsend_pdf(ind,l) = l
         prc%drecv_pdf(ind,l) = l
         prc%dsend_pdf(ind+1,l) = l
         prc%drecv_pdf(ind+1,l) = l
#endif
      enddo



      ind = 5
      ! first direction : z-
      prc%pos_sends(ind,1:3) = (/0,    0,    1    /) 
      prc%pos_send (ind,1:3) = (/lx+1, ly+1, 1    /) 
      prc%pos_recvs(ind,1:3) = (/0,    0,    0    /) 
      prc%pos_recv (ind,1:3) = (/lx+1, ly+1, 0    /) 
      prc%pos_sends(ind+1,1:3) = (/0,    0,    lz       /)
      prc%pos_send (ind+1,1:3) = (/lx+1, ly+1, lz       /)
      prc%pos_recvs(ind+1,1:3) = (/0,    0,    lz+1 /) 
      prc%pos_recv (ind+1,1:3) = (/lx+1, ly+1, lz+1 /) 
      cnt = 0
      do l = 1,nnod
#ifdef COMM_REDUCED
         if(cx(l,3) == -1) then
         cnt = cnt+1
         prc%dsend_pdf(ind,cnt) = l 
         prc%drecv_pdf(ind+1,cnt) = l 
         do ll=1,nnod
            if(cx(ll,1) == cx(l,1) .and. cx(ll,2) == cx(l,2) .and. cx(ll,3) == -cx(l,3)) then
            prc%drecv_pdf(ind,cnt) = ll
            prc%dsend_pdf(ind+1,cnt) = ll
            endif
         enddo
         endif 
#else
         prc%dsend_pdf(ind,l) = l
         prc%drecv_pdf(ind,l) = l
         prc%dsend_pdf(ind+1,l) = l
         prc%drecv_pdf(ind+1,l) = l
#endif
      enddo

#endif /* ORIGINAL */

#endif /* OLD */
!#endif /* BUFFER_EACH */ 




#ifdef D2Q9
      prc%pos_sends(:,3) = 1
      prc%pos_send(:,3)  = 1
      prc%pos_recvs(:,3) = 1
      prc%pos_recv(:,3)  = 1
#endif


end subroutine mpl_set_nodes_to_send
!------------------------------------------------------------------------

#endif /* USE_OLD_MPI */


end module mpl_set


