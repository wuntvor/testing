module function_set
   use lbmodel
   use tools
   use mpl_set
   implicit none
#include "include/replace.h"
contains

#ifdef DEBUG
subroutine dump_matrix(matrix,lb_dom,prc,s_par,description,tstep)
!------------------------------------------------------------------------
!
! Master thread gathers the output from all threads in order to write it out.
!


   type(lb_block),intent(in)       :: lb_dom
   real(R8B) :: matrix(LB_NODE(nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))
   type(sim_parameter),intent(in)  :: s_par
   type(mpl_var)      :: prc
   integer,intent(in) :: tstep
   integer            :: i,j,k,l,lx,ly,lz
   character(*)       :: description
 
   lx = lb_dom%lx(1)
   ly = lb_dom%lx(2)
   lz = lb_dom%lx(3)

   open(71,file='dump.mtr',position='append')
   write(71,'(a9,i4)') "Timestep ",tstep 
   write(71,*) description 
   write(71,*) " " 
   do k=0,lb_dom%lx(3)
   do j=0,lb_dom%lx(2)
   do i=0,lb_dom%lx(1)
   do l=1,nnod
      write(71,'(5i3,f11.6)') tstep,i,j,k,l,matrix(LB_NODE(l,i,j,k))
   enddo
   enddo
   enddo
   enddo

   write(71,*) " " 
   write(71,*) "-----------------------------------------------------------------------------" 
   write(71,*) " "
   close(71) 

end subroutine dump_matrix 
!------------------------------------------------------------------------
#endif




subroutine gather_output(lb_dom,countOut,tStep,s_par,prc)
!------------------------------------------------------------------------
!
! Master thread gathers the output from all threads in order to write it out.
!


   type(lb_block),intent(inout)    :: lb_dom
   type(sim_parameter),intent(inout)  :: s_par
   type(mpl_var)      :: prc
   integer :: countOut,tStep
   integer :: lx,ly,lz
   integer :: b_l(3),b_u(3)
#ifdef INIT_WITH_ROOT
   integer :: length
   integer :: i,j,k
#endif
   integer :: ab(prc%np(1),prc%np(2),prc%np(3),3),ae(prc%np(1),prc%np(2),prc%np(3),3)
#ifdef USE_CAF
   real(R8B),allocatable       :: u_buf(:,:,:,:)[:,:,:]
   real(R8B),allocatable       :: rho_buf(:,:,:)[:,:,:]
   integer,allocatable         :: state_buf(:,:,:)[:,:,:]
#endif 

   lx = lb_dom%lx(1)
   ly = lb_dom%lx(2)
   lz = lb_dom%lx(3)


   ab = prc%bnd(:,:,:,:,1)
   ae = prc%bnd(:,:,:,:,2)


#ifdef INIT_WITH_ROOT
#ifdef USE_CAF
   allocate(rho_buf(cobnd(1),cobnd(2),cobnd(3))[prc%np(1),prc%np(2),*])
   allocate(u_buf(LB_NODE(NDIM,cobnd(1),cobnd(2),cobnd(3)))[prc%np(1),prc%np(2),*])
   allocate(state_buf(cobnd(1),cobnd(2),cobnd(3))[prc%np(1),prc%np(2),*])
   ! for caf worker threads: have to copy local data to buffers
    if (prc%rk /= 0 ) then ! copy local arrays to coarray buffers
      u_buf(LB_NODE(1:NDIM,1:lx,1:ly,1:lz))   = lb_dom%u(LB_NODE(1:NDIM,1:lx,1:ly,1:lz))
      state_buf(1:lx,1:ly,1:lz)      = lb_dom%state(1:lx,1:ly,1:lz)
      rho_buf(1:lx,1:ly,1:lz)      = lb_dom%rho(1:lx,1:ly,1:lz)
   end if
#endif 

   call mpl_barrier()

  ! Master Thread gathers all data
  ! and executes output
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
                   k == (prc%root_crd(3)+1)) then ! is master thread
               ! copy values of local arrays to global array


                  lb_dom%grho(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                     &        prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                     &        prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)) &
                     &  = lb_dom%rho(1:lx,1:ly, 1:lz)

               lb_dom%gu(&
               LB_NODE(:,ab(i,j,k,1):ae(i,j,k,1),ab(i,j,k,2):ae(i,j,k,2),ab(i,j,k,3):ae(i,j,k,3))) &
                       = lb_dom%u(LB_NODE(:,1:lx,1:ly,1:lz))

#ifdef SPONGE1
                 lb_dom%gomega(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                     &        prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                     &        prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)) &
                     & =  lb_dom%omega(1:lx,1:ly,1:lz)
#endif /* SPONGE1 */
                 lb_dom%gstate(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                     &        prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                     &        prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)) &
                     & =  lb_dom%state(1:lx,1:ly,1:lz)

                       s_par%g_x(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2)) =  lb_dom%x(1:lx)
                       s_par%g_y(prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2)) =  lb_dom%y(1:ly)
                       s_par%g_z(prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)) =  lb_dom%z(1:lz)
                else ! if master thread's loop is at another crd than 0,0,0 -> recv from other threads 
#ifdef USE_MPI
                call MPI_CART_RANK(prc%cart_comm,(/(i-1),(j-1),(k-1)/),target_rank,prc%ierr)
                call MPI_RECV(lb_dom%gu(LB_NODE(1:NDIM,b_l(1):b_u(1),b_l(2):b_u(2),b_l(3):b_u(3))),&
                                    length*NDIM,MPI_DOUBLE_PRECISION,&
                                    target_rank,tag_rho,prc%cart_comm,prc%stat,prc%ierr)
                call MPI_RECV(lb_dom%grho(b_l(1):b_u(1),b_l(2):b_u(2),b_l(3):b_u(3)),length,mpi_double_precision,&
                                    target_rank,tag_rho,prc%cart_comm,prc%stat,prc%ierr)
                call MPI_RECV(lb_dom%gstate(b_l(1):b_u(1),b_l(2):b_u(2),b_l(3):b_u(3)),length,MPI_INTEGER,&
                                    target_rank,tag_rho,prc%cart_comm,prc%stat,prc%ierr)
#ifdef SPONGE1
                call MPI_RECV(lb_dom%gomega(b_l(1):b_u(1),b_l(2):b_u(2),b_l(3):b_u(3)),length,mpi_double_precision,&
                                    target_rank,tag_rho,prc%cart_comm,prc%stat,prc%ierr)
#endif /* SPONGE1 */
#endif /* USE_MPI */
#ifdef USE_CAF
                lb_dom%gu(&
            LB_NODE(:,ab(i,j,k,1):ae(i,j,k,1),ab(i,j,k,2):ae(i,j,k,2),ab(i,j,k,3):ae(i,j,k,3))) &
                       =     u_buf(LB_NODE(:,1:lx,1:ly,1:lz))[i,j,k]

                 lb_dom%grho(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                              prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                              prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)) &
                       =      rho_buf(1:lx,1:ly,1:lz)[i,j,k]


                 lb_dom%gstate(prc%bnd(i,j,k,1,1):prc%bnd(i,j,k,1,2), & 
                              prc%bnd(i,j,k,2,1):prc%bnd(i,j,k,2,2), &
                              prc%bnd(i,j,k,3,1):prc%bnd(i,j,k,3,2)) &
                       =      state_buf(1:lx,1:ly,1:lz)[i,j,k]


#endif 
                end if
            end do  
         end do  
      end do  
#ifdef TEC_OUT
      call write_tec_output(lb_dom,countOut,tStep,s_par,prc)
#endif
 else ! If current rank is not master -> send local array values to master 
#ifdef USE_MPI
       ! all other threads: 
       ! wait until mpi_recv is done
      length = lx*ly*lz !(b_u(1)-b_l(1) + 1) * (b_u(2)-b_l(2)+1) * (b_u(3)-b_l(3)+1)
      call MPI_SEND(lb_dom%u(LB_NODE(1:NDIM,1:lx,1:ly,1:lz)),length*NDIM,           &
                     MPI_DOUBLE_PRECISION,prc%root_th,tag_rho,prc%cart_comm,prc%ierr)
      call MPI_SEND(lb_dom%rho(1:lx,1:ly,1:lz),length,                              &
                     MPI_DOUBLE_PRECISION,prc%root_th,tag_rho,prc%cart_comm,prc%ierr)
      call MPI_SEND(lb_dom%state(1:lx,1:ly,1:lz),length,MPI_INTEGER,                &
                     prc%root_th,tag_rho,prc%cart_comm,prc%ierr)
#ifdef SPONGE1
      call MPI_SEND(lb_dom%omega(1:lx,1:ly,1:lz),length,mpi_double_precision,       &
                     prc%root_th,tag_rho,prc%cart_comm,prc%ierr)
#endif /* SPONGE1 */
#endif /* USE_MPI */
   end if
#ifdef USE_CAF
   deallocate(rho_buf)
   deallocate(u_buf)
   deallocate(state_buf)
#endif
#else /* INIT_WITH_ROOT */
#ifdef TEC_OUT
   call write_tec_output(lb_dom,countOut,tstep,s_par,prc)
#endif
#endif /* INIT_WITH_ROOT */
end subroutine gather_output

!------------------------------------------------------------------------

#ifdef TEC_OUT






subroutine write_tec_output(lb_dom,countOut,tstep,s_par,prc)

   !------------------------------------------------------------------------
   !
   ! output results to file
   !


   INCLUDE 'tecio.f90'
   type(lb_block),intent(inout)  :: lb_dom
   type(sim_parameter),intent(inout):: s_par 
   type(mpl_var)   :: prc
   integer, intent(inout) :: countOut, tstep
   character(16)              :: cLx,cRe,cumax,cLy,cLz,cRk,cStep,clchar
   character(64)              :: filename 
   Integer IMax,JMax,KMax
   character*1 NULLCHR
   Integer   Debug,III   !,NPts,NElm
   real(R4B),allocatable :: X(:,:,:), Y(:,:,:), Z(:,:,:), P(:,:,:)
   Real(R8B)    SolTime,offx(3)
   Integer VIsDouble, FileType
   Integer ZoneType,StrandID,ParentZn,IsBlock
   Integer ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
   pointer (NullPtr,Null)
   Integer Null(*)
   Integer I,J,K



   write(unit=cLx,fmt='(i4)') s_par%gx(1)
   write(unit=cRe,fmt='(f8.5)') s_par%omega
   write(unit=cumax,fmt='(f5.3)') s_par%umax
   write(unit=cLy,fmt='(i4)') s_par%gx(2)
   write(unit=cLz,fmt='(i4)') s_par%gx(3)
   write(unit=cRk,fmt='(i4)') prc%rk
   write(unit=cStep,fmt='(i6)') tstep
   write(unit=clchar,fmt='(f6.3)') s_par%obst_l

      NULLCHR = CHAR(0)
      NullPtr = 0
      Debug   = 0
      FileType = 0
      VIsDouble = 0
#ifdef INIT_WITH_ROOT
      IMax    = s_par%gx(1)
      JMax    = s_par%gx(2)
      KMax    = s_par%gx(3)
      offx    = 0.0d0
#else
      IMax    = lb_dom%lx(1)
      JMax    = lb_dom%lx(2)
      KMax    = lb_dom%lx(3)
      offx(1) = prc%crd(1)*real(lb_dom%lx(1))
      offx(2) = prc%crd(2)*real(lb_dom%lx(2))
      offx(3) = prc%crd(3)*real(lb_dom%lx(3))
#endif

      ZoneType = 0
      SolTime = real(tstep)
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0
      allocate(X(IMax,JMax,KMax))
      allocate(Y(IMax,JMax,KMax))
      allocate(Z(IMax,JMax,KMax))
      allocate(P(IMax,JMax,KMax))



if(countOut == 1) then
   s_par%problem_name = adjustl(s_par%problem_name)
   filename =  trim(s_par%problem_name)//"_"//&
         trim(adjustl(s_par%modelname))//'_'//&
         'Re'//trim(adjustl(cRe))//'_'//&
         'Lx_'// trim(adjustl(cLx))//'_'//&
         'u'// trim(adjustl(cumax))//'_'//&
#ifdef CRVP
         'l'// trim(adjustl(clchar))//'_'//&
#endif
         trim(adjustl(cRk))//'.plt'//NULLCHR

      ! only write header information if the file has not been opened yet
      !
      !... Open the file and write the tecplot datafile 
      !... header information.
      !
      I = TecIni112('SIMPLE DATASET'//NULLCHR, &
#ifdef D2Q9
         'X Y Z RHO UX UY STATE'//NULLCHR, &
#endif
#ifdef D3Q19
         'X Y Z RHO UX UY UZ STATE'//NULLCHR, &
#endif
         filename//NULLCHR,&
         '.'//NULLCHR, &
         FileType, &
         Debug, &
         VIsDouble)
endif

!
!... Write the zone header information.
!
      I = TecZne112(trim(adjustl(cStep))//'tStep'//NULLCHR, &
                    ZoneType, &
                    IMax,  JMax, KMax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0,  0,  0,  Null,  Null,  Null, &
                    ShrConn)
!
!... Write out the field data.
!

   do k=1,KMax
      do j=1,JMax
         do i=1,IMax
            X(i,j,k)=real(i,4)+offx(1)
            Y(i,j,k)=real(j,4)+offx(2)
            Z(i,j,k)=real(k,4)+offx(3)
#ifdef INIT_WITH_ROOT
#ifdef D2Q9
            P(i,j,k)=real(lb_dom%gu(LB_NODE(1,i,j,k))**2 + lb_dom%gu(LB_NODE(2,i,j,k))**2 ,4)
#endif
#ifdef D3Q19
            P(i,j,k)=real(lb_dom%gu(LB_NODE(1,i,j,k))**2 + lb_dom%gu(LB_NODE(2,i,j,k))**2 +lb_dom%gu(LB_NODE(3,i,j,k))**2     ,4)
#endif

#else /* INIT_WITH_ROOT */
#ifdef D2Q9
            P(i,j,k)=real(lb_dom%u(LB_NODE(1,i,j,k))**2 + lb_dom%u(LB_NODE(2,i,j,k))**2 ,4)
#endif
#ifdef D3Q19
            P(i,j,k)=real(lb_dom%u(LB_NODE(1,i,j,k))**2 + lb_dom%u(LB_NODE(2,i,j,k))**2 +lb_dom%u(LB_NODE(3,i,j,k))**2     ,4)
#endif
#endif /* INIT_WITH_ROOT */
         enddo
      enddo
   enddo
      III = IMax*JMax*KMax
      I   = TecDat112(III,X,0)
      I   = TecDat112(III,Y,0)
      I   = TecDat112(III,Z,0)
!      I   = TecDat112(III,P,0)
#ifdef INIT_WITH_ROOT
      I   = TecDat112(III,real(lb_dom%grho(1:IMax,1:JMax,1:KMax),4),0)
      I   = TecDat112(III,real(lb_dom%gu(LB_NODE(1,1:IMax,1:JMax,1:KMax)),4),0)
      I   = TecDat112(III,real(lb_dom%gu(LB_NODE(2,1:IMax,1:JMax,1:KMax)),4),0)
#ifdef D3Q19
      I   = TecDat112(III,real(lb_dom%gu(LB_NODE(3,1:IMax,1:JMax,1:KMax)),4),0)
#endif
#ifdef SPONGE
     I   = TecDat112(III,real(lb_dom%gomega(1:IMax,1:JMax,1:KMax),4),0)
#else
     I   = TecDat112(III,real(lb_dom%gstate(1:IMax,1:JMax,1:KMax),4),0)
#endif /* SPONGE */
#else
      I   = TecDat112(III,real(lb_dom%rho(1:IMax,1:JMax,1:KMax),4),0)
      I   = TecDat112(III,real(lb_dom%u(LB_NODE(1,1:IMax,1:JMax,1:KMax)),4),0)
      I   = TecDat112(III,real(lb_dom%u(LB_NODE(2,1:IMax,1:JMax,1:KMax)),4),0)
#ifdef D3Q19
      I   = TecDat112(III,real(lb_dom%u(LB_NODE(3,1:IMax,1:JMax,1:KMax)),4),0)
#endif
     I   = TecDat112(III,real(lb_dom%state(1:IMax,1:JMax,1:KMax),4),0)
!     I   = TecDat112(III,real(lb_dom%omega(1:IMax,1:JMax,1:KMax),4),0)
#endif /* INIT_WITH_ROOT */
      deallocate(Y)
      deallocate(Z)
      deallocate(P)

if(s_par%goend .eqv. .true.) then
      I = TecEnd112()
endif

end subroutine 
#endif /* TEC_OUT */










   subroutine read_params(s_par,filename,prc)


   !------------------------------------------------------------------------
   !
   ! Simulation parameters are read from the file lbc.params
   ! These include
   ! umax         the max velocity, here: of the lid
   ! Re           Reynolds number 
   ! tMax         max. timesteps), tOut (interval for writing out results
   ! problem      1=cylinder, 2=cavity
   ! s_par%gx(1),s_par%gx(2),s_par%gx(3)  global domain size
   ! s_par%setObstacles set small obstacles in equal distance in the whole domain
   ! prc%np(1:3)   number of processes/threads in each direction (cartesian topolos_par%gx(2))
   !

      type(sim_parameter),intent(inout)   :: s_par
      type(mpl_var)   :: prc
      character(10), intent(in)           :: filename
      character(10)           :: temp_name 
      logical :: file_exists     

      !------------------------------------
      ! Reset some variables and set Name variables for output

      !   Set rho0 
      s_par%rho0 = 1.0d0 

      s_par%goend        = .false.
      s_par%initial      = .false.
      s_par%comm_method =  ""//&
#if defined SENDRECV_BUF
     "sendrecv_buf" &
#elif defined ISEND_IRECV_BUF
 "isend_irecv_buf" &
#elif defined ISEND_IRECV
 "isend_irecv" &
#elif defined SENDRECV
 "sendrecv" &
#elif defined MPI_SUBARRAY 
 "subarray" &
#elif defined USE_ADCL
 "adcl" &
#elif defined USE_CAF
 "caf" &
#elif defined PERSISTENT 
 "persistent" &
#elif defined PERSISTENT2
 "persistent" &
#else
 "serial" &
#endif
#if defined MPI_OLD
//"_old"&
#endif
#if defined BUFFER_EACH  
//"B" &
#endif
#if defined COMM_REDUCED 
//"R" &
#endif
//""

#if defined D2Q9
      s_par%modelname = "d2q9"
#elif defined D3Q19
      s_par%modelname = "d3q19"
#endif
      temp_name       = s_par%modelname
#ifdef MRT 
      s_par%modelname = trim(temp_name)//'_mrt'
#else
#ifdef TRT 
      s_par%modelname = trim(temp_name)//'_trt'
#endif
#endif
#ifdef SPONGE
      s_par%modelname = "spg_"//s_par%modelname
#endif
#if defined LAYOUT_LIJK
      s_par%layoutname = "lijk"
#elif defined LAYOUT_IJKL
      s_par%layoutname = "ijkl"
#endif

      !------------------------------------
      ! Read in parameter file

#ifdef USE_MPI
      if(prc%rk == 0) then
#endif
        ! Root thread reads in simulation parameters
         inquire(file=filename,exist=file_exists)
         if(file_exists .eqv. .false.) then
            write(*,*) "Error:"
            write(*,*) "Parameter file '",filename,"' does not exist. "
            write(*,*) "Exiting ... "  
            stop
         endif
         open(1,file=filename)
         read(1,*) s_par%umax                               ! boundary condition, lid speed (0.01<u<0.1)
         read(1,*) s_par%omega                              ! Reynolds number Re<5000?
         read(1,*) s_par%tMax,s_par%tOut,s_par%tRamp,s_par%tOutBeg,s_par%tInitial   ! max timesteps, output interval,ramp tmax, tOutBegin, tInitial
         read(1,*) s_par%problem                            
! problem: 1 cylinder in channel, 2 cavity, 3 shock,gaussian=4, cavity_same=5 
         read(1,*) s_par%gx(1),s_par%gx(2),s_par%gx(3)      ! global domain size x, y, z
         read(1,*) s_par%setObstacles                       ! set additional obstacles in domain t / f
         read(1,*) prc%np(1),prc%np(2),prc%np(3)            ! # processes in each direction
         read(1,*) s_par%save_output                        ! save output t / f
         read(1,*) s_par%obst_l                             ! Additional obstacle Parameter. For CRVP: Distance of vortices 
         close(1)

      !------------------------------------
      ! set characteristic lengths and problem lengths (cylinder radius etc.) 

         if(NDIM==2) s_par%gx(3) = 1
     ! Set position for obstacle (only problem cylinder)
         if(s_par%problem==cylinder) then
         s_par%obst_x=s_par%gx(1)/5+1
         s_par%obst_y=int(real(s_par%gx(2))*(7.d0/12.d0))
         s_par%obst_z=int(real(s_par%gx(3))*(1.d0/2.d0))
         s_par%obst_r=int(real(s_par%gx(2))*(1.d0/15.d0))+1
         elseif(s_par%problem==flute) then
         s_par%obst_x=s_par%gx(1)/5
         s_par%obst_r=int(real(s_par%gx(2))/90.)
         s_par%obst_y=s_par%gx(2)/4-s_par%obst_r
         s_par%obst_z=0
         elseif(s_par%problem==corotating_vortex) then
         s_par%obst_r = s_par%gx(1)/100+1
         elseif(s_par%problem==gaussian .or. s_par%problem==gauss_convect) then
            if(s_par%obst_l < 1.0) then
               s_par%obst_r = s_par%gx(1)/100 +1
            else 
               s_par%obst_r = int(s_par%obst_l) 
            endif
         endif    
     ! Set sponge layer thickness
         s_par%sponge_size = 20
 
     ! set characteristic length lchar
         if(s_par%problem==cylinder .or. s_par%problem==flute) then
            s_par%lchar = real(s_par%obst_r)*2.
         else
         if(NDIM==2) then
            s_par%lchar = real(MAX(s_par%gx(1),s_par%gx(2)))
         else
            s_par%lchar = real(MAX(s_par%gx(1),s_par%gx(2),s_par%gx(3)))
         endif
      endif


      !------------------------------------
      ! set problem names for output 

      if(s_par%problem == cylinder) then
         s_par%problem_name = "cylinder"
      elseif(s_par%problem == cavity) then
         s_par%problem_name = "cavity"
      elseif(s_par%problem == shock) then
         s_par%problem_name = "shock"
      elseif(s_par%problem == gaussian) then
         s_par%problem_name = "gaussian"
      elseif(s_par%problem == gauss_convect) then
         s_par%problem_name = "gauss_conv"
      elseif(s_par%problem == cavity_same) then
         s_par%problem_name = "cavity_same"
      elseif(s_par%problem == flute ) then
         s_par%problem_name = "flute"
      elseif(s_par%problem == corotating_vortex ) then
         s_par%problem_name = "crvp"
      elseif(s_par%problem == cylinder_impulse ) then
         s_par%problem_name = "cyl_imp"
      elseif(s_par%problem == plate_trail ) then
         s_par%problem_name = "plate_trail"
      elseif(s_par%problem == taylor_vortex) then
         s_par%problem_name = "taylor_vortex"
      elseif(s_par%problem == planar_standing_wave) then
         s_par%problem_name = "plane_wave"
      else
         s_par%problem_name = ""
      endif

      ! Calculate Viscosity
       s_par%nu = (2./s_par%omega-1)/6.0
      ! Calculate Re
       s_par%Re = s_par%lchar*s_par%umax/s_par%nu 



      !------------------------------------
      ! Output information about parameters 

      if(prc%rk == 0) then
         ! print out some info 
         write(*,*) 'PARAMETERS    '
         write(*,*) '-----------------------------'
         write(*,'(a12,a)') 'Problem:     ',  s_par%problem_name
         write(*,'(a12,i3,i3,i3)') 'prc%np       ',  prc%np(1),prc%np(2),prc%np(3)
         write(*,'(a12,f6.4)')     'umax        ',   s_par%umax
         write(*,'(a12,f10.4)')     'Re          ',   s_par%Re
         write(*,'(a12,f15.4)')     'Re_acc      ',   s_par%lchar*cs/s_par%nu
         write(*,'(a12,f6.4)')     'nu          ',   s_par%nu
         write(*,'(a12,f8.5)')     'omega       ',   s_par%omega
         write(*,'(a12,i4,i4,i4)') 'Problem size',   s_par%gx(1),s_par%gx(2),s_par%gx(3)
         write(*,'(a12,i8)')     'Max Tstep   ',   s_par%tMax
         write(*,'(a12,i8)')     'Out Tstep   ',   s_par%tOut
         write(*,'(a12,i8)')     'Ramp until  ',   s_par%tRamp
         write(*,'(a12,i8)')     'Beg Out at  ',   s_par%tOutBeg
         write(*,'(a12,i8)')     'Init until  ',   s_par%tInitial
         if(s_par%setObstacles .EQV. .TRUE.) write(*,*) 'Insert regular obstacles ', s_par%setObstacles 
#ifdef D2Q9
         write(*,*) "Model:    D2Q9"
#endif
#ifdef D3Q19
         write(*,*) "Model:    D3Q19"
#endif
         write(*,*) "Relaxation is: ", &
#ifdef MRT
"mrt"
#else
#ifdef TRT
"trt"
#else
"bgk"
#endif
#endif

         write(*,*) "----------------------------------------------"
         write(*,*) 
         if(s_par%tInitial>2) then
            write(*,*) "omega for initial relaxation (u fixed): ",1.1," ...  ",s_par%omega
            write(*,*) "----------------------------------------------"
            write(*,*) 
         endif
     end if
#ifdef USE_MPI
  end if

      !------------------------------------
      ! Broadcast input variables 

   call MPI_BCAST(prc%np,3,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%omega,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%nu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%Re,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%umax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%obst_x,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%obst_y,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%obst_z,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%obst_r,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%gx,3,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%setObstacles,1,MPI_LOGICAL,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%save_output ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%tMax,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%tOut,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%tRamp,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%tOutBeg,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%tInitial,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%problem,1,MPI_INTEGER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%problem_name,16,MPI_CHARACTER,0,MPI_COMM_WORLD,prc%ierr)
   call MPI_BCAST(s_par%modelname,16,MPI_CHARACTER,0,MPI_COMM_WORLD,prc%ierr)
#endif /* USE_MPI */
   end subroutine read_params
!------------------------------------------------------------------------








subroutine write_performance_results(result_state,s_par,prc,meas)
!
!  Performance results (MLUPs) are written to a file for statistical review
!  The file-format is simple can be evaluated with gnuplot
!
   implicit none
   type(sim_parameter)  :: s_par
   type(mpl_var)  :: prc  
   type(measure)  :: meas
   integer  :: result_state,lxyztot
   logical       :: file_exists
   character(15) :: filename

  filename='performance.res'
  inquire(file=filename,exist=file_exists)
  open(11,file=filename,position='append') !TRIM(ADJUSTL(filename)))
  if(file_exists .EQV. .false.) then
     write(11,*) '#loglx*ly*lz lx*ly*lz MLUPs  lx   ly   lz   nprc  tTot   tComm  %Comm mpi_meth    layo  resok'
  end if
   lxyztot = s_par%gx(1)*s_par%gx(2)*s_par%gx(3)
  write(11,'(f6.2,i11,f8.2,4i5,3f8.2,a1,a14,a5,a6,i3)') log10(dble(lxyztot)),       &
            & int(s_par%gx(1)*s_par%gx(2)*s_par%gx(3),4),real(meas%mlups,4),    &
            & int(s_par%gx(1),2),int(s_par%gx(2),2), int(s_par%gx(3),2),                                             &
            & int(prc%size,2),meas%mean_total,meas%mean_comm,meas%mean_comm/meas%mean_total*100.,' ',            &
            & s_par%comm_method,s_par%layoutname,s_par%modelname,result_state
  close(11)
end subroutine write_performance_results
!------------------------------------------------------------------------






subroutine write_temp_data(tStep,nr_unit,size_x,size_y,size_z,size_nnod,fIn,limit_x,limit_y,limit_z,sx,sy,sz)

   integer, intent(in) :: nr_unit,size_x,size_y,size_z,size_nnod,tStep,limit_x,limit_y,limit_z,sx,sy,sz
   real(R8B), intent(in)    :: fIn(LB_NODE(size_nnod,0:size_x+1,0:size_y+1,0:size_z+1))
   integer             :: i,j,k,l
   write(nr_unit,*) "Timestep",tStep
   do k=sz,limit_z
      do j=sy,limit_y
         do i=sx,limit_x
            do l=1,size_nnod
               write(nr_unit,'(3i3,f23.10)') i+1-sx,j+1-sy,k+1-sz,fIn(LB_NODE(l,i,j,k))
            end do 
         end do 
      end do
   end do
end subroutine write_temp_data







subroutine  dealloc_mem(lb_dom,s_par)

!------------------------------------------------------------------------
! Deallocate allocated arrays
!


   type(lb_block)             :: lb_dom
   type(sim_parameter)        :: s_par 
   integer                    :: alloc_stat
#ifdef INIT_WITH_ROOT
   deallocate(lb_dom%gfIn,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%gfIn"
   deallocate(lb_dom%gu,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%gu"
   deallocate(lb_dom%grho,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%grho"
   deallocate(lb_dom%gstate,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%g_state"
#ifdef SPONGE
   deallocate(lb_dom%gomega,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%gsponge"
#endif /* SPONGE */
   deallocate( s_par%g_x,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%g_x"
   deallocate( s_par%g_y,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%g_y"
   deallocate( s_par%g_z,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%g_z"
#endif /* INIT_WITH_ROOT */
   deallocate(lb_dom%obs,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating obs       "
!   deallocate(lb_dom%obs_sponge,stat=alloc_stat)
!      if(alloc_stat /=0) write(*,*) "Error deallocating obs_sponge"
!   deallocate(lb_dom%obs_reset,stat=alloc_stat)
!      if(alloc_stat /=0) write(*,*) "Error deallocating obs_reset "
   deallocate(lb_dom%obs_nrwall,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating obs_nrwall"
   deallocate(lb_dom%nrwall_0val,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating nr_wall0val"

   deallocate(lb_dom%fIn,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%fin "
   deallocate(lb_dom%fOut,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%fOut"
   deallocate(lb_dom%u,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%u  "
if(s_par%tInitial > 0) then
   deallocate(lb_dom%u0,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%u0  "
endif
   deallocate(lb_dom%rho,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%rho"
   deallocate(lb_dom%state,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%state"
#ifdef SPONGE
   deallocate(lb_dom%omega,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%omega"
#endif /* SPONGE */
   deallocate(lb_dom%x,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%x  "
   deallocate(lb_dom%y,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%y  "
   deallocate(lb_dom%z,stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) "Error deallocating lb_dom%z  "
#ifdef SPARSE
   deallocate(lb_dom%fsIn,stat=alloc_stat)
#endif

end subroutine dealloc_mem






#ifdef WRITE_TSTATUS
   subroutine write_tStat(tStep,duration,total,lb_dom,s_par,prc)


      implicit none
      type(mpl_var)   :: prc
      type(lb_block)   :: lb_dom
      type(sim_parameter),intent(inout) :: s_par
      real(R8B) :: duration,total 
      integer  :: tStep
      character(len=8)   :: datechar
      character(len=10)  :: timechar
      character(len=100) :: filename

   call mpl_barrier()

   if(prc%rk == 0) then
      if(tStep==WRITE_TSTATUS) then
         call date_and_time(datechar, timechar)
   write(unit=filename,fmt='(4a)') "tstep_",trim(s_par%comm_method),"_",timechar
         open(42,file=trim(filename),position='append')
         write(42,*) "#"
         write(42,*) "# LBC run lx",lb_dom%lx
         write(42,*) "# Date ccyymmdd   ",datechar,"   at Time hhmmss  ",timechar 
         write(42,*) "# Comm.-Method   ",s_par%comm_method
         write(42,*) "# "
         write(42,*) "# tStep        Duration ts        Duration total"
      endif
      write(42,'(i12,2f16.6)') tStep,duration,total
   end if

   call mpl_barrier()

   end subroutine write_tStat 

#endif /* WRITE_TSTATUS */





   subroutine write_status(meas,result_state,cur_tstep,s_par,prc)

      implicit none
#ifdef USE_CAF
      type(measure),intent(inout) :: meas[*]
#else
      type(measure),intent(inout) :: meas
#endif
      type(sim_parameter),intent(inout) :: s_par
      type(mpl_var)   :: prc
      integer,intent(in)  :: result_state,cur_tstep
      real(R8B) :: prop_min,prop_max,prop_sum,comm_min,comm_max,comm_sum
      real(R8B) :: sbnd_min,sbnd_max,sbnd_sum,coll_min,coll_max,coll_sum
      real(R8B) :: mpl_cp_buf_max,mpl_cp_buf_min,mpl_cp_buf_sum,mpl_cp_buf_back_max,mpl_cp_buf_back_sum
      real(R8B) :: mpl_cp_buf_back_min,mpl_exch_max,mpl_exch_min,mpl_exch_sum
      real(R8B) :: mpl_sync_max,mpl_sync_min,mpl_sync_sum
      real(R8B) :: totl_min,totl_max,totl_sum
      real(R8B) :: mpl_comm_size_min,mpl_comm_size_max,mpl_comm_size_sum 
#if defined USE_CAF
      integer :: i
#endif

      comm_min             = 1.e37 
      mpl_cp_buf_min       = 1.e37 
      mpl_exch_min         = 1.e37 
      mpl_cp_buf_back_min  = 1.e37 
      mpl_exch_min         = 1.e37 
      mpl_sync_min         = 1.e37 
      mpl_comm_size_min    = 1.e37 
      prop_min             = 1.e37 
      sbnd_min             = 1.e37 
      coll_min             = 1.e37 

      comm_max             = 0.
      mpl_cp_buf_max       = 0.
      mpl_exch_max         = 0.
      mpl_cp_buf_back_max  = 0.
      mpl_exch_max         = 0.
      mpl_sync_max         = 0.
      mpl_comm_size_max    = 0.
      prop_max             = 0.
      sbnd_max             = 0.
      coll_max             = 0.

      !-------------------------------
      ! get the min and max values of measured communication times

#ifdef USE_MPI
      call mpi_reduce(meas%duration,totl_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%duration,totl_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%duration,totl_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%sbnd_duration,sbnd_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%sbnd_duration,sbnd_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%sbnd_duration,sbnd_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%coll_duration,coll_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%coll_duration,coll_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%coll_duration,coll_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%prop_duration,prop_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%prop_duration,prop_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%prop_duration,prop_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%comm_duration,comm_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%comm_duration,comm_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%comm_duration,comm_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf,mpl_cp_buf_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf,mpl_cp_buf_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf,mpl_cp_buf_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf_back,mpl_cp_buf_back_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf_back,mpl_cp_buf_back_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_cp_buf_back,mpl_cp_buf_back_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_exch,mpl_exch_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_exch,mpl_exch_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_exch,mpl_exch_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_sync,mpl_sync_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_sync,mpl_sync_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_sync,mpl_sync_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_comm_size,mpl_comm_size_min,1,mpi_double_precision,mpi_min,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_comm_size,mpl_comm_size_max,1,mpi_double_precision,mpi_max,0,prc%cart_comm,prc%ierr)
      call mpi_reduce(meas%mpl_comm_size,mpl_comm_size_sum,1,mpi_double_precision,mpi_sum,0,prc%cart_comm,prc%ierr)


#elif defined USE_CAF
      do i = 1,prc%size
         if(meas[i]%sbnd_duration < sbnd_min) sbnd_min =  meas[i]%sbnd_duration
         if(meas[i]%sbnd_duration > sbnd_max) sbnd_max =  meas[i]%sbnd_duration
         if(meas[i]%sbnd_duration > 0.d0    ) sbnd_sum =  meas[i]%sbnd_duration + sbnd_sum
         if(meas[i]%coll_duration < coll_min) coll_min =  meas[i]%coll_duration
         if(meas[i]%coll_duration > coll_max) coll_max =  meas[i]%coll_duration
         if(meas[i]%coll_duration > 0.d0    ) coll_sum =  meas[i]%coll_duration + coll_sum
         if(meas[i]%prop_duration < prop_min) prop_min =  meas[i]%prop_duration
         if(meas[i]%prop_duration > prop_max) prop_max =  meas[i]%prop_duration
         if(meas[i]%prop_duration > 0.d0    ) prop_sum =  meas[i]%prop_duration + prop_sum
         if(meas[i]%comm_duration < comm_min) comm_min =  meas[i]%comm_duration
         if(meas[i]%comm_duration > comm_max) comm_max =  meas[i]%comm_duration
         if(meas[i]%comm_duration > 0.d0    ) comm_sum =  meas[i]%comm_duration + comm_sum
         if(meas[i]%mpl_cp_buf    < mpl_cp_buf_min) mpl_cp_buf_min =  meas[i]%mpl_cp_buf
         if(meas[i]%mpl_cp_buf    > mpl_cp_buf_max) mpl_cp_buf_max =  meas[i]%mpl_cp_buf
         if(meas[i]%mpl_cp_buf    > 0.d0    ) mpl_cp_buf_sum =  meas[i]%mpl_cp_buf + mpl_cp_buf_sum
         if(meas[i]%mpl_cp_buf_back < mpl_cp_buf_back_min) mpl_cp_buf_back_min =  meas[i]%mpl_cp_buf_back
         if(meas[i]%mpl_cp_buf_back > mpl_cp_buf_back_max) mpl_cp_buf_back_max =  meas[i]%mpl_cp_buf_back
         if(meas[i]%mpl_cp_buf_back > 0.d0    ) mpl_cp_buf_back_sum =  meas[i]%mpl_cp_buf_back + mpl_cp_buf_back_sum
         if(meas[i]%mpl_exch      < mpl_exch_min) mpl_exch_min =  meas[i]%mpl_exch
         if(meas[i]%mpl_exch      > mpl_exch_max) mpl_exch_max =  meas[i]%mpl_exch
         if(meas[i]%mpl_exch      > 0.d0    ) mpl_exch_sum =  meas[i]%mpl_exch + mpl_exch_sum
         if(meas[i]%mpl_sync      < mpl_sync_min) mpl_sync_min =  meas[i]%mpl_sync
         if(meas[i]%mpl_sync      > mpl_sync_max) mpl_sync_max =  meas[i]%mpl_sync
         if(meas[i]%mpl_sync      > 0.d0    ) mpl_sync_sum =  meas[i]%mpl_sync      + mpl_sync_sum
         if(meas[i]%mpl_comm_size < mpl_comm_size_min) mpl_comm_size_min =  meas[i]%mpl_comm_size
         if(meas[i]%mpl_comm_size > mpl_comm_size_max) mpl_comm_size_max =  meas[i]%mpl_comm_size
         if(meas[i]%mpl_comm_size > 0.d0    ) mpl_comm_size_sum =  meas[i]%mpl_comm_size + mpl_comm_size_sum
      enddo
#endif

! Mean values
   if(prc%rk == 0) then
      meas%mean_total = totl_sum/real(prc%size)
      meas%mean_sbnd  = sbnd_sum/real(prc%size)
      meas%mean_coll  = coll_sum/real(prc%size)
      meas%mean_prop  = prop_sum/real(prc%size)
      meas%mean_comm  = comm_sum/real(prc%size)
      meas%mean_mpl_cp_buf  = mpl_cp_buf_sum/real(prc%size)
      meas%mean_mpl_cp_buf_back  = mpl_cp_buf_back_sum/real(prc%size)
      meas%mean_mpl_exch  = mpl_exch_sum/real(prc%size)
      meas%mean_mpl_sync  = mpl_sync_sum/real(prc%size)
      meas%mean_mpl_comm_size  = mpl_comm_size_sum/real(prc%size)
      meas%max_total = totl_max

      meas%mlups = real(cur_tstep,8)*real(s_par%gx(1),8)*real(s_par%gx(2),8)*real(s_par%gx(3),8)/ meas%mean_total/10**6
      write(*,*) 
      write(*,*) "  Statistics" 
      write(*,*) 
      write(*,'(a35,a32)')          "Communication method:           ",s_par%comm_method
      write(*,'(a35,a32)')          "Memory layout of densities:     ",s_par%layoutname
      write(*,*) 
      write(*,'(a35,e10.3)')        "Max. Total wall time:           ", meas%max_total
      write(*,*) 
      write(*,'(a35,f8.3,a12,f8.3,a2)') "Mean Total duration:            ", &
      & meas%mean_total," s     (th0: " , meas%duration, " )" !meas%comm_duration+meas%sbnd_duration+meas%cceq_duration+ &
    !  & meas%coll_duration+meas%prop_duration+meas%ccmv_duration+meas%bncb_duration
      write(*,'(a35,e8.3,a2)') "initialization duration         ", meas%init_duration," s       min     max"
#if defined USE_MPI || defined USE_CAF
      if(prc%np(1)*prc%np(2)*prc%np(3) > 1) then
      write(*,'(a35,f8.3,a2)') "Communication duration          ", &
      & meas%mean_comm," s"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "   buffer copy                  ", &
      & meas%mean_mpl_cp_buf   ," s", meas%mean_mpl_cp_buf/meas%mean_comm*100.,"%    (", & 
      & mpl_cp_buf_min,mpl_cp_buf_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "   buffer exchange              ", &
      & meas%mean_mpl_exch     ," s", meas%mean_mpl_exch/meas%mean_comm*100.,"%    (",mpl_exch_min,mpl_exch_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "   synchronization              ", &
      & meas%mean_mpl_sync     ," s", meas%mean_mpl_sync/meas%mean_comm*100.,"%    (",mpl_sync_min,mpl_sync_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "   buffer copy back             ", &
      & meas%mean_mpl_cp_buf_back  ," s", meas%mean_mpl_cp_buf_back/meas%mean_comm*100.,"%    (", & 
      & mpl_cp_buf_back_min,mpl_cp_buf_back_max,")"

      write(*,*) 
      write(*,'(a35,f12.6,a15,f6.3,a4,f6.3)') "   Data size                    ", &
      & meas%mean_mpl_comm_size/1000000.d0,"MB      min,",mpl_comm_size_min/1000000.," max",mpl_comm_size_max/1000000.
      write(*,'(a35,f12.6,a14,f12.6,a5)') "   Data rate                    ", &
      & meas%mean_mpl_comm_size/mpl_exch_min/1000000.,"MB/s     mean ",meas%mean_mpl_comm_size/meas%mean_mpl_exch/1000000.,"MB/s"
      end if
#endif
      write(*,*)
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "Communication duration          ", &
      & meas%mean_comm," s", meas%mean_comm/meas%mean_total*100.,"%    (", & 
      & comm_min,comm_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "Set boundary  duration          ",  &
& meas%mean_sbnd," s", meas%mean_sbnd/meas%mean_total*100.,"%    (",sbnd_min,sbnd_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a2)') "Calc fEq      duration          ",  &
& meas%mean_cceq," s", meas%mean_cceq/meas%mean_total*100.,"%"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "Collision     duration          ",  &
& meas%mean_coll," s", meas%mean_coll/meas%mean_total*100.,"%    (",coll_min,coll_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a7,2e9.2,a2)') "Propagation   duration          ",  &
& meas%mean_prop," s", meas%mean_prop/meas%mean_total*100.,"%    (",prop_min,prop_max,")"
      write(*,'(a35,f8.3,a2,f8.3,a2)') "Bounceback    duration          ",  &
& meas%mean_bncb," s", meas%mean_bncb/meas%mean_total*100.,"%"
      write(*,'(a35,f8.3,a2,f8.3,a2)') "Calcmacr vals duration          ",  &
& meas%mean_ccmv," s", meas%mean_ccmv/meas%mean_total*100.,"%"
         write(*,*)
      write(*,'(a35,f8.3,a6)') "Performance [MLUPs]             ",meas%mlups," MLUPs"

      call write_performance_results(result_state,s_par,prc,meas)

#ifdef VERBOSE
      write(*,*)
      write(*,*) "Timings of all processes"
      write(*,*)
      write(*,*) "Rk,     set_bnd,     coll,       comm,       prop,       bncb"
#endif /* VERBOSE */
   end if

#ifdef VERBOSE
   call mpl_barrier()
   write(*,'(i4,5f12.8)') prc%rk, meas%sbnd_duration, meas%coll_duration, meas%comm_duration, meas%prop_duration, meas%bncb_duration
#endif /* VERBOSE */

   end subroutine write_status 






subroutine alloc_mem(lb_dom,s_par,prc)

   type(lb_block)        :: lb_dom
   type(mpl_var)         :: prc
   type(sim_parameter)   :: s_par
   integer               :: alloc_stat,lx,ly,lz
#ifdef USE_CAF
#ifdef CHECK_CAF
   real   ,dimension(:,:,:)[:,:,:],allocatable :: caf_test
#endif
#endif

   lx = lb_dom%lx(1) 
   ly = lb_dom%lx(2) 
   lz = lb_dom%lx(3) 

#ifdef INIT_WITH_ROOT
   if(prc%rk == 0) then
      allocate(lb_dom%gu(LB_NODE(NDIM,s_par%gx(1),s_par%gx(2),s_par%gx(3))),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%gu"
      allocate(lb_dom%gfIn(LB_NODE(nnod,s_par%gx(1),s_par%gx(2),s_par%gx(3))),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%gfIn"
      allocate(lb_dom%grho(s_par%gx(1),s_par%gx(2),s_par%gx(3)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%grho"
      allocate(lb_dom%gstate(s_par%gx(1),s_par%gx(2),s_par%gx(3)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%g_state"
      allocate(s_par%g_x(s_par%gx(1)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%g_x"
      allocate(s_par%g_y(s_par%gx(2)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%g_y"
      allocate(s_par%g_z(s_par%gx(3)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%g_z"
#ifdef SPONGE
      allocate(lb_dom%gomega(s_par%gx(1),s_par%gx(2),s_par%gx(3)))
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%gsponge"
#endif /* SPONGE */
   else 
      allocate(lb_dom%gu(LB_NODE(1,1,1,1)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%gu"
      allocate(lb_dom%gfIn(LB_NODE(1,1,1,1)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%gfIn"
      allocate(lb_dom%grho(1,1,1),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%grho"
      allocate(lb_dom%gstate(1,1,1),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%g_state"
#ifdef SPONGE
      allocate(lb_dom%gomega(1,1,1))
#endif /* SPONGE */
      allocate( s_par%g_x(s_par%gx(1)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%g_x"
      allocate( s_par%g_y(s_par%gx(2)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%g_y"
      allocate( s_par%g_z(s_par%gx(3)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%g_z"
   end if
#endif

   if(prc%rk==0) call calc_mem(prc,lb_dom,s_par)
   allocate(lb_dom%fOut(LB_NODE(nnod,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%fOut"

!   if(associated(lb_dom%fOut)) write(*,*) "fOut OK"

if(s_par%tInitial > 0) then
   allocate(lb_dom%u0(  LB_NODE(NDIM,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%u0 "
endif

   allocate(lb_dom%u(  LB_NODE(NDIM,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%u  "
!   if(allocated(lb_dom%u)) write(*,*) "u OK"
   allocate(lb_dom%rho(0:(lx+1),0:(ly+1),0:(lz+1)),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%rho"
   allocate(lb_dom%state(lx,ly,lz),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%state"
#ifdef SPONGE
   allocate(lb_dom%omega(lx,ly,lz),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%omega"
#endif /* SPONGE */
   allocate(lb_dom%x(lx),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%x  "
   allocate(lb_dom%y(ly),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%y  "
   allocate(lb_dom%z(lz),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%z  "
   allocate(lb_dom%fIn(LB_NODE(nnod,0:(lx+1),0:(ly+1),0:(lz+1))),stat=alloc_stat)
      if(alloc_stat /=0) write(*,*) prc%rk,"Error allocating lb_dom%fIn"
end subroutine alloc_mem 

end module function_set
!------------------------------------------------------------------------


