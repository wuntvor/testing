! ------------------------------------------------------------------------
!
!         oooo   .o8                       
!          888   888                       
!          888   888oooo.   .ooooo.        
!          888   d88   88b d88    Y8       
!          888   888   888 888             
!          888   888   888 888   .o8       
!         o888o  `Y8bod8P  `Y8bod8P       
!
! This is the Lattice Boltzmann Kernel LBC for testing purposes only.
! Loosely based on J. Bernsdorfs anb-Lattice Boltzmann-Code
! 2D Lattice Boltzmann calculations with boundary and obstacle implementation.
!
!------------------------------------------------------------------------
!
!
! - MRT Model parameter was tuned, so that bulk viscosity is set equal to bgk. Watch parameter s1 s_mrt(2) 
!
!========================================================================

PROGRAM lbmain
   use mpl_lib
   use lbmodel
   use lbm_functions
   use lb_bc
   use mpl_set
   use function_set
   use tools
   use timing

#include "include/replace.h"
   implicit none
#ifdef USE_CAF
   type(measure)        :: meas[*]
#else
   type(measure)        :: meas
#endif
   type(lb_block)       :: lb_dom
   type(sim_parameter)  :: s_par 
   type(mpl_var)        :: prc
   real(R8B)            :: total_start,total_end
   integer              :: tStep,countOut,result_state
   integer              :: i,j,k
   integer,pointer      :: tstep_cnt 
   logical              :: stop_request    
   real(R8B),pointer,dimension(:,:,:,:) :: pnt_cur,pnt_nxt

   countOut             = 0
   result_state         = 0
   tstep_cnt => gtstep_cur

   stop_request         = .false.
  
!------------------------------------------------------------------------
! Initialize environment

   call mpl_init(meas,prc)
   call cpu_time_measure(total_start)
   call cpu_time_measure(meas%tStart)
   call read_params(s_par,'lbc.params',prc)

   call mpl_init_domain(lb_dom,s_par,prc)

   call alloc_mem(lb_dom,s_par,prc)
#ifdef USE_ADCL
   call mpl_adcl_init(lb_dom,s_par,prc)
#endif /*USE_ADCL*/


!------------------------------------------------------------------------
! Initialize data and distribute to processes
!
#ifdef D2Q9
   lb_dom%fIn(NDX(1,:,:,:))=0.4444444444444d0
   lb_dom%fIn(NDX(2:5,:,:,:))=0.111111111111d0
   lb_dom%fIn(NDX(6:9,:,:,:))=0.02777777777777d0
#else
   lb_dom%fIn(NDX(1,:,:,:))=0.33333333333333d0
   lb_dom%fIn(NDX(2:7,:,:,:))=0.0555555555555d0
   lb_dom%fIn(NDX(8:19,:,:,:))=0.02777777777777d0
#endif
!   lb_dom%fIn  = -10.d0
   lb_dom%fOut = lb_dom%fIn


   !---------------------------
   ! Master thread initializes initial field and geometry

#ifdef INIT_WITH_ROOT
   if(prc%rk==0) then
      call get_geo(lb_dom,s_par,prc)
      call init_field(lb_dom%gfIn,lb_dom%gu,lb_dom%grho,lb_dom%gstate,s_par,lb_dom,prc)
      call calc_fEq_global(lb_dom%grho,lb_dom%gu,lb_dom%gfIn,s_par) ! reset complete domain to eq
   end if

   !---------------------------
   ! Decompose domain 
   call mpl_decompose_domain(lb_dom,s_par,prc)
   
#else
   call get_geo_each(lb_dom,s_par,prc)
   call init_field_each(lb_dom,s_par,prc)
   call calc_fEq(lb_dom,lb_dom%rho,lb_dom%u,lb_dom%fIn,0,0,0,0,0,0,s_par)
#endif
   !---------------------------
   ! Initialize the obstacle and boundary vectors 
   call init_vector(lb_dom,s_par,prc)

   call ramp_u(s_par,0)

   call set_bnd(lb_dom,lb_dom%state,lb_dom%fIn,lb_dom%rho,lb_dom%u,s_par,prc,meas,0)


   call cpu_time_measure(meas%tEnd)
   meas%init_duration = real(meas%tEnd,8) - real(meas%tStart,8)
!------------------------------------------------------------------------
! Main Loop 
!
   if(prc%rk==0) then  
     write(*,*); write(*,*) "  --------------------------------------"
     write(*,*) "  Starting main time loop"; write(*,*)
   endif

   tstep_cnt = 0 
   do tStep = 1,s_par%tMax
   pnt_nxt => lb_dom%fIn

      call cpu_time_measure(meas%tStart)
      !---------------------------
      ! Set Boundaries on fIn  
      call treat_initial(lb_dom,s_par,tstep,prc)

      tstep_cnt = tstep_cnt+1

      call ramp_u(s_par,tstep)
#ifdef COMBINE_STREAM_COLLIDE
      if(mod(tStep,2)==0) then
      pnt_nxt => lb_dom%fIn
      pnt_cur => lb_dom%fOut
#ifdef USE_ADCL
      prc%adcl_active_req=1
#endif /*USE_ADCL*/
      else
      pnt_cur => lb_dom%fIn
      pnt_nxt => lb_dom%fOut
#ifdef USE_ADCL
      prc%adcl_active_req=0
#endif /*USE_ADCL*/
      endif
#else
      pnt_cur => lb_dom%fIn
#ifdef USE_ADCL
      prc%adcl_active_req=0
#endif /*USE_ADCL*/
      pnt_nxt => lb_dom%fOut
#endif

      !---------------------------
      ! Exchange between domains 
      if(prc%size > 1) then
         call mpl_exchange(lb_dom,pnt_cur,prc,meas)
      end if
!if(s_par%initial ) then
!#ifndef COMBINE_STREAM_COLLIDE
!      call propagate(lb_dom,pnt_nxt,pnt_cur,lb_dom%state,meas)
!
!      pnt_nxt => lb_dom%fIn
!      pnt_cur => lb_dom%fOut
!      call bounceback(lb_dom,pnt_nxt,pnt_cur,meas)
!#endif

!else
#ifdef COMBINE_STREAM_COLLIDE
      !---------------------------
      ! Combined Stream-collide Routine
#ifdef TRT
      call stream_collide_trt(lb_dom,pnt_nxt,pnt_cur,s_par%omega)
#else
      call stream_collide_bgk( lb_dom,pnt_nxt,pnt_cur,s_par,meas)
#endif
#else
      !---------------------------
      ! Propagation fIn => fOut 
      call propagate(lb_dom,pnt_nxt,pnt_cur,lb_dom%state,meas)

      ! pnt_nxt : fOut: post-propagation, pnt_cur: fIn: post-collision

      call set_nrbc_antibb(lb_dom,pnt_nxt,pnt_cur,s_par,prc,meas)
      call set_bnd(lb_dom,lb_dom%state,pnt_nxt,lb_dom%rho,lb_dom%u,s_par,prc,meas,tStep)

      pnt_nxt => lb_dom%fIn
      pnt_cur => lb_dom%fOut

      !---------------------------
      ! Collision routine fOut => fIn, excluding solid nodes 
#ifdef TRT
      call collide_trt(lb_dom,pnt_nxt,pnt_cur,s_par%omega)
#else
#ifdef MRT
      call collide_mrt( lb_dom,pnt_nxt,pnt_cur,s_par,meas)
#else
      call collide_bgk( lb_dom,pnt_nxt,pnt_cur,s_par,meas)
#endif
#endif
#endif

!endif

      call cpu_time_measure(meas%tEnd)
      meas%duration = meas%duration + real(meas%tEnd,8) - real(meas%tStart,8)

      !---------------------------
      ! Output and error check 
      if(MOD(tStep,s_par%tOut)==0 .and. tStep >= s_par%tOutBeg) then
         ! calc macroscopic values (is also done in do_kinetic)
         call calc_macr_vals(lb_dom,lb_dom%u,lb_dom%rho,pnt_nxt,meas)
         countOut=countOut+1
         call calc_rho_ges(lb_dom,pnt_nxt,s_par%rhoges,tStep,prc)
         if(tStep > s_par%tMax-s_par%tOut  ) s_par%goend=.true. 
         if(result_state == -1) s_par%goend=.true. 
         if(stop_request .eqv. .true.) s_par%goend=.true. 
      
         inquire(FILE='stop',exist=stop_request)
         if(s_par%save_output) then
            call gather_output(lb_dom,countOut,tStep,s_par,prc)
         endif

         ! Check for errors. If errors true, then exit after next iteration
         if(result_state == 0) then
            if(s_par%initial .eqv. .false.) &
            call check_density(lb_dom,tStep,result_state,prc%rk)
            call comm_res(result_state,prc)
         elseif(result_state ==-1 .and. s_par%goend .eqv. .true.) then
            write(*,*) "Result is errorous. Exiting..."
            exit 
         endif
     
         if(stop_request .and. s_par%goend ) then
            write(*,*) "Stop has been requested by stopfile 'Stop'"
            exit
         endif

      endif

#ifdef WRITE_TSTATUS
      if(MOD(tStep,WRITE_TSTATUS)==0) then
         ! Output duration of current timestep
         call write_tStat(tStep,meas%tEnd-meas%tStart,meas%duration,lb_dom,s_par,prc)
      endif
#endif
   end do

!------------------------------------------------------------------------
! Done with Main loop 

   if(prc%rk==0) then  
     write(*,*); write(*,*) "  --------------------------------------"
     write(*,*) "  "; write(*,*)
   endif
!------------------------------------------------------------------------
! Finish up 
!
   call cpu_time_measure(total_end)
   meas%total_duration = total_end - total_start

   call write_status(meas,result_state,tstep_cnt,s_par,prc)

   call dealloc_mem(lb_dom,s_par)

   call mpl_finish(prc)





END PROGRAM lbmain
!========================================================================





