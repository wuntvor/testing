module tools 
   use nrt_lib
   use lbmodel
#include "include/replace.h"
contains




   subroutine pdf_to_macro(fpdf,rho,u)
   !------------------------------------------------------------------------
   !
   ! calculate the macroscopic values from the microscopic distributions 
   !
     implicit none
     type(measure )  :: meas   
     real(R8B)       :: fpdf(nnod)
     real(R8B)       :: rho,u(3) 

     intent(in)      :: fpdf
     intent(out)     :: rho,u 

#ifdef D2Q9
               rho =  fpdf(I__0) + fpdf(I__E)  + &
                      fpdf(I__W) + fpdf(I__N)  + &
                      fpdf(I__S) + fpdf(I_NE)  + &
                      fpdf(I_SW) + fpdf(I_SE)  + & 
                      fpdf(I_NW) 
               u(1) = fpdf(I__E) - fpdf(I__W)  + &
                      fpdf(I_NE) - fpdf(I_SW)  + &
                      fpdf(I_SE) - fpdf(I_NW)     
               u(2) = fpdf(I__N) - fpdf(I__S)  + &
                      fpdf(I_NE) - fpdf(I_SW)  + &
                      fpdf(I_NW) - fpdf(I_SE)     
               u(3) = 0.d0
#endif
#ifdef D3Q19
               rho =  fpdf(I__0) + &
                      fpdf(I__E) + fpdf(I__W) + &
                      fpdf(I__N) + fpdf(I__S) + &
                      fpdf(I__T) + fpdf(I__B) + &
                      fpdf(I_NE) + fpdf(I_SW) + &
                      fpdf(I_SE) + fpdf(I_NW) + &
                      fpdf(I_TE) + fpdf(I_BW) + &
                      fpdf(I_BE) + fpdf(I_TW) + &
                      fpdf(I_TN) + fpdf(I_BS) + &
                      fpdf(I_BN) + fpdf(I_TS)
               u(1) = fpdf(I__E) - fpdf(I__W) + &
                      fpdf(I_NE) - fpdf(I_SW) + &
                      fpdf(I_SE) - fpdf(I_NW) + &
                      fpdf(I_TE) - fpdf(I_BW) + &
                      fpdf(I_BE) - fpdf(I_TW)  
               u(2) = fpdf(I__N) - fpdf(I__S) + &
                      fpdf(I_NE) - fpdf(I_SW) + &
                      fpdf(I_NW) - fpdf(I_SE) + & 
                      fpdf(I_TN) - fpdf(I_BS) + &
                      fpdf(I_BN) - fpdf(I_TS)
               u(3) = fpdf(I__T) - fpdf(I__B)  + &
                      fpdf(I_TE) - fpdf(I_BW)  + &
                      fpdf(I_TW) - fpdf(I_BE)  + &
                      fpdf(I_TN) - fpdf(I_BS)  + &
                      fpdf(I_TS) - fpdf(I_BN)
#endif
               u=u/rho
   end subroutine pdf_to_macro
   !------------------------------------------------------------------------



   function lb_feq(pos,rho,u,rho0)
   !------------------------------------------------------------------------
   !
   ! calculate the macroscopic values from the microscopic distributions 
   !
     implicit none
     real(R8B)       :: lb_feq    
     real(R8B)       :: rho,rho0,drho0,u(3),usq,ucx 
     integer         :: pos

     if(NDIM==2) u(3)=0.d0
     usq = u(1)*u(1) + u(2)*u(2) + u(3)*u(3)
     ucx = cx(pos,1)*u(1) + cx(pos,2)*u(2) + cx(pos,3)*u(3) 
#ifdef INCOMPRESSIBLE
     drho0 = rho0
#else
     drho0 = rho
#endif 

     lb_feq = t(pos)*(rho + drho0*(  &
     +  ucx*cs2inv & 
     +  ucx*ucx*cs2inv*cs2inv*0.5d0 & 
     -  usq*0.5d0*cs2inv))

     return

   end function lb_feq      
   !------------------------------------------------------------------------


#ifdef CALC_FEQ_SINGLE
   function lb_feq_all(rho,u,fEq,s_par)
   !------------------------------------------------------------------------
   !
   ! calculate the Maxwellian / equlibrium distribution function
   !

      type(sim_parameter)           :: s_par
      type(lb_block) :: lb_dom
      real(R8B)      :: usq,cs2inv,t2cs4inv,t2cs2inv,rholoc,rho0
      integer        :: i,j,k,l,err

      intent(in)     :: u,rho
      intent(inout)  :: fEq,lb_dom

      ! define constants 
      cs2inv  = 1._R8B/cs**2
      t2cs4inv = 1._R8B/(2*cs**4)
      t2cs2inv = 1._R8B/(2*cs**2)

#ifdef D2Q9    
             rholoc = rho
#ifdef INCOMPRESSIBLE
             rho0 = s_par%rho0
#else
             rho0 = rholoc
#endif 
             usq =  (u(NDX(1,i,j,k))**2 + u(NDX(2,i,j,k))**2)*t2cs2inv
             
             fEq(NDX(1,i,j,k)) = t(1)*(rholoc - rho0*usq)

             fEq(NDX(2,i,j,k)) = t(2)*(rholoc + rho0*((u(NDX(1,i,j,k)))*cs2inv &
            + (u(NDX(1,i,j,k)))**2*t2cs4inv   - usq))

             fEq(NDX(3,i,j,k)) = t(3)*(rholoc + rho0*((u(NDX(2,i,j,k)))*cs2inv &
            + (u(NDX(2,i,j,k)))**2*t2cs4inv   - usq))

             fEq(NDX(4,i,j,k)) = t(4)*(rholoc + rho0*((-u(NDX(1,i,j,k)))*cs2inv &
           + (-u(NDX(1,i,j,k)))**2*t2cs4inv  - usq))

             fEq(NDX(5,i,j,k)) = t(5)*(rholoc + rho0*((-u(NDX(2,i,j,k)))*cs2inv &
           + (-u(NDX(2,i,j,k)))**2*t2cs4inv   - usq))

             fEq(NDX(6,i,j,k)) = t(6)*(rholoc  &
            + rho0*((u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))*cs2inv &
            + (u(NDX(1,i,j,k))+u(NDX(2,i,j,k)))**2*t2cs4inv  - usq))

             fEq(NDX(7,i,j,k)) = t(7)*(rholoc  &
           + rho0*((-u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))*cs2inv &
           + (-u(NDX(1,i,j,k))+u(NDX(2,i,j,k)))**2*t2cs4inv  - usq))

             fEq(NDX(8,i,j,k)) = t(8)*(rholoc  &
           + rho0*((-u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))*cs2inv &
           + (-u(NDX(1,i,j,k))-u(NDX(2,i,j,k)))**2*t2cs4inv    - usq))

             fEq(NDX(9,i,j,k)) = t(9)*(rholoc  &
            + rho0*((u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))*cs2inv &
            + (u(NDX(1,i,j,k))-u(NDX(2,i,j,k)))**2*t2cs4inv  - usq))


#endif
#ifdef D3Q19

          usq =  (u(NDX(1,i,j,k))**2 + u(NDX(2,i,j,k))**2 + u(NDX(3,i,j,k))**2)*t2cs2inv
          rholoc   = rho(i,j,k) 
#ifdef INCOMPRESSIBLE
             rho0 = s_par%rho0
#else
             rho0 = rholoc
#endif 

           fEq(NDX(1,i,j,k)) = t(1)*(rholoc  &
                - rho0*usq) 
           fEq(NDX(2,i,j,k)) = t(2)*(rholoc + rho0*((u(NDX(1,i,j,k)))*cs2inv &
                + (u(NDX(1,i,j,k)))**2*t2cs4inv  &
                - usq))
           fEq(NDX(3,i,j,k)) = t(3)*(rholoc + rho0*((-u(NDX(1,i,j,k)))*cs2inv &
                + u(NDX(1,i,j,k))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(4,i,j,k)) = t(4)*(rholoc + rho0*((u(NDX(2,i,j,k)))*cs2inv &
                + (u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(5,i,j,k)) = t(5)*(rholoc + rho0*((-u(NDX(2,i,j,k)))*cs2inv &
                + (-u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(6,i,j,k)) = t(6)*(rholoc + rho0*((u(NDX(3,i,j,k)))*cs2inv &
                + (u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(7,i,j,k)) = t(7)*(rholoc + rho0*((-u(NDX(3,i,j,k)))*cs2inv &
                + (-u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(8,i,j,k)) = t(8)*(rholoc  &
                + rho0*((u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))*cs2inv &
                + (u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(9,i,j,k)) = t(9)*(rholoc  &
                + rho0*((-u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))*cs2inv &
                + (-u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(10,i,j,k)) = t(10)*(rholoc  &
                + rho0*((u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))*cs2inv &
                + (u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(11,i,j,k)) = t(11)*(rholoc  &
                + rho0*((-u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))*cs2inv &
                + (-u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(12,i,j,k)) = t(12)*(rholoc  &
                + rho0*((u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                + (u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(13,i,j,k)) = t(13)*(rholoc  &
                + rho0*((-u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                + (-u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(14,i,j,k)) = t(14)*(rholoc  &
                + rho0*((u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                + (u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(15,i,j,k)) = t(15)*(rholoc  &
                + rho0*((-u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                + (-u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(16,i,j,k)) = t(16)*(rholoc  &
                + rho0*((u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                + (u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(17,i,j,k)) = t(17)*(rholoc  &
                + rho0*((-u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                + (-u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(18,i,j,k)) = t(18)*(rholoc  &
                + rho0*((u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                + (u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
           fEq(NDX(19,i,j,k)) = t(19)*(rholoc  &
                + rho0*((-u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                + (-u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                - usq) )
               
#endif    
      return
  end subroutine lb_feq_all
  !------------------------------------------------------------------------
#endif /* calc_feq_single */







   subroutine check_density(lb_dom,tStep,result_state,rank)
   use nrt_lib
   use lbmodel
   !
   ! check densities if below zero or similar errors occur
   !
      implicit none
      type(lb_block),intent(inout)  :: lb_dom
      integer, intent(in)           :: tStep,rank
      integer                       :: i,j,k,l,count_error,result_state

      count_error = 0

      if(result_state == 0) then
         do k=1,lb_dom%lx(3)
            do j=1,lb_dom%lx(2)
               do i=1,lb_dom%lx(1)
               if(count_error > 10) exit
               if(btest(lb_dom%state(i,j,k),wall)    .eqv. .false.) then ! .and. &
               if(btest(lb_dom%state(i,j,k),nr_wall) .eqv. .false.) then 
               if(btest(lb_dom%state(i,j,k),nr_wall_in) .eqv. .false.) then 
                  do l= 1,nnod
                     if(lb_dom%fIn(NDX(l,i,j,k)) == 0.) then
                        write(*,*) rank,"zero at ",i,j,k,l,lb_dom%fIn(NDX(l,i,j,k))
                     endif
                     if(lb_dom%fIn(NDX(l,i,j,k)) .ne. -100.d0) then
                        if(lb_dom%fIn(NDX(l,i,j,k)) < 0 .or. &
                     &lb_dom%fIn(NDX(l,i,j,k)) .ne. lb_dom%fIn(NDX(l,i,j,k))) then
                           write(*,'(2i3,a20,4i4,a2,f10.4)') rank,tStep,&
                           &"Error at (i,j,k,l)",i,j,k,l,":",lb_dom%fIn(NDX(l,i,j,k))
                           count_error = count_error + 1
                           result_state = -1
                        endif
                     end if
                  end do
               endif
               endif
               endif
               end do
            end do
         end do
      end if

   end subroutine check_density

   subroutine calc_mem(prc,lb_dom,s_par)
      implicit none
      type(lb_block),intent(inout) :: lb_dom
      type(mpl_var) :: prc
      type(measure) :: meas
      type(sim_parameter) :: s_par
      integer  :: gx,gy,gz  
      gx = s_par%gx(1)
      gy = s_par%gx(2)
      gz = s_par%gx(3)
      meas%mem_master = 0.d0
#ifdef INIT_WITH_ROOT
      meas%mem_master = real(R8B*(int(gx)*int(gy)*int(gz)*int(1+NDIM+nnod) + &
            &  int(gx) + int(gy) + int(gz)) + int(I4B)*(int(gx)*int(gy)*int(gz)),R8B)/10**6
#endif
      meas%mem_thread = real(R8B*(int(lb_dom%lx(1)+2)*int(lb_dom%lx(2)+2)*int(2+lb_dom%lx(3))*int(1+NDIM+nnod) + & 
            & int(lb_dom%lx(1)+2) + int(lb_dom%lx(2)+2) + int(lb_dom%lx(3)+2)) + &
            & int(I4B)*(int(2+gx)*int(2+gy)*int(2+gz)),R8B)/10**6
#ifdef USE_ADCL
#ifdef D2Q9
      meas%mem_thread = meas%mem_thread + real((lb_dom%lx(1)+2)*(lb_dom%lx(2)+2)*R8B*nnod,R8B)/10**6
#endif
#ifdef D3Q19
      meas%mem_thread = meas%mem_thread + real((lb_dom%lx(1)+2)*(lb_dom%lx(2)+2)*(lb_dom%lx(3)+2)*R8B*nnod,R8B)/10**6
#endif
#endif
      write(*,'(a35,f28.3)') "Memory allocated per thread  [MB]",meas%mem_thread
      write(*,'(a35,f28.3)') "Additional Memory for master [MB]",meas%mem_master
      write(*,'(a35,f28.3)') "Total Memory needed          [MB]",meas%mem_thread*real(prc%size)+meas%mem_master

   end subroutine calc_mem 

   subroutine ramp_u(s_par,tstep)
      type(sim_parameter),intent(inout) :: s_par
      integer,intent(in)                :: tstep
      if(tstep == 1) then
         s_par%umaxRamp = s_par%umax 
      endif
      if(tstep <= s_par%tRamp) then
         s_par%umax = s_par%umaxRamp*real(tstep)/real(s_par%tRamp)
      endif
   end subroutine ramp_u




end module tools
