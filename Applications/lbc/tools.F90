module tools 
   use nrt_lib
   use lbmodel
#include "include/replace.h"
contains
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

#ifdef VECTOR
      logical :: printout
 
      count_error = 0

      if (result_state .eq. 0) then
         do k=2,lb_dom%lx(3)-1
            do j=2,lb_dom%lx(2)-1
               do i=2,lb_dom%lx(1)-1
                  do l= 1,nnod
                     if(lb_dom%fIn(LB_NODE(l,i,j,k)) == 0.) then
                        printout = .true. 
                     endif
                     if(lb_dom%fIn(LB_NODE(l,i,j,k)) .ne. -100.d0) then
                        if(lb_dom%fIn(LB_NODE(l,i,j,k)) < 0 .or. &
                     &lb_dom%fIn(LB_NODE(l,i,j,k)) .ne. lb_dom%fIn(LB_NODE(l,i,j,k))) then
                           result_state = -1
                        endif
                     end if
                  end do
               end do
            end do
         end do
      end if
      if ( printout .or. result_state .eq. -1 ) then 
        ! do the non-vectorizing anaylsis
#endif
      count_error = 0

      if(result_state == 0) then
         do k=1,lb_dom%lx(3)
            do j=1,lb_dom%lx(2)
               do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
!if(lb_dom%state(i,j,k) /= wall) then
                  do l= 1,nnod
                     if(lb_dom%fIn(LB_NODE(l,i,j,k)) == 0.) then
                        write(*,*) rank,"zero at ",i,j,k,l,lb_dom%fIn(LB_NODE(l,i,j,k))
                     endif
                     if(lb_dom%fIn(LB_NODE(l,i,j,k)) .ne. -100.d0) then
                        if(lb_dom%fIn(LB_NODE(l,i,j,k)) < 0 .or. &
                     &lb_dom%fIn(LB_NODE(l,i,j,k)) .ne. lb_dom%fIn(LB_NODE(l,i,j,k))) then
                           write(*,'(2i3,a20,4i4,a2,f10.4)') rank,tStep,&
                           &"Error at (i,j,k,l)",i,j,k,l,":",lb_dom%fIn(LB_NODE(l,i,j,k))
                           count_error = count_error + 1
                           result_state = -1
                        endif
                     end if
                  end do
endif
               end do
            end do
         end do
      end if
#ifdef VECTOR
   endif 
#endif
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
      if(tstep == 0) then
         s_par%umaxRamp = s_par%umax 
      endif
      if(tstep <= s_par%tRamp) then
         s_par%umax = s_par%umaxRamp*real(tstep)/real(s_par%tRamp)
      endif
   end subroutine ramp_u




end module tools
