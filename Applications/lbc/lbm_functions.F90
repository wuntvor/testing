module lbm_functions
   use mpl_lib
   use lbmodel
   use function_set
   use timing
   use lb_init
   use lb_bc
   use lb_geo 
#include "include/replace.h"
contains



#ifndef COMBINE_STREAM_COLLIDE
   subroutine propagate(lb_dom,fOut,fIn,state,meas)

!------------------------------------------------------------------------
! propagation step: move the densities to the neighbor nodes
! performs PULL of densities


   implicit none
   type(lb_block),intent(inout)  :: lb_dom
   type(measure )            :: meas  
   real(R8B), intent(in)     :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
   real(R8B), intent(inout)  :: fOut(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
   integer, intent(in)       :: state(lb_dom%lx(1),lb_dom%lx(2),lb_dom%lx(3))
   integer                   :: i,j,k,l,count_err
   real(R8B)                 :: ftmp(nnod)

      call cpu_time_measure(meas%tSt_comm) 
      count_err=0


      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
                  fOut(NDX(I__0,i,j,k))=fIn(NDX(I__0,i-cx(I__0,1),j-cx(I__0,2),k-cx(I__0,3)))
                  fOut(NDX(I__E,i,j,k))=fIn(NDX(I__E,i-cx(I__E,1),j-cx(I__E,2),k-cx(I__E,3)))
                  fOut(NDX(I__W,i,j,k))=fIn(NDX(I__W,i-cx(I__W,1),j-cx(I__W,2),k-cx(I__W,3)))
                  fOut(NDX(I__N,i,j,k))=fIn(NDX(I__N,i-cx(I__N,1),j-cx(I__N,2),k-cx(I__N,3)))
                  fOut(NDX(I__S,i,j,k))=fIn(NDX(I__S,i-cx(I__S,1),j-cx(I__S,2),k-cx(I__S,3)))
                  fOut(NDX(I_NE,i,j,k))=fIn(NDX(I_NE,i-cx(I_NE,1),j-cx(I_NE,2),k-cx(I_NE,3)))
                  fOut(NDX(I_SW,i,j,k))=fIn(NDX(I_SW,i-cx(I_SW,1),j-cx(I_SW,2),k-cx(I_SW,3)))
                  fOut(NDX(I_SE,i,j,k))=fIn(NDX(I_SE,i-cx(I_SE,1),j-cx(I_SE,2),k-cx(I_SE,3)))
                  fOut(NDX(I_NW,i,j,k))=fIn(NDX(I_NW,i-cx(I_NW,1),j-cx(I_NW,2),k-cx(I_NW,3)))
#ifdef D3Q19
                  fOut(NDX(I__T,i,j,k))=fIn(NDX(I__T,i-cx(I__T,1),j-cx(I__T,2),k-cx(I__T,3)))
                  fOut(NDX(I__B,i,j,k))=fIn(NDX(I__B,i-cx(I__B,1),j-cx(I__B,2),k-cx(I__B,3)))
                  fOut(NDX(I_TE,i,j,k))=fIn(NDX(I_TE,i-cx(I_TE,1),j-cx(I_TE,2),k-cx(I_TE,3)))
                  fOut(NDX(I_BW,i,j,k))=fIn(NDX(I_BW,i-cx(I_BW,1),j-cx(I_BW,2),k-cx(I_BW,3)))
                  fOut(NDX(I_BE,i,j,k))=fIn(NDX(I_BE,i-cx(I_BE,1),j-cx(I_BE,2),k-cx(I_BE,3)))
                  fOut(NDX(I_TW,i,j,k))=fIn(NDX(I_TW,i-cx(I_TW,1),j-cx(I_TW,2),k-cx(I_TW,3)))
                  fOut(NDX(I_TN,i,j,k))=fIn(NDX(I_TN,i-cx(I_TN,1),j-cx(I_TN,2),k-cx(I_TN,3)))
                  fOut(NDX(I_BS,i,j,k))=fIn(NDX(I_BS,i-cx(I_BS,1),j-cx(I_BS,2),k-cx(I_BS,3)))
                  fOut(NDX(I_BN,i,j,k))=fIn(NDX(I_BN,i-cx(I_BN,1),j-cx(I_BN,2),k-cx(I_BN,3)))
                  fOut(NDX(I_TS,i,j,k))=fIn(NDX(I_TS,i-cx(I_TS,1),j-cx(I_TS,2),k-cx(I_TS,3)))
#endif
   if(btest(lb_dom%state(i,j,k),wall) .eqv. .true.) then
                  ftmp = fOut(NDX(:,i,j,k)) 
                  fOut(NDX(opp(I__0),i,j,k))=ftmp(I__0)
                  fOut(NDX(opp(I__E),i,j,k))=ftmp(I__E)
                  fOut(NDX(opp(I__W),i,j,k))=ftmp(I__W)
                  fOut(NDX(opp(I__N),i,j,k))=ftmp(I__N)
                  fOut(NDX(opp(I__S),i,j,k))=ftmp(I__S)
                  fOut(NDX(opp(I_NE),i,j,k))=ftmp(I_NE)
                  fOut(NDX(opp(I_SW),i,j,k))=ftmp(I_SW)
                  fOut(NDX(opp(I_SE),i,j,k))=ftmp(I_SE)
                  fOut(NDX(opp(I_NW),i,j,k))=ftmp(I_NW)
#ifdef D3Q19
                  fOut(NDX(opp(I__T),i,j,k))=ftmp(I__T)
                  fOut(NDX(opp(I__B),i,j,k))=ftmp(I__B)
                  fOut(NDX(opp(I_TE),i,j,k))=ftmp(I_TE)
                  fOut(NDX(opp(I_BW),i,j,k))=ftmp(I_BW)
                  fOut(NDX(opp(I_BE),i,j,k))=ftmp(I_BE)
                  fOut(NDX(opp(I_TW),i,j,k))=ftmp(I_TW)
                  fOut(NDX(opp(I_TN),i,j,k))=ftmp(I_TN)
                  fOut(NDX(opp(I_BS),i,j,k))=ftmp(I_BS)
                  fOut(NDX(opp(I_BN),i,j,k))=ftmp(I_BN)
                  fOut(NDX(opp(I_TS),i,j,k))=ftmp(I_TS)
#endif
endif
            enddo
         enddo
      enddo

      call cpu_time_measure(meas%tEnd_comm)
      meas%prop_duration = meas%tEnd_comm - meas%tSt_comm + meas%prop_duration

   end subroutine propagate
  !------------------------------------------------------------------------
#endif /*COMBINE_STREAM_COLLIDE */






#ifdef COMBINE_STREAM_COLLIDE
#ifdef TRT

  subroutine stream_collide_trt(lb_dom,fOut,fIn,omega)


!------------------------------------------------------------------------
! Two-relaxation time collision routine



    implicit none
    type(lb_block),intent(inout)  :: lb_dom
    real(R8B) , intent(in)  :: omega
    real(R8B), intent(out)     :: fOut(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
    real(R8B), intent(inout)   :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))

    integer i_ct,i,j,k,lx(3),l
    real(R8B) ftmp(nnod) 

    ! relaxation parameter for the antisymmetric directions

    ! square speed of sound
    real(R8B), parameter :: c_squ = 1.d0 / 3.d0
    ! speed of sound related factor for equilibrium distribution functon
    real(R8B), parameter :: inv2csq2 = 1.d0 / (2.d0 * c_squ*c_squ)

#ifdef D3Q19
    ! common prefactors for calculating the equilibrium parts (D3Q19 model)
    real(R8B), parameter :: t0_0 = 1.d0/3.d0
    real(R8B), parameter :: t1x2_0 = 1.d0/18.d0 * 2.d0 ! twice the usual value due
    real(R8B), parameter :: t2x2_0 = 1.d0/36.d0 * 2.d0 ! to the way of calculations
#endif
#ifdef D2Q9
    ! common prefactors for calculating the equilibrium parts (D2Q9 model)
    real(R8B), parameter :: t0_0 = 4.d0/9.d0
    real(R8B), parameter :: t1x2_0 = 1.d0/9.d0 * 2.d0 ! twice the usual value due
    real(R8B), parameter :: t2x2_0 = 1.d0/36.d0 * 2.d0 ! to the way of calculations
#endif

    ! multiples of the relaxation parameter
    real(R8B) omega_h,asym_omega_h

    ! common part of the equilibrium distribution function for all directions
    real(R8B)  feq_common

    ! components of the velocity for calculating the (next two) eq. distributions
    real(R8B)  ui

    ! symmetric / antisymmetric part of the equilibrium distribution function
    real(R8B)  sym,asym

    ! in case of the default model, the actual values depend on density
    real(R8B)  t0,t1x2,t2x2,fac1,fac2

    ! local x,y,z velocities
    real(R8B)  u_x,u_y,u_z
    ! local density
    real(R8B)  loc_dens,inv_loc_dens

    ! inverse square velocity
    real(R8B) usq


    ! some global variables
    omega_h      = 0.5d0 * omega
    ! check that BGK is recovered; o.k. works fine
     asym_omega_h = 0.5d0 * omega
    ! magic Irina for fixed walls
    !asym_omega_h = 0.5d0 * 8.d0 * (2.d0-omega)/(8.d0-omega)
      lx = lb_dom%lx
#ifdef VECTOR
#include 'lbm_functions_vector.F90'
#else
    do k=1,lx(3)
      do j=1,lx(2)
         do i=1,lx(1)


            !--------------------------------------------
            ! Do Streaming process ( Pull values) 

                  ftmp(1)=fIn(NDX(1,i-cx(1,1),j-cx(1,2),k-cx(1,3)))
                  ftmp(2)=fIn(NDX(2,i-cx(2,1),j-cx(2,2),k-cx(2,3)))
                  ftmp(3)=fIn(NDX(3,i-cx(3,1),j-cx(3,2),k-cx(3,3)))
                  ftmp(4)=fIn(NDX(4,i-cx(4,1),j-cx(4,2),k-cx(4,3)))
                  ftmp(5)=fIn(NDX(5,i-cx(5,1),j-cx(5,2),k-cx(5,3)))
                  ftmp(6)=fIn(NDX(6,i-cx(6,1),j-cx(6,2),k-cx(6,3)))
                  ftmp(7)=fIn(NDX(7,i-cx(7,1),j-cx(7,2),k-cx(7,3)))
                  ftmp(8)=fIn(NDX(8,i-cx(8,1),j-cx(8,2),k-cx(8,3)))
                  ftmp(9)=fIn(NDX(9,i-cx(9,1),j-cx(9,2),k-cx(9,3)))
#ifdef D3Q19
                  ftmp(10)=fIn(NDX(10,i-cx(10,1),j-cx(10,2),k-cx(10,3)))
                  ftmp(11)=fIn(NDX(11,i-cx(11,1),j-cx(11,2),k-cx(11,3)))
                  ftmp(12)=fIn(NDX(12,i-cx(12,1),j-cx(12,2),k-cx(12,3)))
                  ftmp(13)=fIn(NDX(13,i-cx(13,1),j-cx(13,2),k-cx(13,3)))
                  ftmp(14)=fIn(NDX(14,i-cx(14,1),j-cx(14,2),k-cx(14,3)))
                  ftmp(15)=fIn(NDX(15,i-cx(15,1),j-cx(15,2),k-cx(15,3)))
                  ftmp(16)=fIn(NDX(16,i-cx(16,1),j-cx(16,2),k-cx(16,3)))
                  ftmp(17)=fIn(NDX(17,i-cx(17,1),j-cx(17,2),k-cx(17,3)))
                  ftmp(18)=fIn(NDX(18,i-cx(18,1),j-cx(18,2),k-cx(18,3)))
                  ftmp(19)=fIn(NDX(19,i-cx(19,1),j-cx(19,2),k-cx(19,3)))
#endif





               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
!            if(lb_dom%state(i,j,k) /= wall ) then


         !--------------------------------------------
         ! Calculate macroscopic variables 

#ifdef D3Q19
       ! local density
       loc_dens =  ftmp(1)  + ftmp(2)  + &
                   ftmp(3)  + ftmp(4)  + ftmp(5) + &
                   ftmp(6)  + ftmp(7)  + &
                   ftmp(8)  + ftmp(9)  + ftmp(10) + & 
                   ftmp(11) + ftmp(12) + &
                   ftmp(13) + ftmp(14) + ftmp(15) + & 
                   ftmp(16) + ftmp(17) + &
                   ftmp(18) + ftmp(19)

       ! local x-velocity
       u_x = (ftmp(2)  - ftmp(3)  + &
              ftmp(8)  - ftmp(9)  + &
              ftmp(10) - ftmp(11) + & 
              ftmp(12) - ftmp(13) + &
              ftmp(14) - ftmp(15)  )

       ! local y-velocity
       u_y = (ftmp(4)  - ftmp(5) + &
             ftmp(8)  - ftmp(9)  - ftmp(10) + & 
             ftmp(11) + & 
             ftmp(16) - ftmp(17) + &
             ftmp(18) - ftmp(19))

       ! local z-velocity
       u_z = (ftmp(6)  - ftmp(7)  + &
              ftmp(12) - fIn(NDX(13,i,j,k) )- ftmp(14)  + & 
              ftmp(15) + ftmp(16) - &
              ftmp(17) - ftmp(18) + &
              ftmp(19))
#endif
#ifdef D2Q9
       ! local density
       loc_dens =  ftmp(1)  + ftmp(2)  + &
                   ftmp(3)  + ftmp(4)  + ftmp(5) + &
                   ftmp(6)  + ftmp(7)  + &
                   ftmp(8)  + ftmp(9)  

       ! local x-velocity
       u_x =  ftmp(2)  - ftmp(4)  + &
              ftmp(6)  - ftmp(7)  + &
              ftmp(9)  - ftmp(8)  

       ! local y-velocity
       u_y = ftmp(3)  - ftmp(5)  + &
             ftmp(6)  + ftmp(7)  - & 
             ftmp(8)  - ftmp(9) 

       ! local z-velocity
       u_z = 0.0d0 
#endif
       ! inverse
       inv_loc_dens = 1.d0 / loc_dens

       ! transfer moments to velocities
       u_x = u_x * inv_loc_dens
       u_y = u_y * inv_loc_dens
       u_z = u_z * inv_loc_dens
       ! square velocity and derived constants
       usq  = u_x*u_x + u_y*u_y + u_z*u_z

       ! common part of F^eq (however depends on default/incompressible model)
       ! firt part and last part of the F^eq expression
       feq_common = 1.d0  - 1.5d0 * usq
       t0 = t0_0 * loc_dens

       !--------------------------------------------
       ! Relaxation Process 

       fOut(NDX(1,i,j,k)) = ftmp(1)*(1.d0-omega) + omega*t0*feq_common

#ifdef D3Q19
       t2x2 = t2x2_0 * loc_dens
       fac2 = t2x2 * inv2csq2

       ui   = u_x + u_y
       sym  = omega_h*(ftmp(8) + ftmp(9) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(8) - ftmp(9) - 3.d0*t2x2*ui )
       fOut(NDX(8,i,j,k))  = ftmp(8) - sym - asym
       fOut(NDX(9,i,j,k))  = ftmp(9) - sym + asym

       ui   = u_x - u_y
       sym  = omega_h*(ftmp(10) + ftmp(11) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(10) - ftmp(11) - 3.d0*t2x2*ui )
       fOut(NDX(10,i,j,k)) = ftmp(10) - sym - asym
       fOut(NDX(11,i,j,k)) = ftmp(11) - sym + asym

       ui   = u_x + u_z
       sym  = omega_h*(ftmp(12) + ftmp(13) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(12) - ftmp(13) - 3.d0*t2x2*ui )
       fOut(NDX(12,i,j,k)) = ftmp(12) - sym - asym
       fOut(NDX(13,i,j,k)) = ftmp(13) - sym + asym

       ui   = u_x - u_z
       sym  = omega_h*(ftmp(14) + ftmp(15) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(14) - ftmp(15) - 3.d0*t2x2*ui )
       fOut(NDX(14,i,j,k)) = ftmp(14) - sym - asym
       fOut(NDX(15,i,j,k)) = ftmp(15) - sym + asym

       ui   = u_y + u_z
       sym  = omega_h*(ftmp(16) + ftmp(17) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(16) - ftmp(17) - 3.d0*t2x2*ui )
       fOut(NDX(16,i,j,k)) = ftmp(16) - sym - asym
       fOut(NDX(17,i,j,k)) = ftmp(17) - sym + asym

       ui   = u_y - u_z
       sym  = omega_h*(ftmp(18) + ftmp(19) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(18) - ftmp(19) - 3.d0*t2x2*ui )
       fOut(NDX(18,i,j,k)) = ftmp(18) - sym - asym
       fOut(NDX(19,i,j,k)) = ftmp(19) - sym + asym

       t1x2 = t1x2_0 * loc_dens
       fac1 = t1x2 * inv2csq2

       ui   = u_y
       sym  = omega_h*(ftmp(4) + ftmp(5) - fac1*ui*ui - t1x2*feq_common)
       asym = asym_omega_h*( ftmp(4) - ftmp(5) - 3.d0*t1x2*ui )
       fOut(NDX(4,i,j,k)) = ftmp(4) - sym - asym
       fOut(NDX(5,i,j,k)) = ftmp(5) - sym + asym

       ui   = u_x
       sym  = omega_h*(ftmp(2) + ftmp(3) - fac1*ui*ui - t1x2*feq_common)
       asym = asym_omega_h*( ftmp(2) - ftmp(3) - 3.d0*t1x2*ui )
       fOut(NDX(2,i,j,k)) = ftmp(2) - sym - asym
       fOut(NDX(3,i,j,k)) = ftmp(3) - sym + asym

       ui   = u_z
       sym  = omega_h*(ftmp(6) + ftmp(7)  - fac1*ui*ui - t1x2*feq_common)
       asym = asym_omega_h*( ftmp(6) - ftmp(7) - 3.d0*t1x2*ui )
       fOut(NDX(6,i,j,k)) = ftmp(6) - sym - asym
       fOut(NDX(7,i,j,k)) = ftmp(7) - sym + asym
#endif
#ifdef D2Q9
       t1x2 = t1x2_0 * loc_dens
       fac1 = t1x2 * inv2csq2

       ! y-direction densities
       ui   = u_y
       sym  = omega_h*(ftmp(3) + ftmp(5) - fac1*ui*ui - t1x2*feq_common)
       asym = asym_omega_h*( ftmp(3) - ftmp(5) - 3.d0*t1x2*ui )
       fOut(NDX(3,i,j,k)) = ftmp(3) - sym - asym
       fOut(NDX(5,i,j,k)) = ftmp(5) - sym + asym

       ! x-direction densities
       ui   = u_x
       sym  = omega_h*(ftmp(2) + ftmp(4) - fac1*ui*ui - t1x2*feq_common)
       asym = asym_omega_h*( ftmp(2) - ftmp(4) - 3.d0*t1x2*ui )
       fOut(NDX(2,i,j,k)) = ftmp(2) - sym - asym
       fOut(NDX(4,i,j,k)) = ftmp(4) - sym + asym

       ! xy-direction densities
       t2x2 = t2x2_0 * loc_dens
       fac2 = t2x2 * inv2csq2

       ui   = u_x + u_y
       sym  = omega_h*(ftmp(6) + ftmp(8) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(6) - ftmp(8) - 3.d0*t2x2*ui )
       fOut(NDX(6,i,j,k))  = ftmp(6) - sym - asym
       fOut(NDX(8,i,j,k))  = ftmp(8) - sym + asym

       ui   = u_x - u_y
       sym  = omega_h*(ftmp(9) + ftmp(7) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(9) - ftmp(7) - 3.d0*t2x2*ui )
       fOut(NDX(9,i,j,k)) = ftmp(9) - sym - asym
       fOut(NDX(7,i,j,k)) = ftmp(7) - sym + asym
#endif
      else

         !--------------------------------------------
         ! Do Bounce Back 

         fOut(NDX(1,i,j,k)) = ftmp(opp(1))
         fOut(NDX(2,i,j,k)) = ftmp(opp(2))
         fOut(NDX(3,i,j,k)) = ftmp(opp(3))
         fOut(NDX(4,i,j,k)) = ftmp(opp(4))
         fOut(NDX(5,i,j,k)) = ftmp(opp(5))
         fOut(NDX(6,i,j,k)) = ftmp(opp(6))
         fOut(NDX(7,i,j,k)) = ftmp(opp(7))
         fOut(NDX(8,i,j,k)) = ftmp(opp(8))
         fOut(NDX(9,i,j,k)) = ftmp(opp(9))
#ifdef D3Q19
         fOut(NDX(10,i,j,k)) = ftmp(opp(10))
         fOut(NDX(11,i,j,k)) = ftmp(opp(11))
         fOut(NDX(12,i,j,k)) = ftmp(opp(12))
         fOut(NDX(13,i,j,k)) = ftmp(opp(13))
         fOut(NDX(14,i,j,k)) = ftmp(opp(14))
         fOut(NDX(15,i,j,k)) = ftmp(opp(15))
         fOut(NDX(16,i,j,k)) = ftmp(opp(16))
         fOut(NDX(17,i,j,k)) = ftmp(opp(17))
         fOut(NDX(18,i,j,k)) = ftmp(opp(18))
         fOut(NDX(19,i,j,k)) = ftmp(opp(19))
#endif

      endif
          end do
       end do
    end do
#endif /* VECTOR */

end subroutine stream_collide_trt
#endif /* TRT */
#endif /* COMBINE_STREAM_COLLIDE */



#ifdef MRT
#ifdef D3Q19
   subroutine collide_mrt_safe(lb_dom,fOut,fIn,s_par,meas)

!------------------------------------------------------------------------
! This is the MRT collision step
! This routine is intended for tests. It runs very slow, but is stable.
! Tested with lid-driven cavity on 2010-05-20.
! higher stability than bgk.
! Max achieved at lx=10**3 was u=0.03, om=1.95 (bgk crashed)
! 

      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure)              :: meas  
      real(R8B)                  :: usq,omega
      real(R8B)                  :: meq(nnod),mneq(nnod),mout(nnod),s_mrt(nnod)
      real(R8B)                  :: M_tr(nnod,nnod),M_inv(nnod,nnod),fOut_tmp(nnod)
      real(R8B)                  :: fOut_2d(nnod)
      real(R8B)                  :: rhot,ut(3)
      real(R8B)                  :: weps,wepsj,wxx,rholoc_inv 
      real(R8B)                  :: ux(3),rholoc,fEq(nnod),t2cs4inv,t2cs2inv,jx(3)
      real(R8B)                  :: rho0,inv_rho0 
      real(R8B), intent(out)     :: fOut(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)   :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer                    :: i,j,k,l,ll
#ifdef MRT_MASS_MOMENTUM 
      real(R8B)  :: check_rho,check_jx,check_jy,check_jz 
#endif /* LES_SMAGORINSKY */

M_tr(:, 1) = (/1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.         /)
M_tr(:, 2) = (/-1.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.        /)
M_tr(:, 3) = (/1.,-2.,-2.,-2.,-2.,-2.,-2.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.   /)
M_tr(:, 4) = (/0.,1.,-1.,0.,0.,0.,0.,1.,-1.,1.,-1.,1.,-1.,1.,-1.,0.,0.,0.,0.    /)
M_tr(:, 5) = (/0.,-2.,2.,0.,0.,0.,0.,1.,-1.,1.,-1.,1.,-1.,1.,-1.,0.,0.,0.,0.    /)
M_tr(:, 6) = (/0.,0.,0.,1.,-1.,0.,0.,1.,-1.,-1.,1.,0.,0.,0.,0.,1.,-1.,1.,-1.    /)
M_tr(:, 7) = (/0.,0.,0.,-2.,2.,0.,0.,1.,-1.,-1.,1.,0.,0.,0.,0.,1.,-1.,1.,-1.    /)
M_tr(:, 8) = (/0.,0.,0.,0.,0.,1.,-1.,0.,0.,0.,0.,1.,-1.,-1.,1.,1.,-1.,-1.,1.    /)
M_tr(:, 9) = (/0.,0.,0.,0.,0.,-2.,2.,0.,0.,0.,0.,1.,-1.,-1.,1.,1.,-1.,-1.,1.    /)
M_tr(:,10) = (/0.,2.,2.,-1.,-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,-2.,-2.,-2.,-2. /)
M_tr(:,11) = (/0.,-2.,-2.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,-2.,-2.,-2.,-2.   /)
M_tr(:,12) = (/0.,0.,0.,1.,1.,-1.,-1.,1.,1.,1.,1.,-1.,-1.,-1.,-1.,0.,0.,0.,0.   /)
M_tr(:,13) = (/0.,0.,0.,-1.,-1.,1.,1.,1.,1.,1.,1.,-1.,-1.,-1.,-1.,0.,0.,0.,0.   /)
M_tr(:,14) = (/0.,0.,0.,0.,0.,0.,0.,1.,1.,-1.,-1.,0.,0.,0.,0.,0.,0.,0.,0.       /)
M_tr(:,15) = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,-1.,-1.       /)
M_tr(:,16) = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,-1.,-1.,0.,0.,0.,0.       /)
M_tr(:,17) = (/0.,0.,0.,0.,0.,0.,0.,1.,-1.,1.,-1.,-1.,1.,-1.,1.,0.,0.,0.,0.     /)
M_tr(:,18) = (/0.,0.,0.,0.,0.,0.,0.,-1.,1.,1.,-1.,0.,0.,0.,0.,1.,-1.,1.,-1.     /)
M_tr(:,19) = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,-1.,-1.,1.,-1.,1.,1.,-1.     /)

#include "include/mrt_inv.incl"


      call cpu_time_measure(meas%tSt_comm) 

      omega = s_par%omega
        s_mrt(:) = 1.0_R8B
        s_mrt(2) = 2.0_R8B/(9.0_R8B * s_par%bulkViscosity + 1.0_R8B) !Bulk viscosity
        s_mrt(10)= omega
        s_mrt(12)= omega
        s_mrt(14:16)= omega


#ifdef CHECK_MRT
!     s_mrt(1:nnod) = omega 
#endif /* CHECK_MRT */


      ! Multiple Relaxation time Collision routine
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then

#ifdef D3Q19

#ifdef CHECK_MRT
! test data
fIn(NDX(1,i,j,k)) = 0.332333333503913_R8B
fIn(NDX(2,i,j,k)) = 0.060909350699940_R8B
fIn(NDX(3,i,j,k)) = 0.049868427134697_R8B
fIn(NDX(4,i,j,k)) = 0.101307711335621_R8B
fIn(NDX(5,i,j,k)) = 0.055555555555500_R8B
fIn(NDX(6,i,j,k)) = 0.060909350699940_R8B
fIn(NDX(7,i,j,k)) = 0.050368426964117_R8B
fIn(NDX(8,i,j,k)) = 0.018918438141908_R8B
fIn(NDX(9,i,j,k)) = 0.027777777777770_R8B
fIn(NDX(10,i,j,k)) =0.027777777777770_R8B
fIn(NDX(11,i,j,k)) =0.010882810381340_R8B
fIn(NDX(12,i,j,k)) =0.033464906155991_R8B
fIn(NDX(13,i,j,k)) =0.022923982420168_R8B
fIn(NDX(14,i,j,k)) =0.027694444458659_R8B
fIn(NDX(15,i,j,k)) =0.027694444458659_R8B
fIn(NDX(16,i,j,k)) =0.018818739533272_R8B
fIn(NDX(17,i,j,k)) =0.027777777777770_R8B
fIn(NDX(18,i,j,k)) =0.010900357542234_R8B
fIn(NDX(19,i,j,k)) =0.027777777777770_R8B

#endif /* CHECK_MRT */



call pdf_to_macro(fIn(NDX(:,i,j,k)),rholoc,jx) 
rho0 = rholoc
ux = jx
jx = jx*rho0

rholoc_inv = 1.d0/rholoc
inv_rho0 = 1.d0/rho0

#ifdef CHECK_MRT
write(78,*) "fIn   "
write(78,*) fIn(NDX(:,i,j,k))
write(78,*) "rholoc",rholoc-0.993661390097039_R8B
write(78,*) "ux    ",ux(1:3)
#endif          

#ifdef TEST_MRT
!if(gtstep_cur==0) then
  rholoc = 1.1d0 
  rho0   = 1.1d0
  ux(1) =  0.01
  ux(2) = -0.01
  ux(3) =  0.d0
!endif
#endif 

                  meq(1:nnod) =  0.0d0 

                  meq( 1) = rholoc 
! e_eq
                  meq( 2) =   rho0*((ux(1))*(ux(1)) + &
                          (ux(2))*(ux(2))+(ux(3))*(ux(3)))
! eps_eq
! qx_eq, qy_eq, qz_eq
                  meq( 4) = ux(1)*rho0 
                  meq( 6) = ux(2)*rho0 
                  meq( 8) = ux(3)*rho0 
                  meq(10) = rho0*(2.*ux(1)*ux(1) - (ux(2)*ux(2) + ux(3)*ux(3)))
!p_xy,p_yz,p_xz
                  meq(12) =  rho0*ux(1)*ux(2)
                  meq(14) =  rho0*ux(2)*ux(3)
                  meq(16) =  rho0*ux(1)*ux(3)


! do the back transformation once again

               do l=1,nnod
                 mneq(l) = 0.0d0
                  ! Berechnung der  nicht-gg-Werte im Momentenraum
                  do ll=1,nnod
                     mneq(l) = mneq(l) + fIn(NDX(ll,i,j,k))*M_tr(ll,l)
                  end do
               end do
#ifdef TEST_MRT
!  mneq = 0.0d0 
#endif 



               ! relaxation
               do l=1,nnod
                  !mout(l) = mneq(l) - s_mrt(l)*(mneq(l) - meq(l)) 
                  mout(l) = s_mrt(l)*(mneq(l) - meq(l)) 
                  !write(45,*) l, " mout",mout(l)
                  !write(45,*) l, " meq",meq(l)
               end do 
               ! Rücktransformation
               do l=1,nnod
                  fOut_2d(l)  = 0.0d0
                  ! Berechnung der  nicht-gg-Werte im Momentenraum
                  do ll=1,nnod
                     fOut_2d(l)  = fOut_2d(l)  + mout(ll)*M_inv(l,ll)
                  end do
         
               end do
                !write(45,*) l," mout(ll)",fOut_2d(:)
               fOut(NDX(:,i,j,k))= fIn(NDX(:,i,j,k)) - fOut_2d(:) 


#ifdef TEST_MRT
write(45,*) "mneq",mneq
write(45,*) "meq",meq

write(45,*) fOut(NDX(1,i,j,k)) - 0.346866666666666768
write(45,*) fOut(NDX(2,i,j,k)) - 0.0810944444444444135
write(45,*) fOut(NDX(3,i,j,k)) - 0.0444277777777778021
write(45,*) fOut(NDX(4,i,j,k)) - 0.0394777777777777575
write(45,*) fOut(NDX(5,i,j,k)) - 0.0761444444444444313
write(45,*) fOut(NDX(6,i,j,k)) - 0.0627611111111111147
write(45,*) fOut(NDX(7,i,j,k)) - 0.0627611111111111147
write(45,*) fOut(NDX(8,i,j,k)) - 0.0313805555555555574
write(45,*) fOut(NDX(9,i,j,k)) - 0.0313805555555555712
write(45,*) fOut(NDX(10,i,j,k)) - 0.0497138888888888839
write(45,*) fOut(NDX(11,i,j,k)) - 0.0130472222222222222
write(45,*) fOut(NDX(12,i,j,k)) - 0.0430222222222222325
write(45,*) fOut(NDX(13,i,j,k)) - 0.0246888888888888887
write(45,*) fOut(NDX(14,i,j,k)) - 0.0430222222222222256
write(45,*) fOut(NDX(15,i,j,k)) - 0.0246888888888888887
write(45,*) fOut(NDX(16,i,j,k)) - 0.0222138888888888837
write(45,*) fOut(NDX(17,i,j,k)) - 0.0405472222222222276
write(45,*) fOut(NDX(18,i,j,k)) - 0.0222138888888888907
write(45,*) fOut(NDX(18,i,j,k)) - 0.0405472222222222276

!write(45,*) fOut(NDX(:,i,j,k)) 
call pdf_to_macro(fOut(NDX(:,i,j,k)),rhot,ut)
write(*,*) "with",rhot,ut 
stop
#endif 

#endif /* D3Q19 */



#ifdef MRT_MASS_MOMENTUM
call pdf_to_macro(fOut(NDX(:,i,j,k)),rhot,ut)
        if(abs(check_rho) > 0.00001) then
            write(*,*) "error in rho",i,j,k," checkrho",rhot-rholoc
            stop
         endif
        if(abs(check_jx ) > 0.00001) then
            write(*,*) "error in jx ",i,j,k," check_jx",ut(1)-jx(1)
            stop
         endif
        if(abs(check_jy ) > 0.00001) then
            write(*,*) "error in jy ",i,j,k," checkrho",ut(2)-jx(2)
            stop
         endif
        if(abs(check_jz ) > 0.00001) then
            write(*,*) "error in jz ",i,j,k," checkrho",ut(3)-jx(3)
            stop
         endif

#endif
#ifdef CHECK_MRT
write(78,*) "mEq   "
write(78,*) meq(:)
write(78,*) 
write(78,*) "mnEq   "
write(78,*) mneq(:)
write(78,*) 
write(78,*) "mout   "
write(78,*) mout(:)
write(78,*) 
write(78,*) "fOut   "
write(78,*) fOut(NDX(:,i,j,k))
write(78,*) 
write(78,*) "Out-momentums rho,jx  "
write(78,*) rholoc,jx 


   write(*,*) "rho",rholoc
   write(*,*) "jx ",jx
   write(*,*) 


write(*,*) fOut(NDX( 1,i,j,k)) - 0.309663695742359302  
write(*,*) fOut(NDX( 2,i,j,k)) - 0.0616057383160554112 
write(*,*) fOut(NDX( 3,i,j,k)) - 0.0517332466260554055 
write(*,*) fOut(NDX( 4,i,j,k)) - 0.0515416601071256447 
write(*,*) fOut(NDX( 5,i,j,k)) - 0.0534878633704589587 
write(*,*) fOut(NDX( 6,i,j,k)) - 0.0612593490684856523 
write(*,*) fOut(NDX( 7,i,j,k)) - 0.0515926059118189811 
write(*,*) fOut(NDX( 8,i,j,k)) - 0.0290845851724469764 
write(*,*) fOut(NDX( 9,i,j,k)) - 0.0251214409591136410 
write(*,*) fOut(NDX(10,i,j,k)) - 0.0334252989253408447 
write(*,*) fOut(NDX(11,i,j,k)) - 0.0275159514486741780 
write(*,*) fOut(NDX(12,i,j,k)) - 0.0358162124613323821 
write(*,*) fOut(NDX(13,i,j,k)) - 0.0260465950379990471 
write(*,*) fOut(NDX(14,i,j,k)) - 0.0306048873878154472 
write(*,*) fOut(NDX(15,i,j,k)) - 0.0305020131211487834 
write(*,*) fOut(NDX(16,i,j,k)) - 0.0290115202110680925 
write(*,*) fOut(NDX(17,i,j,k)) - 0.0251512502644014278 
write(*,*) fOut(NDX(18,i,j,k)) - 0.0273455014291499598 
write(*,*) fOut(NDX(19,i,j,k)) - 0.0331519746391499628 

stop
#endif /* CHECK_MRT */
         else
             ! copy the boundary nodes!!
             fOut(NDX(:,i,j,k)) = fIn(NDX(:,i,j,k))
         endif
            end do
         end do
      end do
      
      call cpu_time_measure(meas%tEnd_comm)
#ifdef DEBUG_LES
      write(66,*) gtstep_cur,omega_min
#endif /* LES_SMAGORINSKY */


      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration


   end subroutine collide_mrt_safe
  !------------------------------------------------------------------------
#endif /* D3Q19 */












#ifdef D3Q19
   subroutine collide_mrt(lb_dom,fOut,fIn,s_par,meas)

!------------------------------------------------------------------------
! This is the BGK collision step
! used by default



      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure)              :: meas  
      real(R8B)                  :: usq,omega
      real(R8B)                  :: ux(3),rholoc,fEq(nnod),t2cs4inv,t2cs2inv,rholoc_inv,jx(3)
      real(R8B)                  :: rho0 ,inv_rho0 
      real(R8B), intent(out)     :: fOut(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)   :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer                    :: i,j,k,l
      ! MRT Variables
      real(R8B)                  :: meq(nnod),mneq(nnod),mout(nnod),s_mrt(nnod)
      real(R8B)                  :: M_tr(nnod,nnod),M_inv(nnod,nnod),fOut_tmp(nnod)
      real(R8B)                  :: fOut_2d(nnod)
      real(R8B)                  :: rhot,ut(3)
#ifdef MRT_MASS_MOMENTUM 
      real(R8B)  :: check_rho,check_jx,check_jy,check_jz 
#endif /* LES_SMAGORINSKY */


      call cpu_time_measure(meas%tSt_comm) 

#ifdef CHECK_MRT
      omega = 1.8 
#else /* CHECK_MRT */
      omega = s_par%omega
#endif/* CHECK_MRT */

        s_mrt(:) = 1.0_R8B
        s_mrt(2) = 2.0_R8B/(9.0_R8B * s_par%bulkViscosity + 1.0_R8B) !Bulk viscosity
        s_mrt(10)= omega
        s_mrt(12)= omega
        s_mrt(14:16)= omega


      ! Multiple Relaxation time Collision routine

      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then

#ifdef SPONGE
        omega = lb_dom%omega(i,j,k)
        s_mrt(:) = 1.0_R8B
        s_mrt(2) = 2.0_R8B/(9.0_R8B * s_par%bulkViscosity + 1.0_R8B) !Bulk viscosity
        s_mrt(10)= omega
        s_mrt(12)= omega
        s_mrt(14:16)= omega
#endif /* SPONGE */

! macroscopic variables
   rholoc   =  fIn(NDX(I__0,i,j,k)) + fIn(NDX(I__E,i,j,k)) + &
               fIn(NDX(I__W,i,j,k)) + fIn(NDX(I__N,i,j,k)) + &
               fIn(NDX(I__S,i,j,k)) + fIn(NDX(I__T,i,j,k)) + &
               fIn(NDX(I__B,i,j,k)) + fIn(NDX(I_NE,i,j,k)) + &
               fIn(NDX(I_SW,i,j,k)) + fIn(NDX(I_SE,i,j,k)) + & 
               fIn(NDX(I_NW,i,j,k)) + fIn(NDX(I_TE,i,j,k)) + &
               fIn(NDX(I_BW,i,j,k)) + fIn(NDX(I_BE,i,j,k)) + &
               fIn(NDX(I_TW,i,j,k)) + fIn(NDX(I_TN,i,j,k)) + &
               fIn(NDX(I_BS,i,j,k)) + fIn(NDX(I_BN,i,j,k)) + &
               fIn(NDX(I_TS,i,j,k))
#ifdef INCOMPRESSIBLE
rho0=s_par%rho0
#else
rho0=rholoc
#endif

      if(s_par%initial .eqv. .false.) then
   jx(1)    = fIn(NDX(I__E,i,j,k)) - fIn(NDX(I__W,i,j,k)) + &
              fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
              fIn(NDX(I_SE,i,j,k)) - fIn(NDX(I_NW,i,j,k)) + &
              fIn(NDX(I_TE,i,j,k)) - fIn(NDX(I_BW,i,j,k)) + &
              fIn(NDX(I_BE,i,j,k)) - fIn(NDX(I_TW,i,j,k))  
   jx(2)    = fIn(NDX(I__N,i,j,k)) - fIn(NDX(I__S,i,j,k)) + &
              fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
              fIn(NDX(I_NW,i,j,k)) - fIn(NDX(I_SE,i,j,k)) + & 
              fIn(NDX(I_TN,i,j,k)) - fIn(NDX(I_BS,i,j,k)) + &
              fIn(NDX(I_BN,i,j,k)) - fIn(NDX(I_TS,i,j,k))
   jx(3)    = fIn(NDX(I__T,i,j,k)) - fIn(NDX(I__B,i,j,k))  + &
              fIn(NDX(I_TE,i,j,k)) - fIn(NDX(I_BW,i,j,k))  + &
              fIn(NDX(I_TW,i,j,k)) - fIn(NDX(I_BE,i,j,k))  + &
              fIn(NDX(I_TN,i,j,k)) - fIn(NDX(I_BS,i,j,k))  + &
              fIn(NDX(I_TS,i,j,k)) - fIn(NDX(I_BN,i,j,k))

       else
         jx(1) = lb_dom%u0(NDX(1,i,j,k))*rho0
         jx(2) = lb_dom%u0(NDX(2,i,j,k))*rho0
         jx(3) = lb_dom%u0(NDX(3,i,j,k))*rho0
       endif

       rholoc_inv = 1.d0/rholoc
       inv_rho0 = 1.d0/rho0

       ux(:) = jx(:)/rho0

 
! equilibrium moments

                 meq(1:nnod) =  0.0d0 

                  meq( 1) = rholoc 
! e_eq
                  meq( 2) =   rho0*((ux(1))*(ux(1)) + &
                          (ux(2))*(ux(2))+(ux(3))*(ux(3)))
! eps_eq
! qx_eq, qy_eq, qz_eq
                  meq( 4) = ux(1)*rho0 
                  meq( 6) = ux(2)*rho0 
                  meq( 8) = ux(3)*rho0 
                  meq(10) = rho0*(2.*ux(1)*ux(1) - (ux(2)*ux(2) + ux(3)*ux(3)))
!p_xy,p_yz,p_xz
                  meq(12) =  rho0*ux(1)*ux(2)
                  meq(14) =  rho0*ux(2)*ux(3)
                  meq(16) =  rho0*ux(1)*ux(3)


! transformation

! rho
mneq(1) = rholoc
! e
mneq(2) = -fIn(NDX( 1,i,j,k)) + fIn(NDX( 8,i,j,k))+ fIn(NDX( 9,i,j,k))+ fIn(NDX(10,i,j,k)) &
         + fIn(NDX(11,i,j,k))+ fIn(NDX(12,i,j,k))+ fIn(NDX(13,i,j,k))+ fIn(NDX(14,i,j,k))  &
         + fIn(NDX(15,i,j,k))+ fIn(NDX(16,i,j,k))+ fIn(NDX(17,i,j,k))+ fIn(NDX(18,i,j,k))  &
         + fIn(NDX(19,i,j,k))
! epsilon
mneq(3) =  fIn(NDX( 1,i,j,k))- 2.d0*(fIn(NDX( 2,i,j,k))+ fIn(NDX( 3,i,j,k))+ fIn(NDX( 4,i,j,k)) &
         + fIn(NDX( 5,i,j,k))+ fIn(NDX( 6,i,j,k))+ fIn(NDX( 7,i,j,k)))+ fIn(NDX( 8,i,j,k))  &
         + fIn(NDX( 9,i,j,k))+ fIn(NDX(10,i,j,k)) &
         + fIn(NDX(11,i,j,k))+ fIn(NDX(12,i,j,k))+ fIn(NDX(13,i,j,k))+ fIn(NDX(14,i,j,k))  &
         + fIn(NDX(15,i,j,k))+ fIn(NDX(16,i,j,k))+ fIn(NDX(17,i,j,k))+ fIn(NDX(18,i,j,k))  &
         + fIn(NDX(19,i,j,k))
! jx
mneq(4) =                          fIn(NDX( 2,i,j,k))- fIn(NDX( 3,i,j,k))+ fIn(NDX( 8,i,j,k))  &
         - fIn(NDX( 9,i,j,k))+ fIn(NDX(10,i,j,k)) &
         - fIn(NDX(11,i,j,k))+ fIn(NDX(12,i,j,k))- fIn(NDX(13,i,j,k))+ fIn(NDX(14,i,j,k))  &
         - fIn(NDX(15,i,j,k))
! qx
mneq(5) =                    -2.0d0*( fIn(NDX( 2,i,j,k))- fIn(NDX( 3,i,j,k)))+fIn(NDX( 8,i,j,k))  &
         - fIn(NDX( 9,i,j,k))+ fIn(NDX(10,i,j,k)) &
         - fIn(NDX(11,i,j,k))+ fIn(NDX(12,i,j,k))- fIn(NDX(13,i,j,k))+ fIn(NDX(14,i,j,k))  &
         - fIn(NDX(15,i,j,k))
! jy
mneq(6) =                          fIn(NDX( 4,i,j,k))- fIn(NDX( 5,i,j,k))+ fIn(NDX( 8,i,j,k))  &
         - fIn(NDX( 9,i,j,k))- fIn(NDX(10,i,j,k)) &
         + fIn(NDX(11,i,j,k))+ fIn(NDX(16,i,j,k))- fIn(NDX(17,i,j,k))+ fIn(NDX(18,i,j,k))  &
         - fIn(NDX(19,i,j,k))
! qy
mneq(7) =                    -2.0d0*(fIn(NDX( 4,i,j,k))- fIn(NDX( 5,i,j,k)))+ fIn(NDX( 8,i,j,k))  &
         - fIn(NDX( 9,i,j,k))- fIn(NDX(10,i,j,k)) &
         + fIn(NDX(11,i,j,k))+ fIn(NDX(16,i,j,k))- fIn(NDX(17,i,j,k))+ fIn(NDX(18,i,j,k))  &
         - fIn(NDX(19,i,j,k))
! jz
mneq(8) =                          fIn(NDX( 6,i,j,k))- fIn(NDX( 7,i,j,k)) &
         + fIn(NDX(12,i,j,k))- fIn(NDX(13,i,j,k))- fIn(NDX(14,i,j,k))  + fIn(NDX(15,i,j,k))&
         + fIn(NDX(16,i,j,k))- fIn(NDX(17,i,j,k))- fIn(NDX(18,i,j,k))  + fIn(NDX(19,i,j,k))
! qz
mneq(9) =                    -2.0d0*(fIn(NDX( 6,i,j,k))- fIn(NDX( 7,i,j,k))) &
         + fIn(NDX(12,i,j,k))- fIn(NDX(13,i,j,k))- fIn(NDX(14,i,j,k))  + fIn(NDX(15,i,j,k))&
         + fIn(NDX(16,i,j,k))- fIn(NDX(17,i,j,k))- fIn(NDX(18,i,j,k))  + fIn(NDX(19,i,j,k))
! 3pxx
mneq(10) =                         2.d0*(fIn(NDX( 2,i,j,k))+ fIn(NDX( 3,i,j,k)))-(fIn(NDX( 4,i,j,k)) &
         + fIn(NDX( 5,i,j,k))+ fIn(NDX( 6,i,j,k))+ fIn(NDX( 7,i,j,k))) & 
         + fIn(NDX( 8,i,j,k))+ fIn(NDX( 9,i,j,k))+ fIn(NDX(10,i,j,k))+ fIn(NDX(11,i,j,k)) &
         + fIn(NDX(12,i,j,k))+ fIn(NDX(13,i,j,k))+ fIn(NDX(14,i,j,k))+ fIn(NDX(15,i,j,k)) &
   -2.0*( fIn(NDX(16,i,j,k))+ fIn(NDX(17,i,j,k))+ fIn(NDX(18,i,j,k))+ fIn(NDX(19,i,j,k)))
! 3 pixx
mneq(11) =                        -2.d0*(fIn(NDX( 2,i,j,k))+ fIn(NDX( 3,i,j,k)))+ fIn(NDX( 4,i,j,k)) &
         + fIn(NDX( 5,i,j,k))+ fIn(NDX( 6,i,j,k))+ fIn(NDX( 7,i,j,k))  & 
         + fIn(NDX( 8,i,j,k))+ fIn(NDX( 9,i,j,k))+ fIn(NDX(10,i,j,k))+ fIn(NDX(11,i,j,k)) &
         + fIn(NDX(12,i,j,k))+ fIn(NDX(13,i,j,k))+ fIn(NDX(14,i,j,k))+ fIn(NDX(15,i,j,k)) &
   -2.0*( fIn(NDX(16,i,j,k))+ fIn(NDX(17,i,j,k))+ fIn(NDX(18,i,j,k))+ fIn(NDX(19,i,j,k)))
! pww
mneq(12) =                                                                                fIn(NDX( 4,i,j,k)) &
         + fIn(NDX( 5,i,j,k))- fIn(NDX( 6,i,j,k))- fIn(NDX( 7,i,j,k))  & 
         + fIn(NDX( 8,i,j,k))+ fIn(NDX( 9,i,j,k))+ fIn(NDX(10,i,j,k))+ fIn(NDX(11,i,j,k)) &
         -(fIn(NDX(12,i,j,k))+ fIn(NDX(13,i,j,k))+ fIn(NDX(14,i,j,k))+ fIn(NDX(15,i,j,k)))
! piww
mneq(13) =                                                                              - fIn(NDX( 4,i,j,k)) &
         - fIn(NDX( 5,i,j,k))+ fIn(NDX( 6,i,j,k))+ fIn(NDX( 7,i,j,k))  & 
         + fIn(NDX( 8,i,j,k))+ fIn(NDX( 9,i,j,k))+ fIn(NDX(10,i,j,k))+ fIn(NDX(11,i,j,k)) &
         -(fIn(NDX(12,i,j,k))+ fIn(NDX(13,i,j,k))+ fIn(NDX(14,i,j,k))+ fIn(NDX(15,i,j,k)))
! pxy
mneq(14) =  fIn(NDX( 8,i,j,k))+ fIn(NDX( 9,i,j,k))- fIn(NDX(10,i,j,k))- fIn(NDX(11,i,j,k))
! pyz
mneq(15) =  fIn(NDX(16,i,j,k))+ fIn(NDX(17,i,j,k))- fIn(NDX(18,i,j,k))- fIn(NDX(19,i,j,k))
! pxz
mneq(16) =  fIn(NDX(12,i,j,k))+ fIn(NDX(13,i,j,k))- fIn(NDX(14,i,j,k))- fIn(NDX(15,i,j,k))
! mx
mneq(17) = fIn(NDX( 8,i,j,k))- fIn(NDX( 9,i,j,k))+ fIn(NDX(10,i,j,k))- fIn(NDX(11,i,j,k)) &
          -fIn(NDX(12,i,j,k))+ fIn(NDX(13,i,j,k))- fIn(NDX(14,i,j,k))+ fIn(NDX(15,i,j,k))
! my
mneq(18) =-fIn(NDX( 8,i,j,k))+ fIn(NDX( 9,i,j,k))+ fIn(NDX(10,i,j,k))- fIn(NDX(11,i,j,k)) &
          +fIn(NDX(16,i,j,k))- fIn(NDX(17,i,j,k))+ fIn(NDX(18,i,j,k))- fIn(NDX(19,i,j,k))
! mz
mneq(19) =+fIn(NDX(12,i,j,k))- fIn(NDX(13,i,j,k))- fIn(NDX(14,i,j,k))+ fIn(NDX(15,i,j,k)) &
          -fIn(NDX(16,i,j,k))+ fIn(NDX(17,i,j,k))+ fIn(NDX(18,i,j,k))- fIn(NDX(19,i,j,k))

!relaxation 


mout( 1) = s_mrt( 1)*(mneq( 1) - meq( 1))
mout( 2) = s_mrt( 2)*(mneq( 2) - meq( 2))
mout( 3) = s_mrt( 3)*(mneq( 3) - meq( 3))
mout( 4) = s_mrt( 4)*(mneq( 4) - meq( 4))
mout( 5) = s_mrt( 5)*(mneq( 5) - meq( 5)) 
mout( 6) = s_mrt( 6)*(mneq( 6) - meq( 6))
mout( 7) = s_mrt( 7)*(mneq( 7) - meq( 7)) 
mout( 8) = s_mrt( 8)*(mneq( 8) - meq( 8))
mout( 9) = s_mrt( 9)*(mneq( 9) - meq( 9)) 
mout(10) = s_mrt(10)*(mneq(10) - meq(10)) 
mout(11) = s_mrt(11)*(mneq(11) - meq(11)) 
mout(12) = s_mrt(12)*(mneq(12) - meq(12)) 
mout(13) = s_mrt(13)*(mneq(13) - meq(13)) 
mout(14) = s_mrt(14)*(mneq(14) - meq(14)) 
mout(15) = s_mrt(15)*(mneq(15) - meq(15)) 
mout(16) = s_mrt(16)*(mneq(16) - meq(16)) 
mout(17) = s_mrt(17)*(mneq(17)) 
mout(18) = s_mrt(18)*(mneq(18)) 
mout(19) = s_mrt(19)*(mneq(19)) 


! Rücktransformation

fOut(NDX( 1,i,j,k))= fIn(NDX( 1,i,j,k)) -(1.d0/3.d0 *mout(1) &
 -1.d0/2.d0*mout(2)+  1.d0/6.d0*mout(3))
fOut(NDX( 2,i,j,k))= fIn(NDX( 2,i,j,k)) -(1.d0/18.d0*mout(1)&
 - 1.d0/18.d0*mout(3) +1.d0/6.d0 *mout( 4) - 1.d0/6.d0 *mout( 5) &
 +1.d0/12.d0*mout(10) - 1.d0/12.d0*mout(11))
fOut(NDX( 3,i,j,k))= fIn(NDX( 3,i,j,k)) -(1.d0/18.d0*mout(1)&
 - 1.d0/18.d0*mout(3) -1.d0/6.d0 *mout( 4) + 1.d0/6.d0 *mout( 5) &
 +1.d0/12.d0*mout(10) - 1.d0/12.d0*mout(11))
fOut(NDX( 4,i,j,k))= fIn(NDX( 4,i,j,k)) -(1.d0/18.d0*mout(1)&
  -1.d0/18.d0*mout(3)+   1.d0/6.d0*mout(6) -1.d0/6.d0*mout(7)-1.d0/24.d0*mout(10)&
 + 1.d0/24.d0*mout(11)+  1.d0/8.d0*mout(12) -1.d0/8.d0*mout(13))
fOut(NDX( 5,i,j,k))= fIn(NDX( 5,i,j,k)) -(1.d0/18.d0*mout(1)&
 -1.d0/18.d0*mout(3) -1.d0/6.d0*mout(6)+  1.d0/6.d0*mout(7)&
 -1.d0/24.d0*mout(10)+ 1.d0/24.d0*mout(11)+  1.d0/8.d0*mout(12) -1.d0/8.d0*mout(13))
fOut(NDX( 6,i,j,k))= fIn(NDX( 6,i,j,k)) -(1.d0/18.d0*mout(1)&
 -1.d0/18.d0*mout(3)+        1.d0/6.d0*mout(8) -1.d0/6.d0*mout(9) &
 -1.d0/24.d0*mout(10)+ 1.d0/24.d0*mout(11) -1.d0/8.d0*mout(12)+  1.d0/8.d0*mout(13))
fOut(NDX( 7,i,j,k))= fIn(NDX( 7,i,j,k)) -(1.d0/18.d0*mout(1)&
 -1.d0/18.d0*mout(3)+       -1.d0/6.d0*mout(8)+  1.d0/6.d0*mout(9) &
 -1.d0/24.d0*mout(10)+ 1.d0/24.d0*mout(11) -1.d0/8.d0*mout(12)+  1.d0/8.d0*mout(13))
fOut(NDX( 8,i,j,k))= fIn(NDX( 8,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3)+ 1.d0/12.d0*mout(4)&
 + 1.d0/24.d0*mout(5)+ 1.d0/12.d0*mout(6)+ 1.d0/24.d0*mout(7)&
 + 1.d0/48.d0*mout(10)+ 1.d0/48.d0*mout(11)+ 1.d0/16.d0*mout(12)&
 + 1.d0/16.d0*mout(13)+  1.d0/4.d0*mout(14)+  1.d0/8.d0*mout(17) -1.d0/8.d0*mout(18))
fOut(NDX( 9,i,j,k))= fIn(NDX( 9,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3) -1.d0/12.d0*mout(4) &
 -1.d0/24.d0*mout(5) -1.d0/12.d0*mout(6) -1.d0/24.d0*mout(7)+ &
 1.d0/48.d0*mout(10)+ 1.d0/48.d0*mout(11)+ 1.d0/16.d0*mout(12)&
 + 1.d0/16.d0*mout(13)+  1.d0/4.d0*mout(14) -1.d0/8.d0*mout(17)+  1.d0/8.d0*mout(18))
fOut(NDX(10,i,j,k))= fIn(NDX(10,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3)+ 1.d0/12.d0*mout(4)&
 + 1.d0/24.d0*mout(5) -1.d0/12.d0*mout(6) -1.d0/24.d0*mout(7)&
 + 1.d0/48.d0*mout(10)+ 1.d0/48.d0*mout(11)+ 1.d0/16.d0*mout(12)&
 + 1.d0/16.d0*mout(13) -1.d0/4.d0*mout(14)+  1.d0/8.d0*mout(17)+  1.d0/8.d0*mout(18))
fOut(NDX(11,i,j,k))= fIn(NDX(11,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3) -1.d0/12.d0*mout(4) &
 -1.d0/24.d0*mout(5)+ 1.d0/12.d0*mout(6)+ 1.d0/24.d0*mout(7)&
 + 1.d0/48.d0*mout(10)+ 1.d0/48.d0*mout(11)+ 1.d0/16.d0*mout(12)+&
  1.d0/16.d0*mout(13) -1.d0/4.d0*mout(14) -1.d0/8.d0*mout(17) -1.d0/8.d0*mout(18))
fOut(NDX(12,i,j,k))= fIn(NDX(12,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3)+ 1.d0/12.d0*mout(4)&
 + 1.d0/24.d0*mout(5)+ 1.d0/12.d0*mout(8)+ 1.d0/24.d0*mout(9)&
 + 1.d0/48.d0*mout(10)+ 1.d0/48.d0*mout(11) -1.d0/16.d0*mout(12) &
 -1.d0/16.d0*mout(13)+  1.d0/4.d0*mout(16) -1.d0/8.d0*mout(17)    +  1.d0/8.d0*mout(19))
fOut(NDX(13,i,j,k))= fIn(NDX(13,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3) -1.d0/12.d0*mout(4) &
 -1.d0/24.d0*mout(5) -1.d0/12.d0*mout(8) -1.d0/24.d0*mout(9)&
 + 1.d0/48.d0*mout(10)+ 1.d0/48.d0*mout(11) -1.d0/16.d0*mout(12)&
  -1.d0/16.d0*mout(13)+  1.d0/4.d0*mout(16)+  1.d0/8.d0*mout(17)     -1.d0/8.d0*mout(19))
fOut(NDX(14,i,j,k))= fIn(NDX(14,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3)+ 1.d0/12.d0*mout(4)+ &
 1.d0/24.d0*mout(5) -1.d0/12.d0*mout(8) -1.d0/24.d0*mout(9)+ &
 1.d0/48.d0*mout(10)+ 1.d0/48.d0*mout(11) -1.d0/16.d0*mout(12) &
 -1.d0/16.d0*mout(13) -1.d0/4.d0*mout(16) -1.d0/8.d0*mout(17)     -1.d0/8.d0*mout(19))
fOut(NDX(15,i,j,k))= fIn(NDX(15,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3) -1.d0/12.d0*mout(4) &
 -1.d0/24.d0*mout(5)+ 1.d0/12.d0*mout(8)+ 1.d0/24.d0*mout(9)&
 + 1.d0/48.d0*mout(10)+ 1.d0/48.d0*mout(11) -1.d0/16.d0*mout(12) &
 -1.d0/16.d0*mout(13) -1.d0/4.d0*mout(16)+  1.d0/8.d0*mout(17)    +  1.d0/8.d0*mout(19))
fOut(NDX(16,i,j,k))= fIn(NDX(16,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3)+ 1.d0/12.d0*mout(6)&
 + 1.d0/24.d0*mout(7)+ 1.d0/12.d0*mout(8)+ 1.d0/24.d0*mout(9) &
 -1.d0/24.d0*mout(10) -1.d0/24.d0*mout(11)+      1.d0/4.d0*mout(15)+  1.d0/8.d0*mout(18) -1.d0/8.d0*mout(19))
fOut(NDX(17,i,j,k))= fIn(NDX(17,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3) -1.d0/12.d0*mout(6) &
 -1.d0/24.d0*mout(7) -1.d0/12.d0*mout(8) -1.d0/24.d0*mout(9) &
 -1.d0/24.d0*mout(10) -1.d0/24.d0*mout(11)+      1.d0/4.d0*mout(15) -1.d0/8.d0*mout(18)+  1.d0/8.d0*mout(19))
fOut(NDX(18,i,j,k))= fIn(NDX(18,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3)+ 1.d0/12.d0*mout(6)&
 + 1.d0/24.d0*mout(7) -1.d0/12.d0*mout(8) -1.d0/24.d0*mout(9) &
 -1.d0/24.d0*mout(10) -1.d0/24.d0*mout(11) -1.d0/4.d0*mout(15)+  1.d0/8.d0*mout(18)+  1.d0/8.d0*mout(19))
fOut(NDX(19,i,j,k))= fIn(NDX(19,i,j,k)) -(1.d0/36.d0*mout(1)&
 + 1.d0/24.d0*mout(2)+ 1.d0/72.d0*mout(3) -1.d0/12.d0*mout(6) &
 -1.d0/24.d0*mout(7)+ 1.d0/12.d0*mout(8)+ 1.d0/24.d0*mout(9) &
 -1.d0/24.d0*mout(10) -1.d0/24.d0*mout(11) -1.d0/4.d0*mout(15) -1.d0/8.d0*mout(18) -1.d0/8.d0*mout(19))



#ifdef MRT_MASS_MOMENTUM
call pdf_to_macro(fOut(NDX(:,i,j,k)),rhot,ut)
if(abs(check_rho) > 0.00001) then
    write(*,*) "error in rho",i,j,k," checkrho",rhot-rholoc
    stop
 endif
if(abs(check_jx ) > 0.00001) then
    write(*,*) "error in jx ",i,j,k," check_jx",ut(1)-jx(1)
    stop
 endif
if(abs(check_jy ) > 0.00001) then
    write(*,*) "error in jy ",i,j,k," checkrho",ut(2)-jx(2)
    stop
 endif
if(abs(check_jz ) > 0.00001) then
    write(*,*) "error in jz ",i,j,k," checkrho",ut(3)-jx(3)
    stop
 endif
#endif




else 
   ! simply copy the unrelaxed variables
   fOut(NDX(:,i,j,k)) = fIn(NDX(:,i,j,k))
endif
            end do
         end do
      end do
      
      call cpu_time_measure(meas%tEnd_comm)

      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration


   end subroutine collide_mrt 


#endif /* D3Q19 */



#ifdef D2Q9 
   subroutine collide_mrt(lb_dom,fOut,fIn,s_par,meas)

!------------------------------------------------------------------------
! This is the MRT collision step for the D2Q9 model
! This routine is not optimized for speed yet.
! Bulk viscosity is fixed. (relaxation factor 1.63)
! The LES SMagorinsky model is not working properly yet

      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure)              :: meas  
      real(R8B)                  :: usq,omega
      real(R8B)                  :: meq(nnod),mneq(nnod),mout(nnod),s_mrt(nnod)
      real(R8B)                  :: M_tr(nnod,nnod),M_inv(nnod,nnod),fOut_tmp(nnod)
      real(R8B)                  :: fOut_2d
      real(R8B)                  :: rhot,ut(3)
      real(R8B)                  :: weps,wepsj,wxx,rholoc_inv 
      real(R8B)                  :: ux(3),rholoc,fEq(nnod),t2cs4inv,t2cs2inv,jx(3)
      real(R8B)                  :: rho0,inv_rho0 
      real(R8B), intent(out)     :: fOut(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)   :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer                    :: i,j,k,l,ll
#ifdef MRT_MASS_MOMENTUM 
      real(R8B)  :: check_rho,check_jx,check_jy,check_jz 
#endif /* LES_SMAGORINSKY */



      call cpu_time_measure(meas%tSt_comm) 

      omega = s_par%omega

      if(s_par%initial) omega = 1.0
     s_mrt(2)      = 1.63 !Damped too much!!
!     s_mrt(2)      = omega 
!      s_mrt(2) = 2.0_R8B/(9.0_R8B * s_par%bulkViscosity + 1.0_R8B) !Bulk viscosity
      s_mrt(3)      = 1.14_R8B
      s_mrt(5)      = 1.92_R8B
      s_mrt(7)      = 1.92_R8B
      s_mrt(8)      = omega
      s_mrt(9)      = omega
      if(s_par%initial) then
         s_mrt(4)      = 1._R8B
         s_mrt(6)      = 1._R8B
      endif

M_tr(1:9,1) = (/ 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0 /)
M_tr(1:9,2) = (/-4.d0,-1.d0,-1.d0,-1.d0,-1.d0, 2.d0, 2.d0, 2.d0, 2.d0 /)
M_tr(1:9,3) = (/ 4.d0,-2.d0,-2.d0,-2.d0,-2.d0, 1.d0, 1.d0, 1.d0, 1.d0 /)
M_tr(1:9,4) = (/ 0.d0, 1.d0, 0.d0,-1.d0, 0.d0, 1.d0,-1.d0,-1.d0, 1.d0 /)
M_tr(1:9,5) = (/ 0.d0,-2.d0, 0.d0, 2.d0, 0.d0, 1.d0,-1.d0,-1.d0, 1.d0 /)
M_tr(1:9,6) = (/ 0.d0, 0.d0, 1.d0, 0.d0,-1.d0, 1.d0, 1.d0,-1.d0,-1.d0 /)
M_tr(1:9,7) = (/ 0.d0, 0.d0,-2.d0, 0.d0, 2.d0, 1.d0, 1.d0,-1.d0,-1.d0 /)
M_tr(1:9,8) = (/ 0.d0, 1.d0,-1.d0, 1.d0,-1.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_tr(1:9,9) = (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0,-1.d0, 1.d0,-1.d0 /)

M_inv(1:9,1) = (/ 1.d0/9.d0, -1.d0/9.d0,   1.d0/9.d0,  0.d0,     0.d0,    0.d0,     0.d0,     0.d0,    0.d0 /)
M_inv(1:9,2) = (/ 1.d0/9.d0, -1.d0/36.d0, -1.d0/18.d0, 1.d0/6.d0, -1.d0/6.d0,    0.d0,     0.d0,  1.d0/4.d0,    0.d0 /)
M_inv(1:9,3) = (/ 1.d0/9.d0, -1.d0/36.d0, -1.d0/18.d0, 0.d0   , 0.d0,     1.d0/6.d0, -1.d0/6.d0, -1.d0/4.d0,    0.d0 /)
M_inv(1:9,4) = (/ 1.d0/9.d0, -1.d0/36.d0, -1.d0/18.d0,-1.d0/6.d0, 1.d0/6.d0,     0.d0,     0.d0,  1.d0/4.d0,    0.d0 /)
M_inv(1:9,5) = (/ 1.d0/9.d0, -1.d0/36.d0, -1.d0/18.d0, 0.d0,    0.d0,    -1.d0/6.d0,  1.d0/6.d0, -1.d0/4.d0,    0.d0 /)
M_inv(1:9,6) = (/ 1.d0/9.d0,  1.d0/18.d0,  1.d0/36.d0, 1.d0/6.d0, 1.d0/12.d0, 1.d0/6.d0,  1.d0/12.d0,    0.d0, 1.d0/4.d0 /)
M_inv(1:9,7) = (/ 1.d0/9.d0,  1.d0/18.d0,  1.d0/36.d0,-1.d0/6.d0,-1.d0/12.d0, 1.d0/6.d0,  1.d0/12.d0,    0.d0,-1.d0/4.d0 /)
M_inv(1:9,8) = (/ 1.d0/9.d0,  1.d0/18.d0,  1.d0/36.d0,-1.d0/6.d0,-1.d0/12.d0,-1.d0/6.d0, -1.d0/12.d0,    0.d0, 1.d0/4.d0 /)
M_inv(1:9,9) = (/ 1.d0/9.d0,  1.d0/18.d0,  1.d0/36.d0, 1.d0/6.d0, 1.d0/12.d0,-1.d0/6.d0, -1.d0/12.d0,    0.d0,-1.d0/4.d0 /)
 

      ! Multiple Relaxation time Collision routine
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
          rholoc   =  fIn(NDX(I__0,i,j,k)) + fIn(NDX(I__E,i,j,k)) + &
                      fIn(NDX(I__W,i,j,k)) + fIn(NDX(I__N,i,j,k)) + &
                      fIn(NDX(I__S,i,j,k)) + fIn(NDX(I_NE,i,j,k)) + &
                      fIn(NDX(I_SW,i,j,k)) + fIn(NDX(I_SE,i,j,k)) + & 
                      fIn(NDX(I_NW,i,j,k)) 


#ifdef INCOMPRESSIBLE
rho0=s_par%rho0
#else
rho0=rholoc
#endif
! Am I currently relaxing with fixed velocity field?

              if(s_par%initial .eqv. .false.) then
! - No: Calculate macroscopic Variables (momentum and density) from pdfs
           jx(1)    = fIn(NDX(I__E,i,j,k)) - fIn(NDX(I__W,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_SE,i,j,k)) - fIn(NDX(I_NW,i,j,k))
           jx(2)    = fIn(NDX(I__N,i,j,k)) - fIn(NDX(I__S,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_NW,i,j,k)) - fIn(NDX(I_SE,i,j,k))
               else
! - Yes: Read Velocity field from Initial Field
                 jx(1) = lb_dom%u0(NDX(1,i,j,k))/rholoc
                 jx(2) = lb_dom%u0(NDX(2,i,j,k))/rholoc
               endif
! Incompressible assumption

! rho
                  meq(1) =  rholoc
! e
                  meq(2) =  -2.*rholoc + 3.*(jx(1)*jx(1) + jx(2)*jx(2))
! eps
                  meq(3) =  rholoc - 3.*(jx(1)*jx(1) + jx(2)*jx(2))
! jx
                  meq(4) =  jx(1)
! qx
                  meq(5) = -jx(1)
! jy
                  meq(6) =  jx(2)
! qy
                  meq(7) = -jx(2)
! pxx
                  meq(8) =  jx(1)*jx(1) - jx(2)*jx(2)
! pxy
                  meq(9) =  jx(1)*jx(2)

               do l=1,nnod
                  mneq(l) = 0.0d0
                  ! Berechnung der  nicht-gg-Werte im Momentenraum
                  do ll=1,nnod
                     mneq(l) = mneq(l) + fIn(NDX(ll,i,j,k))*M_tr(ll,l)
                  end do
               end do

#ifdef LES_SMAGORINSKY
! calculate additional turbulent viscosity with local non-equilibrium parts 
! acutally, neighbor information is used (after propagation step!)

drho = rholoc - rho0
!FIXME this calculation is not correct. 

s11=0.25d0*s(1)*mneq_1(i,j)+0.75d0*s(7)*mneq_7(i,j)
s12=s21=1.5d0*s(8)*mneq_8(i,j)
s22=0.25d0*s(1)*mneq_1(i,j)-0.75d0*s(7)*mneq_7(i,j) 

q_mrt(1,1) = 1.d0/3.d0*drho + jx(1)*jx(1) - 1.d0/3.d0*((mneq(2)+2.d0*drho) + mneq( 8))
q_mrt(2,2) = 1.d0/3.d0*drho + jx(2)*jx(2) - 1.d0/3.d0*((mneq(2)+2.d0*drho) - mneq( 8))

q_mrt(1,2) = mneq( 9)               
q_mrt(2,1) = mneq( 9)               
write(*,*) "smagorinsky model for 2d is not ready"
stop

q_sum = sqrt(2.d0*(q_mrt(1,1)*q_mrt(1,1) & 
            + q_mrt(2,1)*q_mrt(2,1) & 
            + q_mrt(1,2)*q_mrt(1,2) & 
            + q_mrt(2,2)*q_mrt(2,2) & 
           ) )

tau_turb = 0.5d0 *(sqrt(tau0*tau0 + 2.d0*c_smag*c_smag/rho0*cs4inv*q_sum)-tau0)
omega = 1._R8B/(tau0+tau_turb)

s_mrt( 8) = omega
s_mrt( 9) = omega
#ifdef DEBUG_LES
if(omega < omega_min) omega_min = omega
#endif /* DEBUG_LES */


#endif /* LES_SMAGORINSKY */
#ifdef SPONGE 
    s_mrt(8:9) = lb_dom%omega(i,j,k)
#endif /* SPONGE */

               ! relaxation
mout( 1) = s_mrt( 1)*(mneq( 1) - meq( 1))
mout( 2) = s_mrt( 2)*(mneq( 2) - meq( 2))
mout( 3) = s_mrt( 3)*(mneq( 3) - meq( 3))
mout( 4) = s_mrt( 4)*(mneq( 4) - meq( 4))
mout( 5) = s_mrt( 5)*(mneq( 5) - meq( 5)) 
mout( 6) = s_mrt( 6)*(mneq( 6) - meq( 6))
mout( 7) = s_mrt( 7)*(mneq( 7) - meq( 7)) 
mout( 8) = s_mrt( 8)*(mneq( 8) - meq( 8))
mout( 9) = s_mrt( 9)*(mneq( 9) - meq( 9)) 

               ! Rücktransformation
               do l=1,nnod
                  fOut_2d  = 0.0d0
                  ! Berechnung der  nicht-gg-Werte im Momentenraum
                  do ll=1,nnod
                     fOut_2d  = fOut_2d  + mout(ll)*M_inv(ll,l)
                  end do
                  fOut(NDX(l,i,j,k))= fIn(NDX(l,i,j,k)) - fOut_2d 
               end do


#ifdef MRT_MASS_MOMENTUM
call pdf_to_macro(fOut(NDX(:,i,j,k)),rhot,ut)
        if(abs(check_rho) > 0.00001) then
            write(*,*) "error in rho",i,j,k," checkrho",rhot-rholoc
            stop
         endif
        if(abs(check_jx ) > 0.00001) then
            write(*,*) "error in jx ",i,j,k," check_jx",ut(1)-jx(1)
            stop
         endif
        if(abs(check_jy ) > 0.00001) then
            write(*,*) "error in jy ",i,j,k," checkrho",ut(2)-jx(2)
            stop
         endif
        if(abs(check_jz ) > 0.00001) then
            write(*,*) "error in jz ",i,j,k," checkrho",ut(3)-jx(3)
            stop
         endif

#endif
         else
             ! copy the boundary nodes!!
             fOut(NDX(:,i,j,k)) = fIn(NDX(:,i,j,k))
         endif
            end do
         end do
      end do
      
      call cpu_time_measure(meas%tEnd_comm)
      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration


   end subroutine collide_mrt
  !------------------------------------------------------------------------
#endif /* D2Q9  */

#endif /* MRT */
















#ifndef COMBINE_STREAM_COLLIDE
   subroutine collide_bgk(lb_dom,fOut,fIn,s_par,meas)

!------------------------------------------------------------------------
! This is the BGK collision step
! used by default



      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure)              :: meas  
      real(R8B)                  :: usq,omega
      real(R8B)                  :: u(3),rholoc,fEq(nnod),t2cs4inv,t2cs2inv
      real(R8B)                  :: rho0 
      real(R8B), intent(out)     :: fOut(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)   :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer                    :: i,j,k,l
#ifdef LES_SMAGORINSKY
      integer                    :: is,js,ks,ls
      real(R8B)  ::  tau0,q_bgk(NDIM,NDIM),q_sum,tau_turb,cs4inv
      real(R8B),parameter :: sqrt2  = 1.4142135623731d0 
      real(R8B),parameter :: c_smag = 0.16d0 
#ifdef DEBUG_LES
      real(R8B)  ::  omega_min 
      omega_min = 3.0 
#endif /* LES_SMAGORINSKY */

#endif /* LES_SMAGORINSKY */


      call cpu_time_measure(meas%tSt_comm) 
      omega = s_par%omega
      rho0=s_par%rho0

      t2cs4inv = 1._R8B/(2._R8B*cs*cs*cs*cs)
      t2cs2inv = 1._R8B/(2._R8B*cs*cs)
#ifdef LES_SMAGORINSKY
      tau0 = 1.d0/omega
      cs4inv = t2cs4inv*2._R8B 
#endif /* LES_SMAGORINSKY */



      ! Single Relaxation time Collision routine
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
        !       if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
#ifdef SPONGE
! In the sponge layer, the viscosity is increased up to the border to reduce reflections
! Not working very well so far....
            omega = lb_dom%omega(i,j,k)
#endif



! Calculate current Density 
#ifdef D3Q19
       rholoc = fIn(NDX(I__0,i,j,k)) + &
                fIn(NDX(I__E,i,j,k)) + fIn(NDX(I__W,i,j,k)) + &
                fIn(NDX(I__N,i,j,k)) + fIn(NDX(I__S,i,j,k)) + &
                fIn(NDX(I__T,i,j,k)) + fIn(NDX(I__B,i,j,k)) + &
                fIn(NDX(I_NE,i,j,k)) + fIn(NDX(I_SW,i,j,k)) + &
                fIn(NDX(I_SE,i,j,k)) + fIn(NDX(I_NW,i,j,k)) + &
                fIn(NDX(I_TE,i,j,k)) + fIn(NDX(I_BW,i,j,k)) + &
                fIn(NDX(I_BE,i,j,k)) + fIn(NDX(I_TW,i,j,k)) + &
                fIn(NDX(I_TN,i,j,k)) + fIn(NDX(I_BS,i,j,k)) + &
                fIn(NDX(I_BN,i,j,k)) + fIn(NDX(I_TS,i,j,k))
#else /* D2Q9 */
       rholoc = fIn(NDX(I__0,i,j,k)) + fIn(NDX(I__E,i,j,k))  + &
                fIn(NDX(I__W,i,j,k)) + fIn(NDX(I__N,i,j,k))  + &
                fIn(NDX(I__S,i,j,k)) + fIn(NDX(I_NE,i,j,k))  + &
                fIn(NDX(I_SW,i,j,k)) + fIn(NDX(I_SE,i,j,k))  + &
                fIn(NDX(I_NW,i,j,k))
#endif /* D3Q9 */


! Choose which reference density is used. 
! Constant density rho0 for incompressible model

#ifdef INCOMPRESSIBLE
   rho0 = s_par%rho0 
#else
   rho0 = rholoc
#endif
  
! Calculate Momentum

#ifdef D3Q19
              if(s_par%initial .eqv. .false.) then
    u(1) = fIn(NDX(I__E,i,j,k)) - fIn(NDX(I__W,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_SE,i,j,k)) - fIn(NDX(I_NW,i,j,k)) + &
                      fIn(NDX(I_TE,i,j,k)) - fIn(NDX(I_BW,i,j,k)) + &
                      fIn(NDX(I_BE,i,j,k)) - fIn(NDX(I_TW,i,j,k))  
    u(2) = fIn(NDX(I__N,i,j,k)) - fIn(NDX(I__S,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_NW,i,j,k)) - fIn(NDX(I_SE,i,j,k)) + & 
                      fIn(NDX(I_TN,i,j,k)) - fIn(NDX(I_BS,i,j,k)) + &
                      fIn(NDX(I_BN,i,j,k)) - fIn(NDX(I_TS,i,j,k))
    u(3) = fIn(NDX(I__T,i,j,k)) - fIn(NDX(I__B,i,j,k))  + &
                      fIn(NDX(I_TE,i,j,k)) - fIn(NDX(I_BW,i,j,k))  + &
                      fIn(NDX(I_TW,i,j,k)) - fIn(NDX(I_BE,i,j,k))  + &
                      fIn(NDX(I_TN,i,j,k)) - fIn(NDX(I_BS,i,j,k))  + &
                      fIn(NDX(I_TS,i,j,k)) - fIn(NDX(I_BN,i,j,k))
   u=u/rho0
               else
         ! Decide, if a constant velocity field is used for initial relaxation process
                  u(1) = lb_dom%u0(NDX(1,i,j,k))
                  u(2) = lb_dom%u0(NDX(2,i,j,k))
                  u(3) = lb_dom%u0(NDX(3,i,j,k))
               endif
 
               usq =  (u(1)**2 + u(2)**2 + u(3)**2)*t2cs2inv

#else /* D3Q19 -> D2Q9*/ 

              if(s_par%initial .eqv. .false.) then
         u(1) =  fIn(NDX(I__E,i,j,k)) - fIn(NDX(I__W,i,j,k))  + &
                 fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k))  + &
                 fIn(NDX(I_SE,i,j,k)) - fIn(NDX(I_NW,i,j,k))    
         u(2) =  fIn(NDX(I__N,i,j,k)) - fIn(NDX(I__S,i,j,k))  + &
                 fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k))  + &  
                 fIn(NDX(I_NW,i,j,k)) - fIn(NDX(I_SE,i,j,k) )    
         u = u/rho0
               else
         ! Decide, if a constant velocity field is used for initial relaxation process
                  u(1) = lb_dom%u0(NDX(1,i,j,k))
                  u(2) = lb_dom%u0(NDX(2,i,j,k))
               endif

                usq =  (u(1)*u(1) + u(2)*u(2))*t2cs2inv

#endif /* D3Q19 */

! Calculate equilibrium distributions

#ifdef D3Q19
               fEq(1) = t(1)*(rholoc -rho0*usq)
 
               fEq(2) = t(2)*(rholoc + rho0*((u(1))*cs2inv &
                      + (u(1))**2*t2cs4inv  &
                      - usq))
               fEq(3) = t(3)*(rholoc + rho0*((-u(1))*cs2inv &
                      + (-u(1))**2*t2cs4inv  &
                      - usq) )
               fEq(4) = t(4)*(rholoc + rho0*((u(2))*cs2inv &
                      + (u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(5) = t(5)*(rholoc + rho0*((-u(2))*cs2inv &
                      + (-u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(6) = t(6)*(rholoc + rho0*((u(3))*cs2inv &
                      + (u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(7) = t(7)*(rholoc + rho0*((-u(3))*cs2inv &
                      + (-u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(8) = t(8)*(rholoc + rho0*((u(1) + u(2))*cs2inv &
                      + (u(1) + u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(9) = t(9)*(rholoc + rho0*((-u(1) -u(2))*cs2inv &
                      + (-u(1) -u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(10) = t(10)*(rholoc + rho0*((u(1) -u(2))*cs2inv &
                      + (u(1) -u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(11) = t(11)*(rholoc + rho0*((-u(1) + u(2))*cs2inv &
                      + (-u(1) + u(2))**2*t2cs4inv  &
                      - usq)) 
               fEq(12) = t(12)*(rholoc + rho0*((u(1) + u(3))*cs2inv &
                      + (u(1) + u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(13) = t(13)*(rholoc + rho0*((-u(1) -u(3))*cs2inv &
                      + (-u(1) -u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(14) = t(14)*(rholoc + rho0*((u(1) -u(3))*cs2inv &
                      + (u(1) -u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(15) = t(15)*(rholoc + rho0*((-u(1) + u(3))*cs2inv &
                      + (-u(1) + u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(16) = t(16)*(rholoc + rho0*((u(2) + u(3))*cs2inv &
                      + (u(2) + u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(17) = t(17)*(rholoc + rho0*((-u(2) -u(3))*cs2inv &
                      + (-u(2) -u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(18) = t(18)*(rholoc + rho0*((u(2) -u(3))*cs2inv &
                      + (u(2)  -u(3))**2*t2cs4inv  &
                      - usq)) 
               fEq(19) = t(19)*(rholoc + rho0*((-u(2) + u(3))*cs2inv &
                      + (-u(2) + u(3))**2*t2cs4inv  &
                      - usq)) 
#else /* D2Q9 */
                fEq(1) = t(1)*(rholoc  - rho0*usq)
                fEq(2) = t(2)*(rholoc  + rho0*((u(1))*cs2inv &
                      + (u(1))**2*t2cs4inv & 
                      - usq))
                fEq(3) = t(3)*(rholoc + rho0*((u(2))*cs2inv &
                      + (u(2))**2*t2cs4inv  &
                      - usq))
                fEq(4) = t(4)*(rholoc + rho0*((-u(1) )*cs2inv &
                      + (-u(1))**2*t2cs4inv  &
                      - usq))
                fEq(5) = t(5)*(rholoc + rho0*((-u(2))*cs2inv &
                      + (-u(2))**2*t2cs4inv  &
                      - usq))
                fEq(6) = t(6)*(rholoc + rho0*((u(1) + u(2))*cs2inv &
                      + (u(1)+u(2))**2*t2cs4inv  &
                      - usq))
                fEq(7) = t(7)*(rholoc + rho0*((-u(1) + u(2))*cs2inv &
                      + (-u(1)+u(2))**2*t2cs4inv  &
                      - usq))
                fEq(8) = t(8)*(rholoc + rho0*((-u(1) -u(2))*cs2inv &
                      + (-u(1)-u(2))**2*t2cs4inv  &
                      - usq))
                fEq(9) = t(9)*(rholoc + rho0*((u(1) -u(2))*cs2inv &
                      + (u(1)-u(2))**2*t2cs4inv  &
                      - usq))

#endif /* D3Q19 */

#ifdef LES_SMAGORINSKY
! calculate additional turbulent viscosity with local non-equilibrium parts 
! acutally, neighbor information is used (after propagation step!)
q_bgk = 0.d0
do js = 1,NDIM
do is = 1,NDIM
do ls = 1,nnod
q_bgk(is,js)=q_bgk(is,js) + cx(ls,is)*cx(ls,js)*(fIn(NDX(ls,i,j,k))-fEq(ls)) 
enddo
enddo
enddo

#ifdef UNROLL
q_bgk(1,1) = cx( 1,1)*cx( 1,1)*(fIn(NDX( 1,i,j,k))-fEq( 1)) &
           + cx( 2,1)*cx( 2,1)*(fIn(NDX( 2,i,j,k))-fEq( 2)) &
           + cx( 3,1)*cx( 3,1)*(fIn(NDX( 3,i,j,k))-fEq( 3)) &
           + cx( 4,1)*cx( 4,1)*(fIn(NDX( 4,i,j,k))-fEq( 4)) &
           + cx( 5,1)*cx( 5,1)*(fIn(NDX( 5,i,j,k))-fEq( 5)) &
           + cx( 6,1)*cx( 6,1)*(fIn(NDX( 6,i,j,k))-fEq( 6)) &
           + cx( 7,1)*cx( 7,1)*(fIn(NDX( 7,i,j,k))-fEq( 7)) &
           + cx( 8,1)*cx( 8,1)*(fIn(NDX( 8,i,j,k))-fEq( 8)) &
           + cx( 9,1)*cx( 9,1)*(fIn(NDX( 9,i,j,k))-fEq( 9)) &
#ifdef D3Q19
           + cx(10,1)*cx(10,1)*(fIn(NDX(10,i,j,k))-fEq(10)) &
           + cx(11,1)*cx(11,1)*(fIn(NDX(11,i,j,k))-fEq(11)) &
           + cx(12,1)*cx(12,1)*(fIn(NDX(12,i,j,k))-fEq(12)) &
           + cx(13,1)*cx(13,1)*(fIn(NDX(13,i,j,k))-fEq(13)) &
           + cx(14,1)*cx(14,1)*(fIn(NDX(14,i,j,k))-fEq(14)) &
           + cx(15,1)*cx(15,1)*(fIn(NDX(15,i,j,k))-fEq(15)) &
           + cx(16,1)*cx(16,1)*(fIn(NDX(16,i,j,k))-fEq(16)) &
           + cx(17,1)*cx(17,1)*(fIn(NDX(17,i,j,k))-fEq(17)) &
           + cx(18,1)*cx(18,1)*(fIn(NDX(18,i,j,k))-fEq(18)) &
           + cx(19,1)*cx(19,1)*(fIn(NDX(19,i,j,k))-fEq(19)) 
#else 
           + 0.d0
#endif
#endif /* UNROLL */
!q_sum = 0.d0
!do js = 1,NDIM
!do is = 1,NDIM
!q_sum=q_sum + 2.d0*q_bgk(is,js)*q_bgk(is,js)
!enddo
!enddo
q_sum = sqrt(2.d0*(q_bgk(1,1)*q_bgk(1,1) & 
            + q_bgk(2,1)*q_bgk(2,1) & 
            + q_bgk(1,2)*q_bgk(1,2) & 
            + q_bgk(2,2)*q_bgk(2,2) & 
#ifdef D3Q19
            + q_bgk(3,1)*q_bgk(3,1) & 
            + q_bgk(3,2)*q_bgk(3,2) & 
            + q_bgk(1,3)*q_bgk(1,3) & 
            + q_bgk(2,3)*q_bgk(2,3) & 
            + q_bgk(3,3)*q_bgk(3,3) ) )
#else /* D2Q9 */
           ) )
#endif /* D3Q19 */
!q_sum = sqrt(q_sum)

tau_turb = 0.5d0 *(sqrt(tau0*tau0 + 2.d0*c_smag*c_smag/rho0*cs4inv*q_sum)-tau0)
omega = 1._R8B/(tau0+tau_turb)

#ifdef DEBUG_LES
if(omega < omega_min) omega_min = omega
#endif /* DEBUG_LES */

#endif /* LES_SMAGORINSKY */

! Finally, relax to local equilibrium

               fOut(NDX(1,i,j,k)) = fIn(NDX(1,i,j,k)) - omega*(fIn(NDX(1,i,j,k)) - fEq(1))
               fOut(NDX(2,i,j,k)) = fIn(NDX(2,i,j,k)) - omega*(fIn(NDX(2,i,j,k)) - fEq(2))
               fOut(NDX(3,i,j,k)) = fIn(NDX(3,i,j,k)) - omega*(fIn(NDX(3,i,j,k)) - fEq(3))
               fOut(NDX(4,i,j,k)) = fIn(NDX(4,i,j,k)) - omega*(fIn(NDX(4,i,j,k)) - fEq(4))
               fOut(NDX(5,i,j,k)) = fIn(NDX(5,i,j,k)) - omega*(fIn(NDX(5,i,j,k)) - fEq(5))
               fOut(NDX(6,i,j,k)) = fIn(NDX(6,i,j,k)) - omega*(fIn(NDX(6,i,j,k)) - fEq(6))
               fOut(NDX(7,i,j,k)) = fIn(NDX(7,i,j,k)) - omega*(fIn(NDX(7,i,j,k)) - fEq(7))
               fOut(NDX(8,i,j,k)) = fIn(NDX(8,i,j,k)) - omega*(fIn(NDX(8,i,j,k)) - fEq(8))
               fOut(NDX(9,i,j,k)) = fIn(NDX(9,i,j,k)) - omega*(fIn(NDX(9,i,j,k)) - fEq(9))
#ifdef D3Q19
               fOut(NDX(10,i,j,k)) = fIn(NDX(10,i,j,k)) - omega*(fIn(NDX(10,i,j,k)) - fEq(10))
               fOut(NDX(11,i,j,k)) = fIn(NDX(11,i,j,k)) - omega*(fIn(NDX(11,i,j,k)) - fEq(11))
               fOut(NDX(12,i,j,k)) = fIn(NDX(12,i,j,k)) - omega*(fIn(NDX(12,i,j,k)) - fEq(12))
               fOut(NDX(13,i,j,k)) = fIn(NDX(13,i,j,k)) - omega*(fIn(NDX(13,i,j,k)) - fEq(13))
               fOut(NDX(14,i,j,k)) = fIn(NDX(14,i,j,k)) - omega*(fIn(NDX(14,i,j,k)) - fEq(14))
               fOut(NDX(15,i,j,k)) = fIn(NDX(15,i,j,k)) - omega*(fIn(NDX(15,i,j,k)) - fEq(15))
               fOut(NDX(16,i,j,k)) = fIn(NDX(16,i,j,k)) - omega*(fIn(NDX(16,i,j,k)) - fEq(16))
               fOut(NDX(17,i,j,k)) = fIn(NDX(17,i,j,k)) - omega*(fIn(NDX(17,i,j,k)) - fEq(17))
               fOut(NDX(18,i,j,k)) = fIn(NDX(18,i,j,k)) - omega*(fIn(NDX(18,i,j,k)) - fEq(18))
               fOut(NDX(19,i,j,k)) = fIn(NDX(19,i,j,k)) - omega*(fIn(NDX(19,i,j,k)) - fEq(19))
#endif /* D3Q19 */

#ifdef DEBUG
if(i==5 .and. j==8 .and. k==5) then
write(*,*) "u",u
write(*,*) "rholoc",rholoc
write(*,*) "fEq"
do l=1,nnod
write(*,*) fEq(l) 
enddo
write(*,*) "fIn"
do l=1,nnod
   write(*,*) fIn(NDX(l,5,8,5)) 
enddo
write(*,*) "fOut"
do l=1,nnod
   write(*,*) fOut(NDX(l,5,8,5)) 
enddo
endif
#endif /* DEBUG */

!else 
   ! simply copy the unrelaxed variables
!   fOut(NDX(:,i,j,k)) = fIn(NDX(:,i,j,k))
!endif
            end do
         end do
      end do
      
      call cpu_time_measure(meas%tEnd_comm)

      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration
#ifdef DEBUG_LES
      write(66,*) gtstep_cur,omega_min
#endif /* LES_SMAGORINSKY */


   end subroutine collide_bgk 
  !------------------------------------------------------------------------
#endif /* COMBINE_STREAM_COLLIDE */


#ifndef TRT
   subroutine stream_collide_bgk(lb_dom,fOut,fIn,s_par,meas)

!------------------------------------------------------------------------
! This is the BGK collision step
! with combined stream collide 


      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure)              :: meas  
      real (R8B)                 :: t2cs4inv,t2cs2inv
      real(R8B), intent(out)     :: fOut(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)   ::  fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer                    :: i,j,k,l
#ifdef SPLITLOOPS
#define INDEXVAR i
#define SPLITDIMENSION lb_dom%lx(1)
#else
#define INDEXVAR 1
#define SPLITDIMENSION 1 
#endif
      real(R8B),dimension(SPLITDIMENSION) :: rholoc,u1,u2,u3,omega,usq
      real(R8B),dimension(SPLITDIMENSION) :: ftmp1,ftmp2,ftmp3,ftmp4,ftmp5,ftmp6,ftmp7,ftmp8,ftmp9
      real(R8B)                           :: fEq1,fEq2,fEq3,fEq4,fEq5,fEq6,fEq7,fEq8,fEq9
#ifdef D3Q19
      real(R8B),dimension(SPLITDIMENSION) :: ftmp10,ftmp11,ftmp12,ftmp13,ftmp14,ftmp15,ftmp16,ftmp17,ftmp18,ftmp19
      real(R8B)                           :: fEq10,fEq11,fEq12,fEq13,fEq14,fEq15,fEq16,fEq17,fEq18,fEq19
#endif
      call cpu_time_measure(meas%tSt_comm) 

#ifndef SPONGE
      omega = s_par%omega
#endif
      t2cs4inv = 1._R8B/(2._R8B*cs*cs*cs*cs)
      t2cs2inv = 1._R8B/(2._R8B*cs*cs)

      ! Single Relaxation time Collision routine
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)

            !--------------------------------------------
            ! Do Streaming process ( Pull values) 

                  ftmp1(INDEXVAR)=fIn(NDX(1,i-cx(1,1),j-cx(1,2),k-cx(1,3)))
                  ftmp2(INDEXVAR)=fIn(NDX(2,i-cx(2,1),j-cx(2,2),k-cx(2,3)))
                  ftmp3(INDEXVAR)=fIn(NDX(3,i-cx(3,1),j-cx(3,2),k-cx(3,3)))
                  ftmp4(INDEXVAR)=fIn(NDX(4,i-cx(4,1),j-cx(4,2),k-cx(4,3)))
                  ftmp5(INDEXVAR)=fIn(NDX(5,i-cx(5,1),j-cx(5,2),k-cx(5,3)))
                  ftmp6(INDEXVAR)=fIn(NDX(6,i-cx(6,1),j-cx(6,2),k-cx(6,3)))
                  ftmp7(INDEXVAR)=fIn(NDX(7,i-cx(7,1),j-cx(7,2),k-cx(7,3)))
                  ftmp8(INDEXVAR)=fIn(NDX(8,i-cx(8,1),j-cx(8,2),k-cx(8,3)))
                  ftmp9(INDEXVAR)=fIn(NDX(9,i-cx(9,1),j-cx(9,2),k-cx(9,3)))
#ifdef D3Q19
                  ftmp10(INDEXVAR)=fIn(NDX(10,i-cx(10,1),j-cx(10,2),k-cx(10,3)))
                  ftmp11(INDEXVAR)=fIn(NDX(11,i-cx(11,1),j-cx(11,2),k-cx(11,3)))
                  ftmp12(INDEXVAR)=fIn(NDX(12,i-cx(12,1),j-cx(12,2),k-cx(12,3)))
                  ftmp13(INDEXVAR)=fIn(NDX(13,i-cx(13,1),j-cx(13,2),k-cx(13,3)))
                  ftmp14(INDEXVAR)=fIn(NDX(14,i-cx(14,1),j-cx(14,2),k-cx(14,3)))
                  ftmp15(INDEXVAR)=fIn(NDX(15,i-cx(15,1),j-cx(15,2),k-cx(15,3)))
                  ftmp16(INDEXVAR)=fIn(NDX(16,i-cx(16,1),j-cx(16,2),k-cx(16,3)))
                  ftmp17(INDEXVAR)=fIn(NDX(17,i-cx(17,1),j-cx(17,2),k-cx(17,3)))
                  ftmp18(INDEXVAR)=fIn(NDX(18,i-cx(18,1),j-cx(18,2),k-cx(18,3)))
                  ftmp19(INDEXVAR)=fIn(NDX(19,i-cx(19,1),j-cx(19,2),k-cx(19,3)))
#endif

            ! Split loops. (Thomas: 3 und 5 splitloops)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
#ifdef SPONGE
               omega(INDEXVAR) = lb_dom%omega(i,j,k)
#endif


         !--------------------------------------------
         ! Calculate macroscopic variables 

#ifdef D2Q9
              rholoc(INDEXVAR) =   ftmp1(INDEXVAR) + ftmp2(INDEXVAR) + & 
                         ftmp3(INDEXVAR) + ftmp4(INDEXVAR) + &
                         ftmp5(INDEXVAR) + ftmp6(INDEXVAR) + & 
                         ftmp7(INDEXVAR) + ftmp8(INDEXVAR) + ftmp9(INDEXVAR)
               u1(INDEXVAR) = 1/rholoc(INDEXVAR) * (ftmp2(INDEXVAR)  - ftmp4(INDEXVAR) + &
                         ftmp6(INDEXVAR) - ftmp7(INDEXVAR) - & 
                         ftmp8(INDEXVAR) + ftmp9(INDEXVAR))
               u2(INDEXVAR) = 1/rholoc(INDEXVAR) * (ftmp3(INDEXVAR)  - ftmp5(INDEXVAR) + &
                        ftmp6(INDEXVAR) + ftmp7(INDEXVAR) - & 
                        ftmp8(INDEXVAR) - ftmp9(INDEXVAR))
               u3(INDEXVAR) = 0.d0
#endif
#ifdef D3Q19
               rholoc(INDEXVAR)  = ftmp1(INDEXVAR)  + ftmp2(INDEXVAR)  + &
                         ftmp3(INDEXVAR)  + ftmp4(INDEXVAR)  + ftmp5(INDEXVAR) + &
                         ftmp6(INDEXVAR)  + ftmp7(INDEXVAR)  + &
                         ftmp8(INDEXVAR)  + ftmp9(INDEXVAR)  + ftmp10(INDEXVAR) + & 
                         ftmp11(INDEXVAR) + ftmp12(INDEXVAR) + &
                         ftmp13(INDEXVAR) + ftmp14(INDEXVAR) + ftmp15(INDEXVAR) + & 
                         ftmp16(INDEXVAR) + ftmp17(INDEXVAR) + &
                         ftmp18(INDEXVAR) + ftmp19(INDEXVAR)

               u1(INDEXVAR) = (  ftmp2(INDEXVAR)  - ftmp3(INDEXVAR)  + &
                         ftmp8(INDEXVAR)  - ftmp9(INDEXVAR)  + &
                         ftmp10(INDEXVAR) - ftmp11(INDEXVAR) + & 
                         ftmp12(INDEXVAR) - ftmp13(INDEXVAR) + &
                         ftmp14(INDEXVAR) - ftmp15(INDEXVAR)  & 
                             )/rholoc(INDEXVAR)    

               u2(INDEXVAR) = (  ftmp4(INDEXVAR)  - ftmp5(INDEXVAR) + &
                         ftmp8(INDEXVAR)  - ftmp9(INDEXVAR)  - ftmp10(INDEXVAR) + & 
                         ftmp11(INDEXVAR) + & 
                         ftmp16(INDEXVAR) - ftmp17(INDEXVAR) + &
                         ftmp18(INDEXVAR) - ftmp19(INDEXVAR))/rholoc(INDEXVAR)    

               u3(INDEXVAR) = (  ftmp6(INDEXVAR)  - ftmp7(INDEXVAR)  + &
                         ftmp12(INDEXVAR) - ftmp13(INDEXVAR)- ftmp14(INDEXVAR)  + & 
                         ftmp15(INDEXVAR) + ftmp16(INDEXVAR) - &
                         ftmp17(INDEXVAR) - ftmp18(INDEXVAR) + &
                         ftmp19(INDEXVAR))/rholoc(INDEXVAR)    
#endif
 
               usq(INDEXVAR)=  (u1(INDEXVAR)**2 + u2(INDEXVAR)**2 + u3(INDEXVAR)**2)*t2cs2inv


         !--------------------------------------------
         ! Calculate equilibrium distribution


#ifdef D2Q9
                fEq1 = t(1)*rholoc(INDEXVAR) *(1._R8B - usq(INDEXVAR))
                fEq2 = t(2)*rholoc(INDEXVAR) *(1._R8B + (u1(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR))**2*t2cs4inv & 
                      - usq(INDEXVAR))
                fEq3 = t(3)*rholoc(INDEXVAR)*(1._R8B + (u2(INDEXVAR))*cs2inv &
                      + (u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR))
                fEq4 = t(4)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) )*cs2inv &
                      + (-u1(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR))
                fEq5 = t(5)*rholoc(INDEXVAR)*(1._R8B + (-u2(INDEXVAR))*cs2inv &
                      + (-u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR))


               fOut(NDX(1,i,j,k)) = ftmp1(INDEXVAR) - omega(INDEXVAR)*(ftmp1(INDEXVAR) - fEq1)
               fOut(NDX(2,i,j,k)) = ftmp2(INDEXVAR) - omega(INDEXVAR)*(ftmp2(INDEXVAR) - fEq2)
               fOut(NDX(3,i,j,k)) = ftmp3(INDEXVAR) - omega(INDEXVAR)*(ftmp3(INDEXVAR) - fEq3)
               fOut(NDX(4,i,j,k)) = ftmp4(INDEXVAR) - omega(INDEXVAR)*(ftmp4(INDEXVAR) - fEq4)
               fOut(NDX(5,i,j,k)) = ftmp5(INDEXVAR) - omega(INDEXVAR)*(ftmp5(INDEXVAR) - fEq5)

#ifdef SPLITLOOPS
            endif ! not wall 
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
#endif
                fEq6 = t(6)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR) + u2(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR)+u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR))
                fEq7 = t(7)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) + u2(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR)+u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR))
                fEq8 = t(8)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) -u2(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR)-u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR))
                fEq9 = t(9)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR) -u2(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR)-u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR))
               fOut(NDX(6,i,j,k)) = ftmp6(INDEXVAR) - omega(INDEXVAR)*(ftmp6(INDEXVAR) - fEq6)
               fOut(NDX(7,i,j,k)) = ftmp7(INDEXVAR) - omega(INDEXVAR)*(ftmp7(INDEXVAR) - fEq7)
               fOut(NDX(8,i,j,k)) = ftmp8(INDEXVAR) - omega(INDEXVAR)*(ftmp8(INDEXVAR) - fEq8)
               fOut(NDX(9,i,j,k)) = ftmp9(INDEXVAR) - omega(INDEXVAR)*(ftmp9(INDEXVAR) - fEq9)
#endif
#ifdef D3Q19
               fEq1 = t(1)*rholoc(INDEXVAR)*(1._R8B  &
                      - usq(INDEXVAR)) 
               fEq2 = t(2)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq3 = t(3)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq4 = t(4)*rholoc(INDEXVAR)*(1._R8B + (u2(INDEXVAR))*cs2inv &
                      + (u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq5 = t(5)*rholoc(INDEXVAR)*(1._R8B + (-u2(INDEXVAR))*cs2inv &
                      + (-u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq6 = t(6)*rholoc(INDEXVAR)*(1._R8B + (u3(INDEXVAR))*cs2inv &
                      + (u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq7 = t(7)*rholoc(INDEXVAR)*(1._R8B + (-u3(INDEXVAR))*cs2inv &
                      + (-u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fOut(NDX(1,i,j,k)) = ftmp1(INDEXVAR) - omega(INDEXVAR)*(ftmp1(INDEXVAR) - fEq1)
               fOut(NDX(2,i,j,k)) = ftmp2(INDEXVAR) - omega(INDEXVAR)*(ftmp2(INDEXVAR) - fEq2)
               fOut(NDX(3,i,j,k)) = ftmp3(INDEXVAR) - omega(INDEXVAR)*(ftmp3(INDEXVAR) - fEq3)
               fOut(NDX(4,i,j,k)) = ftmp4(INDEXVAR) - omega(INDEXVAR)*(ftmp4(INDEXVAR) - fEq4)
               fOut(NDX(5,i,j,k)) = ftmp5(INDEXVAR) - omega(INDEXVAR)*(ftmp5(INDEXVAR) - fEq5)
               fOut(NDX(6,i,j,k)) = ftmp6(INDEXVAR) - omega(INDEXVAR)*(ftmp6(INDEXVAR) - fEq6)
               fOut(NDX(7,i,j,k)) = ftmp7(INDEXVAR) - omega(INDEXVAR)*(ftmp7(INDEXVAR) - fEq7)

#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
#endif



               fEq8 = t(8)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR) + u2(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR) + u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq9 = t(9)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) -u2(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR) -u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq10 = t(10)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR) -u2(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR) -u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq11 = t(11)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) + u2(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR) + u2(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq12 = t(12)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR) + u3(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR) + u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq13 = t(13)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) -u3(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR) -u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fOut(NDX(8,i,j,k))  = ftmp8(INDEXVAR)  - omega(INDEXVAR)*(ftmp8(INDEXVAR)  -  fEq8)
               fOut(NDX(9,i,j,k))  = ftmp9(INDEXVAR)  - omega(INDEXVAR)*(ftmp9(INDEXVAR)  -  fEq9)
               fOut(NDX(10,i,j,k)) = ftmp10(INDEXVAR) - omega(INDEXVAR)*(ftmp10(INDEXVAR) - fEq10)
               fOut(NDX(11,i,j,k)) = ftmp11(INDEXVAR) - omega(INDEXVAR)*(ftmp11(INDEXVAR) - fEq11)
               fOut(NDX(12,i,j,k)) = ftmp12(INDEXVAR) - omega(INDEXVAR)*(ftmp12(INDEXVAR) - fEq12)
               fOut(NDX(13,i,j,k)) = ftmp13(INDEXVAR) - omega(INDEXVAR)*(ftmp13(INDEXVAR) - fEq13)

#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
#endif


               fEq14 = t(14)*rholoc(INDEXVAR)*(1._R8B + (u1(INDEXVAR) -u3(INDEXVAR))*cs2inv &
                      + (u1(INDEXVAR) -u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq15 = t(15)*rholoc(INDEXVAR)*(1._R8B + (-u1(INDEXVAR) + u3(INDEXVAR))*cs2inv &
                      + (-u1(INDEXVAR) + u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq16 = t(16)*rholoc(INDEXVAR)*(1._R8B + (u2(INDEXVAR) + u3(INDEXVAR))*cs2inv &
                      + (u2(INDEXVAR) + u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq17 = t(17)*rholoc(INDEXVAR)*(1._R8B + (-u2(INDEXVAR) -u3(INDEXVAR))*cs2inv &
                      + (-u2(INDEXVAR) -u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq18 = t(18)*rholoc(INDEXVAR)*(1._R8B + (u2(INDEXVAR) -u3(INDEXVAR))*cs2inv &
                      + (u2(INDEXVAR)  -u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fEq19 = t(19)*rholoc(INDEXVAR)*(1._R8B + (-u2(INDEXVAR) + u3(INDEXVAR))*cs2inv &
                      + (-u2(INDEXVAR) + u3(INDEXVAR))**2*t2cs4inv  &
                      - usq(INDEXVAR)) 
               fOut(NDX(14,i,j,k)) = ftmp14(INDEXVAR) - omega(INDEXVAR)*(ftmp14(INDEXVAR) - fEq14)
               fOut(NDX(15,i,j,k)) = ftmp15(INDEXVAR) - omega(INDEXVAR)*(ftmp15(INDEXVAR) - fEq15)
               fOut(NDX(16,i,j,k)) = ftmp16(INDEXVAR) - omega(INDEXVAR)*(ftmp16(INDEXVAR) - fEq16)
               fOut(NDX(17,i,j,k)) = ftmp17(INDEXVAR) - omega(INDEXVAR)*(ftmp17(INDEXVAR) - fEq17)
               fOut(NDX(18,i,j,k)) = ftmp18(INDEXVAR) - omega(INDEXVAR)*(ftmp18(INDEXVAR) - fEq18)
               fOut(NDX(19,i,j,k)) = ftmp19(INDEXVAR) - omega(INDEXVAR)*(ftmp19(INDEXVAR) - fEq19)



#endif

else

         !--------------------------------------------
         ! Do Bounce Back 
#ifdef D2Q9
         fOut(NDX(opp(1),i,j,k)) = ftmp1(INDEXVAR)
         fOut(NDX(opp(2),i,j,k)) = ftmp2(INDEXVAR)
         fOut(NDX(opp(3),i,j,k)) = ftmp3(INDEXVAR)
         fOut(NDX(opp(4),i,j,k)) = ftmp4(INDEXVAR)
         fOut(NDX(opp(5),i,j,k)) = ftmp5(INDEXVAR)

#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
#endif


         fOut(NDX(opp(6),i,j,k)) = ftmp6(INDEXVAR)
         fOut(NDX(opp(7),i,j,k)) = ftmp7(INDEXVAR)
         fOut(NDX(opp(8),i,j,k)) = ftmp8(INDEXVAR)
         fOut(NDX(opp(9),i,j,k)) = ftmp9(INDEXVAR)
#endif

#ifdef D3Q19
         fOut(NDX(opp(1),i,j,k)) = ftmp1(INDEXVAR)
         fOut(NDX(opp(2),i,j,k)) = ftmp2(INDEXVAR)
         fOut(NDX(opp(3),i,j,k)) = ftmp3(INDEXVAR)
         fOut(NDX(opp(4),i,j,k)) = ftmp4(INDEXVAR)
         fOut(NDX(opp(5),i,j,k)) = ftmp5(INDEXVAR)
         fOut(NDX(opp(6),i,j,k)) = ftmp6(INDEXVAR)
         fOut(NDX(opp(7),i,j,k)) = ftmp7(INDEXVAR)
#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
#endif
         fOut(NDX(opp(8),i,j,k)) = ftmp8(INDEXVAR)
         fOut(NDX(opp(9),i,j,k)) = ftmp9(INDEXVAR)
         fOut(NDX(opp(10),i,j,k)) = ftmp10(INDEXVAR)
         fOut(NDX(opp(11),i,j,k)) = ftmp11(INDEXVAR)
         fOut(NDX(opp(12),i,j,k)) = ftmp12(INDEXVAR)
         fOut(NDX(opp(13),i,j,k)) = ftmp13(INDEXVAR)
#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
#endif
         fOut(NDX(opp(14),i,j,k)) = ftmp14(INDEXVAR)
         fOut(NDX(opp(15),i,j,k)) = ftmp15(INDEXVAR)
         fOut(NDX(opp(16),i,j,k)) = ftmp16(INDEXVAR)
         fOut(NDX(opp(17),i,j,k)) = ftmp17(INDEXVAR)
         fOut(NDX(opp(18),i,j,k)) = ftmp18(INDEXVAR)
         fOut(NDX(opp(19),i,j,k)) = ftmp19(INDEXVAR)
#endif
      endif ! not wall
            end do
         end do
      end do
      
      call cpu_time_measure(meas%tEnd_comm)

      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration


   end subroutine stream_collide_bgk 
  !------------------------------------------------------------------------

#endif /* No TRT */





!------------------------------------------------------------------------
subroutine bounceback(lb_dom,fOut,fIn,meas,s_par,prc) !
! calculate bounce-back routine for the treatment of wall boundaries
!
   type(sim_parameter),intent(inout)  :: s_par
   type(mpl_var),intent(inout)  :: prc
   type(lb_block),intent(inout)  :: lb_dom
   real(R8B), intent(in)     :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
   real(R8B), intent(inout)  :: fOut(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
   type(measure )                :: meas  
   integer                 :: ii
   integer                 :: cdx(3)

   call cpu_time_measure(meas%tSt_comm) 

   do ii=1,lb_dom%nobs
      cdx(1:3) = lb_dom%obs(ii,1:3)
      fOut(NDX(1,cdx(1),cdx(2),cdx(3))) = fIn(NDX(opp(1),cdx(1),cdx(2),cdx(3)))
      fOut(NDX(2,cdx(1),cdx(2),cdx(3))) = fIn(NDX(opp(2),cdx(1),cdx(2),cdx(3)))
      fOut(NDX(3,cdx(1),cdx(2),cdx(3))) = fIn(NDX(opp(3),cdx(1),cdx(2),cdx(3)))
      fOut(NDX(4,cdx(1),cdx(2),cdx(3))) = fIn(NDX(opp(4),cdx(1),cdx(2),cdx(3)))
      fOut(NDX(5,cdx(1),cdx(2),cdx(3))) = fIn(NDX(opp(5),cdx(1),cdx(2),cdx(3)))
      fOut(NDX(6,cdx(1),cdx(2),cdx(3))) = fIn(NDX(opp(6),cdx(1),cdx(2),cdx(3)))
      fOut(NDX(7,cdx(1),cdx(2),cdx(3))) = fIn(NDX(opp(7),cdx(1),cdx(2),cdx(3)))
      fOut(NDX(8,cdx(1),cdx(2),cdx(3))) = fIn(NDX(opp(8),cdx(1),cdx(2),cdx(3)))
      fOut(NDX(9,cdx(1),cdx(2),cdx(3))) = fIn(NDX(opp(9),cdx(1),cdx(2),cdx(3)))
#ifdef D3Q19
      fOut(NDX(10,cdx(1),cdx(2),cdx(3))) =fIn(NDX(opp(10),cdx(1),cdx(2),cdx(3))) 
      fOut(NDX(11,cdx(1),cdx(2),cdx(3))) =fIn(NDX(opp(11),cdx(1),cdx(2),cdx(3))) 
      fOut(NDX(12,cdx(1),cdx(2),cdx(3))) =fIn(NDX(opp(12),cdx(1),cdx(2),cdx(3))) 
      fOut(NDX(13,cdx(1),cdx(2),cdx(3))) =fIn(NDX(opp(13),cdx(1),cdx(2),cdx(3))) 
      fOut(NDX(14,cdx(1),cdx(2),cdx(3))) =fIn(NDX(opp(14),cdx(1),cdx(2),cdx(3))) 
      fOut(NDX(15,cdx(1),cdx(2),cdx(3))) =fIn(NDX(opp(15),cdx(1),cdx(2),cdx(3))) 
      fOut(NDX(16,cdx(1),cdx(2),cdx(3))) =fIn(NDX(opp(16),cdx(1),cdx(2),cdx(3))) 
      fOut(NDX(17,cdx(1),cdx(2),cdx(3))) =fIn(NDX(opp(17),cdx(1),cdx(2),cdx(3))) 
      fOut(NDX(18,cdx(1),cdx(2),cdx(3))) =fIn(NDX(opp(18),cdx(1),cdx(2),cdx(3))) 
      fOut(NDX(19,cdx(1),cdx(2),cdx(3))) =fIn(NDX(opp(19),cdx(1),cdx(2),cdx(3))) 
#endif

   end do

      call cpu_time_measure(meas%tEnd_comm)
      meas%bncb_duration = meas%tEnd_comm - meas%tSt_comm + meas%bncb_duration

end subroutine bounceback
!------------------------------------------------------------------------







end module lbm_functions


