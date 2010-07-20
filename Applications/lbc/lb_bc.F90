module lb_bc
   use mpl_lib
   use lbmodel
   use function_set
   use timing
   use lb_init
   use lb_geo 
#include "include/replace.h"
contains





   subroutine set_periodic(lb_dom,fIn,s_par,prc,meas)
   !------------------------------------------------------------------------
   !
   ! Treat periodic boundaries according to init_vector
   ! Tested and ok
   !
      implicit none
      type(lb_block)       :: lb_dom
      type(sim_parameter)  :: s_par 
      type(mpl_var)        :: prc
      type(measure)        :: meas
      integer              :: ii,ll
      real(R8B),dimension(NDX(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
      real(R8B)    :: diff


      if(s_par%periodic(1).eqv..true.) then
      if(prc%size > 1) then; write(*,*) "Periodic boundaries not yet implemented for np>1"; stop; endif 
      !write(*,*) "x periodic"
      fIn(NDX(:,1,:,:)) = fIn(NDX(:,lb_dom%lx(1)-1,:,:)) 
      fIn(NDX(:,lb_dom%lx(1),:,:)) = fIn(NDX(:,2,:,:)) 
      endif      

      if(s_par%periodic(2).eqv..true.) then
      if(prc%size > 1) then; write(*,*) "Periodic boundaries not yet implemented for np>1"; stop; endif 
      !write(*,*) "y periodic"
      fIn(NDX(:,:,1,:)) = fIn(NDX(:,:,lb_dom%lx(2)-1,:)) 
      fIn(NDX(:,:,lb_dom%lx(2),:)) = fIn(NDX(:,:,2,:)) 
      endif     

      if(s_par%periodic(3).eqv..true.) then
      if(prc%size > 1) then; write(*,*) "Periodic boundaries not yet implemented for np>1"; stop; endif 
      !write(*,*) "z periodic"
      fIn(NDX(:,:,:,1)) = fIn(NDX(:,:,:,lb_dom%lx(3)-1)) 
      fIn(NDX(:,:,:,lb_dom%lx(3))) = fIn(NDX(:,:,:,2)) 
      endif     

   end subroutine set_periodic
   !------------------------------------------------------------------------






   subroutine set_nrbc_antibb(lb_dom,fOut,fIn,s_par,prc,meas)
   !------------------------------------------------------------------------
   !
   ! fOut Propagated densities
   ! fIn  Collided densities
   !
      implicit none
      type(lb_block)       :: lb_dom
      type(sim_parameter)  :: s_par 
      type(mpl_var)        :: prc
      type(measure)        :: meas
      integer              :: i,j,k,l,jj
      integer,dimension(-1:3) :: xc,yc,zc
      integer              :: border,ind 
      real(R8B),dimension(NDX(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn,fOut
      real(R8B)  :: dist(2),m,rho0,u0(3),rhod,ud(3),rhod1,ud1(3),rhod2,ud2(3),rhod3,ud3(3)
      real(R8B),dimension(0:3)    :: c1,c2,c3,c4,c5
      real(R8B),dimension(0:3)    :: rhoc
      real(R8B),dimension(3,0:3)  :: uc
      real(R8B)  :: ftmp(nnod)
      real(R8B)  :: ulodi(3),rholodi
      real(R8B)  :: cs_mod, K_mod,kappa 
      real(R8B)  :: usq,usqc
      real(R8B)  :: drhodx,dudx,dvdx,dwdx 
      real(R8B)  :: snu,feqplus(nnod)
      real(R8B)  :: L1,L2,L3,L4,L5 
      real(R8B)  :: Ma,swap,uw,sigma 

      intent(in)    :: fIn
      intent(inout) :: fOut
     
      ! find the neighbors for the nrbc nodes
      if(gtstep_cur == 1 .and. lb_dom%nnr_wall > 0) then
         call nrbc_find_neighbors(lb_dom,s_par,prc)
      endif

      ! Izquierdo constants
      kappa = 1.0d0
      K_mod=0.05d0
      cs_mod = cs*sqrt(kappa)

      ! Do the loop over all nrbc nodes
      do  i=1,lb_dom%nnr_wall  

         ! get coordinates
         xc(0) =lb_dom%obs_nrwall(i,1) 
         yc(0) =lb_dom%obs_nrwall(i,2) 
         zc(0) =lb_dom%obs_nrwall(i,3) 
         xc(1) =lb_dom%obs_nrwall(i,4) 
         yc(1) =lb_dom%obs_nrwall(i,5) 
         zc(1) =lb_dom%obs_nrwall(i,6) 
         xc(2) =lb_dom%obs_nrwall(i,7) 
         yc(2) =lb_dom%obs_nrwall(i,8) 
         zc(2) =lb_dom%obs_nrwall(i,9) 

         ! get current macroscopic variables
        rholodi   = lb_dom%nrwall_prev(i,0) 
        ulodi(1:3)= lb_dom%nrwall_prev(i,1:3)

         call pdf_to_macro(fIn(NDX(:,xc(0),yc(0),zc(0))),rhoc(0),uc(1:3,0))
         call pdf_to_macro(fIn(NDX(:,xc(1),yc(1),zc(1))),rhoc(1),uc(1:3,1))
         call pdf_to_macro(fIn(NDX(:,xc(2),yc(2),zc(2))),rhoc(2),uc(1:3,2))
         ! set zero values
         rho0          = lb_dom%nrwall_0val(i,0) 
         u0(1:3)       = lb_dom%nrwall_0val(i,1:3)

         ! check on which border the current node lies.
         ! the bits correspond to 
         !   1: right border  x+ !   2: left border   x-
         !   3: top border    y+ !   4: bottom border y-
         !   5: front border  z+ !   6: back  border  z-

         if(btest(lb_dom%obs_nrwall(i,0),0)) then
            border = inlet
         else
            border = outlet 
         endif

#ifdef CHECK_FOR_INOUT
         ! determine, if inlet or outlet
         if(btest(lb_dom%obs_nrwall(i,0),dWEST) .or.btest(lb_dom%obs_nrwall(i,0),dEAST) ) then
            if(btest(lb_dom%obs_nrwall(i,0),dEAST)) then
               swap = +1
            else
               swap = -1
            endif
            if(swap*uc(1,0) .ge. 0) then
               border = outlet
            else
               border = inlet
            endif
         else ! y / z border
            if(btest(lb_dom%obs_nrwall(i,0),dNORTH) .or.btest(lb_dom%obs_nrwall(i,0),dSOUTH) ) then
               if(btest(lb_dom%obs_nrwall(i,0),dNORTH)) then
                  swap = +1
               else
                  swap = -1
               endif
               if(swap*uc(2,0) .ge. 0) then
                  border = outlet
               else
                  border = inlet
               endif
            else ! z border 
            endif
         endif
#endif
         if(border == inlet) then
            sigma = 1.0d0
         else
            sigma = 0.0001d0
         endif


         ! --- left and right border  x+ x- 
         !     inlet (or outlet), upstream extrapolation
         if(btest(lb_dom%obs_nrwall(i,0),dWEST) .or.btest(lb_dom%obs_nrwall(i,0),dEAST) ) then
            if(btest(lb_dom%obs_nrwall(i,0),dEAST) ) then  ! x+
              swap = +1.d0
            else
              swap = -1.d0
            endif 
            Ma   = abs(uc(1,0))/cs

          K_mod = swap*sigma*(1.d0 - Ma*Ma)*cs/(s_par%gx(1))    !lodi_length= distance outlet-inlet
!            write(26,*) xc(0),yc(0),zc(0),K_mod
          if(border == inlet) then
              swap = -1.d0*swap
          endif
          ! calculate derivatives
               drhodx =  swap/3.d0*(8.d0*rholodi  - 9.d0*rhoc(1) + rhoc(2))
               dudx   =  swap/3.d0*(8.d0*ulodi(1) - 9.d0*uc(1,1) + uc(1,2))
               dvdx   =  swap/3.d0*(8.d0*ulodi(2) - 9.d0*uc(2,1) + uc(2,2))
               dwdx   =  swap/3.d0*(8.d0*ulodi(3) - 9.d0*uc(3,1) + uc(3,2))
          if(border == inlet) then
               uw = 1.5d0*uc(1,1) - 0.5d0*uc(1,2)
               L1 = (ulodi(1)+swap*cs_mod)*(cs_mod*cs_mod*drhodx/kappa + rholodi*swap*cs_mod*dudx)
               L2 = 0.d0 
               L3 = 0.d0
               L4 = (ulodi(1))*cs_mod*cs_mod*drhodx*(1.-1./kappa)
               L5 = K_mod*(uw-u0(1)) ! Model outgoing wave
          else 
               ! wave amplitudes
               L1 = K_mod*cs_mod*cs_mod*(rholodi-rho0)
               L2 = (ulodi(1))*dvdx
               L3 = (ulodi(1))*dwdx
               L4 = (ulodi(1))*cs_mod*cs_mod*drhodx*(1.-1./kappa)
               L5 = (ulodi(1)+swap*cs_mod)*(cs_mod*cs_mod*drhodx/kappa + rholodi*swap*cs_mod*dudx)
          endif !x+ in / out

               ! LODI Equations
               rholodi  = rholodi  - (L4+0.5d0*(L5+L1))/cs_mod/cs_mod
               ulodi(1) = ulodi(1) - ((L5-L1))/(rho0*swap*cs_mod)*0.5d0
               ulodi(2) = ulodi(2) - L2 
               ulodi(3) = ulodi(3) - L3 

         lb_dom%nrwall_prev(i,0)      =  rholodi       
         lb_dom%nrwall_prev(i,1:NDIM) =  ulodi(1:NDIM) 
         do ind=1,nnod
             fOut(NDX(ind,xc(0),yc(0),zc(0))) = lb_feq(ind,rholodi,ulodi,rho0) 
         enddo

         else  ! check for top / down up/downstream

         ! --- bottom border y-
         !     inlet (or outlet), downstream extrapolation

         if(btest(lb_dom%obs_nrwall(i,0),dNORTH) .or.btest(lb_dom%obs_nrwall(i,0),dSOUTH) ) then
            if(btest(lb_dom%obs_nrwall(i,0),dNORTH) ) then  ! y+
              swap = +1.d0
            else
              swap = -1.d0
            endif 
            Ma   = abs(uc(2,0))/cs

            K_mod = swap*sigma*(1.d0 - Ma*Ma)*cs/(s_par%gx(2))  !lodi_length= distance outlet-inlet
          if(border == inlet) then
              swap = -1.d0*swap
          endif
         ! Derivatives
               drhodx =  swap/3.d0*(8.d0*rholodi  - 9.d0*rhoc(1) + rhoc(2))
               dudx   =  swap/3.d0*(8.d0*ulodi(1) - 9.d0*uc(1,1) + uc(1,2))
               dvdx   =  swap/3.d0*(8.d0*ulodi(2) - 9.d0*uc(2,1) + uc(2,2))
               dwdx   =  swap/3.d0*(8.d0*ulodi(3) - 9.d0*uc(3,1) + uc(3,2))

          if(border == inlet) then
               uw = 1.5d0*uc(2,1) - 0.5d0*uc(2,2)
               L1 = (ulodi(2)+swap*cs_mod)*(cs_mod*cs_mod*drhodx/kappa + rholodi*swap*cs_mod*dvdx)
               L2 = 0.d0 
               L3 = 0.d0
               L4 = (ulodi(2))*cs_mod*cs_mod*drhodx*(1.-1./kappa)
               L5 = K_mod*(uw-u0(2)) ! Model outgoing wave

          else 
               ! wave amplitudes
               L1 = K_mod*cs_mod*cs_mod*(rholodi-rho0)
               L2 = (ulodi(2))*dudx
               L3 = (ulodi(2))*dwdx
               L4 = (ulodi(2))*cs_mod*cs_mod*drhodx*(1.d0-1.d0/kappa)
               L5 = (ulodi(2)+swap*cs_mod)*(cs_mod*cs_mod*drhodx/kappa + rholodi*swap*cs_mod*dvdx)
          endif

         ! LODI Equations
         rholodi  = rholodi - (L4+0.5d0*(L5+L1))/cs_mod/cs_mod
         ulodi(1) = ulodi(1) - L2 
         ulodi(2) = ulodi(2) - ((L5-L1))/(rho0*swap*cs_mod)*0.5d0
         ulodi(3) = ulodi(3)  - L3 
         lb_dom%nrwall_prev(i,0)      =  rholodi       
         lb_dom%nrwall_prev(i,1:NDIM) =  ulodi(1:NDIM) 
               

         do ind=1,nnod
             fOut(NDX(ind,xc(0),yc(0),zc(0))) = lb_feq(ind,rholodi,ulodi,rho0) 
         enddo

            else ! if third dimension has to be extrapolated
#ifdef D3Q19

         !z-
         if(btest(lb_dom%obs_nrwall(i,0),6)) then
          if(border == inlet) then
          else  ! z- outlet
          endif  ! z- in/out

         ! --- front border !z+ 
         !     inlet (or outlet), downstream extrapolation

         elseif(btest(lb_dom%obs_nrwall(i,0),5)) then 
          if(border == inlet) then
          else ! z+ outlet
            
          endif !z- in / out
          endif 
#endif
         end if


         end if

if(s_par%initial) then
         do ind=1,nnod
             fOut(NDX(ind,xc(0),yc(0),zc(0))) = lb_feq(ind,rho0,u0,rho0) 
         enddo
endif

      end do

   end subroutine set_nrbc_antibb
   !------------------------------------------------------------------------







   subroutine set_nrbc(lb_dom,fIn,s_par,prc,meas)
   !------------------------------------------------------------------------
   !
   ! Non-reflecting boundary implementation
   ! Based on work by Giles 
   ! set boundaries as in Andreas Babucke's diss
   ! inlet: set entropic, vortictiy and incoming wave characteristics to 0
   ! outlet: set incoming wave characteristic to 0
! Literatur:   "Non-Reflecting Boundary Conditions for the Euler Equations", Michael Giles,  !
!      CFDL Report 88-1, MIT Dept. of Aero. and Astro., 1988            !
! Download:     http://web.comlab.ox.ac.uk/oucl/work/mike.giles/psfiles/bcs.ps.gz      !
!                                   !
! Bezeichnung der charakteristischen Variablen im Code:                 !
!  c1 Entropie Stoerung                      !
!  c2 Wirbelstaerke Stoerung (in z)                   !
!  c3 Wirbelstaerke Stoerung durch w (in y oder x, je nach Rand)        !
!  c4 Schall Stoerung stromauf                     !
!  c5 Scahll Stoerung stromab                      !
   ! 
   ! this routine should be rechecked. I'm not sure, why there is a problem with 
   ! non-0 velocites imposed on the boundaries
   !
   !
      implicit none
      type(lb_block)       :: lb_dom
      type(sim_parameter)  :: s_par 
      type(mpl_var)        :: prc
      type(measure)        :: meas
      integer              :: i,j,k,l,jj
      integer,dimension(-1:3) :: xc,yc,zc
      integer              :: border 
      real(R8B),dimension(NDX(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
      real(R8B)  :: dist(2),m,rho0,u0(3),rhod,ud(3),rhod1,ud1(3),rhod2,ud2(3),rhod3,ud3(3)
      real(R8B),dimension(0:3)  :: c1,c2,c3,c4,c5,rhoc
      real(R8B),dimension(3,0:3)  :: uc
      real(R8B)  :: factoru(1:3),damping 
      real(R8B)  :: usq,cs2inv,t2cs4inv,t2cs2inv

!------------------
! If NRBC_LINEAR_EXTRAPOLATION is defined, use linear instead of parabolic extrapolation
!#!define NRBC_LINEAR_EXTRAPOLATION
! Note: 
! parabolic interpolation is supposed to work better (higher accuracy) 
! but seems to be instable with the corotating vortex pair problem. BUG???!
! Use linear interpolation until bug is found

! if(s_par%initial ) then ! Not doing anything, because during initialization, no nrbc treatment is required ! endif

if(gtstep_cur == 1) then
      if (lb_dom%nnr_wall > 0) then
         call nrbc_find_neighbors(lb_dom,s_par,prc)
      end if
endif

      do  i=1,lb_dom%nnr_wall  

         xc(-1)=lb_dom%obs_nrwall(i,1)  
         yc(-1)=lb_dom%obs_nrwall(i,2)  
         zc(-1)=lb_dom%obs_nrwall(i,3)  
         xc(0) =lb_dom%obs_nrwall(i,4)  
         yc(0) =lb_dom%obs_nrwall(i,5)  
         zc(0) =lb_dom%obs_nrwall(i,6)  
         xc(1) =lb_dom%obs_nrwall(i,7)  
         yc(1) =lb_dom%obs_nrwall(i,8)  
         zc(1) =lb_dom%obs_nrwall(i,9)  
         xc(2) =lb_dom%obs_nrwall(i,10) 
         yc(2) =lb_dom%obs_nrwall(i,11) 
         zc(2) =lb_dom%obs_nrwall(i,12) 
         xc(3) =lb_dom%obs_nrwall(i,13) 
         yc(3) =lb_dom%obs_nrwall(i,14) 
         zc(3) =lb_dom%obs_nrwall(i,15) 

!         xc(0) =lb_dom%obs_nrwall(i,1) 
!         yc(0) =lb_dom%obs_nrwall(i,2) 
!         zc(0) =lb_dom%obs_nrwall(i,3) 
!         xc(1) =lb_dom%obs_nrwall(i,4) 
!         yc(1) =lb_dom%obs_nrwall(i,5) 
!         zc(1) =lb_dom%obs_nrwall(i,6) 
!         xc(2) =lb_dom%obs_nrwall(i,7) 
!         yc(2) =lb_dom%obs_nrwall(i,8) 
!         zc(2) =lb_dom%obs_nrwall(i,9) 
!         xc(3) =lb_dom%obs_nrwall(i,10)
!         yc(3) =lb_dom%obs_nrwall(i,11)
!         zc(3) =lb_dom%obs_nrwall(i,12)


         ! calculate current macroscopic variables
         do jj = 0,3
#if defined D2Q9
        rhoc(jj) = fIn(NDX(1,xc(jj),yc(jj),zc(jj)))  + fIn(NDX(2,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(3,xc(jj),yc(jj),zc(jj)))  + fIn(NDX(4,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(5,xc(jj),yc(jj),zc(jj)))  +                                         &
                   fIn(NDX(6,xc(jj),yc(jj),zc(jj)))  + fIn(NDX(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(8,xc(jj),yc(jj),zc(jj)))  + fIn(NDX(9,xc(jj),yc(jj),zc(jj))) 

       uc(1,jj) = (fIn(NDX(2,xc(jj),yc(jj),zc(jj)))  - fIn(NDX(4,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(6,xc(jj),yc(jj),zc(jj)))  - fIn(NDX(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(9,xc(jj),yc(jj),zc(jj)))  - fIn(NDX(8,xc(jj),yc(jj),zc(jj)))    & 
                       )/rhoc(jj)    

       uc(2,jj) = (fIn(NDX(3,xc(jj),yc(jj),zc(jj)))  - fIn(NDX(5,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(6,xc(jj),yc(jj),zc(jj)))  - fIn(NDX(9,xc(jj),yc(jj),zc(jj)))  - &
                   fIn(NDX(7,xc(jj),yc(jj),zc(jj)))  + fIn(NDX(8,xc(jj),yc(jj),zc(jj)))    & 
                      ) /rhoc(jj)  
       uc(3,jj) = 0.0d0

#elif defined D3Q19 /**/
        rhoc(jj) = fIn(NDX(1,xc(jj),yc(jj),zc(jj)))  + fIn(NDX(2,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(3,xc(jj),yc(jj),zc(jj)))  + fIn(NDX(4,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(5,xc(jj),yc(jj),zc(jj)))  +                                         &
                   fIn(NDX(6,xc(jj),yc(jj),zc(jj)))  + fIn(NDX(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(8,xc(jj),yc(jj),zc(jj)))  + fIn(NDX(9,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(10,xc(jj),yc(jj),zc(jj))) +                                         & 
                   fIn(NDX(11,xc(jj),yc(jj),zc(jj))) + fIn(NDX(12,xc(jj),yc(jj),zc(jj))) + &
                   fIn(NDX(13,xc(jj),yc(jj),zc(jj))) + fIn(NDX(14,xc(jj),yc(jj),zc(jj))) + &
                   fIn(NDX(15,xc(jj),yc(jj),zc(jj))) +                                         & 
                   fIn(NDX(16,xc(jj),yc(jj),zc(jj))) + fIn(NDX(17,xc(jj),yc(jj),zc(jj))) + &
                   fIn(NDX(18,xc(jj),yc(jj),zc(jj))) + fIn(NDX(19,xc(jj),yc(jj),zc(jj)))

       uc(1,jj) = (fIn(NDX(2,xc(jj),yc(jj),zc(jj)))  - fIn(NDX(3,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(8,xc(jj),yc(jj),zc(jj)))  - fIn(NDX(9,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(10,xc(jj),yc(jj),zc(jj))) - fIn(NDX(11,xc(jj),yc(jj),zc(jj))) + & 
                   fIn(NDX(12,xc(jj),yc(jj),zc(jj))) - fIn(NDX(13,xc(jj),yc(jj),zc(jj))) + &
                   fIn(NDX(14,xc(jj),yc(jj),zc(jj))) - fIn(NDX(15,xc(jj),yc(jj),zc(jj)))   & 
                       )/rhoc(jj)    

       uc(2,jj) = (fIn(NDX(4,xc(jj),yc(jj),zc(jj)))  - fIn(NDX(5,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(8,xc(jj),yc(jj),zc(jj)))  - fIn(NDX(9,xc(jj),yc(jj),zc(jj)))  - &
                   fIn(NDX(10,xc(jj),yc(jj),zc(jj))) + fIn(NDX(11,xc(jj),yc(jj),zc(jj))) + & 
                   fIn(NDX(16,xc(jj),yc(jj),zc(jj))) - fIn(NDX(17,xc(jj),yc(jj),zc(jj))) + &
                   fIn(NDX(18,xc(jj),yc(jj),zc(jj))) - fIn(NDX(19,xc(jj),yc(jj),zc(jj)))   & 
                      ) /rhoc(jj)     

       uc(3,jj) = (fIn(NDX(6,xc(jj),yc(jj),zc(jj)))  - fIn(NDX(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(NDX(12,xc(jj),yc(jj),zc(jj))) - fIn(NDX(13,xc(jj),yc(jj),zc(jj))) - &
                   fIn(NDX(14,xc(jj),yc(jj),zc(jj)))  + & 
                   fIn(NDX(15,xc(jj),yc(jj),zc(jj))) + fIn(NDX(16,xc(jj),yc(jj),zc(jj))) - &
                   fIn(NDX(17,xc(jj),yc(jj),zc(jj))) - fIn(NDX(18,xc(jj),yc(jj),zc(jj))) + &
                   fIn(NDX(19,xc(jj),yc(jj),zc(jj))))/rhoc(jj)
#endif
         enddo
         ! Calculate the distance of the nodes from the neighbor nodes
         dist(1) = sqrt(real(yc(1) - yc(0))**2 + real(xc(1)-xc(0))**2 + real(zc(1)-zc(0))**2)
         dist(2) = sqrt(real(yc(2) - yc(0))**2 + real(xc(2)-xc(0))**2 + real(zc(2)-zc(0))**2)

         ! set zero values
         rho0       = lb_dom%nrwall_0val(i,0) 
         u0(1:NDIM) = lb_dom%nrwall_0val(i,1:NDIM)

         ud(1) = uc(1,0) - u0(1)
         ud(2) = uc(2,0) - u0(2)
#ifdef D3Q19
         ud(3) = uc(3,0) - u0(3)
#else /* D3Q19 */
         ud(3) = 0.d0 
#endif /* D3Q19 */
         rhod = rhoc(0) - rho0


         ! calculate the fluctuating part of the first neighbor
         ud1(1) = uc(1,1) - u0(1)
         ud1(2) = uc(2,1) - u0(2)
#ifdef D3Q19
         ud1(3) = uc(3,1) - u0(3)
#else /* D3Q19 */
         ud1(3) = 0.d0 
#endif /* D3Q19 */
         rhod1 = rhoc(1) - rho0

         ! calculate the fluctuating part of the second neighbor
         ud2(1) = uc(1,2) - u0(1)
         ud2(2) = uc(2,2) - u0(2)
#ifdef D3Q19             
         ud2(3) = uc(3,2) - u0(3)
#else /* D3Q19 */
         ud2(3) = 0.d0 
#endif /* D3Q19 */
         rhod2 = rhoc(2) - rho0

         ! calculate the fluctuating part of the third neighbor
         ud3(1) = uc(1,3) - u0(1)
         ud3(2) = uc(2,3) - u0(2)
#ifdef D3Q19             
         ud3(3) = uc(3,3) - u0(3)
#else /* D3Q19 */
         ud3(3) = 0.d0 
#endif /* D3Q19 */
         rhod3 = rhoc(3) - rho0


         ! check on which border the current node lies.
         ! the bits correspond to 
         !   1: right border  x+
         !   2: left border   x-
         !   3: top border    y+
         !   4: bottom border y-
         !   5: front border  z+
         !   6: back  border  z-

         if(btest(lb_dom%obs_nrwall(i,0),0)) then
            border = inlet
         else          
            border = outlet
         endif

         ! --- left  border  x- 
         !     inlet (or outlet), upstream extrapolation

         if(btest(lb_dom%obs_nrwall(i,0),3)) then
          if(border == inlet) then
            c5(1) = -ud1(1)*rho0*cs + 1./3.*rhod1
            c5(2) = -ud2(1)*rho0*cs + 1./3.*rhod2
            c5(3) = -ud3(1)*rho0*cs + 1./3.*rhod3
#ifdef NRBC_LINEAR_EXTRAPOLATION
            ! linear extrapolation     
            m = (c5(2) - c5(1))/(dist(2)-dist(1))
            c5(0) = m*0 + (c5(1) - dist(1)*m)    
#else
           ! parabolic extrapolation
            c5(0) = 3._R8B*c5(1) - 3._R8B*c5(2) + c5(3)
#endif
            c4(0) = 0.d0
            c2=0.d0            
            c3=0.d0            
            factoru(1)=c4(0)-c5(0)
            factoru(2)=0.d0
            factoru(3)=0.d0
          else  ! x- outlet
               c4(0) = 0.d0     ! upstream perturbation is zero

               c2(1) = ud1(2)*rho0*cs
               c2(2) = ud2(2)*rho0*cs
               c2(3) = ud3(2)*rho0*cs
               c3(1) = ud1(3)*rho0*cs
               c3(2) = ud2(3)*rho0*cs
               c3(3) = ud3(3)*rho0*cs
               c5(1) = -ud1(1)*rho0*cs + 1./3.*rhod1
               c5(2) = -ud2(1)*rho0*cs + 1./3.*rhod2
               c5(3) = -ud3(1)*rho0*cs + 1./3.*rhod3

#ifdef NRBC_LINEAR_EXTRAPOLATION
               ! linear extrapolation
               m = (c2(2) - c2(1))/(dist(2)-dist(1))
               c2(0) = m*0 + (c2(1) - dist(1)*m)    
               m = (c3(2) - c3(1))/(dist(2)-dist(1))
               c3(0) = m*0 + (c3(1) - dist(1)*m)    
               m = (c5(2) - c5(1))/(dist(2)-dist(1))
               c5(0) = m*0 + (c5(1) - dist(1)*m)    
#else
               ! parabolic extrapolation
               c2(0) = 3._R8B*c2(1) - 3._R8B*c2(2) + c2(3)
               c3(0) = 3._R8B*c3(1) - 3._R8B*c3(2) + c3(3)
               c5(0) = 3._R8B*c5(1) - 3._R8B*c5(2) + c5(3)
#endif
               factoru(1)=c4(0)-c5(0)
               factoru(2)=c2(0)
               factoru(3)=c3(0)
          endif  ! x- in/out

         ! --- right border !x+ 
         !     inlet (or outlet), downstream extrapolation

         elseif(btest(lb_dom%obs_nrwall(i,0),1)) then 
          if(border == inlet) then
            c4(1) = ud1(1)*rho0*cs + 1./3.*rhod1
            c4(2) = ud2(1)*rho0*cs + 1./3.*rhod2
            c4(3) = ud3(1)*rho0*cs + 1./3.*rhod3
#ifdef NRBC_LINEAR_EXTRAPOLATION
            ! linear extrapolation
            m = (c4(2) - c4(1))/(dist(2)-dist(1))
            c4(0) = m*0 + (c4(1) - dist(1)*m) 
#else
            ! parabolic extrapolation
            c4(0) = 3._R8B*c4(1) - 3._R8B*c4(2) + c4(3)
#endif
            c5(0) = 0.d0
            c1=0.d0            
            c2=0.d0            
            factoru(1)=c4(0)-c5(0)
            factoru(2)=0.d0
            factoru(3)=0.d0
          else ! x- outlet

               c5(0) = 0.d0     ! upstream perturbation is zero

               c2(1) = ud1(2)*rho0*cs
               c2(2) = ud2(2)*rho0*cs
               c2(3) = ud3(2)*rho0*cs
               c3(1) = ud1(3)*rho0*cs
               c3(2) = ud2(3)*rho0*cs
               c3(3) = ud3(3)*rho0*cs
               c4(1) = ud1(1)*rho0*cs + 1./3.*rhod1
               c4(2) = ud2(1)*rho0*cs + 1./3.*rhod2
               c4(3) = ud3(1)*rho0*cs + 1./3.*rhod3

#ifdef NRBC_LINEAR_EXTRAPOLATION
               ! linear extrapolation
               m = (c2(2) - c2(1))/(dist(2)-dist(1))
               c2(0) = m*0 + (c2(1) - dist(1)*m)    
               m = (c3(2) - c3(1))/(dist(2)-dist(1))
               c3(0) = m*0 + (c3(1) - dist(1)*m)    
               m = (c4(2) - c4(1))/(dist(2)-dist(1))
               c4(0) = m*0 + (c4(1) - dist(1)*m)    
#else
               ! parabolic extrapolation
               c2(0) = 3._R8B*c2(1) - 3._R8B*c2(2) + c2(3)
               c3(0) = 3._R8B*c3(1) - 3._R8B*c3(2) + c3(3)
               c4(0) = 3._R8B*c4(1) - 3._R8B*c4(2) + c4(3)
#endif
               factoru(1)=c4(0)-c5(0)
               factoru(2)=c2(0)
               factoru(3)=c3(0)
            
          endif !x- in / out

         else  ! check for top / down up/downstream

         ! --- bottom border y-
         !     inlet (or outlet), downstream extrapolation

            if(btest(lb_dom%obs_nrwall(i,0),4)) then 
            if(border == inlet) then ! y- inlet 
               c2=0.d0       ! vorticity fluctuation y is zero        
               c3=0.d0       !                       z 
               c4(0) = 0.d0  ! upstream perturbation is zero
               c5(1) = -ud1(2)*rho0*cs + 1./3.*rhod1
               c5(2) = -ud2(2)*rho0*cs + 1./3.*rhod2
               c5(3) = -ud2(2)*rho0*cs + 1./3.*rhod3
               m = (c5(2) - c5(1))/(dist(2)-dist(1))
#ifdef NRBC_LINEAR_EXTRAPOLATION
               ! linear extrapolation
               c5(0) = m*0 + (c5(1) - dist(1)*m)       
#else
               ! parabolic extrapolation
               c5(0) = 3._R8B*c5(1) - 3._R8B*c5(2) + c5(3)
#endif       
               factoru(1)=0.d0
               factoru(2)=c4(0) - c5(0) 
               factoru(3)=0.d0
            else ! y- outlet 
               c4(0) = 0.d0     ! upstream perturbation is zero

               c2(1) =  ud1(1)*rho0*cs
               c2(2) =  ud2(1)*rho0*cs
               c2(3) =  ud3(1)*rho0*cs
               c3(1) =  ud1(3)*rho0*cs
               c3(2) =  ud2(3)*rho0*cs
               c3(3) =  ud3(3)*rho0*cs
               c5(1) = -ud1(2)*rho0*cs + 1./3.*rhod1
               c5(2) = -ud2(2)*rho0*cs + 1./3.*rhod2
               c5(3) = -ud3(2)*rho0*cs + 1./3.*rhod3

#ifdef NRBC_LINEAR_EXTRAPOLATION
!               m = (c2(2) - c2(1))/(dist(2)-dist(1))
!FIXME               c2(0) = m*0 + (c2(1) - dist(1)*m)    
!               c2(0) = (dist(1)*c2(1)+c2(1)-c2(2))/dist(1)    
               c2(0) = 2.0d0*c2(1)-c2(2)
!               m = (c3(2) - c3(1))/(dist(2)-dist(1))
!FIXME               c3(0) = m*0 + (c3(1) - dist(1)*m)    
!               c3(0) = (dist(1)*c3(1)+c3(1)-c3(2))/dist(1)    
               c3(0) = 2.0d0*c3(1)-c3(2)
!               m = (c5(2) - c5(1))/(dist(2)-dist(1))
!FIXME               c5(0) = m*0 + (c5(1) - dist(1)*m)    
!               c5(0) = (dist(1)*c5(1)+c5(1)-c5(2))/dist(1)    
               c5(0) = 2.0d0*c5(1)-c5(2)
#else
               ! parabolic extrapolation
               c2(0) = 3._R8B*c2(1) - 3._R8B*c2(2) + c2(3)
               c3(0) = 3._R8B*c3(1) - 3._R8B*c3(2) + c3(3)
               c5(0) = 3._R8B*c5(1) - 3._R8B*c5(2) + c5(3)
#endif
               factoru(1)=c2(0)
               factoru(2)=c4(0)-c5(0)
               factoru(3)=c3(0)
            endif ! y- inlet/outlet

         ! --- top border y+
         !     inlet (or outlet), upstream extrapolation

            elseif(btest(lb_dom%obs_nrwall(i,0),2)) then 
            if(border == inlet) then ! y+
               c2(0) = 0.d0          ! vorticity fluctuation y is zero  
               c3(0) = 0.d0          !                       z         
               c5(0) = 0.d0          ! upstream perturbation is zero
               c4(1) = ud1(2)*rho0*cs + 1./3.*rhod1
               c4(2) = ud2(2)*rho0*cs + 1./3.*rhod2
               c4(3) = ud3(2)*rho0*cs + 1./3.*rhod3

#ifdef NRBC_LINEAR_EXTRAPOLATION
               m = (c4(2) - c4(1))/(dist(2)-dist(1))
               c4(0) = m*0 + (c4(1) - dist(1)*m) 
#else
               ! parabolic extrapolation
               c4(0) = 3._R8B*c4(1) - 3._R8B*c4(2) + c4(3)
#endif
               factoru(1)=c2(0)
               factoru(2)=c4(0)-c5(0)
               factoru(3)=c3(0)
            else ! y+ outlet
               c5(0) = 0.d0     ! upstream perturbation is zero

               c2(1) = ud1(1)*rho0*cs
               c2(2) = ud2(1)*rho0*cs
               c2(3) = ud3(1)*rho0*cs
               c3(1) = ud1(3)*rho0*cs
               c3(2) = ud2(3)*rho0*cs
               c3(3) = ud3(3)*rho0*cs
               c4(1) = ud1(2)*rho0*cs + 1./3.*rhod1
               c4(2) = ud2(2)*rho0*cs + 1./3.*rhod2
               c4(3) = ud3(2)*rho0*cs + 1./3.*rhod3

#ifdef NRBC_LINEAR_EXTRAPOLATION
               ! linear extrapolation
               m = (c2(2) - c2(1))/(dist(2)-dist(1))
               c2(0) = m*0 + (c2(1) - dist(1)*m)
               m = (c3(2) - c3(1))/(dist(2)-dist(1))
               c3(0) = m*0 + (c3(1) - dist(1)*m)
               m = (c4(2) - c4(1))/(dist(2)-dist(1))
               c4(0) = m*0 + (c4(1) - dist(1)*m)
#else
               ! parabolic extrapolation
               c2(0) = 3._R8B*c2(1) - 3._R8B*c2(2) + c2(3)
               c3(0) = 3._R8B*c3(1) - 3._R8B*c3(2) + c3(3)
               c4(0) = 3._R8B*c4(1) - 3._R8B*c4(2) + c4(3)
#endif

               factoru(1)=c2(0)
               factoru(2)=c4(0)-c5(0)
               factoru(3)=c3(0)
            endif ! y+ inlet/outlet
            else ! if third dimension has to be extrapolated
#ifdef D3Q19

         !z-
         if(btest(lb_dom%obs_nrwall(i,0),6)) then
          if(border == inlet) then
            c5(1) = -ud1(3)*rho0*cs + 1./3.*rhod1
            c5(2) = -ud2(3)*rho0*cs + 1./3.*rhod2
            c5(3) = -ud3(3)*rho0*cs + 1./3.*rhod3
#ifdef NRBC_LINEAR_EXTRAPOLATION
            ! linear extrapolation     
            m = (c5(2) - c5(1))/(dist(2)-dist(1))
            c5(0) = m*0 + (c5(1) - dist(1)*m)    
#else
           ! parabolic extrapolation
            c5(0) = 3._R8B*c5(1) - 3._R8B*c5(2) + c5(3)
#endif
            c4(0) = 0.d0
            factoru(1)=0.d0
            factoru(2)=0.d0
            factoru(3)=c4(0)-c5(0)
          else  ! z- outlet
               c4(0) = 0.d0     ! upstream perturbation is zero

               c2(1) = ud1(1)*rho0*cs
               c2(2) = ud2(1)*rho0*cs
               c2(3) = ud3(1)*rho0*cs
               c3(1) = ud1(2)*rho0*cs
               c3(2) = ud2(2)*rho0*cs
               c3(3) = ud3(2)*rho0*cs
               c5(1) = -ud1(3)*rho0*cs + 1./3.*rhod1
               c5(2) = -ud2(3)*rho0*cs + 1./3.*rhod2
               c5(3) = -ud3(3)*rho0*cs + 1./3.*rhod3

#ifdef NRBC_LINEAR_EXTRAPOLATION
               ! linear extrapolation
               m = (c2(2) - c2(1))/(dist(2)-dist(1))
               c2(0) = m*0 + (c2(1) - dist(1)*m)    
               m = (c3(2) - c3(1))/(dist(2)-dist(1))
               c3(0) = m*0 + (c3(1) - dist(1)*m)    
               m = (c5(2) - c5(1))/(dist(2)-dist(1))
               c5(0) = m*0 + (c5(1) - dist(1)*m)    
#else
               ! parabolic extrapolation
               c2(0) = 3._R8B*c2(1) - 3._R8B*c2(2) + c2(3)
               c3(0) = 3._R8B*c3(1) - 3._R8B*c3(2) + c3(3)
               c5(0) = 3._R8B*c5(1) - 3._R8B*c5(2) + c5(3)
#endif
               factoru(1)=c2(0)
               factoru(2)=c3(0)
               factoru(3)=c4(0)-c5(0)
          endif  ! z- in/out

         ! --- front border !z+ 
         !     inlet (or outlet), downstream extrapolation

         elseif(btest(lb_dom%obs_nrwall(i,0),5)) then 
          if(border == inlet) then
            c4(1) = ud1(3)*rho0*cs + 1./3.*rhod1
            c4(2) = ud2(3)*rho0*cs + 1./3.*rhod2
            c4(3) = ud3(3)*rho0*cs + 1./3.*rhod3
#ifdef NRBC_LINEAR_EXTRAPOLATION
            ! linear extrapolation
            m = (c4(2) - c4(1))/(dist(2)-dist(1))
            c4(0) = m*0 + (c4(1) - dist(1)*m) 
#else
            ! parabolic extrapolation
            c4(0) = 3._R8B*c4(1) - 3._R8B*c4(2) + c4(3)
#endif
            c5(0) = 0.d0
            c1=0.d0            
            c2=0.d0            
            factoru(1)=0.d0
            factoru(2)=0.d0
            factoru(3)=c4(0)-c5(0)

          else ! z+ outlet

               c5(0) = 0.d0     ! upstream perturbation is zero

               c2(1) = ud1(1)*rho0*cs
               c2(2) = ud2(1)*rho0*cs
               c2(3) = ud3(1)*rho0*cs
               c3(1) = ud1(2)*rho0*cs
               c3(2) = ud2(2)*rho0*cs
               c3(3) = ud3(2)*rho0*cs
               c4(1) = ud1(3)*rho0*cs + 1./3.*rhod1
               c4(2) = ud2(3)*rho0*cs + 1./3.*rhod2
               c4(3) = ud3(3)*rho0*cs + 1./3.*rhod3

#ifdef NRBC_LINEAR_EXTRAPOLATION
               ! linear extrapolation
               m = (c2(2) - c2(1))/(dist(2)-dist(1))
               c2(0) = m*0 + (c2(1) - dist(1)*m)    
               m = (c3(2) - c3(1))/(dist(2)-dist(1))
               c3(0) = m*0 + (c3(1) - dist(1)*m)    
               m = (c4(2) - c4(1))/(dist(2)-dist(1))
               c4(0) = m*0 + (c4(1) - dist(1)*m)    
#else
               ! parabolic extrapolation
               c2(0) = 3._R8B*c2(1) - 3._R8B*c2(2) + c2(3)
               c3(0) = 3._R8B*c3(1) - 3._R8B*c3(2) + c3(3)
               c4(0) = 3._R8B*c4(1) - 3._R8B*c4(2) + c4(3)
#endif
               factoru(1)=c2(0)
               factoru(2)=c3(0)
               factoru(3)=c4(0)-c5(0)
            
          endif !z- in / out
          endif 




#endif
            end if


         end if

         ! if the nrbc macr vals should be damped, define
         damping = 1.00

         ! calculate new macroscopic values

         cs2inv  = 1._R8B/cs**2
         t2cs4inv = 1._R8B/(2*cs**4)
         t2cs2inv = 1._R8B/(2*cs**2)

if(s_par%initial) then
   damping = 0.0d0
endif

         uc(1:NDIM,0) = damping* (factoru(1:NDIM))/(2*rho0*cs) + u0(1:NDIM) 
         rhoc(0)      = rho0 + damping*(cs2inv*0.5*(c4(0) + c5(0)))

         ! now calculate the equilibrium distribution of the 
         ! calculated macroscopic values
!#define NRBC_SPEEDUP
#ifndef NRBC_SPEEDUP


         lb_dom%u(NDX(1:NDIM,xc(0),yc(0),zc(0))) = damping* (factoru(1:NDIM))/(2*rho0*cs) + u0(1:NDIM) 
         lb_dom%rho(xc(0),yc(0),zc(0))      = rho0 + damping*(cs2inv*0.5*(c4(0) + c5(0)))
         call calc_fEq(lb_dom,lb_dom%rho,lb_dom%u,fIn,xc(0),xc(0),yc(0) ,yc(0) ,zc(0),zc(0),s_par) 
         
         lb_dom%u(NDX(1:NDIM,xc(-1),yc(-1),zc(-1))) = damping* (factoru(1:NDIM))/(2*rho0*cs) + u0(1:NDIM) 
         lb_dom%rho(xc(-1),yc(-1),zc(-1))      = rho0 + damping*(cs2inv*0.5*(c4(0) + c5(0)))
         call calc_fEq(lb_dom,lb_dom%rho,lb_dom%u,fIn,xc(-1),xc(-1),yc(-1) ,yc(-1) ,zc(-1),zc(-1),s_par) 
 
 

#endif /* NRBC_SPEEDUP*/
#ifdef NRBC_SPEEDUP
         uc(1:NDIM,0) = damping* (factoru(1:NDIM))/(2*rho0*cs) + u0(1:NDIM) 
         rhoc(0)      = rho0 + damping*(cs2inv*0.5*(c4(0) + c5(0)))
         ! calculate  equilibrium distribution
         usq =  (uc(1,0)*uc(1,0) + uc(2,0)*uc(2,0) + uc(3,0)*uc(3,0))*t2cs2inv
#if defined D2Q9
             fIn(NDX(1,xc(0),yc(0),zc(0))) = t(1)*rhoc(0) *(1._R8B - usq)
             fIn(NDX(2,xc(0),yc(0),zc(0))) = t(2)*rhoc(0) *(1._R8B + (uc(1,0))*cs2inv &
            + (uc(1,0))**2*t2cs4inv   - usq)
             fIn(NDX(3,xc(0),yc(0),zc(0))) = t(3)*rhoc(0)*(1._R8B - (-uc(2,0))*cs2inv &
            + (uc(2,0))**2*t2cs4inv   - usq)
             fIn(NDX(4,xc(0),yc(0),zc(0))) = t(4)*rhoc(0)*(1._R8B + (-uc(1,0))*cs2inv &
           + (-uc(1,0))**2*t2cs4inv  - usq)
             fIn(NDX(5,xc(0),yc(0),zc(0))) = t(5)*rhoc(0)*(1._R8B + (-uc(2,0))*cs2inv &
           + (-uc(2,0))**2*t2cs4inv   - usq)
             fIn(NDX(6,xc(0),yc(0),zc(0))) = t(6)*rhoc(0)*(1._R8B  &
            + (uc(1,0) + uc(2,0))*cs2inv &
            + (uc(1,0)+uc(2,0))**2*t2cs4inv  - usq)
             fIn(NDX(7,xc(0),yc(0),zc(0))) = t(7)*rhoc(0)*(1._R8B  &
           + (-uc(1,0) + uc(2,0))*cs2inv &
           + (-uc(1,0)+uc(2,0))**2*t2cs4inv  - usq)
             fIn(NDX(8,xc(0),yc(0),zc(0))) = t(8)*rhoc(0)*(1._R8B  &
           + (-uc(1,0) -uc(2,0))*cs2inv &
           + (-uc(1,0)-uc(2,0))**2*t2cs4inv    - usq)
             fIn(NDX(9,xc(0),yc(0),zc(0))) = t(9)*rhoc(0)*(1._R8B  &
            + (uc(1,0) -uc(2,0))*cs2inv &
            + (uc(1,0)-uc(2,0))**2*t2cs4inv  - usq)
#elif defined D3Q19
           fIn(NDX(1,xc(0),yc(0),zc(0))) = t(1)*rhoc(0)*(1  &
                - usq) 
           fIn(NDX(2,xc(0),yc(0),zc(0))) = t(2)*rhoc(0)*(1._R8B + (uc(1,0))*cs2inv &
                + (uc(1,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(3,xc(0),yc(0),zc(0))) = t(3)*rhoc(0)*(1._R8B + (-uc(1,0))*cs2inv &
                -uc(1,0)**2*t2cs4inv  &
                - usq) 
           fIn(NDX(4,xc(0),yc(0),zc(0))) = t(4)*rhoc(0)*(1._R8B + uc(2,0)*cs2inv &
                + (uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(5,xc(0),yc(0),zc(0))) = t(5)*rhoc(0)*(1._R8B - uc(2,0)*cs2inv &
                + (-uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(6,xc(0),yc(0),zc(0))) = t(6)*rhoc(0)*(1._R8B + (uc(3,0))*cs2inv &
                + (uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(7,xc(0),yc(0),zc(0))) = t(7)*rhoc(0)*(1._R8B + (-uc(3,0))*cs2inv &
                + (-uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(8,xc(0),yc(0),zc(0))) = t(8)*rhoc(0)*(1._R8B  &
                + (uc(1,0) + uc(2,0))*cs2inv &
                + (uc(1,0) + uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(9,xc(0),yc(0),zc(0))) = t(9)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) -uc(2,0))*cs2inv &
                + (-uc(1,0) -uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(10,xc(0),yc(0),zc(0))) = t(10)*rhoc(0)*(1._R8B  &
                + (uc(1,0) -uc(2,0))*cs2inv &
                + (uc(1,0) -uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(11,xc(0),yc(0),zc(0))) = t(11)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) + uc(2,0))*cs2inv &
                + (-uc(1,0) + uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(12,xc(0),yc(0),zc(0))) = t(12)*rhoc(0)*(1._R8B  &
                + (uc(1,0) + uc(3,0))*cs2inv &
                + (uc(1,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(13,xc(0),yc(0),zc(0))) = t(13)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) -uc(3,0))*cs2inv &
                + (-uc(1,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(14,xc(0),yc(0),zc(0))) = t(14)*rhoc(0)*(1._R8B  &
                + (uc(1,0) -uc(3,0))*cs2inv &
                + (uc(1,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(15,xc(0),yc(0),zc(0))) = t(15)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) + uc(3,0))*cs2inv &
                + (-uc(1,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(16,xc(0),yc(0),zc(0))) = t(16)*rhoc(0)*(1._R8B  &
                + (uc(2,0) + uc(3,0))*cs2inv &
                + (uc(2,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(17,xc(0),yc(0),zc(0))) = t(17)*rhoc(0)*(1._R8B  &
                + (-uc(2,0) -uc(3,0))*cs2inv &
                + (-uc(2,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(18,xc(0),yc(0),zc(0))) = t(18)*rhoc(0)*(1._R8B  &
                + (uc(2,0) -uc(3,0))*cs2inv &
                + (uc(2,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(NDX(19,xc(0),yc(0),zc(0))) = t(19)*rhoc(0)*(1._R8B  &
                + (-uc(2,0) + uc(3,0))*cs2inv &
                + (-uc(2,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
#endif /* D3Q19 */
#endif /* NRBC_SPEEDUP */     
      end do

   end subroutine set_nrbc
   !------------------------------------------------------------------------





   subroutine microphone(lb_dom,pos,fIn,s_par,prc,tstep,nr)
   !------------------------------------------------------------------------
   !
   ! Get current density at a certain place 
   ! 

      implicit none
      type(lb_block)       :: lb_dom
      type(sim_parameter)  :: s_par 
      type(mpl_var)        :: prc
      integer              :: tstep,ll,nr
      integer              :: xmin(3),xmax(3) 
      integer,intent(in)   :: pos(3)
      real(R8B)            :: rholoc
      real(R8B),dimension(NDX(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
      character(len=8) :: ch_rank,ch_posx,ch_posy,ch_posz
      logical :: record

      xmin = lb_dom%lx*prc%crd+1
      xmax = lb_dom%lx*(prc%crd+1)

      record = .false.

      if(xmin(1) .le. pos(1) .and. xmax(1) .ge. pos(1) .and. &
         xmin(2) .le. pos(2) .and. xmax(2) .ge. pos(2).and. &
         xmin(3) .le. pos(3) .and. xmax(3) .ge. pos(3)      &
         ) record = .true. 

      if(record .eqv. .true.) then
         rholoc  = 0.0d0
         do ll=1,nnod
            if(record .eqv. .true. ) &
            rholoc = rholoc+fIn(NDX(ll,pos(1)-xmin(1)+1,pos(2)-xmin(2)+1,pos(3)-xmin(3)+1))
         enddo
      if(tstep==1) then
      write(ch_rank,'(i7)') prc%rk
      write(ch_posx,'(i7)') pos(1) 
      write(ch_posy,'(i7)') pos(2)
      write(ch_posz,'(i7)') pos(3)
      ch_rank = adjustl(ch_rank)
      ch_posx = adjustl(ch_posx)
      ch_posy = adjustl(ch_posy)
      ch_posz = adjustl(ch_posz)

      open(44+nr,file='mic_'//trim(ch_posx)//'_'//trim(ch_posy)//'_'//trim(ch_posz)//'_'//trim(ch_rank)//'.res',position='append')
         write(44+nr,*) 
         write(44+nr,*)
         write(44+nr,*) "# Density Fluctuation for ",trim(s_par%problem_name)
         if(record .eqv. .true. ) then
         write(44+nr,*) "# Sensing pos 1    was x: ",pos(1),pos(1)-xmin(1)+1 
         write(44+nr,*) "#                      y: ",pos(2),pos(2)-xmin(2)+1 
         write(44+nr,*) "#                      z: ",pos(3),pos(3)-xmin(3)+1 
         endif
      write(44+nr,*) "#                  Omega: ",s_par%omega
      write(44+nr,*) "#                   umax: ",s_par%umax 
      write(44+nr,*) "#                     dx: ",s_par%gx(1)
      write(44+nr,*) "#                   dist: ",s_par%obst_l
      write(44+nr,*) "#       Relaxation model: ",s_par%modelname
      write(44+nr,*) "# tstep     dens   " 
   endif
   if(tstep>0) then
      write(44+nr,'(i10,3e18.10)') tstep,rholoc
   endif
      if(s_par%goend .eqv. .true. .or. gtstep_cur == s_par%tMax) then
         close(44+nr)
      endif
   endif

end subroutine microphone





subroutine set_bnd(lb_dom,state,fIn,rho,u,s_par,prc,meas,tStep)

!------------------------------------------------------------------------
!
! set the boundary conditions in each timestep.
!


   implicit none
   type(lb_block)       :: lb_dom
   type(sim_parameter)  :: s_par 
   type(mpl_var)        :: prc
   type(measure)        :: meas
   integer              :: tStep
   integer              :: state(lb_dom%lx(1),lb_dom%lx(2),lb_dom%lx(3))
   integer              :: i,j,k,l,x_start,x_end,y_start,y_end,z_start,z_end
   integer              :: incr_x,incr_y,incr_z,lx,ly,lz
   integer              :: pos(3) 
   real(R8B),dimension(NDX(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
   real(R8B),dimension(NDX(1:NDIM,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: u
   real(R8B),dimension(0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1)    :: rho

   intent(inout) :: fIn,rho,u


   call cpu_time_measure(meas%tSt_comm) 


   lx = lb_dom%lx(1)
   ly = lb_dom%lx(2)
   lz = lb_dom%lx(3)


   ! Set the different starting coordinates:
   ! If crd == 1     ->  first cell is boundary!   -> start bc on 2 
   ! If crd == np-1  ->  last  cell is boundary!   -> start bc on lx-1 

   if(prc%crd(1) == 0)          then; x_start = 2;    else; x_start=1;  end if
   if(prc%crd(1) == prc%np(1)-1) then; x_end   = lx-1; else; x_end  =lx; end if
   if(prc%crd(2) == 0)          then; y_start = 2;    else; y_start=1;  end if
   if(prc%crd(2) == prc%np(2)-1) then; y_end   = ly-1; else; y_end  =ly; end if
#ifdef D2Q9
   z_start=1; z_end=1
#endif
#ifdef D3Q19
   if(prc%crd(3) == 0)          then; z_start = 2;    else; z_start=1;  end if
   if(prc%crd(3) == prc%np(3)-1) then; z_end   = lz-1; else; z_end  =lz; end if
#endif

   ! Set the delta increase of coordinates for process coordinates > 1 
   incr_x = prc%crd(1)*lx
   incr_y = prc%crd(2)*ly
   incr_z = prc%crd(3)*lz

!--------------------------------------------------------------------
! Treat non-reflecting and periodic boundaries 

!if (gtstep_cur==0) then
!   call nrbc_init(lb_dom,fIn,s_par,prc)
!  do ind=1,nnod
!     lb_dom%fplus(i,ind) = 0.5*(fIn(NDX(ind,xc(0),yc(0),zc(0))) + fIn(NDX(opp(ind),xc(0),yc(0),zc(0))))
!  enddo
!endif
 

!   call set_nrbc_antibb(lb_dom,pnt_nxt,pnt_cur,s_par,prc,meas)
!   call set_nrbc(lb_dom,fIn,s_par,prc,meas)
   call set_periodic(lb_dom,fIn,s_par,prc,meas)
 
! calculate initial fplus for nrbc
      ! calculate correction term fplus

 
! Treat inlet / outlet 
 

!--------------------------------------------------------------------
! Select defined problem  

      select case(s_par%problem)

      case(plate_trail)

      case(gaussian)

      case(gauss_convect)

      case(corotating_vortex)

#ifdef CRVP
        ! check if pressure averaging should be active
        if(gtstep_cur == 0) then 
        ! reset density sum field
           lb_dom%rho_inc = 0.d0 
        endif

        if(s_par%calc_rho_inc .eqv. .true. ) then
           ! Copy the density field to current density sum field
           if(gtstep_cur == 5*s_par%crvp_period + s_par%tInitial) then 
        ! this should be the last iteration of averaging
         write(*,*) gtstep_cur,'last averaging iteration!!'
         write(*,*) '         averaged over',gtstep_cur-s_par%crvp_av_start-1,'time steps'
            ! how long was being averaged? gtstep_cur-s_par%crvp_av_start 
             lb_dom%rho_inc = lb_dom%rho_inc/(gtstep_cur-s_par%crvp_av_start-1)
            else
              lb_dom%rho_inc = lb_dom%rho_inc + lb_dom%rho

           endif
        endif
        if(gtstep_cur > 3*s_par%crvp_period + s_par%tInitial .and. &
           gtstep_cur < 5*s_par%crvp_period + s_par%tInitial ) then 
           if(s_par%calc_rho_inc .eqv. .false. ) s_par%crvp_av_start = gtstep_cur 
            s_par%calc_rho_inc = .true. 
        else
            s_par%calc_rho_inc = .false.
        endif
#endif /* CRVP */

        ! read out pressure at stationary microphone
        pos = (/ int(5.d0/8.d0*real(s_par%gx(1))), s_par%gx(2)/2, s_par%gx(3)/2+1    /)
        call microphone(lb_dom,pos,fIn,s_par,prc,tstep,1)
        pos = (/ int(6.d0/8.d0*real(s_par%gx(1))), s_par%gx(2)/2, s_par%gx(3)/2+1    /)
        call microphone(lb_dom,pos,fIn,s_par,prc,tstep,2)
        pos = (/ int(7.d0/8.d0*real(s_par%gx(1))), s_par%gx(2)/2, s_par%gx(3)/2+1    /)
        call microphone(lb_dom,pos,fIn,s_par,prc,tstep,3)


!----------------------------------------
! Cylinder in channel test with eq boundaries 

      case(channel)
! unstable boundary condition!! better use zou
         do k=z_start,z_end
         if(prc%crd(1)== 0) then
         do j=y_start,y_end
            i=1
            u(NDX(2,i,j,k)) = 0._R8B 
#if defined D2Q9
            u(NDX(1,i,j,k)) = s_par%umax*get_parabolic2d(real(j+incr_y,R8B),real(s_par%gx(2)+1,R8B))
#elif defined D3Q19
!FIXME errorous u calculation
            u(NDX(1,i,j,k)) = s_par%umax*get_parabolic3d(real(j+incr_y,R8B) &
            ,real(s_par%gx(2)+1,R8B),real(k+incr_z,R8B),real(s_par%gx(3)+1,R8B))
            u(NDX(3,i,j,k)) = 0._R8B 
#endif
            rho(i,j,k) = 1._R8B 
            call calc_fEq(lb_dom,rho,u,fIn,i,i,j,j,k,k,s_par)
         end do
         end if
            if(prc%crd(1)== prc%np(1) -1) then
         do j=y_start,y_end
            i=lx
            u(NDX(2,i,j,k)) = 0._R8B 
#if defined D2Q9
            u(NDX(1,i,j,k)) = s_par%umax*get_parabolic2d(real(j+incr_y,R8B),real(s_par%gx(2)+1,R8B))
#elif defined D3Q19
            u(NDX(1,i,j,k)) = s_par%umax*get_parabolic3d(real(j+incr_y,R8B)        &
         ,real(s_par%gx(2)+1,R8B),real(k+incr_z,R8B),real(s_par%gx(3)+1,R8B))
            u(NDX(3,i,j,k)) = 0._R8B 
#endif
            rho(i,j,k) = 1._R8B 
            call calc_fEq(lb_dom,rho,u,fIn,i,i,j,j,k,k,s_par)
         end do
         end if
         enddo          

   

!----------------------------------------
! Cylinder in channel test with eq boundaries 

      case(cylinder)
! unstable boundary condition!! better use zou
         do k=z_start,z_end
         if(prc%crd(1)== 0) then
         do j=y_start,y_end
            i=1
            u(NDX(2,i,j,k)) = 0._R8B 
#if defined D2Q9
            u(NDX(1,i,j,k)) = s_par%umax*get_parabolic2d(real(j+incr_y,R8B),real(s_par%gx(2)+1,R8B))
#elif defined D3Q19
!FIXME errorous u calculation
            u(NDX(1,i,j,k)) = s_par%umax*get_parabolic3d(real(j+incr_y,R8B) &
            ,real(s_par%gx(2)+1,R8B),real(k+incr_z,R8B),real(s_par%gx(3)+1,R8B))
            u(NDX(3,i,j,k)) = 0._R8B 
#endif
            rho(i,j,k) = 1._R8B 
            call calc_fEq(lb_dom,rho,u,fIn,i,i,j,j,k,k,s_par)
         end do
         end if
            if(prc%crd(1)== prc%np(1) -1) then
         do j=y_start,y_end
            i=lx
            u(NDX(2,i,j,k)) = 0._R8B 
#if defined D2Q9
            u(NDX(1,i,j,k)) = s_par%umax*get_parabolic2d(real(j+incr_y,R8B),real(s_par%gx(2)+1,R8B))
#elif defined D3Q19
            u(NDX(1,i,j,k)) = s_par%umax*get_parabolic3d(real(j+incr_y,R8B)        &
         ,real(s_par%gx(2)+1,R8B),real(k+incr_z,R8B),real(s_par%gx(3)+1,R8B))
            u(NDX(3,i,j,k)) = 0._R8B 
#endif
            rho(i,j,k) = 1._R8B 
            call calc_fEq(lb_dom,rho,u,fIn,i,i,j,j,k,k,s_par)
         end do
         end if
         enddo          



!----------------------------------------
! Cylinder in channel

      case(100)

   !---------------------------
   ! Standard Zou / He boundary conditions (only working for 2d so far) 

#ifdef D2Q9

        ! ---- macroscopic boundary conditions
         k=1
!         Length             = real(s_par%gx(2))
!         y_phys(:)     = y(:)
         do j=y_start,y_end
            if(prc%crd(1) == 0) then 
               ! INLET: reset to poisseuille Profile
               do i=1,1
               u(NDX(1,i,j,k))   =  s_par%umax*get_parabolic2d(real(j+incr_y,R8B),real(s_par%gx(2)+1,R8B))!4._R8B*s_par%umax/real(s_par%gx(2)+1)**2*&
               !      & ((real(j+incr_y))*real(s_par%gx(2)+1)-real(j+incr_y)*real(j+incr_y))
               u(NDX(2,i,j,k))   = 0._R8B
               rho(i,j,k)   = 1._R8B / (1-u(NDX(1,i,j,k)))*(fIn(NDX(1,i,j,k)) &
            + fIn(NDX(3,i,j,k)) + fIn(NDX(5,i,j,k)) &
            + 2*(fIn(NDX(4,i,j,k)) + fIn(NDX(7,i,j,k)) + fIn(NDX(8,i,j,k))))
 !              u(NDX(:,i-1,j,k)) = u(NDX(:,i,j,k))
 !              rho(i-1,j,k) = rho(i,j,k)
!               write(*,*) i,j,k,"dichte",rho(i,j,k),"geschwindigkeit",u(NDX(:,i,j,k)) 
               enddo
            end if
            ! OUTLET: reset to constant pressure
            if(prc%crd(1)== prc%np(1) -1) then
         do i=lx,lx
               rho(i,j,k) = 1._R8B
               u(NDX(1,i,j,k)) = -1 + 1._R8B/rho(i,j,k)*(fIn(NDX(1,i,j,k))&
                  +fIn(NDX(3,i,j,k))+fIn(NDX(5,i,j,k)) &
                  + 2*(fIn(NDX(2,i,j,k))+fIn(NDX(6,i,j,k))+fIn(NDX(9,i,j,k))))
               u(NDX(2,i,j,k))   = 0._R8B
!            u(NDX(:,i+1,j,k)) = u(NDX(:,i,j,k))
!            rho(i+1,j,k) = rho(i,j,k)
!            write(*,*) i,j,k,"dichte",rho(i,j,k),"geschwindigkeit",u(NDX(:,i,j,k)) 
         end do
            end if
         end do
         ! ---- microscopic boundary condition
         do j=y_start,y_end
            ! INLET: Zou/He BC
            if(prc%crd(1) == 0) then 
         do i=1,1
            fIn(NDX(2,i,j,k)) = fIn(NDX(4,i,j,k)) + 2._R8B/3._R8B*rho(i,j,k)*u(NDX(1,i,j,k))
            fIn(NDX(6,i,j,k)) = fIn(NDX(8,i,j,k)) + 1._R8B/2._R8B*(fIn(NDX(5,i,j,k))-fIn(NDX(3,i,j,k))) &
                 &                  + 1._R8B/2._R8B*rho(i,j,k)*u(NDX(2,i,j,k) )&
                 &                  + 1._R8B/6._R8B*rho(i,j,k)*u(NDX(1,i,j,k)) 
            fIn(NDX(9,i,j,k)) = fIn(NDX(7,i,j,k)) + 1._R8B/2._R8B*(fIn(NDX(3,i,j,k))-fIn(NDX(5,i,j,k))) &
                &                  - 1._R8B/2._R8B*rho(i,j,k)*u(NDX(2,i,j,k)) &
                &                  + 1._R8B/6._R8B*rho(i,j,k)*u(NDX(1,i,j,k)) 
!            fIn(NDX(:,i-1,:,:)) = fIn(NDX(:,i,:,:))
         end do
            end if
           ! OUTLET:Zou/He BC
            if(prc%crd(1)== prc%np(1) -1) then
         do i=lx,lx
           fIn(NDX(4,i,j,k)) = fIn(NDX(2,i,j,k) )- 2._R8B/3._R8B*rho(i,j,k)*u(NDX(1,i,j,k))
           fIn(NDX(8,i,j,k)) = fIn(NDX(6,i,j,k)) + 1._R8B/2._R8B*(fIn(NDX(3,i,j,k))-fIn(NDX(5,i,j,k))) &
                &                  - 1._R8B/2._R8B*rho(i,j,k)*u(NDX(2,i,j,k)) &
                &                  - 1._R8B/6._R8B*rho(i,j,k)*u(NDX(1,i,j,k) )
           fIn(NDX(7,i,j,k)) = fIn(NDX(9,i,j,k)) + 1._R8B/2._R8B*(fIn(NDX(5,i,j,k))-fIn(NDX(3,i,j,k))) &
                &                  + 1._R8B/2._R8B*rho(i,j,k)*u(NDX(2,i,j,k)) &
                &                  - 1._R8B/6._R8B*rho(i,j,k)*u(NDX(1,i,j,k)) 
        end do
            end if
        end do
#endif

!----------------------------------------
! Flute 

      case(flute)

        pos = (/ lb_dom%lx(1)/2+1, lb_dom%lx(2)-10, lb_dom%lx(3)/2+1    /)
        call microphone(lb_dom,pos,fIn,s_par,prc,tstep,1)

#ifdef OLD
#ifdef D2Q9
        ! ---- macroscopic boundary conditions
         k=1
!         Length             = real(s_par%gx(2))
!         y_phys(:)     = y(:)
         do j=y_start,y_end
            do i=1,5
            if(lb_dom%state(i,j,1)==inlet) then 
               ! INLET: reset to poisseuille Profile
!               u(NDX(1,i,j,k))   = 4._R8B*s_par%umax/real(lb_dom%num_in(2)+1)**2*&
!                     & ((real(j))*real(lb_dom%num_in(2)+1)-real(j)*real(j))
!FIXME this is my own poisseuille profile
               u(NDX(1,i,j,k)) = s_par%umax*get_parabolic2d(real(j+1,R8B),real(2*s_par%obst_r+1,R8B))  
               u(NDX(2,i,j,k))   = 0._R8B
               rho(i,j,k)   = 1._R8B / (1-u(NDX(1,i,j,k)))*(fIn(NDX(1,i,j,k)) &
               + fIn(NDX(3,i,j,k)) + fIn(NDX(5,i,j,k)) &
               + 2*(fIn(NDX(4,i,j,k)) + fIn(NDX(7,i,j,k)) + fIn(NDX(8,i,j,k))))
            end if
            end do
            ! OUTLET: reset to constant pressure
            i=lx
            if(lb_dom%state(i,j,1) == outlet) then
               rho(i,j,k) = 1._R8B
               u(NDX(1,i,j,k)) = 0.d0 !-1 + 1._R8B/rho(i,j,k)*(fIn(NDX(1,i,j,k))+fIn(NDX(3,i,j,k))+fIn(NDX(5,i,j,k)) &
                                           !  & + 2*(fIn(NDX(2,i,j,k))+fIn(NDX(6,i,j,k))+fIn(NDX(9,i,j,k))))
               u(NDX(2,i,j,k))   = 0._R8B
                       
            end if
         end do

         

         ! ---- microscopic boundary condition
         do j=y_start,y_end
            ! INLET: Zou/He BC
            i = 1
            if(lb_dom%state(i,j,k) == inlet) then 
            fIn(NDX(2,i,j,k)) = fIn(NDX(4,i,j,k)) + 2._R8B/3._R8B*rho(i,j,k)*u(NDX(1,i,j,k))
            fIn(NDX(6,i,j,k)) = fIn(NDX(8,i,j,k)) + 1._R8B/2._R8B*(fIn(NDX(5,i,j,k))-fIn(NDX(3,i,j,k))) &
                 &                  + 1._R8B/2._R8B*rho(i,j,k)*u(NDX(2,i,j,k) )&
                 &                  + 1._R8B/6._R8B*rho(i,j,k)*u(NDX(1,i,j,k)) 
            fIn(NDX(9,i,j,k)) = fIn(NDX(7,i,j,k)) + 1._R8B/2._R8B*(fIn(NDX(3,i,j,k))-fIn(NDX(5,i,j,k))) &
                &                  - 1._R8B/2._R8B*rho(i,j,k)*u(NDX(2,i,j,k)) &
                &                  + 1._R8B/6._R8B*rho(i,j,k)*u(NDX(1,i,j,k)) 
            end if
           ! OUTLET:Zou/He BC
            if(prc%crd(1)== prc%np(1) -1) then
           i = lx
           fIn(NDX(4,i,j,k)) = fIn(NDX(2,i,j,k) )- 2._R8B/3._R8B*rho(i,j,k)*u(NDX(1,i,j,k))
           fIn(NDX(8,i,j,k)) = fIn(NDX(6,i,j,k)) + 1._R8B/2._R8B*(fIn(NDX(3,i,j,k))-fIn(NDX(5,i,j,k))) &
                &                  - 1._R8B/2._R8B*rho(i,j,k)*u(NDX(2,i,j,k)) &
                &                  - 1._R8B/6._R8B*rho(i,j,k)*u(NDX(1,i,j,k) )
           fIn(NDX(7,i,j,k)) = fIn(NDX(9,i,j,k)) + 1._R8B/2._R8B*(fIn(NDX(5,i,j,k))-fIn(NDX(3,i,j,k))) &
                &                  + 1._R8B/2._R8B*rho(i,j,k)*u(NDX(2,i,j,k)) &
                &                  - 1._R8B/6._R8B*rho(i,j,k)*u(NDX(1,i,j,k)) 
            end if
        end do
      rho(:,lb_dom%lx(2),:) = 1.0d0
      u(NDX(:,:,lb_dom%lx(2),:)) = 0.d0
      call calc_feq(lb_dom,rho,u,lb_dom%fIn,0,0,lb_dom%lx(2),lb_dom%lx(2),0,0)
#endif
#endif


!----------------------------------------
! Lid-driven cavity

     case(cavity)
        ! Reset boundary 
        ! macroscopic Variables
#ifdef D3Q19
        j = ly
        if(prc%crd(2) == prc%np(2)-1) then !wenn Block in oberster Reihe
         do k=z_start,z_end
            do i=x_start,x_end
              ! top lid
               rho(i,j,k)   = 1._R8B
               u(NDX(1,i,j,k))   = s_par%umax/sqrt(2.d0)
               u(NDX(2,i,j,k))   = 0._R8B
               u(NDX(3,i,j,k))   = s_par%umax/sqrt(2.d0)
            end do
         end do
         end if

         if(prc%crd(2) == prc%np(2)-1) then
            call calc_fEq(lb_dom,rho,u,fIn,x_start,x_end,ly,ly,z_start,z_end,s_par) !top lid
         end if
#endif
#ifdef D2Q9
         k = 1
         j = ly
         if(prc%crd(2) == prc%np(2)-1) then
            do i=x_start,x_end
               ! top lid
               rho(i,j,k)   = 1._R8B
               u(NDX(1,i,j,k))   = s_par%umax
               u(NDX(2,i,j,k))   = 0._R8B
            end do
        call calc_fEq(lb_dom,rho,u,fIn,x_start,x_end,ly,ly,0,0,s_par) ! top wall
         end if

#endif





!----------------------------------------
! Lid-driven cavity on all processes the same calculation
! this was used for testing different execution times on Nehalem 
! for the same amount of work

     case(cavity_same)
        ! Reset boundary 
        ! macroscopic Variables
#ifdef D3Q19
        j = ly
         do k=2,lz-1
            do i=2,lx-1       
              ! top lid
               rho(i,j,k)   = 1._R8B
               u(NDX(1,i,j,k))   = sqrt(s_par%umax)
               u(NDX(2,i,j,k))   = 0._R8B
               u(NDX(3,i,j,k))   = sqrt(s_par%umax)
            end do
         end do

        call calc_fEq(lb_dom,rho,u,fIn,2      ,lx-1 ,ly,ly,2      ,lz-1,s_par ) !top lid
#endif
#ifdef D2Q9
         k = 1
         j = ly
         do i=2, lx-1      
               ! top lid
               rho(i,j,k)   = 1._R8B
               u(NDX(1,i,j,k))   = s_par%umax
               u(NDX(2,i,j,k))   = 0._R8B
         end do

        call calc_fEq(lb_dom,rho,u,fIn,2      ,lx-1 ,ly,ly,0,0,s_par) ! top wall
#endif



!----------------------------------------
! Default case 

     
   case (planar_standing_wave)
!        pos = (/ int(1.d0/4.d0*real(s_par%gx(1))), s_par%gx(2)/2, s_par%gx(3)/2+1    /)
!        pos = (/ int(0.25*real(s_par%gx(1))+1.), s_par%gx(2)/2, s_par%gx(3)/2+1    /)
        pos = (/ 1+(s_par%gx(1)-2)/4, s_par%gx(2)/2, s_par%gx(3)/2+1    /)
        call microphone(lb_dom,pos,fIn,s_par,prc,tstep,1)

   end select


   call cpu_time_measure(meas%tEnd_comm)

   meas%sbnd_duration = meas%tEnd_comm - meas%tSt_comm + meas%sbnd_duration



end  subroutine set_bnd
!------------------------------------------------------------------------








end module lb_bc        


