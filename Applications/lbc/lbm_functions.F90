module lbm_functions
   use mpl_lib
   use lbmodel
   use function_set
   use timing
#include "include/replace.h"
contains







subroutine init_field(fIn,u,rho,state,s_par)
!------------------------------------------------------------------------
!
! initialization routine of initial values, done by master process only.
! values have to be distributed to slave processes
!
   implicit none
   type(sim_parameter)              :: s_par
   real(R8B),dimension(s_par%gx(2)) :: y_phys
   real(R8B) :: fIn(LB_NODE(nnod,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
   real(R8B) :: rho(s_par%gx(1),s_par%gx(2),s_par%gx(3))
   real(R8B) :: u(LB_NODE(NDIM,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
   real(R8B) :: radius 
   integer   :: state(s_par%gx(1),s_par%gx(2),s_par%gx(3))
   integer                 :: ii,jj,kk,i,j,k,l,N
   integer                 :: xmax,ymin,ymax 
   intent(inout) :: u,rho,fIn,state

   ! assign physical coordinates
   do j=1,s_par%gx(2) 
     y_phys(j) = s_par%g_y(j)
   end do
   ! reset complete domain to infinite solution
      rho = 1.0_R8B
      u = 0.0_R8B



!--------------------------------------------------------------------
! now assign initial conditions depending on the problem 


   select case(s_par%problem)

!----------------------------------------
! Single Taylor Vortex decay


   case(taylor_vortex)
      ! vortex core radius
      
      call init_taylor_vortex(u,s_par)


!----------------------------------------
! Single Taylor Vortex decay


   case(planar_standing_wave)
      ! reset complete domain to eq
      do kk =1,s_par%gx(3) 
      do jj =1,s_par%gx(2) 
      do ii =1,s_par%gx(1) 
         rho(ii,jj,kk) = 1.0d0 + s_par%umax/100.d0*sin(2.d0*pi*real(ii-1)/real(s_par%gx(1)-2))*cs2inv
         u(LB_NODE(1,ii,jj,kk)) = (rho(ii,jj,kk)-1.d0)*cs 
         u(LB_NODE(2:NDIM,ii,jj,kk)) = 0.d0 
   ! it's a standing wave, if you set u(LB_NODE(1))=0.
      enddo
      enddo
      enddo

!----------------------------------------
! Flat plate with trailing edge

   case(plate_trail)
      u(LB_NODE(1,1,:,:)) = s_par%umax
!      do ii=1,s_par%gx(1)
!!      rho(ii,:,:) = rho(ii,:,:)*(1.15 - 0.15*ii/s_par%gax
!      enddo
!----------------------------------------
! Cylinder Impulse   

   case(flute)
      ! cylinder with travelling impulse from close end

      xmax = s_par%obst_x 
      ymin = s_par%obst_y
      ymax = s_par%obst_y+2*s_par%obst_r

       do j=1,2*s_par%obst_r+1
          ! inlet links 
          u(LB_NODE(1,1:xmax ,j+ymin-1,:)) = s_par%umax*get_parabolic2d(real(j,R8B),real(2*s_par%obst_r+2,R8B)) 
          u(LB_NODE(2,1:xmax ,j+ymin-1,:)) = 0.0 
          rho(1:xmax ,j+ymin-1,:) = rho(1:xmax ,j+ymin-1,:)*1.15
       end do



!----------------------------------------
! Cylinder Impulse   

   case(cylinder_impulse)
      ! cylinder with travelling impulse from close end
      N=10   ! Impulsweite
         do i=1,N
            rho(1+i,s_par%gx(2)/2-5:s_par%gx(2)/2+5,:) = & 
               & 1.d0+0.001*(0.5+0.5*cos(2*pi*i/N+pi))
            u(LB_NODE(1,1+i,s_par%gx(2)/2-5:s_par%gx(2)/2+5,:)) = 0.001/cs*(0.5+0.5*cos(2*pi*i/N+pi))
         enddo 
 




!----------------------------------------
! Cylinder in channel

    case(cylinder)
       ! Poisseuille Profile for Flow around Cylinder
      do k=1,s_par%gx(3)
        do j=1,s_par%gx(2)
          do i=1,s_par%gx(1)
#ifdef D2Q9
               u(LB_NODE(1,i,j,k)) = s_par%umax*get_parabolic2d(real(j,R8B),real(s_par%gx(2)+1,R8B))
               if(state(i,j,k)==wall) u(LB_NODE(1,i,j,k))=0._R8B
               u(LB_NODE(2,i,j,k))   = 0._R8B
#endif
#ifdef D3Q19
               u(LB_NODE(1,i,j,k)) = s_par%umax*get_parabolic3d(real(j,R8B),real(s_par%gx(2)+1,R8B), &
                                       real(k,R8B),real(s_par%gx(3)+1,R8B))
               u(LB_NODE(2,i,j,k))   = 0._R8B
               u(LB_NODE(3,i,j,k))   = 0._R8B
#endif
               rho(i,j,k)   = 1._R8B
             end do
          end do
       end do



!----------------------------------------
! Lid-driven Cavity  


    case(cavity)
       ! Driven Cavity
       do k=1,s_par%gx(3)
          do j=1,s_par%gx(2)
             do i=1,s_par%gx(1)
               ! In Benchmark of of Hou, Chen and Doolen (1995) rho = 2.7
               rho(i,j,k)   = 1.0_R8B !2.7
               u(LB_NODE(1,i,j,k))   = 0._R8B
               u(LB_NODE(2,i,j,k))   = 0._R8B
#ifdef D3Q19
               u(LB_NODE(3,i,j,k))   = 0._R8B
#endif
             end do
          end do
       end do



!----------------------------------------
! Density shock 

 
    case(shock)
       ! Shock propagation
       do k=1,s_par%gx(3)
          do j=1,s_par%gx(2)
             do i=1,s_par%gx(1)
               ! In Benchmark of of Hou, Chen and Doolen (1995) rho = 2.7
               rho(i,j,k)   = 1.0_R8B 
               u(LB_NODE(1,i,j,k))   = 0._R8B
               u(LB_NODE(2,i,j,k))   = 0._R8B
#ifdef D3Q19
               u(LB_NODE(3,i,j,k))   = 0._R8B
#endif
             end do
          end do
          end do
#ifdef D3Q19
      do k=s_par%gx(3)/2-2,s_par%gx(3)/2+2
#else
      k=1
#endif
          do j=s_par%gx(2)/2-1 ,s_par%gx(2)/2+1
             do i=s_par%gx(1)/2-1,s_par%gx(1)/2+1
               rho(i,j,k)   = s_par%umax*rho(i,j,k)
             end do
          end do
#ifdef D3Q19
       end do
#endif



!----------------------------------------
! Gaussian Pulse 
 
    case(gaussian)
       radius=real(s_par%obst_r)
       do k=1,s_par%gx(3)
          do j=1,s_par%gx(2)
             do i=1,s_par%gx(1)
        rho(i,j,k)   = s_par%umax*exp(-log(2.)*real(i-s_par%gx(1)/2)*real(i-s_par%gx(1)/2)/radius/radius)  &
                                 *exp(-log(2.)*real(j-s_par%gx(2)/2)*real(j-s_par%gx(2)/2)/radius/radius)  &
#ifdef D3Q19
                                 *exp(-log(2.)*real(k-s_par%gx(3)/2)*real(k-s_par%gx(3)/2)/radius/radius)  &
#endif
                                 + rho(i,j,k)
             enddo
         enddo
      enddo




!----------------------------------------
! Convected Gaussian Pulse 
 
    case(gauss_convect)
       radius=10. !real(s_par%obst_r)
       do k=1,s_par%gx(3)
          do j=1,s_par%gx(2)
             do i=1,s_par%gx(1)
          rho(i,j,k)   = s_par%umax/cs2inv*exp(-log(2.)*real(i-s_par%gx(1)/2)*real(i-s_par%gx(1)/2)/radius/radius)  &
#ifdef PULSE_3d
                                   *exp(-log(2.)*real(j-s_par%gx(2)/2)*real(j-s_par%gx(2)/2)/radius/radius)  &
#ifdef D3Q19
                                   *exp(-log(2.)*real(k-s_par%gx(3)/2)*real(k-s_par%gx(3)/2)/radius/radius)  &
#endif
#endif              
                     + rho(i,j,k)
               u(LB_NODE(1,i,j,k)) = (rho(i,j,k) - 1.0)/cs         !s_par%umax
               u(LB_NODE(2:NDIM,i,j,k)) =  0.0 !s_par%umax !0.0
             enddo
         enddo
      enddo



!----------------------------------------
! Corotating Vortex 


   case(corotating_vortex)
      ! reset complete domain to eq
      rho      = 1.0d0
      u        = 0.0d0
      
      call init_vortex_macr(u,s_par)

!----------------------------------------
! Default case 


   case default
      rho = 1.0d0
      u   = 0.0d0
      write(*,*) " Default case: Set complete domain to rho=1, u=0" 
   end select




  end subroutine init_field
  !------------------------------------------------------------------------


  subroutine init_taylor_vortex(u,s_par)
   !-----------------------------------------------------------------------
   ! Initialization routine for Taylor Vortex
   ! x0=50
   ! y0=50
   ! adius=10.
   ! a=radius*sqrt(2.)
   ! phi0=0.0005
   ! set xrange[-100:100]
   ! set yrange[-100:100]
   ! splot 2*(x-x0)/a/a*phi0*exp(-((x-x0)*(x-x0)+(y-y0)*(y-y0))/a/a), \
   ! -2*(y-y0)/a/a*phi0*exp(-((x-x0)*(x-x0)+(y-y0)*(y-y0))/a/a)
   ! pause -1

      implicit none
      type(sim_parameter)           :: s_par
      real(R8B) :: u(LB_NODE(NDIM,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
      integer   :: i,j,k,l
      integer   :: x0,y0,z0 

      ! corotating vortex
      real(R8B) :: vor_rc, vor_pos1(3),vor_pos2(3),vor_r
      real(R8B) :: delta_x,delta_y,delta_z
      real(R8B) :: cutoff_radius 
      real(R8B) :: phi0,vor_phi,pos_r

      double precision MachRot, MachCore
      double precision p0, rho0, c0, gamma
      double precision u_0, u_c, omega
      
      double precision Vortex_Gamma, r_0, r_c, alpha_Vortex
      double precision Vortex_Strength, rholoc, norm_factor
         
 
      phi0 = 20*s_par%umax

      x0    = real(s_par%gx(1))/2.0d0+0.5
      y0    = real(s_par%gx(2))/2.0d0+0.5
      z0    = real(s_par%gx(3))/2.0d0+0.5
      vor_r = real(s_par%gx(1))/20.0d0*sqrt(2.0)
      cutoff_radius = real(s_par%gx(1))/2.0d0-1.0d0

      do k=1,s_par%gx(3)
         do j=1,s_par%gx(2)
            do i=1,s_par%gx(1)


               vor_phi = phi0*exp(-(dble(i-x0)*dble(i-x0) + dble(j-y0)*dble(j-y0))/(vor_r*vor_r))
               pos_r    = sqrt((dble(i)-x0 )**2 + (dble(j)-y0)**2 +(dble(k)-z0)**2)
               if(pos_r  < cutoff_radius) then
                  u(LB_NODE(1,i,j,k)) = (-2.d0*dble(j-y0)/vor_r/vor_r*vor_phi)
                  u(LB_NODE(2,i,j,k)) = ( 2.d0*dble(i-x0)/vor_r/vor_r*vor_phi)
               endif
            enddo
         enddo
      enddo


  end subroutine init_taylor_vortex
  !------------------------------------------------------------------------



  subroutine init_vortex_macr(u,s_par)
   !-----------------------------------------------------------------------
   ! Initialization routine for Corotating Vortex Pair. See Dissertation Roller 2003 or Schwartzkopf 2006, 
   ! M S Howe Theory of Vortex Sound, P. 121 
   ! Gnuplot script is located in scripts/gnuplot/ini_crvp.gpl 

      implicit none
      type(sim_parameter)           :: s_par
      real(R8B) :: u(LB_NODE(NDIM,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
      integer   :: i,j,k,kk,l

      real(R8B) :: x0,y0,z0,x1 

      ! corotating vortex
      real(R8B) :: vor_rc, vor_pos1(3),vor_pos2(3),vor_r
      real(R8B) :: delta_x,delta_y,delta_z
      real(R8B) :: cutoff_radius,distance,factor 
      real(R8B) :: phi0,vor_phi,pos_r,prer

      double precision MachRot, MachCore
      double precision p0, rho0, c0, gamma
      double precision u_0, u_c, omega, u_max,u_cur

      double precision Vortex_Gamma, r_0, r_c, alpha_Vortex
      double precision Vortex_Strength, radius, rholoc, norm_factor
         
      u_max     = 0.0        ! 

#ifndef TAYLOR_CRVP 
                                ! Fluid dimensionslos
                                ! -------------------
      distance    = dble(s_par%obst_r)

      gamma       = 1.4d0
      p0          = 1.d0
      rho0        = 1.d0
      c0          = sqrt(gamma)
  
                                ! Vortex dimensionslos
                                ! --------------------
      MachRot     = sqrt(gamma)*0.08d0 ! RotationsMach-Zahl
      r_0         = 1.d0        ! r_0 = halber Wirbelabstand
      r_c         = 0.3d0*r_0       ! Core-Radius, Wirbelmodell 
      alpha_vortex= 1.256431d0  ! Parameter fuer Homentropic Vortex
                                ! dann max. Geschw.Betr. in r_c
      u_0         = s_par%umax       ! Ind. Geschw.Betrag in Entfernung2*r_0
      omega       = 1.d0        ! Winkelgeschw.
                                ! Homentropic Vortex
                                ! Zirkulation in Entfernung r_c
      Vortex_Gamma= 4.d0*Pi/(1.d0-exp(-4.d0*alpha_Vortex/r_c**2)) 
      Vortex_Strength=Vortex_Gamma/2.d0/Pi
      u_c         = Vortex_Gamma/2.d0/Pi/r_c *(1.d0-exp(-alpha_Vortex))
      MachCore    = u_c/c0*MachRot ! Mach-Zahl am Radius r_c

      x0    = dble(s_par%gx(1))/2.0d0 + distance*r_0 + 0.5
      x1    = dble(s_par%gx(1))/2.0d0 - distance*r_0 + 0.5
      y0    = dble(s_par%gx(2))/2.0d0 + 0.5
      z0    = dble(s_par%gx(3))/2.0d0 + 0.5

      cutoff_radius  = abs(dble(s_par%gx(1))/2.0d0 - distance*r_0)*0.95
      
                                !
                                !---Dokumentation
                                !
      print*,'  Wirbelpaar CRVP:'
      print*,'  -----------------------------------------------------'
      print*,'  Wirbelabstand/2        r_0    [-] :', r_0*distance
      print*,'  Core Radius            r_c    [-] :', r_c*distance
      print*,'  Rot.-Mach-Zahl in 2r_0 MaRot  [-] :', MachRot
      print*,'  Max.-Mach-Zahl in  r_c MaCore [-] :', MachCore
      print*,'  Zirkulation            Gamma  [-] :', Vortex_Gamma
      print*,'  GeschwindigkeitsBetrag in 2r_0[-] :', u_0
      print*,'  GeschwindigkeitsBetrag in  r_c[-] :', u_c
      print*,'  Zeit fuer einen Wirbel-Umlauf [-] :', 2.d0*Pi
      print*,'  Akustische Wellenlaenge       [-] :', Pi*c0/MachRot
      print*,'  Cutoff-Radius                 [-] :', cutoff_radius
      print*,'  -----------------------------------------------------'
      print*,'  Vortex_Strength               [-] :', Vortex_Strength
      print*,'  -----------------------------------------------------'
      
      do k=1,s_par%gx(3)
         do j=1,s_par%gx(2)
            do i=1,s_par%gx(1)
            radius = sqrt((dble(i)-x0)**2 + (dble(j)-y0)**2)
                                ! Homentropic Vortex 
            if(radius .ne. 0.0d0 .and. radius .le. cutoff_radius) &
            u(LB_NODE(1,i,j,k))=u(LB_NODE(1,i,j,k)) + Vortex_Gamma/2.d0/Pi/radius**2 &
                *(1.d0-exp(-alpha_Vortex*(radius/(r_c*distance))**2)) &
                *(-1.0d0)*(-dble(j)+y0)

            radius = sqrt((dble(i)-x1)**2 + (dble(j)-y0)**2)
            if(radius .ne. 0.0d0 .and. radius .le. cutoff_radius) &
            u(LB_NODE(1,i,j,k))=u(LB_NODE(1,i,j,k)) + Vortex_Gamma/2.d0/Pi/radius**2 &
                *(1.d0-exp(-alpha_Vortex*(radius/(r_c*distance))**2)) &
                *(-1.0d0)*(-dble(j)+y0)

            radius = sqrt((dble(i)-x0)**2 + (dble(j)-y0)**2)
                                ! Homentropic Vortex 
            if(radius .ne. 0.0d0 .and. radius .le. cutoff_radius) &
            u(LB_NODE(2,i,j,k))=u(LB_NODE(2,i,j,k)) + Vortex_Gamma/2.d0/Pi/radius**2 &
                *(1.d0-exp(-alpha_Vortex*(radius/(r_c*distance))**2)) &
                *(-dble(i)+x0)

            radius = sqrt((dble(i)-x1)**2 + (dble(j)-y0)**2)
            if(radius .ne. 0.0d0 .and. radius .le. cutoff_radius) &
            u(LB_NODE(2,i,j,k))=u(LB_NODE(2,i,j,k)) + Vortex_Gamma/2.d0/Pi/radius**2 &
                *(1.d0-exp(-alpha_Vortex*(radius/(r_c*distance))**2)) &
                *(-dble(i)+x1)

            if(u(LB_NODE(1,i,j,k)) > u_max) u_max=u(LB_NODE(1,i,j,k))
            if(u(LB_NODE(2,i,j,k)) > u_max) u_max=u(LB_NODE(2,i,j,k))
            enddo
         enddo
      enddo

      print*,'  Maximal Velocity before rescale   :', u_max 

      factor = s_par%umax/u_max 
      u = u*factor

      u_max=0.0d0      

      ! Re-check maximum velocity
      do k=1,s_par%gx(3)
         do j=1,s_par%gx(2)
            do i=1,s_par%gx(1)
            if(u(LB_NODE(1,i,j,k)) > u_max) u_max=u(LB_NODE(1,i,j,k))
            if(u(LB_NODE(2,i,j,k)) > u_max) u_max=u(LB_NODE(2,i,j,k))
            enddo
         enddo
      enddo
      print*,'  Maximal Velocity after rescale    :', u_max 

! Done initializing velocities. Pressure has to be found via the initial relaxation with fixed velocities




#else /* TAYLOR_CRVP */
 
      phi0 = s_par%umax
!      distance = 3.3/2.    ! vortices are diverging
      distance = 3.3        ! standard distance. Vortices do NOT move at all
      distance = s_par%obst_l !

      vor_r = real(s_par%obst_r)*sqrt(2.0)
      cutoff_radius = real(s_par%gx(1))/2.0d0-1.0d0

      write(*,*) 
      write(*,*) "Initializing vortices "
      write(*,*) "Cutoff Radius:   ",int(cutoff_radius)
      write(*,*) "Distance:        ",int(distance*vor_r)
      write(*,*) "Characteristic l:",s_par%obst_r
      write(*,*) "Factor:          ",s_par%obst_l
! Vortex 1

      do kk=1,2
      prer  = (-1.)**kk
      x0    = real(s_par%gx(1))/2.0d0 + prer*distance*vor_r + 0.5
      y0    = real(s_par%gx(2))/2.0d0 + 0.5
      z0    = real(s_par%gx(3))/2.0d0 + 0.5

      write(*,*) "Vortex ",kk,": ",int(x0),int(y0),int(z0)

      do k=1,s_par%gx(3)
         do j=1,s_par%gx(2)
            do i=1,s_par%gx(1)

               vor_phi = phi0*exp(-((dble(i)-x0)*(dble(i)-x0) + (dble(j)-y0)*(dble(j)-y0))/(vor_r*vor_r))
               pos_r    = sqrt((dble(i)-x0 )**2 + (dble(j)-y0)**2 +(dble(k)-z0)**2)
               if(pos_r  < cutoff_radius) then
                  u(LB_NODE(1,i,j,k)) = u(LB_NODE(1,i,j,k)) + (-2.d0*(dble(j)-y0)/vor_r/vor_r*vor_phi)
                  u(LB_NODE(2,i,j,k)) = u(LB_NODE(2,i,j,k)) + ( 2.d0*(dble(i)-x0)/vor_r/vor_r*vor_phi)
                  u_cur = sqrt(u(LB_NODE(1,i,j,k)**2 + u(LB_NODE(2,i,j,k))))
                  if(u_cur > u_max) u_max = u_cur
               endif
            enddo
         enddo
      enddo
   enddo

   write(*,*) "Maximum Velocity: ", u_max 
   write(*,*) 


#endif
  end subroutine init_vortex_macr  
  !------------------------------------------------------------------------




   subroutine init_field_each(lb_dom,s_par,prc)
   !------------------------------------------------------------------------
   !
   ! initialization routine for initial conditions of domains
   ! this routine is intended for usage of more than one processes
   ! centralized initialization can be done by using init_field routine

      implicit none
      type(lb_block),intent(inout)  :: lb_dom
      type(sim_parameter)           :: s_par
      type(mpl_var)                 :: prc
      integer                 :: i,j,k,l,lx(3)

      lx(:) = lb_dom%lx(:)
       

      ! set complete domain to density 1, u=0
       do k=1,lb_dom%lx(3)
          do j=1,lb_dom%lx(2)
             do i=1,lb_dom%lx(1)
               ! In Benchmark of of Hou, Chen and Doolen (1995) rho = 2.7
               lb_dom%rho(i,j,k)   = 1.0_R8B 

               lb_dom%u(LB_NODE(1,i,j,k))   = 0._R8B
               lb_dom%u(LB_NODE(2,i,j,k))   = 0._R8B
#ifdef D3Q19
               lb_dom%u(LB_NODE(3,i,j,k))   = 0._R8B
#endif
             end do
          end do
       end do


!--------------------------------------------------------------------
! now assign initial conditions depending on the problem 

      if(s_par%problem == cylinder) then
       ! Poisseuille Profile for Flow around Cylinder
      do k=1,lx(3)
        do j=1,lx(2)
          do i=1,lx(1)
#ifdef D2Q9
               lb_dom%u(LB_NODE(1,i,j,k)) = s_par%umax*get_parabolic2d(real(   &
                        j+lb_dom%lx(1)*prc%crd(1)  &
                        ,R8B),real(s_par%gx(2)+1,R8B))
               if(lb_dom%state(i,j,k)==wall) lb_dom%u(LB_NODE(1,i,j,k))=0._R8B
               lb_dom%u(LB_NODE(2,i,j,k))   = 0._R8B
#endif
#ifdef D3Q19
               lb_dom%u(LB_NODE(1,i,j,k)) = s_par%umax*get_parabolic3d( &
                  real(j+lb_dom%lx(2)*prc%crd(2),R8B),real(s_par%gx(2)+1,R8B), &
                  real(k+lb_dom%lx(3)*prc%crd(3),R8B),real(s_par%gx(3)+1,R8B))
               lb_dom%u(LB_NODE(2,i,j,k))   = 0._R8B
               lb_dom%u(LB_NODE(3,i,j,k))   = 0._R8B
#endif
               lb_dom%rho(i,j,k)   = 1._R8B
             end do
          end do
       end do
#ifdef MURKS
       do k=1,s_par%gx(3)
        do j=1,s_par%gx(2)
          do i=1,s_par%gx(1)
#ifdef D2Q9
               lb_dom%u(LB_NODE(1,i,j,k))   = 4._R8B*s_par%umax/real(s_par%gx(2)+1)**2*(real(j*s_par%gx(2)+1)-real(j)*real(j))
               if(lb_dom%state(i,j,k)==wall) lb_dom%u(LB_NODE(1,i,j,k))=0._R8B
               lb_dom%u(LB_NODE(2,i,j,k))   = 0._R8B
#endif
#ifdef D3Q19
               lb_dom%u(LB_NODE(1,i,j,k))   = 4._R8B*s_par%umax/real(s_par%gx(2)+1)**2  &
         & *(real(j*s_par%gx(2)+1)-real(j)*real(j))*(real(k*s_par%gx(3)+1)-real(k)*real(k))
               lb_dom%u(LB_NODE(2,i,j,k))   = 0._R8B
               lb_dom%u(LB_NODE(3,i,j,k))   = 0._R8B
#endif
               lb_dom%rho(i,j,k)   = 1._R8B

             end do
          end do
       end do
#endif
    else if(s_par%problem==cavity) then
       ! Driven Cavity
       do k=1,lb_dom%lx(3)
          do j=1,lb_dom%lx(2)
             do i=1,lb_dom%lx(1)
               ! In Benchmark of of Hou, Chen and Doolen (1995) rho = 2.7
               lb_dom%rho(i,j,k)   = 1.0_R8B 
               lb_dom%u(LB_NODE(1,i,j,k))   = 0._R8B
               lb_dom%u(LB_NODE(2,i,j,k))   = 0._R8B
#ifdef D3Q19
               lb_dom%u(LB_NODE(3,i,j,k))   = 0._R8B
#endif
             end do
          end do
       end do
      
    else if(s_par%problem==cavity_same) then
       ! Driven Cavity
       do k=1,lb_dom%lx(3)
          do j=1,lb_dom%lx(2)
             do i=1,lb_dom%lx(1)
               ! In Benchmark of of Hou, Chen and Doolen (1995) rho = 2.7
               lb_dom%rho(i,j,k)            = 1.0_R8B 
               lb_dom%u(LB_NODE(1,i,j,k))   = 0._R8B
               lb_dom%u(LB_NODE(2,i,j,k))   = 0._R8B
#ifdef D3Q19
               lb_dom%u(LB_NODE(3,i,j,k))   = 0._R8B
#endif
             end do
          end do
       end do
    else if(s_par%problem==shock) then
         write(*,*) "not implemented yet"

   end if
  end subroutine init_field_each
  !------------------------------------------------------------------------





  subroutine nrbc_find_neighbors(lb_dom,s_par,prc)
!------------------------------------------------------------------------
!
! subroutine to find next and next-next neighbor for non-reflecting boundary
!

   implicit none
   type(lb_block     )                 :: lb_dom
   type(sim_parameter)                 :: s_par
   type(mpl_var      )                 :: prc  
   integer                        :: i
   integer                        :: x,y,z 

   do i=1,lb_dom%nnr_wall   
      x = lb_dom%obs_nrwall(i,1)
      y = lb_dom%obs_nrwall(i,2)
      z = lb_dom%obs_nrwall(i,3)
      lb_dom%obs_nrwall(i,0) = 0  ! integer with directions (set by setting the bits)

      ! Check: Is Inlet Set?? If not, outlet. If yes, inlet
      if(btest(lb_dom%state(x,y,z),nr_wall_in) .eqv. .true.) then
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),0)
      endif
      ! x coordinate
      if(x == 1 .and. prc%crd(1)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,4) = x+1 
         lb_dom%obs_nrwall(i,7) = x+2
         lb_dom%obs_nrwall(i,10) = x+3
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),3)
      elseif(x == lb_dom%lx(1) .and. prc%crd(1)==prc%np(1)-1) then ! if on right border, use left neighbors
         lb_dom%obs_nrwall(i,4) = x-1
         lb_dom%obs_nrwall(i,7) = x-2
         lb_dom%obs_nrwall(i,10) = x-3
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),1)
      else ! if not, use same coordinate
         lb_dom%obs_nrwall(i,4) = x
         lb_dom%obs_nrwall(i,7) = x
         lb_dom%obs_nrwall(i,10) = x
      endif
      ! y coordinate
!FIXME original: not diagonal!      if(y == 1 .and. x>1  .and. x<lb_dom%lx(1)) then ! if on left boarder, use right neighbors
      if(y == 1 .and. prc%crd(2)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,5) = y+1 
         lb_dom%obs_nrwall(i,8) = y+2 
         lb_dom%obs_nrwall(i,11) = y+3 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),4)
!FIXME original: not diagnol      elseif(y == lb_dom%lx(2) .and. x>1 .and. x<lb_dom%lx(1)) then! if on right border, use left neighbors
      elseif(y == lb_dom%lx(2) .and. prc%crd(2)==prc%np(2)-1) then! if on right border, use left neighbors
         lb_dom%obs_nrwall(i,5) = y-1
         lb_dom%obs_nrwall(i,8) = y-2
         lb_dom%obs_nrwall(i,11) = y-3 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),2)
      else ! if not, use same coordinate
         lb_dom%obs_nrwall(i,5) = y
         lb_dom%obs_nrwall(i,8) = y
         lb_dom%obs_nrwall(i,11) = y
      endif
! this is for testing, if extrapolation from the same node is useful for stability reasons
!      if(x==1 .and. y==1) then
!         lb_dom%obs_nrwall(i,4)  = x
!         lb_dom%obs_nrwall(i,7)  = x
!         lb_dom%obs_nrwall(i,10) = x
!         lb_dom%obs_nrwall(i,5)  = y
!         lb_dom%obs_nrwall(i,8)  = y
!         lb_dom%obs_nrwall(i,11) = y
!      endif
#ifdef D3Q19
      ! z coordinate
!FIXME      if(z == 1.and. x>1 .and. x<lb_dom%lx(1).and. y>1 .and. y<lb_dom%lx(2)) then ! if on left boarder, use right neighbors
      if(z == 1.and.  prc%crd(3)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,6) = z+1 
         lb_dom%obs_nrwall(i,9) = z+2 
         lb_dom%obs_nrwall(i,12) = z+3 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),6)
!FIXME      elseif(z == lb_dom%lx(3).and. x>1 .and. x<lb_dom%lx(1).and. y>1 .and. y<lb_dom%lx(2)) then
      elseif(z == lb_dom%lx(3).and. prc%crd(3)==prc%np(3)-1 ) then
         lb_dom%obs_nrwall(i,6) = z-1
         lb_dom%obs_nrwall(i,9) = z-2
         lb_dom%obs_nrwall(i,12) = z-3 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),5)
      else ! if not, use same coordinate
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
      endif
#else
      ! if 2d, set third coord to 1
      lb_dom%obs_nrwall(i,6) = 1
      lb_dom%obs_nrwall(i,9) = 1
      lb_dom%obs_nrwall(i,12) = 1 
#endif
      ! assign zero values
      lb_dom%nrwall_0val(i,0)      = lb_dom%rho(x,y,z)
      lb_dom%nrwall_0val(i,1:NDIM) = lb_dom%u(LB_NODE(1:NDIM,x,y,z))
   end do
   end subroutine nrbc_find_neighbors
  !------------------------------------------------------------------------





  subroutine per_find_neighbors(lb_dom,s_par,prc)
!------------------------------------------------------------------------
!
! subroutine to find next and next-next neighbor for non-reflecting boundary
!

   implicit none
   type(lb_block     )      :: lb_dom
   type(sim_parameter)      :: s_par
   type(mpl_var      )      :: prc  
   integer                  :: i,j,k,l
   integer                  :: x,y,z 

   do i=1,lb_dom%nper_wall   
      x = lb_dom%obs_per(i,1)
      y = lb_dom%obs_per(i,2)
      z = lb_dom%obs_per(i,3)
      lb_dom%obs_per(i,0) = 0  ! integer with directions (set by setting the bits)
      ! x coordinate
      if(x == 1) then ! if on left boarder, use right neighbors
         lb_dom%obs_per(i,4) = lb_dom%lx(1)-1 
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),3)
      elseif(x == lb_dom%lx(1)) then ! if on right border, use left neighbors
         lb_dom%obs_per(i,4) = 2
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),1)
      else ! if not, use same coordinate
         lb_dom%obs_per(i,4) = x
      endif
      ! y coordinate
      if(y == 1) then ! if on left boarder, use right neighbors
         lb_dom%obs_per(i,5) = lb_dom%lx(2)-1 
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),4)
      elseif(y == lb_dom%lx(2)) then! if on right border, use left neighbors
         lb_dom%obs_per(i,5) = 2  
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),2)
      else ! if not, use same coordinate
         lb_dom%obs_per(i,5) = y
      endif
#ifdef D3Q19
      ! z coordinate
      if(z == 1) then ! if on left boarder, use right neighbors
         lb_dom%obs_per(i,6) = lb_dom%lx(3)-1 
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),6)
      elseif(z == lb_dom%lx(3)) then
         lb_dom%obs_per(i,6) = 2  
         lb_dom%obs_per(i,0) = ibset(lb_dom%obs_per(i,0),5)
      else ! if not, use same coordinate
         lb_dom%obs_per(i,6) = z
      endif
#else
      ! if 2d, set third coord to 1
      lb_dom%obs_per(i,6) = 1
#endif
   end do
   end subroutine per_find_neighbors
  !------------------------------------------------------------------------







   subroutine init_vector(lb_dom,s_par,prc)

!------------------------------------------------------------------------
! create a vector for the boundary conditions containing the coordinates 
! of the boundary node (and additionally first and second neighbors for nrbc)
!


      implicit none
      type(lb_block)                :: lb_dom
      type(sim_parameter)           :: s_par
      type(mpl_var)                 :: prc  
      integer      :: i,j,k,l
      integer      :: sum_nnr,sum_nper
      integer      :: counter,counter_sponge,counter_reset,counter_nrwall,counter_perwall

      lb_dom%nobs = 0
!      lb_dom%ninlet= 0
!      lb_dom%noutlet= 0
!      lb_dom%nobs_sponge = 0
!      lb_dom%nobs_reset = 0
      lb_dom%nnr_wall = 0
      lb_dom%nper_wall = 0

      ! count how many of each boundaries are present
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
                 lb_dom%nobs = lb_dom%nobs + 1
               end if
!               if(btest(lb_dom%state(i,j,k),sponge_wall)) then
!                 lb_dom%nobs_sponge = lb_dom%nobs_sponge + 1
!               end if
!               if(btest(lb_dom%state(i,j,k),reset_wall)) then
!                 lb_dom%nobs_reset = lb_dom%nobs_reset + 1
!               end if
               if(btest(lb_dom%state(i,j,k),nr_wall)) then
                 lb_dom%nnr_wall = lb_dom%nnr_wall + 1
               end if
               if(btest(lb_dom%state(i,j,k),per_wall)) then
                 lb_dom%nper_wall = lb_dom%nper_wall + 1
               end if
!               if(btest(lb_dom%state(i,j,k),inlet ) ) then
!                 lb_dom%ninlet  = lb_dom%ninlet   + 1
!               end if
!               if(btest(lb_dom%state(i,j,k),outlet) ) then
!                 lb_dom%noutlet  = lb_dom%noutlet  + 1
!               end if


            end do
         end do
      end do

      allocate(lb_dom%obs(lb_dom%nobs,3))
!      allocate(lb_dom%obs_sponge(lb_dom%nobs_sponge,9))
!      allocate(lb_dom%obs_reset(lb_dom%nobs_reset,9))
      allocate(lb_dom%obs_per(lb_dom%nper_wall,0:6))
!      allocate(lb_dom%per_val(lb_dom%nper_wall,nnod))
      allocate(lb_dom%obs_nrwall(lb_dom%nnr_wall,0:12))
               ! 0     ... holds bitmap information, which border
               ! 1..3  ... coordinates of nrbc node
               ! 4..12 ... hold 1st, 2nd, 3rd neighbor node coordinates 
      allocate(lb_dom%nrwall_0val(lb_dom%nnr_wall,0:3))  ! holds values of rho and u of bc
!      allocate(lb_dom%inlet(lb_dom%ninlet,3))
!      allocate(lb_dom%outlet(lb_dom%noutlet,3))

      counter        = 0
      counter_sponge = 0
      counter_reset  = 0
      counter_nrwall = 0
      counter_perwall= 0
!      counter_ninlet= 0
!      counter_noutlet=0
      
      ! then assign the boundaries with their coordinates
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
                  counter = counter+1
                  lb_dom%obs(counter,1:3) = (/i,j,k/)
               end if
!               if(btest(lb_dom%state(i,j,k),sponge_wall)) then
!                  counter_sponge = counter_sponge+1
!                  lb_dom%obs_sponge(counter_sponge,1:3) = (/i,j,k/)
!               end if
!               if(btest(lb_dom%state(i,j,k),reset_wall)) then
!                  counter_reset = counter_reset+1
!                  lb_dom%obs_reset(counter_reset,1:3) = (/i,j,k/)
!               end if
               if(btest(lb_dom%state(i,j,k),nr_wall)) then
                  counter_nrwall = counter_nrwall +1
                  lb_dom%obs_nrwall(counter_nrwall,1:3) = (/i,j,k/)
               end if
               if(btest(lb_dom%state(i,j,k),per_wall)) then
                  counter_perwall = counter_perwall +1
                  lb_dom%obs_per(counter_perwall,1:3) = (/i,j,k/)
               end if
!               if(btest(lb_dom%state(i,j,k),inlet ) ) then
!                  counter_ninlet= counter_ninlet+1
!                  lb_dom%inlet(counter_ninlet,1:3) = (/i,j,k/)
!               end if

            end do
         end do
      end do

#ifdef USE_MPI
       call MPI_REDUCE(lb_dom%nnr_wall, sum_nnr, 1,mpi_integer,MPI_SUM,prc%root_th,prc%cart_comm,prc%ierr)
       call MPI_REDUCE(lb_dom%nper_wall,sum_nper,1,mpi_integer,MPI_SUM,prc%root_th,prc%cart_comm,prc%ierr)
#else
       sum_nnr  = lb_dom%nnr_wall
       sum_nper = lb_dom%nper_wall
#endif

       if(prc%rk==0) then
          write(*,*) 
          write(*,*) "Setting Boundary conditions resulted in"
          write(*,'(a10,i6,a14,i6)') " # NRBC:  ",sum_nnr,"   # Periodic: ",sum_nper
          write(*,*) 
       endif
      if (lb_dom%nnr_wall > 0) then
         call nrbc_find_neighbors(lb_dom,s_par,prc)
      end if

      if (lb_dom%nper_wall > 0) then
         if(prc%size > 1) then 
            write(*,*) "Periodic functions are only working with serial execution!"
            stop
         endif
         call per_find_neighbors(lb_dom,s_par,prc)
      end if




#ifdef SPARSE /* SPARSE */
! this is a preparation for the usage of indirect addressing. 
! not used until now
   ! allocate the number of (nodes + ghost nodes - solid nodes)*nnod
   allocate(lb_dom%fsIn(((lb_dom%lx(1)+2*prc%n_ghost_nodes)*&
      & (lb_dom%lx(1)+2*prc%n_ghost_nodes)*(lb_dom%lx(1)+2*prc%n_ghost_nodes)-lb_dom%nobs)*nnod)
   ! assign the densities one by one
      counter = 0
      do k=1-prc%n_ghost_nodes,lb_dom%lx(3)+prc%n_ghost_nodes
         do j=1-prc%n_ghost_nodes,lb_dom%lx(2)+prc%n_ghost_nodes
            do i=1-prc%n_ghost_nodes,lb_dom%lx(1)+prc%n_ghost_nodes
               if(btest(lb_dom%state(i,j,k),fluid)) then
                  counter = counter+1
                   
               end if
            end do
         end do
      end do

   ! assign adjacency list
!FIXME ghost nodes also have to be defined as fluid nodes
      ! create the target vector
      counter = 0
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               counter = counter+1
               do l=1,nnod
                  lb_dom%f_sparse((counter*nnod-1)+l) = lb_dom%fIn(LB_NODE(l,i,j,k))
               end do
            end do
         end do
      end do
#endif

  end subroutine init_vector
  !------------------------------------------------------------------------
  





   subroutine calc_rho_ges(lb_dom,fIn,rho_ges,tStep,prc)
!------------------------------------------------------------------------
! calculate the total density of the domain
! This can be used for convergence checks



      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(mpl_var)              :: prc
      real(R8B),intent(in)       :: fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B),intent(out)      :: rho_ges
      integer,intent(in)         :: tStep  
      integer                    :: i,j,k,l
      real(R8B)                  :: sum_rho,uc(3),umax(3)
#ifdef USE_CAF
      real(R8B)                  :: rho_tmp[*]
#endif
      sum_rho  = 0.0d0
      rho_ges  = 0.0d0
      umax     = 0.0d0


      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
!               if(state(i,j,k) /= wall) then
               uc = 0.0d0
               do l=1,nnod
                  rho_ges=rho_ges+fIn(LB_NODE(l,i,j,k))
               enddo
#ifdef CALC_MAX_VELOCITY
                  uc(1) = cx(l,1)*fIn(LB_NODE(l,i,j,k)) 
                  uc(2) = cx(l,2)*fIn(LB_NODE(l,i,j,k))
                  uc(3) = cx(l,3)*fIn(LB_NODE(l,i,j,k)) 
               enddo
               if(abs(uc(1)) > umax(1)) then
                  umax(1)=abs(uc(1))
               endif
               if(abs(uc(2)) > umax(2)) then
                  umax(2)=abs(uc(2))
               endif
               if(abs(uc(3)) > umax(3)) then
                  umax(3)=abs(uc(3))
               endif
#endif /* CALC_MAX_VELOCITY */
               endif
            enddo
         enddo
      enddo
#ifdef USE_CAF
      rho_tmp = rho_ges 
#endif /* USE_CAF */

      if(prc%size > 0) then
#ifdef USE_MPI
         ! if parallel execution, sum up the calculated values and sen to process 0
         call MPI_REDUCE(rho_ges,sum_rho,1,mpi_double_precision,MPI_SUM,prc%root_th,prc%cart_comm,prc%ierr)
         rho_ges=sum_rho
#endif /* MPI */
#ifdef USE_CAF
      sync all
      if(prc%rk==0) then
      do i=1,prc%size
        sum_rho = sum_rho + rho_tmp[i] 
      enddo
      endif
#endif /* USE_CAF */
      end if
      if(prc%rk == 0) then
#ifdef CALC_MAX_VELOCITY
         write(*,'(a16,f10.5,a4,f10.5,a4,f10.5)') "max velocity  x:",umax(1),"  y:",umax(2),"  z:",umax(3)
#endif /* CALC_MAX_VELOCITY */
         write(*,'(a16,f16.6,a12,i7)') "Total density:  ",rho_ges," at timestep ",tStep
         open(43,file='dens.inf',position='append')
         write(43,*) tStep,rho_ges
         close(43)
      end if


   end subroutine calc_rho_ges
  !------------------------------------------------------------------------
  







   subroutine propagate(lb_dom,fOut,fIn,state,meas)

!------------------------------------------------------------------------
! propagation step: move the densities to the neighbor nodes
! performs PULL of densities


      implicit none
      type(lb_block),intent(inout)  :: lb_dom
      type(measure )            :: meas  
      real(R8B), intent(in)     :: fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)  :: fOut(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer, intent(in)       :: state(lb_dom%lx(1),lb_dom%lx(2),lb_dom%lx(3))
      integer                   :: i,j,k,l,count_err

      call cpu_time_measure(meas%tSt_comm) 
      count_err=0


      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)

                  fOut(LB_NODE(1,i,j,k))=fIn(LB_NODE(1,i-cx(1,1),j-cx(1,2),k-cx(1,3)))
                  fOut(LB_NODE(2,i,j,k))=fIn(LB_NODE(2,i-cx(2,1),j-cx(2,2),k-cx(2,3)))
                  fOut(LB_NODE(3,i,j,k))=fIn(LB_NODE(3,i-cx(3,1),j-cx(3,2),k-cx(3,3)))
                  fOut(LB_NODE(4,i,j,k))=fIn(LB_NODE(4,i-cx(4,1),j-cx(4,2),k-cx(4,3)))
                  fOut(LB_NODE(5,i,j,k))=fIn(LB_NODE(5,i-cx(5,1),j-cx(5,2),k-cx(5,3)))
                  fOut(LB_NODE(6,i,j,k))=fIn(LB_NODE(6,i-cx(6,1),j-cx(6,2),k-cx(6,3)))
                  fOut(LB_NODE(7,i,j,k))=fIn(LB_NODE(7,i-cx(7,1),j-cx(7,2),k-cx(7,3)))
                  fOut(LB_NODE(8,i,j,k))=fIn(LB_NODE(8,i-cx(8,1),j-cx(8,2),k-cx(8,3)))
                  fOut(LB_NODE(9,i,j,k))=fIn(LB_NODE(9,i-cx(9,1),j-cx(9,2),k-cx(9,3)))
#ifdef D3Q19
                  fOut(LB_NODE(10,i,j,k))=fIn(LB_NODE(10,i-cx(10,1),j-cx(10,2),k-cx(10,3)))
                  fOut(LB_NODE(11,i,j,k))=fIn(LB_NODE(11,i-cx(11,1),j-cx(11,2),k-cx(11,3)))
                  fOut(LB_NODE(12,i,j,k))=fIn(LB_NODE(12,i-cx(12,1),j-cx(12,2),k-cx(12,3)))
                  fOut(LB_NODE(13,i,j,k))=fIn(LB_NODE(13,i-cx(13,1),j-cx(13,2),k-cx(13,3)))
                  fOut(LB_NODE(14,i,j,k))=fIn(LB_NODE(14,i-cx(14,1),j-cx(14,2),k-cx(14,3)))
                  fOut(LB_NODE(15,i,j,k))=fIn(LB_NODE(15,i-cx(15,1),j-cx(15,2),k-cx(15,3)))
                  fOut(LB_NODE(16,i,j,k))=fIn(LB_NODE(16,i-cx(16,1),j-cx(16,2),k-cx(16,3)))
                  fOut(LB_NODE(17,i,j,k))=fIn(LB_NODE(17,i-cx(17,1),j-cx(17,2),k-cx(17,3)))
                  fOut(LB_NODE(18,i,j,k))=fIn(LB_NODE(18,i-cx(18,1),j-cx(18,2),k-cx(18,3)))
                  fOut(LB_NODE(19,i,j,k))=fIn(LB_NODE(19,i-cx(19,1),j-cx(19,2),k-cx(19,3)))
#endif
            enddo
         enddo
      enddo

      call cpu_time_measure(meas%tEnd_comm)
      meas%prop_duration = meas%tEnd_comm - meas%tSt_comm + meas%prop_duration

   end subroutine propagate
  !------------------------------------------------------------------------







#ifdef TRT

  subroutine stream_collide_trt(lb_dom,fOut,fIn,omega)


!------------------------------------------------------------------------
! Two-relaxation time collision routine



    implicit none
    type(lb_block),intent(inout)  :: lb_dom
    real(R8B) , intent(in)  :: omega
    real(R8B), intent(out)     :: fOut(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
    real(R8B), intent(inout)   :: fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))

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

                  ftmp(1)=fIn(LB_NODE(1,i-cx(1,1),j-cx(1,2),k-cx(1,3)))
                  ftmp(2)=fIn(LB_NODE(2,i-cx(2,1),j-cx(2,2),k-cx(2,3)))
                  ftmp(3)=fIn(LB_NODE(3,i-cx(3,1),j-cx(3,2),k-cx(3,3)))
                  ftmp(4)=fIn(LB_NODE(4,i-cx(4,1),j-cx(4,2),k-cx(4,3)))
                  ftmp(5)=fIn(LB_NODE(5,i-cx(5,1),j-cx(5,2),k-cx(5,3)))
                  ftmp(6)=fIn(LB_NODE(6,i-cx(6,1),j-cx(6,2),k-cx(6,3)))
                  ftmp(7)=fIn(LB_NODE(7,i-cx(7,1),j-cx(7,2),k-cx(7,3)))
                  ftmp(8)=fIn(LB_NODE(8,i-cx(8,1),j-cx(8,2),k-cx(8,3)))
                  ftmp(9)=fIn(LB_NODE(9,i-cx(9,1),j-cx(9,2),k-cx(9,3)))
#ifdef D3Q19
                  ftmp(10)=fIn(LB_NODE(10,i-cx(10,1),j-cx(10,2),k-cx(10,3)))
                  ftmp(11)=fIn(LB_NODE(11,i-cx(11,1),j-cx(11,2),k-cx(11,3)))
                  ftmp(12)=fIn(LB_NODE(12,i-cx(12,1),j-cx(12,2),k-cx(12,3)))
                  ftmp(13)=fIn(LB_NODE(13,i-cx(13,1),j-cx(13,2),k-cx(13,3)))
                  ftmp(14)=fIn(LB_NODE(14,i-cx(14,1),j-cx(14,2),k-cx(14,3)))
                  ftmp(15)=fIn(LB_NODE(15,i-cx(15,1),j-cx(15,2),k-cx(15,3)))
                  ftmp(16)=fIn(LB_NODE(16,i-cx(16,1),j-cx(16,2),k-cx(16,3)))
                  ftmp(17)=fIn(LB_NODE(17,i-cx(17,1),j-cx(17,2),k-cx(17,3)))
                  ftmp(18)=fIn(LB_NODE(18,i-cx(18,1),j-cx(18,2),k-cx(18,3)))
                  ftmp(19)=fIn(LB_NODE(19,i-cx(19,1),j-cx(19,2),k-cx(19,3)))
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
              ftmp(12) - fIn(LB_NODE(13,i,j,k) )- ftmp(14)  + & 
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

       fOut(LB_NODE(1,i,j,k)) = ftmp(1)*(1.d0-omega) + omega*t0*feq_common

#ifdef D3Q19
       t2x2 = t2x2_0 * loc_dens
       fac2 = t2x2 * inv2csq2

       ui   = u_x + u_y
       sym  = omega_h*(ftmp(8) + ftmp(9) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(8) - ftmp(9) - 3.d0*t2x2*ui )
       fOut(LB_NODE(8,i,j,k))  = ftmp(8) - sym - asym
       fOut(LB_NODE(9,i,j,k))  = ftmp(9) - sym + asym

       ui   = u_x - u_y
       sym  = omega_h*(ftmp(10) + ftmp(11) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(10) - ftmp(11) - 3.d0*t2x2*ui )
       fOut(LB_NODE(10,i,j,k)) = ftmp(10) - sym - asym
       fOut(LB_NODE(11,i,j,k)) = ftmp(11) - sym + asym

       ui   = u_x + u_z
       sym  = omega_h*(ftmp(12) + ftmp(13) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(12) - ftmp(13) - 3.d0*t2x2*ui )
       fOut(LB_NODE(12,i,j,k)) = ftmp(12) - sym - asym
       fOut(LB_NODE(13,i,j,k)) = ftmp(13) - sym + asym

       ui   = u_x - u_z
       sym  = omega_h*(ftmp(14) + ftmp(15) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(14) - ftmp(15) - 3.d0*t2x2*ui )
       fOut(LB_NODE(14,i,j,k)) = ftmp(14) - sym - asym
       fOut(LB_NODE(15,i,j,k)) = ftmp(15) - sym + asym

       ui   = u_y + u_z
       sym  = omega_h*(ftmp(16) + ftmp(17) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(16) - ftmp(17) - 3.d0*t2x2*ui )
       fOut(LB_NODE(16,i,j,k)) = ftmp(16) - sym - asym
       fOut(LB_NODE(17,i,j,k)) = ftmp(17) - sym + asym

       ui   = u_y - u_z
       sym  = omega_h*(ftmp(18) + ftmp(19) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(18) - ftmp(19) - 3.d0*t2x2*ui )
       fOut(LB_NODE(18,i,j,k)) = ftmp(18) - sym - asym
       fOut(LB_NODE(19,i,j,k)) = ftmp(19) - sym + asym

       t1x2 = t1x2_0 * loc_dens
       fac1 = t1x2 * inv2csq2

       ui   = u_y
       sym  = omega_h*(ftmp(4) + ftmp(5) - fac1*ui*ui - t1x2*feq_common)
       asym = asym_omega_h*( ftmp(4) - ftmp(5) - 3.d0*t1x2*ui )
       fOut(LB_NODE(4,i,j,k)) = ftmp(4) - sym - asym
       fOut(LB_NODE(5,i,j,k)) = ftmp(5) - sym + asym

       ui   = u_x
       sym  = omega_h*(ftmp(2) + ftmp(3) - fac1*ui*ui - t1x2*feq_common)
       asym = asym_omega_h*( ftmp(2) - ftmp(3) - 3.d0*t1x2*ui )
       fOut(LB_NODE(2,i,j,k)) = ftmp(2) - sym - asym
       fOut(LB_NODE(3,i,j,k)) = ftmp(3) - sym + asym

       ui   = u_z
       sym  = omega_h*(ftmp(6) + ftmp(7)  - fac1*ui*ui - t1x2*feq_common)
       asym = asym_omega_h*( ftmp(6) - ftmp(7) - 3.d0*t1x2*ui )
       fOut(LB_NODE(6,i,j,k)) = ftmp(6) - sym - asym
       fOut(LB_NODE(7,i,j,k)) = ftmp(7) - sym + asym
#endif
#ifdef D2Q9
       t1x2 = t1x2_0 * loc_dens
       fac1 = t1x2 * inv2csq2

       ! y-direction densities
       ui   = u_y
       sym  = omega_h*(ftmp(3) + ftmp(5) - fac1*ui*ui - t1x2*feq_common)
       asym = asym_omega_h*( ftmp(3) - ftmp(5) - 3.d0*t1x2*ui )
       fOut(LB_NODE(3,i,j,k)) = ftmp(3) - sym - asym
       fOut(LB_NODE(5,i,j,k)) = ftmp(5) - sym + asym

       ! x-direction densities
       ui   = u_x
       sym  = omega_h*(ftmp(2) + ftmp(4) - fac1*ui*ui - t1x2*feq_common)
       asym = asym_omega_h*( ftmp(2) - ftmp(4) - 3.d0*t1x2*ui )
       fOut(LB_NODE(2,i,j,k)) = ftmp(2) - sym - asym
       fOut(LB_NODE(4,i,j,k)) = ftmp(4) - sym + asym

       ! xy-direction densities
       t2x2 = t2x2_0 * loc_dens
       fac2 = t2x2 * inv2csq2

       ui   = u_x + u_y
       sym  = omega_h*(ftmp(6) + ftmp(8) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(6) - ftmp(8) - 3.d0*t2x2*ui )
       fOut(LB_NODE(6,i,j,k))  = ftmp(6) - sym - asym
       fOut(LB_NODE(8,i,j,k))  = ftmp(8) - sym + asym

       ui   = u_x - u_y
       sym  = omega_h*(ftmp(9) + ftmp(7) - fac2*ui*ui - t2x2*feq_common)
       asym = asym_omega_h*( ftmp(9) - ftmp(7) - 3.d0*t2x2*ui )
       fOut(LB_NODE(9,i,j,k)) = ftmp(9) - sym - asym
       fOut(LB_NODE(7,i,j,k)) = ftmp(7) - sym + asym
#endif
      else

         !--------------------------------------------
         ! Do Bounce Back 

         fOut(LB_NODE(1,i,j,k)) = ftmp(opp(1))
         fOut(LB_NODE(2,i,j,k)) = ftmp(opp(2))
         fOut(LB_NODE(3,i,j,k)) = ftmp(opp(3))
         fOut(LB_NODE(4,i,j,k)) = ftmp(opp(4))
         fOut(LB_NODE(5,i,j,k)) = ftmp(opp(5))
         fOut(LB_NODE(6,i,j,k)) = ftmp(opp(6))
         fOut(LB_NODE(7,i,j,k)) = ftmp(opp(7))
         fOut(LB_NODE(8,i,j,k)) = ftmp(opp(8))
         fOut(LB_NODE(9,i,j,k)) = ftmp(opp(9))
#ifdef D3Q19
         fOut(LB_NODE(10,i,j,k)) = ftmp(opp(10))
         fOut(LB_NODE(11,i,j,k)) = ftmp(opp(11))
         fOut(LB_NODE(12,i,j,k)) = ftmp(opp(12))
         fOut(LB_NODE(13,i,j,k)) = ftmp(opp(13))
         fOut(LB_NODE(14,i,j,k)) = ftmp(opp(14))
         fOut(LB_NODE(15,i,j,k)) = ftmp(opp(15))
         fOut(LB_NODE(16,i,j,k)) = ftmp(opp(16))
         fOut(LB_NODE(17,i,j,k)) = ftmp(opp(17))
         fOut(LB_NODE(18,i,j,k)) = ftmp(opp(18))
         fOut(LB_NODE(19,i,j,k)) = ftmp(opp(19))
#endif

      endif
          end do
       end do
    end do
#endif /* VECTOR */

end subroutine stream_collide_trt
#endif /* TRT */

#ifdef MRT_NOT_READY_YET
   subroutine collide_mrt(lb_dom,fOut,fIn,s_par,meas)

!------------------------------------------------------------------------
! This is the MRT collision step
! Note:
! In the paper of Mei and Luo 2006, drho = rho-rho0 is used instead of rho to reduce roundoff errors.
! This should be implemented here also.
!
! !!! This routine is NOT optimized for performance
! See EXCEL-Sheet BGK_MRT.xls 

      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure)              :: meas  
      real(R8B)                  :: usq,omega
      real(R8B)                  :: meq(nnod),mneq(nnod),mout(nnod),s_mrt(nnod)
      real(R8B)                  :: M_tr(nnod,nnod),M_inv(nnod,nnod),fOut_tmp
      real(R8B)                  :: weps,wepsj,wxx,rholoc_inv 
      real(R8B)                  :: u(3),rholoc,fEq(nnod),t2cs4inv,t2cs2inv,jx(3)
      real(R8B)                  :: rho0 
      real(R8B), intent(out)     :: fOut(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)   :: fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer                    :: i,j,k,l,ll

    real*8   dd_Q19_0,dd_Q19_NE,dd_Q19_N,dd_Q19_NW,dd_Q19_W,        &
             dd_Q19_SW,dd_Q19_S,dd_Q19_SE,dd_Q19_E,dd_Q19_T,        &
             dd_Q19_TE,dd_Q19_TN,dd_Q19_TW,dd_Q19_TS,dd_Q19_B,      &
             dd_Q19_BE,dd_Q19_BN,dd_Q19_BW,dd_Q19_BS
    real*8 u_x,u_y,u_z
    real*8 j_x,j_y,j_z
    real*8 loc_dens, inv_loc_dens

    real*8, parameter :: one_over__2 = 0.5d0
    real*8, parameter :: one_over__3 = 0.333333333333333333333d0
    real*8, parameter :: two_over__3 = 0.666666666666666666666d0
    real*8, parameter :: one_over__4 = 0.25d0
    real*8, parameter :: one_over__6 = 0.166666666666666666667d0
    real*8, parameter :: one_over__8 = 0.125d0
    real*8, parameter :: one_over__9 = 0.111111111111111111111d0
    real*8, parameter :: one_over_12 = 0.083333333333333333333d0
    real*8, parameter :: one_over_16 = 0.0625d0
    real*8, parameter :: one_over_18 = 0.055555555555555555556d0
    real*8, parameter :: one_over_19 = 0.052631578947368421053d0
    real*8, parameter :: one_over_21 = 0.047619047619047619048d0
    real*8, parameter :: one_over_24 = 0.041666666666666666667d0
    real*8, parameter :: one_over_36 = 0.027777777777777777778d0
    real*8, parameter :: one_over_48 = 0.020833333333333333333d0
    real*8, parameter :: one_over_72 = 0.013888888888888888889d0

      call cpu_time_measure(meas%tSt_comm) 

      omega = s_par%omega
      if(s_par%initial) omega = 1.1

  s_mrt(1) = 0._R8B
  s_mrt(2) = 1.19_R8B
  s_mrt(3) = 1.4_R8B
  s_mrt(4) = 0._R8B
  s_mrt(5) = 1.2_R8B
  s_mrt(6) = 0._R8B
  s_mrt(7) = 1.2_R8B
  s_mrt(8) = 0._R8B
  s_mrt(9) = 1.2_R8B
  s_mrt(10)= omega
  s_mrt(11)= 1.4_R8B
  s_mrt(12)= omega
  s_mrt(13)= 1.4_R8B
  s_mrt(14)= omega
  s_mrt(15)= omega
  s_mrt(16)= omega
  s_mrt(17)= 1.98_R8B
  s_mrt(18)= 1.98_R8B
  s_mrt(19)= 1.98_R8B

      if(s_par%initial) then
         s_mrt(4) = 1._R8B
         s_mrt(6) = 1._R8B
         s_mrt(8) = 1._R8B
      endif

!#ifdef CHECK_MRT
!FIXME These are BGK parameters. With optimized parameters, strange effects in cavity. 
!  See also Excel sheet with these parameters: rest-denstiy somehow decreases strongly in comparison to other parameter settings
      weps  = 3._R8B 
      wepsj = -5.5_R8B
      wxx   = -0.5_R8B
!#else
!      weps  = 0._R8B 
!      wepsj = -475._R8B/63._R8B
!      wxx   = 0._R8B
!#endif /* CHECK_MRT */


#ifdef CHECK_MRT
!     s_mrt(1:nnod) = omega 
#endif /* CHECK_MRT */


 

      ! Multiple Single Relaxation time Collision routine
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then

dd_Q19_0   =  fIn(LB_NODE(1 ,i,j,k))
dd_Q19_E   =  fIn(LB_NODE(2 ,i,j,k))
dd_Q19_W   =  fIn(LB_NODE(3 ,i,j,k))
dd_Q19_N   =  fIn(LB_NODE(4 ,i,j,k))
dd_Q19_S   =  fIn(LB_NODE(5 ,i,j,k))
dd_Q19_T   =  fIn(LB_NODE(6 ,i,j,k))
dd_Q19_B   =  fIn(LB_NODE(7 ,i,j,k))
dd_Q19_NE  =  fIn(LB_NODE(8 ,i,j,k))
dd_Q19_SW  =  fIn(LB_NODE(9 ,i,j,k))
dd_Q19_SE  =  fIn(LB_NODE(10,i,j,k))
dd_Q19_NW  =  fIn(LB_NODE(11,i,j,k))
dd_Q19_TE  =  fIn(LB_NODE(12,i,j,k))
dd_Q19_BW  =  fIn(LB_NODE(13,i,j,k))
dd_Q19_BE  =  fIn(LB_NODE(14,i,j,k))
dd_Q19_TW  =  fIn(LB_NODE(15,i,j,k))
dd_Q19_TN  =  fIn(LB_NODE(16,i,j,k))
dd_Q19_BS  =  fIn(LB_NODE(17,i,j,k))
dd_Q19_BN  =  fIn(LB_NODE(18,i,j,k))
dd_Q19_TS  =  fIn(LB_NODE(19,i,j,k))



       loc_dens = dd_Q19_0                                                    &
            + dd_Q19_NE  + dd_Q19_N  + dd_Q19_NW  + dd_Q19_W                     &
            + dd_Q19_SW  + dd_Q19_S  + dd_Q19_SE  + dd_Q19_E                     &
            + dd_Q19_T  + dd_Q19_TE + dd_Q19_TN + dd_Q19_TW                    &
            + dd_Q19_TS + dd_Q19_B + dd_Q19_BE + dd_Q19_BN                    &
            + dd_Q19_BW + dd_Q19_BS
              if(s_par%initial .eqv. .false.) then

       ! local x-momentum (velocity)
       j_x = (dd_Q19_NE+dd_Q19_SE+dd_Q19_E+dd_Q19_TE+dd_Q19_BE                   &
            -dd_Q19_NW-dd_Q19_W-dd_Q19_SW-dd_Q19_TW-dd_Q19_BW)

       ! local y-momentum (velocity)
       j_y = (dd_Q19_NE+dd_Q19_N+dd_Q19_NW+dd_Q19_BN+dd_Q19_TN                   &
            -dd_Q19_SW-dd_Q19_S-dd_Q19_SE-dd_Q19_TS-dd_Q19_BS)

       ! local z-momentum (velocity)
       j_z = (dd_Q19_T +dd_Q19_TE+dd_Q19_TN+dd_Q19_TW+dd_Q19_TS                &
            -dd_Q19_B-dd_Q19_BE-dd_Q19_BN-dd_Q19_BW-dd_Q19_BS)



               else
                 j_x = lb_dom%u0(LB_NODE(1,i,j,k))*loc_dens
                 j_y = lb_dom%u0(LB_NODE(2,i,j,k))*loc_dens
                 j_z = lb_dom%u0(LB_NODE(3,i,j,k))*loc_dens
               endif

               inv_loc_dens = 1.d0/loc_dens
   
       u_x = j_x * inv_loc_dens
       u_y = j_y * inv_loc_dens
       u_z = j_z * inv_loc_dens

       ! square velocity and derived constants
       usq  = u_x*u_x + u_y*u_y + u_z*u_z




! This is the multiplication m_neq = M*f
! rho
mneq(1) = loc_dens 
mneq(2) = -30.d0*dd_Q19_0 -11.0d0*( dd_Q19_E + dd_Q19_W   + dd_Q19_N + dd_Q19_S +dd_Q19_T + dd_Q19_B) &
                           + 8.0d0*( dd_Q19_NE + dd_Q19_SW + dd_Q19_SE + dd_Q19_NW & 
                                   + dd_Q19_TE + dd_Q19_BW + dd_Q19_BE + dd_Q19_TW & 
                                   + dd_Q19_TN + dd_Q19_BS + dd_Q19_BN + dd_Q19_TS) 
mneq(3) =  12.d0*dd_Q19_0  -4.0d0*( dd_Q19_E + dd_Q19_W   + dd_Q19_N + dd_Q19_S +dd_Q19_T + dd_Q19_B) &
                                   + dd_Q19_NE + dd_Q19_SW + dd_Q19_SE + dd_Q19_NW & 
                                   + dd_Q19_TE + dd_Q19_BW + dd_Q19_BE + dd_Q19_TW & 
                                   + dd_Q19_TN + dd_Q19_BS + dd_Q19_BN + dd_Q19_TS 
! j_x 
mneq(4) = j_x 
! q_x 
mneq(5) =                    4.0d0*(-dd_Q19_E + dd_Q19_W)                           &
                                   + dd_Q19_NE - dd_Q19_SW + dd_Q19_SE - dd_Q19_NW & 
                                   + dd_Q19_TE - dd_Q19_BW + dd_Q19_BE - dd_Q19_TW
! j_y 
mneq(6) = j_y 
! q_y 
mneq(7) =                    4.0d0*(-dd_Q19_N + dd_Q19_S) &
                                   + dd_Q19_NE - dd_Q19_SW - dd_Q19_SE + dd_Q19_NW & 
                                   + dd_Q19_TN - dd_Q19_BS + dd_Q19_BN - dd_Q19_TS 

! j_z 
mneq(8) = j_z
! q_z 
mneq(9) =                    4.0d0*(-dd_Q19_T + dd_Q19_B) &
                                   + dd_Q19_TE - dd_Q19_BW - dd_Q19_BE + dd_Q19_TW &
                                   + dd_Q19_TN - dd_Q19_BS - dd_Q19_BN + dd_Q19_TS 
! p_xx
mneq(10)=                    2.0d0*( dd_Q19_E + dd_Q19_W) - dd_Q19_N  - dd_Q19_S - dd_Q19_T - dd_Q19_B &
                                   + dd_Q19_NE + dd_Q19_SW + dd_Q19_SE + dd_Q19_NW & 
                                   + dd_Q19_TE + dd_Q19_BW + dd_Q19_BE + dd_Q19_TW &
                            -2.0d0*( dd_Q19_TN + dd_Q19_BS + dd_Q19_BN + dd_Q19_TS)
! pi_xx
mneq(11)=                   -4.0d0*( dd_Q19_E + dd_Q19_W) +2.0d0*(dd_Q19_N  + dd_Q19_S + dd_Q19_T + dd_Q19_B) &
                                   + dd_Q19_NE + dd_Q19_SW + dd_Q19_SE + dd_Q19_NW & 
                                   + dd_Q19_TE + dd_Q19_BW + dd_Q19_BE + dd_Q19_TW &
                            -2.0d0*( dd_Q19_TN + dd_Q19_BS + dd_Q19_BN + dd_Q19_TS)

! p_ww
mneq(12)=                            dd_Q19_N  + dd_Q19_S  - dd_Q19_T  - dd_Q19_B &
                                   + dd_Q19_NE + dd_Q19_SW + dd_Q19_SE + dd_Q19_NW & 
                                   - dd_Q19_TE - dd_Q19_BW - dd_Q19_BE - dd_Q19_TW  

! pi_ww
mneq(13)=                    2.0d0*(-dd_Q19_N  - dd_Q19_S  + dd_Q19_T  + dd_Q19_B) &
                                   + dd_Q19_NE + dd_Q19_SW + dd_Q19_SE + dd_Q19_NW & 
                                   - dd_Q19_TE - dd_Q19_BW - dd_Q19_BE - dd_Q19_TW  

! p_xy 
mneq(14)=               & 
                                     dd_Q19_NE + dd_Q19_SW - dd_Q19_SE - dd_Q19_NW 
! p_yz 
mneq(15)=               & 
                                     dd_Q19_TN + dd_Q19_BS - dd_Q19_BN - dd_Q19_TS 
! p_xz 
mneq(16)=               & 
                                     dd_Q19_TE + dd_Q19_BW - dd_Q19_BE - dd_Q19_TW
! m_x 
mneq(17)=               & 
                                     dd_Q19_NE - dd_Q19_SW + dd_Q19_SE - dd_Q19_NW & 
                                   - dd_Q19_TE + dd_Q19_BW - dd_Q19_BE + dd_Q19_TW
! m_y 
mneq(18)=               & 
                                   - dd_Q19_NE + dd_Q19_SW + dd_Q19_SE - dd_Q19_NW & 
                                   + dd_Q19_TN - dd_Q19_BS + dd_Q19_BN - dd_Q19_TS 
! m_z 
mneq(19)=               & 
                                   + dd_Q19_TE - dd_Q19_BW - dd_Q19_BE + dd_Q19_TW &
                                   - dd_Q19_TN + dd_Q19_BS + dd_Q19_BN - dd_Q19_TS 






! For compressible assumption, rho0 = loc_dens

                  meq( 1) = loc_dens
                  meq( 2) = -11.0d0*loc_dens + 19.0d0*inv_loc_dens*(j_x*j_x+j_y*j_y+j_z*j_z)
                  meq( 3) = weps*loc_dens+wepsj*inv_loc_dens*((j_x)*(j_x) + &
                                            (j_y)*(j_y)+(j_z)*(j_z))
                  meq( 4) =  j_x
                  meq( 6) =  j_y
                  meq( 8) =  j_z
! qx_eq, qy_eq, qz_eq
                  meq( 5) =  - two_over__3*j_x
                  meq( 7) =  - two_over__3*j_y
                  meq( 9) =  - two_over__3*j_z
! p_xx, pi_xx
                  meq(10) = inv_loc_dens*(2.d0*j_x*j_x - j_y*j_y - j_z*j_z)
                  meq(11) = wxx*meq(10) 
! p_ww, pi_ww
                  meq(12) = inv_loc_dens*(j_x*j_x - j_z*j_z)
                  meq(13) = wxx*meq(12) 
!p_xy,p_yz,p_xz
                  meq(14) =  inv_loc_dens*j_x*j_y
                  meq(15) =  inv_loc_dens*j_y*j_z
                  meq(16) =  inv_loc_dens*j_x*j_z
                  meq(17) =  0.0d0 
                  meq(18) =  0.0d0
                  meq(19) =  0.0d0

! Relax to eq
mout( 1) = mneq( 1)                                        !   mout( 1) = mneq( 1) - s_mrt( 1)*(mneq( 1) - meq( 1))   !  
mout( 2) = mneq( 2) - s_mrt( 2)*(mneq( 2) - meq( 2))       !   mout( 2) = mneq( 2) - s_mrt( 2)*(mneq( 2) - meq( 2))   !  
mout( 3) = mneq( 3) - s_mrt( 3)*(mneq( 3) - meq( 3))       !   mout( 3) = mneq( 3) - s_mrt( 3)*(mneq( 3) - meq( 3))   !  
mout( 4) = mneq( 4)                                        !   mout( 4) = mneq( 4) - s_mrt( 4)*(mneq( 4) - meq( 4))   !  
mout( 5) = mneq( 5) - s_mrt( 5)*(mneq( 5) - meq( 5))       !   mout( 5) = mneq( 5) - s_mrt( 5)*(mneq( 5) - meq( 5))   !  
mout( 6) = mneq( 6)                                        !   mout( 6) = mneq( 6) - s_mrt( 6)*(mneq( 6) - meq( 6))   !  
mout( 7) = mneq( 7) - s_mrt( 7)*(mneq( 7) - meq( 7))       !   mout( 7) = mneq( 7) - s_mrt( 7)*(mneq( 7) - meq( 7))   !  
mout( 8) = mneq( 8)                                        !   mout( 8) = mneq( 8) - s_mrt( 8)*(mneq( 8) - meq( 8))   !  
mout( 9) = mneq( 9) - s_mrt( 9)*(mneq( 9) - meq( 9))       !   mout( 9) = mneq( 9) - s_mrt( 9)*(mneq( 9) - meq( 9))   !  
mout(10) = mneq(10) - s_mrt(10)*(mneq(10) - meq(10))       !   mout(10) = mneq(10) - s_mrt(10)*(mneq(10) - meq(10))   !  
mout(11) = mneq(11) - s_mrt(11)*(mneq(11) - meq(11))       !   mout(11) = mneq(11) - s_mrt(11)*(mneq(11) - meq(11))   !  
mout(12) = mneq(12) - s_mrt(12)*(mneq(12) - meq(12))       !   mout(12) = mneq(12) - s_mrt(12)*(mneq(12) - meq(12))   !  
mout(13) = mneq(13) - s_mrt(13)*(mneq(13) - meq(13))       !   mout(13) = mneq(13) - s_mrt(13)*(mneq(13) - meq(13))   !  
mout(14) = mneq(14) - s_mrt(14)*(mneq(14) - meq(14))       !   mout(14) = mneq(14) - s_mrt(14)*(mneq(14) - meq(14))   !  
mout(15) = mneq(15) - s_mrt(15)*(mneq(15) - meq(15))       !   mout(15) = mneq(15) - s_mrt(15)*(mneq(15) - meq(15))   !  
mout(16) = mneq(16) - s_mrt(16)*(mneq(16) - meq(16))       !   mout(16) = mneq(16) - s_mrt(16)*(mneq(16) - meq(16))   !  
mout(17) = mneq(17) - s_mrt(17)*(mneq(17)          )       !   mout(17) = mneq(17) - s_mrt(17)*(mneq(17) - meq(17))   !  
mout(18) = mneq(18) - s_mrt(18)*(mneq(18)          )       !   mout(18) = mneq(18) - s_mrt(18)*(mneq(18) - meq(18))   !  
mout(19) = mneq(19) - s_mrt(19)*(mneq(19)          )       !   mout(19) = mneq(19) - s_mrt(19)*(mneq(19) - meq(19))   !  

! Transform back to pdf-space
      dd_Q19_0  = one_over_19*mout(1) -  5.d0/399.d0*mout(2) +one_over_21*mout(3)

      dd_Q19_E  = one_over_19*mout(1) - 11.d0/2394.d0*mout(2) -1.d0/63.d0*mout(3)         &
            +1.d0/10.d0*mout(4) -1.d0/10.d0*mout(5) + 1.d0/18.d0*mout(10) - 1.d0/18.d0*mout(11)

      dd_Q19_W  = one_over_19*mout(1) - 11.d0/2394.d0*mout(2) -1.d0/63.d0*mout(3)         &
            -1.d0/10.d0*mout(4) +1.d0/10.d0*mout(5) + 1.d0/18.d0*mout(10) - 1.d0/18.d0*mout(11)

      dd_Q19_N  = one_over_19*mout(1) - 11.d0/2394.d0*mout(2) -1.d0/63.d0*mout(3)         &
            +1.d0/10.d0*mout(6)  - 1.d0/10.d0*mout(7) - 1.d0/36.d0*mout(10) + 1.d0/36.d0*mout(11)  &
            +1.d0/12.d0*mout(12) - 1.d0/12.d0*mout(13)

      dd_Q19_S  = one_over_19*mout(1) - 11.d0/2394.d0*mout(2) -1.d0/63.d0*mout(3)         &
            -1.d0/10.d0*mout(6)  + 1.d0/10.d0*mout(7) - 1.d0/36.d0*mout(10) + 1.d0/36.d0*mout(11)  &
            +1.d0/12.d0*mout(12) - 1.d0/12.d0*mout(13)

      dd_Q19_T  = one_over_19*mout(1) - 11.d0/2394.d0*mout(2) -1.d0/63.d0*mout(3)         &
            +1.d0/10.d0*mout(8)  - 1.d0/10.d0*mout(9) - 1.d0/36.d0*mout(10) + 1.d0/36.d0*mout(11)  &
            -1.d0/12.d0*mout(12) + 1.d0/12.d0*mout(13)

      dd_Q19_B  = one_over_19*mout(1) - 11.d0/2394.d0*mout(2) -1.d0/63.d0*mout(3)         &
            -1.d0/10.d0*mout(8)  + 1.d0/10.d0*mout(9) - 1.d0/36.d0*mout(10) + 1.d0/36.d0*mout(11)  &
            -1.d0/12.d0*mout(12) + 1.d0/12.d0*mout(13)

      dd_Q19_NE = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
+ 1.d0/10.d0*mout(4) +1.d0/40.d0*mout(5) +1.d0/10.d0*mout(6) +1.d0/40.d0*mout(7)           &
+ 1.d0/36.d0*mout(10) + 1.d0/72.d0*mout(11) + 1.d0/12.d0*mout(12) + 1.d0/24.d0*mout(13)    &
+ 1.d0/4.d0*mout(14) + 1.d0/ 8.d0*mout(17) - 1.d0/ 8.d0*mout(18)

      dd_Q19_SW = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
- 1.d0/10.d0*mout(4) -1.d0/40.d0*mout(5) -1.d0/10.d0*mout(6) -1.d0/40.d0*mout(7)           &
+ 1.d0/36.d0*mout(10) + 1.d0/72.d0*mout(11) + 1.d0/12.d0*mout(12) + 1.d0/24.d0*mout(13)    &
+ 1.d0/4.d0*mout(14) - 1.d0/ 8.d0*mout(17) + 1.d0/ 8.d0*mout(18)

      dd_Q19_SE = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
+ 1.d0/10.d0*mout(4) +1.d0/40.d0*mout(5) -1.d0/10.d0*mout(6) -1.d0/40.d0*mout(7)           &
+ 1.d0/36.d0*mout(10) + 1.d0/72.d0*mout(11) + 1.d0/12.d0*mout(12) + 1.d0/24.d0*mout(13)    &
- 1.d0/4.d0*mout(14) + 1.d0/ 8.d0*mout(17)  + 1.d0/ 8.d0*mout(18)

      dd_Q19_NW = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
- 1.d0/10.d0*mout(4) -1.d0/40.d0*mout(5) +1.d0/10.d0*mout(6) +1.d0/40.d0*mout(7)           &
+ 1.d0/36.d0*mout(10) + 1.d0/72.d0*mout(11) + 1.d0/12.d0*mout(12) + 1.d0/24.d0*mout(13)    &
- 1.d0/4.d0*mout(14) - 1.d0/ 8.d0*mout(17) - 1.d0/ 8.d0*mout(18)

      dd_Q19_TE = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
+ 1.d0/10.d0*mout(4) +1.d0/40.d0*mout(5) +1.d0/10.d0*mout(8) +1.d0/40.d0*mout(9)           &
+ 1.d0/36.d0*mout(10) + 1.d0/72.d0*mout(11) - 1.d0/12.d0*mout(12) - 1.d0/24.d0*mout(13)    &
+ 1.d0/4.d0*mout(16) - 1.d0/ 8.d0*mout(17) + 1.d0/ 8.d0*mout(19)

      dd_Q19_BW = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
- 1.d0/10.d0*mout(4) -1.d0/40.d0*mout(5) -1.d0/10.d0*mout(8) -1.d0/40.d0*mout(9)           &
+ 1.d0/36.d0*mout(10) + 1.d0/72.d0*mout(11) - 1.d0/12.d0*mout(12) - 1.d0/24.d0*mout(13)    &
+ 1.d0/4.d0*mout(16) + 1.d0/ 8.d0*mout(17) - 1.d0/ 8.d0*mout(19)

      dd_Q19_BE = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
+ 1.d0/10.d0*mout(4) +1.d0/40.d0*mout(5) -1.d0/10.d0*mout(8) -1.d0/40.d0*mout(9)           &
+ 1.d0/36.d0*mout(10) + 1.d0/72.d0*mout(11) - 1.d0/12.d0*mout(12) - 1.d0/24.d0*mout(13)    &
- 1.d0/4.d0*mout(16) - 1.d0/ 8.d0*mout(17) - 1.d0/ 8.d0*mout(19)

      dd_Q19_TW = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
- 1.d0/10.d0*mout(4) -1.d0/40.d0*mout(5) +1.d0/10.d0*mout(8) +1.d0/40.d0*mout(9)           &
+ 1.d0/36.d0*mout(10) + 1.d0/72.d0*mout(11) - 1.d0/12.d0*mout(12) - 1.d0/24.d0*mout(13)    &
- 1.d0/4.d0*mout(16) + 1.d0/ 8.d0*mout(17) + 1.d0/ 8.d0*mout(19)

      dd_Q19_TN = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
+ 1.d0/10.d0*mout(6) +1.d0/40.d0*mout(7) +1.d0/10.d0*mout(8) +1.d0/40.d0*mout(9)           &
- 1.d0/18.d0*mout(10) - 1.d0/36.d0*mout(11)                                                &
+ 1.d0/4.d0*mout(15) + 1.d0/ 8.d0*mout(18) - 1.d0/ 8.d0*mout(19)

      dd_Q19_BS = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
- 1.d0/10.d0*mout(6) -1.d0/40.d0*mout(7) -1.d0/10.d0*mout(8) -1.d0/40.d0*mout(9)           &
- 1.d0/18.d0*mout(10) - 1.d0/36.d0*mout(11)                                                &
+ 1.d0/4.d0*mout(15) - 1.d0/ 8.d0*mout(18) + 1.d0/ 8.d0*mout(19)

      dd_Q19_BN = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
+ 1.d0/10.d0*mout(6) +1.d0/40.d0*mout(7) -1.d0/10.d0*mout(8) -1.d0/40.d0*mout(9)           &
- 1.d0/18.d0*mout(10) - 1.d0/36.d0*mout(11)                                                &
- 1.d0/4.d0*mout(15) + 1.d0/ 8.d0*mout(18) + 1.d0/ 8.d0*mout(19)

      dd_Q19_TS = one_over_19*mout(1) +  4.d0/1197.d0*mout(2) +1.d0/252.d0*mout(3)        &
- 1.d0/10.d0*mout(6) -1.d0/40.d0*mout(7) +1.d0/10.d0*mout(8) +1.d0/40.d0*mout(9)           &
- 1.d0/18.d0*mout(10) - 1.d0/36.d0*mout(11)                                                &
- 1.d0/4.d0*mout(15) - 1.d0/ 8.d0*mout(18) - 1.d0/ 8.d0*mout(19)


fOut(LB_NODE(1 ,i,j,k)) = dd_Q19_0  
fOut(LB_NODE(2 ,i,j,k)) = dd_Q19_E  
fOut(LB_NODE(3 ,i,j,k)) = dd_Q19_W  
fOut(LB_NODE(4 ,i,j,k)) = dd_Q19_N  
fOut(LB_NODE(5 ,i,j,k)) = dd_Q19_S  
fOut(LB_NODE(6 ,i,j,k)) = dd_Q19_T  
fOut(LB_NODE(7 ,i,j,k)) = dd_Q19_B  
fOut(LB_NODE(8 ,i,j,k)) = dd_Q19_NE 
fOut(LB_NODE(9 ,i,j,k)) = dd_Q19_SW 
fOut(LB_NODE(10,i,j,k)) = dd_Q19_SE 
fOut(LB_NODE(11,i,j,k)) = dd_Q19_NW 
fOut(LB_NODE(12,i,j,k)) = dd_Q19_TE 
fOut(LB_NODE(13,i,j,k)) = dd_Q19_BW 
fOut(LB_NODE(14,i,j,k)) = dd_Q19_BE 
fOut(LB_NODE(15,i,j,k)) = dd_Q19_TW 
fOut(LB_NODE(16,i,j,k)) = dd_Q19_TN 
fOut(LB_NODE(17,i,j,k)) = dd_Q19_BS 
fOut(LB_NODE(18,i,j,k)) = dd_Q19_BN 
fOut(LB_NODE(19,i,j,k)) = dd_Q19_TS 




#ifdef CHECK_MRT
if(i==5.and.j==10.and.k==5) then
   write(*,*) "rho",loc_dens
   write(*,*) "jx ",j_x,j_y,j_z
write(*,*) "meq",meq
write(*,*) "mneq",mneq
write(*,*) "mout",mout
write(*,*) "fIn ",fIn(LB_NODE(:,i,j,k))
write(*,*) "fOut",fOut(LB_NODE(:,i,j,k))
endif
#endif /* CHECK_MRT */
         endif
            end do
         end do
      end do
      
      call cpu_time_measure(meas%tEnd_comm)

      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration


   end subroutine collide_mrt 
  !------------------------------------------------------------------------
#endif /* MRT */ 

#ifdef MRT
#ifndef MRT_SIRANN
   subroutine collide_mrt(lb_dom,fOut,fIn,s_par,meas)

!------------------------------------------------------------------------
! This is the MRT collision step
! Note:
! In the paper of Mei and Luo 2006, drho = rho-rho0 is used instead of rho to reduce roundoff errors.
! This should be implemented here also.
!
! !!! This routine is NOT optimized for performance
! See EXCEL-Sheet BGK_MRT.xls 

      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure)              :: meas  
      real(R8B)                  :: usq,omega
      real(R8B)                  :: meq(nnod),mneq(nnod),mout(nnod),s_mrt(nnod)
      real(R8B)                  :: M_tr(nnod,nnod),M_inv(nnod,nnod),fOut_tmp
      real(R8B)                  :: weps,wepsj,wxx,rholoc_inv 
      real(R8B)                  :: u(3),rholoc,fEq(nnod),t2cs4inv,t2cs2inv,jx(3)
      real(R8B)                  :: rho0 
      real(R8B), intent(out)     :: fOut(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)   :: fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer                    :: i,j,k,l,ll

    real*8, parameter :: one_over__2 = 0.5d0
    real*8, parameter :: one_over__3 = 0.333333333333333333333d0
    real*8, parameter :: two_over__3 = 0.666666666666666666666d0
    real*8, parameter :: one_over__4 = 0.25d0
    real*8, parameter :: one_over__6 = 0.166666666666666666667d0
    real*8, parameter :: one_over__8 = 0.125d0
    real*8, parameter :: one_over__9 = 0.111111111111111111111d0
    real*8, parameter :: one_over_12 = 0.083333333333333333333d0
    real*8, parameter :: one_over_16 = 0.0625d0
    real*8, parameter :: one_over_18 = 0.055555555555555555556d0
    real*8, parameter :: one_over_19 = 0.052631578947368421053d0
    real*8, parameter :: one_over_21 = 0.047619047619047619048d0
    real*8, parameter :: one_over_24 = 0.041666666666666666667d0
    real*8, parameter :: one_over_36 = 0.027777777777777777778d0
    real*8, parameter :: one_over_48 = 0.020833333333333333333d0
    real*8, parameter :: one_over_72 = 0.013888888888888888889d0

      call cpu_time_measure(meas%tSt_comm) 

      omega = s_par%omega
      if(s_par%initial) omega = 1.1

      s_mrt(1:nnod) = 0.d0
#ifdef D2Q9
      s_mrt(2)      = 1.63 !Damps  too much!!
!      s_mrt(2)      = omega 
      s_mrt(3)      = 1.14
      s_mrt(5)      = 1.92
      s_mrt(7)      = 1.92
      s_mrt(8)      = omega
      s_mrt(9)      = omega
      if(s_par%initial) then
         s_mrt(4)      = 1.
         s_mrt(6)      = 1.
      endif
#endif
#ifdef D3Q19
  s_mrt(1) = 0._R8B
!FIXME  s_mrt(2) = 1.19_R8B ! bulk viscosity is damped too much!
  s_mrt(2) = omega   
  s_mrt(3) = 1.4_R8B
  s_mrt(4) = 0._R8B
  s_mrt(5) = 1.2_R8B
  s_mrt(6) = 0._R8B
  s_mrt(7) = 1.2_R8B
  s_mrt(8) = 0._R8B
  s_mrt(9) = 1.2_R8B
  s_mrt(10)= omega
  s_mrt(11)= 1.4_R8B
  s_mrt(12)= omega
  s_mrt(13)= 1.4_R8B
  s_mrt(14)= omega
  s_mrt(15)= omega
  s_mrt(16)= omega
  s_mrt(17)= 1.98_R8B
  s_mrt(18)= 1.98_R8B
  s_mrt(19)= 1.98_R8B

      if(s_par%initial) then
         s_mrt(4) = 1._R8B
         s_mrt(6) = 1._R8B
         s_mrt(8) = 1._R8B
      endif

#ifdef CHECK_MRT
      weps  = 3._R8B 
      wepsj = -5.5_R8B
      wxx   = -0.5_R8B
#else
!FIXME: These parameteres should be taken. But! Strange behaviour of flow.
!      weps  = 0._R8B 
!      wepsj = -475._R8B/63._R8B
!      wxx   = 0._R8B
      weps  = 3._R8B 
      wepsj = -5.5_R8B
      wxx   = -0.5_R8B
#endif /* CHECK_MRT */

#endif

#ifdef CHECK_MRT
     s_mrt(1:nnod) = omega 
#endif /* CHECK_MRT */


#ifdef D2Q9
      M_tr(1:9,1) = (/ 1., 1., 1., 1., 1., 1., 1., 1., 1. /)
      M_tr(1:9,2) = (/-4.,-1.,-1.,-1.,-1., 2., 2., 2., 2. /)
      M_tr(1:9,3) = (/ 4.,-2.,-2.,-2.,-2., 1., 1., 1., 1. /)
      M_tr(1:9,4) = (/ 0., 1., 0.,-1., 0., 1.,-1.,-1., 1. /)
      M_tr(1:9,5) = (/ 0.,-2., 0., 2., 0., 1.,-1.,-1., 1. /)
      M_tr(1:9,6) = (/ 0., 0., 1., 0.,-1., 1., 1.,-1.,-1. /)
      M_tr(1:9,7) = (/ 0., 0.,-2., 0., 2., 1., 1.,-1.,-1. /)
      M_tr(1:9,8) = (/ 0., 1.,-1., 1.,-1., 0., 0., 0., 0. /)
      M_tr(1:9,9) = (/ 0., 0., 0., 0., 0., 1.,-1., 1.,-1. /)

      M_inv(1:9,1) = (/ 1./9., -1./9.,    1./9.,    0.,     0.,    0.,     0.,     0.,    0. /)
      M_inv(1:9,2) = (/ 1./9., -1./36., -1./18., 1./6., -1./6.,    0.,     0.,  1./4.,    0. /)
      M_inv(1:9,3) = (/ 1./9., -1./36., -1./18., 0.   , 0.,     1./6., -1./6., -1./4.,    0. /)
      M_inv(1:9,4) = (/ 1./9., -1./36., -1./18.,-1./6., 1./6.,     0.,     0.,  1./4.,    0. /)
      M_inv(1:9,5) = (/ 1./9., -1./36., -1./18., 0.,    0.,    -1./6.,  1./6., -1./4.,    0. /)
      M_inv(1:9,6) = (/ 1./9.,  1./18.,  1./36., 1./6., 1./12., 1./6.,  1./12.,    0., 1./4. /)
      M_inv(1:9,7) = (/ 1./9.,  1./18.,  1./36.,-1./6.,-1./12., 1./6.,  1./12.,    0.,-1./4. /)
      M_inv(1:9,8) = (/ 1./9.,  1./18.,  1./36.,-1./6.,-1./12.,-1./6., -1./12.,    0., 1./4. /)
      M_inv(1:9,9) = (/ 1./9.,  1./18.,  1./36., 1./6., 1./12.,-1./6., -1./12.,    0.,-1./4. /)
#endif

#ifdef D3Q19

M_tr(1:19, 1) = (/1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0 /)
M_tr(1:19, 2) = (/-30.d0, -11.d0, -11.d0, -11.d0, -11.d0, -11.d0, -11.d0, &
18.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0, 8.d0 /)
M_tr(1:19, 3) = (/12.d0, -4.d0, -4.d0, -4.d0, -4.d0, -4.d0, -4.d0, &
1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0 /)
M_tr(1:19, 4) = (/0.d0, 1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_tr(1:19, 5) = (/0.d0, -4.d0, 4.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0, 1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_tr(1:19, 6) = (/0.d0, 0.d0, 0.d0, 1.d0, -1.d0, 0.d0, 0.d0, &
1.d0, -1.d0, -1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, 1.d0, -1.d0 /)
M_tr(1:19, 7) = (/0.d0, 0.d0, 0.d0, -4.d0, 4.d0, 0.d0, 0.d0, &
1.d0, -1.d0, -1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, 1.d0, -1.d0 /)
M_tr(1:19, 8) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, &
0.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, -1.d0, 1.d0, 1.d0, -1.d0, -1.d0, 1.d0 /)
M_tr(1:19, 9) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -4.d0, 4.d0, &
0.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, -1.d0, 1.d0, 1.d0, -1.d0, -1.d0, 1.d0 /)
M_tr(1:19,10) = (/0.d0, 2.d0, 2.d0, -1.d0, -1.d0, -1.d0, -1.d0, &
1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, -2.d0, -2.d0, -2.d0, -2.d0 /)
M_tr(1:19,11) = (/0.d0, -4.d0, -4.d0, 2.d0, 2.d0, 2.d0, 2.d0, &
1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, -2.d0, -2.d0, -2.d0, -2.d0 /)
M_tr(1:19,12) = (/0.d0, 0.d0, 0.d0, 1.d0, 1.d0, -1.d0, -1.d0, &
1.d0, 1.d0, 1.d0, 1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_tr(1:19,13) = (/0.d0, 0.d0, 0.d0, -2.d0, -2.d0, 2.d0, 2.d0, &
1.d0, 1.d0, 1.d0, 1.d0, -1.d0, -1.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_tr(1:19,14) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
1.d0, 1.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_tr(1:19,15) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 1.d0, -1.d0, -1.d0 /)
M_tr(1:19,16) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 1.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_tr(1:19,17) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
1.d0, -1.d0, 1.d0, -1.d0, -1.d0, 1.d0, -1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_tr(1:19,18) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
-1.d0, 1.d0, 1.d0, -1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, 1.d0, -1.d0 /)
M_tr(1:19,19) = (/0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
0.d0, 0.d0, 0.d0, 0.d0, 1.d0, -1.d0, -1.d0, 1.d0, -1.d0, 1.d0, 1.d0, -1.d0 /)


M_inv(1:19,1)=(/ 1.d0/19.d0, -5.d0/399.d0, 1.d0/21.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_inv(1:19,2)=(/ 1.d0/19.d0, -11.d0/2394.d0, -1.d0/63.d0, 1.d0/10.d0, -1.d0/10.d0, &
0.d0, 0.d0, 0.d0, 0.d0, 1.d0/18.d0, -1.d0/18.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_inv(1:19,3)=(/ 1.d0/19.d0, -11.d0/2394.d0, -1.d0/63.d0, -1.d0/10.d0, 1.d0/10.d0, &
0.d0, 0.d0, 0.d0, 0.d0, 1.d0/18.d0, -1.d0/18.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_inv(1:19,4)=(/ 1.d0/19.d0, -11.d0/2394.d0, -1.d0/63.d0, 0.d0, 0.d0, 1.d0/10.d0, &
-1.d0/10.d0, 0.d0, 0.d0, -1.d0/36.d0, 1.d0/36.d0, 1.d0/12.d0, &
-1.d0/12.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_inv(1:19,5)=(/ 1.d0/19.d0, -11.d0/2394.d0, -1.d0/63.d0, 0.d0, 0.d0, -1.d0/10.d0, &
1.d0/10.d0, 0.d0, 0.d0, -1.d0/36.d0, 1.d0/36.d0, 1.d0/12.d0, &
-1.d0/12.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_inv(1:19,6)=(/ 1.d0/19.d0, -11.d0/2394.d0, -1.d0/63.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
1.d0/10.d0, -1.d0/10.d0, -1.d0/36.d0, 1.d0/36.d0, -1.d0/12.d0, &
1.d0/12.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_inv(1:19,7)=(/ 1.d0/19.d0, -11.d0/2394.d0, -1.d0/63.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
-1.d0/10.d0, 1.d0/10.d0, -1.d0/36.d0, 1.d0/36.d0, -1.d0/12.d0, &
1.d0/12.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)
M_inv(1:19,8)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, 1.d0/10.d0, 1.d0/40.d0, &
1.d0/10.d0, 1.d0/40.d0, 0.d0, 0.d0, 1.d0/36.d0, 1.d0/72.d0, 1.d0/12.d0, &
1.d0/24.d0, 1.d0/4.d0, 0.d0, 0.d0, 1.d0/8.d0, -1.d0/8.d0, 0.d0 /)
M_inv(1:19,9)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, -1.d0/10.d0, &
-1.d0/40.d0, -1.d0/10.d0, -1.d0/40.d0, 0.d0, 0.d0, 1.d0/36.d0, 1.d0/72.d0, &
1.d0/12.d0, 1.d0/24.d0, 1.d0/4.d0, 0.d0, 0.d0, -1.d0/8.d0, 1.d0/8.d0, 0.d0 /)
M_inv(1:19,10)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, 1.d0/10.d0, 1.d0/40.d0, &
-1.d0/10.d0, -1.d0/40.d0, 0.d0, 0.d0, 1.d0/36.d0, 1.d0/72.d0, 1.d0/12.d0, &
1.d0/24.d0, -1.d0/4.d0, 0.d0, 0.d0, 1.d0/8.d0, 1.d0/8.d0, 0.d0 /)
M_inv(1:19,11)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, -1.d0/10.d0, -1.d0/40.d0, &
1.d0/10.d0, 1.d0/40.d0, 0.d0, 0.d0, 1.d0/36.d0, 1.d0/72.d0, 1.d0/12.d0, 1.d0/24.d0,&
 -1.d0/4.d0, 0.d0, 0.d0, -1.d0/8.d0, -1.d0/8.d0, 0.d0 /)
M_inv(1:19,12)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, 1.d0/10.d0, 1.d0/40.d0, &
0.d0, 0.d0, 1.d0/10.d0, 1.d0/40.d0, 1.d0/36.d0, 1.d0/72.d0, -1.d0/12.d0, &
-1.d0/24.d0, 0.d0, 0.d0, 1.d0/4.d0, -1.d0/8.d0, 0.d0, 1.d0/8.d0 /)
M_inv(1:19,13)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, -1.d0/10.d0, -1.d0/40.d0, &
0.d0, 0.d0, -1.d0/10.d0, -1.d0/40.d0, 1.d0/36.d0, 1.d0/72.d0, -1.d0/12.d0, &
-1.d0/24.d0, 0.d0, 0.d0, 1.d0/4.d0, 1.d0/8.d0, 0.d0, -1.d0/8.d0 /)
M_inv(1:19,14)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, 1.d0/10.d0, 1.d0/40.d0, &
0.d0, 0.d0, -1.d0/10.d0, -1.d0/40.d0, 1.d0/36.d0, 1.d0/72.d0, -1.d0/12.d0, &
-1.d0/24.d0, 0.d0, 0.d0, -1.d0/4.d0, -1.d0/8.d0, 0.d0, -1.d0/8.d0 /)
M_inv(1:19,15)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, -1.d0/10.d0, -1.d0/40.d0, &
0.d0, 0.d0, 1.d0/10.d0, 1.d0/40.d0, 1.d0/36.d0, 1.d0/72.d0, -1.d0/12.d0, -1.d0/24.d0, &
0.d0, 0.d0, -1.d0/4.d0, 1.d0/8.d0, 0.d0, 1.d0/8.d0 /)
M_inv(1:19,16)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, 0.d0, 0.d0, 1.d0/10.d0, &
1.d0/40.d0, 1.d0/10.d0, 1.d0/40.d0, -1.d0/18.d0, -1.d0/36.d0, 0.d0, 0.d0, 0.d0, &
1.d0/4.d0, 0.d0, 0.d0, 1.d0/8.d0, -1.d0/8.d0 /)
M_inv(1:19,17)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, 0.d0, 0.d0, -1.d0/10.d0, &
-1.d0/40.d0, -1.d0/10.d0, -1.d0/40.d0, -1.d0/18.d0, -1.d0/36.d0, 0.d0, 0.d0, 0.d0, &
1.d0/4.d0, 0.d0, 0.d0, -1.d0/8.d0, 1.d0/8.d0 /)
M_inv(1:19,18)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, 0.d0, 0.d0, 1.d0/10.d0, &
1.d0/40.d0, -1.d0/10.d0, -1.d0/40.d0, -1.d0/18.d0, -1.d0/36.d0, 0.d0, 0.d0, 0.d0, &
-1.d0/4.d0, 0.d0, 0.d0, 1.d0/8.d0, 1.d0/8.d0 /)
M_inv(1:19,19)=(/ 1.d0/19.d0, 4.d0/1197.d0, 1.d0/252.d0, 0.d0, 0.d0, -1.d0/10.d0, &
-1.d0/40.d0, 1.d0/10.d0, 1.d0/40.d0, -1.d0/18.d0, -1.d0/36.d0, 0.d0, 0.d0, 0.d0, &
-1.d0/4.d0, 0.d0, 0.d0, -1.d0/8.d0, -1.d0/8.d0 /)
 


#endif

 

      ! Multiple Single Relaxation time Collision routine
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
! How shall rho be used? rho' and rho0 or simply rholoc?
#ifdef D2Q9
              rholoc =   fIn(LB_NODE(1,i,j,k)) + fIn(LB_NODE(2,i,j,k)) + & 
                         fIn(LB_NODE(3,i,j,k) )+ fIn(LB_NODE(4,i,j,k)) + &
                         fIn(LB_NODE(5,i,j,k)) + fIn(LB_NODE(6,i,j,k)) + & 
                         fIn(LB_NODE(7,i,j,k)) + fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(9,i,j,k) )

! Am I currently relaxing with fixed velocity field?

              if(s_par%initial .eqv. .false.) then
! - No: Calculate macroscopic Variables (momentum and density) from pdfs
              jx(1) =  (fIn(LB_NODE(2,i,j,k))  - fIn(LB_NODE(4,i,j,k)) + &
                         fIn(LB_NODE(6,i,j,k)) - fIn(LB_NODE(7,i,j,k)) - & 
                         fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(9,i,j,k)))
              jx(2) =  (fIn(LB_NODE(3,i,j,k))  - fIn(LB_NODE(5,i,j,k)) + &
                        fIn(LB_NODE(6,i,j,k)) + fIn(LB_NODE(7,i,j,k)) - & 
                        fIn(LB_NODE(8,i,j,k)) - fIn(LB_NODE(9,i,j,k)))
               else
! - Yes: Read Velocity field from Initial Field
                 jx(1) = lb_dom%u0(LB_NODE(1,i,j,k))*rholoc
                 jx(2) = lb_dom%u0(LB_NODE(2,i,j,k))*rholoc
               endif
! Incompressible assumption

!            rholoc = rholoc - s_par%rho0

! rho
                  meq(1) =  rholoc
! e
                  meq(2) =  -2.*rholoc + 3.*(jx(1)*jx(1) + jx(2)*jx(2))
! eps
                  meq(3) =  rholoc - 3.*(jx(1)*jx(1) + jx(2)*jx(2))
                  meq(4) =  jx(1)
                  meq(5) = -jx(1)
                  meq(6) =  jx(2)
                  meq(7) = -jx(2)
                  meq(8) =  jx(1)*jx(1) - jx(2)*jx(2)
                  meq(9) =  jx(1)*jx(2)
#endif
#ifdef D3Q19
               rholoc  = fIn(LB_NODE(1,i,j,k))  + fIn(LB_NODE(2,i,j,k))  + &
                         fIn(LB_NODE(3,i,j,k))  + fIn(LB_NODE(4,i,j,k))  + fIn(LB_NODE(5,i,j,k)) + &
                         fIn(LB_NODE(6,i,j,k))  + fIn(LB_NODE(7,i,j,k))  + &
                         fIn(LB_NODE(8,i,j,k))  + fIn(LB_NODE(9,i,j,k))  + fIn(LB_NODE(10,i,j,k)) + & 
                         fIn(LB_NODE(11,i,j,k)) + fIn(LB_NODE(12,i,j,k)) + &
                         fIn(LB_NODE(13,i,j,k)) + fIn(LB_NODE(14,i,j,k)) + fIn(LB_NODE(15,i,j,k)) + & 
                         fIn(LB_NODE(16,i,j,k)) + fIn(LB_NODE(17,i,j,k)) + &
                         fIn(LB_NODE(18,i,j,k)) + fIn(LB_NODE(19,i,j,k))
              if(s_par%initial .eqv. .false.) then
              jx(1) = (  fIn(LB_NODE(2,i,j,k))  - fIn(LB_NODE(3,i,j,k))  + &
                         fIn(LB_NODE(8,i,j,k))  - fIn(LB_NODE(9,i,j,k))  + &
                         fIn(LB_NODE(10,i,j,k)) - fIn(LB_NODE(11,i,j,k)) + & 
                         fIn(LB_NODE(12,i,j,k)) - fIn(LB_NODE(13,i,j,k)) + &
                         fIn(LB_NODE(14,i,j,k)) - fIn(LB_NODE(15,i,j,k))  & 
                             )           

              jx(2) = (  fIn(LB_NODE(4,i,j,k))  - fIn(LB_NODE(5,i,j,k)) + &
                         fIn(LB_NODE(8,i,j,k))  - fIn(LB_NODE(9,i,j,k))  - fIn(LB_NODE(10,i,j,k)) + & 
                         fIn(LB_NODE(11,i,j,k)) + & 
                         fIn(LB_NODE(16,i,j,k)) - fIn(LB_NODE(17,i,j,k)) + &
                         fIn(LB_NODE(18,i,j,k)) - fIn(LB_NODE(19,i,j,k)))

              jx(3) = (  fIn(LB_NODE(6,i,j,k))  - fIn(LB_NODE(7,i,j,k))  + &
                         fIn(LB_NODE(12,i,j,k)) - fIn(LB_NODE(13,i,j,k))- fIn(LB_NODE(14,i,j,k))  + & 
                         fIn(LB_NODE(15,i,j,k)) + fIn(LB_NODE(16,i,j,k)) - &
                         fIn(LB_NODE(17,i,j,k)) - fIn(LB_NODE(18,i,j,k)) + &
                         fIn(LB_NODE(19,i,j,k)))


               else
                 jx(1) = lb_dom%u0(LB_NODE(1,i,j,k))*rholoc
                 jx(2) = lb_dom%u0(LB_NODE(2,i,j,k))*rholoc
                 jx(3) = lb_dom%u0(LB_NODE(3,i,j,k))*rholoc
               endif

               rholoc_inv = 1.d0/rholoc
   
                  meq( 1) = rholoc 
! e_eq
                  meq( 2) =  -11.*rholoc + 19.*rholoc_inv*((jx(1))*(jx(1)) + &
                          (jx(2))*(jx(2))+(jx(3))*(jx(3)))
! eps_eq
                  meq( 3) =  weps*rholoc+wepsj*rholoc_inv*((jx(1))*(jx(1)) + &
                          (jx(2))*(jx(2))+(jx(3))*(jx(3)))
                  meq( 4) =  jx(1)
                  meq( 6) =  jx(2)
                  meq( 8) =  jx(3)
! qx_eq, qy_eq, qz_eq
                  meq( 5) =  -2./3.*jx(1)
                  meq( 7) =  -2./3.*jx(2)
                  meq( 9) =  -2./3.*jx(3)
! p_xx, pi_xx
                  meq(10) = 1./3.*rholoc_inv*(2.*jx(1)*jx(1) - (jx(2)*jx(2) + jx(3)*jx(3)))
                  meq(11) = wxx*meq(10) 
! p_ww, pi_ww
                  meq(12) = rholoc_inv*(jx(1)*jx(1) - jx(3)*jx(3))
                  meq(13) = wxx*meq(12) 
!p_xy,p_yz,p_xz
                  meq(14) =  rholoc_inv*jx(1)*jx(2)
                  meq(15) =  rholoc_inv*jx(2)*jx(3)
                  meq(16) =  rholoc_inv*jx(1)*jx(3)
                  meq(17) =  0.0d0 
                  meq(18) =  0.0d0
                  meq(19) =  0.0d0
#endif
               do l=1,nnod
                  mneq(l) = 0.0d0
                  ! Berechnung der  nicht-gg-Werte im Momentenraum
                  do ll=1,nnod
                     mneq(l) = mneq(l) + fIn(LB_NODE(ll,i,j,k))*M_tr(ll,l)
                  end do
               end do
               ! relaxation
               do l=1,nnod
                  mout(l) = mneq(l) - s_mrt(l)*(mneq(l) - meq(l)) 
               end do 
               ! Ruecktransformation
               do l=1,nnod
                  fOut_tmp = 0.0d0
                  ! Berechnung der  nicht-gg-Werte im Momentenraum
                  do ll=1,nnod
                     fOut_tmp = fOut_tmp + mout(ll)*M_inv(ll,l)
                  end do
                  fOut(LB_NODE(l,i,j,k))= fOut_tmp
               end do
#ifdef CHECK_MRT
if(i==5.and.j==10.and.k==5) then
   write(*,*) "rho",rholoc
   write(*,*) "jx ",jx
write(*,*) "meq",meq
write(*,*) "mneq",mneq
write(*,*) "mout",mout
write(*,*) "fIn ",fIn(LB_NODE(:,i,j,k))
write(*,*) "fOut",fOut(LB_NODE(:,i,j,k))
endif
#endif /* CHECK_MRT */
         endif
            end do
         end do
      end do
      
      call cpu_time_measure(meas%tEnd_comm)

      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration


   end subroutine collide_mrt
  !------------------------------------------------------------------------



#else /* SIRANN MRT */


   subroutine collide_mrt(lb_dom,fOut,fIn,s_par,meas)

!------------------------------------------------------------------------
! This is the MRT collision step
! Note:
! In the paper of Mei and Luo 2006, drho = rho-rho0 is used instead of rho to reduce roundoff errors.
! This should be implemented here also.
!
! !!! This routine is NOT optimized for performance
! See EXCEL-Sheet BGK_MRT.xls 

      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure)              :: meas  
      real(R8B)                  :: usq,omega
      real(R8B)                  :: meq(nnod),mneq(nnod),mout(nnod),s_mrt(nnod)
      real(R8B)                  :: M_tr(nnod,nnod),M_inv(nnod,nnod),fOut_tmp
      real(R8B)                  :: weps,wepsj,wxx,rholoc_inv 
      real(R8B)                  :: ux(3),rholoc,fEq(nnod),t2cs4inv,t2cs2inv,jx(3)
      real(R8B)                  :: rho0 
      real(R8B), intent(out)     :: fOut(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)   :: fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer                    :: i,j,k,l,ll


      call cpu_time_measure(meas%tSt_comm) 

      omega = s_par%omega
      if(s_par%initial) omega = 1.1

      s_mrt(1:nnod) = 0.d0
#ifdef D2Q9
!      s_mrt(2)      = 1.63 Damped too much!!
      s_mrt(2)      = omega 
      s_mrt(3)      = 1.14
      s_mrt(5)      = 1.92
      s_mrt(7)      = 1.92
      s_mrt(8)      = omega
      s_mrt(9)      = omega
      if(s_par%initial) then
         s_mrt(4)      = 1.
         s_mrt(6)      = 1.
      endif
#endif
#ifdef D3Q19
  s_mrt(1) = 0._R8B
!FIXME  s_mrt(2) = 1.19_R8B ! bulk viscosity is damped too much!
  s_mrt(2) = omega   
  s_mrt(3) = 1.4_R8B
  s_mrt(4) = 0._R8B
  s_mrt(5) = 1.2_R8B
  s_mrt(6) = 0._R8B
  s_mrt(7) = 1.2_R8B
  s_mrt(8) = 0._R8B
  s_mrt(9) = 1.2_R8B
  s_mrt(10)= omega
  s_mrt(11)= 1.4_R8B
  s_mrt(12)= omega
  s_mrt(13)= 1.4_R8B
  s_mrt(14)= omega
  s_mrt(15)= omega
  s_mrt(16)= omega
  s_mrt(17)= 1.98_R8B
  s_mrt(18)= 1.98_R8B
  s_mrt(19)= 1.98_R8B

      if(s_par%initial) then
         s_mrt(4) = 1._R8B
         s_mrt(6) = 1._R8B
         s_mrt(8) = 1._R8B
      endif

#ifdef CHECK_MRT
      weps  = 3._R8B 
      wepsj = -5.5_R8B
      wxx   = -0.5_R8B
#else
!FIXME: These parameteres should be taken. But! Strange behaviour of flow.
!      weps  = 0._R8B 
!      wepsj = -475._R8B/63._R8B
!      wxx   = 0._R8B
      weps  = 3._R8B 
      wepsj = -5.5_R8B
      wxx   = -0.5_R8B
#endif /* CHECK_MRT */

#endif

#ifdef CHECK_MRT
     s_mrt(1:nnod) = omega 
#endif /* CHECK_MRT */


#ifdef D2Q9
      M_tr(1:9,1) = (/ 1., 1., 1., 1., 1., 1., 1., 1., 1. /)
      M_tr(1:9,2) = (/-4.,-1.,-1.,-1.,-1., 2., 2., 2., 2. /)
      M_tr(1:9,3) = (/ 4.,-2.,-2.,-2.,-2., 1., 1., 1., 1. /)
      M_tr(1:9,4) = (/ 0., 1., 0.,-1., 0., 1.,-1.,-1., 1. /)
      M_tr(1:9,5) = (/ 0.,-2., 0., 2., 0., 1.,-1.,-1., 1. /)
      M_tr(1:9,6) = (/ 0., 0., 1., 0.,-1., 1., 1.,-1.,-1. /)
      M_tr(1:9,7) = (/ 0., 0.,-2., 0., 2., 1., 1.,-1.,-1. /)
      M_tr(1:9,8) = (/ 0., 1.,-1., 1.,-1., 0., 0., 0., 0. /)
      M_tr(1:9,9) = (/ 0., 0., 0., 0., 0., 1.,-1., 1.,-1. /)

      M_inv(1:9,1) = (/ 1./9., -1./9.,    1./9.,    0.,     0.,    0.,     0.,     0.,    0. /)
      M_inv(1:9,2) = (/ 1./9., -1./36., -1./18., 1./6., -1./6.,    0.,     0.,  1./4.,    0. /)
      M_inv(1:9,3) = (/ 1./9., -1./36., -1./18., 0.   , 0.,     1./6., -1./6., -1./4.,    0. /)
      M_inv(1:9,4) = (/ 1./9., -1./36., -1./18.,-1./6., 1./6.,     0.,     0.,  1./4.,    0. /)
      M_inv(1:9,5) = (/ 1./9., -1./36., -1./18., 0.,    0.,    -1./6.,  1./6., -1./4.,    0. /)
      M_inv(1:9,6) = (/ 1./9.,  1./18.,  1./36., 1./6., 1./12., 1./6.,  1./12.,    0., 1./4. /)
      M_inv(1:9,7) = (/ 1./9.,  1./18.,  1./36.,-1./6.,-1./12., 1./6.,  1./12.,    0.,-1./4. /)
      M_inv(1:9,8) = (/ 1./9.,  1./18.,  1./36.,-1./6.,-1./12.,-1./6., -1./12.,    0., 1./4. /)
      M_inv(1:9,9) = (/ 1./9.,  1./18.,  1./36., 1./6., 1./12.,-1./6., -1./12.,    0.,-1./4. /)
#endif

#ifdef D3Q19

M_tr(1:19, 1) = (/1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1. /)
M_tr(1:19, 2) = (/-1.,0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./)
M_tr(1:19, 3) = (/1.,-2.,-2.,-2.,-2.,-2.,-2.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./)
M_tr(1:19, 4) = (/0.,1.,-1.,0.,0.,0.,0.,1.,-1.,1.,-1.,1.,-1.,1.,-1.,0.,0.,0.,0./)
M_tr(1:19, 5) = (/0.,-2.,2.,0.,0.,0.,0.,1.,-1.,1.,-1.,1.,-1.,1.,-1.,0.,0.,0.,0./)
M_tr(1:19, 6) = (/0.,0.,0.,1.,-1.,0.,0.,1.,-1.,-1.,1.,0.,0.,0.,0.,1.,-1.,1.,-1./)
M_tr(1:19, 7) = (/0.,0.,0.,-2.,2.,0.,0.,1.,-1.,-1.,1.,0.,0.,0.,0.,1.,-1.,1.,-1./)
M_tr(1:19, 8) = (/0.,0.,0.,0.,0.,1.,-1.,0.,0.,0.,0.,1.,-1.,-1.,1.,1.,-1.,-1.,1./)
M_tr(1:19, 9) = (/0.,0.,0.,0.,0.,-2.,2.,0.,0.,0.,0.,1.,-1.,-1.,1.,1.,-1.,-1.,1./)
M_tr(1:19,10) = (/0.,2.,2.,-1.,-1.,-1.,-1.,1.,1.,1.,1.,1.,1.,1.,1.,-2.,-2.,-2.,-2./)
M_tr(1:19,11) = (/0.,-2.,-2.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,-2.,-2.,-2.,-2./)
M_tr(1:19,12) = (/0.,0.,0.,1.,1.,-1.,-1.,1.,1.,1.,1.,-1.,-1.,-1.,-1.,0.,0.,0.,0./)
M_tr(1:19,13) = (/0.,0.,0.,-1.,-1.,1.,1.,1.,1.,1.,1.,-1.,-1.,-1.,-1.,0.,0.,0.,0./)
M_tr(1:19,14) = (/0.,0.,0.,0.,0.,0.,0.,1.,1.,-1.,-1.,0.,0.,0.,0.,0.,0.,0.,0./)
M_tr(1:19,15) = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,-1.,-1./)
M_tr(1:19,16) = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,1.,-1.,-1.,0.,0.,0.,0./)
M_tr(1:19,17) = (/0.,0.,0.,0.,0.,0.,0.,1.,-1.,1.,-1.,-1.,1.,-1.,1.,0.,0.,0.,0./)
M_tr(1:19,18) = (/0.,0.,0.,0.,0.,0.,0.,-1.,1.,1.,-1.,0.,0.,0.,0.,1.,-1.,1.,-1./)
M_tr(1:19,19) = (/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1.,-1.,-1.,1.,-1.,1.,1.,-1.  /) 


M_inv(1:19,1) =(/1./3.,  -1./2.,   1./7.,     0.,     0.,     0.,     0.,     0.,  &
    0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0.    /)  
M_inv(1:19,2) =(/1./18.,     0., -1./18.,   1./6.,  -1./6.,     0.,     0.,     0.,&
     0.,  1./12., -1./12.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0. /) 
M_inv(1:19,3) =(/1./18.,     0., -1./18.,  -1./6.,   1./6.,     0.,     0.,     0.,&
     0.,  1./12., -1./12.,     0.,     0.,     0.,     0.,     0.,     0.,     0.,     0. /) 
M_inv(1:19,4) =(/1./18.,     0., -1./18.,     0.,     0.,   1./6.,  -1./6.,     0.,&
     0., -1./24.,  1./24.,   1./8.,  -1./8.,     0.,     0.,     0.,     0.,     0.,     0. /) 
M_inv(1:19,5) =(/1./18.,     0., -1./18.,     0.,     0.,  -1./6.,   1./6.,     0.,&
     0., -1./24.,  1./24.,   1./8.,  -1./8.,     0.,     0.,     0.,     0.,     0.,     0. /) 
M_inv(1:19,6) =(/1./18.,     0., -1./18.,     0.,     0.,     0.,     0.,   1./6.,&
  -1./6., -1./24.,  1./24.,  -1./8.,   1./8.,     0.,     0.,     0.,     0.,     0.,     0. /) 
M_inv(1:19,7) =(/1./18.,     0., -1./18.,     0.,     0.,     0.,     0.,  -1./6.,&
   1./6., -1./24.,  1./24.,  -1./8.,   1./8.,     0.,     0.,     0.,     0.,     0.,     0. /) 
M_inv(1:19,8) =(/1./36.,  1./24.,  1./72.,  1./12.,  1./24.,  1./12.,  1./24.,   &
  0.,     0.,  1./48.,  1./48.,  1./16.,  1./16.,   1./4.,     0.,     0.,   1./8.,  -1./8.,0. /) 
M_inv(1:19,9) =(/1./36.,  1./24.,  1./72., -1./12., -1./24., -1./12., -1./24.,   &
  0.,     0.,  1./48.,  1./48.,  1./16.,  1./16.,   1./4.,     0.,     0.,  -1./8.,   1./8.,0. /) 
M_inv(1:19,10)=(/1./36.,  1./24.,  1./72.,  1./12.,  1./24., -1./12., -1./24.,   &
  0.,     0.,  1./48.,  1./48.,  1./16.,  1./16.,  -1./4.,     0.,     0.,   1./8.,   1./8.,0. /) 
M_inv(1:19,11)=(/1./36.,  1./24.,  1./72., -1./12., -1./24.,  1./12.,  1./24.,   &
  0.,     0.,  1./48.,  1./48.,  1./16.,  1./16.,  -1./4.,     0.,     0.,  -1./8.,  -1./8.,0. /) 
M_inv(1:19,12)=(/1./36.,  1./24.,  1./72.,  1./12.,  1./24.,     0.,     0.,  1./12.,&
  1./24.,  1./48.,  1./48., -1./16., -1./16.,     0.,     0.,   1./4.,  -1./8.,     0.,   1./8./) 
M_inv(1:19,13)=(/1./36.,  1./24.,  1./72., -1./12., -1./24.,     0.,     0., -1./12.,&
 -1./24.,  1./48.,  1./48., -1./16., -1./16.,     0.,     0.,   1./4.,   1./8.,     0.,  -1./8./) 
M_inv(1:19,14)=(/1./36.,  1./24.,  1./72.,  1./12.,  1./24.,     0.,     0., -1./12.,&
 -1./24.,  1./48.,  1./48., -1./16., -1./16.,     0.,     0.,  -1./4.,  -1./8.,     0.,  -1./8./) 
M_inv(1:19,15)=(/1./36.,  1./24.,  1./72., -1./12., -1./24.,     0.,     0.,  1./12.,&
  1./24.,  1./48.,  1./48., -1./16., -1./16.,     0.,     0.,  -1./4.,   1./8.,     0.,   1./8./) 
M_inv(1:19,16)=(/1./36.,  1./24.,  1./72.,     0.,     0.,  1./12.,  1./24.,  1./12.,&
  1./24., -1./24., -1./24.,     0.,     0.,     0.,   1./4.,     0.,     0.,   1./8.,  -1./8./) 
M_inv(1:19,17)=(/1./36.,  1./24.,  1./72.,     0.,     0., -1./12., -1./24., -1./12.,&
 -1./24., -1./24., -1./24.,     0.,     0.,     0.,   1./4.,     0.,     0.,  -1./8.,   1./8./) 
M_inv(1:19,18)=(/1./36.,  1./24.,  1./72.,     0.,     0.,  1./12.,  1./24., -1./12.,&
 -1./24., -1./24., -1./24.,     0.,     0.,     0.,  -1./4.,     0.,     0.,   1./8.,   1./8./) 
M_inv(1:19,19)=(/1./36.,  1./24.,  1./72.,     0.,     0., -1./12., -1./24.,  1./12.,&
  1./24., -1./24., -1./24.,     0.,     0.,     0.,  -1./4.,     0.,     0.,  -1./8.,  -1./8.     /) 
 


#endif

 

      ! Multiple Single Relaxation time Collision routine
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
! How shall rho be used? rho' and rho0 or simply rholoc?
#ifdef D2Q9
              rholoc =   fIn(LB_NODE(1,i,j,k)) + fIn(LB_NODE(2,i,j,k)) + & 
                         fIn(LB_NODE(3,i,j,k) )+ fIn(LB_NODE(4,i,j,k)) + &
                         fIn(LB_NODE(5,i,j,k)) + fIn(LB_NODE(6,i,j,k)) + & 
                         fIn(LB_NODE(7,i,j,k)) + fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(9,i,j,k) )

! Am I currently relaxing with fixed velocity field?

              if(s_par%initial .eqv. .false.) then
! - No: Calculate macroscopic Variables (momentum and density) from pdfs
              jx(1) =  (fIn(LB_NODE(2,i,j,k))  - fIn(LB_NODE(4,i,j,k)) + &
                         fIn(LB_NODE(6,i,j,k)) - fIn(LB_NODE(7,i,j,k)) - & 
                         fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(9,i,j,k)))
              jx(2) =  (fIn(LB_NODE(3,i,j,k))  - fIn(LB_NODE(5,i,j,k)) + &
                        fIn(LB_NODE(6,i,j,k)) + fIn(LB_NODE(7,i,j,k)) - & 
                        fIn(LB_NODE(8,i,j,k)) - fIn(LB_NODE(9,i,j,k)))
               else
! - Yes: Read Velocity field from Initial Field
                 jx(1) = lb_dom%u0(LB_NODE(1,i,j,k))*rholoc
                 jx(2) = lb_dom%u0(LB_NODE(2,i,j,k))*rholoc
               endif
! Incompressible assumption

!            rholoc = rholoc - s_par%rho0
rho0=1.0
! rho
                  meq(1) =  rholoc
! e
                  meq(2) =  (jx(1)*jx(1) + jx(2)*jx(2))
! eps
                  meq(3) =  rholoc - 3.*(jx(1)*jx(1) + jx(2)*jx(2))
                  meq(4) =  jx(1)
                  meq(5) = -jx(1)
                  meq(6) =  jx(2)
                  meq(7) = -jx(2)
                  meq(8) =  jx(1)*jx(1) - jx(2)*jx(2)
                  meq(9) =  jx(1)*jx(2)
#endif
#ifdef D3Q19
               rholoc  = fIn(LB_NODE(1,i,j,k))  + fIn(LB_NODE(2,i,j,k))  + &
                         fIn(LB_NODE(3,i,j,k))  + fIn(LB_NODE(4,i,j,k))  + fIn(LB_NODE(5,i,j,k)) + &
                         fIn(LB_NODE(6,i,j,k))  + fIn(LB_NODE(7,i,j,k))  + &
                         fIn(LB_NODE(8,i,j,k))  + fIn(LB_NODE(9,i,j,k))  + fIn(LB_NODE(10,i,j,k)) + & 
                         fIn(LB_NODE(11,i,j,k)) + fIn(LB_NODE(12,i,j,k)) + &
                         fIn(LB_NODE(13,i,j,k)) + fIn(LB_NODE(14,i,j,k)) + fIn(LB_NODE(15,i,j,k)) + & 
                         fIn(LB_NODE(16,i,j,k)) + fIn(LB_NODE(17,i,j,k)) + &
                         fIn(LB_NODE(18,i,j,k)) + fIn(LB_NODE(19,i,j,k))
              if(s_par%initial .eqv. .false.) then
              jx(1) = (  fIn(LB_NODE(2,i,j,k))  - fIn(LB_NODE(3,i,j,k))  + &
                         fIn(LB_NODE(8,i,j,k))  - fIn(LB_NODE(9,i,j,k))  + &
                         fIn(LB_NODE(10,i,j,k)) - fIn(LB_NODE(11,i,j,k)) + & 
                         fIn(LB_NODE(12,i,j,k)) - fIn(LB_NODE(13,i,j,k)) + &
                         fIn(LB_NODE(14,i,j,k)) - fIn(LB_NODE(15,i,j,k))  & 
                             )           

              jx(2) = (  fIn(LB_NODE(4,i,j,k))  - fIn(LB_NODE(5,i,j,k)) + &
                         fIn(LB_NODE(8,i,j,k))  - fIn(LB_NODE(9,i,j,k))  - fIn(LB_NODE(10,i,j,k)) + & 
                         fIn(LB_NODE(11,i,j,k)) + & 
                         fIn(LB_NODE(16,i,j,k)) - fIn(LB_NODE(17,i,j,k)) + &
                         fIn(LB_NODE(18,i,j,k)) - fIn(LB_NODE(19,i,j,k)))

              jx(3) = (  fIn(LB_NODE(6,i,j,k))  - fIn(LB_NODE(7,i,j,k))  + &
                         fIn(LB_NODE(12,i,j,k)) - fIn(LB_NODE(13,i,j,k))- fIn(LB_NODE(14,i,j,k))  + & 
                         fIn(LB_NODE(15,i,j,k)) + fIn(LB_NODE(16,i,j,k)) - &
                         fIn(LB_NODE(17,i,j,k)) - fIn(LB_NODE(18,i,j,k)) + &
                         fIn(LB_NODE(19,i,j,k)))


               else
                 jx(1) = lb_dom%u0(LB_NODE(1,i,j,k))*rholoc
                 jx(2) = lb_dom%u0(LB_NODE(2,i,j,k))*rholoc
                 jx(3) = lb_dom%u0(LB_NODE(3,i,j,k))*rholoc
               endif

               rholoc_inv = 1.d0/rholoc

                 ux(1) = jx(1)*rholoc_inv
                 ux(2) = jx(2)*rholoc_inv
                 ux(3) = jx(3)*rholoc_inv

                 rho0=1.0 
                  meq =  0.0d0 

                  meq( 1) = rholoc 
! e_eq
                  meq( 2) =   rho0*((ux(1))*(ux(1)) + &
                          (ux(2))*(ux(2))+(ux(3))*(ux(3)))
! eps_eq
                  meq( 4) =  jx(1)
                  meq( 6) =  jx(2)
                  meq( 8) =  jx(3)
! qx_eq, qy_eq, qz_eq
                  meq( 3) = ux(1)*rho0 
                  meq( 5) = ux(2)*rho0 
                  meq( 7) = ux(3)*rho0 
                  meq( 9) = rho0*(2.*ux(1)*ux(1) - (ux(2)*ux(2) + ux(3)*ux(3)))
! p_xx, pi_xx
                  meq(11) = rho0*(ux(2)*ux(2) - ux(3)*ux(3))
!p_xy,p_yz,p_xz
                  meq(13) =  rho0*ux(1)*ux(2)
                  meq(14) =  rho0*ux(2)*ux(3)
                  meq(15) =  rho0*ux(1)*ux(3)
#endif
               do l=1,nnod
                  mneq(l) = 0.0d0
                  ! Berechnung der  nicht-gg-Werte im Momentenraum
                  do ll=1,nnod
                     mneq(l) = mneq(l) + fIn(LB_NODE(ll,i,j,k))*M_tr(ll,l)
                  end do
               end do
               ! relaxation
               do l=1,nnod
                  mout(l) = mneq(l) - s_mrt(l)*(mneq(l) - meq(l)) 
               end do 
               ! Rücktransformation
               do l=1,nnod
                  fOut_tmp = 0.0d0
                  ! Berechnung der  nicht-gg-Werte im Momentenraum
                  do ll=1,nnod
                     fOut_tmp = fOut_tmp + mout(ll)*M_inv(ll,l)
                  end do
                  fOut(LB_NODE(l,i,j,k))= fOut_tmp
               end do
#ifdef CHECK_MRT
if(i==5.and.j==10.and.k==5) then
   write(*,*) "rho",rholoc
   write(*,*) "jx ",jx
write(*,*) "meq",meq
write(*,*) "mneq",mneq
write(*,*) "mout",mout
write(*,*) "fIn ",fIn(LB_NODE(:,i,j,k))
write(*,*) "fOut",fOut(LB_NODE(:,i,j,k))
endif
#endif /* CHECK_MRT */
         endif
            end do
         end do
      end do
      
      call cpu_time_measure(meas%tEnd_comm)

      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration


   end subroutine collide_mrt
  !------------------------------------------------------------------------


#endif /* MRT_SIRANN*/
#endif /* MRT */








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
      real(R8B), intent(out)     :: fOut(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)   :: fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer                    :: i,j,k,l


      call cpu_time_measure(meas%tSt_comm) 
!      om0 = 1.5
      omega = s_par%omega
      t2cs4inv = 1._R8B/(2._R8B*cs*cs*cs*cs)
      t2cs2inv = 1._R8B/(2._R8B*cs*cs)
      rho0=s_par%rho0
      if(s_par%initial) then
!         ramp = real(gtstep_cur)/real(s_par%tInitial)
!      om_ch = om0 !+ (s_par%omega-om0)/(real(s_par%tInitial))*real(gtstep_cur)
!      write(*,*) "current omega",om_ch
      endif
      ! Single Relaxation time Collision routine
      do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
#ifdef SPONGE
            omega = lb_dom%omega(i,j,k)
#endif
#ifdef D3Q19
               rholoc  = fIn(LB_NODE(1,i,j,k))  + fIn(LB_NODE(2,i,j,k))  + &
                         fIn(LB_NODE(3,i,j,k))  + fIn(LB_NODE(4,i,j,k))  + fIn(LB_NODE(5,i,j,k)) + &
                         fIn(LB_NODE(6,i,j,k))  + fIn(LB_NODE(7,i,j,k))  + &
                         fIn(LB_NODE(8,i,j,k))  + fIn(LB_NODE(9,i,j,k))  + fIn(LB_NODE(10,i,j,k)) + & 
                         fIn(LB_NODE(11,i,j,k)) + fIn(LB_NODE(12,i,j,k)) + &
                         fIn(LB_NODE(13,i,j,k)) + fIn(LB_NODE(14,i,j,k)) + fIn(LB_NODE(15,i,j,k)) + & 
                         fIn(LB_NODE(16,i,j,k)) + fIn(LB_NODE(17,i,j,k)) + &
                         fIn(LB_NODE(18,i,j,k)) + fIn(LB_NODE(19,i,j,k))

               u(1) = (  fIn(LB_NODE(2,i,j,k))  - fIn(LB_NODE(3,i,j,k))  + &
                         fIn(LB_NODE(8,i,j,k))  - fIn(LB_NODE(9,i,j,k))  + &
                         fIn(LB_NODE(10,i,j,k)) - fIn(LB_NODE(11,i,j,k)) + & 
                         fIn(LB_NODE(12,i,j,k)) - fIn(LB_NODE(13,i,j,k)) + &
                         fIn(LB_NODE(14,i,j,k)) - fIn(LB_NODE(15,i,j,k))  & 
                             )/rho0

               u(2) = (  fIn(LB_NODE(4,i,j,k))  - fIn(LB_NODE(5,i,j,k)) + &
                         fIn(LB_NODE(8,i,j,k))  - fIn(LB_NODE(9,i,j,k))  - fIn(LB_NODE(10,i,j,k)) + & 
                         fIn(LB_NODE(11,i,j,k)) + & 
                         fIn(LB_NODE(16,i,j,k)) - fIn(LB_NODE(17,i,j,k)) + &
                         fIn(LB_NODE(18,i,j,k)) - fIn(LB_NODE(19,i,j,k)))/rho0    

               u(3) = (  fIn(LB_NODE(6,i,j,k))  - fIn(LB_NODE(7,i,j,k))  + &
                         fIn(LB_NODE(12,i,j,k)) - fIn(LB_NODE(13,i,j,k))- fIn(LB_NODE(14,i,j,k))  + & 
                         fIn(LB_NODE(15,i,j,k)) + fIn(LB_NODE(16,i,j,k)) - &
                         fIn(LB_NODE(17,i,j,k)) - fIn(LB_NODE(18,i,j,k)) + &
                         fIn(LB_NODE(19,i,j,k)))/rho0
 
               usq =  (u(1)**2 + u(2)**2 + u(3)**2)*t2cs2inv

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

               fOut(LB_NODE(1,i,j,k)) = fIn(LB_NODE(1,i,j,k)) - omega*(fIn(LB_NODE(1,i,j,k)) - fEq(1))
               fOut(LB_NODE(2,i,j,k)) = fIn(LB_NODE(2,i,j,k)) - omega*(fIn(LB_NODE(2,i,j,k)) - fEq(2))
               fOut(LB_NODE(3,i,j,k)) = fIn(LB_NODE(3,i,j,k)) - omega*(fIn(LB_NODE(3,i,j,k)) - fEq(3))
               fOut(LB_NODE(4,i,j,k)) = fIn(LB_NODE(4,i,j,k)) - omega*(fIn(LB_NODE(4,i,j,k)) - fEq(4))
               fOut(LB_NODE(5,i,j,k)) = fIn(LB_NODE(5,i,j,k)) - omega*(fIn(LB_NODE(5,i,j,k)) - fEq(5))
               fOut(LB_NODE(6,i,j,k)) = fIn(LB_NODE(6,i,j,k)) - omega*(fIn(LB_NODE(6,i,j,k)) - fEq(6))
               fOut(LB_NODE(7,i,j,k)) = fIn(LB_NODE(7,i,j,k)) - omega*(fIn(LB_NODE(7,i,j,k)) - fEq(7))
               fOut(LB_NODE(8,i,j,k)) = fIn(LB_NODE(8,i,j,k)) - omega*(fIn(LB_NODE(8,i,j,k)) - fEq(8))
               fOut(LB_NODE(9,i,j,k)) = fIn(LB_NODE(9,i,j,k)) - omega*(fIn(LB_NODE(9,i,j,k)) - fEq(9))
               fOut(LB_NODE(10,i,j,k)) = fIn(LB_NODE(10,i,j,k)) - omega*(fIn(LB_NODE(10,i,j,k)) - fEq(10))
               fOut(LB_NODE(11,i,j,k)) = fIn(LB_NODE(11,i,j,k)) - omega*(fIn(LB_NODE(11,i,j,k)) - fEq(11))
               fOut(LB_NODE(12,i,j,k)) = fIn(LB_NODE(12,i,j,k)) - omega*(fIn(LB_NODE(12,i,j,k)) - fEq(12))
               fOut(LB_NODE(13,i,j,k)) = fIn(LB_NODE(13,i,j,k)) - omega*(fIn(LB_NODE(13,i,j,k)) - fEq(13))
               fOut(LB_NODE(14,i,j,k)) = fIn(LB_NODE(14,i,j,k)) - omega*(fIn(LB_NODE(14,i,j,k)) - fEq(14))
               fOut(LB_NODE(15,i,j,k)) = fIn(LB_NODE(15,i,j,k)) - omega*(fIn(LB_NODE(15,i,j,k)) - fEq(15))
               fOut(LB_NODE(16,i,j,k)) = fIn(LB_NODE(16,i,j,k)) - omega*(fIn(LB_NODE(16,i,j,k)) - fEq(16))
               fOut(LB_NODE(17,i,j,k)) = fIn(LB_NODE(17,i,j,k)) - omega*(fIn(LB_NODE(17,i,j,k)) - fEq(17))
               fOut(LB_NODE(18,i,j,k)) = fIn(LB_NODE(18,i,j,k)) - omega*(fIn(LB_NODE(18,i,j,k)) - fEq(18))
               fOut(LB_NODE(19,i,j,k)) = fIn(LB_NODE(19,i,j,k)) - omega*(fIn(LB_NODE(19,i,j,k)) - fEq(19))
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
   write(*,*) fIn(LB_NODE(l,5,8,5)) 
enddo
write(*,*) "fOut"
do l=1,nnod
   write(*,*) fOut(LB_NODE(l,5,8,5)) 
enddo
endif
#endif /* DEBUG */
#endif
#ifdef D2Q9
              rholoc =   fIn(LB_NODE(1,i,j,k)) + fIn(LB_NODE(2,i,j,k)) + & 
                         fIn(LB_NODE(3,i,j,k) )+ fIn(LB_NODE(4,i,j,k)) + &
                         fIn(LB_NODE(5,i,j,k)) + fIn(LB_NODE(6,i,j,k)) + & 
                         fIn(LB_NODE(7,i,j,k)) + fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(9,i,j,k) )
              if(s_par%initial .eqv. .false.) then
               u(1) = 1._R8B/rholoc * (fIn(LB_NODE(2,i,j,k))  - fIn(LB_NODE(4,i,j,k)) + &
                         fIn(LB_NODE(6,i,j,k)) - fIn(LB_NODE(7,i,j,k)) - & 
                         fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(9,i,j,k)))
               u(2) = 1._R8B/rholoc * (fIn(LB_NODE(3,i,j,k))  - fIn(LB_NODE(5,i,j,k)) + &
                        fIn(LB_NODE(6,i,j,k)) + fIn(LB_NODE(7,i,j,k)) - & 
                        fIn(LB_NODE(8,i,j,k)) - fIn(LB_NODE(9,i,j,k)))
               else
                  u(1) = lb_dom%u0(LB_NODE(1,i,j,k))
                  u(2) = lb_dom%u0(LB_NODE(2,i,j,k))
               endif

                usq =  (u(1)**2 + u(2)**2)*t2cs2inv


!              if(s_par%initial .eqv. .false.) then
                fEq(1) = t(1)*rholoc *(1._R8B - usq)

                fEq(2) = t(2)*rholoc *(1._R8B + (u(1))*cs2inv &
                      + (u(1))**2*t2cs4inv & 
                      - usq)
                fEq(3) = t(3)*rholoc*(1._R8B + (u(2))*cs2inv &
                      + (u(2))**2*t2cs4inv  &
                      - usq)
                fEq(4) = t(4)*rholoc*(1._R8B + (-u(1) )*cs2inv &
                      + (-u(1))**2*t2cs4inv  &
                      - usq)
                fEq(5) = t(5)*rholoc*(1._R8B + (-u(2))*cs2inv &
                      + (-u(2))**2*t2cs4inv  &
                      - usq)
                fEq(6) = t(6)*rholoc*(1._R8B + (u(1) + u(2))*cs2inv &
                      + (u(1)+u(2))**2*t2cs4inv  &
                      - usq)
                fEq(7) = t(7)*rholoc*(1._R8B + (-u(1) + u(2))*cs2inv &
                      + (-u(1)+u(2))**2*t2cs4inv  &
                      - usq)
                fEq(8) = t(8)*rholoc*(1._R8B + (-u(1) -u(2))*cs2inv &
                      + (-u(1)-u(2))**2*t2cs4inv  &
                      - usq)
                fEq(9) = t(9)*rholoc*(1._R8B + (u(1) -u(2))*cs2inv &
                      + (u(1)-u(2))**2*t2cs4inv  &
                      - usq)
!               else
!                do l=1,nnod
!                  fEq(l) = t(l)*(rholoc+cs2inv*rho0*(real(cx(l,1))*u(1) + real(cx(l,2))*u(2)) &
!                              + t2cs4inv*rho0*((real(cx(l,1))*u(1) + real(cx(l,2))*u(2))**2 - cs2*(u(1)*u(1) + u(2)*u(2))))
!                enddo 
!               endif

               fOut(LB_NODE(1,i,j,k)) = fIn(LB_NODE(1,i,j,k)) - omega*(fIn(LB_NODE(1,i,j,k)) - fEq(1))
               fOut(LB_NODE(2,i,j,k)) = fIn(LB_NODE(2,i,j,k)) - omega*(fIn(LB_NODE(2,i,j,k)) - fEq(2))
               fOut(LB_NODE(3,i,j,k)) = fIn(LB_NODE(3,i,j,k)) - omega*(fIn(LB_NODE(3,i,j,k)) - fEq(3))
               fOut(LB_NODE(4,i,j,k)) = fIn(LB_NODE(4,i,j,k)) - omega*(fIn(LB_NODE(4,i,j,k)) - fEq(4))
               fOut(LB_NODE(5,i,j,k)) = fIn(LB_NODE(5,i,j,k)) - omega*(fIn(LB_NODE(5,i,j,k)) - fEq(5))
               fOut(LB_NODE(6,i,j,k)) = fIn(LB_NODE(6,i,j,k)) - omega*(fIn(LB_NODE(6,i,j,k)) - fEq(6))
               fOut(LB_NODE(7,i,j,k)) = fIn(LB_NODE(7,i,j,k)) - omega*(fIn(LB_NODE(7,i,j,k)) - fEq(7))
               fOut(LB_NODE(8,i,j,k)) = fIn(LB_NODE(8,i,j,k)) - omega*(fIn(LB_NODE(8,i,j,k)) - fEq(8))
               fOut(LB_NODE(9,i,j,k)) = fIn(LB_NODE(9,i,j,k)) - omega*(fIn(LB_NODE(9,i,j,k)) - fEq(9))
 

!if(j==53 .and. i==58)   then
!   rholoc=0
!   do l=1,nnod
!      rholoc=rholoc+fOut(LB_NODE(l,i,j,k))
!   enddo
!      ux   =            (fOut(LB_NODE(2,i,j,k))  - fOut(LB_NODE(4,i,j,k)) + &
!                         fOut(LB_NODE(6,i,j,k)) - fOut(LB_NODE(7,i,j,k)) - & 
!                         fOut(LB_NODE(8,i,j,k)) + fOut(LB_NODE(9,i,j,k)))
!      uy   =            (fOut(LB_NODE(3,i,j,k))  - fOut(LB_NODE(5,i,j,k)) + &
!                        fOut(LB_NODE(6,i,j,k)) + fOut(LB_NODE(7,i,j,k)) - & 
!                        fOut(LB_NODE(8,i,j,k)) - fOut(LB_NODE(9,i,j,k)))
!  write(47,'(2i3,a,3e15.6)') i,j,"aft feq",1.0d0-rholoc,ux,lb_dom%u0(LB_NODE(1,i,j,k))
!endif
#endif
endif
            end do
         end do
      end do
      
      call cpu_time_measure(meas%tEnd_comm)

      meas%coll_duration = meas%tEnd_comm - meas%tSt_comm + meas%coll_duration


   end subroutine collide_bgk 
  !------------------------------------------------------------------------


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
      real(R8B), intent(out)     :: fOut(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(inout)   ::  fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
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

                  ftmp1(INDEXVAR)=fIn(LB_NODE(1,i-cx(1,1),j-cx(1,2),k-cx(1,3)))
                  ftmp2(INDEXVAR)=fIn(LB_NODE(2,i-cx(2,1),j-cx(2,2),k-cx(2,3)))
                  ftmp3(INDEXVAR)=fIn(LB_NODE(3,i-cx(3,1),j-cx(3,2),k-cx(3,3)))
                  ftmp4(INDEXVAR)=fIn(LB_NODE(4,i-cx(4,1),j-cx(4,2),k-cx(4,3)))
                  ftmp5(INDEXVAR)=fIn(LB_NODE(5,i-cx(5,1),j-cx(5,2),k-cx(5,3)))
                  ftmp6(INDEXVAR)=fIn(LB_NODE(6,i-cx(6,1),j-cx(6,2),k-cx(6,3)))
                  ftmp7(INDEXVAR)=fIn(LB_NODE(7,i-cx(7,1),j-cx(7,2),k-cx(7,3)))
                  ftmp8(INDEXVAR)=fIn(LB_NODE(8,i-cx(8,1),j-cx(8,2),k-cx(8,3)))
                  ftmp9(INDEXVAR)=fIn(LB_NODE(9,i-cx(9,1),j-cx(9,2),k-cx(9,3)))
#ifdef D3Q19
                  ftmp10(INDEXVAR)=fIn(LB_NODE(10,i-cx(10,1),j-cx(10,2),k-cx(10,3)))
                  ftmp11(INDEXVAR)=fIn(LB_NODE(11,i-cx(11,1),j-cx(11,2),k-cx(11,3)))
                  ftmp12(INDEXVAR)=fIn(LB_NODE(12,i-cx(12,1),j-cx(12,2),k-cx(12,3)))
                  ftmp13(INDEXVAR)=fIn(LB_NODE(13,i-cx(13,1),j-cx(13,2),k-cx(13,3)))
                  ftmp14(INDEXVAR)=fIn(LB_NODE(14,i-cx(14,1),j-cx(14,2),k-cx(14,3)))
                  ftmp15(INDEXVAR)=fIn(LB_NODE(15,i-cx(15,1),j-cx(15,2),k-cx(15,3)))
                  ftmp16(INDEXVAR)=fIn(LB_NODE(16,i-cx(16,1),j-cx(16,2),k-cx(16,3)))
                  ftmp17(INDEXVAR)=fIn(LB_NODE(17,i-cx(17,1),j-cx(17,2),k-cx(17,3)))
                  ftmp18(INDEXVAR)=fIn(LB_NODE(18,i-cx(18,1),j-cx(18,2),k-cx(18,3)))
                  ftmp19(INDEXVAR)=fIn(LB_NODE(19,i-cx(19,1),j-cx(19,2),k-cx(19,3)))
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


               fOut(LB_NODE(1,i,j,k)) = ftmp1(INDEXVAR) - omega(INDEXVAR)*(ftmp1(INDEXVAR) - fEq1)
               fOut(LB_NODE(2,i,j,k)) = ftmp2(INDEXVAR) - omega(INDEXVAR)*(ftmp2(INDEXVAR) - fEq2)
               fOut(LB_NODE(3,i,j,k)) = ftmp3(INDEXVAR) - omega(INDEXVAR)*(ftmp3(INDEXVAR) - fEq3)
               fOut(LB_NODE(4,i,j,k)) = ftmp4(INDEXVAR) - omega(INDEXVAR)*(ftmp4(INDEXVAR) - fEq4)
               fOut(LB_NODE(5,i,j,k)) = ftmp5(INDEXVAR) - omega(INDEXVAR)*(ftmp5(INDEXVAR) - fEq5)

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
               fOut(LB_NODE(6,i,j,k)) = ftmp6(INDEXVAR) - omega(INDEXVAR)*(ftmp6(INDEXVAR) - fEq6)
               fOut(LB_NODE(7,i,j,k)) = ftmp7(INDEXVAR) - omega(INDEXVAR)*(ftmp7(INDEXVAR) - fEq7)
               fOut(LB_NODE(8,i,j,k)) = ftmp8(INDEXVAR) - omega(INDEXVAR)*(ftmp8(INDEXVAR) - fEq8)
               fOut(LB_NODE(9,i,j,k)) = ftmp9(INDEXVAR) - omega(INDEXVAR)*(ftmp9(INDEXVAR) - fEq9)
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
               fOut(LB_NODE(1,i,j,k)) = ftmp1(INDEXVAR) - omega(INDEXVAR)*(ftmp1(INDEXVAR) - fEq1)
               fOut(LB_NODE(2,i,j,k)) = ftmp2(INDEXVAR) - omega(INDEXVAR)*(ftmp2(INDEXVAR) - fEq2)
               fOut(LB_NODE(3,i,j,k)) = ftmp3(INDEXVAR) - omega(INDEXVAR)*(ftmp3(INDEXVAR) - fEq3)
               fOut(LB_NODE(4,i,j,k)) = ftmp4(INDEXVAR) - omega(INDEXVAR)*(ftmp4(INDEXVAR) - fEq4)
               fOut(LB_NODE(5,i,j,k)) = ftmp5(INDEXVAR) - omega(INDEXVAR)*(ftmp5(INDEXVAR) - fEq5)
               fOut(LB_NODE(6,i,j,k)) = ftmp6(INDEXVAR) - omega(INDEXVAR)*(ftmp6(INDEXVAR) - fEq6)
               fOut(LB_NODE(7,i,j,k)) = ftmp7(INDEXVAR) - omega(INDEXVAR)*(ftmp7(INDEXVAR) - fEq7)

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
               fOut(LB_NODE(8,i,j,k))  = ftmp8(INDEXVAR)  - omega(INDEXVAR)*(ftmp8(INDEXVAR)  -  fEq8)
               fOut(LB_NODE(9,i,j,k))  = ftmp9(INDEXVAR)  - omega(INDEXVAR)*(ftmp9(INDEXVAR)  -  fEq9)
               fOut(LB_NODE(10,i,j,k)) = ftmp10(INDEXVAR) - omega(INDEXVAR)*(ftmp10(INDEXVAR) - fEq10)
               fOut(LB_NODE(11,i,j,k)) = ftmp11(INDEXVAR) - omega(INDEXVAR)*(ftmp11(INDEXVAR) - fEq11)
               fOut(LB_NODE(12,i,j,k)) = ftmp12(INDEXVAR) - omega(INDEXVAR)*(ftmp12(INDEXVAR) - fEq12)
               fOut(LB_NODE(13,i,j,k)) = ftmp13(INDEXVAR) - omega(INDEXVAR)*(ftmp13(INDEXVAR) - fEq13)

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
               fOut(LB_NODE(14,i,j,k)) = ftmp14(INDEXVAR) - omega(INDEXVAR)*(ftmp14(INDEXVAR) - fEq14)
               fOut(LB_NODE(15,i,j,k)) = ftmp15(INDEXVAR) - omega(INDEXVAR)*(ftmp15(INDEXVAR) - fEq15)
               fOut(LB_NODE(16,i,j,k)) = ftmp16(INDEXVAR) - omega(INDEXVAR)*(ftmp16(INDEXVAR) - fEq16)
               fOut(LB_NODE(17,i,j,k)) = ftmp17(INDEXVAR) - omega(INDEXVAR)*(ftmp17(INDEXVAR) - fEq17)
               fOut(LB_NODE(18,i,j,k)) = ftmp18(INDEXVAR) - omega(INDEXVAR)*(ftmp18(INDEXVAR) - fEq18)
               fOut(LB_NODE(19,i,j,k)) = ftmp19(INDEXVAR) - omega(INDEXVAR)*(ftmp19(INDEXVAR) - fEq19)



#endif

#ifdef OLD
         !--------------------------------------------
         ! Relaxation Process 

               fOut(LB_NODE(1,i,j,k)) = ftmp1(INDEXVAR) - omega(INDEXVAR)*(ftmp1(INDEXVAR) - fEq1)
               fOut(LB_NODE(2,i,j,k)) = ftmp2(INDEXVAR) - omega(INDEXVAR)*(ftmp2(INDEXVAR) - fEq2)
               fOut(LB_NODE(3,i,j,k)) = ftmp3(INDEXVAR) - omega(INDEXVAR)*(ftmp3(INDEXVAR) - fEq3)
               fOut(LB_NODE(4,i,j,k)) = ftmp4(INDEXVAR) - omega(INDEXVAR)*(ftmp4(INDEXVAR) - fEq4)
               fOut(LB_NODE(5,i,j,k)) = ftmp5(INDEXVAR) - omega(INDEXVAR)*(ftmp5(INDEXVAR) - fEq5)
               fOut(LB_NODE(6,i,j,k)) = ftmp6(INDEXVAR) - omega(INDEXVAR)*(ftmp6(INDEXVAR) - fEq6)
               fOut(LB_NODE(7,i,j,k)) = ftmp7(INDEXVAR) - omega(INDEXVAR)*(ftmp7(INDEXVAR) - fEq7)
               fOut(LB_NODE(8,i,j,k)) = ftmp8(INDEXVAR) - omega(INDEXVAR)*(ftmp8(INDEXVAR) - fEq8)
               fOut(LB_NODE(9,i,j,k)) = ftmp9(INDEXVAR) - omega(INDEXVAR)*(ftmp9(INDEXVAR) - fEq9)
#ifdef D3Q19
               fOut(LB_NODE(10,i,j,k)) = ftmp10(INDEXVAR) - omega(INDEXVAR)*(ftmp10(INDEXVAR) - fEq10)
               fOut(LB_NODE(11,i,j,k)) = ftmp11(INDEXVAR) - omega(INDEXVAR)*(ftmp11(INDEXVAR) - fEq11)
               fOut(LB_NODE(12,i,j,k)) = ftmp12(INDEXVAR) - omega(INDEXVAR)*(ftmp12(INDEXVAR) - fEq12)
               fOut(LB_NODE(13,i,j,k)) = ftmp13(INDEXVAR) - omega(INDEXVAR)*(ftmp13(INDEXVAR) - fEq13)
               fOut(LB_NODE(14,i,j,k)) = ftmp14(INDEXVAR) - omega(INDEXVAR)*(ftmp14(INDEXVAR) - fEq14)
               fOut(LB_NODE(15,i,j,k)) = ftmp15(INDEXVAR) - omega(INDEXVAR)*(ftmp15(INDEXVAR) - fEq15)
               fOut(LB_NODE(16,i,j,k)) = ftmp16(INDEXVAR) - omega(INDEXVAR)*(ftmp16(INDEXVAR) - fEq16)
               fOut(LB_NODE(17,i,j,k)) = ftmp17(INDEXVAR) - omega(INDEXVAR)*(ftmp17(INDEXVAR) - fEq17)
               fOut(LB_NODE(18,i,j,k)) = ftmp18(INDEXVAR) - omega(INDEXVAR)*(ftmp18(INDEXVAR) - fEq18)
               fOut(LB_NODE(19,i,j,k)) = ftmp19(INDEXVAR) - omega(INDEXVAR)*(ftmp19(INDEXVAR) - fEq19)
#endif
#endif /* OLD */
else

         !--------------------------------------------
         ! Do Bounce Back 
#ifdef D2Q9
         fOut(LB_NODE(opp(1),i,j,k)) = ftmp1(INDEXVAR)
         fOut(LB_NODE(opp(2),i,j,k)) = ftmp2(INDEXVAR)
         fOut(LB_NODE(opp(3),i,j,k)) = ftmp3(INDEXVAR)
         fOut(LB_NODE(opp(4),i,j,k)) = ftmp4(INDEXVAR)
         fOut(LB_NODE(opp(5),i,j,k)) = ftmp5(INDEXVAR)

#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
#endif


         fOut(LB_NODE(opp(6),i,j,k)) = ftmp6(INDEXVAR)
         fOut(LB_NODE(opp(7),i,j,k)) = ftmp7(INDEXVAR)
         fOut(LB_NODE(opp(8),i,j,k)) = ftmp8(INDEXVAR)
         fOut(LB_NODE(opp(9),i,j,k)) = ftmp9(INDEXVAR)
#endif

#ifdef D3Q19
         fOut(LB_NODE(opp(1),i,j,k)) = ftmp1(INDEXVAR)
         fOut(LB_NODE(opp(2),i,j,k)) = ftmp2(INDEXVAR)
         fOut(LB_NODE(opp(3),i,j,k)) = ftmp3(INDEXVAR)
         fOut(LB_NODE(opp(4),i,j,k)) = ftmp4(INDEXVAR)
         fOut(LB_NODE(opp(5),i,j,k)) = ftmp5(INDEXVAR)
         fOut(LB_NODE(opp(6),i,j,k)) = ftmp6(INDEXVAR)
         fOut(LB_NODE(opp(7),i,j,k)) = ftmp7(INDEXVAR)
#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
#endif
         fOut(LB_NODE(opp(8),i,j,k)) = ftmp8(INDEXVAR)
         fOut(LB_NODE(opp(9),i,j,k)) = ftmp9(INDEXVAR)
         fOut(LB_NODE(opp(10),i,j,k)) = ftmp10(INDEXVAR)
         fOut(LB_NODE(opp(11),i,j,k)) = ftmp11(INDEXVAR)
         fOut(LB_NODE(opp(12),i,j,k)) = ftmp12(INDEXVAR)
         fOut(LB_NODE(opp(13),i,j,k)) = ftmp13(INDEXVAR)
#ifdef SPLITLOOPS
            endif
         enddo
         do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall)) then
#endif
         fOut(LB_NODE(opp(14),i,j,k)) = ftmp14(INDEXVAR)
         fOut(LB_NODE(opp(15),i,j,k)) = ftmp15(INDEXVAR)
         fOut(LB_NODE(opp(16),i,j,k)) = ftmp16(INDEXVAR)
         fOut(LB_NODE(opp(17),i,j,k)) = ftmp17(INDEXVAR)
         fOut(LB_NODE(opp(18),i,j,k)) = ftmp18(INDEXVAR)
         fOut(LB_NODE(opp(19),i,j,k)) = ftmp19(INDEXVAR)
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







#ifdef MRT_ILBDC


   subroutine collide(lb_dom,fIn,fOut,s_par,meas)



!------------------------------------------------------------------------
! collision step MRT model
! not working yet!



      implicit none
      type(lb_block),intent(in)  :: lb_dom
      type(sim_parameter),intent(in)  :: s_par
      type(measure      )             :: meas 
      real(R8B)   :: usq,om,M_tr(nnod,nnod),mrt_m(nnod),temp,M_inv(nnod,nnod),s_mrt(nnod),mrt_out(nnod),mrt_m0(nnod)
      real (R8B) :: d_tmp,id_tmp,d1,d2,d3,d4,d5,d6
   real (R8B) :: m_rho,m_e,m_eps,m_jx,m_qx,m_jy,m_qy,m_jz,m_qz,m_pxx,m_pixx,m_pww,m_piww,m_pxy,m_pyz,m_pxz,m_mx,m_my,m_mz
  real (R8B) :: meq_rho,meq_e,meq_eps,meq_jx,meq_qx,meq_jy,meq_qy,meq_jz,meq_qz,meq_pxx,meq_pixx,meq_pww,meq_piww,meq_pxy
  real (R8B) :: meq_pyz,meq_pxz,meq_mx,meq_my,meq_mz 
  real (R8B) :: s_rho,s_e,s_eps,s_jx,s_qx,s_jy,s_qy,s_jz,s_qz,s_pxx,s_pixx,s_pww,s_piww,s_pxy,s_pyz,s_pxz,s_mx,s_my,s_mz
  real (R8B), parameter :: w_eps  = 0._R8B
  real (R8B), parameter :: w_epsj = - 475._R8B/63._R8B
  real (R8B), parameter :: w_xx   = 0._R8B
  real(R8B), intent(out)     :: fOut(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
  real(R8B), intent(inout)   :: fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer               :: i,j,k,l,ll
      om = s_par%omega
!#!define COLL_MRT
#ifdef COLL_MRT
      ! Multi-relaxation time collision routine
#ifdef D2Q9
      M_tr(1:9,1) = (/ 1., 1., 1., 1., 1., 1., 1., 1., 1. /)
      M_tr(1:9,2) = (/-4.,-1.,-1.,-1.,-1., 2., 2., 2., 2. /)
      M_tr(1:9,3) = (/ 4.,-2.,-2.,-2.,-2., 1., 1., 1., 1. /)
      M_tr(1:9,4) = (/ 0., 1., 0.,-1., 0., 1.,-1.,-1., 1. /)
      M_tr(1:9,5) = (/ 0.,-2., 0., 2., 0., 1.,-1.,-1., 1. /)
      M_tr(1:9,6) = (/ 0., 0., 1., 0.,-1., 1., 1.,-1.,-1. /)
      M_tr(1:9,7) = (/ 0., 0.,-2., 0., 2., 1., 1.,-1.,-1. /)
      M_tr(1:9,8) = (/ 0., 1.,-1., 1.,-1., 0., 0., 0., 0. /)
      M_tr(1:9,9) = (/ 0., 0., 0., 0., 0., 1.,-1., 1.,-1. /)

      M_inv(1:9,1) = (/ 1./9., -1./9., 1./9.,0.,0.,0.,0.,0.,0. /)
      M_inv(1:9,2) = (/ 1./9., -1./36., -1./18., 1./6., -1./6., 0., 0., 1./4., 0. /)
      M_inv(1:9,3) = (/ 1./9., -1./36., -1./18., 0., 0., 1./6., -1./6., -1./4., 0. /)
      M_inv(1:9,4) = (/ 1./9., -1./36., -1./18.,-1./6.,  1./6., 0., 0., 1./4., 0. /)
      M_inv(1:9,5) = (/ 1./9., -1./36., -1./18., 0., 0.,-1./6.,  1./6., -1./4., 0. /)
      M_inv(1:9,6) = (/ 1./9., 1./18., 1./36., 1./6., 1./12., 1./6., 1./12., 0., 1./4. /)
      M_inv(1:9,7) = (/ 1./9., 1./18., 1./36.,-1./6.,-1./12., 1./6., 1./12., 0.,-1./4. /)
      M_inv(1:9,8) = (/ 1./9., 1./18., 1./36.,-1./6.,-1./12.,-1./6.,-1./12., 0., 1./4. /)
      M_inv(1:9,9) = (/ 1./9., 1./18., 1./36., 1./6., 1./12.,-1./6.,-1./12., 0.,-1./4. /)

      s_mrt(1:nnod) = 0.d0
      s_mrt(2)      = 1.63
      s_mrt(3)      = 1.14
      s_mrt(5)      = 1.92
      s_mrt(7)      = 1.92
      s_mrt(8)      = om
      s_mrt(9)      = om
     do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)

write(*,*) "the m is ",mrt_m
               do l=1,nnod 
                  mrt_m(l) = 0.d0
                   
                  do ll=1,nnod
                     mrt_m(l) = mrt_m(l) + M_inv(l,ll)*M_tr(l,ll) * fIn(LB_NODE(ll,i,j,k)) 
                  enddo
               enddo
write(*,*) "the m is ",mrt_m
               mrt_m0(:) = 0.d0
               mrt_m0(2) = -2.d0*mrt_m(1) + 3.d0*(mrt_m(4)**2 + mrt_m(6)**2)
               mrt_m0(3) =       mrt_m(1) - 3.d0*(mrt_m(4)**2 + mrt_m(6)**2)
               mrt_m0(5) =                   - mrt_m(4)
               mrt_m0(7) =                   - mrt_m(6)
               mrt_m0(8) =   mrt_m(4)**2 - mrt_m(6)**2 
               mrt_m0(9) =   mrt_m(4) * mrt_m(6)
write(*,*) "M0-mrt is ",mrt_m0
               do l=1,nnod 
               mrt_out(l) = mrt_m(l) - s_mrt(l) * (mrt_m(l) - mrt_m0(l)  )
                  temp = 0.d0
                   
                  do ll=1,nnod
                     temp = temp + M_inv(l,ll) * mrt_out(l) 
                  enddo
                  fOut(LB_NODE(l,i,j,k)) =  temp
               enddo
write(*,*) "RESULT FOR fOUT is",fOut(LB_NODE(:,i,j,k))
if(i==3) stop
            end do
         end do
      end do
#endif /* D2Q9 */
#ifdef D3Q19
!FIXME check values! have to be same as BGK with omega in every s-spot.
  s_rho  = 0._R8B
  s_e    = 1.19_R8B
  s_eps  = 1.4_R8B
  s_jx   = 0._R8B
  s_qx   = 1.2_R8B
  s_jy   = 0._R8B
  s_qy   = 1.2_R8B
  s_jz   = 0._R8B
  s_qz   = 1.2_R8B
  s_pxx  = om
  s_pixx = 1.4_R8B
  s_pww  = om
  s_piww = 1.4_R8B
  s_pxy  = om
  s_pyz  = om
  s_pxz  = om
  s_mx   = 1.98_R8B
  s_my   = 1.98_R8B
  s_mz   = 1.98_R8B
     do k=1,lb_dom%lx(3)
         do j=1,lb_dom%lx(2)
            do i=1,lb_dom%lx(1)

           d1 = fIn(LB_NODE(4,i,j,k)) + fIn(LB_NODE(3,i,j,k)) + fIn(LB_NODE(5,i,j,k)) &
              + fIn(LB_NODE(2,i,j,k)) + fIn(LB_NODE(6,i,j,k)) + fIn(LB_NODE(7,i,j,k)) 
           d2 = fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(11,i,j,k)) + fIn(LB_NODE(9,i,j,k)) &
              + fIn(LB_NODE(10,i,j,k)) + fIn(LB_NODE(12,i,j,k)) + fIn(LB_NODE(16,i,j,k)) &
              + fIn(LB_NODE(15,i,j,k)) + fIn(LB_NODE(19,i,j,k)) + fIn(LB_NODE(14,i,j,k)) &
              + fIn(LB_NODE(18,i,j,k)) + fIn(LB_NODE(13,i,j,k)) + fIn(LB_NODE(17,i,j,k))
           !
           m_rho =            fIn(LB_NODE(1,i,j,k)) +          d1 +         d2 
           m_e   = - 30._R8B * fIn(LB_NODE(1,i,j,k)) - 11._R8B * d1 + 8._R8B * d2
           m_eps =   12._R8B * fIn(LB_NODE(1,i,j,k)) -  4._R8B * d1 +         d2
           !
           d1 = fIn(LB_NODE(8,i,j,k)) - fIn(LB_NODE(11,i,j,k)) - fIn(LB_NODE(9,i,j,k)) &
              + fIn(LB_NODE(10,i,j,k)) + fIn(LB_NODE(12,i,j,k)) - fIn(LB_NODE(15,i,j,k)) &
              + fIn(LB_NODE(14,i,j,k)) - fIn(LB_NODE(13,i,j,k))
           !
           m_jx =         - fIn(LB_NODE(3,i,j,k)) + fIn(LB_NODE(2,i,j,k))   + d1 
           m_qx = 4._R8B * ( fIn(LB_NODE(3,i,j,k)) - fIn(LB_NODE(2,i,j,k)) ) + d1
           !
           d1 = fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(11,i,j,k)) - fIn(LB_NODE(9,i,j,k)) &
              - fIn(LB_NODE(10,i,j,k)) + fIn(LB_NODE(16,i,j,k)) - fIn(LB_NODE(19,i,j,k)) &
              + fIn(LB_NODE(18,i,j,k)) - fIn(LB_NODE(17,i,j,k))
           !
           m_jy =             fIn(LB_NODE(4,i,j,k)) - fIn(LB_NODE(5,i,j,k))   + d1 
           m_qy = - 4._R8B * ( fIn(LB_NODE(4,i,j,k)) - fIn(LB_NODE(5,i,j,k)) ) + d1
           !
           d1 = fIn(LB_NODE(12,i,j,k)) + fIn(LB_NODE(16,i,j,k)) + fIn(LB_NODE(15,i,j,k)) &
              + fIn(LB_NODE(19,i,j,k)) - fIn(LB_NODE(14,i,j,k)) - fIn(LB_NODE(18,i,j,k)) &
              - fIn(LB_NODE(13,i,j,k)) - fIn(LB_NODE(17,i,j,k))
           !
           m_jz =             fIn(LB_NODE(6,i,j,k)) - fIn(LB_NODE(7,i,j,k))   + d1
           m_qz = - 4._R8B * ( fIn(LB_NODE(6,i,j,k)) - fIn(LB_NODE(7,i,j,k)) ) + d1
           !
           d1 = - 2._R8B * ( fIn(LB_NODE(3,i,j,k)) + fIn(LB_NODE(2,i,j,k)) )      &
              + fIn(LB_NODE(4,i,j,k)) + fIn(LB_NODE(5,i,j,k)) + fIn(LB_NODE(6,i,j,k)) &
              + fIn(LB_NODE(7,i,j,k))
           d2 = fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(11,i,j,k)) + fIn(LB_NODE(9,i,j,k)) &
              + fIn(LB_NODE(10,i,j,k)) 
           d3 = fIn(LB_NODE(12,i,j,k)) + fIn(LB_NODE(15,i,j,k)) + fIn(LB_NODE(14,i,j,k)) &
              + fIn(LB_NODE(13,i,j,k))
           d4 = - 2._R8B * ( fIn(LB_NODE(16,i,j,k)) + fIn(LB_NODE(19,i,j,k)) &
                          + fIn(LB_NODE(18,i,j,k)) + fIn(LB_NODE(17,i,j,k)) )
           !
           m_pxx  = -         d1 + d2 + d3 + d4
           m_pixx = + 2._R8B * d1 + d2 + d3 + d4 
           !
           d1 = fIn(LB_NODE(4,i,j,k)) + fIn(LB_NODE(5,i,j,k)) - fIn(LB_NODE(6,i,j,k)) &
              - fIn(LB_NODE(7,i,j,k))
           !
           m_pww  =           d1 + d2 - d3
           m_piww = - 2._R8B * d1 + d2 - d3
           !
           m_pxy = fIn(LB_NODE(8,i,j,k)) - fIn(LB_NODE(11,i,j,k)) + fIn(LB_NODE(9,i,j,k)) &
                 - fIn(LB_NODE(10,i,j,k)) 
           m_pyz = fIn(LB_NODE(16,i,j,k)) - fIn(LB_NODE(19,i,j,k)) - fIn(LB_NODE(18,i,j,k)) &
                 + fIn(LB_NODE(17,i,j,k)) 
           m_pxz = fIn(LB_NODE(12,i,j,k)) - fIn(LB_NODE(15,i,j,k)) - fIn(LB_NODE(14,i,j,k)) &
                 + fIn(LB_NODE(13,i,j,k)) 
           m_mx  = fIn(LB_NODE(8,i,j,k)) - fIn(LB_NODE(11,i,j,k)) - fIn(LB_NODE(9,i,j,k)) &
                 + fIn(LB_NODE(10,i,j,k)) - fIn(LB_NODE(12,i,j,k)) + fIn(LB_NODE(15,i,j,k)) &
                 - fIn(LB_NODE(14,i,j,k)) + fIn(LB_NODE(13,i,j,k)) 
           m_my  = fIn(LB_NODE(8,i,j,k)) - fIn(LB_NODE(11,i,j,k)) + fIn(LB_NODE(9,i,j,k)) &
                 + fIn(LB_NODE(10,i,j,k)) + fIn(LB_NODE(16,i,j,k)) - fIn(LB_NODE(19,i,j,k)) &
                 + fIn(LB_NODE(18,i,j,k)) - fIn(LB_NODE(17,i,j,k)) 
           m_mz  = fIn(LB_NODE(12,i,j,k)) - fIn(LB_NODE(16,i,j,k)) + fIn(LB_NODE(15,i,j,k)) &
                 - fIn(LB_NODE(19,i,j,k)) - fIn(LB_NODE(14,i,j,k)) + fIn(LB_NODE(18,i,j,k)) &
                 - fIn(LB_NODE(13,i,j,k)) + fIn(LB_NODE(17,i,j,k)) 
           !
           !
           !
           m_jx = m_jx 
           m_jy = m_jy
           m_jz = m_jz
           !
           usq = ( m_jx*m_jx + m_jy*m_jy + m_jz*m_jz ) / m_rho
           !
           meq_e    = -11._R8B * m_rho + 19._R8B * usq
           meq_eps  =   w_eps * m_rho + w_epsj * usq
           meq_qx   = - 2.d0/3.d0 * m_jx
           meq_qy   = - 2.d0/3.d0 * m_jy
           meq_qz   = - 2.d0/3.d0 * m_jz
           meq_pxx  = (2*m_jx*m_jx-m_jy*m_jy-m_jz*m_jz)/m_rho
           meq_pixx = w_xx*(2*m_jx*m_jx-m_jy*m_jy-m_jz*m_jz)/m_rho
           meq_pww  = (m_jy*m_jy-m_jz*m_jz)/m_rho
           meq_piww = w_xx*(m_jy*m_jy-m_jz*m_jz)/m_rho
           meq_pxy  = m_jx * m_jy / m_rho
           meq_pyz  = m_jy * m_jz / m_rho
           meq_pxz  = m_jx * m_jz / m_rho
           meq_mx   = 0._R8B
           meq_my   = 0._R8B
           meq_mz   = 0._R8B
           !
           m_e    = ( m_e    - s_e    * ( m_e    - meq_e    ) ) * 4._R8B/1197._R8B
           m_eps  = ( m_eps  - s_eps  * ( m_eps  - meq_eps  ) ) / 252._R8B
           m_qx   = ( m_qx   - s_qx   * ( m_qx   - meq_qx   ) ) / 40._R8B
           m_qy   = ( m_qy   - s_qy   * ( m_qy   - meq_qy   ) ) / 40._R8B
           m_qz   = ( m_qz   - s_qz   * ( m_qz   - meq_qz   ) ) / 40._R8B
           m_pxx  = ( m_pxx  - s_pxx  * ( m_pxx  - meq_pxx  ) ) / 36._R8B
           m_pixx = ( m_pixx - s_pixx * ( m_pixx - meq_pixx ) ) / 72._R8B
           m_pww  = ( m_pww  - s_pww  * ( m_pww  - meq_pww  ) ) / 12._R8B
           m_piww = ( m_piww - s_piww * ( m_piww - meq_piww ) ) / 24._R8B
           m_pxy  = ( m_pxy  - s_pxy  * ( m_pxy  - meq_pxy  ) ) / 4._R8B
           m_pyz  = ( m_pyz  - s_pyz  * ( m_pyz  - meq_pyz  ) ) / 4._R8B
           m_pxz  = ( m_pxz  - s_pxz  * ( m_pxz  - meq_pxz  ) ) / 4._R8B
           m_mx   = ( m_mx   - s_mx   * ( m_mx   - meq_mx   ) ) / 8._R8B
           m_my   = ( m_my   - s_my   * ( m_my   - meq_my   ) ) / 8._R8B
           m_mz   = ( m_mz   - s_mz   * ( m_mz   - meq_mz   ) ) / 8._R8B
           !
           ! 61 flops
           !

           !
           m_rho = m_rho / 19._R8B           
           m_jx  = m_jx  / 10._R8B
           m_jy  = m_jy  / 10._R8B
           m_jz  = m_jz  / 10._R8B
           !
           ! 4 flops
           !
           fOut(LB_NODE(1,i,j,k)) = m_rho - 15._R8B/4._R8B * m_e + 12._R8B * m_eps
           !
           d1 = - 11._R8B/8._R8B * m_e - 4._R8B * m_eps
           d2 = m_jx  - 4._R8B * m_qx
           d3 = m_jy  - 4._R8B * m_qy
           d4 = m_jz  - 4._R8B * m_qz 
           d5 = m_pxx - 2._R8B * m_pixx 
           d6 = m_pww - 2._R8B * m_piww
           !
           ! 18 flops
           !
           fOut(LB_NODE(2,i,j,k)) = m_rho + d1 + d2 + 2._R8B * d5      
           fOut(LB_NODE(3,i,j,k)) = m_rho + d1 - d2 + 2._R8B * d5     
           fOut(LB_NODE(4,i,j,k)) = m_rho + d1 + d3 -         d5 + d6 
           fOut(LB_NODE(5,i,j,k)) = m_rho + d1 - d3 -         d5 + d6
           fOut(LB_NODE(6,i,j,k)) = m_rho + d1 + d4 -         d5 - d6
           fOut(LB_NODE(7,i,j,k)) = m_rho + d1 - d4 -         d5 - d6
           !
           d1 = m_e   + m_eps
           d2 = m_jx  + m_qx
           d3 = m_jy  + m_qy
           d4 = m_jz  + m_qz 
           d5 = m_pxx + m_pixx
           d6 = m_pww + m_piww
           !
           ! 36 flops
           !
           fOut(LB_NODE(8,i,j,k)) = m_rho + d1 + d2 + d3 +         d5 + d6 + m_pxy + m_mx - m_my 
           fOut(LB_NODE(11,i,j,k)) = m_rho + d1 - d2 + d3 +         d5 + d6 - m_pxy - m_mx - m_my 
           fOut(LB_NODE(9,i,j,k)) = m_rho + d1 - d2 - d3 +         d5 + d6 + m_pxy - m_mx + m_my 
           fOut(LB_NODE(10,i,j,k)) = m_rho + d1 + d2 - d3 +         d5 + d6 - m_pxy + m_mx + m_my 
           fOut(LB_NODE(12,i,j,k)) = m_rho + d1 + d2 + d4 +         d5 - d6 + m_pxz - m_mx + m_mz 
           fOut(LB_NODE(15,i,j,k)) = m_rho + d1 - d2 + d4 +         d5 - d6 - m_pxz + m_mx + m_mz 
           fOut(LB_NODE(14,i,j,k)) = m_rho + d1 + d2 - d4 +         d5 - d6 - m_pxz - m_mx - m_mz 
           fOut(LB_NODE(13,i,j,k)) = m_rho + d1 - d2 - d4 +         d5 - d6 + m_pxz + m_mx - m_mz 
           !
           fOut(LB_NODE(16,i,j,k)) = m_rho + d1 + d3 + d4 - 2._R8B * d5      + m_pyz + m_my - m_mz 
           fOut(LB_NODE(19,i,j,k)) = m_rho + d1 - d3 + d4 - 2._R8B * d5      - m_pyz - m_my - m_mz 
           fOut(LB_NODE(18,i,j,k)) = m_rho + d1 + d3 - d4 - 2._R8B * d5      - m_pyz + m_my + m_mz 
           fOut(LB_NODE(17,i,j,k)) = m_rho + d1 - d3 - d4 - 2._R8B * d5      + m_pyz - m_my + m_mz 

            end do
         end do
      end do

 write(*,*) "MRT model has errors!!"
#endif /* D3Q19 */
  
#else
#endif /* COLL_MRT */
   end subroutine collide
  !------------------------------------------------------------------------
#endif /* MRT */





#ifdef INIT_WITH_ROOT  




   subroutine calc_fEq_global(rho,u,fEq,s_par)
   !------------------------------------------------------------------------
   ! same as calc_fEq, but with global limits.
   !


      implicit none
      real(R8B)                  :: usq,cs2inv,t2cs4inv,t2cs2inv
      type(sim_parameter)        :: s_par
      real(R8B), intent(in)      :: u(LB_NODE(NDIM,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
      real(R8B), intent(in)      :: rho(s_par%gx(1),s_par%gx(2),s_par%gx(3))
      real(R8B), intent(out)     :: fEq(LB_NODE(nnod,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
      integer               :: i,j,k,l

      cs2inv  = 1._R8B/cs**2
      t2cs4inv = 1._R8B/(2._R8B*cs**4)
      t2cs2inv = 1._R8B/(2._R8B*cs**2)
#ifdef D2Q9    
! ersetze diese rho durch lokale berechnungen und schauen, was sich aendert
      do k=1,s_par%gx(3)
         do j=1,s_par%gx(2)
            do i=1,s_par%gx(1)
                usq =  (u(LB_NODE(1,i,j,k))**2 + u(LB_NODE(2,i,j,k))**2)*t2cs2inv
                fEq(LB_NODE(1,i,j,k)) = t(1)*rho(i,j,k) *(1._R8B - usq)
                fEq(LB_NODE(2,i,j,k)) = t(2)*rho(i,j,k) *(1._R8B + (u(LB_NODE(1,i,j,k)))*cs2inv &
                      + (u(LB_NODE(1,i,j,k)))**2*t2cs4inv   - usq)
                fEq(LB_NODE(3,i,j,k)) = t(3)*rho(i,j,k)*(1._R8B - (-u(LB_NODE(2,i,j,k)))*cs2inv &
                      + (u(LB_NODE(2,i,j,k)))**2*t2cs4inv   - usq)
                fEq(LB_NODE(4,i,j,k)) = t(4)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(1,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(1,i,j,k)))**2*t2cs4inv  - usq)
                fEq(LB_NODE(5,i,j,k)) = t(5)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(2,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(2,i,j,k)))**2*t2cs4inv   - usq)
                fEq(LB_NODE(6,i,j,k)) = t(6)*rho(i,j,k)*(1._R8B + (u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))*cs2inv &
                      + (u(LB_NODE(1,i,j,k))+u(LB_NODE(2,i,j,k)))**2*t2cs4inv  - usq)
                fEq(LB_NODE(7,i,j,k)) = t(7)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(1,i,j,k))+u(LB_NODE(2,i,j,k)))**2*t2cs4inv  - usq)
                fEq(LB_NODE(8,i,j,k)) = t(8)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(1,i,j,k))-u(LB_NODE(2,i,j,k)))**2*t2cs4inv    - usq)
                fEq(LB_NODE(9,i,j,k)) = t(9)*rho(i,j,k)*(1._R8B + (u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))*cs2inv &
                      + (u(LB_NODE(1,i,j,k))-u(LB_NODE(2,i,j,k)))**2*t2cs4inv  - usq)
            end do
         end do
      end do
#endif
#ifdef D3Q19
      do k=1,s_par%gx(3)
         do j=1,s_par%gx(2)
            do i=1,s_par%gx(1)
                usq =  (u(LB_NODE(1,i,j,k))**2 + u(LB_NODE(2,i,j,k))**2 + u(LB_NODE(3,i,j,k))**2)*t2cs2inv
                 fEq(LB_NODE(1,i,j,k)) = t(1)*rho(i,j,k)*(1  &
                      - usq) 
                 fEq(LB_NODE(2,i,j,k)) = t(2)*rho(i,j,k)*(1._R8B + (u(LB_NODE(1,i,j,k)))*cs2inv &
                      + (u(LB_NODE(1,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(3,i,j,k)) = t(3)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(1,i,j,k)))*cs2inv &
                      -u(LB_NODE(1,i,j,k))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(4,i,j,k)) = t(4)*rho(i,j,k)*(1._R8B + u(LB_NODE(2,i,j,k))*cs2inv &
                      + (u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(5,i,j,k)) = t(5)*rho(i,j,k)*(1._R8B - u(LB_NODE(2,i,j,k))*cs2inv &
                      + (-u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(6,i,j,k)) = t(6)*rho(i,j,k)*(1._R8B + (u(LB_NODE(3,i,j,k)))*cs2inv &
                      + (u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(7,i,j,k)) = t(7)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(3,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(8,i,j,k)) = t(8)*rho(i,j,k)*(1._R8B + (u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))*cs2inv &
                      + (u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(9,i,j,k)) = t(9)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(10,i,j,k)) = t(10)*rho(i,j,k)*(1._R8B + (u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))*cs2inv &
                      + (u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(11,i,j,k)) = t(11)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(12,i,j,k)) = t(12)*rho(i,j,k)*(1._R8B + (u(LB_NODE(1,i,j,k)) + u(LB_NODE(3,i,j,k)))*cs2inv &
                      + (u(LB_NODE(1,i,j,k)) + u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(13,i,j,k)) = t(13)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(1,i,j,k)) -u(LB_NODE(3,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(1,i,j,k)) -u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(14,i,j,k)) = t(14)*rho(i,j,k)*(1._R8B + (u(LB_NODE(1,i,j,k)) -u(LB_NODE(3,i,j,k)))*cs2inv &
                      + (u(LB_NODE(1,i,j,k)) -u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(15,i,j,k)) = t(15)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(1,i,j,k)) + u(LB_NODE(3,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(1,i,j,k)) + u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(16,i,j,k)) = t(16)*rho(i,j,k)*(1._R8B + (u(LB_NODE(2,i,j,k)) + u(LB_NODE(3,i,j,k)))*cs2inv &
                      + (u(LB_NODE(2,i,j,k)) + u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(17,i,j,k)) = t(17)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(2,i,j,k)) -u(LB_NODE(3,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(2,i,j,k)) -u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(18,i,j,k)) = t(18)*rho(i,j,k)*(1._R8B + (u(LB_NODE(2,i,j,k)) -u(LB_NODE(3,i,j,k)))*cs2inv &
                      + (u(LB_NODE(2,i,j,k)) -u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(LB_NODE(19,i,j,k)) = t(19)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(2,i,j,k)) + u(LB_NODE(3,i,j,k)))*cs2inv &
                      + (-u(LB_NODE(2,i,j,k)) + u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
            end do
         end do
      end do
#endif    
   end subroutine calc_fEq_global
!------------------------------------------------------------------------
#endif  /*INIT_WITH_ROOT  */






   subroutine calc_fEq(lb_dom,rho,u,fEq,xmin,xmax,ymin,ymax,zmin,zmax)
   !------------------------------------------------------------------------
   !
   ! calculate the Maxwellian / equlibrium distribution function
   !




      type(lb_block), intent(inout) :: lb_dom
      real(R8B)                  :: usq,cs2inv,t2cs4inv,t2cs2inv,rholoc
      real(R8B), intent(in)      :: u(LB_NODE(NDIM,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B), intent(in)      :: rho(0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1))
      real(R8B), intent(inout)   :: fEq(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer                    :: i,j,k
      integer, intent(in)        :: xmin,xmax,ymin,ymax,zmin,zmax
      integer                    :: x_l,x_u,y_l,y_u,z_l,z_u


      !---------------------------------------------------
      ! define upper and lower limits
      ! if no limit given (0), use lower: 1 or upper: lb_dom%lx(i)
      if(xmin==0) then;  x_l=1;  else; x_l=xmin; end if
      if(xmax==0) then;  x_u=lb_dom%lx(1); else; x_u=xmax; end if
      if(ymin==0) then;  y_l=1;  else; y_l=ymin; end if
      if(ymax==0) then;  y_u=lb_dom%lx(2); else; y_u=ymax; end if
      if(zmin==0) then;  z_l=1;  else; z_l=zmin; end if
      if(zmax==0) then;  z_u=lb_dom%lx(3); else; z_u=zmax; end if
      err = 0

      !---------------------------------------------------
      ! define constants 
      cs2inv  = 1._R8B/cs**2
      t2cs4inv = 1._R8B/(2*cs**4)
      t2cs2inv = 1._R8B/(2*cs**2)

#ifdef D2Q9    
   do k=z_l,z_u
       do j=y_l,y_u
          do i=x_l,x_u
             rholoc = rho(i,j,k)
 
             usq =  (u(LB_NODE(1,i,j,k))**2 + u(LB_NODE(2,i,j,k))**2)*t2cs2inv

             fEq(LB_NODE(1,i,j,k)) = t(1)*rho(i,j,k) *(1._R8B - usq)
             fEq(LB_NODE(2,i,j,k)) = t(2)*rho(i,j,k) *(1._R8B + (u(LB_NODE(1,i,j,k)))*cs2inv &
            + (u(LB_NODE(1,i,j,k)))**2*t2cs4inv   - usq)
             fEq(LB_NODE(3,i,j,k)) = t(3)*rho(i,j,k)*(1._R8B - (-u(LB_NODE(2,i,j,k)))*cs2inv &
            + (u(LB_NODE(2,i,j,k)))**2*t2cs4inv   - usq)
             fEq(LB_NODE(4,i,j,k)) = t(4)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(1,i,j,k)))*cs2inv &
           + (-u(LB_NODE(1,i,j,k)))**2*t2cs4inv  - usq)
             fEq(LB_NODE(5,i,j,k)) = t(5)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(2,i,j,k)))*cs2inv &
           + (-u(LB_NODE(2,i,j,k)))**2*t2cs4inv   - usq)
             fEq(LB_NODE(6,i,j,k)) = t(6)*rho(i,j,k)*(1._R8B  &
            + (u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))*cs2inv &
            + (u(LB_NODE(1,i,j,k))+u(LB_NODE(2,i,j,k)))**2*t2cs4inv  - usq)
             fEq(LB_NODE(7,i,j,k)) = t(7)*rho(i,j,k)*(1._R8B  &
           + (-u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))*cs2inv &
           + (-u(LB_NODE(1,i,j,k))+u(LB_NODE(2,i,j,k)))**2*t2cs4inv  - usq)
             fEq(LB_NODE(8,i,j,k)) = t(8)*rho(i,j,k)*(1._R8B  &
           + (-u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))*cs2inv &
           + (-u(LB_NODE(1,i,j,k))-u(LB_NODE(2,i,j,k)))**2*t2cs4inv    - usq)
             fEq(LB_NODE(9,i,j,k)) = t(9)*rho(i,j,k)*(1._R8B  &
            + (u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))*cs2inv &
            + (u(LB_NODE(1,i,j,k))-u(LB_NODE(2,i,j,k)))**2*t2cs4inv  - usq)


          end do
       end do
    end do
#endif
#ifdef D3Q19
    do k=z_l,z_u
       do j=y_l,y_u
          do i=x_l,x_u
          usq =  (u(LB_NODE(1,i,j,k))**2 + u(LB_NODE(2,i,j,k))**2 + u(LB_NODE(3,i,j,k))**2)*t2cs2inv
          rholoc   = rho(i,j,k) 
           fEq(LB_NODE(1,i,j,k)) = t(1)*rho(i,j,k)*(1  &
                - usq) 
           fEq(LB_NODE(2,i,j,k)) = t(2)*rho(i,j,k)*(1._R8B + (u(LB_NODE(1,i,j,k)))*cs2inv &
                + (u(LB_NODE(1,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(3,i,j,k)) = t(3)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(1,i,j,k)))*cs2inv &
                -u(LB_NODE(1,i,j,k))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(4,i,j,k)) = t(4)*rho(i,j,k)*(1._R8B + u(LB_NODE(2,i,j,k))*cs2inv &
                + (u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(5,i,j,k)) = t(5)*rho(i,j,k)*(1._R8B - u(LB_NODE(2,i,j,k))*cs2inv &
                + (-u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(6,i,j,k)) = t(6)*rho(i,j,k)*(1._R8B + (u(LB_NODE(3,i,j,k)))*cs2inv &
                + (u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(7,i,j,k)) = t(7)*rho(i,j,k)*(1._R8B + (-u(LB_NODE(3,i,j,k)))*cs2inv &
                + (-u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(8,i,j,k)) = t(8)*rho(i,j,k)*(1._R8B  &
                + (u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))*cs2inv &
                + (u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(9,i,j,k)) = t(9)*rho(i,j,k)*(1._R8B  &
                + (-u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))*cs2inv &
                + (-u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(10,i,j,k)) = t(10)*rho(i,j,k)*(1._R8B  &
                + (u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))*cs2inv &
                + (u(LB_NODE(1,i,j,k)) -u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(11,i,j,k)) = t(11)*rho(i,j,k)*(1._R8B  &
                + (-u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))*cs2inv &
                + (-u(LB_NODE(1,i,j,k)) + u(LB_NODE(2,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(12,i,j,k)) = t(12)*rho(i,j,k)*(1._R8B  &
                + (u(LB_NODE(1,i,j,k)) + u(LB_NODE(3,i,j,k)))*cs2inv &
                + (u(LB_NODE(1,i,j,k)) + u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(13,i,j,k)) = t(13)*rho(i,j,k)*(1._R8B  &
                + (-u(LB_NODE(1,i,j,k)) -u(LB_NODE(3,i,j,k)))*cs2inv &
                + (-u(LB_NODE(1,i,j,k)) -u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(14,i,j,k)) = t(14)*rho(i,j,k)*(1._R8B  &
                + (u(LB_NODE(1,i,j,k)) -u(LB_NODE(3,i,j,k)))*cs2inv &
                + (u(LB_NODE(1,i,j,k)) -u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(15,i,j,k)) = t(15)*rho(i,j,k)*(1._R8B  &
                + (-u(LB_NODE(1,i,j,k)) + u(LB_NODE(3,i,j,k)))*cs2inv &
                + (-u(LB_NODE(1,i,j,k)) + u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(16,i,j,k)) = t(16)*rho(i,j,k)*(1._R8B  &
                + (u(LB_NODE(2,i,j,k)) + u(LB_NODE(3,i,j,k)))*cs2inv &
                + (u(LB_NODE(2,i,j,k)) + u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(17,i,j,k)) = t(17)*rho(i,j,k)*(1._R8B  &
                + (-u(LB_NODE(2,i,j,k)) -u(LB_NODE(3,i,j,k)))*cs2inv &
                + (-u(LB_NODE(2,i,j,k)) -u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(18,i,j,k)) = t(18)*rho(i,j,k)*(1._R8B  &
                + (u(LB_NODE(2,i,j,k)) -u(LB_NODE(3,i,j,k)))*cs2inv &
                + (u(LB_NODE(2,i,j,k)) -u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                - usq) 
           fEq(LB_NODE(19,i,j,k)) = t(19)*rho(i,j,k)*(1._R8B  &
                + (-u(LB_NODE(2,i,j,k)) + u(LB_NODE(3,i,j,k)))*cs2inv &
                + (-u(LB_NODE(2,i,j,k)) + u(LB_NODE(3,i,j,k)))**2*t2cs4inv  &
                - usq) 
               
          end do
       end do
    end do
#endif    
  end subroutine calc_fEq
  !------------------------------------------------------------------------







  subroutine get_geo_each(lb_dom,s_par,prc)
    !------------------------------------------------------------------------
    !
    ! set up the geometry, each process by itself
    !
      implicit none
      type(lb_block),intent(inout)      :: lb_dom 
      type(sim_parameter),intent(inout) :: s_par
      type(mpl_var),intent(inout)       :: prc  
      integer        :: i,j,k,l
      integer        :: x_start,x_end,y_start,y_end,z_start,z_end,lx,ly,lz
      real   (R8B)   :: x_pos, xs,xe
      integer        :: ys,ye 
      ! Fill up the Coordinates for the Grid
      lx = lb_dom%lx(1)
      ly = lb_dom%lx(2)
      lz = lb_dom%lx(3)

      if(prc%crd(1) == 0)          then; x_start = 2;    else; x_start=1;  end if
      if(prc%crd(1) == prc%np(1)-1) then; x_end   = lx-1; else; x_end  =lx; end if
      if(prc%crd(2) == 0)          then; y_start = 2;    else; y_start=1;  end if
      if(prc%crd(2) == prc%np(2)-1) then; y_end   = ly-1; else; y_end  =ly; end if
#ifdef D3Q19
      if(prc%crd(3) == 0)          then; z_start = 2;    else; z_start=1;  end if
      if(prc%crd(3) == prc%np(3)-1) then; z_end   = lz-1; else; z_end  =lz; end if
#endif
#ifdef D2Q9
      z_start = 1
      z_end   = 1
#endif

      do i=1,prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1,2)-prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1,1)+1
         lb_dom%x(i)  = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,1,1)+real(i,8) - 2
      end do
      do j=1,prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,2,2)-prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,2,1)+1
         lb_dom%y(j)  = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,2,1)+real(j,8) - 2
      end do
      do k=1,prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,3,2)-prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,3,1)+1
         lb_dom%z(k)  = prc%bnd(prc%crd(1)+1,prc%crd(2)+1,prc%crd(3)+1,3,1)+real(k,8) - 2
      end do    

      ! Set all nodes to fluid
      lb_dom%state(:,:,:)   = 0 
      lb_dom%state = ibset(lb_dom%state,fluid) 


#ifdef SPONGE
   write(*,*) "calculation of omega for sponge layer not implemented yet."
   stop
#endif



!--------------------------------------------------------------------
! Select defined problem  


      select case(s_par%problem)



!----------------------------------------
! Cylinder in channel

       case(cylinder)
         ! wall for z=0 and z=lz has to be implemented
         do j=y_start,y_end

            if(prc%crd(1)==0) then
            ! inlet left  
!            lb_dom%state(1,j,:)    = nr_wall ! inlet  !
            lb_dom%state(1,j,:)    = ibset(lb_dom%state(1,j,:),nr_wall)
!            do kk=1,size(lb_dom%state,3)
!            lb_dom%state(1,j,:)    = ibset(lb_dom%state(1,j,:),nr_wall)
!            enddo
            endif
            if(prc%crd(1)==prc%np(1)-1) then
            ! outlet right  
!            lb_dom%state(lx,j,:)   = nr_wall ! outlet !
            lb_dom%state(lx,j,:)   = ibset(lb_dom%state(lx,j,:),nr_wall) 
            endif
         end do
         ! wall boundaries
         if(prc%crd(2)==0)           &
            lb_dom%state(x_start:x_end,1,:)  = ibset(lb_dom%state(x_start:x_end,1,:), wall)
         if(prc%crd(2)==prc%np(1)-1) &
            lb_dom%state(x_start:x_end,ly,:) = ibset(lb_dom%state(x_start:x_end,ly,:),wall)

#ifdef D3Q19
         ! wall boundaries
         if(prc%crd(3) ==0)            lb_dom%state(:,:,1)  = ibset(lb_dom%state(:,:,1), wall)
         if(prc%crd(3) ==prc%np(3)-1)  lb_dom%state(:,:,lz) = ibset(lb_dom%state(:,:,lz),wall)
#endif
         ! cylinder
         write(*,*) "todo: implement cylinder"


      case(flute   )
         lb_dom%state(:,lb_dom%lx(2),:)      = ibset(lb_dom%state(:,lb_dom%lx(2),:),outlet)
         lb_dom%state(1:lb_dom%lx(1)/5,:,:)  = ibset(lb_dom%state(1:lb_dom%lx(1)/5,:,:),wall)
         lb_dom%state(lb_dom%lx(1)/10:lb_dom%lx(1)/5, lb_dom%lx(2)/4+8:lb_dom%lx(2),:) = & 
         ibset(lb_dom%state(lb_dom%lx(1)/10:lb_dom%lx(1)/5, lb_dom%lx(2)/4+8:lb_dom%lx(2),:),fluid)
         lb_dom%state(lb_dom%lx(1),:,:)      = ibset(lb_dom%state(lb_dom%lx(1),:,:),wall)
      do j=lb_dom%lx(2)/4+2,lb_dom%lx(2)
          ! outlet rechts
          lb_dom%state(lb_dom%lx(1),j,:)     = ibset(lb_dom%state(lb_dom%lx(1),j,:),nr_wall)
       end do

       do j=lb_dom%lx(2)/4-4,lb_dom%lx(2)/4+4
          ! inlet links 
          lb_dom%state(1,j,:)                 = ibset(lb_dom%state(1,j,:),inlet)
          lb_dom%state(2:lb_dom%lx(1)/5 ,j,:) = ibset(lb_dom%state(2:lb_dom%lx(1)/5 ,j,:),fluid )
       end do
       ! wall boundaries
       lb_dom%state(1:lb_dom%lx(1),1,:)        = ibset(lb_dom%state(1:lb_dom%lx(1),1,:),wall)
       lb_dom%state(1:lb_dom%lx(1),lb_dom%lx(2),:)  = ibset(lb_dom%state(1:lb_dom%lx(1),lb_dom%lx(2),:),wall)

       ! flute labium 
       xs = lb_dom%lx(1)/10*2.4d0
       xe = lb_dom%lx(1)/10*3.1
       ys = lb_dom%lx(2)/4-2+2
       ye = lb_dom%lx(2)/4+8-1
       do j=int(ys),int(ye)
         x_pos = xs + (xe-xs)/(ye-ys)*(j-ys) 
         lb_dom%state(int(x_pos):lb_dom%lx(1),j,1) = ibset(lb_dom%state(int(x_pos):lb_dom%lx(1),j,1),wall)
       end do


      case(cavity)
         ! where is the lid? y has to be ly, x and z have to be 2 or 1, depending on where the prc%crd is
         if(prc%crd(2) == prc%np(2) - 1) then 
             lb_dom%state(  x_start:x_end,  ly,  z_start:z_end  )  = &
         ibset(lb_dom%state(  x_start:x_end,  ly,  z_start:z_end  ),lid)
         end if
         ! left wall (at y- border)
!FIXME
         if(prc%crd(1) == 0)           lb_dom%state(1,:,:)     = ibset(lb_dom%state(1,:,:),wall)
         if(prc%crd(1) == prc%np(1)-1) lb_dom%state(lx,:,:)    = ibset(lb_dom%state(lx,:,:),wall)
         ! bottom wall
         if(prc%crd(2) == 0)           lb_dom%state(:,1,:)     = ibset(lb_dom%state(:,1,:) ,wall)
         if(NDIM==3) then
           ! set also front and back wall to wall condition
            if(prc%crd(3) == 0)        lb_dom%state(:,:,1)     = ibset(lb_dom%state(:,:,1),wall)
            if(prc%crd(3) == prc%np(3)-1) lb_dom%state(:,:,lz) = ibset(lb_dom%state(:,:,lz),wall )        
         end if

!----------------------------------------
! Same Lid-driven cavity on all partitions

      case(cavity_same)
         lb_dom%state(2:(lx-1),ly,:)  = ibset(lb_dom%state(2:(lx-1),ly,:),lid)
         lb_dom%state(1:lx,1,:)       = ibset(lb_dom%state(1:lx,1,:),wall)
         lb_dom%state(1,1:ly,:)       = ibset(lb_dom%state(1,1:ly,:),wall)
         lb_dom%state(lx,1:ly,:)      = ibset(lb_dom%state(lx,1:ly,:),wall)
         if(NDIM==3) then
            ! set also front and back wall to wall condition
            lb_dom%state(:,:,1)       = ibset(lb_dom%state(:,:,1),wall)
            lb_dom%state(:,:,lz)      = ibset(lb_dom%state(:,:,lz),wall)
        end if


   case(shock)
         write(*,*) "not implemented yet"
         stop
   case(gaussian)
         ! left wall (at y- border)
         if(prc%crd(1) == 0)              lb_dom%state(1,:,:)     = ibset(lb_dom%state(1,:,:), nr_wall)
         if(prc%crd(1) == prc%np(1)-1)    lb_dom%state(lx,:,:)    = ibset(lb_dom%state(lx,:,:),nr_wall)
         ! bottom wall
         if(prc%crd(2) == 0)              lb_dom%state(:,1,:)     = ibset(lb_dom%state(:,1,:), nr_wall)
         if(prc%crd(2) == prc%np(2)-1)    lb_dom%state(:,ly,:)    = ibset(lb_dom%state(:,ly,:),nr_wall)
         if(NDIM==3) then
           ! set also front and back wall to wall condition
            if(prc%crd(3) == 0)           lb_dom%state(:,:,1)     = ibset(lb_dom%state(:,:,1), nr_wall)
            if(prc%crd(3) == prc%np(3)-1) lb_dom%state(:,:,lz)    = ibset(lb_dom%state(:,:,lz),nr_wall)        
         end if
   end select

! get number of inlet nodes in directions
#ifdef D2Q9
lb_dom%num_in(2) = 0
do j=1,lb_dom%lx(2)
   if(lb_dom%state(1,j,1) == inlet) then
      lb_dom%num_in(2) = lb_dom%num_in(2) + 1
   endif
end do
#endif
   end subroutine get_geo_each
   !------------------------------------------------------------------------

  







  subroutine get_geo(lb_dom,s_par,prc)


    !------------------------------------------------------------------------
    !
    ! set up the geometry
    !


      implicit none
      type(lb_block),intent(inout)  :: lb_dom 
      type(sim_parameter),intent(inout)  :: s_par
      type(mpl_var),intent(in)     :: prc  
      integer               :: i,j,k,l,lx,ly,lz
      real(R8B)             :: xs,xe,ys,ye,x_pos 
      integer               :: in_xmax,right_wall_min,posbc 
#ifdef SPONGE
      real(R8B)             :: distance,dist_max,outer_omega
      integer               :: sponge_radius,lchar
#endif /* SPONGE */
      ! Fill up the Coordinates for the Grid
      lx=s_par%gx(1)
      ly=s_par%gx(2)
      lz=s_par%gx(3)

      do i=1,s_par%gx(1)
         s_par%g_x(i)  = real(i,8) - 1.
      end do
      do j=1,s_par%gx(2)
         s_par%g_y(j)  = real(j,8) - 1.
      end do
      do k=1,s_par%gx(3)
         s_par%g_z(k)  = real(k,8) - 1.
      end do    

      ! Reset 
      lb_dom%gstate = 0
      ! Set all nodes to fluid
      lb_dom%gstate = ibset(lb_dom%gstate,fluid) 

#ifdef SPONGE
!--------------------------------------------------------------------
! Spongelayer: assign viscosity to domain 


         ! set the radius, from where the sponge layer starts, starting from the middle of the domain
         if(NDIM==2) then
            lchar = min(lx,ly)
         else
            lchar = min(lx,ly,lz)
         endif
         sponge_radius =int(real(lchar)/5.d0)
         ! set the minimal omega
         outer_omega = 0.1 

         ! assign the viscosity to each grid point
         do k=1,lz
            do j=1,ly
               do i=1,lx

                  distance = sqrt((real(i)-real(lx)/2.)**2 + (real(j)-real(ly)/2.)**2 + (real(k)-real(lz)/2.)**2 ) 
                  dist_max = real(lchar)/2. 

                  if(distance .ge. dist_max) then
                     lb_dom%gomega(i,j,k) =outer_omega 
                  elseif(distance .ge. sponge_radius .and. distance .lt.dist_max) then
                     lb_dom%gomega(i,j,k) = s_par%omega + (outer_omega - s_par%omega)  &
                                       /(dist_max-sponge_radius)*(distance-sponge_radius)
                  else 
                     lb_dom%gomega(i,j,k) = s_par%omega
                  endif
               enddo
            enddo
         enddo
#endif /* SPONGE */



!--------------------------------------------------------------------
! Select defined problem  


      select case(s_par%problem)
         


!----------------------------------------
! Single Taylor Vortex decay

   case(taylor_vortex)

         posbc = nr_wall
         posbc = per_wall 

         ! use non-reflecting boundary
         lb_dom%gstate(1:lx,ly,:)  = ibset(lb_dom%gstate(1:lx,ly,:),posbc) 
         lb_dom%gstate(1:lx,1,:)   = ibset(lb_dom%gstate(1:lx,1,:), posbc)
         lb_dom%gstate(1,1:ly,:)   = ibset(lb_dom%gstate(1,1:ly,:), posbc)
         lb_dom%gstate(lx,1:ly,:)  = ibset(lb_dom%gstate(lx,1:ly,:),posbc)
         if(NDIM==3) then
            ! set also front and back wall to wall condition
            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), posbc)
            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),posbc)
          endif


!----------------------------------------
! Single Taylor Vortex decay

   case(planar_standing_wave)

         posbc = nr_wall
         posbc = per_wall 

         ! use non-reflecting boundary
         lb_dom%gstate(1:lx,ly,:)  = ibset(lb_dom%gstate(1:lx,ly,:),posbc) 
         lb_dom%gstate(1:lx,1,:)   = ibset(lb_dom%gstate(1:lx,1,:), posbc)
         lb_dom%gstate(1,1:ly,:)   = ibset(lb_dom%gstate(1,1:ly,:), posbc)
         lb_dom%gstate(lx,1:ly,:)  = ibset(lb_dom%gstate(lx,1:ly,:),posbc)
         if(NDIM==3) then
            ! set also front and back wall to wall condition
            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), posbc)
            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),posbc)
          endif




!----------------------------------------
! Cylinder in channel

       case(cylinder)
         ! wall for z=0 and z=lz has to be implemented
         do j=2,ly-1

            ! inlet left  
!            lb_dom%gstate(1,j,:)    = nr_wall ! inlet  !
            lb_dom%gstate(1,j,:)    = ibset(lb_dom%gstate(1,j,:),nr_wall)
         
            ! outlet right  
!           lb_dom%gstate(lx,j,:)   = nr_wall ! outlet !
            lb_dom%gstate(lx,j,:)   = ibset(lb_dom%gstate(lx,j,:),nr_wall) 


         end do
         ! wall boundaries
!         lb_dom%gstate(1:lx,1,:)    = wall
!         lb_dom%gstate(1:lx,ly,:)   = wall
         lb_dom%gstate(1:lx,1,:)    = ibset(lb_dom%gstate(1:lx,1,:), wall)
         lb_dom%gstate(1:lx,ly,:)   = ibset(lb_dom%gstate(1:lx,ly,:),wall)

#ifdef D3Q19
         ! wall boundaries
!         lb_dom%gstate(:,:,1)    = wall
!         lb_dom%gstate(:,:,lz)   = wall
         lb_dom%gstate(:,:,1)    = ibset(lb_dom%gstate(:,:,1),wall)
         lb_dom%gstate(:,:,lz)   = ibset(lb_dom%gstate(:,:,lz),wall)
#endif
         ! cylinder
         do i=1,lx
            do j=1,ly
               if((real(i)-real(s_par%obst_x))**2. + (real(j)-real(s_par%obst_y))**2. <= real(s_par%obst_r)**2.) then
                  lb_dom%gstate(i,j,:) = ibset(lb_dom%gstate(i,j,:),wall)
               end if
            end do
         end do



!----------------------------------------
! Flute

      case(flute   )
         right_wall_min = s_par%obst_y + s_par%obst_r
         in_xmax=s_par%obst_x
         ! upper outlet 
         lb_dom%gstate(:,s_par%gx(2),:)  = ibset(lb_dom%gstate(:,s_par%gx(2),:),nr_wall)
         ! left solid block
       do i=1,in_xmax
         lb_dom%gstate(i,:,:)    = ibset(lb_dom%gstate(i,:,:),wall)
       end do
 
       do j=right_wall_min,s_par%gx(2)
          ! outlet rechts
          lb_dom%gstate(s_par%gx(1),j,:) = ibset(lb_dom%gstate(s_par%gx(1),j,:),nr_wall)
       end do
       do j=1,right_wall_min
       lb_dom%gstate(s_par%gx(1),j,:)     = ibset(lb_dom%gstate(s_par%gx(1),j,:),wall)
       end do

       do j=s_par%obst_y,s_par%obst_y+2*s_par%obst_r
          ! inlet links 
          lb_dom%gstate(1,j,:)           = ibset(lb_dom%gstate(1,j,:),nr_wall)
          lb_dom%gstate(1,j,:)           = ibclr(lb_dom%gstate(1,j,:),wall)
          do i=2,in_xmax
             lb_dom%gstate(i ,j,:)  = ibclr(lb_dom%gstate(i ,j,:),wall)
          end do
       end do

       ! lower and upper wall boundaries
       lb_dom%gstate(1:s_par%gx(1),1,:)           = ibset(lb_dom%gstate(1:s_par%gx(1),1,:),wall)
       lb_dom%gstate(1:s_par%gx(1),s_par%gx(2),:) = ibset(lb_dom%gstate(1:s_par%gx(1),s_par%gx(2),:),wall)

       ! flute labium 
       xs = s_par%gx(1)/10*2.4d0
       xe = s_par%gx(1)/10*3.1
       ys = s_par%gx(2)/4
       ye = s_par%gx(2)/4+s_par%gx(2)/50

       do j=int(ys),int(ye)
         x_pos = xs + (xe-xs)/(ye-ys)*(j-ys) 
         lb_dom%gstate(int(x_pos):s_par%gx(1),j,1) = ibset(lb_dom%gstate(int(x_pos):s_par%gx(1),j,1),wall)
       end do



!----------------------------------------
! CAA Cavity      

      case(caa_cavity)
      if(NDIM==3) then
         write(*,*) "not implemented"
         stop
      endif
      
      ! flat plate      
      lb_dom%gstate(1:s_par%gx(1)/2,s_par%gx(2)/2-1:s_par%gx(2)/2+1,:) = &
ibset(lb_dom%gstate(1:s_par%gx(1)/2,s_par%gx(2)/2-1:s_par%gx(2)/2+1,:),wall)

      ! left, right borders (nrbc)
      lb_dom%gstate(1,:,:) = ibset(lb_dom%gstate(1,:,:),nr_wall)
      lb_dom%gstate(s_par%gx(1),:,:) = ibset(lb_dom%gstate(s_par%gx(1),:,:),nr_wall)

      ! top, bottom borders (nrbc)
      lb_dom%gstate(:,1,:) = ibset(lb_dom%gstate(:,1,:),nr_wall)
      lb_dom%gstate(:,s_par%gx(2),:) = ibset(lb_dom%gstate(:,s_par%gx(2),:),nr_wall)

!----------------------------------------
! CAA Cavity      

       case(plate_trail)
         


!----------------------------------------
! Lid-driven cavity

       case(cavity)
         lb_dom%gstate(2:(lx-1),ly,:)  = ibset(lb_dom%gstate(2:(lx-1),ly,:),lid)
         lb_dom%gstate(1:lx,1,:)       = ibset(lb_dom%gstate(1:lx,1,:),wall)
         lb_dom%gstate(1,1:ly,:)       = ibset(lb_dom%gstate(1,1:ly,:),wall)
         lb_dom%gstate(lx,1:ly,:)      = ibset(lb_dom%gstate(lx,1:ly,:),wall)
         if(NDIM==3) then
            ! set also front and back wall to wall condition
            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1),wall)
            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),wall)         ! insert some obstacles
            if(s_par%setObstacles .EQV. .TRUE.) then
            do k = 4,lz-4,4
               do j = 4,ly-4,4
                  do i = 4,lx-4,4
                     lb_dom%gstate(i,j,k) = ibset(lb_dom%gstate(i,j,k),wall)
                  end do 
               end do 
            end do 
           end if
        end if




!----------------------------------------
! Density shock 

   case(shock)
        lb_dom%gstate(2:(lx-1),ly,:)  = ibset(lb_dom%gstate(2:(lx-1),ly,:),wall) 
        lb_dom%gstate(1:lx,1,:)       = ibset(lb_dom%gstate(1:lx,1,:), wall)
        lb_dom%gstate(1,1:ly,:)       = ibset(lb_dom%gstate(1,1:ly,:), wall)
        lb_dom%gstate(lx,1:ly,:)      = ibset(lb_dom%gstate(lx,1:ly,:),wall)
        if(NDIM==3) then
           ! set also front and back wall to wall condition
           lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), wall)
           lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),wall)
           ! insert some obstacles
           if(s_par%setObstacles .EQV. .TRUE.) then
            do k = 4,lz-4,4
               do j = 4,ly-4,4
                  do i = 4,lx-4,4
                     lb_dom%gstate(i,j,k) = ibset(lb_dom%gstate(i,j,k),wall)
                  end do 
               end do 
            end do 
           end if
        end if




!----------------------------------------
! Gaussian Pulse (one source) 


   case(gaussian)
         posbc = per_wall
         posbc = nr_wall
         ! use non-reflecting boundary
         lb_dom%gstate(1,1:ly,:)   = ibset(lb_dom%gstate(1,1:ly,:), posbc) 
         lb_dom%gstate(lx,1:ly,:)  = ibset(lb_dom%gstate(lx,1:ly,:),posbc) 
         lb_dom%gstate(1:lx,ly,:)  = ibset(lb_dom%gstate(1:lx,ly,:),posbc) 
         lb_dom%gstate(1:lx,1,:)   = ibset(lb_dom%gstate(1:lx,1,:), posbc) 
          if(NDIM==3) then
            ! set also front and back wall to wall condition
            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), posbc)
            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),posbc)
          endif

   case(gauss_convect)
         ! use non-reflecting boundary
         posbc = nr_wall
        posbc = per_wall 
         lb_dom%gstate(1,1:ly,:)   = ibset(lb_dom%gstate(1,1:ly,:), posbc) 
         ! if inlet, then set separately
         posbc = nr_wall
        posbc = per_wall 
         lb_dom%gstate(1,1:ly,:)   = ibset(lb_dom%gstate(1,1:ly,:), posbc) 
         posbc = nr_wall
        posbc = per_wall 
         lb_dom%gstate(lx,1:ly,:)  = ibset(lb_dom%gstate(lx,1:ly,:),posbc) 
         posbc = nr_wall 
        posbc = per_wall 
         lb_dom%gstate(1:lx,ly,:)  = ibset(lb_dom%gstate(1:lx,ly,:),posbc) 
         lb_dom%gstate(1:lx,1,:)   = ibset(lb_dom%gstate(1:lx,1,:), posbc) 
          if(NDIM==3) then
            ! set also front and back wall to wall condition
            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), posbc)
            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),posbc)
          endif




!----------------------------------------
! Corotating Vortex with Spongelayer

   case(corotating_vortex)
#ifdef SPONGE
         posbc = nr_wall
         ! use sponge layer and reset wall
         lb_dom%gstate(1:lx,ly,:)  = ibset(lb_dom%gstate(1:lx,ly,:),posbc) 
         lb_dom%gstate(1:lx,1,:)   = ibset(lb_dom%gstate(1:lx,1,:), posbc)
         lb_dom%gstate(1,1:ly,:)   = ibset(lb_dom%gstate(1,1:ly,:), posbc)
         lb_dom%gstate(lx,1:ly,:)  = ibset(lb_dom%gstate(lx,1:ly,:),posbc)
         if(NDIM==3) then
           ! set also front and back wall to wall condition
            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), posbc)
            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),posbc)
         endif
         ! assign the viscosity
         
         ! set the radius, from where the sponge layer starts, starting from the middle of the domain
         if(NDIM==2) then
            lchar = min(lx,ly)
         else
            lchar = min(lx,ly,lz)
         endif
         sponge_radius =int(real(lchar)/3.d0)

         outer_omega = 0.1 
         ! the max distance is half the domain diameter minus sponge radius
         dist_max = real(lchar)/2.d0-sponge_radius
 !sqrt((real(lx/2))**2 + (real(ly/2))**2 + (real(lz/2))**2)-sponge_radius

         ! assign the viscosity to each grid point
         do k=1,lz
            do j=1,ly
               do i=1,lx
                  ! distance from domain middle
                  distance = sqrt((real(i)-real(lx)/2.)**2 + (real(j)-real(ly)/2.)**2 + (real(k)-real(lz)/2.)**2 ) 
                  if(distance .ge. sponge_radius .and. distance .lt. dist_max+sponge_radius) then
                     ! if inside the sponge area, that means outside the inner circle defined by sponge_radius
                     lb_dom%gomega(i,j,k) = s_par%omega + (outer_omega - s_par%omega)/(dist_max)*(distance-sponge_radius)
                  elseif(distance .ge. dist_max+sponge_radius) then
                     lb_dom%gomega(i,j,k) = outer_omega
                  else 
                     lb_dom%gomega(i,j,k) = s_par%omega
                  endif
               enddo
            enddo
         enddo
#else /* not SPONGE -> use nrbc*/


!----------------------------------------
! Corotating Vortex without Sponge Layer
         posbc = per_wall
!FIXME careful, using periodic boundaries, not nrbc

         do i=1,2 !3   If until 3 -> 
            if(i==1)  posbc = nr_wall
            if(i==2)  posbc = wall
            if(i==3)  posbc = nr_wall_in
         ! use non-reflecting boundary
         lb_dom%gstate(1:lx,ly,:)  = ibset(lb_dom%gstate(1:lx,ly,:),posbc) 
         lb_dom%gstate(1:lx,1,:)   = ibset(lb_dom%gstate(1:lx,1,:), posbc)
         lb_dom%gstate(1,1:ly,:)   = ibset(lb_dom%gstate(1,1:ly,:), posbc)
         lb_dom%gstate(lx,1:ly,:)  = ibset(lb_dom%gstate(lx,1:ly,:),posbc)
         if(NDIM==3) then
            ! set also front and back wall to wall condition
            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), posbc)
            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),posbc)
          endif
         enddo

#endif





!----------------------------------------
! Cylinder Impulse

   case(cylinder_impulse)
        lb_dom%gstate(2:(lx-1),ly,:)  = ibset(lb_dom%gstate(2:(lx-1),ly,:),wall) 
        lb_dom%gstate(1:lx,1,:)       = ibset(lb_dom%gstate(1:lx,1,:),wall)
        lb_dom%gstate(1,1:ly,:)       = ibset(lb_dom%gstate(1,1:ly,:),wall)
        lb_dom%gstate(lx,1:ly,:)      = ibset(lb_dom%gstate(lx,1:ly,:),wall)
        if(NDIM==3) then
           ! set also front and back wall to wall condition
           lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1) ,wall)
           lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),wall)         ! insert some obstacles
        end if
        ! set cylinder walls
        lb_dom%gstate(1:(lx/2),ly/2+6,:)  = ibset(lb_dom%gstate(1:(lx/2),ly/2+6,:),wall) 
        lb_dom%gstate(1:(lx/2),ly/2-6,:)  = ibset(lb_dom%gstate(1:(lx/2),ly/2-6,:),wall) 
        

   end select

   end subroutine get_geo
   !------------------------------------------------------------------------






   subroutine calc_macr_vals(lb_dom,u,rho,fIn,meas)
   !------------------------------------------------------------------------
   !
   ! calculate the macroscopic values from the microscopic distributions 
   !
     implicit none
     type(lb_block),intent(inout)  :: lb_dom
     type(measure )                :: meas   
     real(R8B), intent(out)    :: u(LB_NODE(NDIM,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
     real(R8B), intent(out)    :: rho(0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1))
     real(R8B), intent(in)     :: fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
     integer              :: i,j,k,l
     
         call cpu_time_measure(meas%tSt_comm) 

     do k=1,lb_dom%lx(3)
        do j=1,lb_dom%lx(2)
           do i=1,lb_dom%lx(1)
               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
!if(lb_dom%state(i,j,k) /= wall) then
#ifdef D2Q9
               rho(i,j,k) = fIn(LB_NODE(1,i,j,k)) + fIn(LB_NODE(2,i,j,k)) + & 
                            fIn(LB_NODE(3,i,j,k)) + fIn(LB_NODE(4,i,j,k)) + &
                            fIn(LB_NODE(5,i,j,k)) + fIn(LB_NODE(6,i,j,k)) + &
                            fIn(LB_NODE(7,i,j,k)) + fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(9,i,j,k)) 
               u(LB_NODE(1,i,j,k)) = 1/rho(i,j,k) * (fIn(LB_NODE(2,i,j,k))  - &
                            fIn(LB_NODE(4,i,j,k)) + &
                            fIn(LB_NODE(6,i,j,k)) - fIn(LB_NODE(7,i,j,k)) - &
                            fIn(LB_NODE(8,i,j,k)) + fIn(LB_NODE(9,i,j,k)))
               u(LB_NODE(2,i,j,k)) = 1/rho(i,j,k) * (fIn(LB_NODE(3,i,j,k))  - fIn(LB_NODE(5,i,j,k)) + &
                            fIn(LB_NODE(6,i,j,k)) + fIn(LB_NODE(7,i,j,k)) - &
                            fIn(LB_NODE(8,i,j,k)) - fIn(LB_NODE(9,i,j,k)))
#endif
#ifdef D3Q19
               rho(i,j,k) = fIn(LB_NODE(1,i,j,k))  + fIn(LB_NODE(2,i,j,k))  + fIn(LB_NODE(3,i,j,k))  + &
                            fIn(LB_NODE(4,i,j,k))  + fIn(LB_NODE(5,i,j,k)) + &
                            fIn(LB_NODE(6,i,j,k))  + fIn(LB_NODE(7,i,j,k))  + &
                            fIn(LB_NODE(8,i,j,k))  + fIn(LB_NODE(9,i,j,k))  + fIn(LB_NODE(10,i,j,k)) + & 
                            fIn(LB_NODE(11,i,j,k)) + fIn(LB_NODE(12,i,j,k)) + &
                            fIn(LB_NODE(13,i,j,k)) + fIn(LB_NODE(14,i,j,k)) + fIn(LB_NODE(15,i,j,k)) + & 
                            fIn(LB_NODE(16,i,j,k)) + fIn(LB_NODE(17,i,j,k)) + &
                            fIn(LB_NODE(18,i,j,k)) + fIn(LB_NODE(19,i,j,k))
               u(LB_NODE(1,i,j,k)) = (fIn(LB_NODE(2,i,j,k))  - fIn(LB_NODE(3,i,j,k))  + &
                            fIn(LB_NODE(8,i,j,k))  - fIn(LB_NODE(9,i,j,k))  + &
                            fIn(LB_NODE(10,i,j,k)) - & 
                            fIn(LB_NODE(11,i,j,k)) + fIn(LB_NODE(12,i,j,k)) - &
                            fIn(LB_NODE(13,i,j,k)) + fIn(LB_NODE(14,i,j,k)) - fIn(LB_NODE(15,i,j,k))  & 
                            )/rho(i,j,k)
               u(LB_NODE(2,i,j,k)) = (fIn(LB_NODE(4,i,j,k))  - fIn(LB_NODE(5,i,j,k)) + &
                            fIn(LB_NODE(8,i,j,k))  - fIn(LB_NODE(9,i,j,k))  - &
                            fIn(LB_NODE(10,i,j,k)) + & 
                            fIn(LB_NODE(11,i,j,k)) + & 
                            fIn(LB_NODE(16,i,j,k)) - fIn(LB_NODE(17,i,j,k)) + &
                            fIn(LB_NODE(18,i,j,k)) - fIn(LB_NODE(19,i,j,k)))/rho(i,j,k)
               u(LB_NODE(3,i,j,k)) = (fIn(LB_NODE(6,i,j,k))  - fIn(LB_NODE(7,i,j,k))  + &
                            fIn(LB_NODE(12,i,j,k)) - fIn(LB_NODE(13,i,j,k) )- fIn(LB_NODE(14,i,j,k))  + & 
                            fIn(LB_NODE(15,i,j,k)) + fIn(LB_NODE(16,i,j,k)) - &
                            fIn(LB_NODE(17,i,j,k)) - fIn(LB_NODE(18,i,j,k)) + fIn(LB_NODE(19,i,j,k)))/rho(i,j,k)
#endif
endif

            enddo
         enddo
      enddo
         call cpu_time_measure(meas%tEnd_comm)
         meas%ccmv_duration = meas%tEnd_comm - meas%tSt_comm + meas%ccmv_duration
   end subroutine calc_macr_vals
   !------------------------------------------------------------------------









!   subroutine set_reset_wall(lb_dom,state,fIn,rho,u,s_par,prc,meas)
   !------------------------------------------------------------------------
   !
   ! set the boundary conditions in each timestep.
   ! equilibrium distribution is assumed at the boundaries
   !
!      implicit none
!      type(lb_block)       :: lb_dom
!      type(sim_parameter)  :: s_par 
!      type(mpl_var)        :: prc
!      type(measure)        :: meas
!      integer              :: state(lb_dom%lx(1),lb_dom%lx(2),lb_dom%lx(3))
!      integer              :: i,j,k,l
!      integer              :: x,y,z,x1,y1,z1
!      real(R8B),dimension(LB_NODE(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
!      real(R8B),dimension(LB_NODE(1:NDIM,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: u
!      real(R8B),dimension(0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1)    :: rho
!
!      do  i=1,lb_dom%nobs_reset
!         x = lb_dom%obs_reset(i,1)
!         y = lb_dom%obs_reset(i,2)
!         z = lb_dom%obs_reset(i,3)
!         rho(x,y,z) = 1.0d0
!         u(LB_NODE(:,x,y,z)) = 0.0d0
!      ! reset boundary nodes to equilibrium 
!      !FIXME repair the coordinates (use x_start, x_end etc.)
!         call calc_fEq(lb_dom,rho,u,fIn,x,x,y ,y ,z,z) ! top wall
!         lb_dom%fIn(LB_NODE(:,x,y,z)) = lb_dom%fIn(LB_NODE(:,x1,y1,z1))
!      end do
!   end subroutine set_reset_wall
!   !------------------------------------------------------------------------





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
      integer              :: ii
      real(R8B),dimension(LB_NODE(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn

      do ii=1,lb_dom%nper_wall
!         lb_dom%per_val(ii,:) = &
         fIn(LB_NODE(:,lb_dom%obs_per(ii,1),lb_dom%obs_per(ii,2),lb_dom%obs_per(ii,3))) = &
         fIn(LB_NODE(:,lb_dom%obs_per(ii,4),lb_dom%obs_per(ii,5),lb_dom%obs_per(ii,6)))
      enddo
      
   end subroutine set_periodic
   !------------------------------------------------------------------------




   subroutine set_nrbc(lb_dom,fIn,s_par,prc,meas)
   !------------------------------------------------------------------------
   !
   ! Non-reflecting boundary implementation
   ! Based on work by Giles 
   ! set boundaries as in Andreas Babucke's diss
   ! inlet: set entropic, vortictiy and incoming wave characteristics to 0
   ! outlet: set incoming wave characteristic to 0
   !
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
      integer,dimension(0:3) :: xc,yc,zc
      integer              :: border 
      real(R8B),dimension(LB_NODE(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
      real(R8B)  :: dist(2),m,rho0,u0(3),rhod,ud(3),rhod1,ud1(3),rhod2,ud2(3),rhod3,ud3(3)
      real(R8B),dimension(0:3)  :: c1,c2,c3,c4,c5,rhoc
      real(R8B),dimension(3,0:3)  :: uc
      real(R8B)  :: factoru(1:3),damping 
      real(R8B)  :: usq,cs2inv,t2cs4inv,t2cs2inv

!------------------
! If NRBC_LINEAR_EXTRAPOLATION is defined, use linear instead of parabolic extrapolation
#define NRBC_LINEAR_EXTRAPOLATION
! Note: 
! parabolic interpolation is supposed to work better (higher accuracy) 
! but seems to be instable with the corotating vortex pair problem. BUG???!
! Use linear interpolation until bug is found

! if(s_par%initial ) then ! Not doing anything, because during initialization, no nrbc treatment is required ! endif

      do  i=1,lb_dom%nnr_wall  

         xc(0) = lb_dom%obs_nrwall(i,1)
         yc(0) = lb_dom%obs_nrwall(i,2)
         zc(0) = lb_dom%obs_nrwall(i,3)
         xc(1) = lb_dom%obs_nrwall(i,4)
         yc(1) = lb_dom%obs_nrwall(i,5)
         zc(1) = lb_dom%obs_nrwall(i,6)
         xc(2) = lb_dom%obs_nrwall(i,7)
         yc(2) = lb_dom%obs_nrwall(i,8)
         zc(2) = lb_dom%obs_nrwall(i,9)
         xc(3) = lb_dom%obs_nrwall(i,10)
         yc(3) = lb_dom%obs_nrwall(i,11)
         zc(3) = lb_dom%obs_nrwall(i,12)

         ! calculate current macroscopic variables
         do jj = 0,3
#if defined D2Q9
        rhoc(jj) = fIn(LB_NODE(1,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(2,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(3,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(4,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(5,xc(jj),yc(jj),zc(jj)))  +                                         &
                   fIn(LB_NODE(6,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj))) 

       uc(1,jj) = (fIn(LB_NODE(2,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(4,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(6,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))    & 
                       )/rhoc(jj)    

       uc(2,jj) = (fIn(LB_NODE(3,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(5,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(6,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj)))  - &
                   fIn(LB_NODE(7,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))    & 
                      ) /rhoc(jj)  
       uc(3,jj) = 0.0d0

#elif defined D3Q19 /**/
        rhoc(jj) = fIn(LB_NODE(1,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(2,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(3,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(4,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(5,xc(jj),yc(jj),zc(jj)))  +                                         &
                   fIn(LB_NODE(6,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(10,xc(jj),yc(jj),zc(jj))) +                                         & 
                   fIn(LB_NODE(11,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(12,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(13,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(14,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(15,xc(jj),yc(jj),zc(jj))) +                                         & 
                   fIn(LB_NODE(16,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(17,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(18,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(19,xc(jj),yc(jj),zc(jj)))

       uc(1,jj) = (fIn(LB_NODE(2,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(3,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(10,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(11,xc(jj),yc(jj),zc(jj))) + & 
                   fIn(LB_NODE(12,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(13,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(14,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(15,xc(jj),yc(jj),zc(jj)))   & 
                       )/rhoc(jj)    

       uc(2,jj) = (fIn(LB_NODE(4,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(5,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj)))  - &
                   fIn(LB_NODE(10,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(11,xc(jj),yc(jj),zc(jj))) + & 
                   fIn(LB_NODE(16,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(17,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(18,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(19,xc(jj),yc(jj),zc(jj)))   & 
                      ) /rhoc(jj)     

       uc(3,jj) = (fIn(LB_NODE(6,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(12,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(13,xc(jj),yc(jj),zc(jj))) - &
                   fIn(LB_NODE(14,xc(jj),yc(jj),zc(jj)))  + & 
                   fIn(LB_NODE(15,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(16,xc(jj),yc(jj),zc(jj))) - &
                   fIn(LB_NODE(17,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(18,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(19,xc(jj),yc(jj),zc(jj))))/rhoc(jj)
#endif
         enddo
         ! Calculate the distance of the nodes from the neighbor nodes
         dist(1) = sqrt(real(yc(1) - yc(0))**2 + real(xc(1)-xc(0))**2 + real(zc(1)-zc(0))**2)
         dist(2) = sqrt(real(yc(2) - yc(0))**2 + real(xc(2)-xc(0))**2 + real(zc(2)-zc(0))**2)

         ! set zero values
         rho0       = lb_dom%nrwall_0val(i,0) 
         u0(1:NDIM) = 0.0 ! FIXME lb_dom%nrwall_0val(i,1:NDIM)

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
         !   5: back border   z+
         !   6: front border  z-

         if(btest(lb_dom%obs_nrwall(i,0),0)) then
            border = inlet
         else          
            border = outlet
         endif

         ! --- left  border
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
          else
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
          endif

         ! --- right border
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
          else

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
            
          endif
         else  ! check for top / down up/downstream

         ! --- bottom border
         !     inlet (or outlet), downstream extrapolation

            if(btest(lb_dom%obs_nrwall(i,0),4)) then 
            if(border == inlet) then
               c2=0.d0       ! vorticity fluctuation y is zero        
               c3=0.d0       !                       z 
               c4(0) = 0.d0  ! upstream perturbation is zero
               c5(1) = -ud1(2)*rho0*cs + 1./3.*rhod1
               c5(2) = -ud2(2)*rho0*cs + 1./3.*rhod2
               c5(3) = -ud2(2)*rho0*cs + 1./3.*rhod3
               m = (c4(2) - c4(1))/(dist(2)-dist(1))
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
            else ! border is outlet!
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
            endif !inlet/outlet

         ! --- top border
         !     inlet (or outlet), upstream extrapolation

            elseif(btest(lb_dom%obs_nrwall(i,0),2)) then 
            if(border == inlet) then
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
            else ! border is outlet!
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
            endif !inlet/outlet
            else ! if third dimension has to be extrapolated
#ifdef D3Q19

#endif
            end if


         end if

         ! if the nrbc macr vals should be damped, define
#ifdef NRBC_DAMPING
         damping = NRBC_DAMPING 
#else
         damping = 1.00
#endif

         ! calculate new macroscopic values

         cs2inv  = 1._R8B/cs**2
         t2cs4inv = 1._R8B/(2*cs**4)
         t2cs2inv = 1._R8B/(2*cs**2)

if(s_par%initial) then
   damping = 0.0d0
endif

#define DEBUG_NRBC
#ifdef DEBUG_NRBC
if(xc(0)==10.and. yc(0) ==1) then
if(gtstep_cur <=1) write(76,*)  "#  tstep   rho   u1   u2  u3  " 
write(76,*) gtstep_cur,rhoc(0),uc(1:2,0)
if(gtstep_cur <=1) write(75,*)  "#  tstep   rho   rho1 rho2 rho3"
write(75,*) "# ",gtstep_cur 
write(75,*) rhoc(0)-1.0, uc(1:2,0)
write(75,*) rhoc(1)-1.0, uc(1:2,1)
write(75,*) rhoc(2)-1.0, uc(1:2,2)
write(75,*) rhoc(3)-1.0, uc(1:2,3)
write(75,*) 
write(75,*) 
endif
#endif /* DEBUG_NRBC */
         uc(1:NDIM,0) = damping* (factoru(1:NDIM))/(2*rho0*cs) + u0(1:NDIM) 
         rhoc(0)      = rho0 + damping*(cs2inv*0.5*(c4(0) + c5(0)))
#ifdef DEBUG_NRBC
if(xc(0)==10.and. yc(0) ==1) then
write(77,*) gtstep_cur,rhoc(0),uc(1:2,0)
write(74,*) gtstep_cur,rhoc(0:3)
endif
#endif /* DEBUG_NRBC */

         ! now calculate the equilibrium distribution of the 
         ! calculated macroscopic values

#ifndef NRBC_SPEEDUP
         lb_dom%u(LB_NODE(1:NDIM,xc(0),yc(0),zc(0))) = damping* (factoru(1:NDIM))/(2*rho0*cs) + u0(1:NDIM) 
         lb_dom%rho(xc(0),yc(0),zc(0))      = rho0 + damping*(cs2inv*0.5*(c4(0) + c5(0)))
         call calc_fEq(lb_dom,lb_dom%rho,lb_dom%u,fIn,xc(0),xc(0),yc(0) ,yc(0) ,zc(0),zc(0)) 


#endif /* NRBC_SPEEDUP*/
#ifdef NRBC_SPEEDUP
         uc(1:NDIM,0) = damping* (factoru(1:NDIM))/(2*rho0*cs) + u0(1:NDIM) 
         rhoc(0)      = rho0 + damping*(cs2inv*0.5*(c4(0) + c5(0)))
         ! calculate  equilibrium distribution
         usq =  (uc(1,0)*uc(1,0) + uc(2,0)*uc(2,0) + uc(3,0)*uc(3,0))*t2cs2inv
#if defined D2Q9
             fIn(LB_NODE(1,xc(0),yc(0),zc(0))) = t(1)*rhoc(0) *(1._R8B - usq)
             fIn(LB_NODE(2,xc(0),yc(0),zc(0))) = t(2)*rhoc(0) *(1._R8B + (uc(1,0))*cs2inv &
            + (uc(1,0))**2*t2cs4inv   - usq)
             fIn(LB_NODE(3,xc(0),yc(0),zc(0))) = t(3)*rhoc(0)*(1._R8B - (-uc(2,0))*cs2inv &
            + (uc(2,0))**2*t2cs4inv   - usq)
             fIn(LB_NODE(4,xc(0),yc(0),zc(0))) = t(4)*rhoc(0)*(1._R8B + (-uc(1,0))*cs2inv &
           + (-uc(1,0))**2*t2cs4inv  - usq)
             fIn(LB_NODE(5,xc(0),yc(0),zc(0))) = t(5)*rhoc(0)*(1._R8B + (-uc(2,0))*cs2inv &
           + (-uc(2,0))**2*t2cs4inv   - usq)
             fIn(LB_NODE(6,xc(0),yc(0),zc(0))) = t(6)*rhoc(0)*(1._R8B  &
            + (uc(1,0) + uc(2,0))*cs2inv &
            + (uc(1,0)+uc(2,0))**2*t2cs4inv  - usq)
             fIn(LB_NODE(7,xc(0),yc(0),zc(0))) = t(7)*rhoc(0)*(1._R8B  &
           + (-uc(1,0) + uc(2,0))*cs2inv &
           + (-uc(1,0)+uc(2,0))**2*t2cs4inv  - usq)
             fIn(LB_NODE(8,xc(0),yc(0),zc(0))) = t(8)*rhoc(0)*(1._R8B  &
           + (-uc(1,0) -uc(2,0))*cs2inv &
           + (-uc(1,0)-uc(2,0))**2*t2cs4inv    - usq)
             fIn(LB_NODE(9,xc(0),yc(0),zc(0))) = t(9)*rhoc(0)*(1._R8B  &
            + (uc(1,0) -uc(2,0))*cs2inv &
            + (uc(1,0)-uc(2,0))**2*t2cs4inv  - usq)
#elif defined D3Q19
           fIn(LB_NODE(1,xc(0),yc(0),zc(0))) = t(1)*rhoc(0)*(1  &
                - usq) 
           fIn(LB_NODE(2,xc(0),yc(0),zc(0))) = t(2)*rhoc(0)*(1._R8B + (uc(1,0))*cs2inv &
                + (uc(1,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(3,xc(0),yc(0),zc(0))) = t(3)*rhoc(0)*(1._R8B + (-uc(1,0))*cs2inv &
                -uc(1,0)**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(4,xc(0),yc(0),zc(0))) = t(4)*rhoc(0)*(1._R8B + uc(2,0)*cs2inv &
                + (uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(5,xc(0),yc(0),zc(0))) = t(5)*rhoc(0)*(1._R8B - uc(2,0)*cs2inv &
                + (-uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(6,xc(0),yc(0),zc(0))) = t(6)*rhoc(0)*(1._R8B + (uc(3,0))*cs2inv &
                + (uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(7,xc(0),yc(0),zc(0))) = t(7)*rhoc(0)*(1._R8B + (-uc(3,0))*cs2inv &
                + (-uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(8,xc(0),yc(0),zc(0))) = t(8)*rhoc(0)*(1._R8B  &
                + (uc(1,0) + uc(2,0))*cs2inv &
                + (uc(1,0) + uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(9,xc(0),yc(0),zc(0))) = t(9)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) -uc(2,0))*cs2inv &
                + (-uc(1,0) -uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(10,xc(0),yc(0),zc(0))) = t(10)*rhoc(0)*(1._R8B  &
                + (uc(1,0) -uc(2,0))*cs2inv &
                + (uc(1,0) -uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(11,xc(0),yc(0),zc(0))) = t(11)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) + uc(2,0))*cs2inv &
                + (-uc(1,0) + uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(12,xc(0),yc(0),zc(0))) = t(12)*rhoc(0)*(1._R8B  &
                + (uc(1,0) + uc(3,0))*cs2inv &
                + (uc(1,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(13,xc(0),yc(0),zc(0))) = t(13)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) -uc(3,0))*cs2inv &
                + (-uc(1,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(14,xc(0),yc(0),zc(0))) = t(14)*rhoc(0)*(1._R8B  &
                + (uc(1,0) -uc(3,0))*cs2inv &
                + (uc(1,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(15,xc(0),yc(0),zc(0))) = t(15)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) + uc(3,0))*cs2inv &
                + (-uc(1,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(16,xc(0),yc(0),zc(0))) = t(16)*rhoc(0)*(1._R8B  &
                + (uc(2,0) + uc(3,0))*cs2inv &
                + (uc(2,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(17,xc(0),yc(0),zc(0))) = t(17)*rhoc(0)*(1._R8B  &
                + (-uc(2,0) -uc(3,0))*cs2inv &
                + (-uc(2,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(18,xc(0),yc(0),zc(0))) = t(18)*rhoc(0)*(1._R8B  &
                + (uc(2,0) -uc(3,0))*cs2inv &
                + (uc(2,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(19,xc(0),yc(0),zc(0))) = t(19)*rhoc(0)*(1._R8B  &
                + (-uc(2,0) + uc(3,0))*cs2inv &
                + (-uc(2,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
#endif /* D3Q19 */
#endif /* NRBC_SPEEDUP */     
      end do

   end subroutine set_nrbc
   !------------------------------------------------------------------------






#ifdef ORG
   subroutine set_nrbc(lb_dom,state,fIn,rho,u,s_par,prc,meas)
   !------------------------------------------------------------------------
   !
   ! Non-reflecting boundary implementation
   ! Based on work by Giles 
   ! set boundaries as in Andreas Babucke's diss
   ! inlet: set entropic, vortictiy and incoming wave characteristics to 0
   ! outlet: set incoming wave characteristic to 0
   !

      implicit none
      type(lb_block)       :: lb_dom
      type(sim_parameter)  :: s_par 
      type(mpl_var)        :: prc
      type(measure)        :: meas
      integer              :: state(lb_dom%lx(1),lb_dom%lx(2),lb_dom%lx(3))
      integer              :: i,j,k,l,jj
!      integer              :: x_start,x_end,y_start,y_end,z_start,z_end,incr_y,incr_z,lx,ly,lz
      integer,dimension(-1:3) :: xc,yc,zc
      integer              :: border 
      real(R8B),dimension(LB_NODE(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
      real(R8B),dimension(LB_NODE(1:NDIM,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: u
      real(R8B),dimension(0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1)    :: rho
      real(R8B)  :: dist(2),m,rho0,u0(3),rhod,ud(3),rhod1,ud1(3),rhod2,ud2(3),rhod3,ud3(3)
      real(R8B),dimension(0:3)  :: c1,c2,c3,c4,c5,rhoc
      real(R8B),dimension(3,0:3)  :: uc
      real(R8B)  :: factoru(1:3),damping 
      real(R8B)  :: usq,cs2inv,t2cs4inv,t2cs2inv
      real(R8B)  :: rholoc 

!------------------
! If NRBC_LINEAR_EXTRAPOLATION is defined, use linear instead of parabolic extrapolation
#define NRBC_LINEAR_EXTRAPOLATION
! Note: 
! parabolic interpolation is supposed to work better (higher accuracy) 
! but seems to be instable with the corotating vortex pair problem. BUG???!
! Use linear interpolation until bug is found
  
      do  i=1,lb_dom%nnr_wall  
         ! Values should be all copied to the ghost nodes because propagation
         ! pulls wrong values without 
         ! FIXME: The corners are not copied yet!
         xc(-1) = lb_dom%obs_nrwall(i,1)-(lb_dom%obs_nrwall(i,4)-lb_dom%obs_nrwall(i,1))
         yc(-1) = lb_dom%obs_nrwall(i,2)-(lb_dom%obs_nrwall(i,5)-lb_dom%obs_nrwall(i,2))
         zc(-1) = lb_dom%obs_nrwall(i,3)-(lb_dom%obs_nrwall(i,6)-lb_dom%obs_nrwall(i,3))

         xc(0) = lb_dom%obs_nrwall(i,1)
         yc(0) = lb_dom%obs_nrwall(i,2)
         zc(0) = lb_dom%obs_nrwall(i,3)
         xc(1) = lb_dom%obs_nrwall(i,4)
         yc(1) = lb_dom%obs_nrwall(i,5)
         zc(1) = lb_dom%obs_nrwall(i,6)
         xc(2) = lb_dom%obs_nrwall(i,7)
         yc(2) = lb_dom%obs_nrwall(i,8)
         zc(2) = lb_dom%obs_nrwall(i,9)
         xc(3) = lb_dom%obs_nrwall(i,10)
         yc(3) = lb_dom%obs_nrwall(i,11)
         zc(3) = lb_dom%obs_nrwall(i,12)

         ! First, copy all the values of the nrbc node to the ghost node for
         ! propagation
!        fIn(LB_NODE(:,xc(-1),yc(-1),zc(-1))) =  -0.001 !fIn(LB_NODE(:,xc(0),yc(0),zc(0))) 
        fIn(LB_NODE(:,xc(-1),yc(-1),zc(-1))) =  fIn(LB_NODE(:,xc(0),yc(0),zc(0))) 
         ! calculate current macroscopic variables
         do jj = 0,3
#if defined D2Q9
        rhoc(jj) = fIn(LB_NODE(1,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(2,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(3,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(4,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(5,xc(jj),yc(jj),zc(jj)))  +                                         &
                   fIn(LB_NODE(6,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj))) 

       uc(1,jj) = (fIn(LB_NODE(2,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(4,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(6,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))    & 
                       )/rhoc(jj)    

       uc(2,jj) = (fIn(LB_NODE(3,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(5,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(6,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj)))  - &
                   fIn(LB_NODE(7,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))    & 
                      ) /rhoc(jj)  
       uc(3,jj) = 0.0d0

#elif defined D3Q19 /**/
        rhoc(jj) = fIn(LB_NODE(1,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(2,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(3,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(4,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(5,xc(jj),yc(jj),zc(jj)))  +                                         &
                   fIn(LB_NODE(6,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))  + fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(10,xc(jj),yc(jj),zc(jj))) +                                         & 
                   fIn(LB_NODE(11,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(12,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(13,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(14,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(15,xc(jj),yc(jj),zc(jj))) +                                         & 
                   fIn(LB_NODE(16,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(17,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(18,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(19,xc(jj),yc(jj),zc(jj)))

       uc(1,jj) = (fIn(LB_NODE(2,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(3,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(10,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(11,xc(jj),yc(jj),zc(jj))) + & 
                   fIn(LB_NODE(12,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(13,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(14,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(15,xc(jj),yc(jj),zc(jj)))   & 
                       )/rhoc(jj)    

       uc(2,jj) = (fIn(LB_NODE(4,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(5,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(8,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(9,xc(jj),yc(jj),zc(jj)))  - &
                   fIn(LB_NODE(10,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(11,xc(jj),yc(jj),zc(jj))) + & 
                   fIn(LB_NODE(16,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(17,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(18,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(19,xc(jj),yc(jj),zc(jj)))   & 
                      ) /rhoc(jj)     

       uc(3,jj) = (fIn(LB_NODE(6,xc(jj),yc(jj),zc(jj)))  - fIn(LB_NODE(7,xc(jj),yc(jj),zc(jj)))  + &
                   fIn(LB_NODE(12,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(13,xc(jj),yc(jj),zc(jj))) - &
                   fIn(LB_NODE(14,xc(jj),yc(jj),zc(jj)))  + & 
                   fIn(LB_NODE(15,xc(jj),yc(jj),zc(jj))) + fIn(LB_NODE(16,xc(jj),yc(jj),zc(jj))) - &
                   fIn(LB_NODE(17,xc(jj),yc(jj),zc(jj))) - fIn(LB_NODE(18,xc(jj),yc(jj),zc(jj))) + &
                   fIn(LB_NODE(19,xc(jj),yc(jj),zc(jj))))/rhoc(jj)
#endif
         enddo
         ! Calculate the distance of the nodes from the neighbor nodes
         dist(1) = sqrt(real(yc(1) - yc(0))**2 + real(xc(1)-xc(0))**2 + real(zc(1)-zc(0))**2)
         dist(2) = sqrt(real(yc(2) - yc(0))**2 + real(xc(2)-xc(0))**2 + real(zc(2)-zc(0))**2)

         ! set zero values
         rho0       = lb_dom%nrwall_0val(i,0) 
!         rho0       = rhoc(0) 
         u0(1:NDIM) = lb_dom%nrwall_0val(i,1:NDIM)
         

         !FIXME REMOVE!!!
!         rhoc(0) = rhoc(1)
!         uc(:,0) = uc(:,1)
         
!         if(btest(lb_dom%state(xc(0),yc(0),zc(0)),nr_wall) .and. yc(0)==1 &
!         .and.  xc(0)==3) then 
!           write(*,*) gtstep_cur,"nrbc detected",xc(0),yc(0)
!           write(*,*) "fIn ",fIn(LB_NODE(:,xc(0),yc(0),zc(0)))
!           write(*,*) "rho0",rho0
!           write(*,*) "rhoc",rhoc(0)
!           write(*,*) 
!           write(*,*) "u0  ",u0(1:2)
!           write(*,*) "uc  ",uc(1:2,0)
!stop
!               endif


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
         !   5: back border   z+
         !   6: front border  z-

         border = outlet
         border = inlet

         ! --- left  border
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
          else
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
          endif

         ! --- right border
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
          else

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
            
          endif
         else  ! check for top / down up/downstream

         ! --- bottom border
         !     inlet (or outlet), downstream extrapolation

            if(btest(lb_dom%obs_nrwall(i,0),4)) then 
            if(border == inlet) then
               c2=0.d0       ! vorticity fluctuation y is zero        
               c3=0.d0       !                       z 
               c4(0) = 0.d0  ! upstream perturbation is zero
               c5(1) = -ud1(2)*rho0*cs + 1./3.*rhod1
               c5(2) = -ud2(2)*rho0*cs + 1./3.*rhod2
               c5(3) = -ud2(2)*rho0*cs + 1./3.*rhod3
               m = (c4(2) - c4(1))/(dist(2)-dist(1))
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
            else ! border is outlet!
               c4(0) = 0.d0     ! upstream perturbation is zero

               c2(1) = ud1(1)*rho0*cs
               c2(2) = ud2(1)*rho0*cs
               c2(3) = ud3(1)*rho0*cs
               c3(1) = ud1(3)*rho0*cs
               c3(2) = ud2(3)*rho0*cs
               c3(3) = ud3(3)*rho0*cs
               c5(1) = -ud1(2)*rho0*cs + 1./3.*rhod1
               c5(2) = -ud2(2)*rho0*cs + 1./3.*rhod2
               c5(3) = -ud3(2)*rho0*cs + 1./3.*rhod3

#ifdef NRBC_LINEAR_EXTRAPOLATION
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
               factoru(2)=c4(0)-c5(0)
               factoru(3)=c3(0)
            endif !inlet/outlet

         ! --- top border
         !     inlet (or outlet), upstream extrapolation

            elseif(btest(lb_dom%obs_nrwall(i,0),2)) then 
            if(border == inlet) then
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
            else ! border is outlet!
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
            endif !inlet/outlet
            else ! if third dimension has to be extrapolated
#ifdef D3Q19

#endif
            end if


         end if

         ! if the nrbc macr vals should be damped, define
#ifdef NRBC_DAMPING
         damping = NRBC_DAMPING 
#else
         damping = 1.00
#endif

         ! calculate new macroscopic values

         cs2inv  = 1._R8B/cs**2
         t2cs4inv = 1._R8B/(2*cs**4)
         t2cs2inv = 1._R8B/(2*cs**2)

         u(LB_NODE(1:NDIM,xc(0),yc(0),zc(0))) =damping* (factoru(1:NDIM))/(2*rho0*cs) + u0(1:NDIM) 
         uc(1:NDIM,0) = u(LB_NODE(1:NDIM,xc(0),yc(0),zc(0)))
         rho(xc(0),yc(0),zc(0)) = rho0 + damping*(cs2inv*0.5*(c4(0) + c5(0)))
         rhoc(0) = rho(xc(0),yc(0),zc(0))
!         if(btest(lb_dom%state(xc(0),yc(0),zc(0)),nr_wall) .and. yc(0)==1 &
!         .and.  xc(0)==3) then 
!            write(*,*) "new value rhoc",rhoc(0)
!            write(*,*) "            uc",uc(1:2,0)
!            write(*,*) 
!            write(*,*) 
!         endif

         ! now calculate the equilibrium distribution of the 
         ! calculated macroscopic values
!         call calc_fEq(lb_dom,rho,u,fIn,xc(0),xc(0),yc(0) ,yc(0) ,zc(0),zc(0)) 


         ! calculate  equilibrium distribution

         usq =  (uc(1,0)*uc(1,0) + uc(2,0)*uc(2,0) + uc(3,0)*uc(3,0))*t2cs2inv
#if defined D2Q9
             fIn(LB_NODE(1,xc(0),yc(0),zc(0))) = t(1)*rhoc(0) *(1._R8B - usq)
             fIn(LB_NODE(2,xc(0),yc(0),zc(0))) = t(2)*rhoc(0) *(1._R8B + (uc(1,0))*cs2inv &
            + (uc(1,0))**2*t2cs4inv   - usq)
             fIn(LB_NODE(3,xc(0),yc(0),zc(0))) = t(3)*rhoc(0)*(1._R8B - (-uc(2,0))*cs2inv &
            + (uc(2,0))**2*t2cs4inv   - usq)
             fIn(LB_NODE(4,xc(0),yc(0),zc(0))) = t(4)*rhoc(0)*(1._R8B + (-uc(1,0))*cs2inv &
           + (-uc(1,0))**2*t2cs4inv  - usq)
             fIn(LB_NODE(5,xc(0),yc(0),zc(0))) = t(5)*rhoc(0)*(1._R8B + (-uc(2,0))*cs2inv &
           + (-uc(2,0))**2*t2cs4inv   - usq)
             fIn(LB_NODE(6,xc(0),yc(0),zc(0))) = t(6)*rhoc(0)*(1._R8B  &
            + (uc(1,0) + uc(2,0))*cs2inv &
            + (uc(1,0)+uc(2,0))**2*t2cs4inv  - usq)
             fIn(LB_NODE(7,xc(0),yc(0),zc(0))) = t(7)*rhoc(0)*(1._R8B  &
           + (-uc(1,0) + uc(2,0))*cs2inv &
           + (-uc(1,0)+uc(2,0))**2*t2cs4inv  - usq)
             fIn(LB_NODE(8,xc(0),yc(0),zc(0))) = t(8)*rhoc(0)*(1._R8B  &
           + (-uc(1,0) -uc(2,0))*cs2inv &
           + (-uc(1,0)-uc(2,0))**2*t2cs4inv    - usq)
             fIn(LB_NODE(9,xc(0),yc(0),zc(0))) = t(9)*rhoc(0)*(1._R8B  &
            + (uc(1,0) -uc(2,0))*cs2inv &
            + (uc(1,0)-uc(2,0))**2*t2cs4inv  - usq)
#elif defined D3Q19
           fIn(LB_NODE(1,xc(0),yc(0),zc(0))) = t(1)*rhoc(0)*(1  &
                - usq) 
           fIn(LB_NODE(2,xc(0),yc(0),zc(0))) = t(2)*rhoc(0)*(1._R8B + (uc(1,0))*cs2inv &
                + (uc(1,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(3,xc(0),yc(0),zc(0))) = t(3)*rhoc(0)*(1._R8B + (-uc(1,0))*cs2inv &
                -uc(1,0)**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(4,xc(0),yc(0),zc(0))) = t(4)*rhoc(0)*(1._R8B + uc(2,0)*cs2inv &
                + (uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(5,xc(0),yc(0),zc(0))) = t(5)*rhoc(0)*(1._R8B - uc(2,0)*cs2inv &
                + (-uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(6,xc(0),yc(0),zc(0))) = t(6)*rhoc(0)*(1._R8B + (uc(3,0))*cs2inv &
                + (uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(7,xc(0),yc(0),zc(0))) = t(7)*rhoc(0)*(1._R8B + (-uc(3,0))*cs2inv &
                + (-uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(8,xc(0),yc(0),zc(0))) = t(8)*rhoc(0)*(1._R8B  &
                + (uc(1,0) + uc(2,0))*cs2inv &
                + (uc(1,0) + uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(9,xc(0),yc(0),zc(0))) = t(9)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) -uc(2,0))*cs2inv &
                + (-uc(1,0) -uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(10,xc(0),yc(0),zc(0))) = t(10)*rhoc(0)*(1._R8B  &
                + (uc(1,0) -uc(2,0))*cs2inv &
                + (uc(1,0) -uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(11,xc(0),yc(0),zc(0))) = t(11)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) + uc(2,0))*cs2inv &
                + (-uc(1,0) + uc(2,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(12,xc(0),yc(0),zc(0))) = t(12)*rhoc(0)*(1._R8B  &
                + (uc(1,0) + uc(3,0))*cs2inv &
                + (uc(1,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(13,xc(0),yc(0),zc(0))) = t(13)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) -uc(3,0))*cs2inv &
                + (-uc(1,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(14,xc(0),yc(0),zc(0))) = t(14)*rhoc(0)*(1._R8B  &
                + (uc(1,0) -uc(3,0))*cs2inv &
                + (uc(1,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(15,xc(0),yc(0),zc(0))) = t(15)*rhoc(0)*(1._R8B  &
                + (-uc(1,0) + uc(3,0))*cs2inv &
                + (-uc(1,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(16,xc(0),yc(0),zc(0))) = t(16)*rhoc(0)*(1._R8B  &
                + (uc(2,0) + uc(3,0))*cs2inv &
                + (uc(2,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(17,xc(0),yc(0),zc(0))) = t(17)*rhoc(0)*(1._R8B  &
                + (-uc(2,0) -uc(3,0))*cs2inv &
                + (-uc(2,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(18,xc(0),yc(0),zc(0))) = t(18)*rhoc(0)*(1._R8B  &
                + (uc(2,0) -uc(3,0))*cs2inv &
                + (uc(2,0) -uc(3,0))**2*t2cs4inv  &
                - usq) 
           fIn(LB_NODE(19,xc(0),yc(0),zc(0))) = t(19)*rhoc(0)*(1._R8B  &
                + (-uc(2,0) + uc(3,0))*cs2inv &
                + (-uc(2,0) + uc(3,0))**2*t2cs4inv  &
                - usq) 
#endif /* D3Q19 */
     
      end do

   end subroutine set_nrbc
   !------------------------------------------------------------------------
#endif /* ORG */

   subroutine microphone(lb_dom,pos,fIn,s_par,prc,tstep)
   !------------------------------------------------------------------------
   !
   ! Get current density at a certain place 
   ! 

      implicit none
      type(lb_block)       :: lb_dom
      type(sim_parameter)  :: s_par 
      type(mpl_var)        :: prc
      integer              :: tstep,ll
      integer              :: xmin(3),xmax(3) 
      integer,intent(in)   :: pos(3)
      real(R8B)            :: rholoc
      real(R8B),dimension(LB_NODE(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
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
            rholoc = rholoc+fIn(LB_NODE(ll,pos(1)-xmin(1)+1,pos(2)-xmin(2)+1,pos(3)-xmin(3)+1))
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

      open(44,file='mic_'//trim(ch_posx)//'_'//trim(ch_posy)//'_'//trim(ch_posz)//'_'//trim(ch_rank)//'.res',position='append')
         write(44,*) 
         write(44,*)
         write(44,*) "# Density Fluctuation for ",trim(s_par%problem_name)
         if(record .eqv. .true. ) then
         write(44,*) "# Sensing pos 1    was x: ",pos(1),pos(1)-xmin(1)+1 
         write(44,*) "#                      y: ",pos(2),pos(2)-xmin(2)+1 
         write(44,*) "#                      z: ",pos(3),pos(3)-xmin(3)+1 
         endif
      write(44,*) "#                  Omega: ",s_par%omega
      write(44,*) "#                   umax: ",s_par%umax 
      write(44,*) "#                     dx: ",s_par%gx(1)
      write(44,*) "#                   dist: ",s_par%obst_l
      write(44,*) "#       Relaxation model: ",s_par%modelname
      write(44,*) "# tstep     dens   " 
   endif
   if(tstep>0) then
      write(44,'(i10,3e15.7)') tstep,rholoc
   endif
      if(s_par%goend .eqv. .true. .or. gtstep_cur == s_par%tMax) then
         close(44)
      endif
   endif

end subroutine microphone





subroutine set_bnd(lb_dom,state,fIn,rho,u,s_par,prc,meas,tStep)


!------------------------------------------------------------------------
!
! set the boundary conditions in each timestep.
! equilibrium distribution is assumed at the boundaries
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
   real(R8B),dimension(LB_NODE(1:nnod,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: fIn
   real(R8B),dimension(LB_NODE(1:NDIM,0:lb_dom%lx(1)+1,0:lb_dom%lx(2)+1,0:lb_dom%lx(3)+1))    :: u
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

   call set_nrbc(lb_dom,fIn,s_par,prc,meas)
   call set_periodic(lb_dom,fIn,s_par,prc,meas)
   

!--------------------------------------------------------------------
! Select defined problem  

      select case(s_par%problem)



!----------------------------------------
! Flat plate with trailing edge

      case(plate_trail)
!         call set_nrbc(lb_dom,state,fIn,rho,u,s_par,prc,meas)
 

!----------------------------------------
! Gaussian Pulse (one source) 

      case(gaussian)
         ! set period of sine for gaussian pulse
!         t_per=20

         ! set sinusodial distribution in middle.
!         x_b = (/1*lx/2-1,1*ly/2-1,lz/2-1/)
!         x_e = (/1*lx/2+1,1*ly/2+1,lz/2+1/)
!         x_b = (/1*lx/3-1,1*ly/3-1,lz/2+1/)
!         x_e = (/1*lx/3+1,1*ly/3+1,lz/2+1/)
#ifdef D2Q9
!         x_b(3) = 1
!         x_e(3) = 1
#endif
           
!         if(tStep <= t_per) then
!            if(prc%crd(1) == 0 .and. prc%crd(2) ==prc%np(2)-1 .and. prc%crd(3) == 0) then
!            rho(x_b(1):x_e(1),x_b(2):x_e(2),x_b(3):x_e(3)) = 1.+(s_par%umax/100.d0)*sin(pi/dble(t_per)*dble(tStep))
!            write(49,*) tstep, rho(x_b(1),x_b(2),x_b(3))  
!            call calc_fEq(lb_dom,rho,u,fIn,x_b(1),x_e(1),x_b(2),x_e(2),x_b(3),x_e(3)) 
!             end if
!          end if

      case(gauss_convect)
         ! set period of sine for gaussian pulse
!         t_per=20

         ! set sinusodial distribution in middle.
!         x_b = (/1*lx/2-1,1*ly/2-1,lz/2-1/)
!         x_e = (/1*lx/2+1,1*ly/2+1,lz/2+1/)
#ifdef D2Q9
!         x_b(3) = 1
!         x_e(3) = 1
#endif
           
!         if(tStep <= t_per) then
!            rho(x_b(1):x_e(1),x_b(2):x_e(2),x_b(3):x_e(3)) = 1.+0.1*sin(pi/dble(t_per)*dble(tStep))
!            call calc_fEq(lb_dom,rho,u,fIn,x_b(1),x_e(1),x_b(2),x_e(2),x_b(3),x_e(3)) 
!          end if

!----------------------------------------
! Corotating vortex 

      case(corotating_vortex)


        pos = (/ int(5.d0/8.d0*real(s_par%gx(1))), s_par%gx(2)/2, s_par%gx(3)/2+1    /)
        call microphone(lb_dom,pos,fIn,s_par,prc,tstep)
        pos = (/ int(6.d0/8.d0*real(s_par%gx(1))), s_par%gx(2)/2, s_par%gx(3)/2+1    /)
        call microphone(lb_dom,pos,fIn,s_par,prc,tstep)
        pos = (/ int(7.d0/8.d0*real(s_par%gx(1))), s_par%gx(2)/2, s_par%gx(3)/2+1    /)
        call microphone(lb_dom,pos,fIn,s_par,prc,tstep)



   

!----------------------------------------
! Cylinder in channel test with eq boundaries 

      case(100)
! unstable boundary condition!! better use zou
         do k=z_start,z_end
         if(prc%crd(1)== 0) then
         do j=y_start,y_end
            i=1
            u(LB_NODE(2,i,j,k)) = 0._R8B 
#if defined D2Q9
            u(LB_NODE(1,i,j,k)) = s_par%umax*get_parabolic2d(real(j+incr_y,R8B),real(s_par%gx(2)+1,R8B))
#elif defined D3Q19
!FIXME errorous u calculation
            u(LB_NODE(1,i,j,k)) = s_par%umax*get_parabolic3d(real(j+incr_y,R8B) &
            ,real(s_par%gx(2)+1,R8B),real(k+incr_z,R8B),real(s_par%gx(3)+1,R8B))
            u(LB_NODE(3,i,j,k)) = 0._R8B 
#endif
            rho(i,j,k) = 1._R8B 
            call calc_fEq(lb_dom,rho,u,fIn,i,i,j,j,k,k)
         end do
         end if
            if(prc%crd(1)== prc%np(1) -1) then
         do j=y_start,y_end
            i=lx
            u(LB_NODE(2,i,j,k)) = 0._R8B 
#if defined D2Q9
            u(LB_NODE(1,i,j,k)) = s_par%umax*get_parabolic2d(real(j+incr_y,R8B),real(s_par%gx(2)+1,R8B))
#elif defined D3Q19
            u(LB_NODE(1,i,j,k)) = s_par%umax*get_parabolic3d(real(j+incr_y,R8B)        &
         ,real(s_par%gx(2)+1,R8B),real(k+incr_z,R8B),real(s_par%gx(3)+1,R8B))
            u(LB_NODE(3,i,j,k)) = 0._R8B 
#endif
            rho(i,j,k) = 1._R8B 
            call calc_fEq(lb_dom,rho,u,fIn,i,i,j,j,k,k)
         end do
         end if
         enddo          



!----------------------------------------
! Cylinder in channel

      case(cylinder)
#define CYLINDER_NRBC
#ifdef  CYLINDER_NRBC

   !---------------------------
   ! Non-reflective boundary condition treatment at open endings

!          call set_nrbc(lb_dom,state,fIn,rho,u,s_par,prc,meas)

#else

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
               u(LB_NODE(1,i,j,k))   =  s_par%umax*get_parabolic2d(real(j+incr_y,R8B),real(s_par%gx(2)+1,R8B))!4._R8B*s_par%umax/real(s_par%gx(2)+1)**2*&
               !      & ((real(j+incr_y))*real(s_par%gx(2)+1)-real(j+incr_y)*real(j+incr_y))
               u(LB_NODE(2,i,j,k))   = 0._R8B
               rho(i,j,k)   = 1._R8B / (1-u(LB_NODE(1,i,j,k)))*(fIn(LB_NODE(1,i,j,k)) &
            + fIn(LB_NODE(3,i,j,k)) + fIn(LB_NODE(5,i,j,k)) &
            + 2*(fIn(LB_NODE(4,i,j,k)) + fIn(LB_NODE(7,i,j,k)) + fIn(LB_NODE(8,i,j,k))))
 !              u(LB_NODE(:,i-1,j,k)) = u(LB_NODE(:,i,j,k))
 !              rho(i-1,j,k) = rho(i,j,k)
!               write(*,*) i,j,k,"dichte",rho(i,j,k),"geschwindigkeit",u(LB_NODE(:,i,j,k)) 
               enddo
            end if
            ! OUTLET: reset to constant pressure
            if(prc%crd(1)== prc%np(1) -1) then
         do i=lx,lx
               rho(i,j,k) = 1._R8B
               u(LB_NODE(1,i,j,k)) = -1 + 1._R8B/rho(i,j,k)*(fIn(LB_NODE(1,i,j,k))&
                  +fIn(LB_NODE(3,i,j,k))+fIn(LB_NODE(5,i,j,k)) &
                  + 2*(fIn(LB_NODE(2,i,j,k))+fIn(LB_NODE(6,i,j,k))+fIn(LB_NODE(9,i,j,k))))
               u(LB_NODE(2,i,j,k))   = 0._R8B
!            u(LB_NODE(:,i+1,j,k)) = u(LB_NODE(:,i,j,k))
!            rho(i+1,j,k) = rho(i,j,k)
!            write(*,*) i,j,k,"dichte",rho(i,j,k),"geschwindigkeit",u(LB_NODE(:,i,j,k)) 
         end do
            end if
         end do
         ! ---- microscopic boundary condition
         do j=y_start,y_end
            ! INLET: Zou/He BC
            if(prc%crd(1) == 0) then 
         do i=1,1
            fIn(LB_NODE(2,i,j,k)) = fIn(LB_NODE(4,i,j,k)) + 2._R8B/3._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k))
            fIn(LB_NODE(6,i,j,k)) = fIn(LB_NODE(8,i,j,k)) + 1._R8B/2._R8B*(fIn(LB_NODE(5,i,j,k))-fIn(LB_NODE(3,i,j,k))) &
                 &                  + 1._R8B/2._R8B*rho(i,j,k)*u(LB_NODE(2,i,j,k) )&
                 &                  + 1._R8B/6._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k)) 
            fIn(LB_NODE(9,i,j,k)) = fIn(LB_NODE(7,i,j,k)) + 1._R8B/2._R8B*(fIn(LB_NODE(3,i,j,k))-fIn(LB_NODE(5,i,j,k))) &
                &                  - 1._R8B/2._R8B*rho(i,j,k)*u(LB_NODE(2,i,j,k)) &
                &                  + 1._R8B/6._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k)) 
!            fIn(LB_NODE(:,i-1,:,:)) = fIn(LB_NODE(:,i,:,:))
         end do
            end if
           ! OUTLET:Zou/He BC
            if(prc%crd(1)== prc%np(1) -1) then
         do i=lx,lx
           fIn(LB_NODE(4,i,j,k)) = fIn(LB_NODE(2,i,j,k) )- 2._R8B/3._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k))
           fIn(LB_NODE(8,i,j,k)) = fIn(LB_NODE(6,i,j,k)) + 1._R8B/2._R8B*(fIn(LB_NODE(3,i,j,k))-fIn(LB_NODE(5,i,j,k))) &
                &                  - 1._R8B/2._R8B*rho(i,j,k)*u(LB_NODE(2,i,j,k)) &
                &                  - 1._R8B/6._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k) )
           fIn(LB_NODE(7,i,j,k)) = fIn(LB_NODE(9,i,j,k)) + 1._R8B/2._R8B*(fIn(LB_NODE(5,i,j,k))-fIn(LB_NODE(3,i,j,k))) &
                &                  + 1._R8B/2._R8B*rho(i,j,k)*u(LB_NODE(2,i,j,k)) &
                &                  - 1._R8B/6._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k)) 
        end do
            end if
        end do
#endif
#endif /* NRBC */

!----------------------------------------
! Flute 

      case(flute)

        pos = (/ lb_dom%lx(1)/2+1, lb_dom%lx(2)-10, lb_dom%lx(3)/2+1    /)
        call microphone(lb_dom,pos,fIn,s_par,prc,tstep)

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
!               u(LB_NODE(1,i,j,k))   = 4._R8B*s_par%umax/real(lb_dom%num_in(2)+1)**2*&
!                     & ((real(j))*real(lb_dom%num_in(2)+1)-real(j)*real(j))
!FIXME this is my own poisseuille profile
               u(LB_NODE(1,i,j,k)) = s_par%umax*get_parabolic2d(real(j+1,R8B),real(2*s_par%obst_r+1,R8B))  
               u(LB_NODE(2,i,j,k))   = 0._R8B
               rho(i,j,k)   = 1._R8B / (1-u(LB_NODE(1,i,j,k)))*(fIn(LB_NODE(1,i,j,k)) &
               + fIn(LB_NODE(3,i,j,k)) + fIn(LB_NODE(5,i,j,k)) &
               + 2*(fIn(LB_NODE(4,i,j,k)) + fIn(LB_NODE(7,i,j,k)) + fIn(LB_NODE(8,i,j,k))))
            end if
            end do
            ! OUTLET: reset to constant pressure
            i=lx
            if(lb_dom%state(i,j,1) == outlet) then
               rho(i,j,k) = 1._R8B
               u(LB_NODE(1,i,j,k)) = 0.d0 !-1 + 1._R8B/rho(i,j,k)*(fIn(LB_NODE(1,i,j,k))+fIn(LB_NODE(3,i,j,k))+fIn(LB_NODE(5,i,j,k)) &
                                           !  & + 2*(fIn(LB_NODE(2,i,j,k))+fIn(LB_NODE(6,i,j,k))+fIn(LB_NODE(9,i,j,k))))
               u(LB_NODE(2,i,j,k))   = 0._R8B
                       
            end if
         end do

         

         ! ---- microscopic boundary condition
         do j=y_start,y_end
            ! INLET: Zou/He BC
            i = 1
            if(lb_dom%state(i,j,k) == inlet) then 
            fIn(LB_NODE(2,i,j,k)) = fIn(LB_NODE(4,i,j,k)) + 2._R8B/3._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k))
            fIn(LB_NODE(6,i,j,k)) = fIn(LB_NODE(8,i,j,k)) + 1._R8B/2._R8B*(fIn(LB_NODE(5,i,j,k))-fIn(LB_NODE(3,i,j,k))) &
                 &                  + 1._R8B/2._R8B*rho(i,j,k)*u(LB_NODE(2,i,j,k) )&
                 &                  + 1._R8B/6._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k)) 
            fIn(LB_NODE(9,i,j,k)) = fIn(LB_NODE(7,i,j,k)) + 1._R8B/2._R8B*(fIn(LB_NODE(3,i,j,k))-fIn(LB_NODE(5,i,j,k))) &
                &                  - 1._R8B/2._R8B*rho(i,j,k)*u(LB_NODE(2,i,j,k)) &
                &                  + 1._R8B/6._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k)) 
            end if
           ! OUTLET:Zou/He BC
            if(prc%crd(1)== prc%np(1) -1) then
           i = lx
           fIn(LB_NODE(4,i,j,k)) = fIn(LB_NODE(2,i,j,k) )- 2._R8B/3._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k))
           fIn(LB_NODE(8,i,j,k)) = fIn(LB_NODE(6,i,j,k)) + 1._R8B/2._R8B*(fIn(LB_NODE(3,i,j,k))-fIn(LB_NODE(5,i,j,k))) &
                &                  - 1._R8B/2._R8B*rho(i,j,k)*u(LB_NODE(2,i,j,k)) &
                &                  - 1._R8B/6._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k) )
           fIn(LB_NODE(7,i,j,k)) = fIn(LB_NODE(9,i,j,k)) + 1._R8B/2._R8B*(fIn(LB_NODE(5,i,j,k))-fIn(LB_NODE(3,i,j,k))) &
                &                  + 1._R8B/2._R8B*rho(i,j,k)*u(LB_NODE(2,i,j,k)) &
                &                  - 1._R8B/6._R8B*rho(i,j,k)*u(LB_NODE(1,i,j,k)) 
            end if
        end do
      rho(:,lb_dom%lx(2),:) = 1.0d0
      u(LB_NODE(:,:,lb_dom%lx(2),:)) = 0.d0
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
               u(LB_NODE(1,i,j,k))   = sqrt(s_par%umax)
               u(LB_NODE(2,i,j,k))   = 0._R8B
               u(LB_NODE(3,i,j,k))   = sqrt(s_par%umax)
            end do
         end do
         end if

         if(prc%crd(2) == prc%np(2)-1) then
            call calc_fEq(lb_dom,rho,u,fIn,x_start,x_end,ly,ly,z_start,z_end) !top lid
         end if
#endif
#ifdef D2Q9
         k = 1
         j = ly
         if(prc%crd(2) == prc%np(2)-1) then
            do i=x_start,x_end
               ! top lid
               rho(i,j,k)   = 1._R8B
               u(LB_NODE(1,i,j,k))   = s_par%umax
               u(LB_NODE(2,i,j,k))   = 0._R8B
            end do
        call calc_fEq(lb_dom,rho,u,fIn,x_start,x_end,ly,ly,0,0) ! top wall
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
               u(LB_NODE(1,i,j,k))   = sqrt(s_par%umax)
               u(LB_NODE(2,i,j,k))   = 0._R8B
               u(LB_NODE(3,i,j,k))   = sqrt(s_par%umax)
            end do
         end do

        call calc_fEq(lb_dom,rho,u,fIn,2      ,lx-1 ,ly,ly,2      ,lz-1 ) !top lid
#endif
#ifdef D2Q9
         k = 1
         j = ly
         do i=2, lx-1      
               ! top lid
               rho(i,j,k)   = 1._R8B
               u(LB_NODE(1,i,j,k))   = s_par%umax
               u(LB_NODE(2,i,j,k))   = 0._R8B
         end do

        call calc_fEq(lb_dom,rho,u,fIn,2      ,lx-1 ,ly,ly,0,0) ! top wall
#endif



!----------------------------------------
! Default case 

     
   case (planar_standing_wave)
!        pos = (/ int(1.d0/4.d0*real(s_par%gx(1))), s_par%gx(2)/2, s_par%gx(3)/2+1    /)
!        pos = (/ int(0.25*real(s_par%gx(1))+1.), s_par%gx(2)/2, s_par%gx(3)/2+1    /)
        pos = (/ 1, s_par%gx(2)/2, s_par%gx(3)/2+1    /)
        call microphone(lb_dom,pos,fIn,s_par,prc,tstep)

   case default

   end select


   call cpu_time_measure(meas%tEnd_comm)

   meas%sbnd_duration = meas%tEnd_comm - meas%tSt_comm + meas%sbnd_duration



end  subroutine set_bnd
!------------------------------------------------------------------------




      function get_parabolic2d(r,r0)
         ! 
         ! output magnitude for vortex depending on the radius r and the vortex core radius r0 
         ! Note: this function has been error checked with gnuplot

         real(R8B)               :: get_parabolic2d
         real(R8B),intent(in)    :: r,r0 

         get_parabolic2d = -4./r0/r0*(r-r0/2.)*(r-r0/2.) + 1. 
         
      end function



      function get_parabolic3d(r,r0,t,t0)
         ! 
         ! output magnitude for vortex depending on the radius r and the vortex core radius r0 
         ! Note: this function has been error checked with gnuplot
         !
         real(R8B)               :: get_parabolic3d
         real(R8B),intent(in)    :: r,r0,t,t0 

         get_parabolic3d = -2./r0/t0*((r-r0/2.)*(r-r0/2.)+(t-t0/2.)*(t-t0/2)) + 1. 
      end function







!------------------------------------------------------------------------
subroutine bounceback(lb_dom,fOut,fIn,meas) !
! calculate bounce-back routine for the treatment of wall boundaries
!
   type(lb_block),intent(inout)  :: lb_dom
   real(R8B), intent(in)     :: fIn(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
   real(R8B), intent(inout)  :: fOut(LB_NODE(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
   type(measure )                :: meas  
   integer                 :: ii
   integer                 :: coordx(3)

   call cpu_time_measure(meas%tSt_comm) 

   do ii=1,lb_dom%nobs
      coordx(1:3) = lb_dom%obs(ii,1:3)
      fOut(LB_NODE(1,coordx(1),coordx(2),coordx(3))) = fIn(LB_NODE(opp(1),coordx(1),coordx(2),coordx(3)))
      fOut(LB_NODE(2,coordx(1),coordx(2),coordx(3))) = fIn(LB_NODE(opp(2),coordx(1),coordx(2),coordx(3)))
      fOut(LB_NODE(3,coordx(1),coordx(2),coordx(3))) = fIn(LB_NODE(opp(3),coordx(1),coordx(2),coordx(3)))
      fOut(LB_NODE(4,coordx(1),coordx(2),coordx(3))) = fIn(LB_NODE(opp(4),coordx(1),coordx(2),coordx(3)))
      fOut(LB_NODE(5,coordx(1),coordx(2),coordx(3))) = fIn(LB_NODE(opp(5),coordx(1),coordx(2),coordx(3)))
      fOut(LB_NODE(6,coordx(1),coordx(2),coordx(3))) = fIn(LB_NODE(opp(6),coordx(1),coordx(2),coordx(3)))
      fOut(LB_NODE(7,coordx(1),coordx(2),coordx(3))) = fIn(LB_NODE(opp(7),coordx(1),coordx(2),coordx(3)))
      fOut(LB_NODE(8,coordx(1),coordx(2),coordx(3))) = fIn(LB_NODE(opp(8),coordx(1),coordx(2),coordx(3)))
      fOut(LB_NODE(9,coordx(1),coordx(2),coordx(3))) = fIn(LB_NODE(opp(9),coordx(1),coordx(2),coordx(3)))
#ifdef D3Q19
      fOut(LB_NODE(10,coordx(1),coordx(2),coordx(3))) =fIn(LB_NODE(opp(10),coordx(1),coordx(2),coordx(3))) 
      fOut(LB_NODE(11,coordx(1),coordx(2),coordx(3))) =fIn(LB_NODE(opp(11),coordx(1),coordx(2),coordx(3))) 
      fOut(LB_NODE(12,coordx(1),coordx(2),coordx(3))) =fIn(LB_NODE(opp(12),coordx(1),coordx(2),coordx(3))) 
      fOut(LB_NODE(13,coordx(1),coordx(2),coordx(3))) =fIn(LB_NODE(opp(13),coordx(1),coordx(2),coordx(3))) 
      fOut(LB_NODE(14,coordx(1),coordx(2),coordx(3))) =fIn(LB_NODE(opp(14),coordx(1),coordx(2),coordx(3))) 
      fOut(LB_NODE(15,coordx(1),coordx(2),coordx(3))) =fIn(LB_NODE(opp(15),coordx(1),coordx(2),coordx(3))) 
      fOut(LB_NODE(16,coordx(1),coordx(2),coordx(3))) =fIn(LB_NODE(opp(16),coordx(1),coordx(2),coordx(3))) 
      fOut(LB_NODE(17,coordx(1),coordx(2),coordx(3))) =fIn(LB_NODE(opp(17),coordx(1),coordx(2),coordx(3))) 
      fOut(LB_NODE(18,coordx(1),coordx(2),coordx(3))) =fIn(LB_NODE(opp(18),coordx(1),coordx(2),coordx(3))) 
      fOut(LB_NODE(19,coordx(1),coordx(2),coordx(3))) =fIn(LB_NODE(opp(19),coordx(1),coordx(2),coordx(3))) 
#endif

   end do

      call cpu_time_measure(meas%tEnd_comm)
      meas%bncb_duration = meas%tEnd_comm - meas%tSt_comm + meas%bncb_duration

end subroutine bounceback
!------------------------------------------------------------------------







  subroutine treat_initial(lb_dom,s_par,tstep,prc)
!------------------------------------------------------------------------
!
! subroutine to find next and next-next neighbor for non-reflecting boundary
!

   implicit none
   type(lb_block)       :: lb_dom
   type(sim_parameter)  :: s_par
   integer, intent(in)  :: tstep 
   type(mpl_var)        :: prc  
   integer              :: i,j
   integer              :: x,y,z 
   logical              :: current_status,status_changed


   current_status = s_par%initial 
   status_changed = .false.
   if(tstep <= s_par%tInitial ) then 
      s_par%initial=.true. 
   else 
      s_par%initial=.false.
   endif 

   
if(s_par%initial .eqv. .true.) then
if(tstep ==1) then
   write(*,*) "assigning zero field..."
   if( &
#ifdef F2003_ALLOC_DERIVED_TYPE
   allocated(  &
#else
   associated(  &
#endif
lb_dom%u0)) then
      lb_dom%u0   = lb_dom%u
   else
      write(*,*) "array u0 is not allocated! Stopping..."
      stop
   endif
   write(*,*) "Done."
endif
endif

   if(s_par%initial .eqv. .false.) then 
   if(current_status .eqv. .true.) then 
      status_changed = .true.
   endif; endif

if(status_changed) then
if(prc%rk == 0) then
   write(*,*) "------------------------ "
   write(*,*) "Initial Relaxation done. "
   write(*,*) "Here, a Norm calculation should be implemented to check, if pressure field has converged"
   write(*,*) "This routine must be called before the set_bnd of the current iteration."
   write(*,*) 
   write(*,*) "Setting Zero-values for nrbc to current values..."
endif
if(lb_dom%nnr_wall > 0) then
   do i=1,lb_dom%nnr_wall   
      x = lb_dom%obs_nrwall(i,1)
      y = lb_dom%obs_nrwall(i,2)
      z = lb_dom%obs_nrwall(i,3)

      ! assign zero values
      lb_dom%nrwall_0val(i,0)      = lb_dom%rho(x,y,z)
!      lb_dom%nrwall_0val(i,1:NDIM) = lb_dom%u(LB_NODE(1:NDIM,x,y,z))
   end do
   write(*,*) prc%rk," done."
   write(*,*) 
endif
endif
   end subroutine treat_initial 
  !------------------------------------------------------------------------





end module lbm_functions


