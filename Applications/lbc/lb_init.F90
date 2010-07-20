module lb_init      
   use mpl_lib
   use lbmodel
   use function_set
   use timing
   use lb_geo
   use tools  
#include "include/replace.h"
contains






#ifdef INIT_WITH_ROOT

subroutine init_field(fIn,u,rho,state,s_par,lb_dom,prc)
!------------------------------------------------------------------------
!
! initialization routine of initial values, done by master process only.
! values have to be distributed to slave processes
!
   implicit none
   type(sim_parameter)              :: s_par
   type(lb_block     )              :: lb_dom
   type(mpl_var      )              :: prc   
   real(R8B),dimension(s_par%gx(2)) :: y_phys
   real(R8B) :: fIn(NDX(nnod,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
   real(R8B) :: rho(s_par%gx(1),s_par%gx(2),s_par%gx(3))
   real(R8B) :: u(NDX(NDIM,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
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

#ifdef SPONGE
   lb_dom%gomega = s_par%omega
#endif

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
         rho(ii,jj,kk) = 1.0d0 + s_par%umax*sin(2.d0*pi*real(ii-1)/real(s_par%gx(1)-2)) !*cs2inv
         u(NDX(1,ii,jj,kk)) = (rho(ii,jj,kk)-1.d0)*cs 
         u(NDX(2:NDIM,ii,jj,kk)) = 0.d0 
   ! it's a standing wave, if you set u(NDX(1))=0.
      enddo
      enddo
      enddo

!----------------------------------------
! Flat plate with trailing edge

   case(plate_trail)
      u(NDX(1,1,:,:)) = s_par%umax
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
          u(NDX(1,1:xmax ,j+ymin-1,:)) = s_par%umax*get_parabolic2d(real(j,R8B),real(2*s_par%obst_r+2,R8B)) 
          u(NDX(2,1:xmax ,j+ymin-1,:)) = 0.0 
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
            u(NDX(1,1+i,s_par%gx(2)/2-5:s_par%gx(2)/2+5,:)) = 0.001/cs*(0.5+0.5*cos(2*pi*i/N+pi))
         enddo 
 




!----------------------------------------
! Cylinder in channel

    case(cylinder)
       ! Poisseuille Profile for Flow around Cylinder
      do k=1,s_par%gx(3)
        do j=1,s_par%gx(2)
          do i=1,s_par%gx(1)
#ifdef D2Q9
               u(NDX(1,i,j,k)) = s_par%umax*get_parabolic2d(real(j,R8B),real(s_par%gx(2)+1,R8B))
               if(state(i,j,k)==wall) u(NDX(1,i,j,k))=0._R8B
               u(NDX(2,i,j,k))   = 0._R8B
#endif
#ifdef D3Q19
               u(NDX(1,i,j,k)) = s_par%umax*get_parabolic3d(real(j,R8B),real(s_par%gx(2)+1,R8B), &
                                       real(k,R8B),real(s_par%gx(3)+1,R8B))
               u(NDX(2,i,j,k))   = 0._R8B
               u(NDX(3,i,j,k))   = 0._R8B
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
               u(NDX(1,i,j,k))   = 0._R8B
               u(NDX(2,i,j,k))   = 0._R8B
#ifdef D3Q19
               u(NDX(3,i,j,k))   = 0._R8B
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
               u(NDX(1,i,j,k))   = 0._R8B
               u(NDX(2,i,j,k))   = 0._R8B
#ifdef D3Q19
               u(NDX(3,i,j,k))   = 0._R8B
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
       radius=  3. !real(s_par%obst_r)
       do k=1,s_par%gx(3)
          do j=1,s_par%gx(2)
             do i=1,s_par%gx(1)
        rho(i,j,k)   = s_par%umax      &
*exp(-log(2.)*real(i-s_par%gx(1)/2)*real(i-s_par%gx(1)/2)/radius/radius)  &
*exp(-log(2.)*real(j-s_par%gx(2)/2)*real(j-s_par%gx(2)/2)/radius/radius)  &
#ifdef D3Q19
!                                 *exp(-log(2.)*real(k-s_par%gx(3)/2)*real(k-s_par%gx(3)/2)/radius/radius)  &
#endif
                                + rho(i,j,k)
!               u(NDX(1,i,j,k)) = s_par%umax
!               u(NDX(2,i,j,k)) = s_par%umax
             enddo
         enddo
      enddo

    case(gauss1d)
       radius=  3. !real(s_par%obst_r)
       do k=1,s_par%gx(3)
          do j=1,s_par%gx(2)
             do i=1,s_par%gx(1)
        rho(i,j,k)   = s_par%umax      &
*exp(-log(2.)*real(i-s_par%gx(1)/2)*real(i-s_par%gx(1)/2)/radius/radius)  &
!*exp(-log(2.)*real(j-s_par%gx(2)/2)*real(j-s_par%gx(2)/2)/radius/radius)  &
#ifdef D3Q19
!                                 *exp(-log(2.)*real(k-s_par%gx(3)/2)*real(k-s_par%gx(3)/2)/radius/radius)  &
#endif
                                + rho(i,j,k)
!               u(NDX(1,i,j,k)) = -s_par%umax
             enddo
         enddo
      enddo



!----------------------------------------
! Convected Gaussian Pulse 
 
    case(gauss_convect)
       radius=real(s_par%obst_r)
       do k=1,s_par%gx(3)  
          do j=1,s_par%gx(2)  
             do i=1,s_par%gx(1)  
          rho(i,j,k)   = s_par%umax/cs2inv*exp(-log(2.)*real(i-s_par%gx(1)/2)*real(i-s_par%gx(1)/2)/radius/radius)  &
                                   *exp(-log(2.)*real(j-s_par%gx(2)/2)*real(j-s_par%gx(2)/2)/radius/radius)  &
#ifdef D3Q19
                                   *exp(-log(2.)*real(k-s_par%gx(3)/2)*real(k-s_par%gx(3)/2)/radius/radius)  &
#endif
                     + rho(i,j,k)

               u(NDX(1:NDIM,i,j,k)) =  0.0 !s_par%umax !0.0
#ifdef D3Q19
               u(NDX(3,i,j,k)) = s_par%umax  !(rho(i,j,k) - 1.0)/cs         !s_par%umax
#else 
               u(NDX(1,i,j,k)) = s_par%umax  !(rho(i,j,k) - 1.0)/cs         !s_par%umax
#endif
             enddo
         enddo
      enddo



!----------------------------------------
! Corotating Vortex 


   case(corotating_vortex)
      ! reset complete domain to eq
      rho      = 1.0d0
      u        = 0.0d0
      
      call init_vortex_macr(u,s_par,prc,lb_dom,(/0,0,0/))


!----------------------------------------
! Default case 

   case default
      rho = 1.0d0
      u   = 0.0d0
      write(*,*) " Default case: Set complete domain to rho=1, u=0" 
   end select




  end subroutine init_field
#endif /*INIT_WITH_ROOT */






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
      real(R8B) :: u(NDX(NDIM,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
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
                  u(NDX(1,i,j,k)) = (-2.d0*dble(j-y0)/vor_r/vor_r*vor_phi)
                  u(NDX(2,i,j,k)) = ( 2.d0*dble(i-x0)/vor_r/vor_r*vor_phi)
               endif
            enddo
         enddo
      enddo


  end subroutine init_taylor_vortex
  !------------------------------------------------------------------------



  subroutine init_vortex_macr(u,s_par,prc,lb_dom,offset)
   !-----------------------------------------------------------------------
   ! Initialization routine for Corotating Vortex Pair. See Dissertation Roller 2003 or Schwartzkopf 2006, 
   ! M S Howe Theory of Vortex Sound, P. 121 
   ! Gnuplot script is located in scripts/gnuplot/ini_crvp.gpl 

      implicit none
      type(sim_parameter)  :: s_par
      type(mpl_var)        :: prc
      type(lb_block)       :: lb_dom
#ifdef INIT_WITH_ROOT
      real(R8B) :: u(NDX(NDIM,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
#else
   real(R8B) :: u(NDX(NDIM,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
#endif
      integer   :: i,j,k,kk,l
      integer   :: xmin,xmax,ymin,ymax,zmin,zmax 
      integer   :: offset(3),p1,p2,p3 
      real(R8B) :: x0,y0,z0,x1,swap,umax_tmp 

      ! corotating vortex
      real(R8B) :: vor_rc, vor_pos1(3),vor_pos2(3),vor_r
      real(R8B) :: delta_x,delta_y,delta_z
      real(R8B) :: cutoff_radius,distance,factor 
      real(R8B) :: phi0,vor_phi,pos_r,prer

      double precision MachRot, MachCore
      double precision p0, rho0, c0, gamma
      double precision u_0, u_c, omega, u_max,u_cur

      double precision Vortex_Gamma, r0, rc, alpha_Vortex
      double precision Vortex_Strength, radius, rholoc, norm_factor
         
      u_max     = 0.0        ! 
      distance    = dble(s_par%obst_r)

#ifdef INIT_WITH_ROOT
     xmin = 1; xmax=s_par%gx(1)
     ymin = 1; ymax=s_par%gx(2)
     zmin = 1; zmax=s_par%gx(3)
#else
     xmin = 1+offset(1); xmax=lb_dom%lx(1)+offset(1)
     ymin = 1+offset(2); ymax=lb_dom%lx(2)+offset(2)
     zmin = 1+offset(3); zmax=lb_dom%lx(3)+offset(3)
#endif
#ifndef TAYLOR_CRVP 
                                ! Fluid dimensionslos
                                ! -------------------

      gamma       = 1.4d0
      p0          = 1.d0*3.d0
      rho0        = 1.d0
      c0          = cs 
  
                                ! Vortex dimensionslos
                                ! --------------------
      MachRot     = sqrt(gamma)*0.08d0 ! RotationsMach-Zahl
      r0         = distance        ! r_0 = halber Wirbelabstand
      rc         = 0.3d0*r0       ! Core-Radius, Wirbelmodell 
      alpha_vortex= 1.256431d0  ! Parameter fuer Homentropic Vortex
                                ! dann max. Geschw.Betr. in r_c
      u_0         = s_par%umax       ! Ind. Geschw.Betrag in Entfernung2*r_0
      omega       = 1.d0        ! Winkelgeschw.
                                ! Homentropic Vortex
                                ! Zirkulation in Entfernung r_c
! was before      Vortex_Gamma= s_par%umax*4.d0*Pi/(1.d0-exp(-4.d0*alpha_Vortex/r_c**2)) 
      Vortex_Gamma= s_par%umax*4.d0*Pi/(1.d0-exp(-4.d0*alpha_Vortex/(rc/r0)**2)) 
      Vortex_Strength=Vortex_Gamma/2.d0/Pi
      u_c         = Vortex_Gamma/2.d0/Pi/(rc/r0) *(1.d0-exp(-alpha_Vortex))
      MachCore    = u_c/c0*MachRot ! Mach-Zahl am Radius rc

      x0    = dble(s_par%gx(1))/2.0d0 + r0 
      x1    = dble(s_par%gx(1))/2.0d0 - r0 
      y0    = dble(s_par%gx(2))/2.0d0 
      z0    = dble(s_par%gx(3))/2.0d0 + 0.5

      cutoff_radius  = 1.e37 !abs(dble(s_par%gx(1))/2.0d0 - r0)*0.95
      cutoff_radius  = abs(dble(s_par%gx(1))/2.0d0 - r0)*0.95
      

     do k=zmin,zmax
        do j=ymin,ymax
           do i=xmin,xmax
         p1 = i - offset(1)
         p2 = j - offset(2)
         p3 = k - offset(3)

!      do k=1,s_par%gx(3)
!         do j=1,s_par%gx(2)
!            do i=1,s_par%gx(1)
            radius = sqrt((dble(i)-x0)**2 + (dble(j)-y0)**2)
            factor = calc_cutoff_ratio(radius,cutoff_radius)
                                ! Homentropic Vortex 
!            u(NDX(1,i,j,k)) = s_par%umax !*sin(2.d0*Pi*dble(i/s_par%gx(1)))
!            u(NDX(2,i,j,k)) = s_par%umax !*sin(2.d0*Pi*dble(j/s_par%gx(2)))
            if(radius .ne. 0.0d0) & !.and. radius .le. cutoff_radius)  &
            u(NDX(1,p1,p2,p3))=u(NDX(1,p1,p2,p3)) + factor*Vortex_Gamma/2.d0/Pi/radius**2*r0 &
                *(1.d0-exp(-alpha_Vortex*(radius/(rc))**2)) &
                *(-1.0d0)*(-dble(j)+y0)

            radius = sqrt((dble(i)-x1)**2 + (dble(j)-y0)**2)
            factor = calc_cutoff_ratio(radius,cutoff_radius)
            if(radius .ne. 0.0d0 .and. radius .le. cutoff_radius) &
            u(NDX(1,p1,p2,p3))=u(NDX(1,p1,p2,p3)) + factor*Vortex_Gamma/2.d0/Pi/radius**2*r0 &
                *(1.d0-exp(-alpha_Vortex*(radius/(rc))**2)) &
                *(-1.0d0)*(-dble(j)+y0)

            radius = sqrt((dble(i)-x0)**2 + (dble(j)-y0)**2)
            factor = calc_cutoff_ratio(radius,cutoff_radius)
                                ! Homentropic Vortex 
            if(radius .ne. 0.0d0) & ! .and. radius .le. cutoff_radius) &
            u(NDX(2,p1,p2,p3))=u(NDX(2,p1,p2,p3)) + factor*Vortex_Gamma/2.d0/Pi/radius**2*r0 &
                *(1.d0-exp(-alpha_Vortex*(radius/(rc))**2)) &
                *(-dble(i)+x0)

            radius = sqrt((dble(i)-x1)**2 + (dble(j)-y0)**2)
            factor = calc_cutoff_ratio(radius,cutoff_radius)
            if(radius .ne. 0.0d0 .and. radius .le. cutoff_radius) &
            u(NDX(2,p1,p2,p3))=u(NDX(2,p1,p2,p3)) + factor*Vortex_Gamma/2.d0/Pi/radius**2*r0 &
                *(1.d0-exp(-alpha_Vortex*(radius/(rc))**2)) &
                *(-dble(i)+x1)

            if(u(NDX(1,p1,p2,p3)) > u_max) u_max=u(NDX(1,p1,p2,p3))
            if(u(NDX(2,p1,p2,p3)) > u_max) u_max=u(NDX(2,p1,p2,p3))
            enddo
         enddo
      enddo

#else /* TAYLOR_CRVP */
 
      phi0 = s_par%umax

      vor_r = real(s_par%obst_r)*sqrt(2.0)
      cutoff_radius = real(s_par%gx(1))/2.0d0-1.0d0

      write(*,*) 
      write(*,*) "Initializing vortices "
      write(*,*) "Cutoff Radius:   ",int(cutoff_radius)
      write(*,*) "Distance:        ",int(distance*vor_r)
      write(*,*) "Characteristic l:",s_par%obst_r
      write(*,*) "Factor:          ",s_par%obst_l
! Vortex 1

      phi0 = 20*s_par%umax

      do kk=1,2
      if(kk == 1) then
         swap =-1.d0
      else
         swap = 1.d0
      endif


      phi0 = 20*s_par%umax

      y0    = real(s_par%gx(2))/2.0d0+0.5
      z0    = real(s_par%gx(3))/2.0d0+0.5
      vor_r = real(s_par%gx(1))/20.0d0*sqrt(2.0)
      cutoff_radius = real(s_par%gx(1))/2.0d0-1.0d0
      x0    = real(s_par%gx(1))/2.0d0+0.5+swap*vor_r

      do k=zmin,zmax
         do j=ymin,ymax
            do i=xmin,xmax
               vor_phi = phi0*exp(-(dble(i-x0)*dble(i-x0) + dble(j-y0)*dble(j-y0))/(vor_r*vor_r))
               pos_r    = sqrt((dble(i)-x0 )**2 + (dble(j)-y0)**2 +(dble(k)-z0)**2)
               if(pos_r  < cutoff_radius) then
                  u(NDX(1,i,j,k)) = u(NDX(1,i,j,k)) + (-2.d0*dble(j-y0)/vor_r/vor_r*vor_phi)
                  u(NDX(2,i,j,k)) = u(NDX(2,i,j,k)) + ( 2.d0*dble(i-x0)/vor_r/vor_r*vor_phi)
               endif
            if(u(NDX(1,i,j,k)) > u_max) u_max=u(NDX(1,i,j,k))
            if(u(NDX(2,i,j,k)) > u_max) u_max=u(NDX(2,i,j,k))
            enddo
         enddo
      enddo
   enddo
#endif

      print*,'  Maximal Velocity before rescale   :', u_max 

#ifndef INIT_WITH_ROOT
   umax_tmp = u_max
   if(prc%size > 1) then
      call mpi_allreduce(umax_tmp,u_max,1,mpi_double_precision,mpi_max,prc%cart_comm,prc%ierr)
      else
   endif
#endif

#ifdef RESCALE
     factor = s_par%umax/u_max 
     u = u*factor

     u_max=0.0d0      

    ! Re-check maximum velocity
     do k=zmin,zmax
        do j=ymin,ymax
           do i=xmin,xmax
         p1 = i - offset(1)
         p2 = j - offset(2)
         p3 = k - offset(3)
           if(u(NDX(1,p1,p2,p3)) > u_max) u_max=u(NDX(1,p1,p2,p3))
           if(u(NDX(2,p1,p2,p3)) > u_max) u_max=u(NDX(2,p1,p2,p3))
           enddo
        enddo
     enddo
#endif
 print*,'  Maximal Velocity after rescale    :', u_max 
                                !
                                !---Dokumentation
                                !
  !    print*,'  Wirbelpaar CRVP:'
  !    print*,'  -----------------------------------------------------'
      print*,'  Wirbelabstand/2        r_0    [-] :', r0
      print*,'  Core Radius            rc    [-] :', rc
      print*,'  MaRot               r0        [-] :', u_0/cs
      print*,'  Max.-Mach-Zahl in  r_c MaCore [-] :', MachCore
      print*,'  Zirkulation            Gamma  [-] :', Vortex_Gamma
      print*,'  GeschwindigkeitsBetrag in 2r_0[-] :', u_0
      print*,'  GeschwindigkeitsBetrag in  r_c[-] :', u_c
      print*,'  Zeit fuer einen Wirbel-Umlauf [-] :', 2.d0*r0*Pi/u_0
      print*,'  Akustische Wellenlaenge       [-] :', Pi*c0/MachRot
      print*,'  Cutoff-Radius                 [-] :', cutoff_radius
   !   print*,'  -----------------------------------------------------'
   !   print*,'  Vortex_Strength               [-] :', Vortex_Strength
   !   print*,'  -----------------------------------------------------'
    print*,'  Vortex position one               :', s_par%gx(1)/2+r0,&
                                                    s_par%gx(2)/2+1,   s_par%gx(3)/2+1 
!print*,'  Velocity in vortex core at the beginning :', &
!u(NDX(2,s_par%gx(1)/2+distance,s_par%gx(2)/2+1,s_par%gx(3)/2+1 ))
!
!s_par%crvp_period=2.d0*Pi*distance/&
!u(NDX(2,s_par%gx(1)/2+distance,s_par%gx(2)/2+1,s_par%gx(3)/2+1))
!print*,'  Revolution time      :',s_par%crvp_period 
#ifdef CRVP
!print*,'  Start averaging after 3 periods for 2 periods...'
#else
!print*,'  -------------- Averaging is deactivated. -----------------' 
#endif
! Done initializing velocities. Pressure has to be found via the initial relaxation with fixed velocities



  contains

      function calc_cutoff_ratio(radius,cutoff_radius)
         real(R8B) :: calc_cutoff_ratio   
         real(R8B) :: pos_rel
         real(R8B),parameter :: limit=0.8
         real(R8B) :: radius,cutoff_radius
         calc_cutoff_ratio = 0.d0
         pos_rel = radius / cutoff_radius
         ! check if inside cutoff radius
         if(pos_rel < 1.0) then
            if(pos_rel < limit) then
               ! inside: 100%
               calc_cutoff_ratio = 1.0d0
            else
               ! border region: linear degression
              calc_cutoff_ratio = 1.d0-(pos_rel-limit)/(1.0d0-limit)
            endif
         endif
         calc_cutoff_ratio = -1.d0*calc_cutoff_ratio
         return
      end function calc_cutoff_ratio

      function expint(n,x)
       implicit none
      INTEGER n,MAXIT
      double precision expint,x,EPS,FPMIN,EULER
      PARAMETER (MAXIT=100,EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649d0)
      INTEGER i,ii,nm1
      double precision a,b,c,d,del,fact,h,psi
      nm1=n-1
      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
         print*,'bad arguments in expint'
         stop
      else if(n.eq.0)then
         expint=exp(-x)/x
      else if(x.eq.0.)then
         expint=1./nm1
      else if(x.gt.1.)then
         b=x+n
         c=1./FPMIN
         d=1./b
         h=d
         do 11 i=1,MAXIT
            a=-i*(nm1+i)
            b=b+2.
            d=1./(a*d+b)
            c=b+a/c
            del=c*d
            h=h*del
            if(abs(del-1.).lt.EPS)then
               expint=h*exp(-x)
               return
            endif
 11      continue
         print*,'continued fraction failed in expint'
         stop
      else
         if(nm1.ne.0)then
            expint=1./nm1
         else
            expint=-log(x)-EULER
         endif
         fact=1.
         do 13 i=1,MAXIT
            fact=-fact*x/i
            if(i.ne.nm1)then
               del=-fact/(i-nm1)
            else
               psi=-EULER
               do 12 ii=1,nm1
                  psi=psi+1./ii
 12            continue
               del=fact*(-log(x)+psi)
            endif
            expint=expint+del
            if(abs(del).lt.abs(expint)*EPS) return
 13      continue
         print*,'series failed in expint'
         stop
      endif
      return
      end function expint





  end subroutine init_vortex_macr  
  !------------------------------------------------------------------------



#ifndef INIT_WITH_ROOT
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
      integer                 :: i,j,k,l
      integer                 :: lx(3),offset(3)
      real(R8B)               :: radius           

      lx(:) = lb_dom%lx(:)
      offset(:) = lb_dom%lx(:)*prc%crd(:) 


      ! set complete domain to density 1, u=0
       do k=1,lb_dom%lx(3)
          do j=1,lb_dom%lx(2)
             do i=1,lb_dom%lx(1)
               ! In Benchmark of of Hou, Chen and Doolen (1995) rho = 2.7
               lb_dom%rho(i,j,k)   = 1.0_R8B 

               lb_dom%u(NDX(1,i,j,k))   = 0._R8B
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
#ifdef D3Q19
               lb_dom%u(NDX(3,i,j,k))   = 0._R8B
#endif
             end do
          end do
       end do


   select case(s_par%problem)


!----------------------------------------
! Corotating Vortex 


   case(corotating_vortex)
      ! reset complete domain to eq
      lb_dom%rho      = 1.0d0
      lb_dom%u        = 0.0d0
      
      
      call init_vortex_macr(lb_dom%u,s_par,prc,lb_dom,offset)


!----------------------------------------
! Gaussian Pulse 
 
    case(gaussian)
       radius=  real(s_par%obst_r)
       do k=1,lb_dom%lx(3)
          do j=1,lb_dom%lx(2)
             do i=1,lb_dom%lx(1)
        lb_dom%rho(i,j,k)   = s_par%umax      &
*exp(-log(2.)*real(i-lb_dom%lx(1)/2)*real(i-lb_dom%lx(1)/2)/radius/radius)  &
*exp(-log(2.)*real(j-lb_dom%lx(2)/2)*real(j-lb_dom%lx(2)/2)/radius/radius)  &
#ifdef D3Q19
!                                 *exp(-log(2.)*real(k-s_par%gx(3)/2)*real(k-s_par%gx(3)/2)/radius/radius)  &
#endif
                                + lb_dom%rho(i,j,k)
!               u(NDX(1,i,j,k)) = s_par%umax
!               u(NDX(2,i,j,k)) = s_par%umax
             enddo
         enddo
      enddo



!----------------------------------------
! Default case 

   case default
      lb_dom%rho = 1.0d0
      lb_dom%u   = 0.0d0
      write(*,*) " Default case: Set complete domain to rho=1, u=0" 
   end select




  end subroutine init_field_each




#ifdef OLD_INIT_EACH


!--------------------------------------------------------------------
! now assign initial conditions depending on the problem 

      if(s_par%problem == cylinder) then
       ! Poisseuille Profile for Flow around Cylinder
      do k=1,lx(3)
        do j=1,lx(2)
          do i=1,lx(1)
#ifdef D2Q9
               lb_dom%u(NDX(1,i,j,k)) = s_par%umax*get_parabolic2d(real(   &
                        j+lb_dom%lx(1)*prc%crd(1)  &
                        ,R8B),real(s_par%gx(2)+1,R8B))
               if(lb_dom%state(i,j,k)==wall) lb_dom%u(NDX(1,i,j,k))=0._R8B
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
#endif
#ifdef D3Q19
               lb_dom%u(NDX(1,i,j,k)) = s_par%umax*get_parabolic3d( &
                  real(j+lb_dom%lx(2)*prc%crd(2),R8B),real(s_par%gx(2)+1,R8B), &
                  real(k+lb_dom%lx(3)*prc%crd(3),R8B),real(s_par%gx(3)+1,R8B))
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
               lb_dom%u(NDX(3,i,j,k))   = 0._R8B
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
               lb_dom%u(NDX(1,i,j,k))   = 4._R8B*s_par%umax/real(s_par%gx(2)+1)**2*(real(j*s_par%gx(2)+1)-real(j)*real(j))
               if(lb_dom%state(i,j,k)==wall) lb_dom%u(NDX(1,i,j,k))=0._R8B
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
#endif
#ifdef D3Q19
               lb_dom%u(NDX(1,i,j,k))   = 4._R8B*s_par%umax/real(s_par%gx(2)+1)**2  &
         & *(real(j*s_par%gx(2)+1)-real(j)*real(j))*(real(k*s_par%gx(3)+1)-real(k)*real(k))
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
               lb_dom%u(NDX(3,i,j,k))   = 0._R8B
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
               lb_dom%u(NDX(1,i,j,k))   = 0._R8B
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
#ifdef D3Q19
               lb_dom%u(NDX(3,i,j,k))   = 0._R8B
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
               lb_dom%rho(i,j,k)        = 1.0_R8B 
               lb_dom%u(NDX(1,i,j,k))   = 0._R8B
               lb_dom%u(NDX(2,i,j,k))   = 0._R8B
#ifdef D3Q19
               lb_dom%u(NDX(3,i,j,k))   = 0._R8B
#endif
             end do
          end do
       end do
    else if(s_par%problem==shock) then
         write(*,*) "not implemented yet"

   end if
!----------------------------------------
! Corotating Vortex 


   case(corotating_vortex)
      ! reset complete domain to eq
      rho      = 1.0d0
      u        = 0.0d0
      
      call init_vortex_macr(u,s_par,prc,lb_dom,(/0,0,0/))

  end subroutine init_field_each
  !------------------------------------------------------------------------
#endif /* OLD_INIT_EACH */

#endif /*INIT_WITH_ROOT */


  subroutine nrbc_init(lb_dom,fIn,s_par,prc)
!------------------------------------------------------------------------
!
! subroutine to find next and next-next neighbor for non-reflecting boundary
! this is non-diagonal for the corners

   implicit none
   type(lb_block     )    :: lb_dom
   type(sim_parameter)    :: s_par
   type(mpl_var      )    :: prc  
   real(R8B),intent(in)       :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
   integer                :: i,ind
   integer                :: x,y,z ,x1,y1,z1,x2,y2,z2  

!   allocate(lb_dom%nrfplus(lb_dom%nnr_wall,nnod))
   do i=1,lb_dom%nnr_wall   
      x = lb_dom%obs_nrwall(i,1)
      y = lb_dom%obs_nrwall(i,2)
      z = lb_dom%obs_nrwall(i,3)
      x1= lb_dom%obs_nrwall(i,4)
      y1= lb_dom%obs_nrwall(i,5)
      z1= lb_dom%obs_nrwall(i,6)
      x2= lb_dom%obs_nrwall(i,7)
      y2= lb_dom%obs_nrwall(i,8)
      z2= lb_dom%obs_nrwall(i,9)
      
!  do ind=1,nnod
!     lb_dom%nrfplus(i,ind) = 0.5*(fIn(NDX(ind,x,y,z)) + fIn(NDX(opp(ind),x,y,z)))
!  enddo

   end do
   end subroutine nrbc_init
  !------------------------------------------------------------------------





  subroutine nrbc_find_neighbors(lb_dom,s_par,prc)
!------------------------------------------------------------------------
!
! subroutine to find next and next-next neighbor for non-reflecting boundary
! this is non-diagonal for the corners

   implicit none
   type(lb_block     )    :: lb_dom
   type(sim_parameter)    :: s_par
   type(mpl_var      )    :: prc  
   integer                :: i
   integer                :: x,y,z ,x1,y1,z1,x2,y2,z2,x3,y3,z3  
   integer                :: nsp 
   real(R8B)              :: outer_omega,om_sponge(0:3)

#ifdef SPONGE
   !----------------------------
   ! calculate the sponge viscosities
   outer_omega = 1.2 
   nsp = 3
   do i = 0,nsp
      om_sponge(i) = outer_omega + (s_par%omega - outer_omega)/dble(nsp+1)*(dble(i))
   enddo
#endif

   do i=1,lb_dom%nnr_wall   
      x = lb_dom%obs_nrwall(i,1)
      y = lb_dom%obs_nrwall(i,2)
      z = lb_dom%obs_nrwall(i,3)

      lb_dom%obs_nrwall(i,0) = 0  ! integer with directions (set by setting the bits)

      ! Check: Is Inlet Set?? If not, outlet. If yes, inlet
      if(btest(lb_dom%state(x,y,z),nr_wall_in) .eqv. .true.) then
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),0)
      endif

!#!define NRBC_DIAGONAL
#ifdef NRBC_DIAGONAL
      if(x == 1 .and. prc%crd(1)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,4) = x+1 
         lb_dom%obs_nrwall(i,7) = x+2
         lb_dom%obs_nrwall(i,10) = x+3
         lb_dom%obs_nrwall(i,13) = x+4
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dWEST)
      elseif(x == lb_dom%lx(1) .and. prc%crd(1)==prc%np(1)-1) then ! if on right border, use left neighbors
         lb_dom%obs_nrwall(i,4) = x-1
         lb_dom%obs_nrwall(i,7) = x-2
         lb_dom%obs_nrwall(i,10) = x-3
         lb_dom%obs_nrwall(i,13) = x-4
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dEAST)
      else
         lb_dom%obs_nrwall(i,4) = x
         lb_dom%obs_nrwall(i,7) = x
         lb_dom%obs_nrwall(i,10) = x
         lb_dom%obs_nrwall(i,13) = x
      endif

      if(y == 1 .and. prc%crd(2)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,5) = y+1 
         lb_dom%obs_nrwall(i,8) = y+2 
         lb_dom%obs_nrwall(i,11) = y+3 
         lb_dom%obs_nrwall(i,14) = y+4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dSOUTH)
      elseif(y == lb_dom%lx(2) .and. prc%crd(2)==prc%np(2)-1) then! if on right border, use left neighbors
         lb_dom%obs_nrwall(i,5) = y-1
         lb_dom%obs_nrwall(i,8) = y-2
         lb_dom%obs_nrwall(i,11) = y-3 
         lb_dom%obs_nrwall(i,14) = y-4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dNORTH)
      else
         lb_dom%obs_nrwall(i,5) = y
         lb_dom%obs_nrwall(i,8) = y
         lb_dom%obs_nrwall(i,11) = y
         lb_dom%obs_nrwall(i,14) = y 
      endif 

#ifdef D3Q19
      if(z == 1.and.  prc%crd(3)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,6) = z+1 
         lb_dom%obs_nrwall(i,9) = z+2 
         lb_dom%obs_nrwall(i,12) = z+3 
         lb_dom%obs_nrwall(i,15) = z+4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dBOTTOM)
      elseif(z == lb_dom%lx(3).and. prc%crd(3)==prc%np(3)-1 ) then
         lb_dom%obs_nrwall(i,6) = z-1
         lb_dom%obs_nrwall(i,9) = z-2
         lb_dom%obs_nrwall(i,12) = z-3 
         lb_dom%obs_nrwall(i,15) = z-4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dTOP)
      else ! if not, use same coordinate
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      endif  ! z
#else
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
#endif


#else /* NRBC_DIAGONAL */
!----------------------------
! x: check for x border
      if(x == 1 .and. prc%crd(1)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,4) = x+1 
         lb_dom%obs_nrwall(i,7) = x+2
         lb_dom%obs_nrwall(i,10) = x+3
         lb_dom%obs_nrwall(i,13) = x+4
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dWEST)
         lb_dom%obs_nrwall(i,5) = y
         lb_dom%obs_nrwall(i,8) = y
         lb_dom%obs_nrwall(i,11) = y
         lb_dom%obs_nrwall(i,14) = y 
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      elseif(x == lb_dom%lx(1) .and. prc%crd(1)==prc%np(1)-1) then ! if on right border, use left neighbors
         lb_dom%obs_nrwall(i,4) = x-1
         lb_dom%obs_nrwall(i,7) = x-2
         lb_dom%obs_nrwall(i,10) = x-3
         lb_dom%obs_nrwall(i,13) = x-4
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dEAST)
         lb_dom%obs_nrwall(i,5) = y
         lb_dom%obs_nrwall(i,8) = y
         lb_dom%obs_nrwall(i,11) = y
         lb_dom%obs_nrwall(i,14) = y 
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      else ! if not, use same coordinate

!----------------------------
! y if not x-border check y coordinate
         lb_dom%obs_nrwall(i,4) = x
         lb_dom%obs_nrwall(i,7) = x
         lb_dom%obs_nrwall(i,10) = x
         lb_dom%obs_nrwall(i,13) = x
      if(y == 1 .and. prc%crd(2)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,5) = y+1 
         lb_dom%obs_nrwall(i,8) = y+2 
         lb_dom%obs_nrwall(i,11) = y+3 
         lb_dom%obs_nrwall(i,14) = y+4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dSOUTH)
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      elseif(y == lb_dom%lx(2) .and. prc%crd(2)==prc%np(2)-1) then! if on right border, use left neighbors
         lb_dom%obs_nrwall(i,5) = y-1
         lb_dom%obs_nrwall(i,8) = y-2
         lb_dom%obs_nrwall(i,11) = y-3 
         lb_dom%obs_nrwall(i,14) = y-4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dNORTH)
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      else

!----------------------------
! z:if neither on x- nor on y-border, check z coordinate
         lb_dom%obs_nrwall(i,5) = y
         lb_dom%obs_nrwall(i,8) = y
         lb_dom%obs_nrwall(i,11) = y
         lb_dom%obs_nrwall(i,14) = y 
      if(z == 1.and.  prc%crd(3)==0) then ! if on left boarder, use right neighbors
         lb_dom%obs_nrwall(i,6) = z+1 
         lb_dom%obs_nrwall(i,9) = z+2 
         lb_dom%obs_nrwall(i,12) = z+3 
         lb_dom%obs_nrwall(i,15) = z+4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dBOTTOM)
      elseif(z == lb_dom%lx(3).and. prc%crd(3)==prc%np(3)-1 ) then
         lb_dom%obs_nrwall(i,6) = z-1
         lb_dom%obs_nrwall(i,9) = z-2
         lb_dom%obs_nrwall(i,12) = z-3 
         lb_dom%obs_nrwall(i,15) = z-4 
         lb_dom%obs_nrwall(i,0) = ibset(lb_dom%obs_nrwall(i,0),dTOP)
      else ! if not, use same coordinate
         lb_dom%obs_nrwall(i,6) = z
         lb_dom%obs_nrwall(i,9) = z
         lb_dom%obs_nrwall(i,12) = z 
         lb_dom%obs_nrwall(i,15) = z 
      endif  ! z
      endif  ! y
      endif  ! x
#endif      


!      write(31,*) x,y,z ,"  >  ",lb_dom%obs_nrwall(i,4),lb_dom%obs_nrwall(i,5),lb_dom%obs_nrwall(i,6)
      x1= lb_dom%obs_nrwall(i,4)
      y1= lb_dom%obs_nrwall(i,5)
      z1= lb_dom%obs_nrwall(i,6)
      x2= lb_dom%obs_nrwall(i,7)
      y2= lb_dom%obs_nrwall(i,8)
      z2= lb_dom%obs_nrwall(i,9)
      x3= lb_dom%obs_nrwall(i,10)
      y3= lb_dom%obs_nrwall(i,11)
      z3= lb_dom%obs_nrwall(i,12)

      ! assign zero values
      lb_dom%nrwall_0val(i,0)          = lb_dom%rho(x,y,z)
      lb_dom%nrwall_0val(i,1:NDIM    ) = lb_dom%u(NDX(1:NDIM,x,y,z))
      lb_dom%nrwall_prev(i,0)          = lb_dom%rho(x,y,z)
      lb_dom%nrwall_prev(i,1:NDIM    ) = lb_dom%u(NDX(1:NDIM,x,y,z))
      lb_dom%nrwall_prev(i,4)          = lb_dom%rho(x1,y1,z1)
      lb_dom%nrwall_prev(i,5:4+NDIM  ) = lb_dom%u(NDX(1:NDIM,x1,y1,z1))
      lb_dom%nrwall_prev(i,8)          = lb_dom%rho(x2,y2,z2)
      lb_dom%nrwall_prev(i,9:8+NDIM  ) = lb_dom%u(NDX(1:NDIM,x2,y2,z2))
#ifdef D2Q9
      lb_dom%nrwall_0val(i,3)  = 0.d0
      lb_dom%nrwall_prev(i,3)  = 0.d0
      lb_dom%nrwall_prev(i,7)  = 0.d0
      lb_dom%nrwall_prev(i,11) = 0.d0
#endif
 
#ifdef SPONGE_NRBC
   !----------------------------------------
   ! Increase viscosity on border nodes to increase stability
   
   lb_dom%omega(x, y, z ) = om_sponge(0)
   lb_dom%omega(x1,y1,z1) = om_sponge(1)
   lb_dom%omega(x2,y2,z2) = om_sponge(2)
   lb_dom%omega(x3,y3,z3) = om_sponge(3)

#endif /* SPONGE_NRBC */
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

#ifdef CRVP
      ! Communicate crvp_period
   call MPI_BCAST(s_par%crvp_period,1,MPI_INTEGER,prc%root_th,prc%cart_comm,prc%ierr)
#endif
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
               if(btest(lb_dom%state(i,j,k),nr_wall) .or. &
                  btest(lb_dom%state(i,j,k),nr_wall_in) ) then
                  if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
                    lb_dom%nnr_wall = lb_dom%nnr_wall + 1
                  end if
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
      allocate(lb_dom%obs_per(lb_dom%nper_wall,0:6))
      allocate(lb_dom%obs_nrwall(lb_dom%nnr_wall,0:15))
               ! 0     ... holds bitmap information, which border
               ! 1..3  ... coordinates of nrbc node
               ! 4..12 ... hold 1st, 2nd, 3rd neighbor node coordinates 
               ! 13.15 ... hold neighboring ghost node                  
      allocate(lb_dom%nrwall_0val(lb_dom%nnr_wall,0:3))  ! holds values of rho and u of bc
      allocate(lb_dom%nrwall_prev(lb_dom%nnr_wall,0:11))  ! holds lodi values of rho and u 
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
               if(btest(lb_dom%state(i,j,k),nr_wall) .or. &
                  btest(lb_dom%state(i,j,k),nr_wall_in)) then
                  if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
                     counter_nrwall = counter_nrwall +1
                     lb_dom%obs_nrwall(counter_nrwall,1:3) = (/i,j,k/)
                  end if
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
                  lb_dom%f_sparse((counter*nnod-1)+l) = lb_dom%fIn(NDX(l,i,j,k))
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
      real(R8B),intent(in)       :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
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
               if(btest(lb_dom%state(i,j,k),lid) .eqv. .false.) then
               uc = 0.0d0
               do l=1,nnod
                  rho_ges=rho_ges+fIn(NDX(l,i,j,k))
               !   if(fIn(NDX(l,i,j,k)) .le. 0.) write(*,*) 'found l 0  at',i,j,k,l
               enddo
#ifdef CALC_MAX_VELOCITY
                  uc(1) = cx(l,1)*fIn(NDX(l,i,j,k)) 
                  uc(2) = cx(l,2)*fIn(NDX(l,i,j,k))
                  uc(3) = cx(l,3)*fIn(NDX(l,i,j,k)) 
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
  


#ifdef INIT_WITH_ROOT  




   subroutine calc_fEq_global(rho,u,fEq,s_par)
   !------------------------------------------------------------------------
   ! same as calc_fEq, but with global limits.
   !


      implicit none
      real(R8B)                  :: usq,cs2inv,t2cs4inv,t2cs2inv,rholoc,rho0
      type(sim_parameter)        :: s_par
      real(R8B), intent(in)      :: u(NDX(NDIM,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
      real(R8B), intent(in)      :: rho(s_par%gx(1),s_par%gx(2),s_par%gx(3))
      real(R8B), intent(out)     :: fEq(NDX(nnod,s_par%gx(1),s_par%gx(2),s_par%gx(3)))
      integer               :: i,j,k,l

      cs2inv  = 1._R8B/cs**2
      t2cs4inv = 1._R8B/(2._R8B*cs**4)
      t2cs2inv = 1._R8B/(2._R8B*cs**2)
#ifdef D2Q9    
! ersetze diese rho durch lokale berechnungen und schauen, was sich aendert
      do k=1,s_par%gx(3)
         do j=1,s_par%gx(2)
            do i=1,s_par%gx(1)
             rholoc = rho(i,j,k)
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

!do l =1,nnod
!   if(abs(lb_feq(l,rholoc,(/u(NDX(1,i,j,k)),u(NDX(2,i,j,k)),0.d0/),rho0) - fEq(NDX(l,i,j,k))) .ge. 0.000001) then
!      write(*,*) 'feq error at ',i,j,k, 'pos',l
!      stop
!   endif
!enddo
            end do
         end do
      end do
#endif
#ifdef D3Q19
      do k=1,s_par%gx(3)
         do j=1,s_par%gx(2)
            do i=1,s_par%gx(1)

          usq =  (u(NDX(1,i,j,k))**2 + u(NDX(2,i,j,k))**2 + u(NDX(3,i,j,k))**2)*t2cs2inv
          rholoc   = rho(i,j,k) 
#ifdef INCOMPRESSIBLE
             rho0 = s_par%rho0
#else
             rho0 = rholoc
#endif 


           fEq(NDX(1,i,j,k)) = t(1)*(rholoc - rho0*usq) 
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


#ifdef OLD
                usq =  (u(NDX(1,i,j,k))**2 + u(NDX(2,i,j,k))**2 + u(NDX(3,i,j,k))**2)*t2cs2inv
                 fEq(NDX(1,i,j,k)) = t(1)*rho(i,j,k)*(1  &
                      - usq) 
                 fEq(NDX(2,i,j,k)) = t(2)*rho(i,j,k)*(1._R8B + (u(NDX(1,i,j,k)))*cs2inv &
                      + (u(NDX(1,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(3,i,j,k)) = t(3)*rho(i,j,k)*(1._R8B + (-u(NDX(1,i,j,k)))*cs2inv &
                      -u(NDX(1,i,j,k))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(4,i,j,k)) = t(4)*rho(i,j,k)*(1._R8B + u(NDX(2,i,j,k))*cs2inv &
                      + (u(NDX(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(5,i,j,k)) = t(5)*rho(i,j,k)*(1._R8B - u(NDX(2,i,j,k))*cs2inv &
                      + (-u(NDX(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(6,i,j,k)) = t(6)*rho(i,j,k)*(1._R8B + (u(NDX(3,i,j,k)))*cs2inv &
                      + (u(NDX(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(7,i,j,k)) = t(7)*rho(i,j,k)*(1._R8B + (-u(NDX(3,i,j,k)))*cs2inv &
                      + (-u(NDX(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(8,i,j,k)) = t(8)*rho(i,j,k)*(1._R8B + (u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))*cs2inv &
                      + (u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(9,i,j,k)) = t(9)*rho(i,j,k)*(1._R8B + (-u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))*cs2inv &
                      + (-u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(10,i,j,k)) = t(10)*rho(i,j,k)*(1._R8B + (u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))*cs2inv &
                      + (u(NDX(1,i,j,k)) -u(NDX(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(11,i,j,k)) = t(11)*rho(i,j,k)*(1._R8B + (-u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))*cs2inv &
                      + (-u(NDX(1,i,j,k)) + u(NDX(2,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(12,i,j,k)) = t(12)*rho(i,j,k)*(1._R8B + (u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                      + (u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(13,i,j,k)) = t(13)*rho(i,j,k)*(1._R8B + (-u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                      + (-u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(14,i,j,k)) = t(14)*rho(i,j,k)*(1._R8B + (u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                      + (u(NDX(1,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(15,i,j,k)) = t(15)*rho(i,j,k)*(1._R8B + (-u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                      + (-u(NDX(1,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(16,i,j,k)) = t(16)*rho(i,j,k)*(1._R8B + (u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                      + (u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(17,i,j,k)) = t(17)*rho(i,j,k)*(1._R8B + (-u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                      + (-u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(18,i,j,k)) = t(18)*rho(i,j,k)*(1._R8B + (u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))*cs2inv &
                      + (u(NDX(2,i,j,k)) -u(NDX(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
                 fEq(NDX(19,i,j,k)) = t(19)*rho(i,j,k)*(1._R8B + (-u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))*cs2inv &
                      + (-u(NDX(2,i,j,k)) + u(NDX(3,i,j,k)))**2*t2cs4inv  &
                      - usq) 
#endif
            end do
         end do
      end do
#endif    
   end subroutine calc_fEq_global
!------------------------------------------------------------------------
#endif  /*INIT_WITH_ROOT  */




   subroutine calc_fEq(lb_dom,rho,u,fEq,xmin,xmax,ymin,ymax,zmin,zmax,s_par)
   !------------------------------------------------------------------------
   !
   ! calculate the Maxwellian / equlibrium distribution function
   !

      type(sim_parameter)           :: s_par
      type(lb_block) :: lb_dom
      real(R8B)      :: usq,cs2inv,t2cs4inv,t2cs2inv,rholoc,rho0
      real(R8B)      :: u(NDX(NDIM,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      real(R8B)      :: rho(0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1))
      real(R8B)      :: fEq(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
      integer        :: i,j,k,l,err
      integer        :: xmin,xmax,ymin,ymax,zmin,zmax
      integer        :: x_l,x_u,y_l,y_u,z_l,z_u

      intent(in)     :: u,rho
      intent(inout)  :: fEq,lb_dom

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
#ifdef DEBUG_FEQ
            do l=1,nnod
             fEq(NDX(l,i,j,k)) = lb_feq(l,rho(i,j,k),(/u(NDX(1,i,j,k)),u(NDX(2,i,j,k)),0.0d0/),rho0)
            enddo
#else
             rholoc = rho(i,j,k)
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
          end do
       end do
    end do
#endif
#ifdef D3Q19
    do k=z_l,z_u
       do j=y_l,y_u
          do i=x_l,x_u

          usq =  (u(NDX(1,i,j,k))**2 + u(NDX(2,i,j,k))**2 + u(NDX(3,i,j,k))**2)*t2cs2inv
          rholoc   = rho(i,j,k) 
#ifdef INCOMPRESSIBLE
             rho0 = s_par%rho0
#else
             rho0 = rholoc
#endif 

           fEq(NDX(1,i,j,k)) = t(1)*rholoc*(1._R8B  &
                - usq) 
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
               
          end do
       end do
    end do
#endif    
  end subroutine calc_fEq
  !------------------------------------------------------------------------



#ifndef INIT_WITH_ROOT
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
      integer        :: x_start,x_end,y_start,y_end,z_start,z_end
      integer        :: lx,ly,lz,gx,gy,gz
      real   (R8B)   :: x_pos, xs,xe,radius
      integer        :: ys,ye,posbc,offset(3)
#ifdef SPONGE
      real(R8B)             :: distance,dist_max,outer_omega
      integer               :: sponge_radius,lchar
#endif /* SPONGE */
      ! Fill up the Coordinates for the Grid
      lx = lb_dom%lx(1)
      ly = lb_dom%lx(2)
      lz = lb_dom%lx(3)
      gx = s_par%gx(1)
      gy = s_par%gx(2)
      gz = s_par%gx(3)
      offset(:) = lb_dom%lx(:)*prc%crd(:) 

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
!   write(*,*) "calculation of omega for sponge layer not implemented yet."
!   stop


#endif

#ifdef SPONGE
!--------------------------------------------------------------------
! Spongelayer: assign viscosity to domain 

      lb_dom%omega = s_par%omega

         ! set the radius, from where the sponge layer starts, starting from the middle of the domain
         if(NDIM==2) then
            lchar = min(gx,gy)
         else
            lchar = min(gx,gy,gz)
         endif
         sponge_radius =int(real(lchar)*0.41d0)

         ! set the minimal omega
         outer_omega = 0.5 
         dist_max = real(lchar)*0.59d0

         ! assign the viscosity to each grid point
         do k=1,lz
            do j=1,ly
               do i=1,lx

                  distance = sqrt((real(i+offset(1))-real(gx)/2.)**2 + &
                                  (real(j+offset(2))-real(gy)/2.)**2 + &
                                  (real(k+offset(3))-real(gz)/2.)**2 ) 

                  if(distance .ge. dist_max) then
                     lb_dom%omega(i,j,k) =outer_omega 
                  elseif(distance .ge. sponge_radius .and. distance .lt.dist_max) then
                     lb_dom%omega(i,j,k) = s_par%omega + (outer_omega - s_par%omega)  &
                                       /(dist_max-sponge_radius)*(distance-sponge_radius)
                  else 
                     lb_dom%omega(i,j,k) = s_par%omega
                  endif
               enddo
            enddo
         enddo
#endif /* SPONGE_GLOBAl */




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

!----------------------------------------
! Lid driven Cavity 

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
      case(corotating_vortex)
         ! where is the lid? y has to be ly, x and z have to be 2 or 1, depending on where the prc%crd is
         posbc = nr_wall

         if(prc%crd(2) == prc%np(2) - 1) then 
             lb_dom%state(  x_start:x_end,  ly,  z_start:z_end  )  = &
         ibset(lb_dom%state(  x_start:x_end,  ly,  z_start:z_end  ),posbc)
         end if
         ! left posbc (at y- border)
!FIXME
         if(prc%crd(1) == 0)           lb_dom%state(1,:,:)     = ibset(lb_dom%state(1,:,:),posbc)
         if(prc%crd(1) == prc%np(1)-1) lb_dom%state(lx,:,:)    = ibset(lb_dom%state(lx,:,:),posbc)
         ! bottom posbc
         if(prc%crd(2) == 0)           lb_dom%state(:,1,:)     = ibset(lb_dom%state(:,1,:) ,posbc)
         if(NDIM==3) then
           ! set also front and back posbc to posbc condition
            if(prc%crd(3) == 0)        lb_dom%state(:,:,1)     = ibset(lb_dom%state(:,:,1),posbc)
            if(prc%crd(3) == prc%np(3)-1) lb_dom%state(:,:,lz) = ibset(lb_dom%state(:,:,lz),posbc )        
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

  

#else /*INIT_WITH_ROOT  */






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
      real(R8B)             :: xs,xe,ys,ye,x_pos,radius 
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

! read in geometry from file
      if(s_par%read_obs .eqv. .true.) call read_obstacles(lb_dom,s_par,prc)


#ifdef SPONGE_GLOBAL
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
         outer_omega = 0.5 

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
#endif /* SPONGE_GLOBAl */



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

       case(channel)
         ! wall for z=0 and z=lz has to be implemented
         do j=2,ly-1

            ! inlet left  
            posbc = inlet    
            lb_dom%gstate(1,j,:)    = ibset(lb_dom%gstate(1,j,:),posbc)
         
            ! outlet right  
            posbc = outlet   
            lb_dom%gstate(lx,j,:)   = ibset(lb_dom%gstate(lx,j,:),posbc) 


         end do
         ! wall boundaries
         lb_dom%gstate(1:lx,1,:)    = ibset(lb_dom%gstate(1:lx,1,:), wall)
         lb_dom%gstate(1:lx,ly,:)   = ibset(lb_dom%gstate(1:lx,ly,:),wall)

#ifdef D3Q19
         ! wall boundaries
         lb_dom%gstate(:,:,1)    = ibset(lb_dom%gstate(:,:,1),wall)
         lb_dom%gstate(:,:,lz)   = ibset(lb_dom%gstate(:,:,lz),wall)
#endif
         radius= min(s_par%gx(2)/2,s_par%gx(3)/2)
         do k=1,lz
            do j=1,ly
               if((real(j)-real(s_par%gx(2))/2.d0)**2. + (real(k)-real(s_par%gx(3))/2.d0)**2. >= real(radius)**2.) then
!               if((real(i)-real(s_par%obst_x))**2. + (real(j)-real(s_par%obst_y))**2. <= real(s_par%obst_r)**2.) then
                  lb_dom%gstate(:,j,k) = ibset(lb_dom%gstate(:,j,k),wall)
               end if
            end do
         end do



       case(cylinder)
         ! wall for z=0 and z=lz has to be implemented
         do j=2,ly-1

            ! inlet left  
            posbc = inlet    
            lb_dom%gstate(1,j,:)    = ibset(lb_dom%gstate(1,j,:),posbc)
         
            ! outlet right  
            posbc = outlet   
            lb_dom%gstate(lx,j,:)   = ibset(lb_dom%gstate(lx,j,:),posbc) 


         end do
         ! wall boundaries
         lb_dom%gstate(1:lx,1,:)    = ibset(lb_dom%gstate(1:lx,1,:), wall)
         lb_dom%gstate(1:lx,ly,:)   = ibset(lb_dom%gstate(1:lx,ly,:),wall)

#ifdef D3Q19
         ! wall boundaries
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
         posbc = nr_wall
         ! use non-reflecting boundary
         lb_dom%gstate(1,1:ly,:)   = ibset(lb_dom%gstate(1,1:ly-0,:), posbc) 
         posbc = nr_wall
         lb_dom%gstate(lx,1:ly,:)  = ibset(lb_dom%gstate(lx,1:ly,:),posbc) 
         posbc = nr_wall
         lb_dom%gstate(1:lx-0,1,:)   = ibset(lb_dom%gstate(1:lx-0,1,:), posbc) 
         posbc = nr_wall
         lb_dom%gstate(1:lx-0,ly,:)   = ibset(lb_dom%gstate(1:lx-0,ly,:), posbc) 
!         posbc = fluid 
!         lb_dom%gstate(1:2,1:2,:)  = ibset(lb_dom%gstate(1:2,1:2,:),posbc) 
!         lb_dom%gstate(lx-1:lx,1:2,:)  = ibset(lb_dom%gstate(lx-1:lx,1:2,:),posbc) 
!         lb_dom%gstate(1:2,ly:ly-1,:)  = ibset(lb_dom%gstate(1:2,ly:ly-1,:),posbc) 
!         lb_dom%gstate(lx-1:lx,ly-1:ly,:)  = ibset(lb_dom%gstate(lx-1:lx,ly-1:ly,:),posbc) 
           if(NDIM==3) then
            ! set also front and back wall to wall condition
            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), posbc)
            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),posbc)
          endif

   case(gauss1d)
         posbc = wall
         ! use non-reflecting boundary
         lb_dom%gstate(1,1:ly-0,:)   = ibset(lb_dom%gstate(1,1:ly-0,:), posbc) 
         posbc = nr_wall
         lb_dom%gstate(lx,1:ly-0,:)   = ibset(lb_dom%gstate(lx,1:ly-0,:), posbc) 
         posbc = wall
         lb_dom%gstate(1:lx-0,1,:)   = ibset(lb_dom%gstate(1:lx-0,1,:), posbc) 
         posbc = wall
         lb_dom%gstate(1:lx-0,ly,:)  = ibset(lb_dom%gstate(1:lx-0,ly,:),posbc) 
          if(NDIM==3) then
            ! set also front and back wall to wall condition
            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), posbc)
            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),posbc)
          endif

   case(gauss_convect)
         ! use non-reflecting boundary
         ! if inlet, then set separately
        posbc = per_wall 
         posbc = nr_wall
         lb_dom%gstate(1,1:ly,:)   = ibset(lb_dom%gstate(1,1:ly,:), posbc) 
         posbc = nr_wall
         lb_dom%gstate(lx,1:ly,:)  = ibset(lb_dom%gstate(lx,1:ly,:),posbc) 
!         lb_dom%gstate(1:lx,ly,:)  = ibset(lb_dom%gstate(1:lx,ly,:),posbc) 
         posbc = nr_wall
!         lb_dom%gstate(1:lx,1,:)   = ibset(lb_dom%gstate(1:lx,1,:), posbc) 
          if(NDIM==3) then
            ! set also front and back wall to wall condition
            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), posbc)
         posbc = nr_wall
            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),posbc)
          endif

!         posbc = wall
!         lb_dom%gstate(1,1:ly,:)   = ibset(lb_dom%gstate(1,1:ly,:), posbc) 
!         lb_dom%gstate(lx,1:ly,:)  = ibset(lb_dom%gstate(lx,1:ly,:),posbc) 
!         lb_dom%gstate(1:lx,ly,:)  = ibset(lb_dom%gstate(1:lx,ly,:),posbc) 
!         lb_dom%gstate(1:lx,1,:)   = ibset(lb_dom%gstate(1:lx,1,:), posbc) 
!          if(NDIM==3) then
!            ! set also front and back wall to wall condition
!            lb_dom%gstate(:,:,1)       = ibset(lb_dom%gstate(:,:,1), posbc)
!            lb_dom%gstate(:,:,lz)      = ibset(lb_dom%gstate(:,:,lz),posbc)
!          endif




!----------------------------------------
! Corotating Vortex with Spongelayer

   case(corotating_vortex)
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


  subroutine read_obstacles(lb_dom,s_par,prc)


    !------------------------------------------------------------------------
    !
    ! set up the geometry
    !

      implicit none
      type(lb_block),intent(inout)  :: lb_dom 
      type(sim_parameter),intent(inout)  :: s_par
      type(mpl_var),intent(in)     :: prc  
      character(7)                :: filename
      integer :: x,y,z,iostatus
  character*72 dummy_c

      filename = 'lbc.obs'

      open(10,file=filename)
  read(10,*) dummy_c

  ! read obstacle coordinates
    do

    read(10,*,IOSTAT=iostatus) x,y,z
    if(iostatus < 0) then
       exit
    endif

!    ! add offset
!    x = x + x_off
!    y = y + y_off
!    z = z + z_off

    ! check if obstacle inside domain boundaries, skip otherwise
    if ( x .le. s_par%gx(1) .and. x .gt. 0 .and.                                        &
         y .le. s_par%gx(2) .and. y .gt. 0 .and.                                        &
         z .le. s_par%gx(3) .and. z .gt. 0      ) then

          ! define obstacle
         lb_dom%gstate(x,y,z)  = ibset(lb_dom%gstate(x,y,z),wall) 
    else
!       write(6,*) '!!! obstacle out of range, skipped'
!       write(6,*) '!!! lx = ', x, ' , ly = ', y,' , lz = ', z
    end if

  enddo
      close(10)

   end subroutine read_obstacles
#endif  /*INIT_WITH_ROOT  */



   subroutine calc_macr_vals(lb_dom,u,rho,fIn,meas)
   !------------------------------------------------------------------------
   !
   ! calculate the macroscopic values from the microscopic distributions 
   !
     implicit none
     type(lb_block),intent(inout)  :: lb_dom
     type(measure )                :: meas   
     real(R8B), intent(out)    :: u(NDX(NDIM,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
     real(R8B), intent(out)    :: rho(0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1))
     real(R8B), intent(in)     :: fIn(NDX(nnod,0:(lb_dom%lx(1)+1),0:(lb_dom%lx(2)+1),0:(lb_dom%lx(3)+1)))
     integer              :: i,j,k,l
     
         call cpu_time_measure(meas%tSt_comm) 

     do k=1,lb_dom%lx(3)
        do j=1,lb_dom%lx(2)
           do i=1,lb_dom%lx(1)
!               if(btest(lb_dom%state(i,j,k),wall) .eqv. .false.) then
!if(lb_dom%state(i,j,k) /= wall) then
#ifdef D2Q9
               rho(i,j,k) = fIn(NDX(I__0,i,j,k)) + fIn(NDX(I__E,i,j,k))  + &
                            fIn(NDX(I__W,i,j,k)) + fIn(NDX(I__N,i,j,k))  + &
                            fIn(NDX(I__S,i,j,k)) + fIn(NDX(I_NE,i,j,k))  + &
                            fIn(NDX(I_SW,i,j,k)) + fIn(NDX(I_SE,i,j,k))  + &
                            fIn(NDX(I_NW,i,j,k))                 
                                                           
         u(NDX(1,i,j,k)) =  fIn(NDX(I__E,i,j,k)) - fIn(NDX(I__W,i,j,k))  + &
                            fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k))  + &
                            fIn(NDX(I_SE,i,j,k)) - fIn(NDX(I_NW,i,j,k))    
         u(NDX(2,i,j,k)) =  fIn(NDX(I__N,i,j,k)) - fIn(NDX(I__S,i,j,k))  + &
                            fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k))  + &  
                            fIn(NDX(I_NW,i,j,k)) - fIn(NDX(I_SE,i,j,k))       
#endif
#ifdef D3Q19
       rho(i,j,k)  =  fIn(NDX(I__0,i,j,k)) + &
                      fIn(NDX(I__E,i,j,k)) + fIn(NDX(I__W,i,j,k)) + &
                      fIn(NDX(I__N,i,j,k)) + fIn(NDX(I__S,i,j,k)) + &
                      fIn(NDX(I__T,i,j,k)) + fIn(NDX(I__B,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) + fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_SE,i,j,k)) + fIn(NDX(I_NW,i,j,k)) + &
                      fIn(NDX(I_TE,i,j,k)) + fIn(NDX(I_BW,i,j,k)) + &
                      fIn(NDX(I_BE,i,j,k)) + fIn(NDX(I_TW,i,j,k)) + &
                      fIn(NDX(I_TN,i,j,k)) + fIn(NDX(I_BS,i,j,k)) + &
                      fIn(NDX(I_BN,i,j,k)) + fIn(NDX(I_TS,i,j,k))
    u(NDX(1,i,j,k)) = fIn(NDX(I__E,i,j,k)) - fIn(NDX(I__W,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_SE,i,j,k)) - fIn(NDX(I_NW,i,j,k)) + &
                      fIn(NDX(I_TE,i,j,k)) - fIn(NDX(I_BW,i,j,k)) + &
                      fIn(NDX(I_BE,i,j,k)) - fIn(NDX(I_TW,i,j,k))  
    u(NDX(2,i,j,k)) = fIn(NDX(I__N,i,j,k)) - fIn(NDX(I__S,i,j,k)) + &
                      fIn(NDX(I_NE,i,j,k)) - fIn(NDX(I_SW,i,j,k)) + &
                      fIn(NDX(I_NW,i,j,k)) - fIn(NDX(I_SE,i,j,k)) + & 
                      fIn(NDX(I_TN,i,j,k)) - fIn(NDX(I_BS,i,j,k)) + &
                      fIn(NDX(I_BN,i,j,k)) - fIn(NDX(I_TS,i,j,k))
    u(NDX(3,i,j,k)) = fIn(NDX(I__T,i,j,k)) - fIn(NDX(I__B,i,j,k))  + &
                      fIn(NDX(I_TE,i,j,k)) - fIn(NDX(I_BW,i,j,k))  + &
                      fIn(NDX(I_TW,i,j,k)) - fIn(NDX(I_BE,i,j,k))  + &
                      fIn(NDX(I_TN,i,j,k)) - fIn(NDX(I_BS,i,j,k))  + &
                      fIn(NDX(I_TS,i,j,k)) - fIn(NDX(I_BN,i,j,k))


#endif
               u(NDX(:,i,j,k))=u(NDX(:,i,j,k))/rho(i,j,k)
!endif

            enddo
         enddo
      enddo
         call cpu_time_measure(meas%tEnd_comm)
         meas%ccmv_duration = meas%tEnd_comm - meas%tSt_comm + meas%ccmv_duration
   end subroutine calc_macr_vals
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
      lb_dom%nrwall_0val(i,1:NDIM) = lb_dom%u(NDX(1:NDIM,x,y,z))
   end do
   write(*,*) prc%rk," done."
   write(*,*) 
endif
endif
   end subroutine treat_initial 
  !------------------------------------------------------------------------





end module lb_init     


