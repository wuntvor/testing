module lb_geo       
   use lbmodel
#include "include/replace.h"
contains


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
         real(R8B)               :: get_parabolic3d,distance,radius
         real(R8B)               :: t0_2,r0_2,r0sq,t0sq 
         real(R8B),intent(in)    :: r,r0,t,t0 

         r0_2 = r0/2.d0
         t0_2 = t0/2.d0
         r0_sq= (r-r0_2)**2
         t0_sq= (t-t0_2)**2

         distance = r0_sq + t0_sq 
         radius   = min(r0_2**2,t0_2**2)
 
         if(distance < radius) then
!             get_parabolic3d = -2./r0/t0*((r-r0/2.)*(r-r0/2.)+(t-t0/2.)*(t-t0/2)) + 1. 
            get_parabolic3d = - r0_sq/(r0_2)**2 - t0_sq/(t0_2)**2 + 1.d0 
         else
            get_parabolic3d = 0.0d0
         endif
!         get_ux = u_max*16.d0 / ( ((gy-2.d0)*(gy-2.d0))*((gz-2.d0)*(gz-2.d0)) )&
!                * (y-1.5d0)*((gy-0.5d0)-y) * (z-1.5d0)*((gz-0.5d0)-z)
      end function





end module lb_geo      


