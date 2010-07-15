c-----------------------------------------------------------------------
cTime-stamp: <03/07/19 23:50:30 roller@linux.local>
c-----------------------------------------------------------------------
        subroutine ini_vortex_crvp
!
!          Sprache: erweitertes FORTRAN77
!
!           Inhalt: Setzt die Anfangsbedingungen fuer 
!                   Co-Rotating Vortex Pair
!
!      Subroutinen: dieses File enthaelt die Subroutinen       
!                   ini_vortex_crvp:     Co-Rotating-Vortex-Pair
!                   FUNCTION expint(n,x)
!
!                   Mittels C-Preprocessor-defines werden die Faelle
!                   unterschieden, und der jeweils 'richtige' Quellcode
!                   erzeugt.  
!                   Diese Datei muss deshalb vor uebersetzen mit einem 
!                   FORTRAN-Compiler mit dem C-Preprocessor bearbeitet 
!                   werden
!
! 
c-----------------------------------------------------------------------

#include "definitions.h"
       implicit none

        include 'parameter.h'
        include 'grid.h'
        include 'stroemung.h'
        include 'mpvdruck.h'
        include 'rand_werte.h'
        include 'abgeleitet.h'

        

                                !
                                ! Lokale Variablen, die immer gebraucht
                                ! werden 
                                !
      integer i, j
      double precision Pi, MachRot, MachCore
      double precision p0, rho0, c0, gamma
      double precision u_0, u_c, omega
      
      double precision Vortex_Gamma, r_0, r_c, alpha_Vortex
      double precision Vortex_Strength, radius
      double precision expint
      double precision x_z, y_z
      
      Pi          = 3.14159265d0




#ifdef HINDERNIS
        write(*,*) '----------------- FEHLER !'
        write(*,*) 'Es wurde in "definitions.h" als Anfangsbedingung '
        write(*,*) '"VORTEX" gewaehlt.'
        write(*,*) 'Es ist aber nicht moeglich, dieses Beispiel mit'
        write(*,*) 'einem Hindernis zu rechnen.'
        write(*,*) 'Korrigieren Sie den Eintrag in "definitions.h"'
        write(*,*) ' ---------------- STOP in subroutine startup '
        stop
#endif /* HINDERNIS */
                                !
                                ! Definition der Wandgeschwindigkeiten
                                ! und Geschwindigkeitsprofile am Einlass
                                !
                                ! Setzen der Wandgeschwindigkeiten
                                ! oben/unten 
                                !
        uwallo = 0.d0
        uwallu = 0.d0
                                ! Setzen der Wandgeschwindigkeiten
                                ! rechts/links
                                !
        vwallr = 0.d0
        vwalll = 0.d0
                                !
                                ! Setzen der Bedingungen am Ein-/Auslass
                                ! 
                                ! links/rechts:
        do j = 0, npy
           uinl(j)   = 1.d0
           vinl(j)   = 0.d0
           rhoinl(j) = 1.d0
           poutr(j)  = 1.d0
        enddo
                                !
                                ! in y-Richtung laeuft v von -1:npy
                                !
        vinl(-1) = 0.d0
                                !
                                ! Anfangsbedingungen setzen:
                                ! Inkompressibel, deshalb Dichte und
                                ! Druck konstant, normiert auf 1
                                !
      
      
                                ! Fluid dimensionslos
                                ! -------------------
      gamma       = 1.4d0
      p0          = 1.d0
      rho0        = 1.d0
      c0          = sqrt(gamma)
      
      
                                ! Vortex dimensionslos
                                ! --------------------
      MachRot     = sqrt(gamma)*0.08d0 ! RotationsMach-Zahl
      r_0         = 1.d0        ! r_0 = halber Wirbelabstand
      r_c         = 0.3d0       ! Core-Radius, Wirbelmodell 
      alpha_vortex= 1.256431d0  ! Parameter fuer Homentropic Vortex
                                ! dann max. Geschw.Betr. in r_c
      u_0         = 1.d0        ! Ind. Geschw.Betrag in Entfernung2*r_0
      omega       = 1.d0        ! Winkelgeschw.
                                ! Homentropic Vortex
                                ! Zirkulation in Entfernung r_c
      Vortex_Gamma= 4.d0*Pi/(1.d0-exp(-4.d0*alpha_Vortex/r_c**2)) 
      Vortex_Strength=Vortex_Gamma/2.d0/Pi
      u_c         = Vortex_Gamma/2.d0/Pi/r_c *(1.d0-exp(-alpha_Vortex))
      MachCore    = u_c/c0*MachRot ! Mach-Zahl am Radius r_c
      
                                !
                                !---Dokumentation
                                !
      print*,'  Wirbelpaar CRVP:'
      print*,'  -----------------------------------------------------'
      print*,'  Wirbelabstand/2        r_0    [-] :', r_0
      print*,'  Core Radius            r_c    [-] :', r_c
      print*,'  Rot.-Mach-Zahl in 2r_0 MaRot  [-] :', MachRot
      print*,'  Max.-Mach-Zahl in  r_c MaCore [-] :', MachCore
      print*,'  Zirkulation            Gamma  [-] :', Vortex_Gamma
      print*,'  GeschwindigkeitsBetrag in 2r_0[-] :', u_0
      print*,'  GeschwindigkeitsBetrag in  r_c[-] :', u_c
      print*,'  Zeit fuer einen Wirbel-Umlauf [-] :', 2.d0*Pi
      print*,'  Akustische Wellenlaenge       [-] :', Pi*c0/MachRot
      print*,'  -----------------------------------------------------'
      print*,'  Vortex_Strength               [-] :', Vortex_Strength
      print*,'  -----------------------------------------------------'
      
      
      do j =  0, nraumy
         do i= 0, nraumx  
                                !
                                ! Fall: Randdefinierte Groessen
                                !       Geschw. u
                                !
                                ! 1. Wirbel mit Wirbelzentrum
            x_z= 1.d0
            y_z= 0.d0
                                ! Koordinate 
            radius = sqrt((x(i)-x_z+dx(i)/2.0d0)**2 + (y(j)-y_z)**2)
                                ! Homentropic Vortex 
            u(i,j)=Vortex_Gamma/2.d0/Pi/radius**2
     &           *(1.d0-exp(-alpha_Vortex*(radius/r_c)**2))
     &           *(-y(j)+y_z)
                                ! 2. Wirbel mit Wirbelzentrum
            x_z= -1.d0
            y_z= 0.d0
                                ! Koordinate 
            radius = sqrt((x(i)-x_z+dx(i)/2.0d0)**2 + (y(j)-y_z)**2)
                                ! Homentropic Vortex 
            u(i,j)= u(i,j) + Vortex_Gamma/2.d0/Pi/radius**2
     &                     *(1.d0-exp(-alpha_Vortex*(radius/r_c)**2))
     &                     *(-y(j)+y_z)            
                                !
                                ! Fall: Randdefinierte Groessen
                                !       Geschw. v
                                !
                                ! 1. Wirbel mit Wirbelzentrum
            x_z= 1.d0
            y_z= 0.d0
                                ! Koordinate
            radius = sqrt((x(i)-x_z)**2 + (y(j)-y_z+dx(i)/2.0d0)**2)
                                ! Homentropic Vortex 
            v(i,j)=Vortex_Gamma/2.d0/Pi/radius**2
     &               *(1.d0-exp(-alpha_Vortex*(radius/r_c)**2))
     &               *(x(i)-x_z)          
                                ! 2. Wirbel mit Wirbelzentrum
            x_z= -1.d0
            y_z= 0.d0
                                ! Koordinate
            radius = sqrt((x(i)-x_z)**2 + (y(j)-y_z+dx(i)/2.0d0)**2)
                                ! Homentropic Vortex 
            v(i,j)=v(i,j) + Vortex_Gamma/2.d0/Pi/radius**2
     &               *(1.d0-exp(-alpha_Vortex*(radius/r_c)**2))
     &               *(x(i)-x_z)    
                                !
                                ! Fall: Zellzentrierte Groesse
                                !
                                !
                                !       Druck und Dichte
                                !
                                ! 1. Wirbel mit Wirbelzentrum
            x_z= 1.d0
            y_z= 0.d0
            radius = sqrt((x(i)-x_z)**2 + (y(j)-y_z)**2)     
            if (radius .eq. 0.0) then 
               print*,'radius =0'
            endif
            p(i,j)=rho0*c0**2/gamma*
     &               (
     &               1.d0-MachRot**2*(gamma-1.d0)*Vortex_Gamma**2/
     &               (4.d0*Pi**2*c0**2*radius**2)*
     &                (
     &           0.5d0-exp(-alpha_Vortex*(radius/r_c)**2)+
     &           0.5d0*exp(-2.d0*alpha_Vortex*(radius/r_c)**2)-
     &           alpha_Vortex*(radius/r_c)**2*
     &           expint(1,2.d0*alpha_Vortex*(radius/r_c)**2)+
     &           alpha_Vortex*(radius/r_c)**2*
     &           expint(1,alpha_Vortex*(radius/r_c)**2)
     &                 )
     &                )**(gamma/(gamma-1.d0))
                                ! 2. Wirbel mit Wirbelzentrum: 
                                ! ACHTUNG Behandlung Wirbel 1,2
            x_z= -1.d0
            y_z= 0.d0
            radius = sqrt((x(i)-x_z)**2 + (y(j)-y_z)**2)     
            p(i,j)= p(i,j) + rho0*c0**2/gamma*
     &               (
     &               1.d0-MachRot**2*(gamma-1.d0)*Vortex_Gamma**2/
     &               (4.d0*Pi**2*c0**2*radius**2)*
     &                (
     &           0.5d0-exp(-alpha_Vortex*(radius/r_c)**2)+
     &           0.5d0*exp(-2.d0*alpha_Vortex*(radius/r_c)**2)-
     &           alpha_Vortex*(radius/r_c)**2*
     &           expint(1,2.d0*alpha_Vortex*(radius/r_c)**2)+
     &           alpha_Vortex*(radius/r_c)**2*
     &           expint(1,alpha_Vortex*(radius/r_c)**2)
     &                )
     &                )**(gamma/(gamma-1.d0)) - rho0*c0**2/gamma
            
            
            
         enddo
      enddo  

      pnull= 1.d0
      peins= 0.d0
      pzwei= (p-pnull)/MachRot**2
      rho= p**(1.d0/gamma) 
      
      return
      end
      


      

      
      
c-----------------------------------------------------------------------
c     
      FUNCTION expint(n,x)
c     
c-----------------------------------------------------------------------
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
      END
      
