    do k=2,lx(3)-1
      do j=2,lx(2)-1
         do i=2,lx(1)-1
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
         end do
      end do
   end do

! for walls
    do k=1,lx(3),lx(3)-1
      do j=1,lx(2),lx(2)-1
         do i=1,lx(1),lx(1)-1

            !--------------------------------------------
            ! Do Bounce Back 
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
         end do
      end do
   end do
