// Boundary conditions are computed after the collision step and
// before the propagation one
// Corrected!
// (IMPORTANT: in the PRE paper wave amplitudes are called: L1, L4, L2,
// (L4 but in the code these are: L1, L2, L4, L5
// ((ie. when a L4 appears in the source code it should be replaced by
// (L2 in order to have the same nomenclature used in the paper).

IF_BOUNDARY_U_IN{
 G=1.0;            // \kappa in the article
 Cs=C*sqrt(G); // modified speed of sound
 vy=0.0;           // v velocity component
 u_=0.1;          // u velocity component
 sigma_5=0.0; // to be used for modeling of the inlet amplitude wave
 //
 // first order forward stencil is used here instead of Eq.(21)
 //
 drhodx=  -(rho_prev[y][x+1] - rho_prev[y][x+2]);
 dudx  =  -(uprev[y][x+1] - uprev[y][x+2]) ;
 //
 // Extrapolation of the velocity value to the boundary
 //
 uw=1.5*ufield[y][x+1]-0.5*ufield[y][x+2];
 //
 // wave amplitudes  Eq.(9)
 //
 L5=sigma_5*(uw-u_0);
 L1=(uprev[y][x]-Cs)*(Cs*Cs*drhodx - rho_0*Cs*dudx);
 L4=uprev[y][x]*(Cs*Cs*drhodx - Cs*Cs*drhodx/G);
 //
 // new rho and vx from LODI equations; Eq.(10) with a first order
time derivative stencil as in Eq.(18)
 //
 d_rho=rho_prev[y][x]-(L4+0.5*(L5+L1))/(Cs*Cs);
 vx   = uprev[y][x]-(L5-L1)/(2.0*rho_0*Cs);
 //
 // macroscopic variables at this t to be used the next time step
 //
 rho_prev[y][x]  =d_rho;
 rho_prev[y][x+1]=f[0][y][x+1]+f[1][y][x+1]+f[2][y][x+1]+f[3][y][x+1]+f[4][y][x+1]+f[5][y][x+1]+f[6][y][x+1]+f[7][y][x+1]+f[8][y][x+1];
 rho_prev[y][x+2]=f[0][y][x+2]+f[1][y][x+2]+f[2][y][x+2]+f[3][y][x+2]+f[4][y][x+2]+f[5][y][x+2]+f[6][y][x+2]+f[7][y][x+2]+f[8][y][x+2];
 uprev[y][x]=vx;
 uprev[y][x+1]=f[1][y][x+1]+f[5][y][x+1]+f[8][y][x+1]-f[3][y][x+1]-f[6][y][x+1]-f[7][y][x+1];
 uprev[y][x+2]=f[1][y][x+2]+f[5][y][x+2]+f[8][y][x+2]-f[3][y][x+2]-
f[6][y][x+2]-f[7][y][x+2];
 //
 // BC for f's using macroscopic variables computed using LODI; Eq.(19)
 //
 f[1][y][x]=f[3][y][x+1]   + 0.666666666666666666666667*vx*d_rho;
 f[5][y][x]=f[7][y+1][x+1] + 0.166666666666666666666667*vx*d_rho;
 f[8][y][x]=f[6][y-1][x+1] + 0.166666666666666666666667*vx*d_rho;
}

//##############################################################################################################

IF_BOUNDARY_RHO_OUT{
 rho_0=1.0;      // boundary condition for rho
 d_rho=rho_0;  //
 G=1.0;           // \kappa in the paper
 K=0.0;           // to be used for modeling inlet amplitude wave
 Cs=C*sqrt(G);// modified speed of sound
 //
 // derivatives; Eq.(15) !!!!!!!!!!!!!!!!
 //
 drhodx= 0.33333333333*(8.0*rho_prev[y][x] - 9.0*rho_prev[y][x-1] +
rho_prev[y][x-2]);
 dudx  = 0.33333333333*(8.0*uprev[y][x] - 9.0*uprev[y][x-1] + uprev[y][x-2]) ;
 dvdx  = 0.33333333333*(8.0*vprev[y][x] - 9.0*vprev[y][x-1] + vprev[y][x-2]);
 //
 // wave amplitudes Eq.(9)
 //
 L1=K*Cs*Cs*(rho_prev[y][x]-rho_0);//rho_prev[y][x]
 L2=uprev[y][x]*dvdx;
 L4=uprev[y][x]*(Cs*Cs*drhodx - Cs*Cs*drhodx/G);
 L5=(uprev[y][x]+Cs)*(Cs*Cs*drhodx/G + rho_prev[y][x]*Cs*dudx);
 //
 // LODI equations; Eq(18)
 //
 d_rho=rho_prev[y][x]-(L4+0.5*(L5+L1))/(Cs*Cs);
 vx   = uprev[y][x]-(L5-L1)/(2.0*rho_0*Cs);
 vy   = vprev[y][x]-L2;
 //
 // storing macroscopic variables at this t to be used the next time step
 //
 rho_prev[y][x]  =d_rho;
 rho_prev[y][x-1]=f[0][y][x-1]+f[1][y][x-1]+f[2][y][x-1]+f[3][y][x-1]+f[4][y][x-1]+f[5][y][x-1]+f[6][y][x-1]+f[7][y][x-1]+f[8][y][x-1];
 rho_prev[y][x-2]=f[0][y][x-2]+f[1][y][x-2]+f[2][y][x-2]+f[3][y][x-2]+f[4][y][x-2]+
f[5][y][x-2]+f[6][y][x-2]+f[7][y][x-2]+f[8][y][x-2];
 uprev[y][x]=vx;
 uprev[y][x-1]=f[1][y][x-1]+f[5][y][x-1]+f[8][y][x-1]-f[3][y][x-1]-f[6][y][x-1]-f[7][y][x-1];
 uprev[y][x-2]=f[1][y][x-2]+f[5][y][x-2]+f[8][y][x-2]-f[3][y][x-2]-f[6][y][x-2]-f[7][y][x-2];
 vprev[y][x]=vy;
 vprev[y][x-1]=f[2][y][x-1]+f[5][y][x-1]+f[6][y][x-1]-f[4][y][x-1]-f[7][y][x-1]-f[8][y][x-1];
 vprev[y][x-2]=f[2][y][x-2]+f[5][y][x-2]+f[6][y][x-2]-f[4][y][x-2]-f[7][y][x-2]-f[8][y][x-2];
 //
 // Eq.(11)
 //
 square =1.5*(vx * vx +vy *vy );
 f[3][y][x]=-f[1][y][x-1]   + 0.2222222222222222222*(d_rho +
4.5*vx*vx - square) + (2.0-s8)*(fplus[1][y]-eplus[1][y]);
 f[6][y][x]=-f[8][y+1][x-1] + 0.0555555555555555556*(d_rho + 4.5
*(vx-vy) * (vx-vy) - square)+(2.0-s8)*(fplus[8][y]-eplus[8][y]);
 f[7][y][x]=-f[5][y-1][x-1] + 0.0555555555555555556*(d_rho + 4.5
*(vx+vy)*(vx+vy) - square)+(2.0-s8)*(fplus[5][y]-eplus[5][y]);
}

