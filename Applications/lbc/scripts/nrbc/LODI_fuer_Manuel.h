//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//

f     = distributions->getStartAdressOfSortedArray(x1,x2,x3,0);
ftemp = tempdistributions->getStartAdressOfSortedArray(x1,x2,x3,0);

//////////////////////////////////////////////////////////////////////////
// LODI
//////////////////////////////////////////////////////////////////////////
nx1 = nnx1 = x1;
nx2 = nnx2 = x2;
nx3 = nnx3 = x3;

double Ma,swap;
calcMacroVals( f, rho, vx1, vx2, vx3 );

// determine local Mach number and "mirror" derivates if neccessary (swap = +- 1)
// East border heisst, Rand zeigt in Ost-Richtung. genau umgekehrt wie bei mir.
if     ( bc->hasDensityBoundaryFlag(D3Q19System::E) ) {nx1 += 1; nnx1 += 2; Ma = vx1; swap = -1.0; }
else if( bc->hasDensityBoundaryFlag(D3Q19System::W) ) {nx1 -= 1; nnx1 -= 2; Ma = vx1; swap =  1.0; }
else if( bc->hasDensityBoundaryFlag(D3Q19System::N) ) {nx2 += 1; nnx2 += 2; Ma = vx2; swap = -1.0; }
else if( bc->hasDensityBoundaryFlag(D3Q19System::S) ) {nx2 -= 1; nnx2 -= 2; Ma = vx2; swap =  1.0; }
else if( bc->hasDensityBoundaryFlag(D3Q19System::T) ) {nx3 += 1; nnx3 += 2; Ma = vx3; swap = -1.0; }
else if( bc->hasDensityBoundaryFlag(D3Q19System::B) ) {nx3 -= 1; nnx3 -= 2; Ma = vx3; swap =  1.0; }
else throw UbException(UB_EXARGS,"Danger...kein orthogonales BC-Flag am Dichterand...");

if( nnx1<0 || nnx1>=matrixLX1 ) throw UbException(UB_EXARGS,"nnx1<0 || nnx1>=lengthX1");
if( nnx2<0 || nnx2>=matrixLX2 ) throw UbException(UB_EXARGS,"nnx2<0 || nnx2>=lengthX2");
if( nnx3<0 || nnx3>=matrixLX3 ) throw UbException(UB_EXARGS,"nnx3<0 || nnx3>=lengthX3");

const double c_s    = std::sqrt( 1./3. );
Ma = std::fabs(Ma) / c_s;

// Parameter k1 bestimmen
const double sigma1 = 1.0; //chose!
const double k1 = swap * sigma1 * ( 1.0 - Ma*Ma ) * c_s / bc->getDensityLodiLength();

// makroskopische groessen an den nachbarpunkten berechnen
double nrho,  nvx1,  nvx2,  nvx3;
double nnrho, nnvx1, nnvx2, nnvx3;
// Das hier muessten doch post-collision verteilungen sein.
calcMacroVals( distributions->getStartAdressOfSortedArray(nx1 ,nx2 ,nx3 ,0), nrho,  nvx1,  nvx2,  nvx3  );
calcMacroVals( distributions->getStartAdressOfSortedArray(nnx1,nnx2,nnx3,0), nnrho, nnvx1, nnvx2, nnvx3 );

//derivatives
const double dpdx  = swap * UbMath::c1o3 * c_s*c_s* ( 8.0*bc->getDensityLodiDensity()    - 9.0*nrho + 1.0*nnrho );
const double duxdx = swap * UbMath::c1o3 *          ( 8.0*bc->getDensityLodiVelocityX1() - 9.0*nvx1 + 1.0*nnvx1 );
const double duydx = swap * UbMath::c1o3 *          ( 8.0*bc->getDensityLodiVelocityX2() - 9.0*nvx2 + 1.0*nnvx2 );
const double duzdx = swap * UbMath::c1o3 *          ( 8.0*bc->getDensityLodiVelocityX3() - 9.0*nvx3 + 1.0*nnvx3 );

//L2,L4,L5
// Kann ich hier nicht einfach anstatt densityLodi die Dichte aus den Post-Collisions nehmen?
// Wie setze ich anfangs rho_lodi?? 
const double L1 = k1*( bc->getDensityLodiDensity() - bc->getBoundaryDensity() )*c_s*c_s;

double L2,L4,L5;

if(  bc->hasDensityBoundaryFlag( D3Q19System::W )||bc->hasDensityBoundaryFlag( D3Q19System::E ) )
{
   //y,z,x
   L2 = duydx*bc->densityLodiVelocityX1();
   L5 = duzdx*bc->densityLodiVelocityX1();
   L4 = ( bc->densityLodiVelocityX1() + (swap*c_s) )*( dpdx + ( 1.0 + bc->getBoundaryDensity() )*(swap*c_s)*duxdx );

   bc->densityLodiDensity()    -= (float)(( L4+L1 ) / ( 2.0*c_s*c_s ));

   bc->densityLodiVelocityX1() -= (float)(( L4-L1 ) / ( 2.0*( 1.0+bc->getBoundaryDensity() )*(swap*c_s) ));
   bc->densityLodiVelocityX2() -= (float)L2;
   bc->densityLodiVelocityX3() -= (float)L5;
}
else if(  bc->hasDensityBoundaryFlag( D3Q19System::S )||bc->hasDensityBoundaryFlag( D3Q19System::N ) )
{
   //z,x,y
   L2 = duzdx*bc->densityLodiVelocityX2();
   L5 = duxdx*bc->densityLodiVelocityX2();
   L4 = ( bc->densityLodiVelocityX2() + (swap*c_s) )*( dpdx + ( 1.0 + bc->getBoundaryDensity() )*(swap*c_s)*duydx );

   bc->densityLodiDensity()    -= (float)(( L4+L1 ) / ( 2.0*c_s*c_s ));

   bc->densityLodiVelocityX2() -= (float)(( L4-L1 ) / ( 2.0*( 1.0+bc->getBoundaryDensity() )*(swap*c_s) ));
   bc->densityLodiVelocityX3() -= (float)L2;
   bc->densityLodiVelocityX1() -= (float)L5;
}
else if(  bc->hasDensityBoundaryFlag( D3Q19System::B )||bc->hasDensityBoundaryFlag( D3Q19System::T ) )
{
   //x,y,z
   L2 = duxdx*bc->densityLodiVelocityX3();
   L5 = duydx*bc->densityLodiVelocityX3();
   L4 = ( bc->densityLodiVelocityX3() + (swap*c_s) )*( dpdx + ( 1.0 + bc->getBoundaryDensity() )*(swap*c_s)*duzdx );

   bc->densityLodiDensity()    -= (float)(( L4+L1 ) / ( 2.0*c_s*c_s ));

   bc->densityLodiVelocityX3() -= (float)(( L4-L1 ) / ( 2.0*( 1.0+bc->getBoundaryDensity() )*(swap*c_s) ));
   bc->densityLodiVelocityX1() -= (float)L2;
   bc->densityLodiVelocityX2() -= (float)L5;
}


// Nach dem Berechnen der LODI-Größen wird damit dann Anti-BB gemacht
