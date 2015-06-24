/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Shearing Box setup.

  Set initial conditions for a shearing box model in two or 
  three dimensions:
  \f[ 
    \rho = \left\{\begin{array}{ll}
      1              & \quad{\rm(unstratified)} \\ \noalign{\medskip}
     \exp(-z^2/H^2) & \quad{\rm(stratified)}
   \end{array}\right.
    ;\quad 
    v_y = -q\Omega x\,; \quad
    p = \rho c_s^2\,; \quad
    \vec{B} = \left\{\begin{array}{ll}
     (0,0,B_z)          & \quad{\rm(net\,flux)} \\ \noalign{\medskip}
     (0,0,B_z\sin 2\pi x/L_x) &  \quad{\rm(no\,net\,flux)} \\ \noalign{\medskip}
    \end{array}\right.
   \f]
  where, by default, the angular rotation frequency is \f$\Omega = 1\f$ 
  while the shear parameter is \f$ q = 3/2 \f$ (Keplerian rotation).
  In the stratified case the scale height is given by 
  \f$ H = \sqrt{2}c_s/\Omega \f$. 
  The parameters controlling the model are the elements of the 
  array \c g_inputParam:
 
  - <tt>g_inputParam[BETA]</tt>: sets the plasma beta (midplane),
     \f$ B_z = \sqrt{2c_s^2/\beta} \f$;
  - <tt>g_inputParam[CSOUND]</tt>: sets the speed of sound \f$ c_s \f$.
 
  The additional constants \c NET_FLUX and \c STRATIFICATION are 
  defined inside \c definitions.h and are used to:
 
  - \c NET_FLUX (\c YES/\c NO): set the magnetic field configuration 
     corresponding to a net magnetic flux (\c YES) or zero-net flux
     (\c NO);
  - \c STRATIFICATION (\c YES/\c NO): enable or disable stratified shearing 
    box models. With stratification, density and vertical gravity are changed.
    
  The \c Shearing_Box/ directory contains several configurations 
  for different purposes:

  <CENTER>
  Conf.|DIM|NET_FLUX|STRAT.|T. STEPPING|INTERPOLATION| Ref
  -----|---|--------|------|---------- | -- ---------|--------
   #01 | 2 | YES    | NO   | HANCOCK   | LINEAR      |   -
   #02 | 2 | YES    | NO   | ChTr      | PARABOLIC   |    -
   #03 | 3 | YES    | NO   | RK2       | LINEAR      | [Bod08]
   #04 | 3 | YES    | NO   | RK2       | LINEAR      | [Bod08]
   #05 | 3 | YES    | NO   | ChTr      | PARABOLIC   | [Bod08]
   #06 | 3 | YES    | NO   | ChTr      | PARABOLIC   | [Bod08]
   #07 | 3 | NO     | NO   | RK2       | LINEAR      | [Mig12]
   #08 | 3 | NO     | NO   |RK2        | LINEAR      | [Mig12] (*)
   #09 | 3 | NO     | NO   |ChTr       | PARABOLIC   | [Mig12]
   #10 | 3 | NO     | NO   |ChTr       | PARABOLIC   | [Mig12] (*)
   #11 | 3 | NO     | YES  |HANCOCK    | LINEAR      | [Bod14] 
  </CENTER>

  (*) used with FARGO to be compared to the previous configurations.
      Only difference: we use Omega = cs = 1 instead of Omega = cs = 1.e-3.

  All configurations employs an \c ISOTHERMAL equation of state.

  \author A. Mignone (mignone@ph.unito.it)
  \date   May 23, 2014

  \b References:
     - [Bod08]: "Aspect ratio dependence in magnetorotational instability 
                 shearing box simulations" Bodo et al. A&A (2008) 487, 1
     - [Bod14]: "On the Convergence of Magnetorotational Turbulence in 
                 Stratified Isothermal Shearing Boxes" Bodo et al. (2014) 787 L13
     - [Mig12]: "A conservative orbital advection scheme for simulations of
                 magnetized shear flows with the PLUTO code",
                 Mignone et al., A&A (2012) 545, Sec. 3.3.2
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

double sb_Omega = 1.0; /* Orbital frequency (global variable) */
double sb_q     = 1.5; /* Shear parameter   (global variable) */

/* ********************************************************************* */
void Init (double *v, double x, double y, double z)
/*
 *
 *********************************************************************** */
{
  static int first_call = 1;
  double Bz0, rnd, dvy, cs, H;
  double Lx, Ly, Lz;
  double kx, ky, kz;

  #ifndef SHEARINGBOX
   print1 ("! ShearingBox module has not been included.\n");
   print1 ("! Cannot continue.\n");
   QUIT_PLUTO(1);
  #endif

/* -- compute domain sizes -- */

  Lx = g_domEnd[IDIR] - g_domBeg[IDIR]; kx = 2.0*CONST_PI/Lx;
  Ly = g_domEnd[JDIR] - g_domBeg[JDIR]; ky = 2.0*CONST_PI/Ly;
  Lz = g_domEnd[KDIR] - g_domBeg[KDIR]; kz = 2.0*CONST_PI/Lz;

/* -- get sound speed and pressure scale height -- */

  cs = g_inputParam[CSOUND];    /* sound speed */
  H  = sqrt(2.0)*cs/sb_Omega;   /* pressure scale height */

/* -- seed random numbers differently for each processor -- */

  if (first_call == 1){
    srand(time(NULL) + prank);
    first_call = 0;
  }
  rnd = (double)(rand())/((double)RAND_MAX + 1.0);

/* -- set velocity perturbation [1]: random noise -- */

/*  dvy = 0.01*cs*rnd; */

/* -- set velocity perturbation [2], two harmonics for each direction -- */

  dvy  = sin(kx*x + 0.20) + sin(2.0*kx*x - 0.37);
  dvy *= sin(ky*y + 0.13) + sin(2.0*ky*y + 0.04);
  dvy *= sin(kz*z + 0.56) + sin(2.0*kz*z + 0.62);
  dvy *= 0.01*cs/8.0;

/* -- in 2D we don't use any perturbation -- */

  #if DIMENSIONS == 2
   dvy = 0.0;
  #endif

/* -- set initial condition -- */

  #if STRATIFICATION == YES
   v[RHO] = exp(-z*z/(H*H));
  #else
   v[RHO] = 1.0;
  #endif

  v[VX1] = 0.0;
  v[VX2] = -sb_q*sb_Omega*x + dvy;
  v[VX3] = 0.0;
  #if EOS == IDEAL
   v[PRS] = cs*cs*v[RHO];
  #elif EOS == ISOTHERMAL
   g_isoSoundSpeed = cs;
  #endif

  v[TRC] = 0.0;

/* ----------------------------------------------------------------
    The magnetic field amplitude is set by the parameter
    beta = 2p/B^2 = 2*rho*c^2/b0^2   -->   b0 = c*sqrt(2/beta)
    where it is assumed that rho = 1 (midplane).
   ---------------------------------------------------------------- */

  #if PHYSICS == MHD 
   Bz0 = cs*sqrt(2.0/g_inputParam[BETA]);

   #if  NET_FLUX  == YES  /* -- Net flux, constant vertical field Bz = B0 -- */

    v[BX1] = 0.0;
    v[BX2] = 0.0;
    v[BX3] = Bz0;

    v[AX1] = 0.0;
    v[AX2] = Bz0*x;
    v[AX3] = 0.0;

   #else /* -- Zero net flux, Bz = B0*sin(2*pi*x/Lx) -- */

   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = Bz0*sin(kx*x);

   v[AX1] = 0.0;
   v[AX2] = -Bz0*cos(kx*x)/kx;
   v[AX3] = 0.0;

   #endif

   #if DIMENSIONS == 2  /* 2D Case only for testing */
    v[BX1] = Bz0*sin(ky*y);
    v[BX2] = 0.0;
    v[BX3] = 0.0;

    v[AX1] = 0.0;
    v[AX2] = 0.0;
    v[AX3] = -Bz0*cos(ky*y)/ky;
   #endif
  #endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 * Compute volume-integrated pressure, Maxwell and Reynolds stresses.
 *
 *********************************************************************** */
{
  int    i,j,k;
  static int first_call=1;
  double *dx, *dy, *dz;
  double *x, *y, *z;
  double Lx, Ly, Lz, scrh;
  double pm, Mxy, Rxy, tot_vol, aM, aR, dV;
  FILE *fp;

/* -- compute domain sizes and box volume -- */

  Lx = g_domEnd[IDIR] - g_domBeg[IDIR];
  Ly = g_domEnd[JDIR] - g_domBeg[JDIR];
  Lz = g_domEnd[KDIR] - g_domBeg[KDIR];

  tot_vol = Lx*Ly*Lz;

/* -- pointers to mesh spacing and coordinates -- */

  dx = grid[IDIR].dx; dy = grid[JDIR].dx; dz = grid[KDIR].dx;
   x = grid[IDIR].x;   y = grid[JDIR].x;   z = grid[KDIR].x;

/* ------------------------------------------------------------
    Main analysis loop.
    Compute volume-integrated magnetic pressure and stresses
   ------------------------------------------------------------ */

  pm = Mxy = Rxy = 0.0;
  DOM_LOOP(k,j,i){
    dV  = dx[i]*dy[j]*dz[k];

  /* -- magnetic pressure -- */

    pm += 0.5*(EXPAND(  d->Vc[BX1][k][j][i]*d->Vc[BX1][k][j][i],
                      + d->Vc[BX2][k][j][i]*d->Vc[BX2][k][j][i],
                      + d->Vc[BX3][k][j][i]*d->Vc[BX3][k][j][i]))*dV;

  /* -- Maxwell and Reynolds stresses -- */

    aM = d->Vc[BX1][k][j][i]*d->Vc[BX2][k][j][i];
    aR = d->Vc[RHO][k][j][i]*d->Vc[VX1][k][j][i]*
        (d->Vc[VX2][k][j][i] + sb_q*sb_Omega*x[i]);

    Mxy   += aM*dV;
    Rxy   += aR*dV;
  }

/* -- divide summations by total volume -- */

  pm  /= tot_vol;
  Mxy /= tot_vol;
  Rxy /= tot_vol;

/* -- parallel reduce -- */

  #ifdef PARALLEL
   MPI_Allreduce (&pm, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   pm = scrh;

   MPI_Allreduce (&Mxy, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   Mxy = scrh;

   MPI_Allreduce (&Rxy, &scrh, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   Rxy = scrh;

   MPI_Barrier (MPI_COMM_WORLD);
  #endif

/* -- only proc #0 writes ascii data-file to disk -- */

  if (prank == 0){
    static double tpos = -1.0;
    if (g_stepNumber == 0){   /* -- open for writing if initial step -- */
      fp = fopen("averages.dat","w");
      fprintf (fp,"# %4s  %12s  %12s  %12s  %12s\n",
                   "step","  time  "," <B^2/2> "," <Bx*By> ","<rho*ux*duy>");
      first_call = 0;
    }else{                 
      if (tpos < 0.0){  /* obtain time coordinate of last written line */
        char   sline[512];
        fp = fopen("averages.dat","r");
        while (fgets(sline, 512, fp))  {}
        sscanf(sline, "%lf\n",&tpos);
        fclose(fp);
      }
      fp = fopen("averages.dat","a");
    }
    if (g_time > tpos){ /* -- write -- */
      fprintf (fp, "%4d    %12.6e  %12.6e  %12.6e  %12.6e\n", 
               g_stepNumber, g_time, pm, Mxy, Rxy);
    }
    fclose(fp);
  }
}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/* 
 *
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x, *y, *z, *zi, *dz;
  double ***rho, ***prs;
  double ***Bx, ***By, ***Bz;
  double ***vx, ***vy, ***vz;
  double ***Bxs, ***Bys, ***Bzs;
  double H, cs;

  cs = g_inputParam[CSOUND];    /* sound speed */
  H  = sqrt(2.0)*cs/sb_Omega;   /* pressure scale height */

/* -- set pointer shortcuts -- */

  rho = d->Vc[RHO];
  #if HAVE_PRESSURE
   prs = d->Vc[PRS];
  #endif
  EXPAND(vx  = d->Vc[VX1]; Bx = d->Vc[BX1];  ,
         vy  = d->Vc[VX2]; By = d->Vc[BX2];  ,
         vz  = d->Vc[VX3]; Bz = d->Vc[BX3];)
  #ifdef STAGGERED_MHD
   EXPAND(Bxs  = d->Vs[BX1s];  ,
          Bys  = d->Vs[BX2s];  ,
          Bzs  = d->Vs[BX3s];)
  #endif

  x  = grid[IDIR].x; 
  y  = grid[JDIR].x;
  z  = grid[KDIR].x;
  zi = grid[KDIR].xr;
  dz = grid[KDIR].dx;

  if (side == 0) {    /* -- Density threshold -- */
    DOM_LOOP(k,j,i){
      if (d->Vc[RHO][k][j][i] < 1.e-5) d->Vc[RHO][k][j][i] = 1.e-5;
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    double gz, a;
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        gz = -g_domBeg[KDIR]*sb_Omega*sb_Omega;  /* > 0 */
        a  = 0.5*dz[k]*gz/(g_isoSoundSpeed*g_isoSoundSpeed);
        rho[k][j][i] = rho[KBEG][j][i]*(1.0 - a)/(1.0 + a);
        vx[k][j][i]  =  vx[2*KBEG-k-1][j][i];
        vy[k][j][i]  =  vy[2*KBEG-k-1][j][i];
        vz[k][j][i]  = -vz[2*KBEG-k-1][j][i];

        Bx[k][j][i]  = -Bx[2*KBEG-k-1][j][i];
        By[k][j][i]  = -By[2*KBEG-k-1][j][i];
        Bz[k][j][i]  =  Bz[2*KBEG-k-1][j][i];
      }
    }
    #ifdef STAGGERED_MHD
     if (box->vpos == X1FACE){
       BOX_LOOP(box,k,j,i) Bxs[k][j][i] = -Bxs[2*KBEG-k-1][j][i];
     }else if (box->vpos == X2FACE){
       BOX_LOOP(box,k,j,i) Bys[k][j][i] = -Bys[2*KBEG-k-1][j][i];
     }
    #endif
  }

  if (side == X3_END) {  /* -- X3_END boundary -- */
    double gz, a;
    if (box->vpos == CENTER){
      BOX_LOOP(box,k,j,i){
        gz = -g_domEnd[KDIR]*sb_Omega*sb_Omega;  /* < 0 */
        a  = 0.5*dz[k]*gz/(g_isoSoundSpeed*g_isoSoundSpeed);
        rho[k][j][i] = rho[KEND][j][i]*(1.0 + a)/(1.0 - a);
        vx[k][j][i]  =  vx[2*KEND-k+1][j][i];
        vy[k][j][i]  =  vy[2*KEND-k+1][j][i];
        vz[k][j][i]  = -vz[2*KEND-k+1][j][i];

        Bx[k][j][i]  = -Bx[2*KEND-k+1][j][i];
        By[k][j][i]  = -By[2*KEND-k+1][j][i];
        Bz[k][j][i]  =  Bz[2*KEND-k+1][j][i];
      }
    }
    #ifdef STAGGERED_MHD
     if (box->vpos == X1FACE){
       BOX_LOOP(box,k,j,i) Bxs[k][j][i] = -Bxs[2*KEND-k+1][j][i];
     }else if (box->vpos == X2FACE){
       BOX_LOOP(box,k,j,i) Bys[k][j][i] = -Bys[2*KEND-k+1][j][i];
     }
    #endif
  }
	
}

/* ********************************************************************* */
void BodyForceVector (double *v, double *g, double x, double y, double z)
/*
 *  Include gravitational force in the shearing box module.
 *  Coriolis terms are included elsewhere.
 *  
 *********************************************************************** */
{
  double om2  = sb_Omega*sb_Omega;
  double zb  = g_domBeg[KDIR], ze  = g_domEnd[KDIR];

  #ifdef FARGO
   g[IDIR] = 0.0;
   g[JDIR] = sb_q*sb_Omega*v[VX1];
  #else
   g[IDIR] = om2*2.0*sb_q*x;
   g[JDIR] = 0.0;
  #endif

  #if STRATIFICATION == YES
   g[KDIR] = -om2*z;
  #else
   g[KDIR] = 0.0;
  #endif

/* -----------------------------------------------------------
    The following modification should be turned on only
    when a reflective vertical boundary is adopted, since
    the predictor step of CTU will not preserve the symmetry
    or anti-symmetry.
   ----------------------------------------------------------- */

  #if STRATIFICATION == YES && (defined CTU)
   if (z > ze) g[KDIR] += 2.0*om2*ze;
   if (z < zb) g[KDIR] += 2.0*om2*zb;
  #endif

}
/* ********************************************************************* */
double BodyForcePotential(double x, double y, double z)
/*
 * Return the gravitational potential as function of the coordinates.
 *
 *********************************************************************** */
{
  double om2  = sb_Omega*sb_Omega;
  double zb  = g_domBeg[KDIR], ze  = g_domEnd[KDIR];
  double psi;
 
  #if STRATIFICATION == YES
   psi = om2*(0.5*z*z - sb_q*x*x);
   #ifdef CTU /* see note in BodyForceVector() */
    if (z > ze) psi -= 2.0*om2*ze*(z - ze);
    if (z < zb) psi -= 2.0*om2*zb*(z - zb); 
   #endif
   return psi;
  #else
   return -om2*sb_q*x*x;
  #endif
}

/* ************************************************************** */
double FARGO_SetVelocity(double x1, double x2)
/*
 *
 * PURPOSE
 *
 *   Compute the shear angular velocity to be subtracted from 
 *   the HD or MHD equations.
 * 
 **************************************************************** */
{
  return -sb_q*sb_Omega*x1;
}
