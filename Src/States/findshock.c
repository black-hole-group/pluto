#include"pluto.h"

#ifndef EPS_P
 #define EPS_P 5.0
#endif

#if SHOCK_FLATTENING == MULTID
/* *************************************************************** */
void FindShock (const Data *d, Grid *grid)
/*
 *
 *
 * PURPOSE
 *
 *  Check if strong shocks are present in the computational domain.  
 *  If so, revert to minmod limiter in the interpolation 
 *  routines and to HLL in the Riemann solver.
 *
 *  Example: if "X" marks the shock, then flag zones as follows:
 * 
 *  |---|---|---|---|---|---| 
 *            X
 *        M   M   M
 *          H   H
 *
 *  M = FLAG_MINMOD
 *  H = FLAG_HLL
 *
 *  In more than 1D, the flag strategy is extended to neighbor
 *  zones.
 *  
 ***************************************************************** */
{
  int  i, j, k, nv;
  int  ibeg, jbeg, kbeg;
  int  iend, jend, kend;
  int  ip, jp, kp;
  real dpx, dvx, pxmin;
  real dpy, dvy, pymin;
  real dpz, dvz, pzmin;
  real divv, pmin, gradp;
  real ***rho, ***vx, ***vy, ***vz, ***pr;
  unsigned char ***flag;
  double *dVx, *dVy, *dVz;
  double *Ar, *Ath, *r, *th, s;
  
/* -----------------------------------------------
     Define pointers to variables and geometrical 
     factors.
   ----------------------------------------------- */

  rho = pr = d->Vc[RHO];
  EXPAND(vx = d->Vc[VX1];  ,
         vy = d->Vc[VX2];  ,
         vz = d->Vc[VX3];)
  #if HAVE_ENERGY
   pr = d->Vc[PRS];
  #endif

  flag = d->flag;

  dVx = grid[IDIR].dV;
  dVy = grid[JDIR].dV; 
  dVz = grid[KDIR].dV; 
  
  Ar  = grid[IDIR].A;
  r   = grid[IDIR].x;
  Ath = grid[JDIR].A;
  th  = grid[JDIR].x;

/* ----------------------------------------------
       define maximum allowable stencil 
   ---------------------------------------------- */

  jbeg = jend = kbeg = kend = 0;

  D_EXPAND(ibeg = 1; iend = NX1_TOT - 2;  ,   
           jbeg = 1; jend = NX2_TOT - 2;  ,  
           kbeg = 1; kend = NX3_TOT - 2;)


  for (k = kbeg; k <= kend; k++){ 
  for (j = jbeg; j <= jend; j++){ 
  for (i = ibeg; i <= iend; i++){

    #if GEOMETRY == CARTESIAN

     D_EXPAND(dvx = vx[k][j][i + 1] - vx[k][j][i - 1];   ,
              dvy = vy[k][j + 1][i] - vy[k][j - 1][i];   ,
              dvz = vz[k + 1][j][i] - vz[k - 1][j][i];)

    #elif GEOMETRY == CYLINDRICAL

     D_EXPAND(dvx =   Ar[i]  *(vx[k][j][i + 1] + vx[k][j][i])
                    - Ar[i-1]*(vx[k][j][i - 1] + vx[k][j][i]);   ,
              dvy = vy[k][j + 1][i] - vy[k][j - 1][i];               , 
              dvz = vz[k + 1][j][i] - vz[k - 1][j][i];)

    #elif GEOMETRY == POLAR

     D_EXPAND(dvx =  Ar[i]  *(vx[k][j][i + 1] + vx[k][j][i])
                   - Ar[i-1]*(vx[k][j][i - 1] + vx[k][j][i]);  ,
              dvy = (vy[k][j + 1][i] - vy[k][j - 1][i])/r[i];      ,
              dvz =  vz[k + 1][j][i] - vz[k - 1][j][i];)

    #elif GEOMETRY == SPHERICAL

     D_EXPAND(dvx =   Ar[i]  *(vx[k][j][i + 1] + vx[k][j][i])
                    - Ar[i-1]*(vx[k][j][i - 1] + vx[k][j][i]);               ,
              dvy = (  Ath[j]  *(vy[k][j + 1][i] + vy[k][j][i])
                     - Ath[j-1]*(vy[k][j - 1][i] + vy[k][j][i]))/fabs(r[i]); ,
              s   = th[j];
              dvz = (vz[k + 1][j][i] - vz[k - 1][j][i])/(r[i]*sin(s));)

    #endif

    divv = D_EXPAND(dvx, + dvy*dVx[i]/dVy[j], + dvz*dVx[i]/dVz[k]);
 
    if (divv < 0.0){
      D_EXPAND(pxmin = MIN(pr[k][j][i + 1], pr[k][j][i - 1]);  ,  
               pymin = MIN(pr[k][j + 1][i], pr[k][j - 1][i]);  , 
               pzmin = MIN(pr[k + 1][j][i], pr[k - 1][j][i]);)   

      D_EXPAND(dpx = fabs(pr[k][j][i + 1] - pr[k][j][i - 1])/pxmin;  ,  
               dpy = fabs(pr[k][j + 1][i] - pr[k][j - 1][i])/pymin;  , 
               dpz = fabs(pr[k + 1][j][i] - pr[k - 1][j][i])/pzmin;)   
                
      D_EXPAND(gradp  = dpx;   ,
               gradp += dpy;   ,
               gradp += dpz;)

      if (gradp > EPS_P) {

        flag[k][j][i]   |= FLAG_MINMOD;

        D_EXPAND(
          flag[k][j][i+1] |= FLAG_MINMOD;
          flag[k][j][i-1] |= FLAG_MINMOD;  ,
          flag[k][j-1][i] |= FLAG_MINMOD;  
          flag[k][j+1][i] |= FLAG_MINMOD;  ,
          flag[k-1][j][i] |= FLAG_MINMOD;
          flag[k+1][j][i] |= FLAG_MINMOD;)

          flag[k][j][i]   |= FLAG_HLL;
      }
    }
      
  }}}

  #ifdef PARALLEL
   AL_Exchange (flag[0][0], SZ_char);
  #endif
/*
{
  double t0, t1, i0, i1;
  FILE *fp;

  t0 = g_time;
  t1 = t0 + g_dt;
  i0 = (int)(t0/0.08);
  i1 = (int)(t1/0.08);

  if ((i0 != i1 && g_intStage == 1) || g_stepNumber == 0) {
    if (g_stepNumber == 0) fp = fopen("flat.out","wb");
    else {fp = fopen("flat.out","r+b");
          fseek(fp, 0, SEEK_END);
         }
    DOM_LOOP(k,j,i){
      dpx = (double) flag[k][j][i];
      fwrite (&dpx, sizeof(double), 1, fp);
    }
    fclose(fp);
  } 
}
*/

}
#endif
#undef EPS_P

