#include "pluto.h"

static void JETVAL (real x1, real x2, real x3, double *vj);
static double Profile (double, double, double);
static void Perturbations (double phi, double *vj);
static double BTOR (double r);
static double YVECPOT (double x1, double x2, double x3);

static double pmin, a, bm;

/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  int    nv;
  double r, z, p0, vjet[256], vamb[256];

  g_gamma = 5./3.;
  JETVAL (x1, x2, x3, vjet);

/* -- set ambient values -- */

  vamb[RHO] = g_inputParam[RHOA];
  EXPAND(vamb[VX1] = 0.0;   ,
         vamb[VX2] = 0.0;   ,
         vamb[VX3] = 0.0;)
  vamb[PRS] = pmin;
  #if PHYSICS == RMHD
   EXPAND(vamb[BX1] = 0.0;       ,
          vamb[BX2] = vjet[BX2];  ,
          vamb[BX3] = 0.0;)
   vamb[AX1] = 0.0;
   vamb[AX2] = 0.0;
   vamb[AX3] = vjet[AX3];
  #endif
  #if NTRACER > 0
   vamb[TR] = 0.0;
  #endif
  #ifdef PSI_GLM 
   vamb[PSI_GLM] = 0.0;
  #endif

  for (nv = 0; nv < 2*NVAR; nv++){  /* 2*NVAR: include vec pot as well */
    us[nv] = vamb[nv] - (vamb[nv] - vjet[nv])*Profile(x1,x2,x3);
  }
  g_smallPressure = 1.e-5;
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in/out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies on which side boundary conditions need 
 *                    to be assigned. side can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3;
  double *x1r, *x2r, *x3r;
  double *dx1, *dx2, *dx3;
  double ***bxs, ***bys, ***bzs;
  double r, phi, prof, vjet[256], vout[NVAR];
  double vp[256], vm[256];

  x1  = grid[IDIR].xgc;  x2  = grid[JDIR].xgc;  x3  = grid[KDIR].xgc;
  dx1 = grid[IDIR].dx; dx2 = grid[JDIR].dx; dx3 = grid[KDIR].dx;
  x1r = grid[IDIR].xr; x2r = grid[JDIR].xr; x3r = grid[KDIR].xr;

  #ifdef STAGGERED_MHD
   D_EXPAND(bxs = d->Vs[BX1s];  , 
            bys = d->Vs[BX2s];  , 
            bzs = d->Vs[BX3s];)
  #endif

  if (side == 0) {    /* -- check solution inside domain -- */
    DOM_LOOP(k,j,i){};
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER){    /* -- cell-centered boundary conditions -- */
      BOX_LOOP(box,k,j,i){
 
   /* -- obtain beam values -- */

        JETVAL (x1[i], x2[j], x3[k], vjet);
/*
      #if (PHYSICS == MHD || PHYSICS == RMHD) && DIMENSIONS == 3
       JETVAL (x1[i], x2[j], x3[k] + 0.5*dx3[k], vp);
       JETVAL (x1[i], x2[j], x3[k] - 0.5*dx3[k], vm);
       vjet[BX1] = -(vp[AX2] - vm[AX2])/dx3[k];

       JETVAL (x1[i] + 0.5*dx1[i], x2[j], x3[k], vp);
       JETVAL (x1[i] - 0.5*dx1[i], x2[j], x3[k], vm);
       vjet[BX3] = (vp[AX2] - vm[AX2])/dx1[i];
     #endif
*/
   /* -- add perturbations (3D only) ---- */

        #if GEOMETRY == CARTESIAN && DIMENSIONS == 3
         phi = atan2(x3[k],x1[i]);
         Perturbations(phi, vjet);
        #endif

   /* -- ambient values -- */

        for (nv = NVAR; nv--; ) vout[nv] = d->Vc[nv][k][2*JBEG - j - 1][i];
        vout[VX2] *= -1.0; 
        #if PHYSICS == RMHD
         EXPAND(vout[BX1] *= -1.0;,  , vout[BX3] *= -1.0;)
         #ifdef PSI_GLM
          vjet[PSI_GLM]  = d->Vc[PSI_GLM][k][JBEG][i] - d->Vc[BX2][k][JBEG][i];
          vout[PSI_GLM] *= -1.0;
         #endif
        #endif

        prof = Profile(x1[i],x2[j],x3[k]);
        for (nv = NVAR; nv--; ){
          d->Vc[nv][k][j][i] = vout[nv] - (vout[nv] - vjet[nv])*prof;
        }
      }

    }else if (box->vpos == X1FACE){  /* -- staggered fields -- */

      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i){
         #if DIMENSIONS == 2
          vout[BX1] = -bxs[k][2*JBEG - j - 1][i];
          prof = Profile(x1r[i],x2[j],x3[k]);
          bxs[k][j][i] = vout[BX1] - (vout[BX1] - vjet[BX1])*prof;
         #elif DIMENSIONS == 3
          vout[BX1] = -bxs[k][2*JBEG - j - 1][i];                                                 
       
          vjet[BX1] = -(  YVECPOT(x1r[i], x2[j], x3[k] + 0.5*dx3[k])
                       - YVECPOT(x1r[i], x2[j], x3[k] - 0.5*dx3[k]))/dx3[k];

          prof = Profile(x1r[i],x2[j],x3[k]);
          bxs[k][j][i] = vout[BX1] - (vout[BX1] - vjet[BX1])*prof;  
         #endif
       }
      #endif

    }else if (box->vpos == X3FACE){

      #ifdef STAGGERED_MHD
       BOX_LOOP(box,k,j,i){
         vout[BX3] = -bzs[k][2*JBEG - j - 1][i];
         vjet[BX3] = (  YVECPOT(x1[i] + 0.5*dx1[i], x2[j], x3r[k])
                     - YVECPOT(x1[i] - 0.5*dx1[i], x2[j], x3r[k]))/dx1[i];
         prof = Profile(x1[i],x2[j],x3r[k]);
         bzs[k][j][i] = vout[BX3] - (vout[BX3] - vjet[BX3])*prof;
       }
      #endif
    }
  }
}

/* **************************************************************** */
void JETVAL (real x1, real x2, real x3, double *vj)
/*
 *
 *
 *
 * 
 ****************************************************************** */
{
  static int  first_call = 1;
  double b2_av, bphi;
  double r, z, phi, x, scrh;

  a = 0.5;
  #if GEOMETRY == CYLINDRICAL
   r   = x1;
   phi = 0.0;
  #elif GEOMETRY == CARTESIAN
   r   = D_EXPAND(x1*x1,   , + x3*x3);
   r   = sqrt(r);
   phi = atan2(x3,x1);
   #if DIMENSIONS == 2
    r = x1; phi = 0.0;
   #endif
  #endif
  z   = x2;

  if (fabs(r) < 1.e-9) r = 1.e-9;
  if (fabs(z) < 1.e-9) z = 1.e-9;

  x = r/a;
  vj[RHO] = g_inputParam[RHOJ];

  EXPAND(vj[VX1] = 0.0;                                        ,
         vj[VX2] = sqrt(1.0 - 1.0/g_inputParam[LORENTZ]/g_inputParam[LORENTZ]);  ,
         vj[VX3] = 0.0;)

  scrh = g_inputParam[MACH]/vj[VX2];
  pmin = vj[RHO]/(g_gamma*scrh*scrh - g_gamma/(g_gamma - 1.0));

  bm  = 4.0*pmin*g_inputParam[SIGMA_TOR];
  bm /= a*a*(1.0 - 4.0*log(a));
  bm  = sqrt(bm);

  scrh   = MIN(x*x, 1.0);
  vj[PRS] = pmin + bm*bm*(1.0 - scrh);

  bphi = BTOR(r);

  #if PHYSICS == RMHD
   #if GEOMETRY == CYLINDRICAL
    EXPAND(vj[iBR]   = 0.0;                            ,
           vj[iBZ]   = sqrt(2.0*g_inputParam[SIGMA_POL]*pmin);  ,
           vj[iBPHI] = bphi;)
    vj[AX1] = 0.0;
    vj[AX2] = 0.0;
    vj[AX3] = 0.5*r*vj[iBZ];
   #elif GEOMETRY == CARTESIAN
    EXPAND(vj[BX1] = -bphi*sin(phi);                  ,
           vj[BX2] =  sqrt(2.0*g_inputParam[SIGMA_POL]*pmin);  ,
           vj[BX3] =  bphi*cos(phi);)
    vj[AX1] = 0.0; 
    vj[AX2] = YVECPOT(x1, x2, x3);
    vj[AX3] = -vj[BX2]*x1;
   #endif
   
  #endif

  if (first_call){
    print1 ("pmin  = %f\n",pmin);
    first_call = 0;
  }

/* -- set tracer value -- */

  #if NTRACER > 0
   vj[TR] = 1.0;
  #endif
  #ifdef PSI_GLM 
   vj[PSI_GLM] = 0.0;
  #endif
}

/* ******************************************************* */
double Profile (double x1, double x2, double x3)
/*
 *
 *
 ********************************************************* */
{
  double r, z, scrh;

  #if GEOMETRY == CARTESIAN
   r = sqrt(x1*x1 + x3*x3);
  #elif GEOMETRY == CYLINDRICAL
   r = x1; 
  #endif
  z = x2;

  scrh = pow(r, 14.0);
  scrh = 1.0/cosh(scrh)*(z < 2.0);
//  return(scrh);
  return((fabs(r) <= 1.0)*(z < 0.0));
}

/* ****************************************************** */
void Perturbations (double phi, double *vj)
/*
 *
 *
 ******************************************************** */
{
  int m, l;
  double scrh, amp, om1, tom1;
  double arg = 0.0;
  const double oml[8] = {0.5, 1.0, 2.0, 3.0, 0.03, 0.06, 0.12, 0.25};
  const double bl[8]  = {0.2, 1.6, 2.3, 4.9, 5.9 , 5.4 , 0.7 , 3.1};

  scrh = 1.0 + g_inputParam[EPS];
  amp  = sqrt(scrh*scrh - 1.0)/(g_inputParam[LORENTZ]*scrh);
  om1  = sqrt(1.0 - 1.0/g_inputParam[LORENTZ]/g_inputParam[LORENTZ])/g_inputParam[MACH];
  tom1 = g_time*om1;
  for(m = 0; m <= 2; m++){
  for(l = 0; l < 8; l++){
    arg += cos(m*phi + oml[l]*tom1 + bl[l]);
  }}
  arg *= amp/24.0;
  vj[VX1] += arg*cos(phi);
  vj[VX3] += arg*sin(phi);
}


/* ********************************************************* */
double BTOR (double r)
/*
 *
 *
 *********************************************************** */
{
  double scrh, bphi,x;

  x   = r/a;
  scrh  = MIN(x*x, 1.0);
  bphi  = bm*(fabs(x) < 1.0 ? x: 1.0/x);
  #if PHYSICS == RMHD
   bphi *= g_inputParam[LORENTZ];
  #endif

  bphi *= (r <= 1.0);
//  bphi *= 1.0/cosh(pow(r,10.0));
  return (bphi);
}

/* *********************************************************** */
double YVECPOT (double x1, double x2, double x3)
/*
 *
 *  Return the vector potential associated with 
 *  the toroidal component of mag. field
 *
 ************************************************************* */
{
  int    i, nr = 2048;
  double r, rmax, dr, ay;
  static double *psi, *rpos;
  

  rmax = 1.5;
  dr   = rmax/(double)nr;
  if (psi == NULL){
    psi  = ARRAY_1D(nr, double);
    rpos = ARRAY_1D(nr, double);
    psi[0]  = dr*BTOR(dr*0.5);
    rpos[0] = 0.0;
    for (i = 1; i < nr; i++){
      rpos[i] = i*dr; 
      psi[i]  = psi[i - 1] + dr*BTOR((i+.5)*dr);
    }

  /* -- make ay continuous at r = rmax -- */

    for (i = 0; i < nr; i++) psi[i] -= psi[nr-1];

  }
  
  r  = sqrt(x1*x1 + x3*x3);
  i  = floor(r/dr);
  if (i > (nr-2)) return(0.0);
  ay = (psi[i]*(rpos[i+1] - r) + psi[i+1]*(r - rpos[i]))/dr;

  return (ay);
}
