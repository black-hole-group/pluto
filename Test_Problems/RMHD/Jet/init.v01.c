#include "pluto.h"

static void   GetJetValues (double x1, double x2, double x3, double *vj);
static double Profile (double, double, double);

static double pmin, a;

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
  GetJetValues (x1, x2, x3, vjet);

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
/* 
 *
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double *x1, *x2, *x3;
  double ***bxs, ***bys, ***bzs;
  double r, phi, prof, vjet[256], vout[NVAR];

  #ifdef STAGGERED_MHD
   D_EXPAND(bxs = d->Vs[BX1s];  , 
            bys = d->Vs[BX2s];  , 
            bzs = d->Vs[BX3s];)
  #endif

  x1  = grid[IDIR].xgc;  
  x2  = grid[JDIR].xgc;  
  x3  = grid[KDIR].xgc;

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER){    /* -- cell-centered boundary conditions -- */
      BOX_LOOP(box,k,j,i){
 
      /* -- obtain beam values -- */

        GetJetValues (x1[i], x2[j], x3[k], vjet);

      /* -- ambient values -- */

        VAR_LOOP(nv) vout[nv] = d->Vc[nv][k][2*JBEG - j - 1][i];
        vout[VX2] *= -1.0; 
        EXPAND(vout[BX1] *= -1.0;  ,  
                                ;  , 
               vout[BX3] *= -1.0;)
        #ifdef PSI_GLM
         vjet[PSI_GLM]  = d->Vc[PSI_GLM][k][JBEG][i] - d->Vc[BX2][k][JBEG][i];
         vout[PSI_GLM] *= -1.0;
        #endif

        prof = Profile(x1[i],x2[j],x3[k]);
        VAR_LOOP(nv) d->Vc[nv][k][j][i] = vout[nv] - (vout[nv] - vjet[nv])*prof;
      }

    }else if (box->vpos == X1FACE){  /* -- staggered fields -- */
      #ifdef STAGGERED_MHD
       x1 = grid[IDIR].xr;
       BOX_LOOP(box,k,j,i){
          vout[BX1] = -bxs[k][2*JBEG - j - 1][i];
          prof = Profile(x1[i],x2[j],x3[k]);
          bxs[k][j][i] = vout[BX1] - (vout[BX1] - vjet[BX1])*prof;
       }
      #endif
    }
  }
}

/* **************************************************************** */
void GetJetValues (double x1, double x2, double x3, double *vj)
/*
 *
 *
 *
 * 
 ****************************************************************** */
{
  static int  first_call = 1;
  double bm, b2_av, bphi;
  double r, z, phi, x, scrh;

  a = 0.5;
  r = x1;
  z = x2;

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

  scrh    = MIN(x*x, 1.0);
  vj[PRS] = pmin + bm*bm*(1.0 - scrh);

  bphi  = bm*(fabs(x) < 1.0 ? x: 1.0/x);
  bphi *= g_inputParam[LORENTZ]*(r <= 1.0 ? 1.0:0.0);

  #if GEOMETRY == CYLINDRICAL
   EXPAND(vj[iBR]   = 0.0;                            ,
          vj[iBZ]   = sqrt(2.0*g_inputParam[SIGMA_POL]*pmin);  ,
          vj[iBPHI] = bphi;)
   vj[AX1] = 0.0;
   vj[AX2] = 0.0;
   vj[AX3] = 0.5*r*vj[iBZ];
  #endif

  if (first_call){
    print1 ("pmin  = %f\n",pmin);
    first_call = 0;
  }

/* -- set tracer value -- */

  #if NTRACER > 0
   vj[TRC] = 1.0;
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

  r = x1; 
  z = x2;

  scrh = pow(r, 14.0);
  scrh = 1.0/cosh(scrh)*(z < 2.0);
  return((fabs(r) <= 1.0)*(z < 0.0));
}

