#include "pluto.h"

/* ***************************************************************** */
void Flux (double **ucons, double **uprim, double *h, 
           double **fx, real *pr, int beg, int end)
/*
 *
 *
 *
 ******************************************************************* */
{
  int    nv, i;
  double vB, b2, wtg2, Bmag2, pt; 
  double *u, *v, b[4], g, g2, g_2;

  for (i = beg; i <= end; i++) {

    u = ucons[i];
    v = uprim[i];

    g     = u[RHO]/v[RHO];
    g2    = g*g;
    g_2   = 1.0/g2;

    vB    = EXPAND(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
    Bmag2 = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);    

    EXPAND(b[IDIR] = g*(v[BX1]*g_2 + vB*v[VX1]); , 
           b[JDIR] = g*(v[BX2]*g_2 + vB*v[VX2]); ,
           b[KDIR] = g*(v[BX3]*g_2 + vB*v[VX3]);)

    b2 = Bmag2*g_2 + vB*vB;
   
    pt = v[PRS] + 0.5*b2;
    
    wtg2 = (v[RHO]*h[i] + b2)*g2;
    
    fx[i][RHO]  = u[RHO]*v[VXn];
    EXPAND(fx[i][MX1] = wtg2*v[VX1]*v[VXn] - b[IDIR]*b[g_dir];  ,
           fx[i][MX2] = wtg2*v[VX2]*v[VXn] - b[JDIR]*b[g_dir];  ,
           fx[i][MX3] = wtg2*v[VX3]*v[VXn] - b[KDIR]*b[g_dir];)
    EXPAND(fx[i][BXn] = 0.0;             ,
           fx[i][BXt] = v[VXn]*v[BXt] - v[BXn]*v[VXt];   ,
           fx[i][BXb] = v[VXn]*v[BXb] - v[BXn]*v[VXb]; )
    #if RESISTIVE_RMHD == YES
     fx[i][EC] = ??;  
     EXPAND(fx[i][E1] = ??;    ,
            fx[i][E2] = ??;    ,
            fx[i][E3] = ??; )
    #endif

    fx[i][ENG] = u[MXn];
    #if SUBTRACT_DENSITY
     fx[i][ENG] -= fx[i][RHO];
    #endif

    pr[i] = pt;

    #ifdef GLM_MHD
     fx[i][BXn]      = v[PSI_GLM];
     fx[i][PSI_GLM] = glm_ch*glm_ch*v[BXn];
    #endif

  }
}
