#include "pluto.h"

#define VISCOSITY_COEFF 0.1

/* *********************************************************** */
void VISC_FLUX (const State_1D *state, int beg, int end, Grid *grid)
/* 
 *
 * PURPOSE
 *
 *   Add artificial viscosity term to the fluxes
 *
 *
 ************************************************************* */
{
  int   i, nv;
  double  scrh;
  double  **flux, **v, **u;
  double *uL, *uR, *vL, *vR;

  flux = state->flux;

 /* -- OLD version used cell-centered viscosity -- */
/*
  v    = state->v;
  u    = state->u;
  for (i = beg; i <= end; i++) {
    scrh = VISCOSITY_COEFF*MAX(0.0, v[i][VXn] - v[i+1][VXn]);
    for (nv = 0; nv < NFLX; nv++) {
      flux[i][nv] += scrh*(u[i][nv] - u[i + 1][nv]);
    }
  }
*/
  for (i = beg; i <= end; i++) {
    uL = state->uL[i]; uR = state->uR[i];
    vL = state->vL[i]; vR = state->vR[i];
    scrh = VISCOSITY_COEFF*(MAX(0.0, vL[VXn] - vR[VXn]));
    for (nv = 0; nv < NFLX; nv++) {
      flux[i][nv] += scrh*(uL[nv] - uR[nv]);
    }
  }

}

#undef VISCOSITY_COEFF
