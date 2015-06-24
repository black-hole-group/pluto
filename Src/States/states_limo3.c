#include "pluto.h"
static double LIM_FUNC (double dvp, double dvm, double);

/* ************************************************************* */
void States (const State_1D *state, int beg, int end, Grid *grid)
/* 
 *
 * PURPOSE
 *
 *   Provide a three-point stencil, third-order 
 *   reconstruction algorithm based on the limiter 
 *   function introduced  by:
 *
 *   "Compact third-order limiter functions for finite volume
 *    methods", Cada & Torrilhon, JCP (2009) 228, 4118.
 *
 * LAST MODIFIED
 *
 *   Nov 3rd 2009
 *   by A. Mignone (mignone@ph.unito.it)
 *
 **************************************************************** */
{
  int    k, nv, i;
  double dmm;
  double **v, *vp, *vm, *dvp, *dvm, *dx;
  static double **dv;
  double **L, **R, *lambda;
  double dwp[NVAR], dwp_lim[NVAR];
  double dwm[NVAR], dwm_lim[NVAR];
  double dvpR, dvmR;

  if (dv == NULL) {
    dv = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  v  = state->v;
  dx = grid[g_dir].dx;

/* ----------------------------------------------------
    compute slopes and left and right interface values 
   ---------------------------------------------------- */

  #if CHAR_LIMITING == NO  /* ----------------------------------------
                                 Limiter on primitive variables
                              ----------------------------------------  */
   for (i = beg-1; i <= end; i++){
   for (nv = NVAR; nv--; ){
     dv[i][nv] = v[i+1][nv] - v[i][nv];
   }}

   for (i = beg; i <= end; i++){
     dvp = dv[i];   vp = state->vp[i];
     dvm = dv[i-1]; vm = state->vm[i];

     #if SHOCK_FLATTENING == MULTID    
      if (CheckZone (i, FLAG_MINMOD)) {
        for (nv = NVAR; nv--;    ) {
          dmm = MINMOD(dvp[nv], dvm[nv]);
          vp[nv] = v[i][nv] + 0.5*dmm;
          vm[nv] = v[i][nv] - 0.5*dmm;
        }
        continue;
      }
     #endif

     for (nv = NVAR; nv--; ){
       vp[nv] = v[i][nv] + 0.5*dvp[nv]*LIM_FUNC(dvp[nv], dvm[nv], dx[i]);
       vm[nv] = v[i][nv] - 0.5*dvm[nv]*LIM_FUNC(dvm[nv], dvp[nv], dx[i]);
     }
   }       
  
  #else       /* --------------------------------------------
                    Limiter on characteristic variables
                 --------------------------------------------  */

   SoundSpeed2 (state->v, state->a2, state->h, beg, end, CELL_CENTER, grid);
   i = beg-1;
   for (nv = NVAR; nv--;   ) dv[i][nv] = v[i+1][nv] - v[i][nv];

   for (i = beg; i <= end; i++){

     for (nv = NVAR; nv--; ) dv[i][nv] = v[i+1][nv] - v[i][nv];
     L      = state->Lp[i];
     R      = state->Rp[i];
     lambda = state->lambda[i];
     dvp = dv[i];   vp = state->vp[i];
     dvm = dv[i-1]; vm = state->vm[i];
    
  /* -------------------------------
      project undivided differences 
      onto characteristic space
     ------------------------------- */
     
     PrimEigenvectors (state->v[i], state->a2[i], state->h[i], lambda, L, R);
     PrimToChar(L, dvp, dwp);
     PrimToChar(L, dvm, dwm);

     #if SHOCK_FLATTENING == MULTID    
      if (CheckZone (i, FLAG_MINMOD)) {
        for (nv = NVAR; nv--;    ) {
          dmm = MINMOD(dvp[nv], dvm[nv]);
          vp[nv] = v[i][nv] + 0.5*dmm;
          vm[nv] = v[i][nv] - 0.5*dmm;
        }
        continue;
      }
     #endif

  /* -----------------------------
      limit undivided differences
     ----------------------------- */

     for (k = NFLX; k--; ){
       dwp_lim[k] = dwp[k]*LIM_FUNC(dwp[k], dwm[k], dx[i]);
       dwm_lim[k] = dwm[k]*LIM_FUNC(dwm[k], dwp[k], dx[i]);
     }
     for (nv = NFLX; nv--;   ){
       dvpR = dvmR = 0.0;
       #ifdef STAGGERED_MHD
        if (nv == BXn) continue;
       #endif
       for (k = NFLX; k--; ){
         dvpR += dwp_lim[k]*R[nv][k];
         dvmR += dwm_lim[k]*R[nv][k];
       }
       vp[nv] = v[i][nv] + 0.5*dvpR;
       vm[nv] = v[i][nv] - 0.5*dvmR;
     }

  /* -------------------------------------- 
      Compute limited slopes for tracers
      exploiting the simple characteristic 
      structure, L=R=diag(1).
     -------------------------------------- */            

    #if NFLX != NVAR
     for (nv = NFLX; nv < NVAR; nv++ ){
       dvpR = dvp[nv]*LIM_FUNC(dvp[nv], dvm[nv], dx[i]);
       dvmR = dvm[nv]*LIM_FUNC(dvm[nv], dvp[nv], dx[i]);
       vp[nv] = v[i][nv] + 0.5*dvpR;
       vm[nv] = v[i][nv] - 0.5*dvmR;
     }
    #endif
   }
   
  #endif /* CHAR_LIMITING == YES */
   
/* --------------------------------------------------
    Since the third-order limiter is not TVD, we 
    need to ensure positivity of density and pressure 
   -------------------------------------------------- */

  for (i = beg; i <= end; i++){
    dvp = dv[i];   vp = state->vp[i];
    dvm = dv[i-1]; vm = state->vm[i];
    if (vp[RHO] < 0.0 || vm[RHO] < 0.0){
      dmm = MINMOD(dvp[RHO], dvm[RHO]);
      vp[RHO] = v[i][RHO] + 0.5*dmm;
      vm[RHO] = v[i][RHO] - 0.5*dmm;
    }

    #if HAVE_ENERGY
     if (vp[PRS] < 0.0 || vm[PRS] < 0.0){
       dmm = MINMOD(dvp[PRS], dvm[PRS]);
       vp[PRS] = v[i][PRS] + 0.5*dmm;
       vm[PRS] = v[i][PRS] - 0.5*dmm;
     }
    #endif
    #if ENTROPY_SWITCH == YES
     if (vp[ENTR] < 0.0 || vm[ENTR] < 0.0){
       dmm = MINMOD(dvp[ENTR], dvm[ENTR]);
       vp[ENTR] = v[i][ENTR] + 0.5*dmm;
       vm[ENTR] = v[i][ENTR] - 0.5*dmm;
     }
    #endif      

  /* -- relativistic limiter --*/

    #if PHYSICS == RHD || PHYSICS == RMHD
     VelocityLimiter (v[i], vp, vm);
    #endif
  }

/*  -------------------------------------------
        Shock flattening 
    -------------------------------------------  */

  #if SHOCK_FLATTENING == ONED 
   Flatten (state, beg, end, grid);
  #endif

/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg - 1; i <= end; i++) {
     state->vR[i][BXn] = state->vL[i][BXn] = state->bn[i];
   }
  #endif

/* --------------------------------------------------------
      evolve center values by dt/2
   -------------------------------------------------------- */

  #if TIME_STEPPING == CHARACTERISTIC_TRACING
   CharTracingStep(state, beg, end, grid);
  #elif TIME_STEPPING == HANCOCK
   HancockStep(state, beg, end, grid);
  #endif

/* -------------------------------------------
    compute states in conservative variables
   ------------------------------------------- */

  PrimToCons (state->vp, state->up, beg, end);
  PrimToCons (state->vm, state->um, beg, end);

}

/* ***************************************************** */
double LIM_FUNC (double dvp, double dvm, double dx)
/*
 *
 *  Implement the 3rd-order limiter function, 
 *  Eq. 3.22
 *
 *
 ******************************************************* */
{
  double r = 0.1;
  double a,b,c,q, th, lim;
  double eta, psi, eps = 1.e-12;

  th  = dvm/(dvp + 1.e-16);

  q = (2.0 + th)/3.0;

  a = MIN(1.5,2.0*th);
  a = MIN(q,a);
  b = MAX(-0.5*th,a);
  c = MIN(q,b);
  psi = MAX(0.0,c);

  eta = r*dx;
  eta = (dvm*dvm + dvp*dvp)/(eta*eta);
  if ( eta <= 1.0 - eps) {
    lim = q;
  }else if (eta >= 1.0 + eps){
    lim = psi;
  }else{
    psi =   (1.0 - (eta - 1.0)/eps)*q
          + (1.0 + (eta - 1.0)/eps)*psi;
    lim = 0.5*psi;
  }
  return (lim);
}



