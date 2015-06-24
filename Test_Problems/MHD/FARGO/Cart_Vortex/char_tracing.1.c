#include "pluto.h"

/* ------------------------------------------
    Set PARABOLIC_LIM to 0,1,2 to apply the
    parabolic limiter of PPM to either
    characteristic (0), primitive (1) or 
    both variables (2). Default is 0.
   ------------------------------------------ */

#define PARABOLIC_LIM  1

/* ----------------------------------------------
    Set REF_STATE to 1,2,3 to use cell centered
    value (1), interpolated states (2) or 
    fastest wave (3, the usual PPM prescription)
    Default is 1.
   ---------------------------------------------- */

#define REF_STATE     3

/* *************************************************************************** */
void STATES (const State_1D *state, int beg, int end, real dt, Grid *grid)
/*
 *
 * PURPOSE:
 * 
 *   Compute 1D left and right interface states using the characteristic
 *   decomposition of the quasi-linear form of the equations.
 *   This is done by first extrapolating the cell center value to the 
 *   interface using piecewise limited reconstruction (LINEAR, WENO or PPM) 
 *   on the characteristic variables.
 *   Left and right states are then evolved for the half time step 
 *   using characteristic tracing.
 *   It should be faster than the standard PLUTO char_tracing since everything
 *   (interpolation + time extrapolation) is done is in the same function.
 *
 *
 * LAST MODIFIED
 *
 *   May 21, 2011 by A. Mignone (mignone@ph.unito.it)
 *
 *************************************************************************** */
{
  int    i, j, k, nv, S=1;
  double dtdx, dx, dx2;
  double dwh, d2w;
  double dwp[NVAR], dwm[NVAR];
  double dwm2[NVAR], dwm1[NVAR], dwp1[NVAR], dwp2[NVAR];
  double Smm, Smp, Spp, dvp[NVAR], dvm[NVAR];
  double nu[NVAR], nup=1.0, num=-1.0;
  double *vp, *vm, *v;
  double **L, **R, *lambda;
  double tau, a0, a1, w0, w1;
  double kstp[NVAR];  
  const double one_sixth = 1.0/6.0;
  static double *a2, *h, **src, **dv;

/* --------------------------------------------
     allocate some memory here
   -------------------------------------------- */

  if (dv == NULL){
    dv  = Array_2D(NMAX_POINT, NVAR, double);
    a2  = Array_1D(NMAX_POINT, double);
    h   = Array_1D(NMAX_POINT, double);
    src = Array_2D(NMAX_POINT, NVAR, double);
  } 

/* ---------------------------------------------
    define some useful quantities, compute
    source term and undivided differences
   --------------------------------------------- */

  #if INTERPOLATION == PARABOLIC
   S = 2;
  #endif
  dx   = grid[DIR].dx[beg];
  dx2  = dx*dx;
  dtdx = delta_t/dx;
  SOUND_SPEED2 (state->v, a2, h, beg, end);
  PRIM_SOURCE  (state, beg, end, a2, h, src, grid);

  for (i = beg-S; i <= end+S-1; i++){    
    for (nv = NVAR; nv--;   ) 
      dv[i][nv] = state->v[i+1][nv] - state->v[i][nv];
  }

/* ------------------------------------------
    set the amount of steepening for each 
    characteristic family (only for LINEAR).
    Default is 2, but nonlinear fields may 
    be safely set to 1 for strongly 
    nonlinear problems.
   ------------------------------------------ */

  #if INTERPOLATION == LINEAR
   for (k = NVAR; k--;  ) kstp[k] = 2.0;
/*
   #if PHYSICS == MHD
    kstp[KFASTP] = kstp[KFASTM] = 1.0;
    #if COMPONENTS > 1
     kstp[KSLOWP] = kstp[KSLOWM] = 1.0; 
    #endif
*/
  #endif
/* --------------------------------------------------------------
                    main spatial loop
   -------------------------------------------------------------- */

  for (i = beg; i <= end; i++){    

    v      = state->v[i]; 
    vp     = state->vp[i];
    vm     = state->vm[i];
    L      = state->Lp[i];
    R      = state->Rp[i];
    lambda = state->lambda[i];

    PRIM_EIGENVECTORS(v, a2[i], h[i], lambda, L, R);
    #if NVAR != NFLX
     for (k = NFLX; k < NVAR; k++) lambda[k] = v[V1]; 
    #endif
    for (k = 0; k < NVAR; k++) nu[k] = dtdx*lambda[k];

    #if SHOCK_FLATTENING == MULTID    
     if (CHECK_ZONE (i, FLAG_MINMOD)) {
       for (nv = 0; nv < NVAR; nv++){  
         dwh = MINMOD(dv[i][nv], dv[i-1][nv]);
         dvp[nv] = 0.5*dwh;
         vp[nv] = v[nv] + dvp[nv];
         vm[nv] = v[nv] - dvp[nv];
       }
       LpPROJECT(L, dvp, dwp);

       for (k = 0; k < NFLX; k++){
         if (nu[k] >= 0.0) 
           for (nv = 0; nv < NFLX; nv++) vp[nv] -= dwp[k]*nu[k]*R[nv][k];
         else     
           for (nv = 0; nv < NFLX; nv++) vm[nv] -= dwp[k]*nu[k]*R[nv][k];
       }
       #if NVAR != NFLX
        if (nu[NFLX] >= 0.0){  /* -- scalars all move at the flow speed -- */
          for (k = NFLX; k < NVAR; k++) vp[k] -= dwp[k]*nu[k];
        }else {
          for (k = NFLX; k < NVAR; k++) vm[k] -= dwp[k]*nu[k];
       }
       #endif 

       for (nv = NFLX; nv--;   ) {
         dwh = 0.5*delta_t*src[i][nv];
         vp[nv] += dwh;
         vm[nv] += dwh;
       }
       continue;
     }
    #endif  /* SHOCK_FLATTENING == MULTID */

   /* --------------------------------------------------------
       [Interpolation Step]

       Provide piecewise polynomial interpolation at time 
       level t(n) using characteristic variable differences.
       Interpolation can be done using LINEAR, WENO3 or 
       PARABOLIC algorithms.
       Here dwp and dwm are such that

         wp = w + dwp
         wm = w + dwm
      -------------------------------------------------------- */

    #if INTERPOLATION == LINEAR

   /* -- compute characteristic variable differences -- */

     LpPROJECT(L, dv[i-1], dwm1);
     LpPROJECT(L, dv[i  ], dwp1);

   /* -- slope limiter on chracteristics -- */

     for (k = 0; k < NVAR; k++){
       dwh = 0.0;
       if (dwp1[k]*dwm1[k] > 0.0) {
         dwh = 0.5*(dwm1[k] + dwp1[k]);
         d2w = kstp[k]*(fabs(dwp1[k]) < fabs(dwm1[k]) ? dwp1[k]:dwm1[k]);
         dwh = fabs(d2w) < fabs(dwh) ? d2w:dwh;
       }
       dwp[k] = 0.5*dwh;
     }

   /* -- project on right eigenvectors to obtain primite slopes -- */

     for (nv = 0; nv < NFLX; nv++) dvp[nv] = 0.0;
     for (k = 0; k < NFLX; k++){
       for (nv = 0; nv < NFLX; nv++) dvp[nv] += dwp[k]*R[nv][k];
     }
     #if NVAR != NFLX
      for (nv = NFLX; nv < NVAR; nv++) dvp[nv] = dwp[nv];
     #endif 

    /* -- compute L/R states at opposite interfaces at time t(n) -- */

     for (nv = 0; nv < NVAR; nv++) {
       vp[nv] = v[nv] + dvp[nv];
       vm[nv] = v[nv] - dvp[nv];
     }

    /* -- check positivity of density & pressure -- */

     if (vp[RHO] < 0.0 || vm[RHO] < 0.0) {
       dvp[RHO] = 0.5*(MINMOD(dv[i][RHO], dv[i-1][RHO]));
       vp[RHO] = v[RHO] + dvp[RHO];
       vm[RHO] = v[RHO] - dvp[RHO];
       LpPROJECT(L, dvp, dwp);
     }
     if (vp[PRS] < 0.0 || vm[PRS] < 0.0) {
       dvp[PRS] = 0.5*(MINMOD(dv[i][PRS], dv[i-1][PRS]));
       vp[PRS] = v[PRS] + dvp[PRS];
       vm[PRS] = v[PRS] - dvp[PRS];
       LpPROJECT(L, dvp, dwp);
     }

   /* -- Evolve by half time step using characteristic tracing -- */

     for (k = 0; k < NFLX; k++){
       if (nu[k] >= 0.0) 
         for (nv = 0; nv < NFLX; nv++) vp[nv] -= dwp[k]*nu[k]*R[nv][k];
       else     
         for (nv = 0; nv < NFLX; nv++) vm[nv] -= dwp[k]*nu[k]*R[nv][k];
     }
     #if NVAR != NFLX
      if (nu[NFLX] >= 0.0){  /* -- scalars all move at the flow speed -- */
        for (k = NFLX; k < NVAR; k++) vp[k] -= dwp[k]*nu[k];
      }else {
        for (k = NFLX; k < NVAR; k++) vm[k] -= dwp[k]*nu[k];
      }
     #endif 
     for (nv = NFLX; nv--;   ) {
       dwh = 0.5*delta_t*src[i][nv];
       vp[nv] += dwh;
       vm[nv] += dwh;
     }
     continue;

    #endif     /* INTERPOLATION == LINEAR */

    #if INTERPOLATION == WENO3

   /* -- compute undivided differences and 
         reconstruct characteristic fields -- */

     LpPROJECT(L, dv[i-1], dwm1);
     LpPROJECT(L, dv[i  ], dwp1);
     for (k = 0; k < NVAR; k++){
       tau = (dwp1[k] - dwm1[k]); 
       tau = tau*tau;

       a0 = 1.0 + tau/(dx2 + dwp1[k]*dwp1[k]);
       a1 = 1.0 + tau/(dx2 + dwm1[k]*dwm1[k]);

       dwp[k] =  (a0*dwp1[k] + 0.5*a1*dwm1[k])/(2.0*a0 + a1);
       dwm[k] = -(a1*dwm1[k] + 0.5*a0*dwp1[k])/(2.0*a1 + a0);
     }

    #endif     /* INTERPOLATION == WENO3 */

    #if INTERPOLATION == PARABOLIC

   /* -- compute undivided differences and 
         reconstruct characteristic fields -- */

     LpPROJECT(L, dv[i-2], dwm2);
     LpPROJECT(L, dv[i-1], dwm1);
     LpPROJECT(L, dv[i  ], dwp1);
     LpPROJECT(L, dv[i+1], dwp2);

     for (k = 0; k < NVAR; k++){  
       Smm = MC(dwm2[k], dwm1[k]);
       Smp = MC(dwm1[k], dwp1[k]);
       Spp = MC(dwp1[k], dwp2[k]);
       dwp[k] =  0.5*dwp1[k] - (Spp - Smp)*one_sixth;
       dwm[k] = -0.5*dwm1[k] - (Smp - Smm)*one_sixth;

      /* -- parabolic limiter-- */

       #if PARABOLIC_LIM == 0 || PARABOLIC_LIM == 2
        if (dwp[k]*dwm[k] >= 0.0) dwm[k] = dwp[k] = 0.0;
        else{
          if      (fabs(dwp[k]) >= 2.0*fabs(dwm[k])) dwp[k] = -2.0*dwm[k];
          else if (fabs(dwm[k]) >= 2.0*fabs(dwp[k])) dwm[k] = -2.0*dwp[k];
        }
       #endif
     }

    #endif   /* INTERPOLATION == PARABOLIC */

   /* ----------------------------------------------
       [Positivity Step]

       Check that left and right states at time
       t^n are physically admissible. 
       If not, use linear reconstruction on density 
       and pressure.
      ---------------------------------------------- */

   /* -- compute total jump in primitive variables -- */

     for (nv = 0; nv < NFLX; nv++) dvp[nv] = dvm[nv] = 0.0;
     for (k = 0; k < NFLX; k++){
       for (nv = 0; nv < NFLX; nv++){
         dvp[nv] += dwp[k]*R[nv][k];
         dvm[nv] += dwm[k]*R[nv][k];
       }
     }
     #if NVAR != NFLX
      for (nv = NFLX; nv < NVAR; nv++){
        dvp[nv] = dwp[nv];
        dvm[nv] = dwm[nv];
      }
     #endif 

    /* -- compute L/R states at opposite interfaces at time t(n) -- */

     for (nv = 0; nv < NVAR; nv++) {
       #if PARABOLIC_LIM == 1 || PARABOLIC_LIM == 2
        if (dvp[nv]*dvm[nv] >= 0.0) dvp[nv] = dvm[nv] = 0.0;
        else {
          if      (fabs(dvp[nv]) >= 2.0*fabs(dvm[nv])) dvp[nv] = -2.0*dvm[nv];
          else if (fabs(dvm[nv]) >= 2.0*fabs(dvp[nv])) dvm[nv] = -2.0*dvp[nv];
        }
       #endif
       vp[nv] = v[nv] + dvp[nv];
       vm[nv] = v[nv] + dvm[nv];
     }

     #if PARABOLIC_LIM == 1 || PARABOLIC_LIM == 2
      LpPROJECT(L, dvp, dwp);
      LpPROJECT(L, dvm, dwm);
     #endif

    /* -- check positivity of density & pressure -- */

     if (vp[RHO] < 0.0 || vm[RHO] < 0.0) {
       dvp[RHO] = 0.5*(MINMOD(dv[i][RHO], dv[i-1][RHO]));
       dvm[RHO] = - dvp[RHO]; 
       vp[RHO] = v[RHO] + dvp[RHO];
       vm[RHO] = v[RHO] + dvm[RHO];
       LpPROJECT(L, dvp, dwp);
       LpPROJECT(L, dvm, dwm);
     }
     if (vp[PRS] < 0.0 || vm[PRS] < 0.0) {
       dvp[PRS] = 0.5*(MINMOD(dv[i][PRS], dv[i-1][PRS]));
       dvm[PRS] = - dvp[PRS];       
       vp[PRS] = v[PRS] + dvp[PRS];
       vm[PRS] = v[PRS] + dvm[PRS];
       LpPROJECT(L, dvp, dwp);
       LpPROJECT(L, dvm, dwm);
     }

   /* ----------------------------------------------------------
       [Characteristic Tracing Step]
  
       evolve vp and vm previously defined by dt/2: project
       the characteristic differences dwp and dwm onto right 
       eigenvectors and discard contributions from waves not
       reaching the interface [only for WENO3 & PARABOLIC].
       The reference state is somewhat arbitrary and 
       can be set using REF_STATE:

       REF_STATE==1: use cell-center value;
       REF_STATE==2: interpolated value at base time level 
       REF_STATE==3: traditional PPM reference state (fastest
                     wave), minimize the size of the term 
                     subject to characteristic limiting. 

       Passive scalars use always REF_STATE == 2.
      ---------------------------------------------------------- */

    #if REF_STATE == 1
     for (nv = NFLX; nv--;   ) vp[nv] = vm[nv] = v[nv];
    #elif REF_STATE == 3
     nup = MAX(nu[1], 0.0); num = MIN(nu[0], 0.0);
     for (nv = NFLX; nv--;   ){
       dwh = vp[nv] - vm[nv];
       d2w = vp[nv] + vm[nv] - 2.0*v[nv];
       vp[nv] -= 0.5*nup*(dwh + d2w*(3.0 - 2.0*nup));
       vm[nv] -= 0.5*num*(dwh - d2w*(3.0 + 2.0*num));
     }
    #endif

    for (k = 0; k < NFLX; k++){  /* -- characteristic tracing step -- */
      dwh = dwp[k] - dwm[k];
      d2w = dwp[k] + dwm[k]; 
      if (nu[k] >= 0.0) {
        #if REF_STATE == 1
         Spp = dwp[k] - 0.5*nu[k]*(dwh + d2w*(3.0 - 2.0*nu[k]));
        #elif REF_STATE == 2
         Spp =        - 0.5*nu[k]*(dwh + d2w*(3.0 - 2.0*nu[k]));
        #elif REF_STATE == 3
         Spp = -0.5*(nu[k]-nup)*(dwh + 3.0*d2w) + (nu[k]*nu[k] - nup*nup)*d2w;
        #endif
        for (nv = NFLX; nv--;   ) vp[nv] += Spp*R[nv][k];
      } else {
        #if REF_STATE == 1
         Smm = dwm[k] - 0.5*nu[k]*(dwh - d2w*(3.0 + 2.0*nu[k]));
        #elif REF_STATE == 2
         Smm =        - 0.5*nu[k]*(dwh - d2w*(3.0 + 2.0*nu[k]));
        #elif REF_STATE == 3
         Smm = -0.5*(nu[k]-num)*(dwh - 3.0*d2w) + (nu[k]*nu[k] - num*num)*d2w;
        #endif
        for (nv = NFLX; nv--;   ) vm[nv] += Smm*R[nv][k];
      }
    }

    /* -- add source term -- */

    for (nv = NFLX; nv--;   ){   
      dwh = 0.5*delta_t*src[i][nv];
      vp[nv] += dwh;
      vm[nv] += dwh;
    }
    #if NVAR != NFLX
     if (nu[NFLX] >= 0.0) {   /* -- scalars all move at the flow speed -- */
       for (k = NFLX; k < NVAR; k++){
         dwh = dwp[k] - dwm[k];
         d2w = dwp[k] + dwm[k]; 
         vp[k] -= 0.5*nu[k]*(dwh + d2w*(3.0 - 2.0*nu[k]));
       }
     }else{
       for (k = NFLX; k < NVAR; k++){
         dwh = dwp[k] - dwm[k];
         d2w = dwp[k] + dwm[k]; 
         vm[k] -= 0.5*nu[k]*(dwh - d2w*(3.0 + 2.0*nu[k]));
       }
     }
    #endif 
  }  /* -- end main loop on grid points -- */


/*  -------------------------------------------
      Assign face-centered magnetic field
    -------------------------------------------  */

  #ifdef STAGGERED_MHD
   for (i = beg-1; i <= end; i++) {
     state->vR[i][B1] = state->vL[i][B1] = state->bn[i];
   }
  #endif

  CHECK_PRIM_STATES (state->vm, state->vp, state->v, beg, end);

/* -----------------------------------
      evolve center value by dt/2 
   ----------------------------------- */

  for (i = beg; i <= end; i++) {
    vp = state->vp[i];
    vm = state->vm[i];
    for (nv = NVAR; nv--; ) {
      state->vh[i][nv] = 0.5*(vp[nv] + vm[nv]);
    }
  }

  PRIMTOCON (state->vh, state->uh, beg, end);
  PRIMTOCON (state->vp, state->up, beg, end);
  PRIMTOCON (state->vm, state->um, beg, end);
}
#undef PARABOLIC_LIM
#undef REF_STATE 
