/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute the right hand side of the conservative 
         HD/MHD equations.

  This function constructs the one-dimensional right hand side of 
  the conservative MHD or HD equations in the direction given by ::g_dir 
  in different geometries.
  The right hand side is computed as a two-point flux difference 
  term plus a source term:
  
  \f[ \mathrm{RHS}_i = 
      \frac{\Delta t}{\Delta V_i}\
      \Big(A_{i+1/2}F_{i+1/2} - A_{i-1/2}F_{i-1/2}\Big)  + \Delta t S_i  \f]
 
   where 
 
   - \f$ A_{i\pm 1/2} \f$ : interface areas
   - \f$ dV_i         \f$ : cell volume
   - \f$ F_{i\pm 1/2} \f$ : interface fluxes
   - \f$ dt           \f$ : time step
   - \f$ S_i          \f$ : source term including geometrical terms and 
                            body forces.
  
  See also the \ref RHS_page.
  The right hand side is assembled through the following steps:
 
   - If either one of FARGO, ROTATION or gravitational potential is 
      used, fluxes are combined to enforce conservation of total angular
      momentum and/or energy, see TotalFlux()
   - initialize rhs with flux differences
   - add geometrical source terms
   - enforce conservation of total angular  momentum and/or energy;
   - add gravity
   - add dissipative effects (viscosity and thermal conduction) to
     entropy equation. Ohmic dissipation is included in a separate step.

  For the entropy equation, the dissipative contributions can be recovered
  from the internal energy equation, for which (Boyd, Eq. [3.27]):
 
  \f[
      \pd{p}{t} + \vec{v}\cdot\nabla p + \Gamma p \nabla\cdot\vec{v}
    = (\Gamma-1)\left[\nabla\cdot\left(\kappa\nabla T\right)
      + \eta\vec{J}\cdot\vec{J} + \mu\left(\pd{v_i}{x_j}+\pd{v_j}{x_i}
        - \frac{2}{3}\delta_{ij}\nabla\cdot\vec{v}\right)\pd{v_i}{x_j}\right]
  \f]
 
  To obtain the corresponding contribution to the (conservative form of)
  entropy equation, just divide the previous one by 
   \f$ \rho^{\Gamma-1}\f$:
 
  \f[
     \pd{(\rho s)}{t} + \nabla\cdot(\rho s\vec{v}) = 
     (\Gamma-1)\rho^{1-\Gamma}\left[\nabla\cdot\left(\kappa\nabla T\right)
       + \eta\vec{J}^2 + \mu \Pi_{ij}\pd{v_i}{x_j}\right]
  \f]
      
  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 5, 2014

  \b References
    - Goedbloed & Poedts, "Principle of MHD", page 165-166
    - Boyd, "The Physics of Plasmas", page 57, Eq. 3.27
  

*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

static void TotalFlux (const State_1D *, double *, int, int, Grid *);

#ifdef CH_SPACEDIM   /*  implies Chombo is being used  */
 #define USE_PR_GRADIENT  YES   
#else
 #ifdef FINITE_DIFFERENCE 
  #define USE_PR_GRADIENT  NO   /* -- default for Finite Difference schemes -- */
 #else
  #define USE_PR_GRADIENT  YES   /* -- default, do not change!! -- */
 #endif
#endif

#if defined FARGO && !defined SHEARINGBOX
 #define IF_FARGO(a)  a
#else 
 #define IF_FARGO(a)  
#endif

#if ROTATING_FRAME == YES
 #define IF_ROTATION(a)  a
#else 
 #define IF_ROTATION(a)  
#endif

#if PHYSICS == MHD
#if BACKGROUND_FIELD == YES
 #define TotBB(v, b0, a, b) (v[a]*(v[b] + b0[b]) + b0[a]*v[b])
#else
 #define TotBB(v, b0, a, b) (v[a]*v[b])
#endif
#endif

/* *********************************************************************** */
void RightHandSide (const State_1D *state, Time_Step *Dts, 
                    int beg, int end, double dt, Grid *grid)
/*! 
 *
 * \param [in,out]  state  pointer to State_1D structure
 * \param [in]      Dts    pointer to time step structure
 * \param [in]      beg    initial index of computation
 * \param [in]      end    final   index of computation
 * \param [in]      dt     time increment
 * \param [in]      grid  pointer to Grid structure
 *
 * \return This function has no return value.
 * \note    --
 * \todo    --
 ************************************************************************* */
{
  int    i, j, k, nv;
  double dtdx, dtdV, scrh, rhog;
  double *x1, *x1p, *dx1, *dV1;
  double *x2, *x2p, *dx2, *dV2;
  double *x3, *x3p, *dx3, *dV3;
  double **vh, **vp, **vm;
  double **flux, **rhs, *p;
  double *A, *dV;
  double cl;
  double **Bg0, **wA, w, wp, vphi, phi_c;
  double g[3];
  static double **fA, *phi_p;
  #if ENTROPY_SWITCH == YES
   double rhs_entr;
   double **visc_flux, **visc_src, **tc_flux, **res_flux;
   static double **fvA;
  #endif

  #if GEOMETRY != CARTESIAN
   if (fA == NULL) fA = ARRAY_2D(NMAX_POINT, NVAR, double);
  #endif
  #ifdef FARGO
   wA = FARGO_GetVelocity();
  #endif

  if (phi_p == NULL) phi_p = ARRAY_1D(NMAX_POINT, double);

  #if (defined SHEARINGBOX) && (defined FARGO) && (EOS == IDEAL)
   print1 ("! ShearingBox+Fargo+Ideal EoS not properly implemented\n");
   QUIT_PLUTO(1);
  #endif

/* --------------------------------------------------
             Compute passive scalar fluxes
   -------------------------------------------------- */

  #if NSCL > 0
   AdvectFlux (state, beg - 1, end, grid);
  #endif

  #if (PHYSICS == MHD) && (BACKGROUND_FIELD == YES)
   Bg0 = GetBackgroundField (beg, end, CELL_CENTER, grid);
  #endif

/* --------------------------
      pointer shortcuts
   -------------------------- */

  rhs    = state->rhs;
  flux   = state->flux;
  p      = state->press;
  vh     = state->vh;
  vp     = state->vp;
  vm     = state->vm;
  A      = grid[g_dir].A;
  dV     = grid[g_dir].dV;
  
  #if ENTROPY_SWITCH == YES
   visc_flux = state->visc_flux; 
   visc_src  = state->visc_src; 
   tc_flux   = state->tc_flux; 
   res_flux  = state->res_flux;
   if (fvA == NULL) fvA = ARRAY_2D(NMAX_POINT, NVAR, double);
  #endif 

  x1  = grid[IDIR].x;  x2  = grid[JDIR].x;  x3  = grid[KDIR].x;
  x1p = grid[IDIR].xr; x2p = grid[JDIR].xr; x3p = grid[KDIR].xr;
  dx1 = grid[IDIR].dx; dx2 = grid[JDIR].dx; dx3 = grid[KDIR].dx; 
  dV1 = grid[IDIR].dV; dV2 = grid[JDIR].dV; dV3 = grid[KDIR].dV;

  i = g_i;  /* will be redefined during x1-sweep */
  j = g_j;  /* will be redefined during x2-sweep */
  k = g_k;  /* will be redefined during x3-sweep */
  
/* ------------------------------------------------
     Add pressure to normal component of 
     momentum flux if necessary.
   ------------------------------------------------ */

  #if USE_PR_GRADIENT == NO
   for (i = beg - 1; i <= end; i++) flux[i][MXn] += p[i];
  #endif

/* -------------------------------------------------
    compute gravitational potential and add 
    contribution to energy flux.
   ------------------------------------------------- */

  #if (defined FARGO && !defined SHEARINGBOX) ||\
      (ROTATING_FRAME == YES) || (BODY_FORCE & POTENTIAL)
   TotalFlux(state, phi_p, beg-1, end, grid);
  #endif

#if GEOMETRY == CARTESIAN
/* ***********************************************************
   
        Compute right-hand side of the MHD/HD 
        equations in CARTESIAN geometry.
    
   *********************************************************** */
{
  double x, y, z, *vc;

  if (g_dir == IDIR){

  /* ****************************************************
      Cartesian x-direction,

       - initialize rhs with flux differences (I1)
       - enforce conservation of total angular
         momentum and/or energy               (I3)
       - add gravity                          (I4)
     **************************************************** */

    y = x2[j];
    z = x3[k];
    for (i = beg; i <= end; i++) {
      x    = x1[i];
      dtdx = dt/dx1[i];

    /* -----------------------------------------------
       I1. initialize rhs with flux difference
       ----------------------------------------------- */

      VAR_LOOP(nv) rhs[i][nv] = -dtdx*(flux[i][nv] - flux[i-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[i][MX1] -= dtdx*(p[i] - p[i-1]);
      #endif

    /* ----------------------------------------------------
       I3. modify rhs to achieve conservation
       ---------------------------------------------------- */

      #if (defined FARGO && !defined SHEARINGBOX) 
       w = wA[g_k][i];
       rhs[i][MX2] -= w*rhs[i][RHO];
       #if HAVE_ENERGY
        rhs[i][ENG] -= w*(rhs[i][MX2] + 0.5*w*rhs[i][RHO]);
       #endif
      #endif

    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      vc = vh[i];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x, y, z);
       rhs[i][MX1] += dt*vc[RHO]*g[IDIR];
       #if DIMENSIONS == 1 && COMPONENTS > 1 /* In 1D and when COMPONENTS == 2  */
                                             /* or 3, add gravity contributions */
                                             /* along non-existing dimensions. */
        EXPAND(                                     ,
               rhs[i][MX2] += dt*vc[RHO]*g[JDIR];   ,
               rhs[i][MX3] += dt*vc[RHO]*g[KDIR];)
       #endif
       
       #if HAVE_ENERGY
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[IDIR];
        #if DIMENSIONS == 1 && COMPONENTS > 1
         rhs[i][ENG] += dt*(EXPAND(0.0, + vc[RHO]*vc[VX2]*g[JDIR],
                                        + vc[RHO]*vc[VX3]*g[KDIR]));
        #endif
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][MX1] -= dtdx*vc[RHO]*(phi_p[i] - phi_p[i-1]);
       #if HAVE_ENERGY
        phi_c        = BodyForcePotential(x, y, z);
        rhs[i][ENG] -= phi_c*rhs[i][RHO];
       #endif
      #endif

      #ifdef SHEARINGBOX
       rhs[i][MX1] += dt*2.0*vc[RHO]*vc[VX2]*sb_Omega;
      #endif

    /* ----------------------------------------------------
       I5. Add dissipative terms to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[i][ENG] - visc_flux[i-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(visc_flux[i][MX1] - visc_flux[i-1][MX1])  ,
                           + vc[VX2]*(visc_flux[i][MX2] - visc_flux[i-1][MX2])  ,
                           + vc[VX3]*(visc_flux[i][MX3] - visc_flux[i-1][MX3]));
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[i][ENG] - tc_flux[i-1][ENG]);       
       #endif
       rhs[i][ENTR] += rhs_entr*dtdx*rhog;              
      #endif

    }
  } else if (g_dir == JDIR){

  /* ****************************************************
      Cartesian y-direction,

       - initialize rhs with flux differences (J1)
       - add gravity                          (J4)
     **************************************************** */

    x = x1[i];
    z = x3[k];
    for (j = beg; j <= end; j++) {
      y    = x2[j];
      dtdx = dt/dx2[j];

    /* -----------------------------------------------
       J1. initialize rhs with flux difference
       ----------------------------------------------- */

      VAR_LOOP(nv) rhs[j][nv] = -dtdx*(flux[j][nv] - flux[j-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[j][MX2] -= dtdx*(p[j] - p[j-1]);
      #endif

    /* ----------------------------------------------------
       J4. Include gravity
       ---------------------------------------------------- */

      vc = vh[j];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[j], x3[k]);
       rhs[j][MX2] += dt*vc[RHO]*g[JDIR];
       #if DIMENSIONS == 2 && COMPONENTS == 3
        rhs[j][MX3] += dt*vc[RHO]*g[KDIR];
       #endif

       #if HAVE_ENERGY
        rhs[j][ENG] += dt*0.5*(flux[j][RHO] + flux[j-1][RHO])*g[JDIR];
        #if DIMENSIONS == 2 && COMPONENTS == 3
         rhs[j][ENG] += dt*vc[RHO]*vc[VX3]*g[KDIR];
        #endif
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[j][MX2] -= dtdx*vc[RHO]*(phi_p[j] - phi_p[j-1]);
       #if HAVE_ENERGY
        phi_c       = BodyForcePotential(x, y, z);
        rhs[j][ENG] -= phi_c*rhs[j][RHO];
       #endif
      #endif

      #ifdef SHEARINGBOX
       rhs[j][MX2] -= dt*2.0*vc[RHO]*vc[VX1]*sb_Omega;
      #endif

    /* ----------------------------------------------------
       J5. Add dissipative terms to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[j][ENG] - visc_flux[j-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(visc_flux[j][MX1] - visc_flux[j-1][MX1])  ,
                           + vc[VX2]*(visc_flux[j][MX2] - visc_flux[j-1][MX2])  ,
                           + vc[VX3]*(visc_flux[j][MX3] - visc_flux[j-1][MX3]));
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[j][ENG] - tc_flux[j-1][ENG]);       
       #endif
       rhs[j][ENTR] += rhs_entr*dtdx*rhog;              
      #endif

    }

  }else if (g_dir == KDIR){

  /* ****************************************************
      Cartesian z-direction,

       - initialize rhs with flux differences (K1)
       - enforce conservation of total angular
         momentum and/or energy               (K3)
       - add gravity                          (K4)
     **************************************************** */

    x = x1[i];
    y = x2[j];
    for (k = beg; k <= end; k++) {
      z    = x3[k];
      dtdx = dt/dx3[k];

    /* -----------------------------------------------
       K1. initialize rhs with flux difference
       ----------------------------------------------- */

      VAR_LOOP(nv) rhs[k][nv] = -dtdx*(flux[k][nv] - flux[k-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[k][MX3] -= dtdx*(p[k] - p[k-1]);
      #endif

    /* ----------------------------------------------------
       K3. modify rhs to achieve conservation
       ---------------------------------------------------- */

      #if (defined FARGO && !defined SHEARINGBOX) 
       w = wA[k][i];
       rhs[k][MX2] -= w*rhs[k][RHO];
       #if HAVE_ENERGY
        rhs[k][ENG] -= w*(rhs[k][MX2] + 0.5*w*rhs[k][RHO]);
       #endif
      #endif

    /* ----------------------------------------------------
       K4. Include gravity
       ---------------------------------------------------- */

      vc = vh[k];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[j], x3[k]);
       rhs[k][MX3] += dt*vc[RHO]*g[KDIR];
       #if HAVE_ENERGY
        rhs[k][ENG] += dt*0.5*(flux[k][RHO] + flux[k-1][RHO])*g[KDIR];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[k][MX3] -= dtdx*vc[RHO]*(phi_p[k] - phi_p[k-1]);
       #if HAVE_ENERGY
        phi_c        = BodyForcePotential(x, y, z);
        rhs[k][ENG] -= phi_c*rhs[k][RHO];
       #endif
      #endif

    /* ----------------------------------------------------
       K5. Add dissipative terms to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[k][ENG] - visc_flux[k-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(visc_flux[k][MX1] - visc_flux[k-1][MX1])  ,
                           + vc[VX2]*(visc_flux[k][MX2] - visc_flux[k-1][MX2])  ,
                           + vc[VX3]*(visc_flux[k][MX3] - visc_flux[k-1][MX3]));
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[k][ENG] - tc_flux[k-1][ENG]);
       #endif
       rhs[k][ENTR] += rhs_entr*dtdx*rhog;              
      #endif
    }
  }
}
#elif GEOMETRY == CYLINDRICAL

/* ***********************************************************
   
        Compute right-hand side of the MHD/HD 
        equations in CYLINDRICAL geometry.
    
   *********************************************************** */
{
  double R, z, phi, R_1; 

  if (g_dir == IDIR) {  
    double vc[NVAR];

  /* ****************************************************
      Cylindrical radial direction:
      multiply fluxes times interface area
     **************************************************** */

    z   = x2[j];
    phi = 0.0;
    for (i = beg - 1; i <= end; i++){ 
      R = grid[IDIR].A[i];

      fA[i][RHO] = flux[i][RHO]*R;
      EXPAND(fA[i][iMR]   = flux[i][iMR]*R;     ,
             fA[i][iMZ]   = flux[i][iMZ]*R;     ,
             fA[i][iMPHI] = flux[i][iMPHI]*R*R;)       
      #if (ENTROPY_SWITCH == YES) && (VISCOSITY == EXPLICIT)
       EXPAND(fvA[i][iMR]   = visc_flux[i][iMR]*R; ,
              fvA[i][iMZ]   = visc_flux[i][iMZ]*R; ,
              fvA[i][iMPHI] = visc_flux[i][iMPHI]*R*R;)
       #if HAVE_ENERGY
        fvA[i][ENG] = visc_flux[i][ENG]*R;
       #endif
      #endif
      #if PHYSICS == MHD
       EXPAND(fA[i][iBR]   = flux[i][iBR]*R;  ,
              fA[i][iBZ]   = flux[i][iBZ]*R;  ,
              fA[i][iBPHI] = flux[i][iBPHI]*R;)
      #endif
      #if HAVE_ENERGY
       fA[i][ENG] = flux[i][ENG]*R;
      #endif
      #ifdef GLM_MHD
       fA[i][PSI_GLM] = flux[i][PSI_GLM]*R;
      #endif
      for (nv = NFLX; nv < NVAR; nv++) fA[i][nv] = flux[i][nv]*R;
    }

  /* ****************************************************
      Cylindrical radial direction,

       - initialize rhs with flux differences (I1)
       - add source terms                     (I2)
       - add gravity                          (I4)
     **************************************************** */

    for (i = beg; i <= end; i++){ 
      R    = x1[i];
      dtdV = dt/dV1[i];
      dtdx = dt/dx1[i];
      R_1  = 1.0/R;

    /* ---------------------------------------------------------------
       I1. initialize rhs with flux difference

           Note: when there's no explicit resistivity, we use the 
                 formulation with source terms since it seems to have 
                 better stability properties.
       ---------------------------------------------------------------- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(rhs[i][iMR]   = - dtdV*(fA[i][iMR]   - fA[i-1][iMR]);  ,
             rhs[i][iMZ]   = - dtdV*(fA[i][iMZ]   - fA[i-1][iMZ]);  ,
             rhs[i][iMPHI] = - dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*fabs(R_1);)
      #if USE_PR_GRADIENT == YES
       rhs[i][iMR] -= dtdx*(p[i] - p[i-1]);  
      #endif
      #if PHYSICS == MHD
       EXPAND(rhs[i][iBR]   = - dtdV*(fA[i][iBR]   - fA[i-1][iBR]);  ,
              rhs[i][iBZ]   = - dtdV*(fA[i][iBZ]   - fA[i-1][iBZ]);  ,
              rhs[i][iBPHI] = - dtdx*(flux[i][iBPHI] - flux[i-1][iBPHI]);)
       #ifdef GLM_MHD
        rhs[i][iBR]     = - dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
        rhs[i][PSI_GLM] = - dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      #if HAVE_ENERGY
       rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);
      #endif

      for (nv = NFLX; nv < NVAR; nv++) {
        rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      }

    /* ----------------------------------------------------
       I2. Add source terms
       ---------------------------------------------------- */

      for (nv = NVAR; nv--;  ) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]); 
/* for (nv = NVAR; nv--;  ) vc[nv] = vh[i][nv]; */

      #if COMPONENTS == 3
       vphi = vc[iVPHI];
       #if ROTATING_FRAME == YES
        w     = g_OmegaZ*R;
        vphi += w;
       #endif

       #if PHYSICS == HD
        rhs[i][iMR] += dt*vc[RHO]*vphi*vphi*R_1;
       #elif PHYSICS == MHD
        rhs[i][iMR] += dt*(vc[RHO]*vphi*vphi - TotBB(vc, Bg0[i], iBPHI, iBPHI))*R_1;
       #endif /* PHYSICS == MHD */
      #endif  /* COMPONENTS == 3 */

    /* ----------------------------------------------------
       I3. modify rhs to achieve conservation
       ---------------------------------------------------- */

      #if ROTATING_FRAME == YES
       rhs[i][iMPHI] -= w*rhs[i][RHO];
       #if HAVE_ENERGY
        rhs[i][ENG]  -= w*(rhs[i][iMPHI] + 0.5*w*rhs[i][RHO]);
       #endif
      #endif

    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[j], x3[k]);
       rhs[i][iMR] += dt*vc[RHO]*g[IDIR];
       #if HAVE_ENERGY
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[IDIR];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMR] -= dtdx*vc[RHO]*(phi_p[i] - phi_p[i-1]);
       #if HAVE_ENERGY
        phi_c        = BodyForcePotential(R, z, phi);
        rhs[i][ENG] -= phi_c*rhs[i][RHO];
       #endif
      #endif

    /* ----------------------------------------------------
       I5. Add dissipative terms to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (fvA[i][ENG] - fvA[i-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(fvA[i][MX1] - fvA[i-1][MX1])      ,
                           + vc[VX2]*(fvA[i][MX2] - fvA[i-1][MX2])      ,
                           + vc[VX3]*(fvA[i][MX3] - fvA[i-1][MX3])*fabs(R_1));
                    
        rhs_entr -= EXPAND(  vc[VX1]*visc_src[i][MX1], 
                           + vc[VX2]*visc_src[i][MX2],
                           + vc[VX3]*visc_src[i][MX3]);
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (A[i]*tc_flux[i][ENG] - A[i-1]*tc_flux[i-1][ENG]);       
       #endif
       rhs[i][ENTR] += rhs_entr*dtdV*rhog;              
      #endif

    }
     
  } else if (g_dir == JDIR) { 
    double *vc;

  /* ****************************************************
      Cylindrical vertical direction:

       - initialize rhs with flux differences (J1)
       - add gravity                          (J4)
     **************************************************** */

    R   = x1[i];
    phi = 0.0;
    for (j = beg; j <= end; j++){ 
      z    = x2[j];   
      dtdx = dt/dx2[j];

    /* -----------------------------------------------
       J1. initialize rhs with flux difference
       ----------------------------------------------- */

      VAR_LOOP(nv) rhs[j][nv] = -dtdx*(flux[j][nv] - flux[j-1][nv]);
      #if USE_PR_GRADIENT == YES
       rhs[j][iMZ] += - dtdx*(p[j] - p[j-1]);
      #endif

    /* ----------------------------------------------------
       J4. Include gravity
       ---------------------------------------------------- */

      vc = vh[j];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[j], x3[k]);
       rhs[j][iMZ] += dt*vc[RHO]*g[JDIR];
       #if HAVE_ENERGY
        rhs[j][ENG] += dt*0.5*(flux[j][RHO] + flux[j-1][RHO])*g[JDIR];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[j][iMZ] += -dtdx*vc[RHO]*(phi_p[j] - phi_p[j-1]);
       #if HAVE_ENERGY
        phi_c       = BodyForcePotential(R, z, phi);
        rhs[j][ENG] -= phi_c*rhs[j][RHO];
       #endif
      #endif

    /* ----------------------------------------------------
       J5. Add dissipative terms to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[j][ENG] - visc_flux[j-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(visc_flux[j][MX1] - visc_flux[j-1][MX1])  ,
                           + vc[VX2]*(visc_flux[j][MX2] - visc_flux[j-1][MX2])  ,
                           + vc[VX3]*(visc_flux[j][MX3] - visc_flux[j-1][MX3]));
       #endif
       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[j][ENG] - tc_flux[j-1][ENG]);
       #endif
       rhs[j][ENTR] += rhs_entr*dtdx*rhog;              
      #endif

    }
  }
}

#elif GEOMETRY == POLAR

/* ***********************************************************
   
        Compute right-hand side of the MHD/HD 
        equations in POLAR geometry.
    
   *********************************************************** */
{
  double R, phi, z; 
  double R_1;
   
  if (g_dir == IDIR) { 
    double vc[NVAR];

  /* ****************************************************
      Polar radial direction:
      multiply fluxes times interface area
     **************************************************** */

    phi = x2[j];
    z   = x3[k];
    for (i = beg - 1; i <= end; i++) { 
      R = grid[IDIR].A[i];

      fA[i][RHO] = flux[i][RHO]*R;
      EXPAND(fA[i][iMR]   = flux[i][iMR]*R;      ,
             fA[i][iMPHI] = flux[i][iMPHI]*R*R;  ,
             fA[i][iMZ]   = flux[i][iMZ]*R;)
      #if (ENTROPY_SWITCH == YES) && (VISCOSITY == EXPLICIT)
       EXPAND(fvA[i][iMR]   = visc_flux[i][iMR]*R;     ,
              fvA[i][iMPHI] = visc_flux[i][iMPHI]*R*R; ,
              fvA[i][iMZ]   = visc_flux[i][iMZ]*R;)
       #if HAVE_ENERGY
        fvA[i][ENG] = visc_flux[i][ENG]*R;
       #endif
      #endif
      #if PHYSICS == MHD 
       EXPAND(fA[i][iBR]   = flux[i][iBR]*R;    ,
              fA[i][iBPHI] = flux[i][iBPHI];    ,
              fA[i][iBZ]   = flux[i][iBZ]*R;)
      #endif
      #if HAVE_ENERGY
       fA[i][ENG] = flux[i][ENG]*R;
      #endif
      #ifdef GLM_MHD
       fA[i][PSI_GLM] = flux[i][PSI_GLM]*R;
      #endif
      for (nv = NFLX; nv < NVAR; nv++) fA[i][nv] = flux[i][nv]*R;
    }

  /* ****************************************************
      Polar radial direction,

       - initialize rhs with flux differences (I1)
       - add source terms                     (I2)
       - enforce conservation of total angular
         momentum and/or energy               (I3)
       - add gravity                          (I4)
     **************************************************** */

    for (i = beg; i <= end; i++) {
      R    = x1[i];
      dtdV = dt/dV1[i];
      dtdx = dt/dx1[i];
      R_1  = grid[IDIR].r_1[i];

    /* -----------------------------------------------
       I1. initialize rhs with flux difference
       ----------------------------------------------- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(rhs[i][iMR]   = - dtdV*(fA[i][iMR]   - fA[i-1][iMR])
                             - dtdx*(p[i] - p[i-1]);                      ,      
             rhs[i][iMPHI] = - dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*R_1;  ,
             rhs[i][iMZ]   = - dtdV*(fA[i][iMZ]   - fA[i-1][iMZ]);)
      #if PHYSICS == MHD 
       EXPAND(rhs[i][iBR]   = -dtdV*(fA[i][iBR]   - fA[i-1][iBR]);    ,
              rhs[i][iBPHI] = -dtdx*(fA[i][iBPHI] - fA[i-1][iBPHI]);  ,
              rhs[i][iBZ]   = -dtdV*(fA[i][iBZ]   - fA[i-1][iBZ]);)
       #ifdef GLM_MHD
        rhs[i][iBR]     = -dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
        rhs[i][PSI_GLM] = -dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      #if HAVE_ENERGY
       rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);
      #endif
      for (nv = NFLX; nv < NVAR; nv++) {
        rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      }

    /* ----------------------------------------------------
       I2. Add source terms
       ---------------------------------------------------- */

      for (nv = NVAR; nv--;  ) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]); 
/* for (nv = NVAR; nv--;  ) vc[nv] = vh[i][nv]; */
      vphi = vc[iVPHI];
      #if (defined FARGO) || (ROTATING_FRAME == YES)
       w = 0.0;
       IF_FARGO   (w += wA[g_k][i];)
       IF_ROTATION(w += g_OmegaZ*R;)
       vphi += w;
      #endif
      #if PHYSICS == HD
       rhs[i][iMR] += dt*vc[RHO]*vphi*vphi*R_1;
      #elif PHYSICS == MHD
       rhs[i][iMR] += dt*(  vc[RHO]*vphi*vphi 
                          - TotBB(vc, Bg0[i], iBPHI, iBPHI))*R_1;
      #endif

    /* ----------------------------------------------------
       I3. modify rhs to achieve conservation
       ---------------------------------------------------- */

      #if (defined FARGO) || (ROTATING_FRAME == YES)
       rhs[i][iMPHI] -= w*rhs[i][RHO];
       #if HAVE_ENERGY
        rhs[i][ENG]   -= w*(rhs[i][iMPHI] + 0.5*w*rhs[i][RHO]);
       #endif
      #endif
      
    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[j], x3[k]);
       rhs[i][iMR] += dt*vc[RHO]*g[IDIR];
       #if HAVE_ENERGY
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[IDIR];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMR] -= dtdx*vc[RHO]*(phi_p[i] - phi_p[i-1]);
       #if HAVE_ENERGY
        phi_c       = BodyForcePotential(R, phi, z);
        rhs[i][ENG] -= phi_c*rhs[i][RHO];
       #endif
      #endif

    /* ----------------------------------------------------
       I5. Add dissipative terms to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (fvA[i][ENG] - fvA[i-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(fvA[i][iMR]   - fvA[i-1][iMR])        ,
                           + vc[VX2]*(fvA[i][iMPHI] - fvA[i-1][iMPHI])*R_1  ,
                           + vc[VX3]*(fvA[i][iMZ]   - fvA[i-1][iMZ]));
                    
        rhs_entr -= EXPAND(  vc[VX1]*visc_src[i][MX1], 
                           + vc[VX2]*visc_src[i][MX2],
                           + vc[VX3]*visc_src[i][MX3]);
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (A[i]*tc_flux[i][ENG] - A[i-1]*tc_flux[i-1][ENG]);       
       #endif
       rhs[i][ENTR] += rhs_entr*dtdV*rhog;              
      #endif

    }
     
  } else if (g_dir == JDIR) {
    double *vc;

  /* ****************************************************
      Polar azimuthal direction:

       - initialize rhs with flux differences (J1)
       - add gravity                          (J4)
     **************************************************** */

    R = x1[i];
    z = x3[k];
    scrh = dt/R;
    for (j = beg; j <= end; j++){ 
      phi  = x2[j];
      dtdx = scrh/dx2[j];

    /* ------------------------------------------------
       J1. Compute equations rhs for phi-contributions
       ------------------------------------------------ */

      VAR_LOOP(nv) rhs[j][nv] = -dtdx*(flux[j][nv] - flux[j-1][nv]);
      rhs[j][iMPHI] -= dtdx*(p[j] - p[j-1]);

    /* -------------------------------------------------------
       J4. Include gravity
       ------------------------------------------------------- */

      vc = vh[j];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[j], x3[k]);
       rhs[j][iMPHI] += dt*vc[RHO]*g[JDIR];
       #if HAVE_ENERGY
        rhs[j][ENG] += dt*0.5*(flux[j][RHO] + flux[j-1][RHO])*g[JDIR];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[j][iMPHI] -= dtdx*vc[RHO]*(phi_p[j] - phi_p[j-1]);
       #if HAVE_ENERGY
        phi_c        = BodyForcePotential(R, phi, z);
        rhs[j][ENG] -= phi_c*rhs[j][RHO];
       #endif
      #endif

    /* ----------------------------------------------------
       J5. Add dissipative terms to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[j][ENG] - visc_flux[j-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(visc_flux[j][MX1] - visc_flux[j-1][MX1]) ,
                           + vc[VX2]*(visc_flux[j][MX2] - visc_flux[j-1][MX2]) ,
                           + vc[VX3]*(visc_flux[j][MX3] - visc_flux[j-1][MX3]));
                    
        rhs_entr -= EXPAND(  vc[VX1]*visc_src[j][MX1], 
                           + vc[VX2]*visc_src[j][MX2],
                           + vc[VX3]*visc_src[j][MX3]);
       #endif
       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[j][ENG] - tc_flux[j-1][ENG]);
       #endif
       rhs[j][ENTR] += rhs_entr*dtdx*rhog;              
      #endif
    }

  } else if (g_dir == KDIR) { 
    double *vc;

  /* ****************************************************
      Polar vertical direction:

       - initialize rhs with flux differences (K1)
       - enforce conservation of total angular
         momentum and/or energy               (K3)
       - add gravity                          (K4)
     **************************************************** */

    R   = x1[i];
    phi = x2[j];
    for (k = beg; k <= end; k++){ 
      z    = x3[k];
      dtdx = dt/dx3[k];

    /* -----------------------------------------------
       K1. initialize rhs with flux difference
       ----------------------------------------------- */

      VAR_LOOP(nv) rhs[k][nv] = -dtdx*(flux[k][nv] - flux[k-1][nv]);
      rhs[k][iMZ] -= dtdx*(p[k] - p[k-1]);

    /* ------------------------------------------------------
       K3. modify rhs to enforce conservation (FARGO only)
           (solid body rotations are not included the
            velocity depends on the cylindrical radius only)
       ------------------------------------------------------ */

      #ifdef FARGO 
       w = wA[k][i];
       rhs[k][iMPHI] -= w*rhs[k][RHO];
       #if HAVE_ENERGY
        rhs[k][ENG]  -= w*(rhs[k][iMPHI] + 0.5*w*rhs[k][RHO]);       
       #endif
      #endif

    /* ----------------------------------------------------
       K4. Include gravity
       ---------------------------------------------------- */

      vc = vh[k];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[j], x3[k]); 
       rhs[k][iMZ] += dt*vc[RHO]*g[KDIR];
       #if HAVE_ENERGY
        rhs[k][ENG] += dt*0.5*(flux[k][RHO] + flux[k-1][RHO])*g[KDIR];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[k][iMZ] += -dtdx*vc[RHO]*(phi_p[k] - phi_p[k-1]);
       #if HAVE_ENERGY
        phi_c       = BodyForcePotential(R, phi, z);
        rhs[k][ENG] -= phi_c*rhs[k][RHO];
       #endif
      #endif

    /* ----------------------------------------------------
       K5. Add dissipative terms to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[k][ENG] - visc_flux[k-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(visc_flux[k][MX1] - visc_flux[k-1][MX1])  ,
                           + vc[VX2]*(visc_flux[k][MX2] - visc_flux[k-1][MX2])  ,
                           + vc[VX3]*(visc_flux[k][MX3] - visc_flux[k-1][MX3]));
       #endif

       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[k][ENG] - tc_flux[k-1][ENG]);
       #endif
       rhs[k][ENTR] += rhs_entr*dtdx*rhog;              
      #endif

    }
  }
}
#elif GEOMETRY == SPHERICAL

/* ***********************************************************
   
        Compute right-hand side of the MHD/HD 
        equations in SPHERICAL geometry.
    
   *********************************************************** */
{
  double r, th, phi;
  double r2, r3, r_1;
  double s, s2, ct, s_1;

  if (g_dir == IDIR) { 
    double Sm, vc[NVAR];

  /* ****************************************************
      Spherical radial direction: 
      multiply fluxes by interface area 
     **************************************************** */

    th  = x2[j]; 
    s   = sin(th);
    phi = x3[k];
    for (i = beg - 1; i <= end; i++){
      r  = x1p[i];
      r2 = r*r; 
      r3 = r2*r;

      fA[i][RHO] = flux[i][RHO]*r2;
      EXPAND(fA[i][iMR]   = flux[i][iMR]*r2;   ,
             fA[i][iMTH]  = flux[i][iMTH]*r2;  ,
             fA[i][iMPHI] = flux[i][iMPHI]*r3;)
      #if (ENTROPY_SWITCH == YES) && (VISCOSITY == EXPLICIT)
       EXPAND(fvA[i][iMR]   = visc_flux[i][iMR]*r2;    ,
              fvA[i][iMTH]  = visc_flux[i][iMTH]*r2;   ,
              fvA[i][iMPHI] = visc_flux[i][iMPHI]*r3;)
       #if HAVE_ENERGY
        fvA[i][ENG] = visc_flux[i][ENG]*r2;
       #endif
      #endif
      #if PHYSICS == MHD
       EXPAND(fA[i][iBR]   = flux[i][iBR]*r2;   ,
              fA[i][iBTH]  = flux[i][iBTH]*r;  ,
              fA[i][iBPHI] = flux[i][iBPHI]*r;)
      #endif
      #if HAVE_ENERGY
       fA[i][ENG] = flux[i][ENG]*r2;
      #endif
      #ifdef GLM_MHD
       fA[i][PSI_GLM] = flux[i][PSI_GLM]*r2;
      #endif
      for (nv = NFLX; nv < NVAR; nv++) fA[i][nv] = flux[i][nv]*r2;
    } 

  /* ****************************************************
      Spherical radial direction:

       - initialize rhs with flux differences (I1)
       - add source terms                     (I2)
       - enforce conservation of total angular 
         momentum and/or energy               (I3)
       - add gravity                          (I4)
     **************************************************** */

    for (i = beg; i <= end; i++) { 
      r    = x1[i];
      dtdV = dt/dV1[i];
      dtdx = dt/dx1[i];
      r_1  = grid[IDIR].r_1[i];

    /* -----------------------------------------------
       I1. initialize rhs with flux difference
       ----------------------------------------------- */

      rhs[i][RHO] = -dtdV*(fA[i][RHO] - fA[i-1][RHO]);
      EXPAND(
        rhs[i][iMR]   = - dtdV*(fA[i][iMR] - fA[i-1][iMR])
                        - dtdx*(p[i] - p[i-1]);                    ,
        rhs[i][iMTH]  = -dtdV*(fA[i][iMTH]  - fA[i-1][iMTH]);      ,
        rhs[i][iMPHI] = -dtdV*(fA[i][iMPHI] - fA[i-1][iMPHI])*r_1; 
      )
      #if PHYSICS == MHD
       EXPAND(                                                     
         rhs[i][iBR]   = -dtdV*(fA[i][iBR]   - fA[i-1][iBR]);       ,
         rhs[i][iBTH]  = -dtdx*(fA[i][iBTH]  - fA[i-1][iBTH])*r_1;  ,
         rhs[i][iBPHI] = -dtdx*(fA[i][iBPHI] - fA[i-1][iBPHI])*r_1;
       )
       #ifdef GLM_MHD
        rhs[i][iBR]     = -dtdx*(flux[i][iBR]   - flux[i-1][iBR]);
        rhs[i][PSI_GLM] = -dtdV*(fA[i][PSI_GLM] - fA[i-1][PSI_GLM]);
       #endif
      #endif
      #if HAVE_ENERGY
       rhs[i][ENG] = -dtdV*(fA[i][ENG] - fA[i-1][ENG]);
      #endif

      for (nv = NFLX; nv < NVAR; nv++){
        rhs[i][nv] = -dtdV*(fA[i][nv] - fA[i-1][nv]);
      }

    /* ----------------------------------------------------
       I2. Add source terms 
       ---------------------------------------------------- */
  
      for (nv = NVAR; nv--;  ) vc[nv] = 0.5*(vp[i][nv] + vm[i][nv]);
/*  for (nv = NVAR; nv--;  ) vc[nv] = vh[i][nv]; */

      vphi = SELECT(0.0, 0.0, vc[iVPHI]);
      #if (defined FARGO) || (ROTATING_FRAME == YES)
       w = 0.0;
       IF_FARGO   (w += wA[g_j][i];)
       IF_ROTATION(w += g_OmegaZ*r*s;)
       vphi += w;
      #endif
      Sm = vc[RHO]*(EXPAND(  0.0, + vc[iVTH]*vc[iVTH], + vphi*vphi));
      #if PHYSICS == MHD
       Sm += EXPAND(  0.0, - TotBB(vc, Bg0[i], iBTH, iBTH), 
                           - TotBB(vc, Bg0[i], iBPHI,iBPHI));
      #endif
      rhs[i][iMR] += dt*Sm*r_1;

    /* ----------------------------------------------------
       I3. modify rhs to enforce conservation
       ---------------------------------------------------- */

      #if (defined FARGO) || (ROTATING_FRAME == YES)
       rhs[i][iMPHI] -= w*rhs[i][RHO];
       #if HAVE_ENERGY
        rhs[i][ENG]   -= w*(rhs[i][iMPHI] + 0.5*w*rhs[i][RHO]);
       #endif
      #endif

    /* ----------------------------------------------------
       I4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[j], x3[k]); 
       rhs[i][iMR] += dt*vc[RHO]*g[IDIR];
       #if HAVE_ENERGY
        rhs[i][ENG] += dt*0.5*(flux[i][RHO] + flux[i-1][RHO])*g[IDIR]; 
       #endif
      #endif
      
      #if (BODY_FORCE & POTENTIAL)
       rhs[i][iMR] -= dtdx*vc[RHO]*(phi_p[i] - phi_p[i-1]);
       #if HAVE_ENERGY
        phi_c       = BodyForcePotential(r, th, phi); 
        rhs[i][ENG] -= phi_c*rhs[i][RHO];
       #endif
      #endif

    /* ----------------------------------------------------
       I5. Add dissipative terms to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (fvA[i][ENG] - fvA[i-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(fvA[i][iMR]   - fvA[i-1][iMR])        ,
                           + vc[VX2]*(fvA[i][iMTH]  - fvA[i-1][iMTH])       ,
                           + vc[VX3]*(fvA[i][iMPHI] - fvA[i-1][iMPHI])*r_1);
                    
        rhs_entr -= EXPAND(  vc[VX1]*visc_src[i][MX1], 
                           + vc[VX2]*visc_src[i][MX2],
                           + vc[VX3]*visc_src[i][MX3]);
       #endif
       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (A[i]*tc_flux[i][ENG] - A[i-1]*tc_flux[i-1][ENG]);       
       #endif
       rhs[i][ENTR] += rhs_entr*dtdV*rhog;              
      #endif
    }

  } else if (g_dir == JDIR) {

    double Sm, *vc;

  /* ****************************************************
      Spherical meridional direction:
      multiply fluxes by zone-interface area
     **************************************************** */

    r   = x1[i];
    phi = x3[k];
    
    for (j = beg - 1; j <= end; j++){ 
      s  = grid[JDIR].A[j];
      s2 = s*s;

      fA[j][RHO] = flux[j][RHO]*s;
      EXPAND(fA[j][iMR]   = flux[j][iMR]*s;   ,
             fA[j][iMTH]  = flux[j][iMTH]*s;  ,
             fA[j][iMPHI] = flux[j][iMPHI]*s2;) 

      #if (ENTROPY_SWITCH == YES) && (VISCOSITY == EXPLICIT)
       EXPAND(fvA[j][iMR]   = visc_flux[j][iMR]*s;    ,
              fvA[j][iMTH]  = visc_flux[j][iMTH]*s;   ,
              fvA[j][iMPHI] = visc_flux[j][iMPHI]*s2;)
       #if HAVE_ENERGY
        fvA[j][ENG] = visc_flux[j][ENG]*s;
       #endif
      #endif

      #if PHYSICS == MHD
       EXPAND(fA[j][iBR]   = flux[j][iBR]*s;   ,
              fA[j][iBTH]  = flux[j][iBTH]*s;  ,
              fA[j][iBPHI] = flux[j][iBPHI];)
      #endif
      #if HAVE_ENERGY
       fA[j][ENG] = flux[j][ENG]*s;
      #endif
      #ifdef GLM_MHD
       fA[j][PSI_GLM] = flux[j][PSI_GLM]*s;
      #endif
      for (nv = NFLX; nv < NVAR; nv++) fA[j][nv] = flux[j][nv]*s;
    }

  /* ****************************************************
      Spherical meridional direction:

       - initialize rhs with flux differences (J1)
       - add source terms                     (J2)
       - enforce conservation of total angular
         momentum and/or energy               (J3)
       - add gravity                          (J4)
     **************************************************** */
    
//    r_1 = grid[IDIR].r_1[i];
    r_1 = 0.5*(x1p[i]*x1p[i] - x1p[i-1]*x1p[i-1])/dV1[i];
    
    for (j = beg; j <= end; j++){
      th   = x2[j];
      dtdV = dt/dV2[j]*r_1;
      dtdx = dt/dx2[j]*r_1;      
      s    = sin(th);
      s_1  = 1.0/s;   
      ct   = grid[JDIR].ct[j];         /* = cot(theta)  */

    /* -----------------------------------------------
       J1. initialize rhs with flux difference
       ----------------------------------------------- */

      rhs[j][RHO] = -dtdV*(fA[j][RHO] - fA[j-1][RHO]);
      EXPAND(
        rhs[j][iMR]   = - dtdV*(fA[j][iMR] - fA[j-1][iMR]);  , 
        rhs[j][iMTH]  = - dtdV*(fA[j][iMTH] - fA[j-1][iMTH])
                        - dtdx*(p[j] - p[j-1]);              , 
        rhs[j][iMPHI] = - dtdV*(fA[j][iMPHI] - fA[j-1][iMPHI])*fabs(s_1);
      )       
      #if PHYSICS == MHD
       EXPAND(                                                     
         rhs[j][iBR]   = -dtdV*(fA[j][iBR]   - fA[j-1][iBR]);  ,
         rhs[j][iBTH]  = -dtdV*(fA[j][iBTH]  - fA[j-1][iBTH]);  ,
         rhs[j][iBPHI] = -dtdx*(fA[j][iBPHI] - fA[j-1][iBPHI]);
       )
       #ifdef GLM_MHD
        rhs[j][iBTH]    = -dtdx*(flux[j][iBTH] - flux[j-1][iBTH]);
        rhs[j][PSI_GLM] = -dtdV*(fA[j][PSI_GLM] - fA[j-1][PSI_GLM]);
       #endif
      #endif
      #if HAVE_ENERGY
       rhs[j][ENG] = -dtdV*(fA[j][ENG] - fA[j-1][ENG]);
      #endif
      for (nv = NFLX; nv < NVAR; nv++){
        rhs[j][nv] = -dtdV*(fA[j][nv] - fA[j-1][nv]);
      }

    /* ----------------------------------------------------
       J2. Add source terms
       ---------------------------------------------------- */
       
      vc = vh[j];
      vphi = SELECT(0.0, 0.0, vc[iVPHI]);
      #if (defined FARGO) || (ROTATING_FRAME == YES)
       w = 0.0; 
       IF_FARGO   (w += wA[j][i];)
       IF_ROTATION(w += g_OmegaZ*r*s;)
       vphi += w;
      #endif
      Sm = vc[RHO]*(EXPAND(  0.0, - vc[iVTH]*vc[iVR], + ct*vphi*vphi));
      #if PHYSICS == MHD
       Sm += EXPAND(0.0, +    TotBB(vc, Bg0[j], iBTH, iBR), 
                         - ct*TotBB(vc, Bg0[j], iBPHI, iBPHI));
      #endif
      rhs[j][iMTH] += dt*Sm*r_1;

    /* ----------------------------------------------------
       J3. modify rhs to enforce conservation
       ---------------------------------------------------- */

      #if (defined FARGO) || (ROTATING_FRAME == YES)
       rhs[j][iMPHI] -= w*rhs[j][RHO];
       #if HAVE_ENERGY
        rhs[j][ENG]  -= w*(rhs[j][iMPHI] + 0.5*w*rhs[j][RHO]);
       #endif
      #endif

    /* ----------------------------------------------------
       J4. Include gravity
       ---------------------------------------------------- */

      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[j], x3[k]);
       rhs[j][iMTH] += dt*vc[RHO]*g[JDIR];
       #if HAVE_ENERGY
        rhs[j][ENG] += dt*0.5*(flux[j][RHO] + flux[j-1][RHO])*g[JDIR];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[j][iMTH] -= dtdx*vc[RHO]*(phi_p[j] - phi_p[j-1]);
       #if HAVE_ENERGY
        phi_c        = BodyForcePotential(r, th, phi); 
        rhs[j][ENG] -= phi_c*rhs[j][RHO];
       #endif
      #endif

    /* ----------------------------------------------------
       J5. Add TC dissipative term to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (fvA[j][ENG] - fvA[j-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(fvA[j][iMR]   - fvA[j-1][iMR])        ,
                           + vc[VX2]*(fvA[j][iMTH]  - fvA[j-1][iMTH])       ,
                           + vc[VX3]*(fvA[j][iMPHI] - fvA[j-1][iMPHI])*fabs(s_1));
                    
        rhs_entr -= EXPAND(  vc[VX1]*visc_src[j][MX1], 
                           + vc[VX2]*visc_src[j][MX2],
                           + vc[VX3]*visc_src[j][MX3]);
       #endif
       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (A[j]*tc_flux[j][ENG] - A[j-1]*tc_flux[j-1][ENG]);
       #endif
       rhs[j][ENTR] += rhs_entr*dtdV*rhog;              
      #endif

    }

  } else if (g_dir == KDIR) {

    double Sm, *vc;

  /* ****************************************************
      Spherical azimuthal direction:

       - initialize rhs with flux differences (K1)
       - add gravity                          (K4)
     **************************************************** */

    r    = x1[i];
    th   = x2[j];

    r_1  = 0.5*(x1p[i]*x1p[i] - x1p[i-1]*x1p[i-1])/dV1[i];
    scrh = dt*r_1*dx2[j]/dV2[j];
    for (k = beg; k <= end; k++) {
      phi  = x3[k];
      dtdx = scrh/dx3[k];

    /* ------------------------------------------------
       K1.  initialize rhs with flux difference
       ------------------------------------------------ */

      VAR_LOOP(nv) rhs[k][nv] = -dtdx*(flux[k][nv] - flux[k-1][nv]);
      rhs[k][iMPHI] -= dtdx*(p[k] - p[k-1]); 

    /* -------------------------------------------------------
       K4. Include gravity
       ------------------------------------------------------- */

      vc = vh[k];
      #if (BODY_FORCE & VECTOR)
       BodyForceVector(vc, g, x1[i], x2[j], x3[k]);
       rhs[k][iMPHI] += dt*vc[RHO]*g[KDIR];
       #if HAVE_ENERGY
        rhs[k][ENG] += dt*0.5*(flux[k][RHO] + flux[k-1][RHO])*g[KDIR];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       rhs[k][iMPHI] -= dtdx*vc[RHO]*(phi_p[k] - phi_p[k-1]);
       #if HAVE_ENERGY
        phi_c        = BodyForcePotential(r, th, phi); 
        rhs[k][ENG] -= phi_c*rhs[k][RHO];
       #endif
      #endif

    /* ----------------------------------------------------
       K5. Add dissipative terms to entropy equation
       ---------------------------------------------------- */

      #if (ENTROPY_SWITCH == YES) && (PARABOLIC_FLUX & EXPLICIT)
       rhog = vc[RHO];
       rhog = (g_gamma - 1.0)*pow(rhog,1.0-g_gamma);
 
       rhs_entr = 0.0;
       #if VISCOSITY == EXPLICIT
        rhs_entr += (visc_flux[k][ENG] - visc_flux[k-1][ENG]);
        rhs_entr -= EXPAND(  vc[VX1]*(visc_flux[k][MX1] - visc_flux[k-1][MX1])  ,
                           + vc[VX2]*(visc_flux[k][MX2] - visc_flux[k-1][MX2])  ,
                           + vc[VX3]*(visc_flux[k][MX3] - visc_flux[k-1][MX3]));
       #endif
       #if THERMAL_CONDUCTION == EXPLICIT
        rhs_entr += (tc_flux[k][ENG] - tc_flux[k-1][ENG]);
       #endif
       rhs[k][ENTR] += rhs_entr*dtdx*rhog;              
      #endif
    }
  }
}
#endif  /* GEOMETRY == SPHERICAL */

/* ---------------------------------------------------------------
    Source terms coming from tensor discretazion of parabolic
    terms in curvilinear coordinates (only for viscosity)
  ---------------------------------------------------------------- */

  #if GEOMETRY != CARTESIAN
   #if VISCOSITY == EXPLICIT
    for (i = beg; i <= end; i++) {
      EXPAND(rhs[i][MX1] += dt*state->visc_src[i][MX1];  ,
             rhs[i][MX2] += dt*state->visc_src[i][MX2];  ,
             rhs[i][MX3] += dt*state->visc_src[i][MX3];)
    }
   #endif
  #endif

/* --------------------------------------------------
              Powell's source terms
   -------------------------------------------------- */

  #if PHYSICS == MHD 
   #if MHD_FORMULATION == EIGHT_WAVES
    for (i = beg; i <= end; i++) {
      EXPAND(rhs[i][MX1] += dt*state->src[i][MX1];  ,
             rhs[i][MX2] += dt*state->src[i][MX2];  ,
             rhs[i][MX3] += dt*state->src[i][MX3];)

      EXPAND(rhs[i][BX1] += dt*state->src[i][BX1];  ,
             rhs[i][BX2] += dt*state->src[i][BX2];  ,
             rhs[i][BX3] += dt*state->src[i][BX3];)
      #if HAVE_ENERGY
       rhs[i][ENG] += dt*state->src[i][ENG];
      #endif
    }
   #endif
  #endif

/* -------------------------------------------------
            Extended GLM source terms
   ------------------------------------------------- */

  #if (defined GLM_MHD) && (EGLM == YES)
   EGLM_Source (state, dt, beg, end, grid);
  #endif

/* --------------------------------------------------
    Reset right hand side in internal boundary zones
   -------------------------------------------------- */
   
  #if INTERNAL_BOUNDARY == YES
   InternalBoundaryReset(state, Dts, beg, end, grid);
  #endif
  
/* --------------------------------------------------
           Time step determination
   -------------------------------------------------- */

#if !GET_MAX_DT
return;
#endif

  cl = 0.0;
  for (i = beg-1; i <= end; i++) {
    scrh = Dts->cmax[i]*grid[g_dir].inv_dxi[i];
    cl = MAX(cl, scrh);   
  }
  #if GEOMETRY == POLAR || GEOMETRY == SPHERICAL
   if (g_dir == JDIR) {
     cl /= fabs(grid[IDIR].xgc[g_i]);
   }
   #if GEOMETRY == SPHERICAL
    if (g_dir == KDIR){
      cl /= fabs(grid[IDIR].xgc[g_i])*sin(grid[JDIR].xgc[g_j]);
    }
   #endif
  #endif
  Dts->inv_dta = MAX(cl, Dts->inv_dta);  
}

/* ********************************************************************* */
void TotalFlux (const State_1D *state, double *phi_p,
                int beg, int end, Grid *grid)
/*!
 *  Compute the total flux in order to enforce conservation of 
 *  angular momentum and energy in presence of FARGO source 
 *  terms, rotation or gravitational potential.
 *
 * \param [in]     state pointer to State_1D structure;
 * \param [in,out] phi_p  1D array defining the gravitational potential;
 * \param [in]     beg    initial index of computation; 
 * \param [in]     end    final   index of computation;
 * \param [in]     grid  pointer to Grid structure;
 *********************************************************************** */
#ifndef iMPHI
 #define iMPHI MY  /* -- for Cartesian coordinates -- */
#endif
{
  int i;
  double wp, R;
  double **flux, *vp;
  double *x1,  *x2,  *x3;
  double *x1p, *x2p, *x3p;
  #ifdef FARGO
   double **wA;
   wA = FARGO_GetVelocity();
  #endif

  flux = state->flux;
  x1  = grid[IDIR].x;  x1p = grid[IDIR].xr;
  x2  = grid[JDIR].x;  x2p = grid[JDIR].xr;
  x3  = grid[KDIR].x;  x3p = grid[KDIR].xr;

  if (g_dir == IDIR){ 
    for (i = beg; i <= end; i++){

    /* ----------------------------------------------------
        include flux contributions from FARGO or Rotation 
        Note: ShearingBox terms are not included here but
              only within the BodyForce function.
       ---------------------------------------------------- */

      #if (defined FARGO && !defined SHEARINGBOX) || (ROTATING_FRAME == YES)
       wp = 0.0;
       #if GEOMETRY == SPHERICAL
        IF_FARGO(wp = 0.5*(wA[g_j][i] + wA[g_j][i+1]);)
        R = x1p[i]*sin(x2[g_j]);  /* -- cylindrical radius -- */
       #else
        IF_FARGO(wp = 0.5*(wA[g_k][i] + wA[g_k][i+1]);)
        R = x1p[i];                   /* -- cylindrical radius -- */
       #endif
       IF_ROTATION(wp += g_OmegaZ*R;)
       #if HAVE_ENERGY
        flux[i][ENG] += wp*(0.5*wp*flux[i][RHO] + flux[i][iMPHI]);
       #endif
       flux[i][iMPHI] += wp*flux[i][RHO];
      #endif

    /* -- gravitational potential -- */

      #if (BODY_FORCE & POTENTIAL)
       phi_p[i] = BodyForcePotential(x1p[i], x2[g_j], x3[g_k]);
       #if HAVE_ENERGY
        flux[i][ENG] += flux[i][RHO]*phi_p[i];                          
       #endif
      #endif
    }
  }else if (g_dir == JDIR){
    for (i = beg; i <= end; i++){ 

    /* ----------------------------------------------------
        include flux contributions from FARGO and Rotation 
       ---------------------------------------------------- */

      #if GEOMETRY == SPHERICAL
       #if defined FARGO || (ROTATING_FRAME == YES)
        wp = 0.0;
        R  = x1[g_i]*sin(x2p[i]);
        IF_FARGO   (wp += 0.5*(wA[i][g_i] + wA[i+1][g_i]);)
        IF_ROTATION(wp += g_OmegaZ*R;)
        #if HAVE_ENERGY
         flux[i][ENG] += wp*(0.5*wp*flux[i][RHO] + flux[i][iMPHI]);
        #endif
        flux[i][iMPHI] += wp*flux[i][RHO];
       #endif
      #endif

      #if (BODY_FORCE & POTENTIAL)
       phi_p[i] = BodyForcePotential(x1[g_i], x2p[i], x3[g_k]);
       #if HAVE_ENERGY
        flux[i][ENG] += flux[i][RHO]*phi_p[i];
       #endif
      #endif      

    }
  }else if (g_dir == KDIR){
    R = x1[g_i];
    for (i = beg; i <= end; i++) {

    /* ----------------------------------------------------
        include flux contributions from FARGO
        (polar/cartesian geometries only)
       ---------------------------------------------------- */

      #if (GEOMETRY != SPHERICAL) && (defined FARGO) && (!defined SHEARINGBOX)
       wp = 0.5*(wA[i][g_i] + wA[i+1][g_i]);
       #if HAVE_ENERGY
        flux[i][ENG] += wp*(0.5*wp*flux[i][RHO] + flux[i][iMPHI]);
       #endif
       flux[i][iMPHI] += wp*flux[i][RHO];
      #endif

      #if (BODY_FORCE & POTENTIAL)
       phi_p[i] = BodyForcePotential(x1[g_i], x2[g_j], x3p[i]); 
       #if HAVE_ENERGY
        flux[i][ENG] += flux[i][RHO]*phi_p[i];                          
       #endif
      #endif
    }
  }
}
#undef TotBB
#undef IF_FARGO
#undef IF_ROTATION
