/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Compute the thermal conduction flux.

  Compute the thermal conduction flux along one row of computational 
  zones for the HD and MHD  modules according to Spitzer (1962):
  \f[
     F_c = F_{\rm class}/(|F_{\rm class}| + q)
  \f]
  where \f$ F_{\rm class} \f$ is the classical thermal conduction flux, 
  \f$ q \f$ is the saturated flux.
  Since the first term is purely parabolic, it is discretized using standard
  finite difference. 
  The saturated flux is treated in an upwind manner following the guidelines
  given in Balsara (2008) (see also Mignone et al. 2012)

  The classical MHD flux further splits into 2 components, along and 
  across the magnetic field lines.
  This function also computes the inverse of the time
  step and return its maximum over the current sweep.
  
  \b References
     - "The PLUTO Code for Adaptive Mesh Computations in Astrophysical 
        Fluid Dynamics" \n
       Mignone et al, ApJS (2012) 198, 7M
       
  \authors A. Mignone (mignone@ph.unito.it)\n
           T. Matsakos
  \date   Sep 1, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/*! When set to YES, saturated flux is computed using an upwind 
    selection rule.
    When set to NO, staurated flux is treated in the same manner as 
    the  conduction flux.                                          */ 
#define HYPERBOLIC_SAT_FLUX YES 

/* ********************************************************************* */
void TC_Flux (double ***T, const State_1D *state,  
              double **dcoeff, int beg, int end, Grid *grid)
/*! 
 * Compute the thermal conduction flux, state->par_flx.
 *
 * \param [in]     T       3D array containing the dimensionless 
 *                         temperature
 * \param [in,out] state   pointer to a State_1D structure
 * \param [out]    dcoeff  the diffusion coefficient needed for computing
 *                         the time step.
 * \param [in]      beg   initial index of computation
 * \param [in]      end   final   index of computation
 * \param [in]      grid  pointer to an array of Grid structures
 *
 * \return This function has no return value.
 *                       
 *********************************************************************** */
{
  int  i, j, k, nv;
  double bgradT, Bmag, dTmag;
  double Fc, Fcmag, Fsat, g1;
  double x1, x2, x3;
  double alpha, uL, uR, suL, suR, bn;
  double vi[NVAR], kpar=0.0, knor=0.0, phi;
  double bck_fld[3];
  static double **gradT;

  if (gradT == NULL) {
    gradT = ARRAY_2D(NMAX_POINT, 3, double);
  }

  GetGradient (T, gradT, beg, end, grid);
  D_EXPAND(x1 = grid[IDIR].x[g_i]; ,
           x2 = grid[JDIR].x[g_j]; ,
           x3 = grid[KDIR].x[g_k]; )

/* ----------------------------------------------- 
    Compute Thermal Conduction Flux (tcflx).
   ----------------------------------------------- */

  for (i = beg; i <= end; i++){
    
    for (nv = 0; nv < NVAR; nv++) {
      vi[nv]  = 0.5*(state->vh[i][nv] + state->vh[i+1][nv]);
    }

  /* ---------------------------------------------
      obtain the thermal conduction coefficients
      along (kpar) and across (knor) the field
      lines. This is done at the cell interface.
     --------------------------------------------- */

    if (g_dir == IDIR) x1 = grid[IDIR].xr[i];
    if (g_dir == JDIR) x2 = grid[JDIR].xr[i];
    if (g_dir == KDIR) x3 = grid[KDIR].xr[i];

    #if PHYSICS == MHD && BACKGROUND_FIELD == YES
     BackgroundField(x1,x2,x3, bck_fld);
     EXPAND(vi[BX1] += bck_fld[0];  ,
            vi[BX2] += bck_fld[1];  ,
            vi[BX3] += bck_fld[2];)
    #endif
    
    TC_kappa(vi, x1, x2, x3, &kpar, &knor, &phi);
    dTmag  = D_EXPAND(  gradT[i][0]*gradT[i][0], 
                      + gradT[i][1]*gradT[i][1], 
                      + gradT[i][2]*gradT[i][2]);
    dTmag = sqrt(dTmag) + 1.e-12;

  /* ---------------------------------------------------------
      compute magnitude of saturated flux using a Roe average
     --------------------------------------------------------- */

    #if HYPERBOLIC_SAT_FLUX == YES
/*
     uL = state->vL[i][PRS]/g1; suL = sqrt(uL);
     uR = state->vR[i][PRS]/g1; suR = sqrt(uR);
     scrh = 5.0*phi*g1*sqrt(g1/vi[RHO]);
     csat[i] = scrh*(uR + uL + suR*suL)/(suL + suR);
     Fsat    = 0.5*scrh*(uL*suL + uR*suR);
     if      (gradT[i][g_dir] > 0.0) Fsat += 0.5*csat[i]*(uR - uL);
     else if (gradT[i][g_dir] < 0.0) Fsat -= 0.5*csat[i]*(uR - uL);

*/
     uL = state->vL[i][PRS];
     uR = state->vR[i][PRS];
     if      (gradT[i][g_dir] > 0.0) Fsat = 5.0*phi/sqrt(vi[RHO])*uR*sqrt(uR);
     else if (gradT[i][g_dir] < 0.0) Fsat = 5.0*phi/sqrt(vi[RHO])*uL*sqrt(uL);
     else                          Fsat = 5.0*phi*vi[PRS]*sqrt(vi[PRS]/vi[RHO]);
    #else
     Fsat = 5.0*phi*vi[PRS]*sqrt(vi[PRS]/vi[RHO]);
    #endif

    #if PHYSICS == HD 

     Fc = kpar*gradT[i][g_dir];  /* -- classical thermal conduction flux -- */
/*
{
 double x = kpar*dTmag/Fsat;
 double flx_lim  = exp(-x*x);
 state->par_flx[i][ENG] =   flx_lim*Fc 
                          + (1.0 - flx_lim)*Fsat*gradT[i][g_dir]/dTmag;
 dcoeff[i][ENG] = flx_lim*kpar/vi[RHO]*(g_gamma - 1.0);
}
*/
     alpha                  = Fsat/(Fsat + kpar*dTmag);
     state->tc_flux[i][ENG] = alpha*Fc;
     dcoeff[i][ENG]         = fabs(alpha*kpar/vi[RHO])*(g_gamma - 1.0);

    #elif PHYSICS == MHD

     Bmag = EXPAND(vi[BX1]*vi[BX1], + vi[BX2]*vi[BX2], + vi[BX3]*vi[BX3]);
     Bmag = sqrt(Bmag) + 1.e-12;

     bgradT = D_EXPAND(  vi[BX1]*gradT[i][0], 
                       + vi[BX2]*gradT[i][1], 
                       + vi[BX3]*gradT[i][2]);
     bgradT /= Bmag;

  /* ------------------------------------------------------
      compute the classical MHD thermal conduction flux Fc
     ------------------------------------------------------ */

     bn    = vi[BX + g_dir]/Bmag; /* -- unit vector component -- */
     Fc    = kpar*bgradT*bn + knor*(gradT[i][g_dir] - bn*bgradT);
     Fcmag = sqrt((kpar*kpar - knor*knor)*bgradT*bgradT + 
                   knor*knor*dTmag*dTmag);     

     alpha                  = Fsat/(Fsat + Fcmag);
     state->tc_flux[i][ENG] = alpha*Fc;
     dcoeff[i][ENG]  = fabs(Fcmag/dTmag*bn*alpha/vi[RHO])*(g_gamma - 1.0);
    #endif
  }

}
#undef HYPERBOLIC_SAT_FLUX 

