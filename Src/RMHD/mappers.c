/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to 
  an array of conservative variables for the RMHD equations.
  
  The ConsToPrim() converts an array of conservative quantities to 
  an array of primitive quantities.
  During the conversion, pressure is normally recovered from total 
  energy using the algorithm outlined in
  - "Equation of state in relativistic magnetohydrodynamics: variable versus
     constant adiabatic index"\n
     Mignone \& Mc Kinney, MNRAS (2007) 378, 1118.

  However, if the zone has been tagged with FLAG_ENTROPY, primitive
  variables are recovered by using the conserved entropy rather than
  total energy.
  
  In other words:
  
      if (FLAG_ENTROPY is TRUE)  --> p = p(S)
      else                       --> p = p(E)
  
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   Oct 4, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void PrimToCons (double **uprim, double **ucons, int ibeg, int iend)
/*!
 * Convert primitive variables to conservative variables. 
 *
 * \param [in]  uprim array of primitive variables
 * \param [out] ucons array of conservative variables
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final index of computation
 *
 *********************************************************************** */
{
  int   i, nv;
  double  vel2, vB, Bmag2;
  double  g, g2, wt;
  double  *u, *v;
  static double *h;
  #if EOS == IDEAL
   double gmmr = g_gamma/(g_gamma - 1.0);
  #endif

  if (h == NULL) h = ARRAY_1D(NMAX_POINT, double);

  Enthalpy(uprim, h, ibeg, iend);

  for (i = ibeg; i <= iend; i++) {

    v = uprim[i];
    u = ucons[i];

    vel2  = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
    vB    = EXPAND(v[VX1]*v[BX1], + v[VX2]*v[BX2], + v[VX3]*v[BX3]);
    Bmag2 = EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);

    g2 = 1.0/(1.0 - vel2);
    g  = sqrt(g2);

    wt  = v[RHO]*h[i]*g2 + Bmag2;

  /* -------------------------------------------------------
       Convert from primitive (v) to conservative (u)   
     ------------------------------------------------------- */

    u[RHO] = g*v[RHO];
    EXPAND (u[MX1] = wt*v[VX1] - vB*v[BX1];  ,
            u[MX2] = wt*v[VX2] - vB*v[BX2];  ,
            u[MX3] = wt*v[VX3] - vB*v[BX3];)

    EXPAND (u[BX1] = v[BX1];  ,
            u[BX2] = v[BX2];  ,
            u[BX3] = v[BX3];)

    #if RESISTIVE_RMHD == YES
     v[EC] = u[EC];
     EXPAND (u[EX] = v[EX];  ,
             u[EY] = v[EY];  ,
             u[EZ] = v[EZ];)
    #endif

    #if SUBTRACT_DENSITY == YES
     #if EOS == IDEAL
      u[ENG]  = v[PRS]*(g2*gmmr - 1.0) + u[RHO]*g2*vel2/(g + 1.0)
               + 0.5*(Bmag2*(1.0 + vel2) - vB*vB);
     #elif EOS == TAUB
      wt    = v[PRS]/v[RHO];
      u[ENG] = v[PRS]*(g2*2.5 - 1.0) 
               + u[RHO]*g2*(2.25*wt*wt + vel2)/(g*sqrt(1.0 + 2.25*wt*wt) + 1.0)
               + 0.5*(Bmag2*(1.0 + vel2) - vB*vB);
     #endif
    #else
     u[ENG]  = v[RHO]*h[i]*g2 - v[PRS] + 0.5*(Bmag2*(1.0 + vel2) - vB*vB);
    #endif

    for (nv = NFLX; nv < (NFLX + NSCL); nv++) u[nv] = u[RHO]*v[nv];
    #ifdef GLM_MHD
     u[PSI_GLM] = v[PSI_GLM]; 
    #endif
  }
}

/* ********************************************************************* */
int ConsToPrim (double **ucons, double **uprim, int ibeg, int iend, 
                unsigned char *flag)
/*!
 * Convert from conservative to primitive variables.
 *
 * \param [in]  ucons  array of conservative variables
 * \param [out] uprim  array of primitive variables
 * \param [in]  beg    starting index of computation
 * \param [in]  end    final index of computation
 * \param [out] flag   array of flags tagging zones where conversion
 *                     went wrong.
 * 
 * \return Return (0) if conversion was succesful in every zone 
 *         [ibeg,iend]. 
 *         Otherwise, return a non-zero integer number giving the bit 
 *         flag(s) turned on during the conversion process.
 *         In this case, flag contains the failure codes of those
 *         zones where where conversion did not go through.
 *
 *********************************************************************** */
{
  int    i, nv, status=0;
  int    use_energy;
  double *u, *v, scrh, w_1;
  Map_param par;

  for (i = ibeg; i <= iend; i++) {

    flag[i] = 0;
    u = ucons[i];
    v = uprim[i];

/* ------------------------------------------------------------
      Define the input parameters of the parameter structure
   ------------------------------------------------------------ */

    par.D  = u[DN];
    par.E  = u[EN];
    par.S  = EXPAND(u[MX1]*u[BX1], + u[MX2]*u[BX2], + u[MX3]*u[BX3]);
    par.m2 = EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]);
    par.B2 = EXPAND(u[BX1]*u[BX1], + u[BX2]*u[BX2], + u[BX3]*u[BX3]); 
    par.S2 = par.S*par.S;

  /* -------------------------------------------
        Check density and energy positivity 
     ------------------------------------------- */
  
    if (u[RHO] < 0.0) {
      print("! ConsToPrim: negative density (%8.2e), ", u[RHO]);
      Where (i, NULL);
      u[RHO]   = g_smallDensity;
      flag[i] |= RHO_FAIL;
    }

    if (u[ENG] < 0.0) {
      WARNING(
        print("! ConsToPrim: negative energy (%8.2e), ", u[ENG]);
        Where (i, NULL);
      )
      u[ENG]   = 1.e-5;
      flag[i] |= ENG_FAIL;
    }

  /* --------------------------------
      recover pressure by inverting 
      the energy or entropy equation.
     -------------------------------- */

    use_energy = 1; /* -- default -- */
    #if ENTROPY_SWITCH == YES
     par.sigma_c = u[ENTR];
     #ifdef CH_SPACEDIM
      if (g_intStage > 0)   /* -- hot fix to be used with Chombo: avoid calling 
                             CheckZone when writing file to disk          -- */
     #endif
     if (CheckZone(i,FLAG_ENTROPY)) {
       use_energy = 0;
       if (EntropySolve(&par) != 0) {
         flag[i] |= ENG_FAIL;
         WARNING(Where (i, NULL);)
         if (PressureFix(&par) != 0) flag[i] |= PRS_FAIL;
         u[ENG] = par.E;  /* redefine total energy */
       }
     } 
    #endif
 
    if (use_energy){
      if (EnergySolve (&par) != 0){
        WARNING(Where(i,NULL);)
        flag[i] |= ENG_FAIL;
        if (PressureFix(&par) != 0) flag[i] |= PRS_FAIL;
/* printf ("! PressureFix: p = %12.6e, lor = %12.6e\n", par.p, par.lor); */
        u[ENG] = par.E;  /* redefine total energy */
      }
    }

  /* --------------------------------
      quit if a consistent pressure 
      value could not be found.  
     -------------------------------- */
    
    if (flag[i] & PRS_FAIL){
      Where(i,NULL);
      QUIT_PLUTO(1);
    }

 /* -----------------------------------------------------
      W, p and lor have been found. 
      Now complete conversion 
    ----------------------------------------------------- */

    v[RHO] = u[RHO]/par.lor;
    v[PRS] = par.prs;

    w_1  = 1.0/(par.W + par.B2);
    scrh = par.S/par.W;
    EXPAND(v[VX1] = w_1*(u[MX1] + scrh*u[BX1]);  ,
           v[VX2] = w_1*(u[MX2] + scrh*u[BX2]);  ,
           v[VX3] = w_1*(u[MX3] + scrh*u[BX3]);)

    scrh = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
    if (scrh >= 1.0){
      print ("! v^2 = %f > 1 in mappers (p = %12.6e)\n", scrh, par.prs);
      QUIT_PLUTO(1);
      flag[i] |= RHO_FAIL;
    }

    EXPAND(v[BX1] = u[BX1];  ,
           v[BX2] = u[BX2];  ,
           v[BX3] = u[BX3];)

    for (nv = NFLX; nv < NVAR; nv++) v[nv] = u[nv]/u[RHO];
    #ifdef GLM_MHD
     v[PSI_GLM] = u[PSI_GLM]; 
    #endif

    status |= flag[i];
  }
  return(status);
}
