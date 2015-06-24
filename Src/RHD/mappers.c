/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to 
  an array of conservative variables for the RHD equations.
  
  The ConsToPrim() converts an array of conservative quantities to 
  an array of primitive quantities.
  During the conversion, pressure is normally recovered from total 
  energy unless zone has been tagged with FLAG_ENTROPY.
  In this case we recover pressure from conserved entropy:
  
      if (FLAG_ENTROPY is TRUE)  --> p = p(S)
      else                       --> p = p(E)
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   Oct 4, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void PrimToCons  (double *uprim[], double *ucons[],
                 int ibeg, int iend)
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
  int     nv, ii;
  double  rhoh_g2, scrh, g;
  double  beta_fix = 0.9999;
  double  *u, *v;
  static double  *h;

  if (h == NULL){
    h = ARRAY_1D(NMAX_POINT, double);
  }

  Enthalpy (uprim, h, ibeg, iend);

  for (ii = ibeg; ii <= iend; ii++) {
   
    u = ucons[ii];
    v = uprim[ii];

    g = EXPAND(v[VX1]*v[VX1], +v[VX2]*v[VX2], +v[VX3]*v[VX3]);

    #if USE_FOUR_VELOCITY == YES
     g       = sqrt(1.0 + g);
     scrh    = v[RHO]*h[ii]*g;
     rhoh_g2 = scrh*g;
    #else

     if (g >= 1.0){
       WARNING( 
         print ("! u^2 > 1 (%f) in PrimToCons\n", scrh);
         Where (ii, NULL);
       )

       g = beta_fix/sqrt(g);
       EXPAND(v[VX1] *= g;  ,
              v[VX2] *= g;  ,
              v[VX3] *= g;)

       g = beta_fix*beta_fix;
     }
     g    = 1.0/(1.0 - g);
     scrh = rhoh_g2 = v[RHO]*h[ii]*g;
     g    = sqrt(g);

    #endif

    u[RHO] = v[RHO]*g;
    EXPAND(u[MX1] = scrh*v[VX1];  ,
           u[MX2] = scrh*v[VX2];  ,
           u[MX3] = scrh*v[VX3];)

    ucons[ii][ENG] = rhoh_g2 - v[PRS];

    for (nv = NFLX; nv < (NFLX + NSCL); nv++) u[nv] = v[nv]*u[RHO];

    #if CHECK_CONSERVATIVE_VAR == YES
     m2 = EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]);
     g  = u[ENG] - sqrt(m2);
     if (g <= 0.0){
       printf ("! E - m < 0 in PrimToCons (%12.6e)\n",g);
       QUIT_PLUTO(1);
     }
     g = u[RHO]/u[ENG] - sqrt(1.0 - m2/u[ENG]/u[ENG]);
     if (g >=0){
        print ("! g > 0.0 in PrimToCons(2) (%12.6e)\n",g);
        Show(ucons,ii);
        QUIT_PLUTO(1);
     }
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
  int     nv, i, status=0, use_energy;
  double  scrh, m, g;
  double *u, *v;
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
    par.m2 = EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]);

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
      u[ENG]   = sqrt(1.e-8 + par.m2 + u[RHO]*u[RHO]);
      flag[i] |= ENG_FAIL;
    }

/*
    scrh = sqrt(m2 + u[RHO]*u[RHO]);
    if (u[ENG] < scrh){
      WARNING( 
        print ("! ConsToPrim: E*E < m*m + D*D, ");
        Where(i, NULL);
      )
      u[ENG]   = sqrt(scrh*scrh + 1.e-9); 
      u[ENG]   = sqrt(scrh*scrh + g_smallPressure*g_smallPressure);
      flag[i] |= ENG_FAIL;
    }
*/
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
         if (PressureFix(&par) != 0) flag[i] = PRS_FAIL;
         u[ENG] = par.E;
       }
     } 
    #endif
 
    if (use_energy){
      if (EnergySolve (&par) != 0){
        WARNING(Where(i,NULL);)
        flag[i] |= ENG_FAIL;
        if (PressureFix(&par) != 0) flag[i] |= PRS_FAIL;
        u[ENG] = par.E;
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

  /*  ------------------------------------------
              complete conversion 
      ------------------------------------------ */

    v[PRS] = par.prs;
    scrh   = 1.0/(u[ENG] + v[PRS]);  /* = 1 / W */

    EXPAND(v[VX1] = u[MX1]*scrh; ,
           v[VX2] = u[MX2]*scrh; ,
           v[VX3] = u[MX3]*scrh;)

    g = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
    g = 1.0/sqrt(1.0 - g);
    #if USE_FOUR_VELOCITY == YES
     EXPAND(v[VX1] *= g;  ,
            v[VX2] *= g;  ,
            v[VX3] *= g;)
    #endif
    v[RHO] = u[RHO]/g;

    for (nv = NFLX; nv < NVAR; nv++) v[nv] = u[nv]/u[RHO];
    status |= flag[i];
  }

  return(status);
}
