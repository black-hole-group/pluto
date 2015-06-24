/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Convert between primitive and conservative variables.

  The PrimToCons() converts an array of primitive quantities to 
  an array of conservative variables for the MHD equations.
  
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
void PrimToCons (double **uprim, double **ucons, int ibeg, int iend)
/*!
 * Convert primitive variables in conservative variables. 
 *
 * \param [in]  uprim array of primitive variables
 * \param [out] ucons array of conservative variables
 * \param [in]  beg   starting index of computation
 * \param [in]  end   final index of computation
 *
 *********************************************************************** */
{
  int  i, nv, status;
  double *v, *u;
  double rhoe, kin, T, gmm1;

  #if EOS == IDEAL
   gmm1 = g_gamma - 1.0;
  #endif
  for (i = ibeg; i <= iend; i++) {
  
    v = uprim[i];
    u = ucons[i];

    u[RHO] = v[RHO];
        
    EXPAND (u[MX1] = v[RHO]*v[VX1];  ,
            u[MX2] = v[RHO]*v[VX2];  ,
            u[MX3] = v[RHO]*v[VX3];)

    EXPAND (u[BX1] = v[BX1];  ,
            u[BX2] = v[BX2];  ,
            u[BX3] = v[BX3];)

    kin   = EXPAND(v[VX1]*v[VX1], + v[VX2]*v[VX2], + v[VX3]*v[VX3]);
    kin   = v[RHO]*kin + EXPAND(v[BX1]*v[BX1], + v[BX2]*v[BX2], + v[BX3]*v[BX3]);
    kin  *= 0.5;

    #if EOS == IDEAL
     u[ENG] = kin + v[PRS]/gmm1;
    #elif EOS == PVTE_LAW
     status = GetPV_Temperature(v, &T);
     if (status != 0){
       T      = T_CUT_RHOE;
       v[PRS] = Pressure(v, T);
     }
     rhoe   = InternalEnergy(v, T);
     u[ENG] = rhoe + kin;

     if (u[ENG] != u[ENG]){
       print("! PrimToCons: KE:%12.6e uRHO : %12.6e, m2 : %12.6e \n",rhoe,v[RHO],u[ENG]);
       QUIT_PLUTO(1);
     }
    #endif
    
    #ifdef GLM_MHD
     u[PSI_GLM] = v[PSI_GLM]; 
    #endif

    for (nv = NFLX; nv < (NFLX + NSCL); nv++) u[nv] = v[RHO]*v[nv];
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
  int  i, nv, status=0, use_energy;
  double tau, rho, gmm1, rhoe, T;
  double b2, m2, kin, rhog1;
  double *u, *v;

  #if EOS == IDEAL
   gmm1 = g_gamma - 1.0;
  #endif

  for (i = ibeg; i <= iend; i++) {

    flag[i] = 0;
    u = ucons[i];
    v = uprim[i];

    m2 = EXPAND(u[MX1]*u[MX1], + u[MX2]*u[MX2], + u[MX3]*u[MX3]);
    b2 = EXPAND(u[BX1]*u[BX1], + u[BX2]*u[BX2], + u[BX3]*u[BX3]);

  /* -------------------------------------------
           Check density positivity 
     ------------------------------------------- */
  
    if (u[RHO] < 0.0) {
      print("! ConsToPrim: negative density (%8.2e), ", u[RHO]);
      Where (i, NULL);
      u[RHO]   = g_smallDensity;
      flag[i] |= RHO_FAIL;
    }
/*
   #if COOLING != NO
    if (u[X_HI] < 0.0) {
      print("! ConsToPrim: negative fractions (%8.2e), ", u[X_HI]);
      Where (i, NULL);
      Show(ucons,i);
      QUIT_PLUTO(1);
    }
  #endif
*/
    v[RHO] = rho = u[RHO];
    tau = 1.0/u[RHO];
    EXPAND(v[VX1] = u[MX1]*tau;  ,
           v[VX2] = u[MX2]*tau;  ,
           v[VX3] = u[MX3]*tau;)

    EXPAND(v[BX1] = u[BX1];  ,
           v[BX2] = u[BX2];  ,
           v[BX3] = u[BX3];)

  /* --------------------------------------------
       now try to recover pressure from 
       total energy or entropy
     -------------------------------------------- */

    #if EOS == IDEAL
     kin = 0.5*(m2*tau + b2);
     use_energy = 1;     
    
     #if ENTROPY_SWITCH == YES
      #ifdef CH_SPACEDIM
       if (g_intStage > 0)   /* -- HOT FIX used with Chombo: avoid calling 
                              CheckZone when writing file to disk          -- */
      #endif
       if (CheckZone(i,FLAG_ENTROPY)){
         use_energy = 0;
         rhog1 = pow(rho, g_gamma - 1.0);
         v[PRS] = u[ENTR]*rhog1; 
         if (v[PRS] < 0.0){
            WARNING(
              print("! ConsToPrim: negative p(S) (%8.2e, %8.2e), ", v[PRS], u[ENTR]);
              Where (i, NULL);
            )
          v[PRS]   = g_smallPressure;
          flag[i] |= PRS_FAIL;
         }
         u[ENG] = v[PRS]/gmm1 + kin; /* -- redefine energy -- */
       }
     #endif  /* ENTROPY_SWITCH == YES */

     if (use_energy){
       if (u[ENG] < 0.0) {
         WARNING(
           print("! ConsToPrim: negative energy (%8.2e), ", u[ENG]);
           Where (i, NULL);
         )
         v[PRS]    = g_smallPressure;
         u[ENG]    = v[PRS]/gmm1 + kin;
         flag[i] |= ENG_FAIL;
       }else{
         v[PRS] = gmm1*(u[ENG] - kin);
         if (v[PRS] < 0.0){
           WARNING(
             print("! ConsToPrim: negative p(E) (%8.2e), ", v[PRS]);
             Where (i, NULL);
           )
           v[PRS]    = g_smallPressure;
           flag[i] |= PRS_FAIL;
           u[ENG]    = v[PRS]/gmm1 + kin; /* -- redefine energy -- */
         }
       }
     }

    #endif  /* EOS == IDEAL */

  /* -- do remaning variables -- */

    for (nv = NFLX; nv < NVAR; nv++) v[nv] = u[nv]*tau;

/* --------------------------------------------
      Recover pressure from energy or entropy:
      2. PVTE_LAW Equation of state
     -------------------------------------------- */

    #if EOS == PVTE_LAW
     if (u[ENG] != u[ENG]){
       print("! ConsToPrim: NaN found\n");
       Show(ucons,i);
       QUIT_PLUTO(1);
     }
     kin  = 0.5*(m2*tau + b2);      
     rhoe = u[ENG] - kin; 

     status = GetEV_Temperature (rhoe, v, &T);
     if (status != 0){  /* If something went wrong while retrieving the  */
                        /* temperature, we floor \c T to \c T_CUT_RHOE,  */
                        /* recompute internal and total energies.        */
       T = T_CUT_RHOE;
       WARNING(  
         print ("! ConsToPrim: rhoe < 0 or T < T_CUT_RHOE; ");
         Where(i,NULL);
       )
       rhoe     = InternalEnergy(v, T);
       u[ENG]   = rhoe + kin; /* -- redefine total energy -- */
       flag[i] |= PRS_FAIL;
     }

     v[PRS] = Pressure(v, T);
    #endif  /* EOS == PVTE_LAW */

    #ifdef GLM_MHD
     v[PSI_GLM] = u[PSI_GLM]; 
    #endif

    status |= flag[i];
  }

  return(status);
}
