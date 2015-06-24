#include "pluto.h"

/* ****************************************************************** */
void SetDefaultVarNames(Output *output)
/*
 *
 *  PURPOSE
 *
 *    Set file names for I/O
 *
 *
 ******************************************************************** */
{
  int nv;

/* ----------------------------------------------
    Physics module file names; 
    these pertain to the physics module ONLY
   ---------------------------------------------- */

  output->var_name[RHO] = "rho";
  EXPAND(output->var_name[VX1] = "vx1";  ,
         output->var_name[VX2] = "vx2";  ,
         output->var_name[VX3] = "vx3";)
  #if HAVE_PRESSURE
   output->var_name[PRS] = "prs";
  #endif

  #if PHYSICS == MHD || PHYSICS == RMHD
   EXPAND(output->var_name[BX1] = "bx1";  ,
          output->var_name[BX2] = "bx2";  ,
          output->var_name[BX3] = "bx3";)
  #endif
  
  /* (staggered field names are set in SetOutput) */

  #ifdef GLM_MHD
   output->var_name[PSI_GLM] = "psi_glm";
  #endif
  
/* ------------------------------------------------
                   Tracers 
   ------------------------------------------------ */

  for (nv = TRC; nv < TRC + NTRACER; nv++){
    sprintf (output->var_name[nv],"tr%d",nv - TRC + 1);
  } 

  #if ENTROPY_SWITCH == YES
   sprintf (output->var_name[ENTR],"entropy");
  #endif

/* ------------------------------------------------
               Cooling vars
   ------------------------------------------------ */

  #if COOLING == MINEq
  {
   static char *ion_name[] = {"X_HI", "X_HeI", "X_HeII" 
                       C_EXPAND("X_CI","fCII", "X_CIII", "X_CIV", "X_CV")
                       N_EXPAND("X_NI","fNII", "X_NIII", "X_NIV", "X_NV")
                       O_EXPAND("X_OI","fOII", "X_OIII", "X_OIV", "X_OV")
                      Ne_EXPAND("X_NeI","X_NeII", "X_NeIII", "X_NeIV", "X_NeV")
                       S_EXPAND("X_SI","X_SII", "X_SIII", "X_SIV", "X_SV")
                      Fe_EXPAND("X_FeI", "X_FeII", "X_FeIII")};

   for (nv = 0; nv < NIONS; nv++) output->var_name[NFLX+nv] = ion_name[nv];  

  }
  #elif COOLING == SNEq

   output->var_name[X_HI] = "X_HI";

  #elif COOLING == H2_COOL

   static char *molnames[] = {"X_HI", "X_H2", "X_HII"};
   
   for (nv = 0; nv < NIONS; nv++) output->var_name[NFLX+nv] = molnames[nv];

  #endif

}
