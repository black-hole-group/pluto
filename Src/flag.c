/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Flag database.

  The FlagReset() function is simply used to reset all the bits 
  ot the flag array to zero (immediately before starting integration) and
  keep a static copy of its memory address.
  
  The CheckZone() function is used to check whether a zone has been 
  tagged with a particular flag using the Bit flag labels.

  \authors A. Mignone (mignone@ph.unito.it)
  \date    Oct 29, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#if (defined CH_SPACEDIM) && (TIME_STEPPING == RK2)
 unsigned char ***flag;
#else
 static unsigned char ***flag;
#endif

/* ********************************************************************* */
void FlagReset (const Data *d)
/*!
 *
 *  Set d->flag array to zero and keep a static reference.
 *
 *********************************************************************** */
{
  int  i, j, k;
  
  #if (defined CHOMBO) && (TIME_STEPPING == RK2)
   flag = d->flag;
   if (g_intStage == 1) TOT_LOOP(k,j,i) flag[k][j][i] = 0;
  #else
   if (flag == NULL) flag = d->flag;
   TOT_LOOP(k,j,i) flag[k][j][i] = 0;
  #endif

}

/* ********************************************************************* */
unsigned char CheckZone (int z, int bit)
/*!
 * Check if bit "bit" is turned on in zone "z".
 *
 *********************************************************************** */
{
  if (flag == NULL){

  /* ----------------------------------------------
      the flag array is not defined when Chombo 
      is converting conservative variables to 
      primitive ones.
     ---------------------------------------------- */  
    #ifdef CH_SPACEDIM
     return (0);
    #endif
    print1 ("! CheckZone: NULL pointer in CheckZone\n");
    QUIT_PLUTO(1);
  }

  if (g_dir == IDIR) return (flag[g_k][g_j][z] & bit);
  if (g_dir == JDIR) return (flag[g_k][z][g_i] & bit);
  if (g_dir == KDIR) return (flag[z][g_j][g_i] & bit);
  else{
    print1 ("! CheckZone: wrong direction\n");
    QUIT_PLUTO(1);
  }

}

