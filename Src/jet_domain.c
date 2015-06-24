/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Adjust the effective domain size when the -xNjet option is given.

  The SetJetDomain() function shortens the domain integration index in the 
  direction specified by either '-x1jet', '-x2jet' or '-x3jet' to save 
  computational time.
  This is done by setting the final integration index in the direction
  of jet propagation to that of the first zone (counting from the top) with 
  non-zero pressure gradient (GetRightmostIndex()).
  The function UnsetJetDomain() restores the initial computational domain
  size.
  Useful for problems involving jet propagation.
  \note In parallel, the domain shall not be decomposed along the
        propagation direction.
  
  \author A. Mignone (mignone@ph.unito.it)
  \date   Oct 3, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */

#include "pluto.h"

static int NBEG, NEND, NPT, NPT_TOT, rbound;
static int GetRightmostIndex(int, double ***);

/* ********************************************************************* */
void SetJetDomain (const Data *d, int dir, int log_freq, Grid *grid)
/*!
 *
 * Adjust the size of computational domain by reducing the final
 * computationa index in the direction of propagation.
 *
 * \param [in]         d   pointer to Data structure
 * \param [in]       dir   the direction of propagation
 * \param [in]  log_freq   the output log frequency
 * \param [in,out]  grid   pointer to array of Grid structures
 * 
 *********************************************************************** */
{
  int i, j, k, ngh;
  int n, n_glob;
  static int first_call = 1;
  double ***pr, ***dn, dp;

  dn = d->Vc[RHO];
  #if HAVE_ENERGY
   pr = d->Vc[PRS];
  #else
   pr = d->Vc[RHO];
  #endif

  ngh = grid[dir].nghost;

/* -- save original domain offsets 
      and return for the first time -- */

  if (first_call){
    if (dir == IDIR) {
      NBEG = IBEG; NEND    = IEND; 
      NPT  = NX1;  NPT_TOT = NX1_TOT;
    }else if (dir == JDIR){
      NBEG = JBEG; NEND    = JEND; 
      NPT  = NX2;  NPT_TOT = NX2_TOT;
    }else if (dir == KDIR){
      NBEG = KBEG; NEND    = KEND; 
      NPT  = NX3;  NPT_TOT = NX3_TOT;
    }

    rbound = grid[dir].rbound;
    first_call = 0;
    return;
  }

/* --------------------------------------
    find where grad(p) is not zero and 
    add safety guard cells 
   -------------------------------------- */

  n = GetRightmostIndex(dir, pr) + 2*ngh;

  if (n < NBEG + ngh) n = NBEG + ngh;
  if (n > NEND)       n = NEND;

  #ifdef PARALLEL
   MPI_Allreduce (&n, &n_glob, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   n = n_glob;
  #endif

  if (g_stepNumber%log_freq==0){
/*    print1 ("- SetJetDomain: index %d / %d\n",n,NEND); */
  }

/* ------------------------------------------------------------------
    Change global integer variables giving the rightmost integration 
    indexes. Also, suppress the rightmost physics boundary (rbound) 
    if solution has not yet reached it. 
   ------------------------------------------------------------------ */

  if (dir == IDIR) {
    IEND = grid[IDIR].lend = n;
    NX1_TOT = IEND + ngh;
    NX1     = IEND - ngh;

    grid[IDIR].rbound = rbound;
    if (IEND != NEND) grid[IDIR].rbound = 0;  
    
  } else if (dir == JDIR){

    JEND = grid[JDIR].lend = n;
    NX2_TOT = JEND + ngh;
    NX2     = JEND - ngh;
    grid[JDIR].rbound = rbound;
    if (JEND != NEND) grid[JDIR].rbound = 0;  

  }else if (dir == KDIR){

    KEND = grid[KDIR].lend = n;
    NX3_TOT = KEND + ngh;
    NX3     = KEND - ngh;

    grid[KDIR].rbound = rbound;
    if (KEND != NEND) grid[KDIR].rbound = 0;  
  }
  
/*
  print1 ("Box = %d/%d\n", n, NEND);
  {
    double t0, t1, i0, i1;
    FILE *fp;
    t0 = g_time;
    t1 = t0 + g_dt;
    i0 = (int)(t0/0.5);
    i1 = (int)(t1/0.5);

    if (i0 != i1){
      fp = fopen("box.out","w");
      fprintf (fp,"%f  %f  %f  %f\n",0.0, grid[IDIR].x[IEND], 
                                     0.0, grid[JDIR].x[JEND]);
      fclose(fp);
    } 
  }
*/
}

/* ********************************************************************* */
void UnsetJetDomain (const Data *d, int dir, Grid *grid)
/*!
 *  Restore original (full domain) indexes.
 *
 *********************************************************************** */
{

  if (dir == IDIR){
    IBEG = grid[IDIR].lbeg = NBEG;
    IEND = grid[IDIR].lend = NEND;
    NX1_TOT = NPT_TOT;
    NX1     = NPT;
    grid[IDIR].rbound = rbound;
  }

  if (dir == JDIR){
    JBEG = grid[JDIR].lbeg = NBEG;
    JEND = grid[JDIR].lend = NEND;
    NX2_TOT = NPT_TOT;
    NX2     = NPT;
    grid[JDIR].rbound = rbound;
  }

  if (dir == KDIR){
    KBEG = grid[KDIR].lbeg = NBEG;
    KEND = grid[KDIR].lend = NEND;
    NX3_TOT = NPT_TOT;
    NX3     = NPT;
    grid[KDIR].rbound = rbound;
  }

}

/* ********************************************************************* */
int GetRightmostIndex (int dir, double ***q)
/*!
 *
 * Find the local index j where grad(p) exceeds
 * a ceratin threshold
 *
 *********************************************************************** */
{
  int    i, j, k;
  double dp;

/* -- backward sweep on x-axis -- */

  if (dir == IDIR){
    for (i = NEND; i >= NBEG; i--){
    KDOM_LOOP(k){
    JDOM_LOOP(j){

      dp  = q[k][j][i+1] - q[k][j][i-1];
      dp /= q[k][j][i+1] + q[k][j][i-1];

      if (fabs(dp) > 1.e-5) return(i);
    }}}
  }

/* -- backward sweep on y-axis -- */

  if (dir == JDIR){
    for (j = NEND; j >= NBEG; j--){
    KDOM_LOOP(k){
    IDOM_LOOP(i){

      dp  = q[k][j+1][i] - q[k][j-1][i];
      dp /= q[k][j+1][i] + q[k][j-1][i];

      if (fabs(dp) > 1.e-5) return(j);
    }}}
  }

/* -- backward sweep on z-axis -- */

  if (dir == KDIR){
    for (k = NEND; k >= NBEG; k--){
    JDOM_LOOP(j){
    IDOM_LOOP(i){

      dp  = q[k+1][j][i] - q[k-1][j][i];
      dp /= q[k+1][j][i] + q[k-1][j][i];

      if (fabs(dp) > 1.e-5) return(k);
    }}}
  }

  return(-1);
}
