/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief 3D wrapper for conservative/primitive conversion.

  Provide 3D wrappers to the standard 1D conversion functions
  ConsToPrim() and PrimToCons().

  \authors A. Mignone (mignone@ph.unito.it)
  \date    March 10, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void ConsToPrim3D (Data_Arr U, Data_Arr V, Grid *grid)
/*!
 *  Convert a 3D array of conservative variables \c U to
 *  an array of primitive variables \c V.
 *  Note that <tt>[nv]</tt> is the fastest running index for \c U 
 *  while it is the slowest running index for \c V.
 *
 * \param [in]   U      pointer to 3D array of conserved variables,
 *                      with array indexing <tt>[k][j][i][nv]</tt>
 * \param [out]  V      pointer to 3D array of primitive variables,
 *                      with array indexing <tt>[nv][k][j][i]</tt>
 * \param [in]   grid   pointer to an array of Grid structures
 *********************************************************************** */
{
  int   i, j, k, nv;
  int   current_dir;
  static unsigned char *flag;
  static double **v;

  if (v == NULL){
    v    = ARRAY_2D(NMAX_POINT, NVAR, double);
    flag = ARRAY_1D(NMAX_POINT, unsigned char);   
  }

  current_dir = g_dir;  /* save current direction */
  g_dir = IDIR;
  KDOM_LOOP(k) { g_k = k;
  JDOM_LOOP(j) { g_j = j;
    ConsToPrim (U[k][j], v, IBEG, IEND, flag);
    IDOM_LOOP(i) VAR_LOOP(nv) V[nv][k][j][i] = v[i][nv];
  }}
  g_dir = current_dir;  /* restore current direction */
}
/* ********************************************************************* */
void PrimToCons3D (Data_Arr V, Data_Arr U, Grid *grid)
/*!
 *  Convert a 3D array of primitive variables \c V  to
 *  an array of conservative variables \c U.
 *  Note that <tt>[nv]</tt> is the fastest running index for \c U 
 *  while it is the slowest running index for \c V.
 *
 * \param [in]    V     pointer to 3D array of primitive variables,
 *                      with array indexing <tt>[nv][k][j][i]</tt>
 * \param [out]   U     pointer to 3D array of conserved variables,
 *                      with array indexing <tt>[k][j][i][nv]</tt>
 * \param [in]   grid   pointer to an array of Grid structures
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  int   current_dir;
  static double **v;

  if (v == NULL) v = ARRAY_2D(NMAX_POINT, NVAR, double);

/* ------------------------------------------------------------
     Convert solution vector from primitive to conservative
   ------------------------------------------------------------ */

  current_dir = g_dir; /* save current direction */
  g_dir = IDIR;
  KDOM_LOOP(k) { g_k = k;
  JDOM_LOOP(j) { g_j = j;
    IDOM_LOOP(i)  VAR_LOOP(nv) v[i][nv] = V[nv][k][j][i];
    PrimToCons (v, U[k][j], IBEG, IEND);
  }}
  g_dir = current_dir; /* restore current direction */
}
