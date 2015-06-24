/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Return the number of ghost zones.

  Return the number of ghost zones employed by the selected numerical
  algorithm. 
  The minimum number for a 2nd-order algorithm is 2.
  Higher-order interpolation scheme may require more zones.
  
  \authors A. Mignone (mignone@ph.unito.it)
  \date    Sep 20, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
int GetNghost (Input *ini)
/*! 
 * Compute the number of ghost zones, depending on the selected
 * scheme.
 *
 *********************************************************************** */
{
  int nv, nghost ;
  Limiter *limiter[NVAR];

  #if INTERPOLATION == FLAT

   nghost = 2;

  #elif INTERPOLATION == LINEAR || INTERPOLATION == LimO3\
                                || INTERPOLATION == WENO3
   #if LIMITER == FOURTH_ORDER_LIM
    nghost = 3;
   #else
    nghost = 2;
   #endif

  #elif INTERPOLATION == LINEAR_MULTID

   nghost = 2;

  #elif INTERPOLATION == PARABOLIC

   #if PHYSICS == HD || PHYSICS == RHD
    nghost = 4;   /* -- since HD makes use of contact steepener -- */
   #else
    nghost = 3;
   #endif

  #elif INTERPOLATION == WENO3_FD || INTERPOLATION == LIMO3_FD

   nghost = 2;

  #elif INTERPOLATION == WENOZ_FD  || \
        INTERPOLATION == MP5_FD

   nghost = 3;

  #endif

  #if SHOCK_FLATTENING == ONED

   nghost = MAX(4, nghost);

  #elif SHOCK_FLATTENING == MULTID

/* ------------------------------------------------------
    The MULTID shock flattening only need 2 ghost zones.
    However for axisymmetric simulations with CTU 
    3 zones will ensure that flag[][][] will remain 
    symmetric around the axis.
   ------------------------------------------------------ */

   nghost = MAX(3, nghost);
  #endif

/* ----------------------------------------------------
    The following should operate on the static grid 
    version of the code. Add an extra row of boundary
    zones if CTU+CT is selected.
    At least 3 ghost zones.
   ---------------------------------------------------- */

  #ifdef CTU
   #ifdef STAGGERED_MHD
    nghost++;
   #endif
  #endif

/* ----------------------------------------
    FARGO PPM needs at least 3 ghost zones 
   ---------------------------------------- */

  #ifdef FARGO
   #if FARGO_ORDER == 3
    nghost = MAX(3,nghost);
   #endif
  #endif

  #if (defined CH_SPACEDIM) && (TIME_STEPPING == RK_MIDPOINT) 
   nghost++;  /* AMR + RK_MIDPOINT */
  #endif  

  return (nghost);
}

