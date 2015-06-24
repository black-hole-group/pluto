/* ///////////////////////////////////////////////////////////////////// */
/*! 
 \file
 \brief Compute the flux for the relativistic hydro (RHD) equations.

  Compute the flux of the conservative RHD equations in the direction 
  given by ::g_dir.
  This function defines the component of the hyperbolic flux tensor 
  of the standard RHD equations.\n
  In what follows:
  - \c VXn, \c MXn are the velocity, momentum components in the direction 
    given by ::g_dir (normal, \c "n")
  - \c VXt, \c MXt and \c VXb, \c MXb are the transverse components
    (tangent \c "t" and bi-tangent \c "b").

 \author A. Mignone (mignone@ph.unito.it)
 \date   Aug 16, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Flux (double **u, double **v, double *a2, double **fx, 
           double *p, int beg, int end)
/*
 *
 *
 *
 *********************************************************************** */
{
  int    nv, i;
  double vn;

  for (i = beg ; i <= end; i++) {
    #if USE_FOUR_VELOCITY == YES
     vn = EXPAND( v[i][VXn]*v[i][VXn], 
                + v[i][VXt]*v[i][VXt],
                + v[i][VXb]*v[i][VXb]);

     vn = v[i][VXn]/sqrt(1.0 + vn);
    #else  
     vn = v[i][VXn];
    #endif

    fx[i][RHO]  = u[i][RHO]*vn;
    EXPAND(fx[i][MX1] = u[i][MX1]*vn;  ,
           fx[i][MX2] = u[i][MX2]*vn;  ,
           fx[i][MX3] = u[i][MX3]*vn;)
    fx[i][ENG] = u[i][MXn];
    p[i] = v[i][PRS];
  }
}
