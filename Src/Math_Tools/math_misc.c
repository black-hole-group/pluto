/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief  Miscellaneous math functions.
  \author A. Mignone (mignone@ph.unito.it)
  \date   July 6, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
double BesselJ0(double x)
/*!
 * Returns the Bessel function J0(x) for any real x.
 * (Adapted from Numerical Recipes)
 * 
 *********************************************************************** */
{
  double ax,z;
  double xx,y,ans,ans1,ans2;      /* Accumulate pol.in double prec. */
   
  if ((ax=fabs(x)) < 8.0) {       /* Direct rational function fit. */
    y   = x*x;
    ans1 =  57568490574.0+y*(-13362590354.0+y*(651619640.7
          + y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
    ans2 =  57568490411.0+y*(1029532985.0+y*(9494680.718
          + y*(59272.64853+y*(267.8532712+y*1.0))));
    ans = ans1/ans2;
  } else {                                  /*Fitting function (6.5.9).*/
    z    = 8.0/ax;
    y    = z*z;
    xx   = ax-0.785398164;
    ans1 = 1.0 + y*(-0.1098628627e-2+y*(0.2734510407e-4
           + y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 = - 0.1562499995e-1+y*(0.1430488765e-3
           + y*(-0.6911147651e-5+y*(0.7621095161e-6
           - y*0.934945152e-7)));
    ans = sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}

/* ********************************************************************* */
double BesselJ1(double x)
/*!
 * Returns the Bessel function J1(x) for any real x.
 * (Adapted from Numerical Recipes)
 * 
 *********************************************************************** */
{
  double ax,z;
  double xx,y,ans,ans1,ans2;    /* Accumulate polynomials in double precision.*/
  if ((ax=fabs(x)) < 8.0) {     /* Direct rational approximation.*/
  
    y=x*x;
    ans1 =   x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
           + y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
    ans2 =  144725228442.0+y*(2300535178.0+y*(18583304.74
          + y*(99447.43394+y*(376.9991397+y*1.0))));
    ans  = ans1/ans2;
  } else {                                  /*Fitting function (6.5.9).*/
    z    = 8.0/ax;
    y    = z*z;
    xx   = ax-2.356194491;
    ans1 = 1.0 +y*(0.183105e-2+y*(-0.3516396496e-4
           + y*(0.2457520174e-5+y*(-0.240337019e-6))));
     ans2 = 0.04687499995+y*(-0.2002690873e-3
           + y*(0.8449199096e-5+y*(-0.88228987e-6
           + y*0.105787412e-6)));
     ans  = sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
    if (x < 0.0) ans = -ans;
  }
  return ans;
}


/* ********************************************************************* */
double RandomNumber (double rmin, double rmax)
/*!
 * Generate and return a random number between [rmin, rmax] 
 *********************************************************************** */
{
  static int first_call = 1;
  double eps, rnd;

  if (first_call == 1){
    srand(time(NULL) + prank);
    first_call = 0;    
  }
  eps = (double)(rand())/((double)RAND_MAX + 1.0);  /* 0 < eps < 1 */
  rnd = rmin + eps*(rmax-rmin); /* Map random number between [rmin,rmax] */
  return rnd;
}

/* ********************************************************************* */
void VectorCartesianComponents(double *v, double x1, double x2, double x3)
/*!
 * Transform the vector \f$ v = (v_{x1}, v_{x2}, v_{x3})\f$ from
 * the chosen coordinate system (= GEOMETRY) to Cartesian components.
 *  
 * \param [in,out]  *v  the original vector.
 *                      On output \c v is replaced by the three Cartesian
 *                      components.
 * \param [in] x1,x2,x3  the grid coordinates
 *
 *********************************************************************** */
{
  double vx1, vx2, vx3;

  vx1 = v[0];
  vx2 = v[1];
  vx3 = v[2];

  #if GEOMETRY == CARTESIAN

  /* nothing to do */

  #elif GEOMETRY == POLAR

   EXPAND(v[0] = vx1*cos(x2) - vx2*sin(x2);   ,
          v[1] = vx1*sin(x2) + vx2*cos(x2);   ,
          v[2] = vx3;)
  
  #elif GEOMETRY == SPHERICAL

   #if DIMENSIONS == 2
    EXPAND(v[0] = vx1*sin(x2) + vx2*cos(x2);   ,
           v[1] = vx1*cos(x2) - vx2*sin(x2);   ,
           v[2] = vx3;)
   #elif DIMENSIONS == 3 
    v[0] = (vx1*sin(x2) + vx2*cos(x2))*cos(x3) - vx3*sin(x3);
    v[1] = (vx1*sin(x2) + vx2*cos(x2))*sin(x3) + vx3*cos(x3);
    v[2] = (vx1*cos(x2) - vx2*sin(x2));
   #endif
  #endif
}
    


