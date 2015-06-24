/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Solve quartic and cubic equations-
  

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 10, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     These two functions may be moved into tools.c
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

#define swap(x,y) f = x; x = y; y = f;

/* ********************************************************************* */
int QuarticSolve (double b, double c, double d, double e, double *z)
/*!
 * Solve a quartic equation in the form 
 * \f[
 *      z^4 + bz^3 + cz^2 + dz + e = 0
 * \f]
 *  For its purpose, it is assumed that \b ALL roots are double. 
 *  This makes things faster.
 *
 * \param [in] b   coefficient of the quartic
 * \param [in] c   coefficient of the quartic
 * \param [in] d   coefficient of the quartic
 * \param [in] e   coefficient of the quartic
 * \param [out] z  vector containing the 
 *                 (double) roots of the quartic
 *   
 * \b Reference:
 *
 *   http://www.1728.com/quartic2.htm 
 * 
 * \return Return 0 on success.
 *
 *********************************************************************** */
{
  int  n, j, ifail;
  double b2, f, g, h;
  double a2, a1, a0, u[4];
  double p, q, r, s;
  static double three_256 = 3.0/256.0;
  static double one_64 = 1.0/64.0;
  
  b2 = b*b;

  f = c - b2*0.375;
  g = d + b2*b*0.125 - b*c*0.5;
  h = e - b2*b2*three_256 + 0.0625*b2*c - 0.25*b*d;
  
  a2 = 0.5*f;
  a1 = (f*f - 4.0*h)*0.0625;
  a0 = -g*g*one_64;

  ifail = CubicSolve(a2, a1, a0, u);

  if (ifail)return(1);

  if (u[1] < 1.e-14){

    p = sqrt(u[2]);
    s = 0.25*b;
    z[0] = z[2] = - p - s;
    z[1] = z[3] = + p - s;
    
  }else{
  
    p = sqrt(u[1]);
    q = sqrt(u[2]);
  
    r = -0.125*g/(p*q);
    s =  0.25*b;
     
    z[0] = - p - q + r - s;
    z[1] =   p - q - r - s;
    z[2] = - p + q - r - s;
    z[3] =   p + q + r - s;

  }  

/* ----------------------------------------
    Sort roots in ascending order.
    Since z[0] < z[3], we first order

     z[0] < z[1] < z[3]

    and then put z[2] into the sequence.
   ---------------------------------------- */

  if      (z[1] < z[0]) {swap(z[1],z[0]);}
  else if (z[1] > z[3]) {swap(z[1],z[3]);}
 
 
  if (z[2] < z[0]) {
    swap(z[2],z[0]);
    swap(z[2],z[1]);
  }else if (z[2] < z[1]){
    swap(z[2],z[1]);
  }else if (z[2] > z[3]){
    swap(z[2],z[3]);
  }

/*
  for (j = 1; j < 4; j++){
    f = z[j];
    n = j - 1;
    while(n >= 0 && z[n] > f) {
       z[n+1] = z[n];
       n--;
    }
    z[n+1] = f;
  }
*/

  /* ----------------------------------------------
       verify that cmax and cmin satisfy original 
       equation
     ---------------------------------------------- */  

  for (n = 0; n < 4; n++){
    s = e + z[n]*(d + z[n]*(c + z[n]*(b + z[n])));
    if (s != s) {
      print ("! QuarticSolve: NaN found.\n");
      return(1);
    }
    if (fabs(s) > 1.e-6) {
      print ("! QuarticSolve: solution does not satisfy f(z) = 0; f(z) = %12.6e\n",s);
      return(1);
    }
  }

  return(0);
/*  
  printf (" z: %f ; %f ; %f ; %f\n",z[0], z[1], z[2], z[3]);
  */
}
#undef swap
/* ********************************************************************* */
int CubicSolve (double b, double c, double d, double z[])
/*!
 *  Solve a cubic equation in the form 
 *  \f[
 *      z^3 + bz^2 + cz + d = 0
 *  \f]
 *  For its purpose, it is assumed that \b ALL roots are double. 
 *  This makes things faster.
 *
 * \param [in] b coefficient of the cubic
 * \param [in] c coefficient of the cubic
 * \param [in] d coefficient of the cubic
 * \param [out] z vector containing the roots of the cubic.
 *                Roots should be sorted in increasing order.
 *   
 * \b Reference:  http://www.1728.com/cubic2.htm 
 *
 * \return Return 0 on success.
 *
 *********************************************************************** */
{
  double b2, g2;
  double f, g, h;
  double i, i2, j, k, m, n, p;
  static double one_3 = 1.0/3.0, one_27=1.0/27.0;

  b2 = b*b;
 
/*  ----------------------------------------------
     the expression for f should be 
     f = c - b*b/3.0; however, to avoid negative
     round-off making h > 0.0 or g^2/4 - h < 0.0
     we let c --> c(1- 1.1e-16)
    ---------------------------------------------- */

  f  = c*(1.0 - 1.e-16) - b2*one_3;
  g  = b*(2.0*b2 - 9.0*c)*one_27 + d; 
  g2 = g*g;
  i2 = -f*f*f*one_27;
  h  = g2*0.25 - i2;

/* --------------------------------------------
     double roots are possible only when 
   
               h <= 0 
   -------------------------------------------- */

  if (h > 1.e-12){
    printf ("Only one double root (%12.6e)!\n", h);
  }
  if (i2 < 0.0){
/*
    printf ("i2 < 0.0 %12.6e\n",i2);
    return(1);
*/
    i2 = 0.0;
  }

/* --------------------------------------
       i^2 must be >= g2*0.25
   -------------------------------------- */
  
  i = sqrt(i2);       /*  > 0   */
  j = pow(i, one_3);  /*  > 0   */
  k = -0.5*g/i;

/*  this is to prevent unpleseant situation 
    where both g and i are close to zero       */

  k = (k < -1.0 ? -1.0:k);
  k = (k >  1.0 ?  1.0:k);
  
  k = acos(k)*one_3;       /*  pi/3 < k < 0 */
 
  m = cos(k);              /*   > 0   */
  n = sqrt(3.0)*sin(k);    /*   > 0   */
  p = -b*one_3;

  z[0] = -j*(m + n) + p;
  z[1] = -j*(m - n) + p;
  z[2] =  2.0*j*m + p;

/* ------------------------------------------------------
    Since j, m, n > 0, it should follow that from
    
      z0 = -jm - jn + p
      z1 = -jm + jn + p
      z2 = 2jm + p
      
    z2 is the greatest of the roots, while z0 is the 
    smallest one.
   ------------------------------------------------------ */
      
  return(0);
}

