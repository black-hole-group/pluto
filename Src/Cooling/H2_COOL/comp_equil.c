/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Compute equilibrium fractions for the H2_COOL module.

  Compute the equilibrium fractions of ions for a given density
  and temperature.
  This is accomplished by solving one of the rate equations 
  with (almost) fixed temperature and replacing the other two equations 
  with the normalization condition (HI + HII + 2*H2 = 1) and 
  the fact that HI = cr/ci*HII. 

  \authors A. Mignone (mignone@ph.unito.it)\n
           B. Vaidya 

  \date    Sept 10, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
static double T_;
static void RateFunc (int, double *, double *);

void RateFunc (int n, double *val, double *fval)
{
  double T, st, tev;
  double cr, ci, kr1, kr2, kr3, kr4;
  double frac_Z, x;
  
  T = T_;
  cr  = 2.6e-11/sqrt(T);
  ci  = 1.08e-8*sqrt(T)*exp(-157890.0/T)/(13.6*13.6);
  frac_Z = ((1.0 - H_MASS_FRAC - He_MASS_FRAC)/CONST_AZ)*(CONST_AH/H_MASS_FRAC);

  st  = sqrt(T);
  tev = T*CONST_kB/CONST_eV;
  
  kr1 = 3.e-18*st/(1.0 + 0.04*st + 2e-3*T + 8.e-6*T*T); /* fa ~ 1 and Tg << T  */
  kr2 = 1.067e-10*pow(tev,2.012)*exp(-(4.463/tev)*pow((1.0 + 0.2472),3.512));
  kr3 = 1.e-8*exp(-84100.0/T); 
  kr4 = 4.4e-10*pow(T,0.35)*exp(-102000.0/T);

  //double gn = (1.0 - val[0] - val[1])*0.5;
  x = (val[2] + 0.5*CONST_AZ*frac_Z);  /* -- electron number density, in cm^{-3} -- */
  /*
  fval[0] = (kr1*val[0]*val[0] - gn*(kr2*val[0] + kr3*gn + kr4*x))/kr1;
  fval[1] = (cr*val[1] - ci*val[0])/MAX(ci,cr);
  */
     
  fval[0] = (kr1*val[0]*val[0] - val[1]*(kr2*val[0] + kr3*val[1] + kr4*x))/kr1;
  fval[1] = val[1] - (1.0 - val[0] - val[2])*0.5;
  fval[2] = (cr*val[2] - ci*val[0])/MAX(ci,cr);
  

  /*fval[1] = val[1]*(kr2*val[0] + kr3*val[1] + kr4*x) + cr*val[2]*x - val[0]*(ci*x + kr1*val[0]);
    fval[1] /= kr1;*/
}


/* ********************************************************************* */
void CompEquil(double n, double T, double *v)
/*!
 *  \param [in]      n  the particle number density (not needed)
 *  \param [in]      T  the temperature (in K) for which equilibrium must 
 *                      be found.
 *  \param [in/,out] v  an array of primitive variables. 
 *                      On input, only density needs to be defined.
 *                      On output, fractions will be updated to the 
 *                      equilibrium values.
 *
 *********************************************************************** */
{
  int i, j, sing, check, neq = 3;
  double x[3] = {0.1, 0.4, 0.1};
  T_  = T;
  Broyden(x, neq, &check, RateFunc);
  v[X_HI] = x[0];
  v[X_H2] = x[1];
  v[X_HII] = x[2];
}




