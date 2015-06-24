/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Inversion scheme for RMHD using a pressure fix.

  Fix p to a small value, solve for the square of velocity by using
  secant algorithm applied to Eq (A3).
  This step involved re-computing W at each step of the iteration.
  Once the root has been found, we recompute total energy E.
  Return 0 if succesful, 1 otherwise.

  \authors A. Mignone \n
           C. Zanni

  \date Oct 5, 2012
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

#define MAX_ITER 50
static double VelocitySquareFunc(double, Map_param *);
/* ********************************************************************* */
int PressureFix(Map_param *par)
/*!
 *
 *********************************************************************** */
{
  int    k, done=0;
  double v2, v2c, fc, f, dW, S2_W2;
  double fmin, fmax, v2min, v2max;
  
  par->prs = g_smallPressure; 

  v2max = 1.0-1.e-8;
  v2c = 0.95;
  fc  = VelocitySquareFunc(v2c, par);
  v2  = 0.96;
  for (k = 1; k < MAX_ITER; k++){
    f   = VelocitySquareFunc(v2, par);
    if (done == 1) break;
    dW  = (v2 - v2c)/(f - fc)*f;
    v2c = v2; fc = f;
    v2 -= dW;
    v2 = MIN(v2max,v2);
    v2 = MAX(v2, 0.0);
    if (fabs(f) < 1.e-9) done = 1;
  }
  if (v2 >= 1.0 || k >= MAX_ITER) {
    print ("! PressureFix: too many iter while fixing p , v^2 = %f\n", v2);
    return (1);
  }
/*
v2min = 0.0;
for (k = 1; k < MAX_ITER; k++){
  v2c = 0.5*(v2min + v2max);
  fc  = FUNV2(v2c, p, u, m2, S2, Bmag2); 
  if (fc*fmin > 0.0){
    v2min = v2c; fmin  = fc;
  }else{
    v2max = v2c; fmax  = fc;
  }
  if (fabs(fc) < 1.e-9) break;
}    
if (fabs (v2c-v2) > 1.e-8) {
 print ("! Solution mismatch\n");
 QUIT_PLUTO(1);
}
*/

/* -----------------------------------------------------
         redefine energy and proper density 
   ----------------------------------------------------- */
  
  S2_W2 = par->S2/(par->W*par->W); 
  #if SUBTRACT_DENSITY == YES
   par->E = par->W - par->D - par->prs + 0.5*(1.0 + v2)*par->B2 - 0.5*S2_W2;
  #else
   par->E = par->W - par->prs + 0.5*(1.0 + v2)*par->B2 - 0.5*S2_W2;
  #endif
  par->rho = par->D/par->lor;

  return(0);  /* -- success -- */
} 

/* ****************************************************************** */
double VelocitySquareFunc (double v2, Map_param *par)
/*!
 * 
 * Implement Eq (A3) of Mignone \& McKinney (2007).
 * 
 ******************************************************************** */
{
  double lor2, W2, f, pg;

  lor2     = 1.0/(1.0 - v2);
  par->lor = sqrt(lor2);
  pg       = par->prs*par->lor;
  #if EOS == IDEAL
   par->W = (par->D + pg*g_gamma/(g_gamma - 1.0))*par->lor;
  #elif EOS == TAUB
   par->W = (2.5*pg + sqrt(2.25*pg*pg + par->D*par->D))*par->lor;
  #endif
  
  W2  = par->W*par->W;

  f  =  par->S2*(2.0*par->W + par->B2) + par->m2*W2;
  f /= (par->W + par->B2)*(par->W + par->B2)*W2;
  f -= v2;

  return f;
}
#undef MAX_ITER