#include "pluto.h"



/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  double c, s, v0=sqrt(5.0);
  double A=1.e-3, R=0.3,r;

  r = sqrt(x1*x1 + x2*x2);
  c = 2.0/sqrt(5.0);
  s = 1.0/sqrt(5.0);
/*
  c = 1.0/sqrt(2.0);
  s = 1.0/sqrt(2.0);
*/
  us[RHO] = 1.0;
  us[VX1] = v0*c;
  us[VX2] = v0*s;
  us[VX3] = 0.0;
  us[PRS] = 1.0;
  us[TR] = 0.0;

  #if PHYSICS == MHD || PHYSICS == RMHD

   us[BX1] = -A*x2/r*(r <= R);  /*  =   dAz/dy  */
   us[BX2] =  A*x1/r*(r <= R);  /*  = - dAz/dx  */
   us[BX3] = 0.0;

   us[AX1] = 0.0;
   us[AX2] = 0.0;
   us[AX3] = A*(R - r)*(r <= R);

  #endif
}
/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in/out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies on which side boundary conditions need 
 *                    to be assigned. side can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{ }

